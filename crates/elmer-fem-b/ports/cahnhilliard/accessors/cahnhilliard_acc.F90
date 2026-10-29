! T022-A CahnHilliard port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/CahnHilliard.F90, compiled against the built
! .mod files of the target configure (A2 shipped mechanism, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, list-directed timing reports,
! whole-array updates — stays inside Fortran, so the Rust side is pure C
! calls. The wrapper module mirrors CahnHilliardSolver's own locals (the
! current element, the LocalMatrix Nodes/IP locals, the entry-level
! Newton/Nonlin state) with identical lifetimes. All floating-point work
! (the 2x2-block Galerkin assembly, SUM reductions, SQRT, the WRITE+Info
! reports) and every Elmer library call stay in Fortran (same compiler),
! so bit-identity holds by construction wherever the call sequence
! matches. The `NewtonLinearization` flag the original sets-but-never-
! reads is replicated as module state with identical updates (no
! observable effect, documented in lib.rs).
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module cahnhilliard_acc
  use iso_c_binding
  use defutils
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(valuelist_t), pointer, save :: acc_params => null()
  type(element_t), pointer, save :: acc_elem => null()
  ! LocalMatrix context (mirrors its Nodes local and IP local, which are
  ! per-call automatics in the original; consumed before the next fill).
  type(nodes_t), target, save :: acc_nodes
  type(gaussintegrationpoints_t), save :: acc_ip
  ! Entry-level state (mirrors the solver body's locals).
  integer, save :: acc_newton_iter = 0
  integer, save :: acc_nonlin_iter = 1
  real(kind=dp), save :: acc_newton_tol = 0.0_dp
  real(kind=dp), save :: acc_nonlin_tol = 0.0_dp
  real(kind=dp), save :: acc_relchange = 0.0_dp
  logical, save :: acc_newton_lin = .false.
  logical, save :: acc_transient = .false.

contains

  ! -- entry state ------------------------------------------------------
  subroutine ch_begin(h_solver, transient) bind(c, name="ch_begin")
    type(c_ptr), value :: h_solver
    integer(c_int), value :: transient
    call c_f_pointer(h_solver, acc_solver)
    acc_transient = (transient /= 0)
  end subroutine ch_begin

  ! -- solver-params reads (exact statement sequence of the entry head) --
  subroutine ch_read_params() bind(c, name="ch_read_params")
    logical :: found
    acc_params => GetSolverParams()
    acc_newton_iter = GetInteger(acc_params, &
        'Nonlinear System Newton After Iterations', found)
    acc_newton_tol = GetConstReal(acc_params, &
        'Nonlinear System Newton After Tolerance', found)
    acc_nonlin_iter = GetInteger(acc_params, &
        'Nonlinear System Max Iterations', found)
    if (.not. found) acc_nonlin_iter = 1
    acc_nonlin_tol = GetConstReal(acc_params, &
        'Nonlinear System Convergence Tolerance', found)
    acc_newton_lin = .false.
  end subroutine ch_read_params

  function ch_nonlin_iter() result(n) bind(c, name="ch_nonlin_iter")
    integer(c_int) :: n
    n = acc_nonlin_iter
  end function ch_nonlin_iter

  ! -- assembly driver ---------------------------------------------------
  subroutine ch_default_initialize() bind(c, name="ch_default_initialize")
    call DefaultInitialize()
  end subroutine ch_default_initialize

  function ch_nof_active() result(n) bind(c, name="ch_nof_active")
    integer(c_int) :: n
    n = GetNOFActive()
  end function ch_nof_active

  ! -- BulkAssembly element body (exact transcription of LocalMatrix) ----
  !
  ! Two layers: the BIND(C) entry resolves the element and its sizes,
  ! then the inner body takes n/nd as dummies so the 2*nd automatic
  ! arrays are bound on entry exactly as in the original internal
  ! procedure (automatic bounds evaluate on entry; n/nd must be dummies).
  subroutine ch_local_elem(t) bind(c, name="ch_local_elem")
    integer(c_int), value :: t
    integer :: n, nd
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes(acc_elem)
    nd = GetElementNOFDOFs(acc_elem)
    call ch_local_body(n, nd)
  end subroutine ch_local_elem

  subroutine ch_local_body(n, nd)
    integer :: n, nd
    real(kind=dp), target :: mass(2*nd,2*nd), stiff(2*nd,2*nd)
    real(kind=dp) :: force(2*nd), load(2,nd)
    real(kind=dp), pointer :: a(:,:), m(:,:)

    real(kind=dp) :: basis(nd), dbasisdx(nd,3), detj, sol(2,nd), phi, s
    logical :: stat, isscalar, found, found_dens, found_scoeff, found_mobil
    integer :: i, j, p, q, t, dim
    type(gaussintegrationpoints_t) :: ip

    real(kind=dp) :: velo(3), h, mixenergydensity, mobility
    real(kind=dp) :: velo_n(3,n), h_n(n), mixenergydensity_n(n), &
        mobility_n(n), st_n(n)

    type(valuelist_t), pointer :: material

    type(nodes_t) :: nodes
    ! -- LocalMatrix transcription starts here --
    call GetElementNodes(nodes, acc_elem)

    material => GetMaterial(acc_elem)

    dim = CoordinateSystemDimension()

    velo_n = 0._dp
    do i = 1, dim
      velo_n(i,:) = GetReal(material, 'Convection Velocity '//i2s(i), found)
    end do

    h_n = GetReal(material, 'Interface Thickness', found)
    if (.not. found) then
      h_n = ElementDiameter(acc_elem, nodes) / 2
    end if

    mobility_n = GetReal(material, 'Mobility', found_mobil, acc_elem)

    mixenergydensity_n = GetReal(material, 'Mixing Energy Density', found_dens)
    if (.not. found_dens) then
      st_n = GetReal(material, 'Surface Tension Coefficient', found_scoeff)
    end if

    call GetVectorLocalSolution(sol, UElement=acc_elem)

    mass = 0.0_dp
    stiff = 0.0_dp
    force = 0.0_dp

    ! Numerical integration:
    ! ----------------------
    ip = GaussPoints(acc_elem)
    do t = 1, ip % n
      stat = ElementInfo(acc_elem, nodes, ip % u(t), ip % v(t), &
          ip % w(t), detj, basis, dbasisdx)

      velo = 0._dp
      do i = 1, dim
        velo(i) = sum(velo_n(i,:)*basis(1:n))
      end do

      h = sum(h_n*basis(1:n))

      if (found_mobil) then
        mobility = sum(mobility_n*basis(1:n))
      else
        mobility = h**2
      end if

      mixenergydensity = 1._dp
      if (found_dens) then
        mixenergydensity = sum(mixenergydensity_n*basis(1:n))
      else if (found_scoeff) then
        mixenergydensity = 3*h/2/sqrt(2._dp)*sum(st_n*basis(1:n))
      end if

      phi = sum(sol(1,1:nd)*basis(1:nd))

      s = ip % s(t) * detj

      do p = 1, nd
        i = 2*(p-1)
        do q = 1, nd
          j = 2*(q-1)
          m => mass(i+1:i+2,j+1:j+2)
          a => stiff(i+1:i+2,j+1:j+2)

          ! Phi-eq: DPhi/Dt - div(m*l/h**2*grad(theta)) = 0

          ! Time derivative
          ! ---------------
          m(1,1) = m(1,1) + s*basis(q)*basis(p)

          ! Convection
          ! ---------------
          a(1,1) = a(1,1) + s*sum(velo*dbasisdx(q,:))*basis(p)

          ! diffusion
          ! -----------------
          a(1,2) = a(1,2) + s*mobility*mixenergydensity/h**2* &
              sum(dbasisdx(q,:)*dbasisdx(p,:))

          ! Theta-eq: Theta + div(h**2*grad(Phi)) + Phi - Phi**3 = 0

          a(2,2) = a(2,2) + s*basis(q)*basis(p)
          a(2,1) = a(2,1) - s*h**2*sum(dbasisdx(q,:)*dbasisdx(p,:))
          a(2,1) = a(2,1) + s*(1-3*phi**2)*basis(q)*basis(p)
        end do

        force(i+2) = force(i+2) - s * 2*phi**3 * basis(p)
      end do
    end do

    if (acc_transient) then
      call Default1stOrderTime(mass, stiff, force, UElement=acc_elem)
    end if
    call DefaultUpdateEquations(stiff, force, UElement=acc_elem)
  end subroutine ch_local_body

  ! -- BoundaryAssembly (BoundaryLocalMatrix is an empty stub) -----------
  function ch_nof_boundary() result(n) bind(c, name="ch_nof_boundary")
    integer(c_int) :: n
    n = GetNOFBoundaryElements()
  end function ch_nof_boundary

  subroutine ch_boundary_elem(t) bind(c, name="ch_boundary_elem")
    integer(c_int), value :: t
    integer :: n, nd
    acc_elem => GetBoundaryElement(t)
    if (.not. ActiveBoundaryElement(acc_elem)) return
    n = GetElementNOFNodes(acc_elem)
    nd = GetElementNOFDOFs(acc_elem)
    ! The original's BoundaryLocalMatrix body is empty: no further action.
  end subroutine ch_boundary_elem

  ! -- shared assembly tail ----------------------------------------------
  subroutine ch_default_finish_bulk() bind(c, name="ch_default_finish_bulk")
    call DefaultFinishBulkAssembly()
  end subroutine ch_default_finish_bulk

  subroutine ch_default_finish() bind(c, name="ch_default_finish")
    call DefaultFinishAssembly()
  end subroutine ch_default_finish

  subroutine ch_default_dirichlet() bind(c, name="ch_default_dirichlet")
    call DefaultDirichletBCs()
  end subroutine ch_default_dirichlet

  ! -- solve + reports (exact statement sequence of the loop tail) -------
  function ch_solve() result(norm) bind(c, name="ch_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
    acc_relchange = acc_solver % variable % nonlinchange
    write(message, *) 'Result Norm   : ', norm
    call Info('Cahn-Hilliard interface Solver', message, level=4)
    write(message, *) 'Relative Change : ', acc_relchange
    call Info('Cahn-Hilliard interface Solver', message, level=4)
  end function ch_solve

  subroutine ch_newton_update(iter) bind(c, name="ch_newton_update")
    integer(c_int), value :: iter
    if (acc_relchange < acc_newton_tol .or. iter >= acc_newton_iter) &
        acc_newton_lin = .true.
  end subroutine ch_newton_update

  function ch_converged() result(r) bind(c, name="ch_converged")
    integer(c_int) :: r
    r = merge(1, 0, acc_relchange < acc_nonlin_tol)
  end function ch_converged

end module cahnhilliard_acc
