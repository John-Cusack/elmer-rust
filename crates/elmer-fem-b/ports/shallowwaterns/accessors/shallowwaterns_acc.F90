! T021-B ShallowWaterNS port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/ShallowWaterNS.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module mirrors the original's locals (current element, Nodes SAVE
! state per internal procedure, IP shallow copies, material/body-force
! handles, per-element work arrays) with identical lifetimes.
!
! The whole LocalMatrix / BoundaryLocalMatrix numerical cores (SUM and
! MATMUL reductions, whole-array updates, pointer-aliased block
! assembly, the dim==1 bubble fixup) are transcribed statement for
! statement into swa_local_assemble / swa_boundary_assemble so
! bit-identity holds by construction wherever the call sequence
! matches. Clean-room code: statement-for-statement wrappers around
! Elmer library calls, no upstream text.
module shallowwaterns_acc
  use iso_c_binding
  use types
  use defutils
  use integration
  use matrixassembly
  use eleminfo
  use coordinatesystems
  use elementgeometry
  use messages
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(valuelist_t), pointer, save :: acc_prm => null()
  type(valuelist_t), pointer, save :: acc_mat => null()
  type(valuelist_t), pointer, save :: acc_bf => null()
  type(valuelist_t), pointer, save :: acc_bc => null()
  ! Shared element context (mirrors the internal procedures' SAVE Nodes
  ! and IP locals; consumed before the next fill, as in the original).
  type(nodes_t), target, save :: acc_nodes
  type(nodes_t), target, save :: acc_pnodes
  type(gaussintegrationpoints_t), save :: acc_ip
  ! LocalMatrix / BoundaryLocalMatrix work state (mirrors the automatic
  ! arrays, sized per element; TARGET for the M=>/A=> block aliasing).
  real(c_double), allocatable, target, save :: acc_sol(:, :)     ! (3, nd)
  real(c_double), allocatable, target, save :: acc_depth(:)      ! (nd)
  real(c_double), allocatable, target, save :: acc_uwindn(:, :)  ! (2, n)
  real(c_double), allocatable, target, save :: acc_diffn(:)      ! (n)
  real(c_double), allocatable, target, save :: acc_mass(:, :)   ! (3*nd, 3*nd)
  real(c_double), allocatable, target, save :: acc_stiff(:, :)  ! (3*nd, 3*nd)
  real(c_double), allocatable, save :: acc_force(:)      ! (3*nd)
  real(c_double), allocatable, target, save :: acc_basis(:)      ! (nd)
  real(c_double), allocatable, target, save :: acc_dbasis(:, :)  ! (nd, 3)
contains

  ! -- entry state ----------------------------------------------------
  subroutine swa_begin(h_solver) bind(c, name="swa_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
    acc_prm => GetSolverParams(acc_solver)
  end subroutine swa_begin

  function swa_newton_iter() result(n) bind(c, name="swa_newton_iter")
    integer(c_int) :: n
    logical :: found
    n = GetInteger(acc_prm, 'Nonlinear System Newton After Iterations', found)
  end function swa_newton_iter

  function swa_newton_tol() result(x) bind(c, name="swa_newton_tol")
    real(c_double) :: x
    logical :: found
    x = GetConstReal(acc_prm, 'Nonlinear System Newton After Tolerance', found)
  end function swa_newton_tol

  function swa_nonlin_iter() result(n) bind(c, name="swa_nonlin_iter")
    integer(c_int) :: n
    logical :: found
    n = GetInteger(acc_prm, 'Nonlinear System Max Iterations', found)
    if (.not. found) n = 1
  end function swa_nonlin_iter

  function swa_nonlin_tol() result(x) bind(c, name="swa_nonlin_tol")
    real(c_double) :: x
    logical :: found
    x = GetConstReal(acc_prm, 'Nonlinear System Convergence Tolerance', found)
  end function swa_nonlin_tol

  function swa_massflux() result(r) bind(c, name="swa_massflux")
    integer(c_int) :: r
    logical :: found
    r = merge(1, 0, GetLogical(acc_prm, 'Mass Flux Discretization', found))
  end function swa_massflux

  ! -- assembly driver -------------------------------------------------
  function swa_nof_active() result(n) bind(c, name="swa_nof_active")
    integer(c_int) :: n
    n = GetNOFActive()
  end function swa_nof_active

  function swa_nof_boundary() result(n) bind(c, name="swa_nof_boundary")
    integer(c_int) :: n
    n = GetNOFBoundaryElements()
  end function swa_nof_boundary

  subroutine swa_default_initialize() bind(c, name="swa_default_initialize")
    call DefaultInitialize()
  end subroutine swa_default_initialize

  subroutine swa_default_finish_bulk() bind(c, name="swa_default_finish_bulk")
    call DefaultFinishBulkAssembly()
  end subroutine swa_default_finish_bulk

  subroutine swa_default_finish_bnd() bind(c, name="swa_default_finish_bnd")
    call DefaultFinishBoundaryAssembly()
  end subroutine swa_default_finish_bnd

  subroutine swa_default_finish() bind(c, name="swa_default_finish")
    call DefaultFinishAssembly()
  end subroutine swa_default_finish

  subroutine swa_default_dirichlet() bind(c, name="swa_default_dirichlet")
    call DefaultDirichletBCs()
  end subroutine swa_default_dirichlet

  function swa_default_solve() result(norm) bind(c, name="swa_default_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function swa_default_solve

  function swa_nonlin_change() result(r) bind(c, name="swa_nonlin_change")
    real(c_double) :: r
    r = acc_solver % variable % nonlinchange
  end function swa_nonlin_change

  subroutine swa_info(norm, rel) bind(c, name="swa_info")
    real(c_double), value :: norm, rel
    ! Exact transcription of the two WRITE/Info pairs (same statement,
    ! same compiler, same list-directed formatting as the original).
    write(message, *) 'Result Norm   : ', norm
    call info('ShallowWaterNSSolver', message, level=4)
    write(message, *) 'Relative Change : ', rel
    call info('ShallowWaterNSSolver', message, level=4)
  end subroutine swa_info

  ! -- bulk-element context --------------------------------------------
  subroutine swa_begin_active(t, n, nd) bind(c, name="swa_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n, nd
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes(acc_elem)
    nd = GetElementNOFDOFs(acc_elem)
  end subroutine swa_begin_active

  subroutine swa_local_ensure(n, nd) bind(c, name="swa_local_ensure")
    integer(c_int), value :: n, nd
    integer :: istat
    if (.not. allocated(acc_sol) .or. size(acc_sol, 2) < nd) then
      if (allocated(acc_sol)) deallocate(acc_sol)
      allocate(acc_sol(3, nd), stat=istat)
      if (istat /= 0) call fatal('ShallowWaterNSSolver', 'Memory allocation error.')
    end if
    if (.not. allocated(acc_depth) .or. size(acc_depth) < nd) then
      if (allocated(acc_depth)) deallocate(acc_depth)
      allocate(acc_depth(nd), stat=istat)
      if (istat /= 0) call fatal('ShallowWaterNSSolver', 'Memory allocation error.')
    end if
    if (.not. allocated(acc_uwindn) .or. size(acc_uwindn, 2) < n) then
      if (allocated(acc_uwindn)) deallocate(acc_uwindn)
      allocate(acc_uwindn(2, n), stat=istat)
      if (istat /= 0) call fatal('ShallowWaterNSSolver', 'Memory allocation error.')
    end if
    if (.not. allocated(acc_diffn) .or. size(acc_diffn) < n) then
      if (allocated(acc_diffn)) deallocate(acc_diffn)
      allocate(acc_diffn(n), stat=istat)
      if (istat /= 0) call fatal('ShallowWaterNSSolver', 'Memory allocation error.')
    end if
    if (.not. allocated(acc_mass) .or. size(acc_mass, 1) < 3 * nd) then
      if (allocated(acc_mass)) deallocate(acc_mass)
      if (allocated(acc_stiff)) deallocate(acc_stiff)
      if (allocated(acc_force)) deallocate(acc_force)
      allocate(acc_mass(3 * nd, 3 * nd), stat=istat)
      if (istat /= 0) call fatal('ShallowWaterNSSolver', 'Memory allocation error.')
      allocate(acc_stiff(3 * nd, 3 * nd), stat=istat)
      if (istat /= 0) call fatal('ShallowWaterNSSolver', 'Memory allocation error.')
      allocate(acc_force(3 * nd), stat=istat)
      if (istat /= 0) call fatal('ShallowWaterNSSolver', 'Memory allocation error.')
    end if
    if (.not. allocated(acc_basis) .or. size(acc_basis) < nd) then
      if (allocated(acc_basis)) deallocate(acc_basis)
      if (allocated(acc_dbasis)) deallocate(acc_dbasis)
      allocate(acc_basis(nd), stat=istat)
      if (istat /= 0) call fatal('ShallowWaterNSSolver', 'Memory allocation error.')
      allocate(acc_dbasis(nd, 3), stat=istat)
      if (istat /= 0) call fatal('ShallowWaterNSSolver', 'Memory allocation error.')
    end if
  end subroutine swa_local_ensure

  subroutine swa_local_fetch(n, nd, g, icz, wcoeff, fcor) &
      bind(c, name="swa_local_fetch")
    integer(c_int), value :: n, nd
    real(c_double) :: g, icz, wcoeff, fcor
    logical :: found
    call GetElementNodes(acc_nodes, acc_elem)
    acc_mat => GetMaterial(acc_elem)
    ! Mirrors the original's empty ASSOCIATED(BF) block: the lookup
    ! itself, kept as a touch call with no observable effect.
    acc_bf => GetBodyForce(acc_elem)
    if (associated(acc_bf)) then
    end if
    call GetVectorLocalSolution(acc_sol(1:3, 1:nd), UElement=acc_elem)
    call GetScalarLocalSolution(acc_depth(1:nd), 'Ground', UElement=acc_elem)
    acc_depth(1:nd) = -acc_depth(1:nd)
    g = GetCReal(acc_mat, 'Graviational Acceleration', found)
    icz = GetCReal(acc_mat, 'Inverse Chezy Coefficient', found)
    wcoeff = GetCReal(acc_mat, 'Wind-Stress Coefficient', found)
    acc_uwindn(1, 1:n) = GetReal(acc_mat, 'Wind 1', found, acc_elem)
    acc_uwindn(2, 1:n) = GetReal(acc_mat, 'Wind 2', found, acc_elem)
    fcor = GetCReal(acc_mat, 'Coriolis Coefficient', found)
    acc_diffn(1:n) = GetReal(acc_mat, 'Viscosity Coefficient', found, acc_elem)
  end subroutine swa_local_fetch

  ! -- LocalMatrix core (exact transcription) --------------------------
  subroutine swa_local_assemble(n, nd, newton, massflux, g, icz, wcoeff, fcor) &
      bind(c, name="swa_local_assemble")
    integer(c_int), value :: n, nd, newton, massflux
    real(c_double), value :: g, icz, wcoeff, fcor
    ! Mirrors: MASS(3*nd,3*nd), STIFF(3*nd,3*nd), FORCE(3*nd),
    ! SOL(3,nd), Depth(nd), Basis(nd), dBasisdx(nd,3), Uwind_n(2,n),
    ! diff_n(n).
    real(c_double), pointer :: m(:, :), a(:, :)
    real(c_double), pointer :: sol(:, :), depth(:), uwn(:, :), diffn(:)
    real(c_double), pointer :: basis(:), dbasis(:, :)
    real(c_double) :: s, detj, diff, speed, h, u(3), uwind(2)
    real(c_double) :: gradu(2, 2), gradf(2), gradh(2)
    logical :: stat
    integer :: i, j, k, p, q, t, dim, l
    sol => acc_sol(1:3, 1:nd)
    depth => acc_depth(1:nd)
    uwn => acc_uwindn(1:2, 1:n)
    diffn => acc_diffn(1:n)
    basis => acc_basis(1:nd)
    dbasis => acc_dbasis(1:nd, 1:3)
    dim = CoordinateSystemDimension()
    l = dim + 1
    acc_mass(1:3 * nd, 1:3 * nd) = 0.0_c_double
    acc_stiff(1:3 * nd, 1:3 * nd) = 0.0_c_double
    acc_force(1:3 * nd) = 0.0_c_double
    uwind(1) = 0.0_c_double
    uwind(2) = 0.0_c_double
    ! Numerical integration:
    ! ----------------------
    acc_ip = GaussPoints(acc_elem)
    do t = 1, acc_ip % n
      stat = ElementInfo(acc_elem, acc_nodes, acc_ip % u(t), acc_ip % v(t), &
          acc_ip % w(t), detj, basis, dbasis)
      s = acc_ip % s(t) * detj
      do i = 1, l
        u(i) = sum(sol(i, 1:nd) * basis(1:nd))
      end do
      gradu(1:dim, 1:dim) = matmul(sol(1:dim, 1:nd), dbasis(1:nd, 1:dim))
      do i = 1, dim
        gradf(i) = sum(sol(l, 1:nd) * dbasis(1:nd, i))
      end do
      speed = icz**2 * sqrt(sum(u(1:dim)**2))
      h = sum(depth(1:n) * basis(1:n))
      do i = 1, dim
        gradh(i) = sum(depth(1:n) * dbasis(1:n, i))
      end do
      do i = 1, dim
        uwind(i) = sum(basis(1:n) * uwn(i, 1:n))
      end do
      diff = sum(basis(1:n) * diffn(1:n))
      do p = 1, nd
        i = l * (p - 1)
        do q = 1, nd
          j = l * (q - 1)
          m => acc_mass(i + 1:i + l, j + 1:j + l)
          a => acc_stiff(i + 1:i + l, j + 1:j + l)
          do j = 1, l
            m(j, j) = m(j, j) + s * basis(q) * basis(p)
          end do
          do j = 1, dim
            a(j, j) = a(j, j) + s * diff * sum(dbasis(q, :) * dbasis(p, :))
            do k = 1, dim
              a(j, j) = a(j, j) + s * u(k) * dbasis(q, k) * basis(p)
              if (newton /= 0) then
                a(j, k) = a(j, k) + s * basis(q) * gradu(j, k) * basis(p)
              end if
            end do
            a(j, j) = a(j, j) + s * g * speed * basis(q) * basis(p)
            if (massflux /= 0) then
              a(j, l) = a(j, l) + s * g * dbasis(q, j) * basis(p)
            else
              a(j, l) = a(j, l) - s * g * basis(q) * dbasis(p, j)
            end if
            a(j, j) = a(j, j) + s * wcoeff * basis(q) / (h + u(l)) * basis(p)
            if (massflux /= 0) then
              a(l, j) = a(l, j) - s * basis(q) * (h + u(l)) * dbasis(p, j)
              if (newton /= 0) &
                  a(l, l) = a(l, l) - s * u(j) * basis(q) * dbasis(p, j)
            else
              a(l, j) = a(l, j) + s * (dbasis(q, j) * (h + u(l)) + &
                  basis(q) * (gradh(j) + gradf(j))) * basis(p)
              if (newton /= 0) &
                  a(l, l) = a(l, l) + s * (gradu(j, j) * basis(q) + &
                  u(j) * dbasis(q, j)) * basis(p)
            end if
          end do
          if (dim == 2) then
            a(1, 2) = a(1, 2) + s * fcor * basis(q) * basis(p)
            a(2, 1) = a(2, 1) + s * fcor * basis(q) * basis(p)
          end if
        end do
        do j = 1, dim
          if (newton /= 0) then
            if (massflux /= 0) then
              acc_force(i + l) = acc_force(i + l) - s * u(j) * u(l) * dbasis(p, j)
            else
              acc_force(i + l) = acc_force(i + l) + s * &
                  (u(l) * gradu(j, j) + u(j) * gradf(j)) * basis(p)
            end if
            do k = 1, dim
              acc_force(i + j) = acc_force(i + j) + s * u(k) * gradu(j, k) * basis(p)
            end do
          end if
        end do
        do j = 1, dim
          acc_force(i + j) = acc_force(i + j) + s * wcoeff * uwind(j) / &
              (h + u(l)) * basis(p)
        end do
      end do
    end do
    j = n + 1
    if (dim == 1) j = max(j, nd - 2)
    do i = j, nd
      q = l * (i - 1) + l
      acc_stiff(q, :) = 0.0_c_double
      acc_stiff(:, q) = 0.0_c_double
      acc_mass(q, :) = 0.0_c_double
      acc_mass(:, q) = 0.0_c_double
      acc_stiff(q, q) = 1.0_c_double
      acc_force(q) = 0.0_c_double
    end do
  end subroutine swa_local_assemble

  subroutine swa_local_update(nd, transient) bind(c, name="swa_local_update")
    integer(c_int), value :: nd, transient
    if (transient /= 0) then
      call Default1stOrderTime(acc_mass(1:3 * nd, 1:3 * nd), &
          acc_stiff(1:3 * nd, 1:3 * nd), acc_force(1:3 * nd), &
          UElement=acc_elem, USolver=acc_solver)
    end if
    call DefaultUpdateEquations(acc_stiff(1:3 * nd, 1:3 * nd), &
        acc_force(1:3 * nd), UElement=acc_elem, USolver=acc_solver)
  end subroutine swa_local_update

  ! -- boundary-element context -----------------------------------------
  subroutine swa_begin_boundary(t, n, nd, active) bind(c, name="swa_begin_boundary")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, active
    acc_elem => GetBoundaryElement(t)
    active = merge(1, 0, ActiveBoundaryElement(acc_elem))
    n = GetElementNOFNodes(acc_elem)
    nd = GetElementNOFDOFs(acc_elem)
  end subroutine swa_begin_boundary

  function swa_bc() result(r) bind(c, name="swa_bc")
    integer(c_int) :: r
    acc_bc => GetBC(acc_elem)
    r = merge(1, 0, associated(acc_bc))
  end function swa_bc

  ! -- BoundaryLocalMatrix core (exact transcription) -------------------
  subroutine swa_boundary_assemble(n, nd, transient) &
      bind(c, name="swa_boundary_assemble")
    integer(c_int), value :: n, nd, transient
    real(c_double), pointer :: a(:, :)
    real(c_double), pointer :: sol(:, :), depth(:), basis(:), dbasis(:, :)
    type(element_t), pointer :: parent
    real(c_double) :: s, detj
    real(c_double) :: u(3), udot, h, nrm(3)
    logical :: stat, rout
    integer :: i, j, p, q, t, pind(2), dim, l
    sol => acc_sol(1:3, 1:nd)
    depth => acc_depth(1:nd)
    basis => acc_basis(1:nd)
    dbasis => acc_dbasis(1:nd, 1:3)
    dim = CoordinateSystemDimension()
    l = dim + 1
    call GetElementNodes(acc_nodes, acc_elem)
    rout = .false.
    if (dim == 1) then
      parent => acc_elem % BoundaryInfo % Left
      call GetElementNodes(acc_pnodes, parent)
      pind = [1, 2]
      if (acc_elem % NodeIndexes(1) /= parent % NodeIndexes(1)) pind = [2, 1]
      rout = acc_pnodes % x(pind(1)) > acc_pnodes % x(pind(2))
    end if
    acc_mass(1:3 * nd, 1:3 * nd) = 0.0_c_double
    acc_stiff(1:3 * nd, 1:3 * nd) = 0.0_c_double
    acc_force(1:3 * nd) = 0.0_c_double
    call GetVectorLocalSolution(sol, UElement=acc_elem)
    call GetScalarLocalSolution(depth, 'Ground', UElement=acc_elem)
    depth(1:nd) = -depth(1:nd)
    ! Numerical integration:
    ! ----------------------
    acc_ip = GaussPoints(acc_elem)
    do t = 1, acc_ip % n
      stat = ElementInfo(acc_elem, acc_nodes, acc_ip % u(t), acc_ip % v(t), &
          acc_ip % w(t), detj, basis, dbasis)
      s = acc_ip % s(t) * detj
      do i = 1, dim + 1
        u(i) = sum(sol(i, 1:nd) * basis(1:nd))
      end do
      if (dim > 1) then
        nrm = NormalVector(acc_elem, acc_nodes, acc_ip % u(t), acc_ip % v(t), .true.)
        udot = sum(nrm(1:dim) * u(1:dim))
      else if (rout) then
        udot = u(1)
      else
        udot = -u(1)
      end if
      h = sum(basis(1:n) * depth(1:n))
      do p = 1, nd
        i = l * (p - 1)
        do q = 1, nd
          j = l * (q - 1)
          a => acc_stiff(i + 1:i + l, j + 1:j + l)
          a(l, l) = a(l, l) + s * udot * basis(q) * basis(p)
        end do
        acc_force(i + l) = acc_force(i + l) - s * udot * h * basis(p)
      end do
    end do
    if (transient /= 0) then
      call Default1stOrderTime(acc_mass(1:3 * nd, 1:3 * nd), &
          acc_stiff(1:3 * nd, 1:3 * nd), acc_force(1:3 * nd), &
          UElement=acc_elem)
    end if
    call DefaultUpdateEquations(acc_stiff(1:3 * nd, 1:3 * nd), &
        acc_force(1:3 * nd), UElement=acc_elem)
  end subroutine swa_boundary_assemble
end module shallowwaterns_acc
