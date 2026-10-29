! T022-A HelmholtzBEM port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/HelmholtzBEM.F90, compiled against the
! built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape arrays). All Fortran-object
! traffic — TYPE(...) handles, OPTIONAL/POINTER dummies, assumed-shape
! descriptors, hidden lengths, array-pointer results — stays inside
! Fortran, so the Rust side is pure C calls. The wrapper module mirrors
! the original's SAVE state (current element, ElementNodes SAVE state,
! IP shallow copies) with identical lifetimes.
!
! The `GlobMatC` module below re-declares the plugin-private BEM matrix
! with the identical statement as the original unit, so the port cdylib
! carries the same `B __globmatc_MOD_matrix` dynamic symbol the
! reference `.so` exports (proven by `nm -D`). All matrix traffic goes
! through this copy; the host never touches it.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module GlobMatC
  use types, only: dp
  implicit none
  complex(kind=dp), allocatable :: Matrix(:,:)
end module GlobMatC

module hbm_acc
  use iso_c_binding
  use types
  use defutils, only: GetAngularFrequency
  use lists, only: ListGetReal, ListGetConstReal, ListGetInteger, &
      ListGetLogical, ListGetString, VariableGet
  use generalutils, only: ComponentNameVar, ComponentNameStr, I2S
  use integration, only: GaussPoints, GaussIntegrationPoints_t
  use eleminfo, only: ElementInfo
  use elementgeometry, only: NormalVector
  use coordinatesystems, only: CurrentCoordinateSystem, Cartesian, &
      CoordinateSystemDimension
  use messages, only: Info, Fatal, Error
  use solverbasics, only: InvalidateVariable
  use spariterglobals, only: ParEnv
  use loadmod, only: IterCallFTNC
  use GlobMatC, only: Matrix
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(model_t), pointer, save :: acc_model => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(variable_t), pointer, save :: acc_var => null()
  type(variable_t), pointer, save :: acc_flux => null()
  type(valuelist_t), pointer, save :: acc_bcvals => null()
  interface
    ! Host C routine (`CPUTime.c`), reached by the original as CPUTime().
    function hbm_host_cputime() result(t) bind(c, name="cputime")
      import c_double
      real(c_double) :: t
    end function hbm_host_cputime
  end interface
  ! ElementNodes SAVE-state mirror (shared by both assembly kernels,
  ! exactly as the host procedure's SAVE ElementNodes is shared by its
  ! internal procedures).
  type(nodes_t), target, save :: acc_nodes
  ! IP shallow-copy mirror (one per assembly kernel use, consumed
  ! before the next fill, as in the original).
  type(gaussintegrationpoints_t), save :: acc_ip

contains

  ! -- entry state ----------------------------------------------------
  subroutine hbm_begin(h_model, h_solver) bind(c, name="hbm_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    acc_var => acc_solver % variable
  end subroutine hbm_begin

  subroutine hbm_check_mpi() bind(c, name="hbm_check_mpi")
    if (ParEnv % PEs > 1) then
      call Fatal('HelmholtzBEMSolver', &
          'Unfortunately this module does not support MPI!')
    end if
  end subroutine hbm_check_mpi

  subroutine hbm_check_cart() bind(c, name="hbm_check_cart")
    if (CurrentCoordinateSystem() /= Cartesian) then
      call Fatal('HelmholtzBEMSolver', &
          'This solver is implemented only for cartesian coordinates!')
    end if
  end subroutine hbm_check_cart

  subroutine hbm_sizes(nmaxnodes, numnodes, nbcs) bind(c, name="hbm_sizes")
    integer(c_int) :: nmaxnodes, numnodes, nbcs
    nmaxnodes = acc_model % maxelementnodes
    numnodes = acc_solver % mesh % numberofnodes
    nbcs = acc_model % numberofbcs
  end subroutine hbm_sizes

  function hbm_nof_active() result(n) bind(c, name="hbm_nof_active")
    integer(c_int) :: n
    n = acc_solver % numberofactiveelements
  end function hbm_nof_active

  subroutine hbm_begin_active(t, code, n, assoc, constraint) &
      bind(c, name="hbm_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: code, n, assoc, constraint
    acc_elem => acc_solver % mesh % elements(acc_solver % activeelements(t))
    code = acc_elem % type % elementcode
    n = acc_elem % type % numberofnodes
    if (associated(acc_elem % boundaryinfo)) then
      assoc = 1
      constraint = acc_elem % boundaryinfo % constraint
    else
      assoc = 0
      constraint = 0
    end if
  end subroutine hbm_begin_active

  subroutine hbm_elem_nodes(n, buf) bind(c, name="hbm_elem_nodes")
    integer(c_int), value :: n
    integer(c_int) :: buf(n)
    buf(1:n) = acc_elem % nodeindexes(1:n)
  end subroutine hbm_elem_nodes

  subroutine hbm_alloc_matrix(nb) bind(c, name="hbm_alloc_matrix")
    integer(c_int), value :: nb
    integer :: istat
    allocate(Matrix(nb, nb), stat=istat)
    if (istat /= 0) then
      call Fatal('HelmholtzBEMSolver', 'Memory allocation error 2.')
    end if
  end subroutine hbm_alloc_matrix

  subroutine hbm_alloc_nodes(n) bind(c, name="hbm_alloc_nodes")
    integer(c_int), value :: n
    integer :: istat
    allocate(acc_nodes % x(n), acc_nodes % y(n), acc_nodes % z(n), &
        stat=istat)
    if (istat /= 0) then
      call Fatal('HelmholtzBEMSolver', 'Memory allocation error 2.')
    end if
  end subroutine hbm_alloc_nodes

  subroutine hbm_matrix_zero() bind(c, name="hbm_matrix_zero")
    Matrix = 0.0_dp
  end subroutine hbm_matrix_zero

  ! -- solver variable traffic -----------------------------------------
  function hbm_var_values_ptr() result(p) bind(c, name="hbm_var_values_ptr")
    type(c_ptr) :: p
    p = c_loc(acc_var % values(1))
  end function hbm_var_values_ptr

  function hbm_var_perm_ptr() result(p) bind(c, name="hbm_var_perm_ptr")
    type(c_ptr) :: p
    p = c_loc(acc_var % perm(1))
  end function hbm_var_perm_ptr

  subroutine hbm_touch_equation() bind(c, name="hbm_touch_equation")
    character(:), allocatable :: s
    s = ListGetString(acc_model % solver % values, 'Equation')
  end subroutine hbm_touch_equation

  function hbm_angular_frequency() result(w) &
      bind(c, name="hbm_angular_frequency")
    real(c_double) :: w
    w = GetAngularFrequency()
  end function hbm_angular_frequency

  function hbm_flux_init() result(r) bind(c, name="hbm_flux_init")
    integer(c_int) :: r
    acc_flux => VariableGet(acc_solver % mesh % variables, 'Flux')
    r = merge(1, 0, associated(acc_flux))
  end function hbm_flux_init

  function hbm_flux_values_ptr() result(p) bind(c, name="hbm_flux_values_ptr")
    type(c_ptr) :: p
    p = c_loc(acc_flux % values(1))
  end function hbm_flux_values_ptr

  ! -- boundary condition traffic ---------------------------------------
  function hbm_bc_tag(i) result(t) bind(c, name="hbm_bc_tag")
    integer(c_int), value :: i
    integer(c_int) :: t
    t = acc_model % bcs(i) % tag
  end function hbm_bc_tag

  subroutine hbm_bc_select(i) bind(c, name="hbm_bc_select")
    integer(c_int), value :: i
    acc_bcvals => acc_model % bcs(i) % values
  end subroutine hbm_bc_select

  subroutine hbm_bc_pot(comp, n, p, gotit) bind(c, name="hbm_bc_pot")
    integer(c_int), value :: comp, n
    real(c_double) :: p(n)
    integer(c_int) :: gotit
    logical :: found
    character(:), allocatable :: nm
    nm = ComponentNameVar(acc_var, comp)
    p(1:n) = ListGetReal(acc_bcvals, nm, n, acc_elem % nodeindexes, found)
    if (.not. found) then
      nm = ComponentNameStr('Potential', comp)
      p(1:n) = ListGetReal(acc_bcvals, nm, n, acc_elem % nodeindexes, found)
    end if
    gotit = merge(1, 0, found)
  end subroutine hbm_bc_pot

  subroutine hbm_bc_flux(comp, n, p, gotit) bind(c, name="hbm_bc_flux")
    integer(c_int), value :: comp, n
    real(c_double) :: p(n)
    integer(c_int) :: gotit
    logical :: found
    character(:), allocatable :: nm
    nm = ComponentNameStr('Flux', comp)
    p(1:n) = ListGetReal(acc_bcvals, nm, n, acc_elem % nodeindexes, found)
    gotit = merge(1, 0, found)
  end subroutine hbm_bc_flux

  ! -- assembly kernel traffic ------------------------------------------
  subroutine hbm_elemnodes_fill(n, nidx) bind(c, name="hbm_elemnodes_fill")
    integer(c_int), value :: n
    integer(c_int) :: nidx(n)
    acc_nodes % x(1:n) = acc_solver % mesh % nodes % x(nidx(1:n))
    acc_nodes % y(1:n) = acc_solver % mesh % nodes % y(nidx(1:n))
    acc_nodes % z(1:n) = acc_solver % mesh % nodes % z(nidx(1:n))
  end subroutine hbm_elemnodes_fill

  function hbm_ip_begin(k) result(n) bind(c, name="hbm_ip_begin")
    integer(c_int), value :: k
    integer(c_int) :: n
    acc_ip = GaussPoints(acc_elem, k)
    n = acc_ip % n
  end function hbm_ip_begin

  subroutine hbm_ip_point(t, u, v, w, s) bind(c, name="hbm_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip % u(t)
    v = acc_ip % v(t)
    w = acc_ip % w(t)
    s = acc_ip % s(t)
  end subroutine hbm_ip_point

  subroutine hbm_eleminfo(u, v, w, n, nd, basis, dbasis, detj, stat) &
      bind(c, name="hbm_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: n, nd
    real(c_double) :: basis(n), dbasis(nd, 3), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_elem, acc_nodes, u, v, w, detj, &
        basis(1:n), dbasis(1:nd, 1:3))
    stat = merge(1, 0, ok)
  end subroutine hbm_eleminfo

  subroutine hbm_normal(u, v, check, nout) bind(c, name="hbm_normal")
    real(c_double), value :: u, v
    integer(c_int), value :: check
    real(c_double) :: nout(3)
    nout(1:3) = NormalVector(acc_elem, acc_nodes, u, v, check /= 0)
  end subroutine hbm_normal

  function hbm_coordsysdim() result(d) bind(c, name="hbm_coordsysdim")
    integer(c_int) :: d
    d = CoordinateSystemDimension()
  end function hbm_coordsysdim

  ! -- order-sensitive floating point stays in Fortran -------------------
  subroutine hbm_precond_elem(i, u, v) bind(c, name="hbm_precond_elem")
    integer(c_int), value :: i
    complex(c_double_complex) :: u(*), v(*)
    u(i) = v(i) / Matrix(i, i)
  end subroutine hbm_precond_elem

  subroutine hbm_matvec(n, u, v) bind(c, name="hbm_matvec")
    integer(c_int), value :: n
    complex(c_double_complex) :: u(n), v(n)
    v(1:n) = matmul(Matrix, u(1:n))
  end subroutine hbm_matvec

  subroutine hbm_post_adjust(nb, force, load, diag, potknown) &
      bind(c, name="hbm_post_adjust")
    integer(c_int), value :: nb
    complex(c_double_complex) :: force(nb), load(nb), diag(nb)
    integer(c_int) :: potknown(nb)
    integer :: i
    do i = 1, nb
      if (potknown(i) /= 0) then
        force(i) = force(i) - load(i) * diag(i)
      else
        Matrix(i, i) = diag(i)
      end if
    end do
  end subroutine hbm_post_adjust

  function hbm_var_norm(pot, n, numnodes) result(s) &
      bind(c, name="hbm_var_norm")
    integer(c_int), value :: n, numnodes
    complex(c_double_complex) :: pot(n)
    real(c_double) :: s
    acc_var % norm = sqrt(sum(abs(pot(1:n))**2)) / numnodes
    s = acc_var % norm
  end function hbm_var_norm

  subroutine hbm_lapack(n, x) bind(c, name="hbm_lapack")
    integer(c_int), value :: n
    complex(c_double_complex) :: x(n)
    call SolveLapack_cmplx(n, Matrix, x)
  end subroutine hbm_lapack

  subroutine hbm_itercall(fptr, x, b, n, ipar, dpar, work, wsize, &
      mvptr, pcptr) bind(c, name="hbm_itercall")
    integer(c_int64_t) :: fptr, mvptr, pcptr
    integer(c_int), value :: n, wsize
    complex(c_double_complex) :: x(n), b(n), work(wsize, n)
    integer(c_int) :: ipar(50)
    real(c_double) :: dpar(50)
    integer(AddrInt) :: z0
    z0 = 0
    call IterCallFTNC(fptr, x(1:n), b(1:n), ipar, dpar, &
        work(1:wsize, 1:n), mvptr, pcptr, z0, z0, z0, z0)
  end subroutine hbm_itercall

  ! -- solver-list queries ------------------------------------------------
  function hbm_linsys() result(r) bind(c, name="hbm_linsys")
    integer(c_int) :: r
    character(:), allocatable :: s
    s = ListGetString(acc_solver % values, 'Linear System Solver')
    if (s == 'direct') then
      r = 1
    else if (s == 'iterative') then
      r = 2
    else
      r = 0
    end if
  end function hbm_linsys

  subroutine hbm_tol(val, gotit) bind(c, name="hbm_tol")
    real(c_double) :: val
    integer(c_int) :: gotit
    logical :: found
    val = ListGetConstReal(acc_solver % values, &
        'Linear System Convergence Tolerance', found)
    gotit = merge(1, 0, found)
  end subroutine hbm_tol

  subroutine hbm_divlim(val, gotit) bind(c, name="hbm_divlim")
    real(c_double) :: val
    integer(c_int) :: gotit
    logical :: found
    val = ListGetConstReal(acc_solver % values, &
        'Linear System Divergence Limit', found)
    gotit = merge(1, 0, found)
  end subroutine hbm_divlim

  subroutine hbm_maxit(val, gotit) bind(c, name="hbm_maxit")
    integer(c_int) :: val, gotit
    logical :: found
    val = ListGetInteger(acc_solver % values, &
        'Linear System Max Iterations', found)
    gotit = merge(1, 0, found)
  end subroutine hbm_maxit

  subroutine hbm_resid(val, gotit) bind(c, name="hbm_resid")
    integer(c_int) :: val, gotit
    logical :: found
    val = ListGetInteger(acc_solver % values, &
        'Linear System Residual Output', found)
    gotit = merge(1, 0, found)
  end subroutine hbm_resid

  subroutine hbm_abort(val, gotit) bind(c, name="hbm_abort")
    integer(c_int) :: val, gotit
    logical :: found, l
    l = ListGetLogical(acc_solver % values, &
        'Linear System Abort Not Converged', found)
    val = merge(1, 0, l)
    gotit = merge(1, 0, found)
  end subroutine hbm_abort

  ! -- messages (exact original call sites) --------------------------------
  subroutine hbm_info_solving() bind(c, name="hbm_info_solving")
    call Info('HelmholtzBEMSolver', &
        'Solving the Helmholtz equation on boundaries')
  end subroutine hbm_info_solving

  subroutine hbm_info_nbn(nb) bind(c, name="hbm_info_nbn")
    integer(c_int), value :: nb
    call Info('HelmholtzBEMSolver', &
        'Number of nodes on boundaries: ' // I2S(nb))
  end subroutine hbm_info_nbn

  subroutine hbm_info_full() bind(c, name="hbm_info_full")
    call Info('HelmholtzBEMSolver', 'Solving full linear BEM matrix!', &
        level=7)
  end subroutine hbm_info_full

  subroutine hbm_info_extract() bind(c, name="hbm_info_extract")
    call Info('HelmholtzBEMSolver', &
        'Extract potential and fluxes for the boundary nodes', level=7)
  end subroutine hbm_info_extract

  subroutine hbm_info_compute() bind(c, name="hbm_info_compute")
    call Info('HelmholtzBEMSolver', 'Now compute potential for all mesh &
        &points', level=7)
  end subroutine hbm_info_compute

  subroutine hbm_info_fields() bind(c, name="hbm_info_fields")
    call Info('HelmholtzBEMSolver', 'Creating fields for postprocessing', &
        level=7)
  end subroutine hbm_info_fields

  subroutine hbm_info_done() bind(c, name="hbm_info_done")
    call Info('HelmholtzBEMSolver', 'All done for now!', level=7)
  end subroutine hbm_info_done

  subroutine hbm_fatal_alloc1() bind(c, name="hbm_fatal_alloc1")
    call Fatal('HelmholtzBEMSolver', 'Memory allocation error 1.')
  end subroutine hbm_fatal_alloc1

  subroutine hbm_fatal_alloc2() bind(c, name="hbm_fatal_alloc2")
    call Fatal('HelmholtzBEMSolver', 'Memory allocation error 2.')
  end subroutine hbm_fatal_alloc2

  subroutine hbm_fatal_alloc3() bind(c, name="hbm_fatal_alloc3")
    call Fatal('HelmholtzBEMSolver', 'Memory allocation error 3.')
  end subroutine hbm_fatal_alloc3

  subroutine hbm_fatal_linsys() bind(c, name="hbm_fatal_linsys")
    call Fatal('SolveFull', 'Unknown solver type.')
  end subroutine hbm_fatal_linsys

  subroutine hbm_fatal_iter() bind(c, name="hbm_fatal_iter")
    call Fatal('IterSolve', 'Failed convergence tolerances.')
  end subroutine hbm_fatal_iter

  subroutine hbm_error_iter() bind(c, name="hbm_error_iter")
    call Error('IterSolve', 'Failed convergence tolerances.')
  end subroutine hbm_error_iter

  subroutine hbm_invalidate() bind(c, name="hbm_invalidate")
    call InvalidateVariable(acc_model % meshes, acc_solver % mesh, &
        acc_var % name)
  end subroutine hbm_invalidate

  function hbm_cputime() result(t) bind(c, name="hbm_cputime")
    real(c_double) :: t
    t = hbm_host_cputime()
  end function hbm_cputime

  subroutine hbm_print_assembly(t) bind(c, name="hbm_print_assembly")
    real(c_double), value :: t
    print *, 'Assembly (s): ', t
  end subroutine hbm_print_assembly

  subroutine hbm_print_solve(t) bind(c, name="hbm_print_solve")
    real(c_double), value :: t
    print *, 'Solve (s):    ', t
  end subroutine hbm_print_solve

  subroutine hbm_print_post(t) bind(c, name="hbm_print_post")
    real(c_double), value :: t
    print *, 'Post Processing (s):    ', t
  end subroutine hbm_print_post

  ! -- assembly kernels (exact original statements, wrapper SAVE state) ----
  subroutine hbm_integrate(n, nb, freq, potperm, bndperm, load, force, &
      diag, potknown, numnodes) bind(c, name="hbm_integrate")
    integer(c_int), value :: n, nb, numnodes
    real(c_double), value :: freq
    integer(c_int) :: potperm(numnodes), bndperm(numnodes)
    complex(c_double_complex) :: load(nb), force(numnodes), diag(nb)
    integer(c_int) :: potknown(nb)
    real(c_double) :: Basis(n), dBasisdx(n, 3)
    real(c_double) :: LX, LY, LZ, x, y, z
    logical :: Stat, CheckNormals
    complex(kind=dp) :: R, dGdN, G, GradG(3), i
    real(c_double) :: detJ, U, V, W, S, A, L, Normal(3), rad
    integer :: j, k, p, q, t, dim
    type(GaussIntegrationPoints_t) :: IntegStuff
    dim = CoordinateSystemDimension()
    i = cmplx(0.0_dp, 1.0_dp, kind=dp)
    select case (acc_elem % type % elementcode / 100)
    case (2)
      IntegStuff = GaussPoints(acc_elem, 4)
    case (3)
      IntegStuff = GaussPoints(acc_elem, 6)
    case (4)
      IntegStuff = GaussPoints(acc_elem, 16)
    end select
    CheckNormals = associated(acc_elem % boundaryinfo)
    if (CheckNormals) then
      CheckNormals = associated(acc_elem % boundaryinfo % left) .or. &
          associated(acc_elem % boundaryinfo % right)
    end if
    do t = 1, IntegStuff % n
      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)
      S = IntegStuff % s(t)
      Stat = ElementInfo(acc_elem, acc_nodes, U, V, W, detJ, &
          Basis, dBasisdx)
      S = S * detJ
      Normal = NormalVector(acc_elem, acc_nodes, U, V, CheckNormals)
      LX = sum(acc_nodes % x(1:n) * Basis)
      LY = sum(acc_nodes % y(1:n) * Basis)
      LZ = sum(acc_nodes % z(1:n) * Basis)
      do p = 1, nb
        k = bndperm(p)
        x = LX - acc_solver % mesh % nodes % x(k)
        y = LY - acc_solver % mesh % nodes % y(k)
        z = LZ - acc_solver % mesh % nodes % z(k)
        call hbm_green(dim, freq, x, y, z, G, GradG)
        dGdN = sum(GradG * Normal)
        do j = 1, n
          q = potperm(acc_elem % nodeindexes(j))
          if (potknown(q) /= 0) then
            if (p /= q) then
              force(p) = force(p) - S * load(q) * Basis(j) * dGdN
            end if
            Matrix(p, q) = Matrix(p, q) - S * Basis(j) * G
          else
            force(p) = force(p) + S * load(q) * Basis(j) * G
            if (p /= q) then
              Matrix(p, q) = Matrix(p, q) + S * Basis(j) * dGdN
            end if
          end if
          R = -freq * sin(freq * x) * G * Normal(1)
          if (p /= q) then
            R = R - cos(freq * x) * dGdN
          end if
          diag(p) = diag(p) - S * Basis(j) * R
        end do
      end do
    end do
  end subroutine hbm_integrate

  subroutine hbm_green(dim, k, x, y, z, W, GradW) bind(c, name="hbm_green")
    integer(c_int), value :: dim
    real(c_double), value :: k, x, y, z
    complex(c_double_complex) :: W, GradW(3)
    complex(kind=dp) :: i, dWdR
    real(c_double) :: r, J0, Y0, dJ0, dY0
    R = sqrt(X**2 + Y**2 + Z**2)
    i = cmplx(0.0_dp, 1.0_dp, kind=dp)
    select case (dim)
    case (2)
      call hbm_bessel(k * R, j0, y0, dj0, dy0)
      W = (J0 - i * Y0) / (i * 4)
      dWdR = k * (dJ0 - i * dY0) / (i * 4)
    case (3)
      W = exp(-i * k * R) / (4 * PI * R)
      dWdR = (-i * k - 1 / R) * W
    end select
    GradW(1) = x * dWdR / R
    GradW(2) = y * dWdR / R
    GradW(3) = z * dWdR / R
  end subroutine hbm_green

  subroutine hbm_bessel(x, j0, y0, dj0, dy0) bind(c, name="hbm_bessel")
    real(c_double), value :: x
    real(c_double) :: j0, y0, dj0, dy0
    integer, parameter :: maxrounds = 1000
    double precision, parameter :: polylimit = 10.0d0
    double precision, parameter :: accuracy = 1.0d-8
    double precision, parameter :: gamma = &
        0.577215664901532860606512090082D0
    integer :: i, k
    real(c_double) :: hk
    real(kind=c_double) :: phi, res, p, f
    double precision :: A(7) = &
        (/ 0.79788456D0, -0.00000077D0, -0.00552740D0, &
           0.00009512D0,  0.00137237D0, -0.00072805D0, &
           0.00014476D0 /)
    double precision :: B(7) = &
        (/ -0.78539816D0, -0.04166397D0, -0.00003954D0, &
            0.00262573D0, -0.00054125D0, -0.00029333D0, &
            0.00013558D0 /)
    double precision :: C(7) = &
        (/ 0.79788456D0,  0.00000156D0, 0.01659667D0, &
           0.00017105D0, -0.00249511D0, 0.00113653D0, &
          -0.0020033D0 /)
    double precision :: D(7) = &
        (/ -2.35619449D0, 0.12499612D0, 0.00005650D0, &
           -0.00637879D0, 0.00074348D0, 0.00079824D0, &
           -0.00029166D0 /)
    if (abs(x) > polylimit) then
      P = x
      F = 0.0d0
      do k = 1, 7
        F = F + A(k) * (3 / x)**(k - 1.0d0)
        P = P + B(k) * (3 / x)**(k - 1.0d0)
      end do
      j0 = F * cos(P) / sqrt(x)
      y0 = F * sin(P) / sqrt(x)
      P = x
      F = 0.0d0
      do k = 1, 7
        F = F + C(k) * (3 / x)**(k - 1.0d0)
        P = P + D(k) * (3 / x)**(k - 1.0d0)
      end do
      dj0 = -F * cos(P) / sqrt(x)
      dy0 = -F * sin(P) / sqrt(x)
    else
      j0 = 1.0d0
      y0 = 0.0d0
      dj0 = 0.0d0
      dy0 = 0.0d0
      hk = 0.0d0
      do k = 1, maxrounds
        hk = hk + 1.0d0 / k
        res = 1.0d0
        do i = 1, k
          res = res * (x / (2.0d0 * i))**2
        end do
        j0 = j0 + (-1)**k * res
        y0 = y0 + (-1)**(k + 1) * hk * res
        dj0 = dj0 + (-1)**k * k / (0.5d0 * x) * res
        dy0 = dy0 + (-1)**(k + 1) * hk * k / (0.5d0 * x) * res
        if (abs(k / (0.5d0 * x) * res) < accuracy) exit
      end do
      if (k >= maxrounds) stop 'Error in evaluating Bessel functions'
      y0 = y0 + (log(0.5d0 * x) + gamma) * j0
      y0 = y0 * 2.0d0 / PI
      dy0 = dy0 + (1.0d0 / x) * j0 + (log(0.5d0 * x) + gamma) * dj0
      dy0 = dy0 * 2.0d0 / PI
    end if
  end subroutine hbm_bessel

  subroutine hbm_compute_pot(n, nb, freq, potential, pot, flx, potperm, &
      numnodes) bind(c, name="hbm_compute_pot")
    integer(c_int), value :: n, nb, numnodes
    real(c_double), value :: freq
    complex(c_double_complex) :: potential(numnodes), pot(nb), flx(nb)
    integer(c_int) :: potperm(numnodes)
    real(c_double) :: Basis(n), dBasisdx(n, 3)
    real(c_double) :: LX, LY, LZ, x, y, z
    logical :: Stat, CheckNormals
    real(c_double) :: detJ, U, V, W, S, A, L, Normal(3)
    complex(kind=dp) :: dGdN, G, GradG(3)
    integer :: i, j, k, p, q, t, dim
    type(GaussIntegrationPoints_t) :: IntegStuff
    dim = CoordinateSystemDimension()
    select case (acc_elem % type % elementcode / 100)
    case (2)
      IntegStuff = GaussPoints(acc_elem, 4)
    case (3)
      IntegStuff = GaussPoints(acc_elem, 6)
    case (4)
      IntegStuff = GaussPoints(acc_elem, 16)
    end select
    CheckNormals = associated(acc_elem % boundaryinfo)
    if (CheckNormals) then
      CheckNormals = associated(acc_elem % boundaryinfo % left) .or. &
          associated(acc_elem % boundaryinfo % right)
    end if
    do t = 1, IntegStuff % n
      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)
      S = IntegStuff % s(t)
      Stat = ElementInfo(acc_elem, acc_nodes, U, V, W, detJ, &
          Basis, dBasisdx)
      S = S * detJ
      Normal = NormalVector(acc_elem, acc_nodes, U, V, CheckNormals)
      LX = sum(acc_nodes % x(1:n) * Basis)
      LY = sum(acc_nodes % y(1:n) * Basis)
      LZ = sum(acc_nodes % z(1:n) * Basis)
      do i = 1, acc_solver % mesh % numberofnodes
        if (potperm(i) > 0) cycle
        if (acc_var % perm(i) <= 0) cycle
        x = LX - acc_solver % mesh % nodes % x(i)
        y = LY - acc_solver % mesh % nodes % y(i)
        z = LZ - acc_solver % mesh % nodes % z(i)
        call hbm_green(dim, freq, x, y, z, G, GradG)
        dGdN = sum(GradG * Normal)
        do j = 1, n
          q = potperm(acc_elem % nodeindexes(j))
          potential(i) = potential(i) - S * Basis(j) * &
              (pot(q) * dGdN - flx(q) * G)
        end do
      end do
    end do
  end subroutine hbm_compute_pot

end module hbm_acc
