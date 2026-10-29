! T021-B ScalarPotentialSolver port: statement-level BIND(C) wrappers for
! the Elmer calls in fem/src/modules/ScalarPotentialSolver.F90, compiled
! against the built .mod files of the target configure.
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module mirrors the original's locals (current element, Nodes SAVE state
! of the BulkAssembly internal procedure, the Conductivity SAVE pointer,
! and the fetched solver-parameter snapshot) with identical lifetimes.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module scalarpotentialsolver_acc
  use iso_c_binding
  use types
  use defutils
  use coordinatesystems
  use integration
  use matrixassembly
  use eleminfo
  use messages
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(element_t), pointer, save :: acc_elem => null()
  ! BulkAssembly context (mirrors its SAVE Nodes and SAVE Conductivity).
  type(nodes_t), target, save :: acc_nodes
  real(kind=dp), pointer, save :: acc_cond(:,:,:) => null()
  type(gaussintegrationpoints_t), save :: acc_ip
  ! Solver-parameter snapshot taken by sps_fetch, mirroring the entry's
  ! locals (VarName, CondName, ConstantBulkMatrix, CSymmetry, dim).
  character(len=MAX_NAME_LEN), save :: acc_varname = ' '
  character(len=MAX_NAME_LEN), save :: acc_condname = ' '
  logical, save :: acc_constbulk = .false.
  logical, save :: acc_inuse = .false.
  integer, save :: acc_dim = 0
  logical, save :: acc_csym = .false.

contains

  ! -- init entry -------------------------------------------------------
  subroutine sps_init() bind(c, name="sps_init")
    type(valuelist_t), pointer :: params
    integer :: dim
    params => GetSolverParams()
    dim = CoordinateSystemDimension()
    call ListAddNewString(params, 'Variable', 'Scalar Potential')
    call ListAddInteger(params, 'Time derivative order', 0)
    if (.not. ListCheckPresent(params, 'Linear System Solver')) then
      call Info('ScalarPotentialSolver_init', &
          'Setting defaults for linear system solver')
      call ListAddNewString(params, 'Linear System Solver', 'Iterative')
      call ListAddNewString(params, 'Linear System Iterative Method', 'cg')
      call ListAddNewString(params, 'Linear System Preconditioning', 'ILU0')
      call ListAddNewInteger(params, 'Linear System Max Iterations', 500)
      call ListAddNewInteger(params, 'Linear System Residual Output', 10)
      call ListAddNewConstReal(params, &
          'Linear System Convergence Tolerance', 1.0e-10_dp)
    end if
  end subroutine sps_init

  ! -- entry state ------------------------------------------------------
  subroutine sps_begin(h_solver) bind(c, name="sps_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine sps_begin

  function sps_assoc_matrix() result(r) bind(c, name="sps_assoc_matrix")
    integer(c_int) :: r
    r = merge(1, 0, associated(acc_solver % matrix))
  end function sps_assoc_matrix

  function sps_active_dofs() result(r) bind(c, name="sps_active_dofs")
    integer(c_int) :: r
    r = merge(1, 0, count(acc_solver % variable % perm > 0) > 0)
  end function sps_active_dofs

  subroutine sps_fetch() bind(c, name="sps_fetch")
    type(valuelist_t), pointer :: params
    logical :: got
    params => GetSolverParams()
    acc_dim = CoordinateSystemDimension()
    acc_varname = GetString(params, 'Flux Variable')
    acc_condname = ListGetString(params, 'Flux Coefficient', got)
    if (.not. got) acc_condname = 'none'
    acc_constbulk = GetLogical(params, 'Constant Bulk Matrix', got)
    acc_inuse = acc_constbulk .and. &
        associated(acc_solver % matrix % bulkvalues)
    acc_csym = CurrentCoordinateSystem() == AxisSymmetric .or. &
        CurrentCoordinateSystem() == CylindricSymmetric
  end subroutine sps_fetch

  function sps_dim() result(n) bind(c, name="sps_dim")
    integer(c_int) :: n
    n = acc_dim
  end function sps_dim

  function sps_inuse() result(r) bind(c, name="sps_inuse")
    integer(c_int) :: r
    r = merge(1, 0, acc_inuse)
  end function sps_inuse

  function sps_max_elem() result(n) bind(c, name="sps_max_elem")
    integer(c_int) :: n
    n = max(acc_solver % mesh % maxelementdofs, &
        acc_solver % mesh % maxelementnodes)
  end function sps_max_elem

  function sps_nof_active() result(n) bind(c, name="sps_nof_active")
    integer(c_int) :: n
    n = acc_solver % NumberOfActiveElements
  end function sps_nof_active

  ! -- driver calls ------------------------------------------------------
  subroutine sps_info_header() bind(c, name="sps_info_header")
    call Info('ScalarPotentialSolver', &
        '-------------------------------------', level=4)
    call Info('ScalarPotentialSolver', &
        'Computing scalar potential of a vector field', level=4)
    call Info('ScalarPotentialSolver', &
        '-------------------------------------', level=4)
  end subroutine sps_info_header

  subroutine sps_default_initialize() bind(c, name="sps_default_initialize")
    call DefaultInitialize(acc_solver, acc_inuse)
  end subroutine sps_default_initialize

  subroutine sps_finish_bulk() bind(c, name="sps_finish_bulk")
    if (acc_constbulk) then
      call DefaultFinishBulkAssembly(bulkupdate = .not. acc_inuse, &
          rhsupdate = .false.)
    else
      call DefaultFinishBulkAssembly()
    end if
  end subroutine sps_finish_bulk

  subroutine sps_finish() bind(c, name="sps_finish")
    call DefaultFinishAssembly()
  end subroutine sps_finish

  subroutine sps_dirichlet() bind(c, name="sps_dirichlet")
    call DefaultDirichletBCs()
  end subroutine sps_dirichlet

  function sps_solve() result(norm) bind(c, name="sps_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function sps_solve

  function sps_realtime() result(t) bind(c, name="sps_realtime")
    real(c_double) :: t
    t = RealTime()
  end function sps_realtime

  subroutine sps_info_assembly(t0, t1) bind(c, name="sps_info_assembly")
    real(c_double), value :: t0, t1
    write(message, *) 'Assembly Time: ', t1 - t0
    call Info('ScalarPotentialSolver', message, level=5)
  end subroutine sps_info_assembly

  subroutine sps_info_solution(t1, t2) bind(c, name="sps_info_solution")
    real(c_double), value :: t1, t2
    write(message, *) 'Solution Time: ', t2 - t1
    call Info('ScalarPotentialSolver', message, level=5)
  end subroutine sps_info_solution

  subroutine sps_info_result(norm) bind(c, name="sps_info_result")
    real(c_double), value :: norm
    write(message, *) 'Result Norm: ', norm
    call Info('ScalarPotentialSolver', message, level=4)
  end subroutine sps_info_result

  ! -- BulkAssembly context ----------------------------------------------
  subroutine sps_begin_elem(t, nd, n) bind(c, name="sps_begin_elem")
    integer(c_int), value :: t
    integer(c_int) :: nd, n
    acc_elem => GetActiveElement(t)
    call GetElementNodes(acc_nodes)
    nd = GetElementNOFDOFs()
    n = GetElementNOFNodes()
  end subroutine sps_begin_elem

  subroutine sps_fetch_cond(got, rank) bind(c, name="sps_fetch_cond")
    integer(c_int) :: got, rank
    logical :: found
    call GetRealArray(GetMaterial(), acc_cond, acc_condname, found)
    got = merge(1, 0, found)
    rank = 0
    if (found) then
      ! Exact transcription of the internal GetTensorRank function.
      if (size(acc_cond, 1) == 1) then
        rank = 1
      else if (size(acc_cond, 2) == 1) then
        rank = 2
      else
        rank = 3
      end if
    end if
  end subroutine sps_fetch_cond

  subroutine sps_vector_field(ndim, ld, vf) bind(c, name="sps_vector_field")
    integer(c_int), value :: ndim, ld
    real(c_double) :: vf(3, ld)
    integer :: i
    do i = 1, ndim
      call GetScalarLocalSolution(vf(i, :), ComponentName(acc_varname, i))
    end do
  end subroutine sps_vector_field

  function sps_ip_begin() result(n) bind(c, name="sps_ip_begin")
    integer(c_int) :: n
    acc_ip = GaussPoints(acc_elem)
    n = acc_ip % n
  end function sps_ip_begin

  subroutine sps_ip_point(t, u, v, w, s) bind(c, name="sps_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip % u(t)
    v = acc_ip % v(t)
    w = acc_ip % w(t)
    s = acc_ip % s(t)
  end subroutine sps_ip_point

  subroutine sps_eleminfo(u, v, w, ld, basis, dbasis, detj) &
      bind(c, name="sps_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: ld
    real(c_double) :: basis(ld), dbasis(ld, 3), detj
    logical :: ok
    ok = ElementInfo(acc_elem, acc_nodes, u, v, w, detj, basis, dbasis)
  end subroutine sps_eleminfo

  function sps_weight(s, detj, nelem, ld, basis) result(w) &
      bind(c, name="sps_weight")
    real(c_double), value :: s, detj
    integer(c_int), value :: nelem, ld
    real(c_double) :: basis(ld)
    real(c_double) :: w
    w = s * detj
    if (acc_csym) w = w * sum(basis(1:nelem) * acc_nodes % x(1:nelem))
  end function sps_weight

  subroutine sps_stiff_laplace(nd, ndim, ld, dbasis, weight, stiff) &
      bind(c, name="sps_stiff_laplace")
    integer(c_int), value :: nd, ndim, ld
    real(c_double) :: dbasis(ld, 3), stiff(ld, ld)
    real(c_double), value :: weight
    integer :: p, q
    ! Exact transcription of the no-coefficient update (same compiler,
    ! same evaluation order as the original internal procedure).
    do p = 1, nd
      do q = 1, nd
        stiff(p, q) = stiff(p, q) + weight * &
            sum(dbasis(q, 1:ndim) * dbasis(p, 1:ndim))
      end do
    end do
  end subroutine sps_stiff_laplace

  subroutine sps_stiff_cond(rank, nd, nelem, ndim, ld, basis, dbasis, &
      weight, stiff) bind(c, name="sps_stiff_cond")
    integer(c_int), value :: rank, nd, nelem, ndim, ld
    real(c_double) :: basis(ld), dbasis(ld, 3), stiff(ld, ld)
    real(c_double), value :: weight
    real(kind=dp) :: C(3, 3)
    integer :: p, q, i, j
    ! Exact transcription of the coefficient branch (same compiler,
    ! same evaluation order as the original internal procedure).
    C = 0.0_dp
    select case (rank)
    case (1)
      do i = 1, ndim
        C(i, i) = sum(basis(1:nelem) * acc_cond(1, 1, 1:nelem))
      end do
    case (2)
      do i = 1, ndim
        C(i, i) = sum(basis(1:nelem) * acc_cond(i, 1, 1:nelem))
      end do
    case default
      do i = 1, ndim
        do j = 1, ndim
          C(i, j) = sum(basis(1:nelem) * acc_cond(i, j, 1:nelem))
        end do
      end do
    end select
    do p = 1, nd
      do q = 1, nd
        do i = 1, ndim
          do j = 1, ndim
            stiff(p, q) = stiff(p, q) - weight * &
                dbasis(q, i) * C(i, j) * dbasis(p, j)
          end do
        end do
      end do
    end do
  end subroutine sps_stiff_cond

  subroutine sps_force_update(nd, nelem, ndim, ld, vf, basis, dbasis, &
      weight, force) bind(c, name="sps_force_update")
    integer(c_int), value :: nd, nelem, ndim, ld
    real(c_double) :: vf(3, ld), basis(ld), dbasis(ld, 3), force(ld)
    real(c_double), value :: weight
    integer :: p, i
    ! Exact transcription of the FORCE update (same compiler, same
    ! evaluation order as the original internal procedure).
    do p = 1, nd
      do i = 1, ndim
        force(p) = force(p) + weight * dbasis(p, i) * &
            sum(basis(1:nd) * vf(i, 1:nd))
      end do
    end do
  end subroutine sps_force_update

  ! -- assembly tail ------------------------------------------------------
  subroutine sps_update_eq(ld, nd, stiff, force) bind(c, name="sps_update_eq")
    integer(c_int), value :: ld, nd
    real(c_double) :: stiff(ld, ld), force(ld)
    call DefaultUpdateEquations(stiff(1:ld, 1:ld), force(1:nd))
  end subroutine sps_update_eq

  subroutine sps_update_force(nd, force) bind(c, name="sps_update_force")
    integer(c_int), value :: nd
    real(c_double) :: force(nd)
    call DefaultUpdateForce(force(1:nd))
  end subroutine sps_update_force

  subroutine sps_fatal_alloc() bind(c, name="sps_fatal_alloc")
    call Fatal('ScalarPotentialSolver', 'Memory allocation error.')
  end subroutine sps_fatal_alloc

end module scalarpotentialsolver_acc
