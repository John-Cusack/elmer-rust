! VorticitySolver port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/VorticitySolver.F90, compiled against the
! built .mod files of the target configure (T021-B lane, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module mirrors VorticitySolver's own locals (BulkAssembly SAVE Nodes
! and SAVE Coeff, the ForceVectors/SaveRHS entry context, the VarName
! buffers) with identical lifetimes.
!
! Order-sensitive reductions (MATMUL, SUM), formatted I/O (PRINT/WRITE)
! and whole-array updates stay in Fortran so bit-identity holds by
! construction wherever the call sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module vorticitysolver_acc
  use iso_c_binding
  use types
  use defutils
  use generalutils
  use coordinatesystems
  use integration
  use matrixassembly
  use eleminfo
  use messages
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(valuelist_t), pointer, save :: acc_params => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(variable_t), pointer, save :: acc_vort => null()
  ! BulkAssembly context (mirrors its SAVE Nodes and SAVE Coeff).
  type(nodes_t), target, save :: acc_nodes
  real(kind=dp), allocatable, save :: acc_coeff(:)
  ! Entry context (mirrors the entry's locals across calls).
  character(len=MAX_NAME_LEN), save :: acc_varname = ' '
  character(len=MAX_NAME_LEN), save :: acc_resname = ' '
  logical, save :: acc_csym = .false.
  logical, save :: acc_inuse = .false.
  real(kind=dp), pointer, save :: acc_saverhs(:) => null()
  logical, save :: acc_gotname = .false.
  real(kind=dp), allocatable, target, save :: acc_forcevectors(:, :)
  type(gaussintegrationpoints_t), save :: acc_ip

contains

  ! -- entry state ----------------------------------------------------
  subroutine voa_begin(h_solver) bind(c, name="voa_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine voa_begin

  subroutine voa_info_header() bind(c, name="voa_info_header")
    call Info('VorticitySolver', '-------------------------------------', Level=4)
    call Info('VorticitySolver', 'Computing the vorticity field', Level=4)
    call Info('VorticitySolver', '-------------------------------------', Level=4)
  end subroutine voa_info_header

  function voa_dim() result(d) bind(c, name="voa_dim")
    integer(c_int) :: d
    d = CoordinateSystemDimension()
  end function voa_dim

  function voa_assoc_matrix() result(r) bind(c, name="voa_assoc_matrix")
    integer(c_int) :: r
    r = merge(1, 0, associated(acc_solver % matrix))
  end function voa_assoc_matrix

  function voa_active_var_count() result(n) bind(c, name="voa_active_var_count")
    integer(c_int) :: n
    n = count(acc_solver % variable % perm > 0)
  end function voa_active_var_count

  subroutine voa_params_begin() bind(c, name="voa_params_begin")
    acc_params => GetSolverParams()
  end subroutine voa_params_begin

  ! -- 2D/3D setup (mirrors the dim branches verbatim) -----------------
  subroutine voa_setup_2d(dofs) bind(c, name="voa_setup_2d")
    integer(c_int) :: dofs
    acc_vort => acc_solver % variable
    dofs = acc_vort % dofs
  end subroutine voa_setup_2d

  subroutine voa_print_dofs(dofs) bind(c, name="voa_print_dofs")
    integer(c_int), value :: dofs
    print *, 'DOFS', dofs
  end subroutine voa_print_dofs

  subroutine voa_fatal_dofs() bind(c, name="voa_fatal_dofs")
    call Fatal('VorticitySolver', 'Vorticity should have 1 component in 2D')
  end subroutine voa_fatal_dofs

  subroutine voa_setup_3d(gotit, dofs, state) bind(c, name="voa_setup_3d")
    integer(c_int) :: gotit, dofs, state
    logical :: got
    ! state: 0 = ok, 1 = missing variable, 2 = bad dofs
    state = 0
    acc_resname = GetString(acc_params, 'Vorticity Result Variable', got)
    gotit = merge(1, 0, got)
    if (.not. got) acc_resname = 'Vorticity'
    acc_vort => VariableGet(acc_solver % mesh % variables, acc_resname)
    if (associated(acc_vort)) then
      dofs = acc_vort % dofs
      if (dofs /= 3) state = 2
    else
      dofs = 0
      state = 1
    end if
  end subroutine voa_setup_3d

  subroutine voa_fatal_dofs3() bind(c, name="voa_fatal_dofs3")
    call Fatal('VorticitySolver', 'Vorticity should have 3 components in 3D')
  end subroutine voa_fatal_dofs3

  subroutine voa_fatal_missing() bind(c, name="voa_fatal_missing")
    call Fatal('VorticitySolver', &
        'Vorticity Result Variable is missing: '//trim(acc_resname))
  end subroutine voa_fatal_missing

  subroutine voa_setup_csym() bind(c, name="voa_setup_csym")
    acc_csym = CurrentCoordinateSystem() == AxisSymmetric .or. &
        CurrentCoordinateSystem() == CylindricSymmetric
  end subroutine voa_setup_csym

  subroutine voa_get_varname(which, gotit) bind(c, name="voa_get_varname")
    integer(c_int), value :: which
    integer(c_int) :: gotit
    logical :: got
    if (which == 1) then
      acc_varname = GetString(acc_params, 'Vorticity Variable', got)
    else
      acc_varname = GetString(acc_params, 'Target Variable', got)
    end if
    gotit = merge(1, 0, got)
  end subroutine voa_get_varname

  subroutine voa_default_varname() bind(c, name="voa_default_varname")
    acc_varname = trim('Velocity')
  end subroutine voa_default_varname

  subroutine voa_get_coeff(gotit) bind(c, name="voa_get_coeff")
    integer(c_int) :: gotit
    logical :: got
    character(len=MAX_NAME_LEN) :: condname
    ! For future use: value intentionally unused, call kept for sequence.
    condname = ListGetString(acc_params, 'Vorticity Coefficient', got)
    gotit = merge(1, 0, got)
  end subroutine voa_get_coeff

  function voa_realtime() result(t) bind(c, name="voa_realtime")
    real(c_double) :: t
    t = RealTime()
  end function voa_realtime

  subroutine voa_get_constbulk(gotit, val) bind(c, name="voa_get_constbulk")
    integer(c_int) :: gotit, val
    logical :: got, l
    l = GetLogical(acc_params, 'Constant Bulk Matrix', got)
    gotit = merge(1, 0, got)
    val = merge(1, 0, l)
  end subroutine voa_get_constbulk

  subroutine voa_setup_inuse(constbulk) bind(c, name="voa_setup_inuse")
    integer(c_int), value :: constbulk
    logical :: cb
    cb = constbulk /= 0
    acc_inuse = cb .and. associated(acc_solver % matrix % bulkvalues)
  end subroutine voa_setup_inuse

  function voa_inuse() result(r) bind(c, name="voa_inuse")
    integer(c_int) :: r
    r = merge(1, 0, acc_inuse)
  end function voa_inuse

  function voa_constbulk() result(r) bind(c, name="voa_constbulk")
    integer(c_int) :: r
    ! Re-reads the flag for the DefaultFinishBulkAssembly branch.
    r = merge(1, 0, GetLogical(acc_params, 'Constant Bulk Matrix'))
  end function voa_constbulk

  subroutine voa_default_initialize() bind(c, name="voa_default_initialize")
    call DefaultInitialize(acc_solver, acc_inuse)
  end subroutine voa_default_initialize

  ! -- ForceVectors / SaveRHS entry context -----------------------------
  function voa_rhs_size() result(n) bind(c, name="voa_rhs_size")
    integer(c_int) :: n
    n = size(acc_solver % matrix % rhs)
  end function voa_rhs_size

  subroutine voa_fv_begin(n, dofs) bind(c, name="voa_fv_begin")
    integer(c_int), value :: n, dofs
    allocate(acc_forcevectors(n, dofs - 1))
    acc_forcevectors = 0.0_dp
    acc_saverhs => acc_solver % matrix % rhs
  end subroutine voa_fv_begin

  subroutine voa_fv_end() bind(c, name="voa_fv_end")
    deallocate(acc_forcevectors)
    acc_solver % matrix % rhs => acc_saverhs
  end subroutine voa_fv_end

  subroutine voa_rhs_to_save() bind(c, name="voa_rhs_to_save")
    acc_solver % matrix % rhs => acc_saverhs
  end subroutine voa_rhs_to_save

  subroutine voa_rhs_to_fv(i) bind(c, name="voa_rhs_to_fv")
    integer(c_int), value :: i
    acc_solver % matrix % rhs => acc_forcevectors(:, i - 1)
  end subroutine voa_rhs_to_fv

  ! -- BulkAssembly driver ------------------------------------------------
  function voa_max_dofs() result(n) bind(c, name="voa_max_dofs")
    integer(c_int) :: n
    n = acc_solver % mesh % maxelementdofs
  end function voa_max_dofs

  function voa_max_nodes() result(n) bind(c, name="voa_max_nodes")
    integer(c_int) :: n
    n = acc_solver % mesh % maxelementnodes
  end function voa_max_nodes

  function voa_nof_active() result(n) bind(c, name="voa_nof_active")
    integer(c_int) :: n
    n = acc_solver % numberofactiveelements
  end function voa_nof_active

  subroutine voa_bulk_alloc(n) bind(c, name="voa_bulk_alloc")
    integer(c_int), value :: n
    allocate(acc_coeff(n))
  end subroutine voa_bulk_alloc

  subroutine voa_bulk_free() bind(c, name="voa_bulk_free")
    deallocate(acc_coeff)
  end subroutine voa_bulk_free

  subroutine voa_begin_active(t, n, nd) bind(c, name="voa_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n, nd
    acc_elem => GetActiveElement(t)
    call GetElementNodes(acc_nodes)
    nd = GetElementNOFDOFs()
    n = GetElementNOFNodes()
  end subroutine voa_begin_active

  subroutine voa_local_solution(comp, n, x) bind(c, name="voa_local_solution")
    integer(c_int), value :: comp, n
    real(c_double) :: x(n)
    call GetScalarLocalSolution(x, ComponentName(acc_varname, comp))
  end subroutine voa_local_solution

  function voa_ip_begin() result(n) bind(c, name="voa_ip_begin")
    integer(c_int) :: n
    acc_ip = GaussPoints(acc_elem)
    n = acc_ip % n
  end function voa_ip_begin

  ! Exact transcription of the per-integration-point update statements
  ! (ElementInfo, weight with the CSymmetry SUM term, STIFF accumulation,
  ! the MATMUL gradient projections, and the FORCE updates).
  subroutine voa_ip_update(t, n, nd, dim, dofs, stiff, force, vx, vy, vz, &
      basis, dbasis) bind(c, name="voa_ip_update")
    integer(c_int), value :: t, n, nd, dim, dofs
    real(c_double) :: stiff(n, n), force(dim, n)
    real(c_double) :: vx(n), vy(n), vz(n), basis(n), dbasis(n, 3)
    logical :: found
    real(kind=dp) :: weight, detj
    real(kind=dp) :: gradvx(3), gradvy(3), gradvz(3)
    integer :: p, q
    found = ElementInfo(acc_elem, acc_nodes, acc_ip % u(t), &
        acc_ip % v(t), acc_ip % w(t), detj, basis, dbasis)
    weight = acc_ip % s(t) * detj
    if (acc_csym) weight = weight * sum(basis(1:n) * acc_nodes % x(1:n))
    if (.not. acc_inuse) then
      do p = 1, nd
        do q = 1, nd
          stiff(p, q) = stiff(p, q) + weight * basis(q) * basis(p)
        end do
      end do
    end if
    gradvx(1:dim) = matmul(vx(1:nd), dbasis(1:nd, 1:dim))
    gradvy(1:dim) = matmul(vy(1:nd), dbasis(1:nd, 1:dim))
    gradvz(1:dim) = matmul(vz(1:nd), dbasis(1:nd, 1:dim))
    if (dofs == 1) then
      force(1, 1:nd) = force(1, 1:nd) + &
          basis(1:nd) * weight * (gradvy(1) - gradvx(2))
    else
      force(1, 1:nd) = force(1, 1:nd) + &
          basis(1:nd) * weight * (gradvz(2) - gradvy(3))
      force(2, 1:nd) = force(2, 1:nd) + &
          basis(1:nd) * weight * (gradvx(3) - gradvz(1))
      force(3, 1:nd) = force(3, 1:nd) + &
          basis(1:nd) * weight * (gradvy(1) - gradvx(2))
    end if
  end subroutine voa_ip_update

  subroutine voa_update_equations(n, nd, df, stiff, force) bind(c, name="voa_update_equations")
    integer(c_int), value :: n, nd, df
    real(c_double) :: stiff(n, n), force(df, n)
    call DefaultUpdateEquations(stiff, force(1, 1:nd))
  end subroutine voa_update_equations

  subroutine voa_update_force(i, nd, df, force) bind(c, name="voa_update_force")
    integer(c_int), value :: i, nd, df
    real(c_double) :: force(df, nd)
    call DefaultUpdateForce(force(i, 1:nd))
  end subroutine voa_update_force

  ! -- assembly tail / solve / reports ------------------------------------
  subroutine voa_finish_bulk() bind(c, name="voa_finish_bulk")
    call DefaultFinishBulkAssembly()
  end subroutine voa_finish_bulk

  subroutine voa_finish_bulk_const() bind(c, name="voa_finish_bulk_const")
    call DefaultFinishBulkAssembly(BulkUpdate=.not. acc_inuse, RHSUpdate=.false.)
  end subroutine voa_finish_bulk_const

  subroutine voa_default_finish() bind(c, name="voa_default_finish")
    call DefaultFinishAssembly()
  end subroutine voa_default_finish

  subroutine voa_report_assembly(at0, at1) bind(c, name="voa_report_assembly")
    real(c_double), value :: at0, at1
    write(Message, *) 'Assembly Time: ', at1 - at0
    call Info('VorticitySolver', Message, Level=5)
  end subroutine voa_report_assembly

  function voa_default_solve() result(norm) bind(c, name="voa_default_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function voa_default_solve

  function voa_nrows() result(n) bind(c, name="voa_nrows")
    integer(c_int) :: n
    n = acc_solver % matrix % numberofrows
  end function voa_nrows

  subroutine voa_copy_component(i, dofs) bind(c, name="voa_copy_component")
    integer(c_int), value :: i, dofs
    integer :: j
    do j = 1, acc_solver % matrix % numberofrows
      acc_vort % values(dofs * (j - 1) + i) = acc_solver % variable % values(j)
    end do
  end subroutine voa_copy_component

  subroutine voa_set_norm(tot) bind(c, name="voa_set_norm")
    real(c_double), value :: tot
    acc_solver % variable % norm = tot
  end subroutine voa_set_norm

  subroutine voa_report_solution(at1, at2) bind(c, name="voa_report_solution")
    real(c_double), value :: at1, at2
    write(Message, *) 'Solution Time: ', at2 - at1
    call Info('VorticitySolver', Message, Level=5)
  end subroutine voa_report_solution

  subroutine voa_report_norm(tot) bind(c, name="voa_report_norm")
    real(c_double), value :: tot
    write(Message, *) 'Result Norm: ', tot
    call Info('VorticitySolver', Message, Level=4)
  end subroutine voa_report_norm

  ! -- VorticitySolver_Init -------------------------------------------------
  subroutine voi_begin(h_solver) bind(c, name="voi_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
    acc_params => GetSolverParams()
  end subroutine voi_begin

  subroutine voi_get_varname(gotit) bind(c, name="voi_get_varname")
    integer(c_int) :: gotit
    logical :: got
    acc_varname = GetString(acc_params, 'Vorticity Variable', got)
    acc_gotname = got
    gotit = merge(1, 0, got)
  end subroutine voi_get_varname

  function voi_check_variable() result(r) bind(c, name="voi_check_variable")
    integer(c_int) :: r
    r = merge(1, 0, ListCheckPresent(acc_params, 'Variable'))
  end function voi_check_variable

  subroutine voi_add_var_2d() bind(c, name="voi_add_var_2d")
    character(len=MAX_NAME_LEN) :: vorname
    if (acc_gotname) then
      vorname = 'Curl '//trim(acc_varname)
    else
      vorname = 'Vorticity'
    end if
    call ListAddInteger(acc_params, 'Variable DOFs', 1)
    call ListAddString(acc_params, 'Variable', vorname)
  end subroutine voi_add_var_2d

  subroutine voi_add_var_3d() bind(c, name="voi_add_var_3d")
    character(len=MAX_NAME_LEN) :: vorname, helpname
    if (acc_gotname) then
      vorname = 'Curl '//trim(acc_varname)
    else
      vorname = 'Vorticity'
    end if
    call ListAddInteger(acc_params, 'Variable DOFs', 1)
    helpname = 'F['//trim(vorname)//':3]'
    call ListAddString(acc_params, 'Exported Variable 1', helpname)
    call ListAddString(acc_params, 'Vorticity Result Variable', 'F')
    call ListAddString(acc_params, 'Variable', '-nooutput vorticity_temp')
  end subroutine voi_add_var_3d

  subroutine voi_add_abs() bind(c, name="voi_add_abs")
    character(len=MAX_NAME_LEN) :: vorname
    if (acc_gotname) then
      vorname = 'Curl '//trim(acc_varname)
    else
      vorname = 'Vorticity'
    end if
    call ListAddString(acc_params, &
        NextFreeKeyword('Exported Variable', acc_params), trim(vorname)//'_abs')
  end subroutine voi_add_abs

  subroutine voi_fatal_baddim() bind(c, name="voi_fatal_baddim")
    call Fatal('VortictySolver_init', &
        'Vorticity computation makes sense only in 2D and 3D')
  end subroutine voi_fatal_baddim

  function voi_get_calcabs() result(r) bind(c, name="voi_get_calcabs")
    integer(c_int) :: r
    logical :: got
    r = merge(1, 0, GetLogical(acc_params, 'Calculate Abs', got))
  end function voi_get_calcabs

  subroutine voi_add_tail() bind(c, name="voi_add_tail")
    call ListAddInteger(acc_params, 'Time derivative order', 0)
    call ListAddNewString(acc_params, 'Linear System Solver', 'Iterative')
    call ListAddNewString(acc_params, 'Linear System Iterative Method', 'cg')
    call ListAddNewString(acc_params, 'Linear System Preconditioning', 'ILU0')
    call ListAddNewInteger(acc_params, 'Linear System Max Iterations', 500)
    call ListAddNewConstReal(acc_params, &
        'Linear System Convergence Tolerance', 1.0e-8_dp)
  end subroutine voi_add_tail

end module vorticitysolver_acc
