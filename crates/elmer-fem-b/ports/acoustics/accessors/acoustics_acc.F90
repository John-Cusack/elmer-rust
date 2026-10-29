! Port of `fem/src/modules/Acoustics.F90` (Elmer @ c0b4bc6, the campaign
! fork): statement-level BIND(C) wrappers for the Elmer calls in the
! `AcousticsSolver` driver, compiled against the built .mod files of the
! target configure (A2 shipped mechanism, ADR-057 S3).
!
! Every wrapper takes only C-interoperable arguments (TYPE(C_PTR) handles,
! VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All Fortran-object
! traffic (TYPE(...) handles, OPTIONAL/POINTER dummies, assumed-shape
! descriptors, hidden lengths) stays inside Fortran, so the Rust side is
! pure C calls. The wrapper module mirrors the original's host-associated
! element context (current element/parent, material and BC lists, Nodes
! SAVE state per local-matrix context, Gauss-point state) with identical
! lifetimes.
!
! Complex arithmetic (complex multiply/divide, whose libgcc lowering the
! Rust side must not second-guess), order-sensitive reductions
! (SUM/MATMUL), formatted/list-directed I/O (WRITE/PRINT) and whole-array
! updates stay in Fortran wrappers, so bit-identity holds by construction
! wherever the call sequence matches. The BlockPreconditioning-only
! iterative suite and the BEM/Compa/slip/fluidic-force-only kernels are
! preserved verbatim as hidden module procedures (same statements, same
! compiler) and entered only through thin BIND(C) shims under the exact
! original branches; they never execute in the proof configs.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text beyond the ported solver's own statements.
module acoustics_acc
  use iso_c_binding
  use types
  use defutils
  use lists
  use integration
  use matrixassembly
  use eleminfo
  use elementutils
  use solverbasics
  use boundaryconditionutils
  use crsmatrix
  use coordinatesystems
  use elementgeometry
  use generalutils
  use linearalgebra
  use messages
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(element_t), pointer, save :: acc_parent => null()
  type(valuelist_t), pointer, save :: acc_mat => null()
  type(valuelist_t), pointer, save :: acc_bc => null()
  type(valuelist_t), pointer, save :: acc_bf => null()
  integer, pointer, save :: acc_nidx(:) => null()
  integer, save :: acc_nn = 0
  ! Local-matrix context (mirrors the internal-procedure SAVE Nodes/IP).
  type(nodes_t), target, save :: acc_nodes
  type(nodes_t), target, save :: acc_pnodes
  real(kind=dp), pointer, save :: acc_h(:) => null()
  integer, pointer, save :: acc_hp(:) => null()
  type(gaussintegrationpoints_t), save :: acc_ip
  character(len=MAX_NAME_LEN), save :: acc_message
  character(len=MAX_NAME_LEN), save :: acc_eqname = ''
  ! BEM-coupling SAVE mirror (original SAVE lifetimes, Fortran-owned).
  integer, allocatable, save :: acc_bem_elem(:)
  integer, allocatable, save :: acc_ainodes(:)
  real(kind=dp), allocatable, save :: acc_aires(:,:)
  logical, save :: acc_bem_created = .false.
  ! BlockPreconditioning SAVE mirror (original SAVE lifetimes, Fortran-owned).
  type(matrix_t), pointer, save :: acc_amatrix => null()
  type(matrix_t), pointer, save :: acc_smatrix => null()
  type(matrix_t), pointer, save :: acc_mmatrix => null()
  type(matrix_t), pointer, save :: acc_lmatrix => null()
  real(kind=dp), allocatable, save :: acc_slocal(:,:), acc_alocal(:,:), acc_anlocal(:,:)

contains

  ! -- init entry ------------------------------------------------------
  subroutine aca_init_begin(h_solver) bind(c, name="aca_init_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine aca_init_begin

  subroutine aca_info_init() bind(c, name="aca_info_init")
    call Info('AcousticsSolver','Initialization the solver')
  end subroutine aca_info_init

  subroutine aca_init_add_complex() bind(c, name="aca_init_add_complex")
    type(valuelist_t), pointer :: params
    params => GetSolverParams()
    call ListAddNewLogical(params,'Linear System Complex',.TRUE.)
  end subroutine aca_init_add_complex

  subroutine aca_init_add_variable(is2d) bind(c, name="aca_init_add_variable")
    integer(c_int), value :: is2d
    type(valuelist_t), pointer :: params
    params => GetSolverParams()
    if (is2d /= 0) then
      call ListAddNewString(params,'Variable',&
          'Flow[Re Velocity 1:1 Im Velocity 1:1 Re Velocity 2:1 Im Velocity 2:1 '&
          //' Re Temperature:1 Im Temperature:1 Re Pressure:1 Im Pressure:1]')
    else
      call ListAddNewString(params,'Variable',&
          'Flow[Re Velocity 1:1 Im Velocity 1:1 Re Velocity 2:1 Im Velocity 2:1 Re Velocity 3:1 Im Velocity 3:1 '&
          //' Re Temperature:1 Im Temperature:1 Re Pressure:1 Im Pressure:1]')
    end if
  end subroutine aca_init_add_variable

  subroutine aca_begin(h_model, h_solver) bind(c, name="aca_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine aca_begin

  subroutine aca_alloc_nodes(nn) bind(c, name="aca_alloc_nodes")
    integer(c_int), value :: nn
    integer :: istat
    if (associated(acc_nodes % x)) deallocate(acc_nodes % x, acc_nodes % y, acc_nodes % z, &
        acc_pnodes % x, acc_pnodes % y, acc_pnodes % z, stat=istat)
    allocate(acc_nodes % x(nn), acc_nodes % y(nn), acc_nodes % z(nn), &
        acc_pnodes % x(nn), acc_pnodes % y(nn), acc_pnodes % z(nn), stat=istat)
    if (istat /= 0) call Fatal('AcousticsSolver', 'Memory allocation error.')
    acc_nodes % x = 0.0d0
    acc_nodes % y = 0.0d0
    acc_nodes % z = 0.0d0
  end subroutine aca_alloc_nodes

  function aca_coord_dim() result(n) bind(c, name="aca_coord_dim")
    integer(c_int) :: n
    n = CoordinateSystemDimension()
  end function aca_coord_dim

  ! -- driver guards and model queries -----------------------------------

  function aca_assoc_matrix() result(r) bind(c, name="aca_assoc_matrix")
    integer(c_int) :: r
    r = merge(1, 0, associated(acc_solver % matrix))
  end function aca_assoc_matrix

  function aca_flow_vals() result(p) bind(c, name="aca_flow_vals")
    type(c_ptr) :: p
    p = c_loc(acc_solver % variable % values(1))
  end function aca_flow_vals

  function aca_flow_perm() result(p) bind(c, name="aca_flow_perm")
    type(c_ptr) :: p
    p = c_loc(acc_solver % variable % perm(1))
  end function aca_flow_perm
  function aca_model_num_bcs() result(n) bind(c, name="aca_model_num_bcs")
    integer(c_int) :: n
    n = acc_model % numberofbcs
  end function aca_model_num_bcs

  function aca_model_num_nodes() result(n) bind(c, name="aca_model_num_nodes")
    integer(c_int) :: n
    n = acc_model % numberofnodes
  end function aca_model_num_nodes
  function aca_flow_dofs() result(n) bind(c, name="aca_flow_dofs")
    integer(c_int) :: n
    n = acc_solver % variable % dofs
  end function aca_flow_dofs

  function aca_var_norm() result(x) bind(c, name="aca_var_norm")
    real(c_double) :: x
    x = acc_solver % variable % norm
  end function aca_var_norm

  subroutine aca_var_norm_set(x) bind(c, name="aca_var_norm_set")
    real(c_double), value :: x
    acc_solver % variable % norm = x
  end subroutine aca_var_norm_set

  function aca_local_nodes() result(n) bind(c, name="aca_local_nodes")
    integer(c_int) :: n
    n = count(acc_solver % variable % perm > 0)
  end function aca_local_nodes

  function aca_num_nodes() result(n) bind(c, name="aca_num_nodes")
    integer(c_int) :: n
    n = acc_solver % mesh % numberofnodes
  end function aca_num_nodes


  function aca_mesh_changed() result(r) bind(c, name="aca_mesh_changed")
    integer(c_int) :: r
    r = merge(1, 0, acc_solver % meshchanged)
  end function aca_mesh_changed

  function aca_max_elem_dofs() result(n) bind(c, name="aca_max_elem_dofs")
    integer(c_int) :: n
    n = acc_solver % mesh % maxelementdofs
  end function aca_max_elem_dofs

  function aca_num_active() result(n) bind(c, name="aca_num_active")
    integer(c_int) :: n
    n = acc_solver % numberofactiveelements
  end function aca_num_active

  function aca_num_bulk() result(n) bind(c, name="aca_num_bulk")
    integer(c_int) :: n
    n = acc_solver % mesh % numberofbulkelements
  end function aca_num_bulk

  function aca_num_bnd_elems() result(n) bind(c, name="aca_num_bnd_elems")
    integer(c_int) :: n
    n = acc_solver % mesh % numberofboundaryelements
  end function aca_num_bnd_elems
  ! -- simulation/solution option reads ------------------------------------
  subroutine aca_sim_type(s, n, got) bind(c, name="aca_sim_type")
    integer(c_int), value :: n
    character(kind=c_char) :: s(n)
    integer(c_int) :: got
    character(len=MAX_NAME_LEN) :: str
    logical :: f
    integer :: i, m
    str = ListGetString(acc_model % simulation, 'Simulation Type', f)
    got = merge(1, 0, f)
    m = len_trim(str)
    do i = 1, n
      if (i <= m) then
        s(i) = str(i:i)
      else if (i == m + 1) then
        s(i) = c_null_char
      else
        s(i) = ' '
      end if
    end do
  end subroutine aca_sim_type

  subroutine aca_time_var(found, v1) bind(c, name="aca_time_var")
    integer(c_int) :: found
    real(c_double) :: v1
    type(variable_t), pointer :: tvar
    tvar => VariableGet(acc_model % mesh % variables, 'Time')
    if (.not. associated(tvar)) then
      found = 0
      v1 = 0.0_c_double
    else
      found = 1
      v1 = tvar % values(1)
    end if
  end subroutine aca_time_var

  subroutine aca_fatal_timevar() bind(c, name="aca_fatal_timevar")
    call Fatal('AcousticsSolver', 'Error in reading time variable.')
  end subroutine aca_fatal_timevar

  subroutine aca_sim_logical(which, val, got) bind(c, name="aca_sim_logical")
    integer(c_int), value :: which
    integer(c_int) :: val, got
    logical :: l, f
    if (which == 1) then
      l = ListGetLogical(acc_model % simulation, 'BEM Coupling', f)
    else
      l = ListGetLogical(acc_model % simulation, 'DD Preconditioning', f)
    end if
    val = merge(1, 0, l)
    got = merge(1, 0, f)
  end subroutine aca_sim_logical

  subroutine aca_find_pres(found) bind(c, name="aca_find_pres")
    integer(c_int) :: found
    type(variable_t), pointer :: hsol, tmp
    hsol => null()
    tmp => acc_solver % mesh % variables
    do while (associated(tmp))
      if (tmp % namelen == 4) then
        if (tmp % name == 'pres') then
          if (tmp % valid) then
            hsol => tmp
          end if
          exit
        end if
      end if
      tmp => tmp % next
    end do
    found = merge(1, 0, associated(hsol))
    if (associated(hsol)) then
      acc_h => hsol % values
      acc_hp => hsol % perm
    else
      nullify(acc_h)
      nullify(acc_hp)
    end if
  end subroutine aca_find_pres
  subroutine aca_fatal_helmholtz_node() bind(c, name="aca_fatal_helmholtz_node")
    call Fatal('AcousticsSolver', 'Helmholtz solution is not available on the node')
  end subroutine aca_fatal_helmholtz_node

  subroutine aca_fatal_helmholtz_bnd() bind(c, name="aca_fatal_helmholtz_bnd")
    call Fatal('AcousticsSolver', 'Helmholtz solution is not available on boundary')
  end subroutine aca_fatal_helmholtz_bnd
  subroutine aca_var_len(nv, np) bind(c, name="aca_var_len")
    integer(c_int) :: nv, np
    nv = size(acc_solver % variable % values)
    np = size(acc_solver % variable % perm)
  end subroutine aca_var_len

  subroutine aca_info_no_pres() bind(c, name="aca_info_no_pres")
    call Info('AcousticsSolver', 'Helmholtz solution is not available')
    call Info('AcousticsSolver', 'Acoustic interface boundary conditions cannot be applied')
  end subroutine aca_info_no_pres

  subroutine aca_fatal_no_pres() bind(c, name="aca_fatal_no_pres")
    call Fatal('AcousticsSolver', 'Helmholtz solution Pres is not available')
  end subroutine aca_fatal_no_pres

  subroutine aca_solver_logical(which, val, got) bind(c, name="aca_solver_logical")
    integer(c_int), value :: which
    integer(c_int) :: val, got
    logical :: l, f
    select case (which)
    case (1)
      l = ListGetLogical(acc_solver % values, 'Block Preconditioning', f)
    case (2)
      l = ListGetLogical(acc_solver % values, 'Bubbles', f)
    case (3)
      l = ListGetLogical(acc_solver % values, 'Mini Bubbles', f)
    case (4)
      l = ListGetLogical(acc_solver % values, 'Utilize Previous Solution', f)
    case (5)
      l = ListGetLogical(acc_solver % values, 'Velocity Assembly', f)
    case default
      l = ListGetLogical(acc_solver % values, 'Optimize Bandwidth', f)
    end select
    val = merge(1, 0, l)
    got = merge(1, 0, f)
  end subroutine aca_solver_logical

  subroutine aca_solver_string(which, s, n, got) bind(c, name="aca_solver_string")
    integer(c_int), value :: which, n
    character(kind=c_char) :: s(n)
    integer(c_int) :: got
    character(len=MAX_NAME_LEN) :: str
    logical :: f
    integer :: i, m
    if (which == 1) then
      str = ListGetString(acc_solver % values, 'Outer Iteration Method', f)
    else
      str = ListGetString(acc_solver % values, 'Equation', f)
    end if
    got = merge(1, 0, f)
    m = len_trim(str)
    do i = 1, n
      if (i <= m) then
        s(i) = str(i:i)
      else if (i == m + 1) then
        s(i) = c_null_char
      else
        s(i) = ' '
      end if
    end do
  end subroutine aca_solver_string

  subroutine aca_set_eqname() bind(c, name="aca_set_eqname")
    logical :: f
    acc_eqname = ListGetString(acc_solver % values, 'Equation', f)
  end subroutine aca_set_eqname

  subroutine aca_solver_creal(which, v, got) bind(c, name="aca_solver_creal")
    integer(c_int), value :: which
    real(c_double) :: v
    integer(c_int) :: got
    logical :: f
    if (which == 1) then
      v = ListGetConstReal(acc_solver % values, 'Ratio of Convergence Tolerances', f)
    else if (which == 2) then
      v = ListGetConstReal(acc_solver % values, 'Linear System Convergence Tolerance', f)
    else
      v = ListGetConstReal(acc_solver % values, 'Steady State Convergence Tolerance', f)
    end if
    got = merge(1, 0, f)
  end subroutine aca_solver_creal

  subroutine aca_solver_int(which, v, got) bind(c, name="aca_solver_int")
    integer(c_int), value :: which
    integer(c_int) :: v, got
    logical :: f
    if (which == 1) then
      v = ListGetInteger(acc_solver % values, 'Relative Integration Order', f)
    else if (which == 2) then
      v = ListGetInteger(acc_solver % values, 'Max Outer Iterations', f)
    else if (which == 3) then
      v = ListGetInteger(acc_solver % values, 'Linear System Max Iterations', f)
    else
      v = ListGetInteger(acc_solver % values, 'BiCGStab Polynomial Degree', f)
    end if
    got = merge(1, 0, f)
  end subroutine aca_solver_int

  subroutine aca_angular_freq(af, got_af, got_f) bind(c, name="aca_angular_freq")
    real(c_double) :: af
    integer(c_int) :: got_af, got_f
    logical :: f1, f2
    real(kind=dp), parameter :: pi = 3.14159265358979323846264338327950288_dp
    af = GetCReal(acc_model % simulation, 'Angular Frequency', f1)
    got_af = merge(1, 0, f1)
    if (.not. f1) then
      af = 2 * pi * GetCReal(acc_model % simulation, 'Frequency', f2)
      got_f = merge(1, 0, f2)
    else
      got_f = 0
    end if
  end subroutine aca_angular_freq

  subroutine aca_check_coordsys() bind(c, name="aca_check_coordsys")
    integer :: cs
    cs = CurrentCoordinateSystem()
    if (.not. (cs == Cartesian .or. cs == AxisSymmetric)) then
      call Fatal('AcousticsSolver', &
          'Currently only Cartesian or cylindrical coordinates are allowed')
    end if
  end subroutine aca_check_coordsys

  function aca_cur_coordsys() result(n) bind(c, name="aca_cur_coordsys")
    integer(c_int) :: n
    n = CurrentCoordinateSystem()
  end function aca_cur_coordsys

  subroutine aca_relorder_diag(relorder) bind(c, name="aca_relorder_diag")
    integer(c_int), value :: relorder
    type(element_t), pointer :: el
    type(gaussintegrationpoints_t) :: ip
    integer :: i, j
    if (ParEnv % MyPe == 0) then
      el => acc_solver % mesh % elements(acc_solver % activeelements(1))
      ip = GaussPoints(el)
      i = ip % n
      ip = GaussPoints(el, RelOrder=relorder)
      j = ip % n
      call Info('AcousticsSolver', 'Number of Gauss Points: ' &
          //I2S(j)//' (vs. '//I2S(i)//')', Level=6)
    end if
  end subroutine aca_relorder_diag
  ! -- bulk/boundary element context ---------------------------------------
  subroutine aca_begin_bulk(t, n, nd, nb, bodyid) bind(c, name="aca_begin_bulk")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, nb, bodyid
    acc_elem => acc_solver % mesh % elements(acc_solver % activeelements(t))
    acc_model % currentelement => acc_elem
    CurrentModel % currentelement => acc_elem
    acc_nidx => acc_elem % nodeindexes
    acc_nn = acc_elem % type % numberofnodes
    n = acc_nn
    nd = GetElementNOFDOFs()
    nb = GetElementNOFBDOFs()
    bodyid = acc_elem % bodyid
    acc_nidx => acc_elem % nodeindexes
  end subroutine aca_begin_bulk

  subroutine aca_elem_nodes() bind(c, name="aca_elem_nodes")
    integer :: n
    n = acc_nn
    acc_nodes % x(1:n) = acc_solver % mesh % nodes % x(acc_nidx(1:n))
    acc_nodes % y(1:n) = acc_solver % mesh % nodes % y(acc_nidx(1:n))
    acc_nodes % z(1:n) = acc_solver % mesh % nodes % z(acc_nidx(1:n))
  end subroutine aca_elem_nodes

  function aca_nidx_val(j) result(v) bind(c, name="aca_nidx_val")
    integer(c_int), value :: j
    integer(c_int) :: v
    v = acc_nidx(j)
  end function aca_nidx_val

  function aca_perm_of(node) result(v) bind(c, name="aca_perm_of")
    integer(c_int), value :: node
    integer(c_int) :: v
    v = acc_solver % variable % perm(node)
  end function aca_perm_of

  function aca_perm_has_zero(n) result(r) bind(c, name="aca_perm_has_zero")
    integer(c_int), value :: n
    integer(c_int) :: r
    if (any(acc_solver % variable % perm(acc_nidx(1:n)) == 0)) then
      r = 1
    else
      r = 0
    end if
  end function aca_perm_has_zero

  subroutine aca_begin_bnd(t, n, active, constraint) bind(c, name="aca_begin_bnd")
    integer(c_int), value :: t
    integer(c_int) :: n, active, constraint
    acc_elem => acc_solver % mesh % elements(t)
    acc_model % currentelement => acc_elem
    CurrentModel % currentelement => acc_elem
    acc_nidx => acc_elem % nodeindexes
    acc_nn = acc_elem % type % numberofnodes
    n = acc_nn
    if (ActiveBoundaryElement(acc_elem, CurrentModel % solver)) then
      active = 1
    else
      active = 0
    end if
    constraint = acc_elem % boundaryinfo % constraint
  end subroutine aca_begin_bnd

  function aca_bc_tag(i) result(v) bind(c, name="aca_bc_tag")
    integer(c_int), value :: i
    integer(c_int) :: v
    v = acc_model % bcs(i) % tag
  end function aca_bc_tag

  subroutine aca_use_bc(i) bind(c, name="aca_use_bc")
    integer(c_int), value :: i
    acc_bc => acc_model % bcs(i) % values
  end subroutine aca_use_bc

  function aca_check_eqn() result(r) bind(c, name="aca_check_eqn")
    integer(c_int) :: r
    if (CheckElementEquation(acc_model, acc_elem, acc_eqname)) then
      r = 1
    else
      r = 0
    end if
  end function aca_check_eqn

  function aca_check_parent_eqn() result(r) bind(c, name="aca_check_parent_eqn")
    integer(c_int) :: r
    if (CheckElementEquation(acc_model, acc_parent, acc_eqname)) then
      r = 1
    else
      r = 0
    end if
  end function aca_check_parent_eqn

  subroutine aca_parent_left(ok) bind(c, name="aca_parent_left")
    integer(c_int) :: ok
    logical :: stat
    acc_parent => acc_elem % boundaryinfo % left
    stat = associated(acc_parent)
    if (stat) stat = all(acc_solver % variable % perm( &
        acc_parent % nodeindexes(1:acc_parent % type % numberofnodes)) > 0)
    ok = merge(1, 0, stat)
  end subroutine aca_parent_left

  subroutine aca_parent_right(ok) bind(c, name="aca_parent_right")
    integer(c_int) :: ok
    logical :: stat
    acc_parent => acc_elem % boundaryinfo % right
    stat = associated(acc_parent)
    if (stat) stat = all(acc_solver % variable % perm( &
        acc_parent % nodeindexes(1:acc_parent % type % numberofnodes)) > 0)
    ok = merge(1, 0, stat)
  end subroutine aca_parent_right

  subroutine aca_fatal_no_parent() bind(c, name="aca_fatal_no_parent")
    call Fatal('AcousticsSolver', 'No parent element can be found for given boundary element')
  end subroutine aca_fatal_no_parent

  ! -- material / body-force / BC property reads -----------------------------
  function aca_body_material() result(k) bind(c, name="aca_body_material")
    integer(c_int) :: k
    logical :: f
    k = ListGetInteger(acc_model % bodies(acc_elem % bodyid) % values, 'Material')
  end function aca_body_material

  function aca_parent_material() result(k) bind(c, name="aca_parent_material")
    integer(c_int) :: k
    k = ListGetInteger(acc_model % bodies(acc_parent % bodyid) % values, 'Material')
  end function aca_parent_material

  subroutine aca_use_material(k) bind(c, name="aca_use_material")
    integer(c_int), value :: k
    acc_mat => acc_model % materials(k) % values
  end subroutine aca_use_material

  subroutine aca_mat_real(which, n, out, got) bind(c, name="aca_mat_real")
    integer(c_int), value :: which, n
    real(c_double) :: out(n)
    integer(c_int) :: got
    logical :: f
    f = .false.
    select case (which)
    case (1)
      out(1:n) = ListGetReal(acc_mat, 'Specific Heat', n, acc_nidx)
    case (2)
      out(1:n) = ListGetReal(acc_mat, 'Specific Heat Ratio', n, acc_nidx)
    case (3)
      out(1:n) = ListGetReal(acc_mat, 'Equilibrium Density', n, acc_nidx)
    case (4)
      out(1:n) = ListGetReal(acc_mat, 'Equilibrium Temperature', n, acc_nidx)
    case (5)
      out(1:n) = ListGetReal(acc_mat, 'Heat Conductivity', n, acc_nidx)
    case (6)
      out(1:n) = ListGetReal(acc_mat, 'Viscosity', n, acc_nidx)
    case (7)
      out(1:n) = ListGetReal(acc_mat, ' Bulk Viscosity', n, acc_nidx, f)
    case default
      out(1:n) = ListGetReal(acc_mat, 'Bulk Viscosity', n, acc_nidx, f)
    end select
    got = merge(1, 0, f)
  end subroutine aca_mat_real

  subroutine aca_body_bfcount(k, got) bind(c, name="aca_body_bfcount")
    integer(c_int) :: k, got
    logical :: f
    k = ListGetInteger(acc_model % bodies(acc_elem % bodyid) % values, 'Body Force', f)
    got = merge(1, 0, f)
  end subroutine aca_body_bfcount

  subroutine aca_use_bodyforce(k) bind(c, name="aca_use_bodyforce")
    integer(c_int), value :: k
    acc_bf => acc_model % bodyforces(k) % values
  end subroutine aca_use_bodyforce

  subroutine aca_bf_real(which, n, out, got) bind(c, name="aca_bf_real")
    integer(c_int), value :: which, n
    real(c_double) :: out(n)
    integer(c_int) :: got
    logical :: f
    select case (which)
    case (1)
      out(1:n) = ListGetReal(acc_bf, 'Re Heat Source', n, acc_nidx, f)
    case (2)
      out(1:n) = ListGetReal(acc_bf, 'Im Heat Source', n, acc_nidx, f)
    case (3)
      out(1:n) = ListGetReal(acc_bf, 'Re Body Force 1', n, acc_nidx, f)
    case (4)
      out(1:n) = ListGetReal(acc_bf, 'Im Body Force 1', n, acc_nidx, f)
    case (5)
      out(1:n) = ListGetReal(acc_bf, 'Re Body Force 2', n, acc_nidx, f)
    case (6)
      out(1:n) = ListGetReal(acc_bf, 'Im Body Force 2', n, acc_nidx, f)
    case (7)
      out(1:n) = ListGetReal(acc_bf, 'Re Body Force 3', n, acc_nidx, f)
    case default
      out(1:n) = ListGetReal(acc_bf, 'Im Body Force 3', n, acc_nidx, f)
    end select
    got = merge(1, 0, f)
  end subroutine aca_bf_real

  subroutine aca_bc_logical(which, val, got) bind(c, name="aca_bc_logical")
    integer(c_int), value :: which
    integer(c_int) :: val, got
    logical :: l, f
    select case (which)
    case (1)
      l = ListGetLogical(acc_bc, 'Acoustic Interface', f)
    case (2)
      l = ListGetLogical(acc_bc, 'Impedance Update', f)
    case (3)
      l = ListGetLogical(acc_bc, 'Potential Flow BC', f)
    case (4)
      l = ListGetLogical(acc_bc, 'Slip Boundary', f)
    case (5)
      l = ListGetLogical(acc_bc, 'Traction Boundary', f)
    case (6)
      l = ListGetLogical(acc_bc, 'BEM Boundary', f)
    case (7)
      l = ListGetLogical(acc_bc, 'Calculate Fluidic Force', f)
    case (8)
      l = ListGetLogical(acc_bc, 'Calculate Acoustic Impedance', f)
    case default
      l = ListGetLogical(acc_bc, 'Impedance Target Boundary', f)
    end select
    val = merge(1, 0, l)
    got = merge(1, 0, f)
  end subroutine aca_bc_logical

  subroutine aca_bc_real(which, n, out, got) bind(c, name="aca_bc_real")
    integer(c_int), value :: which, n
    real(c_double) :: out(n)
    integer(c_int) :: got
    logical :: f
    select case (which)
    case (1)
      out(1:n) = ListGetReal(acc_bc, 'Re Specific Acoustic Impedance', n, acc_nidx, f)
    case (2)
      out(1:n) = ListGetReal(acc_bc, 'Im Specific Acoustic Impedance', n, acc_nidx, f)
    case (3)
      out(1:n) = ListGetReal(acc_bc, 'Re Specific Thermal Impedance', n, acc_nidx, f)
    case (4)
      out(1:n) = ListGetReal(acc_bc, 'Im Specific Thermal Impedance', n, acc_nidx, f)
    case (5)
      out(1:n) = ListGetReal(acc_bc, 'Re Surface Traction 1', n, acc_nidx, f)
    case (6)
      out(1:n) = ListGetReal(acc_bc, 'Im Surface Traction 1', n, acc_nidx, f)
    case (7)
      out(1:n) = ListGetReal(acc_bc, 'Re Surface Traction 2', n, acc_nidx, f)
    case (8)
      out(1:n) = ListGetReal(acc_bc, 'Im Surface Traction 2', n, acc_nidx, f)
    case (9)
      out(1:n) = ListGetReal(acc_bc, 'Re Surface Traction 3', n, acc_nidx, f)
    case (10)
      out(1:n) = ListGetReal(acc_bc, 'Im Surface Traction 3', n, acc_nidx, f)
    case (11)
      out(1:n) = ListGetReal(acc_bc, 'Reference Wall Temperature', n, acc_nidx, f)
    case (12)
      out(1:n) = ListGetReal(acc_bc, 'Re Reference Wall Velocity 1', n, acc_nidx, f)
    case (13)
      out(1:n) = ListGetReal(acc_bc, 'Im Reference Wall Velocity 1', n, acc_nidx, f)
    case (14)
      out(1:n) = ListGetReal(acc_bc, 'Re Reference Wall Velocity 2', n, acc_nidx, f)
    case (15)
      out(1:n) = ListGetReal(acc_bc, 'Im Reference Wall Velocity 2', n, acc_nidx, f)
    case (16)
      out(1:n) = ListGetReal(acc_bc, 'Re Reference Wall Velocity 3', n, acc_nidx, f)
    case default
      out(1:n) = ListGetReal(acc_bc, 'Im Reference Wall Velocity 3', n, acc_nidx, f)
    end select
    got = merge(1, 0, f)
  end subroutine aca_bc_real

  subroutine aca_bc_constreal(which, v, got) bind(c, name="aca_bc_constreal")
    integer(c_int), value :: which
    real(c_double) :: v
    integer(c_int) :: got
    logical :: f
    select case (which)
    case (1)
      v = ListGetConstReal(acc_bc, 'Moment About 1', f)
    case (2)
      v = ListGetConstReal(acc_bc, 'Moment About 2', f)
    case (3)
      v = ListGetConstReal(acc_bc, 'Moment About 3', f)
    case (4)
      v = ListGetConstReal(acc_bc, 'Tangential Momentum Accommodation Coefficient', f)
    case (5)
      v = ListGetConstReal(acc_bc, 'Energy Accommodation Coefficient', f)
    case default
      v = ListGetConstReal(acc_bc, 'Normal Momentum Accommodation Coefficient', f)
    end select
    got = merge(1, 0, f)
  end subroutine aca_bc_constreal
  ! -- assembly actions --------------------------------------------------------
  subroutine aca_begin_direct(t, n, bodyid) bind(c, name="aca_begin_direct")
    integer(c_int), value :: t
    integer(c_int) :: n, bodyid
    acc_elem => acc_solver % mesh % elements(t)
    acc_model % currentelement => acc_elem
    CurrentModel % currentelement => acc_elem
    acc_nidx => acc_elem % nodeindexes
    acc_nn = acc_elem % type % numberofnodes
    n = acc_nn
    bodyid = acc_elem % bodyid
  end subroutine aca_begin_direct

  function aca_gauss_begin(which, relorder) result(n) bind(c, name="aca_gauss_begin")
    integer(c_int), value :: which, relorder
    integer(c_int) :: n
    if (which == 1) then
      acc_ip = GaussPoints(acc_elem, acc_elem % type % gausspoints2)
    else if (which == 2) then
      acc_ip = GaussPoints(acc_elem, RelOrder=relorder)
    else
      acc_ip = GaussPoints(acc_elem)
    end if
    n = acc_ip % n
  end function aca_gauss_begin

  subroutine aca_gauss_point(t, u, v, w, s) bind(c, name="aca_gauss_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip % u(t)
    v = acc_ip % v(t)
    w = acc_ip % w(t)
    s = acc_ip % s(t)
  end subroutine aca_gauss_point

  subroutine aca_eleminfo(which, u, v, w, n, nd, basis, dbasis, ddbasis, &
      detj, bubbles, stat) bind(c, name="aca_eleminfo")
    integer(c_int), value :: which, n, nd, bubbles
    real(c_double), value :: u, v, w
    real(c_double) :: basis(2*n), dbasis(2*n,3), ddbasis(n,3,3), detj
    integer(c_int) :: stat
    logical :: ok, bub
    bub = (bubbles /= 0)
    if (which == 0) then
      ok = ElementInfo(acc_elem, acc_nodes, u, v, w, detj, &
          basis(1:2*n), dbasis(1:2*n,1:3))
    else
      ok = ElementInfo(acc_elem, acc_nodes, u, v, w, detj, &
          basis(1:2*n), dbasis(1:2*n,1:3), ddbasis(1:n,1:3,1:3), .FALSE., bub)
    end if
    stat = merge(1, 0, ok)
  end subroutine aca_eleminfo

  function aca_elem_area(n) result(a) bind(c, name="aca_elem_area")
    integer(c_int), value :: n
    real(c_double) :: a
    a = ElementArea(acc_solver % mesh, acc_elem, n)
  end function aca_elem_area

  subroutine aca_normal(u, v, nrm) bind(c, name="aca_normal")
    real(c_double), value :: u, v
    real(c_double) :: nrm(3)
    nrm = NormalVector(acc_elem, acc_nodes, u, v, .TRUE.)
  end subroutine aca_normal

  function aca_dot(n, a, b) result(s) bind(c, name="aca_dot")
    integer(c_int), value :: n
    real(c_double) :: a(n), b(n), s
    s = sum(a(1:n) * b(1:n))
  end function aca_dot

  function aca_nodes_dotx(n, basis) result(s) bind(c, name="aca_nodes_dotx")
    integer(c_int), value :: n
    real(c_double) :: basis(n), s
    s = sum(basis(1:n) * acc_nodes % x(1:n))
  end function aca_nodes_dotx

  function aca_sumv(n, a) result(s) bind(c, name="aca_sumv")
    integer(c_int), value :: n
    real(c_double) :: a(n), s
    s = sum(a(1:n))
  end function aca_sumv

  subroutine aca_cmul(ar, ai, br, bi, rr, ri) bind(c, name="aca_cmul")
    real(c_double), value :: ar, ai, br, bi
    real(c_double) :: rr, ri
    complex(kind=dp) :: c
    c = cmplx(ar, ai, kind=dp) * cmplx(br, bi, kind=dp)
    rr = real(c, kind=dp)
    ri = aimag(c)
  end subroutine aca_cmul

  subroutine aca_cdiv(ar, ai, br, bi, rr, ri) bind(c, name="aca_cdiv")
    real(c_double), value :: ar, ai, br, bi
    real(c_double) :: rr, ri
    complex(kind=dp) :: c
    c = cmplx(ar, ai, kind=dp) / cmplx(br, bi, kind=dp)
    rr = real(c, kind=dp)
    ri = aimag(c)
  end subroutine aca_cdiv

  subroutine aca_cmatmul(m, k, n, a, b, c) bind(c, name="aca_cmatmul")
    integer(c_int), value :: m, k, n
    complex(c_double_complex) :: a(m,k), b(k,n), c(m,n)
    c(1:m,1:n) = matmul(a(1:m,1:k), b(1:k,1:n))
  end subroutine aca_cmatmul

  subroutine aca_cinv(n, k) bind(c, name="aca_cinv")
    integer(c_int), value :: n
    complex(c_double_complex) :: k(n,n)
    call ComplexInvertMatrix(k, n)
  end subroutine aca_cinv

  subroutine aca_update_global(lsm, lf, nn, n, dofs) bind(c, name="aca_update_global")
    integer(c_int), value :: nn, n, dofs
    real(c_double) :: lsm(dofs*nn,dofs*nn), lf(dofs*nn)
    call UpdateGlobalEquations(acc_solver % matrix, lsm, &
        acc_solver % matrix % rhs, lf, n, dofs, &
        acc_solver % variable % perm(acc_nidx(1:n)))
  end subroutine aca_update_global

  subroutine aca_setdir_velo(comp, dim) bind(c, name="aca_setdir_velo")
    integer(c_int), value :: comp, dim
    character(len=MAX_NAME_LEN) :: vn
    write(vn,'(A,A,I1)') 'Re Velocity',' ', comp
    call SetDirichletBoundaries(acc_model, acc_solver % matrix, &
        acc_solver % matrix % rhs, vn, (comp-1)*2+1, dim*2+4, &
        acc_solver % variable % perm)
    write(vn,'(A,A,I1)') 'Im Velocity',' ', comp
    call SetDirichletBoundaries(acc_model, acc_solver % matrix, &
        acc_solver % matrix % rhs, vn, (comp-1)*2+2, dim*2+4, &
        acc_solver % variable % perm)
  end subroutine aca_setdir_velo

  subroutine aca_setdir_temp(which, dofs) bind(c, name="aca_setdir_temp")
    integer(c_int), value :: which, dofs
    if (which == 0) then
      call SetDirichletBoundaries(acc_model, acc_solver % matrix, &
          acc_solver % matrix % rhs, 'Re Temperature', dofs-3, dofs, &
          acc_solver % variable % perm)
    else
      call SetDirichletBoundaries(acc_model, acc_solver % matrix, &
          acc_solver % matrix % rhs, 'Im Temperature', dofs-2, dofs, &
          acc_solver % variable % perm)
    end if
  end subroutine aca_setdir_temp

  subroutine aca_zerorow(row) bind(c, name="aca_zerorow")
    integer(c_int), value :: row
    call ZeroRow(acc_solver % matrix, row)
  end subroutine aca_zerorow

  subroutine aca_setmat(row, col, val) bind(c, name="aca_setmat")
    integer(c_int), value :: row, col
    real(c_double), value :: val
    call SetMatrixElement(acc_solver % matrix, row, col, val)
  end subroutine aca_setmat

  subroutine aca_rhs_set(row, val) bind(c, name="aca_rhs_set")
    integer(c_int), value :: row
    real(c_double), value :: val
    acc_solver % matrix % rhs(row) = val
  end subroutine aca_rhs_set

  subroutine aca_default_start() bind(c, name="aca_default_start")
    call DefaultStart()
  end subroutine aca_default_start

  subroutine aca_default_initialize() bind(c, name="aca_default_initialize")
    call DefaultInitialize()
  end subroutine aca_default_initialize

  subroutine aca_default_finish_bulk() bind(c, name="aca_default_finish_bulk")
    call DefaultFinishBulkAssembly()
  end subroutine aca_default_finish_bulk

  subroutine aca_default_dirichlet() bind(c, name="aca_default_dirichlet")
    call DefaultDirichletBCs()
  end subroutine aca_default_dirichlet

  function aca_default_solve() result(norm) bind(c, name="aca_default_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function aca_default_solve

  function aca_cputime() result(t) bind(c, name="aca_cputime")
    real(c_double) :: t
    t = CPUTime()
  end function aca_cputime

  function aca_realtime() result(t) bind(c, name="aca_realtime")
    real(c_double) :: t
    t = RealTime()
  end function aca_realtime
  ! -- formatted output / messages (list-directed and F8.2 writes) ---------------
  subroutine aca_info_blank() bind(c, name="aca_info_blank")
    call Info('AcousticsSolver', ' ', Level=4)
  end subroutine aca_info_blank

  subroutine aca_info_dashes() bind(c, name="aca_info_dashes")
    call Info('AcousticsSolver', '-------------------------------------', Level=4)
  end subroutine aca_info_dashes

  subroutine aca_info_starting() bind(c, name="aca_info_starting")
    call Info('AcousticsSolver', 'Starting Assembly', Level=4)
  end subroutine aca_info_starting

  subroutine aca_info_done() bind(c, name="aca_info_done")
    call Info('AcousticsSolver', 'Assembly done', Level=4)
  end subroutine aca_info_done

  subroutine aca_info_exiting() bind(c, name="aca_info_exiting")
    call Info('AcousticsSolver', 'Exiting the solver...', Level=4)
  end subroutine aca_info_exiting

  subroutine aca_info_freq(val) bind(c, name="aca_info_freq")
    real(c_double), value :: val
    write(acc_message, *) 'Frequency (Hz): ', val
    call Info('AcousticsSolver', acc_message, Level=4)
  end subroutine aca_info_freq

  subroutine aca_info_result(norm) bind(c, name="aca_info_result")
    real(c_double), value :: norm
    write(acc_message, *) 'Result Norm    : ', norm
    call Info('AcousticsSolver', acc_message, Level=4)
  end subroutine aca_info_result

  subroutine aca_info_relchange(rc) bind(c, name="aca_info_relchange")
    real(c_double), value :: rc
    write(acc_message, *) 'Relative Change: ', rc
    call Info('AcousticsSolver', acc_message, Level=4)
  end subroutine aca_info_relchange

  subroutine aca_info_asmtime(at) bind(c, name="aca_info_asmtime")
    real(c_double), value :: at
    write(acc_message,'(a,F8.2)') ' Assembly: (s)', at
    call Info('AcousticsSolver', acc_message, Level=4)
  end subroutine aca_info_asmtime

  subroutine aca_info_soltime(st) bind(c, name="aca_info_soltime")
    real(c_double), value :: st
    write(acc_message,'(a,F8.2)') ' Solve:    (s)', st
    call Info('AcousticsSolver', acc_message, Level=4)
  end subroutine aca_info_soltime

  subroutine aca_info_asm_percent(pct) bind(c, name="aca_info_asm_percent")
    integer(c_int), value :: pct
    write(acc_message,'(a,i3,a)') '   Assembly: ', pct, ' % done'
    call Info('AcousticsSolver', acc_message, Level=5)
  end subroutine aca_info_asm_percent

  subroutine aca_info_asm_est(est) bind(c, name="aca_info_asm_est")
    real(c_double), value :: est
    write(acc_message,'(a,F8.2)') 'Estimated assembly time (s): ', est
    call Info('AcousticsSolver', acc_message, Level=5)
    end subroutine aca_info_asm_est

  subroutine aca_info_block() bind(c, name="aca_info_block")
    call Info('AcousticsSolver', 'Block preconditioning will be used.')
  end subroutine aca_info_block

  subroutine aca_warn_dofs() bind(c, name="aca_warn_dofs")
    call Warn('AcousticsSolver', 'Inconsistent number of Variable Dofs')
  end subroutine aca_warn_dofs

  subroutine aca_fatal_alloc() bind(c, name="aca_fatal_alloc")
    call Fatal('AcousticsSolver', 'Memory allocation error.')
  end subroutine aca_fatal_alloc

  subroutine aca_fatal_bubbles_global() bind(c, name="aca_fatal_bubbles_global")
    write(acc_message, *) 'The keyword Bubbles in Global System should have the value False'
    call Fatal('AcousticsSolver', acc_message)
  end subroutine aca_fatal_bubbles_global

  subroutine aca_error_minibubbles(code) bind(c, name="aca_error_minibubbles")
    integer(c_int), value :: code
    write(acc_message, '(a,i4,a)') 'Mini-bubbles for element ', code, ' are not implemented.'
    call Error('Acoustics', acc_message)
  end subroutine aca_error_minibubbles

  subroutine aca_error_2ndderiv(code) bind(c, name="aca_error_2ndderiv")
    integer(c_int), value :: code
    write(acc_message, '(a,i4,a)') 'Second derivatives for element ', code, ' are not implemented.'
    call Error('Acoustics', acc_message)
  end subroutine aca_error_2ndderiv

  subroutine aca_list_add_freq(val) bind(c, name="aca_list_add_freq")
    real(c_double), value :: val
    call ListAddConstReal(acc_model % simulation, 'res: frequency', val)
  end subroutine aca_list_add_freq

  ! -- per-node complex rescale statements (verbatim arithmetic) -----------------
  subroutine aca_pre_rescale(hr, sh, dens, temp, visc, lam, af, &
      f1, f2, f3, f4, g1, g2, g3, g4) bind(c, name="aca_pre_rescale")
    real(c_double), value :: hr, sh, dens, temp, visc, lam, af, f1, f2, f3, f4
    real(c_double) :: g1, g2, g3, g4
    complex(kind=dp) :: a
    a = (hr-1.0d0) * sh / af * cmplx(f1, f2, kind=dp)
    g1 = real(a, kind=dp)
    g2 = aimag(a)
    a = cmplx(1.0d0, lam*af/((hr-1.0d0)*sh*dens*temp), kind=dp) * &
        (1.0d0/(af*dens) * cmplx(f3, f4, kind=dp) - cmplx(f1, f2, kind=dp))
    g3 = real(a, kind=dp)
    g4 = aimag(a)
  end subroutine aca_pre_rescale

  subroutine aca_post_rescale(dens, temp, lam, press, af, &
      f1, f2, f3, f4, g1, g2, g3, g4) bind(c, name="aca_post_rescale")
    real(c_double), value :: dens, temp, lam, press, af, f1, f2, f3, f4
    real(c_double) :: g1, g2, g3, g4
    complex(kind=dp) :: a

    a = af * dens * (cmplx(f1, f2, kind=dp) + &
        1.0d0/cmplx(1.0d0, lam*af/press, kind=dp) * cmplx(f3, f4, kind=dp))
    g3 = real(a, kind=dp)
    g4 = aimag(a)
    a = af * dens * temp/press * cmplx(f1, f2, kind=dp)
    g1 = real(a, kind=dp)
    g2 = aimag(a)
  end subroutine aca_post_rescale
  ! -- BlockPreconditioning iterative suite (verbatim preservation) ----------------
  ! The outer/inner GCR and BiCGStab(l) iterations below are the original
  ! internal procedures of AcousticsSolver, preserved statement for
  ! statement as internals of this container. Host association (`Solver`,
  ! `Message`, `UtilizePreviousSolution`) is re-established by the
  ! container from explicit arguments; SAVE lifetimes are unchanged.
  subroutine aca_block_solve(which, maxit, tol, cdim, jdeg, utilprev, norm) &
      bind(c, name="aca_block_solve")
    integer(c_int), value :: which, maxit, cdim, jdeg, utilprev
    real(c_double), value :: tol
    real(c_double) :: norm
    type(solver_t), pointer :: solver
    type(matrix_t), pointer :: a, pa, ps, m, l
    real(kind=dp), pointer :: px(:), pb(:)
    integer :: nrows, mrows
    real(kind=dp), allocatable :: xt(:), bt(:)
    character(len=MAX_NAME_LEN) :: message
    logical :: UtilizePreviousSolution
    solver => acc_solver
    a => acc_solver % matrix
    pa => acc_amatrix
    ps => acc_smatrix
    m => acc_mmatrix
    l => acc_lmatrix
    nrows = a % numberofrows
    allocate(xt(nrows), bt(nrows))
    px => acc_solver % variable % values
    pb => a % rhs
    xt = px(1:nrows)
    bt = pb(1:nrows)
    UtilizePreviousSolution = (utilprev /= 0)
    mrows = m % numberofrows
    select case (which)
    case (1)
      call GCROuterIteration(nrows, a, mrows, pa, ps, m, l, xt, bt, &
          maxit, tol, cdim, norm)
    case (2)
      call BiCGStabOuterIteration(nrows, a, mrows, pa, ps, m, l, xt, bt, &
          maxit, tol, cdim, norm)
    case (3)
      call BiCGStablOuterIteration(jdeg, nrows, a, mrows, pa, ps, m, l, xt, bt, &
          maxit, tol, cdim, norm)
    case default
      call InnerOuterIteration(nrows, a, mrows, pa, ps, m, l, xt, bt, maxit, &
          tol, cdim, norm)
    end select
    px(1:nrows) = xt
    pb(1:nrows) = bt
    deallocate(xt, bt)
  contains
!------------------------------------------------------------------------------
  SUBROUTINE OptimalScaling( n, A, s )
!------------------------------------------------------------------------------
!   This subroutine equilibrates the rows of the complex-coefficient matrix A
!   to minimize the condition number. Only the coefficients necessary for obtaining
!   Re(A*z), with z a complex vector, are modified. The vector s will contain
!   the complex scaling factors.
!-----------------------------------------------------------------------------
    INTEGER :: n
    TYPE(Matrix_t), POINTER :: A
    COMPLEX(kind=dp) :: s(n)
!-----------------------------------------------------------------------------
    INTEGER :: i, j
    REAL(kind=dp) :: norm, tmp
    INTEGER, POINTER :: Cols(:), Rows(:)
    REAL(KIND=dp), POINTER :: Values(:)
    !-------------------------------------------------------------------------
    Rows   => A % Rows
    Cols   => A % Cols
    Values => A % Values

    norm = 0.0d0
    DO i=1,n
      tmp = 0.0d0
      DO j=Rows(2*i-1),Rows(2*i)-1,2
        tmp = tmp + CDABS( CMPLX( Values(j), -Values(j+1), kind=dp ) )
      END DO
      IF (tmp > norm) norm = tmp
      s(i) = CMPLX( 1.0d0,0.0d0, kind=dp) / tmp
    END DO

!------------------------------------------------------------------------------
    END SUBROUTINE OptimalScaling
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE OptimalMatrixScaling( n, A, s )
!------------------------------------------------------------------------------
!   This subroutine equilibrates the rows of the complex-coefficient matrix A
!   to minimize the condition number. Only the coefficients necessary for obtaining
!   Re(A*z), with z a complex vector, are modified. The vector s will contain
!   the complex scaling factors.
!-----------------------------------------------------------------------------
    INTEGER :: n
    TYPE(Matrix_t), POINTER :: A
    COMPLEX(kind=dp) :: s(n)
!-----------------------------------------------------------------------------
    INTEGER :: i, j
    REAL(kind=dp) :: norm, tmp
    INTEGER, POINTER :: Cols(:), Rows(:)
    REAL(KIND=dp), POINTER :: Values(:)
    !-------------------------------------------------------------------------
    Rows   => A % Rows
    Cols   => A % Cols
    Values => A % Values

    norm = 0.0d0
    DO i=1,n
      tmp = 0.0d0
      DO j=Rows(2*i-1),Rows(2*i)-1,2
        tmp = tmp + CDABS( CMPLX( Values(j), -Values(j+1), kind=dp ) )
      END DO
      IF (tmp > norm) norm = tmp
      s(i) = CMPLX( 1.0d0,0.0d0, kind=dp) / tmp
      DO j=Rows(2*i-1),Rows(2*i)-1,2
        Values(j) = Values(j) * s(i)
        Values(j+1) = Values(j+1) * s(i)
      END DO
    END DO
    !--------------------------------------------------------------
    WRITE( Message, * ) 'Unscaled matrix norm: ', norm
    CALL Info( 'AcousticsSolver', Message, Level=5 )
    !--------------------------------------------------------------
    norm = 0.0d0
    DO i=1,n
      tmp = 0.0d0
      DO j=Rows(2*i-1),Rows(2*i)-1,2
        tmp = tmp + CDABS( CMPLX( Values(j), -Values(j+1), kind=dp ) )
      END DO
      IF (tmp > norm) norm = tmp
    END DO
    !---------------------------------------------------------------
    WRITE( Message, * ) 'Scaled matrix norm: ', norm
    CALL Info( 'AcousticsSolver', Message, Level=5 )
    !---------------------------------------------------------------

    ! recovering original system
    !DO i=1,n
    !  b(i) = b(i)/s(i)
    !  DO j=Rows(2*i-1),Rows(2*i)-1,2
    !    Values(j) = Values(j)/s(i)
    !    Values(j+1) = Values(j+1)/s(i)
    !  END DO
    !END DO


!------------------------------------------------------------------------------
    END SUBROUTINE OptimalMatrixScaling
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  SUBROUTINE GCROuterIteration( n, A, q, PA, PS, MMatrix, LMatrix, &
      x, b, Rounds, TOL, dim, Norm )
!------------------------------------------------------------------------------
!    This is the preconditioned GCR iteration for the complex linear system Ax=b.
!    The Schur complement preconditioning strategy is used.
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: A, PA, PS, MMatrix, LMatrix
    INTEGER :: n, q, Rounds, dim
    REAL(KIND=dp), TARGET :: x(n), b(n)
    REAL(KIND=dp) :: TOL, Norm
!-------------------------------------------------------------------------------
    INTEGER :: i, j, k, m, InnerRounds, IluOrder, VelocityPrecond = 1, &
        MaxRestarts
    LOGICAL :: Condition, GotIt, SystemScaling, ConvergedSol = .FALSE.
    REAL(KIND=dp) :: res, tottime, res0, const, stime, InnerTol, alpha
    COMPLEX(KIND=dp) :: r(n/2),T1(n/2),T2(n/2), &
        S(n/2,Rounds), V(n/2,Rounds), e(n/2), &
        Vel(dim*q/2), VelRhs(dim*q/2), da(n/2), dps(q), dpa(dim*q/2)
    ! x and b already hold the complex vectors as consecutive (Re,Im) pairs, so
    ! alias them rather than copying into complex temporaries and back.
    COMPLEX(KIND=dp), POINTER :: y(:),f(:)
    COMPLEX(KIND=dp) :: beta
!------------------------------------------------------------------------------

    tottime = CPUTime()
    !-----------------------------------------------------------------------------
    ! Scaling to minimize the condition number
    !-----------------------------------------------------------------------------
    SystemScaling = ListGetLogical( Solver % Values, 'Linear System Scaling', GotIt)
    IF (.NOT. GotIt) SystemScaling = .FALSE.
    IF (SystemScaling) THEN
      CALL OptimalMatrixScaling( n/2, A, da )
      CALL OptimalMatrixScaling( q, PS, dps )
      CALL OptimalMatrixScaling( dim*q/2, PA, dpa )
    END IF

    !-------------------------------------------------------------------------------
    ! Compute ILU factorizations for the preconditioner matrices.
    ! This needs to be done only once.
    !-------------------------------------------------------------------------------
    IluOrder = ListGetInteger( Solver % Values, 'ILU Order for Schur Complement')

    CALL Info( 'AcousticsSolver', ' ', Level=4)
    CALL Info( 'AcousticsSolver', 'ILU factorization for the Schur complement preconditioner', &
        Level=4)
    CALL Info( 'AcousticsSolver', ' ', Level=4)
    Condition = CRS_ComplexIncompleteLU( PS, IluOrder )

    IluOrder = ListGetInteger( Solver % Values, 'ILU Order for Velocities', GotIt)

    IF (GotIt) THEN
      CALL Info( 'AcousticsSolver', ' ', Level=4)
      CALL Info( 'AcousticsSolver', 'ILU factorization for the Velocity preconditioner', &
          Level=4)
      CALL Info( 'AcousticsSolver', ' ', Level=4)
      Condition = CRS_ComplexIncompleteLU( PA, IluOrder )
    ELSE
      VelocityPrecond = 0
    END IF

    !--------------------------------------------------------------------------------
    !   Some initializations
    !--------------------------------------------------------------------------------
    InnerTol = ListGetConstReal( Solver % Values, &
        'Linear System Convergence Tolerance' )
    InnerRounds = ListGetInteger( Solver % Values, &
        'Linear System Max Iterations')
    MaxRestarts = ListGetInteger( Solver % Values, &
        'Max GCR Restarts', GotIt)
    IF (.NOT. GotIt) MaxRestarts = 1

    !--------------------------------------------------------------------------------------------
    ! The solution of an initial guess, the previous solution cannot be used as an initial guess
    !--------------------------------------------------------------------------------------------
    x(1:n) = 0.0d0

    Vel(1:dim*q/2) = CMPLX( 0.0d0,0.0d0, kind=dp )
    DO i=1,dim
      DO j=1,q/2
        VelRhs((j-1)*dim+i) = CMPLX( b(2*(dim+2)*(j-1)+2*i-1), &
            b(2*(dim+2)*(j-1)+2*i), kind=dp )
      END DO
    END DO

    IF ( ANY( VelRhs /= CMPLX(0.0d0, 0.0d0, kind=dp) ) ) THEN
      !-----------------------------------
      ! Scale the right-hand side vector
      !-----------------------------------
      IF (SystemScaling) THEN
        DO i=1,dim*q/2
          VelRhs(i) = dpa(i) * VelRhs(i)
        END DO
      END IF
      !-----------------------------------------------------------------
      CALL Info( 'AcousticsSolver', ' ', Level=4)
      CALL Info( 'AcousticsSolver', &
          'Solving initial guess for velocities', Level=4)
      CALL Info( 'AcousticsSolver', ' ', Level=4)
      !------------------------------------------------------------------
      CALL ComplexBiCGStab( q*dim, PA, Vel, VelRhs, InnerRounds, 1.0d-7, 0 )

      DO j=1,q/2
        DO i=1,dim
          x(2*(dim+2)*(j-1)+2*i-1) = REAL(Vel((j-1)*dim+i))
          x(2*(dim+2)*(j-1)+2*i) = AIMAG(Vel((j-1)*dim+i))
        END DO
      END DO
    END IF

    !--------------------------------------------------------------------------------
    ! The start of the GCR iteration...
    !--------------------------------------------------------------------------------
    m = n/2
    !----------------------------------------------------------------------------
    ! View the solution vector x and the right-hand side vector b as the
    ! complex-valued vectors y and f. Nothing is copied, so the solution needs
    ! no transforming back either.
    !---------------------------------------------------------------------------
    y => ComplexValues( x, m )
    f => ComplexValues( b, m )

    CALL ComplexMatrixVectorProduct( A, y, r )
    !--------------------------------------------------------
    ! Scale the right-hand side vector and form the residual
    !--------------------------------------------------------
    IF (SystemScaling) THEN
      DO i=1,m
        f(i) = da(i) * f(i)
      END DO
    END IF
    r(1:m) = f(1:m) - r(1:m)
    res0 = ComplexNorm(m,f)
    res = ComplexNorm(m,r)/res0

    ! PRINT *,'OuterIteration ',0, ComplexNorm(m,r)/res0, StoppingCriterion(m,A,y,f,r), &
    !     CPUTime() - tottime
    WRITE(*,'(a,I4,ES12.3,ES12.3)') 'OuterIteration residual for iterate', &
        0, res, CPUTime() - tottime

    DO j=1,MaxRestarts
      IF (ConvergedSol) EXIT
      V(1:m,1:Rounds) = CMPLX( 0.0d0, 0.0d0, kind=dp)
      S(1:m,1:Rounds) = CMPLX( 0.0d0, 0.0d0, kind=dp)
      DO k=1,Rounds
        !----------------------------------------------------------
        ! Perform the preconditioning...
        !---------------------------------------------------------------
        T1(1:m) = r(1:m)
        CALL PreconditioningIteration(n, A, q, PA, PS, T1, dim, &
            da, dpa, dps, SystemScaling, MMatrix, LMatrix)
        CALL ComplexMatrixVectorProduct( A, T1, T2 )

        !--------------------------------------------------------------
        ! Perform the orthogonalisation of the search directions...
        !--------------------------------------------------------------
        DO i=1,k-1
          beta = ComplexDotProduct( m, V(1:m,i), T2(1:m) )
          T1(1:m) = T1(1:m) - beta * S(1:m,i)
          T2(1:m) = T2(1:m) - beta * V(1:m,i)
        END DO
        alpha = ComplexNorm(m,T2)
        T1(1:m) = CMPLX( 1.0d0, 0.0d0, kind=dp)/CMPLX( alpha, 0.0d0, kind=dp) * T1(1:m)
        T2(1:m) = CMPLX( 1.0d0, 0.0d0, kind=dp)/CMPLX( alpha, 0.0d0, kind=dp) * T2(1:m)

        !-------------------------------------------------------------
        ! The update of the solution and save the search data...
        !-------------------------------------------------------------
        beta = ComplexDotProduct(m, T2, r)
        y(1:m) = y(1:m) + beta * T1(1:m)
        r(1:m) = r(1:m) - beta * T2(1:m)
        S(1:m,k) = T1(1:m)
        V(1:m,k) = T2(1:m)

        !----------------------------------------------------------------
        ! Check the accuracy of the residual, if desired...
        !----------------------------------------------------------------
        IF (.FALSE.) THEN
          e(1:m) = CMPLX( 0.0d0,0.0d0, kind=dp )
          CALL ComplexMatrixVectorProduct( A, y, e )
          e(1:m) = f(1:m) - e(1:m)
          norm = ComplexNorm(m,e(1:m)-r(1:m))/ComplexNorm(m,e(1:m))
          PRINT *, 'Relative error of the residual: ', norm
        END IF

        !----------------------------------------------------
        ! Check whether the convergence criterion is met
        !----------------------------------------------------
        res = ComplexNorm(m,r)/res0

        !res = StoppingCriterion( m, A, y, f, r )
        !PRINT *,'OuterIteration ',i,res, StoppingCriterion(m,A,y,f,r),CPUTime() - tottime

        WRITE(*,'(a,I4,ES12.3,ES12.3)') 'OuterIteration residual for iterate', &
            k + (j-1) * Rounds, res, CPUTime() - tottime
        ConvergedSol = ( res < TOL)
        IF (ConvergedSol) EXIT
      END DO
    END DO

    WRITE(*,'(a,ES12.3)') 'An approximate lower bound for the condition number: ', &
        ConditionEstimate( m, A, y, f)

    !----------------------------------------------
    ! No transforming back: y aliases x.
    !----------------------------------------------
    Norm = SQRT(DOT_PRODUCT( x(1:2*m), x(1:2*m) )/(2*m))

!------------------------------------------------------------------------------
  END SUBROUTINE GCROuterIteration
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  SUBROUTINE BiCGStabOuterIteration( n, A, q, PA, PS, MMatrix, LMatrix, &
      x, b, Rounds, TOL, dim, Norm )
!------------------------------------------------------------------------------
!   This is the preconditioned BiCGStab iteration for the complex linear
!   system Ax=b. The Schur complement preconditioning strategy is used.
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: A, PA, PS, MMatrix, LMatrix
    INTEGER :: n, q, Rounds, dim
    REAL(KIND=dp), TARGET :: x(n), b(n)
    REAL(KIND=dp) :: TOL, Norm
!-------------------------------------------------------------------------------
    INTEGER :: i, j, m, InnerRounds, IluOrder, VelocityPrecond = 1
    LOGICAL :: Condition, GotIt, SystemScaling
    REAL(KIND=dp) :: res, tottime, res0, InnerTol
    COMPLEX(KIND=dp) :: r(n/2),Ri(n/2),P(n/2),V(n/2),T(n/2),T1(n/2),T2(n/2),S(n/2), &
        Vel(dim*q/2), VelRhs(dim*q/2), da(n/2), dps(q), dpa(dim*q/2)
    ! x and b already hold the complex vectors as consecutive (Re,Im) pairs, so
    ! alias them rather than copying into complex temporaries and back.
    COMPLEX(KIND=dp), POINTER :: y(:),f(:)
    COMPLEX(KIND=dp) :: alpha, beta, omega, rho, oldrho
!------------------------------------------------------------------------------

    tottime = CPUTime()
    !-----------------------------------------------------------------------------
    ! Scaling to minimize the condition number
    !-----------------------------------------------------------------------------
    SystemScaling = ListGetLogical( Solver % Values, 'Linear System Scaling', GotIt)
    IF (.NOT. GotIt) SystemScaling = .FALSE.
    IF (SystemScaling) THEN
      CALL OptimalMatrixScaling( n/2, A, da )
      CALL OptimalMatrixScaling( q, PS, dps )
      CALL OptimalMatrixScaling( dim*q/2, PA, dpa )
    END IF

    !-------------------------------------------------------------------------------
    ! Compute ILU factorizations for the preconditioner matrices.
    ! This needs to be done only once.
    !-------------------------------------------------------------------------------
    IluOrder = ListGetInteger( Solver % Values, 'ILU Order for Schur Complement')

    CALL Info( 'AcousticsSolver', ' ', Level=4)
    CALL Info( 'AcousticsSolver', 'ILU factorization for the Schur complement preconditioner', &
        Level=4)
    CALL Info( 'AcousticsSolver', ' ', Level=4)
    Condition = CRS_ComplexIncompleteLU( PS, IluOrder )

    IluOrder = ListGetInteger( Solver % Values, 'ILU Order for Velocities', GotIt)

    IF (GotIt) THEN
      CALL Info( 'AcousticsSolver', ' ', Level=4)
      CALL Info( 'AcousticsSolver', 'ILU factorization for the Velocity preconditioner', &
          Level=4)
      CALL Info( 'AcousticsSolver', ' ', Level=4)
      Condition = CRS_ComplexIncompleteLU( PA, IluOrder )
    ELSE
      VelocityPrecond = 0
    END IF

    !--------------------------------------------------------------------------------
    !   Some initializations
    !--------------------------------------------------------------------------------
    InnerTol = ListGetConstReal( Solver % Values, &
        'Linear System Convergence Tolerance' )
    InnerRounds = ListGetInteger( Solver % Values, &
        'Linear System Max Iterations')

    !--------------------------------------------------------------------------------------------
    ! The solution of an initial guess, the previous solution cannot be used as an initial guess
    !--------------------------------------------------------------------------------------------
    x(1:n) = 0.0d0

    Vel(1:dim*q/2) = CMPLX( 0.0d0,0.0d0, kind=dp )
    DO i=1,dim
      DO j=1,q/2
        VelRhs((j-1)*dim+i) = CMPLX( b(2*(dim+2)*(j-1)+2*i-1), &
            b(2*(dim+2)*(j-1)+2*i), kind=dp )
      END DO
    END DO

    IF ( ANY( VelRhs /= CMPLX(0.0d0, 0.0d0, kind=dp) ) ) THEN
      !-----------------------------------
      ! Scale the right-hand side vector
      !-----------------------------------
      IF (SystemScaling) THEN
        DO i=1,dim*q/2
          VelRhs(i) = dpa(i) * VelRhs(i)
        END DO
      END IF
      !-----------------------------------------------------------------
      CALL Info( 'AcousticsSolver', ' ', Level=4)
      CALL Info( 'AcousticsSolver', &
          'Solving initial guess for velocities', Level=4)
      CALL Info( 'AcousticsSolver', ' ', Level=4)
      !------------------------------------------------------------------
      CALL ComplexBiCGStab( q*dim, PA, Vel, VelRhs, InnerRounds, 1.0d-7, 0 )

      DO j=1,q/2
        DO i=1,dim
          x(2*(dim+2)*(j-1)+2*i-1) = REAL(Vel((j-1)*dim+i))
          x(2*(dim+2)*(j-1)+2*i) = AIMAG(Vel((j-1)*dim+i))
        END DO
      END DO
    END IF

    !--------------------------------------------------------------------------------
    ! The start of the BiCGStab iteration...
    !--------------------------------------------------------------------------------
    m = n/2
    !----------------------------------------------------------------------------
    ! View the solution vector x and the right-hand side vector b as the
    ! complex-valued vectors y and f. Nothing is copied, so the solution needs
    ! no transforming back either.
    !---------------------------------------------------------------------------
    y => ComplexValues( x, m )
    f => ComplexValues( b, m )

    CALL ComplexMatrixVectorProduct( A, y, r )
    !--------------------------------------------------------
    ! Scale the right-hand side vector and form the residual
    !--------------------------------------------------------
    IF (SystemScaling) THEN
      DO i=1,m
        f(i) = da(i) * f(i)
      END DO
    END IF
    r(1:m) = f(1:m) - r(1:m)
    res0 = ComplexNorm(m,f)
    res = ComplexNorm(m,r)/res0

    ! PRINT *,'OuterIteration ',0, ComplexNorm(m,r)/res0, StoppingCriterion(m,A,y,f,r), &
    !     CPUTime() - tottime

    WRITE(*,'(a,I4,ES12.3,ES12.3)') 'OuterIteration residual for iterate', &
        0, res, CPUTime() - tottime

    Ri(1:m) = r(1:m)
    P(1:m) = CMPLX( 0.0d0, 0.0d0, kind=dp)
    V(1:m) = CMPLX( 0.0d0, 0.0d0, kind=dp)
    omega  = CMPLX( 1.0d0, 0.0d0, kind=dp)
    alpha  = CMPLX( 0.0d0, 0.0d0, kind=dp)
    oldrho = CMPLX( 1.0d0, 0.0d0, kind=dp)

    DO i=1,Rounds
      rho = ComplexDotProduct( m, Ri, r )
      beta = alpha * rho / ( oldrho * omega )
      P(1:m) = r(1:m) + beta * (P(1:m) - omega*V(1:m))
      V(1:m) = P(1:m)
      !----------------------------------------------------------
      ! Perform the preconditioning...
      !---------------------------------------------------------------
      CALL PreconditioningIteration(n, A, q, PA, PS, V, dim, &
          da, dpa, dps, SystemScaling, MMatrix, LMatrix)
      !--------------------------------------------------------------
      T1(1:m) = V(1:m)
      CALL ComplexMatrixVectorProduct( A, T1, V )
      alpha = rho / ComplexDotProduct( m, Ri, V )
      S(1:m) = r(1:m) - alpha * V(1:m)

      !---------------------------------------------------------------------------------
      ! The update of the solution and the computation of the residual-based error indicator
      !---------------------------------------------------------------------------------
      y(1:m) = y(1:m) + alpha*T1(1:m)

      res = ComplexNorm(m,S)/res0
      !res = StoppingCriterion( m, A, y, f, S )
      IF ( res < TOL ) THEN
        WRITE(*,'(a,I4,ES12.3,ES12.3)') 'OuterIteration residual for iterate', i, res, &
            StoppingCriterion( m, A, y, f, S )
        WRITE(*,'(a,ES12.3)') 'An approximate lower bound for the condition number: ', &
            ConditionEstimate( m, A, y, f)
        EXIT
      END IF

      T(1:m) = S(1:m)
      !----------------------------------------------------------
      ! Perform the preconditioning...
      !-----------------------------------------------------------------
      CALL PreconditioningIteration(n, A, q, PA, PS, T, dim, &
          da, dpa, dps, SystemScaling, MMatrix, LMatrix)
      !-----------------------------------------------------------------
      T2(1:m) = T(1:m)
      CALL ComplexMatrixVectorProduct( A, T2, T )
      omega = ComplexDotProduct( m,T,S ) / ComplexDotProduct( m,T,T )
      oldrho = rho
      r(1:m) = S(1:m) - omega*T(1:m)
      y(1:m) = y(1:m) + omega*T2(1:m)

      res = ComplexNorm(m,r)/res0
      ! res = StoppingCriterion( m, A, y, f, r )

      WRITE(*,'(a,I4,ES12.3,ES12.3)') 'OuterIteration residual for iterate', i, res, &
          StoppingCriterion( m, A, y, f, r )

      IF ( res < TOL .OR. i==Rounds) THEN
        WRITE(*,'(a,ES12.3)') 'An approximate lower bound for the condition number: ', &
            ConditionEstimate( m, A, y, f)
        EXIT
      END IF
    END DO

    !----------------------------------------------
    ! No transforming back: y aliases x.
    !----------------------------------------------
    Norm = SQRT(DOT_PRODUCT( x(1:2*m), x(1:2*m) )/(2*m))


!------------------------------------------------------------------------------
  END SUBROUTINE BiCGStabOuterIteration
!------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
  SUBROUTINE BiCGStablOuterIteration( l, m, A, q, PA, PS, MMatrix, LMatrix, &
      v, f, MaxRounds, Tol, dim, Norm)
!-----------------------------------------------------------------------------------
!  This is the preconditioned BiCGStab(l) iteration for the complex linear
!  system Ax=b. The Schur complement preconditioning strategy is used.
!------------------------------------------------------------------------------
   INTEGER :: l, m, q, MaxRounds, dim
   TYPE(Matrix_t), POINTER :: A, PA, PS, MMatrix, LMatrix
   REAL(KIND=dp), TARGET :: v(m), f(m)
   REAL(KIND=dp) :: Tol, Norm
!-----------------------------------------------------------------------------------
   LOGICAL :: SystemScaling, GotIt, Condition
   COMPLEX(KIND=dp) :: da(m/2), dps(q), dpa(dim*q/2), Vel(dim*q/2), VelRhs(dim*q/2)
   ! v and f already hold the complex vectors as consecutive (Re,Im) pairs, so
   ! alias them rather than copying into complex temporaries and back.
   COMPLEX(KIND=dp), POINTER :: x(:), b(:)
   INTEGER :: IluOrder, InnerRounds, VelocityPrecond = 1
   REAL(KIND=dp) :: InnerTol, tottime

   COMPLEX(KIND=dp) :: zzero, zone, t(m/2), kappa0, kappal
   REAL(KIND=dp) :: rnrm0, bnrm, rnrm, mxnrmx, mxnrmr, errorind, &
       delta = 1.0d-2
   INTEGER :: i, j, rr, r, u, xp, bp, z, zz, y0, yl, y, k, iwork(l-1), stat, Round
   INTEGER :: n
   COMPLEX(KIND=dp) :: work(m/2,3+2*(l+1)), rwork(l+1,3+2*(l+1)), &
       alpha, beta, omega, rho0, rho1, sigma, varrho, hatgamma
   LOGICAL rcmp, xpdt

!------------------------------------------------------------------------------

   tottime = CPUTime()
   !-----------------------------------------------------------------------------
   ! Scaling to minimize the condition number
   !-----------------------------------------------------------------------------
   SystemScaling = ListGetLogical( Solver % Values, 'Linear System Scaling', GotIt)
   IF (.NOT. GotIt) SystemScaling = .FALSE.
   IF (SystemScaling) THEN
     CALL OptimalMatrixScaling( m/2, A, da )
     CALL OptimalMatrixScaling( q, PS, dps )
     CALL OptimalMatrixScaling( dim*q/2, PA, dpa )
   END IF

   !-------------------------------------------------------------------------------
   ! Compute ILU factorizations for the preconditioner matrices.
   ! This needs to be done only once.
   !-------------------------------------------------------------------------------
   IluOrder = ListGetInteger( Solver % Values, 'ILU Order for Schur Complement')

   CALL Info( 'AcousticsSolver', ' ', Level=4)
   CALL Info( 'AcousticsSolver', 'ILU factorization for the Schur complement preconditioner', &
       Level=4)
   CALL Info( 'AcousticsSolver', ' ', Level=4)
   Condition = CRS_ComplexIncompleteLU( PS, IluOrder )

   IluOrder = ListGetInteger( Solver % Values, 'ILU Order for Velocities', GotIt)

   IF (GotIt) THEN
     CALL Info( 'AcousticsSolver', ' ', Level=4)
     CALL Info( 'AcousticsSolver', 'ILU factorization for the Velocity preconditioner', &
         Level=4)
     CALL Info( 'AcousticsSolver', ' ', Level=4)
     Condition = CRS_ComplexIncompleteLU( PA, IluOrder )
   ELSE
     VelocityPrecond = 0
   END IF

   !--------------------------------------------------------------------------------
   !   Some initializations
   !--------------------------------------------------------------------------------
   InnerTol = ListGetConstReal( Solver % Values, &
       'Linear System Convergence Tolerance' )
   InnerRounds = ListGetInteger( Solver % Values, &
       'Linear System Max Iterations')

   !--------------------------------------------------------------------------------------------
   ! The solution of an initial guess, the previous solution cannot be used as an initial guess
   !--------------------------------------------------------------------------------------------
   v(1:m) = 0.0d0

   Vel(1:dim*q/2) = CMPLX( 0.0d0,0.0d0, kind=dp )
   DO i=1,dim
     DO j=1,q/2
       VelRhs((j-1)*dim+i) = CMPLX( f(2*(dim+2)*(j-1)+2*i-1), &
           f(2*(dim+2)*(j-1)+2*i), kind=dp )
     END DO
   END DO

   IF ( ANY( VelRhs /= CMPLX(0.0d0, 0.0d0, kind=dp) ) ) THEN
     !-----------------------------------
     ! Scale the right-hand side vector
     !-----------------------------------
     IF (SystemScaling) THEN
       DO i=1,dim*q/2
         VelRhs(i) = dpa(i) * VelRhs(i)
       END DO
     END IF
     !-----------------------------------------------------------------
     CALL Info( 'AcousticsSolver', ' ', Level=4)
     CALL Info( 'AcousticsSolver', &
         'Solving initial guess for velocities', Level=4)
     CALL Info( 'AcousticsSolver', ' ', Level=4)
     !------------------------------------------------------------------
     CALL ComplexBiCGStab( q*dim, PA, Vel, VelRhs, InnerRounds, 1.0d-7, 0 )

     DO j=1,q/2
       DO i=1,dim
         v(2*(dim+2)*(j-1)+2*i-1) = REAL(Vel((j-1)*dim+i))
         v(2*(dim+2)*(j-1)+2*i) = AIMAG(Vel((j-1)*dim+i))
       END DO
     END DO
   END IF

   !--------------------------------------------------------------------------------
   ! The start of the BiCGStabl iteration...
   !--------------------------------------------------------------------------------
   n = m/2
   !----------------------------------------------------------------------------
   ! View the solution vector v and the right-hand side vector f as the
   ! complex-valued vectors x and b. Nothing is copied, so the solution needs
   ! no transforming back either.
   !---------------------------------------------------------------------------
   x => ComplexValues( v, n )
   b => ComplexValues( f, n )

   zzero = CMPLX( 0.0d0,0.0d0, kind=dp)
   zone =  CMPLX( 1.0d0,0.0d0, kind=dp)
   work = CMPLX( 0.0d0, 0.0d0, kind=dp )
   rwork = CMPLX( 0.0d0, 0.0d0, kind=dp )

   rr = 1
   r = rr+1
   u = r+(l+1)
   xp = u+(l+1)
   bp = xp+1

   z = 1
   zz = z+(l+1)
   y0 = zz+(l+1)
   yl = y0+1
   y = yl+1

   CALL ComplexMatrixVectorProduct( A, x, work(1:n,r) )
   work(1:n,r) = b(1:n) - work(1:n,r)
   bnrm = dznrm2(n, b(1:n), 1)

   work(1:n,rr) = work(1:n,r)
   work(1:n,bp) = work(1:n,r)
   work(1:n,xp) = x(1:n)
   x(1:n) = zzero
   rnrm0 = dznrm2(n, work(1:n,r), 1)
   rnrm = rnrm0
   mxnrmx = rnrm0
   mxnrmr = rnrm0

   alpha = zzero
   omega = zone
   sigma = zone
   rho0 = zone

   Round = 0
   errorind = rnrm/bnrm

   WRITE(*,'(a,I4,ES12.3)') 'OuterIteration residual for iterate', &
        0, errorind

   DO WHILE ( errorind > Tol .AND. Round < MaxRounds)
     Round = Round + 1
     !-------------------------
     ! --- The BiCG part ---
     !-------------------------
     rho0 = -omega*rho0

     DO k=1,l
       rho1 = zdotc(n, work(1:n,rr), 1, work(1:n,r+k-1), 1)
       IF (rho0 == zzero) THEN
         CALL Fatal( 'ComplexBiCGStab(l)', 'Breakdown error.' )
       ENDIF
       beta = alpha*(rho1/rho0)
       rho0 = rho1
       DO j=0,k-1
         work(1:n,u+j) = work(1:n,r+j) - beta*work(1:n,u+j)
       ENDDO
       t(1:n) = work(1:n,u+k-1)
       !----------------------------------------------------------
       ! Perform the preconditioning...
       !---------------------------------------------------------------
       CALL PreconditioningIteration(m, A, q, PA, PS, t, dim, &
          da, dpa, dps, SystemScaling, MMatrix, LMatrix)
       !--------------------------------------------------------------

       CALL ComplexMatrixVectorProduct( A, t, work(1:n,u+k) )
       sigma = zdotc(n, work(1:n,rr), 1, work(1:n,u+k), 1)
       IF (sigma == zzero) THEN
         CALL Fatal( 'ComplexBiCGStab(l)', 'Breakdown error.' )
       ENDIF
       alpha = rho1/sigma
       x(1:n) = x(1:n) + alpha * work(1:n,u)
       DO j=0,k-1
         work(1:n,r+j) = work(1:n,r+j) - alpha * work(1:n,u+j+1)
       ENDDO
       t(1:n) = work(1:n,r+k-1)
       !-----------------------------------------------------------------
       ! Preconditioning...
       !-----------------------------------------------------------------
       CALL PreconditioningIteration(m, A, q, PA, PS, t, dim, &
           da, dpa, dps, SystemScaling, MMatrix, LMatrix)
       !------------------------------------------------------------------
       CALL ComplexMatrixVectorProduct( A, t, work(1:n,r+k) )
       rnrm = dznrm2(n, work(1:n,r), 1)
       mxnrmx = MAX (mxnrmx, rnrm)
       mxnrmr = MAX (mxnrmr, rnrm)
     ENDDO

     !--------------------------------------
     ! --- The convex polynomial part ---
     !--------------------------------------

     DO i=1,l+1
       DO j=1,i
         rwork(i,j) = zdotc(n, work(1:n,r+i-1), 1, work(1:n,r+j-1),1 )
       END DO
     END DO
     DO j=2,l+1
       rwork(1:j-1,j) = CONJG( rwork(j,1:j-1) )
     END DO

     rwork(1:l+1,zz:zz+l) = rwork(1:l+1,z:z+l)
     CALL zgetrf (l-1, l-1, rwork(2:l,zz+1:zz+l-1), l-1, &
         iwork, stat)

     ! --- tilde r0 and tilde rl (small vectors)

     rwork(1,y0) = -zone
     rwork(2:l,y0) = rwork(2:l,z)
     CALL zgetrs('n', l-1, 1, rwork(2:l,zz+1:zz+l-1), l-1, iwork, &
         rwork(2:l,y0), l-1, stat)
     rwork(l+1,y0) = zzero

     rwork(1,yl) = zzero
     rwork(2:l,yl) = rwork(2:l,z+l)
     CALL zgetrs ('n', l-1, 1, rwork(2:l,zz+1:zz+l-1), l-1, iwork, &
         rwork(2:l,yl), l-1, stat)
     rwork(l+1,yl) = -zone

     ! --- Convex combination

     CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
         rwork(1:l+1,y0), 1, zzero, rwork(1:l+1,y), 1)
     kappa0 = SQRT( ABS(zdotc(l+1, rwork(1:l+1,y0), 1, rwork(1:l+1,y), 1)) )
     CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
         rwork(1:l+1,yl), 1, zzero, rwork(1:l+1,y), 1)
     kappal = SQRT( ABS(zdotc(l+1, rwork(1:l+1,yl), 1, rwork(1:l+1,y), 1)) )
     CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
         rwork(1:l+1,y0), 1, zzero, rwork(1:l+1,y), 1)
     varrho = zdotc(l+1, rwork(1:l+1,yl), 1, rwork(1:l+1,y), 1) / &
         (kappa0*kappal)
     hatgamma = varrho/ABS(varrho) * MAX(ABS(varrho),7d-1) * &
         kappa0/kappal
     rwork(1:l+1,y0) = rwork(1:l+1,y0) - hatgamma * rwork(1:l+1,yl)

     !  --- Update

     omega = rwork(l+1,y0)
     DO j=1,l
       work(1:n,u) = work(1:n,u) - rwork(j+1,y0) * work(1:n,u+j)
       x(1:n) = x(1:n) + rwork(j+1,y0) * work(1:n,r+j-1)
       work(1:n,r) = work(1:n,r) - rwork(j+1,y0) * work(1:n,r+j)
     ENDDO


     CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
         rwork(1:l+1,y0), 1, zzero, rwork(1:l+1,y), 1)
     rnrm = SQRT( ABS(zdotc(l+1, rwork(1:l+1,y0), 1, rwork(1:l+1,y), 1)) )

      !---------------------------------------
      !  --- The reliable update part ---
      !---------------------------------------
     IF (.FALSE.) THEN
       IF (.TRUE.) THEN
         mxnrmx = MAX (mxnrmx, rnrm)
         mxnrmr = MAX (mxnrmr, rnrm)
         xpdt = (rnrm < delta*rnrm0 .AND. rnrm0 < mxnrmx)
         rcmp = ((rnrm < delta*mxnrmr .AND. rnrm0 < mxnrmr) .OR. xpdt)
         IF (rcmp) THEN
           PRINT *, 'Performing residual update...'
           t(1:n) = x(1:n)
           CALL CRS_ComplexLUSolve2( n, A, t )
           CALL ComplexMatrixVectorProduct( A, t, work(1:n,r) )
           work(1:n,r) = work(1:n,bp) - work(1:n,r)
           mxnrmr = rnrm
           IF (xpdt) THEN
             PRINT *, 'Performing solution update...'
             !CALL CRS_ComplexLUSolve2( n, A, x )
             work(1:n,xp) = work(1:n,xp) + t(1:n)
             x(1:n) = zzero
             work(1:n,bp) = work(1:n,r)
             mxnrmx = rnrm
           ENDIF
         ENDIF
       ENDIF

       IF (rcmp) THEN
         IF (xpdt) THEN
           t(1:n) = work(1:n,xp)
         ELSE
           t(1:n) = t(1:n) + work(1:n,xp)
         END IF
       ELSE
         t(1:n) = x(1:n)
         CALL CRS_ComplexLUSolve2( n, A, t )
         t(1:n) =  t(1:n) + work(1:n,xp)
       END IF
     END IF

     !errorind = StoppingCriterion( n, A, t, b, work(1:n,r) )

     errorind = rnrm/bnrm
     WRITE(*,'(a,I4,ES12.3)') 'OuterIteration residual for iterate', Round, errorind

   END DO

   !------------------------------------------------------------
   ! We have solved z = P*x, so finally solve the true unknown x
   !------------------------------------------------------------
   !
   ! This is a nontrivial computation ...
   !
   !CALL CRS_ComplexLUSolve2( n, A, x )
   !x(1:n) = x(1:n) + work(1:n,xp)
   STOP

   WRITE(*,'(a,ES12.3)') 'An approximate lower bound for the condition number: ', &
       ConditionEstimate( m, A, x, b)
   !----------------------------------------------
   ! No transforming back: x aliases v.
   !----------------------------------------------
   Norm = SQRT(DOT_PRODUCT( v(1:2*n), v(1:2*n) )/(2*n))

!------------------------------------------------------------------------------
   END SUBROUTINE BiCGStablOuterIteration
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  SUBROUTINE InnerOuterIteration( n, A, q, PA, PS, MMatrix, LMatrix, &
      x, b, Rounds, TOL, dim, Norm )
!------------------------------------------------------------------------------
!   This is the nested GCR iteration for the complex linear system Ax=b.
!   The new search direction is solved from the residual equation As = r.
!   The residual equations is solved using the block-preconditioned GCR(m) method
!   ideally with a small m.
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: A, PA, PS, MMatrix, LMatrix
    INTEGER :: n, q, Rounds, dim
    REAL(KIND=dp), TARGET :: x(n), b(n)
    REAL(KIND=dp) :: TOL, Norm
!-------------------------------------------------------------------------------
    INTEGER :: i, j, k, m, InnerRounds, IluOrder, VelocityPrecond = 1, InnerRestart, &
        MaxRestarts
    LOGICAL :: Condition, GotIt, SystemScaling, Truncation, &
        ConvergedSol, FirstVisit = .TRUE.
    REAL(KIND=dp) :: res, tottime, res0, InnerTol, alpha, &
        ResidualReductionRatio, bw_error
    COMPLEX(KIND=dp) :: r(n/2),T1(n/2),T2(n/2), &
        S(n/2,Rounds), V(n/2,Rounds), &
        Vel(dim*q/2), VelRhs(dim*q/2), Sol(n/2), da(n/2), dps(q), dpa(dim*q/2)
    ! x and b already hold the complex vectors as consecutive (Re,Im) pairs, so
    ! alias them rather than copying into complex temporaries and back.
    COMPLEX(KIND=dp), POINTER :: y(:),f(:)
    COMPLEX(KIND=dp) :: beta
    SAVE FirstVisit
!------------------------------------------------------------------------------
    ConvergedSol = .FALSE.
    tottime = CPUTime()
    !-----------------------------------------------------------------------------
    ! Scaling to minimize the condition number
    !-----------------------------------------------------------------------------
    SystemScaling = ListGetLogical( Solver % Values, 'Linear System Scaling', GotIt)

    IF (.NOT. GotIt) SystemScaling = .TRUE.
    IF (SystemScaling) THEN
      WRITE( Message, * ) 'Scaling the system matrix...'
      CALL Info( 'AcousticsSolver', Message, Level=5 )
      CALL OptimalMatrixScaling( n/2, A, da )
      WRITE( Message, * ) 'Scaling the preconditioning matrix for Schur complement...'
      CALL Info( 'AcousticsSolver', Message, Level=5 )
      CALL OptimalMatrixScaling( q, PS, dps )
      IF ( ASSOCIATED(PA) ) THEN
        WRITE( Message, * ) 'Scaling the preconditioning matrix for velocities...'
        CALL Info( 'AcousticsSolver', Message, Level=5 )
        CALL OptimalMatrixScaling( dim*q/2, PA, dpa )
      END IF
    END IF

    !-------------------------------------------------------------------------------
    ! Compute ILU factorizations for the preconditioner matrices.
    ! This needs to be done only once.
    !-------------------------------------------------------------------------------
    IluOrder = ListGetInteger( Solver % Values, 'ILU Order for Schur Complement')

    CALL Info( 'AcousticsSolver', ' ', Level=4)
    CALL Info( 'AcousticsSolver', 'ILU factorization for the Schur complement preconditioner', &
        Level=4)
    CALL Info( 'AcousticsSolver', ' ', Level=4)
    Condition = CRS_ComplexIncompleteLU( PS, IluOrder )
    !Condition = CRS_ComplexILUT(PS,1.0d-2)

    IluOrder = ListGetInteger( Solver % Values, 'ILU Order for Velocities', GotIt)

    IF (GotIt .AND. ASSOCIATED(PA) ) THEN
      CALL Info( 'AcousticsSolver', ' ', Level=4)
      CALL Info( 'AcousticsSolver', 'ILU factorization for the Velocity preconditioner', &
          Level=4)
      CALL Info( 'AcousticsSolver', ' ', Level=4)
      Condition = CRS_ComplexIncompleteLU( PA, IluOrder )
    ELSE
      VelocityPrecond = 0
    END IF

    !--------------------------------------------------------------------------------
    !   Some initializations
    !--------------------------------------------------------------------------------
    InnerTol = ListGetConstReal( Solver % Values, &
        'Linear System Convergence Tolerance' )
    InnerRounds = ListGetInteger( Solver % Values, &
        'Linear System Max Iterations')
    MaxRestarts = ListGetInteger( Solver % Values, &
        'Max Outer GCR Cycles', GotIt)
    IF (.NOT. GotIt) MaxRestarts = 1
    InnerRestart = ListGetInteger( Solver % Values, 'Max Inner GCR Iterations', GotIt)
    IF (.NOT. GotIt) InnerRestart = 5
    ResidualReductionRatio = ListGetConstReal( Solver % Values, 'Residual Reduction Ratio', GotIt)
    IF (.NOT. GotIt) ResidualReductionRatio = 1.0d-1
    Truncation = ListGetLogical( Solver % Values, &
        'Use Truncation', GotIt)
    IF (.NOT. GotIt) Truncation = .FALSE.



    !--------------------------------------------------------------------------------------------
    ! The solution of an initial guess, the previous solution can be used as an initial guess
    !--------------------------------------------------------------------------------------------
    IF ( FirstVisit .OR. (.NOT. UtilizePreviousSolution) ) THEN
      x(1:n) = 0.0d0

      Vel(1:dim*q/2) = CMPLX( 0.0d0,0.0d0, kind=dp )
      DO i=1,dim
        DO j=1,q/2
          VelRhs((j-1)*dim+i) = CMPLX( b(2*(dim+2)*(j-1)+2*i-1), &
              b(2*(dim+2)*(j-1)+2*i), kind=dp )
        END DO
      END DO

      IF ( ANY( VelRhs /= CMPLX(0.0d0, 0.0d0, kind=dp) ) ) THEN
        IF ( ASSOCIATED(PA) ) THEN
          !-----------------------------------
          ! Scale the right-hand side vector
          !-----------------------------------
          IF (SystemScaling) THEN
            DO i=1,dim*q/2
              VelRhs(i) = dpa(i) * VelRhs(i)
            END DO
          END IF
          !-----------------------------------------------------------------
          CALL Info( 'AcousticsSolver', ' ', Level=4)
          CALL Info( 'AcousticsSolver', &
              'Solving initial guess for velocities', Level=4)
          CALL Info( 'AcousticsSolver', ' ', Level=4)
          !------------------------------------------------------------------
          !CALL ComplexBiCGStab( q*dim, PA, Vel, VelRhs, InnerRounds, 1.0d-7, 0 )
          CALL ComplexBiCGStabl( 4, q*dim/2, PA, Vel, VelRhs, InnerRounds, 1.0d-7, 0)
        ELSE
          !-----------------------------------
          ! Scale the right-hand side vector
          !-----------------------------------
          IF (SystemScaling) THEN
            DO i=1,dim
              DO j=1,q/2
                VelRhs((j-1)*dim+i) = da((j-1)*(dim+2)+i) * VelRhs((j-1)*dim+i)
              END DO
            END DO
          END IF
          !-----------------------------------------------------------------
          CALL Info( 'AcousticsSolver', ' ', Level=4)
          CALL Info( 'AcousticsSolver', &
              'Solving initial guess for velocities', Level=4)
          CALL Info( 'AcousticsSolver', ' ', Level=4)
          !------------------------------------------------------------------
          CALL VelocitySolve( 4, q*dim/2, A, Vel, VelRhs, dim, InnerRounds, 1.0d-7, 0)
        END IF

        DO j=1,q/2
          DO i=1,dim
            x(2*(dim+2)*(j-1)+2*i-1) = REAL(Vel((j-1)*dim+i))
            x(2*(dim+2)*(j-1)+2*i) = AIMAG(Vel((j-1)*dim+i))
          END DO
        END DO
      END IF
      FirstVisit = .FALSE.
    END IF

    !--------------------------------------------------------------------------------
    ! The start of the GCR iteration...
    !--------------------------------------------------------------------------------
    m = n/2
    !----------------------------------------------------------------------------
    ! View the solution vector x and the right-hand side vector b as the
    ! complex-valued vectors y and f. Nothing is copied, so the solution needs
    ! no transforming back either.
    !---------------------------------------------------------------------------
    y => ComplexValues( x, m )
    f => ComplexValues( b, m )

    CALL ComplexMatrixVectorProduct( A, y, r )
    !--------------------------------------------------------
    ! Scale the right-hand side vector and form the residual
    !--------------------------------------------------------
    IF (SystemScaling) THEN
      DO i=1,m
        f(i) = da(i) * f(i)
      END DO
    END IF
    r(1:m) = f(1:m) - r(1:m)
    res0 = ComplexNorm(m,f)
    res = ComplexNorm(m,r)/res0

    WRITE(*,'(a,I4,ES12.3,ES12.3,ES12.3)') 'OuterIteration residual for iterate', &
        0, res, StoppingCriterion(m,A,y,f,r), CPUTime() - tottime

    DO j=1,MaxRestarts
      IF (ConvergedSol) EXIT
      V(1:m,1:Rounds) = CMPLX( 0.0d0, 0.0d0, kind=dp)
      S(1:m,1:Rounds) = CMPLX( 0.0d0, 0.0d0, kind=dp)
      DO k=1,Rounds
        !----------------------------------------------------------
        ! Perform the preconditioning...
        !---------------------------------------------------------------
        T1(1:m) = r(1:m)
        Sol(1:m) = CMPLX(0.0d0, 0.0d0, kind=dp)
        CALL InnerIteration( n, A, q, PA, PS, Sol, T1, &
            InnerRestart, ResidualReductionRatio, dim, MMatrix, LMatrix, &
            SystemScaling, da, dpa, dps)
        T1(1:m) = Sol(1:m)
        CALL ComplexMatrixVectorProduct( A, T1, T2 )

        !--------------------------------------------------------------
        ! Perform the orthogonalisation of the search directions....
        !--------------------------------------------------------------
        DO i=1,k-1
          beta = ComplexDotProduct( m, V(1:m,i), T2(1:m) )
          T1(1:m) = T1(1:m) - beta * S(1:m,i)
          T2(1:m) = T2(1:m) - beta * V(1:m,i)
        END DO
        alpha = ComplexNorm(m,T2)
        T1(1:m) = CMPLX( 1.0d0, 0.0d0, kind=dp)/CMPLX( alpha, 0.0d0, kind=dp) * T1(1:m)
        T2(1:m) = CMPLX( 1.0d0, 0.0d0, kind=dp)/CMPLX( alpha, 0.0d0, kind=dp) * T2(1:m)

        !-------------------------------------------------------------
        ! The update of the solution and save the search data...
        !-------------------------------------------------------------
        beta = ComplexDotProduct(m, T2, r)
        y(1:m) = y(1:m) + beta * T1(1:m)
        r(1:m) = r(1:m) - beta * T2(1:m)
        S(1:m,k) = T1(1:m)
        V(1:m,k) = T2(1:m)

        !--------------------------------------------------------------
        ! Check whether the convergence criterion is met
        !--------------------------------------------------------------
        res = ComplexNorm(m,r)/res0

        bw_error = StoppingCriterion( m, A, y, f, r )
        !PRINT *,'OuterIteration residual',i,res, StoppingCriterion(m,A,y,f,r),CPUTime() - tottime

        WRITE(*,'(a,I4,ES12.3,ES12.3,ES12.3)') 'OuterIteration residual for iterate', &
            k + (j-1) * Rounds, res, bw_error, CPUTime() - tottime
        ConvergedSol = ( bw_error < TOL)
        IF (ConvergedSol) EXIT
      END DO
    END DO

    WRITE(*,'(a,ES12.3)') 'An approximate lower bound for the condition number: ', &
        ConditionEstimate( m, A, y, f )

    !----------------------------------------------
    ! No transforming back: y aliases x.
    !----------------------------------------------
    Norm = SQRT(DOT_PRODUCT( x(1:2*m), x(1:2*m) )/(2*m))

!------------------------------------------------------------------------------
   END SUBROUTINE InnerOuterIteration
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  SUBROUTINE InnerIteration( n, A, q, PA, PS, x, b, Rounds, TOL, dim, &
      MMatrix, LMatrix, SystemScaling, da, dpa, dps)
!------------------------------------------------------------------------------
!   This is the preconditioned GCR(Rounds) iteration for the complex linear system Ax=b.
!   The Schur complement preconditioning strategy is used.
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: A, PA, PS, MMatrix, LMatrix
    INTEGER :: n, q, Rounds, dim
    REAL(KIND=dp) :: TOL
    COMPLEX(KIND=dp) :: x(n/2), b(n/2)
    LOGICAL, OPTIONAL :: SystemScaling
    COMPLEX(KIND=dp), OPTIONAL :: da(n/2), dpa(dim*q/2), dps(q)
!-------------------------------------------------------------------------------
    INTEGER :: k, m
    INTEGER :: i
    REAL(KIND=dp) :: res, tottime, res0, alpha
    COMPLEX(KIND=dp) :: r(n/2), T1(n/2), T2(n/2), S(n/2,Rounds), V(n/2,Rounds)
    COMPLEX(KIND=dp) :: beta
!------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
!   The start of the GCR iteration...
!--------------------------------------------------------------------------------
    tottime = CPUTime()
    m = n/2

    CALL ComplexMatrixVectorProduct( A, x, r )
    r(1:m) = b(1:m) - r(1:m)
    res0 = ComplexNorm(m,r)

    V(1:m,1:Rounds) = CMPLX( 0.0d0, 0.0d0, kind=dp)
    S(1:m,1:Rounds) = CMPLX( 0.0d0, 0.0d0, kind=dp)

    DO k=1,Rounds
      !----------------------------------------------------------
      ! Perform the preconditioning...
      !---------------------------------------------------------------
      T1(1:m) = r(1:m)
      CALL PreconditioningIteration( n, A, q, PA, PS, T1, dim, &
          da, dpa, dps, SystemScaling, MMatrix, LMatrix)
      CALL ComplexMatrixVectorProduct( A, T1, T2 )

      !--------------------------------------------------------------
      ! Perform the orthogonalisation of the search directions....
      !--------------------------------------------------------------
      DO i=1,k-1
        beta = ComplexDotProduct( m, V(1:m,i), T2(1:m) )
        T1(1:m) = T1(1:m) - beta * S(1:m,i)
        T2(1:m) = T2(1:m) - beta * V(1:m,i)
      END DO
      alpha = ComplexNorm(m,T2)
      T1(1:m) = CMPLX( 1.0d0, 0.0d0, kind=dp)/CMPLX( alpha, 0.0d0, kind=dp) * T1(1:m)
      T2(1:m) = CMPLX( 1.0d0, 0.0d0, kind=dp)/CMPLX( alpha, 0.0d0, kind=dp) * T2(1:m)

      !-------------------------------------------------------------
      ! The update of the solution and save the search data...
      !-------------------------------------------------------------
      beta = ComplexDotProduct(m, T2, r)
      x(1:m) = x(1:m) + beta * T1(1:m)
      r(1:m) = r(1:m) - beta * T2(1:m)
      S(1:m,k) = T1(1:m)
      V(1:m,k) = T2(1:m)

      !--------------------------------------------------------------
      ! Check whether the convergence criterion is met
      !--------------------------------------------------------------
      res = ComplexNorm(m,r)/res0
      WRITE(*,'(a,I4,ES12.3,ES12.3)') 'InnerIteration residual for iterate', &
          k, res, CPUTime() - tottime
      IF ( res < TOL) EXIT
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE InnerIteration
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
     SUBROUTINE PreconditioningIteration( n, A, q, PA, PS, V, dim, &
         da, dpa, dps, Scaling, MMatrix, LMatrix)
!------------------------------------------------------------------------------
!     This subroutine solves iteratively the upper triangular preconditioning
!     system P*z = V. The vector V is overwritten by the solution z.
!------------------------------------------------------------------------------
      TYPE(Matrix_t), POINTER :: A, PA, PS
      INTEGER :: n, q, dim
      COMPLEX(KIND=dp) :: V(n/2)
      COMPLEX(KIND=dp), OPTIONAL :: da(n/2), dpa(dim*q/2), dps(q)
      LOGICAL, OPTIONAL :: Scaling
      TYPE(Matrix_t), POINTER, OPTIONAL :: MMatrix, LMatrix
!--------------------------------------------------------------------------------
      REAL(kind=dp) :: InnerTol, VelocityTol, SchurTol
      INTEGER :: i, j, InnerRounds
      LOGICAL :: SystemScaling, ConsistentSplitting, VelocityCriterion, &
          SchurCriterion, GotIt
      COMPLEX(kind=dp) :: z(n/2), y(q/2), f(q/2), w(q), g(q), &
          Vel(dim*q/2), VelRhs(dim*q/2)
!------------------------------------------------------------------------------
      InnerTol = ListGetConstReal( Solver % Values, &
          'Linear System Convergence Tolerance', GotIt )
      InnerRounds = ListGetInteger( Solver % Values, &
          'Linear System Max Iterations', GotIt)
      VelocityTol =  ListGetConstReal( Solver % Values, &
          'Velocity Convergence Tolerance', VelocityCriterion)
      SchurTol = ListGetConstReal( Solver % Values, &
          'Schur Complement Convergence Tolerance', SchurCriterion)

      ConsistentSplitting = PRESENT(MMatrix) .AND. PRESENT(LMatrix)

      IF ( PRESENT(Scaling) ) THEN
        SystemScaling = Scaling
      ELSE
        SystemScaling = .FALSE.
      END IF

      z(1:n/2) = CMPLX( 0.0d0,0.0d0, kind=dp )
      !-------------------------------------------------------------------------
      ! If the consistent splitting approach is used,
      ! compute the continuous approximation of the continuity equation residual
      !-------------------------------------------------------------------------
      IF ( ConsistentSplitting ) THEN
        y(1:q/2) = CMPLX( 0.0d0,0.0d0, kind=dp )
        IF (SystemScaling) THEN
          !-------------------------------
          ! Recover the unscaled residual
          !--------------------------------
          DO j=1,q/2
            f(j) = 1.0d0/da( (dim+2)*j ) * V( (dim+2)*j )
          END DO
        ELSE
          DO j=1,q/2
            f(j) = V( (dim+2)*j )
          END DO
        END IF
        !-----------------------------------------------------------------
        CALL Info( 'AcousticsSolver', ' ', Level=4)
        CALL Info( 'AcousticsSolver', &
            'Iteration for the continuity equation residual', Level=4)
        CALL Info( 'AcousticsSolver', ' ', Level=4)
        !------------------------------------------------------------------
        CALL ComplexBiCGStab( q, MMatrix, y, f, InnerRounds, 1.0d-7, 0 )
        CALL ComplexMatrixVectorProduct( LMatrix, y, f )
      END IF

      !-----------------------------------------------------------------
      ! The solution of the Schur complement equation...
      !------------------------------------------------------------------
      w(1:q) = CMPLX( 0.0d0,0.0d0, kind=dp )
      !------------------------------------------------------------------
      ! Construct the right-hand side vector g for the Schur complement
      ! system. The vector g consists of the unscaled residuals and the
      ! modification term which arises from the consistent splitting.
      !------------------------------------------------------------------
      IF ( ConsistentSplitting ) THEN

        IF ( SystemScaling ) THEN
          DO j=1,q/2
            g(2*j-1) = 1.0d0 / da( (dim+2)*j-1 ) * V((dim+2)*j-1)
            g(2*j) = 1.0d0 / da( (dim+2)*j ) * V((dim+2)*j) + f(j)
          END DO
        ELSE
          DO j=1,q/2
            g(2*j-1) = V((dim+2)*j-1)
            g(2*j) = V((dim+2)*j) + f(j)
          END DO
        END IF

      ELSE

        IF ( SystemScaling ) THEN
          g(2*j-1) = 1.0d0 / da( (dim+2)*j-1 ) * V((dim+2)*j-1)
          g(2*j) = 1.0d0 / da( (dim+2)*j ) * V((dim+2)*j)
        ELSE
          DO j=1,q/2
            g(2*j-1) = V((dim+2)*j-1)
            g(2*j) = V((dim+2)*j)
          END DO
        END IF

      END IF

      !--------------------------------------------------------------------
      ! If scaling is used, we must scale g according to the scaling of
      ! the Schur complement matrix
      !--------------------------------------------------------------------
      IF ( SystemScaling ) THEN
        DO j=1,q/2
          g(2*j-1) = dps( 2*j-1 ) * g(2*j-1)
          g(2*j) = dps( 2*j ) * g(2*j)
        END DO
      END IF
      !-------------------------------------------------------------------
      ! Solve the Schur complement system...
      !-------------------------------------------------------------------
      IF ( ANY( g /= CMPLX( 0.0d0,0.0d0, KIND=dp ) ) ) THEN
        !-----------------------------------------------------------------
        CALL Info( 'AcousticsSolver', ' ', Level=4)
        CALL Info( 'AcousticsSolver', &
            'Preconditioning iteration for the Schur complement system', Level=4)
        CALL Info( 'AcousticsSolver', ' ', Level=4)
        !-------------------------------------------------------------------
        IF (SchurCriterion) THEN
          CALL ComplexBiCGStabl( 4, q, PS, w, g, InnerRounds, SchurTol, 0 )
        ELSE
          CALL ComplexBiCGStabl( 4, q, PS, w, g, InnerRounds, InnerTol, 0 )
          !CALL ComplexBiCGStabl( 2, q, PS, w, g, 5, 1.0d-6, 0 )
        END IF
      END IF

      DO j=1,q/2
        V((dim+2)*j-1) = w(2*j-1)
        V((dim+2)*j) = w(2*j)
      END DO

      !---------------------------------------------------------------------------------
      ! The computation of the specific matrix vector product related to preconditining
      !---------------------------------------------------------------------------------
      CALL ComplexMatrixVectorProduct2( A, V, z, dim )

      !--------------------------------------------------------
      ! The solution of the velocity preconditioning equation
      !--------------------------------------------------------
      Vel(1:dim*q/2) = CMPLX( 0.0d0,0.0d0, kind=dp )
      DO i=1,dim
        DO j=1,q/2
          VelRhs((j-1)*dim+i) = z((dim+2)*(j-1)+i)
        END DO
      END DO


      IF ( ASSOCIATED(PA) ) THEN
        IF ( SystemScaling ) THEN
          DO i=1,dim
            DO j=1,q/2
              VelRhs((j-1)*dim+i) = dpa( (j-1)*dim+i ) / da( (dim+2)*(j-1)+i ) * &
                  VelRhs((j-1)*dim+i)
            END DO
          END DO
        END IF

        IF ( ANY( VelRhs /= CMPLX( 0.0d0,0.0d0, KIND=dp ) ) ) THEN
          !-----------------------------------------------------------------
          CALL Info( 'AcousticsSolver', ' ', Level=4)
          CALL Info( 'AcousticsSolver', &
              'Preconditioning iteration for velocities', Level=4)
          CALL Info( 'AcousticsSolver', ' ', Level=4)
          !------------------------------------------------------------------
          IF (VelocityCriterion) THEN
            CALL ComplexBiCGStabl( 2, q*dim/2, PA, Vel, VelRhs, InnerRounds, VelocityTol, 0)
          ELSE
            CALL ComplexBiCGStabl( 2, q*dim/2, PA, Vel, VelRhs, InnerRounds, InnerTol, 0)
          END IF
        END IF

      ELSE

        IF ( ANY( VelRhs /= CMPLX( 0.0d0,0.0d0, KIND=dp ) ) ) THEN
          !-----------------------------------------------------------------
          CALL Info( 'AcousticsSolver', ' ', Level=4)
          CALL Info( 'AcousticsSolver', &
              'Preconditioning iteration for velocities', Level=4)
          CALL Info( 'AcousticsSolver', ' ', Level=4)
          !------------------------------------------------------------------
          IF (VelocityCriterion) THEN
            CALL VelocitySolve( 2, q*dim/2, A, Vel, VelRhs, dim, InnerRounds, VelocityTol, 0)
          ELSE
            CALL VelocitySolve( 2, q*dim/2, A, Vel, VelRhs, dim, InnerRounds, InnerTol, 0)
          END IF
        END IF
      END IF

      DO j=1,q/2
        DO i=1,dim
          z((dim+2)*(j-1)+i) = Vel((j-1)*dim+i)
        END DO
      END DO

      V(1:n/2)=z(1:n/2)
!------------------------------------------------------------------------------
    END SUBROUTINE PreconditioningIteration
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  FUNCTION StoppingCriterion( n, A, x, b, res, CriterionType ) RESULT(err)
!------------------------------------------------------------------------------
    INTEGER :: n
    TYPE(Matrix_t), POINTER :: A
    COMPLEX(KIND=dp) :: x(n), b(n), res(n)
    INTEGER, OPTIONAL :: CriterionType
    REAL(kind=dp) :: err
!-----------------------------------------------------------------------------
    REAL(kind=dp) :: norm, tmp, normb, normres, normx
    INTEGER :: i, j
    INTEGER, POINTER :: Cols(:), Rows(:)
    REAL(KIND=dp), POINTER :: Values(:)
!--------------------------------------------------------------------------
    Rows   => A % Rows
    Cols   => A % Cols
    Values => A % Values

    norm = 0.0d0
    normb = 0.0d0
    normres = 0.0d0
    normx = 0.0d0

    IF (PRESENT(CriterionType)) THEN
      DO i=1,n
        normb = MAX( normb, CDABS(b(i)) )
        normres = MAX( normres, CDABS(res(i)) )
      END DO
      err = normres/normb
    ELSE
      DO i=1,n
        tmp = 0.0d0
        DO j=Rows(2*i-1),Rows(2*i)-1,2
          tmp = tmp + CDABS( CMPLX( Values(j), -Values(j+1), kind=dp ) )
        END DO
        IF (tmp > norm) norm = tmp
        normb = MAX( normb, CDABS(b(i)) )
        normres = MAX( normres, CDABS(res(i)) )
        normx = MAX( normx, CDABS(x(i)) )
      END DO
      err = normres / (norm * normx + normb)
    END IF

!------------------------------------------------------------------------------
  END FUNCTION StoppingCriterion
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  FUNCTION ConditionEstimate( n, A, x, b ) RESULT(err)
!------------------------------------------------------------------------------
    INTEGER :: n
    TYPE(Matrix_t), POINTER :: A
    COMPLEX(KIND=dp) :: x(n), b(n)
    REAL(kind=dp) :: err
!------------------------------------------------------------------------------
    INTEGER :: i, j
    REAL(kind=dp) :: norm, tmp, normb, normx
    INTEGER, POINTER :: Cols(:), Rows(:)
    REAL(KIND=dp), POINTER :: Values(:)
!-------------------------------------------------------------------------------
    Rows   => A % Rows
    Cols   => A % Cols
    Values => A % Values

    norm = 0.0d0
    normb = 0.0d0
    normx = 0.0d0
    DO i=1,n
      tmp = 0.0d0
      DO j=Rows(2*i-1),Rows(2*i)-1,2
        tmp = tmp + CDABS( CMPLX( Values(j), -Values(j+1), kind=dp ) )
      END DO
      IF (tmp > norm) norm = tmp
      normb = MAX( normb, CDABS(b(i)) )
      normx = MAX( normx, CDABS(x(i)) )
    END DO

    err = norm*normx/normb
!------------------------------------------------------------------------------
  END FUNCTION ConditionEstimate
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE ComplexMatrixVectorProduct( A,u,v )
!------------------------------------------------------------------------------
!   Matrix vector product (v = Au) for a matrix given in CRS format.
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp), DIMENSION(*) :: u,v
    TYPE(Matrix_t), POINTER :: A
!------------------------------------------------------------------------------
    INTEGER, POINTER :: Cols(:),Rows(:)
    REAL(KIND=dp), POINTER :: Values(:)

    INTEGER :: i,j,n
    COMPLEX(KIND=dp) :: s
!------------------------------------------------------------------------------

    n = A % NumberOfRows / 2
    Rows   => A % Rows
    Cols   => A % Cols
    Values => A % Values

    v(1:n) = CMPLX( 0.0d0, 0.0d0, kind=dp )
    DO i=1,n
       DO j=Rows(2*i-1),Rows(2*i)-1,2
          s = CMPLX( Values(j), -Values(j+1), kind=dp )
          v(i) = v(i) + s * u((Cols(j)+1)/2)
       END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE ComplexMatrixVectorProduct
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  FUNCTION ComplexInfinityNorm( m, x ) RESULT(s)
!------------------------------------------------------------------------------
    INTEGER :: m
    COMPLEX(KIND=dp) :: x(m)
    REAL(KIND=dp) :: s
    INTEGER :: i
!--------------------------------------------------------------------
    s = 0.0d0
    DO i=1,m
      s = MAX( s, ABS(x(i)) )
    END DO
!--------------------------------------
 END FUNCTION ComplexInfinityNorm
!--------------------------------------




!------------------------------------------------------------------------------
    FUNCTION ComplexNorm( n, x ) RESULT(s)
!------------------------------------------------------------------------------
       INTEGER :: n, i
       REAL(KIND=dp) :: s
       COMPLEX(KIND=dp) :: r, x(:)
!------------------------------------------------------------------------------
       s =  SQRT( REAL( DOT_PRODUCT( x(1:n), x(1:n) ) ) )
!------------------------------------------------------------------------------
    END FUNCTION ComplexNorm
!------------------------------------------------------------------------------





!------------------------------------------------------------------------------
    FUNCTION ComplexDotProduct( n, x, y ) RESULT(s)
!------------------------------------------------------------------------------
       INTEGER :: n
       COMPLEX(KIND=dp) :: s, x(:), y(:)
!------------------------------------------------------------------------------
       s = DOT_PRODUCT( x(1:n), y(1:n) )
!------------------------------------------------------------------------------
    END FUNCTION ComplexDotProduct
!------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
  SUBROUTINE ComplexBiCGStabl( l, n, A, x, b, MaxRounds, Tol, StoppingCriterionType )
!-----------------------------------------------------------------------------------
!   This subroutine solves complex linear systems by using the BiCGStab(l) algorithm
!   with l >= 2. It has been developed by using as a starting point the work of D.R. Fokkema
!   (subroutine zbistbl v1.1 1998). Dr. Fokkema has given the right to distribute
!   the derived work under GPL and hence the original more conservative
!   copyright notice of the subroutine has been removed accordingly.
!
!   This version uses right-oriented ILU(n) preconditioning.
!-----------------------------------------------------------------------------------
    INTEGER :: l   ! polynomial degree
    INTEGER :: n, MaxRounds
    TYPE(Matrix_t), POINTER :: A
    COMPLEX(KIND=dp) :: x(n), b(n)
    REAL(KIND=dp) :: Tol
    INTEGER, OPTIONAL :: StoppingCriterionType
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp) :: zzero, zone, t(n), kappa0, kappal
    REAL(KIND=dp) :: rnrm0, rnrm, mxnrmx, mxnrmr, errorind, &
        delta = 1.0d-2, bnrm, bw_errorind, tottime
    INTEGER :: i, j, rr, r, u, xp, bp, z, zz, y0, yl, y, k, iwork(l-1), stat, Round, &
        IluOrder
    COMPLEX(KIND=dp) :: work(n,3+2*(l+1)), rwork(l+1,3+2*(l+1)), &
        alpha, beta, omega, rho0, rho1, sigma, varrho, hatgamma
    LOGICAL rcmp, xpdt, GotIt, BackwardError
    CHARACTER(LEN=MAX_NAME_LEN) :: str
!------------------------------------------------------------------------------
    tottime = CPUTime()
    BackwardError = .TRUE.
    IF ( PRESENT(StoppingCriterionType) ) THEN
      IF (StoppingCriterionType == 0) BackwardError = .FALSE.
    END IF
    IF ( ALL(x == CMPLX(0.0d0,0.0d0,kind=dp)) ) x = b

    zzero = dcmplx( 0.0d0,0.0d0)
    zone =  dcmplx( 1.0d0,0.0d0)
    work = dcmplx( 0.0d0, 0.0d0 )
    rwork = dcmplx( 0.0d0, 0.0d0 )

    rr = 1
    r = rr+1
    u = r+(l+1)
    xp = u+(l+1)
    bp = xp+1

    z = 1
    zz = z+(l+1)
    y0 = zz+(l+1)
    yl = y0+1
    y = yl+1

    CALL ComplexMatrixVectorProduct( A, x, work(1:n,r) )
    work(1:n,r) = b(1:n) - work(1:n,r)
    bnrm = dznrm2(n, b(1:n), 1)

    work(1:n,rr) = work(1:n,r)
    work(1:n,bp) = work(1:n,r)
    work(1:n,xp) = x(1:n)
    x(1:n) = zzero
    rnrm0 = dznrm2(n, work(1:n,r), 1)
    rnrm = rnrm0
    mxnrmx = rnrm0
    mxnrmr = rnrm0

    alpha = zzero
    omega = zone
    sigma = zone
    rho0 = zone

    Round = 0
    errorind = 1.0d0
    DO WHILE ( errorind > Tol .AND. Round < MaxRounds)
      Round = Round + 1
      !-------------------------
      ! --- The BiCG part ---
      !-------------------------
      rho0 = -omega*rho0

      DO k=1,l
        rho1 = zdotc(n, work(1:n,rr), 1, work(1:n,r+k-1), 1)
        IF (rho0 == zzero) THEN
          CALL Fatal( 'ComplexBiCGStab(l)', 'Breakdown error.' )
        ENDIF
        beta = alpha*(rho1/rho0)
        rho0 = rho1
        DO j=0,k-1
          work(1:n,u+j) = work(1:n,r+j) - beta*work(1:n,u+j)
        ENDDO
        t(1:n) = work(1:n,u+k-1)
        CALL CRS_ComplexLUSolve2( n, A, t )
        CALL ComplexMatrixVectorProduct( A, t, work(1:n,u+k) )
        sigma = zdotc(n, work(1:n,rr), 1, work(1:n,u+k), 1)
        IF (sigma == zzero) THEN
          CALL Fatal( 'ComplexBiCGStab(l)', 'Breakdown error.' )
        ENDIF
        alpha = rho1/sigma
        x(1:n) = x(1:n) + alpha * work(1:n,u)
        DO j=0,k-1
          work(1:n,r+j) = work(1:n,r+j) - alpha * work(1:n,u+j+1)
        ENDDO
        t(1:n) = work(1:n,r+k-1)
        CALL CRS_ComplexLUSolve2( n, A, t )
        CALL ComplexMatrixVectorProduct( A, t, work(1:n,r+k) )
        rnrm = dznrm2(n, work(1:n,r), 1)
        mxnrmx = MAX (mxnrmx, rnrm)
        mxnrmr = MAX (mxnrmr, rnrm)
      ENDDO

      !--------------------------------------
      ! --- The convex polynomial part ---
      !--------------------------------------

      DO i=1,l+1
        DO j=1,i
          rwork(i,j) = zdotc(n, work(1:n,r+i-1), 1, work(1:n,r+j-1),1 )
        END DO
      END DO
      DO j=2,l+1
        rwork(1:j-1,j) = CONJG( rwork(j,1:j-1) )
      END DO

      rwork(1:l+1,zz:zz+l) = rwork(1:l+1,z:z+l)
      CALL zgetrf (l-1, l-1, rwork(2:l,zz+1:zz+l-1), l-1, &
          iwork, stat)

      ! --- tilde r0 and tilde rl (small vectors)

      rwork(1,y0) = -zone
      rwork(2:l,y0) = rwork(2:l,z)
      CALL zgetrs('n', l-1, 1, rwork(2:l,zz+1:zz+l-1), l-1, iwork, &
          rwork(2:l,y0), l-1, stat)
      rwork(l+1,y0) = zzero

      rwork(1,yl) = zzero
      rwork(2:l,yl) = rwork(2:l,z+l)
      CALL zgetrs ('n', l-1, 1, rwork(2:l,zz+1:zz+l-1), l-1, iwork, &
          rwork(2:l,yl), l-1, stat)
      rwork(l+1,yl) = -zone

      ! --- Convex combination

      CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
          rwork(1:l+1,y0), 1, zzero, rwork(1:l+1,y), 1)
      kappa0 = SQRT( ABS(zdotc(l+1, rwork(1:l+1,y0), 1, rwork(1:l+1,y), 1)) )
      CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
          rwork(1:l+1,yl), 1, zzero, rwork(1:l+1,y), 1)
      kappal = SQRT( ABS(zdotc(l+1, rwork(1:l+1,yl), 1, rwork(1:l+1,y), 1)) )
      CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
          rwork(1:l+1,y0), 1, zzero, rwork(1:l+1,y), 1)
      varrho = zdotc(l+1, rwork(1:l+1,yl), 1, rwork(1:l+1,y), 1) / &
          (kappa0*kappal)
      hatgamma = varrho/ABS(varrho) * MAX(ABS(varrho),7d-1) * &
          kappa0/kappal
      rwork(1:l+1,y0) = rwork(1:l+1,y0) - hatgamma * rwork(1:l+1,yl)

      !  --- Update

      omega = rwork(l+1,y0)
      DO j=1,l
        work(1:n,u) = work(1:n,u) - rwork(j+1,y0) * work(1:n,u+j)
        x(1:n) = x(1:n) + rwork(j+1,y0) * work(1:n,r+j-1)
        work(1:n,r) = work(1:n,r) - rwork(j+1,y0) * work(1:n,r+j)
      ENDDO


      CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
          rwork(1:l+1,y0), 1, zzero, rwork(1:l+1,y), 1)
      rnrm = SQRT( ABS(zdotc(l+1, rwork(1:l+1,y0), 1, rwork(1:l+1,y), 1)) )

      !---------------------------------------
      !  --- The reliable update part ---
      !---------------------------------------

      IF (.TRUE.) THEN
        mxnrmx = MAX (mxnrmx, rnrm)
        mxnrmr = MAX (mxnrmr, rnrm)
        xpdt = (rnrm < delta*rnrm0 .AND. rnrm0 < mxnrmx)
        rcmp = ((rnrm < delta*mxnrmr .AND. rnrm0 < mxnrmr) .OR. xpdt)
        IF (rcmp) THEN
          ! PRINT *, 'Performing residual update...'
          t(1:n) = x(1:n)
          CALL CRS_ComplexLUSolve2( n, A, t )
          CALL ComplexMatrixVectorProduct( A, t, work(1:n,r) )
          work(1:n,r) = work(1:n,bp) - work(1:n,r)
          mxnrmr = rnrm
          IF (xpdt) THEN
            ! PRINT *, 'Performing solution update...'
            work(1:n,xp) = work(1:n,xp) + t(1:n)
            x(1:n) = zzero
            work(1:n,bp) = work(1:n,r)
            mxnrmx = rnrm
          ENDIF
        ENDIF
      ENDIF

      IF (rcmp) THEN
        IF (xpdt) THEN
          t(1:n) = work(1:n,xp)
        ELSE
          t(1:n) = t(1:n) + work(1:n,xp)
        END IF
      ELSE
        t(1:n) = x(1:n)
        CALL CRS_ComplexLUSolve2( n, A, t )
        t(1:n) =  t(1:n) + work(1:n,xp)
      END IF

      bw_errorind = StoppingCriterion( n, A, t, b, work(1:n,r) )
      errorind = rnrm/bnrm
      WRITE(*,'(I4,ES12.3,ES12.3)') Round, errorind, bw_errorind

      IF (BackwardError) errorind = bw_errorind

    END DO

    !------------------------------------------------------------
    ! We have solved z = P*x, so finally solve the true unknown x
    !------------------------------------------------------------
    CALL CRS_ComplexLUSolve2( n, A, x )
    x(1:n) = x(1:n) + work(1:n,xp)

    WRITE(*,'(a,ES12.3)') 'An approximate lower bound for the condition number ', &
        ConditionEstimate( n, A, x, b )
    WRITE(*,'(a,ES12.3)') 'The 2-norm of the solution: ', &
        ComplexNorm( n, x )

!------------------------------------------------------------------------------
  END SUBROUTINE ComplexBiCGStabl
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
    SUBROUTINE ComplexBiCGStab( n, A, x, b, Rounds, TOL, StoppingCriterionType)
!------------------------------------------------------------------------------
!   This is the ILU or Jacobi preconditioned BiCGStab method for the complex
!   linear system Ax = b. The optional argument StoppingCriterionType can be
!   used to define the type of the stopping criterion.
!------------------------------------------------------------------------------
      TYPE(Matrix_t), POINTER :: A
      INTEGER :: n, Rounds
      REAL(KIND=dp) :: TOL
      COMPLEX(kind=dp) :: x(n/2), b(n/2)
      INTEGER, OPTIONAL :: StoppingCriterionType
!------------------------------------------------------------------------------
      INTEGER :: i, m, k
      LOGICAL :: BackwardError
      REAL(KIND=dp) :: res, tottime, res0, const
      COMPLEX(KIND=dp) :: r(n/2),Ri(n/2),P(n/2),V(n/2),T(n/2),T1(n/2),T2(n/2),&
          S(n/2)
      COMPLEX(KIND=dp) :: alpha,beta,omega,rho,oldrho
!------------------------------------------------------------------------------

      BackwardError = .TRUE.
      IF ( PRESENT(StoppingCriterionType) ) THEN
        IF (StoppingCriterionType == 0) BackwardError = .FALSE.
      END IF

      IF ( ALL(x == CMPLX(0.0d0,0.0d0,kind=dp)) ) x = b
      m = n/2

      CALL ComplexMatrixVectorProduct( A, x, r )
      r(1:m) = b(1:m) - r(1:m)
      IF ( .NOT. BackwardError) res0 = ComplexNorm(m,b)

      Ri(1:m) = r(1:m)
      P(1:m) = CMPLX( 0.0d0, 0.0d0, kind=dp)
      V(1:m) = CMPLX( 0.0d0, 0.0d0, kind=dp)
      omega  = CMPLX( 1.0d0, 0.0d0, kind=dp)
      alpha  = CMPLX( 0.0d0, 0.0d0, kind=dp)
      oldrho = CMPLX( 1.0d0, 0.0d0, kind=dp)
      tottime = CPUTime()

      DO i=1,Rounds
        rho = ComplexDotProduct( m, Ri, r )
        beta = alpha * rho / ( oldrho * omega )
        P(1:m) = r(1:m) + beta * (P(1:m) - omega*V(1:m))
        V(1:m) = P(1:m)

        CALL CRS_ComplexLUSolve2( m, A, V )

        T1(1:m) = V(1:m)
        CALL ComplexMatrixVectorProduct( A, T1, V )
        alpha = rho / ComplexDotProduct( m, Ri, V )
        S(1:m) = r(1:m) - alpha * V(1:m)
        x(1:m) = x(1:m) + alpha*T1(1:m)

        IF (BackwardError) THEN
          res = StoppingCriterion( m, A, x, b, S)
        ELSE
          res = ComplexNorm(m,S)/res0
        END IF

        IF ( res < TOL ) THEN
          WRITE(*,'(I4,ES12.3)') i, res
          WRITE(*,'(a,F8.2)') 'Solution time (s):    ', CPUTime() - tottime
          WRITE(*,'(a,ES12.3)') 'An approximate lower bound for the condition number ', &
              ConditionEstimate( m, A, x, b )
          WRITE(*,'(a,ES12.3)') 'The 2-norm of the solution: ', &
              ComplexNorm( m, x )
          EXIT
        END IF

        T(1:m) = S(1:m)

        CALL CRS_ComplexLUSolve2( m, A, T )

        T2(1:m) = T(1:m)
        CALL ComplexMatrixVectorProduct( A, T2, T )
        omega = ComplexDotProduct( m,T,S ) / ComplexDotProduct( m,T,T )
        oldrho = rho
        r(1:m) = S(1:m) - omega*T(1:m)
        x(1:m) = x(1:m) + omega*T2(1:m)

        IF (BackwardError) THEN
          res = StoppingCriterion( m, A, x, b, r )
        ELSE
          res = ComplexNorm(m,r)/res0
        END IF

        WRITE(*,'(I4,ES12.3)') i, res
        IF ( res < TOL ) THEN
          WRITE(*,'(a,F8.2)') 'Solution time (s):    ', CPUTime() - tottime
          WRITE(*,'(a,ES12.3)') 'An approximate lower bound for the condition number ', &
              ConditionEstimate( m, A, x, b )
          WRITE(*,'(a,ES12.3)') 'The 2-norm of the solution: ', &
              ComplexNorm( m, x )
          EXIT
        END IF
      END DO

!------------------------------------------------------------------------------
    END SUBROUTINE ComplexBiCGStab
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE CRS_ComplexLUSolve2( N,A,b )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  DESCRIPTION:
!    Solve a (complex9 system (Ax=b) after factorization A=LUD has been
!    done. This routine is meant as a part of  a preconditioner for an
!    iterative solver.
!
!  ARGUMENTS:
!
!    INTEGER :: N
!      INPUT: Size of the system
!
!    TYPE(Matrix_t) :: A
!      INPUT: Structure holding input matrix
!
!    DOUBLE PRECISION :: b
!      INOUT: on entry the RHS vector, on exit the solution vector.
!
!******************************************************************************
!------------------------------------------------------------------------------

    TYPE(Matrix_t), POINTER :: A
    INTEGER :: N
    COMPLEX(KIND=dp) :: b(N)

!------------------------------------------------------------------------------

    COMPLEX(KIND=dp), POINTER :: Values(:)
    INTEGER :: i,j
    COMPLEX(KIND=dp) :: x, s
    INTEGER, POINTER :: Cols(:),Rows(:),Diag(:)

!------------------------------------------------------------------------------

    Diag => A % ILUDiag
    Rows => A % ILURows
    Cols => A % ILUCols
    Values => A % CILUValues

!
!   if no ilu provided do diagonal solve:
!   -------------------------------------
    IF ( .NOT. ASSOCIATED( Values ) ) THEN
       Diag => A % Diag

       !DO i=1,n/2    ! THIS IS PROBABLY ERRATIC IN CVS-VERSION
       DO i=1,n
          x = CMPLX( A % Values(Diag(2*i-1)), -A % Values(Diag(2*i-1)+1), kind=dp )
          b(i) = b(i) / x
       END DO
       RETURN
    END IF

    IF (.FALSE.) THEN
      CALL ComplexLUSolve( n,SIZE(Cols),Rows,Cols,Diag,Values,b )
    ELSE
      ! Forward substitute
      DO i=1,n
        s = b(i)
        DO j=Rows(i),Diag(i)-1
          s = s - Values(j) * b(Cols(j))
        END DO
        b(i) = s
      END DO
      !
      ! Backward substitute
      DO i=n,1,-1
        s = b(i)
        DO j=Diag(i)+1,Rows(i+1)-1
          s = s - Values(j) * b(Cols(j))
        END DO
        b(i) = Values(Diag(i)) * s
      END DO
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE CRS_ComplexLUSolve2
!------------------------------------------------------------------------------

!-----------------------------------------------------------------------------
  SUBROUTINE ComplexLUSolve( n,m,Rows,Cols,Diag,Values,b )
!-----------------------------------------------------------------------------
    INTEGER :: n,m,Rows(n+1),Cols(m),Diag(n)
    COMPLEX(KIND=dp) :: Values(m),b(n)
    INTEGER :: i,j


    ! Forward substitute
    DO i=1,n
      DO j=Rows(i),Diag(i)-1
        b(i) = b(i) - Values(j) * b(Cols(j))
      END DO
    END DO

    ! Backward substitute
    DO i=n,1,-1
      DO j=Diag(i)+1,Rows(i+1)-1
        b(i) = b(i) - Values(j) * b(Cols(j))
      END DO
      b(i) = Values(Diag(i)) * b(i)
    END DO
!-----------------------------------------------------------------------------
  END SUBROUTINE ComplexLUSolve
!------------------------------------------------------------------------------




!-------------------------------------------------------------------------------
  SUBROUTINE ComplexMatrixVectorProduct2( A,u,v,dim )
!------------------------------------------------------------------------------
!
!   The computation of a specific matrix-vector product for preconditioning
!
!------------------------------------------------------------------------------

    COMPLEX(KIND=dp), DIMENSION(*) :: u,v
    TYPE(Matrix_t), POINTER :: A
    INTEGER :: dim
!------------------------------------------------------------------------------
    INTEGER, POINTER :: Cols(:), Rows(:), Diag(:)
    REAL(KIND=dp), POINTER :: Values(:)

    INTEGER :: i,j,k,n,p,q
    COMPLEX(KIND=dp) :: s
!------------------------------------------------------------------------------

    n = A % NumberOfRows / 2
    q = n/(dim+2)

    Rows   => A % Rows
    Cols   => A % Cols
    Diag   => A % Diag
    Values => A % Values

    v(1:n) = u(1:n)

    DO k=1,q
      DO p=1,dim
        i = (k-1)*(dim+2)+p
        DO j = Rows(2*i-1)+2*dim, Rows(2*i)-1, 2*(dim+2)
          s = CMPLX( Values(j), -Values(j+1), kind=dp )
          v(i) = v(i) - s * u((Cols(j)+1)/2)
        END DO
        DO j = Rows(2*i-1)+2*(dim+1), Rows(2*i)-1, 2*(dim+2)
          s = CMPLX( Values(j), -Values(j+1), kind=dp )
          v(i) = v(i) - s * u((Cols(j)+1)/2)
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------------
   END SUBROUTINE ComplexMatrixVectorProduct2
!------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
  SUBROUTINE VelocitySolve( l, n, A, x, b, dim, MaxRounds, Tol, StoppingCriterionType )
!-----------------------------------------------------------------------------------
!   This subroutine solves the velocity preconditioning system without requiring
!   that the assembly has been made for the preconditioning system. The required
!   matrix-vector products are performed by extracting the required entries from
!   the primary coefficient matrix A. This version uses diagonal preconditioning.
!-----------------------------------------------------------------------------------
    INTEGER :: l   ! polynomial degree
    INTEGER :: n, MaxRounds, dim
    TYPE(Matrix_t), POINTER :: A
    COMPLEX(KIND=dp) :: x(n), b(n)
    REAL(KIND=dp) :: Tol
    INTEGER, OPTIONAL :: StoppingCriterionType
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp) :: zzero, zone, t(n), kappa0, kappal
    REAL(KIND=dp) :: rnrm0, rnrm, mxnrmx, mxnrmr, errorind, &
        delta = 1.0d-2, bnrm, bw_errorind, tottime
    INTEGER :: i, j, rr, r, u, xp, bp, z, zz, y0, yl, y, k, iwork(l-1), stat, Round, &
        IluOrder
    COMPLEX(KIND=dp) :: work(n,3+2*(l+1)), rwork(l+1,3+2*(l+1)), &
        alpha, beta, omega, rho0, rho1, sigma, varrho, hatgamma
    LOGICAL rcmp, xpdt, GotIt, BackwardError
    CHARACTER(LEN=MAX_NAME_LEN) :: str
!------------------------------------------------------------------------------
    tottime = CPUTime()
    BackwardError = .FALSE.
    IF ( PRESENT(StoppingCriterionType) ) THEN
      IF (StoppingCriterionType == 0) BackwardError = .FALSE.
    END IF
    IF ( ALL(x == CMPLX(0.0d0,0.0d0,kind=dp)) ) x = b

    zzero = dcmplx( 0.0d0,0.0d0)
    zone =  dcmplx( 1.0d0,0.0d0)
    work = dcmplx( 0.0d0, 0.0d0 )
    rwork = dcmplx( 0.0d0, 0.0d0 )

    rr = 1
    r = rr+1
    u = r+(l+1)
    xp = u+(l+1)
    bp = xp+1

    z = 1
    zz = z+(l+1)
    y0 = zz+(l+1)
    yl = y0+1
    y = yl+1

    CALL ComplexMatrixVelocityVectorProduct( A, x, work(1:n,r), dim )
    work(1:n,r) = b(1:n) - work(1:n,r)
    bnrm = dznrm2(n, b(1:n), 1)

    work(1:n,rr) = work(1:n,r)
    work(1:n,bp) = work(1:n,r)
    work(1:n,xp) = x(1:n)
    x(1:n) = zzero
    rnrm0 = dznrm2(n, work(1:n,r), 1)
    rnrm = rnrm0
    mxnrmx = rnrm0
    mxnrmr = rnrm0

    alpha = zzero
    omega = zone
    sigma = zone
    rho0 = zone

    Round = 0
    errorind = 1.0d0
    DO WHILE ( errorind > Tol .AND. Round < MaxRounds)
      Round = Round + 1
      !-------------------------
      ! --- The BiCG part ---
      !-------------------------
      rho0 = -omega*rho0

      DO k=1,l
        rho1 = zdotc(n, work(1:n,rr), 1, work(1:n,r+k-1), 1)
        IF (rho0 == zzero) THEN
          CALL Fatal( 'ComplexBiCGStab(l)', 'Breakdown error.' )
        ENDIF
        beta = alpha*(rho1/rho0)
        rho0 = rho1
        DO j=0,k-1
          work(1:n,u+j) = work(1:n,r+j) - beta*work(1:n,u+j)
        ENDDO
        t(1:n) = work(1:n,u+k-1)
        !CALL CRS_ComplexLUSolve2( n, A, t )
        CALL DiagonalVelocityPreconditioning( A, t, dim)
        !CALL ComplexMatrixVectorProduct( A, t, work(1:n,u+k) )
        CALL ComplexMatrixVelocityVectorProduct( A, t, work(1:n,u+k), dim )
        sigma = zdotc(n, work(1:n,rr), 1, work(1:n,u+k), 1)
        IF (sigma == zzero) THEN
          CALL Fatal( 'ComplexBiCGStab(l)', 'Breakdown error.' )
        ENDIF
        alpha = rho1/sigma
        x(1:n) = x(1:n) + alpha * work(1:n,u)
        DO j=0,k-1
          work(1:n,r+j) = work(1:n,r+j) - alpha * work(1:n,u+j+1)
        ENDDO
        t(1:n) = work(1:n,r+k-1)
        !CALL CRS_ComplexLUSolve2( n, A, t )
        CALL DiagonalVelocityPreconditioning( A, t, dim)
        !CALL ComplexMatrixVectorProduct( A, t, work(1:n,r+k) )
        CALL ComplexMatrixVelocityVectorProduct( A, t, work(1:n,r+k), dim )
        rnrm = dznrm2(n, work(1:n,r), 1)
        mxnrmx = MAX (mxnrmx, rnrm)
        mxnrmr = MAX (mxnrmr, rnrm)
      ENDDO

      !--------------------------------------
      ! --- The convex polynomial part ---
      !--------------------------------------

      DO i=1,l+1
        DO j=1,i
          rwork(i,j) = zdotc(n, work(1:n,r+i-1), 1, work(1:n,r+j-1),1 )
        END DO
      END DO
      DO j=2,l+1
        rwork(1:j-1,j) = CONJG( rwork(j,1:j-1) )
      END DO

      rwork(1:l+1,zz:zz+l) = rwork(1:l+1,z:z+l)
      CALL zgetrf (l-1, l-1, rwork(2:l,zz+1:zz+l-1), l-1, &
          iwork, stat)

      ! --- tilde r0 and tilde rl (small vectors)

      rwork(1,y0) = -zone
      rwork(2:l,y0) = rwork(2:l,z)
      CALL zgetrs('n', l-1, 1, rwork(2:l,zz+1:zz+l-1), l-1, iwork, &
          rwork(2:l,y0), l-1, stat)
      rwork(l+1,y0) = zzero

      rwork(1,yl) = zzero
      rwork(2:l,yl) = rwork(2:l,z+l)
      CALL zgetrs ('n', l-1, 1, rwork(2:l,zz+1:zz+l-1), l-1, iwork, &
          rwork(2:l,yl), l-1, stat)
      rwork(l+1,yl) = -zone

      ! --- Convex combination

      CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
          rwork(1:l+1,y0), 1, zzero, rwork(1:l+1,y), 1)
      kappa0 = SQRT( ABS(zdotc(l+1, rwork(1:l+1,y0), 1, rwork(1:l+1,y), 1)) )
      CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
          rwork(1:l+1,yl), 1, zzero, rwork(1:l+1,y), 1)
      kappal = SQRT( ABS(zdotc(l+1, rwork(1:l+1,yl), 1, rwork(1:l+1,y), 1)) )
      CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
          rwork(1:l+1,y0), 1, zzero, rwork(1:l+1,y), 1)
      varrho = zdotc(l+1, rwork(1:l+1,yl), 1, rwork(1:l+1,y), 1) / &
          (kappa0*kappal)
      hatgamma = varrho/ABS(varrho) * MAX(ABS(varrho),7d-1) * &
          kappa0/kappal
      rwork(1:l+1,y0) = rwork(1:l+1,y0) - hatgamma * rwork(1:l+1,yl)

      !  --- Update

      omega = rwork(l+1,y0)
      DO j=1,l
        work(1:n,u) = work(1:n,u) - rwork(j+1,y0) * work(1:n,u+j)
        x(1:n) = x(1:n) + rwork(j+1,y0) * work(1:n,r+j-1)
        work(1:n,r) = work(1:n,r) - rwork(j+1,y0) * work(1:n,r+j)
      ENDDO


      CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
          rwork(1:l+1,y0), 1, zzero, rwork(1:l+1,y), 1)
      rnrm = SQRT( ABS(zdotc(l+1, rwork(1:l+1,y0), 1, rwork(1:l+1,y), 1)) )

      !---------------------------------------
      !  --- The reliable update part ---
      !---------------------------------------

      IF (.TRUE.) THEN
        mxnrmx = MAX (mxnrmx, rnrm)
        mxnrmr = MAX (mxnrmr, rnrm)
        xpdt = (rnrm < delta*rnrm0 .AND. rnrm0 < mxnrmx)
        rcmp = ((rnrm < delta*mxnrmr .AND. rnrm0 < mxnrmr) .OR. xpdt)
        IF (rcmp) THEN
          ! PRINT *, 'Performing residual update...'
          t(1:n) = x(1:n)
          !CALL CRS_ComplexLUSolve2( n, A, t )
          CALL DiagonalVelocityPreconditioning( A, t, dim)
          !CALL ComplexMatrixVectorProduct( A, t, work(1:n,r) )
          CALL ComplexMatrixVelocityVectorProduct( A, t, work(1:n,r), dim )
          work(1:n,r) = work(1:n,bp) - work(1:n,r)
          mxnrmr = rnrm
          IF (xpdt) THEN
            ! PRINT *, 'Performing solution update...'
            work(1:n,xp) = work(1:n,xp) + t(1:n)
            x(1:n) = zzero
            work(1:n,bp) = work(1:n,r)
            mxnrmx = rnrm
          ENDIF
        ENDIF
      ENDIF

      IF (rcmp) THEN
        IF (xpdt) THEN
          t(1:n) = work(1:n,xp)
        ELSE
          t(1:n) = t(1:n) + work(1:n,xp)
        END IF
      ELSE
        t(1:n) = x(1:n)
        !CALL CRS_ComplexLUSolve2( n, A, t )
        CALL DiagonalVelocityPreconditioning( A, t, dim)
        t(1:n) =  t(1:n) + work(1:n,xp)
      END IF

      ! bw_errorind = StoppingCriterion( n, A, t, b, work(1:n,r) )
      errorind = rnrm/bnrm
      WRITE(*,'(I4,ES12.3,ES12.3)') Round, errorind

      !IF (BackwardError) errorind = bw_errorind

    END DO

    !------------------------------------------------------------
    ! We have solved z = P*x, so finally solve the true unknown x
    !------------------------------------------------------------
    !CALL CRS_ComplexLUSolve2( n, A, x )
    CALL DiagonalVelocityPreconditioning( A, x, dim)
    x(1:n) = x(1:n) + work(1:n,xp)

    !WRITE(*,'(a,ES12.3)') 'An approximate lower bound for the condition number ', &
    !    ConditionEstimate( n, A, x, b )
    WRITE(*,'(a,ES12.3)') 'The 2-norm of the solution: ', &
        ComplexNorm( n, x )

!------------------------------------------------------------------------------
  END SUBROUTINE VelocitySolve
!------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  SUBROUTINE ComplexMatrixVelocityVectorProduct( K,u,v,dim )
!------------------------------------------------------------------------------
!
!   The computation of a specific matrix-vector product needed in preconditioning.
!   The subroutine computes the matrix-vector product v = Au where A is
!   the coefficient matrix for unknown velocities occupying
!   the (1,1) block of the coefficient matrix K.
!
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp), DIMENSION(*) :: u, v
    TYPE(Matrix_t), POINTER :: K
    INTEGER :: dim
!------------------------------------------------------------------------------
    INTEGER, POINTER :: Cols(:), Rows(:)
    REAL(KIND=dp), POINTER :: Values(:)
    INTEGER :: i, j, l, m, n, p, q, t
    COMPLEX(KIND=dp) :: s
!------------------------------------------------------------------------------

    n = K % NumberOfRows / 2
    q = n/(dim+2)

    Rows   => K % Rows
    Cols   => K % Cols
    Values => K % Values

    v(1:dim*q) = CMPLX( 0.0d0, 0.0d0, kind=dp)

    DO m=1,q
      DO p=1,dim
        i = (m-1)*(dim+2)+p
        DO l = 1,dim
          DO j = Rows(2*i-1)+2*(l-1), Rows(2*i)-1, 2*(dim+2)
            s = CMPLX( Values(j), -Values(j+1), kind=dp )
            t = (Cols(j)+1)/2
            v((m-1)*dim+p) = v((m-1)*dim+p) + s * u( (t-l)/(dim+2)*dim+p )
          END DO
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------------
  END SUBROUTINE ComplexMatrixVelocityVectorProduct
!------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
  SUBROUTINE DiagonalVelocityPreconditioning( K, v, dim )
!------------------------------------------------------------------------------
!
!   This subroutine performs a specific diagonal preconditioning needed in
!   solving the velocities. The subroutine computes the matrix-vector product
!   v = inv(P)* v where P is the diagonal matrix consisting of the diagonal
!   entries of the coefficient matrix for unknown velocities occupying
!   the (1,1) block of the coefficient matrix K.
!
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp), DIMENSION(*) :: v
    TYPE(Matrix_t), POINTER :: K
    INTEGER :: dim
!------------------------------------------------------------------------------
    INTEGER, POINTER :: Diag(:)
    REAL(KIND=dp), POINTER :: Values(:)
    INTEGER :: i, j, n, p, q
    COMPLEX(KIND=dp) :: s
!------------------------------------------------------------------------------

    n = K % NumberOfRows / 2
    q = n/(dim+2)

    Diag   => K % Diag
    Values => K % Values

    DO j=1,q
      DO p=1,dim
        i = (j-1)*(dim+2)+p
        s = CMPLX( Values(Diag(2*i-1)), -Values(Diag(2*i-1)+1), kind=dp )
        v((j-1)*dim+p) = v((j-1)*dim+p) / s
      END DO
    END DO

!-----------------------------------------------------------------------------
  END SUBROUTINE DiagonalVelocityPreconditioning
!------------------------------------------------------------------------------
  END SUBROUTINE aca_block_solve
!------------------------------------------------------------------------------
! Verbatim-kept assembly kernels for the BlockPreconditioning / slip /
! interface-only paths (same statements as the original internals; the
! author's host variables travel as explicit trailing arguments).
!------------------------------------------------------------------------------
   SUBROUTINE CoupledVelocityMatrix( StiffMatrix, Viscosity, AngularFrequency, &
       Density, Element, n, dim, CoordSysIn)
!------------------------------------------------------------------------------
!    This subroutine computes the element matrix for the preconditioning
!    equation of velocities. The coordinate system can be either orthogonal
!    Cartesian or axially symmetric.
!-------------------------------------------------------------------------------
     REAL(KIND=dp), TARGET :: StiffMatrix(:,:)
     REAL(KIND=dp) :: Viscosity(:), AngularFrequency, Density(:)
     INTEGER :: dim, n, CoordSysIn
     TYPE(Element_t), TARGET :: Element
!------------------------------------------------------------------------------
     COMPLEX(kind=dp) :: CStiff(dim*n,dim*n)
     REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), DetJ, r, &
         mu, rho0, s
     LOGICAL :: Stat
     INTEGER :: t, i, j, p, q
     TYPE(GaussIntegrationPoints_t) :: IP
     TYPE(Nodes_t) :: Nodes
     SAVE Nodes
!------------------------------------------------------------------------------
     CALL GetElementNodes( Nodes )
     StiffMatrix = 0.0d0
     CStiff = CMPLX(0.0d0,0.0d0, KIND=dp)
     !-------------------------
     ! Numerical integration:
     !-------------------------
     IP = GaussPoints( Element,n)
     DO t=1,IP % n
       !--------------------------------------------------------------
       ! Basis function values & derivatives at the integration point:
       !--------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t), detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

       s = IP % s(t) * detJ
       IF (CoordSysIn == AxisSymmetric) THEN
         r = SUM( Basis * Nodes % x(1:n) )
         s = r * s
       END IF

       !-----------------------------------------------
       ! Material parameters at the integration point:
       !----------------------------------------------
       mu  = SUM( Basis(1:n) * Viscosity(1:n) )
       rho0 = SUM( Density(1:n) * Basis(1:n) )
       !---------------------------------------------
       ! the stiffness matrix...
       !---------------------------------------------
       DO i=1,dim
         DO p=1,n
           DO q=1,n
             DO j = 1,dim
               CStiff((p-1)*dim+i,(q-1)*dim+i) = CStiff((p-1)*dim+i,(q-1)*dim+i) + &
                   s * CMPLX( 0.0d0, -mu/(AngularFrequency*rho0), KIND=dp ) * &
                   dBasisdx(q,j) * dBasisdx(p,j)
               CStiff((p-1)*dim+i, (q-1)*dim+j) = &
                   CStiff((p-1)*dim+i, (q-1)*dim+j) + &
                   CMPLX( 0.0d0, -mu/(AngularFrequency*rho0), KIND=dp ) * &
                   dBasisdx(q,i) * dBasisdx(p,j) * s
             END DO

             IF ( (i==1) .AND. (CoordSysIn == AxisSymmetric) ) THEN
               CStiff((p-1)*dim+i, (q-1)*dim+i) = &
                   CStiff((p-1)*dim+i, (q-1)*dim+i) + &
                   CMPLX( 0.0d0, -2*mu/(AngularFrequency*rho0), KIND=dp ) * 1/r**2 * &
                   Basis(q) * Basis(p) * s
             END IF

             CStiff((p-1)*dim+i,(q-1)*dim+i) = CStiff((p-1)*dim+i,(q-1)*dim+i) + &
                 s * CMPLX(1.0d0, 0.0d0, KIND=dp ) * Basis(q) * Basis(p)
           END DO
         END DO
       END DO
     END DO


     DO p=1,n
       DO i=1,DIM
         DO q=1,n
           DO j=1,DIM
             StiffMatrix( 2*DIM*(p-1)+2*i-1, 2*DIM*(q-1)+2*j-1 ) =  &
                 REAL( CSTIFF(DIM*(p-1)+i,DIM*(q-1)+j) )
             StiffMatrix( 2*DIM*(p-1)+2*i-1, 2*DIM*(q-1)+2*j ) =  &
                 -AIMAG( CSTIFF(DIM*(p-1)+i,DIM*(q-1)+j) )
             StiffMatrix( 2*DIM*(p-1)+2*i, 2*DIM*(q-1)+2*j-1 ) =  &
                 AIMAG( CSTIFF(DIM*(p-1)+i,DIM*(q-1)+j) )
             StiffMatrix( 2*DIM*(p-1)+2*i, 2*DIM*(q-1)+2*j ) =  &
                 REAL( CSTIFF(DIM*(p-1)+i,DIM*(q-1)+j) )
           END DO
         END DO
       END DO
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE CoupledVelocityMatrix
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
   SUBROUTINE PressureLaplaceMatrix( StiffMatrix, Viscosity, AngularFrequency, &
       Density, Element, n, dim, CoordSysIn)
!------------------------------------------------------------------------------
!   This subroutine computes the element matrix for the Laplacian type
!   term arising in the consistent splitting approach. The coordinate system can
!   be either orthogonal Cartesian or axially symmetric.
!-------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: StiffMatrix(:,:)
    REAL(KIND=dp) :: Viscosity(:), AngularFrequency, Density(:)
    INTEGER :: dim, n, CoordSysIn
    TYPE(Element_t), TARGET :: Element
!------------------------------------------------------------------------------
    COMPLEX(kind=dp) :: CStiff(n,n)
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), DetJ, r, &
        mu, rho0, s
    LOGICAL :: Stat
    INTEGER :: t, i, j, p, q
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    StiffMatrix = 0.0d0
    CStiff = CMPLX(0.0d0,0.0d0, KIND=dp)
    !-------------------------
    ! Numerical integration:
    !-------------------------
    IP = GaussPoints( Element,n)
    DO t=1,IP % n
      !--------------------------------------------------------------
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

      s = IP % s(t) * detJ
      IF (CoordSysIn == AxisSymmetric) THEN
        r = SUM( Basis * Nodes % x(1:n) )
        s = r * s
      END IF

      !-----------------------------------------------
      ! Material parameters at the integration point:
      !----------------------------------------------
      mu  = SUM( Basis(1:n) * Viscosity(1:n) )
      rho0 = SUM( Density(1:n) * Basis(1:n) )
      !---------------------------------------------
      ! the stiffness matrix...
      !---------------------------------------------
      DO p=1,n
        DO q=1,n
          DO j = 1,dim
            CStiff(p,q) = CStiff(p,q) - &
                s * CMPLX( 0.0d0, 2.0d0*mu/(AngularFrequency*rho0), KIND=dp ) * &
                dBasisdx(q,j) * dBasisdx(p,j)
          END DO
        END DO
      END DO
    END DO

    DO i=1,n
      DO j=1,n
        StiffMatrix( 2*(i-1)+1, 2*(j-1)+1 ) =  REAL( CStiff(i,j) )
        StiffMatrix( 2*(i-1)+1, 2*(j-1)+2 ) = -AIMAG( CStiff(i,j) )
        StiffMatrix( 2*(i-1)+2, 2*(j-1)+1 ) =  AIMAG( CStiff(i,j) )
        StiffMatrix( 2*(i-1)+2, 2*(j-1)+2 ) =  REAL( CStiff(i,j) )
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE PressureLaplaceMatrix
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE PressureMassMatrix( StiffMatrix, Element, n, dim, CoordSysIn)
!------------------------------------------------------------------------------
!   This subroutine computes the element matrix for the mass matrix which is
!   needed in the consistent splitting approach to compute the continuous
!   approximation of residual. The coordinate system can be either orthogonal
!   Cartesian or axially symmetric.
!-------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: StiffMatrix(:,:)
    INTEGER :: dim, n, CoordSysIn
    TYPE(Element_t), TARGET :: Element
!------------------------------------------------------------------------------
    COMPLEX(kind=dp) :: CStiff(n,n)
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), DetJ, r, s
    LOGICAL :: Stat
    INTEGER :: t, i, j, p, q
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    StiffMatrix = 0.0d0
    CStiff = CMPLX(0.0d0,0.0d0, KIND=dp)
    !-------------------------
    ! Numerical integration:
    !-------------------------
    IP = GaussPoints( Element,n)
    DO t=1,IP % n
      !--------------------------------------------------------------
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

      s = IP % s(t) * detJ
      IF (CoordSysIn == AxisSymmetric) THEN
        r = SUM( Basis * Nodes % x(1:n) )
        s = r * s
      END IF

      !---------------------------------------------
      ! the stiffness matrix...
      !---------------------------------------------
      DO p=1,n
        DO q=1,n
          CStiff(p,q) = CStiff(p,q) + s * CMPLX(1.0d0, 0.0d0, KIND=dp) * &
              Basis(q) * Basis(p)
        END DO
      END DO
    END DO

    DO i=1,n
      DO j=1,n
        StiffMatrix( 2*(i-1)+1, 2*(j-1)+1 ) =  REAL( CStiff(i,j) )
        StiffMatrix( 2*(i-1)+1, 2*(j-1)+2 ) = -AIMAG( CStiff(i,j) )
        StiffMatrix( 2*(i-1)+2, 2*(j-1)+1 ) =  AIMAG( CStiff(i,j) )
        StiffMatrix( 2*(i-1)+2, 2*(j-1)+2 ) =  REAL( CStiff(i,j) )
      END DO
    END DO
!------------------------------------------------------------------------------
    END SUBROUTINE PressureMassMatrix
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
    SUBROUTINE SchurComplementMatrix( StiffMatrix, AngularFrequency , &
        SpecificHeat, HeatRatio, Density, Pressure,               &
        Temperature, Conductivity, Viscosity, Lambda,             &
        Element, n, dim, CoordSysIn)
!------------------------------------------------------------------------------
!     This subroutine computes the element matrix for the preconditioning
!     equation of temperature and pressure. The coordinate system can be either
!     orthogonal Cartesian or axially symmetric.
!-------------------------------------------------------------------------------
      REAL(KIND=dp), TARGET :: StiffMatrix(:,:)
      REAL(KIND=dp) :: AngularFrequency, SpecificHeat(:), HeatRatio(:), Density(:), &
          Pressure(:),  Temperature(:), Conductivity(:), Viscosity(:), Lambda(:)
      TYPE(Element_t), TARGET :: Element
      INTEGER :: n, dim, CoordSysIn
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: Basis(n),dBasisdx(n,3), ddBasisddx(n,3,3), DetJ, r, &
          CV, kappa, mu, rho0, P0, gamma, la, K1, K2, s
      COMPLEX(kind=dp) :: CStiff(2*n,2*n), C1, C2, C3, C4
      LOGICAL :: Stat
      INTEGER :: t, i, j, p, q
      TYPE(GaussIntegrationPoints_t) :: IP
      TYPE(Nodes_t) :: Nodes
      SAVE Nodes
!------------------------------------------------------------------------------
      CALL GetElementNodes( Nodes )
      StiffMatrix = 0.0d0
      CStiff = CMPLX( 0.0d0,0.0d0, kind=dp )
      !-------------------------
      ! Numerical integration:
      !-------------------------
      IP = GaussPoints( Element,n)
      DO t=1,IP % n
        !--------------------------------------------------------------
        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

        s = IP % s(t) * detJ
        IF (CoordSysIn == AxisSymmetric) THEN
          r = SUM( Basis * Nodes % x(1:n) )
          s = r * s
        END IF

        !-----------------------------------------------
        ! Material parameters at the integration point:
        !----------------------------------------------
        CV = SUM( SpecificHeat(1:n) * Basis(1:n) )
        kappa = SUM( Conductivity(1:n) * Basis(1:n) )
        mu  = SUM( Basis(1:n) * Viscosity(1:n) )
        rho0  = SUM( Basis(1:n) * Density(1:n) )
        P0  = SUM( Basis(1:n) * Pressure(1:n) )
        gamma  = SUM( Basis(1:n) * HeatRatio(1:n) )
        la = SUM( Basis(1:n) * Lambda(1:n) )

        C1 = CMPLX( 1.0d0,AngularFrequency/P0*(2.0d0*mu+la), kind=dp ) / &
            CMPLX( 1.0d0,AngularFrequency/P0*(la), kind=dp )

        C2 = CMPLX( rho0*AngularFrequency,0.0d0, kind=dp) / &
            CMPLX( P0/AngularFrequency, la, kind=dp )
        C3 = CMPLX( 1.0d0, 0.0d0, kind=dp)
        C4 = CMPLX( 1.0d0,0.0d0, kind=dp) / &
            CMPLX( 1.0d0, AngularFrequency/P0*(la), kind=dp )

        K1 = kappa/(rho0*AngularFrequency*(gamma-1.0d0)*CV)
        K2 = 1.0d0/(gamma-1.0d0)

        !---------------------------------------------
        ! the stiffness matrix...
        !----------------------------------------
        DO p=1,n
          DO q=1,n
            DO i=1,dim

              CStiff((p-1)*2+1,(q-1)*2+1) = CStiff((p-1)*2+1,(q-1)*2+1) + &
                  s * CMPLX( 0.0d0, K1, kind=dp) * dBasisdx(q,i) * dBasisdx(p,i)

              CStiff((p-1)*2+2,(q-1)*2+1) = CStiff((p-1)*2+2,(q-1)*2+1) + &
                  s * C3 * dBasisdx(q,i) * dBasisdx(p,i)

              CStiff((p-1)*2+2,(q-1)*2+2) = CStiff((p-1)*2+2,(q-1)*2+2) + &
                  s * C1 * dBasisdx(q,i) * dBasisdx(p,i)
            END DO
            CStiff((p-1)*2+1,(q-1)*2+1) = CStiff((p-1)*2+1,(q-1)*2+1) + &
                s * CMPLX(-K2, 0.0d0, kind=dp) * Basis(q) * Basis(p)

            CStiff((p-1)*2+1,(q-1)*2+2) = CStiff((p-1)*2+1,(q-1)*2+2) + &
                s * C4 * Basis(q) * Basis(p)

            CStiff((p-1)*2+2,(q-1)*2+2) = CStiff((p-1)*2+2,(q-1)*2+2) - &
                s * C2 * Basis(q) * Basis(p)
          END DO
        END DO
      END DO

      DO p=1,n
        DO i=1,2
          DO q=1,n
            DO j=1,2
              StiffMatrix( 4*(p-1)+2*i-1, 4*(q-1)+2*j-1 ) =  &
                  REAL( CStiff(2*(p-1)+i,2*(q-1)+j) )
              StiffMatrix( 4*(p-1)+2*i-1, 4*(q-1)+2*j ) =  &
                  -AIMAG( CStiff(2*(p-1)+i,2*(q-1)+j) )
              StiffMatrix( 4*(p-1)+2*i, 4*(q-1)+2*j-1 ) =  &
                  AIMAG( CStiff(2*(p-1)+i,2*(q-1)+j) )
              StiffMatrix( 4*(p-1)+2*i, 4*(q-1)+2*j ) =  &
                  REAL( CStiff(2*(p-1)+i,2*(q-1)+j) )
            END DO
          END DO
        END DO
      END DO
!------------------------------------------------------------------------------
    END SUBROUTINE SchurComplementMatrix
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  SUBROUTINE VelocityImpedanceMatrix( StiffMatrix, AngularFrequency, &
      Density, Impedance, Element, n, Nodes, dim, CoordSysIn)
!------------------------------------------------------------------------------
!   This subroutine computes the contribution of the impedance boundary
!   condition to the preconditioning equation of velocities. The coordinate
!   system can be either orthogonal Cartesian or axially symmetric.
!-------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: StiffMatrix(:,:)
    REAL(KIND=dp) :: AngularFrequency, Density(:), Impedance(:,:)
    INTEGER :: dim, n, CoordSysIn
    TYPE(Element_t), TARGET :: Element
    TYPE(Nodes_t) :: Nodes
!------------------------------------------------------------------------------
    COMPLEX(kind=dp) :: CStiff(dim*n,dim*n)
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), DetJ, r, &
        Normal(3), Impedance1, Impedance2
    LOGICAL :: Stat
    INTEGER :: t, i, j, k, l, p, q
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: rho0, s
!------------------------------------------------------------------------------
    StiffMatrix = 0.0d0
    CStiff = CMPLX(0.0d0,0.0d0, kind=dp)
    !-------------------------
    ! Numerical integration:
    !-------------------------
    IP = GaussPoints( Element,n)
    DO t=1,IP % n
      !--------------------------------------------------------------
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )
      s = IP % s(t) * detJ
      IF (CoordSysIn == AxisSymmetric) THEN
        r = SUM( Basis * Nodes % x(1:n) )
        s = r * s
      END IF
      !-----------------------------------------------
      ! Material parameters etc. at the integration point:
      !----------------------------------------------
      rho0 = SUM( Density(1:n) * Basis(1:n) )
      Normal = Normalvector(Element, Nodes, IP % U(t), IP % V(t), .TRUE.)
      Impedance1 = 1.0d0/(AngularFrequency*rho0) * SUM( Impedance(1,1:n) * Basis(1:n) )
      Impedance2 = 1.0d0/(AngularFrequency*rho0) * SUM( Impedance(2,1:n) * Basis(1:n) )
      !---------------------------------------------
      ! the stiffness matrix...
      !---------------------------------------------
      DO p=1,n
        DO i=1,dim
          DO q=1,n
            DO j=1,dim
              CStiff( (p-1)*DIM+i, (q-1)*DIM+j) = &
                  CStiff( (p-1)*DIM+i, (q-1)*DIM+j) + &
                  CMPLX(-Impedance2, Impedance1, kind=dp) * &
                  Basis(q) * Normal(j) * Basis(p) * Normal(i) * s
            END DO
          END DO
        END DO
      END DO
    END DO   ! Loop over integration points

    DO p=1,n
      DO i=1,DIM
        DO q=1,n
          DO j=1,DIM
            StiffMatrix( 2*DIM*(p-1)+2*i-1, 2*DIM*(q-1)+2*j-1 ) =  &
                REAL( CSTIFF(DIM*(p-1)+i,DIM*(q-1)+j) )
            StiffMatrix( 2*DIM*(p-1)+2*i-1, 2*DIM*(q-1)+2*j ) =  &
                -AIMAG( CSTIFF(DIM*(p-1)+i,DIM*(q-1)+j) )
            StiffMatrix( 2*DIM*(p-1)+2*i, 2*DIM*(q-1)+2*j-1 ) =  &
                AIMAG( CSTIFF(DIM*(p-1)+i,DIM*(q-1)+j) )
            StiffMatrix( 2*DIM*(p-1)+2*i, 2*DIM*(q-1)+2*j ) =  &
                REAL( CSTIFF(DIM*(p-1)+i,DIM*(q-1)+j) )
          END DO
        END DO
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE VelocityImpedanceMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE SchurComplementImpedanceMatrix( StiffMatrix, Impedance, &
      AngularFrequency, SpecificHeat, HeatRatio, Density, Conductivity, &
      Element, n, Nodes, dim, CoordSysIn)
!------------------------------------------------------------------------------
!   This subroutine computes the contribution of the impedance boundary
!   conditions to the preconditioning equation of temperature and pressure.
!   The coordinate system can be either orthogonal Cartesian or axially symmetric.
!-------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: StiffMatrix(:,:)
    REAL(KIND=dp) :: Impedance(:,:), AngularFrequency, SpecificHeat(:), &
        HeatRatio(:), Density(:), Conductivity(:)
    INTEGER :: dim, n, CoordSysIn
    TYPE(Element_t), TARGET :: Element
    TYPE(Nodes_t) :: Nodes
!------------------------------------------------------------------------------
    COMPLEX(kind=dp) :: CStiff(2*n,2*n), K1, ZT
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), DetJ, r, &
        Normal(3), Impedance1, Impedance2
    LOGICAL :: Stat
    INTEGER :: t, i, j, k, l, p, q
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: CV, kappa, rho0, gamma, s
!------------------------------------------------------------------------------
    StiffMatrix = 0.0d0
    CStiff = CMPLX(0.0d0,0.0d0, kind=dp)
    !-------------------------
    ! Numerical integration:
    !-------------------------
    IP = GaussPoints( Element,n)
    DO t=1,IP % n
      !--------------------------------------------------------------
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

      s = IP % s(t) * detJ
      IF (CoordSysIn == AxisSymmetric) THEN
        r = SUM( Basis * Nodes % x(1:n) )
        s = r * s
      END IF

      !-----------------------------------------------
      ! Material parameters at the integration point:
      !----------------------------------------------
      CV = SUM( SpecificHeat(1:n) * Basis(1:n) )
      kappa = SUM( Conductivity(1:n) * Basis(1:n) )
      rho0  = SUM( Basis(1:n) * Density(1:n) )
      gamma  = SUM( Basis(1:n) * HeatRatio(1:n) )

      K1 = CMPLX( 0.0d0, kappa/(rho0*AngularFrequency*(gamma-1.0d0)*CV), kind=dp )
      ZT = CMPLX( SUM( Impedance(3,1:n) * Basis(1:n) ), SUM( Impedance(4,1:n) * Basis(1:n) ), kind=dp )

      DO p=1,n
        DO q=1,n
          CStiff((p-1)*2+1,(q-1)*2+1) = CStiff((p-1)*2+1,(q-1)*2+1) - &
              s * K1 * ZT * Basis(q) * Basis(p)
        END DO
      END DO
    END DO

    DO p=1,n
      DO i=1,2
        DO q=1,n
          DO j=1,2
            StiffMatrix( 4*(p-1)+2*i-1, 4*(q-1)+2*j-1 ) =  &
                REAL( CStiff(2*(p-1)+i,2*(q-1)+j) )
            StiffMatrix( 4*(p-1)+2*i-1, 4*(q-1)+2*j ) =  &
                -AIMAG( CStiff(2*(p-1)+i,2*(q-1)+j) )
            StiffMatrix( 4*(p-1)+2*i, 4*(q-1)+2*j-1 ) =  &
                AIMAG( CStiff(2*(p-1)+i,2*(q-1)+j) )
            StiffMatrix( 4*(p-1)+2*i, 4*(q-1)+2*j ) =  &
                REAL( CStiff(2*(p-1)+i,2*(q-1)+j) )
          END DO
        END DO
      END DO
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE SchurComplementImpedanceMatrix
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  SUBROUTINE VelocitySlipMatrix( StiffMatrix, SpecificHeat, &
      HeatRatio, Density, Temperature, &
      AngularFrequency, WallTemperature, SlipCoefficient1, &
      Element, n, Nodes, dim, CoordSysIn)
!------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: StiffMatrix(:,:)
    REAL(KIND=dp) :: SpecificHeat(:), HeatRatio(:), Density(:), &
        Temperature(:), AngularFrequency, WallTemperature(:), SlipCoefficient1
    INTEGER :: dim, n, CoordSysIn
    TYPE(Element_t), TARGET :: Element
    TYPE(Nodes_t) :: Nodes
!------------------------------------------------------------------------------
    COMPLEX(kind=dp) :: CStiff(dim*n,dim*n)
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), DetJ, r, &
        CV, gamma, rho0, T0, WallT0, C1, s, &
        Normal(3), Tangent1(3), Tangent2(3)
    LOGICAL :: Stat
    INTEGER :: t, i, j, k, l, p, q
    TYPE(GaussIntegrationPoints_t) :: IP
!------------------------------------------------------------------------------
    StiffMatrix = 0.0d0
    CStiff = CMPLX(0.0d0,0.0d0,kind=dp)
    !-------------------------
    ! Numerical integration:
    !-------------------------
    IP = GaussPoints( Element,n)
    DO t=1,IP % n
      !--------------------------------------------------------------
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )
      s = IP % s(t) * detJ
      IF (CoordSysIn == AxisSymmetric) THEN
        r = SUM( Basis * Nodes % x(1:n) )
        s = r * s
      END IF
      !-----------------------------------------------
      ! Material parameters etc. at the integration point:
      !----------------------------------------------
      CV = SUM( SpecificHeat(1:n) * Basis(1:n) )
      gamma = SUM( HeatRatio(1:n) * Basis(1:n) )
      rho0 = SUM( Density(1:n) * Basis(1:n) )
      T0 =  SUM( Temperature(1:n) * Basis(1:n) )
      WallT0 = SUM( WallTemperature(1:n) * Basis(1:n) )

      C1 = -SlipCoefficient1/(2.0d0-SlipCoefficient1) * &
          rho0*SQRT(2.0d0*(gamma-1.0d0)*CV* &
          (T0+WallT0)/PI)

      Normal = Normalvector(Element, Nodes, IP % U(t), IP % V(t), .TRUE.)
      CALL TangentDirections(Normal, Tangent1, Tangent2)

      DO p=1,n
        DO i=1,dim
          DO j=1,dim
            DO q=1,n
              CStiff( (p-1)*DIM+i, (q-1)*DIM+j) = &
                  CStiff( (p-1)*DIM+i, (q-1)*DIM+j) + &
                  CMPLX(0.0d0, C1, kind=dp) * 1.0d0/(AngularFrequency*rho0) * &
                  Basis(q) * Tangent1(j) * Basis(p) * Tangent1(i) * s
            END DO
          END DO
        END DO
      END DO

      IF (dim > 2) THEN
        DO p=1,n
          DO i=1,dim
            DO j=1,dim
              DO q=1,n
                CStiff( (p-1)*DIM+i, (q-1)*DIM+j) = &
                    CStiff( (p-1)*DIM+i, (q-1)*DIM+j) + &
                    CMPLX(0.0d0, C1, kind=dp) * 1.0d0/(AngularFrequency*rho0) * &
                    Basis(q) * Tangent2(j) * Basis(p) * Tangent2(i) * s
              END DO
            END DO
          END DO
        END DO
      END IF
    END DO

    DO p=1,n
      DO i=1,DIM
        DO q=1,n
          DO j=1,DIM
            StiffMatrix( 2*DIM*(p-1)+2*i-1, 2*DIM*(q-1)+2*j-1 ) =  &
                REAL( CSTIFF(DIM*(p-1)+i,DIM*(q-1)+j) )
            StiffMatrix( 2*DIM*(p-1)+2*i-1, 2*DIM*(q-1)+2*j ) =  &
                -AIMAG( CSTIFF(DIM*(p-1)+i,DIM*(q-1)+j) )
            StiffMatrix( 2*DIM*(p-1)+2*i, 2*DIM*(q-1)+2*j-1 ) =  &
                AIMAG( CSTIFF(DIM*(p-1)+i,DIM*(q-1)+j) )
            StiffMatrix( 2*DIM*(p-1)+2*i, 2*DIM*(q-1)+2*j ) =  &
                REAL( CSTIFF(DIM*(p-1)+i,DIM*(q-1)+j) )
          END DO
        END DO
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE VelocitySlipMatrix
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  SUBROUTINE SchurComplementSlipMatrix( StiffMatrix, SpecificHeat, &
      HeatRatio, Density, Temperature, AngularFrequency, Conductivity, &
      WallTemperature, SlipCoefficient2, &
      Element, n, Nodes, dim, CoordSysIn)
!------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: StiffMatrix(:,:)
    REAL(KIND=dp) :: SpecificHeat(:), HeatRatio(:), Density(:), Temperature(:), &
        AngularFrequency, Conductivity(:), WallTemperature(:), SlipCoefficient2
    INTEGER :: dim, n, CoordSysIn
    TYPE(Element_t), TARGET :: Element
    TYPE(Nodes_t) :: Nodes
!------------------------------------------------------------------------------
    COMPLEX(kind=dp) :: CStiff(2*n,2*n)
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), DetJ, r, &
        Normal(3), Impedance1, Impedance2
    LOGICAL :: Stat
    INTEGER :: t, i, j, k, l, p, q
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: CV, kappa, rho0, gamma, T0, WallT0, K1, C2, s
!------------------------------------------------------------------------------
    StiffMatrix = 0.0d0
    CStiff = CMPLX(0.0d0,0.0d0, kind=dp)
    !-------------------------
    ! Numerical integration:
    !-------------------------
    IP = GaussPoints( Element,n)
    DO t=1,IP % n
      !--------------------------------------------------------------
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

      s = IP % s(t) * detJ
      IF (CoordSysIn == AxisSymmetric) THEN
        r = SUM( Basis * Nodes % x(1:n) )
        s = r * s
      END IF

      !-----------------------------------------------
      ! Material parameters at the integration point:
      !----------------------------------------------
      CV = SUM( SpecificHeat(1:n) * Basis(1:n) )
      gamma = SUM( HeatRatio(1:n) * Basis(1:n) )
      rho0 = SUM( Density(1:n) * Basis(1:n) )
      kappa = SUM( Conductivity(1:n) * Basis(1:n) )
      T0 =  SUM( Temperature(1:n) * Basis(1:n) )
      WallT0 = SUM( WallTemperature(1:n) * Basis(1:n) )

      K1 = kappa/(rho0*AngularFrequency*(gamma-1.0d0)*CV)
      C2 = 1/kappa*SlipCoefficient2/(2.0d0-SlipCoefficient2) * &
          (gamma+1.0d0)/2.0d0 * rho0 * CV * &
          SQRT(2.0d0*(gamma-1.0d0)*CV*(T0+WallT0)/PI)

      DO p=1,n
        DO q=1,n
          CStiff((p-1)*2+1,(q-1)*2+1) = CStiff((p-1)*2+1,(q-1)*2+1) + &
              s * CMPLX(0.0d0, K1, kind=dp) * CMPLX(C2, 0.0d0, kind=dp) * Basis(q) * Basis(p)
        END DO
      END DO
    END DO

    DO p=1,n
      DO i=1,2
        DO q=1,n
          DO j=1,2
            StiffMatrix( 4*(p-1)+2*i-1, 4*(q-1)+2*j-1 ) =  &
                REAL( CStiff(2*(p-1)+i,2*(q-1)+j) )
            StiffMatrix( 4*(p-1)+2*i-1, 4*(q-1)+2*j ) =  &
                -AIMAG( CStiff(2*(p-1)+i,2*(q-1)+j) )
            StiffMatrix( 4*(p-1)+2*i, 4*(q-1)+2*j-1 ) =  &
                AIMAG( CStiff(2*(p-1)+i,2*(q-1)+j) )
            StiffMatrix( 4*(p-1)+2*i, 4*(q-1)+2*j ) =  &
                REAL( CStiff(2*(p-1)+i,2*(q-1)+j) )
          END DO
        END DO
      END DO
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE SchurComplementSlipMatrix
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  SUBROUTINE LocalInterfaceMatrix( StiffMatrix, Force, AngularFrequency, &
      SpecificHeat, HeatRatio, Density, Pressure,               &
      Temperature, Conductivity, Impedance, Load, &
      Element, n, Nodes, Dofs )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: StiffMatrix(:,:), Force(:), AngularFrequency, &
        SpecificHeat(:), HeatRatio(:), Density(:), Pressure(:),    &
        Temperature(:), Conductivity(:), Impedance(:,:), &
        Load(:,:)
    INTEGER :: n, Dofs
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), TARGET :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: SqrtElementMetric, U, V, W, S, Impedance1, Impedance2, &
        Impedance3, Impedance4, CV, gamma, rho0, P0, T0, kappa, K1, L(6),   &
        Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), X, Y, Z, Normal(3), r
    COMPLEX(KIND=dp) :: LSTIFF(n*(Dofs-2),n*(Dofs-2)), LFORCE(n*(Dofs-2))
    LOGICAL :: Stat
    INTEGER :: i, j, p, q, t, DIM, cSys
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
!------------------------------------------------------------------------------
    DIM = CoordinateSystemDimension()
    cSys = CurrentCoordinateSystem()
    LSTIFF = CMPLX(0.0d0,0.0d0, kind=dp)
    LFORCE = CMPLX(0.0d0,0.0d0, kind=dp)
!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------
    IntegStuff = GaussPoints( Element )
!------------------------------------------------------------------------------
    DO t=1, IntegStuff % n
      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)
      S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, U, V, W, SqrtElementMetric, &
          Basis, dBasisdx, ddBasisddx, .FALSE., .FALSE. )
      s = s * SqrtElementMetric

      IF (cSys == AxisSymmetric) THEN
        r = SUM( Basis * Nodes % x(1:n) )
        s = r * s
      END IF
!------------------------------------------------------------------------------
!     Problem parameters at the integration point
!------------------------------------------------------------------------------
      CV = SUM( SpecificHeat(1:n) * Basis(1:n) )
      gamma = SUM( HeatRatio(1:n) * Basis(1:n) )
      rho0 = SUM( Density(1:n) * Basis(1:n) )
      P0 = SUM( Pressure(1:n) * Basis(1:n) )
      T0 = SUM( Temperature(1:n) * Basis(1:n) )
      kappa = SUM( Conductivity(1:n) * Basis(1:n) )
!      K1 = kappa*AngularFrequency/(CV*(gamma-1.0d0)*P0)
      K1 = kappa/(rho0*AngularFrequency*(gamma-1.0d0)*CV)

      Normal = Normalvector(Element, Nodes, U, V, .TRUE.)

      DO j=1,dim
        L(2*j-1) = -1.0d0/(AngularFrequency*rho0) * SUM( Load(1,1:n) * Basis(1:n) ) * Normal(j)
        L(2*j) = -1.0d0/(AngularFrequency*rho0) * SUM( Load(2,1:n) * Basis(1:n) ) * Normal(j)
      END DO

      DO p=1,n
        DO i=1,dim
          LFORCE( (p-1)*(DIM+2)+i) = LFORCE( (p-1)*(DIM+2)+i ) - &
              s * Basis(p) * CMPLX( -L((i-1)*2+2), L((i-1)*2+1), kind=dp  )
        END DO
      END DO

    END DO   ! Loop over integration points

    DO p=1,n
      DO i=1,DIM+2
        Force( 2*(DIM+2)*(p-1)+2*i-1 ) = REAL(LFORCE( (p-1)*(DIM+2)+i ))
        Force( 2*(DIM+2)*(p-1)+2*i ) = AIMAG(LFORCE( (p-1)*(DIM+2)+i ))
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalInterfaceMatrix
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  SUBROUTINE SlipMatrix( StiffMatrix, Force, SpecificHeat, HeatRatio, &
      Density, Conductivity, Pressure, Temperature, &
      AngularFrequency, WallTemperature, &
      WallVelocity, SlipCoefficient1, SlipCoefficient2, SlipCoefficient3, &
      Element, n, Nodes, Dofs )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: StiffMatrix(:,:), Force(:), &
        SpecificHeat(:), HeatRatio(:), Density(:), Conductivity(:), &
        Pressure(:), Temperature(:), AngularFrequency, WallTemperature(:), &
        WallVelocity(:,:), SlipCoefficient1, SlipCoefficient2, SlipCoefficient3
    INTEGER :: n, Dofs
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), TARGET :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: SqrtElementMetric, U, V, W, S, Impedance1, Impedance2, &
        Impedance3, Impedance4, CV, gamma, rho0, P0, T0, ReV0(3), ImV0(3), &
        WallT0, kappa, K1, L(6),C1, C2, C3, &
        Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), X, Y, Z, Normal(3), r,  &
        Tangent1(3), Tangent2(3)
    COMPLEX(KIND=dp) :: LSTIFF(n*(Dofs-2),n*(Dofs-2)), LFORCE(n*(Dofs-2))
    LOGICAL :: Stat
    INTEGER :: i, j, p, q, t, DIM, cSys
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
!------------------------------------------------------------------------------
    DIM = CoordinateSystemDimension()
    cSys = CurrentCoordinateSystem()
    LSTIFF = CMPLX(0.0d0,0.0d0, kind=dp)
    LFORCE = CMPLX(0.0d0,0.0d0, kind=dp)
!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------
    IntegStuff = GaussPoints( Element )
!------------------------------------------------------------------------------
    DO t=1, IntegStuff % n
      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)
      S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, U, V, W, SqrtElementMetric, &
          Basis, dBasisdx, ddBasisddx, .FALSE., .FALSE. )
      IF (cSys == Cartesian) THEN
        s = s * SqrtElementMetric
      ELSE
        r = SUM( Basis * Nodes % x(1:n) )
        s = r * s * SqrtElementMetric
      END IF
!------------------------------------------------------------------------------
!     Problem parameters at the integration point
!------------------------------------------------------------------------------
      CV = SUM( SpecificHeat(1:n) * Basis(1:n) )
      gamma = SUM( HeatRatio(1:n) * Basis(1:n) )
      rho0 = SUM( Density(1:n) * Basis(1:n) )
      kappa = SUM( Conductivity(1:n) * Basis(1:n) )
      P0 =  SUM( Pressure(1:n) * Basis(1:n) )
      T0 =  SUM( Temperature(1:n) * Basis(1:n) )
      WallT0 = SUM( WallTemperature(1:n) * Basis(1:n) )

      DO i=1,3
        ReV0(i) = SUM( WallVelocity((i-1)*2+1,1:n) * Basis(1:n) )
        ImV0(i) = SUM( WallVelocity((i-1)*2+2,1:n) * Basis(1:n) )
      END DO

      K1 = kappa/(rho0*AngularFrequency*(gamma-1.0d0)*CV)

      C1 = -SlipCoefficient1/(2.0d0-SlipCoefficient1) * &
          rho0*SQRT(2.0d0*(gamma-1.0d0)*CV* &
          (T0+WallT0)/PI)

      C2 = 1/kappa*SlipCoefficient2/(2.0d0-SlipCoefficient2) * &
          (gamma+1.0d0)/2.0d0 * rho0 * CV * &
          SQRT(2.0d0*(gamma-1.0d0)*CV*(T0+WallT0)/PI)

      C3 = -SlipCoefficient3/(2.0d0-SlipCoefficient3) * &
          rho0*SQRT(2.0d0*(gamma-1.0d0)*CV* &
          (T0+WallT0)/PI)

      Normal = Normalvector(Element, Nodes, U, V, .TRUE.)

      CALL TangentDirections(Normal, Tangent1, Tangent2)

      DO p=1,n
        DO i=1,dim
          DO j=1,dim
            DO q=1,n
              LSTIFF( (p-1)*(DIM+2)+i, (q-1)*(DIM+2)+j) = &
                  LSTIFF( (p-1)*(DIM+2)+i, (q-1)*(DIM+2)+j) + &
                  CMPLX(0.0d0, C1, kind=dp) * 1.0d0/(AngularFrequency*rho0) * &
                  Basis(q) * Tangent1(j) * Basis(p) * Tangent1(i) * s
            END DO

            LFORCE( (p-1)*(DIM+2)+i) = LFORCE( (p-1)*(DIM+2)+i) + &
                CMPLX(0.0d0, C1, kind=dp) * 1.0d0/(AngularFrequency*rho0) * &
                CMPLX(ReV0(j),ImV0(j),kind=dp) * &
                Tangent1(j) * Basis(p) * Tangent1(i) * s
          END DO
        END DO
      END DO

      IF (dim > 2) THEN
        DO p=1,n
          DO i=1,dim
            DO j=1,dim
              DO q=1,n
                LSTIFF( (p-1)*(DIM+2)+i, (q-1)*(DIM+2)+j) = &
                    LSTIFF( (p-1)*(DIM+2)+i, (q-1)*(DIM+2)+j) + &
                    CMPLX(0.0d0, C1, kind=dp) * 1.0d0/(AngularFrequency*rho0) * &
                    Basis(q) * Tangent2(j) * Basis(p) * Tangent2(i) * s
              END DO

              LFORCE( (p-1)*(DIM+2)+i) = LFORCE( (p-1)*(DIM+2)+i) + &
                  CMPLX(0.0d0, C1, kind=dp) * 1.0d0/(AngularFrequency*rho0) * &
                  CMPLX(ReV0(j),ImV0(j),kind=dp) * &
                  Tangent2(j) * Basis(p) * Tangent2(i) * s
            END DO
          END DO
        END DO
      END IF

      !Normal direction

      DO p=1,n
         DO i=1,dim
            DO j=1,dim
               DO q=1,n
                  LSTIFF( (p-1)*(DIM+2)+i, (q-1)*(DIM+2)+j) = &
                       LSTIFF( (p-1)*(DIM+2)+i, (q-1)*(DIM+2)+j) + &
                       CMPLX(0.0d0, C3, kind=dp) * 1.0d0/(AngularFrequency*rho0) * &
                       Basis(q) * Normal(j) * Basis(p) * Normal(i) * s
               END DO

               LFORCE( (p-1)*(DIM+2)+i) = LFORCE( (p-1)*(DIM+2)+i) + &
                    CMPLX(0.0d0, C3, kind=dp) * 1.0d0/(AngularFrequency*rho0) * &
                    CMPLX(ReV0(j),ImV0(j),kind=dp) * &
                    Normal(j) * Basis(p) * Normal(i) * s
            END DO
         END DO
      END DO




      DO p=1,n
        DO q=1,n
          LSTIFF( (p-1)*(DIM+2)+DIM+1, (q-1)*(DIM+2)+DIM+1) = &
              LSTIFF( (p-1)*(DIM+2)+DIM+1, (q-1)*(DIM+2)+DIM+1) + &
              CMPLX(0.0d0, K1, kind=dp) * CMPLX(C2, 0.0d0, kind=dp) * &
              Basis(q) * Basis(p) * s
        END DO
        LFORCE( (p-1)*(DIM+2)+DIM+1) = LFORCE( (p-1)*(DIM+2)+DIM+1) + &
            CMPLX(0.0d0, K1, kind=dp) * CMPLX(C2, 0.0d0, kind=dp) * &
            WallT0 * Basis(p) * s
      END DO

    END DO  ! Loop over integration points

    DO p=1,n
      DO i=1,DIM+2
        Force( 2*(DIM+2)*(p-1)+2*i-1 ) = REAL(LFORCE( (p-1)*(DIM+2)+i ))
        Force( 2*(DIM+2)*(p-1)+2*i ) = AIMAG(LFORCE( (p-1)*(DIM+2)+i ))
        DO q=1,n
          DO j=1,DIM+2
            StiffMatrix( 2*(DIM+2)*(p-1)+2*i-1, 2*(DIM+2)*(q-1)+2*j-1 ) =  &
                REAL( LSTIFF((DIM+2)*(p-1)+i,(DIM+2)*(q-1)+j) )
            StiffMatrix( 2*(DIM+2)*(p-1)+2*i-1, 2*(DIM+2)*(q-1)+2*j ) =  &
                -AIMAG( LSTIFF((DIM+2)*(p-1)+i,(DIM+2)*(q-1)+j) )
            StiffMatrix( 2*(DIM+2)*(p-1)+2*i, 2*(DIM+2)*(q-1)+2*j-1 ) =  &
                AIMAG( LSTIFF((DIM+2)*(p-1)+i,(DIM+2)*(q-1)+j) )
            StiffMatrix( 2*(DIM+2)*(p-1)+2*i, 2*(DIM+2)*(q-1)+2*j ) =  &
                REAL( LSTIFF((DIM+2)*(p-1)+i,(DIM+2)*(q-1)+j) )
          END DO
        END DO
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE SlipMatrix
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  SUBROUTINE UpdateGlobalPreconditioner( StiffMatrix, LocalStiffMatrix, &
      n, NDOFs, NodeIndexes )
!------------------------------------------------------------------------------
!
! Add element matrices to global matrices
!
!------------------------------------------------------------------------------
     TYPE(Matrix_t), POINTER :: StiffMatrix

     REAL(KIND=dp) :: LocalStiffMatrix(:,:)

     INTEGER :: n, NDOFs, NodeIndexes(:)
!------------------------------------------------------------------------------
     INTEGER :: i,j,k
!------------------------------------------------------------------------------
!    Update global matrix .
!------------------------------------------------------------------------------
     SELECT CASE( StiffMatrix % FORMAT )
     CASE( MATRIX_CRS )
       CALL CRS_GlueLocalMatrix( StiffMatrix,n,NDOFs, NodeIndexes, &
           LocalStiffMatrix )
     END SELECT

!------------------------------------------------------------------------------
   END SUBROUTINE UpdateGlobalPreconditioner
!------------------------------------------------------------------------------






!------------------------------------------------------------------------------
   SUBROUTINE SetBoundaryConditions( Model, StiffMatrix, &
                   Name, DOF, NDOFs, Perm )
!------------------------------------------------------------------------------
!
! Set dirichlet boundary condition for given dof
!
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(Matrix_t), POINTER :: StiffMatrix

    CHARACTER(LEN=*) :: Name
    INTEGER :: DOF, NDOFs, Perm(:)
!------------------------------------------------------------------------------

    TYPE(Element_t), POINTER :: CurrentElement
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER :: i,j,k,n,t,k1,k2
    LOGICAL :: GotIt, periodic
    REAL(KIND=dp) :: Work(Model % MaxElementNodes),s

!------------------------------------------------------------------------------

    DO t = Model % NumberOfBulkElements + 1, &
        Model % NumberOfBulkElements + Model % NumberOfBoundaryElements

      CurrentElement => Model % Elements(t)
!------------------------------------------------------------------------------
!      Set the current element pointer in the model structure to
!      reflect the element being processed
!------------------------------------------------------------------------------
      Model % CurrentElement => Model % Elements(t)
!------------------------------------------------------------------------------
      n = CurrentElement % TYPE % NumberOfNodes
      NodeIndexes => CurrentElement % NodeIndexes(1:n)

      DO i=1,Model % NumberOfBCs
        IF ( CurrentElement % BoundaryInfo % Constraint == &
            Model % BCs(i) % Tag ) THEN

          Work(1:n) = ListGetReal( Model % BCs(i) % Values, &
              Name,n,NodeIndexes, gotIt )
          IF ( gotIt ) THEN
            DO j=1,n
              k = Perm(NodeIndexes(j))
              IF ( k > 0 ) THEN
                k = NDOFs * (k-1) + DOF
                s = 1.0d0
                CALL ZeroRow( StiffMatrix,k )
                CALL SetMatrixElement( StiffMatrix,k,k, 1.0d0 * s )
              END IF
            END DO
          END IF
        END IF
      END DO
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE SetBoundaryConditions
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  SUBROUTINE ComputeAcousticImpedance( AcousticI, Bndries, AcImpedances, &
      Model, Solver, Flow, FlowPerm, Dofs, VelocityComponents, ElementNodes )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: Flow(:)
    INTEGER :: FlowPerm(:), Dofs, VelocityComponents
    TYPE(Nodes_t) :: ElementNodes
    INTEGER :: AcousticI, Bndries(:)
    REAL(KIND=dp) :: AcImpedances(:,:)
!------------------------------------------------------------------------------

    TYPE(Element_t), POINTER :: CurrentElement
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
    REAL(KIND=dp), ALLOCATABLE :: PressureData(:,:), VelocityData(:)
    REAL(KIND=dp) :: ElementPresRe, ElementPresIm
    REAL(KIND=dp) :: ElementVeloRe, ElementVeloIm
    REAL(KIND=dp) :: u, v, w, s, xpos, ypos, zpos
    REAL(KIND=dp) :: SqrtMetric, SqrtElementMetric
    REAL(KIND=dp) :: Basis(Model % MaxElementNodes), Metric(3,3), &
         dBasisdx(Model % MaxElementNodes,3), Symb(3,3,3), &
         ddBasisddx(Model % MaxElementNodes,3,3), dSymb(3,3,3,3)
    REAL(KIND=dp) :: Normal(3), VReal(3), VIm(3), Preal, PIm
    REAL(KIND=dp) :: AbsNormalVelo, InPhasePres, OutPhasePres
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER :: i, j, k, n, t
    LOGICAL :: stat
!------------------------------------------------------------------------------

    ALLOCATE( PressureData( AcousticI, 2 ) )
    PressureData = 0.0d0
    ALLOCATE( VelocityData( 2 ) )
    VelocityData = 0.0d0

    VReal = 0.0d0
    VIm = 0.0d0

    DO t = Model % NumberOfBulkElements +1, &
         Model % NumberOfBulkElements + Model % NumberOfBoundaryElements

       CurrentElement => Solver % Mesh % Elements(t)
       n = CurrentElement % TYPE % NumberOfNodes
       NodeIndexes => CurrentElement % NodeIndexes

       DO i = 1, AcousticI
          IF ( CurrentElement % BoundaryInfo % Constraint == &
               Model % BCs(Bndries(i)) % Tag ) THEN

             ElementPresRe = 0.0d0
             ElementPresIm = 0.0d0
             ElementVeloRe = 0.0d0
             ElementVeloIm = 0.0d0

             ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
             ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
             ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)
!------------------------------------------------------------------------------
!      Numerical integration
!------------------------------------------------------------------------------

             IntegStuff = GaussPoints( CurrentElement )
             DO j = 1, IntegStuff % n
                u = IntegStuff % u(j)
                v = IntegStuff % v(j)
                w = IntegStuff % w(j)

!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------

                stat = ElementInfo( CurrentElement, ElementNodes, u, v, w, &
                     SqrtElementMetric, Basis, dBasisdx, ddBasisddx, &
                     .FALSE., .FALSE. )

!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
                s = 1.0d0
                IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
                   xpos = SUM( ElementNodes % x(1:n) * Basis(1:n) )
                   ypos = SUM( ElementNodes % y(1:n) * Basis(1:n) )
                   zpos = SUM( ElementNodes % z(1:n) * Basis(1:n) )
                   s = 2*PI
                END IF

                CALL CoordinateSystemInfo( Metric, SqrtMetric, Symb, dSymb, &
                     xpos, ypos, zpos)

                s = s * SqrtMetric * SqrtElementMetric * IntegStuff % s(j)

!------------------------------------------------------------------------------

                PReal = SUM( Flow( DOFs * ( FlowPerm(NodeIndexes) ) &
                     - 1 ) * Basis(1:n) )
                PIm = SUM( Flow( DOFs * ( FlowPerm(NodeIndexes) ) ) &
                     * Basis(1:n) )


                ElementPresRe = ElementPresRe + s * PReal
                ElementPresIm = ElementPresIm + s * PIm

!------------------------------------------------------------------------------
!   The normal velocity is needed only for the transmitting boundary i = 1
!------------------------------------------------------------------------------

                IF ( i == 1 ) THEN
                   Normal = Normalvector( CurrentElement, ElementNodes, &
                        u, v, .TRUE. )
                   Normal = -Normal   ! use inward normal vector

                   DO k = 1, VelocityComponents
                      VReal(k) = SUM( Flow( DOFs * ( FlowPerm(NodeIndexes)-1 ) &
                           + 2 * k - 1 ) * Basis(1:n) )
                      VIm(k) = SUM( Flow( DOFs * ( FlowPerm(NodeIndexes)-1 ) &
                           + 2 * k ) * Basis(1:n) )
                   END DO

                   ElementVeloRe = ElementVeloRe + s * &
                        SUM( VReal(1:3) * Normal(1:3) )
                   ElementVeloIm = ElementVeloIm + s * &
                        SUM( VIm(1:3) * Normal(1:3) )
                END IF

!------------------------------------------------------------------------------

             END DO

          PressureData(i,1) = PressureData(i,1) + ElementPresRe
          PressureData(i,2) = PressureData(i,2) + ElementPresIm
          VelocityData(1) = VelocityData(1) + ElementVeloRe
          VelocityData(2) = VelocityData(2) + ElementVeloIm

          END IF

       END DO
    END DO

!------------------------------------------------------------------------------

    DO i = 1, AcousticI
       AcImpedances(i,1) = ( PressureData(i,1) * VelocityData(1) + &
            PressureData(i,2) * VelocityData(2) ) / &
            SUM( VelocityData**2 )
       AcImpedances(i,2) = ( PressureData(i,2) * VelocityData(1) - &
            PressureData(i,1) * VelocityData(2) ) / &
            SUM( VelocityData**2 )
    END DO

    DEALLOCATE( PressureData )
    DEALLOCATE( VelocityData )

!------------------------------------------------------------------------------
  END SUBROUTINE ComputeAcousticImpedance
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  SUBROUTINE ComputeAverageVelocity( Element, n, WallVelo, res, Solver )
!------------------------------------------------------------------------------
    TYPE(Element_t), TARGET :: Element
    TYPE(Solver_t) :: Solver
    INTEGER :: n
    REAL(KIND=dp) :: WallVelo(:,:)
    COMPLEX(KIND=dp) :: res
!------------------------------------------------------------------------------
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    INTEGER :: t
    LOGICAL :: stat
    REAL(KIND=dp) :: AK, SqrtElementMetric,U,V,W,S
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3)
    REAL(KIND=dp) :: Normal(3)
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
    COMPLEX(KIND=dp) :: Velo(3), NormVelo
!------------------------------------------------------------------------------
    res = CMPLX(0.0d0,0.0d0, kind=dp)
    AK = ElementArea(Solver % Mesh, Element, n)
!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    IntegStuff = GaussPoints( Element )
!------------------------------------------------------------------------------
    DO t=1,IntegStuff % n
      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)
      S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, U, V, W, SqrtElementMetric, &
          Basis, dBasisdx, ddBasisddx, .FALSE., .FALSE. )

      s = s * SqrtElementMetric

      Normal = Normalvector(Element, Nodes, U, V, .TRUE.)
      !PRINT *, Normal

      Velo(1) = CMPLX( SUM( WallVelo(1,1:n) * Basis(1:n) ), SUM( WallVelo(2,1:n) * Basis(1:n) ), kind=dp )
      Velo(2) = CMPLX( SUM( WallVelo(3,1:n) * Basis(1:n) ), SUM( WallVelo(4,1:n) * Basis(1:n) ), kind=dp )
      Velo(3) = CMPLX( SUM( WallVelo(5,1:n) * Basis(1:n) ), SUM( WallVelo(6,1:n) * Basis(1:n) ), kind=dp )
      NormVelo = Velo(1) * CMPLX(Normal(1), 0.0d0, kind=dp) + Velo(2) * CMPLX(Normal(2), 0.0d0, kind=dp) + &
          Velo(3) * CMPLX(Normal(3), 0.0d0, kind=dp)
      res = res + s * NormVelo
!------------------------------------------------------------------------------
    END DO
    res = res / CMPLX(AK, 0.0d0, kind=dp)
!------------------------------------------------------------------------------
  END SUBROUTINE ComputeAverageVelocity
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  SUBROUTINE ExplicitStabilisationMatrix( StiffMatrix, Force, AngularFrequency , &
      SpecificHeat, HeatRatio, Density,               &
      Temperature, Conductivity, Viscosity, Lambda,             &
      HeatSource, Load, Bubbles, Element, n, Nodes, Dofs, NodeOnBoundary )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: StiffMatrix(:,:), Force(:), AngularFrequency, &
        SpecificHeat(:), HeatRatio(:), Density(:),    &
        Temperature(:), Conductivity(:), Viscosity(:), Lambda(:),  &
        HeatSource(:,:), Load(:,:)
    LOGICAL :: Bubbles, NodeOnBoundary(:)
    INTEGER :: n, Dofs
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), TARGET :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(2*n), dBasisdx(2*n,3), ddBasisddx(n,3,3)
    REAL(KIND=dp) :: SqrtElementMetric, U, V, W, S, L(6), &
        CV, gamma, rho0, P0, T0, kappa, mu, la, f1, f2, K1, K2, K3, k
    COMPLEX(KIND=dp) :: LSTIFF(n*(Dofs-2),n*(Dofs-2)), LFORCE(n*(Dofs-2)), A, &
        SchurConst, C1, C2, C3

    INTEGER :: i, j, p, q, t, DIM, NBasis, CoordSys, VelocityDofs, &
        VelocityComponents
    LOGICAL :: stat
    TYPE(GaussIntegrationPoints_t) :: IntegStuff

    REAL(KIND=dp) :: X, Y, Z, Metric(3,3), SqrtMetric, Symb(3,3,3), &
        dSymb(3,3,3,3), StabTerms(n), AK
!------------------------------------------------------------------------------
    DIM = CoordinateSystemDimension()

    Metric = 0.0d0
    Metric(1,1) = 1.0d0
    Metric(2,2) = 1.0d0
    Metric(3,3) = 1.0d0

    LSTIFF = CMPLX(0.0d0,0.0d0, kind=dp)
    LFORCE = CMPLX(0.0d0,0.0d0, kind=dp)
!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------
    NBasis = n
    IntegStuff = GaussPoints( Element )
!------------------------------------------------------------------------------
    DO t=1,IntegStuff % n
      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)
      S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, U, V, W, SqrtElementMetric, &
          Basis, dBasisdx, ddBasisddx, .FALSE., Bubbles )
      s = s * SqrtElementMetric
!------------------------------------------------------------------------------
!     Problem parameters and the real and imaginary part of the
!     load at the integration point
!------------------------------------------------------------------------------
      CV = SUM( SpecificHeat(1:n) * Basis(1:n) )
      gamma = SUM( HeatRatio(1:n) * Basis(1:n) )
      rho0 = SUM( Density(1:n) * Basis(1:n) )
      T0 = SUM( Temperature(1:n) * Basis(1:n) )
      kappa = SUM( Conductivity(1:n) * Basis(1:n) )
      mu = SUM( Viscosity(1:n) * Basis(1:n) )
      la = SUM( Lambda(1:n) * Basis(1:n) )

      P0 = (gamma-1.0d0)*CV*rho0*T0
      C1 = CMPLX( 1.0d0,AngularFrequency/P0*(2.0d0*mu+la), kind=dp ) / &
          CMPLX( 1.0d0,AngularFrequency/P0*(la), kind=dp )
      C2 = CMPLX( rho0*AngularFrequency,0.0d0, kind=dp) / &
          CMPLX( P0/AngularFrequency, la, kind=dp )
      C3 = CMPLX( 1.0d0, 0.0d0, kind=dp)

      DO p=1,N
        IF ( NodeOnBoundary(p) ) CYCLE
        DO q=1,N
          DO i=1,dim

            LSTIFF( p*(DIM+2), (q-1)*(DIM+2)+DIM+1 ) = &
                LSTIFF( p*(DIM+2), (q-1)*(DIM+2)+DIM+1 ) - &
                s * C3 * dBasisdx(q,i) * dBasisdx(p,i)

            LSTIFF( p*(DIM+2), q*(DIM+2) ) = &
                LSTIFF( p*(DIM+2), q*(DIM+2) ) - &
                s * C1 * dBasisdx(q,i) * dBasisdx(p,i)

          END DO

          LSTIFF( p*(DIM+2), q*(DIM+2) ) = &
              LSTIFF( p*(DIM+2), q*(DIM+2) ) + &
              s * C2 * Basis(q) * Basis(p)
        END DO
      END DO

    END DO


    DO p=1,n
      DO i=1,DIM+2
        DO q=1,n
          DO j=1,DIM+2
            StiffMatrix( 2*(DIM+2)*(p-1)+2*i-1, 2*(DIM+2)*(q-1)+2*j-1 ) =  &
                REAL( LSTIFF((DIM+2)*(p-1)+i,(DIM+2)*(q-1)+j) )
            StiffMatrix( 2*(DIM+2)*(p-1)+2*i-1, 2*(DIM+2)*(q-1)+2*j ) =  &
                -AIMAG( LSTIFF((DIM+2)*(p-1)+i,(DIM+2)*(q-1)+j) )
            StiffMatrix( 2*(DIM+2)*(p-1)+2*i, 2*(DIM+2)*(q-1)+2*j-1 ) =  &
                AIMAG( LSTIFF((DIM+2)*(p-1)+i,(DIM+2)*(q-1)+j) )
            StiffMatrix( 2*(DIM+2)*(p-1)+2*i, 2*(DIM+2)*(q-1)+2*j ) =  &
                REAL( LSTIFF((DIM+2)*(p-1)+i,(DIM+2)*(q-1)+j) )
          END DO
        END DO
      END DO
    END DO
    Force = 0.0d0
!--------------------------------------------------------------------------------
  END SUBROUTINE ExplicitStabilisationMatrix
!--------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
  SUBROUTINE SurfaceForceIntegration(Element, Parent, Traction, Moment, MomentAbout, Area, &
     CalculateMoment, Velo, Pressure, Viscosity, BulkViscosity, Nodes, ParentNodes, np, &
     Solver, CoordSys)
!----------------------------------------------------------------------------------
  TYPE(Element_t), TARGET :: Element, Parent
  TYPE(Solver_t) :: Solver
  INTEGER :: CoordSys
  REAL(kind=dp) :: Traction(6), Moment(6), MomentAbout(3), Area, Velo(:,:), &
      Pressure(:,:), Viscosity(:), BulkViscosity(:)
  LOGICAL :: CalculateMoment
  TYPE(Nodes_t) :: Nodes, ParentNodes
  INTEGER :: np
!----------------------------------------------------------------------------------
  LOGICAL :: stat
  INTEGER :: N_Integ, t, i, dim
  REAL(kind=dp) :: u, v, w, detJ, s, Basis(np), dBasisdx(np,3), Normal(3), &
      Visc, Lambda, ReGrad(3,3), ImGrad(3,3), ReDiv, ImDiv, ReD(3,3), ImD(3,3), &
      ReP, ImP, tmpmat(3,np), tmp = 0.0d0, r1, r2, r3, r
  REAL(KIND=dp), POINTER :: U_Integ(:), V_Integ(:), W_Integ(:), S_Integ(:)
  TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
  SAVE tmp
!----------------------------------------------------------------------------------
  dBasisdx = 0.0d0
  ReDiv = 0.0d0
  ImDiv = 0.0d0
  Traction = 0.0d0
  Area = 0.0d0
  dim = CoordinateSystemDimension()
!----------------------------------------------------------------------------------
! The normal vector to the boundary element
!----------------------------------------------------------------------------------
  SELECT CASE( Element % TYPE % NumberOfNodes )
  CASE(2)
    Normal = Normalvector(Element, Nodes, 0.0d0, 0.0d0, .TRUE.)
  CASE( 3 )
    Normal = Normalvector(Element, Nodes, 0.3d0, 0.3d0, .TRUE.)
  CASE( 4 )
    Normal = Normalvector(Element, Nodes, 0.0d0, 0.0d0, .TRUE.)
  END SELECT


!----------------------------------------------------------------------------------
! Compute field derivatives at the integration points of the parent element
!----------------------------------------------------------------------------------
  IntegStuff = Gausspoints(Parent)
  U_Integ => IntegStuff % u
  V_Integ => IntegStuff % v
  W_Integ => IntegStuff % w
  S_Integ => IntegStuff % s
  N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
  DO t=1,N_Integ
!------------------------------------------------------------------------------
    u = U_Integ(t)
    v = V_Integ(t)
    w = W_Integ(t)
!------------------------------------------------------------------------------
!   Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
    stat = ElementInfo( Parent, ParentNodes, u, v, w, &
        detJ, Basis, dBasisdx )

    s = detJ * S_Integ(t)

    IF (CoordSys == AxisSymmetric) THEN
       r = SUM( Basis(1:np) * ParentNodes % x(1:np) )
       s = r * s
    END IF

    Visc = SUM( Viscosity(1:np) * Basis(1:np) )
!    Visc = 0.0d0
    Lambda = SUM( BulkViscosity(1:np) * Basis(1:np) )
!    Lambda = 0.0d0
    ReP = SUM(Pressure(1,1:np) * Basis(1:np))
    ImP = SUM(Pressure(2,1:np) * Basis(1:np))
    tmpmat(1,1:np) = Velo(1,1:np)
    tmpmat(2,1:np) = Velo(3,1:np)
    tmpmat(3,1:np) = Velo(5,1:np)
    ReGrad = MATMUL( tmpmat, dBasisdx )
    tmpmat(1,1:np) = Velo(2,1:np)
    tmpmat(2,1:np) = Velo(4,1:np)
    tmpmat(3,1:np) = Velo(6,1:np)
    ImGrad = MATMUL( tmpmat, dBasisdx )

    ReD = Visc * ( ReGrad + TRANSPOSE(ReGrad) )
    ImD = Visc * ( ImGrad + TRANSPOSE(ImGrad) )

    ReDiv = 0.0d0
    ImDiv = 0.0d0
    DO i = 1,dim
      ReDiv = ReDiv + ReGrad(i,i)
      ImDiv = ImDiv + ImGrad(i,i)
    END DO

    IF (CoordSys == AxisSymmetric) THEN
       ReDiv = ReDiv + 1/r * SUM( Velo(1,1:np) * Basis(1:np) )
       ImDiv = ImDiv + 1/r * SUM( Velo(2,1:np) * Basis(1:np) )
    END IF

    Traction(1:3) = Traction(1:3) - ReP * Normal(1:3) + Lambda * &
        ReDiv * Normal(1:3) + MATMUL( ReD, Normal)

    Traction(4:6) = Traction(4:6) - ImP * Normal(1:3) + Lambda * &
        ImDiv * Normal(1:3) + MATMUL( ImD, Normal)

  END DO

  Area = ElementArea(Solver % Mesh, Element, Element % TYPE % NumberOfNodes)

  IF (CoordSys == AxisSymmetric) Area = 2*pi*Area

  Traction = -Area/N_Integ * Traction


  IF (CalculateMoment) THEN
    r1 = SUM( Nodes % x(1:Element % TYPE % NumberOfNodes) )/Element % TYPE % NumberOfNodes - MomentAbout(1)
    r2 = SUM( Nodes % y(1:Element % TYPE % NumberOfNodes) )/Element % TYPE % NumberOfNodes - MomentAbout(2)
    r3 = SUM( Nodes % z(1:Element % TYPE % NumberOfNodes) )/Element % TYPE % NumberOfNodes - MomentAbout(3)
    Moment(1) = r2 * Traction(3) - r3 * Traction(2)
    Moment(2) = r3 * Traction(1) - r1 * Traction(3)
    Moment(3) = r1 * Traction(2) - r2 * Traction(1)
    Moment(4) = r2 * Traction(6) - r3 * Traction(5)
    Moment(5) = r3 * Traction(4) - r1 * Traction(6)
    Moment(6) = r1 * Traction(5) - r2 * Traction(4)
  END IF

!----------------------------------------------------------------------------------
END SUBROUTINE SurfaceForceIntegration
!-----------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
SUBROUTINE SurfaceImpedanceIntegration(Element, Velo, Pressure, Nodes, n, Impedance1, Impedance2, &
     C3, C4, Area, CoordSys)
!---------------------------------------------------------------------------------------
  TYPE(Element_t), TARGET :: Element
  INTEGER :: CoordSys
  REAL(kind=dp) :: Velo(:,:), Pressure(:,:), Area
  TYPE(Nodes_t) :: Nodes
  INTEGER :: n
  COMPLEX(kind=dp) :: Impedance1, Impedance2, C3, C4
!------------------------------------------------------------------------------
  REAL(KIND=dp) :: SqrtElementMetric, U, V, W, S, &
       Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), Normal(3), r, ReP, ImP, &
       u1, v1, u2, v2, u3, v3, un, vn

  LOGICAL :: Stat
  INTEGER :: i, j, p, q, t, DIM
  TYPE(GaussIntegrationPoints_t) :: IntegStuff
!------------------------------------------------------------------------------
!  Impedance1 = CMPLX(0.0d0,0.0d0, kind=dp)
!  Impedance2 = CMPLX(0.0d0,0.0d0, kind=dp)
!------------------------------------------------------------------------------
! Numerical integration
!------------------------------------------------------------------------------
  IntegStuff = GaussPoints( Element )
!------------------------------------------------------------------------------
  DO t=1, IntegStuff % n
     U = IntegStuff % u(t)
     V = IntegStuff % v(t)
     W = IntegStuff % w(t)
     S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!    Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
     stat = ElementInfo( Element, Nodes, U, V, W, SqrtElementMetric, &
          Basis, dBasisdx, ddBasisddx, .FALSE., .FALSE. )
     s = s * SqrtElementMetric
     IF (CoordSys == AxisSymmetric) THEN
        r = SUM( Basis * Nodes % x(1:n) )
        s = r * s
     END IF

     Normal = Normalvector(Element, Nodes, U, V, .TRUE.)

     ReP = SUM(Pressure(1,1:n) * Basis(1:n))
     ImP = SUM(Pressure(2,1:n) * Basis(1:n))

     u1 = SUM(Velo(1,1:n) * Basis(1:n))
     v1 = SUM(Velo(2,1:n) * Basis(1:n))

     u2 = SUM(Velo(3,1:n) * Basis(1:n))
     v2 = SUM(Velo(4,1:n) * Basis(1:n))

     u3 = SUM(Velo(5,1:n) * Basis(1:n))
     v3 = SUM(Velo(6,1:n) * Basis(1:n))

     un = u1 * Normal(1) + u2 * Normal(2) + u3 * Normal(3)
     vn = v1 * Normal(1) + v2 * Normal(2) + v3 * Normal(3)

     Impedance1 = Impedance1 + CMPLX(ReP, ImP, kind=dp) * CMPLX(ReP, ImP, kind=dp) * s
     Impedance2 = Impedance2 - CMPLX(ReP, ImP, kind=dp) * CMPLX(un, vn, kind=dp) * s

     C3 = C3 - CMPLX(ReP, ImP, kind=dp) * CMPLX(un, vn, kind=dp) * s
     C4 = C4 + CMPLX(un, vn, kind=dp) * CMPLX(un, vn, kind=dp) * s


     Area = Area + s
  END DO

!-----------------------------------------------------------------------------
END SUBROUTINE SurfaceImpedanceIntegration
!-----------------------------------------------------------------------------



!------------------------------------------------------------------------------
FUNCTION DiscontIndexes( Solver, Element, n ) RESULT(GapIndexes)
!------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    TYPE(Element_t), TARGET :: Element
    INTEGER :: n
    INTEGER :: GapIndexes(n)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Parent, Left, Right
    INTEGER :: i,j,k,mini
    REAL(KIND=dp) :: x0,y0,z0,dx,dy,dz,dist,mindist
!------------------------------------------------------------------------------
    Left  => Element % BoundaryInfo % Left
    Right => Element % BoundaryInfo % Right

    GapIndexes(1:n) = 0
    IF ( .NOT. ASSOCIATED(Left) .OR. .NOT. ASSOCIATED(Right) ) RETURN

    DO i=1,n
      Parent => Left
      k = Element % NodeIndexes(i)
      IF ( ANY( Parent % NodeIndexes == k ) ) Parent => Right

      x0 = Solver % Mesh % Nodes % x(k)
      y0 = Solver % Mesh % Nodes % y(k)
      z0 = Solver % Mesh % Nodes % z(k)

      mini = 0
      mindist = HUGE(mindist)
      DO j=1,Parent % TYPE % NumberOfNodes
        k = Parent % NodeIndexes(j)
        dx = Solver % Mesh % Nodes % x(k) - x0
        dy = Solver % Mesh % Nodes % y(k) - y0
        dz = Solver % Mesh % Nodes % z(k) - z0
        dist = dx**2 + dy**2 + dz**2
        IF( dist < mindist) THEN
          mini = k
          mindist = dist
        END IF
      END DO

      IF(mini /= Element % NodeIndexes(i)) GapIndexes(i) = mini
    END DO
!------------------------------------------------------------------------------
  END FUNCTION DiscontIndexes
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE SetDirichletPoints( StiffMatrix, ForceVector, DOF, NDOFs, &
      Perm, n, NodeIndexes, NodeValues)
!------------------------------------------------------------------------------

    INTEGER :: n
    TYPE(Matrix_t), POINTER :: StiffMatrix
    REAL(KIND=dp) :: ForceVector(:), NodeValues(:)
    INTEGER :: DOF, NDOFs, Perm(:), NodeIndexes(:)
!------------------------------------------------------------------------------

    INTEGER :: i,PermIndex
    REAL(KIND=dp) :: s

!------------------------------------------------------------------------------

    DO i=1,n
      PermIndex = Perm(NodeIndexes(i))
      IF ( PermIndex == 0 ) CYCLE

      PermIndex = NDOFs * (PermIndex-1) + DOF

      IF ( StiffMatrix % FORMAT == MATRIX_SBAND ) THEN

        CALL SBand_SetDirichlet( StiffMatrix,ForceVector,PermIndex,NodeValues(i) )

      ELSE IF ( StiffMatrix % FORMAT == MATRIX_CRS .AND. &
          StiffMatrix % Symmetric ) THEN

        CALL CRS_SetSymmDirichlet(StiffMatrix,ForceVector,PermIndex,NodeValues(i) )

      ELSE

        s = StiffMatrix % Values(StiffMatrix % Diag(PermIndex))
        ForceVector(PermIndex) = NodeValues(i) * s
        CALL ZeroRow( StiffMatrix,PermIndex )
        CALL SetMatrixElement( StiffMatrix,PermIndex,PermIndex,1.0d0*s )

      END IF
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE SetDirichletPoints
!------------------------------------------------------------------------------
SUBROUTINE AcousticShellInterface(Model, Solver, Flow, FlowPerm, GapIndexes, &
    VelocityComponents, Dofs, AngularFrequency, dim, StiffMatrix, ForceVector)

  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: Flow(:)
  INTEGER :: FlowPerm(:), GapIndexes(:), VelocityComponents, Dofs, dim
  REAL(KIND=dp) :: AngularFrequency
  TYPE(Matrix_t), POINTER :: StiffMatrix
  REAL(KIND=dp) :: ForceVector(:)
  INTEGER :: BC,t,n
  INTEGER :: i
  TYPE(Variable_t), POINTER :: Var
  REAL (KIND=dp), POINTER :: Deflection(:)
  INTEGER, POINTER :: DeflectionPerm(:)
  INTEGER :: Indexes(MAX_ELEMENT_NODES)
  REAL (KIND=dp) :: NodeValues(MAX_ELEMENT_NODES)
  INTEGER :: sign
  TYPE(Element_t), POINTER :: Element
  LOGICAL :: AcousticShell, AcousticNoslip, AcousticIsothermal, DoneDiscont, GotIt

  Var => VariableGet( Model % Variables, 'U' )
  IF(.NOT. ASSOCIATED(Var)) THEN
    CALL Warn('Acoustics','Could not find variable U for FSI simulation')
    RETURN
  END IF
  Deflection => Var % Values
  DeflectionPerm => Var % Perm

  DO BC=1,Model % NumberOfBCs

    AcousticShell = GetLogical( Model % BCs(BC) % Values,'Solid Interaction',GotIt)
    AcousticNoslip = GetLogical( Model % BCs(BC) % Values,'Acoustic Wall Noslip',GotIt)
    AcousticIsothermal = GetLogical( Model % BCs(BC) % Values,'Acoustic Wall Isothermal',GotIt)

    IF(.NOT. (AcousticShell .OR. AcousticNoslip .OR. AcousticIsothermal)) CYCLE

    DO t = Model % NumberOfBulkElements + 1, &
        Model % NumberOfBulkElements + Model % NumberOfBoundaryElements

      Element => Model % Elements(t)
      IF ( Element % BoundaryInfo % Constraint /= Model % BCs(BC) % Tag ) CYCLE

      Model % CurrentElement => Element
      n = Element % TYPE % NumberOfNodes
      Indexes(1:n) = Element % NodeIndexes

      IF( AcousticShell) THEN
        GapIndexes(1:n) = Indexes(1:n)
        !DoneDiscont = .FALSE.

        DO i=1,VelocityComponents
          ! Real part
          NodeValues(1:n) = -AngularFrequency * Deflection( 6 * (DeflectionPerm(Indexes(1:n)) - 1) + 2*i)
          CALL SetDirichletPoints( StiffMatrix, ForceVector, 2*i-1, 10, FlowPerm, n, GapIndexes, NodeValues)
          ! Im part
          NodeValues(1:n) = AngularFrequency * Deflection( 6 * (DeflectionPerm(Indexes(1:n)) - 1) + 2*i-1)
          CALL SetDirichletPoints( StiffMatrix, ForceVector, 2*i, 10, FlowPerm, n, GapIndexes, NodeValues)
        END DO

        IF(AcousticIsothermal) THEN
          NodeValues(1:n) = 0.0
          CALL SetDirichletPoints( StiffMatrix, ForceVector, Dofs-3, Dofs, FlowPerm, n, GapIndexes, NodeValues)
          CALL SetDirichletPoints( StiffMatrix, ForceVector, Dofs-2, Dofs, FlowPerm, n, GapIndexes, NodeValues)
        END IF

        !IF( .NOT. DoneDiscont ) THEN
        !  DoneDiscont = .FALSE.
        !  GapIndexes(1:n) = DiscontIndexes(Solver, Element, n)
        !  IF( ALL(GapIndexes(1:n) /= 0)) THEN
        !    IF(.FALSE.) PRINT *,'found discontinuity',n,GapIndexes(1:n)
        !    GOTO 100
        !  END IF
        !END IF

      ELSE IF (AcousticNoslip) THEN
        NodeValues(1:n) = 0.0

        DO i=1, VelocityComponents
          CALL SetDirichletPoints( StiffMatrix, ForceVector, 2*i-1, Dofs, FlowPerm, n, Indexes, NodeValues)
          CALL SetDirichletPoints( StiffMatrix, ForceVector, 2*i, Dofs, FlowPerm, n, Indexes, NodeValues)
        END DO
        IF(AcousticIsothermal) THEN
          CALL SetDirichletPoints( StiffMatrix, ForceVector, Dofs-3, Dofs, FlowPerm, n, Indexes, NodeValues)
          CALL SetDirichletPoints( StiffMatrix, ForceVector, Dofs-2, Dofs, FlowPerm, n, Indexes, NodeValues)
        END IF
      END IF
    END DO
  END DO

END SUBROUTINE AcousticShellInterface



!---------------------------------------------------------------------------------------
SUBROUTINE FSIIntegration(Element, Velo, Pressure, Nodes, n, uf, wf, Area, CoordSys, AngularFrequency)
!---------------------------------------------------------------------------------------
  TYPE(Element_t), TARGET :: Element
  INTEGER :: CoordSys
  REAL(kind=dp) :: AngularFrequency
  REAL(kind=dp) :: Velo(:,:), Pressure(:,:), Area
  TYPE(Nodes_t) :: Nodes
  INTEGER :: n
  COMPLEX(kind=dp) :: uf, wf
!------------------------------------------------------------------------------
  REAL(KIND=dp) :: SqrtElementMetric, U, V, W, S, &
       Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), Normal(3), r, ReP, ImP, &
       u1, v1, u2, v2, u3, v3, un, vn

  LOGICAL :: Stat
  INTEGER :: i, j, p, q, t, DIM
  TYPE(GaussIntegrationPoints_t) :: IntegStuff
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Numerical integration
!------------------------------------------------------------------------------
  IntegStuff = GaussPoints( Element )
!------------------------------------------------------------------------------
  DO t=1, IntegStuff % n
     U = IntegStuff % u(t)
     V = IntegStuff % v(t)
     W = IntegStuff % w(t)
     S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!    Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
     stat = ElementInfo( Element, Nodes, U, V, W, SqrtElementMetric, &
          Basis, dBasisdx, ddBasisddx, .FALSE., .FALSE. )
     s = s * SqrtElementMetric
     IF (CoordSys == AxisSymmetric) THEN
        r = SUM( Basis * Nodes % x(1:n) )
        s = r * s
     END IF

     Normal = Normalvector(Element, Nodes, U, V, .TRUE.)

     ReP = SUM(Pressure(1,1:n) * Basis(1:n))
     ImP = SUM(Pressure(2,1:n) * Basis(1:n))

     u1 = SUM(Velo(1,1:n) * Basis(1:n))
     v1 = SUM(Velo(2,1:n) * Basis(1:n))

     u2 = SUM(Velo(3,1:n) * Basis(1:n))
     v2 = SUM(Velo(4,1:n) * Basis(1:n))

     u3 = SUM(Velo(5,1:n) * Basis(1:n))
     v3 = SUM(Velo(6,1:n) * Basis(1:n))

     un = u1 * Normal(1) + u2 * Normal(2) + u3 * Normal(3)
     vn = v1 * Normal(1) + v2 * Normal(2) + v3 * Normal(3)


     ! U = \int (u.u*)
     uf = uf + (s / AngularFrequency ** 2) * ( CMPLX(u1, v1, kind=dp) * CMPLX(u1, -v1, kind=dp) + &
          CMPLX(u2, v2, kind=dp) * CMPLX(u2, -v2, kind=dp) + &
          CMPLX(u3, v3, kind=dp) * CMPLX(u3, -v3, kind=dp) )

     ! Kf = \int p.u*
     wf = wf + CMPLX(0.0d0, 1.0d0, kind=dp) * (s / AngularFrequency) * CMPLX(ReP, ImP, kind=dp) * &
          CMPLX(un, -vn, kind=dp)

     Area = Area + s
  END DO

!-----------------------------------------------------------------------------
END SUBROUTINE FSIIntegration
!-----------------------------------------------------------------------------
  ! -- BEM-coupling bodies (verbatim driver regions; Fortran-owned SAVE) --------
  subroutine aca_bem_pre(firstvisit) bind(c, name="aca_bem_pre")
    integer(c_int), value :: firstvisit
    character(len=MAX_NAME_LEN) :: ElmerStatusFile, CompaStatusFile, CompaMeshPrefix
    character(30) :: str2
    logical :: GotIt
    integer :: i, istat
    ElmerStatusFile = ListGetString(acc_model % simulation, 'Elmer Status File', GotIt)
    IF ( .NOT. GotIt ) ElmerStatusFile = 'ELMER_STATUS'
    CompaStatusFile = ListGetString(acc_model % simulation, 'Compa Status File', GotIt)
    IF ( .NOT. GotIt ) CompaStatusFile = 'COMPA_STATUS'
    CompaMeshPrefix = ListGetString(acc_model % simulation, 'Compa Mesh File Prefix', GotIt)
    IF ( .NOT. GotIt ) CompaMeshPrefix = 'mesh'
    IF (firstvisit == 0) THEN
      OPEN( 10, FILE = CompaStatusFile, status='OLD', IOSTAT = istat )
      IF ( istat /= 0) &
          CALL Fatal( 'AcousticsSolver', 'Cannot open Compa Status File' )
      READ(10,'(A)') str2
      CLOSE(10)
      i = INDEX(str2,'$') + 1
      IF ( str2(i:i+3) == 'STOP' .OR. str2(i:i+4) == 'PANIC' ) &
          CALL Fatal( 'AcousticsSolver', 'CompaSolver is not willing to continue' )
    END IF
    OPEN( 10, FILE = ElmerStatusFile, status='REPLACE', IOSTAT = istat )
    IF ( istat /= 0) THEN
      CALL Fatal( 'AcousticsSolver', 'Cannot open Elmer Status File' )
    ELSE
      WRITE( 10, '(A)', ADVANCE='NO') '$WAIT'
    END IF
    CLOSE(10)
  end subroutine aca_bem_pre
  subroutine aca_bem_nodes() bind(c, name="aca_bem_nodes")
    character(len=MAX_NAME_LEN) :: DataFile, CompaMeshPrefix
    logical :: GotIt
    integer :: i, j, k, m, t, istat
    real(kind=dp) :: ReP, ImP
    CompaMeshPrefix = ListGetString(acc_model % simulation, 'Compa Mesh File Prefix', GotIt)
    IF ( .NOT. GotIt ) CompaMeshPrefix = 'mesh'
    WRITE(DataFile,'(A,A)') TRIM(CompaMeshPrefix), '.P'
    OPEN( 10, FILE = DataFile, status='OLD', IOSTAT = istat )
    IF ( istat /= 0) THEN
      WRITE( acc_message, * ) 'Cannot open file ', DataFile
      CALL Fatal( 'AcousticsSolver', acc_message )
    END IF
    j = 0
    DO
      READ( 10, *, IOSTAT = istat) i
      IF (istat == 0) THEN
        j = j + 1
      ELSE
        EXIT
      END IF
    END DO
    IF (.NOT. acc_bem_created) THEN
      ALLOCATE( acc_ainodes(j), acc_aires(j,2), STAT=istat )
      IF ( istat /= 0 ) CALL Fatal( 'AcousticsSolver', 'Memory allocation error.' )
      acc_ainodes(1:j) = 0
    END IF
    REWIND 10
    OPEN( 20, FILE = 'mesh.nodes', status='OLD', IOSTAT = istat )
    IF ( istat /= 0) &
        CALL Fatal( 'AcousticsSolver', 'Cannot open file mesh.nodes' )
    DO i = 1, j
      READ( 10, *) t, ReP, ImP
      k = 0
      DO
        READ( 20, *, IOSTAT = istat) m
        IF (istat == 0) THEN
          k = k + 1
          IF ( m==t ) THEN
            acc_ainodes(i) = k
            acc_aires(i,1) = ReP
            acc_aires(i,2) = -1.0d0*ImP
            EXIT
          END IF
        ELSE
          WRITE( acc_message, * ) 'Inconsistent node numbering in the file ', DataFile
          CALL Fatal( 'AcousticsSolver', acc_message )
        END IF
      END DO
      REWIND 20
    END DO
    CLOSE(20)
    CLOSE(10)
    acc_bem_created = .TRUE.
  end subroutine aca_bem_nodes
  subroutine aca_bem_write(firstvisit, angularfreq, cdim, WallVelocity, density, Load, wsz) &
      bind(c, name="aca_bem_write")
    integer(c_int), value :: firstvisit, cdim, wsz
    real(c_double), value :: angularfreq
    real(c_double) :: WallVelocity(6,wsz), density(wsz), Load(6,wsz)
    character(len=MAX_NAME_LEN) :: DataFile, CompaMeshPrefix
    character(len=MAX_NAME_LEN) :: EquationName
    type(element_t), pointer :: CurrentElement, Parent
    type(valuelist_t), pointer :: Material
    integer, pointer :: NodeIndexes(:)
    type(nodes_t) :: ElementNodes
    logical :: GotIt, stat, BEMBoundary, AcousticInterface
    integer :: i, j, k, n, t, istat
    complex(kind=dp) :: AverVel
    logical :: firstv
    firstv = (firstvisit /= 0)
    EquationName = acc_eqname
    CALL Info('AcousticsSolver', 'Writing input file for BEM solver')
    IF (firstv) THEN
      ALLOCATE( acc_bem_elem(acc_solver % mesh % numberofboundaryelements), STAT=istat )
      IF ( istat /= 0 ) CALL Fatal( 'AcousticsSolver', 'Memory allocation error.' )
      OPEN( 10, FILE = 'mesh.boundary', status='OLD')
      DO t = 1, acc_solver % mesh % numberofboundaryelements
        READ( 10, *) j
        acc_bem_elem(t) = j
      END DO
      CLOSE(10)
    END IF
    CompaMeshPrefix = ListGetString(acc_model % simulation, 'Compa Mesh File Prefix', GotIt)
    IF ( .NOT. GotIt ) CompaMeshPrefix = 'mesh'
    WRITE(DataFile,'(A,A,A)') TRIM(CompaMeshPrefix), '.', 'param'
    OPEN( 10, FILE = DataFile, status='REPLACE')
    DO t = acc_solver % mesh % numberofbulkelements + 1, &
        acc_solver % mesh % numberofbulkelements + acc_solver % mesh % numberofboundaryelements
      CurrentElement => acc_solver % mesh % elements(t)
      acc_model % currentelement => CurrentElement
      CurrentModel % currentelement => CurrentElement
      n = CurrentElement % TYPE % NumberOfNodes
      NodeIndexes => CurrentElement % NodeIndexes
      DO i=1,acc_model % numberofbcs
        IF ( CurrentElement % BoundaryInfo % Constraint == &
            acc_model % BCs(i) % Tag ) THEN
          BEMBoundary = ListGetLogical( acc_model % BCs(i) % Values, &
              'BEM Boundary', GotIt )
          AcousticInterface = ListGetLogical( acc_model % BCs(i) % Values, &
              'Acoustic Interface', GotIt )
          IF ( BEMBoundary .OR. AcousticInterface) THEN
            IF (AcousticInterface) THEN
              ElementNodes % x(1:n) = acc_solver % mesh % nodes % x(NodeIndexes)
              ElementNodes % y(1:n) = acc_solver % mesh % nodes % y(NodeIndexes)
              ElementNodes % z(1:n) = acc_solver % mesh % nodes % z(NodeIndexes)
              Parent => CurrentELement % BoundaryInfo % Left
              stat = ASSOCIATED( Parent )
              IF (stat) stat = ALL(acc_solver % variable % perm(Parent % NodeIndexes(1:n)) > 0)
              IF ( .NOT. stat) THEN
                Parent => CurrentELement % BoundaryInfo % Right
                stat = ASSOCIATED( Parent )
                IF (stat) stat = ALL(acc_solver % variable % perm(Parent % NodeIndexes(1:n)) > 0)
                IF ( .NOT. stat )  CALL Fatal( 'AcousticsSolver', &
                    'No parent element can be found for given boundary element' )
              END IF
              IF ( .NOT. CheckElementEquation( acc_model, Parent, EquationName ) ) CYCLE
              k = ListGetInteger( acc_model % Bodies(Parent % Bodyid) % Values, &
                  'Material' )
              Material => acc_model % Materials(k) % Values
              Density(1:n) = ListGetReal( Material, 'Equilibrium Density', &
                  n, NodeIndexes )
              DO j=1,n
                k = acc_solver % variable % perm( NodeIndexes(j) )
                WallVelocity(1,j) =  acc_solver % variable % values( (k-1)*(cdim*2+4)+1 )
                WallVelocity(2,j) =  acc_solver % variable % values( (k-1)*(cdim*2+4)+2 )
                WallVelocity(3,j) =  acc_solver % variable % values( (k-1)*(cdim*2+4)+3 )
                WallVelocity(4,j) =  acc_solver % variable % values( (k-1)*(cdim*2+4)+4 )
                IF (cdim > 2) THEN
                  WallVelocity(5,j) =  acc_solver % variable % values( (k-1)*(cdim*2+4)+5 )
                  WallVelocity(6,j) =  acc_solver % variable % values( (k-1)*(cdim*2+4)+6 )
                ELSE
                  WallVelocity(5,j) =  0.0d0
                  WallVelocity(6,j) =  0.0d0
                END IF
              END DO
              CALL ComputeAverageVelocity( CurrentElement, n, WallVelocity, AverVel, acc_solver)
              AverVel = AverVel * CMPLX(0.0d0, -1.0d0 * angularfreq * Density(1), kind=dp )
              AverVel = CONJG(AverVel)
              WRITE( 10, '(I9,6e23.15)',ADVANCE='NO') acc_bem_elem( CurrentElement % ElementIndex - &
                  acc_solver % mesh % NumberOfBulkElements ), 0.0d+0, 0.0d+0, 1.0d+0, 0.0d0, &
                  REAL(AverVel), AIMAG(AverVel)
              WRITE( 10,* ) ''
            ELSE
              Load(1,1:n) = ListGetReal( acc_model % BCs(i) % Values, &
                  'Re a', n, NodeIndexes, GotIt )
              Load(2,1:n) = ListGetReal( acc_model % BCs(i) % Values, &
                  'Im a', n, NodeIndexes, GotIt )
              Load(3,1:n) = ListGetReal( acc_model % BCs(i) % Values, &
                  'Re b', n, NodeIndexes, GotIt )
              Load(4,1:n) = ListGetReal( acc_model % BCs(i) % Values, &
                  'Im b', n, NodeIndexes, GotIt )
              Load(5,1:n) = ListGetReal( acc_model % BCs(i) % Values, &
                  'Re c', n, NodeIndexes, GotIt )
              Load(6,1:n) = ListGetReal( acc_model % BCs(i) % Values, &
                  'Im c', n, NodeIndexes, GotIt )
              WRITE( 10, '(I9,6e23.15)',ADVANCE='NO') acc_bem_elem( CurrentElement % ElementIndex - &
                  acc_solver % mesh % NumberOfBulkElements ), Load(1,1), Load(2,1), Load(3,1), Load(4,1), &
                  Load(5,1), Load(6,1)
              WRITE( 10,* ) ''
            END IF
          END IF
        END IF
      END DO
    END DO
    CLOSE(10)
  end subroutine aca_bem_write
  subroutine aca_bem_post(relchange, calls) bind(c, name="aca_bem_post")
    real(c_double), value :: relchange
    integer(c_int), value :: calls
    character(len=MAX_NAME_LEN) :: ElmerStatusFile
    character(3) :: str
    logical :: GotIt
    integer :: istat, MaxCoupledIterations
    real(kind=dp) :: CoupledTolerance
    ElmerStatusFile = ListGetString(acc_model % simulation, 'Elmer Status File', GotIt)
    IF ( .NOT. GotIt ) ElmerStatusFile = 'ELMER_STATUS'
    OPEN( 10, FILE = ElmerStatusFile, status='REPLACE', IOSTAT = istat )
    IF ( istat /= 0) THEN
      CALL Fatal( 'AcousticsSolver', 'Cannot open Elmer Status File' )
    ELSE
      CoupledTolerance = ListGetConstReal( acc_solver % values, &
          'Steady State Convergence Tolerance' )
      MaxCoupledIterations =  ListGetInteger( acc_model % simulation, &
          'Steady State Max Iterations' )
      IF ( (relchange < CoupledTolerance) .OR. (calls >= MaxCoupledIterations) ) THEN
        WRITE( 10, '(A)', ADVANCE='NO') '$STOP'
      ELSE
        WRITE(str(1:3),'(I3)') calls
        WRITE( 10, '(A)', ADVANCE='NO') '$CONTINUE #' // ADJUSTL(str(1:3))
      END IF
    END IF
    CLOSE(10)
  end subroutine aca_bem_post

  subroutine aca_info_bem_writing() bind(c, name="aca_info_bem_writing")
    call Info('AcousticsSolver', 'Writing input file for BEM solver')
  end subroutine aca_info_bem_writing
  subroutine aca_fluidic_all(cdim, visc, lam, bulkv, wallvelo, lod, &
      totforce, totmoment, totarea, momout, nn, nbc) &
      bind(c, name="aca_fluidic_all")
    integer(c_int), value :: cdim, nn, nbc
    real(c_double) :: visc(nn), lam(nn), bulkv(nn), wallvelo(6,nn), lod(6,nn)
    real(c_double) :: totforce(nbc,6), totmoment(nbc,6), totarea(nbc)
    integer(c_int) :: momout(nbc)
    type(element_t), pointer :: CurrentElement, Parent
    type(valuelist_t), pointer :: Material
    integer, pointer :: NodeIndexes(:)
    type(nodes_t) :: ElementNodes, ParentNodes
    logical :: GotIt, stat, CalculateMoment
    integer :: i, j, k, n, t, np
    real(kind=dp) :: MomentAbout(3), Traction(6), Moment(6), Area
    type(valuelist_t), pointer :: bcv
    character(len=MAX_NAME_LEN) :: BoundaryName
    integer :: nlen
    totforce = 0.0d0
    totmoment = 0.0d0
    totarea = 0.0d0
    momout = 0
    DO t = acc_solver % mesh % numberofbulkelements + 1, &
        acc_solver % mesh % numberofbulkelements + acc_solver % mesh % numberofboundaryelements
      CurrentElement => acc_solver % mesh % elements(t)
      acc_model % currentelement => CurrentElement
      CurrentModel % currentelement => CurrentElement
      IF ( .NOT. ActiveBoundaryElement(CurrentElement, CurrentModel % solver ) ) CYCLE
      n = CurrentElement % TYPE % NumberOfNodes
      NodeIndexes => CurrentElement % NodeIndexes
      IF (ANY(acc_solver % variable % perm(NodeIndexes(1:n)) == 0)) CYCLE
      DO i=1,acc_model % numberofbcs
        IF ( CurrentElement % BoundaryInfo % Constraint == &
            acc_model % BCs(i) % Tag ) THEN
          bcv => acc_model % bcs(i) % values
          IF ( .NOT. ListGetLogical( bcv, 'Calculate Fluidic Force', GotIt ) ) CYCLE
          ElementNodes % x(1:n) = acc_solver % mesh % nodes % x(NodeIndexes)
          ElementNodes % y(1:n) = acc_solver % mesh % nodes % y(NodeIndexes)
          ElementNodes % z(1:n) = acc_solver % mesh % nodes % z(NodeIndexes)
          Parent => CurrentELement % BoundaryInfo % Left
          stat = ASSOCIATED( Parent )
          IF (stat) stat = ALL( acc_solver % variable % perm( &
              Parent % NodeIndexes(1: Parent % TYPE % NumberOfNodes)) > 0 )
          IF ( .NOT. stat) THEN
            Parent => CurrentELement % BoundaryInfo % Right
            stat = ASSOCIATED( Parent )
            IF (stat) stat = ALL( acc_solver % variable % perm( &
                Parent % NodeIndexes(1: Parent % TYPE % NumberOfNodes)) > 0 )
            IF ( .NOT. stat )  CALL Fatal( 'AcousticsSolver', &
                'No parent element can be found for given boundary element' )
          END IF
          np = Parent % TYPE % NumberOfNodes
          ParentNodes % x(1:np) = acc_solver % mesh % nodes % x(Parent % NodeIndexes)
          ParentNodes % y(1:np) = acc_solver % mesh % nodes % y(Parent % NodeIndexes)
          ParentNodes % z(1:np) = acc_solver % mesh % nodes % z(Parent % NodeIndexes)
          k = ListGetInteger( acc_model % Bodies(Parent % Bodyid) % Values, &
              'Material' )
          Material => acc_model % Materials(k) % Values
          visc(1:np) = ListGetReal( Material, 'Viscosity', &
              np, Parent % NodeIndexes )
          lam(1:np) = -2.0d0/3.0d0 * visc(1:np)
          bulkv(1:np) = ListGetReal( Material, 'Bulk Viscosity', &
              np, Parent % NodeIndexes, GotIt )
          IF (GotIt) lam(1:np) = bulkv(1:np) - 2.0d0/3.0d0 * visc(1:np)
          MomentAbout(1) = ListGetConstReal( bcv, 'Moment About 1', GotIt )
          MomentAbout(2) = ListGetConstReal( bcv, 'Moment About 2', CalculateMoment )
          CalculateMoment = GotIt .OR. CalculateMoment
          MomentAbout(3) = ListGetConstReal( bcv, 'Moment About 3', GotIt)
          CalculateMoment = GotIt .OR. CalculateMoment
          IF (CalculateMoment) THEN
            momout(i) = 1
          ELSE
            momout(i) = 0
          END IF
          DO j=1,np
            k = acc_solver % variable % perm( Parent % NodeIndexes(j) )
            wallvelo(1,j) = acc_solver % variable % values( (k-1)*(cdim*2+4)+1 )
            wallvelo(2,j) = acc_solver % variable % values( (k-1)*(cdim*2+4)+2 )
            wallvelo(3,j) = acc_solver % variable % values( (k-1)*(cdim*2+4)+3 )
            wallvelo(4,j) = acc_solver % variable % values( (k-1)*(cdim*2+4)+4 )
            IF (cdim > 2) THEN
              wallvelo(5,j) = acc_solver % variable % values( (k-1)*(cdim*2+4)+5 )
              wallvelo(6,j) = acc_solver % variable % values( (k-1)*(cdim*2+4)+6 )
              lod(1,j) = acc_solver % variable % values( (k-1)*(cdim*2+4)+9 )
              lod(2,j) = acc_solver % variable % values( (k-1)*(cdim*2+4)+10 )
            ELSE
              wallvelo(5,j) = 0.0d0
              wallvelo(6,j) = 0.0d0
              lod(1,j) =  acc_solver % variable % values( (k-1)*(cdim*2+4)+7 )
              lod(2,j) =  acc_solver % variable % values( (k-1)*(cdim*2+4)+8 )
            END IF
          END DO
          Traction = 0.0d0
          Moment = 0.0d0
          CALL SurfaceForceIntegration(CurrentElement, Parent, Traction, Moment, &
              MomentAbout, Area, CalculateMoment, wallvelo, lod, visc, &
              bulkv, ElementNodes, ParentNodes, np, acc_solver, CurrentCoordinateSystem())
          totforce(i,1:6) = totforce(i,1:6) + Traction(1:6)
          totarea(i) = totarea(i) + Area
          IF (CalculateMoment) THEN
            totmoment(i,1:6) = totmoment(i,1:6) + Moment(1:6)
          END IF
        END IF
      END DO
    END DO
    DO k=1, acc_model % numberofbcs
      bcv => acc_model % bcs(k) % values
      IF ( .NOT. ListGetLogical(bcv, 'Calculate Fluidic Force', stat ) ) CYCLE
      IF(acc_model % numberofbcs < 10) THEN
        WRITE( BoundaryName, '("bc ",I1)') k
      ELSE
        WRITE( BoundaryName, '("bc ",I2)') k
      END IF
      nlen = LEN_TRIM(BoundaryName)
      CALL Info('ForceCompute','Forces on Boundary '//BoundaryName(1:nlen),Level=4 )
      WRITE( acc_message, '("Real Fluidic Force (X,Y,Z):", 3ES17.6E2)') totforce(k,1:3)
      CALL Info( 'ForceCompute', acc_message, Level=4 )
      WRITE( acc_message, '("Imaginary Fluidic Force (X,Y,Z):", 3ES17.6E2)') totforce(k,4:6)
      CALL Info( 'ForceCompute', acc_message, Level=4 )
      WRITE( acc_message, '("Contact Area:   ", ES17.6E2)') totarea(k)
      CALL Info( 'ForceCompute', acc_message, Level=4 )
      CALL ListAddConstReal( acc_model % simulation, &
          'res: re contact force 1 '//BoundaryName(1:nlen), totforce(k,1) )
      CALL ListAddConstReal( acc_model % simulation, &
          'res: re contact force 2 '//BoundaryName(1:nlen), totforce(k,2) )
      IF ( cdim > 2 )  CALL ListAddConstReal( acc_model % simulation, &
          'res: re contact force 3 '//BoundaryName(1:nlen), totforce(k,3) )
      CALL ListAddConstReal( acc_model % simulation, &
          'res: im contact force 1 '//BoundaryName(1:nlen), totforce(k,4) )
      CALL ListAddConstReal( acc_model % simulation, &
          'res: im contact force 2 '//BoundaryName(1:nlen), totforce(k,5) )
      IF ( cdim > 2 )  CALL ListAddConstReal( acc_model % simulation, &
          'res: im contact force 3 '//BoundaryName(1:nlen), totforce(k,6) )
      IF ( momout(k) /= 0 ) THEN
        WRITE( acc_message, &
            '("Real Moment about (",ES9.3E1,",",ES9.3E1,",",ES9.3E1,"):",3Es14.6E2)') &
            MomentAbout(1:3), totmoment(k,1:3)
        CALL Info( 'ForceCompute', acc_message, Level=4 )
        WRITE( acc_message, &
            '("Imaginary Moment about (",ES9.3E1,",",ES9.3E1,",",ES9.3E1,"):",3Es14.6E2)') &
            MomentAbout(1:3), totmoment(k,4:6)
        CALL Info( 'ForceCompute', acc_message, Level=4 )
        CALL ListAddConstReal( acc_model % simulation, &
            'res: re contact moment 1 '//BoundaryName(1:nlen), totmoment(k,1) )
        CALL ListAddConstReal( acc_model % simulation, &
            'res: re contact moment 2 '//BoundaryName(1:nlen), totmoment(k,2) )
        CALL ListAddConstReal( acc_model % simulation, &
            'res: re contact moment 3 '//BoundaryName(1:nlen), totmoment(k,3) )
        CALL ListAddConstReal( acc_model % simulation, &
            'res: im contact moment 1 '//BoundaryName(1:nlen), totmoment(k,4) )
        CALL ListAddConstReal( acc_model % simulation, &
            'res: im contact moment 2 '//BoundaryName(1:nlen), totmoment(k,5) )
        CALL ListAddConstReal( acc_model % simulation, &
            'res: im contact moment 3 '//BoundaryName(1:nlen), totmoment(k,6) )
      END IF
      CALL ListAddConstReal( acc_model % simulation, &
          'res: contact force area '//BoundaryName(1:nlen), totarea(k) )
    END DO
  end subroutine aca_fluidic_all
  subroutine aca_impedance_pass(cdim) bind(c, name="aca_impedance_pass")
    integer(c_int), value :: cdim
    integer, allocatable :: Bndries(:)
    integer :: AcousticI, i, j
    real(kind=dp), allocatable :: AcImpedances(:,:)
    type(nodes_t) :: kEN
    logical :: GotIt
    ALLOCATE( Bndries( acc_model % numberofbcs ) )
    Bndries = 0
    AcousticI = 0
    DO i = 1, acc_model % numberofbcs
      IF ( ListGetLogical( acc_model % BCs(i) % values, &
          'Calculate Acoustic Impedance', GotIt ) ) THEN
        AcousticI = 1
        Bndries(1) = i
        EXIT
      END IF
    END DO
    IF ( AcousticI > 0 ) THEN
      j = 1
      DO i = 1, acc_model % numberofbcs
        IF ( ListGetLogical( acc_model % BCs(i) % values, &
            'Impedance Target Boundary', GotIt ) ) THEN
          AcousticI = AcousticI + 1
          j = j + 1
          Bndries(j) = i
        END IF
      END DO
      ALLOCATE( AcImpedances( AcousticI, 2 ) )
      CALL ComputeAcousticImpedance( AcousticI, Bndries, AcImpedances, &
          acc_model, acc_solver, acc_solver % variable % values, &
          acc_solver % variable % perm, acc_solver % variable % dofs, cdim, kEN )
      WRITE( acc_message, * ) 'Self specific acoustic impedance on bc ', &
          Bndries(1)
      CALL INFO( 'Acoustics', acc_message, LEVEL=5 )
      WRITE( acc_message, * ) '  In-phase with velocity:     ', &
          AcImpedances(1,1)
      CALL INFO( 'Acoustics', acc_message, LEVEL=5 )
      WRITE( acc_message, * ) '  Out-of-phase with velocity: ', &
          AcImpedances(1,2)
      CALL INFO( 'Acoustics', acc_message, LEVEL=5 )
      DO i = 2, AcousticI
        WRITE( acc_message, * ) 'Cross specific acoustic impedance, bcs ', &
            Bndries(i), ' , ', Bndries(1)
        CALL INFO( 'Acoustics', acc_message, LEVEL=5 )
        WRITE( acc_message, * ) '  In-phase with velocity:     ', &
            AcImpedances(i,1)
        CALL INFO( 'Acoustics', acc_message, LEVEL=5 )
        WRITE( acc_message, * ) '  Out-of-phase with velocity: ', &
            AcImpedances(i,2)
        CALL INFO( 'Acoustics', acc_message, LEVEL=5 )
      END DO
      DO i = 2, AcousticI
        WRITE( acc_message, '(A,I1,A,I1)' ) &
            'res: out-of-phase cross acoustic impedance ', &
            Bndries(AcousticI - i + 2), ' , ', Bndries(1)
        CALL ListAddConstReal( acc_model % simulation, acc_message, &
            AcImpedances(AcousticI - i + 2,2) )
        WRITE( acc_message, '(A,I1,A,I1)' ) &
            'res: in-phase cross acoustic impedance ', &
            Bndries(AcousticI - i + 2), ' , ', Bndries(1)
        CALL ListAddConstReal( acc_model % simulation, acc_message, &
            AcImpedances(AcousticI - i + 2,1) )
      END DO
      WRITE( acc_message, '(A,I1)' ) &
          'res: out-of-phase self acoustic impedance ', &
          Bndries(1)
      CALL ListAddConstReal( acc_model % simulation, acc_message, &
          AcImpedances(1,2) )
      WRITE( acc_message, '(A,I1)' ) &
          'res: in-phase self acoustic impedance ', &
          Bndries(1)
      CALL ListAddConstReal( acc_model % simulation, acc_message, &
          AcImpedances(1,1) )
      DEALLOCATE( AcImpedances )
    END IF
    DEALLOCATE( Bndries )
  end subroutine aca_impedance_pass
  ! -- BlockPreconditioning coarse bodies (verbatim driver regions) ---------------
  subroutine aca_block_create() bind(c, name="aca_block_create")
    logical :: GotIt, OptimizeBW, VelocityAssembly
    integer :: dim, n, istat
    dim = CoordinateSystemDimension()
    OptimizeBW = ListGetLogical(acc_solver % values,'Optimize Bandwidth',gotIt)
    IF(.NOT. GotIt) OptimizeBW = .TRUE.
    VelocityAssembly = ListGetLogical(acc_solver % values, 'Velocity Assembly', GotIt)
    IF (.NOT. GotIt) VelocityAssembly = .FALSE.
    IF (VelocityAssembly) &
        acc_amatrix => CreateMatrix( acc_model, acc_solver, acc_solver % mesh, &
        acc_solver % variable % perm, 2*dim, &
        MATRIX_CRS, OptimizeBW, ListGetString( acc_solver % values, 'Equation' ) )
    acc_smatrix => CreateMatrix( acc_model, acc_solver, acc_solver % mesh, &
        acc_solver % variable % perm, &
        4, MATRIX_CRS, OptimizeBW, ListGetString( acc_solver % values, 'Equation' ) )
    acc_mmatrix => CreateMatrix( acc_model, acc_solver, acc_solver % mesh, &
        acc_solver % variable % perm, &
        2, MATRIX_CRS, OptimizeBW, ListGetString( acc_solver % values, 'Equation' ) )
    acc_lmatrix => CreateMatrix( acc_model, acc_solver, acc_solver % mesh, &
        acc_solver % variable % perm, &
        2, MATRIX_CRS, OptimizeBW, ListGetString( acc_solver % values, 'Equation' ) )
    n = acc_solver % mesh % maxelementdofs
    ALLOCATE( acc_slocal(4*n,4*n), acc_alocal(2*dim*n,2*dim*n), acc_anlocal(2*n,2*n), &
        STAT=istat )
    IF ( istat /= 0 ) THEN
      CALL Fatal( 'AcousticsSolver', 'Memory allocation error.' )
    END IF
  end subroutine aca_block_create

  subroutine aca_block_zero(velasm) bind(c, name="aca_block_zero")
    integer(c_int), value :: velasm
    IF (velasm /= 0) CALL CRS_ZeroMatrix( acc_amatrix )
    CALL CRS_ZeroMatrix( acc_smatrix )
    CALL CRS_ZeroMatrix( acc_mmatrix )
    CALL CRS_ZeroMatrix( acc_lmatrix )
  end subroutine aca_block_zero

  subroutine aca_block_bulk(visc, dens, sh, hr, temp, cond, lam, press, angfreq, &
      nn, n, cdim, coordsys, velasm) bind(c, name="aca_block_bulk")
    integer(c_int), value :: nn, n, cdim, coordsys, velasm
    real(c_double), value :: angfreq
    real(c_double) :: visc(nn), dens(nn), sh(nn), hr(nn), temp(nn), cond(nn), &
        lam(nn), press(nn)
    IF (velasm /= 0) THEN
      CALL CoupledVelocityMatrix( acc_alocal, visc, angfreq, dens, &
          acc_elem, n, cdim, coordsys)
      CALL UpdateGlobalPreconditioner( acc_amatrix, acc_alocal, n, 2*cdim, &
          acc_solver % variable % perm( acc_elem % nodeindexes ) )
    END IF
    CALL PressureLaplaceMatrix( acc_anlocal, visc, angfreq, dens, &
        acc_elem, n, cdim, coordsys)
    CALL UpdateGlobalPreconditioner( acc_lmatrix, acc_anlocal, n, 2, &
        acc_solver % variable % perm( acc_elem % nodeindexes ) )
    CALL PressureMassMatrix( acc_anlocal, acc_elem, n, cdim, coordsys)
    CALL UpdateGlobalPreconditioner( acc_mmatrix, acc_anlocal, n, 2, &
        acc_solver % variable % perm( acc_elem % nodeindexes ) )
    CALL SchurComplementMatrix( acc_slocal, angfreq, sh, hr, dens, press, &
        temp, cond, visc, lam, acc_elem, n, cdim, coordsys)
    CALL UpdateGlobalPreconditioner( acc_smatrix, acc_slocal, n, 4, &
        acc_solver % variable % perm( acc_elem % nodeindexes ) )
  end subroutine aca_block_bulk
  subroutine aca_block_bnd_imp(imp, sh, hr, dens, cond, angfreq, nn, n, cdim, coordsys, velasm) &
      bind(c, name="aca_block_bnd_imp")
    integer(c_int), value :: nn, n, cdim, coordsys, velasm
    real(c_double), value :: angfreq
    real(c_double) :: imp(4,nn), sh(nn), hr(nn), dens(nn), cond(nn)
    IF (velasm /= 0) THEN
      CALL VelocityImpedanceMatrix( acc_alocal, angfreq, dens, imp, &
          acc_elem, n, acc_nodes, cdim, coordsys)
      CALL UpdateGlobalPreconditioner( acc_amatrix, acc_alocal, n, 2*cdim, &
          acc_solver % variable % perm( acc_elem % nodeindexes ) )
    END IF
    CALL SchurComplementImpedanceMatrix( acc_slocal, imp, angfreq, sh, hr, &
        dens, cond, acc_elem, n, acc_nodes, cdim, coordsys)
    CALL UpdateGlobalPreconditioner( acc_smatrix, acc_slocal, n, 4, &
        acc_solver % variable % perm( acc_elem % nodeindexes ) )
  end subroutine aca_block_bnd_imp
  subroutine aca_block_bnd_slip(sh, hr, dens, temp, angfreq, cond, wallt, c1, c2, &
      nn, n, cdim, coordsys, velasm) bind(c, name="aca_block_bnd_slip")
    integer(c_int), value :: nn, n, cdim, coordsys, velasm
    real(c_double), value :: angfreq, c1, c2
    real(c_double) :: sh(nn), hr(nn), dens(nn), temp(nn), cond(nn), wallt(nn)
    IF (velasm /= 0) THEN
      CALL VelocitySlipMatrix( acc_alocal, sh, hr, dens, temp, angfreq, &
          wallt, c1, acc_elem, n, acc_nodes, cdim, coordsys)
      CALL UpdateGlobalPreconditioner( acc_amatrix, acc_alocal, n, 2*cdim, &
          acc_solver % variable % perm( acc_elem % nodeindexes ) )
    END IF
    CALL SchurComplementSlipMatrix( acc_slocal, sh, hr, dens, temp, &
        angfreq, cond, wallt, c2, acc_elem, n, acc_nodes, cdim, coordsys)
    CALL UpdateGlobalPreconditioner( acc_smatrix, acc_slocal, n, 4, &
        acc_solver % variable % perm( acc_elem % nodeindexes ) )
  end subroutine aca_block_bnd_slip
  subroutine aca_block_setbcs(cdim, velasm) bind(c, name="aca_block_setbcs")
    integer(c_int), value :: cdim, velasm
    IF (velasm /= 0) THEN
      CALL SetBoundaryConditions(acc_model, acc_amatrix, 'Re Velocity 1', 1, cdim*2, &
          acc_solver % variable % perm)
      CALL SetBoundaryConditions(acc_model, acc_amatrix, 'Im Velocity 1', 2, cdim*2,  &
          acc_solver % variable % perm)
      CALL SetBoundaryConditions(acc_model, acc_amatrix, 'Re Velocity 2', 3, cdim*2,  &
          acc_solver % variable % perm)
      CALL SetBoundaryConditions(acc_model, acc_amatrix, 'Im Velocity 2', 4, cdim*2,  &
          acc_solver % variable % perm)
      IF (cdim > 2) THEN
        CALL SetBoundaryConditions(acc_model, acc_amatrix, 'Re Velocity 3', 5, cdim*2,  &
            acc_solver % variable % perm)
        CALL SetBoundaryConditions(acc_model, acc_amatrix, 'Im Velocity 3', 6, cdim*2,  &
            acc_solver % variable % perm)
      END IF
    END IF
    CALL SetBoundaryConditions(acc_model, acc_smatrix, 'Re Temperature', 1, 4,  &
        acc_solver % variable % perm)
    CALL SetBoundaryConditions(acc_model, acc_smatrix, 'Im Temperature', 2, 4,  &
        acc_solver % variable % perm)
  end subroutine aca_block_setbcs

  subroutine aca_traction_precond(cdim, angfreq, sh, hr, dens, temp, cond, visc, lam, bulk, press, nn) &
      bind(c, name="aca_traction_precond")
    integer(c_int), value :: cdim, nn
    real(c_double), value :: angfreq
    real(c_double) :: sh(nn), hr(nn), dens(nn), temp(nn), cond(nn), visc(nn), &
        lam(nn), bulk(nn), press(nn)
    type(element_t), pointer :: CurrentElement, Parent
    type(valuelist_t), pointer :: Material
    integer, pointer :: NodeIndexes(:)
    logical :: stat, GotIt, SlipBoundary
    integer :: i, j, k, n, t, m
    real(kind=dp) :: reb, imb, rec, imc
    complex(kind=dp) :: C1, C2
    DO t = acc_solver % mesh % numberofbulkelements + 1, &
         acc_solver % mesh % numberofbulkelements + &
         acc_solver % mesh % numberofboundaryelements
       CurrentElement => acc_solver % mesh % elements(t)
       acc_model % currentelement => CurrentElement
       CurrentModel % currentelement => CurrentElement
       IF ( .NOT. ActiveBoundaryElement(CurrentElement, CurrentModel % solver ) ) CYCLE
       n = CurrentElement % TYPE % NumberOfNodes
       NodeIndexes => CurrentElement % NodeIndexes
       IF (ANY(acc_solver % variable % perm(NodeIndexes(1:n)) == 0)) CYCLE
       DO i=1,acc_model % numberofbcs
          IF ( CurrentElement % BoundaryInfo % Constraint == &
               acc_model % BCs(i) % Tag ) THEN
             SlipBoundary = ListGetLogical( acc_model % BCs(i) % values, &
                  'Traction Boundary', GotIt )
             IF ( .NOT. SlipBoundary) CYCLE
             Parent => CurrentELement % BoundaryInfo % Left
             stat = ASSOCIATED( Parent )
             IF (stat) stat = ALL( acc_solver % variable % perm( &
                 Parent % NodeIndexes(1: Parent % TYPE % NumberOfNodes)) > 0 )
             IF ( .NOT. stat) THEN
                Parent => CurrentELement % BoundaryInfo % Right
                stat = ASSOCIATED( Parent )
                IF (stat) stat = ALL( acc_solver % variable % perm( &
                    Parent % NodeIndexes(1: Parent % TYPE % NumberOfNodes)) > 0 )
                IF ( .NOT. stat )  CALL Fatal( 'AcousticsSolver', &
                     'No parent element can be found for given boundary element' )
             END IF
             k = ListGetInteger( acc_model % Bodies(Parent % Bodyid) % Values, &
                  'Material' )
             Material => acc_model % Materials(k) % Values
             sh(1:n) = ListGetReal( Material, 'Specific Heat', &
                  n, NodeIndexes )
             hr(1:n) = ListGetReal( Material, 'Specific Heat Ratio', &
                  n, NodeIndexes )
             dens(1:n) = ListGetReal( Material, 'Equilibrium Density', &
                  n, NodeIndexes )
             temp(1:n) = ListGetReal( Material, 'Equilibrium Temperature', &
                  n, NodeIndexes )
             cond(1:n) = ListGetReal( Material, 'Heat Conductivity', &
                  n, NodeIndexes )
             visc(1:n) = ListGetReal( Material, 'Viscosity', &
                  n, NodeIndexes )
             lam(1:n) = -2.0d0/3.0d0 * visc(1:n)
             bulk(1:n) = ListGetReal( Material, ' Bulk Viscosity', &
                  n, NodeIndexes, GotIt )
             IF (GotIt) lam(1:n) = bulk(1:n) - 2.0d0/3.0d0 * visc(1:n)
             press(1:n) = (hr(1:n)-1.0d0)* sh(1:n) * dens(1:n) * temp(1:n)
             C1 = CMPLX(0.0d0, -0.5d0*dens(1)*angfreq/visc(1), kind=dp)
             C2 = C1 * CMPLX( 1.0d0,angfreq/press(1)*(2.0d0*visc(1)+lam(1)), kind=dp ) / &
                  CMPLX( 1.0d0,angfreq/press(1)*lam(1), kind=dp )
             reb = REAL(C1)
             imb = AIMAG(C1)
             rec = REAL(C2)
             imc = AIMAG(C2)
             DO j=1,n
                m = acc_solver % variable % perm( NodeIndexes(j) )
                CALL ZeroRow( acc_solver % matrix, (m-1)*(cdim*2+4)+cdim*2+3 )
                CALL ZeroRow( acc_solver % matrix, (m-1)*(cdim*2+4)+cdim*2+4 )
                CALL SetMatrixElement( acc_solver % matrix, (m-1)*(cdim*2+4)+cdim*2+3, &
                  (m-1)*(cdim*2+4)+cdim*2+1, reb )
                CALL SetMatrixElement( acc_solver % matrix, (m-1)*(cdim*2+4)+cdim*2+3, &
                  (m-1)*(cdim*2+4)+cdim*2+2, -imb )
                CALL SetMatrixElement( acc_solver % matrix, (m-1)*(cdim*2+4)+cdim*2+3, &
                  (m-1)*(cdim*2+4)+cdim*2+3, rec )
                CALL SetMatrixElement( acc_solver % matrix, (m-1)*(cdim*2+4)+cdim*2+3, &
                  (m-1)*(cdim*2+4)+cdim*2+4, -imc )
                CALL SetMatrixElement( acc_solver % matrix, (m-1)*(cdim*2+4)+cdim*2+4, &
                  (m-1)*(cdim*2+4)+cdim*2+1, imb )
                CALL SetMatrixElement( acc_solver % matrix, (m-1)*(cdim*2+4)+cdim*2+4, &
                  (m-1)*(cdim*2+4)+cdim*2+2, reb )
                CALL SetMatrixElement( acc_solver % matrix, (m-1)*(cdim*2+4)+cdim*2+4, &
                  (m-1)*(cdim*2+4)+cdim*2+3, imc )
                CALL SetMatrixElement( acc_solver % matrix, (m-1)*(cdim*2+4)+cdim*2+4, &
                  (m-1)*(cdim*2+4)+cdim*2+4, rec )
             END DO
          END IF
       END DO
    END DO
  end subroutine aca_traction_precond
  ! -- shims into the verbatim-kept kernels -------------------------------------
  subroutine aca_parent_nodes() bind(c, name="aca_parent_nodes")
    integer :: n
    n = acc_parent % type % numberofnodes
    acc_pnodes % x(1:n) = acc_solver % mesh % nodes % x(acc_parent % nodeindexes)
    acc_pnodes % y(1:n) = acc_solver % mesh % nodes % y(acc_parent % nodeindexes)
    acc_pnodes % z(1:n) = acc_solver % mesh % nodes % z(acc_parent % nodeindexes)
  end subroutine aca_parent_nodes

  subroutine aca_shell_interface(gap, cdim, dofs, angfreq, nn) &
      bind(c, name="aca_shell_interface")
    integer(c_int), value :: cdim, dofs, nn
    real(c_double), value :: angfreq
    integer(c_int) :: gap(nn)
    call AcousticShellInterface(acc_model, acc_solver, &
        acc_solver % variable % values, acc_solver % variable % perm, gap, &
        cdim, dofs, angfreq, cdim, acc_solver % matrix, acc_solver % matrix % rhs)
  end subroutine aca_shell_interface

  subroutine aca_local_interface(lsm, lf, angfreq, sh, hr, dens, press, temp, &
      cond, imp, lod, nn, n, dofs) bind(c, name="aca_local_interface")
    integer(c_int), value :: nn, n, dofs
    real(c_double), value :: angfreq
    real(c_double) :: lsm(dofs*nn,dofs*nn), lf(dofs*nn)
    real(c_double) :: sh(nn), hr(nn), dens(nn), press(nn), temp(nn), &
        cond(nn), imp(4,nn), lod(6,nn)
    call LocalInterfaceMatrix(lsm, lf, angfreq, sh, hr, dens, press, temp, &
        cond, imp, lod, acc_elem, n, acc_nodes, dofs)
  end subroutine aca_local_interface

  subroutine aca_slip_matrix(lsm, lf, angfreq, sh, hr, dens, cond, press, &
      temp, wallt, wallv, c1, c2, c3, nn, n, dofs) bind(c, name="aca_slip_matrix")
    integer(c_int), value :: nn, n, dofs
    real(c_double), value :: angfreq, c1, c2, c3
    real(c_double) :: lsm(dofs*nn,dofs*nn), lf(dofs*nn)
    real(c_double) :: sh(nn), hr(nn), dens(nn), cond(nn), press(nn), &
        temp(nn), wallt(nn), wallv(6,nn)
    call SlipMatrix(lsm, lf, sh, hr, dens, cond, press, temp, angfreq, &
        wallt, wallv, c1, c2, c3, acc_elem, n, acc_nodes, dofs)
  end subroutine aca_slip_matrix

  subroutine aca_avg_velocity(n, WallVelo, wsz, re, im) bind(c, name="aca_avg_velocity")
    integer(c_int), value :: n, wsz
    real(c_double) :: WallVelo(6,wsz)
    real(c_double) :: re, im
    complex(kind=dp) :: res
    call ComputeAverageVelocity(acc_elem, n, WallVelo, res, acc_solver)
    re = real(res, kind=dp)
    im = aimag(res)
  end subroutine aca_avg_velocity

  subroutine aca_surface_force(traction, moment, momentabout, area, calc, &
      wallvelo, lod, visc, bulk, npn, nn) bind(c, name="aca_surface_force")
    integer(c_int), value :: calc, npn, nn
    real(c_double) :: traction(6), moment(6), momentabout(3), area
    real(c_double) :: wallvelo(6,nn), lod(2,nn), visc(nn), bulk(nn)
    logical :: lcalc
    lcalc = (calc /= 0)
    call SurfaceForceIntegration(acc_elem, acc_parent, traction, moment, &
        momentabout, area, lcalc, wallvelo, lod, visc, bulk, acc_nodes, &
        acc_pnodes, npn, acc_solver, CurrentCoordinateSystem())
  end subroutine aca_surface_force

  subroutine aca_setdir_points(dof, ndofs, n, gap, vals, nn) &
      bind(c, name="aca_setdir_points")
    integer(c_int), value :: dof, ndofs, n, nn
    integer(c_int) :: gap(nn)
    real(c_double) :: vals(nn)
    call SetDirichletPoints(acc_solver % matrix, acc_solver % matrix % rhs, &
        dof, ndofs, acc_solver % variable % perm, n, gap, vals)
  end subroutine aca_setdir_points

  ! -- Helmholtz-solution and misc queries --------------------------------------
  subroutine aca_dd_pres() bind(c, name="aca_dd_pres")
    type(variable_t), pointer :: hsol
    hsol => VariableGet(acc_model % variables, 'Pres')
    if (associated(hsol)) then
      acc_h => hsol % values
      acc_hp => hsol % perm
    else
      call Fatal('AcousticsSolver', 'Helmholtz solution Pres is not available')
    end if
  end subroutine aca_dd_pres

  subroutine aca_h_at(k, re, im) bind(c, name="aca_h_at")
    integer(c_int), value :: k
    real(c_double) :: re, im
    re = acc_h(2*k-1)
    im = acc_h(2*k)
  end subroutine aca_h_at

  function aca_hp_at(node) result(v) bind(c, name="aca_hp_at")
    integer(c_int), value :: node
    integer(c_int) :: v
    v = acc_hp(node)
  end function aca_hp_at

  function aca_ain_size() result(n) bind(c, name="aca_ain_size")
    integer(c_int) :: n
    if (allocated(acc_ainodes)) then
      n = size(acc_ainodes)
    else
      n = 0
    end if
  end function aca_ain_size

  subroutine aca_ain_at(k, node, r1, r2) bind(c, name="aca_ain_at")
    integer(c_int), value :: k
    integer(c_int) :: node
    real(c_double) :: r1, r2
    node = acc_ainodes(k)
    r1 = acc_aires(k,1)
    r2 = acc_aires(k,2)
  end subroutine aca_ain_at

  function aca_elem_code() result(c) bind(c, name="aca_elem_code")
    integer(c_int) :: c
    c = acc_elem % type % elementcode
  end function aca_elem_code

  function aca_h_found() result(r) bind(c, name="aca_h_found")
    integer(c_int) :: r
    r = merge(1, 0, associated(acc_h))
  end function aca_h_found

end module acoustics_acc
