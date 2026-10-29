! StructuredFlowLine port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/StructuredFlowLine.F90, compiled against the
! built .mod files of the target configure.
!
! Every wrapper takes only C-interoperable arguments (TYPE(C_PTR), VALUE
! scalars, explicit-shape REAL(C_DOUBLE)/INTEGER(C_INT) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module mirrors the original's locals and SAVE state with identical
! lifetimes:
!   - outer-procedure SAVE: VisitedTimes, Initialized, UnitVector, Coord,
!     TopPointer, BotPointer, DownPointer, UpPointer, Dim, HardDisp,
!     DispPerm, MaskPerm, MaskExist, ActiveDirection (TopMode, BotMode,
!     TopField, BotField, TopPerm, BotPerm, Field, Surface and DOFs are
!     declared+SAVE'd in the original but never used — NULLIFY(Field)
!     aside — so they need no mirror);
!   - internal-procedure SAVE of ComputeVeloAtPoint: Visited, Elements,
!     RootQuadrant, CurrentElement, ElementNodes, Basis, dBasisdx;
!   - timing locals at0/at1 (CPUTime) feeding only the wall-clock Info
!     lines (filtered by the proof norm).
! Order-sensitive reductions (SUM, MATMUL-free here), whole-array updates
! and every formatted WRITE/Info stay in Fortran wrappers, so
! bit-identity holds by construction wherever the call sequence matches.
!
! Dead code elided (all guarded by IF(.FALSE.) in the original and never
! executed): the three bare PRINT markers and the debug PRINT block in
! ComputeVeloAtPoint.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module sfl_acc
  use iso_c_binding
  use types
  use lists
  use defutils
  use meshbasics
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(mesh_t), pointer, save :: acc_mesh => null()
  type(valuelist_t), pointer, save :: acc_params => null()
  type(element_t), pointer, save :: acc_element => null()
  ! Mirrors the outer procedure's `Nodes` SAVE (GetElementNodes target).
  type(nodes_t), target, save :: acc_nodes

  ! -- outer-procedure SAVE mirrors ------------------------------------
  integer(c_int), save :: s_visited = 0
  logical, save :: s_initialized = .false.
  real(kind=dp), save :: s_unit(3) = 0.0_dp
  real(kind=dp), pointer, save :: s_coord(:) => null()
  integer, pointer, save :: s_top(:) => null(), s_up(:) => null(), &
      s_bot(:) => null(), s_down(:) => null(), s_maskperm(:) => null()
  logical, save :: s_maskexist = .false.
  integer(c_int), save :: s_actdir = 0, s_dim = 0
  ! Per-call associations (re-associated on every entry, as original).
  real(kind=dp), pointer, save :: s_flow(:) => null(), &
      s_hard(:) => null(), s_hitf(:) => null()
  integer, pointer, save :: s_flowperm(:) => null(), &
      s_dispperm(:) => null(), s_hitperm(:) => null()
  integer(c_int), save :: s_flowdofs = 0, s_dispdofs = 0
  logical, save :: s_gothit = .false.
  real(kind=dp), save :: s_at0 = 0.0_dp

  ! -- ComputeVeloAtPoint SAVE mirrors ----------------------------------
  type sfl_elemptr
    type(element_t), pointer :: el => null()
  end type sfl_elemptr
  logical, save :: v_visited = .false.
  type(sfl_elemptr), allocatable, save :: v_elems(:)
  type(quadrant_t), pointer, save :: v_root => null()
  type(element_t), pointer, save :: v_cur => null()
  type(nodes_t), save :: v_nodes
  real(kind=dp), allocatable, save :: v_basis(:), v_dbasis(:,:)

  character(*), parameter :: Caller = 'StructuredFlowLine'

contains

  ! -- entry state ------------------------------------------------------
  subroutine sfl_begin(h_model, h_solver) bind(c, name="sfl_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    acc_mesh => acc_solver % mesh
    acc_params => GetSolverParams()
  end subroutine sfl_begin

  subroutine sfl_info_header() bind(c, name="sfl_info_header")
    call Info(Caller, '---------------------------------------', Level=4)
    call Info(Caller, 'Setting edges to follow the streamlines ', Level=4)
    call Info(Caller, '---------------------------------------', Level=4)
  end subroutine sfl_info_header

  subroutine sfl_free_quadrant() bind(c, name="sfl_free_quadrant")
    call FreeQuadrantTree(acc_mesh % RootQuadrant)
    acc_solver % Mesh % RootQuadrant => null()
  end subroutine sfl_free_quadrant

  function sfl_get_true_iters() result(i) bind(c, name="sfl_get_true_iters")
    integer(c_int) :: i
    i = GetInteger(acc_params, 'True Flow Line Iterations')
  end function sfl_get_true_iters

  function sfl_visited() result(v) bind(c, name="sfl_visited")
    integer(c_int) :: v
    v = s_visited
  end function sfl_visited

  subroutine sfl_mark_visited() bind(c, name="sfl_mark_visited")
    s_visited = s_visited + 1
  end subroutine sfl_mark_visited

  function sfl_get_ave_order(tl) result(o) bind(c, name="sfl_get_ave_order")
    integer(c_int), value :: tl
    integer(c_int) :: o
    logical :: gotit
    o = GetInteger(acc_params, 'Averaging Order', gotit)
    if (tl /= 0) o = max(1, o)
  end function sfl_get_ave_order

  function sfl_get_ave_mode() result(m) bind(c, name="sfl_get_ave_mode")
    integer(c_int) :: m
    logical :: gotit
    m = GetInteger(acc_params, 'Averaging Method', gotit)
  end function sfl_get_ave_mode

  function sfl_initialized() result(r) bind(c, name="sfl_initialized")
    integer(c_int) :: r
    r = merge(1, 0, s_initialized)
  end function sfl_initialized

  subroutine sfl_mark_initialized() bind(c, name="sfl_mark_initialized")
    s_initialized = .true.
  end subroutine sfl_mark_initialized

  function sfl_enter_init() result(d) bind(c, name="sfl_enter_init")
    integer(c_int) :: d
    s_dim = CoordinateSystemDimension()
    d = s_dim
  end function sfl_enter_init

  function sfl_get_dim() result(d) bind(c, name="sfl_get_dim")
    integer(c_int) :: d
    d = s_dim
  end function sfl_get_dim

  function sfl_get_actdir() result(a) bind(c, name="sfl_get_actdir")
    integer(c_int) :: a
    a = s_actdir
  end function sfl_get_actdir

  function sfl_is_mask() result(r) bind(c, name="sfl_is_mask")
    integer(c_int) :: r
    r = merge(1, 0, s_maskexist)
  end function sfl_is_mask

  function sfl_init_coord() result(ad) bind(c, name="sfl_init_coord")
    integer(c_int) :: ad
    integer :: iact
    iact = GetInteger(acc_params, 'Active Coordinate')
    if (iact == 1) then
      s_coord => acc_mesh % Nodes % x
    else if (iact == 2) then
      s_coord => acc_mesh % Nodes % y
    else if (iact == 3) then
      s_coord => acc_mesh % Nodes % z
    else
      call Fatal(Caller, 'Invalid value for Active Coordinate')
    end if
    s_actdir = iact
    ad = iact
  end function sfl_init_coord

  subroutine sfl_init_unit() bind(c, name="sfl_init_unit")
    s_unit = 0.0_dp
    s_unit(s_actdir) = 1.0_dp
  end subroutine sfl_init_unit

  function sfl_init_mask() result(r) bind(c, name="sfl_init_mask")
    integer(c_int) :: r
    character(len=MAX_NAME_LEN) :: maskname
    integer :: n
    maskname = 'Flow Line'
    n = acc_mesh % NumberOfNodes
    allocate(s_maskperm(n))
    n = 0
    s_maskperm = 0
    call MakePermUsingMask(acc_model, acc_solver, acc_mesh, &
        trim(maskname), .false., s_maskperm, n)
    if (n == 0) then
      s_maskexist = .false.
      deallocate(s_maskperm)
    else
      s_maskexist = .true.
      write(message, '(A,I6)') 'Number of nodes on '//trim(maskname), n
      call Info(Caller, message)
    end if
    r = merge(1, 0, s_maskexist)
  end function sfl_init_mask

  subroutine sfl_alloc_pointers() bind(c, name="sfl_alloc_pointers")
    integer :: i, nsize
    nsize = size(s_coord)
    allocate(s_top(nsize), s_up(nsize))
    do i = 1, nsize
      s_top(i) = i
      s_up(i) = i
    end do
    allocate(s_bot(nsize), s_down(nsize))
    do i = 1, nsize
      s_bot(i) = i
      s_down(i) = i
    end do
  end subroutine sfl_alloc_pointers

  function sfl_get_eps() result(e) bind(c, name="sfl_get_eps")
    real(c_double) :: e
    logical :: gotit
    e = GetConstReal(acc_params, 'Dot Product Tolerance', gotit)
    if (.not. gotit) e = 1.0e-4
  end function sfl_get_eps

  ! -- init element-loop primitives --------------------------------------
  function sfl_num_bulk() result(n) bind(c, name="sfl_num_bulk")
    integer(c_int) :: n
    n = acc_mesh % NumberOfBulkElements
  end function sfl_num_bulk

  function sfl_max_elem_nodes() result(n) bind(c, name="sfl_max_elem_nodes")
    integer(c_int) :: n
    n = acc_mesh % MaxElementNodes
  end function sfl_max_elem_nodes

  function sfl_elem_begin(elem) result(n) bind(c, name="sfl_elem_begin")
    integer(c_int), value :: elem
    integer(c_int) :: n
    acc_element => acc_mesh % elements(elem)
    acc_model % CurrentElement => acc_element
    call GetElementNodes(acc_nodes)
    n = GetElementNOFNodes()
  end function sfl_elem_begin

  subroutine sfl_elem_nodeidx(n, buf) bind(c, name="sfl_elem_nodeidx")
    integer(c_int), value :: n
    integer(c_int) :: buf(n)
    buf(1:n) = acc_element % NodeIndexes(1:n)
  end subroutine sfl_elem_nodeidx

  subroutine sfl_elem_coords(n, xs, ys, zs) bind(c, name="sfl_elem_coords")
    integer(c_int), value :: n
    real(c_double) :: xs(n), ys(n), zs(n)
    xs(1:n) = acc_nodes % x(1:n)
    ys(1:n) = acc_nodes % y(1:n)
    zs(1:n) = acc_nodes % z(1:n)
  end subroutine sfl_elem_coords

  function sfl_mask_at(i) result(m) bind(c, name="sfl_mask_at")
    integer(c_int), value :: i
    integer(c_int) :: m
    m = s_maskperm(i)
  end function sfl_mask_at

  function sfl_vec_len(ev) result(l) bind(c, name="sfl_vec_len")
    real(c_double) :: ev(3)
    real(c_double) :: l
    l = sqrt(sum(ev * ev))
  end function sfl_vec_len

  function sfl_vec_dot(ev) result(d) bind(c, name="sfl_vec_dot")
    real(c_double) :: ev(3)
    real(c_double) :: d
    d = sum(ev * s_unit)
  end function sfl_vec_dot

  subroutine sfl_set_updown(ii, jj, dir) bind(c, name="sfl_set_updown")
    integer(c_int), value :: ii, jj, dir
    if (dir > 0) then
      s_up(ii) = jj
      s_down(jj) = ii
    else
      s_down(ii) = jj
      s_up(jj) = ii
    end if
  end subroutine sfl_set_updown

  ! Verbatim top/bottom pointer propagation plus its two Info reports
  ! (the timing line is wall-clock, filtered by the proof norm).
  function sfl_propagate() result(rounds) bind(c, name="sfl_propagate")
    integer(c_int) :: rounds
    logical :: found
    integer :: i, j, nsize
    real(kind=dp) :: at1
    nsize = size(s_coord)
    found = .true.
    rounds = 0
    do while (found)
      found = .false.
      do i = 1, nsize
        if (s_maskexist) then
          if (s_maskperm(i) == 0) cycle
        end if
        j = s_up(i)
        if (s_top(i) /= s_top(j)) then
          found = .true.
          s_top(i) = s_top(j)
        end if
        j = s_down(i)
        if (s_bot(i) /= s_bot(j)) then
          found = .true.
          s_bot(i) = s_bot(j)
        end if
      end do
      if (found) rounds = rounds + 1
    end do
    at1 = CPUTime()
    write(message, *) 'Top and bottom pointer init rounds: ', rounds
    call Info(Caller, message)
    write(message, '(A,F8.3)') 'Top and bottom pointer init time: ', at1 - s_at0
    call Info(Caller, message)
  end function sfl_propagate

  subroutine sfl_tic() bind(c, name="sfl_tic")
    s_at0 = CPUTime()
  end subroutine sfl_tic

  ! -- per-call variable lookups ------------------------------------------
  function sfl_get_flow() result(dofs) bind(c, name="sfl_get_flow")
    integer(c_int) :: dofs
    character(len=MAX_NAME_LEN) :: varname
    type(variable_t), pointer :: var
    logical :: gotit
    varname = GetString(acc_params, 'Velocity Variable Name', gotit)
    if (.not. gotit) varname = 'Flow Solution'
    var => VariableGet(acc_mesh % Variables, varname)
    if (associated(var)) then
      s_flow => var % Values
      s_flowperm => var % Perm
      s_flowdofs = var % DOFs
    else
      call Fatal(Caller, 'Velocity variable is missing: '//trim(varname))
    end if
    dofs = s_flowdofs
  end function sfl_get_flow

  function sfl_get_disp() result(dofs) bind(c, name="sfl_get_disp")
    integer(c_int) :: dofs
    character(len=MAX_NAME_LEN) :: varname
    type(variable_t), pointer :: var
    logical :: gotit
    varname = GetString(acc_params, 'Flow Line Displacement Name', gotit)
    if (.not. gotit) varname = 'FlowLineDisp'
    nullify(var)
    var => VariableGet(acc_mesh % Variables, varname, ThisOnly=.true.)
    if (associated(var)) then
      s_hard => var % Values
      s_dispperm => var % Perm
      s_dispdofs = var % DOFs
    else
      call Fatal(Caller, 'Flow line displacement variable does not exist!')
    end if
    if (s_dispdofs /= s_dim) then
      call Fatal(Caller, 'Displacement variable of wrong size')
    end if
    dofs = s_dispdofs
  end function sfl_get_disp

  function sfl_get_hitstyle() result(r) bind(c, name="sfl_get_hitstyle")
    integer(c_int) :: r
    character(len=MAX_NAME_LEN) :: varname
    type(variable_t), pointer :: var
    logical :: gotit
    varname = GetString(acc_params, 'Hit Style Name', gotit)
    if (.not. gotit) varname = 'HitStyle'
    var => VariableGet(acc_mesh % Variables, varname)
    s_gothit = associated(var)
    if (s_gothit) then
      s_hitf => var % Values
      s_hitperm => var % Perm
    end if
    r = merge(1, 0, s_gothit)
  end function sfl_get_hitstyle

  function sfl_get_dispmode() result(r) bind(c, name="sfl_get_dispmode")
    integer(c_int) :: r
    logical :: found, l
    l = GetLogical(acc_params, 'Displacement Mode', found)
    r = merge(1, 0, l)
  end function sfl_get_dispmode

  subroutine sfl_coord_limits(cmin, gotmin, cmax, gotmax) &
      bind(c, name="sfl_coord_limits")
    real(c_double) :: cmin, cmax
    integer(c_int) :: gotmin, gotmax
    logical :: f1, f2
    cmin = ListGetConstReal(acc_solver % Values, 'Min Coordinate', f1)
    if (.not. f1) cmin = -huge(cmin)
    cmax = ListGetConstReal(acc_solver % Values, 'Max Coordinate', f2)
    if (.not. f2) cmax = huge(cmax)
    gotmin = merge(1, 0, f1)
    gotmax = merge(1, 0, f2)
  end subroutine sfl_coord_limits

  function sfl_nsize() result(n) bind(c, name="sfl_nsize")
    integer(c_int) :: n
    n = size(s_coord)
  end function sfl_nsize

  ! -- streamline-loop primitives ------------------------------------------
  function sfl_disp_perm(i) result(p) bind(c, name="sfl_disp_perm")
    integer(c_int), value :: i
    integer(c_int) :: p
    p = s_dispperm(i)
  end function sfl_disp_perm

  function sfl_top_at(i) result(t) bind(c, name="sfl_top_at")
    integer(c_int), value :: i
    integer(c_int) :: t
    t = s_top(i)
  end function sfl_top_at

  function sfl_down_at(k) result(d) bind(c, name="sfl_down_at")
    integer(c_int), value :: k
    integer(c_int) :: d
    d = s_down(k)
  end function sfl_down_at

  function sfl_bot_at(k) result(b) bind(c, name="sfl_bot_at")
    integer(c_int), value :: k
    integer(c_int) :: b
    b = s_bot(k)
  end function sfl_bot_at

  subroutine sfl_flow_at(k, dim, buf) bind(c, name="sfl_flow_at")
    integer(c_int), value :: k, dim
    real(c_double) :: buf(3)
    integer :: i
    do i = 1, dim
      buf(i) = s_flow(s_flowdofs * (s_flowperm(k) - 1) + i)
    end do
  end subroutine sfl_flow_at

  subroutine sfl_node_xyz(k, buf) bind(c, name="sfl_node_xyz")
    integer(c_int), value :: k
    real(c_double) :: buf(3)
    buf(1) = acc_mesh % Nodes % x(k)
    buf(2) = acc_mesh % Nodes % y(k)
    buf(3) = acc_mesh % Nodes % z(k)
  end subroutine sfl_node_xyz

  subroutine sfl_mark_hit(k, val) bind(c, name="sfl_mark_hit")
    integer(c_int), value :: k
    real(c_double), value :: val
    s_hitf(s_hitperm(k)) = val
  end subroutine sfl_mark_hit

  function sfl_norm3(a) result(s) bind(c, name="sfl_norm3")
    real(c_double) :: a(3)
    real(c_double) :: s
    s = sqrt(sum(a ** 2))
  end function sfl_norm3

  function sfl_norm3_sum(a, b) result(s) bind(c, name="sfl_norm3_sum")
    real(c_double) :: a(3), b(3)
    real(c_double) :: s
    s = sqrt(sum((a + b) ** 2))
  end function sfl_norm3_sum

  subroutine sfl_harddisp_set(k, dim, s) bind(c, name="sfl_harddisp_set")
    integer(c_int), value :: k, dim
    real(c_double) :: s(3)
    integer :: l
    do l = 1, dim
      s_hard(s_dispdofs * (s_dispperm(k) - 1) + l) = s(l)
    end do
  end subroutine sfl_harddisp_set

  subroutine sfl_report_norm() bind(c, name="sfl_report_norm")
    real(kind=dp) :: norm
    norm = sqrt(sum(s_hard ** 2) / maxval(s_dispperm))
    write(message, *) 'Result Norm: ', norm
    call Info(Caller, message)
  end subroutine sfl_report_norm

  subroutine sfl_relax_apply() bind(c, name="sfl_relax_apply")
    real(kind=dp) :: relax
    logical :: gotit
    relax = GetCReal(acc_solver % Values, &
        'Nonlinear System Relaxation Factor', gotit)
    if (gotit .and. abs(relax - 1.0_dp) > tiny(relax)) then
      s_hard = relax * s_hard
    end if
  end subroutine sfl_relax_apply

  ! Verbatim displacement block, including the reference's leftover-`i`
  ! quirk: the mask guard reads the outer i-loop's post-loop value
  ! (nsize+1, out of bounds) instead of the k-loop index. The caller
  ! passes that leftover value through; the proof keeps defined paths
  ! (Displacement Mode off, or no mask) for the differentials.
  subroutine sfl_move_nodes(leftover_i) bind(c, name="sfl_move_nodes")
    integer(c_int), value :: leftover_i
    integer :: k, i
    i = leftover_i
    call Info(Caller, 'Moving the nodes')
    do k = 1, size(s_coord)
      if (s_maskexist) then
        if (s_maskperm(i) == 0) cycle
      end if
      acc_mesh % Nodes % x(k) = acc_mesh % Nodes % x(k) + &
          s_hard(s_dispdofs * (s_dispperm(k) - 1) + 1)
      acc_mesh % Nodes % y(k) = acc_mesh % Nodes % y(k) + &
          s_hard(s_dispdofs * (s_dispperm(k) - 1) + 2)
      if (s_dim == 3) acc_mesh % Nodes % z(k) = acc_mesh % Nodes % z(k) + &
          s_hard(s_dispdofs * (s_dispperm(k) - 1) + 3)
    end do
  end subroutine sfl_move_nodes

  subroutine sfl_report_shots(shots, hits) bind(c, name="sfl_report_shots")
    integer(c_int), value :: shots, hits
    write(message, *) 'Velocity estimations, total: ', shots
    call Info(Caller, message)
    write(message, *) 'Velocity estimations, inside: ', hits
    call Info(Caller, message)
  end subroutine sfl_report_shots

  subroutine sfl_report_maptime() bind(c, name="sfl_report_maptime")
    real(kind=dp) :: at1
    at1 = CPUTime()
    write(message, '(A,F8.3)') 'Streamline mapping time: ', at1 - s_at0
    call Info(Caller, message)
  end subroutine sfl_report_maptime

  ! -- ComputeVeloAtPoint (verbatim internal procedure) ---------------------
  ! Called through sfl_velo_at_point below; WithGrad selects the PRESENT
  ! shape of the optional GradVeloAtPoint dummy.
  subroutine sfl_velo_core(GlobalCoords, node, OutOfScope, VeloAtPoint, &
      GradVeloAtPoint, WithGrad)
    real(kind=dp) :: GlobalCoords(3)
    integer :: node
    logical :: OutOfScope
    real(kind=dp) :: VeloAtPoint(3), GradVeloAtPoint(3, 3)
    logical :: WithGrad
    logical :: Stat, Hit, DummySearch
    integer :: i, j, k, n
    real(kind=dp) :: u, v, w, LocalCoords(3), SqrtElementMetric
    real(kind=dp) :: GradientAtPoint(3)
    integer, pointer :: NodeIndexes(:)
    type(element_t), pointer :: CurrentElement
    type(quadrant_t), pointer :: LeafQuadrant
    real(kind=dp) :: BoundingBox(6), eps2, eps1

    eps1 = 1.0e-3
    if (.not. v_visited) then
      n = acc_solver % mesh % MaxElementNodes
      allocate(v_nodes % x(n), v_nodes % y(n), v_nodes % z(n))
      allocate(v_basis(n), v_dbasis(n, 3))
      allocate(v_elems(acc_mesh % NumberOfNodes))
      nullify(v_cur)
      v_visited = .true.
    end if

    GradientAtPoint = 0.0_dp
    VeloAtPoint = 0.0_dp
    Hit = .false.
    DummySearch = .false.

    CurrentElement => v_elems(node) % el
    if (associated(CurrentElement)) then
      n = CurrentElement % TYPE % NumberOfNodes
      NodeIndexes => CurrentElement % NodeIndexes
      v_nodes % x(1:n) = acc_mesh % Nodes % x(NodeIndexes)
      v_nodes % y(1:n) = acc_mesh % Nodes % y(NodeIndexes)
      v_nodes % z(1:n) = acc_mesh % Nodes % z(NodeIndexes)
      if (PointInElement(CurrentElement, v_nodes, &
          GlobalCoords, LocalCoords)) then
        Hit = .true.
      end if
    end if

    if (.not. Hit) then
      if (s_dim == 2 .or. DummySearch) then
        do k = 1, acc_solver % mesh % NumberOfBulkElements
          CurrentElement => acc_mesh % Elements(k)
          n = CurrentElement % TYPE % NumberOfNodes
          NodeIndexes => CurrentElement % NodeIndexes
          v_nodes % x(1:n) = acc_mesh % Nodes % x(NodeIndexes)
          v_nodes % y(1:n) = acc_mesh % Nodes % y(NodeIndexes)
          v_nodes % z(1:n) = acc_mesh % Nodes % z(NodeIndexes)
          if (PointInElement(CurrentElement, v_nodes, &
              GlobalCoords, LocalCoords)) then
            Hit = .true.
            exit
          end if
        end do
      else
        if (.not. associated(acc_mesh % RootQuadrant)) then
          BoundingBox(1) = minval(acc_mesh % Nodes % x)
          BoundingBox(2) = minval(acc_mesh % Nodes % y)
          BoundingBox(3) = minval(acc_mesh % Nodes % z)
          BoundingBox(4) = maxval(acc_mesh % Nodes % x)
          BoundingBox(5) = maxval(acc_mesh % Nodes % y)
          BoundingBox(6) = maxval(acc_mesh % Nodes % z)
          eps2 = eps1 * maxval(BoundingBox(4:6) - BoundingBox(1:3))
          BoundingBox(1:3) = BoundingBox(1:3) - eps2
          BoundingBox(4:6) = BoundingBox(4:6) + eps2
          call BuildQuadrantTree(acc_mesh, BoundingBox, &
              acc_mesh % RootQuadrant)
        end if
        v_root => acc_mesh % RootQuadrant
        if (associated(v_root)) then
          nullify(CurrentElement)
          call FindLeafElements(GlobalCoords, s_dim, v_root, LeafQuadrant)
          if (associated(LeafQuadrant)) then
            do k = 1, LeafQuadrant % NElemsInQuadrant
              CurrentElement => acc_mesh % Elements(LeafQuadrant % Elements(k))
              n = CurrentElement % TYPE % NumberOfNodes
              NodeIndexes => CurrentElement % NodeIndexes
              v_nodes % x(1:n) = acc_mesh % Nodes % x(NodeIndexes)
              v_nodes % y(1:n) = acc_mesh % Nodes % y(NodeIndexes)
              v_nodes % z(1:n) = acc_mesh % Nodes % z(NodeIndexes)
              if (PointInElement(CurrentElement, v_nodes, &
                  GlobalCoords, LocalCoords)) then
                Hit = .true.
                exit
              end if
            end do
          end if
        end if
      end if
    end if

    if (Hit) then
      v_elems(node) % el => CurrentElement
      u = LocalCoords(1)
      v = LocalCoords(2)
      w = LocalCoords(3)
      Stat = ElementInfo(CurrentElement, v_nodes, U, V, W, &
          SqrtElementMetric, v_basis, v_dbasis)
      if (all(s_flowperm(NodeIndexes) > 0)) then
        do i = 1, s_dim
          VeloAtPoint(i) = sum(v_basis(1:n) * &
              s_flow(s_flowdofs * (s_flowperm(NodeIndexes) - 1) + i))
          if (WithGrad) then
            do j = 1, s_dim
              GradVeloAtPoint(i, j) = sum(v_dbasis(1:n, j) * &
                  s_flow(s_flowdofs * (s_flowperm(NodeIndexes) - 1) + i))
            end do
          end if
        end do
      end if
    end if

    OutOfScope = .not. Hit
  end subroutine sfl_velo_core

  subroutine sfl_velo_at_point(g, node, with_grad, velo, grad, outmesh) &
      bind(c, name="sfl_velo_at_point")
    real(c_double) :: g(3), velo(3), grad(9)
    integer(c_int), value :: node, with_grad
    integer(c_int) :: outmesh
    real(kind=dp) :: gc(3), vv(3), gg(3, 3)
    logical :: oos
    gc = g
    call sfl_velo_core(gc, node, oos, vv, gg, with_grad /= 0)
    velo = vv
    if (with_grad /= 0) grad = reshape(gg, [9])
    outmesh = merge(1, 0, oos)
  end subroutine sfl_velo_at_point

  ! -- StructuredFlowLine_init entry state -----------------------------------
  subroutine sfl_init_begin(h_model, h_solver, dim) &
      bind(c, name="sfl_init_begin")
    type(c_ptr), value :: h_model, h_solver
    integer(c_int) :: dim
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    acc_mesh => acc_solver % mesh
    acc_params => GetSolverParams()
    dim = CoordinateSystemDimension()
  end subroutine sfl_init_begin

  subroutine sfl_init_nomatrix() bind(c, name="sfl_init_nomatrix")
    call ListAddLogical(acc_params, 'No Matrix', .true.)
  end subroutine sfl_init_nomatrix

  subroutine sfl_init_var(dim) bind(c, name="sfl_init_var")
    integer(c_int), value :: dim
    character(len=MAX_NAME_LEN) :: varname
    logical :: gotit
    varname = GetString(acc_params, 'Flow Line Displacement Name', gotit)
    if (.not. gotit) varname = 'FlowLineDisp'
    call ListAddString(acc_params, 'Variable', varname)
    call ListAddInteger(acc_params, 'Variable Dofs', dim)
  end subroutine sfl_init_var

end module sfl_acc
