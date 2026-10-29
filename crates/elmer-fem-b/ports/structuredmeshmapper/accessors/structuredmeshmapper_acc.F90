! StructuredMeshMapper port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/StructuredMeshMapper.F90, compiled against the
! built .mod files of the target configure (A2 shipped mechanism, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module mirrors the original's SAVE Fortran-object state (pointer
! associations, ALLOCATABLE arrays, variable handles) with identical
! lifetimes; plain SAVE scalars (Visited, Initialized, TopMode, BotMode,
! nsize, nnodes, Cnt, ...) live in Rust.
!
! Order-sensitive reductions (MATMUL in the multilayer projector,
! ParallelReduction over the tangled/limited counters), formatted I/O
! (PRINT/WRITE, Message construction incl. I2S/E11.4 edits), and
! whole-array updates (OrigCoord/Coord sync, VeloVar zeroing, stride
! commit/WHERE blocks, BC field gathers) stay in Fortran wrappers, so
! bit-identity holds by construction wherever the call sequence matches.
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module structuredmeshmapper_acc
  use iso_c_binding
  use types
  use defutils
  use meshtransform
  use messages
  use generalutils
  use parallelutils
  implicit none

  character(*), parameter :: Caller = 'StructuredMeshMapper'

  type(solver_t), pointer, save :: acc_solver => null()
  type(model_t), pointer, save :: acc_model => null()
  type(valuelist_t), pointer, save :: acc_params => null()
  type(mesh_t), pointer, save :: acc_mesh => null()
  type(variable_t), pointer, save :: acc_var => null()
  type(variable_t), pointer, save :: acc_velovar => null()
  type(variable_t), pointer, save :: acc_updatevar => null()
  type(variable_t), pointer, save :: acc_tangledmaskvar => null()
  type(variable_t), pointer, save :: acc_basevar => null()
  type(variable_t), pointer, save :: acc_fixedvar => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(valuelist_t), pointer, save :: acc_bc => null()
  real(kind=dp), pointer, save :: acc_coord(:) => null()
  real(kind=dp), pointer, save :: acc_topfield(:) => null()
  real(kind=dp), pointer, save :: acc_botfield(:) => null()
  real(kind=dp), pointer, save :: acc_tangledmask(:) => null()
  real(kind=dp), pointer, save :: acc_coordp(:) => null()
  integer, pointer, save :: acc_maskperm(:) => null()
  integer, pointer, save :: acc_topperm(:) => null()
  integer, pointer, save :: acc_botperm(:) => null()
  integer, pointer, save :: acc_tangledmaskperm(:) => null()
  integer, pointer, save :: acc_toppointer(:) => null()
  integer, pointer, save :: acc_botpointer(:) => null()
  integer, pointer, save :: acc_midpointer(:) => null()
  integer, pointer, save :: acc_uppointer(:) => null()
  integer, pointer, save :: acc_downpointer(:) => null()
  integer, pointer, save :: acc_nodelayer(:) => null()
  integer, pointer, save :: acc_fixed(:) => null()
  integer, pointer, save :: acc_tmpperm(:) => null()
  real(kind=dp), allocatable, save :: acc_origcoord(:)
  real(kind=dp), allocatable, save :: acc_field(:)
  real(kind=dp), allocatable, save :: acc_surface(:)

contains

  ! Copy a C buffer (stops at NUL) into a Fortran string.
  subroutine sma_c2f(cstr, n, fstr)
    integer, intent(in) :: n
    character(kind=c_char), intent(in) :: cstr(n)
    character(len=*), intent(out) :: fstr
    integer :: k, m
    fstr = ' '
    m = min(n, len(fstr))
    do k = 1, m
      if (cstr(k) == c_null_char) exit
      fstr(k:k) = cstr(k)
    end do
  end subroutine sma_c2f

  ! -- entry state ----------------------------------------------------
  subroutine sma_begin(h_model, h_solver) bind(c, name="sma_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    acc_mesh => acc_solver % mesh
    acc_params => GetSolverParams()
  end subroutine sma_begin

  subroutine sma_add_nomatrix() bind(c, name="sma_add_nomatrix")
    call ListAddNewLogical(acc_params, 'No Matrix', .true.)
  end subroutine sma_add_nomatrix

  subroutine sma_info_banner() bind(c, name="sma_info_banner")
    call Info(Caller, '---------------------------------------', Level=4)
    call Info(Caller, 'Performing mapping on a structured mesh ', Level=4)
    call Info(Caller, '---------------------------------------', Level=4)
  end subroutine sma_info_banner

  function sma_mesh_dim() result(d) bind(c, name="sma_mesh_dim")
    integer(c_int) :: d
    d = acc_mesh % MeshDim
  end function sma_mesh_dim

  subroutine sma_reinit_flag(reinit, found) bind(c, name="sma_reinit_flag")
    integer(c_int) :: reinit, found
    logical :: f, v
    v = ListGetLogical(acc_params, 'Always Detect Structure', f)
    reinit = merge(1, 0, v)
    found = merge(1, 0, f)
  end subroutine sma_reinit_flag

  subroutine sma_reinit_dealloc() bind(c, name="sma_reinit_dealloc")
    if (allocated(acc_field)) deallocate(acc_field)
    if (allocated(acc_surface)) deallocate(acc_surface)
  end subroutine sma_reinit_dealloc

  subroutine sma_recomp_setup(stab, interval) bind(c, name="sma_recomp_setup")
    integer(c_int) :: stab, interval
    logical :: f, f2, v
    integer :: iv
    v = ListGetLogical(acc_params, 'Recompute Stabilization', f)
    if (.not. f) then
      call Info(Caller, 'Defaulting "Recompute Stabilization" to True.', Level=8)
      v = .true.
    end if
    stab = merge(1, 0, v)
    if (v) then
      iv = ListGetInteger(acc_params, 'Recompute Stabilization Interval', f2)
      if (.not. f2) iv = 1
      interval = iv
    else
      interval = 1
    end if
  end subroutine sma_recomp_setup

  subroutine sma_fixed_layers(n, present) bind(c, name="sma_fixed_layers")
    integer(c_int) :: n, present
    logical :: f
    acc_fixed => ListGetIntegerArray(acc_params, 'Fixed Layer Indexes', f)
    if (associated(acc_fixed)) then
      present = 1
      n = size(acc_fixed)
    else
      present = 0
      n = 0
    end if
  end subroutine sma_fixed_layers

  function sma_fixed_get(k) result(v) bind(c, name="sma_fixed_get")
    integer(c_int), value :: k
    integer(c_int) :: v
    v = acc_fixed(k)
  end function sma_fixed_get

  subroutine sma_botproj_flag(val, found) bind(c, name="sma_botproj_flag")
    integer(c_int) :: val, found
    logical :: f, v
    v = ListGetLogical(acc_params, 'Project To Bottom', f)
    val = merge(1, 0, v)
    found = merge(1, 0, f)
  end subroutine sma_botproj_flag

  subroutine sma_ptr_dealloc() bind(c, name="sma_ptr_dealloc")
    if (associated(acc_botpointer)) deallocate(acc_botpointer)
    if (associated(acc_toppointer)) deallocate(acc_toppointer)
    if (associated(acc_uppointer)) deallocate(acc_uppointer)
    if (associated(acc_downpointer)) deallocate(acc_downpointer)
    if (associated(acc_nodelayer)) deallocate(acc_nodelayer)
  end subroutine sma_ptr_dealloc

  ! Multilayer structure detection + fixed-layer enforcement (exact
  ! transcription of the MultiLayer branch of the init block).
  function sma_detect_multi(nfixed) result(nlayers) bind(c, name="sma_detect_multi")
    integer(c_int), value :: nfixed
    integer(c_int) :: nlayers
    integer :: i, nl
    call DetectExtrudedStructure(acc_mesh, acc_solver, ExtVar=acc_var, &
        TopNodePointer=acc_toppointer, BotNodePointer=acc_botpointer, &
        UpNodePointer=acc_uppointer, DownNodePointer=acc_downpointer, &
        NumberOfLayers=nl, NodeLayer=acc_nodelayer)
    nl = nl + 1
    i = acc_fixed(1)
    if (i /= 1) then
      call Warn(Caller, 'Enforcing first fixed layer to: 1 (was '//I2S(i)//')')
      acc_fixed(1) = 1
    end if
    i = acc_fixed(nfixed)
    if (i /= nl) then
      call Warn(Caller, 'Enforcing last fixed layer to: '//I2S(nl)//' (was '//I2S(i)//')')
      acc_fixed(nfixed) = nl
    end if
    nlayers = nl
  end function sma_detect_multi

  ! Binary structure detection (exact transcription of the ELSE branch).
  function sma_detect_binary() result(midex) bind(c, name="sma_detect_binary")
    integer(c_int) :: midex
    logical :: m
    call DetectExtrudedStructure(acc_mesh, acc_solver, ExtVar=acc_var, &
        TopNodePointer=acc_toppointer, BotNodePointer=acc_botpointer, &
        MidNodePointer=acc_midpointer, MidLayerExists=m)
    midex = merge(1, 0, m)
  end function sma_detect_binary

  ! Post-detection sizing, mask setup, OrigCoord allocation (exact
  ! transcription of the tail of the init block).
  subroutine sma_post_detect(nsize, nnodes, mask) bind(c, name="sma_post_detect")
    integer(c_int) :: nsize, nnodes, mask
    integer :: i, istat
    acc_coord => acc_var % Values
    nnodes = acc_mesh % NumberOfNodes
    nsize = min(size(acc_coord), acc_mesh % NumberOfNodes)
    if (associated(acc_var % Perm)) then
      mask = 1
      acc_maskperm => acc_var % Perm
    else
      mask = 0
      if (.not. associated(acc_maskperm)) then
        allocate(acc_maskperm(nsize))
        do i = 1, nsize
          acc_maskperm(i) = i
        end do
      end if
    end if
    if (allocated(acc_origcoord)) deallocate(acc_origcoord)
    allocate(acc_origcoord(nsize), stat=istat)
    if (istat /= 0) then
      call Fatal(Caller, 'Memory allocation error')
    end if
  end subroutine sma_post_detect

  subroutine sma_orig_sync(nsize) bind(c, name="sma_orig_sync")
    integer(c_int), value :: nsize
    acc_origcoord(1:nsize) = acc_coord(1:nsize)
  end subroutine sma_orig_sync

  function sma_cputime() result(t) bind(c, name="sma_cputime")
    real(c_double) :: t
    t = CPUTime()
  end function sma_cputime

  ! -- ValueList string reads ------------------------------------------
  ! which: 1 Base Displacement Variable, 2 Mesh Velocity Variable,
  ! 3 Mesh Update Variable, 4 Top Surface Variable Name,
  ! 5 Bottom Surface Variable Name, 6 Correct Surface Mask.
  function sma_get_string(which, buf, nbuf) result(found) bind(c, name="sma_get_string")
    integer(c_int), value :: which, nbuf
    character(kind=c_char) :: buf(nbuf)
    integer(c_int) :: found
    character(:), allocatable :: s
    logical :: f
    integer :: k
    character(len=64) :: key
    select case (which)
    case (1)
      key = 'Base Displacement Variable'
    case (2)
      key = 'Mesh Velocity Variable'
    case (3)
      key = 'Mesh Update Variable'
    case (4)
      key = 'Top Surface Variable Name'
    case (5)
      key = 'Bottom Surface Variable Name'
    case (6)
      key = 'Correct Surface Mask'
    case default
      key = ''
    end select
    s = GetString(acc_params, trim(key), f)
    found = merge(1, 0, f)
    do k = 1, nbuf
      buf(k) = c_null_char
    end do
    if (f) then
      do k = 1, min(nbuf - 1, len(s))
        buf(k) = s(k:k)
      end do
    end if
  end function sma_get_string

  subroutine sma_base_resolve(cname, namelen) bind(c, name="sma_base_resolve")
    integer(c_int), value :: namelen
    character(kind=c_char) :: cname(namelen)
    character(len=2048) :: fname
    call sma_c2f(cname, namelen, fname)
    acc_basevar => VariableGet(acc_mesh % Variables, trim(fname), ThisOnly=.true.)
    if (.not. associated(acc_basevar)) then
      call Fatal(Caller, 'The variable does not exist: '//trim(fname))
    end if
  end subroutine sma_base_resolve

  subroutine sma_base_first(first, found) bind(c, name="sma_base_first")
    integer(c_int) :: first, found
    logical :: f, v
    v = GetLogical(acc_params, 'Base Displacement First', f)
    first = merge(1, 0, v)
    found = merge(1, 0, f)
    if (v) then
      call Info(Caller, 'Applying base displacement before structural mapping')
    end if
  end subroutine sma_base_first

  subroutine sma_velo_setup(cname, namelen, dim, nsize) bind(c, name="sma_velo_setup")
    integer(c_int), value :: namelen, dim, nsize
    character(kind=c_char) :: cname(namelen)
    character(len=2048) :: fname
    integer :: n
    call sma_c2f(cname, namelen, fname)
    acc_velovar => VariableGet(acc_mesh % Variables, trim(fname))
    if (associated(acc_velovar)) then
      if (acc_velovar % Dofs /= 1) then
        call Fatal(Caller, 'The size of mesh velocity must be one')
      end if
      if (size(acc_velovar % Values) /= nsize) then
        call Fatal(Caller, 'Sizes must agree: '//trim(fname))
      end if
    else
      n = len_trim(fname)
      if (verify(fname(n:n), '123') == 0) then
        call DefaultVariableAdd(fname(1:n-1), DOFs=dim, Perm=acc_maskperm, USolver=acc_solver)
        acc_velovar => VariableGet(acc_mesh % Variables, trim(fname), ThisOnly=.true.)
      else
        call DefaultVariableAdd(trim(fname), Perm=acc_maskperm, Var=acc_velovar, USolver=acc_solver)
      end if
    end if
  end subroutine sma_velo_setup

  subroutine sma_update_setup(cname, namelen, dim, nsize) bind(c, name="sma_update_setup")
    integer(c_int), value :: namelen, dim, nsize
    character(kind=c_char) :: cname(namelen)
    character(len=2048) :: fname
    integer :: n
    call sma_c2f(cname, namelen, fname)
    acc_updatevar => VariableGet(acc_mesh % Variables, trim(fname), ThisOnly=.true.)
    if (associated(acc_updatevar)) then
      if (acc_updatevar % Dofs /= 1) then
        call Fatal(Caller, 'The size of mesh update must be one')
      end if
      if (size(acc_updatevar % Values) /= nsize) then
        call Fatal(Caller, 'Sizes must agree: '//trim(fname))
      end if
    else
      n = len_trim(fname)
      if (verify(fname(n:n), '123') == 0) then
        call DefaultVariableAdd(fname(1:n-1), DOFs=dim, Perm=acc_maskperm, USolver=acc_solver)
        acc_updatevar => VariableGet(acc_mesh % Variables, trim(fname), ThisOnly=.true.)
      else
        call DefaultVariableAdd(trim(fname), Perm=acc_maskperm, Var=acc_updatevar, USolver=acc_solver)
      end if
    end if
  end subroutine sma_update_setup

  subroutine sma_displacement_mode(val, found) bind(c, name="sma_displacement_mode")
    integer(c_int) :: val, found
    logical :: f, v
    v = GetLogical(acc_params, 'Displacement Mode', f)
    val = merge(1, 0, v)
    found = merge(1, 0, f)
  end subroutine sma_displacement_mode

  function sma_minheight() result(mh) bind(c, name="sma_minheight")
    real(c_double) :: mh
    logical :: gotit
    real(kind=dp) :: v
    v = GetCReal(acc_params, 'Minimum Mesh Height', gotit)
    if (.not. gotit) v = GetCReal(acc_params, 'Minimum Height', gotit)
    if (.not. gotit) v = epsilon(v)
    mh = v
  end function sma_minheight

  subroutine sma_info_minheight(mh) bind(c, name="sma_info_minheight")
    real(c_double), value :: mh
    write(Message, '(A,E11.4)') 'Adjusting upper surface to maintain minimum height to:', mh
    call Info(Caller, Message, Level=6)
  end subroutine sma_info_minheight

  ! -- BinaryLayerMapper setup ------------------------------------------
  subroutine sma_detangle_info(detangle, mapheight) bind(c, name="sma_detangle_info")
    integer(c_int) :: detangle, mapheight
    logical :: gotit, v
    v = GetLogical(acc_params, 'Correct Surface', gotit)
    detangle = merge(1, 0, v)
    if (v) then
      call Info(Caller, &
          '> Correct Surface < in case of intersecting upper and lower surface', Level=4)
      if (ListCheckPresent(acc_params, 'Mesh Height Map') ) then
        mapheight = 1
        call Info(Caller, 'Using function to map heights', Level=5)
      else
        mapheight = 0
      end if
    else
      mapheight = 0
    end if
  end subroutine sma_detangle_info

  ! Correct-surface mask setup; returns ComputeTangledMask (the Found flag
  ! of 'Correct Surface Mask', exactly as in the original).
  function sma_tangled_setup() result(compute) bind(c, name="sma_tangled_setup")
    integer(c_int) :: compute
    character(len=2048) :: fname
    logical :: f
    integer :: i
    fname = GetString(acc_params, 'Correct Surface Mask', f)
    compute = merge(1, 0, f)
    if (.not. f) return
    acc_tangledmaskvar => VariableGet(acc_mesh % Variables, trim(fname), ThisOnly=.true.)
    if (.not. associated(acc_tangledmaskvar)) then
      call Info(Caller, &
          'Given > Correct Surface Mask < variable not present, creating it.')
      allocate(acc_tmpperm(acc_mesh % NumberOfNodes))
      do i = 1, acc_mesh % NumberOfNodes
        acc_tmpperm(i) = i
      end do
      call DefaultVariableAdd(trim(fname), Perm=acc_tmpperm, Var=acc_tangledmaskvar, &
          USolver=acc_solver)
      nullify(acc_tmpperm)
    end if
    if (acc_tangledmaskvar % DOFs /= 1) then
      call Fatal(Caller, '> Correct Surface Mask < variable should have only 1 dof')
    end if
    acc_tangledmask => acc_tangledmaskvar % Values
    acc_tangledmask = 1.0_dp
    acc_tangledmaskperm => acc_tangledmaskvar % Perm
    write(Message, '(A,A)')
    call Info(Caller, &
        'Output of > Correct Surface Mask < to: '//trim(fname), Level=6)
  end function sma_tangled_setup

  subroutine sma_top_setup(mode, topval0) bind(c, name="sma_top_setup")
    integer(c_int) :: mode
    real(c_double) :: topval0
    type(variable_t), pointer :: Var
    logical :: gotit
    real(kind=dp) :: v
    character(:), allocatable :: s
    v = GetCReal(acc_params, 'Top Surface Level', gotit)
    if (gotit) then
      mode = 1
      topval0 = v
    else
      s = GetString(acc_params, 'Top Surface Variable Name', gotit)
      if (gotit) then
        Var => VariableGet(acc_mesh % Variables, s, ThisOnly=.true.)
        if (associated(Var)) then
          if (Var % DOFs /= 1) then
            call Fatal(Caller, 'Top surface variable should have only 1 dof')
          else
            acc_topfield => Var % Values
            acc_topperm => Var % Perm
            mode = 2
          end if
          if (InfoActive(20)) then
            print *, 'TopField range:', minval(acc_topfield), maxval(acc_topfield), &
                sum(acc_topfield) / size(acc_topfield)
          end if
        else
          call Fatal(Caller, 'Top surface variable is missing: '//trim(s))
        end if
      else
        mode = 0
        topval0 = 0.0_c_double
      end if
    end if
  end subroutine sma_top_setup

  function sma_top_present() result(p) bind(c, name="sma_top_present")
    integer(c_int) :: p
    p = merge(1, 0, ListCheckPresentAnyBC(acc_model, 'Top Surface'))
  end function sma_top_present

  subroutine sma_fatal_topbc() bind(c, name="sma_fatal_topbc")
    call Fatal(Caller, 'Top surface BC entry is missing')
  end subroutine sma_fatal_topbc

  ! Exact transcription of the TopMode==3 BC gather loop.
  subroutine sma_top_bc(nsize, mask) bind(c, name="sma_top_bc")
    integer(c_int), value :: nsize, mask
    type(Element_t), pointer :: Element
    type(ValueList_t), pointer :: BC
    integer, pointer :: NodeIndexes(:)
    logical :: Found
    integer :: elem, n
    n = acc_mesh % MaxElementNodes
    if (.not. allocated(acc_field)) then
      allocate(acc_field(nsize), acc_surface(n))
      acc_field = 0.0_dp
      acc_surface = 0.0_dp
    end if
    do elem = 1, acc_mesh % NumberOfBoundaryElements
      Element => GetBoundaryElement(elem)
      BC => GetBC()
      if (.not. associated(BC)) cycle
      NodeIndexes => Element % NodeIndexes
      n = GetElementNOFNodes()
      acc_surface(1:n) = GetReal(BC, 'Top Surface', Found)
      if (.not. Found) cycle
      if (mask /= 0) then
        if (all(acc_maskperm(NodeIndexes(1:n)) > 0)) then
          acc_field(acc_maskperm(NodeIndexes(1:n))) = acc_surface(1:n)
        end if
      else
        acc_field(NodeIndexes(1:n)) = acc_surface(1:n)
      end if
    end do
  end subroutine sma_top_bc

  subroutine sma_bot_setup(mode, botval0) bind(c, name="sma_bot_setup")
    integer(c_int) :: mode
    real(c_double) :: botval0
    type(variable_t), pointer :: Var
    logical :: gotit
    real(kind=dp) :: v
    character(:), allocatable :: s
    v = GetCReal(acc_params, 'Bottom Surface Level', gotit)
    if (gotit) then
      mode = 1
      botval0 = v
    else
      s = GetString(acc_params, 'Bottom Surface Variable Name', gotit)
      if (gotit) then
        Var => VariableGet(acc_mesh % Variables, s, ThisOnly=.true.)
        if (associated(Var)) then
          if (Var % DOFs /= 1) then
            call Fatal(Caller, 'Bottom surface variable should have only 1 dof')
          else
            acc_botfield => Var % Values
            acc_botperm => Var % Perm
            mode = 2
          end if
          if (InfoActive(20)) then
            print *, 'BotField range:', minval(acc_botfield), maxval(acc_botfield), &
                sum(acc_botfield) / size(acc_botfield)
          end if
        else
          call Fatal(Caller, 'Bottom surface variable is missing: '//trim(s))
        end if
      else
        mode = 0
        botval0 = 0.0_c_double
      end if
    end if
  end subroutine sma_bot_setup

  function sma_bot_present() result(p) bind(c, name="sma_bot_present")
    integer(c_int) :: p
    p = merge(1, 0, ListCheckPresentAnyBC(acc_model, 'Bottom Surface'))
  end function sma_bot_present

  ! Exact transcription of the BotMode==3 BC gather loop.
  subroutine sma_bot_bc(nsize, mask) bind(c, name="sma_bot_bc")
    integer(c_int), value :: nsize, mask
    type(Element_t), pointer :: Element
    type(ValueList_t), pointer :: BC
    integer, pointer :: NodeIndexes(:)
    logical :: Found
    integer :: elem, n
    n = acc_mesh % MaxElementNodes
    if (.not. allocated(acc_field)) then
      allocate(acc_field(nsize), acc_surface(n))
      acc_field = 0.0_dp
      acc_surface = 0.0_dp
    end if
    do elem = 1, acc_mesh % NumberOfBoundaryElements
      Element => GetBoundaryElement(elem)
      BC => GetBC()
      if (.not. associated(BC)) cycle
      NodeIndexes => Element % NodeIndexes
      n = GetElementNOFNodes()
      acc_surface(1:n) = GetReal(BC, 'Bottom Surface', Found)
      if (.not. Found) cycle
      if (mask /= 0) then
        if (all(acc_maskperm(NodeIndexes(1:n)) > 0)) then
          acc_field(acc_maskperm(NodeIndexes(1:n))) = acc_surface(1:n)
        end if
      else
        acc_field(NodeIndexes(1:n)) = acc_surface(1:n)
      end if
    end do
  end subroutine sma_bot_bc

  ! Exact transcription of the mid-surface BC gather loop.
  subroutine sma_mid_bc(nsize) bind(c, name="sma_mid_bc")
    integer(c_int), value :: nsize
    type(Element_t), pointer :: Element
    type(ValueList_t), pointer :: BC
    integer, pointer :: NodeIndexes(:)
    logical :: Found
    integer :: elem, n
    n = acc_mesh % MaxElementNodes
    if (.not. allocated(acc_field)) then
      allocate(acc_field(nsize), acc_surface(n))
      acc_field = 0.0_dp
      acc_surface = 0.0_dp
    end if
    do elem = 1, acc_mesh % NumberOfBoundaryElements
      Element => GetBoundaryElement(elem)
      BC => GetBC()
      if (associated(BC)) then
        NodeIndexes => Element % NodeIndexes
        n = GetElementNOFNodes()
        acc_surface(1:n) = GetReal(BC, 'Mid Surface', Found)
        if (Found) acc_field(NodeIndexes(1:n)) = acc_surface(1:n)
      end if
    end do
  end subroutine sma_mid_bc

  ! -- per-node data access (binary mapper) ------------------------------
  function sma_maskperm_get(i) result(v) bind(c, name="sma_maskperm_get")
    integer(c_int), value :: i
    integer(c_int) :: v
    v = acc_maskperm(i)
  end function sma_maskperm_get

  function sma_topptr_get(j) result(v) bind(c, name="sma_topptr_get")
    integer(c_int), value :: j
    integer(c_int) :: v
    v = acc_toppointer(j)
  end function sma_topptr_get

  function sma_botptr_get(j) result(v) bind(c, name="sma_botptr_get")
    integer(c_int), value :: j
    integer(c_int) :: v
    v = acc_botpointer(j)
  end function sma_botptr_get

  function sma_midptr_get(j) result(v) bind(c, name="sma_midptr_get")
    integer(c_int), value :: j
    integer(c_int) :: v
    v = acc_midpointer(j)
  end function sma_midptr_get

  function sma_upptr_get(j) result(v) bind(c, name="sma_upptr_get")
    integer(c_int), value :: j
    integer(c_int) :: v
    v = acc_uppointer(j)
  end function sma_upptr_get

  function sma_orig_get(k) result(v) bind(c, name="sma_orig_get")
    integer(c_int), value :: k
    real(c_double) :: v
    v = acc_origcoord(k)
  end function sma_orig_get

  function sma_field_get(k) result(v) bind(c, name="sma_field_get")
    integer(c_int), value :: k
    real(c_double) :: v
    v = acc_field(k)
  end function sma_field_get

  function sma_topperm_get(node) result(v) bind(c, name="sma_topperm_get")
    integer(c_int), value :: node
    integer(c_int) :: v
    v = acc_topperm(node)
  end function sma_topperm_get

  function sma_topfield_at(p) result(v) bind(c, name="sma_topfield_at")
    integer(c_int), value :: p
    real(c_double) :: v
    v = acc_topfield(p)
  end function sma_topfield_at

  function sma_botperm_get(node) result(v) bind(c, name="sma_botperm_get")
    integer(c_int), value :: node
    integer(c_int) :: v
    v = acc_botperm(node)
  end function sma_botperm_get

  function sma_botfield_at(p) result(v) bind(c, name="sma_botfield_at")
    integer(c_int), value :: p
    real(c_double) :: v
    v = acc_botfield(p)
  end function sma_botfield_at

  function sma_coord_get(j) result(v) bind(c, name="sma_coord_get")
    integer(c_int), value :: j
    real(c_double) :: v
    v = acc_coord(j)
  end function sma_coord_get

  subroutine sma_coord_set(j, v) bind(c, name="sma_coord_set")
    integer(c_int), value :: j
    real(c_double), value :: v
    acc_coord(j) = v
  end subroutine sma_coord_set

  subroutine sma_tangled_mark(i) bind(c, name="sma_tangled_mark")
    integer(c_int), value :: i
    acc_tangledmask(acc_tangledmaskperm(i)) = -1.0_dp
  end subroutine sma_tangled_mark

  subroutine sma_warn_tangle_cancel() bind(c, name="sma_warn_tangle_cancel")
    call Warn(Caller, 'Cancelling tanglement when mask exists!')
  end subroutine sma_warn_tangle_cancel

  function sma_fatal_topperm() result(v) bind(c, name="sma_fatal_topperm")
    integer(c_int) :: v
    call Fatal(Caller, 'Top surface variable perm is zero!')
    v = 0
  end function sma_fatal_topperm

  function sma_fatal_botperm() result(v) bind(c, name="sma_fatal_botperm")
    integer(c_int) :: v
    call Fatal(Caller, 'Bottom surface variable perm is zero!')
    v = 0
  end function sma_fatal_botperm

  function sma_listfun(x) result(f) bind(c, name="sma_listfun")
    real(c_double), value :: x
    real(c_double) :: f
    real(kind=dp) :: xx
    xx = x
    f = ListGetFun(acc_params, 'Mesh Height Map', xx)
  end function sma_listfun

  subroutine sma_velo_update(i, v) bind(c, name="sma_velo_update")
    integer(c_int), value :: i
    real(c_double), value :: v
    if (acc_velovar % Perm(i) > 0) acc_velovar % Values(acc_velovar % Perm(i)) = v
  end subroutine sma_velo_update

  subroutine sma_update_update(i, v) bind(c, name="sma_update_update")
    integer(c_int), value :: i
    real(c_double), value :: v
    acc_updatevar % Values(acc_updatevar % Perm(i)) = v
  end subroutine sma_update_update

  ! -- tail of the main body ---------------------------------------------
  function sma_parreduce(v) result(r) bind(c, name="sma_parreduce")
    integer(c_int), value :: v
    integer(c_int) :: r
    r = ParallelReduction(v)
  end function sma_parreduce

  subroutine sma_info_limited(lim, nsize) bind(c, name="sma_info_limited")
    integer(c_int), value :: lim, nsize
    call Info(Caller, 'There seems to be '//I2S(lim)//' (out of '//I2S(nsize)// &
        ') limited heights!', Level=6)
  end subroutine sma_info_limited

  subroutine sma_info_tangled(tan, nsize) bind(c, name="sma_info_tangled")
    integer(c_int), value :: tan, nsize
    call Info(Caller, 'There seems to be '//I2S(tan)//' (out of '//I2S(nsize)// &
        ') tangled nodes!', Level=5)
  end subroutine sma_info_tangled

  function sma_passive(nsize) result(p) bind(c, name="sma_passive")
    integer(c_int), value :: nsize
    integer(c_int) :: p
    logical :: f
    if (ListGetLogical(acc_params, 'Mesh Mapping Passive', f)) then
      call Info(Caller, 'Taking back the suggested mapping!', Level=5)
      acc_coord(1:nsize) = acc_origcoord(1:nsize)
      p = 1
    else
      p = 0
    end if
  end function sma_passive

  function sma_actdir() result(d) bind(c, name="sma_actdir")
    integer(c_int) :: d
    d = ListGetInteger(acc_solver % Values, 'Active Coordinate')
  end function sma_actdir

  subroutine sma_coordp_sel(dir) bind(c, name="sma_coordp_sel")
    integer(c_int), value :: dir
    if (dir == 1) then
      acc_coordp => acc_mesh % Nodes % x
    else if (dir == 2) then
      acc_coordp => acc_mesh % Nodes % y
    else if (dir == 3) then
      acc_coordp => acc_mesh % Nodes % z
    else
      call Fatal(Caller, 'Unknown active coordinate!')
    end if
  end subroutine sma_coordp_sel

  subroutine sma_coordp_set(i, v) bind(c, name="sma_coordp_set")
    integer(c_int), value :: i
    real(c_double), value :: v
    acc_coordp(i) = v
  end subroutine sma_coordp_set

  function sma_velo_firstzero() result(z) bind(c, name="sma_velo_firstzero")
    integer(c_int) :: z
    logical :: f
    z = merge(1, 0, GetLogical(acc_params, 'Mesh Velocity First Zero', f))
  end function sma_velo_firstzero

  subroutine sma_velo_zero() bind(c, name="sma_velo_zero")
    acc_velovar % Values = 0.0_dp
  end subroutine sma_velo_zero

  subroutine sma_info_maptime(at0) bind(c, name="sma_info_maptime")
    real(c_double), value :: at0
    real(kind=dp) :: at1
    at1 = CPUTime()
    write(Message, *) 'Active coordinate mapping time: ', at1 - at0
    call Info(Caller, Message)
  end subroutine sma_info_maptime

  subroutine sma_maybe_write() bind(c, name="sma_maybe_write")
    character(len=2048) :: mname
    logical :: w
    mname = GetString(acc_params, 'Mapped Mesh Name', w)
    if (w) then
      call WriteMeshToDisk(acc_mesh, trim(mname))
    end if
  end subroutine sma_maybe_write

  subroutine sma_meshstab() bind(c, name="sma_meshstab")
    call MeshStabParams(acc_mesh)
  end subroutine sma_meshstab

  ! -- MultiLayerMapper ----------------------------------------------------
  subroutine sma_info_multilay_using(nfixed) bind(c, name="sma_info_multilay_using")
    integer(c_int), value :: nfixed
    call Info(Caller, 'Mapping using '//I2S(nfixed)//' fixed layers', Level=6)
  end subroutine sma_info_multilay_using

  subroutine sma_fatal_multimask() bind(c, name="sma_fatal_multimask")
    call Fatal(Caller, 'Mask not available yet for multiple layers!')
  end subroutine sma_fatal_multimask

  function sma_detangle_flag() result(d) bind(c, name="sma_detangle_flag")
    integer(c_int) :: d
    logical :: gotit
    d = merge(1, 0, GetLogical(acc_params, 'Correct Surface', gotit))
  end function sma_detangle_flag

  ! Fixed Layer Variable resolve + DOFs check (name via UnfoundFatal,
  ! exactly as the original).
  subroutine sma_fixedvar_setup(nfixed) bind(c, name="sma_fixedvar_setup")
    integer(c_int), value :: nfixed
    character(:), allocatable :: s
    s = ListGetString(acc_params, 'Fixed Layer Variable', UnfoundFatal=.true.)
    acc_fixedvar => VariableGet(acc_mesh % Variables, s)
    if (.not. associated(acc_fixedvar)) then
      call Fatal(Caller, 'Could not find variable: '//trim(s))
    end if
    if (acc_fixedvar % Dofs /= nfixed) then
      call Fatal(Caller, 'Invalid number of components in fixed layer variable:'// &
          I2S(acc_fixedvar % Dofs))
    end if
  end subroutine sma_fixedvar_setup

  function sma_fixedvar_perm(node) result(p) bind(c, name="sma_fixedvar_perm")
    integer(c_int), value :: node
    integer(c_int) :: p
    p = acc_fixedvar % Perm(node)
  end function sma_fixedvar_perm

  subroutine sma_fixedvar_vals(k, nf, buf) bind(c, name="sma_fixedvar_vals")
    integer(c_int), value :: k, nf
    real(c_double) :: buf(nf)
    buf(1:nf) = acc_fixedvar % Values(nf*(k-1)+1:nf*k)
  end subroutine sma_fixedvar_vals

  subroutine sma_fatal_nofield() bind(c, name="sma_fatal_nofield")
    call Fatal(Caller, 'Could not find fixed fields at top or bottom!')
  end subroutine sma_fatal_nofield

  function sma_fatal_nolayer(k) result(v) bind(c, name="sma_fatal_nolayer")
    integer(c_int), value :: k
    integer(c_int) :: v
    call Fatal(Caller, 'Could not find mapping for layer: '//I2S(k))
    v = 0
  end function sma_fatal_nolayer

  ! Order-sensitive projector application (MATMUL stays in Fortran).
  ! proj is column-major (nl rows), matching the dummy proj(nl,nf).
  subroutine sma_stride_project(proj, fixed, stride, nl, nf, displ) &
      bind(c, name="sma_stride_project")
    integer(c_int), value :: nl, nf, displ
    real(c_double) :: proj(nl, nf), fixed(nf), stride(nl)
    if (displ /= 0) then
      stride(1:nl) = stride(1:nl) + matmul(proj(1:nl, 1:nf), fixed(1:nf))
    else
      stride(1:nl) = matmul(proj(1:nl, 1:nf), fixed(1:nf))
    end if
  end subroutine sma_stride_project

  ! Indexed whole-array commit + conditional velocity/update blocks
  ! (exact transcription of the stride commit and WHERE blocks).
  subroutine sma_stride_commit(ind, vals, n) bind(c, name="sma_stride_commit")
    integer(c_int), value :: n
    integer(c_int) :: ind(n)
    real(c_double) :: vals(n)
    integer :: k
    do k = 1, n
      acc_coord(ind(k)) = vals(k)
    end do
  end subroutine sma_stride_commit

  subroutine sma_stride_velo(ind, vals, n, dt) bind(c, name="sma_stride_velo")
    integer(c_int), value :: n
    integer(c_int) :: ind(n)
    real(c_double) :: vals(n)
    real(c_double), value :: dt
    integer :: k, p
    do k = 1, n
      p = acc_velovar % Perm(ind(k))
      if (p > 0) acc_velovar % Values(p) = (vals(k) - acc_origcoord(ind(k))) / dt
    end do
  end subroutine sma_stride_velo

  subroutine sma_stride_update(ind, vals, n) bind(c, name="sma_stride_update")
    integer(c_int), value :: n
    integer(c_int) :: ind(n)
    real(c_double) :: vals(n)
    integer :: k
    do k = 1, n
      acc_updatevar % Values(acc_updatevar % Perm(ind(k))) = &
          vals(k) - acc_origcoord(ind(k))
    end do
  end subroutine sma_stride_update

  subroutine sma_info_multilay_done() bind(c, name="sma_info_multilay_done")
    call Info(Caller, 'Finished multilayer mapping', Level=8)
  end subroutine sma_info_multilay_done

  ! -- BaseVarDisplace -------------------------------------------------------
  subroutine sma_info_base_adding() bind(c, name="sma_info_base_adding")
    call Info(Caller, 'Adding base displacement to the displacements!')
  end subroutine sma_info_base_adding

  function sma_basedofs() result(d) bind(c, name="sma_basedofs")
    integer(c_int) :: d
    integer :: dofs
    dofs = acc_basevar % Dofs
    if (dofs /= 2 .and. dofs /= 3) then
      call Fatal(Caller, 'Invalid base displacement dimension: '//I2S(dofs))
    end if
    d = dofs
  end function sma_basedofs

  subroutine sma_basevals(ibot, dofs, bx, by, bz) bind(c, name="sma_basevals")
    integer(c_int), value :: ibot, dofs
    real(c_double) :: bx, by, bz
    integer :: pb
    pb = acc_basevar % Perm(ibot)
    bx = acc_basevar % Values(dofs*(pb-1)+1)
    by = acc_basevar % Values(dofs*(pb-1)+2)
    bz = 0.0_c_double
    if (dofs == 3) then
      bz = acc_basevar % Values(dofs*pb)
    end if
  end subroutine sma_basevals

  subroutine sma_node_add(i, dx, dy, dz, dofs) bind(c, name="sma_node_add")
    integer(c_int), value :: i, dofs
    real(c_double), value :: dx, dy, dz
    acc_mesh % Nodes % x(i) = acc_mesh % Nodes % x(i) + dx
    acc_mesh % Nodes % y(i) = acc_mesh % Nodes % y(i) + dy
    if (dofs == 3) then
      acc_mesh % Nodes % z(i) = acc_mesh % Nodes % z(i) + dz
    end if
  end subroutine sma_node_add

  subroutine sma_fatal_novar(cname, namelen) bind(c, name="sma_fatal_novar")
    integer(c_int), value :: namelen
    character(kind=c_char) :: cname(namelen)
    character(len=2048) :: fname
    call sma_c2f(cname, namelen, fname)
    call Fatal(Caller, 'The variable does not exist: '//trim(fname))
  end subroutine sma_fatal_novar

end module structuredmeshmapper_acc
