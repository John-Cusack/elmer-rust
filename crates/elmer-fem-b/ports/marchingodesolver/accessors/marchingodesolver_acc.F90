! Port of fem/src/modules/MarchingODESolver.F90 (Elmer fork) for lane T021-B,
! at the solver-plugin boundary (whole `.so` -> Rust cdylib).
!
! Statement-level BIND(C) wrappers for the Elmer calls in MarchingODESolver.F90,
! compiled against the built .mod files of the target configure. Every
! wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) out-args). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls.
!
! The four CONTAINS internal procedures (FirstExtrudedMaterial,
! SetLayerValues, GetLayerValues, GetCoefficients) are transcribed verbatim
! as module procedures over mirrored module state (same compiler, same
! statement order): component access, Info/Fatal I/O, I2S string building,
! list-directed PRINT, bare STOP, and the whole-array updates with their
! SUM reductions all stay in Fortran, so bit-identity holds by construction
! wherever the call sequence matches. The entry-level control flow
! (initialization dispatch, Beta branch, transient dtn logic, layer /
! substep / nonlinear-iteration loops, the GOTO-1 marching-step loop) is
! driven from Rust through the moa_* wrappers below, each a verbatim
! statement block of the original body.
!
! SAVE semantics: every name in the original's SAVE list (plus the =init
! implied-SAVE locals MaskExist, ParabolicModel, RequireBC, DoTransient,
! AnyDG, VectorSource) lives here with identical lifetime — set once at
! initialization, never reset per call. The GetCoefficients SAVE state
! (Visited, Source_h, xloc, cdofs) likewise persists across invocations.
! Per-call locals (Model context, Params, Mesh, Var3D, Material, dofs,
! n, m, layer, subt, Norm, Change, dtime, HaveF/HaveC/HaveR) are module
! variables re-set on every entry call, exactly as the original
! re-executes those statements.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module marchingodesolver_acc
  use iso_c_binding
  use types
  use defutils
  use messages
  use meshtransform, only : DetectExtrudedStructure
  use meshbasics, only : MarkBCNodes
  implicit none

  character(*), parameter :: Caller = 'MarchingODESolver'

  ! -- mirrored SAVE state: the original's SAVE list --
  logical, save :: Initialized = .false.
  integer, pointer, save :: BotPointer(:) => null(), UpPointer(:) => null(), &
      BotPerm(:) => null(), InvPerm(:) => null(), PrevInvPerm(:) => null(), &
      MaskPerm(:) => null(), SingleIndex(:) => null()
  integer, allocatable, save :: ParentElem(:)
  type(variable_t), pointer, save :: ExtVar => null()
  integer, save :: NumberOfLayers = 0, BotNodes = 0
  logical, save :: MaskExist = .false.
  character(len=MAX_NAME_LEN), save :: TimeMethod = ''
  logical, save :: RequireBC = .false., UseInternalVals = .false., &
      AnyDG = .false.
  type(variable_t), pointer, save :: AddVar => null()
  real(kind=dp), save :: HParam = 0.0_dp, mincons = 0.0_dp
  logical, save :: SetMin = .false., SetSum = .false.
  real(kind=dp), pointer, save :: Coord(:) => null(), xvec(:) => null(), &
      xivec(:) => null(), dxvec(:) => null(), x0vec(:) => null(), &
      fvec(:) => null(), rvec(:) => null(), cvec(:) => null(), &
      f0vec(:) => null(), r0vec(:) => null(), c0vec(:) => null()
  ! -- implied SAVE via =init at declaration (Fortran SAVE semantics) --
  logical, save :: ParabolicModel = .false., VectorSource = .false., &
      DoTransient = .false.
  ! -- GetCoefficients internal-procedure SAVE state --
  logical, save :: Visited = .false.
  type(ValueHandle_t), pointer, save :: Source_h(:) => null()
  real(kind=dp), allocatable, save :: xloc(:)
  integer, save :: cdofs = 0

  ! -- per-call entry context (mirrors the entry's un-saved locals) --
  type(model_t), pointer, save :: Model => null()
  type(solver_t), pointer, save :: PSolver => null()
  type(valuelist_t), pointer, save :: Params => null()
  type(mesh_t), pointer, save :: Mesh => null()
  type(variable_t), pointer, save :: Var3D => null()
  type(valuelist_t), pointer, save :: Material => null()
  character(len=MAX_NAME_LEN), save :: VarName = ''
  integer, save :: dofs = 0, n = 0, m = 0
  integer, save :: layer = 0, subt = 0, maxsubt = 1
  real(kind=dp), save :: velo = 0.0_dp, dtime = 0.0_dp
  real(kind=dp), save :: Norm = 0.0_dp, Change = 0.0_dp
  logical, save :: HaveF = .false., HaveC = .false., HaveR = .false.

contains

  ! -- MarchingODESolver_init body (verbatim) ---------------------------
  subroutine moa_init_params() bind(c, name="moa_init_params")
    type(valuelist_t), pointer :: lparams
    lparams => GetSolverParams()
    ! By construction time derivative order must be one!
    call ListAddNewInteger( lparams, 'Time derivative order', 1 )
    ! Global mass matrix allows us to use library routines for time integration
    call ListAddNewLogical( lparams, 'Use Global Mass Matrix', .TRUE. )
    ! We initialize our own matrix structures
    call ListAddNewLogical( lparams, 'No Matrix', .TRUE. )
  end subroutine moa_init_params

  ! -- per-call prologue: Info headers, Params/Mesh/Var3D/dofs (verbatim) --
  subroutine moa_begin(h_model, h_solver) bind(c, name="moa_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, Model)
    call c_f_pointer(h_solver, PSolver)
    call Info(Caller,'-----------------------------------------------------',Level=6)
    call Info(Caller,'Solving ODE on moving coordinates in structured mesh',Level=4)
    Params => GetSolverParams()
    Mesh => PSolver % Mesh
    ! In principle we could use the same solver to advect many fields on
    ! same extruded mesh. Hence these are not saved.
    Var3D => PSolver % Variable
    if (.not. associated(Var3D)) then
      call Fatal(Caller,'Solver variable not associated!')
    end if
    VarName = trim( Var3D % Name )
    call Info(Caller,'Working with variable: '//trim(VarName),Level=7)
    dofs = Var3D % dofs
    call Info(Caller,'Number of components in variable: '//I2S(dofs),Level=7)
    if (dofs < 1) call Fatal(Caller,'Invalid number of components!')
  end subroutine moa_begin

  function moa_initialized() result(r) bind(c, name="moa_initialized")
    integer(c_int) :: r
    r = merge(1, 0, Initialized)
  end function moa_initialized

  ! -- IF (.NOT. Initialized) block, transcribed verbatim ----------------
  subroutine moa_initialize() bind(c, name="moa_initialize")
    integer :: NoActive, t, i, j, NoBCNodes
    type(element_t), pointer :: Element => null()
    logical, allocatable :: BCNode(:)
    character(len=MAX_NAME_LEN) :: str
    logical :: Found
    real(kind=dp) :: dz

    ! Just check what type of elements we have
    NoActive = GetNOFActive()
    do t=1,NoActive
      Element => GetActiveElement(t)
      i = Element % TYPE % ElementCode
      if (t == 1) then
        j = i
      else if (i /= j) then
        exit
      end if
    end do
    if ( i /= j ) then
      call Info(Caller,'There are at least elements of type: '//I2S(i)//' and '//I2S(j))
    else
      call Info(Caller,'All elements are of type: '//I2S(i),Level=12)
    end if

    call Info(Caller,'Initializing structured mesh and ODE structures',Level=6)
    call Info(Caller,'Solving for variable: '//trim(VarName),Level=6)

    call DetectExtrudedStructure( Mesh, PSolver, ExtVar = ExtVar, &
        BotNodePointer = BotPointer, UpNodePointer = UpPointer, &
        NumberOfLayers = NumberOfLayers, MaskVar = Var3D )

    call Info(Caller,'Number of element layers: '//I2S(NumberOfLayers),Level=7)

    MaskExist = associated( ExtVar % Perm )
    if ( MaskExist ) then
      call Info(Caller,'We have a mask of size:'//I2S(maxval(ExtVar % Perm)),Level=7)
    else
      call Info(Caller,'No mask associated to solver',Level=20)
    end if

    if ( MaskExist ) MaskPerm => ExtVar % Perm
    Coord => ExtVar % Values

    if (.not. associated(Coord) ) then
      call Fatal(Caller,'Coord is not associated!')
    end if

    write(Message,'(A,2ES12.3)') 'Extruded Coordinate Range:',minval(Coord),maxval(Coord)
    call Info(Caller,Message,Level=8)

    dz = maxval(Coord)-minval(Coord)
    if ( dz < epsilon(dz) ) then
      call Fatal(Caller,'We cannot march when dz is ~zero!')
    end if

    AnyDG = ListGetLogicalAnySolver( Model,'Discontinuous Galerkin')

    UseInternalVals = ListGetLogical( Params,'Use Internal Values', Found )
    str = ListGetString( Params,'Additional Internal Variable', Found )
    if (Found) then
      AddVar => VariableGet( Mesh % Variables, str, UnfoundFatal = .TRUE. )
      call Info(Caller,'Using additional internal variable: '//trim(AddVar % Name))
      if (InfoActive(20)) then
        call VectorValuesRange(AddVar % Values,size(AddVar % Values),'AddVar')
      end if
    end if

    mincons = ListGetCReal( Params,'Minimum Cons',SetMin )
    SetSum = ListGetLogical( Params,'Enforce Unity Sum',Found )

    ! It is not trivial to know to which element a node belongs to.
    ! This structure is needed when we want to know the DG field value of a given node.
    ! Only if we also have "Discontinuous Bodies" within this active domain will this be
    ! uniquely defined.
    !-----------------------------------------------------------------------------------
    if ( MaskExist .or. AnyDG .or. UseInternalVals ) then
      call Info(Caller,'Creating inverse node parent look-up table',Level=7)
      allocate( ParentElem(Mesh % NumberOfNodes) )
      ParentElem = 0

      NoActive = GetNOFActive()
      do t=1,NoActive
        Element => GetActiveElement(t)
        n = Element % TYPE % NumberOfNodes
        if ( MaskExist ) then
          if ( any( MaskPerm(Element % NodeIndexes) == 0 ) ) cycle
        end if
        do i=1,n
          j = Element % NodeIndexes(i)
          if ( ParentElem(j) == 0 ) ParentElem(j) = Element % ElementIndex
        end do
      end do
    end if

    Hparam = ( maxval( Coord ) - minval( Coord ) ) / NumberOfLayers
    write(Message,'(A,ES12.3)') 'Constant mesh parameter: ',Hparam
    call Info(Caller,Message,Level=7)

    ! We may choose only to apply the ODE to BC nodes
    RequireBC = ListGetLogical( Params,'Apply BCs Only',Found )
    if ( RequireBC ) then
      call MarkBCNodes( Mesh,BCNode,NoBCNodes)
      if (NoBCNodes == 0 ) RequireBC = .FALSE.
      call Info(Caller,'Number of BC nodes: '//I2S(NoBCNodes),Level=7)
    end if

    ! Create the permutation using the bottom layer
    BotNodes = 0
    allocate( BotPerm( Mesh % NumberOfNodes ) )
    BotPerm = 0
    do i=1,Mesh % NumberOfNodes
      ! The variable to be marched does not exist at the 1st layer
      if ( Var3D % Perm(i) == 0 ) cycle

      if ( RequireBC ) then
        if ( .not. BcNode(i) ) cycle
      end if

      j = i
      if ( MaskExist ) then
        j = MaskPerm(i)
        if ( j == 0 ) cycle
      end if

      ! This is not at the bottom
      if (BotPointer(j) /= i) cycle

      ! Ok, check that also the next node would be at BC
      ! The 1st layer may be an exception, check for 2nd too

      ! The variable to be marched does not exist
      if ( Var3D % Perm(UpPointer(j)) == 0 ) cycle

      ! We want BC node but don't have one
      if ( RequireBC ) then
        if ( .not. BcNode(UpPointer(j)) ) cycle
      end if

      BotNodes = BotNodes + 1
      BotPerm(i) = BotNodes
    end do
    n = BotNodes
    allocate( InvPerm(n), PrevInvPerm(n), SingleIndex(1) )

    call Info(Caller,'Number of bottom nodes: '//I2S(n),Level=7)

    ! Allocate some vectors to study convergence
    m = dofs * n
    allocate( xvec(m), fvec(m), rvec(m), cvec(m), f0vec(m), r0vec(m), &
        c0vec(m), xivec(m), dxvec(m), x0vec(m) )

    if ( RequireBC ) deallocate(BcNode)

    Initialized = .TRUE.
    call Info(Caller,'Initialization done',Level=10)
  end subroutine moa_initialize

  ! -- The variable on the layer: n = BotNodes; m = dofs * n (verbatim) --
  subroutine moa_set_dims() bind(c, name="moa_set_dims")
    n = BotNodes
    m = dofs * n
  end subroutine moa_set_dims

  ! -- Beta / TimeMethod branch (verbatim) --------------------------------
  function moa_beta() result(Beta) bind(c, name="moa_beta")
    real(c_double) :: Beta
    logical :: Found
    ! We just use one parameter to define the timestepping.
    ! This defines how the coefficients are to be evaluated.
    Beta = ListGetCReal( Params,'Newmark Beta',Found )
    if (.not. Found ) then
      ! Default timestepping is impicit euler
      Beta = 1.0_dp
      TimeMethod = ListGetString( Params, 'Timestepping Method',Found )
      if ( Found ) then
        if ( TimeMethod == 'implicit euler' ) then
          Beta = 1.0_dp
        else if ( TimeMethod == 'explicit euler' ) then
          Beta = 0.0_dp
        else if ( TimeMethod == 'crank-nicolson' ) then
          Beta = 0.5_dp
        else if ( TimeMethod == 'newmark' ) then
          Beta = ListGetCReal( Params,'Newmark Beta',UnfoundFatal=.TRUE. )
        end if
      end if
    end if
  end function moa_beta

  ! -- ParabolicModel flag (verbatim) --------------------------------------
  subroutine moa_parabolic() bind(c, name="moa_parabolic")
    logical :: Found
    ParabolicModel = ListGetLogical( Params,'Parabolic Model',Found )
    if ( ParabolicModel ) then
      call Info(Caller,'Using parabolic growth model',Level=7)
    end if
  end subroutine moa_parabolic

  ! -- Draw Velocity (verbatim) ---------------------------------------------
  subroutine moa_velo() bind(c, name="moa_velo")
    logical :: Found
    velo = GetCReal( Model % Simulation,'Draw Velocity',Found )
    if (.not. Found ) velo = GetCReal( Params,'Draw Velocity',Found )
    if (.not. Found ) then
      call Fatal(Caller,'>Draw Velocity< is needed for marching solver!')
    end if
  end subroutine moa_velo

  ! -- Nonlinear iteration count (verbatim) ---------------------------------
  function moa_maxiter() result(MaxIter) bind(c, name="moa_maxiter")
    integer(c_int) :: MaxIter
    logical :: Found
    MaxIter = GetInteger( Params,'Nonlinear System Max Iterations',Found )
    if (.not. Found) MaxIter = 1
    call Info(Caller,'Number of nonlinear iterations set to: '//I2S(MaxIter),Level=7)
  end function moa_maxiter

  ! -- Nonlinear convergence tolerance (verbatim) ------------------------------
  function moa_nonlintol() result(NonLinTol) bind(c, name="moa_nonlintol")
    real(c_double) :: NonLinTol
    logical :: Found
    NonLinTol = GetCReal( Params,'Nonlinear System Convergence Tolerance',Found )
  end function moa_nonlintol

  ! -- Material, zeroed vectors, source flags (verbatim) ------------------------
  subroutine moa_material_setup() bind(c, name="moa_material_setup")
    logical :: Found
    Material => FirstExtrudedMaterial()

    xvec = 0.0_dp; fvec = 0.0_dp; rvec = 0.0_dp; cvec = 0.0_dp
    x0vec = 0.0_dp; f0vec = 0.0_dp; r0vec = 0.0_dp; c0vec = 0.0_dp
    xivec = 0.0_dp; dxvec = 0.0_dp

    HaveF = ListCheckPresent( Material,trim(VarName)//': Source')
    HaveR = ListCheckPresent( Material,trim(VarName)//': Reaction Coefficient')
    HaveC = ListCheckPresent( Material,trim(VarName)//': Time Derivative Coefficient')

    VectorSource = .FALSE.
    if (dofs > 1 ) then
      if (HaveF .or. HaveR .or. HaveC ) then
        VectorSource = .TRUE.
      else
        HaveF = ListCheckSuffix( Material,': Source')
        HaveR = ListCheckSuffix( Material,': Reaction Coefficient')
        HaveC = ListCheckSuffix( Material,': Time Derivative Coefficient')
      end if
    end if

    if ( HaveR ) then
      call Fatal(Caller,'Code some more to account for reaction term!')
    else if ( .not. HaveC ) then
      call Info(Caller,'By default "Time Derivative Coefficient" one will be used',Level=7)
    end if
  end subroutine moa_material_setup

  ! -- Transient dtn logic (verbatim) ----------------------------------------------
  function moa_transient(dt_in, trans_in) result(dtn) bind(c, name="moa_transient")
    real(c_double), value :: dt_in
    integer(c_int), value :: trans_in
    integer(c_int) :: dtn
    logical :: trans
    real(kind=dp) :: dth
    trans = (trans_in /= 0)
    dtn = 0
    if ( trans ) then
      ! timestep defined by velocity & mesh parameter
      dth = Hparam / velo
      dtn = nint(dt_in / dth )
      if ( dtn >= NumberOfLayers ) then
        call Info(Caller,'Timestep so large than we can use steady algo!')
        dtn = 0
      else
        if ( abs( dt_in/dth - dtn ) > 0.01 ) then
          print *,'Mesh parameter:',Hparam
          print *,'Draw Velocity:',velo
          print *,'Suggested timesteps:',dt_in,dth
          call Fatal(Caller,'Timesteps are not matching')
        else
          call Info(Caller,'Number of marching steps for each timestep: '//I2S(dtn),Level=5)
        end if
      end if
    end if
    DoTransient = ( dtn > 0 )
  end function moa_transient

  ! -- Passage head at label 1: rebuild InvPerm, read layer, Set0 (verbatim) --
  function moa_passage_begin() result(msub) bind(c, name="moa_passage_begin")
    integer(c_int) :: msub
    integer :: i, j
    logical :: Found
    ! First layer (= 0) is determined by the initial conditions (=boundary conditions)
    !-----------------------------------------------------------------------------------
    do i=1,Mesh % NumberOfNodes
      j = BotPerm(i)
      if ( j > 0 ) then
        InvPerm(j) = i
      end if
    end do
    PrevInvPerm = InvPerm

    if ( any(InvPerm == 0 ) ) then
      call Fatal(Caller,'Number of nodes has InvPerm undefined: '//I2S(count(InvPerm==0)))
    end if

    call GetLayerValues( xvec )

    maxsubt = ListGetInteger( Params,'Timestep Divisions',Found)
    if (.not. Found) maxsubt = 1
    if ( maxsubt > 1 .and. .not. UseInternalVals ) then
      call Fatal(caller,'We can only have substeps when we use internal values!')
    end if

    ! 0-values at values at the previous layer
    ! The 1st layer cannot really change since it is the BC.
    layer = 0
    call GetCoefficients(Set0=.TRUE., xlayer = xvec)
    x0vec = xvec
    msub = maxsubt
  end function moa_passage_begin

  function moa_nlayers() result(r) bind(c, name="moa_nlayers")
    integer(c_int) :: r
    r = NumberOfLayers
  end function moa_nlayers

  ! -- Advance to the next layer: InvPerm, dz/dtime, GetLayerValues (verbatim) --
  subroutine moa_advance_layer(layer_in, dtime_out) bind(c, name="moa_advance_layer")
    integer(c_int), value :: layer_in
    real(c_double) :: dtime_out
    real(kind=dp) :: dz
    layer = layer_in
    call Info(Caller,'Solving for layer: '//I2S(layer),Level=12)
    ! Find the next level of nodes, and remember the previous one.
    PrevInvPerm = InvPerm
    if ( MaskExist ) then
      InvPerm = UpPointer(MaskPerm(PrevInvPerm))
    else
      InvPerm = UpPointer(PrevInvPerm)
    end if
    if ( InvPerm(1) == PrevInvPerm(1) ) then
      call Fatal(Caller,'InvPerm is the same on different layers: '//I2S(InvPerm(1)))
    end if
    ! This sets the timestep assuming that all nodes are extruded equally.
    ! Hence this only applied to cartesian drawing.
    if ( MaskExist ) then
      dz = Coord(MaskPerm(InvPerm(1))) - Coord(MaskPerm(PrevInvPerm(1)))
    else
      dz = Coord(InvPerm(1)) - Coord(PrevInvPerm(1))
    end if
    dtime = dz / velo
    write(Message,'(A,2ES12.3)') 'Layer thickness and timestep: ',dz,dtime
    call Info(Caller,Message,Level=12)
    if (dtime < epsilon(dtime) ) then
      call Fatal(Caller,'Cannot march if timestep is ~zero')
    end if
    call GetLayerValues( xivec )
    dtime_out = dtime
  end subroutine moa_advance_layer


  ! -- Sub-timestep head (verbatim) ----------------------------------------------
  subroutine moa_sub_head(subt_in) bind(c, name="moa_sub_head")
    integer(c_int), value :: subt_in
    subt = subt_in
    if ( maxsubt > 1 ) then
      call Info(Caller,'Solving for sub timestep: '//I2S(subt),Level=15)
    end if
  end subroutine moa_sub_head

  ! -- Nonlinear-iteration head (verbatim) -------------------------------------------
  subroutine moa_iter_head(iter_in, maxiter_in) bind(c, name="moa_iter_head")
    integer(c_int), value :: iter_in, maxiter_in
    ! We may have iteration if the ODE is nonlinear.
    ! However, more typical could be to iterate over coupled systems.
    if ( maxiter_in > 1 ) then
      call Info(Caller,'Nonlinear iteration: '//I2S(iter_in),Level=20)
      xivec = xvec
    end if
  end subroutine moa_iter_head

  ! -- One nonlinear iteration: coefficients, update, norms (verbatim) --
  subroutine moa_solve_iter(q_beta, change_out) bind(c, name="moa_solve_iter")
    real(c_double), value :: q_beta
    real(c_double) :: change_out
    integer :: i
    real(kind=dp) :: q, sumcons, corr
    q = q_beta
    call GetCoefficients(q, xlayer = xvec )
    ! We don't really need any linear solver for this as there is no coupling among dofs.
    if ( HaveC ) then
      xvec = x0vec + ( dtime / maxsubt ) * fvec / cvec
    else
      xvec = x0vec + ( dtime / maxsubt ) * fvec
    end if
    if (SetMin) then
      xvec = max(xvec,mincons)
    end if
    if (SetSum) then
      do i=1,n
        sumcons = sum(xvec(dofs*(i-1)+1:dofs*i))
        ! The scaling get a little bit more complex when we want to maintain the minimum cuts.
        if (SetMin) then
          corr = (1-dofs*mincons) / (sumcons - dofs*mincons)
          xvec(dofs*(i-1)+1:dofs*i) = corr * ( xvec(dofs*(i-1)+1:dofs*i) - mincons) + mincons
        else
          xvec(dofs*(i-1)+1:dofs*i) = xvec(dofs*(i-1)+1:dofs*i) / sumcons
        end if
      end do
    end if
    ! For the 1st iteration the error corresponds to the error with respect to previous solution.
    ! For 2nd etc. iteration the error is of the nonlinear iteration.
    Norm = sqrt(sum(xvec*xvec))
    dxvec = xvec-xivec
    Change = sqrt(sum(dxvec*dxvec)) / Norm
    change_out = Change
    ! This must be in the loop since we may have dependence on some field value that has changed!
    if (.not. UseInternalVals ) then
      call SetLayerValues( xvec )
    end if
  end subroutine moa_solve_iter

  ! -- Sub-timestep tail: internal-values write, PRINT, x0 advance (verbatim) --
  subroutine moa_sub_tail() bind(c, name="moa_sub_tail")
    if (UseInternalVals .and. subt == maxsubt ) then
      call SetLayerValues( xvec )
    end if
    if ( InfoActive(20) ) then
      print *,'Layer:',layer,dtime,Norm,Change
    end if
    ! 0-values at values at the previous layer
    if ( layer < NumberOfLayers .or. subt < maxsubt ) then
      if ( DoTransient ) then
        ! The only way to have transient simulation is that the timestep is such
        ! that we take exactly one extruded layer. Then the initial value is the
        ! next starting value of the next layer. We have to do some back-and-forth
        ! stuff to have true previous timestep starting values for the coefficients.
        if (.not. UseInternalVals ) then
          call SetLayerValues( xivec )
        end if
        ! Set the coefficients at the start at the timestep, hence Set0 = .TRUE.
        call GetCoefficients(Set0=.TRUE.,xlayer=xivec)
        x0vec = xivec
        if (.not. UseInternalVals ) then
          call SetLayerValues( xvec )
        end if
      else
        ! For steady state the initial value is the final value of this layer.
        call GetCoefficients(Set0=.TRUE.,xlayer=xvec)
        x0vec = xvec
      end if
    end if
  end subroutine moa_sub_tail

  subroutine moa_info_step(dti_in) bind(c, name="moa_info_step")
    integer(c_int), value :: dti_in
    call Info(Caller,'Taking marching step: '//I2S(dti_in))
  end subroutine moa_info_step

  subroutine moa_info_done() bind(c, name="moa_info_done")
    call Info(Caller,'All done',Level=5)
    call Info(Caller,'-----------------------------------------------------',Level=6)
  end subroutine moa_info_done

  ! Finds pointer to the extruded material.
  ! Note that we assume that this is the only material being extruded.
  !--------------------------------------------------------------------------------
  function FirstExtrudedMaterial() result ( Material0 )
    type(valuelist_t), pointer :: Material0
    type(element_t), pointer :: Element => null()
    Material0 => null()
    Element => GetActiveElement(1)
    if ( associated( Element ) ) then
      Material0 => GetMaterial(Element)
    end if
    if (.not. associated( Material0 ) ) then
      call Fatal(Caller,'Could not set material for extruded domain!')
    end if
  end function FirstExtrudedMaterial

  ! Back-substitute values on one single layer back to the distributed field on the
  ! finite element mesh
  !-----------------------------------------------------------------------------------
  subroutine SetLayerValues( xlayer )
    real(kind=dp) :: xlayer(:)
    integer :: j
    if (dofs == 1 ) then
      if ( ParabolicModel ) then
        Var3D % Values(Var3D % Perm(InvPerm)) = sqrt( 2 * xlayer(1:m) )
      else
        Var3D % Values(Var3D % Perm(InvPerm)) = xlayer(1:m)
      end if
    else
      do j=1,dofs
        if ( ParabolicModel ) then
          Var3D % Values(dofs*(Var3D % Perm(InvPerm)-1)+j) = sqrt( 2 * xlayer(j:m:dofs))
        else
          Var3D % Values(dofs*(Var3D % Perm(InvPerm)-1)+j) = xlayer(j:m:dofs)
        end if
      end do
    end if
  end subroutine SetLayerValues

  ! Reverse of the previous
  !-----------------------------------------------------------------------------------
  subroutine GetLayerValues( xlayer )
    real(kind=dp) :: xlayer(:)
    integer :: j
    if ( dofs == 1 ) then
      if ( ParabolicModel ) then
        xlayer(1:n) = 0.5_dp * Var3D % Values(Var3D % Perm(InvPerm))**2
      else
        xlayer(1:n) = Var3D % Values(Var3D % Perm(InvPerm))
      end if
    else
      do j=1,dofs
        if ( ParabolicModel ) then
          xlayer(j:m:dofs) = 0.5_dp * Var3D % Values(dofs*(Var3D % Perm(InvPerm)-1)+j)**2
        else
          xlayer(j:m:dofs) = Var3D % Values(dofs*(Var3D % Perm(InvPerm)-1)+j)
        end if
      end do
    end if
  end subroutine GetLayerValues


  ! Creates the local matrix equation for the ODY before time integration.
  ! We use namespace in order to allow the same solver to be used for
  ! several fields. The equation is of type:
  ! c*du/dt + r*u = f.
  !-----------------------------------------------------------------------
  subroutine GetCoefficients( q, Set0, xlayer )
    real(kind=DP), optional :: q
    logical, optional :: Set0
    real(kind=dp), optional :: xlayer(:)
    integer :: i,j,k,l
    type(element_t), pointer :: Element => null()
    real(kind=dp), pointer :: farray(:,:) => null()
    integer :: interp
    real(kind=dp) :: addv(2), qadd
    logical :: Found

    if ( UseInternalVals ) then
      if (.not. present(xlayer) ) then
        call Fatal(Caller,'"xlayer" must be provided as a parameter!')
      end if

      if (.not. Visited) then
        cdofs = dofs
        if (associated(AddVar)) then
          cdofs = cdofs + 1
        end if
        allocate(xloc(cdofs))
        call Info(Caller,'Get ODE coefficients using '//I2S(cdofs)//' internal variables!')
        if ( dofs == 1 .or. VectorSource ) then
          allocate(Source_h(1))
          call ListInitElementKeyword( Source_h(1),'Material',trim(VarName)//': Source', &
              DummyCount = cdofs )
        else
          allocate(Source_h(dofs))
          do j=1,dofs
            call ListInitElementKeyword( Source_h(j),'Material',trim(VarName)//' '//I2S(j)//': Source', &
                DummyCount = cdofs )
          end do
        end if
        Visited = .TRUE.
      end if

      do i=1,n
        j = InvPerm(i)
        if (j==0) call Fatal('GetCoefficients','We should have positive index!')
        Element => Mesh % Elements( ParentElem(j) )
        Model % CurrentElement => Element

        ! field values of the additional variable (most likely temperature)
        if (associated(AddVar)) then
          addv = 0.0_dp
          do interp=1,2
            if ( AddVar % TYPE == Variable_on_nodes_on_elements ) then
              k = 0
              do l=1,Element % TYPE % NumberOfNodes
                if (Element % NodeIndexes(l) == j) k = Element % DGIndexes(l)
              end do
              if (k==0) then
                call Fatal('GetCoefficients','Could not define DG index!')
              end if
            else
              k = j
            end if

            addv(interp) = AddVar % Values(AddVar % Perm(k))
            if (maxsubt == 1 .or. layer == 0) exit

            if (interp == 1) then
              j = PrevInvPerm(i)
            end if
          end do

          if (layer == 0 .or. maxsubt == 1 ) then
            xloc(1) = addv(1)
          else
            qadd = 1.0_dp * subt / maxsubt
            xloc(1) = qadd * addv(1) + (1-qadd) * addv(2)
          end if

          k = 1
        else
          k = 0
        end if

        ! the field values of this solver
        do j=1,dofs
          xloc(k+j) = xlayer(dofs*(i-1)+j)
        end do

        if ( VectorSource ) then
          j = 0
          qadd = ListGetElementReal( Source_h(1), Element = Element, Found = Found, DummyVals = xloc, Rdim=j, Rtensor=farray )
          do j=1,dofs
            l = dofs*(i-1)+j
            fvec(l:l) = farray(j,1)
          end do
        else
          do j=1,dofs
            l = dofs*(i-1)+j
            fvec(l:l) = ListGetElementReal( Source_h(j), Element = Element, Found = Found, DummyVals = xloc )
          end do
        end if
      end do

    else if (MaskExist .or. AnyDG ) then
      k = 1
      do i=1,n
        j = InvPerm(i)
        SingleIndex(1) = j

        if (j==0) call Fatal('GetCoefficients','We should have positive index!')
        Element => Mesh % Elements( ParentElem(j) )

        Model % CurrentElement => Element

        if (dofs == 1 ) then
          if ( HaveF ) then
            fvec(i:i) = ListGetReal( Material,&
                trim(VarName)//': Source',k,SingleIndex(1:1))
          end if
          if ( HaveR ) then
            rvec(i:i) = ListGetReal( Material,&
                trim(VarName)//': Reaction Coefficient',k,SingleIndex(1:1) )
          end if
          if ( HaveC ) then
            cvec(i:i) = ListGetReal( Material,&
                trim(VarName)//': Time Derivative Coefficient',k,SingleIndex(1:1) )
          end if
        else
          do j=1,dofs
            l = dofs*(i-1)+j
            if ( HaveF ) then
              fvec(l:l) = ListGetReal( Material,&
                  trim(VarName)//' '//I2S(j)//': Source',k,SingleIndex(1:1), Found )
            end if
            if ( HaveR ) then
              rvec(l:l) = ListGetReal( Material,&
                  trim(VarName)//' '//I2S(j)//': Reaction Coefficient',k,SingleIndex(1:1), Found )
            end if
            if ( HaveC ) then
              cvec(l:l) = ListGetReal( Material,&
                  trim(VarName)//' '//I2S(j)//': Time Derivative Coefficient',k,SingleIndex(1:1), Found)
              if (.not. Found) cvec(l:l) = 1.0_dp
            end if
          end do
        end if

        if (all(SingleIndex(1) /= Element % NodeIndexes ) ) stop
      end do
    else
      if ( dofs == 1 ) then
        if ( HaveF ) then
          fvec(1:n) = ListGetReal( Material,&
              trim(VarName)//': Source',n,InvPerm )
        end if
        if ( HaveR ) then
          rvec(1:n) = ListGetReal( Material,&
              trim(VarName)//': Reaction Coefficient',n,InvPerm )
        end if
        if ( HaveC ) then
          cvec(1:n) = ListGetReal( Material,&
              trim(VarName)//': Time Derivative Coefficient',n,invPerm )
        end if
      else
        do j=1,dofs
          if ( HaveF ) then
            fvec(j:m:dofs) = ListGetReal( Material,&
                trim(VarName)//' '//I2S(j)//': Source',n,InvPerm, Found )
          end if
          if ( HaveR ) then
            rvec(j:m:dofs) = ListGetReal( Material,&
                trim(VarName)//' '//I2S(j)//': Reaction Coefficient',n,InvPerm, Found )
          end if
          if ( HaveC ) then
            cvec(j:m:dofs) = ListGetReal( Material,&
                trim(VarName)//' '//I2S(j)//': Time Derivative Coefficient',n,InvPerm,Found)
            if (.not. Found) cvec(j:m:dofs) = 1.0_dp
          end if
        end do
      end if

    end if

    ! When using different integration we may need to access the
    ! value of coefficients at previous mesh layer.
    if ( present( Set0 ) ) then
      if ( Set0 ) then
        if ( HaveF ) f0vec = fvec
        if ( HaveR ) r0vec = rvec
        if ( HaveC ) c0vec = cvec
        return
      end if
    end if

    if (.not. present( q ) ) return
    if ( abs(q-1.0_dp) < epsilon(q) ) return

    if ( HaveF ) fvec = q * fvec + (1-q) * f0vec
    if ( HaveR ) rvec = q * rvec + (1-q) * r0vec
    if ( HaveC ) cvec = q * cvec + (1-q) * c0vec

  end subroutine GetCoefficients

end module marchingodesolver_acc
