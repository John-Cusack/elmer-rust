! T021-B port of fem/src/modules/PhaseChangeSolve.F90: statement-level
! BIND(C) wrappers compiled against the built .mod files of the target
! configure (ADR-057 §3, T020 Poisson pattern).
!
! Every wrapper takes only C-interoperable arguments. All Fortran-object
! traffic stays inside Fortran. SAVE state mirrors the original's SAVE
! statement with identical lifetimes: first-mesh heap and latched scalars
! live as module variables (implicitly SAVE); per-call associations are
! refreshed by pca_begin / the bind wrappers on every entry, exactly
! where the original re-associates its locals. The three CONTAINS
! procedures (VelocityLocalMatrix, SurfaceLocalMatrix, FindPullBoundary)
! are kept as Fortran module procedures transcribed statement for
! statement, so all order-sensitive reductions (SUM/MAXVAL), formatted
! I/O (WRITE to Message, PRINT *) and whole-array updates hold
! bit-identity by construction under the same compiler and flags.
!
! Dead state omitted with zero observable effect (verified against the
! source: assigned but never read, or declared but never referenced):
! PNodes (host copy, never touched), Weights, ForceVector (host SAVE,
! shadowed per call), FluxCorrect, AllocationsDone (set, never read),
! PrevUpull, NodalTemp host values (allocated for Fatal parity but never
! filled, as in the original), and the unused host scalars/logicals.
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module pca_acc
  use iso_c_binding
  use DefUtils
  implicit none

  ! -- entry context (per-call: refreshed by pca_begin each entry) --------
  type(Model_t), pointer :: Model => null()
  type(Solver_t), pointer :: Solver => null()
  type(Solver_t), pointer :: PSolver => null()
  type(Mesh_t), pointer :: Mesh => null()
  type(ValueList_t), pointer :: Params => null()
  type(ValueList_t), pointer :: Material => null()
  type(Variable_t), pointer :: SurfSol => null()
  type(Variable_t), pointer :: TempSol => null()
  type(Variable_t), pointer :: HelpSol => null()
  type(Variable_t), pointer :: HelpSol2 => null()
  type(Matrix_t), pointer :: StiffMatrix => null()
  real(kind=dp), pointer :: Surface(:), Temperature(:)
  real(kind=dp), pointer :: ForceVector(:)
  integer, pointer :: SurfPerm(:), TempPerm(:), NodeIndexes(:)
  type(Element_t), pointer :: CurrentElement => null()

  ! -- SAVE-latched state (mirrors the SAVE statement lifetimes) ---------
  type(Nodes_t) :: Nodes
  real(kind=dp), pointer :: x(:), y(:), z(:), Basis(:), dBasisdx(:,:)
  real(kind=dp), pointer :: NodalTemp(:), Conductivity(:), LatentHeat(:), TempDiff(:)
  real(kind=dp), pointer :: LocalStiffMatrix(:,:), LocalForceVector(:), LocalMassMatrix(:,:)
  real(kind=dp), pointer :: PrevTemp(:), SurfaceMove(:), SurfaceMoveAve(:)
  real(kind=dp), pointer :: Normals(:)
  integer, pointer :: NormalsPerm(:)
  logical, pointer :: NodeDone(:)
  logical :: AverageNormal = .false.
  integer :: NoBNodes = 0, Trip_node = 0, Axis_node = 0
  integer :: NormalDir = 0, TangentDirection = 0
  real(kind=dp) :: CoordMax(3) = 0.0_dp, CoordMin(3) = 0.0_dp
  integer :: CoordMaxi(3) = 0, CoordMini(3) = 0
  real(kind=dp) :: UPull(3) = 0.0_dp
  logical :: PullVelocitySet = .false.
  logical :: SurfaceVelocitySet = .false.
  ! Per-call config mirrors (re-read every entry, as in the original):
  integer :: m_dim = 0
  real(kind=dp) :: m_density = 0.0_dp
  real(kind=dp) :: m_meltpoint = 0.0_dp
  logical :: m_pullcontrol = .false., m_triplefixed = .false., m_useloads = .false.
  ! Steady per-entry accumulators + cross-call latches:
  real(kind=dp) :: m_area = 0.0_dp, m_volume = 0.0_dp, m_tave = 0.0_dp
  real(kind=dp) :: m_tabs = 0.0_dp, m_volabs = 0.0_dp, m_maxtempdiff = 0.0_dp
  real(kind=dp) :: m_dpos = 0.0_dp, m_pos0 = 0.0_dp, m_prevpos0 = 0.0_dp
  integer :: m_visited = -1
  real(kind=dp) :: m_ccum = 1.0_dp
  real(kind=dp) :: m_prevvolume = 0.0_dp, m_prevtave = 0.0_dp
  real(kind=dp) :: m_prevtabs = 0.0_dp, m_prevvolabs = 0.0_dp
  real(kind=dp) :: m_newton_after_tol = 0.0_dp
  integer :: m_nelems = 0
  real(kind=dp) :: m_xmin = 0.0_dp, m_xmax = 0.0_dp
  real(kind=dp), pointer :: m_taverage(:), m_taverage2(:)
  real(kind=dp), allocatable :: m_isosurf(:,:)
  logical :: m_isoalloc = .false.
  ! Host-associated mirrors for the CONTAINS procedures:
  real(kind=dp) :: Normal(3) = 0.0_dp
  real(kind=dp) :: SpeedUp = 1.0_dp, m_dt = 0.0_dp
  logical :: m_transient = .false.
  ! MeltingHeat SAVE mirror (Nodes buffer allocated once; the lookup
  ! itself re-runs every call, exactly as in the original).
  type(Nodes_t), save :: mh_nodes
  logical, save :: mh_nodes_alloc = .false.

contains

  ! -- entry -----------------------------------------------------------
  subroutine pca_begin(h_model, h_solver) bind(c, name="pca_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, Model)
    call c_f_pointer(h_solver, Solver)
    PSolver => Solver
    Mesh => Solver % Mesh
    Params => GetSolverParams()
  end subroutine pca_begin

  subroutine pca_info_start() bind(c, name="pca_info_start")
    call Info('PhaseChangeSolve','--------------------------------------------')
  end subroutine pca_info_start

  ! Lines 141-151: algorithm decision + the two branch Info lines.
  function pca_algo(dtv, transient) result(steady) bind(c, name="pca_algo")
    real(c_double), value :: dtv
    integer(c_int), value :: transient
    integer(c_int) :: steady
    logical :: stat, transalgo
    real(kind=dp) :: steadydt
    transalgo = (transient /= 0)
    m_transient = transalgo
    m_dt = dtv
    if (transalgo) then
      steadydt = ListGetConstReal(Solver % Values,'Steady Transition Timestep',stat)
      if (stat) transalgo = (dtv < steadydt)
    end if
    if (transalgo) then
      call Info('PhaseChangeSolve','Using transient algorithm for surface update')
    else
      call Info('PhaseChangeSolve','Using steady algorithm to find the isotherm')
    end if
    call Info('PhaseChangeSolve','--------------------------------------------')
    steady = merge(1, 0, .not. transalgo)
  end function pca_algo

  function pca_coord_dim() result(d) bind(c, name="pca_coord_dim")
    integer(c_int) :: d
    m_dim = CoordinateSystemDimension()
    d = m_dim
  end function pca_coord_dim

  subroutine pca_fatal_2d() bind(c, name="pca_fatal_2d")
    call Fatal('PhaseChangeSolve','Implemented only in 2D')
  end subroutine pca_fatal_2d

  ! -- primary / temperature bindings -----------------------------------
  subroutine pca_bind_primary() bind(c, name="pca_bind_primary")
    SurfSol => Solver % Variable
    Surface => SurfSol % Values
    SurfPerm => SurfSol % Perm
    if (.not. associated(Surface) .or. all(SurfPerm <= 0)) then
      call Fatal('PhaseChangeSolve','Surface field needed for Phase Change')
    end if
  end subroutine pca_bind_primary

  subroutine pca_bind_temp() bind(c, name="pca_bind_temp")
    logical :: stat
    character(len=MAX_NAME_LEN) :: vname
    vname = ListGetString(Solver % Values, 'Phase Change Variable', stat)
    if (stat) then
      TempSol => VariableGet(Solver % Mesh % Variables, trim(vname))
    else
      TempSol => VariableGet(Solver % Mesh % Variables, 'Temperature')
    end if
    TempPerm => TempSol % Perm
    Temperature => TempSol % Values
    if (.not. associated(Temperature) .or. all(TempPerm <= 0)) then
      call Fatal('PhaseChangeSolve','Temperature field needed for Phase Change')
    end if
  end subroutine pca_bind_temp

  ! -- per-call solver parameters (lines 184-194) -------------------------
  subroutine pca_cfg(relax, niter, ntol) bind(c, name="pca_cfg")
    real(c_double) :: relax, ntol
    integer(c_int) :: niter
    logical :: stat
    relax = GetCReal(Solver % Values, 'Nonlinear System Relaxation Factor', stat)
    if (.not. stat) relax = 1.0d0
    niter = ListGetInteger(Solver % Values, 'Nonlinear System Max Iterations', stat)
    if (.not. stat) niter = 1
    ntol = ListGetConstReal(Solver % Values, 'Nonlinear System Convergence Tolerance', stat)
    m_pullcontrol = ListGetLogical(Solver % Values, 'Pull Rate Control', stat)
    m_triplefixed = ListGetLogical(Solver % Values, 'Triple Point Fixed', stat)
  end subroutine pca_cfg

  ! -- first-time setup (lines 201-313) -----------------------------------
  subroutine pca_firsttime(steady) bind(c, name="pca_firsttime")
    integer(c_int), value :: steady
    logical :: stat
    character(len=MAX_NAME_LEN) :: vname
    integer :: i, j, k, n, istat
    real(kind=dp) :: xx
    UPull = 0.0d0
    NoBNodes = 0
    CoordMax = -huge(CoordMax)
    CoordMin = huge(CoordMin)
    do k = 1, Model % Mesh % NumberOfNodes
      if (SurfPerm(k) <= 0) cycle
      NoBNodes = NoBNodes + 1
      do j = 1, m_dim
        if (j == 1) xx = Model % Mesh % Nodes % x(k)
        if (j == 2) xx = Model % Mesh % Nodes % y(k)
        if (j == 3) xx = Model % Mesh % Nodes % z(k)
        if (xx > CoordMax(j)) then
          CoordMax(j) = xx
          CoordMaxi(j) = k
        end if
        if (xx < CoordMin(j)) then
          CoordMin(j) = xx
          CoordMini(j) = k
        end if
      end do
    end do
    j = 1
    do i = 1, m_dim
      if (CoordMax(i) - CoordMin(i) < CoordMax(j) - CoordMin(j)) j = i
    end do
    NormalDir = j
    j = 1
    do i = 1, m_dim
      if (CoordMax(i) - CoordMin(i) > CoordMax(j) - CoordMin(j)) j = i
    end do
    TangentDirection = j
    Trip_node = CoordMaxi(TangentDirection)
    Axis_node = CoordMini(TangentDirection)
    n = Solver % Mesh % MaxElementNodes
    allocate(Nodes % x(n), Nodes % y(n), Nodes % z(n), &
        x(n), y(n), z(n), Basis(n), dBasisdx(n,3), NodalTemp(n), &
        Conductivity(n), LatentHeat(n), TempDiff(n), &
        LocalStiffMatrix(n,n), LocalForceVector(n), LocalMassMatrix(n,n), &
        PrevTemp(NoBNodes), stat=istat)
    if (istat /= 0) call Fatal('PhaseChangeSolver', 'Memory allocation error 1.')
    Nodes % x = 0.0d0
    Nodes % y = 0.0d0
    Nodes % z = 0.0d0
    PrevTemp = 0.0d0
    if (steady /= 0) then
      allocate(NodeDone(NoBNodes), stat=istat)
      if (istat /= 0) call Fatal('PhaseChangeSolver', 'Memory allocation error 3.')
    end if
    vname = ListGetString(Solver % Values, 'Normal Variable', stat)
    if (stat) then
      HelpSol => VariableGet(Solver % Mesh % Variables, trim(vname), ThisOnly=.true.)
    else
      HelpSol => VariableGet(Solver % Mesh % Variables, 'Normals', ThisOnly=.true.)
    end if
    if (associated(HelpSol)) then
      Normals => HelpSol % Values
      NormalsPerm => HelpSol % Perm
      AverageNormal = .true.
    else
      AverageNormal = .false.
    end if
    HelpSol => VariableGet(Solver % Mesh % Variables, &
        trim(ComponentName(Solver % Variable))//'Move', stat)
    if (associated(HelpSol)) then
      SurfaceMove => HelpSol % Values
    else
      allocate(SurfaceMove(size(Surface)), stat=istat)
      if (istat /= 0) call Fatal('PhaseChangeSolver', 'Memory allocation error.')
      SurfaceMove = 0.0d0
      call VariableAdd(Solver % Mesh % Variables, Solver % Mesh, &
          PSolver, trim(ComponentName(Solver % Variable))//'Move', 1, &
          SurfaceMove, SurfPerm, Output=.false.)
    end if
    if (ListGetLogical(Solver % Values, 'Use Average Velocity', stat)) then
      HelpSol => VariableGet(Solver % Mesh % Variables, &
          trim(ComponentName(Solver % Variable))//'MoveAve', stat)
      if (associated(HelpSol)) then
        SurfaceMoveAve => HelpSol % Values
      else
        allocate(SurfaceMoveAve(size(Surface)), stat=istat)
        if (istat /= 0) call Fatal('PhaseChangeSolver', 'Memory allocation error.')
        SurfaceMoveAve = 0.0d0
        call VariableAdd(Solver % Mesh % Variables, Solver % Mesh, &
            PSolver, trim(ComponentName(Solver % Variable))//'MoveAve', &
            1, SurfaceMoveAve, SurfPerm, Output=.false.)
      end if
    end if
  end subroutine pca_firsttime

  function pca_trip_temp() result(t) bind(c, name="pca_trip_temp")
    real(c_double) :: t
    t = Temperature(TempPerm(Trip_node))
  end function pca_trip_temp

  ! -- Passive Steps short-circuit (lines 321-322) -------------------------
  function pca_passive(visit) result(skip) bind(c, name="pca_passive")
    integer(c_int), value :: visit
    integer(c_int) :: skip
    logical :: stat
    integer :: i
    i = ListGetInteger(Solver % Values, 'Passive Steps', stat)
    skip = merge(1, 0, i >= visit)
  end function pca_passive

  ! -- material parameters + pull velocity latch (lines 327-341) -----------
  subroutine pca_material() bind(c, name="pca_material")
    logical :: stat
    integer :: k, n
    CurrentElement => Solver % Mesh % Elements(Solver % ActiveElements(1))
    n = CurrentElement % TYPE % NumberOfNodes
    NodeIndexes => CurrentElement % NodeIndexes
    k = ListGetInteger(Model % Bodies(CurrentElement % BodyId) % Values, 'Material')
    Material => Model % Materials(k) % Values
    m_meltpoint = ListGetConstReal(Material, 'Melting Point')
    m_density = ListGetConstReal(Material, 'Density')
    if (.not. PullVelocitySet) then
      UPull(1) = ListGetConstReal(Material, 'Convection Velocity 1', stat)
      UPull(2) = ListGetConstReal(Material, 'Convection Velocity 2', stat)
      UPull(3) = 0.0d0
      PullVelocitySet = .true.
    end if
  end subroutine pca_material


  ! -- transient branch --------------------------------------------------
  ! Lines 351-361: matrix association, heat-load binding, PrevTemp snap.
  function pca_trans_begin() result(ul) bind(c, name="pca_trans_begin")
    integer(c_int) :: ul
    logical :: stat
    StiffMatrix => Solver % Matrix
    ForceVector => Solver % Matrix % RHS
    m_useloads = ListGetLogical(Solver % Values, 'Use Heat Load', stat)
    if (m_useloads) then
      HelpSol => VariableGet(Solver % Mesh % Variables, 'Nodal Heat Load')
    end if
    PrevTemp = SurfaceMove
    ul = merge(1, 0, m_useloads)
  end function pca_trans_begin

  subroutine pca_vel_zero() bind(c, name="pca_vel_zero")
    call InitializeToZero(StiffMatrix, ForceVector)
  end subroutine pca_vel_zero

  ! Lines 369-375: direct heat-load fill (no nodal-loads machinery here).
  subroutine pca_vel_loads() bind(c, name="pca_vel_loads")
    integer :: t, i
    do t = 1, Solver % Mesh % NumberOfNodes
      i = SurfPerm(t)
      if (i <= 0) cycle
      ForceVector(i) = HelpSol % Values(HelpSol % Perm(t))
    end do
  end subroutine pca_vel_loads

  function pca_nactive() result(n) bind(c, name="pca_nactive")
    integer(c_int) :: n
    n = Solver % NumberOfActiveElements
  end function pca_nactive

  ! Lines 377-394: per-element gather + VelocityLocalMatrix + update.
  subroutine pca_vel_elem(t) bind(c, name="pca_vel_elem")
    integer(c_int), value :: t
    logical :: stat
    integer :: k, n
    CurrentElement => Solver % Mesh % Elements(Solver % ActiveElements(t))
    n = CurrentElement % TYPE % NumberOfNodes
    NodeIndexes => CurrentElement % NodeIndexes
    k = ListGetInteger(Model % Bodies(CurrentElement % BodyId) % Values, 'Material')
    Material => Model % Materials(k) % Values
    LatentHeat(1:n) = ListGetReal(Material, 'Latent Heat', n, NodeIndexes)
    Conductivity(1:n) = ListGetReal(Material, 'Heat Conductivity', n, NodeIndexes)
    TempDiff(1:n) = Temperature(TempPerm(NodeIndexes)) - m_meltpoint
    call VelocityLocalMatrix(LocalStiffMatrix, LocalMassMatrix, LocalForceVector, &
        CurrentElement, n, Nodes)
    call UpdateGlobalEquations(StiffMatrix, LocalStiffMatrix, &
        ForceVector, LocalForceVector, n, 1, SurfPerm(NodeIndexes))
  end subroutine pca_vel_elem

  subroutine pca_vel_asm() bind(c, name="pca_vel_asm")
    call FinishAssembly(Solver, ForceVector)
  end subroutine pca_vel_asm

  ! Lines 402-409: solve + the two Level-4 Info reports.
  subroutine pca_vel_solve(norm, relchange) bind(c, name="pca_vel_solve")
    real(c_double) :: norm, relchange
    real(kind=dp) :: lnorm
    call SolveSystem(StiffMatrix, ParMatrix, ForceVector, SurfaceMove, lnorm, 1, Solver)
    relchange = Solver % Variable % NonlinChange
    write(Message, *) 'Result Norm     : ', lnorm
    call Info('PhaseChangeSolve', Message, level=4)
    write(Message, *) 'Relative Change : ', relchange
    call Info('PhaseChangeSolve', Message, level=4)
    norm = lnorm
  end subroutine pca_vel_solve

  ! Lines 414-429: velocity averaging / relaxation selection.
  subroutine pca_vel_average() bind(c, name="pca_vel_average")
    logical :: stat
    real(kind=dp) :: localrelax
    if (ListGetLogical(Solver % Values, 'Use Average Velocity', stat)) then
      if (SurfaceVelocitySet) then
        localrelax = GetCReal(Solver % Values, 'Velocity Averaging Factor')
        SurfaceMoveAve = localrelax * SurfaceMove + (1.0d0 - localrelax) * SurfaceMoveAve
      else
        SurfaceMoveAve = SurfaceMove
        SurfaceVelocitySet = .true.
      end if
      localrelax = GetCReal(Solver % Values, 'Velocity Relaxation Factor', stat)
      if (.not. stat) localrelax = 1.0d0
      SurfaceMove = localrelax * SurfaceMove + (1.0d0 - localrelax) * SurfaceMoveAve
    else
      localrelax = GetCReal(Solver % Values, 'Velocity Relaxation Factor', stat)
      if (.not. stat) localrelax = 1.0d0
      SurfaceMove = localrelax * SurfaceMove + (1.0d0 - localrelax) * PrevTemp
    end if
  end subroutine pca_vel_average

  ! Lines 431-437: pull-velocity feedback from the triple point.
  subroutine pca_pull_update() bind(c, name="pca_pull_update")
    if (m_pullcontrol .or. m_triplefixed) then
      UPull(NormalDir) = -SurfaceMove(SurfPerm(Trip_node))
      if (m_pullcontrol) then
        write(Message, *) 'Pull velocity: ', UPull(NormalDir)
        call Info('PhaseChangeSolve', Message)
      end if
    end if
  end subroutine pca_pull_update

  function pca_speedup() result(s) bind(c, name="pca_speedup")
    real(c_double) :: s
    logical :: stat
    SpeedUp = ListGetConstReal(Solver % Values, 'Transient SpeedUp', stat)
    if (.not. stat) SpeedUp = 1.0d0
    s = SpeedUp
  end function pca_speedup

  ! Line 447: the smoothing-factor read is fatal-if-missing, verbatim.
  function pca_smooth() result(s) bind(c, name="pca_smooth")
    real(c_double) :: s
    s = ListGetConstReal(Solver % Values, 'Surface Smoothing Factor')
  end function pca_smooth

  ! Line 448: direct surface update (no smoothing).
  subroutine pca_simple_update(dtv) bind(c, name="pca_simple_update")
    real(c_double), value :: dtv
    Surface = Surface + SpeedUp * dtv * (SurfaceMove + UPull(NormalDir))
  end subroutine pca_simple_update

  subroutine pca_sur_begin() bind(c, name="pca_sur_begin")
    call InitializeToZero(StiffMatrix, ForceVector)
  end subroutine pca_sur_begin

  ! Lines 452-464: per-element gather + SurfaceLocalMatrix + 1st-order time.
  subroutine pca_sur_elem(t, dtv) bind(c, name="pca_sur_elem")
    integer(c_int), value :: t
    real(c_double), value :: dtv
    integer :: n
    CurrentElement => Solver % Mesh % Elements(Solver % ActiveElements(t))
    n = CurrentElement % TYPE % NumberOfNodes
    NodeIndexes => CurrentElement % NodeIndexes
    call SurfaceLocalMatrix(LocalStiffMatrix, LocalMassMatrix, LocalForceVector, &
        CurrentElement, n, Nodes)
    call Add1stOrderTime(LocalMassMatrix, LocalStiffMatrix, &
        LocalForceVector, dtv, n, 1, SurfPerm(NodeIndexes), Solver)
    call UpdateGlobalEquations(StiffMatrix, LocalStiffMatrix, &
        ForceVector, LocalForceVector, n, 1, SurfPerm(NodeIndexes))
  end subroutine pca_sur_elem

  function pca_sur_solve() result(norm) bind(c, name="pca_sur_solve")
    real(c_double) :: norm
    real(kind=dp) :: lnorm
    call FinishAssembly(Solver, ForceVector)
    call SolveSystem(StiffMatrix, ParMatrix, ForceVector, Surface, lnorm, 1, Solver)
    norm = lnorm
  end function pca_sur_solve

  ! Lines 471-489: pulled-distance integration + FindPullBoundary.
  subroutine pca_pullpos(dtv) bind(c, name="pca_pullpos")
    real(c_double), value :: dtv
    if (m_pullcontrol .or. m_triplefixed) then
      if (m_visited /= Solver % DoneTime) then
        if (Solver % DoneTime == 1) then
          m_prevpos0 = 0.0d0
        else
          m_prevpos0 = m_pos0
        end if
        m_visited = Solver % DoneTime
      end if
      m_pos0 = m_prevpos0 + UPull(NormalDir) * dtv
      if (m_pullcontrol) call FindPullBoundary()
    end if
  end subroutine pca_pullpos

  ! -- steady branch -----------------------------------------------------
  ! Lines 505-545: optional temperature averaging (Taverage/Taverage Slow).
  subroutine pca_taverage(dtv) bind(c, name="pca_taverage")
    real(c_double), value :: dtv
    logical :: stat
    real(kind=dp) :: averagingdt, averagerelax, averagerelax2
    integer :: istat
    if (.not. ListGetLogical(Solver % Values, 'Average Temperature', stat)) return
    averagingdt = ListGetConstReal(Solver % Values, 'Steady Averaging Timestep', stat)
    HelpSol => VariableGet(Solver % Mesh % Variables, 'Taverage')
    if (.not. associated(HelpSol)) then
      allocate(m_taverage(Model % NumberOfNodes), stat=istat)
      if (istat /= 0) call Fatal('PhaseChangeSolver', 'Memory allocation error.')
      call VariableAdd(Solver % Mesh % Variables, Solver % Mesh, &
          PSolver, 'Taverage', 1, m_taverage, TempPerm)
      HelpSol => VariableGet(Solver % Mesh % Variables, 'Taverage')
    end if
    nullify(HelpSol2)
    averagerelax2 = GetCReal(Solver % Values, 'Temperature Averaging Factor', stat)
    if (stat) then
      HelpSol2 => VariableGet(Solver % Mesh % Variables, 'Taverage Slow')
      if (.not. associated(HelpSol2)) then
        allocate(m_taverage2(Model % NumberOfNodes), stat=istat)
        if (istat /= 0) call Fatal('PhaseChangeSolver', 'Memory allocation error.')
        call VariableAdd(Solver % Mesh % Variables, Solver % Mesh, &
            PSolver, 'Taverage Slow', 1, m_taverage2, TempPerm)
        HelpSol2 => VariableGet(Solver % Mesh % Variables, 'Taverage Slow')
      end if
    end if
    if (dtv >= averagingdt) then
      HelpSol % Values = TempSol % Values
      if (associated(HelpSol2)) HelpSol2 % Values = TempSol % Values
    else
      averagerelax = GetCReal(Solver % Values, 'Temperature Relaxation Factor')
      if (.not. associated(HelpSol2)) then
        HelpSol % Values = (1.0d0 - averagerelax) * HelpSol % Values + averagerelax * TempSol % Values
      else
        HelpSol2 % Values = (1.0d0 - averagerelax2) * HelpSol2 % Values + averagerelax2 * TempSol % Values
        HelpSol % Values = (1.0d0 - averagerelax) * HelpSol2 % Values + averagerelax * TempSol % Values
      end if
      print *, 'temperature set to average temperature'
      Temperature => HelpSol % Values
    end if
  end subroutine pca_taverage

  ! Lines 547-558: Newton formulation latch + branch Info lines.
  function pca_newton_begin(visit, newton_in) result(nw) bind(c, name="pca_newton_begin")
    integer(c_int), value :: visit, newton_in
    integer(c_int) :: nw
    logical :: stat
    integer :: newtonafteriter
    logical :: newton
    newton = (newton_in /= 0)
    m_newton_after_tol = 0.0d0
    newtonafteriter = ListGetInteger(Solver % Values, &
        'Nonlinear System Newton After Iterations', stat)
    if (stat .and. visit > newtonafteriter) newton = .true.
    m_newton_after_tol = ListGetConstReal(Solver % Values, &
        'Nonlinear System Newton After Tolerance', stat)
    if (newton) then
      call Info('PhaseChangeSolve', 'Steady state newton formulation', level=4)
    else
      call Info('PhaseChangeSolve', 'Steady state isoterm formulation', level=4)
    end if
    nw = merge(1, 0, newton)
  end function pca_newton_begin

  ! Lines 563-567: triple-point melting point.
  subroutine pca_triple_melt(trip_temp) bind(c, name="pca_triple_melt")
    real(c_double), value :: trip_temp
    if (m_triplefixed) then
      m_meltpoint = Temperature(TempPerm(Trip_node))
      write(Message, *) 'Melting point set to triple point temperature', m_meltpoint, Trip_node
      call Info('PhaseChangeSolve', Message)
    end if
  end subroutine pca_triple_melt

  ! Lines 573-667: isotherm construction (two-pass with GOTO 100).
  subroutine pca_createisotherm() bind(c, name="pca_createisotherm")
    logical :: stat
    integer :: t, n, nn, next, vertex, elementcode, imin, k
    integer :: nelems
    real(kind=dp) :: temp1, temp2, temppi, xmin, xmax
    m_isoalloc = .false.
    xmin = huge(xmin)
    xmax = -huge(xmax)
100 nelems = 0
    do t = 1, Solver % Mesh % NumberOfBulkElements
      CurrentElement => Solver % Mesh % Elements(t)
      n = CurrentElement % TYPE % NumberOfNodes
      NodeIndexes => CurrentElement % NodeIndexes
      elementcode = CurrentElement % TYPE % ElementCode
      k = ListGetInteger(Model % Bodies(CurrentElement % BodyId) % Values, 'Material')
      if (.not. (ListGetLogical(Model % Materials(k) % Values, 'Solid', stat) .or. &
          ListGetLogical(Model % Materials(k) % Values, 'Liquid', stat))) cycle
      if (elementcode < 300) cycle
      if (any(TempPerm(NodeIndexes(1:n)) <= 0)) cycle
      TempDiff(1:n) = Temperature(TempPerm(NodeIndexes(1:n))) - m_meltpoint
      if (all(TempDiff(1:n) < 0.0d0)) cycle
      if (all(TempDiff(1:n) > 0.0d0)) cycle
      Nodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
      Nodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
      vertex = elementcode / 100
      n = 0
      do nn = 1, vertex
        next = modulo(nn, vertex) + 1
        temp1 = TempDiff(nn)
        temp2 = TempDiff(next)
        if (((temp1 < 0.0d0) .and. (0.0d0 <= temp2)) .or. &
            ((temp2 <= 0.0d0) .and. (0.0d0 < temp1))) then
          n = n + 1
          if (n <= 2) then
            nelems = nelems + 1
            if (m_isoalloc) then
              m_isosurf(nelems, 1) = Nodes % x(nn) + &
                  temp1 * ((Nodes % x(next) - Nodes % x(nn)) / (temp1 - temp2))
              m_isosurf(nelems, 2) = Nodes % y(nn) + &
                  temp1 * ((Nodes % y(next) - Nodes % y(nn)) / (temp1 - temp2))
              xmin = min(m_isosurf(nelems, 1), xmin)
              xmax = max(m_isosurf(nelems, 1), xmax)
            end if
          else
            call Warn('PhaseChangeSolve', 'Wiggly Isotherm')
            write(Message, *) 'nodeindexes', NodeIndexes(1:vertex)
            call Warn('PhaseChangeSolve', Message)
            write(Message, *) 'temperature', Temperature(TempPerm(NodeIndexes))
            call Warn('PhaseChangeSolve', Message)
            print *, 'Nodes % x', Nodes % x(n)
            print *, 'Nodes % y', Nodes % y(n)
          end if
        end if
      end do
      if (n == 1) then
        nelems = nelems - 1
        cycle
      end if
      if (m_isoalloc .and. n == 2) then
        if (m_isosurf(nelems - 1, TangentDirection) > m_isosurf(nelems, TangentDirection)) then
          temppi = m_isosurf(nelems - 1, 1)
          m_isosurf(nelems - 1, 1) = m_isosurf(nelems, 1)
          m_isosurf(nelems, 1) = temppi
          temppi = m_isosurf(nelems - 1, 2)
          m_isosurf(nelems - 1, 2) = m_isosurf(nelems, 2)
          m_isosurf(nelems, 2) = temppi
        end if
      end if
    end do
    if (nelems == 0) &
        call Fatal('PhaseChangeSolve', 'Isotherm is empty thus cannot map phase change surface')
    if (.not. m_isoalloc) then
      allocate(m_isosurf(nelems + 1, 2))
      m_isoalloc = .true.
      write(Message, *) 'Isotherm created with number of segments', nelems
      call Info('PhaseChangeSolve', Message)
      goto 100
    end if
    m_nelems = nelems
    m_xmin = xmin
    m_xmax = xmax
  end subroutine pca_createisotherm

  ! Lines 669-826: accumulator init + free-surface mapping + integration.
  subroutine pca_steady_map(newton) bind(c, name="pca_steady_map")
    integer(c_int), value :: newton
    logical :: stat
    integer :: t, n, nn, i, imin, elementcode, k
    real(kind=dp) :: xx, yy, eps, dxmin, dymin, d, update, ttemp, dtdz
    real(kind=dp) :: u, v, w, s, detj
    type(GaussIntegrationPoints_t) :: integstuff
    m_area = 0.0d0
    m_volume = 0.0d0
    m_tave = 0.0d0
    m_tabs = 0.0d0
    m_volabs = 0.0d0
    NodeDone = .false.
    m_maxtempdiff = 0.0d0
    do t = 1, Solver % NumberOfActiveElements
      CurrentElement => Solver % Mesh % Elements(Solver % ActiveElements(t))
      n = CurrentElement % TYPE % NumberOfNodes
      NodeIndexes => CurrentElement % NodeIndexes
      elementcode = CurrentElement % TYPE % ElementCode
      if (elementcode < 200) cycle
      if (elementcode > 203) then
        call Fatal('PhaseChangeSolve', 'Implemented only for elements 202 and 203!')
        cycle
      end if
      Nodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
      Nodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
      Nodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)
      TempDiff(1:n) = Temperature(TempPerm(NodeIndexes(1:n))) - m_meltpoint
      m_maxtempdiff = max(m_maxtempdiff, maxval(abs(TempDiff(1:n))))
      do nn = 1, n
        k = SurfPerm(NodeIndexes(nn))
        if (NodeDone(k)) cycle
        NodeDone(k) = .true.
        if (m_triplefixed .and. NodeIndexes(nn) == Trip_node) then
          SurfaceMove(k) = 0.0d0
          print *, 'triple point fixed by construction'
          cycle
        end if
        if (nn == 3) then
          Surface(k) = 0.5d0 * sum(Surface(SurfPerm(NodeIndexes(1:2))))
          update = 0.0d0
        end if
        if (TangentDirection == 1) then
          xx = Nodes % x(nn)
          yy = Nodes % y(nn)
        else
          xx = Nodes % y(nn)
          yy = Nodes % x(nn)
        end if
        if (newton == 0) then
          eps = 1.0d-6 * (m_xmax - m_xmin)
          dxmin = huge(dxmin)
          dymin = huge(dymin)
          stat = .false.
          do i = 1, m_nelems - 1, 2
            if ((xx > m_isosurf(i, TangentDirection) - eps) .and. &
                (xx < m_isosurf(i + 1, TangentDirection) + eps)) then
              dxmin = 0.0d0
              d = min(abs(yy - m_isosurf(i, NormalDir)), &
                  abs(yy - m_isosurf(i + 1, NormalDir)))
              d = d + max(0.0d0, m_isosurf(i, TangentDirection) - xx)
              d = d + max(0.0d0, xx - m_isosurf(i + 1, TangentDirection))
              if (d <= dymin) then
                stat = .true.
                dymin = d
                imin = i
              end if
            end if
            if (.not. stat) then
              d = min(abs(xx - m_isosurf(i, TangentDirection)), &
                  abs(xx - m_isosurf(i + 1, TangentDirection)))
              if (d <= dxmin) then
                dxmin = d
                imin = i
              end if
            end if
          end do
          i = imin
          if (.not. stat) then
            if (dxmin > 1.0d-2 * abs(m_isosurf(i, TangentDirection) - &
                m_isosurf(i + 1, TangentDirection))) then
              call Warn('PhaseChangeSolve', 'Isotherm error?')
              write(Message, *) 'Nodeindexes', NodeIndexes(nn)
              call Warn('PhaseChangeSolve', Message)
              write(Message, *) 'x:', xx, ' y:', yy
              call Warn('PhaseChangeSolve', Message)
              write(Message, *) 'dxmin:', dxmin, ' dymin:', dymin
              call Warn('PhaseChangeSolve', Message)
            end if
          end if
          if (abs(m_isosurf(i + 1, TangentDirection) - &
              m_isosurf(i, TangentDirection)) > AEPS) then
            update = m_isosurf(i, NormalDir) + (xx - m_isosurf(i, TangentDirection)) * &
                (m_isosurf(i + 1, NormalDir) - m_isosurf(i, NormalDir)) / &
                (m_isosurf(i + 1, TangentDirection) - m_isosurf(i, TangentDirection)) - yy
          else
            update = 0.5d0 * (m_isosurf(i, NormalDir) + &
                m_isosurf(i + 1, NormalDir)) - yy
          end if
        end if
        ttemp = Temperature(TempPerm(NodeIndexes(nn)))
        if (newton /= 0) then
          dtdz = ttemp - PrevTemp(k)
          if (abs(dtdz) < AEPS) then
            call Warn('PhaseChangeSolve', 'Very small temperature update.')
            dtdz = 1.0d0
          end if
          update = SurfaceMove(k) * (m_meltpoint - ttemp) / dtdz
        end if
        PrevTemp(k) = ttemp
        SurfaceMove(k) = update
      end do
      integstuff = GaussPoints(CurrentElement)
      do i = 1, integstuff % n
        u = integstuff % u(i)
        v = integstuff % v(i)
        w = integstuff % w(i)
        stat = ElementInfo(CurrentElement, Nodes, u, v, w, detJ, Basis, dBasisdx)
        s = integstuff % s(i) * detJ
        if (CurrentCoordinateSystem() /= Cartesian) then
          s = s * sum(Basis(1:n) * Nodes % x(1:n)) * 2.0d0 * PI
        end if
        m_area = m_area + s
        m_volume = m_volume + s * sum(Basis(1:n) * SurfaceMove(SurfPerm(NodeIndexes(1:n))))
        m_tave = m_tave + s * sum(Basis(1:n) * TempDiff(1:n))
        m_volabs = m_volabs + s * sum(Basis(1:n) * abs(SurfaceMove(SurfPerm(NodeIndexes(1:n)))))
        m_tabs = m_tabs + s * sum(Basis(1:n) * abs(TempDiff(1:n)))
      end do
    end do
  end subroutine pca_steady_map

  ! Lines 835-906: normalization, lumped Newton, surface update, reports.
  function pca_steady_update(relax, visit, newton_in) result(nw) bind(c, name="pca_steady_update")
    real(c_double), value :: relax
    integer(c_int), value :: visit, newton_in
    integer(c_int) :: nw
    logical :: stat, newton
    integer :: i, j
    real(kind=dp) :: localrelax, cvol, clim, maxupdate, lnorm
    newton = (newton_in /= 0)
    m_tave = m_tave / m_area
    m_tabs = m_tabs / m_area
    m_volume = m_volume / m_area
    m_volabs = m_volabs / m_area
    i = ListGetInteger(Solver % Values, 'Lumped Newton After Iterations', stat)
    if (stat .and. visit > i) then
      j = ListGetInteger(Solver % Values, 'Lumped Newton Mode', stat)
      select case (j)
      case (1)
        cvol = 0.5d0 * (m_prevtave + m_tave) / (m_prevtave - m_tave)
      case (2)
        cvol = 0.5d0 * (m_prevvolabs + m_volabs) / (m_prevvolabs - m_volabs)
      case (3)
        cvol = 0.5d0 * (m_prevtabs + m_tabs) / (m_prevtabs - m_tabs)
      case default
        cvol = 0.5d0 * (m_prevvolume + m_volume) / (m_prevvolume - m_volume)
      end select
      if (cvol < 0.0d0) then
        cvol = 1.0d0
        m_ccum = 1.0d0
      end if
      clim = ListGetConstReal(Solver % Values, 'Lumped Newton Limiter', stat)
      if (.not. stat) clim = 100.0d0
      cvol = min(clim, cvol)
      cvol = max(1.0d0 / clim, cvol)
      m_ccum = m_ccum * cvol
      write(Message, *) 'Lumped Newton relaxation: ', m_ccum
      call Info('PhaseChangeSolve', Message)
      localrelax = relax * m_ccum
    else
      localrelax = relax
    end if
    SurfaceMove = localrelax * SurfaceMove
    Surface = Surface + SurfaceMove
    m_dpos = SurfaceMove(SurfPerm(Trip_node))
    maxupdate = maxval(abs(SurfaceMove)) / maxval(abs(Surface))
    if (abs(maxupdate) < m_newton_after_tol) newton = .true.
    write(Message, *) 'Maximum surface update: ', maxupdate
    call Info('PhaseChangeSolve', Message)
    write(Message, *) 'Maximum temperature difference: ', m_maxtempdiff
    call Info('PhaseChangeSolve', Message)
    lnorm = sqrt(sum(Surface**2) / size(Surface))
    write(Message, *) 'Result Norm     : ', lnorm
    call Info('PhaseChangeSolve', Message, level=4)
    Solver % Variable % Norm = lnorm
    m_prevvolume = m_volume
    m_prevtave = m_tave
    m_prevtabs = m_tabs
    m_prevvolabs = m_volabs
    if (m_isoalloc) then
      m_isoalloc = .false.
      deallocate(m_isosurf)
    end if
    nw = merge(1, 0, newton)
  end function pca_steady_update

  ! Lines 913-923: label-200 tail (simulation result reports).
  subroutine pca_finish(trip_temp) bind(c, name="pca_finish")
    real(c_double), value :: trip_temp
    call ListAddConstReal(Model % Simulation, 'res: Triple point temperature', trip_temp)
    call ListAddConstReal(Model % Simulation, 'res: triple point movement', m_dpos)
    call ListAddConstReal(Model % Simulation, 'res: Pull Position', m_pos0)
    if (m_pullcontrol) then
      if (NormalDir == 1) then
        call ListAddConstReal(Model % Simulation, 'res: Pull Velocity 1', UPull(NormalDir))
      else
        call ListAddConstReal(Model % Simulation, 'res: Pull Velocity 2', UPull(NormalDir))
      end if
    end if
  end subroutine pca_finish

  ! ==========================================================================
  ! Transcribed CONTAINS procedures (statement for statement; host
  ! association now resolves to the module mirrors above with identical
  ! lifetimes: NormalDir/TangentDirection/Trip_node/AverageNormal/Normals
  ! latched at first time, Density/MeltPoint per call, dt/SpeedUp per
  ! entry, m_useloads per entry).
  ! ==========================================================================

  subroutine VelocityLocalMatrix(StiffMatrix, MassMatrix, ForceVector, &
      Element, nCoord, Nodes)

    ! external variables:
    real(kind=dp) :: StiffMatrix(:,:), MassMatrix(:,:), ForceVector(:)
    integer :: nCoord
    type(Nodes_t) :: Nodes
    type(Element_t), target :: Element

    ! internal variables:
    type(Nodes_t) :: PNodes
    type(Element_t), pointer :: Parent
    real(kind=dp) :: Basis(3*nCoord),dBasisdx(3*nCoord,3), &
        X,Y,Z,U,V,W,S,detJ, TGrad(3,3),Flux,pu,pv,pw,pull, LocalHeat, &
        NodalTemp(3*nCoord), xx(10),yy(10),zz(10),  NodalSurf(nCoord), NodalNormal(2,nCoord)
    real(kind=dp) :: Velo, Ny, Normal(3), StabFactor, StabCoeff, DerSurf, xcoord
    logical :: Stat
    integer :: i,j,k,l,t,p,q, n, pn
    type(GaussIntegrationPoints_t) :: IntegStuff

!------------------------------------------------------------------------------

    ForceVector = 0.0d0
    StiffMatrix = 0.0d0
    MassMatrix  = 0.0d0
    n = nCoord

    StabFactor = ListGetConstReal(Solver % Values,'Velocity Smoothing Factor')

    Nodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
    Nodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
    NodalSurf(1:n) = Surface( SurfPerm(NodeIndexes) )

    if (AverageNormal) then
      NodalNormal(1,1:n) = Normals(2*NormalsPerm(NodeIndexes(1:n))-1)
      NodalNormal(2,1:n) = Normals(2*NormalsPerm(NodeIndexes(1:n)))
    end if

    if (m_transient) then
      allocate(PNodes % x(10), PNodes % y(10), PNodes % z(10))
    end if

!      Numerical integration:
!      ----------------------

    IntegStuff = GaussPoints( Element )

    do t = 1,IntegStuff % n

      u = IntegStuff % u(t)
      v = IntegStuff % v(t)
      w = IntegStuff % w(t)
      s = IntegStuff % s(t)

!        Basis function values & derivatives at the integration point:
!        -------------------------------------------------------------
      stat = ElementInfo( Element,Nodes,U,V,W,detJ,Basis,dBasisdx)
      s = s * detJ
      xcoord = sum( Nodes % x(1:nCoord) * Basis(1:nCoord) )

!      IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
!        s = s * xcoord
!      END IF

      if (AverageNormal) then
        Normal(1) = sum( Basis(1:n) * NodalNormal(1,1:n))
        Normal(2) = sum( Basis(1:n) * NodalNormal(2,1:n))
        Normal(3) = 0.0d0
      else
        Normal = NormalVector( Element, Nodes, u, v, .TRUE. )
      end if

      LocalHeat = sum(Basis(1:n) * LatentHeat(1:n))
      StabCoeff = StabFactor * LocalHeat * m_density

      if (m_useloads) then
        ! do nothing, loads already computed
      else
        ! Compute the flux from normal derivatives
        TGrad = 0.0d0
        l = 0
        do i=1,2

          if ( i == 1) then
            Parent => Element % BoundaryInfo % Left
          else
            Parent => Element % BoundaryInfo % Right
          end if

          k = ListGetInteger(Model % Bodies(Parent % BodyId) % Values,'Material')
          if (ListGetLogical(Model % Materials(k) % Values,'Solid',stat)) then
            if (l == 2) call Fatal('PhaseChangeSolve','Both materials cannot be solid!')
            l = 2
          else
            if (l == 1) call Fatal('PhaseChangeSolve','Both materials cannot be liquid!')
            l = 1
          end if

          pn = Parent % TYPE % NumberOfNodes
          k = ListGetInteger(Model % Bodies(Parent % BodyId) % Values,'Material')

          Conductivity(1:pn) = ListGetReal( Model % Materials(k) % Values, &
              'Heat Conductivity', pn, Parent % NodeIndexes )
          NodalTemp(1:pn) = Temperature( TempPerm(Parent % NodeIndexes) )

          stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis, dBasisdx)

          !           Calculate the basis functions for the parent element:
          !           -----------------------------------------------------
          do j = 1,n
            do k = 1,pn
              if ( NodeIndexes(j) == Parent % NodeIndexes(k) ) then
                xx(j) = Parent % TYPE % NodeU(k)
                yy(j) = Parent % TYPE % NodeV(k)
                zz(j) = Parent % TYPE % NodeW(k)
                exit
              end if
            end do
          end do

          pu = sum( Basis(1:n) * xx(1:n) )
          pv = sum( Basis(1:n) * yy(1:n) )
          pw = sum( Basis(1:n) * zz(1:n) )

          PNodes % x(1:pn) = Solver % Mesh % Nodes % x(Parent % NodeIndexes)
          PNodes % y(1:pn) = Solver % Mesh % Nodes % y(Parent % NodeIndexes)
          PNodes % z(1:pn) = Solver % Mesh % Nodes % z(Parent % NodeIndexes)

          stat = ElementInfo( Parent, PNodes, pu, pv, pw, detJ, Basis, dBasisdx )

          do j=1,m_dim
            TGrad(l,j) = sum( Conductivity(1:pn) * Basis(1:pn) ) * &
                sum( dBasisdx(1:pn,j) * NodalTemp(1:pn) )
          end do
        end do

        Flux = sum( (TGrad(1,:) - TGrad(2,:)) *  Normal)
        stat = ElementInfo( Element,Nodes,U, V, W, detJ, Basis, dBasisdx )
      end if

      ! Assembly the matrix
      do p=1,n
        do q=1,n
          StiffMatrix(p,q) = StiffMatrix(p,q) + s * Basis(p) * Basis(q) * &
              Normal(NormalDir) * m_density * LocalHeat
          if (NodeIndexes(p) /= Trip_node) then
            StiffMatrix(p,q) = StiffMatrix(p,q) + &
                s * StabCoeff * dBasisdx(q,TangentDirection) * dBasisdx(p,TangentDirection)
          end if
        end do

        ! transient part of heat flux
        if (.not. (m_useloads)) then
          ForceVector(p) = ForceVector(p) - s * Basis(p) * Flux
        end if

      end do

    end do

    if (m_transient) deallocate( PNodes % x, PNodes % y, PNodes % z )

!------------------------------------------------------------------------------
  end subroutine VelocityLocalMatrix
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
  subroutine SurfaceLocalMatrix( StiffMatrix, MassMatrix, ForceVector,&
      Element, nCoord, Nodes )

    ! external variables:
    real(kind=dp) :: StiffMatrix(:,:), MassMatrix(:,:), ForceVector(:)
    integer :: nCoord
    type(Nodes_t) :: Nodes
    type(Element_t), pointer :: Element

    ! internal variables:
    real(kind=dp) :: Basis(3*nCoord),dBasisdx(3*nCoord,3), &
        X,Y,Z,U,V,W,S,detJ,pull, NodalVelo(nCoord), xcoord, NodalNormal(2,nCoord)
    real(kind=dp) :: Velo, StabFactor, StabCoeff
    logical :: Stat
    integer :: i,j,k,l,t,p,q, n, pn
    type(GaussIntegrationPoints_t) :: IntegStuff

!------------------------------------------------------------------------------

    ForceVector = 0.0d0
    StiffMatrix = 0.0d0
    MassMatrix  = 0.0d0
    n = nCoord

    StabFactor = ListGetConstReal(Solver % Values,'Surface Smoothing Factor')
    StabCoeff = StabFactor / m_dt

    Nodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
    Nodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
    NodalVelo(1:n) = SurfaceMove( SurfPerm(NodeIndexes) )

    if (AverageNormal) then
      NodalNormal(1,1:n) = Normals(2*NormalsPerm(NodeIndexes(1:n))-1)
      NodalNormal(2,1:n) = Normals(2*NormalsPerm(NodeIndexes(1:n)))
    end if


!      Numerical integration:
!      ----------------------

    IntegStuff = GaussPoints( Element )

    do t = 1,IntegStuff % n

      u = IntegStuff % u(t)
      v = IntegStuff % v(t)
      w = IntegStuff % w(t)
      s = IntegStuff % s(t)

!        Basis function values & derivatives at the integration point:
!        -------------------------------------------------------------
      stat = ElementInfo( Element,Nodes,U,V,W,detJ,Basis,dBasisdx)
      s = s * detJ
      xcoord = sum( Nodes % x(1:nCoord) * Basis(1:nCoord) )

      if (AverageNormal) then
        Normal(1) = sum( Basis(1:n) * NodalNormal(1,1:n) )
        Normal(2) = sum( Basis(1:n) * NodalNormal(2,1:n) )
        Normal(3) = 0.0d0
      else
        Normal = NormalVector( Element, Nodes, u, v, .TRUE. )
      end if

      Velo = sum( Basis(1:n) * NodalVelo(1:n)) + UPull(NormalDir)
      Velo = SpeedUp * Velo

      do p=1,n
        do q=1,n
          MassMatrix(p,q) = MassMatrix(p,q) + s * Basis(p) * Basis(q)
          StiffMatrix(p,q) = StiffMatrix(p,q) + &
              s * StabCoeff * dBasisdx(q,TangentDirection) * dBasisdx(p,TangentDirection)

          ! BC for tipple node
          if (NodeIndexes(p) == Trip_node) then
            StiffMatrix(p,q) = StiffMatrix(p,q) - &
                StabCoeff * Basis(p) * dBasisdx(q,TangentDirection)
          end if
        end do
        ForceVector(p) = ForceVector(p) + s * Basis(p) * Velo
      end do

    end do

!------------------------------------------------------------------------------
  end subroutine SurfaceLocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    subroutine FindPullBoundary()

      real(kind=dp) :: Ybot, Ymax, Ytop, Xtrip, x, y
      integer :: t,k,i,n
      logical :: PullBoundary
      ! Undeclared in the original (implicit REAL under CONTAINS scoping).
      ! Declared LOGICAL here: the value is never read, only passed as the
      ! FOUND dummy, so the 4-byte write lands identically either way.
      logical :: stat

      if (TangentDirection /= 1) call Fatal('PhaseChangeSolve',&
          'FindPullBoundary implemented only for lateral boundaries!')

      PullBoundary = .false.
      Ybot = Solver % Mesh % Nodes % y(Trip_node)
      Ymax = maxval(Solver % Mesh % Nodes % y)
      Ytop = Ymax
      Xtrip = Solver % Mesh % Nodes % x(Trip_node)

      do t = Solver % Mesh % NumberOfBulkElements + 1, &
          Solver % Mesh % NumberOfBulkElements + Solver % Mesh % NumberOfBoundaryElements

        CurrentElement => Solver % Mesh % Elements(t)
        Model % CurrentElement => CurrentElement
        n = CurrentElement % TYPE % NumberOfNodes
        NodeIndexes => CurrentElement % NodeIndexes

        do k=1, Model % NumberOfBCs
          if ( Model % BCs(k) % Tag /= CurrentElement % BoundaryInfo % Constraint ) cycle
          if ( ListGetLogical(Model % BCs(k) % Values,'Pull Boundary',stat ) ) then
            PullBoundary = .true.
            do i = 1,n
              x = Solver % Mesh % Nodes % x(NodeIndexes(i))
              y = Solver % Mesh % Nodes % y(NodeIndexes(i))
              if (y > Ybot .and. abs(x-Xtrip) > 1.0e-6 * (Ymax-Ybot)) Ytop = y
            end do
          end if
        end do
      end do

      if (PullBoundary) then
        call ListAddConstReal(Model % Simulation,'res: Triple point position',Ybot)
        call ListAddConstReal(Model % Simulation,'res: Full pull position',Ytop)
      end if

    end subroutine FindPullBoundary


!------------------------------------------------------------------------------
  ! Init entry: PhaseChangeSolve_Init (lines 43-66).
  subroutine pca_init_begin(h_model, h_solver) bind(c, name="pca_init_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, Model)
    call c_f_pointer(h_solver, Solver)
    Mesh => Solver % Mesh
  end subroutine pca_init_begin

  subroutine pca_init_body() bind(c, name="pca_init_body")
    logical :: found
    character(len=MAX_NAME_LEN) :: variablename
    Params => GetSolverParams()
    variablename = GetString(Params, 'Variable')
    call ListAddString(Params, NextFreeKeyword('Exported Variable ', Params), &
        '-nooutput '//trim(ComponentName(variablename))//'Move')
    if (ListGetLogical(Params, 'Use Average Velocity', found)) &
      call ListAddString(Params, NextFreeKeyword('Exported Variable ', Params), &
          '-nooutput '//trim(ComponentName(variablename))//'MoveAve')
  end subroutine pca_init_body

!------------------------------------------------------------------------------
  ! MeltingHeat co-entry: the Heat Flux UDF called by HeatSolve. The
  ! original keeps SAVE(NormalExist, Normals, NormalsPerm, Nodes) but its
  ! Visited flag is a fresh local every call, so the NormalSol lookup and
  ! the NormalExist/Normals/NormalsPerm (re-)association run on every
  ! call; the Nodes buffer is allocated once. Mirrored exactly: locals
  ! for the per-call state, one module SAVE buffer for Nodes.
  function pca_mh_eval(h_model, node, tval) result(flux) bind(c, name="pca_mh_eval")
    type(c_ptr), value :: h_model
    integer(c_int), value :: node
    real(c_double), value :: tval
    real(c_double) :: flux
    type(Model_t), pointer :: mh_model
    type(Variable_t), pointer :: normalsol
    type(Element_t), pointer :: currentelement, parent
    type(Mesh_t), pointer :: mh_mesh
    integer, pointer :: mh_nodeindexes(:), mh_normalsperm(:)
    real(kind=dp), pointer :: mh_normals(:)
    real(kind=dp) :: upull(3), normal(3), normalpull
    real(kind=dp) :: nodelatentheat, density, elemlatentheat(4)
    real(kind=dp) :: u, v
    integer :: k, n, i
    logical :: stat, normalexist
    real(kind=dp) :: tt
    tt = tval
    call c_f_pointer(h_model, mh_model)
    mh_mesh => mh_model % Mesh
    normalsol => VariableGet(mh_model % Variables, 'Normals', ThisOnly=.true.)
    if (associated(normalsol)) then
      mh_normalsperm => normalsol % Perm
      mh_normals => normalsol % Values
      normalexist = .true.
    else
      ! The original re-runs this lookup on every call (its Visited flag
      ! is a fresh local, so NormalExist is re-set every call): mirror
      ! that by assigning the flag on both paths.
      normalexist = .false.
      n = mh_mesh % MaxElementNodes
      if (.not. mh_nodes_alloc) then
        allocate(mh_nodes % x(n), mh_nodes % y(n), mh_nodes % z(n))
        mh_nodes_alloc = .true.
      end if
    end if
    currentelement => mh_model % CurrentElement
    mh_nodeindexes => currentelement % NodeIndexes
    n = currentelement % TYPE % NumberOfNodes
    k = ListGetInteger(mh_model % Bodies(currentelement % BodyId) % Values, 'Material')
    elemlatentheat(1:n) = ListGetReal(mh_model % Materials(k) % Values, &
        'Latent Heat', n, mh_nodeindexes)
    do i = 1, n
      if (mh_nodeindexes(i) == node) exit
    end do
    if (mh_nodeindexes(i) /= node) call Fatal('PhaseProcs', 'Node not found')
    nodelatentheat = elemlatentheat(i)
    upull = 0.0d0
    upull(1) = ListGetConstReal(mh_model % Simulation, 'res: Pull Velocity 1', stat)
    if (.not. stat) upull(1) = ListGetConstReal(mh_model % Materials(k) % Values, &
        'Convection Velocity 1', stat)
    upull(2) = ListGetConstReal(mh_model % Simulation, 'res: Pull Velocity 2', stat)
    if (.not. stat) upull(2) = ListGetConstReal(mh_model % Materials(k) % Values, &
        'Convection Velocity 2', stat)
    parent => currentelement % BoundaryInfo % Left
    k = ListGetInteger(mh_model % Bodies(parent % BodyId) % Values, 'Material')
    if (ListGetLogical(mh_model % Materials(k) % Values, 'Solid', stat)) then
      density = ListGetConstReal(mh_model % Materials(k) % Values, 'Density')
    else
      parent => currentelement % BoundaryInfo % Right
      k = ListGetInteger(mh_model % Bodies(parent % BodyId) % Values, 'Material')
      density = ListGetConstReal(mh_model % Materials(k) % Values, 'Density')
    end if
    if (.not. associated(parent)) call Fatal('PhaseProcs', 'Parent not found')
    if (normalexist) then
      normal(1) = mh_normals(2 * mh_normalsperm(node) - 1)
      normal(2) = mh_normals(2 * mh_normalsperm(node))
      normal(3) = 0.0d0
    else
      mh_nodes % x(1:n) = mh_model % Nodes % x(mh_nodeindexes)
      mh_nodes % y(1:n) = mh_model % Nodes % y(mh_nodeindexes)
      mh_nodes % z(1:n) = mh_model % Nodes % z(mh_nodeindexes)
      u = 0.0d0
      v = 0.0d0
      normal = NormalVector(currentelement, mh_nodes, u, v, .true.)
    end if
    normalpull = sum(normal * upull)
    flux = nodelatentheat * density * normalpull
  end function pca_mh_eval

!------------------------------------------------------------------------------
end module pca_acc
