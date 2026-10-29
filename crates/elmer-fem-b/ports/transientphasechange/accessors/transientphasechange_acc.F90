! T021-B port of fem/src/modules/TransientPhaseChange.F90: statement-level
! BIND(C) wrappers compiled against the built .mod files of the target
! configure (ADR-057 §3, T020 Poisson pattern).
!
! Every wrapper takes only C-interoperable arguments. All Fortran-object
! traffic stays inside Fortran. SAVE state mirrors the original's SAVE
! statement with identical lifetimes: latched scalars/arrays live as
! module variables (implicitly SAVE); per-call associations are refreshed
! by tpa_begin / the bind wrappers on every entry, exactly where the
! original re-associates its locals. The three CONTAINS procedures
! (VelocityLocalMatrix, SurfaceLocalMatrix, FindPullBoundary) are kept as
! Fortran module procedures transcribed statement for statement, so all
! order-sensitive reductions (SUM/MINVAL/MAXVAL), formatted I/O
! (WRITE ... ES12.5 / I2S) and whole-array updates hold bit-identity by
! construction under the same compiler and flags.
!
! Dead state omitted with zero observable effect (verified against the
! source: assigned but never read, or declared but never referenced):
! PSolver, AllocationsDone, RelativeChange, NonlinearTol, PrevUpull,
! host scalars d/s/u/v/w/yy/ds, and the never-associated Normals pointer.
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module tpc_acc
  use iso_c_binding
  use DefUtils
  implicit none

  ! -- entry context (per-call: refreshed by tpa_begin each entry) --------
  type(Model_t), pointer :: Model => null()
  type(Solver_t), pointer :: Solver => null()
  type(Mesh_t), pointer :: Mesh => null()
  type(ValueList_t), pointer :: Params => null()
  type(ValueList_t), pointer :: Material => null()
  type(ValueList_t), pointer :: ValueList => null()
  type(Variable_t), pointer :: Var => null()
  type(Variable_t), pointer :: SurfSol => null()
  type(Variable_t), pointer :: TempSol => null()
  type(Variable_t), pointer :: HelpSol => null()
  type(Variable_t), pointer :: LoadsSol => null()
  type(Matrix_t), pointer :: StiffMatrix => null()
  real(kind=dp), pointer :: Surface(:), PrevSurface(:), Temperature(:)
  real(kind=dp), pointer :: ForceVector(:)
  integer, pointer :: SurfPerm(:), TempPerm(:), NodeIndexes(:)
  type(Element_t), pointer :: CurrentElement => null()

  ! -- SAVE-latched state (mirrors the SAVE statement lifetimes) ---------
  real(kind=dp), pointer :: SurfaceVelo(:), PrevSurfaceVelo(:)
  real(kind=dp), pointer :: CurrentLoads(:), PrevLoads(:)
  real(kind=dp), pointer :: x(:), y(:), z(:), Basis(:), dBasisdx(:,:)
  real(kind=dp), pointer :: NodalTemp(:), Conductivity(:), LatentHeat(:)
  real(kind=dp), pointer :: Density(:)
  real(kind=dp), pointer :: LocalStiffMatrix(:,:), LocalForceVector(:)
  real(kind=dp), pointer :: LocalMassMatrix(:,:)
  logical, pointer :: IsBoundaryNode(:)
  type(Nodes_t) :: Nodes
  type(Variable_t), pointer :: NrmSol => null()
  logical :: AverageNormal = .false.
  integer :: NormalDir = 0
  ! Per-call mirrors for the SurfaceLocalMatrix host association.
  real(kind=dp) :: m_speedup = 0.0_dp, m_upull = 0.0_dp, m_dt = 0.0_dp
  real(kind=dp) :: m_norm = 0.0_dp
  integer :: m_dim = 0
  logical :: m_useloads = .false.
  character(len=MAX_NAME_LEN) :: TemperatureName = ''

contains

  ! -- entry -----------------------------------------------------------
  subroutine tpa_begin(h_model, h_solver) bind(c, name="tpa_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, Model)
    call c_f_pointer(h_solver, Solver)
    Mesh => Solver % Mesh
    Params => GetSolverParams()
  end subroutine tpa_begin

  subroutine tpa_info_header() bind(c, name="tpa_info_header")
    call Info('TransientPhaseChange','--------------------------------------------')
    call Info('TransientPhaseChange','Using transient algorithm for surface update')
    call Info('TransientPhaseChange','--------------------------------------------')
  end subroutine tpa_info_header

  subroutine tpa_fatal_notransient() bind(c, name="tpa_fatal_notransient")
    call Fatal('TransientPhaseChange','This only makes sense in a transient setting')
  end subroutine tpa_fatal_notransient

  function tpa_coord_dim() result(d) bind(c, name="tpa_coord_dim")
    integer(c_int) :: d
    m_dim = CoordinateSystemDimension()
    d = m_dim
  end function tpa_coord_dim

  ! -- primary / temperature bindings -----------------------------------
  subroutine tpa_bind_primary() bind(c, name="tpa_bind_primary")
    SurfSol => Solver % Variable
    Surface => SurfSol % Values
    PrevSurface => SurfSol % PrevValues(:,1)
    SurfPerm => SurfSol % Perm
    if (.not. associated(Surface) .or. all(SurfPerm <= 0)) then
      call Fatal('TransientPhaseChange','Surface field needed for Phase Change')
    end if
  end subroutine tpa_bind_primary

  subroutine tpa_bind_temp() bind(c, name="tpa_bind_temp")
    logical :: stat
    TemperatureName = ListGetString(Params, 'Phase Change Variable', stat)
    if (.not. stat) TemperatureName = 'Temperature'
    TempSol => VariableGet(Mesh % Variables, trim(TemperatureName))
    TempPerm => TempSol % Perm
    Temperature => TempSol % Values
    if (.not. associated(Temperature) .or. all(TempPerm <= 0)) then
      call Fatal('TransientPhaseChange','Temperature field needed for Phase Change')
    end if
  end subroutine tpa_bind_temp

  ! -- solver params ------------------------------------------------------
  function tpa_nonlin_iter() result(n) bind(c, name="tpa_nonlin_iter")
    integer(c_int) :: n
    logical :: stat
    n = ListGetInteger(Params, 'Nonlinear System Max Iterations', stat)
    if (.not. stat) n = 1
  end function tpa_nonlin_iter

  subroutine tpa_nonlintol_touch() bind(c, name="tpa_nonlintol_touch")
    logical :: stat
    real(kind=dp) :: tol
    tol = ListGetConstReal(Params, 'Nonlinear System Convergence Tolerance', stat)
  end subroutine tpa_nonlintol_touch

  function tpa_pullcontrol(dim) result(v) bind(c, name="tpa_pullcontrol")
    integer(c_int), value :: dim
    integer(c_int) :: v
    logical :: stat, pc
    pc = ListGetLogical(Params, 'Pull Rate Control', stat)
    if (pc) then
      if (dim == 3) then
        call Fatal('TransientPhaseChange','Pull rate control not implemented in 3D')
      end if
      call Info('TransientPhaseChange','Using pull control for the phase change', level=7)
    end if
    v = merge(1, 0, pc)
  end function tpa_pullcontrol

  function tpa_triplefixed(dim) result(v) bind(c, name="tpa_triplefixed")
    integer(c_int), value :: dim
    integer(c_int) :: v
    logical :: stat, tf
    tf = ListGetLogical(Params, 'Triple Point Fixed', stat)
    if (tf) then
      if (dim == 3) then
        call Fatal('TransientPhaseChange','Fixed triple point position not implemented in 3D')
      end if
      call Info('TransientPhaseChange','Fixing triple point position', level=7)
    end if
    v = merge(1, 0, tf)
  end function tpa_triplefixed

  function tpa_coupled_iter() result(n) bind(c, name="tpa_coupled_iter")
    integer(c_int) :: n
    HelpSol => VariableGet(Mesh % Variables, 'coupled iter')
    n = nint(HelpSol % Values(1))
  end function tpa_coupled_iter

  function tpa_timestep() result(n) bind(c, name="tpa_timestep")
    integer(c_int) :: n
    HelpSol => VariableGet(Mesh % Variables, 'timestep')
    n = nint(HelpSol % Values(1))
  end function tpa_timestep

  ! -- first-time initializations ------------------------------------------
  function tpa_nof_nodes() result(n) bind(c, name="tpa_nof_nodes")
    integer(c_int) :: n
    n = Model % Mesh % NumberOfNodes
  end function tpa_nof_nodes

  function tpa_surfperm(node) result(p) bind(c, name="tpa_surfperm")
    integer(c_int), value :: node
    integer(c_int) :: p
    p = SurfPerm(node)
  end function tpa_surfperm

  function tpa_nodecoord(node, j1) result(c) bind(c, name="tpa_nodecoord")
    integer(c_int), value :: node, j1
    real(c_double) :: c
    if (j1 == 1) then
      c = Model % Mesh % Nodes % x(node)
    else if (j1 == 2) then
      c = Model % Mesh % Nodes % y(node)
    else
      c = Model % Mesh % Nodes % z(node)
    end if
  end function tpa_nodecoord

  function tpa_normaldir(jdet) result(nd) bind(c, name="tpa_normaldir")
    integer(c_int), value :: jdet
    integer(c_int) :: nd
    logical :: stat
    NormalDir = ListGetInteger(Params, 'Normal Direction', stat)
    if (.not. stat) then
      NormalDir = jdet
      call Info('TransientPhaseChange','Normal coordinate set to: '//I2S(NormalDir), level=7)
    end if
    nd = NormalDir
  end function tpa_normaldir

  subroutine tpa_alloc_bndflags() bind(c, name="tpa_alloc_bndflags")
    allocate(IsBoundaryNode(size(Surface)))
    IsBoundaryNode = .false.
  end subroutine tpa_alloc_bndflags

  subroutine tpa_set_bndflag(idx) bind(c, name="tpa_set_bndflag")
    integer(c_int), value :: idx
    IsBoundaryNode(idx) = .true.
  end subroutine tpa_set_bndflag

  subroutine tpa_bnd_scan(t, n, has) bind(c, name="tpa_bnd_scan")
    integer(c_int), value :: t
    integer(c_int) :: n, has
    logical :: stat, side
    CurrentElement => GetBoundaryElement(t)
    n = GetElementNOFNodes()
    NodeIndexes => CurrentElement % NodeIndexes
    ValueList => GetBC()
    side = ListGetLogical(ValueList, 'Phase Change Side', stat)
    if (.not. stat) side = .false.
    has = merge(1, 0, side)
  end subroutine tpa_bnd_scan

  function tpa_bnd_nodeidx(i) result(node) bind(c, name="tpa_bnd_nodeidx")
    integer(c_int), value :: i
    integer(c_int) :: node
    node = NodeIndexes(i)
  end function tpa_bnd_nodeidx

  function tpa_num_bnd() result(n) bind(c, name="tpa_num_bnd")
    integer(c_int) :: n
    n = Mesh % NumberOfBoundaryElements
  end function tpa_num_bnd

  subroutine tpa_report_bndcount() bind(c, name="tpa_report_bndcount")
    integer :: n
    n = count(IsBoundaryNode)
    call Info('TransientPhaseChange','Number of boundary nodes: '//I2S(n), level=7)
  end subroutine tpa_report_bndcount

  subroutine tpa_alloc_first() bind(c, name="tpa_alloc_first")
    integer :: n, istat
    n = Mesh % MaxElementNodes
    allocate(Nodes % x(n), Nodes % y(n), Nodes % z(n), &
        x(n), y(n), z(n), Basis(n), dBasisdx(n,3), NodalTemp(n), &
        Conductivity(n), LatentHeat(n), Density(n), &
        LocalStiffMatrix(n,n), LocalForceVector(n), LocalMassMatrix(n,n), &
        stat=istat)
    if (istat /= 0) call Fatal('TransientPhaseChange', 'Memory allocation error 1.')
  end subroutine tpa_alloc_first

  subroutine tpa_bind_normals(nliter) bind(c, name="tpa_bind_normals")
    integer(c_int), value :: nliter
    logical :: stat
    character(len=MAX_NAME_LEN) :: vname
    vname = ListGetString(Params, 'Normal Variable', stat)
    if (stat) then
      NrmSol => VariableGet(Mesh % Variables, trim(vname), ThisOnly=.true.)
    else
      NrmSol => VariableGet(Mesh % Variables, 'Normals', ThisOnly=.true.)
    end if
    AverageNormal = associated(NrmSol)
    if (AverageNormal .and. nliter > 1) then
      call Warn('TransientPhaseChange','With external normal field there is no nonlinearity to iterate!')
    end if
  end subroutine tpa_bind_normals

  subroutine tpa_bind_velo() bind(c, name="tpa_bind_velo")
    character(len=MAX_NAME_LEN) :: vname
    vname = trim(ComponentName(Solver % Variable))//' Velo'
    HelpSol => VariableGet(Mesh % Variables, trim(vname))
    if (.not. associated(HelpSol)) then
      call Fatal('TransientPhaseChange','Surface Velo field should exist!')
    end if
    SurfaceVelo => HelpSol % Values
  end subroutine tpa_bind_velo

  subroutine tpa_velrelax(val, found) bind(c, name="tpa_velrelax")
    real(c_double) :: val
    integer(c_int) :: found
    logical :: f
    val = GetCReal(Params, 'Velocity Relaxation Factor', f)
    found = merge(1, 0, f)
  end subroutine tpa_velrelax

  subroutine tpa_alloc_prevvelo() bind(c, name="tpa_alloc_prevvelo")
    integer :: istat
    allocate(PrevSurfaceVelo(size(Surface)), stat=istat)
    if (istat /= 0) call Fatal('TransientPhaseChange', 'Memory allocation error 2.')
    PrevSurfaceVelo = 0.0_dp
  end subroutine tpa_alloc_prevvelo

  ! -- passive steps / pull velocity / loads --------------------------------
  subroutine tpa_passive_steps(val, found) bind(c, name="tpa_passive_steps")
    integer(c_int) :: val, found
    logical :: f
    val = ListGetInteger(Params, 'Passive Steps', f)
    found = merge(1, 0, f)
  end subroutine tpa_passive_steps

  function tpa_pull_velo() result(v) bind(c, name="tpa_pull_velo")
    real(c_double) :: v
    logical :: stat
    integer :: k
    character(len=MAX_NAME_LEN) :: str
    CurrentElement => GetActiveElement(1)
    k = GetMaterialId()
    Material => Model % Materials(k) % Values
    write (str,'(A,I2)') 'Convection Velocity',NormalDir
    v = ListGetConstReal(Material, str, stat)
  end function tpa_pull_velo

  function tpa_use_loads() result(v) bind(c, name="tpa_use_loads")
    integer(c_int) :: v
    logical :: stat
    v = merge(1, 0, ListGetLogical(Params, 'Use Nodal Loads', stat))
  end function tpa_use_loads

  subroutine tpa_bind_loads(sign, usefirst) bind(c, name="tpa_bind_loads")
    integer(c_int) :: sign, usefirst
    logical :: stat
    integer :: istat
    LoadsSol => VariableGet(Mesh % Variables, trim(TemperatureName)//' Loads')
    if (.not. associated(LoadsSol)) then
      call Fatal('TransientPhaseChange','Loads are requested to be used but missing!')
    else
      call Info('TransientPhaseChange','Using nodal loads to determine phase change', level=6)
    end if
    if (ListGetLogical(Params, 'Nodal Loads Negative', stat)) then
      call Info('TransientPhaseChange','Assigning negative sign to nodal loads!', level=6)
      sign = -1
    else
      sign = 1
    end if
    if (.not. associated(PrevLoads)) then
      call Info('TransientPhaseChange','Allocating for previous loads', level=12)
      allocate(PrevLoads(size(Surface)), stat=istat)
      PrevLoads = 0.0_dp
    end if
    if (.not. associated(CurrentLoads)) then
      call Info('TransientPhaseChange','Allocating for current loads', level=12)
      allocate(CurrentLoads(size(Surface)), stat=istat)
      CurrentLoads = 0.0_dp
    end if
    if (ListGetLogical(Params, 'Use First Loads', stat)) then
      usefirst = 1
    else
      usefirst = 0
    end if
  end subroutine tpa_bind_loads

  subroutine tpa_loads_snapshot() bind(c, name="tpa_loads_snapshot")
    PrevLoads = CurrentLoads
  end subroutine tpa_loads_snapshot

  function tpa_loadsperm(node) result(p) bind(c, name="tpa_loadsperm")
    integer(c_int), value :: node
    integer(c_int) :: p
    p = LoadsSol % Perm(node)
  end function tpa_loadsperm

  function tpa_loadsval(j) result(v) bind(c, name="tpa_loadsval")
    integer(c_int), value :: j
    real(c_double) :: v
    v = LoadsSol % Values(j)
  end function tpa_loadsval

  subroutine tpa_curloads_set(i, v) bind(c, name="tpa_curloads_set")
    integer(c_int), value :: i
    real(c_double), value :: v
    CurrentLoads(i) = v
  end subroutine tpa_curloads_set

  subroutine tpa_report_loads() bind(c, name="tpa_report_loads")
    write(Message,'(A,ES12.5)') 'Minimum nodal load at interface: ', minval(CurrentLoads)
    call Info('TransientPhaseChange', Message)
    write(Message,'(A,ES12.5)') 'Maximum nodal load at interface: ', maxval(CurrentLoads)
    call Info('TransientPhaseChange', Message)
  end subroutine tpa_report_loads

  subroutine tpa_snap_prevvelo() bind(c, name="tpa_snap_prevvelo")
    PrevSurfaceVelo = SurfaceVelo
  end subroutine tpa_snap_prevvelo

  subroutine tpa_assoc_mat() bind(c, name="tpa_assoc_mat")
    StiffMatrix => Solver % Matrix
    ForceVector => Solver % Matrix % RHS
  end subroutine tpa_assoc_mat

  subroutine tpa_loadsrelax(val, useavg) bind(c, name="tpa_loadsrelax")
    real(c_double) :: val
    integer(c_int) :: useavg
    logical :: f
    val = GetCReal(Params, 'Loads Relaxation Factor', f)
    useavg = merge(1, 0, f)
  end subroutine tpa_loadsrelax

  subroutine tpa_fill_force(order, relax) bind(c, name="tpa_fill_force")
    integer(c_int), value :: order
    real(c_double), value :: relax
    integer :: t, i
    do t = 1, Mesh % NumberOfNodes
      i = SurfPerm(t)
      if (i <= 0) cycle
      if (order == 2) then
        ForceVector(i) = relax * CurrentLoads(i) + (1.0_dp - relax) * PrevLoads(i)
      else if (order == 1) then
        ForceVector(i) = relax * CurrentLoads(i) + (1.0_dp - relax) * PrevLoads(i)
      else
        ForceVector(i) = CurrentLoads(i)
      end if
    end do
  end subroutine tpa_fill_force

  ! -- assembly drivers ------------------------------------------------------
  subroutine tpa_default_init() bind(c, name="tpa_default_init")
    call DefaultInitialize()
  end subroutine tpa_default_init

  function tpa_nactive() result(n) bind(c, name="tpa_nactive")
    integer(c_int) :: n
    n = Solver % NumberOfActiveElements
  end function tpa_nactive

  subroutine tpa_ebegin(t, n, nd) bind(c, name="tpa_ebegin")
    integer(c_int), value :: t
    integer(c_int) :: n, nd
    integer :: k
    CurrentElement => GetActiveElement(t)
    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
    k = GetMaterialId()
    Material => Model % Materials(k) % Values
    NodeIndexes => CurrentElement % NodeIndexes
  end subroutine tpa_ebegin

  subroutine tpa_egather(n) bind(c, name="tpa_egather")
    integer(c_int), value :: n
    LatentHeat(1:n) = ListGetReal(Material, 'Latent Heat', n, NodeIndexes)
    Density(1:n) = ListGetReal(Material, 'Density', n, NodeIndexes)
    Conductivity(1:n) = ListGetReal(Material, 'Heat Conductivity', n, NodeIndexes)
    call GetElementNodes(Nodes)
  end subroutine tpa_egather

  subroutine tpa_evelo(n, nd, ul) bind(c, name="tpa_evelo")
    integer(c_int), value :: n, nd, ul
    m_useloads = (ul /= 0)
    call VelocityLocalMatrix(LocalStiffMatrix, LocalMassMatrix, LocalForceVector, &
        CurrentElement, n, nd, Nodes)
  end subroutine tpa_evelo

  subroutine tpa_eupdate() bind(c, name="tpa_eupdate")
    call DefaultUpdateEquations(LocalStiffMatrix, LocalForceVector)
  end subroutine tpa_eupdate

  subroutine tpa_finish_asm() bind(c, name="tpa_finish_asm")
    call DefaultFinishBulkAssembly()
    call DefaultFinishBoundaryAssembly()
    call DefaultFinishAssembly()
  end subroutine tpa_finish_asm

  subroutine tpa_dirichlet() bind(c, name="tpa_dirichlet")
    call DefaultDirichletBCs()
  end subroutine tpa_dirichlet

  function tpa_solve_velo() result(nrm) bind(c, name="tpa_solve_velo")
    real(c_double) :: nrm
    real(kind=dp) :: norm
    call SolveSystem(StiffMatrix, ParMatrix, ForceVector, SurfaceVelo, norm, 1, Solver)
    nrm = norm
  end function tpa_solve_velo

  subroutine tpa_nonlin(change, conv) bind(c, name="tpa_nonlin")
    real(c_double) :: change
    integer(c_int) :: conv
    change = Solver % Variable % NonlinChange
    conv = Solver % Variable % NonlinConverged
  end subroutine tpa_nonlin

  function tpa_velrelax2() result(v) bind(c, name="tpa_velrelax2")
    real(c_double) :: v
    logical :: stat
    v = GetCReal(Params, 'Velocity Relaxation Factor', stat)
  end function tpa_velrelax2

  subroutine tpa_relax_velo(vr) bind(c, name="tpa_relax_velo")
    real(c_double), value :: vr
    SurfaceVelo = vr * SurfaceVelo + (1.0_dp - vr) * PrevSurfaceVelo
  end subroutine tpa_relax_velo

  subroutine tpa_report_velo() bind(c, name="tpa_report_velo")
    write(Message,'(A,ES12.5)') 'Minimum surface velocity: ', minval(SurfaceVelo)
    call Info('TransientPhaseChange', Message, level=12)
    write(Message,'(A,ES12.5)') 'Maximum surface velocity: ', maxval(SurfaceVelo)
    call Info('TransientPhaseChange', Message, level=12)
  end subroutine tpa_report_velo

  function tpa_trip_velo(trip) result(v) bind(c, name="tpa_trip_velo")
    integer(c_int), value :: trip
    real(c_double) :: v
    v = -SurfaceVelo(SurfPerm(trip))
  end function tpa_trip_velo

  subroutine tpa_info_pull(upull) bind(c, name="tpa_info_pull")
    real(c_double), value :: upull
    write(Message,'(A,ES12.5)') 'Pull velocity: ', upull
    call Info('TransientPhaseChange', Message)
  end subroutine tpa_info_pull

  ! -- displacement update ----------------------------------------------------
  subroutine tpa_speedup(val, found) bind(c, name="tpa_speedup")
    real(c_double) :: val
    integer(c_int) :: found
    logical :: f
    val = ListGetConstReal(Solver % Values, 'Transient SpeedUp', f)
    found = merge(1, 0, f)
  end subroutine tpa_speedup

  subroutine tpa_disprelax(val, stat) bind(c, name="tpa_disprelax")
    real(c_double) :: val
    integer(c_int) :: stat
    logical :: f
    val = ListGetCReal(Params, 'Surface Smoothing Factor', f)
    stat = merge(1, 0, f)
  end subroutine tpa_disprelax

  function tpa_solver_order() result(o) bind(c, name="tpa_solver_order")
    integer(c_int) :: o
    o = Solver % Order
  end function tpa_solver_order

  subroutine tpa_simple_update(spdup, dtv, upullv) bind(c, name="tpa_simple_update")
    real(c_double), value :: spdup, dtv, upullv
    call Info('TransientPhaseChange','Updating surface by simple update', level=7)
    Surface = PrevSurface + spdup * dtv * (SurfaceVelo + upullv)
  end subroutine tpa_simple_update

  subroutine tpa_pde_info() bind(c, name="tpa_pde_info")
    call Info('TransientPhaseChange','Updating surfface by solving a PDE', level=7)
  end subroutine tpa_pde_info

  subroutine tpa_sbegin(t, n, nd) bind(c, name="tpa_sbegin")
    integer(c_int), value :: t
    integer(c_int) :: n, nd
    CurrentElement => GetActiveElement(t)
    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
    NodeIndexes => CurrentElement % NodeIndexes
  end subroutine tpa_sbegin

  subroutine tpa_sgather() bind(c, name="tpa_sgather")
    call GetElementNodes(Nodes)
  end subroutine tpa_sgather

  subroutine tpa_surf_compute(nc, spdup, upullv, dtv) bind(c, name="tpa_surf_compute")
    integer(c_int), value :: nc
    real(c_double), value :: spdup, upullv, dtv
    m_speedup = spdup
    m_upull = upullv
    m_dt = dtv
    call SurfaceLocalMatrix(LocalStiffMatrix, LocalMassMatrix, LocalForceVector, &
        CurrentElement, nc, Nodes)
  end subroutine tpa_surf_compute

  subroutine tpa_first_order() bind(c, name="tpa_first_order")
    call Default1stOrderTime(LocalMassMatrix, LocalStiffMatrix, LocalForceVector)
  end subroutine tpa_first_order

  function tpa_solve_surf() result(nrm) bind(c, name="tpa_solve_surf")
    real(c_double) :: nrm
    real(kind=dp) :: norm
    call SolveSystem(StiffMatrix, ParMatrix, ForceVector, Surface, norm, 1, Solver)
    nrm = norm
  end function tpa_solve_surf

  ! -- pull position bookkeeping -----------------------------------------------
  function tpa_donetime() result(n) bind(c, name="tpa_donetime")
    integer(c_int) :: n
    n = Solver % DoneTime
  end function tpa_donetime

  subroutine tpa_find_pull(trip) bind(c, name="tpa_find_pull")
    integer(c_int), value :: trip
    integer :: Trip_save
    Trip_save = trip
    call FindPullBoundary(Trip_save)
  end subroutine tpa_find_pull

  subroutine tpa_res_pull(pos0v, upullv, dopull) bind(c, name="tpa_res_pull")
    real(c_double), value :: pos0v, upullv
    integer(c_int), value :: dopull
    call ListAddConstReal(Model % Simulation, 'res: Pull Position', pos0v)
    call ListAddConstReal(Model % Simulation, 'res: Pull Velocity', upullv)
    if (dopull /= 0) then
      Var => VariableGet(Mesh % Variables, 'pull velocity')
      Var % Values(1) = upullv
      Var => VariableGet(Mesh % Variables, 'pull position')
      Var % Values(1) = pos0v
    end if
  end subroutine tpa_res_pull

  ! -- init entry ----------------------------------------------------------------
  subroutine tpa_init_begin(h_model, h_solver) bind(c, name="tpa_init_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, Model)
    call c_f_pointer(h_solver, Solver)
    Mesh => Solver % Mesh
  end subroutine tpa_init_begin

  subroutine tpa_init_body() bind(c, name="tpa_init_body")
    logical :: found
    character(len=MAX_NAME_LEN) :: vname
    Params => GetSolverParams()
    vname = GetString(Params, 'Variable')
    call ListAddString(Params, NextFreeKeyword('Exported Variable ', Params), &
        '-nooutput '//trim(ComponentName(vname))//' Velo')
    call ListAddNewInteger(Params, 'Time Derivative Order', 1)
    if (ListGetLogical(Params, 'Pull Rate Control', found)) then
      call ListAddString(Params, NextFreeKeyword('Exported Variable ', Params), &
          '-global pull velocity')
      call ListAddString(Params, NextFreeKeyword('Exported Variable ', Params), &
          '-global pull position')
    end if
  end subroutine tpa_init_body

  ! ==========================================================================
  ! Transcribed CONTAINS procedures (statement for statement; host
  ! association now resolves to the module mirrors above with identical
  ! lifetimes: NormalDir/AverageNormal/NrmSol latched at first time,
  ! UseLoads/SpeedUp/Upull/dt refreshed per entry, m_dim per entry).
  ! ==========================================================================

  subroutine VelocityLocalMatrix(StiffMatrix, MassMatrix, ForceVector, &
      Element, n, nd, Nodes)

    ! external variables:
    real(kind=dp) :: StiffMatrix(:,:), MassMatrix(:,:), ForceVector(:)
    integer :: n, nd
    type(Nodes_t) :: Nodes
    type(Element_t), target :: Element

    ! internal variables:
    type(Nodes_t) :: PNodes
    type(Element_t), pointer :: Parent
    real(kind=dp) :: Basis(3*n),dBasisdx(3*n,3), &
        U,V,W,S,detJ, TGrad(3,3),Flux,pu,pv,pw,LocalHeat, LocalDens, &
        NodalTemp(3*n), xx(10),yy(10),zz(10),  NodalSurf(n)
    real(kind=dp) :: Normal(3), StabFactor, StabCoeff, xcoord
    logical :: Stat, Found
    integer :: i,j,k,l,t,p,q, pn, NBasis
    type(GaussIntegrationPoints_t) :: IntegStuff

    !------------------------------------------------------------------------------

    ForceVector = 0.0d0
    StiffMatrix = 0.0d0
    MassMatrix  = 0.0d0

    StabFactor = ListGetConstReal(Params,'Velocity Smoothing Factor',Stat )

    NodalSurf(1:n) = Surface( SurfPerm(NodeIndexes) )

    allocate( PNodes % x(10), PNodes % y(10), PNodes % z(10) )

    !      Numerical integration:
    !      ----------------------

    NBasis = nd
    IntegStuff = GaussPoints( Element )

    do t = 1,IntegStuff % n

      u = IntegStuff % u(t)
      v = IntegStuff % v(t)
      w = IntegStuff % w(t)
      s = IntegStuff % s(t)


      !        Basis function values & derivatives at the integration point:
      !        -------------------------------------------------------------
      stat = ElementInfo( Element,Nodes,U,V,W,detJ,Basis,dBasisdx )
      s = s * detJ
      xcoord = sum( Nodes % x(1:n) * Basis(1:n) )

      ! This is not really a physical equation and hence weighing with the radius is not necessary
      ! However, do it when using loads since the loads have the weights in!
      if ( m_useloads ) then
        if ( CurrentCoordinateSystem() /= Cartesian ) then
          s = s * xcoord
        end if
      end if

      Found = .false.
      if (AverageNormal) then
        Normal = ConsistentNormalVector( CurrentModel % Solver, NrmSol, Element, Found, Basis = Basis )
      end if
      if (.not. Found) then
        Normal = NormalVector( Element, Nodes, u, v, .true. )
      end if

      LocalHeat = sum(Basis(1:n) * LatentHeat(1:n))
      LocalDens = sum(Basis(1:n) * Density(1:n))
      StabCoeff = StabFactor * LocalHeat * LocalDens

      if (m_useloads) then
        ! do nothing, loads already inserted

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
            if (l == 2) call Fatal('TransientPhaseChange','Both materials cannot be solid!')
            l = 2
          else
            if (l == 1) call Fatal('TransientPhaseChange','Both materials cannot be liquid!')
            l = 1
          end if

          pn = Parent % TYPE % NumberOfNodes
          k = ListGetInteger(Model % Bodies(Parent % BodyId) % Values,'Material')

          Conductivity(1:pn) = ListGetReal( Model % Materials(k) % Values, &
              'Heat Conductivity', pn, Parent % NodeIndexes )
          NodalTemp(1:pn) = Temperature( TempPerm(Parent % NodeIndexes) )

          stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis, dBasisdx)

          ! Calculate the basis functions for the parent element:
          !-----------------------------------------------------
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

          PNodes % x(1:pn) = Mesh % Nodes % x(Parent % NodeIndexes)
          PNodes % y(1:pn) = Mesh % Nodes % y(Parent % NodeIndexes)
          PNodes % z(1:pn) = Mesh % Nodes % z(Parent % NodeIndexes)

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
      do p=1,NBasis

        do q=1,NBasis
          StiffMatrix(p,q) = StiffMatrix(p,q) + s * Basis(p) * Basis(q) * &
              Normal(NormalDir) * LocalDens * LocalHeat

          ! Jump over the special nodes (when defined)
          if ( IsBoundaryNode( SurfPerm(NodeIndexes(p))) ) cycle

          do i=1,m_dim
            if ( i == NormalDir ) cycle
            StiffMatrix(p,q) = StiffMatrix(p,q) + &
                s * StabCoeff * dBasisdx(q,i) * dBasisdx(p,i)
          end do
        end do

        ! transient part of heat flux
        if (.not. m_useloads) then
          ForceVector(p) = ForceVector(p) - s * Basis(p) * Flux
        end if


      end do

    end do

    deallocate( PNodes % x, PNodes % y, PNodes % z )

    !------------------------------------------------------------------------------
  end subroutine VelocityLocalMatrix
  !------------------------------------------------------------------------------


  !------------------------------------------------------------------------------
  subroutine SurfaceLocalMatrix(StiffMatrix, MassMatrix, ForceVector, &
      Element, nCoord, Nodes)

    ! external variables:
    real(kind=dp) :: StiffMatrix(:,:), MassMatrix(:,:), ForceVector(:)
    integer :: nCoord
    type(Nodes_t) :: Nodes
    type(Element_t), target :: Element

    ! internal variables:
    real(kind=dp) :: Basis(3*nCoord),dBasisdx(3*nCoord,3), &
        U,V,W,S,detJ,NodalVelo(nCoord), xcoord
    real(kind=dp) :: Velo, StabFactor, StabCoeff
    logical :: Stat
    integer :: i,j,k,l,t,p,q, n
    type(GaussIntegrationPoints_t) :: IntegStuff

    !------------------------------------------------------------------------------

    ForceVector = 0.0d0
    StiffMatrix = 0.0d0
    MassMatrix  = 0.0d0
    n = nCoord

    StabFactor = ListGetConstReal(Params,'Surface Smoothing Factor')
    StabCoeff = StabFactor / m_dt

    NodalVelo(1:n) = SurfaceVelo( SurfPerm(NodeIndexes) )


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

      Velo = m_speedup * ( sum( Basis(1:n) * NodalVelo(1:n)) + m_upull )

      do p=1,n
        do q=1,n
          MassMatrix(p,q) = MassMatrix(p,q) + s * Basis(p) * Basis(q)

          ! Jump over the special nodes (when defined)
          if ( IsBoundaryNode( SurfPerm(NodeIndexes(p))) ) cycle

          do i=1,m_dim
            if ( NormalDir == i ) cycle
            StiffMatrix(p,q) = StiffMatrix(p,q) + &
                s * StabCoeff * dBasisdx(q,i) * dBasisdx(p,i)
          end do
        end do
        ForceVector(p) = ForceVector(p) + s * Basis(p) * Velo
      end do

    end do

    !------------------------------------------------------------------------------
  end subroutine SurfaceLocalMatrix
  !------------------------------------------------------------------------------


  !------------------------------------------------------------------------------
  subroutine FindPullBoundary(Trip_arg)

    integer :: Trip_arg

    real(kind=dp) :: Ybot, Ymax, Ytop, Xtrip, x, y
    integer :: t,k,i,n
    logical :: PullBoundary
    ! Undeclared in the original (implicit REAL under CONTAINS scoping).
    ! Declared LOGICAL here: the value is never read, only passed as the
    ! FOUND dummy, so the 4-byte write lands identically either way.
    logical :: stat

    if (NormalDir /= 2) call Fatal('TransientPhaseChange', &
        'FindPullBoundary implemented only for lateral boundaries!')

    PullBoundary = .false.
    Ybot = Mesh % Nodes % y(Trip_arg)
    Ymax = maxval(Mesh % Nodes % y)
    Ytop = Ymax
    Xtrip = Mesh % Nodes % x(Trip_arg)

    do t = Mesh % NumberOfBulkElements + 1, &
        Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

      CurrentElement => Mesh % Elements(t)
      Model % CurrentElement => CurrentElement
      n = CurrentElement % TYPE % NumberOfNodes
      NodeIndexes => CurrentElement % NodeIndexes

      do k=1, Model % NumberOfBCs
        if ( Model % BCs(k) % Tag /= CurrentElement % BoundaryInfo % Constraint ) cycle
        if ( ListGetLogical(Model % BCs(k) % Values,'Pull Boundary',stat ) ) then
          PullBoundary = .true.
          do i = 1,n
            x = Mesh % Nodes % x(NodeIndexes(i))
            y = Mesh % Nodes % y(NodeIndexes(i))
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
end module tpc_acc
