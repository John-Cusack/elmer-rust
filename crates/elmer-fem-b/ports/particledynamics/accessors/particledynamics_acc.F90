! Port of fem/src/modules/ParticleDynamics.F90 (Elmer @ 9f6af2f85):
! statement-level BIND(C) wrappers for the Elmer calls in
! ParticleDynamics.F90, compiled against the built .mod files of the
! target configure. Imitates the T020 Poisson port shape (ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays).
! All Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER
! dummies, assumed-shape descriptors, hidden lengths, array-pointer
! results — stays inside Fortran, so the Rust side is pure C calls.
! The wrapper module mirrors the original's SAVE state (per-entry
! groups with identical lifetimes): first-visit sizing/lookup results
! are cached wrapper-side, consumed before the next fill, as in the
! original.
!
! Order-sensitive reductions (SUM/MATMUL), formatted I/O (PRINT/WRITE)
! and whole-array updates stay in Fortran wrappers so bit-identity
! holds by construction where the call sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer
! library calls, no upstream text.
module particledynamics_acc
  use iso_c_binding
  use types
  use defutils
  use interpolation
  use meshbasics
  use elementutils
  use particleutils
  use lists
  use messages
  use modeldescription
  use parallelutils
  use loadmod
  use generalutils
  implicit none

  ! -- shared entry context (mirrors main-entry locals) -----------------
  type(solver_t), pointer, save :: acc_solver => null()
  type(valuelist_t), pointer, save :: acc_params => null()
  type(mesh_t), pointer, save :: acc_mesh => null()
  type(particle_t), pointer, save :: acc_particles => null()
  integer, save :: acc_dim = 0
  type(valuelistentry_t), pointer, save :: acc_cond_ptr => null()

  ! -- ParticleFieldInteraction SAVE mirror ----------------------------
  type(mesh_t), pointer, save :: pfi_mesh => null()
  integer, save :: pfi_visited = 0
  integer, save :: pfi_dim = 0, pfi_timeorder = 0, pfi_maxfield = 0, pfi_nogroups = 0
  real(kind=dp), allocatable, save :: pfi_basis(:)
  real(kind=dp), allocatable, save :: pfi_dbasisdx(:,:)
  real(kind=dp), save :: pfi_timedecay = 0.0_dp, pfi_distdecay = 0.0_dp
  real(kind=dp), save :: pfi_gravity(3) = 0.0_dp, pfi_damping = 0.0_dp
  real(kind=dp), save :: pfi_particlevolume = 0.0_dp, pfi_prevdtime = 0.0_dp
  integer, save :: pfi_usegradvelo = 0, pfi_gotfieldmode = 0, pfi_gotfieldweight = 0
  integer, save :: pfi_gotgravity = 0, pfi_gotdamping = 0, pfi_gottimedecay = 0
  integer, save :: pfi_gotdistdecay = 0, pfi_gotpot = 0, pfi_gotb = 0, pfi_gotpot2 = 0
  integer, save :: pfi_gotvelo = 0, pfi_velocond = 0, pfi_coordcond = 0
  integer, save :: pfi_gotbuoyancy = 0, pfi_gotfield = 0, pfi_gotmass = 0
  integer, save :: pfi_gotcharge = 0, pfi_gotdrag = 0, pfi_gotrad = 0
  type(variable_t), pointer, save :: pfi_potvar => null()
  type(variable_t), pointer, save :: pfi_potvar2 => null()
  type(variable_t), pointer, save :: pfi_bvar => null()
  type(variable_t), pointer, save :: pfi_velovar => null()
  type(variable_t), pointer, save :: pfi_velocondvar => null()
  type(variable_t), pointer, save :: pfi_coordcondvar => null()
  type(variable_t), pointer, save :: pfi_distvar => null()
  character(len=MAX_NAME_LEN), save :: pfi_fieldmode = '', pfi_fieldweight = ''
  character(len=MAX_NAME_LEN), save :: pfi_densityname = 'Density'
  type VarPointerAcc_t
    type(variable_t), pointer :: var => null()
  end type VarPointerAcc_t
  type(VarPointerAcc_t), save :: pfi_activevars(20), pfi_particlevars(20)
  integer, save :: pfi_activeopers(20) = 0, pfi_activegroups(20) = 0
  integer, save :: pfi_normvars(20) = 0
  real(kind=dp), pointer, save :: pfi_massv(:) => null()
  real(kind=dp), pointer, save :: pfi_dampingv(:) => null()
  real(kind=dp), pointer, save :: pfi_radv(:) => null()
  real(kind=dp), pointer, save :: pfi_chargev(:) => null()
  real(kind=dp), pointer, save :: pfi_dragcoeffv(:) => null()
  type(element_t), pointer, save :: pfi_bulkelement => null()
  type(element_t), pointer, save :: pfi_bulkelement2 => null()
  real(kind=dp), save :: pfi_force(3) = 0.0_dp
  real(kind=dp), save :: pfi_potatpoint = 0.0_dp
  real(kind=dp), save :: pfi_veloatpoint(3) = 0.0_dp
  real(kind=dp), save :: pfi_coord(3) = 0.0_dp, pfi_velo(3) = 0.0_dp
  integer, save :: pfi_currgroup = -1, pfi_prevgroup = -1
  character(len=MAX_NAME_LEN), save :: pfi_varname = ''

  ! -- ParticleWallContact SAVE mirror ---------------------------------
  type(mesh_t), pointer, save :: pwc_mesh => null()
  type(solver_t), pointer, save :: pwc_solver => null()
  type(valuelist_t), pointer, save :: pwc_params => null()
  integer, save :: pwc_visited = 0
  integer, save :: pwc_dim = 0, pwc_anyinteraction = 0
  real(kind=dp), save :: pwc_rad = 0.0_dp, pwc_mass = 0.0_dp
  real(kind=dp), save :: pwc_coeff = 1.0_dp, pwc_spring = 0.0_dp
  real(kind=dp), save :: pwc_wallvelo(3) = 0.0_dp
  real(kind=dp), save :: pwc_velo(3) = 0.0_dp, pwc_coord(3) = 0.0_dp
  integer, save :: pwc_truecollision = 0, pwc_gotvelovar = 0, pwc_walltrace = 0
  type(variable_t), pointer, save :: pwc_velovar => null()
  type(variable_t), pointer, save :: pwc_wallvar => null()
  real(kind=dp), allocatable, save :: pwc_basis(:)
  real(kind=dp), allocatable, save :: pwc_dbasisdx(:,:)
  type(nodes_t), target, save :: pwc_boundarynodes
  type(element_t), pointer, save :: pwc_bulkelement => null()
  type(element_t), pointer, save :: pwc_bdryelement => null()
  type(valuelist_t), pointer, save :: pwc_bc => null()
  integer, save :: pwc_nnodes = 0
  real(kind=dp), save :: pwc_val = 0.0_dp
  real(kind=dp), save :: pwc_normal(3) = 0.0_dp
  real(kind=dp), save :: pwc_gradvelo(3, 3) = 0.0_dp
  integer, pointer :: pwc_wallperm(:) => null()

  ! -- ParticleParticleContact SAVE mirror ------------------------------
  type(valuelist_t), pointer, save :: ppc_params => null()
  integer, save :: ppc_visited = 0
  real(kind=dp), save :: ppc_rad = 0.0_dp, ppc_mass = 0.0_dp
  real(kind=dp), save :: ppc_spring = 0.0_dp, ppc_damping = 0.0_dp
  real(kind=dp), save :: ppc_friction = 0.0_dp

  ! -- ParticleParticleCollision SAVE mirror ----------------------------
  type(valuelist_t), pointer, save :: ppk_params => null()
  integer, save :: ppk_visited = 0
  real(kind=dp), save :: ppk_rad1 = 0.0_dp, ppk_rad2 = 0.0_dp
  real(kind=dp), save :: ppk_mass1 = 0.0_dp, ppk_mass2 = 0.0_dp
  real(kind=dp), save :: ppk_coeff = 1.0_dp
  integer, save :: ppk_similar = 1, ppk_truecollision = 0

  ! -- ParticleWallProc SAVE mirror -------------------------------------
  type(valuelist_t), pointer, save :: pwp_params => null()
  type(solver_t), pointer, save :: pwp_solver => null()
  type(mesh_t), pointer, save :: pwp_mesh => null()
  type(particle_t), pointer, save :: pwp_particles => null()
  integer, save :: pwp_visited = 0
  integer, save :: pwp_meshdim = 0, pwp_particlebunch = 0
  type(nodes_t), target, save :: pwp_elementnodes
  type(element_t), pointer, save :: pwp_faceelement => null()
  real(kind=dp), allocatable, save :: pwp_basis(:), pwp_reflectivity(:)
  integer, allocatable, save :: pwp_indexes(:)
  type(variable_t), pointer, save :: pwp_bc_var => null()
  real(kind=dp), pointer, save :: pwp_particleenergy(:) => null()
  real(kind=dp), pointer, save :: pwp_particlecollisions(:) => null()
  real(kind=dp), pointer, save :: pwp_collisionenergy(:) => null()
  integer, pointer, save :: pwp_collisionperm(:) => null()

  ! -- TimeInfo estimate mirror -----------------------------------------
  real(kind=dp), save :: acc_cput1 = 0.0_dp
  integer, save :: acc_estindexes(6) = 0

contains

  ! -- main entry context ------------------------------------------
  subroutine pda_m_begin(h_solver, dim) bind(c, name="pda_m_begin")
    type(c_ptr), value :: h_solver
    integer(c_int) :: dim
    call c_f_pointer(h_solver, acc_solver)
    acc_params => GetSolverParams()
    acc_mesh => acc_solver % mesh
    acc_dim = CoordinateSystemDimension()
    dim = acc_dim
    acc_particles => GlobalParticles
  end subroutine pda_m_begin

  subroutine pda_m_get_int(which, val, found) bind(c, name="pda_m_get_int")
    integer(c_int), value :: which
    integer(c_int) :: val, found
    logical :: f
    integer :: v
    f = .false.; v = 0
    select case (which)
    case (1)
      v = GetInteger(acc_params, 'Time Order', f)
    case (2)
      v = GetInteger(acc_params, 'Number Of Particle Groups', f)
    case (3)
      v = GetInteger(acc_params, 'Output Interval', f)
    case (4)
      v = GetInteger(acc_params, 'Max Timestep Intervals', f)
    case (5)
      v = ListGetInteger(acc_params, 'Default Group', f)
    end select
    val = v
    found = merge(1, 0, f)
  end subroutine pda_m_get_int

  subroutine pda_m_get_log(which, val, found) bind(c, name="pda_m_get_log")
    integer(c_int), value :: which
    integer(c_int) :: val, found
    logical :: f, l
    f = .false.; l = .false.
    select case (which)
    case (1)
      l = GetLogical(acc_params, 'Particle To Field', f)
    case (2)
      l = GetLogical(acc_params, 'Particle Particle Collision', f)
    case (3)
      l = GetLogical(acc_params, 'Particle Particle Contact', f)
    case (4)
      l = GetLogical(acc_params, 'Box Particle Collision', f)
    case (5)
      l = GetLogical(acc_params, 'Box Particle Contact', f)
    case (6)
      l = GetLogical(acc_params, 'Wall Particle Contact', f)
    case (7)
      l = GetLogical(acc_params, 'Particle Scattering', f)
    case (8)
      l = GetLogical(acc_params, 'Statistical Info', f)
    case (9)
      l = GetLogical(acc_params, 'Particle Info', f)
    case (10)
      l = GetLogical(acc_params, 'Timing Info', f)
    case (11)
      l = GetLogical(acc_params, 'Particle To Field Reset', f)
    case (12)
      l = GetLogical(acc_params, 'Reinitialize Particles', f)
    case (13)
      l = ListGetLogical(acc_params, 'Set Particle Group By Domain', f)
    case (14)
      l = ListGetLogical(acc_params, 'Set Particle Group By Condition', f)
    case (15)
      l = GetLogical(acc_params, 'Particle Distance', f)
    end select
    val = merge(1, 0, l)
    found = merge(1, 0, f)
  end subroutine pda_m_get_log

  subroutine pda_m_add_truecollision() bind(c, name="pda_m_add_truecollision")
    logical :: f
    f = .false.
    call ListAddLogical(acc_params, 'True Collision Mode', f)
  end subroutine pda_m_add_truecollision

  subroutine pda_m_set_prelim(timeorder) bind(c, name="pda_m_set_prelim")
    integer(c_int), value :: timeorder
    call SetParticlePreliminaries(acc_particles, acc_dim, timeorder)
  end subroutine pda_m_set_prelim

  function pda_m_ngroups_get() result(n) bind(c, name="pda_m_ngroups_get")
    integer(c_int) :: n
    n = acc_particles % NumberOfGroups
  end function pda_m_ngroups_get

  subroutine pda_m_ngroups_set(n) bind(c, name="pda_m_ngroups_set")
    integer(c_int), value :: n
    acc_particles % NumberOfGroups = n
  end subroutine pda_m_ngroups_set

  subroutine pda_m_init_particles(group, append) bind(c, name="pda_m_init_particles")
    integer(c_int), value :: group, append
    logical :: doappend
    doappend = (append /= 0)
    if (group >= 0) then
      if (doappend) then
        call InitializeParticles(acc_particles, AppendParticles=.true., Group=group)
      else
        call InitializeParticles(acc_particles, Group=group)
      end if
    else
      call InitializeParticles(acc_particles)
    end if
  end subroutine pda_m_init_particles

  subroutine pda_m_push_ns(group) bind(c, name="pda_m_push_ns")
    integer(c_int), value :: group
    call ListPushNameSpace('group'//I2S(group)//':')
  end subroutine pda_m_push_ns

  subroutine pda_m_pop_ns() bind(c, name="pda_m_pop_ns")
    call ListPopNameSpace()
  end subroutine pda_m_pop_ns

  subroutine pda_m_info_which(which, a) bind(c, name="pda_m_info_which")
    integer(c_int), value :: which, a
    character(*), parameter :: Caller = 'ParticleDynamics'
    select case (which)
    case (1)
      call Info(Caller, 'Setting particle group by domain', Level=9)
    case (2)
      call Info(Caller, 'Initializing particles in group '//I2S(a), Level=5)
    case (3)
      call Info(Caller, 'Group '//I2S(a), Level=5)
    case (4)
      call Info(Caller, 'Setting particle group by condition', Level=9)
    end select
  end subroutine pda_m_info_which

  subroutine pda_m_info_groupcount(group, cnt) bind(c, name="pda_m_info_groupcount")
    integer(c_int), value :: group, cnt
    character(*), parameter :: Caller = 'ParticleDynamics'
    call Info(Caller, 'Group '//I2S(group)//' particles: '//I2S(cnt), Level=5)
  end subroutine pda_m_info_groupcount

  function pda_m_group_count(group) result(n) bind(c, name="pda_m_group_count")
    integer(c_int), value :: group
    integer(c_int) :: n
    n = COUNT(acc_particles % Group == group)
  end function pda_m_group_count

  subroutine pda_m_group_body(i, bodyid) bind(c, name="pda_m_group_body")
    integer(c_int), value :: i
    integer(c_int) :: bodyid
    integer :: j
    j = acc_particles % ElementIndex(i)
    if (j == 0) then
      call LocateParticles(acc_particles, pda_tramp_wall)
      j = acc_particles % ElementIndex(i)
    end if
    bodyid = acc_mesh % Elements(j) % BodyId
  end subroutine pda_m_group_body

  subroutine pda_m_group_set(i, g) bind(c, name="pda_m_group_set")
    integer(c_int), value :: i, g
    acc_particles % Group(i) = g
  end subroutine pda_m_group_set

  subroutine pda_m_group_elem(i, elemidx) bind(c, name="pda_m_group_elem")
    integer(c_int), value :: i
    integer(c_int) :: elemidx
    elemidx = acc_particles % ElementIndex(i)
  end subroutine pda_m_group_elem

  function pda_m_groupcond_begin(group) result(has) bind(c, name="pda_m_groupcond_begin")
    integer(c_int), value :: group
    integer(c_int) :: has
    character(len=MAX_NAME_LEN) :: str
    logical :: found
    str = 'Group Condition '//I2S(group)
    acc_cond_ptr => ListFind(acc_params, str, found)
    if (.not. associated(acc_cond_ptr)) then
      has = 0
    else
      has = 1
    end if
  end function pda_m_groupcond_begin

  function pda_m_groupcond_eval(i, coord) result(cond) bind(c, name="pda_m_groupcond_eval")
    integer(c_int), value :: i
    real(c_double) :: coord(3)
    real(c_double) :: cond
    cond = ExecRealFunction(acc_cond_ptr % PROCEDURE, CurrentModel, i, coord)
  end function pda_m_groupcond_eval

  subroutine pda_m_getcoord(i, coord) bind(c, name="pda_m_getcoord")
    integer(c_int), value :: i
    real(c_double) :: coord(3)
    coord = GetParticleCoord(acc_particles, i)
  end subroutine pda_m_getcoord

  subroutine pda_m_pdist_create() bind(c, name="pda_m_pdist_create")
    call ParticleVariableCreate(acc_particles, 'Particle Distance')
  end subroutine pda_m_pdist_create

  subroutine pda_m_release_waiting() bind(c, name="pda_m_release_waiting")
    call ReleaseWaitingParticles(acc_particles)
  end subroutine pda_m_release_waiting

  subroutine pda_m_statistics(which) bind(c, name="pda_m_statistics")
    integer(c_int), value :: which
    call ParticleStatistics(acc_particles, which)
  end subroutine pda_m_statistics

  subroutine pda_m_timestep(isfirst, dtime) bind(c, name="pda_m_timestep")
    integer(c_int), value :: isfirst
    real(c_double) :: dtime
    logical :: first
    first = (isfirst /= 0)
    dtime = GetParticleTimeStep(acc_particles, first)
  end subroutine pda_m_timestep

  subroutine pda_m_init_time() bind(c, name="pda_m_init_time")
    call ParticleInitializeTime(acc_particles)
  end subroutine pda_m_init_time

  subroutine pda_m_locate() bind(c, name="pda_m_locate")
    call LocateParticles(acc_particles, pda_tramp_wall)
  end subroutine pda_m_locate

  subroutine pda_m_eliminate() bind(c, name="pda_m_eliminate")
    call EliminateExitingParticles(acc_particles)
  end subroutine pda_m_eliminate

  function pda_m_nparticles() result(n) bind(c, name="pda_m_nparticles")
    integer(c_int) :: n
    n = acc_particles % NumberOfParticles
  end function pda_m_nparticles

  subroutine pda_m_box_contact() bind(c, name="pda_m_box_contact")
    call ParticleBoxContact(acc_particles)
  end subroutine pda_m_box_contact

  subroutine pda_m_neigh_create() bind(c, name="pda_m_neigh_create")
    call CreateNeighbourList(acc_particles)
  end subroutine pda_m_neigh_create

  subroutine pda_m_interact_collision(dtime) bind(c, name="pda_m_interact_collision")
    real(c_double), value :: dtime
    call ParticleParticleInteraction(acc_particles, dtime, .true., pda_tramp_collision)
  end subroutine pda_m_interact_collision

  subroutine pda_m_interact_contact(dtime) bind(c, name="pda_m_interact_contact")
    real(c_double), value :: dtime
    call ParticleParticleInteraction(acc_particles, dtime, .false., pda_tramp_contact)
  end subroutine pda_m_interact_contact

  subroutine pda_m_destroy_ghosts() bind(c, name="pda_m_destroy_ghosts")
    call DestroyGhostParticles(acc_particles)
  end subroutine pda_m_destroy_ghosts

  subroutine pda_m_advance() bind(c, name="pda_m_advance")
    call ParticleAdvanceTimestep(acc_particles)
  end subroutine pda_m_advance

  subroutine pda_m_box_periodic() bind(c, name="pda_m_box_periodic")
    call ParticleBoxPeriodic(acc_particles)
  end subroutine pda_m_box_periodic

  subroutine pda_m_delete_lost() bind(c, name="pda_m_delete_lost")
    call DeleteLostParticles(acc_particles)
  end subroutine pda_m_delete_lost

  subroutine pda_m_path_integral() bind(c, name="pda_m_path_integral")
    call ParticlePathIntegral(acc_particles)
  end subroutine pda_m_path_integral

  subroutine pda_m_save(h_model, dt, transient) bind(c, name="pda_m_save")
    type(c_ptr), value :: h_model
    real(c_double), value :: dt
    integer(c_int), value :: transient
    type(model_t), pointer :: model
    real(kind=dp) :: ldt
    logical :: ltrans
    call c_f_pointer(h_model, model)
    ldt = dt
    ltrans = (transient /= 0)
    call SaveParticleData(model, acc_solver, ldt, ltrans)
  end subroutine pda_m_save

  subroutine pda_m_time_init(nstep) bind(c, name="pda_m_time_init")
    integer(c_int), value :: nstep
    integer :: i
    acc_cput1 = CPUTime()
    acc_estindexes(1) = nstep / 2
    do i = 1, 5
      acc_estindexes(i+1) = acc_estindexes(i) / 10
    end do
  end subroutine pda_m_time_init

  subroutine pda_m_time_check(i, nstep) bind(c, name="pda_m_time_check")
    integer(c_int), value :: i, nstep
    character(*), parameter :: Caller = 'ParticleDynamics'
    real(kind=dp) :: cput2, dcput
    if (any(acc_estindexes == i)) then
      cput2 = CPUTime()
      dcput = cput2 - acc_cput1
      if (dcput > 0.5_dp) then
        write (Message, '(A,F8.3)') 'Fraction computed (s)  :', (100.0_dp)*i/nstep
        call Info(Caller, Message)
        write (Message, '(A,F8.3)') 'Consumed time (s)  :', dcput
        call Info(Caller, Message)
        write (Message, '(A,F8.3)') 'Remaining time (s) :', dcput*(nstep-i)/i
        call Info(Caller, Message)
      end if
    end if
  end subroutine pda_m_time_check

  subroutine pda_m_info_alldone() bind(c, name="pda_m_info_alldone")
    character(*), parameter :: Caller = 'ParticleDynamics'
    call Info(Caller, 'All done', Level=4)
    call Info(Caller, '-----------------------------------------', Level=4)
  end subroutine pda_m_info_alldone

  subroutine pda_m_info_start() bind(c, name="pda_m_info_start")
    character(*), parameter :: Caller = 'ParticleDynamics'
    call Info(Caller, '-----------------------------------------', Level=4)
    call Info(Caller, 'Following the path of the particles', Level=4)
  end subroutine pda_m_info_start

  function pda_m_pconst(which) result(v) bind(c, name="pda_m_pconst")
    integer(c_int), value :: which
    integer(c_int) :: v
    select case (which)
    case (1)
      v = PARTICLE_LOST
    case (2)
      v = PARTICLE_INITIATED
    case (3)
      v = PARTICLE_WALLBOUNDARY
    case (4)
      v = PARTICLE_FIXEDVELO
    case (5)
      v = PARTICLE_FIXEDCOORD
    case (6)
      v = PARTICLE_GHOST
    case default
      v = 0
    end select
  end function pda_m_pconst

  function pda_m_tol_tiny() result(t) bind(c, name="pda_m_tol_tiny")
    real(c_double) :: t
    t = tiny(1.0_dp)
  end function pda_m_tol_tiny

  ! -- PseudoNorm wrappers -----------------------------------------
  subroutine pda_n_v2sum(dim, v2sum, n) bind(c, name="pda_n_v2sum")
    integer(c_int), value :: dim
    real(c_double) :: v2sum
    integer(c_int) :: n
    integer :: i
    v2sum = 0.0_dp
    n = acc_particles % NumberOfParticles
    do i = 1, n
      v2sum = v2sum + SUM(acc_particles % Velocity(i,1:dim)**2)
    end do
  end subroutine pda_n_v2sum

  subroutine pda_n_reduce(n, v2sum) bind(c, name="pda_n_reduce")
    integer(c_int) :: n
    real(c_double) :: v2sum
    n = ParallelReduction(n)
    v2sum = ParallelReduction(v2sum)
  end subroutine pda_n_reduce

  subroutine pda_n_finish(prevnrm, nrm) bind(c, name="pda_n_finish")
    real(c_double) :: prevnrm, nrm
    prevnrm = acc_solver % Variable % Norm
    acc_solver % Variable % Norm = nrm
    acc_solver % Variable % Values = nrm
  end subroutine pda_n_finish

  ! -- ParticleScattering wrappers ----------------------------------
  subroutine pda_s_first(timeorder, dim, mfp) bind(c, name="pda_s_first")
    integer(c_int) :: timeorder, dim
    real(c_double) :: mfp
    type(valuelist_t), pointer :: params
    logical :: found
    params => GetSolverParams()
    timeorder = acc_particles % TimeOrder
    dim = acc_particles % dim
    if (.not. acc_particles % DtConstant) then
      call Fatal('ParticleScattering', 'Timestep should not be variable in this routine!')
    end if
    mfp = ListGetConstReal(params, 'Particle MFP', found)
    if (.not. found) then
      call Fatal('ParticleScattering', 'Keyword > Particle MFP < not given!')
    end if
  end subroutine pda_s_first

  function pda_s_dt() result(dt) bind(c, name="pda_s_dt")
    real(c_double) :: dt
    dt = acc_particles % dTime
  end function pda_s_dt

  subroutine pda_s_try(no, dt, dim, mfp, timeorder, scattered) bind(c, name="pda_s_try")
    integer(c_int), value :: no, dim, timeorder
    real(c_double), value :: dt, mfp
    integer(c_int) :: scattered
    real(kind=dp) :: velo(3), veloabs, pscatter, ptest, angle
    scattered = 0
    velo(1:dim) = acc_particles % Velocity(no, :)
    veloabs = SQRT(SUM(velo(1:dim)**2))
    pscatter = 1.0_dp - EXP(-(dt*veloabs)/mfp)
    ptest = EvenRandom()
    if (ptest < pscatter) then
      angle = 2*PI*EvenRandom()
      acc_particles % Velocity(no,1) = COS(angle)*velo(1) - SIN(angle)*velo(2)
      acc_particles % Velocity(no,2) = COS(angle)*velo(2) + SIN(angle)*velo(1)
      if (timeorder > 1) then
        acc_particles % PrevVelocity(no, :) = acc_particles % Velocity(no, :)
      end if
      scattered = 1
    end if
  end subroutine pda_s_try

  function pda_s_status(no) result(s) bind(c, name="pda_s_status")
    integer(c_int), value :: no
    integer(c_int) :: s
    s = acc_particles % Status(no)
  end function pda_s_status

  subroutine pda_s_finish(noscattered, h_model) bind(c, name="pda_s_finish")
    integer(c_int), value :: noscattered
    type(c_ptr), value :: h_model
    type(model_t), pointer :: model
    call c_f_pointer(h_model, model)
    write (Message, '(A,I0)') 'Number of Scattered particles: ', noscattered
    call Info('ParticleScattering', Message, Level=9)
    call ListAddConstReal(model % Simulation, 'res: Scattered particles', 1.0_dp*noscattered)
  end subroutine pda_s_finish


  ! -- ParticleFieldInteraction wrappers ------------------------------
  subroutine pda_f_touch() bind(c, name="pda_f_touch")
    type(valuelist_t), pointer :: params
    params => GetSolverParams()
  end subroutine pda_f_touch

  subroutine pda_f_first_alloc(dim, timeorder) bind(c, name="pda_f_first_alloc")
    integer(c_int) :: dim, timeorder
    integer :: i, n
    pfi_mesh => GetMesh()
    pfi_dim = pfi_mesh % MeshDim
    dim = pfi_dim
    n = pfi_mesh % MaxElementNodes
    allocate (pfi_basis(n), pfi_dbasisdx(n, 3))
    do i = 1, 20
      nullify (pfi_activevars(i) % var)
      nullify (pfi_particlevars(i) % var)
    end do
    pfi_timeorder = acc_particles % TimeOrder
    timeorder = pfi_timeorder
  end subroutine pda_f_first_alloc

  subroutine pda_f_first_gravity() bind(c, name="pda_f_first_gravity")
    character(*), parameter :: Caller = 'ParticleFieldInteraction'
    type(valuelist_t), pointer :: params
    real(kind=dp), pointer :: gwork(:,:)
    logical :: found
    params => GetSolverParams()
    pfi_gotbuoyancy = merge(1, 0, GetLogical(params, 'Particle Lift', found))
    pfi_gotgravity = 0
    if (pfi_gotbuoyancy == 1 .or. ListGetLogical(params, 'Particle Gravity', found)) then
      pfi_gotgravity = 1
      gwork => ListGetConstRealArray(CurrentModel % Constants, 'Gravity', found)
      if (associated(gwork)) then
        pfi_gravity = gwork(4, 1)*gwork(1:3, 1)
      else
        call Fatal(Caller, 'Gravity and Lift requires gravity!')
      end if
    end if
  end subroutine pda_f_first_gravity

  subroutine pda_f_first_var(which) bind(c, name="pda_f_first_var")
    integer(c_int), value :: which
    character(*), parameter :: Caller = 'ParticleFieldInteraction'
    type(valuelist_t), pointer :: params
    character(len=MAX_NAME_LEN) :: variablename
    logical :: got, found
    params => GetSolverParams()
    select case (which)
    case (1)
      variablename = ListGetString(params, 'Potential Variable Name', got)
      if (.not. got) variablename = ListGetString(params, 'Potential Variable', got)
      pfi_varname = variablename
      if (got) then
        pfi_potvar => VariableGet(pfi_mesh % Variables, trim(variablename))
        if (.not. associated(pfi_potvar)) then
          call Fatal(Caller, 'Potential field variable does not exist: '//trim(variablename))
        end if
      end if
      pfi_gotpot = merge(1, 0, got)
    case (2)
      variablename = ListGetString(params, 'Secondary Potential Variable Name', got)
      if (.not. got) variablename = ListGetString(params, 'Secondary Potential Variable', got)
      pfi_varname = variablename
      if (got) then
        pfi_potvar2 => VariableGet(pfi_mesh % Variables, trim(variablename))
        if (.not. associated(pfi_potvar2)) then
          call Fatal(Caller, 'Potential field variable does not exist: '//trim(variablename))
        end if
      end if
      pfi_gotpot2 = merge(1, 0, got)
    case (3)
      variablename = ListGetString(params, 'Magnetic Field Variable Name', got)
      if (.not. got) variablename = ListGetString(params, 'Magnetic Field Variable', got)
      pfi_varname = variablename
      if (got) then
        pfi_bvar => VariableGet(pfi_mesh % Variables, trim(variablename))
        if (.not. associated(pfi_bvar)) then
          call Fatal(Caller, 'Magnetic field variable does not exist: '//trim(variablename))
        end if
      end if
      pfi_gotb = merge(1, 0, got)
    case (4)
      variablename = ListGetString(params, 'Velocity Variable Name', got)
      if (.not. got) variablename = ListGetString(params, 'Velocity Variable', got)
      pfi_varname = variablename
      if (got) then
        pfi_velovar => VariableGet(pfi_mesh % Variables, trim(variablename))
        if (.not. associated(pfi_velovar)) then
          call Fatal(Caller, 'Velocity field variable does not exist: '//trim(variablename))
        end if
        pfi_usegradvelo = merge(1, 0, GetLogical(params, 'Velocity Gradient Correction', found))
      else
        pfi_usegradvelo = 0
      end if
      pfi_gotvelo = merge(1, 0, got)
    case (5)
      variablename = ListGetString(params, 'Velocity Condition Variable Name', got)
      if (.not. got) variablename = ListGetString(params, 'Velocity Condition Variable', got)
      pfi_varname = variablename
      if (got) then
        pfi_velocondvar => VariableGet(pfi_mesh % Variables, trim(variablename))
        if (.not. associated(pfi_velocondvar)) then
          call Fatal(Caller, 'Velocity condition field variable does not exist: '//trim(variablename))
        end if
      end if
      pfi_velocond = merge(1, 0, got)
    case (6)
      variablename = ListGetString(params, 'Coordinate Condition Variable Name', got)
      if (.not. got) variablename = ListGetString(params, 'Coordinate Condition Variable', got)
      pfi_varname = variablename
      if (got) then
        pfi_coordcondvar => VariableGet(pfi_mesh % Variables, trim(variablename))
        if (.not. associated(pfi_coordcondvar)) then
          call Fatal(Caller, 'Coordinate condition field variable does not exist: '//trim(variablename))
        end if
      end if
      pfi_coordcond = merge(1, 0, got)
    end select
  end subroutine pda_f_first_var

  subroutine pda_f_first_fields() bind(c, name="pda_f_first_fields")
    character(*), parameter :: Caller = 'ParticleFieldInteraction'
    type(valuelist_t), pointer :: params
    type(variable_t), pointer :: var
    character(len=MAX_NAME_LEN) :: str, groupname, variablename
    logical :: gotit, found, gotfw, gotdd, gottd, gotfm
    integer :: i, j, k, nfield
    params => GetSolverParams()
    pfi_fieldweight = GetString(params, 'Particle To Field Weight', gotfw)
    pfi_gotfieldweight = merge(1, 0, gotfw)
    pfi_distdecay = GetCReal(params, 'Particle Decay Distance', gotdd)
    pfi_gotdistdecay = merge(1, 0, gotdd)
    pfi_timedecay = GetCReal(params, 'Particle Decay Time', gottd)
    pfi_gottimedecay = merge(1, 0, gottd)
    pfi_fieldmode = GetString(params, 'Particle To Field Mode', gotfm)
    pfi_gotfieldmode = merge(1, 0, gotfm)
    if (gotfm) then
      call Fatal(Caller, '> Particle to Field Mode < is obsolete. Use > Field 1 < instead!')
    end if
    pfi_activeopers = 0
    pfi_activegroups = 0
    pfi_normvars = 0
    i = 0
    do i = 1, 20
      write (str, '(A,I0)') 'Field ', i
      pfi_fieldmode = ListGetString(params, str, gotit)
      if (.not. gotit) then
        pfi_maxfield = i - 1
        exit
      end if
      select case (pfi_fieldmode)
      case ('weight')
        j = 1
      case ('energy')
        j = 2
      case ('kinetic energy')
        j = 3
      case ('potential energy')
        j = 4
      case ('electrostatic energy')
        j = 5
      case ('charge')
        j = 6
      case ('speed')
        j = 7
      case ('force')
        j = 8
      case default
        if (pfi_fieldmode(1:8) == 'particle') then
          var => ParticleVariableGet(acc_particles, pfi_fieldmode)
          if (associated(var)) then
            call Info(Caller, 'Dealing with particle variable: '//trim(pfi_fieldmode))
            j = 9
            pfi_particlevars(i) % var => var
          else
            call Warn(Caller, 'Field modes does not exist as operator or variable: '//trim(pfi_varname))
          end if
        else
          call Fatal(Caller, 'Unknown field mode: '//trim(pfi_fieldmode))
        end if
      end select
      pfi_activeopers(i) = j
      k = 0
      if (acc_particles % NumberOfGroups > 0) then
        write (str, '(A,I0)') 'Group ', i
        k = ListGetInteger(params, str, gotit)
        pfi_activegroups(i) = k
      end if
      if (k > 0) then
        write (str, '(A,I0)') 'Group Name ', i
        groupname = ListGetString(params, str, gotit)
        if (.not. gotit) call Fatal(Caller, 'Unfound keyword: '//trim(str))
        variablename = 'Particle '//trim(pfi_fieldmode)//' '//trim(groupname)
      else
        var => pfi_particlevars(i) % var
        if (associated(var)) then
          variablename = trim(var % Name)
        else
          variablename = 'Particle '//trim(pfi_fieldmode)
        end if
      end if
      var => VariableGet(pfi_mesh % Variables, variablename)
      if (.not. associated(var)) then
        call Info(Caller, 'Creating variable: '//variablename)
        call VariableAddVector(pfi_mesh % Variables, pfi_mesh, acc_solver, variablename)
        var => VariableGet(pfi_mesh % Variables, variablename)
      end if
      pfi_activevars(i) % var => var
      write (str, '(A,I0)') 'Field Normalize ', i
      pfi_normvars(i) = merge(1, 0, ListGetLogical(params, str, gotit))
    end do
    call Info(Caller, 'Number of particle fields: '//I2S(pfi_maxfield), Level=6)
  end subroutine pda_f_first_fields

  subroutine pda_f_first_tail(particletofield) bind(c, name="pda_f_first_tail")
    integer(c_int), value :: particletofield
    if (pfi_gotdistdecay == 1) then
      pfi_distvar => ParticleVariableGet(acc_particles, 'Particle Distance')
      if (.not. associated(pfi_distvar)) then
        call Fatal('ParticleFieldInteraction', '> Particle Distance < should exist as variable!')
      end if
    end if
    pfi_gotfield = 0
    if (particletofield == 1 .or. pfi_gotvelo == 1 .or. pfi_gotpot == 1 .or. &
        pfi_velocond == 1 .or. pfi_coordcond == 1 .or. pfi_gotbuoyancy == 1 .or. &
        pfi_gotb == 1) pfi_gotfield = 1
    pfi_visited = 1
  end subroutine pda_f_first_tail

  subroutine pda_f_reset(dtime, fieldreset, usercoeff) bind(c, name="pda_f_reset")
    real(c_double), value :: dtime
    integer(c_int), value :: fieldreset
    real(c_double) :: usercoeff
    type(valuelist_t), pointer :: params
    logical :: found
    real(kind=dp) :: tchar, cchar
    integer :: i
    params => GetSolverParams()
    usercoeff = GetCReal(params, 'Particle To Field Coefficient', found)
    if (.not. found) usercoeff = 1.0_dp
    tchar = GetCReal(params, 'Field Decay Time', found)
    if (found) then
      cchar = EXP(-dtime/tchar)
      do i = 1, pfi_maxfield
        pfi_activevars(i) % var % Values = cchar*pfi_activevars(i) % var % Values
      end do
    end if
    if (fieldreset == 1) then
      do i = 1, pfi_maxfield
        pfi_activevars(i) % var % Values = 0.0_dp
      end do
    end if
  end subroutine pda_f_reset

  subroutine pda_f_counts(nop, nogrp) bind(c, name="pda_f_counts")
    integer(c_int) :: nop, nogrp
    nop = acc_particles % NumberOfParticles
    nogrp = acc_particles % NumberOfGroups
    pfi_nogroups = nogrp
  end subroutine pda_f_counts

  subroutine pda_f_debug() bind(c, name="pda_f_debug")
    type(variable_t), pointer :: var
    integer :: i
    if (InfoActive(20)) then
      do i = 1, 10
        var => pfi_particlevars(i) % var
        if (associated(var)) then
          print *, 'Range:', trim(var % Name), MINVAL(var % Values), MAXVAL(var % Values)
        end if
      end do
    end if
  end subroutine pda_f_debug

  subroutine pda_f_scalars(mass, damping, rad, charge, drag, gm, gd, gr, gc, gdrag) &
      bind(c, name="pda_f_scalars")
    real(c_double) :: mass, damping, rad, charge, drag
    integer(c_int) :: gm, gd, gr, gc, gdrag
    type(valuelist_t), pointer :: params
    logical :: fm, fd, fr, fc, fg
    params => GetSolverParams()
    mass = GetCReal(params, 'Particle Mass', fm)
    damping = GetCReal(params, 'Particle Damping', fd)
    rad = GetCReal(params, 'Particle Radius', fr)
    charge = GetCReal(params, 'Particle Charge', fc)
    drag = GetCReal(params, 'Particle Drag Coefficient', fg)
    gm = merge(1, 0, fm); gd = merge(1, 0, fd); gr = merge(1, 0, fr)
    gc = merge(1, 0, fc); gdrag = merge(1, 0, fg)
    pfi_gotmass = gm; pfi_gotdamping = gd; pfi_gotrad = gr
    pfi_gotcharge = gc; pfi_gotdrag = gdrag
  end subroutine pda_f_scalars

  subroutine pda_f_groups_fetch() bind(c, name="pda_f_groups_fetch")
    type(valuelist_t), pointer :: params
    logical :: fm, fd, fr, fc, fg
    params => GetSolverParams()
    pfi_massv => ListGetConstRealArray1(params, 'Particle Mass', fm)
    pfi_dampingv => ListGetConstRealArray1(params, 'Particle Damping', fd)
    pfi_radv => ListGetConstRealArray1(params, 'Particle Radius', fr)
    pfi_chargev => ListGetConstRealArray1(params, 'Particle Charge', fc)
    pfi_dragcoeffv => ListGetConstRealArray1(params, 'Particle Drag Coefficient', fg)
    pfi_gotmass = merge(1, 0, fm); pfi_gotdamping = merge(1, 0, fd)
    pfi_gotrad = merge(1, 0, fr); pfi_gotcharge = merge(1, 0, fc)
    pfi_gotdrag = merge(1, 0, fg)
  end subroutine pda_f_groups_fetch

  subroutine pda_f_checks() bind(c, name="pda_f_checks")
    character(*), parameter :: Caller = 'ParticleFieldInteraction'
    if (pfi_gotbuoyancy == 1 .and. pfi_gotrad == 0) then
      call Fatal(Caller, '> Particle Radius < is needed for buoyancy!')
    end if
    if (pfi_gotgravity == 1 .and. pfi_gotmass == 0) then
      call Warn(Caller, '> Particle Mass < is needed by gravity!')
    end if
    if ((pfi_gotpot == 1 .or. pfi_gotb == 1) .and. pfi_gotcharge == 0) then
      call Fatal(Caller, '> Particle Charge < is needed by electric and magnetic fields!')
    end if
    if (pfi_gotvelo == 1 .and. pfi_gotdrag == 0) then
      call Fatal(Caller, '> Particle Drag Coefficient < required with velocity!')
    end if
  end subroutine pda_f_checks

  subroutine pda_f_buoyvol(rad) bind(c, name="pda_f_buoyvol")
    real(c_double), value :: rad
    if (pfi_dim == 2) then
      pfi_particlevolume = PI*rad**2
    else
      pfi_particlevolume = (4.0_dp/3)*PI*rad**3
    end if
  end subroutine pda_f_buoyvol

  subroutine pda_f_state(s) bind(c, name="pda_f_state")
    integer(c_int) :: s(18)
    s(1) = pfi_dim; s(2) = pfi_timeorder; s(3) = pfi_maxfield
    s(4) = pfi_nogroups; s(5) = pfi_gotfieldweight; s(6) = pfi_gottimedecay
    s(7) = pfi_gotdistdecay; s(8) = pfi_usegradvelo; s(9) = pfi_gotpot
    s(10) = pfi_gotb; s(11) = pfi_gotpot2; s(12) = pfi_gotvelo
    s(13) = pfi_velocond; s(14) = pfi_coordcond; s(15) = pfi_gotgravity
    s(16) = pfi_gotbuoyancy; s(17) = pfi_gotfield; s(18) = pfi_visited
  end subroutine pda_f_state

  subroutine pda_f_begin_no(no, status, elemidx, group, proceed) bind(c, name="pda_f_begin_no")
    integer(c_int), value :: no
    integer(c_int) :: status, elemidx, group, proceed
    integer :: idx
    proceed = 0
    status = GetParticleStatus(acc_particles, no)
    if (status >= PARTICLE_LOST) return
    if (status <= PARTICLE_INITIATED) return
    idx = GetParticleElement(acc_particles, no)
    elemidx = idx
    if (idx < 1) return
    pfi_bulkelement => pfi_mesh % Elements(idx)
    if (.not. associated(pfi_bulkelement)) return
    pfi_coord = GetParticleCoord(acc_particles, no)
    pfi_velo = GetParticleVelo(acc_particles, no)
    pfi_force = 0.0_dp
    group = GetParticleGroup(acc_particles, no)
    proceed = 1
  end subroutine pda_f_begin_no

  subroutine pda_f_get_coordvelo(coord, velo) bind(c, name="pda_f_get_coordvelo")
    real(c_double) :: coord(3), velo(3)
    coord = pfi_coord
    velo = pfi_velo
  end subroutine pda_f_get_coordvelo

  subroutine pda_f_group_update(group, mass, damping, rad, charge, drag) &
      bind(c, name="pda_f_group_update")
    integer(c_int), value :: group
    real(c_double) :: mass, damping, rad, charge, drag
    pfi_currgroup = group
    if (group /= pfi_prevgroup) then
      if (pfi_gotmass == 1) mass = pfi_massv(MIN(SIZE(pfi_massv), group))
      if (pfi_gotdamping == 1) damping = pfi_dampingv(MIN(SIZE(pfi_dampingv), group))
      if (pfi_gotrad == 1) rad = pfi_radv(MIN(SIZE(pfi_radv), group))
      if (pfi_gotcharge == 1) charge = pfi_chargev(MIN(SIZE(pfi_chargev), group))
      if (pfi_gotdrag == 1) drag = pfi_dragcoeffv(MIN(SIZE(pfi_dragcoeffv), group))
      if (pfi_gotbuoyancy == 1) then
        if (pfi_dim == 2) then
          pfi_particlevolume = PI*rad**2
        else
          pfi_particlevolume = (4.0_dp/3)*PI*rad**3
        end if
      end if
      pfi_prevgroup = group
    end if
  end subroutine pda_f_group_update

  subroutine pda_f_force(dtime, mass, damping, charge, drag, rad, status, proceed) &
      bind(c, name="pda_f_force")
    character(*), parameter :: Caller = 'ParticleFieldInteraction'
    real(c_double), value :: dtime, mass, damping, charge, drag, rad
    integer(c_int) :: status, proceed
    real(kind=dp) :: sqrtele, weight, dist, fluiddensity, volfrac, val
    real(kind=dp) :: potatpoint, gradpotatpoint(3), veloatpoint(3)
    real(kind=dp) :: gradveloatpoint(3, 3), batpoint(3)
    type(element_t), pointer :: be2
    logical :: stat
    integer :: i
    proceed = 1
    if (pfi_gotgravity == 1) then
      pfi_force = pfi_force + pfi_gravity*mass
    end if
    if (pfi_gotdamping == 1 .and. pfi_timeorder == 2) then
      pfi_force = pfi_force - damping*pfi_velo*mass
    end if
    if (pfi_gotfield == 1) then
      if (pfi_gotpot == 1 .or. pfi_usegradvelo == 1) then
        stat = ParticleElementInfo(pfi_bulkelement, pfi_coord, &
            sqrtele, pfi_basis, pfi_dbasisdx)
      else
        stat = ParticleElementInfo(pfi_bulkelement, pfi_coord, &
            sqrtele, pfi_basis)
      end if
      if (.not. stat) then
        call Warn(Caller, 'Particle not in element')
        proceed = 0
        return
      end if
      if (pfi_velocond == 1) then
        call GetScalarFieldInMesh(pfi_velocondvar, pfi_bulkelement, pfi_basis, val)
        if (val > tiny(val)) status = PARTICLE_FIXEDVELO
      end if
      if (pfi_coordcond == 1) then
        call GetScalarFieldInMesh(pfi_coordcondvar, pfi_bulkelement, pfi_basis, val)
        if (val > tiny(val)) status = PARTICLE_FIXEDCOORD
      end if
      if (pfi_gotvelo == 1) then
        if (pfi_usegradvelo == 1) then
          call GetVectorFieldInMesh(pfi_velovar, pfi_bulkelement, pfi_basis, veloatpoint, &
              pfi_dbasisdx, gradveloatpoint)
          do i = 1, pfi_dim
            veloatpoint(i) = veloatpoint(i) + &
                0.5_dp*SUM(gradveloatpoint(i,1:pfi_dim)*pfi_velo(1:pfi_dim))*dtime
          end do
        else
          call GetVectorFieldInMesh(pfi_velovar, pfi_bulkelement, pfi_basis, veloatpoint)
        end if
        if (pfi_timeorder == 2) then
          pfi_force = pfi_force + drag*(veloatpoint - pfi_velo)
        else
          pfi_force = pfi_force + drag*veloatpoint
        end if
      else
        veloatpoint = 0.0_dp
      end if
      pfi_veloatpoint = veloatpoint
      if (pfi_gotpot == 1) then
        call GetScalarFieldInMesh(pfi_potvar, pfi_bulkelement, pfi_basis, potatpoint, &
            pfi_dbasisdx, gradpotatpoint)
        pfi_force = pfi_force - charge*gradpotatpoint
        pfi_potatpoint = potatpoint
      end if
      if (pfi_gotb == 1) then
        call GetVectorFieldInMesh(pfi_bvar, pfi_bulkelement, pfi_basis, batpoint)
        pfi_force = pfi_force + charge*CrossProduct(pfi_velo, batpoint)
        if (pfi_dim < 3) pfi_force(3) = 0.0_dp
      end if
      if (pfi_gotpot2 == 1) then
        call GetScalarFieldInMesh(pfi_potvar2, pfi_bulkelement, pfi_basis, potatpoint, &
            pfi_dbasisdx, gradpotatpoint)
        pfi_force = pfi_force - charge*gradpotatpoint
        pfi_potatpoint = potatpoint
      end if
      if (pfi_gotbuoyancy == 1) then
        if (GetParticleElementIntersection(acc_particles, pfi_bulkelement, pfi_basis, &
            pfi_coord, rad, be2, volfrac)) then
          fluiddensity = GetMaterialPropertyInMesh(pfi_densityname, pfi_bulkelement, &
              pfi_basis, be2, volfrac)
        else
          fluiddensity = GetMaterialPropertyInMesh(pfi_densityname, pfi_bulkelement, &
              pfi_basis)
        end if
        pfi_force = pfi_force - pfi_gravity*pfi_particlevolume*fluiddensity
      end if
    else
      pfi_veloatpoint = 0.0_dp
    end if
  end subroutine pda_f_force

  subroutine pda_f_accum(no, dtime, tottime, usercoeff, mass, charge, group) &
      bind(c, name="pda_f_accum")
    integer(c_int), value :: no, group
    real(c_double), value :: dtime, tottime, usercoeff, mass, charge
    character(*), parameter :: Caller = 'ParticleFieldInteraction'
    real(kind=dp) :: val, valcoeff, weight, dist, sumf, sumw
    real(kind=dp), pointer :: forcevector(:)
    real(kind=dp), pointer :: weightvar_vals(:)
    integer, pointer :: forceperm(:)
    type(variable_t), pointer :: var, weightvar
    integer, pointer :: nodeindexes(:)
    type(variable_t), pointer :: wvar
    character(len=MAX_NAME_LEN) :: wname
    logical :: gotit
    integer :: i, j, k, l, n
    if (pfi_gotfieldweight == 1) then
      select case (pfi_fieldweight)
      case ('distance')
        valcoeff = SQRT(SUM(pfi_velo**2))*pfi_prevdtime
      case ('time')
        valcoeff = pfi_prevdtime
      case ('speed')
        valcoeff = SQRT(SUM(pfi_velo**2))
      case default
        call Fatal(Caller, 'Unknown field weight: '//trim(pfi_fieldmode))
      end select
    else
      valcoeff = 1.0_dp
    end if
    if (pfi_gottimedecay == 1) then
      valcoeff = valcoeff*EXP(-tottime/pfi_timedecay)
    end if
    if (pfi_gotdistdecay == 1) then
      dist = pfi_distvar % Values(no)
      valcoeff = valcoeff*EXP(-dist/pfi_distdecay)
    end if
    valcoeff = usercoeff*valcoeff
    do k = 1, pfi_maxfield
      l = pfi_activegroups(k)
      if (l > 0) then
        if (group /= l) cycle
      end if
      select case (pfi_activeopers(k))
      case (1)
        val = 1.0_dp
      case (2)
        val = 0.5_dp*mass*SUM(pfi_velo(1:pfi_dim)**2) + &
            mass*SUM(pfi_gravity(1:pfi_dim)*pfi_coord(1:pfi_dim)) + &
            charge*pfi_potatpoint
      case (3)
        val = 0.5_dp*mass*SUM(pfi_velo(1:pfi_dim)**2)
      case (4)
        val = mass*SUM(pfi_gravity(1:pfi_dim)*pfi_coord(1:pfi_dim))
      case (5)
        val = charge*pfi_potatpoint
      case (6)
        val = charge
      case (7)
        val = SQRT(SUM(pfi_veloatpoint(1:pfi_dim)**2))
      case (8)
        val = SQRT(SUM(pfi_force(1:pfi_dim)**2))
      case (9)
        val = pfi_particlevars(k) % var % Values(no)
      end select
      val = val*valcoeff
      var => pfi_activevars(k) % var
      forcevector => var % Values
      if (associated(var % Perm)) then
        forceperm => var % Perm
      else
        nullify (forceperm)
      end if
      nodeindexes => pfi_bulkelement % NodeIndexes
      n = pfi_bulkelement % Type % NumberOfNodes
      do i = 1, n
        weight = pfi_basis(i)
        j = nodeindexes(i)
        if (associated(forceperm)) then
          j = forceperm(j)
          if (j == 0) cycle
        end if
        forcevector(j) = forcevector(j) + weight*val
      end do
    end do
  end subroutine pda_f_accum

  subroutine pda_f_get_force(force) bind(c, name="pda_f_get_force")
    real(c_double) :: force(3)
    force = pfi_force
  end subroutine pda_f_get_force

  subroutine pda_f_setparticles(no, status) bind(c, name="pda_f_setparticles")
    integer(c_int), value :: no, status
    call AddParticleForce(acc_particles, no, pfi_force)
    call SetParticleStatus(acc_particles, no, status)
  end subroutine pda_f_setparticles

  subroutine pda_f_normalize() bind(c, name="pda_f_normalize")
    character(*), parameter :: Caller = 'ParticleFieldInteraction'
    type(variable_t), pointer :: weightvar
    real(kind=dp), pointer :: forcevector(:)
    character(len=MAX_NAME_LEN) :: wname
    real(kind=dp) :: sumf, sumw
    logical :: gotit
    integer :: j
    wname = ListGetString(acc_params, 'Particle Nodal Weights', gotit)
    if (gotit) then
      call Info(Caller, 'Setting the average to zero', Level=9)
      weightvar => VariableGet(pfi_mesh % Variables, trim(wname))
      sumw = SUM(weightvar % Values)
      do j = 1, 20
        if (pfi_normvars(j) == 1) then
          forcevector => pfi_activevars(j) % var % Values
          sumf = SUM(forcevector)
          if (SIZE(weightvar % Values) /= SIZE(forcevector)) then
            call Fatal(Caller, 'Sizes are assumed to be the same')
          end if
          if (ABS(sumf) > tiny(sumf)) then
            sumw = SUM(weightvar % Values)
            forcevector = forcevector - (sumf/sumw)*weightvar % Values
          end if
        end if
      end do
    end if
  end subroutine pda_f_normalize

  subroutine pda_f_set_prevdtime(dtime) bind(c, name="pda_f_set_prevdtime")
    real(c_double), value :: dtime
    pfi_prevdtime = dtime
  end subroutine pda_f_set_prevdtime





  ! -- ParticleWallContact wrappers ----------------------------------
  subroutine pda_t_first(h_model, dim, anyinter) bind(c, name="pda_t_first")
    type(c_ptr), value :: h_model
    integer(c_int) :: dim, anyinter
    type(valuelist_t), pointer :: params, bc
    type(solver_t), pointer :: solver
    type(model_t), pointer :: model
    logical :: found, collision, contact, movingwall, accumulation
    logical :: accumulationlimit, walltrace
    integer :: k, n, wallnodes
    character(len=MAX_NAME_LEN) :: variablename
    call c_f_pointer(h_model, model)
    pwc_mesh => GetMesh()
    params => GetSolverParams()
    solver => GetSolver()
    pwc_params => params
    pwc_solver => solver
    pwc_dim = pwc_mesh % Meshdim
    dim = pwc_dim
    n = pwc_mesh % MaxElementNodes
    allocate (pwc_basis(n), pwc_dbasisdx(n, 3))
    pwc_rad = GetCReal(params, 'Wall Particle Radius', found)
    if (.not. found) pwc_rad = GetCReal(params, 'Particle Radius', found)
    if (.not. found) then
      call Fatal('ParticleWallContact', '> Particle Radius < needed!')
    end if
    collision = .false.
    contact = .false.
    movingwall = .false.
    accumulation = .false.
    accumulationlimit = .false.
    walltrace = .false.
    do k = 1, CurrentModel % NumberOfBCs
      bc => CurrentModel % BCs(k) % Values
      collision = collision .or. GetLogical(bc, 'Wall Particle Collision', found)
      contact = contact .or. GetLogical(bc, 'Wall Particle Contact', found)
      movingwall = movingwall .or. GetLogical(bc, 'Moving Wall', found)
      accumulation = accumulation .or. ListGetLogical(bc, 'Particle Accumulation', found)
      accumulationlimit = accumulationlimit .or. &
          ListCheckPresent(bc, 'Particle Accumulation Max Shear')
      walltrace = walltrace .or. ListGetLogical(bc, 'Particle Trace', found)
    end do
    if (contact) then
      pwc_spring = GetCReal(params, 'Wall Particle Spring', found)
      if (.not. found) then
        call Fatal('ParticleWallContact', '> Wall Particle Spring < needed!')
      end if
    end if
    if (collision) then
      pwc_mass = GetCReal(params, 'Particle Mass', found)
      if (.not. found) then
        call Fatal('ParticleWallContact', '> Particle Mass < needed!')
      end if
      pwc_coeff = GetCReal(params, 'Wall Particle Bounciness', found)
      if (.not. found) pwc_coeff = 1.0_dp
      pwc_truecollision = merge(1, 0, GetLogical(params, 'True Collision Mode', found))
    end if
    pwc_gotvelovar = 0
    if (movingwall .or. accumulationlimit) then
      if (movingwall) then
        if (collision) then
          variablename = ListGetString(params, 'Wall Velocity Variable Name', found)
          if (.not. found) variablename = ListGetString(params, 'Wall Velocity Variable', found)
          if (.not. found) then
            call Fatal('ParticleWallContact', 'Moving wall needs > Wall Velocity Variable <')
          end if
        else
          call Fatal('ParticleWallContact', 'Moving Wall assumes > Wall Particle Collision <')
        end if
      end if
      if (accumulationlimit) then
        if (contact) then
          variablename = ListGetString(params, 'Velocity Variable Name', found)
          if (.not. found) variablename = ListGetString(params, 'Velocity Variable', found)
          if (.not. found) then
            call Fatal('ParticleWallContact', 'Particle Accumulation needs > Velocity Variable <')
          end if
        else
          call Warn('ParticleWallContact', 'Particle Accumulation assumes > Wall Particle Collision <')
        end if
      end if
      pwc_gotvelovar = 1
      pwc_velovar => VariableGet(pwc_mesh % Variables, trim(variablename))
      if (.not. associated(pwc_velovar)) then
        call Fatal('ParticleWallContact', 'Velocity field variable does not exist: '//trim(variablename))
      end if
    end if
    if (walltrace) then
      variablename = 'Particle Trace'
      pwc_wallvar => VariableGet(pwc_mesh % Variables, variablename)
      if (.not. associated(pwc_wallvar)) then
        allocate (pwc_wallperm(pwc_mesh % NumberOfNodes))
        pwc_wallperm = 0
        call MakePermUsingMask(model, pwc_solver, pwc_mesh, variablename, .false., &
            pwc_wallperm, wallnodes)
        if (wallnodes > 0) then
          call VariableAddVector(pwc_mesh % Variables, pwc_mesh, pwc_solver, &
              variablename, Perm=pwc_wallperm)
        end if
      end if
      pwc_wallvar => VariableGet(pwc_mesh % Variables, variablename)
    end if
    pwc_anyinteraction = 0
    if (contact .or. collision) pwc_anyinteraction = 1
    anyinter = pwc_anyinteraction
    pwc_walltrace = merge(1, 0, walltrace)
    pwc_gradvelo = 0.0_dp
    pwc_wallvelo = 0.0_dp
    pwc_velo = 0.0_dp
    pwc_coord = 0.0_dp
    pwc_visited = 1
  end subroutine pda_t_first

  subroutine pda_t_state(s) bind(c, name="pda_t_state")
    integer(c_int) :: s(6)
    s(1) = pwc_dim; s(2) = pwc_anyinteraction; s(3) = pwc_visited
    s(4) = pwc_truecollision; s(5) = pwc_gotvelovar; s(6) = pwc_walltrace
  end subroutine pda_t_state

  subroutine pda_t_rad(mass, coeff, spring, rad) bind(c, name="pda_t_rad")
    real(c_double) :: mass, coeff, spring, rad
    mass = pwc_mass; coeff = pwc_coeff; spring = pwc_spring; rad = pwc_rad
  end subroutine pda_t_rad

  subroutine pda_t_begin_no(no, status, elemidx, proceed) bind(c, name="pda_t_begin_no")
    integer(c_int), value :: no
    integer(c_int) :: status, elemidx, proceed
    integer :: eidx
    proceed = 0
    status = acc_particles % Status(no)
    if (status >= PARTICLE_LOST) return
    if (status < PARTICLE_INITIATED) return
    eidx = acc_particles % ElementIndex(no)
    elemidx = eidx
    if (eidx == 0) return
    if (acc_particles % InternalElements(eidx)) return
    pwc_bulkelement => pwc_mesh % Elements(eidx)
    proceed = 1
  end subroutine pda_t_begin_no

  function pda_t_imax() result(m) bind(c, name="pda_t_imax")
    integer(c_int) :: m
    if (pwc_bulkelement % Type % DIMENSION == 3) then
      m = pwc_bulkelement % Type % NumberOfFaces
    else
      m = pwc_bulkelement % Type % NumberOfEdges
    end if
  end function pda_t_imax

  subroutine pda_t_face(face_i, found, collision, contact) bind(c, name="pda_t_face")
    integer(c_int), value :: face_i
    integer(c_int) :: found, collision, contact
    type(valuelist_t), pointer :: bc
    logical :: fnd, coll, cont, lfound
    integer :: j, jj
    found = 0; collision = 0; contact = 0
    if (pwc_bulkelement % Type % DIMENSION == 3) then
      j = pwc_bulkelement % FaceIndexes(face_i)
      pwc_bdryelement => pwc_mesh % Faces(j)
    else
      j = pwc_bulkelement % EdgeIndexes(face_i)
      pwc_bdryelement => pwc_mesh % Edges(j)
    end if
    fnd = .false.
    do jj = 1, CurrentModel % NumberOfBCs
      if (.not. associated(pwc_bdryelement % BoundaryInfo)) cycle
      if (pwc_bdryelement % BoundaryInfo % Constraint == CurrentModel % BCs(jj) % Tag) then
        fnd = .true.
        exit
      end if
    end do
    if (.not. fnd) return
    found = 1
    pwc_bc => CurrentModel % BCs(jj) % Values
    bc => pwc_bc
    coll = GetLogical(bc, 'Wall Particle Collision', lfound)
    cont = GetLogical(bc, 'Wall Particle Contact', lfound)
    collision = merge(1, 0, coll)
    contact = merge(1, 0, cont)
  end subroutine pda_t_face

  subroutine pda_t_geom(no, dim, coord, normal, dist) bind(c, name="pda_t_geom")
    integer(c_int), value :: no, dim
    real(c_double) :: coord(3), normal(3), dist
    coord(1:dim) = acc_particles % Coordinate(no, 1:dim)
    call GetElementNodes(pwc_boundarynodes, pwc_bdryelement)
    normal = NormalVector(pwc_bdryelement, pwc_boundarynodes)
    pwc_normal = normal
    pwc_coord(1:dim) = coord(1:dim)
    dist = PointFaceDistance(pwc_bdryelement, pwc_boundarynodes, coord, normal)
  end subroutine pda_t_geom

  subroutine pda_t_contact(no, dt, dim, dist, hit, exitface) bind(c, name="pda_t_contact")
    integer(c_int), value :: no, dim
    real(c_double), value :: dt, dist
    integer(c_int) :: hit, exitface
    real(kind=dp) :: eta, force(3), sqrtele, accumshear, shearrate, s
    real(kind=dp) :: speed
    logical :: stat, accumulation, found
    integer :: k, l
    hit = 0; exitface = 0
    eta = dist - pwc_rad
    if (eta > 0.0_dp) return
    force = eta*pwc_spring*pwc_normal
    acc_particles % Force(no, 1:dim) = acc_particles % Force(no, 1:dim) + force(1:dim)
    accumulation = GetLogical(pwc_bc, 'Particle Accumulation', found)
    if (accumulation) then
      stat = ParticleElementInfo(pwc_bulkelement, pwc_coord, &
          sqrtele, pwc_basis, pwc_dbasisdx)
      accumshear = GetCReal(pwc_bc, 'Particle Accumulation Max Speed', found)
      if (found) then
        pwc_velo(1:dim) = acc_particles % Velocity(no, 1:dim)
        speed = SQRT(SUM(pwc_velo(1:dim)**2))
        if (speed > accumshear) accumulation = .false.
      end if
      accumshear = GetCReal(pwc_bc, 'Particle Accumulation Max Shear', found)
      if (found .and. pwc_gotvelovar == 1) then
        call GetVectorFieldInMesh(pwc_velovar, pwc_bulkelement, pwc_basis, &
            pwc_wallvelo, pwc_dbasisdx, pwc_gradvelo)
        shearrate = 0.0_dp
        do k = 1, dim
          do l = 1, dim
            s = 0.5_dp*(pwc_gradvelo(k, l) + pwc_gradvelo(l, k))
            shearrate = shearrate + s*s
          end do
        end do
        shearrate = SQRT(shearrate)
        if (shearrate > accumshear) accumulation = .false.
      end if
      if (accumulation) then
        acc_particles % Status(no) = PARTICLE_LOST
        hit = 1
        pwc_walltrace = merge(1, 0, GetLogical(pwc_bc, 'Particle Trace', found))
        pwc_val = 1.0_dp
        exitface = 1
      end if
    end if
  end subroutine pda_t_contact

  subroutine pda_t_collision(no, dt, dim, dist, hit, exitface, skipface, walltrace) &
      bind(c, name="pda_t_collision")
    integer(c_int), value :: no, dim
    real(c_double), value :: dt, dist
    integer(c_int) :: hit, exitface, skipface, walltrace
    real(kind=dp) :: vn, dta, dtb, rn(3), v1na, v2na, v1nb, sqrtele
    real(kind=dp) :: coord(3), velo(3), force(3)
    logical :: stat, movingwall, found
    hit = 0; exitface = 0; skipface = 0
    walltrace = pwc_walltrace
    pwc_velo(1:dim) = acc_particles % Velocity(no, 1:dim)
    velo = pwc_velo
    vn = SUM(pwc_normal(1:dim)*pwc_velo(1:dim))
    movingwall = GetLogical(pwc_bc, 'Moving Wall', found)
    if (movingwall) then
      stat = ParticleElementInfo(pwc_bulkelement, pwc_coord, &
          sqrtele, pwc_basis, pwc_dbasisdx)
      call GetVectorFieldInMesh(pwc_velovar, pwc_bulkelement, pwc_basis, pwc_wallvelo)
      vn = vn - SUM(pwc_normal(1:dim)*pwc_wallvelo(1:dim))
    end if
    if (ABS(vn) < tiny(vn)) then
      skipface = 1
      return
    end if
    dta = (dist - pwc_rad)/vn
    if (dta >= dt) then
      skipface = 1
      return
    end if
    if (dta < 0.0_dp .and. dist - pwc_rad > 0.0_dp) then
      skipface = 1
      return
    end if
    coord = pwc_coord
    if (dta < 0.0_dp) then
      print *, 'Coord:', coord(1:dim)
      print *, 'Velo:', velo(1:dim)
      print *, 'Dist', dist, pwc_rad
      print *, 'vn', vn, dta, dt
    end if
    dtb = dt - dta
    rn = -pwc_normal
    v1na = SUM(velo(1:dim)*rn(1:dim))
    v2na = SUM(pwc_wallvelo(1:dim)*rn(1:dim))
    v1nb = pwc_coeff*(v2na - v1na) + v2na
    if (pwc_truecollision == 1) then
      coord = coord + dta*velo
      velo = velo + (v1nb - v1na)*rn
      coord = coord + dtb*velo
    else
      coord = coord + (v1na - v1nb)*rn*dta
      force = pwc_mass*(v1nb - v1na)*rn/dt
    end if
    if (dta < 0.0_dp) then
      print *, 'dtb', dtb
      print *, 'rn', rn
      print *, 'na', v1na, v2na, v1nb
      print *, 'Coord2', coord
      print *, 'Force2', force
    end if
    if (pwc_truecollision == 1) then
      acc_particles % Coordinate(no, 1:dim) = coord(1:dim)
      acc_particles % Velocity(no, 1:dim) = velo(1:dim)
      pwc_velo = velo
    else
      acc_particles % Coordinate(no, 1:dim) = coord(1:dim)
      acc_particles % Force(no, 1:dim) = acc_particles % Force(no, 1:dim) + force(1:dim)
    end if
    pwc_coord = coord
    hit = 1
    walltrace = merge(1, 0, GetLogical(pwc_params, 'Particle Trace', found))
    pwc_walltrace = walltrace
    if (walltrace == 1) then
      if (.not. movingwall) then
        stat = ParticleElementInfo(pwc_bulkelement, coord, &
            sqrtele, pwc_basis, pwc_dbasisdx)
      end if
      pwc_val = pwc_mass*(v1nb - v1na)/dt
    end if
    exitface = 1
  end subroutine pda_t_collision

  subroutine pda_t_trace(no, dim) bind(c, name="pda_t_trace")
    integer(c_int), value :: no, dim
    real(kind=dp) :: sumbasis, weight, val
    integer, pointer :: nodeindexes(:)
    integer :: i, j, n
    n = pwc_bulkelement % Type % NumberOfNodes
    nodeindexes => pwc_bulkelement % NodeIndexes
    sumbasis = 0.0_dp
    do i = 1, n
      j = pwc_wallvar % Perm(nodeindexes(i))
      if (j > 0) sumbasis = sumbasis + pwc_basis(i)
    end do
    if (sumbasis > tiny(sumbasis)) then
      val = pwc_val/sumbasis
      do i = 1, n
        j = pwc_wallvar % Perm(nodeindexes(i))
        if (j > 0) then
          weight = pwc_basis(i)
          pwc_wallvar % Values(j) = pwc_wallvar % Values(j) + weight*val
        end if
      end do
    end if
  end subroutine pda_t_trace

  ! -- ParticleParticleContact wrappers -------------------------------
  subroutine pda_c_first(rad, mass, spring, damping, friction) bind(c, name="pda_c_first")
    real(c_double) :: rad, mass, spring, damping, friction
    logical :: found
    ppc_params => GetSolverParams()
    if (GlobalParticles % NumberOfGroups > 1) then
      call Fatal('ParticleParticleContact', 'Implemented only for one particle type!')
    end if
    rad = GetCReal(ppc_params, 'Particle Radius', found)
    mass = GetCReal(ppc_params, 'Particle Mass', found)
    spring = GetCReal(ppc_params, 'Particle Spring', found)
    damping = GetCReal(ppc_params, 'Particle Damping', found)
    friction = GetCReal(ppc_params, 'Particle Friction', found)
    ppc_rad = rad; ppc_mass = mass; ppc_spring = spring
    ppc_damping = damping; ppc_friction = friction
    ppc_visited = 1
  end subroutine pda_c_first

  function pda_c_dot(a, b) result(s) bind(c, name="pda_c_dot")
    real(c_double) :: a(3), b(3)
    real(c_double) :: s
    s = SUM(a*b)
  end function pda_c_dot

  subroutine pda_c_warn_samepoint() bind(c, name="pda_c_warn_samepoint")
    call Warn('ParticleParticleContact', 'Particles are at same point!')
  end subroutine pda_c_warn_samepoint

  subroutine pda_c_state(visited, rad, mass, spring, damping, friction) &
      bind(c, name="pda_c_state")
    integer(c_int) :: visited
    real(c_double) :: rad, mass, spring, damping, friction
    visited = ppc_visited
    rad = ppc_rad; mass = ppc_mass; spring = ppc_spring
    damping = ppc_damping; friction = ppc_friction
  end subroutine pda_c_state

  ! -- ParticleParticleCollision wrappers -----------------------------
  subroutine pda_k_first(rad1, mass1, coeff, truecollision) bind(c, name="pda_k_first")
    real(c_double) :: rad1, mass1, coeff
    integer(c_int) :: truecollision
    logical :: found
    ppk_params => GetSolverParams()
    if (GlobalParticles % NumberOfGroups > 1) then
      call Fatal('ParticleParticleCollision', 'Implemented only for one particle type')
    end if
    rad1 = GetCReal(ppk_params, 'Particle Radius', found)
    if (.not. found) then
      call Fatal('ParticleParticleCollision', '> Particle Radius < needed!')
    end if
    coeff = GetCReal(ppk_params, 'Particle Bounciness', found)
    if (.not. found) coeff = 1.0_dp
    mass1 = GetCReal(ppk_params, 'Particle Mass', found)
    if (.not. found) then
      call Fatal('ParticleParticleCollision', '> Particle Mass < needed!')
    end if
    truecollision = merge(1, 0, GetLogical(ppk_params, 'True Collision Mode'))
    ppk_rad1 = rad1; ppk_mass1 = mass1; ppk_coeff = coeff
    ppk_truecollision = truecollision
    ppk_visited = 1
  end subroutine pda_k_first

  function pda_k_dot(a, b) result(s) bind(c, name="pda_k_dot")
    real(c_double) :: a(3), b(3)
    real(c_double) :: s
    s = SUM(a*b)
  end function pda_k_dot

  subroutine pda_k_state(visited, rad1, mass1, coeff, truecollision) &
      bind(c, name="pda_k_state")
    integer(c_int) :: visited, truecollision
    real(c_double) :: rad1, mass1, coeff
    visited = ppk_visited
    rad1 = ppk_rad1; mass1 = ppk_mass1; coeff = ppk_coeff
    truecollision = ppk_truecollision
  end subroutine pda_k_state





  ! -- ParticleWallProc wrappers ---------------------------------------
  subroutine pda_w_first(meshdim, particlebunch) bind(c, name="pda_w_first")
    integer(c_int) :: meshdim, particlebunch
    type(valuelist_t), pointer :: params
    type(solver_t), pointer :: solver
    type(variable_t), pointer :: var
    type(mesh_t), pointer :: mesh
    logical :: found
    real(kind=dp) :: totenergy
    integer :: noparticles, n, istat
    params => GetSolverParams()
    solver => GetSolver()
    mesh => GetMesh()
    pwp_params => params
    pwp_solver => solver
    pwp_mesh => mesh
    pwp_particles => GlobalParticles
    pwp_meshdim = mesh % MeshDim
    meshdim = pwp_meshdim
    particlebunch = merge(1, 0, GetLogical(params, 'Particle Bunch', found))
    pwp_particlebunch = particlebunch
    if (particlebunch == 1) then
      totenergy = ListGetConstReal(params, 'Initial Total Energy', found)
      if (.not. found) totenergy = 1.0_dp
      noparticles = ListGetInteger(params, 'Number of Particles', found)
      if (.not. found) noparticles = 1
      call ParticleVariableCreate(pwp_particles, 'particle energy')
      var => ParticleVariableGet(pwp_particles, 'particle energy')
      pwp_particleenergy => var % Values
      pwp_particleenergy = totenergy/noparticles
      call ParticleVariableCreate(pwp_particles, 'particle collisions')
      var => ParticleVariableGet(pwp_particles, 'particle collisions')
      pwp_particlecollisions => var % Values
      pwp_particlecollisions = 0
      var => VariableGet(mesh % Variables, 'Collision Energy')
      if (.not. associated(var)) then
        call Info('ParticleDynamics', 'Creating variable > Collision Energy <')
        call VariableAddVector(mesh % Variables, mesh, solver, 'Collision Energy')
        var => VariableGet(mesh % Variables, 'Collision Energy')
      end if
      pwp_collisionenergy => var % Values
      pwp_collisionperm => var % Perm
      n = mesh % MaxElementNodes
      allocate (pwp_reflectivity(n), pwp_basis(n), pwp_indexes(n), STAT=istat)
      pwp_reflectivity = 0.0_dp
      if (istat /= 0) then
        call Fatal('ParticleUtils', 'PhotonBunchivity, Memory allocation error')
      end if
    end if
    pwp_visited = 1
  end subroutine pda_w_first

  subroutine pda_w_state(visited, meshdim, particlebunch) bind(c, name="pda_w_state")
    integer(c_int) :: visited, meshdim, particlebunch
    visited = pwp_visited
    meshdim = pwp_meshdim
    particlebunch = pwp_particlebunch
  end subroutine pda_w_state

  subroutine pda_w_face(faceindex) bind(c, name="pda_w_face")
    integer(c_int), value :: faceindex
    if (pwp_meshdim == 2) then
      pwp_faceelement => pwp_mesh % Edges(faceindex)
    else
      pwp_faceelement => pwp_mesh % Faces(faceindex)
    end if
  end subroutine pda_w_face

  subroutine pda_w_advance(lambda, rinit, rfin) bind(c, name="pda_w_advance")
    real(c_double), value :: lambda
    real(c_double) :: rinit(3), rfin(3)
    rinit = rinit + lambda*(rfin - rinit)
  end subroutine pda_w_advance

  subroutine pda_w_reflect(rinit, rfin, vinit, vfin) bind(c, name="pda_w_reflect")
    real(c_double) :: rinit(3), rfin(3), vinit(3), vfin(3)
    real(kind=dp) :: normal(3), unitvector(3)
    call GetElementNodes(pwp_elementnodes, pwp_faceelement)
    normal = NormalVector(pwp_faceelement, pwp_elementnodes)
    rfin = rfin - 2*SUM((rfin - rinit)*normal)*normal
    unitvector = rfin - rinit
    unitvector = unitvector/SQRT(SUM(unitvector**2))
    vfin = unitvector*SQRT(SUM(vinit**2))
  end subroutine pda_w_reflect

  subroutine pda_w_bunch(no, rinit) bind(c, name="pda_w_bunch")
    integer(c_int), value :: no
    real(c_double) :: rinit(3)
    type(valuelist_t), pointer :: bc
    real(kind=dp) :: u, v, w, detj, localreflectivity, denergy
    integer, pointer :: nodeindexes(:)
    logical :: stat
    integer :: n, i
    n = pwp_mesh % MaxElementNodes
    if (.not. allocated(pwp_basis)) allocate (pwp_basis(n))
    call GlobalToLocal(u, v, w, rinit(1), rinit(2), rinit(3), &
        pwp_faceelement, pwp_elementnodes)
    stat = ElementInfo(pwp_faceelement, pwp_elementnodes, u, v, w, detj, pwp_basis)
    bc => GetBC()
    n = pwp_faceelement % Type % NumberOfNodes
    nodeindexes => pwp_faceelement % NodeIndexes
    pwp_reflectivity(1:n) = GetReal(bc, 'Particle Reflectivity', stat)
    localreflectivity = SUM(pwp_basis(1:n)*pwp_reflectivity(1:n))
    pwp_particlecollisions(no) = pwp_particlecollisions(no) + 1
    denergy = (1 - localreflectivity)*pwp_particleenergy(no)
    pwp_particleenergy(no) = pwp_particleenergy(no) - denergy
    if (associated(pwp_collisionperm)) then
      pwp_indexes(1:n) = pwp_collisionperm(nodeindexes(1:n))
    else
      pwp_indexes(1:n) = nodeindexes(1:n)
    end if
    pwp_collisionenergy(pwp_indexes(1:n)) = pwp_collisionenergy(pwp_indexes(1:n)) + denergy
  end subroutine pda_w_bunch

  ! -- procedure-argument trampolines (match the library interfaces) ---
  subroutine pda_tramp_wall(no, rinit, rfin, vinit, vfin, lambda, faceindex, pstatus)
    integer :: no
    real(kind=dp) :: rinit(3)
    real(kind=dp) :: rfin(3)
    real(kind=dp) :: vinit(3)
    real(kind=dp) :: vfin(3)
    real(kind=dp) :: lambda
    integer :: faceindex
    integer :: pstatus
    interface
      subroutine pda_rust_wallproc(no, rinit, rfin, vinit, vfin, lambda, &
          faceindex, pstatus) bind(c, name="__particledynamicsstuff_MOD_particlewallproc")
        use iso_c_binding
        integer(c_int) :: no, faceindex, pstatus
        real(c_double) :: rinit(3), rfin(3), vinit(3), vfin(3), lambda
      end subroutine pda_rust_wallproc
    end interface
    call pda_rust_wallproc(no, rinit, rfin, vinit, vfin, lambda, faceindex, pstatus)
  end subroutine pda_tramp_wall

  subroutine pda_tramp_contact(dt, coord, coord2, velo, velo2, force, force2, contact)
    real(kind=dp) :: dt, coord(3), coord2(3), velo(3), velo2(3), force(3), force2(3)
    logical :: contact
    interface
      subroutine pda_rust_contact(dt, coord, coord2, velo, velo2, force, force2, &
          contact) bind(c, name="__particledynamicsstuff_MOD_particleparticlecontact")
        use iso_c_binding
        real(c_double) :: dt, coord(3), coord2(3), velo(3), velo2(3), force(3), force2(3)
        integer(c_int) :: contact
      end subroutine pda_rust_contact
    end interface
    integer(c_int) :: c_contact
    call pda_rust_contact(dt, coord, coord2, velo, velo2, force, force2, c_contact)
    contact = (c_contact /= 0)
  end subroutine pda_tramp_contact

  subroutine pda_tramp_collision(dt, coord, coord2, velo, velo2, force, force2, collision)
    real(kind=dp) :: dt, coord(3), coord2(3), velo(3), velo2(3), force(3), force2(3)
    logical :: collision
    interface
      subroutine pda_rust_collision(dt, coord, coord2, velo, velo2, force, force2, &
          collision) bind(c, name="__particledynamicsstuff_MOD_particleparticlecollision")
        use iso_c_binding
        real(c_double) :: dt, coord(3), coord2(3), velo(3), velo2(3), force(3), force2(3)
        integer(c_int) :: collision
      end subroutine pda_rust_collision
    end interface
    integer(c_int) :: c_collision
    call pda_rust_collision(dt, coord, coord2, velo, velo2, force, force2, c_collision)
    collision = (c_collision /= 0)
  end subroutine pda_tramp_collision

  subroutine pda_m_particle_info(steps, timesteps, tottime) bind(c, name="pda_m_particle_info")
    integer(c_int), value :: steps, timesteps
    real(c_double), value :: tottime
    call ParticleInformation(acc_particles, steps, timesteps, tottime)
  end subroutine pda_m_particle_info

  ! -- Init entry ------------------------------------------------------
  subroutine pda_i_init(h_solver) bind(c, name="pda_i_init")
    type(c_ptr), value :: h_solver
    type(solver_t), pointer :: solver
    type(valuelist_t), pointer :: params
    call c_f_pointer(h_solver, solver)
    params => solver % Values
    call ListAddNewString(params, 'Variable', '-nooutput -global ParticleDynamics_var')
  end subroutine pda_i_init

end module particledynamics_acc


