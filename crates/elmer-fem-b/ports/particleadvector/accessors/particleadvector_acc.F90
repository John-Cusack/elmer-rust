! Port of fem/src/modules/ParticleAdvector.F90 (Elmer @ 9f6af2f85): statement-level
! BIND(C) wrappers for the Elmer calls in ParticleAdvector, compiled against the
! built .mod files of the target configure (T020 Poisson pattern, ADR-057 S3).
!
! Every wrapper takes only C-interoperable arguments (TYPE(C_PTR), VALUE
! scalars, explicit-shape REAL(C_DOUBLE) arrays). All Fortran-object traffic --
! TYPE(...) handles, OPTIONAL/POINTER dummies, assumed-shape descriptors, hidden
! lengths, array-pointer results -- stays inside Fortran, so the Rust side is
! pure C calls. Internal-procedure SAVE state (SetParticleVelocities,
! SetAdvectedField) lives in wrapper-side SAVE locals with identical lifetimes;
! OrigCoordinate (an automatic array of the entry, host-associated into
! SetAdvectedField) is mirrored as a SAVE allocatable with matching
! allocate-per-invocation / free-at-return lifetime.
!
! Clean-room code: statement-for-statement wrappers around Elmer library calls,
! no upstream text.
module particleadvector_acc
  use iso_c_binding
  use DefUtils
  use Interpolation
  use MeshBasics
  use ElementUtils
  use ParticleUtils
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(model_t), pointer, save :: acc_model => null()
  type(particle_t), pointer, save :: pav_p => null()
  real(kind=dp), allocatable, save :: acc_orig(:, :)
  character(*), parameter :: Caller = 'ParticleAdvector'

contains

  ! -- entry state ------------------------------------------------------
  subroutine pav_begin(h_model, h_solver) bind(c, name="pav_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    pav_p => GlobalParticles
  end subroutine pav_begin

  function pav_dim() result(r) bind(c, name="pav_dim")
    integer(c_int) :: r
    r = CoordinateSystemDimension()
  end function pav_dim

  subroutine pav_info_banner() bind(c, name="pav_info_banner")
    call Info(Caller, '-----------------------------------------', Level=4)
    call Info(Caller, 'Advecting fields using particle tracking', Level=4)
  end subroutine pav_info_banner

  subroutine pav_done() bind(c, name="pav_done")
    call Info(Caller, 'All done', Level=4)
    call Info(Caller, '-----------------------------------------', Level=4)
  end subroutine pav_done

  ! -- top-level ValueList queries ----------------------------------------
  subroutine pav_reverse_query(val, got) bind(c, name="pav_reverse_query")
    integer(c_int) :: val, got
    type(ValueList_t), pointer :: params
    logical :: v, g
    params => GetSolverParams()
    v = ListGetLogical(params, 'Particle time reverse', g)
    val = merge(1, 0, v)
    got = merge(1, 0, g)
  end subroutine pav_reverse_query

  function pav_reinit_query() result(r) bind(c, name="pav_reinit_query")
    integer(c_int) :: r
    type(ValueList_t), pointer :: params
    logical :: gotit
    params => GetSolverParams()
    r = merge(1, 0, GetLogical(params, 'Reinitialize Particles', gotit))
  end function pav_reinit_query

  function pav_reverse_conflict(reverse, gotit, reinit) result(r) &
      bind(c, name="pav_reverse_conflict")
    integer(c_int), value :: reverse, gotit, reinit
    integer(c_int) :: r
    r = reinit
    if ((reverse /= 0) .and. (reinit == 0)) then
      if (gotit /= 0) then
        call Fatal(Caller, 'Reverse time and no reinitialization are in conflict!')
      else
        call Info(Caller, 'Enforcing reinitialization because of reversed time!')
        r = 1
      end if
    end if
  end function pav_reverse_conflict

  function pav_fix_tangent_query() result(r) bind(c, name="pav_fix_tangent_query")
    integer(c_int) :: r
    type(ValueList_t), pointer :: params
    logical :: gotit
    params => GetSolverParams()
    r = merge(1, 0, ListGetLogical(params, 'Fix Tangent Velocity', gotit))
  end function pav_fix_tangent_query

  ! -- first-visit initialization ------------------------------------------
  subroutine pav_first_init(torder, nstep, pinfo) bind(c, name="pav_first_init")
    integer(c_int) :: torder, nstep, pinfo
    type(ValueList_t), pointer :: params
    logical :: found
    integer :: dim
    pav_p => GlobalParticles
    params => GetSolverParams()
    torder = GetInteger(params, 'Time Order', found)
    dim = CoordinateSystemDimension()
    call SetParticlePreliminaries(pav_p, dim, torder)
    nstep = GetInteger(params, 'Max Timestep Intervals', found)
    if (.not. found) nstep = 1
    if (GetLogical(params, 'Particle Info', found)) then
      pinfo = 1
    else
      pinfo = 0
    end if
  end subroutine pav_first_init

  subroutine pav_init_particles() bind(c, name="pav_init_particles")
    pav_p => GlobalParticles
    call InitializeParticles(pav_p, IsAdvector=.TRUE.)
    call ReleaseWaitingParticles(pav_p)
    pav_p % Status = PARTICLE_LOCATED
  end subroutine pav_init_particles

  subroutine pav_orig_capture() bind(c, name="pav_orig_capture")
    integer :: n, dim, i, j
    pav_p => GlobalParticles
    dim = CoordinateSystemDimension()
    n = pav_p % NumberOfParticles
    if (allocated(acc_orig)) deallocate(acc_orig)
    allocate(acc_orig(n, dim))
    do i = 1, n
      do j = 1, dim
        acc_orig(i, j) = pav_p % Coordinate(i, j)
      end do
    end do
  end subroutine pav_orig_capture

  subroutine pav_orig_free() bind(c, name="pav_orig_free")
    if (allocated(acc_orig)) deallocate(acc_orig)
  end subroutine pav_orig_free

  subroutine pav_particle_create() bind(c, name="pav_particle_create")
    pav_p => GlobalParticles
    call ParticleVariableCreate(pav_p, 'particle time')
    call ParticleVariableCreate(pav_p, 'particle distance')
  end subroutine pav_particle_create

  subroutine pav_particle_zero() bind(c, name="pav_particle_zero")
    type(Variable_t), pointer :: ptimevar, var
    pav_p => GlobalParticles
    ptimevar => ParticleVariableGet(pav_p, 'particle time')
    if (associated(ptimevar)) then
      ptimevar % Values = 0.0_dp
    end if
    var => ParticleVariableGet(pav_p, 'particle distance')
    if (associated(var)) then
      var % Values = 0.0_dp
    end if
  end subroutine pav_particle_zero

  function pav_iorder() result(r) bind(c, name="pav_iorder")
    integer(c_int) :: r
    pav_p => GlobalParticles
    if (pav_p % RK2) then
      r = 2
    else
      r = 1
    end if
  end function pav_iorder

  subroutine pav_dtsign_neg() bind(c, name="pav_dtsign_neg")
    pav_p => GlobalParticles
    pav_p % DtSign = -1
  end subroutine pav_dtsign_neg

  ! -- timestep driver ------------------------------------------------------
  function pav_maxdt(initts) result(m) bind(c, name="pav_maxdt")
    integer(c_int), value :: initts
    real(c_double) :: m
    logical :: initl
    pav_p => GlobalParticles
    initl = (initts /= 0)
    m = GetParticleTimeStep(pav_p, initl)
  end function pav_maxdt

  function pav_is_tiny(m) result(r) bind(c, name="pav_is_tiny")
    real(c_double), value :: m
    integer(c_int) :: r
    r = merge(1, 0, abs(m) < tiny(m))
  end function pav_is_tiny

  subroutine pav_tiny_exit(im1) bind(c, name="pav_tiny_exit")
    integer(c_int), value :: im1
    write (Message, '(A,I0)') 'Number of steps used: ', im1
    call Info(Caller, Message, Level=6)
  end subroutine pav_tiny_exit

  function pav_num_particles() result(r) bind(c, name="pav_num_particles")
    integer(c_int) :: r
    pav_p => GlobalParticles
    r = pav_p % NumberOfParticles
  end function pav_num_particles

  subroutine pav_advance(istep) bind(c, name="pav_advance")
    integer(c_int), value :: istep
    pav_p => GlobalParticles
    call ParticleAdvanceTimestep(pav_p, istep)
  end subroutine pav_advance

  function pav_info_active(level) result(r) bind(c, name="pav_info_active")
    integer(c_int), value :: level
    integer(c_int) :: r
    r = merge(1, 0, InfoActive(level))
  end function pav_info_active

  subroutine pav_status_count() bind(c, name="pav_status_count")
    pav_p => GlobalParticles
    call ParticleStatusCount(pav_p)
  end subroutine pav_status_count

  subroutine pav_locate() bind(c, name="pav_locate")
    pav_p => GlobalParticles
    call LocateParticles(pav_p)
  end subroutine pav_locate

  subroutine pav_path_integral(reverse, istep) bind(c, name="pav_path_integral")
    integer(c_int), value :: reverse, istep
    pav_p => GlobalParticles
    if (.not. (reverse /= 0) .or. pav_p % DtSign == 1) then
      call ParticlePathIntegral(pav_p, istep)
    end if
  end subroutine pav_path_integral

  function pav_nomoving() result(r) bind(c, name="pav_nomoving")
    integer(c_int) :: r
    integer :: m
    pav_p => GlobalParticles
    m = pav_p % NumberOfMovingParticles
    m = ParallelReduction(m)
    r = m
  end function pav_nomoving

  subroutine pav_info_step(i, nomoving) bind(c, name="pav_info_step")
    integer(c_int), value :: i, nomoving
    write (Message, '(A,I0,A,I0,A)') 'Timestep ', i, ' with ', nomoving, ' moving particles'
    call Info(Caller, Message, Level=6)
  end subroutine pav_info_step

  subroutine pav_particle_info(psteps, tsteps, tottime) bind(c, name="pav_particle_info")
    integer(c_int), value :: psteps, tsteps
    real(c_double), value :: tottime
    pav_p => GlobalParticles
    call ParticleInformation(pav_p, psteps, tsteps, tottime)
  end subroutine pav_particle_info

  function pav_reverse_turn(reverse) result(r) bind(c, name="pav_reverse_turn")
    integer(c_int), value :: reverse
    integer(c_int) :: r
    type(Variable_t), pointer :: ptimevar
    integer :: i
    pav_p => GlobalParticles
    r = 0
    if (reverse == 0) return
    if (pav_p % DtSign == -1) then
      call Info(Caller, 'Reversing time and restarting integration', Level=7)
      pav_p % DtSign = 1
      ptimevar => ParticleVariableGet(pav_p, 'particle time')
      if (associated(ptimevar)) then
        ptimevar % Values = 0.0_dp
      end if
      do i = 1, pav_p % NumberOfParticles
        call SetParticleStatus(pav_p, i, PARTICLE_MOVING)
      end do
      r = 1
    else
      call Info(Caller, 'Time reversal finished!', Level=7)
    end if
  end function pav_reverse_turn

  !> Eliminate particles that sit at fixed boundaries.
  !-------------------------------------------------------------------------
  subroutine pav_set_fixed() bind(c, name="pav_set_fixed")
    integer :: i, j, k
    type(ValueList_t), pointer :: bc, bodyforce
    type(Element_t), pointer :: element
    real(kind=dp), allocatable :: ptime(:), pcond(:)
    real(kind=dp) :: ptimeconst
    logical :: found, gottime, somebodyforce, somebodyc
    integer :: fixedcount
    type(ValueList_t), pointer :: params
    type(Variable_t), pointer :: ptimevar
    type(Mesh_t), pointer :: mesh
    integer :: n

    pav_p => GlobalParticles
    params => GetSolverParams()
    mesh => acc_solver % Mesh

    somebodyc = ListCheckPresentAnyBC(acc_model, 'Particle Fixed Condition')
    somebodyforce = ListCheckPresentAnyBodyForce(acc_model, 'Particle Fixed Condition')

    if (.not. (somebodyc .or. somebodyforce)) return

    i = acc_model % MaxElementNodes
    allocate(ptime(i), pcond(i))

    ptimevar => ParticleVariableGet(pav_p, 'particle time')
    gottime = associated(ptimevar)

    ptimeconst = ListGetCReal(params, 'Fixed Particle Time', found)
    fixedcount = 0

    if (somebodyc) then
      do j = mesh % NumberOfBulkElements + 1, &
          mesh % NumberOfBulkElements + mesh % NumberOfBoundaryElements
        element => mesh % Elements(j)
        acc_model % CurrentElement => element
        bc => GetBC(element)
        n = GetElementNOFNodes()

        pcond(1:n) = GetReal(bc, 'Particle Fixed Condition', found)
        if (.not. found) cycle

        ptime(1:n) = GetReal(bc, 'Particle Time', found)

        do i = 1, n
          if (pcond(i) < 0.0_dp) cycle

          k = element % NodeIndexes(i)
          if (pav_p % Status(k) == PARTICLE_FIXEDCOORD) cycle

          pav_p % Status(k) = PARTICLE_FIXEDCOORD
          fixedcount = fixedcount + 1
          if (.not. gottime) cycle

          if (found) then
            ptimevar % Values(k) = ptime(i)
          else
            ptimevar % Values(k) = ptimeconst
          end if
        end do
      end do
    end if

    if (somebodyforce) then
      do j = 1, mesh % NumberOfBulkElements
        element => mesh % Elements(j)
        acc_model % CurrentElement => element
        bodyforce => GetBodyForce(element)
        n = GetElementNOFNodes()

        pcond(1:n) = GetReal(bodyforce, 'Particle Fixed Condition', found)
        if (.not. found) cycle

        ptime(1:n) = GetReal(bodyforce, 'Particle Time', found)

        do i = 1, n
          if (pcond(i) < 0.0_dp) cycle

          k = element % NodeIndexes(i)
          if (pav_p % Status(k) == PARTICLE_FIXEDCOORD) cycle

          pav_p % Status(k) = PARTICLE_FIXEDCOORD
          fixedcount = fixedcount + 1
          if (.not. gottime) cycle

          if (found) then
            ptimevar % Values(k) = ptime(i)
          else
            ptimevar % Values(k) = ptimeconst
          end if
        end do
      end do
    end if
    deallocate(ptime, pcond)

    call Info('SetFixedParticles', 'Number of new fixed particles: ' &
        // I2S(fixedcount), Level=7)
  end subroutine pav_set_fixed

  !> Eliminate particles that are too old.
  !-------------------------------------------------------------------------
  subroutine pav_retired() bind(c, name="pav_retired")
    integer :: i, j, k
    real(kind=dp) :: maxintegtime
    logical :: found, gotmaxintegtime
    integer :: fixedcount
    type(ValueList_t), pointer :: params
    type(Variable_t), pointer :: ptimevar

    pav_p => GlobalParticles
    params => GetSolverParams()

    maxintegtime = ListGetCReal(params, 'Max Integration Time', gotmaxintegtime)

    if (.not. gotmaxintegtime) return

    ptimevar => ParticleVariableGet(pav_p, 'particle time')
    if (.not. associated(ptimevar)) then
      call Fatal('SetRetiredParticles', 'Cannot retire particles without time!')
    end if

    fixedcount = 0
    do k = 1, size(ptimevar % Values)
      if (pav_p % Status(k) == PARTICLE_FIXEDCOORD) cycle
      if (ptimevar % Values(k) >= maxintegtime) then
        fixedcount = fixedcount + 1
        pav_p % Status(k) = PARTICLE_FIXEDCOORD
      end if
    end do

    if (fixedcount > 0) then
      call Info('SetRetiredParticles', 'Number of new retired particles: ' &
          // I2S(fixedcount), Level=7)
    end if
  end subroutine pav_retired

  subroutine do_tangent_fix(undofix)
    logical :: undofix
    integer :: sgn, i, jv, jf, vdofs
    type(Variable_t), pointer :: velovar, fixvar
    character(len=MAX_NAME_LEN) :: variablename
    type(ValueList_t), pointer :: params
    logical :: found

    sgn = -1
    if (undofix) sgn = 1

    params => GetSolverParams()
    variablename = ListGetString(params, 'Velocity Variable Name', found)
    if (.not. found) variablename = 'flow solution'
    velovar => VariableGet(acc_solver % Mesh % Variables, trim(variablename))

    fixvar => VariableGet(acc_solver % Mesh % Variables, 'fixvelo')
    if (.not. associated(fixvar)) call Fatal('ApplyTangentFix', '"FixVelo" not found!')

    vdofs = velovar % dofs

    do i = 1, acc_solver % Mesh % NumberOfNodes
      jf = fixvar % Perm(i)
      jv = velovar % Perm(i)
      if (jf == 0 .or. jv == 0) cycle

      velovar % Values(vdofs * (jv - 1) + 3) = velovar % Values(vdofs * (jv - 1) + 3) + &
          sgn * fixvar % Values(jf)
    end do
  end subroutine do_tangent_fix

  subroutine pav_apply_tangent(undo) bind(c, name="pav_apply_tangent")
    integer(c_int), value :: undo
    call do_tangent_fix(undo /= 0)
  end subroutine pav_apply_tangent

  subroutine pav_fix_undo() bind(c, name="pav_fix_undo")
    type(ValueList_t), pointer :: params
    logical :: gotit
    params => GetSolverParams()
    if (ListGetLogical(params, 'Fix Tangent Velocity', gotit)) then
      call do_tangent_fix(.TRUE.)
    end if
  end subroutine pav_fix_undo

  !------------------------------------------------------------------------
  !> Compute field values at the given points in the FE mesh.
  !-------------------------------------------------------------------------
  subroutine pav_setvel(first) bind(c, name="pav_setvel")
    integer(c_int), value :: first
    logical :: firststep

    type(Element_t), pointer :: bulkelement
    type(Element_t), pointer :: bulkelement2
    integer :: no, status
    real(kind=dp) :: coord(3), velo(3), gradvelo(3, 3)
    type(Mesh_t), pointer :: mesh
    type(ValueList_t), pointer :: params
    real(kind=dp) :: veloatpoint(3), gradveloatpoint(3, 3), dtime
    logical :: stat, usegradvelo, visited = .FALSE.
    integer :: i, j, k, l, n, timeorder, newlost(3), oldlost, fixedlost
    integer :: elementindex, pdim
    integer, pointer :: nodeindexes(:), fieldperm(:), fieldindexes(:)
    real(kind=dp) :: sqrtelementmetric, weight, speed, speedmin
    real(kind=dp), pointer :: basis(:), dbasisdx(:, :), coordinate(:, :), velocity(:, :)
    logical :: gotit, skipzerotime
    character(len=MAX_NAME_LEN) :: variablename
    type(Variable_t), pointer :: velovar
    type(Variable_t), pointer :: dtvar

    save :: visited, mesh, basis, dbasisdx, params, velovar, usegradvelo, dtvar, &
        speedmin, newlost

    firststep = (first /= 0)
    pav_p => GlobalParticles

    if (.not. visited) then
      mesh => GetMesh()
      n = mesh % MaxElementNodes
      allocate(basis(n), dbasisdx(n, 3))

      params => GetSolverParams()

      variablename = ListGetString(params, 'Velocity Variable Name', gotit)
      if (.not. gotit) variablename = 'flow solution'
      velovar => VariableGet(mesh % Variables, trim(variablename))
      if (.not. associated(velovar)) then
        call Fatal('ParticleFieldInteraction', 'Velocity field variable does not exist: ' &
            // trim(variablename))
      end if
      usegradvelo = GetLogical(params, 'Velocity Gradient Correction', gotit)

      if (usegradvelo .and. pav_p % RK2) then
        call Warn('ParticlePathIntegral', 'Quadratic source correction incompatibe with Runge-Kutta')
        usegradvelo = .FALSE.
      end if

      if (.not. pav_p % DtConstant) then
        dtvar => ParticleVariableGet(pav_p, 'particle dt')
        if (.not. associated(dtvar)) then
          call Fatal('SetParticleVelocities', 'Required field > particle dt < not present!')
        end if
      end if

      speedmin = ListGetConstReal(params, 'Particle Min Speed', gotit)
      if (.not. gotit) speedmin = epsilon(speedmin)

      newlost = 0

      visited = .TRUE.
    end if

    coordinate => pav_p % Coordinate
    velocity => pav_p % Velocity
    coord = 0.0_dp
    velo = 0.0_dp

    oldlost = 0
    fixedlost = 0

    if (pav_p % DtConstant) then
      dtime = pav_p % DtSign * pav_p % dtime
    end if

    skipzerotime = .not. (pav_p % DtConstant .or. firststep)
    pdim = CoordinateSystemDimension()

    do no = 1, pav_p % NumberOfParticles
      status = GetParticleStatus(pav_p, no)
      if (status >= PARTICLE_LOST .or. &
          status <= PARTICLE_INITIATED .or. &
          status == PARTICLE_FIXEDCOORD .or. &
          status == PARTICLE_WALLBOUNDARY) then
        oldlost = oldlost + 1
        cycle
      end if

      elementindex = GetParticleElement(pav_p, no)
      if (elementindex == 0) then
        pav_p % Status(no) = PARTICLE_LOST
        newlost(1) = newlost(1) + 1
        cycle
      end if

      ! If the particle has not moved then it cannot have
      ! any change in the velocity.
      if (skipzerotime) then
        if (abs(dtvar % Values(no)) < tiny(dtime)) cycle
      end if

      bulkelement => mesh % Elements(elementindex)

      coord(1:pdim) = coordinate(no, 1:pdim)

      !-------------------------------------------------------------------------
      ! Get velocity from mesh
      !-------------------------------------------------------------------------
      if (usegradvelo) then
        stat = ParticleElementInfo(bulkelement, coord, &
            sqrtelementmetric, basis, dbasisdx)
      else
        stat = ParticleElementInfo(bulkelement, coord, &
            sqrtelementmetric, basis)
      end if

      if (.not. stat) then
        pav_p % Status(no) = PARTICLE_LOST
        newlost(2) = newlost(2) + 1
        cycle
      end if

      if (usegradvelo) then
        call GetVectorFieldInMesh(velovar, bulkelement, basis, veloatpoint, &
            dbasisdx, gradveloatpoint)
        if (.not. pav_p % DtConstant) then
          dtime = pav_p % DtSign * dtvar % Values(no)
        end if
        do i = 1, pdim
          velo(i) = veloatpoint(i) + &
              0.5_dp * sum(gradveloatpoint(i, 1:pdim) * veloatpoint(1:pdim)) * dtime
        end do
      else
        call GetVectorFieldInMesh(velovar, bulkelement, basis, veloatpoint)
        velo(1:pdim) = veloatpoint(1:pdim)
      end if
      speed = sqrt(sum(velo(1:pdim)**2))

      if (speed < speedmin) then
        pav_p % Status(no) = PARTICLE_FIXEDCOORD
        velocity(no, 1:pdim) = 0.0_dp
        fixedlost = fixedlost + 1
      else
        velocity(no, 1:pdim) = velo(1:pdim)
      end if
    end do

    if (.FALSE.) then
      if (newlost(1) > 0) call Warn('SetParticleVelocities', 'Some particles could not be located')
      print *, 'Total number of particles:', pav_p % NumberOfParticles
      print *, 'Passive particles:', oldlost
      print *, 'New lost particles:', newlost
      print *, 'New fixed velo particles:', fixedlost
    end if
  end subroutine pav_setvel

  !------------------------------------------------------------------------
  !> Compute field values at the given points in the FE mesh.
  !-------------------------------------------------------------------------
  subroutine pav_set_advected(dertime, dim, reversetime) bind(c, name="pav_set_advected")
    real(c_double), value :: dertime
    integer(c_int), value :: dim
    integer(c_int), value :: reversetime

    type(Element_t), pointer :: bulkelement
    integer :: no, status, noparticles
    real(kind=dp) :: dtime, coord(3), velo(3), val, vals(10)
    type(Mesh_t), pointer :: mesh
    type(ValueList_t), pointer :: params
    logical :: stat, visited = .FALSE.
    integer :: i, j, k, l, n, nsize, wallcount, novar, nonorm, dofs, maxdim, vartype
    integer :: elementindex2
    integer, pointer :: nodeindexes(:), pperm(:)
    real(kind=dp) :: sqrtelementmetric, norm, prevnorm = 0.0_dp, change
    real(kind=dp), pointer :: basis(:)
    logical :: gotit, difference, cumulative, derivative, gotvar, gotres, gotoper, &
        useperm, internalvariable, initiated, parallel, initintegral, revtime
    character(len=MAX_NAME_LEN) :: variablename, resultname, opername, name
    type(Variable_t), pointer :: targetvar, resultvar, datavar, var
    type(Variable_t), pointer :: particlevar
    real(kind=dp), pointer :: tmpvalues(:), nodevalues(:), newvalues(:)
    integer, pointer :: tmpperm(:), unitperm(:)
    real(kind=dp) :: x, y, z

    save :: visited, prevnorm, unitperm

    pav_p => GlobalParticles
    revtime = (reversetime /= 0)
    dtime = dertime

    if (revtime) then
      if (pav_p % DtSign == -1) then
        call Info(Caller, 'Setting the advected fields at turning point', Level=10)
      else
        call Info(Caller, 'Setting the advected fields after reverse tracking', Level=10)
      end if
    else
      call Info(Caller, 'Setting the advected fields', Level=10)
    end if

    mesh => GetMesh()
    maxdim = dim
    n = mesh % MaxElementNodes
    allocate(basis(n))
    coord = 0.0_dp
    velo = 0.0_dp

    params => GetSolverParams()
    nonorm = GetInteger(params, 'Norm Variable Index', gotit)
    noparticles = pav_p % NumberOfParticles

    datavar => VariableGet(mesh % Variables, 'AdvectorData')
    if (associated(datavar)) then
      nsize = size(datavar % Values)
      datavar % Output = .FALSE.
    else
      nsize = mesh % NumberOfNodes
    end if

    parallel = (ParEnv % PEs > 1)

    initiated = .FALSE.
100 novar = 0
    do while (.TRUE.)
      novar = novar + 1

      write (name, '(A,I0)') 'Variable ', novar
      variablename = GetString(params, name, gotvar)
      if (.not. gotvar) exit

      write (name, '(A,I0)') 'Operator ', novar
      opername = GetString(params, name, gotoper)

      initintegral = .FALSE.
      if (revtime) then
        if (variablename == 'particle disp' .or. variablename == 'particle disp abs') then
          ! This is done only at end of forward cycle
          if (pav_p % DtSign == -1) cycle
        else if (variablename == 'particle time integral' .or. &
            variablename == 'particle distance integral') then
          ! These are done both at forward and backward cycle
          ! End of backward cycle, initialize the integral
          ! End of forward cycle, project particle info to nodal info
          initintegral = (pav_p % DtSign == -1)
        else
          ! This is done at end of backward cycle
          if (pav_p % DtSign == 1) cycle
        end if
      end if

      call Info(Caller, 'Setting field ' // I2S(novar) // ' for variable: ' // trim(variablename), Level=15)

      ! Get the target variables
      ! Variables starting with 'particle' as associated with particles
      !----------------------------------------------------------------
      targetvar => null()
      resultvar => null()

      if (variablename == 'particle coordinate' .or. &
          variablename == 'particle velocity' .or. &
          variablename == 'particle force' .or. &
          variablename == 'particle disp') then

        dofs = dim
        internalvariable = .TRUE.
      else if (SEQL(variablename, 'particle')) then
        dofs = 1
        internalvariable = .TRUE.
      else
        targetvar => VariableGet(mesh % Variables, trim(variablename))
        if (.not. associated(targetvar)) then
          call Fatal(Caller, 'Could not obtain target variable:' &
              // trim(variablename))
          cycle
        end if
        dofs = targetvar % dofs
        if (dofs > 1) then
          call Fatal(Caller, 'Advection implemented so far only for scalars')
        end if
        internalvariable = .FALSE.
        maxdim = max(dofs, maxdim)
      end if

      if (.not. initiated) cycle

      ! For internal variables the target name is the field name
      !---------------------------------------------------------
      difference = .FALSE.
      derivative = .FALSE.
      cumulative = .FALSE.

      if (gotoper) then
        if (opername == 'difference') then
          difference = .TRUE.
        else if (opername == 'derivative') then
          derivative = .TRUE.
        else if (opername == 'cumulative' .or. opername == 'age') then
          cumulative = .TRUE.
        else
          call Fatal('SetAdvectedField', 'Unknown operator: ' // trim(opername))
        end if
      else
        opername = 'adv'
      end if

      call Info(Caller, 'Using operator ' // I2S(novar) // ' for variable: ' // trim(opername), Level=15)

      write (name, '(A,I0)') 'Result Variable ', novar
      resultname = GetString(params, name, gotres)
      if (.not. gotres) then
        resultname = trim(opername) // ' ' // trim(variablename)
      end if

      ! Create variables if they do not exist
      !---------------------------------------------------------
      resultvar => VariableGet(mesh % Variables, trim(resultname))
      if (associated(resultvar)) then
        if (associated(datavar)) then
          if (datavar % TYPE /= resultvar % TYPE) then
            call Fatal(Caller, 'ResultVar is of wrong type, use new name for result variable!')
          end if
          if (size(datavar % Values) /= size(resultvar % Values)) then
            call Fatal(Caller, 'ResultVar is of wrong size, use new name for result variable!')
          end if
        end if
        call Info(Caller, 'Found a pre-existing result variable: ' // trim(resultname), Level=20)
      else
        gotit = .FALSE.
        pperm => null()
        if (associated(datavar)) then
          call Info(Caller, 'Using non-nodal given permutation for data', Level=15)
          pperm => datavar % Perm
          vartype = datavar % TYPE
        else if (associated(targetvar)) then
          call Info(Caller, 'Using inherited permutation for data', Level=15)
          pperm => targetvar % Perm
          vartype = targetvar % TYPE
        end if

        if (.not. associated(pperm)) then
          if (.not. associated(unitperm)) then
            call Info(Caller, 'Creating unity permutation for data', Level=15)
            allocate(unitperm(nsize))
            do i = 1, nsize
              unitperm(i) = i
            end do
          end if
          pperm => unitperm
          vartype = 0
        end if

        call VariableAddVector(mesh % Variables, mesh, acc_solver, resultname, dofs, &
            Perm=pperm, VarType=vartype)

        if (dofs == 1) then
          call Info(Caller, 'Created a scalar variable: ' // trim(resultname))
        else
          call Info(Caller, 'Created a vector variable: ' // trim(resultname))
        end if
        resultvar => VariableGet(mesh % Variables, trim(resultname))
        if (.not. associated(resultvar)) call Fatal(Caller, 'Problems in VariableAdd')
      end if

      ! Finally, set the values
      !---------------------------------------------------------
      if (internalvariable .and. .not. initintegral) then
        call Info(Caller, 'Setting particle variable "' // trim(variablename) // '" to fields', Level=15)

        if (variablename == 'particle disp abs') then
          if (resultvar % Dofs /= 1) then
            call Fatal(Caller, 'Variable should have one dof: ' // trim(variablename))
          end if
          do i = 1, noparticles
            newvalues(i) = sqrt(sum((pav_p % Coordinate(i, :) - acc_orig(i, :))**2))
          end do

        else if (variablename == 'particle disp') then
          if (resultvar % Dofs /= dim) then
            call Fatal(Caller, 'Variable should have dim dofs: ' // trim(variablename))
          end if
          do i = 1, noparticles
            do j = 1, dim
              newvalues(dim * (i - 1) + j) = pav_p % Coordinate(i, j) - acc_orig(i, j)
            end do
          end do

        else if (variablename == 'particle coordinate') then
          if (resultvar % Dofs /= dim) then
            call Fatal(Caller, 'Variable should have dim dofs: ' // trim(variablename))
          end if
          do j = 1, dim
            do i = 1, noparticles
              newvalues(dim * (i - 1) + j) = pav_p % Coordinate(i, j)
            end do
          end do
        else if (variablename == 'particle coordinate_abs') then
          do i = 1, noparticles
            val = 0.0_dp
            do j = 1, dim
              val = val + pav_p % Coordinate(i, j)**2
            end do
            newvalues(i) = sqrt(val)
          end do
        else if (variablename == 'particle velocity') then
          if (resultvar % Dofs /= dim) then
            call Fatal(Caller, 'Variable should have dim dofs: ' // trim(variablename))
          end if
          do i = 1, noparticles
            do j = 1, dim
              newvalues(dim * (i - 1) + j) = pav_p % Velocity(i, j)
            end do
          end do

        else if (variablename == 'particle velocity_abs') then
          do i = 1, noparticles
            val = 0.0_dp
            do j = 1, dim
              val = val + pav_p % Velocity(i, j)**2
            end do
            newvalues(i) = sqrt(val)
          end do

        else if (variablename == 'particle force') then
          if (resultvar % Dofs /= dim) then
            call Fatal(Caller, 'Variable should have dim dofs: ' // trim(variablename))
          end if
          do i = 1, noparticles
            do j = 1, dim
              newvalues(dim * (i - 1) + j) = pav_p % Force(i, j)
            end do
          end do

        else if (variablename == 'particle status') then
          newvalues(1:noparticles) = 1.0_dp * pav_p % Status(1:noparticles)

        else if (variablename == 'particle partition') then
          if (associated(pav_p % Partition)) then
            newvalues(1:noparticles) = ParEnv % MyPe + 1.0_dp
          end if

        else if (variablename == 'particle number') then
          do i = 1, noparticles
            newvalues(i) = 1.0_dp * i
          end do

        else if (variablename == 'particle index') then
          newvalues(1:noparticles) = 1.0_dp * pav_p % NodeIndex(1:noparticles)

        else if (SEQL(variablename, 'particle')) then
          particlevar => ParticleVariableGet(pav_p, variablename)
          if (associated(particlevar)) then
            newvalues(1:noparticles) = particlevar % Values(1:noparticles)
          else
            call Warn(Caller, 'Field does not exist: ' // trim(variablename))
          end if
        end if

      else
        call Info(Caller, 'Setting field variable "' // trim(variablename) // '" to advected fields', Level=15)

        if (initintegral) then
          ! This is a cludge as the path integrals are somewhat special
          targetvar => resultvar
        end if

        do i = 1, noparticles
          status = GetParticleStatus(pav_p, i)

          if (status >= PARTICLE_LOST) cycle
          if (status <= PARTICLE_INITIATED) cycle

          elementindex2 = GetParticleElement(pav_p, i)

          if (elementindex2 == 0) cycle

          bulkelement => mesh % Elements(elementindex2)

          coord(1:dim) = pav_p % Coordinate(i, 1:dim)
          velo(1:dim) = pav_p % Velocity(i, 1:dim)

          stat = ParticleElementInfo(bulkelement, coord, sqrtelementmetric, basis)
          if (.not. stat) cycle

          if (dofs == 1) then
            call GetScalarFieldInMesh(targetvar, bulkelement, basis, val)
            newvalues(i) = val
          else
            call GetVectorFieldInMesh(targetvar, bulkelement, basis, vals)
            do j = 1, dofs
              newvalues(dofs * (i - 1) + j) = vals(j)
            end do
          end if
        end do

        ! These are initializations related to path integrals.
        ! No parallel communication required for these.
        if (initintegral) then
          if (associated(targetvar)) then
            particlevar => ParticleVariableGet(pav_p, variablename)
            if (associated(particlevar)) then
              call Info(Caller, 'Initialize "' // trim(particlevar % Name) // '" with "' &
                  // trim(targetvar % Name) // '"')
              particlevar % Values = newvalues(1:noparticles)
            end if
          end if
          targetvar => null()
          cycle
        end if

      end if

      ! In a serial case the nodes and particles are directly associated.
      ! In a parallel case we need to transfer the values from particles in
      ! different partitions to nodes.
      !---------------------------------------------------------------------
      if (parallel) then
        nodevalues = 0._dp
        call ParticleAdvectParallel(pav_p, newvalues, nodevalues, dofs)
      end if

      ! Finally move the nodal values to the target variable
      !---------------------------------------------------------------------
      if (associated(datavar)) then
        n = size(resultvar % Values)
        if (difference .or. derivative) then
          resultvar % Values = nodevalues(1:n) - targetvar % Values
        else if (cumulative) then
          resultvar % Values = nodevalues(1:n) + resultvar % Values
        else
          resultvar % Values = nodevalues(1:n)
        end if
      else
        do i = 1, nsize
          k = i
          if (associated(resultvar % Perm)) k = resultvar % Perm(k)
          if (k == 0) cycle
          dofs = resultvar % Dofs
          do j = 1, dofs
            if (difference .or. derivative) then
              resultvar % Values(dofs * (k - 1) + j) = &
                  nodevalues(dofs * (i - 1) + j) - targetvar % Values(dofs * (k - 1) + j)
            else if (cumulative) then
              resultvar % Values(dofs * (k - 1) + j) = &
                  nodevalues(dofs * (i - 1) + j) + resultvar % Values(dofs * (k - 1) + j)
            else
              resultvar % Values(dofs * (k - 1) + j) = nodevalues(dofs * (i - 1) + j)
            end if
          end do
        end do
      end if

      if (derivative) resultvar % Values = resultvar % Values / dertime

      block
        integer :: t, localperm(10)
        real(kind=dp) :: cval
        type(Element_t), pointer :: element
        real(kind=dp) :: dgscale
        logical :: gotscale

        if (resultvar % TYPE == variable_on_nodes_on_elements) then

          dgscale = ListGetCReal(acc_solver % Values, 'DG Nodes Scale', gotscale)
          if (.not. gotscale) dgscale = 1.0 / sqrt(3.0_dp)
          gotscale = (abs(dgscale - 1.0_dp) > tiny(dgscale))

          if (gotscale) then
            call Info(Caller, 'Expanding shrinked DG field', Level=12)
            do t = 1, mesh % NumberOfBulkElements
              element => mesh % Elements(t)
              n = element % TYPE % NumberOfNodes
              localperm(1:n) = resultvar % Perm(element % DGIndexes)
              if (any(localperm(1:n) == 0)) cycle
              cval = sum(resultvar % Values(localperm(1:n))) / n
              do i = 1, n
                j = localperm(i)
                resultvar % Values(j) = cval + &
                    (resultvar % Values(j) - cval) * (1.0_dp / dgscale)
              end do
            end do
          end if
        end if
      end block

      ! To allow computation of change in the standard manner the Variable
      ! is set to point to the one of interest. This is mainly used in the
      ! tests, or could be used in for convergence monitoring also.
      !---------------------------------------------------------------
      if (novar == nonorm) then
        n = size(resultvar % Values)
        norm = sum(resultvar % Values**2)

        n = ParallelReduction(n)
        norm = ParallelReduction(norm)
        norm = sqrt(norm / n)

        change = 2.0 * abs(norm - prevnorm) / (norm + prevnorm)
        prevnorm = norm

        acc_solver % Variable % Norm = norm
        acc_solver % Variable % NonlinChange = change
        acc_solver % Variable % Values = norm

        ! Here the name is ComputeChange in order to get the change also to ElmerGUI
        ! albeit in a very dirt style. One could also edit ElmerGUI....
        write (Message, '(a,g15.8,g15.8,a)') &
            'SS (ITER=1) (NRM,RELC): (', norm, change, ' ) :: ' // trim(resultname)
        call Info('ComputeChange', Message, Level=3)
      end if
    end do

    ! Allocate the local new temporal vectors
    if (.not. initiated) then
      call Info(Caller, 'Allocating for temporal value vectors', Level=15)
      novar = novar - 1
      if (novar < 1) then
        call Fatal(Caller, 'No target and result variables exist!')
      end if
      allocate(newvalues(maxdim * pav_p % NumberOfParticles))
      newvalues = 0.0_dp
      if (parallel) then
        allocate(nodevalues(maxdim * nsize))
        nodevalues = 0.0_dp
      else
        nodevalues => newvalues
      end if
      initiated = .TRUE.
      goto 100
    end if

    deallocate(newvalues)
    if (parallel) deallocate(nodevalues)
    visited = .TRUE.
  end subroutine pav_set_advected

  !------------------------------------------------------------------------------
  !> Initialization for the primary solver: ParticleAdvector
  !------------------------------------------------------------------------------
  subroutine pav_init_body() bind(c, name="pav_init_body")
    type(ValueList_t), pointer :: params
    logical :: found, advelem, advdg, advip

    params => GetSolverParams()

    ! These are default setting that make the operation of the advection solver
    ! possible. There should always be one passive particle for each active node.
    !---------------------------------------------------------------------------
    advelem = ListGetLogical(params, 'Advect Elemental', found)
    advdg = ListGetLogical(params, 'Advect DG', found)
    advip = ListGetLogical(params, 'Advect Ip', found)

    if (advelem) then
      call ListAddString(params, NextFreeKeyword('Exported Variable', params), '-elem AdvectorData')
    else if (advdg) then
      call ListAddString(params, NextFreeKeyword('Exported Variable', params), '-dg AdvectorData')
    else if (advip) then
      call ListAddString(params, NextFreeKeyword('Exported Variable', params), '-ip AdvectorData')
    end if

    call ListAddInteger(params, 'Time Order', 0)
    call ListAddNewLogical(params, 'Particle Accurate At Face', .FALSE.)

    ! If we want to show a pseudonorm add a variable for which the norm
    ! is associated with.
    call ListAddNewString(acc_solver % Values, 'Variable', '-nooutput -global particleadvector_var')
    call ListAddNewLogical(params, 'No Matrix', .TRUE.)
  end subroutine pav_init_body

end module particleadvector_acc
