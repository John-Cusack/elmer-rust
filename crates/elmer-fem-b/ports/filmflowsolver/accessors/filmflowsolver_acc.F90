! T022-A FilmFlowSolver port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/FilmFlowSolver.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 section 3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! hidden lengths, assumed-shape descriptor dummies, formatted I/O —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module mirrors FilmFlowSolver's own SAVE state (entry-level SAVE
! arrays and flags plus the CONTAINS internals' SAVE Nodes/Visited with
! identical lifetimes). All floating-point work (the bulk/boundary
! Galerkin assembly, MATMUL/SUM reductions, WHERE clamps) and all
! formatted/whole-array I/O stay in Fortran (same compiler options as
! the module build), so bit-identity holds by construction wherever the
! call sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module filmflowsolver_acc
  use iso_c_binding
  use defutils
  use parallelutils, only: ParallelSumNodalVector
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(model_t), pointer, save :: acc_model => null()
  type(valuelist_t), pointer, save :: acc_params => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(mesh_t), pointer, save :: acc_mesh => null()
  type(valuelist_t), pointer, save :: acc_bf => null()
  type(valuelist_t), pointer, save :: acc_mat => null()
  type(valuelist_t), pointer, save :: acc_bc => null()
  type(variable_t), pointer, save :: acc_thisvar => null()
  type(variable_t), pointer, save :: acc_pvar => null()
  type(variable_t), pointer, save :: acc_hvar => null()
  real(kind=dp), pointer, save :: acc_heatingenergy(:) => null()
  real(kind=dp), pointer, save :: acc_fricflux(:) => null()
  real(kind=dp), pointer, save :: acc_presflux(:) => null()

  ! Entry-level host scalars (mirrors the solver body locals).
  integer, save :: acc_mdim = 0, acc_dim = 0, acc_maxiter = 1, acc_active = 0
  integer, save :: acc_gapdir = 0, acc_fricmodel = 0, acc_coupled = 0
  integer, save :: acc_itime = 0, acc_itime0 = -1
  logical, save :: acc_gradp = .false., acc_lateral = .false.
  logical, save :: acc_gotac = .false., acc_surfac = .false.
  logical, save :: acc_useprevgap = .false., acc_convect = .false.
  logical, save :: acc_newton = .false., acc_transient = .false.
  logical, save :: acc_cs = .false., acc_gotgrav = .false.
  logical, save :: acc_gotheight = .false.
  logical, save :: acc_calcfric = .false., acc_calcpres = .false.
  logical, save :: acc_calcheat = .false., acc_useheat = .false.
  logical, save :: acc_doneweight = .false., acc_firstround = .true.
  logical, save :: acc_allocdone = .false.
  real(kind=dp), save :: acc_dt = 0.0_dp, acc_mingap = 0.0_dp
  real(kind=dp), save :: acc_gapfactor = 0.0_dp, acc_grav = 0.0_dp
  real(kind=dp), save :: acc_fsimult = 1.0_dp, acc_nm = 0.0_dp
  real(kind=dp), save :: acc_meltheat = 0.0_dp, acc_cp = 0.0_dp
  real(kind=dp), save :: acc_ct = 0.0_dp, acc_norm = 0.0_dp
  character(:), allocatable, save :: acc_densityname, acc_viscname

  ! Entry-level SAVE arrays (allocated once, first-mesh sizing).
  real(kind=dp), allocatable, save :: acc_stiff(:,:), acc_mass(:,:)
  real(kind=dp), allocatable, save :: acc_force(:), acc_load(:,:)
  real(kind=dp), allocatable, save :: acc_rho(:), acc_ac(:), acc_gap(:)
  real(kind=dp), allocatable, save :: acc_gap0(:), acc_height(:), acc_mu(:)
  real(kind=dp), allocatable, save :: acc_acpres(:), acc_velo(:,:)
  real(kind=dp), allocatable, save :: acc_acprevp(:), acc_fsirhs(:,:)
  real(kind=dp), allocatable, save :: acc_prevgap(:), acc_heatingw(:)

  ! CONTAINS-internal SAVE mirrors (identical lifetimes).
  type(nodes_t), target, save :: acc_nodes_bulk, acc_nodes_bnd
  logical, save :: acc_bvisited = .false.
  real(kind=dp), save :: acc_minpres = 0.0_dp, acc_minspeed = 0.0_dp

contains

  ! -- FilmFlowSolver_Init0 body --------------------------------------
  subroutine ffw_init0_begin(h_model, h_solver) bind(c, name="ffw_init0_begin")
    implicit none
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine ffw_init0_begin

  subroutine ffw_init0_body() bind(c, name="ffw_init0_body")
    implicit none
    type(valuelist_t), pointer :: params
    logical :: found, serendipity
    params => GetSolverParams()
    serendipity = GetLogical(GetSimulation(), 'Serendipity P Elements', found)
    if (.not. found) serendipity = .true.
    call ListAddNewInteger(params, 'Time derivative Order', 1)
    if (serendipity) then
      call ListAddNewString(params, 'Element', 'p:1 -line b:1 -tri b:1 -quad b:3')
    else
      call ListAddNewString(params, 'Element', 'p:1 -line b:1 -tri b:1 -quad b:4')
    end if
  end subroutine ffw_init0_body

  ! -- FilmFlowSolver_Init body ---------------------------------------
  subroutine ffw_init_begin(h_model, h_solver) bind(c, name="ffw_init_begin")
    implicit none
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine ffw_init_begin

  subroutine ffw_init_body() bind(c, name="ffw_init_body")
    implicit none
    type(valuelist_t), pointer :: params
    logical :: found, found2
    integer :: mdim
    params => GetSolverParams()
    mdim = ListGetInteger(params, 'Model Dimension', found)
    if (.not. found) then
      call Fatal('FilmFlowSolver_init', &
          'Give "Model Dimension" i.e. the dimension of N-S equation!')
    end if
    if (mdim == 2) then
      call ListAddNewString(params, 'Variable', &
          'Flow[FilmVelocity:2 FilmPressure:1]')
    else if (mdim == 1) then
      call ListAddNewString(params, 'Variable', &
          'Flow[FilmSpeed:1 FilmPressure:1]')
    else
      call Fatal('FilmFlowSolver_init', &
          'This module does not make sense in dim: '//I2S(mdim))
    end if
    call ListAddNewInteger(params, 'Nonlinear System Norm DOFs', mdim)
    call ListAddNewLogical(params, 'Relative Pressure Relaxation', .true.)
    call ListAddNewLogical(params, 'Bubbles in Global System', .false.)
    call ListAddNewLogical(params, 'Global Mass Matrix', .true.)
    if (ListGetLogical(params, 'Calculate Friction Heating', found)) then
      call ListAddString(params, NextFreeKeyword('Exported Variable ', params), &
          'FilmFriction HeatFlux')
    end if
    if (ListGetLogical(params, 'Calculate Pressure Heating', found2)) then
      call ListAddString(params, NextFreeKeyword('Exported Variable ', params), &
          'FilmPressure HeatFlux')
    end if
    if (found .or. found2) then
      call ListAddString(params, NextFreeKeyword('Exported Variable ', params), &
          'Heating Energy')
    end if
    if (ListGetLogical(params, 'Calculate Sensitivity', found)) then
      call ListAddString(params, NextFreeKeyword('Exported Variable ', params), &
          '-dofs '//I2S(mdim+1)//' Flow Sensitivity')
      call ListAddString(params, 'Sensitivity Variable', 'Flow Sensitivity')
    end if
  end subroutine ffw_init_body

  ! -- sweep setup (reads + first-mesh allocation) ---------------------
  ! Returns 1 when the original takes the Skip-First-Solution early
  ! RETURN path, 0 to continue.
  subroutine ffw_begin(h_model, h_solver, dt, transient, code) &
      bind(c, name="ffw_begin")
    implicit none
    type(c_ptr), value :: h_model, h_solver
    real(c_double), value :: dt
    integer(c_int), value :: transient
    integer(c_int) :: code
    type(valuelist_t), pointer :: params
    type(variable_t), pointer :: hvar
    real(kind=dp), pointer :: gwork(:,:)
    logical :: found, gotgrav
    integer :: n, istat
    character(:), allocatable :: str

    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    acc_dt = dt
    acc_transient = (transient /= 0)

    call Info('FilmFlowSolver', &
        'Computing reduced dimensional Navier-Stokes equations!')

    acc_mesh => GetMesh()
    acc_elem => GetActiveElement(1)

    acc_dim = CoordinateSystemDimension()
    params => GetSolverParams()
    acc_params => params
    acc_thisvar => acc_solver % Variable

    acc_mdim = ListGetInteger(params, 'Model Dimension', UnFoundFatal=.true.)
    acc_gradp = GetLogical(params, 'GradP Discretization', found)
    acc_lateral = GetLogical(params, 'Lateral Strain', found)
    acc_mingap = ListGetCReal(params, 'Min Gap Height', found)
    if (.not. found) acc_mingap = 1.0e-20_dp
    acc_gotac = ListCheckPresentAnyMaterial(acc_model, 'Artificial Compressibility')
    acc_useprevgap = ListGetLogical(params, 'Use Gap Average', found)
    acc_convect = GetLogical(params, 'Convect', found)
    acc_coupled = GetCoupledIter()

    acc_densityname = ListGetString(params, 'Density Name', found)
    if (.not. found) acc_densityname = 'Density'
    acc_viscname = ListGetString(params, 'Viscosity Name', found)
    if (.not. found) acc_viscname = 'Viscosity'

    acc_gapdir = 0
    acc_gapfactor = ListGetCReal(params, 'Gap Addition Factor', found)
    if (found) then
      acc_gapdir = acc_mdim + 1
      if (abs(acc_gapfactor) > 1.0_dp) then
        call Warn('FilmFlowSolver', &
            '"Gap Addition Factor" greater to unity does not make sense!')
      end if
    end if

    acc_cs = ListGetLogical(params, 'Axi Symmetric', found)
    if (.not. found) then
      acc_cs = (CurrentCoordinateSystem() == AxisSymmetric .or. &
          CurrentCoordinateSystem() == CylindricSymmetric)
    end if

    acc_fricmodel = 0
    str = ListGetString(params, 'Friction model', found)
    if (found) then
      if (str == 'laminar') then
        acc_fricmodel = 0
      else if (str == 'darcy') then
        acc_fricmodel = 1
      else if (str == 'darcy2') then
        acc_fricmodel = 2
      else if (str == 'manning') then
        acc_fricmodel = 3
      else if (str == 'manning2') then
        acc_fricmodel = 4
      else
        call Fatal('FilmFlowSolver', 'Uknown friction model: '//trim(str))
      end if
      call Info('FilmFlowSolver', 'Using friction model: '//trim(str), level=7)
    end if

    block
      logical :: usegravity
      usegravity = ListGetLogical(params, 'Use Gravity', found)
      gwork => ListGetConstRealArray(CurrentModel % Constants, 'Gravity', gotgrav)
      acc_grav = 0.0_dp
      acc_gotgrav = .false.
      if (gotgrav) then
        acc_gotgrav = .true.
        acc_grav = abs(gwork(size(gwork, 1), 1))
        if (.not. found) usegravity = .true.
      else
        if (usegravity) call Fatal('FilmFlowSolver', &
            'Gravity requested but not given as constant!')
      end if
    end block


    if (any(acc_fricmodel == [3, 4])) then
      if (.not. acc_gotgrav) call Fatal('FilmFlowSolver', &
          'Manning equation not possible without gravity!')
      if (acc_cs) call Fatal('FilmFlowSolver', &
          'Manning equation not applicable to axial symmetry!')
    end if

    acc_calcfric = ListGetLogical(params, 'Calculate Friction Heating', found)
    acc_calcpres = ListGetLogical(params, 'Calculate Pressure Heating', found)
    acc_calcheat = acc_calcfric .or. acc_calcpres
    acc_useheat = ListGetLogical(params, 'Use Heating Source', found)

    if (.not. acc_allocdone) then
      call Info('FilmFlowSolver', &
          'Dimension of Navier-Stokes equation: '//I2S(acc_mdim))
      call Info('FilmFlowSolver', &
          'Dimension of coordinate system: '//I2S(acc_dim))
      n = (acc_mdim + 1) * (acc_mesh % MaxElementDOFs + 4 * 1)
      allocate (acc_force(n), acc_load(acc_mdim + 4, n), acc_stiff(n, n), &
          acc_mass(n, n), acc_rho(n), acc_ac(n), acc_gap(n), acc_gap0(n), &
          acc_height(n), acc_mu(n), acc_acpres(n), acc_velo(acc_mdim + 1, n), &
          stat=istat)
      acc_velo = 0.0_dp
      if (istat /= 0) then
        call Fatal('FilmFlowSolver', 'Memory allocation error.')
      end if
      if (acc_gradp) then
        call Info('FilmFlowSolver', '"Gradp Discretization" is set True', level=10)
      end if
      acc_pvar => VariableGet(acc_mesh % Variables, 'FilmPressure')
      if (.not. associated(acc_pvar)) then
        call Fatal('FilmFlowSolver', 'Could not find required field "FilmPressure"!')
      end if
      n = size(acc_pvar % Values)
      if (acc_gotac) then
        call Info('FilmFlowSolver', &
            'Using artificial compressibility for FSI emulation!')
        allocate (acc_acprevp(n), acc_fsirhs(2, n))
        acc_acprevp = 0.0_dp
        acc_fsirhs = 0.0_dp
        acc_surfac = ListGetLogical(params, 'Surface Compressibility', found)
      end if
      if (acc_useprevgap) then
        allocate (acc_prevgap(n))
        acc_prevgap = 0.0_dp
      end if
      if (acc_calcfric) then
        call Info('FilmFlowSolver', 'Registering friction heat flux', level=7)
        hvar => VariableGet(acc_mesh % Variables, 'FilmFriction HeatFlux')
        acc_hvar => hvar
        acc_fricflux => hvar % Values
      end if
      if (acc_calcpres) then
        call Info('FilmFlowSolver', 'Registering convected heat flux', level=7)
        hvar => VariableGet(acc_mesh % Variables, 'FilmPressure HeatFlux')
        acc_hvar => hvar
        acc_presflux => hvar % Values
      end if
      if (acc_calcheat) then
        call Info('FilmFlowSolver', 'Registering total energy production', level=7)
        hvar => VariableGet(acc_mesh % Variables, 'Heating Energy')
        acc_hvar => hvar
        acc_heatingenergy => hvar % Values
        allocate (acc_heatingw(size(acc_heatingenergy)))
      end if
      acc_allocdone = .true.
      if (ListGetLogical(params, 'Skip First Solution', found)) then
        call Info('FilmFlowSolver', &
            'Skipping first solution phase completely!', level=5)
        code = 1
        return
      end if
    end if
    code = 0
  end subroutine ffw_begin

  ! -- pre-sweep heating correction + timestep -------------------------
  subroutine ffw_heating_pre() bind(c, name="ffw_heating_pre")
    implicit none
    acc_itime = GetTimestep()
    if (acc_calcheat) then
      if (acc_itime - acc_itime0 == 0) then
        if (acc_calcfric .and. acc_calcpres) then
          acc_heatingenergy = acc_heatingenergy - acc_dt * &
              max(acc_fricflux + acc_presflux, 0.0_dp)
        else
          if (acc_calcfric) then
            acc_heatingenergy = acc_heatingenergy - acc_dt * acc_fricflux
          end if
          if (acc_calcpres) then
            acc_heatingenergy = acc_heatingenergy - acc_dt * &
                max(acc_presflux, 0.0_dp)
          end if
        end if
      end if
      if (acc_calcfric) acc_fricflux = 0.0_dp
      if (acc_calcpres) acc_presflux = 0.0_dp
      if (.not. acc_doneweight) acc_heatingw = 0.0_dp
      acc_itime0 = acc_itime
    end if
  end subroutine ffw_heating_pre

  ! -- FSI pressure memorize --------------------------------------------
  subroutine ffw_ac_begin() bind(c, name="ffw_ac_begin")
    implicit none
    if (acc_gotac) then
      acc_acprevp = acc_pvar % Values
      if (InfoActive(12)) then
        write (message, '(A,2E15.4)') 'PrevPressure: ', &
            sum(acc_acprevp) / size(acc_acprevp), maxval(acc_acprevp)
        call Info('FilmFlowSolver', message, level=5)
      end if
      acc_fsirhs = 0.0_dp
    end if
  end subroutine ffw_ac_begin

  ! -- nonlinear-loop setup: maxiter + FSI multiplier -------------------
  function ffw_maxiter() result(m) bind(c, name="ffw_maxiter")
    implicit none
    integer(c_int) :: m
    logical :: found
    m = ListGetInteger(acc_params, 'Nonlinear System Max Iterations', &
        found, minv=1)
    if (.not. found) m = 1
    if (acc_calcheat) m = max(2, m)
    acc_maxiter = m
    acc_fsimult = ListGetCReal(acc_params, 'Fsi Velocity Multiplier', found)
    if (.not. found) acc_fsimult = 1.0_dp
    if (acc_itime == 1) then
      if (ListGetLogical(acc_params, 'Steady State Fsi Start', found)) then
        call Info('FilmFlowSolver', &
            'Ignoring FSI velocity for the 1st timestep', level=7)
        acc_fsimult = 0.0_dp
      end if
    end if
  end function ffw_maxiter

  ! -- per-iteration setup; returns GetNOFActive() ----------------------
  function ffw_iter_begin() result(active) bind(c, name="ffw_iter_begin")
    implicit none
    integer(c_int) :: active
    logical :: found
    real(kind=dp) :: s
    call DefaultInitialize()
    acc_newton = GetNewtonActive()
    s = ListGetCReal(acc_params, 'Convect Condition', found)
    if (found) acc_convect = (s > 0.0_dp)
    active = GetNOFActive()
    acc_active = active
  end function ffw_iter_begin

  ! -- bulk element body (exact transcription) ---------------------------
  subroutine ffw_bulk_elem(t) bind(c, name="ffw_bulk_elem")
    implicit none
    integer(c_int), value :: t
    integer :: n, nd, nb
    logical :: found
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
    nb = GetElementNOFBDOFs()
    acc_bf => GetBodyForce()
    acc_load = 0.0d0
    if (associated(acc_bf)) then
      acc_load(1, 1:n) = GetReal(acc_bf, 'FilmFlow Bodyforce 1', found)
      if (acc_mdim > 1) &
          acc_load(2, 1:n) = GetReal(acc_bf, 'FilmFlow Bodyforce 2', found)
      acc_load(acc_mdim + 1, 1:n) = GetReal(acc_bf, 'Normal Velocity', found)
      acc_load(acc_mdim + 2, 1:n) = &
          acc_fsimult * GetReal(acc_bf, 'Fsi Velocity', found)
      acc_load(acc_mdim + 3, 1:n) = GetReal(acc_bf, 'Flow Admittance', found)
      acc_load(acc_mdim + 4, 1:n) = &
          GetReal(acc_bf, 'External FilmPressure', found)
    end if
    acc_mat => GetMaterial()
    acc_rho(1:n) = GetReal(acc_mat, acc_densityname)
    acc_mu(1:n) = GetReal(acc_mat, acc_viscname)
    acc_gap(1:n) = GetReal(acc_mat, 'Gap Height')
    acc_height(1:n) = GetReal(acc_mat, 'Bedrock Height', acc_gotheight)
    if (any(acc_fricmodel == [1, 2])) then
      acc_nm = ListGetCReal(acc_mat, 'Darcy Roughness', UnfoundFatal=.true.)
    else if (any(acc_fricmodel == [3, 4])) then
      acc_nm = ListGetCReal(acc_mat, 'Manning coefficient', UnfoundFatal=.true.)
    end if
    if (acc_useheat) then
      acc_meltheat = &
          ListGetConstReal(acc_mat, 'Latent Heat', UnFoundFatal=.true.)
    end if
    if (acc_calcpres) then
      acc_cp = ListGetConstReal(acc_mat, 'Water Heat Capacity', &
          UnFoundFatal=.true.)
      acc_ct = ListGetConstReal(acc_mat, 'Pressure Melting Coefficient', &
          UnFoundFatal=.true.)
    end if
    where (acc_gap(1:n) < acc_mingap)
      acc_gap(1:n) = acc_mingap
    end where
    if (acc_useprevgap) then
      if (acc_coupled == 1) then
        acc_prevgap(acc_pvar % Perm(acc_elem % NodeIndexes)) = acc_gap(1:n)
      end if
      acc_gap0(1:n) = acc_prevgap(acc_pvar % Perm(acc_elem % NodeIndexes))
    end if
    if (acc_gotac) then
      acc_ac(1:n) = GetReal(acc_mat, 'Artificial Compressibility', found)
      acc_acpres(1:n) = acc_acprevp(acc_pvar % Perm(acc_elem % NodeIndexes))
    end if
    call GetVectorLocalSolution(acc_velo)
    call LocalBulkMatrix(acc_mass, acc_stiff, acc_force, acc_load, acc_rho, &
        acc_gap, acc_gap0, acc_height, acc_mu, acc_ac, acc_velo, acc_acpres, &
        acc_elem, n, nd, nd + nb, acc_dim, acc_mdim, acc_firstround)
    if (nb > 0) then
      call NSCondensate(nd, nb, acc_mdim, acc_stiff, acc_force)
    end if
    if (acc_transient) then
      call Default1stOrderTime(acc_mass, acc_stiff, acc_force)
    end if
    call DefaultUpdateEquations(acc_stiff, acc_force)
  end subroutine ffw_bulk_elem

  subroutine ffw_finish_bulk() bind(c, name="ffw_finish_bulk")
    implicit none
    call DefaultFinishBulkAssembly()
  end subroutine ffw_finish_bulk

  ! -- FSI flux report ----------------------------------------------------
  subroutine ffw_ac_report() bind(c, name="ffw_ac_report")
    implicit none
    if (acc_gotac) then
      block
        real(kind=dp) :: sorig, sfsi, coeff
        sorig = sum(acc_fsirhs(1, :))
        sfsi = sum(acc_fsirhs(2, :))
        coeff = 1.0_dp
        if (sfsi /= 0.0) coeff = sorig / sfsi
        if (sfsi < sorig) coeff = 1.0_dp
        if (InfoActive(20)) then
          print *, 'RHSComp:', sorig, sfsi, coeff
        end if
      end block
    end if
  end subroutine ffw_ac_report

  function ffw_nbnd() result(n) bind(c, name="ffw_nbnd")
    implicit none
    integer(c_int) :: n
    n = acc_solver % Mesh % NumberOfBoundaryElements
  end function ffw_nbnd

  ! -- boundary element body (exact transcription) ------------------------
  subroutine ffw_boundary_elem(t) bind(c, name="ffw_boundary_elem")
    implicit none
    integer(c_int), value :: t
    integer :: n, i
    logical :: found
    acc_elem => GetBoundaryElement(t)
    if (.not. ActiveBoundaryElement()) return
    n = GetElementNOFNodes()
    acc_bc => GetBC()
    if (.not. associated(acc_bc)) return
    acc_rho(1:n) = GetParentMatProp(acc_densityname, acc_elem, found)
    acc_mu(1:n) = GetParentMatProp(acc_viscname, acc_elem, found)
    acc_gap(1:n) = GetParentMatProp('Gap Height', acc_elem, found)
    where (acc_gap(1:n) < acc_mingap)
      acc_gap(1:n) = acc_mingap
    end where
    if (acc_useprevgap) then
      acc_gap0(1:n) = acc_prevgap(acc_pvar % Perm(acc_elem % NodeIndexes))
    end if
    do i = 1, acc_mdim
      acc_load(i, 1:n) = GetReal(acc_bc, 'FilmPressure '//I2S(i), found)
    end do
    acc_load(acc_mdim + 1, 1:n) = GetReal(acc_bc, 'Mass Flux', found)
    call LocalBoundaryMatrix(acc_mass, acc_stiff, acc_force, acc_load, &
        acc_rho, acc_gap, acc_mu, acc_elem, n, acc_dim, acc_mdim)
    call DefaultUpdateEquations(acc_stiff, acc_force)
  end subroutine ffw_boundary_elem

  ! -- iteration tail: finish assembly + solve ----------------------------
  function ffw_iter_end() result(norm) bind(c, name="ffw_iter_end")
    implicit none
    real(c_double) :: norm
    call DefaultFinishBoundaryAssembly()
    call DefaultFinishAssembly()
    call DefaultDirichletBCs()
    acc_norm = DefaultSolve()
    norm = acc_norm
  end function ffw_iter_end

  function ffw_nonlin_converged() result(r) bind(c, name="ffw_nonlin_converged")
    implicit none
    integer(c_int) :: r
    r = merge(1, 0, acc_solver % Variable % NonlinConverged == 1)
  end function ffw_nonlin_converged

  subroutine ffw_firstround_done() bind(c, name="ffw_firstround_done")
    implicit none
    acc_firstround = .false.
  end subroutine ffw_firstround_done

  function ffw_sensitivity() result(r) bind(c, name="ffw_sensitivity")
    implicit none
    integer(c_int) :: r
    r = merge(1, 0, DefaultSensitivity())
  end function ffw_sensitivity

  subroutine ffw_finish() bind(c, name="ffw_finish")
    implicit none
    call DefaultFinish()
  end subroutine ffw_finish

  ! -- post-sweep heating (exact transcription) ------------------------------
  subroutine ffw_heating_post() bind(c, name="ffw_heating_post")
    implicit none
    logical :: found
    if (acc_calcheat) then
      block
        real(kind=dp) :: totflux, area, cfix
        totflux = 0.0_dp
        if (acc_calcfric) &
            totflux = ParallelReduction(sum(acc_fricflux))
        if (acc_calcpres) &
            totflux = totflux + ParallelReduction(sum(acc_fricflux))
        if (ParEnv % PEs > 1) then
          if (acc_calcfric) &
              call ParallelSumNodalVector(acc_mesh, acc_fricflux, &
              acc_hvar % Perm)
          if (acc_calcpres) &
              call ParallelSumNodalVector(acc_mesh, acc_presflux, &
              acc_hvar % Perm)
          if (.not. acc_doneweight) &
              call ParallelSumNodalVector(acc_mesh, acc_heatingw, &
              acc_hvar % Perm)
        end if
        acc_doneweight = .true.
        if (acc_calcfric) then
          where (acc_heatingw > 0.0_dp)
            acc_fricflux = acc_fricflux / acc_heatingw
          end where
        end if
        if (acc_calcpres) then
          where (acc_heatingw > 0.0_dp)
            acc_presflux = acc_presflux / acc_heatingw
          end where
        end if
        if (ListGetLogical(acc_params, 'Enforce Zero Heating', found)) then
          area = ParallelReduction(sum(acc_heatingw))
          cfix = -totflux / area
          acc_fricflux = acc_fricflux + cfix
        end if
        if (acc_calcfric .and. acc_calcpres) then
          acc_heatingenergy = acc_heatingenergy + acc_dt * &
              max(acc_fricflux + acc_presflux, 0.0_dp)
        else
          if (acc_calcfric) then
            acc_heatingenergy = acc_heatingenergy + acc_dt * acc_fricflux
          end if
          if (acc_calcpres) then
            acc_heatingenergy = acc_heatingenergy + acc_dt * &
                max(acc_presflux, 0.0_dp)
          end if
        end if
        write (message, '(A,ES12.5)') 'Total heating power: ', totflux
        call Info('FilmFlowSolver', message, level=7)
      end block
    end if
  end subroutine ffw_heating_post

  ! -- diagnostics -------------------------------------------------------------
  subroutine ffw_diagnostics() bind(c, name="ffw_diagnostics")
    implicit none
    integer :: i, n
    real(kind=dp), pointer :: comp(:)
    if (InfoActive(12)) then
      n = size(acc_pvar % Values)
      do i = 1, acc_solver % Variable % dofs
        comp => acc_solver % Variable % Values(i::acc_solver % Variable % Dofs)
        call VectorValuesRange(comp, n, 'Velocity '//I2S(i))
      end do
      call VectorValuesRange(acc_pvar % Values, n, 'Pressure')
      if (acc_gotac) then
        call VectorValuesRange(acc_acprevp, n, 'Pressure0')
        call VectorValuesRange(acc_pvar % Values - acc_acprevp, n, &
            'PressureDiff '//I2S(acc_coupled))
      end if
    end if
    call Info('FilmFlowSolver', 'All done', level=12)
  end subroutine ffw_diagnostics

  !---------------------------------------------------------------------------!
  ! CONTAINS internals, transcribed with host association redirected at the   !
  ! module SAVE mirrors above.                                                !
  !---------------------------------------------------------------------------!

  ! Eq (29) in Praks & Brkic (2020); the original's SAVE Visited flag is
  ! set-but-never-read and is not replicated (no observable state).
  function FrictionLawPraks(v, rho, nu, D, eps) result(f)
    implicit none
    real(kind=dp) :: v, rho, nu, D, eps, f
    real(kind=dp) :: Re, A, B, C, x
    real(kind=dp) :: lambda
    Re = v * (D / 2) * rho / nu
    A = Re * eps / 8.0897_dp
    B = log(Re) - 0.779626_dp
    x = A + B
    if (x <= AEPS) then
      f = 64._dp / Re
      return
    end if
    C = log(x)
    lambda = Re * (6.94871_dp * (B - C + C / (x - 0.5588_dp * C + 1.2079_dp)))**(-2.0_dp)
    lambda = max(1.0_dp, lambda)
    f = 64.0_dp * lambda / Re
  end function FrictionLawPraks

  subroutine LocalBulkMatrix(MASS, STIFF, FORCE, LOAD, Nodalrho, NodalGap, &
      NodalGap0, NodalH, Nodalmu, NodalAC, NodalVelo, NodalAcPres, &
      Element, n, nd, ntot, dim, mdim, FirstRound)
    implicit none
    real(kind=dp), target :: MASS(:,:), STIFF(:,:), FORCE(:), LOAD(:,:)
    real(kind=dp) :: Nodalmu(:), NodalAC(:), Nodalrho(:), &
        NodalGap(:), NodalGap0(:), NodalH(:), NodalAcPres(:), NodalVelo(:,:)
    integer :: dim, mdim, n, nd, ntot
    type(Element_t), target :: Element
    logical :: FirstRound
    real(kind=dp) :: Basis(ntot), dBasisdx(ntot, 3)
    real(kind=dp) :: DetJ, LoadAtIP(mdim + 4), Velo(mdim), VeloGrad(mdim, mdim), &
        gapGrad(mdim), hGrad(mdim), presGrad(mdim)
    real(kind=dp) :: NodalPres(n), NodalPrevPres(n)
    real(kind=dp), pointer :: A(:,:), F(:), M(:,:)
    logical :: Stat, Found
    integer :: t, i, j, p, q, geomc
    type(GaussIntegrationPoints_t) :: IP
    real(kind=dp) :: mu = 1.0d0, rho = 1.0d0, AcPres, gap, gap0, gap2, gapi, &
        ac, s, s0, s1, MinPres, MuCoeff, MinSpeed, Speed, h, q_p, q_f, Pres, &
        PrevPres, FlowAdm


    call GetElementNodes(acc_nodes_bulk)

    if (acc_gapdir > 0) then
      select case (acc_gapdir)
      case (1)
        acc_nodes_bulk % x(1:n) = &
            acc_nodes_bulk % x(1:n) + acc_gapfactor * NodalGap(1:n)
      case (2)
        acc_nodes_bulk % y(1:n) = &
            acc_nodes_bulk % y(1:n) + acc_gapfactor * NodalGap(1:n)
      case (3)
        acc_nodes_bulk % z(1:n) = &
            acc_nodes_bulk % z(1:n) + acc_gapfactor * NodalGap(1:n)
      end select
    end if

    call GetLocalSolution(NodalPres, UElement=Element, UVariable=acc_pvar)
    call GetLocalSolution(NodalPrevPres, UElement=Element, &
        UVariable=acc_pvar, tStep=-1)

    STIFF = 0.0d0
    MASS = 0.0d0
    FORCE = 0.0d0
    gapGrad = 0.0_dp
    hGrad = 0.0_dp
    presGrad = 0.0_dp
    q_f = 0.0_dp
    q_p = 0.0_dp

    if (acc_cs) then
      geomc = 2
    else
      geomc = 1
    end if

    IP = GaussPointsAdapt(Element, PReferenceElement=.true.)

    if (.not. acc_bvisited) then
      call Info('FilmFlowSolver', &
          'Number of integration points: '//I2S(IP % n))
      MinPres = ListGetConstReal(acc_params, 'Min FilmPressure', Found)
      if (.not. Found) MinPres = -huge(MinPres)
      MinSpeed = ListGetConstReal(acc_params, 'Min Speed', Found)
      if (.not. Found) MinSpeed = 1.0e-6
      acc_minpres = MinPres
      acc_minspeed = MinSpeed
      acc_bvisited = .true.
    end if
    MinPres = acc_minpres
    MinSpeed = acc_minspeed

    do t = 1, IP % n
      Stat = ElementInfo(Element, acc_nodes_bulk, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx)
      s = IP % s(t) * detJ
      s1 = s
      s0 = s
      mu = sum(Basis(1:n) * Nodalmu(1:n))
      rho = sum(Basis(1:n) * Nodalrho(1:n))
      gap = sum(Basis(1:n) * NodalGap(1:n))
      gap0 = sum(Basis(1:n) * NodalGap0(1:n))
      AcPres = max(MinPres, sum(NodalAcPres(1:n) * Basis(1:n)))
      Pres = sum(NodalPres(1:n) * Basis(1:n))
      do i = 1, mdim
        gapGrad(i) = sum(NodalGap(1:nd) * dBasisdx(1:nd, i))
        presGrad(i) = sum(NodalPres(1:nd) * dBasisdx(1:nd, i))
      end do
      if (acc_gotheight) then
        h = sum(NodalH(1:nd) * Basis(1:nd))
        do i = 1, mdim
          hGrad(i) = sum(NodalH(1:nd) * dBasisdx(1:nd, i))
        end do
      else
        if (mdim == 1) then
          h = sum(acc_nodes_bulk % y(1:n) * Basis(1:n))
          hGrad(1) = sum(acc_nodes_bulk % y(1:n) * dBasisdx(1:n, 1))
        else
          h = sum(acc_nodes_bulk % z(1:n) * Basis(1:n))
          hGrad(1) = sum(acc_nodes_bulk % z(1:n) * dBasisdx(1:n, 1))
          hGrad(2) = sum(acc_nodes_bulk % z(1:n) * dBasisdx(1:n, 2))
        end if
      end if
      Velo = matmul(NodalVelo(1:mdim, 1:nd), Basis(1:nd))
      VeloGrad = matmul(NodalVelo(1:mdim, 1:nd), dBasisdx(1:nd, 1:mdim))
      Speed = sqrt(sum(Velo(1:mdim)**2))
      if (acc_gotac) then
        ac = sum(NodalAC(1:n) * Basis(1:n)) / acc_dt
      end if
      do i = 1, mdim + 4
        LoadAtIP(i) = sum(Basis(1:n) * LOAD(i, 1:n))
      end do
      if (acc_convect .and. acc_newton) then
        LoadAtIp(1:mdim) = LoadAtIp(1:mdim) + &
            rho * matmul(VeloGrad(1:mdim, 1:mdim), Velo(1:mdim))
      end if
      LoadAtIp(mdim + 1:mdim + 2) = geomc * LoadAtIp(mdim + 1:mdim + 2)
      if (acc_useprevgap) then
        gap2 = gap * gap0
        gapi = sqrt(gap2)
      else
        gap2 = gap**2
        gapi = gap
      end if
      select case (acc_fricmodel)
      case (1, 2)
        block
          real(kind=dp) :: D, R, fd, GradZphi2
          Speed = max(MinSpeed, Speed)
          D = 2 * gapi
          fd = FrictionLawPraks(Speed, rho, mu, D, acc_nm)
          if (acc_fricmodel == 1) then
            MuCoeff = fd * rho * Speed / (2 * D)
          else
            GradZphi2 = max(sum((hGrad(1:mdim) + &
                presGrad(1:mdim) / (rho * acc_grav))**2), 1.0E-09)
            MuCoeff = rho * sqrt(fd * acc_grav) * (2 * gap)**(-1.0 / 2.0) * &
                GradZphi2**(1.0 / 4.0)
          end if
        end block
      case (3)
        block
          real(kind=dp) :: GradZphi2
          GradZphi2 = max(sum((hGrad(1:mdim) + &
              presGrad(1:mdim) / (rho * acc_grav))**2), 1.0E-09)
          MuCoeff = acc_nm * rho * acc_grav * (gapi / 2)**(-2.0 / 3) * &
              GradZphi2**(1.0 / 4.0)
        end block
      case (4)
        block
          MuCoeff = rho * acc_grav * acc_nm**2 * (gapi / 2)**(-4.0 / 3) * Speed
        end block
      case default
        if (acc_cs) then
          MuCoeff = 8 * mu / gap2
        else
          MuCoeff = 12 * mu / gap2
        end if
      end select
      if (acc_calcfric) then
        q_f = MuCoeff * gapi * Speed**2
      end if
      if (acc_calcpres) then
        PrevPres = sum(NodalPrevPres(1:n) * Basis(1:n))
        q_p = rho * gapi * acc_cp * acc_ct * ((Pres - PrevPres) / acc_dt + &
            sum(Velo(1:mdim) * presGrad(1:mdim)))
      end if
      if (FirstRound) then
        q_f = 0.0_dp
        q_p = 0.0_dp
      end if
      do p = 1, ntot
        do q = 1, ntot
          i = (mdim + 1) * (p - 1) + 1
          j = (mdim + 1) * (q - 1) + 1
          A => STIFF(i:i + mdim, j:j + mdim)
          M => MASS(i:i + mdim, j:j + mdim)
          do i = 1, mdim
            if (acc_transient) then
              M(i, i) = M(i, i) + s * rho * Basis(q) * Basis(p)
            end if
            A(i, i) = A(i, i) + s * MuCoeff * Basis(q) * Basis(p)
            do j = 1, mdim
              if (acc_lateral) then
                A(i, i) = A(i, i) + s * mu * dBasisdx(q, j) * dBasisdx(p, j)
                A(i, j) = A(i, j) + s * mu * dBasisdx(q, i) * dBasisdx(p, j)
              end if
              if (acc_convect) then
                A(i, i) = A(i, i) + s * rho * Velo(j) * dBasisdx(q, j) * Basis(p)
                if (acc_newton) then
                  A(i, j) = A(i, j) + s * rho * VeloGrad(i, j) * &
                      Basis(q) * Basis(p)
                end if
              end if
            end do
            if (acc_gradp) then
              A(i, mdim + 1) = A(i, mdim + 1) + s * dBasisdx(q, i) * Basis(p)
              A(mdim + 1, i) = A(mdim + 1, i) - s * gap * rho * &
                  Basis(q) * dBasisdx(p, i)
            else
              A(i, mdim + 1) = A(i, mdim + 1) - s * Basis(q) * dBasisdx(p, i)
              A(mdim + 1, i) = A(mdim + 1, i) + s * gap * rho * &
                  dBasisdx(q, i) * Basis(p) + &
                  geomc * s * rho * Basis(q) * gapGrad(i) * Basis(p)
            end if
          end do
          if (acc_gotac) A(mdim + 1, mdim + 1) = A(mdim + 1, mdim + 1) + &
              ac * s * rho * Basis(q) * Basis(p)
          A(mdim + 1, mdim + 1) = A(mdim + 1, mdim + 1) + s * rho * &
              LoadAtIP(mdim + 3) * Basis(q) * Basis(p)
        end do
        i = (mdim + 1) * (p - 1) + 1
        F => FORCE(i:i + mdim)
        F(1:mdim + 1) = F(1:mdim + 1) - s * rho * Basis(p) * LoadAtIp(1:mdim + 1)
        F(mdim + 1) = F(mdim + 1) - s * rho * Basis(p) * LoadAtIp(mdim + 2)
        if (acc_gotgrav) then
          F(1:mdim) = F(1:mdim) - s * rho * acc_grav * Basis(p) * hGrad(1:mdim)
        end if
        if (acc_gotac) F(mdim + 1) = F(mdim + 1) + &
            ac * s * rho * Basis(p) * AcPres
        F(mdim + 1) = F(mdim + 1) + s * rho * Basis(p) * &
            LoadAtIp(mdim + 3) * LoadAtIP(mdim + 4)
        if (acc_useheat) then
          F(1:mdim + 1) = F(1:mdim + 1) + s * Basis(p) * &
              max(q_f + q_p, 0.0_dp) / (acc_meltheat * rho)
        end if
      end do
      if (acc_calcfric) then
        acc_fricflux(acc_hvar % Perm(Element % NodeIndexes)) = &
            acc_fricflux(acc_hvar % Perm(Element % NodeIndexes)) + &
            s * Basis(1:n) * q_f
      end if
      if (acc_calcpres) then
        acc_presflux(acc_hvar % Perm(Element % NodeIndexes)) = &
            acc_presflux(acc_hvar % Perm(Element % NodeIndexes)) + &
            s * Basis(1:n) * q_p
      end if
      if (acc_calcheat .and. .not. acc_doneweight) then
        acc_heatingw(acc_hvar % Perm(Element % NodeIndexes)) = &
            acc_heatingw(acc_hvar % Perm(Element % NodeIndexes)) + &
            s * Basis(1:n)
      end if
      if (acc_gotac) then
        acc_fsirhs(1, acc_thisvar % Perm(Element % NodeIndexes)) = &
            acc_fsirhs(1, acc_thisvar % Perm(Element % NodeIndexes)) + &
            s * rho * LoadAtIp(mdim + 1) * Basis(1:n)
        acc_fsirhs(2, acc_thisvar % Perm(Element % NodeIndexes)) = &
            acc_fsirhs(2, acc_thisvar % Perm(Element % NodeIndexes)) + &
            s * rho * LoadAtIp(mdim + 2) * Basis(1:n)
      end if
    end do
    do p = n + 1, ntot
      i = (mdim + 1) * p
      FORCE(i) = 0.0d0
      MASS(:, i) = 0.0d0
      MASS(i, :) = 0.0d0
      STIFF(i, :) = 0.0d0
      STIFF(:, i) = 0.0d0
      STIFF(i, i) = 1.0d0
    end do
  end subroutine LocalBulkMatrix

  subroutine LocalBoundaryMatrix(MASS, STIFF, FORCE, NodalLoad, Nodalrho, &
      NodalGap, Nodalmu, Element, n, dim, mdim)
    implicit none
    real(kind=dp), target :: MASS(:,:), STIFF(:,:), FORCE(:), NodalLoad(:,:)
    real(kind=dp) :: Nodalmu(:), Nodalrho(:), NodalGap(:)
    integer :: dim, mdim, n
    type(Element_t), target :: Element
    real(kind=dp) :: Basis(n), dBasisdx(n, 3), DetJ, Load(mdim + 1)
    real(kind=dp), pointer :: F(:)
    logical :: Stat
    integer :: t, i, p, geomc
    type(GaussIntegrationPoints_t) :: IP
    real(kind=dp) :: mu = 1.0d0, rho = 1.0d0, gap, s
    call GetElementNodes(acc_nodes_bnd)
    STIFF = 0.0d0
    MASS = 0.0d0
    FORCE = 0.0d0
    if (acc_cs) then
      geomc = 2
    else
      geomc = 1
    end if
    IP = GaussPoints(Element)
    do t = 1, IP % n
      Stat = ElementInfo(Element, acc_nodes_bnd, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx)
      s = IP % s(t) * detJ
      mu = sum(Basis(1:n) * Nodalmu(1:n))
      rho = sum(Basis(1:n) * Nodalrho(1:n))
      gap = sum(Basis(1:n) * NodalGap(1:n))
      do i = 1, mdim + 1
        Load(i) = sum(Basis(1:n) * NodalLoad(i, 1:n))
      end do
      do p = 1, n
        i = (mdim + 1) * (p - 1) + 1
        F => FORCE(i:i + mdim)
        F = F + s * Basis(p) * Load
      end do
    end do
  end subroutine LocalBoundaryMatrix

end module filmflowsolver_acc
