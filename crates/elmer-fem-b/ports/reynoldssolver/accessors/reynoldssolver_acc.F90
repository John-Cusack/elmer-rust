! T021-B port of fem/src/modules/ReynoldsSolver.F90: statement-level
! BIND(C) phase wrappers for the Elmer calls in the unit, compiled against
! the built .mod files of the target configure (the Poisson-template
! mechanism, ADR-057 section 3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape out integers). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL dummies,
! assumed-shape arrays, array-pointer results, keyword actuals
! (UVariable=, Ux=, tstep=), hidden lengths — stays inside Fortran, so
! the Rust side is pure C calls. The wrapper modules mirror the
! original's SAVE locals with identical lifetimes (module SAVE state,
! allocated once at first-mesh sizing, never resized).
!
! Shape of the port (bit-identity by construction where the call
! sequence matches):
! - `reynolds_solve_acc`: mirrors ReynoldsSolver. The whole-bulk and
!   whole-boundary assemblies (GlobalBulkAssembly, LocalBulkMatrix,
!   GlobalBoundaryAssemby, LocalBoundaryMatrix) are transcribed
!   statement-for-statement from the original's CONTAINS section, so
!   order-sensitive SUM reductions, integer-power/exponential
!   evaluations and whole-array updates run under the same compiler in
!   the same order. Rust owns the entry-level control flow (nonlinear
!   iteration loop, linesearch redo, sensitivity redo, gap-sensitivity
!   gate) and calls one phase wrapper per statement block.
! - `reynolds_post_acc`: mirrors ReynoldsPostprocess the same way; the
!   per-component element loop and LocalMatrix kernel are transcribed
!   verbatim (including the formatted WRITE/Info force-moment-heating
!   tails and the strided `Values(Component::Components)` stores).
! - The two `_init` entries are pure list manipulation whose exact
!   keyword/I2S behaviour is kept verbatim in single coarse wrappers
!   (no floating point involved).
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls plus transcriptions of the unit's own kernels; no upstream text
! beyond the call/shape facts any caller must use.
module reynolds_solve_acc
  use iso_c_binding
  use DefUtils
  use Types
  use Lists
  use Integration
  use ElementDescription
  use SolverBasics
  implicit none

  ! -- mirrored SAVE state of ReynoldsSolver (identical lifetimes) --
  type(Nodes_t) :: ElementNodes
  type(Element_t), pointer :: Element => null(), Parent => null()
  type(ValueList_t), pointer :: Params => null(), Material => null(), &
      Equation => null(), BC => null(), BodyForce => null()
  type(Variable_t), pointer :: SensVar => null(), PresVar => null()
  type(Model_t), pointer :: acc_model => null()
  type(Solver_t), pointer :: acc_solver => null()
  integer :: iter, i, j, k, l, n, nd, t, istat, mat_id, ent_id, body_id, &
      mat_idold, NoIterations, ViscosityType, CompressibilityType, itime
  integer, pointer :: NodeIndexes(:), PressurePerm(:)
  logical :: GotIt, GotIt2, GotIt3, stat, &
      AllocationsDone = .FALSE., SubroutineVisited = .FALSE., &
      UseVelocity, Bubbles, ApplyLimiter, LinearModel, ManningModel, &
      GotMinGap, OpenSide, GotExt, GotFlux, GotVelo, AnyBC, &
      GotPseudoPressure, Converged
  real(kind=dp), pointer :: Pressure(:)
  real(kind=dp) :: Norm, ReferencePressure, HeatRatio, BulkModulus, &
      mfp0, Pres, Dens, ManningCoeff, GravityCoeff, MinGap, MinGradPres, &
      ACScale, FSIMult
  real(kind=dp), allocatable :: STIFF(:,:), MASS(:,:), FORCE(:), &
      TimeForce(:), Viscosity(:), GapHeight(:), NormalVelocity(:), &
      Velocity(:,:), Admittance(:), Impedance(:), ElemPressure(:), &
      PrevElemPressure(:), ExtPressure(:), ElemDensity(:), ElemArtif(:), &
      ExtPres(:), FluxPres(:), VeloPres(:), CoeffPres(:), &
      ElemPseudoPressure(:), PseudoPressure(:)
  character(len=MAX_NAME_LEN) :: ViscosityModel, CompressibilityModel, varname
  character(*), parameter :: Caller = 'ReynoldsSolver'
  integer, parameter :: Compressibility_None = 1, Compressibility_Weak = 2, &
      Compressibility_GasIsothermal = 3, Compressibility_GasAdiabatic = 4, &
      Compressibility_Artificial = 5
  integer, parameter :: Viscosity_Newtonian = 1, Viscosity_Rarefied = 2
  real(kind=dp) :: dt_mod
  logical :: transient_mod

contains

  subroutine rrs_frame_begin(h_model, h_solver, dt, transient) &
      bind(c, name="rrs_frame_begin")
    implicit none
    type(c_ptr), value :: h_model, h_solver
    real(c_double), value :: dt
    integer(c_int), value :: transient
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    dt_mod = dt
    transient_mod = (transient /= 0)
  end subroutine rrs_frame_begin

  subroutine rrs_frame_end() bind(c, name="rrs_frame_end")
    implicit none
    acc_model => null()
    acc_solver => null()
    Element => null()
    Parent => null()
    Params => null()
    Material => null()
    Equation => null()
    BC => null()
    BodyForce => null()
    SensVar => null()
    PresVar => null()
    Pressure => null()
    PressurePerm => null()
    NodeIndexes => null()
  end subroutine rrs_frame_end

  ! Lines 99-153 of the original: header Info, DefaultStart, parameter
  ! fetch, early-return gates. Copies NoIterations out for the Rust
  ! iteration loop; proceed=0 reproduces the bare RETURNs.
  subroutine rrs_preamble(noiter, proceed) bind(c, name="rrs_preamble")
    implicit none
    integer(c_int) :: noiter, proceed
    CALL Info(Caller,'---------------------------------------',Level=5)
    IF(transient_mod) THEN
      CALL Info(Caller,'Solving the transient Reynolds equation',Level=5)
    ELSE
      CALL Info(Caller,'Solving the steady-state Reynolds equation',Level=5)
    END IF
    CALL Info(Caller,'---------------------------------------',Level=5)
    CALL DefaultStart()
    Params => GetSolverParams()
    Bubbles = GetLogical( Params, 'Bubbles', GotIt )
    IF ( .NOT. ASSOCIATED( acc_solver % Matrix ) ) THEN
      proceed = 0
      noiter = 0
      RETURN
    END IF
    IF(acc_solver % Variable % Dofs /= 1) THEN
      CALL Fatal(Caller,'Impossible number of dofs! (should be 1)')
    END IF
    PresVar => acc_solver % Variable
    Pressure => PresVar % Values
    PressurePerm => PresVar % Perm
    Varname = TRIM(PresVar % Name)
    IF( COUNT( PressurePerm > 0 ) <= 0) THEN
      proceed = 0
      noiter = 0
      RETURN
    END IF
    ManningModel = GetLogical( Params,'Manning Model',GotIt)
    IF( ManningModel ) THEN
      GravityCoeff = GetCReal( CurrentModel % Constants,'Gravity Coefficient',GotIt)
      IF(.NOT. GotIt) GravityCoeff = 9.81_dp
    END IF
    AnyBC = ListGetLogicalAnyBC( acc_model, 'Open Side') .OR. &
        ListCheckPresentAnyBC( acc_model, 'Filmpressure Flux') .OR. &
        ListCheckPresentAnyBC( acc_model, 'Filmpressure Velocity') .OR. &
        ListCheckPresentAnyBC( acc_model, 'Filmpressure Transfer Coefficient')
    MinGap = ListGetCReal( Params,'Min Gap Height',GotMinGap)
    MinGradPres = ListGetCReal( Params,'Initial Pressure Gradient',GotIt)
    IF(.NOT. GotIt .OR. AllocationsDone ) THEN
      MinGradPres = EPSILON( MinGradPres )
    END IF
    NoIterations = GetInteger( Params,'Nonlinear System Max Iterations',GotIt)
    IF(.NOT. GotIt) NoIterations = 1
    noiter = NoIterations
    proceed = 1
  end subroutine rrs_preamble

  ! Lines 157-203: first-mesh permanent storage + Skip First Solution.
  subroutine rrs_allocate(skipfirst) bind(c, name="rrs_allocate")
    implicit none
    integer(c_int) :: skipfirst
    skipfirst = 0
    IF ( .NOT. AllocationsDone  ) THEN
      n = acc_solver % Mesh % MaxElementNodes
      ALLOCATE(ElementNodes % x(n),  &
          ElementNodes % y(n),       &
          ElementNodes % z(n),       &
          Viscosity(n),              &
          GapHeight(n),          &
          ExtPres(n), &
          FluxPres(n), &
          VeloPres(n), &
          CoeffPres(n), &
          ElemArtif(n), &
          ElemDensity(n), &
          Velocity(3,N),         &
          NormalVelocity(n),     &
          Admittance(n),         &
          FORCE( 2*N ),           &
          STIFF( 2*N, 2*N ), &
          MASS( 2*N, 2*N ), &
          TimeForce( 2*N ), &
          ElemPressure(n), &
          PrevElemPressure(n), &
          ExtPressure(n), &
          ElemPseudoPressure(n), &
          STAT=istat )
      IF ( istat /= 0 ) CALL FATAL(Caller,'Memory allocation error')
      GotPseudoPressure = .FALSE.
      DO k=1,acc_model % NumberOfMaterials
        Material => acc_model % Materials(k) % Values
        CompressibilityModel = ListGetString( Material,'Compressibility Model', GotIt)
        IF (.NOT. GotIt ) CYCLE
        IF( CompressibilityModel == 'artificial compressible') THEN
          GotPseudoPressure = .TRUE.
          ALLOCATE( PseudoPressure(SIZE(Pressure)),STAT=istat)
          EXIT
        END IF
      END DO
      AllocationsDone = .TRUE.
      IF(ListGetLogical( Params,'Skip First Solution', GotIt ) ) THEN
        CALL Info(Caller,'Skipping first solution phase completely!',Level=5)
        skipfirst = 1
        RETURN
      END IF
    END IF
  end subroutine rrs_allocate

  ! Lines 205-230: per-frame setup (pseudo-pressure, FSI multiplier).
  subroutine rrs_frame_setup() bind(c, name="rrs_frame_setup")
    implicit none
    IF(GotPseudoPressure) THEN
      PseudoPressure = Pressure
      ACScale = ListGetConstReal( acc_model % Simulation, &
          'Artificial Compressibility Scaling',GotIt)
      IF(.NOT.GotIt) ACScale = 1.0
    END IF
    FsiMult = ListGetCReal( Params,'Fsi Velocity Multiplier',GotIt )
    IF(.NOT. GotIt) FsiMult = 1.0_dp
    itime = GetTimestep()
    IF(itime == 1) THEN
      IF( ListGetLogical( Params,'Steady State Fsi Start', GotIt )  ) THEN
        CALL Info(Caller,'Ignoring FSI velocity for the 1st timestep',Level=7)
        FsiMult = 0.0_dp
      END IF
    END IF
    mat_idold = 0
    CALL Info(Caller,'-------------------------------------------------',Level=5)
  end subroutine rrs_frame_setup

  ! Lines 234-238: iteration header (LinearModel feeds the assembly).
  subroutine rrs_iter_head(iter_in) bind(c, name="rrs_iter_head")
    implicit none
    integer(c_int), value :: iter_in
    iter = iter_in
    LinearModel = ( iter == 1 ) .AND. ListGetLogical( Params, 'Linear First Iteration',GotIt)
    WRITE(Message,'(A,T35,I5)') 'Reynolds iteration:',iter
    CALL Info(Caller,Message,Level=5)
  end subroutine rrs_iter_head

  subroutine rrs_initialize() bind(c, name="rrs_initialize")
    implicit none
    CALL DefaultInitialize()
  end subroutine rrs_initialize

  subroutine rrs_bulk(sens) bind(c, name="rrs_bulk")
    implicit none
    integer(c_int), value :: sens
    CALL GlobalBulkAssembly( sens )
    CALL DefaultFinishBulkAssembly( )
  end subroutine rrs_bulk

  subroutine rrs_boundary() bind(c, name="rrs_boundary")
    implicit none
    IF(AnyBC) THEN
      CALL GlobalBoundaryAssemby()
    END IF
    CALL DefaultFinishBoundaryAssembly( )
  end subroutine rrs_boundary

  subroutine rrs_finish_dirichlet() bind(c, name="rrs_finish_dirichlet")
    implicit none
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()
  end subroutine rrs_finish_dirichlet

  subroutine rrs_linesearch(redo, conv) bind(c, name="rrs_linesearch")
    implicit none
    integer(c_int) :: redo, conv
    IF( DefaultLinesearch( Converged ) ) THEN
      redo = 1
    ELSE
      redo = 0
    END IF
    conv = MERGE(1, 0, Converged)
  end subroutine rrs_linesearch

  subroutine rrs_solve() bind(c, name="rrs_solve")
    implicit none
    Norm = DefaultSolve()
  end subroutine rrs_solve

  subroutine rrs_nonlin(nc) bind(c, name="rrs_nonlin")
    implicit none
    integer(c_int) :: nc
    nc = acc_solver % Variable % NonlinConverged
  end subroutine rrs_nonlin

  subroutine rrs_sensitivity(redo) bind(c, name="rrs_sensitivity")
    implicit none
    integer(c_int) :: redo
    IF(DefaultSensitivity()) THEN
      NoIterations = 1
      redo = 1
    ELSE
      redo = 0
    END IF
  end subroutine rrs_sensitivity

  subroutine rrs_gap_gate(doit) bind(c, name="rrs_gap_gate")
    implicit none
    integer(c_int) :: doit
    IF( ListGetLogical( Params,'Gap Sensitivity', GotIt ) ) THEN
      doit = 1
    ELSE
      doit = 0
    END IF
  end subroutine rrs_gap_gate

  ! Lines 278-308: gap-sensitivity block (variable swap + re-assembly).
  subroutine rrs_gap_sensitivity() bind(c, name="rrs_gap_sensitivity")
    implicit none
    CALL Info(Caller,'Computing FilmPressure sentivity to gap height',Level=5)
    SensVar => VariableGet( acc_model % Variables,TRIM(Varname)//' Gap Sensitivity')
    IF( .NOT. ASSOCIATED( SensVar ) ) THEN
      CALL Fatal(Caller,'> '//TRIM(Varname)//' gap sensitivity < should exist!')
    END IF
    CALL ListAddLogical(Params,'Skip Compute Nonlinear Change',.TRUE.)
    ApplyLimiter = ListGetLogical( Params,'Apply Limiter', GotIt )
    IF( ApplyLimiter ) CALL ListAddLogical( Params,'Apply Limiter', .FALSE. )
    acc_solver % Variable => SensVar
    SensVar % Values = 0.0_dp
    CALL DefaultInitialize()
    CALL GlobalBulkAssembly( 1 )
    CALL DefaultFinishBulkAssembly( )
    CALL DefaultFinishBoundaryAssembly( )
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs( Ux = SensVar )
    Norm = DefaultSolve()
    CALL ListAddLogical(Params,'Skip Compute Nonlinear Change',.FALSE.)
    acc_solver % Variable => PresVar
    IF( ApplyLimiter ) CALL ListAddLogical( Params,'Apply Limiter', .TRUE. )
  end subroutine rrs_gap_sensitivity

  subroutine rrs_finish() bind(c, name="rrs_finish")
    implicit none
    CALL DefaultFinish()
    CALL Info(Caller,'-------------------------------------------------',Level=5)
  end subroutine rrs_finish

  ! Verbatim ReynoldsSolver_init body (pure list manipulation).
  subroutine rrs_init_wrap() bind(c, name="rrs_init_wrap")
    implicit none
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL :: Found
    CHARACTER(LEN=MAX_NAME_LEN) :: Varname
    Params => GetSolverParams()
    VarName = ListGetString( Params,'Variable',Found)
    IF(.NOT. Found ) THEN
      Varname = 'FilmPressure'
      CALL ListAddString( Params, 'Variable', VarName )
    END IF
    CALL ListAddLogical( Params,'Limiter Load Sign Negative',.TRUE.)
    IF( ListGetLogical( Params,'Gap Sensitivity', Found ) ) THEN
      CALL ListAddStrinG( Params,NextFreeKeyword('Exported Variable',Params),&
          TRIM(VarName)//' Gap Sensitivity')
    END IF
    IF( ListGetLogical( Params,'Calculate Sensitivity', Found ) ) THEN
      CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
          'Filmpressure Sensitivity' )
      CALL ListAddString(Params,'Sensitivity Variable','Filmpressure Sensitivity')
    END IF
  end subroutine rrs_init_wrap

  ! --- transcribed CONTAINS kernels (statement-for-statement) ---

  SUBROUTINE GlobalBulkAssembly( SensitivityMode )
    implicit none
    INTEGER, OPTIONAL :: SensitivityMode
    INTEGER :: SensMode
    SensMode = 0
    IF( PRESENT(SensitivityMode) ) SensMode = SensitivityMode
    mat_idold = -1
    DO t=1,acc_solver % NumberOfActiveElements
      Element => GetActiveElement(t)
      IF( Element % TYPE % ElementCode > 500 ) THEN
        CALL Fatal(Caller,'This is a reduced dimensional solver for 1D and 2D only!')
      END IF
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      CALL GetElementNodes( ElementNodes )
      CALL GetScalarLocalSolution( ElemPressure, UVariable = PresVar)
      IF( SensMode > 0 ) THEN
        IF( transient_mod ) THEN
          CALL GetScalarLocalSolution( PrevElemPressure, tstep = -1, UVariable = PresVar )
        END IF
      END IF
      IF( GotPseudoPressure ) THEN
        ElemPseudoPressure(1:n) = PseudoPressure( PressurePerm(Element % NodeIndexes))
      END IF
      body_id =  Element % Bodyid
      ent_id = ListGetInteger( acc_model % Bodies(body_id) % Values,'Equation')
      Equation => acc_model % Equations(ent_id) % Values
      mat_id = GetInteger( acc_model % Bodies( body_id ) % Values, 'Material')
      Material => acc_model % Materials(mat_id) % Values
      ent_id = GetInteger( acc_model % Bodies( ent_id ) % Values, 'Body Force',GotIt)
      IF(ent_id>0) THEN
        BodyForce => acc_model % BodyForces(ent_id) % Values
      ELSE
        BodyForce => NULL()
      END IF
      Velocity = 0.0_dp
      UseVelocity = .FALSE.
      GotIt = .FALSE.; GotIt2 = .FALSE.; GotIt3 = .FALSE.
      IF( ListCheckPrefix( BodyForce,'Surface Velocity') ) THEN
        Velocity(1,1:n) = GetReal(BodyForce,'Surface Velocity 1',GotIt)
        Velocity(2,1:n) = GetReal(BodyForce,'Surface Velocity 2',GotIt2)
        Velocity(3,1:n) = GetReal(BodyForce,'Surface Velocity 3',GotIt3)
        UseVelocity = GotIt .OR. GotIt2 .OR. GotIt3
      END IF
      IF(.NOT. UseVelocity) THEN
        IF( ListCheckPrefix( Material,'Surface Velocity') ) THEN
          Velocity(1,1:n) = GetReal(Material,'Surface Velocity 1',GotIt)
          Velocity(2,1:n) = GetReal(Material,'Surface Velocity 2',GotIt2)
          Velocity(3,1:n) = GetReal(Material,'Surface Velocity 3',GotIt3)
          UseVelocity = GotIt .OR. GotIt2 .OR. GotIt3
        END IF
      END IF
      IF(.NOT. UseVelocity) THEN
        IF( ListCheckPrefix( BodyForce,'Tangent Velocity') ) THEN
          Velocity(1,1:n) = GetReal(BodyForce,'Tangent Velocity 1',GotIt)
          Velocity(2,1:n) = GetReal(BodyForce,'Tangent Velocity 2',GotIt2)
          Velocity(3,1:n) = GetReal(BodyForce,'Tangent Velocity 3',GotIt3)
        END IF
        IF(.NOT. (GotIt .OR. GotIt2 .OR. GotIt3)) THEN
          IF( ListCheckPrefix( Material,'Tangent Velocity') ) THEN
            Velocity(1,1:n) = GetReal(Material,'Tangent Velocity 1',GotIt)
            Velocity(2,1:n) = GetReal(Material,'Tangent Velocity 2',GotIt2)
            Velocity(3,1:n) = GetReal(Material,'Tangent Velocity 3',GotIt3)
          END IF
        END IF
      END IF
      NormalVelocity(1:n) = GetReal(BodyForce,'Normal Velocity',GotIt)
      IF(.NOT. GotIt) NormalVelocity(1:n) = GetReal(Material,'Normal Velocity',GotIt)
      NormalVelocity(1:n) = NormalVelocity(1:n) + &
          FsiMult * GetReal(BodyForce,'Fsi Velocity',GotIt)
      IF(.NOT. GotIt) NormalVelocity(1:n) = NormalVelocity(1:n) + &
          FsiMult * GetReal(Material,'Fsi Velocity',GotIt)
      IF( ManningModel ) THEN
        ElemDensity(1:) = GetReal( Material,'Density')
      END IF
      GapHeight(1:n) = GetReal( Material,'Gap Height')
      IF(GotMinGap) GapHeight(1:n) = MAX(GapHeight(1:n),MinGap)
      Admittance(1:n) = GetReal( Material, 'Flow Admittance', GotIt)
      IF(ASSOCIATED(BodyForce)) THEN
        ExtPressure(1:n) = GetReal( BodyForce, 'External FilmPressure', GotIt)
      ELSE
        ExtPressure(1:n) = 0.0_dp
      END IF
      Viscosity(1:n) = GetReal( Material, 'Viscosity')
      IF(mat_id /= mat_idold) THEN
        mat_idold = mat_id
        IF( ManningModel ) THEN
          ManningCoeff = GetCReal(Material,'Manning Coefficient')
        END IF
        ReferencePressure = GetCReal( Material,'Reference Pressure', GotIt )
        ViscosityModel = GetString(Material,'Viscosity Model',GotIt)
        IF(GotIt) THEN
          IF( ViscosityModel == 'newtonian') THEN
            ViscosityType = Viscosity_Newtonian
          ELSE IF( ViscosityModel == 'rarefied') THEN
            ViscosityType = Viscosity_Rarefied
            mfp0 = GetCReal(Material,'Mean Free Path')
          ELSE
            CALL Warn(Caller,'Unknown viscosity model')
          END IF
        ELSE
          ViscosityType = Viscosity_Newtonian
        END IF
        CompressibilityType = Compressibility_None
        IF( .NOT. LinearModel ) THEN
          CompressibilityModel = GetString(Material,'Compressibility Model',GotIt)
          IF(GotIt) THEN
            IF(CompressibilityModel == 'incompressible') THEN
              CompressibilityType = Compressibility_None
            ELSE IF(CompressibilityModel == 'weakly compressible') THEN
              CompressibilityType = Compressibility_Weak
              BulkModulus = GetCReal( Material, 'Bulk Modulus')
            ELSE IF(CompressibilityModel == 'isothermal ideal gas') THEN
              CompressibilityType = Compressibility_GasIsothermal
            ELSE IF(CompressibilityModel == 'adiabatic ideal gas') THEN
              CompressibilityType = Compressibility_GasAdiabatic
              HeatRatio = GetCReal( Material, 'Specific Heat Ratio')
            ELSE IF( CompressibilityModel == 'artificial compressible') THEN
              CompressibilityType = Compressibility_Artificial
            ELSE
              CALL Fatal(Caller,'Unknown compressibility model: '//TRIM(CompressibilityModel))
            END IF
          END IF
        END IF
      END IF
      IF( CompressibilityType == Compressibility_Artificial ) THEN
        ElemArtif(1:n) = GetReal( Material,'Artificial Compressibility',GotIt)
      END IF
      STIFF = 0.0_dp
      MASS = 0.0_dp
      FORCE = 0.0_dp
      CALL LocalBulkMatrix( MASS, STIFF, FORCE, Element, n, nd, ElementNodes, SensMode )
      TimeForce  = 0.0_dp
      IF ( transient_mod ) THEN
        CALL Default1stOrderTime( MASS,STIFF,FORCE)
      END IF
      IF (  Bubbles ) THEN
        CALL Condensate( N, STIFF, FORCE, TimeForce )
      END IF
      CALL DefaultUpdateEquations( STIFF, FORCE )
    END DO
  END SUBROUTINE GlobalBulkAssembly

  SUBROUTINE LocalBulkMatrix(MassMatrix, StiffMatrix, ForceVector, &
      Element, n, nd, Nodes, SensMode )
    implicit none
    REAL(KIND=dp) :: MassMatrix(:,:), StiffMatrix(:,:), ForceVector(:)
    INTEGER :: n, nd
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), TARGET :: Element
    INTEGER :: SensMode
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3), detJ
    REAL(KIND=dp) :: x,y,z,Metric(3,3),SqrtMetric,Symb(3,3,3),dSymb(3,3,3,3)
    REAL(KIND=dp) :: U, V, W, S, MS, MM, MA, L, A, B, HR, SL(3), SLR, SLL(3), F, AS
    REAL(KIND=dp) :: Normal(3), Velo(3), NormalVelo, TangentVelo(3), Damp, Pres, PrevPres, &
        TotPres, GradPres(3), AbsGradPres, dPdt, Gap, Visc, mfp, Kn, Density, DensityDer, &
        PseudoPres, ExtPres
    LOGICAL :: Stat, GotAC
    INTEGER :: i,p,q,t,DIM, NBasis, CoordSys
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
    DIM = CoordinateSystemDimension()
    CoordSys = CurrentCoordinateSystem()
    GotAC = .FALSE.
    Metric = 0.0_dp
    Metric(1,1) = 1.0_dp
    Metric(2,2) = 1.0_dp
    Metric(3,3) = 1.0_dp
    NBasis = n
    IF ( Bubbles ) THEN
      NBasis = 2*n
      IntegStuff = GaussPoints( Element, Element % TYPE % GaussPoints2 )
    ELSE
      NBasis = nd
      IntegStuff = GaussPoints( Element )
    END IF
    DO t=1,IntegStuff % n
      u = IntegStuff % u(t)
      v = IntegStuff % v(t)
      w = IntegStuff % w(t)
      s = IntegStuff % s(t)
      stat = ElementInfo( Element, Nodes, u, v, w, DetJ, &
               Basis, dBasisdx, Bubbles = Bubbles)
      s = s * DetJ
      IF ( CoordSys /= Cartesian ) THEN
        x = SUM( Nodes % x(1:n) * Basis(1:n) )
        y = SUM( Nodes % y(1:n) * Basis(1:n) )
        z = SUM( Nodes % z(1:n) * Basis(1:n) )
        CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,X,Y,Z )
        s = s * SqrtMetric
      END IF
      NormalVelo = SUM( Basis(1:n) * NormalVelocity(1:n) )
      IF(UseVelocity) THEN
        IF( ASSOCIATED( Element % BoundaryInfo)) THEN
          Normal = NormalVector( Element,Nodes,u,v,.TRUE. )
        ELSE
          Normal = NormalVector( Element,Nodes,u,v,.FALSE. )
        END IF
        DO i=1,3
          Velo(i) = SUM( Basis(1:n) * Velocity(i,1:n) )
        END DO
        NormalVelo = NormalVelo + SUM( Normal * Velo)
        TangentVelo = Velo - NormalVelo * Normal
      ELSE
        DO i=1,3
          TangentVelo(i) = SUM( Basis(1:n) * Velocity(i,1:n) )
        END DO
      END IF
      Damp = SUM(Basis(1:n) * Admittance(1:n))
      ExtPres = SUM(Basis(1:n) * ExtPressure(1:n))
      Pres = SUM(Basis(1:n) * ElemPressure(1:n))
      Gap = SUM(Basis(1:n) * GapHeight(1:n))
      TotPres = ReferencePressure + Pres
      IF( SensMode > 0 ) THEN
        IF( transient_mod ) THEN
          PrevPres = SUM( Basis(1:n) * PrevElemPressure(1:n) )
          dPdt = ( Pres - PrevPres ) / dt_mod
        ELSE
          dPdt = 0.0_dp
        END IF
        DO i = 1,3
          GradPres(i) = SUM( dBasisdx(1:n,i) * ElemPressure(1:n) )
        END DO
      END IF
      SELECT CASE (ViscosityType)
      CASE (Viscosity_Newtonian)
        Visc = SUM(Basis(1:n) * Viscosity(1:n))
      CASE (Viscosity_Rarefied)
        Visc = SUM(Basis(1:n) * Viscosity(1:n))
        mfp = mfp0 * ReferencePressure / TotPres
        Kn = mfp / ABS(Gap)
        Visc = Visc / (1+9.638_dp*Kn**1.159_dp)
      END SELECT
      SELECT CASE (CompressibilityType)
      CASE (Compressibility_None)
        Density = 1.0d0
        DensityDer = 0.0d0
      CASE (Compressibility_Weak)
        Density = EXP(TotPres/BulkModulus)
        DensityDer = Density / BulkModulus
      CASE (Compressibility_GasIsothermal)
        Density = TotPres
        DensityDer = 1.0d0
      CASE(Compressibility_GasAdiabatic)
        Density = TotPres ** (1.0_dp/HeatRatio)
        DensityDer = (1/HeatRatio) * TotPres ** (1.0_dp/HeatRatio - 1.0_dp)
      CASE (Compressibility_Artificial )
        Density = 1.0d0
        DensityDer = 0.0d0
        GotAC = .TRUE.
      END SELECT
      IF( ManningModel ) THEN
        Density = Density * SUM( Basis(1:n) * ElemDensity(1:n) )
        DensityDer = DensityDer * SUM( Basis(1:n) * ElemDensity(1:n) )
        DO i = 1,3
          GradPres(i) = SUM( dBasisdx(1:n,i) * ElemPressure(1:n) )
        END DO
        AbsGradPres = SQRT( SUM( GradPres**2 ) )
        AbsGradPres = MAX( AbsGradPres, MinGradPres )
        MS = -SQRT(Density/(GravityCoeff*AbsGradPres)) * Gap**(5.0/3)  / (2**(2.0/3) * ManningCoeff)
      ELSE
        MS = -Density * Gap**3 / (12 * Visc)
      END IF
      IF(CompressibilityType == Compressibility_GasIsothermal) THEN
        HR = -Damp * ( Density + ExtPres ) / 2
      ELSE IF(CompressibilityType == Compressibility_GasAdiabatic) THEN
        HR = -Damp * ( Density  + ExtPres ** (1.0_dp/HeatRatio) ) / 2
      ELSE
        HR = -Damp * Density
      END IF
      MM = -DensityDer * Gap
      MA = 0.0_dp
      IF( GotAC ) THEN
        PseudoPres = SUM(Basis(1:n) * ElemPseudoPressure(1:n) )
        MA = ( -Density / dt_mod ) * SUM( ElemArtif(1:n) * Basis(1:n) )
        MA = ACScale * MA
      END IF
      L = Density * NormalVelo
      SLR = 0.0_dp
      SLL = 0.0_dp
      DO i=1,dim
        SLR = SLR + 0.5_dp * Density * SUM( dBasisdx(1:n,i) * Velocity(i,1:n) * GapHeight(1:n))
      END DO
      SLL = -0.5_dp* DensityDer * Gap * TangentVelo
      DO p=1,NBasis
        F = 0.0_dp
        DO q=1,NBasis
          A = (MA + HR) * Basis(q) * Basis(p)
          DO i=1,DIM
            DO j=1,DIM
              A = A + MS * Metric(i,j) * dBasisdx(q,i) * dBasisdx(p,j)
              A = A + SLL(j) * Metric(i,j) * dBasisdx(q,i) * Basis(p)
            END DO
          END DO
          IF( SensMode == 1 ) THEN
            A = A + Damp * ExtPres * ( ExtPres/TotPres - 1.0_dp) / 2
          END IF
          StiffMatrix(p,q) = StiffMatrix(p,q) + s * A
          IF( transient_mod ) THEN
            B = MM * Basis(q) * Basis(p)
            MassMatrix(p,q)  = MassMatrix(p,q)  + s * B
          END IF
        END DO
        IF( SensMode == 0 ) THEN
          F = L + SLR + HR * ExtPres
          IF(GotAC) F = F + MA * PseudoPres
        ELSE IF( SensMode == 1 ) THEN
          IF( transient_mod ) THEN
            F = -2.0 * DensityDer * dPdt
          END IF
          F = F - DensityDer * SUM( TangentVelo * GradPres )
          DO i = 1,dim
            F = F - 1.5_dp * ( Density / Gap ) * SUM( dBasisdx(1:n,i) * Velocity(i,1:n) * GapHeight(1:n) )
          END DO
          F = F - Density * NormalVelo * (3 / Gap )
          F = F + HR * ( ExtPres - Pres) * ( 3 / Gap )
        END IF
        ForceVector(p) = ForceVector(p) + s * Basis(p) * F
      END DO
    END DO
  END SUBROUTINE LocalBulkMatrix

   SUBROUTINE GlobalBoundaryAssemby()
    implicit none
    DO t=1, acc_solver % Mesh % NumberOfBoundaryElements
      Element => GetBoundaryElement(t)
      IF ( .NOT. ActiveBoundaryElement() ) CYCLE
      BC => GetBC()
      IF ( .NOT. ASSOCIATED( BC ) ) CYCLE
      n  = GetElementNOFNodes()
      OpenSide = GetLogical(BC,'Open Side',gotIt)
      FluxPres(1:n) = GetReal(BC,'Filmpressure Flux',GotFlux)
      VeloPres(1:n) = GetReal(BC,'Filmpressure Velocity',GotVelo)
      CoeffPres(1:n) = GetReal(BC,'Filmpressure Transfer Coefficient',GotIt)
      ExtPres(1:n) = GetReal(BC,'External FilmPressure',GotExt)
      IF(GotExt .NEQV. GotIt) THEN
        CALL Fatal(Caller,'Give neither or both keywords for Robin BC!')
      END IF
      IF(.NOT. (OpenSide .OR. GotExt .OR. GotFlux .OR. GotVelo ) ) CYCLE
      nd = GetElementNOFDOFs()
      IF ( GetElementFamily() == 1 ) CYCLE
      NodeIndexes => Element % NodeIndexes
      IF ( ANY( PressurePerm(NodeIndexes(1:n)) == 0 ) ) CYCLE
      Parent => Element % BoundaryInfo % Left
      stat = ASSOCIATED( Parent )
      IF ( stat ) stat = ALL(PressurePerm(Parent % NodeIndexes) > 0)
      IF(.NOT. stat) THEN
        Parent => ELement % BoundaryInfo % Right
        stat = ASSOCIATED( Parent )
        IF ( stat ) stat = ALL(PressurePerm(Parent % NodeIndexes) > 0)
        IF ( .NOT. stat )  CALL Fatal( Caller, &
            'No proper parent element available for specified boundary' )
      END IF
      CALL GetElementNodes( ElementNodes )
      mat_id = GetInteger( acc_model % Bodies(Parent % BodyId) % Values,'Material')
      Material => acc_model % Materials(mat_id) % Values
      GapHeight(1:n) = GetReal(Material,'Gap Height')
      IF(GotMinGap) GapHeight(1:n) = MAX( GapHeight(1:n), MinGap )
      Viscosity(1:n) = GetReal( Material, 'Viscosity')
      IF( ManningModel ) THEN
        ElemDensity(1:) = GetReal( Material,'Density')
      END IF
      STIFF = 0.0d0
      MASS = 0.0d0
      FORCE = 0.0d0
      CALL LocalBoundaryMatrix( MASS, STIFF, FORCE, Element, n, ElementNodes )
      IF ( transient_mod ) THEN
        MASS = 0.d0
        CALL Default1stOrderTime( MASS, STIFF, FORCE )
      END IF
      CALL DefaultUpdateEquations( STIFF, FORCE )
    END DO
  END SUBROUTINE GlobalBoundaryAssemby

  SUBROUTINE LocalBoundaryMatrix(MASS, STIFF, FORCE, Element, n, Nodes)
    implicit none
    REAL(KIND=dp) :: MASS(:,:), STIFF(:,:), FORCE(:)
    INTEGER :: n
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), TARGET :: Element
    REAL(KIND=dp) :: DetJ,U,V,W,S
    REAL(KIND=dp) :: Basis(n)
    REAL(KIND=dp) :: Visc, dl, mfp, Kn, Damp, TotPres, Pres, Density, Gap, A, B
    LOGICAL :: Stat
    INTEGER :: i,p,q,t,DIM,CoordSys
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
    DIM = CoordinateSystemDimension()
    CoordSys = CurrentCoordinateSystem()
    IntegStuff = GaussPoints( Element )
    DO t=1,IntegStuff % n
      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)
      S = IntegStuff % s(t)
      stat = ElementInfo( Element, Nodes, U, V, W, DetJ, Basis )
      s = s * DetJ
      Gap = SUM( GapHeight(1:n) * Basis(1:n) )
      Pres = SUM( ElemPressure(1:n) * Basis(1:n))
      TotPres = ReferencePressure + Pres
      SELECT CASE (ViscosityType)
      CASE (Viscosity_Newtonian)
        Visc = SUM(Basis(1:n) * Viscosity(1:n))
        dl = 0.8488_dp
      CASE (Viscosity_Rarefied)
        Visc = SUM(Basis(1:n) * Viscosity(1:n))
        mfp = mfp0 * ReferencePressure / TotPres
        Kn = mfp / ABS(Gap)
        Visc = Visc / (1+9.638_dp*Kn**1.159_dp)
        dl = 0.8488_dp*(1+2.676_dp*Kn**0.659_dp)
      END SELECT
      Damp = Gap / (12 * dl * Visc)
      SELECT CASE (CompressibilityType)
      CASE (Compressibility_None)
        Density = 1.0d0
      CASE (Compressibility_Weak)
        Density = EXP(TotPres/BulkModulus)
      CASE (Compressibility_GasIsothermal)
        Density = TotPres
      CASE(Compressibility_GasAdiabatic)
        Density = TotPres ** (1.0/HeatRatio)
      END SELECT
      IF(ManningModel) THEN
        Density = Density * SUM( Basis(1:n) * ElemDensity(1:n) )
      END IF
      A = 0.0_dp
      B = 0.0_dp
      IF( GotExt ) THEN
        A = -SUM(Basis(1:n) * CoeffPres(1:n) )
        B = A * SUM( Basis(1:n) * ExtPres(1:n) )
      END IF
      IF( GotFlux ) THEN
        B = B - SUM( Basis(1:n) * FluxPres(1:n) )
      END IF
      IF( GotVelo ) THEN
        B = B - Gap * SUM( Basis(1:n) * VeloPres(1:n) )
      END IF
      IF(OpenSide) THEN
        A = A - Damp * Gap * Density
      END IF
      DO p=1,n
        FORCE(p) = FORCE(p) + s * Basis(p) * B
        DO q=1,n
          STIFF(p,q) = STIFF(p,q) + s * Basis(q) * Basis(p) * A
        END DO
      END DO
     END DO
  END SUBROUTINE LocalBoundaryMatrix

end module reynolds_solve_acc

! ---------------------------------------------------------------------------
! Mirror of ReynoldsPostprocess: same wrapper discipline as above. The
! per-component element loop and the LocalMatrix kernel are transcribed
! statement-for-statement (all SUM reductions, the strided
! Values(Component::Components) store, and the formatted force / moment /
! heating tails stay in Fortran).
! ---------------------------------------------------------------------------
module reynolds_post_acc
  use iso_c_binding
  use DefUtils
  use Types
  use Lists
  use Integration
  use ElementDescription
  use SolverBasics
  implicit none

  ! -- mirrored SAVE state of ReynoldsPostprocess (identical lifetimes) --
  type(Variable_t), pointer :: PressureVar => null(), VarResult => null(), &
      SolverVar => null(), SensVar => null()
  type(Nodes_t) :: ElementNodes
  type(Element_t), pointer :: Element => null()
  type(ValueList_t), pointer :: Params => null(), Material => null(), &
      Equation => null(), BodyForce => null()
  type(Model_t), pointer :: acc_model => null()
  type(Solver_t), pointer :: acc_solver => null()
  integer :: iter, i, j, k, l, n, nd, t, istat, ent_id, body_id, mat_id, &
      mat_idold, ViscosityType, dim
  integer :: Component, Components, Mode
  integer, pointer :: NodeIndexes(:), PressurePerm(:)
  real(kind=dp), pointer :: mWork(:,:)
  logical :: GotIt, GotIt2, GotIt3, stat, UseVelocity, &
      AllocationsDone = .FALSE., OpposingWall, CalculateMoment, &
      ManningModel, GotMinGap
  real(kind=dp), pointer :: Pressure(:), gWork(:,:)
  real(kind=dp) :: Norm, ReferencePressure, mfp0, HeatSlide, HeatPres, &
      HeatTotal, Pforce(3), Vforce(3), Sforce(3), TotForce, Moment(3), &
      MomentAbout(3), AmbientPres, ManningCoeff, GravityCoeff, MinGap
  real(kind=dp), allocatable :: STIFF(:,:), FORCE(:), Viscosity(:), &
      GapHeight(:), Velocity(:,:), ElemPressure(:), SensPressure(:), &
      BotHeight(:), ElemDensity(:)
  character(len=MAX_NAME_LEN) :: ViscosityModel, PressureName
  character(*), parameter :: Caller = 'ReynoldsPostprocess'
  integer, parameter :: Viscosity_Newtonian = 1, Viscosity_Rarefied = 2
  real(kind=dp) :: dt_mod
  logical :: transient_mod

contains

  subroutine rpp_frame_begin(h_model, h_solver, dt, transient) &
      bind(c, name="rpp_frame_begin")
    implicit none
    type(c_ptr), value :: h_model, h_solver
    real(c_double), value :: dt
    integer(c_int), value :: transient
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    dt_mod = dt
    transient_mod = (transient /= 0)
  end subroutine rpp_frame_begin

  subroutine rpp_frame_end() bind(c, name="rpp_frame_end")
    implicit none
    acc_model => null()
    acc_solver => null()
    Element => null()
    Params => null()
    Material => null()
    Equation => null()
    BodyForce => null()
    PressureVar => null()
    VarResult => null()
    SolverVar => null()
    SensVar => null()
    Pressure => null()
    PressurePerm => null()
    NodeIndexes => null()
    mWork => null()
    gWork => null()
  end subroutine rpp_frame_end

  ! Original lines 1085-1191: checks, parameter fetch, first-mesh
  ! storage, moment setup; proceed=0 reproduces the bare RETURNs.
  subroutine rpp_preamble(proceed) bind(c, name="rpp_preamble")
    implicit none
    integer(c_int) :: proceed
    proceed = 0
    CALL Info(Caller,'--------------------------------------------------',Level=5)
    CALL Info(Caller,'Computing postprocessing fields from film pressure',Level=5)
    CALL Info(Caller,'--------------------------------------------------',Level=5)
    Params => GetSolverParams()
    IF ( .NOT. ASSOCIATED( acc_solver % Matrix ) ) RETURN
    IF(acc_solver % Variable % Dofs /= 1) THEN
      CALL Fatal(Caller,'Impossible number of dofs! (should be 1)')
    END IF
    ManningModel = GetLogical(Params,'Manning Model',GotIt)
    PressureName = GetString(Params,'Reynolds Pressure Variable Name',GotIt)
    IF(.NOT. GotIt) PressureName = 'FilmPressure'
    PressureVar => VariableGet( acc_solver % Mesh % Variables, PressureName)
    IF(.NOT. ASSOCIATED(PressureVar)) THEN
      CALL Info(Caller,'Give pressure variable name with: "Reynolds Pressure Variable Name"',Level=3)
      CALL Fatal(Caller,'Could not find primary variable: '//TRIM(PressureName))
    END IF
    SensVar => VariableGet( acc_solver % Mesh % Variables,TRIM(PressureName)//' Gap Sensitivity')
    IF(ASSOCIATED(SensVar)) THEN
      CALL Info(Caller,'We got pressure sensitivity available!')
    END IF
    IF( ManningModel ) THEN
      gWork => ListGetConstRealArray( CurrentModel % Constants,'Gravity',GotIt)
      IF ( GotIt ) THEN
        GravityCoeff = gWork(4,1)
      ELSE
        GravityCoeff = GetCReal( CurrentModel % Constants,'Gravity Coefficient',GotIt)
      END IF
      IF(.NOT. GotIt) GravityCoeff = 9.81_dp
    END IF
    MinGap = ListGetCReal( Params,'Min Gap Height',GotMinGap)
    AmbientPres = ListGetCReal( Params,'Ambient Pressure',GotIt)
    Pressure => PressureVar % Values
    PressurePerm => PressureVar % Perm
    IF( COUNT( PressurePerm > 0 ) <= 0) RETURN
    DIM = CoordinateSystemDimension()
    ReferencePressure = 0.0_dp
    IF ( .NOT. AllocationsDone  ) THEN
      N = acc_solver % Mesh % MaxElementNodes
      ALLOCATE(ElementNodes % x(n),  &
          ElementNodes % y(n),       &
          ElementNodes % z(n),       &
          Viscosity(n),              &
          GapHeight(n),          &
          ElemDensity(n), &
          BotHeight(n), &
          Velocity(3,n),         &
          FORCE(n),           &
          STIFF(n,n), &
          ElemPressure(n), &
          SensPressure(n), &
          STAT=istat )
      IF ( istat /= 0 ) CALL FATAL(Caller,'Memory allocation error')
      AllocationsDone = .TRUE.
    END IF
    mWork => ListGetConstRealArray( Params,'Moment About',CalculateMoment)
    IF( CalculateMoment ) THEN
      MomentAbout(1:dim) = mWork(1:dim,1)
    ELSE
      MomentAbout(1) = ListGetCReal( Params,'Moment About 1', CalculateMoment )
      MomentAbout(2) = ListGetCReal( Params,'Moment About 2', stat )
      CalculateMoment = stat .OR. CalculateMoment
      MomentAbout(3) = ListGetCReal( Params,'Moment About 3', stat )
      CalculateMoment = stat .OR. CalculateMoment
    END IF
    mat_idold = 0
    HeatSlide = 0.0_dp
    HeatPres = 0.0_dp
    Pforce = 0.0_dp
    Sforce = 0.0_dp
    Vforce = 0.0_dp
    Moment = 0.0_dp
    SolverVar => acc_solver % Variable
    IF( .NOT. ASSOCIATED( SolverVar ) ) THEN
      CALL Fatal(Caller,'Solver Variable not associated')
    END IF
    CALL Info(Caller,'Primary variable name: '//TRIM( acc_solver % variable % Name) )
    proceed = 1
  end subroutine rpp_preamble

  ! Original mode dispatch (lines 1193-1221 head): SELECT CASE over Mode
  ! plus the unassociated-VariableGet CYCLE. active=0 means CYCLE.
  subroutine rpp_mode_open(mode_in, comps, active) bind(c, name="rpp_mode_open")
    implicit none
    integer(c_int), value :: mode_in
    integer(c_int) :: comps, active
    Mode = mode_in
    active = 1
    comps = 0
    SELECT CASE( Mode )
    CASE( 0 )
      IF( .NOT. ManningModel ) THEN
        active = 0
        RETURN
      END IF
      VarResult => VariableGet( acc_solver % Mesh % Variables,TRIM(PressureName)//' Corrected')
    CASE( 1 )
      VarResult => VariableGet( acc_solver % Mesh % Variables,TRIM(PressureName)//' Force')
    CASE( 2 )
      VarResult => VariableGet( acc_solver % Mesh % Variables,TRIM(PressureName)//' Flux')
    CASE( 3 )
      VarResult => VariableGet( acc_solver % Mesh % Variables,TRIM(PressureName)//' Mean Velocity')
    CASE( 4 )
      VarResult => VariableGet( acc_solver % Mesh % Variables,TRIM(PressureName)//' Heating')
    CASE DEFAULT
      CALL Fatal(Caller,'Unknow Mode for operation:'//I2S(Mode))
    END SELECT
    IF(.NOT. ASSOCIATED(VarResult)) THEN
      active = 0
      RETURN
    END IF
    Components = VarResult % Dofs
    comps = Components
  end subroutine rpp_mode_open

  ! Original lines 1221-1344: DefaultInitialize plus the per-element
  ! bulk assembly for one (Mode, Component) pair (verbatim).
  subroutine rpp_assemble(mode_in, comp_in) bind(c, name="rpp_assemble")
    implicit none
    integer(c_int), value :: mode_in, comp_in
    Mode = mode_in
    Component = comp_in
    CALL DefaultInitialize()
    DO t=1,acc_solver % NumberOfActiveElements
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      CALL GetElementNodes( ElementNodes )
      NodeIndexes => Element % NodeIndexes
      ElemPressure(1:n) = Pressure(PressurePerm(NodeIndexes(1:n)))
      IF(ASSOCIATED(SensVar)) THEN
        SensPressure(1:n) = SensVar % Values(SensVar % Perm(NodeIndexes(1:n)))
      END IF
      body_id =  Element % Bodyid
      IF( body_id <= 0 ) CYCLE
      ent_id = ListGetInteger( acc_model % Bodies(body_id) % Values,'Equation')
      Equation => acc_model % Equations(ent_id) % Values
      mat_id = GetInteger( acc_model % Bodies( body_id ) % Values, 'Material')
      Material => acc_model % Materials(mat_id) % Values
      ent_id = GetInteger( acc_model % Bodies( ent_id ) % Values, 'Body Force',GotIt)
      IF(ent_id>0) THEN
        BodyForce => acc_model % BodyForces(ent_id) % Values
      ELSE
        BodyForce => NULL()
      END IF
      Velocity = 0.0_dp
      UseVelocity = .FALSE.
      GotIt = .FALSE.; GotIt2 = .FALSE.; GotIt3 = .FALSE.
      IF( ListCheckPrefix( BodyForce,'Surface Velocity') ) THEN
        Velocity(1,1:n) = GetReal(BodyForce,'Surface Velocity 1',GotIt)
        Velocity(2,1:n) = GetReal(BodyForce,'Surface Velocity 2',GotIt2)
        Velocity(3,1:n) = GetReal(BodyForce,'Surface Velocity 3',GotIt3)
        UseVelocity = GotIt .OR. GotIt2 .OR. GotIt3
      END IF
      IF(.NOT. UseVelocity) THEN
        IF( ListCheckPrefix( Material,'Surface Velocity') ) THEN
          Velocity(1,1:n) = GetReal(Material,'Surface Velocity 1',GotIt)
          Velocity(2,1:n) = GetReal(Material,'Surface Velocity 2',GotIt2)
          Velocity(3,1:n) = GetReal(Material,'Surface Velocity 3',GotIt3)
          UseVelocity = GotIt .OR. GotIt2 .OR. GotIt3
        END IF
      END IF
      IF(.NOT. UseVelocity) THEN
        IF( ListCheckPrefix( BodyForce,'Tangent Velocity') ) THEN
          Velocity(1,1:n) = GetReal(BodyForce,'Tangent Velocity 1',GotIt)
          Velocity(2,1:n) = GetReal(BodyForce,'Tangent Velocity 2',GotIt2)
          Velocity(3,1:n) = GetReal(BodyForce,'Tangent Velocity 3',GotIt3)
        END IF
        IF(.NOT. (GotIt .OR. GotIt2 .OR. GotIt3 )) THEN
          IF( ListCheckPrefix( Material,'Tangent Velocity') ) THEN
            Velocity(1,1:n) = GetReal(Material,'Tangent Velocity 1',GotIt)
            Velocity(2,1:n) = GetReal(Material,'Tangent Velocity 2',GotIt2)
            Velocity(3,1:n) = GetReal(Material,'Tangent Velocity 3',GotIt3)
          END IF
        END IF
      END IF
      GapHeight(1:n) = GetReal( Material,'Gap Height')
      IF(GotMinGap) GapHeight(1:n) = MAX(GapHeight(1:n), MinGap)
      Viscosity(1:n) = GetReal( Material, 'Viscosity')
      IF(mat_id /= mat_idold) THEN
        mat_idold = mat_id
        ReferencePressure = GetCReal( Material,'Reference Pressure', GotIt )
        ViscosityModel = GetString(Material,'Viscosity Model',GotIt)
        IF(GotIt) THEN
          IF( ViscosityModel == 'newtonian') THEN
            ViscosityType = Viscosity_Newtonian
          ELSE IF( ViscosityModel == 'rarefied') THEN
            ViscosityType = Viscosity_Rarefied
            mfp0 = GetCReal(Material,'Mean Free Path')
          ELSE
            CALL Fatal(Caller,'Unknown viscosity model: '//TRIM(ViscosityModel))
          END IF
        ELSE
          ViscosityType = Viscosity_Newtonian
        END IF
        OpposingWall = ListGetLogical( Material,'Opposing Wall',GotIt)
      END IF
      IF( ManningModel ) THEN
        ManningCoeff = ListGetCReal( Material,'Manning Coefficient')
      END IF
      IF( Mode == 0 ) THEN
        ElemDensity(1:n) = GetReal( Material,'Density')
        BotHeight(1:n) = GetReal( Material,'Bedrock Elevation')
        VarResult % Values(VarResult % Perm(NodeIndexes)) = Pressure(PressurePerm(NodeIndexes)) &
            - GravityCoeff * ElemDensity(1:n) * ( BotHeight(1:n) + GapHeight(1:n) )
        CYCLE
      END IF
      STIFF = 0.0d0
      FORCE = 0.0d0
      CALL LocalMatrix( STIFF, FORCE, Element, n, nd, ElementNodes)
      CALL DefaultUpdateEquations( STIFF, FORCE )
    END DO
  end subroutine rpp_assemble

  ! Original lines 1348-1362: assembly finish, Dirichlet BCs, solve and
  ! the per-component store (Mode 0 never reaches here: the caller
  ! replicates the line-1346 CYCLE).
  subroutine rpp_solve_store() bind(c, name="rpp_solve_store")
    implicit none
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()
    CALL Info( Caller, 'Dirichlet conditions done', Level=4 )
    Norm = DefaultSolve()
    IF( Mode /= 4 ) THEN
      VarResult % Values(Component::Components) = acc_solver % Variable % Values
    END IF
  end subroutine rpp_solve_store

  ! Original lines 1365-1420: per-mode result reporting (formatted
  ! WRITE + Info + ListAddConstReal kept verbatim in Fortran).
  subroutine rpp_mode_tail() bind(c, name="rpp_mode_tail")
    implicit none
    IF( Mode == 0) THEN
      CONTINUE
    ELSE IF( Mode == 1 ) THEN
      DO i=1,3
        WRITE(Message,'(A,I1,A,T35,ES15.4)') 'Pressure force ',i,' (N):',Pforce(i)
        CALL Info(Caller,Message,Level=5)
        CALL ListAddConstReal( acc_model % Simulation,'res: Pressure force '&
            //I2S(i),Pforce(i))
      END DO
      DO i=1,3
        WRITE(Message,'(A,I1,A,T35,ES15.4)') 'Sliding force ',i,' (N):',Vforce(i)
        CALL Info(Caller,Message,Level=5)
        CALL ListAddConstReal( acc_model % Simulation,'res: Sliding force '&
            //I2S(i),Vforce(i))
      END DO
      TotForce = SQRT( SUM((Pforce + Vforce)**2) )
      WRITE(Message,'(A,T35,ES15.4)') 'Reynolds force (N): ',TotForce
      CALL Info(Caller,Message,Level=5)
      CALL ListAddConstReal( acc_model % Simulation,'res: Reynolds force',TotForce)
      IF(ASSOCIATED(SensVar)) THEN
        DO i=1,3
          WRITE(Message,'(A,I1,A,T35,ES15.4)') 'Pressure spring ',i,' (N/m):',Sforce(i)
          CALL Info(Caller,Message,Level=5)
          CALL ListAddConstReal( acc_model % Simulation,'res: Pressure spring '&
              //I2S(i),Sforce(i))
        END DO
      END IF
      IF( CalculateMoment ) THEN
        DO i=1,3
          WRITE(Message,'(A,I1,A,T35,ES15.4)') 'Reynolds moment ',i,' (Nm):',Moment(i)
          CALL Info(Caller,Message,Level=5)
          CALL ListAddConstReal( acc_model % Simulation,'res: Reynolds moment '&
              //I2S(i),Moment(i))
        END DO
      END IF
    ELSE IF( Mode == 2 ) THEN
      CONTINUE
    ELSE IF( Mode == 3 ) THEN
      CONTINUE
    ELSE IF( Mode == 4 ) THEN
      HeatTotal = HeatPres + HeatSlide
      WRITE(Message,'(A,T35,ES15.4)') 'Pressure heating (W): ',HeatPres
      CALL Info(Caller,Message,Level=5)
      WRITE(Message,'(A,T35,ES15.4)') 'Sliding heating (W): ',HeatSlide
      CALL Info(Caller,Message,Level=5)
      WRITE(Message,'(A,T35,ES15.4)') 'Reynolds heating (W): ',HeatTotal
      CALL Info(Caller,Message,Level=5)
      CALL ListAddConstReal( acc_model % Simulation,'res: Reynolds heating',HeatTotal)
    END IF
  end subroutine rpp_mode_tail

  subroutine rpp_finish() bind(c, name="rpp_finish")
    implicit none
    CALL Info(Caller,'Finished computing postprocessing fields',Level=8)
    CALL Info(Caller,'--------------------------------------------------',Level=8)
  end subroutine rpp_finish

  ! Verbatim ReynoldsPostprocess_Init body (pure list manipulation).
  subroutine rpp_init_wrap() bind(c, name="rpp_init_wrap")
    implicit none
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL :: Found, Calculate, ManningModel
    INTEGER :: Dim,GivenDim,dofs
    CHARACTER(LEN=MAX_NAME_LEN) :: PressureName
    CHARACTER(*), PARAMETER :: CName = 'ReynoldsPostprocess_init'
    Params => GetSolverParams()
    Dim = CoordinateSystemDimension()
    PressureName = GetString(Params,'Reynolds Pressure Variable Name',Found)
    IF(.NOT. Found) PressureName = 'FilmPressure'
    Calculate = ListGetLogical(Params,'Calculate Heating',Found)
    IF( Calculate ) THEN
      CALL ListAddString( Params,'Variable', &
          TRIM(PressureName)//' Heating' )
    ELSE IF( .NOT. ListCheckPresent( Params,'Variable') ) THEN
      CALL Info(CName,'Defaulting field name to: ReynoldsPost')
      CALL ListAddString( Params,'Variable', &
          '-nooutput ReynoldsPost' )
    END IF
    ManningModel = ListGetLogical( Params,'Manning Model',Found )
    IF( ManningModel ) THEN
      CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),&
          TRIM(PressureName)//' corrected' )
    END IF
    Calculate = ListGetLogical(Params,'Calculate Force',Found)
    IF( Calculate ) THEN
      GivenDim = ListGetInteger(Params,'Calculate Force Dim',Found)
      IF( Dim == 1 .OR. GivenDim == 2 ) THEN
        dofs = 2
      ELSE
        dofs = 3
      END IF
      CALL Info(CName,'Creating "'//TRIM(PressureName)//'" Force with '&
          //I2S(dofs)//' components',Level=12)
      CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params), &
          '-dofs '//I2S(dofs)//' '//TRIM(PressureName)//' Force' )
    END IF
    Calculate = ListGetLogical(Params,'Calculate Flux',Found)
    IF( Calculate ) THEN
      GivenDim = ListGetInteger(Params,'Calculate Flux Dim',Found)
      IF( Dim == 1 .OR. GivenDim == 2 ) THEN
        dofs = 2
      ELSE
        dofs = 3
      END IF
      CALL Info(CName,'Creating "'//TRIM(PressureName)//' Flux" with '&
          //I2S(dofs)//' components',Level=12)
      CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params), &
          '-dofs '//I2S(dofs)//' '//TRIM(PressureName)//' Flux' )
    END IF
    Calculate = ListGetLogical(Params,'Calculate Mean Velocity',Found)
    IF( Calculate ) THEN
      GivenDim = ListGetInteger(Params,'Calculate Mean Velocity Dim',Found)
      IF( Dim == 1 .OR. GivenDim == 2 ) THEN
        dofs = 2
      ELSE
        dofs = 3
      END IF
      CALL Info(CName,'Creating "'//TRIM(PressureName)//' mean velocity" with '&
          //I2S(dofs)//' components',Level=12)
      CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params), &
          '-dofs '//I2S(dofs)//' '//TRIM(PressureName)//' Mean Velocity' )
    END IF
    CALL ListAddInteger( Params, 'Time derivative order', 0 )
    CALL ListAddNewString(Params,'Linear System Solver','Iterative')
    CALL ListAddNewString(Params,'Linear System Iterative Method','cg')
    CALL ListAddNewString(Params,'Linear System Preconditioning','ILU0')
    CALL ListAddNewInteger(Params,'Linear System Max Iterations',500)
    CALL ListAddNewInteger(Params,'Linear System Residual Output',10)
    CALL ListAddNewConstReal(Params,'Linear System Convergence Tolerance',1.0e-10_dp)
  end subroutine rpp_init_wrap

  ! --- transcribed kernel (statement-for-statement) ---

  SUBROUTINE LocalMatrix(STIFF, FORCE, Element, n, nd, Nodes)
    implicit none
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:)
    INTEGER :: n, nd
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), TARGET :: Element
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3), DetJ
    REAL(KIND=dp) :: x,y,z,Metric(3,3),SqrtMetric,Symb(3,3,3),dSymb(3,3,3,3)
    REAL(KIND=dp) :: U, V, W, S, A
    REAL(KIND=dp) :: Spres, Sslide
    REAL(KIND=dp) :: Normal(3), Velo(3), TangentVelo(3), Pres, Gap, GradPres(3), &
        Visc, mfp, Kn, TotPres, source, Radius(3)
    LOGICAL :: Stat
    INTEGER :: i,p,q,t,NBasis, CoordSys
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
    CoordSys = CurrentCoordinateSystem()
    Metric = 0.0_dp
    Metric(1,1) = 1.0_dp
    Metric(2,2) = 1.0_dp
    Metric(3,3) = 1.0_dp
    NBasis = nd
    IntegStuff = GaussPoints( Element )
    DO t=1,IntegStuff % n
      u = IntegStuff % u(t)
      v = IntegStuff % v(t)
      w = IntegStuff % w(t)
      s = IntegStuff % s(t)
      stat = ElementInfo(Element, Nodes, u, v, w, DetJ, Basis, dBasisdx)
      s = s * DetJ
      IF ( CoordSys /= Cartesian .OR. CalculateMoment ) THEN
        x = SUM( Nodes % x(1:n) * Basis(1:n) )
        y = SUM( Nodes % y(1:n) * Basis(1:n) )
        z = SUM( Nodes % z(1:n) * Basis(1:n) )
      END IF
      IF ( CoordSys /= Cartesian ) THEN
        CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,X,Y,Z )
        s = s * SqrtMetric
      END IF
      IF( CalculateMoment ) THEN
        Radius(1) = X - MomentAbout(1)
        Radius(2) = Y - MomentAbout(2)
        Radius(3) = Z - MomentAbout(3)
      END IF
      IF( ASSOCIATED( Element % BoundaryInfo)) THEN
        Normal = NormalVector( Element,Nodes,u,v,.TRUE. )
      ELSE
        Normal = NormalVector( Element,Nodes,u,v,.FALSE. )
      END IF
      IF( UseVelocity ) THEN
        DO i=1,3
          Velo(i) = SUM( Basis(1:n) * Velocity(i,1:n) )
        END DO
        TangentVelo = Velo - SUM(Normal * Velo) * Normal
      ELSE
        DO i=1,3
          TangentVelo(i) = SUM( Basis(1:n) * Velocity(i,1:n) )
        END DO
      END IF
      Pres = SUM(Basis(1:n) * ElemPressure(1:n))
      TotPres = Pres + ReferencePressure
      Gap = SUM(Basis(1:n) * GapHeight(1:n))
      DO i=1,3
        GradPres(i) = SUM(dBasisdx(1:n,i) * ElemPressure(1:n))
      END DO
      SELECT CASE (ViscosityType)
      CASE (Viscosity_Newtonian)
        Visc = SUM(Basis(1:n) * Viscosity(1:n))
      CASE (Viscosity_Rarefied)
        Visc = SUM(Basis(1:n) * Viscosity(1:n))
        mfp = mfp0 * ReferencePressure / TotPres
        Kn = mfp / ABS(Gap)
        Visc = Visc / (1+9.638_dp*Kn**1.159_dp)
      END SELECT
      TotPres = TotPres - AmbientPres
      Sslide = 0.0_dp
      Spres = 0.0_dp
      SELECT CASE( Mode )
      CASE( 1 )
        Spres = -TotPres * Normal( Component )
        IF( OpposingWall ) Spres = -Spres
        Spres = Spres - Gap * GradPres( Component ) / 2
        Sslide = Visc * TangentVelo(Component ) / Gap
        source = Spres + Sslide
        Pforce(Component) = Pforce(Component) + s * Spres
        Vforce(Component) = Vforce(Component) + s * Sslide
        IF( ASSOCIATED(SensVar) ) THEN
          Sforce(Component) = Sforce(Component) + s * &
              SUM(Basis(1:n) * SensPressure(1:n)) * Normal( Component )
        END IF
        IF( CalculateMoment ) THEN
          IF( Component == 1 ) THEN
            Moment(2) = Moment(2) + Radius(3) * s * source
            Moment(3) = Moment(3) - Radius(2) * s * source
          ELSE IF ( Component == 2 ) THEN
            Moment(1) = Moment(1) - Radius(3) * s * source
            Moment(3) = Moment(3) + Radius(1) * s * source
          ELSE IF ( Component == 3 ) THEN
            Moment(1) = Moment(1) + Radius(2) * s * source
            Moment(2) = Moment(2) - Radius(1) * s * source
          END IF
        END IF
      CASE( 2 )
        Spres = - (Gap**3 / (12 * Visc) ) * GradPres(Component)
        Sslide = Gap * TangentVelo(Component) / 2
        source = Spres + Sslide
      CASE( 3 )
        Spres = - (Gap**2 / (12 * Visc) ) * GradPres(Component)
        Sslide = TangentVelo(Component) / 2
        source = Spres + Sslide
      CASE( 4 )
        Spres = (Gap**3 / (12 * Visc) ) * SUM(GradPres *GradPres )
        Sslide = (Visc / Gap) * SUM(TangentVelo * TangentVelo )
        source = Spres + Sslide
        HeatPres = HeatPres + s * Spres
        HeatSlide = HeatSlide + s * Sslide
      CASE DEFAULT
        CALL Fatal(Caller,'Unknow mode: '//I2S(Mode))
      END SELECT
      DO p=1,NBasis
        DO q=1,NBasis
          A = Basis(q) * Basis(p)
          STIFF(p,q) = STIFF(p,q) + s * A
        END DO
        FORCE(p) = FORCE(p) + s * Basis(p) * source
      END DO
    END DO
  END SUBROUTINE LocalMatrix

end module reynolds_post_acc
