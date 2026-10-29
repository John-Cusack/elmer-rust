! T022-A AdvectionDiffusion port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/AdvectionDiffusion.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays, C_INT
! codes). All Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER
! dummies, assumed-shape/descriptor dummies, hidden lengths, list-directed
! timing reports — stays inside Fortran, so the Rust side is pure C calls.
! The wrapper module mirrors AdvectionDiffusionSolver's own entry locals
! (the SAVE'd assembly arrays plus every per-invocation local, which must
! persist across the per-element wrapper calls within one invocation) with
! identical lifetimes, and carries the six CONTAINS internal procedures as
! verbatim module procedures (their host-associated Diffusivity resolves to
! the same module SAVE array). All floating-point work (the Galerkin
! assembly, SUM reductions, the CPUTime-based reports) and every formatted
! WRITE stay in Fortran (same compiler), so bit-identity holds by
! construction wherever the call sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module advectiondiffusion_acc
  use iso_c_binding
  use SolverBasics
  use Differentials
  use DefUtils
  use MaterialModels
  use Integration
  use messages
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()

  ! -- entry locals mirrored as module SAVE (original SAVE statement plus
  ! -- every other entry local, which must survive across wrapper calls) --
  integer, save :: i,j,k,m,n,pn,t,tmax,istat,bf_id,body_id,lbody,rbody,eq_id,MDOFs,NSDOFs
  integer, save :: NonlinearIter,CoordinateSystem,outbody
  integer, save :: CompressibilityModel

  type(matrix_t),pointer,save  :: StiffMatrix => null()
  type(nodes_t), target, save  :: ElementNodes, ParentNodes
  type(element_t),pointer,save :: CurrentElement => null(), Parent => null()
  type(mesh_t), pointer, save  :: Mesh => null()
  real(kind=dp), save :: Norm, RelativeChange
  real(kind=dp), pointer, save :: Species(:) => null(), ForceVector(:) => null(), &
      MeshVelocity(:) => null(), Hwrk(:,:,:) => null(), Temperature(:) => null(), &
      EPotential(:) => null(), FlowSolution(:) => null()
  integer, pointer, save :: NodeIndexes(:) => null()
  logical, save :: Stabilize = .FALSE., Bubbles = .FALSE., GotIt = .FALSE., &
      GotIt2 = .FALSE., AbsoluteMass = .FALSE.
  logical, save :: AllocationsDone = .FALSE., ScaledToSolubility = .FALSE.
  logical, save :: ErrorWritten = .FALSE.
  logical, save :: TransientSim = .FALSE.

  type(valuelist_t), pointer, save :: Material => null(), BC => null(), Eq => null()

  integer, pointer, save :: SpeciesPerm(:) => null(), MeshPerm(:) => null()
  real(kind=dp), allocatable, save :: LocalMassMatrix(:,:), SoretDiffusivity(:), &
    NEConst(:), LocalStiffMatrix(:,:),Load(:),Diffusivity(:,:,:), &
              C0(:),C1(:),CT(:),C2(:,:,:),LocalForce(:), TimeForce(:)
  character(len=MAX_NAME_LEN), save :: ConvectionFlag = '', HeatSolName = '', &
      PhiSolName = '', ConvectName = ''
  type(variable_t), pointer, save :: TempSol => null(),PhiSol => null(), &
      FlowSol => null(),MeshSol => null()

  integer, pointer, save :: TempPerm(:) => null(),PotentialPerm(:) => null(), &
      FlowPerm(:) => null()
  real(kind=dp), save :: Relax = 1.0_dp, dt = 0.0_dp
  real(kind=dp), save :: SpecificHeatRatio, ReferencePressure, Ratio, MaxSol, IonCharge
  real(kind=dp), allocatable, save :: U(:), V(:), W(:), LocalTemperature(:), &
      LocalPotential(:), Density(:), HeatCapacity(:), Pressure(:),GasConstant(:), &
      SpeciesTransferCoeff(:), SExt(:), MU(:), MV(:), MW(:)

  real(kind=dp), save :: at,totat,st,totst,at0
  real(kind=dp), save :: Nrm(3), Nrm2(3)
  real(kind=dp), save :: Mass, PreviousMass

  character(len=MAX_NAME_LEN), save :: CompressibilityFlag = '', StabilizeFlag = ''
  character(len=MAX_NAME_LEN), save :: EquationName = ''
  character(len=MAX_NAME_LEN), save :: ConcentrationUnits = ''

contains

  ! -- AdvectionDiffusionSolver_init body ---------------------------------
  subroutine ada_init_begin(h_solver) bind(c, name="ada_init_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine ada_init_begin

  subroutine ada_init_body() bind(c, name="ada_init_body")
    type(ValueList_t), pointer :: Params
    Params => GetSolverParams()
    call ListAddInteger( Params,'Time Derivative Order', 1 )
  end subroutine ada_init_body

  ! -- entry state ---------------------------------------------------------
  subroutine ada_begin(h_model, h_solver, dt_in, transient_in) &
      bind(c, name="ada_begin")
    type(c_ptr), value :: h_model, h_solver
    real(c_double), value :: dt_in
    integer(c_int), value :: transient_in
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    dt = dt_in
    TransientSim = (transient_in /= 0)
  end subroutine ada_begin

  ! -- preamble: matrix guard, headers, variable setup, allocation, method --
  ! -- flags; returns 1 when the body must RETURN early --------------------
  function ada_preamble() result(code) bind(c, name="ada_preamble")
    integer(c_int) :: code
    code = 0
     IF ( .NOT. ASSOCIATED( acc_solver % Matrix ) ) THEN
       code = 1
       RETURN
     END IF

     CALL Info( 'AdvectionDiffusion', ' ', Level=4 )
     CALL Info( 'AdvectionDiffusion', &
         '-------------------------------------', Level=4 )
     CALL Info( 'AdvectionDiffusion', 'Solving for species transport', Level=4 )
     CALL Info( 'AdvectionDiffusion', &
         '-------------------------------------', Level=4 )

     CoordinateSystem = CurrentCoordinateSystem()
     Mesh => acc_solver % Mesh

     Species     => acc_solver % Variable % Values
     SpeciesPerm => acc_solver % Variable % Perm

     IF ( ALL( SpeciesPerm == 0) ) THEN
       code = 1
       RETURN
     END IF

     TempSol => VariableGet( Mesh % Variables, 'Temperature' )
     IF ( ASSOCIATED( TempSol ) ) THEN
       TempPerm    => TempSol % Perm
       Temperature => TempSol % Values
     END IF

     PhiSol => VariableGet( Mesh % Variables, 'Potential' )
     IF ( ASSOCIATED( PhiSol ) ) THEN
       PotentialPerm => PhiSol % Perm
       EPotential    => PhiSol % Values
     END IF

     FlowSol => VariableGet( Mesh % Variables, 'Flow Solution' )
     IF ( ASSOCIATED( FlowSol ) ) THEN
       FlowPerm     => FlowSol % Perm
       NSDOFs       =  FlowSol % DOFs
       FlowSolution => FlowSol % Values
     END IF

     MeshSol => VariableGet( Mesh % Variables, 'Mesh Velocity' )
     NULLIFY( MeshVelocity )
     IF ( ASSOCIATED(MeshSol ) ) THEN
       MDOFs    =  MeshSol % DOFs
       MeshPerm => MeshSol % Perm
       MeshVelocity => MeshSol % Values
     END IF

     StiffMatrix => acc_solver % Matrix
     ForceVector => StiffMatrix % RHS

     Norm = acc_solver % Variable % Norm
     IF ( .NOT. AllocationsDone ) THEN
       N = Mesh % MaxElementNodes
       ALLOCATE( U( N ),  V( N ), W( N ), &
                 MU( N ),MV( N ),MW( N ), &
                 LocalTemperature( N ),   &
                 LocalPotential( N ),     &
                 Pressure( N ),           &
                 Density( N ),            &
                 ElementNodes % x( N ),   &
                 ElementNodes % y( N ),   &
                 ElementNodes % z( N ),   &
                 ParentNodes % x( N ),   &
                 ParentNodes % y( N ),   &
                 ParentNodes % z( N ),   &
                 LocalForce( 2*N ),         &
                 TimeForce( 2*N ),         &
                 LocalMassMatrix( 2*N,2*N ),  &
                 LocalStiffMatrix( 2*N,2*N ), &
                 Load( N ), Diffusivity( 3,3,N ), &
                 SoretDiffusivity( N ),   &
                 NEConst( N ),            &
                 HeatCapacity( N ),       &
                 GasConstant( N ),        &
                 SpeciesTransferCoeff( N ), &
                 SExt( N ),              &
                 C0( N ), C1( N ), CT( N ), C2( 3,3,N ),STAT=istat )

       NULLIFY( HWrk)

       IF ( istat /= 0 ) THEN
         CALL Fatal( 'AdvectionDiffusion', 'Memory allocation error.' )
       END IF

       Mass = 0.0d0

       AllocationsDone = .TRUE.
     END IF
     Bubbles = .FALSE.
     Stabilize = .FALSE.

     StabilizeFlag = GetString( acc_solver % Values, 'Stabilization Method', GotIt )
     IF( GotIt ) THEN
       SELECT CASE(StabilizeFlag)
       CASE('stabilized')
         Stabilize = .TRUE.
       CASE('bubbles')
         Bubbles = .TRUE.
       CASE DEFAULT
         CALL Fatal('AdvectionDiffusion','Unknown stabilization method: '//TRIM(StabilizeFlag))
       END SELECT
     ELSE
       Stabilize = ListGetLogical( acc_solver % Values,'Stabilize',GotIt )
       Bubbles = ListGetLogical( acc_solver % Values,'Bubbles',GotIt2 )

       IF( .NOT. ( GotIt .OR. GotIt2 ) ) THEN
         CALL Info('AdvectionDiffusion','Defaulting stabilization to bubbles',Level=10)
         Bubbles = .TRUE.
       ELSE IF( Stabilize .AND. Bubbles ) THEN
         CALL Fatal('AdvectionDiffusion','Choose either Bubbles or Stabilize!')
       END IF
     END IF

     NonlinearIter = ListGetInteger( acc_solver % Values, &
        'Nonlinear System Max Iterations',GotIt )
     IF ( .NOT.GotIt ) NonlinearIter = 1

     Relax = GetCReal( acc_solver % Values, &
         'Nonlinear System Relaxation Factor',GotIt )
     IF ( .NOT.GotIt ) Relax = 1

     EquationName = ListGetString( acc_solver % Values, 'Equation' )

     totat = 0.0d0
     totst = 0.0d0
  end function ada_preamble

  function ada_nonlinear_max() result(nm) bind(c, name="ada_nonlinear_max")
    integer(c_int) :: nm
    nm = NonlinearIter
  end function ada_nonlinear_max

  subroutine ada_iter_start(iter) bind(c, name="ada_iter_start")
    integer(c_int), value :: iter
       at  = CPUTime()
       at0 = RealTime()

       CALL Info( 'AdvectionDiffusion', &
             '-------------------------------------', Level=4 )
       WRITE( Message, * ) 'Iteration step: ', iter
       CALL Info( 'AdvectionDiffusion', Message, Level=4 )
       CALL Info( 'AdvectionDiffusion', &
             '-------------------------------------', Level=4 )

       CALL DefaultInitialize()
       body_id = -1
       NULLIFY(Material)

       tmax = acc_solver % NumberOfActiveElements
       CALL Info( 'AdvectionDiffusion','Bulk Assembly')
  end subroutine ada_iter_start

  function ada_tmax() result(nm) bind(c, name="ada_tmax")
    integer(c_int) :: nm
    nm = tmax
  end function ada_tmax


  ! -- one bulk element: full per-element body incl. update ------------------
  subroutine ada_bulk(t) bind(c, name="ada_bulk")
    integer(c_int), value :: t
         IF ( RealTime() - at0 > 1.0 ) THEN
           WRITE(Message,'(a,i3,a)' ) '   Assembly: ', &
               INT(100.0 - 100.0 * (tmax-t) / (1.0*tmax)), ' % done'
           CALL Info( 'AdvectionDiffusion', Message, Level=5 )
           at0 = RealTime()
         END IF

         CurrentElement => GetActiveElement(t)
         IF ( CurrentElement % BodyId /= body_id ) THEN
           body_id = CurrentElement % Bodyid
           eq_id = ListGetInteger( acc_model % Bodies(body_id) % Values,'Equation')
           Eq => acc_model % Equations(eq_id) % Values

           ConvectionFlag = ListGetString( Eq, &
               TRIM(ComponentName(acc_solver % Variable)) //' Convection', GotIt )
           IF ( .NOT. GotIt ) &
             ConvectionFlag = ListGetString( Eq,'Convection', GotIt )

           ScaledToSolubility = .FALSE.
           ConcentrationUnits = ListGetString( Eq, 'Concentration Units', GotIt )
           IF ( .NOT.GotIt ) AbsoluteMass = .FALSE.
           IF (ConcentrationUnits == 'absolute mass') THEN
              AbsoluteMass = .TRUE.
           ELSE IF (ConcentrationUnits == 'mass to max solubility' ) THEN
              AbsoluteMass = .TRUE.
              ScaledToSolubility = .TRUE.
           ELSE
              AbsoluteMass = .FALSE.
           END IF

           k = ListGetInteger( acc_model % Bodies( body_id) % Values, 'Material' )
           Material => acc_model % Materials(k) % Values

           HeatSolName = ListGetString( Material, &
               'Temperature Field Variable', GotIt )
           IF ( Gotit ) THEN
             TempSol => VariableGet( Mesh % Variables, &
                 TRIM( HeatSolName ) )
             IF ( ASSOCIATED( TempSol ) ) THEN
               TempPerm     => TempSol % Perm
               Temperature  => TempSol % Values
             ELSE
               WRITE( Message, * ) 'No temperature variable ' &
                   // TRIM( HeatSolName ) // ' available'
               CALL Fatal( 'AdvectionDiffusion', Message )
             END IF
           END IF

           PhiSolName = ListGetString( Material, &
               TRIM(ComponentName(acc_solver % Variable)) // &
               ' EPotential Field Variable', GotIt )
           IF ( .NOT. GotIt ) &
             PhiSolName = ListGetString( Material, &
                 'EPotential Field Variable', GotIt )
           IF ( Gotit ) THEN
             PhiSol => VariableGet( Mesh % Variables, &
                 TRIM( PhiSolName ) )
             IF ( ASSOCIATED( PhiSol ) ) THEN
               PotentialPerm => PhiSol % Perm
               EPotential    => PhiSol % Values
             ELSE
               WRITE( Message, * ) 'No electrical potential variable ' &
                   // TRIM( PhiSolName ) // ' available'
               CALL Fatal( 'AdvectionDiffusion', Message )
             END IF
           END IF

           ConvectName = ListGetString( Eq, &
               TRIM(ComponentName(acc_solver % Variable)) // &
               ' Convection Field Variable', GotIt )
           IF ( .NOT. GotIt ) &
             ConvectName = ListGetString( Material, &
                 'Convection Field Variable', GotIt )
           IF ( GotIt ) THEN
             FlowSol => VariableGet( Mesh % Variables, &
                 TRIM( ConvectName ) )
             IF ( ASSOCIATED( FlowSol ) ) THEN
               FlowPerm     => FlowSol % Perm
               NSDOFs       =  FlowSol % DOFs
               FlowSolution => FlowSol % Values
             ELSE
               WRITE( Message, * ) 'No convection  variable ' // &
                   TRIM( ConvectName ) // ' available'
               CALL Fatal( 'AdvectionDiffusion', Message )
             END IF
           END IF

           CompressibilityFlag = ListGetString( Material, &
               'Compressibility Model', GotIt)
           IF ( .NOT.GotIt ) CompressibilityModel = Incompressible

           SELECT CASE( CompressibilityFlag )

             CASE( 'incompressible' )
             CompressibilityModel = Incompressible

             CASE( 'user defined 1' )
             CompressibilityModel = UserDefined1

             CASE( 'user defined 2' )
             CompressibilityModel = UserDefined2

             CASE( 'perfect gas equation 1' )
             CompressibilityModel = PerfectGas1

             CASE( 'perfect gas equation 2' )
             CompressibilityModel = PerfectGas2

             CASE( 'perfect gas equation 3' )
             CompressibilityModel = PerfectGas3

           CASE DEFAULT
             CompressibilityModel = Incompressible
           END SELECT
         END IF

         n = CurrentElement % TYPE % NumberOfNodes
         NodeIndexes => CurrentElement % NodeIndexes

         ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
         ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
         ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)
         IF ( ASSOCIATED( TempSol ) ) THEN
           LocalTemperature(1:n) = Temperature( TempPerm(NodeIndexes) )
         ELSE
           LocalTemperature(1:n) = ListGetReal(  Material, &
                'Reference Temperature',n,NodeIndexes,GotIt )
         ENDIF
         IF ( ASSOCIATED( PhiSol ) ) THEN
           LocalPotential(1:n) = EPotential( PotentialPerm(NodeIndexes) )
         ELSE
           LocalPotential(1:n) = ListGetReal(  Material, &
                'Reference EPotential',n,NodeIndexes,GotIt )
         ENDIF

         CALL ListGetRealArray( Material,  &
             TRIM(ComponentName(acc_solver % Variable)) // &
                     ' Diffusivity', Hwrk, n, NodeIndexes )

         Diffusivity = 0.0d0
         IF ( SIZE(Hwrk,1) == 1 ) THEN
           DO i=1,3
             Diffusivity( i,i,1:n ) = Hwrk( 1,1,1:n )
           END DO
         ELSE IF ( SIZE(Hwrk,2) == 1 ) THEN
           DO i=1,MIN(3,SIZE(Hwrk,1))
             Diffusivity(i,i,1:n) = Hwrk(i,1,1:n)
           END DO
         ELSE
           DO i=1,MIN(3,SIZE(Hwrk,1))
             DO j=1,MIN(3,SIZE(Hwrk,2))
               Diffusivity( i,j,1:n ) = Hwrk(i,j,1:n)
             END DO
           END DO
         END IF

         SoretDiffusivity(1:n) = ListGetReal( Material, &
              TRIM(ComponentName(acc_solver % Variable)) // &
              ' Soret Diffusivity', n, NodeIndexes, GotIt )
         IF ( .NOT. GotIt )  SoretDiffusivity = 0.0d0

         IonCharge = ListGetConstReal( Material, &
              TRIM(ComponentName(acc_solver % Variable)) // &
              ' Ion Charge', GotIt )
         IF ( .NOT. GotIt ) THEN
           NEConst   = 0._dp
           IonCharge = 0._dp
         ELSE
           NEConst(1:n) = IonCharge / LocalTemperature(1:n) * &
               ListGetConstReal( acc_model % Constants,'Unit Charge') / &
               ListGetConstReal( acc_model % Constants,'Boltzmann Constant')
         END IF

         HeatCapacity(1:n) = ListGetReal( Material,'Heat Capacity', &
             n,NodeIndexes,GotIt )

         IF ( CompressibilityModel >= PerfectGas1 .AND. &
             CompressibilityModel <= PerfectGas3 ) THEN
           SpecificHeatRatio = ListGetConstReal( Material, &
               'Specific Heat Ratio', GotIt )
           IF ( .NOT.GotIt ) SpecificHeatRatio = 1.4d0
           GasConstant(1:n) = ( SpecificHeatRatio - 1.d0 ) * &
               HeatCapacity(1:n) / SpecificHeatRatio
         ELSE
           Density(1:n) = ListGetReal( Material,'Density',n,NodeIndexes )
         END IF
         ReferencePressure = 0.0_dp
         IF ( CompressibilityModel /= Incompressible ) THEN
           ReferencePressure = ListGetConstReal( Material, &
               'Reference Pressure', GotIt)
         END IF

         Load = 0.d0
         Pressure = 0.d0

         IF ( ConvectionFlag == 'constant' ) THEN
           U = ListGetReal( Eq, TRIM(ComponentName(acc_solver % Variable)) // &
               ' Convection Velocity 1',n,NodeIndexes,GotIt )
           IF ( .NOT. GotIt ) &
               U = ListGetReal( Material,'Convection Velocity 1',n,NodeIndexes, GotIt)
           V = ListGetReal( Eq, TRIM(ComponentName(acc_solver % Variable)) // &
               ' Convection Velocity 2',n,NodeIndexes,GotIt )
           IF ( .NOT. GotIt ) &
             V = ListGetReal( Material,'Convection Velocity 2',n,NodeIndexes, GotIt)
           W = ListGetReal( Eq, TRIM(ComponentName(acc_solver % Variable)) // &
                ' Convection Velocity 3',n,NodeIndexes,GotIt )
           IF ( .NOT. GotIt ) &
             W = ListGetReal( Material,'Convection Velocity 3',n,NodeIndexes, GotIt)
         ELSE IF ( ConvectionFlag == 'computed' ) THEN

           IF( .NOT. ASSOCIATED( FlowSol ) ) THEN
             CALL Fatal('AdvectionDiffusion','Give > Convection Field Variable <')
           END IF

           DO i=1,n
             k = FlowPerm(NodeIndexes(i))
             IF ( k > 0 ) THEN

               SELECT CASE( CompressibilityModel )
                 CASE( PerfectGas1,PerfectGas2,PerfectGas3 )
                 Pressure(i) = FlowSolution(NSDOFs*k) + ReferencePressure
                 Density(i)  = Pressure(i) / &
                       ( GasConstant(i) * LocalTemperature(i) )
               END SELECT

               SELECT CASE( NSDOFs )
                 CASE(3)
                 U(i) = FlowSolution( NSDOFs*k-2 )
                 V(i) = FlowSolution( NSDOFs*k-1 )
                 W(i) = 0.0D0

               CASE(4)
                 U(i) = FlowSolution( NSDOFs*k-3 )
                 V(i) = FlowSolution( NSDOFs*k-2 )
                 W(i) = FlowSolution( NSDOFs*k-1 )
               END SELECT

             ELSE
               U(i) = 0.0d0
               V(i) = 0.0d0
               W(i) = 0.0d0
             END IF
           END DO

         ELSE
           U = 0.0d0
           V = 0.0d0
           W = 0.0d0
         END IF

         IF ( ConvectionFlag == 'constant' .OR. ConvectionFlag == 'computed' &
              .OR. IonCharge /= 0.0d0 ) THEN
            IF (.NOT. AbsoluteMass) THEN
               C1 = Density
            ELSE
               C1 = 1.d0
            END IF
         ELSE
           C1 = 0.0d0
         END IF

         IF (.NOT. AbsoluteMass) THEN
            CT = Density
            DO i=1,3
               DO j=1,3
                  C2(i,j,:) = Density * Diffusivity(i,j,:)
               END DO
            END DO
            WHERE ( Density(1:n) /= 0.0d0 )
               NEConst(1:n) = NEConst(1:n) / Density(1:n)
            ELSEWHERE
               NEConst(1:n) = 0.0d0
            END WHERE
         ELSE
            CT = 1.0d0
            DO i=1,3
               DO j=1,3
                  C2(i,j,:) = Diffusivity(i,j,:)
               END DO
            END DO
         END IF

         MU  = 0.0d0
         MV  = 0.0d0
         MW  = 0.0d0
         IF ( ASSOCIATED( MeshVelocity ) ) THEN
            DO i=1,n
              IF ( MeshPerm( NodeIndexes(i) ) > 0 ) THEN
                 MU(i) = MeshVelocity( MDOFs*(MeshPerm(NodeIndexes(i))-1)+1 )
                 MV(i) = MeshVelocity( MDOFs*(MeshPerm(NodeIndexes(i))-1)+2 )
                 IF ( MDOFs > 2 ) THEN
                    MW(i) = MeshVelocity( MDOFs*(MeshPerm(NodeIndexes(i))-1)+3 )
                 END IF
              END IF
            END DO
         END IF

         IF (.NOT. AbsoluteMass) THEN
            C0 = Density
         ELSE
            C0 = 1.0d0
         END IF
         IF ( ScaledToSolubility ) THEN
            MaxSol = ListGetConstReal( Material, &
              TRIM(ComponentName(acc_solver % Variable)) // &
                 ' Maximum Solubility', GotIt )
            IF ( .NOT. GotIT ) THEN
               WRITE( Message, * ) 'Maximum solubility not defined in body : ', &
                    CurrentElement % BodyId
               CALL Fatal( 'AdvectionDiffusion', Message )
            END IF
            C0 = C0 / MaxSol
         END IF

         bf_id = ListGetInteger( acc_model % Bodies(CurrentElement % BodyId) % &
           Values, 'Body Force', GotIt, 1, acc_model % NumberOfBodyForces )

         IF ( GotIt ) THEN

           IF ( ScaledToSolubility .AND. ListGetLogical( acc_model % &
               BodyForces(bf_id) % Values, 'Physical Units', GotIt ) ) THEN

             Ratio = ListGetConstReal( Material, &
              TRIM(ComponentName(acc_solver % Variable)) // &
                 ' Maximum Solubility', GotIt )
             IF ( .NOT. GotIT ) THEN
               WRITE( Message, * ) 'Maximum solubility not defined in body : ', &
                   CurrentElement % BodyId
               CALL Fatal( 'AdvectionDiffusion', Message )
             END IF
           ELSE
             Ratio = 1.0d0
           END IF
           Load(1:n) = Load(1:n) + &
               ListGetReal( acc_model % BodyForces(bf_id) % Values,  &
               TRIM(ComponentName(acc_solver % Variable)) // &
                 ' Diffusion Source',n,NodeIndexes,gotIt ) / Ratio
         END IF

         IF ( CoordinateSystem == Cartesian ) THEN
           CALL DiffuseConvectiveCompose( &
               LocalMassMatrix, LocalStiffMatrix, LocalForce, Load, &
               CT, C0, C1, C2, LocalTemperature, LocalPotential, U(1:n), V(1:n), W(1:n), MU, MV, MW, &
               SoretDiffusivity, NEConst, (CompressibilityModel /= Incompressible), &
               AbsoluteMass,Stabilize, Bubbles, CurrentElement, n, ElementNodes )
         ELSE
           CALL DiffuseConvectiveGenCompose( &
               LocalMassMatrix, LocalStiffMatrix, LocalForce, Load, &
               CT, C0, C1, C2, LocalTemperature, LocalPotential, U, V, W, MU, MV, MW, &
               SoretDiffusivity, NEConst, (CompressibilityModel /= Incompressible), &
               AbsoluteMass,Stabilize,CurrentElement, n, ElementNodes )
         END IF

         TimeForce  = 0.0_dp
         IF ( TransientSim ) THEN
           CALL Default1stOrderTime( LocalMassMatrix, LocalStiffMatrix, &
             LocalForce )
         END IF
         IF ( Bubbles .AND. ( ConvectionFlag == 'computed' .OR. &
              ConvectionFlag == 'constant' ) ) THEN
           CALL Condensate( N, LocalStiffMatrix,  LocalForce, TimeForce )
         END IF

         CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce )
  end subroutine ada_bulk

  subroutine ada_bulk_finish() bind(c, name="ada_bulk_finish")
    call DefaultFinishBulkAssembly()
  end subroutine ada_bulk_finish

  ! -- mixed bulk-boundary assembly (no-op unless ScaledToSolubility) --------
  subroutine ada_mixed() bind(c, name="ada_mixed")
       IF ( ScaledToSolubility ) THEN
         CALL Info( 'AdvectionDiffusion', 'Mixed bulk-boundary assembly')

         DO t=Mesh % NumberOfBulkElements + 1, &
             Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

           CurrentElement => Mesh % Elements(t)
           GotIt = .FALSE.
           DO i=1,acc_model % NumberOfBCs
             GotIt = CurrentElement % BoundaryInfo % Constraint == &
                 acc_model % BCs(i) % Tag
             IF( GotIt ) EXIT
           END DO
           IF( .NOT. GotIt ) CYCLE

           BC => acc_model % BCs(i) % Values
           IF ( .NOT. ListGetLogical( BC, TRIM(ComponentName(acc_solver % Variable)) &
               // ' Solubility Change Boundary', GotIt ) )  CYCLE

           acc_model % CurrentElement => Mesh % Elements(t)
           n = CurrentElement % TYPE % NumberOfNodes
           NodeIndexes => CurrentElement % NodeIndexes

           ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
           ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
           ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)

           body_id = ListGetInteger( BC,'Normal Target Body', GotIt )

           lbody = 0
           IF ( ASSOCIATED( CurrentElement % BoundaryInfo % Left ) ) &
               lbody = CurrentElement % BoundaryInfo % Left % BodyId
           rbody = 0
           IF ( ASSOCIATED( CurrentElement % BoundaryInfo % Right ) ) &
               rbody = CurrentElement % BoundaryInfo % Right % BodyId

           IF ( .NOT. GotIt ) THEN
             body_id = lbody
             outbody = CurrentElement % BoundaryInfo % OutBody
             CurrentElement % BoundaryInfo % OutBody = body_id

             Nrm = NormalVector( CurrentElement, ElementNodes, &
                 CurrentElement % TYPE % NodeU(1), CurrentElement % TYPE % NodeV(1), &
                 .TRUE. )
             CurrentElement % BoundaryInfo % OutBody = outbody
             Nrm2 = NormalVector( CurrentElement, ElementNodes, &
                 CurrentElement % TYPE % NodeU(1), CurrentElement % TYPE % NodeV(1), &
                 .TRUE. )
             IF ( SUM( Nrm(1:3) * Nrm2(1:3) ) < 0 ) body_id = rbody
           END IF

           k = ListGetInteger( acc_model % Bodies( body_id ) % Values, &
               'Material', minv=1, maxv=acc_model % NumberOfMaterials )
           Material => acc_model % Materials(k) % Values

           Ratio = ListGetConstReal( Material, &
               TRIM(ComponentName(acc_solver % Variable)) // &
               ' Maximum Solubility', GotIt )

           IF ( .NOT. GotIT ) THEN
             WRITE( Message, * ) 'No maximum solubility defined for material : ', k
             CALL Fatal( 'AdvectionDiffusion', Message )
           END IF

           IF ( lbody == body_id ) THEN
             k = ListGetInteger( acc_model % Bodies( rbody ) % Values, 'Material', &
                 minv=1,maxv=acc_model % NumberOfMaterials )
             Parent => CurrentElement % BoundaryInfo % Right
           ELSE
             k = ListGetInteger( acc_model % Bodies(lbody) % Values, 'Material', &
                 minv=1,maxv=acc_model % NumberOfMaterials )
             Parent => CurrentElement % BoundaryInfo % Left
           END IF
           Material => acc_model % Materials(k) % Values
           Ratio = ListGetConstReal( Material, &
               TRIM(ComponentName(acc_solver % Variable)) // &
               ' Maximum Solubility', GotIt ) / Ratio

           IF ( .NOT. GotIT ) THEN
             WRITE( Message, * ) 'No maximum solubility defined for material : ', k
             CALL Fatal( 'AdvectionDiffusion', Message )
           END IF
           Ratio = Ratio - 1.0d0

           CALL ListGetRealArray( Material,  &
               TRIM(ComponentName(acc_solver % Variable)) // &
               ' Diffusivity', Hwrk, n, NodeIndexes )

           Diffusivity = 0.0d0
           IF ( SIZE(Hwrk,1) == 1 ) THEN
             DO m=1,3
               Diffusivity( m,m,1:n ) = Hwrk( 1,1,1:n )
             END DO
           ELSE IF ( SIZE(Hwrk,2) == 1 ) THEN
             DO m=1,MIN(3,SIZE(Hwrk,1))
               Diffusivity(m,m,1:n) = Hwrk(m,1,1:n)
             END DO
           ELSE
             DO m=1,MIN(3,SIZE(Hwrk,1))
               DO j=1,MIN(3,SIZE(Hwrk,2))
                 Diffusivity( m,j,1:n ) = Hwrk(m,j,1:n)
               END DO
             END DO
           END IF

           pn = Parent % TYPE % NumberOfNodes

           ParentNodes % x(1:pn) = Mesh % Nodes % x(Parent % NodeIndexes)
           ParentNodes % y(1:pn) = Mesh % Nodes % y(Parent % NodeIndexes)
           ParentNodes % z(1:pn) = Mesh % Nodes % z(Parent % NodeIndexes)
           IF ( CoordinateSystem == Cartesian ) THEN
             CALL DiffuseConvectiveBBoundary( LocalStiffMatrix, Parent, pn, &
                 ParentNodes, Ratio, CurrentElement, n, ElementNodes )
           ELSE
             CALL DiffuseConvectiveGenBBoundary(LocalStiffMatrix, Parent, &
                 pn, ParentNodes, Ratio, CurrentElement,n , ElementNodes )
           END IF
           IF ( TransientSim ) THEN
             LocalMassMatrix = 0.0d0
             LocalForce = 0.0d0
             CALL Add1stOrderTime( LocalMassMatrix, LocalStiffMatrix, &
                 LocalForce,dt,pn,1,SpeciesPerm(Parent % NodeIndexes),acc_solver )
           END IF

           CALL UpdateGlobalEquations( StiffMatrix, LocalStiffMatrix, &
               ForceVector, LocalForce, pn, 1, SpeciesPerm(Parent % NodeIndexes) )
         END DO
       END IF
  end subroutine ada_mixed

  subroutine ada_boundary_start() bind(c, name="ada_boundary_start")
       CALL Info( 'AdvectionDiffusion','Boundary Assembly')
       ErrorWritten = .FALSE.
  end subroutine ada_boundary_start

  function ada_nbulk() result(nb) bind(c, name="ada_nbulk")
    integer(c_int) :: nb
    nb = Mesh % NumberOfBulkElements
  end function ada_nbulk

  function ada_nboundary() result(nb) bind(c, name="ada_nboundary")
    integer(c_int) :: nb
    nb = Mesh % NumberOfBoundaryElements
  end function ada_nboundary

  ! -- one boundary element (global index t): full body incl. update --------
  subroutine ada_boundary(t) bind(c, name="ada_boundary")
    integer(c_int), value :: t
         CurrentElement => Mesh % Elements(t)

         DO i=1,acc_model % NumberOfBCs
           IF ( CurrentElement % BoundaryInfo % Constraint == &
               acc_model % BCs(i) % Tag ) THEN

              acc_model % CurrentElement => Mesh % Elements(t)
              n = CurrentElement % TYPE % NumberOfNodes
              NodeIndexes => CurrentElement % NodeIndexes

              IF ( ANY( SpeciesPerm( NodeIndexes ) <= 0 ) ) CYCLE

              ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
              ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
              ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)

              SpeciesTransferCoeff = 0.0D0
              SExt = 0.0D0
              Load = 0.0D0

              BC => acc_model % BCs(i) % Values

              SpeciesTransferCoeff(1:n) = ListGetReal( BC, &
                  'Mass Transfer Coefficient', n, NodeIndexes, GotIt )
              IF ( GotIt ) THEN
                SExt(1:n) = ListGetReal( BC, &
                    'External Concentration', n, NodeIndexes, GotIt )

                IF ( .NOT. AbsoluteMass .OR. ScaledToSolubility ) THEN
                  IF ( .NOT. ErrorWritten ) THEN
                    CALL Error( 'AdvectionDiffusion', '--------------------' )
                    CALL Error( 'AdvectionDiffusion', &
                        'Mass transfer coefficient possible to use only with absolute mass concentrations' )
                    CALL Error( 'AdvectionDiffusion', &
                        'Ignoring mass transfer BC' )
                    CALL Error( 'AdvectionDiffusion', '--------------------' )
                    ErrorWritten = .TRUE.
                  END IF
                  SExt = 0.0d0
                  SpeciesTransferCoeff = 0.0d0
                END IF
              ELSE
                SExt(1:n) = 0.0d0
              END IF

              DO j=1,n
                Load(j) = Load(j) + SpeciesTransferCoeff(j) * SExt(j)
              END DO

              IF ( ScaledToSolubility .AND. &
                  ListGetLogical( BC, 'Physical Units', GotIt ) ) THEN

                Ratio = ListGetConstReal( Material, &
                    TRIM(ComponentName(acc_solver % Variable)) // &
                    ' Maximum Solubility', GotIt )
                IF ( .NOT. GotIT ) THEN
                  WRITE( Message, * ) 'No maximum solubility defined in body : ', &
                      CurrentElement % BodyId
                  CALL Fatal( 'AdvectionDiffusion', Message )
                END IF
              ELSE
                Ratio = 1.0d0
              END IF
              Load(1:n) = Load(1:n) + &
                  ListGetReal( BC, &
                  TRIM(ComponentName(acc_solver % Variable)) // &
                  ' Flux', n,NodeIndexes,gotIt ) / Ratio
              IF ( CoordinateSystem == Cartesian ) THEN
                CALL DiffuseConvectiveBoundary( LocalStiffMatrix,LocalForce, &
                    Load,SpeciesTransferCoeff,CurrentElement,n,ElementNodes )
              ELSE
                CALL DiffuseConvectiveGenBoundary(LocalStiffMatrix,LocalForce,&
                    Load,SpeciesTransferCoeff,CurrentElement,n,ElementNodes )
              END IF
              IF ( TransientSim ) THEN
                LocalMassMatrix = 0.0d0
                CALL Add1stOrderTime( LocalMassMatrix, LocalStiffMatrix, &
                    LocalForce,dt,n,1,SpeciesPerm(NodeIndexes),acc_solver )
              END IF

              CALL UpdateGlobalEquations( StiffMatrix, LocalStiffMatrix, &
                  ForceVector, LocalForce, n, 1, SpeciesPerm(NodeIndexes) )
          END IF
        END DO
  end subroutine ada_boundary

  subroutine ada_assembly_finish() bind(c, name="ada_assembly_finish")
      CALL DefaultFinishBoundaryAssembly()
      CALL DefaultFinishAssembly()
      CALL DefaultDirichletBCs()
      CALL Info( 'AdvectionDiffusion', 'Assembly done')
  end subroutine ada_assembly_finish

  ! -- solve + iteration tail; returns 1 when the loop must EXIT ------------
  function ada_solve_tail(iter, maxiter) result(code) bind(c, name="ada_solve_tail")
    integer(c_int), value :: iter, maxiter
    integer(c_int) :: code
      code = 0
      at = CPUTime() - at
      st = CPUTime()

      Norm = DefaultSolve()

      st = CPUTIme()-st
      totat = totat + at
      totst = totst + st
      WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Assembly: (s)', at, totat
      CALL Info( 'AdvectionDiffusion', Message, Level=5 )
      WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Solve:    (s)', st, totst
      CALL Info( 'AdvectionDiffusion', Message, Level=5 )
      RelativeChange = acc_solver % Variable % NonlinChange

      WRITE( Message, * ) 'Result Norm   : ',Norm
      CALL Info( 'AdvectionDiffusion', Message, Level=4 )
      WRITE( Message, * ) 'Relative Change : ',RelativeChange
      CALL Info( 'AdvectionDiffusion', Message, Level=4 )

      IF ( acc_solver % Variable % NonlinConverged == 1 ) code = 1
  end function ada_solve_tail

  ! -- mass-conservation check tail ------------------------------------------
  subroutine ada_mass_check() bind(c, name="ada_mass_check")
    IF( ListCheckPresent( acc_model % Simulation, 'Species Density') ) THEN
      PreviousMass = Mass
      Mass = VolumeIntegrate( acc_model, acc_solver % ActiveElements, 'Species Density' )
      PRINT *,'Species Mass: ',Mass
      IF ( TransientSim )  THEN
        PRINT *,'Mass Gain: ',Mass - PreviousMass
      END IF
    END IF
  end subroutine ada_mass_check

  ! -- CONTAINS internals (verbatim transcriptions, module procedures) --------

   SUBROUTINE DiffuseConvectiveCompose( MassMatrix,StiffMatrix,ForceVector,  &
       LoadVector,NodalCT,NodalC0,NodalC1,NodalC2,Temperature,EPotential, &
       Ux,Uy,Uz,MUx,MUy,MUz,SoretD,NEC,Compressible,AbsoluteMass, &
       Stabilize,UseBubbles,Element,n,Nodes )
     REAL(KIND=dp), DIMENSION(:)   :: ForceVector,Ux,Uy,Uz,MUx,MUy,MUz,LoadVector
     REAL(KIND=dp), DIMENSION(:,:) :: MassMatrix,StiffMatrix
     REAL(KIND=dp) :: Temperature(:), SoretD(:)
     REAL(KIND=dp) :: EPotential(:), NEC(:)
     REAL(KIND=dp) :: NodalC0(:),NodalC1(:),NodalCT(:),NodalC2(:,:,:),dT

     LOGICAL :: Stabilize,UseBubbles,Compressible,AbsoluteMass

     INTEGER :: n

     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t), TARGET :: Element

     REAL(KIND=dp) :: ddBasisddx(n,3,3)
     REAL(KIND=dp) :: Basis(2*n)
     REAL(KIND=dp) :: dBasisdx(2*n,3),SqrtElementMetric

     REAL(KIND=dp) :: Velo(3),dVelodx(3,3),Force

     REAL(KIND=dp) :: A,M
     REAL(KIND=dp) :: Load

     REAL(KIND=dp) :: VNorm,hK,mK
     REAL(KIND=dp) :: Lambda=1.0,Pe,Pe1,Pe2,Tau,x,y,z

     REAL(KIND=dp) :: SorD, GradTemp(3), SoretForce
     REAL(KIND=dp) :: GradPhi(3), NEConst

     INTEGER :: i,j,k,c,p,q,t,dim,N_Integ,NBasis

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     REAL(KIND=dp) :: s,u,v,w,DivVelo

     REAL(KIND=dp) :: C0,C00,C1,CT,C2(3,3),dC2dx(3,3,3),SU(n),SW(n)
     REAL(KIND=dp) :: NodalCThermal(n), CThermal
     REAL(KIND=dp) :: NodalCNP(n), CNP

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     LOGICAL :: stat,Convection,ConvectOrNPAndStabilize,Bubbles,ThermalDiffusion,NPDiffusion

     dim = CoordinateSystemDimension()
     c = dim + 1

     ForceVector = 0.0D0
     StiffMatrix = 0.0D0
     MassMatrix  = 0.0D0
     Load = 0.0D0
     Convection =  ANY( NodalC1 /= 0.0d0 )
     NBasis = n
     Bubbles = .FALSE.
     IF ( Convection .AND. .NOT. Stabilize .AND. UseBubbles ) THEN
        NBasis = 2*n
        Bubbles = .TRUE.
     END IF

     ThermalDiffusion = .FALSE.
     IF ( ANY( ABS( SoretD(1:n) ) > AEPS ) ) THEN
        ThermalDiffusion = .TRUE.
        NodalCThermal = NodalC0(1:n)
     END IF
     NPDiffusion = .FALSE.
     IF ( ANY( ABS( NEC(1:n) ) > AEPS ) ) THEN
        NPDiffusion = .TRUE.
        NodalCNP = NodalC0(1:n)
     END IF
     NodalC0 = 0.0d0

     IF ( Bubbles ) THEN
        IntegStuff = GaussPoints( element, Element % TYPE % GaussPoints2 )
     ELSE
        IntegStuff = GaussPoints( element )
     END IF
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

     ConvectOrNPAndStabilize = .FALSE.
     IF ( Stabilize .AND. ( Convection .OR. NPDiffusion ) ) THEN
       ConvectOrNPAndStabilize = .TRUE.
       hK = element % hK
       mK = element % StabilizationMK
     END IF

     DO t=1,N_Integ

       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)

       stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
             Basis,dBasisdx,ddBasisddx,ConvectOrNPAndStabilize,Bubbles )

       s = SqrtElementMetric * S_Integ(t)
       C0 = SUM( NodalC0(1:n) * Basis(1:n) )
       C1 = SUM( NodalC1(1:n) * Basis(1:n) )
       CT = SUM( NodalCT(1:n) * Basis(1:n) )

       DO i=1,dim
         DO j=1,dim
           C2(i,j) = SUM( NodalC2(i,j,1:n) * Basis(1:n) )
         END DO
       END DO
       Convection = .FALSE.
       IF ( C1 /= 0.0D0 ) THEN
          Convection = .TRUE.
          Velo = 0.0D0
          Velo(1) = SUM( (Ux(1:n)-MUx(1:n))*Basis(1:n) )
          Velo(2) = SUM( (Uy(1:n)-MUy(1:n))*Basis(1:n) )
          IF ( dim > 2 ) Velo(3) = SUM( (Uz(1:n)-MUz(1:n))*Basis(1:n) )

          IF ( Compressible .AND. AbsoluteMass ) THEN
            dVelodx = 0.0D0
            DO i=1,3
              dVelodx(1,i) = SUM( Ux(1:n)*dBasisdx(1:n,i) )
              dVelodx(2,i) = SUM( Uy(1:n)*dBasisdx(1:n,i) )
              IF ( dim > 2 ) dVelodx(3,i) = SUM( Uz(1:n)*dBasisdx(1:n,i) )
            END DO

            DivVelo = 0.0D0
            DO i=1,dim
              DivVelo = DivVelo + dVelodx(i,i)
            END DO
            C0 = DivVelo
          END IF

          IF ( NPDiffusion ) THEN
             NEConst = SUM ( NEC(1:n)*Basis(1:n) )
             GradPhi = 0.0d0
             DO i = 1, dim
                GradPhi(i) = SUM( dBasisdx(1:n,i) * EPotential(1:n) )
             END DO
             DO i = 1, dim
                DO j = 1, dim
                   Velo(i) = Velo(i) - C2(i,j) * NEConst * GradPhi(j)
                END DO
             END DO
          END IF

          IF ( Stabilize ) THEN
            VNorm = SQRT( SUM(Velo(1:dim)**2) )

            Pe  = MIN( 1.0D0, mK*hK*C1*VNorm/(2*ABS(C2(1,1))) )

            Tau = 0.0D0
            IF ( VNorm /= 0.0 ) THEN
               Tau = hK * Pe / (2 * C1 * VNorm)
            END IF

            DO i=1,dim
              DO j=1,dim
                DO k=1,dim
                  dC2dx(i,j,k) = SUM( NodalC2(i,j,1:n)*dBasisdx(1:n,k) )
                END DO
              END DO
            END DO

            DO p=1,N
              SU(p) = C0 * Basis(p)
              DO i = 1,dim
                SU(p) = SU(p) + C1 * dBasisdx(p,i) * Velo(i)
                DO j=1,dim
                  SU(p) = SU(p) - C2(i,j) * ddBasisddx(p,i,j)
                  SU(p) = SU(p) - dC2dx(i,j,j) * dBasisdx(p,i)
                END DO
              END DO

              SW(p) = C0 * Basis(p)
              DO i = 1,dim
                SW(p) = SW(p) + C1 * dBasisdx(p,i) * Velo(i)
                DO j=1,dim
                  SW(p) = SW(p) - C2(i,j) * ddBasisddx(p,i,j)
                  SW(p) = SW(p) - dC2dx(i,j,j) * dBasisdx(p,i)
                END DO
              END DO
            END DO
          END IF
        END IF

        DO p=1,NBasis
        DO q=1,NBasis
          M = CT * Basis(q) * Basis(p)
          A = C0 * Basis(q) * Basis(p)
          DO i=1,dim
            DO j=1,dim
              A = A + C2(i,j) * dBasisdx(q,i) * dBasisdx(p,j)
            END DO
          END DO

          IF ( Convection ) THEN
            DO i=1,dim
              A = A + C1 * Velo(i) * dBasisdx(q,i) * Basis(p)
            END DO
            IF ( Stabilize ) THEN
              A = A + Tau * SU(q) * SW(p)
              M = M + Tau * CT * Basis(q) * SW(p)
            END IF
          END IF

          StiffMatrix(p,q) = StiffMatrix(p,q) + s * A
          MassMatrix(p,q)  = MassMatrix(p,q)  + s * M
        END DO
        END DO

        Force = SUM( LoadVector(1:n)*Basis(1:n) )

        DO p=1,NBasis
          Load = Basis(p)
          IF ( ConvectOrNPAndStabilize ) Load = Load + Tau * SW(p)
          ForceVector(p) = ForceVector(p) + s * Force * Load
        END DO

        IF ( ThermalDiffusion ) THEN

           CThermal = SUM( NodalCThermal(1:n) * Basis(1:n) )

           GradTemp = 0.0d0
           DO i = 1, dim
              GradTemp(i) = SUM( dBasisdx(1:n,i) * Temperature(1:n) )
           END DO
           SorD = SUM( Basis(1:n) * SoretD(1:n) )

           DO p=1,NBasis
              IF ( ConvectOrNPAndStabilize ) THEN
                 Load = Tau * SW(p)
              ELSE
                 Load = 1.0d0
              END IF

              SoretForce = CThermal * SorD * SUM( GradTemp(1:dim) * dBasisdx(p,1:dim) )
              ForceVector(p) = ForceVector(p) - s * SoretForce * Load
           END DO

        END IF

     END DO

   END SUBROUTINE DiffuseConvectiveCompose

   SUBROUTINE DiffuseConvectiveBBoundary( BoundaryMatrix, Parent, &
               pn, ParentNodes, Ratio, Element, n, Nodes )
     TYPE(Nodes_t) :: Nodes, ParentNodes
     TYPE(Element_t), TARGET :: Element, Parent
     REAL(KIND=dp) :: BoundaryMatrix(:,:), Ratio
     INTEGER :: n, pn

     REAL(KIND=dp) :: ddBasisddx(n,3,3), ParentdBasisdx(pn,3)
     REAL(KIND=dp) :: Basis(n), ParentBasis(pn)
     REAL(KIND=dp) :: dBasisdx(n,3), SqrtElementMetric

     REAL(KIND=dp) :: Diff(3,3)
     REAL(KIND=dp) :: u, v, w, s, x(n), y(n), z(n), Normal(3), FluxVector(3)
     REAL(KIND=dp), POINTER :: U_Integ(:), V_Integ(:), W_Integ(:), S_Integ(:)

     INTEGER :: ParentNodeIndexes(n)
     INTEGER :: i, t, q, p, N_Integ, j

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: stat

     BoundaryMatrix = 0.0D0
     IntegStuff = GaussPoints( Element )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

     DO t=1,N_Integ
       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)
       stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
           Basis,dBasisdx,ddBasisddx,.FALSE. )

       s = SqrtElementMetric * S_Integ(t)

       Normal = Normalvector( Element, Nodes, u, v, .TRUE. )
       DO i = 1,n
         DO j = 1,pn
           IF ( Element % NodeIndexes(i) == Parent % NodeIndexes(j) ) THEN
             x(i) = Parent % TYPE % NodeU(j)
             y(i) = Parent % TYPE % NodeV(j)
             z(i) = Parent % TYPE % NodeW(j)
             ParentNodeIndexes(i) = j
             EXIT
           END IF
         END DO
       END DO

       u = SUM( Basis(1:n) * x(1:n) )
       v = SUM( Basis(1:n) * y(1:n) )
       w = SUM( Basis(1:n) * z(1:n) )

       stat = ElementInfo( Parent, ParentNodes,u, v, w, SqrtElementMetric, &
           ParentBasis, ParentdBasisdx, ddBasisddx, .FALSE. )

       FluxVector = 0.0d0

       DO i = 1, 3
         DO j = 1, 3
           Diff(i,j) = SUM( Diffusivity(i,j,1:n) * Basis(1:n) )
         END DO
       END DO

       DO q = 1, pn
         DO j = 1, 3
           FluxVector(j) = SUM( Diff(j,1:3) * ParentdBasisdx(q,1:3) )
         END DO
         DO i = 1, n
           p = ParentNodeIndexes(i)
           BoundaryMatrix(p,q) = BoundaryMatrix(p,q) + Ratio * &
               s * Basis(i) * SUM( FluxVector(1:3) * Normal(1:3) )
         END DO
       END DO

     END DO

   END SUBROUTINE DiffuseConvectiveBBoundary

   SUBROUTINE DiffuseConvectiveBoundary( BoundaryMatrix,BoundaryVector, &
               LoadVector,NodalAlpha,Element,n,Nodes )
     REAL(KIND=dp) :: BoundaryMatrix(:,:),BoundaryVector(:), &
                    LoadVector(:),NodalAlpha(:)

     TYPE(Nodes_t)   :: Nodes
     TYPE(Element_t) :: Element

     INTEGER :: n

     REAL(KIND=dp) :: ddBasisddx(n,3,3)
     REAL(KIND=dp) :: Basis(n)
     REAL(KIND=dp) :: dBasisdx(n,3),SqrtElementMetric

     REAL(KIND=dp) :: u,v,w,s
     REAL(KIND=dp) :: Force,Alpha
     REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)

     INTEGER :: i,t,q,p,N_Integ

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: stat

     BoundaryVector = 0.0D0
     BoundaryMatrix = 0.0D0
     IntegStuff = GaussPoints( Element )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

     DO t=1,N_Integ
       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)
       stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
                  Basis,dBasisdx,ddBasisddx,.FALSE. )

       s = SqrtElementMetric * S_Integ(t)
       Force = SUM( LoadVector(1:n)*Basis )
       Alpha = SUM( NodalAlpha(1:n)*Basis )

       DO p=1,N
         DO q=1,N
           BoundaryMatrix(p,q) = BoundaryMatrix(p,q) + &
              s * Alpha * Basis(q) * Basis(p)
         END DO
       END DO

       DO q=1,N
         BoundaryVector(q) = BoundaryVector(q) + s * Basis(q) * Force
       END DO
     END DO
   END SUBROUTINE DiffuseConvectiveBoundary

   SUBROUTINE DiffuseConvectiveGenCompose( MassMatrix,StiffMatrix,ForceVector, &
       LoadVector,NodalCT,NodalC0,NodalC1,NodalC2,Temperature,EPotential, &
         Ux,Uy,Uz,MUx,MUy, MUz,SoretD,NEC, Compressible,AbsoluteMass, &
             Stabilize,Element,n,Nodes )
     REAL(KIND=dp), DIMENSION(:) :: ForceVector,Ux,Uy,Uz,MUx,MUy,MUz,LoadVector
     REAL(KIND=dp), DIMENSION(:,:) :: MassMatrix,StiffMatrix
     REAL(KIND=dp) :: NodalC0(:),NodalC1(:),NodalCT(:),NodalC2(:,:,:)
     REAL(KIND=dp) :: Temperature(:), SoretD(:), dT
     REAL(KIND=dp) :: EPotential(:), NEC(:), dPhi

     LOGICAL :: Stabilize,Compressible,AbsoluteMass

     INTEGER :: n

     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t), TARGET :: Element

     REAL(KIND=dp) :: ddBasisddx(n,3,3)
     REAL(KIND=dp) :: Basis(2*n)
     REAL(KIND=dp) :: dBasisdx(2*n,3),SqrtElementMetric

     REAL(KIND=dp) :: Velo(3),Force

     REAL(KIND=dp) :: A,M
     REAL(KIND=dp) :: Load

     REAL(KIND=dp) :: VNorm,hK,mK
     REAL(KIND=dp) :: Lambda=1.0,Pe,Pe1,Pe2,C00,Tau,Delta,x,y,z

     REAL(KIND=dp) :: SorD, GradTemp(3), SoretForce
     REAL(KIND=dp) :: GradPhi(3), NEConst

     INTEGER :: i,j,k,c,p,q,t,dim,N_Integ,NBasis

     REAL(KIND=dp) :: s,u,v,w,DivVelo,dVelodx(3,3)

     REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3)

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     REAL(KIND=dp) :: C0,CT,C1,C2(3,3),dC2dx(3,3,3),SU(n),SW(n)
     REAL(KIND=dp) :: NodalCThermal(n), CThermal
     REAL(KIND=dp) :: NodalCNP(n), CNP

     LOGICAL :: stat,CylindricSymmetry,Convection,ConvectOrNPAndStabilize,Bubbles
     LOGICAL :: ThermalDiffusion,NPDiffusion

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     CylindricSymmetry = (CurrentCoordinateSystem() == CylindricSymmetric .OR. &
                  CurrentCoordinateSystem() == AxisSymmetric)

     IF ( CylindricSymmetry ) THEN
       dim = 3
     ELSE
       dim = CoordinateSystemDimension()
     END IF
     n = element % TYPE % NumberOfNodes

     ForceVector = 0.0D0
     StiffMatrix = 0.0D0
     MassMatrix  = 0.0D0
     Load = 0.0D0

     Convection =  ANY( NodalC1 /= 0.0d0 )
     NBasis = n
     Bubbles = .FALSE.
     IF ( Convection .AND. .NOT. Stabilize ) THEN
        NBasis = 2*n
        Bubbles = .TRUE.
     END IF

     ThermalDiffusion = .FALSE.
     IF ( ANY( ABS( SoretD(1:n) ) > AEPS ) ) THEN
        ThermalDiffusion = .TRUE.
        NodalCThermal = NodalC0(1:n)
     END IF
     NPDiffusion = .FALSE.
     IF ( ANY( ABS( NEC(1:n) ) > AEPS ) ) THEN
        NPDiffusion = .TRUE.
        NodalCNP = NodalC0(1:n)
     END IF
     NodalC0 = 0.0d0

     IF ( Bubbles ) THEN
        IntegStuff = GaussPoints( element, Element % TYPE % GaussPoints2 )
     ELSE
        IntegStuff = GaussPoints( element )
     END IF
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

     ConvectOrNPAndStabilize = .FALSE.
     IF ( Stabilize .AND. (ANY(NodalC1 /= 0.0D0) .OR. NPDiffusion) ) THEN
       ConvectOrNPAndStabilize = .TRUE.
       hK = element % hK
       mK = element % StabilizationMK
     END IF

     DO t=1,N_Integ

       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)
      stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
             Basis,dBasisdx,ddBasisddx,ConvectOrNPAndStabilize,Bubbles )

       IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
         x = SUM( nodes % x(1:n)*Basis(1:n) )
         y = SUM( nodes % y(1:n)*Basis(1:n) )
         z = SUM( nodes % z(1:n)*Basis(1:n) )
       END IF

       CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )

       s = SqrtMetric * SqrtElementMetric * S_Integ(t)
       C0 = SUM( NodalC0(1:n)*Basis(1:n) )
       CT = SUM( NodalCT(1:n)*Basis(1:n) )
       C1 = SUM( NodalC1(1:n)*Basis(1:n) )

       DO i=1,dim
         DO j=1,dim
           C2(i,j) = SQRT(Metric(i,i)) * SQRT(Metric(j,j)) * &
                SUM( NodalC2(i,j,1:n) * Basis(1:n) )
         END DO
       END DO

       Convection = .FALSE.
       IF ( C1 /= 0.0D0 ) THEN
         Convection = .TRUE.
         Velo = 0.0D0
         Velo(1) = SUM( (Ux(1:n)-MUx(1:n))*Basis(1:n) )
         Velo(2) = SUM( (Uy(1:n)-MUy(1:n))*Basis(1:n) )
         IF ( dim > 2 .AND. CurrentCoordinateSystem() /= AxisSymmetric ) THEN
           Velo(3) = SUM( (Uz(1:n)-MUz(1:n))*Basis(1:n) )
         END IF

         IF ( Compressible .AND. AbsoluteMass ) THEN

           dVelodx = 0.0D0
           DO i=1,3
             dVelodx(1,i) = SUM( Ux(1:n)*dBasisdx(1:n,i) )
             dVelodx(2,i) = SUM( Uy(1:n)*dBasisdx(1:n,i) )
             IF ( dim > 2 .AND. CurrentCoordinateSystem() /= AxisSymmetric ) &
               dVelodx(3,i) = SUM( Uz(1:n)*dBasisdx(1:n,i) )
           END DO

           DivVelo = 0.0D0
           DO i=1,dim
             DivVelo = DivVelo + dVelodx(i,i)
           END DO
           IF ( CurrentCoordinateSystem() >= Cylindric .AND. &
                CurrentCoordinateSystem() <= AxisSymmetric ) THEN
             DivVelo = DivVelo + Velo(1)/x
           ELSE
             DO i=1,dim
               DO j=i,dim
                 DivVelo = DivVelo + Velo(j)*Symb(i,j,i)
               END DO
             END DO
           END IF
           C0 = DivVelo
         END IF

          IF ( NPDiffusion ) THEN
             NEConst = SUM ( NEC(1:n)*Basis(1:n) )
             GradPhi = 0.0d0
             DO i = 1, dim
                GradPhi(i) = SUM( dBasisdx(1:n,i) * EPotential(1:n) )
             END DO
             DO i = 1, dim
                DO j = 1, dim
                   Velo(i) = Velo(i) - C2(i,j) * NEConst * GradPhi(j)
                END DO
             END DO
          END IF

         IF ( Stabilize ) THEN

           Vnorm = 0.0D0
           DO i=1,dim
              Vnorm = Vnorm + Velo(i)*Velo(i) / Metric(i,i)
           END DO
           Vnorm = SQRT( Vnorm )

           Pe = MIN(1.0D0,mK*hK*C1*VNorm/(2*ABS(C2(1,1))))

           Tau = 0.0D0
           IF ( VNorm /= 0.0D0 ) THEN
             Tau = hK * Pe / (2 * C1 * VNorm)
           END IF

           DO i=1,dim
             DO j=1,dim
               DO k=1,3
                 dC2dx(i,j,k) = SQRT(Metric(i,i))*SQRT(Metric(j,j))* &
                      SUM(NodalC2(i,j,1:n)*dBasisdx(1:n,k))
               END DO
             END DO
           END DO
           DO p=1,n
             SU(p) = C0 * Basis(p)
             DO i = 1,dim
               SU(p) = SU(p) + C1 * dBasisdx(p,i) * Velo(i)
               IF ( Element % TYPE % BasisFunctionDegree <= 1 ) CYCLE

               DO j=1,dim
                 SU(p) = SU(p) - C2(i,j) * ddBasisddx(p,i,j)
                 SU(p) = SU(p) - dC2dx(i,j,j) * dBasisdx(p,i)
                 DO k=1,dim
                   SU(p) = SU(p) + C2(i,j) * Symb(i,j,k) * dBasisdx(p,k)
                   SU(p) = SU(p) - C2(i,k) * Symb(k,j,j) * dBasisdx(p,i)
                   SU(p) = SU(p) - C2(k,j) * Symb(k,j,i) * dBasisdx(p,i)
                 END DO
               END DO
             END DO

             SW(p) = C0 * Basis(p)

             DO i = 1,dim
               SW(p) = SW(p) + C1 * dBasisdx(p,i) * Velo(i)
               IF ( Element % TYPE % BasisFunctionDegree <= 1 ) CYCLE

               DO j=1,dim
                 SW(p) = SW(p) - C2(i,j) * ddBasisddx(p,i,j)
                 SW(p) = SW(p) - dC2dx(i,j,j) * dBasisdx(p,i)
                 DO k=1,dim
                   SW(p) = SW(p) + C2(i,j) * Symb(i,j,k) * dBasisdx(p,k)
                   SW(p) = SW(p) - C2(i,k) * Symb(k,j,j) * dBasisdx(p,i)
                   SW(p) = SW(p) - C2(k,j) * Symb(k,j,i) * dBasisdx(p,i)
                 END DO
               END DO
             END DO
           END DO
         END IF
       END IF
       DO p=1,NBasis
       DO q=1,NBasis
         M = CT * Basis(q) * Basis(p)
         A = C0 * Basis(q) * Basis(p)
         DO i=1,dim
           DO j=1,dim
             A = A + C2(i,j) * dBasisdx(q,i) * dBasisdx(p,j)
           END DO
         END DO

         IF ( Convection ) THEN
           DO i=1,dim
             A = A + C1 * Velo(i) * dBasisdx(q,i) * Basis(p)
           END DO

           IF ( Stabilize ) THEN
             A = A + Tau * SU(q) * SW(p)
             M = M + Tau * CT * Basis(q) * SW(p)
           END IF
         END IF

         StiffMatrix(p,q) = StiffMatrix(p,q) + s * A
         MassMatrix(p,q)  = MassMatrix(p,q)  + s * M
       END DO
       END DO

       Force = SUM( LoadVector(1:n)*Basis(1:n) )

       DO p=1,NBasis
         Load = Basis(p)

         IF ( ConvectOrNPAndStabilize ) THEN
           Load = Load + Tau * SW(p)
         END IF

         ForceVector(p) = ForceVector(p) + s * Load * Force
       END DO

        IF ( ThermalDiffusion ) THEN

           CThermal = SUM( NodalCThermal(1:n) * Basis(1:n) )

           GradTemp = 0.0d0
           IF ( CurrentCoordinateSystem() >= Cylindric .AND. &
                CurrentCoordinateSystem() <= AxisSymmetric ) THEN

             DO i = 1, dim
               GradTemp(i) = SUM( dBasisdx(1:n,i) * Temperature(1:n) )
             END DO
           ELSE
             CALL Error( 'AdvectionDiffusion', &
                 'Thermal diffusion not implemented for this coordinate system. Ignoring it.' )
           END IF

           SorD = SUM( Basis(1:n) * SoretD(1:n) )

           DO p=1,NBasis
              IF ( ConvectOrNPAndStabilize ) THEN
                 Load = Tau * SW(p)
              ELSE
                 Load = 1.0d0
              END IF

              SoretForce = CThermal * SorD * SUM( GradTemp(1:dim) * dBasisdx(p,1:dim) )
              ForceVector(p) = ForceVector(p) - s * SoretForce * Load
           END DO

        END IF

     END DO
   END SUBROUTINE DiffuseConvectiveGenCompose

   SUBROUTINE DiffuseConvectiveGenBBoundary( BoundaryMatrix, Parent, &
               pn, ParentNodes, Ratio, Element, n, Nodes )
     TYPE(Nodes_t)   :: Nodes, ParentNodes
     TYPE(Element_t), TARGET :: Element, Parent
     REAL(KIND=dp) :: BoundaryMatrix(:,:), Ratio
     INTEGER :: n, pn

     REAL(KIND=dp) :: ddBasisddx(n,3,3), ParentdBasisdx(pn,3)
     REAL(KIND=dp) :: Basis(n), ParentBasis(pn)
     REAL(KIND=dp) :: dBasisdx(n,3), SqrtElementMetric

     REAL(KIND=dp) :: Diff(3,3)
     REAL(KIND=dp) :: xpos, ypos, zpos
     REAL(KIND=dp) :: u, v, w, s, x(n), y(n), z(n), Normal(3), FluxVector(3)
     REAL(KIND=dp), POINTER :: U_Integ(:), V_Integ(:), W_Integ(:), S_Integ(:)

     INTEGER :: ParentNodeIndexes(n)
     INTEGER :: i, t, q, p, N_Integ, j

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: stat

     BoundaryMatrix = 0.0D0
     IntegStuff = GaussPoints( Element )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

     DO t=1,N_Integ
       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)
       stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
           Basis,dBasisdx,ddBasisddx,.FALSE. )

       s = SqrtElementMetric * S_Integ(t)
       IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
         xpos = SUM( Nodes % x(1:n)*Basis )
         ypos = SUM( Nodes % y(1:n)*Basis )
         zpos = SUM( Nodes % z(1:n)*Basis )
         s = s * CoordinateSqrtMetric( xpos,ypos,zpos )
       END IF

       Normal = Normalvector( Element, Nodes, u, v, .TRUE. )
       DO i = 1,n
         DO j = 1,pn
           IF ( Element % NodeIndexes(i) == Parent % NodeIndexes(j) ) THEN
             x(i) = Parent % TYPE % NodeU(j)
             y(i) = Parent % TYPE % NodeV(j)
             z(i) = Parent % TYPE % NodeW(j)
             ParentNodeIndexes(i) = j
             EXIT
           END IF
         END DO
       END DO

       u = SUM( Basis(1:n) * x(1:n) )
       v = SUM( Basis(1:n) * y(1:n) )
       w = SUM( Basis(1:n) * z(1:n) )

       stat = ElementInfo( Parent, ParentNodes,u, v, w, SqrtElementMetric, &
           ParentBasis, ParentdBasisdx, ddBasisddx, .FALSE. )

       FluxVector = 0.0d0
       DO i = 1, 3
         DO j = 1, 3
           Diff(i,j) = SUM( Diffusivity(i,j,1:n) * Basis(1:n) )
         END DO
       END DO

       DO q = 1, pn
         DO j = 1, 3
           FluxVector(j) = SUM( Diff(j,1:3) * ParentdBasisdx(q,1:3) )
         END DO
         DO i = 1, n
           p = ParentNodeIndexes(i)
           BoundaryMatrix(p,q) = BoundaryMatrix(p,q) + Ratio * &
               s * Basis(i) * SUM( FluxVector(1:3) * Normal(1:3) )
         END DO
       END DO

     END DO

   END SUBROUTINE DiffuseConvectiveGenBBoundary

   SUBROUTINE DiffuseConvectiveGenBoundary( BoundaryMatrix,BoundaryVector, &
              LoadVector,NodalAlpha,Element,n,Nodes)
     REAL(KIND=dp) :: BoundaryMatrix(:,:),BoundaryVector(:)
     REAL(KIND=dp) :: LoadVector(:),NodalAlpha(:)
     TYPE(Nodes_t)    :: Nodes
     TYPE(Element_t), TARGET  :: Element

     INTEGER :: n

     REAL(KIND=dp) :: ddBasisddx(n,3,3)
     REAL(KIND=dp) :: Basis(n)
     REAL(KIND=dp) :: dBasisdx(n,3),SqrtElementMetric

     REAL(KIND=dp) :: u,v,w,s,x,y,z
     REAL(KIND=dp) :: Force,Alpha
     REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)

     REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3)

     INTEGER :: i,t,q,p,N_Integ

     LOGICAL :: stat,CylindricSymmetry

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     BoundaryVector = 0.0D0
     BoundaryMatrix = 0.0D0

     IntegStuff = GaussPoints( element )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

     DO t=1,N_Integ
       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)

       stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
                  Basis,dBasisdx,ddBasisddx,.FALSE. )

       s =  S_Integ(t) * SqrtElementMetric
       IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
         x = SUM( Nodes % x(1:n)*Basis )
         y = SUM( Nodes % y(1:n)*Basis )
         z = SUM( Nodes % z(1:n)*Basis )
         s = s * CoordinateSqrtMetric( x,y,z )
       END IF
       Alpha = SUM( NodalAlpha(1:n)*Basis )
       Force = SUM( LoadVector(1:n)*Basis )

       DO p=1,N
         DO q=1,N
           BoundaryMatrix(p,q) = BoundaryMatrix(p,q) + &
               s * Alpha * Basis(q) * Basis(p)
         END DO
       END DO

       DO q=1,N
         BoundaryVector(q) = BoundaryVector(q) + s * Basis(q) * Force
       END DO
     END DO
  END SUBROUTINE DiffuseConvectiveGenBoundary

end module advectiondiffusion_acc
