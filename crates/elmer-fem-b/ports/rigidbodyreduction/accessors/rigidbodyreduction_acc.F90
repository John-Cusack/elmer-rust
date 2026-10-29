! T021-B RigidBodyReduction port: section-level BIND(C) wrappers for the
! Before-Linsolve entry `rigidbody_` in fem/src/modules/RigidBodyReduction.F90,
! compiled against the built .mod files of the target configure (A2 shipped
! mechanism, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(KIND=dp) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results — stays
! inside Fortran, so the Rust side is pure C calls plus the INTEGER result.
!
! The entry body is split into 11 ordered sections; the Rust entry calls
! them in the original statement order with the original branch structure
! (the only Rust-visible branch is the NbrRigids==0 early return). Each
! section body is a verbatim copy of the corresponding original lines, with
! only these documented adaptations:
!  1. Entry dummies (Model, Solver, Matrix POINTER, b/x arrays, n, DOFs,
!     Norm) arrive as C handles and are converted with C_F_POINTER to
!     locals named exactly as the original dummies (model/solver/A/b/x/
!     n/Norm), so section bodies read verbatim. The Matrix POINTER dummy
!     arrives double-indirected (address of the pointer variable — the
!     gfortran convention under the explicit POINTER interface used by
!     LoadMod's execlinsolveprocs, proven by probe); the Rust entry
!     dereferences it once and passes the Matrix_t address as h_a.
!  2. The original SAVE list (C, C_Trans, RedA, AllocationsDone, u,
!     RealNodeTypes, Permutation, CenterOfRigidBody, NbrRigids, RigidIndex,
!     f, NodeTypes, Verbose) lives as module SAVE state with identical
!     lifetimes. Cross-section non-SAVE locals (Perm, BMatrix, FVector,
!     XVector, LocalEigenVectors, TotTime, RowsInUnityBlock, NEigen, DOFs,
!     istat) live as module state too; every one is written before it is
!     read on every path that reaches the read, exactly as in the
!     original. (RowsInUnityBlock is read by the rigid-C fill on
!     MeshChanged re-entry without recompute — there it keeps the
!     first-call value; the original leaves it undefined on that path.
!     The proof never takes MeshChanged, so the proven paths are exact.)
!  3. The CONTAINS internal procedures are copied verbatim as module
!     procedures with two one-word adaptations for lost host
!     association: ComputeMassCenter gains a local `LOGICAL :: Rigid`
!     (host scratch, always written before read), and
!     MatrixBandWidthOptimize reads the calling solver through the
!     module pointer g_solver (set by rbr_optimize before the call).
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text beyond the transcribed entry body.
module rbr_acc
  use iso_c_binding
  use types
  use lists
  use solverbasics
  use crsmatrix
  use generalutils
  use integration
  use defutils
  use messages
  implicit none

  ! -- mirror of the RigidBody SAVE list (identical lifetimes) --
  type(matrix_t), pointer, save :: C => null()
  type(matrix_t), pointer, save :: C_Trans => null()
  type(matrix_t), pointer, save :: RedA => null()
  logical, save :: AllocationsDone = .false.
  real(kind=dp), pointer, save :: u(:) => null()
  integer, allocatable, save :: RealNodeTypes(:)
  integer, pointer, save :: Permutation(:) => null()
  real(kind=dp), pointer, save :: CenterOfRigidBody(:,:) => null()
  integer, save :: NbrRigids = 0
  integer, allocatable, save :: RigidIndex(:)
  real(kind=dp), pointer CONTIG, save :: f(:) => null()
  real(kind=dp), pointer, save :: NodeTypes(:) => null()
  logical, save :: Verbose = .false.
  ! -- first-call scalars read across sections (non-SAVE in the original,
  ! -- set in section 2 and read in section 5 within the first call only)
  integer, save :: RigidNodes = 0
  integer, save :: Dirichlet = 0

  ! -- cross-section locals (written before read on every taken path) --
  type(solver_t), pointer, save :: g_solver => null()
  type(matrix_t), pointer, save :: BMatrix => null()
  integer, pointer, save :: Perm(:) => null()
  real(kind=dp), pointer, save :: FVector(:) => null()
  real(kind=dp), pointer, save :: XVector(:) => null()
  real(kind=dp), pointer, save :: LocalEigenVectors(:,:) => null()
  integer, save :: DOFs = 0
  integer, save :: istat = 0
  integer, save :: NEigen = 0
  integer, save :: RowsInUnityBlock = 0
  real(kind=dp), save :: TotTime = 0.0_dp

contains

  ! -- section 1: rigid-body scan (original lines 114-140) -------------
  subroutine rbr_scan(h_model, h_solver, nrig) bind(c, name="rbr_scan")
    type(c_ptr), value :: h_model, h_solver
    integer(c_int) :: nrig
    type(model_t), pointer :: model
    type(solver_t), pointer :: solver
    integer :: t
    logical :: Rigid, GotIt
    call c_f_pointer(h_model, model)
    call c_f_pointer(h_solver, solver)

    TotTime = CPUtime()

    Verbose = .FALSE.

    IF ( .NOT. AllocationsDone ) THEN

       ALLOCATE( RigidIndex( Model % NumberOfBodies ) )
       RigidIndex = 0
       NbrRigids = 0
       DO t = 1, Model % NumberOfBodies
          Rigid = ListGetLogical( Model % Bodies(t) % Values, 'Rigid Body', gotIt )
          IF ( Rigid ) THEN
             NbrRigids = NbrRigids + 1
             RigidIndex(NbrRigids) = t
          END IF
       END DO

       IF ( NbrRigids == 0 ) THEN
          nrig = 0
          CALL Warn( 'RigidBody', 'No Rigid bodies defined.' )
          RETURN
       END IF

       WRITE( Message, * ) 'Found ', NbrRigids, ' rigid blocks'
       CALL Info( 'RigidBody', Message, Level=5 )

       Verbose = ListGetLogical( Solver % Values, 'Additional Info', GotIt )
    END IF
    nrig = NbrRigids
  end subroutine rbr_scan

  ! -- section 2: node classification + Dirichlet elimination --
  ! -- (original lines 141-497) ----------------------------------------
  subroutine rbr_setup(h_model, h_solver, h_a, b, n, ndofs) bind(c, name="rbr_setup")
    type(c_ptr), value :: h_model, h_solver, h_a
    integer(c_int), value :: n, ndofs
    real(kind=dp) :: b(n)
    type(model_t), pointer :: Model
    type(solver_t), pointer :: Solver
    type(matrix_t), pointer :: A
    type(element_t), pointer :: CurrentElement
    type(solver_t), pointer :: DummySolver
    integer, pointer :: NodeIndexes(:)
    integer, pointer :: DCount(:)
    integer, allocatable :: FixedIndex(:)
    logical, allocatable :: Visited(:)
    real(kind=dp) :: s, val
    logical :: GotIt, Rigid
    logical :: Mass, Damp
    ! (RigidNodes, Dirichlet bind to module state: set here on the first
    ! call, read by section 5 within the same first call — identical values)
    integer :: NbrFixed, LBody, RBody
    integer :: i, j, k, m, t, p, l
    INTRINSIC PACK, MERGE
    call c_f_pointer(h_model, Model)
    call c_f_pointer(h_solver, Solver)
    call c_f_pointer(h_a, A)
    DOFs = ndofs

    IF ( .NOT. AllocationsDone ) THEN

       Perm => Solver % Variable % Perm

       ALLOCATE( Visited( Model % NumberOfNodes ), &
            RealNodeTypes( Model % NumberOfNodes ), STAT=istat )

       IF ( istat /= 0 )  CALL Fatal( 'RigidBody', 'Memory allocation error 1' )
       Visited = .FALSE.
       RealNodeTypes = 0

       DO t = 1, Solver % NumberOfActiveElements

          CurrentElement => Solver % Mesh % Elements( Solver % ActiveElements(t) )
          NodeIndexes => CurrentElement % NodeIndexes
          m = CurrentElement % TYPE % NumberOfNodes

          IF ( .NOT. ALL( Visited( NodeIndexes(1:m) ) ) ) THEN
             DO i = 1, NbrRigids
                IF ( RigidIndex(i) == CurrentElement % BodyId ) THEN
                   RealNodeTypes( NodeIndexes(1:m) ) = i
                   EXIT
                END IF
             END DO
             Visited( NodeIndexes(1:m) ) = .TRUE.
          END IF
       END DO

       IF ( ListGetLogical( Solver % Values, 'Extend Elastic Region', GotIt )) THEN

          p = ListGetInteger( Solver % Values, 'Extend Elastic Layers', GotIt )
          IF ( .NOT. GotIt )  p = 1

          IF ( Verbose ) THEN
             RigidNodes = 0
             DO t = 1, NbrRigids
                RigidNodes = Rigidnodes + COUNT( RealNodeTypes == t )
                WRITE ( Message, * ) 'Rigid nodes in body ', RigidIndex(t), &
                     ': ', COUNT( RealNodeTypes == t )
                CALL Info( 'RigidBody', Message, Level=24 )
             END DO
          END IF
          WRITE ( Message, * ) 'Extending elastic regions by  ', p , &
               ' layers'
          CALL Info( 'RigidBody', Message, Level=16 )

          DO i = 1, p

             Visited = .FALSE.
             DO t = 1, Solver % NumberOfActiveElements

                CurrentElement => &
                     Solver % Mesh % Elements( Solver % ActiveElements(t) )
                NodeIndexes => CurrentElement % NodeIndexes
                m = CurrentElement % TYPE % NumberOfNodes

                IF ( .NOT. ALL( Visited( NodeIndexes(1:m) ) ) ) THEN

                   IF ( ANY( RealNodeTypes( NodeIndexes(1:m) ) == 0 ) .AND. &
                        ANY( RealNodeTypes( NodeIndexes(1:m) ) > 0 ) ) THEN

                      RealNodeTypes( NodeIndexes(1:m) ) = &
                           MERGE( -1, RealNodeTypes( NodeIndexes(1:m) ), &
                           RealNodeTypes( NodeIndexes(1:m) ) > 0 )

                   END IF
                   Visited( NodeIndexes(1:m) ) = .TRUE.
                END IF
             END DO

             RealNodeTypes = MERGE( 0, RealNodeTypes, RealNodeTypes < 0 )
          END DO

          IF ( Verbose ) THEN
             RigidNodes = 0
             DO t = 1, NbrRigids
                RigidNodes = Rigidnodes + COUNT( RealNodeTypes == t )
                WRITE ( Message, * ) 'Rigid nodes in body ', RigidIndex(t), &
                     ': ', COUNT( RealNodeTypes == t )
                CALL Info( 'RigidBody', Message, Level=24 )
             END DO
          END IF

       END IF

       DEALLOCATE( Visited )

       ALLOCATE( Visited( Model % NumberOfBCs ) )
       ALLOCATE( FixedIndex( NbrRigids ) )

       Visited = .FALSE.
       FixedIndex = 0

       NbrFixed = 0

       DO t = Model % NumberOfBulkElements + 1, Model % NumberOfBulkElements &
            + Model % NumberOfBoundaryElements

          CurrentElement => Solver % Mesh % Elements(t)

          IF ( CurrentElement % BoundaryInfo % Constraint == 0 )  CYCLE
          IF ( Visited( CurrentElement % BoundaryInfo % Constraint ) )  CYCLE
          Visited( CurrentElement % BoundaryInfo % Constraint ) = .TRUE.

          NodeIndexes => CurrentElement % NodeIndexes
          m = CurrentElement % TYPE % NumberOfNodes

          DO i = 1, Model % NumberOfBCs
             IF ( CurrentElement % BoundaryInfo % Constraint == &
                  Model % BCs(i) % Tag ) THEN

                GotIt = ListCheckPresent( Model % BCs(i) % Values, &
                     ComponentName(Solver % Variable,1))
                IF( .NOT. GotIt ) &
                   GotIt = ListCheckPresent( Model % BCs(i) % Values, &
                     ComponentName(Solver % Variable,2) )
                IF ( .NOT. GotIt )  &
                   GotIt = ListCheckPresent( Model % BCs(i) % Values, &
                     ComponentName(Solver % Variable,3) )
                IF ( .NOT. GotIt )  CYCLE


                Lbody = -1
                IF ( ASSOCIATED( CurrentElement % BoundaryInfo % Left ) ) &
                   LBody = CurrentElement % BoundaryInfo % Left % BodyID

                Rbody = -1
                IF ( ASSOCIATED( CurrentElement % BoundaryInfo % Right ) ) &
                   RBody = CurrentElement % BoundaryInfo % Right % BodyID


                IF ( ANY( RigidIndex == LBody ) .OR. ANY( RigidIndex == RBody )) THEN
                   NbrRigids = NbrRigids - 1
                   NbrFixed = NbrFixed + 1
                   DO j = 1, NbrRigids + 1
                      IF ( RigidIndex(j) == LBody .OR. &
                           RigidIndex(j) == RBody ) THEN
                         FixedIndex( NbrFixed ) = RigidIndex(j)
                         RigidIndex(j) = 0
                      END IF
                   END DO
                   RigidIndex = PACK( RigidIndex, RigidIndex > 0 )
                   RigidIndex( NbrRigids + 1: Model % NumberOfBodies ) = 0
                END IF

             END IF
          END DO

          IF ( NbrRigids == 0 )  EXIT


       END DO

       DEALLOCATE( Visited )

       IF ( NbrFixed > 0 ) THEN

          IF ( NbrFixed == 1 ) THEN
             WRITE( Message, * ) 'Setting ', NbrFixed, &
                  ' rigid block as fixed (body ', FixedIndex(1), ')'
          ELSE
             WRITE( Message, * ) 'Setting ', NbrFixed, &
                  ' rigid blocks as fixed (bodies ', FixedIndex(1:NbrFixed), ')'
          END IF
          CALL Info( 'RigidBody', Message, Level=8 )

          DO t = 1, Solver % NumberOfActiveElements

             CurrentElement => &
                  Solver % Mesh % Elements( Solver % ActiveElements(t))

             IF ( ANY( FixedIndex == CurrentElement % BodyId ) ) THEN

                m = CurrentElement % TYPE % NumberOfNodes
                NodeIndexes => CurrentElement % NodeIndexes(1:m)
                val = 0.0_dp

                DO j = 1, m

                   IF ( RealNodeTypes( NodeIndexes( j ) ) == 0 )  CYCLE

                   RealNodeTypes( NodeIndexes( j ) ) = -1
                   k = Perm(NodeIndexes(j))

                   IF ( k > 0 ) THEN
                      DO i = 1, DOFs
                        CALL SetDirichletPoint( A, b, i, DOFs, Perm, NodeIndexes(j), val )
                      END DO
                   END IF

                END DO
             END IF
          END DO
       END IF


       DEALLOCATE( FixedIndex )


       IF ( ListGetLogical( Solver % Values, 'Output Node Types', GotIt ) ) THEN
          ALLOCATE( NodeTypes( Model % NumberOfNodes ) )
          NodeTypes = 0.0d0

          NodeTypes = REAL( RealNodeTypes )
          DummySolver => Solver
          CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, &
               DummySolver, 'Node Type', 1, NodeTypes )
       END IF

       BMatrix => A % EMatrix

       Mass = .FALSE.
       IF ( ASSOCIATED( A % MassValues ) ) THEN
          IF ( SIZE( A % MassValues ) == SIZE( A % Values ) ) Mass = .TRUE.
       END IF

       Damp = .FALSE.
       IF ( ASSOCIATED( A % DampValues ) ) THEN
          IF ( SIZE( A % DampValues ) == SIZE( A % Values ) ) Damp = .TRUE.
       END IF

       IF ( .NOT. ASSOCIATED( BMatrix ) ) THEN

          ALLOCATE( Permutation( A % NumberOfRows ), STAT=istat )
          IF ( istat /= 0 )  &
               CALL Fatal( 'RigidBody', 'Memory allocation error 1.1' )

          Permutation = 0
          Dirichlet = 0

          l = 0
          DO i=1,A % NumberOFRows
             s = A % Values( A % Diag(i) )
             A % Values( A % Diag(i) ) = 0.0d0
             j = A % Rows(i)
             k = A % Rows(i+1)-1
             IF ( ALL( A % Values(j:k) == 0.0d0 ) ) THEN
                Dirichlet = Dirichlet + 1
             ELSE
                l = l + 1
                Permutation(i) = l
             END IF
             A % Values( A % Diag(i) ) = s
          END DO

          IF ( Verbose ) THEN
             WRITE( Message, * ) 'Will eliminate ', Dirichlet, ' Dirichlet DOFs'
             CALL Info( 'RigidBody', Message, Level=16 )
          END IF

          Bmatrix => AllocateMatrix()

          Bmatrix % NumberOFRows = A % NumberOFRows - Dirichlet
          ALLOCATE( Bmatrix % Rows( Bmatrix % NumberOfRows + 1 ), &
               Bmatrix % Diag( Bmatrix % NumberOfRows ) )

          j = 0
          k = 1
          DO i=1, A % NumberOFRows
             IF ( Permutation(i) /= 0 ) THEN
                j = j + 1
                Bmatrix % Rows(j) = k
                DO l = A % Rows(i),A % Rows(i+1)-1
                   t = A % Cols(l)
                   IF ( Permutation(t) /= 0 ) k = k + 1
                END DO
             END IF
          END DO
          Bmatrix % Rows(Bmatrix % NumberOFRows+1) = k

          IF ( Mass )  ALLOCATE( Bmatrix % MassValues(k-1) )
          IF ( Damp )  ALLOCATE( Bmatrix % DampValues(k-1) )

          ALLOCATE( Bmatrix % Values(k-1), Bmatrix % Cols(k-1), &
               BMatrix % RHS( Bmatrix % NumberOfRows ) )

          ALLOCATE( Dcount( A % NumberOfRows ) )
          IF ( Permutation(1) == 0 ) THEN
             DCount(1) = 1
          ELSE
             DCount(1) = 0
          END IF

          DO i=2, A % NumberOFRows
             IF ( Permutation(i) /= 0 ) THEN
                DCount(i) = DCount(i-1)
             ELSE
                DCount(i) = DCount(i-1) + 1
             END IF
          END DO

          j = 0
          k = 1
          DO i=1, A % NumberOFRows
             IF ( Permutation(i) /= 0 ) THEN
                j = j + 1
                DO l = A % Rows(i),A % Rows(i+1)-1
                   t = A % Cols(l)
                   IF ( Permutation( t ) /= 0 ) THEN
                      Bmatrix % Cols(k) = t - DCount(t)
                      k = k + 1
                   END IF
                END DO
             END IF
          END DO

          DEALLOCATE( DCount )

          DO i=1,Bmatrix % NumberOfRows
             DO j=Bmatrix % Rows(i), Bmatrix % Rows(i+1)-1
                IF ( Bmatrix % Cols(j) == i ) THEN
                   Bmatrix % Diag(i) = j
                END IF
             END DO
          END DO
          CALL CRS_SortMatrix( BMatrix )

          A % EMatrix => BMatrix
          A % Perm => Permutation
       END IF

       WRITE( Message, * ) 'Total nodes:              ', Model % NumberOfNodes
       CALL Info( 'RigidBody', Message, Level=8 )
       RigidNodes = 0
       DO t = 1, NbrRigids
          RigidNodes = Rigidnodes + COUNT( RealNodeTypes == t )
          WRITE( Message, * ) 'Rigid nodes in body ', RigidIndex(t), ': ', &
               COUNT( RealNodeTypes == t )
          CALL Info( 'RigidBody', Message, Level=8 )
       END DO

       ALLOCATE( CenterOfRigidBody( NbrRigids, 3 ) )
    END IF
  end subroutine rbr_setup

  ! -- section 3: Dirichlet-eliminated fill (original lines 499-537) ----
  subroutine rbr_fill(h_a, b, n) bind(c, name="rbr_fill")
    type(c_ptr), value :: h_a
    integer(c_int), value :: n
    real(kind=dp) :: b(n)
    type(matrix_t), pointer :: A
    integer :: j, k, l, t, i
    logical :: Mass, Damp
    call c_f_pointer(h_a, A)

    BMatrix => A % EMatrix

    Mass = .FALSE.
    IF ( ASSOCIATED( A % MassValues ) ) THEN
       IF ( SIZE( A % MassValues ) == SIZE( A % Values ) ) Mass = .TRUE.
    END IF

    Damp = .FALSE.
    IF ( ASSOCIATED( A % DampValues ) ) THEN
       IF ( SIZE( A % DampValues ) == SIZE( A % Values ) ) Damp = .TRUE.
    END IF

    Permutation => A % Perm
    FVector => BMatrix % RHS

    j = 0
    k = 1
    DO i=1, A % NumberOFRows
       IF ( Permutation(i) /= 0 ) THEN
          j = j + 1
          FVector(j) = B(i)
          DO l = A % Rows(i), A % Rows(i+1)-1
             t = A % Cols(l)
             IF ( Permutation(t) == 0 ) THEN
                FVector(j) = FVector(j) - &
                     A % Values(l) * b(t) / A % Values( A % Diag(t) )
             ELSE
                Bmatrix % Values(k) = A % Values(l)
                IF ( Mass )  Bmatrix % MassValues(k) = A % MassValues(l)
                IF ( Damp )  Bmatrix % DampValues(k) = A % DampValues(l)
                k = k + 1
             END IF
          END DO
       END IF
    END DO
  end subroutine rbr_fill

  ! -- section 4: rigid mass centers (original lines 544-549) -----------
  subroutine rbr_centers(h_model, h_solver) bind(c, name="rbr_centers")
    type(c_ptr), value :: h_model, h_solver
    type(model_t), pointer :: Model
    type(solver_t), pointer :: Solver
    call c_f_pointer(h_model, Model)
    call c_f_pointer(h_solver, Solver)

    IF ( .NOT. AllocationsDone .OR. Solver % MeshChanged ) THEN

       CALL ComputeMassCenter( CenterOfRigidBody, Model, Solver, &
            RigidIndex, NbrRigids )

    END IF
  end subroutine rbr_centers

  ! -- section 5: projection-matrix structure (original lines 555-626) --
  subroutine rbr_cstruct(h_model, h_solver, h_a, ndofs) bind(c, name="rbr_cstruct")
    type(c_ptr), value :: h_model, h_solver, h_a
    integer(c_int), value :: ndofs
    type(model_t), pointer :: Model
    type(solver_t), pointer :: Solver
    type(matrix_t), pointer :: A
    integer :: NumberOfRows, NumberOfValues
    integer :: i, j, k, t
    call c_f_pointer(h_model, Model)
    call c_f_pointer(h_solver, Solver)
    call c_f_pointer(h_a, A)
    DOFs = ndofs

    IF ( .NOT. AllocationsDone ) THEN

       NumberOfRows = DOFs *(Model % NumberOfNodes -RigidNodes) - Dirichlet &
            + NbrRigids * 3 * (DOFs -1)

       NumberOfValues = NumberOfRows - NbrRigids * 3 * (DOFs-1) + &
            DOFs*DOFs*( RigidNodes )

       C => AllocateMatrix()
       C % NumberOfRows = NumberOfRows
       ALLOCATE( C % Rows( NumberOfRows+1 ), C % Cols( NumberOfValues ), &
            C % Values( NumberOfValues ), C % Diag( NumberOfRows ), STAT=istat )
       IF ( istat /= 0 )  CALL Fatal( 'RigidBody', 'Memory allocation error 2' )

       ALLOCATE( u( NumberOfRows ) )
       u = 0.0d0

       RowsInUnityBlock = NumberOfRows - NbrRigids * 3 * (DOFs-1)

       C % Rows(1:RowsInUnityBlock+1) = (/ (i, i=1,RowsInUnityBlock+1) /)

       k = RowsInUnityBlock + 1
       DO j = 1, NbrRigids

          DO i = 1, DOFs
             C % Rows(k + i) = C % Rows(k + i - 1) + &
                  COUNT( RealNodeTypes == j )
          END DO

          DO i = 1, 2*DOFs-3
             C % Rows(k + DOFs + i) = C % Rows(k + DOFs + i - 1) + &
                  2 * COUNT( RealNodeTypes == j )
          END DO
          k = k + 3*(DOFs-1)

       END DO

       k = 0
       DO t = 1, Model % NumberOfNodes
          IF ( RealNodeTypes(t) /= 0 )  CYCLE
          DO i = 1, DOFs
             IF ( Permutation(DOFs * Perm(t) -DOFs + i ) == 0 )  CYCLE
             k = k + 1
             C % Cols(k) = Permutation( DOFs * Perm(t) - DOFs + i )
             C % Values(k) = 1.0d0
          END DO
       END DO

       DO j = 1, NbrRigids
          k = 0
          DO t = 1, Model % NumberOfNodes
             IF ( RealNodeTypes(t) /= j )  CYCLE
             DO i = 1, DOFs
                C % Cols( C % Rows( RowsInUnityBlock + (j-1)*3*(DOFs-1) + i ) + k ) &
                     = Permutation( DOFs * Perm(t) - DOFs + i )
                C % Values( C % Rows( RowsInUnityBlock + (j-1)*3*(DOFs-1) +i ) +k ) &
                     = 1.0d0
             END DO
             k = k + 1
          END DO
       END DO

    END IF
  end subroutine rbr_cstruct

  ! -- section 6: rigid-C values (original lines 632-712) ---------------
  subroutine rbr_cvals(h_model, h_solver, h_a, ndofs) bind(c, name="rbr_cvals")
    type(c_ptr), value :: h_model, h_solver, h_a
    integer(c_int), value :: ndofs
    type(model_t), pointer :: Model
    type(solver_t), pointer :: Solver
    type(matrix_t), pointer :: A
    integer :: j, k, t
    call c_f_pointer(h_model, Model)
    call c_f_pointer(h_solver, Solver)
    call c_f_pointer(h_a, A)
    DOFs = ndofs

    IF ( .NOT. AllocationsDone .OR. Solver % MeshChanged ) THEN

       DO j = 1, NbrRigids

          k = 0
          IF ( DOFs == 2 ) THEN
             DO t = 1, Model % NumberOfNodes
                IF ( RealNodeTypes(t) /= j )  CYCLE
                C % Cols( C % Rows( RowsInUnityBlock + (j-1)*3*(DOFs-1) &
                     + DOFs + 1 ) + k ) = &
                     Permutation( DOFs * Perm(t) - DOFs + 1 )
                C % Values( C % Rows( RowsInUnityBlock + (j-1)*3*(DOFs-1) &
                     + DOFs + 1 ) + k ) = &
                     CenterOfRigidBody(j,2) - Model % Mesh % Nodes % y(t)
                k = k + 1
                C % Cols( C % Rows( RowsInUnityBlock + (j-1)*3*(DOFs-1) &
                     + DOFs + 1 ) + k ) = &
                     Permutation( DOFs * Perm(t) - DOFs + 2 )
                C % Values( C % Rows( RowsInUnityBlock + (j-1)*3*(DOFs-1) &
                     + DOFs + 1 ) + k ) = &
                     Model % Mesh % Nodes % x(t) - CenterOfRigidBody(j,1)
                k = k + 1
             END DO
          ELSE
             DO t = 1, Model % NumberOfNodes
                IF ( RealNodeTypes(t) /= j ) CYCLE
                C % Cols( C % Rows( RowsInUnityBlock + (j-1)*3*(DOFs-1) &
                     + DOFs + 1 ) + k ) = &
                     Permutation( DOFs * Perm(t) - DOFs + 2 )
                C % Values( C % Rows( RowsInUnityBlock + (j-1)*3*(DOFs-1) &
                     + DOFs + 1 ) + k ) = &
                     Model % Mesh % Nodes % z(t) - CenterOfRigidBody(j,3)
                C % Cols( C % Rows( RowsInUnityBlock + (j-1)*3*(DOFs-1) &
                     + DOFs + 2 ) + k ) = &
                     Permutation( DOFs * Perm(t) - DOFs + 1 )
                C % Values( C % Rows( RowsInUnityBlock + (j-1)*3*(DOFs-1) &
                     + DOFs + 2 ) + k ) = &
                     CenterOfRigidBody(j,3) - Model % Mesh % Nodes % z(t)
                C % Cols( C % Rows( RowsInUnityBlock + (j-1)*3*(DOFs-1) &
                     + DOFs + 3 ) + k ) = &
                     Permutation( DOFs * Perm(t) - DOFs + 1 )
                C % Values( C % Rows( RowsInUnityBlock + (j-1)*3*(DOFs-1) &
                     + DOFs + 3 ) + k ) = &
                     Model % Mesh % Nodes % y(t) - CenterOfRigidBody(j,2)
                k = k + 1

                C % Cols( C % Rows( RowsInUnityBlock + (j-1)*3*(DOFs-1) &
                     + DOFs + 1 ) + k ) = &
                     Permutation( DOFs * Perm(t) - DOFs + 3 )
                C % Values( C % Rows( RowsInUnityBlock + (j-1)*3*(DOFs-1) &
                     + DOFs + 1 ) + k ) = &
                     CenterOfRigidBody(j,2) - Model % Mesh % Nodes % y(t)
                C % Cols( C % Rows( RowsInUnityBlock + (j-1)*3*(DOFs-1) &
                     + DOFs + 2 ) + k ) = &
                     Permutation( DOFs * Perm(t) - DOFs + 3 )
                C % Values( C % Rows( RowsInUnityBlock + (j-1)*3*(DOFs-1) &
                     + DOFs + 2 ) + k ) = &
                     Model % Mesh % Nodes % x(t) - CenterOfRigidBody(j,1)
                C % Cols( C % Rows( RowsInUnityBlock + (j-1)*3*(DOFs-1) &
                     + DOFs + 3 ) + k ) = &
                     Permutation( DOFs * Perm(t) - DOFs + 2 )
                C % Values( C % Rows( RowsInUnityBlock + (j-1)*3*(DOFs-1) &
                     + DOFs + 3 ) + k ) = &
                     CenterOfRigidBody(j,1) - Model % Mesh % Nodes % x(t)
                k = k + 1
             END DO
          END IF
       END DO


       C % Ordered = .FALSE.
       C % Diag = 0
       CALL CRS_SortMatrixValues( C )

       ALLOCATE( f ( C % NumberOfRows ) )

       NULLIFY( RedA )
       C_Trans => CRS_TransposeLocal( C )
       AllocationsDone = .TRUE.

    END IF
  end subroutine rbr_cvals

  ! -- section 7: eigen setup (original lines 718-733) ------------------
  subroutine rbr_eigsetup(h_solver, eig, neigen) bind(c, name="rbr_eigsetup")
    type(c_ptr), value :: h_solver
    integer(c_int) :: eig, neigen
    type(solver_t), pointer :: Solver
    logical :: EigAnal, GotIt
    call c_f_pointer(h_solver, Solver)

    EigAnal = ListGetLogical( Solver % Values, 'Eigen Analysis', GotIt )
    eig = 0
    neigen = 0
    IF ( EigAnal ) THEN
       NEigen = ListGetInteger( Solver % Values,  &
            'Eigen System Values', GotIt )

       IF ( GotIt .AND. NEigen > 0 ) THEN

          Solver % NOFEigenValues = NEigen

          Solver % Variable % EigenValues  = 0.0d0
          Solver % Variable % EigenVectors = 0.0d0

          ALLOCATE( LocalEigenVectors( NEigen, C % NumberOfRows ) )
          eig = 1

       END IF
    END IF
    neigen = NEigen
  end subroutine rbr_eigsetup

  ! -- section 8: reduced-matrix assembly (original lines 739-771) -----
  subroutine rbr_multiply(h_a, eig) bind(c, name="rbr_multiply")
    type(c_ptr), value :: h_a
    integer(c_int), value :: eig
    type(matrix_t), pointer :: A
    type(matrix_t), pointer :: A_Trans, BB, BB_Trans
    logical :: EigAnal
    call c_f_pointer(h_a, A)
    EigAnal = (eig /= 0)

    WRITE( Message, * ) 'Constructing the reduced matrix ...'
    CALL Info( 'RigidBody', Message, Level=5 )

    A_Trans => CRS_TransposeLocal( BMatrix, EigAnal, ASSOCIATED( BMatrix % DampValues ) )

    NULLIFY( BB )
    CALL CRS_MatrixMatrixMultiply( BB, A_Trans, C_Trans, EigAnal, &
         ASSOCIATED( A_Trans % DampValues ) )

    CALL LocalFreeMatrix( A_Trans )
    NULLIFY( A_Trans )

    BB_Trans => CRS_TransposeLocal( BB, EigAnal, ASSOCIATED( BMatrix % DampValues ) )

    CALL LocalFreeMatrix( BB )
    NULLIFY( BB )

    CALL Info( 'RigidBody', Message, Level=5 )
    CALL CRS_MatrixMatrixMultiply( RedA, BB_Trans, C_Trans, EigAnal, &
         ASSOCIATED( BB_Trans % DampValues ) )

    CALL LocalFreeMatrix( BB_Trans )

    FVector => BMatrix % RHS

    ALLOCATE( XVector( BMatrix % NumberOfRows ), &
         STAT=istat  )
    IF ( istat /= 0 )  CALL Fatal( 'RigidBody', 'Memory allocation error 3' )
    f = 0.0d0

    IF ( .NOT. EigAnal )  CALL CRS_MatrixVectorMultiply( C, FVector, f )
  end subroutine rbr_multiply

  ! -- section 9: bandwidth optimize (original lines 777-837) -----------
  subroutine rbr_optimize(h_solver, h_a) bind(c, name="rbr_optimize")
    type(c_ptr), value :: h_solver, h_a
    type(solver_t), pointer :: Solver
    type(matrix_t), pointer :: A
    integer :: i, j, HalfBandWidth
    logical :: GotIt
    call c_f_pointer(h_solver, Solver)
    call c_f_pointer(h_a, A)
    g_solver => Solver

    IF ( Verbose ) THEN

       HalfBandWidth = 0
       DO i = 1, A % NumberOfRows
          j = MAXVAL( ABS( A % Cols( A % Rows(i) : A % Rows(i+1)-1 ) - i ) )
          HalfBandWidth = MAX( HalfBandWidth, j )
       END DO
       WRITE ( Message, '( "Original matrix number of rows:      ", I6 )' ) &
            A % NumberOfRows
       CALL Info( 'RigidBody', Message, Level=24 )
       WRITE ( Message, '( "Half bandwidth of original matrix:   ", I6 )' ) &
            HalfBandWidth
       CALL Info( 'RigidBody', Message, Level=31 )

       HalfBandWidth = 0
       DO i = 1, BMatrix % NumberOfRows
          j = MAXVAL( ABS( BMatrix % Cols( BMatrix % Rows(i) : &
               BMatrix % Rows(i+1)-1 ) - i ) )
          HalfBandWidth = MAX( HalfBandWidth, j )
       END DO
       WRITE ( Message, '( "Eliminated matrix number of rows:    ", I6 )' ) &
            BMatrix % NumberOfRows
       CALL Info( 'RigidBody', Message, Level=24 )
       WRITE ( Message, '( "Half bandwidth of eliminated matrix: ", I6 )' ) &
            HalfBandWidth
       CALL Info( 'RigidBody', Message, Level=31 )

    END IF

    IF ( ListGetLogical( Solver % Values, 'Optimize Matrix Structure', GotIt ) ) THEN

       CALL Info( 'RigidBody', 'Optimizing matrix bandwidth', Level=5 )

       HalfBandWidth = 0
       DO i = 1, RedA % NumberOfRows
          j = MAXVAL( ABS( RedA % Cols( RedA % Rows(i) : RedA % Rows(i+1)-1 ) - i ) )
          HalfBandWidth = MAX( HalfBandWidth, j )
       END DO

       RedA % RHS => f
       CALL MatrixBandWidthOptimize( RedA, GotIt )
       f => RedA % RHS

       IF ( .NOT. GotIt ) THEN
         WRITE ( Message, '( "Reduced matrix number of rows:       ", I6 )' ) &
             RedA % NumberOfRows
         CALL Info( 'RigidBody', Message, Level=16 )
         WRITE ( Message, '( "Half bandwidth of reduced matrix:    ", I6 )' ) &
             HalfBandWidth
         CALL Info( 'RigidBody', Message, Level=16 )

         HalfBandWidth = 0
         DO i = 1, RedA % NumberOfRows
           j = MAXVAL( ABS( RedA % Cols( RedA % Rows(i) : RedA % Rows(i+1)-1 ) - i ) )
           HalfBandWidth = MAX( HalfBandWidth, j )
         END DO
         WRITE ( Message, '( "Half bandwidth of optimized matrix:  ", I6 )' ) &
             HalfBandWidth
         CALL Info( 'RigidBody', Message, Level=16 )
       END IF
     END IF
  end subroutine rbr_optimize

  ! -- section 10: reduced solve (original lines 843-848) ---------------
  subroutine rbr_solve(h_solver, h_a, Norm, ndofs) bind(c, name="rbr_solve")
    type(c_ptr), value :: h_solver, h_a
    integer(c_int), value :: ndofs
    real(kind=dp) :: Norm
    type(solver_t), pointer :: Solver
    type(matrix_t), pointer :: A
    real(kind=dp) :: at
    call c_f_pointer(h_solver, Solver)
    call c_f_pointer(h_a, A)
    DOFs = ndofs

    u = 0.0d0
    at = CPUtime()
    CALL SolveLinearSystem( RedA, f, u, Norm, DOFs, Solver )
    at = CPUtime() - at
    WRITE( Message, * ) 'Time used in SolveLinearSystem (CPU): ', at
    CALL Info( 'RigidBody', Message, Level=5 )
  end subroutine rbr_solve

  ! -- section 11: map back + teardown (original lines 853-1006) --------
  subroutine rbr_mapback(h_model, h_solver, h_a, b, x, n, ndofs, eig) bind(c, name="rbr_mapback")
    type(c_ptr), value :: h_model, h_solver, h_a
    integer(c_int), value :: n, ndofs, eig
    real(kind=dp) :: b(n), x(n)
    type(model_t), pointer :: Model
    type(solver_t), pointer :: Solver
    type(matrix_t), pointer :: A
    real(kind=dp), pointer :: ElimDiriEigenVectors(:,:)
    real(kind=dp), pointer :: u2(:)
    logical :: EigAnal, GotIt
    integer :: i, j, k
    call c_f_pointer(h_model, Model)
    call c_f_pointer(h_solver, Solver)
    call c_f_pointer(h_a, A)
    DOFs = ndofs
    EigAnal = (eig /= 0)

    IF ( EigAnal )  LocalEigenVectors(:, 1:RedA % NumberOfRows ) &
         = REAL( Solver % Variable % EigenVectors(:, 1:RedA % NumberOfRows) )

    IF ( ListGetLogical( Solver % Values, 'Optimize Matrix Structure', GotIt ) ) THEN
       IF ( EigAnal ) THEN
          LocalEigenVectors = 0.0d0
          DO i = 1, RedA % NumberOfRows
             LocalEigenVectors(:,i) = &
                  REAL( Solver % Variable % EigenVectors(:, RedA % Perm(i) ) )
          END DO

       ELSE
          ALLOCATE( u2( RedA % NumberOfRows ) )
          u2 = u
          DO i = 1, RedA % NumberOfRows
             u(i) = u2( RedA % Perm( i ) )
          END DO
          DEALLOCATE( u2 )
       END IF
    END IF

    IF ( EigAnal ) THEN
      Solver % Variable % EigenVectors = 0.0d0

      ALLOCATE( ElimDiriEigenVectors( NEigen, BMatrix % NumberOfRows ) )
      DO i = 1, NEigen
         CALL CRS_MatrixVectorMultiply( C_Trans, &
              LocalEigenVectors( i, 1: C % NumberOfRows ), &
              ElimDiriEigenVectors(i, 1: BMatrix % NumberOfRows) )
      END DO

      j = 0
      DO i = 1, A % NumberOFRows
         IF ( Permutation(i) == 0 ) THEN
            DO k = 1, NEigen
               Solver % Variable % EigenVectors(k,i) = b(i) / A % Values( A % Diag(i))
            END DO
         ELSE
            j = j + 1
            DO k = 1, NEigen
               Solver % Variable % EigenVectors(k,j) = ElimDiriEigenVectors(k,j)
            END DO
         END IF
      END DO

      IF ( .NOT. ASSOCIATED( BMatrix % DampValues ) ) THEN
         CALL Info( 'RigidBody', 'Eigenvalues  [Hz]:', Level=8 )
         CALL Info( 'RigidBody', ' ', Level=8 )

         DO i = 1, NEigen
            WRITE( Message, * ) '   ', i, &
                 SQRT( REAL(Solver % Variable % EigenValues(i) ) ) / (2.0d0 * PI )
            CALL Info( 'RigidBody', Message, Level=8 )
         END DO
      END IF

      DEALLOCATE( LocalEigenVectors )
      LocalEigenVectors => null()
      DEALLOCATE( ElimDiriEigenVectors )
    ELSE
       CALL Info( 'RigidBody', ' ', Level=5 )
       DO i = 1, NbrRigids
         IF ( DOFs == 2 ) THEN
            WRITE( Message, '("Rigid block ", I2, " (body ", I2, ") :" )' ) &
                 i, RigidIndex(i)
            CALL Info( 'RigidBody', Message )
            WRITE( Message, '("   Rigid body displacements: ", 2ES17.6E2 )' ) &
                 u( RedA % NumberOfRows - 3*NbrRigids + 3*i - 2 ), &
                 u( RedA % NumberOfRows - 3*NbrRigids + 3*i - 1 )
            CALL Info( 'RigidBody', Message )
            WRITE( Message, '("   Rigid body rotation:      ", ES17.6E2 )' ) &
                 u( RedA % NumberOfRows - 3*NbrRigids + 3*i )
            CALL Info( 'RigidBody', Message )

            WRITE( Message, '("res: Rigid block ", I2, " rotation 1" )' ) i
            CALL ListAddConstReal( Model % Simulation, &
                 Message, u( RedA % NumberOfRows - 3*NbrRigids + 3*i ) )
            WRITE( Message, '("res: Rigid block ", I2, " displacement 2" )' ) i
            CALL ListAddConstReal( Model % Simulation, &
                 Message, u( RedA % NumberOfRows - 3*NbrRigids + 3*i - 1) )
            WRITE( Message, '("res: Rigid block ", I2, " displacement 1" )' ) i
            CALL ListAddConstReal( Model % Simulation, &
                 Message, u( RedA % NumberOfRows - 3*NbrRigids + 3*i - 2) )
         ELSE
            WRITE( Message, '("Rigid block ", I2, " (body ", I2, ") :" )' ) &
                 i, RigidIndex(i)
            CALL Info( 'RigidBody', Message )
            WRITE( Message, '("Rigid body displacements: ", 3ES17.6E2 )' ) &
                 u( RedA % NumberOfRows - 6*NbrRigids + 6*i - 5), &
                 u( RedA % NumberOfRows - 6*NbrRigids + 6*i - 4), &
                 u( RedA % NumberOfRows - 6*NbrRigids + 6*i - 3)
            CALL Info( 'RigidBody', Message )
            WRITE( Message, '("Rigid body rotations:     ", 3ES17.6E2 )' ) &
                 u( RedA % NumberOfRows - 6*NbrRigids + 6*i - 2), &
                 u( RedA % NumberOfRows - 6*NbrRigids + 6*i - 1), &
                 u( RedA % NumberOfRows - 6*NbrRigids + 6*i)
            CALL Info( 'RigidBody', Message )

            WRITE( Message, '("res: Rigid block ", I2, " rotation 3" )' ) i
            CALL ListAddConstReal( Model % Simulation, &
                 Message, u( RedA % NumberOfRows - 6*NbrRigids + 6*i ) )
            WRITE( Message, '("res: Rigid block ", I2, " rotation 2" )' ) i
            CALL ListAddConstReal( Model % Simulation, &
                 Message, u( RedA % NumberOfRows - 6*NbrRigids + 6*i - 1) )
            WRITE( Message, '("res: Rigid block ", I2, " rotation 1" )' ) i
            CALL ListAddConstReal( Model % Simulation, &
                 Message, u( RedA % NumberOfRows - 6*NbrRigids + 6*i - 2) )

            WRITE( Message, '("res: Rigid block ", I2, " displacement 3" )' ) i
            CALL ListAddConstReal( Model % Simulation, &
                 Message, u( RedA % NumberOfRows - 6*NbrRigids + 6*i - 3) )
            WRITE( Message, '("res: Rigid block ", I2, " displacement 2" )' ) i
            CALL ListAddConstReal( Model % Simulation, &
                 Message, u( RedA % NumberOfRows - 6*NbrRigids + 6*i - 4) )
            WRITE( Message, '("res: Rigid block ", I2, " displacement 1" )' ) i
            CALL ListAddConstReal( Model % Simulation, &
                 Message, u( RedA % NumberOfRows - 6*NbrRigids + 6*i - 5) )
         END IF
      END DO
      CALL CRS_MatrixVectorMultiply( C_Trans, U, XVector )
    END IF

    IF ( .NOT. EigAnal ) THEN
      j = 0
      DO i=1, A % NumberOFRows
        IF ( Permutation(i) == 0 ) THEN
          x(i) = b(i) / A % Values( A % Diag(i) )
        ELSE
          j = j + 1
          X(i) = XVector(j)
        END IF
      END DO
    END IF

    DEALLOCATE( XVector )
    XVector => null()

    DEALLOCATE( RedA % Rows )
    DEALLOCATE( RedA % Cols )
    DEALLOCATE( RedA % Values )
    IF ( ASSOCIATED( RedA % Diag ) )  DEALLOCATE( RedA % Diag )
    IF ( ASSOCIATED( RedA % MassValues ) )  DEALLOCATE( RedA % MassValues )
    IF ( ASSOCIATED( RedA % DampValues ) )  DEALLOCATE( RedA % DampValues )

    TotTime = CPUtime() - TotTime

    WRITE( Message, * ) 'Total time spent in RigidBodyReduction (CPU): ', TotTime
    CALL Info( 'RigidBody', ' ', Level=5 )
    CALL Info( 'RigidBody', Message, Level=5 )
    CALL Info( 'RigidBody', ' ', Level=5 )
  end subroutine rbr_mapback

  ! -- CONTAINS internals copied verbatim as module procedures -----------
  ! -- (original lines 1010-1833; see header for the two adaptations) -----
  SUBROUTINE ComputeMassCenter( Centers, Model, Solver, Indx, n )
    USE Lists
    USE Integration
    IMPLICIT NONE

    TYPE(Model_t) :: Model
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: Centers(:,:)
    INTEGER :: n, Indx( Model % NumberOfBodies )

    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    TYPE(Element_t), POINTER :: CurrentElement
    TYPE(Nodes_t) :: ElementNodes
    REAL(KIND=dp), POINTER :: Basis(:), dBasisdx(:,:)
    REAL(KIND=dp), POINTER :: Density(:)
    REAL(KIND=dp) :: detJ, xpos, ypos, zpos, uu, vv, ww, s
    REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)
    REAL(KIND=dp) :: ElementMassCenter(3), ElementMass, Dens, TotalMass(n)
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER :: i, j, m, t, k, p, istat
    LOGICAL :: stat
    LOGICAL :: Rigid
    LOGICAL :: GotIt

    Centers = 0.0d0

    m = Model % MaxElementNodes

    ALLOCATE( ElementNodes % x( m ), &
         ElementNodes % y( m ), &
         ElementNodes % z( m ), &
         Density( m ),          &
         Basis( m ),            &
         dBasisdx( m,3 ),       &
         STAT=istat )
    IF ( istat /= 0 )  CALL Fatal( 'ComputeMassCenter', 'Memory allocation error' )

    TotalMass = 0.0d0

    DO t = 1, Solver % NumberOfActiveElements

       CurrentElement => Solver % Mesh % Elements( Solver % ActiveElements(t) )
       NodeIndexes => CurrentElement % NodeIndexes
       m = CurrentElement % TYPE % NumberOfNodes

       Rigid = .FALSE.
       DO i = 1, Model % NumberOfBodies
          IF ( CurrentElement % BodyId == Indx(i) ) THEN
             Rigid = .TRUE.
             EXIT
          END IF
       END DO
       IF ( .NOT. Rigid )  CYCLE

       k = ListGetInteger( Model % Bodies( CurrentElement % BodyId ) &
            % Values, 'Material' )

       Density(1:m) = ListGetReal( Model % Materials(k) % Values, &
            'Density', m, NodeIndexes, gotIt )
       IF ( .NOT. GotIt ) Density(1:m) = 1.0d0

       ElementNodes % x(1:m) = Solver % Mesh % Nodes % x(NodeIndexes)
       ElementNodes % y(1:m) = Solver % Mesh % Nodes % y(NodeIndexes)
       ElementNodes % z(1:m) = Solver % Mesh % Nodes % z(NodeIndexes)

       ElementMassCenter = 0.0d0
       ElementMass = 0.0d0

       IntegStuff = GaussPoints( CurrentElement )

       DO p = 1, IntegStuff % n
          uu = IntegStuff % u(p)
          vv = IntegStuff % v(p)
          ww = IntegStuff % w(p)
          s = IntegStuff % s(p)
          stat = ElementInfo( CurrentElement, ElementNodes, uu, vv, ww, detJ, &
               Basis, dBasisdx )

          xpos = SUM( ElementNodes % x(1:m)*Basis(1:m) )
          ypos = SUM( ElementNodes % y(1:m)*Basis(1:m) )
          zpos = SUM( ElementNodes % z(1:m)*Basis(1:m) )
          IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
             s = s * 2 * pi
          END IF

          CALL CoordinateSystemInfo( Metric, SqrtMetric, Symb, dSymb, &
               xpos, ypos, zpos )

          S = S * detJ * SqrtMetric

          Dens = SUM( Density(1:m) * Basis(1:m) )

          ElementMassCenter(1) = ElementMassCenter(1) + s * xpos * Dens
          ElementMassCenter(2) = ElementMassCenter(2) + s * ypos * Dens
          IF ( DOFs == 3 )  &
               ElementMassCenter(3) = ElementMassCenter(3) + s * zpos * Dens

          ElementMass = ElementMass + s * Dens

       END DO

       TotalMass(i) = TotalMass(i) + ElementMass
       DO j = 1, 3
          CenterOfRigidBody(i,j) = CenterOfRigidBody(i,j) + ElementMassCenter(j)
       END DO

    END DO

    CALL Info( 'RigidBody', ' ', Level=12 )
    DO j = 1, n
       CenterOfRigidBody(j,:) = CenterOfRigidBody(j,:) / TotalMass(j)

       WRITE( Message, '(" Center of Rigid body:", 3ES17.6E2)') &
            CenterOfRigidBody(j,:)
       CALL Info( 'RigidBody', Message, Level=12 )
       WRITE( Message, '(" Mass of rigid body:  ", ES17.6E2)') TotalMass(j)
       CALL Info( 'RigidBody', Message, Level=24 )
    END DO
    CALL Info( 'RigidBody', ' ', Level=12 )

    DEALLOCATE( ElementNodes % x, ElementNodes % y, ElementNodes % z, &
         Density, Basis, dBasisdx )
  END SUBROUTINE ComputeMassCenter
  SUBROUTINE CRS_SortMatrixValues( A )
    USE CRSMatrix
    USE GeneralUtils
    IMPLICIT NONE

    TYPE(Matrix_t), POINTER :: A

    LOGICAL :: Mass, Damp

    INTEGER :: i, j, n, k

    INTEGER, ALLOCATABLE :: RowPerm(:)
    INTEGER, POINTER :: Cols(:), Rows(:), Diag(:)
    REAL(KIND=dp), POINTER CONTIG :: Values(:), MassValues(:), DampValues(:)


    Diag   => A % Diag
    Rows   => A % Rows
    Cols   => A % Cols

    n = A % NumberOfRows

    Mass = ASSOCIATED( A % MassValues )
    Damp = ASSOCIATED( A % DampValues )

    IF ( .NOT. A % Ordered ) THEN
       IF ( Mass .OR. Damp ) THEN

          ALLOCATE( Values( A % Rows( A % NumberOfRows + 1 ) - 1 ) )
          IF ( Mass )  ALLOCATE( MassValues( A % Rows( A % NumberOfRows + 1 ) - 1 ) )
          IF ( Damp )  ALLOCATE( DampValues( A % Rows( A % NumberOfRows + 1 ) - 1 ) )

          j = 0
          DO i = 1, n
             j = MAX( j, Rows(i+1) - Rows(i) )
          END DO
          ALLOCATE( RowPerm( j ) )

          DO i = 1, n
             RowPerm = 0
             j = Rows(i+1) - Rows(i)
             RowPerm = (/ (k, k=1,j) /)
             CALL SortI( Rows(i+1)-Rows(i), Cols(Rows(i):Rows(i+1)-1), &
                  RowPerm(1:j) )
             DO k = 1, j
                Values(Rows(i) + k - 1) = A % Values(Rows(i) + RowPerm(k) - 1)
                IF ( Mass )  MassValues(Rows(i) + k - 1) = &
                     A % MassValues(Rows(i) + RowPerm(k) - 1)
                IF ( Damp )  DampValues(Rows(i) + k - 1) = &
                     A % DampValues(Rows(i) + RowPerm(k) - 1)
             END DO
          END DO

          DEALLOCATE( RowPerm )
          DEALLOCATE( A % Values )
          A % Values => Values

          IF ( Mass )  THEN
             DEALLOCATE( A % MassValues )
             A % MassValues => MassValues
          END IF

          IF ( Damp ) THEN
             DEALLOCATE( A % DampValues )
             A % DampValues => DampValues
          END IF

       ELSE

          Values => A % Values
          DO i = 1, n
             CALL SortF( Rows(i+1)-Rows(i), Cols(Rows(i):Rows(i+1)-1), &
                  Values(Rows(i):Rows(i+1)-1) )
          END DO

       END IF

       DO i = 1, n
          DO j = Rows(i), Rows(i+1)-1
             IF ( Cols(j) == i ) THEN
                Diag(i) = j
                EXIT
             END IF
          END DO
       END DO

       A % Ordered = .TRUE.
    END IF

  END SUBROUTINE CRS_SortMatrixValues

  FUNCTION CRS_TransposeLocal( A, MVal, DVal ) RESULT(B)
    USE CRSMatrix
    IMPLICIT NONE

    TYPE(Matrix_t), POINTER :: A, B
    LOGICAL, OPTIONAL :: MVal, DVal

    INTEGER, ALLOCATABLE :: Row(:)
    INTEGER :: NVals
    INTEGER :: i,j,k
    LOGICAL :: Mass, Damp

    Mass = .FALSE.
    Damp = .FALSE.
    IF ( PRESENT( MVal ) )  Mass = MVal
    IF ( PRESENT( DVal ) )  Damp = DVal

    B => AllocateMatrix()

    NVals = SIZE( A % Values )
    B % NumberOfRows = MAXVAL( A % Cols )

    ALLOCATE( B % Rows( B % NumberOfRows +1 ), B % Cols( NVals ), &
        B % Values( Nvals ), B % Diag( B % NumberOfRows ), STAT=istat )
    IF ( istat /= 0 )  CALL Fatal( 'CRS_TransposeLocal', &
        'Memory allocation error.' )
    IF ( Mass )  ALLOCATE( B % MassValues( NVals ) )
    IF ( Damp )  ALLOCATE( B % DampValues( NVals ) )

    B % Diag = 0

    ALLOCATE( Row( B % NumberOfRows ) )
    Row = 0

    DO i = 1, NVals
      Row( A % Cols(i) ) = Row( A % Cols(i) ) + 1
    END DO

    B % Rows(1) = 1
    DO i = 1, B % NumberOfRows
      B % Rows(i+1) = B % Rows(i) + Row(i)
    END DO
    B % Cols = 0

    DO i = 1, B % NumberOfRows
       Row(i) = B % Rows(i)
    END DO
    IF ( Mass ) THEN
       DO i = 1, A % NumberOfRows
          DO j = A % Rows(i), A % Rows(i+1) - 1
             k = A % Cols(j)
             IF ( Row(k) < B % Rows(k+1) ) THEN
                B % Cols( Row(k) ) = i
                B % Values( Row(k) ) = A % Values(j)
                B % MassValues( Row(k) ) = A % MassValues(j)
                IF ( Damp )  B % DampValues( Row(k) ) = A % DampValues(j)
                Row(k) = Row(k) + 1
             ELSE
                WRITE( Message, * ) 'Trying to access non-existent column', i,k
                CALL Error( 'CRS_TransposeLocal', Message )
                RETURN
             END IF
          END DO
       END DO
    ELSE
       DO i = 1, A % NumberOfRows
          DO j = A % Rows(i), A % Rows(i+1) - 1
             k = A % Cols(j)
             IF ( Row(k) < B % Rows(k+1) ) THEN
                B % Cols( Row(k) ) = i
                B % Values( Row(k) ) = A % Values(j)
                Row(k) = Row(k) + 1
             ELSE
                WRITE( Message, * ) 'Trying to access non-existent column', i,k
                CALL Error( 'CRS_TransposeLocal', Message )
                RETURN
             END IF
          END DO
       END DO
    END IF

    DEALLOCATE( Row )

  END FUNCTION CRS_TransposeLocal

  SUBROUTINE CRS_MatrixMatrixMultiply( C, A, B, MVal, DVal )
    USE CRSMatrix
    IMPLICIT NONE

    TYPE(Matrix_t), POINTER :: A, B, C
    LOGICAL :: Mval
    LOGICAL, OPTIONAL :: DVal

    REAL(KIND=dp) , POINTER :: AMassVals(:), BMassVals(:)
    REAL(KIND=dp) , POINTER :: ADampVals(:), BDampVals(:)
    LOGICAL, ALLOCATABLE :: Row(:)
    LOGICAL :: Damp, Allocated
    INTEGER, ALLOCATABLE :: Ind(:)
    INTEGER :: ki, NCols
    INTEGER :: i,j,k,l,m,n,Rownonzeros( A % NumberOfRows ), TotalNonzeros

    Damp = .FALSE.
    IF ( PRESENT( DVal ) )  Damp = DVal

    IF ( .NOT. ASSOCIATED( C ) ) THEN
       Allocated = .FALSE.
    ELSE IF ( .NOT. ASSOCIATED( C % Values ) ) THEN
       Allocated = .FALSE.
    ELSE
       Allocated = .TRUE.
    END IF


    IF ( .NOT. Allocated ) THEN

       Rownonzeros   = 0
       TotalNonzeros = 0

       NCols = MAXVAL( B % Cols )
       ALLOCATE( Row( NCols ), Ind( NCols ) )

       Row = .FALSE.
       DO i=1,A % NumberOfRows
          ki = 0
          DO l=A % Rows(i),A % Rows(i+1)-1
             k = A % Cols(l)
             DO m=B % Rows(k), B % Rows(k+1)-1
                j = B % Cols(m)

                IF ( .NOT. Row(j) ) THEN
                   ki = ki + 1
                   Ind(ki) = j
                   Row(j) = .TRUE.
                END IF
             END DO
          END DO

          DO j=1,ki
             Row(Ind(j)) = .FALSE.
          END DO
          RowNonzeros(i) = ki
          TotalNonzeros  = TotalNonzeros + ki
       END DO

       IF ( .NOT. ASSOCIATED( C ) )  C => AllocateMatrix()

       C % NumberOfRows = A % NumberOFRows
       ALLOCATE( C % Cols( TotalNonzeros ), C % Values( TotalNonzeros ) )
       ALLOCATE( C % Rows( C % NumberOfRows + 1 ), &
            C % Diag( C % NumberOfRows ) )
       IF ( MVal )  ALLOCATE( C % MassValues( TotalNonzeros ) )
       IF ( Damp )  ALLOCATE( C % DampValues( TotalNonzeros ) )

       C % Diag = 0
       C % Rows(1) = 1
       DO i=1, A % NumberOfRows
          C % Rows(i+1) = C % Rows(i) + Rownonzeros(i)
       END DO

       C % Cols = 0
       Row = .FALSE.
       DO i=1,A % NumberOfRows
          ki = 0
          DO l=A % Rows(i),A % Rows(i+1)-1
             k = A % Cols(l)
             DO m=B % Rows(k), B % Rows(k+1)-1
                j = B % Cols(m)

                IF ( .NOT. Row(j) ) THEN
                   ki = ki + 1
                   Ind(ki) = j
                   Row(j) = .TRUE.
                   CALL CRS_MakeMatrixIndex( C, i, j )
                END IF
             END DO
          END DO

          DO j=1,ki
             Row(Ind(j)) = .FALSE.
          END DO
       END DO

       CALL CRS_SortMatrix( C )

       DEALLOCATE( Row, Ind )

    END IF

    IF ( MVal ) THEN
       IF ( ASSOCIATED( A % MassValues ) ) THEN
          AMassVals => A % MassValues
       ELSE
          AMassVals => A % Values
       END IF

       IF ( ASSOCIATED( B % MassValues ) ) THEN
          BMassVals => B % MassValues
       ELSE
          BMassVals => B % Values
       END IF

       IF ( Damp ) THEN
          IF ( ASSOCIATED( A % DampValues ) ) THEN
             ADampVals => A % DampValues
          ELSE
             ADampVals => A % Values
          END IF

          IF ( ASSOCIATED( B % DampValues ) ) THEN
             BDampVals => B % DampValues
          ELSE
             BDampVals => B % Values
          END IF
          C % DampValues = 0.0d0
       END IF

       C % Values = 0.0d0
       C % MassValues = 0.0d0

       DO i=1,A % NumberOfRows
          DO l=A % Rows(i),A % Rows(i+1)-1
             k = A % Cols(l)
             DO m=B % Rows(k), B % Rows(k+1)-1
                j = B % Cols(m)
                DO n=C % Rows(i), C % Rows(i+1)-1
                   IF ( C % Cols(n) == j ) THEN
                      C % Values(n) = C % Values(n) + &
                           A % Values(l) * B % Values(m)
                      C % MassValues(n) = C % MassValues(n) + &
                           AMassVals(l) * BMassVals(m)
                      IF ( Damp )  C % DampValues(n) = C % DampValues(n) + &
                           ADampVals(l) * BDampVals(m)
                      EXIT
                   END IF
                END DO
             END DO
          END DO
       END DO
    ELSE
       C % Values = 0.0d0
       DO i=1,A % NumberOfRows
          DO l=A % Rows(i),A % Rows(i+1)-1
             k = A % Cols(l)
             DO m=B % Rows(k), B % Rows(k+1)-1
                j = B % Cols(m)
                DO n=C % Rows(i), C % Rows(i+1)-1
                   IF ( C % Cols(n) == j ) THEN
                      C % Values(n) = C % Values(n) + &
                           A % Values(l) * B % Values(m)
                      EXIT
                   END IF
                END DO
             END DO
          END DO
       END DO
    END IF

  END SUBROUTINE CRS_MatrixMatrixMultiply

  SUBROUTINE LocalFreeMatrix( A )
    USE GeneralUtils
    IMPLICIT NONE

    TYPE(matrix_t), POINTER :: A

    DEALLOCATE( A % Cols )
    DEALLOCATE( A % Rows )
    DEALLOCATE( A % Diag )
    DEALLOCATE( A % Values )
    IF ( ASSOCIATED( A % MassValues ) )  DEALLOCATE( A % MassValues )
    IF ( ASSOCIATED( A % DampValues ) )  DEALLOCATE( A % DampValues )
    IF ( ASSOCIATED( A % RHS ) )  DEALLOCATE( A % RHS )
    IF ( ASSOCIATED( A % Perm ) )  DEALLOCATE( A % Perm )
    DEALLOCATE( A )
  END SUBROUTINE LocalFreeMatrix
  SUBROUTINE MatrixBandWidthOptimize( A, DoneAlready )
    USE CRSMatrix
    IMPLICIT NONE

    TYPE(Matrix_t), POINTER :: A
    LOGICAL, OPTIONAL :: DoneAlready

    REAL(KIND=dp), POINTER CONTIG :: Vals(:)
    REAL(KIND=dp), POINTER CONTIG :: MassVals(:), DampVals(:), b(:)
    INTEGER, POINTER CONTIG :: Rows(:), Cols(:), Diag(:), RowEntrys(:)
    INTEGER, POINTER :: NewPerm(:), NbrNeighbors(:,:), TempCol(:)
    INTEGER :: RowCompleted, RowAdded, Lag, TempMin
    INTEGER :: i, j, k, n, istat
    LOGICAL, POINTER :: Visited(:)
    LOGICAL :: Mass, Damp, GotIt

    IF ( PRESENT( DoneAlready) )  DoneAlready = .TRUE.

    IF ( .NOT. ASSOCIATED( A % Perm ) ) THEN

       IF ( PRESENT( DoneAlready ) )  DoneAlready = .FALSE.

       ALLOCATE( NewPerm( A % NumberOfRows ), &
            NbrNeighbors( A % NumberOfRows, 2 ), &
            Cols( A % Rows( A % NumberOfRows + 1 ) - 1 ), &
            STAT = istat )
       IF ( istat /= 0 )  CALL Fatal( 'MatrixBandWidthOptimize', &
            'Memory allocation error 1' )
       NewPerm = 0
       Cols = A % Cols

       TempMin = 100
       DO i = 1, A % NumberOfRows
          NbrNeighbors(i,1) = i
          NbrNeighbors(i,2) = A % Rows(i+1) - A % Rows(i)
          IF ( NbrNeighbors(i,2) < TempMin ) THEN
             TempMin = NbrNeighbors(i,2)
             j = i
          END IF
       END DO

       ALLOCATE( TempCol( MAXVAL( NbrNeighbors(:,2) ) ) )

       DO i = 1, A % NumberOfRows
          n = NbrNeighbors(i,2)
          TempCol(1:n) = NbrNeighbors( Cols( A % Rows(i) : A % Rows(i+1) - 1 ), 2 )

          CALL SortI( n, TempCol(1:n), &
               Cols(A % Rows(i) : A % Rows(i+1) - 1 ) )
       END DO

       ALLOCATE( Visited( A % NumberOfRows ) )
       Visited = .FALSE.

       NewPerm(1) = NbrNeighbors(j,1)
       Visited(NbrNeighbors(j,1)) = .TRUE.

       RowCompleted = 0
       RowAdded = 1

       DO WHILE ( RowAdded < A % NumberOfRows )

          Lag = RowAdded - RowCompleted
          IF ( Lag == 0 ) THEN
             PRINT*, 'No success'
             DO k = 1, A % NumberOfRows
                IF ( .NOT. Visited(k) )  EXIT
             END DO
             NewPerm( RowAdded + 1 ) = k
             RowAdded = RowAdded + 1
             Lag = RowAdded - RowCompleted
          END IF

          DO k = 1, Lag

             i = NewPerm( RowCompleted + k )

             n = 0
             DO j = 1, NbrNeighbors(i,2)
                IF ( Visited( Cols( A % Rows(i) + j - 1 ) ) ) CYCLE
                n = n + 1
                NewPerm( RowAdded + n ) = Cols( A % Rows(i) + j - 1 )
                Visited( Cols( A % Rows(i) + j - 1 ) ) = .TRUE.
             END DO

             RowAdded = RowAdded + n

             IF ( RowAdded == A % NumberOfRows )  EXIT
          END DO
          RowCompleted = RowCompleted + Lag

       END DO

       DEALLOCATE( Cols )
       DEALLOCATE( Visited )
       DEALLOCATE( NbrNeighbors )
       DEALLOCATE( TempCol )

       ALLOCATE( TempCol( A % NumberOfRows ), STAT=istat )
       IF ( istat /= 0 )  CALL Fatal( 'MatrixBandWidthOptimize', &
            'Memory allocation error 2' )

       IF ( ListGetLogical( g_solver % Values, 'Reverse Ordering', GotIt ) ) THEN

          TempCol = NewPerm
          DO i = 1, A % NumberOfRows
             NewPerm(i) = TempCol( A % NumberOfRows + 1 - i )
          END DO

          IF ( Verbose )  &
               CALL Info( 'RigidBody', 'Matrix ordering reversed', Level=16 )

       END IF

       TempCol = NewPerm

       DO i = 1, A % NumberOfRows
          NewPerm( TempCol(i) ) = i
       END DO

       DEALLOCATE( TempCol )

       A % Perm => NewPerm

    END IF

    NewPerm => A % Perm

    ALLOCATE( Rows( A % NumberOfRows + 1 ), &
              Diag( A % NumberOfRows ), &
              Cols( A % Rows( A % NumberOfRows + 1 ) - 1 ), &
              Vals( A % Rows( A % NumberOfRows + 1 ) - 1 ), &
              STAT=istat )
    IF ( istat /= 0 )  CALL Fatal( 'MatrixBandWidthOptimize', &
         'Memory allocation error 3' )

    Mass = .FALSE.
    IF ( ASSOCIATED( A % MassValues ) ) THEN
       Mass = .TRUE.
       ALLOCATE( MassVals( A % Rows( A % NumberOfRows + 1 ) - 1 ), &
            STAT=istat )
       IF ( istat /= 0 )  CALL Fatal( 'MatrixBandWidthOptimize', &
            'Rare memory allocation error' )
       MassVals = 0.0d0
    END IF
    Damp = .FALSE.
    IF ( ASSOCIATED( A % DampValues ) ) THEN
       Damp = .TRUE.
       ALLOCATE( DampVals( A % Rows( A % NumberOfRows + 1 ) - 1 ), &
            STAT=istat )
       IF ( istat /= 0 )  CALL Fatal( 'MatrixBandWidthOptimize', &
            'Unusual memory allocation error' )
       DampVals = 0.0d0
    END IF
    Vals = 0.0d0

    ALLOCATE( RowEntrys( A % NumberOfRows ), STAT=istat )
    IF ( istat /= 0 )  CALL Fatal( 'MatrixBandWidthOptimize', &
         'Memory allocation error 2.1' )

    RowEntrys = 0
    DO i = 1, A % NumberOfRows
       RowEntrys( NewPerm(i) ) = A % Rows(i+1) - A % Rows(i)
    END DO

    Rows(1) = 1
    DO i = 1, A % NumberOfRows
       Rows(i+1) = Rows(i) + RowEntrys(i)
    END DO

    DEALLOCATE( RowEntrys )

    DO i = 1, A % NumberOfRows
       DO j = 1, Rows(NewPerm(i)+1) - Rows(NewPerm(i))
          Cols(Rows(NewPerm(i)) + j - 1) = &
               NewPerm( A % Cols( A % Rows(i) + j - 1 ) )
          Vals(Rows(NewPerm(i)) + j - 1) = A % Values( A % Rows(i) + j - 1 )
          IF ( Mass )  MassVals(Rows(NewPerm(i)) + j - 1) = &
               A % MassValues( A % Rows(i) + j - 1 )
          IF ( Damp )  DampVals(Rows(NewPerm(i)) + j - 1) = &
               A % DampValues( A % Rows(i) + j - 1 )
          IF ( Cols(Rows(NewPerm(i)) + j - 1) == NewPerm(i) ) &
               Diag(NewPerm(i)) = Rows( NewPerm(i)) + j - 1
       END DO
    END DO

    IF ( ASSOCIATED( A % Rows ) )  DEALLOCATE( A % Rows )
    IF ( ASSOCIATED( A % Cols ) )  DEALLOCATE( A % Cols )
    IF ( ASSOCIATED( A % Values ) )  DEALLOCATE( A % Values )
    IF ( ASSOCIATED( A % Diag ) )  DEALLOCATE( A % Diag )

    A % Rows => Rows
    A % Cols => Cols
    A % Diag => Diag
    A % Values => Vals
    A % Ordered = .FALSE.

    IF ( Mass ) THEN
       DEALLOCATE( A % MassValues )
       A % MassValues => MassVals
    END IF

    IF ( Damp ) THEN
       DEALLOCATE( A % DampValues )
       A % DampValues => DampVals
    END IF

    A % Ordered = .FALSE.

    CALL CRS_SortMatrixValues( A )

    IF ( ASSOCIATED( A % RHS ) ) THEN
       ALLOCATE( b( A % NumberOfRows ), STAT=istat )
       IF ( istat /= 0 )  CALL Fatal( 'MatrixBandWidthOptimize', &
            'Memory allocation error 4' )

       DO i = 1, A % NumberOfRows
          b(NewPerm(i)) = A % RHS(i)
       END DO

       DEALLOCATE( A % RHS )
       A % RHS => b
    END IF

  END SUBROUTINE MatrixBandWidthOptimize

end module rbr_acc
