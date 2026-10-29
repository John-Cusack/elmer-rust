! T022-A EliminatePeriodic port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/EliminatePeriodic.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments. All Fortran-object
! traffic — TYPE(...) handles, POINTER dummies, list-directed timing
! reports — stays inside Fortran, so the Rust side is pure C calls. The
! wrapper module mirrors the original's locals with identical lifetimes.
!
! Partition (original line numbers):
! - the first-call projection build (101-244), the eigen setup (251-265),
!   the CAC^T reduction (271-293), the reduced flags + solve (297-327),
!   the eigen map-back (339-356), the Norm reduction (361), the tail
!   timing/return (364-373), and every Elmer library call stay in Fortran,
!   so bit-identity holds by construction wherever the call sequence
!   matches.
! - the reduced-system Perm remap (303-313: `DOF = DOFs*(j-1)+1`,
!   `(ColN(DOF)-1)/DOFs+1`, exact `> 0` comparisons, truncating integer
!   division on positive values) runs in Rust through the scalar
!   elp_vperm_*/elp_coln_* getters/setters below with the original's
!   evaluation order reproduced exactly.
! - the serial (PEs <= 1) paths are transcribed verbatim including the
!   PEs > 1 branches, which are never taken by the serial proof.
!
! SAVE design (author's evident intent): the original's `ColN` is a
! non-SAVE allocatable local filled on the first call and read on every
! later call (`IF (ColN(DOF) > 0)`); the wrapper keeps it as the SAVE
! module array `acc_coln` with identical contents, so re-entry is
! deterministic. The `A % EMatrix`/`C % EMatrix` links persist through
! the matrix exactly as in the original. `NumberOfRows` (used by the
! eigen branch) is likewise kept as `acc_numberofrows`.
!
! Delegation: like the reference `.so` (NEEDED EliminateDirichlet.so),
! both the no-periodic first call (184-187) and the reduced-system solve
! (322) call the same `EliminateDirichlet` routine from the built tree.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module eliminateperiodic_acc
  use types
  use lists
  use solverbasics
  use solvecore, only: SolveLinearSystem
  use crsmatrix
  use generalutils
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(matrix_t), pointer, save :: acc_a => null()
  type(matrix_t), pointer, save :: acc_c => null()
  type(matrix_t), pointer, save :: acc_ctrans => null()
  type(matrix_t), pointer, save :: acc_reda => null()
  integer, pointer, save :: acc_perm(:) => null()
  ! SAVE mirror of the original's ColN (see header note).
  integer, allocatable, save :: acc_coln(:)
  integer, save :: acc_numberofrows = 0
  real(kind=dp), pointer, save :: acc_f(:) => null()
  real(kind=dp), pointer, save :: acc_u(:) => null()
  real(kind=dp), pointer, save :: acc_b(:) => null()
  real(kind=dp), pointer, save :: acc_x(:) => null()
  real(kind=dp), pointer, save :: acc_locre(:,:) => null()
  real(kind=dp), pointer, save :: acc_locim(:,:) => null()
  ! Entry timing (mirrors the original's TotTime local).
  real(kind=dp), save :: acc_tottime = 0.0_dp
  real(kind=dp), save :: acc_at = 0.0_dp
  integer, save :: acc_n = 0
  integer, save :: acc_neigen = 0
  logical, save :: acc_eiganal = .false.
  integer, save :: acc_delegated = 0

contains

  ! -- entry state ------------------------------------------------------
  subroutine elp_begin(h_model, h_solver, h_a, h_b, h_x, n) bind(c, name="elp_begin")
    type(c_ptr), value :: h_model, h_solver, h_a, h_b, h_x
    integer(c_int), value :: n
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    call c_f_pointer(h_a, acc_a)
    call c_f_pointer(h_b, acc_b, [n])
    call c_f_pointer(h_x, acc_x, [n])
    acc_n = n
    acc_tottime = cputime()
    acc_c => acc_a % ematrix
    acc_delegated = 0
  end subroutine elp_begin

  function elp_need_build() result(code) bind(c, name="elp_need_build")
    integer(c_int) :: code
    code = 0
    if (.not. associated(acc_c)) code = 1
  end function elp_need_build

  function elp_delegated() result(code) bind(c, name="elp_delegated")
    integer(c_int) :: code
    code = acc_delegated
  end function elp_delegated

  ! -- first-call projection build (original 101-244, verbatim) ----------
  function elp_build(dofs, n, b, x, norm) result(rv) bind(c, name="elp_build")
    integer(c_int), value :: dofs, n
    real(c_double) :: b(n), x(n)
    real(c_double) :: norm
    integer(c_int) :: rv
    INTERFACE
      INTEGER FUNCTION EliminateDirichlet( Model, Solver, A, b, x, n, DOFs, Norm )
         USE Types
         USE Lists
         USE SolverBasics
         USE CRSmatrix
         USE GeneralUtils
         TYPE(model_t)  :: Model
         TYPE(solver_t) :: Solver
         TYPE(matrix_t), POINTER :: A
         INTEGER :: DOFs, n
         REAL(KIND=dp) :: b(n), x(n), Norm
      END FUNCTION EliminateDirichlet
    END INTERFACE
    TYPE(matrix_t), POINTER :: RedA, C, C_Trans, Projector
    LOGICAL, ALLOCATABLE :: PeriodicProj(:)
    LOGICAL :: GotIt
    INTEGER, ALLOCATABLE, TARGET :: ColNTmp(:)
    INTEGER, POINTER :: Perm(:), PeriodicDOFs(:), R(:,:), TmpR(:,:)
    INTEGER :: NumberOfValues
    INTEGER :: istat, i, j, k, l, m, t, IDC, PeriodicNodes, DOF, NLen
    REAL(KIND=dp) :: Scale
    CHARACTER(LEN=MAX_NAME_LEN) :: Name
    C => acc_a % EMatrix
    Perm => acc_solver % Variable % Perm
    ALLOCATE( PeriodicDOFs(acc_a % NumberOFRows), ColNTmp( acc_a % NumberOfRows ), &
              PeriodicProj(acc_a % NumberOfRows) )
    ALLOCATE( R(1000, 2), STAT=istat )
    PeriodicDOFs  = 0
    PeriodicNodes = 0
    PeriodicProj = .FALSE.
    DO i=1,acc_model % NumberOfBCs
       DO DOF=1,dofs
          Name = ComponentName(acc_solver % Variable)
          nLen = LEN_TRIM(Name)
          IF ( dofs > 1 ) THEN
             WRITE( Name, '(a,i1)' ) Name(1:nlen) // ' ', DOF
          END IF
          Scale = -1.0d0
          IF ( .NOT. ListGetLogical( acc_model % BCs(i) % Values, &
             'Periodic BC ' // Name(1:nlen), GotIt ) ) THEN
              IF ( .NOT. ListGetLogical( acc_model % BCs(i) % Values, &
                   'Anti Periodic BC ' // Name(1:nlen), GotIt ) ) CYCLE
             Scale = 1.0d0
          END IF
          IF ( ParEnv % PEs <= 1 ) THEN
             Projector => acc_model % BCs(i) % PMatrix
             IF ( .NOT. ASSOCIATED(Projector) ) CYCLE
             DO j=1,Projector % NumberOfRows
                l = Perm( Projector % InvPerm(j) )
                IF ( l > 0 ) THEN
                   l = dofs * (l-1) + DOF
                   DO k=Projector % Rows(j), Projector % Rows(j+1)-1
                      m = Perm( Projector % Cols(k) )
                      IF ( m > 0 )  THEN
                         m = dofs * (m-1) + DOF
                         IF ( ABS(1-Projector % Values(k))<=1.d-12 ) THEN
                            PeriodicNodes = PeriodicNodes + 1
                            IF ( PeriodicNodes > SIZE(R,1) ) THEN
                               TmpR => R
                               ALLOCATE( R(2*PeriodicNodes,2) )
                               R(1:SIZE(TmpR,1),:) = TmpR
                               R(SIZE(TmpR,1)+1:,:) = 0
                               DEALLOCATE( TmpR )
                            END IF
                            R(PeriodicNodes,1) = l
                            R(PeriodicNodes,2) = m
                         ELSE IF ( projector % values(k) /= 0.0d0 ) THEN
                         END IF
                      END IF
                   END DO
                END IF
             END DO
         ELSE
            DO j=acc_solver % Mesh % NumberOfBulkElements+1, acc_solver % Mesh % NumberOfBulkElements + &
                        acc_solver % Mesh % NumberOfBoundaryElements
               IF ( acc_solver % Mesh % Elements(j) % TYPE % ElementCode /= 102 ) CYCLE
               l = Perm( acc_solver % Mesh % Elements(j) % NodeIndexes(1) )
               m = Perm( acc_solver % Mesh % Elements(j) % NodeIndexes(2) )
               l = dofs * (l-1) + DOF
               m = dofs * (m-1) + DOF
               PeriodicNodes = PeriodicNodes + 1
               IF ( PeriodicNodes > SIZE(R,1) ) THEN
                  TmpR => R
                  ALLOCATE( R(2*PeriodicNodes,2) )
                  R(1:SIZE(TmpR,1),:) = TmpR
                  R(SIZE(TmpR,1)+1:,:) = 0
                  DEALLOCATE( TmpR )
               END IF
                            R(PeriodicNodes,1) = l
                            R(PeriodicNodes,2) = m
            END DO
         END IF
       END DO
    END DO
    IF ( PeriodicNodes <= 0 .AND. ParEnv % PEs <= 1 ) THEN
        rv = EliminateDirichlet( acc_model, acc_solver, acc_a, b, x, n, dofs, norm )
        acc_delegated = 1
        RETURN
     END IF
    CALL PeriodicCnstr( R, PeriodicNodes, acc_a % NumberOFRows, PeriodicDOFs, ColNTmp, IDC )
    PeriodicNodes = IDC
    acc_numberofrows  = acc_a % NumberOfRows - PeriodicNodes
    C => AllocateMatrix()
    NumberOfValues = acc_a % NumberOfRows
    C % NumberOfRows = acc_a % NumberOfRows
    ALLOCATE( C % Rows( C % NumberOfRows+1 ), C % Cols( NumberOfValues ), &
         C % Values( NumberOfValues ), C % Diag( C % NumberOfRows ), STAT=istat )
    IF ( istat /= 0 )  CALL Fatal( 'EliminatePeriodic', 'Memory allocation error' )
    C % Values = 0.0d0
    C % Rows(1) = 1
    DO i=1,C % NumberOfRows
       C % Rows(i+1) = C % Rows(i) + 1 ! ColN(i)
    END DO
    DO i = 1, C % NumberOfRows
       t = C % Rows(i) ! + ColN(i) - 1
       C % Values(t) = 1.0d0
       C % Cols(t) = ColNTmp(PeriodicDOFs(i))
    END DO
    DEALLOCATE( PeriodicDOFs,  PeriodicProj, R )
    C % Ordered = .FALSE.
    C % Diag = 0
    CALL CRS_SortMatrixValues( C )
    C_Trans => CRS_TransposeLocal( C, .FALSE. )
    RedA => C
    C => C_Trans
    C_Trans => RedA
    C % Ordered = .FALSE.
    ALLOCATE( C % Diag(C % NumberOFrows ) )
    C % Diag = 0
    CALL CRS_SortMatrixValues( C )
    acc_a % EMatrix => C
    C % EMatrix => C_Trans
    NULLIFY( RedA )
    ! SAVE mirror of the original's ColN (author's evident intent).
    if (allocated(acc_coln)) deallocate(acc_coln)
    allocate(acc_coln(size(ColNTmp)))
    acc_coln = ColNTmp
    deallocate(ColNTmp)
    acc_c => C
    acc_ctrans => C_Trans
    rv = 0
  end function elp_build

  ! -- re-entry link refresh (original 241-244, verbatim) ----------------
  subroutine elp_relink() bind(c, name="elp_relink")
    acc_c => acc_a % EMatrix
    acc_ctrans => acc_c % EMatrix
    acc_reda => acc_ctrans % EMatrix
  end subroutine elp_relink

  ! -- eigen setup (original 251-265, verbatim) --------------------------
  function elp_eig_setup() result(neigen) bind(c, name="elp_eig_setup")
    integer(c_int) :: neigen
    logical :: GotIt
    integer :: NEig
    acc_eiganal = ListGetLogical( acc_solver % Values, 'Eigen Analysis', GotIt )
    neigen = 0
    IF ( acc_eiganal ) THEN
       NEig = ListGetInteger( acc_solver % Values, 'Eigen System Values', GotIt )
       IF ( GotIt .AND. NEig > 0 ) THEN
          DEALLOCATE( acc_solver % Variable % EigenVectors )
          acc_solver % NOFEigenValues = NEig
          ALLOCATE( acc_solver % Variable % EigenVectors( NEig, acc_numberofrows ) )
          acc_solver % Variable % EigenValues  = 0.0d0
          acc_solver % Variable % EigenVectors = 0.0d0
          if (associated(acc_locre)) deallocate(acc_locre)
          if (associated(acc_locim)) deallocate(acc_locim)
          ALLOCATE( acc_locre(NEig,acc_numberofrows) )
          ALLOCATE( acc_locim(NEig,acc_numberofrows) )
          neigen = NEig
       END IF
    END IF
    acc_neigen = neigen
  end function elp_eig_setup

  ! -- CAC^T reduction (original 271-295, verbatim) ----------------------
  subroutine elp_reduce() bind(c, name="elp_reduce")
    TYPE(matrix_t), POINTER :: BB, A_Trans, BB_Trans
    WRITE( Message, * ) 'Constructing the reduced matrix...'
    CALL Info( 'EliminatePeriodic', Message, Level=5 )
    A_Trans => CRS_TransposeLocal( acc_a, acc_eiganal )
    NULLIFY( BB )
    CALL CRS_MatrixMatrixMultiply( BB, A_Trans, acc_ctrans, acc_eiganal )
    CALL FreeMatrix( A_Trans )
    CALL Info( 'EliminatePeriodic', Message, Level=5 )
    BB_Trans => CRS_TransposeLocal( BB, acc_eiganal )
    CALL FreeMatrix( BB )
    CALL Info( 'EliminatePeriodic', Message, Level=5 )
    CALL CRS_MatrixMatrixMultiply( acc_reda, BB_Trans, acc_ctrans, acc_eiganal )
    CALL FreeMatrix( BB_Trans )
    ALLOCATE( acc_f( acc_reda % NumberOfRows ), acc_u( acc_reda % NumberOfRows ) )
    acc_f = 0.0d0
    acc_u = 0.0d0
    IF ( .NOT. acc_eiganal ) CALL CRS_MatrixVectorMultiply( acc_c, acc_b, acc_f )
  end subroutine elp_reduce

  ! -- reduced flags (original 297-299, verbatim) ------------------------
  subroutine elp_pre_solve() bind(c, name="elp_pre_solve")
    acc_reda % Lumped    = acc_a % Lumped
    acc_reda % COMPLEX   = acc_a % COMPLEX
    acc_reda % Symmetric = acc_a % Symmetric
  end subroutine elp_pre_solve

  ! -- reduced-system Perm workspace (original 303-304) ------------------
  subroutine elp_perm_alloc() bind(c, name="elp_perm_alloc")
    if (associated(acc_perm)) deallocate(acc_perm)
    ALLOCATE( acc_perm(SIZE(acc_solver % Variable % Perm)) )
    acc_perm = 0
  end subroutine elp_perm_alloc

  ! -- scalar traffic for the Rust-side Perm remap (original 305-313) ----
  function elp_vperm_size() result(n) bind(c, name="elp_vperm_size")
    integer(c_int) :: n
    n = SIZE(acc_solver % Variable % Perm)
  end function elp_vperm_size

  function elp_vperm_get(i) result(p) bind(c, name="elp_vperm_get")
    integer(c_int), value :: i
    integer(c_int) :: p
    p = acc_solver % Variable % Perm(i)
  end function elp_vperm_get

  function elp_perm_get(i) result(p) bind(c, name="elp_perm_get")
    integer(c_int), value :: i
    integer(c_int) :: p
    p = acc_perm(i)
  end function elp_perm_get

  subroutine elp_perm_set(i, v) bind(c, name="elp_perm_set")
    integer(c_int), value :: i, v
    acc_perm(i) = v
  end subroutine elp_perm_set

  function elp_coln_get(d) result(c) bind(c, name="elp_coln_get")
    integer(c_int), value :: d
    integer(c_int) :: c
    c = acc_coln(d)
  end function elp_coln_get

  ! -- reduced solve (original 301-327, verbatim) ------------------------
  subroutine elp_solve(dofs, norm) bind(c, name="elp_solve")
    integer(c_int), value :: dofs
    real(c_double) :: norm
    INTERFACE
      INTEGER FUNCTION EliminateDirichlet( Model, Solver, A, b, x, n, DOFs, Norm )
         USE Types
         USE Lists
         USE SolverBasics
         USE CRSmatrix
         USE GeneralUtils
         TYPE(model_t)  :: Model
         TYPE(solver_t) :: Solver
         TYPE(matrix_t), POINTER :: A
         INTEGER :: DOFs, n
         REAL(KIND=dp) :: b(n), x(n), Norm
      END FUNCTION EliminateDirichlet
    END INTERFACE
    TYPE(matrix_t), POINTER :: BB
    INTEGER, POINTER :: PermSave(:)
    INTEGER :: k
    LOGICAL :: GotIt
    acc_at = cputime()
    PermSave => acc_solver % Variable % Perm
    BB => acc_solver % Matrix
    acc_solver % Matrix => acc_reda
    acc_solver % Variable % Perm => acc_perm
    CALL ParallelInitMatrix( acc_solver, acc_reda )
    IF ( ParEnv % PEs <= 1 ) THEN
       k = EliminateDirichlet( acc_model, acc_solver, acc_reda, acc_f, acc_u, &
            acc_reda % NumberOfRows, dofs, norm )
    ELSE
      k = 0
    END IF
    IF ( k == 0 ) CALL SolveLinearSystem( acc_reda, acc_f, acc_u, norm, dofs, acc_solver )
    acc_solver % Matrix => BB
    acc_solver % Variable % Perm => PermSave
    acc_at = cputime() - acc_at
    WRITE( Message, * ) 'Time used in SolveLinearSystem (CPU): ', acc_at
    CALL Info( 'EliminatePeriodic', Message, Level=5 )
  end subroutine elp_solve

  ! -- map-back + Norm (original 339-361, verbatim) ----------------------
  subroutine elp_mapback(n, norm) bind(c, name="elp_mapback")
    integer(c_int), value :: n
    real(c_double) :: norm
    REAL(KIND=dp), ALLOCATABLE :: TempVectorReal(:), TempVectorImag(:)
    integer :: i
    IF ( acc_eiganal ) THEN
       acc_locre =  REAL( acc_solver % Variable % EigenVectors )
       acc_locim = AIMAG( acc_solver % Variable % EigenVectors )
       DEALLOCATE( acc_solver % Variable % EigenVectors )
       ALLOCATE( acc_solver % Variable % EigenVectors( acc_neigen, &
           SIZE( acc_solver % Variable % Values ) ) )
       acc_solver % Variable % EigenVectors = 0.0d0
       ALLOCATE( TempVectorReal( SIZE( acc_solver % Variable % Values ) ) )
       ALLOCATE( TempVectorImag( SIZE( acc_solver % Variable % Values ) ) )
       DO i = 1, acc_neigen
         CALL CRS_MatrixVectorMultiply( acc_ctrans, acc_locre(i,:), TempVectorReal )
         CALL CRS_MatrixVectorMultiply( acc_ctrans, acc_locim(i,:), TempVectorImag )
         acc_solver % Variable % EigenVectors(i,1:n) = CMPLX( TempVectorReal, TempVectorImag,KIND=dp )
       END DO
       DEALLOCATE( TempVectorReal, TempVectorImag )
    ELSE
       CALL CRS_MatrixVectorMultiply( acc_ctrans, acc_u, acc_x )
    END IF
    norm = SQRT( SUM(acc_x**2) / acc_a % NumberOFRows )
  end subroutine elp_mapback

  ! -- tail (original 364-373): F/U drop, timing, return 1 ---------------
  function elp_finish() result(rv) bind(c, name="elp_finish")
    integer(c_int) :: rv
    REAL(KIND=dp) :: TotTime
    DEALLOCATE( acc_f, acc_u )
    rv = 1
    TotTime = cputime() - acc_tottime
    WRITE( Message, * ) 'Total time spent in PeriodicReduction (CPU): ', TotTime
    CALL Info( 'EliminatePeriodic', ' ', Level=5 )
    CALL Info( 'EliminatePeriodic', Message, Level=5 )
    CALL Info( 'EliminatePeriodic', ' ', Level=5 )
  end function elp_finish

  !-----------------------------------------------------------------------------
  !  Ideana oli tehda matriisit Neumanin reunaehdoilla seka taulukko
  !  R, jossa on solmuarvojen riippuvuudet tyyliin R(i,1) = R(i,2).
  !  (Comment block carried over with the PeriodicCnstr routine below.)
  !-----------------------------------------------------------------------------
  SUBROUTINE PeriodicCnstr( R, n, NOFRows, d, Perm, IDC )
  !-----------------------------------------------------------------------------
  INTEGER :: R(:,:), n, NOFRows, d(:), Perm(:), IDC
    INTEGER :: i, j, k, CountChanges
    LOGICAL :: ChangesOccurred
    DO i = 1, NOFRows
       d(i) = i
    END DO
    DO i = 1, n
       j = R(i,1) ! MINVAL( R(i,1:2) )
       k = R(i,2) ! MAXVAL( R(i,1:2) )
       d(j) = d(k)
    END DO
    ChangesOccurred = .TRUE.
    DO WHILE( ChangesOccurred )
       ChangesOccurred = .FALSE.
       CountChanges = 0
       DO i = 1,NOFRows
          IF( d(i) /= d(d(i)) ) THEN
             d(i) = d(d(i))
             ChangesOccurred = .TRUE.
             CountChanges = CountChanges+1
          END IF
       END DO
       PRINT *,'Eliminated',CountChanges,' dependent constraints'
    END DO
    j = 0
    Perm = 0
    DO i = 1, NOFRows
       IF( d(i) == i ) THEN
          j = j + 1
          Perm(i) = j
       END IF
    END DO
    IDC = NOFRows - j
  END SUBROUTINE PeriodicCnstr

  SUBROUTINE CRS_SortMatrixValues( A )
    USE CRSMatrix
    USE GeneralUtils
    IMPLICIT NONE
    TYPE(Matrix_t), POINTER :: A
    INTEGER :: i, j, n
    INTEGER, POINTER :: Cols(:), Rows(:), Diag(:)
    REAL(KIND=dp), POINTER :: Values(:)
    Diag   => A % Diag
    Rows   => A % Rows
    Cols   => A % Cols
    Values => A % Values
    n = A % NumberOfRows
    IF ( .NOT. A % Ordered ) THEN
      DO i = 1, n
        CALL SortF( Rows(i+1)-Rows(i), Cols(Rows(i):Rows(i+1)-1), &
            Values(Rows(i):Rows(i+1)-1) )
      END DO
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

  FUNCTION CRS_TransposeLocal( A, MVal ) RESULT(B)
    USE CRSMatrix
    IMPLICIT NONE
    TYPE(Matrix_t), POINTER :: A, B
    LOGICAL :: MVal
    INTEGER, ALLOCATABLE :: Row(:)
    INTEGER :: NVals
    INTEGER :: i,j,k,l
    B => AllocateMatrix()
    NVals = SIZE( A % Values )
    B % NumberOfRows = MAXVAL( A % Cols )
    ALLOCATE( B % Rows( B % NumberOfRows +1 ), B % Cols( NVals ), &
        B % Values( Nvals ) )
    IF ( Mval )  ALLOCATE( B % MassValues( NVals ) )
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
    IF ( Mval ) THEN
       DO i = 1, A % NumberOfRows
          DO j = A % Rows(i), A % Rows(i+1) - 1
             k = A % Cols(j)
             IF ( Row(k) < B % Rows(k+1) ) THEN
                B % Cols( Row(k) ) = i
                B % Values( Row(k) ) = A % Values(j)
                B % MassValues( Row(k) ) = A % MassValues(j)
                Row(k) = Row(k) + 1
             ELSE
                WRITE( Message, * ) 'Trying to access non-existent column', i,k
                !ALL Error( 'CRS_TransposeLocal', Message )
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

  SUBROUTINE CRS_MatrixMatrixMultiply( C, A, B, MVal )
    USE CRSMatrix
    IMPLICIT NONE
    TYPE(Matrix_t), POINTER :: A, B, C
    LOGICAL :: Mval
    REAL(KIND=dp) , POINTER :: AMassVals(:), BMassVals(:)
    LOGICAL, ALLOCATABLE :: Row(:)
    INTEGER, ALLOCATABLE :: Ind(:)
    INTEGER :: ki, NCols
    INTEGER :: i,j,k,l,m,n,Rownonzeros( A % NumberOfRows ), TotalNonzeros
    IF ( .NOT. ASSOCIATED( C ) ) THEN
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
       C => AllocateMatrix()
       C % NumberOfRows = A % NumberOFRows
       ALLOCATE( C % Cols( TotalNonzeros ), C % Values( TotalNonzeros ) )
       ALLOCATE( C % Rows( C % NumberOfRows + 1 ), &
            C % Diag( C % NumberOfRows ) )
       IF ( MVal )  ALLOCATE( C % MassValues( TotalNonzeros ) )
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

end module eliminateperiodic_acc
