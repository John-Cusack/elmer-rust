! T021-B ModelMixedPoisson port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/ModelMixedPoisson.F90, compiled against
! the built .mod files of the target configure (Poisson-in-T020 pattern).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module mirrors the original's own locals (SAVE Nodes per internal
! procedure, SAVE ValueHandles, first-call init flags) with identical
! lifetimes. Order-sensitive reductions (MATMUL, SUM) and whole-array
! updates stay in Fortran, so bit-identity holds by construction wherever
! the call sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module modelmixedpoisson_acc
  use iso_c_binding
  use DefUtils
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(model_t), pointer, save :: acc_model => null()

  ! LocalMatrix context (mirrors its SAVE Nodes and SAVE ValueHandles).
  type(nodes_t), target, save :: lm_nodes
  type(ValueHandle_t), save :: lm_source_h, lm_matpar_h, lm_mattensor_h, lm_conv_h
  logical, save :: lm_inited = .FALSE.

  ! LocalMatrixBC context (mirrors its SAVE Nodes and SAVE ValueHandle).
  type(nodes_t), target, save :: lbc_nodes
  type(ValueHandle_t), save :: lbc_scalar_h
  logical, save :: lbc_inited = .FALSE.

  ! MixedPoisson_post GetFlux context (mirrors its SAVE Nodes).
  type(nodes_t), target, save :: post_nodes
  type(element_t), pointer, save :: post_elem => null()
  logical, save :: post_secondfamily = .FALSE.
  real(kind=dp), allocatable, save :: post_fx(:), post_fy(:), post_fz(:)

contains

  ! -- entry state ----------------------------------------------------
  subroutine mmp_begin(h_solver) bind(c, name="mmp_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine mmp_begin

  subroutine mmp_init0_begin(h_model) bind(c, name="mmp_init0_begin")
    type(c_ptr), value :: h_model
    call c_f_pointer(h_model, acc_model)
  end subroutine mmp_init0_begin

  ! -- MixedPoisson_Init0 (verbatim body; Model/dt/transient unused) --
  subroutine mmp_init0_go() bind(c, name="mmp_init0_go")
    TYPE(ValueList_t), POINTER :: SolverPars
    LOGICAL :: Found, SecondFamily
    INTEGER :: dim
    CHARACTER(LEN=MAX_NAME_LEN) :: csys, VarName

    SolverPars => GetSolverParams()
    SecondFamily = GetLogical(SolverPars, 'Second Kind Basis', Found)

    csys = ListGetString(acc_model % Simulation, 'Coordinate System', Found)
    IF (.NOT. Found) THEN
      IF (.NOT. ListCheckPresent(SolverPars, 'Element')) &
          CALL Fatal('MixedPoisson_Init0', 'The keyword Element should be specified')
    ELSE
      SELECT CASE (csys)
      CASE('cartesian 2d')
        IF (SecondFamily) THEN
          CALL ListAddNewString(SolverPars, "Element", "n:0 e:2 b:1")
        ELSE
          CALL ListAddNewString(SolverPars, "Element", "n:0 e:1 -tri b:1 -quad b:3")
        END IF
      CASE('cartesian 3d')
        IF (SecondFamily) THEN
          CALL ListAddNewString(SolverPars, "Element", &
              "n:0 -tetra b:1 -brick b:25 -quad_face b:4 -tri_face b:3")
        ELSE
          CALL ListAddNewString(SolverPars, "Element", &
              "n:0 -tetra b:1 -brick b:25 -quad_face b:4 -tri_face b:1")
        END IF
      CASE DEFAULT
        IF (.NOT. ListCheckPresent(SolverPars, 'Element')) &
            CALL Fatal('MixedPoisson_Init0', 'The keyword Element should be specified')
      END SELECT
    END IF

    CALL ListAddNewLogical(SolverPars, 'Bubbles in Global System', .TRUE.)

    CALL ListAddNewString(SolverPars,'Potential Variable','mixedpot' )
    VarName = ListGetString(SolverPars,'Potential Variable')

    CALL ListAddString( SolverPars,NextFreeKeyword( &
        'Exported Variable',SolverPars),'-elem '//TRIM(VarName))

    CALL ListAddString( SolverPars,NextFreeKeyword( &
        'Exported Variable',SolverPars), TRIM(VarName))

    CALL ListAddNewString(SolverPars,'Flux Variable','mixedflux' )
    VarName = ListGetString(SolverPars,'Flux Variable')

    CALL ListAddString( SolverPars,NextFreeKeyword( &
        'Exported Variable',SolverPars),'-dofs 3 -dg '//TRIM(VarName))

    CALL ListAddString( SolverPars,NextFreeKeyword( &
        'Exported Variable',SolverPars), TRIM(VarName))

    CALL ListAddNewString( SolverPars,'Element Integration Points', &
        '-tri 3 -quad 9 -tetra 4 -brick 64')

    CALL ListAddLogical( SolverPars,'PostSolver Active',.TRUE.)
  end subroutine mmp_init0_go

  ! -- assembly driver -------------------------------------------------
  subroutine mmp_default_start() bind(c, name="mmp_default_start")
    CALL DefaultStart()
  end subroutine mmp_default_start

  function mmp_coord_dim() result(d) bind(c, name="mmp_coord_dim")
    integer(c_int) :: d
    d = CoordinateSystemDimension()
  end function mmp_coord_dim

  function mmp_mesh_maxdofs() result(n) bind(c, name="mmp_mesh_maxdofs")
    integer(c_int) :: n
    TYPE(Mesh_t), POINTER :: Mesh
    Mesh => GetMesh()
    n = Mesh % MaxElementDOFs
  end function mmp_mesh_maxdofs

  function mmp_nof_active() result(n) bind(c, name="mmp_nof_active")
    integer(c_int) :: n
    n = GetNOFActive()
  end function mmp_nof_active

  subroutine mmp_default_initialize() bind(c, name="mmp_default_initialize")
    CALL DefaultInitialize()
  end subroutine mmp_default_initialize

  function mmp_second_family() result(r) bind(c, name="mmp_second_family")
    integer(c_int) :: r
    LOGICAL :: Found, SecondFamily
    SecondFamily = GetLogical(GetSolverParams(), 'Second Kind Basis', Found)
    r = merge(1, 0, SecondFamily)
  end function mmp_second_family

  subroutine mmp_begin_active(t, n, nd, nb) bind(c, name="mmp_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, nb
    acc_elem => GetActiveElement(t)
    n  = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
    nb = SIZE(acc_elem % BubbleIndexes(:))
  end subroutine mmp_begin_active

  subroutine mmp_lm_reset() bind(c, name="mmp_lm_reset")
    lm_inited = .FALSE.
  end subroutine mmp_lm_reset

  ! -- LocalMatrix (verbatim body over the (ld,ld) SAVE sections) -----
  subroutine mmp_lm_go(ld, stiff, mass, force, n, nd, nb, dim, &
      second_family, need_mass) bind(c, name="mmp_lm_go")
    USE MGDynMaterialUtils, ONLY : GetTensor
    integer(c_int), value :: ld, n, nd, nb, dim, second_family, need_mass
    real(c_double) :: stiff(ld,ld), mass(ld,ld), force(ld)
    LOGICAL :: SecondFamily, NeedMass
    INTEGER, PARAMETER :: MaxFaceBasisDim = 48
    TYPE(GaussIntegrationPoints_t) :: IP
    LOGICAL :: Stat, Found
    INTEGER :: t, i, j, p, q, np, mu_rank
    REAL(KIND=dp), POINTER :: mu_tensor(:,:)
    REAL(KIND=dp) :: a, f, v(3), tmp(dim), mu
    REAL(KIND=dp) :: FaceBasis(MaxFaceBasisDim,3), DivFaceBasis(MaxFaceBasisDim)
    REAL(KIND=dp) :: Basis(nd), DetJ, s
    REAL(KIND=dp) :: StiffW(nd,nd), MassW(nd,nd), ForceW(nd)

    SecondFamily = (second_family /= 0)
    NeedMass = (need_mass /= 0)

    IF (.NOT. lm_inited) THEN
      CALL ListInitElementKeyword( lm_source_h, 'Body Force','Source Field')

      CALL ListInitElementKeyword( lm_matpar_h,'Material','Material Parameter',DefRValue=1.0_dp)
      CALL ListInitElementKeyword( lm_mattensor_h,'Material','Material Tensor')

      CALL ListInitElementKeyword( lm_conv_h,'Material','Convection Velocity',InitVec3D=.TRUE.)
      lm_inited = .TRUE.
    END IF

    CALL GetElementNodes( lm_nodes )

    ! Whole-array zeroing, exactly as the original's `Stiff = 0` etc.;
    ! the (1:nd,1:nd) accumulation below then matches element for element.
    stiff = 0.0d0
    mass = 0.0d0
    force = 0.0d0

    StiffW = 0.0d0
    MassW = 0.0d0
    ForceW = 0.0d0

    IF( .FALSE. ) THEN
      IP = GaussPointsAdapt( acc_elem, PReferenceElement = .TRUE. )
    ELSE
      SELECT CASE( GetElementFamily(acc_elem) )
      CASE(3)
        IP = GaussPointsTriangle(3, PReferenceElement=.TRUE.)
      CASE(4)
        IP = GaussPointsQuad(9)
      CASE(5)
        IP = GaussPointsTetra(4, PReferenceElement=.TRUE.)
      CASE(8)
        IP = GaussPointsBrick(64)
      CASE DEFAULT
        CALL Fatal('MixedPoisson', 'A non-supported element type')
      END SELECT
    END IF

    np = n * acc_solver % Def_Dofs(GetElementFamily(acc_elem), acc_elem % BodyId, 1)

    DO t=1,IP % n
      stat = FaceElementInfo(acc_elem, lm_nodes, IP % U(t), IP % V(t), &
          IP % W(t), detF=detJ, Basis=Basis, FBasis=FaceBasis, &
          DivFBasis=DivFaceBasis, BDM=SecondFamily, ApplyPiolaTransform=.TRUE.)

      s = detJ * IP % s(t)

      IF (np > 0) THEN
        DO p = 1,n
          DO q = 1,n
            StiffW(p,q) = StiffW(p,q) + Basis(p) * Basis(q) * s
          END DO

          DO q = nd,nd
            StiffW(p,q) = StiffW(p,q) - Basis(p) * 1.0d0 * s
          END DO
        END DO
      END IF

      mu = ListGetElementReal( lm_mattensor_h, Basis, acc_elem, Found, GaussPoint = t, &
          Rdim = mu_rank, Rtensor = mu_tensor)
      IF( Found ) THEN
        IF( mu_rank == 0 ) THEN
          DO p = 1,nd-np-1
            i = np + p
            DO q = 1,nd-np-1
              j = np + q
              StiffW(i,j) = StiffW(i,j) + mu * &
                  SUM( FaceBasis(q,1:dim) * FaceBasis(p,1:dim) ) * s
            END DO
          END DO
        ELSE IF (mu_rank == 1 ) THEN
          DO p = 1,nd-np-1
            i = np + p
            DO q = 1,nd-np-1
              j = np + q
              tmp = SUM(mu_tensor(1:dim,1)*FaceBasis(q,1:dim))
              StiffW(i,j) = StiffW(i,j) + &
                  SUM( tmp(1:dim) * FaceBasis(p,1:dim) ) * s
            END DO
          END DO
        ELSE IF (mu_rank == 2 ) THEN
          DO p = 1,nd-np-1
            i = np + p
            DO q = 1,nd-np-1
              j = np + q
              tmp = MATMUL(mu_tensor(1:dim,1:dim),FaceBasis(q,1:dim))
              StiffW(i,j) = StiffW(i,j) + &
                  SUM( tmp(1:dim) * FaceBasis(p,1:dim) ) * s
            END DO
          END DO
        END IF
      ELSE
        a = ListGetElementReal( lm_matpar_h, Basis, acc_elem, Found, GaussPoint = t )
        DO p = 1,nd-np-1
          i = np + p
          DO q = 1,nd-np-1
            j = np + q
            StiffW(i,j) = StiffW(i,j) + (1.0_dp / a ) * &
                SUM( FaceBasis(q,1:dim) * FaceBasis(p,1:dim) ) * s
          END DO
        END DO
      END IF

      DO p = 1,nd-np-1
        i = np + p
        DO q = nd,nd
          StiffW(i,q) = StiffW(i,q) + 1.0d0 * DivFaceBasis(p) * s
        END DO
      END DO

      DO p = nd,nd
        DO q = 1,nd-np-1
          j = np + q
          StiffW(p,j) = StiffW(p,j) + 1.0d0 * DivFaceBasis(q) * s
        END DO
      END DO

      v = ListGetElementReal3D( lm_conv_h, Basis, acc_elem, Found, GaussPoint = t )
      IF ( Found ) THEN
        DO p = nd,nd
          DO q = 1,nd-np-1
            j = np + q
            StiffW(p,j) = StiffW(p,j) - SUM(FaceBasis(q,1:dim) * v(1:dim)) * 1.0d0 * s
          END DO
        END DO
      END IF

      f = ListGetElementReal( lm_source_h, Basis, acc_elem, Found, GaussPoint = t )
      IF ( Found ) THEN
        DO p = nd,nd
          ForceW(p) = ForceW(p) - f * 1.0d0 * s
        END DO
      END IF

      IF (NeedMass) THEN
        MassW(nd,nd) = MassW(nd,nd) - s
      END IF

    END DO

    stiff(1:nd,1:nd) = StiffW(1:nd,1:nd)
    mass(1:nd,1:nd) = MassW(1:nd,1:nd)
    force(1:nd) = ForceW(1:nd)
  end subroutine mmp_lm_go

  subroutine mmp_first_order_time(ld, stiff, mass, force) &
      bind(c, name="mmp_first_order_time")
    integer(c_int), value :: ld
    real(c_double) :: stiff(ld,ld), mass(ld,ld), force(ld)
    CALL Default1stOrderTime(mass(1:ld,1:ld), stiff(1:ld,1:ld), force(1:ld))
  end subroutine mmp_first_order_time

  subroutine mmp_default_update(ld, stiff, force) bind(c, name="mmp_default_update")
    integer(c_int), value :: ld
    real(c_double) :: stiff(ld,ld), force(ld)
    CALL DefaultUpdateEquations(stiff(1:ld,1:ld), force(1:ld))
  end subroutine mmp_default_update

  subroutine mmp_default_finish_bulk() bind(c, name="mmp_default_finish_bulk")
    CALL DefaultFinishBulkAssembly()
  end subroutine mmp_default_finish_bulk

  ! -- boundary assembly ----------------------------------------------
  function mmp_nof_bnd() result(n) bind(c, name="mmp_nof_bnd")
    integer(c_int) :: n
    n = GetNOFBoundaryElements()
  end function mmp_nof_bnd

  subroutine mmp_begin_bnd(t, n, nd, active) bind(c, name="mmp_begin_bnd")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, active
    acc_elem => GetBoundaryElement(t)
    IF (ActiveBoundaryElement()) THEN
      active = 1
      n  = GetElementNOFNodes(acc_elem)
      nd = GetElementNOFDOFs(acc_elem)
    ELSE
      active = 0
    END IF
  end subroutine mmp_begin_bnd

  subroutine mmp_bnd_cond(n, cond) bind(c, name="mmp_bnd_cond")
    integer(c_int), value :: n
    real(c_double) :: cond
    LOGICAL :: Found
    cond = SUM(GetReal(GetBC(), GetVarName(acc_solver % Variable)//' Condition', Found))/n
  end subroutine mmp_bnd_cond

  subroutine mmp_lbc_reset() bind(c, name="mmp_lbc_reset")
    lbc_inited = .FALSE.
  end subroutine mmp_lbc_reset

  ! -- LocalMatrixBC (verbatim body; Mesh via GetMesh) ----------------
  subroutine mmp_lbc_go(n, nd, second_family) bind(c, name="mmp_lbc_go")
    USE ElementDescription, ONLY : PickActiveFace
    integer(c_int), value :: n, nd, second_family
    LOGICAL :: SecondFamily
    TYPE(Mesh_t), POINTER :: parMesh
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Element_t), POINTER :: Parent, Face
    TYPE(GaussIntegrationPoints_t) :: IP
    LOGICAL :: ReverseSign(6), stat, AssembleForce, OrientationsMatch
    LOGICAL :: RevertIndices
    INTEGER, POINTER :: FaceMap(:,:)
    INTEGER, TARGET :: TetraFaceMap(4,3), BrickFaceMap(6,4)
    INTEGER :: FDofMap(6,4)
    INTEGER :: OriginalIndices(4)
    INTEGER :: j, k, l, t, np, p, ActiveFaceId, Family, FDOFs
    INTEGER :: ParentFamily
    REAL(KIND=dp) :: Force(nd), Basis(n), TraceBasis(nd), WorkTrace(nd)
    REAL(KIND=dp) :: detJ, s, u, v, w, g

    SecondFamily = (second_family /= 0)
    parMesh => GetMesh()

    Family = GetElementFamily(acc_elem)
    IF (Family == 1) RETURN

    BC => GetBC()
    IF (.NOT. ASSOCIATED(BC)) RETURN

    IF (.NOT. lbc_inited) THEN
      CALL ListInitElementKeyword(lbc_scalar_h, 'Boundary Condition', &
          'Scalar Field')
      lbc_inited = .TRUE.
    END IF
    IF (lbc_scalar_h % NotPresentAnywhere) RETURN

    Parent => acc_elem % BoundaryInfo % Left
    IF (.NOT. ASSOCIATED(Parent)) THEN
      Parent => acc_elem % BoundaryInfo % Right
    END IF
    IF (.NOT. ASSOCIATED(Parent)) RETURN
    ParentFamily = GetElementFamily(Parent)

    CALL PickActiveFace(parMesh, Parent, acc_elem, Face, ActiveFaceId)
    IF (ActiveFaceId == 0) RETURN
    FDOFs = Face % BDOFs
    IF (FDOFs < 1) RETURN

    CALL FaceElementOrientation(Parent, ReverseSign, ActiveFaceId)

    RevertIndices = .FALSE.
    IF ((ParentFamily == 3 .OR. ParentFamily == 5) .AND. SecondFamily) THEN
      SELECT CASE(Family)
      CASE(2)
        FaceMap => GetEdgeMap(GetElementFamily(Parent))
        IF (ReverseSign(ActiveFaceId)) THEN
          OrientationsMatch = ( acc_elem % NodeIndexes(1) == Parent % NodeIndexes(FaceMap(ActiveFaceId,2)) )
        ELSE
          OrientationsMatch = ( acc_elem % NodeIndexes(1) == Parent % NodeIndexes(FaceMap(ActiveFaceId,1)) )
        END IF

      CASE(3)
        IF (FDOFs /= 3) CALL Fatal('ModelMixedPoisson', '3-DOF faces expected')
        TetraFaceMap(1,:) = (/ 2, 1, 3 /)
        TetraFaceMap(2,:) = (/ 1, 2, 4 /)
        TetraFaceMap(3,:) = (/ 2, 3, 4 /)
        TetraFaceMap(4,:) = (/ 3, 1, 4 /)

        CALL FaceElementBasisOrdering(Parent, FDofMap, ActiveFaceId)

        IF (ANY(acc_elem % NodeIndexes(1:3) /= Parent % NodeIndexes(TetraFaceMap(ActiveFaceId,1:3)))) THEN
          OriginalIndices(1:3) = acc_elem % NodeIndexes(1:3)
          acc_elem % NodeIndexes(1:3) = Parent % NodeIndexes(TetraFaceMap(ActiveFaceId,1:3))
          RevertIndices = .TRUE.
        END IF

      END SELECT
    ELSE
      SELECT CASE(GetElementFamily(Parent))
      CASE(8)
        IF (FDOFs /= 4) CALL Fatal('ModelMixedPoisson', '4-DOF faces expected')

        BrickFaceMap(1,:) = (/ 2, 1, 4, 3 /)
        BrickFaceMap(2,:) = (/ 5, 6, 7, 8 /)
        BrickFaceMap(3,:) = (/ 1, 2, 6, 5 /)
        BrickFaceMap(4,:) = (/ 2, 3, 7, 6 /)
        BrickFaceMap(5,:) = (/ 3, 4, 8, 7 /)
        BrickFaceMap(6,:) = (/ 4, 1, 5, 8 /)

        CALL FaceElementBasisOrdering(Parent, FDofMap, ActiveFaceId)

        IF (ANY(acc_elem % NodeIndexes(1:4) /= Parent % NodeIndexes(BrickFaceMap(ActiveFaceId,1:4)))) THEN
          OriginalIndices(1:4) = acc_elem % NodeIndexes(1:4)
          acc_elem % NodeIndexes(1:4) = Parent % NodeIndexes(BrickFaceMap(ActiveFaceId,1:4))
          RevertIndices = .TRUE.
        END IF
      END SELECT
    END IF

    np = n * acc_solver % Def_Dofs(GetElementFamily(Parent), Parent % BodyId, 1)

    IF (ReverseSign(ActiveFaceId)) THEN
      s = -1.0d0
    ELSE
      s = 1.0d0
    END IF

    CALL GetElementNodes(lbc_nodes)
    IF (Family == 3) THEN
      IP = GaussPointsTriangle(3, PReferenceElement=.TRUE.)
    ELSE
      IP = GaussPoints(acc_elem)
    END IF

    Force = 0.0d0
    DO t=1,IP % n
      SELECT CASE(Family)
      CASE(2)
        stat = ElementInfo(acc_elem, lbc_nodes, IP % U(t), IP % V(t), &
            IP % W(t), DetJ, Basis)
        IF (SecondFamily) THEN
          u = IP % U(t)
          IF (OrientationsMatch) THEN
            TraceBasis(1) = s * 0.5d0 * (1.0d0 - sqrt(3.0d0)*u)
            TraceBasis(2) = s * 0.5d0 * (1.0d0 + sqrt(3.0d0)*u)
          ELSE
            TraceBasis(1) = s * 0.5d0 * (1.0d0 + sqrt(3.0d0)*u)
            TraceBasis(2) = s * 0.5d0 * (1.0d0 - sqrt(3.0d0)*u)
          END IF
        ELSE
          TraceBasis(1) = s * 0.5d0
        END IF
      CASE(3)
        u = IP % U(t)
        v = IP % V(t)
        DO k=1,n
          Basis(k) = TriangleNodalPBasis(k, u, v)
        END DO

        IF (SecondFamily) THEN
          WorkTrace(1) = -1.0d0*(-8 - 12*u + 4*Sqrt(3.0d0)*v)/12.0d0
          WorkTrace(2) = -1.0d0*(u + (-8 + 4*Sqrt(3.0d0)*v)/12.0d0)
          WorkTrace(3) = -1.0d0*(4.0d0 - 8*Sqrt(3.0d0)*v)/12.0d0
          DO j=1,FDOFs
            k = FDofMap(ActiveFaceId,j)
            TraceBasis(j) = s * WorkTrace(k)
          END DO
        ELSE
          TraceBasis(1) = s / SQRT(3.0d0)
        END IF
      CASE(4)
        stat = ElementInfo(acc_elem, lbc_nodes, IP % U(t), IP % V(t), &
            IP % W(t), DetJ, Basis)

        DO j=1,FDOFs
          k = FDofMap(ActiveFaceId,j)
          TraceBasis(j) = s * Basis(k)
        END DO
      END SELECT

      w = IP % s(t)
      g = ListGetElementReal(lbc_scalar_h, Basis, acc_elem, AssembleForce)

      IF (AssembleForce) THEN
        DO p = 1,nd-np
          j = np + p
          Force(j) = Force(j) + g * TraceBasis(p) * w
        END DO
      END IF
    END DO

    IF (AssembleForce) CALL DefaultUpdateForce(Force)

    IF (RevertIndices) THEN
      SELECT CASE(Family)
      CASE(3)
        acc_elem % NodeIndexes(1:3) = OriginalIndices(1:3)
      CASE(4)
        acc_elem % NodeIndexes(1:4) = OriginalIndices(1:4)
      END SELECT
    END IF
  end subroutine mmp_lbc_go

  subroutine mmp_default_finish_bnd() bind(c, name="mmp_default_finish_bnd")
    CALL DefaultFinishBoundaryAssembly()
  end subroutine mmp_default_finish_bnd

  subroutine mmp_default_finish() bind(c, name="mmp_default_finish")
    CALL DefaultFinishAssembly()
  end subroutine mmp_default_finish

  subroutine mmp_default_dirichlet() bind(c, name="mmp_default_dirichlet")
    CALL DefaultDirichletBCs()
  end subroutine mmp_default_dirichlet

  function mmp_default_solve() result(norm) bind(c, name="mmp_default_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function mmp_default_solve

  subroutine mmp_fatal_alloc() bind(c, name="mmp_fatal_alloc")
    CALL Fatal('MixedPoisson', 'Memory allocation error.')
  end subroutine mmp_fatal_alloc

  ! -- MixedPoisson_post (verbatim body) ------------------------------
  subroutine mmp_post_go() bind(c, name="mmp_post_go")
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Mesh_t), POINTER :: Mesh
    CHARACTER(LEN=MAX_NAME_LEN) :: VarName
    LOGICAL :: Found
    TYPE(Variable_t), POINTER :: pVar, Var, fVar
    TYPE(Element_t), POINTER :: Element
    INTEGER, ALLOCATABLE :: Indexes(:)
    INTEGER :: dim, active, t, n, nd, nb, np
    REAL(KIND=dp) :: val

    Params => GetSolverParams()

    post_secondfamily = GetLogical(Params, 'Second Kind Basis', Found)

    Mesh => GetMesh()
    Var => acc_solver % Variable
    dim = CoordinateSystemDimension()

    n = acc_solver % Mesh % MaxElementDOFs
    ALLOCATE( Indexes(n) )
    ALLOCATE( post_fx(n), post_fy(n), post_fz(n))

    VarName = ListGetString(Params,'Potential Variable')
    pVar => VariableGet( Mesh % Variables, VarName )

    VarName = ListGetString(Params,'Flux Variable')
    fVar => VariableGet( Mesh % Variables, VarName )

    active = GetNOFActive()

    DO t=1,active
      Element => GetActiveElement(t)
      post_elem => Element
      n  = GetElementNOFNodes()
      nd = GetElementDOFs( Indexes, Element )
      nb = SIZE(Element % BubbleIndexes(:))

      np = n * acc_solver % Def_Dofs(GetElementFamily(Element), Element % BodyId, 1)

      val = Var % Values( Var % Perm(Indexes(nd)) )

      pVar % Values(pVar % Perm(Element % ElementIndex)) = val

      CALL post_getflux( n,nd-np-1,Var % Values(Var % Perm(Indexes(np+1:nd-1))) )

      fVar % Values(3*(fVar % Perm(Element % dGIndexes(1:n))-1)+1) = post_fx(1:n)
      fVar % Values(3*(fVar % Perm(Element % dGIndexes(1:n))-1)+2) = post_fy(1:n)
      fVar % Values(3*(fVar % Perm(Element % dGIndexes(1:n))-1)+3) = post_fz(1:n)

    END DO

    DEALLOCATE( Indexes, post_fx, post_fy, post_fz )
  end subroutine mmp_post_go

  ! -- GetFlux (verbatim body over module state) ----------------------
  SUBROUTINE post_getflux( n, nval, vals )
    INTEGER :: n, nval
    REAL(KIND=dp) :: vals(:)
    REAL(KIND=dp) :: MASS(n,n), FORCE_x(n), FORCE_y(n), FORCE_z(n)
    INTEGER, PARAMETER :: MaxFaceBasisDim = 48
    TYPE(GaussIntegrationPoints_t) :: IP
    LOGICAL :: Stat, Found
    INTEGER :: t, i, j, p, q, np
    REAL(KIND=dp) :: FaceBasis(MaxFaceBasisDim,3), DivFaceBasis(MaxFaceBasisDim)
    REAL(KIND=dp) :: Basis(n), DetJ, s

    SELECT CASE( GetElementFamily(post_elem) )
    CASE(3)
      IP = GaussPointsTriangle(3, PReferenceElement=.TRUE.)
    CASE(4)
      IP = GaussPointsQuad(9)
    CASE(5)
      IP = GaussPointsTetra(4, PReferenceElement=.TRUE.)
    CASE(8)
      IP = GaussPointsBrick(64)
    CASE DEFAULT
      CALL Fatal('MixedPoisson', 'A non-supported element type')
    END SELECT

    CALL GetElementNodes(post_nodes)

    MASS = 0._dp
    FORCE_x = 0._dp
    FORCE_y = 0._dp
    FORCE_z = 0._dp
    np = n * acc_solver % Def_Dofs(GetElementFamily(post_elem), post_elem % BodyId, 1)

    DO t=1,IP % n
      stat = FaceElementInfo(post_elem, post_nodes, IP % U(t), IP % V(t), &
          IP % W(t), detF=detJ, Basis=Basis, FBasis=FaceBasis, &
          DivFBasis=DivFaceBasis, BDM=post_secondfamily, ApplyPiolaTransform=.TRUE.)

      s = IP % s(t) * DetJ

      DO i=1,n
        DO j=1,n
          MASS(i,j) = MASS(i,j) + s * Basis(i)*Basis(j)
        END DO
        FORCE_x(i) = FORCE_x(i) + s * SUM(FaceBasis(1:nval,1) * vals(1:nval)) * Basis(i)
        FORCE_y(i) = FORCE_y(i) + s * SUM(FaceBasis(1:nval,2) * vals(1:nval)) * Basis(i)
        FORCE_z(i) = FORCE_z(i) + s * SUM(FaceBasis(1:nval,3) * vals(1:nval)) * Basis(i)
      END DO
    END DO

    CALL InvertMatrix(MASS,n)
    post_fx(1:n) = MATMUL(MASS, FORCE_x)
    post_fy(1:n) = MATMUL(MASS, FORCE_y)
    post_fz(1:n) = MATMUL(MASS, FORCE_z)
  END SUBROUTINE post_getflux

end module modelmixedpoisson_acc
