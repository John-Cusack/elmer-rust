! T022-A EMWaveSolver port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/EMWaveSolver.F90, compiled against the
! built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, list-directed timing
! reports — stays inside Fortran, so the Rust side is pure C calls. The
! wrapper module mirrors EMWaveSolver's own locals (current element,
! Nodes SAVE state per internal procedure, ValueHandle SAVE state, the
! SAVE-by-initialization ConstantBulkInUse flag; the entry SAVE work
! arrays are Rust-side) with identical lifetimes.
!
! Order-sensitive floating-point work (the SUM/MATMUL Galerkin assembly,
! the whole-array fills, the strided eigen copies, the LU solves) stays
! in Fortran wrappers transcribed statement for statement, so
! bit-identity holds by construction wherever the call sequence matches.
! The commented-out ListSetCoefficients block of the original `_Init0`
! has no observable effect and is documented as skipped rather than
! replicated.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module emwavesolver_acc
  use iso_c_binding
  use types
  use defutils
  use integration
  use matrixassembly
  use eleminfo
  use messages
  use meshbasics, only: CalculateBodyAverage
  implicit none

  ! -- shared entry state ------------------------------------------------
  type(solver_t), pointer, save :: acc_solver => null()
  type(model_t), pointer, save :: acc_model => null()
  type(valuelist_t), pointer, save :: acc_params => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(mesh_t), pointer, save :: acc_mesh => null()
  type(solver_t), pointer, save :: acc_psolver => null()
  real(c_double), save :: acc_eps0 = 0.0_c_double
  real(c_double), save :: acc_mu0 = 0.0_c_double
  ! Default-kind integers (converted at the BIND(C) boundary) so internal
  ! calls match the library dummies exactly.
  integer, save :: acc_timeorder = 2
  integer, save :: acc_edgedeg = 1
  logical, save :: acc_piola = .false.
  ! DoBulkAssembly SAVE-by-initialization mirror: the original declares
  ! `ConstantBulkInUse = .FALSE.` (implicit SAVE), reads it on entry and
  ! refreshes it at the end, so it persists across invocations.
  logical, save :: acc_cbulkinuse = .false.

  ! LocalMatrix context (mirrors its SAVE Nodes and SAVE ValueHandles).
  type(nodes_t), target, save :: acc_nodes_lm
  type(gaussintegrationpoints_t), save :: acc_ip_lm
  type(valuehandle_t), save :: acc_cdh(3), acc_muh_lm, acc_epsh_lm, acc_condh

  ! LocalMatrixBC context (mirrors its SAVE Nodes, SAVE ValueHandles and
  ! SAVE Visited flag, which is set-but-never-read in the original).
  type(nodes_t), target, save :: acc_nodes_bc
  type(gaussintegrationpoints_t), save :: acc_ip_bc
  type(valuehandle_t), save :: acc_blh(3), acc_damph, acc_muh_bc, acc_epsh_bc
  logical, save :: acc_visited = .false.

  ! EMWaveCalcFields setup state (mirrors the entry locals that the
  ! CONTAINS procedures host-associate).
  type(variable_t), pointer, save :: acc_pvar => null()
  type(solver_t), pointer, save :: acc_cf_psolver => null()
  type(mesh_t), pointer, save :: acc_cfmesh => null()
  type(variable_t), pointer, save :: acc_ef => null()
  type(variable_t), pointer, save :: acc_def => null()
  type(variable_t), pointer, save :: acc_ddef => null()
  type(variable_t), pointer, save :: acc_ef_e => null()
  type(variable_t), pointer, save :: acc_def_e => null()
  type(variable_t), pointer, save :: acc_ddef_e => null()
  real(c_double), pointer, save :: acc_fsave(:) => null()
  logical, save :: acc_eigen = .false.
  logical, save :: acc_anytimeder = .false.
  logical, save :: acc_nodal = .false.
  logical, save :: acc_elemental = .false.
  logical, save :: acc_doave = .false.
  logical, save :: acc_cf_cbulkmatrix = .false.
  logical, save :: acc_cf_cbulkinuse = .false.
  integer, save :: acc_nofeigen = 1
  integer, save :: acc_cdofs = 1
  integer, save :: acc_cf_dofs = 0
  integer, save :: acc_ieigen = 1

contains

  ! -- EMWaveSolver_Init0 --------------------------------------------------
  subroutine emw_i0_begin(h_model, h_solver) bind(c, name="emw_i0_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine emw_i0_begin

  subroutine emw_i0_body() bind(c, name="emw_i0_body")
    type(valuelist_t), pointer :: solverparams
    logical :: found, secondorder, secondkind, piolaversion
    real(c_double) :: mu0, eps0
    solverparams => GetSolverParams()
    if (.not. ListCheckPresent(solverparams, "Element")) then
      call EdgeElementStyle(solverparams, piolaversion, secondkind, &
          secondorder, check = .true.)
      if (secondorder) then
        call ListAddString(solverparams, "Element", &
            "n:0 e:2 -tri b:2 -quad b:4 -brick b:6 -pyramid b:3 -prism b:2" // &
            " -quad_face b:4 -tri_face b:2")
      else if (secondkind) then
        call ListAddString(solverparams, "Element", "n:0 e:2")
      else if (piolaversion) then
        call ListAddString(solverparams, "Element", &
            "n:0 e:1 -quad b:2 -brick b:3 -quad_face b:2")
      else
        call ListAddString(solverparams, "Element", "n:0 e:1")
      end if
    end if
    call ListAddNewLogical(solverparams, 'Hcurl Basis', .true.)
    if (ListGetLogical(solverparams, 'Constant Bulk Matrix', found) .or. &
        ListGetLogical(solverparams, 'Eigen Analysis', found)) then
      call ListAddNewLogical(solverparams, 'Use Global Mass Matrix', .true.)
    end if
    call ListAddNewLogical(solverparams, 'Variable Output', .false.)
    call ListAddNewString(solverparams, 'Variable', 'E')
    call ListAddNewLogical(solverparams, 'Linear System Complex', .false.)
    call ListAddNewInteger(solverparams, 'Time derivative order', 2)
    eps0 = GetConstReal(acc_model % constants, 'Permittivity of Vacuum', found)
    if (.not. found) call Fatal('EMWaveSolver_Init0', &
        '> Permittivity of Vacuum < is required')
    mu0 = GetConstReal(acc_model % constants, 'Permeability of Vacuum', found)
    if (.not. found) call Fatal('EMWaveSolver_Init0', &
        '> Permeability of Vacuum < is required')
  end subroutine emw_i0_body

  ! -- EMWaveSolver entry --------------------------------------------------
  subroutine emw_begin(h_model, h_solver) bind(c, name="emw_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine emw_begin

  subroutine emw_prologue() bind(c, name="emw_prologue")
    type(valuelist_t), pointer :: solverparams
    integer, pointer :: perm(:)
    integer :: n, dofs, basideg
    logical :: piola
    solverparams => GetSolverParams()
    acc_params => solverparams
    call EdgeElementStyle(solverparams, piola, basisdegree = basideg)
    acc_piola = piola
    acc_edgedeg = basideg
    if (CoordinateSystemDimension() == 2 .and. .not. piola) then
      call Fatal('EMWaveSolver', 'A 2D model needs Use Piola Transform = True')
    end if
    acc_timeorder = ListGetInteger(solverparams, 'Time derivative order')
    dofs = acc_solver % variable % dofs
    acc_eps0 = GetConstReal(acc_model % constants, 'Permittivity of Vacuum')
    acc_mu0 = GetConstReal(acc_model % constants, 'Permeability of Vacuum')
    acc_mesh => GetMesh()
    n = acc_mesh % numberofnodes
    perm => acc_solver % variable % perm
    acc_psolver => acc_solver
    if (dofs /= 1) call Fatal('EMWaveSolver', &
        'Invalid variable size:' // I2S(dofs))
  end subroutine emw_prologue

  function emw_maxdofs() result(n) bind(c, name="emw_maxdofs")
    integer(c_int) :: n
    n = acc_mesh % maxelementdofs
  end function emw_maxdofs

  function emw_noitmax() result(m) bind(c, name="emw_noitmax")
    integer(c_int) :: m
    logical :: found
    m = GetInteger(acc_params, 'Nonlinear System Max Iterations', found)
    if (.not. found) m = 1
  end function emw_noitmax

  function emw_edgebasis() result(r) bind(c, name="emw_edgebasis")
    integer(c_int) :: r
    logical :: found
    if (.not. ListCheckPresent(acc_params, 'Linear System Refactorize') &
        .and. GetLogical(acc_params, 'Edge Basis', found)) then
      r = 1
    else
      r = 0
    end if
  end function emw_edgebasis

  subroutine emw_default_start() bind(c, name="emw_default_start")
    call DefaultStart()
  end subroutine emw_default_start

  subroutine emw_bulk(stiff, damp, mass, force, n) bind(c, name="emw_bulk")
    integer(c_int), value :: n
    real(c_double) :: stiff(n, n), damp(n, n), mass(n, n), force(n)
    integer :: active, t, nd, nn
    logical :: found
    if (acc_cbulkinuse) then
      call DefaultInitialize(UseConstantBulk = acc_cbulkinuse)
      return
    end if
    call DefaultInitialize()
    active = GetNOFActive()
    do t = 1, active
      acc_elem => GetActiveElement(t)
      nn = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      call LocalMatrix(mass, damp, stiff, force, &
          acc_elem, nn, nd, acc_piola, t == 1)
      if (acc_timeorder == 2) then
        call Default2ndOrderTime(mass, damp, stiff, force)
      else if (acc_timeorder == 1) then
        call Default1stOrderTime(damp, stiff, force)
      end if
      call DefaultUpdateEquations(stiff, force)
    end do
    call DefaultFinishBulkAssembly()
    acc_cbulkinuse = ListGetLogical(acc_params, 'Constant Bulk Matrix', found)
  end subroutine emw_bulk

  subroutine emw_bnd(stiff, damp, mass, force, n) bind(c, name="emw_bnd")
    integer(c_int), value :: n
    real(c_double) :: stiff(n, n), damp(n, n), mass(n, n), force(n)
    integer :: active, t, nd, nn, k
    logical :: init_handles
    type(element_t), pointer :: dummy_element
    type(valuelist_t), pointer :: bc
    init_handles = .true.
    active = GetNOFBoundaryElements()
    do t = 1, active
      acc_elem => GetBoundaryElement(t)
      bc => GetBC()
      if (.not. associated(bc)) cycle
      select case (GetElementFamily())
      case (1)
        cycle
      case (2)
        k = GetBoundaryEdgeIndex(acc_elem, 1)
        acc_elem => acc_mesh % edges(k)
      case (3, 4)
        k = GetBoundaryFaceIndex(acc_elem)
        acc_elem => acc_mesh % faces(k)
      end select
      if (.not. ActiveBoundaryElement(acc_elem)) cycle
      nd = GetElementNOFDOFs(acc_elem)
      nn = GetElementNOFNodes(acc_elem)
      dummy_element => SetCurrentElement(acc_elem)
      call LocalMatrixBC(mass, damp, stiff, force, &
          acc_elem, nn, nd, acc_piola, init_handles)
      if (acc_timeorder == 2) then
        call Default2ndOrderTimeR(mass, damp, stiff, force(1:nd), &
            UElement = acc_elem)
      else if (acc_timeorder == 1) then
        call Default1stOrderTimeR(damp, stiff, force(1:nd), &
            UElement = acc_elem)
      end if
      call DefaultUpdateEquationsR(stiff, force(1:nd), UElement = acc_elem)
      init_handles = .false.
    end do
    call DefaultFinishBoundaryAssembly()
    call DefaultFinishAssembly()
    call DefaultDirichletBCs()
  end subroutine emw_bnd

  function emw_solve() result(norm) bind(c, name="emw_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function emw_solve

  function emw_converged() result(r) bind(c, name="emw_converged")
    integer(c_int) :: r
    r = merge(1, 0, DefaultConverged())
  end function emw_converged

  subroutine emw_norefactor() bind(c, name="emw_norefactor")
    call ListAddLogical(acc_params, 'Linear System Refactorize', .false.)
  end subroutine emw_norefactor

  subroutine emw_rmrefactor() bind(c, name="emw_rmrefactor")
    call ListRemove(acc_params, 'Linear System Refactorize')
  end subroutine emw_rmrefactor

  subroutine emw_finish() bind(c, name="emw_finish")
    call DefaultFinish()
  end subroutine emw_finish

  subroutine emw_done() bind(c, name="emw_done")
    call Info('EMWaveSolver', 'All done', level = 10)
  end subroutine emw_done

  subroutine emw_info_start() bind(c, name="emw_info_start")
    call Info('EMWaveSolver', 'Solving electromagnetic waves in time', level = 5)
  end subroutine emw_info_start

  subroutine emw_fatal_alloc() bind(c, name="emw_fatal_alloc")
    call Fatal('EMWaveSolver', 'Memory allocation error.')
  end subroutine emw_fatal_alloc

  subroutine LocalMatrix(MASS, DAMP, STIFF, FORCE, &
      Element, n, nd, PiolaVersion, InitHandles)
    real(c_double) :: STIFF(:, :), DAMP(:, :), FORCE(:), MASS(:, :)
    integer :: n, nd
    type(element_t), pointer :: Element
    logical :: PiolaVersion, InitHandles
    real(c_double) :: DetJ, weight, eps, mu, muinv, L(3), cond
    real(c_double) :: Basis(n), dBasisdx(n, 3), WBasis(nd, 3), RotWBasis(nd, 3)
    logical :: Stat
    integer :: t, i, j
    type(GaussIntegrationPoints_t) :: IP
    type(Nodes_t), save :: Nodes
    type(ValueHandle_t) :: CD_h(3), Mu_h, Eps_h, Cond_h
    logical :: AllocationsDone = .false.
    save Cd_h, Mu_h, Eps_h, Cond_h, AllocationsDone
    logical :: Found
    if (InitHandles) then
      call ListInitElementKeyword(Cd_h(1), 'Body Force', 'Current Density Rate 1')
      call ListInitElementKeyword(Cd_h(2), 'Body Force', 'Current Density Rate 2')
      call ListInitElementKeyword(Cd_h(3), 'Body Force', 'Current Density Rate 3')
      call ListInitElementKeyword(Mu_h, 'Material', 'Relative Permeability')
      call ListInitElementKeyword(Eps_h, 'Material', 'Relative Permittivity')
      call ListInitElementKeyword(Cond_h, 'Material', 'Electric Conductivity')
    end if
    call GetElementNodes(Nodes, Element)
    STIFF = 0.0_dp
    DAMP = 0.0_dp
    FORCE = 0.0_dp
    MASS = 0.0_dp
    IP = GaussPoints(Element, EdgeBasis = .true., PReferenceElement = PiolaVersion)
    do t = 1, IP % n
      stat = ElementInfo(Element, Nodes, IP % u(t), IP % v(t), IP % w(t), detJ, &
          Basis, dBasisdx, EdgeBasis = Wbasis, RotBasis = RotWBasis, &
          USolver = acc_psolver)
      weight = detJ * IP % s(t)
      eps = acc_eps0 * ListGetElementReal(Eps_h, Basis, Element, Found)
      mu = acc_mu0 * ListGetElementReal(Mu_h, Basis, Element, Found)
      muinv = 1.0_dp / mu
      cond = ListGetElementReal(Cond_h, Basis, Element, Found)
      L(1) = ListGetElementReal(CD_h(1), Basis, Element, Found)
      L(2) = ListGetElementReal(CD_h(2), Basis, Element, Found)
      L(3) = ListGetElementReal(CD_h(3), Basis, Element, Found)
      do i = 1, nd
        FORCE(i) = FORCE(i) - sum(L * WBasis(i, :)) * weight
        do j = 1, nd
          STIFF(i, j) = STIFF(i, j) + muinv * &
              sum(RotWBasis(i, :) * RotWBasis(j, :)) * weight
          MASS(i, j) = MASS(i, j) + &
              eps * sum(WBasis(j, :) * WBasis(i, :)) * weight
          DAMP(i, j) = DAMP(i, j) + &
              cond * sum(WBasis(j, :) * WBasis(i, :)) * weight
        end do
      end do
    end do
  end subroutine LocalMatrix

  subroutine LocalMatrixBC(MASS, DAMP, STIFF, FORCE, Element, n, nd, &
      PiolaVersion, InitHandles)
    real(c_double) :: MASS(:, :), DAMP(:, :), STIFF(:, :), FORCE(:)
    integer :: n, nd
    type(element_t), pointer :: Element
    logical :: PiolaVersion, InitHandles
    real(c_double) :: DetJ, Tem(n)
    real(c_double) :: B, L(3), muinv, mu, weight, tanWBasis(3)
    real(c_double) :: Basis(n), WBasis(nd, 3), RotWBasis(nd, 3), dBasisdx(n, 3)
    logical :: Stat
    type(GaussIntegrationPoints_t) :: IP
    integer :: t, i, j, p, q
    type(Nodes_t), save :: Nodes
    type(ValueHandle_t) :: BL_h(3), Damp_h, Mu_h, Eps_h
    logical :: Visited = .false., GotTem
    type(element_t), pointer :: Parent
    logical :: AllocationsDone = .false.
    save Visited, BL_h, Damp_h, Mu_h, Eps_h
    logical :: Found
    type(valuelist_t), pointer :: bc
    if (InitHandles) then
      call ListInitElementKeyword(Bl_h(1), 'Boundary Condition', 'Magnetic Boundary Load 1')
      call ListInitElementKeyword(Bl_h(2), 'Boundary Condition', 'Magnetic Boundary Load 2')
      call ListInitElementKeyword(Bl_h(3), 'Boundary Condition', 'Magnetic Boundary Load 3')
      call ListInitElementKeyword(Damp_h, 'Boundary Condition', 'Electric Damping Coefficient')
      call ListInitElementKeyword(Mu_h, 'Material', 'Relative Permeability')
      call ListInitElementKeyword(Eps_h, 'Material', 'Relative Permittivity')
    end if
    bc => GetBC()
    Parent => GetBulkElementAtBoundary(Element)
    TEM(1:n) = GetReal(bc, 'TEM Potential', GotTem)
    call GetElementNodes(Nodes, Element)
    DAMP = 0.0_dp
    STIFF = 0.0_dp
    FORCE = 0.0_dp
    MASS = 0.0_dp
    IP = GaussPoints(Element, EdgeBasis = .true., PReferenceElement = PiolaVersion)
    do t = 1, IP % n
      if (GetElementFamily(Element) == 2) then
        stat = EdgeElementInfo(Element, Nodes, IP % u(t), IP % v(t), IP % w(t), &
            detF = detJ, Basis = Basis, EdgeBasis = Wbasis, &
            dBasisdx = dBasisdx, BasisDegree = acc_edgedeg, &
            ApplyPiolaTransform = .true.)
      else
        stat = ElementInfo(Element, Nodes, IP % u(t), IP % v(t), IP % w(t), detJ, &
            Basis, dBasisdx, EdgeBasis = Wbasis, RotBasis = RotWBasis, &
            USolver = acc_psolver)
      end if
      weight = detJ * IP % s(t)
      mu = acc_mu0 * ListGetElementReal(Mu_h, Basis, Parent)
      muinv = 1.0_dp / mu
      L(1) = ListGetElementReal(Bl_h(1), Basis, Element, Found)
      L(2) = ListGetElementReal(Bl_h(2), Basis, Element, Found)
      L(3) = ListGetElementReal(Bl_h(3), Basis, Element, Found)
      if (GotTem) then
        L = L + matmul(Tem(1:n), dBasisdx(1:n, 1:3))
      end if
      do i = 1, nd
        FORCE(i) = FORCE(i) - muinv * sum(L(1:3) * WBasis(i, 1:3)) * weight
      end do
      B = ListGetElementReal(Damp_h, Basis, Element, Found)
      if (Found) then
        do i = 1, nd
          do j = 1, nd
            DAMP(i, j) = DAMP(i, j) + muinv * B * &
                sum(WBasis(i, :) * WBasis(j, :)) * weight
          end do
        end do
      end if
    end do
  end subroutine LocalMatrixBC

  ! -- EMWaveCalcFields_Init0 ------------------------------------------------
  subroutine emw_cfi0_begin(h_model, h_solver) bind(c, name="emw_cfi0_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine emw_cfi0_begin

  subroutine emw_cfi0_body() bind(c, name="emw_cfi0_body")
    character(len=MAX_NAME_LEN) :: sname, pname
    logical :: found, elementfields, eigenanalysis
    integer, pointer :: active(:)
    integer :: mysolver, i, j, n, m, soln
    type(valuelist_t), pointer :: solverparams
    type(solver_t), pointer :: solvers(:), psolver
    solverparams => GetSolverParams()
    soln = 0
    do i = 1, acc_model % numberofsolvers
      sname = GetString(acc_model % solvers(i) % values, 'Procedure', found)
      j = index(sname, 'EMWaveSolver')
      if (j > 0) then
        soln = i
        exit
      end if
    end do
    if (soln == 0) then
      call Fatal('EMWaveCalcFields_Init0', 'Cannot locate the primary solver: ' // I2S(soln))
    else
      call Info('EMWaveCalcFields_Init0', 'The primary solver index is: ' // I2S(soln), level = 12)
      call ListAddInteger(solverparams, 'Primary Solver Index', soln)
    end if
    eigenanalysis = ListGetLogical(acc_model % solvers(soln) % values, &
        'Eigen Analysis', found)
    if (eigenanalysis) then
      call ListAddNewLogical(solverparams, 'Eigen Analysis', .true.)
      call ListAddNewLogical(solverparams, 'Constant Bulk Matrix', .true.)
    end if
    if (GetLogical(solverparams, 'Discontinuous Galerkin', found)) return
    elementfields = GetLogical(solverparams, 'Calculate Elemental Fields', found)
    if (found .and. .not. elementfields) return
    psolver => acc_solver
    do mysolver = 1, acc_model % numberofsolvers
      if (associated(psolver, acc_model % solvers(mysolver))) exit
    end do
    n = acc_model % numberofsolvers
    do i = 1, acc_model % numberofequations
      active => ListGetIntegerArray(acc_model % equations(i) % values, &
          'Active Solvers', found)
      m = size(active)
      if (any(active == mysolver)) &
          call ListAddIntegerArray(acc_model % equations(i) % values, &
          'Active Solvers', m + 1, [active, n + 1])
    end do
    allocate (solvers(n + 1))
    solvers(1:n) = acc_model % solvers
    solverparams => null()
    call ListAddLogical(solverparams, 'Discontinuous Galerkin', .true.)
    solvers(n + 1) % dg = .true.
    solvers(n + 1) % values => solverparams
    solvers(n + 1) % procedure = C_NULL_FUNPTR
    solvers(n + 1) % activeelements => null()
    call ListAddString(solverparams, 'Exec Solver', 'never')
    call ListAddLogical(solverparams, 'No Matrix', .true.)
    call ListAddString(solverparams, 'Equation', 'never')
    call ListAddString(solverparams, 'Procedure', &
        'AllocateSolver AllocateSolver', .false.)
    call ListAddString(solverparams, 'Variable', '-nooutput cf_dummy')
    pname = ListGetString(acc_model % solvers(soln) % values, 'Mesh', found)
    if (found) then
      call ListAddString(solverparams, 'Mesh', pname)
    end if
    call ListAddString(solverparams, &
        NextFreeKeyword('Exported Variable', solverparams), &
        "Elfield E[Elfield E:3]")
    if (.not. eigenanalysis) then
      if (GetLogical(solverparams, 'Calculate Electric field derivatives', found)) then
        call ListAddString(solverparams, &
            NextFreeKeyword('Exported Variable', solverparams), &
            "dEdt E[dEdt E:3]")
        call ListAddString(solverparams, &
            NextFreeKeyword('Exported Variable', solverparams), &
            "ddEddt E[ddEddt E:3]")
      end if
    end if
    deallocate (acc_model % solvers)
    acc_model % solvers => solvers
    acc_model % numberofsolvers = n + 1
  end subroutine emw_cfi0_body

  ! -- EMWaveCalcFields_Init ---------------------------------------------------
  subroutine emw_cfi_begin(h_solver) bind(c, name="emw_cfi_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine emw_cfi_begin

  subroutine emw_cfi_body() bind(c, name="emw_cfi_body")
    logical :: found, nodalfields, eigenanalysis
    type(valuelist_t), pointer :: solverparams
    solverparams => GetSolverParams()
    call ListAddString(solverparams, 'Variable', '-nooutput hr_dummy')
    call ListAddLogical(solverparams, 'Linear System refactorize', .false.)
    call ListAddNewLogical(solverparams, 'Skip Compute Nonlinear Change', .true.)
    nodalfields = GetLogical(solverparams, 'Calculate Nodal Fields', found)
    if (found .and. .not. nodalfields) return
    eigenanalysis = ListGetLogical(solverparams, 'Eigen Analysis', found)
    call ListAddString(solverparams, &
        NextFreeKeyword('Exported Variable', solverparams), &
        "Elfield[Elfield:3]")
    if (.not. eigenanalysis) then
      if (GetLogical(solverparams, 'Calculate Electric field derivatives', found)) then
        call ListAddString(solverparams, &
            NextFreeKeyword('Exported Variable', solverparams), &
            "dEdt[dEdt:3]")
        call ListAddString(solverparams, &
            NextFreeKeyword('Exported Variable', solverparams), &
            "ddEddt[ddEddt:3]")
      end if
    end if
  end subroutine emw_cfi_body

  ! -- EMWaveCalcFields entry --------------------------------------------------
  subroutine emw_cf_begin(h_model, h_solver) bind(c, name="emw_cf_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine emw_cf_begin

  subroutine emw_cf_setup() bind(c, name="emw_cf_setup")
    type(valuelist_t), pointer :: solverparams
    character(len=MAX_NAME_LEN) :: pname
    logical :: found
    integer :: i, j, n, dofs, vdofs
    call Info('EMWaveCalcFields', 'Computing postprocessing fields')
    solverparams => GetSolverParams()
    acc_params => solverparams
    acc_eps0 = GetConstReal(acc_model % constants, 'Permittivity of Vacuum')
    acc_mu0 = GetConstReal(acc_model % constants, 'Permeability of Vacuum')
    i = ListGetInteger(solverparams, 'Primary Solver Index', found)
    if (i == 0) then
      call Fatal('EMWaveCalcFields', 'We should know > Primary Solver Index <')
    end if
    acc_cf_psolver => acc_model % solvers(i)
    acc_pvar => acc_cf_psolver % variable
    pname = getVarName(acc_pvar)
    call Info('EMWaveCalcFields', 'Name of potential variable: ' // trim(pname), level = 10)
    vdofs = acc_pvar % dofs
    if (vdofs /= 1) then
      call Fatal('EMWaveCalcFields', 'Primary variable should have 1 dofs: ' // I2S(vdofs))
    end if
    dofs = 3
    call EdgeElementStyle(acc_cf_psolver % values, acc_piola)
    if (acc_piola) call Info('EMWaveCalcFields', &
        'Using Piola transformed finite elements', level = 5)
    acc_nofeigen = 1
    acc_eigen = GetLogical(acc_cf_psolver % values, 'Eigen Analysis', found)
    if (acc_eigen) then
      acc_nofeigen = size(acc_cf_psolver % variable % eigenvalues)
      call Info('EMWaveCalcFields', &
          'Computing fields for ' // I2S(acc_nofeigen) // ' eigen vectors', level = 7)
    end if
    acc_cfmesh => GetMesh()
    acc_ef => VariableGet(acc_cfmesh % variables, 'Elfield')
    acc_ef_e => VariableGet(acc_cfmesh % variables, 'Elfield E')
    acc_cdofs = 1
    if (acc_eigen) then
      do i = 1, 2
        if (i == 1) then
          acc_def => acc_ef
        else if (i == 2) then
          acc_def => acc_ef_e
        end if
        if (associated(acc_def)) then
          n = size(acc_def % values)
          if (.not. associated(acc_def % eigenvectors)) then
            allocate (acc_def % eigenvectors(acc_nofeigen, n))
            acc_def % eigenvectors = 0.0_dp
            allocate (acc_def % eigenvalues(acc_nofeigen))
            acc_def % eigenvalues = acc_pvar % eigenvalues
          end if
        end if
      end do
      acc_def => null()
      acc_def_e => null()
      acc_ddef => null()
      acc_ddef_e => null()
      acc_anytimeder = .false.
      acc_cdofs = 2
    else
      acc_def => VariableGet(acc_cfmesh % variables, 'dEdt')
      acc_def_e => VariableGet(acc_cfmesh % variables, 'dEdt E')
      acc_ddef => VariableGet(acc_cfmesh % variables, 'ddEddt')
      acc_ddef_e => VariableGet(acc_cfmesh % variables, 'ddEddt E')
    end if
    acc_anytimeder = associated(acc_def) .or. associated(acc_ddef) .or. &
        associated(acc_def_e) .or. associated(acc_ddef_e)
    i = 0
    if (associated(acc_ef)) i = i + acc_cdofs * 3
    if (associated(acc_def)) i = i + 3
    if (associated(acc_ddef)) i = i + 3
    acc_nodal = (i > 0)
    j = 0
    if (associated(acc_ef_e)) j = j + acc_cdofs * 3
    if (associated(acc_def_e)) j = j + 3
    if (associated(acc_ddef_e)) j = j + 3
    acc_elemental = (j > 0)
    dofs = max(i, j)
    acc_cf_dofs = dofs
    n = acc_cfmesh % maxelementdofs
    acc_cf_cbulkmatrix = ListGetLogical(solverparams, 'Constant Bulk Matrix', found)
    acc_cf_cbulkinuse = associated(acc_solver % matrix % bulkvalues) .and. &
        acc_cf_cbulkmatrix
    if (acc_nodal) then
      call DefaultInitialize(UseConstantBulk = acc_cf_cbulkinuse)
    else
      call DefaultInitialize()
    end if
  end subroutine emw_cf_setup

  subroutine emw_cf_counts(nofeigen, dofs, nmax, nrhs, ncols) &
      bind(c, name="emw_cf_counts")
    integer(c_int) :: nofeigen, dofs, nmax, nrhs, ncols
    integer :: i
    nofeigen = acc_nofeigen
    dofs = acc_cf_dofs
    nmax = acc_cfmesh % maxelementdofs
    if (acc_nodal) then
      nrhs = size(acc_solver % matrix % rhs)
      i = 0
      if (associated(acc_ef)) i = i + acc_cdofs * 3
      if (associated(acc_def)) i = i + 3
      if (associated(acc_ddef)) i = i + 3
      ncols = i
    else
      nrhs = 0
      ncols = 0
    end if
  end subroutine emw_cf_counts

  function emw_cf_nactive() result(n) bind(c, name="emw_cf_nactive")
    integer(c_int) :: n
    n = GetNOFActive()
  end function emw_cf_nactive

  subroutine emw_cf_eigenhead(ieigen) bind(c, name="emw_cf_eigenhead")
    integer(c_int), value :: ieigen
    acc_ieigen = ieigen
    if (acc_eigen) then
      call Info('EMwaveCalcFields', 'Computing fields for eigen vector: ' // &
          I2S(acc_ieigen), level = 10)
    end if
  end subroutine emw_cf_eigenhead

  subroutine emw_cf_elem(t, ieigen, mass, force, gf, pivot, n0, d0, nr, nc) &
      bind(c, name="emw_cf_elem")
    integer(c_int), value :: t, ieigen, n0, d0, nr, nc
    real(c_double) :: mass(n0, n0), force(n0, d0)
    real(c_double), target :: gf(nr, nc)
    integer(c_int) :: pivot(n0)
    integer :: n, nd, l, dofcount
    logical :: erroneous
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs(uSolver = acc_cf_psolver)
    call LocalAssembly(t == 1 .and. ieigen == 1, &
        acc_elemental .or. .not. acc_cf_cbulkinuse, mass, force, n0, d0, &
        nd, int(ieigen))
    if (acc_nodal) then
      if (.not. acc_cf_cbulkinuse) then
        call DefaultUpdateEquations(mass, force(:, 1))
      end if
      acc_fsave => acc_solver % matrix % rhs
      do l = 1, acc_cf_dofs
        acc_solver % matrix % rhs => gf(:, l)
        call DefaultUpdateForce(force(:, l))
      end do
      acc_solver % matrix % rhs => acc_fsave
    end if
    if (acc_elemental) then
      dofcount = 0
      call LUdecomp(mass, n, pivot, erroneous)
      if (erroneous) call Fatal('EMWaveCalcFields', 'LU-decomposition fails')
      call LocalSol(acc_ef_e, acc_cdofs * 3, n, mass, force, pivot, dofcount)
      call LocalSol(acc_def_e, 3, n, mass, force, pivot, dofcount)
      call LocalSol(acc_ddef_e, 3, n, mass, force, pivot, dofcount)
    end if
  end subroutine emw_cf_elem

  subroutine emw_cf_post(ieigen, mass, force, gf, n0, d0, nr, nc) &
      bind(c, name="emw_cf_post")
    integer(c_int), value :: ieigen, n0, d0, nr, nc
    real(c_double) :: mass(n0, n0), force(n0, d0)
    real(c_double), target :: gf(nr, nc)
    integer :: dofcount
    logical :: found
    acc_doave = GetLogical(acc_params, 'Average Within Materials', found)
    if (acc_doave) then
      if (GetLogical(acc_params, 'Discontinuous Galerkin', found)) then
        if (.not. acc_cf_cbulkinuse) then
          force = 0.0_c_double
          call AddLocalFaceTerms(mass, force(:, 1))
        end if
      end if
    end if
    if (acc_nodal .or. acc_doave) then
      if (acc_cf_cbulkmatrix .and. .not. acc_cf_cbulkinuse) then
        call Info('EMWaveCalcFields', 'Saving the system matrix', level = 6)
        call CopyBulkMatrix(acc_solver % matrix, BulkRHS = .false.)
      end if
      acc_fsave => null()
      if (associated(acc_solver % matrix)) acc_fsave => acc_solver % matrix % rhs
      dofcount = 0
      call GlobalSol(acc_ef, acc_cdofs * 3, gf, dofcount, acc_ef_e)
      call GlobalSol(acc_def, 3, gf, dofcount, acc_def_e)
      call GlobalSol(acc_ddef, 3, gf, dofcount, acc_ddef_e)
      if (associated(acc_fsave)) acc_solver % matrix % rhs => acc_fsave
    end if
    acc_cf_cbulkinuse = associated(acc_solver % matrix % bulkvalues)
  end subroutine emw_cf_post

  subroutine emw_cf_done() bind(c, name="emw_cf_done")
    call Info('EMCalcFields', 'Done computing postprocessed fields!', level = 10)
  end subroutine emw_cf_done

  subroutine emw_cf_fatal_alloc() bind(c, name="emw_cf_fatal_alloc")
    call Fatal('EMWaveCalcFields', 'Memory allocation error.')
  end subroutine emw_cf_fatal_alloc

  subroutine LocalAssembly(InitHandles, IntegrateMass, MASS, FORCE, n0, d0, nd, ieigen)
    logical :: InitHandles, IntegrateMass
    integer :: n0, d0, nd, ieigen
    real(c_double) :: MASS(n0, n0), FORCE(n0, d0)
    type(GaussIntegrationPoints_t) :: IP
    type(Nodes_t), save :: Nodes
    real(c_double), allocatable :: WBasis(:, :), SOL(:), dsol(:), ddsol(:)
    real(c_double), allocatable :: RotWBasis(:, :), Basis(:), dBasisdx(:, :)
    integer, allocatable :: Indexes(:)
    real(c_double) :: detJ, s, u, v, w
    real(c_double) :: EF_ip(3), dEF_ip(3), ddEF_ip(3)
    logical :: AllocationsDone = .false.
    type(ValueHandle_t) :: Mu_h, Eps_h, Cd_h(3)
    integer :: i, j, k, n, p, q
    logical :: stat
    save Cd_h, Mu_h, Eps_h, sol, dsol, ddsol, &
        Indexes, WBasis, RotWBasis, Basis, dBasisdx
    if (.not. AllocationsDone) then
      n = acc_cfmesh % maxelementdofs
      allocate (SOL(n), dsol(n), ddsol(n), &
          Indexes(n), WBasis(n, 3), RotWBasis(n, 3), Basis(n), dBasisdx(n, 3))
      sol = 0.0_dp
      dsol = 0.0_dp
      ddsol = 0.0_dp
      AllocationsDone = .true.
    end if
    if (InitHandles) then
      call ListInitElementKeyword(Mu_h, 'Material', 'Relative Permeability')
      call ListInitElementKeyword(Eps_h, 'Material', 'Relative Permittivity')
      call ListInitElementKeyword(Cd_h(1), 'Body Force', 'Current Density 1')
      call ListInitElementKeyword(Cd_h(2), 'Body Force', 'Current Density 2')
      call ListInitElementKeyword(Cd_h(3), 'Body Force', 'Current Density 3')
    end if
    call GetElementNodes(Nodes)
    n = GetElementDOFs(Indexes, acc_elem, acc_cf_psolver)
    sol(1:n) = 0.0_dp
    dsol(1:n) = 0.0_dp
    ddsol(1:n) = 0.0_dp
    if (acc_eigen) then
      do i = 1, n
        j = acc_pvar % perm(Indexes(i))
        if (j > 0) then
          sol(i) = real(acc_pvar % eigenvectors(ieigen, j))
          dsol(i) = aimag(acc_pvar % eigenvectors(ieigen, j))
        end if
      end do
    else
      do i = 1, n
        j = acc_pvar % perm(Indexes(i))
        if (j > 0) then
          sol(i) = acc_pvar % values(j)
          if (acc_anytimeder) then
            dsol(i) = acc_pvar % prevvalues(j, 1)
            ddsol(i) = acc_pvar % prevvalues(j, 2)
          end if
        end if
      end do
    end if
    IP = GaussPoints(acc_elem, EdgeBasis = .true., PReferenceElement = acc_piola)
    MASS = 0._dp
    FORCE = 0._dp
    do j = 1, IP % n
      u = IP % U(j)
      v = IP % V(j)
      w = IP % W(j)
      stat = ElementInfo(acc_elem, Nodes, u, v, w, detJ, Basis, dBasisdx, &
          EdgeBasis = Wbasis, RotBasis = RotWBasis, USolver = acc_cf_psolver)
      EF_ip = matmul(sol(1:nd), WBasis(1:nd, :))
      if (acc_anytimeder .or. acc_eigen) then
        dEF_ip = matmul(dsol(1:nd), WBasis(1:nd, :))
        ddEF_ip = matmul(ddsol(1:nd), WBasis(1:nd, :))
      end if
      s = IP % s(j) * detJ
      if (IntegrateMass) then
        do p = 1, n
          do q = 1, n
            MASS(p, q) = MASS(p, q) + s * Basis(p) * Basis(q)
          end do
        end do
      end if
      do p = 1, n
        k = 0
        if (associated(acc_ef) .or. associated(acc_ef_e)) then
          FORCE(p, k + 1:k + 3) = FORCE(p, k + 1:k + 3) + s * (EF_ip) * Basis(p)
          k = k + 3
        end if
        if (associated(acc_def) .or. associated(acc_def_e) .or. acc_eigen) then
          FORCE(p, k + 1:k + 3) = FORCE(p, k + 1:k + 3) + s * (dEF_ip) * Basis(p)
          k = k + 3
        end if
        if (associated(acc_ddef) .or. associated(acc_ddef_e)) then
          FORCE(p, k + 1:k + 3) = FORCE(p, k + 1:k + 3) + s * (ddEF_ip) * Basis(p)
          k = k + 3
        end if
      end do
    end do
  end subroutine LocalAssembly

  subroutine GlobalSol(Var, m, b, dofs, EL_Var)
    real(c_double), target CONTIG :: b(:, :)
    integer :: m, dofs
    type(variable_t), pointer :: Var
    type(variable_t), pointer, optional :: EL_Var
    real(c_double) :: Norm
    integer :: i
    if (present(EL_Var)) then
      if (associated(El_Var)) then
        El_Var % DgAveraged = .false.
        if (acc_doave) then
          call Info('EMWaveSolver', 'Averaging for field: ' // trim(El_Var % Name), level = 10)
          call CalculateBodyAverage(acc_cfmesh, El_Var, .false.)
        end if
        if (.not. (associated(var) .and. acc_nodal)) then
          dofs = dofs + m
          return
        end if
      end if
    end if
    if (.not. associated(Var)) return
    call Info('EMWaveSolver', 'Solving for field: ' // trim(Var % Name), level = 6)
    do i = 1, m
      dofs = dofs + 1
      acc_solver % matrix % rhs => b(:, dofs)
      acc_solver % variable % values = 0.0_dp
      Norm = DefaultSolve()
      if (acc_eigen) then
        if (i <= 3) then
          Var % EigenVectors(acc_ieigen, i::m / 2) = acc_solver % variable % values
        else
          Var % EIgenVectors(acc_ieigen, i - 3::m / 2) = &
              cmplx(real(Var % EIgenVectors(acc_ieigen, i - 3::m / 2)), &
              acc_solver % variable % values, kind = dp)
        end if
      else
        var % values(i::m) = acc_solver % variable % values
      end if
    end do
  end subroutine GlobalSol

  subroutine LocalSol(Var, m, n, A, b, pivot, dofs)
    type(variable_t), pointer :: Var
    integer :: pivot(:), m, n, dofs
    real(c_double) :: b(:, :), A(:, :)
    integer :: ind(n), i
    real(c_double) :: x(n)
    if (.not. associated(var)) return
    ind = Var % dofs * (Var % Perm(acc_elem % DGIndexes(1:n)) - 1)
    do i = 1, m
      dofs = dofs + 1
      x = b(1:n, dofs)
      call LUSolve(n, A, x, pivot)
      if (acc_eigen) then
        if (i <= 3) then
          Var % EigenVectors(acc_ieigen, ind(1:n) + i) = x
        else
          Var % EIgenVectors(acc_ieigen, ind(1:n) + i - 3) = &
              cmplx(real(Var % EigenVectors(acc_ieigen, ind(1:n) + i - 3)), x, kind = dp)
        end if
      else
        Var % Values(ind(1:n) + i) = x
      end if
    end do
  end subroutine LocalSol

  subroutine AddLocalFaceTerms(STIFF, FORCE)
    real(c_double) :: STIFF(:, :), FORCE(:)
    type(element_t), pointer :: P1, P2, Face, Faces(:)
    integer :: t, n, n1, n2, NumberOfFaces, dim
    dim = CoordinateSystemDimension()
    if (dim == 2) then
      Faces => acc_solver % mesh % edges
      NumberOfFaces = acc_solver % mesh % numberofedges
    else
      Faces => acc_solver % mesh % faces
      NumberOfFaces = acc_solver % mesh % numberoffaces
    end if
    do t = 1, NumberOfFaces
      Face => Faces(t)
      if (.not. ActiveBoundaryElement(Face)) cycle
      P1 => Face % BoundaryInfo % Left
      P2 => Face % BoundaryInfo % Right
      if (associated(P2) .and. associated(P1)) then
        if (.not. associated(GetMaterial(P1), GetMaterial(P2))) cycle
        n = GetElementNOFNodes(Face)
        n1 = GetElementNOFNodes(P1)
        n2 = GetElementNOFNodes(P2)
        call LocalJumps(STIFF, Face, n, P1, n1, P2, n2)
        call DefaultUpdateEquations(STIFF, FORCE, Face)
      end if
    end do
  end subroutine AddLocalFaceTerms

  subroutine LocalJumps(STIFF, Face, n, P1, n1, P2, n2)
    implicit none
    real(c_double) :: STIFF(:, :)
    integer :: n, n1, n2
    type(element_t), target :: Face, P1, P2
    real(c_double) :: FaceBasis(n), P1Basis(n1), P2Basis(n2)
    real(c_double) :: Jump(n1 + n2), detJ, U, V, W, S
    logical :: Stat
    integer :: p, q, t
    type(GaussIntegrationPoints_t) :: IntegStuff
    type(Nodes_t) :: FaceNodes, P1Nodes, P2Nodes
    save FaceNodes, P1Nodes, P2Nodes
    STIFF = 0._dp
    call GetElementNodes(FaceNodes, Face)
    call GetElementNodes(P1Nodes, P1)
    call GetElementNodes(P2Nodes, P2)
    IntegStuff = GaussPoints(Face)
    do t = 1, IntegStuff % n
      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)
      S = IntegStuff % s(t)
      stat = ElementInfo(Face, FaceNodes, U, V, W, detJ, FaceBasis)
      S = S * detJ
      call GetParentUVW(Face, n, P1, n1, U, V, W, FaceBasis)
      stat = ElementInfo(P1, P1Nodes, U, V, W, detJ, P1Basis)
      call GetParentUVW(Face, n, P2, n2, U, V, W, FaceBasis)
      stat = ElementInfo(P2, P2Nodes, U, V, W, detJ, P2Basis)
      Jump(1:n1) = P1Basis(1:n1)
      Jump(n1 + 1:n1 + n2) = -P2Basis(1:n2)
      do p = 1, n1 + n2
        do q = 1, n1 + n2
          STIFF(p, q) = STIFF(p, q) + s * Jump(q) * Jump(p)
        end do
      end do
    end do
  end subroutine LocalJumps

end module emwavesolver_acc
