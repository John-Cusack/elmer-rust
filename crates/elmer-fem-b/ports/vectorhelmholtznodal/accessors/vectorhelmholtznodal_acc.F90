! T021-B VectorHelmholtzNodal port: statement-level BIND(C) wrappers for
! the Elmer calls in fem/src/modules/VectorHelmholtzNodal.F90, compiled
! against the built .mod files of the target configure (T020 Poisson
! pattern, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape arrays). All Fortran-object
! traffic — TYPE(...) handles, OPTIONAL/POINTER dummies, assumed-shape
! descriptors, hidden lengths, array-pointer results — stays inside
! Fortran, so the Rust side is pure C calls. The wrapper module mirrors
! the original's SAVE state with identical lifetimes:
! - the entry SAVE `Proj` matrix pointer;
! - the `LocalMatrix` internal-procedure SAVE (Nodes, the four
!   ValueHandles, Basis/dBasisdx/STIFF/FORCE allocatables);
! - the `LocalMatrixBC` internal-procedure SAVE (Nodes, the nine
!   ValueHandles, STIFF/FORCE allocatables, GotPort — declared without
!   an initializer exactly as in the original);
! - the per-element locals that the split across BIND(C) calls forces to
!   persist between a head call and its IP calls (Normal, NT,
!   GoodConductor/Absorb flags, Parent, the boundary ValuesList, the IP
!   object, the trailing Found, kappa) — set by the head wrapper and
!   consumed by the IP wrappers before the next head call, exactly the
!   original lifetimes (Found included: the Parent-guarded muinv fetch
!   may read the head's trailing Found on the first IP pass, and later
!   passes read the previous pass's trailing Found, as in the original
!   single-procedure variable).
! The entry-level automatics (Norm, host mirrors of dim/dofs/compn/
! flags/coefficients) live here as module state, set during the entry
! prologue and consumed within the same invocation; Norm is zeroed at
! entry to mirror its per-invocation automatic lifetime (every read
! follows the per-component DefaultSolve write, as in the original).
!
! Order-sensitive reductions (MATMUL), formatted I/O (PRINT), and
! whole-array/strided updates stay in these wrappers, so bit-identity
! holds by construction wherever the call sequence matches. The `im`
! unit-imaginary parameter comes from VectorHelmholtzUtils, exactly as
! in the original unit.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module vhn_acc
  use iso_c_binding
  use types
  use defutils
  use meshbasics, only : FollowCurvedBoundary
  use crsmatrix, only : CRS_TransposeMatrixVectorMultiply
  use vectorhelmholtzutils
  implicit none

  character(*), parameter :: caller = 'VectorHelmholtzNodal'

  ! Entry handles (re-pointed on every invocation).
  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(mesh_t), pointer, save :: acc_mesh => null()
  type(valuelist_t), pointer, save :: acc_params => null()
  type(variable_t), pointer, save :: acc_ei => null()
  type(variable_t), pointer, save :: acc_ef => null()
  type(variable_t), pointer, save :: acc_edgesol => null()
  type(variable_t), pointer, save :: acc_edgeres => null()
  type(valuelist_t), pointer, save :: acc_edgeparams => null()
  type(element_t), pointer, save :: acc_elem => null()
  ! The entry SAVE matrix pointer (NULL until the first PrecUse build).
  type(matrix_t), pointer, save :: acc_proj => null()

  ! Host-variable mirrors (entry locals in the original).
  integer, save :: acc_dim = 0, acc_dofs = 0, acc_compn = 0
  integer, save :: acc_relorder = 0
  logical, save :: acc_monolithic = .false., acc_segregated = .false.
  logical, save :: acc_precuse = .false., acc_curlcurl = .false.
  logical, save :: acc_hasprecdamp = .false.
  complex(kind=dp), save :: acc_precdamp = (0.0_dp, 0.0_dp)
  real(kind=dp), save :: acc_mu0inv = 0.0_dp, acc_eps0 = 0.0_dp
  real(kind=dp), save :: acc_omega = 0.0_dp, acc_rob0 = 0.0_dp
  real(kind=dp), save :: acc_norm(3) = 0.0_dp

  ! LocalMatrix SAVE state.
  type(nodes_t), target, save :: acc_nodes_lm
  type(valuehandle_t), save :: acc_lm_cond, acc_lm_eps, acc_lm_mu, acc_lm_curr
  real(kind=dp), allocatable, save :: acc_lm_basis(:), acc_lm_dbasis(:,:)
  complex(kind=dp), allocatable, save :: acc_lm_stiff(:,:,:), acc_lm_force(:,:)
  complex(kind=dp), save :: acc_kappa = (1.0_dp, 0.0_dp)
  type(gaussintegrationpoints_t), save :: acc_ip_lm

  ! LocalMatrixBC SAVE state.
  type(nodes_t), target, save :: acc_nodes_lbc
  logical, save :: acc_gotport
  type(valuehandle_t), save :: acc_lbc_robin, acc_lbc_magload, acc_lbc_absorb, &
      acc_lbc_goodcond, acc_lbc_temre, acc_lbc_temim, acc_lbc_mu, &
      acc_lbc_cond, acc_lbc_relnu
  complex(kind=dp), allocatable, save :: acc_lbc_stiff(:,:,:), acc_lbc_force(:,:)
  real(kind=dp), save :: acc_normal(3)
  logical, save :: acc_nt = .false.
  logical, save :: acc_goodcond = .false., acc_absorb = .false.
  logical, save :: acc_found_bc = .false.
  type(valuelist_t), pointer, save :: acc_bc => null()
  type(gaussintegrationpoints_t), save :: acc_ip_lbc
  type(element_t), pointer, save :: acc_parent => null()

contains

  ! -- VectorHelmholtzNodal_init entry ---------------------------------
  subroutine vhn_init_begin(h_model, h_solver) bind(c, name="vhn_init_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    acc_params => GetSolverParams()
  end subroutine vhn_init_begin

  function vhn_init_precuse() result(r) bind(c, name="vhn_init_precuse")
    integer(c_int) :: r
    logical :: found, v
    v = ListGetLogical(acc_params, 'Preconditioning Solver', found)
    r = merge(1, 0, v)
  end function vhn_init_precuse

  function vhn_init_monolithic() result(r) bind(c, name="vhn_init_monolithic")
    integer(c_int) :: r
    logical :: found, v
    v = ListGetLogical(acc_params, 'Monolithic Solver', found)
    r = merge(1, 0, v)
  end function vhn_init_monolithic

  function vhn_init_curlcurl() result(r) bind(c, name="vhn_init_curlcurl")
    integer(c_int) :: r
    logical :: found, v
    v = ListGetLogical(acc_params, 'curl-curl Form', found)
    r = merge(1, 0, v)
  end function vhn_init_curlcurl

  subroutine vhn_init_set_monolithic() bind(c, name="vhn_init_set_monolithic")
    call ListAddLogical(acc_params, 'Monolithic Solver', .TRUE.)
  end subroutine vhn_init_set_monolithic

  subroutine vhn_init_var_prec() bind(c, name="vhn_init_var_prec")
    call ListAddNewString(acc_params, 'Variable', &
        "Prec Elfield[PEX re:1 PEX im:1 PEY re:1 PEY im:1 PEZ re:1 PEZ im:1]")
  end subroutine vhn_init_var_prec

  subroutine vhn_init_var_field() bind(c, name="vhn_init_var_field")
    call ListAddNewString(acc_params, 'Variable', &
        "Elfield[EX re:1 EX im:1 EY re:1 EY im:1 EZ re:1 EZ im:1]")
  end subroutine vhn_init_var_field

  subroutine vhn_init_var_etmp() bind(c, name="vhn_init_var_etmp")
    call ListAddNewLogical(acc_params, 'Variable Output', .FALSE.)
    call ListAddNewString(acc_params, 'Variable', 'Etmp[Etmp re:1 Etmp im:1]')
  end subroutine vhn_init_var_etmp

  subroutine vhn_init_export_prec() bind(c, name="vhn_init_export_prec")
    call ListAddString(acc_params, &
        NextFreeKeyword('Exported Variable', acc_params), &
        "Prec Elfield[PEX re:1 PEX im:1 PEY re:1 PEY im:1 PEZ re:1 PEZ im:1]")
  end subroutine vhn_init_export_prec

  subroutine vhn_init_export_field() bind(c, name="vhn_init_export_field")
    call ListAddString(acc_params, &
        NextFreeKeyword('Exported Variable', acc_params), &
        "Elfield[EX re:1 EX im:1 EY re:1 EY im:1 EZ re:1 EZ im:1]")
  end subroutine vhn_init_export_field

  subroutine vhn_init_complex_time() bind(c, name="vhn_init_complex_time")
    call ListAddNewLogical(acc_params, "Linear System Complex", .TRUE.)
    call ListAddInteger(acc_params, 'Time Derivative Order', 0)
  end subroutine vhn_init_complex_time

  function vhn_init_findeigen() result(r) bind(c, name="vhn_init_findeigen")
    integer(c_int) :: r
    logical :: found, findeigen
    integer :: i
    findeigen = .FALSE.
    do i = 1, acc_model % NumberOfBCs
      if (ListGetString(acc_model % BCs(i) % Values, 'Port Type', found) == 'eigenmode') then
        findeigen = .TRUE.
        exit
      end if
    end do
    r = merge(1, 0, findeigen)
  end function vhn_init_findeigen

  subroutine vhn_init_eigensolver() bind(c, name="vhn_init_eigensolver")
    logical :: found
    character(len=MAX_NAME_LEN) :: sname
    integer :: soln, i, j
    soln = 0
    do i = 1, acc_model % NumberOfSolvers
      sname = GetString(acc_model % Solvers(i) % Values, 'Procedure', found)
      j = INDEX(sname, 'EMPortSolver')
      if (j > 0) then
        soln = i
        exit
      end if
    end do
    if (soln == 0) then
      call Fatal('VectorHelmholtzNodal_Init', &
          'Eigenfunction BC given without solving a port model')
    else
      call Info('VectorHelmholtzNodal_Init', &
          'The eigensolver index is: '//I2S(soln), Level=12)
      call ListAddInteger(acc_params, 'Eigensolver Index', soln)
    end if
  end subroutine vhn_init_eigensolver

  ! -- VectorHelmholtzNodal entry prologue ------------------------------
  subroutine vhn_begin(h_model, h_solver) bind(c, name="vhn_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    acc_norm = 0.0_dp
  end subroutine vhn_begin

  subroutine vhn_info_header() bind(c, name="vhn_info_header")
    call Info(caller, '', Level=8)
    call Info(caller, '------------------------------------------------', Level=6)
    call Info(caller, 'Solving harmonic electric waves using nodal basis!')
  end subroutine vhn_info_header

  function vhn_dim() result(n) bind(c, name="vhn_dim")
    integer(c_int) :: n
    acc_dim = CoordinateSystemDimension()
    if (CurrentCoordinateSystem() /= Cartesian) then
      call Fatal(caller, 'Only implemented for Cartesian problems!')
    end if
    n = acc_dim
  end function vhn_dim

  subroutine vhn_mesh_params() bind(c, name="vhn_mesh_params")
    acc_mesh => GetMesh()
    acc_params => GetSolverParams()
  end subroutine vhn_mesh_params

  subroutine vhn_dofs_branch(compn, segregated) bind(c, name="vhn_dofs_branch")
    integer(c_int) :: compn, segregated
    logical :: found
    acc_monolithic = ListGetLogical(acc_params, 'Monolithic Solver', found)
    acc_ei => acc_solver % Variable
    acc_dofs = acc_ei % Dofs / 2
    if (acc_dofs == 1) then
      if (acc_monolithic) call Fatal(caller, &
          'Variable DOFs incompatible with the segregated solution')
      call Info(caller, 'Treating the equation in segregated manner!')
      acc_compn = acc_dim
    else if (acc_dofs == acc_dim) then
      if (.NOT. acc_monolithic) call Fatal(caller, &
          'Variable DOFs incompatible with the monolithic solution')
      call Info(caller, 'Treating the equation in monolithic manner!')
      acc_compn = 1
    else
      call Fatal(caller, &
          'Invalid number of dofs for the solver variable: '//I2S(acc_dofs))
    end if
    acc_segregated = .NOT. acc_monolithic
    compn = acc_compn
    segregated = merge(1, 0, acc_segregated)
  end subroutine vhn_dofs_branch

  function vhn_field_setup() result(precuse) bind(c, name="vhn_field_setup")
    integer(c_int) :: precuse
    logical :: found, piolaversion, secondfamily, secondorder
    character(len=MAX_NAME_LEN) :: sname
    acc_precuse = ListGetLogical(acc_params, 'Preconditioning Solver', found)
    acc_curlcurl = ListGetLogical(acc_params, 'curl-curl Form', found)
    if (acc_precuse) then
      if (.NOT. acc_monolithic) call Fatal(caller, &
          'The use as a preconditioner needs Monolithic Solver = True')
      acc_ef => VariableGet(acc_mesh % Variables, 'Prec ElField')
      acc_edgesol => NULL()
      sname = ListGetString(acc_params, 'Edge Update Name', found)
      if (found) then
        acc_edgesol => VariableGet(acc_mesh % Variables, sname)
        if (.NOT. ASSOCIATED(acc_edgesol)) &
            call Fatal(caller, 'Could not found field: '//TRIM(sname))
      else
        call Warn(caller, &
            'Give Edge Update Name to enable the use as a preconditioner')
        acc_precuse = .FALSE.
      end if
      acc_edgeres => NULL()
      sname = ListGetString(acc_params, 'Edge Residual Name', found)
      if (found) then
        acc_edgeres => VariableGet(acc_mesh % Variables, sname)
        if (.NOT. ASSOCIATED(acc_edgeres)) &
            call Fatal(caller, 'Could not find field: '//TRIM(sname))
        acc_edgeparams => GetSolverParams(acc_edgeres % Solver)
        call EdgeElementStyle(acc_edgeparams, piolaversion, secondfamily, &
            secondorder, Check = .TRUE.)
        if (secondorder) call Fatal(caller, 'The lowest-order edge basis assumed')
      else
        call Warn(caller, &
            'Give Edge Residual Name to enable the use as a preconditioner')
        acc_precuse = .FALSE.
      end if
      if (acc_precuse) then
        acc_precdamp = GetCReal(acc_params, &
            'Linear System Preconditioning Damp Coefficient', acc_hasprecdamp)
        acc_precdamp = CMPLX(REAL(acc_precdamp), GetCReal(acc_params, &
            'Linear System Preconditioning Damp Coefficient im', found), KIND=dp)
        acc_hasprecdamp = acc_hasprecdamp .OR. found
      end if
    else
      acc_ef => VariableGet(acc_mesh % Variables, 'ElField')
    end if
    if (.NOT. ASSOCIATED(acc_ef)) then
      call Fatal(caller, 'Variable for Electric field not found!')
    end if
    precuse = merge(1, 0, acc_precuse)
  end function vhn_field_setup

  function vhn_follow() result(r) bind(c, name="vhn_follow")
    integer(c_int) :: r
    logical :: found, v
    v = ListGetLogical(acc_params, 'Follow P Curvature', found)
    r = merge(1, 0, v)
  end function vhn_follow

  subroutine vhn_follow_apply() bind(c, name="vhn_follow_apply")
    call FollowCurvedBoundary(acc_model, acc_mesh, .TRUE.)
  end subroutine vhn_follow_apply

  subroutine vhn_default_start() bind(c, name="vhn_default_start")
    call DefaultStart()
  end subroutine vhn_default_start

  function vhn_maxiter() result(n) bind(c, name="vhn_maxiter")
    integer(c_int) :: n
    logical :: found
    n = ListGetInteger(acc_params, 'Nonlinear System Max Iterations', found, minv=1)
    if (.NOT. found) n = 1
  end function vhn_maxiter

  function vhn_relorder() result(n) bind(c, name="vhn_relorder")
    integer(c_int) :: n
    logical :: found
    acc_relorder = GetInteger(acc_params, 'Relative Integration Order', found)
    n = acc_relorder
  end function vhn_relorder

  subroutine vhn_init_stuff() bind(c, name="vhn_init_stuff")
    logical :: found
    found = .FALSE.
    if (ASSOCIATED(acc_model % Constants)) then
      acc_mu0inv = 1.0_dp / GetConstReal(acc_model % Constants, &
          'Permeability of Vacuum', found)
    end if
    if (.NOT. found) acc_mu0inv = 1.0_dp / (PI * 4.0d-7)
    found = .FALSE.
    if (ASSOCIATED(acc_model % Constants)) then
      acc_eps0 = GetConstReal(acc_model % Constants, &
          'Permittivity of Vacuum', found)
    end if
    if (.NOT. found) acc_eps0 = 8.854187817d-12
    acc_omega = GetAngularFrequency(Found=found)
    if (.NOT. found) call Fatal(caller, 'We should have Omega!')
    acc_rob0 = acc_omega * SQRT(acc_eps0 / acc_mu0inv)
  end subroutine vhn_init_stuff

  subroutine vhn_proj() bind(c, name="vhn_proj")
    if (acc_precuse .AND. .NOT. ASSOCIATED(acc_proj)) then
      call Info(caller, 'Creating projection matrix to map a nodal '// &
          'solution into vector element space', Level=6)
      call NodalToNedelecInterpolation_GlobalMatrix(acc_mesh, acc_ef, &
          acc_edgesol, acc_proj, cdim=3)
    end if
  end subroutine vhn_proj

  ! -- compi-loop phases ------------------------------------------------
  subroutine vhn_info_compi(compi) bind(c, name="vhn_info_compi")
    integer(c_int), value :: compi
    call Info(caller, 'Solving for component '//I2S(compi), Level=6)
  end subroutine vhn_info_compi

  subroutine vhn_default_initialize() bind(c, name="vhn_default_initialize")
    call DefaultInitialize()
  end subroutine vhn_default_initialize

  function vhn_nof_active() result(n) bind(c, name="vhn_nof_active")
    integer(c_int) :: n
    n = GetNOFActive(acc_solver)
  end function vhn_nof_active

  function vhn_nof_boundary() result(n) bind(c, name="vhn_nof_boundary")
    integer(c_int) :: n
    n = GetNOFBoundaryActive(acc_solver)
  end function vhn_nof_boundary

  subroutine vhn_begin_active(t, n, nd, nb) bind(c, name="vhn_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, nb
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes(acc_elem)
    nd = GetElementNOFDOFs(acc_elem)
    nb = GetElementNOFBDOFs(acc_elem)
  end subroutine vhn_begin_active

  subroutine vhn_begin_boundary(t, n, nd, active) bind(c, name="vhn_begin_boundary")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, active
    acc_elem => GetBoundaryElement(t)
    active = merge(1, 0, ActiveBoundaryElement(acc_elem))
    n = GetElementNOFNodes(acc_elem)
    nd = GetElementNOFDOFs(acc_elem)
  end subroutine vhn_begin_boundary

  ! LocalMatrix head: InitHandles block, Gauss points, first-mesh SAVE
  ! allocation, GetElementNodes, kappa, and the zeroing.
  subroutine vhn_lm_setup(h) bind(c, name="vhn_lm_setup")
    integer(c_int) :: h
    integer :: m, k, q, allocstat
    if (h /= 0) then
      call ListInitElementKeyword(acc_lm_cond, 'Material', 'Electric Conductivity')
      call ListInitElementKeyword(acc_lm_eps, 'Material', &
          'Relative Permittivity', InitIm=.TRUE.)
      call ListInitElementKeyword(acc_lm_mu, 'Material', &
          'Relative Reluctivity', InitIm=.TRUE.)
      call ListInitElementKeyword(acc_lm_curr, 'Body Force', &
          'Current Density', InitIm=.TRUE., InitVec3D=.TRUE.)
      h = 0
    end if
    if (acc_relorder /= 0) then
      acc_ip_lm = GaussPoints(acc_elem, RelOrder=acc_relorder)
    else
      acc_ip_lm = GaussPoints(acc_elem)
    end if
    if (.NOT. ALLOCATED(acc_lm_basis)) then
      m = acc_mesh % MaxElementDofs
      if (acc_curlcurl) then
        k = 1
        q = 3*m
      else
        k = 3
        q = m
      end if
      allocate(acc_lm_basis(m), acc_lm_dbasis(m, 3), &
          acc_lm_stiff(q, q, k), acc_lm_force(q, k), STAT=allocstat)
      if (allocstat /= 0) call Fatal(caller, 'Local storage allocation failed')
    end if
    call GetElementNodes(acc_nodes_lm, UElement=acc_elem)
    acc_kappa = CMPLX(1.0_dp, 0.0_dp, kind=dp)
    if (acc_precuse) then
      if (acc_hasprecdamp) acc_kappa = acc_kappa - acc_precdamp
    end if
    acc_lm_stiff = 0._dp
    acc_lm_force = 0._dp
  end subroutine vhn_lm_setup

  function vhn_lm_nip() result(n) bind(c, name="vhn_lm_nip")
    integer(c_int) :: n
    n = acc_ip_lm % n
  end function vhn_lm_nip

  ! LocalMatrix per-IP body (ElementInfo, coefficient fetches, then the
  ! curl-curl or component-wise STIFF/FORCE accumulation with the
  ! original MATMUL reduction). `nd` is the combined nd+nb the entry
  ! passes to the internal procedure.
  subroutine vhn_lm_ip(t, compi, nd) bind(c, name="vhn_lm_ip")
    integer(c_int), value :: t, compi, nd
    real(kind=dp) :: weight, detj, condatip
    complex(kind=dp) :: muinvatip, epsatip
    logical :: stat, found, withconductivity
    stat = ElementInfo(acc_elem, acc_nodes_lm, acc_ip_lm % U(t), &
        acc_ip_lm % V(t), acc_ip_lm % W(t), detj, acc_lm_basis, acc_lm_dbasis)
    weight = acc_ip_lm % s(t) * detj
    muinvatip = ListGetElementComplex(acc_lm_mu, acc_lm_basis, acc_elem, &
        found, GaussPoint=t)
    if (found) then
      muinvatip = muinvatip * acc_mu0inv
    else
      muinvatip = acc_mu0inv
    end if
    epsatip = ListGetElementComplex(acc_lm_eps, acc_lm_basis, acc_elem, &
        found, GaussPoint=t)
    if (found) then
      epsatip = epsatip * acc_eps0
    else
      epsatip = acc_eps0
    end if
    condatip = ListGetElementReal(acc_lm_cond, acc_lm_basis, acc_elem, &
        withconductivity, GaussPoint=t)
    if (acc_curlcurl) then
      call vhn_lm_ip_curl(nd, weight, muinvatip, epsatip, &
          condatip, withconductivity)
    else
      call vhn_lm_ip_scalar(nd, compi, weight, muinvatip, epsatip, &
          condatip, withconductivity)
    end if
  end subroutine vhn_lm_ip

  ! Curl-curl STIFF accumulation (verbatim p/j/q nest over 1:nd with
  ! dim-strided destination rows).
  subroutine vhn_lm_ip_curl(nd, weight, muinvatip, epsatip, condatip, withcond)
    integer :: nd
    real(kind=dp) :: weight, condatip
    complex(kind=dp) :: muinvatip, epsatip
    logical :: withcond
    integer :: j, p, q
    do p = 1, nd
      do j = 1, acc_dim
        do q = 1, nd
          acc_lm_stiff(acc_dim*(p-1)+j, acc_dim*(q-1)+j, 1) = &
              acc_lm_stiff(acc_dim*(p-1)+j, acc_dim*(q-1)+j, 1) - &
              weight * acc_kappa * acc_omega**2 * epsatip * &
              acc_lm_basis(q) * acc_lm_basis(p)
          if (withcond) then
            acc_lm_stiff(acc_dim*(p-1)+j, acc_dim*(q-1)+j, 1) = &
                acc_lm_stiff(acc_dim*(p-1)+j, acc_dim*(q-1)+j, 1) - &
                im * weight * acc_omega * condatip * &
                acc_lm_basis(q) * acc_lm_basis(p)
          end if
          if (j == 1) then
            acc_lm_stiff(acc_dim*(p-1)+j, acc_dim*(q-1)+1, 1) = &
                acc_lm_stiff(acc_dim*(p-1)+j, acc_dim*(q-1)+1, 1) + &
                weight * muinvatip * ( &
                acc_lm_dbasis(q, 3) * acc_lm_dbasis(p, 3) + &
                acc_lm_dbasis(q, 2) * acc_lm_dbasis(p, 2))
            acc_lm_stiff(acc_dim*(p-1)+j, acc_dim*(q-1)+2, 1) = &
                acc_lm_stiff(acc_dim*(p-1)+j, acc_dim*(q-1)+2, 1) - &
                weight * muinvatip * &
                acc_lm_dbasis(q, 1) * acc_lm_dbasis(p, 2)
            acc_lm_stiff(acc_dim*(p-1)+j, acc_dim*(q-1)+3, 1) = &
                acc_lm_stiff(acc_dim*(p-1)+j, acc_dim*(q-1)+3, 1) - &
                weight * muinvatip * &
                acc_lm_dbasis(q, 1) * acc_lm_dbasis(p, 3)
          end if
          if (j == 2) then
            acc_lm_stiff(acc_dim*(p-1)+j, acc_dim*(q-1)+1, 1) = &
                acc_lm_stiff(acc_dim*(p-1)+j, acc_dim*(q-1)+1, 1) - &
                weight * muinvatip * &
                acc_lm_dbasis(q, 2) * acc_lm_dbasis(p, 1)
            acc_lm_stiff(acc_dim*(p-1)+j, acc_dim*(q-1)+2, 1) = &
                acc_lm_stiff(acc_dim*(p-1)+j, acc_dim*(q-1)+2, 1) + &
                weight * muinvatip * ( &
                acc_lm_dbasis(q, 3) * acc_lm_dbasis(p, 3) + &
                acc_lm_dbasis(q, 1) * acc_lm_dbasis(p, 1))
            acc_lm_stiff(acc_dim*(p-1)+j, acc_dim*(q-1)+3, 1) = &
                acc_lm_stiff(acc_dim*(p-1)+j, acc_dim*(q-1)+3, 1) - &
                weight * muinvatip * &
                acc_lm_dbasis(q, 2) * acc_lm_dbasis(p, 3)
          end if
          if (j == 3) then
            acc_lm_stiff(acc_dim*(p-1)+j, acc_dim*(q-1)+1, 1) = &
                acc_lm_stiff(acc_dim*(p-1)+j, acc_dim*(q-1)+1, 1) - &
                weight * muinvatip * &
                acc_lm_dbasis(q, 3) * acc_lm_dbasis(p, 1)
            acc_lm_stiff(acc_dim*(p-1)+j, acc_dim*(q-1)+2, 1) = &
                acc_lm_stiff(acc_dim*(p-1)+j, acc_dim*(q-1)+2, 1) - &
                weight * muinvatip * &
                acc_lm_dbasis(q, 3) * acc_lm_dbasis(p, 2)
            acc_lm_stiff(acc_dim*(p-1)+j, acc_dim*(q-1)+3, 1) = &
                acc_lm_stiff(acc_dim*(p-1)+j, acc_dim*(q-1)+3, 1) + &
                weight * muinvatip * ( &
                acc_lm_dbasis(q, 2) * acc_lm_dbasis(p, 2) + &
                acc_lm_dbasis(q, 1) * acc_lm_dbasis(p, 1))
          end if
        end do
      end do
    end do
  end subroutine vhn_lm_ip_curl

  ! Component-wise STIFF/FORCE accumulation (verbatim scalar nest with
  ! the original MATMUL diffusion reduction).
  subroutine vhn_lm_ip_scalar(nd, compi, weight, muinvatip, epsatip, &
      condatip, withcond)
    integer :: nd, compi
    real(kind=dp) :: weight, condatip
    complex(kind=dp) :: muinvatip, epsatip
    logical :: withcond
    logical :: found
    complex(kind=dp) :: curratip(3)
    integer :: i, p
    acc_lm_stiff(1:nd, 1:nd, 1) = acc_lm_stiff(1:nd, 1:nd, 1) + weight * &
        muinvatip * MATMUL(acc_lm_dbasis(1:nd, :), &
        TRANSPOSE(acc_lm_dbasis(1:nd, :)))
    if (withcond) then
      do p = 1, nd
        acc_lm_stiff(p, 1:nd, 1) = acc_lm_stiff(p, 1:nd, 1) - &
            im * weight * acc_omega * condatip * acc_lm_basis(p) * &
            acc_lm_basis(1:nd)
      end do
    end if
    do p = 1, nd
      acc_lm_stiff(p, 1:nd, 1) = acc_lm_stiff(p, 1:nd, 1) - &
          weight * acc_kappa * acc_omega**2 * epsatip * acc_lm_basis(p) * &
          acc_lm_basis(1:nd)
    end do
    if (.NOT. acc_precuse) then
      curratip = ListGetElementComplex3D(acc_lm_curr, acc_lm_basis, &
          acc_elem, found)
      if (found) then
        if (acc_monolithic) then
          do i = 1, acc_dofs
            acc_lm_force(1:nd, i) = acc_lm_force(1:nd, i) + &
                weight * curratip(i) * acc_lm_basis(1:nd)
          end do
        else
          acc_lm_force(1:nd, 1) = acc_lm_force(1:nd, 1) + &
              weight * curratip(compi) * acc_lm_basis(1:nd)
        end if
      end if
    end if
  end subroutine vhn_lm_ip_scalar

  ! LocalMatrix tail: DefaultUpdateEquations routing (verbatim).
  subroutine vhn_lm_update(compi, nd) bind(c, name="vhn_lm_update")
    integer(c_int), value :: compi, nd
    integer :: i
    if (acc_curlcurl) then
      call DefaultUpdateEquations( &
          acc_lm_stiff(1:acc_dim*nd, 1:acc_dim*nd, 1), &
          acc_lm_force(1:acc_dim*nd, 1), UElement=acc_elem)
    else
      if (acc_monolithic) then
        do i = 2, acc_dofs
          acc_lm_stiff(1:nd, 1:nd, i) = acc_lm_stiff(1:nd, 1:nd, 1)
        end do
        call DefaultUpdateEquations(acc_lm_stiff, acc_lm_force, UElement=acc_elem)
      else
        call DefaultUpdateEquations(acc_lm_stiff(:, :, 1), &
            acc_lm_force(:, 1), UElement=acc_elem)
      end if
    end if
  end subroutine vhn_lm_update

  subroutine vhn_finish_bulk() bind(c, name="vhn_finish_bulk")
    call DefaultFinishBulkAssembly()
  end subroutine vhn_finish_bulk

  subroutine vhn_apply_residual() bind(c, name="vhn_apply_residual")
    if (acc_precuse) then
      if (acc_monolithic) then
        call Info(caller, &
            'Using Transposed Projection Matrix: H(curl) -> H1', Level=6)
        call CRS_TransposeMatrixVectorMultiply(acc_proj, &
            acc_edgeres % Values, acc_solver % Matrix % rhs)
      else
        ! TO DO: Add the transformation of the residual for the
        ! component-wise wave equation
      end if
    end if
  end subroutine vhn_apply_residual

  subroutine vhn_info_range() bind(c, name="vhn_info_range")
    if (InfoActive(20)) then
      call VectorValuesRange(acc_solver % Matrix % Values, &
          SIZE(acc_solver % Matrix % Values), 'A0')
    end if
  end subroutine vhn_info_range

  ! LocalMatrixBC head: BC fetch, InitHandles block (+ElectricPortModel
  ! phase 1), GetElementNodes, first-mesh SAVE allocation, zeroing,
  ! Normal/NT/GoodConductor/Absorb fetches, ElectricPortModel phase 2,
  ! Gauss points and the bulk parent.
  subroutine vhn_lbc_head(h, proceed, nip) bind(c, name="vhn_lbc_head")
    integer(c_int) :: h, proceed, nip
    integer :: m, k, q, allocstat
    logical :: found
    acc_bc => GetBC(acc_elem)
    if (.NOT. ASSOCIATED(acc_bc)) then
      proceed = 0
      nip = 0
      return
    end if
    proceed = 1
    if (h /= 0) then
      call ListInitElementKeyword(acc_lbc_robin, 'Boundary Condition', &
          'Electric Robin Coefficient', InitIm=.TRUE.)
      call ListInitElementKeyword(acc_lbc_magload, 'Boundary Condition', &
          'Magnetic Boundary Load', InitIm=.TRUE., InitVec3D=.TRUE.)
      call ListInitElementKeyword(acc_lbc_absorb, 'Boundary Condition', &
          'Absorbing BC')
      call ListInitElementKeyword(acc_lbc_goodcond, 'Boundary Condition', &
          'Good Conductor BC')
      call ListInitElementKeyword(acc_lbc_temre, 'Boundary Condition', &
          'TEM Potential')
      call ListInitElementKeyword(acc_lbc_temim, 'Boundary Condition', &
          'TEM Potential Im')
      call ListInitElementKeyword(acc_lbc_mu, 'Material', &
          'Relative Reluctivity', InitIm=.TRUE.)
      call ListInitElementKeyword(acc_lbc_cond, 'Boundary Condition', &
          'Layer Electric Conductivity', InitIm=.TRUE.)
      call ListInitElementKeyword(acc_lbc_relnu, 'Boundary Condition', &
          'Layer Relative Reluctivity', InitIm=.TRUE.)
      call ElectricPortModel(1, acc_solver)
      h = 0
    end if
    call GetElementNodes(acc_nodes_lbc, UElement=acc_elem)
    if (.NOT. ALLOCATED(acc_lbc_stiff)) then
      m = acc_mesh % MaxElementDofs
      if (acc_curlcurl) then
        k = 1
        q = 3*m
      else
        k = 3
        q = m
      end if
      allocate(acc_lbc_stiff(q, q, k), acc_lbc_force(q, k), STAT=allocstat)
      if (allocstat /= 0) call Fatal(caller, 'Local storage allocation failed')
    end if
    acc_lbc_stiff = 0._dp
    acc_lbc_force = 0._dp
    if (.NOT. acc_curlcurl) acc_normal = NormalVector(acc_elem, acc_nodes_lbc)
    acc_nt = ListGetLogical(acc_bc, &
        'Normal-Tangential '//GetVarName(acc_solver % Variable), found)
    if (acc_nt .AND. .NOT. acc_monolithic) then
      call Fatal(caller, &
          'Normal-tangential conditions require monolithic solver!')
    end if
    acc_goodcond = ListGetElementLogical(acc_lbc_goodcond, acc_elem, found)
    acc_absorb = ListGetElementLogical(acc_lbc_absorb, acc_elem, found)
    acc_found_bc = found
    call ElectricPortModel(2, acc_solver, acc_elem, acc_gotport)
    acc_ip_lbc = GaussPoints(acc_elem)
    acc_parent => GetBulkElementAtBoundary(acc_elem)
    nip = acc_ip_lbc % n
  end subroutine vhn_lbc_head

  ! LocalMatrixBC per-IP body (ElementInfo, coefficient fetches, FORCE
  ! update, curl-curl or component-wise STIFF update), verbatim. The
  ! original's procedure-level Found is the module acc_found_bc, primed
  ! by the head wrapper with the Absorb-fetch Found.
  subroutine vhn_lbc_ip(t, nd) bind(c, name="vhn_lbc_ip")
    integer(c_int), value :: t, nd
    real(kind=dp) :: basis(nd), dbasisdx(nd, 3)
    real(kind=dp) :: F, C, Ext, weight, detj, coeff
    complex(kind=dp) :: muinvatip, muinv, cond, surfimp, temgrad(3), l(3), b
    real(kind=dp) :: testvec(3), trialvec(3)
    logical :: stat, gotsome
    integer :: i, p, q, j
    stat = ElementInfo(acc_elem, acc_nodes_lbc, acc_ip_lbc % U(t), &
        acc_ip_lbc % V(t), acc_ip_lbc % W(t), detj, basis, dbasisdx)
    weight = acc_ip_lbc % s(t) * detj
    if (ASSOCIATED(acc_parent)) then
      muinvatip = ListGetElementComplex(acc_lbc_mu, basis, acc_parent, &
          acc_found_bc, GaussPoint=t)
    end if
    if (acc_found_bc) then
      muinvatip = muinvatip * acc_mu0inv
    else
      muinvatip = acc_mu0inv
    end if
    l = 0.0_dp
    gotsome = .TRUE.
    if (ListGetElementLogical(acc_lbc_absorb, acc_elem, acc_found_bc)) then
      b = CMPLX(0.0_dp, acc_rob0, KIND=dp)
    else if (acc_goodcond) then
      cond = ListGetElementComplex(acc_lbc_cond, basis, acc_elem, &
          acc_found_bc, GaussPoint=t)
      muinv = ListGetElementComplex(acc_lbc_relnu, basis, acc_elem, &
          acc_found_bc, GaussPoint=t)
      if (acc_found_bc) then
        muinv = muinv * acc_mu0inv
      else
        muinv = acc_mu0inv
      end if
      surfimp = CMPLX(1.0_dp, -1.0_dp, KIND=dp) * &
          SQRT(acc_omega / (2.0_dp * cond * muinv))
      b = 1.0_dp / surfimp
    else if (acc_gotport) then
      if (acc_precuse) then
        call ElectricPortModel(3, acc_solver, acc_elem, acc_gotport, b, &
            Basis=basis, dBasisdx=dbasisdx)
      else
        call ElectricPortModel(3, acc_solver, acc_elem, acc_gotport, b, l, &
            Basis=basis, dBasisdx=dbasisdx)
      end if
    else
      if (.NOT. acc_precuse) then
        l = ListGetElementComplex3D(acc_lbc_magload, basis, acc_elem, &
            acc_found_bc, GaussPoint=t)
        temgrad = CMPLX(ListGetElementRealGrad(acc_lbc_temre, dbasisdx, &
            acc_elem, acc_found_bc), ListGetElementRealGrad(acc_lbc_temim, &
            dbasisdx, acc_elem, acc_found_bc), KIND=dp)
        l = l + temgrad
      end if
      b = ListGetElementComplex(acc_lbc_robin, basis, acc_elem, &
          acc_found_bc, GaussPoint=t)
      gotsome = acc_found_bc
    end if
    if (gotsome) then
      if (.NOT. acc_precuse) then
        do i = 1, acc_dim
          acc_lbc_force(1:nd, i) = acc_lbc_force(1:nd, i) - &
              muinvatip * l(i) * basis(1:nd) * weight
        end do
      end if
      if (acc_curlcurl) then
        acc_normal = NormalVector(acc_elem, acc_nodes_lbc, &
            acc_ip_lbc % U(t), acc_ip_lbc % V(t), .TRUE.)
        do p = 1, nd
          do j = 1, acc_dim
            testvec = 0.0d0
            testvec(j) = basis(p)
            testvec = CrossProduct(testvec, acc_normal)
            do q = 1, nd
              do i = 1, acc_dim
                trialvec = 0.0d0
                trialvec(i) = basis(q)
                trialvec = CrossProduct(trialvec, acc_normal)
                acc_lbc_stiff(acc_dim*(p-1)+j, acc_dim*(q-1)+i, 1) = &
                    acc_lbc_stiff(acc_dim*(p-1)+j, acc_dim*(q-1)+i, 1) - &
                    muinvatip * b * SUM(testvec(:) * trialvec(:)) * weight
              end do
            end do
          end do
        end do
      else
        do i = 1, acc_dim
          if (acc_nt) then
            if (i == 1) cycle
            coeff = 1.0_dp
          else
            coeff = 1.0_dp - acc_normal(i)**2
            if (coeff <= 0.0_dp) then
              coeff = 0.0_dp
            else
              coeff = SQRT(coeff)
            end if
          end if
          do p = 1, nd
            acc_lbc_stiff(p, 1:nd, i) = acc_lbc_stiff(p, 1:nd, i) - &
                coeff * muinvatip * b * basis(p) * basis(1:nd) * weight
          end do
        end do
      end if
    end if
  end subroutine vhn_lbc_ip

  ! LocalMatrixBC tail: DefaultUpdateEquations routing (verbatim).
  subroutine vhn_lbc_update(compi, nd) bind(c, name="vhn_lbc_update")
    integer(c_int), value :: compi, nd
    if (acc_curlcurl) then
      call DefaultUpdateEquations( &
          acc_lbc_stiff(1:acc_dim*nd, 1:acc_dim*nd, 1), &
          acc_lbc_force(1:acc_dim*nd, 1), UElement=acc_elem)
    else
      if (acc_monolithic) then
        call DefaultUpdateEquations(acc_lbc_stiff, acc_lbc_force, &
            UElement=acc_elem)
      else
        call DefaultUpdateEquations(acc_lbc_stiff(:, :, compi), &
            acc_lbc_force(:, compi), UElement=acc_elem)
      end if
    end if
  end subroutine vhn_lbc_update

  subroutine vhn_finish_bnd() bind(c, name="vhn_finish_bnd")
    call DefaultFinishBoundaryAssembly()
  end subroutine vhn_finish_bnd

  subroutine vhn_finish_assem() bind(c, name="vhn_finish_assem")
    call DefaultFinishAssembly()
  end subroutine vhn_finish_assem

  ! Dirichlet routing: monolithic DefaultDirichletBCs, else the
  ! per-component SetDirichletBoundaries pair (verbatim).
  subroutine vhn_dirichlet(compi) bind(c, name="vhn_dirichlet")
    integer(c_int), value :: compi
    character(len=MAX_NAME_LEN) :: sname
    integer :: i
    if (acc_monolithic) then
      call DefaultDirichletBCs()
    else
      do i = 1, 2
        sname = ComponentName(acc_ef, 2*(compi-1)+i)
        call SetDirichletBoundaries(CurrentModel, acc_solver % Matrix, &
            acc_solver % Matrix % rhs, sname, i, 2, acc_solver % Variable % Perm)
      end do
      call EnforceDirichletConditions(acc_solver, acc_solver % Matrix, &
          acc_solver % Matrix % rhs)
    end if
  end subroutine vhn_dirichlet

  ! Segregated component gather (verbatim strided copies).
  subroutine vhn_gather(compi) bind(c, name="vhn_gather")
    integer(c_int), value :: compi
    if (acc_segregated) then
      acc_ei % Values(1::2) = acc_ef % Values(2*compi-1::2*acc_dim)
      acc_ei % Values(2::2) = acc_ef % Values(2*compi::2*acc_dim)
    end if
  end subroutine vhn_gather

  subroutine vhn_solve(compi) bind(c, name="vhn_solve")
    integer(c_int), value :: compi
    acc_norm(compi) = DefaultSolve()
  end subroutine vhn_solve

  ! Segregated component scatter (verbatim strided copies + the
  ! InfoActive(25) diagnostics with their PRINT).
  subroutine vhn_scatter(compi) bind(c, name="vhn_scatter")
    integer(c_int), value :: compi
    if (acc_segregated) then
      if (InfoActive(25)) then
        call VectorValuesRange(acc_ei % Values, SIZE(acc_ei % Values), &
            'E'//I2S(compi))
        print *, 'Component Norm:', acc_norm(compi)
      end if
      acc_ef % Values(2*compi-1::2*acc_dim) = acc_ei % Values(1::2)
      acc_ef % Values(2*compi::2*acc_dim) = acc_ei % Values(2::2)
    end if
  end subroutine vhn_scatter

  subroutine vhn_project() bind(c, name="vhn_project")
    if (acc_precuse) then
      call Info(caller, &
          'Projecting nodal solution to vector element space', Level=6)
      call Info(caller, 'Using Projection Matrix: H1 -> H(curl)', Level=6)
      call CRS_MatrixVectorMultiply(acc_proj, acc_ef % Values, &
          acc_edgesol % Values)
    end if
  end subroutine vhn_project

  subroutine vhn_finish() bind(c, name="vhn_finish")
    call DefaultFinish()
  end subroutine vhn_finish

  subroutine vhn_final_norm() bind(c, name="vhn_final_norm")
    if (acc_segregated) then
      acc_solver % Variable % Norm = &
          SQRT(SUM(acc_norm(1:acc_compn)**2) / acc_compn)
    end if
  end subroutine vhn_final_norm

  subroutine vhn_info_done() bind(c, name="vhn_info_done")
    call Info(caller, 'All done', Level=12)
  end subroutine vhn_info_done

end module vhn_acc
