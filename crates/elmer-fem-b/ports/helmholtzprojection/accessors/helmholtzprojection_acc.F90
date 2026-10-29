! T022-A HelmholtzProjection port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/HelmholtzProjection.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments. All Fortran-object
! traffic — TYPE(...) handles, OPTIONAL/POINTER dummies, hidden lengths,
! assumed-shape descriptors, whole-array updates — stays inside Fortran,
! so the Rust side is pure C calls. Each entry's SAVE state (first-mesh
! sizing, Nodes/IP context, potential-solver handles, Piola/SecondOrder
! flags, ConstantBulkMatrix flag) lives wrapper-side with identical
! lifetimes. All floating-point work (the edge-element Galerkin assembly,
! MATMUL/SUM reductions, the A-vector accumulation) stays in Fortran
! (same compiler), so bit-identity holds by construction wherever the
! call sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module helmholtzprojection_acc
  use iso_c_binding
  use defutils
  implicit none

  ! -- HelmholtzProjector entry state (mirrors its SAVE + locals) ---------
  type(model_t), pointer, save :: hp_model => null()
  type(solver_t), pointer, save :: hp_solver => null()
  type(solver_t), pointer, save :: hp_pot => null()
  type(valuelist_t), pointer, save :: hp_params => null()
  type(element_t), pointer, save :: hp_elem => null()
  type(nodes_t), target, save :: hp_nodes
  type(gaussintegrationpoints_t), save :: hp_ip
  character(len=MAX_NAME_LEN), save :: hp_potname = ''
  logical, save :: hp_piola = .false.
  logical, save :: hp_second = .false.
  integer, save :: hp_dimv = 0
  logical, save :: hp_cbflag = .false.

  ! -- RemoveKernelComponent entry state (mirrors its SAVE + locals) ------
  type(model_t), pointer, save :: rk_model => null()
  type(solver_t), pointer, save :: rk_solver => null()
  type(solver_t), pointer, save :: rk_pot => null()
  type(solver_t), pointer, save :: rk_ker => null()
  type(valuelist_t), pointer, save :: rk_params => null()
  type(element_t), pointer, save :: rk_elem => null()
  type(nodes_t), target, save :: rk_nodes
  type(gaussintegrationpoints_t), save :: rk_ip
  character(len=MAX_NAME_LEN), save :: rk_potname = ''
  character(len=MAX_NAME_LEN), save :: rk_kername = ''
  logical, save :: rk_piola = .false.
  logical, save :: rk_second = .false.
  integer, save :: rk_dimv = 0
  logical, save :: rk_cbflag = .false.

contains

  ! -- HelmholtzProjector entry ----------------------------------------
  subroutine hp_begin(h_model, h_solver) bind(c, name="hp_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, hp_model)
    call c_f_pointer(h_solver, hp_solver)
  end subroutine hp_begin

  subroutine hp_default_start() bind(c, name="hp_default_start")
    call DefaultStart()
  end subroutine hp_default_start

  function hp_dim() result(d) bind(c, name="hp_dim")
    integer(c_int) :: d
    hp_dimv = CoordinateSystemDimension()
    d = hp_dimv
  end function hp_dim

  function hp_params_begin() result(n) bind(c, name="hp_params_begin")
    integer(c_int) :: n
    type(mesh_t), pointer :: m
    hp_params => GetSolverParams()
    m => GetMesh()
    n = m % maxelementdofs
  end function hp_params_begin

  function hp_constbulk() result(v) bind(c, name="hp_constbulk")
    integer(c_int) :: v
    logical :: found
    hp_cbflag = GetLogical(hp_params, 'Constant Bulk Matrix', found)
    v = merge(1, 0, hp_cbflag)
  end function hp_constbulk

  function hp_ready() result(v) bind(c, name="hp_ready")
    integer(c_int) :: v
    v = merge(1, 0, associated(hp_solver % matrix % bulkvalues) .and. hp_cbflag)
  end function hp_ready

  ! Potential-variable solver lookup + FE-style inheritance (Fatal paths
  ! inside, exactly as in the original).
  subroutine hp_lookup() bind(c, name="hp_lookup")
    integer :: i, potdofs
    logical :: found
    hp_potname = GetString(hp_params, 'Potential Variable', found)
    if (.not. found) hp_potname = 'av'
    found = .false.
    do i = 1, hp_model % numberofsolvers
      hp_pot => hp_model % solvers(i)
      if (hp_potname == GetVarName(hp_pot % variable)) then
        found = .true.
        exit
      end if
    end do
    if (found) then
      call Info('HelmholtzProjection', 'Solver inherits potential ' &
          //trim(hp_potname)//' from solver: '//I2S(i), level=7)
    else
      call Fatal('HelmholtzProjection', 'Solver associated with potential variable > ' &
          //trim(hp_potname)//' < not found!')
    end if
    potdofs = hp_pot % variable % dofs
    if (potdofs > 1) call Fatal('HelmholtzProjection', 'A real-valued potential expected')
    call EdgeElementStyle(hp_pot % values, hp_piola, &
        QuadraticApproximation=hp_second)
    if (hp_piola) call Info('HelmholtzProjection', &
        'Using Piola-transformed finite elements', level=5)
  end subroutine hp_lookup

  subroutine hp_get_potsol(buf, nb) bind(c, name="hp_get_potsol")
    integer(c_int), value :: nb
    real(c_double) :: buf(1, nb)
    call GetVectorLocalSolution(buf(1:1, 1:nb), hp_potname, USolver=hp_pot)
  end subroutine hp_get_potsol

  function hp_nof_active() result(n) bind(c, name="hp_nof_active")
    integer(c_int) :: n
    n = GetNOFActive()
  end function hp_nof_active

  subroutine hp_default_initialize(ready) bind(c, name="hp_default_initialize")
    integer(c_int), value :: ready
    call DefaultInitialize(hp_solver, ready /= 0)
  end subroutine hp_default_initialize

  subroutine hp_elem_begin(t, n, npot, ndpot) bind(c, name="hp_elem_begin")
    integer(c_int), value :: t
    integer(c_int) :: n, npot, ndpot
    hp_elem => GetActiveElement(t)
    n = GetElementNOFNodes()
    npot = n * hp_pot % def_dofs(GetElementFamily(hp_elem), hp_elem % bodyid, 1)
    ndpot = GetElementNOFDOFs(USolver=hp_pot)
  end subroutine hp_elem_begin


  subroutine hp_nodes_touch() bind(c, name="hp_nodes_touch")
    call GetElementNodes(hp_nodes)
  end subroutine hp_nodes_touch

  function hp_ip_begin() result(n) bind(c, name="hp_ip_begin")
    integer(c_int) :: n
    integer :: degree
    if (hp_second) then
      degree = 2
      hp_ip = GaussPoints(hp_elem, EdgeBasis=.true., &
          PReferenceElement=hp_piola, EdgeBasisDegree=degree)
    else
      degree = 1
      hp_ip = GaussPoints(hp_elem, EdgeBasis=.true., &
          PReferenceElement=hp_piola)
    end if
    n = hp_ip % n
  end function hp_ip_begin

  subroutine hp_check_2d() bind(c, name="hp_check_2d")
    if (hp_dimv == 2 .and. .not. hp_piola) &
        call Fatal('HelmholtzProjection', '"Use Piola Transform = True" needed in 2D')
  end subroutine hp_check_2d

  subroutine hp_ip_point(t, u, v, w, s) bind(c, name="hp_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = hp_ip % u(t)
    v = hp_ip % v(t)
    w = hp_ip % w(t)
    s = hp_ip % s(t)
  end subroutine hp_ip_point

  subroutine hp_eleminfo(u, v, w, n, nd, ne, basis, dbasis, wbasis, curlw, &
      detj, stat) bind(c, name="hp_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: n, nd, ne
    real(c_double) :: basis(n), dbasis(n, 3), wbasis(ne, 3), curlw(ne, 3)
    real(c_double) :: detj
    integer(c_int) :: stat
    stat = merge(1, 0, ElementInfo(hp_elem, hp_nodes, u, v, w, detj, &
        basis(1:n), dbasis(1:n, 1:3), EdgeBasis=wbasis(1:ne, 1:3), &
        RotBasis=curlw(1:ne, 1:3), USolver=hp_pot))
  end subroutine hp_eleminfo

  ! One IP contribution (exact transcription of the update statements).
  subroutine hp_proj_point(ne, dim, npot, ndpot, potsol, sdet, n, ld, &
      dbasis, wbasis, stiff, force, ready) bind(c, name="hp_proj_point")
    integer(c_int), value :: ne, dim, npot, ndpot, n, ld, ready
    real(c_double) :: potsol(ndpot), dbasis(n, 3), wbasis(ne, 3)
    real(c_double) :: stiff(ld, ld), force(ld)
    real(c_double), value :: sdet
    real(c_double) :: a(1, 3)
    integer :: p, q
    a(1, 1:dim) = matmul(potsol(npot+1:ndpot), wbasis(1:ne, 1:dim))
    if (ready == 0) then
      do p = 1, n
        do q = 1, n
          stiff(p, q) = stiff(p, q) + sum(dbasis(q, 1:dim) * dbasis(p, 1:dim)) * sdet
        end do
      end do
    end if
    do p = 1, n
      force(p) = force(p) + sum(a(1, 1:dim) * dbasis(p, 1:dim)) * sdet
    end do
  end subroutine hp_proj_point

  subroutine hp_default_update_equations(ld, stiff, force) &
      bind(c, name="hp_default_update_equations")
    integer(c_int), value :: ld
    real(c_double) :: stiff(ld, ld), force(ld)
    call DefaultUpdateEquations(stiff(1:ld, 1:ld), force(1:ld))
  end subroutine hp_default_update_equations

  subroutine hp_default_update_force(ld, force) bind(c, name="hp_default_update_force")
    integer(c_int), value :: ld
    real(c_double) :: force(ld)
    call DefaultUpdateForce(force(1:ld))
  end subroutine hp_default_update_force

  subroutine hp_finish_bulk(cb, ready) bind(c, name="hp_finish_bulk")
    integer(c_int), value :: cb, ready
    if (cb /= 0) then
      call DefaultFinishBulkAssembly(bulkupdate=.not. (ready /= 0), rhsupdate=.false.)
    else
      call DefaultFinishBulkAssembly()
    end if
  end subroutine hp_finish_bulk

  subroutine hp_default_dirichlet() bind(c, name="hp_default_dirichlet")
    call DefaultDirichletBCs()
  end subroutine hp_default_dirichlet

  subroutine hp_zero_values() bind(c, name="hp_zero_values")
    hp_solver % variable % values = 0.0d0
  end subroutine hp_zero_values

  function hp_default_solve() result(norm) bind(c, name="hp_default_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function hp_default_solve

  subroutine hp_fatal_alloc() bind(c, name="hp_fatal_alloc")
    call Fatal('HelmholtzProjection', 'Memory allocation error.')
  end subroutine hp_fatal_alloc

  ! -- HelmholtzProjector_Init entry -----------------------------------
  subroutine hp_init_begin(h_solver) bind(c, name="hp_init_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, hp_solver)
  end subroutine hp_init_begin

  subroutine hp_init_body() bind(c, name="hp_init_body")
    type(valuelist_t), pointer :: solverpars
    solverpars => GetSolverParams()
    call ListAddLogical(solverpars, 'Linear System Refactorize', .false.)
  end subroutine hp_init_body

  ! -- RemoveKernelComponent entry -------------------------------------
  subroutine rk_begin(h_model, h_solver) bind(c, name="rk_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, rk_model)
    call c_f_pointer(h_solver, rk_solver)
  end subroutine rk_begin

  subroutine rk_default_start() bind(c, name="rk_default_start")
    call DefaultStart()
  end subroutine rk_default_start

  function rk_dim() result(d) bind(c, name="rk_dim")
    integer(c_int) :: d
    rk_dimv = CoordinateSystemDimension()
    d = rk_dimv
  end function rk_dim

  function rk_params_begin() result(n) bind(c, name="rk_params_begin")
    integer(c_int) :: n
    type(mesh_t), pointer :: m
    rk_params => GetSolverParams()
    m => GetMesh()
    n = m % maxelementdofs
  end function rk_params_begin

  function rk_constbulk() result(v) bind(c, name="rk_constbulk")
    integer(c_int) :: v
    logical :: found
    rk_cbflag = GetLogical(rk_params, 'Constant Bulk Matrix', found)
    v = merge(1, 0, rk_cbflag)
  end function rk_constbulk

  function rk_ready() result(v) bind(c, name="rk_ready")
    integer(c_int) :: v
    v = merge(1, 0, associated(rk_solver % matrix % bulkvalues) .and. rk_cbflag)
  end function rk_ready

  ! Both variable-solver lookups + FE-style inheritance (Fatal paths
  ! inside, exactly as in the original).
  subroutine rk_lookup() bind(c, name="rk_lookup")
    integer :: i, potdofs
    logical :: found
    rk_potname = GetString(rk_params, 'Potential Variable', found)
    if (.not. found) rk_potname = 'av'
    found = .false.
    do i = 1, rk_model % numberofsolvers
      rk_pot => rk_model % solvers(i)
      if (rk_potname == GetVarName(rk_pot % variable)) then
        found = .true.
        exit
      end if
    end do
    if (.not. found) then
      call Fatal('RemoveKernelComponent', 'Solver associated with potential variable > ' &
          //trim(rk_potname)//' < not found!')
    end if
    potdofs = rk_pot % variable % dofs
    if (potdofs > 1) call Fatal('RemoveKernelComponent', 'A real-valued potential expected')
    rk_kername = GetString(rk_params, 'Kernel Variable', found)
    if (.not. found) rk_kername = 'phi'
    found = .false.
    do i = 1, rk_model % numberofsolvers
      rk_ker => rk_model % solvers(i)
      if (rk_kername == GetVarName(rk_ker % variable)) then
        found = .true.
        exit
      end if
    end do
    if (.not. found) then
      call Fatal('RemoveKernelComponent', 'Solver associated with kernel variable > ' &
          //trim(rk_kername)//' < not found!')
    end if
    if (rk_ker % variable % dofs > 1) call Fatal('RemoveKernelComponent', &
        'A real-valued potential expected')
    call EdgeElementStyle(rk_pot % values, rk_piola, &
        QuadraticApproximation=rk_second)
    if (rk_piola) call Info('RemoveKernelComponent', &
        'Using Piola-transformed finite elements', level=5)
  end subroutine rk_lookup

  subroutine rk_get_potsol(buf, nb) bind(c, name="rk_get_potsol")
    integer(c_int), value :: nb
    real(c_double) :: buf(1, nb)
    call GetVectorLocalSolution(buf(1:1, 1:nb), rk_potname, USolver=rk_pot)
  end subroutine rk_get_potsol

  function rk_nof_active() result(n) bind(c, name="rk_nof_active")
    integer(c_int) :: n
    n = GetNOFActive()
  end function rk_nof_active

  subroutine rk_default_initialize(ready) bind(c, name="rk_default_initialize")
    integer(c_int), value :: ready
    call DefaultInitialize(rk_solver, ready /= 0)
  end subroutine rk_default_initialize

  subroutine rk_elem_begin(t, n, nd, npot, ndpot) bind(c, name="rk_elem_begin")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, npot, ndpot
    rk_elem => GetActiveElement(t)
    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
    npot = n * rk_pot % def_dofs(GetElementFamily(rk_elem), rk_elem % bodyid, 1)
    ndpot = GetElementNOFDOFs(USolver=rk_pot)
  end subroutine rk_elem_begin

  subroutine rk_fatal_dofcount() bind(c, name="rk_fatal_dofcount")
    call Fatal('RemoveKernelComponent', &
        'Potential variable DOFs count /= the solver DOFs count')
  end subroutine rk_fatal_dofcount


  subroutine rk_get_phisol(buf, nb) bind(c, name="rk_get_phisol")
    integer(c_int), value :: nb
    real(c_double) :: buf(nb)
    call GetScalarLocalSolution(buf(1:nb), rk_kername, USolver=rk_ker)
  end subroutine rk_get_phisol

  subroutine rk_nodes_touch() bind(c, name="rk_nodes_touch")
    call GetElementNodes(rk_nodes)
  end subroutine rk_nodes_touch

  function rk_ip_begin() result(n) bind(c, name="rk_ip_begin")
    integer(c_int) :: n
    if (rk_second) then
      rk_ip = GaussPoints(rk_elem, EdgeBasis=.true., &
          PReferenceElement=rk_piola, EdgeBasisDegree=2)
    else
      rk_ip = GaussPoints(rk_elem, EdgeBasis=.true., &
          PReferenceElement=rk_piola, EdgeBasisDegree=1)
    end if
    n = rk_ip % n
  end function rk_ip_begin

  subroutine rk_check_2d() bind(c, name="rk_check_2d")
    if (rk_dimv == 2 .and. .not. rk_piola) &
        call Fatal('RemoveKernelComponent', '"Use Piola Transform = True" needed in 2D')
  end subroutine rk_check_2d

  subroutine rk_ip_point(t, u, v, w, s) bind(c, name="rk_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = rk_ip % u(t)
    v = rk_ip % v(t)
    w = rk_ip % w(t)
    s = rk_ip % s(t)
  end subroutine rk_ip_point

  subroutine rk_eleminfo(u, v, w, n, nd, basis, dbasis, wbasis, curlw, &
      detj, stat) bind(c, name="rk_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: n, nd
    real(c_double) :: basis(n), dbasis(n, 3), wbasis(nd, 3), curlw(nd, 3)
    real(c_double) :: detj
    integer(c_int) :: stat
    stat = merge(1, 0, ElementInfo(rk_elem, rk_nodes, u, v, w, detj, &
        basis(1:n), dbasis(1:n, 1:3), EdgeBasis=wbasis(1:nd, 1:3), &
        RotBasis=curlw(1:nd, 1:3), USolver=rk_pot))
  end subroutine rk_eleminfo

  ! One IP contribution (exact transcription of the update statements).
  subroutine rk_point(nd, dim, n, phisol, sdet, dbasis, wbasis, stiff, &
      force, ready, ld) bind(c, name="rk_point")
    integer(c_int), value :: nd, dim, n, ready, ld
    real(c_double) :: phisol(n), dbasis(n, 3), wbasis(nd, 3)
    real(c_double) :: stiff(ld, ld), force(ld)
    real(c_double), value :: sdet
    real(c_double) :: a(3)
    integer :: i, p, q
    a = 0.0d0
    do i = 1, n
      a(1:dim) = a(1:dim) + phisol(i) * dbasis(i, 1:dim)
    end do
    if (ready == 0) then
      do p = 1, nd
        do q = 1, nd
          stiff(p, q) = stiff(p, q) + sum(wbasis(q, 1:dim) * wbasis(p, 1:dim)) * sdet
        end do
      end do
    end if
    do p = 1, nd
      force(p) = force(p) + sum(a(1:dim) * wbasis(p, 1:dim)) * sdet
    end do
  end subroutine rk_point

  subroutine rk_default_update_equations(ld, stiff, force) &
      bind(c, name="rk_default_update_equations")
    integer(c_int), value :: ld
    real(c_double) :: stiff(ld, ld), force(ld)
    call DefaultUpdateEquations(stiff(1:ld, 1:ld), force(1:ld))
  end subroutine rk_default_update_equations

  subroutine rk_default_update_force(ld, force) bind(c, name="rk_default_update_force")
    integer(c_int), value :: ld
    real(c_double) :: force(ld)
    call DefaultUpdateForce(force(1:ld))
  end subroutine rk_default_update_force

  subroutine rk_finish_bulk(cb, ready) bind(c, name="rk_finish_bulk")
    integer(c_int), value :: cb, ready
    if (cb /= 0) then
      call DefaultFinishBulkAssembly(bulkupdate=.not. (ready /= 0), rhsupdate=.false.)
    else
      call DefaultFinishBulkAssembly()
    end if
  end subroutine rk_finish_bulk

  subroutine rk_default_finish() bind(c, name="rk_default_finish")
    call DefaultFinishAssembly()
  end subroutine rk_default_finish

  subroutine rk_default_dirichlet() bind(c, name="rk_default_dirichlet")
    call DefaultDirichletBCs()
  end subroutine rk_default_dirichlet

  subroutine rk_zero_values() bind(c, name="rk_zero_values")
    rk_solver % variable % values = 0.0d0
  end subroutine rk_zero_values

  function rk_default_solve() result(norm) bind(c, name="rk_default_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function rk_default_solve

  ! -- final projection: redefines the potential variable ----------------
  function rk_perm_size(which) result(n) bind(c, name="rk_perm_size")
    integer(c_int), value :: which
    integer(c_int) :: n
    if (which == 1) then
      n = size(rk_solver % variable % perm(:))
    else
      n = size(rk_pot % variable % perm(:))
    end if
  end function rk_perm_size

  function rk_perm_val(which, i) result(v) bind(c, name="rk_perm_val")
    integer(c_int), value :: which, i
    integer(c_int) :: v
    if (which == 1) then
      v = rk_solver % variable % perm(i)
    else
      v = rk_pot % variable % perm(i)
    end if
  end function rk_perm_val

  function rk_values_get(which, j) result(v) bind(c, name="rk_values_get")
    integer(c_int), value :: which, j
    real(c_double) :: v
    if (which == 1) then
      v = rk_solver % variable % values(j)
    else
      v = rk_pot % variable % values(j)
    end if
  end function rk_values_get

  subroutine rk_values_set(which, j, v) bind(c, name="rk_values_set")
    integer(c_int), value :: which, j
    real(c_double), value :: v
    if (which == 1) then
      rk_solver % variable % values(j) = v
    else
      rk_pot % variable % values(j) = v
    end if
  end subroutine rk_values_set

  subroutine rk_fatal_nonmatch() bind(c, name="rk_fatal_nonmatch")
    call Fatal('RemoveKernelComponent', &
        'The variable and potential permutations are nonmatching?')
  end subroutine rk_fatal_nonmatch

  subroutine rk_fatal_differ() bind(c, name="rk_fatal_differ")
    call Fatal('RemoveKernelComponent', 'The variable and potential permutations differ')
  end subroutine rk_fatal_differ

  subroutine rk_fatal_alloc() bind(c, name="rk_fatal_alloc")
    call Fatal('RemoveKernelComponent', 'Memory allocation error.')
  end subroutine rk_fatal_alloc

  ! -- RemoveKernelComponent_Init0 entry --------------------------------
  subroutine rk_init0_begin(h_solver) bind(c, name="rk_init0_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, rk_solver)
  end subroutine rk_init0_begin

  subroutine rk_init0_body() bind(c, name="rk_init0_body")
    type(valuelist_t), pointer :: solverpars
    logical :: piolaversion, secondkind, secondorder
    solverpars => GetSolverParams()
    call ListAddLogical(solverpars, 'Linear System Refactorize', .false.)
    if (.not. ListCheckPresent(solverpars, "Element")) then
      call EdgeElementStyle(solverpars, piolaversion, secondkind, secondorder, &
          Check=.true.)
      if (secondorder) then
        call ListAddString(solverpars, "Element", &
            "n:0 e:2 -brick b:6 -pyramid b:3 -prism b:2 -quad_face b:4 -tri_face b:2")
      else if (secondkind) then
        call ListAddString(solverpars, "Element", "n:0 e:2")
      else if (piolaversion) then
        call ListAddString(solverpars, "Element", &
            "n:0 e:1 -brick b:3 -quad_face b:2")
      else
        call ListAddString(solverpars, "Element", "n:0 e:1")
      end if
    end if
  end subroutine rk_init0_body

end module helmholtzprojection_acc
