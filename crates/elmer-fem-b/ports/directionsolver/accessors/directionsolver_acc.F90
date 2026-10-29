! T022-A DirectionSolver port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/DirectionSolver.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, whole-model surgery
! (ALLOCATE/pointer reassignment of Model % Solvers) — stays inside
! Fortran, so the Rust side is pure C calls. The wrapper module mirrors
! DirectionSolver's own locals (current element, LocalMatrix and
! BoundaryCondition SAVE Nodes and IP locals, SaveElementSolution SAVE
! state) with identical lifetimes. Order-sensitive reductions (MATMUL,
! SUM) and the DG-indexed save stay in Fortran (same compiler), so
! bit-identity holds by construction wherever the call sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module directionsolver_acc
  use iso_c_binding
  use types
  use defutils
  use integration
  use matrixassembly
  use eleminfo
  use messages
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(mesh_t), pointer, save :: acc_mesh => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(valuelist_t), pointer, save :: acc_bf => null()
  type(valuelist_t), pointer, save :: acc_bc => null()
  ! LocalMatrix context (mirrors its SAVE Nodes and IP local).
  type(nodes_t), target, save :: acc_nodes_lm
  type(gaussintegrationpoints_t), save :: acc_ip_lm
  ! BoundaryCondition context (mirrors its SAVE Nodes and IP local).
  type(nodes_t), target, save :: acc_nodes_bc
  type(gaussintegrationpoints_t), save :: acc_ip_bc
  ! SaveElementSolution context (mirrors its SAVE First/varname/
  ! directionvar). acc_savename starts undefined, exactly like the
  ! original's SAVE varname, so the first-visit compare behaves the same.
  type(variable_t), pointer, save :: acc_dirvar => null()
  logical, save :: acc_first = .true.
  character(len=MAX_NAME_LEN), save :: acc_savename

contains

  ! -- entry state ----------------------------------------------------
  subroutine dsa_begin(h_model, h_solver) bind(c, name="dsa_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine dsa_begin

  function dsa_assoc_matrix() result(r) bind(c, name="dsa_assoc_matrix")
    integer(c_int) :: r
    r = merge(1, 0, associated(acc_solver % matrix))
  end function dsa_assoc_matrix

  subroutine dsa_getmesh_touch() bind(c, name="dsa_getmesh_touch")
    acc_mesh => GetMesh()
  end subroutine dsa_getmesh_touch

  function dsa_max_elem_dofs() result(n) bind(c, name="dsa_max_elem_dofs")
    integer(c_int) :: n
    n = acc_solver % mesh % maxelementdofs
  end function dsa_max_elem_dofs

  function dsa_number_of_bodies() result(n) &
      bind(c, name="dsa_number_of_bodies")
    integer(c_int) :: n
    n = acc_model % numberofbodies
  end function dsa_number_of_bodies

  subroutine dsa_fatal_alloc() bind(c, name="dsa_fatal_alloc")
    call Fatal('DirectionSolver', 'Memory allocation error.')
  end subroutine dsa_fatal_alloc

  ! -- DirectionSolver_Init0 body --------------------------------------
  ! DG guard: mirrors `IF (GetLogical(GetSolverParams(),
  ! 'Discontinuous Galerkin',Found)) RETURN`.
  function dsa_init0_dg() result(r) bind(c, name="dsa_init0_dg")
    integer(c_int) :: r
    logical :: found
    r = merge(1, 0, GetLogical(GetSolverParams(), &
        'Discontinuous Galerkin', found))
  end function dsa_init0_dg

  ! Model surgery: append the never-executed DG dummy solver and widen
  ! the equation Active Solvers lists. Mirrors the Init0 statements in
  ! order (same List call sequence, same intrinsic array copy, same
  ! DEALLOCATE/pointer reassignment). The write-only `visited` counter
  ! of the original is not replicated (no observable state).
  subroutine dsa_init0_body() bind(c, name="dsa_init0_body")
    type(solver_t), pointer :: ps, newsolvers(:)
    type(valuelist_t), pointer :: spar
    integer, pointer :: act(:)
    character(len=MAX_NAME_LEN) :: vname
    integer :: mine, i, m, n
    logical :: found
    ps => acc_solver
    do mine = 1, acc_model % NumberOfSolvers
      if (associated(ps, acc_model % Solvers(mine))) exit
    end do
    vname = GetString(GetSolverParams(), 'Variable', found)
    n = acc_model % NumberOfSolvers
    do i = 1, acc_model % NumberOfEquations
      act => ListGetIntegerArray(acc_model % Equations(i) % Values, &
          'Active Solvers', found)
      m = size(act)
      if (any(act == mine)) call ListAddIntegerArray( &
          acc_model % Equations(i) % Values, &
          'Active Solvers', m + 1, [act, n + 1])
    end do
    allocate(newsolvers(n + 1))
    newsolvers(1:n) = acc_model % Solvers
    newsolvers(n + 1) % Values => ListAllocate()
    spar => newsolvers(n + 1) % Values
    call ListAddLogical(spar, 'Discontinuous Galerkin', .true.)
    newsolvers(n + 1) % DG = .true.
    newsolvers(n + 1) % PROCEDURE = C_NULL_FUNPTR
    newsolvers(n + 1) % ActiveElements => null()
    call ListAddString(spar, 'Exec Solver', 'never')
    call ListAddLogical(spar, 'No Matrix', .true.)
    call ListAddLogical(spar, 'Optimize Bandwidth', .false.)
    call ListAddString(spar, 'Equation', 'elementaladd'//trim(vname))
    call ListAddString(spar, 'Procedure', &
        'DirectionSolver DirectionSolver_Dummy', .false.)
    call ListAddString(spar, 'Variable', &
        '-nooutput '//trim(vname)//'_dummy')
    i = 1
    do while (.true.)
      if (ListCheckPresent(spar, 'Exported Variable '//i2s(i))) then
        i = i + 1
      else
        exit
      end if
    end do
    call ListAddString(spar, 'Exported Variable '//i2s(i), &
        trim(vname)//' Direction')
    deallocate(acc_model % Solvers)
    acc_model % Solvers => newsolvers
    acc_model % NumberOfSolvers = n + 1
  end subroutine dsa_init0_body

  ! -- per-body namespace prologue/epilogue ----------------------------
  subroutine dsa_ns_setup(ns) bind(c, name="dsa_ns_setup")
    integer(c_int), value :: ns
    character(len=MAX_NAME_LEN) :: nsspace, vnwithns, vname
    logical :: found
    vname = GetString(GetSolverParams(), 'Variable', found)
    nsspace = 'body '//i2s(ns)//':'
    vnwithns = trim(nsspace)//' '//trim(vname)
    if (ListCheckPresentAnyBC(acc_model, vnwithns)) &
        call ListSetNameSpace(nsspace)
  end subroutine dsa_ns_setup

  subroutine dsa_ns_reset() bind(c, name="dsa_ns_reset")
    call ListSetNameSpace('')
  end subroutine dsa_ns_reset

  ! -- assembly driver -------------------------------------------------
  function dsa_nof_active() result(n) bind(c, name="dsa_nof_active")
    integer(c_int) :: n
    n = GetNOFActive()
  end function dsa_nof_active

  subroutine dsa_default_initialize() bind(c, name="dsa_default_initialize")
    call DefaultInitialize()
  end subroutine dsa_default_initialize

  ! -- active-element context ------------------------------------------
  subroutine dsa_begin_active(t, n, nd, nb) bind(c, name="dsa_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, nb
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
    nb = GetElementNOFBDOFs()
  end subroutine dsa_begin_active

  function dsa_bodyforce() result(r) bind(c, name="dsa_bodyforce")
    integer(c_int) :: r
    acc_bf => GetBodyForce()
    r = merge(1, 0, associated(acc_bf))
  end function dsa_bodyforce

  subroutine dsa_load_source(n, load, found) bind(c, name="dsa_load_source")
    integer(c_int), value :: n
    real(c_double) :: load(n)
    integer(c_int) :: found
    logical :: f
    load(1:n) = GetReal(acc_bf, 'Source', f)
    found = merge(1, 0, f)
  end subroutine dsa_load_source

  ! -- boundary-element context ----------------------------------------
  function dsa_num_bnd_elems() result(n) bind(c, name="dsa_num_bnd_elems")
    integer(c_int) :: n
    n = acc_solver % mesh % numberofboundaryelements
  end function dsa_num_bnd_elems

  subroutine dsa_begin_boundary(t, n, active, family) &
      bind(c, name="dsa_begin_boundary")
    integer(c_int), value :: t
    integer(c_int) :: n, active, family
    acc_elem => GetBoundaryElement(t)
    active = merge(1, 0, ActiveBoundaryElement())
    n = GetElementNOFNodes()
    family = GetElementFamily()
  end subroutine dsa_begin_boundary

  subroutine dsa_bc() bind(c, name="dsa_bc")
    acc_bc => GetBC()
  end subroutine dsa_bc

  subroutine dsa_bc_logical(which, val, found) bind(c, name="dsa_bc_logical")
    integer(c_int), value :: which
    integer(c_int) :: val, found
    logical :: l, f
    if (which == 1) then
      l = GetLogical(acc_bc, 'Positive Electrode', f)
    else
      l = GetLogical(acc_bc, 'Negative Electrode', f)
    end if
    val = merge(1, 0, l)
    found = merge(1, 0, f)
  end subroutine dsa_bc_logical

  function dsa_var_is_w() result(r) bind(c, name="dsa_var_is_w")
    integer(c_int) :: r
    r = merge(1, 0, acc_solver % variable % name == 'w')
  end function dsa_var_is_w

  ! -- LocalMatrix context ----------------------------------------------
  subroutine dsa_lm_nodes() bind(c, name="dsa_lm_nodes")
    call GetElementNodes(acc_nodes_lm)
  end subroutine dsa_lm_nodes

  function dsa_lm_ip_begin() result(n) bind(c, name="dsa_lm_ip_begin")
    integer(c_int) :: n
    acc_ip_lm = GaussPoints(acc_elem)
    n = acc_ip_lm % n
  end function dsa_lm_ip_begin

  subroutine dsa_lm_ip_point(t, u, v, w, s) bind(c, name="dsa_lm_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip_lm % u(t)
    v = acc_ip_lm % v(t)
    w = acc_ip_lm % w(t)
    s = acc_ip_lm % s(t)
  end subroutine dsa_lm_ip_point

  subroutine dsa_lm_eleminfo(u, v, w, n, nd, basis, dbasis, detj, stat) &
      bind(c, name="dsa_lm_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: n, nd
    real(c_double) :: basis(n), dbasis(nd, 3), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_elem, acc_nodes_lm, u, v, w, detj, &
        basis(1:n), dbasis(1:nd, 1:3))
    stat = merge(1, 0, ok)
  end subroutine dsa_lm_eleminfo

  function dsa_dot(n, a, b) result(s) bind(c, name="dsa_dot")
    integer(c_int), value :: n
    real(c_double) :: a(n), b(n), s
    s = sum(a(1:n) * b(1:n))
  end function dsa_dot

  subroutine dsa_lm_assemble(nd, dbasis, sdet, basis, loadatip, stiff, force) &
      bind(c, name="dsa_lm_assemble")
    integer(c_int), value :: nd
    real(c_double) :: dbasis(nd, 3), basis(nd), stiff(nd, nd), force(nd)
    real(c_double), value :: sdet, loadatip
    ! Exact transcription of the two update statements (same compiler,
    ! same evaluation order as the original internal procedure).
    stiff(1:nd, 1:nd) = stiff(1:nd, 1:nd) + sdet * &
        matmul(dbasis, transpose(dbasis))
    force(1:nd) = force(1:nd) + sdet * loadatip * basis(1:nd)
  end subroutine dsa_lm_assemble

  ! -- BoundaryCondition context ----------------------------------------
  subroutine dsa_bc_nodes() bind(c, name="dsa_bc_nodes")
    call GetElementNodes(acc_nodes_bc)
  end subroutine dsa_bc_nodes

  function dsa_bc_ip_begin() result(n) bind(c, name="dsa_bc_ip_begin")
    integer(c_int) :: n
    acc_ip_bc = GaussPoints(acc_elem)
    n = acc_ip_bc % n
  end function dsa_bc_ip_begin

  subroutine dsa_bc_ip_point(t, u, v, w, s) bind(c, name="dsa_bc_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip_bc % u(t)
    v = acc_ip_bc % v(t)
    w = acc_ip_bc % w(t)
    s = acc_ip_bc % s(t)
  end subroutine dsa_bc_ip_point

  subroutine dsa_bc_eleminfo(u, v, w, n, nd, basis, dbasis, detj, stat) &
      bind(c, name="dsa_bc_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: n, nd
    real(c_double) :: basis(n), dbasis(nd, 3), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_elem, acc_nodes_bc, u, v, w, detj, &
        basis(1:n), dbasis(1:nd, 1:3))
    stat = merge(1, 0, ok)
  end subroutine dsa_bc_eleminfo

  subroutine dsa_bc_force_update(n, sdet, loadatip, basis, force) &
      bind(c, name="dsa_bc_force_update")
    integer(c_int), value :: n
    real(c_double) :: basis(n), force(n)
    real(c_double), value :: sdet, loadatip
    integer :: j
    ! Exact transcription of the per-node update loop.
    do j = 1, n
      force(j) = force(j) + sdet * loadatip * basis(j)
    end do
  end subroutine dsa_bc_force_update

  ! -- shared assembly tail ---------------------------------------------
  subroutine dsa_condensate(nd, nb, ld, k, f) bind(c, name="dsa_condensate")
    integer(c_int), value :: nd, nb, ld
    real(c_double) :: k(ld, ld), f(ld)
    call CondensateP(nd, nb, k(1:ld, 1:ld), f(1:ld))
  end subroutine dsa_condensate

  subroutine dsa_default_update(ld, g, f) bind(c, name="dsa_default_update")
    integer(c_int), value :: ld
    real(c_double) :: g(ld, ld), f(ld)
    call DefaultUpdateEquations(g(1:ld, 1:ld), f(1:ld))
  end subroutine dsa_default_update

  subroutine dsa_default_finish() bind(c, name="dsa_default_finish")
    call DefaultFinishAssembly()
  end subroutine dsa_default_finish

  subroutine dsa_default_dirichlet() bind(c, name="dsa_default_dirichlet")
    call DefaultDirichletBCs()
  end subroutine dsa_default_dirichlet

  function dsa_default_solve() result(norm) bind(c, name="dsa_default_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function dsa_default_solve

  ! -- per-body save (SaveSolutionWithBodyMethod + SaveElementSolution) --
  subroutine dsa_save_sweep(ns) bind(c, name="dsa_save_sweep")
    integer(c_int), value :: ns
    integer :: active, t, nn
    active = GetNOFActive()
    do t = 1, active
      acc_elem => GetActiveElement(t)
      if (ns /= acc_elem % BodyId) cycle
      nn = GetElementNOFNodes()
      call dsa_save_elem(nn)
    end do
  end subroutine dsa_save_sweep

  subroutine dsa_save_elem(nn) bind(c, name="dsa_save_elem")
    integer(c_int), value :: nn
    real(c_double) :: direction(nn)
    integer :: j, k
    logical :: found
    if (acc_savename /= GetString(GetSolverParams(), 'Variable', found)) &
        acc_first = .true.
    acc_savename = GetString(GetSolverParams(), 'Variable', found)
    if (acc_first) then
      acc_first = .false.
      acc_dirvar => VariableGet(acc_mesh % Variables, &
          trim(acc_savename)//' Direction')
      if (.not. associated(acc_dirvar)) &
          call Fatal('SaveElementSolution()', 'Direction variable not found')
    end if
    call GetScalarLocalSolution(direction, acc_savename)
    do j = 1, nn
      if (associated(acc_dirvar)) then
        do k = 1, acc_dirvar % DOFs
          acc_dirvar % Values(acc_dirvar % DOFs*(acc_dirvar % Perm( &
              acc_elem % DGIndexes(j)) - 1) + k) = direction(j)
        end do
      end if
    end do
  end subroutine dsa_save_elem

end module directionsolver_acc
