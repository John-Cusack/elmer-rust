! T020 Poisson pilot: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/Poisson.F90, compiled against the built .mod
! files of the target configure (A2 shipped mechanism, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module mirrors Poisson's own locals (current element, Nodes SAVE state
! per internal procedure, IP shallow copies) with identical lifetimes.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module poisson_acc
  use iso_c_binding
  use types
  use defutils
  use integration
  use matrixassembly
  use eleminfo
  use messages
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(valuelist_t), pointer, save :: acc_bf => null()
  type(valuelist_t), pointer, save :: acc_bc => null()
  ! LocalMatrix context (mirrors its SAVE Nodes and IP local).
  type(nodes_t), target, save :: acc_nodes_lm
  type(gaussintegrationpoints_t), save :: acc_ip_lm
  ! BoundaryCondition context (mirrors its SAVE Nodes and IP local).
  type(nodes_t), target, save :: acc_nodes_bc
  type(gaussintegrationpoints_t), save :: acc_ip_bc

contains

  ! -- entry state ----------------------------------------------------
  subroutine poa_begin(h_solver) bind(c, name="poa_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine poa_begin

  function poa_assoc_matrix() result(r) bind(c, name="poa_assoc_matrix")
    integer(c_int) :: r
    r = merge(1, 0, associated(acc_solver % matrix))
  end function poa_assoc_matrix

  function poa_max_elem_dofs() result(n) bind(c, name="poa_max_elem_dofs")
    integer(c_int) :: n
    n = acc_solver % mesh % maxelementdofs
  end function poa_max_elem_dofs

  function poa_num_bnd_elems() result(n) bind(c, name="poa_num_bnd_elems")
    integer(c_int) :: n
    n = acc_solver % mesh % numberofboundaryelements
  end function poa_num_bnd_elems

  function poa_var_is_w() result(r) bind(c, name="poa_var_is_w")
    integer(c_int) :: r
    r = merge(1, 0, acc_solver % variable % name == 'w')
  end function poa_var_is_w

  ! -- assembly driver -------------------------------------------------
  function poa_nof_active() result(n) bind(c, name="poa_nof_active")
    integer(c_int) :: n
    n = GetNOFActive()
  end function poa_nof_active

  subroutine poa_default_initialize() bind(c, name="poa_default_initialize")
    call DefaultInitialize()
  end subroutine poa_default_initialize

  subroutine poa_getmesh_touch() bind(c, name="poa_getmesh_touch")
    type(mesh_t), pointer :: m
    m => GetMesh()
  end subroutine poa_getmesh_touch

  ! -- active-element context ------------------------------------------
  subroutine poa_begin_active(t, n, nd, nb) bind(c, name="poa_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, nb
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
    nb = GetElementNOFBDOFs()
  end subroutine poa_begin_active

  function poa_bodyforce() result(r) bind(c, name="poa_bodyforce")
    integer(c_int) :: r
    acc_bf => GetBodyForce()
    r = merge(1, 0, associated(acc_bf))
  end function poa_bodyforce

  subroutine poa_load_source(n, load, found) bind(c, name="poa_load_source")
    integer(c_int), value :: n
    real(c_double) :: load(n)
    integer(c_int) :: found
    logical :: f
    load(1:n) = GetReal(acc_bf, 'Source', f)
    found = merge(1, 0, f)
  end subroutine poa_load_source

  ! -- boundary-element context ----------------------------------------
  subroutine poa_begin_boundary(t, n, active, family) bind(c, name="poa_begin_boundary")
    integer(c_int), value :: t
    integer(c_int) :: n, active, family
    acc_elem => GetBoundaryElement(t)
    active = merge(1, 0, ActiveBoundaryElement())
    n = GetElementNOFNodes()
    family = GetElementFamily()
  end subroutine poa_begin_boundary

  subroutine poa_bc() bind(c, name="poa_bc")
    acc_bc => GetBC()
  end subroutine poa_bc

  subroutine poa_bc_logical(which, val, found) bind(c, name="poa_bc_logical")
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
  end subroutine poa_bc_logical

  ! -- LocalMatrix context ----------------------------------------------
  subroutine poa_lm_nodes() bind(c, name="poa_lm_nodes")
    call GetElementNodes(acc_nodes_lm)
  end subroutine poa_lm_nodes

  function poa_lm_ip_begin() result(n) bind(c, name="poa_lm_ip_begin")
    integer(c_int) :: n
    acc_ip_lm = GaussPoints(acc_elem)
    n = acc_ip_lm % n
  end function poa_lm_ip_begin

  subroutine poa_lm_ip_point(t, u, v, w, s) bind(c, name="poa_lm_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip_lm % u(t)
    v = acc_ip_lm % v(t)
    w = acc_ip_lm % w(t)
    s = acc_ip_lm % s(t)
  end subroutine poa_lm_ip_point

  subroutine poa_lm_eleminfo(u, v, w, n, nd, basis, dbasis, detj, stat) &
      bind(c, name="poa_lm_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: n, nd
    real(c_double) :: basis(n), dbasis(nd, 3), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_elem, acc_nodes_lm, u, v, w, detj, &
        basis(1:n), dbasis(1:nd, 1:3))
    stat = merge(1, 0, ok)
  end subroutine poa_lm_eleminfo

  function poa_dot(n, a, b) result(s) bind(c, name="poa_dot")
    integer(c_int), value :: n
    real(c_double) :: a(n), b(n), s
    s = sum(a(1:n) * b(1:n))
  end function poa_dot

  subroutine poa_lm_assemble(nd, dbasis, sdet, basis, loadatip, stiff, force) &
      bind(c, name="poa_lm_assemble")
    integer(c_int), value :: nd
    real(c_double) :: dbasis(nd, 3), basis(nd), stiff(nd, nd), force(nd)
    real(c_double), value :: sdet, loadatip
    ! Exact transcription of the two update statements (same compiler,
    ! same evaluation order as the original internal procedure).
    stiff(1:nd, 1:nd) = stiff(1:nd, 1:nd) + sdet * &
        matmul(dbasis, transpose(dbasis))
    force(1:nd) = force(1:nd) + sdet * loadatip * basis(1:nd)
  end subroutine poa_lm_assemble

  ! -- BoundaryCondition context ----------------------------------------
  subroutine poa_bc_nodes() bind(c, name="poa_bc_nodes")
    call GetElementNodes(acc_nodes_bc)
  end subroutine poa_bc_nodes

  function poa_bc_ip_begin() result(n) bind(c, name="poa_bc_ip_begin")
    integer(c_int) :: n
    acc_ip_bc = GaussPoints(acc_elem)
    n = acc_ip_bc % n
  end function poa_bc_ip_begin

  subroutine poa_bc_ip_point(t, u, v, w, s) bind(c, name="poa_bc_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip_bc % u(t)
    v = acc_ip_bc % v(t)
    w = acc_ip_bc % w(t)
    s = acc_ip_bc % s(t)
  end subroutine poa_bc_ip_point

  subroutine poa_bc_eleminfo(u, v, w, n, nd, basis, dbasis, detj, stat) &
      bind(c, name="poa_bc_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: n, nd
    real(c_double) :: basis(n), dbasis(nd, 3), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_elem, acc_nodes_bc, u, v, w, detj, &
        basis(1:n), dbasis(1:nd, 1:3))
    stat = merge(1, 0, ok)
  end subroutine poa_bc_eleminfo

  subroutine poa_bc_force_update(n, sdet, loadatip, basis, force) &
      bind(c, name="poa_bc_force_update")
    integer(c_int), value :: n
    real(c_double) :: basis(n), force(n)
    real(c_double), value :: sdet, loadatip
    integer :: j
    ! Exact transcription of the per-node update loop.
    do j = 1, n
      force(j) = force(j) + sdet * loadatip * basis(j)
    end do
  end subroutine poa_bc_force_update

  ! -- shared assembly tail ---------------------------------------------
  subroutine poa_condensate(nd, nb, ld, k, f) bind(c, name="poa_condensate")
    integer(c_int), value :: nd, nb, ld
    real(c_double) :: k(ld, ld), f(ld)
    call CondensateP(nd, nb, k(1:ld, 1:ld), f(1:ld))
  end subroutine poa_condensate

  subroutine poa_default_update(ld, g, f) bind(c, name="poa_default_update")
    integer(c_int), value :: ld
    real(c_double) :: g(ld, ld), f(ld)
    call DefaultUpdateEquations(g(1:ld, 1:ld), f(1:ld))
  end subroutine poa_default_update

  subroutine poa_default_finish() bind(c, name="poa_default_finish")
    call DefaultFinishAssembly()
  end subroutine poa_default_finish

  subroutine poa_default_dirichlet() bind(c, name="poa_default_dirichlet")
    call DefaultDirichletBCs()
  end subroutine poa_default_dirichlet

  function poa_default_solve() result(norm) bind(c, name="poa_default_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function poa_default_solve

  subroutine poa_fatal_alloc() bind(c, name="poa_fatal_alloc")
    call Fatal('PoissonSolve', 'Memory allocation error.')
  end subroutine poa_fatal_alloc

end module poisson_acc
