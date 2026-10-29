! T022-A DistributeSource port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/DistributeSource.F90, compiled against
! the built .mod files of the target configure (ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) /
! INTEGER(C_INT) arrays). All Fortran-object traffic — TYPE(...) handles,
! OPTIONAL/POINTER dummies, assumed-shape descriptors, hidden lengths,
! array-pointer results — stays inside Fortran, so the Rust side is pure
! C calls. The wrapper module mirrors the original's SAVE'd Nodes state
! (ElementNodes, ElementNodesB) with identical lifetimes, plus the
! internal-procedure locals (integration points, element/leaf handles,
! solver/variable pointers) the Rust driver threads through.
!
! Order-sensitive reductions (SUM, MATMUL), formatted I/O (PRINT) and
! whole-array updates stay in Fortran wrappers, so bit-identity holds by
! construction wherever the call sequence matches. No upstream text:
! clean-room wrappers around Elmer library calls.
module distributesource_acc
  use iso_c_binding
  use types
  use defutils
  use integration
  use messages
  implicit none

  character(*), parameter :: caller = 'IntegrateSource'
  character(*), parameter :: src_name = 'Local Source'
  character(*), parameter :: rhs_name = 'g'
  character(*), parameter :: fix_name = 'g fix'

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(solver_t), pointer, save :: acc_gsolver => null()
  type(variable_t), pointer, save :: acc_rhs => null()
  type(variable_t), pointer, save :: acc_fix => null()
  type(mesh_t), pointer, save :: acc_local => null()
  type(mesh_t), pointer, save :: acc_global => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(element_t), pointer, save :: acc_elemb => null()
  type(quadrant_t), pointer, save :: acc_leaf => null()
  type(valuelist_t), pointer, save :: acc_bf => null()
  ! SAVE mirrors of the internal procedure's SAVE Nodes state. Allocated
  ! once at first-mesh sizing and grown only if a later mesh needs more;
  ! the leading sections used below are unaffected, so valid inputs never
  ! observe the growth.
  type(nodes_t), save :: acc_nodes
  type(nodes_t), save :: acc_nodesb
  type(gaussintegrationpoints_t), save :: acc_ip
  integer, save :: acc_degree = 0
  logical, save :: acc_piola = .false.

contains

  ! -- entry state -----------------------------------------------------
  subroutine dsa_begin(h_model, h_solver) bind(c, name="dsa_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine dsa_begin

  function dsa_meshes_differ() result(r) bind(c, name="dsa_meshes_differ")
    integer(c_int) :: r
    acc_local => GetMesh()
    acc_global => acc_model % meshes
    r = merge(0, 1, associated(acc_local, acc_global))
  end function dsa_meshes_differ

  subroutine dsa_fatal_same_mesh() bind(c, name="dsa_fatal_same_mesh")
    call Fatal('DistributeSource', 'Meshes are the same!')
  end subroutine dsa_fatal_same_mesh

  ! -- setup ------------------------------------------------------------
  subroutine dsa_info_projecting() bind(c, name="dsa_info_projecting")
    call Info(caller, 'Projecting current sources between meshes')
  end subroutine dsa_info_projecting

  subroutine dsa_edge_param(val, found) bind(c, name="dsa_edge_param")
    integer(c_int) :: val, found
    type(valuelist_t), pointer :: params
    logical :: l, f
    params => GetSolverParams(acc_solver)
    l = GetLogical(params, 'Edge Basis', f)
    val = merge(1, 0, l)
    found = merge(1, 0, f)
  end subroutine dsa_edge_param

  function dsa_dim() result(d) bind(c, name="dsa_dim")
    integer(c_int) :: d
    d = CoordinateSystemDimension()
  end function dsa_dim

  function dsa_ensure_qtree() result(t) bind(c, name="dsa_ensure_qtree")
    integer(c_int) :: t
    type(quadrant_t), pointer :: root
    real(kind=dp) :: bbox(6), eps2
    logical :: useqtree
    useqtree = .true.
    root => acc_global % rootquadrant
    if (.not. associated(root)) then
      bbox(1) = minval(acc_global % nodes % x)
      bbox(2) = minval(acc_global % nodes % y)
      bbox(3) = minval(acc_global % nodes % z)
      bbox(4) = maxval(acc_global % nodes % x)
      bbox(5) = maxval(acc_global % nodes % y)
      bbox(6) = maxval(acc_global % nodes % z)
      eps2 = 0.1_dp * maxval(bbox(4:6) - bbox(1:3))
      bbox(1:3) = bbox(1:3) - eps2
      bbox(4:6) = bbox(4:6) + eps2
      call BuildQuadrantTree(acc_global, bbox, acc_global % rootquadrant)
      root => acc_global % rootquadrant
      call Info(caller, 'Quadrant tree build for interpolation', level=10)
    end if
    t = merge(1, 0, associated(root) .and. useqtree)
  end function dsa_ensure_qtree

  function dsa_max_dofs() result(n) bind(c, name="dsa_max_dofs")
    integer(c_int) :: n
    n = max(acc_global % maxelementdofs, acc_local % maxelementdofs)
  end function dsa_max_dofs

  subroutine dsa_alloc_nodes(n) bind(c, name="dsa_alloc_nodes")
    integer(c_int), value :: n
    if (associated(acc_nodes % x)) then
      if (size(acc_nodes % x) >= n) return
      deallocate(acc_nodes % x, acc_nodes % y, acc_nodes % z)
      deallocate(acc_nodesb % x, acc_nodesb % y, acc_nodesb % z)
    end if
    allocate(acc_nodes % x(n), acc_nodes % y(n), acc_nodes % z(n))
    allocate(acc_nodesb % x(n), acc_nodesb % y(n), acc_nodesb % z(n))
  end subroutine dsa_alloc_nodes

  subroutine dsa_sim_params(geps, leps, tries, neps) bind(c, name="dsa_sim_params")
    real(c_double) :: geps, leps, neps
    integer(c_int) :: tries
    logical :: stat
    geps = ListGetConstReal(acc_model % simulation, 'Interpolation Global Epsilon', stat)
    if (.not. stat) geps = 2.0d-10
    leps = ListGetConstReal(acc_model % simulation, 'Interpolation Local Epsilon', stat)
    if (.not. stat) leps = 1.0d-10
    tries = ListGetInteger(acc_model % simulation, 'Interpolation Max Iterations', stat)
    if (.not. stat) tries = 12
    neps = ListGetConstReal(acc_model % simulation, 'Interpolation Numeric Epsilon', stat)
    if (.not. stat) neps = 1.0e-10
  end subroutine dsa_sim_params

  function dsa_find_rhs() result(r) bind(c, name="dsa_find_rhs")
    integer(c_int) :: r
    nullify(acc_gsolver)
    nullify(acc_fix)
    acc_rhs => VariableGet(acc_global % variables, rhs_name, .true.)
    r = merge(1, 0, associated(acc_rhs))
  end function dsa_find_rhs

  subroutine dsa_fatal_no_rhs() bind(c, name="dsa_fatal_no_rhs")
    call Fatal(caller, 'Could not find variable: ' // trim(rhs_name))
  end subroutine dsa_fatal_no_rhs

  function dsa_assoc_gsolver() result(r) bind(c, name="dsa_assoc_gsolver")
    integer(c_int) :: r
    integer :: i
    nullify(acc_gsolver)
    do i = 1, acc_model % numberofsolvers
      if (associated(acc_rhs % solver, acc_model % solvers(i))) then
        call Info(caller, 'Source will be projected to solver id: ' // i2s(i), level=8)
        acc_gsolver => acc_model % solvers(i)
        exit
      end if
    end do
    r = merge(1, 0, associated(acc_gsolver))
  end function dsa_assoc_gsolver

  subroutine dsa_fatal_no_gsolver() bind(c, name="dsa_fatal_no_gsolver")
    call Fatal(caller, 'Could not associate r.h.s. vector to any solver!')
  end subroutine dsa_fatal_no_gsolver

  function dsa_edge_detect() result(r) bind(c, name="dsa_edge_detect")
    integer(c_int) :: r
    integer :: n, m
    logical :: edge
    edge = .false.
    n = acc_global % numberofnodes
    m = size(acc_rhs % perm)
    if (m > n) edge = any(acc_rhs % perm(n + 1:m) > 0)
    r = merge(1, 0, edge)
  end function dsa_edge_detect

  subroutine dsa_info_edge() bind(c, name="dsa_info_edge")
    call Info(caller, 'Projecting source to edge basis functions')
  end subroutine dsa_info_edge

  function dsa_find_fix() result(r) bind(c, name="dsa_find_fix")
    integer(c_int) :: r
    nullify(acc_fix)
    acc_fix => VariableGet(acc_global % variables, fix_name, .true.)
    r = merge(1, 0, associated(acc_fix))
  end function dsa_find_fix

  subroutine dsa_warn_no_fix() bind(c, name="dsa_warn_no_fix")
    call Warn(caller, 'Could not find variable: ' // trim(fix_name))
  end subroutine dsa_warn_no_fix

  subroutine dsa_edge_style(degree) bind(c, name="dsa_edge_style")
    integer(c_int) :: degree
    acc_piola = .false.
    call EdgeElementStyle(acc_gsolver % values, acc_piola, basisdegree=acc_degree)
    degree = acc_degree
  end subroutine dsa_edge_style

  subroutine dsa_fatal_no_edge() bind(c, name="dsa_fatal_no_edge")
    call Fatal(caller, 'Currently we assume edge basis!')
  end subroutine dsa_fatal_no_edge

  subroutine dsa_zero_values() bind(c, name="dsa_zero_values")
    acc_rhs % values = 0.0_dp
    if (associated(acc_fix)) acc_fix % values = 0.0_dp
  end subroutine dsa_zero_values

  subroutine dsa_info_integrating() bind(c, name="dsa_info_integrating")
    call Info(caller, 'Integrating over bulk elements of local mesh', level=8)
  end subroutine dsa_info_integrating

  ! -- local-mesh traversal ----------------------------------------------
  function dsa_local_num_bulk() result(n) bind(c, name="dsa_local_num_bulk")
    integer(c_int) :: n
    n = acc_local % numberofbulkelements
  end function dsa_local_num_bulk

  subroutine dsa_begin_local(i, n) bind(c, name="dsa_begin_local")
    integer(c_int), value :: i
    integer(c_int) :: n
    acc_elem => acc_local % elements(i)
    acc_model % currentelement => acc_elem
    n = acc_elem % type % numberofnodes
  end subroutine dsa_begin_local

  subroutine dsa_gather_local(n) bind(c, name="dsa_gather_local")
    integer(c_int), value :: n
    integer, pointer :: ni(:)
    ni => acc_elem % nodeindexes
    acc_nodes % x(1:n) = acc_local % nodes % x(ni(1:n))
    acc_nodes % y(1:n) = acc_local % nodes % y(ni(1:n))
    acc_nodes % z(1:n) = acc_local % nodes % z(ni(1:n))
  end subroutine dsa_gather_local

  function dsa_bodyforce() result(r) bind(c, name="dsa_bodyforce")
    integer(c_int) :: r
    acc_bf => GetBodyForce(acc_elem)
    r = merge(1, 0, associated(acc_bf))
  end function dsa_bodyforce

  subroutine dsa_curr_dens(n, a, found) bind(c, name="dsa_curr_dens")
    integer(c_int), value :: n
    real(c_double) :: a(3, n)
    integer(c_int) :: found
    logical :: f
    call GetRealVector(acc_bf, a, 'Local Current Density', f)
    found = merge(1, 0, f)
  end subroutine dsa_curr_dens

  ! -- quadrature ----------------------------------------------------------
  function dsa_ip_begin(edge) result(n) bind(c, name="dsa_ip_begin")
    integer(c_int), value :: edge
    integer(c_int) :: n
    if (edge /= 0) then
      acc_ip = GaussPoints(acc_elem, EdgeBasis=.true., &
          PReferenceElement=acc_piola, EdgeBasisDegree=acc_degree)
    else
      acc_ip = GaussPoints(acc_elem)
    end if
    n = acc_ip % n
  end function dsa_ip_begin

  subroutine dsa_ip_point(t, u, v, w, s) bind(c, name="dsa_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip % u(t)
    v = acc_ip % v(t)
    w = acc_ip % w(t)
    s = acc_ip % s(t)
  end subroutine dsa_ip_point

  subroutine dsa_eleminfo_local(u, v, w, n, basis, ddx, detj, stat) &
      bind(c, name="dsa_eleminfo_local")
    real(c_double), value :: u, v, w
    integer(c_int), value :: n
    real(c_double) :: basis(n), ddx(n, 3), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_elem, acc_nodes, u, v, w, detj, &
        basis(1:n), ddx(1:n, 1:3))
    stat = merge(1, 0, ok)
  end subroutine dsa_eleminfo_local

  subroutine dsa_ip_eval(n, basis, curr, point, currip) bind(c, name="dsa_ip_eval")
    integer(c_int), value :: n
    real(c_double) :: basis(n), curr(3, n), point(3), currip(3)
    point(1) = sum(basis(1:n) * acc_nodes % x(1:n))
    point(2) = sum(basis(1:n) * acc_nodes % y(1:n))
    point(3) = sum(basis(1:n) * acc_nodes % z(1:n))
    currip(1) = sum(basis(1:n) * curr(1, 1:n))
    currip(2) = sum(basis(1:n) * curr(2, 1:n))
    currip(3) = sum(basis(1:n) * curr(3, 1:n))
  end subroutine dsa_ip_eval

  ! -- global-mesh search ---------------------------------------------------
  function dsa_find_leaf(px, py, pz, dim) result(r) bind(c, name="dsa_find_leaf")
    real(c_double), value :: px, py, pz
    integer(c_int), value :: dim
    integer(c_int) :: r
    real(kind=dp) :: point(3)
    nullify(acc_leaf)
    point = [px, py, pz]
    call FindLeafElements(point, dim, acc_global % rootquadrant, acc_leaf)
    r = merge(1, 0, associated(acc_leaf))
  end function dsa_find_leaf

  function dsa_leaf_count() result(n) bind(c, name="dsa_leaf_count")
    integer(c_int) :: n
    n = acc_leaf % nelemsinquadrant
  end function dsa_leaf_count

  subroutine dsa_begin_leaf(k, nb, ne) bind(c, name="dsa_begin_leaf")
    integer(c_int), value :: k
    integer(c_int) :: nb, ne
    acc_elemb => acc_global % elements(acc_leaf % elements(k))
    nb = acc_elemb % type % numberofnodes
    ne = acc_elemb % type % numberofedges
  end subroutine dsa_begin_leaf

  function dsa_global_num_bulk() result(n) bind(c, name="dsa_global_num_bulk")
    integer(c_int) :: n
    n = acc_global % numberofbulkelements
  end function dsa_global_num_bulk

  subroutine dsa_begin_global(k, nb, ne) bind(c, name="dsa_begin_global")
    integer(c_int), value :: k
    integer(c_int) :: nb, ne
    acc_elemb => acc_global % elements(k)
    nb = acc_elemb % type % numberofnodes
    ne = acc_elemb % type % numberofedges
  end subroutine dsa_begin_global

  subroutine dsa_edge_idx(ne, idx) bind(c, name="dsa_edge_idx")
    integer(c_int), value :: ne
    integer(c_int) :: idx(ne)
    idx(1:ne) = acc_elemb % edgeindexes + acc_global % numberofnodes
  end subroutine dsa_edge_idx

  subroutine dsa_gather_global(nb) bind(c, name="dsa_gather_global")
    integer(c_int), value :: nb
    integer, pointer :: ni(:)
    ni => acc_elemb % nodeindexes
    acc_nodesb % x(1:nb) = acc_global % nodes % x(ni(1:nb))
    acc_nodesb % y(1:nb) = acc_global % nodes % y(ni(1:nb))
    acc_nodesb % z(1:nb) = acc_global % nodes % z(ni(1:nb))
  end subroutine dsa_gather_global

  function dsa_point_in_elem(point, lc) result(r) bind(c, name="dsa_point_in_elem")
    real(c_double) :: point(3), lc(3)
    integer(c_int) :: r
    real(kind=dp) :: local(3)
    logical :: ok
    ok = PointInElement(acc_elemb, acc_nodesb, point, local)
    lc = local
    r = merge(1, 0, ok)
  end function dsa_point_in_elem

  subroutine dsa_info_noelem(i) bind(c, name="dsa_info_noelem")
    integer(c_int), value :: i
    call Info(caller, 'Point ' // i2s(i) // ' was not found in any of the elements!', level=15)
  end subroutine dsa_info_noelem

  ! -- assembly --------------------------------------------------------------
  subroutine dsa_rhs_ind(ne, idx, rhs) bind(c, name="dsa_rhs_ind")
    integer(c_int), value :: ne
    integer(c_int) :: idx(ne), rhs(ne)
    rhs(1:ne) = acc_rhs % perm(idx(1:ne))
  end subroutine dsa_rhs_ind

  subroutine dsa_add_rhs(ne, rhs, wb, currip, detj, s) bind(c, name="dsa_add_rhs")
    integer(c_int), value :: ne
    integer(c_int) :: rhs(ne)
    real(c_double) :: wb(ne, 3), currip(3)
    real(c_double), value :: detj, s
    if (all(rhs(1:ne) > 0)) then
      acc_rhs % values(rhs(1:ne)) = acc_rhs % values(rhs(1:ne)) + &
          matmul(wb(1:ne, :), currip) * detj * s
    end if
  end subroutine dsa_add_rhs

  function dsa_fix_assoc() result(r) bind(c, name="dsa_fix_assoc")
    integer(c_int) :: r
    r = merge(1, 0, associated(acc_fix))
  end function dsa_fix_assoc

  subroutine dsa_add_fix(nb, ddxb, currip, detj, s) bind(c, name="dsa_add_fix")
    integer(c_int), value :: nb
    real(c_double) :: ddxb(nb, 3), currip(3)
    real(c_double), value :: detj, s
    integer :: fix(nb)
    fix(1:nb) = acc_fix % perm(acc_elemb % nodeindexes)
    if (all(fix(1:nb) > 0)) then
      acc_fix % values(fix(1:nb)) = acc_fix % values(fix(1:nb)) + &
          matmul(ddxb(1:nb, :), currip) * detj * s
    end if
  end subroutine dsa_add_fix

  subroutine dsa_eleminfo_global(ub, vb, wb, nb, ne, basisb, ddxb, wb0, rotw, detjb) &
      bind(c, name="dsa_eleminfo_global")
    real(c_double), value :: ub, vb, wb
    integer(c_int), value :: nb, ne
    real(c_double) :: basisb(nb), ddxb(nb, 3), wb0(ne, 3), rotw(ne, 3), detjb
    logical :: stat
    stat = ElementInfo(acc_elemb, acc_nodesb, ub, vb, wb, detjb, &
        basisb(1:nb), ddxb(1:nb, 1:3), edgebasis=wb0, rotbasis=rotw, &
        usolver=acc_gsolver)
  end subroutine dsa_eleminfo_global

  ! -- tail --------------------------------------------------------------------
  subroutine dsa_finish(nofails, nofound) bind(c, name="dsa_finish")
    integer(c_int), value :: nofails, nofound
    call Info(caller, 'Integration finished', level=15)
    if (nofails == 0) then
      call Info(caller, 'Found all nodes ' // i2s(nofound) // &
          ' in the target mesh', level=6)
    else
      call Warn(caller, 'Points not found: ' // i2s(nofails) // ' (found ' &
          // i2s(nofound) // ')')
    end if
    print *, 'Global source sum:', sum(acc_rhs % values), sum(abs(acc_rhs % values))
    if (associated(acc_fix)) then
      print *, 'Global fixing sum:', sum(acc_fix % values), sum(abs(acc_fix % values))
    end if
  end subroutine dsa_finish

  subroutine dsa_fatal_alloc() bind(c, name="dsa_fatal_alloc")
    call Fatal(caller, 'Memory allocation error.')
  end subroutine dsa_fatal_alloc

end module distributesource_acc
