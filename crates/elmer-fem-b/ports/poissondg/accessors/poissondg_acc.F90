! T021-B PoissonDG port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/PoissonDG.F90, compiled against the built .mod
! files of the target configure (ADR-057 §3, T020 Poisson pattern).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module mirrors PoissonDG's own locals (current bulk/edge elements,
! edge parents, Nodes SAVE state per internal procedure, IP shallow
! copies) with identical lifetimes.
!
! Bit-identity notes: every order-sensitive reduction (SUM), MATMUL
! update, formatted I/O (Warn/Fatal inside FindParentUVW), and
! whole-array zero/fill stays in Fortran with the original statement
! order, so results match by construction where the call sequence
! matches. LocalJumps and LocalJumpsDiscontBC are textually identical in
! the original (names only differ); each gets its own wrapper set with
! its own SAVE Nodes, mirroring the original's separate SAVE lifetimes.
!
! LocalMatrixBoundary is dead code in the original (never called from the
! entry); it has no wrappers. EpsilonBoundary/LocalSolution are
! allocated but never read on any live path (see lib.rs).
module poissondg_acc
  use iso_c_binding
  use defutils
  use elementdescription
  use messages
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(mesh_t), pointer, save :: acc_mesh => null()
  type(element_t), pointer, save :: acc_elem => null()
  ! LeftOut scratch: the original keeps LeftOut as an automatic in each
  ! internal procedure (recomputed per call, consumed within the IP loop).
  ! One SAVE copy per variant with identical lifetime: filled by *_setup,
  ! consumed by *_point before the next fill.
  real(c_double), save :: acc_leftout_lj(3) = 0.0d0
  real(c_double), save :: acc_leftout_ld(3) = 0.0d0
  type(element_t), pointer, save :: acc_edge => null()
  type(element_t), pointer, save :: acc_left => null()
  type(element_t), pointer, save :: acc_right => null()
  type(valuelist_t), pointer, save :: acc_bf => null()
  type(valuelist_t), pointer, save :: acc_bc => null()
  logical, save :: acc_use_faces = .false.
  ! LocalMatrix context (mirrors its SAVE Nodes and IP local).
  type(nodes_t), target, save :: acc_nodes_lm
  type(gaussintegrationpoints_t), save :: acc_ip_lm
  ! LocalJumps context (mirrors its SAVE Nodes x3 and IP local).
  type(nodes_t), target, save :: acc_nodes_lj, acc_leftnodes_lj, acc_rightnodes_lj
  type(gaussintegrationpoints_t), save :: acc_ip_lj
  ! LocalJumpsDiscontBC context (mirrors its SAVE Nodes x3 and IP local).
  type(nodes_t), target, save :: acc_nodes_ld, acc_leftnodes_ld, acc_rightnodes_ld
  type(gaussintegrationpoints_t), save :: acc_ip_ld

contains

  ! -- entry state ----------------------------------------------------
  subroutine pda_begin(h_solver) bind(c, name="pda_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine pda_begin

  subroutine pda_touch_mesh() bind(c, name="pda_touch_mesh")
    acc_mesh => GetMesh()
  end subroutine pda_touch_mesh

  ! N = 2 * Mesh % MaxElementDOFs (first-mesh SAVE sizing).
  function pda_first_n() result(n) bind(c, name="pda_first_n")
    integer(c_int) :: n
    n = 2 * acc_mesh % maxelementdofs
  end function pda_first_n

  ! CoordinateSystemDimension + edge/face count (Edges for 2D, Faces else).
  subroutine pda_edge_counts(dim, np) bind(c, name="pda_edge_counts")
    integer(c_int) :: dim, np
    dim = CoordinateSystemDimension()
    if (dim == 2) then
      acc_use_faces = .false.
      np = acc_mesh % numberofedges
    else
      acc_use_faces = .true.
      np = acc_mesh % numberoffaces
    end if
  end subroutine pda_edge_counts

  function pda_nof_active() result(n) bind(c, name="pda_nof_active")
    integer(c_int) :: n
    n = GetNOFActive()
  end function pda_nof_active

  subroutine pda_default_initialize() bind(c, name="pda_default_initialize")
    call DefaultInitialize()
  end subroutine pda_default_initialize

  function pda_gamma(which) result(g) bind(c, name="pda_gamma")
    integer(c_int), value :: which
    real(c_double) :: g
    logical :: gotit
    real(c_double) :: v
    if (which == 1) then
      v = ListGetConstReal(acc_solver % values, 'gamma1', gotit)
      if (.not. gotit) v = 1.0d3
    else
      v = ListGetConstReal(acc_solver % values, 'gamma2', gotit)
      if (.not. gotit) v = 1.0d-3
    end if
    g = v
  end function pda_gamma

  function pda_beta(gamma) result(b) bind(c, name="pda_beta")
    real(c_double), value :: gamma
    real(c_double) :: b
    logical :: gotit
    real(c_double) :: v
    v = ListGetConstReal(acc_solver % values, 'beta', gotit)
    if (.not. gotit) v = gamma
    b = v
  end function pda_beta

  ! -- bulk-element context --------------------------------------------
  subroutine pda_begin_active(t, n) bind(c, name="pda_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes(acc_elem)
  end subroutine pda_begin_active

  subroutine pda_bodyforce() bind(c, name="pda_bodyforce")
    acc_bf => GetBodyForce(acc_elem)
  end subroutine pda_bodyforce

  subroutine pda_load_source(n, load, found) bind(c, name="pda_load_source")
    integer(c_int), value :: n
    real(c_double) :: load(n)
    integer(c_int) :: found
    logical :: f
    load(1:n) = GetReal(acc_bf, 'Source', f)
    found = merge(1, 0, f)
  end subroutine pda_load_source

  ! -- LocalMatrix context ----------------------------------------------
  subroutine pda_lm_nodes() bind(c, name="pda_lm_nodes")
    call GetElementNodes(acc_nodes_lm, acc_elem)
  end subroutine pda_lm_nodes

  function pda_lm_ip_begin() result(n) bind(c, name="pda_lm_ip_begin")
    integer(c_int) :: n
    acc_ip_lm = GaussPoints(acc_elem)
    n = acc_ip_lm % n
  end function pda_lm_ip_begin

  subroutine pda_lm_ip_point(t, u, v, w, s) bind(c, name="pda_lm_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip_lm % u(t)
    v = acc_ip_lm % v(t)
    w = acc_ip_lm % w(t)
    s = acc_ip_lm % s(t)
  end subroutine pda_lm_ip_point

  subroutine pda_lm_eleminfo(u, v, w, n, basis, dbasis, detj, stat) &
      bind(c, name="pda_lm_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: n
    real(c_double) :: basis(n), dbasis(n, 3), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_elem, acc_nodes_lm, u, v, w, detj, &
        basis(1:n), dbasis(1:n, 1:3))
    stat = merge(1, 0, ok)
  end subroutine pda_lm_eleminfo

  function pda_dot(n, a, b) result(s) bind(c, name="pda_dot")
    integer(c_int), value :: n
    real(c_double) :: a(n), b(n), s
    s = sum(a(1:n) * b(1:n))
  end function pda_dot

  subroutine pda_lm_assemble(n, dbasis, sdet, basis, loadatip, ld, stiff, &
      lf, force) bind(c, name="pda_lm_assemble")
    integer(c_int), value :: n, ld, lf
    real(c_double) :: dbasis(n, 3), basis(n), stiff(ld, ld), force(lf)
    real(c_double), value :: sdet, loadatip
    ! Exact transcription of the two update statements (same compiler,
    ! same evaluation order as the original internal procedure).
    stiff(1:n, 1:n) = stiff(1:n, 1:n) + sdet * &
        matmul(dbasis(1:n, 1:3), transpose(dbasis(1:n, 1:3)))
    force(1:n) = force(1:n) + sdet * loadatip * basis(1:n)
  end subroutine pda_lm_assemble

  ! -- edge/face context -------------------------------------------------
  subroutine pda_begin_edge(t, n, nl, nr, active, hasright) &
      bind(c, name="pda_begin_edge")
    integer(c_int), value :: t
    integer(c_int) :: n, nl, nr, active, hasright
    logical :: l
    if (acc_use_faces) then
      acc_edge => acc_mesh % faces(t)
    else
      acc_edge => acc_mesh % edges(t)
    end if
    l = ActiveBoundaryElement(acc_edge, DGBoundary=.TRUE.)
    active = merge(1, 0, l)
    if (active == 0) return
    acc_left => acc_edge % boundaryinfo % left
    acc_right => acc_edge % boundaryinfo % right
    hasright = merge(1, 0, associated(acc_right))
    if (hasright == 0) return
    n = GetElementNOFNodes(acc_edge)
    nl = GetElementNOFDOFs(acc_left)
    nr = GetElementNOFDOFs(acc_right)
  end subroutine pda_begin_edge

  ! GetBC + 'Discont BC' logical (the LocalJumps vs LocalJumpsDiscontBC
  ! branch). Returns 1 when the discont-BC path must be taken.
  function pda_edge_discont() result(r) bind(c, name="pda_edge_discont")
    integer(c_int) :: r
    logical :: found, l
    acc_bc => GetBC(acc_edge)
    r = 0
    if (associated(acc_bc)) then
      l = GetLogical(acc_bc, 'Discont BC', found)
      if (l) r = 1
    end if
  end function pda_edge_discont

  ! -- LocalJumps context -------------------------------------------------
  subroutine pda_lj_setup(n, nl, nr, he) bind(c, name="pda_lj_setup")
    integer(c_int), value :: n, nl, nr
    real(c_double) :: he
    call GetElementNodes(acc_nodes_lj, acc_edge)
    call GetElementNodes(acc_leftnodes_lj, acc_left)
    call GetElementNodes(acc_rightnodes_lj, acc_right)
    he = ElementDiameter(acc_edge, acc_nodes_lj)
    acc_leftout_lj(1) = sum(acc_leftnodes_lj % x(1:nl)) / nl
    acc_leftout_lj(2) = sum(acc_leftnodes_lj % y(1:nl)) / nl
    acc_leftout_lj(3) = sum(acc_leftnodes_lj % z(1:nl)) / nl
    acc_leftout_lj(1) = sum(acc_nodes_lj % x(1:n)) / n - acc_leftout_lj(1)
    acc_leftout_lj(2) = sum(acc_nodes_lj % y(1:n)) / n - acc_leftout_lj(2)
    acc_leftout_lj(3) = sum(acc_nodes_lj % z(1:n)) / n - acc_leftout_lj(3)
  end subroutine pda_lj_setup

  function pda_lj_ip_begin() result(n) bind(c, name="pda_lj_ip_begin")
    integer(c_int) :: n
    acc_ip_lj = GaussPoints(acc_edge)
    n = acc_ip_lj % n
  end function pda_lj_ip_begin

  subroutine pda_lj_ip_point(t, u, v, w, s) bind(c, name="pda_lj_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip_lj % u(t)
    v = acc_ip_lj % v(t)
    w = acc_ip_lj % w(t)
    s = acc_ip_lj % s(t)
  end subroutine pda_lj_ip_point

  ! One integration point: edge ElementInfo, normal (+flip),
  ! FindParentUVW + ElementInfo for each parent. U/V/W/S are in/out,
  ! exactly as in the original loop body.
  subroutine pda_lj_point(u, v, w, s, n, nl, nr, basis, dbasis, &
      leftbasis, leftdbasis, rightbasis, rightdbasis, normal, detj, stat) &
      bind(c, name="pda_lj_point")
    real(c_double) :: u, v, w, s
    integer(c_int), value :: n, nl, nr
    real(c_double) :: basis(n), dbasis(n, 3)
    real(c_double) :: leftbasis(nl), leftdbasis(nl, 3)
    real(c_double) :: rightbasis(nr), rightdbasis(nr, 3)
    real(c_double) :: normal(3), detj
    integer(c_int) :: stat
    logical :: ok
    real(c_double) :: nor(3)
    ok = ElementInfo(acc_edge, acc_nodes_lj, u, v, w, detj, &
        basis(1:n), dbasis(1:n, 1:3))
    s = s * detj
    nor = NormalVector(acc_edge, acc_nodes_lj, u, v, .FALSE.)
    if (sum(acc_leftout_lj * nor) < 0) nor = -nor
    normal = nor
    call FindParentUVW(acc_edge, n, acc_left, nl, u, v, w, basis)
    ok = ElementInfo(acc_left, acc_leftnodes_lj, u, v, w, detj, &
        leftbasis(1:nl), leftdbasis(1:nl, 1:3))
    call FindParentUVW(acc_edge, n, acc_right, nr, u, v, w, basis)
    ok = ElementInfo(acc_right, acc_rightnodes_lj, u, v, w, detj, &
        rightbasis(1:nr), rightdbasis(1:nr, 1:3))
    stat = merge(1, 0, ok)
  end subroutine pda_lj_point

  ! Jump/AverageFlux + (nl+nr)^2 STIFF update, transcribed exactly.
  subroutine pda_lj_assemble(nl, nr, gamma, he, leftbasis, leftdbasis, &
      rightbasis, rightdbasis, normal, s, ld, stiff) &
      bind(c, name="pda_lj_assemble")
    integer(c_int), value :: nl, nr, ld
    real(c_double), value :: gamma, he, s
    real(c_double) :: leftbasis(nl), leftdbasis(nl, 3)
    real(c_double) :: rightbasis(nr), rightdbasis(nr, 3)
    real(c_double) :: normal(3), stiff(ld, ld)
    integer :: i, p, q, m
    real(c_double) :: jump(nl+nr), averageflux(nl+nr)
    real(c_double) :: leftdn(nl), rightdn(nr)
    m = nl + nr
    jump(1:nl) = -leftbasis(1:nl)
    jump(nl+1:m) = rightbasis(1:nr)
    do i = 1, nl
      leftdn(i) = sum(leftdbasis(i, :) * normal(:))
    end do
    do i = 1, nr
      rightdn(i) = sum(rightdbasis(i, :) * normal(:))
    end do
    averageflux(1:nl) = leftdn(1:nl) / 2.0d0
    averageflux(nl+1:m) = rightdn(1:nr) / 2.0d0
    do p = 1, m
      do q = 1, m
        stiff(p, q) = stiff(p, q) + (gamma / he) * jump(p) * jump(q) * s
        stiff(p, q) = stiff(p, q) - averageflux(p) * jump(q) * s
        stiff(p, q) = stiff(p, q) + jump(p) * averageflux(q) * s
      end do
    end do
  end subroutine pda_lj_assemble

  ! -- LocalJumpsDiscontBC context (identical body, own SAVE state) --------
  subroutine pda_ld_setup(n, nl, nr, he) bind(c, name="pda_ld_setup")
    integer(c_int), value :: n, nl, nr
    real(c_double) :: he
    call GetElementNodes(acc_nodes_ld, acc_edge)
    call GetElementNodes(acc_leftnodes_ld, acc_left)
    call GetElementNodes(acc_rightnodes_ld, acc_right)
    he = ElementDiameter(acc_edge, acc_nodes_ld)
    acc_leftout_ld(1) = sum(acc_leftnodes_ld % x(1:nl)) / nl
    acc_leftout_ld(2) = sum(acc_leftnodes_ld % y(1:nl)) / nl
    acc_leftout_ld(3) = sum(acc_leftnodes_ld % z(1:nl)) / nl
    acc_leftout_ld(1) = sum(acc_nodes_ld % x(1:n)) / n - acc_leftout_ld(1)
    acc_leftout_ld(2) = sum(acc_nodes_ld % y(1:n)) / n - acc_leftout_ld(2)
    acc_leftout_ld(3) = sum(acc_nodes_ld % z(1:n)) / n - acc_leftout_ld(3)
  end subroutine pda_ld_setup

  function pda_ld_ip_begin() result(n) bind(c, name="pda_ld_ip_begin")
    integer(c_int) :: n
    acc_ip_ld = GaussPoints(acc_edge)
    n = acc_ip_ld % n
  end function pda_ld_ip_begin

  subroutine pda_ld_ip_point(t, u, v, w, s) bind(c, name="pda_ld_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip_ld % u(t)
    v = acc_ip_ld % v(t)
    w = acc_ip_ld % w(t)
    s = acc_ip_ld % s(t)
  end subroutine pda_ld_ip_point

  subroutine pda_ld_point(u, v, w, s, n, nl, nr, basis, dbasis, &
      leftbasis, leftdbasis, rightbasis, rightdbasis, normal, detj, stat) &
      bind(c, name="pda_ld_point")
    real(c_double) :: u, v, w, s
    integer(c_int), value :: n, nl, nr
    real(c_double) :: basis(n), dbasis(n, 3)
    real(c_double) :: leftbasis(nl), leftdbasis(nl, 3)
    real(c_double) :: rightbasis(nr), rightdbasis(nr, 3)
    real(c_double) :: normal(3), detj
    integer(c_int) :: stat
    logical :: ok
    real(c_double) :: nor(3)
    ok = ElementInfo(acc_edge, acc_nodes_ld, u, v, w, detj, &
        basis(1:n), dbasis(1:n, 1:3))
    s = s * detj
    nor = NormalVector(acc_edge, acc_nodes_ld, u, v, .FALSE.)
    if (sum(acc_leftout_ld * nor) < 0) nor = -nor
    normal = nor
    call FindParentUVW(acc_edge, n, acc_left, nl, u, v, w, basis)
    ok = ElementInfo(acc_left, acc_leftnodes_ld, u, v, w, detj, &
        leftbasis(1:nl), leftdbasis(1:nl, 1:3))
    call FindParentUVW(acc_edge, n, acc_right, nr, u, v, w, basis)
    ok = ElementInfo(acc_right, acc_rightnodes_ld, u, v, w, detj, &
        rightbasis(1:nr), rightdbasis(1:nr, 1:3))
    stat = merge(1, 0, ok)
  end subroutine pda_ld_point

  subroutine pda_ld_assemble(nl, nr, gamma, he, leftbasis, leftdbasis, &
      rightbasis, rightdbasis, normal, s, ld, stiff) &
      bind(c, name="pda_ld_assemble")
    integer(c_int), value :: nl, nr, ld
    real(c_double), value :: gamma, he, s
    real(c_double) :: leftbasis(nl), leftdbasis(nl, 3)
    real(c_double) :: rightbasis(nr), rightdbasis(nr, 3)
    real(c_double) :: normal(3), stiff(ld, ld)
    integer :: i, p, q, m
    real(c_double) :: jump(nl+nr), averageflux(nl+nr)
    real(c_double) :: leftdn(nl), rightdn(nr)
    m = nl + nr
    jump(1:nl) = -leftbasis(1:nl)
    jump(nl+1:m) = rightbasis(1:nr)
    do i = 1, nl
      leftdn(i) = sum(leftdbasis(i, :) * normal(:))
    end do
    do i = 1, nr
      rightdn(i) = sum(rightdbasis(i, :) * normal(:))
    end do
    averageflux(1:nl) = leftdn(1:nl) / 2.0d0
    averageflux(nl+1:m) = rightdn(1:nr) / 2.0d0
    do p = 1, m
      do q = 1, m
        stiff(p, q) = stiff(p, q) + (gamma / he) * jump(p) * jump(q) * s
        stiff(p, q) = stiff(p, q) - averageflux(p) * jump(q) * s
        stiff(p, q) = stiff(p, q) + jump(p) * averageflux(q) * s
      end do
    end do
  end subroutine pda_ld_assemble

  ! -- shared assembly tail ----------------------------------------------
  subroutine pda_default_update(ld, stiff, lf, force) &
      bind(c, name="pda_default_update")
    integer(c_int), value :: ld, lf
    real(c_double) :: stiff(ld, ld), force(lf)
    call DefaultUpdateEquations(stiff(1:ld, 1:ld), force(1:lf))
  end subroutine pda_default_update

  subroutine pda_default_update_edge(ld, stiff, lf, force) &
      bind(c, name="pda_default_update_edge")
    integer(c_int), value :: ld, lf
    real(c_double) :: stiff(ld, ld), force(lf)
    call DefaultUpdateEquations(stiff(1:ld, 1:ld), force(1:lf), acc_edge)
  end subroutine pda_default_update_edge

  subroutine pda_default_finish() bind(c, name="pda_default_finish")
    call DefaultFinishAssembly()
  end subroutine pda_default_finish

  subroutine pda_default_dirichlet() bind(c, name="pda_default_dirichlet")
    call DefaultDirichletBCs()
  end subroutine pda_default_dirichlet

  function pda_default_solve() result(norm) bind(c, name="pda_default_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function pda_default_solve

  subroutine pda_fatal_alloc() bind(c, name="pda_fatal_alloc")
    call Fatal('PoissonDG', 'Memory allocation error.')
  end subroutine pda_fatal_alloc

end module poissondg_acc
