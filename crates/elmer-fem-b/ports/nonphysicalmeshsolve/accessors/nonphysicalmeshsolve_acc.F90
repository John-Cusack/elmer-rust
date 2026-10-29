! T021-B NonphysicalMeshSolve port: statement-level BIND(C) wrappers for
! the Elmer calls in fem/src/modules/NonphysicalMeshSolve.F90, compiled
! against the built .mod files of the target configure (A2 shipped
! mechanism, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module mirrors MeshSolver's own SAVE state (the entry SAVE arrays that
! must alias the mesh node storage, VisitedTimes/dim, the target-field
! handles, the per-procedure SAVE Nodes/IP states, the ParentBasis SAVE
! of MeshBoundary, the NodalDisplacementPenalty Visited flag) with
! identical lifetimes.
!
! Bit-identity rules followed throughout:
! - Whole-array fills, SUM reductions, the (p,q,i,j) element assembly
!   nest, the symmetric fill, the bubble-DOF clamp, NormalVector,
!   GetParentUVW, the (TrueX/OrigX) mesh-pointer swaps, the MovingMesh
!   coordinate update, Info output, and the CRS matrix penalty update
!   stay in Fortran wrappers with the exact statement order of the
!   original, so bit-identity holds by construction wherever the call
!   sequence matches.
! - The Rust side owns the numeric SAVE buffers (STIFF, FORCE, LOAD,
!   ElasticModulus, PoissonRatio, Alpha, Beta, Gamma, RefSurface) and
!   passes (ld,ld) / exact-shape dummies so no (nd,nd)-over-(ld,ld)
!   stride mismatch can occur; wrappers update (1:ld,1:ld) sections
!   exactly as the original updates whole arrays.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module nonphysicalmeshsolve_acc
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
  type(element_t), pointer, save :: acc_elem => null()
  type(valuelist_t), pointer, save :: acc_mat => null()
  type(valuelist_t), pointer, save :: acc_bc => null()
  ! Frame state re-pointed on every entry (mirrors the MeshSolver
  ! locals MeshSol/MeshPerm/MeshUpdate, plus dim and the node pointers).
  type(variable_t), pointer, save :: acc_meshsol => null()
  integer, pointer, save :: acc_meshperm(:) => null()
  real(c_double), pointer, save :: acc_meshupdate(:) => null()
  type(variable_t), pointer, save :: acc_targetsol => null()
  type(variable_t), pointer, save :: acc_gradsol => null()
  ! Entry SAVE mirrors: pointer swaps (CONTIG spelling matches the unit,
  ! compiled with -DCONTIG= like the Elmer tree itself).
  real(kind=dp), pointer CONTIG, save :: acc_truex(:) => null()
  real(kind=dp), pointer CONTIG, save :: acc_truey(:) => null()
  real(kind=dp), pointer CONTIG, save :: acc_truez(:) => null()
  real(kind=dp), pointer CONTIG, save :: acc_origx(:) => null()
  real(kind=dp), pointer CONTIG, save :: acc_origy(:) => null()
  real(kind=dp), pointer CONTIG, save :: acc_origz(:) => null()
  real(c_double), allocatable, save :: acc_prev(:)
  integer, save :: acc_visited = 0
  integer, save :: acc_dim = 0
  logical, save :: acc_cumulative = .false.
  logical, save :: acc_moving = .false.
  logical, save :: acc_got_field = .false.
  logical, save :: acc_got_surface = .false.
  logical, save :: acc_got_grad = .false.
  ! LocalMatrix context (mirrors its SAVE Nodes; IP reassigned per call
  ! exactly as the original's local is).
  type(nodes_t), target, save :: acc_nodes_lm
  type(gaussintegrationpoints_t), save :: acc_ip_lm
  ! MeshBoundary context (mirrors its SAVE Nodes/ParentNodes, the
  ! AllocationsDone flag, and the ParentBasis SAVE pair).
  type(nodes_t), target, save :: acc_nodes_mb
  type(nodes_t), target, save :: acc_nodes_parent
  type(gaussintegrationpoints_t), save :: acc_ip_mb
  logical, save :: acc_pb_alloc = .false.
  integer, save :: acc_np_save = 0
  real(c_double), allocatable, save :: acc_parentbasis(:)
  real(c_double), allocatable, save :: acc_dparentbasis(:, :)
  ! NodalDisplacementPenalty context (mirrors its SAVE Visited flag).
  logical, save :: acc_pen_visited = .false.

contains

  ! -- entry state ----------------------------------------------------
  subroutine nma_begin(h_model, h_solver) bind(c, name="nma_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine nma_begin

  subroutine nma_info_start() bind(c, name="nma_info_start")
    call info('NonPhysicalMeshSolve', '-------------------------------------', level=4)
    call info('NonPhysicalMeshSolve', 'Nonphysical Mesh Solver:', level=4)
    call info('NonPhysicalMeshSolve', '-------------------------------------', level=4)
  end subroutine nma_info_start

  function nma_assoc_matrix() result(r) bind(c, name="nma_assoc_matrix")
    integer(c_int) :: r
    r = merge(1, 0, associated(acc_solver % matrix))
  end function nma_assoc_matrix

  ! -- frame: pointer capture, dim, counts ------------------------------
  subroutine nma_frame() bind(c, name="nma_frame")
    acc_meshsol => acc_solver % variable
    acc_meshperm => acc_meshsol % perm
    acc_meshupdate => acc_meshsol % values
    acc_dim = CoordinateSystemDimension()
  end subroutine nma_frame

  function nma_stdofs() result(n) bind(c, name="nma_stdofs")
    integer(c_int) :: n
    n = acc_meshsol % dofs
  end function nma_stdofs

  function nma_dim() result(n) bind(c, name="nma_dim")
    integer(c_int) :: n
    n = acc_dim
  end function nma_dim

  function nma_local_nodes() result(n) bind(c, name="nma_local_nodes")
    integer(c_int) :: n
    n = count(acc_meshperm > 0)
  end function nma_local_nodes

  function nma_max_elem_dofs() result(n) bind(c, name="nma_max_elem_dofs")
    integer(c_int) :: n
    n = acc_solver % mesh % maxelementdofs
  end function nma_max_elem_dofs

  function nma_cumulative() result(r) bind(c, name="nma_cumulative")
    integer(c_int) :: r
    type(valuelist_t), pointer :: params
    logical :: found
    params => GetSolverParams()
    acc_cumulative = ListGetLogical(params, 'Cumulative Displacement', found)
    r = merge(1, 0, acc_cumulative)
  end function nma_cumulative

  ! -- first-call permanent storage (Orig node coordinates) --------------
  subroutine nma_first_alloc() bind(c, name="nma_first_alloc")
    integer :: n, istat
    if (.not. acc_cumulative) then
      n = size(acc_solver % mesh % nodes % x)
      allocate(acc_origx(n), acc_origy(n), acc_origz(n), stat=istat)
      if (istat /= 0) then
        call fatal('NonPhysicalMeshSolve', 'Memory allocation error.')
      end if
      acc_origx = acc_solver % mesh % nodes % x
      acc_origy = acc_solver % mesh % nodes % y
      acc_origz = acc_solver % mesh % nodes % z
    end if
  end subroutine nma_first_alloc

  subroutine nma_visit() bind(c, name="nma_visit")
    acc_visited = acc_visited + 1
  end subroutine nma_visit

  subroutine nma_nodes_capture() bind(c, name="nma_nodes_capture")
    acc_truex => acc_solver % mesh % nodes % x
    acc_truey => acc_solver % mesh % nodes % y
    acc_truez => acc_solver % mesh % nodes % z
  end subroutine nma_nodes_capture

  subroutine nma_nodes_to_orig() bind(c, name="nma_nodes_to_orig")
    if (.not. acc_cumulative) then
      acc_solver % mesh % nodes % x => acc_origx
      acc_solver % mesh % nodes % y => acc_origy
      acc_solver % mesh % nodes % z => acc_origz
    end if
  end subroutine nma_nodes_to_orig

  subroutine nma_nodes_restore() bind(c, name="nma_nodes_restore")
    if (.not. acc_cumulative) then
      acc_solver % mesh % nodes % x => acc_truex
      acc_solver % mesh % nodes % y => acc_truey
      acc_solver % mesh % nodes % z => acc_truez
    end if
  end subroutine nma_nodes_restore

  ! -- MovingMesh block (VisitedTimes-driven PrevMeshUpdate) --------------
  subroutine nma_moving_update() bind(c, name="nma_moving_update")
    type(valuelist_t), pointer :: params
    logical :: found, moving
    integer :: n
    params => GetSolverParams()
    moving = ListGetLogical(params, 'Moving Mesh', found)
    if (.not. found) moving = .true.
    moving = moving .and. (acc_visited > 1)
    acc_moving = moving
    if (acc_moving) then
      if (acc_visited == 2) then
        n = size(acc_meshupdate)
        allocate(acc_prev(n))
        acc_prev = 0.0_c_double
      end if
      acc_prev = acc_meshupdate
    end if
  end subroutine nma_moving_update

  ! -- permanent storage failure -----------------------------------------
  subroutine nma_fatal_alloc() bind(c, name="nma_fatal_alloc")
    call fatal('NonPhysicalMeshSolve', 'Memory allocation error.')
  end subroutine nma_fatal_alloc

  ! -- assembly driver -----------------------------------------------------
  function nma_assembly_begin() result(noactive) bind(c, name="nma_assembly_begin")
    integer(c_int) :: noactive
    real(c_double) :: at, at0
    at = CPUTime()
    at0 = RealTime()
    call DefaultInitialize()
    call StartAdvanceOutput('NonphysicalMeshSolve', 'Assembly: ')
    noactive = GetNOFActive()
  end function nma_assembly_begin

  subroutine nma_advance(t, noactive) bind(c, name="nma_advance")
    integer(c_int), value :: t, noactive
    call AdvanceOutput(t, noactive)
  end subroutine nma_advance

  ! -- bulk-element context --------------------------------------------------
  subroutine nma_begin_active(t, n, nd, nb) bind(c, name="nma_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, nb
    acc_elem => GetActiveElement(t)
    nd = GetElementNOFDOFs()
    nb = GetElementNOFBDOFs()
    n = GetElementNOFNodes()
  end subroutine nma_begin_active

  subroutine nma_mat_props(n, buf_n, emod, pr) bind(c, name="nma_mat_props")
    integer(c_int), value :: n, buf_n
    real(c_double) :: emod(buf_n), pr(buf_n)
    logical :: found
    acc_mat => GetMaterial()
    emod(1:n) = GetReal(acc_mat, 'Mesh Elastic Modulus', found)
    if (.not. found) emod(1:n) = GetReal(acc_mat, 'Youngs Modulus', found)
    if (.not. found) emod(1:n) = 1.0d0
    pr(1:n) = GetReal(acc_mat, 'Mesh Poisson Ratio', found)
    if (.not. found) pr(1:n) = GetReal(acc_mat, 'Poisson Ratio', found)
    if (.not. found) pr(1:n) = 0.25d0
  end subroutine nma_mat_props

  ! -- LocalMatrix (nodal Lame + per-IP assembly + finish) ----------------------
  subroutine nma_lm_nodes() bind(c, name="nma_lm_nodes")
    call GetElementNodes(acc_nodes_lm)
  end subroutine nma_lm_nodes

  subroutine nma_lame(n, buf_n, emod, pr, lame1, lame2, plane) bind(c, name="nma_lame")
    integer(c_int), value :: n, buf_n, plane
    real(c_double) :: emod(buf_n), pr(buf_n), lame1(n), lame2(n)
    if (plane /= 0) then
      lame1(1:n) = emod(1:n) * pr(1:n) / ((1.0d0 - pr(1:n)**2))
    else
      lame1(1:n) = emod(1:n) * pr(1:n) / &
          ((1.0d0 + pr(1:n)) * (1.0d0 - 2.0d0 * pr(1:n)))
    end if
    lame2(1:n) = emod(1:n) / (2 * (1.0d0 + pr(1:n)))
  end subroutine nma_lame

  function nma_lm_coeff() result(c) bind(c, name="nma_lm_coeff")
    real(c_double) :: c
    type(valuelist_t), pointer :: params
    logical :: stat
    params => GetSolverParams()
    c = ListGetConstReal(params, 'Mass Coefficient', stat)
  end function nma_lm_coeff

  function nma_lm_ip_begin() result(n) bind(c, name="nma_lm_ip_begin")
    integer(c_int) :: n
    acc_ip_lm = GaussPoints(acc_elem)
    n = acc_ip_lm % n
  end function nma_lm_ip_begin

  subroutine nma_lm_ip_point(t, u, v, w, s) bind(c, name="nma_lm_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip_lm % u(t)
    v = acc_ip_lm % v(t)
    w = acc_ip_lm % w(t)
    s = acc_ip_lm % s(t)
  end subroutine nma_lm_ip_point

  subroutine nma_lm_eleminfo(u, v, w, n, nd, basis, dbasis, detj, stat) &
      bind(c, name="nma_lm_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: n, nd
    real(c_double) :: basis(nd), dbasis(nd, 3), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_elem, acc_nodes_lm, u, v, w, detj, &
        basis(1:nd), dbasis(1:nd, 1:3))
    stat = merge(1, 0, ok)
  end subroutine nma_lm_eleminfo

  function nma_dot(n, a, b) result(s) bind(c, name="nma_dot")
    integer(c_int), value :: n
    real(c_double) :: a(n), b(n), s
    s = sum(a(1:n) * b(1:n))
  end function nma_dot

  subroutine nma_lm_assemble(ld, sdim, nd, s, l1, l2, coeff, dbasis, stiff) &
      bind(c, name="nma_lm_assemble")
    integer(c_int), value :: ld, sdim, nd
    real(c_double), value :: s, l1, l2, coeff
    real(c_double) :: dbasis(nd, 3)
    real(c_double), target :: stiff(ld, ld)
    real(c_double), pointer :: a(:, :)
    integer :: i, j, p, q
    ! Exact transcription of the (p,q,i,j) update nest.
    do p = 1, nd
      do q = p, nd
        a => stiff(sdim * (p - 1) + 1:sdim * p, sdim * (q - 1) + 1:sdim * q)
        do i = 1, sdim
          do j = 1, sdim
            a(i, j) = a(i, j) + s * l1 * dbasis(q, j) * dbasis(p, i)
            a(i, i) = a(i, i) + s * l2 * dbasis(q, j) * dbasis(p, j)
            a(i, j) = a(i, j) + s * l2 * dbasis(q, i) * dbasis(p, j)
          end do
          a(i, i) = a(i, i) + s * coeff * dbasis(q, i) * dbasis(p, i)
        end do
      end do
    end do
  end subroutine nma_lm_assemble

  ! Symmetric fill + bubble-DOF clamp (shared by LocalMatrix/MeshBoundary).
  subroutine nma_bubble(ld, sdim, nd, n, nb, stiff, force) bind(c, name="nma_bubble")
    integer(c_int), value :: ld, sdim, nd, n, nb
    real(c_double) :: stiff(ld, ld), force(ld)
    integer :: i, j, p, q
    do p = 1, sdim * nd
      do q = 1, p - 1
        stiff(p, q) = stiff(q, p)
      end do
    end do
    if (nb == 0) then
      do p = nd - acc_elem % bdofs + 1, nd
        do i = 1, sdim
          j = (p - 1) * sdim + i
          stiff(j, :) = 0.0d0
          stiff(:, j) = 0.0d0
          stiff(j, j) = 1.0d0
          force(j) = 0.0d0
        end do
      end do
    end if
  end subroutine nma_bubble

  ! -- MeshBoundary parent setup (once per call, before the IP loop;
  ! mirrors the early-RETURN block: 0 means "no parent, skip element").
  function nma_mb_parent_setup() result(ok) bind(c, name="nma_mb_parent_setup")
    integer(c_int) :: ok
    type(element_t), pointer :: parent
    integer :: np
    ok = 1
    if (acc_got_surface .and. .not. acc_got_grad) then
      if (.not. acc_pb_alloc) then
        np = acc_solver % mesh % maxelementnodes
        allocate(acc_parentbasis(np), acc_dparentbasis(np, 3))
        acc_pb_alloc = .true.
      end if
      parent => acc_elem % boundaryinfo % left
      if (.not. associated(parent)) &
          parent => acc_elem % boundaryinfo % right
      if (.not. associated(parent)) then
        ok = 0
        return
      end if
      np = GetElementNOFDOFs(parent)
      call GetElementNodes(acc_nodes_parent, parent)
      acc_np_save = np
    end if
  end function nma_mb_parent_setup

  ! -- shared assembly tail ------------------------------------------------------

  subroutine nma_default_update(ld, g, f) bind(c, name="nma_default_update")
    integer(c_int), value :: ld
    real(c_double) :: g(ld, ld), f(ld)
    call DefaultUpdateEquations(g(1:ld, 1:ld), f(1:ld))
  end subroutine nma_default_update

  subroutine nma_finish_bulk() bind(c, name="nma_finish_bulk")
    call DefaultFinishBulkAssembly()
  end subroutine nma_finish_bulk

  ! -- target field/surface setup --------------------------------------------------
  subroutine nma_target_setup(got_field, got_surface, got_grad) bind(c, name="nma_target_setup")
    integer(c_int) :: got_field, got_surface, got_grad
    type(valuelist_t), pointer :: params
    character(len=MAX_NAME_LEN) :: acc_name
    logical :: found
    params => GetSolverParams()
    acc_name = GetString(params, 'Target Field', acc_got_field)
    if (acc_got_field) then
      acc_targetsol => VariableGet(acc_model % variables, acc_name)
      if (.not. associated(acc_targetsol)) then
        call fatal('NonphysicalMeshSolve', &
            'Given > Target Field < does not exist: '//trim(acc_name))
      end if
    end if
    acc_name = GetString(params, 'Target Surface', acc_got_surface)
    acc_got_grad = .false.
    if (acc_got_surface) then
      if (acc_got_field) then
        call fatal('NonphysicalMeshSolve', &
            'Cannot have > Target Field < and > Target Surface < at same time!')
      end if
      acc_targetsol => VariableGet(acc_model % variables, acc_name)
      if (.not. associated(acc_targetsol)) then
        call fatal('NonphysicalMeshSolve', &
            'Given > Target Surface < does not exist: '//trim(acc_name))
      end if
      acc_name = GetString(params, 'Grad Surface', acc_got_grad)
      if (acc_got_grad) then
        acc_gradsol => VariableGet(acc_model % variables, acc_name)
        if (.not. associated(acc_targetsol)) then
          call fatal('NonphysicalMeshSolve', &
              'Given > Grad Surface < does not exist: '//trim(acc_name))
        end if
      end if
    end if
    got_field = merge(1, 0, acc_got_field)
    got_surface = merge(1, 0, acc_got_surface)
    got_grad = merge(1, 0, acc_got_grad)
  end subroutine nma_target_setup

  ! -- boundary loop -------------------------------------------------------------
  function nma_nbnd() result(n) bind(c, name="nma_nbnd")
    integer(c_int) :: n
    n = acc_solver % mesh % numberofboundaryelements
  end function nma_nbnd

  subroutine nma_begin_boundary(t, n, nd, nb, active) bind(c, name="nma_begin_boundary")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, nb, active
    acc_elem => GetBoundaryElement(t)
    active = merge(1, 0, ActiveBoundaryElement())
    nd = GetElementNOFDOFs()
    n = GetElementNOFNodes()
    nb = GetElementNOFBDOFs()
  end subroutine nma_begin_boundary

  function nma_bc() result(r) bind(c, name="nma_bc")
    integer(c_int) :: r
    acc_bc => GetBC()
    r = merge(1, 0, associated(acc_bc))
  end function nma_bc

  function nma_bc_vals(n, buf_n, load, alpha, beta, gamma, refsurf) result(got) &
      bind(c, name="nma_bc_vals")
    integer(c_int), value :: n, buf_n
    real(c_double) :: load(4, buf_n), alpha(3, buf_n), beta(buf_n)
    real(c_double) :: gamma(buf_n), refsurf(buf_n)
    integer(c_int) :: got
    logical :: found
    got = 0
    alpha(1, 1:n) = GetReal(acc_bc, 'Mesh Coefficient 1', found)
    if (found) got = 1
    alpha(2, 1:n) = GetReal(acc_bc, 'Mesh Coefficient 2', found)
    if (found) got = 1
    alpha(3, 1:n) = GetReal(acc_bc, 'Mesh Coefficient 3', found)
    if (found) got = 1
    load(1, 1:n) = GetReal(acc_bc, 'Mesh Force 1', found)
    if (found) got = 1
    load(2, 1:n) = GetReal(acc_bc, 'Mesh Force 2', found)
    if (found) got = 1
    load(3, 1:n) = GetReal(acc_bc, 'Mesh Force 3', found)
    if (found) got = 1
    beta(1:n) = GetReal(acc_bc, 'Mesh Normal Force', found)
    if (found) got = 1
    gamma(1:n) = GetReal(acc_bc, 'Mesh Penalty Factor', found)
    if (found) got = 1
    if (acc_got_surface) then
      refsurf(1:n) = GetReal(acc_bc, 'Reference Surface', found)
    end if
  end function nma_bc_vals

  ! -- MeshBoundary (per-IP assembly; Nodes/ParentNodes SAVE + IP mirror) -----------
  subroutine nma_mb_nodes() bind(c, name="nma_mb_nodes")
    call GetElementNodes(acc_nodes_mb)
  end subroutine nma_mb_nodes

  function nma_mb_sdim() result(n) bind(c, name="nma_mb_sdim")
    integer(c_int) :: n
    n = acc_elem % type % dimension + 1
  end function nma_mb_sdim

  function nma_mb_ip_begin() result(n) bind(c, name="nma_mb_ip_begin")
    integer(c_int) :: n
    acc_ip_mb = GaussPoints(acc_elem)
    n = acc_ip_mb % n
  end function nma_mb_ip_begin

  subroutine nma_mb_ip_point(t, u, v, w, s) bind(c, name="nma_mb_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip_mb % u(t)
    v = acc_ip_mb % v(t)
    w = acc_ip_mb % w(t)
    s = acc_ip_mb % s(t)
  end subroutine nma_mb_ip_point

  subroutine nma_mb_eleminfo(u, v, w, n, nd, basis, dbasis, detj, stat) &
      bind(c, name="nma_mb_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: n, nd
    real(c_double) :: basis(nd), dbasis(nd, 3), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_elem, acc_nodes_mb, u, v, w, detj, &
        basis(1:nd), dbasis(1:nd, 1:3))
    stat = merge(1, 0, ok)
  end subroutine nma_mb_eleminfo

  subroutine nma_mb_assemble(ld, n, nd, nb, u, v, w, s, load, alpha, beta, &
      gamma, refsurf, basis, stiff, force) bind(c, name="nma_mb_assemble")
    integer(c_int), value :: ld, n, nd, nb
    real(c_double), value :: u, v, w, s
    real(c_double) :: load(4, n), alpha(3, n), beta(n), gamma(n)
    real(c_double) :: refsurf(n), basis(nd), stiff(ld, ld), force(ld)
    real(c_double) :: loadatip(3), avec(3), normal(3), extdisp(3)
    real(c_double) :: refsurface, surface, gradsurface(3), absgrad2, g
    real(c_double) :: parentdetj
    type(element_t), pointer :: parent
    integer :: i, j, p, q, sdim, np
    logical :: ok
    sdim = acc_elem % type % dimension + 1
    loadatip = 0.0d0
    do i = 1, sdim
      loadatip(i) = sum(load(i, 1:n) * basis(1:n))
      avec(i) = sum(alpha(i, 1:n) * basis(1:n))
    end do
    normal = NormalVector(acc_elem, acc_nodes_mb, u, v, .true.)
    loadatip = loadatip + sum(beta(1:n) * basis(1:n)) * normal
    if (acc_got_surface .or. acc_got_field) then
      if (acc_got_field) then
        do i = 1, sdim
          extdisp(i) = sum(basis(1:n) * acc_targetsol % values( &
              sdim * (acc_targetsol % perm(acc_elem % nodeindexes(1:n)) - 1) + i))
        end do
      else
        refsurface = sum(basis(1:n) * refsurf(1:n))
        surface = sum(basis(1:n) * acc_targetsol % values(acc_elem % nodeindexes(1:n)))
        if (acc_got_grad) then
          do i = 1, sdim
            gradsurface(i) = &
                sum(basis(1:n) * acc_gradsol % values(sdim * &
                (acc_elem % nodeindexes(1:n) - 1) + i))
          end do
        else
          parent => acc_elem % boundaryinfo % left
          if (.not. associated(parent)) &
              parent => acc_elem % boundaryinfo % right
          np = acc_np_save
          call GetParentUVW(acc_elem, n, parent, np, u, v, w, basis)
          ok = ElementInfo(parent, acc_nodes_parent, u, v, w, parentdetj, &
              acc_parentbasis, acc_dparentbasis)
          do i = 1, sdim
            gradsurface(i) = &
                sum(acc_dparentbasis(1:np, i) * &
                acc_targetsol % values(parent % nodeindexes(1:np)))
          end do
        end if
        absgrad2 = sum(gradsurface(1:sdim)**2)
        do i = 1, sdim
          extdisp(i) = (refsurface - surface) * gradsurface(i) / absgrad2
        end do
      end if
      g = sum(gamma(1:n) * basis(1:n))
      loadatip = loadatip + extdisp * g
    else
      g = 0.0d0
    end if
    do p = 1, nd
      do q = 1, nd
        do i = 1, sdim
          stiff((p - 1) * sdim + i, (q - 1) * sdim + i) = &
              stiff((p - 1) * sdim + i, (q - 1) * sdim + i) + &
              s * (avec(i) + g) * basis(q) * basis(p)
        end do
      end do
    end do
    do q = 1, nd
      do i = 1, sdim
        force((q - 1) * sdim + i) = force((q - 1) * sdim + i) + &
            s * basis(q) * loadatip(i)
      end do
    end do
    if (nb == 0) then
      do p = nd - acc_elem % bdofs + 1, nd
        do i = 1, sdim
          j = (p - 1) * sdim + i
          stiff(j, :) = 0.0d0
          stiff(:, j) = 0.0d0
          stiff(j, j) = 1.0d0
          force(j) = 0.0d0
        end do
      end do
    end if
  end subroutine nma_mb_assemble

  ! -- solve + tail ------------------------------------------------------------------
  subroutine nma_finish_assembly() bind(c, name="nma_finish_assembly")
    call DefaultFinishAssembly()
    call info('NonPhysicalMeshSolve', 'Assembly done', level=4)
  end subroutine nma_finish_assembly

  subroutine nma_nodal_penalty() bind(c, name="nma_nodal_penalty")
    integer :: i, j, j2, k, ind
    real(c_double) :: extdisp(3), gradsurface(3), absgrad2, surface, &
        refsurface, targetcoeff, coeff
    logical :: localpenalty, gotweight, found
    integer, pointer :: penaltyperm(:), weightperm(:)
    type(variable_t), pointer :: weightsol
    type(valuelist_t), pointer :: params
    if (.not. (acc_got_surface .or. acc_got_field)) return
    params => GetSolverParams()
    targetcoeff = GetCReal(params, 'Nodal Penalty Factor', found)
    if (abs(targetcoeff) < tiny(targetcoeff)) return
    gotweight = .false.
    if (GetLogical(params, 'Use Boundary Weights', found)) then
      if (.not. acc_pen_visited) then
        call CalculateNodalWeights(acc_solver, .true., VarName='Boundary Weights')
        weightsol => VariableGet(acc_model % variables, 'Boundary Weights')
        weightsol % output = .false.
      else
        weightsol => VariableGet(acc_model % variables, 'Boundary Weights')
      end if
      gotweight = associated(weightsol)
      if (gotweight) then
        call info('NonPhysicalMeshSolve', &
            'Using Boundary Weights for setting nodal penalty', level=8)
      else
        call fatal('NonPhysicalMeshSolve', 'Boundary Weights not present!')
      end if
    end if
    localpenalty = ListCheckPresentAnyBC(acc_model, 'Apply Nodal Penalty')
    if (.not. localpenalty) then
      localpenalty = ListCheckPresentAnyBodyForce(acc_model, 'Apply Nodal Penalty')
    end if
    if (localpenalty) then
      call info('NonPhysicalMeshSolve', 'Setting permutation for local penalties')
      allocate(penaltyperm(acc_solver % mesh % numberofnodes))
      call MakePermUsingMask(acc_model, acc_solver, acc_solver % mesh, &
          'Apply Nodal Penalty', .false., penaltyperm, i)
    end if
    if (.not. acc_got_field) then
      refsurface = ListGetCReal(params, 'Reference Surface')
    end if
    do i = 1, acc_solver % mesh % numberofnodes
      j = acc_meshperm(i)
      if (j == 0) cycle
      j2 = acc_targetsol % perm(i)
      if (j2 == 0) cycle
      if (localpenalty) then
        if (penaltyperm(i) == 0) cycle
      end if
      if (acc_got_field) then
        do k = 1, acc_dim
          extdisp(k) = acc_targetsol % values(acc_dim * (j2 - 1) + k)
        end do
      else
        surface = acc_targetsol % values(j2)
        if (acc_got_grad) then
          j2 = acc_gradsol % perm(i)
          do k = 1, acc_dim
            gradsurface(k) = acc_gradsol % values(acc_dim * (j2 - 1) + k)
          end do
        else
          call fatal('NonphysicalMeshSolve', 'You should provide GradSol!')
        end if
        absgrad2 = sum(gradsurface(1:acc_dim)**2)
        do k = 1, acc_dim
          extdisp(k) = (refsurface - surface) * gradsurface(k) / absgrad2
        end do
      end if
      do k = 1, acc_dim
        ind = 3 * (j - 1) + k
        coeff = targetcoeff
        if (gotweight) coeff = coeff * weightsol % values(weightsol % perm(i))
        acc_solver % matrix % rhs(ind) = acc_solver % matrix % rhs(ind) + &
            targetcoeff * extdisp(k)
        call CRS_AddToMatrixElement(acc_solver % matrix, ind, ind, targetcoeff)
      end do
    end do
    if (localpenalty) deallocate(penaltyperm)
    acc_pen_visited = .true.
  end subroutine nma_nodal_penalty

  subroutine nma_dirichlet() bind(c, name="nma_dirichlet")
    call DefaultDirichletBCs()
    call info('NonPhysicalMeshSolve', 'Set boundaries done', level=4)
  end subroutine nma_dirichlet

  function nma_default_solve() result(norm) bind(c, name="nma_default_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function nma_default_solve

  subroutine nma_nonlin_relax() bind(c, name="nma_nonlin_relax")
    type(valuelist_t), pointer :: params
    logical :: found
    real(c_double) :: relax
    params => GetSolverParams()
    relax = ListGetCReal(params, 'Nonlinear System Relaxation Factor', found)
    if (found) acc_meshupdate = relax * acc_meshupdate
  end subroutine nma_nonlin_relax

  function nma_mesh_relax() result(r) bind(c, name="nma_mesh_relax")
    real(c_double) :: r
    type(valuelist_t), pointer :: params
    logical :: found
    params => GetSolverParams()
    r = ListGetCReal(params, 'Mesh Relaxation Factor', found)
    if (.not. found) r = 1.0_c_double
  end function nma_mesh_relax

  ! -- mesh coordinate update + node-pointer restore -----------------------------------
  subroutine nma_displace(relax) bind(c, name="nma_displace")
    real(c_double), value :: relax
    integer :: i, j, n
    n = size(acc_meshperm)
    if (acc_moving) then
      do i = 1, n
        j = acc_meshperm(i)
        if (j == 0) cycle
        if (acc_dim == 2) then
          acc_truex(i) = acc_truex(i) + &
              relax * (acc_meshupdate(2 * j - 1) - acc_prev(2 * j - 1))
          acc_truey(i) = acc_truey(i) + &
              relax * (acc_meshupdate(2 * j) - acc_prev(2 * j))
        else
          acc_truex(i) = acc_truex(i) + &
              relax * (acc_meshupdate(3 * j - 2) - acc_prev(3 * j - 2))
          acc_truey(i) = acc_truey(i) + &
              relax * (acc_meshupdate(3 * j - 1) - acc_prev(3 * j - 1))
          acc_truez(i) = acc_truez(i) + &
              relax * (acc_meshupdate(3 * j) - acc_prev(3 * j))
        end if
      end do
    else
      do i = 1, n
        j = acc_meshperm(i)
        if (j == 0) cycle
        if (acc_dim == 2) then
          acc_truex(i) = acc_truex(i) + relax * acc_meshupdate(2 * j - 1)
          acc_truey(i) = acc_truey(i) + relax * acc_meshupdate(2 * j)
        else
          acc_truex(i) = acc_truex(i) + relax * acc_meshupdate(3 * j - 2)
          acc_truey(i) = acc_truey(i) + relax * acc_meshupdate(3 * j - 1)
          acc_truez(i) = acc_truez(i) + relax * acc_meshupdate(3 * j)
        end if
      end do
    end if
    if (.not. acc_cumulative) then
      acc_solver % mesh % nodes % x => acc_truex
      acc_solver % mesh % nodes % y => acc_truey
      acc_solver % mesh % nodes % z => acc_truez
    end if
  end subroutine nma_displace

end module nonphysicalmeshsolve_acc
