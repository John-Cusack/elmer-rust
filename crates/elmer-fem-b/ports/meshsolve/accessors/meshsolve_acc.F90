! T021-B MeshSolve port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/MeshSolve.F90, compiled against the built .mod
! files of the target configure (A2 shipped mechanism, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module mirrors MeshSolver's own SAVE state (velocities, perms, SAVE
! scalars, Nodes/IP contexts of both internal procedures) with identical
! lifetimes.
!
! Bit-identity rules followed throughout:
! - Whole-array fills, SUM reductions, the (p,q,i,j) element assembly
!   updates, the symmetric fill, the bubble-DOF zeroing, NormalVector,
!   formatted WRITE(Message,...)/Info output, and the timestep velocity
!   SELECT CASE stay in Fortran wrappers with the exact statement order
!   of the original, so bit-identity holds by construction wherever the
!   call sequence matches.
! - The Rust side owns the SAVE buffers (STIFF, FORCE, LOAD,
!   ElasticModulus, PoissonRatio, Alpha, Beta, TPerm) and passes (ld,ld)
!   / exact-shape dummies so no (nd,nd)-over-(ld,ld) stride mismatch can
!   occur; wrappers update (1:ld,1:ld) sections exactly as the original
!   updates whole arrays.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module meshsolve_acc
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
  ! Frame pointers re-pointed on every entry (mirrors the MeshSolver
  ! locals, including the SAVE alias pair MeshVelocity/MeshVeloPerm).
  type(variable_t), pointer, save :: acc_meshsol => null()
  type(variable_t), pointer, save :: acc_stresssol => null()
  real(c_double), pointer, save :: acc_meshupdate(:) => null()
  real(c_double), pointer, save :: acc_displacement(:) => null()
  real(c_double), pointer, save :: acc_meshvelocity(:) => null()
  integer, pointer, save :: acc_meshperm(:) => null()
  integer, pointer, save :: acc_stressperm(:) => null()
  integer, pointer, save :: acc_meshveloperm(:) => null()
  integer, save :: acc_stdofs = 0
  integer, save :: acc_msdofs = 0
  integer, save :: acc_dim = 0
  real(c_double), save :: acc_dt = 0.0_c_double
  integer, save :: acc_transient = 0
  ! SAVE scalars (mirror the original's SAVE list lifetimes).
  logical, save :: acc_firsttime = .TRUE.
  logical, save :: acc_skipfirst = .FALSE.
  real(c_double), save :: acc_at0 = 0.0_c_double
  ! LocalMatrix context (mirrors its SAVE Nodes; IP is reassigned per
  ! call exactly as the original's local is).
  type(nodes_t), target, save :: acc_nodes_lm
  type(gaussintegrationpoints_t), save :: acc_ip_lm
  ! MeshBoundary context (mirrors its SAVE Nodes the same way).
  type(nodes_t), target, save :: acc_nodes_mb
  type(gaussintegrationpoints_t), save :: acc_ip_mb

contains

  ! -- entry state ----------------------------------------------------
  subroutine ms_begin(h_model, h_solver, dt, transient) bind(c, name="ms_begin")
    type(c_ptr), value :: h_model, h_solver
    real(c_double), value :: dt
    integer(c_int), value :: transient
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    acc_dt = dt
    acc_transient = transient
  end subroutine ms_begin

  function ms_assoc_matrix() result(r) bind(c, name="ms_assoc_matrix")
    integer(c_int) :: r
    r = merge(1, 0, associated(acc_solver % matrix))
  end function ms_assoc_matrix

  ! -- init entry body (MeshSolver_Init transcribed) -------------------
  subroutine ms_init_body() bind(c, name="ms_init_body")
    type(valuelist_t), pointer :: params
    integer :: dim
    logical :: found, calculate
    params => GetSolverParams()
    dim = CoordinateSystemDimension()
    calculate = ListGetLogical(params, 'Compute Mesh Velocity', found)
    if (.not. found) calculate = .true.
    if (calculate) then
      if (acc_transient /= 0) then
        if (dim == 2) then
          call ListAddString(params, &
              NextFreeKeyword('Exported Variable', params), &
              '-dofs 2 Mesh Velocity')
        else
          call ListAddString(params, &
              NextFreeKeyword('Exported Variable', params), &
              '-dofs 3 Mesh Velocity')
        end if
      end if
    end if
  end subroutine ms_init_body

  ! -- frame: pointer capture, counts, dim adjust, stress/displace flags
  subroutine ms_frame() bind(c, name="ms_frame")
    type(variable_t), pointer :: meshes
    logical :: found
    nullify(acc_meshvelocity)
    if (acc_transient /= 0) then
      meshes => VariableGet(acc_solver % mesh % variables, 'Mesh Velocity')
      if (associated(meshes)) then
        acc_meshvelocity => meshes % values
        acc_meshveloperm => meshes % perm
      end if
    end if
    acc_meshsol => acc_solver % variable
    acc_meshperm => acc_meshsol % perm
    acc_stdofs = acc_meshsol % dofs
    acc_msdofs = acc_meshsol % dofs
    acc_meshupdate => acc_meshsol % values
    acc_dim = CoordinateSystemDimension()
    if (acc_stdofs < acc_dim) then
      call warn('MeshSolve', &
          'Displacement dofs smaller than dim, assuming reduced dimension!')
      acc_dim = acc_stdofs
    end if
    acc_stresssol => VariableGet(acc_solver % mesh % variables, 'Displacement')
    if (ListGetLogical(acc_solver % values, 'Ignore Displacement', found)) then
      nullify(acc_stresssol)
    end if
  end subroutine ms_frame

  function ms_local_nodes() result(n) bind(c, name="ms_local_nodes")
    integer(c_int) :: n
    n = count(acc_meshperm > 0)
  end function ms_local_nodes

  function ms_stdofs() result(n) bind(c, name="ms_stdofs")
    integer(c_int) :: n
    if (associated(acc_stresssol)) acc_stdofs = acc_stresssol % dofs
    n = acc_stdofs
  end function ms_stdofs

  function ms_dim() result(n) bind(c, name="ms_dim")
    integer(c_int) :: n
    n = acc_dim
  end function ms_dim

  function ms_stress_assoc() result(r) bind(c, name="ms_stress_assoc")
    integer(c_int) :: r
    r = merge(1, 0, associated(acc_stresssol))
  end function ms_stress_assoc

  function ms_displace_first() result(r) bind(c, name="ms_displace_first")
    integer(c_int) :: r
    logical :: found
    r = merge(1, 0, ListGetLogical(acc_solver % values, 'First Time Non-Zero', found))
  end function ms_displace_first

  subroutine ms_skip_displace(r) bind(c, name="ms_skip_displace")
    integer(c_int) :: r
    logical :: found
    if (ListGetLogical(acc_solver % values, 'Skip Displace Mesh', found)) then
      call info('MeshSolver', 'Skipping the displacement of mesh!', level=5)
      r = 1
    else
      r = 0
    end if
  end subroutine ms_skip_displace

  function ms_mesh_changed() result(r) bind(c, name="ms_mesh_changed")
    integer(c_int) :: r
    r = merge(1, 0, acc_solver % meshchanged)
  end function ms_mesh_changed

  function ms_nperm() result(n) bind(c, name="ms_nperm")
    integer(c_int) :: n
    n = size(acc_meshperm)
  end function ms_nperm

  function ms_max_elem_dofs() result(n) bind(c, name="ms_max_elem_dofs")
    integer(c_int) :: n
    n = acc_solver % mesh % maxelementdofs
  end function ms_max_elem_dofs

  ! -- TPerm setup (TPerm = MeshPerm, then zero shared nodes) -----------
  subroutine ms_tperm_fill(tperm, n) bind(c, name="ms_tperm_fill")
    integer(c_int), value :: n
    integer(c_int) :: tperm(n)
    integer :: i
    acc_stressperm => acc_stresssol % perm
    tperm(1:n) = acc_meshperm(1:n)
    do i = 1, n
      if (acc_stressperm(i) /= 0 .and. acc_meshperm(i) /= 0) tperm(i) = 0
    end do
  end subroutine ms_tperm_fill

  ! -- DisplaceMesh call shapes ------------------------------------------
  subroutine ms_disp_mu_tperm(tperm, n, dir) bind(c, name="ms_disp_mu_tperm")
    integer(c_int), value :: n, dir
    integer(c_int) :: tperm(n)
    call DisplaceMesh(acc_solver % mesh, acc_meshupdate, dir, tperm, acc_stdofs)
  end subroutine ms_disp_mu_tperm

  subroutine ms_disp_disp_stress(dir) bind(c, name="ms_disp_disp_stress")
    integer(c_int), value :: dir
    call DisplaceMesh(acc_solver % mesh, acc_displacement, dir, &
        acc_stressperm, acc_stdofs)
  end subroutine ms_disp_disp_stress

  subroutine ms_disp_disp_stress_norhs() bind(c, name="ms_disp_disp_stress_norhs")
    acc_displacement => acc_stresssol % values
    call DisplaceMesh(acc_solver % mesh, acc_displacement, 1, &
        acc_stressperm, acc_stdofs, .FALSE.)
  end subroutine ms_disp_disp_stress_norhs

  subroutine ms_disp_mu_meshperm(dir) bind(c, name="ms_disp_mu_meshperm")
    integer(c_int), value :: dir
    call DisplaceMesh(acc_solver % mesh, acc_meshupdate, dir, &
        acc_meshperm, acc_stdofs)
  end subroutine ms_disp_mu_meshperm

  subroutine ms_stress_vals() bind(c, name="ms_stress_vals")
    acc_stressperm => acc_stresssol % perm
    acc_displacement => acc_stresssol % values
  end subroutine ms_stress_vals

  ! -- permanent storage failure -----------------------------------------
  subroutine ms_fatal_alloc() bind(c, name="ms_fatal_alloc")
    call Fatal('MeshSolve', 'Memory allocation error.')
  end subroutine ms_fatal_alloc

  ! -- assembly driver -----------------------------------------------------
  subroutine ms_tick_begin() bind(c, name="ms_tick_begin")
    acc_at0 = RealTime()
    call info('MeshSolve', ' ', level=4)
    call info('MeshSolve', '-------------------------------------', level=4)
    call info('MeshSolve', 'MESH UPDATE SOLVER:', level=4)
    call info('MeshSolve', '-------------------------------------', level=4)
    call info('MeshSolve', ' ', level=4)
    call info('MeshSolve', 'Starting assembly...', level=4)
  end subroutine ms_tick_begin

  function ms_max_iter() result(n) bind(c, name="ms_max_iter")
    integer(c_int) :: n
    logical :: found
    n = GetInteger(acc_solver % values, 'Nonlinear System Max Iterations', found)
    if (.not. found) n = 1
  end function ms_max_iter

  subroutine ms_default_initialize() bind(c, name="ms_default_initialize")
    call DefaultInitialize()
  end subroutine ms_default_initialize

  function ms_nof_active() result(n) bind(c, name="ms_nof_active")
    integer(c_int) :: n
    n = acc_solver % numberofactiveelements
  end function ms_nof_active

  subroutine ms_tick_maybe(t, nactive) bind(c, name="ms_tick_maybe")
    integer(c_int), value :: t, nactive
    if (RealTime() - acc_at0 > 1.0) then
      write(message, '(a,i3,a)') '   Assembly: ', int(100.0 - 100.0 * &
          (nactive - t) / (1.0 * nactive)), ' % done'
      call info('MeshSolve', message, level=5)
      acc_at0 = RealTime()
    end if
  end subroutine ms_tick_maybe

  ! -- bulk-element context --------------------------------------------------
  subroutine ms_begin_active(t, n, nd, nb) bind(c, name="ms_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, nb
    acc_elem => GetActiveElement(t)
    nd = GetElementNOFDOFs()
    nb = GetElementNOFBDOFs()
    n = GetElementNOFNodes()
  end subroutine ms_begin_active

  subroutine ms_mat_props(n, nn, emod, pr) bind(c, name="ms_mat_props")
    integer(c_int), value :: n, nn
    real(c_double) :: emod(6, 6, nn), pr(nn)
    logical :: found
    acc_mat => GetMaterial()
    emod(1, 1, 1:n) = GetReal(acc_mat, 'Mesh Elastic Modulus', found)
    if (.not. found) emod(1, 1, 1:n) = GetReal(acc_mat, 'Youngs Modulus', found)
    if (.not. found) emod(1, 1, 1:n) = 1.0d0
    pr(1:n) = GetReal(acc_mat, 'Mesh Poisson Ratio', found)
    if (.not. found) pr(1:n) = GetReal(acc_mat, 'Poisson Ratio', found)
    if (.not. found) pr(1:n) = 0.25d0
  end subroutine ms_mat_props

  ! -- LocalMatrix (nodal Lame + per-IP assembly + finish) ----------------------
  subroutine ms_lm_nodes() bind(c, name="ms_lm_nodes")
    call GetElementNodes(acc_nodes_lm)
  end subroutine ms_lm_nodes

  subroutine ms_lame(n, nn, emod, pr, lame1, lame2, plane) bind(c, name="ms_lame")
    integer(c_int), value :: n, nn, plane
    real(c_double) :: emod(6, 6, nn), pr(nn), lame1(n), lame2(n)
    if (plane /= 0) then
      lame1(1:n) = emod(1, 1, 1:n) * pr(1:n) / ((1.0d0 - pr(1:n)**2))
    else
      lame1(1:n) = emod(1, 1, 1:n) * pr(1:n) / &
          ((1.0d0 + pr(1:n)) * (1.0d0 - 2.0d0 * pr(1:n)))
    end if
    lame2(1:n) = emod(1, 1, 1:n) / (2 * (1.0d0 + pr(1:n)))
  end subroutine ms_lame

  function ms_lm_ip_begin() result(n) bind(c, name="ms_lm_ip_begin")
    integer(c_int) :: n
    acc_ip_lm = GaussPoints(acc_elem)
    n = acc_ip_lm % n
  end function ms_lm_ip_begin

  subroutine ms_lm_ip_point(t, u, v, w, s) bind(c, name="ms_lm_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip_lm % u(t)
    v = acc_ip_lm % v(t)
    w = acc_ip_lm % w(t)
    s = acc_ip_lm % s(t)
  end subroutine ms_lm_ip_point

  subroutine ms_lm_eleminfo(u, v, w, n, nd, basis, dbasis, detj, stat) &
      bind(c, name="ms_lm_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: n, nd
    real(c_double) :: basis(nd), dbasis(nd, 3), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_elem, acc_nodes_lm, u, v, w, detj, &
        basis(1:nd), dbasis(1:nd, 1:3))
    stat = merge(1, 0, ok)
  end subroutine ms_lm_eleminfo

  function ms_dot(n, a, b) result(s) bind(c, name="ms_dot")
    integer(c_int), value :: n
    real(c_double) :: a(n), b(n), s
    s = sum(a(1:n) * b(1:n))
  end function ms_dot

  subroutine ms_lm_assemble(ld, dim, nd, s, l1, l2, dbasis, stiff) &
      bind(c, name="ms_lm_assemble")
    integer(c_int), value :: ld, dim, nd
    real(c_double), value :: s, l1, l2
    real(c_double) :: dbasis(nd, 3)
    real(c_double), target :: stiff(ld, ld)
    real(c_double), pointer :: a(:, :)
    integer :: i, j, p, q
    ! Exact transcription of the (p,q,i,j) update nest.
    do p = 1, nd
      do q = p, nd
        a => stiff(dim * (p - 1) + 1:dim * p, dim * (q - 1) + 1:dim * q)
        do i = 1, dim
          do j = 1, dim
            a(i, j) = a(i, j) + s * l1 * dbasis(q, j) * dbasis(p, i)
            a(i, i) = a(i, i) + s * l2 * dbasis(q, j) * dbasis(p, j)
            a(i, j) = a(i, j) + s * l2 * dbasis(q, i) * dbasis(p, j)
          end do
        end do
      end do
    end do
  end subroutine ms_lm_assemble

  ! Symmetric fill + bubble-DOF clamp (shared by LocalMatrix/MeshBoundary).
  subroutine ms_bubble(ld, dim, nd, n, nb, stiff, force) bind(c, name="ms_bubble")
    integer(c_int), value :: ld, dim, nd, n, nb
    real(c_double) :: stiff(ld, ld), force(ld)
    integer :: i, j, p, q
    do p = 1, dim * nd
      do q = 1, p - 1
        stiff(p, q) = stiff(q, p)
      end do
    end do
    if (nb == 0) then
      do p = max(n + 1, nd - acc_elem % bdofs + 1), nd
        do i = 1, dim
          j = (p - 1) * dim + i
          stiff(j, :) = 0.0d0
          stiff(:, j) = 0.0d0
          stiff(j, j) = 1.0d0
          force(j) = 0.0d0
        end do
      end do
    end if
  end subroutine ms_bubble

  subroutine ms_default_update(ld, g, f) bind(c, name="ms_default_update")
    integer(c_int), value :: ld
    real(c_double) :: g(ld, ld), f(ld)
    call DefaultUpdateEquations(g(1:ld, 1:ld), f(1:ld))
  end subroutine ms_default_update

  subroutine ms_finish_bulk() bind(c, name="ms_finish_bulk")
    call DefaultFinishBulkAssembly()
  end subroutine ms_finish_bulk

  ! -- boundary loop -------------------------------------------------------------
  function ms_doit_bc() result(r) bind(c, name="ms_doit_bc")
    integer(c_int) :: r
    logical :: doit
    doit = ListCheckPresentAnyBC(acc_model, 'Mesh Normal Force') .or. &
        ListCheckPrefixAnyBC(acc_model, 'Mesh Force')
    r = merge(1, 0, doit)
  end function ms_doit_bc

  function ms_nbnd() result(n) bind(c, name="ms_nbnd")
    integer(c_int) :: n
    n = acc_solver % mesh % numberofboundaryelements
  end function ms_nbnd

  subroutine ms_begin_boundary(t, n, nd, nb, active) bind(c, name="ms_begin_boundary")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, nb, active
    acc_elem => GetBoundaryElement(t)
    active = merge(1, 0, ActiveBoundaryElement())
    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
    nb = GetElementNOFBDOFs()
  end subroutine ms_begin_boundary

  function ms_bc_assoc() result(r) bind(c, name="ms_bc_assoc")
    integer(c_int) :: r
    acc_bc => GetBC()
    r = merge(1, 0, associated(acc_bc))
  end function ms_bc_assoc

  subroutine ms_bc_penalty(i, n, nn, alpha) bind(c, name="ms_bc_penalty")
    integer(c_int), value :: i, n, nn
    real(c_double) :: alpha(3, nn)
    logical :: found
    alpha(i, 1:n) = ListGetReal(acc_bc, &
        'Mesh Penalty Factor ' // trim(i2s(i)), n, &
        acc_elem % nodeindexes, found)
    if (found) then
      write(message, *) 'Mesh Penalty Factor ' // trim(i2s(i)) // ' =', alpha(i, 1)
      call info("MeshSolve", message, level=20)
    end if
  end subroutine ms_bc_penalty

  function ms_bc_loads(n, nn, load, beta) result(got) bind(c, name="ms_bc_loads")
    integer(c_int), value :: n, nn
    real(c_double) :: load(4, nn), beta(nn)
    integer(c_int) :: got
    logical :: found
    got = 0
    load(1, 1:n) = GetReal(acc_bc, 'Mesh Force 1', found)
    if (found) got = 1
    load(2, 1:n) = GetReal(acc_bc, 'Mesh Force 2', found)
    if (found) got = 1
    load(3, 1:n) = GetReal(acc_bc, 'Mesh Force 3', found)
    if (found) got = 1
    beta(1:n) = GetReal(acc_bc, 'Mesh Normal Force', found)
    if (found) got = 1
  end function ms_bc_loads

  ! -- MeshBoundary (per-IP assembly; Nodes SAVE + IP mirror) -----------------------
  subroutine ms_mb_nodes() bind(c, name="ms_mb_nodes")
    call GetElementNodes(acc_nodes_mb)
  end subroutine ms_mb_nodes

  function ms_mb_ip_begin() result(n) bind(c, name="ms_mb_ip_begin")
    integer(c_int) :: n
    acc_ip_mb = GaussPoints(acc_elem)
    n = acc_ip_mb % n
  end function ms_mb_ip_begin

  subroutine ms_mb_ip_point(t, u, v, w, s) bind(c, name="ms_mb_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip_mb % u(t)
    v = acc_ip_mb % v(t)
    w = acc_ip_mb % w(t)
    s = acc_ip_mb % s(t)
  end subroutine ms_mb_ip_point

  subroutine ms_mb_eleminfo(u, v, w, n, nd, basis, dbasis, detj, stat) &
      bind(c, name="ms_mb_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: n, nd
    real(c_double) :: basis(nd), dbasis(nd, 3), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_elem, acc_nodes_mb, u, v, w, detj, &
        basis(1:nd), dbasis(1:nd, 1:3))
    stat = merge(1, 0, ok)
  end subroutine ms_mb_eleminfo

  subroutine ms_mb_assemble(ld, dim, n, nd, u, v, w, s, load, alpha, beta, &
      basis, stiff, force) bind(c, name="ms_mb_assemble")
    integer(c_int), value :: ld, dim, n, nd
    real(c_double), value :: u, v, w, s
    real(c_double) :: load(4, n), alpha(3, n), beta(n)
    real(c_double) :: basis(nd), stiff(ld, ld), force(ld)
    real(c_double) :: loadatip(3), avec(3), normal(3)
    integer :: i, p, q
    loadatip = 0.0d0
    do i = 1, dim
      loadatip(i) = sum(load(i, 1:n) * basis(1:n))
      avec(i) = sum(alpha(i, 1:n) * basis(1:n))
    end do
    normal = NormalVector(acc_elem, acc_nodes_mb, u, v, .TRUE.)
    loadatip = loadatip + sum(beta(1:n) * basis(1:n)) * normal
    do p = 1, nd
      do q = 1, nd
        do i = 1, dim
          stiff((p - 1) * dim + i, (q - 1) * dim + i) = &
              stiff((p - 1) * dim + i, (q - 1) * dim + i) + &
              s * avec(i) * basis(q) * basis(p)
        end do
      end do
    end do
    do q = 1, nd
      do i = 1, dim
        force((q - 1) * dim + i) = force((q - 1) * dim + i) + &
            s * basis(q) * loadatip(i)
      end do
    end do
  end subroutine ms_mb_assemble

  ! -- solve + tail ------------------------------------------------------------------
  subroutine ms_finish_assembly() bind(c, name="ms_finish_assembly")
    call DefaultFinishBoundaryAssembly()
    call DefaultFinishAssembly()
    call info('MeshSolve', 'Assembly done', level=4)
  end subroutine ms_finish_assembly

  subroutine ms_dirichlet() bind(c, name="ms_dirichlet")
    call DefaultDirichletBCs()
    call info('MeshSolve', 'Set boundaries done', level=4)
  end subroutine ms_dirichlet

  function ms_default_solve() result(norm) bind(c, name="ms_default_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function ms_default_solve

  function ms_converged() result(r) bind(c, name="ms_converged")
    integer(c_int) :: r
    r = merge(1, 0, DefaultConverged())
  end function ms_converged

  ! -- transient mesh velocity (whole block transcribed; SAVE lifetimes mirror) ---------
  subroutine ms_velocity() bind(c, name="ms_velocity")
    logical :: compute, found
    integer :: k, j, m, i, dofs, l
    compute = ListGetLogical(acc_solver % values, 'Compute Mesh Velocity', found)
    if (.not. found) compute = .TRUE.
    acc_skipfirst = .FALSE.
    if (compute .and. acc_firsttime) then
      acc_skipfirst = ListGetLogical(acc_solver % values, 'Skip First Mesh Velocity', found)
      if (.not. found) then
        acc_skipfirst = .FALSE.
      else
        call info('MeshSolve', 'Skipping computation of initial Mesh Velocity', level=3)
      end if
      acc_firsttime = .FALSE.
    end if
    if (compute .and. (.not. (acc_skipfirst))) then
      k = min(size(acc_solver % variable % prevvalues, 2), acc_solver % donetime)
      j = ListGetInteger(acc_solver % values, 'Compute Mesh Velocity Order', found)
      if (found) then
        k = min(k, j)
      else
        k = 1
      end if
      dofs = acc_msdofs
      do i = 1, acc_solver % mesh % numberofnodes
        if (acc_meshveloperm(i) <= 0) cycle
        j = acc_meshperm(i) * dofs
        l = acc_meshveloperm(i) * dofs
        select case (k)
        case (0)
          acc_meshvelocity = 0._dp
        case (1)
          do m = 0, dofs - 1
            acc_meshvelocity(l - m) = (acc_meshupdate(j - m) - &
                acc_solver % variable % prevvalues(j - m, 1)) / acc_dt
          end do
        case (2)
          do m = 0, dofs - 1
            acc_meshvelocity(l - m) = ( &
                acc_meshupdate(j - m) - (4.0d0 / 3.0d0) * acc_solver % variable % prevvalues(j - m, 1) &
                + (1.0d0 / 3.0d0) * acc_solver % variable % prevvalues(j - m, 2)) / acc_dt
          end do
        case default
          do m = 0, dofs - 1
            acc_meshvelocity(l - m) = ( &
                acc_meshupdate(j - m) - (18.0d0 / 11.0d0) * acc_solver % variable % prevvalues(j - m, 1) &
                + (9.0d0 / 11.0d0) * acc_solver % variable % prevvalues(j - m, 2) &
                - (2.0d0 / 11.0d0) * acc_solver % variable % prevvalues(j - m, 3)) / acc_dt
          end do
        end select
      end do
    else if (associated(acc_meshvelocity) .and. (.not. (acc_skipfirst))) then
      acc_meshvelocity = 0.0d0
    end if
  end subroutine ms_velocity

end module meshsolve_acc
