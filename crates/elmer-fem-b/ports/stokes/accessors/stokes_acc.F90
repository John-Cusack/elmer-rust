! T021-B Stokes port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/Stokes.F90, compiled against the built .mod
! files of the target configure (Poisson-pattern accessor, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module mirrors StokesSolver's own locals (current element, Nodes SAVE
! state, Indexes SAVE array, TractionBCIndeces SAVE array, block-solver
! handles) with identical lifetimes.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module stokes_acc
  use iso_c_binding
  use types
  use defutils
  use lists
  use integration
  use matrixassembly
  use eleminfo
  use elementbasis
  use elementgeometry
  use coordinatesystems
  use solvecore
  use crsmatrix
  use bandmatrix
  use generalutils
  use messages
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(model_t), pointer, save :: acc_model => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(valuelist_t), pointer, save :: acc_bf => null()
  type(valuelist_t), pointer, save :: acc_mat => null()
  type(valuelist_t), pointer, save :: acc_bc => null()
  type(matrix_t), pointer, save :: acc_amat => null()
  type(matrix_t), pointer, save :: acc_mmat => null()
  type(matrix_t), pointer, save :: acc_pmat => null()
  type(solver_t), pointer, save :: acc_vsolver => null()
  type(solver_t), pointer, save :: acc_jsolver => null()
  type(solver_t), pointer, save :: acc_psolver => null()
  ! Element-integration context (mirrors the SAVE Nodes + IP locals of
  ! LocalMatrix / LocalMatrixBoundary / SchurComplementMatrix; the three
  ! never overlap in time, so one shared pair has identical lifetimes).
  type(nodes_t), target, save :: acc_nodes
  type(gaussintegrationpoints_t), save :: acc_ip
  ! SAVE Indexes array (allocated once at first sizing).
  integer, allocatable, save :: acc_idx(:)
  ! SAVE TractionBCIndeces array (block path only; integer node numbers).
  integer, allocatable, save :: acc_traction(:)

contains

  subroutine sta_key_to_str(key, keylen, s)
    character(kind=c_char), intent(in) :: key(*)
    integer(c_int), value, intent(in) :: keylen
    character(len=*), intent(out) :: s
    integer :: i, m
    s = ' '
    m = min(keylen, len(s))
    do i = 1, m
      if (key(i) == c_null_char) exit
      s(i:i) = key(i)
    end do
  end subroutine sta_key_to_str

  subroutine sta_ensure_idx()
    integer :: pmax
    pmax = acc_solver % mesh % maxelementdofs
    if (.not. allocated(acc_idx)) allocate(acc_idx(pmax))
  end subroutine sta_ensure_idx

  ! -- entry state ----------------------------------------------------
  subroutine sta_begin(h_model, h_solver) bind(c, name="sta_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine sta_begin

  function sta_done_time() result(n) bind(c, name="sta_done_time")
    integer(c_int) :: n
    n = acc_solver % donetime
  end function sta_done_time

  function sta_dim() result(n) bind(c, name="sta_dim")
    integer(c_int) :: n
    n = CoordinateSystemDimension()
  end function sta_dim

  function sta_max_elem_dofs() result(n) bind(c, name="sta_max_elem_dofs")
    integer(c_int) :: n
    n = acc_solver % mesh % maxelementdofs
  end function sta_max_elem_dofs

  function sta_num_rows() result(n) bind(c, name="sta_num_rows")
    integer(c_int) :: n
    n = acc_solver % matrix % numberofrows
  end function sta_num_rows

  function sta_solver_dt() result(x) bind(c, name="sta_solver_dt")
    real(c_double) :: x
    x = acc_solver % dt
  end function sta_solver_dt

  function sta_nof_active() result(n) bind(c, name="sta_nof_active")
    integer(c_int) :: n
    n = GetNOFActive()
  end function sta_nof_active

  function sta_num_bnd_elems() result(n) bind(c, name="sta_num_bnd_elems")
    integer(c_int) :: n
    n = acc_solver % mesh % numberofboundaryelements
  end function sta_num_bnd_elems

  ! -- ValueList getters (list: 0 = Solver % Values, 1 = SolverParams,
  ! -- 2 = current BC) ------------------------------------------------
  function sta_list_logical(which, key, keylen, found) result(r) &
      bind(c, name="sta_list_logical")
    integer(c_int), value :: which, keylen
    character(kind=c_char) :: key(*)
    integer(c_int) :: found, r
    character(len=128) :: k
    logical :: f, v
    call sta_key_to_str(key, keylen, k)
    if (which == 0) then
      v = ListGetLogical(acc_solver % values, trim(k), f)
    else if (which == 1) then
      v = GetLogical(GetSolverParams(acc_solver), trim(k), f)
    else
      v = GetLogical(acc_bc, trim(k), f)
    end if
    r = merge(1, 0, v)
    found = merge(1, 0, f)
  end function sta_list_logical

  function sta_list_const_real(key, keylen, found, use_minv, minv, with_found) result(r) &
      bind(c, name="sta_list_const_real")
    character(kind=c_char) :: key(*)
    integer(c_int), value :: keylen, use_minv, with_found
    real(c_double), value :: minv
    integer(c_int) :: found
    real(c_double) :: r
    character(len=128) :: k
    logical :: f
    type(valuelistentry_t), pointer :: ptr
    call sta_key_to_str(key, keylen, k)
    if (with_found /= 0) then
      if (use_minv /= 0) then
        r = ListGetConstReal(acc_solver % values, trim(k), f, minv=minv)
      else
        r = ListGetConstReal(acc_solver % values, trim(k), f)
      end if
    else
      ! Mirror the no-Found call: the ListFind warning is emitted, then
      ! the found-status is probed silently for the caller's branch.
      if (use_minv /= 0) then
        r = ListGetConstReal(acc_solver % values, trim(k), minv=minv)
      else
        r = ListGetConstReal(acc_solver % values, trim(k))
      end if
      ptr => ListFind(acc_solver % values, trim(k), f)
      f = associated(ptr)
    end if
    found = merge(1, 0, f)
  end function sta_list_const_real

  function sta_list_integer(key, keylen, found, use_minv, minv, with_found) result(r) &
      bind(c, name="sta_list_integer")
    character(kind=c_char) :: key(*)
    integer(c_int), value :: keylen, use_minv, minv, with_found
    integer(c_int) :: found, r
    character(len=128) :: k
    logical :: f
    type(valuelistentry_t), pointer :: ptr
    call sta_key_to_str(key, keylen, k)
    if (with_found /= 0) then
      if (use_minv /= 0) then
        r = ListGetInteger(acc_solver % values, trim(k), f, minv=minv)
      else
        r = ListGetInteger(acc_solver % values, trim(k), f)
      end if
    else
      ! Mirror the no-Found call: the ListFind warning is emitted, then
      ! the found-status is probed silently for the caller's branch.
      if (use_minv /= 0) then
        r = ListGetInteger(acc_solver % values, trim(k), minv=minv)
      else
        r = ListGetInteger(acc_solver % values, trim(k))
      end if
      ptr => ListFind(acc_solver % values, trim(k), f)
      f = associated(ptr)
    end if
    found = merge(1, 0, f)
  end function sta_list_integer

  function sta_list_creal(key, keylen, found) result(r) &
      bind(c, name="sta_list_creal")
    character(kind=c_char) :: key(*)
    integer(c_int), value :: keylen
    integer(c_int) :: found
    real(c_double) :: r
    character(len=128) :: k
    logical :: f
    call sta_key_to_str(key, keylen, k)
    r = ListGetCReal(acc_solver % values, trim(k), f)
    found = merge(1, 0, f)
  end function sta_list_creal

  ! -- string getters (which: 0 = Nonlinear Iteration Method,
  ! -- 1 = Equation of Model % Solvers(idx)) ---------------------------
  subroutine sta_get_string(which, idx, buf, buflen, found) &
      bind(c, name="sta_get_string")
    integer(c_int), value :: which, idx, buflen
    character(kind=c_char) :: buf(*)
    integer(c_int) :: found
    character(len=max_name_len) :: s
    logical :: f
    integer :: i, m
    if (which == 0) then
      s = ListGetString(acc_solver % values, 'Nonlinear Iteration Method', f)
    else
      s = ListGetString(acc_model % solvers(idx) % values, 'Equation', f)
    end if
    found = merge(1, 0, f)
    do i = 1, buflen
      buf(i) = ' '
    end do
    m = min(buflen, len_trim(s))
    do i = 1, m
      buf(i) = s(i:i)
    end do
  end subroutine sta_get_string

  ! -- fatals ----------------------------------------------------------
  subroutine sta_fatal_alloc() bind(c, name="sta_fatal_alloc")
    call Fatal('NavierStokesSolver', 'Memory allocation error.')
  end subroutine sta_fatal_alloc

  subroutine sta_fatal_msg(which) bind(c, name="sta_fatal_msg")
    integer(c_int), value :: which
    if (which == 0) then
      call Fatal('NavierStokesSolver', &
          'Pressure Preconditioning matrix (S) was not found')
    else if (which == 1) then
      call Fatal('NavierStokesSolver', &
          'Pressure Mass matrix (M) was not found')
    else if (which == 2) then
      call Fatal('NavierStokesSolver', &
          'Velocity Preconditioning matrix (A) was not found')
    else if (which == 3) then
      call Fatal('NavierStokesSolver', &
          'Incompatible permutations for pressure preconditioning and primary unknown')
    else if (which == 4) then
      call Fatal('NavierStokesSolver', &
          'Incompatible permutations for divergence projection and primary unknown')
    else
      call Fatal('NavierStokesSolver', &
          'Incompatible permutations for divergence projection and primary unknown')
    end if
  end subroutine sta_fatal_msg

  ! -- block-preconditioning setup -------------------------------------
  subroutine sta_block_find() bind(c, name="sta_block_find")
    integer :: i, nlen, n
    character(len=max_name_len) :: eq
    logical :: found
    nullify(acc_amat, acc_mmat, acc_pmat)
    nullify(acc_vsolver, acc_jsolver, acc_psolver)
    do i = 1, acc_model % numberofsolvers
      eq = ListGetString(acc_model % solvers(i) % values, 'Equation', found)
      if (found) then
        nlen = len_trim(eq)
        if (eq(1:nlen) == 'pressure preconditioning') then
          acc_psolver => acc_model % solvers(i)
          acc_pmat => acc_psolver % matrix
        else if (eq(1:nlen) == 'divergence projection') then
          acc_jsolver => acc_model % solvers(i)
          acc_mmat => acc_jsolver % matrix
        else if (eq(1:nlen) == 'velocity preconditioning') then
          acc_vsolver => acc_model % solvers(i)
          acc_amat => acc_vsolver % matrix
        end if
      end if
    end do
    if (.not. associated(acc_pmat)) &
        call Fatal('NavierStokesSolver', &
        'Pressure Preconditioning matrix (S) was not found')
    if (.not. associated(acc_mmat)) &
        call Fatal('NavierStokesSolver', &
        'Pressure Mass matrix (M) was not found')
    if (.not. associated(acc_amat)) &
        call Fatal('NavierStokesSolver', &
        'Velocity Preconditioning matrix (A) was not found')
    n = size(acc_solver % variable % perm)
    found = all(acc_psolver % variable % perm(1:n) == &
        acc_solver % variable % perm(1:n))
    if (.not. found) call Fatal('NavierStokesSolver', &
        'Incompatible permutations for pressure preconditioning and primary unknown')
    found = all(acc_jsolver % variable % perm(1:n) == &
        acc_solver % variable % perm(1:n))
    if (.not. found) call Fatal('NavierStokesSolver', &
        'Incompatible permutations for divergence projection and primary unknown')
    found = all(acc_vsolver % variable % perm(1:n) == &
        acc_solver % variable % perm(1:n))
    if (.not. found) call Fatal('NavierStokesSolver', &
        'Incompatible permutations for divergence projection and primary unknown')
  end subroutine sta_block_find

  function sta_block_alloc() result(cnt) bind(c, name="sta_block_alloc")
    integer(c_int) :: cnt
    type(element_t), pointer :: element
    type(valuelist_t), pointer :: bc
    integer :: t, n, nd, i, nn, istat
    logical :: found, outflowbc
    call sta_ensure_idx()
    nn = 0
    i = 0
    do t = 1, acc_solver % mesh % numberofboundaryelements
      element => GetBoundaryElement(t)
      if (.not. ActiveBoundaryElement()) cycle
      n = GetElementNOFNodes()
      nd = GetElementDOFs(acc_idx)
      if (GetElementFamily() == 1) cycle
      bc => GetBC()
      if (associated(bc)) then
        outflowbc = ListGetLogical(bc, 'Outflow Boundary', found)
        if (.not. outflowbc) cycle
        if (any(acc_solver % variable % perm(acc_idx(1:nd)) == 0)) cycle
        i = i + 1
        nn = nn + nd
      end if
    end do
    if (nn > 0) then
      allocate(acc_traction(nn), stat=istat)
      if (istat /= 0) &
          call Fatal('NavierStokesSolver', 'Memory allocation error.')
    end if
    cnt = nn
  end function sta_block_alloc

  subroutine sta_block_fill() bind(c, name="sta_block_fill")
    type(element_t), pointer :: element
    type(valuelist_t), pointer :: bc
    integer :: t, n, nd, i, m, k
    logical :: found, outflowbc
    m = 0
    do t = 1, acc_solver % mesh % numberofboundaryelements
      element => GetBoundaryElement(t)
      if (.not. ActiveBoundaryElement()) cycle
      n = GetElementNOFNodes()
      nd = GetElementDOFs(acc_idx)
      if (GetElementFamily() == 1) cycle
      bc => GetBC()
      if (associated(bc)) then
        outflowbc = ListGetLogical(bc, 'Outflow Boundary', found)
        if (.not. outflowbc) cycle
        if (any(acc_solver % variable % perm(acc_idx(1:nd)) == 0)) cycle
        do i = 1, nd
          k = acc_solver % variable % perm(acc_idx(i))
          acc_traction(m+i) = k
          call ZeroRow(acc_pmat, k)
          call SetMatrixElement(acc_pmat, k, k, 1.0_dp)
        end do
        m = m + nd
      end if
    end do
  end subroutine sta_block_fill

  subroutine sta_zero_crs(which) bind(c, name="sta_zero_crs")
    integer(c_int), value :: which
    if (which == 0) then
      call CRS_ZeroMatrix(acc_amat)
    else if (which == 1) then
      call CRS_ZeroMatrix(acc_mmat)
    else
      call CRS_ZeroMatrix(acc_pmat)
    end if
  end subroutine sta_zero_crs

  ! -- timers -----------------------------------------------------------
  function sta_cpu_time() result(t) bind(c, name="sta_cpu_time")
    real(c_double) :: t
    t = CPUTime()
  end function sta_cpu_time

  function sta_real_time() result(t) bind(c, name="sta_real_time")
    real(c_double) :: t
    t = RealTime()
  end function sta_real_time

  ! -- assembly driver --------------------------------------------------
  subroutine sta_default_initialize() bind(c, name="sta_default_initialize")
    call DefaultInitialize()
  end subroutine sta_default_initialize

  subroutine sta_info_progress(t) bind(c, name="sta_info_progress")
    integer(c_int), value :: t
    integer :: nact
    nact = acc_solver % numberofactiveelements
    write(message, '(a,i3,a)') '   Assembly: ', int(100.0_dp - 100.0_dp * &
        (nact - t) / (nact)), ' % done'
    call Info('NavierStokesSolver', message, level=5)
  end subroutine sta_info_progress

  subroutine sta_info_nonlin(res) bind(c, name="sta_info_nonlin")
    real(c_double), value :: res
    write(message, '(a,E12.4)') 'Nonlinear iteration residual: ', res
    call Info('NavierStokesSolver', message, level=4)
  end subroutine sta_info_nonlin

  subroutine sta_info_snorm(norm) bind(c, name="sta_info_snorm")
    real(c_double), value :: norm
    write(message, '(a,E12.4)') 'Solution norm: ', norm
    call Info('NavierStokesSolver', message, level=4)
  end subroutine sta_info_snorm

  subroutine sta_info_atime(atime) bind(c, name="sta_info_atime")
    real(c_double), value :: atime
    write(message, '(a,F8.2)') ' Assembly:  (s)', atime
    call Info('NavierStokesSolver', message, level=4)
  end subroutine sta_info_atime

  subroutine sta_info_stime(stime) bind(c, name="sta_info_stime")
    real(c_double), value :: stime
    write(message, '(a,F8.2)') ' Solution:  (s)', stime
    call Info('NavierStokesSolver', message, level=4)
  end subroutine sta_info_stime

  subroutine sta_outer_timelevel(cur) bind(c, name="sta_outer_timelevel")
    integer(c_int), value :: cur
    print *, 'OuterIteration residuals for time level', cur
  end subroutine sta_outer_timelevel

  subroutine sta_print_msg(which) bind(c, name="sta_print_msg")
    integer(c_int), value :: which
    if (which == 0) then
      print *, 'Velocity solve...'
    else if (which == 1) then
      print *, 'Pressure mass matrix solve...'
    else
      print *, 'Pressure Laplace matrix solve...'
    end if
  end subroutine sta_print_msg

  subroutine sta_gcr_info0(res) bind(c, name="sta_gcr_info0")
    real(c_double), value :: res
    write(message, '(a,I4,ES12.3)') 'OuterIteration residual for iterate', 0, res
    call Info('NavierStokesSolver', message, level=4)
  end subroutine sta_gcr_info0

  subroutine sta_gcr_info(k, res) bind(c, name="sta_gcr_info")
    integer(c_int), value :: k
    real(c_double), value :: res
    write(message, '(a,I4,ES12.3)') 'OuterIteration residual for iterate', k, res
    call Info('NavierStokesSolver', message, level=4)
  end subroutine sta_gcr_info

  subroutine sta_bicg_write(i, res) bind(c, name="sta_bicg_write")
    integer(c_int), value :: i
    real(c_double), value :: res
    write(*, '(I4,ES12.3)') i, res
  end subroutine sta_bicg_write

  ! -- active-element context -------------------------------------------
  subroutine sta_begin_active(t, n, nd, nb) bind(c, name="sta_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, nb
    call sta_ensure_idx()
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes()
    nd = GetElementDOFs(acc_idx)
    nb = GetElementNOFBDOFs()
  end subroutine sta_begin_active

  function sta_bodyforce() result(r) bind(c, name="sta_bodyforce")
    integer(c_int) :: r
    acc_bf => GetBodyForce()
    r = merge(1, 0, associated(acc_bf))
  end function sta_bodyforce

  subroutine sta_load_force(k, n, load, found) bind(c, name="sta_load_force")
    integer(c_int), value :: k, n
    real(c_double) :: load(n)
    integer(c_int) :: found
    character(len=16) :: fname
    logical :: f
    write(fname, '(a,i1)') 'Force ', k
    load(1:n) = GetReal(acc_bf, fname, f)
    found = merge(1, 0, f)
  end subroutine sta_load_force

  subroutine sta_mat_rho_mu(n, rho, mu) bind(c, name="sta_mat_rho_mu")
    integer(c_int), value :: n
    real(c_double) :: rho(n), mu(n)
    acc_mat => GetMaterial()
    rho(1:n) = GetReal(acc_mat, 'Density')
    mu(1:n) = GetReal(acc_mat, 'Viscosity')
  end subroutine sta_mat_rho_mu

  ! -- velocity gather (which: 0 = PrevValues, 1 = Values) --------------
  subroutine sta_gather_vel(which, dim, ldv, nd, vel) &
      bind(c, name="sta_gather_vel")
    integer(c_int), value :: which, dim, ldv, nd
    real(c_double) :: vel(ldv, *)
    integer :: i, j, dofs
    dofs = acc_solver % variable % dofs
    if (which == 0) then
      do j = 1, nd
        do i = 1, dim
          vel(i, j) = acc_solver % variable % prevvalues( &
              dofs*(acc_solver % variable % perm(acc_idx(j))-1)+i, 1)
        end do
      end do
    else
      do j = 1, nd
        do i = 1, dim
          vel(i, j) = acc_solver % variable % values( &
              dofs*(acc_solver % variable % perm(acc_idx(j))-1)+i)
        end do
      end do
    end if
  end subroutine sta_gather_vel

  ! -- element integration context (shared Nodes/IP) --------------------
  subroutine sta_ip_begin(nip, ch) bind(c, name="sta_ip_begin")
    integer(c_int) :: nip
    real(c_double) :: ch
    call GetElementNodes(acc_nodes)
    acc_ip = GaussPoints(acc_elem)
    nip = acc_ip % n
    ch = ElementDiameter(acc_elem, acc_nodes)
  end subroutine sta_ip_begin

  subroutine sta_schur_begin(nip) bind(c, name="sta_schur_begin")
    integer(c_int) :: nip
    call GetElementNodes(acc_nodes)
    acc_ip = GaussPoints(acc_elem)
    nip = acc_ip % n
  end subroutine sta_schur_begin

  subroutine sta_ip_point(t, u, v, w, s) bind(c, name="sta_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip % u(t)
    v = acc_ip % v(t)
    w = acc_ip % w(t)
    s = acc_ip % s(t)
  end subroutine sta_ip_point

  subroutine sta_lm_eleminfo(u, v, w, n, nd, basis, dbasis, detj, stat) &
      bind(c, name="sta_lm_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: n, nd
    real(c_double) :: basis(nd), dbasis(nd, 3), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_elem, acc_nodes, u, v, w, detj, &
        basis(1:nd), dbasis(1:nd, 1:3))
    stat = merge(1, 0, ok)
  end subroutine sta_lm_eleminfo

  subroutine sta_lm_eleminfo2(u, v, w, n, nd, basis, dbasis, ddb, detj, stat) &
      bind(c, name="sta_lm_eleminfo2")
    real(c_double), value :: u, v, w
    integer(c_int), value :: n, nd
    real(c_double) :: basis(nd), dbasis(nd, 3), ddb(nd, 3, 3), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_elem, acc_nodes, u, v, w, detj, &
        basis(1:nd), dbasis(1:nd, 1:3), ddb(1:nd, 1:3, 1:3), .false.)
    stat = merge(1, 0, ok)
  end subroutine sta_lm_eleminfo2

  ! -- order-sensitive reductions (built at -O0, same evaluation order
  ! -- as the reference module build) -----------------------------------
  subroutine sta_lm_velo(dim, ldv, nd, vel, basis, out) &
      bind(c, name="sta_lm_velo")
    integer(c_int), value :: dim, ldv, nd
    real(c_double) :: vel(ldv, *), basis(nd), out(dim)
    out(1:dim) = matmul(vel(1:dim, 1:nd), basis(1:nd))
  end subroutine sta_lm_velo

  subroutine sta_lm_w(need12, ldv, nd, vel, dd, w) bind(c, name="sta_lm_w")
    integer(c_int), value :: need12, ldv, nd
    real(c_double) :: vel(ldv, *), dd(nd, *), w(3)
    w(3) = sum(vel(2, 1:nd) * dd(1:nd, 1)) - &
        sum(vel(1, 1:nd) * dd(1:nd, 2))
    if (need12 /= 0) then
      w(1) = sum(vel(3, 1:nd) * dd(1:nd, 2)) - &
          sum(vel(2, 1:nd) * dd(1:nd, 3))
      w(2) = sum(vel(1, 1:nd) * dd(1:nd, 3)) - &
          sum(vel(3, 1:nd) * dd(1:nd, 1))
    end if
  end subroutine sta_lm_w

  subroutine sta_lm_grad(ldv, nd, vel, dd, g) bind(c, name="sta_lm_grad")
    integer(c_int), value :: ldv, nd
    real(c_double) :: vel(ldv, *), dd(nd, *), g(3, 3)
    g(1, 1) = sum(vel(1, 1:nd) * dd(1:nd, 1))
    g(1, 2) = sum(vel(1, 1:nd) * dd(1:nd, 2))
    g(2, 1) = sum(vel(2, 1:nd) * dd(1:nd, 1))
    g(2, 2) = sum(vel(2, 1:nd) * dd(1:nd, 2))
  end subroutine sta_lm_grad

  subroutine sta_lm_murho(n, basis, nodalmu, nodalrho, mu, rho) &
      bind(c, name="sta_lm_murho")
    integer(c_int), value :: n
    real(c_double) :: basis(n), nodalmu(n), nodalrho(n), mu, rho
    mu = sum(basis(1:n) * nodalmu(1:n))
    rho = sum(basis(1:n) * nodalrho(1:n))
  end subroutine sta_lm_murho

  subroutine sta_lm_loadatip(n, dim, ldl, load, basis, out) &
      bind(c, name="sta_lm_loadatip")
    integer(c_int), value :: n, dim, ldl
    real(c_double) :: load(ldl, 4), basis(n), out(4)
    out(1:dim+1) = matmul(basis(1:n), load(1:n, 1:dim+1))
  end subroutine sta_lm_loadatip

  subroutine sta_lm_hk_mk(hk, mk) bind(c, name="sta_lm_hk_mk")
    real(c_double) :: hk, mk
    hk = acc_elem % hk
    mk = acc_elem % stabilizationmk
  end subroutine sta_lm_hk_mk

  subroutine sta_lb_normal(u, v, nl) bind(c, name="sta_lb_normal")
    real(c_double), value :: u, v
    real(c_double) :: nl(3)
    nl(1:3) = Normalvector(acc_elem, acc_nodes, u, v, .true.)
  end subroutine sta_lb_normal

  ! -- shared assembly tail ----------------------------------------------
  subroutine sta_default_1st_order(ld, m, a, f) &
      bind(c, name="sta_default_1st_order")
    integer(c_int), value :: ld
    real(c_double) :: m(ld, ld), a(ld, ld), f(ld)
    call Default1stOrderTime(m(1:ld, 1:ld), a(1:ld, 1:ld), f(1:ld))
  end subroutine sta_default_1st_order

  subroutine sta_default_update(ld, g, f) bind(c, name="sta_default_update")
    integer(c_int), value :: ld
    real(c_double) :: g(ld, ld), f(ld)
    call DefaultUpdateEquations(g(1:ld, 1:ld), f(1:ld))
  end subroutine sta_default_update

  subroutine sta_update_precond(which, ld, nd, ndofs, local) &
      bind(c, name="sta_update_precond")
    integer(c_int), value :: which, ld, nd, ndofs
    real(c_double) :: local(ld, *)
    type(matrix_t), pointer :: sm
    integer :: ni(nd)
    if (which == 0) then
      sm => acc_amat
    else if (which == 1) then
      sm => acc_mmat
    else
      sm => acc_pmat
    end if
    ni(1:nd) = acc_solver % variable % perm(acc_idx(1:nd))
    select case (sm % format)
    case (matrix_crs)
      call CRS_GlueLocalMatrix(sm, nd, ndofs, ni, &
          local(1:nd*ndofs, 1:nd*ndofs))
    case (matrix_band, matrix_sband)
      call Band_GlueLocalMatrix(sm, nd, ndofs, ni, &
          local(1:nd*ndofs, 1:nd*ndofs))
    end select
  end subroutine sta_update_precond

  ! -- boundary-element context ------------------------------------------
  subroutine sta_begin_boundary(t, n, nd, active, family) &
      bind(c, name="sta_begin_boundary")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, active, family
    call sta_ensure_idx()
    acc_elem => GetBoundaryElement(t)
    active = merge(1, 0, ActiveBoundaryElement())
    n = GetElementNOFNodes()
    nd = GetElementDOFs(acc_idx)
    family = GetElementFamily()
  end subroutine sta_begin_boundary

  subroutine sta_bc() bind(c, name="sta_bc")
    acc_bc => GetBC()
  end subroutine sta_bc

  function sta_bc_associated() result(r) bind(c, name="sta_bc_associated")
    integer(c_int) :: r
    r = merge(1, 0, associated(acc_bc))
  end function sta_bc_associated

  subroutine sta_default_finish_bulk() bind(c, name="sta_default_finish_bulk")
    call DefaultFinishBulkAssembly()
  end subroutine sta_default_finish_bulk

  subroutine sta_default_finish_assembly() &
      bind(c, name="sta_default_finish_assembly")
    call DefaultFinishAssembly()
  end subroutine sta_default_finish_assembly

  ! -- internal Dirichlet handling (replicates SetBoundaryConditions) ---
  subroutine sta_set_bcs(mat, dof, ndofs, with_rhs) &
      bind(c, name="sta_set_bcs")
    integer(c_int), value :: mat, dof, ndofs, with_rhs
    type(matrix_t), pointer :: sm
    type(element_t), pointer :: curelm
    integer, pointer :: nodeindexes(:)
    integer :: i, j, k, t, n, nd
    logical :: gotit
    real(kind=dp) :: work(acc_model % maxelementnodes)
    character(:), allocatable :: nm
    call sta_ensure_idx()
    if (mat == 0) then
      sm => acc_solver % matrix
    else
      sm => acc_amat
    end if
    nm = ComponentName(acc_solver % variable % name, dof)
    do t = acc_model % numberofbulkelements + 1, &
        acc_model % numberofbulkelements + acc_model % numberofboundaryelements
      curelm => acc_model % elements(t)
      acc_model % currentelement => acc_model % elements(t)
      n = GetElementNOFNodes()
      nd = GetElementDOFs(acc_idx)
      nodeindexes => curelm % nodeindexes(1:n)
      do i = 1, acc_model % numberofbcs
        if (curelm % boundaryinfo % constraint == acc_model % bcs(i) % tag) then
          work(1:n) = ListGetReal(acc_model % bcs(i) % values, &
              nm, n, nodeindexes, gotit)
          if (gotit) then
            do j = 1, n
              k = acc_solver % variable % perm(nodeindexes(j))
              if (k > 0) then
                k = ndofs * (k-1) + dof
                call ZeroRow(sm, k)
                call SetMatrixElement(sm, k, k, 1.0_dp)
                if (with_rhs /= 0) sm % rhs(k) = work(j)
              end if
            end do
            do j = n + 1, nd
              k = acc_solver % variable % perm(acc_idx(j))
              k = ndofs * (k-1) + dof
              call ZeroRow(sm, k)
              call SetMatrixElement(sm, k, k, 1.0_dp)
              if (with_rhs /= 0) sm % rhs(k) = 0.0_dp
            end do
          end if
        end if
      end do
    end do
  end subroutine sta_set_bcs

  ! -- solve phase -------------------------------------------------------
  subroutine sta_prev_sol(which, n, prev) bind(c, name="sta_prev_sol")
    integer(c_int), value :: which, n
    real(c_double) :: prev(n)
    if (which == 0) then
      prev(1:n) = acc_solver % variable % prevvalues(1:n, 1)
    else
      prev(1:n) = acc_solver % variable % values(1:n)
    end if
  end subroutine sta_prev_sol

  subroutine sta_matvec_rhs(n, prev, res) bind(c, name="sta_matvec_rhs")
    integer(c_int), value :: n
    real(c_double) :: prev(n), res(n)
    call CRS_MatrixVectorMultiply(acc_solver % matrix, prev(1:n), res(1:n))
    res(1:n) = acc_solver % matrix % rhs(1:n) - res(1:n)
  end subroutine sta_matvec_rhs

  function sta_pnorm_vec(n, x) result(s) bind(c, name="sta_pnorm_vec")
    integer(c_int), value :: n
    real(c_double) :: x(n), s
    s = sqrt(dot_product(x(1:n), x(1:n)))
  end function sta_pnorm_vec

  function sta_rhs_pnorm() result(s) bind(c, name="sta_rhs_pnorm")
    real(c_double) :: s
    integer :: n
    n = acc_solver % matrix % numberofrows
    s = sqrt(dot_product(acc_solver % matrix % rhs(1:n), &
        acc_solver % matrix % rhs(1:n)))
  end function sta_rhs_pnorm

  function sta_values_pnorm(n) result(s) bind(c, name="sta_values_pnorm")
    integer(c_int), value :: n
    real(c_double) :: s
    s = sqrt(dot_product(acc_solver % variable % values(1:n), &
        acc_solver % variable % values(1:n)))
  end function sta_values_pnorm

  function sta_pdot(n, x, y) result(s) bind(c, name="sta_pdot")
    integer(c_int), value :: n
    real(c_double) :: x(n), y(n), s
    s = dot_product(x(1:n), y(1:n))
  end function sta_pdot

  subroutine sta_relax(n, prev, factor) bind(c, name="sta_relax")
    integer(c_int), value :: n
    real(c_double) :: prev(n)
    real(c_double), value :: factor
    acc_solver % variable % values(1:n) = (1.0_dp - factor) * prev(1:n) + &
        factor * acc_solver % variable % values(1:n)
  end subroutine sta_relax

  function sta_default_solve() result(norm) bind(c, name="sta_default_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function sta_default_solve

  ! -- GCR data movement --------------------------------------------------
  subroutine sta_values_copy(x, n) bind(c, name="sta_values_copy")
    integer(c_int), value :: n
    real(c_double) :: x(n)
    x(1:n) = acc_solver % variable % values(1:n)
  end subroutine sta_values_copy

  subroutine sta_rhs_copy(b, n) bind(c, name="sta_rhs_copy")
    integer(c_int), value :: n
    real(c_double) :: b(n)
    b(1:n) = acc_solver % matrix % rhs(1:n)
  end subroutine sta_rhs_copy

  subroutine sta_values_set(x, n) bind(c, name="sta_values_set")
    integer(c_int), value :: n
    real(c_double) :: x(n)
    acc_solver % variable % values(1:n) = x(1:n)
  end subroutine sta_values_set

  subroutine sta_pmv_main(n, x, b) bind(c, name="sta_pmv_main")
    integer(c_int), value :: n
    real(c_double) :: x(n), b(n)
    call CRS_MatrixVectorMultiply(acc_solver % matrix, x(1:n), b(1:n))
  end subroutine sta_pmv_main

  subroutine sta_crs_matvec_pp(n, x, b) bind(c, name="sta_crs_matvec_pp")
    integer(c_int), value :: n
    real(c_double) :: x(n), b(n)
    call CRS_MatrixVectorMultiply(acc_pmat, x(1:n), b(1:n))
  end subroutine sta_crs_matvec_pp

  subroutine sta_crs_lusolve(which, n, b) bind(c, name="sta_crs_lusolve")
    integer(c_int), value :: which, n
    real(c_double) :: b(n)
    if (which == 0) then
      call CRS_LUSolve(n, acc_solver % matrix, b(1:n))
    else if (which == 1) then
      call CRS_LUSolve(n, acc_amat, b(1:n))
    else if (which == 2) then
      call CRS_LUSolve(n, acc_mmat, b(1:n))
    else
      call CRS_LUSolve(n, acc_pmat, b(1:n))
    end if
  end subroutine sta_crs_lusolve

  ! Transcription of CRS_MatrixDivMultiply (internal, block path only).
  subroutine sta_crs_divmult(dim, u, v) bind(c, name="sta_crs_divmult")
    integer(c_int), value :: dim
    real(c_double) :: u(*), v(*)
    type(matrix_t), pointer :: a
    integer, pointer :: cols(:), rows(:), diag(:)
    real(kind=dp), pointer :: values(:)
    integer :: i, j, k, n
    real(kind=dp) :: s
    a => acc_solver % matrix
    n = a % numberofrows
    rows => a % rows
    cols => a % cols
    diag => a % diag
    values => a % values
    do i = 1, n
      k = diag(i)
      if (mod(i, dim+1) == 0) then
        s = 0.0_dp
        do j = rows(i), rows(i+1)-1, dim+1
          s = s - u(cols(j)) * values(j)
        end do
        do j = rows(i)+1, rows(i+1)-1, dim+1
          s = s - u(cols(j)) * values(j)
        end do
        if (dim > 2) then
          do j = rows(i)+2, rows(i+1)-1, dim+1
            s = s - u(cols(j)) * values(j)
          end do
        end if
        v(i/(dim+1)) = s
      end if
    end do
  end subroutine sta_crs_divmult

  subroutine sta_solve_lin(which, n, f, y, dofs, norm) &
      bind(c, name="sta_solve_lin")
    integer(c_int), value :: which, n, dofs
    real(c_double) :: f(n), y(n), norm
    type(matrix_t), pointer :: pa
    type(solver_t), pointer :: ps
    if (which == 0) then
      pa => acc_amat
      ps => acc_vsolver
    else if (which == 1) then
      pa => acc_mmat
      ps => acc_jsolver
    else
      pa => acc_pmat
      ps => acc_psolver
    end if
    currentmodel % solver => ps
    call SolveLinearSystem(pa, f(1:n), y(1:n), norm, dofs, ps)
    currentmodel % solver => acc_solver
  end subroutine sta_solve_lin

  subroutine sta_prec_traction(fvec, dvec) bind(c, name="sta_prec_traction")
    real(c_double) :: fvec(*), dvec(*)
    integer :: j, m, i
    j = size(acc_traction)
    do m = 1, j
      i = acc_traction(m)
      fvec(i) = -dvec(i)
    end do
  end subroutine sta_prec_traction

  subroutine sta_precomp_toggle(which, val) bind(c, name="sta_precomp_toggle")
    integer(c_int), value :: which, val
    type(solver_t), pointer :: ps
    character(len=max_name_len) :: str
    logical :: gotit
    if (which == 0) then
      ps => acc_vsolver
    else if (which == 1) then
      ps => acc_jsolver
    else
      ps => acc_psolver
    end if
    str = ListGetString(ps % values, 'Linear System Preconditioning', gotit)
    if (gotit .and. seql(str, 'ilu')) &
        call ListAddLogical(ps % values, 'No Precondition Recompute', val /= 0)
  end subroutine sta_precomp_toggle

  subroutine sta_get_precond(which, buf, buflen, found) &
      bind(c, name="sta_get_precond")
    integer(c_int), value :: which, buflen
    character(kind=c_char) :: buf(*)
    integer(c_int) :: found
    type(solver_t), pointer :: ps
    character(len=max_name_len) :: str
    logical :: gotit
    integer :: i, m
    if (which == 0) then
      ps => acc_vsolver
    else if (which == 1) then
      ps => acc_jsolver
    else
      ps => acc_psolver
    end if
    str = ListGetString(ps % values, 'Linear System Preconditioning', gotit)
    found = merge(1, 0, gotit)
    do i = 1, buflen
      buf(i) = ' '
    end do
    m = min(buflen, len_trim(str))
    do i = 1, m
      buf(i) = str(i:i)
    end do
  end subroutine sta_get_precond

  function sta_get_max_restarts(found) result(r) &
      bind(c, name="sta_get_max_restarts")
    integer(c_int) :: found, r
    logical :: gotit
    r = ListGetInteger(acc_solver % values, 'Max GCR Restarts', gotit)
    if (.not. gotit) r = 1
    found = merge(1, 0, gotit)
  end function sta_get_max_restarts
  function sta_block_rows(which) result(n) bind(c, name="sta_block_rows")
    integer(c_int), value :: which
    integer(c_int) :: n
    if (which == 0) then
      n = acc_amat % numberofrows
    else
      n = acc_mmat % numberofrows
    end if
  end function sta_block_rows

  subroutine sta_outer_nonlin_iter(iter) bind(c, name="sta_outer_nonlin_iter")
    integer(c_int), value :: iter
    print *, 'OuterIteration residuals for nonlinear iteration: ', iter
  end subroutine sta_outer_nonlin_iter

  subroutine sta_crs_matvec_which(which, n, x, b) &
      bind(c, name="sta_crs_matvec_which")
    integer(c_int), value :: which, n
    real(c_double) :: x(n), b(n)
    if (which == 0) then
      call CRS_MatrixVectorMultiply(acc_solver % matrix, x(1:n), b(1:n))
    else if (which == 1) then
      call CRS_MatrixVectorMultiply(acc_amat, x(1:n), b(1:n))
    else if (which == 2) then
      call CRS_MatrixVectorMultiply(acc_mmat, x(1:n), b(1:n))
    else
      call CRS_MatrixVectorMultiply(acc_pmat, x(1:n), b(1:n))
    end if
  end subroutine sta_crs_matvec_which


end module stokes_acc
