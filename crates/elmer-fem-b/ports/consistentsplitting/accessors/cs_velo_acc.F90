! ConsistentSplitting pilot: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/ConsistentSplitting/VelocitySolver.F90.
! Same contract as cs_comp_acc (C-only arguments; SAVE state mirrored in
! Fortran; order-sensitive reductions stay in Fortran).
!
! The rotational-form outflow-boundary block and the serial
! SetBoundaryConditions calls run in wrappers off the captured handles.
! The parallel arms are transcribed and compile-checked (never taken on
! the serial tree).

module cs_velo_acc
  use iso_c_binding
  use defutils
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(model_t), pointer, save :: acc_model => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(variable_t), pointer, save :: acc_flowvar => null()
  type(variable_t), pointer, save :: acc_presvar => null()
  type(valuelist_t), pointer, save :: acc_bodyforce => null()
  type(valuelist_t), pointer, save :: acc_bc => null()
  integer, allocatable, save :: acc_indexes(:)
  ! LocalMatrix context (mirrors its SAVE Nodes and IP local).
  type(nodes_t), target, save :: acc_nodes_lm
  type(gaussintegrationpoints_t), save :: acc_ip_lm
  ! LocalMatrixBoundary context.
  type(nodes_t), target, save :: acc_nodes_lb
  type(gaussintegrationpoints_t), save :: acc_ip_lb
  ! Entry SAVE: timestep tracking.
  integer, save :: acc_current_done_time = 0

contains

  ! -- entry context ---------------------------------------------------
  subroutine csv_begin(h_model, h_solver) bind(c, name="csv_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine csv_begin

  function csv_parallel() result(r) bind(c, name="csv_parallel")
    integer(c_int) :: r
    r = merge(1, 0, parenv % pes > 1)
  end function csv_parallel

  subroutine csv_conv_flags(conv, block, found_c, found_b) &
      bind(c, name="csv_conv_flags")
    integer(c_int) :: conv, block, found_c, found_b
    type(valuelist_t), pointer :: sp
    logical :: v, f
    sp => GetSolverParams()
    v = GetLogical(sp, 'Convective', f)
    conv = merge(1, 0, v)
    found_c = merge(1, 0, f)
    v = GetLogical(sp, 'Block Preconditioning', f)
    block = merge(1, 0, v)
    found_b = merge(1, 0, f)
  end subroutine csv_conv_flags

  function csv_new_timestep() result(r) bind(c, name="csv_new_timestep")
    integer(c_int) :: r
    if (acc_current_done_time /= acc_solver % donetime) then
      acc_current_done_time = acc_current_done_time + 1
      r = 1
    else
      r = 0
    end if
  end function csv_new_timestep

  subroutine csv_const_bulk(cb, cb_in_use, newts) &
      bind(c, name="csv_const_bulk")
    integer(c_int) :: cb, cb_in_use
    integer(c_int), value :: newts
    type(valuelist_t), pointer :: sp
    logical :: v, f
    sp => GetSolverParams()
    v = GetLogical(sp, 'Constant Bulk Matrix', f)
    cb = merge(1, 0, v)
    cb_in_use = merge(1, 0, v .and. &
        associated(acc_solver % matrix % bulkvalues) .and. (newts == 0))
  end subroutine csv_const_bulk

  subroutine csv_ilu_recompute(newts) bind(c, name="csv_ilu_recompute")
    integer(c_int), value :: newts
    character(:), allocatable :: str
    logical :: found
    str = ListGetString(acc_solver % values, 'Linear System Preconditioning', found)
    if (found) then
      if (SEQL(str, 'ilu')) then
        if (newts /= 0) then
          call ListAddLogical(acc_solver % values, 'No Precondition Recompute', .false.)
        else
          call ListAddLogical(acc_solver % values, 'No Precondition Recompute', .true.)
        end if
      end if
    end if
  end subroutine csv_ilu_recompute

  subroutine csv_convection(block, picard, newton, hybrid, form, nliter) &
      bind(c, name="csv_convection")
    integer(c_int), value :: block
    integer(c_int) :: picard, newton, hybrid, form, nliter
    character(:), allocatable :: s
    logical :: found
    real(kind=dp) :: tol
    ! form: 0 = standard, 1 = rotational. Mirrors the source decision tree;
    ! a bad Convection Form is Fatal, exactly as in the source.
    picard = 0
    newton = 0
    hybrid = 0
    form = 0
    nliter = 1
    s = ListGetString(acc_solver % values, 'Convection Form', found)
    if (.not. found) s = 'standard'
    if (s /= 'rotational') then
      if (s == 'standard') then
        picard = 1
      else
        call Fatal('VelocitySolver', 'Convection Form must be rotational or standard')
      end if
    else
      form = 1
      s = ListGetString(acc_solver % values, 'Nonlinear Iteration Method', found)
      if (.not. found) s = 'picard'
      if (s == 'hybrid') then
        newton = 1
        hybrid = 1
      else
        if (s /= 'picard') then
          newton = 1
        else
          picard = 1
        end if
      end if
    end if
    if (block == 0) then
      nliter = ListGetInteger(acc_solver % values, 'Nonlinear System Max Iterations', minv=0)
      tol = ListGetConstReal(acc_solver % values, &
          'Nonlinear System Convergence Tolerance', minv=0.0d0)
    end if
  end subroutine csv_convection

  function csv_dim() result(d) bind(c, name="csv_dim")
    integer(c_int) :: d
    d = CoordinateSystemDimension()
  end function csv_dim

  function csv_mesh_max_dofs() result(n) bind(c, name="csv_mesh_max_dofs")
    integer(c_int) :: n
    n = acc_solver % mesh % maxelementdofs
  end function csv_mesh_max_dofs

  function csv_matrix_rows() result(n) bind(c, name="csv_matrix_rows")
    integer(c_int) :: n
    n = acc_solver % matrix % numberofrows
  end function csv_matrix_rows

  subroutine csv_bind_vars(block, ok_flow, ok_pres) bind(c, name="csv_bind_vars")
    integer(c_int), value :: block
    integer(c_int) :: ok_flow, ok_pres
    type(mesh_t), pointer :: m
    m => GetMesh()
    if (block /= 0) then
      acc_flowvar => VariableGet(m % variables, 'Flow')
      ok_flow = merge(1, 0, associated(acc_flowvar))
      ok_pres = 1
    else
      acc_presvar => VariableGet(m % variables, 'TotalPressure')
      ok_pres = merge(1, 0, associated(acc_presvar))
      ok_flow = 1
    end if
  end subroutine csv_bind_vars

  subroutine csv_fatal_flow() bind(c, name="csv_fatal_flow")
    call Fatal('VelocitySolver', 'The coupled flow variable Flow was not found')
  end subroutine csv_fatal_flow

  subroutine csv_fatal_pres() bind(c, name="csv_fatal_pres")
    call Fatal('VelocitySolver', 'The pressure variable TotalPressure was not found')
  end subroutine csv_fatal_pres

  subroutine csv_fatal_alloc() bind(c, name="csv_fatal_alloc")
    call Fatal('VelocitySolver', 'Memory allocation error.')
  end subroutine csv_fatal_alloc

  subroutine csv_alloc(n) bind(c, name="csv_alloc")
    integer(c_int), value :: n
    if (allocated(acc_indexes)) deallocate(acc_indexes)
    allocate(acc_indexes(n))
  end subroutine csv_alloc

  subroutine csv_default_init() bind(c, name="csv_default_init")
    call DefaultInitialize()
  end subroutine csv_default_init

  subroutine csv_bulk_restore() bind(c, name="csv_bulk_restore")
    acc_solver % matrix % values = acc_solver % matrix % bulkvalues
    acc_solver % matrix % rhs = 0.0_dp
  end subroutine csv_bulk_restore

  ! -- active-element driver -------------------------------------------
  function csv_nof_active() result(n) bind(c, name="csv_nof_active")
    integer(c_int) :: n
    n = GetNOFActive()
  end function csv_nof_active

  subroutine csv_begin_active(t, n, nb, nd) bind(c, name="csv_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n, nb, nd
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes()
    nb = GetElementNOFBDOFs()
    nd = GetElementDOFs(acc_indexes)
  end subroutine csv_begin_active

  subroutine csv_bodyforce(n, load) bind(c, name="csv_bodyforce")
    integer(c_int), value :: n
    real(c_double) :: load(n, 4)
    logical :: found
    integer :: has
    load = 0.0d0
    acc_bodyforce => GetBodyForce()
    has = merge(1, 0, associated(acc_bodyforce))
    if (has /= 0) then
      load(1:n, 1) = GetReal(acc_bodyforce, 'Force 1', found)
      load(1:n, 2) = GetReal(acc_bodyforce, 'Force 2', found)
      load(1:n, 3) = GetReal(acc_bodyforce, 'Force 3', found)
    end if
  end subroutine csv_bodyforce

  subroutine csv_material(n, rho, mu) bind(c, name="csv_material")
    integer(c_int), value :: n
    real(c_double) :: rho(n), mu(n)
    type(valuelist_t), pointer :: mat
    mat => GetMaterial()
    rho(1:n) = GetReal(mat, 'Density')
    mu(1:n) = GetReal(mat, 'Viscosity')
  end subroutine csv_material

  subroutine csv_gather(block, newts, iter, nd, vel, pres) &
      bind(c, name="csv_gather")
    integer(c_int), value :: block, newts, iter, nd
    real(c_double) :: vel(*), pres(*)
    ! Exact transcription of the previous-iterate gather. vel is
    ! (dim,n)-major flat (dim slots per node), pres is (1,n) flat.
    integer :: i, k, d
    d = CoordinateSystemDimension()
    if (block /= 0) then
      do i = 1, d
        do k = 1, nd
          vel((k - 1) * d + i) = acc_flowvar % prevvalues( &
              acc_flowvar % dofs * (acc_flowvar % perm(acc_indexes(k)) - 1) + i, 1)
        end do
      end do
      do k = 1, nd
        pres(k) = acc_flowvar % values( &
            acc_flowvar % dofs * (acc_flowvar % perm(acc_indexes(k)) - 1) + d + 1)
      end do
    else
      if (newts /= 0) then
        do k = 1, nd
          pres(k) = acc_presvar % prevvalues(acc_presvar % perm(acc_indexes(k)), 1)
        end do
        if (iter == 1) then
          do i = 1, d
            do k = 1, nd
              vel((k - 1) * d + i) = acc_solver % variable % prevvalues( &
                  acc_solver % variable % dofs * &
                  (acc_solver % variable % perm(acc_indexes(k)) - 1) + i, 1)
            end do
          end do
        else
          do i = 1, d
            do k = 1, nd
              vel((k - 1) * d + i) = acc_solver % variable % values( &
                  acc_solver % variable % dofs * &
                  (acc_solver % variable % perm(acc_indexes(k)) - 1) + i)
            end do
          end do
        end if
      else
        do k = 1, nd
          pres(k) = acc_presvar % values(acc_presvar % perm(acc_indexes(k)))
        end do
        do i = 1, d
          do k = 1, nd
            vel((k - 1) * d + i) = acc_solver % variable % values( &
                acc_solver % variable % dofs * &
                (acc_solver % variable % perm(acc_indexes(k)) - 1) + i)
          end do
        end do
      end if
    end if
  end subroutine csv_gather

  ! -- LocalMatrix context ----------------------------------------------
  subroutine csv_lm_nodes() bind(c, name="csv_lm_nodes")
    call GetElementNodes(acc_nodes_lm)
  end subroutine csv_lm_nodes

  function csv_lm_ip_begin() result(n) bind(c, name="csv_lm_ip_begin")
    integer(c_int) :: n
    acc_ip_lm = GaussPoints(acc_elem)
    n = acc_ip_lm % n
  end function csv_lm_ip_begin

  subroutine csv_lm_ip_point(t, u, v, w, s) bind(c, name="csv_lm_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip_lm % u(t)
    v = acc_ip_lm % v(t)
    w = acc_ip_lm % w(t)
    s = acc_ip_lm % s(t)
  end subroutine csv_lm_ip_point

  subroutine csv_lm_eleminfo(u, v, w, ntot, basis, dbasis, detj, stat) &
      bind(c, name="csv_lm_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: ntot
    real(c_double) :: basis(ntot), dbasis(ntot, 3), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_elem, acc_nodes_lm, u, v, w, detj, &
        basis(1:ntot), dbasis(1:ntot, 1:3))
    stat = merge(1, 0, ok)
  end subroutine csv_lm_eleminfo

  subroutine csv_lm_point(n, nd, ntot, dim, convect, form, newton, bulk, &
      sdet, load, rho_n, mu_n, vel_n, pres_n, basis, dbasis, stiff, force, mass) &
      bind(c, name="csv_lm_point")
    integer(c_int), value :: n, nd, ntot, dim, convect, form, newton, bulk
    real(c_double), value :: sdet
    real(c_double) :: load(ntot, 4), rho_n(ntot), mu_n(ntot)
    real(c_double) :: vel_n(dim, ntot), pres_n(ntot)
    real(c_double) :: basis(ntot), dbasis(ntot, 3)
    real(c_double), target :: stiff(dim * ntot, dim * ntot)
    real(c_double), target :: force(dim * ntot), mass(dim * ntot, dim * ntot)
    ! Exact transcription of the LocalMatrix per-IP body (reductions, the
    ! p-loop with all branch arms, the Newton force correction). Same
    ! compiler, same evaluation order as the original internal procedure.
    ! The StabilizedMethod arm is dead upstream (forced .FALSE.); it is
    ! omitted here — the dummy is accepted and ignored — and documented
    ! in the Rust caller.
    real(c_double) :: mu, rho, pressure, velo(dim), grad(dim, dim)
    real(c_double) :: loadatip(dim), w1, w2, w3
    integer :: i, j, p, q
    real(c_double), pointer :: a(:, :), m(:, :), f(:)
    mu = sum(basis(1:n) * mu_n(1:n))
    rho = sum(basis(1:n) * rho_n(1:n))
    pressure = sum(pres_n(1:n) * basis(1:n))
    velo = matmul(vel_n(1:dim, 1:nd), basis(1:nd))
    grad = matmul(vel_n(1:dim, 1:nd), dbasis(1:nd, 1:dim))
    if (dim > 2) then
      w1 = sum(vel_n(3, 1:nd) * dbasis(1:nd, 2)) - &
          sum(vel_n(2, 1:nd) * dbasis(1:nd, 3))
      w2 = sum(vel_n(1, 1:nd) * dbasis(1:nd, 3)) - &
          sum(vel_n(3, 1:nd) * dbasis(1:nd, 1))
      w3 = sum(vel_n(2, 1:nd) * dbasis(1:nd, 1)) - &
          sum(vel_n(1, 1:nd) * dbasis(1:nd, 2))
    else
      w3 = sum(vel_n(2, 1:nd) * dbasis(1:nd, 1)) - &
          sum(vel_n(1, 1:nd) * dbasis(1:nd, 2))
    end if
    loadatip = matmul(basis(1:n), load(1:n, 1:dim))
    do p = 1, ntot
      if (bulk == 0) then
        do q = 1, ntot
          i = dim * (p - 1)
          j = dim * (q - 1)
          a => stiff(i + 1:i + dim, j + 1:j + dim)
          m => mass(i + 1:i + dim, j + 1:j + dim)
          do i = 1, dim
            m(i, i) = m(i, i) + sdet * rho * basis(q) * basis(p)
            do j = 1, dim
              a(i, i) = a(i, i) + sdet * mu * dbasis(q, j) * dbasis(p, j)
              a(i, j) = a(i, j) + sdet * mu * dbasis(q, i) * dbasis(p, j)
            end do
          end do
          if (convect /= 0) then
            if (form == 0) then
              call csv_convect_standard(dim, sdet, rho, velo, basis, &
                  dbasis, p, q, stiff, dim * ntot, ntot)
            else
              call csv_convect_rotational(dim, sdet, rho, velo, w1, w2, w3, &
                  newton, basis, dbasis, p, q, stiff, dim * ntot, ntot)
            end if
          end if
        end do
      else
        do q = 1, ntot
          i = dim * (p - 1)
          j = dim * (q - 1)
          m => mass(i + 1:i + dim, j + 1:j + dim)
          do i = 1, dim
            m(i, i) = m(i, i) + sdet * rho * basis(q) * basis(p)
          end do
        end do
      end if
      i = dim * (p - 1)
      f => force(i + 1:i + dim)
      f(1:dim) = f(1:dim) + &
          sdet * (loadatip(1:dim) * basis(p) + pressure * dbasis(p, 1:dim))
      if ((convect /= 0) .and. (newton /= 0)) then
        f(1) = f(1) - rho * sdet * basis(p) * velo(2) * w3
        f(2) = f(2) + rho * sdet * basis(p) * velo(1) * w3
        if (dim > 2) then
          f(1) = f(1) + rho * sdet * basis(p) * velo(3) * w2
          f(2) = f(2) - rho * sdet * basis(p) * velo(3) * w1
          f(3) = f(3) - rho * sdet * basis(p) * velo(1) * w2
          f(3) = f(3) + rho * sdet * basis(p) * velo(2) * w1
        end if
      end if
    end do
  end subroutine csv_lm_point

  subroutine csv_convect_standard(dim, sdet, rho, velo, basis, dbasis, &
      p, q, stiff, ld, ldb) bind(c)
    integer(c_int), value :: dim, p, q, ld, ldb
    real(c_double), value :: sdet, rho
    real(c_double) :: velo(dim), basis(*), dbasis(ldb, 3), stiff(ld, ld)
    integer :: r1, c1
    r1 = dim * (p - 1) + 1
    c1 = dim * (q - 1) + 1
    stiff(r1, c1) = stiff(r1, c1) + rho * sdet * velo(2) * dbasis(q, 2) * basis(p)
    stiff(r1, c1) = stiff(r1, c1) + rho * sdet * velo(1) * dbasis(q, 1) * basis(p)
    stiff(r1 + 1, c1 + 1) = stiff(r1 + 1, c1 + 1) + &
        rho * sdet * velo(1) * dbasis(q, 1) * basis(p)
    stiff(r1 + 1, c1 + 1) = stiff(r1 + 1, c1 + 1) + &
        rho * sdet * velo(2) * dbasis(q, 2) * basis(p)
    if (dim > 2) then
      stiff(r1, c1) = stiff(r1, c1) + rho * sdet * velo(3) * dbasis(q, 3) * basis(p)
      stiff(r1 + 1, c1 + 1) = stiff(r1 + 1, c1 + 1) + &
          rho * sdet * velo(3) * dbasis(q, 3) * basis(p)
      stiff(r1 + 2, c1 + 2) = stiff(r1 + 2, c1 + 2) + &
          rho * sdet * velo(1) * dbasis(q, 1) * basis(p)
      stiff(r1 + 2, c1 + 2) = stiff(r1 + 2, c1 + 2) + &
          rho * sdet * velo(2) * dbasis(q, 2) * basis(p)
      stiff(r1 + 2, c1 + 2) = stiff(r1 + 2, c1 + 2) + &
          rho * sdet * velo(3) * dbasis(q, 3) * basis(p)
    end if
  end subroutine csv_convect_standard

  subroutine csv_convect_rotational(dim, sdet, rho, velo, w1, w2, w3, &
      newton, basis, dbasis, p, q, stiff, ld, ldb) bind(c)
    integer(c_int), value :: dim, newton, p, q, ld, ldb
    real(c_double), value :: sdet, rho, w1, w2, w3
    real(c_double) :: velo(dim), basis(*), dbasis(ldb, 3), stiff(ld, ld)
    integer :: r1, c1
    r1 = dim * (p - 1) + 1
    c1 = dim * (q - 1) + 1
    if (dim > 2) then
      stiff(r1, c1) = stiff(r1, c1) + rho * sdet * velo(2) * dbasis(q, 2) * basis(p) + &
          rho * sdet * velo(3) * dbasis(q, 3) * basis(p)
      stiff(r1, c1 + 1) = stiff(r1, c1 + 1) - &
          rho * sdet * velo(2) * dbasis(q, 1) * basis(p)
      stiff(r1, c1 + 2) = stiff(r1, c1 + 2) - &
          rho * sdet * velo(3) * dbasis(q, 1) * basis(p)
      stiff(r1 + 1, c1) = stiff(r1 + 1, c1) - &
          rho * sdet * velo(1) * dbasis(q, 2) * basis(p)
      stiff(r1 + 1, c1 + 1) = stiff(r1 + 1, c1 + 1) + &
          rho * sdet * velo(1) * dbasis(q, 1) * basis(p) + &
          rho * sdet * velo(3) * dbasis(q, 3) * basis(p)
      stiff(r1 + 1, c1 + 2) = stiff(r1 + 1, c1 + 2) - &
          rho * sdet * velo(3) * dbasis(q, 2) * basis(p)
      stiff(r1 + 2, c1) = stiff(r1 + 2, c1) - &
          rho * sdet * velo(1) * dbasis(q, 3) * basis(p)
      stiff(r1 + 2, c1 + 1) = stiff(r1 + 2, c1 + 1) - &
          rho * sdet * velo(2) * dbasis(q, 3) * basis(p)
      stiff(r1 + 2, c1 + 2) = stiff(r1 + 2, c1 + 2) + &
          rho * sdet * velo(2) * dbasis(q, 2) * basis(p) + &
          rho * sdet * velo(1) * dbasis(q, 1) * basis(p)
      if (newton /= 0) then
        stiff(r1, c1 + 1) = stiff(r1, c1 + 1) - &
            rho * sdet * w3 * basis(q) * basis(p)
        stiff(r1, c1 + 2) = stiff(r1, c1 + 2) + &
            rho * sdet * w2 * basis(q) * basis(p)
        stiff(r1 + 1, c1) = stiff(r1 + 1, c1) + &
            rho * sdet * w3 * basis(q) * basis(p)
        stiff(r1 + 1, c1 + 2) = stiff(r1 + 1, c1 + 2) - &
            rho * sdet * w1 * basis(q) * basis(p)
        stiff(r1 + 2, c1) = stiff(r1 + 2, c1) - &
            rho * sdet * w2 * basis(q) * basis(p)
        stiff(r1 + 2, c1 + 1) = stiff(r1 + 2, c1 + 1) + &
            rho * sdet * w1 * basis(q) * basis(p)
      end if
    else
      stiff(r1, c1) = stiff(r1, c1) + &
          rho * sdet * velo(2) * dbasis(q, 2) * basis(p)
      stiff(r1, c1 + 1) = stiff(r1, c1 + 1) - &
          rho * sdet * velo(2) * dbasis(q, 1) * basis(p)
      stiff(r1 + 1, c1) = stiff(r1 + 1, c1) - &
          rho * sdet * velo(1) * dbasis(q, 2) * basis(p)
      stiff(r1 + 1, c1 + 1) = stiff(r1 + 1, c1 + 1) + &
          rho * sdet * velo(1) * dbasis(q, 1) * basis(p)
      if (newton /= 0) then
        stiff(r1, c1 + 1) = stiff(r1, c1 + 1) - &
            rho * sdet * w3 * basis(q) * basis(p)
        stiff(r1 + 1, c1) = stiff(r1 + 1, c1) + &
            rho * sdet * w3 * basis(q) * basis(p)
      end if
    end if
  end subroutine csv_convect_rotational

  ! -- shared assembly tail ---------------------------------------------
  subroutine csv_default_time(m_ld, mass, stiff, force) &
      bind(c, name="csv_default_time")
    integer(c_int), value :: m_ld
    real(c_double) :: mass(m_ld, m_ld), stiff(m_ld, m_ld), force(m_ld)
    call Default1stOrderTime(mass(1:m_ld, 1:m_ld), stiff(1:m_ld, 1:m_ld), &
        force(1:m_ld))
  end subroutine csv_default_time

  subroutine csv_update_eq(ld, g, f) bind(c, name="csv_update_eq")
    integer(c_int), value :: ld
    real(c_double) :: g(ld, ld), f(ld)
    call DefaultUpdateEquations(g(1:ld, 1:ld), f(1:ld))
  end subroutine csv_update_eq

  subroutine csv_update_force(ld, f) bind(c, name="csv_update_force")
    integer(c_int), value :: ld
    real(c_double) :: f(ld)
    call DefaultUpdateForce(f(1:ld))
  end subroutine csv_update_force

  subroutine csv_finish_bulk_cond(in_use) bind(c, name="csv_finish_bulk_cond")
    integer(c_int), value :: in_use
    if (in_use == 0) call DefaultFinishBulkAssembly()
  end subroutine csv_finish_bulk_cond

  subroutine csv_finish() bind(c, name="csv_finish")
    call DefaultFinishAssembly()
  end subroutine csv_finish

  subroutine csv_dirichlet() bind(c, name="csv_dirichlet")
    call DefaultDirichletBCs()
  end subroutine csv_dirichlet

  function csv_solve() result(norm) bind(c, name="csv_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function csv_solve

  ! -- rotational-form outflow boundary ----------------------------------
  function csv_nof_boundary() result(n) bind(c, name="csv_nof_boundary")
    integer(c_int) :: n
    n = GetNOFBoundaryElements()
  end function csv_nof_boundary

  subroutine csv_begin_boundary(t, n, nd, active, family) &
      bind(c, name="csv_begin_boundary")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, active, family
    acc_elem => GetBoundaryElement(t)
    active = merge(1, 0, ActiveBoundaryElement())
    n = GetElementNOFNodes()
    nd = GetElementDOFs(acc_indexes)
    family = GetElementFamily()
  end subroutine csv_begin_boundary

  subroutine csv_bc() bind(c, name="csv_bc")
    acc_bc => GetBC()
  end subroutine csv_bc

  function csv_bc_associated() result(r) bind(c, name="csv_bc_associated")
    integer(c_int) :: r
    r = merge(1, 0, associated(acc_bc))
  end function csv_bc_associated

  subroutine csv_bc_outflow(val, found) bind(c, name="csv_bc_outflow")
    integer(c_int) :: val, found
    logical :: v, f
    v = ListGetLogical(acc_bc, 'Outflow Boundary', f)
    val = merge(1, 0, v)
    found = merge(1, 0, f)
  end subroutine csv_bc_outflow

  subroutine csv_gather_outflow(newts, iter, nd, vel) &
      bind(c, name="csv_gather_outflow")
    integer(c_int), value :: newts, iter, nd
    real(c_double) :: vel(*)
    integer :: i, k, d
    d = CoordinateSystemDimension()
    if ((newts /= 0) .and. (iter == 1)) then
      do i = 1, d
        do k = 1, nd
          vel((k - 1) * d + i) = acc_solver % variable % prevvalues( &
              acc_solver % variable % dofs * &
              (acc_solver % variable % perm(acc_indexes(k)) - 1) + i, 1)
        end do
      end do
    else
      do i = 1, d
        do k = 1, nd
          vel((k - 1) * d + i) = acc_solver % variable % values( &
              acc_solver % variable % dofs * &
              (acc_solver % variable % perm(acc_indexes(k)) - 1) + i)
        end do
      end do
    end if
  end subroutine csv_gather_outflow

  subroutine csv_lb_nodes() bind(c, name="csv_lb_nodes")
    call GetElementNodes(acc_nodes_lb)
  end subroutine csv_lb_nodes

  function csv_lb_ip_begin() result(n) bind(c, name="csv_lb_ip_begin")
    integer(c_int) :: n
    acc_ip_lb = GaussPoints(acc_elem)
    n = acc_ip_lb % n
  end function csv_lb_ip_begin

  subroutine csv_lb_ip_point(t, u, v, w, s) bind(c, name="csv_lb_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip_lb % u(t)
    v = acc_ip_lb % v(t)
    w = acc_ip_lb % w(t)
    s = acc_ip_lb % s(t)
  end subroutine csv_lb_ip_point

  subroutine csv_lb_eleminfo(u, v, w, nd, basis, dbasis, detj, stat) &
      bind(c, name="csv_lb_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: nd
    real(c_double) :: basis(nd), dbasis(nd, 3), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_elem, acc_nodes_lb, u, v, w, detj, &
        basis(1:nd), dbasis(1:nd, 1:3))
    stat = merge(1, 0, ok)
  end subroutine csv_lb_eleminfo

  subroutine csv_lb_point(nd, dim, newton, rho, sdet, vel_n, basis, &
      dbasis, normal, stiff, force) bind(c, name="csv_lb_point")
    integer(c_int), value :: nd, dim, newton
    real(c_double), value :: rho, sdet
    real(c_double) :: vel_n(dim, nd), basis(nd), dbasis(nd, 3)
    real(c_double) :: normal(dim)
    real(c_double) :: stiff(dim * nd, dim * nd), force(dim * nd)
    ! Exact transcription of the LocalMatrixBoundary per-IP body (the
    ! commented-out Picard block stays commented, as upstream).
    real(c_double) :: velo(dim), squaredvelo, c
    integer :: p, q, i
    velo = matmul(vel_n(1:dim, 1:nd), basis(1:nd))
    squaredvelo = rho * sum(velo(1:dim) * velo(1:dim))
    if (newton /= 0) then
      c = 1.0d0 * rho
    else
      c = 0.5d0 * rho
    end if
    do p = 1, nd
      do i = 1, dim
        if (newton /= 0) &
            force(dim * (p - 1) + i) = force(dim * (p - 1) + i) + &
            sdet * normal(i) * basis(p) * 0.5d0 * squaredvelo
        do q = 1, nd
          stiff(dim * (p - 1) + i, dim * (q - 1) + 1) = &
              stiff(dim * (p - 1) + i, dim * (q - 1) + 1) + &
              sdet * c * basis(p) * normal(i) * velo(1) * basis(q)
          stiff(dim * (p - 1) + i, dim * (q - 1) + 2) = &
              stiff(dim * (p - 1) + i, dim * (q - 1) + 2) + &
              sdet * c * basis(p) * normal(i) * velo(2) * basis(q)
          if (dim > 2) &
              stiff(dim * (p - 1) + i, dim * (q - 1) + 3) = &
              stiff(dim * (p - 1) + i, dim * (q - 1) + 3) + &
              sdet * c * basis(p) * normal(i) * velo(3) * basis(q)
        end do
      end do
    end do
  end subroutine csv_lb_point

  subroutine csv_lb_normal(u, v, normal) bind(c, name="csv_lb_normal")
    real(c_double), value :: u, v
    real(c_double) :: normal(*)
    real(c_double) :: n3(3)
    integer :: d
    d = CoordinateSystemDimension()
    n3 = Normalvector(acc_elem, acc_nodes_lb, u, v, .true.)
    normal(1:d) = n3(1:d)
  end subroutine csv_lb_normal

  ! -- serial Dirichlet handling (nonstandard, see source comment) -------
  subroutine csv_set_bc(dof, ndofs) bind(c, name="csv_set_bc")
    integer(c_int), value :: dof, ndofs
    character(len=64) :: compname
    compname = ComponentName(acc_solver % variable % name, dof)
    call csv_set_boundary_conditions(acc_model, acc_solver % matrix, compname, &
        dof, ndofs, acc_solver % variable % perm, acc_solver % matrix % rhs)
  end subroutine csv_set_bc

  subroutine csv_set_boundary_conditions(model, stiffmatrix, name, dof, &
      ndofs, perm, rhs)
    type(model_t) :: model
    type(matrix_t), pointer :: stiffmatrix
    character(len=*) :: name
    integer :: dof, ndofs, perm(:)
    real(kind=dp), optional :: rhs(:)
    ! Verbatim transcription of the entry's internal SetBoundaryConditions
    ! (host-associated Indexes shared via the module SAVE).
    type(element_t), pointer :: currentelement
    integer, pointer :: nodeindexes(:)
    integer :: i, j, k, n, t, nd
    logical :: gotit
    real(kind=dp) :: work(model % maxelementnodes)
    do t = model % numberofbulkelements + 1, &
        model % numberofbulkelements + model % numberofboundaryelements
      currentelement => model % elements(t)
      model % currentelement => model % elements(t)
      n = GetElementNOFNodes()
      nd = GetElementDOFs(acc_indexes)
      nodeindexes => currentelement % nodeindexes(1:n)
      do i = 1, model % numberofbcs
        if (currentelement % boundaryinfo % constraint == &
            model % bcs(i) % tag) then
          work(1:n) = ListGetReal(model % bcs(i) % values, &
              name, n, nodeindexes, gotit)
          if (gotit) then
            do j = 1, n
              k = perm(nodeindexes(j))
              if (k > 0) then
                k = ndofs * (k - 1) + dof
                call ZeroRow(stiffmatrix, k)
                call SetMatrixElement(stiffmatrix, k, k, 1.0d0)
                if (present(rhs)) rhs(k) = work(j)
              end if
            end do
            do j = n + 1, nd
              k = perm(acc_indexes(j))
              k = ndofs * (k - 1) + dof
              call ZeroRow(stiffmatrix, k)
              call SetMatrixElement(stiffmatrix, k, k, 1.0d0)
              if (present(rhs)) rhs(k) = 0.0d0
            end do
          end if
        end if
      end do
    end do
  end subroutine csv_set_boundary_conditions

end module cs_velo_acc
