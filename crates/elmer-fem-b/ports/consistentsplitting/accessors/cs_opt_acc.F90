! ConsistentSplitting pilot: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/ConsistentSplitting/OptimalSolutionUpdate.F90.
! Same contract as cs_comp_acc (C-only arguments; SAVE state mirrored in
! Fortran; order-sensitive reductions stay in Fortran).
!
! Notes: NonlinearIterationMethod is overwritten to 'picard' upstream
! (lines 133-134), so Newton/Hybrid are always false — the port hardcodes
! the resulting flags and transcribes every arm. GradDivParam is read and
! passed but never used by any live statement; the dummy is accepted and
! ignored. Stabilization is forced .FALSE. upstream. The IF(.FALSE.)
! outflow block is dead and omitted. The serial path (Orthogonalize) is
! live; the parallel path (GCRUpdate et al.) is transcribed and
! compile-checked.

module cs_opt_acc
  use iso_c_binding
  use defutils
  use parallelutils
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(model_t), pointer, save :: acc_model => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(variable_t), pointer, save :: acc_flowsol => null()
  type(valuelist_t), pointer, save :: acc_bodyforce => null()
  integer, allocatable, save :: acc_indexes(:)
  ! LocalMatrix context (mirrors its SAVE Nodes and IP local).
  type(nodes_t), target, save :: acc_nodes_lm
  type(gaussintegrationpoints_t), save :: acc_ip_lm
  ! Entry SAVE: timestep/round tracking.
  integer, save :: acc_current_done_time = 0
  integer, save :: acc_round = 0
  character(len=512), save :: acc_message

contains

  ! -- entry context ---------------------------------------------------
  subroutine cso_begin(h_model, h_solver) bind(c, name="cso_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine cso_begin

  ! -- assembly driver ----------------------------------------------------
  function cso_num_active_elems() result(n) bind(c, name="cso_num_active_elems")
    integer(c_int) :: n
    n = acc_solver % numberofactiveelements
  end function cso_num_active_elems

  function cso_nof_active() result(n) bind(c, name="cso_nof_active")
    integer(c_int) :: n
    n = GetNOFActive()
  end function cso_nof_active

  subroutine cso_dirichlet() bind(c, name="cso_dirichlet")
    call DefaultDirichletBCs()
  end subroutine cso_dirichlet

  subroutine cso_new_timestep(newts) bind(c, name="cso_new_timestep")
    integer(c_int) :: newts
    if (acc_current_done_time /= acc_solver % donetime) then
      acc_current_done_time = acc_current_done_time + 1
      newts = 1
      acc_round = 0
    else
      newts = 0
    end if
    acc_round = acc_round + 1
  end subroutine cso_new_timestep

  function cso_round() result(r) bind(c, name="cso_round")
    integer(c_int) :: r
    r = acc_round
  end function cso_round

  subroutine cso_const_system(cb_in_use) bind(c, name="cso_const_system")
    integer(c_int) :: cb_in_use
    type(valuelist_t), pointer :: sp
    logical :: v, f
    sp => GetSolverParams()
    v = GetLogical(sp, 'Constant System', f)
    cb_in_use = merge(1, 0, v .and. &
        associated(acc_solver % matrix % bulkvalues) .and. &
        (acc_current_done_time == acc_solver % donetime) .and. &
        associated(acc_solver % matrix % bulkrhs))
  end subroutine cso_const_system

  function cso_dim() result(d) bind(c, name="cso_dim")
    integer(c_int) :: d
    d = CoordinateSystemDimension()
  end function cso_dim

  subroutine cso_graddiv(gd_param) bind(c, name="cso_graddiv")
    real(c_double) :: gd_param
    logical :: gotit, stab
    real(kind=dp) :: p
    stab = ListGetLogical(acc_solver % values, 'Grad-Div Stabilization', gotit)
    if (gotit .and. stab) then
      p = ListGetConstReal(acc_solver % values, &
          'Grad-Div Stabilization Parameter', gotit)
      if (.not. gotit) p = 0.0d0
    else
      p = 0.0d0
    end if
    gd_param = p
  end subroutine cso_graddiv

  subroutine cso_block_precond(val, found) bind(c, name="cso_block_precond")
    integer(c_int) :: val, found
    logical :: v, f
    v = ListGetLogical(acc_solver % values, 'Block Preconditioning', f)
    if (.not. f) v = .true.
    val = merge(1, 0, v)
    found = merge(1, 0, f)
  end subroutine cso_block_precond

  function cso_max_iter() result(n) bind(c, name="cso_max_iter")
    integer(c_int) :: n
    n = ListGetInteger(acc_model % simulation, 'Steady State Max Iterations')
  end function cso_max_iter

  function cso_mesh_max_dofs() result(n) bind(c, name="cso_mesh_max_dofs")
    integer(c_int) :: n
    n = acc_solver % mesh % maxelementdofs
  end function cso_mesh_max_dofs

  function cso_matrix_rows() result(n) bind(c, name="cso_matrix_rows")
    integer(c_int) :: n
    n = acc_solver % matrix % numberofrows
  end function cso_matrix_rows

  subroutine cso_conv_flags(conv, found) bind(c, name="cso_conv_flags")
    integer(c_int) :: conv, found
    type(valuelist_t), pointer :: sp
    logical :: v, f
    sp => GetSolverParams()
    v = GetLogical(sp, 'Convective', f)
    conv = merge(1, 0, v)
    found = merge(1, 0, f)
  end subroutine cso_conv_flags

  subroutine cso_alloc(n) bind(c, name="cso_alloc")
    integer(c_int), value :: n
    if (allocated(acc_indexes)) deallocate(acc_indexes)
    allocate(acc_indexes(n))
  end subroutine cso_alloc

  subroutine cso_fatal_steady() bind(c, name="cso_fatal_steady")
    call Fatal('NavierStokesSolver', 'Steady state problems cannot be handled yet')
  end subroutine cso_fatal_steady

  subroutine cso_fatal_velo() bind(c, name="cso_fatal_velo")
    call Fatal('CoupledNSUpdate', 'Could not find velocity solution VelocityTot')
  end subroutine cso_fatal_velo

  subroutine cso_fatal_pres() bind(c, name="cso_fatal_pres")
    call Fatal('CoupledNSUpdate', 'Could not pressure solution TotalPressure')
  end subroutine cso_fatal_pres

  subroutine cso_fatal_alloc() bind(c, name="cso_fatal_alloc")
    call Fatal('NavierStokesSolver', 'Memory allocation error.')
  end subroutine cso_fatal_alloc

  subroutine cso_default_init() bind(c, name="cso_default_init")
    call DefaultInitialize()
  end subroutine cso_default_init

  subroutine cso_bulk_restore() bind(c, name="cso_bulk_restore")
    acc_solver % matrix % values = acc_solver % matrix % bulkvalues
    acc_solver % matrix % rhs = acc_solver % matrix % bulkrhs
  end subroutine cso_bulk_restore

  function cso_progress_due(at0) result(r) bind(c, name="cso_progress_due")
    real(c_double), value :: at0
    integer(c_int) :: r
    r = merge(1, 0, RealTime() - at0 > 1.0)
  end function cso_progress_due

  subroutine cso_progress(nactive, t, at0) bind(c, name="cso_progress")
    integer(c_int), value :: nactive, t
    real(c_double) :: at0
    character(len=512) :: msg
    write (msg, '(a,i3,a)') '   Assembly: ', &
        int(100.0 - 100.0 * (nactive - t) / nactive), ' % done'
    call Info('NavierStokesSolver', msg, level=5)
    at0 = RealTime()
  end subroutine cso_progress

  function cso_realtime() result(r) bind(c, name="cso_realtime")
    real(c_double) :: r
    r = RealTime()
  end function cso_realtime

  ! -- active-element driver -------------------------------------------
  subroutine cso_begin_active(t, n, nd, nb) bind(c, name="cso_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, nb
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes()
    nd = GetElementDOFs(acc_indexes)
    nb = GetElementNOFBDOFs()
  end subroutine cso_begin_active

  subroutine cso_bodyforce(n, load) bind(c, name="cso_bodyforce")
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
  end subroutine cso_bodyforce

  subroutine cso_material(n, rho, mu) bind(c, name="cso_material")
    integer(c_int), value :: n
    real(c_double) :: rho(n), mu(n)
    type(valuelist_t), pointer :: mat
    mat => GetMaterial()
    rho(1:n) = GetReal(mat, 'Density')
    mu(1:n) = GetReal(mat, 'Viscosity')
  end subroutine cso_material

  subroutine cso_lm_point(n, nd, dim, stab, convect, newton, hybrid, &
      sdet, load, rho_n, mu_n, vel_n, basis, dbasis, stiff, force, mass) &
      bind(c, name="cso_lm_point")
    integer(c_int), value :: n, nd, dim, stab, convect, newton, hybrid
    real(c_double), value :: sdet
    real(c_double) :: load(n, 4), rho_n(n), mu_n(n), vel_n(dim + 1, n)
    real(c_double) :: basis(nd), dbasis(nd, 3)
    real(c_double), target :: stiff((dim + 1) * n, (dim + 1) * n)
    real(c_double), target :: force((dim + 1) * n)
    real(c_double), target :: mass((dim + 1) * n, (dim + 1) * n)
    ! Exact transcription of the LocalMatrix per-IP body. Dead upstream
    ! values (AK/rotterm/StabTerms accumulation, all read only inside
    ! IF(.FALSE.) or commented code) are omitted; every live statement
    ! keeps its order. Newton/Hybrid/Stabilization arms are transcribed
    ! and selected by the flags (all false on the picard path).
    real(c_double) :: mu, rho, velo(dim), grad(dim, dim)
    real(c_double) :: loadatip(dim + 1), w1, w2, w3
    real(c_double), pointer :: a(:, :), m(:, :), f(:)
    integer :: i, j, p, q
    mu = sum(basis(1:n) * mu_n(1:n))
    rho = sum(basis(1:n) * rho_n(1:n))
    velo = matmul(vel_n(1:dim, 1:nd), basis(1:nd))
    grad = matmul(vel_n(1:dim, 1:nd), dbasis(1:nd, 1:dim))
    w1 = 0.0d0
    w2 = 0.0d0
    w3 = 0.0d0
    loadatip(1:dim + 1) = matmul(basis(1:n), load(1:n, 1:dim + 1))
    do p = 1, nd
      do q = 1, nd
        i = (dim + 1) * (p - 1) + 1
        j = (dim + 1) * (q - 1) + 1
        a => stiff(i:i + dim, j:j + dim)
        m => mass(i:i + dim, j:j + dim)
        do i = 1, dim
          do j = 1, dim
            a(i, i) = a(i, i) + sdet * mu * dbasis(q, j) * dbasis(p, j)
            a(i, j) = a(i, j) + sdet * mu * dbasis(q, i) * dbasis(p, j)
          end do
          m(i, i) = m(i, i) + sdet * rho * basis(p) * basis(q)
          if (stab /= 0) then
            a(i, dim + 1) = a(i, dim + 1) - sdet * basis(q) * dbasis(p, i)
          else
            if (q <= n) &
                a(i, dim + 1) = a(i, dim + 1) - sdet * basis(q) * dbasis(p, i)
          end if
          if (p <= n) &
              a(dim + 1, i) = a(dim + 1, i) - sdet * dbasis(q, i) * basis(p)
        end do
        if (convect /= 0) then
          call cso_convect_standard(dim, sdet, rho, velo, newton, hybrid, &
              grad, basis, dbasis, p, q, n, stiff, (dim + 1) * n, nd)
        end if
      end do
      i = (dim + 1) * (p - 1) + 1
      f => force(i:i + dim)
      f = f + sdet * loadatip * basis(p)
      if (newton /= 0) then
        f(1) = f(1) - sdet * rho * basis(p) * velo(2) * w3
        f(2) = f(2) + sdet * rho * basis(p) * velo(1) * w3
        if ((hybrid /= 0) .and. (dim == 2)) then
          f(1) = f(1) - sdet * rho * basis(p) * &
              (grad(1, 1) * velo(1) + grad(1, 2) * velo(2))
          f(2) = f(2) - sdet * rho * basis(p) * &
              (grad(2, 1) * velo(1) + grad(2, 2) * velo(2))
        end if
        if (dim > 2) then
          f(1) = f(1) + sdet * rho * basis(p) * velo(3) * w2
          f(2) = f(2) - sdet * rho * basis(p) * velo(3) * w1
          f(3) = f(3) - sdet * rho * basis(p) * velo(1) * w2
          f(3) = f(3) + sdet * rho * basis(p) * velo(2) * w1
        end if
      end if
    end do
    do p = n + 1, nd
      i = (dim + 1) * p
      force(i) = 0.0d0
      stiff(i, :) = 0.0d0
      if (stab == 0) stiff(i, i) = 1.0d0
    end do
  end subroutine cso_lm_point

  subroutine cso_convect_standard(dim, sdet, rho, velo, newton, hybrid, &
      grad, basis, dbasis, p, q, n, stiff, ld, ldb) bind(c)
    ! The IF(.TRUE.) standard-form arm (the rotation ELSE is dead: the
    ! .TRUE. literal always selects this path — transcribed as taken).
    integer(c_int), value :: dim, newton, hybrid, p, q, n, ld, ldb
    real(c_double), value :: sdet, rho
    real(c_double) :: velo(dim), grad(dim, dim), basis(*), dbasis(ldb, 3)
    real(c_double) :: stiff(ld, ld)
    integer :: r1, c1
    r1 = (dim + 1) * (p - 1) + 1
    c1 = (dim + 1) * (q - 1) + 1
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
  end subroutine cso_convect_standard

  ! -- shared assembly tail ---------------------------------------------
  subroutine cso_default_time(m_ld, mass, stiff, force) &
      bind(c, name="cso_default_time")
    integer(c_int), value :: m_ld
    real(c_double) :: mass(m_ld, m_ld), stiff(m_ld, m_ld), force(m_ld)
    call Default1stOrderTime(mass(1:m_ld, 1:m_ld), stiff(1:m_ld, 1:m_ld), &
        force(1:m_ld))
  end subroutine cso_default_time

  subroutine cso_update_eq(ld, g, f) bind(c, name="cso_update_eq")
    integer(c_int), value :: ld
    real(c_double) :: g(ld, ld), f(ld)
    call DefaultUpdateEquations(g(1:ld, 1:ld), f(1:ld))
  end subroutine cso_update_eq

  subroutine cso_finish_bulk() bind(c, name="cso_finish_bulk")
    call DefaultFinishBulkAssembly()
  end subroutine cso_finish_bulk

  subroutine cso_finish() bind(c, name="cso_finish")
    call DefaultFinishAssembly()
  end subroutine cso_finish

  function cso_parallel() result(r) bind(c, name="cso_parallel")
    integer(c_int) :: r
    r = merge(1, 0, parenv % pes > 1)
  end function cso_parallel

  subroutine cso_info_assembly_done() bind(c, name="cso_info_assembly_done")
    call Info('CoupledNSUpdate', 'Assembly done', level=4)
  end subroutine cso_info_assembly_done

  ! -- serial Dirichlet handling ------------------------------------------
  subroutine cso_set_bc(dof, ndofs) bind(c, name="cso_set_bc")
    integer(c_int), value :: dof, ndofs
    character(len=64) :: compname
    compname = ComponentName(acc_solver % variable % name, dof)
    call cso_set_boundary_conditions(acc_model, acc_solver % matrix, compname, &
        dof, ndofs, acc_solver % variable % perm, acc_solver % matrix % rhs)
  end subroutine cso_set_bc

  subroutine cso_set_boundary_conditions(model, stiffmatrix, name, dof, &
      ndofs, perm, rhs)
    type(model_t) :: model
    type(matrix_t), pointer :: stiffmatrix
    character(len=*) :: name
    integer :: dof, ndofs, perm(:)
    real(kind=dp), optional :: rhs(:)
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
  end subroutine cso_set_boundary_conditions

  ! -- search-direction update --------------------------------------------
  subroutine cso_bind_flow(ok) bind(c, name="cso_bind_flow")
    integer(c_int) :: ok
    acc_flowsol => VariableGet(acc_model % variables, 'VelocityTot')
    ok = merge(1, 0, associated(acc_flowsol))
  end subroutine cso_bind_flow

  subroutine cso_gather_snew(active, nrows, snew) bind(c, name="cso_gather_snew")
    integer(c_int), value :: active, nrows
    real(c_double) :: snew(nrows)
    integer :: t, i, j, k, nd, d, p
    d = CoordinateSystemDimension()
    do t = 1, active
      acc_elem => GetActiveElement(t)
      nd = GetElementDOFs(acc_indexes)
      do i = 1, nd
        j = acc_solver % variable % perm(acc_indexes(i))
        k = acc_flowsol % perm(acc_indexes(i))
        do p = 1, d
          snew((d + 1) * (j - 1) + p) = &
              acc_flowsol % values(d * (k - 1) + p)
        end do
      end do
    end do
    acc_flowsol => VariableGet(acc_model % variables, 'TotalPressure')
    if (.not. associated(acc_flowsol)) &
        call Fatal('CoupledNSUpdate', 'Could not pressure solution TotalPressure')
    do t = 1, active
      acc_elem => GetActiveElement(t)
      nd = GetElementDOFs(acc_indexes)
      do i = 1, nd
        j = acc_solver % variable % perm(acc_indexes(i))
        k = acc_flowsol % perm(acc_indexes(i))
        snew((d + 1) * j) = acc_flowsol % values(k)
      end do
    end do
  end subroutine cso_gather_snew

  subroutine cso_orthogonalize(n, maxrounds, round, snew, s, v, r) &
      bind(c, name="cso_orthogonalize")
    integer(c_int), value :: n, maxrounds, round
    real(c_double) :: snew(n), s(n, maxrounds), v(n, maxrounds), r(n)
    ! Serial path of the source (Orthogonalize): x/b/r through the
    ! captured handles, S/V/R/Snew through the Rust-owned arrays.
    real(c_double) :: y(n), f(n), t1(n), t2(n)
    real(c_double) :: res, alpha, beta
    integer :: i
    type(matrix_t), pointer :: a
    real(c_double), pointer :: x(:), b(:)
    a => acc_solver % matrix
    x => acc_solver % variable % values
    b => a % rhs
    if (round == 1) then
      call cso_mymv(a, x, r)
      r(1:n) = b(1:n) - r(1:n)
      res = cso_mynorm(n, r) / cso_mynorm(n, b)
      print *, 'Residual norm for initial guess = ', res
    end if
    t1(1:n) = snew(1:n) - x(1:n)
    call cso_mymv(a, t1, t2)
    do i = 1, round - 1
      beta = cso_mydot(n, v(1:n, i), t2(1:n))
      t1(1:n) = t1(1:n) - beta * s(1:n, i)
      t2(1:n) = t2(1:n) - beta * v(1:n, i)
    end do
    alpha = cso_mynorm(n, t2)
    t1(1:n) = 1.0d0 / alpha * t1(1:n)
    t2(1:n) = 1.0d0 / alpha * t2(1:n)
    beta = cso_mydot(n, t2, r)
    x(1:n) = x(1:n) + beta * t1(1:n)
    r(1:n) = r(1:n) - beta * t2(1:n)
    s(1:n, round) = t1(1:n)
    v(1:n, round) = t2(1:n)
    res = cso_mynorm(n, r) / cso_mynorm(n, b)
    write (acc_message, '(a,I4,ES12.3)') 'OuterIteration residual for iterate', &
        round, res
    call Info('NavierStokesSolver', acc_message, level=4)
  end subroutine cso_orthogonalize

  subroutine cso_gather_vel(nd, vel) bind(c, name="cso_gather_vel")
    integer(c_int), value :: nd
    real(c_double) :: vel(*)
    ! Velocity array is (dim+1)-major like the source's Velocity(dim+1,p).
    integer :: i, k, d
    d = CoordinateSystemDimension()
    do i = 1, d
      do k = 1, nd
        vel((k - 1) * (d + 1) + i) = acc_solver % variable % prevvalues( &
            acc_solver % variable % dofs * &
            (acc_solver % variable % perm(acc_indexes(k)) - 1) + i, 1)
      end do
    end do
  end subroutine cso_gather_vel

  function cso_mydot(n, x, y) result(s)
    integer(c_int), value :: n
    real(c_double) :: x(n), y(n), s
    if (parenv % pes > 1) then
      s = ParallelDot(n, x, y)
    else
      s = dot_product(x(1:n), y(1:n))
    end if
  end function cso_mydot

  subroutine cso_mymv(a, x, b)
    type(matrix_t), pointer :: a
    real(c_double) :: x(:), b(:)
    if (parenv % pes > 1) then
      call ParallelMatrixVector(a, x, b)
    else
      call CRS_MatrixVectorMultiply(a, x, b)
    end if
  end subroutine cso_mymv

  function cso_mynorm(n, x) result(s)
    integer(c_int), value :: n
    real(c_double) :: x(n), s
    if (parenv % pes > 1) then
      s = ParallelNorm(n, x)
    else
      s = sqrt(dot_product(x(1:n), x(1:n)))
    end if
  end function cso_mynorm

  subroutine cso_parallel_update(n, maxrounds, round, snew, s, v, r) &
      bind(c, name="cso_parallel_update")
    integer(c_int), value :: n, maxrounds, round
    real(c_double) :: snew(n), s(n, maxrounds), v(n, maxrounds), r(n)
    ! Parallel-only path (transcribed, compile-checked): mirrors the
    ! source's Parallel arm with GCRUpdate.
    type(matrix_t), pointer :: a, mmat
    real(c_double), pointer :: x(:), b(:), mx(:), mb(:), mr(:)
    real(c_double), allocatable :: residual(:)
    a => acc_solver % matrix
    if (.not. associated(a % parmatrix)) &
        call ParallelInitMatrix(acc_solver, a)
    allocate(residual(n))
    residual = 0.0d0
    x => acc_solver % variable % values
    b => a % rhs
    call ParallelInitSolve(a, x, b, residual)
    call ParallelUpdateSolve(a, x, residual)
    mmat => ParallelMatrix(a, mx, mb, mr)
    call cso_gcr_update(n, a, mmat, mx, mb, mr, snew, s, v, r, round)
    call ParallelUpdateResult(a, x, residual)
    deallocate(residual)
  end subroutine cso_parallel_update

  subroutine cso_gcr_update(n, a, m, x, b, r, snew, s, v, rr, round)
    integer :: n, round
    type(matrix_t), pointer :: a, m
    real(c_double) :: x(:), b(:), r(:), snew(:)
    real(c_double) :: s(:, :), v(:, :), rr(:)
    real(c_double) :: t1(n), t2(n), beta, alpha, res
    integer :: i
    if (parenv % pes > 1) call ParallelVector(a, snew)
    if (round == 1) then
      call cso_mymv(a, x, r)
      r(1:n) = b(1:n) - r(1:n)
    else
      r(1:n) = rr(1:n)
    end if
    t1(1:n) = snew(1:n) - x(1:n)
    call cso_mymv(a, t1, t2)
    do i = 1, round - 1
      beta = cso_mydot(n, v(1:n, i), t2(1:n))
      t1(1:n) = t1(1:n) - beta * s(1:n, i)
      t2(1:n) = t2(1:n) - beta * v(1:n, i)
    end do
    alpha = cso_mynorm(n, t2)
    t1(1:n) = 1.0d0 / alpha * t1(1:n)
    t2(1:n) = 1.0d0 / alpha * t2(1:n)
    beta = cso_mydot(n, t2, r)
    x(1:n) = x(1:n) + beta * t1(1:n)
    r(1:n) = r(1:n) - beta * t2(1:n)
    s(1:n, round) = t1(1:n)
    v(1:n, round) = t2(1:n)
    rr(1:n) = r(1:n)
    res = cso_mynorm(n, r) / cso_mynorm(n, b)
    write (acc_message, '(a,I4,ES12.3)') 'OuterIteration residual for iterate', &
        round, res
    call Info('NavierStokesSolver', acc_message, level=4)
  end subroutine cso_gcr_update

  ! -- LocalMatrix context (IP loop owned by Rust) -------------------------
  subroutine cso_lm_nodes() bind(c, name="cso_lm_nodes")
    call GetElementNodes(acc_nodes_lm)
  end subroutine cso_lm_nodes

  function cso_lm_ip_begin() result(n) bind(c, name="cso_lm_ip_begin")
    integer(c_int) :: n
    acc_ip_lm = GaussPoints(acc_elem)
    n = acc_ip_lm % n
  end function cso_lm_ip_begin

  subroutine cso_lm_ip_point(t, u, v, w, s) bind(c, name="cso_lm_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip_lm % u(t)
    v = acc_ip_lm % v(t)
    w = acc_ip_lm % w(t)
    s = acc_ip_lm % s(t)
  end subroutine cso_lm_ip_point

  subroutine cso_lm_eleminfo(u, v, w, nd, basis, dbasis, detj, stat) &
      bind(c, name="cso_lm_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: nd
    real(c_double) :: basis(nd), dbasis(nd, 3), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_elem, acc_nodes_lm, u, v, w, detj, &
        basis(1:nd), dbasis(1:nd, 1:3))
    stat = merge(1, 0, ok)
  end subroutine cso_lm_eleminfo

end module cs_opt_acc
