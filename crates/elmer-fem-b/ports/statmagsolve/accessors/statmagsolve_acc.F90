! T021-B StatMagSolve port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/StatMagSolve.F90, compiled against the
! built .mod files of the target configure (A2 shipped mechanism, ADR-057
! §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, function-result arrays,
! strided POINTER views (MVP(2::2), MFD(1::3)), whole-array assignments —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module mirrors StatMagSolve's own SAVE locals (ElementNodes xyz, the
! MFD/Br/.../Joule POINTER views, Brim/Bzim, joule scratch, Wang, progress
! clock) with identical lifetimes.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module statmagsolve_acc
  use iso_c_binding
  use types
  use defutils
  use differentials
  use integration
  use matrixassembly
  use eleminfo
  use messages
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(valuelist_t), pointer, save :: acc_mat => null()
  type(nodes_t), target, save :: acc_enodes
  real(dp), pointer, save :: acc_mvp(:) => null()
  integer, pointer, save :: acc_perm(:) => null()
  real(dp), pointer, save :: acc_mfd(:) => null()
  real(dp), pointer, save :: acc_br(:) => null()
  real(dp), pointer, save :: acc_bz(:) => null()
  real(dp), pointer, save :: acc_bp(:) => null()
  real(dp), pointer, save :: acc_babs(:) => null()
  real(dp), pointer, save :: acc_joule(:) => null()
  real(dp), pointer, save :: acc_absj(:) => null()
  real(dp), allocatable, save :: acc_brim(:), acc_bzim(:)
  real(dp), allocatable, save :: acc_ae(:), acc_condj(:)
  real(dp), save :: acc_wang = 0.0_dp
  real(dp), save :: acc_at = 0.0_dp
  real(dp), save :: acc_at0 = 0.0_dp
  integer, save :: acc_n = 0

contains

  ! -- entry state ----------------------------------------------------
  subroutine sma_begin(h_model, h_solver) bind(c, name="sma_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine sma_begin

  subroutine sma_assoc_var() bind(c, name="sma_assoc_var")
    ! Solver % Variable only exists at solve time (never in _init).
    acc_mvp => acc_solver % variable % values
    acc_perm => acc_solver % variable % perm
  end subroutine sma_assoc_var


  function sma_assoc_matrix() result(r) bind(c, name="sma_assoc_matrix")
    integer(c_int) :: r
    r = merge(1, 0, associated(acc_solver % matrix))
  end function sma_assoc_matrix

  subroutine sma_info_start() bind(c, name="sma_info_start")
    call Info('StatMagSolve', ' ', Level=4)
    call Info('StatMagSolve', '------------------------------------------------', Level=4)
    call Info('StatMagSolve', 'Starting simulation', Level=4)
  end subroutine sma_info_start

  function sma_mesh_dim() result(d) bind(c, name="sma_mesh_dim")
    integer(c_int) :: d
    d = acc_solver % mesh % meshdim
  end function sma_mesh_dim

  function sma_perm_vacuum() result(mu0) bind(c, name="sma_perm_vacuum")
    real(c_double) :: mu0
    logical :: got
    mu0 = GetConstReal(acc_model % constants, 'Permeability Of Vacuum', got)
    if (.not. got) mu0 = PI * 4.0d-7
  end function sma_perm_vacuum

  function sma_count_local() result(n) bind(c, name="sma_count_local")
    integer(c_int) :: n
    n = count(acc_perm > 0)
  end function sma_count_local

  function sma_var_dofs() result(d) bind(c, name="sma_var_dofs")
    integer(c_int) :: d
    d = acc_solver % variable % dofs
  end function sma_var_dofs

  function sma_harmonic(found) result(h) bind(c, name="sma_harmonic")
    integer(c_int) :: h, found
    logical :: f
    h = merge(1, 0, ListGetLogical(acc_solver % values, 'Harmonic Simulation', f))
    found = merge(1, 0, f)
  end function sma_harmonic

  subroutine sma_info_mode(h) bind(c, name="sma_info_mode")
    integer(c_int), value :: h
    if (h /= 0) then
      call Info('StatMagSolver', 'Assuming harmonic simulation')
    else
      call Info('StatMagSolver', 'Assuming steady-state simulation')
    end if
  end subroutine sma_info_mode

  subroutine sma_set_complex() bind(c, name="sma_set_complex")
    acc_solver % matrix % complex = .true.
  end subroutine sma_set_complex

  subroutine sma_warn_whitney() bind(c, name="sma_warn_whitney")
    call Warn('StatMagSolver', &
        'This solver does not really fulfill the Coulomb gauge, use Whitney solver in 3D!')
  end subroutine sma_warn_whitney

  function sma_angfreq() result(w) bind(c, name="sma_angfreq")
    real(c_double) :: w
    w = GetAngularFrequency()
    acc_wang = w
  end function sma_angfreq

  subroutine sma_fatal_3d() bind(c, name="sma_fatal_3d")
    call Fatal('StatMagSolver', &
        'In 3D there must be three components for the vector potential!')
  end subroutine sma_fatal_3d

  subroutine sma_fatal_harm() bind(c, name="sma_fatal_harm")
    call Fatal('StatMagSolver', &
        'For harmonic simulation there must be two components for the vector potential!')
  end subroutine sma_fatal_harm

  subroutine sma_fatal_steady() bind(c, name="sma_fatal_steady")
    call Fatal('StatMagSolver', &
        'For 2d steady-state cases there should be just one component for vector potential!')
  end subroutine sma_fatal_steady

  subroutine sma_fatal_alloc() bind(c, name="sma_fatal_alloc")
    call Fatal('StatMagSolve', 'Memory allocation error.')
  end subroutine sma_fatal_alloc

  subroutine sma_warn_joule_only() bind(c, name="sma_warn_joule_only")
    call Warn('StatMagSolver', &
        'Computation of Joule heating is only relevan for harmonic systems!')
  end subroutine sma_warn_joule_only

  subroutine sma_calc_flags(mflux, mabs, mjoule) bind(c, name="sma_calc_flags")
    integer(c_int) :: mflux, mabs, mjoule
    logical :: got
    mflux = merge(1, 0, ListGetLogical(acc_solver % values, 'Calculate Magnetic Flux', got))
    mabs = merge(1, 0, ListGetLogical(acc_solver % values, 'Calculate Magnetic Flux Abs', got))
    mjoule = merge(1, 0, ListGetLogical(acc_solver % values, 'Calculate Joule Heating', got))
  end subroutine sma_calc_flags

  function sma_max_elem_nodes() result(n) bind(c, name="sma_max_elem_nodes")
    integer(c_int) :: n
    n = acc_model % maxelementnodes
  end function sma_max_elem_nodes

  function sma_num_nodes() result(n) bind(c, name="sma_num_nodes")
    integer(c_int) :: n
    n = acc_model % numberofnodes
  end function sma_num_nodes

  ! -- first-call SAVE allocation (mirrors the AllocationsDone block) --
  subroutine sma_allocate(n) bind(c, name="sma_allocate")
    integer(c_int), value :: n
    acc_n = n
    allocate(acc_enodes % x(n), acc_enodes % y(n), acc_enodes % z(n))
  end subroutine sma_allocate

  subroutine sma_assoc_flux(mabs, harm) bind(c, name="sma_assoc_flux")
    integer(c_int), value :: mabs, harm
    type(variable_t), pointer :: tv
    integer :: istat
    tv => VariableGet(acc_solver % mesh % variables, 'Magnetic Flux Density')
    acc_mfd => tv % values
    acc_br => acc_mfd(1::3)
    acc_bz => acc_mfd(2::3)
    acc_bp => acc_mfd(3::3)
    if (mabs /= 0) then
      tv => VariableGet(acc_solver % mesh % variables, 'Magnetic Flux Density_abs')
      acc_babs => tv % values
    end if
    if (harm /= 0) then
      allocate(acc_brim(acc_model % numberofnodes), &
          acc_bzim(acc_model % numberofnodes), stat=istat)
      if (istat /= 0) call Fatal('StatMagSolve', 'Memory allocation error.')
    end if
  end subroutine sma_assoc_flux

  subroutine sma_assoc_joule() bind(c, name="sma_assoc_joule")
    type(variable_t), pointer :: tv
    tv => VariableGet(acc_solver % mesh % variables, 'Joule Heating')
    acc_joule => tv % values
    tv => VariableGet(acc_solver % mesh % variables, 'Joule Field')
    acc_absj => tv % values
    allocate(acc_ae(acc_n), acc_condj(acc_n))
  end subroutine sma_assoc_joule

  function sma_maxiter(found) result(ni) bind(c, name="sma_maxiter")
    integer(c_int) :: ni, found
    logical :: f
    ni = ListGetInteger(acc_solver % values, 'Nonlinear System Max Iterations', f)
    found = merge(1, 0, f)
    if (found == 0) ni = 1
  end function sma_maxiter

  subroutine sma_info_iter(iter) bind(c, name="sma_info_iter")
    integer(c_int), value :: iter
    write(message, '(A,I0)') 'Magnetic Field Iteration: ', iter
    call Info('StatMagSolve', message, Level=4)
  end subroutine sma_info_iter

  subroutine sma_prog_begin() bind(c, name="sma_prog_begin")
    acc_at = CPUTime()
    acc_at0 = RealTime()
  end subroutine sma_prog_begin

  subroutine sma_prog_check(t) bind(c, name="sma_prog_check")
    integer(c_int), value :: t
    if (RealTime() - acc_at0 > 1.0) then
      write(message, '(a,i3,a)') '   Assembly: ', &
          int(100.0 - 100.0 * (acc_model % numberofbulkelements - t) / &
          (1.0 * acc_model % numberofbulkelements)), ' % done'
      call Info('StatMagSolve', message, Level=5)
      acc_at0 = RealTime()
    end if
  end subroutine sma_prog_check

  subroutine sma_default_initialize() bind(c, name="sma_default_initialize")
    call DefaultInitialize()
  end subroutine sma_default_initialize

  function sma_nof_active() result(n) bind(c, name="sma_nof_active")
    integer(c_int) :: n
    n = acc_solver % numberofactiveelements
  end function sma_nof_active

  ! -- per-element context (material lookup + node fill) ---------------
  function sma_elem_begin(t) result(n) bind(c, name="sma_elem_begin")
    integer(c_int), value :: t
    integer(c_int) :: n
    integer :: k, body_id
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes()
    acc_enodes % x(1:n) = acc_model % nodes % x(acc_elem % nodeindexes(1:n))
    acc_enodes % y(1:n) = acc_model % nodes % y(acc_elem % nodeindexes(1:n))
    acc_enodes % z(1:n) = acc_model % nodes % z(acc_elem % nodeindexes(1:n))
    body_id = acc_elem % bodyid
    k = ListGetInteger(acc_model % bodies(body_id) % values, 'Material', &
        minv=1, maxv=acc_model % numberofmaterials)
    acc_mat => acc_model % materials(k) % values
  end function sma_elem_begin

  function sma_bf_id() result(bf) bind(c, name="sma_bf_id")
    integer(c_int) :: bf
    logical :: got
    integer :: body_id
    body_id = acc_elem % bodyid
    bf = ListGetInteger(acc_model % bodies(body_id) % values, &
        'Body Force', got, minv=1, maxv=acc_model % numberofbodyforces)
  end function sma_bf_id

  subroutine sma_perm_fill(n, ld, perm, mu0) bind(c, name="sma_perm_fill")
    integer(c_int), value :: n, ld
    real(c_double) :: perm(ld)
    real(c_double), value :: mu0
    logical :: got
    perm(1:n) = ListGetReal(acc_mat, 'Relative Permeability', n, &
        acc_elem % nodeindexes, got)
    if (got) then
      perm(1:n) = mu0 * perm(1:n)
    else
      perm(1:n) = ListGetReal(acc_mat, 'Magnetic Permeability', n, &
          acc_elem % nodeindexes, got)
      if (.not. got) perm = mu0
    end if
  end subroutine sma_perm_fill

  subroutine sma_reluct(n, reluc, perm) bind(c, name="sma_reluct")
    integer(c_int), value :: n
    real(c_double) :: reluc(*), perm(*)
    reluc(1:n) = 1.0_dp / perm(1:n)
  end subroutine sma_reluct

  subroutine sma_cond_harm(n, cond) bind(c, name="sma_cond_harm")
    integer(c_int), value :: n
    real(c_double) :: cond(*)
    logical :: got
    cond(1:n) = ListGetReal(acc_mat, 'Electrical Conductivity', n, &
        acc_elem % nodeindexes, got)
    if (got) then
      call Warn('StatMagSolve', 'Use > Electric Conductivity < instead of electrial')
    else
      cond(1:n) = ListGetReal(acc_mat, 'Electric Conductivity', n, &
          acc_elem % nodeindexes)
    end if
  end subroutine sma_cond_harm

  subroutine sma_current(n, bf, cur) bind(c, name="sma_current")
    integer(c_int), value :: n, bf
    real(c_double) :: cur(*)
    logical :: got
    if (bf > 0) then
      cur(1:n) = ListGetReal(acc_model % bodyforces(bf) % values, &
          'Current Density', n, acc_elem % nodeindexes, got)
    else
      cur(1:n) = 0.0d0
    end if
  end subroutine sma_current

  subroutine sma_phase(n, bf, phase) bind(c, name="sma_phase")
    integer(c_int), value :: n, bf
    real(c_double) :: phase(*)
    logical :: got
    if (bf > 0) then
      phase(1:n) = ListGetReal(acc_model % bodyforces(bf) % values, &
          'Current Phase Angle', n, acc_elem % nodeindexes, got)
    else
      phase(1:n) = 0.0d0
    end if
  end subroutine sma_phase

  subroutine sma_vecload(n, ldn, bf, vlv) bind(c, name="sma_vecload")
    integer(c_int), value :: n, ldn, bf
    real(c_double) :: vlv(3, ldn)
    logical :: got
    vlv = 0.0_dp
    if (bf > 0) then
      vlv(1, 1:n) = vlv(1, 1:n) + ListGetReal(acc_model % bodyforces(bf) % values, &
          'Current Density 1', n, acc_elem % nodeindexes, got)
      vlv(2, 1:n) = vlv(2, 1:n) + ListGetReal(acc_model % bodyforces(bf) % values, &
          'Current Density 2', n, acc_elem % nodeindexes, got)
      vlv(3, 1:n) = vlv(3, 1:n) + ListGetReal(acc_model % bodyforces(bf) % values, &
          'Current Density 3', n, acc_elem % nodeindexes, got)
    end if
  end subroutine sma_vecload

  ! -- local assemblers (verbatim CONTAINS transcriptions) -------------
  subroutine sma_stat_axis(ld, n, stiff, force, mass, load, reluc) &
      bind(c, name="sma_stat_axis")
    integer(c_int), value :: ld, n
    real(c_double) :: stiff(ld, ld), force(ld), mass(ld, ld), load(ld), reluc(ld)
    real(dp) :: basis(n), dbasisdx(n, 3), ddbasisddx(n, 3, 3)
    real(dp) :: sqrtmetric, fint, r, reluct, mat
    integer :: ddim, t, p, q
    real(dp) :: s, u, v, w
    type(GaussIntegrationPoints_t), target :: integ
    integer :: nip
    real(dp), dimension(:), pointer :: uint, vint, wint, sint
    logical :: stat
    ddim = 2
    force = 0.0d0
    mass = 0.0d0
    stiff = 0.0d0
    integ = GaussPoints(acc_elem)
    uint => integ % u
    vint => integ % v
    wint => integ % w
    sint => integ % s
    nip = integ % n
    do t = 1, nip
      u = uint(t)
      v = vint(t)
      w = wint(t)
      stat = ElementInfo(acc_elem, acc_enodes, u, v, w, sqrtmetric, &
          basis, dbasisdx, ddbasisddx, .false.)
      r = sum(basis(1:n) * acc_enodes % x(1:n))
      s = sqrtmetric * sint(t)
      fint = sum(load(1:n) * basis(1:n))
      reluct = sum(reluc(1:n) * basis(1:n))
      do p = 1, n
        do q = 1, n
          mat = r * dbasisdx(p, 1) * dbasisdx(q, 1) + &
              r * dbasisdx(p, 2) * dbasisdx(q, 2) + &
              basis(p) * dbasisdx(q, 1) + &
              basis(q) * dbasisdx(p, 1) + &
              basis(p) * basis(q) / r
          mat = mat * reluct * s
          stiff(p, q) = stiff(p, q) + mat
        end do
      end do
      do p = 1, n
        force(p) = force(p) + fint * basis(p) * r * s
      end do
    end do
  end subroutine sma_stat_axis

  subroutine sma_harm_axis(ld, n, stiff, force, curdens, phase, reluc, cond, wang) &
      bind(c, name="sma_harm_axis")
    integer(c_int), value :: ld, n
    real(c_double) :: stiff(ld, ld), force(ld)
    real(c_double) :: curdens(ld), phase(ld), reluc(ld), cond(ld)
    real(c_double), value :: wang
    real(dp) :: basis(n), dbasisdx(n, 3), ddbasisddx(n, 3, 3)
    real(dp) :: sqrtmetric, fint, r, a11, a21, a12, a22, reluct, condv, angle
    integer :: ddim, t, p, q
    real(dp) :: s, u, v, w
    type(GaussIntegrationPoints_t), target :: integ
    integer :: nip
    real(dp), dimension(:), pointer :: uint, vint, wint, sint
    logical :: stat
    ddim = 2
    force = 0.0d0
    stiff = 0.0d0
    integ = GaussPoints(acc_elem)
    uint => integ % u
    vint => integ % v
    wint => integ % w
    sint => integ % s
    nip = integ % n
    do t = 1, nip
      u = uint(t)
      v = vint(t)
      w = wint(t)
      stat = ElementInfo(acc_elem, acc_enodes, u, v, w, sqrtmetric, &
          basis, dbasisdx, ddbasisddx, .false.)
      r = sum(basis * acc_enodes % x(1:n))
      s = sqrtmetric * sint(t)
      fint = 0.0d0
      fint = sum(curdens(1:n) * basis)
      reluct = sum(reluc(1:n) * basis(1:n))
      condv = sum(cond(1:n) * basis(1:n))
      angle = (PI / 180.0d0) * sum(phase(1:n) * basis(1:n))
      do p = 1, n
        do q = 1, n
          a11 = dbasisdx(p, 1) * dbasisdx(q, 1) * r + &
              dbasisdx(p, 2) * dbasisdx(q, 2) * r + &
              basis(p) * dbasisdx(q, 1) + &
              basis(q) * dbasisdx(p, 1) + &
              basis(p) * basis(q) / r
          a11 = reluct * s * a11
          a22 = a11
          a21 = -condv * wang * s * r * basis(q) * basis(p)
          a12 = -a21
          stiff(2 * p - 1, 2 * q - 1) = stiff(2 * p - 1, 2 * q - 1) + a11
          stiff(2 * p, 2 * q) = stiff(2 * p, 2 * q) + a22
          stiff(2 * p, 2 * q - 1) = stiff(2 * p, 2 * q - 1) + a21
          stiff(2 * p - 1, 2 * q) = stiff(2 * p - 1, 2 * q) + a12
        end do
      end do
      do p = 1, n
        force(2 * p - 1) = force(2 * p - 1) + fint * cos(angle) * basis(p) * r * s
        force(2 * p) = force(2 * p) + fint * sin(angle) * basis(p) * r * s
      end do
    end do
  end subroutine sma_harm_axis

  subroutine sma_cart_3d(ld, n, ldn, stiff, force, vlv, reluc) &
      bind(c, name="sma_cart_3d")
    integer(c_int), value :: ld, n, ldn
    real(c_double) :: stiff(ld, ld), force(ld), vlv(3, ldn), reluc(ldn)
    real(dp) :: basis(n), dbasisdx(n, 3), ddbasisddx(n, 3, 3)
    real(dp) :: sqrtmetric, fint(3), mat, reluct
    real(dp) :: s, u, v, w
    type(GaussIntegrationPoints_t), target :: integ
    integer :: nip
    real(dp), dimension(:), pointer :: uint, vint, wint, sint
    integer :: t, p, q, i
    logical :: stat
    force = 0.0d0
    stiff = 0.0d0
    integ = GaussPoints(acc_elem)
    uint => integ % u
    vint => integ % v
    wint => integ % w
    sint => integ % s
    nip = integ % n
    do t = 1, nip
      u = uint(t)
      v = vint(t)
      w = wint(t)
      stat = ElementInfo(acc_elem, acc_enodes, u, v, w, sqrtmetric, &
          basis, dbasisdx, ddbasisddx, .false.)
      s = sqrtmetric * sint(t)
      do i = 1, 3
        fint(i) = sum(vlv(i, 1:n) * basis(1:n))
      end do
      reluct = sum(reluc(1:n) * basis(1:n))
      do p = 1, n
        do q = 1, n
          mat = sum(dbasisdx(p, 1:3) * dbasisdx(q, 1:3))
          mat = mat * reluct * s
          stiff(3 * p - 2, 3 * q - 2) = stiff(3 * p - 2, 3 * q - 2) + mat
          stiff(3 * p - 1, 3 * q - 1) = stiff(3 * p - 1, 3 * q - 1) + mat
          stiff(3 * p, 3 * q) = stiff(3 * p, 3 * q) + mat
        end do
      end do
      do p = 1, n
        force(3 * p - 2) = force(3 * p - 2) + fint(1) * basis(p) * s
        force(3 * p - 1) = force(3 * p - 1) + fint(2) * basis(p) * s
        force(3 * p) = force(3 * p) + fint(3) * basis(p) * s
      end do
    end do
  end subroutine sma_cart_3d

  ! -- shared assembly tail ---------------------------------------------
  subroutine sma_default_update(ld, g, f) bind(c, name="sma_default_update")
    integer(c_int), value :: ld
    real(c_double) :: g(ld, ld), f(ld)
    call DefaultUpdateEquations(g(1:ld, 1:ld), f(1:ld))
  end subroutine sma_default_update

  subroutine sma_info_asm_done() bind(c, name="sma_info_asm_done")
    call Info('StatMagSolve', 'Assembly done', Level=4)
  end subroutine sma_info_asm_done

  subroutine sma_default_finish() bind(c, name="sma_default_finish")
    call DefaultFinishAssembly()
  end subroutine sma_default_finish

  subroutine sma_default_dirichlet() bind(c, name="sma_default_dirichlet")
    call DefaultDirichletBCs()
  end subroutine sma_default_dirichlet

  subroutine sma_info_bnd_done() bind(c, name="sma_info_bnd_done")
    call Info('StatMagSolve', 'Set boundaries done', Level=4)
  end subroutine sma_info_bnd_done

  function sma_default_solve() result(norm) bind(c, name="sma_default_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function sma_default_solve

  ! -- magnetic flux density post-processing ------------------------------
  subroutine sma_flux_zero() bind(c, name="sma_flux_zero")
    acc_mfd = 0.0d0
  end subroutine sma_flux_zero

  subroutine sma_flux_harm() bind(c, name="sma_flux_harm")
    real(dp), pointer :: a(:)
    integer :: i, j
    a => acc_mvp(2::2)
    call AxiSCurl(acc_bp, acc_bp, a, acc_brim, acc_bzim, acc_bp, acc_perm)
    a => acc_mvp(1::2)
    call AxiSCurl(acc_bp, acc_bp, a, acc_br, acc_bz, acc_bp, acc_perm)
    do i = 1, acc_model % numberofnodes
      j = acc_perm(i)
      if (j > 0) then
        acc_br(j) = sqrt(acc_br(j)**2 + acc_brim(j)**2)
        acc_bz(j) = sqrt(acc_bz(j)**2 + acc_bzim(j)**2)
        acc_brim(j) = 0.0d0
        acc_bzim(j) = 0.0d0
      end if
    end do
  end subroutine sma_flux_harm

  subroutine sma_flux_steady(dim) bind(c, name="sma_flux_steady")
    integer(c_int), value :: dim
    real(dp), pointer :: a(:), ax(:), ay(:), az(:)
    if (dim < 3) then
      a => acc_mvp
      call AxiSCurl(acc_bp, acc_bp, a, acc_br, acc_bz, acc_bp, acc_perm)
    else
      ax => acc_mvp(1::3)
      ay => acc_mvp(2::3)
      az => acc_mvp(3::3)
      acc_br = 0.0d0
      acc_bp = 0.0d0
      acc_bz = 0.0d0
      call Curl(ax, ay, az, acc_br, acc_bz, acc_bp, acc_perm)
    end if
  end subroutine sma_flux_steady

  subroutine sma_invalidate(which) bind(c, name="sma_invalidate")
    integer(c_int), value :: which
    if (which == 0) then
      call InvalidateVariable(acc_model % meshes, acc_solver % mesh, &
          'Magnetic Flux Density')
    else
      call InvalidateVariable(acc_model % meshes, acc_solver % mesh, &
          'Magnetic Flux Density_abs')
    end if
  end subroutine sma_invalidate

  subroutine sma_babs() bind(c, name="sma_babs")
    integer :: i
    do i = 1, size(acc_babs)
      acc_babs(i) = sqrt(acc_br(i)**2 + acc_bz(i)**2 + acc_bp(i)**2)
    end do
  end subroutine sma_babs

  function sma_nonlin_converged() result(r) bind(c, name="sma_nonlin_converged")
    integer(c_int) :: r
    r = merge(1, 0, acc_solver % variable % nonlinconverged == 1)
  end function sma_nonlin_converged

  subroutine sma_info_converged(iter) bind(c, name="sma_info_converged")
    integer(c_int), value :: iter
    write(message, '(A,I0,A)') 'Convergence after ', iter, ' iterations'
    call Info('StatMagSolve', message, Level=4)
  end subroutine sma_info_converged

  ! -- Joule heating section (whole IF block, verbatim) --------------------
  subroutine sma_joule_section() bind(c, name="sma_joule_section")
    integer :: i, j, n, t, k, body_id
    integer, pointer :: ni(:)
    real(dp) :: jc, jre, jim, totalheating, totalvolume, desired
    type(valuelist_t), pointer :: mat
    logical :: got
    jc = 0.5_dp * acc_wang**2
    do i = 1, acc_model % numberofnodes
      j = acc_perm(i)
      if (j > 0) then
        jre = acc_mvp(2 * j - 1)
        jim = acc_mvp(2 * j)
        acc_absj(j) = jc * (jre * jre + jim * jim)
      end if
    end do
    totalheating = 0.0d0
    totalvolume = 0.0d0
    do t = 1, acc_solver % numberofactiveelements
      acc_elem => GetActiveElement(t)
      ni => acc_elem % nodeindexes
      n = GetElementNOFNodes()
      acc_enodes % x(1:n) = acc_model % nodes % x(ni)
      acc_enodes % y(1:n) = acc_model % nodes % y(ni)
      acc_enodes % z(1:n) = acc_model % nodes % z(ni)
      body_id = acc_elem % bodyid
      k = ListGetInteger(acc_model % bodies(body_id) % values, 'Material', &
          minv=1, maxv=acc_model % numberofmaterials)
      mat => acc_model % materials(k) % values
      acc_condj(1:n) = ListGetReal(mat, 'Electrical Conductivity', n, ni, got)
      if (got) then
        call Warn('StatMagSolve', 'Use electric conductivity instead of electrical')
      else
        acc_condj(1:n) = ListGetReal(mat, 'Electric Conductivity', n, ni)
      end if
      acc_ae(1:n) = acc_absj(acc_perm(ni(1:n)))
      call sma_joule_integrate(n, totalheating, totalvolume)
      do i = 1, n
        j = acc_perm(ni(i))
        if (j > 0) then
          jc = acc_condj(i) * acc_absj(j)
          if (jc > acc_joule(j)) acc_joule(j) = jc
        end if
      end do
    end do
    desired = ListGetConstReal(acc_solver % values, 'Desired Heating Power', got)
    if (.not. got) desired = ListGetConstReal(acc_solver % values, 'Power Control', got)
    if (got .and. totalheating > 0.0d0) then
      acc_absj = (desired / totalheating) * acc_absj
      acc_joule = (desired / totalheating) * acc_joule
    end if
    write(message, '(A,ES15.4)') 'Joule Heating (W): ', totalheating
    call Info('StatMagSolve', message, Level=4)
    call ListAddConstReal(acc_model % simulation, 'res: Joule heating', totalheating)
  end subroutine sma_joule_section

  subroutine sma_joule_integrate(n, totalheating, totalvolume) 
    integer, intent(in) :: n
    real(dp), intent(inout) :: totalheating, totalvolume
    real(dp) :: basis(n), dbasisdx(n, 3), ddbasisddx(n, 3, 3)
    real(dp) :: sqrtmetric, conductivity, field
    integer :: t
    real(dp) :: r, s, u, v, w
    type(GaussIntegrationPoints_t), target :: integ
    integer :: nip, p
    real(dp), dimension(:), pointer :: uint, vint, wint, sint
    logical :: stat
    integ = GaussPoints(acc_elem)
    uint => integ % u
    vint => integ % v
    wint => integ % w
    sint => integ % s
    nip = integ % n
    do t = 1, nip
      u = uint(t)
      v = vint(t)
      w = wint(t)
      stat = ElementInfo(acc_elem, acc_enodes, u, v, w, sqrtmetric, &
          basis, dbasisdx, ddbasisddx, .false.)
      r = sum(basis * acc_enodes % x(1:n))
      s = sqrtmetric * sint(t)
      field = sum(acc_ae(1:n) * basis(1:n))
      conductivity = sum(acc_condj(1:n) * basis(1:n))
      do p = 1, n
        totalvolume = totalvolume + 2.0d0 * PI * r * s * basis(p)
        totalheating = totalheating + 2.0d0 * PI * r * s * basis(p) * field * conductivity
      end do
    end do
  end subroutine sma_joule_integrate

  subroutine sma_info_done() bind(c, name="sma_info_done")
    call Info('StatMagSolve', 'All done, exiting', Level=4)
    call Info('StatMagSolve', '------------------------------------------------', Level=4)
  end subroutine sma_info_done

  ! -- Init entry body (verbatim StatMagSolver_Init) ----------------------
  subroutine sma_init_run() bind(c, name="sma_init_run")
    type(valuelist_t), pointer :: params
    logical :: found, calculate, calculateabs
    integer :: dim, givendim
    params => GetSolverParams()
    dim = CoordinateSystemDimension()
    givendim = dim
    if (GetLogical(params, 'Harmonic Simulation', found)) then
      call ListAddInteger(params, 'Variable Dofs', 2)
    end if
    calculate = ListGetLogical(params, 'Calculate Magnetic Flux', found)
    calculateabs = ListGetLogical(params, 'Calculate Magnetic Flux Abs', found)
    if (calculateabs .and. .not. calculate) then
      call Warn('StatMagSolver_init', 'Cannot compute Abs without computing field')
    end if
    if (calculate) then
      call ListAddString(params, NextFreeKeyword('Exported Variable', params), &
          '-dofs 3 Magnetic Flux Density')
      if (calculateabs) then
        call ListAddString(params, NextFreeKeyword('Exported Variable', params), &
            'Magnetic Flux Density_abs')
      end if
    end if
    calculate = ListGetLogical(params, 'Calculate Joule Heating', found)
    if (calculate) then
      call ListAddString(params, NextFreeKeyword('Exported Variable', params), &
          'Joule Heating')
      call ListAddString(params, NextFreeKeyword('Exported Variable', params), &
          'Joule Field')
    end if
  end subroutine sma_init_run

end module statmagsolve_acc
