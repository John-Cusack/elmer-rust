! T021-B ModelPDE port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/ModelPDE.F90, compiled against the built .mod
! files of the target configure (T020 Poisson pattern, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module mirrors ModelPDE's own locals (current element, Nodes SAVE state
! per internal procedure, material/BC/body-force pointers, IP state) with
! identical lifetimes.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module modelpde_acc
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
  type(valuelist_t), pointer, save :: acc_mat => null()
  type(valuelist_t), pointer, save :: acc_bc => null()
  ! LocalMatrix context (mirrors its SAVE Nodes and IP local).
  type(nodes_t), target, save :: acc_nodes_lm
  type(gaussintegrationpoints_t), save :: acc_ip_lm
  ! LocalMatrixBC context (mirrors its SAVE Nodes and IP local).
  type(nodes_t), target, save :: acc_nodes_bc
  type(gaussintegrationpoints_t), save :: acc_ip_bc

contains

  ! -- entry state ----------------------------------------------------
  subroutine moa_begin(h_solver) bind(c, name="moa_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine moa_begin

  subroutine moa_default_start() bind(c, name="moa_default_start")
    call DefaultStart()
  end subroutine moa_default_start

  function moa_max_iter() result(m) bind(c, name="moa_max_iter")
    integer(c_int) :: m
    logical :: found
    m = ListGetInteger(GetSolverParams(), &
        'Nonlinear System Max Iterations', found, minv=1)
    if (.not. found) m = 1
  end function moa_max_iter

  ! -- assembly driver -------------------------------------------------
  function moa_nof_active() result(n) bind(c, name="moa_nof_active")
    integer(c_int) :: n
    n = GetNOFActive()
  end function moa_nof_active

  subroutine moa_default_initialize() bind(c, name="moa_default_initialize")
    call DefaultInitialize()
  end subroutine moa_default_initialize

  subroutine moa_begin_active(t, n, nd, nb) bind(c, name="moa_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, nb
    acc_elem => GetActiveElement(t)
    n  = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
    nb = GetElementNOFBDOFs()
  end subroutine moa_begin_active

  subroutine moa_finish_bulk() bind(c, name="moa_finish_bulk")
    call DefaultFinishBulkAssembly()
  end subroutine moa_finish_bulk

  function moa_num_bnd() result(n) bind(c, name="moa_num_bnd")
    integer(c_int) :: n
    n = GetNOFBoundaryElements()
  end function moa_num_bnd

  subroutine moa_begin_boundary(t, n, nd, active) bind(c, name="moa_begin_boundary")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, active
    acc_elem => GetBoundaryElement(t)
    active = merge(1, 0, ActiveBoundaryElement())
    n  = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
  end subroutine moa_begin_boundary

  function moa_cut_fem() result(r) bind(c, name="moa_cut_fem")
    integer(c_int) :: r
    r = merge(1, 0, DefaultCutFEM())
  end function moa_cut_fem

  subroutine moa_finish_boundary() bind(c, name="moa_finish_boundary")
    call DefaultFinishBoundaryAssembly()
  end subroutine moa_finish_boundary

  subroutine moa_finish_assembly() bind(c, name="moa_finish_assembly")
    call DefaultFinishAssembly()
  end subroutine moa_finish_assembly

  subroutine moa_dirichlet() bind(c, name="moa_dirichlet")
    call DefaultDirichletBCs()
  end subroutine moa_dirichlet

  function moa_solve() result(norm) bind(c, name="moa_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function moa_solve

  function moa_converged() result(r) bind(c, name="moa_converged")
    integer(c_int) :: r
    r = merge(1, 0, DefaultConverged())
  end function moa_converged

  subroutine moa_default_finish() bind(c, name="moa_default_finish")
    call DefaultFinish()
  end subroutine moa_default_finish

  subroutine moa_report(totarea, totlen, totsrc) bind(c, name="moa_report")
    real(c_double), value :: totarea, totlen, totsrc
    logical :: found
    if (ListGetLogical(GetSolverParams(), 'CutFEM', found) &
        .or. ListGetLogical(GetSolverParams(), 'Integ Test', found)) then
      call ListAddConstReal(CurrentModel % Simulation, 'res: integ total area', totarea)
      call ListAddConstReal(CurrentModel % Simulation, 'res: integ total len', totlen)
      call ListAddConstReal(CurrentModel % Simulation, 'res: integ total src', totsrc)
    end if
  end subroutine moa_report

  ! -- LocalMatrix context ----------------------------------------------
  subroutine moa_lm_nodes() bind(c, name="moa_lm_nodes")
    call GetElementNodes(acc_nodes_lm)
  end subroutine moa_lm_nodes

  function moa_bodyforce() result(r) bind(c, name="moa_bodyforce")
    integer(c_int) :: r
    acc_bf => GetBodyForce()
    r = merge(1, 0, associated(acc_bf))
  end function moa_bodyforce

  subroutine moa_load_source(n, load, found) bind(c, name="moa_load_source")
    integer(c_int), value :: n
    real(c_double) :: load(n)
    integer(c_int) :: found
    logical :: f
    load(1:n) = GetReal(acc_bf, 'field source', f)
    found = merge(1, 0, f)
  end subroutine moa_load_source

  subroutine moa_mat_begin() bind(c, name="moa_mat_begin")
    acc_mat => GetMaterial()
  end subroutine moa_mat_begin

  subroutine moa_mat_coeff(which, n, out, found) bind(c, name="moa_mat_coeff")
    integer(c_int), value :: which, n
    real(c_double) :: out(n)
    integer(c_int) :: found
    logical :: f
    select case (which)
    case (1)
      out(1:n) = GetReal(acc_mat, 'diffusion coefficient', f)
    case (2)
      out(1:n) = GetReal(acc_mat, 'reaction coefficient', f)
    case (3)
      out(1:n) = GetReal(acc_mat, 'convection coefficient', f)
    case default
      out(1:n) = GetReal(acc_mat, 'time derivative coefficient', f)
    end select
    found = merge(1, 0, f)
  end subroutine moa_mat_coeff

  subroutine moa_conv_vel(i, n, velo) bind(c, name="moa_conv_vel")
    integer(c_int), value :: i, n
    real(c_double) :: velo(3, n)
    logical :: f
    velo(i, 1:n) = GetReal(acc_mat, 'convection velocity '//I2S(i), f)
  end subroutine moa_conv_vel

  function moa_coord_dim() result(d) bind(c, name="moa_coord_dim")
    integer(c_int) :: d
    d = CoordinateSystemDimension()
  end function moa_coord_dim

  function moa_lm_ip_begin() result(n) bind(c, name="moa_lm_ip_begin")
    integer(c_int) :: n
    acc_ip_lm = GaussPointsAdapt(acc_elem)
    if (acc_elem % ElementIndex == 1) then
      call Info('AdvDiffSolver', &
          'Integration points in 1st element: '//I2S(acc_ip_lm % n), Level=8)
    end if
    n = acc_ip_lm % n
  end function moa_lm_ip_begin

  subroutine moa_lm_ip_point(t, u, v, w, s) bind(c, name="moa_lm_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip_lm % u(t)
    v = acc_ip_lm % v(t)
    w = acc_ip_lm % w(t)
    s = acc_ip_lm % s(t)
  end subroutine moa_lm_ip_point

  subroutine moa_lm_eleminfo(u, v, w, n, nd, basis, dbasis, detj, stat) &
      bind(c, name="moa_lm_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: n, nd
    real(c_double) :: basis(nd), dbasis(nd, 3), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_elem, acc_nodes_lm, u, v, w, detj, &
        basis(1:nd), dbasis(1:nd, 1:3))
    stat = merge(1, 0, ok)
  end subroutine moa_lm_eleminfo

  ! Exact transcription of the LocalMatrix integration-point body: all
  ! SUM/MATMUL reductions, the diffusion/advection/reaction/time updates,
  ! and the TotArea/TotSrc accumulation keep the original evaluation order.
  subroutine moa_lm_ip(n, nd, dim, basis, dbasis, load, timec, velo, &
      diffc, convc, reactc, weight, stiff, mass, force, loadatip, &
      totarea, totsrc) bind(c, name="moa_lm_ip")
    integer(c_int), value :: n, nd, dim
    real(c_double) :: basis(nd), dbasis(nd, 3)
    real(c_double) :: load(n), timec(n), velo(3, n)
    real(c_double) :: diffc(n), convc(n), reactc(n)
    real(c_double), value :: weight
    real(c_double) :: stiff(nd, nd), mass(nd, nd), force(nd)
    real(c_double) :: loadatip, totarea, totsrc
    real(c_double) :: rho, avec(3), dcoef, ccoef, rcoef
    integer :: p, q
    loadatip = sum(basis(1:n) * load(1:n))
    rho = sum(basis(1:n) * timec(1:n))
    avec = matmul(velo(:, 1:n), basis(1:n))
    dcoef = sum(basis(1:n) * diffc(1:n))
    ccoef = sum(basis(1:n) * convc(1:n))
    rcoef = sum(basis(1:n) * reactc(1:n))
    stiff(1:nd, 1:nd) = stiff(1:nd, 1:nd) + weight * &
        dcoef * matmul(dbasis, transpose(dbasis))
    do p = 1, nd
      do q = 1, nd
        stiff(p, q) = stiff(p, q) + weight * &
            ccoef * sum(avec(1:dim) * dbasis(q, 1:dim)) * basis(p)
        stiff(p, q) = stiff(p, q) + weight * rcoef * basis(q) * basis(p)
        mass(p, q) = mass(p, q) + weight * rho * basis(q) * basis(p)
      end do
    end do
    force(1:nd) = force(1:nd) + weight * loadatip * basis(1:nd)
    totarea = totarea + weight
    totsrc = totsrc + weight * loadatip
  end subroutine moa_lm_ip

  subroutine moa_lm_transient(nd, mass, stiff, force) &
      bind(c, name="moa_lm_transient")
    integer(c_int), value :: nd
    real(c_double) :: mass(nd, nd), stiff(nd, nd), force(nd)
    call Default1stOrderTime(mass, stiff, force)
  end subroutine moa_lm_transient

  subroutine moa_condensate(nd, nb, ld, k, f) bind(c, name="moa_condensate")
    integer(c_int), value :: nd, nb, ld
    real(c_double) :: k(ld, ld), f(ld)
    call CondensateP(nd, nb, k(1:ld, 1:ld), f(1:ld))
  end subroutine moa_condensate

  subroutine moa_default_update(ld, g, f) bind(c, name="moa_default_update")
    integer(c_int), value :: ld
    real(c_double) :: g(ld, ld), f(ld)
    call DefaultUpdateEquations(g(1:ld, 1:ld), f(1:ld))
  end subroutine moa_default_update

  ! -- LocalMatrixBC context ----------------------------------------------
  function moa_bc_begin() result(r) bind(c, name="moa_bc_begin")
    integer(c_int) :: r
    acc_bc => GetBC()
    r = merge(1, 0, associated(acc_bc))
  end function moa_bc_begin

  subroutine moa_bc_nodes() bind(c, name="moa_bc_nodes")
    call GetElementNodes(acc_nodes_bc)
  end subroutine moa_bc_nodes

  subroutine moa_bc_coeffs(n, flux, coeff, ext) bind(c, name="moa_bc_coeffs")
    integer(c_int), value :: n
    real(c_double) :: flux(n), coeff(n), ext(n)
    logical :: f
    flux(1:n)  = GetReal(acc_bc, 'field flux', f)
    coeff(1:n) = GetReal(acc_bc, 'robin coefficient', f)
    ext(1:n)   = GetReal(acc_bc, 'external field', f)
  end subroutine moa_bc_coeffs

  function moa_bc_ip_begin() result(n) bind(c, name="moa_bc_ip_begin")
    integer(c_int) :: n
    acc_ip_bc = GaussPoints(acc_elem)
    n = acc_ip_bc % n
  end function moa_bc_ip_begin

  subroutine moa_bc_ip_point(t, u, v, w, s) bind(c, name="moa_bc_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip_bc % u(t)
    v = acc_ip_bc % v(t)
    w = acc_ip_bc % w(t)
    s = acc_ip_bc % s(t)
  end subroutine moa_bc_ip_point

  subroutine moa_bc_eleminfo(u, v, w, n, nd, basis, dbasis, detj, stat) &
      bind(c, name="moa_bc_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: n, nd
    real(c_double) :: basis(nd), dbasis(nd, 3), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_elem, acc_nodes_bc, u, v, w, detj, &
        basis(1:nd), dbasis(1:nd, 1:3))
    stat = merge(1, 0, ok)
  end subroutine moa_bc_eleminfo

  ! Exact transcription of the LocalMatrixBC integration-point body.
  subroutine moa_bc_ip(n, nd, basis, flux, coeff, extt, weight, stiff, &
      force, totlen) bind(c, name="moa_bc_ip")
    integer(c_int), value :: n, nd
    real(c_double) :: basis(nd)
    real(c_double) :: flux(n), coeff(n), extt(n)
    real(c_double), value :: weight
    real(c_double) :: stiff(nd, nd), force(nd)
    real(c_double) :: totlen
    real(c_double) :: f, c, ext
    integer :: p, q
    f = sum(basis(1:n) * flux(1:n))
    c = sum(basis(1:n) * coeff(1:n))
    ext = sum(basis(1:n) * extt(1:n))
    do p = 1, nd
      do q = 1, nd
        stiff(p, q) = stiff(p, q) + weight * c * basis(q) * basis(p)
      end do
    end do
    force(1:nd) = force(1:nd) + weight * (f + c * ext) * basis(1:nd)
    totlen = totlen + weight
  end subroutine moa_bc_ip

  subroutine moa_fatal_alloc() bind(c, name="moa_fatal_alloc")
    call Fatal('AdvDiffSolver', 'Memory allocation error.')
  end subroutine moa_fatal_alloc

end module modelpde_acc
