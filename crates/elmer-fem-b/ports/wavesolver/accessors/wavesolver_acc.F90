! T021-B WaveSolver port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/WaveSolver.F90, compiled against the built
! .mod files of the target configure (T020 Poisson pattern, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module mirrors WaveSolver's SAVE locals (per-internal-procedure
! ValueHandle sets and Nodes, plus the Gauss-point state) with identical
! lifetimes.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module wavesolver_acc
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
  ! LocalMatrix context (mirrors its SAVE handles and Nodes).
  type(valuehandle_t), save :: acc_lm_load, acc_lm_ss, acc_lm_damp
  type(valuehandle_t), save :: acc_lm_rho, acc_lm_react
  type(nodes_t), target, save :: acc_nodes_lm
  type(gaussintegrationpoints_t), save :: acc_ip_lm
  ! LocalMatrixBC context (mirrors its SAVE handles and Nodes).
  type(valuehandle_t), save :: acc_bc_flux, acc_bc_ss, acc_bc_rho
  type(nodes_t), target, save :: acc_nodes_lbc
  type(gaussintegrationpoints_t), save :: acc_ip_lbc

contains

  ! -- WaveSolver_Init entry ------------------------------------------
  subroutine wva_init() bind(c, name="wva_init")
    type(valuelist_t), pointer :: pars
    pars => GetSolverParams()
    call ListAddInteger(pars, 'Time derivative order', 2)
  end subroutine wva_init

  ! -- entry state ------------------------------------------------------
  subroutine wva_begin(h_solver) bind(c, name="wva_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine wva_begin

  subroutine wva_info_start() bind(c, name="wva_info_start")
    call Info('WaveSolver', 'Solving the compressional pressure wave')
  end subroutine wva_info_start

  function wva_max_iter() result(n) bind(c, name="wva_max_iter")
    integer(c_int) :: n
    type(valuelist_t), pointer :: pars
    logical :: found
    pars => GetSolverParams()
    n = ListGetInteger(pars, 'Nonlinear System Max Iterations', found, minv=1)
    if (.not. found) n = 1
  end function wva_max_iter

  function wva_coord_dim() result(n) bind(c, name="wva_coord_dim")
    integer(c_int) :: n
    n = CoordinateSystemDimension()
  end function wva_coord_dim

  function wva_need_mass() result(r) bind(c, name="wva_need_mass")
    integer(c_int) :: r
    type(valuelist_t), pointer :: pars
    logical :: found, e, h1, h2
    pars => GetSolverParams()
    e = EigenOrHarmonicAnalysis()
    h1 = GetLogical(pars, 'Harmonic Mode', found)
    h2 = GetLogical(pars, 'Harmonic Analysis', found)
    if (e .or. h1 .or. h2) then
      call Info('WaveSolver', 'We have a harmonic or eigenmode system')
      r = 1
    else
      r = 0
    end if
  end function wva_need_mass

  function wva_nonlin_converged() result(n) bind(c, name="wva_nonlin_converged")
    integer(c_int) :: n
    n = acc_solver % variable % nonlinconverged
  end function wva_nonlin_converged

  subroutine wva_avg_zero() bind(c, name="wva_avg_zero")
    logical :: found
    real(kind=dp) :: pave
    if (GetLogical(acc_solver % values, 'Set Average To Zero', found)) then
      pave = sum(acc_solver % variable % values) / &
          size(acc_solver % variable % values)
      acc_solver % variable % values = acc_solver % variable % values - pave
    end if
  end subroutine wva_avg_zero

  ! -- assembly driver --------------------------------------------------
  subroutine wva_default_start() bind(c, name="wva_default_start")
    call DefaultStart()
  end subroutine wva_default_start

  subroutine wva_default_initialize() bind(c, name="wva_default_initialize")
    call DefaultInitialize()
  end subroutine wva_default_initialize

  function wva_nof_active() result(n) bind(c, name="wva_nof_active")
    integer(c_int) :: n
    n = GetNOFActive()
  end function wva_nof_active

  subroutine wva_begin_active(t, n, nd, nb) bind(c, name="wva_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, nb
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
    nb = GetElementNOFBDOFs()
  end subroutine wva_begin_active

  function wva_nof_boundary() result(n) bind(c, name="wva_nof_boundary")
    integer(c_int) :: n
    n = GetNOFBoundaryElements()
  end function wva_nof_boundary

  subroutine wva_begin_boundary(t, n, nd, active) bind(c, name="wva_begin_boundary")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, active
    acc_elem => GetBoundaryElement(t)
    active = merge(1, 0, ActiveBoundaryElement())
    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
  end subroutine wva_begin_boundary

  subroutine wva_default_finish_bulk() bind(c, name="wva_default_finish_bulk")
    call DefaultFinishBulkAssembly()
  end subroutine wva_default_finish_bulk

  subroutine wva_default_finish_bnd() bind(c, name="wva_default_finish_bnd")
    call DefaultFinishBoundaryAssembly()
  end subroutine wva_default_finish_bnd

  subroutine wva_default_finish_assem() bind(c, name="wva_default_finish_assem")
    call DefaultFinishAssembly()
  end subroutine wva_default_finish_assem

  subroutine wva_default_dirichlet() bind(c, name="wva_default_dirichlet")
    call DefaultDirichletBCs()
  end subroutine wva_default_dirichlet

  function wva_default_solve() result(norm) bind(c, name="wva_default_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function wva_default_solve

  subroutine wva_default_finish() bind(c, name="wva_default_finish")
    call DefaultFinish()
  end subroutine wva_default_finish

  ! -- LocalMatrix context ----------------------------------------------
  subroutine wva_lm_init(ph) bind(c, name="wva_lm_init")
    integer(c_int) :: ph
    if (ph == 0) return
    call ListInitElementKeyword(acc_lm_load, 'Body Force', 'Sound source')
    call ListInitElementKeyword(acc_lm_ss, 'Material', 'Sound speed', &
        UnfoundFatal=.true.)
    call ListInitElementKeyword(acc_lm_damp, 'Material', 'Sound damping')
    call ListInitElementKeyword(acc_lm_rho, 'Material', 'Density')
    call ListInitElementKeyword(acc_lm_react, 'Material', &
        'Sound reaction damping')
    ph = 0
  end subroutine wva_lm_init

  subroutine wva_lm_nodes() bind(c, name="wva_lm_nodes")
    call GetElementNodes(acc_nodes_lm)
  end subroutine wva_lm_nodes

  function wva_lm_ip_begin() result(n) bind(c, name="wva_lm_ip_begin")
    integer(c_int) :: n
    acc_ip_lm = GaussPoints(acc_elem)
    n = acc_ip_lm % n
  end function wva_lm_ip_begin

  subroutine wva_lm_ip_point(t, u, v, w, s) bind(c, name="wva_lm_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip_lm % u(t)
    v = acc_ip_lm % v(t)
    w = acc_ip_lm % w(t)
    s = acc_ip_lm % s(t)
  end subroutine wva_lm_ip_point

  subroutine wva_lm_eleminfo(u, v, w, nd, basis, dbasis, detj, stat) &
      bind(c, name="wva_lm_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: nd
    real(c_double) :: basis(nd), dbasis(nd, 3), detj
    integer(c_int) :: stat
    logical :: ok
    ! Whole-array travel, exactly like the original internal procedure
    ! (Basis(nd), dBasisdx(nd,3) automatics passed whole).
    ok = ElementInfo(acc_elem, acc_nodes_lm, u, v, w, detj, &
        basis, dbasis)
    stat = merge(1, 0, ok)
  end subroutine wva_lm_eleminfo

  subroutine wva_lm_ip(nd, ld, basis, dbasisdx, sdet, stiff, damp, mass, force) &
      bind(c, name="wva_lm_ip")
    ! `stiff`/`damp`/`mass` ride on the caller's `ld`-strided scratch (the
    ! Rust-side analogue of the original's per-call automatics); only the
    ! leading `(1:nd,1:nd)`/`(1:nd)` sections are touched, exactly as the
    ! original touches its whole `nd`-sized automatics.
    integer(c_int), value :: nd, ld
    real(c_double) :: basis(nd), dbasisdx(nd, 3)
    real(c_double) :: stiff(ld, ld), damp(ld, ld), mass(ld, ld), force(ld)
    real(c_double), value :: sdet
    real(kind=dp) :: c, att, react, loadatip, rho, weight
    logical :: dampingactive, reactivemedium, assemblesource, found
    integer :: p, q
    ! Exact transcription of the per-IP statement sequence (same compiler,
    ! same evaluation order as the original internal procedure).
    c = ListGetElementReal(acc_lm_ss, basis, acc_elem)
    att = ListGetElementReal(acc_lm_damp, basis, acc_elem, dampingactive)
    react = ListGetElementReal(acc_lm_react, basis, acc_elem, reactivemedium)
    assemblesource = .false.
    if (.not. acc_lm_load % notpresentanywhere) &
        loadatip = ListGetElementReal(acc_lm_load, basis, acc_elem, assemblesource)
    weight = sdet
    rho = ListGetElementReal(acc_lm_rho, basis, acc_elem, found)
    if (found) weight = weight / rho
    stiff(1:nd, 1:nd) = stiff(1:nd, 1:nd) + weight * &
        matmul(dbasisdx, transpose(dbasisdx))
    if (dampingactive) damp(1:nd, 1:nd) = damp(1:nd, 1:nd) + weight * &
        (att / c**2) * matmul(dbasisdx, transpose(dbasisdx))
    do p = 1, nd
      do q = 1, nd
        if (reactivemedium) damp(p, q) = damp(p, q) + weight * &
            (react / c**2) * basis(q) * basis(p)
        mass(p, q) = mass(p, q) + weight * &
            (1.0_dp / c**2) * basis(q) * basis(p)
      end do
    end do
    if (assemblesource) &
        force(1:nd) = force(1:nd) - weight * loadatip * basis(1:nd)
  end subroutine wva_lm_ip

  subroutine wva_lm_tail(nd, nb, ld, mass, damp, stiff, force, transient, needmass) &
      bind(c, name="wva_lm_tail")
    integer(c_int), value :: nd, nb, ld, transient, needmass
    real(c_double) :: mass(ld, ld), damp(ld, ld), stiff(ld, ld), force(ld)
    if (transient /= 0 .or. needmass /= 0) then
      if (transient /= 0) then
        call Default2ndOrderTime(mass(1:nd, 1:nd), damp(1:nd, 1:nd), &
            stiff(1:nd, 1:nd), force(1:nd))
      else
        call DefaultUpdateMass(mass(1:nd, 1:nd))
        call DefaultUpdateDamp(damp(1:nd, 1:nd))
      end if
    end if
    call CondensateP(nd - nb, nb, stiff(1:ld, 1:ld), force(1:ld))
    call DefaultUpdateEquations(stiff(1:ld, 1:ld), force(1:ld))
  end subroutine wva_lm_tail

  ! -- LocalMatrixBC context ---------------------------------------------
  subroutine wva_lbc_head(ph, outflow, proceed) bind(c, name="wva_lbc_head")
    integer(c_int) :: ph, outflow, proceed
    type(valuelist_t), pointer :: bc
    logical :: obc, found
    bc => GetBC()
    if (.not. associated(bc)) then
      outflow = 0
      proceed = 0
      return
    end if
    obc = GetLogical(bc, 'Plane Wave BC', found)
    if (.not. found) obc = GetLogical(bc, 'Outflow Boundary', found)
    outflow = merge(1, 0, obc)
    if (ph /= 0) call ListInitElementKeyword(acc_bc_flux, &
        'Boundary Condition', 'Source Acceleration')
    if ((.not. obc) .and. acc_bc_flux % notpresentanywhere) then
      proceed = 0
      return
    end if
    if (ph /= 0) then
      call ListInitElementKeyword(acc_bc_ss, 'Material', 'Sound speed', &
          UnfoundFatal=.true.)
      call ListInitElementKeyword(acc_bc_rho, 'Material', 'Density', &
          UnfoundFatal=.true.)
      ph = 0
    end if
    proceed = 1
  end subroutine wva_lbc_head

  subroutine wva_lbc_nodes() bind(c, name="wva_lbc_nodes")
    call GetElementNodes(acc_nodes_lbc)
  end subroutine wva_lbc_nodes

  function wva_lbc_ip_begin() result(n) bind(c, name="wva_lbc_ip_begin")
    integer(c_int) :: n
    acc_ip_lbc = GaussPoints(acc_elem)
    n = acc_ip_lbc % n
  end function wva_lbc_ip_begin

  subroutine wva_lbc_ip_point(t, u, v, w, s) bind(c, name="wva_lbc_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip_lbc % u(t)
    v = acc_ip_lbc % v(t)
    w = acc_ip_lbc % w(t)
    s = acc_ip_lbc % s(t)
  end subroutine wva_lbc_ip_point

  subroutine wva_lbc_eleminfo(u, v, w, nd, basis, dbasis, detj, stat) &
      bind(c, name="wva_lbc_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: nd
    real(c_double) :: basis(nd), dbasis(nd, 3), detj
    integer(c_int) :: stat
    logical :: ok
    ! Whole-array travel, exactly like the original internal procedure.
    ok = ElementInfo(acc_elem, acc_nodes_lbc, u, v, w, detj, &
        basis, dbasis)
    stat = merge(1, 0, ok)
  end subroutine wva_lbc_eleminfo

  subroutine wva_lbc_ip(nd, ld, basis, sdet, damp, force, outflow) &
      bind(c, name="wva_lbc_ip")
    ! `damp` rides on the caller's `ld`-strided scratch; only the leading
    ! `(1:nd,1:nd)`/`(1:nd)` sections are touched, exactly as the original
    ! touches its whole `nd`-sized automatics.
    integer(c_int), value :: nd, ld, outflow
    real(c_double) :: basis(nd), damp(ld, ld), force(ld)
    real(c_double), value :: sdet
    real(kind=dp) :: c, g, rho, weight
    logical :: assembleflux, found
    integer :: p, q
    ! Exact transcription of the per-IP statement sequence.
    weight = sdet
    rho = ListGetElementRealParent(acc_bc_rho, basis, acc_elem, found)
    if (found) weight = weight / rho
    if (outflow /= 0) then
      c = ListGetElementRealParent(acc_bc_ss, basis, acc_elem)
      do p = 1, nd
        do q = 1, nd
          damp(p, q) = damp(p, q) + weight * basis(q) * basis(p) / c
        end do
      end do
    end if
    g = ListGetElementReal(acc_bc_flux, basis, acc_elem, assembleflux)
    if (assembleflux) then
      force(1:nd) = force(1:nd) - weight * g * basis(1:nd)
    end if
  end subroutine wva_lbc_ip

  subroutine wva_lbc_tail(nd, ld, mass, damp, stiff, force, transient, needmass) &
      bind(c, name="wva_lbc_tail")
    integer(c_int), value :: nd, ld, transient, needmass
    real(c_double) :: mass(ld, ld), damp(ld, ld), stiff(ld, ld), force(ld)
    if (transient /= 0 .or. needmass /= 0) then
      if (transient /= 0) then
        call Default2ndOrderTime(mass(1:nd, 1:nd), damp(1:nd, 1:nd), &
            stiff(1:nd, 1:nd), force(1:nd))
      else
        call DefaultUpdateMass(mass(1:nd, 1:nd))
        call DefaultUpdateDamp(damp(1:nd, 1:nd))
      end if
    end if
    call DefaultUpdateEquations(stiff(1:ld, 1:ld), force(1:ld))
  end subroutine wva_lbc_tail

end module wavesolver_acc
