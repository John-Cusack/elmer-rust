! T021-B ModelPDEhandle port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/ModelPDEhandle.F90, compiled against the
! built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module mirrors ModelPDEhandle's own locals (current element, Nodes SAVE
! state per internal procedure, IP shallow copies, the twelve
! ValueHandle_t element keywords, the dim snapshot, and the
! host-associated TotArea/TotLen/TotSrc accumulators) with identical
! lifetimes.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module modelpdehandle_acc
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
  type(valuelist_t), pointer, save :: acc_bc => null()
  ! LocalMatrix context (mirrors its SAVE Nodes and IP local).
  type(nodes_t), target, save :: acc_nodes_lm
  type(gaussintegrationpoints_t), save :: acc_ip_lm
  ! LocalMatrixBC context (mirrors its SAVE Nodes and IP local).
  type(nodes_t), target, save :: acc_nodes_bc
  type(gaussintegrationpoints_t), save :: acc_ip_bc
  ! Element-keyword handles (entry locals, initialized once per entry).
  ! (The original's FieldSource_h is declared but never initialized, so
  ! it is omitted here entirely.)
  type(valuehandle_t), save :: acc_load_h, acc_diff_h, acc_react_h, &
      acc_conv_h, acc_time_h, acc_vel1_h, acc_vel2_h, acc_vel3_h, &
      acc_bcflux_h, acc_bccoeff_h, acc_bcext_h
  ! Host-associated accumulators and the dim snapshot.
  real(c_double), save :: acc_totarea = 0.0_c_double
  real(c_double), save :: acc_totlen = 0.0_c_double
  real(c_double), save :: acc_totsrc = 0.0_c_double
  integer(c_int), save :: acc_dim = 0

contains

  ! -- entry state ----------------------------------------------------
  subroutine mdh_begin(h_solver) bind(c, name="mdh_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine mdh_begin

  subroutine mdh_init_handles() bind(c, name="mdh_init_handles")
    call ListInitElementKeyword(acc_load_h, 'Body Force', 'Field Source')
    call ListInitElementKeyword(acc_diff_h, 'Material', 'Diffusion Coefficient')
    call ListInitElementKeyword(acc_react_h, 'Material', 'Reaction Coefficient')
    call ListInitElementKeyword(acc_conv_h, 'Material', 'Convection Coefficient')
    call ListInitElementKeyword(acc_time_h, 'Material', 'Time Derivative Coefficient')
    call ListInitElementKeyword(acc_vel1_h, 'Material', 'Convection Velocity 1')
    call ListInitElementKeyword(acc_vel2_h, 'Material', 'Convection Velocity 2')
    call ListInitElementKeyword(acc_vel3_h, 'Material', 'Convection Velocity 3')
    call ListInitElementKeyword(acc_bcflux_h, 'Boundary Condition', 'Field Flux')
    call ListInitElementKeyword(acc_bccoeff_h, 'Boundary Condition', 'Robin Coefficient')
    call ListInitElementKeyword(acc_bcext_h, 'Boundary Condition', 'External Field')
  end subroutine mdh_init_handles

  subroutine mdh_default_start() bind(c, name="mdh_default_start")
    call DefaultStart()
  end subroutine mdh_default_start

  subroutine mdh_default_initialize() bind(c, name="mdh_default_initialize")
    call DefaultInitialize()
  end subroutine mdh_default_initialize

  function mdh_maxiter() result(m) bind(c, name="mdh_maxiter")
    integer(c_int) :: m
    logical :: found
    m = ListGetInteger(GetSolverParams(), &
        'Nonlinear System Max Iterations', found, minv=1)
    if (.not. found) m = 1
  end function mdh_maxiter

  function mdh_coord_dim() result(d) bind(c, name="mdh_coord_dim")
    integer(c_int) :: d
    d = CoordinateSystemDimension()
    acc_dim = d
  end function mdh_coord_dim

  subroutine mdh_reset_totals() bind(c, name="mdh_reset_totals")
    acc_totarea = 0.0_c_double
    acc_totlen = 0.0_c_double
    acc_totsrc = 0.0_c_double
  end subroutine mdh_reset_totals

  subroutine mdh_report_integ() bind(c, name="mdh_report_integ")
    logical :: found
    if (ListGetLogical(GetSolverParams(), 'CutFEM', found) &
        .or. ListGetLogical(GetSolverParams(), 'Integ Test', found)) then
      call ListAddConstReal(CurrentModel % Simulation, &
          'res: integ total area', acc_totarea)
      call ListAddConstReal(CurrentModel % Simulation, &
          'res: integ total len', acc_totlen)
      call ListAddConstReal(CurrentModel % Simulation, &
          'res: integ total src', acc_totsrc)
    end if
  end subroutine mdh_report_integ

  subroutine mdh_default_finish() bind(c, name="mdh_default_finish")
    call DefaultFinish()
  end subroutine mdh_default_finish

  ! -- assembly driver -------------------------------------------------
  function mdh_nof_active() result(n) bind(c, name="mdh_nof_active")
    integer(c_int) :: n
    n = GetNOFActive()
  end function mdh_nof_active

  subroutine mdh_begin_active(t, n, nd, nb) bind(c, name="mdh_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, nb
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
    nb = GetElementNOFBDOFs()
  end subroutine mdh_begin_active

  function mdh_nof_bnd() result(n) bind(c, name="mdh_nof_bnd")
    integer(c_int) :: n
    n = GetNOFBoundaryElements()
  end function mdh_nof_bnd

  subroutine mdh_begin_boundary(t, n, nd, nb, active) bind(c, name="mdh_begin_boundary")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, nb, active
    acc_elem => GetBoundaryElement(t)
    active = merge(1, 0, ActiveBoundaryElement())
    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
    nb = GetElementNOFBDOFs()
  end subroutine mdh_begin_boundary

  function mdh_default_cutfem() result(r) bind(c, name="mdh_default_cutfem")
    integer(c_int) :: r
    r = merge(1, 0, DefaultCutFEM())
  end function mdh_default_cutfem

  subroutine mdh_default_finish_bulk() bind(c, name="mdh_default_finish_bulk")
    call DefaultFinishBulkAssembly()
  end subroutine mdh_default_finish_bulk

  subroutine mdh_default_finish_bnd() bind(c, name="mdh_default_finish_bnd")
    call DefaultFinishBoundaryAssembly()
  end subroutine mdh_default_finish_bnd

  subroutine mdh_default_finish_asm() bind(c, name="mdh_default_finish_asm")
    call DefaultFinishAssembly()
  end subroutine mdh_default_finish_asm

  subroutine mdh_default_dirichlet() bind(c, name="mdh_default_dirichlet")
    call DefaultDirichletBCs()
  end subroutine mdh_default_dirichlet

  function mdh_default_solve() result(norm) bind(c, name="mdh_default_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function mdh_default_solve

  function mdh_default_converged() result(r) bind(c, name="mdh_default_converged")
    integer(c_int) :: r
    r = merge(1, 0, DefaultConverged())
  end function mdh_default_converged

  ! -- LocalMatrix context ----------------------------------------------
  subroutine mdh_lm_nodes() bind(c, name="mdh_lm_nodes")
    call GetElementNodes(acc_nodes_lm)
  end subroutine mdh_lm_nodes

  function mdh_lm_ip_begin() result(n) bind(c, name="mdh_lm_ip_begin")
    integer(c_int) :: n
    acc_ip_lm = GaussPointsAdapt(acc_elem)
    n = acc_ip_lm % n
  end function mdh_lm_ip_begin
  subroutine mdh_lm_ip_point(t, u, v, w, s) bind(c, name="mdh_lm_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip_lm % u(t)
    v = acc_ip_lm % v(t)
    w = acc_ip_lm % w(t)
    s = acc_ip_lm % s(t)
  end subroutine mdh_lm_ip_point


  subroutine mdh_lm_eleminfo(u, v, w, n, nd, basis, dbasis, detj, stat) &
      bind(c, name="mdh_lm_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: n, nd
    real(c_double) :: basis(nd), dbasis(nd, 3), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_elem, acc_nodes_lm, u, v, w, detj, &
        basis(1:nd), dbasis(1:nd, 1:3))
    stat = merge(1, 0, ok)
  end subroutine mdh_lm_eleminfo

  subroutine mdh_lm_ip_values(n, basis, load, rho, a1, a2, a3, d, c, r) &
      bind(c, name="mdh_lm_ip_values")
    integer(c_int), value :: n
    real(c_double) :: basis(n)
    real(c_double) :: load, rho, a1, a2, a3, d, c, r
    logical :: found
    load = ListGetElementReal(acc_load_h, basis, acc_elem, found)
    rho = ListGetElementReal(acc_time_h, basis, acc_elem, found)
    a1 = ListGetElementReal(acc_vel1_h, basis, acc_elem, found)
    a2 = ListGetElementReal(acc_vel2_h, basis, acc_elem, found)
    if (acc_dim == 3) then
      a3 = ListGetElementReal(acc_vel3_h, basis, acc_elem, found)
    else
      a3 = 0.0_c_double
    end if
    d = ListGetElementReal(acc_diff_h, basis, acc_elem, found)
    c = ListGetElementReal(acc_conv_h, basis, acc_elem, found)
    r = ListGetElementReal(acc_react_h, basis, acc_elem, found)
  end subroutine mdh_lm_ip_values

  subroutine mdh_lm_assemble(nd, dim, dbasis, basis, weight, load, d, c, &
      r, rho, a1, a2, a3, mass, stiff, force) bind(c, name="mdh_lm_assemble")
    integer(c_int), value :: nd, dim
    real(c_double) :: dbasis(nd, 3), basis(nd)
    real(c_double), value :: weight, load, d, c, r, rho, a1, a2, a3
    real(c_double) :: mass(nd, nd), stiff(nd, nd), force(nd)
    real(c_double) :: a(3)
    integer :: p, q
    a(1) = a1
    a(2) = a2
    a(3) = a3
    ! Exact transcription of the update statements (same compiler, same
    ! evaluation order as the original internal procedure).
    stiff(1:nd, 1:nd) = stiff(1:nd, 1:nd) + weight * &
        d * matmul(dbasis, transpose(dbasis))
    do p = 1, nd
      do q = 1, nd
        stiff(p, q) = stiff(p, q) + weight * &
            c * sum(a(1:dim) * dbasis(q, 1:dim)) * basis(p)
        stiff(p, q) = stiff(p, q) + weight * r * basis(q) * basis(p)
        mass(p, q) = mass(p, q) + weight * rho * basis(q) * basis(p)
      end do
    end do
    force(1:nd) = force(1:nd) + weight * load * basis(1:nd)
    acc_totarea = acc_totarea + weight
    acc_totsrc = acc_totsrc + weight * load
  end subroutine mdh_lm_assemble

  ! -- shared assembly tail ---------------------------------------------
  subroutine mdh_time_update(nd, mass, stiff, force) bind(c, name="mdh_time_update")
    integer(c_int), value :: nd
    real(c_double) :: mass(nd, nd), stiff(nd, nd), force(nd)
    call Default1stOrderTime(mass, stiff, force)
  end subroutine mdh_time_update

  subroutine mdh_condensate(nd, nb, ld, k, f) bind(c, name="mdh_condensate")
    integer(c_int), value :: nd, nb, ld
    real(c_double) :: k(ld, ld), f(ld)
    call CondensateP(nd, nb, k(1:ld, 1:ld), f(1:ld))
  end subroutine mdh_condensate

  subroutine mdh_default_update(ld, g, f) bind(c, name="mdh_default_update")
    integer(c_int), value :: ld
    real(c_double) :: g(ld, ld), f(ld)
    call DefaultUpdateEquations(g(1:ld, 1:ld), f(1:ld))
  end subroutine mdh_default_update

  ! -- LocalMatrixBC context --------------------------------------------
  function mdh_bc_assoc() result(r) bind(c, name="mdh_bc_assoc")
    integer(c_int) :: r
    acc_bc => GetBC()
    r = merge(1, 0, associated(acc_bc))
  end function mdh_bc_assoc

  subroutine mdh_bc_nodes() bind(c, name="mdh_bc_nodes")
    call GetElementNodes(acc_nodes_bc)
  end subroutine mdh_bc_nodes

  function mdh_bc_ip_begin() result(n) bind(c, name="mdh_bc_ip_begin")
    integer(c_int) :: n
    acc_ip_bc = GaussPoints(acc_elem)
    n = acc_ip_bc % n
  end function mdh_bc_ip_begin

  subroutine mdh_bc_ip_point(t, u, v, w, s) bind(c, name="mdh_bc_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip_bc % u(t)
    v = acc_ip_bc % v(t)
    w = acc_ip_bc % w(t)
    s = acc_ip_bc % s(t)
  end subroutine mdh_bc_ip_point

  subroutine mdh_bc_eleminfo(u, v, w, n, basis, detj, stat) &
      bind(c, name="mdh_bc_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: n
    real(c_double) :: basis(n), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_elem, acc_nodes_bc, u, v, w, detj, basis(1:n))
    stat = merge(1, 0, ok)
  end subroutine mdh_bc_eleminfo

  subroutine mdh_bc_ip_values(n, basis, f, c, ext) bind(c, name="mdh_bc_ip_values")
    integer(c_int), value :: n
    real(c_double) :: basis(n)
    real(c_double) :: f, c, ext
    logical :: found
    f = ListGetElementReal(acc_bcflux_h, basis, acc_elem, found)
    c = ListGetElementReal(acc_bccoeff_h, basis, acc_elem, found)
    ext = ListGetElementReal(acc_bcext_h, basis, acc_elem, found)
  end subroutine mdh_bc_ip_values

  subroutine mdh_bc_assemble(nd, basis, weight, c, f, ext, stiff, force) &
      bind(c, name="mdh_bc_assemble")
    integer(c_int), value :: nd
    real(c_double) :: basis(nd)
    real(c_double), value :: weight, c, f, ext
    real(c_double) :: stiff(nd, nd), force(nd)
    integer :: p, q
    ! Exact transcription of the per-node update loop.
    do p = 1, nd
      do q = 1, nd
        stiff(p, q) = stiff(p, q) + weight * c * basis(q) * basis(p)
      end do
    end do
    force(1:nd) = force(1:nd) + weight * (f + c * ext) * basis(1:nd)
    acc_totlen = acc_totlen + weight
  end subroutine mdh_bc_assemble

end module modelpdehandle_acc
