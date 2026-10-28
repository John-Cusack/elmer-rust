! T021-A DivergenceSolver port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/DivergenceSolver.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments. All Fortran-object
! traffic — TYPE(...) handles, OPTIONAL/POINTER dummies, hidden lengths,
! list-directed timing reports — stays inside Fortran, so the Rust side is
! pure C calls. The wrapper module mirrors DivergenceSolver's own locals
! (current element, Nodes/IP SAVE state, entry-level flags) with
! identical lifetimes. All floating-point work (the Galerkin assembly,
! `**`, `SQRT`, `MIN`/`MAX` reductions) and the `RealTime`-based reports
! stay in Fortran (same compiler), so bit-identity holds by construction
! wherever the call sequence matches. The `Divergence Coefficient`
! keyword is read in the original but used only by a commented-out line;
! the read is pure (no side effects) and is documented as skipped rather
! than replicated.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module divergencesolver_acc
  use iso_c_binding
  use coordinatesystems
  use defutils
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(valuelist_t), pointer, save :: acc_params => null()
  type(element_t), pointer, save :: acc_elem => null()
  ! BulkAssembly context (mirrors its SAVE Nodes; IP is a fresh local
  ! there, kept module-side with identical per-element lifetime).
  type(nodes_t), target, save :: acc_nodes
  type(gaussintegrationpoints_t), target, save :: acc_ip
  ! Entry-level state (mirrors the solver body's locals).
  integer(c_int), save :: acc_dim = 0
  character(len=MAX_NAME_LEN), save :: acc_varname = ''
  logical, save :: acc_relative = .false.
  logical, save :: acc_csym = .false.
  logical, save :: acc_constbulk = .false.
  logical, save :: acc_inuse = .false.
  real(c_double), save :: acc_at0 = 0.0_c_double

contains

  ! -- entry state ----------------------------------------------------
  subroutine dvs_begin(h_solver) bind(c, name="dvs_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine dvs_begin

  subroutine dvs_info_start() bind(c, name="dvs_info_start")
    call Info('DivergenceSolver', '-------------------------------------', level=4)
    call Info('DivergenceSolver', 'Computing the divergence field', level=4)
    call Info('DivergenceSolver', '-------------------------------------', level=4)
  end subroutine dvs_info_start

  function dvs_var_check() result(code) bind(c, name="dvs_var_check")
    integer(c_int) :: code
    type(variable_t), pointer :: ds
    ds => acc_solver % variable
    if (.not. associated(ds)) then
      code = 4
      return
    end if
    if (count(ds % perm > 0) <= 0) then
      code = 1
      return
    end if
    if (.not. associated(acc_solver % matrix)) then
      code = 2
      return
    end if
    if (ds % dofs /= 1) then
      code = 3
      return
    end if
    code = 0
  end function dvs_var_check

  subroutine dvs_warn_zero() bind(c, name="dvs_warn_zero")
    call Warn('DivergenceSolver', 'Size is zero, nothing to compute')
  end subroutine dvs_warn_zero

  subroutine dvs_warn_nomatrix() bind(c, name="dvs_warn_nomatrix")
    call Warn('DivergenceSolver', 'No matrix associated, exiting')
  end subroutine dvs_warn_nomatrix

  subroutine dvs_fatal_dofs() bind(c, name="dvs_fatal_dofs")
    call Fatal('DivergenceSolver', 'Divergence should have 1 component')
  end subroutine dvs_fatal_dofs

  subroutine dvs_fatal_novar() bind(c, name="dvs_fatal_novar")
    call Fatal('DivergenceSolver', 'Variable does not exist!')
  end subroutine dvs_fatal_novar

  subroutine dvs_fatal_alloc() bind(c, name="dvs_fatal_alloc")
    call Fatal('DivergenceSolver', 'Memory allocation error.')
  end subroutine dvs_fatal_alloc

  ! -- sweep state ------------------------------------------------------
  subroutine dvs_sweep_begin() bind(c, name="dvs_sweep_begin")
    logical :: gotit, constbulk
    acc_dim = CoordinateSystemDimension()
    acc_params => GetSolverParams()
    acc_varname = GetString(acc_params, 'Divergence Variable', gotit)
    if (.not. gotit) acc_varname = GetString(acc_params, 'Target Variable', gotit)
    if (.not. gotit) acc_varname = 'Velocity'
    acc_csym = CurrentCoordinateSystem() == AxisSymmetric .or. &
        CurrentCoordinateSystem() == CylindricSymmetric
    acc_at0 = RealTime()
    constbulk = GetLogical(acc_params, 'Constant Bulk Matrix', gotit)
    acc_constbulk = constbulk
    acc_inuse = constbulk .and. associated(acc_solver % matrix % bulkvalues)
    acc_relative = ListGetLogical(acc_params, 'Relative Divergence', gotit)
    call DefaultInitialize(acc_solver, acc_inuse)
  end subroutine dvs_sweep_begin

  function dvs_alloc_n() result(n) bind(c, name="dvs_alloc_n")
    integer(c_int) :: n
    n = max(acc_solver % mesh % maxelementdofs, &
        acc_solver % mesh % maxelementnodes)
  end function dvs_alloc_n

  function dvs_nof_active_elems() result(n) bind(c, name="dvs_nof_active_elems")
    integer(c_int) :: n
    n = acc_solver % numberofactiveelements
  end function dvs_nof_active_elems

  ! -- BulkAssembly element body (exact transcription) ------------------
  subroutine dvs_elem(elem, nmax, stiff, force, vx, vy, vz, basis, dbasis) &
      bind(c, name="dvs_elem")
    integer(c_int), value :: elem, nmax
    real(c_double) :: stiff(nmax, nmax), force(nmax)
    real(c_double) :: vx(nmax), vy(nmax), vz(nmax)
    real(c_double) :: basis(nmax), dbasis(nmax, 3)
    integer :: n, nd, t, p, q
    real(c_double) :: weight, detj, source, x, vabs, velo(3)
    logical :: found
    acc_elem => GetActiveElement(elem)
    call GetElementNodes(acc_nodes)
    nd = GetElementNOFDOFs()
    n = GetElementNOFNodes()
    call GetScalarLocalSolution(vx, ComponentName(acc_varname, 1))
    call GetScalarLocalSolution(vy, ComponentName(acc_varname, 2))
    if (acc_dim == 3) call GetScalarLocalSolution(vz, ComponentName(acc_varname, 3))
    if (acc_relative) then
      velo(1) = sum(vx(1:nd)) / nd
      velo(2) = sum(vy(1:nd)) / nd
      if (acc_dim == 3) velo(3) = sum(vz(1:nd)) / nd
      vabs = max(sqrt(sum(velo(1:acc_dim)**2)), epsilon(vabs))
    end if
    acc_ip = GaussPoints(acc_elem)
    stiff = 0.0_c_double
    force = 0.0_c_double
    do t = 1, acc_ip % n
      found = ElementInfo(acc_elem, acc_nodes, acc_ip % u(t), &
          acc_ip % v(t), acc_ip % w(t), detj, basis, dbasis)
      weight = acc_ip % s(t) * detj
      if (acc_csym) then
        x = sum(basis(1:n) * acc_nodes % x(1:n))
        weight = weight * x
      end if
      if (.not. acc_inuse) then
        do p = 1, nd
          do q = 1, nd
            stiff(p, q) = stiff(p, q) + weight * basis(q) * basis(p)
          end do
        end do
      end if
      source = sum(dbasis(1:nd, 1) * vx(1:nd))
      source = source + sum(dbasis(1:nd, 2) * vy(1:nd))
      if (acc_dim == 3) source = source + sum(dbasis(1:nd, 3) * vz(1:nd))
      if (acc_csym) then
        source = source + sum(basis(1:nd) * vx(1:nd)) / x
      end if
      if (acc_relative) source = source / vabs
      force(1:nd) = force(1:nd) + basis(1:nd) * weight * source
    end do
    if (.not. acc_inuse) then
      call DefaultUpdateEquations(stiff, force(1:nd))
    else
      call DefaultUpdateForce(force(1:nd))
    end if
  end subroutine dvs_elem

  ! -- shared assembly tail ---------------------------------------------
  subroutine dvs_default_finish_bulk() bind(c, name="dvs_default_finish_bulk")
    if (acc_constbulk) then
      call DefaultFinishBulkAssembly(bulkupdate = .not. acc_inuse, rhsupdate = .false.)
    else
      call DefaultFinishBulkAssembly()
    end if
  end subroutine dvs_default_finish_bulk

  subroutine dvs_default_finish() bind(c, name="dvs_default_finish")
    call DefaultFinishAssembly()
  end subroutine dvs_default_finish

  subroutine dvs_default_dirichlet() bind(c, name="dvs_default_dirichlet")
    call DefaultDirichletBCs()
  end subroutine dvs_default_dirichlet

  function dvs_tail() result(norm) bind(c, name="dvs_tail")
    real(c_double) :: norm, at1, at2
    at1 = RealTime()
    write(message, *) 'Assembly Time: ', at1 - acc_at0
    call Info('DivergenceSolver', message, level=5)
    norm = DefaultSolve()
    at2 = RealTime()
    write(message, *) 'Solution Time: ', at2 - at1
    call Info('DivergenceSolver', message, level=5)
    write(message, *) 'Result Norm: ', norm
    call Info('DivergenceSolver', message, level=4)
  end function dvs_tail

  ! -- DivergenceSolver_Init body -----------------------------------------
  subroutine dvs_init_begin(h_solver) bind(c, name="dvs_init_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine dvs_init_begin

  subroutine dvs_init_body() bind(c, name="dvs_init_body")
    type(valuelist_t), pointer :: solverpars
    integer :: dim
    character(len=MAX_NAME_LEN) :: varname, divname
    logical :: gotit
    solverpars => GetSolverParams()
    dim = CoordinateSystemDimension()
    varname = GetString(solverpars, 'Divergence Variable', gotit)
    if (gotit) then
      divname = 'Div ' // trim(varname)
    else
      divname = 'Divergence'
    end if
    if (.not. ListCheckPresent(solverpars, 'Variable')) then
      call ListAddInteger(solverpars, 'Variable DOFs', 1)
      call ListAddString(solverpars, 'Variable', divname)
    end if
    call ListAddInteger(solverpars, 'Time derivative order', 0)
    call ListAddNewString(solverpars, 'Linear System Solver', 'Iterative')
    call ListAddNewString(solverpars, 'Linear System Iterative Method', 'cg')
    call ListAddNewString(solverpars, 'Linear System Preconditioning', 'ILU0')
    call ListAddNewInteger(solverpars, 'Linear System Max Iterations', 500)
    call ListAddNewConstReal(solverpars, &
        'Linear System Convergence Tolerance', 1.0e-8_dp)
  end subroutine dvs_init_body

end module divergencesolver_acc
