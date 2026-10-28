! T021-A ElementSizeSolver port: statement-level BIND(C) wrappers for
! the Elmer calls in fem/src/modules/ElementSizeSolver.F90, compiled
! against the built .mod files of the target configure (A2 shipped
! mechanism, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments. All Fortran-object
! traffic — TYPE(...) handles, OPTIONAL/POINTER dummies, assumed-shape
! descriptors, hidden lengths, ES12.4 message formatting — stays inside
! Fortran, so the Rust side is pure C calls. The wrapper module mirrors
! ElementSizeSolver's own locals (current element, Nodes SAVE state, IP
! shallow copy, Power/NoWeight/ElemMin/ElemMax entry state) with
! identical lifetimes. All floating-point updates, MIN/MAX reductions,
! `**` powers and formatted reports stay in Fortran (same compiler),
! so bit-identity holds by construction wherever the call sequence
! matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module elementsizesolver_acc
  use iso_c_binding
  use defutils
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(valuelist_t), pointer, save :: acc_params => null()
  type(element_t), pointer, save :: acc_elem => null()
  ! LocalMatrix context (mirrors its SAVE Nodes and IP local).
  type(nodes_t), target, save :: acc_nodes
  type(gaussintegrationpoints_t), save :: acc_ip
  ! Entry-level state (mirrors the solver body's locals).
  real(c_double), save :: acc_power = 1.0_c_double
  logical, save :: acc_noweight = .false.
  real(c_double), save :: acc_elemmin = 0.0_c_double
  real(c_double), save :: acc_elemmax = 0.0_c_double

contains

  ! -- entry state ----------------------------------------------------
  subroutine ees_begin(h_solver) bind(c, name="ees_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine ees_begin

  subroutine ees_info_start() bind(c, name="ees_info_start")
    call Info('ElementSizeSolver', 'Computing nodal element size indicator')
  end subroutine ees_info_start

  function ees_max_elem_nodes() result(n) bind(c, name="ees_max_elem_nodes")
    integer(c_int) :: n
    type(mesh_t), pointer :: m
    m => GetMesh()
    n = m % maxelementnodes
  end function ees_max_elem_nodes

  subroutine ees_params_begin() bind(c, name="ees_params_begin")
    acc_params => GetSolverParams()
  end subroutine ees_params_begin

  subroutine ees_read_exponent(power, gotit) bind(c, name="ees_read_exponent")
    real(c_double) :: power
    integer(c_int) :: gotit
    logical :: found
    power = ListGetCReal(acc_params, 'Element Size Exponent', found)
    gotit = merge(1, 0, found)
  end subroutine ees_read_exponent

  subroutine ees_read_noweight(val, gotit) bind(c, name="ees_read_noweight")
    integer(c_int) :: val, gotit
    logical :: found
    val = merge(1, 0, ListGetLogical(acc_params, 'No Integration Weight', found))
    gotit = merge(1, 0, found)
  end subroutine ees_read_noweight

  ! -- sweep state ------------------------------------------------------
  subroutine ees_sweep_begin(power, noweigh) bind(c, name="ees_sweep_begin")
    real(c_double), value :: power
    integer(c_int), value :: noweigh
    acc_power = power
    acc_noweight = (noweigh /= 0)
    acc_elemmin = huge(acc_elemmin)
    acc_elemmax = -huge(acc_elemmax)
  end subroutine ees_sweep_begin

  ! -- assembly driver --------------------------------------------------
  function ees_nof_active() result(n) bind(c, name="ees_nof_active")
    integer(c_int) :: n
    n = GetNOFActive()
  end function ees_nof_active

  subroutine ees_default_initialize() bind(c, name="ees_default_initialize")
    call DefaultInitialize()
  end subroutine ees_default_initialize

  subroutine ees_begin_active(t, n) bind(c, name="ees_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes()
  end subroutine ees_begin_active

  ! -- LocalMatrix context ----------------------------------------------
  subroutine ees_lm_begin() bind(c, name="ees_lm_begin")
    call GetElementNodes(acc_nodes)
    acc_ip = GaussPoints(acc_elem)
  end subroutine ees_lm_begin

  function ees_lm_ip_n() result(n) bind(c, name="ees_lm_ip_n")
    integer(c_int) :: n
    n = acc_ip % n
  end function ees_lm_ip_n

  subroutine ees_lm_ip(t, n, stiff, force) bind(c, name="ees_lm_ip")
    integer(c_int), value :: t, n
    real(c_double) :: stiff(n, n), force(n)
    ! Exact transcription of the LocalMatrix integration-point body.
    real(c_double) :: basis(n), detj, loadatip, weight
    logical :: stat
    integer :: i, j
    stat = ElementInfo(acc_elem, acc_nodes, acc_ip % u(t), &
        acc_ip % v(t), acc_ip % w(t), detj, basis)
    loadatip = detj ** acc_power
    if (acc_noweight) then
      weight = acc_ip % s(t)
    else
      weight = acc_ip % s(t) * detj
    end if
    acc_elemmin = min(acc_elemmin, loadatip)
    acc_elemmax = max(acc_elemmax, loadatip)
    do i = 1, n
      do j = 1, n
        stiff(i, j) = stiff(i, j) + weight * basis(i) * basis(j)
      end do
      force(i) = force(i) + weight * basis(i) * loadatip
    end do
  end subroutine ees_lm_ip

  ! -- shared assembly tail ---------------------------------------------
  subroutine ees_default_update(nmax, g, f) bind(c, name="ees_default_update")
    integer(c_int), value :: nmax
    real(c_double) :: g(nmax, nmax), f(nmax)
    call DefaultUpdateEquations(g, f)
  end subroutine ees_default_update

  subroutine ees_default_finish_bulk() bind(c, name="ees_default_finish_bulk")
    call DefaultFinishBulkAssembly()
  end subroutine ees_default_finish_bulk

  subroutine ees_default_finish() bind(c, name="ees_default_finish")
    call DefaultFinishAssembly()
  end subroutine ees_default_finish

  subroutine ees_default_dirichlet() bind(c, name="ees_default_dirichlet")
    call DefaultDirichletBCs()
  end subroutine ees_default_dirichlet

  function ees_default_solve() result(norm) bind(c, name="ees_default_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function ees_default_solve

  subroutine ees_report() bind(c, name="ees_report")
    write(message, '(A,ES12.4)') 'Minimum Element Size: ', acc_elemmin
    call Info('ElementSizeSolver', message)
    write(message, '(A,ES12.4)') 'Maximum Element Size: ', acc_elemmax
    call Info('ElementSizeSolver', message)
    write(message, '(A,ES12.4)') 'Element Size Ratio: ', &
        acc_elemmax / acc_elemmin
    call Info('ElementSizeSolver', message)
  end subroutine ees_report

  subroutine ees_fatal_alloc() bind(c, name="ees_fatal_alloc")
    call Fatal('PoissonSolve', 'Memory allocation error.')
  end subroutine ees_fatal_alloc

end module elementsizesolver_acc
