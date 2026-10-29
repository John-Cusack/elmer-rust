! T022-A CyclicConvergence port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/CyclicConvergence.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! hidden lengths, formatted I/O — stays inside Fortran, so the Rust side
! is pure C calls. The wrapper module mirrors CyclicConvergence's own
! SAVE state (Visited, valTable, tsize, prevCycle, prevVal, prevDev,
! prevminv, prevmaxv, prevrangev, iCycle, ConvCycle, Converged, plus the
! internal procedure's SAVE ConvVisited) with identical lifetimes. All
! floating-point work (FLOOR, SUM, SQRT, MINVAL, MAXVAL, the error
! ratios) and every formatted WRITE stay in Fortran (same compiler), so
! bit-identity holds by construction wherever the call sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module cyclicconvergence_acc
  use iso_c_binding
  use defutils
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(mesh_t), pointer, save :: acc_mesh => null()
  type(valuelist_t), pointer, save :: acc_params => null()
  type(variable_t), pointer, save :: acc_percycle => null()

  ! Entry SAVE state (mirrors the original's SAVE statement exactly).
  logical, save :: acc_visited = .false.
  logical, save :: acc_converged = .false.
  real(c_double), pointer, save :: acc_valtable(:) => null()
  integer(c_int), save :: acc_tsize = 500
  integer(c_int), save :: acc_ncycle = 0
  integer(c_int), save :: acc_prevcycle = 0
  integer(c_int), save :: acc_icycle = 0
  integer(c_int), save :: acc_convcycle = 0
  real(c_double), save :: acc_prevval = 0.0_c_double
  real(c_double), save :: acc_prevdev = 0.0_c_double
  real(c_double), save :: acc_prevminv = 0.0_c_double
  real(c_double), save :: acc_prevmaxv = 0.0_c_double
  real(c_double), save :: acc_prevrangev = 0.0_c_double
  ! Cycle statistics (entry locals in the original, consumed by the
  ! internal procedure and the evaluate branch within one invocation).
  real(c_double), save :: acc_errval = 0.0_c_double
  real(c_double), save :: acc_errdev = 0.0_c_double
  real(c_double), save :: acc_errrange = 0.0_c_double
  real(c_double), save :: acc_meanval = 0.0_c_double
  real(c_double), save :: acc_meandev = 0.0_c_double
  ! Internal-procedure SAVE state.
  logical, save :: acc_convvisited = .false.

contains

  ! -- entry state ----------------------------------------------------
  subroutine cca_begin(h_model, h_solver) bind(c, name="cca_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    acc_mesh => acc_solver % mesh
    acc_params => acc_solver % values
  end subroutine cca_begin

  subroutine cca_info_start() bind(c, name="cca_info_start")
    call Info('CyclicConvergence', &
        'Checking for convergence of cyclic system')
  end subroutine cca_info_start

  subroutine cca_percycle(dt) bind(c, name="cca_percycle")
    real(c_double), value :: dt
    acc_percycle => VariableGet(acc_mesh % variables, 'Periodic Cycle')
    if (.not. associated(acc_percycle)) then
      call Fatal('CyclicConvergence', 'Could not define "Periodic Cycle"')
    end if
    ! We apply a small offset so that at the exact time we would
    ! swap the cycle is triggered just at the end of cycle.
    acc_ncycle = floor(acc_percycle % values(1) + 1.0e-3*dt)
  end subroutine cca_percycle

  subroutine cca_first_touch() bind(c, name="cca_first_touch")
    if (.not. acc_visited) then
      acc_prevcycle = acc_ncycle
      allocate(acc_valtable(acc_tsize))
      acc_valtable = 0.0_c_double
      acc_visited = .true.
    end if
  end subroutine cca_first_touch

  subroutine cca_resolve(val) bind(c, name="cca_resolve")
    real(c_double) :: val
    logical :: found
    type(variable_t), pointer :: pvar
    character(MAX_NAME_LEN) :: valstr
    val = ListGetCReal(acc_params, 'Convergence Value', found)
    if (.not. found) then
      valstr = ListGetString(acc_params, 'Convergence Value Name', found)
      if (.not. found) then
        call Fatal('CyclicConvergence', &
            'Give "Convergence Value" or "Convergence Value Name"!')
      end if
      val = ListGetCReal(acc_model % simulation, valstr, found)
      if (.not. found) then
        pvar => VariableGet(acc_mesh % variables, valstr)
        if (associated(pvar)) then
          val = pvar % values(1)
        else
          call Fatal('CyclicConvergence', 'Could not find ' // &
              '"Convergence Value Name" as variable or list entry!')
        end if
      end if
    end if
  end subroutine cca_resolve

  subroutine cca_push(val) bind(c, name="cca_push")
    real(c_double), value :: val
    real(c_double), pointer :: tmptable(:)
    ! Tabulate values
    acc_icycle = acc_icycle + 1
    if (acc_icycle > acc_tsize) then
      call Info('CyclicConvergence', &
          'Increasing size of valTable to '//I2S(2*acc_tsize), level=6)
      allocate(tmptable(2*acc_tsize))
      tmptable(1:acc_tsize) = acc_valtable
      tmptable(acc_tsize+1:2*acc_tsize) = 0.0_c_double
      acc_tsize = 2*acc_tsize
      deallocate(acc_valtable)
      acc_valtable => tmptable
      nullify(tmptable)
    end if
    acc_valtable(acc_icycle) = val
    write(message, '(A,ES12.3)') 'Entry '//I2S(acc_icycle) // &
        ' for convergence criterion: ', val
    call Info('CyclicConvergence', message, level=6)
  end subroutine cca_push

  function cca_same_cycle() result(r) bind(c, name="cca_same_cycle")
    integer(c_int) :: r
    ! Same cycle, nothing to do!
    r = merge(1, 0, acc_ncycle == acc_prevcycle)
  end function cca_same_cycle

  subroutine cca_stats() bind(c, name="cca_stats")
    integer :: n, me
    real(c_double) :: meanval, meandev, errval, errdev, errrange, &
        minv, maxv, rangev
    ! Compute mean errors over the cycle
    n = acc_icycle
    ! We assume here that the value that we follow has been communicated
    ! in parallel as intended so there is no need to communicate things
    ! again.
    meanval = sum(acc_valtable(1:n)) / n
    meandev = sqrt(sum((acc_valtable(1:n)-meanval)**2) / n)
    minv = minval(acc_valtable(1:n))
    maxv = maxval(acc_valtable(1:n))
    rangev = maxv - minv
    errval = 2*abs(meanval-acc_prevval)/(abs(meanval)+abs(acc_prevval))
    errdev = 2*abs(meandev-acc_prevdev)/(abs(meandev)+abs(acc_prevdev))
    errrange = 2*abs(rangev-acc_prevrangev)/(abs(meanval)+abs(acc_prevval))
    acc_meanval = meanval
    acc_meandev = meandev
    acc_errval = errval
    acc_errdev = errdev
    acc_errrange = errrange
    ! Take the smallest entry participating here and let it write the file
    ! This also ensures that we sync different processes.
    me = ParEnv % Mype
    me = ParallelReduction(me, 1)
    if (me == ParEnv % MyPe) call WriteConvergenceData(meanval, meandev, &
        errval, errdev, errrange)
    acc_prevval = meanval
    acc_prevdev = meandev
    acc_prevminv = minv
    acc_prevmaxv = maxv
    acc_prevrangev = rangev
  end subroutine cca_stats

  function cca_min_pending() result(r) bind(c, name="cca_min_pending")
    integer(c_int) :: r
    logical :: found
    integer :: mincycle
    mincycle = ListGetInteger(acc_params, &
        'Cyclic System Min Iterations', found)
    if (acc_ncycle < mincycle) then
      call Info('CyclicConvergence', &
          'Number of cycles smaller than required minimum', level=8)
      r = 1
    else
      r = 0
    end if
  end function cca_min_pending

  function cca_is_converged() result(r) bind(c, name="cca_is_converged")
    integer(c_int) :: r
    r = merge(1, 0, acc_converged)
  end function cca_is_converged

  subroutine cca_evaluate() bind(c, name="cca_evaluate")
    logical :: found
    real(c_double) :: tol
    integer :: nconv, ncnt
    type(variable_t), pointer :: pvar
    ! All the given measure must converge.
    nconv = 0
    ncnt = 0
    tol = ListGetCReal(acc_params, 'Mean Value Tolerance', found)
    if (found) then
      ncnt = ncnt + 1
      if (acc_errval < tol) nconv = nconv + 1
    end if
    tol = ListGetCReal(acc_params, 'Mean Deviation Tolerance', found)
    if (found) then
      ncnt = ncnt + 1
      if (acc_errdev < tol) nconv = nconv + 1
    end if
    tol = ListGetCReal(acc_params, 'Range Tolerance', found)
    if (found) then
      ncnt = ncnt + 1
      if (acc_errrange < tol) nconv = nconv + 1
    end if
    if (ncnt == 0) then
      call Fatal('CyclicConvergence', &
          'No tolerances given to evaluate convergence!')
    end if
    call Info('CyclicConvergence', I2S(nconv)//' measures converged ' // &
        'out of '//I2S(ncnt), level=10)
    acc_converged = (nconv == ncnt)
    if (acc_converged) then
      call Info('CyclicConvergence', 'Cyclic convergence reached after ' // &
          I2S(acc_ncycle)//' cycles!', level=5)
      pvar => VariableGet(acc_mesh % variables, 'produce')
      pvar % values(1) = 1.0_c_double
      acc_convcycle = acc_ncycle
    else
      call Info('CyclicConvergence', &
          'Cyclic convergence not reached after '//I2S(acc_ncycle) // &
          ' cycles!', level=10)
    end if
  end subroutine cca_evaluate

  subroutine cca_produce() bind(c, name="cca_produce")
    logical :: found
    integer :: prodcycles
    type(variable_t), pointer :: pvar
    prodcycles = ListGetInteger(acc_params, &
        'Number of Production Cycles', found)
    if (.not. found) prodcycles = 1
    if (acc_ncycle - acc_convcycle == prodcycles) then
      call Info('CyclicConvergence', 'Finish flag was activated after ' // &
          I2S(acc_ncycle)//' cycles!', level=5)
      pvar => VariableGet(acc_mesh % variables, 'finish')
      pvar % values(1) = 1.0_c_double
    end if
  end subroutine cca_produce

  subroutine cca_reset() bind(c, name="cca_reset")
    ! Restart the counter
    acc_icycle = 0
    acc_prevcycle = acc_ncycle
  end subroutine cca_reset

  ! -- internal procedure (exact transcription) --------------------------
  subroutine WriteConvergenceData(meanval, meandev, errval, errdev, &
      errrange)
    real(c_double) :: meanval, meandev, errval, errdev, errrange
    real(c_double) :: convdata(7)
    integer :: convunit
    character(len=MAX_NAME_LEN) :: convfile
    type(variable_t), pointer :: pvar
    logical :: found
    convfile = ListGetString(acc_params, 'Filename', found)
    if (.not. found) convfile = 'convergence.dat'
    if (acc_convvisited) then
      open(newunit=convunit, file=convfile, status='old', &
          position='append')
    else
      open(newunit=convunit, file=convfile)
      ! Write info in the exact same width as the results
      write(convunit, '(A1,A14,6A15)') '!', 'time', 'nCycle', 'meanVal', &
          'meanDev', 'errVal', 'errDev', 'errRange'
      acc_convvisited = .true.
    end if
    pvar => VariableGet(acc_solver % mesh % variables, 'time')
    convdata(1) = pvar % values(1)
    convdata(2) = acc_percycle % values(1)
    convdata(3) = meanval
    convdata(4) = meandev
    convdata(5) = errval
    convdata(6) = errdev
    convdata(7) = errrange
    write(convunit, '(7ES15.6)') convdata
    close(convunit)
  end subroutine WriteConvergenceData

  ! -- CyclicConvergence_Init body --------------------------------------
  subroutine cci_init_begin(h_solver) bind(c, name="cci_init_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine cci_init_begin

  subroutine cci_init_body() bind(c, name="cci_init_body")
    type(valuelist_t), pointer :: solverpars
    solverpars => acc_solver % values
    call ListAddNewString(solverpars, 'Equation', 'CyclicConvergence')
    ! Introduce variable so that we have a place for pseudonorm
    call ListAddNewString(solverpars, 'Variable', &
        '-nooutput -global convcrit_var')
  end subroutine cci_init_body

end module cyclicconvergence_acc
