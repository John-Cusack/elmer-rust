! SaveData pilot: statement-level BIND(C) wrappers for the Elmer calls in
! fem/src/modules/SaveData/SaveDependence.F90.
module sd_depend_acc
  use iso_c_binding
  use defutils
  use types
  use SaveUtils, only: SolverOutputDirectory
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(model_t), pointer, save :: acc_model => null()
  integer, save :: acc_iounit = -1
  integer, save :: acc_noparse = 0
  real(kind=dp), save :: acc_normsum = 0.0_dp
  real(kind=dp), save :: acc_x0 = 0.0_dp, acc_x1 = 1.0_dp

contains

  subroutine sdd_begin(h_model, h_solver) bind(c, name="sdd_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine sdd_begin

  function sdd_parallel_skip() result(r) bind(c, name="sdd_parallel_skip")
    integer(c_int) :: r
    r = 0
    if (parenv % pes > 1) then
      if (parenv % mype > 0) r = 1
    end if
  end function sdd_parallel_skip

  subroutine sdd_info_start() bind(c, name="sdd_info_start")
    call Info('SaveDependence', '--------------------------------')
    call Info('SaveDependence', 'Saving dependencies in a table')
  end subroutine sdd_info_start

  subroutine sdd_init() bind(c, name="sdd_init")
    integer :: normind
    logical :: gotit
    character(len=4096) :: solvername
    normind = ListGetInteger(acc_solver % values, 'Show Norm Index', gotit)
    if (normind > 0) then
      solvername = ListGetString(acc_solver % values, 'Equation', gotit)
      if (.not. ListCheckPresent(acc_solver % values, 'Variable')) then
        call ListAddString(acc_solver % values, 'Variable', &
            '-nooutput -global ' // trim(solvername) // '_var')
      end if
    end if
  end subroutine sdd_init

  subroutine sdd_setup(n, nopar, normind, takeder) bind(c, name="sdd_setup")
    integer(c_int) :: n, nopar, normind, takeder
    type(valuelist_t), pointer :: params
    character(len=4096) :: filename, parname
    character(:), allocatable :: outdir
    logical :: found, gotit
    integer :: j
    real(kind=dp) :: x0, x1
    params => GetSolverParams()
    filename = ListGetString(params, 'Filename', found)
    if (.not. found) filename = 'dep.dat'
    normind = ListGetInteger(params, 'Show Norm Index', gotit)
    call SolverOutputDirectory(acc_solver, filename, outdir)
    filename = trim(outdir) // '/' // trim(filename)
    if (GetLogical(params, 'Filename Numbering', gotit)) then
      filename = NextFreeFilename(filename)
    end if
    takeder = merge(1, 0, ListGetLogical(params, 'Take Derivative', gotit))
    n = ListGetInteger(params, 'Number of points', minv=2)
    x0 = ListGetCReal(params, 'Lower limit')
    x1 = ListGetCReal(params, 'Upper Limit')
    acc_x0 = x0
    acc_x1 = x1
    acc_noparse = 0
    acc_normsum = 0.0_dp
    nopar = 0
    do j = 1, 100
      write (parname, '(A,I0)') 'Expression ', j
      if (ListCheckPresent(params, trim(parname))) then
        nopar = j
      else
        exit
      end if
    end do
    acc_noparse = nopar
    if (nopar == 0) then
      call Warn('SaveDependence', 'No parameter given!')
      return
    end if
    open (newunit=acc_iounit, file=trim(filename))
  end subroutine sdd_setup

  subroutine sdd_write_row(i, n, normind) bind(c, name="sdd_write_row")
    integer(c_int), value :: i, n, normind
    type(valuelist_t), pointer :: params
    character(len=4096) :: parname
    real(kind=dp) :: f, dfdx, w, xx
    integer :: j
    logical :: takeder
    params => GetSolverParams()
    takeder = ListGetLogical(params, 'Take Derivative')
    w = (1.0_dp * (i - 1)) / (n - 1)
    xx = acc_x0 + w * (acc_x1 - acc_x0)
    write (acc_iounit, '(I6,ES15.6)', advance='NO') i, xx
    do j = 1, acc_noparse
      write (parname, '(A,I0)') 'Expression ', j
      if (takeder) then
        f = ListGetFun(params, trim(parname), xx, dfdx=dfdx)
        write (acc_iounit, '(2ES15.6)', advance='NO') f, dfdx
      else
        f = ListGetFun(params, trim(parname), xx)
        write (acc_iounit, '(ES15.6)', advance='NO') f
      end if
      if (normind == j) acc_normsum = acc_normsum + f * f
    end do
    write (acc_iounit, '(A)') ' '
  end subroutine sdd_write_row

  subroutine sdd_finish(n, normind) bind(c, name="sdd_finish")
    integer(c_int), value :: n, normind
    real(kind=dp) :: norm
    close (acc_iounit)
    acc_iounit = -1
    if (normind > 0) then
      norm = sqrt(acc_normsum / n)
      acc_solver % variable % values = norm
      acc_solver % variable % norm = norm
    end if
    acc_normsum = 0.0_dp
  end subroutine sdd_finish

end module sd_depend_acc
