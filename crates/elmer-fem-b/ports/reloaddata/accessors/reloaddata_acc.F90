! T021-B ReloadData port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/ReloadData.F90, compiled against the built .mod
! files of the target configure (A2 shipped mechanism, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars). All Fortran-object traffic — TYPE(...)
! handles, OPTIONAL/POINTER dummies, hidden lengths, allocatable strings —
! stays inside Fortran, so the Rust side is pure C calls.
!
! The wrapper module mirrors ReloadData's own state: the entry's SAVE
! scalars live in Rust; the per-invocation RestartFile local is staged in
! wrapper SAVE (filled by rda_has_restart, consumed by rda_info_reading /
! rda_load before the next fill, as in the original where it is re-read
! every invocation). Order-sensitive formatted I/O (WRITE(Message,*)) and
! LoadRestartFile stay in Fortran so bit-identity holds by construction
! wherever the call sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module reloaddata_acc
  use iso_c_binding
  use DefUtils
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  character(len=MAX_NAME_LEN), save :: acc_restart = ''

contains

  ! -- entry state ----------------------------------------------------
  subroutine rda_begin(h_solver) bind(c, name="rda_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine rda_begin

  ! -- Visit == 0 init reads ------------------------------------------
  subroutine rda_get_int(which, val, found) bind(c, name="rda_get_int")
    integer(c_int), value :: which
    integer(c_int) :: val, found
    logical :: f
    integer :: v
    select case (which)
    case (1)
      v = GetInteger(acc_solver % values, 'Reload Range Minimum', f)
    case (2)
      v = GetInteger(acc_solver % values, 'Reload Range Maximum', f)
    case (3)
      v = GetInteger(acc_solver % values, 'Reload Starting Position', f)
    case default
      v = GetInteger(acc_solver % values, 'Reload Reading Intervals', f)
    end select
    val = v
    found = merge(1, 0, f)
  end subroutine rda_get_int

  subroutine rda_get_cont(val, found) bind(c, name="rda_get_cont")
    integer(c_int) :: val, found
    logical :: f, l
    l = GetLogical(acc_solver % values, 'Continuous Reading', f)
    val = merge(1, 0, l)
    found = merge(1, 0, f)
  end subroutine rda_get_cont

  subroutine rda_info_cont_disabled() bind(c, name="rda_info_cont_disabled")
    call Info('ReloadSolution', 'Continuous reading disabled', Level=8)
  end subroutine rda_info_cont_disabled

  ! -- per-invocation restart file ------------------------------------
  subroutine rda_has_restart(found) bind(c, name="rda_has_restart")
    integer(c_int) :: found
    logical :: f
    acc_restart = GetString(acc_solver % values, 'Reload Solution File', f)
    found = merge(1, 0, f)
  end subroutine rda_has_restart

  subroutine rda_info_reading() bind(c, name="rda_info_reading")
    write(message, *) 'Reading old solution from file ', trim(acc_restart)
    call Info('ReloadSolution', message, Level=4)
  end subroutine rda_info_reading

  ! -- load step (mesh lookup, OutputName build, Info pair, load) -----
  subroutine rda_load(step, cont) bind(c, name="rda_load")
    integer(c_int), value :: step, cont
    type(mesh_t), pointer :: m
    integer :: k
    character(len=MAX_NAME_LEN) :: outputname
    k = step
    m => acc_solver % mesh
    if (len_trim(m % name) > 0) then
      outputname = trim(m % name) // '/' // trim(acc_restart)
      write(message, *) 'Opening file: ', trim(outputname)
    else
      outputname = trim(acc_restart)
      write(message, *) 'Opening file: ', trim(outputname)
    end if
    call Info('ReloadSolution', message, Level=4)
    if (parenv % pes > 1) outputname = trim(outputname) // '.' // i2s(parenv % mype)
    write(message, *) 'Loading Timestep: ', k
    call Info('ReloadSolution', message, Level=4)
    if (cont /= 0) then
      call LoadRestartFile(outputname, k, m, .true.)
    else
      call LoadRestartFile(outputname, k, m)
    end if
  end subroutine rda_load

end module reloaddata_acc
