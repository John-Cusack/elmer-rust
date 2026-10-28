! T021-B ReloadInput port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/ReloadInput.F90, compiled against the built .mod
! files of the target configure (A2 shipped mechanism, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape results). All Fortran-object
! traffic — TYPE(...) handles, assumed-shape descriptors, hidden lengths —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module owns the entry's locals (ModelName, MeshDir, MeshName, IOUnit)
! with static lifetimes; they are overwritten on every path that observes
! them (exactly as the original's fresh locals are), so persistence across
! invocations is unobservable.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module reloadinput_acc
  use iso_c_binding
  use types
  use modeldescription
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  character(len=MAX_NAME_LEN), save :: acc_name = ' '
  character(len=MAX_NAME_LEN), save :: acc_meshdir = ' '
  character(len=MAX_NAME_LEN), save :: acc_meshname = ' '
  integer, save :: acc_iounit = -1

contains

  ! -- entry state ------------------------------------------------------
  subroutine rli_begin(h_model) bind(c, name="rli_begin")
    type(c_ptr), value :: h_model
    call c_f_pointer(h_model, acc_model)
  end subroutine rli_begin

  ! -- OPEN(NEWUNIT=IOUnit, ..., STATUS='OLD', ERR=10); 0 ok / 1 error --
  function rli_open1() result(rc) bind(c, name="rli_open1")
    integer(c_int) :: rc
    open(newunit=acc_iounit, file='ELMERSOLVER_REREADINFO', status='OLD', err=10)
    rc = 0
    return
10  rc = 1
  end function rli_open1

  ! -- READ(IOUnit,'(a)') ModelName; CLOSE(IOUnit) (main path) ----------
  subroutine rli_read1() bind(c, name="rli_read1")
    read(acc_iounit,'(a)') acc_name
    close(acc_iounit)
  end subroutine rli_read1

  ! -- OPEN(InFileUnit, File=ModelName); LoadInputFile(...); CLOSE ------
  subroutine rli_reload() bind(c, name="rli_reload")
    open(InFileUnit, file=acc_name)
    call LoadInputFile(acc_model, InFileUnit, acc_name, &
        acc_meshdir, acc_meshname, .FALSE., .FALSE.)
    close(InFileUnit)
  end subroutine rli_reload

  ! -- second OPEN (label-10 arm), ERR=20; 0 ok / 1 error ----------------
  function rli_open2() result(rc) bind(c, name="rli_open2")
    integer(c_int) :: rc
    open(newunit=acc_iounit, file='ELMERSOLVER_REREADINFO', status='OLD', err=20)
    rc = 0
    return
20  rc = 1
  end function rli_open2

  ! -- READ+CLOSE on the label-10 arm (falls through to label 20) --------
  subroutine rli_read2() bind(c, name="rli_read2")
    read(acc_iounit,'(a)') acc_name
    close(acc_iounit)
  end subroutine rli_read2

end module reloadinput_acc
