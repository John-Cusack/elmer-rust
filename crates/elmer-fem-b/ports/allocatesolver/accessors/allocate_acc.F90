! T021-A AllocateSolver port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/AllocateSolver.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments. All Fortran-object
! traffic — TYPE(...) handles, OPTIONAL/POINTER dummies, assumed-shape
! descriptors, hidden lengths — stays inside Fortran, so the Rust side is
! pure C calls. Clean-room code: statement-for-statement wrappers around
! Elmer library calls, no upstream text.
module allocate_acc
  use iso_c_binding
  use defutils
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()

contains

  ! -- entry state ----------------------------------------------------
  subroutine ala_begin(h_solver) bind(c, name="ala_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine ala_begin

  ! -- AllocateSolver_init body ---------------------------------------
  subroutine ala_init_exec_never() bind(c, name="ala_init_exec_never")
    call ListAddNewString(acc_solver % Values, 'Exec Solver', 'never')
  end subroutine ala_init_exec_never

  ! -- AllocateSolver body --------------------------------------------
  subroutine ala_info() bind(c, name="ala_info")
    call Info('AllocateSolver', &
        'This solver does nothing, it is just for allocating stuff')
  end subroutine ala_info

end module allocate_acc
