! Port of `fem/src/modules/VelocityPrecond.F90` accessor calls for lane
! T021-B, in the T020 Poisson pattern (ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (none here: the
! original `VelocityPrecond_Init0` body calls `GetSolverParams()` with no
! arguments and then issues twelve unconditional `ListAdd*` statements
! with literal arguments). All Fortran-object traffic — the
! `TYPE(ValueList_t)` pointer, OPTIONAL dummies, hidden CHARACTER
! lengths — stays inside Fortran, so the Rust side is pure C calls. The
! module holds the `SolverParams` pointer exactly as the original holds
! its local (set once per entry call, before the twelve statements).
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module velocityprecond_acc
  use iso_c_binding
  use defutils
  implicit none

  type(valuelist_t), pointer, save :: acc_params => null()

contains

  subroutine vpa_begin() bind(c, name="vpa_begin")
    acc_params => GetSolverParams()
  end subroutine vpa_begin

  subroutine vpa_add_exec_solver() bind(c, name="vpa_add_exec_solver")
    call ListAddString(acc_params, 'Exec Solver', 'before simulation')
  end subroutine vpa_add_exec_solver

  subroutine vpa_add_variable_output() bind(c, name="vpa_add_variable_output")
    call ListAddNewLogical(acc_params, 'Variable Output', .FALSE.)
  end subroutine vpa_add_variable_output

  subroutine vpa_add_bubbles() bind(c, name="vpa_add_bubbles")
    call ListAddNewLogical(acc_params, 'Bubbles in Global System', .FALSE.)
  end subroutine vpa_add_bubbles

  subroutine vpa_add_skip_change() bind(c, name="vpa_add_skip_change")
    call ListAddLogical(acc_params, 'Skip Compute Nonlinear Change', .TRUE.)
  end subroutine vpa_add_skip_change

  subroutine vpa_add_back_rotate() bind(c, name="vpa_add_back_rotate")
    call ListAddLogical(acc_params, 'Back Rotate N-T Solution', .FALSE.)
  end subroutine vpa_add_back_rotate

  subroutine vpa_add_variable() bind(c, name="vpa_add_variable")
    call ListAddNewString(acc_params, 'Variable', 'V')
  end subroutine vpa_add_variable

  subroutine vpa_add_linsys_solver() bind(c, name="vpa_add_linsys_solver")
    call ListAddNewString(acc_params, 'Linear System Solver', 'Iterative')
  end subroutine vpa_add_linsys_solver

  subroutine vpa_add_linsys_method() bind(c, name="vpa_add_linsys_method")
    call ListAddNewString(acc_params, 'Linear System Iterative Method', 'BiCGStabL')
  end subroutine vpa_add_linsys_method

  subroutine vpa_add_linsys_maxiter() bind(c, name="vpa_add_linsys_maxiter")
    call ListAddNewInteger(acc_params, 'Linear System Max Iterations', 1000)
  end subroutine vpa_add_linsys_maxiter

  subroutine vpa_add_linsys_precond() bind(c, name="vpa_add_linsys_precond")
    call ListAddNewString(acc_params, 'Linear System Preconditioning', 'ILU0')
  end subroutine vpa_add_linsys_precond

  subroutine vpa_add_linsys_tol() bind(c, name="vpa_add_linsys_tol")
    call ListAddNewConstReal(acc_params, 'Linear System Convergence Tolerance', 1.0d-4)
  end subroutine vpa_add_linsys_tol

  subroutine vpa_add_linsys_abort() bind(c, name="vpa_add_linsys_abort")
    call ListAddNewLogical(acc_params, 'Linear System Abort Not Converged', .FALSE.)
  end subroutine vpa_add_linsys_abort

end module velocityprecond_acc
