! T021-B PressurePrecond port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/PressurePrecond.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (no arguments at
! all here: the entry dummies are unused by the original, and every
! ListAdd* literal is transcribed verbatim into its own wrapper, so no
! CHARACTER dummy — with its hidden length — ever crosses the boundary).
! All Fortran-object traffic — the SolverParams POINTER from
! GetSolverParams(), the assumed-length literals — stays inside Fortran,
! so the Rust side is pure C calls. The wrapper module mirrors the
! original's SolverParams local with identical lifetime (set by
! ppa_begin on every entry call, consumed by the twelve init wrappers in
! the original's statement order, exactly as in the original).
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module pressureprecond_acc
  use iso_c_binding
  use types
  use defutils
  implicit none

  type(valuelist_t), pointer, save :: acc_params => null()

contains

  ! -- entry state ----------------------------------------------------
  subroutine ppa_begin() bind(c, name="ppa_begin")
    acc_params => GetSolverParams()
  end subroutine ppa_begin

  ! -- PressurePrecond_Init0 body, one wrapper per statement -----------
  subroutine ppa_init_exec_solver() bind(c, name="ppa_init_exec_solver")
    call ListAddString(acc_params, 'Exec Solver', 'before simulation')
  end subroutine ppa_init_exec_solver

  subroutine ppa_init_variable_output() bind(c, name="ppa_init_variable_output")
    call ListAddLogical(acc_params, 'Variable Output', .FALSE.)
  end subroutine ppa_init_variable_output

  subroutine ppa_init_bubbles() bind(c, name="ppa_init_bubbles")
    call ListAddNewLogical(acc_params, 'Bubbles in Global System', .FALSE.)
  end subroutine ppa_init_bubbles

  subroutine ppa_init_skip_change() bind(c, name="ppa_init_skip_change")
    call ListAddLogical(acc_params, 'Skip Compute Nonlinear Change', .TRUE.)
  end subroutine ppa_init_skip_change

  subroutine ppa_init_back_rotate() bind(c, name="ppa_init_back_rotate")
    call ListAddLogical(acc_params, 'Back Rotate N-T Solution', .FALSE.)
  end subroutine ppa_init_back_rotate

  subroutine ppa_init_variable() bind(c, name="ppa_init_variable")
    call ListAddNewString(acc_params, 'Variable', 'P')
  end subroutine ppa_init_variable

  subroutine ppa_init_lin_solver() bind(c, name="ppa_init_lin_solver")
    call ListAddNewString(acc_params, 'Linear System Solver', 'Iterative')
  end subroutine ppa_init_lin_solver

  subroutine ppa_init_lin_method() bind(c, name="ppa_init_lin_method")
    call ListAddNewString(acc_params, 'Linear System Iterative Method', 'BiCGStab2')
  end subroutine ppa_init_lin_method

  subroutine ppa_init_lin_maxiter() bind(c, name="ppa_init_lin_maxiter")
    call ListAddNewInteger(acc_params, 'Linear System Max Iterations', 1000)
  end subroutine ppa_init_lin_maxiter

  subroutine ppa_init_lin_precond() bind(c, name="ppa_init_lin_precond")
    call ListAddNewString(acc_params, 'Linear System Preconditioning', 'Diagonal')
  end subroutine ppa_init_lin_precond

  subroutine ppa_init_lin_tol() bind(c, name="ppa_init_lin_tol")
    call ListAddNewConstReal(acc_params, 'Linear System Convergence Tolerance', 1.0d-6)
  end subroutine ppa_init_lin_tol

  subroutine ppa_init_lin_abort() bind(c, name="ppa_init_lin_abort")
    call ListAddNewLogical(acc_params, 'Linear System Abort Not Converged', .FALSE.)
  end subroutine ppa_init_lin_abort

end module pressureprecond_acc
