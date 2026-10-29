! T022-A BeamSolver3D port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/BeamSolver3D.F90, compiled against the
! built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! hidden lengths, assumed-shape array components passed to DisplaceMesh,
! the BeamStiffnessMatrix optional-arg calls — stays inside Fortran, so
! the Rust side is pure C calls. The wrapper module mirrors
! BeamSolver3D's own locals (current element, entry-level flags) with
! identical lifetimes. All floating-point work (the Timoshenko beam
! assembly in SolidMechanicsUtils, MATMUL reductions, the nonlinear
! iteration control) stays in Fortran (same compiler), so bit-identity
! holds by construction wherever the call sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module beamsolver3d_acc
  use iso_c_binding
  use defutils
  use solidmechanicsutils
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(mesh_t), pointer, save :: acc_mesh => null()

contains

  ! -- entry state ----------------------------------------------------
  subroutine b3d_begin(h_model, h_solver) bind(c, name="b3d_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine b3d_begin

  subroutine b3d_default_start() bind(c, name="b3d_default_start")
    call DefaultStart()
  end subroutine b3d_default_start

  subroutine b3d_check_axes() bind(c, name="b3d_check_axes")
    if (.not. ListCheckPresentAnyMaterial(acc_model, 'Principal Direction 2') .and. &
        .not. ListCheckPresentAnyMaterial(acc_model, 'Director')) then
      call Warn('TimoshenkoSolver', &
          'Principal axes unspecified, assuming a circular cross section')
    end if
  end subroutine b3d_check_axes

  function b3d_maxiter() result(n) bind(c, name="b3d_maxiter")
    integer(c_int) :: n
    type(valuelist_t), pointer :: params
    logical :: found
    params => GetSolverParams()
    n = ListGetInteger(params, 'Nonlinear System Max Iterations', found, minv=1)
    if (.not. found) n = 1
  end function b3d_maxiter

  subroutine b3d_modes(transient, harmonic, mass) bind(c, name="b3d_modes")
    integer(c_int), value :: transient
    integer(c_int) :: harmonic, mass
    type(valuelist_t), pointer :: params
    logical :: found, harm, tsim
    params => GetSolverParams()
    harm = EigenOrHarmonicAnalysis() .or. ListGetLogical(params, 'Harmonic Mode', found)
    tsim = transient /= 0
    harmonic = merge(1, 0, harm)
    mass = merge(1, 0, tsim .or. harm)
  end subroutine b3d_modes

  function b3d_displace_active() result(r) bind(c, name="b3d_displace_active")
    integer(c_int) :: r
    type(valuelist_t), pointer :: params
    logical :: found
    params => GetSolverParams()
    r = merge(1, 0, GetLogical(params, 'Displace Mesh', found))
  end function b3d_displace_active

  subroutine b3d_mesh_return() bind(c, name="b3d_mesh_return")
    acc_mesh => GetMesh()
    call Info('TimoshenkoSolver', 'Returning the mesh to its reference position', level=4)
    call DisplaceMesh(acc_mesh, acc_solver % Variable % Values, -1, &
        acc_solver % Variable % Perm, 6, .false., 3)
  end subroutine b3d_mesh_return

  subroutine b3d_displace_forward() bind(c, name="b3d_displace_forward")
    call Info('TimoshenkoSolver', 'Displacing the mesh with computed displacement field', level=4)
    call DisplaceMesh(acc_mesh, acc_solver % Variable % Values, 1, &
        acc_solver % Variable % Perm, 6, .false., 3)
  end subroutine b3d_displace_forward

  ! -- nonlinear iteration body ---------------------------------------
  subroutine b3d_default_initialize() bind(c, name="b3d_default_initialize")
    call DefaultInitialize()
  end subroutine b3d_default_initialize

  function b3d_nof_active() result(n) bind(c, name="b3d_nof_active")
    integer(c_int) :: n
    n = GetNOFActive()
  end function b3d_nof_active

  ! -- per-element assembly (exact transcription of the loop body) -----
  subroutine b3d_elem(t, transient, mass, harmonic) bind(c, name="b3d_elem")
    integer(c_int), value :: t, transient, mass, harmonic
    integer :: n, nd, nb
    logical :: tsim, massl, harml
    acc_elem => GetActiveElement(t)
    if (.not. (GetElementFamily(acc_elem) == 2)) return
    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
    nb = GetElementNOFBDOFs()
    tsim = transient /= 0
    massl = mass /= 0
    harml = harmonic /= 0
    call BeamStiffnessMatrix(acc_elem, n, nd+nb, nb, tsim, &
        MassAssembly=massl, HarmonicAssembly=harml)
  end subroutine b3d_elem

  ! -- shared assembly tail --------------------------------------------
  subroutine b3d_finish_bulk() bind(c, name="b3d_finish_bulk")
    call DefaultFinishBulkAssembly()
  end subroutine b3d_finish_bulk

  subroutine b3d_finish_boundary() bind(c, name="b3d_finish_boundary")
    call DefaultFinishBoundaryAssembly()
  end subroutine b3d_finish_boundary

  subroutine b3d_finish_assembly() bind(c, name="b3d_finish_assembly")
    call DefaultFinishAssembly()
  end subroutine b3d_finish_assembly

  subroutine b3d_dirichlet() bind(c, name="b3d_dirichlet")
    call DefaultDirichletBCs()
  end subroutine b3d_dirichlet

  function b3d_solve() result(norm) bind(c, name="b3d_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function b3d_solve

  function b3d_converged() result(r) bind(c, name="b3d_converged")
    integer(c_int) :: r
    r = merge(1, 0, DefaultConverged())
  end function b3d_converged

  subroutine b3d_finish() bind(c, name="b3d_finish")
    call DefaultFinish()
  end subroutine b3d_finish

  ! -- TimoshenkoSolver_Init0 body (exact transcription) ----------------
  subroutine b3d_init0(h_solver, transient) bind(c, name="b3d_init0")
    type(c_ptr), value :: h_solver
    integer(c_int), value :: transient
    type(valuelist_t), pointer :: solverpars, simulation
    logical :: found, meshactive
    call c_f_pointer(h_solver, acc_solver)
    solverpars => GetSolverParams()
    call ListAddInteger(solverpars, 'Variable DOFs', 6)
    call ListAddNewString(solverpars, 'Variable', 'Deflection[U:3 Theta:3]')
    call ListAddNewString(solverpars, 'Element', 'p:1 b:1')
    call ListAddNewLogical(solverpars, 'Bubbles in Global System', .false.)
    call ListAddNewLogical(solverpars, 'Use Global Mass Matrix', .true.)
    if (transient /= 0) then
      call ListAddInteger(solverpars, 'Time derivative order', 2)
      call ListAddNewString(solverpars, 'Timestepping Method', 'Bossak')
    end if
    meshactive = GetLogical(solverpars, 'Displace Mesh', found)
    if (meshactive) then
      simulation => GetSimulation()
      call ListAddLogical(simulation, 'Initialize Dirichlet Conditions', .false.)
    end if
    call ListAddNewLogical(solverpars, 'Beam Solver', .true.)
  end subroutine b3d_init0

end module beamsolver3d_acc
