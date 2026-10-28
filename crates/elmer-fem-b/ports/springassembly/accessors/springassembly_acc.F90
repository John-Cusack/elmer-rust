! T021-B SpringAssembly port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/SpringAssembly.F90, compiled against the
! built .mod files of the target configure (ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, function-result arrays —
! stays inside Fortran, so the Rust side is pure C calls.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module springassembly_acc
  use iso_c_binding
  use types
  use defutils
  use eleminfo
  use messages
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(solver_t), pointer, save :: acc_solverptr => null()
  type(valuelist_t), pointer, save :: acc_vlist => null()
  integer, pointer, save :: acc_nodeidx(:) => null()
  ! Point element serving the standard assembly routine (mirrors the
  ! original's call-local TARGET Element with identical lifetime: fully
  ! re-pointed on every entry call).
  type(element_t), target, save :: acc_elem

contains

  ! -- entry state ----------------------------------------------------
  subroutine spa_begin(h_model, h_solver) bind(c, name="spa_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine spa_begin

  ! -- prologue: solver search, HarmonicAssembly, DOFs ------------------
  ! Mirrors the original's statement order exactly, including reuse of
  ! Found across GetString/GetLogical and the non-short-circuit .OR.
  subroutine spa_prologue(dofs, harmonic) bind(c, name="spa_prologue")
    integer(c_int) :: dofs, harmonic
    type(valuelist_t), pointer :: solverparams, valuelist
    character(len=MAX_NAME_LEN) :: uname
    logical :: found, hassemble
    integer :: i

    solverparams => GetSolverParams()

    uname = GetString(solverparams, 'Displacement Variable Name', found)
    if (.not. found) uname = 'displacement'
    found = .false.
    do i = 1, acc_model % numberofsolvers
      acc_solverptr => acc_model % solvers(i)
      if (uname == GetVarName(acc_solverptr % variable)) then
        found = .true.
        exit
      end if
    end do
    if (.not. found) then
      call Fatal('SpringAssembler', &
          'Solver associated with Displacement Variable Name > '// &
          trim(uname)//' < not found!')
    end if

    valuelist => GetSolverParams(acc_solverptr)
    hassemble = EigenOrHarmonicAnalysis(acc_solverptr) .or. &
        GetLogical(valuelist, 'Harmonic Mode', found) .or. &
        GetLogical(valuelist, 'Harmonic Analysis', found)

    dofs = acc_solverptr % variable % dofs
    harmonic = merge(1, 0, hassemble)
  end subroutine spa_prologue

  subroutine spa_elem_init() bind(c, name="spa_elem_init")
    acc_elem % type => GetElementType(101, .false.)
    acc_elem % ndofs = 1
  end subroutine spa_elem_init

  function spa_num_bcs() result(n) bind(c, name="spa_num_bcs")
    integer(c_int) :: n
    n = acc_model % numberofbcs
  end function spa_num_bcs

  ! -- per-BC context ---------------------------------------------------
  function spa_bc_begin(bc) result(r) bind(c, name="spa_bc_begin")
    integer(c_int), value :: bc
    integer(c_int) :: r
    acc_vlist => acc_model % bcs(bc) % values
    r = merge(1, 0, ListCheckPresent(acc_vlist, 'Target Nodes'))
  end function spa_bc_begin

  function spa_bc_numnodes() result(n) bind(c, name="spa_bc_numnodes")
    integer(c_int) :: n
    acc_nodeidx => ListGetIntegerArray(acc_vlist, 'Target Nodes')
    n = size(acc_nodeidx)
  end function spa_bc_numnodes

  ! Scalar-'Spring' probe: exact transcription of
  ! Work(1,1:n) = ListGetReal(... 'Spring' ...) plus the Warn pair.
  subroutine spa_spring_scalar(n, ld, out) bind(c, name="spa_spring_scalar")
    integer(c_int), value :: n, ld
    real(c_double) :: out(ld, n)
    logical :: found
    out(1, 1:n) = ListGetReal(acc_vlist, 'Spring', n, acc_nodeidx, found)
    if (found) then
      call Warn('SpringAssembler', 'Define a spring by components (Spring i = ...)')
      call Warn('SpringAssembler', 'Skipping a definition (Spring = ...)')
    end if
  end subroutine spa_spring_scalar

  ! Per-component row read: exact transcription of
  ! Work(i,1:n) = ListGetReal(... ComponentName('Spring',i) ...) /
  ! MassVals(i,1:n) = ListGetReal(... ComponentName('Mass',i) ...).
  subroutine spa_comp_row(which, comp, n, ld, out, found) &
      bind(c, name="spa_comp_row")
    integer(c_int), value :: which, comp, n, ld
    real(c_double) :: out(ld, n)
    integer(c_int) :: found
    logical :: fl
    if (which == 1) then
      out(comp, 1:n) = ListGetReal(acc_vlist, &
          ComponentName('Spring', comp), n, acc_nodeidx, fl)
    else
      out(comp, 1:n) = ListGetReal(acc_vlist, &
          ComponentName('Mass', comp), n, acc_nodeidx, fl)
    end if
    found = merge(1, 0, fl)
  end subroutine spa_comp_row

  ! -- per-node context ---------------------------------------------------
  subroutine spa_set_elem_node(i) bind(c, name="spa_set_elem_node")
    integer(c_int), value :: i
    acc_elem % nodeindexes => acc_nodeidx(i:i)
  end subroutine spa_set_elem_node

  function spa_node_perm(i) result(p) bind(c, name="spa_node_perm")
    integer(c_int), value :: i
    integer(c_int) :: p
    p = acc_solverptr % variable % perm(acc_nodeidx(i))
  end function spa_node_perm

  ! -- assembly tail ------------------------------------------------------
  subroutine spa_update_equations(dd, k, f) &
      bind(c, name="spa_update_equations")
    integer(c_int), value :: dd
    real(c_double) :: k(dd, dd), f(dd)
    call DefaultUpdateEquations(k(1:dd, 1:dd), f(1:dd), acc_elem, acc_solverptr)
  end subroutine spa_update_equations

  subroutine spa_2ndordertime(dd, m, d, k, f) &
      bind(c, name="spa_2ndordertime")
    integer(c_int), value :: dd
    real(c_double) :: m(dd, dd), d(dd, dd), k(dd, dd), f(dd)
    call Default2ndOrderTime(m(1:dd, 1:dd), d(1:dd, 1:dd), &
        k(1:dd, 1:dd), f(1:dd))
  end subroutine spa_2ndordertime

  subroutine spa_update_mass(dd, m) bind(c, name="spa_update_mass")
    integer(c_int), value :: dd
    real(c_double) :: m(dd, dd)
    call DefaultUpdateMass(m(1:dd, 1:dd), acc_elem, acc_solverptr)
  end subroutine spa_update_mass

  subroutine spa_fatal_alloc() bind(c, name="spa_fatal_alloc")
    call Fatal('SpringAssembler', 'Memory allocation error.')
  end subroutine spa_fatal_alloc

end module springassembly_acc
