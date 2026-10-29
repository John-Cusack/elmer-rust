! T022-A CraigBamptonSolver port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/CraigBamptonSolver.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, list-directed file I/O —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module mirrors CraigBamptonSolver's own host-associated locals
! (NoEigenModes, NoConstraintModes, NoComponentModes, Dofs, FileName,
! SaveValues, MatrixNo, Ahat/Ax extents, the CONTAINS-internal file units)
! with identical lifetimes. The order-sensitive reduction (`SUM(x*Ax)`,
! the Frobenius norm) and every formatted/file WRITE stay in Fortran
! (same compiler, same evaluation order as the original internal
! procedures), so bit-identity holds by construction wherever the call
! sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module craigbamptonsolver_acc
  use iso_c_binding
  use defutils
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(solver_t), pointer, save :: acc_esolver => null()
  type(variable_t), pointer, save :: acc_evar => null()
  type(matrix_t), pointer, save :: acc_a => null()
  type(mesh_t), pointer, save :: acc_mesh => null()
  real(kind=dp), pointer, save :: acc_savevalues(:) => null()
  ! Entry-level state (mirrors the solver body's locals).
  character(len=MAX_NAME_LEN), save :: acc_elasticname = ''
  character(len=MAX_NAME_LEN), save :: acc_filename = ''
  integer, save :: acc_noeigen = 0
  integer, save :: acc_noconstraint = 0
  integer, save :: acc_nocomp = 0
  integer, save :: acc_dofs = 0
  integer, save :: acc_matrixno = 0

contains

  ! -- CraigBamptonSolver_Init body -------------------------------------
  subroutine cba_init_begin(h_solver) bind(c, name="cba_init_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine cba_init_begin

  subroutine cba_init_body() bind(c, name="cba_init_body")
    logical :: found
    character(len=MAX_NAME_LEN) :: solvername
    if (listgetlogical(acc_solver % values, 'Calculate Matrix Norm', found)) then
      solvername = listgetstring(acc_solver % values, 'Equation', found)
      if (.not. listcheckpresent(acc_solver % values, 'Variable')) then
        call listaddstring(acc_solver % values, 'Variable', &
            '-nooutput -global '//trim(solvername)//'_var')
      end if
    end if
  end subroutine cba_init_body

  ! -- entry state ------------------------------------------------------
  subroutine cba_begin(h_model, h_solver) bind(c, name="cba_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine cba_begin

  subroutine cba_info_start() bind(c, name="cba_info_start")
    call info('CraigBamptonSolver', 'Performing Craig-Bampton model reduction')
  end subroutine cba_info_start

  ! Resolve the elasticity solver, swap in its bulk values, and cache the
  ! reduction sizing. Returns 0 on success, else the failing check:
  ! 1 = no 'Elasticity Solver Name', 2 = no solver with that Equation
  ! name, 3 = matrix not associated, 4 = bulk values not associated,
  ! 5 = number of eigenmodes is zero. The caller maps each code to the
  ! matching Fatal wrapper below (same messages as the original).
  function cba_resolve() result(code) bind(c, name="cba_resolve")
    integer(c_int) :: code
    type(valuelist_t), pointer :: params
    character(len=MAX_NAME_LEN) :: str
    integer :: i, solverid
    logical :: found
    params => getsolverparams()
    str = listgetstring(params, 'Elasticity Solver Name', found)
    if (.not. found) then
      code = 1
      return
    end if
    acc_elasticname = str
    solverid = 0
    do i = 1, acc_model % numberofsolvers
      if (listgetstring(acc_model % solvers(i) % values, 'Equation', found) == str) then
        solverid = i
        exit
      end if
    end do
    if (solverid == 0) then
      code = 2
      return
    end if
    acc_esolver => acc_model % solvers(solverid)
    acc_evar => acc_esolver % variable
    acc_a => acc_esolver % matrix
    if (.not. associated(acc_a)) then
      code = 3
      return
    end if
    acc_savevalues => acc_a % values
    if (.not. associated(acc_a % bulkvalues)) then
      code = 4
      return
    end if
    acc_a % values => acc_a % bulkvalues
    acc_dofs = acc_evar % dofs
    acc_mesh => acc_esolver % mesh
    acc_noeigen = acc_esolver % nofeigenvalues
    if (acc_noeigen == 0) then
      code = 5
      return
    end if
    call info('CraigBamptonSolver', 'Number of eigen modes: '//i2s(acc_noeigen), level=7)
    acc_noconstraint = acc_evar % numberofconstraintmodes
    call info('CraigBamptonSolver', 'Number of constraint modes: '//i2s(acc_noconstraint), level=7)
    acc_nocomp = acc_noeigen + acc_noconstraint
    code = 0
  end function cba_resolve

  subroutine cba_fatal_noname() bind(c, name="cba_fatal_noname")
    call fatal('CraigBamptonSolver', '> Elasticity Solver Name < not given!')
  end subroutine cba_fatal_noname

  subroutine cba_fatal_nosolver() bind(c, name="cba_fatal_nosolver")
    call fatal('CraigBamptonSolver', 'No solver with Equation name: '//trim(acc_elasticname))
  end subroutine cba_fatal_nosolver

  subroutine cba_fatal_nomatrix() bind(c, name="cba_fatal_nomatrix")
    call fatal('CraigBamptonSolver', 'Matrix not associated!')
  end subroutine cba_fatal_nomatrix

  subroutine cba_fatal_nobulk() bind(c, name="cba_fatal_nobulk")
    call fatal('CraigBamptonSolver', 'Bulk values of matrix not associated!')
  end subroutine cba_fatal_nobulk

  subroutine cba_fatal_noeigen() bind(c, name="cba_fatal_noeigen")
    call fatal('CraigBamptonSolver', 'Number of eigenmodes is zero!')
  end subroutine cba_fatal_noeigen

  subroutine cba_fatal_alloc() bind(c, name="cba_fatal_alloc")
    call fatal('CraigBamptonSolver', 'Memory allocation error.')
  end subroutine cba_fatal_alloc

  ! -- SaveReductionDofs internal (exact transcription) ------------------
  subroutine cba_save_dofs() bind(c, name="cba_save_dofs")
    integer :: iounit, i, j
    acc_filename = 'ReductionHeader.dat'
    call info('CraigBamptonSolver', &
        'Saving information on the reduction basis to: '//trim(acc_filename), level=7)
    open (newunit=iounit, file=acc_filename)
    write(iounit, '(I0)') acc_noeigen
    write(iounit, '(I0)') acc_noconstraint
    close(iounit)
    acc_filename = 'ReductionIndeces.dat'
    call info('CraigBamptonSolver', &
        'Saving Constraint Modes Indeces to: '//trim(acc_filename), level=7)
    open (newunit=iounit, file=acc_filename)
    do i = 1, acc_a % numberofrows
      j = acc_evar % constraintmodesindeces(i)
      if (j == 0) cycle
      write(iounit, '(I0)') i
    end do
    close(iounit)
  end subroutine cba_save_dofs

  ! -- Save Reduction Basis branch --------------------------------------
  function cba_want_basis() result(r) bind(c, name="cba_want_basis")
    integer(c_int) :: r
    logical :: found
    r = merge(1, 0, listgetlogical(acc_solver % values, 'Save Reduction Basis', found))
  end function cba_want_basis

  ! SaveReductionBasis internal (exact transcription).
  subroutine cba_save_basis() bind(c, name="cba_save_basis")
    real(kind=dp) :: val
    integer :: iounit, i, j, k, l
    acc_filename = 'ReductionModes.dat'
    call info('CraigBamptonSolver', &
        'Saving the actual component modes to: '//trim(acc_filename), level=7)
    open (newunit=iounit, file=acc_filename)
    do i = 1, acc_mesh % numberofnodes
      j = acc_evar % perm(i)
      if (j == 0) cycle
      do k = 1, acc_nocomp
        write(iounit, '(I0)', advance='NO') i
        do l = 1, acc_dofs
          if (k <= acc_noeigen) then
            val = acc_evar % eigenvectors(k, acc_dofs*(j-1)+l)
          else
            val = acc_evar % constraintmodes(k-acc_noeigen, acc_dofs*(j-1)+l)
          end if
          write(iounit, '(ES16.7)', advance='NO') val
        end do
      end do
      write(iounit, '(A)') ' '
    end do
    close(iounit)
  end subroutine cba_save_basis

  ! -- Save Reduction Matrix branch --------------------------------------
  function cba_want_matrix() result(r) bind(c, name="cba_want_matrix")
    integer(c_int) :: r
    logical :: found, savethis
    savethis = listgetlogical(acc_solver % values, 'Save Reduction Matrix', found)
    if (.not. found) savethis = .true.
    r = merge(1, 0, savethis)
  end function cba_want_matrix

  function cba_ncomp() result(n) bind(c, name="cba_ncomp")
    integer(c_int) :: n
    n = acc_nocomp
  end function cba_ncomp

  function cba_nrows() result(n) bind(c, name="cba_nrows")
    integer(c_int) :: n
    n = acc_a % numberofrows
  end function cba_nrows

  ! Per-matrix setup (exact transcription of the MatrixNo branch head).
  ! Returns 1 when the mass pass must CYCLE (bulk mass values absent,
  ! after the same Warn call), else 0.
  function cba_matrix_begin(mno) result(skip) bind(c, name="cba_matrix_begin")
    integer(c_int), value :: mno
    integer(c_int) :: skip
    acc_matrixno = mno
    if (mno == 1) then
      call info('CraigBamptonSolver', 'Generating the reduced stiffness matrix', level=5)
      acc_a % values => acc_a % bulkvalues
      acc_filename = 'ReducedStiff.dat'
      skip = 0
    else
      if (.not. associated(acc_a % bulkmassvalues)) then
        call warn('CraigBamptonSolver', 'Bulk mass values not present!')
        skip = 1
        return
      end if
      call info('CraigBamptonSolver', 'Generating the reduced mass matrix', level=5)
      acc_a % values => acc_a % bulkmassvalues
      acc_filename = 'ReducedMass.dat'
      skip = 0
    end if
  end function cba_matrix_begin

  ! One (k, l) reduction pair (exact transcription of the inner loop
  ! body): fill x with mode k, multiply, dot with mode l. `x` aliases
  ! EVar % Values exactly as in the original; `ahat`/`ax` are the
  ! Rust-held work arrays (same extents the original ALLOCATEs).
  subroutine cba_mode_pair(k, l, ncomp, n, ahat, ax) bind(c, name="cba_mode_pair")
    integer(c_int), value :: k, l, ncomp, n
    real(kind=dp) :: ahat(ncomp, ncomp), ax(n)
    real(kind=dp), pointer :: x(:)
    x => acc_evar % values
    if (k <= acc_noeigen) then
      x = acc_evar % eigenvectors(k, :)
    else
      x = acc_evar % constraintmodes(k-acc_noeigen, :)
    end if
    call matrixvectormultiply(acc_a, x, ax)
    if (l <= acc_noeigen) then
      x = acc_evar % eigenvectors(l, :)
    else
      x = acc_evar % constraintmodes(l-acc_noeigen, :)
    end if
    ahat(l, k) = sum(x * ax)
  end subroutine cba_mode_pair

  ! SaveReductionMatrix internal (exact transcription).
  subroutine cba_save_matrix(ncomp, ahat) bind(c, name="cba_save_matrix")
    integer(c_int), value :: ncomp
    real(kind=dp) :: ahat(ncomp, ncomp)
    integer :: iounit, k
    call info('CraigBamptonSolver', &
        'Saving the reduction matrix to: '//trim(acc_filename), level=7)
    open (newunit=iounit, file=acc_filename)
    do k = 1, acc_nocomp
      write(iounit, *) ahat(k, :)
    end do
    close(iounit)
  end subroutine cba_save_matrix

  ! -- Calculate Matrix Norm branch --------------------------------------
  function cba_want_norm() result(r) bind(c, name="cba_want_norm")
    integer(c_int) :: r
    logical :: found
    r = merge(1, 0, listgetlogical(acc_solver % values, 'Calculate Matrix Norm', found))
  end function cba_want_norm

  ! FrobeniusNorm internal (exact transcription) + the norm reporting
  ! block (exact transcription, including the stiffness-only writeback
  ! to Solver % Variable).
  subroutine cba_norm(mno, ncomp, ahat) bind(c, name="cba_norm")
    integer(c_int), value :: mno, ncomp
    real(kind=dp) :: ahat(ncomp, ncomp)
    real(kind=dp) :: norm
    integer :: i, j
    norm = 0.0_dp
    do i = 1, acc_nocomp
      do j = 1, acc_nocomp
        norm = norm + ahat(i, j)**2
      end do
    end do
    norm = sqrt(norm)
    if (mno == 1) then
      acc_solver % variable % values = norm
      acc_solver % variable % norm = norm
    end if
    if (mno == 1) then
      write(message, '(A,ES15.6)') 'Matrix norm for reduced stiffness matrix: ', norm
    else
      write(message, '(A,ES15.6)') 'Matrix norm for reduced mass matrix: ', norm
    end if
    call info('CraigBamptonSolver', message, level=5)
  end subroutine cba_norm

  ! -- shared tail --------------------------------------------------------
  subroutine cba_restore() bind(c, name="cba_restore")
    acc_a % values => acc_savevalues
  end subroutine cba_restore

  subroutine cba_info_done() bind(c, name="cba_info_done")
    call info('CraigBamptonSolver', 'All done for now', level=5)
  end subroutine cba_info_done

end module craigbamptonsolver_acc
