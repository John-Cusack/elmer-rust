! Port of fem/src/modules/OdeSolver.F90 (Elmer fork) for lane T021-B,
! at the solver-plugin boundary (whole `.so` -> Rust cdylib).
!
! Statement-level BIND(C) wrappers for the Elmer calls in OdeSolver.F90,
! compiled against the built .mod files of the target configure. Every
! wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars). All Fortran-object traffic — TYPE(...)
! handles, OPTIONAL/POINTER dummies, assumed-shape descriptors, hidden
! lengths — stays inside Fortran, so the Rust side is pure C calls.
!
! The OdeMatrix internal procedure body is transcribed verbatim into
! oda_matrix (same compiler, same statement order): component access,
! Info/Fatal I/O, I2S string building, and the CRS matrix fills all stay
! in Fortran, so bit-identity holds by construction wherever the call
! sequence matches. The IF(.FALSE.) PRINT block is dead code in the
! original and is not transcribed.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module odesolver_acc
  use iso_c_binding
  use types
  use defutils
  use messages
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(valuelist_t), pointer, save :: acc_params => null()

contains

  ! -- entry state ----------------------------------------------------
  subroutine oda_begin(h_solver) bind(c, name="oda_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine oda_begin

  ! -- OdeSolver_init body (verbatim) ----------------------------------
  subroutine oda_init_params() bind(c, name="oda_init_params")
    type(valuelist_t), pointer :: params
    params => GetSolverParams()
    call ListAddNewInteger( params, 'Time derivative order', 2 )
    call ListAddLogical( params, 'Use Global Mass Matrix', .TRUE. )
    call ListAddLogical( params, 'Ode Matrix', .TRUE. )
    call ListAddLogical( params, 'Variable Global', .TRUE. )
  end subroutine oda_init_params

  ! -- OdeSolver body --------------------------------------------------
  subroutine oda_info_start() bind(c, name="oda_info_start")
    call Info('OdeSolver','Solving ordinary differential equation',Level=4)
  end subroutine oda_info_start

  subroutine oda_params_begin() bind(c, name="oda_params_begin")
    acc_params => GetSolverParams()
  end subroutine oda_params_begin

  function oda_maxiter() result(m) bind(c, name="oda_maxiter")
    integer(c_int) :: m
    logical :: found
    m = GetInteger( acc_params,'Nonlinear System Max Iterations', found )
    if (.not. found) m = 1
  end function oda_maxiter

  function oda_timeorder() result(t) bind(c, name="oda_timeorder")
    integer(c_int) :: t
    t = GetInteger( acc_params,'Time Derivative Order' )
  end function oda_timeorder

  subroutine oda_info_iter(iter) bind(c, name="oda_info_iter")
    integer(c_int), value :: iter
    call Info('OdeSolver','Nonlinear iteration: '//I2S(iter),Level=5)
  end subroutine oda_info_iter

  subroutine oda_default_initialize() bind(c, name="oda_default_initialize")
    call DefaultInitialize()
  end subroutine oda_default_initialize

  ! OdeMatrix internal procedure, transcribed verbatim (TimeOrder is the
  ! host variable, passed down by value from the entry wrapper state).
  subroutine oda_matrix(TimeOrder) bind(c, name="oda_matrix")
    integer(c_int), value :: TimeOrder
    type(matrix_t), pointer :: A
    integer :: i,j,Dofs
    real(kind=dp), pointer CONTIG :: SaveValues(:)
    real(kind=dp) :: val
    type(valuelist_t), pointer :: OdeList
    integer, pointer :: ActiveComponents(:)
    character(len=MAX_NAME_LEN) :: str
    logical :: Found

    ActiveComponents => ListGetIntegerArray( acc_params, &
        'Active Components', Found )
    if ( Found ) then
      if ( size( ActiveComponents ) > 1 ) then
        call Fatal('OdeSolver','Currently implemented only for one component!')
      end if
      i = ActiveComponents(1)
      if ( i > CurrentModel % NumberOfComponents ) then
        call Fatal('OdeSolver','Active Component index out of range')
      end if
      OdeList => CurrentModel % Components(i) % Values
      call Info('OdeSolver','Using active component: '//I2S(i),Level=10)
    else
      OdeList => acc_params
      call Info('OdeSolver','Using solver section',Level=10)
    end if

    A => acc_solver % Matrix
    Dofs = acc_solver % Variable % Dofs

    if ( .not. associated( A % MassValues ) ) then
      call Info('OdeSolver','Allocating mass matrix',Level=10)
      allocate( A % MassValues(size(A % Values)) )
      A % MassValues = 0.0_dp
    end if

    if ( TimeOrder >= 2 ) then
      if ( .not. associated( A % DampValues ) ) then
        call Info('OdeSolver','Allocating damping matrix',Level=10)
        allocate( A % DampValues(size(A % Values)) )
        A % DampValues = 0.0_dp
      end if
    end if

    ! Set stiffness matrix values
    do i=1,Dofs
      do j=1,Dofs
        str = 'Stiffness Matrix '//I2S(i)//I2S(j)
        val = ListGetCReal( OdeList, str, Found )
        call CRS_SetMatrixElement( A,i,j,val )
      end do
    end do

    ! Set mass matrix values
    SaveValues => A % Values
    A % Values => A % MassValues
    do i=1,Dofs
      do j=1,Dofs
        ! In Elmer the convention is to call the highest time-derivatie always mass
        ! Even for 1st order PDEs. However, in this solver we like to call the coefficients
        ! mass, damping and spring.
        if ( TimeOrder >= 2 ) then
          str = 'Mass Matrix '//I2S(i)//I2S(j)
        else
          str = 'Damping Matrix '//I2S(i)//I2S(j)
        end if
        val = ListGetCReal( OdeList, str, Found )
        call CRS_SetMatrixElement( A,i,j,val )
      end do
    end do
    A % Values => SaveValues

    if ( TimeOrder >= 2 ) then
      ! Set damp matrix values
      SaveValues => A % Values
      A % Values => A % DampValues
      do i=1,Dofs
        do j=1,Dofs
          str = 'Damping Matrix '//I2S(i)//I2S(j)
          val = ListGetCReal( OdeList, str, Found )
          call CRS_SetMatrixElement( A,i,j,val )
        end do
      end do
      A % Values => SaveValues
    end if

    do i=1,Dofs
      str = 'Force '//I2S(i)
      val = ListGetCReal( OdeList, str, Found )
      A % Rhs(i) = val
    end do
  end subroutine oda_matrix

  subroutine oda_finish_bulk() bind(c, name="oda_finish_bulk")
    call DefaultFinishBulkAssembly()
  end subroutine oda_finish_bulk

  subroutine oda_finish_bnd() bind(c, name="oda_finish_bnd")
    call DefaultFinishBoundaryAssembly()
  end subroutine oda_finish_bnd

  subroutine oda_finish_assembly() bind(c, name="oda_finish_assembly")
    call DefaultFinishAssembly()
  end subroutine oda_finish_assembly

  function oda_solve() result(norm) bind(c, name="oda_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function oda_solve

  function oda_converged() result(r) bind(c, name="oda_converged")
    integer(c_int) :: r
    r = merge(1, 0, acc_solver % Variable % NonlinConverged == 1)
  end function oda_converged

  subroutine oda_info_done() bind(c, name="oda_info_done")
    call Info('OdeSolver','All done',Level=5)
  end subroutine oda_info_done

end module odesolver_acc
