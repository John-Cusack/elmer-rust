! RateOfChange port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/RateOfChange.F90, compiled against the built
! .mod files of the target configure (ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module mirrors the original's SAVE state (call counter excluded: the
! Rust side owns ThisSolverCalls) with identical lifetimes.
!
! Order-sensitive reductions (SUM in the local matrix, WRITE/Info timing
! lines) stay in Fortran wrappers so bit-identity holds by construction
! wherever the call sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module rateofchange_acc
  use iso_c_binding
  use types
  use lists
  use integration
  use elementdescription
  use solverbasics
  use solvecore
  use boundaryconditionutils
  use generalutils
  use coordinatesystems
  use messages
  implicit none

  type :: roc_varptr_t
    type(variable_t), pointer :: var => null()
  end type roc_varptr_t

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(matrix_t), pointer, save :: acc_stiff => null()
  type(valuelist_t), pointer, save :: acc_material => null()
  type(nodes_t), target, save :: acc_nodes
  type(gaussintegrationpoints_t), save :: acc_ip
  type(roc_varptr_t), allocatable, save :: acc_u0var(:)
  real(kind=dp), pointer, save :: acc_udot(:) => null()
  integer, pointer, save :: acc_udot_perm(:) => null()
  real(kind=dp), pointer, save :: acc_force(:) => null()
  real(kind=dp), pointer, save :: acc_flow(:) => null()
  integer, pointer, save :: acc_flow_perm(:) => null()
  real(kind=dp), save :: acc_norm = 0.0_dp
  real(kind=dp), save :: acc_at = 0.0_dp
  character(len=MAX_NAME_LEN), save :: acc_adv_var = ''
  character(len=MAX_NAME_LEN), save :: acc_adv_flag = ''
  character(len=MAX_NAME_LEN), save :: acc_eq_name = ''
  integer, save :: acc_dofs = 0
  integer, save :: acc_velcomp = 0
  integer, save :: acc_localnodes = 0

contains

  ! -- entry prologue ---------------------------------------------------
  subroutine roc_info_start() bind(c, name="roc_info_start")
    call Info('RateOfChangeSolver', '')
    call Info('RateOfChangeSolver', 'Starting initialization...')
    call Info('RateOfChangeSolver', '')
  end subroutine roc_info_start

  subroutine roc_begin(h_model, h_solver, nmax, dofs, nbulk) &
      bind(c, name="roc_begin")
    type(c_ptr), value :: h_model, h_solver
    integer(c_int) :: nmax, dofs, nbulk
    integer :: istat
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    acc_udot => acc_solver % variable % values
    acc_udot_perm => acc_solver % variable % perm
    acc_localnodes = acc_model % numberofnodes
    acc_stiff => acc_solver % matrix
    acc_force => acc_stiff % rhs
    acc_norm = acc_solver % variable % norm
    acc_dofs = acc_solver % variable % dofs
    nmax = acc_model % maxelementnodes
    dofs = acc_dofs
    nbulk = acc_model % numberofbulkelements
    if (associated(acc_nodes % x)) then
      if (size(acc_nodes % x) < nmax) &
          deallocate(acc_nodes % x, acc_nodes % y, acc_nodes % z)
    end if
    if (.not. associated(acc_nodes % x)) then
      allocate(acc_nodes % x(nmax), acc_nodes % y(nmax), &
          acc_nodes % z(nmax), stat=istat)
      if (istat /= 0) call Fatal('RateOfChangeSolver', &
          'Memory allocation error, aborting.')
    end if
  end subroutine roc_begin

  ! -- advection-variable resolution (original lines: Advection Variable
  !    block with the scalar/vector branch and Fatals) -------------------
  subroutine roc_setup_u0() bind(c, name="roc_setup_u0")
    integer :: i
    character(len=MAX_NAME_LEN) :: varname
    acc_adv_var = ListGetString(acc_solver % values, 'Advection Variable')
    if (allocated(acc_u0var)) then
      if (size(acc_u0var) /= acc_dofs) deallocate(acc_u0var)
    end if
    if (.not. allocated(acc_u0var)) allocate(acc_u0var(acc_dofs))
    if (acc_dofs == 1) then
      acc_u0var(1) % var => VariableGet(acc_solver % mesh % variables, &
          acc_adv_var)
      if (associated(acc_u0var(1) % var)) then
        if (acc_dofs /= acc_u0var(1) % var % dofs) &
            call Fatal('RateOfChangeSolver', &
            'Solver Variable and Advection Variable Dofs do not equal')
      else
        call Fatal('RateOfChangeSolver', &
            'Specified Advection Variable does not exist')
      end if
    else
      do i = 1, acc_dofs
        write(varname,'(A,A,I1)') trim(acc_adv_var),' ',i
        acc_u0var(i) % var => VariableGet(acc_solver % mesh % variables, &
            varname)
        if (.not. associated(acc_u0var(i) % var)) &
            call Fatal('RateOfChangeSolver', &
            'Specified Advection Variable component does not exist')
      end do
    end if
  end subroutine roc_setup_u0

  ! -- advection-velocity resolution (Advection flag block with the
  !    computed/constant branch, Warn and Fatals) ------------------------
  subroutine roc_setup_vel() bind(c, name="roc_setup_vel")
    type(variable_t), pointer :: flowsol
    acc_adv_flag = ListGetString(acc_solver % values, 'Advection')
    if (acc_adv_flag == 'computed') then
      flowsol => VariableGet(acc_solver % mesh % variables, &
          'Flow Solution')
      if (associated(flowsol)) then
        acc_flow_perm => flowsol % perm
        acc_flow => flowsol % values
        acc_velcomp = flowsol % dofs - 1
        if (acc_velcomp /= CoordinateSystemDimension()) &
            call Warn('RateOfChangeSolver', &
            'Coordinate system and Advection Velocity dimensions unequal')
      else
        call Fatal('RateOfChangeSolver', &
            'Flow Solution does not exist')
      end if
    else
      if (acc_adv_flag /= 'constant') call Fatal('RateOfChangeSolver', &
          'Advection flag should be either "computed" or "constant"')
      acc_velcomp = CoordinateSystemDimension()
    end if
  end subroutine roc_setup_vel

  ! -- assembly prologue (timing start, matrix zero, equation name) ------
  subroutine roc_init_zero() bind(c, name="roc_init_zero")
    acc_at = CPUTime()
    call InitializeToZero(acc_stiff, acc_force)
    acc_eq_name = ListGetString(acc_solver % values, 'Equation')
  end subroutine roc_init_zero

  ! -- one bulk element: gather, local vectors, material, velocity,
  !    LocalMatrix, global update (original element loop body) -----------
  subroutine roc_assemble_elem(t, ld, n, ncomp, stiff, force, localu0, &
      v1, v2, v3, skip) bind(c, name="roc_assemble_elem")
    integer(c_int), value :: t, ld, n, ncomp
    real(c_double) :: stiff(ld, ld), force(ld), localu0(ld)
    real(c_double) :: v1(n), v2(n), v3(n)
    integer(c_int) :: skip
    integer :: i, j, k, p, body_id, nn
    acc_elem => acc_solver % mesh % elements(t)
    if (.not. CheckElementEquation(acc_model, acc_elem, acc_eq_name)) then
      skip = 1
      return
    end if
    skip = 0

    nn = acc_elem % type % numberofnodes
    acc_nodes % x(1:nn) = acc_solver % mesh % nodes % x( &
        acc_elem % nodeindexes)
    acc_nodes % y(1:nn) = acc_solver % mesh % nodes % y( &
        acc_elem % nodeindexes)
    acc_nodes % z(1:nn) = acc_solver % mesh % nodes % z( &
        acc_elem % nodeindexes)

    do p = 1, ncomp
      do i = 1, nn
        j = acc_udot_perm(acc_elem % nodeindexes(i))
        k = acc_u0var(p) % var % perm(acc_elem % nodeindexes(i))
        if (j > 0 .and. k > 0) then
          localu0((i-1)*ncomp + p) = acc_u0var(p) % var % values(k)
        else
          call Fatal('RateOfChangeSolver', &
              'Nonmatching variable permutations')
        end if
      end do
    end do

    body_id = acc_elem % bodyid
    k = ListGetInteger(acc_model % bodies(body_id) % values, 'Material')
    acc_material => acc_model % materials(k) % values

    if (acc_adv_flag == 'constant') then
      select case (acc_velcomp)
        case(1)
        v1(1:nn) = ListGetReal(acc_material, 'Advection Velocity 1', &
            nn, acc_elem % nodeindexes)
        v2 = 0.0_dp
        v3 = 0.0_dp
        case(2)
        v1(1:nn) = ListGetReal(acc_material, 'Advection Velocity 1', &
            nn, acc_elem % nodeindexes)
        v2(1:nn) = ListGetReal(acc_material, 'Advection Velocity 2', &
            nn, acc_elem % nodeindexes)
        v3 = 0.0_dp
        case(3)
        v1(1:nn) = ListGetReal(acc_material, 'Advection Velocity 1', &
            nn, acc_elem % nodeindexes)
        v2(1:nn) = ListGetReal(acc_material, 'Advection Velocity 2', &
            nn, acc_elem % nodeindexes)
        v3(1:nn) = ListGetReal(acc_material, 'Advection Velocity 3', &
            nn, acc_elem % nodeindexes)
      end select
    else
      do i = 1, nn
        j = acc_flow_perm(acc_elem % nodeindexes(i))
        select case (acc_velcomp)
          case(1)
          v1(i) = acc_flow((j-1)*2+1)
          v2(i) = 0.0_dp
          v3(i) = 0.0_dp
          case (2)
          v1(i) = acc_flow((j-1)*3+1)
          v2(i) = acc_flow((j-1)*3+2)
          v3(i) = 0.0_dp
          case (3)
          v1(i) = acc_flow((j-1)*4+1)
          v2(i) = acc_flow((j-1)*4+2)
          v3(i) = acc_flow((j-1)*4+3)
        end select
      end do
    end if

    call roc_local_matrix(ld, nn, ncomp, stiff, force, localu0, v1, v2, v3)

    call UpdateGlobalEquations(acc_stiff, stiff, acc_force, force, &
        nn, ncomp, acc_udot_perm(acc_elem % nodeindexes))
  end subroutine roc_assemble_elem

  ! -- the LocalMatrix internal procedure, statement for statement --------
  subroutine roc_local_matrix(ld, n, ncomp, stiff, force, localu0, &
      v1, v2, v3)
    integer :: ld, n, ncomp
    real(kind=dp) :: stiff(ld, ld), force(ld), localu0(ld)
    real(kind=dp) :: v1(*), v2(*), v3(*)
    real(kind=dp), dimension(3) :: vel
    integer :: i, p, q, tt, dim, ndc
    real(kind=dp) :: basis(n), dbasisdx(n, 3)
    real(kind=dp) :: sqrtele, u, v, w, s, temp
    logical :: stat
    type(gaussintegrationpoints_t) :: integ

    dim = CoordinateSystemDimension()
    ndc = n*ncomp
    force(1:ld) = 0.0_dp
    stiff(1:ld, 1:ld) = 0.0_dp
    integ = GaussPoints(acc_elem)
    do tt = 1, integ % n
      u = integ % u(tt)
      v = integ % v(tt)
      w = integ % w(tt)
      s = integ % s(tt)
      stat = ElementInfo(acc_elem, acc_nodes, u, v, w, sqrtele, &
          basis, dbasisdx)

      s = s * sqrtele
      vel(1) = sum(v1(1:n) * basis(1:n))
      vel(2) = sum(v2(1:n) * basis(1:n))
      vel(3) = sum(v3(1:n) * basis(1:n))

      do i = 1, ncomp
        do p = 1, n
          do q = 1, n
            stiff((p-1)*ncomp+i, (q-1)*ncomp+i) = &
                stiff((p-1)*ncomp+i, (q-1)*ncomp+i) &
                + basis(q) * basis(p) * s

            temp = sum(vel(1:dim) * dbasisdx(q, 1:dim))
            force((p-1)*ncomp + i) = &
                force((p-1)*ncomp + i) - &
                temp * localu0((q-1)*ncomp + i) * basis(p) * s
          end do
        end do
      end do
    end do
  end subroutine roc_local_matrix

  ! -- Dirichlet boundaries over components -------------------------------
  subroutine roc_dirichlet() bind(c, name="roc_dirichlet")
    integer :: i
    character(len=MAX_NAME_LEN) :: varname
    do i = 1, acc_dofs
      varname = ComponentName(acc_solver % variable, i)
      call SetDirichletBoundaries(acc_model, acc_stiff, acc_force, &
          varname, i, acc_dofs, acc_udot_perm)
    end do
  end subroutine roc_dirichlet

  ! -- solve with timing lines --------------------------------------------
  subroutine roc_solve() bind(c, name="roc_solve")
    real(kind=dp) :: st
    acc_at = CPUTime() - acc_at
    st = CPUTime()
    call SolveSystem(acc_stiff, ParMatrix, acc_force, acc_udot, &
        acc_norm, acc_dofs, acc_solver)
    st = CPUTime() - st
    write(Message,'(a,F8.2)') ' Assembly: (s)', acc_at
    call Info('RateOfChangeSolver', Message)
    write(Message,'(a,F8.2)') ' Solve:    (s)', st
    call Info('RateOfChangeSolver', Message)
  end subroutine roc_solve

  subroutine roc_fatal_alloc() bind(c, name="roc_fatal_alloc")
    call Fatal('RateOfChangeSolver', 'Memory allocation error, aborting.')
  end subroutine roc_fatal_alloc

end module rateofchange_acc
