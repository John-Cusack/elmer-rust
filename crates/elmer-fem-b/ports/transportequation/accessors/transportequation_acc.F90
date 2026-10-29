! T021-B TransportEquation port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/TransportEquation.F90, compiled against
! the built .mod files of the target configure (T020 Poisson pattern,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module mirrors TransportEquationSolver's SAVE locals (ThisSolverCalls-
! equivalent counter lives in Rust; LocalStiffMatrix/LocalForce/
! LocalMassMatrix/LocalDampMatrix/V1/V2/V3 ride Rust Vecs; ElementNodes,
! U0Var, Gauss-point states live here) with identical lifetimes.
!
! The internal procedures LocalMatrix/LocalMatrixBoundary/
! Add2ndOrderTime2/AverageAccelerationMethod/ModifyBoundaryValues are
! transcribed here statement for statement (order-sensitive SUM
! reductions and formatted I/O stay in Fortran, so bit-identity holds by
! construction wherever the call sequence matches). The uninitialized
! `Method` SELECT CASE in Add2ndOrderTime2 is transcribed verbatim: its
! only case is DEFAULT, so the garbage value is unobservable.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module transportequation_acc
  use iso_c_binding
  use defutils
  use types
  use lists
  use integration
  use elementdescription
  use solverbasics
  use elementutils
  use messages
  implicit none

  type varptr_t
    type(variable_t), pointer :: var => null()
  end type varptr_t

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(valuelist_t), pointer, save :: acc_mat => null()
  type(varptr_t), save, allocatable :: acc_u0v(:)
  type(variable_t), pointer, save :: acc_udot0 => null()
  type(variable_t), pointer, save :: acc_flow => null()
  integer, pointer, save :: acc_vperm(:) => null()
  real(kind=dp), pointer, save :: acc_flowv(:) => null()
  character(len=MAX_NAME_LEN), save :: acc_advvar = ''
  character(len=MAX_NAME_LEN), save :: acc_eqname = ''
  integer, save :: acc_advflag = 0
  integer, save :: acc_dofs = 1
  integer, save :: acc_velcomps = 0
  logical, save :: acc_allocdone = .false.
  ! ElementNodes SAVE mirror (allocated once at first-mesh sizing).
  type(nodes_t), save :: acc_enodes
  ! LocalMatrix / LocalMatrixBoundary Gauss-point states.
  type(gaussintegrationpoints_t), save :: acc_ip_lm
  type(gaussintegrationpoints_t), save :: acc_ip_lmb
  real(kind=dp), save :: acc_at = 0.0_dp
  real(kind=dp), save :: acc_st = 0.0_dp
  integer, save, allocatable :: acc_di(:)

contains

  ! -- entry state ------------------------------------------------------
  subroutine teq_begin(h_model, h_solver) bind(c, name="teq_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine teq_begin

  subroutine teq_info_init() bind(c, name="teq_info_init")
    call Info('TransportEquationSolver', '')
    call Info('TransportEquationSolver', 'Starting initialization...')
    call Info('TransportEquationSolver', '')
  end subroutine teq_info_init

  function teq_var_dofs() result(n) bind(c, name="teq_var_dofs")
    integer(c_int) :: n
    n = acc_solver % variable % dofs
    acc_dofs = n
  end function teq_var_dofs

  function teq_max_elem_nodes() result(n) bind(c, name="teq_max_elem_nodes")
    integer(c_int) :: n
    n = acc_model % maxelementnodes
  end function teq_max_elem_nodes

  ! -- first-mesh permanent storage + initial conditions (lines 112-236) -
  subroutine teq_init_cond(dt) bind(c, name="teq_init_cond")
    real(c_double), value :: dt
    integer :: i, j, k, p, t, istat
    logical :: gotit
    character(len=MAX_NAME_LEN) :: variablename
    integer, pointer :: uperm(:)
    if (.not. acc_allocdone) then
      allocate(acc_u0v(acc_dofs), stat=istat)
      if (istat /= 0) call Fatal('TransportEquationSolver', &
          'Memory allocation error, aborting.')
      allocate(acc_enodes % x(acc_model % maxelementnodes), &
          acc_enodes % y(acc_model % maxelementnodes), &
          acc_enodes % z(acc_model % maxelementnodes), stat=istat)
      if (istat /= 0) call Fatal('TransportEquationSolver', &
          'Memory allocation error, aborting.')
      acc_allocdone = .true.

    end if
    acc_advvar = ListGetString(acc_solver % values, &
        'Advection Variable', gotit)
    if (acc_dofs == 1) then
      acc_u0v(1) % var => VariableGet(acc_solver % mesh % variables, &
          acc_advvar)
      if (associated(acc_u0v(1) % var)) then
        if (acc_dofs /= acc_u0v(1) % var % dofs) &
            call Fatal('TransportEquationSolver', &
            'Solver Variable and Advection Variable Dofs do not equal')
      else
        call Fatal('TransportEquationSolver', &
            'Specified Advection Variable does not exist')
      end if
    else
      do i = 1, acc_dofs
        write(variablename,'(A,A,I1)') trim(acc_advvar),' ', i
        acc_u0v(i) % var => VariableGet(acc_solver % mesh % variables, &
            variablename)
        if (.not. associated(acc_u0v(i) % var)) &
            call Fatal('TransportEquationSolver', &
            'Specified Advection Variable component does not exist')
      end do
    end if
    variablename = ListGetString(acc_solver % values, &
        'Rate Of Change Equation Variable')
    acc_udot0 => VariableGet(acc_solver % mesh % variables, variablename)
    acc_solver % dt = 2.0_dp * dt
    uperm => acc_solver % variable % perm
    do i = 1, acc_model % numberofnodes
      j = uperm(i)
      k = acc_u0v(1) % var % perm(i)
      p = acc_udot0 % perm(i)
      if (j > 0 .and. k > 0 .and. p > 0) then
        do t = 1, acc_dofs
          acc_solver % variable % values((j-1)*acc_dofs + t) = &
              acc_u0v(t) % var % values(k)
          acc_solver % variable % prevvalues((j-1)*acc_dofs + t,1) = &
              acc_udot0 % values((p-1)*acc_dofs + t)
          acc_solver % variable % prevvalues((i-1)*acc_dofs + t,2) = 0.0_dp
        end do
      else
        call Fatal('TransportEquationSolver', &
            'Nonmatching variable permutations')
      end if
    end do
  end subroutine teq_init_cond

  subroutine teq_at_start() bind(c, name="teq_at_start")
    acc_at = CPUTime()
  end subroutine teq_at_start
  subroutine teq_at_done() bind(c, name="teq_at_done")
    acc_at = CPUTime() - acc_at
  end subroutine teq_at_done


  subroutine teq_init_zero() bind(c, name="teq_init_zero")
    call InitializeToZero(acc_solver % matrix, acc_solver % matrix % rhs)
  end subroutine teq_init_zero

  subroutine teq_init_timestep() bind(c, name="teq_init_timestep")
    call InitializeTimestep(acc_solver)
  end subroutine teq_init_timestep

  subroutine teq_eqname() bind(c, name="teq_eqname")
    acc_eqname = ListGetString(acc_solver % values, 'Equation')
  end subroutine teq_eqname

  ! -- advection setup (lines 184-207) ------------------------------------
  subroutine teq_adv_setup(flag, velcomps) bind(c, name="teq_adv_setup")
    integer(c_int) :: flag, velcomps
    character(len=MAX_NAME_LEN) :: advflag
    advflag = ListGetString(acc_solver % values, 'Advection')
    if (advflag == 'computed') then
      acc_flow => VariableGet(acc_solver % mesh % variables, &
          'Flow Solution')
      if (associated(acc_flow)) then
        acc_vperm => acc_flow % perm
        acc_flowv => acc_flow % values
        velcomps = acc_flow % dofs - 1
        if (velcomps /= CoordinateSystemDimension()) &
            call Warn('TransportEquationSolver', &
            'Coordinate system and Advection Velocity dimensions unequal')
        flag = 2
      else
        call Fatal('TransportEquationSolver', &
            'Flow Solution does not exist')
      end if
    else
      if (advflag /= 'constant') call Fatal('TransportEquationSolver', &
          'Advection flag should be either "computed" or "constant"')
      velcomps = CoordinateSystemDimension()
      flag = 1
    end if
    acc_advflag = flag
    acc_velcomps = velcomps
  end subroutine teq_adv_setup

  function teq_num_bulk() result(n) bind(c, name="teq_num_bulk")
    integer(c_int) :: n
    n = acc_solver % mesh % numberofbulkelements
  end function teq_num_bulk

  function teq_num_bnd() result(n) bind(c, name="teq_num_bnd")
    integer(c_int) :: n
    n = acc_solver % mesh % numberofboundaryelements
  end function teq_num_bnd

  ! -- bulk-element context (lines 249-263) --------------------------------
  subroutine teq_begin_bulk(t, n, proceed) bind(c, name="teq_begin_bulk")
    integer(c_int), value :: t
    integer(c_int) :: n, proceed
    integer :: body_id, k
    acc_elem => acc_solver % mesh % elements(t)
    if (.not. CheckElementEquation(acc_model, acc_elem, acc_eqname)) then
      proceed = 0
      return
    end if
    n = acc_elem % type % numberofnodes
    body_id = acc_elem % bodyid
    k = ListGetInteger(acc_model % bodies(body_id) % values, 'Material')
    acc_mat => acc_model % materials(k) % values
    proceed = 1
  end subroutine teq_begin_bulk

  ! -- boundary-element context (lines 328-353) -----------------------------
  subroutine teq_begin_bnd(t, n, proceed) bind(c, name="teq_begin_bnd")
    integer(c_int), value :: t
    integer(c_int) :: n, proceed
    type(element_t), pointer :: parent
    logical :: stat
    integer :: body_id, k
    acc_elem => acc_solver % mesh % elements(t)
    if (size(acc_elem % nodeindexes) > 1) then
      parent => acc_elem % boundaryinfo % left
      stat = associated(parent)
      if (stat) stat = all(acc_solver % variable % perm(parent % nodeindexes) > 0)
      if (.not. stat) then
        parent => acc_elem % boundaryinfo % right
        stat = associated(parent)
        if (stat) stat = all(acc_solver % variable % perm(parent % nodeindexes) > 0)
        if (.not. stat) call Fatal('TransportEquationSolver', &
            'No parent element can be found for given boundary element')
      end if
      if (.not. CheckElementEquation(acc_model, parent, acc_eqname)) then
        proceed = 0
        return
      end if
      body_id = parent % bodyid
      k = ListGetInteger(acc_model % bodies(body_id) % values, 'Material')
      acc_mat => acc_model % materials(k) % values
      n = acc_elem % type % numberofnodes
      proceed = 1
    else
      proceed = 0
    end if
  end subroutine teq_begin_bnd

  ! -- ElementNodes fill (lines 257-259, 355-357) ---------------------------
  subroutine teq_elem_nodes(n) bind(c, name="teq_elem_nodes")
    integer(c_int), value :: n
    acc_enodes % x(1:n) = acc_solver % mesh % nodes % x(acc_elem % nodeindexes)
    acc_enodes % y(1:n) = acc_solver % mesh % nodes % y(acc_elem % nodeindexes)
    acc_enodes % z(1:n) = acc_solver % mesh % nodes % z(acc_elem % nodeindexes)
  end subroutine teq_elem_nodes

  ! -- velocity fill, bulk and boundary shared (lines 265-304, 359-398) -----
  subroutine teq_fill_vel(n, nn, v1, v2, v3) bind(c, name="teq_fill_vel")
    integer(c_int), value :: n, nn
    real(c_double) :: v1(nn), v2(nn), v3(nn)
    integer :: i, j
    if (acc_advflag == 1) then
      select case (acc_velcomps)
      case (1)
        v1(1:n) = ListGetReal(acc_mat, 'Advection Velocity 1', n, &
            acc_elem % nodeindexes)
        v2 = 0.0_dp
        v3 = 0.0_dp
      case (2)
        v1(1:n) = ListGetReal(acc_mat, 'Advection Velocity 1', n, &
            acc_elem % nodeindexes)
        v2(1:n) = ListGetReal(acc_mat, 'Advection Velocity 2', n, &
            acc_elem % nodeindexes)
        v3 = 0.0_dp
      case (3)
        v1(1:n) = ListGetReal(acc_mat, 'Advection Velocity 1', n, &
            acc_elem % nodeindexes)
        v2(1:n) = ListGetReal(acc_mat, 'Advection Velocity 2', n, &
            acc_elem % nodeindexes)
        v3(1:n) = ListGetReal(acc_mat, 'Advection Velocity 3', n, &
            acc_elem % nodeindexes)
      end select
    else
      do i = 1, n
        j = acc_vperm(acc_elem % nodeindexes(i))
        select case (acc_velcomps)
        case (1)
          v1(i) = acc_flowv((j-1)*2+1)
          v2(i) = 0.0_dp
          v3(i) = 0.0_dp
        case (2)
          v1(i) = acc_flowv((j-1)*3+1)
          v2(i) = acc_flowv((j-1)*3+2)
          v3(i) = 0.0_dp
        case (3)
          v1(i) = acc_flowv((j-1)*4+1)
          v2(i) = acc_flowv((j-1)*4+2)
          v3(i) = acc_flowv((j-1)*4+3)
        end select
      end do
    end if
  end subroutine teq_fill_vel

  ! -- LocalMatrix/LocalMatrixBoundary zeroing (lines 502-505, 571-575) -----
  subroutine teq_lm_zero(ld, stiff, mass, damp, force) bind(c, name="teq_lm_zero")
    integer(c_int), value :: ld
    real(c_double) :: stiff(ld,ld), mass(ld,ld), damp(ld,ld), force(ld)
    force = 0.0_dp
    stiff = 0.0_dp
    damp = 0.0_dp
    mass = 0.0_dp
  end subroutine teq_lm_zero

  ! -- LocalMatrix integration ----------------------------------------------
  function teq_lm_ip_begin() result(n) bind(c, name="teq_lm_ip_begin")
    integer(c_int) :: n
    acc_ip_lm = GaussPoints(acc_elem)
    n = acc_ip_lm % n
  end function teq_lm_ip_begin

  subroutine teq_lm_ip_point(t, u, v, w, s) bind(c, name="teq_lm_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip_lm % u(t)
    v = acc_ip_lm % v(t)
    w = acc_ip_lm % w(t)
    s = acc_ip_lm % s(t)
  end subroutine teq_lm_ip_point

  subroutine teq_lm_eleminfo(u, v, w, n, basis, dbasis, detj, stat) &
      bind(c, name="teq_lm_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: n
    real(c_double) :: basis(n), dbasis(n,3), detj
    integer(c_int) :: stat
    logical :: ok

    ok = ElementInfo(acc_elem, acc_enodes, u, v, w, detj, &
        basis(1:n), dbasis(1:n,1:3))
    stat = merge(1, 0, ok)
  end subroutine teq_lm_eleminfo

  ! Exact transcription of the LocalMatrix per-IP body (lines 522-545).
  subroutine teq_lm_ip(ld, n, vc, dim, basis, dbasis, sdet, v1, v2, v3, &
      stiff, mass) bind(c, name="teq_lm_ip")
    integer(c_int), value :: ld, n, vc, dim
    real(c_double) :: basis(n), dbasis(n,3)
    real(c_double) :: v1(n), v2(n), v3(n)
    real(c_double) :: stiff(ld,ld), mass(ld,ld)
    real(c_double), value :: sdet
    real(kind=dp) :: vel(3), temp1, temp2
    integer :: i, p, q
    vel(1) = sum(v1(1:n) * basis(1:n))
    vel(2) = sum(v2(1:n) * basis(1:n))
    vel(3) = sum(v3(1:n) * basis(1:n))
    do i = 1, vc
      do p = 1, n
        do q = 1, n
          mass((p-1)*vc+i, (q-1)*vc+i) = &
              mass((p-1)*vc+i, (q-1)*vc+i) + basis(p) * basis(q) * sdet
          temp1 = sum(vel(1:dim) * dbasis(q,1:dim))
          temp2 = sum(vel(1:dim) * dbasis(p,1:dim))
          stiff((p-1)*vc+i, (q-1)*vc+i) = &
              stiff((p-1)*vc+i, (q-1)*vc+i) + temp1 * temp2 * sdet
        end do
      end do
    end do
  end subroutine teq_lm_ip

  ! -- LocalMatrixBoundary integration ---------------------------------------
  function teq_lmb_ip_begin() result(n) bind(c, name="teq_lmb_ip_begin")
    integer(c_int) :: n
    acc_ip_lmb = GaussPoints(acc_elem)
    n = acc_ip_lmb % n
  end function teq_lmb_ip_begin

  subroutine teq_lmb_ip_point(t, u, v, w, s) bind(c, name="teq_lmb_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip_lmb % u(t)
    v = acc_ip_lmb % v(t)
    w = acc_ip_lmb % w(t)
    s = acc_ip_lmb % s(t)
  end subroutine teq_lmb_ip_point

  subroutine teq_lmb_eleminfo(u, v, w, n, basis, dbasis, detj, stat) &
      bind(c, name="teq_lmb_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: n
    real(c_double) :: basis(n), dbasis(n,3), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_elem, acc_enodes, u, v, w, detj, &
        basis(1:n), dbasis(1:n,1:3))
    stat = merge(1, 0, ok)
  end subroutine teq_lmb_eleminfo

  ! Exact transcription of the LocalMatrixBoundary per-IP body (594-610).
  subroutine teq_lmb_ip(ld, n, vc, dim, basis, dbasis, sdet, u, v, w, &
      v1, v2, v3, damp) bind(c, name="teq_lmb_ip")
    integer(c_int), value :: ld, n, vc, dim
    real(c_double) :: basis(n), dbasis(n,3)
    real(c_double), value :: sdet, u, v, w
    real(c_double) :: v1(n), v2(n), v3(n)
    real(c_double) :: damp(ld,ld)
    real(kind=dp) :: vel(3), normal(3), normvel
    integer :: i, p, q
    vel(1) = sum(v1(1:n) * basis(1:n))
    vel(2) = sum(v2(1:n) * basis(1:n))
    vel(3) = sum(v3(1:n) * basis(1:n))
    normal = Normalvector(acc_elem, acc_enodes, u, v, .true.)
    normvel = sum(normal(1:dim) * vel(1:dim))
    do i = 1, vc
      do p = 1, n
        do q = 1, n
          damp((p-1)*vc+i, (q-1)*vc+i) = &
              damp((p-1)*vc+i, (q-1)*vc+i) + &
              basis(p) * basis(q) * normvel * sdet
        end do
      end do
    end do
  end subroutine teq_lmb_ip

  ! -- Add2ndOrderTime2 (lines 668-741) + AverageAccelerationMethod ---------
  subroutine teq_add2nd(ld, mass, damp, stiff, force, dt, n, dofs) &
      bind(c, name="teq_add2nd")
    integer(c_int), value :: ld, n, dofs
    real(c_double) :: mass(ld,ld), damp(ld,ld), stiff(ld,ld), force(ld)
    real(c_double), value :: dt
    logical :: gotit
    integer :: i, j, k, l
    character(len=MAX_NAME_LEN) :: method
    real(kind=dp) :: x(dofs*n), v(dofs*n), a(dofs*n)
    if (acc_solver % matrix % lumped) then
      call Fatal('TransportEquationSolver', &
          '"Lumped" option is not available')
    end if
    do i = 1, n
      do j = 1, dofs
        k = dofs * (i-1) + j
        l = dofs * (acc_solver % variable % perm(acc_elem % nodeindexes(i))-1) + j
        select case (method)
        case default
          x(k) = acc_solver % variable % prevvalues(l,3)
          v(k) = acc_solver % variable % prevvalues(l,4)
          a(k) = acc_solver % variable % prevvalues(l,5)
        end select
        acc_solver % matrix % force(l,1) = acc_solver % matrix % force(l,1) + &
            force(k)
      end do
    end do
    call teq_avgacc(n*dofs, dt, mass, damp, stiff, force, x, v, a)
  end subroutine teq_add2nd

  ! Exact transcription of AverageAccelerationMethod (lines 746-778).
  subroutine teq_avgacc(m, dt, mass, damp, stiff, force, x, v, a)
    integer :: m
    real(kind=dp) :: force(:), x(:), v(:), a(:), dt
    real(kind=dp) :: alpha, beta, gamma
    real(kind=dp) :: mass(:,:), damp(:,:), stiff(:,:)
    integer :: i, j
    real(kind=dp) :: s
    alpha = 0.0_dp
    gamma = 0.5_dp - alpha
    beta = (1.0_dp - alpha)**2 / 4.0_dp
    do i = 1, m
      s = 0.0_dp
      do j = 1, m
        s = s + ((1.0_dp - alpha) / (beta*dt**2)) * mass(i,j) * x(j)
        s = s + ((1.0_dp - alpha) / (beta*dt)) * mass(i,j) * v(j)
        s = s + (gamma / (beta*dt)) * damp(i,j) * x(j)
        s = s - stiff(i,j) * x(j)
        stiff(i,j) = stiff(i,j) + &
            ((1.0_dp - alpha) / (beta*dt**2)) * mass(i,j) + &
            (gamma / (beta*dt)) * damp(i,j)
      end do
      force(i) = s
    end do
  end subroutine teq_avgacc

  ! -- global update (lines 320-322, 408-410) --------------------------------
  subroutine teq_update_global(ld, n, dofs, stiff, force) &
      bind(c, name="teq_update_global")
    integer(c_int), value :: ld, n, dofs
    real(c_double) :: stiff(ld,ld), force(ld)
    integer :: istat
    if (.not. allocated(acc_di) .or. size(acc_di) < n) then
      if (allocated(acc_di)) deallocate(acc_di)
      allocate(acc_di(n), stat=istat)
      if (istat /= 0) call Fatal('TransportEquationSolver', &
          'Memory allocation error, aborting.')
    end if
    acc_di(1:n) = acc_solver % variable % perm(acc_elem % nodeindexes)
    call UpdateGlobalEquations(acc_solver % matrix, stiff(1:ld,1:ld), &
        acc_solver % matrix % rhs, force(1:ld), n, dofs, acc_di(1:n))
  end subroutine teq_update_global

  subroutine teq_finish() bind(c, name="teq_finish")
    call FinishAssembly(acc_solver, acc_solver % matrix % rhs)
  end subroutine teq_finish

  subroutine teq_dirichlet() bind(c, name="teq_dirichlet")
    call DefaultDirichletBCs()
  end subroutine teq_dirichlet

  function teq_solve() result(norm) bind(c, name="teq_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function teq_solve

  subroutine teq_st_start() bind(c, name="teq_st_start")
    acc_st = CPUTime()
  end subroutine teq_st_start

  ! Timing report (lines 439-442): CPU-time values, normed out of the proof
  ! compare as wall-clock lines.
  subroutine teq_print_times() bind(c, name="teq_print_times")
    character(len=MAX_NAME_LEN) :: message
    acc_st = CPUTime() - acc_st
    write(message,'(a,F8.2)') ' Assembly: (s)', acc_at
    call Info('TransportEquationSolver', message)
    write(message,'(a,F8.2)') ' Solve:    (s)', acc_st
    call Info('TransportEquationSolver', message)
  end subroutine teq_print_times

  subroutine teq_restore_dt(dt) bind(c, name="teq_restore_dt")
    real(c_double), value :: dt
    acc_solver % dt = dt
  end subroutine teq_restore_dt

  ! -- outflow-BC retrieval + advection-variable update (lines 450-473) -----
  subroutine teq_post() bind(c, name="teq_post")
    integer :: i, j, k, t
    character(len=MAX_NAME_LEN) :: variablename
    real(kind=dp), pointer :: u(:)
    integer, pointer :: uperm(:)
    u => acc_solver % variable % values
    uperm => acc_solver % variable % perm
    select case (acc_dofs)
    case (1)
      call teq_modbc(acc_advvar, 1, u, uperm)
    case default
      do i = 1, acc_dofs
        write(variablename,'(A,A,I1)') trim(acc_advvar),' ', i
        call teq_modbc(variablename, i, u, uperm)
      end do
    end select
    do i = 1, acc_model % numberofnodes
      j = uperm(i)
      k = acc_u0v(1) % var % perm(i)
      if (j > 0 .and. k > 0) then
        do t = 1, acc_dofs
          acc_u0v(t) % var % values(k) = u((j-1)*acc_dofs+t)
        end do
      else
        call Fatal('TransportEquationSolver', &
            'Nonmatching variable permutations')
      end if
    end do
  end subroutine teq_post

  ! Exact transcription of ModifyBoundaryValues (lines 618-658).
  subroutine teq_modbc(name, component, u, perm)
    character(len=*) :: name
    integer :: component, perm(:)
    real(kind=dp) :: u(:)
    type(element_t), pointer :: currentelement
    integer, pointer :: nodeindexes(:)
    integer :: i, j, k, n, t
    logical :: gotit
    real(kind=dp) :: values(acc_model % maxelementnodes)
    do t = acc_model % numberofbulkelements + 1, &
        acc_model % numberofbulkelements + acc_model % numberofboundaryelements
      currentelement => acc_model % elements(t)
      n = currentelement % type % numberofnodes
      nodeindexes => currentelement % nodeindexes
      do i = 1, acc_model % numberofbcs
        if (currentelement % boundaryinfo % constraint == &
            acc_model % bcs(i) % tag) then
          values(1:n) = ListGetReal(acc_model % bcs(i) % values, &
              name, n, nodeindexes, gotit)
          if (gotit) then
            do j = 1, n
              k = perm(nodeindexes(j))
              if (k > 0) then
                k = acc_dofs * (k-1) + component
                u(k) = values(j)
              end if
            end do
          end if
        end if
      end do
    end do
  end subroutine teq_modbc

  subroutine teq_fatal_alloc() bind(c, name="teq_fatal_alloc")

    call Fatal('TransportEquationSolver', 'Memory allocation error, aborting.')
  end subroutine teq_fatal_alloc

end module transportequation_acc
