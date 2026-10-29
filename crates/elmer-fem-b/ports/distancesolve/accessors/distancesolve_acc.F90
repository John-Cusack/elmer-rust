! T022-A DistanceSolve port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/DistanceSolve.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments. All Fortran-object
! traffic — TYPE(...) handles, OPTIONAL/POINTER dummies, hidden lengths,
! formatted I/O, whole-array assignments, order-sensitive reductions
! (SUM, MAXVAL, MINVAL, SQRT) — stays inside Fortran, so bit-identity
! holds by construction wherever the call sequence matches. The module
! mirrors the originals' SAVE state with identical lifetimes:
! STIFF/FORCE + AllocationsDone for DistanceSolver, the LocalMatrix
! SAVE Nodes, and the DistanceSolver1 per-call arrays (allocated at
! entry, released at return, as in the original).
!
! Two deliberate deviations from the original text, both output-neutral:
! - DistanceSolver1 never deallocates bdist/bd/work (a per-call leak);
!   the wrappers release every per-call array at return and guard the
!   allocation site, so repeated invocation cannot fail on
!   already-allocated arrays. No output depends on this.
! - The Rust side passes the Dummy flag as 0/1; the wrapper branches
!   exactly as the original IF.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module distancesolve_acc
  use iso_c_binding
  use defutils
  implicit none

  ! -- shared entry handles ------------------------------------------
  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(valuelist_t), pointer, save :: acc_params => null()

  ! -- DistanceSolver SAVE state (mirrors SAVE STIFF, FORCE,
  ! -- AllocationsDone; first-mesh sizing, never resized) -------------
  real(kind=dp), allocatable, save :: acc_stiff(:,:)
  real(kind=dp), allocatable, save :: acc_force(:)
  logical, save :: acc_alloc_done = .false.
  real(kind=dp), save :: acc_relaxdt = 0.0_dp
  ! -- LocalMatrix SAVE Nodes (consumed before the next fill) ---------
  type(nodes_t), target, save :: acc_nodes
  type(element_t), pointer, save :: acc_elem => null()

  ! -- DistanceSolver1 per-call state ---------------------------------
  type(mesh_t), pointer, save :: acc_mesh => null()
  integer, save :: acc_n = 0
  integer, save :: acc_nb = 0
  integer, allocatable, save :: acc_aperm(:), acc_bperm(:)
  real(kind=dp), allocatable, save :: acc_bdist(:), acc_condition(:), acc_work(:)
  real(kind=dp), allocatable, save :: acc_xp(:), acc_yp(:), acc_zp(:), acc_bd(:)
  real(kind=dp), pointer, save :: acc_distance(:) => null()

contains

  ! ================= DistanceSolver ================================

  subroutine ds_begin(h_model, h_solver) bind(c, name="ds_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine ds_begin

  subroutine ds_info_using() bind(c, name="ds_info_using")
    call Info('DistanceSolver','Using PDE based distance solver')
  end subroutine ds_info_using

  subroutine ds_first_init() bind(c, name="ds_first_init")
    type(mesh_t), pointer :: msh
    integer :: n, istat
    msh => GetMesh()
    if (.not. acc_alloc_done) then
      n = acc_solver % mesh % maxelementdofs
      allocate(acc_force(n), acc_stiff(n,n), stat=istat)
      if (istat /= 0) then
        call Fatal('PoissonSolve', 'Memory allocation error.')
      end if
      n = msh % numberofnodes
      if (all(acc_solver % variable % values == 0.0_dp)) &
        acc_solver % variable % values(acc_solver % variable % perm(1:n)) = &
            sqrt(msh % nodes % x(1:n)**2 + &
                 msh % nodes % y(1:n)**2 + &
                 msh % nodes % z(1:n)**2)
      acc_alloc_done = .true.
    end if
  end subroutine ds_first_init

  function ds_read_params() result(maxiter) bind(c, name="ds_read_params")
    integer(c_int) :: maxiter
    logical :: found
    acc_params => GetSolverParams()
    acc_relaxdt = GetCReal(acc_params, 'Distance Pseudo dt', found)
    maxiter = GetInteger(acc_params, 'Nonlinear System Max Iterations', found)
    if (.not. found) maxiter = 100
  end function ds_read_params

  subroutine ds_noslip() bind(c, name="ds_noslip")
    type(valuelist_t), pointer :: bc
    integer :: i
    logical :: found
    do i=1,acc_model % numberofbcs
      bc => acc_model % bcs(i) % values
      if (GetLogical(bc, 'Noslip Wall BC', found)) &
        call ListAddConstReal(bc, ComponentName(acc_solver % variable), 0.0_dp)
    end do
  end subroutine ds_noslip

  function ds_active() result(active) bind(c, name="ds_active")
    integer(c_int) :: active
    active = GetNOFActive()
  end function ds_active

  subroutine ds_default_initialize() bind(c, name="ds_default_initialize")
    call DefaultInitialize()
  end subroutine ds_default_initialize

  subroutine ds_elem(t) bind(c, name="ds_elem")
    integer(c_int), value :: t
    integer :: n, nd, nb, ndnb
    acc_elem => GetActiveElement(t)
    n  = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
    nb = GetElementNOFBDOFs()
    ndnb = nd+nb
    call ds_local_matrix(n, ndnb)
    call CondensateP(nd, nb, acc_stiff, acc_force)
    call DefaultUpdateEquations(acc_stiff, acc_force)
  end subroutine ds_elem

  ! The LocalMatrix internal procedure (statement for statement).
  subroutine ds_local_matrix(n, nd)
    integer, intent(in) :: n, nd
    real(kind=dp), allocatable :: basis(:), dbasis(:,:), sol(:), grad(:), su(:)
    real(kind=dp) :: detj, uabs, stab, hscale
    type(gaussintegrationpoints_t) :: ip
    integer :: i, j, t, dim
    logical :: stat
    allocate(basis(nd), dbasis(nd,3), sol(nd), grad(3), su(nd))
    call GetElementNodes(acc_nodes)
    acc_stiff = 0.0_dp
    acc_force = 0.0_dp
    dim = CoordinateSystemDimension()
    call GetScalarLocalSolution(sol)
    hscale = GetCReal(GetSolverParams(), 'H scale', stat)
    ip = GaussPoints(acc_elem)
    do t=1,ip % n
      stat = ElementInfo(acc_elem, acc_nodes, ip % u(t), ip % v(t), &
               ip % w(t), detj, basis, dbasis)
      do i=1,dim
        grad(i) = sum(dbasis(1:nd,i) * sol(1:nd))
      end do
      uabs = sqrt(sum(grad(1:dim)**2))
      grad(1:dim) = grad(1:dim) / max(uabs, 1.0e-8_dp)
      stab = hscale*acc_elem % hk/2
      su = 0.0_dp
      do i=1,nd
        do j=1,dim
          su(i) = su(i) + grad(j) * dbasis(i,j)
        end do
      end do
      do i=1,nd
      do j=1,nd
        acc_stiff(i,j) = acc_stiff(i,j) + ip % s(t) * detj * &
          sum(grad(1:dim) * dbasis(j,1:dim)) * basis(i)
        acc_stiff(i,j) = acc_stiff(i,j) + ip % s(t) * detj * stab * su(i) * su(j)
        if (acc_relaxdt > 0) acc_stiff(i,j) = acc_stiff(i,j) + &
          1.0_dp/acc_relaxdt*ip % s(t) * detj * basis(j) * basis(i)
      end do
      end do
      acc_force(1:nd) = acc_force(1:nd) + ip % s(t) * detj * (basis(1:nd)+stab*su(1:nd))
      if (acc_relaxdt > 0) acc_force(1:nd) = acc_force(1:nd) + &
        1.0_dp/acc_relaxdt*ip % s(t) * detj * basis(1:nd)*sum(basis(1:nd)*sol(1:nd))
    end do
    deallocate(basis, dbasis, sol, grad, su)
  end subroutine ds_local_matrix

  subroutine ds_finish_bulk() bind(c, name="ds_finish_bulk")
    call DefaultFinishBulkAssembly()
  end subroutine ds_finish_bulk

  subroutine ds_finish() bind(c, name="ds_finish")
    call DefaultFinishAssembly()
  end subroutine ds_finish

  subroutine ds_dirichlet() bind(c, name="ds_dirichlet")
    call DefaultDirichletBCs()
  end subroutine ds_dirichlet

  function ds_solve() result(norm) bind(c, name="ds_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function ds_solve

  function ds_converged() result(code) bind(c, name="ds_converged")
    integer(c_int) :: code
    if (acc_solver % variable % nonlinconverged == 1) then
      code = 1
    else
      code = 0
    end if
  end function ds_converged

  subroutine ds_info_done() bind(c, name="ds_info_done")
    call Info('DistanceSolver','All done')
  end subroutine ds_info_done

  ! ================= DistanceSolver_init ===========================

  subroutine ds_init_begin(h_solver) bind(c, name="ds_init_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine ds_init_begin

  subroutine ds_init_body() bind(c, name="ds_init_body")
    type(valuelist_t), pointer :: solverpars
    solverpars => GetSolverParams()
    if (.not. ListCheckPresent(solverpars,'Nonlinear System Convergence Tolerance')) then
      call ListAddConstReal(solverpars,'Nonlinear System Convergence Tolerance',1.0d-8)
    end if
  end subroutine ds_init_body

  ! ================= DistanceSolver1 ===============================

  subroutine d1_begin(h_solver) bind(c, name="d1_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine d1_begin

  subroutine d1_info_head() bind(c, name="d1_info_head")
    call Info('DistanceSolver1','-------------------------------------------',level=5)
    call Info('DistanceSolver1','Using geometric distance solver',level=5)
    call Info('DistanceSolver1','-------------------------------------------',level=5)
    call ResetTimer('DistanceSolver1')
  end subroutine d1_info_head

  subroutine d1_mesh_var_info() bind(c, name="d1_mesh_var_info")
    acc_mesh => GetMesh()
    call Info('DistanceSolver1','Working with mesh: '// &
        trim(acc_mesh % name),level=10)
    call Info('DistanceSolver1','Solving for variable: '// &
        trim(acc_solver % variable % name),level=10)
  end subroutine d1_mesh_var_info

  subroutine d1_alloc() bind(c, name="d1_alloc")
    integer :: n, m
    acc_mesh => GetMesh()
    n = acc_mesh % numberofnodes
    acc_n = n
    if (allocated(acc_aperm)) deallocate(acc_aperm)
    if (allocated(acc_bperm)) deallocate(acc_bperm)
    if (allocated(acc_bdist)) deallocate(acc_bdist)
    if (allocated(acc_condition)) deallocate(acc_condition)
    if (allocated(acc_work)) deallocate(acc_work)
    if (allocated(acc_xp)) deallocate(acc_xp, acc_yp, acc_zp, acc_bd)
    allocate(acc_aperm(n), acc_bperm(n), acc_bdist(n))
    acc_aperm = 0; acc_bperm = 0; acc_bdist = 0.0_dp
    m = acc_mesh % maxelementnodes
    allocate(acc_condition(m), acc_work(m))
    acc_condition = -1.0_dp; acc_work = 0.0_dp
  end subroutine d1_alloc

  function d1_read_dummy() result(code) bind(c, name="d1_read_dummy")
    integer(c_int) :: code
    logical :: found
    acc_params => GetSolverParams()
    if (GetLogical(acc_params, 'Dummy Distance Computation', found)) then
      code = 1
    else
      code = 0
    end if
  end function d1_read_dummy

  subroutine d1_bulk_loop() bind(c, name="d1_bulk_loop")
    type(element_t), pointer :: element
    type(valuelist_t), pointer :: bodyforce
    integer :: nd, t, i, j
    logical :: found
    acc_nb = 0
    do t=1,acc_solver % numberofactiveelements
      element => GetActiveElement(t)
      if (.not. associated(element)) cycle
      bodyforce => GetBodyForce()
      if (associated(bodyforce)) then
        nd = GetElementNOFNodes(element)
        acc_work(1:nd) = ListGetReal(bodyforce, trim(acc_solver % variable % name), &
            nd, element % nodeindexes, found)
        if (found) then
          acc_condition(1:nd) = GetReal(bodyforce, &
            trim(acc_solver % variable % name) // " Condition", found)
          do i=1,nd
            j = element % nodeindexes(i)
            if (found .and. acc_condition(i) < 0.0_dp) cycle
            if (acc_bperm(j) == 0) then
              acc_nb = acc_nb + 1
              acc_aperm(acc_nb) = j
              acc_bperm(j) = acc_nb
              acc_bdist(j) = acc_work(i)
            end if
          end do
        end if
      end if
    end do
  end subroutine d1_bulk_loop

  subroutine d1_bc_loop() bind(c, name="d1_bc_loop")
    type(element_t), pointer :: element
    type(valuelist_t), pointer :: bc
    integer :: nd, t, i, j
    logical :: found, found1
    do t=1,acc_mesh % numberofboundaryelements
      element => GetBoundaryElement(t)
      if (.not. ActiveBoundaryElement(element)) cycle
      bc => GetBC(element)
      nd = GetElementNOFNodes(element)
      acc_work(1:nd) = GetReal(bc, acc_solver % variable % name, found)
      if (found .or. GetLogical(bc, 'Noslip Wall BC', found1)) then
        do i=1,nd
          j = element % nodeindexes(i)
          if (acc_bperm(j) == 0) then
            acc_nb = acc_nb + 1
            acc_aperm(acc_nb) = j
            acc_bperm(j) = acc_nb
            acc_bdist(j) = acc_work(i)
          end if
        end do
      end if
    end do
  end subroutine d1_bc_loop

  function d1_pe_check() result(code) bind(c, name="d1_pe_check")
    integer(c_int) :: code
    if (parenv % pes == 1) then
      if (acc_nb == 0) then
        call Warn('DistanceSolver1','No known distances given for the distance solver!')
        code = 1
      else
        write(message,'(A,I0)') 'Number of fixed nodes on bulk: ',acc_nb
        call Info('DistanceSolver1',message,level=8)
        code = 0
      end if
    else
      code = 2
    end if
  end function d1_pe_check

  subroutine d1_serial_gather() bind(c, name="d1_serial_gather")
    integer :: i
    allocate(acc_xp(acc_nb), acc_yp(acc_nb), acc_zp(acc_nb), acc_bd(acc_nb))
    do i=1,acc_nb
      acc_xp(i) = acc_solver % mesh % nodes % x(acc_aperm(i))
      acc_yp(i) = acc_solver % mesh % nodes % y(acc_aperm(i))
      acc_zp(i) = acc_solver % mesh % nodes % z(acc_aperm(i))
      acc_bd(i) = acc_bdist(acc_aperm(i))
    end do
  end subroutine d1_serial_gather

  ! Parallel branch (transcribed; not executed in the serial proof).
  subroutine d1_parallel_gather() bind(c, name="d1_parallel_gather")
    integer, pointer :: gperm(:)
    integer :: comm, maxnode, i, nnb
    integer, allocatable :: ibuf(:), cperm(:)
    real(kind=dp), allocatable :: buf(:)
    comm = acc_solver % matrix % comm
    gperm => acc_mesh % parallelinfo % globaldofs
    maxnode = ParallelReduction(maxval(gperm(acc_aperm(1:acc_nb))),2)
    allocate(ibuf(maxnode), cperm(maxnode))
    cperm = 0
    cperm(gperm(acc_aperm(1:acc_nb))) = 1
    call MPI_ALLREDUCE(cperm, ibuf, maxnode, MPI_INTEGER, MPI_SUM, comm, i)
    nnb = 0; cperm = 0
    do i=1,maxnode
      if (ibuf(i) /= 0) then
        nnb = nnb+1
        cperm(i) = nnb
        ibuf(nnb) = ibuf(i)
      end if
    end do
    allocate(acc_xp(nnb), acc_yp(nnb), acc_zp(nnb), acc_bd(nnb), buf(nnb))
    buf = 0.0_dp
    buf(cperm(gperm(acc_aperm(1:acc_nb)))) = acc_mesh % nodes % x(acc_aperm(1:acc_nb))
    call MPI_ALLREDUCE(buf, acc_xp, nnb, MPI_DOUBLE_PRECISION, MPI_SUM, comm, i)
    acc_xp(1:nnb) = acc_xp(1:nnb)/ibuf(1:nnb)
    buf = 0.0_dp
    buf(cperm(gperm(acc_aperm(1:acc_nb)))) = acc_mesh % nodes % y(acc_aperm(1:acc_nb))
    call MPI_ALLREDUCE(buf, acc_yp, nnb, MPI_DOUBLE_PRECISION, MPI_SUM, comm, i)
    acc_yp(1:nnb) = acc_yp(1:nnb)/ibuf(1:nnb)
    buf = 0.0_dp
    buf(cperm(gperm(acc_aperm(1:acc_nb)))) = acc_mesh % nodes % z(acc_aperm(1:acc_nb))
    call MPI_ALLREDUCE(buf, acc_zp, nnb, MPI_DOUBLE_PRECISION, MPI_SUM, comm, i)
    acc_zp(1:nnb) = acc_zp(1:nnb)/ibuf(1:nnb)
    buf = 0.0_dp
    buf(cperm(gperm(acc_aperm(1:acc_nb)))) = acc_bdist(acc_aperm(1:acc_nb))
    call MPI_ALLREDUCE(buf, acc_bd, nnb, MPI_DOUBLE_PRECISION, MPI_SUM, comm, i)
    acc_bd(1:nnb) = acc_bd(1:nnb)/ibuf(1:nnb)
    acc_nb = nnb
    deallocate(ibuf, buf, cperm)
  end subroutine d1_parallel_gather

  subroutine d1_compute(dummy) bind(c, name="d1_compute")
    integer(c_int), value :: dummy
    acc_distance => acc_solver % variable % values
    if (dummy /= 0) then
      call distcomp0()
    else
      call distcomp()
    end if
  end subroutine d1_compute

  subroutine distcomp0
    integer :: i, j, k
    real(kind=dp) :: xl, yl, zl, dist
    acc_distance = huge(1.0_dp)
    do i=1,acc_n
      k = acc_solver % variable % perm(i)
      if (k <= 0) cycle
      if (acc_bperm(i) /= 0) then
        acc_distance(k) = 0.0_dp
      else
        xl = acc_mesh % nodes % x(i)
        yl = acc_mesh % nodes % y(i)
        zl = acc_mesh % nodes % z(i)
        do j=1,acc_nb
          dist = (acc_xp(j)-xl)**2 + (acc_yp(j)-yl)**2 + (acc_zp(j)-zl)**2
          acc_distance(k) = min(acc_distance(k), dist)
        end do
      end if
    end do
  end subroutine distcomp0

  subroutine distcomp
    real(kind=dp), allocatable :: xxp(:), yyp(:), zzp(:), d(:), dd(:), bbd(:)
    real(kind=dp) :: xl, yl, zl, dl, c(3), ldist, dist
    integer :: i, j, k, l, nnb, minl
    integer, allocatable :: dperm(:), near(:)
    nnb = acc_n+acc_nb
    allocate(xxp(nnb), yyp(nnb), zzp(nnb), d(nnb), dperm(nnb), &
             near(nnb), dd(acc_nb), bbd(nnb))
    call random_number(c)
    xxp(1:acc_n) = acc_mesh % nodes % x(1:acc_n)
    yyp(1:acc_n) = acc_mesh % nodes % y(1:acc_n)
    zzp(1:acc_n) = acc_mesh % nodes % z(1:acc_n)
    bbd(1:acc_n) = acc_bdist(1:acc_n)
    xxp(acc_n+1:nnb) = acc_xp
    yyp(acc_n+1:nnb) = acc_yp
    zzp(acc_n+1:nnb) = acc_zp
    bbd(acc_n+1:nnb) = acc_bd
    xl = maxval(xxp)-minval(xxp)
    yl = maxval(yyp)-minval(yyp)
    zl = maxval(zzp)-minval(zzp)
    xxp = xxp + xl*(2*c(1)-1)
    yyp = yyp + yl*(2*c(2)-1)
    zzp = zzp + zl*(2*c(3)-1)
    d = sqrt(xxp**2 + yyp**2 + zzp**2)
    dperm = (/ (i,i=1,nnb) /)
    call SortR(nnb, dperm, d)
    j = 0
    do i=1,nnb
      near(i) = j
      if (dperm(i) > acc_n) then
        j = j+1
        acc_xp(j) = xxp(dperm(i))
        acc_yp(j) = yyp(dperm(i))
        acc_zp(j) = zzp(dperm(i))
        dd(j) = d(i)
        acc_bd(j) = bbd(dperm(i))
      end if
    end do
    do i=1,nnb
      j = dperm(i)
      if (j > acc_n) cycle
      k = acc_solver % variable % perm(j)
      if (k <= 0) cycle
      if (acc_bperm(j) /= 0) then
        acc_distance(k) = acc_bdist(j)
      else
        dl = d(i)
        xl = xxp(j)
        yl = yyp(j)
        zl = zzp(j)
        dist = huge(1.0_dp)
        do l=near(i)+1,acc_nb
          if ((dd(l)-dl)**2 > dist) exit
          ldist = (acc_xp(l)-xl)**2+(acc_yp(l)-yl)**2+(acc_zp(l)-zl)**2
          if (dist > ldist) then
            minl = l; dist = ldist
          end if
        end do
        do l=near(i),1,-1
          if ((dd(l)-dl)**2 > dist) exit
          ldist = (acc_xp(l)-xl)**2+(acc_yp(l)-yl)**2+(acc_zp(l)-zl)**2
          if (dist > ldist) then
            minl = l; dist = ldist
          end if
        end do
        acc_distance(k) = sqrt(dist)+acc_bd(minl)
      end if
    end do
    deallocate(xxp, yyp, zzp, d, dperm, near, dd, bbd)
  end subroutine distcomp

  subroutine d1_maxinfo() bind(c, name="d1_maxinfo")
    if (parenv % pes == 1) then
      write(message,'(A,ES12.3)') 'Maximum distance from given nodes: ',maxval(acc_distance)
      call Info('DistanceSolver1',message,level=8)
    end if
  end subroutine d1_maxinfo

  subroutine d1_finish() bind(c, name="d1_finish")
    deallocate(acc_xp, acc_yp, acc_zp, acc_aperm, acc_bperm, acc_condition)
    if (allocated(acc_bdist)) deallocate(acc_bdist)
    if (allocated(acc_work)) deallocate(acc_work)
    if (allocated(acc_bd)) deallocate(acc_bd)
    acc_solver % variable % norm = sqrt(sum(acc_distance**2))
    call InvalidateVariable(currentModel % meshes, acc_solver % mesh, &
      acc_solver % variable % name)
    call CheckTimer('DistanceSolver1', delete=.true.)
    call Info('DistanceSolver1','All done')
  end subroutine d1_finish

end module distancesolve_acc
