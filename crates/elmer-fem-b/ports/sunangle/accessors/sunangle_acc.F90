! T021-B SunAngle port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/SunAngle.F90, compiled against the built .mod
! files of the serial configure (ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) / INTEGER(C_INT)
! arrays). All Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER
! dummies, assumed-shape descriptors, hidden lengths, array-pointer results,
! the whole-array fills and the MPI calls — stays inside Fortran, so the
! Rust side is pure C calls. The wrapper module mirrors SunAngle's own
! host-associated handles (Solver, Model, Mesh, SunVar, Params, SunAngle,
! SunPerm, x/y/z) with identical lifetimes; the RidgeData SAVE state lives
! on the Rust side (same allocate/release/visited lifetimes as the original).
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module sunangle_acc
  use iso_c_binding
  use defutils
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(model_t), pointer, save :: acc_model => null()
  type(mesh_t), pointer, save :: acc_mesh => null()
  type(variable_t), pointer, save :: acc_sunvar => null()
  type(valuelist_t), pointer, save :: acc_params => null()
  real(kind=dp), pointer, save :: acc_sun(:) => null()
  real(kind=dp), pointer, save :: acc_x(:) => null()
  real(kind=dp), pointer, save :: acc_y(:) => null()
  real(kind=dp), pointer, save :: acc_z(:) => null()
  integer, pointer, save :: acc_perm(:) => null()

contains

  ! -- entry state ------------------------------------------------------
  subroutine sua_begin(h_model, h_solver) bind(c, name="sua_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    acc_mesh => GetMesh()
    acc_params => GetSolverParams()
    acc_sunvar => acc_solver % variable
    acc_sun => acc_sunvar % values
    acc_perm => acc_sunvar % perm
    acc_x => acc_mesh % nodes % x
    acc_y => acc_mesh % nodes % y
    acc_z => acc_mesh % nodes % z
  end subroutine sua_begin

  ! -- init entry -------------------------------------------------------
  subroutine sua_init(h_solver) bind(c, name="sua_init")
    type(c_ptr), value :: h_solver
    type(solver_t), pointer :: solver
    call c_f_pointer(h_solver, solver)
    call ListAddNewString(solver % values, 'Variable', 'CriticalTan')
    call ListAddLogical(solver % values, 'No Matrix', .TRUE.)
  end subroutine sua_init

  ! -- Info transcripts (exact Caller/string/Level of the original) -----
  subroutine sua_info(which) bind(c, name="sua_info")
    integer(c_int), value :: which
    if (which == 1) then
      call Info('SunAngleSolver', '-------------------------------------', Level=4)
    else if (which == 2) then
      call Info('SunAngleSolver', 'Determining the critical sun angle', Level=4)
    else if (which == 3) then
      call Info('SunAngleSolver', '-------------------------------------', Level=4)
    else if (which == 4) then
      call Info('SunAngleSolver', 'Computing the sun angles for real!', Level=7)
    else if (which == 5) then
      call Info('SunAngleSolver', 'All done', Level=10)
    else if (which == 6) then
      call Info('SunAngleSolver', 'Creating mesh edges', Level=10)
    else if (which == 7) then
      call Info('SunAngleSolver', 'Communicating ridges from different partitions', Level=6)
    end if
  end subroutine sua_info

  function sua_info_active20() result(r) bind(c, name="sua_info_active20")
    integer(c_int) :: r
    r = merge(1, 0, InfoActive(20))
  end function sua_info_active20

  subroutine sua_print_totals(mype, k, l, nmax) bind(c, name="sua_print_totals")
    integer(c_int), value :: mype, k, l, nmax
    print *, 'Total number of edges to send and recieve: ', mype, k, l, nmax
  end subroutine sua_print_totals

  subroutine sua_fatal_release() bind(c, name="sua_fatal_release")
    call Fatal('ReleaseRidge', 'Cannot release an unallocated RidgeData!')
  end subroutine sua_fatal_release

  ! -- mesh / variable geometry -----------------------------------------
  function sua_dim() result(n) bind(c, name="sua_dim")
    integer(c_int) :: n
    n = acc_mesh % meshdim
  end function sua_dim

  function sua_nnodes() result(n) bind(c, name="sua_nnodes")
    integer(c_int) :: n
    n = acc_mesh % numberofnodes
  end function sua_nnodes

  function sua_nedges() result(n) bind(c, name="sua_nedges")
    integer(c_int) :: n
    n = acc_mesh % numberofedges
  end function sua_nedges

  function sua_dofs() result(n) bind(c, name="sua_dofs")
    integer(c_int) :: n
    n = acc_sunvar % dofs
  end function sua_dofs

  function sua_nvals() result(n) bind(c, name="sua_nvals")
    integer(c_int) :: n
    n = size(acc_sun)
  end function sua_nvals

  function sua_nperm() result(n) bind(c, name="sua_nperm")
    integer(c_int) :: n
    n = size(acc_perm)
  end function sua_nperm

  function sua_edges_associated() result(r) bind(c, name="sua_edges_associated")
    integer(c_int) :: r
    r = merge(1, 0, associated(acc_mesh % edges))
  end function sua_edges_associated

  subroutine sua_ensure_edges(dim) bind(c, name="sua_ensure_edges")
    integer(c_int), value :: dim
    if (dim == 2) then
      call FindMeshEdges2D(acc_mesh)
    else
      call FindMeshEdges3D(acc_mesh)
    end if
  end subroutine sua_ensure_edges

  subroutine sua_edge_nodes(t, i0, i1) bind(c, name="sua_edge_nodes")
    integer(c_int), value :: t
    integer(c_int) :: i0, i1
    i0 = acc_mesh % edges(t) % nodeindexes(1)
    i1 = acc_mesh % edges(t) % nodeindexes(2)
  end subroutine sua_edge_nodes

  ! -- snapshots into Rust mirrors --------------------------------------
  subroutine sua_fill_x(buf, n) bind(c, name="sua_fill_x")
    integer(c_int), value :: n
    real(c_double) :: buf(n)
    buf(1:n) = acc_x(1:n)
  end subroutine sua_fill_x

  subroutine sua_fill_y(buf, n) bind(c, name="sua_fill_y")
    integer(c_int), value :: n
    real(c_double) :: buf(n)
    buf(1:n) = acc_y(1:n)
  end subroutine sua_fill_y

  subroutine sua_fill_z(buf, n) bind(c, name="sua_fill_z")
    integer(c_int), value :: n
    real(c_double) :: buf(n)
    buf(1:n) = acc_z(1:n)
  end subroutine sua_fill_z

  subroutine sua_fill_perm(buf, n) bind(c, name="sua_fill_perm")
    integer(c_int), value :: n
    integer(c_int) :: buf(n)
    buf(1:n) = acc_perm(1:n)
  end subroutine sua_fill_perm

  subroutine sua_fill_sun(buf, n) bind(c, name="sua_fill_sun")
    integer(c_int), value :: n
    real(c_double) :: buf(n)
    buf(1:n) = acc_sun(1:n)
  end subroutine sua_fill_sun

  subroutine sua_sun_fill(v) bind(c, name="sua_sun_fill")
    real(c_double), value :: v
    acc_sun = v
  end subroutine sua_sun_fill

  subroutine sua_sun_store(buf, n) bind(c, name="sua_sun_store")
    integer(c_int), value :: n
    real(c_double) :: buf(n)
    acc_sun(1:n) = buf(1:n)
  end subroutine sua_sun_store

  subroutine sua_where_fix() bind(c, name="sua_where_fix")
    where (acc_sun < -99.0)
      acc_sun = 0.0_dp
    end where
  end subroutine sua_where_fix

  ! -- ValueList queries -------------------------------------------------
  function sua_test_height_present() result(r) bind(c, name="sua_test_height_present")
    integer(c_int) :: r
    r = merge(1, 0, ListCheckPresentAnyMaterial(acc_model, 'Test Height'))
  end function sua_test_height_present

  ! Whole SetTestHeight internal procedure, verbatim: element loop with
  ! CurrentElement setup, per-material GetReal broadcast and z update.
  subroutine sua_set_test_height() bind(c, name="sua_set_test_height")
    type(element_t), pointer :: element
    type(valuelist_t), pointer :: material
    real(kind=dp) :: elemheight(4)
    integer :: efirst, elast, t, n
    logical :: found
    if (acc_mesh % elements(1) % type % elementcode > 500) then
      efirst = acc_mesh % numberofbulkelements + 1
      elast = acc_mesh % numberofbulkelements + acc_mesh % numberofboundaryelements
    else
      efirst = 1
      elast = acc_mesh % numberofbulkelements
    end if
    do t = efirst, elast
      element => acc_mesh % elements(t)
      acc_model % currentelement => element
      material => GetMaterial(element)
      n = element % type % numberofnodes
      elemheight(1:n) = GetReal(material, 'Test Height', found)
      if (any(acc_perm(element % nodeindexes) == 0)) cycle
      acc_mesh % nodes % z(element % nodeindexes) = elemheight(1:n)
    end do
  end subroutine sua_set_test_height

  subroutine sua_enforce(which, val, found) bind(c, name="sua_enforce")
    integer(c_int), value :: which
    integer(c_int) :: val, found
    logical :: l, f
    if (which == 1) then
      l = ListGetLogical(acc_params, 'Enforce Parallel', f)
    else
      l = ListGetLogical(acc_params, 'Enforce Serial', f)
    end if
    val = merge(1, 0, l)
    found = merge(1, 0, f)
  end subroutine sua_enforce

  function sua_pes() result(n) bind(c, name="sua_pes")
    integer(c_int) :: n
    n = ParEnv % pes
  end function sua_pes

  function sua_mype() result(n) bind(c, name="sua_mype")
    integer(c_int) :: n
    n = ParEnv % mype
  end function sua_mype

  subroutine sua_maxdist(val, got) bind(c, name="sua_maxdist")
    real(c_double) :: val
    integer(c_int) :: got
    logical :: f
    val = ListGetCReal(acc_params, 'SunAngle search distance', f)
    got = merge(1, 0, f)
  end subroutine sua_maxdist

  function sua_test_angle() result(v) bind(c, name="sua_test_angle")
    real(c_double) :: v
    v = ListGetCReal(acc_params, 'Test Angle')
  end function sua_test_angle

  ! -- MPI ridge exchange (serial tree links libmpi_stubs; these run only
  ! -- when PEs > 1, exactly like the original's guarded calls) -----------
  subroutine sua_mpi_allreduce_bb(buf, pes) bind(c, name="sua_mpi_allreduce_bb")
    integer(c_int), value :: pes
    real(c_double) :: buf(4*pes)
    integer :: comm, ierr
    comm = ELMER_COMM_WORLD
    call MPI_ALLREDUCE(MPI_IN_PLACE, buf, 4*pes, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)
  end subroutine sua_mpi_allreduce_bb

  subroutine sua_mpi_barrier() bind(c, name="sua_mpi_barrier")
    integer :: comm, ierr
    comm = ELMER_COMM_WORLD
    call MPI_BARRIER(comm, ierr)
  end subroutine sua_mpi_barrier

  subroutine sua_mpi_send_count(n, proc) bind(c, name="sua_mpi_send_count")
    integer(c_int), value :: n, proc
    integer :: comm, ierr
    comm = ELMER_COMM_WORLD
    call MPI_SEND(n, 1, MPI_INTEGER, proc, 999, comm, ierr)
  end subroutine sua_mpi_send_count

  subroutine sua_mpi_recv_count(n, proc) bind(c, name="sua_mpi_recv_count")
    integer(c_int) :: n
    integer(c_int), value :: proc
    integer :: comm, ierr, status(MPI_STATUS_SIZE)
    comm = ELMER_COMM_WORLD
    call MPI_RECV(n, 1, MPI_INTEGER, proc, 999, comm, status, ierr)
  end subroutine sua_mpi_recv_count

  subroutine sua_mpi_bsend(buf, n, proc) bind(c, name="sua_mpi_bsend")
    integer(c_int), value :: n, proc
    real(c_double) :: buf(n)
    integer :: comm, ierr
    comm = ELMER_COMM_WORLD
    call MPI_BSEND(buf, n, MPI_DOUBLE_PRECISION, proc, 1001, comm, ierr)
  end subroutine sua_mpi_bsend

  subroutine sua_mpi_recv(buf, n, proc) bind(c, name="sua_mpi_recv")
    integer(c_int), value :: n, proc
    real(c_double) :: buf(n)
    integer :: comm, ierr, status(MPI_STATUS_SIZE)
    comm = ELMER_COMM_WORLD
    call MPI_RECV(buf, n, MPI_DOUBLE_PRECISION, proc, 1001, comm, status, ierr)
  end subroutine sua_mpi_recv

end module sunangle_acc
