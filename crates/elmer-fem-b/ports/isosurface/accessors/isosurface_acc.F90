! T022-A Isosurface port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/Isosurface.F90, compiled against the built .mod
! files of the target configure (A2 shipped mechanism, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape C arrays). All Fortran-object
! traffic — TYPE(...) handles, OPTIONAL/POINTER dummies, assumed-shape
! arrays, hidden lengths, pointer-function results — stays inside Fortran,
! so the Rust side is pure C calls. The wrapper module mirrors
! IsosurfaceSolver's own host state (current mesh handles, level-field
! bindings, per-call temporaries, SAVE state) with identical lifetimes.
!
! Granularity (documented, not incidental): the entry body's control flow
! (phases, branches, GOTO-100 exits, collective reporting) is transcribed
! statement-for-statement into Rust, one BIND(C) phase per statement block
! below. The twelve CONTAINS internals stay Fortran-side as whole module
! procedures because their dummies cannot cross the C boundary at all:
! TYPE(Element_t) by value (CreateSurfaceFromTetra, CreateLineFromTriangle),
! assumed-shape arrays (OrientIsoTriangle, CreateSurfaceFromTetra), and
! pointer-function results threaded through host association. Moving those
! across would need the callee-side descriptor ABI — the new-ABI refusal
! class — so they remain whole units operating on the same module state.
! All floating-point work (edge interpolation, magnitude SQRT sums, the
! orientation cross products) stays in Fortran with the same compiler, so
! bit-identity holds by construction wherever the call sequence matches.
!
! SAVE mapping (mirrors the original SAVE statement plus the implicit SAVE
! of `calls`): acc_fixed, acc_exist, acc_invperm, acc_interpolant,
! acc_isomesh, acc_savepartmp, acc_calls. Everything else is per-call
! automatic state, reset in iso_begin (nullify/deallocate guards) so that
! re-entry matches the original's fresh-automatic semantics. The original
! has no STAT= allocation checks, so there are none to replicate.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module isosurface_acc
  use iso_c_binding
  use CoordinateSystems
  use DefUtils
  use SaveUtils
  implicit none

  integer, parameter :: ACC_MAXINTS = 32

  ! -- entry handles ----------------------------------------------------
  type(solver_t), pointer :: acc_solver => null()
  type(model_t), pointer :: acc_model => null()
  type(mesh_t), pointer :: acc_mesh => null()
  type(mesh_t), pointer :: acc_origmesh => null()
  type(mesh_t), pointer :: acc_isomesh => null()
  type(mesh_t), pointer :: acc_tempmesh => null()
  type(valuelist_t), pointer :: acc_params => null()

  ! -- level-field bindings ----------------------------------------------
  type(variable_t), pointer :: acc_levelvar => null()
  real(kind=dp), pointer :: acc_levelfun(:) => null()
  integer, pointer :: acc_levelperm(:) => null()
  real(kind=dp), pointer :: acc_levelvalues(:,:) => null()
  integer :: acc_leveldofs = 0
  integer :: acc_levelcomps = 0
  integer :: acc_nolevels = 0
  real(kind=dp) :: acc_levelvalue = 0.0_dp

  ! -- per-call assembly state --------------------------------------------
  integer :: acc_partmp(6) = 0
  integer :: acc_parsizes(6) = 0
  integer :: acc_elemfirst = 1
  integer :: acc_elemlast = 0
  logical :: acc_gotomask = .false.
  logical :: acc_movingmesh = .false.
  logical :: acc_fixed = .false.
  logical :: acc_exist = .false.
  logical :: acc_isocreated = .false.
  integer :: acc_dim = 0
  integer :: acc_neweltype = 0
  integer :: acc_newelnodes = 0
  integer :: acc_noorigelems = 0
  integer :: acc_nonewelems = 0
  integer :: acc_noedges = 0
  integer :: acc_noisonodes = 0
  integer :: acc_nosurfaces = 0
  integer :: acc_nomeshnodesowned = 0
  integer :: acc_noedgesowned = 0
  integer :: acc_noisonodesowned = 0
  integer :: acc_calls = 0
  integer :: acc_savepartmp(6) = 0
  real(kind=dp), allocatable :: acc_elemfun(:)
  integer, allocatable :: acc_elemperm(:)
  integer, allocatable :: acc_eperm(:)
  integer, allocatable :: acc_isoedge(:)
  integer, allocatable :: acc_isolevel(:)
  integer, allocatable :: acc_edgeowner(:)
  integer, allocatable :: acc_edgenbrs(:,:)
  integer, allocatable :: acc_nedgenbrs(:)
  type(element_t), pointer :: acc_newelements(:) => null()

  ! -- interpolant tables --------------------------------------------------
  integer :: acc_ints = 0
  character(len=MAX_NAME_LEN) :: acc_vname(ACC_MAXINTS)
  type(variabletable_t) :: acc_vfulltab(ACC_MAXINTS)
  type(variabletable_t) :: acc_visotab(ACC_MAXINTS)
  integer, allocatable :: acc_invperm(:,:)
  real(kind=dp), allocatable :: acc_interpolant(:)

  ! -- gmsh-geo gather state -------------------------------------------------
  integer :: acc_ngeonodes = 0
  integer :: acc_ngeolines = 0
  real(kind=dp), allocatable :: acc_geox(:), acc_geoy(:), acc_geoz(:)
  integer, allocatable :: acc_geoline(:,:)

contains

  ! -- entry state --------------------------------------------------------
  subroutine iso_begin(h_model, h_solver) bind(c, name="iso_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    ! Fresh-automatic semantics for every non-SAVE host variable: the
    ! original re-enters with all of these unassociated/unallocated.
    acc_mesh => null()
    acc_origmesh => null()
    acc_tempmesh => null()
    acc_newelements => null()
    acc_params => null()
    acc_levelvar => null()
    acc_levelfun => null()
    acc_levelperm => null()
    acc_levelvalues => null()
    if (allocated(acc_elemfun)) deallocate(acc_elemfun)
    if (allocated(acc_elemperm)) deallocate(acc_elemperm)
    if (allocated(acc_eperm)) deallocate(acc_eperm)
    if (allocated(acc_isoedge)) deallocate(acc_isoedge, acc_isolevel)
    if (allocated(acc_edgeowner)) deallocate(acc_edgeowner, acc_edgenbrs, acc_nedgenbrs)
    if (allocated(acc_geox)) deallocate(acc_geox, acc_geoy, acc_geoz)
    if (allocated(acc_geoline)) deallocate(acc_geoline)
    acc_partmp = 0
    acc_isocreated = .false.
    acc_noedges = 0
    acc_noisonodes = 0
    acc_nosurfaces = 0
    acc_nomeshnodesowned = 0
    acc_noedgesowned = 0
    acc_noisonodesowned = 0
    acc_nonewelems = 0
    acc_newelements => null()
  end subroutine iso_begin

  subroutine iso_info_start() bind(c, name="iso_info_start")
    call Info('IsosurfaceSolver', '-------------------------------------', level=4)
    call Info('IsosurfaceSolver', 'Determining the isosurface', level=4)
    call Info('IsosurfaceSolver', '-------------------------------------', level=4)
  end subroutine iso_info_start

  ! Lines 121-169: handles, params, level variable lookup.
  ! Returns 1 when the body must branch to label 100.
  function iso_read_var() result(code) bind(c, name="iso_read_var")
    integer(c_int) :: code
    character(len=MAX_NAME_LEN) :: levelvariablename
    logical :: found
    code = 0
    acc_mesh => GetMesh()
    acc_origmesh => acc_mesh
    acc_noorigelems = acc_mesh % NumberOfBulkElements
    ! Every exit falls through to the collective size reduction, so the
    ! reported quantities start zeroed (iso_begin did that above).
    acc_params => GetSolverParams()
    acc_fixed = GetLogical(acc_params, 'Isosurface Fixed', found)
    acc_movingmesh = GetLogical(acc_params, 'Moving Mesh', found)
    levelvariablename = GetString(acc_params, 'Isosurface Variable')
    acc_levelvar => VariableGet(acc_mesh % Variables, levelvariablename)
    if (.not. associated(acc_levelvar)) then
      call Error('Isosurface', 'Missing isosurface variable: ' // trim(levelvariablename))
      code = 1
      return
    end if
    acc_levelfun => acc_levelvar % Values
    acc_levelperm => acc_levelvar % Perm
    acc_leveldofs = acc_levelvar % Dofs
    ! For a vector level function the isosurface is taken of its magnitude,
    ! over the spatial components.
    acc_levelcomps = min(acc_mesh % MeshDim, acc_leveldofs)
  end function iso_read_var

  ! Lines 174-187: level values.
  function iso_read_levels() result(code) bind(c, name="iso_read_levels")
    integer(c_int) :: code
    logical :: found
    code = 0
    acc_nolevels = 0
    acc_levelvalues => ListGetConstRealArray(acc_params, 'Isosurface values', found)
    if (found) then
      acc_nolevels = size(acc_levelvalues, 1)
      acc_levelvalue = acc_levelvalues(1, 1)
    else
      acc_levelvalue = ListGetCReal(acc_params, 'Isosurface value', found)
      if (found) acc_nolevels = 1
    end if
    if (.not. found) then
      call Warn('IsosurfaceSolver', 'Could not determine Isosurface value')
      code = 1
    end if
  end function iso_read_levels

  ! Lines 192-203: fixed-surface fast path.
  function iso_fixed_remap() result(code) bind(c, name="iso_fixed_remap")
    integer(c_int) :: code
    code = 0
    if (acc_fixed .and. acc_exist) then
      call Info('IsosurfaceSolver', 'Remapping the fields')
      call acc_remap()
      acc_partmp = acc_savepartmp
      code = 1
    end if
  end function iso_fixed_remap

  ! Lines 205-248: drop the previous surface, size the per-call element
  ! buffers, install the simplex maps, set the element range and mask.
  subroutine iso_prep() bind(c, name="iso_prep")
    logical :: found
    call acc_releaseisomesh()
    allocate(acc_elemfun(acc_mesh % MaxElementNodes))
    allocate(acc_elemperm(acc_mesh % MaxElementNodes))
    acc_elemfirst = 1
    acc_elemlast = acc_mesh % NumberOfBulkElements
    acc_gotomask = .false.
    if (ListGetLogicalAnyBC(acc_model, 'Create Isosurface')) then
      acc_gotomask = .true.
      acc_elemfirst = acc_mesh % NumberOfBulkElements + 1
      acc_elemlast = acc_mesh % NumberOfBulkElements + acc_mesh % NumberOfBoundaryElements
    else if (ListGetLogicalAnyBodyForce(acc_model, 'Create Isosurface')) then
      acc_gotomask = .true.
    end if
  end subroutine iso_prep

  ! Lines 250-298: leading-dimension detection with the cross-partition
  ! agreement reduction. Returns 1 (after the Warn) when no 2D/3D path.
  function iso_detect_dim() result(code) bind(c, name="iso_detect_dim")
    integer(c_int) :: code
    type(element_t), pointer :: element
    integer, pointer :: indexes(:)
    integer :: i, j, n, dim, ierr
    logical :: found
    code = 0
    dim = 0
    do i = acc_elemfirst, acc_elemlast
      element => acc_mesh % Elements(i)
      acc_model % CurrentElement => element
      j = GetElementFamily(element)
      if (acc_gotomask) then
        if (i <= acc_mesh % NumberOfBulkElements) then
          if (.not. GetLogical(GetBodyForce(element), 'Create Isosurface', found)) cycle
        else
          if (.not. GetLogical(GetBC(element), 'Create Isosurface', found)) cycle
        end if
      end if
      if (j <= 2) then
        cycle
      else if (j <= 4) then
        dim = max(dim, 2)
      else
        dim = max(dim, 3)
      end if
    end do
    ! The leading dimension has to agree in every partition.
    if (ParEnv % PEs > 1) then
      i = dim
      call MPI_ALLREDUCE(i, dim, 1, MPI_INTEGER, MPI_MAX, ELMER_COMM_WORLD, ierr)
    end if
    call Info('IsosurfaceSolver', 'Leading element dimension is ' // i2s(dim))
    acc_dim = dim
    if (dim == 2) then
      acc_neweltype = 303
      acc_newelnodes = 3
    else if (dim == 3) then
      acc_neweltype = 504
      acc_newelnodes = 4
    else
      call Warn('IsoSurfaceSolver', 'Isosurface mesh can be created only for 2D or 3D')
      code = 1
    end if
  end function iso_detect_dim

  ! Lines 314-459: temporal simplex mesh plus edge tables.
  subroutine iso_build_temporal() bind(c, name="iso_build_temporal")
    integer, target :: tetratotetramap(1,4), pyramidtotetramap(2,4), &
        wedgetotetramap(3,4), bricktotetramap(5,4), &
        triangletotrianglemap(1,3), quadtotrianglemap(2,3)
    integer, pointer :: map(:,:), indexes(:)
    type(element_t), pointer :: element
    integer :: i, j, k, l, n, level, isallocated
    logical :: found, crossed
    real(kind=dp) :: levelvalue
    levelvalue = acc_levelvalue
    tetratotetramap(1,:) = [1,2,3,4]
    pyramidtotetramap(1,:) = [3,5,4,1]
    pyramidtotetramap(2,:) = [3,5,2,1]
    wedgetotetramap(1,:) = [5,4,3,1]
    wedgetotetramap(2,:) = [5,3,2,1]
    wedgetotetramap(3,:) = [5,6,4,3]
    bricktotetramap(1,:) = [1,2,4,5]
    bricktotetramap(2,:) = [6,7,2,5]
    bricktotetramap(3,:) = [3,4,2,7]
    bricktotetramap(4,:) = [8,7,5,4]
    bricktotetramap(5,:) = [2,4,5,7]
    triangletotrianglemap(1,:) = [1,2,3]
    quadtotrianglemap(1,:) = [1,2,3]
    quadtotrianglemap(2,:) = [1,3,4]
    call Info('IsoSurfaceSolver', 'Creating a temporal mesh')
    do isallocated = 0, 1
      j = 0
      do i = acc_elemfirst, acc_elemlast
        element => acc_mesh % Elements(i)
        acc_model % CurrentElement => element
        if (acc_gotomask) then
          if (i <= acc_mesh % NumberOfBulkElements) then
            if (.not. GetLogical(GetBodyForce(element), 'Create Isosurface', found)) cycle
          else
            if (.not. GetLogical(GetBC(element), 'Create Isosurface', found)) cycle
          end if
        end if
        indexes => element % NodeIndexes
        n = element % TYPE % NumberOfNodes
        if (associated(acc_levelperm)) then
          acc_elemperm(1:n) = acc_levelperm(indexes)
          if (any(acc_levelperm(indexes) == 0)) cycle
        else
          acc_elemperm(1:n) = indexes(1:n)
        end if
        select case (GetElementFamily(element))
        case (3); map => triangletotrianglemap
        case (4); map => quadtotrianglemap
        case (5); map => tetratotetramap
        case (6); map => pyramidtotetramap
        case (7); map => wedgetotetramap
        case (8); map => bricktotetramap
        end select
        if (acc_leveldofs == 1) then
          acc_elemfun(1:n) = acc_levelfun(acc_elemperm(1:n))
        else
          acc_elemfun(1:n) = 0.0_dp
          do l = 1, acc_levelcomps
            acc_elemfun(1:n) = acc_elemfun(1:n) + &
                acc_levelfun(acc_leveldofs * (acc_elemperm(1:n) - 1) + l)**2
          end do
          acc_elemfun(1:n) = sqrt(acc_elemfun(1:n))
        end if
        crossed = .false.
        do level = 1, acc_nolevels
          if (acc_nolevels > 1) levelvalue = acc_levelvalues(level, 1)
          if (all(acc_elemfun(1:n) - levelvalue < 0.0_dp)) cycle
          if (all(acc_elemfun(1:n) - levelvalue > 0.0_dp)) cycle
          crossed = .true.
          exit
        end do
        if (.not. crossed) cycle
        if (isallocated == 0) then
          j = j + size(map, 1)
          cycle
        end if
        do k = 1, size(map, 1)
          j = j + 1
          acc_newelements(j) = element
          ! Detach the copied index tables before the simplex gets its own
          ! identity: FindMeshEdges reuses and overwrites those in place.
          acc_newelements(j) % EdgeIndexes => null()
          acc_newelements(j) % FaceIndexes => null()
          acc_newelements(j) % BubbleIndexes => null()
          acc_newelements(j) % DGIndexes => null()
          acc_newelements(j) % BoundaryInfo => null()
          acc_newelements(j) % PDefs => null()
          acc_newelements(j) % PropertyData => null()
          acc_newelements(j) % TYPE => GetElementType(acc_neweltype)
          acc_newelements(j) % NDOFs = acc_newelnodes
          acc_newelements(j) % BDOFs = 0
          acc_newelements(j) % DGDOFs = 0
          acc_newelements(j) % ElementIndex = j
          allocate(acc_newelements(j) % NodeIndexes(acc_newelnodes))
          do l = 1, acc_newelnodes
            acc_newelements(j) % NodeIndexes(l) = element % NodeIndexes(map(k, l))
          end do
        end do
      end do
      acc_nonewelems = j
      if (acc_nonewelems == 0) exit
      if (isallocated == 0) then
        allocate(acc_newelements(acc_nonewelems))
      end if
    end do
    ! An empty isosurface mesh is still created, so that Model % Meshes
    ! holds the same meshes in every partition.
    if (acc_nonewelems == 0) then
      call Info('IsosurfaceSolver', 'No potential elements found in this partition', level=7)
    else
      call Info('IsosurfaceSolver', 'Found ' // i2s(acc_nonewelems) // ' potential elements')
    end if
    ! Wrap the simplices into a mesh of their own. Nothing below may touch
    ! the parent mesh any more: FindMeshEdges works in place.
    if (acc_nonewelems > 0) then
      acc_tempmesh => AllocateMesh()
      acc_tempmesh % Name = 'isosurface temporal'
      deallocate(acc_tempmesh % Nodes)
      acc_tempmesh % Nodes => acc_mesh % Nodes
      acc_tempmesh % NodesOrig => acc_mesh % Nodes
      acc_tempmesh % NumberOfNodes = acc_mesh % NumberOfNodes
      acc_tempmesh % MeshDim = acc_dim
      acc_tempmesh % Elements => acc_newelements
      acc_tempmesh % NumberOfBulkElements = acc_nonewelems
      acc_tempmesh % NumberOfBoundaryElements = 0
      acc_tempmesh % MaxElementNodes = acc_newelnodes
      acc_tempmesh % MaxElementDOFs = acc_newelnodes
      call Info('IsosurfaceSolver', 'Creating mesh edges', level=9)
      if (acc_dim == 2) then
        call FindMeshEdges2D(acc_tempmesh)
      else
        call FindMeshEdges3D(acc_tempmesh)
      end if
      acc_noedges = acc_tempmesh % NumberOfEdges
    end if
  end subroutine iso_build_temporal

  ! Lines 461-512: fresh isosurface mesh shell plus output directory.
  subroutine iso_create_mesh() bind(c, name="iso_create_mesh")
    logical :: found
    call Info('IsosurfaceSolver', 'Creating a new mesh for isosurface')
    acc_isomesh => AllocateMesh()
    acc_isomesh % Name = GetString(acc_params, 'Mesh Name', found)
    if (.not. found) then
      if (acc_dim == 2) then
        acc_isomesh % Name = "isoline"
      else
        acc_isomesh % Name = "isosurf"
      end if
    end if
    acc_isomesh % Changed = .true.
    acc_isomesh % NumberOfBulkElements = 0
    acc_isomesh % NumberOfNodes = 0
    acc_exist = .true.
    acc_isomesh % OutputActive = GetLogical(acc_params, 'Isomesh Output Active', found)
    if (.not. found) acc_isomesh % OutputActive = .true.
    ! If requested number the output directories.
    if (GetLogical(acc_params, 'Mesh Numbering', found)) then
      acc_calls = acc_calls + 1
      acc_isomesh % name = trim(acc_isomesh % name) // i2s(acc_calls)
    end if
    call MakeDirectory(trim(acc_isomesh % name) // char(0))
  end subroutine iso_create_mesh

  function iso_create_nodes() result(n) bind(c, name="iso_create_nodes")
    integer(c_int) :: n
    n = acc_createnodes()
    acc_noisonodes = n
  end function iso_create_nodes

  function iso_create_surfaces() result(n) bind(c, name="iso_create_surfaces")
    integer(c_int) :: n
    n = acc_createsurfaces()
    acc_nosurfaces = n
  end function iso_create_surfaces

  ! Lines 522-524: ownership, counts, parallel numbering.
  subroutine iso_parallel_info() bind(c, name="iso_parallel_info")
    call acc_isoedgeowners()
    call acc_countowned()
    call acc_buildisoparallel()
  end subroutine iso_parallel_info

  ! Lines 526-541: drop the temporal mesh, link the new mesh in.
  subroutine iso_release_temp() bind(c, name="iso_release_temp")
    type(mesh_t), pointer :: pmesh
    call acc_releasetempmesh()
    if (allocated(acc_elemfun)) deallocate(acc_elemfun, acc_elemperm)
    pmesh => acc_model % Meshes
    do while (associated(pmesh % Next))
      pmesh => pmesh % Next
    end do
    pmesh % Next => acc_isomesh
    acc_isocreated = .true.
  end subroutine iso_release_temp

  ! Label 100 through the end: collective sizes, reports, savers, norm.
  subroutine iso_finish() bind(c, name="iso_finish")
    logical :: found
    integer :: ierr
    if (acc_isocreated) then
      acc_partmp(1) = acc_nomeshnodesowned
      acc_partmp(2) = acc_noorigelems
      acc_partmp(3) = acc_nonewelems
      acc_partmp(4) = acc_noedgesowned
      acc_partmp(5) = acc_noisonodesowned
      acc_partmp(6) = acc_nosurfaces
      acc_savepartmp = acc_partmp
    end if
    if (ParEnv % PEs > 1) then
      call MPI_ALLREDUCE(acc_partmp, acc_parsizes, 6, MPI_INTEGER, MPI_SUM, &
          ELMER_COMM_WORLD, ierr)
    else
      acc_parsizes = acc_partmp
    end if
    ! Only the partitions that really created a mesh may report on it or
    ! save it, the others just participated in the reduction above.
    if (acc_isocreated) then
      call Info('IsoSurfaceSolver', 'Information on isosurface mesh sizes')
      write (Message, '(A,I0,A)') 'Initial mesh has ', acc_parsizes(1), ' nodes'
      call Info('IsoSurfaceSolver', Message)
      write (Message, '(A,I0,A)') 'Initial mesh has ', acc_parsizes(2), ' elements'
      call Info('IsoSurfaceSolver', Message)
      write (Message, '(A,I0,A)') 'Temporal mesh has ', acc_parsizes(3), ' elements'
      call Info('IsoSurfaceSolver', Message)
      write (Message, '(A,I0,A)') 'Temporal mesh has ', acc_parsizes(4), ' edges'
      call Info('IsoSurfaceSolver', Message)
      write (Message, '(A,I0,A)') 'Isosurface mesh has ', acc_parsizes(5), ' nodes'
      call Info('IsoSurfaceSolver', Message)
      write (Message, '(A,I0,A)') 'Isosurface mesh has ', acc_parsizes(6), ' elements'
      call Info('IsoSurfaceSolver', Message)
      if (GetLogical(acc_params, 'Save Gmsh Geo File', found)) then
        call acc_savegmshgeo2d(acc_isomesh)
      end if
      if (GetLogical(acc_params, 'Save STL File', found)) then
        ! The triangles were pointed up the gradient when they were made,
        ! so the saver must not turn the patches to face away from a centre.
        call SaveSTLSurface(acc_isomesh, acc_params, KeepOrientation=.true.)
      end if
    end if
    ! Just create some variable that will act as a norm.
    if (size(acc_solver % Variable % Values) == 1) then
      acc_solver % Variable % Values = sum(acc_parsizes)
      acc_solver % Variable % Norm = sum(acc_parsizes)
    end if
  end subroutine iso_finish

  ! -- IsosurfaceSolver_Init body -------------------------------------------
  subroutine iso_init_begin(h_solver) bind(c, name="iso_init_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine iso_init_begin

  subroutine iso_init_body() bind(c, name="iso_init_body")
    character(len=MAX_NAME_LEN) :: name
    logical :: gotit
    if (.not. ListCheckPresent(acc_solver % Values, 'Variable')) then
      name = ListGetString(acc_solver % Values, 'Equation', gotit)
      call ListAddString(acc_solver % Values, 'Variable', &
          '-nooutput -global ' // trim(name) // '_var')
    end if
  end subroutine iso_init_body

  ! ======================================================================
  ! Transcribed CONTAINS internals (whole units; see module header).
  ! ======================================================================

  subroutine acc_isoedgeowners()
    implicit none
    type hashetry_t
      integer :: g1, g2, ind
      type(hashetry_t), pointer :: next => null()
    end type hashetry_t
    type hashhead_t
      type(hashetry_t), pointer :: head => null()
    end type hashhead_t
    type(hashhead_t), allocatable :: hash(:)
    type(hashetry_t), pointer :: hp, hp1
    integer, allocatable :: key(:,:), cand(:,:), ncand(:), &
        nsend(:), nrecv(:), senddispl(:), recvdispl(:), sendbuf(:), recvbuf(:)
    integer, pointer :: a(:), b(:)
    integer :: i, j, k, p, n1, n2, g1, g2, nh, h, pe, ns, nr, maxcand, ne, ierr
    if (allocated(acc_edgeowner)) deallocate(acc_edgeowner, acc_edgenbrs, acc_nedgenbrs)
    allocate(acc_edgeowner(max(acc_noedges, 1)), &
        acc_edgenbrs(max(ParEnv % PEs, 1), max(acc_noedges, 1)), &
        acc_nedgenbrs(max(acc_noedges, 1)))
    acc_edgeowner = ParEnv % MyPe
    acc_edgenbrs = 0
    acc_nedgenbrs = 1
    do i = 1, acc_noedges
      acc_edgenbrs(1, i) = ParEnv % MyPe
    end do
    ! The only exit that every partition takes together.
    if (ParEnv % PEs <= 1) return
    ne = acc_noedges
    if (.not. associated(acc_tempmesh)) ne = 0
    if (.not. associated(acc_mesh % ParallelInfo % NeighbourList)) ne = 0
    pe = ParEnv % PEs
    maxcand = pe
    allocate(key(2, max(ne, 1)), cand(maxcand, max(ne, 1)), ncand(max(ne, 1)))
    ncand = 0
    do i = 1, ne
      n1 = acc_tempmesh % Edges(i) % NodeIndexes(1)
      n2 = acc_tempmesh % Edges(i) % NodeIndexes(2)
      g1 = acc_mesh % ParallelInfo % GlobalDOFs(n1)
      g2 = acc_mesh % ParallelInfo % GlobalDOFs(n2)
      key(1, i) = min(g1, g2)
      key(2, i) = max(g1, g2)
      a => acc_mesh % ParallelInfo % NeighbourList(n1) % Neighbours
      b => acc_mesh % ParallelInfo % NeighbourList(n2) % Neighbours
      do j = 1, size(a)
        if (a(j) == ParEnv % MyPe) cycle
        do k = 1, size(b)
          if (a(j) == b(k)) then
            ncand(i) = ncand(i) + 1
            cand(ncand(i), i) = a(j)
            exit
          end if
        end do
      end do
    end do
    ! Hash my own edges by their key, to recognise the ones sent back.
    nh = 2 * ne + 1
    allocate(hash(nh))
    do i = 1, ne
      h = modulo(int(key(1, i), 8) * 73856093_8 + int(key(2, i), 8) * 19349663_8, &
          int(nh, 8)) + 1
      allocate(hp)
      hp % g1 = key(1, i); hp % g2 = key(2, i); hp % ind = i
      hp % next => hash(h) % head
      hash(h) % head => hp
    end do
    allocate(nsend(pe), nrecv(pe), senddispl(pe), recvdispl(pe))
    nsend = 0
    do i = 1, ne
      do j = 1, ncand(i)
        nsend(cand(j, i) + 1) = nsend(cand(j, i) + 1) + 2
      end do
    end do
    call MPI_ALLTOALL(nsend, 1, MPI_INTEGER, nrecv, 1, MPI_INTEGER, &
        ELMER_COMM_WORLD, ierr)
    senddispl(1) = 0; recvdispl(1) = 0
    do p = 2, pe
      senddispl(p) = senddispl(p - 1) + nsend(p - 1)
      recvdispl(p) = recvdispl(p - 1) + nrecv(p - 1)
    end do
    ns = sum(nsend); nr = sum(nrecv)
    allocate(sendbuf(max(ns, 1)), recvbuf(max(nr, 1)))
    nsend = 0
    do i = 1, ne
      do j = 1, ncand(i)
        p = cand(j, i) + 1
        sendbuf(senddispl(p) + nsend(p) + 1) = key(1, i)
        sendbuf(senddispl(p) + nsend(p) + 2) = key(2, i)
        nsend(p) = nsend(p) + 2
      end do
    end do
    call MPI_ALLTOALLV(sendbuf, nsend, senddispl, MPI_INTEGER, &
        recvbuf, nrecv, recvdispl, MPI_INTEGER, ELMER_COMM_WORLD, ierr)
    ! An edge that came back from p is one we both hold.
    do p = 1, pe
      do k = 1, nrecv(p), 2
        g1 = recvbuf(recvdispl(p) + k)
        g2 = recvbuf(recvdispl(p) + k + 1)
        h = modulo(int(g1, 8) * 73856093_8 + int(g2, 8) * 19349663_8, &
            int(nh, 8)) + 1
        hp => hash(h) % head
        do while (associated(hp))
          if (hp % g1 == g1 .and. hp % g2 == g2) then
            i = hp % ind
            acc_edgeowner(i) = min(acc_edgeowner(i), p - 1)
            acc_nedgenbrs(i) = acc_nedgenbrs(i) + 1
            acc_edgenbrs(acc_nedgenbrs(i), i) = p - 1
            exit
          end if
          hp => hp % next
        end do
      end do
    end do
    ! The neighbour lists are expected sorted, the owner first.
    do i = 1, ne
      do j = 2, acc_nedgenbrs(i)
        k = acc_edgenbrs(j, i)
        p = j - 1
        do while (p >= 1)
          if (acc_edgenbrs(p, i) <= k) exit
          acc_edgenbrs(p + 1, i) = acc_edgenbrs(p, i)
          p = p - 1
        end do
        acc_edgenbrs(p + 1, i) = k
      end do
    end do
    do i = 1, nh
      hp => hash(i) % head
      do while (associated(hp))
        hp1 => hp % next; deallocate(hp); hp => hp1
      end do
    end do
    deallocate(hash, key, cand, ncand, nsend, nrecv, senddispl, recvdispl, &
        sendbuf, recvbuf)
  end subroutine acc_isoedgeowners

  subroutine acc_buildisoparallel()
    implicit none
    integer, allocatable :: nsend(:), nrecv(:), senddispl(:), recvdispl(:), &
        sendbuf(:), recvbuf(:)
    integer :: i, j, k, p, e, pe, nown, offset, ns, nr, g1, g2, maxid, niso, ierr
    type(parallelInfo_t), pointer :: pi
    ! The only exit that all partitions take together.
    if (ParEnv % PEs <= 1) return
    niso = acc_noisonodes
    if (.not. associated(acc_isomesh)) niso = 0
    if (.not. allocated(acc_isoedge)) niso = 0
    if (.not. allocated(acc_edgeowner)) niso = 0
    if (.not. associated(acc_mesh % ParallelInfo % NeighbourList)) niso = 0
    pe = ParEnv % PEs
    pi => acc_isomesh % ParallelInfo
    allocate(pi % GlobalDOFs(max(acc_noisonodes, 1)), &
        pi % GInterface(max(acc_noisonodes, 1)), &
        pi % NeighbourList(max(acc_noisonodes, 1)))
    pi % GlobalDOFs = 0
    pi % GInterface = .false.
    do i = 1, max(acc_noisonodes, 1)
      nullify(pi % NeighbourList(i) % Neighbours)
    end do
    nown = 0
    do i = 1, niso
      e = acc_isoedge(i)
      if (e <= 0) cycle
      k = acc_nedgenbrs(e)
      allocate(pi % NeighbourList(i) % Neighbours(k))
      pi % NeighbourList(i) % Neighbours = acc_edgenbrs(1:k, e)
      pi % GInterface(i) = (k > 1)
      if (acc_edgeowner(e) == ParEnv % MyPe) nown = nown + 1
    end do
    ! Where this partition's block of numbers starts.
    call MPI_SCAN(nown, offset, 1, MPI_INTEGER, MPI_SUM, ELMER_COMM_WORLD, ierr)
    offset = offset - nown
    k = 0
    do i = 1, niso
      e = acc_isoedge(i)
      if (e <= 0) cycle
      if (acc_edgeowner(e) /= ParEnv % MyPe) cycle
      k = k + 1
      pi % GlobalDOFs(i) = offset + k
    end do
    ! Tell the sharers the numbers of the nodes they do not own.
    allocate(nsend(pe), nrecv(pe), senddispl(pe), recvdispl(pe))
    nsend = 0
    do i = 1, niso
      e = acc_isoedge(i)
      if (e <= 0) cycle
      if (acc_edgeowner(e) /= ParEnv % MyPe) cycle
      do j = 1, acc_nedgenbrs(e)
        if (acc_edgenbrs(j, e) == ParEnv % MyPe) cycle
        nsend(acc_edgenbrs(j, e) + 1) = nsend(acc_edgenbrs(j, e) + 1) + 4
      end do
    end do
    call MPI_ALLTOALL(nsend, 1, MPI_INTEGER, nrecv, 1, MPI_INTEGER, &
        ELMER_COMM_WORLD, ierr)
    senddispl(1) = 0; recvdispl(1) = 0
    do p = 2, pe
      senddispl(p) = senddispl(p - 1) + nsend(p - 1)
      recvdispl(p) = recvdispl(p - 1) + nrecv(p - 1)
    end do
    ns = sum(nsend); nr = sum(nrecv)
    allocate(sendbuf(max(ns, 1)), recvbuf(max(nr, 1)))
    nsend = 0
    do i = 1, niso
      e = acc_isoedge(i)
      if (e <= 0) cycle
      if (acc_edgeowner(e) /= ParEnv % MyPe) cycle
      g1 = acc_mesh % ParallelInfo % GlobalDOFs(acc_invperm(1, i))
      g2 = acc_mesh % ParallelInfo % GlobalDOFs(acc_invperm(2, i))
      do j = 1, acc_nedgenbrs(e)
        if (acc_edgenbrs(j, e) == ParEnv % MyPe) cycle
        p = acc_edgenbrs(j, e) + 1
        sendbuf(senddispl(p) + nsend(p) + 1) = min(g1, g2)
        sendbuf(senddispl(p) + nsend(p) + 2) = max(g1, g2)
        sendbuf(senddispl(p) + nsend(p) + 3) = acc_isolevel(i)
        sendbuf(senddispl(p) + nsend(p) + 4) = pi % GlobalDOFs(i)
        nsend(p) = nsend(p) + 4
      end do
    end do
    call MPI_ALLTOALLV(sendbuf, nsend, senddispl, MPI_INTEGER, &
        recvbuf, nrecv, recvdispl, MPI_INTEGER, ELMER_COMM_WORLD, ierr)
    do k = 1, nr, 4
      do i = 1, niso
        if (pi % GlobalDOFs(i) /= 0) cycle
        if (acc_isolevel(i) /= recvbuf(k + 2)) cycle
        g1 = acc_mesh % ParallelInfo % GlobalDOFs(acc_invperm(1, i))
        g2 = acc_mesh % ParallelInfo % GlobalDOFs(acc_invperm(2, i))
        if (min(g1, g2) /= recvbuf(k) .or. max(g1, g2) /= recvbuf(k + 1)) cycle
        pi % GlobalDOFs(i) = recvbuf(k + 3)
        exit
      end do
    end do
    deallocate(nsend, nrecv, senddispl, recvdispl, sendbuf, recvbuf)
    ! Every node should have a number by now.
    j = 0
    do i = 1, niso
      if (pi % GlobalDOFs(i) == 0) j = j + 1
    end do
    if (j > 0) then
      call Warn('IsosurfaceSolver', 'Left ' // i2s(j) // ' isosurface nodes unnumbered!')
    end if
    maxid = 0
    if (niso > 0) maxid = maxval(pi % GlobalDOFs(1:niso))
    j = maxid
    call MPI_ALLREDUCE(j, maxid, 1, MPI_INTEGER, MPI_MAX, ELMER_COMM_WORLD, ierr)
    call Info('IsosurfaceSolver', 'Isosurface mesh global nodes: ' // i2s(maxid), level=7)
  end subroutine acc_buildisoparallel

  subroutine acc_countowned()
    implicit none
    integer :: i
    acc_nomeshnodesowned = acc_mesh % NumberOfNodes
    acc_noedgesowned = acc_noedges
    acc_noisonodesowned = acc_noisonodes
    if (ParEnv % PEs <= 1) return
    if (associated(acc_mesh % ParallelInfo % NeighbourList)) then
      acc_nomeshnodesowned = 0
      do i = 1, acc_mesh % NumberOfNodes
        if (acc_mesh % ParallelInfo % NeighbourList(i) % Neighbours(1) == ParEnv % MyPe) &
            acc_nomeshnodesowned = acc_nomeshnodesowned + 1
      end do
    end if
    if (allocated(acc_edgeowner)) then
      acc_noedgesowned = 0
      do i = 1, acc_noedges
        if (acc_edgeowner(i) == ParEnv % MyPe) acc_noedgesowned = acc_noedgesowned + 1
      end do
      if (allocated(acc_isoedge)) then
        acc_noisonodesowned = 0
        do i = 1, acc_noisonodes
          if (acc_isoedge(i) <= 0) cycle
          if (acc_edgeowner(acc_isoedge(i)) == ParEnv % MyPe) &
              acc_noisonodesowned = acc_noisonodesowned + 1
        end do
      end if
    end if
  end subroutine acc_countowned

  subroutine acc_releaseisomesh()
    implicit none
    type(mesh_t), pointer :: pmesh
    if (.not. acc_exist) return
    call Info('IsosurfaceSolver', 'Releasing previous mesh', level=9)
    pmesh => acc_model % Meshes
    do while (associated(pmesh % Next))
      if (associated(pmesh % next, acc_isomesh)) then
        pmesh % next => acc_isomesh % Next
        exit
      end if
      pmesh => pmesh % next
    end do
    call ReleaseMesh(acc_isomesh)
    acc_isomesh => null()
    acc_exist = .false.
  end subroutine acc_releaseisomesh

  subroutine acc_getinterpolants()
    implicit none
    integer :: i
    logical :: found
    acc_ints = 0
    do i = 1, ACC_MAXINTS
      acc_vname(i) = GetString(acc_params, ComponentName('Isosurface interpolant', i), found)
      if (.not. found) acc_vname(i) = GetString(acc_params, &
          ComponentName('Variable', i), found)
      if (.not. found) exit
      acc_vfulltab(i) % Variable => VariableGet(acc_mesh % Variables, acc_vname(i))
      if (.not. associated(acc_vfulltab(i) % Variable)) exit
      acc_ints = i
    end do
    if (acc_ints == ACC_MAXINTS) then
      call Warn('IsosurfaceSolver', 'Only the first ' // i2s(ACC_MAXINTS) // &
          ' interpolants are honoured!')
    end if
  end subroutine acc_getinterpolants

  subroutine acc_getisointerpolants()
    implicit none
    integer :: i
    do i = 1, acc_ints
      acc_visotab(i) % Variable => VariableGet(acc_isomesh % Variables, acc_vname(i), &
          ThisOnly=.true.)
    end do
  end subroutine acc_getisointerpolants

  subroutine acc_releasetempmesh()
    implicit none
    integer :: i
    if (.not. associated(acc_tempmesh)) return
    ! This releases the edge table and the EdgeIndexes of the simplices.
    ! It must be done while the elements are still in place.
    call ReleaseMeshEdgeTables(acc_tempmesh)
    if (associated(acc_tempmesh % Elements)) then
      do i = 1, acc_tempmesh % NumberOfBulkElements
        if (associated(acc_tempmesh % Elements(i) % NodeIndexes)) &
            deallocate(acc_tempmesh % Elements(i) % NodeIndexes)
        acc_tempmesh % Elements(i) % NodeIndexes => null()
      end do
      deallocate(acc_tempmesh % Elements)
    end if
    acc_tempmesh % Elements => null()
    acc_newelements => null()
    ! Owned by the parent mesh!
    acc_tempmesh % Nodes => null()
    acc_tempmesh % NodesOrig => null()
    deallocate(acc_tempmesh)
    acc_tempmesh => null()
  end subroutine acc_releasetempmesh

  function acc_createnodes() result(noisonodes)
    implicit none
    integer :: noisonodes
    type(element_t), pointer :: edge
    integer :: i, j, k, l, n1, n2, m1, m2, level, isallocated
    real(kind=dp) :: t, t1, t2, x1, x2, y1, y2, z1, z2
    real(kind=dp) :: levelvalue
    type(variable_t), pointer :: vfull, viso
    integer, pointer :: vperm(:)
    real(kind=dp), pointer :: vals(:)
    levelvalue = acc_levelvalue
    noisonodes = 0
    call acc_getinterpolants()
    ! The isosurface may well miss this partition altogether. Give it an
    ! empty mesh anyway, complete with the same variable set.
    acc_noedges = 0
    if (associated(acc_tempmesh)) acc_noedges = acc_tempmesh % NumberOfEdges
    if (acc_noedges == 0) then
      call Info('IsoSurfaceSolver', 'Creating an empty isosurface mesh', level=8)
      allocate(acc_isomesh % Nodes % x(0), acc_isomesh % Nodes % y(0), &
          acc_isomesh % Nodes % z(0))
      acc_isomesh % NumberOfNodes = 0
      acc_isomesh % Nodes % NumberOfNodes = 0
      acc_isomesh % MeshDim = acc_dim
      allocate(acc_isomesh % Elements(0))
      acc_isomesh % NumberOfBulkElements = 0
      acc_isomesh % NumberOfBoundaryElements = 0
      acc_isomesh % NumberOfEdges = 0
      acc_isomesh % NumberOfFaces = 0
      call VariableAdd(acc_isomesh % Variables, acc_isomesh, acc_solver, &
          'Coordinate 1', 1, acc_isomesh % Nodes % x)
      call VariableAdd(acc_isomesh % Variables, acc_isomesh, acc_solver, &
          'Coordinate 2', 1, acc_isomesh % Nodes % y)
      call VariableAdd(acc_isomesh % Variables, acc_isomesh, acc_solver, &
          'Coordinate 3', 1, acc_isomesh % Nodes % z)
      vfull => VariableGet(acc_mesh % Variables, 'Time')
      call VariableAdd(acc_isomesh % Variables, acc_isomesh, acc_solver, 'Time', 1, &
          vfull % Values)
      do i = 1, acc_ints
        vfull => acc_vfulltab(i) % Variable
        nullify(vperm, vals)
        allocate(vperm(1), vals(vfull % DOFs))
        vperm = 0
        vals = 0.0_dp
        call Info('IsoSurfaceSolver', 'Creating dummy variable ' // trim(acc_vname(i)), level=8)
        call VariableAddVector(acc_isomesh % Variables, acc_isomesh, acc_solver, &
            trim(acc_vname(i)), vfull % DOFs, vals, vperm)
      end do
      call acc_getisointerpolants()
      return
    end if
    ! Loop over edges and check for different signs.
    do isallocated = 0, 1
      j = 0
      do level = 1, acc_nolevels
        if (acc_nolevels > 1) levelvalue = acc_levelvalues(level, 1)
        do i = 1, acc_noedges
          edge => acc_tempmesh % Edges(i)
          n1 = edge % NodeIndexes(1)
          n2 = edge % NodeIndexes(2)
          if (associated(acc_levelperm)) then
            m1 = acc_levelperm(n1)
            m2 = acc_levelperm(n2)
          else
            m1 = n1
            m2 = n2
          end if
          if (m1 <= 0 .or. m2 <= 0) cycle
          if (acc_leveldofs == 1) then
            t1 = acc_levelfun(m1)
            t2 = acc_levelfun(m2)
          else
            t1 = 0.0_dp
            t2 = 0.0_dp
            do k = 1, acc_levelcomps
              t1 = t1 + acc_levelfun(acc_leveldofs * (m1 - 1) + k)**2
              t2 = t2 + acc_levelfun(acc_leveldofs * (m2 - 1) + k)**2
            end do
            t1 = sqrt(t1)
            t2 = sqrt(t2)
          end if
          t1 = t1 - levelvalue
          t2 = t2 - levelvalue
          if (abs(t1 - t2) < tiny(t1)) cycle
          if (t1 * t2 > 0.0_dp) cycle
          j = j + 1
          if (isallocated == 0) cycle
          t = abs(t1 / (t2 - t1))
          ! The parent nodes the isosurface node sits between.
          acc_interpolant(j) = t
          acc_invperm(1, j) = n1
          acc_invperm(2, j) = n2
          x1 = acc_mesh % Nodes % x(n1)
          x2 = acc_mesh % Nodes % x(n2)
          y1 = acc_mesh % Nodes % y(n1)
          y2 = acc_mesh % Nodes % y(n2)
          z1 = acc_mesh % Nodes % z(n1)
          z2 = acc_mesh % Nodes % z(n2)
          acc_eperm(i) = j
          acc_isoedge(j) = i
          acc_isolevel(j) = level
          acc_isomesh % Nodes % x(j) = (1 - t) * x1 + t * x2
          acc_isomesh % Nodes % y(j) = (1 - t) * y1 + t * y2
          acc_isomesh % Nodes % z(j) = (1 - t) * z1 + t * z2
          do k = 1, acc_ints
            vfull => acc_vfulltab(k) % Variable
            viso => acc_visotab(k) % Variable
            if (.not. associated(viso)) cycle
            if (associated(vfull % Perm)) then
              m1 = vfull % Perm(n1)
              m2 = vfull % Perm(n2)
            else
              m1 = n1
              m2 = n2
            end if
            if (m1 <= 0 .or. m2 <= 0) cycle
            do l = 1, vfull % DOFs
              x1 = vfull % Values(vfull % DOFs * (m1 - 1) + l)
              x2 = vfull % Values(vfull % DOFs * (m2 - 1) + l)
              viso % Values(vfull % DOFs * (j - 1) + l) = (1 - t) * x1 + t * x2
            end do
          end do
        end do
      end do
      if (isallocated == 0) then
        call Info('IsosurfaceSolver', 'Creating ' // i2s(j) // ' nodes for isosurface')
        noisonodes = j
        acc_isomesh % NumberOfNodes = j
        acc_isomesh % Nodes % NumberOfNodes = j
        acc_isomesh % MeshDim = acc_dim
        allocate(acc_isomesh % Nodes % x(j))
        allocate(acc_isomesh % Nodes % y(j))
        allocate(acc_isomesh % Nodes % z(j))
        ! Gives the index of the node sitting on an edge ...
        allocate(acc_eperm(acc_noedges))
        acc_eperm = 0
        ! ... and the other way round, which edge a node sits on.
        if (allocated(acc_isoedge)) deallocate(acc_isoedge, acc_isolevel)
        allocate(acc_isoedge(j), acc_isolevel(j))
        acc_isoedge = 0
        acc_isolevel = 0
        call VariableAdd(acc_isomesh % Variables, acc_isomesh, acc_solver, &
            'Coordinate 1', 1, acc_isomesh % Nodes % x)
        call VariableAdd(acc_isomesh % Variables, acc_isomesh, acc_solver, &
            'Coordinate 2', 1, acc_isomesh % Nodes % y)
        call VariableAdd(acc_isomesh % Variables, acc_isomesh, acc_solver, &
            'Coordinate 3', 1, acc_isomesh % Nodes % z)
        vfull => VariableGet(acc_mesh % Variables, 'Time')
        call VariableAdd(acc_isomesh % Variables, acc_isomesh, acc_solver, 'Time', 1, &
            vfull % Values)
        do k = 1, acc_ints
          vfull => acc_vfulltab(k) % Variable
          nullify(vperm, vals)
          allocate(vperm(j), vals(vfull % DOFs * j))
          vperm = [(i, i = 1, j)]
          vals = 0.0_dp
          call Info('IsoSurfaceSolver', 'Creating variable ' // trim(acc_vname(k)))
          call VariableAddVector(acc_isomesh % Variables, acc_isomesh, acc_solver, &
              trim(acc_vname(k)), vfull % DOFs, vals, vperm)
        end do
        ! Now that they exist, bind the isosurface side of the interpolants.
        call acc_getisointerpolants()
        j = 0
        if (allocated(acc_invperm)) deallocate(acc_invperm, acc_interpolant)
        allocate(acc_invperm(2, noisonodes), acc_interpolant(noisonodes))
        acc_invperm = 0
        acc_interpolant = 0.0_dp
      end if
    end do
  end function acc_createnodes

  subroutine acc_remap()
    implicit none
    integer :: j, k, l, n1, n2, m1, m2
    real(kind=dp) :: t, x1, x2
    type(variable_t), pointer :: vfull, viso
    ! If there were no intersections then the dummy field has already
    ! been created.
    if (.not. allocated(acc_invperm)) return
    ! Resolve the fields once, not per node.
    call acc_getinterpolants()
    call acc_getisointerpolants()
    ! If the mesh is moving then remap also the coordinate values.
    if (acc_movingmesh) then
      do j = 1, acc_isomesh % NumberOfNodes
        n1 = acc_invperm(1, j)
        n2 = acc_invperm(2, j)
        t = acc_interpolant(j)
        acc_isomesh % Nodes % x(j) = &
            (1 - t) * acc_mesh % Nodes % x(n1) + t * acc_mesh % Nodes % x(n2)
        acc_isomesh % Nodes % y(j) = &
            (1 - t) * acc_mesh % Nodes % y(n1) + t * acc_mesh % Nodes % y(n2)
        acc_isomesh % Nodes % z(j) = &
            (1 - t) * acc_mesh % Nodes % z(n1) + t * acc_mesh % Nodes % z(n2)
      end do
    end if
    ! Loop over the predefined nodes, field by field.
    do k = 1, acc_ints
      vfull => acc_vfulltab(k) % Variable
      viso => acc_visotab(k) % Variable
      if (.not. associated(viso)) cycle
      do j = 1, acc_isomesh % NumberOfNodes
        n1 = acc_invperm(1, j)
        n2 = acc_invperm(2, j)
        t = acc_interpolant(j)
        if (associated(vfull % Perm)) then
          m1 = vfull % Perm(n1)
          m2 = vfull % Perm(n2)
        else
          m1 = n1
          m2 = n2
        end if
        if (m1 <= 0 .or. m2 <= 0) cycle
        do l = 1, vfull % DOFs
          x1 = vfull % Values(vfull % DOFs * (m1 - 1) + l)
          x2 = vfull % Values(vfull % DOFs * (m2 - 1) + l)
          viso % Values(viso % DOFs * (j - 1) + l) = (1 - t) * x1 + t * x2
        end do
      end do
    end do
  end subroutine acc_remap

  function acc_createsurfaces() result(nosurfaces)
    implicit none
    integer :: nosurfaces
    type(element_t), pointer :: element, defelem
    real(kind=dp) :: f(4)
    integer :: i, j, k, n, triangles(2, 3), line(2)
    integer :: neweltype, newelnodes
    integer :: level, isallocated
    real(kind=dp) :: levelvalue
    levelvalue = acc_levelvalue
    nosurfaces = 0
    if (acc_noisonodes == 0) return
    if (acc_dim == 3) then
      neweltype = 303
      newelnodes = 3
    else
      neweltype = 202
      newelnodes = 2
    end if
    do isallocated = 0, 1
      k = 0
      do level = 1, acc_nolevels
        if (acc_nolevels > 1) levelvalue = acc_levelvalues(level, 1)
        do i = 1, acc_tempmesh % NumberOfBulkElements
          element => acc_tempmesh % Elements(i)
          n = element % TYPE % NumberOfNodes
          if (associated(acc_levelperm)) then
            acc_elemperm(1:n) = acc_levelperm(element % NodeIndexes)
          else
            acc_elemperm(1:n) = element % NodeIndexes
          end if
          if (acc_leveldofs == 1) then
            f(1:n) = acc_levelfun(acc_elemperm(1:n))
          else
            f(1:n) = 0.0_dp
            do j = 1, acc_levelcomps
              f(1:n) = f(1:n) + acc_levelfun(acc_leveldofs * (acc_elemperm(1:n) - 1) + j)**2
            end do
            f(1:n) = sqrt(f(1:n))
          end if
          f(1:n) = f(1:n) - levelvalue
          if (acc_dim == 2) then
            ! An edge that carries no intersection has eperm zero, which can
            ! happen for degenerate sign patterns. Skip those, as the 3D
            ! branch below does, rather than emit a zero node index.
            if (acc_createlinefromtriangle(element, f, line)) then
              if (all(line > 0)) then
                k = k + 1
                if (isallocated == 1) then
                  acc_isomesh % Elements(k) % NodeIndexes = line
                  acc_isomesh % Elements(k) % BodyId = level
                  acc_isomesh % Elements(k) % PartIndex = element % PartIndex
                end if
              end if
            end if
          else
            n = acc_createsurfacefromtetra(element, f, triangles)
            do j = 1, n
              if (all(triangles(j, :) > 0)) then
                k = k + 1
                if (isallocated == 1) then
                  call acc_orientisotriangle(element, f, triangles(j, :))
                  acc_isomesh % Elements(k) % NodeIndexes = triangles(j, :)
                  acc_isomesh % Elements(k) % BodyId = level
                  acc_isomesh % Elements(k) % PartIndex = element % PartIndex
                end if
              end if
            end do
          end if
        end do
      end do
      if (k == 0) return
      if (isallocated == 0) then
        nosurfaces = k
        allocate(acc_isomesh % Elements(k))
        acc_isomesh % MeshDim = acc_dim
        acc_isomesh % NumberOfBulkElements = k
        acc_isomesh % NumberOfFaces = 0
        acc_isomesh % NumberOfEdges = 0
        acc_isomesh % NumberOfBoundaryElements = 0
        defelem => AllocateElement()
        defelem % TYPE => GetElementType(neweltype)
        ! The savers treat an element whose PartIndex is not their own as
        ! a halo element and skip it by default, so the -1 that
        ! AllocateElement gives would write out empty pieces. Each element
        ! inherits the owner of the parent element it was cut from, below.
        defelem % PartIndex = ParEnv % MyPe
        do i = 1, k
          acc_isomesh % Elements(i) = defelem
          acc_isomesh % Elements(i) % ElementIndex = i
          allocate(acc_isomesh % Elements(i) % NodeIndexes(newelnodes))
        end do
      end if
    end do
    deallocate(defelem)
  end function acc_createsurfaces

  subroutine acc_orientisotriangle(element, f, tri)
    implicit none
    type(element_t), pointer :: element
    real(kind=dp) :: f(:)
    integer :: tri(3)
    real(kind=dp) :: p(3, 3), u(3), v(3), nor(3), xp(3), xm(3)
    integer :: i, np, nm, ni, itmp
    do i = 1, 3
      p(1, i) = acc_isomesh % Nodes % x(tri(i))
      p(2, i) = acc_isomesh % Nodes % y(tri(i))
      p(3, i) = acc_isomesh % Nodes % z(tri(i))
    end do
    u = p(:, 2) - p(:, 1)
    v = p(:, 3) - p(:, 1)
    nor(1) = u(2) * v(3) - u(3) * v(2)
    nor(2) = u(3) * v(1) - u(1) * v(3)
    nor(3) = u(1) * v(2) - u(2) * v(1)
    ! Where the field is above the level, and where below.
    xp = 0.0_dp; xm = 0.0_dp; np = 0; nm = 0
    do i = 1, element % TYPE % NumberOfNodes
      ni = element % NodeIndexes(i)
      if (f(i) > 0.0_dp) then
        np = np + 1
        xp(1) = xp(1) + acc_mesh % Nodes % x(ni)
        xp(2) = xp(2) + acc_mesh % Nodes % y(ni)
        xp(3) = xp(3) + acc_mesh % Nodes % z(ni)
      else
        nm = nm + 1
        xm(1) = xm(1) + acc_mesh % Nodes % x(ni)
        xm(2) = xm(2) + acc_mesh % Nodes % y(ni)
        xm(3) = xm(3) + acc_mesh % Nodes % z(ni)
      end if
    end do
    if (np == 0 .or. nm == 0) return
    if (sum(nor * (xp / np - xm / nm)) < 0.0_dp) then
      itmp = tri(2); tri(2) = tri(3); tri(3) = itmp
    end if
  end subroutine acc_orientisotriangle

  function acc_createsurfacefromtetra(tetra, f, surf) result(scount)
    implicit none
    type(element_t) :: tetra
    real(kind=dp) :: f(4)
    integer :: scount, surf(2, 3)
    logical :: s1, s2, s3, s4
    integer :: s, indexes(6)
    scount = 0
    indexes = acc_eperm(tetra % EdgeIndexes)
    s1 = f(1) > 0.0_dp
    s2 = f(2) > 0.0_dp
    s3 = f(3) > 0.0_dp
    s4 = f(4) > 0.0_dp
    s = 0
    if (s1) s = s + 1
    if (s2) s = s + 1
    if (s3) s = s + 1
    if (s4) s = s + 1
    if (s == 0 .or. s == 4) return
    if (s == 1 .or. s == 3) then
      scount = 1
      if ((s == 1 .and. s1) .or. (s == 3 .and. .not. s1)) then
        surf(1, 1) = indexes(1)
        surf(1, 2) = indexes(3)
        surf(1, 3) = indexes(4)
      else if ((s == 1 .and. s2) .or. (s == 3 .and. .not. s2)) then
        surf(1, 1) = indexes(1)
        surf(1, 2) = indexes(2)
        surf(1, 3) = indexes(5)
      else if ((s == 1 .and. s3) .or. (s == 3 .and. .not. s3)) then
        surf(1, 1) = indexes(2)
        surf(1, 2) = indexes(3)
        surf(1, 3) = indexes(6)
      else if ((s == 1 .and. s4) .or. (s == 3 .and. .not. s4)) then
        surf(1, 1) = indexes(4)
        surf(1, 2) = indexes(5)
        surf(1, 3) = indexes(6)
      else
        scount = 0
      end if
    else
      scount = 2
      if ((s1 .and. s2) .or. (.not. s1 .and. .not. s2)) then
        surf(1, 1) = indexes(3)
        surf(1, 2) = indexes(2)
        surf(1, 3) = indexes(5)
        surf(2, 1) = indexes(3)
        surf(2, 2) = indexes(5)
        surf(2, 3) = indexes(4)
      else if ((s1 .and. s3) .or. (.not. s1 .and. .not. s3)) then
        surf(1, 1) = indexes(1)
        surf(1, 2) = indexes(2)
        surf(1, 3) = indexes(6)
        surf(2, 1) = indexes(1)
        surf(2, 2) = indexes(6)
        surf(2, 3) = indexes(4)
      else if ((s1 .and. s4) .or. (.not. s1 .and. .not. s4)) then
        surf(1, 1) = indexes(1)
        surf(1, 2) = indexes(5)
        surf(1, 3) = indexes(6)
        surf(2, 1) = indexes(1)
        surf(2, 2) = indexes(6)
        surf(2, 3) = indexes(3)
      else
        scount = 0
      end if
    end if
  end function acc_createsurfacefromtetra

  function acc_createlinefromtriangle(triangle, f, line) result(gotline)
    implicit none
    type(element_t) :: triangle
    real(kind=dp) :: f(3)
    integer :: line(2)
    logical :: gotline
    logical :: s1, s2, s3
    integer :: s, indexes(3)
    gotline = .false.
    indexes = acc_eperm(triangle % EdgeIndexes)
    s1 = f(1) > 0.0_dp
    s2 = f(2) > 0.0_dp
    s3 = f(3) > 0.0_dp
    s = 0
    if (s1) s = s + 1
    if (s2) s = s + 1
    if (s3) s = s + 1
    if (s == 0 .or. s == 3) return
    gotline = .true.
    if (s == 1) then
      if (s1) then
        line(1) = indexes(1)
        line(2) = indexes(3)
      else if (s2) then
        line(1) = indexes(1)
        line(2) = indexes(2)
      else
        line(1) = indexes(2)
        line(2) = indexes(3)
      end if
    else
      if (.not. s1) then
        line(1) = indexes(1)
        line(2) = indexes(3)
      else if (.not. s2) then
        line(1) = indexes(1)
        line(2) = indexes(2)
      else
        line(1) = indexes(2)
        line(2) = indexes(3)
      end if
    end if
  end function acc_createlinefromtriangle

  subroutine acc_gathergeolines()
    implicit none
    integer, allocatable :: nrecv(:), displ(:), gid(:), gidall(:), &
        segbuf(:), segall(:)
    real(kind=dp), allocatable :: xyz(:), xyzall(:)
    integer :: i, j, k, p, pe, nown, nloc, ntot, ierr2
    type(element_t), pointer :: elem
    acc_ngeonodes = 0
    acc_ngeolines = 0
    if (allocated(acc_geox)) deallocate(acc_geox, acc_geoy, acc_geoz)
    if (allocated(acc_geoline)) deallocate(acc_geoline)
    if (ParEnv % PEs <= 1) then
      acc_ngeonodes = acc_isomesh % NumberOfNodes
      acc_ngeolines = acc_isomesh % NumberOfBulkElements
      allocate(acc_geox(max(acc_ngeonodes, 1)), acc_geoy(max(acc_ngeonodes, 1)), &
          acc_geoz(max(acc_ngeonodes, 1)), acc_geoline(2, max(acc_ngeolines, 1)))
      do i = 1, acc_ngeonodes
        acc_geox(i) = acc_isomesh % Nodes % x(i)
        acc_geoy(i) = acc_isomesh % Nodes % y(i)
        acc_geoz(i) = acc_isomesh % Nodes % z(i)
      end do
      do i = 1, acc_ngeolines
        elem => acc_isomesh % Elements(i)
        if (elem % TYPE % ElementCode /= 202) then
          call Fatal('SaveGmshGeo2D', 'Only elements of type 202 can be saved in Gmsh geo format!')
        end if
        acc_geoline(1:2, i) = elem % NodeIndexes(1:2)
      end do
      return
    end if
    pe = ParEnv % PEs
    ! The nodes this partition owns, and the total over all of them.
    nown = 0
    if (associated(acc_isomesh % ParallelInfo % NeighbourList)) then
      do i = 1, acc_isomesh % NumberOfNodes
        if (acc_isomesh % ParallelInfo % NeighbourList(i) % Neighbours(1) == &
            ParEnv % MyPe) nown = nown + 1
      end do
    end if
    allocate(gid(max(nown, 1)), xyz(max(3 * nown, 1)))
    k = 0
    if (associated(acc_isomesh % ParallelInfo % NeighbourList)) then
      do i = 1, acc_isomesh % NumberOfNodes
        if (acc_isomesh % ParallelInfo % NeighbourList(i) % Neighbours(1) /= &
            ParEnv % MyPe) cycle
        k = k + 1
        gid(k) = acc_isomesh % ParallelInfo % GlobalDOFs(i)
        xyz(3 * k - 2) = acc_isomesh % Nodes % x(i)
        xyz(3 * k - 1) = acc_isomesh % Nodes % y(i)
        xyz(3 * k) = acc_isomesh % Nodes % z(i)
      end do
    end if
    allocate(nrecv(pe), displ(pe))
    call MPI_ALLGATHER(nown, 1, MPI_INTEGER, nrecv, 1, MPI_INTEGER, &
        ELMER_COMM_WORLD, ierr2)
    ntot = sum(nrecv)
    displ(1) = 0
    do p = 2, pe
      displ(p) = displ(p - 1) + nrecv(p - 1)
    end do
    allocate(gidall(max(ntot, 1)))
    call MPI_GATHERV(gid, nown, MPI_INTEGER, gidall, nrecv, displ, &
        MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr2)
    allocate(xyzall(max(3 * ntot, 1)))
    call MPI_GATHERV(xyz, 3 * nown, MPI_DOUBLE_PRECISION, xyzall, 3 * nrecv, &
        3 * displ, MPI_DOUBLE_PRECISION, 0, ELMER_COMM_WORLD, ierr2)
    ! The segments, each held by one partition only.
    nloc = 0
    do i = 1, acc_isomesh % NumberOfBulkElements
      if (acc_isomesh % Elements(i) % PartIndex == ParEnv % MyPe) nloc = nloc + 1
    end do
    allocate(segbuf(max(2 * nloc, 1)))
    k = 0
    do i = 1, acc_isomesh % NumberOfBulkElements
      elem => acc_isomesh % Elements(i)
      if (elem % PartIndex /= ParEnv % MyPe) cycle
      if (elem % TYPE % ElementCode /= 202) then
        call Fatal('SaveGmshGeo2D', 'Only elements of type 202 can be saved in Gmsh geo format!')
      end if
      k = k + 1
      segbuf(2 * k - 1) = acc_isomesh % ParallelInfo % GlobalDOFs(elem % NodeIndexes(1))
      segbuf(2 * k) = acc_isomesh % ParallelInfo % GlobalDOFs(elem % NodeIndexes(2))
    end do
    call MPI_ALLGATHER(nloc, 1, MPI_INTEGER, nrecv, 1, MPI_INTEGER, &
        ELMER_COMM_WORLD, ierr2)
    acc_ngeolines = sum(nrecv)
    displ(1) = 0
    do p = 2, pe
      displ(p) = displ(p - 1) + nrecv(p - 1)
    end do
    allocate(segall(max(2 * acc_ngeolines, 1)))
    call MPI_GATHERV(segbuf, 2 * nloc, MPI_INTEGER, segall, 2 * nrecv, 2 * displ, &
        MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr2)
    if (ParEnv % MyPe == 0) then
      acc_ngeonodes = 0
      do i = 1, ntot
        acc_ngeonodes = max(acc_ngeonodes, gidall(i))
      end do
      allocate(acc_geox(max(acc_ngeonodes, 1)), acc_geoy(max(acc_ngeonodes, 1)), &
          acc_geoz(max(acc_ngeonodes, 1)), acc_geoline(2, max(acc_ngeolines, 1)))
      acc_geox = 0.0_dp; acc_geoy = 0.0_dp; acc_geoz = 0.0_dp
      do i = 1, ntot
        j = gidall(i)
        acc_geox(j) = xyzall(3 * i - 2)
        acc_geoy(j) = xyzall(3 * i - 1)
        acc_geoz(j) = xyzall(3 * i)
      end do
      do i = 1, acc_ngeolines
        acc_geoline(1, i) = segall(2 * i - 1)
        acc_geoline(2, i) = segall(2 * i)
      end do
      if (ntot /= acc_ngeonodes) then
        call Warn('SaveGmshGeo2D', 'Gathered ' // i2s(ntot) // ' nodes but the ' // &
            'numbers run to ' // i2s(acc_ngeonodes) // '!')
      end if
    end if
    deallocate(gid, xyz, gidall, xyzall, segbuf, segall, nrecv, displ)
  end subroutine acc_gathergeolines

  subroutine acc_savegmshgeo2d(mesh)
    implicit none
    type(mesh_t) :: mesh
    character(len=MAX_NAME_LEN) :: filename
    logical :: found, savenode
    logical, allocatable :: nodeused(:)
    integer :: geounit, iostat
    integer :: i1, i2, i, j, usednodes
    integer, allocatable :: neighbours(:,:), nodeorder(:)
    real(kind=dp) :: dx, dy, dz, meshdiam, meshparam, coord(3), nodeeps, &
        prevcoord(3), dist
    integer, parameter :: maxloops = 20
    integer :: loopoffset(maxloops), loopsize(maxloops), noloop, saveloops, nodeindex
    integer :: n, l, k
    logical :: newloop
    ! Unlike the other formats a geo file is one object: everything is
    ! collected on the master first.
    call acc_gathergeolines()
    if (ParEnv % MyPe /= 0) return
    if (acc_ngeonodes == 0) return
    filename = ListGetString(acc_params, 'Geo Filename', found)
    if (.not. found) filename = 'mesh.geo'
    open (newunit=geounit, file=filename, status='UNKNOWN', iostat=iostat)
    if (iostat /= 0) then
      call Warn('SaveGmshGeo2D', 'Could not open file: ' // trim(filename))
      return
    end if
    n = acc_ngeonodes
    allocate(neighbours(n, 2), nodeused(n), nodeorder(n))
    neighbours = 0
    nodeused = .false.
    nodeorder = 0
    ! Create a list of neighbours.
    ! Each node should have exactly two neighbours.
    do i = 1, acc_ngeolines
      i1 = acc_geoline(1, i)
      i2 = acc_geoline(2, i)
      if (i1 == 0 .or. i2 == 0) then
        call Warn('SaveGmshGeo2D', 'Invalid indexes: ' // i2s(i1) // ' and ' // i2s(i2))
      end if
      if (neighbours(i1, 1) == 0) then
        neighbours(i1, 1) = i2
      else
        neighbours(i1, 2) = i2
      end if
      if (neighbours(i2, 1) == 0) then
        neighbours(i2, 1) = i1
      else
        neighbours(i2, 2) = i1
      end if
    end do
    if (any(neighbours(:, 2) == 0)) then
      call Warn('SaveGmshGeo2D', 'This does not seem to be a closed loop!')
      close (geounit)
      return
    end if
    ! Compute the characteristic size of the bounding box.
    dx = maxval(acc_geox) - minval(acc_geox)
    dy = maxval(acc_geoy) - minval(acc_geoy)
    dz = maxval(acc_geoz) - minval(acc_geoz)
    meshdiam = max(dx, max(dy, dz))
    ! Currently defines a constant mesh parameter!
    meshparam = ListGetCReal(acc_params, 'Mesh Parameter', found)
    if (.not. found) meshparam = meshdiam / 50
    nodeeps = ListGetCReal(acc_params, 'Mesh Node Epsilon', found)
    if (.not. found) nodeeps = 1.0e-3 * meshdiam
    ! Find the continuous loops and neglect points that are redundant.
    j = 1
    usednodes = 0
    noloop = 1
    loopoffset = 0
    loopsize = 0
    newloop = .true.
    do i = 1, acc_ngeonodes
      ! Mark the current node before looking for the next one.
      nodeused(j) = .true.
      coord(1) = acc_geox(j)
      coord(2) = acc_geoy(j)
      coord(3) = acc_geoz(j)
      ! Check that the nodes are not so close to each other that Gmsh
      ! cannot handle the spline.
      savenode = .true.
      if (.not. newloop) then
        dist = sum((coord - prevcoord)**2)
        if (dist < nodeeps**2) then
          write (Message, '(A,ES12.4)') 'Skipping node ' // i2s(j) // ' too close: ', &
              sqrt(dist)
          call Info('IsosurfaceSolver', Message, level=12)
          savenode = .false.
        end if
      end if
      if (savenode) then
        loopsize(noloop) = loopsize(noloop) + 1
        usednodes = usednodes + 1
        nodeorder(usednodes) = j
        prevcoord = coord
      end if
      newloop = .false.
      ! Now find the next candidate node. This assumes closed loop!
      if (.not. nodeused(neighbours(j, 1))) then
        j = neighbours(j, 1)
      else if (.not. nodeused(neighbours(j, 2))) then
        j = neighbours(j, 2)
      else
        found = .false.
        do j = 1, acc_ngeonodes
          if (.not. nodeused(j)) then
            call Info('IsosurfaceSolver', 'Found a new start at node: ' // i2s(j), level=10)
            found = .true.
            exit
          end if
        end do
        if (found) then
          if (noloop == maxloops) then
            call Warn('IsosurfaceSolver', 'The static max number of loops (' // &
                i2s(maxloops) // ') is too small, ignoring the rest!')
            exit
          end if
          noloop = noloop + 1
          newloop = .true.
          loopoffset(noloop) = usednodes
        else
          call Info('IsosurfaceSolver', 'Could not find a new start, all nodes checked', level=10)
          exit
        end if
      end if
    end do
    if (usednodes == acc_ngeonodes) then
      call Info('IsosurfaceSolver', 'All nodes used')
    else
      call Info('IsosurfaceSolver', i2s(usednodes) // ' nodes used')
    end if
    call Info('IsoSurfaceSolver', 'Number of loops found: ' // i2s(noloop))
    call Info('IsosurfaceSolver', 'Writing points in geo file', level=10)
    saveloops = ListGetInteger(acc_params, 'Save Number Of Loops', found)
    if (.not. found) saveloops = noloop
    ! Asking for more than exist would just re-save the last one.
    if (saveloops > noloop) then
      call Info('IsosurfaceSolver', 'Only ' // i2s(noloop) // ' loops to save', level=8)
      saveloops = noloop
    end if
    nodeindex = 0
    do l = 1, saveloops
      ! Find the biggest unsaved loop.
      n = 0
      do k = 1, noloop
        if (loopsize(k) > n) then
          j = k
          n = loopsize(k)
        end if
      end do
      call Info('IsosurfaceSolver', 'Size of spline: ' // i2s(loopsize(j)), level=10)
      do i = 1, loopsize(j)
        nodeindex = nodeindex + 1
        k = nodeorder(loopoffset(j) + i)
        coord(1) = acc_geox(k)
        coord(2) = acc_geoy(k)
        coord(3) = acc_geoz(k)
        write (geounit, '(A,I0,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4,A)') 'Point(', &
            nodeindex, ') = {', coord(1), ',', coord(2), ',', coord(3), ',', &
            meshparam, '};'
      end do
      ! Ok, this was saved, erase the list.
      loopsize(j) = -loopsize(j)
    end do
    ! Revert from the negative values, to be able to use the logic again.
    loopsize = abs(loopsize)
    call Info('IsosurfaceSolver', 'Writing spline in geo file', level=10)
    nodeindex = 0
    do l = 1, saveloops
      ! Find the biggest unsaved loop.
      n = 0
      do k = 1, noloop
        if (loopsize(k) > n) then
          j = k
          n = loopsize(k)
        end if
      end do
      call Info('IsosurfaceSolver', 'Size of spline: ' // i2s(loopsize(j)), level=10)
      write (geounit, '(A)', advance='NO') 'Spline(' // i2s(l) // ') = {'
      do i = 1, loopsize(j)
        write (geounit, '(I0,A)', advance='NO') nodeindex + i, ', '
      end do
      write (geounit, '(I0,A)') nodeindex + 1, '};'
      nodeindex = nodeindex + loopsize(j)
      ! Ok, this was saved, erase the list.
      loopsize(j) = -loopsize(j)
    end do
    do j = 1, saveloops
      k = saveloops + 2 * (j - 1)
      write (geounit, '(A)') 'Line Loop(' // i2s(k + 1) // ') = {' // i2s(j) // '};'
      write (geounit, '(A)') 'Plane Surface(' // i2s(k + 2) // ') = {' // i2s(k + 1) // '};'
    end do
    close (geounit)
    call Info('IsosurfaceSolver', 'Wrote geo file: ' // trim(filename), level=6)
  end subroutine acc_savegmshgeo2d

end module isosurface_acc
