! T021-B OpenFoam2ElmerIO port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/OpenFoam2ElmerIO.F90, compiled against
! the built .mod files of the target configure (ADR-057 §3, T020 Poisson
! pattern).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls. Wrapper module
! state mirrors the original's SAVE locals (OFMesh/NoDir/Visited entry
! SAVE, DiffusionAssembly SAVE Nodes, DataAssembly SAVE ElementNodes)
! with identical lifetimes.
!
! File I/O and order-sensitive reductions stay in Fortran by construction:
! OpenFOAMBlocks / CreateFOAMMesh / ReadFOAMField are transcribed whole
! (same compiler, same call sequence, hence the same bytes including the
! original's defined-path behaviour), MINVAL/MAXVAL/SUM statistics and
! the ES12.5 Message writes are emitted by a wrapper, and the element
! update statements (SUM/MATMUL forms) are transcribed verbatim.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text beyond the transcribed call sequences.
module of2e_acc
  use iso_c_binding
  use DefUtils
  use ParallelUtils
  use ParticleUtils
  use Messages
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(variable_t), pointer, save :: acc_var => null()
  type(mesh_t), pointer, save :: acc_mesh => null()
  type(valuelist_t), pointer, save :: acc_params => null()
  type(element_t), pointer, save :: acc_elem => null()
  ! DataAssembly SAVE ElementNodes.
  type(nodes_t), target, save :: acc_nodes_data
  ! DiffusionAssembly SAVE Nodes.
  type(nodes_t), target, save :: acc_nodes_diff
  ! Entry SAVE OFMesh (+ per-call OFField) and OpenFOAMBlocks SAVE NoDir.
  type(mesh_t), pointer, save :: acc_ofmesh => null()
  real(dp), pointer, save :: acc_offield(:) => null()
  real(dp), pointer, save :: acc_weight(:) => null()
  integer, save :: acc_nodir = 0
  ! DiffusionAssembly SAVE DiffusivityName (with the original default).
  character(len=MAX_NAME_LEN), save :: acc_diffname = 'Heat Conductivity'
  ! Entry strings carried across wrappers.
  character(len=MAX_NAME_LEN), save :: acc_basedir = ''
  character(len=MAX_NAME_LEN), save :: acc_dirname = ''
  character(len=MAX_NAME_LEN), save :: acc_filename = ''
  integer, save :: acc_passive = 0
  ! Shared Gauss-point state (consumed before the next fill, as in the
  ! original where IP is a per-call local).
  type(gaussintegrationpoints_t), save :: acc_ip

contains

  ! -- entry state ----------------------------------------------------
  subroutine ofa_begin(h_solver) bind(c, name="ofa_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
    acc_params => GetSolverParams()
    acc_mesh => GetMesh()
    nullify(acc_var)
  end subroutine ofa_begin

  subroutine ofa_info_head() bind(c, name="ofa_info_head")
    call Info('OpenFoam2ElmerFit','-----------------------------------------', Level=4 )
    call Info('OpenFoam2ElmerFit','Projecting OpenFOAM data to Elmer mesh',Level=4)
  end subroutine ofa_info_head

  subroutine ofa_var_check() bind(c, name="ofa_var_check")
    acc_var => acc_solver % Variable
    if (.not. associated(acc_var)) then
      call Fatal('OpenFoam2ElmerFit','Variable must exist for the solver!')
    end if
  end subroutine ofa_var_check

  function ofa_rhs_len() result(n) bind(c, name="ofa_rhs_len")
    integer(c_int) :: n
    type(matrix_t), pointer :: amat
    amat => acc_solver % Matrix
    n = size(amat % Rhs)
  end function ofa_rhs_len

  ! -- entry SAVE mesh handling (mirrors SAVE OFMesh/NoDir/Visited) ----
  subroutine ofa_release_saved() bind(c, name="ofa_release_saved")
    call ReleaseMesh(acc_ofmesh)
    deallocate(acc_ofmesh)
  end subroutine ofa_release_saved

  subroutine ofa_alloc_mesh() bind(c, name="ofa_alloc_mesh")
    acc_ofmesh => AllocateMesh()
  end subroutine ofa_alloc_mesh

  ! -- scalar parameter queries ----------------------------------------
  subroutine ofa_passive(found) bind(c, name="ofa_passive")
    integer(c_int) :: found
    logical :: f
    acc_passive = ListGetInteger(acc_solver % Values,'Passive OpenFOAM Coordinate',f)
    found = merge(1, 0, f)
  end subroutine ofa_passive

  function ofa_passive_val() result(v) bind(c, name="ofa_passive_val")
    integer(c_int) :: v
    v = acc_passive
  end function ofa_passive_val

  function ofa_cylsymm() result(r) bind(c, name="ofa_cylsymm")
    integer(c_int) :: r
    logical :: c
    c = (CurrentCoordinateSystem() == CylindricSymmetric &
        .or. CurrentCoordinateSystem() == AxisSymmetric)
    r = merge(1, 0, c)
  end function ofa_cylsymm

  function ofa_meshdim() result(d) bind(c, name="ofa_meshdim")
    integer(c_int) :: d
    d = acc_mesh % MeshDim
  end function ofa_meshdim

  subroutine ofa_warn_reduced() bind(c, name="ofa_warn_reduced")
    call Warn('OpenFOAM2Elmer','Dimension of Elmer mesh is reduced, and OpenFOAM not?!')
  end subroutine ofa_warn_reduced

  function ofa_basedir(found) result(r) bind(c, name="ofa_basedir")
    integer(c_int) :: found
    integer(c_int) :: r
    logical :: f
    acc_basedir = GetString(acc_params,'OpenFOAM Directory',f)
    found = merge(1, 0, f)
    if (f) then
      call Info('OpenFoam2ElmerFit','Using given > OpenFOAM Directory < : '//trim(acc_basedir),Level=6)
    else
      call Fatal('OpenFoam2ElmerFit','> OpenFOAM Directory < must exist for the solver!')
    end if
    r = 0
  end function ofa_basedir

  function ofa_have_mesh1() result(r) bind(c, name="ofa_have_mesh1")
    integer(c_int) :: r
    r = merge(1, 0, ListCheckPresent(acc_params,'OpenFOAM Mesh 1'))
  end function ofa_have_mesh1

  ! -- OpenFOAMBlocks transcribed whole (file I/O + SystemCommand) -----
  subroutine ofa_blocks() bind(c, name="ofa_blocks")
    character(len=MAX_NAME_LEN) :: DirCommand
    logical :: FileExists
    integer, parameter :: InFileUnit = 28
    integer :: IOStatus
    acc_nodir = 0
    inquire(file=trim(acc_basedir), exist=FileExists)
    if (.not. FileExists) then
      call Fatal('OpenFoam2ElmerFit','OpenFOAM directory does not exist: '//trim(acc_basedir))
    end if
    acc_dirname = trim(acc_basedir)//'/0/'
    inquire(file=acc_dirname, exist=FileExists)
    if (.not. FileExists) then
      call Fatal('OpenFoam2ElmerFit','OpenFOAM mesh does not exist: '//trim(acc_dirname))
    end if
    acc_filename = trim(acc_dirname)//'C'
    call Info('OpenFoam2ElmerFit','Inquire file: '//trim(acc_filename),Level=12)
    inquire(file=acc_filename, exist=FileExists)
    if (FileExists) then
      call Info('OpenFoam2ElmerFit','Using OpenFOAM centers in: '//trim(acc_filename),Level=10)
      call ListAddString(acc_params, 'OpenFOAM Mesh 1', acc_dirname, .false.)
      acc_nodir = 1
      return
    end if
    DirCommand = 'ls -d '//trim(acc_dirname)//'*/ > OpenFOAMBlocks.txt'
    call Info('OpenFoam2ElmerFit','Performing command: '//trim(DirCommand),Level=12)
    call SystemCommand(DirCommand)
    open(InFileUnit,file='OpenFOAMBlocks.txt',iostat=IOstatus)
    if (IOStatus /= 0) then
      call Fatal('OpenFoam2ElmerFit','Could not open file: OpenFOAMBlocks.txt')
    end if
    do
      read(InFileUnit,'(A)',iostat=IOStatus) acc_dirname
      if (IOStatus /= 0) exit
      acc_filename = trim(acc_dirname)//'C'
      call Info('OpenFoam2ElmerFit','Inquire file: '//trim(acc_filename),Level=12)
      inquire(file=acc_filename, exist=FileExists)
      if (FileExists) then
        acc_nodir = acc_nodir + 1
        call Info('OpenFoam2ElmerFit','Using OpenFOAM centers in: '//trim(acc_filename),Level=10)
        call ListAddString(acc_params, 'OpenFOAM Mesh '//I2S(acc_nodir), acc_dirname, .false.)
      else
        call Info('OpenFoam2ElmerFit','No OpenFOAM center in: '//trim(acc_dirname),Level=12)
      end if
    end do
    close(InFileUnit)
    if (acc_nodir == 0) then
      call Fatal('OpenFoam2ElmerFit','No OpenFOAM mesh blocks found!')
    else
      call Info('OpenFoam2ElmerFit','Number of OpenFOAM blocks: '//I2S(acc_nodir),Level=10)
    end if
  end subroutine ofa_blocks

  subroutine ofa_parallel_reduce() bind(c, name="ofa_parallel_reduce")
    acc_nodir = ParallelReduction(acc_nodir)
  end subroutine ofa_parallel_reduce

  function ofa_nodir() result(n) bind(c, name="ofa_nodir")
    integer(c_int) :: n
    n = acc_nodir
  end function ofa_nodir

  subroutine ofa_info_nblocks() bind(c, name="ofa_info_nblocks")
    call Info('OpenFOAM2ElmerFit','Number of active OpenFOAM blocks: '//I2S(acc_nodir),Level=5)
  end subroutine ofa_info_nblocks

  ! -- assembly driver --------------------------------------------------
  subroutine ofa_default_initialize() bind(c, name="ofa_default_initialize")
    call DefaultInitialize()
  end subroutine ofa_default_initialize

  subroutine ofa_weight_alloc(n) bind(c, name="ofa_weight_alloc")
    integer(c_int), value :: n
    allocate(acc_weight(n))
  end subroutine ofa_weight_alloc

  subroutine ofa_weight_zero() bind(c, name="ofa_weight_zero")
    acc_weight = 0.0_dp
  end subroutine ofa_weight_zero

  subroutine ofa_weight_free() bind(c, name="ofa_weight_free")
    if (associated(acc_weight)) deallocate(acc_weight)
  end subroutine ofa_weight_free

  ! -- per-block file handling ------------------------------------------
  subroutine ofa_mesh_dir(i, found) bind(c, name="ofa_mesh_dir")
    integer(c_int), value :: i
    integer(c_int) :: found
    logical :: f
    acc_dirname = ListGetString(acc_params,'OpenFOAM Mesh '//I2S(i),f)
    found = merge(1, 0, f)
    if (.not. f) call Fatal('OpenFoam2ElmerFit','Could not find keyword: '//trim(acc_dirname))
  end subroutine ofa_mesh_dir

  subroutine ofa_info_block_head(i) bind(c, name="ofa_info_block_head")
    integer(c_int), value :: i
    call Info('OpenFOAM2ElmerFit','****************************',Level=8)
    call Info('OpenFOAM2ElmerFit','Treating OpenFOAM block: '//I2S(i),Level=5)
  end subroutine ofa_info_block_head

  subroutine ofa_info_reading() bind(c, name="ofa_info_reading")
    acc_filename = trim(acc_dirname)//'C'
    call Info('OpenFoam2ElmerFit','Reading OpenFOAM cell center from file: '//trim(acc_filename),Level=5)
  end subroutine ofa_info_reading

  ! -- CreateFOAMMesh transcribed whole ---------------------------------
  subroutine ofa_create_mesh() bind(c, name="ofa_create_mesh")
    integer :: line,i,j,k,n
    real(kind=dp) :: x,y,z
    integer :: NumberOfNodes, IOStatus
    integer, parameter :: InFileUnit = 28
    character(len=:), allocatable :: ReadStr
    logical :: f
    allocate(acc_ofmesh % Nodes)
    acc_ofmesh % NumberOfBulkElements = 0
    acc_ofmesh % NumberOfBoundaryElements = 0
    allocate(character(MAX_STRING_LEN)::ReadStr)
    open(InFileUnit,file=acc_filename, status='old', iostat=IOstatus)
    if (IOStatus /= 0) then
      call Fatal('OpenFoam2ElmerFit','Could not open file for reading: '//trim(acc_filename))
    end if
    call Info('OpenFoam2ElmerFit','Reading data points from file: '//trim(acc_filename),Level=6)
    j = 0
    do line = 1, 100
      read(InFileUnit,'(A)',iostat=IOStatus) ReadStr
      if (IOStatus /= 0) then
        call Warn('OpenFoam2ElmerFit','End of file after '//I2S(line)//' lines')
        exit
      end if
      j = index(ReadStr,'internalField',.true.)
      if (j > 0) exit
    end do
    if (j == 0) then
      call Warn('OpenFoam2ElmerFit','Could not find > internalField < in header!')
    else
      call Info('OpenFoam2ElmerFit','internalField found at line: '//I2S(line),Level=7)
    end if
    j = index(ReadStr,'nonuniform',.true.)
    if (j == 0) then
      call Warn('OpenFoam2ElmerFit','This routine only knows how to read nonuniform lists!')
      return
    end if
    read(InFileUnit,*,iostat=IOStatus) NumberOfNodes
    if (IOStatus /= 0) then
      call Fatal('OpenFoam2ElmerFit','Could not read number of nodes!')
    end if
    call Info('OpenFoam2ElmerFit','Number of OpenFOAM nodes: '//I2S(NumberOfNodes))
    i = ListGetInteger(acc_params,'Number of cells',f)
    if (i > 0 .and. i < NumberOfNodes) then
      NumberOfNodes = i
      call Info('OpenFoam2ElmerFit','Limiting number of OpenFOAM nodes: '//I2S(NumberOfNodes))
    end if
    n = NumberOfNodes
    allocate(acc_ofmesh % Nodes % x(n), acc_ofmesh % Nodes % y(n), acc_ofmesh % Nodes % z(n))
    acc_ofmesh % Nodes % x(1:n) = 0.0_dp
    acc_ofmesh % Nodes % y(1:n) = 0.0_dp
    acc_ofmesh % Nodes % z(1:n) = 0.0_dp
    acc_ofmesh % NumberOfNodes = NumberOfNodes
    read(InFileUnit,'(A)',iostat=IOStatus) ReadStr
    do i=1,n
      read(InFileUnit,'(A)',iostat=IOStatus) ReadStr
      if (IOStatus /= 0) then
        call Fatal('OpenFoam2ElmerFit','Could not read coordinate line: '//I2S(i))
      end if
      j = index(ReadStr,'(',.true.)
      if (j == 0) then
        call Fatal('OpenFoam2ElmerFit','Expecting a parenthesis at the start of OpenFOAM line: '//I2S(i))
      end if
      k = index(ReadStr,')',.true.)
      if (k == 0) then
        call Fatal('OpenFoam2ElmerFit','Expecting a parenthesis at the end of OpenFOAM line: '//I2S(i))
      end if
      read(ReadStr(j+1:k-1),*,iostat=IOStatus) x,y,z
      if (IOStatus /= 0) then
        call Fatal('OpenFoam2ElmerFit','Could not read coordinate values: '//I2S(i))
      end if
      acc_ofmesh % Nodes % x(i) = x
      acc_ofmesh % Nodes % y(i) = y
      acc_ofmesh % Nodes % z(i) = z
    end do
    close(InFileUnit)
    call Info('OpenFoam2ElmerFit','Created temporal OpenFOAM mesh just for nodes',Level=8)
    if (acc_passive == 1) then
      acc_ofmesh % Nodes % x = 0.0_dp
    else if (acc_passive == 2) then
      acc_ofmesh % Nodes % y = 0.0_dp
    else if (acc_passive == 3) then
      acc_ofmesh % Nodes % z = 0.0_dp
    end if
  end subroutine ofa_create_mesh

  subroutine ofa_field_alloc() bind(c, name="ofa_field_alloc")
    allocate(acc_offield(acc_ofmesh % NumberOfNodes))
  end subroutine ofa_field_alloc

  subroutine ofa_field_free() bind(c, name="ofa_field_free")
    if (associated(acc_offield)) deallocate(acc_offield)
  end subroutine ofa_field_free

  subroutine ofa_release_contents() bind(c, name="ofa_release_contents")
    call ReleaseMesh(acc_ofmesh)
  end subroutine ofa_release_contents

  ! -- ReadFOAMField transcribed whole ----------------------------------
  subroutine ofa_read_field(got) bind(c, name="ofa_read_field")
    integer(c_int) :: got
    character(len=MAX_NAME_LEN) :: TFileName, TSuffix
    integer :: line,i,j,k,n,nstep
    real(kind=dp) :: val
    integer :: NumberOfNodes, IOStatus
    integer, parameter :: InFileUnit = 28
    character(len=:), allocatable :: ReadStr
    logical :: f
    got = 0
    allocate(character(MAX_STRING_LEN)::ReadStr)
    TSuffix = ListGetString(acc_params,'OpenFOAM field',f)
    if (.not. f) then
      call Fatal('OpenFoam2ElmerFit','Give OpenFOAM field name!')
    end if
    nstep = ListGetInteger(acc_params,'OpenFOAM Timestep',f)
    if (f) call Info('OpenFoam2ElmerFit','Replacing 0 with timestep '//I2S(nstep))
    j = index(acc_filename,'/0/')
    k = len_trim(acc_filename)
    write(TFileName,'(A)') acc_filename(1:j)//I2S(nstep)//acc_filename(j+2:k-1)//trim(TSuffix)
    open(InFileUnit,file=TFileName, status='old', iostat=IOstatus)
    if (IOStatus /= 0) then
      call Warn('OpenFoam2ElmerFit','Could not open file for reading: '//trim(TFileName))
      return
    end if
    call Info('OpenFoam2ElmerFit','Reading data field from file: '//trim(TFileName),Level=6)
    j = 0
    do line = 1, 100
      read(InFileUnit,'(A)',iostat=IOStatus) ReadStr
      if (IOStatus /= 0) then
        call Warn('OpenFoam2ElmerFit','End of file after '//I2S(line)//' lines')
        goto 10
      end if
      j = index(ReadStr,'internalField',.true.)
      if (j > 0) then
        k = index(ReadStr,'<vector>',.true.)
        if (k > 0) then
          call Fatal('OpenFoam2ElmerFit','Currently implemented only for scalar fields!')
        end if
        exit
      end if
    end do
    if (j == 0) then
      call Warn('OpenFoam2ElmerFit','Could not find > internalField < in header!')
      goto 10
    else
      call Info('OpenFoam2ElmerFit','internalField found at line: '//I2S(line),Level=7)
    end if
    j = index(ReadStr,'nonuniform',.true.)
    if (j == 0) then
      call Warn('OpenFoam2ElmerFit','This routine only knows how to read nonuniform lists!')
      goto 10
    end if
    read(InFileUnit,*,iostat=IOStatus) n
    if (IOStatus /= 0) then
      call Fatal('OpenFoam2ElmerFit','Could not read number of nodes!')
    end if
    call Info('OpenFoam2ElmerFit','Number of OpenFOAM nodes: '//I2S(n),Level=10)
    read(InFileUnit,'(A)',iostat=IOStatus) ReadStr
    do i=1,acc_ofmesh % NumberOfNodes
      read(InFileUnit,'(A)',iostat=IOStatus) ReadStr
      if (IOStatus /= 0) then
        call Fatal('OpenFoam2ElmerFit','Could not read coordinate line: '//I2S(i))
      end if
      read(ReadStr,*,iostat=IOStatus) val
      acc_offield(i) = val
      if (IOStatus /= 0) then
        call Fatal('OpenFoam2ElmerFit','Could not read field values: '//I2S(i))
      end if
    end do
    call Info('OpenFoam2ElmerFit','Read data from OpenFOAM mesh region',Level=7)
    got = 1
10  close(InFileUnit)
  end subroutine ofa_read_field

  ! -- field statistics (order-sensitive reductions + ES12.5 writes) ----
  subroutine ofa_report_stats() bind(c, name="ofa_report_stats")
    real(kind=dp) :: MinF, MaxF, MeanF
    MinF = minval(acc_offield)
    MaxF = maxval(acc_offield)
    MeanF = sum(acc_offield) / size(acc_offield)
    write(Message,'(A,ES12.5)') 'Minimum field value: ',MinF
    call Info('OpenFoam2ElmerFit',Message,Level=6)
    write(Message,'(A,ES12.5)') 'Maximum field value: ',MaxF
    call Info('OpenFoam2ElmerFit',Message,Level=6)
    write(Message,'(A,ES12.5)') 'Average field value: ',MeanF
    call Info('OpenFoam2ElmerFit',Message,Level=6)
  end subroutine ofa_report_stats

  ! -- DataAssembly primitives ------------------------------------------
  function ofa_ofcount() result(n) bind(c, name="ofa_ofcount")
    integer(c_int) :: n
    n = acc_ofmesh % NumberOfNodes
  end function ofa_ofcount

  subroutine ofa_ofcoords(t, cyl, g) bind(c, name="ofa_ofcoords")
    integer(c_int), value :: t, cyl
    real(c_double) :: g(3)
    if (cyl /= 0) then
      g(1) = sqrt(acc_ofmesh % Nodes % x(t)**2 + acc_ofmesh % Nodes % y(t)**2)
      g(2) = acc_ofmesh % Nodes % z(t)
      g(3) = 0.0_dp
    else
      g(1) = acc_ofmesh % Nodes % x(t)
      g(2) = acc_ofmesh % Nodes % y(t)
      g(3) = acc_ofmesh % Nodes % z(t)
    end if
  end subroutine ofa_ofcoords

  function ofa_field_val(t) result(v) bind(c, name="ofa_field_val")
    integer(c_int), value :: t
    real(c_double) :: v
    v = acc_offield(t)
  end function ofa_field_val

  subroutine ofa_locate(g, l, idx) bind(c, name="ofa_locate")
    real(c_double) :: g(3), l(3)
    integer(c_int) :: idx
    call LocateParticleInMeshOctree(idx, g, l)
  end subroutine ofa_locate

  subroutine ofa_data_begin_elem(idx, n) bind(c, name="ofa_data_begin_elem")
    integer(c_int), value :: idx
    integer(c_int) :: n
    acc_elem => acc_mesh % Elements(idx)
    n = acc_elem % TYPE % NumberOfNodes
    call GetElementNodes(acc_nodes_data, acc_elem)
  end subroutine ofa_data_begin_elem

  subroutine ofa_data_eleminfo(u, v, w, n, basis, detj, stat) bind(c, name="ofa_data_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: n
    real(c_double) :: basis(n), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_elem, acc_nodes_data, u, v, w, detj, basis)
    stat = merge(1, 0, ok)
  end subroutine ofa_data_eleminfo

  subroutine ofa_data_scatter(val, n, basis) bind(c, name="ofa_data_scatter")
    real(c_double), value :: val
    integer(c_int), value :: n
    real(c_double) :: basis(n)
    type(matrix_t), pointer :: amat
    integer :: i, j
    real(kind=dp) :: weight
    amat => acc_solver % Matrix
    do i = 1, n
      j = acc_var % Perm(acc_elem % NodeIndexes(i))
      if (j == 0) cycle
      weight = basis(i)
      amat % Rhs(j) = amat % Rhs(j) + weight * val
      acc_weight(j) = acc_weight(j) + weight
    end do
  end subroutine ofa_data_scatter

  ! -- fit scaling (whole-array updates) --------------------------------
  function ofa_fit_coeff(found) result(c) bind(c, name="ofa_fit_coeff")
    integer(c_int) :: found
    real(c_double) :: c
    logical :: f
    c = ListGetCReal(acc_params,'Fit Coefficient',f)
    found = merge(1, 0, f)
  end function ofa_fit_coeff

  subroutine ofa_apply_fit(coeff) bind(c, name="ofa_apply_fit")
    real(c_double), value :: coeff
    type(matrix_t), pointer :: amat
    amat => acc_solver % Matrix
    amat % Rhs = coeff * amat % Rhs
    amat % Values(amat % Diag) = coeff * acc_weight
  end subroutine ofa_apply_fit

  ! -- DiffusionAssembly primitives -------------------------------------
  function ofa_max_elem_nodes() result(n) bind(c, name="ofa_max_elem_nodes")
    integer(c_int) :: n
    n = acc_solver % Mesh % MaxElementNodes
  end function ofa_max_elem_nodes

  function ofa_diff_coeff(glob) result(c) bind(c, name="ofa_diff_coeff")
    integer(c_int) :: glob
    real(c_double) :: c
    logical :: g, f
    c = GetCReal(acc_solver % Values,'Diffusion Coefficient',g)
    glob = merge(1, 0, g)
    acc_diffname = GetString(acc_solver % Values,'Diffusivity Name',f)
    if (.not. f) acc_diffname = 'Heat Conductivity'
  end function ofa_diff_coeff

  function ofa_nof_active() result(n) bind(c, name="ofa_nof_active")
    integer(c_int) :: n
    n = GetNOFActive()
  end function ofa_nof_active

  function ofa_coord_dim() result(d) bind(c, name="ofa_coord_dim")
    integer(c_int) :: d
    d = CoordinateSystemDimension()
  end function ofa_coord_dim

  subroutine ofa_begin_active(t, n) bind(c, name="ofa_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes()
  end subroutine ofa_begin_active

  subroutine ofa_diff_nodes() bind(c, name="ofa_diff_nodes")
    call GetElementNodes(acc_nodes_diff)
  end subroutine ofa_diff_nodes

  function ofa_gauss_begin() result(n) bind(c, name="ofa_gauss_begin")
    integer(c_int) :: n
    acc_ip = GaussPoints(acc_elem)
    n = acc_ip % n
  end function ofa_gauss_begin

  subroutine ofa_gauss_point(k, u, v, w, s) bind(c, name="ofa_gauss_point")
    integer(c_int), value :: k
    real(c_double) :: u, v, w, s
    u = acc_ip % u(k)
    v = acc_ip % v(k)
    w = acc_ip % w(k)
    s = acc_ip % s(k)
  end subroutine ofa_gauss_point

  subroutine ofa_diff_eleminfo(u, v, w, n, basis, dbasis, detj, stat) bind(c, name="ofa_diff_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: n
    real(c_double) :: basis(n), dbasis(n,3), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_elem, acc_nodes_diff, u, v, w, detj, basis, dbasis)
    stat = merge(1, 0, ok)
  end subroutine ofa_diff_eleminfo

  subroutine ofa_diff_fill(ld, n, diffdata, loc) bind(c, name="ofa_diff_fill")
    integer(c_int), value :: ld, n
    real(c_double) :: diffdata(3,3,ld)
    integer(c_int) :: loc
    type(valuelist_t), pointer :: material
    real(kind=dp), pointer :: hwrk(:,:,:) => null()
    integer, pointer :: indexes(:)
    integer :: i, j
    logical :: l
    material => GetMaterial()
    indexes => acc_elem % NodeIndexes
    call ListGetRealArray(material, acc_diffname, hwrk, n, indexes, l)
    loc = merge(1, 0, l)
    if (l) then
      diffdata = 0.0_dp
      if (size(hwrk,1) == 1) then
        do i=1,3
          diffdata(i,i,1:n) = hwrk(1,1,1:n)
        end do
      else if (size(hwrk,2) == 1) then
        do i=1,min(3,size(hwrk,1))
          diffdata(i,i,1:n) = hwrk(i,1,1:n)
        end do
      else
        do i=1,min(3,size(hwrk,1))
          do j=1,min(3,size(hwrk,2))
            diffdata(i,j,1:n) = hwrk(i,j,1:n)
          end do
        end do
      end if
    end if
  end subroutine ofa_diff_fill

  subroutine ofa_diff_ip(glob, loc, coeff, dim, n, ld, dbasis, basis, diffdata, stiff, force, sdet) &
      bind(c, name="ofa_diff_ip")
    integer(c_int), value :: glob, loc, dim, n, ld
    real(c_double), value :: coeff, sdet
    real(c_double) :: dbasis(n,3), basis(n), diffdata(3,3,ld)
    real(c_double) :: stiff(ld,ld), force(ld)
    real(kind=dp) :: DiffMatrix(3,3)
    integer :: i, j, p, q
    if (glob /= 0) then
      do i=1,n
        do j=1,n
          stiff(i,j) = stiff(i,j) + sdet * &
              coeff * sum(dbasis(i,1:dim) * dbasis(j,1:dim))
        end do
      end do
    else if (loc /= 0) then
      do p=1,dim
        do q=1,dim
          DiffMatrix(p,q) = sum(diffdata(p,q,1:n) * basis(1:n))
        end do
      end do
      do i=1,n
        do j=1,n
          stiff(i,j) = stiff(i,j) + sdet * &
              sum(matmul(DiffMatrix(1:dim,1:dim), dbasis(j,1:dim)) * dbasis(i,1:dim))
        end do
      end do
    end if
  end subroutine ofa_diff_ip

  subroutine ofa_default_update(ld, stiff, force) bind(c, name="ofa_default_update")
    integer(c_int), value :: ld
    real(c_double) :: stiff(ld,ld), force(ld)
    call DefaultUpdateEquations(stiff(1:ld,1:ld), force(1:ld))
  end subroutine ofa_default_update

  ! -- shared assembly tail ---------------------------------------------
  subroutine ofa_default_bulk() bind(c, name="ofa_default_bulk")
    call DefaultFinishBulkAssembly()
  end subroutine ofa_default_bulk

  subroutine ofa_default_finish() bind(c, name="ofa_default_finish")
    call DefaultFinishAssembly()
  end subroutine ofa_default_finish

  subroutine ofa_default_dirichlet() bind(c, name="ofa_default_dirichlet")
    call DefaultDirichletBCs()
  end subroutine ofa_default_dirichlet

  function ofa_default_solve() result(norm) bind(c, name="ofa_default_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function ofa_default_solve

  subroutine ofa_info_alldone() bind(c, name="ofa_info_alldone")
    call Info('OpenFoam2ElmerFit','All done', Level=4 )
    call Info('OpenFoam2ElmerFit','-----------------------------------------', Level=4 )
  end subroutine ofa_info_alldone

  subroutine ofa_fatal_alloc() bind(c, name="ofa_fatal_alloc")
    call Fatal('OpenFoam2ElmerFit','Allocation error in DiffusionAssembly!')
  end subroutine ofa_fatal_alloc

end module of2e_acc
