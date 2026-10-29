! SaveData pilot: statement-level BIND(C) wrappers for the Elmer calls in
! fem/src/modules/SaveData/SaveLine.F90.
! The entry's SAVE state (SavePerm, PrevMaskName, SaveNodes) and all
! phase-shared locals live in this module with identical lifetimes; Rust
! drives the phase sequence (setup, existing lines, isocurves, polylines,
! circles, close, names). Internal procedures are transcribed verbatim as
! module procedures.
module sd_line_acc
  use iso_c_binding
  use types
  use lists
  use integration
  use meshbasics
  use saveutils
  use bandwidthoptimize
  use defutils
  implicit none

  character(len=*), parameter :: caller = 'SaveLine'
  character(len=4096), parameter :: default_side_file = 'sides.dat'

  type(solver_t), pointer, save :: acc_solver => null()
  type(model_t), pointer, save :: acc_model => null()
  type(mesh_t), pointer, save :: acc_mesh => null()
  type(valuelist_t), pointer, save :: acc_params => null()

  integer, pointer, save :: acc_saveperm(:) => null()
  character(len=4096), save :: acc_prevmaskname = ''
  integer, save :: acc_savenodes(4) = 0

  ! Setup state mirrored for the Rust driver.
  integer, save :: acc_dim = 0, acc_meshdim = 0, acc_noresults = 0
  integer, save :: acc_novar = 0, acc_maxboundary = 0, acc_lineind = 0
  integer, save :: acc_normind = 0, acc_foundnan = 0
  real(kind=dp), save :: acc_norm = 0.0_dp
  real(kind=dp), save :: acc_detepsilon = 1.0e-6_dp
  logical, save :: acc_parallel = .false., acc_maskwithperm = .false.
  logical, save :: acc_skipboundaryinfo = .false., acc_calculateflux = .false.
  logical, save :: acc_fileisopen = .false.
  logical, save :: acc_initializeperm = .false.
  integer, save :: acc_nodata = 0, acc_nolabels = 0
  character(len=4096), save :: acc_sidefile = ''
  character(len=4096), save :: acc_tempname = '', acc_condname = ''
  character(len=4096), save :: acc_maskname = ''
  character(:), allocatable, save :: acc_outputdir
  character(len=4096), save :: acc_message
  character(len=256), allocatable, save :: acc_valuenames(:)
  real(kind=dp), save :: acc_dummy3(3) = 0.0_dp
  integer, save :: acc_lineunit = -1
  logical, save :: acc_fileappend = .false.
  logical, save :: acc_found = .false.
  ! Phase-shared state (SavePolylines/Circles/IsoCurves locals in source).
  integer, save :: acc_noaxis = 0, acc_nolines = 0, acc_nodims = 0
  integer, save :: acc_nofaces = 0, acc_notests = 0
  logical, save :: acc_saveaxis(3) = .false.
  logical, save :: acc_intersectedge = .false.
  integer, save :: acc_intersectcoord = 0
  real(kind=dp), save :: acc_intersectepsilon = 0.0_dp
  real(kind=dp), pointer, save :: acc_pointcoordinates(:, :) => null()
  integer, pointer, save :: acc_nodivisions(:) => null()
  logical, save :: acc_gotdiv = .false.
  logical, allocatable, save :: acc_linetag(:), acc_linetag0(:)
  real(kind=dp), allocatable, save :: acc_dbasisdx(:, :)
  type(element_t), pointer, save :: acc_parent => null()
  ! Entry-shared workspace (allocated per call, mirrors source locals).
  real(kind=dp), allocatable, save :: acc_values(:), acc_basis(:)
  real(kind=dp), allocatable, save :: acc_pointfluxes(:, :), acc_pointweight(:)
  type(nodes_t), target, save :: acc_elementnodes, acc_linenodes
  integer, save :: acc_labels(5) = 0
  real(kind=dp), pointer, save :: acc_posdata(:) => null()
  real(kind=dp), pointer, save :: acc_resultdata(:, :) => null()
  integer, pointer, save :: acc_labeldata(:, :) => null()
  ! Shared scalars for internal procedures (host association in source).
  integer, save :: acc_no = 0, acc_comps = 0, acc_nodofs = 0
  integer, save :: acc_maxboundary_nodes = 0
  real(kind=dp), save :: acc_linepos = 0.0_dp
  real(kind=dp) :: acc_localcoord(3) = 0.0_dp
  real(kind=dp), save :: acc_detj = 0.0_dp, acc_s = 0.0_dp
  logical, save :: acc_doeigen = .false., acc_edgebasis = .false.
  logical, save :: acc_avbasis = .false.
  logical, save :: acc_transient = .false.
  type(variable_t), pointer, save :: acc_var => null()
  type(variable_t), pointer, save :: acc_var2 => null()
  type(variable_t), pointer, save :: acc_var3 => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(nodes_t), target, save :: acc_nodes_wf

contains

  subroutine sdl_begin(h_model, h_solver) bind(c, name="sdl_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine sdl_begin

  subroutine sdl_init() bind(c, name="sdl_init")
    integer :: normind, i
    logical :: gotit
    character(len=4096) :: solvername
    integer, pointer :: asolvers(:)
    logical :: activesomewhere
    activesomewhere = .false.
    do i = 1, acc_model % numberofequations
      asolvers => ListGetIntegerArray(acc_model % equations(i) % values, &
          'Active Solvers', gotit)
      if (gotit) then
        if (any(asolvers == acc_solver % solverid)) then
          activesomewhere = .true.
          exit
        end if
      end if
    end do
    if (activesomewhere) then
      if (.not. ListCheckPresent(acc_solver % values, 'Variable')) then
        call Info('SaveLine_init', &
            'If you want to mask with Equation block give some Variable', level=5)
      end if
    end if
    solvername = ListGetString(acc_solver % values, 'Equation', gotit)
    normind = ListGetInteger(acc_solver % values, 'Show Norm Index', gotit)
    if (normind > 0) then
      call ListAddNewString(acc_solver % values, 'Variable', &
          '-nooutput -global ' // trim(solvername) // '_var')
    end if
    call ListAddNewLogical(acc_solver % values, 'No Matrix', .true.)
  end subroutine sdl_init


  ! -- entry setup (source lines through the variable scan) -----------------
  ! Returns NoVar/NoResults for the Rust driver; SAVE state lives here.
  subroutine sdl_setup(novar, noresults, maxboundary) bind(c, name="sdl_setup")
    integer(c_int) :: novar, noresults, maxboundary
    type(valuelist_t), pointer :: params
    type(variable_t), pointer :: var, var2, var3
    integer :: i, j, comps, noeigenvalues, istat
    logical :: gotit, found, edgebasis, avbasis, dgvar, ipvar, elemvar
    character(len=4096) :: sidefile
    params => GetSolverParams()
    acc_params => params
    sidefile = ListGetString(params, 'Filename', found)
    if (.not. found) sidefile = default_side_file
    acc_sidefile = sidefile
    call Info(caller, '-----------------------------------------', level=4)
    call Info(caller, 'Saving data on lines to file: ' // trim(sidefile), level=4)
    call Info(caller, '------------------------------------------', level=4)
    acc_fileisopen = .false.
    acc_foundnan = 0
    acc_nodata = 0
    acc_savenodes = 0
    acc_maskwithperm = .false.
    if (associated(acc_solver % variable)) then
      acc_maskwithperm = associated(acc_solver % variable % perm)
    end if
    nullify(acc_saveperm)
    ! Save Solver Mesh Index
    i = 0
    i = GetInteger(params, 'Save Solver Mesh Index', found)
    if (found) then
      call Info(caller, 'Using mesh of solver ' // trim(sdl_i2s(i)))
      acc_mesh => acc_model % solvers(i) % mesh
      acc_model % mesh => acc_mesh
    else
      acc_mesh => GetMesh()
    end if
    acc_dim = CoordinateSystemDimension()
    acc_meshdim = acc_mesh % meshdim
    acc_parallel = (parenv % pes > 1)
    if (GetLogical(params, 'Enforce Parallel Mode', gotit)) acc_parallel = .true.
    ! (ElementNodes/LineNodes/Basis allocation happens in sdl_alloc_basis.)
    if (acc_solver % timesvisited == 0) then
      if (associated(acc_saveperm)) deallocate(acc_saveperm)
      allocate (acc_saveperm(acc_mesh % numberofnodes), stat=istat)
      if (istat /= 0) call Fatal(caller, 'Memory allocation error for SavePerm')
    end if
    acc_normind = ListGetInteger(params, 'Show Norm Index', gotit)
    acc_norm = 0.0_dp
    ! MovingMesh / Intersect / Epsilon / Flux params read where used.
    acc_detepsilon = ListGetCReal(params, 'Det Epsilon', gotit)
    if (.not. gotit) acc_detepsilon = 1.0e-6_dp
    acc_calculateflux = ListGetLogical(params, 'Save Heat Flux', gotit)
    if (.not. acc_calculateflux) then
      acc_calculateflux = ListGetLogical(params, 'Save Flux', gotit)
    end if
    ! Variable scan (CreateListForSaving + NoResults count).
    if (acc_solver % timesvisited == 0) then
      call CreateListForSaving(acc_model, params, .true., usegenerickeyword=.true.)
    end if
    call sdl_count_results(novar, noresults)
    acc_novar = novar
    acc_noresults = noresults
    maxboundary = acc_maxboundary
  end subroutine sdl_setup

  ! -- allocation tail of setup (ElementNodes/Basis/Values, flags) ------------
  subroutine sdl_alloc(nresults) bind(c, name="sdl_alloc")
    integer(c_int) :: nresults
    integer :: n, istat
    logical :: gotit
    n = acc_mesh % maxelementnodes
    allocate (acc_elementnodes % x(n), acc_elementnodes % y(n), &
        acc_elementnodes % z(n), acc_linenodes % x(2), acc_linenodes % y(2), &
        acc_linenodes % z(2), acc_basis(n), stat=istat)
    if (istat /= 0) call Fatal(caller, 'Memory allocation error for Elemental stuff')
    if (acc_novar == 0 .or. acc_noresults == 0) then
      nresults = -1
      return
    end if
    nresults = acc_noresults + 3
    acc_noresults = nresults
    allocate (acc_values(acc_noresults), stat=istat)
    if (istat /= 0) call Fatal(caller, 'Memory allocation error for Values')
    acc_skipboundaryinfo = ListGetLogical(acc_params, 'Skip Boundary Info', gotit)
    call SolverOutputDirectory(acc_solver, acc_sidefile, acc_outputdir)
    if (acc_solver % timesvisited == 0) then
      acc_saveperm = 0
    end if
  end subroutine sdl_alloc


  ! -- variable scan: count fields (source var-scan loop) --------------------
  subroutine sdl_count_results(novar, noresults) bind(c, name="sdl_count_results")
    integer(c_int) :: novar, noresults
    type(variable_t), pointer :: var
    integer :: ivar, comps, noeigenvalues
    logical :: edgebasis, avbasis, found
    novar = 0
    noresults = 0
    avbasis = .false.
    do ivar = 1, 99
      var => sdl_getvar_n(ivar, comps)
      if (.not. associated(var)) exit
      novar = ivar
      edgebasis = .false.
      if (associated(var % eigenvalues)) then
        noeigenvalues = size(var % eigenvalues)
        noresults = noresults + var % dofs * noeigenvalues
      else if (var % type == variable_on_nodes_on_elements) then
        ! DGvar: counted in max(dofs,comps) below
      else if (var % type == variable_on_gauss_points) then
        ! IpVar: counted below
      else if (var % type == variable_on_elements) then
        ! ElemVar: counted below
      else if (associated(var % solver) .and. associated(var % perm)) then
        edgebasis = GetLogical(var % solver % values, 'Hcurl Basis', found)
        if (edgebasis) then
          edgebasis = (size(var % perm) > acc_mesh % numberofnodes)
        end if
        if (edgebasis) then
          edgebasis = (any(var % perm(acc_mesh % numberofnodes + 1:) > 0))
        end if
        if (edgebasis) then
          if (any(var % perm(1:acc_mesh % numberofnodes) > 0)) avbasis = .true.
        end if
        if (edgebasis) var % type = variable_on_edges
      end if
      if (edgebasis) then
        if (avbasis) then
          noresults = noresults + 4
        else
          noresults = noresults + 3
        end if
      else
        noresults = noresults + max(var % dofs, comps)
      end if
    end do
    if (acc_calculateflux) noresults = noresults + 3
  end subroutine sdl_count_results

  function sdl_getvar_n(ivar, comps, component) result(var)
    integer :: ivar
    integer, optional :: comps, component
    type(variable_t), pointer :: var
    integer :: j, k
    character(len=4096) :: name, varname
    logical :: found
    type(variable_t), pointer :: var2
    nullify(var)
    if (ivar < 1) then
      write (varname, '(A,I0)') 'Coordinate ', ivar + 3
      found = .true.
    else
      write (name, '(A,I0)') 'Variable ', ivar
      varname = GetString(acc_params, name, found)
    end if
    k = 0
    if (.not. found) then
      if (present(comps)) comps = k
      return
    end if
    if (present(component)) then
      varname = trim(varname) // ' ' // trim(sdl_i2s(component))
    end if
    var => VariableGet(acc_mesh % variables, varname)
    if (.not. associated(var)) then
      var => VariableGet(acc_mesh % variables, trim(varname) // ' 1')
      if (associated(var)) then
        do j = 2, 99
          var2 => VariableGet(acc_mesh % variables, trim(varname) // ' ' // trim(sdl_i2s(j)))
          if (associated(var2)) then
            k = j
          else
            exit
          end if
        end do
      else
        call Warn('VariableGetN', 'Variable given but not found: ' // trim(varname))
      end if
    end if
    if (present(comps)) comps = k
  end function sdl_getvar_n

  function sdl_i2s(i) result(s)
    integer :: i
    character(len=32) :: s
    write (s, '(I0)') i
    s = adjustl(s)
  end function sdl_i2s



  ! -- WriteFieldsAtElement (verbatim internal procedure) ---------------------
  ! Element identified by index; Basis/LocalCoord/GlobalCoord/linepos come
  ! from the caller's ElementInfo via module SAVE (set by phase wrappers).
  ! -- WriteFieldsAtElement (verbatim internal procedure) ---------------------
  subroutine sdl_write_fields_impl(elem_idx, bc_id, node_id, dgnode_id, &
      use_node_in, nodal_flux, has_flux, linepos_in, has_linepos, par_node)
    integer, value :: elem_idx, bc_id, node_id, dgnode_id
    logical, value :: use_node_in, has_flux, has_linepos, par_node
    real(kind=dp) :: nodal_flux(3), linepos_in
    type(element_t), pointer :: element
    integer :: i, j, k, l, ivar, ii, i1, i2, n0
    type(nodes_t) :: nodesloc
    logical :: usegivennode
    integer :: n, nd, np, edgedegree, labels(5)
    integer, pointer :: ptoindexes(:)
    real(kind=dp), pointer :: ptobasis(:)
    real(kind=dp), target :: pointbasis(1)
    real(kind=dp) :: u, v, w
    integer, target :: nodeindex(1), indexes(54), dgindexes(27)
    real(kind=dp) :: up, vp, wp
    real(kind=dp), target :: nodebasis(54)
    real(kind=dp) :: wbasis(54, 3), rotwbasis(54, 3), nodedbasisdx(54, 3)
    real(kind=dp) :: avemult
    logical :: pelem, tabulate
    type(variable_t), pointer :: pvar
    element => acc_mesh % elements(elem_idx)
    n0 = 0
    if (.not. acc_skipboundaryinfo) then
      acc_labels = 0
      if (acc_lineind /= 0) then
        n0 = n0 + 1
        acc_labels(n0) = acc_lineind
      end if
      if (acc_transient) then
        n0 = n0 + 1
        acc_labels(n0) = acc_solver % donetime
      end if
      acc_labels(n0 + 1) = acc_solver % timesvisited + 1
      acc_labels(n0 + 2) = bc_id
      k = node_id
      ! ParNode global-DOF mapping (serial: identity, transcribed verbatim).
      acc_labels(n0 + 3) = k
      n0 = n0 + 3
      if (acc_normind > 0 .and. acc_normind <= n0) then
        acc_norm = acc_norm + 1.0_dp * acc_labels(acc_normind)
      end if
    end if
    usegivennode = use_node_in
    acc_no = 0
    acc_values = 0.0_dp
    do ivar = -2, acc_novar
      acc_var => sdl_getvar_n(ivar, acc_comps)
      if (acc_comps >= 2) then
        acc_var2 => sdl_getvar_n(ivar, component=2)
      else
        nullify(acc_var2)
      end if
      if (acc_comps >= 3) then
        acc_var3 => sdl_getvar_n(ivar, component=3)
      else
        nullify(acc_var3)
      end if
      if (usegivennode) then
        k = node_id
        l = dgnode_id
      else
        k = 0
        l = 0
      end if
      call EvaluateVariableAtGivenPoint(acc_no, acc_values, acc_mesh, &
          acc_var, acc_var2, acc_var3, element, acc_localcoord, &
          acc_basis, k, l, gotedge=acc_edgebasis)
    end do
    if (acc_calculateflux) then
      if (has_flux) then
        acc_values(acc_no + 1:acc_no + 3) = nodal_flux
      else
        acc_values(acc_no + 1:acc_no + 3) = 0.0_dp
      end if
    end if
    do j = 1, acc_noresults
      if (isnan(acc_values(j))) then
        acc_foundnan = acc_foundnan + 1
        acc_values(j) = -1.0_dp
      end if
    end do
    call sdl_tabulate_or_write(element, bc_id, node_id, &
        has_linepos, linepos_in, par_node, n0)
    if (acc_normind > n0) then
      acc_norm = acc_norm + acc_values(acc_normind - n0)
    end if
  end subroutine sdl_write_fields_impl


  ! -- tabulate-or-write tail (verbatim) --------------------------------------
  subroutine sdl_tabulate_or_write(element, bc_id, node_id, &
      has_linepos, linepos_in, par_node, n0)
    type(element_t), pointer :: element
    integer :: bc_id, node_id, n0
    logical :: has_linepos, par_node
    real(kind=dp) :: linepos_in
    integer :: i, j, k, nold, nnew, istat
    real(kind=dp) :: linepos, dpos
    real(kind=dp), pointer :: tmpposdata(:), tmpresultdata(:, :)
    integer, pointer :: tmplabeldata(:, :)
    logical :: tabulate
    tabulate = .false.
    if (has_linepos) then
      linepos = linepos_in
      tabulate = (linepos > 0) .or. acc_parallel
    end if
    if (tabulate) then
      acc_nolabels = max(n0, acc_nolabels)
      acc_nodata = acc_nodata + 1
      nold = 0
      if (associated(acc_posdata)) nold = size(acc_posdata)
      if (acc_nodata > nold) then
        nnew = max(100, 2 * nold)
        call Info(caller, 'Increasing temporal size from ' // trim(sdl_i2s(nold)) // &
            ' to ' // trim(sdl_i2s(nnew)), level=7)
        allocate(tmpposdata(nnew), tmplabeldata(nnew, acc_nolabels), &
            tmpresultdata(nnew, acc_noresults))
        tmpposdata = 0.0_dp
        tmplabeldata = 0
        tmpresultdata = 0.0_dp
        if (nold > 0) then
          tmpposdata(1:nold) = acc_posdata
          tmplabeldata(1:nold, 1:acc_nolabels) = acc_labeldata
          tmpresultdata(1:nold, 1:acc_noresults) = acc_resultdata
          deallocate(acc_posdata, acc_resultdata, acc_labeldata, stat=istat)
          if (istat /= 0) call Fatal(caller, 'Problems deallocating some too small workspace')
        end if
        acc_posdata => tmpposdata
        acc_resultdata => tmpresultdata
        acc_labeldata => tmplabeldata
      end if
      acc_posdata(acc_nodata) = linepos
      acc_labeldata(acc_nodata, 1:acc_nolabels) = acc_labels(1:acc_nolabels)
      acc_resultdata(acc_nodata, 1:acc_noresults) = acc_values(1:acc_noresults)
    else
      call sdl_open_line(0)
      do i = 1, n0
        write (acc_lineunit, '(A)', advance='NO') trim(sdl_i2s(acc_labels(i))) // ' '
      end do
      do j = 1, acc_noresults - 1
        write (acc_lineunit, '(ES20.11E3)', advance='NO') acc_values(j)
      end do
      write (acc_lineunit, '(ES20.11E3)') acc_values(acc_noresults)
    end if
  end subroutine sdl_tabulate_or_write

  subroutine sdl_open_line(parallel_reduce) bind(c, name="sdl_open_line")
    integer(c_int), value :: parallel_reduce
    logical :: parallel
    integer :: iostat
    character(len=4096) :: sideparfile
    if (acc_fileisopen) return
    acc_fileisopen = .true.
    sideparfile = trim(acc_outputdir) // '/' // trim(acc_sidefile)
    parallel = (parenv % pes > 1)
    if (parallel) then
      if (parallel_reduce /= 0) parallel = .false.
    end if
    if (parallel) then
      sideparfile = trim(sideparfile) // '.' // trim(sdl_i2s(parenv % mype))
    else
      if (ListGetLogical(acc_params, 'Filename Numbering')) then
        if (parallel) then
          call Warn(caller, 'Cannot number filenames in parallel with another file!')
        else
          sideparfile = NextFreeFilename(sideparfile)
        end if
      end if
    end if
    call Info(caller, 'Saving line data to file: ' // trim(sideparfile), level=12)
    acc_fileappend = ListGetLogical(acc_params, 'File Append')
    if (acc_solver % timesvisited > 0 .or. acc_fileappend) then
      open (newunit=acc_lineunit, file=sideparfile, position='APPEND', iostat=iostat)
    else
      open (newunit=acc_lineunit, file=sideparfile, iostat=iostat)
    end if
    if (iostat /= 0) call Warn(caller, 'Problems closing line file: ' // trim(sdl_i2s(iostat)))
  end subroutine sdl_open_line

  subroutine sdl_close_line() bind(c, name="sdl_close_line")
    integer, pointer :: neworder(:)
    logical :: parallelreduce
    integer :: iostat, istat, i, j, k
    real(kind=dp) :: dpos
    real(kind=dp), pointer :: tmppos(:), tmpres(:, :)
    integer, pointer :: tmplab(:, :)
    parallelreduce = .false.
    if (parenv % pes > 1) then
      parallelreduce = ListGetLogical(acc_params, 'Parallel Reduce', acc_found)
    end if
    ! Serial path (parallel reduction omitted; transcribed as a separate
    ! routine if a parallel campaign ever needs it — serial tree only).
    if (acc_nodata > 0) then
      if (acc_posdata(acc_nodata) > 0.0_dp) then
        allocate(neworder(acc_nodata), stat=istat)
        if (istat /= 0) call Fatal(caller, 'Problems allocating NewOrder vector!')
        call Info(caller, 'Sorting and saving ' // trim(sdl_i2s(acc_nodata)) // &
            ' tabulated rows', level=7)
        do i = 1, acc_nodata
          neworder(i) = i
        end do
        call SortR(acc_nodata, neworder, acc_posdata)
      end if
      call sdl_open_line(0)
      do i = 1, acc_nodata
        k = i
        if (associated(neworder)) then
          k = neworder(acc_nodata - i + 1)
        end if
        if (i > 1 .and. acc_posdata(i) > 0.0_dp) then
          dpos = abs(acc_posdata(i) - acc_posdata(i - 1))
          if (dpos < epsilon(dpos)) cycle
        end if
        do j = 1, acc_nolabels
          write (acc_lineunit, '(A)', advance='NO') trim(sdl_i2s(acc_labeldata(k, j))) // ' '
        end do
        do j = 1, acc_noresults - 1
          write (acc_lineunit, '(ES20.11E3)', advance='NO') acc_resultdata(k, j)
        end do
        write (acc_lineunit, '(ES20.11E3)') acc_resultdata(k, acc_noresults)
      end do
      acc_nodata = 0
      if (associated(neworder)) then
        deallocate(neworder, stat=istat)
        if (istat /= 0) call Fatal(caller, 'Problems deallocating NewOrder vector')
      end if
      deallocate(acc_posdata, acc_labeldata, acc_resultdata, stat=istat)
      if (istat /= 0) call Fatal(caller, 'Problems deallocating some temporal workspace')
      nullify(acc_posdata, acc_labeldata, acc_resultdata)
    end if
    if (acc_fileisopen) then
      close (acc_lineunit, iostat=iostat)
      if (iostat /= 0) call Fatal(caller, 'Problems closing line file: ' // trim(sdl_i2s(iostat)))
      acc_fileisopen = .false.
    end if
  end subroutine sdl_close_line

  ! -- BoundaryFlux (verbatim internal procedure) ------------------------------
  subroutine sdl_boundary_flux(elem_idx, node_id, varname, coeffname, &
      f1, f2, fn, weight, maxn) bind(c, name="sdl_boundary_flux")
    integer(c_int), value :: elem_idx, node_id, maxn
    character(kind=c_char) :: varname(*), coeffname(*)
    real(c_double) :: f1, f2, fn, weight
    character(len=4096) :: vname, cname
    type(variable_t), pointer :: tvar
    type(element_t), pointer :: parent, element, oldelement
    type(nodes_t) :: nodesloc
    type(valuelist_t), pointer :: material, bc
    real(kind=dp) :: r, u, v, w, ub, detj, normal(3), flow(3)
    logical :: stat, permutated, gotit
    integer :: dimno, i, j, k, lbody, rbody, n, node, body_id, fluxbody
    real(kind=dp), pointer :: pwrk(:, :, :) => null()
    real(kind=dp), allocatable, save :: b_basis(:), b_dbasisdx(:, :)
    real(kind=dp), allocatable, save :: b_coefftensor(:, :, :)
    logical, save :: b_alloc = .false.
    call sdl_to_fortran(varname, vname)
    call sdl_to_fortran(coeffname, cname)
    node = node_id
    if (.not. b_alloc) then
      n = acc_mesh % maxelementnodes
      allocate (b_basis(n), b_dbasisdx(n, 3), b_coefftensor(3, 3, n))
      b_alloc = .true.
    end if
    tvar => VariableGet(acc_mesh % variables, trim(vname))
    if (.not. associated(tvar)) then
      call Fatal(caller, 'Cannot calculate fluxes without potential field!')
    end if
    permutated = associated(tvar % perm)
    element => acc_mesh % elements(elem_idx)
    bc => GetBC(element)
    fluxbody = ListGetInteger(bc, 'Flux Integrate Body', gotit)
    if (fluxbody > 0) then
      lbody = 0
      if (associated(element % boundaryinfo % left)) &
          lbody = element % boundaryinfo % left % bodyid
      rbody = 0
      if (associated(element % boundaryinfo % right)) &
          rbody = element % boundaryinfo % right % bodyid
      if (lbody == fluxbody) then
        parent => element % boundaryinfo % left
      else if (rbody == fluxbody) then
        parent => element % boundaryinfo % right
      else
        write (acc_message, *) 'No such flux integrate body on bc ', &
            element % boundaryinfo % constraint
        call Fatal(caller, acc_message)
      end if
    else
      parent => element % boundaryinfo % left
      stat = associated(parent)
      if (permutated) then
        if (stat) stat = all(tvar % perm(parent % nodeindexes) > 0)
        if (.not. stat) then
          parent => element % boundaryinfo % right
          stat = associated(parent)
          if (stat) stat = all(tvar % perm(parent % nodeindexes) > 0)
        end if
      end if
      if (.not. stat) call Fatal(caller, &
          'No solution available for specified boundary')
    end if
    n = parent % type % numberofnodes
    allocate(nodesloc % x(n), nodesloc % y(n), nodesloc % z(n))
    nodesloc % x(1:n) = acc_mesh % nodes % x(parent % nodeindexes)
    nodesloc % y(1:n) = acc_mesh % nodes % y(parent % nodeindexes)
    nodesloc % z(1:n) = acc_mesh % nodes % z(parent % nodeindexes)
    k = 0
    do j = 1, n
      if (node == parent % nodeindexes(j)) then
        k = j
        exit
      end if
    end do
    if (k == 0) then
      call Warn(caller, 'Side node not in parent element!')
    end if
    call GlobalToLocal(u, v, w, nodesloc % x(k), nodesloc % y(k), &
        nodesloc % z(k), parent, nodesloc)
    stat = ElementInfo(parent, nodesloc, u, v, w, detj, b_basis, b_dbasisdx)
    weight = detj
    if (.not. associated(acc_elementnodes % x)) then
      allocate(acc_elementnodes % x(maxn), acc_elementnodes % y(maxn), acc_elementnodes % z(maxn))
    end if
    acc_elementnodes % x(1:element % type % numberofnodes) = &
        acc_mesh % nodes % x(element % nodeindexes)
    acc_elementnodes % y(1:element % type % numberofnodes) = &
        acc_mesh % nodes % y(element % nodeindexes)
    acc_elementnodes % z(1:element % type % numberofnodes) = &
        acc_mesh % nodes % z(element % nodeindexes)
    do j = 1, element % type % numberofnodes
      if (node == element % nodeindexes(j)) exit
    end do
    if (j == 1) then
      ub = -1.0_dp
    else if (j == 2) then
      ub = 1.0_dp
    else
      ub = 0.0_dp
    end if
    normal = Normalvector(element, acc_elementnodes, ub, 0.0_dp, .true.)
    body_id = parent % bodyid
    k = ListGetInteger(acc_model % bodies(body_id) % values, 'Material', &
        minv=1, maxv=acc_model % numberofmaterials)
    material => acc_model % materials(k) % values
    if (.not. associated(material)) then
      call Warn(caller, 'Could not find material for flux computation!')
      deallocate(nodesloc % x, nodesloc % y, nodesloc % z)
      return
    end if
    oldelement => acc_model % currentelement
    acc_model % currentelement => parent
    call ListGetRealArray(material, trim(cname), pwrk, n, &
        parent % nodeindexes, stat)
    acc_model % currentelement => oldelement
    if (.not. associated(pwrk)) then
      call Warn(caller, 'Coefficient not present for flux computation!')
      deallocate(nodesloc % x, nodesloc % y, nodesloc % z)
      return
    end if
    b_coefftensor = 0.0_dp
    if (stat) then
      if (size(pwrk, 1) == 1) then
        do i = 1, 3
          b_coefftensor(i, i, 1:n) = pwrk(1, 1, 1:n)
        end do
      else if (size(pwrk, 2) == 1) then
        do i = 1, min(3, size(pwrk, 1))
          b_coefftensor(i, i, 1:n) = pwrk(i, 1, 1:n)
        end do
      else
        do i = 1, min(3, size(pwrk, 1))
          do j = 1, min(3, size(pwrk, 2))
            b_coefftensor(i, j, 1:n) = pwrk(i, j, 1:n)
          end do
        end do
      end if
    end if
    dimno = acc_dim
    flow = 0.0_dp
    do j = 1, dimno
      do k = 1, dimno
        if (permutated) then
          flow(j) = flow(j) + sum(b_coefftensor(j, k, 1:n) * b_basis(1:n)) * &
              sum(b_dbasisdx(1:n, k) * tvar % values(tvar % perm(parent % nodeindexes(1:n))))
        else
          flow(j) = flow(j) + sum(b_coefftensor(j, k, 1:n) * b_basis(1:n)) * &
              sum(b_dbasisdx(1:n, k) * tvar % values(parent % nodeindexes(1:n)))
        end if
      end do
    end do
    f1 = flow(1)
    f2 = flow(2)
    fn = sum(normal(1:dimno) * flow(1:dimno))
    deallocate(nodesloc % x, nodesloc % y, nodesloc % z)
  end subroutine sdl_boundary_flux


  subroutine sdl_to_fortran(cstr, fstr)
    character(kind=c_char) :: cstr(*)
    character(len=*) :: fstr
    integer :: i
    fstr = ' '
    do i = 1, len(fstr)
      if (cstr(i) == c_null_char) exit
      fstr(i:i) = cstr(i)
    end do
  end subroutine sdl_to_fortran

  ! -- SaveExistingLines (verbatim phase) --------------------------------------
  subroutine sdl_save_existing() bind(c, name="sdl_save_existing")
    character(len=4096) :: maskname
    type(valuelist_t), pointer :: vlist
    type(element_t), pointer :: el, parent
    integer :: elno, n, i, j, k, node, dgnode, istat, nsize, t
    integer, allocatable :: invperm(:), boundaryindex(:), nodetolement(:)
    integer, allocatable :: tmpidx(:)
    real(kind=dp) :: center(3), coord(3), coord0(3), linepos
    logical :: found, breakloop, parallelcomm, optimizeorder, dgvar
    ! NB: DGVar is intentionally never assigned here, mirroring the source
    ! (uninitialized local; verified by differentials, serial tree only).
    real(kind=dp) :: f1, f2, fn, weight
    maskname = ListGetString(acc_params, 'Save Mask', found)
    if (.not. found) maskname = 'Save Line'
    acc_maskname = maskname
    if (.not. (ListCheckPresentAnyBC(acc_model, maskname) .or. &
        ListCheckPresentAnyBodyForce(acc_model, maskname))) return
    call Info(caller, 'Saving existing nodes into ascii table', level=8)
    if (acc_solver % timesvisited > 0) then
      acc_initializeperm = (maskname /= acc_prevmaskname)
      acc_initializeperm = acc_initializeperm .or. acc_solver % meshchanged
    else
      acc_initializeperm = .true.
    end if
    if (acc_initializeperm) then
      acc_saveperm = 0
      optimizeorder = ListGetLogical(acc_params, 'Optimize Node Ordering', found)
      if (.not. found) optimizeorder = .not. acc_parallel
      breakloop = ListGetLogical(acc_params, 'Break Line Loop', found)
      if (breakloop) optimizeorder = .true.
      parallelcomm = acc_parallel
      if (acc_mesh % singlemesh) parallelcomm = .false.
      call MakePermUsingMask(acc_model, acc_solver, acc_mesh, maskname, &
          optimizeorder, acc_saveperm, acc_savenodes(1), &
          requirelogical=.true., breakloop=breakloop, parallelcomm=parallelcomm)
      if (acc_savenodes(1) > 0) then
        if (ListGetLogical(acc_params, 'Calculate Weights', found)) then
          call CalculateNodalWeights(acc_solver, .true., acc_saveperm, &
              trim(maskname) // ' Weights')
        end if
        call Info(caller, 'Number of nodes in specified boundary: ' // &
            trim(sdl_i2s(acc_savenodes(1))), level=12)
      end if
    else
      acc_savenodes(1) = 0
      if (associated(acc_saveperm)) acc_savenodes(1) = maxval(acc_saveperm)
    end if
    acc_prevmaskname = maskname
    if (acc_savenodes(1) <= 0) return
    allocate(invperm(acc_savenodes(1)), boundaryindex(acc_savenodes(1)), stat=istat)
    if (istat /= 0) call Fatal(caller, 'Memory allocation error 3')
    invperm = 0
    boundaryindex = 0
    do i = 1, size(acc_saveperm)
      if (acc_saveperm(i) > 0) then
        if (invperm(acc_saveperm(i)) > 0) then
          write (acc_message, *) 'Node multiple times in permutation', i, acc_saveperm(i)
          call Warn(caller, acc_message)
        end if
        invperm(acc_saveperm(i)) = i
      end if
    end do
    allocate(nodetolement(acc_mesh % numberofnodes))
    nodetolement = 0
    do i = 1, acc_mesh % numberofbulkelements + acc_mesh % numberofboundaryelements
      el => acc_mesh % elements(i)
      if (parenv % pes > 1) then
        if (el % partindex /= parenv % mype) cycle
      end if
      nodetolement(el % nodeindexes) = el % elementindex
    end do
    if (acc_calculateflux) then
      call Info(caller, 'Calculating nodal fluxes', level=8)
      allocate(acc_pointfluxes(acc_savenodes(1), 3), &
          acc_pointweight(acc_savenodes(1)), stat=istat)
      if (istat /= 0) call Fatal(caller, 'Memory allocation error 4')
      acc_pointfluxes = 0.0_dp
      acc_pointweight = 0.0_dp
      do i = acc_mesh % numberofbulkelements + 1, &
          acc_mesh % numberofbulkelements + acc_mesh % numberofboundaryelements
        el => acc_mesh % elements(i)
        if (parenv % pes > 1) then
          if (el % partindex /= parenv % mype) cycle
        end if
        acc_model % currentelement => el
        n = el % type % numberofnodes
        if (.not. all(acc_saveperm(el % nodeindexes) > 0)) cycle
        if (i <= acc_mesh % numberofbulkelements) then
          vlist => GetBodyForce()
        else
          vlist => GetBC()
        end if
        if (.not. associated(vlist)) cycle
        if (.not. ListCheckPresent(vlist, maskname)) cycle
        acc_elementnodes % x(1:n) = acc_mesh % nodes % x(el % nodeindexes)
        acc_elementnodes % y(1:n) = acc_mesh % nodes % y(el % nodeindexes)
        acc_elementnodes % z(1:n) = acc_mesh % nodes % z(el % nodeindexes)
        do j = 1, n
          node = el % nodeindexes(j)
          call sdl_boundary_flux(el % elementindex, node, acc_tempname, &
              acc_condname, f1, f2, fn, weight, 0)
          k = acc_saveperm(node)
          if (k == 0) cycle
          acc_pointfluxes(k, 1) = acc_pointfluxes(k, 1) + weight * f1
          acc_pointfluxes(k, 2) = acc_pointfluxes(k, 2) + weight * f2
          acc_pointfluxes(k, 3) = acc_pointfluxes(k, 3) + weight * fn
          acc_pointweight(k) = acc_pointweight(k) + weight
        end do
      end do
      do i = 1, acc_savenodes(1)
        acc_pointfluxes(i, 1) = acc_pointfluxes(i, 1) / acc_pointweight(i)
        acc_pointfluxes(i, 2) = acc_pointfluxes(i, 2) / acc_pointweight(i)
        acc_pointfluxes(i, 3) = acc_pointfluxes(i, 3) / acc_pointweight(i)
      end do
    end if
    acc_maxboundary = 0
    do i = 1, acc_mesh % numberofbulkelements + acc_mesh % numberofboundaryelements
      el => acc_mesh % elements(i)
      if (parenv % pes > 1) then
        if (el % partindex /= parenv % mype) cycle
      end if
      acc_model % currentelement => el
      n = el % type % numberofnodes
      if (.not. all(acc_saveperm(el % nodeindexes) > 0)) cycle
      if (i > acc_mesh % numberofbulkelements) then
        k = GetBCId(el)
        if (k == 0) cycle
        vlist => acc_model % bcs(k) % values
      else
        k = GetBodyForceId(el)
        if (k == 0) cycle
        vlist => acc_model % bodyforces(k) % values
      end if
      if (.not. ListCheckPresent(vlist, maskname)) cycle
      if (dgvar) then
        acc_elementnodes % x(1:n) = acc_mesh % nodes % x(el % nodeindexes)
        acc_elementnodes % y(1:n) = acc_mesh % nodes % y(el % nodeindexes)
        acc_elementnodes % z(1:n) = acc_mesh % nodes % z(el % nodeindexes)
        center = 0.0_dp
        center(1) = sum(acc_elementnodes % x(1:n)) / n
        center(2) = sum(acc_elementnodes % y(1:n)) / n
        center(3) = sum(acc_elementnodes % z(1:n)) / n
        do j = 1, n
          node = el % nodeindexes(j)
          coord(1) = acc_elementnodes % x(j)
          coord(2) = acc_elementnodes % y(j)
          coord(3) = acc_elementnodes % z(j)
          coord0 = coord
          coord = center + 0.9999_dp * (coord - center)
          acc_mesh % nodes % x(node) = coord(1)
          acc_mesh % nodes % y(node) = coord(2)
          if (acc_dim == 3) acc_mesh % nodes % z(node) = coord(3)
          linepos = -1.0_dp
        if (acc_calculateflux) then
          call sdl_write_fields_impl(el % elementindex, k, node, dgnode, &
              .true., acc_pointfluxes(t, :), .true., linepos, .true., &
              acc_parallel)
        else
          call sdl_write_fields_impl(el % elementindex, k, node, dgnode, &
              .true., acc_dummy3, .false., linepos, .true., acc_parallel)
        end if
          acc_mesh % nodes % x(node) = coord0(1)
          acc_mesh % nodes % y(node) = coord0(2)
          if (acc_dim == 3) acc_mesh % nodes % z(node) = coord0(3)
        end do
      else
        boundaryindex(acc_saveperm(el % nodeindexes)) = k
      end if
      acc_maxboundary = max(acc_maxboundary, k)
    end do
    if (.not. dgvar) then
      dgnode = 0
      linepos = -1.0_dp
      do i = 1, acc_savenodes(1)
        node = invperm(i)
        el => acc_mesh % elements(nodetolement(node))
        if (acc_calculateflux) then
          call sdl_write_fields_impl(el % elementindex, boundaryindex(i), &
              node, dgnode, .true., acc_pointfluxes(i, :), .true., linepos, &
              .true., acc_parallel)
        else
          call sdl_write_fields_impl(el % elementindex, boundaryindex(i), &
              node, dgnode, .true., acc_dummy3, .false., linepos, .true., &
              acc_parallel)
        end if
      end do
    end if
    deallocate(invperm, boundaryindex, stat=istat)
    if (istat /= 0) call Fatal(caller, 'Problems deallocating stuff on existing lines')
    if (acc_calculateflux) then
      deallocate(acc_pointfluxes, acc_pointweight, stat=istat)
      if (istat /= 0) call Fatal(caller, 'Problems deallocating fluxes on existing lines')
    end if
  end subroutine sdl_save_existing

end module sd_line_acc

