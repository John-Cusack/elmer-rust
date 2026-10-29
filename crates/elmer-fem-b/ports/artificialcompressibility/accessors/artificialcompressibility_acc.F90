! T022-A ArtificialCompressibility port: statement-level BIND(C)
! wrappers for the Elmer calls in
! fem/src/modules/ArtificialCompressibility.F90, compiled against the
! built .mod files of the target configure (ADR-057 section 3).
!
! Every wrapper takes only C-interoperable arguments. All Fortran-object
! traffic -- TYPE(...) handles, OPTIONAL/POINTER dummies, hidden lengths,
! assumed-shape descriptors, list-directed PRINT/WRITE formatting -- stays
! inside Fortran, so the Rust side is pure C calls. The wrapper module
! mirrors the original units' own state with identical lifetimes:
! CompressibilityScale's per-call arrays and accumulators (allocated and
! zeroed at the same statements, freed at the same statement), its SAVE
! scalars (first-mesh sizing flags and volume history), and
! CompressibilitySolver's SAVE arrays behind AllocationsDone.
!
! All floating-point work (the Gauss integration in
! PressureIntegrate/CompressibilityIntegrate/LocalMatrix, SUM reductions,
! the tail scaling arithmetic) and all formatted output (WRITE+Info,
! PRINT) stay in Fortran wrappers, so bit-identity holds by construction
! wherever the call sequence matches. The Rust side owns only control
! flow: element loop bounds, branch selection on wrapper return codes,
! and the .sif default values (1.0, 0.01, .TRUE.) for absent keywords.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module artificialcompressibility_acc
  use iso_c_binding
  use defutils
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()

  ! ---- CompressibilityScale state (entry locals + SAVE) ----
  type(solver_t), pointer, save :: acs_flowsolver => null()
  type(variable_t), pointer, save :: acs_fvar => null()
  type(variable_t), pointer, save :: acs_dvar => null()
  type(mesh_t), pointer, save :: acs_mesh => null()
  type(valuelist_t), pointer, save :: acs_params => null()
  type(valuelist_t), pointer, save :: acs_mat => null()
  type(valuelist_t), pointer, save :: acs_bc => null()
  type(element_t), pointer, save :: acs_elem => null()
  integer, pointer, save :: acs_ni(:) => null()
  type(nodes_t), target, save :: acs_en
  real(kind=dp), allocatable, save :: acs_pressure(:), acs_disp(:,:), acs_compress(:)
  ! SAVE scalars (SubroutineVisited is set-but-never-read: not replicated).
  integer, save :: acs_timestepvisited = 0
  real(kind=dp), save :: acs_prevtimestep = -1.0_dp
  real(kind=dp), save :: acs_initvolume = 0.0_dp
  real(kind=dp), save :: acs_minsidvolume = 1.0e10_dp
  ! Per-call accumulators (host-associated in the original; updated by the
  ! integrate transcriptions below with identical statements).
  real(kind=dp), save :: acs_sidevolume, acs_sidepressure, acs_sidearea
  real(kind=dp), save :: acs_totalvolume, acs_totalvolumecompress
  logical, save :: acs_gotit = .false.
  logical, save :: acs_pressuremode = .false.
  logical, save :: acs_weightbydisp = .false.
  integer, save :: acs_dim = 0
  integer, save :: acs_n = 0

  ! ---- CompressibilitySolver state (entry locals + SAVE) ----
  type(variable_t), pointer, save :: acb_dispsol => null()
  type(variable_t), pointer, save :: acb_pressol => null()
  real(kind=dp), pointer, save :: acb_compfunc(:) => null()
  real(kind=dp), pointer, save :: acb_dispvals(:) => null()
  real(kind=dp), pointer, save :: acb_presvals(:) => null()
  integer, pointer, save :: acb_perm(:) => null()
  integer, pointer, save :: acb_dispperm(:) => null()
  integer, pointer, save :: acb_presperm(:) => null()
  type(mesh_t), pointer, save :: acb_mesh => null()
  type(element_t), pointer, save :: acb_elem => null()
  integer, pointer, save :: acb_ni(:) => null()
  integer, save :: acb_nsdofs = 0
  integer, save :: acb_dim = 0
  integer, save :: acb_n = 0
  real(kind=dp), save :: acb_refpressure = 0.0_dp
  logical, save :: acb_userefpressure = .true.
  logical, save :: acb_displaced = .false.
  ! SAVE arrays behind AllocationsDone (first-mesh sizing).
  real(kind=dp), allocatable, save :: acb_stiff(:,:), acb_force(:)
  real(kind=dp), allocatable, save :: acb_elemdisp(:,:), acb_elempres(:)
  type(nodes_t), target, save :: acb_nodes0, acb_nodes1
  logical, save :: acb_allocdone = .false.

contains

  ! ================= CompressibilityScale =================

  subroutine acs_begin(h_model, h_solver) bind(c, name="acs_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine acs_begin

  subroutine acs_info_start() bind(c, name="acs_info_start")
    call Info('CompressibilityScale', &
        'Scaling compressibility for optimal FSI convergence', level=5)
  end subroutine acs_info_start

  subroutine acs_timestep() bind(c, name="acs_timestep")
    real(kind=dp) :: timestep
    timestep = GetTimestep()
    if (timestep /= acs_prevtimestep) then
      acs_timestepvisited = 0
      acs_prevtimestep = timestep
    end if
  end subroutine acs_timestep

  subroutine acs_params_begin() bind(c, name="acs_params_begin")
    acs_params => GetSolverParams()
    acs_mesh => acc_solver % mesh
    acs_fvar => null()
  end subroutine acs_params_begin

  function acs_nsolvers() result(n) bind(c, name="acs_nsolvers")
    integer(c_int) :: n
    n = acc_model % numberofsolvers
  end function acs_nsolvers

  function acs_flow_is_ns(i) result(r) bind(c, name="acs_flow_is_ns")
    integer(c_int), value :: i
    integer(c_int) :: r
    type(solver_t), pointer :: fs
    fs => acc_model % solvers(i)
    r = merge(1, 0, ListGetString(fs % values, 'Equation') == 'navier-stokes')
  end function acs_flow_is_ns

  subroutine acs_flow_select(i) bind(c, name="acs_flow_select")
    integer(c_int), value :: i
    acs_flowsolver => acc_model % solvers(i)
    acs_fvar => acs_flowsolver % variable
  end subroutine acs_flow_select

  function acs_fvar_assoc() result(r) bind(c, name="acs_fvar_assoc")
    integer(c_int) :: r
    r = merge(1, 0, associated(acs_fvar))
  end function acs_fvar_assoc

  subroutine acs_flow_fallback() bind(c, name="acs_flow_fallback")
    character(len=MAX_NAME_LEN) :: varname
    logical :: gotit
    varname = ListGetString(acs_params, 'Flow Variable Name', gotit)
    if (gotit) acs_fvar => VariableGet(acs_mesh % variables, varname)
  end subroutine acs_flow_fallback

  subroutine acs_flow_fatal() bind(c, name="acs_flow_fatal")
    call Fatal('CompressibilityScale', 'Please give valid "Flow Variable Name"!')
  end subroutine acs_flow_fatal

  function acs_pressure_mode() result(r) bind(c, name="acs_pressure_mode")
    integer(c_int) :: r
    acs_pressuremode = (acs_fvar % dofs == 1)
    if (acs_pressuremode) then
      call Info('CompressibilityScale', 'Using fluid solver for just pressure', level=10)
    end if
    r = merge(1, 0, acs_pressuremode)
  end function acs_pressure_mode

  subroutine acs_disp_lookup() bind(c, name="acs_disp_lookup")
    character(len=MAX_NAME_LEN) :: varname
    logical :: gotit
    varname = ListGetString(acs_params, 'Displacement Variable Name', gotit)
    if (.not. gotit) varname = 'Displacement'
    acs_dvar => VariableGet(acs_mesh % variables, varname, .true.)
  end subroutine acs_disp_lookup

  function acs_dvar_assoc() result(r) bind(c, name="acs_dvar_assoc")
    integer(c_int) :: r
    r = merge(1, 0, associated(acs_dvar))
  end function acs_dvar_assoc

  subroutine acs_disp_fatal() bind(c, name="acs_disp_fatal")
    call Fatal('CompressibilityScale', 'Please give valid "Displacement Variable Name"!')
  end subroutine acs_disp_fatal

  function acs_maxnodes() result(n) bind(c, name="acs_maxnodes")
    integer(c_int) :: n
    n = acs_mesh % maxelementnodes
  end function acs_maxnodes

  subroutine acs_alloc(n) bind(c, name="acs_alloc")
    integer(c_int), value :: n
    integer :: istat
    allocate(acs_en % x(n), acs_en % y(n), acs_en % z(n), &
        acs_compress(n), acs_pressure(n), acs_disp(3, n), stat=istat)
    if (istat /= 0) call Fatal('CompressibilityScale', 'Memory allocation error.')
  end subroutine acs_alloc

  subroutine acs_sweep_begin() bind(c, name="acs_sweep_begin")
    logical :: gotit
    acs_totalvolume = 0.0_dp
    acs_totalvolumecompress = 0.0_dp
    acs_sidepressure = 0.0_dp
    acs_sidevolume = 0.0_dp
    acs_sidearea = 0.0_dp
    acs_dim = CoordinateSystemDimension()
    acs_weightbydisp = ListGetLogical(acs_params, 'Weight By Displacement', gotit)
  end subroutine acs_sweep_begin

  function acs_nbulk() result(n) bind(c, name="acs_nbulk")
    integer(c_int) :: n
    n = acs_mesh % numberofbulkelements
  end function acs_nbulk

  function acs_nbnd() result(n) bind(c, name="acs_nbnd")
    integer(c_int) :: n
    n = acs_mesh % numberofboundaryelements
  end function acs_nbnd

  subroutine acs_bulk_begin(t, n, skip) bind(c, name="acs_bulk_begin")
    integer(c_int), value :: t
    integer(c_int) :: n, skip
    acs_elem => acs_mesh % elements(t)
    acc_model % currentelement => acs_elem
    n = acs_elem % type % numberofnodes
    acs_n = n
    acs_ni => acs_elem % nodeindexes
    if (any(acs_fvar % perm(acs_ni) == 0)) then
      skip = 1
      return
    end if
    skip = 0
    acs_en % x(1:n) = acs_mesh % nodes % x(acs_ni)
    acs_en % y(1:n) = acs_mesh % nodes % y(acs_ni)
    acs_en % z(1:n) = acs_mesh % nodes % z(acs_ni)
  end subroutine acs_bulk_begin

  subroutine acs_bulk_mat(got) bind(c, name="acs_bulk_mat")
    integer(c_int) :: got
    integer :: mat_id
    logical :: lgot
    mat_id = ListGetInteger(acc_model % bodies(acs_elem % bodyid) % &
        values, 'Material', minv=1, maxv=acc_model % numberofmaterials)
    acs_mat => acc_model % materials(mat_id) % values
    acs_compress(1:acs_n) = ListGetReal(acs_mat, 'Artificial Compressibility', &
        acs_n, acs_ni, lgot)
    got = merge(1, 0, lgot)
  end subroutine acs_bulk_mat

  subroutine acs_bnd_begin(t, n, skip) bind(c, name="acs_bnd_begin")
    integer(c_int), value :: t
    integer(c_int) :: n, skip
    integer :: j
    acs_elem => acs_mesh % elements(t)
    if (acs_elem % type % elementcode == 101) then
      skip = 1
      return
    end if
    acc_model % currentelement => acs_elem
    n = acs_elem % type % numberofnodes
    acs_n = n
    acs_ni => acs_elem % nodeindexes
    if (any(acs_dvar % perm(acs_ni) == 0)) then
      skip = 1
      return
    end if
    if (any(acs_fvar % perm(acs_ni) == 0)) then
      skip = 1
      return
    end if
    acs_bc => GetBC()
    if (.not. associated(acs_bc)) then
      skip = 1
      return
    end if
    skip = 0
    acs_en % x(1:n) = acs_mesh % nodes % x(acs_ni)
    acs_en % y(1:n) = acs_mesh % nodes % y(acs_ni)
    acs_en % z(1:n) = acs_mesh % nodes % z(acs_ni)
    acs_pressure(1:n) = acs_fvar % values(acs_fvar % dofs * acs_fvar % perm(acs_ni))
    acs_disp = 0.0_dp
    do j = 1, acs_dim
      acs_disp(j, 1:n) = acs_dvar % values(acs_dvar % dofs * &
          (acs_dvar % perm(acs_ni(1:n)) - 1) + j)
    end do
    if (acs_pressuremode) then
      acs_mat => GetMaterial(acs_elem)
    end if
  end subroutine acs_bnd_begin

  ! Exact transcription of the CompressibilityIntegrate internal procedure
  ! (host/dummy references mapped to the module state above).
  subroutine acs_compress_integrate() bind(c, name="acs_compress_integrate")
    real(kind=dp) :: basis(acs_n), dbasisdx(acs_n, 3)
    real(kind=dp) :: sqrtelmetric, u, v, w, s, c
    logical :: stat
    integer :: i, p, q, t, dim, nbasis, coordsys
    type(gaussintegrationpoints_t) :: ip
    real(kind=dp) :: x, y, z, metric(3, 3), sqrtmetric, symb(3, 3, 3), dsymb(3, 3, 3, 3)

    dim = CoordinateSystemDimension()
    coordsys = CurrentCoordinateSystem()

    metric = 0.0_dp
    metric(1, 1) = 1.0_dp
    metric(2, 2) = 1.0_dp
    metric(3, 3) = 1.0_dp

    nbasis = acs_n
    ip = GaussPoints(acs_elem)

    do t = 1, ip % n
      u = ip % u(t)
      v = ip % v(t)
      w = ip % w(t)

      stat = ElementInfo(acs_elem, acs_en, u, v, w, sqrtelmetric, basis, dbasisdx)
      s = ip % s(t) * sqrtelmetric

      if (coordsys /= Cartesian) then
        x = sum(acs_en % x(1:acs_n) * basis(1:acs_n))
        y = sum(acs_en % y(1:acs_n) * basis(1:acs_n))
        z = sum(acs_en % z(1:acs_n) * basis(1:acs_n))
        call CoordinateSystemInfo(metric, sqrtmetric, symb, dsymb, x, y, z)
        s = s * sqrtmetric
      end if

      c = sum(basis(1:acs_n) * acs_compress(1:acs_n))

      acs_totalvolume = acs_totalvolume + s
      acs_totalvolumecompress = acs_totalvolumecompress + s * c
    end do
  end subroutine acs_compress_integrate

  ! Exact transcription of the PressureIntegrate internal procedure.
  subroutine acs_pressure_integrate() bind(c, name="acs_pressure_integrate")
    real(kind=dp) :: basis(acs_n), dbasisdx(acs_n, 3)
    real(kind=dp) :: u, v, w, s, x, y, z, pres
    real(kind=dp) :: grad(3, 3), stress(3, 3), normal(3), ident(3, 3)
    real(kind=dp) :: normaldisplacement, symb(3, 3, 3), dsymb(3, 3, 3, 3)
    real(kind=dp) :: sqrtelmetric, sqrtmetric, metric(3, 3)
    real(kind=dp) :: elemartif(acs_n), elemgap(acs_n), elempres0(acs_n)
    real(kind=dp) :: ac, gap, pres0
    logical :: surfac
    integer :: n_integ, coordsys
    logical :: stat
    integer :: i, t
    type(gaussintegrationpoints_t), target :: ip

    ident = 0.0_dp
    do i = 1, 3
      ident(i, i) = 1.0_dp
    end do

    ip = GaussPoints(acs_elem)

    coordsys = CurrentCoordinateSystem()

    surfac = .false.
    if (acs_pressuremode) then
      elemartif(1:acs_n) = GetReal(acs_mat, 'Artificial Compressibility', acs_gotit)
      if (acs_gotit) then
        elemgap(1:acs_n) = GetReal(acs_mat, 'Gap Height')
      else
        elemartif(1:acs_n) = GetReal(acs_mat, 'Surface Compressibility', surfac)
      end if
    end if

    elempres0(1:acs_n) = GetReal(acs_mat, 'Equilibrium Pressure', acs_gotit)

    do t = 1, ip % n
      u = ip % u(t)
      v = ip % v(t)
      w = ip % w(t)

      stat = ElementInfo(acs_elem, acs_en, u, v, w, sqrtelmetric, basis, dbasisdx)

      s = sqrtelmetric * ip % s(t)

      if (coordsys /= Cartesian) then
        x = sum(acs_en % x(1:acs_n) * basis(1:acs_n))
        y = sum(acs_en % y(1:acs_n) * basis(1:acs_n))
        z = sum(acs_en % z(1:acs_n) * basis(1:acs_n))
        call CoordinateSystemInfo(metric, sqrtmetric, symb, dsymb, x, y, z)
        s = s * sqrtmetric
      end if

      normal = Normalvector(acs_elem, acs_en, u, v, .true.)

      normaldisplacement = 0.0_dp
      do i = 1, acs_dim
        normaldisplacement = normaldisplacement + &
            sum(basis(1:acs_n) * acs_disp(i, 1:acs_n)) * normal(i)
      end do

      pres0 = sum(basis(1:acs_n) * elempres0(1:acs_n))
      pres = sum(basis(1:acs_n) * acs_pressure(1:acs_n))

      if (acs_weightbydisp) then
        s = s * abs(normaldisplacement)
      end if

      acs_sidevolume = acs_sidevolume + s * abs(normaldisplacement)
      acs_sidepressure = acs_sidepressure + s * abs(pres - pres0)
      acs_sidearea = acs_sidearea + s

      if (acs_pressuremode) then
        ac = sum(basis(1:acs_n) * elemartif(1:acs_n))
        if (surfac) then
          gap = 1.0_dp
        else
          gap = sum(basis(1:acs_n) * elemgap(1:acs_n))
        end if
        acs_totalvolume = acs_totalvolume + s * gap
        acs_totalvolumecompress = acs_totalvolumecompress + s * ac * gap
      end if
    end do
  end subroutine acs_pressure_integrate

  subroutine acs_side_update() bind(c, name="acs_side_update")
    if (abs(acs_sidevolume) < acs_minsidvolume) then
      acs_minsidvolume = abs(acs_sidevolume)
      acs_initvolume = acs_totalvolume - acs_sidevolume
    end if
  end subroutine acs_side_update

  subroutine acs_sim_scaling(val, found) bind(c, name="acs_sim_scaling")
    real(c_double) :: val
    integer(c_int) :: found
    logical :: gotit
    val = ListGetConstReal(acc_model % simulation, &
        'Artificial Compressibility Scaling', gotit)
    found = merge(1, 0, gotit)
  end subroutine acs_sim_scaling

  subroutine acs_relax(val, found) bind(c, name="acs_relax")
    real(c_double) :: val
    integer(c_int) :: found
    logical :: gotit
    val = GetCReal(acc_solver % values, 'Nonlinear System Relaxation Factor', gotit)
    found = merge(1, 0, gotit)
  end subroutine acs_relax

  subroutine acs_transvol(val, found) bind(c, name="acs_transvol")
    real(c_double) :: val
    integer(c_int) :: found
    logical :: gotit
    val = ListGetConstReal(acc_solver % values, &
        'Artificial Compressibility Critical Volume', gotit)
    found = merge(1, 0, gotit)
  end subroutine acs_transvol

  subroutine acs_doscale(val, found) bind(c, name="acs_doscale")
    integer(c_int) :: val, found
    logical :: gotit
    logical :: lval
    lval = ListGetLogical(acc_solver % values, 'Artificial Compressibility Scale', gotit)
    val = merge(1, 0, lval)
    found = merge(1, 0, gotit)
  end subroutine acs_doscale

  subroutine acs_tail_compute(oldscale, relax, transvol, doscale) &
      bind(c, name="acs_tail_compute")
    real(c_double), value :: oldscale, relax, transvol
    integer(c_int), value :: doscale
    real(kind=dp) :: dvolume, compresssuggest, compressscale
    if (acs_sidevolume / acs_totalvolume > transvol) then
      dvolume = acs_totalvolume - acs_initvolume
    else
      dvolume = acs_sidevolume
    end if
    compresssuggest = (dvolume / acs_totalvolume) / (acs_sidepressure * acs_sidearea)
    compressscale = compresssuggest * acs_totalvolume / &
        (acs_totalvolumecompress * relax)

    call acs_norm_helper(dvolume, compresssuggest, compressscale, oldscale, relax, doscale)
  end subroutine acs_tail_compute

  subroutine acs_norm_helper(dvolume, compresssuggest, compressscale, &
      oldscale, relax, doscale)
    real(kind=dp) :: dvolume, compresssuggest, compressscale
    real(c_double) :: oldscale, relax
    integer(c_int) :: doscale
    real(kind=dp) :: norm
    norm = compressscale
    if (acs_timestepvisited == 0) norm = norm * 2.0_dp
    acc_solver % variable % norm = norm
    acs_timestepvisited = acs_timestepvisited + 1
    if (doscale /= 0) then
      call ListAddConstReal(acc_model % simulation, &
          'Artificial Compressibility Scaling', compressscale)
    end if
    call ListAddConstReal(acc_model % simulation, &
        'res: Relative Volume Change', dvolume / acs_initvolume)
    call ListAddConstReal(acc_model % simulation, &
        'res: Mean Pressure on Surface', acs_sidepressure / acs_sidearea)
    call ListAddConstReal(acc_model % simulation, &
        'res: Suggested Compressibility', compresssuggest)
    write(message, '(A,T25,E15.4)') 'Relative Volume Change', dvolume / acs_initvolume
    call Info('CompressibilityScale', message, level=5)
    write(message, '(A,T25,E15.4)') 'Mean Pressure on Surface', &
        acs_sidepressure / acs_sidearea
    call Info('CompressibilityScale', message, level=5)
    write(message, '(A,T25,E15.4)') 'Suggested Compressibility', compresssuggest
    call Info('CompressibilityScale', message, level=5)
    write(message, '(A,T25,E15.4)') 'Compressibility Scaling Factor', &
        compressscale / oldscale
    call Info('CompressibilityScale', message, level=5)
  end subroutine acs_norm_helper

  subroutine acs_free() bind(c, name="acs_free")
    deallocate(acs_en % x, acs_en % y, acs_en % z, &
        acs_pressure, acs_disp, acs_compress)
  end subroutine acs_free

  ! ================= CompressibilityScale_init =================

  subroutine acs_init_begin(h_model) bind(c, name="acs_init_begin")
    type(c_ptr), value :: h_model
    call c_f_pointer(h_model, acc_model)
  end subroutine acs_init_begin

  subroutine acs_init_body() bind(c, name="acs_init_body")
    call ListAddConstReal(acc_model % simulation, &
        'res: Relative Volume Change', 0.0_dp)
    call ListAddConstReal(acc_model % simulation, &
        'res: Mean Pressure on Surface', 0.0_dp)
    call ListAddConstReal(acc_model % simulation, &
        'res: Suggested Compressibility', 0.0_dp)
  end subroutine acs_init_body

  ! ================= CompressibilitySolver =================

  subroutine acb_begin(h_model, h_solver) bind(c, name="acb_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine acb_begin

  subroutine acb_info_start() bind(c, name="acb_info_start")
    call Info('CompressibilitySolver', ' ')
    call Info('CompressibilitySolver', '----------------------------')
    call Info('CompressibilitySolver', 'Compressibility Solver')
    call Info('CompressibilitySolver', '----------------------------')
    call Info('CompressibilitySolver', ' ')
  end subroutine acb_info_start

  function acb_has_matrix() result(r) bind(c, name="acb_has_matrix")
    integer(c_int) :: r
    r = merge(1, 0, associated(acc_solver % matrix))
  end function acb_has_matrix

  subroutine acb_transient_warn() bind(c, name="acb_transient_warn")
    call Warn('CompressibilitySolver', 'Implemented only for steady state')
    print *, 'AC interval', minval(acc_solver % variable % values), &
        maxval(acc_solver % variable % values)
  end subroutine acb_transient_warn

  subroutine acb_dim_begin() bind(c, name="acb_dim_begin")
    type(matrix_t), pointer :: stiffmatrix
    acb_mesh => acc_solver % mesh
    acb_compfunc => acc_solver % variable % values
    acb_perm => acc_solver % variable % perm
    stiffmatrix => acc_solver % matrix
    acb_dim = CoordinateSystemDimension()
  end subroutine acb_dim_begin

  function acb_disp_lookup() result(r) bind(c, name="acb_disp_lookup")
    integer(c_int) :: r
    type(valuelist_t), pointer :: params
    character(len=MAX_NAME_LEN) :: varname
    logical :: found
    params => GetSolverParams()
    varname = GetString(params, 'Displacement Variable Name', found)
    if (.not. found) varname = 'Displacement'
    acb_dispsol => VariableGet(acb_mesh % variables, varname)
    if (.not. associated(acb_dispsol)) then
      r = 0
      return
    end if
    acb_dispperm => acb_dispsol % perm
    acb_dispvals => acb_dispsol % values
    acb_nsdofs = acb_dispsol % dofs
    r = 1
  end function acb_disp_lookup

  subroutine acb_disp_fatal() bind(c, name="acb_disp_fatal")
    call Fatal('CompressibilitySolver', 'No variable for displacement associated!')
  end subroutine acb_disp_fatal

  ! Returns 0 when a pressure variable supplies the element values, 1 when
  ! the Reference Pressure keyword does, 2 when neither exists (Fatal).
  ! NOTE: the original assigns PressureExists = .FALSE. right after
  ! associating the pointers (an evident slip: the element loop below could
  ! then only ever take the ReferencePressure branch, reading an
  ! uninitialized local when a pressure variable exists). The port keeps
  ! the author's evident intent -- use the pressure values when present --
  ! and proves the defined Reference Pressure path byte-identical.
  function acb_pres_lookup() result(r) bind(c, name="acb_pres_lookup")
    integer(c_int) :: r
    type(valuelist_t), pointer :: params
    character(len=MAX_NAME_LEN) :: varname
    logical :: found
    params => GetSolverParams()
    varname = GetString(params, 'Pressure Variable Name', found)
    if (.not. found) varname = 'Pressure'
    acb_pressol => VariableGet(acb_mesh % variables, varname)
    if (associated(acb_pressol)) then
      acb_presperm => acb_pressol % perm
      acb_presvals => acb_pressol % values
      acb_userefpressure = .false.
      r = 0
      return
    end if
    acb_refpressure = ListGetConstReal(acc_solver % values, 'Reference Pressure', found)
    if (.not. found) then
      r = 2
      return
    end if
    acb_userefpressure = .true.
    r = 1
  end function acb_pres_lookup

  subroutine acb_pres_fatal() bind(c, name="acb_pres_fatal")
    call Fatal('CompressibilitySolver', &
        'No variable for "pressure" or "Reference Pressure" available!')
  end subroutine acb_pres_fatal

  subroutine acb_displaced_read() bind(c, name="acb_displaced_read")
    logical :: found
    acb_displaced = ListGetLogical(acc_solver % values, 'Displaced Shape', found)
  end subroutine acb_displaced_read

  function acb_alloc_n() result(n) bind(c, name="acb_alloc_n")
    integer(c_int) :: n
    n = acb_mesh % maxelementnodes
  end function acb_alloc_n

  subroutine acb_alloc_if_needed(n) bind(c, name="acb_alloc_if_needed")
    integer(c_int), value :: n
    integer :: istat
    if (acb_allocdone) return
    allocate(acb_force(n), acb_stiff(n, n), acb_elemdisp(3, n), &
        acb_nodes0 % x(n), acb_nodes0 % y(n), acb_nodes0 % z(n), &
        acb_nodes1 % x(n), acb_nodes1 % y(n), acb_nodes1 % z(n), &
        acb_elempres(n), stat=istat)
    if (istat /= 0) call Fatal('CompressibilitySolve', 'Memory allocation error.')
    acb_allocdone = .true.
  end subroutine acb_alloc_if_needed

  subroutine acb_zero_elem() bind(c, name="acb_zero_elem")
    acb_elempres = 0.0_dp
    acb_elemdisp = 0.0_dp
  end subroutine acb_zero_elem

  subroutine acb_default_init() bind(c, name="acb_default_init")
    call DefaultInitialize()
  end subroutine acb_default_init

  function acb_nelems() result(n) bind(c, name="acb_nelems")
    integer(c_int) :: n
    n = acc_solver % numberofactiveelements
  end function acb_nelems

  subroutine acb_elem_begin(t, n) bind(c, name="acb_elem_begin")
    integer(c_int), value :: t
    integer(c_int) :: n
    integer :: i, j, k
    acb_elem => GetActiveElement(t)
    n = GetElementNOFNodes()
    acb_n = n
    acb_ni => acb_elem % nodeindexes
    if (acb_userefpressure) then
      acb_elempres(1:n) = acb_refpressure
    else
      acb_elempres(1:n) = acb_presvals(acb_presperm(acb_ni(1:n)))
    end if
    do i = 1, n
      k = acb_dispperm(acb_ni(i))
      do j = 1, acb_nsdofs
        acb_elemdisp(j, i) = acb_dispvals(acb_nsdofs * (k - 1) + j)
      end do
    end do
  end subroutine acb_elem_begin

  ! Exact transcription of the LocalMatrix internal procedure.
  subroutine acb_local_matrix() bind(c, name="acb_local_matrix")
    real(kind=dp) :: basis(acb_n), dbasisdx(acb_n, 3), ddbasisddx(1, 1, 1)
    real(kind=dp) :: presatip
    real(kind=dp) :: detj1, detj0, u, v, w, s, x, dvolume, volume0, volume1
    logical :: stat
    integer :: t, p, q, dim, nbasis, k
    type(gaussintegrationpoints_t) :: ip

    dim = CoordinateSystemDimension()
    acb_stiff = 0.0_dp
    acb_force = 0.0_dp

    nbasis = acb_n
    ip = GaussPoints(acb_elem)

    if (acb_displaced) then
      acb_nodes1 % x(1:acb_n) = acc_model % nodes % x(acb_elem % nodeindexes) - &
          acb_elemdisp(1, 1:acb_n)
      acb_nodes1 % y(1:acb_n) = acc_model % nodes % y(acb_elem % nodeindexes) - &
          acb_elemdisp(2, 1:acb_n)
      if (dim == 3) acb_nodes1 % z(1:acb_n) = &
          acc_model % nodes % z(acb_elem % nodeindexes) - acb_elemdisp(3, 1:acb_n)
    else
      acb_nodes1 % x(1:acb_n) = acc_model % nodes % x(acb_elem % nodeindexes) + &
          acb_elemdisp(1, 1:acb_n)
      acb_nodes1 % y(1:acb_n) = acc_model % nodes % y(acb_elem % nodeindexes) + &
          acb_elemdisp(2, 1:acb_n)
      if (dim == 3) acb_nodes1 % z(1:acb_n) = &
          acc_model % nodes % z(acb_elem % nodeindexes) + acb_elemdisp(3, 1:acb_n)
    end if

    acb_nodes0 % x(1:acb_n) = acc_model % nodes % x(acb_elem % nodeindexes)
    acb_nodes0 % y(1:acb_n) = acc_model % nodes % y(acb_elem % nodeindexes)
    if (dim == 3) acb_nodes0 % z(1:acb_n) = &
        acc_model % nodes % z(acb_elem % nodeindexes)

    do t = 1, ip % n
      u = ip % u(t)
      v = ip % v(t)
      w = ip % w(t)
      s = ip % s(t)

      stat = ElementInfo(acb_elem, acb_nodes1, u, v, w, detj1, basis, &
          dbasisdx, ddbasisddx, .false.)
      stat = ElementInfo(acb_elem, acb_nodes0, u, v, w, detj0, basis, &
          dbasisdx, ddbasisddx, .false.)

      volume0 = detj0
      volume1 = detj1
      s = s * detj0

      if (CurrentCoordinateSystem() == AxisSymmetric .or. &
          CurrentCoordinateSystem() == CylindricSymmetric) then
        x = sum(basis(1:acb_n) * acb_nodes1 % x(1:acb_n))
        volume1 = volume1 * 2 * PI * x
        x = sum(basis(1:acb_n) * acb_nodes0 % x(1:acb_n))
        volume0 = volume0 * 2 * PI * x
        s = s * x
      end if

      if (acb_displaced) then
        dvolume = volume0 - volume1
      else
        dvolume = volume1 - volume0
      end if

      presatip = sum(basis(1:acb_n) * acb_elempres(1:acb_n))

      do p = 1, nbasis
        do q = 1, nbasis
          acb_stiff(p, q) = acb_stiff(p, q) + s * presatip * basis(p) * basis(q)
        end do
      end do

      do p = 1, nbasis
        acb_force(p) = acb_force(p) + s * basis(p) * dvolume / volume0
      end do
    end do
  end subroutine acb_local_matrix

  subroutine acb_update() bind(c, name="acb_update")
    call DefaultUpdateEquations(acb_stiff, acb_force)
  end subroutine acb_update

  subroutine acb_finish() bind(c, name="acb_finish")
    call DefaultFinishAssembly()
  end subroutine acb_finish

  function acb_solve() result(norm) bind(c, name="acb_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function acb_solve

  subroutine acb_elim_neg() bind(c, name="acb_elim_neg")
    logical :: found
    if (ListGetLogical(acc_solver % values, 'Eliminate Negative', found)) then
      acb_compfunc = max(0.0_dp, acb_compfunc)
    end if
  end subroutine acb_elim_neg

end module artificialcompressibility_acc
