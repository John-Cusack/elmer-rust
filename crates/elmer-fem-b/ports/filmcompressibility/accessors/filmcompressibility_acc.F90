! T022-A FilmCompressibility port: statement-level BIND(C) wrappers for
! the Elmer calls in fem/src/modules/FilmCompressibility.F90, compiled
! against the built .mod files of the target configure (A2 shipped
! mechanism, ADR-057 S3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, return codes). All Fortran-object traffic
! — TYPE(...) handles, OPTIONAL/POINTER dummies, assumed-shape
! descriptors, hidden lengths, list-directed timing reports — stays
! inside Fortran, so the Rust side is pure C calls. The wrapper module
! mirrors FilmCompressibility's own SAVE statement (STIFF, FORCE, Nodes,
! ElemPres, ElemHeight, ElemSens, Density, Thickness, EqPres,
! PrevElemHeight, PrevElemPres, ACave, ACinst, AllocationsDone,
! nVisited) with identical lifetimes, plus the entry locals the
! CONTAINS procedures host-associate (Found, GotEqPres, the mode flags,
! grav, Diff, SolvePDE, dp0, ac_min, the Volume/Area/PresInt/H3Int/LSMin
! /dpmax/dhmax/nlift accumulators), threaded in the exact original call
! order — including the author's stale-`Found` read after the
! Hydrostatic Pressure lookup, which replicates automatically.
!
! All floating-point work (the Galerkin assembly, MATMUL/SUM
! reductions, the WHERE relaxation update, the formatted Message
! writes) and all formatted I/O stay in Fortran (same compiler), so
! bit-identity holds by construction wherever the call sequence
! matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module filmcompressibility_acc
  use iso_c_binding
  use coordinatesystems
  use defutils
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(mesh_t), pointer, save :: acc_mesh => null()
  type(valuelist_t), pointer, save :: acc_params => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(variable_t), pointer, save :: acc_pressol => null()
  type(variable_t), pointer, save :: acc_heightsol => null()
  type(variable_t), pointer, save :: acc_senssol => null()
  type(variable_t), pointer, save :: acc_lsvar => null()
  real(kind=dp), pointer, save :: acc_acinst(:) => null()
  real(kind=dp), pointer, save :: acc_acave(:) => null()

  ! Shared SAVE Nodes state (consumed before the next fill, as in the
  ! original); per-internal IP state with identical lifetimes.
  type(nodes_t), target, save :: acc_nodes
  type(gaussintegrationpoints_t), save :: acc_ip_ls, acc_ip_lm

  ! Entry locals mirrored with identical lifetimes.
  logical, save :: acc_found = .false., acc_goteq = .false.
  logical, save :: acc_diffmode = .false., acc_sensmode = .false.
  logical, save :: acc_fullmode = .true., acc_solvepde = .false.
  logical, save :: acc_allocdone = .false., acc_lsinit = .false.
  integer, save :: acc_nlift = 0, acc_nvisited = 0, acc_dim = 0
  real(kind=dp), save :: acc_grav = 0.0_dp, acc_diff = 0.0_dp
  real(kind=dp), save :: acc_dp0 = 0.0_dp, acc_acmin = 0.0_dp
  real(kind=dp), save :: acc_volume = 0.0_dp, acc_area = 0.0_dp
  real(kind=dp), save :: acc_presint = 0.0_dp, acc_h3int = 0.0_dp
  real(kind=dp), save :: acc_lsmin = 0.0_dp
  real(kind=dp), save :: acc_dpmax = 0.0_dp, acc_dhmax = 0.0_dp

  ! Mirrors of the SAVE work arrays (first-mesh sizing,
  ! MaxElementNodes), plus internal-procedure locals sized likewise
  ! (only the leading 1:n sections are ever touched).
  real(kind=dp), allocatable, save :: acc_stiff(:,:), acc_force(:)
  real(kind=dp), allocatable, save :: acc_elemheight(:), acc_elempres(:)
  real(kind=dp), allocatable, save :: acc_prevelemheight(:)
  real(kind=dp), allocatable, save :: acc_prevelempres(:)
  real(kind=dp), allocatable, save :: acc_elemsens(:)
  real(kind=dp), allocatable, save :: acc_density(:), acc_thickness(:)
  real(kind=dp), allocatable, save :: acc_eqpres(:)
  real(kind=dp), allocatable, save :: acc_basis_ls(:), acc_elemls(:)
  real(kind=dp), allocatable, save :: acc_basis_lm(:), acc_dbasis_lm(:,:)

contains

  ! -- entry state ----------------------------------------------------
  subroutine fcm_begin(h_model, h_solver) bind(c, name="fcm_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    acc_mesh => acc_solver % mesh
    acc_params => GetSolverParams()
  end subroutine fcm_begin

  subroutine fcm_info_header() bind(c, name="fcm_info_header")
    call Info('FilmCompressibility',' ')
    call Info('FilmCompressibility','----------------------------------------------')
    call Info('FilmCompressibility','Solving compressibility field for FilmPressure')
    call Info('FilmCompressibility','----------------------------------------------')
  end subroutine fcm_info_header

  subroutine fcm_default_start() bind(c, name="fcm_default_start")
    call DefaultStart()
  end subroutine fcm_default_start

  subroutine fcm_check_transient(tr) bind(c, name="fcm_check_transient")
    integer(c_int), value :: tr
    if (tr == 0) then
      call Fatal('FilmCompressibility','Implemented only for transient problems!')
    end if
  end subroutine fcm_check_transient

  subroutine fcm_dim_touch() bind(c, name="fcm_dim_touch")
    acc_dim = CoordinateSystemDimension()
  end subroutine fcm_dim_touch

  subroutine fcm_assoc_acinst() bind(c, name="fcm_assoc_acinst")
    acc_acinst => acc_solver % variable % values
  end subroutine fcm_assoc_acinst

  subroutine fcm_read_modes() bind(c, name="fcm_read_modes")
    acc_diffmode = ListGetLogical(acc_params,'Differential Mode',acc_found)
    acc_sensmode = ListGetLogical(acc_params,'Sensitivity Mode',acc_found)
    acc_fullmode = .not. (acc_diffmode .or. acc_sensmode)
  end subroutine fcm_read_modes

  function fcm_sensmode() result(r) bind(c, name="fcm_sensmode")
    integer(c_int) :: r
    r = merge(1, 0, acc_sensmode)
  end function fcm_sensmode

  subroutine fcm_find_pressure() bind(c, name="fcm_find_pressure")
    character(len=MAX_NAME_LEN) :: varname
    varname = GetString(acc_params,'Pressure Variable Name',acc_found)
    if (.not. acc_found) varname = 'Pressure'
    acc_pressol => VariableGet(acc_mesh % variables, varname)
    if (.not. associated(acc_pressol)) then
      call Fatal('FilmCompressibility','"Pressure" variable '//trim(varname)//' not found!')
    end if
  end subroutine fcm_find_pressure

  subroutine fcm_find_height() bind(c, name="fcm_find_height")
    character(len=MAX_NAME_LEN) :: varname
    varname = GetString(acc_params,'Height Variable Name',acc_found)
    if (.not. acc_found) varname = 'Height'
    acc_heightsol => VariableGet(acc_mesh % variables, varname)
    if (.not. associated(acc_heightsol)) then
      call Fatal('FilmCompressibility','"Height" variable '//trim(varname)//' not found!')
    end if
  end subroutine fcm_find_height

  subroutine fcm_ranges20() bind(c, name="fcm_ranges20")
    if (InfoActive(20)) then
      call VectorValuesRange(acc_pressol % values,size(acc_pressol % values),acc_pressol % name)
      call VectorValuesRange(acc_heightsol % values,size(acc_heightsol % values),acc_heightsol % name)
    end if
  end subroutine fcm_ranges20

  subroutine fcm_sens_block() bind(c, name="fcm_sens_block")
    character(len=MAX_NAME_LEN) :: varname
    varname = ListGetString(acc_params,'Sensitivity Variable Name',acc_found)
    if (.not. acc_found) then
      call Fatal('FilmCompressibility','Give "Sensitivity Variable Name", e.g. "Displacement Adjoint 3"!')
    end if
    acc_senssol => VariableGet(acc_mesh % variables, varname)
    if (.not. associated(acc_senssol)) then
      call Fatal('FilmCompressibility','"Sensitivity" variable '//trim(varname)//' not found!')
    end if
    if (InfoActive(20)) then
      call VectorValuesRange(acc_senssol % values,size(acc_senssol % values),acc_senssol % name)
    end if
  end subroutine fcm_sens_block

  subroutine fcm_gravity() bind(c, name="fcm_gravity")
    real(kind=dp), pointer :: gwork(:,:)
    acc_grav = 0.0_dp
    if (ListGetLogical(acc_params,'Use Gravity',acc_found)) then
      gwork => ListGetConstRealArray(acc_model % constants,'Gravity',acc_found)
      if (acc_found) then
        acc_grav = abs(gwork(size(gwork,1),1))
      end if
    end if
  end subroutine fcm_gravity

  function fcm_fixcheck() result(r) bind(c, name="fcm_fixcheck")
    integer(c_int) :: r
    type(variable_t), pointer :: var
    integer :: i
    r = 0
    i = ListGetInteger(acc_params,'Fix AC After Iters',acc_found)
    if (i > 0) then
      var => VariableGet(acc_mesh % variables,'coupled iter')
      if (nint(var % values(1)) > i) r = 1
    end if
  end function fcm_fixcheck

  subroutine fcm_alloc() bind(c, name="fcm_alloc")
    type(variable_t), pointer :: var
    integer :: n, istat
    var => acc_solver % variable
    if (.not. acc_allocdone) then
      n = acc_mesh % maxelementnodes
      allocate(acc_force(n), acc_stiff(n,n), acc_elemheight(n), &
          acc_prevelemheight(n), acc_elemsens(n), acc_elempres(n), &
          acc_prevelempres(n), acc_density(n), acc_thickness(n), &
          acc_eqpres(n), acc_nodes % x(n), acc_nodes % y(n), &
          acc_nodes % z(n), acc_basis_ls(n), acc_elemls(n), &
          acc_basis_lm(n), acc_dbasis_lm(n,3), stat=istat)
      if (istat /= 0) call Fatal('FilmCompressibility','Memory allocation error.')
      acc_allocdone = .true.
      acc_lsinit = ListGetLogical(acc_params,'Levelset Initialize',acc_found)
    else
      acc_lsinit = .false.
    end if
  end subroutine fcm_alloc

  subroutine fcm_volreset() bind(c, name="fcm_volreset")
    acc_volume = 0.0_dp
    acc_area = 0.0_dp
    acc_presint = 0.0_dp
  end subroutine fcm_volreset

  function fcm_lsinit() result(r) bind(c, name="fcm_lsinit")
    integer(c_int) :: r
    r = merge(1, 0, acc_lsinit)
  end function fcm_lsinit

  function fcm_nof_active() result(n) bind(c, name="fcm_nof_active")
    integer(c_int) :: n
    n = acc_solver % numberofactiveelements
  end function fcm_nof_active

  ! -- LevelsetInitialize path (LocalLS + SetACInit) --------------------
  subroutine fcm_ls_begin() bind(c, name="fcm_ls_begin")
    character(len=MAX_NAME_LEN) :: varname
    varname = GetString(acc_params,'Levelset Variable Name',acc_found)
    acc_lsvar => VariableGet(acc_mesh % variables, varname, UnfoundFatal=.true.)
    acc_h3int = 0.0_dp
    acc_lsmin = 0.0_dp
  end subroutine fcm_ls_begin

  subroutine fcm_localls(t) bind(c, name="fcm_localls")
    integer(c_int), value :: t
    real(kind=dp) :: detj, s, ls, h
    logical :: stat
    integer :: n, q
    type(valuelist_t), pointer :: material
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes()
    acc_ip_ls = GaussPoints(acc_elem)
    acc_nodes % x(1:n) = acc_model % nodes % x(acc_elem % nodeindexes)
    acc_nodes % y(1:n) = acc_model % nodes % y(acc_elem % nodeindexes)
    acc_nodes % z(1:n) = 0_dp
    call GetScalarLocalSolution(acc_elemls, UVariable=acc_lsvar)
    material => GetMaterial(acc_elem)
    acc_thickness(1:n) = GetReal(material,'Thickness')
    do q=1,acc_ip_ls % n
      stat = ElementInfo(acc_elem, acc_nodes, acc_ip_ls % u(q), &
          acc_ip_ls % v(q), acc_ip_ls % w(q), detj, acc_basis_ls(1:n))
      s = acc_ip_ls % s(q) * detj
      ls = sum(acc_basis_ls(1:n) * acc_elemls(1:n))
      if (ls < 0.0) then
        acc_area = acc_area + s
        h = sum(acc_basis_ls(1:n) * acc_thickness(1:n))
        acc_h3int = acc_h3int + s / h**3
        acc_lsmin = min(acc_lsmin, ls)
      end if
    end do
  end subroutine fcm_localls

  subroutine fcm_ls_reduce() bind(c, name="fcm_ls_reduce")
    acc_area = ParallelReduction(acc_area)
    acc_h3int = ParallelReduction(acc_h3int)
    acc_lsmin = ParallelReduction(acc_lsmin,1)
  end subroutine fcm_ls_reduce

  subroutine fcm_setacinit() bind(c, name="fcm_setacinit")
    real(kind=dp) :: youngs, pois, d, reff, ls, h, heff
    integer :: i, j
    type(valuelist_t), pointer :: material
    do i=1,acc_model % numberofmaterials
      material => acc_model % materials(i) % values
      youngs = ListGetCReal(material,'Youngs Modulus',acc_found)
      pois = ListGetCReal(material,'Poisson Ratio',acc_found)
      if (acc_found) exit
    end do
    heff = (acc_h3int/acc_area)**(-1.0/3)
    d = youngs * heff**3 / (12 * (1-pois**2))
#if 1
    reff = sqrt(acc_area / PI)
#else
    reff = -acc_lsmin
#endif
    do i=1,acc_mesh % numberofnodes
      j = acc_solver % variable % perm(i)
      if (j == 0) cycle
      ls = acc_lsvar % values(acc_lsvar % perm(i))
      if (ls > 0.0_dp) then
        h = 0.0_dp
      else
        h = (reff**4 / 64 * d) * (1 - (ls/acc_lsmin)**2)**2
      end if
      acc_acave(j) = h
    end do
    call Info('FilmCompressibility','Initialized the artificial compressibility field!',level=6)
  end subroutine fcm_setacinit

  ! -- maxima sweep (LocalMaxima) ---------------------------------------
  subroutine fcm_dpmax_reset() bind(c, name="fcm_dpmax_reset")
    acc_dpmax = 0.0_dp
    acc_dhmax = 0.0_dp
    acc_nlift = 0
  end subroutine fcm_dpmax_reset

  subroutine fcm_localmaxima(t) bind(c, name="fcm_localmaxima")
    integer(c_int), value :: t
    real(kind=dp) :: p0, p1, dpres, h0, h1, dh, th, rho
    integer :: n, q, nt
    type(valuelist_t), pointer :: material
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes()
    call GetScalarLocalSolution(acc_elemheight, tStep=0, UVariable=acc_heightsol)
    call GetScalarLocalSolution(acc_elempres, tStep=0, UVariable=acc_pressol)
    if (acc_diffmode) then
      call GetScalarLocalSolution(acc_prevelemheight, tStep=-1, UVariable=acc_heightsol)
      call GetScalarLocalSolution(acc_prevelempres, tStep=-1, UVariable=acc_pressol)
    end if
    material => GetMaterial(acc_elem)
    acc_eqpres(1:n) = GetReal(material,'Hydrostatic Pressure',acc_goteq)
    if (.not. acc_found) then
      acc_density(1:n) = GetReal(material,'Ice Density')
      acc_thickness(1:n) = GetReal(material,'Thickness')
    end if
    nt = 0
    do q=1,n
      if (acc_fullmode) then
        h0 = 0.0_dp
        if (acc_goteq) then
          p0 = acc_eqpres(q)
        else
          rho = acc_density(q)
          th = acc_thickness(q)
          p0 = rho * acc_grav * th
        end if
      else
        h0 = acc_prevelemheight(q)
        p0 = acc_prevelempres(q)
      end if
      h1 = acc_elemheight(q)
      dh = h1 - h0
      p1 = acc_elempres(q)
      dpres = p1 - p0
      acc_dhmax = max(acc_dhmax,h1)
      acc_dpmax = max(acc_dpmax,dpres)
      if (h1 * dpres > 1.0e-20) nt = nt + 1
    end do
    if (nt == n) acc_nlift = acc_nlift + 1
  end subroutine fcm_localmaxima

  subroutine fcm_max_reduce() bind(c, name="fcm_max_reduce")
    acc_dpmax = ParallelReduction(acc_dpmax, 2)
    acc_dhmax = ParallelReduction(acc_dhmax, 2)
  end subroutine fcm_max_reduce

  subroutine fcm_max_report() bind(c, name="fcm_max_report")
    write(message,'(A,ES12.3)') 'Maximum pressure overhead:', acc_dpmax
    call Info('FilmCompressibility',message,level=10)
    write(message,'(A,ES12.3)') 'Maximum displacement:', acc_dhmax
    call Info('FilmCompressibility',message,level=10)
    call Info('FilmCompressibility','Number of lifted elements: '//I2S(acc_nlift),level=10)
  end subroutine fcm_max_report

  ! -- coefficients -------------------------------------------------------
  subroutine fcm_read_diffusion() bind(c, name="fcm_read_diffusion")
    acc_diff = ListGetCReal(acc_params,'AC Diffusion',acc_solvepde)
  end subroutine fcm_read_diffusion

  subroutine fcm_read_coeffs() bind(c, name="fcm_read_coeffs")
    real(kind=dp) :: val
    acc_acmin = ListGetCReal(acc_params,'AC Minimum Value',acc_found)
    acc_dp0 = ListGetCReal(acc_params,'AC Pressure Epsilon',acc_found)
    if (.not. acc_found) then
      val = ListGetCReal(acc_params,'AC Pressure Epsilon Relative',acc_found)
      if (.not. acc_found) val = 0.01_dp
      acc_dp0 = max(epsilon(acc_dp0),val*acc_dpmax)
    end if
  end subroutine fcm_read_coeffs

  ! -- assembly driver ------------------------------------------------------
  subroutine fcm_default_initialize() bind(c, name="fcm_default_initialize")
    call DefaultInitialize()
  end subroutine fcm_default_initialize

  subroutine fcm_assembly_elem(t) bind(c, name="fcm_assembly_elem")
    integer(c_int), value :: t
    real(kind=dp) :: detj, s, p0, p1, dpres, h1, dh, th, ac, rho
    logical :: stat
    integer :: n, q, p, qq, dim
    type(valuelist_t), pointer :: material
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes()
    acc_stiff = 0.0_dp
    acc_force = 0.0_dp
    dim = CoordinateSystemDimension()
    acc_ip_lm = GaussPoints(acc_elem)
    acc_nodes % x(1:n) = acc_model % nodes % x(acc_elem % nodeindexes)
    acc_nodes % y(1:n) = acc_model % nodes % y(acc_elem % nodeindexes)
    acc_nodes % z(1:n) = 0_dp
    call GetScalarLocalSolution(acc_elemheight, tStep=0, UVariable=acc_heightsol)
    call GetScalarLocalSolution(acc_elempres, tStep=0, UVariable=acc_pressol)
    if (acc_diffmode) then
      call GetScalarLocalSolution(acc_prevelemheight, tStep=-1, UVariable=acc_heightsol)
      call GetScalarLocalSolution(acc_prevelempres, tStep=-1, UVariable=acc_pressol)
    end if
    if (acc_sensmode) then
      call GetScalarLocalSolution(acc_elemsens, UVariable=acc_senssol)
    end if
    material => GetMaterial(acc_elem)
    if (acc_fullmode) then
      acc_eqpres(1:n) = GetReal(material,'Hydrostatic Pressure',acc_goteq)
      if (.not. acc_found) then
        acc_density(1:n) = GetReal(material,'Ice Density')
        acc_thickness(1:n) = GetReal(material,'Thickness')
      end if
    end if
    do q=1,acc_ip_lm % n
      stat = ElementInfo(acc_elem, acc_nodes, acc_ip_lm % u(q), &
          acc_ip_lm % v(q), acc_ip_lm % w(q), detj, &
          acc_basis_lm(1:n), acc_dbasis_lm(1:n,1:3))
      s = acc_ip_lm % s(q) * detj
      h1 = sum(acc_elemheight(1:n) * acc_basis_lm(1:n))
      p1 = sum(acc_elempres(1:n) * acc_basis_lm(1:n))
      if (acc_goteq) then
        p0 = sum(acc_basis_lm(1:n) * acc_eqpres(1:n))
      else
        rho = sum(acc_basis_lm(1:n) * acc_density(1:n))
        th = sum(acc_basis_lm(1:n) * acc_thickness(1:n))
        p0 = rho * acc_grav * th
      end if
      if (acc_sensmode) then
        dh = h1
        dpres = p1 - p0
        ac = sum(acc_elemsens(1:n) * acc_basis_lm(1:n))
      else if (acc_fullmode) then
        dh = h1
        dpres = p1 - p0
        ac = dh / max(dpres, acc_dp0)
      else
        dh = h1 - sum(acc_prevelemheight(1:n) * acc_basis_lm(1:n))
        dpres = p1 - sum(acc_prevelempres(1:n) * acc_basis_lm(1:n))
        ac = dh / max(dpres, acc_dp0)
      end if
      if (acc_solvepde) then
        acc_stiff(1:n,1:n) = acc_stiff(1:n,1:n) + &
            s * acc_diff * matmul(acc_dbasis_lm(1:n,1:3), &
            transpose(acc_dbasis_lm(1:n,1:3)))
        do p=1,n
          do qq=1,n
            acc_stiff(p,qq) = acc_stiff(p,qq) + s * acc_basis_lm(p) * acc_basis_lm(qq)
          end do
        end do
      else
        do p=1,n
          acc_stiff(p,p) = acc_stiff(p,p) + s * acc_basis_lm(p)
        end do
      end if
      acc_force(1:n) = acc_force(1:n) + s * ac
      if (h1 > epsilon(h1)) then
        acc_volume = acc_volume + s * dh
        acc_area = acc_area + s
        acc_presint = acc_presint + s * dpres
      end if
    end do
    call DefaultUpdateEquations(acc_stiff, acc_force)
  end subroutine fcm_assembly_elem

  ! -- shared assembly tail ---------------------------------------------------
  subroutine fcm_finish_boundary() bind(c, name="fcm_finish_boundary")
    call DefaultFinishBoundaryAssembly()
  end subroutine fcm_finish_boundary

  subroutine fcm_finish_assembly() bind(c, name="fcm_finish_assembly")
    call DefaultFinishAssembly()
  end subroutine fcm_finish_assembly

  subroutine fcm_dirichlet() bind(c, name="fcm_dirichlet")
    call DefaultDirichletBCs()
  end subroutine fcm_dirichlet

  function fcm_solve() result(norm) bind(c, name="fcm_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function fcm_solve

  subroutine fcm_default_finish() bind(c, name="fcm_default_finish")
    call DefaultFinish()
  end subroutine fcm_default_finish

  subroutine fcm_final_reduce() bind(c, name="fcm_final_reduce")
    acc_volume = ParallelReduction(acc_volume)
    acc_area = ParallelReduction(acc_area)
    acc_presint = ParallelReduction(acc_presint)
  end subroutine fcm_final_reduce

  subroutine fcm_acinst_range() bind(c, name="fcm_acinst_range")
    call VectorValuesRange(acc_acinst,size(acc_acinst),acc_solver % variable % name)
  end subroutine fcm_acinst_range

  subroutine fcm_report() bind(c, name="fcm_report")
    real(kind=dp) :: val
    call ListAddConstReal(acc_model % simulation,'res: AC Volume',acc_volume)
    write(message,'(A,ES12.5)') 'AC Volume: ',acc_volume
    call Info('FilmCompressibility',message)
    call ListAddConstReal(acc_model % simulation,'res: AC Area',acc_area)
    write(message,'(A,ES12.5)') 'AC Area: ',acc_area
    call Info('FilmCompressibility',message)
    val = 0.0_dp
    if (acc_presint > epsilon(val)) val = acc_volume / acc_presint
    call ListAddConstReal(acc_model % simulation,'res: AC Average',val)
    write(message,'(A,ES12.5)') 'AC Average: ',val
    call Info('FilmCompressibility',message)
  end subroutine fcm_report

  subroutine fcm_update() bind(c, name="fcm_update")
    type(variable_t), pointer :: var
    real(kind=dp) :: relaxp, relaxm, relax, val
    logical :: doit
    var => VariableGet(acc_mesh % variables,'AC')
    acc_acave => var % values
    doit = .true.
    val = ListGetCReal(acc_params,'AC Update Condition',acc_found)
    if (acc_found .and. val < 0.0_dp) then
      call Info('FilmCompressibility','Not updating the AC field because of condition!',level=7)
      doit = .false.
    end if
    val = ListGetCReal(acc_params,'AC Critical Volume',acc_found)
    if (acc_found .and. acc_volume < val) then
      call Info('FilmCompressibility','Not updating the AC field because of small volume!',level=7)
      doit = .false.
    end if
    val = ListGetCReal(acc_params,'AC Critical Area',acc_found)
    if (acc_found .and. acc_area < val) then
      call Info('FilmCompressibility','Not updating the AC field because of small area!',level=7)
      doit = .false.
    end if
    if (doit) then
      call Info('FilmCompressibility','Updating AC field for FSI coupling!',level=7)
      relax = ListGetCReal(acc_params,'AC Relaxation Factor',acc_found)
      if (.not. acc_found) relax = 1.0_dp
      relaxm = ListGetCReal(acc_params,'AC Relaxation Factor Negative',acc_found)
      if (.not. acc_found) relaxm = relax
      relaxp = ListGetCReal(acc_params,'AC Relaxation Factor Positive',acc_found)
      if (.not. acc_found) relaxp = relax
      if (acc_nvisited == 0) then
        relaxm = 1.0_dp
        relaxp = 1.0_dp
      end if
      val = ListGetCReal(acc_params,'AC Multiplier',acc_found)
      if (.not. acc_found) val = 1.0_dp
      if (relaxm*relaxp < 1.0_dp - epsilon(relax)) then
        where (val * acc_acinst > acc_acave)
          acc_acave = acc_acave + relaxp * (val*acc_acinst - acc_acave)
        elsewhere
          acc_acave = acc_acave + relaxm * (val*acc_acinst - acc_acave)
        end where
      else
        acc_acave = val * acc_acinst
      end if
      acc_acave = max(acc_acave, acc_acmin)
    end if
    if (InfoActive(15)) then
      call VectorValuesRange(acc_acave,size(acc_acave),var % name)
    end if
    acc_nvisited = acc_nvisited + 1
  end subroutine fcm_update

  subroutine fcm_alldone_info() bind(c, name="fcm_alldone_info")
    call Info('FilmCompressibility','All done for now!',level=12)
  end subroutine fcm_alldone_info

  ! -- FilmCompressibility_Init body --------------------------------------------
  subroutine fcm_init_begin(h_solver) bind(c, name="fcm_init_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine fcm_init_begin

  subroutine fcm_init_body() bind(c, name="fcm_init_body")
    type(valuelist_t), pointer :: solverpars
    solverpars => GetSolverParams()
    call ListAddNewInteger(solverpars,'Time Derivative Order',0)
    call ListAddNewString(solverpars,'Variable','AC inst')
    call ListAddString(solverpars,NextFreeKeyword('Exported Variable ',solverpars),'AC')
  end subroutine fcm_init_body

end module filmcompressibility_acc
