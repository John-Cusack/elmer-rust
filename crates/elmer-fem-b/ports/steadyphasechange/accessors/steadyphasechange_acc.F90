! T021-B SteadyPhaseChange port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/SteadyPhaseChange.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls.
!
! The wrapper module mirrors the original's SAVE heap state with identical
! lifetimes: PrevTemp/NodeDone/Nodes/x/y/z/Basis/NodalTemp/TempDiff/
! BoundaryMarker/PhaseElements (allocated once on the first call, sized by
! first-mesh data exactly as in the original), IsoSurf (+flag, allocated
! and freed around each CreateIsotherm pass pair), and the MeltingHeat
! SAVE state (Initialized/Nodes/NormalExist/NormalSol). Plain SAVE
! scalars that only feed branches (Newton) live Rust-side; every other
! SAVE scalar (MeltPoint, Width, Trip_node, ccum, prevtave, ...) is a
! module scalar read back only inside wrappers. Every floating-point
! update, reduction, and formatted I/O statement below is an exact
! transcription of the original (same compiler, same evaluation order),
! so bit-identity holds by construction wherever the wrapper call
! sequence matches.
!
! Dead SAVE weight is dropped, never faked: ForceVector/FluxCorrect/
! NodalTemp-values/Weights/SurfaceVelocitySet/Visited/Eps-host/Norm-host/
! pos0/prevpos0/maxds/SpeedUp/OrigNorm/PrevNorm/MaxAngle-locals and the
! unused `str` are declared-or-SAVEd in the original but never read on
! any path that affects output (Eps/Norm/MaxAngle are recomputed per use;
! NodalTemp is allocated but never filled — the allocation itself is
! kept for Fatal parity).
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module steadyphasechange_acc
  use iso_c_binding
  use types
  use defutils
  use messages
  use solverbasics
  use elementdescription
  implicit none

  ! -- entry handles (re-pointed on every call, as in the original) -----
  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(solver_t), pointer, save :: acc_psolver => null()
  type(valuelist_t), pointer, save :: acc_params => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(variable_t), pointer, save :: acc_surfsol => null()
  type(variable_t), pointer, save :: acc_tempsol => null()
  type(variable_t), pointer, save :: acc_helpsol => null()
  real(kind=dp), pointer, save :: acc_surface(:) => null()
  real(kind=dp), pointer, save :: acc_surfacemove(:) => null()
  real(kind=dp), pointer, save :: acc_temperature(:) => null()
  integer, pointer, save :: acc_surfperm(:) => null()
  integer, pointer, save :: acc_tempperm(:) => null()
  integer, pointer, save :: acc_indexes(:) => null()
  character(len=MAX_NAME_LEN), save :: acc_surfname = ''

  ! -- SAVE heap mirrors (allocated once on the first call) -------------
  real(kind=dp), allocatable, save :: acc_prevtemp(:)
  logical, allocatable, save :: acc_nodedone(:)
  type(nodes_t), target, save :: acc_nodes
  real(kind=dp), allocatable, save :: acc_x(:), acc_y(:), acc_z(:), &
       acc_basis(:), acc_nodaltemp(:), acc_tempdiff(:)
  logical, allocatable, save :: acc_boundarymarker(:)
  integer, allocatable, save :: acc_phaseelements(:)
  logical, save :: acc_allocationsdone = .false.

  ! -- first-call SAVE scalars (set by sp_firsttime) --------------------
  integer, save :: acc_normaldir = 1, acc_tangentdir = 1
  real(kind=dp), save :: acc_meltpoint = 0.0_dp
  real(kind=dp), save :: acc_width = 0.0_dp
  integer, save :: acc_tripnode = -1, acc_axisnode = -1
  integer, save :: acc_nobnodes = 0, acc_nophaseelements = 0
  logical, save :: acc_gotliquid = .false., acc_gotsolid = .false.
  integer, save :: acc_liquidbody = 0, acc_solidbody = 0
  real(kind=dp), save :: acc_coordmax(3) = 0.0_dp, acc_coordmin(3) = 0.0_dp
  integer, save :: acc_coordmaxi(3) = -1, acc_coordmini(3) = -1

  ! -- per-call accumulators + lumped-acceleration SAVE ------------------
  real(kind=dp), save :: acc_area = 0.0_dp, acc_volume = 0.0_dp
  real(kind=dp), save :: acc_tave = 0.0_dp, acc_tabs = 0.0_dp
  real(kind=dp), save :: acc_volabs = 0.0_dp
  real(kind=dp), save :: acc_maxtempdiff = 0.0_dp
  real(kind=dp), save :: acc_dpos = 0.0_dp
  real(kind=dp), save :: acc_ccum = 1.0_dp
  real(kind=dp), save :: acc_prevtave = 0.0_dp, acc_prevtabs = 0.0_dp
  real(kind=dp), save :: acc_prevvolabs = 0.0_dp, acc_prevvolume = 0.0_dp
  real(kind=dp), save :: acc_newtontol = 0.0_dp
  logical, save :: acc_triplefixed = .false.

  ! -- CreateIsotherm SAVE (allocated/freed around each call pair) ------
  real(kind=dp), allocatable, save :: acc_isosurf(:,:)
  logical, save :: acc_isosurfallocated = .false.
  integer, save :: acc_nelems = 0
  real(kind=dp), save :: acc_xmin = 0.0_dp, acc_xmax = 0.0_dp

  ! -- MeltingHeat SAVE (serial proof tree: single slot, tid always 1) --
  logical, save :: acc_mh_initialized = .false.
  type(nodes_t), allocatable, save :: acc_mh_nodes(:)
  logical, save :: acc_mh_normalexist = .false.
  type(variable_t), pointer, save :: acc_mh_normalsol => null()

contains

  ! -- entry state ------------------------------------------------------
  subroutine sp_begin(h_model, h_solver) bind(c, name="sp_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    acc_psolver => acc_solver
    acc_params => GetSolverParams()
  end subroutine sp_begin

  subroutine sp_info_start() bind(c, name="sp_info_start")
    call Info('SteadyPhaseChange',   '--------------------------------------------')
    call Info('SteadyPhaseChange',   'Using steady algorithm to find the isotherm')
    call Info('SteadyPhaseChange',   '--------------------------------------------')
  end subroutine sp_info_start

  function sp_dim() result(n) bind(c, name="sp_dim")
    integer(c_int) :: n
    n = CoordinateSystemDimension()
  end function sp_dim

  subroutine sp_fatal_2d() bind(c, name="sp_fatal_2d")
    call Fatal('SteadyPhaseChange','Implemented only in 2D')
  end subroutine sp_fatal_2d

  subroutine sp_surface_begin() bind(c, name="sp_surface_begin")
    character(len=MAX_NAME_LEN) :: variablename
    logical :: stat
    acc_surfsol => acc_solver % variable
    if (associated(acc_surfsol)) then
       acc_surface => acc_surfsol % values
       acc_surfperm => acc_surfsol % perm
       if (.not. associated(acc_surface) .or. all(acc_surfperm <= 0)) then
          call Fatal('SteadyPhaseChange','Surface field needed for Phase Change')
       end if
       variablename = ComponentName(acc_solver % variable)
    else
       variablename = ListGetString(acc_params, 'Surface Variable', stat)
       if (.not. stat) variablename = 'Surface'
       acc_surfsol => VariableGet(acc_solver % mesh % variables, variablename)
       acc_surface => acc_surfsol % values
       acc_surfperm => acc_surfsol % perm
       if (.not. associated(acc_surface) .or. all(acc_surfperm <= 0)) then
          call Fatal('SteadyPhaseChange','Surface field needed for Phase Change')
       end if
    end if
    acc_surfname = variablename
  end subroutine sp_surface_begin

  subroutine sp_surfacemove_begin() bind(c, name="sp_surfacemove_begin")
    integer :: istat
    acc_helpsol => VariableGet(acc_solver % mesh % variables, &
         trim(acc_surfname)//'Diff')
    if (associated(acc_helpsol)) then
       acc_surfacemove => acc_helpsol % values
    else
       allocate(acc_surfacemove(size(acc_surface)), stat=istat)
       if (istat /= 0) call Fatal('PhaseChangeSolver','Memory allocation error 1.')
       acc_surfacemove = 0.0d0
       call VariableAdd(acc_solver % mesh % variables, acc_solver % mesh, &
            acc_psolver, trim(acc_surfname)//'Diff', 1, acc_surfacemove, acc_surfperm)
    end if
  end subroutine sp_surfacemove_begin

  subroutine sp_temperature_begin() bind(c, name="sp_temperature_begin")
    character(len=MAX_NAME_LEN) :: variablename
    logical :: stat
    variablename = ListGetString(acc_params, 'Phase Change Variable', stat)
    if (.not. stat) variablename = 'Temperature'
    acc_tempsol => VariableGet(acc_solver % mesh % variables, variablename)
    acc_tempperm => acc_tempsol % perm
    acc_temperature => acc_tempsol % values
    if (.not. associated(acc_temperature)) then
       call Fatal('SteadyPhaseChange','Temperature field needed for Phase Change')
    end if
    if (all(acc_tempperm <= 0)) then
       call Fatal('SteadyPhaseChange','Temperature field has zero Perm vector!')
    end if
  end subroutine sp_temperature_begin

  function sp_get_relax() result(r) bind(c, name="sp_get_relax")
    real(c_double) :: r
    logical :: stat
    r = GetCReal(acc_params, 'Nonlinear System Relaxation Factor', stat)
    if (.not. stat) r = 1.0d0
  end function sp_get_relax

  ! -- first call: mesh sizing and permanent vectors ----------------------
  subroutine sp_firsttime() bind(c, name="sp_firsttime")
    type(element_t), pointer :: element, parent, parent2
    type(valuelist_t), pointer :: bc
    integer :: t, n, i, j, k, istat, dim
    real(kind=dp) :: xx
    logical :: gotit
    dim = CoordinateSystemDimension()
    allocate(acc_boundarymarker(size(acc_surfperm)), stat=istat)
    if (istat /= 0) call Fatal('PhaseChangeSolver','Memory allocation error 2.')
    acc_boundarymarker = .false.
    acc_nophaseelements = 0
    do t = 1, acc_solver % mesh % numberofboundaryelements
       element => GetBoundaryElement(t)
       n = GetElementNOFNodes()
       if (GetElementFamily() == 1) cycle
       bc => GetBC()
       if (.not. GetLogical(bc, 'Phase Change', gotit)) cycle
       acc_boundarymarker(element % nodeindexes) = .true.
       acc_nophaseelements = acc_nophaseelements + 1
    end do
    acc_nobnodes = count(acc_boundarymarker)
    write(message,'(A,T35,I12)') 'Number of interface nodes:',acc_nobnodes
    call Info('SteadyPhaseChange',message)
    if (acc_nobnodes == 0) then
       call Warn('SteadyPhaseChange','There is no Phase Change boundary')
    end if
    write(message,'(A,T35,I12)') 'Number of interface elements:',acc_nophaseelements
    call Info('SteadyPhaseChange',message)
    allocate(acc_phaseelements(acc_nophaseelements), stat=istat)
    if (istat /= 0) call Fatal('PhaseChangeSolver','Memory allocation error 3.')
    acc_phaseelements = 0
    acc_nophaseelements = 0
    do t = 1, acc_solver % mesh % numberofboundaryelements
       element => GetBoundaryElement(t)
       n = GetElementNOFNodes()
       if (GetElementFamily() == 1) cycle
       bc => GetBC()
       if (.not. GetLogical(bc, 'Phase Change', gotit)) cycle
       acc_nophaseelements = acc_nophaseelements + 1
       acc_phaseelements(acc_nophaseelements) = t
    end do
    acc_coordmax = -huge(acc_coordmax)
    acc_coordmin = huge(acc_coordmin)
    do k = 1, acc_model % mesh % numberofnodes
       if (.not. acc_boundarymarker(k)) cycle
       do j = 1, dim
          if (j == 1) xx = acc_model % mesh % nodes % x(k)
          if (j == 2) xx = acc_model % mesh % nodes % y(k)
          if (j == 3) xx = acc_model % mesh % nodes % z(k)
          if (xx > acc_coordmax(j)) then
             acc_coordmax(j) = xx
             acc_coordmaxi(j) = k
          end if
          if (xx < acc_coordmin(j)) then
             acc_coordmin(j) = xx
             acc_coordmini(j) = k
          end if
       end do
    end do
    j = 1
    do i = 1, dim
       if (acc_coordmax(i)-acc_coordmin(i) < acc_coordmax(j)-acc_coordmin(j)) j = i
    end do
    acc_normaldir = j
    j = 1
    do i = 1, dim
       if (acc_coordmax(i)-acc_coordmin(i) > acc_coordmax(j)-acc_coordmin(j)) j = i
    end do
    acc_width = acc_coordmax(j) - acc_coordmin(j)
    acc_tangentdir = j
    write(message,'(A,T35,ES12.4)') 'Width of the interface',acc_width
    call Info('SteadyPhaseChange',message)
    acc_tripnode = acc_coordmaxi(acc_tangentdir)
    acc_axisnode = acc_coordmini(acc_tangentdir)
    write(message,'(A,T35,I12)') 'Index of the triple point: ',acc_tripnode
    call Info('SteadyPhaseChange',message)
    write(message,'(A,T35,I12)') 'Index of the axis point: ',acc_axisnode
    call Info('SteadyPhaseChange',message)
    acc_gotliquid = .false.
    acc_gotsolid = .false.
    do i = 1, acc_model % numberofmaterials
       acc_gotliquid = acc_gotliquid .and. &
            ListGetLogical(acc_model % bodies(i) % values,'Liquid', gotit)
       acc_gotsolid = acc_gotsolid .and. &
            ListGetLogical(acc_model % bodies(i) % values,'Solid', gotit)
    end do
    acc_solidbody = 0
    acc_liquidbody = 0
    if (.not. (acc_gotliquid .and. acc_gotsolid)) then
       do t = 1, acc_solver % mesh % numberofboundaryelements
          element => GetBoundaryElement(t)
          n = GetElementNOFNodes()
          if (GetElementFamily() == 1) cycle
          bc => GetBC()
          if (.not. GetLogical(bc, 'Phase Change', gotit)) cycle
          parent => element % boundaryinfo % left
          parent2 => element % boundaryinfo % right
          if (maxval(acc_model % mesh % nodes % y(parent % nodeindexes)) > &
               maxval(acc_model % mesh % nodes % y(parent2 % nodeindexes))) then
             if (.not. acc_gotsolid) acc_solidbody = parent % bodyid
             if (.not. acc_gotliquid) acc_liquidbody = parent2 % bodyid
          else
             if (.not. acc_gotsolid) acc_solidbody = parent2 % bodyid
             if (.not. acc_gotliquid) acc_liquidbody = parent % bodyid
          end if
          exit
       end do
       write(message,'(A,T35,I12)') 'Body index for solid: ',acc_solidbody
       call Info('SteadyPhaseChange',message)
       write(message,'(A,T35,I12)') 'Body index for liquid: ',acc_liquidbody
       call Info('SteadyPhaseChange',message)
    end if
    n = acc_solver % mesh % maxelementnodes
    allocate(acc_nodes % x(n), acc_nodes % y(n), acc_nodes % z(n), &
         acc_x(n), acc_y(n), acc_z(n), acc_basis(n), acc_nodaltemp(n), &
         acc_tempdiff(n), acc_prevtemp(size(acc_surface)), stat=istat)
    if (istat /= 0) call Fatal('SteadyPhaseChange','Memory allocation error 4.')
    acc_nodes % x = 0.0d0
    acc_nodes % y = 0.0d0
    acc_nodes % z = 0.0d0
    acc_prevtemp = 0.0d0
    allocate(acc_nodedone(size(acc_surface)), stat=istat)
    if (istat /= 0) call Fatal('SteadyPhaseChange','Memory allocation error 5.')
    acc_allocationsdone = .true.
  end subroutine sp_firsttime
  ! -- per-call scalars ---------------------------------------------------
  function sp_trip_temp() result(r) bind(c, name="sp_trip_temp")
    real(c_double) :: r
    r = acc_temperature(acc_tempperm(acc_tripnode))
  end function sp_trip_temp

  function sp_passive_steps() result(i) bind(c, name="sp_passive_steps")
    integer(c_int) :: i
    logical :: stat
    i = ListGetInteger(acc_params,'Passive Steps',stat)
  end function sp_passive_steps

  function sp_newton_update(visit, newton_in) result(nout) bind(c, name="sp_newton_update")
    integer(c_int), value :: visit, newton_in
    integer(c_int) :: nout
    integer :: iter
    logical :: stat, newton
    newton = (newton_in /= 0)
    iter = ListGetInteger(acc_params, 'Nonlinear System Newton After Iterations', stat)
    if (stat .and. visit > iter) newton = .true.
    acc_newtontol = ListGetConstReal(acc_params, &
         'Nonlinear System Newton After Tolerance', stat)
    if (newton) then
       call Info('SteadyPhaseChange','Steady state newton formulation', level=4)
    else
       call Info('SteadyPhaseChange','Steady state isoterm formulation', level=4)
    end if
    nout = merge(1, 0, newton)
  end function sp_newton_update

  subroutine sp_meltpoint(trip_temp) bind(c, name="sp_meltpoint")
    real(c_double), value :: trip_temp
    integer :: k
    logical :: stat, gotit
    acc_triplefixed = ListGetLogical(acc_params,'Triple Point Fixed',stat)
    if (.not. stat) acc_triplefixed = .false.
    if (acc_triplefixed) then
       acc_meltpoint = trip_temp
       write(message,'(A,T35,ES12.4)') 'Melting point set: ',acc_meltpoint
       call Info('SteadyPhaseChange',message)
    else
       do k = 1, acc_model % numberofmaterials
          acc_meltpoint = GetCReal(acc_model % materials(k) % values, &
               'Melting Point', gotit)
          if (gotit) exit
       end do
       if (gotit) then
          write(message,'(A,T35,ES12.4)') 'Melting point found: ',acc_meltpoint
          call Info('SteadyPhaseChange',message)
       else
          call Info('SteadyPhaseChange','Could not find melting point in any material!')
       end if
    end if
  end subroutine sp_meltpoint
  ! -- CreateIsotherm internal procedure (exact transcription, incl. ------
  ! -- the two-pass GOTO 100 counting/allocation scheme) ------------------
  subroutine sp_createisotherm() bind(c, name="sp_createisotherm")
    type(element_t), pointer :: element
    integer, pointer :: indexes(:)
    integer :: t, elementcode, k, n, vertex, nn, next, i, imin, nelems
    real(kind=dp) :: temp1, temp2, temppi
    logical :: active, stat
    acc_isosurfallocated = .false.
    acc_xmin = huge(acc_xmin)
    acc_xmax = -huge(acc_xmax)
100 nelems = 0
    do t = 1, acc_solver % mesh % numberofbulkelements
       element => acc_solver % mesh % elements(t)
       elementcode = element % type % elementcode
       if (elementcode < 300) cycle
       k = ListGetInteger(acc_model % bodies(element % bodyid) % values,'Material')
       active = .false.
       if (acc_gotliquid) then
          if (ListGetLogical(acc_model % materials(k) % values, 'Liquid', stat)) &
               active = .true.
       else
          if (element % bodyid == acc_liquidbody) active = .true.
       end if
       if (acc_gotsolid) then
          if (ListGetLogical(acc_model % materials(k) % values, 'Solid', stat)) &
               active = .true.
       else
          if (element % bodyid == acc_solidbody) active = .true.
       end if
       if (.not. active) cycle
       n = element % type % numberofnodes
       indexes => element % nodeindexes
       if (any(acc_tempperm(indexes(1:n)) <= 0)) cycle
       acc_tempdiff(1:n) = acc_temperature(acc_tempperm(indexes(1:n))) - acc_meltpoint
       if (all(acc_tempdiff(1:n) < 0.0)) cycle
       if (all(acc_tempdiff(1:n) > 0.0)) cycle
       acc_nodes % x(1:n) = acc_solver % mesh % nodes % x(indexes)
       acc_nodes % y(1:n) = acc_solver % mesh % nodes % y(indexes)
       vertex = elementcode / 100
       n = 0
       do nn = 1, vertex
          next = modulo(nn,vertex) + 1
          temp1 = acc_tempdiff(nn)
          temp2 = acc_tempdiff(next)
          if (((temp1 < 0.0) .and. (0.0 <= temp2)) .or. &
               ((temp2 <= 0.0) .and. (0.0 < temp1))) then
             n = n + 1
             if (n <= 2) then
                nelems = nelems + 1
                if (acc_isosurfallocated) then
                   acc_isosurf(nelems,1) = acc_nodes % x(nn) + &
                        temp1 * ((acc_nodes % x(next) - acc_nodes % x(nn)) / (temp1-temp2))
                   acc_isosurf(nelems,2) = acc_nodes % y(nn) + &
                        temp1 * ((acc_nodes % y(next) - acc_nodes % y(nn)) / (temp1-temp2))
                   acc_xmin = min(acc_isosurf(nelems,1), acc_xmin)
                   acc_xmax = max(acc_isosurf(nelems,1), acc_xmax)
                end if
             else
                call Warn('SteadyPhaseChange','Wiggly Isotherm')
                write(message,*) 'nodeindexes',indexes(1:vertex)
                call Warn('SteadyPhaseChange',message)
                write(message,*) 'temperature',acc_temperature(acc_tempperm(indexes))
                call Warn('SteadyPhaseChange',message)
                print *,'Nodes % x',acc_nodes % x(n)
                print *,'Nodes % y',acc_nodes % y(n)
             end if
          end if
       end do
       if (n == 1) then
          if (acc_isosurfallocated) then
             acc_isosurf(nelems,1) = 0.0_dp
             acc_isosurf(nelems,2) = 0.0_dp
          end if
          nelems = nelems - 1
          cycle
       end if
       if (acc_isosurfallocated .and. n == 2) then
          if (acc_isosurf(nelems-1,acc_tangentdir) > acc_isosurf(nelems,acc_tangentdir)) then
             temppi = acc_isosurf(nelems-1,1)
             acc_isosurf(nelems-1,1) = acc_isosurf(nelems,1)
             acc_isosurf(nelems,1) = temppi
             temppi = acc_isosurf(nelems-1,2)
             acc_isosurf(nelems-1,2) = acc_isosurf(nelems,2)
             acc_isosurf(nelems,2) = temppi
          end if
       end if
    end do
    acc_nelems = nelems
    if (nelems == 0) call Fatal('SteadyPhaseChange','Isotherm is empty thus cannot map phase change surface')
    if (.not. acc_isosurfallocated) then
       allocate(acc_isosurf(nelems+1,2))
       acc_isosurfallocated = .true.
       write(message,'(A,T35,I12)') 'Number of isotherm segments:',nelems
       call Info('SteadyPhaseChange',message)
       goto 100
    end if
    if (ListGetLogical(acc_params,'Save Isotherm',stat)) then
       call Info('SteadyPhaseChange','Isotherm saved to file isotherm.dat')
       open(10,file='isotherm.dat')
       do i = 1, nelems
          write(10,*) i,acc_isosurf(i,1),acc_isosurf(i,2)
       end do
       close(10)
    end if
  end subroutine sp_createisotherm
  ! -- main boundary loop ---------------------------------------------------
  subroutine sp_loopinit() bind(c, name="sp_loopinit")
    acc_area = 0.0
    acc_volume = 0.0
    acc_tave = 0.0
    acc_tabs = 0.0
    acc_volabs = 0.0
    acc_nodedone = .false.
    acc_maxtempdiff = 0.0
  end subroutine sp_loopinit

  function sp_nbound() result(n) bind(c, name="sp_nbound")
    integer(c_int) :: n
    n = acc_solver % mesh % numberofboundaryelements
  end function sp_nbound

  subroutine sp_elem_begin(t, n, family, doelem) bind(c, name="sp_elem_begin")
    integer(c_int), value :: t
    integer(c_int) :: n, family, doelem
    type(valuelist_t), pointer :: bc
    integer :: nn, elementcode
    logical :: gotit
    acc_elem => GetBoundaryElement(t)
    n = GetElementNOFNodes()
    family = GetElementFamily()
    if (family == 1) then
       doelem = 0
       return
    end if
    bc => GetBC()
    if (.not. GetLogical(bc, 'Phase Change', gotit)) then
       doelem = 0
       return
    end if
    doelem = 1
    acc_indexes => acc_elem % nodeindexes
    elementcode = acc_elem % type % elementcode
    if (elementcode < 200 .or. elementcode > 203) then
       call Fatal('PhaseChangeSolve','Implemented only for elements 202 and 203!')
    end if
    nn = n
    acc_nodes % x(1:nn) = acc_solver % mesh % nodes % x(acc_indexes)
    acc_nodes % y(1:nn) = acc_solver % mesh % nodes % y(acc_indexes)
    acc_nodes % z(1:nn) = acc_solver % mesh % nodes % z(acc_indexes)
    acc_tempdiff(1:nn) = acc_temperature(acc_tempperm(acc_indexes(1:nn))) - acc_meltpoint
    acc_maxtempdiff = max(acc_maxtempdiff, maxval(abs(acc_tempdiff(1:nn))))
  end subroutine sp_elem_begin

  subroutine sp_node_update(nn, te, newton) bind(c, name="sp_node_update")
    integer(c_int), value :: nn, te, newton
    integer :: k, i, imin
    real(kind=dp) :: xx, yy, x1, x2, y1, y2, epsi, dxmin, dymin, d, ttemp, dtdz, update
    logical :: stat
    k = acc_surfperm(acc_indexes(nn))
    if (acc_nodedone(k)) return
    acc_nodedone(k) = .true.
    if (acc_triplefixed .and. acc_indexes(nn) == acc_tripnode) then
       acc_surfacemove(k) = 0.0d0
       call Info('SteadyPhaseChange','Triple point position fixed')
       return
    end if
    if (nn == 3) then
       acc_surfacemove(k) = 0.5*sum(acc_surfacemove(acc_surfperm(acc_indexes(1:2))))
       update = 0.0
    end if
    ttemp = acc_temperature(acc_tempperm(te))
    if (newton /= 0) then
       dtdz = ttemp - acc_prevtemp(k)
       if (abs(dtdz) < aeps) then
          call Warn('SteadyPhaseChange','Very small temperature update.')
          dtdz = 1
       end if
       update = acc_surfacemove(k) * (acc_meltpoint - ttemp) / dtdz
    else
       epsi = 1.0d-6 * (acc_xmax - acc_xmin)
       if (acc_tangentdir == 1) then
          xx = acc_nodes % x(nn)
          yy = acc_nodes % y(nn)
       else
          xx = acc_nodes % y(nn)
          yy = acc_nodes % x(nn)
       end if
       dxmin = huge(dxmin)
       dymin = huge(dymin)
       stat = .false.
       do i = 1, acc_nelems-1, 2
          x1 = acc_isosurf(i,acc_tangentdir)
          x2 = acc_isosurf(i+1,acc_tangentdir)
          y1 = acc_isosurf(i,acc_normaldir)
          y2 = acc_isosurf(i+1,acc_normaldir)
          if ((xx > x1 - epsi) .and. (xx < x2 + epsi)) then
             dxmin = 0.0
             d = min(abs(yy - y1), abs(yy - y2))
             d = d + max(0.0d0, x1 - xx)
             d = d + max(0.0d0, xx - x2)
             if (d <= dymin) then
                stat = .true.
                dymin = d
                imin = i
             end if
          end if
          if (.not. stat) then
             d = min(abs(xx - x1), abs(xx - x2))
             if (d <= dxmin) then
                dxmin = d
                imin = i
             end if
          end if
       end do
       i = imin
       x1 = acc_isosurf(i,acc_tangentdir)
       x2 = acc_isosurf(i+1,acc_tangentdir)
       y1 = acc_isosurf(i,acc_normaldir)
       y2 = acc_isosurf(i+1,acc_normaldir)
       if (.not. stat) then
          if (dxmin > 1.0d-2*abs(x1 - x2)) then
             call Warn('SteadyPhaseChange','Isotherm error?')
             write(message,*) 'Nodeindexes',acc_indexes(nn)
             call Warn('SteadyPhaseChange',message)
             write(message,*) 'x:',xx,' y:',yy
             call Warn('SteadyPhaseChange',message)
             write(message,*) 'dxmin:',dxmin,' dymin:',dymin
             call Warn('SteadyPhaseChange',message)
          end if
       end if
       if (abs(x2 - x1) > aeps) then
          update = (y1 - yy) + (xx - x1) * (y2 - y1) / (x2 - x1)
       else
          update = 0.5_dp*(y1 + y2) - yy
       end if
    end if
    acc_prevtemp(k) = ttemp
    acc_surfacemove(k) = update
  end subroutine sp_node_update

  subroutine sp_elem_integrate() bind(c, name="sp_elem_integrate")
    type(GaussIntegrationPoints_t) :: ip
    integer :: i, nn
    real(kind=dp) :: u, v, w, s, detj
    logical :: ok
    nn = GetElementNOFNodes()
    ip = GaussPoints(acc_elem)
    do i = 1, ip % n
       u = ip % u(i)
       v = ip % v(i)
       w = ip % w(i)
       ok = ElementInfo(acc_elem, acc_nodes, u, v, w, detj, acc_basis)
       s = ip % s(i) * detj
       if (CurrentCoordinateSystem() /= Cartesian) then
          s = s * sum(acc_basis(1:nn) * acc_nodes % x(1:nn)) * 2.0 * PI
       end if
       acc_area = acc_area + s
       acc_volume = acc_volume + s*sum(acc_basis(1:nn) * &
            acc_surfacemove(acc_surfperm(acc_indexes(1:nn))))
       acc_tave = acc_tave + s*sum(acc_basis(1:nn) * acc_tempdiff(1:nn))
       acc_volabs = acc_volabs + s*sum(acc_basis(1:nn) * &
            abs(acc_surfacemove(acc_surfperm(acc_indexes(1:nn)))))
       acc_tabs = acc_tabs + s*sum(acc_basis(1:nn) * abs(acc_tempdiff(1:nn)))
    end do
  end subroutine sp_elem_integrate
  ! -- relaxation, surface update, norms, latches ----------------------------
  function sp_postloop(relax, visit, newton_in) result(nout) bind(c, name="sp_postloop")
    real(c_double), value :: relax
    integer(c_int), value :: visit, newton_in
    integer(c_int) :: nout
    integer :: i, j, k
    real(kind=dp) :: localrelax, cvol, clim, norm, relativechange
    real(kind=dp) :: maxsurfacemove, maxsurface, maxangle
    logical :: stat, newton
    newton = (newton_in /= 0)
    localrelax = relax
    acc_tave = acc_tave / acc_area
    acc_tabs = acc_tabs / acc_area
    acc_volume = acc_volume / acc_area
    acc_volabs = acc_volabs / acc_area
    i = ListGetInteger(acc_params,'Lumped Acceleration After Iterations', stat)
    if (stat .and. visit > i) then
       j = ListGetInteger(acc_params,'Lumped Acceleration Mode', stat)
       select case (j)
       case (1)
          cvol = 0.5*(acc_prevtave+acc_tave)/(acc_prevtave-acc_tave)
       case (2)
          cvol = 0.5*(acc_prevvolabs+acc_volabs)/(acc_prevvolabs-acc_volabs)
       case (3)
          cvol = 0.5*(acc_prevtabs+acc_tabs)/(acc_prevtabs-acc_tabs)
       case default
          cvol = 0.5*(acc_prevvolume+acc_volume)/(acc_prevvolume-acc_volume)
       end select
       if (cvol < 0.0) then
          cvol = 1.0
          acc_ccum = 1.0
       end if
       clim = ListGetConstReal(acc_params,'Lumped Acceleration Limiter', stat)
       if (.not. stat) clim = 100.0
       cvol = min(clim,cvol)
       cvol = max(1.0/clim,cvol)
       acc_ccum = acc_ccum * cvol
       write(message,'(A,T35,ES12.4)') 'Lumped Acceleration relaxation: ', acc_ccum
       call Info('SteadyPhaseChange',message)
       localrelax = localrelax * acc_ccum
    end if
    acc_surfacemove = localrelax * acc_surfacemove
    acc_surface = acc_surface + acc_surfacemove
    acc_dpos = acc_surfacemove(acc_surfperm(acc_tripnode))
    norm = 0.0d0
    relativechange = 0.0d0
    do k = 1, acc_model % mesh % numberofnodes
       if (.not. acc_boundarymarker(k)) cycle
       norm = norm + acc_surface(acc_surfperm(k)) ** 2.0
       relativechange = relativechange + acc_surfacemove(acc_surfperm(k)) ** 2.0
    end do
    relativechange = sqrt(relativechange / norm)
    norm = sqrt(norm / acc_nobnodes)
    maxsurfacemove = maxval(abs(acc_surfacemove))
    maxsurface = maxval(abs(acc_surface))
    maxangle = maxsurface / acc_width
    if (abs(relativechange) < acc_newtontol) newton = .true.
    write(message,'(A,T35,ES12.4)') 'Result Norm:', norm
    call Info('SteadyPhaseChange',message)
    write(message,'(A,T35,ES12.4)') 'Relative Change: ', relativechange
    call Info('SteadyPhaseChange',message)
    write(message,'(A,T35,ES12.4)') 'Max Height:', maxsurface
    call Info('SteadyPhaseChange',message)
    write(message,'(A,T35,ES12.4)') 'Max Change in Height:', maxsurfacemove
    call Info('SteadyPhaseChange',message)
    write(message,'(A,T35,ES12.4)') 'Max Angle:       ', maxangle
    call Info('SteadyPhaseChange',message)
    write(message,'(A,T35,ES12.4)') 'Maximum temperature difference: ', acc_maxtempdiff
    call Info('SteadyPhaseChange',message)
    acc_solver % variable % norm = norm
    acc_prevvolume = acc_volume
    acc_prevtave = acc_tave
    acc_prevtabs = acc_tabs
    acc_prevvolabs = acc_volabs
    if (acc_isosurfallocated) then
       deallocate(acc_isosurf)
       acc_isosurfallocated = .false.
    end if
    nout = merge(1, 0, newton)
  end function sp_postloop

  function sp_internal_mesh_movement() result(r) bind(c, name="sp_internal_mesh_movement")
    integer(c_int) :: r
    logical :: gotit
    if (ListGetLogical(acc_params,'Internal Mesh Movement',gotit)) then
       r = 1
    else
       r = 0
    end if
  end function sp_internal_mesh_movement

  ! -- BoxMoveMesh internal procedure (exact transcription) -----------------
  subroutine sp_boxmovemesh() bind(c, name="sp_boxmovemesh")
    type(element_t), pointer :: element
    integer, pointer :: indexes(:)
    real(kind=dp), pointer :: newy(:)
    real(kind=dp) :: x0, y0, ytop, ybot, coeff, yc, dy, q, r
    real(kind=dp) :: dx, dxmin, s1, s2
    real(kind=dp) :: xx, yy, x1, x2, y1, y2
    integer :: k, t, n, i, tmin
    logical :: hit
    coeff = 1.0_dp
    x0 = acc_solver % mesh % nodes % x(acc_tripnode)
    y0 = acc_solver % mesh % nodes % y(acc_tripnode)
    ytop = y0 + coeff * x0
    ybot = y0 - coeff * x0
    do k = 1, acc_solver % mesh % numberofnodes
       xx = acc_solver % mesh % nodes % x(k)
       if (xx > x0) cycle
       yy = acc_solver % mesh % nodes % y(k)
       if (yy > ytop) cycle
       if (yy < ybot) cycle
       if (acc_boundarymarker(k)) cycle
       hit = .false.
       do t = 1, acc_nophaseelements
          element => GetBoundaryElement(acc_phaseelements(t))
          n = GetElementNOFNodes()
          indexes => element % nodeindexes
          x1 = acc_solver % mesh % nodes % x(indexes(1))
          x2 = acc_solver % mesh % nodes % x(indexes(2))
          if ((xx - x1) * (x2 - xx) >= 0.0_dp) then
             hit = .true.
             exit
          end if
       end do
       if (.not. hit) then
          tmin = 0
          dxmin = huge(dxmin)
          do t = 1, acc_nophaseelements
             element => GetBoundaryElement(acc_phaseelements(t))
             n = GetElementNOFNodes()
             indexes => element % nodeindexes
             do i = 1, 2
                x1 = acc_solver % mesh % nodes % x(indexes(i))
                dx = abs(x1 - xx)
                if (dx < dxmin) then
                   dxmin = dx
                   tmin = t
                end if
             end do
          end do
          t = tmin
          element => GetBoundaryElement(acc_phaseelements(t))
          n = GetElementNOFNodes()
          indexes => element % nodeindexes
          x1 = acc_solver % mesh % nodes % x(indexes(1))
          x2 = acc_solver % mesh % nodes % x(indexes(2))
       end if
       y1 = acc_solver % mesh % nodes % y(indexes(1))
       y2 = acc_solver % mesh % nodes % y(indexes(2))
       q = (xx - x1) / (x2 - x1)
       yc = (1-q) * y1 + q * y2
       if (yy > yc) then
          r = 1.0_dp - (yy - yc) / (ytop - yc)
       else
          r = 1.0_dp - (yy - yc) / (ybot - yc)
       end if
       s1 = acc_surfacemove(acc_surfperm(indexes(1)))
       s2 = acc_surfacemove(acc_surfperm(indexes(2)))
       dy = (1-q) * s1 + q * s2
       acc_solver % mesh % nodes % y(k) = acc_solver % mesh % nodes % y(k) + r * dy
    end do
    do k = 1, acc_solver % mesh % numberofnodes
       if (.not. acc_boundarymarker(k)) cycle
       dy = acc_surfacemove(acc_surfperm(k))
       acc_solver % mesh % nodes % y(k) = acc_solver % mesh % nodes % y(k) + dy
    end do
  end subroutine sp_boxmovemesh

  ! -- label-200 tail: result listing + mesh-y variable ---------------------
  subroutine sp_finish(trip_temp) bind(c, name="sp_finish")
    real(c_double), value :: trip_temp
    real(kind=dp), pointer :: newy(:)
    call ListAddConstReal(acc_model % simulation,'res: Triple point temperature',trip_temp)
    call ListAddConstReal(acc_model % simulation,'res: triple point movement',acc_dpos)
    acc_helpsol => VariableGet(acc_solver % mesh % variables,'newy', ThisOnly=.TRUE.)
    if (.not. associated(acc_helpsol)) then
       newy => acc_solver % mesh % nodes % y
       call VariableAdd(acc_solver % mesh % variables, acc_solver % mesh, &
            acc_psolver,'NewY',1,newy)
    end if
  end subroutine sp_finish
  ! -- MeltingHeat UDF entry --------------------------------------------------
  ! -- Serial proof tree: no OpenMP, so the !$ lines are comments, tid = 1 --
  ! -- and a single Nodes slot is allocated, exactly as the original ------
  ! -- behaves when compiled without -fopenmp. -------------------------------
  subroutine sp_mh_init(h_model) bind(c, name="sp_mh_init")
    type(c_ptr), value :: h_model
    type(model_t), pointer :: model
    type(variable_t), pointer :: normalsol
    integer :: j, n, nthr
    if (acc_mh_initialized) return
    call c_f_pointer(h_model, model)
    normalsol => VariableGet(model % variables, 'Normals', ThisOnly=.TRUE.)
    acc_mh_normalexist = associated(normalsol)
    acc_mh_normalsol => normalsol
    nthr = 1
    n = model % mesh % maxelementnodes
    allocate(acc_mh_nodes(nthr))
    do j = 1, nthr
       allocate(acc_mh_nodes(j) % x(n), acc_mh_nodes(j) % y(n), acc_mh_nodes(j) % z(n))
    end do
    acc_mh_initialized = .true.
  end subroutine sp_mh_init

  function sp_mh_eval(h_model, node, t) result(flux) bind(c, name="sp_mh_eval")
    type(c_ptr), value :: h_model
    integer(c_int), value :: node
    real(c_double), value :: t
    real(c_double) :: flux
    type(model_t), pointer :: model
    type(element_t), pointer :: currentelement, parent
    integer, pointer :: nodeindexes(:)
    integer :: k, n, i, tid
    real(kind=dp) :: nodelatentheat, density, normalpull, u, v
    real(kind=dp) :: upull(3), normal(3), elemlatentheat(4)
    logical :: stat, found
    call c_f_pointer(h_model, model)
    tid = 1
    currentelement => GetCurrentElement()
    nodeindexes => currentelement % nodeindexes
    n = currentelement % type % numberofnodes
    k = ListGetInteger(model % bodies(currentelement % bodyid) % values,'Material')
    elemlatentheat(1:n) = ListGetReal(model % materials(k) % values, 'Latent Heat', n, nodeindexes)
    do i = 1, n
       if (nodeindexes(i) == node) exit
    end do
    if (nodeindexes(i) /= node) call Fatal('MeltingHeat','Node not found')
    nodelatentheat = elemlatentheat(i)
    upull = 0.0
    upull(1) = ListGetConstReal(model % simulation,'res: Pull Velocity 1',stat)
    if (.not. stat) upull(1) = ListGetConstReal(model % materials(k) % values,'Convection Velocity 1',stat)
    upull(2) = ListGetConstReal(model % simulation,'res: Pull Velocity 2',stat)
    if (.not. stat) upull(2) = ListGetConstReal(model % materials(k) % values,'Convection Velocity 2',stat)
    parent => currentelement % boundaryinfo % left
    k = ListGetInteger(model % bodies(parent % bodyid) % values,'Material')
    if (ListGetLogical(model % materials(k) % values, 'Solid', stat)) then
       density = ListGetConstReal(model % materials(k) % values, 'Density')
    else
       parent => currentelement % boundaryinfo % right
       k = ListGetInteger(model % bodies(parent % bodyid) % values,'Material')
       density = ListGetConstReal(model % materials(k) % values, 'Density')
    end if
    if (.not. associated(parent)) call Fatal('MeltingHeat','Parent not found')
    found = .false.
    if (acc_mh_normalexist) then
       normal = ConsistentNormalVector(currentmodel % solver, acc_mh_normalsol, currentelement, found, node=node)
    end if
    if (.not. found) then
       acc_mh_nodes(tid) % x(1:n) = model % nodes % x(nodeindexes)
       acc_mh_nodes(tid) % y(1:n) = model % nodes % y(nodeindexes)
       acc_mh_nodes(tid) % z(1:n) = model % nodes % z(nodeindexes)
       u = 0.0d0
       v = 0.0d0
       normal = NormalVector(currentelement, acc_mh_nodes(tid), u, v, .true.)
    end if
    normalpull = sum(normal * upull)
    flux = nodelatentheat * density * normalpull
  end function sp_mh_eval

end module steadyphasechange_acc
