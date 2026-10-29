! T022-A ElectricForce port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/ElectricForce.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays, or
! nothing). All Fortran-object traffic — TYPE(...) handles,
! OPTIONAL/POINTER dummies, assumed-shape descriptors, hidden lengths,
! list-directed Message writes, array-pointer results — stays inside
! Fortran, so the Rust side is pure C calls. The wrapper module mirrors
! StatElecForce's own locals (entry state, current element, parent,
! material, nodal arrays, Force/Moment/Area accumulators, entry-level
! flags) with identical lifetimes. All floating-point work (the Maxwell
! stress tensor integration, SUM/MATMUL reductions, SQRT) and all
! formatted/list-directed I/O stay in Fortran (same compiler), so
! bit-identity holds by construction wherever the call sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module electricforce_acc
  use iso_c_binding
  use types
  use defutils
  use integration
  use eleminfo
  use elementgeometry
  use coordinatesystems
  use messages
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(mesh_t), pointer, save :: acc_mesh => null()
  type(element_t), pointer, save :: acc_current => null()
  type(element_t), pointer, save :: acc_parent => null()
  type(valuelist_t), pointer, save :: acc_material => null()
  type(variable_t), pointer, save :: acc_potvar => null()
  type(variable_t), pointer, save :: acc_nodalvar => null()
  type(variable_t), pointer, save :: acc_fieldvar => null()
  type(variable_t), pointer, save :: acc_dvar => null()
  ! Entry locals (mirrors the solver body locals).
  type(nodes_t), target, save :: acc_elemnodes
  type(nodes_t), target, save :: acc_parentnodes
  real(c_double), allocatable, save :: acc_localpot(:)
  real(c_double), allocatable, save :: acc_perm(:,:,:)
  real(c_double), allocatable, save :: acc_nodalweight(:)
  real(c_double), pointer, save :: acc_potential(:) => null()
  real(c_double), pointer, save :: acc_forcedensity(:) => null()
  real(c_double), pointer, save :: acc_nodalforce(:) => null()
  real(c_double), pointer, save :: acc_pwrk(:,:,:) => null()
  integer, pointer, save :: acc_nodeindexes(:) => null()
  integer, pointer, save :: acc_potperm(:) => null()
  character(len=MAX_NAME_LEN), save :: acc_potname = ''
  character(len=MAX_NAME_LEN), save :: acc_varname = ''
  real(c_double), save :: acc_force(3) = 0.0_c_double
  real(c_double), save :: acc_moment(3) = 0.0_c_double
  real(c_double), save :: acc_momentabout(3) = 0.0_c_double
  real(c_double), save :: acc_area = 0.0_c_double
  real(c_double), save :: acc_vacuum = 1.0_c_double
  integer, save :: acc_dim = 0
  integer, save :: acc_n = 0
  integer, save :: acc_pn = 0
  logical, save :: acc_calcmoment = .false.
  logical, save :: acc_calcnodal = .false.
  logical, save :: acc_calcfield = .false.
  logical, save :: acc_dodisplaced = .false.

contains

  ! -- entry state ----------------------------------------------------
  subroutine efe_begin(h_model, h_solver) bind(c, name="efe_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    call Info('StatElecForce','Computing electric force on boundaries')
    acc_dim = CoordinateSystemDimension()
    acc_mesh => acc_solver % mesh
  end subroutine efe_begin

  ! -- potential field lookup (mirrors the original's stat logic exactly,
  ! -- including the custom-name path that skips the lookup) ------------
  subroutine efe_pot_setup() bind(c, name="efe_pot_setup")
    logical :: stat
    acc_potname = ListGetString( acc_solver % values, 'Potential Field Name', stat )
    if (.not. stat) acc_potname = 'Potential'
    if (.not. stat) then
      acc_potvar => VariableGet( acc_mesh % variables, acc_potname )
      if (.not. associated(acc_potvar)) then
        call Fatal('ElectricForce','Potential field does not exist:'//trim(acc_potname))
      end if
      acc_potential => acc_potvar % values
      acc_potperm => acc_potvar % perm
    end if
  end subroutine efe_pot_setup

  ! -- nodal-force variable ----------------------------------------------
  subroutine efe_nodal_setup() bind(c, name="efe_nodal_setup")
    logical :: stat
    acc_calcnodal = ListGetLogical( acc_solver % values, 'Calculate Nodal Force', stat )
    acc_nodalvar => VariableGet( acc_mesh % variables, 'Electric Nodal Force' )
    if (.not. associated(acc_nodalvar) .and. acc_calcnodal) then
      call VariableAddVector( acc_mesh % variables, acc_mesh, acc_solver, &
          'Electric Nodal Force', acc_dim, Perm = acc_potvar % perm )
      acc_nodalvar => VariableGet( acc_mesh % variables, 'Electric Nodal Force' )
    else if (associated(acc_nodalvar)) then
      acc_calcnodal = .true.
    end if
    if (acc_calcnodal) acc_nodalforce => acc_nodalvar % values
  end subroutine efe_nodal_setup

  ! -- force-density variable ---------------------------------------------
  subroutine efe_field_setup() bind(c, name="efe_field_setup")
    logical :: stat
    acc_calcfield = ListGetLogical( acc_solver % values, 'Calculate Force Density', stat )
    if (.not. acc_calcnodal) acc_calcfield = .true.
    acc_fieldvar => VariableGet( acc_mesh % variables, 'Electric Force Density' )
    if (.not. associated(acc_fieldvar) .and. acc_calcfield) then
      call VariableAddVector( acc_mesh % variables, acc_mesh, acc_solver, &
          'Electric Force Density', acc_dim, Perm = acc_potvar % perm )
      acc_fieldvar => VariableGet( acc_mesh % variables, 'Electric Force Density' )
    else if (associated(acc_fieldvar)) then
      acc_calcfield = .true.
    end if
    if (acc_calcfield) acc_forcedensity => acc_fieldvar % values
  end subroutine efe_field_setup

  ! -- per-invocation work arrays (mirrors the entry ALLOCATEs) ------------
  subroutine efe_alloc() bind(c, name="efe_alloc")
    integer :: n
    n = acc_mesh % maxelementnodes
    allocate( acc_elemnodes % x(n), acc_elemnodes % y(n), acc_elemnodes % z(n) )
    allocate( acc_parentnodes % x(n), acc_parentnodes % y(n), acc_parentnodes % z(n) )
    allocate( acc_localpot(n), acc_perm(3, 3, n) )
    allocate( acc_nodalweight(acc_mesh % numberofnodes) )
  end subroutine efe_alloc

  ! -- displaced-boundaries setup ------------------------------------------
  subroutine efe_displaced_setup() bind(c, name="efe_displaced_setup")
    logical :: stat
    acc_dodisplaced = ListGetLogical( acc_solver % values, 'All Displaced Boundaries', stat )
    if (acc_dodisplaced) then
      acc_varname = GetString( acc_solver % values, 'Displacement Field Name', stat )
      if (.not. stat) then
        acc_dvar => VariableGet( acc_mesh % variables, 'Displacement', .true. )
      else
        acc_dvar => VariableGet( acc_mesh % variables, acc_varname, .true. )
      end if
    end if
  end subroutine efe_displaced_setup

  ! -- constants + moment-about setup ---------------------------------------
  subroutine efe_const_setup() bind(c, name="efe_const_setup")
    logical :: stat
    acc_vacuum = ListGetConstReal( acc_model % constants, &
        'Permittivity Of Vacuum', stat )
    if (.not. stat) then
      call Warn('StatElecForce', 'Permittivity of Vacuum not given, using 1')
      acc_vacuum = 1.0d0
    end if
    acc_momentabout(1) = ListGetConstReal( acc_model % simulation, 'Moment About 1', stat )
    acc_calcmoment = stat
    acc_momentabout(2) = ListGetConstReal( acc_model % simulation, 'Moment About 2', stat )
    acc_calcmoment = stat .or. acc_calcmoment
    acc_momentabout(3) = ListGetConstReal( acc_model % simulation, 'Moment About 3', stat )
    acc_calcmoment = stat .or. acc_calcmoment
  end subroutine efe_const_setup

  ! -- zeroing (mirrors the whole-array assignments) -------------------------
  subroutine efe_zero() bind(c, name="efe_zero")
    acc_elemnodes % x = 0.0d0
    acc_elemnodes % y = 0.0d0
    acc_elemnodes % z = 0.0d0
    acc_parentnodes % x = 0.0d0
    acc_parentnodes % y = 0.0d0
    acc_parentnodes % z = 0.0d0
    acc_area = 0.0d0
    acc_force = 0.0d0
    acc_moment = 0.0d0
    acc_forcedensity = 0.0d0
    acc_nodalweight = 0.0_dp
  end subroutine efe_zero

  ! -- boundary sweep counts --------------------------------------------------
  subroutine efe_counts(nbulk, nbnd) bind(c, name="efe_counts")
    integer(c_int) :: nbulk, nbnd
    nbulk = acc_mesh % numberofbulkelements
    nbnd = acc_mesh % numberofboundaryelements
  end subroutine efe_counts

  ! -- per-element context ------------------------------------------------------
  subroutine efe_elem_begin(t, n, is101) bind(c, name="efe_elem_begin")
    integer(c_int), value :: t
    integer(c_int) :: n, is101
    acc_current => acc_mesh % elements(t)
    acc_model % currentelement => acc_mesh % elements(t)
    n = acc_current % type % numberofnodes
    acc_n = n
    acc_nodeindexes => acc_current % nodeindexes
    if (acc_current % type % elementcode == 101) then
      is101 = 1
    else
      is101 = 0
    end if
  end subroutine efe_elem_begin

  ! -- active-boundary test -------------------------------------------------------
  function efe_elem_active() result(r) bind(c, name="efe_elem_active")
    integer(c_int) :: r
    logical :: active, stat
    integer :: k
    active = .false.
    if (acc_dodisplaced) then
      active = all( acc_dvar % perm( acc_nodeindexes ) > 0 )
    else
      do k = 1, acc_model % numberofbcs
        if (acc_model % bcs(k) % tag /= acc_current % boundaryinfo % constraint) cycle
        if (ListGetLogical(acc_model % bcs(k) % values, &
            'Calculate Electric Force', stat)) active = .true.
      end do
    end if
    r = merge(1, 0, active)
  end function efe_elem_active

  ! -- per-element gather: nodes, parent, permittivity, potential -----------------
  subroutine efe_elem_compute() bind(c, name="efe_elem_compute")
    integer :: n, pn, i, j
    logical :: stat
    n = acc_n
    acc_elemnodes % x(1:n) = acc_mesh % nodes % x(acc_nodeindexes)
    acc_elemnodes % y(1:n) = acc_mesh % nodes % y(acc_nodeindexes)
    acc_elemnodes % z(1:n) = acc_mesh % nodes % z(acc_nodeindexes)
    acc_parent => acc_current % boundaryinfo % left
    stat = associated(acc_parent)
    if (stat) then
      i = acc_parent % bodyid
      j = ListGetInteger( acc_model % bodies(i) % values, 'Material', &
          minv=1, maxv=acc_model % numberofmaterials )
      acc_material => acc_model % materials(j) % values
      call ListGetRealArray( acc_material, 'Relative Permittivity', acc_pwrk, n, &
          acc_nodeindexes, stat )
      if (.not. stat) call ListGetRealArray( acc_material, &
          'Permittivity', acc_pwrk, n, acc_nodeindexes, stat )
    end if
    if (.not. stat) then
      acc_parent => acc_current % boundaryinfo % right
      if (.not. associated(acc_parent)) then
        write( message, * ) 'No permittivity found on specified boundary'
        call Fatal('StatElecForce', message)
      end if
      i = acc_parent % bodyid
      j = ListGetInteger( acc_model % bodies(i) % values, 'Material', &
          minv=1, maxv=acc_model % numberofmaterials )
      acc_material => acc_model % materials(j) % values
      call ListGetRealArray( acc_material, 'Relative Permittivity', acc_pwrk, n, &
          acc_nodeindexes, stat )
      if (.not. stat) call ListGetRealArray( acc_material, &
          'Permittivity', acc_pwrk, n, acc_nodeindexes, stat )
      if (.not. stat) then
        write( message, * ) 'No permittivity found on specified boundary'
        call Fatal('StatElecForce', message)
      end if
    end if
    stat = all( acc_potvar % perm( acc_nodeindexes ) > 0 )
    if (.not. stat) then
      write( message, * ) 'No potential available for specified boundary'
      call Fatal('StatElecForce', message)
    end if
    pn = acc_parent % type % numberofnodes
    acc_pn = pn
    acc_parentnodes % x(1:pn) = acc_mesh % nodes % x(acc_parent % nodeindexes)
    acc_parentnodes % y(1:pn) = acc_mesh % nodes % y(acc_parent % nodeindexes)
    acc_parentnodes % z(1:pn) = acc_mesh % nodes % z(acc_parent % nodeindexes)
    acc_perm = 0.0d0
    if (size(acc_pwrk,1) == 1) then
      do i = 1, 3
        acc_perm(i,i,1:n) = acc_pwrk(1,1,1:n)
      end do
    else if (size(acc_pwrk,2) == 1) then
      do i = 1, min(3,size(acc_pwrk,1))
        acc_perm(i,i,1:n) = acc_pwrk(i,1,1:n)
      end do
    else
      do i = 1, min(3,size(acc_pwrk,1))
        do j = 1, min(3,size(acc_pwrk,2))
          acc_perm(i,j,1:n) = acc_pwrk(i,j,1:n)
        end do
      end do
    end if
    acc_localpot = 0.0d0
    acc_localpot(1:pn) = acc_potential( acc_potperm( acc_parent % nodeindexes ) )
    call efe_integrate_work(n, pn)
  end subroutine efe_elem_compute

  ! -- MaxwellStressTensorIntegrate: exact transcription of the internal ------
  ! -- procedure (host association replaced by module SAVE state). ------------
  subroutine efe_integrate_work(n, pn)
    integer :: n, pn
    type(gaussintegrationpoints_t), target :: integstuff
    real(c_double), pointer :: u_integ(:), v_integ(:), w_integ(:), s_integ(:)
    real(c_double) :: basis(n), dbasisdx(n,3), ddbasisddx(n,3,3)
    real(c_double) :: parentbasis(pn), parentdbasisdx(pn,3)
    real(c_double) :: sqrtmetric, metric(3,3), symb(3,3,3), dsymb(3,3,3,3)
    real(c_double) :: u, v, w, s, detj, x(n), y(n), z(n)
    real(c_double) :: tensor(3,3), normal(3), efield(3), dflux(3), lforce(3), &
      lmoment(3), radius(3)
    real(c_double) :: elementweight(n), elementforce(3,n), xpos, ypos, zpos
    integer :: n_integ
    integer :: i, j, l
    logical :: stat

    elementweight = 0.0d0
    elementforce = 0.0d0
    integstuff = GaussPoints( acc_current )
    u_integ => integstuff % u
    v_integ => integstuff % v
    w_integ => integstuff % w
    s_integ => integstuff % s
    n_integ = integstuff % n
    do l = 1, n_integ
      u = u_integ(l)
      v = v_integ(l)
      w = w_integ(l)
      stat = ElementInfo( acc_current, acc_elemnodes, u, v, w, &
          detj, basis, dbasisdx, ddbasisddx, .false., .false. )
      s = 1.0d0
      if (CurrentCoordinateSystem() /= Cartesian) then
        xpos = sum( acc_elemnodes % x(1:n) * basis(1:n) )
        ypos = sum( acc_elemnodes % y(1:n) * basis(1:n) )
        zpos = sum( acc_elemnodes % z(1:n) * basis(1:n) )
        s = 2*PI
      end if
      call CoordinateSystemInfo( metric, sqrtmetric, symb, dsymb, xpos, ypos, zpos )
      s = s * sqrtmetric * detj * s_integ(l)
      normal = NormalVector( acc_current, acc_elemnodes, u, v, .true. )
      do i = 1, n
        do j = 1, pn
          if (acc_current % nodeindexes(i) == acc_parent % nodeindexes(j)) then
            x(i) = acc_parent % type % nodeu(j)
            y(i) = acc_parent % type % nodev(j)
            z(i) = acc_parent % type % nodew(j)
            exit
          end if
        end do
      end do
      u = sum( basis(1:n) * x(1:n) )
      v = sum( basis(1:n) * y(1:n) )
      w = sum( basis(1:n) * z(1:n) )
      stat = ElementInfo( acc_parent, acc_parentnodes, u, v, w, detj, parentbasis, &
          parentdbasisdx, ddbasisddx, .false., .false. )
      efield = 0.0d0
      dflux = 0.0d0
      do i = 1, acc_dim
        efield(i) = sum( parentdbasisdx(1:pn,i) * acc_localpot(1:pn) )
        do j = 1, acc_dim
          dflux(i) = dflux(i) + sum( parentdbasisdx(1:pn,j) * acc_localpot(1:pn) ) &
              * sum( acc_perm(i,j,1:n) * basis(1:n) )
        end do
      end do
      dflux = acc_vacuum * dflux
      tensor = 0.0d0
      do i = 1, acc_dim
        do j = 1, acc_dim
          tensor(i,j) = - dflux(i) * efield(j)
        end do
      end do
      do i = 1, acc_dim
        tensor(i,i) = tensor(i,i) + sum( dflux * efield ) / 2.0d0
      end do
      lforce = matmul( tensor, normal )
      acc_force = acc_force + s * lforce
      acc_area = acc_area + s * sum( basis(1:n) )
      if (acc_calcmoment) then
        radius(1) = sum( acc_elemnodes % x(1:n) * basis(1:n) ) - acc_momentabout(1)
        radius(2) = sum( acc_elemnodes % y(1:n) * basis(1:n) ) - acc_momentabout(2)
        radius(3) = sum( acc_elemnodes % z(1:n) * basis(1:n) ) - acc_momentabout(3)
        lmoment(1) = radius(2) * lforce(3) - radius(3) * lforce(2)
        lmoment(2) = radius(3) * lforce(1) - radius(1) * lforce(3)
        lmoment(3) = radius(1) * lforce(2) - radius(2) * lforce(1)
        acc_moment = acc_moment + s * lmoment
      end if
      do i = 1, n
        elementforce(:,i) = elementforce(:,i) + s * lforce * basis(i)
        elementweight(i) = elementweight(i) + s * basis(i)
      end do
    end do
    do i = 1, n
      do j = 1, acc_dim
        if (acc_calcnodal) then
          acc_nodalforce( acc_dim*(acc_potperm(acc_current % nodeindexes(i))-1)+j ) = &
              acc_nodalforce( acc_dim*(acc_potperm(acc_current % nodeindexes(i))-1)+j ) &
              + elementforce(j,i)
        end if
        if (acc_calcfield) then
          acc_forcedensity( acc_dim*(acc_potperm(acc_current % nodeindexes(i))-1)+j ) = &
              acc_forcedensity( acc_dim*(acc_potperm(acc_current % nodeindexes(i))-1)+j ) &
              + elementforce(j,i)
        end if
      end do
      acc_nodalweight(acc_potperm(acc_current % nodeindexes(i))) = &
          acc_nodalweight(acc_potperm(acc_current % nodeindexes(i))) + elementweight(i)
    end do
  end subroutine efe_integrate_work

  ! -- nodal averaging (mirrors the entry post-loop) ----------------------------
  subroutine efe_normalize() bind(c, name="efe_normalize")
    integer :: i, j
    do i = 1, acc_mesh % numberofnodes
      if (acc_nodalweight(i) > 0) then
        do j = 1, acc_dim
          acc_forcedensity(acc_dim*(i-1)+j) = acc_forcedensity(acc_dim*i-acc_dim+j) &
              / acc_nodalweight(i)
        end do
      end if
    end do
  end subroutine efe_normalize

  ! -- reports + res: results (mirrors the entry tail I/O) -----------------------
  subroutine efe_report() bind(c, name="efe_report")
    real(c_double) :: sf
    call Info('StatElecForce', ' ', level=4)
    if (acc_dim == 2) then
      write( message, '("Net electric force  : ", 3ES15.6E2 )' ) acc_force(1:2)
    else
      write( message, '("Net electric force  : ", 3ES15.6E2 )' ) acc_force
    end if
    call Info('StatElecForce', message, level=4)
    sf = sqrt(sum( acc_force**2 ) )
    write( message, '("Resultant force  : ", ES15.6 )' ) sf
    call Info('StatElecForce', message, level=4)
    if (acc_calcmoment) then
      call ListAddConstReal( acc_model % simulation, &
          'res: Electric Moment 3', acc_moment(3) )
      call ListAddConstReal( acc_model % simulation, &
          'res: Electric Moment 2', acc_moment(2) )
      call ListAddConstReal( acc_model % simulation, &
          'res: Electric Moment 1', acc_moment(1) )
    end if
    if (acc_dim > 2) call ListAddConstReal( acc_model % simulation, &
        'res: Electric Force 3', acc_force(3) )
    call ListAddConstReal( acc_model % simulation, &
        'res: Electric Force 2', acc_force(2) )
    call ListAddConstReal( acc_model % simulation, &
        'res: Electric Force 1', acc_force(1) )
  end subroutine efe_report

  ! -- teardown (mirrors the entry deallocations + invalidations) ------------------
  subroutine efe_finish() bind(c, name="efe_finish")
    if (acc_calcfield) then
      acc_fieldvar % primarymesh => acc_mesh
      call InvalidateVariable( acc_model % meshes, acc_mesh, 'Electric Force Density' )
      acc_fieldvar % valid = .true.
    end if
    if (acc_calcnodal) then
      acc_nodalvar % primarymesh => acc_mesh
      call InvalidateVariable( acc_model % meshes, acc_mesh, 'Electric Force Density' )
      acc_nodalvar % valid = .true.
    end if
    call SetCurrentMesh( acc_model, acc_solver % mesh )
    deallocate( acc_elemnodes % x, acc_elemnodes % y, acc_elemnodes % z )
    deallocate( acc_parentnodes % x, acc_parentnodes % y, acc_parentnodes % z )
    deallocate( acc_localpot, acc_perm )
    deallocate( acc_nodalweight )
  end subroutine efe_finish

end module electricforce_acc
