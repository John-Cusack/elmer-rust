! T022-A Electrokinetics port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/Electrokinetics.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays,
! C_CHAR buffers). All Fortran-object traffic — TYPE(...) handles,
! OPTIONAL/POINTER dummies, hidden lengths, assumed-shape descriptors,
! variable search state — stays inside Fortran, so the Rust side is pure
! C calls. The wrapper module mirrors the original's SAVE state
! (Nodes/DIM/FirstTime/vacuumPerm with identical lifetimes; the
! per-call SkipThisTime assignment is replicated by the caller sequence)
! plus the per-call element/material/equation context.
!
! All floating-point work (NormalVector, SUM reductions, the mobility
! quotient, the tangential projection) and all formatted I/O stay in
! Fortran (same compiler), so bit-identity holds by construction
! wherever the call sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module electrokinetics_acc
  use iso_c_binding
  use types
  use defutils
  use elementgeometry
  use coordinatesystems
  use messages
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(element_t), pointer, save :: acc_belem => null()
  type(element_t), pointer, save :: acc_pelem => null()
  type(valuelist_t), pointer, save :: acc_bc => null()
  type(valuelist_t), pointer, save :: acc_parentmat => null()
  type(valuelist_t), pointer, save :: acc_jmat => null()
  type(nodes_t), target, save :: acc_nodes
  integer, target, save :: acc_nodetarget(1) = 0
  integer, pointer, save :: acc_nodearr(:) => null()
  integer(c_int), save :: acc_dim = 0
  integer(c_int), save :: acc_eq = 0
  integer(c_int), save :: acc_mat = 0
  real(c_double), save :: acc_vacuum = 0.0_c_double
  real(c_double), save :: acc_normal(3) = 0.0_c_double
  logical, save :: acc_first = .true.
  character(len=MAX_NAME_LEN), save :: acc_efmethod = ''

contains

  ! -- shared entry state ---------------------------------------------
  subroutine eka_begin(h_model, node) bind(c, name="eka_begin")
    type(c_ptr), value :: h_model
    integer(c_int), value :: node
    call c_f_pointer(h_model, acc_model)
    acc_nodetarget(1) = node
    acc_nodearr => acc_nodetarget
  end subroutine eka_begin

  function eka_init_phase() result(skip) bind(c, name="eka_init_phase")
    integer(c_int) :: skip
    logical :: gotit, skiptime
    skiptime = ListGetLogical(acc_model % Simulation, &
        'Initializaton Phase', gotit)
    skip = merge(1, 0, skiptime)
  end function eka_init_phase

  ! -- helmholtz_smoluchowski_comp FirstTime block -----------------------
  function eka_setup(direction) result(d) bind(c, name="eka_setup")
    integer(c_int), value :: direction
    integer(c_int) :: d
    integer :: n, istat
    logical :: gotit
    if (acc_first) then
      acc_dim = CoordinateSystemDimension()
      if (acc_dim == 1) then
        call Fatal('electrokinetics (helmholtz_smoluchowski_comp)', &
            'The electrokinetic BCs do not make sense in 1D')
      end if
      if (direction == 0 .and. acc_dim > 2) then
        call Fatal('electrokinetics (helmholtz_smoluchowski_comp)', &
            'The tangential velocity may only be defined in 2D')
      end if
      n = acc_model % MaxElementNodes
      allocate(acc_nodes % x(n), acc_nodes % y(n), acc_nodes % z(n), &
          stat=istat)
      if (istat /= 0) then
        call Fatal('electrokinetics (helmholtz_smoluchowski_comp)', &
            'Allocations failed')
      end if
      acc_vacuum = ListGetConstReal(acc_model % Constants, &
          'Permittivity of Vacuum', gotit)
      if (.not. gotit) then
        call Warn('electrokinetics (helmholtz_smoluchowski_comp)', &
            'No value for >Permittivity of Vacuum< found in section Constants')
        call Warn('electrokinetics (helmholtz_smoluchowski_comp)', &
            '         Using default SI value 8.8542E-12')
        acc_vacuum = 8.8542e-12_c_double
      end if
      acc_first = .false.
    end if
    d = acc_dim
  end function eka_setup

  ! -- boundary/parent element resolution --------------------------------
  subroutine eka_find_parent() bind(c, name="eka_find_parent")
    integer :: other_body_id
    logical :: gotit
    acc_belem => acc_model % CurrentElement
    if (.not. associated(acc_belem)) then
      call Fatal('electrokinetics (helmholtz_smoluchowski_comp)', &
          'No boundary element found')
    end if
    other_body_id = acc_belem % BoundaryInfo % outbody
    if (other_body_id < 1) then
      acc_pelem => acc_belem % BoundaryInfo % Right
      if (.not. associated(acc_pelem)) &
          acc_pelem => acc_belem % BoundaryInfo % Left
    else
      acc_pelem => acc_belem % BoundaryInfo % Right
      if (acc_pelem % BodyId == other_body_id) &
          acc_pelem => acc_belem % BoundaryInfo % Left
    end if
    if (.not. associated(acc_pelem)) then
      call Fatal('electrokinetics (helmholtz_smoluchowski_comp)', &
          'No parent element found for boundary element')
    end if
    acc_mat = ListGetInteger( &
        acc_model % Bodies(acc_pelem % BodyId) % Values, 'Material', gotit)
    acc_parentmat => acc_model % Materials(acc_mat) % Values
    if (.not. associated(acc_parentmat)) then
      call Fatal('electrokinetics (helmholtz_smoluchowski_comp)', &
          'No material values could be found')
    end if
    acc_eq = ListGetInteger( &
        acc_model % Bodies(acc_pelem % BodyId) % Values, 'Equation', &
        minv=1, maxv=acc_model % NumberOfEquations)
  end subroutine eka_find_parent

  ! -- BC + boundary frame + normal --------------------------------------
  function eka_bc_normal(normal) result(bc_id) &
      bind(c, name="eka_bc_normal")
    real(c_double), intent(out) :: normal(3)
    integer(c_int) :: bc_id
    integer :: nboundary, j
    real(c_double) :: u, v
    bc_id = GetBCId(acc_belem)
    acc_bc => GetBC(acc_belem)
    nboundary = acc_belem % TYPE % NumberOfNodes
    do j = 1, nboundary
      if (acc_nodetarget(1) == acc_belem % NodeIndexes(j)) exit
    end do
    u = acc_belem % TYPE % NodeU(j)
    v = acc_belem % TYPE % NodeV(j)
    acc_nodes % x(1:nboundary) = &
        acc_model % Nodes % x(acc_belem % NodeIndexes)
    acc_nodes % y(1:nboundary) = &
        acc_model % Nodes % y(acc_belem % NodeIndexes)
    acc_nodes % z(1:nboundary) = &
        acc_model % Nodes % z(acc_belem % NodeIndexes)
    acc_normal = NormalVector(acc_belem, acc_nodes, u, v, .true.)
    normal = acc_normal
  end function eka_bc_normal

  ! -- zeta / EO mobility --------------------------------------------------
  subroutine eka_zeta(zeta, calc) bind(c, name="eka_zeta")
    real(c_double), intent(out) :: zeta
    integer(c_int), intent(out) :: calc
    real(c_double) :: dummyarray(1)
    logical :: calcmobility
    dummyarray = ListGetReal(acc_bc, 'Zeta Potential', 1, &
        acc_nodearr, calcmobility)
    zeta = dummyarray(1)
    calc = merge(1, 0, calcmobility)
  end subroutine eka_zeta

  subroutine eka_use_parent() bind(c, name="eka_use_parent")
    acc_model % CurrentElement => acc_pelem
  end subroutine eka_use_parent

  subroutine eka_viscosity(visc, got) bind(c, name="eka_viscosity")
    real(c_double), intent(out) :: visc
    integer(c_int), intent(out) :: got
    real(c_double) :: dummyarray(1)
    logical :: gotit
    dummyarray = ListGetReal(acc_parentmat, 'Viscosity', 1, &
        acc_nodearr, gotit)
    got = merge(1, 0, gotit)
    if (.not. gotit) then
      write(message, '(a,i0)') &
          'No viscosity found in material section no. ', acc_mat
      call Fatal('electrokinetics (helmholtz_smoluchowski)', message)
    else
      visc = dummyarray(1)
    end if
  end subroutine eka_viscosity

  subroutine eka_relperm(dc, got) bind(c, name="eka_relperm")
    real(c_double), intent(out) :: dc
    integer(c_int), intent(out) :: got
    real(c_double) :: dummyarray(1)
    logical :: gotit
    dummyarray = ListGetReal(acc_parentmat, 'Relative Permittivity', 1, &
        acc_nodearr, gotit)
    got = merge(1, 0, gotit)
    if (.not. gotit) then
      write(message, '(a,i0)') &
          'No keyword >Relative Permittivity< found in material section no. ', &
          acc_mat
      call Info('electrokinetics (helmholtz_smoluchowski)', message, level=4)
      call Info('electrokinetics (helmholtz_smoluchowski)', &
          'setting to 1', level=4)
      dc = 1.0_c_double
    else
      dc = dummyarray(1)
    end if
  end subroutine eka_relperm

  function eka_mob_from_zeta(zeta, visc, dc) result(mob) &
      bind(c, name="eka_mob_from_zeta")
    real(c_double), value :: zeta, visc, dc
    real(c_double) :: mob
    mob = (zeta * dc * acc_vacuum) / visc
  end function eka_mob_from_zeta

  subroutine eka_eo_direct(eo, got) bind(c, name="eka_eo_direct")
    real(c_double), intent(out) :: eo
    integer(c_int), intent(out) :: got
    real(c_double) :: dummyarray(1)
    logical :: gotit
    dummyarray = ListGetReal(acc_bc, 'EO Mobility', 1, &
        acc_nodearr, gotit)
    got = merge(1, 0, gotit)
    if (.not. gotit) then
      write(message, '(a,i0)') &
          'Neither zeta potential nor EO mobility found in boundary condition no. ', &
          eka_bc_id()
      call Info('electrokinetics (helmholtz_smoluchowski)', message, level=4)
      call Info('electrokinetics (helmholtz_smoluchowski)', &
          'setting EO mobility to 0', level=4)
      eo = 0.0_c_double
    else
      eo = dummyarray(1)
    end if
  end subroutine eka_eo_direct

  function eka_bc_id() result(bc_id)
    integer :: bc_id
    bc_id = GetBCId(acc_belem)
  end function eka_bc_id

  ! -- Electric Field method string ---------------------------------------
  subroutine eka_efmethod(buf, buflen, outlen, got) &
      bind(c, name="eka_efmethod")
    integer(c_int), value :: buflen
    character(kind=c_char), intent(out) :: buf(*)
    integer(c_int), intent(out) :: outlen, got
    character(len=:), allocatable :: s
    logical :: found
    integer :: i, n
    s = ListGetString(acc_model % Equations(acc_eq) % Values, &
        'Electric Field', found)
    got = merge(1, 0, found)
    if (found) acc_efmethod = s
    n = min(len_trim(s), buflen)
    do i = 1, n
      buf(i) = s(i:i)
    end do
    outlen = n
  end subroutine eka_efmethod

  subroutine eka_const_comp(i, val, got) bind(c, name="eka_const_comp")
    integer(c_int), value :: i
    real(c_double), intent(out) :: val
    integer(c_int), intent(out) :: got
    character(len=MAX_NAME_LEN) :: valuename
    real(c_double) :: dummyarray(1)
    logical :: gotit
    write(valuename, '(a,i2)') 'Electric Field', i
    dummyarray = ListGetReal(acc_parentmat, valuename, 1, &
        acc_nodearr, gotit)
    val = dummyarray(1)
    got = merge(1, 0, gotit)
  end subroutine eka_const_comp

  subroutine eka_warn_no_const() bind(c, name="eka_warn_no_const")
    write(message, '(a,i0)') &
        'No component for >Electric Field {1,2,3}< found in Material ', &
        acc_mat, ' although defined as constant'
    call Warn('electrokinetics (helmholtz_smoluchowski)', message)
  end subroutine eka_warn_no_const

  subroutine eka_comp_efield(i, val, found) bind(c, name="eka_comp_efield")
    integer(c_int), value :: i
    real(c_double), intent(out) :: val
    integer(c_int), intent(out) :: found
    character(len=MAX_NAME_LEN) :: valuename
    type(variable_t), pointer :: efieldcomp
    write(valuename, '(a,i2)') 'Electric Field', i
    efieldcomp => VariableGet(acc_model % Variables, valuename)
    if (associated(efieldcomp)) then
      val = efieldcomp % Values(efieldcomp % Perm(acc_nodetarget(1)))
      found = 1
    else
      found = 0
    end if
  end subroutine eka_comp_efield

  subroutine eka_info_no_efield() bind(c, name="eka_info_no_efield")
    write(message, '(a,i0)') &
        'No external electric field defined for Equation no. ', acc_eq
    call Info('electrokinetics (helmholtz_smoluchowski)', message, level=4)
  end subroutine eka_info_no_efield

  subroutine eka_warn_unknown() bind(c, name="eka_warn_unknown")
    write(message, '(a,a,a,i0)') 'Unknown entry, ', acc_efmethod, &
        ',for keyword >Electric Field< for Equation no. ', acc_eq
    call Warn('electrokinetics (helmholtz_smoluchowski)', message)
  end subroutine eka_warn_unknown

  subroutine eka_restore() bind(c, name="eka_restore")
    acc_model % CurrentElement => acc_belem
  end subroutine eka_restore

  ! -- HS velocity assembly (order-sensitive reductions stay here) ---------
  function eka_finish(direction, efield, eomob) result(hs) &
      bind(c, name="eka_finish")
    integer(c_int), value :: direction
    real(c_double), intent(in) :: efield(3)
    real(c_double), value :: eomob
    real(c_double) :: hs
    real(c_double) :: electricfield(3), tangefield(3), hsvel(3)
    real(c_double) :: normef, tangent(3), tangent2(3)
    electricfield = efield
    normef = sum(electricfield(1:3) * acc_normal(1:3))
    tangefield = electricfield - normef * acc_normal
    hsvel = tangefield * eomob
    if (direction == 0) then
      call TangentDirections(acc_normal, tangent, tangent2)
      hs = sum(hsvel * tangent)
    else
      hs = hsvel(direction)
    end if
  end function eka_finish

  ! -- getJouleHeat ----------------------------------------------------------
  subroutine ekj_deprecated() bind(c, name="ekj_deprecated")
    call Warn('GetJouleHeat', 'Deprecated function!')
  end subroutine ekj_deprecated

  function ekj_dim() result(d) bind(c, name="ekj_dim")
    integer(c_int) :: d
    d = CoordinateSystemDimension()
  end function ekj_dim

  function ekj_elem(body, eq, mat) result(code) bind(c, name="ekj_elem")
    integer(c_int), intent(out) :: body, eq, mat
    integer(c_int) :: code
    integer :: body_id, material_id, eq_id
    logical :: gotit
    body_id = acc_model % CurrentElement % BodyId
    eq_id = ListGetInteger( &
        acc_model % Bodies(body_id) % Values, 'Equation')
    material_id = ListGetInteger( &
        acc_model % Bodies(body_id) % Values, 'Material')
    acc_jmat => acc_model % Materials(material_id) % Values
    body = body_id
    eq = eq_id
    mat = material_id
    acc_eq = eq_id
    acc_mat = material_id
    if (.not. associated(acc_jmat)) then
      write(message, '(a,i3)') 'No Material found for body-id ', body_id
      call Warn('electrokinetics (getJouleHeat)', message)
      code = 1
    else
      code = 0
    end if
  end function ekj_elem

  subroutine ekj_conductivity(cond, got) bind(c, name="ekj_conductivity")
    real(c_double), intent(out) :: cond
    integer(c_int), intent(out) :: got
    real(c_double) :: dummyarray(1)
    logical :: gotit
    dummyarray = ListGetReal(acc_jmat, &
        'Electric Conductivity', 1, acc_nodearr, gotit)
    got = merge(1, 0, gotit)
    if (.not. gotit) then
      write(message, '(a,i0)') &
          'No >Electric Conductivity< found in Material section ', acc_mat
      call Warn('electrokinetics (getJouleHeat)', message)
    else
      cond = dummyarray(1)
    end if
  end subroutine ekj_conductivity

  subroutine ekj_density(dens, got) bind(c, name="ekj_density")
    real(c_double), intent(out) :: dens
    integer(c_int), intent(out) :: got
    real(c_double) :: dummyarray(1)
    logical :: gotit
    dummyarray = ListGetReal(acc_jmat, 'Density', 1, acc_nodearr, gotit)
    got = merge(1, 0, gotit)
    if (.not. gotit) then
      write(message, '(a,i0)') &
          'No >Density< found in Material section ', acc_mat
      call Warn('electrokinetics (getJouleHeat)', message)
      call Warn('electrokinetics (getJouleHeat)', &
          'setting reference density to 1')
      dens = 1.0_c_double
    else
      dens = dummyarray(1)
    end if
  end subroutine ekj_density

  subroutine ekj_const_comp(i, val, got) bind(c, name="ekj_const_comp")
    integer(c_int), value :: i
    real(c_double), intent(out) :: val
    integer(c_int), intent(out) :: got
    character(len=MAX_NAME_LEN) :: valuename
    real(c_double) :: dummyarray(1)
    logical :: gotit
    write(valuename, '(a,i2)') 'Electric Field', i
    dummyarray = ListGetReal(acc_jmat, valuename, 1, &
        acc_nodearr, gotit)
    val = dummyarray(1)
    got = merge(1, 0, gotit)
  end subroutine ekj_const_comp

  subroutine ekj_warn_no_const() bind(c, name="ekj_warn_no_const")
    write(message, '(a,i3)') &
        'No component for >Electric Field {1,2,3}< found in Material', &
        acc_mat, ' although defined as constant'
    call Warn('electrokinetics (getJouleHeat)', message)
  end subroutine ekj_warn_no_const

  subroutine ekj_info_no_efield() bind(c, name="ekj_info_no_efield")
    write(message, '(a,i3)') &
        'No external electric field defined for Equation no.', acc_eq
    call Info('electrokinetics (helmholtz_smoluchowski)', message, level=4)
  end subroutine ekj_info_no_efield

  subroutine ekj_warn_unknown() bind(c, name="ekj_warn_unknown")
    write(message, '(a,a,a,i3)') 'Unknown entry, ', acc_efmethod, &
        ',for keyword >Electric Field< for Equation no.', acc_eq
    call Warn('electrokinetics (getJouleHeat)', message)
  end subroutine ekj_warn_unknown

  function ekj_finish(efield, cond, dens) result(joule) &
      bind(c, name="ekj_finish")
    real(c_double), intent(in) :: efield(3)
    real(c_double), value :: cond, dens
    real(c_double) :: joule
    real(c_double) :: efsq
    efsq = sum(efield(1:3) * efield(1:3))
    joule = cond * efsq / dens
  end function ekj_finish

end module electrokinetics_acc
