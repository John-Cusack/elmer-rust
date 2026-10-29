! Clean-room BIND(C) wrappers for the Elmer calls in
! fem/src/modules/VectorHelmholtzLumping.F90, compiled against the built
! .mod files of the target configure. Statement-for-statement wrappers
! around Elmer library calls, no upstream text.
!
! Every wrapper takes only C-interoperable arguments. All Fortran-object
! traffic (derived-type handles, OPTIONAL/POINTER dummies, assumed-shape
! descriptors, hidden lengths) stays inside Fortran, so the Rust side is
! pure C calls.
!
! Lifetime mirroring: the original's internal-procedure SAVE state
! (Nodes, ParentNodes, the twelve ValueHandles, WBasis/RotWBasis/Basis/
! dBasisdx/e_local/phi_local plus AllocationsDone) lives here as
! module SAVE state with identical lifetimes. The main-body SAVE arrays
! (IntPoynt/IntEl/IntCurr/IntPot/IntVolt/IntCenter/IntWeight/IntNorm/
! CenterNode/jMode/NoModes/Visited) live on the Rust side; per-IP
! accumulation addends cross as plain doubles so the Rust loop does
! only program-order adds while every SUM/MATMUL reduction stays in
! Fortran, bit-identical by construction.
module vhl_acc
  use iso_c_binding
  use types
  use defutils
  use integration
  use eleminfo
  use messages
  implicit none

  character(*), parameter :: caller = 'VectorHelmholtzLumping'

  ! Entry handles (re-pointed on every invocation).
  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(solver_t), pointer, save :: acc_psolver => null()
  type(mesh_t), pointer, save :: acc_mesh => null()
  type(valuelist_t), pointer, save :: acc_params => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(valuelist_t), pointer, save :: acc_bc => null()
  type(variable_t), pointer, save :: acc_evar => null()
  type(variable_t), pointer, save :: acc_potvar => null()

  ! Edge-element style cached per entry (host association in original).
  logical, save :: acc_piola = .false.
  integer, save :: acc_degree = 0

  ! LocalIntegBC SAVE state (internal-procedure SAVE in the original).
  type(nodes_t), target, save :: acc_nodes, acc_parent_nodes
  type(element_t), pointer, save :: acc_parent => null()
  type(valuehandle_t), save :: acc_h_magload, acc_h_elrobin, acc_h_mu, &
      acc_h_absorb, acc_h_temre, acc_h_temim, acc_h_transfer, acc_h_elcurrent, &
      acc_h_relnu, acc_h_cond, acc_h_currdens, acc_h_eps
  logical, save :: acc_lib_alloc = .false.
  integer, save :: acc_lib_m = 0, acc_lib_ndofs = 0, acc_nd = 0
  real(kind=dp), allocatable, save :: acc_wbasis(:,:), acc_rotwbasis(:,:), &
      acc_basis(:), acc_dbasisdx(:,:), acc_elocal(:,:), acc_philocal(:,:)
  real(kind=dp), save :: acc_normal(3) = 0.0_dp
  type(gaussintegrationpoints_t), save :: acc_ip

contains

  ! -- entry ----------------------------------------------------------
  subroutine vhl_begin(h_model, h_solver) bind(c, name="vhl_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine vhl_begin

  subroutine vhl_info_header() bind(c, name="vhl_info_header")
    call info(caller, '', level=6)
    call info(caller, '----------------------------------------------------------', level=6)
    call info(caller, 'Computing derived fields for electromagnetic wave equation!', level=4)
  end subroutine vhl_info_header

  subroutine vhl_primary(soln) bind(c, name="vhl_primary")
    integer(c_int) :: soln
    logical :: found
    acc_params => getsolverparams()
    soln = listgetinteger(acc_params, 'Primary Solver Index', found)
    if (soln == 0) then
      call fatal(caller, 'We should know > Primary Solver Index <')
    end if
    acc_psolver => acc_model % solvers(soln)
  end subroutine vhl_primary

  function vhl_omega() result(w) bind(c, name="vhl_omega")
    real(c_double) :: w
    logical :: found
    w = getangularfrequency(acc_psolver % values, found=found)
    if (.not. found) call fatal(caller, 'We need angular frequency!')
  end function vhl_omega

  subroutine vhl_constants(m0, e0) bind(c, name="vhl_constants")
    real(c_double) :: m0, e0
    logical :: found
    found = .false.
    if (associated(acc_model % constants)) then
      m0 = 1.0_dp / getconstreal(acc_model % constants, 'Permeability of Vacuum', found)
    end if
    if (.not. found) m0 = 1.0_dp / (pi * 4.0d-7)
    found = .false.
    if (associated(acc_model % constants)) then
      e0 = getconstreal(acc_model % constants, 'Permittivity of Vacuum', found)
    end if
    if (.not. found) e0 = 8.854187817d-12
  end subroutine vhl_constants

  function vhl_count_modes() result(nm) bind(c, name="vhl_count_modes")
    integer(c_int) :: nm
    integer :: i, j
    logical :: found
    nm = 0
    do i = 1, acc_model % numberofbcs
      j = listgetinteger(acc_model % bcs(i) % values, 'Constraint Mode', found)
      nm = max(nm, j)
    end do
  end function vhl_count_modes
  function vhl_nodal_flag() result(r) bind(c, name="vhl_nodal_flag")
    integer(c_int) :: r
    logical :: v, found
    v = listgetlogical(acc_params, 'Nodal Target Field', found)
    r = merge(1, 0, v)
  end function vhl_nodal_flag

  subroutine vhl_fatal_nofield() bind(c, name="vhl_fatal_nofield")
    call fatal(caller, 'Could not find nodal electric field!')
  end subroutine vhl_fatal_nofield

  function vhl_find_potential() result(r) bind(c, name="vhl_find_potential")
    integer(c_int) :: r
    acc_potvar => variableget(acc_mesh % variables, 'Potential', thisonly=.true.)
    r = merge(1, 0, associated(acc_potvar))
  end function vhl_find_potential

  subroutine vhl_find_enodal(found, dofs) bind(c, name="vhl_find_enodal")
    integer(c_int) :: found, dofs
    acc_evar => variableget(acc_mesh % variables, 'Electric field e', thisonly=.true.)
    if (.not. associated(acc_evar)) then
      acc_evar => variableget(acc_mesh % variables, 'Electric field', thisonly=.true.)
    end if
    found = merge(1, 0, associated(acc_evar))
    if (associated(acc_evar)) then
      dofs = acc_evar % dofs
    else
      dofs = 0
    end if
  end subroutine vhl_find_enodal

  subroutine vhl_fatal_dofs(nd) bind(c, name="vhl_fatal_dofs")
    integer(c_int), value :: nd
    call fatal(caller, 'Nodal mode assumes exactly 6 dofs, not '//i2s(nd))
  end subroutine vhl_fatal_dofs

  subroutine vhl_info_nodal() bind(c, name="vhl_info_nodal")
    call info(caller, 'Assuming electric field living on nodes')
  end subroutine vhl_info_nodal

  subroutine vhl_info_edge() bind(c, name="vhl_info_edge")
    call info(caller, 'Assuming electric field to live in Hcurl')
  end subroutine vhl_info_edge

  subroutine vhl_find_enamed(found) bind(c, name="vhl_find_enamed")
    integer(c_int) :: found
    character(len=MAX_NAME_LEN) :: pname
    logical :: f
    pname = listgetstring(acc_params, 'Target Variable', f)
    found = merge(1, 0, f)
    if (f) then
      acc_evar => variableget(acc_psolver % mesh % variables, pname)
    else
      acc_evar => acc_psolver % variable
      pname = getvarname(acc_psolver % variable)
    end if
    call info(caller, 'Name of target variable: '//trim(pname), level=10)
  end subroutine vhl_find_enamed

  subroutine vhl_check_edofs() bind(c, name="vhl_check_edofs")
    if (acc_evar % dofs /= 2) call fatal(caller, 'Primary field should have two components!')
  end subroutine vhl_check_edofs

  subroutine vhl_edge_style(piola, degree) bind(c, name="vhl_edge_style")
    integer(c_int) :: piola, degree
    logical :: piv
    integer :: deg
    call edgeelementstyle(acc_psolver % values, piv, basisdegree=deg)
    acc_piola = piv
    acc_degree = deg
    piola = merge(1, 0, piv)
    degree = deg
    if (piv) call info(caller, 'Using Piola transformed finite elements', level=8)
  end subroutine vhl_edge_style

  subroutine vhl_getmesh_touch() bind(c, name="vhl_getmesh_touch")
    acc_mesh => getmesh()
  end subroutine vhl_getmesh_touch

  ! -- boundary loop ---------------------------------------------------
  function vhl_nof_bnd() result(n) bind(c, name="vhl_nof_bnd")
    integer(c_int) :: n
    n = getnofboundaryelements()
  end function vhl_nof_bnd

  function vhl_begin_boundary(t) result(family) bind(c, name="vhl_begin_boundary")
    integer(c_int), value :: t
    integer(c_int) :: family
    acc_elem => getboundaryelement(t)
    acc_bc => getbc()
    family = getelementfamily()
  end function vhl_begin_boundary

  subroutine vhl_remap_edge() bind(c, name="vhl_remap_edge")
    integer :: k
    k = getboundaryedgeindex(acc_elem, 1)
    acc_elem => acc_mesh % edges(k)
  end subroutine vhl_remap_edge

  subroutine vhl_remap_face() bind(c, name="vhl_remap_face")
    integer :: k
    k = getboundaryfaceindex(acc_elem)
    acc_elem => acc_mesh % faces(k)
  end subroutine vhl_remap_face

  function vhl_bc_mode() result(m) bind(c, name="vhl_bc_mode")
    integer(c_int) :: m
    logical :: found
    m = listgetinteger(acc_bc, 'Constraint Mode', found)
  end function vhl_bc_mode

  ! -- parallel guards (serial proof never takes the PEs>1 branch) -----
  function vhl_pes() result(n) bind(c, name="vhl_pes")
    integer(c_int) :: n
    n = parenv % pes
  end function vhl_pes

  function vhl_mype() result(n) bind(c, name="vhl_mype")
    integer(c_int) :: n
    n = parenv % mype
  end function vhl_mype

  subroutine vhl_par_z(z) bind(c, name="vhl_par_z")
    complex(c_double_complex) :: z
    complex(kind=dp) :: w
    w = z
    w = parallelreduction(w)
    z = w
  end subroutine vhl_par_z

  subroutine vhl_par_r(x) bind(c, name="vhl_par_r")
    real(c_double) :: x
    real(kind=dp) :: w
    w = x
    w = parallelreduction(w)
    x = w
  end subroutine vhl_par_r

  ! -- LocalIntegBC setup ----------------------------------------------
  subroutine vhl_lib_setup(nodal, has_pot, first, nip, found_e) bind(c, name="vhl_lib_setup")
    integer(c_int), value :: nodal, has_pot, first
    integer(c_int) :: nip, found_e
    logical :: found
    integer :: m, ndofs
    if (first /= 0) then
      call listinitelementkeyword(acc_h_elrobin, 'Boundary Condition', &
          'Electric Robin Coefficient', initim=.true.)
      call listinitelementkeyword(acc_h_magload, 'Boundary Condition', &
          'Magnetic Boundary Load', initim=.true., initvec3d=.true.)
      call listinitelementkeyword(acc_h_absorb, 'Boundary Condition', 'Absorbing BC')
      call listinitelementkeyword(acc_h_temre, 'Boundary Condition', 'TEM Potential')
      call listinitelementkeyword(acc_h_temim, 'Boundary Condition', 'TEM Potential Im')
      call listinitelementkeyword(acc_h_mu, 'Material', 'Relative Reluctivity', initim=.true.)
      call listinitelementkeyword(acc_h_eps, 'Material', 'Relative Permittivity', initim=.true.)
      call listinitelementkeyword(acc_h_transfer, 'Boundary Condition', &
          'Electric Transfer Coefficient', initim=.true.)
      call listinitelementkeyword(acc_h_elcurrent, 'Boundary Condition', &
          'Electric Current Density', initim=.true.)
      call listinitelementkeyword(acc_h_currdens, 'Body Force', 'Current Density', &
          initim=.true., initvec3d=.true.)
      call listinitelementkeyword(acc_h_cond, 'Material', 'Electric Conductivity')
    end if
    if (.not. acc_lib_alloc) then
      m = acc_mesh % maxelementdofs
      ndofs = acc_evar % dofs
      allocate(acc_wbasis(m,3), acc_rotwbasis(m,3), acc_basis(m), &
          acc_dbasisdx(m,3), acc_elocal(ndofs,m), acc_philocal(2,m))
      acc_lib_m = m
      acc_lib_ndofs = ndofs
      acc_lib_alloc = .true.
    end if
    call getelementnodes(acc_nodes, acc_elem)
    acc_parent => getbulkelementatboundary(acc_elem)
    if (.not. associated(acc_parent)) &
        call fatal(caller, 'Model lumping requires parent element!')
    acc_elocal = 0.0_dp
    if (has_pot /= 0) then
      acc_philocal = 0.0_dp
      call getvectorlocalsolution(acc_philocal, uvariable=acc_potvar, uelement=acc_elem)
    end if
    if (nodal /= 0) then
      call getvectorlocalsolution(acc_elocal, uelement=acc_elem, &
          uvariable=acc_evar, found=found)
    else
      call getelementnodes(acc_parent_nodes, acc_parent)
      call getvectorlocalsolution(acc_elocal, uelement=acc_parent, &
          uvariable=acc_evar, usolver=acc_psolver, found=found)
      ! Dead in the original (np is overwritten by 0 before any use):
      !   np = n * MAXVAL(Solver % Def_Dofs(GetElementFamily(Parent),:,1))
      acc_nd = getelementnofdofs(acc_parent, usolver=acc_psolver)
    end if
    found_e = merge(1, 0, found)
    if (found_e == 0) call fatal(caller, 'Could not find field data on boundary!?')
    acc_normal = normalvector(acc_elem, acc_nodes, check=.true.)
    if (nodal /= 0) then
      acc_ip = gausspoints(acc_elem)
    else
      acc_ip = gausspoints(acc_elem, preferenceelement=acc_piola, &
          edgebasisdegree=acc_degree)
    end if
    nip = acc_ip % n
  end subroutine vhl_lib_setup

  ! -- LocalIntegBC integration point: full IP body, returns the
  ! -- thirteen accumulation addends (weight already folded in).
  subroutine vhl_lib_ip(t, nodal, has_pot, omega, m0, e0, add) bind(c, name="vhl_lib_ip")
    integer(c_int), value :: t, nodal, has_pot
    real(c_double), value :: omega, m0, e0
    real(c_double) :: add(13)
    complex(kind=dp) :: b, zs, l(3), muinv, temgrad(3), betapar, jn, eps, &
        e_ip(3), e_ip_norm, e_ip_tan(3), f_ip_tan(3), imu, phi, s
    real(kind=dp) :: weight, detj, normal(3), cond, u, v, w, x, y, z
    real(kind=dp) :: gr(3), gi(3)
    logical :: stat, found
    integer :: i, n
    n = acc_elem % type % numberofnodes
    imu = cmplx(0.0_dp, 1.0_dp, kind=dp)
    if (nodal /= 0) then
      stat = elementinfo(acc_elem, acc_nodes, acc_ip % u(t), acc_ip % v(t), &
          acc_ip % w(t), detj, acc_basis, acc_dbasisdx)
    else
      stat = elementinfo(acc_elem, acc_nodes, acc_ip % u(t), acc_ip % v(t), &
          acc_ip % w(t), detj, acc_basis, acc_dbasisdx, &
          edgebasis=acc_wbasis, usolver=acc_psolver)
    end if
    weight = acc_ip % s(t) * detj
    b = listgetelementcomplex(acc_h_elrobin, acc_basis, acc_elem, found, gausspoint=t)
    muinv = listgetelementcomplex(acc_h_mu, acc_basis, acc_parent, found, gausspoint=t)
    if (found) then
      muinv = muinv * m0
    else
      muinv = m0
    end if
    eps = listgetelementcomplex(acc_h_eps, acc_basis, acc_parent, found, gausspoint=t)
    if (found) then
      eps = eps * e0
    else
      eps = e0
    end if
    cond = listgetelementreal(acc_h_cond, acc_basis, acc_parent, found, gausspoint=t)
    zs = imu * omega / (b * muinv)
    x = sum(acc_basis(1:n) * acc_nodes % x(1:n))
    y = sum(acc_basis(1:n) * acc_nodes % y(1:n))
    z = sum(acc_basis(1:n) * acc_nodes % z(1:n))
    add(1) = weight
    add(2) = weight * x
    add(3) = weight * y
    add(4) = weight * z
    l = listgetelementcomplex3d(acc_h_magload, acc_basis, acc_elem, found, gausspoint=t)
    gr = listgetelementrealgrad(acc_h_temre, acc_dbasisdx, acc_elem, found)
    gi = listgetelementrealgrad(acc_h_temim, acc_dbasisdx, acc_elem, found)
    temgrad = cmplx(gr, gi, kind=dp)
    l = l + temgrad
    b = listgetelementcomplex(acc_h_elrobin, acc_basis, acc_elem, found, gausspoint=t)
    l = l / (2 * sqrt(2.0_dp) * b)
    normal = acc_normal
    if (nodal /= 0) then
      do i = 1, 3
        e_ip(i) = cmplx(sum(acc_basis(1:n) * acc_elocal(i,1:n)), &
            sum(acc_basis(1:n) * acc_elocal(i+3,1:n)), kind=dp)
      end do
    else
      u = 0.0_dp
      v = 0.0_dp
      w = 0.0_dp
      call findparentuvw(acc_elem, n, acc_parent, &
          acc_parent % type % numberofnodes, u, v, w, acc_basis)
      if (getelementfamily(acc_elem) == 2) then
        stat = edgeelementinfo(acc_parent, acc_parent_nodes, u, v, w, &
            detf=detj, basis=acc_basis, edgebasis=acc_wbasis, &
            rotbasis=acc_rotwbasis, dbasisdx=acc_dbasisdx, &
            basisdegree=acc_degree, applypiolatransform=.true.)
      else
        stat = elementinfo(acc_parent, acc_parent_nodes, u, v, w, detj, &
            acc_basis, acc_dbasisdx, edgebasis=acc_wbasis, &
            rotbasis=acc_rotwbasis, usolver=acc_psolver)
      end if
      e_ip(1:3) = cmplx(matmul(acc_elocal(1,1:acc_nd), acc_wbasis(1:acc_nd,1:3)), &
          matmul(acc_elocal(2,1:acc_nd), acc_wbasis(1:acc_nd,1:3)), kind=dp)
    end if
    e_ip_norm = sum(e_ip * normal)
    e_ip_tan = e_ip - e_ip_norm * normal
    s = sum(l * conjg(l))
    add(5) = weight * abs(s)
    s = sum(e_ip_tan * conjg(l))
    s = weight * s
    add(6) = real(s)
    add(7) = aimag(s)
    s = sum(e_ip_tan * conjg(e_ip_tan))
    s = weight * 0.5_dp * s / zs
    add(8) = real(s)
    add(9) = aimag(s)
    s = weight * (imu * omega * eps + cond) * e_ip_norm
    add(10) = real(s)
    add(11) = aimag(s)
    if (has_pot /= 0) then
      call getvectorlocalsolution(acc_philocal, uelement=acc_elem, uvariable=acc_potvar)
      phi = cmplx(sum(acc_basis(1:n) * acc_philocal(1,1:n)), &
          sum(acc_basis(1:n) * acc_philocal(2,1:n)), kind=dp)
      s = weight * phi
      add(12) = real(s)
      add(13) = aimag(s)
    else
      add(12) = 0.0_dp
      add(13) = 0.0_dp
    end if
  end subroutine vhl_lib_ip

  ! -- main-body PRINTs (list-directed, kept in Fortran) ----------------
  subroutine vhl_print_el_before(jm, nm, el, cnorm) bind(c, name="vhl_print_el_before")
    integer(c_int), value :: jm, nm
    complex(c_double_complex) :: el(nm, nm)
    real(c_double) :: cnorm(nm)
    print *, 'Elfield before norm through port:', el(jm+1, :)
    print *, 'Elfield norm:', cnorm
  end subroutine vhl_print_el_before

  subroutine vhl_print_el_after(jm, nm, el, wt, poynt, curr) bind(c, name="vhl_print_el_after")
    integer(c_int), value :: jm, nm
    complex(c_double_complex) :: el(nm, nm), poynt(nm, nm), curr(nm, nm)
    real(c_double) :: wt(nm)
    print *, 'Elfield itse through port:', el(jm+1, :)
    print *, 'Elfield area port:', wt(jm+1)
    print *, 'Elfield Energy through port:', poynt(jm+1, :)
    print *, 'Elfield Current through port:', curr(jm+1, :)
  end subroutine vhl_print_el_after

  subroutine vhl_print_pot(jm, nm, pot) bind(c, name="vhl_print_pot")
    integer(c_int), value :: jm, nm
    complex(c_double_complex) :: pot(nm, nm)
    print *, 'Average pot on port:', pot(jm+1, :)
  end subroutine vhl_print_pot

  ! -- CenterPortLoc (exact transcription incl. PRINTs and Fatal) -------
  subroutine vhl_center(nm, center, wt, cnode) bind(c, name="vhl_center")
    integer(c_int), value :: nm
    real(c_double) :: center(nm, 3), wt(nm)
    integer(c_int) :: cnode(nm)
    integer :: i, imode, mini
    real(kind=dp) :: mindist2, dist2, coord0(3), coord1(3)
    print *, 'IntWeight:', wt
    if (any(wt < epsilon(dist2))) then
      print *, 'IntWeight:', wt
      call fatal(caller, 'Some weight is zero!')
    end if
    do i = 1, 3
      center(:, i) = center(:, i) / wt(:)
    end do
    cnode = 0
    do imode = 1, nm
      coord0 = center(imode, :)
      print *, 'Center:', coord0
      mindist2 = huge(mindist2)
      do i = 1, acc_mesh % numberofnodes
        coord1(1) = acc_mesh % nodes % x(i)
        coord1(2) = acc_mesh % nodes % y(i)
        coord1(3) = acc_mesh % nodes % z(i)
        dist2 = sum((coord0 - coord1)**2)
        if (dist2 < mindist2) then
          mindist2 = dist2
          mini = i
        end if
      end do
      print *, 'Minimum distance:', imode, mini, sqrt(mindist2)
      cnode(imode) = mini
    end do
  end subroutine vhl_center

  subroutine vhl_info_voltage() bind(c, name="vhl_info_voltage")
    call info(caller, 'Computing voltage using "Potential" as path indicator')
  end subroutine vhl_info_voltage

  subroutine vhl_info_cannot() bind(c, name="vhl_info_cannot")
    call info(caller, 'Cannot compute voltage as no "Potential" is present')
  end subroutine vhl_info_cannot

  ! -- EdgeVoltageIntegral (exact transcription) ------------------------
  subroutine vhl_voltage(nm, jm, nodal, inpot, volt, cnode) bind(c, name="vhl_voltage")
    integer(c_int), value :: nm, jm, nodal
    complex(c_double_complex) :: inpot(nm, nm)
    complex(c_double_complex) :: volt(nm, nm)
    integer(c_int) :: cnode(nm)
    type(matrix_t), pointer :: nodegraph
    complex :: gradv(3), circ
    real(kind=dp) :: pot, minpot, edgevector(3), s
    integer :: imode, nsteps, sgn, i, j, k, kmin, imin, i1, i2, j1, j2, l, lp, n0
    integer, pointer :: nodeindexes(:)
    type(element_t), pointer :: edge
    real(kind=dp), pointer :: potvals(:)
    integer, pointer :: potperm(:)
    nodegraph => allocatematrix()
    nodegraph % format = matrix_list
    do i = acc_mesh % numberofedges, 1, -1
      edge => acc_mesh % edges(i)
      do j = 1, edge % type % numberofnodes
        call list_addtomatrixelement(nodegraph % listmatrix, edge % nodeindexes(j), i, 1.0_dp)
      end do
    end do
    call list_tocrsmatrix(nodegraph)
    print *, 'Nonzeros per row NodeGraph:', 1.0_dp * size(nodegraph % values) / &
        nodegraph % numberofrows
    n0 = acc_mesh % numberofnodes
    potvals => acc_potvar % values
    potperm => acc_potvar % perm
    do imode = 1, nm
      i = cnode(imode)
      minpot = potvals(potperm(i))
      nsteps = 0
      circ = 0.0_dp
      print *, 'Starting pot:', imode, minpot, inpot(jm+1, imode)
      do while (.true.)
        kmin = 0
        do j = nodegraph % rows(i), nodegraph % rows(i+1) - 1
          k = nodegraph % cols(j)
          edge => acc_mesh % edges(k)
          nodeindexes => edge % nodeindexes
          do l = 1, 2
            if (nodeindexes(l) == i) cycle
            lp = potperm(nodeindexes(l))
            if (lp == 0) cycle
            pot = potvals(lp)
            if (pot < minpot) then
              kmin = k
              imin = nodeindexes(l)
              minpot = pot
            end if
          end do
        end do
        if (kmin == 0) exit
        nsteps = nsteps + 1
        k = kmin
        edge => acc_mesh % edges(k)
        i1 = edge % nodeindexes(1)
        i2 = edge % nodeindexes(2)
        edgevector(1) = acc_mesh % nodes % x(i2) - acc_mesh % nodes % x(i1)
        edgevector(2) = acc_mesh % nodes % y(i2) - acc_mesh % nodes % y(i1)
        edgevector(3) = acc_mesh % nodes % z(i2) - acc_mesh % nodes % z(i1)
        s = sqrt(sum(edgevector**2))
        sgn = 1
        if (i /= i1) sgn = -sgn
        if (nodal /= 0) then
          j1 = acc_evar % perm(i1)
          j2 = acc_evar % perm(i2)
          do k = 1, 3
            gradv(k) = cmplx(acc_evar % values(6*(j1-1)+k) + &
                acc_evar % values(6*(j2-1)+k), &
                acc_evar % values(6*(j1-1)+3+k) + &
                acc_evar % values(6*(j2-1)+3+k), kind=dp) / 2
          end do
          circ = circ + sgn * sum(gradv * edgevector)
        else
          if (parenv % pes > 1) then
            i1 = acc_mesh % parallelinfo % globaldofs(i1)
            i2 = acc_mesh % parallelinfo % globaldofs(i2)
          end if
          if (i1 < i2) sgn = -sgn
          j = acc_evar % perm(n0 + k)
          circ = circ + s * sgn * cmplx(acc_evar % values(2*j-1), &
              acc_evar % values(2*j), kind=dp)
        end if
        i = imin
      end do
      print *, 'Path integral:', imode, minpot, nsteps, circ
      volt(jm+1, imode) = circ
    end do
    call freematrix(nodegraph)
  end subroutine vhl_voltage

  ! -- result files (exact transcription) -------------------------------
  subroutine vhl_write_files(nm, has_pot, poynt, el, curr, pot, volt) &
      bind(c, name="vhl_write_files")
    integer(c_int), value :: nm, has_pot
    complex(c_double_complex) :: poynt(nm, nm), el(nm, nm), curr(nm, nm), &
        pot(nm, nm), volt(nm, nm)
    integer :: i
    open (10, file="Poynt_re.dat")
    do i = 1, nm
      write (10, *) real(poynt(i, :))
    end do
    close (10)
    open (10, file="Poynt_im.dat")
    do i = 1, nm
      write (10, *) aimag(poynt(i, :))
    end do
    close (10)
    open (10, file="El_re.dat")
    do i = 1, nm
      write (10, *) real(el(i, :))
    end do
    close (10)
    open (10, file="El_im.dat")
    do i = 1, nm
      write (10, *) aimag(el(i, :))
    end do
    close (10)
    open (10, file="Curr_re.dat")
    do i = 1, nm
      write (10, *) real(curr(i, :))
    end do
    close (10)
    open (10, file="Curr_im.dat")
    do i = 1, nm
      write (10, *) aimag(curr(i, :))
    end do
    close (10)
    if (has_pot /= 0) then
      open (10, file="dPot_re.dat")
      do i = 1, nm
        write (10, *) real(pot(i, :))
      end do
      close (10)
      open (10, file="dPot_im.dat")
      do i = 1, nm
        write (10, *) aimag(pot(i, :))
      end do
      close (10)
      open (10, file="Volt_re.dat")
      do i = 1, nm
        write (10, *) real(volt(i, :))
      end do
      close (10)
      open (10, file="Volt_im.dat")
      do i = 1, nm
        write (10, *) aimag(volt(i, :))
      end do
      close (10)
    end if
  end subroutine vhl_write_files

  subroutine vhl_info_final() bind(c, name="vhl_info_final")
    call info(caller, 'Writing results on final visit!')
  end subroutine vhl_info_final

  subroutine vhl_info_done() bind(c, name="vhl_info_done")
    call info(caller, 'All done for now!', level=20)
  end subroutine vhl_info_done

  subroutine vhl_fatal_jmode() bind(c, name="vhl_fatal_jmode")
    call fatal(caller, 'The lumping was already called "NoModes" times!')
  end subroutine vhl_fatal_jmode

  subroutine vhl_fatal_alloc() bind(c, name="vhl_fatal_alloc")
    call fatal(caller, 'Memory allocation error.')
  end subroutine vhl_fatal_alloc

end module vhl_acc
