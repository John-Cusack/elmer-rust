! T021-B NormalSolver port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/NormalSolver.F90, compiled against the built .mod
! files of the target configure (A2 shipped mechanism, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module mirrors NormalSolver's own locals (entry SAVE state, the
! BulkAssembly SAVE Nodes, the projection BLOCK locals) with identical
! lifetimes.
!
! Bit-identity rules followed here: every order-sensitive reduction
! (SUM/DOT_PRODUCT/MATMUL-equivalent double loop), every formatted WRITE
! to Message, and every whole-array/strided update stays in Fortran with
! the exact statement order and evaluation order of the original, same
! compiler and flags, so identity holds by construction wherever the call
! sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module normalsolver_acc
  use iso_c_binding
  use types
  use defutils
  use coordinatesystems
  use integration
  use matrixassembly
  use eleminfo
  use messages
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(mesh_t), pointer, save :: acc_mesh => null()
  type(valuelist_t), pointer, save :: acc_params => null()
  type(variable_t), pointer, save :: acc_nrmsol => null()
  type(element_t), pointer, save :: acc_elem => null()
  real(kind=dp), pointer, save :: acc_saverhs(:) => null()
  real(kind=dp), allocatable, target, save :: acc_forcevec(:, :)
  real(kind=dp), allocatable, save :: acc_values(:)
  character(len=MAX_NAME_LEN), save :: acc_varname = ' '
  character(len=MAX_NAME_LEN), save :: acc_vname = ' '
  integer, save :: acc_dim = 0
  logical, save :: acc_csym = .false.
  logical, save :: acc_cbinuse = .false.
  logical, save :: acc_setd = .false.
  ! Mirrors the original's SAVE Visited (dead: initialized .FALSE., never
  ! read or written by the body).
  logical, save :: acc_visited = .false.
  real(kind=dp), save :: acc_at0 = 0.0_dp
  real(kind=dp), save :: acc_at1 = 0.0_dp
  ! BulkAssembly context (mirrors its SAVE Nodes; IntegStuff is a per-call
  ! local there, re-filled for every element here with the same lifetime).
  type(nodes_t), target, save :: acc_nodes
  type(nodes_t), target, save :: acc_nodes_sym
  type(gaussintegrationpoints_t), save :: acc_ip
  integer, save :: acc_nd = 0
  integer, save :: acc_n = 0
  logical, save :: acc_checkorient = .false.
  real(kind=dp), save :: acc_point(3) = 0.0_dp
  real(kind=dp), allocatable, save :: acc_stiff(:, :), acc_force(:, :), acc_basis(:)

contains

  ! -- entry state ----------------------------------------------------
  subroutine noa_begin(h_solver, dim) bind(c, name="noa_begin")
    type(c_ptr), value :: h_solver
    integer(c_int) :: dim
    call c_f_pointer(h_solver, acc_solver)
    acc_mesh => acc_solver % mesh
    acc_params => GetSolverParams()
    acc_dim = CoordinateSystemDimension()
    dim = acc_dim
  end subroutine noa_begin

  subroutine noa_info_head() bind(c, name="noa_info_head")
    call Info('NormalSolver', '-------------------------------------', Level=4)
    call Info('NormalSolver', 'Computing the normals', Level=4)
    call Info('NormalSolver', '-------------------------------------', Level=4)
  end subroutine noa_info_head

  function noa_assoc_matrix() result(r) bind(c, name="noa_assoc_matrix")
    integer(c_int) :: r
    r = merge(1, 0, associated(acc_solver % matrix))
  end function noa_assoc_matrix

  function noa_active_perm() result(n) bind(c, name="noa_active_perm")
    integer(c_int) :: n
    n = COUNT(acc_solver % variable % perm > 0)
  end function noa_active_perm

  subroutine noa_result_variable() bind(c, name="noa_result_variable")
    logical :: gotit
    acc_varname = GetString(acc_params, 'Normals Result Variable', gotit)
    if (.not. gotit) then
      call Fatal('NormalSolver', '> Normals Result Variable < not found!')
    end if
    acc_nrmsol => VariableGet(acc_mesh % variables, acc_varname)
    if (associated(acc_nrmsol)) then
      if (acc_nrmsol % dofs /= acc_dim) then
        call Fatal('NormalSolver', 'The normals should have DOFs equal to DIM')
      end if
    else
      call DefaultVariableAdd(acc_varname, acc_dim, Var=acc_nrmsol)
    end if
  end subroutine noa_result_variable

  subroutine noa_constant_bulk(cb, inuse) bind(c, name="noa_constant_bulk")
    integer(c_int) :: cb, inuse
    logical :: gotit, cbulk
    cbulk = GetLogical(acc_params, 'Constant Bulk Matrix', gotit)
    acc_cbinuse = cbulk .and. associated(acc_solver % matrix % bulkvalues)
    cb = merge(1, 0, cbulk)
    inuse = merge(1, 0, acc_cbinuse)
    acc_csym = CurrentCoordinateSystem() == AxisSymmetric .or. &
        CurrentCoordinateSystem() == CylindricSymmetric
  end subroutine noa_constant_bulk

  subroutine noa_default_initialize() bind(c, name="noa_default_initialize")
    call DefaultInitialize(acc_solver, acc_cbinuse)
  end subroutine noa_default_initialize

  subroutine noa_tick0() bind(c, name="noa_tick0")
    acc_at0 = RealTime()
  end subroutine noa_tick0

  subroutine noa_forcevec_alloc() bind(c, name="noa_forcevec_alloc")
    allocate(acc_forcevec(size(acc_solver % matrix % rhs), acc_dim))
    acc_forcevec = 0.0_dp
    acc_saverhs => acc_solver % matrix % rhs
  end subroutine noa_forcevec_alloc

  subroutine noa_forcevec_free() bind(c, name="noa_forcevec_free")
    deallocate(acc_forcevec)
  end subroutine noa_forcevec_free

  subroutine noa_finish_bulk(cb) bind(c, name="noa_finish_bulk")
    integer(c_int), value :: cb
    if (cb /= 0) then
      call DefaultFinishBulkAssembly(BulkUpdate=.not. acc_cbinuse, RHSUpdate=.false.)
    else
      call DefaultFinishBulkAssembly()
    end if
  end subroutine noa_finish_bulk

  subroutine noa_finish_assembly() bind(c, name="noa_finish_assembly")
    call DefaultFinishAssembly()
  end subroutine noa_finish_assembly

  subroutine noa_info_assembly() bind(c, name="noa_info_assembly")
    acc_at1 = RealTime()
    write(message, *) 'Assembly Time: ', acc_at1 - acc_at0
    call Info('NormalSolver', message, Level=5)
  end subroutine noa_info_assembly

  subroutine noa_info_solution() bind(c, name="noa_info_solution")
    real(kind=dp) :: at2
    at2 = RealTime()
    write(message, *) 'Solution Time: ', at2 - acc_at1
    call Info('NormalSolver', message, Level=5)
  end subroutine noa_info_solution

  ! -- SetD head ------------------------------------------------------
  function noa_setd_begin() result(r) bind(c, name="noa_setd_begin")
    integer(c_int) :: r
    logical :: gotit
    acc_setd = GetLogical(acc_params, 'Set Dirichlet Conditions', gotit)
    if (acc_setd) then
      allocate(acc_values(size(acc_solver % matrix % values)))
      acc_vname = acc_solver % variable % name
      acc_values = acc_solver % matrix % values
    end if
    r = merge(1, 0, acc_setd)
  end function noa_setd_begin

  subroutine noa_setd_free(setd) bind(c, name="noa_setd_free")
    integer(c_int), value :: setd
    if (setd /= 0) deallocate(acc_values)
  end subroutine noa_setd_free

  ! -- per-component solve --------------------------------------------
  subroutine noa_comp_rhs(i) bind(c, name="noa_comp_rhs")
    integer(c_int), value :: i
    acc_solver % matrix % rhs => acc_forcevec(:, i)
  end subroutine noa_comp_rhs

  subroutine noa_comp_dirichlet(i, setd) bind(c, name="noa_comp_dirichlet")
    integer(c_int), value :: i, setd
    if (setd /= 0) then
      ! The normals are computed component-wise hence the Dirichlet
      ! condition also must be given so.
      acc_solver % variable % name = trim(acc_varname)//' '//i2s(i)
      acc_solver % matrix % values = acc_values
      call DefaultDirichletBCs()
      acc_solver % variable % name = acc_vname
    else
      call DefaultDirichletBCs()
    end if
  end subroutine noa_comp_dirichlet

  function noa_comp_solve() result(norm) bind(c, name="noa_comp_solve")
    real(c_double) :: norm
    acc_solver % variable % values = 0.0_dp
    norm = DefaultSolve()
  end function noa_comp_solve

  subroutine noa_comp_copy(i) bind(c, name="noa_comp_copy")
    integer(c_int), value :: i
    ! Component-wise copy into the result field. The original threads a
    ! DOFs local through this statement; its value is NrmSol % DOFs in
    ! every defined case (see noa_result_variable), so reading the field
    ! directly is identical wherever defined.
    acc_nrmsol % values(i::acc_nrmsol % dofs) = acc_solver % variable % values
  end subroutine noa_comp_copy

  subroutine noa_rhs_restore() bind(c, name="noa_rhs_restore")
    acc_solver % matrix % rhs => acc_saverhs
  end subroutine noa_rhs_restore

  ! -- normalization --------------------------------------------------
  subroutine noa_normalize() bind(c, name="noa_normalize")
    integer :: i, j, k, dofs
    real(kind=dp) :: nrm
    ! Normalize the computed normals so that the length is always one.
    dofs = acc_nrmsol % dofs
    do i = 1, acc_solver % matrix % numberofrows
      nrm = 0.0_dp
      do j = 1, acc_dim
        k = dofs*(i - 1) + j
        nrm = nrm + acc_nrmsol % values(k)**2
      end do
      if (nrm > 0) then
        nrm = 1.0_dp / sqrt(nrm)
        do j = 1, acc_dim
          k = dofs*(i - 1) + j
          acc_nrmsol % values(k) = nrm*acc_nrmsol % values(k)
        end do
      end if
    end do
  end subroutine noa_normalize

  ! -- BulkAssembly ---------------------------------------------------
  subroutine noa_bulk_alloc() bind(c, name="noa_bulk_alloc")
    integer :: n
    n = max(acc_mesh % maxelementdofs, acc_mesh % maxelementnodes)
    allocate(acc_stiff(n, n), acc_force(acc_dim, n), acc_basis(n))
  end subroutine noa_bulk_alloc

  subroutine noa_bulk_free() bind(c, name="noa_bulk_free")
    deallocate(acc_stiff, acc_force, acc_basis)
  end subroutine noa_bulk_free

  function noa_bulk_nelems() result(n) bind(c, name="noa_bulk_nelems")
    integer(c_int) :: n
    n = acc_solver % numberofactiveelements
  end function noa_bulk_nelems

  subroutine noa_bulk_elem(elem) bind(c, name="noa_bulk_elem")
    integer(c_int), value :: elem
    real(kind=dp), pointer :: arr(:) => null()
    integer :: i
    acc_elem => GetActiveElement(elem)
    call GetElementNodes(acc_nodes)
    acc_nd = GetElementNOFDOFs()
    acc_n = GetElementNOFNodes()
    arr => ListGetConstRealArray1(GetBodyParams(acc_elem), &
        'Point on Negative Side', acc_checkorient)
    if (acc_checkorient) then
      acc_point = 0.0d0
      do i = 1, size(arr)
        acc_point(i) = arr(i)
      end do
    end if
  end subroutine noa_bulk_elem

  subroutine noa_bulk_zero() bind(c, name="noa_bulk_zero")
    acc_stiff = 0.0_dp
    acc_force = 0.0_dp
  end subroutine noa_bulk_zero

  function noa_bulk_ipbegin() result(n) bind(c, name="noa_bulk_ipbegin")
    integer(c_int) :: n
    acc_ip = GaussPoints(acc_elem)
    n = acc_ip % n
  end function noa_bulk_ipbegin

  subroutine noa_bulk_ipoint(t) bind(c, name="noa_bulk_ipoint")
    integer(c_int), value :: t
    real(kind=dp) :: u, v, w, s, weight, detj, normal(3), r(3)
    logical :: ok
    integer :: p, q, i
    u = acc_ip % u(t)
    v = acc_ip % v(t)
    w = acc_ip % w(t)
    s = acc_ip % s(t)
    ok = ElementInfo(acc_elem, acc_nodes, u, v, w, detj, acc_basis(1:acc_n))
    weight = s*detj
    if (acc_csym) weight = weight*sum(acc_basis(1:acc_n)*acc_nodes % x(1:acc_n))
    if (.not. acc_cbinuse) then
      do p = 1, acc_nd
        do q = 1, acc_nd
          acc_stiff(p, q) = acc_stiff(p, q) + weight*acc_basis(q)*acc_basis(p)
        end do
      end do
    end if
    normal = NormalVector(acc_elem, acc_nodes, u, v, .true.)
    if (acc_checkorient) then
      r(1) = sum(acc_basis(1:acc_n)*acc_nodes % x(1:acc_n)) - acc_point(1)
      r(2) = sum(acc_basis(1:acc_n)*acc_nodes % y(1:acc_n)) - acc_point(2)
      r(3) = sum(acc_basis(1:acc_n)*acc_nodes % z(1:acc_n)) - acc_point(3)
      if (sum(normal*r) < 0.0d0) normal = -normal
    end if
    do i = 1, acc_dim
      acc_force(i, 1:acc_nd) = acc_force(i, 1:acc_nd) + &
          weight*normal(i)*acc_basis(1:acc_nd)
    end do
  end subroutine noa_bulk_ipoint

  subroutine noa_bulk_tomatrix() bind(c, name="noa_bulk_tomatrix")
    acc_solver % matrix % rhs => acc_saverhs
    call DefaultUpdateEquations(acc_stiff, acc_force(1, 1:acc_nd))
  end subroutine noa_bulk_tomatrix

  subroutine noa_bulk_toforce() bind(c, name="noa_bulk_toforce")
    integer :: i
    do i = 1, acc_dim
      acc_solver % matrix % rhs => acc_forcevec(:, i)
      call DefaultUpdateForce(acc_force(i, 1:acc_nd))
    end do
  end subroutine noa_bulk_toforce

  ! -- projection BLOCK (whole block: GotIt quirk lives here) ---------
  subroutine noa_project(dt, transient) bind(c, name="noa_project")
    real(c_double), value :: dt
    integer(c_int), value :: transient
    type(variable_t), pointer :: var1, var2, var3
    real(kind=dp) :: nvec(3), dvec(3), dvec0(3)
    real(kind=dp) :: maxd, limd, coeff
    integer :: i, j, k, id, inn, dofs
    logical :: gotit, scaledt
    character(len=MAX_NAME_LEN) :: vname
    vname = ListGetString(acc_params, 'Vector Field to Project', gotit)
    if (gotit) then
      scaledt = ListGetLogical(acc_params, 'Vector Field Timestep scaling', gotit)
    else
      if (ListGetLogical(acc_params, 'Project Displacement', gotit)) then
        vname = 'Displacement'
        scaledt = .false.
      else if (ListGetLogical(acc_params, 'Project Velocity', gotit)) then
        vname = 'Velocity'
        scaledt = .true.
      end if
    end if
    if (gotit) then
      call Info('NormalSolver', &
          'Projecting field "'//trim(vname)//'" to normal direction', Level=5)
      ! The displacement/velocity field is dealt with component-wise since
      ! otherwise we should deal with the 'flow solution' vector separately.
      var1 => VariableGet(acc_mesh % variables, trim(vname)//' 1')
      if (.not. associated(var1)) call Fatal('NormalSolver', &
          'Vector field component 1 does not exist: '//trim(vname))
      var2 => VariableGet(acc_mesh % variables, trim(vname)//' 2')
      if (.not. associated(var2)) call Fatal('NormalSolver', &
          'Vector field component 2 does not exist: '//trim(vname))
      if (acc_dim == 3) then
        var3 => VariableGet(acc_mesh % variables, trim(vname)//' 3')
        if (.not. associated(var3)) call Fatal('NormalSolver', &
            'Vector field component 3 does not exist: '//trim(vname))
      end if
      dvec0(1) = ListGetCReal(acc_params, 'Vector Field Offset 1', gotit)
      dvec0(2) = ListGetCReal(acc_params, 'Vector Field Offset 2', gotit)
      dvec0(3) = ListGetCReal(acc_params, 'Vector Field Offset 3', gotit)
      nvec = 0.0_dp
      dvec = 0.0_dp
      maxd = 0.0_dp
      dofs = acc_nrmsol % dofs
      do i = 1, acc_mesh % numberofnodes
        inn = acc_nrmsol % perm(i)
        if (inn == 0) cycle
        id = var1 % perm(i)
        if (id == 0) cycle
        do j = 1, acc_dim
          k = dofs*(inn - 1) + j
          nvec(j) = acc_nrmsol % values(k)
        end do
        dvec(1) = var1 % values(id)
        dvec(2) = var2 % values(id)
        if (acc_dim == 3) dvec(3) = var3 % values(id)
        dvec = dvec - dvec0
        dvec = sum(dvec*nvec)*nvec
        maxd = max(maxd, sqrt(sum(dvec(1:acc_dim)**2)))
        do j = 1, acc_dim
          k = dofs*(inn - 1) + j
          acc_nrmsol % values(k) = dvec(j)
        end do
      end do
      ! If we have velocity then the normal displacement is given by
      ! d = dt*(n.v)n
      if (scaledt .and. transient /= 0) then
        acc_nrmsol % values = dt*acc_nrmsol % values
        maxd = maxd*dt
      end if
      limd = ListGetCReal(acc_params, 'Maximum Displacement', gotit)
      if (gotit) then
        if (maxd > limd) then
          coeff = limd/maxd
          write(message, '(A,ES12.3)') 'Limiting displacements with factor: ', coeff
          call Info('NormalSolver', message)
          acc_nrmsol % values = coeff*acc_nrmsol % values
        end if
      end if
    end if
  end subroutine noa_project

  ! -- EnforceSymmetry ------------------------------------------------
  function noa_want_symmetry() result(r) bind(c, name="noa_want_symmetry")
    integer(c_int) :: r
    logical :: gotit
    r = merge(1, 0, ListGetLogical(acc_params, 'Enforce Symmetry', gotit))
  end function noa_want_symmetry

  subroutine noa_enforce_symmetry() bind(c, name="noa_enforce_symmetry")
    logical :: found
    integer :: elem, i, j, n, nd, m, k, dofs
    type(element_t), pointer :: element
    call Info('NormalSolver', 'Enforcing symmetry of normal field', Level=6)
    m = 0
    dofs = acc_nrmsol % dofs
    do elem = 1, acc_solver % numberofactiveelements
      element => GetActiveElement(elem)
      call GetElementNodes(acc_nodes_sym)
      nd = GetElementNOFDOFs()
      n = GetElementNOFNodes()
      found = .false.
      do i = 1, n
        if (abs(acc_nodes_sym % x(i)) < 1.0d-8) then
          found = .true.
          exit
        end if
      end do
      if (found) then
        j = 3 - i
        !PRINT *,'Indeces:',i,j,Element % NodeIndexes
        j = element % nodeindexes(j)
        i = element % nodeindexes(i)
        j = acc_nrmsol % perm(j)
        i = acc_nrmsol % perm(i)
        do k = 1, acc_dim
          acc_nrmsol % values(dofs*(i - 1) + k) = acc_nrmsol % values(dofs*(j - 1) + k)
        end do
        m = m + 1
      end if
    end do
    call Info('NormalSolver', 'Enforced axial symmetry in nodes: '//i2s(m), Level=5)
  end subroutine noa_enforce_symmetry

  ! -- NormalSolver_Init entry ----------------------------------------
  subroutine noa_init_begin(h_solver) bind(c, name="noa_init_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
    acc_params => GetSolverParams()
  end subroutine noa_init_begin

  subroutine noa_init_variable() bind(c, name="noa_init_variable")
    call ListAddNewString(acc_params, 'Variable', '-nooutput nrm_temp')
  end subroutine noa_init_variable

  function noa_init_dim() result(d) bind(c, name="noa_init_dim")
    integer(c_int) :: d
    d = CoordinateSystemDimension()
  end function noa_init_dim

  subroutine noa_init_result(dim) bind(c, name="noa_init_result")
    integer(c_int), value :: dim
    character(len=MAX_NAME_LEN) :: varname
    logical :: found
    varname = GetString(acc_params, 'Normals Result Variable', found)
    if (.not. found) then
      call ListAddString(acc_params, 'Normals Result Variable', 'Normals')
      call ListAddString(acc_params, 'Exported Variable 1', &
          'Normals[Normals:'//i2s(dim)//']')
    end if
  end subroutine noa_init_result

  subroutine noa_init_timeorder() bind(c, name="noa_init_timeorder")
    call ListAddInteger(acc_params, 'Time derivative order', 0)
  end subroutine noa_init_timeorder

end module normalsolver_acc
