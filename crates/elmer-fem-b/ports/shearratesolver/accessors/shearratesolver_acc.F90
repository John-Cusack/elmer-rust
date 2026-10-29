! ShearrateSolver port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/ShearrateSolver.F90, compiled against the
! built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars). All Fortran-object traffic — TYPE(...)
! handles, OPTIONAL/POINTER dummies, assumed-shape descriptors, hidden
! lengths, array-pointer results — stays inside Fortran, so the Rust side
! is pure C calls. The wrapper module mirrors ShearrateSolver's own
! SAVE-carrying locals (entry-level ViscVector/AllocationsDone, the
! Force/Shearrate/Visc pointer binds, dim/branch flags, and the
! BulkAssembly internal procedure's SAVE Nodes + SAVE arrays) with
! identical lifetimes.
!
! BulkAssembly is transcribed verbatim as sra_bulk_assembly (same
! compiler, same evaluation order, same whole-array update shapes) because
! its string-keyed lookups (ComponentName, VariableGet) and
! MaterialModels calls cannot cross the C boundary; the Rust entry drives
! every branch around it in original order.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module shearratesolver_acc
  use iso_c_binding
  use CoordinateSystems
  use DefUtils
  use MaterialModels
  implicit none
  type(solver_t), pointer, save :: acc_solver => null()
  type(variable_t), pointer, save :: acc_shearratesol => null()
  type(variable_t), pointer, save :: acc_viscsol => null()
  ! Entry-level SAVE mirrors (host scope of ShearrateSolver).
  real(kind=dp), pointer, save :: acc_viscvector(:) => null()
  real(kind=dp), pointer, save :: acc_forcevector(:) => null()
  real(kind=dp), pointer, save :: acc_shearratefield(:) => null()
  real(kind=dp), pointer, save :: acc_viscfield(:) => null()
  logical, save :: acc_allocdone = .false.
  integer, save :: acc_dim = 0
  logical, save :: acc_calcvisc = .false.
  logical, save :: acc_constinuse = .false.
  ! BulkAssembly SAVE mirrors (internal-procedure SAVE).
  type(nodes_t), save :: acc_nodes
  real(kind=dp), allocatable, save :: acc_stiff(:,:), acc_force(:), acc_force2(:)
  real(kind=dp), allocatable, save :: acc_vx(:), acc_vy(:), acc_vz(:)
  real(kind=dp), allocatable, save :: acc_basis(:), acc_dbasisdx(:,:)
  real(kind=dp), allocatable, save :: acc_nodalrho(:), acc_nodalvisc(:)
  logical, save :: acc_ballocdone = .false.

contains

  ! -- init entry state -------------------------------------------------
  subroutine sra_init_begin(h_solver) bind(c, name="sra_init_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine sra_init_begin

  function sra_dim() result(n) bind(c, name="sra_dim")
    integer(c_int) :: n
    n = CoordinateSystemDimension()
  end function sra_dim

  function sra_listcheck_variable() result(r) bind(c, name="sra_listcheck_variable")
    integer(c_int) :: r
    type(valuelist_t), pointer :: sp
    sp => GetSolverParams()
    r = merge(1, 0, ListCheckPresent(sp, 'Variable'))
  end function sra_listcheck_variable

  subroutine sra_add_variable_shearrate() bind(c, name="sra_add_variable_shearrate")
    type(valuelist_t), pointer :: sp
    character(len=MAX_NAME_LEN) :: varname, fieldname
    logical :: gotit
    sp => GetSolverParams()
    varname = GetString(sp, 'Target Variable', gotit)
    if (.not. gotit) varname = trim('Velocity')
    fieldname = 'Shearrate '//trim(varname)
    call ListAddString(sp, 'Variable', fieldname)
  end subroutine sra_add_variable_shearrate

  function sra_get_calcvisc() result(r) bind(c, name="sra_get_calcvisc")
    integer(c_int) :: r
    type(valuelist_t), pointer :: sp
    logical :: gotit
    sp => GetSolverParams()
    acc_calcvisc = GetLogical(sp, 'Calculate Viscosity', gotit)
    r = merge(1, 0, acc_calcvisc)
  end function sra_get_calcvisc

  subroutine sra_add_exported_viscosity() bind(c, name="sra_add_exported_viscosity")
    type(valuelist_t), pointer :: sp
    character(len=MAX_NAME_LEN) :: varname, fieldname
    logical :: gotit
    sp => GetSolverParams()
    varname = GetString(sp, 'Target Variable', gotit)
    if (.not. gotit) varname = trim('Velocity')
    fieldname = 'Viscosity '//trim(varname)
    call ListAddString(sp, NextFreeKeyword('Exported Variable', sp), fieldname)
  end subroutine sra_add_exported_viscosity

  subroutine sra_init_tail() bind(c, name="sra_init_tail")
    type(valuelist_t), pointer :: sp
    sp => GetSolverParams()
    call ListAddInteger(sp, 'Time derivative order', 0)
    call ListAddNewString(sp, 'Linear System Solver', 'Iterative')
    call ListAddNewString(sp, 'Linear System Iterative Method', 'cg')
    call ListAddNewString(sp, 'Linear System Preconditioning', 'ILU0')
    call ListAddNewInteger(sp, 'Linear System Max Iterations', 500)
    call ListAddNewConstReal(sp, 'Linear System Convergence Tolerance', 1.0e-8_dp)
  end subroutine sra_init_tail

  ! -- solver entry state -----------------------------------------------
  subroutine sra_begin(h_solver) bind(c, name="sra_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
    acc_dim = CoordinateSystemDimension()
  end subroutine sra_begin

  subroutine sra_info_banner() bind(c, name="sra_info_banner")
    call Info('ShearrateSolver', '-------------------------------------', level=4)
    call Info('ShearrateSolver', 'Computing the shearrate field', level=4)
    call Info('ShearrateSolver', '-------------------------------------', level=4)
  end subroutine sra_info_banner

  function sra_assoc_matrix() result(r) bind(c, name="sra_assoc_matrix")
    integer(c_int) :: r
    r = merge(1, 0, associated(acc_solver % matrix))
  end function sra_assoc_matrix

  function sra_active_perm_count() result(n) bind(c, name="sra_active_perm_count")
    integer(c_int) :: n
    n = count(acc_solver % variable % perm > 0)
  end function sra_active_perm_count

  function sra_var_dofs() result(d) bind(c, name="sra_var_dofs")
    integer(c_int) :: d
    acc_shearratesol => acc_solver % variable
    d = acc_shearratesol % dofs
  end function sra_var_dofs

  subroutine sra_fatal_dofs() bind(c, name="sra_fatal_dofs")
    call Fatal('ShearrateSolver', 'Shearrate should have just 1 component')
  end subroutine sra_fatal_dofs

  function sra_rhs_size() result(n) bind(c, name="sra_rhs_size")
    integer(c_int) :: n
    n = size(acc_solver % matrix % rhs)
  end function sra_rhs_size

  subroutine sra_visc_alloc_maybe() bind(c, name="sra_visc_alloc_maybe")
    if (.not. acc_allocdone) then
      allocate(acc_viscvector(size(acc_solver % matrix % rhs)))
      acc_allocdone = .true.
    end if
  end subroutine sra_visc_alloc_maybe

  subroutine sra_visc_zero() bind(c, name="sra_visc_zero")
    acc_viscvector = 0.0_dp
  end subroutine sra_visc_zero

  subroutine sra_visc_field_bind() bind(c, name="sra_visc_field_bind")
    type(valuelist_t), pointer :: sp
    character(len=MAX_NAME_LEN) :: varname
    logical :: gotit
    sp => GetSolverParams()
    varname = GetString(sp, 'Target Variable', gotit)
    if (.not. gotit) varname = trim('Velocity')
    acc_viscsol => VariableGet(acc_solver % mesh % variables, &
        'Viscosity '//trim(varname))
    acc_viscfield => acc_viscsol % values
  end subroutine sra_visc_field_bind

  function sra_get_constbulk() result(r) bind(c, name="sra_get_constbulk")
    integer(c_int) :: r
    type(valuelist_t), pointer :: sp
    logical :: gotit, cb
    sp => GetSolverParams()
    cb = GetLogical(sp, 'Constant Bulk Matrix', gotit)
    r = merge(1, 0, cb)
  end function sra_get_constbulk

  function sra_bind_force_field(const_bulk) result(in_use) bind(c, name="sra_bind_force_field")
    integer(c_int), value :: const_bulk
    integer(c_int) :: in_use
    acc_constinuse = (const_bulk /= 0) .and. &
        associated(acc_solver % matrix % bulkvalues)
    acc_forcevector => acc_solver % matrix % rhs
    acc_shearratefield => acc_solver % variable % values
    in_use = merge(1, 0, acc_constinuse)
  end function sra_bind_force_field

  subroutine sra_default_initialize(use_const_bulk) bind(c, name="sra_default_initialize")
    integer(c_int), value :: use_const_bulk
    call DefaultInitialize(acc_solver, use_const_bulk /= 0)
  end subroutine sra_default_initialize

  ! -- BulkAssembly: verbatim transcription; host-associated names read
  ! -- the module SAVE mirrors set by the entry wrappers above.
  subroutine sra_bulk_assembly() bind(c, name="sra_bulk_assembly")
    integer :: elem, t, i, j, p, q, n, nd
    type(GaussIntegrationPoints_t), target :: integstuff
    type(element_t), pointer :: element
    real(kind=dp) :: sqrtmetric, metric(3,3), symb(3,3,3), dsymb(3,3,3,3)
    real(kind=dp) :: detj, shearrate, velo(3), dvelodx(3,3), s, x, y, z, &
        u, v, w, viscatip, rhoatip, mu, muder0
    logical :: found
    type(valuelist_t), pointer :: material, solverparams
    type(variable_t), pointer :: shearratesol
    character(len=MAX_NAME_LEN) :: varname
    logical :: gotit
    integer :: dim
    logical :: calculateviscosity, constantbulkmatrixinuse

    dim = acc_dim
    calculateviscosity = acc_calcvisc
    constantbulkmatrixinuse = acc_constinuse
    solverparams => GetSolverParams()
    shearratesol => acc_solver % variable

    varname = GetString(solverparams, 'Target Variable', gotit)
    if (.not. gotit) varname = trim('Velocity')

    if (.not. acc_ballocdone) then
      n = acc_solver % mesh % maxelementnodes
      allocate(acc_stiff(n,n), acc_force(n), acc_force2(n), acc_vx(n), &
          acc_vy(n), acc_vz(n), acc_basis(n), acc_nodalvisc(n), &
          acc_nodalrho(n), acc_dbasisdx(n,3))
      acc_vx = 0.0_dp
      acc_vy = 0.0_dp
      acc_vz = 0.0_dp
      acc_ballocdone = .true.
    end if

    do elem = 1, acc_solver % numberofactiveelements
      element => GetActiveElement(elem)
      call GetElementNodes(acc_nodes)
      nd = GetElementNOFDOFs()
      n = GetElementNOFNodes()

      if (calculateviscosity) then
        material => GetMaterial()
        acc_nodalvisc(1:n) = GetReal(material, 'Viscosity')
        acc_nodalrho(1:n) = GetReal(material, 'Density')
      end if

      call GetScalarLocalSolution(acc_vx, ComponentName(varname, 1))
      call GetScalarLocalSolution(acc_vy, ComponentName(varname, 2))
      if (dim > 2) then
        call GetScalarLocalSolution(acc_vz, ComponentName(varname, 3))
      end if

      integstuff = GaussPoints(element)
      acc_stiff = 0.0_dp
      acc_force = 0.0_dp
      acc_force2 = 0.0_dp

      do t = 1, integstuff % n
        u = integstuff % u(t)
        v = integstuff % v(t)
        w = integstuff % w(t)

        found = ElementInfo(element, acc_nodes, u, v, w, detj, &
            acc_basis, acc_dbasisdx)

        x = sum(acc_nodes % x(1:n) * acc_basis(1:n))
        y = sum(acc_nodes % y(1:n) * acc_basis(1:n))
        z = sum(acc_nodes % z(1:n) * acc_basis(1:n))
        call CoordinateSystemInfo(metric, sqrtmetric, symb, dsymb, x, y, z)

        s = detj * sqrtmetric * integstuff % s(t)

        if (calculateviscosity) then
          rhoatip = sum(acc_nodalrho(1:n) * acc_basis(1:n))
          viscatip = sum(acc_nodalvisc(1:n) * acc_basis(1:n))
          if (ListCheckPresent(material, 'Viscosity Model')) then
            mu = EffectiveViscosity(viscatip, rhoatip, acc_vx, acc_vy, acc_vz, &
                element, acc_nodes, n, n, u, v, w, muder0, LocalIP=t)
          else
            mu = viscatip
          end if
        end if

        if (.not. constantbulkmatrixinuse) then
          do p = 1, nd
            do q = 1, nd
              acc_stiff(p,q) = acc_stiff(p,q) + s * acc_basis(q) * acc_basis(p)
            end do
          end do
        end if

        do j = 1, 3
          dvelodx(1,j) = sum(acc_vx(1:nd) * acc_dbasisdx(1:nd,j))
          dvelodx(2,j) = sum(acc_vy(1:nd) * acc_dbasisdx(1:nd,j))
          dvelodx(3,j) = sum(acc_vz(1:nd) * acc_dbasisdx(1:nd,j))
        end do

        velo(1) = sum(acc_basis(1:nd) * acc_vx(1:nd))
        velo(2) = sum(acc_basis(1:nd) * acc_vy(1:nd))
        velo(3) = sum(acc_basis(1:nd) * acc_vz(1:nd))

        shearrate = sqrt(SecondInvariant(velo, dvelodx, metric, symb) / 2)

        acc_force(1:nd) = acc_force(1:nd) + acc_basis(1:nd) * s * shearrate

        if (calculateviscosity) then
          acc_force2(1:nd) = acc_force2(1:nd) + acc_basis(1:nd) * s * mu
        end if
      end do

      if (.not. constantbulkmatrixinuse) then
        call DefaultUpdateEquations(acc_stiff, acc_force(1:nd))
      else
        call DefaultUpdateForce(acc_force(1:nd))
      end if

      if (calculateviscosity) then
        acc_solver % matrix % rhs => acc_viscvector
        call DefaultUpdateForce(acc_force2)
        acc_solver % matrix % rhs => acc_forcevector
      end if
    end do
  end subroutine sra_bulk_assembly

  subroutine sra_finish_bulk(const_bulk, in_use) bind(c, name="sra_finish_bulk")
    integer(c_int), value :: const_bulk, in_use
    if (const_bulk /= 0) then
      call DefaultFinishBulkAssembly(BulkUpdate = .not. (in_use /= 0), &
          RHSUpdate = .false.)
    else
      call DefaultFinishBulkAssembly()
    end if
  end subroutine sra_finish_bulk

  subroutine sra_default_finish() bind(c, name="sra_default_finish")
    call DefaultFinishAssembly()
  end subroutine sra_default_finish

  subroutine sra_info_visc() bind(c, name="sra_info_visc")
    call Info('ShearrateSolver', 'Solving for viscosity', level=5)
  end subroutine sra_info_visc

  subroutine sra_info_shear() bind(c, name="sra_info_shear")
    call Info('ShearrateSolver', 'Solving for shearrare', level=5)
  end subroutine sra_info_shear

  subroutine sra_swap_to_visc() bind(c, name="sra_swap_to_visc")
    acc_solver % matrix % rhs => acc_viscvector
    acc_solver % variable % values => acc_viscfield
  end subroutine sra_swap_to_visc

  subroutine sra_swap_back() bind(c, name="sra_swap_back")
    acc_solver % variable % values => acc_shearratefield
    acc_solver % matrix % rhs => acc_forcevector
  end subroutine sra_swap_back

  function sra_default_solve() result(norm) bind(c, name="sra_default_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function sra_default_solve

  subroutine sra_fatal_alloc() bind(c, name="sra_fatal_alloc")
    call Fatal('ShearrateSolver', 'Memory allocation error.')
  end subroutine sra_fatal_alloc

end module shearratesolver_acc
