! T021-B RigidMeshMapper port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/RigidMeshMapper.F90, compiled against
! the built .mod files of the target configure (ADR-057 §3, T020 Poisson
! pattern).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls. Order-sensitive
! reductions (MATMUL, SUM), formatted I/O (WRITE), and whole-array updates
! stay in these wrappers, so bit-identity holds by construction wherever
! the call sequence matches.
!
! SAVE mirrors: the original's SAVE state (Visited-adjacent pointers,
! RelaxField/RelaxPerm, RotorElement sizing, main-procedure Nodes, the
! LocalMatrix internal-procedure SAVE Nodes/Hwrk) lives here with
! identical lifetimes; the call sequence is driven from Rust exactly as
! in the original. Main-procedure 4x4 mapping matrices and per-node
! Angles/Scaling/dCoord/Origin are recomputed per node exactly as the
! original's locals are (overwritten before every use in the same order).
!
! Faithfulness notes (constraint transcriptions, same compiler/flags):
! - `ElemYoungs(i,i:n) = Hwrk(1,1,1:n)` is transcribed verbatim.
! - The original recomputes `GotMatrix = GotRotate .OR. GotScale` after
!   the Mesh Matrix fast path, so a lone `Mesh Matrix` array is a
!   (defined) per-node no-op; the port replicates the reassignment.
! - With `Mesh Scale` but no rotation the original reads its uninit
!   `RotMatrix` local (undefined); proofs keep to defined inputs (scale
!   only together with rotation), documented in the Rust binder.
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module rigidmeshmapper_acc
  use iso_c_binding
  use types
  use defutils
  use lists
  use coordinatesystems
  use meshextrusion
  use meshtransform
  use integration
  use eleminfo
  use generalutils
  use parallelutils
  use modeldescription
  use messages
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(mesh_t), pointer, save :: acc_mesh => null()
  type(valuelist_t), pointer, save :: acc_solverparams => null()
  type(valuelist_t), pointer, save :: acc_vlist => null()
  type(valuelist_t), pointer, save :: acc_prevlist => null()
  type(element_t), pointer, save :: acc_elem => null()
  integer, save :: acc_nidx(1) = 0
  ! Main-procedure SAVE Nodes (rotor-marking centroid).
  type(nodes_t), target, save :: acc_nodes_main
  ! LocalMatrix internal-procedure SAVE Nodes / Hwrk.
  type(nodes_t), target, save :: acc_nodes_lm
  type(gaussintegrationpoints_t), save :: acc_ip_lm
  real(kind=dp), pointer, save :: acc_hwrk(:,:,:) => null()
  ! Main-procedure SAVE PArray / RotorBodies / IntArray (re-pointed
  ! before every use, as in the original).
  real(kind=dp), pointer, save :: acc_parray(:,:) => null()
  integer, pointer, save :: acc_rotorbodies(:) => null()
  logical, save :: acc_have_rotorbodies = .false.
  real(kind=dp), save :: acc_rotorrad = 0.0_dp
  integer, pointer, save :: acc_intarray(:) => null()
  ! Main-procedure SAVE RelaxField / RelaxPerm.
  real(kind=dp), pointer, save :: acc_relaxfield(:) => null()
  integer, pointer, save :: acc_relaxperm(:) => null()
  ! Mesh-velocity field created on demand.
  real(kind=dp), pointer, save :: acc_veloval(:) => null()
  integer, pointer, save :: acc_velperm(:) => null()
  ! Per-node mapping scratch (original's 4x4 locals + vectors).
  real(kind=dp), save :: acc_trf(4,4) = 0.0_dp
  real(kind=dp), save :: acc_rot(4,4) = 0.0_dp
  real(kind=dp), save :: acc_angles(3) = 0.0_dp
  real(kind=dp), save :: acc_scaling(3) = 0.0_dp
  real(kind=dp), save :: acc_dcoord(3) = 0.0_dp
  real(kind=dp), save :: acc_origin(4) = 0.0_dp
  integer, save :: acc_rotorder(3) = (/ 1, 2, 3 /)
  logical, save :: acc_use_orig = .true.

contains

  ! -- entry state ----------------------------------------------------
  subroutine rmm_begin(h_model, h_solver) bind(c, name="rmm_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine rmm_begin

  subroutine rmm_solver_params() bind(c, name="rmm_solver_params")
    acc_solverparams => GetSolverParams()
  end subroutine rmm_solver_params

  subroutine rmm_target_mesh_index(k, found) bind(c, name="rmm_target_mesh_index")
    integer(c_int) :: k, found
    logical :: f
    k = ListGetInteger(acc_solverparams, 'Target Mesh Solver Index', f)
    found = merge(1, 0, f)
  end subroutine rmm_target_mesh_index

  subroutine rmm_select_mesh(use_target, k) bind(c, name="rmm_select_mesh")
    integer(c_int), value :: use_target, k
    if (use_target /= 0) then
      acc_mesh => acc_model % solvers(k) % mesh
    else
      acc_mesh => acc_solver % mesh
    end if
  end subroutine rmm_select_mesh

  function rmm_assoc_mesh() result(r) bind(c, name="rmm_assoc_mesh")
    integer(c_int) :: r
    r = merge(1, 0, associated(acc_mesh))
  end function rmm_assoc_mesh

  subroutine rmm_fatal_no_mesh() bind(c, name="rmm_fatal_no_mesh")
    call Fatal('RigidMeshMapper', 'No mesh associated, cannot continue!')
  end subroutine rmm_fatal_no_mesh

  function rmm_coord_dim() result(d) bind(c, name="rmm_coord_dim")
    integer(c_int) :: d
    d = CoordinateSystemDimension()
  end function rmm_coord_dim

  ! which: 1 Cumulative Displacements, 2 Store Original Coordinates,
  ! 3 Translate Before Rotate, 4 Calculate Mesh Velocity
  ! (all DefUtils GetLogical, as original).
  function rmm_get_logical(which, found) result(l) bind(c, name="rmm_get_logical")
    integer(c_int), value :: which
    integer(c_int) :: found, l
    logical :: f, v
    if (which == 1) then
      v = GetLogical(acc_solverparams, 'Cumulative Displacements', f)
    else if (which == 2) then
      v = GetLogical(acc_solverparams, 'Store Original Coordinates', f)
    else if (which == 3) then
      v = GetLogical(acc_solverparams, 'Translate Before Rotate', f)
    else
      v = GetLogical(acc_solverparams, 'Calculate Mesh Velocity', f)
    end if
    l = merge(1, 0, v)
    found = merge(1, 0, f)
  end function rmm_get_logical

  ! which: 1 Whole Mesh Mode, 2 Rotor Mode, 3 Constant Mapping,
  ! 4 mmg remesh, 5 Mesh Relax Normalize (all ListGetLogical).
  function rmm_list_logical(which) result(l) bind(c, name="rmm_list_logical")
    integer(c_int), value :: which
    integer(c_int) :: l
    logical :: f
    if (which == 1) then
      l = merge(1, 0, ListGetLogical(acc_solverparams, 'Whole Mesh Mode', f))
    else if (which == 2) then
      l = merge(1, 0, ListGetLogical(acc_solverparams, 'Rotor Mode', f))
    else if (which == 3) then
      l = merge(1, 0, ListGetLogical(acc_solverparams, 'Constant Mapping', f))
    else if (which == 4) then
      l = merge(1, 0, ListGetLogical(acc_solverparams, 'mmg remesh', f))
    else
      l = merge(1, 0, ListGetLogical(acc_solverparams, 'Mesh Relax Normalize', f))
    end if
  end function rmm_list_logical

  function rmm_angle_coeff() result(c) bind(c, name="rmm_angle_coeff")
    real(c_double) :: c
    logical :: f
    if (ListGetLogical(CurrentModel % Simulation, 'Rotate in Radians', f)) then
      c = 1.0_dp
    else
      c = PI / 180.0_dp
    end if
  end function rmm_angle_coeff

  function rmm_no_nodes() result(n) bind(c, name="rmm_no_nodes")
    integer(c_int) :: n
    n = acc_mesh % numberofnodes
  end function rmm_no_nodes

  function rmm_num_bulk() result(n) bind(c, name="rmm_num_bulk")
    integer(c_int) :: n
    n = acc_mesh % numberofbulkelements
  end function rmm_num_bulk

  function rmm_max_elem_nodes() result(n) bind(c, name="rmm_max_elem_nodes")
    integer(c_int) :: n
    n = acc_mesh % maxelementnodes
  end function rmm_max_elem_nodes

  subroutine rmm_set_use_orig(flag) bind(c, name="rmm_set_use_orig")
    integer(c_int), value :: flag
    acc_use_orig = (flag /= 0)
  end subroutine rmm_set_use_orig

  ! -- rotor setup ----------------------------------------------------
  function rmm_rotor_bodies() result(r) bind(c, name="rmm_rotor_bodies")
    integer(c_int) :: r
    logical :: f
    acc_rotorbodies => ListGetIntegerArray(acc_solverparams, 'Rotor Bodies', f)
    acc_have_rotorbodies = associated(acc_rotorbodies)
    r = merge(1, 0, acc_have_rotorbodies)
  end function rmm_rotor_bodies

  function rmm_rotor_radius() result(x) bind(c, name="rmm_rotor_radius")
    real(c_double) :: x
    logical :: f
    acc_rotorrad = ListGetCReal(CurrentModel % Simulation, 'Rotor Radius', f)
    if (.not. f) then
      call Fatal('RigidMeshMapper', &
        'In "Rotor Mode" give either "Rotor Radius" or "Rotor Bodies"!')
    end if
    x = acc_rotorrad
  end function rmm_rotor_radius

  subroutine rmm_rotor_angle(a, got) bind(c, name="rmm_rotor_angle")
    real(c_double) :: a
    integer(c_int) :: got
    logical :: f
    a = ListGetCReal(CurrentModel % Simulation, 'Rotor Angle', f)
    got = merge(1, 0, f)
  end subroutine rmm_rotor_angle

  ! -- mesh node coordinate access (1-based node numbers) -------------
  function rmm_get_xnew(i) result(x) bind(c, name="rmm_get_xnew")
    integer(c_int), value :: i
    real(c_double) :: x
    x = acc_mesh % nodes % x(i)
  end function rmm_get_xnew

  function rmm_get_ynew(i) result(x) bind(c, name="rmm_get_ynew")
    integer(c_int), value :: i
    real(c_double) :: x
    x = acc_mesh % nodes % y(i)
  end function rmm_get_ynew

  function rmm_get_znew(i) result(x) bind(c, name="rmm_get_znew")
    integer(c_int), value :: i
    real(c_double) :: x
    x = acc_mesh % nodes % z(i)
  end function rmm_get_znew

  function rmm_get_xorig(i) result(x) bind(c, name="rmm_get_xorig")
    integer(c_int), value :: i
    real(c_double) :: x
    if (acc_use_orig) then
      x = acc_mesh % nodesorig % x(i)
    else
      x = acc_mesh % nodes % x(i)
    end if
  end function rmm_get_xorig

  function rmm_get_yorig(i) result(x) bind(c, name="rmm_get_yorig")
    integer(c_int), value :: i
    real(c_double) :: x
    if (acc_use_orig) then
      x = acc_mesh % nodesorig % y(i)
    else
      x = acc_mesh % nodes % y(i)
    end if
  end function rmm_get_yorig

  function rmm_get_zorig(i) result(x) bind(c, name="rmm_get_zorig")
    integer(c_int), value :: i
    real(c_double) :: x
    if (acc_use_orig) then
      x = acc_mesh % nodesorig % z(i)
    else
      x = acc_mesh % nodes % z(i)
    end if
  end function rmm_get_zorig

  ! -- rotor marking (setup-loop body, verbatim) ----------------------
  subroutine rmm_set_mesh_skew() bind(c, name="rmm_set_mesh_skew")
    call SetMeshSkew(acc_mesh, CurrentModel % Simulation)
  end subroutine rmm_set_mesh_skew

  subroutine rmm_rotor_mark(elem, n, isrotor) bind(c, name="rmm_rotor_mark")
    integer(c_int), value :: elem
    integer(c_int) :: n, isrotor
    integer :: nn
    real(kind=dp) :: coord(3)
    logical :: isr
    acc_elem => acc_mesh % elements(elem)
    nn = GetElementNOFNodes(acc_elem)
    n = nn
    call GetElementNodes(acc_nodes_main, acc_elem)
    coord(1) = SUM(acc_nodes_main % x(1:nn)) / nn
    coord(2) = SUM(acc_nodes_main % y(1:nn)) / nn
    coord(3) = SUM(acc_nodes_main % z(1:nn)) / nn
    if (acc_have_rotorbodies) then
      isr = ANY(acc_rotorbodies == acc_elem % bodyid)
    else
      isr = (coord(1)**2 + coord(2)**2 < acc_rotorrad**2)
    end if
    isrotor = merge(1, 0, isr)
  end subroutine rmm_rotor_mark

  subroutine rmm_store_orig() bind(c, name="rmm_store_orig")
    call Info('RigidMeshMapper', 'Storing original coordinates', Level=7)
    call StoreOriginalCoordinates(acc_mesh)
  end subroutine rmm_store_orig

  ! -- mesh velocity --------------------------------------------------
  subroutine rmm_velocity_setup(dim) bind(c, name="rmm_velocity_setup")
    integer(c_int), value :: dim
    type(variable_t), pointer :: velovar
    integer :: i, nn
    nn = SIZE(acc_mesh % nodes % x)
    velovar => VariableGet(acc_mesh % variables, 'Mesh Velocity')
    if (associated(velovar)) then
      acc_veloval => velovar % values
      acc_velperm => velovar % perm
    else
      call Info('RigidMeshMapper', 'Creating variable for mesh velocity')
      allocate(acc_veloval(dim*nn), acc_velperm(nn))
      acc_veloval = 0.0_dp
      do i = 1, nn
        acc_velperm(i) = i
      end do
      if (dim == 2) then
        call VariableAddVector(acc_mesh % variables, acc_mesh, acc_solver, &
            'Mesh Velocity[Mesh Velocity:2]', dim, acc_veloval, acc_velperm)
      else
        call VariableAddVector(acc_mesh % variables, acc_mesh, acc_solver, &
            'Mesh Velocity[Mesh Velocity:3]', dim, acc_veloval, acc_velperm)
      end if
    end if
    call InvalidateVariable(CurrentModel % Meshes, acc_mesh, 'Mesh Velocity')
  end subroutine rmm_velocity_setup

  subroutine rmm_velocity_node(nodei, dx, dt, dim) bind(c, name="rmm_velocity_node")
    integer(c_int), value :: nodei, dim
    real(c_double) :: dx(3)
    real(c_double), value :: dt
    integer :: k
    k = nodei
    if (associated(acc_velperm)) k = acc_velperm(nodei)
    if (k > 0) then
      if (dim == 2) then
        acc_veloval(2*k-1) = (rmm_get_xorig(nodei) + dx(1) - rmm_get_xnew(nodei)) / dt
        acc_veloval(2*k) = (rmm_get_yorig(nodei) + dx(2) - rmm_get_ynew(nodei)) / dt
      else
        acc_veloval(3*k-2) = (rmm_get_xorig(nodei) + dx(1) - rmm_get_xnew(nodei)) / dt
        acc_veloval(3*k-1) = (rmm_get_yorig(nodei) + dx(2) - rmm_get_ynew(nodei)) / dt
        acc_veloval(3*k) = (rmm_get_zorig(nodei) + dx(3) - rmm_get_znew(nodei)) / dt
      end if
    end if
  end subroutine rmm_velocity_node

  ! -- rotation axis order --------------------------------------------
  subroutine rmm_rotorder_fetch(j, found) bind(c, name="rmm_rotorder_fetch")
    integer(c_int) :: j, found
    logical :: f
    acc_intarray => ListGetIntegerArray(acc_solverparams, 'Mesh Rotation Axis Order', f)
    found = merge(1, 0, f)
    if (f) then
      j = SIZE(acc_intarray)
    else
      j = 0
    end if
  end subroutine rmm_rotorder_fetch

  function rmm_rotorder_val(i) result(v) bind(c, name="rmm_rotorder_val")
    integer(c_int), value :: i
    integer(c_int) :: v
    v = acc_intarray(i)
  end function rmm_rotorder_val

  subroutine rmm_set_rotorder(a, b, c) bind(c, name="rmm_set_rotorder")
    integer(c_int), value :: a, b, c
    acc_rotorder(1) = a
    acc_rotorder(2) = b
    acc_rotorder(3) = c
  end subroutine rmm_set_rotorder

  subroutine rmm_fatal_rotaxis() bind(c, name="rmm_fatal_rotaxis")
    call Fatal('RigidMeshMapper', 'Size of Mesh Rotation Axis Order must match dimension.')
  end subroutine rmm_fatal_rotaxis

  ! -- relaxation assembly driver -------------------------------------
  function rmm_assoc_matrix() result(r) bind(c, name="rmm_assoc_matrix")
    integer(c_int) :: r
    r = merge(1, 0, associated(acc_solver % matrix))
  end function rmm_assoc_matrix

  subroutine rmm_relax_bcs() bind(c, name="rmm_relax_bcs")
    type(valuelist_t), pointer :: bc
    integer :: i
    logical :: found
    do i = 1, acc_model % numberofbcs
      bc => acc_model % bcs(i) % values
      if (GetLogical(bc, 'Moving Boundary', found)) then
        call ListAddConstReal(bc, acc_solver % variable % name, 1.0_dp)
      else if (GetLogical(bc, 'Fixed Boundary', found)) then
        call ListAddConstReal(bc, acc_solver % variable % name, 0.0_dp)
      end if
    end do
    do i = 1, acc_model % numberofbodyforces
      bc => acc_model % bodyforces(i) % values
      if (GetLogical(bc, 'Moving Body', found)) then
        call ListAddConstReal(bc, acc_solver % variable % name, 1.0_dp)
      else if (GetLogical(bc, 'Fixed Body', found)) then
        call ListAddConstReal(bc, acc_solver % variable % name, 0.0_dp)
      end if
    end do
  end subroutine rmm_relax_bcs

  function rmm_nonlin_maxiter() result(m) bind(c, name="rmm_nonlin_maxiter")
    integer(c_int) :: m
    logical :: found
    m = GetInteger(acc_solverparams, 'Nonlinear System Max Iterations', found)
    if (.not. found) m = 1
  end function rmm_nonlin_maxiter

  subroutine rmm_nonlin_coeff(coeff, source) bind(c, name="rmm_nonlin_coeff")
    real(c_double) :: coeff, source
    logical :: found
    coeff = GetCReal(acc_solverparams, 'Nonlinear Conductivity Coefficient', found)
    source = GetCReal(acc_solverparams, 'Mesh Relax Source', found)
  end subroutine rmm_nonlin_coeff

  subroutine rmm_default_initialize() bind(c, name="rmm_default_initialize")
    call DefaultInitialize()
  end subroutine rmm_default_initialize

  function rmm_nof_active() result(n) bind(c, name="rmm_nof_active")
    integer(c_int) :: n
    n = GetNOFActive()
  end function rmm_nof_active

  subroutine rmm_begin_active(t, n) bind(c, name="rmm_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes(acc_elem)
  end subroutine rmm_begin_active

  subroutine rmm_default_update(ld, g, f) bind(c, name="rmm_default_update")
    integer(c_int), value :: ld
    real(c_double) :: g(ld, ld), f(ld)
    call DefaultUpdateEquations(g(1:ld, 1:ld), f(1:ld))
  end subroutine rmm_default_update

  subroutine rmm_default_finish_bulk() bind(c, name="rmm_default_finish_bulk")
    call DefaultFinishBulkAssembly()
  end subroutine rmm_default_finish_bulk

  subroutine rmm_default_finish_assembly() bind(c, name="rmm_default_finish_assembly")
    call DefaultFinishAssembly()
  end subroutine rmm_default_finish_assembly

  subroutine rmm_default_dirichlet() bind(c, name="rmm_default_dirichlet")
    call DefaultDirichletBCs()
  end subroutine rmm_default_dirichlet

  function rmm_default_solve() result(norm) bind(c, name="rmm_default_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function rmm_default_solve

  function rmm_nonlin_converged() result(r) bind(c, name="rmm_nonlin_converged")
    integer(c_int) :: r
    r = merge(1, 0, acc_solver % variable % nonlinconverged == 1)
  end function rmm_nonlin_converged

  subroutine rmm_relax_normalize() bind(c, name="rmm_relax_normalize")
    real(kind=dp) :: maxdeform
    maxdeform = MAXVAL(ABS(acc_solver % variable % values))
    maxdeform = ParallelReduction(maxdeform, 2)
    WRITE(Message, '(A,ES12.3)') 'Normalizing deformation by:', maxdeform
    call Info('RigidMeshMapper', Message, Level=6)
    acc_solver % variable % values = acc_solver % variable % values / maxdeform
  end subroutine rmm_relax_normalize

  function rmm_relax_fetch() result(got) bind(c, name="rmm_relax_fetch")
    integer(c_int) :: got
    type(variable_t), pointer :: relaxvar
    logical :: g
    g = .false.
    relaxvar => acc_solver % variable
    if (associated(relaxvar)) then
      if (associated(relaxvar % values)) then
        if (SIZE(relaxvar % values) > 0) then
          g = .true.
          acc_relaxfield => acc_solver % variable % values
          acc_relaxperm => acc_solver % variable % perm
        end if
      end if
    end if
    got = merge(1, 0, g)
  end function rmm_relax_fetch

  subroutine rmm_fatal_alloc() bind(c, name="rmm_fatal_alloc")
    call Fatal('RigidMeshMapper', 'Memory allocation error.')
  end subroutine rmm_fatal_alloc

  ! -- LocalMatrix internal-procedure micros ---------------------------
  subroutine rmm_lm_nodes() bind(c, name="rmm_lm_nodes")
    call GetElementNodes(acc_nodes_lm, acc_elem)
  end subroutine rmm_lm_nodes

  subroutine rmm_lm_localsol(n, arr) bind(c, name="rmm_lm_localsol")
    integer(c_int), value :: n
    real(c_double) :: arr(n)
    call GetScalarLocalSolution(arr(1:n))
  end subroutine rmm_lm_localsol

  subroutine rmm_lm_youngs_begin(n, got, k) bind(c, name="rmm_lm_youngs_begin")
    integer(c_int), value :: n
    integer(c_int) :: got, k
    type(valuelist_t), pointer :: material
    logical :: f
    got = 0
    k = 0
    material => GetMaterial(acc_elem)
    if (associated(material)) then
      call ListGetRealArray(material, 'Mesh Youngs Modulus', &
          acc_hwrk, n, acc_elem % nodeindexes, f)
      if (f) then
        got = 1
        k = SIZE(acc_hwrk, 1)
      end if
    end if
  end subroutine rmm_lm_youngs_begin

  subroutine rmm_lm_elemyoungs(n, k, ey) bind(c, name="rmm_lm_elemyoungs")
    integer(c_int), value :: n, k
    real(c_double) :: ey(3, n)
    integer :: i
    if (k == 1) then
      do i = 1, 3
        ey(i, i:n) = acc_hwrk(1, 1, 1:n)
      end do
    else
      do i = 1, MIN(3, k)
        ey(i, 1:n) = acc_hwrk(i, 1, 1:n)
      end do
      if (k < 3) ey(3, 1:n) = 0.0_dp
    end if
  end subroutine rmm_lm_elemyoungs

  function rmm_lm_ip_begin() result(n) bind(c, name="rmm_lm_ip_begin")
    integer(c_int) :: n
    acc_ip_lm = GaussPoints(acc_elem)
    n = acc_ip_lm % n
  end function rmm_lm_ip_begin

  subroutine rmm_lm_ip_point(t, u, v, w, s) bind(c, name="rmm_lm_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip_lm % u(t)
    v = acc_ip_lm % v(t)
    w = acc_ip_lm % w(t)
    s = acc_ip_lm % s(t)
  end subroutine rmm_lm_ip_point

  subroutine rmm_lm_eleminfo(u, v, w, n, basis, dbasis, detj, stat) &
      bind(c, name="rmm_lm_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: n
    real(c_double) :: basis(n), dbasis(n, 3), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_elem, acc_nodes_lm, u, v, w, detj, &
        basis(1:n), dbasis(1:n, 1:3))
    stat = merge(1, 0, ok)
  end subroutine rmm_lm_eleminfo

  subroutine rmm_lm_grad(n, dbasis, localsol, grad) bind(c, name="rmm_lm_grad")
    integer(c_int), value :: n
    real(c_double) :: dbasis(n, 3), localsol(n), grad(3)
    integer :: i
    do i = 1, 3
      grad(i) = SUM(dbasis(:, i) * localsol(1:n))
    end do
  end subroutine rmm_lm_grad

  subroutine rmm_lm_youngs_ip(n, k, got, basis, ey, youngs) &
      bind(c, name="rmm_lm_youngs_ip")
    integer(c_int), value :: n, k, got
    real(c_double) :: basis(n), ey(3, n), youngs(3)
    integer :: i
    if (got /= 0) then
      do i = 1, 3
        youngs(i) = SUM(basis(1:n) * ey(i, 1:n))
      end do
    else
      youngs(1:3) = 1.0_dp
    end if
  end subroutine rmm_lm_youngs_ip

  subroutine rmm_lm_assemble(n, grad, youngs, coeff, s, detj, source, &
      basis, dbasis, ld, stiff, force) bind(c, name="rmm_lm_assemble")
    integer(c_int), value :: n, ld
    real(c_double) :: grad(3), youngs(3), basis(n), dbasis(n, 3)
    real(c_double), value :: coeff, s, detj, source
    real(c_double) :: stiff(ld, ld), force(ld)
    real(c_double) :: cond, ss
    integer :: i, p, q
    cond = 1.0_dp + coeff * SQRT(SUM(grad * grad))
    ss = s * detj
    do i = 1, 3
      do p = 1, n
        do q = 1, n
          stiff(p, q) = stiff(p, q) + youngs(i) * cond * ss * dbasis(p, i) * dbasis(q, i)
        end do
      end do
    end do
    force(1:n) = force(1:n) + source * ss * basis(1:n)
  end subroutine rmm_lm_assemble

  ! -- obsolete keyword warn -------------------------------------------
  function rmm_obsolete_displacement() result(r) bind(c, name="rmm_obsolete_displacement")
    integer(c_int) :: r
    logical :: any
    any = ListCheckPrefixAnyBodyForce(acc_model, 'Mesh Displacement')
    if (any) then
      call Info('RigidMeshMapper', '> Mesh Displacement < is an obsolete keyword', Level=3)
      call Warn('RigidMeshMapper', 'Replace with > Mesh Translate < ')
      any = .false.
    end if
    r = merge(1, 0, any)
  end function rmm_obsolete_displacement

  ! -- AnyMesh* presence flags ------------------------------------------
  ! mode: 1 means Whole/Rotor mode (ValueList = SolverParams),
  ! otherwise body-force-wide Model checks (as in the original).
  ! which: 1 Matrix, 2 Rotate, 3 Translate, 4 Scale, 5 Origin, 6 Relax.
  function rmm_anymesh(mode, which) result(r) bind(c, name="rmm_anymesh")
    integer(c_int), value :: mode, which
    integer(c_int) :: r
    logical :: v
    if (mode == 1) then
      if (which == 1) then
        v = ListCheckPresent(acc_vlist, 'Mesh Matrix')
      else if (which == 2) then
        v = ListCheckPrefix(acc_vlist, 'Mesh Rotate')
      else if (which == 3) then
        v = ListCheckPrefix(acc_vlist, 'Mesh Translate')
      else if (which == 4) then
        v = ListCheckPrefix(acc_vlist, 'Mesh Scale')
      else if (which == 5) then
        v = ListCheckPrefix(acc_vlist, 'Mesh Origin')
      else
        v = ListCheckPresent(acc_vlist, 'Mesh Relax')
      end if
    else
      if (which == 1) then
        v = ListCheckPresentAnyBodyForce(acc_model, 'Mesh Matrix')
      else if (which == 2) then
        v = ListCheckPrefixAnyBodyForce(acc_model, 'Mesh Rotate')
      else if (which == 3) then
        v = ListCheckPrefixAnyBodyForce(acc_model, 'Mesh Translate')
      else if (which == 4) then
        v = ListCheckPrefixAnyBodyForce(acc_model, 'Mesh Scale')
      else if (which == 5) then
        v = ListCheckPrefixAnyBodyForce(acc_model, 'Mesh Origin')
      else
        v = ListCheckPresentAnyBodyForce(acc_model, 'Mesh Relax')
      end if
    end if
    r = merge(1, 0, v)
  end function rmm_anymesh

  ! -- per-element ValueList selection -----------------------------------
  subroutine rmm_use_solverparams_vlist() bind(c, name="rmm_use_solverparams_vlist")
    acc_vlist => acc_solverparams
  end subroutine rmm_use_solverparams_vlist

  subroutine rmm_vlist_reset() bind(c, name="rmm_vlist_reset")
    acc_prevlist => NULL()
  end subroutine rmm_vlist_reset

  subroutine rmm_elem_bodyforce(bodyid, bf_id, found) bind(c, name="rmm_elem_bodyforce")
    integer(c_int), value :: bodyid
    integer(c_int) :: bf_id, found
    logical :: f
    bf_id = ListGetInteger(acc_model % bodies(bodyid) % values, 'Body Force', f)
    found = merge(1, 0, f)
    if (f) acc_vlist => acc_model % bodyforces(bf_id) % values
  end subroutine rmm_elem_bodyforce

  function rmm_constmap_step() result(got) bind(c, name="rmm_constmap_step")
    integer(c_int) :: got
    got = merge(1, 0, associated(acc_vlist, acc_prevlist))
    acc_prevlist => acc_vlist
  end function rmm_constmap_step

  subroutine rmm_begin_bulk(elem, n) bind(c, name="rmm_begin_bulk")
    integer(c_int), value :: elem
    integer(c_int) :: n
    acc_elem => acc_mesh % elements(elem)
    acc_model % currentelement => acc_elem
    n = GetElementNOFNodes(acc_elem)
  end subroutine rmm_begin_bulk

  function rmm_elem_bodyid() result(b) bind(c, name="rmm_elem_bodyid")
    integer(c_int) :: b
    b = acc_elem % bodyid
  end function rmm_elem_bodyid

  subroutine rmm_node_index(pos, nodei) bind(c, name="rmm_node_index")
    integer(c_int), value :: pos
    integer(c_int) :: nodei
    acc_nidx(1) = acc_elem % nodeindexes(pos)
    nodei = acc_nidx(1)
  end subroutine rmm_node_index

  ! -- mapping keywords --------------------------------------------------
  function rmm_mesh_matrix() result(got) bind(c, name="rmm_mesh_matrix")
    integer(c_int) :: got
    logical :: f
    integer :: i, j
    acc_parray => ListGetConstRealArray(acc_vlist, 'Mesh Matrix', f)
    got = merge(1, 0, f)
    if (f) then
      do i = 1, SIZE(acc_parray, 1)
        do j = 1, SIZE(acc_parray, 2)
          acc_trf(i, j) = acc_parray(j, i)
        end do
      end do
    end if
  end function rmm_mesh_matrix

  ! Fills acc_angles (scaled by coeff); which assumption mirrors the
  ! original: only called when AnyMeshRotate.
  function rmm_mesh_rotate(coeff) result(got) bind(c, name="rmm_mesh_rotate")
    real(c_double), value :: coeff
    integer(c_int) :: got
    logical :: f, found
    real(kind=dp) :: a(1)
    integer :: i
    got = 0
    acc_angles = 0.0_dp
    acc_parray => ListGetConstRealArray(acc_vlist, 'Mesh Rotate', f)
    if (f) then
      got = 1
      do i = 1, SIZE(acc_parray, 1)
        acc_angles(i) = acc_parray(i, 1)
      end do
    else
      a(1:1) = ListGetReal(acc_vlist, 'Mesh Rotate 1', 1, acc_nidx, found)
      if (found) got = 1
      acc_angles(1:1) = a(1:1)
      a(1:1) = ListGetReal(acc_vlist, 'Mesh Rotate 2', 1, acc_nidx, found)
      if (found) got = 1
      acc_angles(2:2) = a(1:1)
      a(1:1) = ListGetReal(acc_vlist, 'Mesh Rotate 3', 1, acc_nidx, found)
      if (found) got = 1
      acc_angles(3:3) = a(1:1)
    end if
    acc_angles = coeff * acc_angles
  end function rmm_mesh_rotate

  subroutine rmm_rotor_override(coeff, rotorangle) bind(c, name="rmm_rotor_override")
    real(c_double), value :: coeff, rotorangle
    acc_angles(3) = coeff * rotorangle
  end subroutine rmm_rotor_override

  function rmm_mesh_scale() result(got) bind(c, name="rmm_mesh_scale")
    integer(c_int) :: got
    logical :: f, found
    real(kind=dp) :: a(1)
    integer :: i
    got = 0
    acc_parray => ListGetConstRealArray(acc_vlist, 'Mesh Scale', f)
    if (f) then
      got = 1
      acc_scaling = 0.0_dp
      do i = 1, SIZE(acc_parray, 1)
        acc_scaling(i) = acc_parray(i, 1)
      end do
    else
      a(1:1) = ListGetReal(acc_vlist, 'Mesh Scale 1', 1, acc_nidx, f)
      got = merge(1, 0, f)
      if (.not. f) a(1) = 1.0_dp
      acc_scaling(1:1) = a(1:1)
      a(1:1) = ListGetReal(acc_vlist, 'Mesh Scale 2', 1, acc_nidx, found)
      if (.not. found) a(1) = 1.0_dp
      if (found) got = 1
      acc_scaling(2:2) = a(1:1)
      a(1:1) = ListGetReal(acc_vlist, 'Mesh Scale 3', 1, acc_nidx, found)
      if (.not. found) a(1) = 1.0_dp
      if (found) got = 1
      acc_scaling(3:3) = a(1:1)
    end if
  end function rmm_mesh_scale

  function rmm_mesh_translate() result(got) bind(c, name="rmm_mesh_translate")
    integer(c_int) :: got
    logical :: f, found
    real(kind=dp) :: a(1)
    integer :: i
    got = 0
    acc_parray => ListGetConstRealArray(acc_vlist, 'Mesh Translate', f)
    if (f) then
      got = 1
      acc_dcoord = 0.0_dp
      do i = 1, SIZE(acc_parray, 1)
        acc_dcoord(i) = acc_parray(i, 1)
      end do
    else
      a(1:1) = ListGetReal(acc_vlist, 'Mesh Translate 1', 1, acc_nidx, f)
      got = merge(1, 0, f)
      acc_dcoord(1:1) = a(1:1)
      a(1:1) = ListGetReal(acc_vlist, 'Mesh Translate 2', 1, acc_nidx, found)
      if (found) got = 1
      acc_dcoord(2:2) = a(1:1)
      a(1:1) = ListGetReal(acc_vlist, 'Mesh Translate 3', 1, acc_nidx, found)
      if (found) got = 1
      acc_dcoord(3:3) = a(1:1)
    end if
  end function rmm_mesh_translate

  ! Verbatim rotation-matrix assembly (TrsMatrix/SclMatrix identity
  ! bases are built in rmm_build_trf, as in the original).
  subroutine rmm_build_rot(gotrotate) bind(c, name="rmm_build_rot")
    integer(c_int), value :: gotrotate
    integer :: i, j, k
    real(kind=dp) :: alpha
    if (gotrotate /= 0) then
      acc_rot = 0.0_dp
      do k = 1, 4
        acc_rot(k, k) = 1.0_dp
      end do
      do i = 1, 3
        j = acc_rotorder(i)
        alpha = acc_angles(j)
        if (ABS(alpha) < TINY(alpha)) cycle
        acc_trf = 0.0_dp
        do k = 1, 4
          acc_trf(k, k) = 1.0_dp
        end do
        select case (j)
        case (1)
          acc_trf(2, 2) = COS(alpha)
          acc_trf(2, 3) = -SIN(alpha)
          acc_trf(3, 2) = SIN(alpha)
          acc_trf(3, 3) = COS(alpha)
        case (2)
          acc_trf(1, 1) = COS(alpha)
          acc_trf(1, 3) = -SIN(alpha)
          acc_trf(3, 1) = SIN(alpha)
          acc_trf(3, 3) = COS(alpha)
        case (3)
          acc_trf(1, 1) = COS(alpha)
          acc_trf(1, 2) = -SIN(alpha)
          acc_trf(2, 1) = SIN(alpha)
          acc_trf(2, 2) = COS(alpha)
        end select
        acc_rot = MATMUL(acc_rot, acc_trf)
      end do
    end if
  end subroutine rmm_build_rot

  subroutine rmm_build_trf(gottranslate, tbr, gotscale) bind(c, name="rmm_build_trf")
    integer(c_int), value :: gottranslate, tbr, gotscale
    real(kind=dp) :: trs(4, 4), scl(4, 4)
    integer :: i, k
    trs = 0.0_dp
    scl = 0.0_dp
    do k = 1, 4
      trs(k, k) = 1.0_dp
      scl(k, k) = 1.0_dp
    end do
    if (gottranslate /= 0) then
      do i = 1, 3
        trs(i, 4) = acc_dcoord(i)
      end do
    end if
    if (tbr /= 0) then
      acc_trf = MATMUL(acc_rot, trs)
    else
      acc_trf = MATMUL(trs, acc_rot)
    end if
    if (gotscale /= 0) then
      do i = 1, 3
        scl(i, i) = acc_scaling(i)
      end do
      trs = acc_trf
      acc_trf = MATMUL(scl, trs)
    end if
  end subroutine rmm_build_trf

  subroutine rmm_mesh_origin(gotmatrix, anymeshorigin) bind(c, name="rmm_mesh_origin")
    integer(c_int), value :: gotmatrix, anymeshorigin
    logical :: found
    real(kind=dp) :: a(1)
    integer :: i
    acc_origin = 0.0_dp
    if ((gotmatrix /= 0) .and. (anymeshorigin /= 0)) then
      acc_parray => ListGetConstRealArray(acc_vlist, 'Mesh Origin', found)
      if (found) then
        do i = 1, SIZE(acc_parray, 1)
          acc_origin(i) = acc_parray(i, 1)
        end do
      else
        a(1:1) = ListGetReal(acc_vlist, 'Mesh Origin 1', 1, acc_nidx, found)
        acc_origin(1:1) = a(1:1)
        a(1:1) = ListGetReal(acc_vlist, 'Mesh Origin 2', 1, acc_nidx, found)
        acc_origin(2:2) = a(1:1)
        a(1:1) = ListGetReal(acc_vlist, 'Mesh Origin 3', 1, acc_nidx, found)
        acc_origin(3:3) = a(1:1)
      end if
    end if
  end subroutine rmm_mesh_origin

  subroutine rmm_apply_node(nodei, gotmatrix, gottranslate, dx, cycled) &
      bind(c, name="rmm_apply_node")
    integer(c_int), value :: nodei, gotmatrix, gottranslate
    real(c_double) :: dx(3)
    integer(c_int) :: cycled
    real(kind=dp) :: x0(4), x1(4)
    if (gotmatrix /= 0) then
      x0(1) = rmm_get_xorig(nodei)
      x0(2) = rmm_get_yorig(nodei)
      x0(3) = rmm_get_zorig(nodei)
      x0(4) = 1.0_dp
      x1 = MATMUL(acc_trf, x0 - acc_origin) + acc_origin
      dx(1:3) = x1(1:3) / x1(4) - x0(1:3)
      cycled = 0
    else if (gottranslate /= 0) then
      dx(1:3) = acc_dcoord(1:3)
      cycled = 0
    else
      cycled = 1
    end if
  end subroutine rmm_apply_node

  subroutine rmm_relax_node(nodei, dx, gotrelax, anyrelax) &
      bind(c, name="rmm_relax_node")
    integer(c_int), value :: nodei, gotrelax, anyrelax
    real(c_double) :: dx(3)
    real(kind=dp) :: relax(1)
    logical :: found
    if (gotrelax /= 0) then
      if (acc_relaxperm(nodei) > 0) then
        relax(1:1) = acc_relaxfield(acc_relaxperm(nodei))
        dx = relax(1) * dx
      end if
    else if (anyrelax /= 0) then
      relax(1:1) = ListGetReal(acc_vlist, 'Mesh Relax', 1, acc_nidx, found)
      if (found) dx = relax(1) * dx
    end if
  end subroutine rmm_relax_node

  subroutine rmm_write_node(nodei, dx) bind(c, name="rmm_write_node")
    integer(c_int), value :: nodei
    real(c_double) :: dx(3)
    acc_mesh % nodes % x(nodei) = dx(1) + rmm_get_xorig(nodei)
    acc_mesh % nodes % y(nodei) = dx(2) + rmm_get_yorig(nodei)
    acc_mesh % nodes % z(nodei) = dx(3) + rmm_get_zorig(nodei)
  end subroutine rmm_write_node

  ! -- info / warn / fatal ---------------------------------------------
  subroutine rmm_info_start() bind(c, name="rmm_info_start")
    call Info('RigidMeshMapper', '---------------------------------------', Level=4)
    call Info('RigidMeshMapper', 'Performing analytic mesh mapping ', Level=4)
    call Info('RigidMeshMapper', '---------------------------------------', Level=4)
  end subroutine rmm_info_start

  subroutine rmm_info_targetmesh(k) bind(c, name="rmm_info_targetmesh")
    integer(c_int), value :: k
    call Info('RigidMeshMapper', 'Operating on mesh from solver ' // I2S(k))
  end subroutine rmm_info_targetmesh

  subroutine rmm_info_movewhole() bind(c, name="rmm_info_movewhole")
    call Info('RigidMeshMapper', 'Moving the whole mesh with keywords from Solver section!')
  end subroutine rmm_info_movewhole

  subroutine rmm_info_solvecaller() bind(c, name="rmm_info_solvecaller")
    call Info('RigidMeshMapper', &
        'Solving mesh relaxation field: ' // TRIM(acc_solver % variable % name), Level=5)
  end subroutine rmm_info_solvecaller

  subroutine rmm_info_laplace() bind(c, name="rmm_info_laplace")
    call Info('RigidMeshMapper', 'Solving mesh relaxation field using Laplace', Level=6)
  end subroutine rmm_info_laplace

  subroutine rmm_info_mapped(count) bind(c, name="rmm_info_mapped")
    integer(c_int), value :: count
    call Info('RigidMeshMapper', 'Number of nodes mapped: ' // I2S(count))
  end subroutine rmm_info_mapped

  function rmm_cputime() result(t) bind(c, name="rmm_cputime")
    real(c_double) :: t
    interface
      function c_cputime() result(x) bind(c, name="cputime")
        import :: c_double
        real(c_double) :: x
      end function c_cputime
    end interface
    t = c_cputime()
  end function rmm_cputime

  subroutine rmm_info_time(at1, at0) bind(c, name="rmm_info_time")
    real(c_double), value :: at1, at0
    if (at1 - at0 > 0.1_dp) then
      WRITE(Message, *) 'Coordinate mapping time: ', at1 - at0
      call Info('RigidMeshMapper', Message)
    end if
  end subroutine rmm_info_time

  subroutine rmm_info_done() bind(c, name="rmm_info_done")
    call Info('RigidMeshMapper', 'All done')
  end subroutine rmm_info_done

  subroutine rmm_fatal_rotor() bind(c, name="rmm_fatal_rotor")
    call Fatal('RigidMeshMapper', 'Rotor state lost, cannot continue!')
  end subroutine rmm_fatal_rotor

  subroutine rmm_default_finish() bind(c, name="rmm_default_finish")
    call DefaultFinish()
  end subroutine rmm_default_finish

end module rigidmeshmapper_acc
