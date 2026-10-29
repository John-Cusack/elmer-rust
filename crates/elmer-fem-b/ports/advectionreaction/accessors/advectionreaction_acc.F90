! T022-A AdvectionReaction port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/AdvectionReaction.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays, C_INT
! codes). All Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER
! dummies, assumed-shape/descriptor dummies, hidden lengths, list-directed
! timing reports — stays inside Fortran, so the Rust side is pure C calls.
! The wrapper module mirrors AdvectionReactionSolver's own SAVE locals
! (assembly arrays, AllocationsDone, DIM, VariableName, the internal
! procedures' SAVE Nodes with identical lifetimes) plus the entry
! current-element context. All floating-point work (the Galerkin assembly,
! SUM reductions, the CPUTime-based reports) and every formatted WRITE
! stay in Fortran (same compiler), so bit-identity holds by construction
! wherever the call sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module advectionreaction_acc
  use iso_c_binding
  use defutils
  use messages
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(mesh_t), pointer, save :: acc_mesh => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(element_t), pointer, save :: acc_face => null()
  type(element_t), pointer, save :: acc_parent => null()
  type(element_t), pointer, save :: acc_left => null()
  type(element_t), pointer, save :: acc_right => null()
  type(element_t), pointer, save :: acc_faces(:) => null()
  type(element_t), pointer, save :: acc_bc_elem => null()
  type(valuelist_t), pointer, save :: acc_material => null()
  type(valuelist_t), pointer, save :: acc_equation => null()
  type(valuelist_t), pointer, save :: acc_bodyforce => null()
  type(valuelist_t), pointer, save :: acc_bc => null()
  type(valuelist_t), pointer, save :: acc_params => null()
  ! Entry SAVE state (mirrors the solver body's SAVE).
  real(kind=dp), allocatable, save :: acc_mass(:,:), acc_stiff(:,:), &
      acc_load(:), acc_force(:), acc_velo(:,:), acc_meshvelo(:,:), &
      acc_gamma(:), acc_upper(:), acc_lower(:)
  integer, save :: acc_n = 0
  logical, save :: acc_alloc_done = .false.
  integer, save :: acc_dim = 0
  character(len=MAX_NAME_LEN), save :: acc_varname = ''
  ! Full-length solver-name scratch (the ALE internal takes
  ! CHARACTER(LEN=MAX_NAME_LEN); a short literal would work but warns).
  character(len=MAX_NAME_LEN), save :: acc_solvername = 'AdvectionReaction'
  ! Internal-procedure SAVE Nodes (mirrors each internal's SAVE).
  type(nodes_t), target, save :: acc_nodes_lm
  type(nodes_t), target, save :: acc_nodes_jump_face
  type(nodes_t), target, save :: acc_nodes_jump_left
  type(nodes_t), target, save :: acc_nodes_jump_right
  type(nodes_t), target, save :: acc_nodes_bnd
  type(nodes_t), target, save :: acc_nodes_bnd_parent
  ! Iteration timing (mirrors the body's at/at0/st/totat/totst locals,
  ! which persist across the nonlinear loop within one invocation).
  real(kind=dp), save :: acc_at = 0.0_dp, acc_at0 = 0.0_dp
  real(kind=dp), save :: acc_st = 0.0_dp
  real(kind=dp), save :: acc_totat = 0.0_dp, acc_totst = 0.0_dp
  ! Limit counters (mirrors the body's CorrectedUpperLimit /
  ! CorrectedLowerLimit locals; reset at the same statement).
  integer, save :: acc_corr_upper = 0, acc_corr_lower = 0

contains

  ! -- entry state ----------------------------------------------------
  subroutine ara_begin(h_solver) bind(c, name="ara_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine ara_begin

  subroutine ara_info_header() bind(c, name="ara_info_header")
    acc_varname = trim(acc_solver % variable % name)
    write(message, '(A,A)') &
        'AdvectionReactionSolver for variable ', acc_varname
    call Info('AdvectionReaction', message, level=5)
  end subroutine ara_info_header

  function ara_faces() result(nfaces) bind(c, name="ara_faces")
    integer(c_int) :: nfaces
    acc_mesh => GetMesh()
    if (CoordinateSystemDimension() == 2) then
      acc_faces => acc_mesh % edges
      nfaces = acc_mesh % numberofedges
    else
      acc_faces => acc_mesh % faces
      nfaces = acc_mesh % numberoffaces
    end if
  end function ara_faces

  function ara_num_boundary_elems() result(n) bind(c, name="ara_num_boundary_elems")
    integer(c_int) :: n
    n = acc_mesh % numberofboundaryelements
  end function ara_num_boundary_elems

  subroutine ara_alloc() bind(c, name="ara_alloc")
    integer :: n, istat
    if (.not. acc_alloc_done) then
      n = 2 * max(acc_mesh % maxelementdofs, acc_mesh % maxelementnodes)
      allocate(acc_force(n), acc_mass(n, n), acc_stiff(n, n), acc_load(n), &
          acc_velo(3, n), acc_meshvelo(3, n), acc_gamma(n), &
          acc_upper(n), acc_lower(n), stat=istat)
      if (istat /= 0) then
        call Fatal('AdvectionReaction', 'Memory allocation error.')
      else
        call Info('AdvectionReaction', 'Memory allocation done', level=12)
      end if
      acc_n = n
      acc_dim = CoordinateSystemDimension()
      acc_alloc_done = .true.
    end if
  end subroutine ara_alloc

  subroutine ara_read_params(nmin, nmax, limit) bind(c, name="ara_read_params")
    integer(c_int) :: nmin, nmax, limit
    logical :: found, lim
    integer :: mn, mx
    acc_params => GetSolverParams()
    mx = GetInteger(acc_params, 'Nonlinear System Max Iterations', found)
    if (.not. found) then
      call Warn('AdvectionReaction', &
          'No > Nonlinear System Max Iterations < found. Setting 1')
      mx = 1
    end if
    mn = GetInteger(acc_params, 'Nonlinear System Min Iterations', found)
    if (.not. found) then
      call Warn('AdvectionReaction', &
          'No >Nonlinear System Min Iterations< found. Setting 1')
      mn = 1
    else if (mn > mx) then
      call Warn('AdvectionReaction', &
          '>Nonlinear System Min Iterations< is exceeding ' // &
          '>Nonlinear System Max Iterations<.')
      call Warn('AdvectionReaction', &
          'First is being reset to the latter.')
      mn = mx
    end if
    lim = GetLogical(acc_params, 'Limit Solution', found)
    if (.not. found) lim = .false.
    if (lim) then
      call Info('AdvectionReaction', &
          'Keyword > Limit Solution < found. Solution will be limited', level=6)
    else
      call Info('AdvectionReaction', &
          'No keyword > Limit Solution < found. Solution will not be limited', &
          level=6)
    end if
    nmin = mn
    nmax = mx
    limit = merge(1, 0, lim)
  end subroutine ara_read_params

  subroutine ara_iter_start(iter, maxiter) bind(c, name="ara_iter_start")
    integer(c_int), value :: iter, maxiter
    ! totat/totst reset matches the original (`totat = 0; totst = 0`
    ! executes on every entry before the loop; the Rust loop calls this
    ! once per nonlinear iteration, so reset only on the first).
    if (iter == 1) then
      acc_totat = 0.0_dp
      acc_totst = 0.0_dp
    end if
    acc_at = CPUTime()
    acc_at0 = RealTime()
    call Info('AdvectionReaction', ' ', level=4)
    call Info('AdvectionReaction', '-------------------------------------', level=4)
    write(message, '(A,I0,A,I0)') &
        'Nonlinear iteration no.', iter, ' of max', maxiter
    call Info('AdvectionReaction', message, level=4)
    call Info('AdvectionReaction', '-------------------------------------', level=4)
    call Info('AdvectionReaction', ' ', level=6)
    call Info('AdvectionReaction', 'Starting Assembly...', level=6)
  end subroutine ara_iter_start

  subroutine ara_default_initialize() bind(c, name="ara_default_initialize")
    call DefaultInitialize()
  end subroutine ara_default_initialize

  function ara_nof_active() result(n) bind(c, name="ara_nof_active")
    integer(c_int) :: n
    n = GetNOFActive()
  end function ara_nof_active

  ! -- bulk element body (exact transcription of the active loop body) --
  subroutine ara_bulk_assemble(t, active, dolimit) bind(c, name="ara_bulk_assemble")
    integer(c_int), value :: t, active, dolimit
    integer :: n, i, dummyint, indexes(128)
    logical :: found, foundupper, foundlower
    real(kind=dp) :: orig
    if (RealTime() - acc_at0 > 1.0) then
      write(message, '(a,i3,a)') '   Assembly: ', &
          int(100.0 - 100.0 * (active-t) / (1.0 * active)), &
          ' % done'
      call Info('AdvectionReaction', message, level=5)
      acc_at0 = RealTime()
    end if
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes(acc_elem)
    acc_material => GetMaterial(acc_elem)
    acc_bodyforce => GetBodyForce(acc_elem)
    acc_equation => GetEquation()
    dummyint = GetElementDOFs(indexes)
    if (dolimit /= 0) then
      acc_upper(1:n) = GetReal(acc_material, &
          trim(acc_varname) // ' Upper Limit', foundupper)
      acc_lower(1:n) = GetReal(acc_material, &
          trim(acc_varname) // ' Lower Limit', foundlower)
      do i = 1, n
        if (foundupper) then
          orig = acc_solver % variable % values( &
              acc_solver % variable % perm(indexes(i)))
          if (orig > acc_upper(i)) then
            acc_solver % variable % values( &
                acc_solver % variable % perm(indexes(i))) = acc_upper(i)
            acc_corr_upper = acc_corr_upper + 1
          end if
          if (orig < acc_lower(i)) then
            acc_solver % variable % values( &
                acc_solver % variable % perm(indexes(i))) = acc_lower(i)
            acc_corr_lower = acc_corr_lower + 1
          end if
        end if
      end do
    end if
    acc_load(1:n) = GetReal(acc_bodyforce, &
        trim(acc_varname) // ' Source', found)
    call GetLocalALEVelocity(acc_velo, acc_meshvelo, acc_solvername, &
        acc_material, acc_equation, acc_solver, acc_elem)
    acc_gamma(1:n) = GetReal(acc_material, &
        trim(acc_varname) // ' Gamma', found)
    call LocalMatrix(acc_mass, acc_stiff, acc_force, acc_load, &
        acc_velo, acc_meshvelo, acc_gamma, acc_elem, n)
  end subroutine ara_bulk_assemble


  subroutine ara_first_order_time() bind(c, name="ara_first_order_time")
    call Default1stOrderTime(acc_mass, acc_stiff, acc_force)
  end subroutine ara_first_order_time

  subroutine ara_default_update() bind(c, name="ara_default_update")
    call DefaultUpdateEquations(acc_stiff, acc_force)
  end subroutine ara_default_update

  ! -- face loop (interior DG jump terms) ------------------------------
  subroutine ara_face_force_zero() bind(c, name="ara_face_force_zero")
    acc_force = 0.0d0
  end subroutine ara_face_force_zero

  subroutine ara_face_assemble(t) bind(c, name="ara_face_assemble")
    integer(c_int), value :: t
    integer :: n, n1, n2
    acc_face => acc_faces(t)
    if (.not. ActiveBoundaryElement(acc_face)) return
    acc_left => acc_face % boundaryinfo % left
    acc_right => acc_face % boundaryinfo % right
    if (associated(acc_right) .and. associated(acc_left)) then
      n = GetElementNOFNodes(acc_face)
      n1 = GetElementNOFNodes(acc_left)
      n2 = GetElementNOFNodes(acc_right)
      acc_material => GetMaterial(acc_left)
      acc_equation => GetEquation(acc_left)
      call GetLocalALEVelocity(acc_velo, acc_meshvelo, acc_solvername, &
          acc_material, acc_equation, acc_solver, acc_face)
      call LocalJumps(acc_stiff, acc_face, n, acc_left, n1, &
          acc_right, n2, acc_velo)
      call DefaultUpdateEquations(acc_stiff, acc_force, acc_face)
    end if
  end subroutine ara_face_assemble

  subroutine ara_finish_bulk() bind(c, name="ara_finish_bulk")
    call DefaultFinishBulkAssembly()
  end subroutine ara_finish_bulk

  ! -- boundary element body (returns 0 on CYCLE, 1 when assembled) -----
  function ara_boundary_assemble(t) result(code) bind(c, name="ara_boundary_assemble")
    integer(c_int), value :: t
    integer(c_int) :: code
    integer :: n, n1
    logical :: found
    acc_bc_elem => GetBoundaryElement(t)
    if (.not. ActiveBoundaryElement()) then
      code = 0
      return
    end if
    acc_parent => acc_bc_elem % boundaryinfo % left
    if (.not. associated(acc_parent)) &
        acc_parent => acc_bc_elem % boundaryinfo % right
    n = GetElementNOFNodes(acc_bc_elem)
    n1 = GetElementNOFNodes(acc_parent)
    acc_material => GetMaterial(acc_parent)
    acc_equation => GetEquation(acc_parent)
    call GetLocalALEVelocity(acc_velo, acc_meshvelo, acc_solvername, &
        acc_material, acc_equation, acc_solver, acc_bc_elem)
    acc_bc => GetBC()
    acc_load = 0.0d0
    found = .false.
    if (associated(acc_bc)) then
      acc_load(1:n) = GetReal(acc_bc, acc_solver % variable % name, found)
    end if
    acc_mass = 0.0d0
    call LocalMatrixBoundary(acc_stiff, acc_force, acc_load, &
        acc_bc_elem, n, acc_parent, n1, acc_velo, found)
    code = 1
  end function ara_boundary_assemble

  subroutine ara_finish_boundary() bind(c, name="ara_finish_boundary")
    call DefaultFinishBoundaryAssembly()
  end subroutine ara_finish_boundary

  subroutine ara_finish_assembly() bind(c, name="ara_finish_assembly")
    call DefaultFinishAssembly()
  end subroutine ara_finish_assembly

  ! -- solve + convergence tail (returns 1 on converged EXIT) -----------
  function ara_solve_tail(iter, maxiter) result(code) bind(c, name="ara_solve_tail")
    integer(c_int), value :: iter, maxiter
    integer(c_int) :: code
    real(kind=dp) :: norm, relchange
    acc_at = CPUTime() - acc_at
    acc_st = CPUTime()
    norm = DefaultSolve()
    acc_st = CPUTime() - acc_st
    acc_totat = acc_totat + acc_at
    acc_totst = acc_totst + acc_st
    write(message, '(a,i4,a,F8.2,F8.2)') &
        'iter: ', iter, ' Assembly: (s)', acc_at, acc_totat
    call Info('AdvectionReaction', message, level=4)
    write(message, '(a,i4,a,F8.2,F8.2)') &
        'iter: ', iter, ' Solve:    (s)', acc_st, acc_totst
    call Info('AdvectionReaction', message, level=4)
    relchange = acc_solver % variable % nonlinchange
    if (acc_solver % variable % nonlinconverged == 1) then
      write(message, '(A,I0,A,I0,A)') &
          'Nonlinear iteration converged after ', iter, &
          ' out of max ', maxiter, ' iterations'
      call Info('AdvectionReaction', message)
      code = 1
    else if ((iter == maxiter) .and. (maxiter > 1)) then
      call Warn('AdvectionReaction', &
          'Maximum nonlinear iterations reached, but system not converged')
      code = 0
    else
      code = 0
    end if
  end function ara_solve_tail

  ! -- post-nonlinear-loop limit section --------------------------------
  subroutine ara_limit_pass(do_limit) bind(c, name="ara_limit_pass")
    integer(c_int), value :: do_limit
    integer :: n, i, t, active, dummyint, indexes(128)
    logical :: foundupper, foundlower
    real(kind=dp) :: orig
    acc_corr_upper = 0
    acc_corr_lower = 0
    if (do_limit /= 0) then
      active = GetNOFActive()
      do t = 1, active
        acc_elem => GetActiveElement(t)
        n = GetElementNOFNodes(acc_elem)
        acc_material => GetMaterial()
        dummyint = GetElementDOFs(indexes)
        acc_upper(1:n) = GetReal(acc_material, &
            trim(acc_varname) // ' Upper Limit', foundupper)
        acc_lower(1:n) = GetReal(acc_material, &
            trim(acc_varname) // ' Lower Limit', foundlower)
        do i = 1, n
          if (foundupper) then
            orig = acc_solver % variable % values( &
                acc_solver % variable % perm(indexes(i)))
            if (orig > acc_upper(i)) then
              acc_solver % variable % values( &
                  acc_solver % variable % perm(indexes(i))) = &
                  min(acc_upper(i), orig)
              acc_corr_upper = acc_corr_upper + 1
            end if
            if (orig < acc_lower(i)) then
              acc_solver % variable % values( &
                  acc_solver % variable % perm(indexes(i))) = &
                  max(acc_lower(i), orig)
              acc_corr_lower = acc_corr_lower + 1
            end if
          end if
        end do
      end do
      write(message, '(a,i10)') &
          'Limited values for upper limit: ', acc_corr_upper
      call Info('AdvectionReaction', message, level=4)
      write(message, '(a,i10)') &
          'Limited values for lower limit: ', acc_corr_lower
      call Info('AdvectionReaction', message, level=4)
    end if
  end subroutine ara_limit_pass

  ! -- AdvectionReactionSolver_init body ---------------------------------
  subroutine ara_init_begin(h_solver) bind(c, name="ara_init_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine ara_init_begin

  subroutine ara_init_body() bind(c, name="ara_init_body")
    type(valuelist_t), pointer :: solverparams
    logical :: found
    character(len=MAX_NAME_LEN) :: variablename
    solverparams => GetSolverParams()
    call ListAddLogical(solverparams, 'Discontinuous Galerkin', .true.)
    call ListAddNewString(solverparams, 'Variable', 'Tracer')
    variablename = ListGetString(solverparams, 'Variable')
    if (ListGetLogical(solverparams, 'Calculate Nodal Average', found)) then
      call ListAddString(solverparams, &
          NextFreeKeyword('Exported Variable ', solverparams), &
          '-nodal ' // trim(variablename) // '_nodal')
      call ListAddLogical(solverparams, 'PostSolver Active', .true.)
    end if
  end subroutine ara_init_body

  ! -- AdvectionReactionSolver_post body ---------------------------------
  subroutine ara_post_begin(h_solver) bind(c, name="ara_post_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine ara_post_begin

  function ara_post_body() result(code) bind(c, name="ara_post_body")
    integer(c_int) :: code
    type(mesh_t), pointer :: mesh
    type(variable_t), pointer :: var
    character(len=MAX_NAME_LEN) :: variablename
    integer, allocatable :: cnt(:)
    integer :: i, j, k, n, t, active
    mesh => GetMesh()
    variablename = trim(acc_solver % variable % name) // '_nodal'
    var => VariableGet(mesh % variables, trim(variablename))
    if (.not. associated(var)) then
      variablename = 'Nodal Result'
      var => VariableGet(mesh % variables, trim(variablename))
    end if
    if (.not. associated(var)) then
      code = 0
      return
    end if
    call Info('AdvectionReaction_post', &
        'Using "' // trim(variablename) // '" for the nodal average output field', &
        level=7)
    allocate(cnt(size(var % values)))
    cnt = 0
    var % values = 0.0_dp
    active = GetNOFActive()
    do t = 1, active
      acc_elem => GetActiveElement(t)
      n = GetElementNOFNodes()
      do i = 1, n
        j = acc_elem % dgindexes(i)
        k = acc_elem % nodeindexes(i)
        j = acc_solver % variable % perm(j)
        if (j == 0) cycle
        k = var % perm(k)
        if (k == 0) cycle
        var % values(k) = var % values(k) + acc_solver % variable % values(j)
        cnt(k) = cnt(k) + 1
      end do
    end do
    where (cnt > 0)
      var % values = var % values / cnt
    end where
    deallocate(cnt)
    code = 1
  end function ara_post_body

  ! -- CONTAINS internals (exact transcriptions, module procedures) -----
  subroutine LocalMatrix(MASS, STIFF, FORCE, LOAD, Velo, Mvelo, Gamma, &
      Element, n)
    real(kind=dp) :: MASS(:,:), STIFF(:,:), FORCE(:), &
        LOAD(:), Velo(:,:), Gamma(:), Mvelo(:,:)
    integer :: n
    type(element_t), target :: Element
    real(kind=dp) :: Basis(n), dBasisdx(n, 3)
    real(kind=dp) :: detJ, U, V, W, S, A, L, cu(3), g, divMVelo
    logical :: Stat
    integer :: i, p, q, t, dim
    type(GaussIntegrationPoints_t) :: IntegStuff
    type(Nodes_t) :: Nodes
    dim = CoordinateSystemDimension()
    FORCE = 0.0d0
    STIFF = 0.0d0
    MASS = 0.0d0
    call GetElementNodes(Nodes, Element)
    IntegStuff = GaussPoints(Element)
    do t = 1, IntegStuff % n
      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)
      S = IntegStuff % s(t)
      stat = ElementInfo(Element, Nodes, U, V, W, detJ, Basis, dBasisdx)
      S = S * detJ
      L = sum(LOAD(1:n) * Basis(1:n))
      g = sum(Basis(1:n) * Gamma(1:n))
      divMVelo = 0.0_dp
      do i = 1, dim
        divMVelo = divMVelo + sum(Mvelo(i, 1:n) * dBasisdx(1:n, i))
      end do
      cu = 0.0d0
      do i = 1, dim
        cu(i) = sum(Basis(1:n) * Velo(i, 1:n))
      end do
      do p = 1, n
        do q = 1, n
          MASS(p, q) = MASS(p, q) + s * Basis(q) * Basis(p)
          STIFF(p, q) = STIFF(p, q) + s * (g + divMVelo) * Basis(q) * Basis(p)
          do i = 1, dim
            STIFF(p, q) = STIFF(p, q) - s * cu(i) * Basis(q) * dBasisdx(p, i)
          end do
        end do
      end do
      FORCE(1:n) = FORCE(1:n) + s * L * Basis(1:n)
    end do
  end subroutine LocalMatrix

  subroutine LocalJumps(STIFF, Face, n, LeftParent, n1, RightParent, n2, Velo)
    real(kind=dp) :: STIFF(:,:), Velo(:,:)
    integer :: n, n1, n2
    type(element_t), target :: Face, LeftParent, RightParent
    real(kind=dp) :: FaceBasis(n), FacedBasisdx(n, 3)
    real(kind=dp) :: LeftBasis(n1), LeftdBasisdx(n1, 3)
    real(kind=dp) :: RightBasis(n2), RightdBasisdx(n2, 3)
    real(kind=dp) :: Jump(n1 + n2), Average(n1 + n2)
    real(kind=dp) :: detJ, U, V, W, S, Udotn, xx, yy
    logical :: Stat
    integer :: i, j, p, q, dim, t, nFace, nParent
    type(GaussIntegrationPoints_t) :: IntegStuff
    real(kind=dp) :: hE, Normal(3), cu(3), LeftOut(3)
    dim = CoordinateSystemDimension()
    STIFF = 0.0d0
    call GetElementNodes(acc_nodes_jump_face, Face)
    call GetElementNodes(acc_nodes_jump_left, LeftParent)
    call GetElementNodes(acc_nodes_jump_right, RightParent)
    IntegStuff = GaussPoints(Face)
    LeftOut(1) = sum(acc_nodes_jump_left % x(1:n1)) / n1
    LeftOut(2) = sum(acc_nodes_jump_left % y(1:n1)) / n1
    LeftOut(3) = sum(acc_nodes_jump_left % z(1:n1)) / n1
    LeftOut(1) = sum(acc_nodes_jump_face % x(1:n)) / n - LeftOut(1)
    LeftOut(2) = sum(acc_nodes_jump_face % y(1:n)) / n - LeftOut(2)
    LeftOut(3) = sum(acc_nodes_jump_face % z(1:n)) / n - LeftOut(3)
    do t = 1, IntegStuff % n
      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)
      S = IntegStuff % s(t)
      stat = ElementInfo(Face, acc_nodes_jump_face, U, V, W, detJ, &
          FaceBasis, FacedBasisdx)
      S = S * detJ
      Normal = NormalVector(Face, acc_nodes_jump_face, U, V, .false.)
      if (sum(LeftOut * Normal) < 0) Normal = -Normal
      call FindParentUVW(Face, n, LeftParent, n1, U, V, W, FaceBasis)
      stat = ElementInfo(LeftParent, acc_nodes_jump_left, U, V, W, detJ, &
          LeftBasis, LeftdBasisdx)
      call FindParentUVW(Face, n, RightParent, n2, U, V, W, FaceBasis)
      stat = ElementInfo(RightParent, acc_nodes_jump_right, U, V, W, &
          detJ, RightBasis, RightdBasisdx)
      Jump(1:n1) = LeftBasis(1:n1)
      Jump(n1 + 1:n1 + n2) = -RightBasis(1:n2)
      Average(1:n1) = LeftBasis(1:n1) / 2
      Average(n1 + 1:n1 + n2) = RightBasis(1:n2) / 2
      cu = 0.0d0
      do i = 1, dim
        cu(i) = sum(Velo(i, 1:n) * FaceBasis(1:n))
      end do
      Udotn = sum(Normal * cu)
      do p = 1, n1 + n2
        do q = 1, n1 + n2
          STIFF(p, q) = STIFF(p, q) + s * Udotn * Average(q) * Jump(p)
          STIFF(p, q) = STIFF(p, q) + s * abs(Udotn) / 2 * Jump(q) * Jump(p)
        end do
      end do
    end do
  end subroutine LocalJumps

  subroutine LocalMatrixBoundary(STIFF, FORCE, LOAD, &
      Element, n, ParentElement, np, Velo, InFlowBC)
    real(kind=dp) :: STIFF(:,:), FORCE(:), LOAD(:), Velo(:,:)
    integer :: n, np
    logical :: InFlowBC
    type(element_t), target :: Element, ParentElement
    real(kind=dp) :: Basis(n), dBasisdx(n, 3)
    real(kind=dp) :: ParentBasis(np), ParentdBasisdx(np, 3)
    integer :: i, j, p, q, t, dim
    real(kind=dp) :: Normal(3), g, L, Udotn, cu(3), cu1(3), detJ, U, V, W, S
    logical :: Stat, Inflow
    type(GaussIntegrationPoints_t) :: IntegStuff
    dim = CoordinateSystemDimension()
    FORCE = 0.0d0
    STIFF = 0.0d0
    call GetElementNodes(acc_nodes_bnd, Element)
    call GetElementNodes(acc_nodes_bnd_parent, ParentElement)
    Normal = NormalVector(Element, acc_nodes_bnd, 0.0d0, 0.0d0, .true.)
    do i = 1, 3
      cu(i) = sum(Velo(i, 1:n)) / n
    end do
    Inflow = InFlowBC .and. sum(Normal * cu) < 0.0d0
    IntegStuff = GaussPoints(Element)
    do t = 1, IntegStuff % n
      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)
      S = IntegStuff % s(t)
      Normal = NormalVector(Element, acc_nodes_bnd, U, V, .true.)
      stat = ElementInfo(Element, acc_nodes_bnd, U, V, W, detJ, &
          Basis, dBasisdx)
      S = S * detJ
      call FindParentUVW(Element, n, ParentElement, np, U, V, W, Basis)
      stat = ElementInfo(ParentElement, acc_nodes_bnd_parent, U, V, W, &
          detJ, ParentBasis, ParentdBasisdx)
      L = sum(LOAD(1:n) * Basis(1:n))
      cu = 0.0d0
      do i = 1, dim
        cu(i) = sum(Velo(i, 1:n) * Basis(1:n))
      end do
      Udotn = sum(Normal * cu)
      do p = 1, np
        if (Inflow) then
          FORCE(p) = FORCE(p) - s * Udotn * L * ParentBasis(p)
        else
          do q = 1, np
            STIFF(p, q) = STIFF(p, q) + s * Udotn * ParentBasis(q) * ParentBasis(p)
          end do
        end if
      end do
    end do
  end subroutine LocalMatrixBoundary

  subroutine GetLocalALEVelocity(Velo, MeshVelo, SolverName, Material, &
      Equation, Solver, ModelElement)
    real(kind=dp) :: Velo(:,:), MeshVelo(:,:)
    character(len=MAX_NAME_LEN) :: SolverName
    type(valuelist_t), pointer :: Material, Equation
    type(solver_t), pointer :: Solver
    type(element_t), target :: ModelElement
    character(len=MAX_NAME_LEN) :: ConvectionFlag, FlowSolName
    integer :: i, j, k, N, FlowDOFs
    integer, pointer :: FlowPerm(:)
    real(kind=dp), pointer :: FlowSolution(:)
    type(variable_t), pointer :: FlowSol
    logical :: Found
    N = GetElementNOFNodes(ModelElement)
    ConvectionFlag = GetString(Equation, 'Convection', Found)
    if (.not. Found) &
        call Fatal(SolverName, 'No string for keyword >Convection< found in Equation')
    Velo = 0.0d00
    FlowDOFs = 0
    if (ConvectionFlag == 'constant') then
      Velo(1, 1:acc_n) = GetReal(Material, 'Convection Velocity 1', Found, ModelElement)
      Velo(2, 1:acc_n) = GetReal(Material, 'Convection Velocity 2', Found, ModelElement)
      Velo(3, 1:acc_n) = GetReal(Material, 'Convection Velocity 3', Found, ModelElement)
    else if (ConvectionFlag == 'computed') then
      FlowSolName = GetString(Equation, 'Flow Solution Name', Found)
      if (.not. Found) then
        call Warn(SolverName, 'Keyword >Flow Solution Name< not found in section >Equation<')
        call Warn(SolverName, 'Taking default value >Flow Solution<')
        write(FlowSolName, '(A)') 'Flow Solution'
      end if
      FlowSol => VariableGet(Solver % Mesh % Variables, FlowSolName)
      if (associated(FlowSol)) then
        FlowPerm => FlowSol % Perm
        FlowDOFs = FlowSol % DOFs
        FlowSolution => FlowSol % Values
      else
        write(Message, '(A,A,A)') &
            'Convection flag set to >computed<, but no variable >', FlowSolName, '< found'
        call Fatal(SolverName, Message)
      end if
      do i = 1, n
        k = FlowPerm(ModelElement % NodeIndexes(i))
        if (k > 0) then
          select case (FlowDOFs)
          case (2)
            Velo(1, i) = FlowSolution(FlowDOFs * k - 1)
            Velo(2, i) = 0.0d0
            Velo(3, i) = 0.0d0
          case (3)
            Velo(1, i) = FlowSolution(FlowDOFs * k - 2)
            Velo(2, i) = FlowSolution(FlowDOFs * k - 1)
            Velo(3, i) = 0.0d0
          case (4)
            Velo(1, i) = FlowSolution(FlowDOFs * k - 3)
            Velo(2, i) = FlowSolution(FlowDOFs * k - 2)
            Velo(3, i) = FlowSolution(FlowDOFs * k - 1)
          end select
        end if
      end do
    else if (ConvectionFlag == 'none') then
      Velo = 0.0d0
    else
      write(Message, '(A,A,A)') 'Convection flag >', ConvectionFlag, '< not recognised'
      call Fatal(SolverName, Message)
    end if
    MeshVelo = 0.0d0
    call GetVectorLocalSolution(MeshVelo, 'Mesh Velocity', ModelElement)
    if (any(MeshVelo /= 0.0d0)) then
      do i = 1, FlowDOFs
        Velo(i, 1:acc_n) = Velo(i, 1:acc_n) - MeshVelo(i, 1:acc_n)
      end do
    end if
  end subroutine GetLocalALEVelocity

end module advectionreaction_acc
