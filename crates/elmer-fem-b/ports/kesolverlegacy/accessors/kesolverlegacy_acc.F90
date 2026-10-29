! T022-A KESolverLegacy port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/TurbulenceModels/KESolverLegacy.F90,
! compiled against the built .mod files of the target configure (A2
! shipped mechanism, ADR-057 section 3).
!
! Every wrapper takes only C-interoperable arguments. All Fortran-object
! traffic — TYPE(...) handles, OPTIONAL/POINTER dummies, hidden lengths,
! list-directed timing reports — stays inside Fortran, so the Rust side is
! pure C calls. The wrapper module mirrors KESolverLegacy's own locals
! (entry SAVE arrays with identical first-mesh sizing and lifetimes, the
! per-element bubble history bx/bxprev, the EpsilonWall SAVE Nodes) and
! all floating-point work (the LocalMatrix Galerkin assembly with its
! SQRT/SUM reductions, the RNG secant loop, the K/E clip floors) plus
! all formatted I/O (WRITE(Message,*)/Info, the RNG PRINT*) stays in
! Fortran (same compiler), so bit-identity holds by construction
! wherever the call sequence matches.
!
! Faithfulness notes (mirrors KESolverLegacy.F90 statement for
! statement):
! - Dead stores replicated as documented skips, never executed: `at`/`at0`
!   (CPUTime/RealTime assigned, never read), `NewtonLinearization` (set,
!   never read — like DivergenceSolver's `Visited`), `ForceVector`,
!   `KEpsilon`, `ElectricCurrent`, `S`, `C`, `k1`, `k2`, `eq_id`, `bf_id`.
! - `body_id` is compared but never updated in the original (stays -1),
!   so the material block runs for every element; transcribed exactly.
! - The entry-level `C0(DOFs,N)` array is allocated (same STAT/Fatal
!   path) but never read — LocalMatrix has its own local `C0(2)`;
!   replicated exactly.
! - `KEWall` stays an implicit external exactly as in the original (no
!   USE — it resolves from libelmersolver at load, `U kewall_` in both
!   the reference and the port).
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module kesolverlegacy_acc
  use iso_c_binding
  use defutils
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  real(kind=dp), save :: acc_dt = 0.0_dp
  logical, save :: acc_transient = .false.

  ! Entry scalar state (mirrors the solver body locals).
  integer, save :: acc_dofs = 0, acc_nsdofs = 0, acc_localnodes = 0
  logical, save :: acc_bubbles = .false.
  integer, save :: acc_newtoniter = 0, acc_nonlineariter = 1
  integer, save :: acc_body_id = -1
  integer, save :: acc_bxstride = 0, acc_bubbletimestep = -1
  real(kind=dp), save :: acc_newtontol = 0.0_dp
  real(kind=dp), save :: acc_clip = 1.0d-6, acc_v2fcp = 0.05_dp
  real(kind=dp), save :: acc_norm = 0.0_dp, acc_prevnorm = 0.0_dp
  real(kind=dp), save :: acc_relchange = 0.0_dp
  logical, save :: acc_bubblesdefault = .true.
  logical, save :: acc_allocdone = .false.
  character(len=MAX_NAME_LEN), save :: acc_kemodel = '', acc_v2fmodel = ''

  ! Entry SAVE arrays (first-mesh sizing, never resized — like the
  ! original, which would overflow silently on a bigger later mesh).
  real(kind=dp), allocatable, save :: acc_u(:), acc_v(:), acc_w(:)
  real(kind=dp), allocatable, save :: acc_density(:), acc_viscosity(:)
  real(kind=dp), allocatable, save :: acc_effectivevisc(:,:)
  real(kind=dp), allocatable, save :: acc_turbvisc(:)
  real(kind=dp), allocatable, save :: acc_c0(:,:)
  real(kind=dp), allocatable, save :: acc_layer(:), acc_rough(:)
  real(kind=dp), allocatable, save :: acc_kec1(:), acc_kec2(:)
  real(kind=dp), allocatable, save :: acc_kesigk(:), acc_kesige(:), acc_kecmu(:)
  real(kind=dp), allocatable, save :: acc_k(:), acc_e(:), acc_v2(:), acc_v2fct(:)
  real(kind=dp), allocatable, save :: acc_nodaldensity(:), acc_nodalvisc(:), acc_nodalcmu(:)
  real(kind=dp), allocatable, save :: acc_prevk(:), acc_preve(:)
  real(kind=dp), allocatable, save :: acc_xl(:), acc_xlprev(:)
  real(kind=dp), allocatable, save :: acc_mass(:,:), acc_stiff(:,:)
  real(kind=dp), allocatable, save :: acc_load(:,:), acc_force(:), acc_timeforce(:)
  real(kind=dp), allocatable, save :: acc_bx(:), acc_bxprev(:)
  real(kind=dp), allocatable, save :: acc_work(:)

  ! Element-loop pointers (re-acquired per call, like the original's
  ! per-iteration pointer assignments).
  type(element_t), pointer, save :: acc_elem => null()
  type(valuelist_t), pointer, save :: acc_bc => null()
  type(valuelist_t), pointer, save :: acc_material => null()
  type(valuelist_t), pointer, save :: acc_equation => null()
  type(nodes_t), save :: acc_enodes
  ! EpsilonWall SAVE Nodes/ParentNodes (identical lifetimes).
  type(nodes_t), save :: acc_wnodes, acc_wpnodes

contains

  subroutine ke_begin(h_model, h_solver, dt, transient) bind(c, name="ke_begin")
    type(c_ptr), value :: h_model, h_solver
    real(c_double), value :: dt
    integer(c_int), value :: transient
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    acc_dt = dt
    acc_transient = (transient /= 0)
  end subroutine ke_begin

  ! Guards: no-matrix and no-local-nodes silent returns. Flow Solution
  ! is only probed (no return). Returns 0 to proceed, 1 to return.
  function ke_prologue() result(code) bind(c, name="ke_prologue")
    integer(c_int) :: code
    type(variable_t), pointer :: ke, flowsol
    if (.not. associated(acc_solver % matrix)) then
      code = 1
      return
    end if
    ke => acc_solver % variable
    if (associated(ke)) then
      acc_dofs = ke % dofs
      acc_norm = ke % norm
    end if
    acc_localnodes = count(ke % perm > 0)
    if (acc_localnodes <= 0) then
      code = 1
      return
    end if
    flowsol => VariableGet(acc_model % variables, 'Flow Solution')
    if (associated(flowsol)) then
      acc_nsdofs = flowsol % dofs
    end if
    code = 0
  end function ke_prologue

  subroutine ke_params() bind(c, name="ke_params")
    logical :: gotit
    acc_newtontol = ListGetConstReal(acc_solver % values, &
        'Nonlinear System Newton After Tolerance', gotit)
    acc_newtoniter = ListGetInteger(acc_solver % values, &
        'Nonlinear System Newton After Iterations', gotit)
    acc_nonlineariter = ListGetInteger(acc_solver % values, &
        'Nonlinear System Max Iterations', gotit)
    if (.not. gotit) acc_nonlineariter = 1
    acc_bubblesdefault = ListGetLogical(acc_solver % values, 'Bubbles', gotit)
    if (.not. gotit) acc_bubblesdefault = GetString(GetSolverParams(), &
        'Stabilization method', gotit) == 'bubbles'
    if (.not. gotit) acc_bubblesdefault = .true.
  end subroutine ke_params

  subroutine ke_noslip() bind(c, name="ke_noslip")
    integer :: i
    type(valuelist_t), pointer :: bc
    logical :: gotit
    do i = 1, acc_model % numberofbcs
      bc => acc_model % bcs(i) % values
      if (GetLogical(bc, 'Noslip wall BC', gotit)) then
        call ListAddConstReal(bc, 'Kinetic Energy', 0.0_dp)
      end if
    end do
  end subroutine ke_noslip

  subroutine ke_alloc() bind(c, name="ke_alloc")
    integer :: n, istat, dofs
    if (.not. acc_allocdone) then
      n = acc_model % maxelementnodes
      dofs = acc_dofs
      allocate(acc_u(n), acc_v(n), acc_w(n), &
          acc_density(n), acc_work(n), &
          acc_viscosity(n), &
          acc_effectivevisc(2,n), &
          acc_turbvisc(n), acc_c0(dofs,n), &
          acc_layer(n), &
          acc_rough(n), &
          acc_kec1(n), acc_kec2(n), &
          acc_kesigk(n), acc_kesige(n), acc_kecmu(n), &
          acc_k(n), &
          acc_e(n), &
          acc_v2(n), acc_v2fct(n), &
          acc_nodaldensity(acc_solver % mesh % numberofnodes), &
          acc_nodalvisc(acc_solver % mesh % numberofnodes), &
          acc_nodalcmu(acc_solver % mesh % numberofnodes), &
          acc_prevk(n), acc_preve(n), &
          acc_xl(dofs*n), acc_xlprev(dofs*n), &
          acc_mass(2*dofs*n,2*dofs*n), &
          acc_stiff(2*dofs*n,2*dofs*n), acc_load(dofs,n), &
          acc_force(2*dofs*n), acc_timeforce(2*dofs*n), stat=istat)
      if (istat /= 0) then
        call Fatal('KESolver', 'Memory allocation error.')
      end if
      acc_allocdone = .true.
    end if
    if (acc_transient .and. .not. allocated(acc_bx)) then
      acc_bxstride = acc_dofs * max(acc_solver % mesh % maxbdofs, &
          acc_solver % mesh % maxelementnodes)
      allocate(acc_bx(acc_bxstride * (acc_solver % mesh % numberofbulkelements + &
          acc_solver % mesh % numberofboundaryelements)), &
          acc_bxprev(acc_bxstride * (acc_solver % mesh % numberofbulkelements + &
          acc_solver % mesh % numberofboundaryelements)))
      acc_bx = 0.0_dp
      acc_bxprev = 0.0_dp
    end if
    if (acc_transient .and. allocated(acc_bx) .and. &
        GetTimestep() /= acc_bubbletimestep) then
      acc_bxprev = acc_bx
      acc_bubbletimestep = GetTimestep()
    end if
  end subroutine ke_alloc

  function ke_niter() result(n) bind(c, name="ke_niter")
    integer(c_int) :: n
    n = acc_nonlineariter
  end function ke_niter

  ! Iteration head: Info reports, DefaultInitialize, output/advance
  ! setup. Returns GetNOFActive().
  function ke_iter_head(iter) result(noactive) bind(c, name="ke_iter_head")
    integer(c_int), value :: iter
    integer(c_int) :: noactive
    call Info('KESolver', ' ', level=4)
    call Info('KESolver', ' ', level=4)
    call Info('KESolver', '-------------------------------------', level=4)
    write(message, *) 'KEpsilon iteration: ', iter
    call Info('KESolver', message, level=4)
    call Info('KESolver', '-------------------------------------', level=4)
    call Info('KESolver', ' ', level=4)
    call Info('KESolver', 'Starting Assembly...', level=4)
    call DefaultInitialize()
    acc_body_id = -1
    call StartAdvanceOutput('KESolver', 'Assembly')
    noactive = GetNOFActive()
  end function ke_iter_head

  ! One bulk element: AdvanceOutput through DefaultUpdateEquations
  ! (material block, LocalMatrix, bubble condensation).
  subroutine ke_bulk_elem(t) bind(c, name="ke_bulk_elem")
    integer(c_int), value :: t
    integer :: i, n, nd, nb, boff, dofs, capn
    integer, pointer :: nodeindexes(:)
    logical :: bubbles, gotit
    type(variable_t), pointer :: ke
    dofs = acc_dofs
    capn = acc_model % maxelementnodes
    call AdvanceOutput(t, GetNOFActive())
    acc_elem => GetActiveElement(t)
    bubbles = acc_bubblesdefault .and. .not. associated(acc_elem % pdefs)
    acc_bubbles = bubbles
    if (acc_elem % bodyid /= acc_body_id) then
      acc_material => GetMaterial()
      acc_equation => GetEquation()
      acc_clip = GetConstReal(acc_material, 'KE Clip', gotit)
      if (.not. gotit) acc_clip = 1.0d-6
      acc_kemodel = GetString(acc_material, 'KE Model', gotit)
      if (.not. gotit) acc_kemodel = 'standard'
      acc_v2fmodel = GetString(acc_material, 'V2-F Model', gotit)
    end if
    call GetElementNodes(acc_enodes)
    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
    if (bubbles) nd = 2*n
    nb = GetElementNOFBDOFs()
    nodeindexes => acc_elem % nodeindexes
    call GetScalarLocalSolution(acc_v2, 'V2')
    call GetScalarLocalSolution(acc_k, 'Kinetic Energy')
    call GetScalarLocalSolution(acc_e, 'Kinetic Dissipation')
    if (acc_transient .and. &
        (bubbles .or. (nb > 0 .and. .not. acc_solver % globalbubbles))) then
      call GetScalarLocalSolution(acc_prevk, 'Kinetic Energy', tStep=-1)
      call GetScalarLocalSolution(acc_preve, 'Kinetic Dissipation', tStep=-1)
    end if
    call GetScalarLocalSolution(acc_u, 'Velocity 1')
    call GetScalarLocalSolution(acc_v, 'Velocity 2')
    call GetScalarLocalSolution(acc_w, 'Velocity 3')
    acc_kesigk(1:n) = GetReal(acc_material, 'KE SigmaK', gotit)
    if (.not. gotit) then
      acc_kesigk = 1.0d0
      call ListAddConstReal(acc_material, 'KE SigmaK', acc_kesigk(1))
    end if
    acc_kesige(1:n) = GetReal(acc_material, 'KE SigmaE', gotit)
    if (.not. gotit) then
      select case (acc_kemodel)
      case ('standard', 'v2-f')
        acc_kesige = 1.3_dp
      case ('rng')
        acc_kesige = 1.0_dp
      case default
        call Fatal('KESolver', 'Unknown K-Epsilon model')
      end select
      call ListAddConstReal(acc_material, 'KE SigmaE', acc_kesige(1))
    end if
    acc_kecmu(1:n) = ListGetConstReal(acc_material, 'KE Cmu', gotit)
    if (.not. gotit) then
      select case (acc_kemodel)
      case ('standard')
        acc_kecmu = 0.09_dp
      case ('v2-f')
        acc_kecmu = 0.22_dp
      case ('rng')
        acc_kecmu = 0.0845_dp
      case default
        call Fatal('KESolver', 'Unknown K-Epsilon model')
      end select
      call ListAddConstReal(acc_material, 'KE Cmu', acc_kecmu(1))
    end if
    acc_kec1(1:n) = GetReal(acc_material, 'KE C1', gotit)
    if (.not. gotit) then
      select case (acc_kemodel)
      case ('standard')
        acc_kec1 = 1.44_dp
      case ('v2-f')
        acc_kec1(1:n) = 1.4_dp
      case ('rng')
        acc_kec1 = 1.42_dp
      case default
        call Fatal('KESolver', 'Unknown K-Epsilon model')
      end select
      call ListAddConstReal(acc_material, 'KE C1', acc_kec1(1))
    end if
    acc_kec2(1:n) = GetReal(acc_material, 'KE C2', gotit)
    if (.not. gotit) then
      select case (acc_kemodel)
      case ('standard')
        acc_kec2 = 1.92_dp
      case ('v2-f')
        acc_kec2 = 1.9_dp
      case ('rng')
        acc_kec2 = 1.68_dp
      case default
        call Fatal('KESolver', 'Unknown K-Epsilon model')
      end select
      call ListAddConstReal(acc_material, 'KE C2', acc_kec2(1))
    end if
    if (acc_kemodel == 'v2-f') then
      acc_v2fct(1:n) = GetReal(acc_material, 'V2-F CT', gotit)
      if (.not. gotit) then
        acc_v2fct(1:n) = 6.0_dp
        call ListAddConstReal(acc_material, 'V2-F CT', acc_v2fct(1))
      end if
      acc_v2fcp = GetCReal(acc_material, 'V2-F Cp', gotit)
      if (.not. gotit) acc_v2fcp = 0.05_dp
    end if
    acc_density(1:n) = GetReal(acc_material, 'Density')
    acc_viscosity(1:n) = GetReal(acc_material, 'Viscosity')
    acc_nodaldensity(nodeindexes(1:n)) = acc_density(1:n)
    acc_nodalvisc(nodeindexes(1:n)) = acc_viscosity(1:n)
    acc_nodalcmu(nodeindexes(1:n)) = acc_kecmu(1:n)
    call ke_localmatrix(acc_mass, acc_stiff, acc_force, acc_load, &
        acc_u, acc_v, acc_w, acc_elem, n, nd+nb, acc_enodes)
    acc_timeforce = 0.0_dp
    if (bubbles) then
      if (acc_transient) then
        acc_xl(1:2*n-1:2) = acc_k(1:n)
        acc_xl(2:2*n:2) = acc_e(1:n)
        acc_xlprev(1:2*n-1:2) = acc_prevk(1:n)
        acc_xlprev(2:2*n:2) = acc_preve(1:n)
        boff = (acc_elem % elementindex - 1) * acc_bxstride
        call CondensatePTransient(n, n, dofs, acc_dt, acc_mass, acc_stiff, &
            acc_force, acc_xlprev(1:dofs*n), acc_xl(1:dofs*n), &
            acc_bxprev(boff+1:boff+dofs*n), acc_bx(boff+1:boff+dofs*n))
      else
        if (acc_transient) call Default1stOrderTime(acc_mass, acc_stiff, acc_force)
        call Condensate(dofs*capn, acc_stiff, acc_force, acc_timeforce)
      end if
    else if (nb > 0) then
      if (acc_transient .and. .not. acc_solver % globalbubbles) then
        acc_xl(1:2*nd-1:2) = acc_k(1:nd)
        acc_xl(2:2*nd:2) = acc_e(1:nd)
        acc_xlprev(1:2*nd-1:2) = acc_prevk(1:nd)
        acc_xlprev(2:2*nd:2) = acc_preve(1:nd)
        boff = (acc_elem % elementindex - 1) * acc_bxstride
        call CondensatePTransient(nd, nb, dofs, acc_dt, acc_mass, acc_stiff, &
            acc_force, acc_xlprev(1:dofs*nd), acc_xl(1:dofs*nd), &
            acc_bxprev(boff+1:boff+dofs*nb), acc_bx(boff+1:boff+dofs*nb))
      else
        if (acc_transient) call Default1stOrderTime(acc_mass, acc_stiff, acc_force)
        call CondensateP(dofs*nd, dofs*nb, acc_stiff, acc_force, acc_timeforce)
      end if
    else
      if (acc_transient) call Default1stOrderTime(acc_mass, acc_stiff, acc_force)
    end if
    call DefaultUpdateEquations(acc_stiff, acc_force)
  end subroutine ke_bulk_elem

  subroutine ke_bulk_tail() bind(c, name="ke_bulk_tail")
    call DefaultFinishBulkAssembly()
    call Info('KESolver', 'Assembly done', level=4)
  end subroutine ke_bulk_tail

  subroutine ke_wall_zero() bind(c, name="ke_wall_zero")
    integer :: t, n, i, j
    integer, pointer :: kinperm(:)
    type(variable_t), pointer :: ke
    logical :: gotit
    if (ListGetLogicalAnyBC(acc_model, 'Epsilon Wall BC') .or. &
        ListGetLogicalAnyBC(acc_model, 'Noslip Wall BC')) then
      ke => acc_solver % variable
      kinperm => ke % perm
      do t = 1, acc_solver % mesh % numberofboundaryelements
        acc_elem => GetBoundaryElement(t)
        n = GetElementNOFNodes()
        acc_bc => GetBC()
        if (associated(acc_bc)) then
          if (GetLogical(acc_bc, 'Epsilon Wall BC', gotit) .or. &
              GetLogical(acc_bc, 'Noslip Wall BC', gotit)) then
            do i = 1, n
              j = kinperm(acc_elem % nodeindexes(i))
              call ZeroRow(acc_solver % matrix, 2*j)
              acc_solver % matrix % rhs(2*j) = 0.0_dp
            end do
          end if
        end if
      end do
    end if
  end subroutine ke_wall_zero

  subroutine ke_boundary() bind(c, name="ke_boundary")
    integer :: t, n, j, k, dofs
    integer, pointer :: nodeindexes(:)
    integer, pointer :: flowperm(:)
    real(kind=dp), pointer :: flowsolution(:)
    type(variable_t), pointer :: flowsol
    logical :: gotit
    dofs = acc_dofs
    do t = 1, acc_solver % mesh % numberofboundaryelements
      acc_elem => GetBoundaryElement(t)
      if (.not. ActiveBoundaryElement()) cycle
      n = GetElementNOFNodes()
      nodeindexes => acc_elem % nodeindexes
      acc_bc => GetBC()
      if (associated(acc_bc)) then
        if (GetLogical(acc_bc, 'Wall Law', gotit)) then
          acc_density(1:n) = GetParentMatProp('Density')
          acc_viscosity(1:n) = GetParentMatProp('Viscosity')
          acc_rough(1:n) = GetReal(acc_bc, 'Surface Roughness', gotit)
          acc_layer(1:n) = GetReal(acc_bc, 'Boundary Layer Thickness', gotit)
          flowsol => VariableGet(acc_model % variables, 'Flow Solution')
          flowperm => flowsol % perm
          flowsolution => flowsol % values
          do j = 1, n
            k = flowperm(nodeindexes(j))
            if (k > 0) then
              select case (acc_nsdofs)
              case (3)
                acc_u(j) = flowsolution(acc_nsdofs*k-2)
                acc_v(j) = flowsolution(acc_nsdofs*k-1)
                acc_w(j) = 0.0d0
              case (4)
                acc_u(j) = flowsolution(acc_nsdofs*k-3)
                acc_v(j) = flowsolution(acc_nsdofs*k-2)
                acc_w(j) = flowsolution(acc_nsdofs*k-1)
              end select
            else
              acc_u(j) = 0.0d0
              acc_v(j) = 0.0d0
              acc_w(j) = 0.0d0
            end if
          end do
          do j = 1, n
            call KEWall(acc_work(1), acc_work(2), acc_work(3), &
                sqrt(acc_u(j)**2+acc_v(j)**2+acc_w(j)**2), &
                acc_layer(j), acc_rough(j), acc_viscosity(j), &
                acc_density(j))
            k = dofs*(ke_kinperm(nodeindexes(j))-1)
            call UpdateDirichletDof(acc_solver % matrix, k+1, acc_work(1))
            call UpdateDirichletDof(acc_solver % matrix, k+2, acc_work(2))
          end do
        end if
        if (GetLogical(acc_bc, 'Epsilon Wall BC', gotit) .or. &
            GetLogical(acc_bc, 'Noslip Wall BC', gotit)) then
          call ke_epsilonwall(acc_elem, n, acc_stiff, acc_force)
        end if
      end if
    end do
    call DefaultFinishAssembly()
    call DefaultDirichletBCs()
    call Info('KESolver', 'Set boundaries done', level=4)
  end subroutine ke_boundary

  function ke_kinperm(node) result(k)
    integer :: node, k
    type(variable_t), pointer :: ke
    ke => acc_solver % variable
    k = ke % perm(node)
  end function ke_kinperm

  ! Solve tail: DefaultFinishAssembly through the convergence check.
  ! Returns 1 when NonlinConverged == 1 (exit the nonlinear loop).
  function ke_solve_tail() result(exitloop) bind(c, name="ke_solve_tail")
    integer(c_int) :: exitloop
    integer :: i, n, k
    real(kind=dp) :: kmax, emax, kval, eval
    type(variable_t), pointer :: ke
    acc_prevnorm = acc_norm
    acc_norm = DefaultSolve()
    n = acc_solver % mesh % numberofnodes
    kmax = maxval(acc_solver % variable % values(1::2))
    emax = maxval(acc_solver % variable % values(2::2))
    ke => acc_solver % variable
    do i = 1, n
      k = ke % perm(i)
      if (k <= 0) cycle
      kval = acc_solver % variable % values(2*k-1)
      eval = acc_solver % variable % values(2*k-0)
      if (kval < acc_clip*kmax) kval = acc_clip*kmax
      if (eval < acc_clip*emax) then
        eval = max(acc_nodaldensity(i)*acc_nodalcmu(i)*kval**2/ &
            acc_nodalvisc(i), acc_clip*emax)
      end if
      acc_solver % variable % values(2*k-1) = max(kval, 1.0d-10)
      acc_solver % variable % values(2*k-0) = max(eval, 1.0d-10)
    end do
    acc_relchange = acc_solver % variable % nonlinchange
    write(message, *) 'Result Norm   : ', acc_norm
    call Info('KESolver', message, level=4)
    write(message, *) 'Relative Change : ', acc_relchange
    call Info('KESolver', message, level=4)
    if (acc_solver % variable % nonlinconverged == 1) then
      exitloop = 1
    else
      exitloop = 0
    end if
  end function ke_solve_tail

  subroutine ke_finish() bind(c, name="ke_finish")
    integer :: n
    type(variable_t), pointer :: ke
    n = size(acc_solver % variable % values)
    ke => VariableGet(acc_solver % mesh % variables, 'Kinetic Energy')
    if (associated(ke)) ke % values = acc_solver % variable % values(1:n:2)
    ke => VariableGet(acc_solver % mesh % variables, 'Kinetic Dissipation')
    if (associated(ke)) ke % values = acc_solver % variable % values(2:n:2)
  end subroutine ke_finish
  ! The LocalMatrix internal procedure: element local matrices and
  ! right-hand side for the diffusion-convection equations. Transcribed
  ! statement for statement; host data (Model, Material, KEModel,
  ! Bubbles, Density, Viscosity, local solutions, V2FCp) arrives via
  ! module state. Like the original it USEs MaterialModels.
  subroutine ke_localmatrix(mass, stiff, force, load, ux, uy, uz, element, n, nd, nodes)
    use materialmodels
    real(kind=dp), dimension(:) :: force, ux, uy, uz
    real(kind=dp), dimension(:,:) :: mass, stiff, load
    integer :: n, nd
    type(nodes_t) :: nodes
    type(element_t) :: element
    real(kind=dp) :: ddbasisddx(nd,3,3)
    real(kind=dp) :: basis(nd)
    real(kind=dp) :: dbasisdx(nd,3), detj
    real(kind=dp) :: velo(3), dvelodx(3,3)
    real(kind=dp) :: a(2,2), m(2,2), prodk, prode, div, prodtensor(3,3)
    integer :: i, j, c, p, q, t, dim, nbasis
    real(kind=dp) :: loadatip(2), cmu, rho, mu, tmu, effmu(2), timescale, re_t, lc1, lc2, lc3, pr
    real(kind=dp) :: s, u, v, w, k, e, eta, strain(3,3), gradasymm(3,3), nu, cv, &
        alpha, oldalpha, dalpha, err, ww, olderr, derr
    type(gaussintegrationpoints_t), target :: integstuff
    real(kind=dp) :: specificheatratio, pressure(n), referencepressure, &
        sound_speed_sq, mach_number_sq
    real(kind=dp) :: metric(3,3), symb(3,3,3), dsymb(3,3,3,3), &
        sqrtmetric, gravity(3), pr_rho(n), rho_g, c3(n)
    real(kind=dp) :: c0(2), c1, ct, c2(2), dcdx(3), secinv, x, y, z, lv2, lct, sigmak, sigmae
    real(kind=dp), pointer :: gwork(:,:)
    logical :: stat, convection, userngmodel
    dim = CoordinateSystemDimension()
    force = 0.0d0
    stiff = 0.0d0
    mass = 0.0d0
    nbasis = nd
    if (acc_bubbles) nbasis = 2*n
    userngmodel = acc_kemodel == 'rng'
    gwork => ListGetConstRealArray(acc_model % constants, 'Gravity', stat)
    if (stat) then
      gravity = gwork(1:3,1)*gwork(4,1)
    else
      gravity = 0.00_dp
      gravity(2) = -9.81_dp
    end if
    pr_rho(1:n) = GetReal(acc_material, 'Turbulent Prandtl Number', stat)
    if (.not. stat) pr_rho(1:n) = 0.85_dp
    c3(1:n) = GetReal(acc_material, 'Dissipation buoyancy coefficient', stat)
    if (.not. stat) c3(1:n) = 0.0_dp
    call ElementDensity(acc_density, n)
    specificheatratio = GetCReal(acc_material, 'Specific Heat Ratio', stat)
    call getScalarLocalSolution(pressure, 'Pressure')
    referencepressure = GetCReal(acc_material, 'Reference Pressure', stat)
    if (acc_bubbles) then
      integstuff = GaussPoints(element, element % type % gaussPoints2)
    else
      integstuff = GaussPoints(element)
    end if
    do t = 1, integstuff % n
      u = integstuff % u(t)
      v = integstuff % v(t)
      w = integstuff % w(t)
      stat = ElementInfo(element, nodes, u, v, w, detj, basis, dbasisdx, bubbles=acc_bubbles)
      s = detj * integstuff % s(t)
      if (CurrentCoordinateSystem() /= cartesian) then
        x = sum(nodes % x(1:n)*basis(1:n))
        y = sum(nodes % y(1:n)*basis(1:n))
        z = sum(nodes % z(1:n)*basis(1:n))
        call CoordinateSystemInfo(metric, sqrtmetric, symb, dsymb, x, y, z)
        s = s * sqrtmetric
      end if
      velo = 0.0d0
      velo(1) = sum(ux(1:n)*basis(1:n))
      velo(2) = sum(uy(1:n)*basis(1:n))
      velo(3) = sum(uz(1:n)*basis(1:n))
      dvelodx = 0.0d0
      do i = 1, dim
        dvelodx(1,i) = sum(ux(1:n)*dbasisdx(1:n,i))
        dvelodx(2,i) = sum(uy(1:n)*dbasisdx(1:n,i))
        dvelodx(3,i) = sum(uz(1:n)*dbasisdx(1:n,i))
      end do
      if (CurrentCoordinateSystem() == cartesian) then
        strain = 0.5d0 * (dvelodx + transpose(dvelodx))
        secinv = 2 * sum(strain * strain)
      else
        secinv = SecondInvariant(velo, dvelodx, metric, symb) / 2
      end if
      k = sum(acc_k(1:n) * basis(1:n))
      e = sum(acc_e(1:n) * basis(1:n))
      eta = sqrt(secinv) * k / e
      pr = sum(pr_rho(1:n) * basis(1:n))
      rho_g = 0._dp
      do i = 1, dim
        rho_g = rho_g + sum(acc_density(1:n) * dbasisdx(1:n,i)) * gravity(i)
      end do
      mu = sum(acc_viscosity(1:n) * basis(1:n))
      rho = sum(acc_density(1:n) * basis(1:n))
      sigmak = sum(acc_kesigk(1:n) * basis(1:n))
      sigmae = sum(acc_kesige(1:n) * basis(1:n))
      cmu = sum(acc_kecmu(1:n) * basis(1:n))
      lc1 = sum(acc_kec1(1:n) * basis(1:n))
      lc2 = sum(acc_kec2(1:n) * basis(1:n))
      lc3 = sum(c3(1:n) * basis(1:n))
      sound_speed_sq = (sum(basis(1:n)*pressure(1:n))+referencepressure) * &
          specificheatratio / rho
      mach_number_sq = 0._dp
      if (sound_speed_sq > 0._dp) mach_number_sq = k/sound_speed_sq
      if (acc_kemodel == 'v2-f') then
        lv2 = sum(acc_v2(1:n) * basis(1:n))
        lct = sum(acc_v2fct(1:n) * basis(1:n))
        lc1 = 1.4_dp * (1+acc_v2fcp*sqrt(k/lv2))
        timescale = max(k/e, lct*sqrt(mu/rho/e))
        tmu = rho * cmu * lv2 * timescale
      else
        tmu = rho * cmu*k**2 / e
      end if
      prodk = tmu * (secinv-rho_g/(rho*pr)) / rho
      prode = tmu * (secinv-lc3*rho_g/(rho*pr)) / rho
      effmu(1) = mu + tmu / sigmak
      effmu(2) = mu + tmu / sigmae
      c0(1) = rho
      if (acc_kemodel == 'v2-f') then
        c0(2) = rho * lc2 / timescale
      else
        c0(2) = rho * lc2 * e / k
      end if
      c1 = rho
      ct = rho
      alpha = 1.0d0
      if (userngmodel) then
        ww = mu / effmu(1)
        alpha = 1.3929d0
        oldalpha = 1
        olderr = abs((oldalpha-1.3929d0)/(1.0d0-1.3929d0))**0.6321d0
        olderr = olderr * abs((oldalpha+2.3929d0)/(1.0d0+2.3929d0))**0.3679d0
        olderr = olderr - ww
        do i = 1, 100
          err = abs((alpha-1.3929d0)/(1.0d0-1.3929d0))**0.6321d0
          err = err * abs((alpha+2.3929d0)/(1.0d0+2.3929d0))**0.3679d0
          err = err - ww
          derr = olderr - err
          olderr = err
          dalpha = oldalpha - alpha
          oldalpha = alpha
          alpha = alpha - 0.5 * err * dalpha / derr
          if (abs(err) < 1.0d-8) exit
        end do
        if (abs(err) > 1.0d-8) then
          print *, 'huh: ', alpha
          alpha = 1.3929d0
        end if
      end if
      c2(1) = alpha * effmu(1)
      c2(2) = alpha * effmu(2)
      do p = 1, nbasis
        do q = 1, nbasis
          m = 0.0d0
          a = 0.0d0
          m(1,1) = ct * basis(q) * basis(p)
          m(2,2) = ct * basis(q) * basis(p)
          a(1,2) = c0(1) * basis(q) * basis(p)
          a(2,2) = c0(2) * basis(q) * basis(p)
          if (CurrentCoordinateSystem() == cartesian) then
            do i = 1, dim
              a(1,1) = a(1,1) + c2(1) * dbasisdx(q,i) * dbasisdx(p,i)
              a(2,2) = a(2,2) + c2(2) * dbasisdx(q,i) * dbasisdx(p,i)
            end do
          else
            do i = 1, dim
              do j = 1, dim
                a(1,1) = a(1,1) + metric(i,j) * c2(1) * &
                    dbasisdx(q,i) * dbasisdx(p,j)
                a(2,2) = a(2,2) + metric(i,j) * c2(2) * &
                    dbasisdx(q,i) * dbasisdx(p,j)
              end do
            end do
          end if
          do i = 1, dim
            a(1,1) = a(1,1) + c1 * velo(i) * dbasisdx(q,i) * basis(p)
            a(2,2) = a(2,2) + c1 * velo(i) * dbasisdx(q,i) * basis(p)
          end do
          do i = 1, 2
            do j = 1, 2
              stiff(2*(p-1)+i,2*(q-1)+j) = stiff(2*(p-1)+i,2*(q-1)+j)+s*a(i,j)
              mass(2*(p-1)+i,2*(q-1)+j) = mass(2*(p-1)+i,2*(q-1)+j) +s*m(i,j)
            end do
          end do
        end do
      end do
      loadatip(1) = rho*(prodk-2*e*mach_number_sq)
      if (acc_kemodel == 'v2-f') then
        loadatip(2) = rho*lc1*prode/timescale
      else
        loadatip(2) = rho*lc1*prode*e/k
      end if
      if (userngmodel) &
          loadatip(2) = loadatip(2) - cmu*rho*eta**3*(1-eta/4.38d0) / &
          (1.0d0 + 0.012d0*eta**3) * e**2 / k
      do p = 1, nbasis
        force(2*(p-1)+1) = force(2*(p-1)+1)+s*loadatip(1)*basis(p)
        force(2*(p-1)+2) = force(2*(p-1)+2)+s*loadatip(2)*basis(p)
      end do
    end do
  end subroutine ke_localmatrix

  ! The EpsilonWall internal procedure: wall-law epsilon boundary
  ! assembly. Transcribed statement for statement; host data (BC, U, V,
  ! W, Density, Viscosity, TransientSimulation, dt, MASS) arrives via
  ! module state, and its SAVE Nodes/ParentNodes live module-side with
  ! identical lifetimes.
  subroutine ke_epsilonwall(element, n, stiff, force)
    type(element_t) :: element
    integer :: n
    real(kind=dp) :: stiff(:,:), force(:)
    real(kind=dp) :: basisb(32), basis(32), basisk(32), dbasisdx(32,3), &
        evals(32), kvals(32), detj
    integer :: i, j, c, p, q, t, np, dim, n_integ, nbasis
    real(kind=dp) :: s, u, v, w, e, k, kder, x, y, z, normal(3), mu, rho, relax
    type(gaussintegrationpoints_t), target :: integstuff
    logical :: stat
    type(element_t), pointer :: parent
    dim = CoordinateSystemDimension()
    force = 0.0_dp
    stiff = 0.0_dp
    integstuff = GaussPoints(element)
    relax = GetCReal(acc_bc, 'Epsilon Relax', stat)
    if (.not. stat) relax = 1
    parent => element % boundaryinfo % left
    if (.not. associated(parent)) &
        parent => element % boundaryinfo % right
    if (.not. associated(parent)) return
    np = GetElementNOFDOFs(parent)
    call GetElementNodes(acc_wnodes)
    call GetElementNodes(acc_wpnodes, parent)
    acc_density(1:np) = GetReal(GetMaterial(parent), 'Density', UElement=parent)
    acc_viscosity(1:np) = GetReal(GetMaterial(parent), 'Viscosity', UElement=parent)
    call GetScalarLocalSolution(kvals, 'Kinetic Energy', parent)
    call GetScalarLocalSolution(evals, 'Kinetic dissipation', parent)
    do t = 1, integstuff % n
      u = integstuff % u(t)
      v = integstuff % v(t)
      w = integstuff % w(t)
      stat = ElementInfo(element, acc_wnodes, u, v, w, detj, basisb)
      s = detj * integstuff % s(t)
      if (CurrentCoordinateSystem() /= cartesian) then
        x = sum(acc_wnodes % x(1:n)*basisb(1:n))
        y = sum(acc_wnodes % y(1:n)*basisb(1:n))
        z = sum(acc_wnodes % z(1:n)*basisb(1:n))
        s = s * CoordinateSqrtMetric(x, y, z)
      end if
      ! Sequence association (as in the original, which passes the whole
      ! U/V/W arrays to these scalar dummies): only element 1 travels.
      normal = NormalVector(element, acc_wnodes, acc_u(1), acc_v(1), .true.)
      call GetParentUVW(element, n, parent, np, acc_u(1), acc_v(1), acc_w(1), basisb)
      stat = ElementInfo(parent, acc_wpnodes, acc_u(1), acc_v(1), acc_w(1), detj, &
          basis, dbasisdx)
      if (abs(u) > 0.999_dp) then
        u = -u
      else if (abs(v) > 0.999_dp) then
        v = -v
      else if (abs(w) > 0.999_dp) then
        w = -w
      end if
      stat = ElementInfo(parent, acc_wpnodes, acc_u(1), acc_v(1), acc_w(1), detj, basisk)
      rho = sum(basis(1:np) * acc_density(1:np))
      mu = sum(basis(1:np) * acc_viscosity(1:np))
      e = sum(basis(1:np) * evals(1:np))
      k = sum(basisk(1:np) * kvals(1:np))
      kvals(1:np) = sqrt(kvals(1:np))
      kder = 0.0_dp
      do i = 1, 3
        kder = kder+sum(dbasisdx(1:np,i)*kvals(1:np))*normal(i)
      end do
      do p = 1, np
        do q = 1, np
          stiff(2*p,2*q) = stiff(2*p,2*q) + s*basis(q)*basis(p)
          stiff(2*p,2*q-1) = stiff(2*p,2*q-1) - s*relax*2*mu/rho*kder**2*basisk(q)/max(k,aeps)*basis(p)
        end do
        force(2*p) = force(2*p) + s*(1-relax)*e*basis(p)
      end do
    end do
    if (acc_transient) then
      acc_mass = 0.0_dp
      call Default1stOrderTime(acc_mass, stiff, force, UElement=parent)
    end if
    call DefaultUpdateEquations(stiff, force, UElement=parent)
  end subroutine ke_epsilonwall

  ! KESolverLegacy_Init body: solver-list defaults.
  subroutine ke_init_begin(h_solver, transient) bind(c, name="ke_init_begin")
    type(c_ptr), value :: h_solver
    integer(c_int), value :: transient
    call c_f_pointer(h_solver, acc_solver)
    acc_transient = (transient /= 0)
  end subroutine ke_init_begin

  subroutine ke_init_body() bind(c, name="ke_init_body")
    type(valuelist_t), pointer :: solverparams
    logical :: found, pbubble, legacybubbles
    character(len=MAX_NAME_LEN) :: str
    solverparams => GetSolverParams()
    str = GetString(solverparams, 'Variable', found)
    if (.not. found) str = 'K-eps'
    if (index(str, '[') <= 0) then
      call ListAddString(solverparams, 'Variable', &
          trim(str) // '[Kinetic Energy:1 Kinetic Dissipation:1]')
    end if
    str = ListGetString(solverparams, 'Element', found)
    pbubble = .false.
    if (found) pbubble = index(str, 'b:') > 0
    if (pbubble) then
      call ListAddNewLogical(solverparams, 'Bubbles in Global System', .false.)
      if (acc_transient .and. &
          .not. ListGetLogical(solverparams, 'Bubbles in Global System', found)) then
        call ListAddNewInteger(solverparams, 'Nonlinear System Min Iterations', 2)
        call ListAddNewInteger(solverparams, 'Nonlinear System Max Iterations', 2)
      end if
    else if (acc_transient) then
      legacybubbles = ListGetLogical(solverparams, 'Bubbles', found)
      if (.not. found) legacybubbles = ListGetString(solverparams, &
          'Stabilization method', found) == 'bubbles'
      if (.not. found) legacybubbles = .true.
      if (legacybubbles) then
        call ListAddNewInteger(solverparams, 'Nonlinear System Min Iterations', 2)
        call ListAddNewInteger(solverparams, 'Nonlinear System Max Iterations', 2)
      end if
    end if
  end subroutine ke_init_body

end module kesolverlegacy_acc
