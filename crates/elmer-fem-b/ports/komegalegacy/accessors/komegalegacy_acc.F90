! T022-A KomegaLegacy port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/TurbulenceModels/KomegaLegacy.F90,
! compiled against the built .mod files of the target configure (A2
! shipped mechanism, ADR-057 section 3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars). All Fortran-object traffic — TYPE(...)
! handles, OPTIONAL/POINTER dummies, assumed-shape descriptors, hidden
! lengths, list-directed reports — stays inside Fortran, so the Rust
! side is pure C calls. The wrapper module mirrors KomegaLegacy's own
! locals (entry SAVE arrays with identical first-mesh sizing and
! lifetimes, the per-element bubble history bx/bxprev with the same
! timestep-shift guard, the LocalMatrix/OmegaWall SAVE Nodes) and all
! floating-point work (the LocalMatrix Galerkin assembly, the wall
! distance MINVAL, the K/Omega clip floors) plus all formatted I/O
! (WRITE(Message,*)/Info) stays in Fortran (same compiler), so
! bit-identity holds by construction wherever the call sequence
! matches.
!
! Faithfulness notes (mirrors KomegaLegacy.F90 statement for
! statement):
! - Dead stores replicated as documented skips, never executed: `at`/
!   `at0` (CPUTime/RealTime assigned, never read), `NewtonTol`/
!   `NewtonIter` (read, never used), `body_id` (set to -1, never read),
!   `FlowSol`, `StiffMatrix`, `Stabilize`, `NewtonLinearization`,
!   `eq_id`, `bf_id`, `KMax`/`EMax` (MAXVAL reductions whose results
!   are never used — like DivergenceSolver's `Visited`), exactly as in
!   the original where they have no observable effect.
! - `Norm = KE % Norm` is a pure component read (overwritten by
!   DefaultSolve before any use); replicated as a state save.
! - `dt`/`TransientSimulation` cross by value; `dt` feeds only the two
!   CondensatePTransient calls (Default1stOrderTime reads its own
!   state), exactly as in the original. `Model` is received and stored
!   (its components are read in wrappers, as in the original).
! - No Rust threads: the solver entries run serially, so the `Mutex`
!   is only a `static`-sharing guard, never contended across threads.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module komegalegacy_acc
  use iso_c_binding
  use coordinatesystems
  use defutils
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  real(kind=dp), save :: acc_dt = 0.0_dp
  logical, save :: acc_transient = .false.

  ! Entry scalar state (mirrors the solver body locals).
  integer, save :: acc_dofs = 0, acc_localnodes = 0
  logical, save :: acc_bubbles = .false.
  integer, save :: acc_newtoniter = 0, acc_nonlineariter = 1
  integer, save :: acc_n0 = 0
  integer, save :: acc_bxstride = 0, acc_bubbletimestep = -1
  real(kind=dp), save :: acc_newtontol = 0.0_dp
  real(kind=dp), save :: acc_norm = 0.0_dp, acc_relchange = 0.0_dp
  logical, save :: acc_bubblesdefault = .true.
  logical, save :: acc_allocdone = .false.

  ! Entry SAVE arrays (first-mesh sizing, never resized — like the
  ! original, which assumes the first mesh is the biggest).
  real(kind=dp), allocatable, save :: acc_mass(:,:), acc_stiff(:,:)
  real(kind=dp), allocatable, save :: acc_load(:,:), acc_force(:)
  real(kind=dp), allocatable, save :: acc_timeforce(:)
  real(kind=dp), allocatable, save :: acc_lke(:), acc_ldis(:)
  real(kind=dp), allocatable, save :: acc_pke(:), acc_pdis(:)
  real(kind=dp), allocatable, save :: acc_xl(:), acc_xlprev(:)
  real(kind=dp), allocatable, save :: acc_bx(:), acc_bxprev(:)

  ! Element-loop handles (re-acquired per call, like the original's
  ! per-element pointer assignments) and SAVE Nodes state.
  type(element_t), pointer, save :: acc_elem => null()
  type(valuelist_t), pointer, save :: acc_bc => null()
  type(valuelist_t), pointer, save :: acc_material => null()
  type(nodes_t), save :: acc_enodes

contains

  subroutine ko_begin(h_model, h_solver, dt, transient) bind(c, name="ko_begin")
    type(c_ptr), value :: h_model, h_solver
    real(c_double), value :: dt
    integer(c_int), value :: transient
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    acc_dt = dt
    acc_transient = (transient /= 0)
  end subroutine ko_begin

  ! Guards: no-matrix and no-local-nodes silent returns, exactly as in
  ! the original (whose KinPerm reference would fault identically when
  ! Solver % Variable is unassociated). Returns 0 to proceed, 1 to
  ! return.
  function ko_prologue() result(code) bind(c, name="ko_prologue")
    integer(c_int) :: code
    type(variable_t), pointer :: ke
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
    code = 0
  end function ko_prologue

  subroutine ko_params() bind(c, name="ko_params")
    logical :: gotit
    acc_newtontol = ListGetConstReal(acc_solver % values, &
        'Nonlinear System Newton After Tolerance', gotit)
    acc_newtoniter = ListGetInteger(acc_solver % values, &
        'Nonlinear System Newton After Iterations', gotit)
    acc_nonlineariter = ListGetInteger(acc_solver % values, &
        'Nonlinear System Max Iterations', gotit)
    if (.not. gotit) acc_nonlineariter = 1
    acc_bubblesdefault = ListGetLogical(acc_solver % values, 'Bubbles', gotit)
    if (.not. gotit) acc_bubblesdefault = .true.
  end subroutine ko_params

  subroutine ko_noslip() bind(c, name="ko_noslip")
    integer :: i
    type(valuelist_t), pointer :: bc
    logical :: gotit
    do i = 1, acc_model % numberofbcs
      bc => acc_model % bcs(i) % values
      if (GetLogical(bc, 'Noslip wall BC', gotit)) then
        call ListAddConstReal(bc, 'Kinetic Energy', 0.0_dp)
      end if
    end do
  end subroutine ko_noslip

  subroutine ko_alloc() bind(c, name="ko_alloc")
    integer :: n, istat, dofs
    if (.not. acc_allocdone) then
      n = acc_solver % mesh % maxelementdofs
      dofs = acc_dofs
      acc_n0 = n
      allocate(acc_mass(2*dofs*n,2*dofs*n), &
          acc_stiff(2*dofs*n,2*dofs*n), acc_load(dofs,n), &
          acc_force(2*dofs*n), acc_timeforce(2*dofs*n), &
          acc_lke(n), acc_ldis(n), &
          acc_pke(n), acc_pdis(n), &
          acc_xl(dofs*n), acc_xlprev(dofs*n), stat=istat)
      if (istat /= 0) then
        call Fatal('KOmega', 'Memory allocation error.')
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
  end subroutine ko_alloc

  function ko_niter() result(n) bind(c, name="ko_niter")
    integer(c_int) :: n
    n = acc_nonlineariter
  end function ko_niter

  ! Iteration head: Info reports, DefaultInitialize, output/advance
  ! setup. Returns GetNOFActive().
  function ko_iter_head(iter) result(noactive) bind(c, name="ko_iter_head")
    integer(c_int), value :: iter
    integer(c_int) :: noactive
    call Info('KOmega', ' ', level=4)
    call Info('KOmega', ' ', level=4)
    call Info('KOmega', '-------------------------------------', level=4)
    write(message, *) 'Komega iteration: ', iter
    call Info('KOmega', message, level=4)
    call Info('KOmega', '-------------------------------------', level=4)
    call Info('KOmega', ' ', level=4)
    call Info('KOmega', 'Starting Assembly...', level=4)
    call DefaultInitialize()
    call StartAdvanceOutput('Komega', 'Assembly:')
    noactive = GetNOFActive()
  end function ko_iter_head

  ! One bulk element: AdvanceOutput through DefaultUpdateEquations
  ! (LocalMatrix, bubble condensation). Exact statement order of the
  ! original element loop.
  subroutine ko_bulk_elem(t) bind(c, name="ko_bulk_elem")
    integer(c_int), value :: t
    integer :: n, nd, nb, boff, dofs
    logical :: bubbles
    dofs = acc_dofs
    call AdvanceOutput(t, GetNOFActive())
    acc_elem => GetActiveElement(t)
    bubbles = acc_bubblesdefault .and. .not. associated(acc_elem % pdefs)
    acc_bubbles = bubbles
    acc_material => GetMaterial()
    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
    if (bubbles) nd = 2*n
    nb = GetElementNOFBDOFs()
    call GetElementNodes(acc_enodes)
    if (acc_transient .and. &
        (bubbles .or. (nb > 0 .and. .not. acc_solver % globalbubbles))) then
      call GetScalarLocalSolution(acc_lke, 'Kinetic energy')
      call GetScalarLocalSolution(acc_ldis, 'Kinetic Dissipation')
      call GetScalarLocalSolution(acc_pke, 'Kinetic energy', tStep=-1)
      call GetScalarLocalSolution(acc_pdis, 'Kinetic Dissipation', tStep=-1)
    end if
    call ko_localmatrix(n, nd+nb)
    acc_timeforce = 0.0_dp
    if (bubbles) then
      if (acc_transient) then
        acc_xl(1:2*n-1:2) = acc_lke(1:n)
        acc_xl(2:2*n:2) = acc_ldis(1:n)
        acc_xlprev(1:2*n-1:2) = acc_pke(1:n)
        acc_xlprev(2:2*n:2) = acc_pdis(1:n)
        boff = (acc_elem % elementindex - 1) * acc_bxstride
        call CondensatePTransient(n, n, dofs, acc_dt, acc_mass, acc_stiff, &
            acc_force, acc_xlprev(1:dofs*n), acc_xl(1:dofs*n), &
            acc_bxprev(boff+1:boff+dofs*n), acc_bx(boff+1:boff+dofs*n))
      else
        if (acc_transient) call Default1stOrderTime(acc_mass, acc_stiff, acc_force)
        call Condensate(dofs*acc_n0, acc_stiff, acc_force, acc_timeforce)
      end if
    else if (nb > 0) then
      if (acc_transient .and. .not. acc_solver % globalbubbles) then
        acc_xl(1:2*nd-1:2) = acc_lke(1:nd)
        acc_xl(2:2*nd:2) = acc_ldis(1:nd)
        acc_xlprev(1:2*nd-1:2) = acc_pke(1:nd)
        acc_xlprev(2:2*nd:2) = acc_pdis(1:nd)
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
  end subroutine ko_bulk_elem

  subroutine ko_bulk_tail() bind(c, name="ko_bulk_tail")
    call DefaultFinishBulkAssembly()
    call Info('KOmega', 'Assembly done', level=4)
  end subroutine ko_bulk_tail

  ! DefaultFinishAssembly runs before the Dirichlet wall loop in the
  ! original (UpdateDirichletDof patches the assembled matrix), so it
  ! is its own step here.
  subroutine ko_finish_assembly() bind(c, name="ko_finish_assembly")
    call DefaultFinishAssembly()
  end subroutine ko_finish_assembly

  ! Dirichlet boundary loop (OmegaWall where requested).
  subroutine ko_boundary() bind(c, name="ko_boundary")
    integer :: t, n
    logical :: gotit
    do t = 1, acc_solver % mesh % numberofboundaryelements
      acc_elem => GetBoundaryElement(t)
      if (.not. ActiveBoundaryElement()) cycle
      n = GetElementNOFNodes()
      acc_bc => GetBC()
      if (.not. associated(acc_bc)) cycle
      if (GetLogical(acc_bc, 'Omega Wall BC', gotit) .or. &
          GetLogical(acc_bc, 'Noslip Wall BC', gotit)) &
          call ko_omegawall(n)
    end do
  end subroutine ko_boundary

  subroutine ko_dirichlet() bind(c, name="ko_dirichlet")
    call DefaultDirichletBCs()
    call Info('KOmega', 'Set boundaries done', level=4)
  end subroutine ko_dirichlet

  ! Solve tail: DefaultSolve through the convergence check. The
  ! KMax/EMax MAXVAL reductions are pure (their results are never
  ! used) and skipped; the 1.0d-12 clip loop is transcribed exactly.
  ! Returns 1 when NonlinConverged == 1 (exit the nonlinear loop).
  function ko_solve_tail() result(exitloop) bind(c, name="ko_solve_tail")
    integer(c_int) :: exitloop
    integer :: i, k
    real(kind=dp) :: kval, eval
    type(variable_t), pointer :: ke
    acc_norm = DefaultSolve()
    ke => acc_solver % variable
    do i = 1, size(ke % perm)
      k = ke % perm(i)
      if (k <= 0) cycle
      kval = ke % values(2*k-1)
      eval = ke % values(2*k-0)
      ke % values(2*k-1) = max(kval, 1.0d-12)
      ke % values(2*k-0) = max(eval, 1.0d-12)
    end do
    acc_relchange = ke % nonlinchange
    write(message, *) 'Result Norm   : ', acc_norm
    call Info('KOmega', message, level=4)
    write(message, *) 'Relative Change : ', acc_relchange
    call Info('KOmega', message, level=4)
    if (ke % nonlinconverged == 1) then
      exitloop = 1
    else
      exitloop = 0
    end if
  end function ko_solve_tail

  ! The LocalMatrix internal procedure: element local matrices and
  ! right-hand side for the diffusion-convection equations. Transcribed
  ! statement for statement; host data (Material, Bubbles, local
  ! solutions) arrives via module state. Like the original it USEs
  ! MaterialModels.
  subroutine ko_localmatrix(n, nd)
    use materialmodels
    integer :: n, nd
    real(kind=dp) :: ddBasisddx(nd,3,3)
    real(kind=dp) :: Basis(nd)
    real(kind=dp) :: dBasisdx(nd,3), detJ
    real(kind=dp) :: UX(n), UY(n), UZ(n), Velo(3), dVelodx(3,3), Energy(n), &
        Dissipation(n), Distance(n), Density(n), Viscosity(n)
    real(kind=dp) :: A(2,2), M(2,2), Prod, div, ProdTensor(3,3)
    integer :: i, j, c, p, q, t, dim, NBasis
    real(kind=dp) :: LoadatIp(2), Cmu, Rho, mu, Tmu, Effmu(2)
    real(kind=dp) :: s, u, v, w, K, Omega, Strain(3,3), Vorticity(3,3), dist
    real(kind=dp) :: StrainMeasure, VorticityMeasure, X, Y, Z, SigmaK, &
        SigmaO, Beta, CD, F1, F2, F3, F4, rGamma, GradK(3), GradO(3)
    real(kind=dp) :: Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3), SqrtMetric
    logical :: stat
    type(gaussintegrationpoints_t), target :: IntegStuff
    dim = CoordinateSystemDimension()
    Viscosity(1:n) = GetReal(acc_material, 'Viscosity')
    Density(1:n) = GetReal(acc_material, 'Density')
    call GetScalarLocalSolution(UX, 'Velocity 1')
    call GetScalarLocalSolution(UY, 'Velocity 2')
    call GetScalarLocalSolution(UZ, 'Velocity 3')
    call GetScalarLocalSolution(Energy, 'Kinetic energy')
    call GetScalarLocalSolution(Dissipation, 'Kinetic Dissipation')
    acc_force = 0.0D0
    acc_stiff = 0.0D0
    acc_mass = 0.0D0
    NBasis = nd
    if (acc_bubbles) then
      IntegStuff = GaussPoints(acc_elem, acc_elem % type % GaussPoints2)
    else
      IntegStuff = GaussPoints(acc_elem)
    end if
    do t = 1, IntegStuff % n
      u = IntegStuff % u(t)
      v = IntegStuff % v(t)
      w = IntegStuff % w(t)
      stat = ElementInfo(acc_elem, acc_enodes, u, v, w, detJ, &
          Basis, dBasisdx, Bubbles=acc_bubbles)
      s = detJ * IntegStuff % s(t)
      if (CurrentCoordinateSystem() /= Cartesian) then
        X = sum(acc_enodes % x(1:n)*Basis(1:n))
        Y = sum(acc_enodes % y(1:n)*Basis(1:n))
        Z = sum(acc_enodes % z(1:n)*Basis(1:n))
        call CoordinateSystemInfo(Metric, SqrtMetric, Symb, dSymb, X, Y, Z)
        s = s * SqrtMetric
      end if
      Velo = 0.0_dp
      Velo(1) = sum(UX(1:n)*Basis(1:n))
      Velo(2) = sum(UY(1:n)*Basis(1:n))
      Velo(3) = sum(UZ(1:n)*Basis(1:n))
      dVelodx = 0.0_dp
      do i = 1, dim
        dVelodx(1,i) = sum(UX(1:n)*dBasisdx(1:n,i))
        dVelodx(2,i) = sum(UY(1:n)*dBasisdx(1:n,i))
        dVelodx(3,i) = sum(UZ(1:n)*dBasisdx(1:n,i))
      end do
      if (CurrentCoordinateSystem() == Cartesian) then
        Strain = 0.5_dp * (dVelodx + transpose(dVelodx))
        StrainMeasure = max(sqrt(2*sum(Strain*Strain)), 1.0d-10)
        Vorticity = 0.5_dp * (dVelodx - transpose(dVelodx))
        VorticityMeasure = sqrt(2*sum(Vorticity*Vorticity))
      else
        StrainMeasure = sqrt(SecondInvariant(Velo, dVelodx, Metric, Symb)/2)
      end if
      K = sum(Energy(1:n)*Basis(1:n))
      Omega = sum(Dissipation(1:n)*Basis(1:n))
      mu = sum(Viscosity(1:n)*Basis(1:n))
      rho = sum(Density(1:n)*Basis(1:n))
      Beta = 0.075_dp
      SigmaK = 2.000_dp
      SigmaO = 2.000_dp
      rGamma = 5._dp/9._dp
      Tmu = rho*K/max(Omega, 1.0d-10)
      Effmu(1) = mu + Tmu/SigmaK
      Effmu(2) = mu + Tmu/SigmaO
      Prod = 2*Tmu*sum(Strain*dVelodx)
      do p = 1, NBasis
        do q = 1, NBasis
          M = 0.0d0
          A = 0.0d0
          M(1,1) = rho * Basis(q) * Basis(p)
          M(2,2) = rho * Basis(q) * Basis(p)
          A(1,1) = A(1,1) + rho * 0.09_dp * Omega * Basis(q) * Basis(p)
          A(2,2) = A(2,2) + rho * Beta * Omega * Basis(q) * Basis(p)
          if (CurrentCoordinateSystem() == Cartesian) then
            do i = 1, dim
              A(1,1) = A(1,1) + Effmu(1) * dBasisdx(q,i) * dBasisdx(p,i)
              A(2,2) = A(2,2) + Effmu(2) * dBasisdx(q,i) * dBasisdx(p,i)
            end do
          else
            do i = 1, dim
              do j = 1, dim
                A(1,1) = A(1,1) + Metric(i,j) * Effmu(1) * &
                    dBasisdx(q,i) * dBasisdx(p,j)
                A(2,2) = A(2,2) + Metric(i,j) * Effmu(2) * &
                    dBasisdx(q,i) * dBasisdx(p,j)
              end do
            end do
          end if
          do i = 1, dim
            A(1,1) = A(1,1) + rho * Velo(i) * dBasisdx(q,i) * Basis(p)
            A(2,2) = A(2,2) + rho * Velo(i) * dBasisdx(q,i) * Basis(p)
          end do
          do i = 1, 2
            do j = 1, 2
              acc_stiff(2*(p-1)+i,2*(q-1)+j) = acc_stiff(2*(p-1)+i,2*(q-1)+j) + s*A(i,j)
              acc_mass(2*(p-1)+i,2*(q-1)+j) = acc_mass(2*(p-1)+i,2*(q-1)+j) + s*M(i,j)
            end do
          end do
        end do
      end do
      LoadAtIP(1) = Prod
      LoadAtIP(2) = rGamma * Prod * Omega / K
      do p = 1, NBasis
        acc_force(2*(p-1)+1) = acc_force(2*(p-1)+1) + s*LoadAtIp(1)*Basis(p)
        acc_force(2*(p-1)+2) = acc_force(2*(p-1)+2) + s*LoadAtIp(2)*Basis(p)
      end do
    end do
  end subroutine ko_localmatrix

  ! Wall law for the k-Omega turbulence model. Transcribed statement
  ! for statement; the boundary element arrives via module state.
  subroutine ko_omegawall(n)
    integer :: n
    real(kind=dp) :: omega_wall, dist, mu(32), rho(32), x0(n), y0(n), z0(n), x, y, z
    integer :: i, j, np
    type(element_t), pointer :: parent
    parent => acc_elem % boundaryinfo % left
    if (.not. associated(parent)) &
        parent => acc_elem % boundaryinfo % right
    if (.not. associated(parent)) return
    np = GetElementNOFNodes(parent)
    rho(1:np) = GetReal(GetMaterial(parent), 'Density', UElement=parent)
    mu(1:np) = GetReal(GetMaterial(parent), 'Viscosity', UElement=parent)
    x0(1:n) = acc_model % nodes % x(acc_elem % nodeindexes)
    y0(1:n) = acc_model % nodes % y(acc_elem % nodeindexes)
    z0(1:n) = acc_model % nodes % z(acc_elem % nodeindexes)
    omega_wall = 1.d10
    do i = 1, np
      j = parent % nodeindexes(i)
      if (any(j == acc_elem % nodeindexes(1:n))) cycle
      x = acc_model % nodes % x(j)
      y = acc_model % nodes % y(j)
      z = acc_model % nodes % z(j)
      dist = minval((x-x0(1:n))**2 + (y-y0(1:n))**2 + (z-z0(1:n))**2)
      if (dist < AEPS) cycle
      omega_wall = 6*mu(i)/rho(i)/0.075_dp/dist
      j = 2*acc_solver % variable % perm(j)
      call UpdateDirichletDof(acc_solver % matrix, j, omega_wall)
    end do
  end subroutine ko_omegawall

  ! -- KOmegaLegacy_Init body -------------------------------------------
  subroutine ko_init_begin(h_solver, transient) bind(c, name="ko_init_begin")
    type(c_ptr), value :: h_solver
    integer(c_int), value :: transient
    call c_f_pointer(h_solver, acc_solver)
    acc_transient = (transient /= 0)
  end subroutine ko_init_begin

  subroutine ko_init_body() bind(c, name="ko_init_body")
    type(valuelist_t), pointer :: solverparams
    logical :: found, pbubble, legacybubbles
    character(len=MAX_NAME_LEN) :: str
    solverparams => GetSolverParams()
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
      if (.not. found) legacybubbles = .true.
      if (legacybubbles) then
        call ListAddNewInteger(solverparams, 'Nonlinear System Min Iterations', 2)
        call ListAddNewInteger(solverparams, 'Nonlinear System Max Iterations', 2)
      end if
    end if
  end subroutine ko_init_body

end module komegalegacy_acc
