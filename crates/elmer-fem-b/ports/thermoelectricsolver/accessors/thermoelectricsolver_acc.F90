! ThermoElectric solver-entry port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/ThermoElectricSolver.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (TYPE(C_PTR) handles,
! VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All Fortran-object
! traffic — TYPE(...) handles, OPTIONAL/POINTER dummies, assumed-shape
! descriptors, hidden CHARACTER lengths, array-pointer List results — stays
! inside Fortran, so the Rust side is pure C calls. The one routine that
! takes a CHARACTER dummy (InputTensor's Name) is replicated once per call
! site with the name spelled literally, selected by an integer tag, so no
! hidden length ever crosses the boundary.
!
! Rust drives element/boundary/nonlinear-loop control flow and owns the
! SAVE bubble-history arrays (bx/bxprev) with identical lifetimes.
module tea_acc
  use iso_c_binding
  use types
  use defutils
  use integration
  use matrixassembly
  use eleminfo
  use elementgeometry
  use messages
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(model_t), pointer, save :: acc_model => null()
  type(element_t), pointer, save :: acc_elem => null()

contains

  ! -- entry state ----------------------------------------------------
  subroutine tea_begin(h_model, h_solver) bind(c, name="tea_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine tea_begin

  ! -- Init entry body (the whole BLOCK of ThermoElectricSolver_Init) --
  subroutine tea_init_fixup(transient) bind(c, name="tea_init_fixup")
    integer(c_int), value :: transient
    type(valuelist_t), pointer :: params
    character(:), allocatable :: elementstr
    logical :: pbubble, found, gotit

    params => GetSolverParams()

    elementstr = ListGetString(params, 'Element', found)
    pbubble = found .and. index(elementstr, 'b:') > 0

    if (pbubble .and. transient /= 0 .and. &
        .not. ListGetLogical(params, 'Bubbles in Global System', gotit)) then
      call ListAddNewInteger(params, 'Nonlinear System Min Iterations', 2)
      call ListAddNewInteger(params, 'Nonlinear System Max Iterations', 2)
    end if
  end subroutine tea_init_fixup

  ! -- mesh sizing (bubble-history SAVE sizing) ------------------------
  function tea_mesh_maxbdofs() result(n) bind(c, name="tea_mesh_maxbdofs")
    integer(c_int) :: n
    n = acc_solver % mesh % maxbdofs
  end function tea_mesh_maxbdofs

  function tea_mesh_bulk() result(n) bind(c, name="tea_mesh_bulk")
    integer(c_int) :: n
    n = acc_solver % mesh % numberofbulkelements
  end function tea_mesh_bulk

  function tea_mesh_bnd() result(n) bind(c, name="tea_mesh_bnd")
    integer(c_int) :: n
    n = acc_solver % mesh % numberofboundaryelements
  end function tea_mesh_bnd

  function tea_get_timestep() result(n) bind(c, name="tea_get_timestep")
    integer(c_int) :: n
    n = GetTimestep()
  end function tea_get_timestep

  function tea_stabilize() result(r) bind(c, name="tea_stabilize")
    integer(c_int) :: r
    type(valuelist_t), pointer :: params
    params => GetSolverParams()
    r = merge(1, 0, GetStabilizeFlag(params))
  end function tea_stabilize

  ! -- nonlinear-loop parameters (single call, original read order) ---
  subroutine tea_nl_params(ni, ntol, mi, mtol, mi_found, mtol_found) &
      bind(c, name="tea_nl_params")
    integer(c_int) :: ni, mi, mi_found, mtol_found
    real(c_double) :: ntol, mtol
    type(valuelist_t), pointer :: params
    logical :: f
    params => GetSolverParams()
    ni = GetInteger(params, 'Nonlinear System Newton After Iterations', f)
    ntol = GetConstReal(params, 'Nonlinear System Newton After Tolerance', f)
    mi = GetInteger(params, 'Nonlinear System Max Iterations', f)
    mi_found = merge(1, 0, f)
    mtol = GetConstReal(params, 'Nonlinear System Convergence Tolerance', f)
    mtol_found = merge(1, 0, f)
  end subroutine tea_nl_params

  ! -- assembly driver -------------------------------------------------
  subroutine tea_default_initialize() bind(c, name="tea_default_initialize")
    call DefaultInitialize()
  end subroutine tea_default_initialize

  function tea_nof_active() result(n) bind(c, name="tea_nof_active")
    integer(c_int) :: n
    n = GetNOFActive()
  end function tea_nof_active

  subroutine tea_default_finish_bulk() bind(c, name="tea_default_finish_bulk")
    call DefaultFinishBulkAssembly()
  end subroutine tea_default_finish_bulk

  function tea_nof_bnd() result(n) bind(c, name="tea_nof_bnd")
    integer(c_int) :: n
    n = GetNOFBoundaryElements()
  end function tea_nof_bnd

  subroutine tea_default_finish() bind(c, name="tea_default_finish")
    call DefaultFinishAssembly()
  end subroutine tea_default_finish

  subroutine tea_default_dirichlet() bind(c, name="tea_default_dirichlet")
    call DefaultDirichletBCs()
  end subroutine tea_default_dirichlet

  function tea_default_solve() result(norm) bind(c, name="tea_default_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function tea_default_solve

  function tea_nonlin_change() result(r) bind(c, name="tea_nonlin_change")
    real(c_double) :: r
    r = acc_solver % variable % nonlinchange
  end function tea_nonlin_change

  subroutine tea_report(norm, rel) bind(c, name="tea_report")
    real(c_double), value :: norm, rel
    write(message, *) 'Result Norm   : ', norm
    call Info('ThermoElectricSolver', message, level=4)
    write(message, *) 'Relative Change : ', rel
    call Info('ThermoElectricSolver', message, level=4)
  end subroutine tea_report

  subroutine tea_fatal_alloc() bind(c, name="tea_fatal_alloc")
    call Fatal('ThermoElectricSolver', 'Memory allocation error.')
  end subroutine tea_fatal_alloc

  ! -- active-element context ------------------------------------------
  subroutine tea_begin_active(t, n, nd, nb) bind(c, name="tea_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, nb
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes(acc_elem)
    nb = GetElementNOFBDOFs(acc_elem)
    nd = GetElementNOFDOFs(acc_elem) + nb
  end subroutine tea_begin_active

  function tea_elem_index() result(k) bind(c, name="tea_elem_index")
    integer(c_int) :: k
    k = acc_elem % elementindex
  end function tea_elem_index

  function tea_global_bubbles() result(r) bind(c, name="tea_global_bubbles")
    integer(c_int) :: r
    r = merge(1, 0, acc_solver % globalbubbles)
  end function tea_global_bubbles

  ! -- boundary-element context ----------------------------------------
  subroutine tea_begin_bnd(t, n, nd, active) bind(c, name="tea_begin_bnd")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, active
    acc_elem => GetBoundaryElement(t)
    active = merge(1, 0, ActiveBoundaryElement(acc_elem))
    n = GetElementNOFNodes(acc_elem)
    nd = GetElementNOFDOFs(acc_elem)
  end subroutine tea_begin_bnd

  ! -- InputTensor replica (CHARACTER name stays in Fortran) -----------
  subroutine InputTensorCopy(Material, Tensor, IsScalar, which, Element)
    type(valuelist_t), pointer :: Material
    real(kind=dp) :: Tensor(:, :, :)
    logical :: IsScalar
    integer :: which
    type(element_t) :: Element
    integer :: i, j, n
    logical :: stat
    real(kind=dp), pointer :: Hwrk(:, :, :)
    nullify(Hwrk)
    Tensor = 0.0_dp
    IsScalar = .true.
    n = Element % TYPE % NumberOfNodes
    select case (which)
    case (1)
      call ListGetRealArray(Material, 'Heat Conductivity', Hwrk, n, &
          Element % NodeIndexes, stat)
    case (2)
      call ListGetRealArray(Material, 'Seebeck Coefficient', Hwrk, n, &
          Element % NodeIndexes, stat)
    case (3)
      call ListGetRealArray(Material, 'Electric Conductivity', Hwrk, n, &
          Element % NodeIndexes, stat)
    case (4)
      call ListGetRealArray(Material, 'Electric Permittivity', Hwrk, n, &
          Element % NodeIndexes, stat)
    case (5)
      call ListGetRealArray(Material, 'Electric Transfer Coefficient', Hwrk, n, &
          Element % NodeIndexes, stat)
    end select
    if (.not. stat) return
    IsScalar = size(Hwrk, 1) == 1 .and. size(Hwrk, 2) == 1
    if (IsScalar) then
      do i = 1, size(Tensor, 1)
        Tensor(i, i, 1:n) = Hwrk(1, 1, 1:n)
      end do
    else
      if (size(Hwrk, 1) == 1) then
        do i = 1, min(3, size(Hwrk, 2))
          Tensor(i, i, 1:n) = Hwrk(1, i, 1:n)
        end do
      else if (size(Hwrk, 2) == 1) then
        do i = 1, min(3, size(Hwrk, 1))
          Tensor(i, i, 1:n) = Hwrk(i, 1, 1:n)
        end do
      else
        do i = 1, min(3, size(Hwrk, 1))
          do j = 1, min(3, size(Hwrk, 2))
            Tensor(i, j, 1:n) = Hwrk(i, j, 1:n)
          end do
        end do
      end if
    end if
  end subroutine InputTensorCopy

  ! -- LocalMatrix core: exact transcription of the internal procedure
  ! -- body through the end of the IP loop (fills MASS/STIFF/FORCE/SOL).
  subroutine AssembleCore(n, nd, MASS, STIFF, FORCE, SOL, stab, newton)
    integer, intent(in) :: n, nd
    real(kind=dp), target :: MASS(2*nd, 2*nd), STIFF(2*nd, 2*nd), FORCE(2*nd)
    real(kind=dp) :: SOL(2, nd)
    logical, intent(in) :: stab, newton
    real(kind=dp), pointer :: A(:, :), M(:, :)
    real(kind=dp) :: Basis(nd), dBasisdx(nd, 3), DetJ
    logical :: Stat, isScalar, Found
    integer :: i, j, p, q, t
    type(GaussIntegrationPoints_t) :: IP
    real(kind=dp) :: EF(3), TG(3), CD(3), JH, HS, IC, SB(3, 3), s, &
        PeltSB(3, 3), PeltSi(3, 3), &
        sigma(3, 3), alpha(3, 3), pelt(3, 3), hcond(3, 3), &
        epsil(3, 3), rho, c_p, Temp
    real(kind=dp) :: v_eff(3), StreamVec(nd), hK, mK, VNorm, Tau, kappa_eff
    real(kind=dp) :: sigma_n(3, 3, n), alpha_n(3, 3, n), pelt_n(3, 3, n), &
        hcond_n(3, 3, n), epsil_n(3, 3, n), rho_n(n), c_p_n(n), epsil0
    real(kind=dp) :: LOAD(2, nd)
    type(valuelist_t), pointer :: Material, BF
    type(nodes_t) :: Nodes

    call GetElementNodes(Nodes, acc_elem)

    Material => GetMaterial(acc_elem)

    call InputTensorCopy(Material, hcond_n, isScalar, 1, acc_elem)
    call InputTensorCopy(Material, alpha_n, isScalar, 2, acc_elem)
    call InputTensorCopy(Material, sigma_n, isScalar, 3, acc_elem)

    if (ListCheckPresent(Material, 'Relative Permittivity')) then
      epsil0 = GetConstReal(acc_model % Simulation, 'Permittivity of Vacuum', Found)
      if (.not. Found) then
        call Fatal('ThermoElectricSolver', &
            '> Relative Permittivity < requires > Permittivity of Vacuum < !')
      end if
      call InputTensorCopy(Material, epsil_n, isScalar, 4, acc_elem)
    else
      epsil0 = 1.0_dp
      call InputTensorCopy(Material, epsil_n, isScalar, 4, acc_elem)
    end if

    rho_n = GetReal(Material, 'Density', Found, acc_elem)
    c_p_n = GetReal(Material, 'Heat Capacity', Found, acc_elem)

    BF => GetBodyForce(acc_elem)
    if (associated(BF)) then
      LOAD(1, 1:n) = GetReal(BF, 'Heat Source', Found, acc_elem)
      if (Found) then
        LOAD(1, 1:n) = rho_n * LOAD(1, 1:n)
      else
        LOAD(1, 1:n) = GetReal(BF, 'Volumetric Heat Source', Found, acc_elem)
      end if
      LOAD(2, 1:n) = GetReal(BF, 'Current Source', Found, acc_elem)
    end if

    call GetVectorLocalSolution(SOL, UElement=acc_elem)

    MASS = 0.0_dp
    STIFF = 0.0_dp
    FORCE = 0.0_dp

    !Numerical integration:
    !----------------------
    IP = GaussPoints(acc_elem)
    do t = 1, IP % n
      stat = ElementInfo(acc_elem, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx)

      do i = 1, 3
        do j = 1, 3
          hcond(i, j) = SUM(hcond_n(i, j, :) * Basis(1:n))
          alpha(i, j) = SUM(alpha_n(i, j, :) * Basis(1:n))
          sigma(i, j) = SUM(sigma_n(i, j, :) * Basis(1:n))
          epsil(i, j) = epsil0 * SUM(epsil_n(i, j, :) * Basis(1:n))
        end do
      end do
      rho = SUM(rho_n * Basis(1:n))
      c_p = SUM(c_p_n * Basis(1:n))

      IC = SUM(LOAD(2, 1:n) * Basis(1:n))    ! input current
      HS = SUM(LOAD(1, 1:n) * Basis(1:n))    ! input heat source

      EF = -MATMUL(SOL(2, 1:nd), dBasisdx(1:nd, :)) ! electric field from prev.sol.
      TG = MATMUL(SOL(1, 1:nd), dBasisdx(1:nd, :)) ! temperature gradient  ""

      CD = MATMUL(sigma, EF - MATMUL(alpha, TG)) ! current density from prev.sol
      JH = SUM(CD * EF)                        ! joule heating       ""

      SB = MATMUL(sigma, alpha)

      Temp = SUM(Basis(1:nd) * SOL(1, 1:nd)) ! temperature from previous iteration
      Pelt = Temp * alpha                   ! Peltier coefficient
      PeltSB = MATMUL(Pelt, SB)
      PeltSi = MATMUL(Pelt, Sigma)

      ! SUPG stabilization of the temperature equation: the Peltier term
      ! div(T*alpha*J) linearizes into an advective "drift" alpha*J acting on
      ! grad(T), which can dominate over thermal diffusion at high current
      ! density -- treated exactly as an externally given (lagged, Picard)
      ! convection field here, the same way HeatSolveVec treats its velocity.
      ! ---------------------------------------------------------------------
      if (stab) then
        v_eff = MATMUL(alpha, CD)
        VNorm = SQRT(SUM(v_eff**2))

        if (VNorm > 0.0_dp) then
          kappa_eff = SUM(v_eff * MATMUL(hcond, v_eff)) / VNorm**2
        else
          kappa_eff = 0.0_dp
        end if

        hK = acc_elem % hK
        mK = acc_elem % StabilizationMK

        if (VNorm > 0.0_dp .and. kappa_eff /= 0.0_dp) then
          Tau = MIN(1.0_dp, mK * hK * VNorm / (2.0_dp * ABS(kappa_eff)))
          Tau = hK * Tau / (2.0_dp * VNorm)
        else
          Tau = 0.0_dp
        end if

        do p = 1, nd
          StreamVec(p) = SUM(v_eff * dBasisdx(p, :))
        end do
      else
        Tau = 0.0_dp
        StreamVec = 0.0_dp
      end if

      s = IP % s(t) * DetJ

      do p = 1, nd
        i = 2 * (p - 1)
        do q = 1, nd
          j = 2 * (q - 1)
          M => MASS(i + 1:i + 2, j + 1:j + 2)
          A => STIFF(i + 1:i + 2, j + 1:j + 2)

          ! thermal damping
          ! ---------------
          M(1, 1) = M(1, 1) + s * rho * c_p * Basis(q) * Basis(p)
          M(1, 1) = M(1, 1) + s * Tau * rho * c_p * Basis(q) * StreamVec(p)  ! SUPG

          ! thermal diffusion
          ! -----------------
          A(1, 1) = A(1, 1) + s * SUM(MATMUL(hcond, dBasisdx(q, :)) * dBasisdx(p, :))

          ! SUPG streamline-diffusion consistency term for the Peltier drift
          ! ------------------------------------------------------------------
          A(1, 1) = A(1, 1) + s * Tau * StreamVec(q) * StreamVec(p)

          if (newton) then

            ! Newton linarization of: div(pelt*J) (= div(T*alpha*J))
            ! ======================================================

            ! Peltier coefficient implicitly
            ! -------------------------------
            A(1, 1) = A(1, 1) - s * Basis(q) * SUM(MATMUL(alpha, CD) * dBasisdx(p, :))

            ! temperature gradient part of div(pelt_0*J)
            ! ------------------------------------------
            A(1, 1) = A(1, 1) - s * SUM(MATMUL(PeltSB, dBasisdx(q, :)) * dBasisdx(p, :))

            ! electric field part of div(pelt_0*J)
            ! ------------------------------------
            A(1, 2) = A(1, 2) - s * SUM(MATMUL(PeltSi, dBasisdx(q, :)) * dBasisdx(p, :))

            ! Newton linarization of Joule heating=-(J,E) ~
            ! -(J_0,E) + (J,E_0) - (J_0,E_0) (<-- to rhs)
            ! =============================================
            ! -(J_0,E)
            ! --------
            A(1, 2) = A(1, 2) + s * SUM(CD * dBasisdx(q, :)) * Basis(p)

            ! temperature gradient part of (J,E_0)
            ! ------------------------------------
            A(1, 1) = A(1, 1) - s * SUM(MATMUL(SB, dBasisdx(q, :)) * EF) * Basis(p)

            ! electric field part part of (J,E_0)
            ! ------------------------------------
            A(1, 2) = A(1, 2) - s * SUM(MATMUL(sigma, dBasisdx(q, :)) * EF) * Basis(p)
          end if

          ! displacement current
          ! --------------------
          M(2, 2) = M(2, 2) + &
              s * SUM(MATMUL(epsil, dBasisdx(q, :)) * dBasisdx(p, :))

          ! Seebeck
          ! -------
          A(2, 1) = A(2, 1) + &
              s * SUM(MATMUL(SB, dBasisdx(q, :)) * dBasisdx(p, :))

          ! electric diffusion
          ! ------------------
          A(2, 2) = A(2, 2) + &
              s * SUM(MATMUL(sigma, dBasisdx(q, :)) * dBasisdx(p, :))
        end do

        ! heat load
        ! ---------
        FORCE(i + 1) = FORCE(i + 1) + s * (HS + JH) * Basis(p)
        FORCE(i + 1) = FORCE(i + 1) + s * Tau * (HS + JH) * StreamVec(p)  ! SUPG

        ! div(pelt_0*J_0)
        ! ---------------
        FORCE(i + 1) = FORCE(i + 1) + &
            s * SUM(MATMUL(Pelt, CD) * dBasisdx(p, :))

        ! electric load
        ! -------------
        FORCE(i + 2) = FORCE(i + 2) + s * IC * Basis(p)
      end do
    end do
  end subroutine AssembleCore

  ! -- LocalMatrix, steady / condensed-bubble tail ----------------------
  subroutine tea_lm_std(n, nd, nb, stab, newton, transient) &
      bind(c, name="tea_lm_std")
    integer(c_int), value :: n, nd, nb, stab, newton, transient
    real(kind=dp) :: MASS(2*nd, 2*nd), STIFF(2*nd, 2*nd), FORCE(2*nd)
    real(kind=dp) :: SOL(2, nd)

    call AssembleCore(n, nd, MASS, STIFF, FORCE, SOL, &
        stab /= 0, newton /= 0)

    if (transient /= 0) call Default1stOrderTime(MASS, STIFF, FORCE, UElement=acc_elem)
    call CondensateP(2 * (nd - nb), 2 * nb, STIFF, FORCE)

    call DefaultUpdateEquations(STIFF, FORCE, UElement=acc_elem)
  end subroutine tea_lm_std

  ! -- LocalMatrix, transient bubble-condensation tail -------------------
  subroutine tea_lm_bubble(n, nd, nb, stab, newton, dt, bprev, bcur) &
      bind(c, name="tea_lm_bubble")
    integer(c_int), value :: n, nd, nb, stab, newton
    real(c_double), value :: dt
    real(c_double) :: bprev(2*nb), bcur(2*nb)
    real(kind=dp) :: MASS(2*nd, 2*nd), STIFF(2*nd, 2*nd), FORCE(2*nd)
    real(kind=dp) :: SOL(2, nd), SOLprev(2, nd)
    real(kind=dp) :: xvec(2*(nd - nb)), xprevvec(2*(nd - nb))
    integer :: Ncond, pp

    call AssembleCore(n, nd, MASS, STIFF, FORCE, SOL, &
        stab /= 0, newton /= 0)

    Ncond = nd - nb
    call GetVectorLocalSolution(SOLprev, UElement=acc_elem, tStep=-1)

    do pp = 1, Ncond
      xvec(2*(pp - 1) + 1) = SOL(1, pp)
      xvec(2*(pp - 1) + 2) = SOL(2, pp)
      xprevvec(2*(pp - 1) + 1) = SOLprev(1, pp)
      xprevvec(2*(pp - 1) + 2) = SOLprev(2, pp)
    end do

    call CondensatePTransient(Ncond, nb, 2, dt, MASS, STIFF, FORCE, &
        xprevvec, xvec, bprev(1:2*nb), bcur(1:2*nb))

    call DefaultUpdateEquations(STIFF, FORCE, UElement=acc_elem)
  end subroutine tea_lm_bubble

  ! -- BoundaryLocalMatrix: exact transcription --------------------------
  subroutine tea_bnd_assemble(n, nd, transient) bind(c, name="tea_bnd_assemble")
    integer(c_int), value :: n, nd, transient
    real(kind=dp), target :: MASS(2*nd, 2*nd), STIFF(2*nd, 2*nd), FORCE(2*nd)
    real(kind=dp), pointer :: A(:, :)
    real(kind=dp) :: Basis(nd), dBasisdx(nd, 3), s, DetJ
    logical :: Stat, Found, IsScalar, EFVFound, AnisoElecBC
    integer :: i, j, p, q, t
    type(GaussIntegrationPoints_t) :: IP
    real(kind=dp) :: HF_n(n), ET_n(n), HT_n(n), EF_n(n), HF, ET, HT, EF
    real(kind=dp) :: EFV_n(3, n), EFV(3), Normal(3), &
        EP_n(n), EP, ETC_n(3, 3, n), ETC(3, 3), HTe
    type(valuelist_t), pointer :: BC
    type(nodes_t) :: Nodes

    BC => GetBC(acc_elem)
    if (.not. associated(BC)) return

    call GetElementNodes(Nodes, acc_elem)

    HF_n = GetReal(BC, 'Heat Flux', Found, acc_elem)
    ET_n = GetReal(BC, 'External Temperature', Found, acc_elem)
    HT_n = GetReal(BC, 'Heat Transfer Coefficient', Found, acc_elem)

    EF_n = GetReal(BC, 'Electric Flux', Found, acc_elem)

    ! Anisotropic current density vector, projected onto the boundary normal
    ! ------------------------------------------------------------------------
    call GetRealVector(BC, EFV_n, 'Electric Flux Vector', EFVFound, acc_elem)

    ! Anisotropic Robin BC: surface conductance tensor + external potential
    ! ------------------------------------------------------------------------
    EP_n = GetReal(BC, 'External Potential', Found, acc_elem)
    AnisoElecBC = ListCheckPresent(BC, 'Electric Transfer Coefficient')
    if (AnisoElecBC) then
      call InputTensorCopy(BC, ETC_n, IsScalar, 5, acc_elem)
    else
      ETC_n = 0.0_dp
    end if

    MASS = 0.0_dp
    STIFF = 0.0_dp
    FORCE = 0.0_dp

    !Numerical integration:
    !----------------------
    IP = GaussPoints(acc_elem)
    do t = 1, IP % n
      stat = ElementInfo(acc_elem, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx)

      s = IP % s(t) * DetJ

      HF = SUM(HF_n * Basis(1:n)) ! heat flux
      ET = SUM(ET_n * Basis(1:n)) ! external temperature
      HT = SUM(HT_n * Basis(1:n)) ! heat transfer coefficient

      EF = SUM(EF_n * Basis(1:n)) ! electric flux
      EP = SUM(EP_n * Basis(1:n)) ! external potential

      HTe = 0.0_dp
      if (EFVFound .or. AnisoElecBC) then
        Normal = NormalVector(acc_elem, Nodes, IP % U(t), IP % V(t), .true.)
      end if

      if (EFVFound) then
        do i = 1, 3
          EFV(i) = SUM(EFV_n(i, 1:n) * Basis(1:n))
        end do
        EF = EF + SUM(EFV * Normal) ! current density vector dotted with the normal
      end if

      if (AnisoElecBC) then
        do i = 1, 3
          do j = 1, 3
            ETC(i, j) = SUM(ETC_n(i, j, :) * Basis(1:n))
          end do
        end do
        HTe = SUM(Normal * MATMUL(ETC, Normal)) ! normal-normal component of the conductance tensor
      end if

      do p = 1, nd
        i = 2 * (p - 1)
        do q = 1, nd
          j = 2 * (q - 1)
          A => STIFF(i + 1:i + 2, j + 1:j + 2)
          A(1, 1) = A(1, 1) + s * HT * Basis(q) * Basis(p)
          A(2, 2) = A(2, 2) + s * HTe * Basis(q) * Basis(p)
        end do
        FORCE(i + 1) = FORCE(i + 1) + s * (HF + HT * ET) * Basis(p)
        FORCE(i + 2) = FORCE(i + 2) + s * (EF + HTe * EP) * Basis(p)
      end do
    end do

    if (transient /= 0) then
      call Default1stOrderTime(MASS, STIFF, FORCE, UElement=acc_elem)
    end if
    call DefaultUpdateEquations(STIFF, FORCE, UElement=acc_elem)
  end subroutine tea_bnd_assemble

end module tea_acc
