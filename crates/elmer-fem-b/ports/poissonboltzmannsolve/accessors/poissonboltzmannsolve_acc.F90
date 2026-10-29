! T021-B PoissonBoltzmann port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/PoissonBoltzmannSolve.F90, compiled
! against the built .mod files of the target configure (A2 shipped
! mechanism, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls. Order-sensitive
! reductions (SUM), transcendental evaluations (SINH/COSH/EXP),
! whole-array updates, formatted I/O (WRITE/Info) and the Elmer library
! calls all stay in Fortran, so bit-identity holds by construction
! wherever the call sequence matches.
!
! SAVE replication: the module SAVE state below mirrors the original's
! SAVE list (LocalStiffMatrix, Load, LocalForce, LocalPot, ElementNodes,
! CalculateCharge, CalculateEnergy, ElectricField1/2/3, ChargeField,
! Permittivity, NameField/NameEnergy/NameCharge, CalculateField,
! ConstantWeights, NewtonLinearization) with identical lifetimes, plus
! the host-associated scalars the CONTAINS procedures share
! (NonDimensional, SymmetricCharges, Beta/Alpha, PermittivityOfVacuum,
! Cunit, NewtonIter/Tol, TotEnergy, EquationName, DIM). Non-SAVEd host
! locals (FluxBC, bf_id, loop indices) are wrapper locals, recomputed on
! every call exactly as in the original. NewtonLinearization keeps its
! declaration-init (.FALSE.) latch: set TRUE on the same conditions,
! never reset per call. Rust owns only the AllocationsDone flag and the
! DO-loop counter, mirroring the original's first-call and EXIT shape.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module poissonboltzmannsolve_acc
  use iso_c_binding
  use types
  use lists
  use integration
  use elementdescription
  use differentials
  use solverbasics
  use elementutils
  use adaptive
  use defutils
  use coordinatesystems
  use messages
  implicit none

  type(model_t), pointer, save :: s_model => null()
  type(solver_t), pointer, save :: s_solver => null()
  type(matrix_t), pointer, save :: s_stiffm => null()
  type(element_t), pointer, save :: s_curelem => null()
  type(valuelist_t), pointer, save :: s_params => null()

  ! First-mesh SAVE state (original SAVE list, first-call sizing).
  type(nodes_t), save :: s_elemnodes
  real(kind=dp), allocatable, save :: s_perm(:), s_localforce(:), &
      s_localstiff(:,:), s_load(:), s_localpot(:)
  real(kind=dp), pointer, save :: s_ef1(:) => null(), s_ef2(:) => null(), &
      s_ef3(:) => null(), s_charge(:) => null(), s_energy(:) => null()
  real(kind=dp), pointer, save :: s_pot(:) => null()
  integer, pointer, save :: s_potperm(:) => null(), &
      s_fieldperm(:) => null()
  character(len=MAX_NAME_LEN), save :: s_namefield, s_nameenergy, &
      s_namecharge
  character(len=MAX_NAME_LEN), save :: s_equation
  logical, save :: s_calccharge, s_calcenergy, s_calcfield, s_constweights
  logical, save :: s_newton_lin = .FALSE.

  ! Host-associated scalars shared with the CONTAINS procedures.
  logical, save :: s_nondim, s_symm
  real(kind=dp), save :: s_permofvac, s_cunit
  real(kind=dp), save :: s_betapos, s_betaneg, s_alphapos, s_alphaneg
  real(kind=dp), save :: s_newtontol, s_norm
  integer, save :: s_newtoniter, s_dim
  real(kind=dp), save :: s_totenergy

contains

  ! -- entry state ----------------------------------------------------
  subroutine pbb_begin(h_model, h_solver) bind(c, name="pbb_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, s_model)
    call c_f_pointer(h_solver, s_solver)
  end subroutine pbb_begin

  function pbb_local_nodes() result(n) bind(c, name="pbb_local_nodes")
    integer(c_int) :: n
    s_pot => s_solver % variable % values
    s_potperm => s_solver % variable % perm
    n = COUNT(s_potperm > 0)
  end function pbb_local_nodes

  subroutine pbb_pre(nonlin, newtontol) bind(c, name="pbb_pre")
    integer(c_int) :: nonlin
    real(c_double) :: newtontol
    logical :: found, gotit
    real(kind=dp) :: cboltz, ccharge
    s_dim = CoordinateSystemDimension()
    s_params => GetSolverParams()
    s_nondim = .TRUE.
    nonlin = GetInteger(s_params, 'Nonlinear System Max Iterations', found)
    IF (.NOT. found) nonlin = 1
    s_newtontol = GetConstReal(s_params, &
        'Nonlinear System Newton After Tolerance', found)
    s_newtoniter = GetInteger(s_params, &
        'Nonlinear System Newton After Iterations', found)
    IF (s_newtoniter == 0) s_newton_lin = .TRUE.
    s_permofvac = ListGetConstReal(s_model % constants, &
        'Permittivity Of Vacuum', gotit)
    IF (.NOT. gotit) s_permofvac = 1.0d0
    cboltz = ListGetConstReal(s_model % constants, 'Boltzmann Constant')
    ccharge = ListGetConstReal(s_model % constants, 'Unit Charge')
    s_cunit = ccharge / cboltz
    s_equation = ListGetString(s_solver % values, 'Equation')
    newtontol = s_newtontol
  end subroutine pbb_pre

  ! -- first-call permanent storage (original AllocationsDone block) --
  subroutine pbb_first_alloc() bind(c, name="pbb_first_alloc")
    real(kind=dp), pointer :: field(:)
    integer :: n, i, istat
    logical :: gotit
    n = s_model % maxelementnodes
    ALLOCATE(s_elemnodes % x(n), s_elemnodes % y(n), s_elemnodes % z(n), &
        s_perm(n), s_localforce(n), s_localstiff(n, n), s_load(n), &
        s_localpot(n), STAT=istat)
    IF (istat /= 0) THEN
      CALL Fatal('StatElecSolve', 'Memory allocation error 1')
    END IF
    s_calcfield = ListGetLogical(s_solver % values, &
        'Calculate Electric Field', gotit)
    IF (.NOT. gotit) s_calcfield = .TRUE.
    IF (s_calcfield) THEN
      ALLOCATE(s_ef1(s_model % numberofnodes), &
          s_ef2(s_model % numberofnodes), &
          s_ef3(s_model % numberofnodes), STAT=istat)
      IF (istat /= 0) THEN
        CALL Fatal('StatElecSolve', 'Memory allocation error 2')
      END IF
      WRITE (s_namefield, '(A,A)') TRIM(ComponentName(s_solver % variable)), &
          ' Field'
    END IF
    s_calccharge = ListGetLogical(s_solver % values, &
        'Calculate Electric Charge', gotit)
    IF (.NOT. gotit) s_calccharge = .TRUE.
    IF (s_calccharge) THEN
      ALLOCATE(s_charge(s_model % numberofnodes), STAT=istat)
      WRITE (s_namecharge, '(A,A)') TRIM(ComponentName(s_solver % variable)), &
          ' Charge'
      IF (istat /= 0) THEN
        CALL Fatal('StatElecSolve', 'Memory allocation error 3')
      END IF
    END IF
    DO i = 1, s_model % numberofequations
      s_calcenergy = ListGetLogical(s_model % equations(i) % values, &
          'Calculate Electric Energy', gotit)
      IF (gotit) EXIT
    END DO
    IF (.NOT. gotit) s_calcenergy = ListGetLogical(s_solver % values, &
        'Calculate Electric Energy', gotit)
    IF (s_calcenergy) THEN
      ALLOCATE(s_energy(s_model % numberofnodes), STAT=istat)
      WRITE (s_nameenergy, '(A,A)') TRIM(ComponentName(s_solver % variable)), &
          ' Energy'
      IF (istat /= 0) THEN
        CALL Fatal('StatElecSolve', 'Memory allocation error 4')
      END IF
    END IF
    s_constweights = ListGetLogical(s_solver % values, &
        'Constant Weights', gotit)
    s_stiffm => s_solver % matrix
    IF (.NOT. ASSOCIATED(s_stiffm % massvalues)) THEN
      ALLOCATE(s_stiffm % massvalues(s_model % numberofnodes))
      s_stiffm % massvalues = 0.0d0
    END IF
    IF (s_calcfield) THEN
      field => s_ef1
      CALL VariableAdd(s_solver % mesh % variables, s_solver % mesh, &
          s_solver, TRIM(s_namefield) // ' 1', 1, field, s_potperm)
      field => s_ef2
      CALL VariableAdd(s_solver % mesh % variables, s_solver % mesh, &
          s_solver, TRIM(s_namefield) // ' 2', 1, field, s_potperm)
      IF (s_dim == 3) THEN
        field => s_ef3
        CALL VariableAdd(s_solver % mesh % variables, s_solver % mesh, &
            s_solver, TRIM(s_namefield) // ' 3', 1, field, s_potperm)
      END IF
    END IF
    IF (s_calccharge) THEN
      field => s_charge
      CALL VariableAdd(s_solver % mesh % variables, s_solver % mesh, &
          s_solver, s_namecharge, 1, field, s_potperm)
    END IF
    IF (s_calcenergy) THEN
      field => s_energy
      CALL VariableAdd(s_solver % mesh % variables, s_solver % mesh, &
          s_solver, s_nameenergy, 1, field, s_potperm)
    END IF
  end subroutine pbb_first_alloc

  ! -- per-call derived-variable refresh -------------------------------
  subroutine pbb_refresh() bind(c, name="pbb_refresh")
    type(variable_t), pointer :: var
    IF (s_calcfield) THEN
      var => VariableGet(s_model % variables, TRIM(s_namefield) // ' 1')
      s_ef1 => var % values
      var => VariableGet(s_model % variables, TRIM(s_namefield) // ' 2')
      s_ef2 => var % values
      s_ef2 = 0.0d0
      s_ef1 = 0.0d0
      IF (s_dim == 3) THEN
        var => VariableGet(s_model % variables, TRIM(s_namefield) // ' 3')
        s_ef3 => var % values
        s_ef3 = 0.0d0
      END IF
    END IF
    IF (s_calccharge) THEN
      var => VariableGet(s_model % variables, TRIM(s_namecharge))
      s_charge => var % values
      s_fieldperm => var % perm
      s_charge = 0.0d0
    END IF
    IF (s_calcenergy) THEN
      var => VariableGet(s_model % variables, TRIM(s_nameenergy))
      s_energy => var % values
      s_fieldperm => var % perm
      s_energy = 0.0d0
    END IF
  end subroutine pbb_refresh

  ! -- nonlinear-iteration driver --------------------------------------
  subroutine pbb_iter_info(iter) bind(c, name="pbb_iter_info")
    integer(c_int), value :: iter
    integer :: it
    it = iter
    CALL Info('PoissonBoltzmannSolve', &
        '-------------------------------------', Level=4)
    WRITE (Message, *) 'Potential iteration', it
    CALL Info('PoissonBoltzmannSolve', Message, Level=4)
    CALL Info('PoissonBoltzmannSolve', &
        '-------------------------------------', Level=4)
    CALL Info('PoissonBoltzmannSolve', 'Starting Assembly...', Level=4)
  end subroutine pbb_iter_info

  subroutine pbb_newton_update(iter) bind(c, name="pbb_newton_update")
    integer(c_int), value :: iter
    IF (iter > s_newtoniter) s_newton_lin = .TRUE.
  end subroutine pbb_newton_update

  subroutine pbb_default_initialize() bind(c, name="pbb_default_initialize")
    CALL DefaultInitialize()
  end subroutine pbb_default_initialize

  subroutine pbb_bulk_assemble() bind(c, name="pbb_bulk_assemble")
    integer :: t, n, bf_id
    integer, pointer :: nodeindexes(:)
    logical :: gotit
    real(kind=dp) :: at0
    at0 = RealTime()
    DO t = 1, s_solver % numberofactiveelements
      IF (RealTime() - at0 > 1.0) THEN
        WRITE (Message, '(a,i3,a)') '   Assembly: ', INT(100.0 - 100.0 * &
            (s_solver % mesh % numberofbulkelements - t) / &
            (1.0 * s_solver % mesh % numberofbulkelements)), ' % done'
        CALL Info('PoissonBoltzmannSolve', Message, Level=5)
        at0 = RealTime()
      END IF
      s_curelem => s_solver % mesh % elements(s_solver % activeelements(t))
      s_model % currentelement => s_curelem
      nodeindexes => s_curelem % nodeindexes
      n = s_curelem % TYPE % numberofnodes
      s_elemnodes % x(1:n) = s_solver % mesh % nodes % x(nodeindexes)
      s_elemnodes % y(1:n) = s_solver % mesh % nodes % y(nodeindexes)
      s_elemnodes % z(1:n) = s_solver % mesh % nodes % z(nodeindexes)
      s_localpot(1:n) = s_pot(s_potperm(nodeindexes(1:n)))
      bf_id = ListGetInteger( &
          s_model % bodies(s_curelem % bodyid) % values, &
          'Body Force', gotit, minv=1, maxv=s_model % numberofbodyforces)
      IF (gotit) THEN
        s_load(1:n) = ListGetReal(s_model % bodyforces(bf_id) % values, &
            'Charge Density', n, nodeindexes, gotit)
      ELSE
        s_load(1:n) = 0.0d0
      END IF
      CALL pbb_material_params()
      CALL pbb_compose(s_localstiff, s_localforce, s_perm, s_load, &
          s_curelem, n, s_elemnodes)
      CALL DefaultUpdateEquations(s_localstiff, s_localforce)
    END DO
    CALL DefaultFinishBulkAssembly()
  end subroutine pbb_bulk_assemble

  subroutine pbb_neumann() bind(c, name="pbb_neumann")
    integer :: t, i, n
    integer, pointer :: nodeindexes(:)
    logical :: gotit, fluxbc
    DO t = s_solver % mesh % numberofbulkelements + 1, &
        s_solver % mesh % numberofbulkelements + &
        s_solver % mesh % numberofboundaryelements
      s_curelem => s_solver % mesh % elements(t)
      s_model % currentelement => s_curelem
      DO i = 1, s_model % numberofbcs
        IF (s_curelem % boundaryinfo % constraint == &
            s_model % bcs(i) % tag) THEN
          fluxbc = ListGetLogical(s_model % bcs(i) % values, &
              'Electric Flux BC', gotit)
          IF (gotit .AND. .NOT. fluxbc) CYCLE
          s_model % currentelement => s_solver % mesh % elements(t)
          n = s_curelem % TYPE % numberofnodes
          nodeindexes => s_curelem % nodeindexes
          IF (ANY(s_potperm(nodeindexes) <= 0)) CYCLE
          s_load(1:n) = ListGetReal(s_model % bcs(i) % values, &
              'Surface Charge', n, nodeindexes, gotit)
          IF (.NOT. gotit) CYCLE
          s_elemnodes % x(1:n) = s_solver % mesh % nodes % x(nodeindexes)
          s_elemnodes % y(1:n) = s_solver % mesh % nodes % y(nodeindexes)
          s_elemnodes % z(1:n) = s_solver % mesh % nodes % z(nodeindexes)
          CALL pbb_statbc(s_localstiff, s_localforce, s_load, &
              s_curelem, n, s_elemnodes)
          CALL DefaultUpdateEquations(s_localstiff, s_localforce)
        END IF
      END DO
    END DO
  end subroutine pbb_neumann

  subroutine pbb_finish() bind(c, name="pbb_finish")
    CALL DefaultFinishAssembly()
  end subroutine pbb_finish

  subroutine pbb_dirichlet() bind(c, name="pbb_dirichlet")
    CALL DefaultDirichletBCs()
  end subroutine pbb_dirichlet

  ! -- timing + reporting (WRITE/Info stay in Fortran) -------------------
  function pbb_cpu() result(t) bind(c, name="pbb_cpu")
    real(c_double) :: t
    t = CPUTime()
  end function pbb_cpu

  subroutine pbb_info_assembly(at) bind(c, name="pbb_info_assembly")
    real(c_double), value :: at
    real(kind=dp) :: a
    a = at
    WRITE (Message, *) 'Assembly (s)          :', a
    CALL Info('PoissonBoltzmannSolve', Message, Level=4)
  end subroutine pbb_info_assembly

  subroutine pbb_info_solve(st) bind(c, name="pbb_info_solve")
    real(c_double), value :: st
    real(kind=dp) :: a
    a = st
    WRITE (Message, *) 'Solve (s)             :', a
    CALL Info('PoissonBoltzmannSolve', Message, Level=4)
  end subroutine pbb_info_solve

  subroutine pbb_solve_core() bind(c, name="pbb_solve_core")
    s_norm = DefaultSolve()
  end subroutine pbb_solve_core

  subroutine pbb_info_normblock(change, conv) bind(c, name="pbb_info_normblock")
    real(c_double) :: change
    integer(c_int) :: conv
    WRITE (Message, *) 'Result Norm   : ', s_norm
    CALL Info('PoissonBoltzmannSolve', Message, Level=4)
    s_solver % variable % norm = s_norm
    WRITE (Message, *) 'Relative Change : ', s_solver % variable % nonlinchange
    CALL Info('PoissonBoltzmannSolve', Message, Level=4)
    IF (s_solver % variable % nonlinchange < s_newtontol) &
        s_newton_lin = .TRUE.
    change = s_solver % variable % nonlinchange
    conv = merge(1, 0, s_solver % variable % nonlinconverged == 1)
  end subroutine pbb_info_normblock

  ! -- post-solve flux + energy ------------------------------------------
  subroutine pbb_flux_tail() bind(c, name="pbb_flux_tail")
    IF (s_calcfield .OR. s_calccharge .OR. s_calcenergy) THEN
      CALL pbb_flux(s_model, s_pot, s_potperm)
    END IF
    IF (s_calcenergy) THEN
      WRITE (Message, *) 'Tot. Electric Energy  :', s_totenergy
      CALL Info('PoissonBoltzmannSolve', Message, Level=4)
      CALL ListAddConstReal(s_model % simulation, &
          'RES: Electric Energy', s_totenergy)
    END IF
  end subroutine pbb_flux_tail

  ! -- transcribed CONTAINS procedures (1:1, host state via module) -------
  subroutine pbb_material_params()
    integer :: k, zpos, zneg, nloc
    logical :: gotit
    real(kind=dp) :: reft, npos, nneg, ucharge
    k = ListGetInteger( &
        s_model % bodies(s_curelem % bodyid) % values, &
        'Material', minv=1, maxv=s_model % numberofmaterials)
    s_betapos = ListGetConstReal(s_model % materials(k) % values, &
        'Poisson Boltzmann Beta', gotit)
    IF (.NOT. gotit) s_nondim = .FALSE.
    s_alphapos = ListGetConstReal(s_model % materials(k) % values, &
        'Poisson Boltzmann Alpha', gotit)
    IF (.NOT. gotit) s_nondim = .FALSE.
    s_symm = .TRUE.
    IF (.NOT. s_nondim) THEN
      nloc = s_curelem % TYPE % numberofnodes
      s_perm(1:nloc) = ListGetReal(s_model % materials(k) % values, &
          'Relative Permittivity', nloc, s_curelem % nodeindexes)
      zpos = ListGetInteger(s_model % materials(k) % values, &
          'Charge Number', gotit)
      IF (.NOT. gotit) THEN
        zpos = ListGetInteger(s_model % materials(k) % values, &
            'Positive Charge Number')
        zneg = ListGetInteger(s_model % materials(k) % values, &
            'Negative Charge Number')
        zneg = -ABS(zneg)
        s_symm = .FALSE.
      ELSE
        zneg = -zpos
      END IF
      reft = ListGetConstReal(s_model % materials(k) % values, &
          'Reference Temperature')
      s_betapos = s_cunit * zpos / reft
      npos = ListGetConstReal(s_model % materials(k) % values, &
          'Ion Density', gotit)
      IF (.NOT. gotit) THEN
        npos = ListGetConstReal(s_model % materials(k) % values, &
            'Positive Ion Density')
        nneg = ListGetConstReal(s_model % materials(k) % values, &
            'Negative Ion Density')
        s_symm = .FALSE.
      ELSE
        nneg = npos
      END IF
      ucharge = ListGetConstReal(s_model % constants, 'Unit Charge')
      s_alphapos = 2.0d0 * npos * zpos * ucharge
      IF (.NOT. s_symm) THEN
        s_betapos = s_cunit * zpos / reft
        s_betaneg = s_cunit * zneg / reft
        s_alphapos = npos * zpos * ucharge
        s_alphaneg = nneg * zneg * ucharge
      END IF
    END IF
  end subroutine pbb_material_params

  subroutine pbb_compose(stiffmatrix, force, permittivity, load, element, &
      n, nodes)
    real(kind=dp) :: stiffmatrix(:, :), force(:), load(:), permittivity(:)
    integer :: n
    type(nodes_t) :: nodes
    type(element_t), TARGET :: element
    real(kind=dp) :: sqrtmetric, metric(3, 3), symb(3, 3, 3), &
        dsymb(3, 3, 3, 3)
    real(kind=dp) :: basis(n), dbasisdx(n, 3)
    real(kind=dp) :: sqrtel, u, v, w, s, a, b, l, c, x, y, z
    real(kind=dp) :: csinh, ccosh, pot
    logical :: stat
    integer :: i, p, q, t, dim
    type(gaussintegrationpoints_t) :: integstuff
    dim = CoordinateSystemDimension()
    force = 0.0d0
    stiffmatrix = 0.0d0
    integstuff = GaussPoints(element)
    DO t = 1, integstuff % n
      u = integstuff % u(t)
      v = integstuff % v(t)
      w = integstuff % w(t)
      s = integstuff % s(t)
      stat = ElementInfo(element, nodes, u, v, w, sqrtel, basis, dbasisdx)
      IF (CurrentCoordinateSystem() /= Cartesian) THEN
        x = SUM(s_elemnodes % x(1:n) * basis(1:n))
        y = SUM(s_elemnodes % y(1:n) * basis(1:n))
        z = SUM(s_elemnodes % z(1:n) * basis(1:n))
      END IF
      CALL CoordinateSystemInfo(metric, sqrtmetric, symb, dsymb, x, y, z)
      s = s * sqrtel * sqrtmetric
      IF (s_nondim) THEN
        c = 1.0_dp
      ELSE
        c = s_permofvac * SUM(permittivity(1:n) * basis(1:n))
      END IF
      pot = SUM(s_localpot(1:n) * basis(1:n))
      l = SUM(load(1:n) * basis(1:n))
      IF (s_symm) THEN
        csinh = SINH(s_betapos * pot)
        ccosh = COSH(s_betapos * pot)
        l = l - s_alphapos * csinh
        IF (s_newton_lin) THEN
          l = l + s_alphapos * s_betapos * ccosh * pot
        END IF
      ELSE
        l = l + s_alphapos * EXP(-s_betapos * pot) + &
            s_alphaneg * EXP(-s_betaneg * pot)
        IF (s_newton_lin) THEN
          l = l + s_alphapos * s_betapos * EXP(-s_betapos * pot) * pot + &
              s_alphaneg * s_betaneg * EXP(-s_betaneg * pot) * pot
        END IF
      END IF
      DO p = 1, n
        DO q = 1, n
          a = c * SUM(dbasisdx(p, 1:dim) * dbasisdx(q, 1:dim))
          IF (s_newton_lin) THEN
            IF (s_symm) THEN
              a = a + s_alphapos * s_betapos * ccosh * basis(p) * basis(q)
            ELSE
              a = a + basis(p) * basis(q) * &
                  (s_alphapos * s_betapos * EXP(-s_betapos * pot) + &
                  s_alphaneg * s_betaneg * EXP(-s_betaneg * pot))
            END IF
          END IF
          stiffmatrix(p, q) = stiffmatrix(p, q) + s * a
        END DO
        force(p) = force(p) + s * l * basis(p)
      END DO
    END DO
  end subroutine pbb_compose

  subroutine pbb_statbc(boundarymatrix, boundaryvector, loadvector, &
      element, n, nodes)
    real(kind=dp) :: boundarymatrix(:, :), boundaryvector(:), &
        loadvector(:)
    type(nodes_t) :: nodes
    type(element_t) :: element
    integer :: n
    real(kind=dp) :: basis(n)
    real(kind=dp) :: dbasisdx(n, 3), sqrtel
    real(kind=dp) :: sqrtmetric, metric(3, 3), symb(3, 3, 3), &
        dsymb(3, 3, 3, 3)
    real(kind=dp) :: u, v, w, s, x, y, z
    real(kind=dp) :: forceload
    real(kind=dp), pointer :: u_integ(:), v_integ(:), w_integ(:), s_integ(:)
    integer :: t, q, n_integ
    type(gaussintegrationpoints_t), TARGET :: integstuff
    logical :: stat
    boundaryvector = 0.0d0
    boundarymatrix = 0.0d0
    integstuff = GaussPoints(element)
    u_integ => integstuff % u
    v_integ => integstuff % v
    w_integ => integstuff % w
    s_integ => integstuff % s
    n_integ = integstuff % n
    DO t = 1, n_integ
      u = u_integ(t)
      v = v_integ(t)
      w = w_integ(t)
      stat = ElementInfo(element, nodes, u, v, w, sqrtel, basis, dbasisdx)
      IF (CurrentCoordinateSystem() /= Cartesian) THEN
        x = SUM(s_elemnodes % x(1:n) * basis(1:n))
        y = SUM(s_elemnodes % y(1:n) * basis(1:n))
        z = SUM(s_elemnodes % z(1:n) * basis(1:n))
      END IF
      CALL CoordinateSystemInfo(metric, sqrtmetric, symb, dsymb, x, y, z)
      s = s_integ(t) * sqrtel * sqrtmetric
      forceload = SUM(loadvector(1:n) * basis)
      DO q = 1, n
        boundaryvector(q) = boundaryvector(q) + s * basis(q) * forceload
      END DO
    END DO
  end subroutine pbb_statbc

  subroutine pbb_flux(model, potential, reorder)
    type(model_t) :: model
    real(kind=dp) :: potential(:)
    integer :: reorder(:)
    type(element_t), pointer :: element
    type(nodes_t) :: nodes
    type(gaussintegrationpoints_t), TARGET :: integstuff
    real(kind=dp), pointer :: u_integ(:), v_integ(:), w_integ(:), s_integ(:)
    real(kind=dp), allocatable :: sumofweights(:)
    real(kind=dp) :: permofvac
    real(kind=dp) :: basis(model % maxelementnodes)
    real(kind=dp) :: dbasisdx(model % maxelementnodes, 3)
    real(kind=dp) :: sqrtel
    real(kind=dp) :: elementpot(model % maxelementnodes)
    real(kind=dp) :: energydensity
    real(kind=dp) :: flux(3), field(3), elemvol
    real(kind=dp) :: s, ug, vg, wg, grad(3), epsgrad(3)
    real(kind=dp) :: sqrtmetric, metric(3, 3), symb(3, 3, 3), &
        dsymb(3, 3, 3, 3)
    real(kind=dp) :: x, y, z, l
    integer, pointer :: nodeindexes(:)
    integer :: n, n_integ, t, tg, i, j, k, dim, bf_id
    logical :: stat, gotit
    ALLOCATE (nodes % x(model % maxelementnodes))
    ALLOCATE (nodes % y(model % maxelementnodes))
    ALLOCATE (nodes % z(model % maxelementnodes))
    ALLOCATE (sumofweights(model % numberofnodes))
    sumofweights = 0.0d0
    s_totenergy = 0.0d0
    permofvac = ListGetConstReal(model % constants, &
        'Permittivity Of Vacuum', gotit)
    IF (.NOT. gotit) permofvac = 1
    dim = CoordinateSystemDimension()
    DO t = 1, model % numberofbulkelements
      element => model % elements(t)
      s_curelem => element
      nodeindexes => element % nodeindexes
      IF (.NOT. CheckElementEquation(model, element, s_equation)) CYCLE
      n = element % TYPE % numberofnodes
      IF (ANY(reorder(nodeindexes) == 0)) CYCLE
      elementpot(1:n) = potential(reorder(nodeindexes(1:n)))
      bf_id = ListGetInteger(model % bodies(s_curelem % bodyid) % values, &
          'Body Force', gotit, minv=1, maxv=model % numberofbodyforces)
      IF (gotit) THEN
        s_load(1:n) = ListGetReal(model % bodyforces(bf_id) % values, &
            'Charge Density', n, nodeindexes, gotit)
      ELSE
        s_load(1:n) = 0.0d0
      END IF
      CALL pbb_material_params()
      nodes % x(1:n) = model % nodes % x(nodeindexes)
      nodes % y(1:n) = model % nodes % y(nodeindexes)
      nodes % z(1:n) = model % nodes % z(nodeindexes)
      integstuff = GaussPoints(element)
      u_integ => integstuff % u
      v_integ => integstuff % v
      w_integ => integstuff % w
      s_integ => integstuff % s
      n_integ = integstuff % n
      k = ListGetInteger(model % bodies(element % bodyid) % values, &
          'Material', minv=1, maxv=model % numberofmaterials)
      IF (.NOT. s_nondim) THEN
        s_perm(1:n) = ListGetReal(model % materials(k) % values, &
            'Relative Permittivity', n, nodeindexes)
      END IF
      energydensity = 0.0d0
      flux = 0.0d0
      field = 0.0d0
      elemvol = 0.0d0
      DO tg = 1, n_integ
        ug = u_integ(tg)
        vg = v_integ(tg)
        wg = w_integ(tg)
        stat = ElementInfo(element, nodes, ug, vg, wg, sqrtel, basis, &
            dbasisdx)
        s = 1
        IF (CurrentCoordinateSystem() /= Cartesian) THEN
          x = SUM(nodes % x(1:n) * basis(1:n))
          y = SUM(nodes % y(1:n) * basis(1:n))
          z = SUM(nodes % z(1:n) * basis(1:n))
          s = 2 * PI
        END IF
        CALL CoordinateSystemInfo(metric, sqrtmetric, symb, dsymb, x, y, z)
        s = s * sqrtmetric * sqrtel * s_integ(tg)
        epsgrad = 0.0d0
        DO j = 1, dim
          grad(j) = SUM(dbasisdx(1:n, j) * elementpot(1:n))
        END DO
        epsgrad(1:dim) = SUM(s_perm(1:n) * basis(1:n)) * grad(1:dim)
        energydensity = energydensity + s * SUM(grad(1:dim) * epsgrad(1:dim))
        DO j = 1, dim
          field(j) = field(j) - grad(j) * s
        END DO
        elemvol = elemvol + s
      END DO
      IF (s_calcenergy) THEN
        s_totenergy = s_totenergy + energydensity
      END IF
      IF (s_constweights) THEN
        energydensity = energydensity / elemvol
        field(1:dim) = field(1:dim) / elemvol
        sumofweights(reorder(nodeindexes(1:n))) = &
            sumofweights(reorder(nodeindexes(1:n))) + 1
      ELSE
        sumofweights(reorder(nodeindexes(1:n))) = &
            sumofweights(reorder(nodeindexes(1:n))) + elemvol
      END IF
      IF (s_calcenergy) THEN
        s_energy(s_fieldperm(nodeindexes(1:n))) = &
            s_energy(s_fieldperm(nodeindexes(1:n))) + energydensity
      END IF
      IF (s_calcfield) THEN
        s_ef1(s_fieldperm(nodeindexes(1:n))) = &
            s_ef1(s_fieldperm(nodeindexes(1:n))) + field(1)
        s_ef2(s_fieldperm(nodeindexes(1:n))) = &
            s_ef2(s_fieldperm(nodeindexes(1:n))) + field(2)
        IF (dim == 3) THEN
          s_ef3(s_fieldperm(nodeindexes(1:n))) = &
              s_ef3(s_fieldperm(nodeindexes(1:n))) + field(3)
        END IF
      END IF
      IF (s_calccharge) THEN
        s_charge(s_fieldperm(nodeindexes(1:n))) = s_load(1:n)
        DO j = 1, n
          IF (s_symm) THEN
            l = -s_alphapos * SINH(s_betapos * elementpot(j))
          ELSE
            l = s_alphapos * EXP(-s_betapos * elementpot(j)) + &
                s_alphaneg * EXP(-s_betaneg * elementpot(j))
          END IF
          s_charge(s_fieldperm(nodeindexes(j))) = s_load(j) + l
        END DO
      END IF
    END DO
    DO i = 1, model % numberofnodes
      IF (reorder(i) == 0) CYCLE
      IF (ABS(sumofweights(reorder(i))) > AEPS) THEN
        IF (s_calcenergy) s_energy(s_fieldperm(i)) = &
            s_energy(s_fieldperm(i)) / sumofweights(reorder(i))
        IF (s_calcfield) THEN
          s_ef1(s_fieldperm(i)) = s_ef1(s_fieldperm(i)) / &
              sumofweights(reorder(i))
          s_ef2(s_fieldperm(i)) = s_ef2(s_fieldperm(i)) / &
              sumofweights(reorder(i))
          IF (dim == 3) THEN
            s_ef3(s_fieldperm(i)) = s_ef3(s_fieldperm(i)) / &
                sumofweights(reorder(i))
          END IF
        END IF
      END IF
    END DO
    s_totenergy = s_permofvac * s_totenergy / 2.0d0
    DEALLOCATE (nodes % x, nodes % y, nodes % z, sumofweights)
  end subroutine pbb_flux

end module poissonboltzmannsolve_acc
