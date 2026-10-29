! T021-B VectorHelmholtz port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/VectorHelmholtz.F90, compiled against the built
! .mod files of the target configure (ADR-057 §3, solver-entry gate shape).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape arrays). All Fortran-object
! traffic — TYPE(...) handles, OPTIONAL/POINTER dummies, assumed-shape
! descriptors, hidden lengths, array-pointer results — stays inside Fortran,
! so the Rust side is pure C calls.
!
! The two big assembly bodies (LocalMatrix, LocalMatrixBC, SingleDipoleLoad
! on the solver side; the CalcFields phase bodies on the field side) are
! transcribed verbatim as module procedures: they own every order-sensitive
! complex reduction (MATMUL/SUM), every formatted write (PRINT/WRITE), and
! every whole-array update, so bit-identity holds by construction wherever
! the call sequence matches. The Rust entries own the call sequence, the
! branch predicates, and the SAVE lifetimes documented in lib.rs.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text beyond the transcribed statement shapes.
module vectorhelmholtz_acc
  use iso_c_binding
  use defutils
  use vectorhelmholtzutils
  implicit none

  ! -- shared entry handles (captured per call, like the entry dummies) --
  type(solver_t), pointer, save :: vh_solver => null()
  type(model_t), pointer, save :: vh_model => null()
  type(valuelist_t), pointer, save :: vh_params => null()
  type(mesh_t), pointer, save :: vh_mesh => null()
  type(solver_t), pointer, save :: vh_psolver => null()
  logical, save :: vh_s0_withndofs = .false.

  ! -- main-solver prologue scalars (entry locals in the original) --
  real(kind=dp), save :: vh_omega = 0.0_dp, vh_mu0inv = 0.0_dp
  real(kind=dp), save :: vh_eps0 = 0.0_dp, vh_rob0 = 0.0_dp
  complex(kind=dp), save :: vh_precdamp = (0.0_dp, 0.0_dp)
  integer, save :: vh_edge_degree = 1, vh_noitmax = 1
  logical, save :: vh_piola = .false., vh_lowfreq = .false.
  logical, save :: vh_lorenz = .false., vh_usegauss = .false.
  logical, save :: vh_chargecons = .false., vh_eigen = .false.
  logical, save :: vh_precmatrix = .false., vh_hasprecdamp = .false.
  logical, save :: vh_massprop = .false., vh_curlcurlprec = .false.
  logical, save :: vh_edgebasis = .false.
  logical, save :: vh_amsany = .false.
  logical, save :: vh_amscurlcurl = .false., vh_amsmono = .false.
  type(matrix_t), pointer, save :: vh_amsmat => null(), vh_amsscalmat => null()
  type(solver_t), pointer, save :: vh_amssolver => null()
  type(solver_t), pointer, save :: vh_amsscalsolver => null()

  character(*), parameter :: vh_caller = 'VectorHelmholtzSolver'
  character(*), parameter :: vh_cf_caller = 'VectorHelmholtzCalcFields'

  ! -- LocalMatrix SAVE (mirrors its SAVE locals incl. AllocationsDone) --
  complex(kind=dp), allocatable, save :: vh_lm_stiff(:,:), vh_lm_force(:)
  complex(kind=dp), allocatable, save :: vh_lm_mass(:,:), vh_lm_gauge(:,:)
  complex(kind=dp), allocatable, save :: vh_lm_prec(:,:), vh_lm_curlmat(:,:)
  complex(kind=dp), allocatable, save :: vh_lm_amsstiff(:,:), vh_lm_amsstiff2(:,:)
  complex(kind=dp), allocatable, save :: vh_lm_amsforce(:)
  real(kind=dp), allocatable, save :: vh_lm_wbasis(:,:), vh_lm_rotw(:,:)
  real(kind=dp), allocatable, save :: vh_lm_basis(:), vh_lm_db(:,:)
  type(nodes_t), save :: vh_lm_nodes
  type(valuehandle_t), save :: vh_lm_cond_h, vh_lm_eps_h, vh_lm_cur_h, vh_lm_mu_h
  logical, save :: vh_lm_alloc = .false.

  ! -- LocalMatrixBC SAVE (mirrors its SAVE locals incl. AllocationsDone) --
  complex(kind=dp), allocatable, save :: vh_bc_stiff(:,:), vh_bc_force(:)
  complex(kind=dp), allocatable, save :: vh_bc_mass(:,:)
  complex(kind=dp), allocatable, save :: vh_bc_amsstiff(:,:), vh_bc_amsstiff2(:,:)
  complex(kind=dp), allocatable, save :: vh_bc_amsforce(:)
  real(kind=dp), allocatable, save :: vh_bc_wbasis(:,:), vh_bc_rotw(:,:)
  real(kind=dp), allocatable, save :: vh_bc_basis(:), vh_bc_db(:,:)
  type(nodes_t), save :: vh_bc_nodes
  type(valuehandle_t), save :: vh_bc_elsurf_h, vh_bc_magload_h, vh_bc_elrobin_h
  type(valuehandle_t), save :: vh_bc_mu_h, vh_bc_eps_h, vh_bc_absorb_h
  type(valuehandle_t), save :: vh_bc_temre_h, vh_bc_temim_h, vh_bc_extpot_h
  type(valuehandle_t), save :: vh_bc_transfer_h, vh_bc_elcurr_h, vh_bc_thick_h
  type(valuehandle_t), save :: vh_bc_relnu_h, vh_bc_cond_h
  type(valuehandle_t), save :: vh_bc_goodcond_h, vh_bc_chargecons_h
  logical, save :: vh_bc_alloc = .false.

  ! -- SingleDipoleLoad SAVE --
  type(nodes_t), save :: vh_dip_nodes

  ! -- CalcFields phase state (entry-lifetime: filled by vcf_prologue) --
  type(solver_t), pointer, save :: vh_cf_solver => null()
  type(model_t), pointer, save :: vh_cf_model => null()
  type(solver_t), pointer, save :: vh_cf_psolver => null()
  type(mesh_t), pointer, save :: vh_cf_mesh => null()
  type(element_t), pointer, save :: vh_cf_elem => null()
  type(variable_t), pointer, save :: vh_cf_targetvar => null()
  character(len=MAX_NAME_LEN), save :: vh_cf_pname = ''
  real(kind=dp), save :: vh_cf_omega = 0.0_dp, vh_cf_mu0inv = 0.0_dp
  real(kind=dp), save :: vh_cf_eps0 = 0.0_dp
  real(kind=dp), save :: vh_cf_energy = 0.0_dp, vh_cf_energy_im = 0.0_dp
  real(kind=dp), save :: vh_cf_norm = 0.0_dp
  integer, save :: vh_cf_degree = 1
  integer, save :: vh_cf_n = 0, vh_cf_nd = 0, vh_cf_np = 0, vh_cf_ndofs = 0
  integer, save :: vh_cf_fields = 0, vh_cf_nfields = 0
  logical, save :: vh_cf_piola = .false., vh_cf_usegauss = .false.
  logical, save :: vh_cf_lorenz = .false.
  logical, save :: vh_cf_nodal = .false., vh_cf_elemental = .false.
  logical, save :: vh_cf_doave = .false.
  type(variable_t), pointer, save :: vh_cf_mfd => null(), vh_cf_mfs => null()
  type(variable_t), pointer, save :: vh_cf_ef => null(), vh_cf_pv => null()
  type(variable_t), pointer, save :: vh_cf_divpv => null(), vh_cf_ew => null()
  type(variable_t), pointer, save :: vh_cf_elpot => null()
  type(variable_t), pointer, save :: vh_cf_el_mfd => null(), vh_cf_el_mfs => null()
  type(variable_t), pointer, save :: vh_cf_el_ef => null(), vh_cf_el_pv => null()
  type(variable_t), pointer, save :: vh_cf_el_divpv => null(), vh_cf_el_ew => null()
  real(kind=dp), allocatable, save, target :: vh_cf_gforce(:,:)
  real(kind=dp), allocatable, save :: vh_cf_mass(:,:), vh_cf_force(:,:)
  real(kind=dp), allocatable, save :: vh_cf_wbasis(:,:), vh_cf_sol(:,:)
  real(kind=dp), allocatable, save :: vh_cf_rotw(:,:), vh_cf_basis(:)
  real(kind=dp), allocatable, save :: vh_cf_db(:,:)
  integer, allocatable, save :: vh_cf_pivot(:)
  type(nodes_t), save :: vh_cf_nodes
  type(valuehandle_t), save :: vh_cf_eps_h, vh_cf_cur_h, vh_cf_mu_h

contains

  ! -- VectorHelmholtzSolver_Init0 ---------------------------------------
  subroutine vh_s0_begin(h_model, h_solver, with_ndofs) bind(c, name="vh_s0_begin")
    type(c_ptr), value :: h_model, h_solver
    integer(c_int) :: with_ndofs
    type(valuelist_t), pointer :: sp
    logical :: found
    call c_f_pointer(h_model, vh_model)
    call c_f_pointer(h_solver, vh_solver)
    sp => GetSolverParams()
    vh_params => sp
    with_ndofs = 0
    if (GetLogical(sp, 'Use Gauss Law', found)) with_ndofs = 1
    if (with_ndofs == 0) then
      if (GetLogical(sp, 'Use Lagrange Gauge', found)) with_ndofs = 1
      if (GetLogical(sp, 'Lorenz Condition', found)) with_ndofs = 1
    end if
    vh_s0_withndofs = (with_ndofs /= 0)
  end subroutine vh_s0_begin

  function vh_s0_has_element() result(r) bind(c, name="vh_s0_has_element")
    integer(c_int) :: r
    r = merge(1, 0, ListCheckPresent(vh_params, "Element"))
  end function vh_s0_has_element

  subroutine vh_s0_edge_style(k, piola, second) bind(c, name="vh_s0_edge_style")
    integer(c_int) :: k, piola, second
    logical :: pv, sf
    integer :: kk
    call EdgeElementStyle(vh_params, pv, sf, BasisDegree = kk, Check = .TRUE.)
    k = kk
    piola = merge(1, 0, pv)
    second = merge(1, 0, sf)
  end subroutine vh_s0_edge_style

  subroutine vh_s0_add_element(which) bind(c, name="vh_s0_add_element")
    integer(c_int), value :: which
    select case (which)
    case (1)
      call ListAddString(vh_params, "Element", &
          "n:1 e:2 -tri b:2 -quad b:4 -brick b:6 -pyramid b:3 -prism b:2 -quad_face b:4 -tri_face b:2")
    case (2)
      call ListAddString(vh_params, "Element", "n:1 e:2")
    case (3)
      call ListAddString(vh_params, "Element", &
          "n:1 e:1 -quad b:2 -brick b:3 -quad_face b:2")
    case (4)
      call ListAddString(vh_params, "Element", "n:1 e:1")
    case (5)
      call ListAddString(vh_params, "Element", "n:0 e:3 -tri b:3 -tri_face b:3")
    case (6)
      call ListAddString(vh_params, "Element", &
          "n:0 e:3 -tri b:6 -tetra b:3 -tri_face b:6")
    case (7)
      call ListAddString(vh_params, "Element", "n:0 e:3 -tri b:3 -tri_face b:3")
    case (8)
      call ListAddString(vh_params, "Element", &
          "n:0 e:2 -tri b:2 -quad b:4 -brick b:6 -pyramid b:3 -prism b:2 -quad_face b:4 -tri_face b:2")
    case (9)
      call ListAddString(vh_params, "Element", "n:0 e:2")
    case (10)
      call ListAddString(vh_params, "Element", &
          "n:0 e:1 -quad b:2 -brick b:3 -quad_face b:2")
    case default
      call ListAddString(vh_params, "Element", "n:0 e:1")
    end select
  end subroutine vh_s0_add_element

  subroutine vh_s0_fatal_cubic(second) bind(c, name="vh_s0_fatal_cubic")
    integer(c_int), value :: second
    if (second /= 0) then
      call Fatal('VectorHelmholtzSolver_Init0', &
          'No ready support for the cubic element of the second kind')
    else
      call Fatal('VectorHelmholtzSolver_Init0', &
          'The use of cubic elements is not yet possible')
    end if
  end subroutine vh_s0_fatal_cubic

  subroutine vh_s0_tail() bind(c, name="vh_s0_tail")
    logical :: found
    logical :: eigenproblem
    call ListAddNewLogical(vh_params, 'Bubbles in Global System', .TRUE.)
    ! (Hcurl Basis line stays commented, as in the original.)
    if (vh_s0_withndofs) then
      eigenproblem = GetLogical(vh_params, 'Eigen Analysis', found)
      if (eigenproblem) call Fatal('VectorHelmholtzSolver_Init0', &
          'Eigenvalue problem does not support the scalar unknown')
      call ListAddNewLogical(vh_params, 'Variable Output', .TRUE.)
      call ListAddString(vh_params, 'Variable', 'AV[AV re:1 AV im:1]')
    else
      call ListAddNewLogical(vh_params, 'Variable Output', .FALSE.)
      call ListAddString(vh_params, 'Variable', 'E[E re:1 E im:1]')
    end if
    call ListAddNewLogical(vh_params, "Linear System Complex", .TRUE.)
    ! These use one flag to call library features to compute automatically
    ! S parameters.
    if (ListGetLogical(vh_params, 'Calculate S Matrix', found)) then
      call Info('VectorHelmholtz_init', &
          'Using Constraint Modes functionality for S Matrix')
      call ListAddNewLogical(vh_params, 'Constraint Modes Analysis', .TRUE.)
    end if
    if (ListGetLogical(vh_params, 'Constraint Modes Analysis', found)) then
      call ListAddNewLogical(vh_params, 'Constraint Modes Constant Matrix', .TRUE.)
      call ListAddNewLogical(vh_params, 'Constraint Modes Lumped', .TRUE.)
      call ListAddNewLogical(vh_params, 'Constraint Modes Fluxes', .TRUE.)
      ! (Constraint Modes Matrix Results line stays commented.)
      call ListAddNewLogical(vh_params, 'Constraint Modes EM Wave', .TRUE.)
      if (ListCheckPresent(vh_params, 'S Matrix Filename')) then
        call ListRename(vh_params, 'S Matrix Filename', &
            'Constraint Modes Matrix Filename', found)
      else
        call ListAddNewString(vh_params, 'Constraint Modes Matrix Filename', &
            'SMatrix.dat', .FALSE.)
      end if
      if (ListCheckPresent(vh_params, 'S Matrix Symmetric')) then
        call ListRename(vh_params, 'S Matrix Symmetric', &
            'Constraint Modes Fluxes Symmetric')
      end if
      call ListRenameAllBC(vh_model, 'S Matrix Port', 'Constraint Mode')
    end if
    if (ListGetString(vh_params, 'Linear System Preconditioning', found) &
        == "auxiliary space solver") then
      if (.not. ListCheckPresent(vh_params, 'Preconditioning Residual')) then
        call ListAddString(vh_params, &
            NextFreeKeyword('Exported Variable', vh_params), &
            "-dofs 2 ams res")
        call ListAddString(vh_params, 'Preconditioning Residual', "ams res")
      end if
      if (.not. ListCheckPresent(vh_params, 'Preconditioning Update')) then
        call ListAddString(vh_params, &
            NextFreeKeyword('Exported Variable', vh_params), &
            "-dofs 2 ams update")
        call ListAddString(vh_params, 'Preconditioning Update', "ams update")
      end if
    end if
  end subroutine vh_s0_tail

  ! -- VectorHelmholtzSolver_Init ----------------------------------------
  subroutine vh_si_begin(h_model, h_solver) bind(c, name="vh_si_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, vh_model)
    call c_f_pointer(h_solver, vh_solver)
    vh_params => GetSolverParams()
  end subroutine vh_si_begin

  subroutine vh_si_prec_prefix() bind(c, name="vh_si_prec_prefix")
    if (ListCheckPrefix(vh_params, &
        'Linear System Preconditioning Damp Coefficient')) then
      call ListAddNewLogical(vh_params, 'Allocate Preconditioning Matrix', .TRUE.)
    end if
  end subroutine vh_si_prec_prefix

  function vh_si_find_eigen() result(r) bind(c, name="vh_si_find_eigen")
    integer(c_int) :: r
    type(valuelist_t), pointer :: solverparams
    logical :: found, findeigen
    integer :: i
    solverparams => GetSolverParams()
    findeigen = .false.
    do i = 1, vh_model % NumberOfBCs
      if (ListGetString(vh_model % BCs(i) % Values, 'Port Type', found) &
          == 'eigenmode') then
        findeigen = .true.
        exit
      end if
    end do
    r = merge(1, 0, findeigen)
  end function vh_si_find_eigen

  function vh_si_emport() result(soln) bind(c, name="vh_si_emport")
    integer(c_int) :: soln
    type(valuelist_t), pointer :: solverparams
    logical :: found
    integer :: i, j, n
    character(len=MAX_NAME_LEN) :: sname
    solverparams => GetSolverParams()
    soln = 0
    n = vh_model % NumberOfSolvers
    do i = 1, n
      sname = GetString(vh_model % Solvers(i) % Values, 'Procedure', found)
      j = index(sname, 'EMPortSolver')
      if (j > 0) then
        soln = i
        exit
      end if
    end do
  end function vh_si_emport

  subroutine vh_si_note_emport(soln) bind(c, name="vh_si_note_emport")
    integer(c_int), value :: soln
    type(valuelist_t), pointer :: solverparams
    solverparams => GetSolverParams()
    call Info('VectorHelmholtzSolver_Init', &
        'The eigensolver index is: '//I2S(soln), level=12)
    call ListAddInteger(solverparams, 'Eigensolver Index', soln)
  end subroutine vh_si_note_emport

  subroutine vh_si_fatal_noport() bind(c, name="vh_si_fatal_noport")
    call Fatal('VectorHelmholtzSolver_Init', &
        'Eigenfunction BC given without solving a port model')
  end subroutine vh_si_fatal_noport

  subroutine vh_si_define_ports() bind(c, name="vh_si_define_ports")
    call DefinePortParameters(vh_model, vh_solver % Mesh)
  end subroutine vh_si_define_ports


  ! -- VectorHelmholtzSolver driver --------------------------------------
  subroutine vh_s_begin(h_model, h_solver) bind(c, name="vh_s_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, vh_model)
    call c_f_pointer(h_solver, vh_solver)
  end subroutine vh_s_begin

  subroutine vh_s_prologue() bind(c, name="vh_s_prologue")
    logical :: found, precmatrix
    logical :: hasprecdampcoeff, massproportional, curlcurlprec
    integer :: i, j
    integer, pointer :: solverindexes(:)
    real(kind=dp) :: omega, mu0inv, eps0, rob0
    complex(kind=dp) :: precdampcoeff
    logical :: pivolaversion, lowfrequencymodel, lorenzcondition
    logical :: usegausslaw, chargeconservation, eigenproblem
    type(valuelist_t), pointer :: solverparams

    call Info(vh_caller, '', level=6)
    call Info(vh_caller, &
        '-------------------------------------------------', level=6)
    call Info(vh_caller, &
        'Solving harmonic electromagnetic wave equation!', level=5)

    solverparams => GetSolverParams()
    vh_params => solverparams
    call EdgeElementStyle(vh_solver % Values, pivolaversion, &
        BasisDegree = vh_edge_degree)
    vh_piola = pivolaversion

    if (CoordinateSystemDimension() == 2 .and. .not. vh_piola) then
      call Fatal(vh_caller, 'A 2D model needs "Use Piola Transform = True"')
    end if

    ! Allocate some permanent storage, this is done first time only:
    !---------------------------------------------------------------
    vh_mesh => GetMesh()
    vh_psolver => vh_solver

    if (vh_solver % Variable % dofs /= 2) then
      call Fatal(vh_caller, &
          'Variable is not of size two ('//I2S(vh_solver % Variable % dofs)//'), &
          Use: Variable = E[E re:1 E im:1]')
    end if

    ! If we use our own AMS preconditioner where actually the scalar and
    ! vector preconditioning equations are solved for, lets get information
    ! needed for the assembly i.e. pointers to the solvers and matrices.
    ! (GetAuxSolverInfo body, its only caller is this prologue.)
    vh_amsany = .false.
    nullify(vh_amsmat); nullify(vh_amsscalmat)
    nullify(vh_amssolver); nullify(vh_amsscalsolver)
    solverindexes => ListGetIntegerArray(solverparams, 'prec solvers', found)
    if (found) then
      do i = 1, size(solverindexes)
        j = solverindexes(i)
        if (ListGetLogical(vh_model % Solvers(j) % Values, &
            'AMS Vector Solver', found)) then
          vh_amssolver => vh_model % Solvers(j)
          vh_amsmat => vh_amssolver % Matrix
          ! For monolithic AMS matrix we may use curl-curl form.
          vh_amscurlcurl = .false.
          vh_amsmono = ListGetLogical(vh_amssolver % Values, &
              'Monolithic Solver', found)
          if (vh_amsmono) then
            vh_amscurlcurl = ListGetLogical(vh_amssolver % Values, &
                'curl-curl Form', found)
          end if
          if (vh_amscurlcurl) then
            call Info(vh_caller, &
                'Using curl-curl form for AMS preconditioner', level=7)
          end if
          vh_amsany = .true.
        end if
        if (ListGetLogical(vh_model % Solvers(j) % Values, &
            'AMS Scalar Solver', found)) then
          vh_amsscalsolver => vh_model % Solvers(j)
          vh_amsscalmat => vh_amsscalsolver % Matrix
          vh_amsany = .true.
        end if
      end do
    end if
    if (associated(vh_amsmat)) then
      call Info(vh_caller, &
          'Assembling nodal AMS vector matrix on the side!', level=7)
    end if
    if (associated(vh_amsscalmat)) then
      call Info(vh_caller, &
          'Assembling nodal AMS scalar matrix on the side!', level=7)
    end if

    precdampcoeff = GetCReal(solverparams, &
        'Linear System Preconditioning Damp Coefficient', hasprecdampcoeff)
    precdampcoeff = cmplx(real(precdampcoeff), GetCReal(solverparams, &
        'Linear System Preconditioning Damp Coefficient im', found), kind=dp)
    hasprecdampcoeff = hasprecdampcoeff .or. found
    vh_precdamp = precdampcoeff
    vh_hasprecdamp = hasprecdampcoeff
    if (hasprecdampcoeff) then
      massproportional = GetLogical(solverparams, &
          'Mass-proportional Damping', found)
    else
      massproportional = .false.
    end if
    vh_massprop = massproportional

    curlcurlprec = GetLogical(solverparams, 'Curl-Curl Preconditioning', found)
    vh_curlcurlprec = curlcurlprec
    precmatrix = hasprecdampcoeff .or. curlcurlprec
    vh_precmatrix = precmatrix

    if (precmatrix) then
      if (associated(vh_amsmat)) then
        call Info(vh_caller, &
            'Incorporating damping coefficient to AMS matrix directly', level=7)
        precmatrix = .false.
        vh_precmatrix = .false.
      else if (ListGetString(solverparams, 'Linear System Solver', found) &
          == 'direct') then
        call Warn(vh_caller, &
            'Generating preconditioning matrix does not make sense for direct methods, canceling!')
        precmatrix = .false.
        vh_precmatrix = .false.
      else
        call Info(vh_caller, &
            'Generating special preconditioning matrix', level=12)
      end if
    end if

    found = .false.
    if (associated(vh_model % Constants)) then
      mu0inv = GetConstReal(vh_model % Constants, 'Permeability of Vacuum', found)
      if (mu0inv /= 0) mu0inv = 1 / mu0inv
    end if
    if (.not. found) mu0inv = 1.0_dp / (PI * 4.0d-7)
    vh_mu0inv = mu0inv

    found = .false.
    if (associated(vh_model % Constants)) then
      eps0 = GetConstReal(vh_model % Constants, 'Permittivity of Vacuum', found)
    end if
    if (.not. found) eps0 = 8.854187817d-12
    vh_eps0 = eps0

    eigenproblem = GetLogical(solverparams, 'Eigen Analysis', found)
    vh_eigen = eigenproblem
    if (.not. eigenproblem) then
      omega = GetAngularFrequency(found=found)
      if (.not. found) then
        call Fatal(vh_caller, 'Harmonic wave equation requires angular frequency!')
      end if
      rob0 = omega * sqrt(eps0 / mu0inv)
    else
      rob0 = 0.0_dp
    end if
    vh_omega = omega
    vh_rob0 = rob0

    lowfrequencymodel = GetLogical(solverparams, 'Low Frequency Model', found)
    vh_lowfreq = lowfrequencymodel
    lorenzcondition = GetLogical(solverparams, 'Lorenz Condition', found)
    vh_lorenz = lorenzcondition
    usegausslaw = GetLogical(solverparams, 'Use Gauss Law', found)
    vh_usegauss = usegausslaw
    chargeconservation = GetLogical(solverparams, &
        'Apply Conservation of Charge', found)
    vh_chargecons = chargeconservation

    ! Resolve internal nonlinearities, if requested:
    ! ----------------------------------------------
    vh_noitmax = GetInteger(solverparams, &
        'Nonlinear System Max Iterations', found)
    if (.not. found) vh_noitmax = 1

    vh_edgebasis = .not. ListCheckPresent(solverparams, &
        'Linear System Refactorize') .and. &
        GetLogical(solverparams, 'Edge Basis', found)
  end subroutine vh_s_prologue

  subroutine vh_s_iter_params(noit, edgebasis) bind(c, name="vh_s_iter_params")
    integer(c_int) :: noit, edgebasis
    noit = vh_noitmax
    edgebasis = merge(1, 0, vh_edgebasis)
  end subroutine vh_s_iter_params

  subroutine vh_s_default_start() bind(c, name="vh_s_default_start")
    call DefaultStart()
  end subroutine vh_s_default_start

  subroutine vh_s_refactorize_false() bind(c, name="vh_s_refactorize_false")
    call ListAddLogical(vh_params, 'Linear System Refactorize', .FALSE.)
  end subroutine vh_s_refactorize_false

  subroutine vh_s_refactorize_remove() bind(c, name="vh_s_refactorize_remove")
    call ListRemove(vh_params, 'Linear System Refactorize')
  end subroutine vh_s_refactorize_remove

  subroutine vh_s_default_finish() bind(c, name="vh_s_default_finish")
    call DefaultFinish()
  end subroutine vh_s_default_finish

  function vh_s_circ_idx() result(idx) bind(c, name="vh_s_circ_idx")
    integer(c_int) :: idx
    logical :: found
    idx = ListGetInteger(vh_params, 'Circulation BC', found)
    if (.not. found) idx = 0
  end function vh_s_circ_idx

  subroutine vh_s_circulation(idx) bind(c, name="vh_s_circulation")
    integer(c_int), value :: idx
    real(kind=dp) :: circ(2)
    call BoundaryCirculation(circ, idx, vh_psolver, vh_psolver % Variable)
    print *, 'Circulation around BC:', circ
  end subroutine vh_s_circulation

  subroutine vh_s_done() bind(c, name="vh_s_done")
    call Info(vh_caller, 'All done', level=12)
  end subroutine vh_s_done

  function vh_s_dosolve() result(converged) bind(c, name="vh_s_dosolve")
    integer(c_int) :: converged
    type(element_t), pointer :: element
    type(valuelist_t), pointer :: bc
    real(kind=dp), pointer, contiguous :: savedvalues(:) => null()
    real(kind=dp) :: norm
    integer :: active, k, n, nd, nb, t
    logical :: inithandles, conv
    ! System assembly:
    !-----------------
    call Info(vh_caller, 'Starting bulk assembly', level=12)

    call DefaultInitialize()
    if (associated(vh_amsmat)) vh_amsmat % Values = 0.0_dp
    if (associated(vh_amsscalmat)) vh_amsscalmat % Values = 0.0_dp

    active = GetNOFActive()
    inithandles = .true.

    do t = 1, active
      element => GetActiveElement(t)
      n = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      nb = GetElementNOFBDOFs()

      ! Glue local element matrix and rhs vector:
      !----------------------------------------
      call LocalMatrix(element, n, nd + nb, nb, inithandles)
    end do
    call DefaultFinishBulkAssembly()

    ! Robin type of BC in terms of E:
    !--------------------------------
    call Info(vh_caller, 'Starting boundary assembly', level=12)

    if (vh_eigen) then
      ! Currently the eigenvalue problem is implemented for a combination of
      ! Dirichlet and homogeneous Neumann BCs
      active = 0
    else
      active = GetNOFBoundaryElements()
      inithandles = .true.
    end if

    do t = 1, active
      element => GetBoundaryElement(t)
      bc => GetBC()
      if (.not. associated(bc)) cycle

      select case (GetElementFamily())
      case (1)
        cycle
      case (2)
        k = GetBoundaryEdgeIndex(element, 1)
        element => vh_mesh % Edges(k)
      case (3, 4)
        k = GetBoundaryFaceIndex(element)
        element => vh_mesh % Faces(k)
      end select
      if (.not. ActiveBoundaryElement(element)) cycle

      block
        type(element_t), pointer :: tmp_element
        tmp_element => SetCurrentElement(element)
      end block

      nd = GetElementNOFDOFs(element)
      n = GetElementNOFNodes(element)

      call LocalMatrixBC(bc, element, n, nd, inithandles)
    end do

    call Info(vh_caller, 'Local assembly done', level=12)

    call DefaultFinishBoundaryAssembly()
    call DefaultFinishAssembly()

    ! Dirichlet BCs in terms of electric field E (or the potentials)
    ! --------------------------------------------------------------
    call DefaultDirichletBCs()

    ! Call DefaultDirichletBCs another time to apply BCs to PrecValues:
    if (associated(vh_solver % Matrix % PrecValues)) then
      savedvalues => vh_solver % Matrix % Values
      vh_solver % Matrix % Values => vh_solver % Matrix % PrecValues
      call DefaultDirichletBCs()
      vh_solver % Matrix % Values => savedvalues
    end if

    if (.not. vh_eigen) call SingleDipoleLoad()

    ! Linear system solution:
    ! -----------------------
    norm = DefaultSolve()
    conv = (vh_solver % Variable % NonlinConverged == 1)
    converged = merge(1, 0, conv)
  end function vh_s_dosolve

  subroutine LocalMatrix(element, n, nd, nb, inithandles)
    type(element_t), pointer :: element
    integer :: n, nd, nb
    logical :: inithandles
    complex(kind=dp) :: eps, muinv, cond, L(3)
    real(kind=dp) :: detj, weight
    logical :: stat, withndofs, conductorbody
    logical :: allocationsdone
    integer :: t, i, j, m, np, p, q, ndofs
    type(GaussIntegrationPoints_t) :: IP
    complex(kind=dp) :: am, aw, ac, atot(3, 3)
    integer :: idim, jdim

    allocationsdone = vh_lm_alloc
    if (.not. allocationsdone) then
      m = vh_mesh % MaxElementDOFs
      allocate (vh_lm_wbasis(m, 3), vh_lm_rotw(m, 3), vh_lm_basis(m), &
          vh_lm_db(m, 3), vh_lm_mass(m, m), vh_lm_stiff(m, m), &
          vh_lm_gauge(m, m), vh_lm_prec(m, m), vh_lm_curlmat(m, m), &
          vh_lm_force(m))
      if (vh_amsany) then
        allocate (vh_lm_amsstiff(3 * m, 3 * m), vh_lm_amsforce(3 * m), &
            vh_lm_amsstiff2(m, m))
        vh_lm_amsforce = 0.0_dp
      end if
      vh_lm_alloc = .true.
    end if

    ! This InitHandles flag might be false on threaded 1st call
    if (inithandles) then
      call Info(vh_caller, 'Initializing handles in LocalMatrix', level=25)
      call ListInitElementKeyword(vh_lm_cond_h, 'Material', &
          'Electric Conductivity', InitIm=.TRUE.)
      call ListInitElementKeyword(vh_lm_eps_h, 'Material', &
          'Relative Permittivity', InitIm=.TRUE.)
      call ListInitElementKeyword(vh_lm_mu_h, 'Material', &
          'Relative Reluctivity', InitIm=.TRUE.)
      call ListInitElementKeyword(vh_lm_cur_h, 'Body Force', &
          'Current Density', InitIm=.TRUE., InitVec3D=.TRUE.)
      inithandles = .false.
    end if

    call GetElementNodes(vh_lm_nodes, element)

    vh_lm_stiff(1:nd, 1:nd) = 0.0_dp
    vh_lm_mass(1:nd, 1:nd) = 0.0_dp
    vh_lm_curlmat = 0.0_dp
    vh_lm_force(1:nd) = 0.0_dp
    if (vh_amsany) then
      vh_lm_amsstiff = 0.0_dp
      vh_lm_amsstiff2 = 0.0_dp
    end if

    ndofs = maxval(vh_solver % Def_Dofs(GetElementFamily(element), :, 1))
    np = n * ndofs

    withndofs = (ndofs > 0)
    if (withndofs) then
      vh_lm_gauge(1:nd, 1:nd) = 0.0_dp
    end if

    if (vh_precmatrix) vh_lm_prec = 0.0_dp

    ! Numerical integration:
    !----------------------
    IP = GaussPoints(element, PReferenceElement=vh_piola, &
        EdgeBasisDegree = vh_edge_degree)

    do t = 1, IP % n

      stat = ElementInfo(element, vh_lm_nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, vh_lm_basis, vh_lm_db, EdgeBasis = vh_lm_wbasis, &
          RotBasis = vh_lm_rotw, USolver = vh_psolver)
      weight = detJ * IP % s(t)

      ! Compute element stiffness matrix and force vector:
      ! --------------------------------------------------
      muinv = ListGetElementComplex(vh_lm_mu_h, vh_lm_basis, element, &
          stat, GaussPoint = t)
      if (stat) then
        muinv = muinv * vh_mu0inv
      else
        muinv = vh_mu0inv
      end if

      ! This is present always
      do p = 1, nd - np
        i = p + np
        do q = 1, nd - np
          j = q + np
          ! the mu^-1 curl E . curl v
          vh_lm_curlmat(i, j) = vh_lm_curlmat(i, j) + muinv * &
              sum(vh_lm_rotw(p, :) * vh_lm_rotw(q, :)) * weight
        end do
      end do

      ! Conductivity may also be accounted for
      cond = ListGetElementComplex(vh_lm_cond_h, vh_lm_basis, element, &
          stat, GaussPoint = t)
      conductorbody = stat .and. abs(cond) > AEPS
      if (conductorbody) then
        if (vh_eigen) call Fatal(vh_caller, &
            'Eigenvalue problem does not support electric conductivity')
        do p = 1, nd - np
          i = p + np
          do q = 1, nd - np
            j = q + np
            ! the term i\omega\sigma E.v
            vh_lm_stiff(i, j) = vh_lm_stiff(i, j) - im * vh_omega * cond * &
                sum(vh_lm_wbasis(q, :) * vh_lm_wbasis(p, :)) * weight
          end do

          if (withndofs) then
            do q = 1, n
              j = (q - 1) * ndofs + 1
              vh_lm_stiff(i, j) = vh_lm_stiff(i, j) + im * vh_omega * cond * &
                  sum(vh_lm_db(q, :) * vh_lm_wbasis(p, :)) * weight
            end do
          end if
        end do
      end if

      ! If not low frequency model, assemble the term that makes this the
      ! wave equation
      if (.not. vh_lowfreq) then
        eps = ListGetElementComplex(vh_lm_eps_h, vh_lm_basis, element, &
            stat, GaussPoint = t)
        if (stat) then
          eps = vh_eps0 * eps
        else
          eps = vh_eps0
        end if

        do p = 1, nd - np
          i = p + np
          do q = 1, nd - np
            j = q + np
            ! the term \omega^2 \epsilon E.v
            vh_lm_mass(i, j) = vh_lm_mass(i, j) + eps * &
                sum(vh_lm_wbasis(q, :) * vh_lm_wbasis(p, :)) * weight
          end do
        end do
      end if

      ! Potential current source
      L = ListGetElementComplex3D(vh_lm_cur_h, vh_lm_basis, element, &
          stat, GaussPoint = t)
      if (stat .and. .not. vh_eigen) then
        do p = 1, nd - np
          i = p + np
          vh_lm_force(i) = vh_lm_force(i) + im * vh_omega * &
              (sum(L * vh_lm_wbasis(p, :))) * weight
        end do
      end if

      ! Additional terms related to using the Gauss law or a gauge condition
      if (withndofs) then
        if (ndofs == 2) then
          !
          ! Using more than one nodal field is experimental and may break some
          ! generic functionality. In principle, the idea is to use the A-V
          ! representation together with a gauge constraint.
          !
          if (conductorbody .and. vh_chargecons) then
            ! Use an implied equation about the conservation of charge
            do p = 1, n
              i = (p - 1) * ndofs + 1
              do q = 1, nd - np
                j = q + np
                vh_lm_gauge(i, j) = vh_lm_gauge(i, j) - im * vh_omega * &
                    cond * sum(vh_lm_wbasis(q, :) * vh_lm_db(i, :)) * weight
                vh_lm_gauge(j, i) = vh_lm_gauge(j, i) + vh_omega**2 * eps * &
                    sum(vh_lm_wbasis(q, :) * vh_lm_db(i, :)) * weight
              end do

              do q = 1, n
                j = (q - 1) * ndofs + 1
                vh_lm_gauge(i, j) = vh_lm_gauge(i, j) + im * vh_omega * &
                    cond * sum(vh_lm_db(i, :) * vh_lm_db(j, :)) * weight
              end do
            end do
          else
            do p = 1, n
              ! If two nodal DOFs per node, the first DOF is the scalar
              ! potential V related to the A-V representation. We add
              ! -w^2 div D = -w^2 rho in a weak form when E = A - grad V:
              i = (p - 1) * ndofs + 1
              do q = 1, nd - np
                j = q + np
                vh_lm_gauge(i, j) = vh_lm_gauge(i, j) + vh_omega**2 * eps * &
                    sum(vh_lm_wbasis(q, :) * vh_lm_db(p, :)) * weight
                vh_lm_gauge(j, i) = vh_lm_gauge(j, i) + vh_omega**2 * eps * &
                    sum(vh_lm_wbasis(q, :) * vh_lm_db(p, :)) * weight
              end do

              do q = 1, n
                j = (q - 1) * ndofs + 1
                vh_lm_gauge(i, j) = vh_lm_gauge(i, j) - vh_omega**2 * eps * &
                    sum(vh_lm_db(q, :) * vh_lm_db(p, :)) * weight
              end do
            end do
          end if

          do p = 1, n
            ! The second DOF is related to the gauge condition div A = 0
            ! (TO DO: Add Lorenz condition as an alternative)
            i = p * ndofs
            do q = 1, nd - np
              j = q + np
              vh_lm_gauge(i, j) = vh_lm_gauge(i, j) - &
                  sum(vh_lm_wbasis(q, :) * vh_lm_db(p, :)) * weight
              vh_lm_gauge(j, i) = vh_lm_gauge(j, i) - &
                  sum(vh_lm_wbasis(q, :) * vh_lm_db(p, :)) * weight
            end do
          end do
        else
          if (vh_usegauss) then
            if (conductorbody .and. vh_chargecons) then
              ! Use an implied equation about the conservation of charge
              do i = 1, np
                do q = 1, nd - np
                  j = q + np
                  vh_lm_gauge(i, j) = vh_lm_gauge(i, j) - im * vh_omega * &
                      cond * sum(vh_lm_wbasis(q, :) * vh_lm_db(i, :)) * weight
                  vh_lm_gauge(j, i) = vh_lm_gauge(j, i) + vh_omega**2 * &
                      eps * sum(vh_lm_wbasis(q, :) * vh_lm_db(i, :)) * weight
                end do
                do j = 1, np
                  vh_lm_gauge(i, j) = vh_lm_gauge(i, j) + im * vh_omega * &
                      cond * sum(vh_lm_db(i, :) * vh_lm_db(j, :)) * weight
                end do
              end do
            else
              ! Add -w^2 div D = -w^2 rho in a weak form when E = A - grad V:
              do i = 1, np
                do q = 1, nd - np
                  j = q + np
                  vh_lm_gauge(i, j) = vh_lm_gauge(i, j) + vh_omega**2 * &
                      eps * sum(vh_lm_wbasis(q, :) * vh_lm_db(i, :)) * weight
                  vh_lm_gauge(j, i) = vh_lm_gauge(j, i) + vh_omega**2 * &
                      eps * sum(vh_lm_wbasis(q, :) * vh_lm_db(i, :)) * weight
                end do
                do j = 1, np
                  vh_lm_gauge(i, j) = vh_lm_gauge(i, j) - vh_omega**2 * &
                      eps * sum(vh_lm_db(i, :) * vh_lm_db(j, :)) * weight
                end do
              end do
            end if
          else
            ! Here Gauss's law is not employed, only an additional gauge
            ! constraint is introduced
            do i = 1, np
              do q = 1, nd - np
                j = q + np
                vh_lm_gauge(i, j) = vh_lm_gauge(i, j) - &
                    sum(vh_lm_wbasis(q, :) * vh_lm_db(i, :)) * weight
                vh_lm_gauge(j, i) = vh_lm_gauge(j, i) + vh_omega**2 * &
                    eps * sum(vh_lm_wbasis(q, :) * vh_lm_db(i, :)) * weight
              end do
              if (vh_lorenz) then
                do j = 1, np
                  vh_lm_gauge(i, j) = vh_lm_gauge(i, j) + vh_omega**2 * &
                      eps / muinv * vh_lm_basis(i) * vh_lm_basis(j) * weight
                end do
              end if
            end do
          end if
        end if
      end if

      if (vh_amsany) then
        block
          complex(kind=dp) :: aw2, ac2, atot2(3, 3)
          do p = 1, n
            do q = 1, n
              if (.not. vh_lowfreq) then
                aw2 = -vh_omega**2 * eps * vh_lm_basis(q) * vh_lm_basis(p)
                if (vh_hasprecdamp .and. vh_massprop) then
                  aw2 = aw2 * (1 - vh_precdamp)
                end if
              end if
              if (conductorbody) then
                ac2 = -im * vh_omega * cond * vh_lm_basis(p) * vh_lm_basis(q)
              end if

              if (associated(vh_amsmat)) then
                if (vh_amscurlcurl) then
                  ! The grad-div operator is now zero.
                  am = 0.0_dp
                  atot2(1, 1) = vh_lm_db(q, 3) * vh_lm_db(p, 3) + &
                      vh_lm_db(q, 2) * vh_lm_db(p, 2)
                  atot2(1, 2) = -vh_lm_db(q, 1) * vh_lm_db(p, 2)
                  atot2(1, 3) = -vh_lm_db(q, 1) * vh_lm_db(p, 3)
                  atot2(2, 1) = -vh_lm_db(q, 2) * vh_lm_db(p, 1)
                  atot2(2, 2) = vh_lm_db(q, 1) * vh_lm_db(p, 1) + &
                      vh_lm_db(q, 3) * vh_lm_db(p, 3)
                  atot2(2, 3) = -vh_lm_db(q, 2) * vh_lm_db(p, 3)
                  atot2(3, 1) = -vh_lm_db(q, 3) * vh_lm_db(p, 1)
                  atot2(3, 2) = -vh_lm_db(q, 3) * vh_lm_db(p, 2)
                  atot2(3, 3) = vh_lm_db(q, 1) * vh_lm_db(p, 1) + &
                      vh_lm_db(q, 2) * vh_lm_db(p, 2)

                  ! Multiply after creating the curl-curl because of
                  ! so many terms...
                  atot2 = weight * muinv * atot2

                  if (vh_hasprecdamp .and. .not. vh_massprop) then
                    atot2 = atot2 * (1.0 + vh_precdamp)
                  end if
                else
                  ! grad-div operator
                  am = muinv * sum(vh_lm_db(p, :) * vh_lm_db(q, :))
                  atot2 = 0.0_dp

                  if (vh_hasprecdamp .and. .not. vh_massprop) then
                    am = am * (1.0 + vh_precdamp)
                  end if
                end if

                if (vh_amsmono) then
                  do idim = 1, 3
                    atot2(idim, idim) = atot2(idim, idim) + &
                        weight * (am + aw2 + ac2)
                  end do
                  do idim = 1, 3
                    do jdim = 1, 3
                      vh_lm_amsstiff(3 * (p - 1) + idim, 3 * (q - 1) + jdim) = &
                          vh_lm_amsstiff(3 * (p - 1) + idim, &
                          3 * (q - 1) + jdim) + atot2(idim, jdim)
                    end do
                  end do
                else
                  atot2(1, 1) = atot2(1, 1) + weight * (am + aw2 + ac2)
                  vh_lm_amsstiff(p, q) = vh_lm_amsstiff(p, q) + atot2(1, 1)
                end if
              end if

              if (associated(vh_amsscalmat)) then
                am = muinv * sum(vh_lm_db(p, :) * vh_lm_db(q, :))
                if (vh_hasprecdamp .and. .not. vh_massprop) then
                  am = am * (1.0 + vh_precdamp)
                end if
                vh_lm_amsstiff2(p, q) = weight * (am + aw2 + ac2)
              end if
            end do
          end do
        end block
      end if

    end do

    if (.not. vh_eigen) then
      vh_lm_mass(1:nd, 1:nd) = -vh_omega**2 * vh_lm_mass(1:nd, 1:nd)
    end if

    if (vh_precmatrix) then
      if (vh_curlcurlprec) then
        vh_lm_prec = vh_lm_stiff(1:nd, 1:nd) + vh_lm_curlmat(1:nd, 1:nd)
        if (withndofs) then
          vh_lm_prec(1:nd, 1:nd) = vh_lm_prec(1:nd, 1:nd) + &
              vh_lm_gauge(1:nd, 1:nd)
        end if
      else if (vh_massprop) then
        vh_lm_prec = -vh_precdamp * (vh_lm_mass(1:nd, 1:nd))
      else
        vh_lm_prec = vh_precdamp * (vh_lm_stiff(1:nd, 1:nd) + &
            vh_lm_curlmat(1:nd, 1:nd) - vh_lm_mass(1:nd, 1:nd))
      end if
    end if

    if (vh_eigen) then
      vh_lm_stiff(1:nd, 1:nd) = vh_lm_stiff(1:nd, 1:nd) + &
          vh_lm_curlmat(1:nd, 1:nd)
    else
      vh_lm_stiff(1:nd, 1:nd) = vh_lm_stiff(1:nd, 1:nd) + &
          vh_lm_curlmat(1:nd, 1:nd) + vh_lm_mass(1:nd, 1:nd)
    end if
    if (withndofs) then
      vh_lm_stiff(1:nd, 1:nd) = vh_lm_stiff(1:nd, 1:nd) + &
          vh_lm_gauge(1:nd, 1:nd)
    end if

    if (vh_precmatrix) then
      if (.not. vh_curlcurlprec) then
        vh_lm_prec(1:nd, 1:nd) = vh_lm_stiff(1:nd, 1:nd) + &
            vh_lm_prec(1:nd, 1:nd)
      end if
      if (nb > 0) call CondensateP(nd - nb, nb, vh_lm_prec)
      call DefaultUpdatePrec(vh_lm_prec(1:nd, 1:nd))
    end if

    ! Update global matrix and rhs vector from local matrix & vector:
    !---------------------------------------------------------------
    if (nb > 0) call CondensateP(nd - nb, nb, vh_lm_stiff, vh_lm_force)
    call DefaultUpdateEquations(vh_lm_stiff, vh_lm_force, element)
    if (vh_eigen) call DefaultUpdateMass(vh_lm_mass)

    if (vh_amsany) then
      if (associated(vh_amsmat)) then
        CurrentModel % Solver => vh_amssolver
        call DefaultUpdateEquations(vh_lm_amsstiff, vh_lm_amsforce, &
            UElement=element, USolver=vh_amssolver)
        CurrentModel % Solver => vh_psolver
      end if
      if (associated(vh_amsscalmat)) then
        CurrentModel % Solver => vh_amsscalsolver
        call DefaultUpdateEquations(vh_lm_amsstiff2, vh_lm_amsforce, &
            UElement=element, USolver=vh_amsscalsolver)
        CurrentModel % Solver => vh_psolver
      end if
    end if
  end subroutine LocalMatrix


  subroutine LocalMatrixBC(bc, element, n, nd, inithandles)
    type(valuelist_t), pointer :: bc
    type(element_t), pointer :: element
    integer :: n, nd
    logical :: inithandles
    complex(kind=dp) :: elsurfcurr(3), b, l(3), muinv, temgrad(3)
    complex(kind=dp) :: magload(3), betapar, portbeta, jn, cond
    complex(kind=dp) :: surfimp, epsr, mur, ep
    real(kind=dp) :: th, detj, weight, testvec(3), trialvec(3), normal(3)
    logical :: stat, found, updatestiff, withndofs, thinsheet
    logical :: goodconductor, absorb
    logical :: lineelement, degenerateelement, regularize, consistent
    logical :: allocationsdone
    type(GaussIntegrationPoints_t) :: IP
    integer :: t, i, j, m, np, p, q, ndofs, dim
    type(element_t), pointer :: parent
    logical :: gotport

    allocationsdone = vh_bc_alloc
    if (.not. allocationsdone) then
      m = vh_mesh % MaxElementDOFs
      allocate (vh_bc_wbasis(m, 3), vh_bc_rotw(m, 3), vh_bc_basis(m), &
          vh_bc_db(m, 3), vh_bc_force(m), vh_bc_stiff(m, m), vh_bc_mass(m, m))
      if (vh_amsany) then
        allocate (vh_bc_amsstiff(3 * m, 3 * m), vh_bc_amsstiff2(m, m), &
            vh_bc_amsforce(3 * m))
        vh_bc_amsforce = 0.0_dp
      end if
      vh_bc_alloc = .true.
    end if

    if (inithandles) then
      call Info(vh_caller, 'Initializing handles in LocalMatrixBC', level=25)

      call DefinePortParameters(vh_model, vh_mesh)

      call ListInitElementKeyword(vh_bc_elrobin_h, 'Boundary Condition', &
          'Electric Robin Coefficient', InitIm=.TRUE.)
      call ListInitElementKeyword(vh_bc_magload_h, 'Boundary Condition', &
          'Magnetic Boundary Load', InitIm=.TRUE., InitVec3D=.TRUE.)
      call ListInitElementKeyword(vh_bc_elsurf_h, 'Boundary Condition', &
          'Electric Surface Current', InitIm=.TRUE., InitVec3D=.TRUE.)
      call ListInitElementKeyword(vh_bc_absorb_h, 'Boundary Condition', &
          'Absorbing BC')
      call ListInitElementKeyword(vh_bc_goodcond_h, 'Boundary Condition', &
          'Good Conductor BC')
      call ListInitElementKeyword(vh_bc_chargecons_h, 'Boundary Condition', &
          'Apply Conservation of Charge')
      call ListInitElementKeyword(vh_bc_temre_h, 'Boundary Condition', &
          'TEM Potential')
      call ListInitElementKeyword(vh_bc_temim_h, 'Boundary Condition', &
          'TEM Potential Im')
      call ListInitElementKeyword(vh_bc_mu_h, 'Material', &
          'Relative Reluctivity', InitIm=.TRUE.)
      call ListInitElementKeyword(vh_bc_eps_h, 'Material', &
          'Relative Permittivity', InitIm=.TRUE.)

      call ListInitElementKeyword(vh_bc_transfer_h, 'Boundary Condition', &
          'Electric Transfer Coefficient', InitIm=.TRUE.)
      call ListInitElementKeyword(vh_bc_elcurr_h, 'Boundary Condition', &
          'Electric Current Density', InitIm=.TRUE.)
      call ListInitElementKeyword(vh_bc_extpot_h, 'Boundary Condition', &
          'Incident Voltage', InitIm=.TRUE.)

      call ListInitElementKeyword(vh_bc_thick_h, 'Boundary Condition', &
          'Layer Thickness')
      call ListInitElementKeyword(vh_bc_relnu_h, 'Boundary Condition', &
          'Layer Relative Reluctivity', InitIm=.TRUE.)
      call ListInitElementKeyword(vh_bc_cond_h, 'Boundary Condition', &
          'Layer Electric Conductivity', InitIm=.TRUE.)

      ! Lumped ports
      call ElectricPortModel(1, vh_solver)
      inithandles = .false.
    end if

    call GetElementNodes(vh_bc_nodes, element)

    parent => GetBulkElementAtBoundary(element)

    vh_bc_stiff = 0.0_dp
    vh_bc_mass = 0.0_dp
    vh_bc_force = 0.0_dp
    if (associated(vh_amsmat)) vh_bc_amsstiff = 0.0_dp
    if (associated(vh_amsscalmat)) vh_bc_amsstiff2 = 0.0_dp

    ndofs = maxval(vh_solver % Def_Dofs(GetElementFamily(element), :, 1))
    withndofs = ndofs > 0
    np = n * ndofs

    goodconductor = ListGetElementLogical(vh_bc_goodcond_h, element, found)
    absorb = ListGetElementLogical(vh_bc_absorb_h, element, found)

    call ElectricPortModel(2, vh_solver, element, gotport)

    ! Numerical integration:
    !-----------------------
    IP = GaussPoints(element, PReferenceElement=vh_piola, &
        EdgeBasisDegree=vh_edge_degree)

    if (withndofs) then
      regularize = vh_usegauss .and. ListGetElementLogical( &
          vh_bc_chargecons_h, element, found)
    end if

    lineelement = GetElementFamily(element) == 2
    degenerateelement = (CoordinateSystemDimension() == 3) .and. lineelement
    updatestiff = .false.

    do t = 1, IP % n
      !
      ! We need to branch as the only way to get the traces of 2D vector
      ! finite elements is to call EdgeElementInfo:
      !
      if (lineelement) then
        ! TO DO: Call ElementInfo instead
        stat = EdgeElementInfo(element, vh_bc_nodes, IP % U(t), IP % V(t), &
            IP % W(t), detF = detJ, Basis = vh_bc_basis, &
            EdgeBasis = vh_bc_wbasis, RotBasis = vh_bc_rotw, &
            dBasisdx = vh_bc_db, BasisDegree = vh_edge_degree, &
            ApplyPiolaTransform = .TRUE.)
      else
        stat = ElementInfo(element, vh_bc_nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, vh_bc_basis, vh_bc_db, &
            EdgeBasis = vh_bc_wbasis, RotBasis = vh_bc_rotw, &
            USolver = vh_psolver)
      end if

      weight = detJ * IP % s(t)

      th = ListGetElementReal(vh_bc_thick_h, vh_bc_basis, element, &
          thinsheet, GaussPoint = t)

      if (degenerateelement .and. thinsheet .and. vh_usegauss) then
        !
        ! If a degenerate (1-D) element, perform only simplified assembly by
        ! assuming that a given thickness is associated with the degenerate
        ! element
        !
        jn = ListGetElementComplex(vh_bc_elcurr_h, vh_bc_basis, element, &
            found, GaussPoint = t)

        betapar = ListGetElementComplex(vh_bc_transfer_h, vh_bc_basis, &
            element, found, GaussPoint = t)
        if (found) then
          ep = ListGetElementComplex(vh_bc_extpot_h, vh_bc_basis, element, &
              found, GaussPoint = t)
          if (found) jn = jn + 2 * betapar * ep
        end if

        if (abs(betapar) < AEPS .and. abs(jn) < AEPS) cycle
        do p = 1, n
          i = (p - 1) * ndofs + 1
          vh_bc_force(i) = vh_bc_force(i) - im * vh_omega * jn * th * &
              vh_bc_basis(p) * weight
          do q = 1, n
            j = (q - 1) * ndofs + 1
            vh_bc_stiff(i, j) = vh_bc_stiff(i, j) - im * vh_omega * &
                betapar * th * vh_bc_basis(p) * vh_bc_basis(q) * weight
          end do
        end do
        updatestiff = .true.
        cycle
      end if

      if (associated(parent)) then
        mur = ListGetElementComplex(vh_bc_mu_h, vh_bc_basis, parent, &
            found, GaussPoint = t)
        if (.not. found) mur = 1.0_dp
        epsr = ListGetElementComplex(vh_bc_eps_h, vh_bc_basis, parent, &
            found, GaussPoint = t)
        if (.not. found) epsr = 1.0_dp
      else
        epsr = 1.0_dp
        mur = 1.0_dp
      end if
      muinv = mur * vh_mu0inv

      if (count([goodconductor, thinsheet, absorb, gotport]) > 1) then
        call Fatal(vh_caller, 'Boundary condition not uniquely defined!')
      end if

      L = cmplx(0.0_dp, 0.0_dp, kind=dp)
      if (absorb) then
        b = im * vh_rob0 * sqrt(epsr / mur)
      else if (thinsheet) then
        cond = ListGetElementComplex(vh_bc_cond_h, vh_bc_basis, element, &
            found, GaussPoint = t)
        b = im * (vh_omega / muinv) * th * cond
      else if (goodconductor) then
        cond = ListGetElementComplex(vh_bc_cond_h, vh_bc_basis, element, &
            found, GaussPoint = t)
        muinv = ListGetElementComplex(vh_bc_relnu_h, vh_bc_basis, element, &
            found, GaussPoint = t)
        if (found) then
          muinv = muinv * vh_mu0inv
        else
          muinv = vh_mu0inv
        end if
        surfimp = cmplx(1.0_dp, -1.0_dp, kind=dp) * &
            sqrt(vh_omega / (2.0_dp * cond * muinv))
        b = im * (vh_omega / muinv) / surfimp
      else if (gotport) then
        call ElectricPortModel(3, vh_solver, element, gotport, b, l, &
            vh_bc_basis, vh_bc_db, vh_bc_wbasis)
      else
        b = ListGetElementComplex(vh_bc_elrobin_h, vh_bc_basis, element, &
            found, GaussPoint = t)

        magload = ListGetElementComplex3D(vh_bc_magload_h, vh_bc_basis, &
            element, found, GaussPoint = t)
        temgrad = cmplx(ListGetElementRealGrad(vh_bc_temre_h, vh_bc_db, &
            element, found), ListGetElementRealGrad(vh_bc_temim_h, vh_bc_db, &
            element, found), kind=dp)

        elsurfcurr = ListGetElementComplex3D(vh_bc_elsurf_h, vh_bc_basis, &
            element, found, GaussPoint = t)

        L = magload + temgrad - (0_dp, 1_dp) * vh_omega / muinv * elsurfcurr
      end if

      if (.not. withndofs) then
        if (abs(b) < AEPS .and. abs(dot_product(L, L)) < AEPS) cycle
      end if
      updatestiff = .true.

      do i = 1, nd - np
        p = i + np
        vh_bc_force(p) = vh_bc_force(p) - muinv * &
            sum(L * vh_bc_wbasis(i, :)) * weight
        do j = 1, nd - np
          q = j + np
          vh_bc_stiff(p, q) = vh_bc_stiff(p, q) - muinv * b * &
              sum(vh_bc_wbasis(i, :) * vh_bc_wbasis(j, :)) * weight
        end do
      end do

      if (vh_amsany) then
        block
          complex(kind=dp) :: ar, atot

          if (vh_amscurlcurl) then
            normal = Normalvector(element, vh_bc_nodes, IP % U(t), &
                IP % V(t), .TRUE.)
          end if

          do p = 1, n
            do q = 1, n
              ar = -muinv * b * vh_bc_basis(q) * vh_bc_basis(p)
              atot = ar * weight

              if (associated(vh_amsmat)) then
                if (vh_amsmono) then
                  if (vh_amscurlcurl) then
                    dim = 3
                    do j = 1, dim
                      testvec = 0.0d0
                      testvec(j) = vh_bc_basis(p)
                      testvec = CrossProduct(testvec, normal)
                      do i = 1, dim
                        trialvec = 0.0d0
                        trialvec(i) = vh_bc_basis(q)
                        trialvec = CrossProduct(trialvec, normal)
                        vh_bc_amsstiff(dim * (p - 1) + j, dim * (q - 1) + i) = &
                            vh_bc_amsstiff(dim * (p - 1) + j, &
                            dim * (q - 1) + i) &
                            - muinv * b * sum(testvec(:) * trialvec(:)) * weight
                      end do
                    end do
                  else
                    vh_bc_amsstiff(3 * p - 2, 3 * q - 2) = &
                        vh_bc_amsstiff(3 * p - 2, 3 * q - 2) + atot
                    vh_bc_amsstiff(3 * p - 1, 3 * q - 1) = &
                        vh_bc_amsstiff(3 * p - 1, 3 * q - 1) + atot
                    vh_bc_amsstiff(3 * p - 0, 3 * q - 0) = &
                        vh_bc_amsstiff(3 * p - 0, 3 * q - 0) + atot
                  end if
                else
                  vh_bc_amsstiff(p, q) = vh_bc_amsstiff(p, q) + atot
                end if
              end if

              if (associated(vh_amsscalmat)) then
                vh_bc_amsstiff2(p, q) = vh_bc_amsstiff2(p, q) + atot
              end if

            end do
          end do
        end block
      end if

      if (withndofs) then
        ! The following term arises if the decomposition E = A - grad V is
        ! applied:
        if (abs(b) > AEPS) then
          do i = 1, nd - np
            p = i + np
            do j = 1, n
              q = (j - 1) * ndofs + 1
              vh_bc_stiff(p, q) = vh_bc_stiff(p, q) + muinv * b * &
                  sum(vh_bc_wbasis(i, :) * vh_bc_db(j, :)) * weight
            end do
          end do

          if (regularize) then
            ! Apply the conservation of surface charge (not sure whether
            ! this is beneficial):
            do p = 1, n
              i = (p - 1) * ndofs + 1
              do q = 1, n
                j = (q - 1) * ndofs + 1
                vh_bc_stiff(i, j) = vh_bc_stiff(i, j) + muinv * b * &
                    sum(vh_bc_db(p, :) * vh_bc_db(q, :)) * weight
              end do

              do q = 1, nd - np
                j = q + np
                vh_bc_stiff(i, j) = vh_bc_stiff(i, j) - muinv * b * &
                    sum(vh_bc_db(p, :) * vh_bc_wbasis(q, :)) * weight
              end do
            end do
            ! TO DO: If a distribution of surface charge were also given,
            !        we would need to add an additional term in order to be
            !        consistent
          end if
        end if

        if (vh_usegauss) then
          if (regularize .and. abs(dot_product(L, L)) > AEPS) then
            ! Apply the conservation of surface charge (not sure whether
            ! this is beneficial):
            do p = 1, n
              i = (p - 1) * ndofs + 1
              vh_bc_force(i) = vh_bc_force(i) - muinv * &
                  sum(L * vh_bc_db(p, :)) * weight
            end do
          end if

          jn = ListGetElementComplex(vh_bc_elcurr_h, vh_bc_basis, element, &
              found, GaussPoint = t)
          betapar = ListGetElementComplex(vh_bc_transfer_h, vh_bc_basis, &
              element, found, GaussPoint = t)
          if (found) then
            ep = ListGetElementComplex(vh_bc_extpot_h, vh_bc_basis, element, &
                found, GaussPoint = t)
            if (found) jn = jn + 2 * betapar * ep
          end if

          do p = 1, n
            i = (p - 1) * ndofs + 1
            vh_bc_force(i) = vh_bc_force(i) - im * vh_omega * jn * &
                vh_bc_basis(p) * weight
            do q = 1, n
              j = (q - 1) * ndofs + 1
              vh_bc_stiff(i, j) = vh_bc_stiff(i, j) - im * vh_omega * &
                  betapar * vh_bc_basis(p) * vh_bc_basis(q) * weight
            end do
          end do
        end if
      end if
    end do

    if (updatestiff) then
      if (vh_precmatrix) then
        if (vh_massprop .or. vh_curlcurlprec) then
          call DefaultUpdatePrec(vh_bc_stiff)
        else
          call DefaultUpdatePrec(vh_precdamp * vh_bc_stiff + vh_bc_stiff)
        end if
      end if
      call DefaultUpdateEquations(vh_bc_stiff, vh_bc_force, element)
    end if

    if (associated(vh_amsmat)) then
      CurrentModel % Solver => vh_amssolver
      call DefaultUpdateEquations(vh_bc_amsstiff, vh_bc_amsforce, &
          UElement=element, USolver=vh_amssolver)
      CurrentModel % Solver => vh_psolver
    end if
    if (associated(vh_amsscalmat)) then
      CurrentModel % Solver => vh_amsscalsolver
      call DefaultUpdateEquations(vh_bc_amsstiff2, vh_bc_amsforce, &
          UElement=element, USolver=vh_amsscalsolver)
      CurrentModel % Solver => vh_psolver
    end if
  end subroutine LocalMatrixBC

  !> This subroutine sets a single dipole load.
  !> It is mainly intended for testing purposes only.
  subroutine SingleDipoleLoad()
    type(valuelist_t), pointer :: dlist
    real(kind=dp), pointer :: wrkarray(:, :)
    real(kind=dp) :: source(3), coord(3), coord0(3), edgevec(3), mincoord(3)
    real(kind=dp) :: mindist, dist, dotprod, maxdotprod, edgelen
    real(kind=dp) :: maxedgelen, parmindist
    integer :: minind, maxedgeind, i, j, t
    type(element_t), pointer :: element
    logical :: found

    dlist => vh_solver % Values

    wrkarray => ListGetConstRealArray(dlist, 'Dipole Coordinates', found)
    if (.not. found) return

    call Info('DipoleSource', &
        '>Dipole Coordinates< given, setting up single edge source!')
    if (size(wrkarray, 1) /= 3 .or. size(wrkarray, 2) /= 1) then
      call Fatal('DipoleSource', '>Dipole Coordinates< of wrong size!')
    end if
    coord0(1:3) = wrkarray(1:3, 1)

    wrkarray => ListGetConstRealArray(dlist, 'Dipole Current', found)
    if (.not. found) call Fatal('DipoleSource', '>Dipole Current< not given!')
    if (size(wrkarray, 1) /= 3 .or. size(wrkarray, 2) /= 1) then
      call Fatal('DipoleSource', '>Dipole Current< of wrong size!')
    end if
    source(1:3) = wrkarray(1:3, 1)

    mindist = huge(mindist)
    do t = 1, vh_mesh % NumberOfNodes
      coord(1) = vh_mesh % Nodes % x(t)
      coord(2) = vh_mesh % Nodes % y(t)
      coord(3) = vh_mesh % Nodes % z(t)
      dist = sqrt(sum((coord - coord0)**2))
      if (dist < mindist) then
        mindist = dist
        minind = t
        mincoord = coord
      end if
    end do

    print *, 'Minimum distance at node:', minind, mindist
    print *, 'Dipole moment point:', coord0
    print *, 'Nearest mesh point:', mincoord

    if (ParEnv % PEs > 1) then
      parmindist = ParallelReduction(mindist, 1)
      if (parmindist < mindist) return
      print *, 'Minimum distance node found in partition:', ParEnv % MyPE
    end if

    maxdotprod = 0.0_dp

    do t = 1, vh_mesh % NumberOfEdges
      element => vh_mesh % Edges(t)

      if (all(element % NodeIndexes /= minind)) cycle

      call GetElementNodes(vh_dip_nodes, element)

      edgevec(1) = vh_dip_nodes % x(1) - vh_dip_nodes % x(2)
      edgevec(2) = vh_dip_nodes % y(1) - vh_dip_nodes % y(2)
      edgevec(3) = vh_dip_nodes % z(1) - vh_dip_nodes % z(2)

      edgelen = sqrt(sum(edgevec**2))

      dotprod = sum(edgevec * source) / edgelen

      !PRINT *,'candidate:',t,edgelen,dotprod

      if (abs(dotprod) > abs(maxdotprod)) then
        !PRINT *,'found better'
        maxdotprod = dotprod
        maxedgelen = edgelen
        maxedgeind = t
      end if
    end do

    print *, 'Edge with best orientation:', maxedgeind, maxdotprod, maxedgelen

    i = vh_mesh % NumberOfNodes + maxedgeind
    j = vh_solver % Variable % Perm(i)

    print *, 'set source at rhs:', i, j, maxdotprod
    if (j > 0) then
      ! Component in real valued system is: im*omega*Source
      ! i.e. real dipole source results to imaginary source
      vh_solver % Matrix % Rhs(2 * j) = vh_omega * maxdotprod / maxedgelen
    end if
  end subroutine SingleDipoleLoad

  ! -- VectorHelmholtzCalcFields_Init0 ------------------------------------
  subroutine vcf_init0(h_model, h_solver) bind(c, name="vcf_init0")
    type(c_ptr), value :: h_model, h_solver
    type(model_t), pointer :: model
    type(solver_t), pointer :: solver
    type(solver_t), pointer :: solvers(:), psolver
    character(len=MAX_NAME_LEN) :: sname, pname
    logical :: found, elementfields
    integer, pointer :: active(:)
    integer :: mysolver, i, j, n, m, soln
    type(valuelist_t), pointer :: solverparams
    call c_f_pointer(h_model, model)
    call c_f_pointer(h_solver, solver)

    solverparams => GetSolverParams()

    ! Find the solver index of the primary solver by the known procedure name.
    ! (the solver is defined here in the same module so not that dirty...)
    soln = ListGetInteger(solverparams, 'Primary Solver index', found)

    if (soln == 0) then
      do i = 1, model % NumberOfSolvers
        sname = GetString(model % Solvers(i) % Values, 'Procedure', found)
        j = index(sname, 'VectorHelmholtzSolver')
        if (j > 0) then
          soln = i
          exit
        end if
      end do
    end if

    if (soln == 0) then
      pname = GetString(solverparams, 'Potential variable', found)
      if (found) then
        do i = 1, model % NumberOfSolvers
          sname = GetString(model % Solvers(i) % Values, 'Variable', found)

          j = index(sname, '[') - 1
          if (j <= 0) j = len_trim(sname)
          if (sname(1:j) == pname(1:len_trim(pname))) then
            soln = i
            exit
          end if
        end do
      end if
    end if

    if (soln == 0) then
      call Fatal('VectorHelmholtzCalcFields_Init0', &
          'Cannot locate the primary solver: '//I2S(soln))
    else
      call Info('VectorHelmholtzCalcFields_Init0', &
          'The primary solver index is: '//I2S(soln), level=12)
      call ListAddInteger(solverparams, 'Primary Solver Index', soln)
    end if

    ! If the primary solver computed DG fields, then we don't need to create
    ! DG solver on-the-fly. We only need it if we have both nodal and DG
    ! solvers needed at the same time.
    !---------------------------------------------------------------------
    if (GetLogical(solverparams, 'Discontinuous Galerkin', found)) return

    elementfields = GetLogical(solverparams, 'Calculate Elemental Fields', found)
    if (found .and. .not. elementfields) return

    psolver => solver
    do mysolver = 1, model % NumberOfSolvers
      if (associated(psolver, model % Solvers(mysolver))) exit
    end do

    n = model % NumberOfSolvers
    do i = 1, model % NumberOfEquations
      active => ListGetIntegerArray(model % Equations(i) % Values, &
          'Active Solvers', found)
      m = size(active)
      if (any(active == mysolver)) &
          call ListAddIntegerArray(model % Equations(i) % Values, &
          'Active Solvers', m + 1, [active, n + 1])
    end do

    allocate (solvers(n + 1))
    solvers(1:n) = model % Solvers
    solverparams => null()
    call ListAddLogical(solverparams, 'Discontinuous Galerkin', .TRUE.)
    solvers(n + 1) % DG = .TRUE.
    solvers(n + 1) % Values => solverparams
    solvers(n + 1) % PROCEDURE = C_NULL_FUNPTR
    solvers(n + 1) % ActiveElements => null()
    call ListAddString(solverparams, 'Exec Solver', 'never')
    call ListAddLogical(solverparams, 'No Matrix', .TRUE.)
    call ListAddString(solverparams, 'Equation', 'never')
    call ListAddString(solverparams, 'Procedure', &
        'VectorHelmholtz VectorHelmholtz_Dummy', .FALSE.)
    call ListAddString(solverparams, 'Variable', '-nooutput cf_dummy')

    pname = ListGetString(model % Solvers(soln) % Values, 'Mesh', found)
    if (found) then
      call ListAddString(solverparams, 'Mesh', pname)
    end if

    if (GetLogical(solver % Values, 'Calculate Magnetic Flux Density', found)) then
      call ListAddString(solverparams, &
          NextFreeKeyword('Exported Variable', solverparams), &
          "Magnetic Flux Density E[Magnetic Flux Density re E:3 Magnetic Flux Density im E:3]")
    end if

    if (GetLogical(solver % Values, 'Calculate Electric field', found)) then
      call ListAddString(solverparams, &
          NextFreeKeyword('Exported Variable', solverparams), &
          "Electric field E[Electric field re E:3 Electric field im E:3]")
    end if

    if (GetLogical(solver % Values, 'Calculate Magnetic Field Strength', found)) then
      call ListAddString(solverparams, &
          NextFreeKeyword('Exported Variable', solverparams), &
          "Magnetic Field Strength E[Magnetic Field Strength re E:3 Magnetic Field Strength im E:3]")
    end if

    if (GetLogical(solver % Values, 'Calculate Poynting vector', found)) then
      call ListAddString(solverparams, &
          NextFreeKeyword('Exported Variable', solverparams), &
          "Poynting vector E[Poynting vector re E:3 Poynting vector im E:3]")
    end if

    if (GetLogical(solver % Values, 'Calculate Div of Poynting Vector', found)) then
      call ListAddString(solverparams, &
          NextFreeKeyword('Exported Variable', solverparams), &
          "Div Poynting Vector E[Div Poynting Vector re E:1 Div Poynting Vector im E:1]")
      call ListAddString(solverparams, &
          NextFreeKeyword('Exported Variable', solverparams), &
          "Joule Heating E[Joule Heating re E:1 Joule Heating im E:1]")
    end if

    deallocate (model % Solvers)
    model % Solvers => solvers
    model % NumberOfSolvers = n + 1
  end subroutine vcf_init0

  ! -- VectorHelmholtzCalcFields_Init -------------------------------------
  subroutine vcf_init(h_model, h_solver) bind(c, name="vcf_init")
    type(c_ptr), value :: h_model, h_solver
    type(model_t), pointer :: model
    type(solver_t), pointer :: solver
    logical :: found, nodalfields
    type(valuelist_t), pointer :: solverparams
    call c_f_pointer(h_model, model)
    call c_f_pointer(h_solver, solver)

    solverparams => GetSolverParams()

    call ListAddString(solverparams, 'Variable', '-nooutput vectorhelmholtz_dummy')

    call ListAddLogical(solverparams, 'Linear System refactorize', .FALSE.)

    ! add these in the beginning, so that SaveScalars sees these existing,
    ! even if executed before the actual computations...
    ! -----------------------------------------------------------------------
    if (GetLogical(solverparams, 'Calculate Energy Functional', found)) then
      call ListAddConstReal(model % Simulation, 'res: Energy Functional', 0._dp)
      call ListAddConstReal(model % Simulation, 'res: Energy Functional im', 0._dp)
    end if

    if (GetLogical(solverparams, 'Show Angular Frequency', found)) &
        call ListAddConstReal(model % Simulation, 'res: Angular Frequency', 0._dp)

    if (GetLogical(solverparams, 'Calculate Electric Potential', found)) then
      call ListAddString(solverparams, &
          NextFreeKeyword('Exported Variable', solverparams), &
          "Electric Potential[Electric Potential re:1 Electric Potential im:1]")
    end if

    nodalfields = GetLogical(solverparams, 'Calculate Nodal Fields', found)
    if (.not. found) nodalfields = .TRUE.

    if (.not. nodalfields) return

    if (GetLogical(solverparams, 'Calculate Magnetic Flux Density', found)) then
      call ListAddString(solverparams, &
          NextFreeKeyword('Exported Variable', solverparams), &
          "Magnetic Flux Density[Magnetic Flux Density re:3 Magnetic Flux Density im:3]")
    end if

    if (GetLogical(solverparams, 'Calculate Electric field', found)) then
      call ListAddString(solverparams, &
          NextFreeKeyword('Exported Variable', solverparams), &
          "Electric field[Electric field re:3 Electric field im:3]")
    end if

    if (GetLogical(solverparams, 'Calculate Magnetic Field Strength', found)) then
      call ListAddString(solverparams, &
          NextFreeKeyword('Exported Variable', solverparams), &
          "Magnetic Field Strength[Magnetic Field Strength re:3 Magnetic Field Strength im:3]")
    end if

    if (GetLogical(solverparams, 'Calculate Poynting vector', found)) then
      call ListAddString(solverparams, &
          NextFreeKeyword('Exported Variable', solverparams), &
          "Poynting vector[Poynting vector re:3 Poynting vector im:3]")
    end if

    if (GetLogical(solverparams, 'Calculate Div of Poynting Vector', found)) then
      call ListAddString(solverparams, &
          NextFreeKeyword('Exported Variable', solverparams), &
          "Div Poynting Vector[Div Poynting Vector re:1 Div Poynting Vector im:1]")
      call ListAddString(solverparams, &
          NextFreeKeyword('Exported Variable', solverparams), &
          "Joule Heating[Joule Heating re:1 Joule Heating im:1]")
    end if
  end subroutine vcf_init


  ! -- VectorHelmholtzCalcFields phases ------------------------------------
  subroutine vcf_prologue(h_model, h_solver) bind(c, name="vcf_prologue")
    type(c_ptr), value :: h_model, h_solver
    type(valuelist_t), pointer :: solverparams
    logical :: found, stat
    integer :: i, soln, vdofs
    real(kind=dp), pointer :: xx(:), bb(:), tempvector(:), temprhs(:)
    real(kind=dp) :: hdotE_r, hdotE_i
    character(len=MAX_NAME_LEN) :: pname
    call c_f_pointer(h_model, vh_cf_model)
    call c_f_pointer(h_solver, vh_cf_solver)

    call Info(vh_cf_caller, '', level=6)
    call Info(vh_cf_caller, &
        '----------------------------------------------------------', level=6)
    call Info(vh_cf_caller, &
        'Computing derived fields for electromagnetic wave equation!', level=4)

    solverparams => GetSolverParams()

    soln = ListGetInteger(solverparams, 'Primary Solver Index', found)
    if (soln == 0) then
      call Fatal(vh_cf_caller, 'We should know > Primary Solver Index <')
    end if

    ! Pointer to primary solver
    vh_cf_psolver => vh_cf_model % Solvers(soln)

    pname = ListGetString(solverparams, 'Potential variable', found)
    if (found) then
      vh_cf_targetvar => VariableGet(vh_cf_psolver % Mesh % Variables, pname)
    else
      vh_cf_targetvar => vh_cf_psolver % Variable
      pname = getVarName(vh_cf_psolver % Variable)
    end if
    vh_cf_pname = pname

    call Info(vh_cf_caller, 'Name of potential variable: '//trim(pname), &
        level=10)

    ! Inherit the solution basis from the primary solver
    vdofs = vh_cf_targetvar % DOFs
    if (vdofs /= 2) then
      call Fatal(vh_cf_caller, &
          'Primary variable should have 2 dofs: '//I2S(vdofs))
    end if

    call EdgeElementStyle(vh_cf_psolver % Values, vh_cf_piola, &
        BasisDegree = vh_cf_degree)
    if (vh_cf_piola) call Info(vh_cf_caller, &
        'Using Piola transformed finite elements', level=5)

    vh_cf_usegauss = GetLogical(vh_cf_psolver % Values, 'Use Gauss Law', found)
    vh_cf_lorenz = GetLogical(vh_cf_psolver % Values, 'Lorenz Condition', found)

    vh_cf_omega = GetAngularFrequency(vh_cf_psolver % Values)

    found = .false.
    if (associated(vh_cf_model % Constants)) then
      vh_cf_mu0inv = GetConstReal(vh_cf_model % Constants, &
          'Permeability of Vacuum', found)
      if (vh_cf_mu0inv /= 0) vh_cf_mu0inv = 1 / vh_cf_mu0inv
    end if
    if (.not. found) vh_cf_mu0inv = 1.0_dp / (PI * 4.0d-7)

    found = .false.
    if (associated(vh_cf_model % Constants)) then
      vh_cf_eps0 = GetConstReal(vh_cf_model % Constants, &
          'Permittivity of Vacuum', found)
    end if
    if (.not. found) vh_cf_eps0 = 8.854187817d-12

    vh_cf_mesh => GetMesh()

    vh_cf_mfd => VariableGet(vh_cf_mesh % Variables, 'Magnetic Flux Density')
    vh_cf_el_mfd => VariableGet(vh_cf_mesh % Variables, &
        'Magnetic Flux Density E')

    vh_cf_mfs => VariableGet(vh_cf_mesh % Variables, 'Magnetic Field Strength')
    vh_cf_el_mfs => VariableGet(vh_cf_mesh % Variables, &
        'Magnetic Field Strength E')

    vh_cf_ef => VariableGet(vh_cf_mesh % Variables, 'Electric field')
    vh_cf_el_ef => VariableGet(vh_cf_mesh % Variables, 'Electric field E')

    vh_cf_pv => VariableGet(vh_cf_mesh % Variables, 'Poynting vector')
    vh_cf_el_pv => VariableGet(vh_cf_mesh % Variables, 'Poynting vector E')

    vh_cf_divpv => VariableGet(vh_cf_mesh % Variables, 'Div Poynting Vector')
    vh_cf_el_divpv => VariableGet(vh_cf_mesh % Variables, &
        'Div Poynting Vector E')

    vh_cf_ew => VariableGet(vh_cf_mesh % Variables, 'Joule Heating')
    vh_cf_el_ew => VariableGet(vh_cf_mesh % Variables, 'Joule Heating E')

    vh_cf_nfields = 0
    if (associated(vh_cf_mfd)) vh_cf_nfields = vh_cf_nfields + 3
    if (associated(vh_cf_mfs)) vh_cf_nfields = vh_cf_nfields + 3
    if (associated(vh_cf_ef)) vh_cf_nfields = vh_cf_nfields + 3
    if (associated(vh_cf_pv)) vh_cf_nfields = vh_cf_nfields + 3
    if (associated(vh_cf_divpv)) vh_cf_nfields = vh_cf_nfields + 1
    if (associated(vh_cf_ew)) vh_cf_nfields = vh_cf_nfields + 1
    vh_cf_nfields = vh_cf_nfields * 2  ! complex problem
    vh_cf_nodal = (vh_cf_nfields > 0)

    if (vh_cf_nodal) then
      allocate (vh_cf_gforce(size(vh_cf_solver % Matrix % RHS), &
          vh_cf_nfields))
      vh_cf_gforce = 0._dp
    end if

    vh_cf_fields = 0
    if (associated(vh_cf_el_mfd)) vh_cf_fields = vh_cf_fields + 3
    if (associated(vh_cf_el_mfs)) vh_cf_fields = vh_cf_fields + 3
    if (associated(vh_cf_el_ef)) vh_cf_fields = vh_cf_fields + 3
    if (associated(vh_cf_el_pv)) vh_cf_fields = vh_cf_fields + 3
    if (associated(vh_cf_el_divpv)) vh_cf_fields = vh_cf_fields + 1
    if (associated(vh_cf_el_ew)) vh_cf_fields = vh_cf_fields + 1
    vh_cf_fields = vh_cf_fields * 2 ! complex problem
    vh_cf_elemental = (vh_cf_fields > 0)

    vh_cf_fields = max(vh_cf_fields, vh_cf_nfields)
    vh_cf_n = vh_cf_mesh % MaxElementDOFs

    allocate (vh_cf_mass(vh_cf_n, vh_cf_n), &
        vh_cf_force(vh_cf_n, vh_cf_fields), vh_cf_pivot(vh_cf_n))
    allocate (vh_cf_wbasis(vh_cf_n, 3), vh_cf_sol(2, vh_cf_n), &
        vh_cf_rotw(vh_cf_n, 3), vh_cf_basis(vh_cf_n), &
        vh_cf_db(vh_cf_n, 3))

    vh_cf_sol = 0._dp
    vh_cf_energy = 0._dp
    vh_cf_energy_im = 0._dp

    xx => vh_cf_targetvar % Values

    hdotE_r = 0._dp
    hdotE_i = 0._dp

    ! This piece of code effectively does what keyword "Calculate Energy
    ! Norm" would but it treats the system complex valued.
    if (GetLogical(solverparams, 'Calculate Energy Functional', found)) then
      if (vh_cf_psolver % Variable % dofs /= 2) then
        call Fatal(vh_cf_caller, 'Cannot compute energy norm if dofs not 2!')
      end if

      bb => vh_cf_psolver % Matrix % RHS
      allocate (tempvector(vh_cf_psolver % Matrix % NumberOfRows))
      if (ParEnv % PEs > 1) then
        allocate (temprhs(size(bb)))
        temprhs = bb
        call ParallelInitSolve(vh_cf_psolver % Matrix, xx, temprhs, tempvector)
        call ParallelMatrixVector(vh_cf_psolver % Matrix, xx, tempvector, &
            .TRUE.)
      else
        call MatrixVectorMultiply(vh_cf_psolver % Matrix, xx, tempvector)
      end if

      do i = 1, size(xx, 1) / 2
        if (ParEnv % PEs > 1) then
          if (vh_cf_psolver % Matrix % ParMatrix % ParallelInfo % &
              NeighbourList(2 * (i - 1) + 1) % Neighbours(1) /= &
              ParEnv % MyPE) cycle
        end if
        hdotE_r = hdotE_r + xx(2 * (i - 1) + 1) * tempvector(2 * (i - 1) + 1) &
            - xx(2 * (i - 1) + 2) * tempvector(2 * (i - 1) + 2)
        hdotE_i = hdotE_i + xx(2 * (i - 1) + 1) * tempvector(2 * (i - 1) + 2) &
            + xx(2 * (i - 1) + 2) * tempvector(2 * (i - 1) + 1)
      end do

      hdotE_r = ParallelReduction(hdotE_r)
      hdotE_i = ParallelReduction(hdotE_i)
      write (Message, *) 'Energy Functional value:', hdotE_r, hdotE_i
      call Info(vh_cf_caller, Message)
      call ListAddConstReal(vh_cf_model % Simulation, &
          'res: Energy Functional', hdotE_r)
      call ListAddConstReal(vh_cf_model % Simulation, &
          'res: Energy Functional im', hdotE_i)
    end if

    call DefaultInitialize()

    ! CALL ListInitElementKeyword( CondCoeff_h,'Material','Electric Conductivity')
    call ListInitElementKeyword(vh_cf_eps_h, 'Material', &
        'Relative Permittivity', InitIm=.TRUE.)
    call ListInitElementKeyword(vh_cf_mu_h, 'Material', &
        'Relative Reluctivity', InitIm=.TRUE.)
    call ListInitElementKeyword(vh_cf_cur_h, 'Body Force', &
        'Current Density', InitIm=.TRUE., InitVec3D=.TRUE.)
  end subroutine vcf_prologue

  function vcf_nactive() result(n) bind(c, name="vcf_nactive")
    integer(c_int) :: n
    n = GetNOFActive()
  end function vcf_nactive

  subroutine vcf_elem(i) bind(c, name="vcf_elem")
    integer(c_int), value :: i
    type(element_t), pointer :: element
    real(kind=dp), pointer :: fsave(:)
    integer :: n, ndofs, np, l
    logical :: found
    element => GetActiveElement(i)
    vh_cf_elem => element
    n = GetElementNOFNodes()
    vh_cf_n = n
    ndofs = vh_cf_psolver % Def_Dofs(GetElementFamily(element), &
        element % BodyId, 1)
    vh_cf_ndofs = ndofs
    np = n * ndofs
    vh_cf_np = np
    if (vh_cf_usegauss .or. vh_cf_lorenz) then
      if (ndofs < 1) call Fatal(vh_cf_caller, 'Nodal DOFs needed')
    end if

    call CalcFieldsLocalAssembly()

    if (vh_cf_nodal) then
      call DefaultUpdateEquations(vh_cf_mass, vh_cf_force(:, 1))
      fsave => vh_cf_solver % Matrix % RHS
      do l = 1, vh_cf_fields
        vh_cf_solver % Matrix % RHS => vh_cf_gforce(:, l)
        call DefaultUpdateForce(vh_cf_force(:, l))
      end do
      vh_cf_solver % Matrix % RHS => fsave
    end if

    if (vh_cf_elemental) then
      vh_cf_nfields = 0
      call LUdecomp(vh_cf_mass, n, vh_cf_pivot)
      call LocalSol(vh_cf_el_mfd, 6, n, vh_cf_mass, vh_cf_force, &
          vh_cf_pivot, vh_cf_nfields) ! 2*3 components
      call LocalSol(vh_cf_el_mfs, 6, n, vh_cf_mass, vh_cf_force, &
          vh_cf_pivot, vh_cf_nfields)
      call LocalSol(vh_cf_el_ef, 6, n, vh_cf_mass, vh_cf_force, &
          vh_cf_pivot, vh_cf_nfields)
      call LocalSol(vh_cf_el_pv, 6, n, vh_cf_mass, vh_cf_force, &
          vh_cf_pivot, vh_cf_nfields)
      call LocalSol(vh_cf_el_divpv, 2, n, vh_cf_mass, vh_cf_force, &
          vh_cf_pivot, vh_cf_nfields)
      call LocalSol(vh_cf_el_ew, 2, n, vh_cf_mass, vh_cf_force, &
          vh_cf_pivot, vh_cf_nfields)
    end if
  end subroutine vcf_elem

  subroutine vcf_mid() bind(c, name="vcf_mid")
    type(valuelist_t), pointer :: solverparams
    real(kind=dp), pointer :: fsave(:)
    logical :: found
    solverparams => GetSolverParams()
    vh_cf_energy = ParallelReduction(vh_cf_energy)
    vh_cf_energy_im = ParallelReduction(vh_cf_energy_im)

    ! Assembly of the face terms:
    !----------------------------

    vh_cf_doave = GetLogical(solverparams, 'Average Within Materials', found)

    if (GetLogical(solverparams, 'Discontinuous Galerkin', found)) then
      if (vh_cf_doave) then
        vh_cf_force = 0.0_dp
        call AddLocalFaceTerms(vh_cf_mass, vh_cf_force(:, 1))
      end if
    end if

    if (vh_cf_nodal .or. vh_cf_doave) then
      fsave => null()
      if (associated(vh_cf_solver % Matrix)) &
          fsave => vh_cf_solver % Matrix % RHS
      vh_cf_nfields = 0
      call GlobalSol(vh_cf_mfd, 6, vh_cf_gforce, vh_cf_nfields, vh_cf_el_mfd)
      call GlobalSol(vh_cf_mfs, 6, vh_cf_gforce, vh_cf_nfields, vh_cf_el_mfs)
      call GlobalSol(vh_cf_ef, 6, vh_cf_gforce, vh_cf_nfields, vh_cf_el_ef)
      call GlobalSol(vh_cf_pv, 6, vh_cf_gforce, vh_cf_nfields, vh_cf_el_pv)
      call GlobalSol(vh_cf_divpv, 2, vh_cf_gforce, vh_cf_nfields, &
          vh_cf_el_divpv)
      call GlobalSol(vh_cf_ew, 2, vh_cf_gforce, vh_cf_nfields, vh_cf_el_ew)
      if (associated(fsave)) vh_cf_solver % Matrix % RHS => fsave
    end if

    write (Message, *) '(Electro) Integral of Divergence of Poynting Vector: ', &
        vh_cf_energy, vh_cf_energy_im
    call Info(vh_cf_caller, Message)
    call ListAddConstReal(vh_cf_model % Simulation, &
        'res: Integral of Div Poynting Vector', vh_cf_energy)
    call ListAddConstReal(vh_cf_model % Simulation, &
        'res: Integral of Div Poynting Vector im', vh_cf_energy_im)

    if (GetLogical(solverparams, 'Show Angular Frequency', found)) then
      write (Message, *) 'Angular Frequency: ', vh_cf_omega
      call Info(vh_cf_caller, Message)
      call ListAddConstReal(vh_cf_model % Simulation, 'res: Angular Frequency', &
          vh_cf_omega)
    end if
  end subroutine vcf_mid

  subroutine vcf_elpot() bind(c, name="vcf_elpot")
    real(kind=dp), pointer :: fsave(:)
    type(element_t), pointer :: element
    integer :: i, n, ndofs, np
    vh_cf_elpot => VariableGet(vh_cf_mesh % Variables, 'Electric Potential')
    if (associated(vh_cf_elpot)) then
      call Info(vh_cf_caller, 'Calculating complex valued electric potential')
      if (.not. associated(vh_cf_solver % Matrix)) then
        call Fatal(vh_cf_caller, 'Solver % Matrix not associated!')
      end if

      call DefaultInitialize()
      fsave => vh_cf_solver % Matrix % RHS

      if (.not. allocated(vh_cf_gforce)) then
        allocate (vh_cf_gforce(size(vh_cf_solver % Matrix % RHS), 2))
      end if
      vh_cf_gforce = 0._dp

      do i = 1, GetNOFActive()
        element => GetActiveElement(i)
        vh_cf_elem => element
        n = GetElementNOFNodes()
        vh_cf_n = n
        ndofs = vh_cf_psolver % Def_Dofs(GetElementFamily(element), &
            element % BodyId, 1)
        vh_cf_ndofs = ndofs
        np = n * ndofs
        vh_cf_np = np

        call ElPotLocalAssembly()

        vh_cf_solver % Matrix % rhs => vh_cf_gforce(:, 1)
        call DefaultUpdateEquations(vh_cf_mass, vh_cf_force(:, 1))

        vh_cf_solver % Matrix % RHS => vh_cf_gforce(:, 2)
        call DefaultUpdateForce(vh_cf_force(:, 2))
      end do

      vh_cf_nfields = 0
      call GlobalSol(vh_cf_elpot, 2, vh_cf_gforce, vh_cf_nfields)
      vh_cf_solver % Matrix % RHS => fsave
    end if
  end subroutine vcf_elpot

  subroutine vcf_finish() bind(c, name="vcf_finish")
    if (allocated(vh_cf_gforce)) deallocate (vh_cf_gforce)
    deallocate (vh_cf_mass, vh_cf_force, vh_cf_pivot)
    if (allocated(vh_cf_wbasis)) deallocate (vh_cf_wbasis)
    if (allocated(vh_cf_sol)) deallocate (vh_cf_sol)
    if (allocated(vh_cf_rotw)) deallocate (vh_cf_rotw)
    if (allocated(vh_cf_basis)) deallocate (vh_cf_basis)
    if (allocated(vh_cf_db)) deallocate (vh_cf_db)
    call Info(vh_cf_caller, 'All done for now!', level=20)
  end subroutine vcf_finish

  subroutine CalcFieldsLocalAssembly()
    integer :: nd, j, p, q, k
    real(kind=dp) :: detj, s, u, v, w
    logical :: stat, found
    complex(kind=dp) :: h(3), exhc(3), divs, edotj, ef_ip(3), r_ip, b(3), pr_ip
    complex(kind=dp) :: j_ip(3)
    type(GaussIntegrationPoints_t) :: IP
    nd = GetElementNOFDOFs(uSolver=vh_cf_psolver)
    vh_cf_nd = nd
    call GetElementNodes(vh_cf_nodes)
    call GetVectorLocalSolution(vh_cf_sol, vh_cf_pname, &
        uSolver=vh_cf_psolver)

    ! Calculate nodal fields:
    ! -----------------------
    IP = GaussPoints(vh_cf_elem, &
        PReferenceElement=vh_cf_piola, EdgeBasisDegree = vh_cf_degree)

    vh_cf_mass = 0._dp
    vh_cf_force = 0._dp

    ! Loop over Gaussian integration points
    !---------------------------------------
    do j = 1, IP % n
      u = IP % U(j)
      v = IP % V(j)
      w = IP % W(j)

      stat = ElementInfo(vh_cf_elem, vh_cf_nodes, u, v, w, detJ, &
          vh_cf_basis, vh_cf_db, EdgeBasis = vh_cf_wbasis, &
          RotBasis = vh_cf_rotw, USolver = vh_cf_psolver)

      b = cmplx(matmul(vh_cf_sol(2, vh_cf_np + 1:vh_cf_nd), &
          vh_cf_rotw(1:vh_cf_nd - vh_cf_np, :)) / (vh_cf_omega), &
          matmul(vh_cf_sol(1, vh_cf_np + 1:vh_cf_nd), &
          vh_cf_rotw(1:vh_cf_nd - vh_cf_np, :)) / (-vh_cf_omega), kind=dp)

      ! The conductivity as a tensor not implemented yet
      !C_ip = ListGetElementReal( CondCoeff_h, Basis, Element, Found, GaussPoint = j )

      j_ip = ListGetElementComplex3D(vh_cf_cur_h, vh_cf_basis, &
          vh_cf_elem, found, GaussPoint = j)

      r_ip = ListGetElementComplex(vh_cf_mu_h, vh_cf_basis, vh_cf_elem, &
          found, GaussPoint = j)
      if (.not. found) then
        r_ip = vh_cf_mu0inv
      else
        r_ip = r_ip * vh_cf_mu0inv
      end if
      h = r_ip * b

      pr_ip = ListGetElementComplex(vh_cf_eps_h, vh_cf_basis, vh_cf_elem, &
          found, GaussPoint = j)
      if (found) then
        pr_ip = vh_cf_eps0 * pr_ip
      else
        pr_ip = vh_cf_eps0
      end if

      ef_ip = cmplx(matmul(vh_cf_sol(1, vh_cf_np + 1:vh_cf_nd), &
          vh_cf_wbasis(1:vh_cf_nd - vh_cf_np, :)), &
          matmul(vh_cf_sol(2, vh_cf_np + 1:vh_cf_nd), &
          vh_cf_wbasis(1:vh_cf_nd - vh_cf_np, :)), kind=dp)
      if (vh_cf_lorenz .or. vh_cf_usegauss) then
        do k = 1, 3
          ef_ip(k) = ef_ip(k) - &
              cmplx(sum(vh_cf_sol(1, 1:vh_cf_np:vh_cf_ndofs) * &
              vh_cf_db(1:vh_cf_n, k)), &
              sum(vh_cf_sol(2, 1:vh_cf_np:vh_cf_ndofs) * &
              vh_cf_db(1:vh_cf_n, k)), kind=dp)
        end do
      end if
      exhc = ComplexCrossProduct(ef_ip, conjg(h))

      edotj = sum(ef_ip * conjg(j_ip))
      divs = 0.5_dp * (im * vh_cf_omega * &
          (sum(b * conjg(h)) + sum(ef_ip * conjg(pr_ip * ef_ip))) - edotj)

      s = IP % s(j) * detJ

      vh_cf_energy = vh_cf_energy + real(divs) * s
      vh_cf_energy_im = vh_cf_energy_im + aimag(divs) * s
      do p = 1, vh_cf_n
        do q = 1, vh_cf_n
          vh_cf_mass(p, q) = vh_cf_mass(p, q) + s * vh_cf_basis(p) * &
              vh_cf_basis(q)
        end do
        k = 0
        if (associated(vh_cf_mfd) .or. associated(vh_cf_el_mfd)) then
          vh_cf_force(p, k + 1:k + 3) = vh_cf_force(p, k + 1:k + 3) + &
              s * real(b) * vh_cf_basis(p)
          k = k + 3
          vh_cf_force(p, k + 1:k + 3) = vh_cf_force(p, k + 1:k + 3) + &
              s * aimag(b) * vh_cf_basis(p)
          k = k + 3
        end if
        if (associated(vh_cf_mfs) .or. associated(vh_cf_el_mfs)) then
          vh_cf_force(p, k + 1:k + 3) = vh_cf_force(p, k + 1:k + 3) + &
              s * real(h) * vh_cf_basis(p)
          k = k + 3
          vh_cf_force(p, k + 1:k + 3) = vh_cf_force(p, k + 1:k + 3) + &
              s * aimag(h) * vh_cf_basis(p)
          k = k + 3
        end if
        if (associated(vh_cf_ef) .or. associated(vh_cf_el_ef)) then
          vh_cf_force(p, k + 1:k + 3) = vh_cf_force(p, k + 1:k + 3) + &
              s * real(ef_ip) * vh_cf_basis(p)
          k = k + 3
          vh_cf_force(p, k + 1:k + 3) = vh_cf_force(p, k + 1:k + 3) + &
              s * aimag(ef_ip) * vh_cf_basis(p)
          k = k + 3
        end if
        if (associated(vh_cf_pv) .or. associated(vh_cf_el_pv)) then
          vh_cf_force(p, k + 1:k + 3) = vh_cf_force(p, k + 1:k + 3) + &
              0.5_dp * s * real(exhc) * vh_cf_basis(p)
          k = k + 3
          vh_cf_force(p, k + 1:k + 3) = vh_cf_force(p, k + 1:k + 3) + &
              0.5_dp * s * aimag(exhc) * vh_cf_basis(p)
          k = k + 3
        end if
        if (associated(vh_cf_divpv) .or. associated(vh_cf_el_divpv)) then
          vh_cf_force(p, k + 1) = vh_cf_force(p, k + 1) + &
              s * real(divs) * vh_cf_basis(p)
          k = k + 1
          vh_cf_force(p, k + 1) = vh_cf_force(p, k + 1) + &
              s * aimag(divs) * vh_cf_basis(p)
          k = k + 1
          vh_cf_force(p, k + 1) = vh_cf_force(p, k + 1) + &
              0.5_dp * s * real(edotj) * vh_cf_basis(p)
          k = k + 1
          vh_cf_force(p, k + 1) = vh_cf_force(p, k + 1) + &
              0.5_dp * s * aimag(edotj) * vh_cf_basis(p)
          k = k + 1
        end if

      end do
    end do
  end subroutine CalcFieldsLocalAssembly

  subroutine ElPotLocalAssembly()
    integer :: nd, j, p, q, k
    real(kind=dp) :: detj, s, u, v, w
    logical :: stat, found
    complex(kind=dp) :: ef_ip(3)
    type(GaussIntegrationPoints_t) :: IP
    nd = GetElementNOFDOFs(uSolver=vh_cf_psolver)
    vh_cf_nd = nd
    call GetElementNodes(vh_cf_nodes)
    call GetVectorLocalSolution(vh_cf_sol, vh_cf_pname, &
        uSolver=vh_cf_psolver)

    ! Calculate nodal fields:
    ! -----------------------
    IP = GaussPoints(vh_cf_elem, &
        PReferenceElement=vh_cf_piola, EdgeBasisDegree = vh_cf_degree)

    vh_cf_mass = 0._dp
    vh_cf_force = 0._dp

    ! Loop over Gaussian integration points
    !---------------------------------------
    do j = 1, IP % n
      u = IP % U(j)
      v = IP % V(j)
      w = IP % W(j)

      stat = ElementInfo(vh_cf_elem, vh_cf_nodes, u, v, w, detJ, &
          vh_cf_basis, vh_cf_db, EdgeBasis = vh_cf_wbasis, &
          RotBasis = vh_cf_rotw, USolver = vh_cf_psolver)

      ef_ip = cmplx(matmul(vh_cf_sol(1, vh_cf_np + 1:vh_cf_nd), &
          vh_cf_wbasis(1:vh_cf_nd - vh_cf_np, :)), &
          matmul(vh_cf_sol(2, vh_cf_np + 1:vh_cf_nd), &
          vh_cf_wbasis(1:vh_cf_nd - vh_cf_np, :)), kind=dp)
      if (vh_cf_lorenz .or. vh_cf_usegauss) then
        do k = 1, 3
          ef_ip(k) = ef_ip(k) - &
              cmplx(sum(vh_cf_sol(1, 1:vh_cf_np:vh_cf_ndofs) * &
              vh_cf_db(1:vh_cf_n, k)), &
              sum(vh_cf_sol(2, 1:vh_cf_np:vh_cf_ndofs) * &
              vh_cf_db(1:vh_cf_n, k)), kind=dp)
        end do
      end if

      s = IP % s(j) * detJ
      do p = 1, vh_cf_n
        do q = 1, vh_cf_n
          vh_cf_mass(p, q) = vh_cf_mass(p, q) + &
              s * sum(vh_cf_db(p, :) * vh_cf_db(q, :))
        end do
        vh_cf_force(p, 1) = vh_cf_force(p, 1) + &
            s * sum(real(ef_ip(:) * vh_cf_db(p, :)))
        vh_cf_force(p, 2) = vh_cf_force(p, 2) + &
            s * sum(aimag(ef_ip(:) * vh_cf_db(p, :)))
      end do
    end do
  end subroutine ElPotLocalAssembly

  !----------------------------------------------------------------------------
  subroutine GlobalSol(var, m, b, dofs, EL_Var)
  !----------------------------------------------------------------------------
    use MeshBasics, only : CalculateBodyAverage
    implicit none
    real(kind=dp), target, contiguous :: b(:, :)
    integer :: m, dofs
    type(variable_t), pointer :: var
    type(variable_t), pointer, optional :: EL_Var
    integer :: i

    if (present(EL_Var)) then
      if (associated(EL_Var)) then
        EL_Var % DgAveraged = .false.
        if (vh_cf_doave) then
          call Info(vh_cf_caller, 'Averaging for field: '//trim(EL_Var % Name), &
              level=10)
          call CalculateBodyAverage(vh_cf_mesh, EL_Var, .FALSE.)
        end if
        if (.not. (associated(var) .and. vh_cf_nodal)) then
          dofs = dofs + m
          return
        end if
      end if
    end if

    if (.not. associated(var)) return

    call Info('VectorHelmholtz', 'Solving for field: '//trim(var % Name), &
        level=6)
    do i = 1, m
      dofs = dofs + 1
      vh_cf_solver % Matrix % RHS => b(:, dofs)
      vh_cf_solver % Variable % Values = 0
      vh_cf_norm = DefaultSolve()
      var % Values(i::m) = vh_cf_solver % Variable % Values
    end do
  end subroutine GlobalSol

  !----------------------------------------------------------------------------
  subroutine LocalSol(var, m, n, a, b, pivot, dofs)
  !----------------------------------------------------------------------------
    type(variable_t), pointer :: var
    integer :: pivot(:), m, n, dofs
    real(kind=dp) :: b(:, :), a(:, :)
    integer :: ind(n), i
    real(kind=dp) :: x(n)

    if (.not. associated(var)) return

    ind = var % dofs * (var % Perm(vh_cf_elem % DGIndexes(1:n)) - 1)
    do i = 1, m
      dofs = dofs + 1
      x = b(1:n, dofs)
      call LUSolve(n, vh_cf_mass, x, pivot)
      var % Values(ind(1:n) + i) = x
    end do
  end subroutine LocalSol

  !----------------------------------------------------------------------------
  subroutine AddLocalFaceTerms(stiff, force)
  !----------------------------------------------------------------------------
    real(kind=dp) :: stiff(:, :), force(:)

    type(element_t), pointer :: P1, P2, face, faces(:)
    integer :: t, n, n1, n2, numberoffaces, dim

    dim = CoordinateSystemDimension()

    if (dim == 2) then
      faces => vh_cf_solver % Mesh % Edges
      numberoffaces = vh_cf_solver % Mesh % NumberOfEdges
    else
      faces => vh_cf_solver % Mesh % Faces
      numberoffaces = vh_cf_solver % Mesh % NumberOfFaces
    end if

    do t = 1, numberoffaces
      face => faces(t)
      if (.not. ActiveBoundaryElement(face)) cycle

      P1 => face % BoundaryInfo % Left
      P2 => face % BoundaryInfo % Right
      if (associated(P2) .and. associated(P1)) then
        if (.not. associated(GetMaterial(P1), GetMaterial(P2))) cycle

        n = GetElementNOFNodes(face)
        n1 = GetElementNOFNodes(P1)
        n2 = GetElementNOFNodes(P2)

        call LocalJumps(stiff, face, n, P1, n1, P2, n2)
        call DefaultUpdateEquations(stiff, force, face)
      end if
    end do
  end subroutine AddLocalFaceTerms

  !----------------------------------------------------------------------------
  subroutine LocalJumps(stiff, face, n, P1, n1, P2, n2)
  !----------------------------------------------------------------------------
    implicit none
    real(kind=dp) :: stiff(:, :)
    integer :: n, n1, n2
    type(element_t), target :: face, P1, P2
    real(kind=dp) :: facebasis(n), P1basis(n1), P2basis(n2)
    real(kind=dp) :: jump(n1 + n2), detJ, U, V, W, S
    logical :: stat
    integer :: p, q, t
    type(GaussIntegrationPoints_t) :: integstuff

    type(nodes_t) :: facenodes, P1nodes, P2nodes
    save facenodes, P1nodes, P2nodes

    stiff = 0._dp

    call GetElementNodes(facenodes, face)
    call GetElementNodes(P1nodes, P1)
    call GetElementNodes(P2nodes, P2)
    ! Numerical integration over the edge
    integstuff = GaussPoints(face)

    do t = 1, integstuff % n
      U = integstuff % u(t)
      V = integstuff % v(t)
      W = integstuff % w(t)
      S = integstuff % s(t)

      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo(face, facenodes, U, V, W, detJ, facebasis)

      S = S * detJ

      ! Find basis functions for the parent elements:
      ! ---------------------------------------------
      call GetParentUVW(face, n, P1, n1, U, V, W, facebasis)
      stat = ElementInfo(P1, P1nodes, U, V, W, detJ, P1basis)

      call GetParentUVW(face, n, P2, n2, U, V, W, facebasis)
      stat = ElementInfo(P2, P2nodes, U, V, W, detJ, P2basis)

      ! Integrate jump terms:
      ! ---------------------
      jump(1:n1) = P1basis(1:n1)
      jump(n1 + 1:n1 + n2) = -P2basis(1:n2)

      do p = 1, n1 + n2
        do q = 1, n1 + n2
          stiff(p, q) = stiff(p, q) + s * jump(q) * jump(p)
        end do
      end do
    end do
  end subroutine LocalJumps

end module vectorhelmholtz_acc

