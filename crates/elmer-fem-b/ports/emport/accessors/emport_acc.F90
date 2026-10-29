! T022-A EMPort port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/EMPort.F90, compiled against the built .mod
! files of the target configure (A2 shipped mechanism, ADR-057 §3).
!
! Every wrapper takes only C-interoperable addresses (TYPE(C_PTR), VALUE
! scalars, explicit-shape REAL(C_DOUBLE) arrays). All Fortran-object
! traffic — TYPE(...) handles, OPTIONAL/POINTER dummies, assumed-shape
! descriptors, hidden lengths, array-pointer results, COMPLEX arithmetic —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module mirrors the original's SAVE state (entry-level ValueHandle_t
! handles, per-internal-procedure SAVE Nodes/IP/basis arrays, SAVE
! scalars, the EMPortPost SAVE set) with identical lifetimes.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module emport_acc
  use iso_c_binding
  use defutils
  use eleminfo
  use integration
  use elementgeometry
  use coordinatesystems
  use parallelutils, only: parallelcdot, parallelreduction, &
      parallelactive, partitionvector, parallelinitmatrix
  use generalutils, only: complexvariablevalues, complexvalues, i2s
  use messages
  use types
  use spariterglobals
  implicit none

  complex(c_double), parameter :: acc_im = (0.0_c_double, 1.0_c_double)

  ! -- shared entry state (mirrors EMPortSolver locals) -------------------
  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(mesh_t), pointer, save :: acc_mesh => null()
  type(valuelist_t), pointer, save :: acc_params => null()
  type(valuelist_t), pointer, save :: acc_bc => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(variable_t), pointer, save :: acc_emvar => null()
  type(valuehandle_t), save :: acc_nucoeff_h, acc_epscoeff_h
  logical, save :: acc_piola = .false.
  logical, save :: acc_usev = .false.
  integer, save :: acc_edge_deg = 1
  integer, save :: acc_maxport = 0
  integer, save :: acc_port_save = 1
  integer, save :: acc_mode = 1
  real(c_double), save :: acc_maxmu = 0.0_c_double
  real(c_double), save :: acc_maxeps = 0.0_c_double
  real(c_double), save :: acc_omega = 0.0_c_double
  real(c_double), save :: acc_mu0inv = 0.0_c_double
  real(c_double), save :: acc_eps0 = 0.0_c_double
  real(c_double), save :: acc_betasum = 0.0_c_double
  complex(c_double), save :: acc_beta = (0.0_c_double, 0.0_c_double)
  complex(c_double), save :: acc_power = (0.0_c_double, 0.0_c_double)
  complex(c_double), save :: acc_e2 = (0.0_c_double, 0.0_c_double)
  logical, save :: acc_output_z = .false.
  integer, allocatable, save :: acc_saveperm(:)
  complex(c_double), pointer, save :: acc_saveev(:, :) => null()
  real(c_double), allocatable, save :: acc_upart(:), acc_vpart(:)

  ! -- LocalMatrix (no potentials) SAVE set --------------------------------
  type(nodes_t), target, save :: acc_nodes_lm
  type(gaussintegrationpoints_t), save :: acc_ip_lm
  real(c_double), allocatable, save :: acc_basis_lm(:), acc_dbasis_lm(:, :), &
      acc_wbasis_lm(:, :), acc_curlw_lm(:, :)
  complex(c_double), allocatable, save :: acc_stiff_lm(:, :), &
      acc_mass_lm(:, :), acc_force_lm(:)

  ! -- LocalMatrix_With_Potentials SAVE set ---------------------------------
  type(nodes_t), target, save :: acc_nodes_lp
  type(gaussintegrationpoints_t), save :: acc_ip_lp
  real(c_double), allocatable, save :: acc_basis_lp(:), acc_dbasis_lp(:, :), &
      acc_wbasis_lp(:, :), acc_curlw_lp(:, :)
  complex(c_double), allocatable, save :: acc_stiff_lp(:, :), &
      acc_mass_lp(:, :), acc_force_lp(:)

  ! -- CalculatePortPower SAVE set (shared by both call sites, as in the
  ! -- original where it is a single internal procedure) -------------------
  type(nodes_t), target, save :: acc_nodes_pw
  type(gaussintegrationpoints_t), save :: acc_ip_pw
  real(c_double), allocatable, save :: acc_wbasis_pw(:, :), &
      acc_curlw_pw(:, :), acc_basis_pw(:), acc_dbasis_pw(:, :)
  real(c_double), allocatable, save :: acc_re_pw(:), acc_im_pw(:)

  ! -- EMPortPost SAVE set (mirrors its SAVE statement exactly) ------------
  type(solver_t), pointer, save :: acc_post => null()
  real(c_double), allocatable, target, save :: acc_postmass(:, :), &
      acc_postlf(:, :), acc_postgf(:, :)
  real(c_double), allocatable, save :: acc_postwb(:, :), acc_postcb(:, :), &
      acc_postbasis(:), acc_postdbasis(:, :), acc_postre(:), acc_postim(:)
  integer, allocatable, save :: acc_postperm(:)
  integer, save :: acc_postdofs = 0
  type(variable_t), pointer, save :: acc_postef => null()
  real(c_double), pointer, save :: acc_postfsave(:) => null()
  integer, pointer, save :: acc_postnperm(:) => null()
  logical, save :: acc_postalloc = .false.
  type(nodes_t), target, save :: acc_nodes_post
  type(gaussintegrationpoints_t), save :: acc_ip_post

  ! -- EMPortPotential solver SAVE set --------------------------------------
  type(model_t), pointer, save :: acc_pmodel => null()
  type(solver_t), pointer, save :: acc_psolver => null()
  type(mesh_t), pointer, save :: acc_pmesh => null()
  type(element_t), pointer, save :: acc_pelem => null()
  type(element_t), pointer, save :: acc_pparent => null()
  real(c_double), allocatable, save :: acc_pbasis(:), acc_pdbasis(:, :), &
      acc_pparentbasis(:)
  real(c_double), allocatable, save :: acc_pstiff(:, :), acc_pforce(:)
  type(nodes_t), target, save :: acc_nodes_pt
  type(gaussintegrationpoints_t), save :: acc_ip_pt
  type(valuehandle_t), save :: acc_psource_h, acc_peps_h
  real(c_double), save :: acc_peps0 = 0.0_c_double

contains

  ! -- EMPortSolver_Init0 ---------------------------------------------------
  subroutine emi0_begin(h_model, h_solver) bind(c, name="emi0_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine emi0_begin

  subroutine emi0_body() bind(c, name="emi0_body")
    type(valuelist_t), pointer :: params, bc, primaryparams
    logical :: found, piolaversion, secondfamily, secondorder
    logical :: usev, ground, unitvoltage
    real(c_double) :: c
    character(:), allocatable :: sname
    character(*), parameter :: caller = 'EMPortSolver_Init0'
    integer :: i, j, soln

    params => GetSolverParams()

    usev = GetLogical(params, 'Use Potential', found)

    call ListAddNewLogical(params, 'Linear System Complex', .true.)
    call ListAddNewInteger(params, 'Variable DOFs', 2)
    call ListAddNewLogical(params, 'Eigen Analysis', .true.)
    call ListAddNewInteger(params, 'Nonlinear System Max Iterations', 1)

    soln = ListGetInteger(params, 'Primary Solver index', found)
    if (soln == 0) then
      do i = 1, acc_model % numberofsolvers
        sname = GetString(acc_model % solvers(i) % values, 'Procedure', found)
        j = index(sname, 'VectorHelmholtzSolver')
        if (j > 0) then
          soln = i
          exit
        end if
      end do
    end if
    if (soln > 0) then
      call Info(caller, 'Copying edge element definitions from primary solver ' // &
          i2s(soln), level=8)
      primaryparams => acc_model % solvers(soln) % values
      call ListCompareAndCopy(primaryparams, params, 'Use Piola Transform')
      call ListCompareAndCopy(primaryparams, params, 'Quadratic Approximation')
      call ListCompareAndCopy(primaryparams, params, 'Second Kind Basis')
      call ListCompareAndCopy(primaryparams, params, 'Gradient Basis Functions')
    end if

    if (.not. ListCheckPresent(params, "Element")) then
      call EdgeElementStyle(params, piolaversion, secondfamily, secondorder, &
          check=.true.)

      if (usev) then
        if (piolaversion) then
          sname = "n:2 e:1 -quad_face b:2 -quad b:2 -brick b:3"
        else
          sname = "n:2 e:1"
        end if
      else
        if (secondfamily) then
          call Warn(caller, 'The formulation for Second Kind Basis seems numerically unstable')
        end if
        if (secondorder) then
          if (secondfamily) then
            sname = "n:1 e:3 -tri b:3 -tri_face b:3"
          else
            sname = "n:1 e:2 -tri b:2 -quad b:4 -brick b:6 -pyramid b:3 -prism b:2 -quad_face b:4 -tri_face b:2"
          end if
        else if (secondfamily) then
          sname = "n:1 e:2"
        else if (piolaversion) then
          sname = "n:1 e:1 -quad_face b:2 -quad b:2 -brick b:3"
        else
          sname = "n:1 e:1"
        end if
      end if
      call Info(caller, 'Setting element type: ' // trim(sname), level=5)
      call ListAddString(params, "Element", trim(sname))
    end if

    do i = 1, acc_model % numberofbcs
      bc => acc_model % bcs(i) % values
      ground = ListGetLogical(bc, "Port Ground", found)
      unitvoltage = ListGetLogical(bc, "Unit Voltage", found)
      if (ground .or. unitvoltage) then
        if (ground) then
          call Info(caller, 'Enforcing zero tangential electric field where "Port Ground" is True', &
              level=10)
          if (usev) then
            call Info(caller, 'Enforcing zero scalar potential where "Port Ground" is True', &
                level=10)
            c = 0.0_c_double
          end if
        else
          call Info(caller, 'Enforcing zero tangential electric field where "Unit Voltage" is True', &
              level=10)
          call Info(caller, 'Enforcing unit scalar potential where "Unit Voltage" is True', &
              level=10)
          c = 1.0_c_double
        end if

        if (usev) then
          call ListAddConstReal(bc, 'Eport re {n} 1', 0.0_c_double)
          call ListAddConstReal(bc, 'Eport re {n} 2', c)
          call ListAddConstReal(bc, 'Eport im {n} 1', 0.0_c_double)
          call ListAddConstReal(bc, 'Eport im {n} 2', 0.0_c_double)
        else
          call ListAddConstReal(bc, 'Eport re', 0.0_c_double)
          call ListAddConstReal(bc, 'Eport im', 0.0_c_double)
        end if
        call ListAddConstReal(bc, 'Eport re {e}', 0.0_c_double)
        call ListAddConstReal(bc, 'Eport im {e}', 0.0_c_double)
      end if
    end do

    call ListAddNewString(params, 'Variable', 'Eport[Eport re:1 Eport im:1]')
    call ListAddLogical(params, 'Linear System refactorize', .true.)

    call ListAddNewLogical(params, 'Skip Compute Nonlinear Change', .true.)
    call ListAddNewLogical(params, 'Skip Compute Steady State Change', .true.)
    call ListAddNewLogical(params, 'Variable Output', .false.)
    call ListAddNewLogical(params, 'post: Skip Compute Nonlinear Change', .true.)
    call ListAddNewLogical(params, 'post: Linear System Complex', .false.)
    call ListAddNewLogical(params, 'post: Variable Output', .false.)

    call Info(caller, 'Setting default sorting and normalization for eigenmodes!', level=7)
    call ListAddNewString(params, 'Eigen System Sorting', 'smallest real part')
    call ListAddNewLogical(params, 'Eigenvector Normalization by Power', .true.)
    if (.not. ListCheckPresent(params, 'Eigen System Shift')) &
        call ListAddNewLogical(params, 'Eigen System Shift Automatic', .true.)

    call ListAddNewLogical(params, 'Linear System Nullify Guess', .true.)
  end subroutine emi0_body

  ! -- EMPortSolver ---------------------------------------------------------
  subroutine ema_begin(h_model, h_solver) bind(c, name="ema_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    acc_betasum = 0.0_c_double
  end subroutine ema_begin

  subroutine ema_start_info() bind(c, name="ema_start_info")
    character(*), parameter :: caller = 'EMPortSolver'
    call Info(caller, '', level=6)
    call Info(caller, '-----------------------------------------------------', level=6)
    call Info(caller, 'Solving electromagnetic port equations over a surface', level=4)
    call Info(caller, '-----------------------------------------------------', level=6)
  end subroutine ema_start_info

  subroutine ema_init_handles() bind(c, name="ema_init_handles")
    call ListInitElementKeyword(acc_nucoeff_h, 'Material', 'Relative Reluctivity', initim=.true.)
    call ListInitElementKeyword(acc_epscoeff_h, 'Material', 'Relative Permittivity', initim=.true.)
  end subroutine ema_init_handles

  subroutine ema_check_cartesian() bind(c, name="ema_check_cartesian")
    character(*), parameter :: caller = 'EMPortSolver'
    if (CurrentCoordinateSystem() /= Cartesian) then
      call Fatal(caller, 'Implemented only for Cartesian problems!')
    end if
  end subroutine ema_check_cartesian

  subroutine ema_check_dofs() bind(c, name="ema_check_dofs")
    character(*), parameter :: caller = 'EMPortSolver'
    integer :: dofs
    dofs = acc_solver % variable % dofs
    if (dofs /= 2) then
      call Fatal(caller, 'Complex field, specify two DOFs instead of ' // i2s(dofs))
    end if
  end subroutine ema_check_dofs

  function ema_use_v() result(r) bind(c, name="ema_use_v")
    integer(c_int) :: r
    logical :: found
    acc_usev = GetLogical(acc_params_check(), 'Use Potential', found)
    r = merge(1, 0, acc_usev)
  end function ema_use_v

  function acc_params_check() result(p)
    type(valuelist_t), pointer :: p
    p => GetSolverParams()
    acc_params => p
  end function acc_params_check

  function ema_scan_ports() result(m) bind(c, name="ema_scan_ports")
    integer(c_int) :: m
    type(valuelist_t), pointer :: bc
    logical :: found
    integer :: i, j, maxport
    character(:), allocatable :: ptype
    maxport = 0
    do i = 1, acc_model % numberofbcs
      bc => acc_model % bcs(i) % values
      j = ListGetInteger(bc, "Port Index", found)
      if (.not. found) then
        ptype = ListGetString(bc, 'Port Type', found)
        if (ptype == 'eigenmode') then
          j = ListGetInteger(bc, "Constraint Mode", found)
          if (j > 0) call ListAddInteger(bc, "Port Index", j)
        end if
      end if
      if (j > 0) then
        call ListAddLogical(acc_model % bcs(i) % values, "Port Label " // i2s(j), .true.)
        maxport = max(maxport, j)
      end if
    end do
    acc_maxport = maxport
    m = maxport
  end function ema_scan_ports

  subroutine ema_fetch_emvar() bind(c, name="ema_fetch_emvar")
    acc_emvar => acc_solver % variable
  end subroutine ema_fetch_emvar

  ! MaxPort > 1 setup: separate matrices per port. Returns 1 when done.
  function ema_multiport_setup(nrows, ncols) result(r) bind(c, name="ema_multiport_setup")
    integer(c_int) :: nrows, ncols
    integer(c_int) :: r
    character(*), parameter :: caller = 'EMPortSolver'
    integer :: n, m
    r = 0
    if (acc_maxport > 1) then
      call Info(caller, 'Creating separate matrices for each of ' // i2s(acc_maxport) // &
          ' ports!', level=5)
      call FreeMatrix(acc_solver % matrix)
      acc_solver % matrix => null()
      allocate (acc_saveperm(size(acc_emvar % perm)))
      acc_saveperm = acc_emvar % perm
      n = size(acc_emvar % eigenvectors, 1)
      m = size(acc_emvar % eigenvectors, 2)
      allocate (acc_saveev(n, m))
      acc_saveev = 0.0_c_double
      nrows = n
      ncols = m
      r = 1
    end if
  end function ema_multiport_setup

  subroutine ema_edge_style() bind(c, name="ema_edge_style")
    call EdgeElementStyle(acc_params, acc_piola, basisdegree=acc_edge_deg)
  end subroutine ema_edge_style

  ! The original assigns EigenProblem = EigenOrHarmonicAnalysis(Solver);
  ! the value is never read afterwards. The call is kept so the call
  ! sequence matches, the dead result discarded.
  subroutine ema_eigen_query() bind(c, name="ema_eigen_query")
    logical :: eigenproblem
    eigenproblem = EigenOrHarmonicAnalysis(acc_solver)
  end subroutine ema_eigen_query

  subroutine ema_default_start() bind(c, name="ema_default_start")
    call DefaultStart()
  end subroutine ema_default_start

  subroutine ema_init_stuff() bind(c, name="ema_init_stuff")
    character(*), parameter :: caller = 'EMPortSolver'
    logical :: found
    found = .false.
    if (associated(acc_model % constants)) then
      if (ListCheckPresent(acc_model % constants, 'Permeability of Vacuum')) &
          acc_mu0inv = 1.0_c_double / GetConstReal(acc_model % constants, &
          'Permeability of Vacuum', found)
    end if
    if (.not. found) acc_mu0inv = 1.0_c_double / (3.14159265358979323846_c_double * 4.0d-7)
    found = .false.
    if (associated(acc_model % constants)) then
      if (ListCheckPresent(acc_model % constants, 'Permittivity of Vacuum')) &
          acc_eps0 = GetConstReal(acc_model % constants, 'Permittivity of Vacuum', found)
    end if
    if (.not. found) acc_eps0 = 8.854187817d-12
    acc_omega = GetAngularFrequency(found=found)
    if (.not. found) call Fatal(caller, 'Angular frequency required')
  end subroutine ema_init_stuff

  subroutine ema_fetch_mesh() bind(c, name="ema_fetch_mesh")
    acc_mesh => GetMesh()
  end subroutine ema_fetch_mesh

  function ema_active() result(n) bind(c, name="ema_active")
    integer(c_int) :: n
    n = GetNOFActive(acc_solver)
  end function ema_active

  ! Port-loop prologue: Info + (MaxPort > 1) matrix surgery + DefaultInitialize
  ! + maxmu/maxeps reset. Returns 1 when the port must be CYCLEd.
  function ema_port_begin(portind) result(r) bind(c, name="ema_port_begin")
    integer(c_int), value :: portind
    integer(c_int) :: r
    character(*), parameter :: caller = 'EMPortSolver'
    logical :: meactive
    integer :: n
    r = 0
    call Info(caller, 'Solving for port: ' // i2s(portind), level=10)
    if (acc_maxport > 1) then
      acc_emvar % perm = 0
      acc_solver % matrix => CreateMatrix(acc_model, acc_solver, acc_solver % mesh, &
          acc_emvar % perm, acc_emvar % dofs, matrix_crs, .false., &
          "Port Label " // i2s(portind), globalbubbles=acc_solver % globalbubbles, bcmode=.true.)
      if (parenv % pes > 1) then
        meactive = .false.
        if (associated(acc_solver % matrix)) then
          meactive = (acc_solver % matrix % numberofrows > 0)
        end if
        if (meactive) then
          allocate (acc_solver % matrix % massvalues(size(acc_solver % matrix % values)))
          acc_solver % matrix % massvalues = 0
        end if
        call ema_parallel_subset(meactive)
        if (.not. meactive) then
          r = 1
          return
        end if
      end if
      n = maxval(acc_emvar % perm) * acc_emvar % dofs
      allocate (acc_solver % matrix % rhs(n))
      acc_solver % matrix % rhs = 0.0_c_double
    end if
    call DefaultInitialize()
    acc_maxmu = 0.0_c_double
    acc_maxeps = 0.0_c_double
  end function ema_port_begin

  function ema_default_mode() result(m) bind(c, name="ema_default_mode")
    integer(c_int) :: m
    logical :: found
    if (acc_maxport == 0) then
      acc_mode = ListGetInteger(acc_params, 'Eigenfunction Index', found)
      if (.not. found) acc_mode = 1
    end if
    m = acc_mode
  end function ema_default_mode

  ! Per-element head: GetActiveElement + optional port filter + family + sizes.
  ! skip = 1 means CYCLE (port mismatch or EFamily > 4).
  subroutine ema_elem_begin(t, n, nd, efamily, skip) bind(c, name="ema_elem_begin")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, efamily, skip
    logical :: found
    skip = 0
    acc_elem => GetActiveElement(t, acc_solver)
    if (acc_maxport > 0) then
      acc_bc => GetBC(acc_elem)
      if (ListGetInteger(acc_bc, 'Port Index', found) /= acc_port_save) then
        skip = 1
        return
      end if
      acc_mode = ListGetInteger(acc_bc, 'Eigenfunction Index', found)
      if (.not. found) acc_mode = 1
    end if
    efamily = GetElementFamily(acc_elem)
    if (efamily > 4) then
      skip = 1
      return
    end if
    n = GetElementNOFNodes(acc_elem)
    nd = GetElementNOFDOFs(acc_elem)
  end subroutine ema_elem_begin

  subroutine ema_set_portind(portind) bind(c, name="ema_set_portind")
    integer(c_int), value :: portind
    acc_port_save = portind
  end subroutine ema_set_portind

  subroutine ema_check_edge(n, efamily) bind(c, name="ema_check_edge")
    integer(c_int), value :: n, efamily
    character(*), parameter :: caller = 'EMPortSolver'
    if (acc_edge_deg == 1) then
      if (n /= efamily) call Fatal(caller, 'A background mesh must have linear elements!')
    else
      select case (efamily)
      case (3)
        if (n < 6) call Fatal(caller, 'A background mesh needs 6-node triangles')
      case (4)
        if (n < 9) call Fatal(caller, 'A background mesh needs 9-node quads')
      end select
    end if
  end subroutine ema_check_edge

  ! -- LocalMatrix (no potentials): exact transcription ----------------------
  subroutine ema_lm(n, nd) bind(c, name="ema_lm")
    integer(c_int), value :: n, nd
    character(*), parameter :: caller = 'EMPortSolver'
    integer :: m, allocstat, t
    integer :: i, j, p, q, vdofs
    logical :: stat, found, gotnu, goteps
    real(c_double) :: weight, detj
    complex(c_double) :: nu, eps

    acc_ip_lm = GaussPoints(acc_elem, edgebasis=.true., &
        PReferenceElement=acc_piola, edgebasisdegree=acc_edge_deg)

    if (.not. allocated(acc_basis_lm)) then
      m = acc_mesh % maxelementdofs
      allocate (acc_wbasis_lm(m, 3), acc_curlw_lm(m, 3), acc_basis_lm(m), &
          acc_dbasis_lm(m, 3), acc_stiff_lm(m, m), acc_mass_lm(m, m), &
          acc_force_lm(m), stat=allocstat)
      if (allocstat /= 0) call Fatal(caller, 'Local storage allocation failed')
    end if

    call GetElementNodes(acc_nodes_lm, acc_elem)

    acc_stiff_lm = cmplx(0.0_c_double, 0.0_c_double, kind=c_double)
    acc_mass_lm = cmplx(0.0_c_double, 0.0_c_double, kind=c_double)
    acc_force_lm = cmplx(0.0_c_double, 0.0_c_double, kind=c_double)
    gotnu = .false.
    goteps = .false.

    vdofs = nd - n

    do t = 1, acc_ip_lm % n
      stat = ElementInfo(acc_elem, acc_nodes_lm, acc_ip_lm % u(t), &
          acc_ip_lm % v(t), acc_ip_lm % w(t), detj, acc_basis_lm, acc_dbasis_lm, &
          edgebasis=acc_wbasis_lm, rotbasis=acc_curlw_lm, usolver=acc_solver)
      weight = acc_ip_lm % s(t) * detj

      if (t == 1 .or. gotnu) then
        nu = ListGetElementComplex(acc_nucoeff_h, acc_basis_lm, acc_elem, gotnu, &
            gausspoint=t)
      end if
      if (.not. gotnu) nu = ListGetElementRealParent(acc_nucoeff_h, acc_basis_lm, &
          acc_elem, found)
      if (gotnu .or. found) then
        nu = acc_mu0inv * nu
      else
        nu = acc_mu0inv
      end if

      if (t == 1 .or. goteps) then
        eps = ListGetElementComplex(acc_epscoeff_h, acc_basis_lm, acc_elem, goteps, &
            gausspoint=t)
      end if
      if (.not. goteps) eps = ListGetElementRealParent(acc_epscoeff_h, acc_basis_lm, &
          acc_elem, found)
      if (goteps .or. found) then
        eps = acc_eps0 * eps
      else
        eps = acc_eps0
      end if

      acc_maxmu = max(acc_maxmu, real(1 / nu))
      acc_maxeps = max(acc_maxeps, real(eps))

      do p = 1, n
        do q = 1, n
          acc_stiff_lm(p, q) = acc_stiff_lm(p, q) - weight * eps * &
              acc_basis_lm(p) * acc_basis_lm(q)
        end do

        i = p
        do q = 1, vdofs
          j = n + q
          acc_stiff_lm(i, j) = acc_stiff_lm(i, j) + eps * &
              sum(acc_wbasis_lm(q, :) * acc_dbasis_lm(p, :)) * weight
          acc_stiff_lm(j, i) = acc_stiff_lm(j, i) + nu * &
              sum(acc_wbasis_lm(q, :) * acc_dbasis_lm(p, :)) * weight
        end do
      end do

      do p = 1, vdofs
        i = n + p
        do q = 1, vdofs
          j = n + q
          acc_stiff_lm(i, j) = acc_stiff_lm(i, j) + weight * (nu * &
              sum(acc_curlw_lm(q, :) * acc_curlw_lm(p, :)) - &
              acc_omega**2 * eps * sum(acc_wbasis_lm(q, :) * acc_wbasis_lm(p, :)))
          acc_mass_lm(i, j) = acc_mass_lm(i, j) + weight * nu * &
              sum(acc_wbasis_lm(q, :) * acc_wbasis_lm(p, :))
        end do
      end do
    end do

    call DefaultUpdateEquations(acc_stiff_lm, acc_force_lm)
    call DefaultUpdateMass(acc_mass_lm)
  end subroutine ema_lm

  ! -- LocalMatrix_With_Potentials: exact transcription ----------------------
  subroutine ema_lm_pot(n, nd) bind(c, name="ema_lm_pot")
    integer(c_int), value :: n, nd
    character(*), parameter :: caller = 'EMPortSolver'
    integer :: m, allocstat, t
    integer :: i, j, p, q, vdofs
    logical :: stat, found, gotnu, goteps
    real(c_double) :: weight, detj
    complex(c_double) :: nu, eps
    integer :: ndofs, np

    acc_ip_lp = GaussPoints(acc_elem, edgebasis=.true., &
        PReferenceElement=acc_piola, edgebasisdegree=acc_edge_deg)

    if (.not. allocated(acc_basis_lp)) then
      m = acc_mesh % maxelementdofs
      allocate (acc_wbasis_lp(m, 3), acc_curlw_lp(m, 3), acc_basis_lp(m), &
          acc_dbasis_lp(m, 3), acc_stiff_lp(m, m), acc_mass_lp(m, m), &
          acc_force_lp(m), stat=allocstat)
      if (allocstat /= 0) call Fatal(caller, 'Local storage allocation failed')
    end if

    call GetElementNodes(acc_nodes_lp, acc_elem)

    acc_stiff_lp = cmplx(0.0_c_double, 0.0_c_double, kind=c_double)
    acc_mass_lp = cmplx(0.0_c_double, 0.0_c_double, kind=c_double)
    acc_force_lp = cmplx(0.0_c_double, 0.0_c_double, kind=c_double)
    gotnu = .false.
    goteps = .false.

    ndofs = maxval(acc_solver % def_dofs(GetElementFamily(acc_elem), :, 1))
    np = n * ndofs

    vdofs = nd - np

    do t = 1, acc_ip_lp % n
      stat = ElementInfo(acc_elem, acc_nodes_lp, acc_ip_lp % u(t), &
          acc_ip_lp % v(t), acc_ip_lp % w(t), detj, acc_basis_lp, acc_dbasis_lp, &
          edgebasis=acc_wbasis_lp, rotbasis=acc_curlw_lp, usolver=acc_solver)
      weight = acc_ip_lp % s(t) * detj

      if (t == 1 .or. gotnu) then
        nu = ListGetElementComplex(acc_nucoeff_h, acc_basis_lp, acc_elem, gotnu, &
            gausspoint=t)
      end if
      if (.not. gotnu) nu = ListGetElementRealParent(acc_nucoeff_h, acc_basis_lp, &
          acc_elem, found)
      if (gotnu .or. found) then
        nu = acc_mu0inv * nu
      else
        nu = acc_mu0inv
      end if

      if (t == 1 .or. goteps) then
        eps = ListGetElementComplex(acc_epscoeff_h, acc_basis_lp, acc_elem, goteps, &
            gausspoint=t)
      end if
      if (.not. goteps) eps = ListGetElementRealParent(acc_epscoeff_h, acc_basis_lp, &
          acc_elem, found)
      if (goteps .or. found) then
        eps = acc_eps0 * eps
      else
        eps = acc_eps0
      end if

      acc_maxmu = max(acc_maxmu, real(1 / nu))
      acc_maxeps = max(acc_maxeps, real(eps))

      do p = 1, n
        i = (p - 1) * ndofs + 1
        do q = 1, n
          j = (q - 1) * ndofs + 1
          acc_stiff_lp(i, j) = acc_stiff_lp(i, j) + (nu * &
              sum(acc_dbasis_lp(q, :) * acc_dbasis_lp(p, :)) - &
              acc_omega**2 * eps * acc_basis_lp(q) * acc_basis_lp(p)) * weight
        end do

        do q = 1, vdofs
          j = np + q
          acc_stiff_lp(i, j) = acc_stiff_lp(i, j) - nu * &
              sum(acc_wbasis_lp(q, :) * acc_dbasis_lp(p, :)) * weight
        end do

        do q = 1, n
          j = (q - 1) * ndofs + 2
          acc_stiff_lp(i, j) = acc_stiff_lp(i, j) + nu * &
              sum(acc_dbasis_lp(q, :) * acc_dbasis_lp(p, :)) * weight
        end do

        i = (p - 1) * ndofs + 2
        do q = 1, n
          j = (q - 1) * ndofs + 2
          acc_stiff_lp(i, j) = acc_stiff_lp(i, j) - weight * acc_omega**2 * eps / nu * &
              acc_basis_lp(p) * acc_basis_lp(q)
        end do

        do q = 1, vdofs
          j = np + q
          acc_stiff_lp(i, j) = acc_stiff_lp(i, j) + &
              sum(acc_wbasis_lp(q, :) * acc_dbasis_lp(p, :)) * weight
        end do
      end do

      do p = 1, vdofs
        i = np + p
        do q = 1, vdofs
          j = np + q
          acc_stiff_lp(i, j) = acc_stiff_lp(i, j) + weight * (nu * &
              sum(acc_curlw_lp(q, :) * acc_curlw_lp(p, :)) - &
              acc_omega**2 * eps * sum(acc_wbasis_lp(q, :) * acc_wbasis_lp(p, :)))
          acc_mass_lp(i, j) = acc_mass_lp(i, j) + weight * nu * &
              sum(acc_wbasis_lp(q, :) * acc_wbasis_lp(p, :))
        end do

        do q = 1, n
          j = (q - 1) * ndofs + 1
          acc_mass_lp(i, j) = acc_mass_lp(i, j) - weight * nu * &
              sum(acc_dbasis_lp(q, :) * acc_wbasis_lp(p, :))

          j = (q - 1) * ndofs + 2
          acc_stiff_lp(i, j) = acc_stiff_lp(i, j) + weight * acc_omega**2 * eps * &
              sum(acc_dbasis_lp(q, :) * acc_wbasis_lp(p, :))
          acc_mass_lp(i, j) = acc_mass_lp(i, j) - weight * nu * &
              sum(acc_dbasis_lp(q, :) * acc_wbasis_lp(p, :))
        end do
      end do
    end do

    call DefaultUpdateEquations(acc_stiff_lp, acc_force_lp)
    call DefaultUpdateMass(acc_mass_lp)
  end subroutine ema_lm_pot

  subroutine ema_finish_bulk() bind(c, name="ema_finish_bulk")
    call DefaultFinishBulkAssembly()
  end subroutine ema_finish_bulk

  subroutine ema_finish_assembly() bind(c, name="ema_finish_assembly")
    call DefaultFinishAssembly()
  end subroutine ema_finish_assembly

  subroutine ema_dirichlet() bind(c, name="ema_dirichlet")
    call DefaultDirichletBCs()
  end subroutine ema_dirichlet

  subroutine ema_shift() bind(c, name="ema_shift")
    character(*), parameter :: caller = 'EMPortSolver'
    logical :: found
    real(c_double) :: betalim
    if (ListGetLogical(acc_params, 'Eigen System Shift Automatic', found)) then
      acc_maxeps = ParallelReduction(acc_maxeps, 2)
      acc_maxmu = ParallelReduction(acc_maxmu, 2)
      betalim = acc_omega * sqrt(acc_maxeps * acc_maxmu)
      call ListAddConstReal(acc_params, 'Eigen System Shift', -betalim**2)
      write (message, '(A,ES15.6)') 'Eigen System Shift set to ', -betalim**2
      call Info(caller, message, level=7)
    end if
  end subroutine ema_shift

  function ema_solve() result(n) bind(c, name="ema_solve")
    real(c_double) :: n
    n = DefaultSolve()
  end function ema_solve

  subroutine ema_beta() bind(c, name="ema_beta")
    character(*), parameter :: caller = 'EMPortSolver'
    acc_beta = sqrt(-acc_solver % variable % eigenvalues(acc_mode))
    write (message, '(A,2ES15.6)') 'Propagation constant beta: ', real(acc_beta), &
        aimag(acc_beta)
    call Info(caller, message, level=5)
    call ListAddConstReal(acc_model % simulation, 'res: Port Beta ' // i2s(acc_port_save), &
        real(acc_beta))
    acc_betasum = acc_betasum + real(acc_beta)
  end subroutine ema_beta

  subroutine ema_beta_to_bcs() bind(c, name="ema_beta_to_bcs")
    type(valuelist_t), pointer :: bc
    logical :: found
    integer :: i, j
    character(:), allocatable :: ptype
    do i = 1, acc_model % numberofbcs
      bc => acc_model % bcs(i) % values
      ptype = ListGetString(bc, 'Port Type', found)
      if (ptype == 'eigenmode') then
        j = ListGetInteger(bc, "Port Index", found)
        if (j == acc_port_save .or. acc_maxport == 0) then
          call ListAddConstReal(bc, 'Port Beta', real(acc_beta))
          call ListAddConstReal(bc, 'Port Beta Im', aimag(acc_beta))
        end if
      end if
      j = ListGetInteger(bc, "Port Beta Parent", found)
      if (found .and. j == acc_port_save) then
        call ListAddConstReal(bc, 'Port Beta', real(acc_beta))
        call ListAddConstReal(bc, 'Port Beta Im', aimag(acc_beta))
      end if
    end do
  end subroutine ema_beta_to_bcs

  ! Power-normalization prologue. Returns 1 when the block runs.
  function ema_norm_begin() result(r) bind(c, name="ema_norm_begin")
    integer(c_int) :: r
    logical :: found
    r = 0
    if (ListGetLogical(acc_params, 'Eigenvector Normalization by Power', found)) then
      acc_power = cmplx(0.0_c_double, 0.0_c_double, kind=c_double)
      r = 1
    end if
  end function ema_norm_begin

  subroutine ema_norm_elem(t) bind(c, name="ema_norm_elem")
    integer(c_int), value :: t
    integer :: efamily
    logical :: found
    acc_elem => GetActiveElement(t, acc_solver)
    efamily = GetElementFamily(acc_elem)
    if (efamily > 4) return
    if (CoordinateSystemDimension() == 2) then
      continue
    else
      if (acc_maxport > 0) then
        acc_bc => GetBC(acc_elem)
        if (ListGetInteger(acc_bc, 'Port Index', found) /= acc_port_save) return
      end if
    end if
    call acc_power_core(.false.)
  end subroutine ema_norm_elem

  subroutine ema_norm_finish() bind(c, name="ema_norm_finish")
    character(*), parameter :: caller = 'EMPortSolver'
    acc_power = ParallelReduction(acc_power)
    write (message, '(A,2ES15.6)') 'Scaling by (average) port power: ', &
        0.5_c_double * real(acc_power)
    call Info(caller, message, level=5)
    acc_solver % variable % eigenvectors(acc_mode, :) = &
        acc_solver % variable % eigenvectors(acc_mode, :) / &
        cmplx(sqrt(real(acc_power) / 2.0_c_double), 0.0_c_double, kind=c_double)
  end subroutine ema_norm_finish

  ! Shared CalculatePortPower core (single internal procedure in the
  ! original, with one set of SAVE arrays used by both call sites).
  subroutine acc_power_core(use_e2)
    logical, intent(in) :: use_e2
    integer :: m, allocstat, t
    integer :: i, j, vdofs, np, ndofs, n
    logical :: stat, found, gotnu
    real(c_double) :: weight, detj
    complex(c_double) :: nu, ef(3), gradez(3), v
    n = GetElementNOFNodes(acc_elem)
    acc_ip_pw = GaussPoints(acc_elem, edgebasis=.true., &
        PReferenceElement=acc_piola, edgebasisdegree=acc_edge_deg)
    if (.not. allocated(acc_basis_pw)) then
      m = acc_mesh % maxelementdofs
      allocate (acc_wbasis_pw(m, 3), acc_curlw_pw(m, 3), acc_basis_pw(m), &
          acc_dbasis_pw(m, 3), acc_re_pw(m), acc_im_pw(m), stat=allocstat)
      if (allocstat /= 0) call Fatal('EMPortSolver', 'Local storage allocation failed')
    end if
    call GetElementNodes(acc_nodes_pw, acc_elem)
    ndofs = maxval(acc_solver % def_dofs(GetElementFamily(acc_elem), :, 1))
    np = n * ndofs
    vdofs = GetElementNOFDOFs(acc_elem) - np
    call GetScalarLocalEigenmode(acc_re_pw, uelement=acc_elem, usolver=acc_solver, &
        noeigen=acc_mode, complexpart=.false.)
    call GetScalarLocalEigenmode(acc_im_pw, uelement=acc_elem, usolver=acc_solver, &
        noeigen=acc_mode, complexpart=.true.)
    gotnu = .false.
    do t = 1, acc_ip_pw % n
      stat = ElementInfo(acc_elem, acc_nodes_pw, acc_ip_pw % u(t), acc_ip_pw % v(t), &
          acc_ip_pw % w(t), detj, acc_basis_pw, acc_dbasis_pw, &
          edgebasis=acc_wbasis_pw, rotbasis=acc_curlw_pw, usolver=acc_solver)
      weight = acc_ip_pw % s(t) * detj
      if (t == 1 .or. gotnu) then
        nu = ListGetElementComplex(acc_nucoeff_h, acc_basis_pw, acc_elem, gotnu, &
            gausspoint=t)
      end if
      if (.not. gotnu) nu = ListGetElementRealParent(acc_nucoeff_h, acc_basis_pw, &
          acc_elem, found)
      if (gotnu .or. found) then
        nu = acc_mu0inv * nu
      else
        nu = acc_mu0inv
      end if
      ef = cmplx(0.0_c_double, 0.0_c_double, kind=c_double)
      do i = 1, vdofs
        ef(:) = ef(:) + cmplx(acc_re_pw(np + i), acc_im_pw(np + i), kind=c_double) * &
            acc_wbasis_pw(i, :)
      end do
      if (acc_usev) then
        do j = 1, 3
          ef(j) = ef(j) - cmplx(sum(acc_re_pw(2:np:ndofs) * acc_dbasis_pw(1:n, j)), &
              sum(acc_im_pw(2:np:ndofs) * acc_dbasis_pw(1:n, j)), kind=c_double)
        end do
        v = cmplx(sum(acc_re_pw(2:np:ndofs) * acc_basis_pw(1:n)), &
            sum(acc_im_pw(2:np:ndofs) * acc_basis_pw(1:n)), kind=c_double)
      end if
      gradez = cmplx(0.0_c_double, 0.0_c_double, kind=c_double)
      if (acc_usev) then
        do j = 1, 3
          gradez(j) = gradez(j) + cmplx(sum(acc_re_pw(1:np:ndofs) * acc_dbasis_pw(1:n, j)), &
              sum(acc_im_pw(1:np:ndofs) * acc_dbasis_pw(1:n, j)), kind=c_double)
        end do
        gradez(:) = gradez(:) * acc_im * acc_beta
      else
        do i = 1, n
          gradez(:) = gradez(:) + cmplx(acc_re_pw(i), acc_im_pw(i), &
              kind=c_double) * acc_dbasis_pw(i, :)
        end do
        gradez(:) = gradez(:) / (acc_im * acc_beta)
      end if
      acc_power = acc_power + nu / acc_omega * (acc_beta * sum(ef * conjg(ef)) - &
          acc_im * sum(ef * conjg(gradez))) * weight
      if (use_e2) acc_e2 = acc_e2 + sum(ef * conjg(ef)) * weight
    end do
  end subroutine acc_power_core

  ! Nonhomogeneous-modes prologue. Returns the mode count (0 = skip).
  function ema_nonhomog_begin() result(r) bind(c, name="ema_nonhomog_begin")
    integer(c_int) :: r
    character(*), parameter :: caller = 'EMPortSolver'
    logical :: found
    integer :: i, m
    i = ListGetInteger(acc_params, 'Number of Nonhomogeneous Modes', found, &
        maxv=size(acc_solver % variable % eigenvalues))
    if (i > 0) then
      call Info(caller, 'Solving an additional component to satisfy nonhomogeneous BCs', &
          level=5)
      call ListAddLogical(acc_params, 'Eigen Analysis', .false.)
      if (parenv % pes > 1) then
        m = acc_solver % matrix % numberofrows
        allocate (acc_upart(m), acc_vpart(m))
      end if
    end if
    r = i
  end function ema_nonhomog_begin

  subroutine ema_nonhomog(j) bind(c, name="ema_nonhomog")
    integer(c_int), value :: j
    complex(c_double), pointer :: cvalues(:), cu(:), cv(:)
    complex(c_double) :: udotv, udotu
    real(c_double) :: norm
    integer :: k, m
    acc_solver % matrix % values = acc_solver % matrix % values - &
        acc_solver % variable % eigenvalues(j) * acc_solver % matrix % massvalues
    acc_solver % variable % values = 0.0_c_double
    norm = DefaultSolve()
    cvalues => ComplexVariableValues(acc_solver % variable)
    if (parenv % pes > 1) then
      m = acc_solver % matrix % numberofrows / 2
      do k = 1, m
        acc_vpart(2 * k - 1) = real(acc_solver % variable % eigenvectors(j, k))
        acc_vpart(2 * k) = aimag(acc_solver % variable % eigenvectors(j, k))
      end do
      call PartitionVector(acc_solver % matrix, acc_upart, acc_vpart)
      cu => ComplexValues(acc_upart, m)
      call PartitionVector(acc_solver % matrix, acc_vpart, acc_solver % variable % values)
      cv => ComplexValues(acc_vpart, m)
      udotv = ParallelCdot(m, cu, cv)
      udotu = ParallelCdot(m, cu, cu)
    else
      m = acc_solver % matrix % numberofrows / 2
      udotv = sum(conjg(acc_solver % variable % eigenvectors(j, 1:m)) * cvalues(1:m))
      udotu = sum(conjg(acc_solver % variable % eigenvectors(j, 1:m)) * &
          acc_solver % variable % eigenvectors(j, 1:m))
    end if
    acc_solver % variable % eigenvectors(j, 1:m) = cvalues(1:m) - &
        udotv / udotu * acc_solver % variable % eigenvectors(j, 1:m)
    acc_solver % matrix % values = acc_solver % matrix % values + &
        acc_solver % variable % eigenvalues(j) * acc_solver % matrix % massvalues
  end subroutine ema_nonhomog

  subroutine ema_nonhomog_end() bind(c, name="ema_nonhomog_end")
    if (parenv % pes > 1) then
      deallocate (acc_upart, acc_vpart)
    end if
    call ListAddLogical(acc_params, 'Eigen Analysis', .true.)
  end subroutine ema_nonhomog_end

  ! Impedance prologue. Returns 1 when the block runs.
  function ema_impedance_check() result(r) bind(c, name="ema_impedance_check")
    integer(c_int) :: r
    logical :: found
    r = 0
    if (ListCheckPresentAnyBC(acc_model, 'Calculate Impedance') .or. &
        (CoordinateSystemDimension() == 2 .and. &
        ListGetLogical(acc_params, 'Calculate Impedance', found))) r = 1
  end function ema_impedance_check

  subroutine ema_impedance_begin() bind(c, name="ema_impedance_begin")
    acc_e2 = cmplx(0.0_c_double, 0.0_c_double, kind=c_double)
    acc_power = cmplx(0.0_c_double, 0.0_c_double, kind=c_double)
    acc_output_z = .false.
  end subroutine ema_impedance_begin

  subroutine ema_impedance_elem(t) bind(c, name="ema_impedance_elem")
    integer(c_int), value :: t
    integer :: efamily
    logical :: found
    acc_elem => GetActiveElement(t, acc_solver)
    efamily = GetElementFamily(acc_elem)
    if (efamily > 4) return
    if (CoordinateSystemDimension() == 2) then
      continue
    else
      if (acc_maxport > 0) then
        acc_bc => GetBC(acc_elem)
        if (ListGetInteger(acc_bc, 'Port Index', found) /= acc_port_save) return
      end if
      if (.not. ListGetLogical(acc_bc, 'Calculate Impedance', found)) return
    end if
    call acc_power_core(.true.)
    if (.not. acc_output_z) acc_output_z = .true.
  end subroutine ema_impedance_elem

  subroutine ema_impedance_finish() bind(c, name="ema_impedance_finish")
    character(*), parameter :: caller = 'EMPortSolver'
    acc_power = ParallelReduction(acc_power)
    acc_e2 = ParallelReduction(acc_e2)
    if (.not. acc_output_z) return
    write (message, '(A,2ES15.6)') 'Port power (average): ', 0.5_c_double * real(acc_power)
    call Info(caller, message, level=5)
    write (message, '(A,2ES15.6)') 'Port power Im (average): ', 0.5_c_double * aimag(acc_power)
    call Info(caller, message, level=5)
    write (message, '(A,2ES15.6)') 'Modal (wave) impedance: ', real(acc_e2) / real(acc_power)
    call Info(caller, message, level=5)
    write (message, '(A,2ES15.6)') 'Port impedance (Re): ', real(1.0_c_double / conjg(acc_power))
    call Info(caller, message, level=5)
    write (message, '(A,2ES15.6)') 'Port impedance (Im): ', aimag(1.0_c_double / conjg(acc_power))
    call Info(caller, message, level=5)
    call ListAddConstReal(acc_model % simulation, 'res: Port Power ' // i2s(acc_port_save), &
        0.5_c_double * real(acc_power))
    call ListAddConstReal(acc_model % simulation, 'res: Port Power Im ' // i2s(acc_port_save), &
        0.5_c_double * aimag(acc_power))
    call ListAddConstReal(acc_model % simulation, 'res: Port Impedance ' // i2s(acc_port_save), &
        real(1.0_c_double / conjg(acc_power)))
    call ListAddConstReal(acc_model % simulation, 'res: Port Impedance Im ' // i2s(acc_port_save), &
        aimag(1.0_c_double / conjg(acc_power)))
  end subroutine ema_impedance_finish

  function ema_calc_nodal() result(r) bind(c, name="ema_calc_nodal")
    integer(c_int) :: r
    logical :: found
    r = merge(1, 0, ListGetLogical(acc_params, 'Calculate Nodal Field', found))
  end function ema_calc_nodal

  subroutine ema_port_collect() bind(c, name="ema_port_collect")
    integer :: i, n
    if (acc_maxport > 1) then
      call FreeMatrix(acc_solver % matrix)
      acc_solver % matrix => null()
      n = size(acc_emvar % perm)
      do i = 1, n
        if (acc_saveperm(i) > 0 .and. acc_emvar % perm(i) > 0) then
          acc_saveev(:, acc_saveperm(i)) = acc_emvar % eigenvectors(:, acc_emvar % perm(i))
        end if
      end do
    end if
  end subroutine ema_port_collect

  subroutine ema_parallel_reactivate() bind(c, name="ema_parallel_reactivate")
    if (parenv % pes > 1) call ParallelActive(.true.)
  end subroutine ema_parallel_reactivate

  subroutine ema_default_finish() bind(c, name="ema_default_finish")
    call DefaultFinish()
  end subroutine ema_default_finish

  subroutine ema_multiport_restore() bind(c, name="ema_multiport_restore")
    if (acc_maxport > 1) then
      acc_emvar % perm = acc_saveperm
      deallocate (acc_emvar % eigenvectors)
      acc_emvar % eigenvectors => acc_saveev
      nullify (acc_saveev)
      deallocate (acc_saveperm)
    end if
  end subroutine ema_multiport_restore

  subroutine ema_set_norm() bind(c, name="ema_set_norm")
    acc_solver % variable % norm = acc_betasum
  end subroutine ema_set_norm

  subroutine ema_done_info() bind(c, name="ema_done_info")
    character(*), parameter :: caller = 'EMPortSolver'
    call Info(caller, 'All done', level=12)
  end subroutine ema_done_info

  ! ParallelActiveSubset: exact transcription of the CONTAINS internal.
  ! In the serial proof (PEs == 1) this returns immediately, mirroring
  ! the original's guard.
  subroutine ema_parallel_subset(meactive) bind(c, name="ema_parallel_subset")
    logical, intent(inout) :: meactive
    integer :: n
    integer, allocatable :: memb(:)
    type(matrix_t), pointer :: m
    integer :: comm_active, group_active, group_world, ierr
    integer :: i
    if (parenv % pes == 1) return
    call ParallelActive(meactive)
    n = count(parenv % active)
    m => acc_solver % matrix
    if (n > 0 .and. n < parenv % pes) then
      if (associated(acc_solver % matrix)) then
        if (acc_solver % matrix % comm /= elmer_comm_world .and. &
            acc_solver % matrix % comm /= mpi_comm_null) &
            call mpi_comm_free(acc_solver % matrix % comm, ierr)
      end if
      call mpi_comm_group(elmer_comm_world, group_world, ierr)
      allocate (memb(n))
      n = 0
      do i = 1, parenv % pes
        if (parenv % active(i)) then
          n = n + 1
          memb(n) = i - 1
        end if
      end do
      call mpi_group_incl(group_world, n, memb, group_active, ierr)
      deallocate (memb)
      call mpi_comm_create(elmer_comm_world, group_active, comm_active, ierr)
      m => acc_solver % matrix
      do while (associated(m))
        m % comm = comm_active
        m => m % parent
      end do
      parenv % activecomm = comm_active
      if (any(parenv % active(minoutputpe + 1:min(maxoutputpe + 1, parenv % pes)))) then
        if (parenv % mype >= minoutputpe .and. parenv % mype <= maxoutputpe) then
          outputpe = parenv % mype
        else
          outputpe = -1
        end if
      else
        do i = 1, parenv % pes
          if (parenv % active(i)) exit
        end do
        outputpe = -1
        if (i - 1 == parenv % mype) then
          outputpe = i - 1
        else if (i > parenv % pes .and. parenv % mype == 0) then
          outputpe = 0
        end if
      end if
    else
      m => acc_solver % matrix
      do while (associated(m))
        m % comm = elmer_comm_world
        m => m % parent
      end do
      if (.not. associated(acc_solver % matrix)) parenv % active = .true.
      if (parenv % mype >= minoutputpe .and. parenv % mype <= maxoutputpe) then
        outputpe = parenv % mype
      else
        outputpe = -1
      end if
    end if
    if (associated(acc_solver % matrix)) then
      if (acc_solver % parallel .and. meactive) then
        if (associated(acc_solver % mesh % parallelInfo % GInterface)) then
          parenv % activecomm = acc_solver % matrix % comm
          if (.not. associated(acc_solver % matrix % parmatrix)) then
            call ParallelInitMatrix(acc_solver, acc_solver % matrix)
          end if
          parenv % activecomm = acc_solver % matrix % comm
        end if
      end if
    end if
  end subroutine ema_parallel_subset

  ! EMPortPost: exact transcription of the CONTAINS internal procedure.
  ! portind/maxport select the element subset and gate the final assembly.
  subroutine ema_post(portind, maxport) bind(c, name="ema_post")
    integer(c_int), value :: portind, maxport
    character(*), parameter :: caller = 'EMPortSolver'
    type(variable_t), pointer :: ef
    type(element_t), pointer :: element
    integer :: i, j, k, n, p, q, nd, normal_ind(1), dofs, vdofs, np, ndofs
    integer :: active
    real(c_double) :: u, v, w, detj, s, norm, totnorm
    real(c_double) :: reez, imez, ree(3), ime(3), normal(3), rel(3), iml(3)
    real(c_double) :: rev, imv
    complex(c_double) :: ez
    complex(c_double), parameter :: im = (0.0_c_double, 1.0_c_double)
    logical :: found, stat
    type(solver_t), pointer :: psolver
    character(:), allocatable :: eqname
    logical :: skip_realloc
    psolver => acc_solver
    if (.not. acc_postalloc) then
      allocate (acc_post)
      call ListCopyPrefixedKeywords(acc_solver % values, acc_post % values, 'post:')
      acc_post % mesh => acc_mesh
      i = size(acc_solver % def_dofs, 1)
      j = size(acc_solver % def_dofs, 2)
      k = size(acc_solver % def_dofs, 3)
      allocate (acc_post % def_dofs(i, j, k))
      acc_post % def_dofs = 0
      acc_post % def_dofs(:, :, 1) = 1
      n = acc_mesh % maxelementdofs
      if (acc_usev) then
        dofs = 8
      else
        dofs = 6
      end if
      acc_postdofs = dofs
      allocate (acc_postmass(n, n), acc_postlf(n, dofs), acc_postwb(n, 3), &
          acc_postcb(n, 3), acc_postbasis(n), acc_postdbasis(n, 3), &
          acc_postperm(n), acc_postre(n), acc_postim(n))
    end if
    if (acc_postalloc) then
      if (size(acc_postnperm) == size(acc_solver % variable % perm)) then
        skip_realloc = .true.
      else
        deallocate (acc_postnperm)
        call FreeMatrix(acc_post % matrix)
        skip_realloc = .false.
      end if
    else
      skip_realloc = .false.
    end if
    acc_postalloc = .true.
    if (.not. skip_realloc) then
      allocate (acc_postnperm(size(acc_solver % variable % perm)))
      acc_postnperm = 0
      eqname = ListGetString(acc_params, 'Equation')
      call ListAddString(acc_post % values, 'Equation', trim(eqname) // '_post')
      acc_post % matrix => CreateMatrix(acc_model, acc_post, acc_mesh, acc_postnperm, &
          1, matrix_crs, .false., eqname, nodaldofsonly=.true.)
      acc_post % matrix % values = 0.0_c_double
      call VariableAddVector(acc_mesh % variables, acc_mesh, acc_post, &
          'EM2D tmp', 1, perm=acc_postnperm, output=.false.)
      acc_post % variable => VariableGet(acc_mesh % variables, 'EM2D tmp')
      if (.not. associated(acc_post % variable)) call Fatal(caller, 'Post solver field not found!')
      if (acc_usev) then
        call VariableAddVector(acc_mesh % variables, acc_mesh, acc_post, &
            'EF2D[EF2D Re:3 EF2D Im:3 Re V:1 Im V:1]', 8, perm=acc_postnperm, &
            secondary=.true., output=.true.)
      else
        call VariableAddVector(acc_mesh % variables, acc_mesh, acc_post, &
            'EF2D[EF2D Re:3 EF2D Im:3]', 6, perm=acc_postnperm, secondary=.true., &
            output=.true.)
      end if
      ef => VariableGet(acc_mesh % variables, 'EF2D')
      if (.not. associated(ef)) call Fatal(caller, 'Could not find field: EF2D!')
      acc_postef => ef
      n = acc_post % matrix % numberofrows
      allocate (acc_post % matrix % rhs(n))
      acc_postfsave => acc_post % matrix % rhs
      acc_post % matrix % rhs = 0.0_c_double
      acc_post % matrix % comm = acc_solver % matrix % comm
    end if
    active = GetNOFActive(acc_solver)
    n = acc_post % matrix % numberofrows
    if (.not. allocated(acc_postgf)) then
      allocate (acc_postgf(n, acc_postdofs - 1))
      acc_postgf = 0.0_c_double
    end if
    do k = 1, active
      element => GetActiveElement(k, acc_solver)
      if (maxport > 0) then
        acc_bc => GetBC(element)
        if (ListGetInteger(acc_bc, 'Port Index', found) /= portind) cycle
      end if
      n = GetElementNOFNodes(element)
      nd = GetElementNOFDOFs(usolver=acc_solver)
      ndofs = maxval(acc_solver % def_dofs(GetElementFamily(element), :, 1))
      np = n * ndofs
      vdofs = nd - np
      call GetElementNodes(acc_nodes_post, element, acc_solver)
      normal = NormalVector(element, acc_nodes_post)
      normal_ind = maxloc(abs(normal))
      if (normal(normal_ind(1)) < 0.0_c_double) normal = -normal
      call GetScalarLocalEigenmode(acc_postre, uelement=element, &
          usolver=acc_solver, noeigen=acc_mode, complexpart=.false.)
      call GetScalarLocalEigenmode(acc_postim, uelement=element, &
          usolver=acc_solver, noeigen=acc_mode, complexpart=.true.)
      acc_postmass = 0.0_c_double
      acc_postlf = 0.0_c_double
      acc_ip_post = GaussPoints(element, edgebasis=.true., &
          PReferenceElement=acc_piola, edgebasisdegree=acc_edge_deg)
      do i = 1, acc_ip_post % n
        u = acc_ip_post % u(i)
        v = acc_ip_post % v(i)
        w = acc_ip_post % w(i)
        stat = ElementInfo(element, acc_nodes_post, u, v, w, detj, acc_postbasis, &
            acc_postdbasis, edgebasis=acc_postwb, rotbasis=acc_postcb, usolver=psolver)
        s = acc_ip_post % s(i) * detj
        reez = sum(acc_postre(1:np:ndofs) * acc_postbasis(1:n))
        imez = sum(acc_postim(1:np:ndofs) * acc_postbasis(1:n))
        if (acc_usev) then
          ez = cmplx(reez, imez, kind=c_double) * im * acc_beta
        else
          ez = cmplx(reez, imez, kind=c_double) / (im * acc_beta)
        end if
        reez = real(ez)
        imez = aimag(ez)
        ree(:) = 0.0_c_double
        ime(:) = 0.0_c_double
        do p = 1, vdofs
          ree(:) = ree(:) + acc_postre(np + p) * acc_postwb(p, :)
          ime(:) = ime(:) + acc_postim(np + p) * acc_postwb(p, :)
        end do
        if (acc_usev) then
          do j = 1, 3
            ree(j) = ree(j) - sum(acc_postre(2:np:ndofs) * acc_postdbasis(1:n, j))
            ime(j) = ime(j) - sum(acc_postim(2:np:ndofs) * acc_postdbasis(1:n, j))
          end do
          rev = sum(acc_postre(2:np:ndofs) * acc_postbasis(1:n))
          imv = sum(acc_postim(2:np:ndofs) * acc_postbasis(1:n))
        end if
        rel = ree + normal * reez
        iml = ime + normal * imez
        do p = 1, n
          do q = 1, n
            acc_postmass(p, q) = acc_postmass(p, q) + s * acc_postbasis(p) * acc_postbasis(q)
          end do
          acc_postlf(p, 1:3) = acc_postlf(p, 1:3) + s * rel * acc_postbasis(p)
          acc_postlf(p, 4:6) = acc_postlf(p, 4:6) + s * iml * acc_postbasis(p)
          if (acc_usev) then
            acc_postlf(p, 7) = acc_postlf(p, 7) + s * rev * acc_postbasis(p)
            acc_postlf(p, 8) = acc_postlf(p, 8) + s * imv * acc_postbasis(p)
          end if
        end do
      end do
      if (acc_usev) then
        acc_postperm(1:n) = acc_post % variable % perm(2 * element % nodeindexes - 1)
      else
        acc_postperm(1:n) = acc_post % variable % perm(element % nodeindexes)
      end if
      call UpdateGlobalEquations(acc_post % matrix, acc_postmass, acc_post % matrix % rhs, &
          acc_postlf(1:n, 1), n, 1, acc_postperm(1:n), uelement=element)
      do j = 2, acc_postdofs
        call UpdateGlobalForce(acc_postgf(:, j - 1), &
            acc_postlf(1:n, j), n, 1, acc_postperm(1:n), uelement=element)
      end do
    end do
    if (portind < maxport) return
    totnorm = 0.0_c_double
    do j = 1, acc_postdofs
      call Info(caller, 'Solving for component: ' // i2s(j), level=10)
      if (j == 1) then
        acc_postfsave => acc_post % matrix % rhs
      else
        acc_post % matrix % rhs => acc_postgf(:, j - 1)
      end if
      acc_post % variable % values = 0.0_c_double
      norm = DefaultSolve(acc_post)
      totnorm = totnorm + norm**2
      acc_postef % values(j::acc_postdofs) = acc_post % variable % values
      call ListAddLogical(acc_post % values, 'Linear System refactorize', .false.)
    end do
    acc_post % variable % norm = sqrt(totnorm)
    acc_post % matrix % rhs => acc_postfsave
    totnorm = sqrt(totnorm)
    acc_post % variable % norm = totnorm
    acc_post % matrix % rhs => acc_postfsave
    deallocate (acc_postgf)
  end subroutine ema_post

  ! -- EMPortPotential_init -------------------------------------------------
  subroutine emp_init_begin(h_model, h_solver) bind(c, name="emp_init_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_pmodel)
    call c_f_pointer(h_solver, acc_psolver)
  end subroutine emp_init_begin

  subroutine emp_init_body() bind(c, name="emp_init_body")
    character(*), parameter :: caller = 'EMPortPotential_init'
    type(valuelist_t), pointer :: params, bc
    logical :: found
    integer :: i
    character(:), allocatable :: varname
    params => GetSolverParams()
    call ListAddNewString(params, 'Variable', 'Port Potential')
    varname = ListGetString(params, 'Variable', found)
    do i = 1, acc_pmodel % numberofbcs
      bc => acc_pmodel % bcs(i) % values
      if (ListGetLogical(bc, "Port Ground", found)) then
        call Info(caller, 'Setting "Port Potential" to zero where "Port Ground" is True', &
            level=10)
        call ListAddConstReal(bc, trim(varname), 0.0_c_double)
      else if (ListGetLogical(bc, "Port Feed", found)) then
        call Info(caller, 'Setting "Port Potential" to one where "Port Feed" is True', &
            level=10)
        call ListAddConstReal(bc, trim(varname), 1.0_c_double)
      end if
    end do
  end subroutine emp_init_body

  ! -- EMPortPotential -------------------------------------------------------
  subroutine emp_solve_begin(h_model, h_solver) bind(c, name="emp_solve_begin")
    type(c_ptr), value :: h_model, h_solver
    character(*), parameter :: caller = 'EMPortPotential'
    call c_f_pointer(h_model, acc_pmodel)
    call c_f_pointer(h_solver, acc_psolver)
    call Info(caller, '------------------------------------------------', level=7)
    call Info(caller, 'Solving static electric field for port boundaries', level=4)
    acc_pmesh => GetMesh()
  end subroutine emp_solve_begin

  subroutine emp_default_start() bind(c, name="emp_default_start")
    call DefaultStart()
  end subroutine emp_default_start

  subroutine emp_default_initialize() bind(c, name="emp_default_initialize")
    call DefaultInitialize()
  end subroutine emp_default_initialize

  function emp_active() result(n) bind(c, name="emp_active")
    integer(c_int) :: n
    n = GetNOFActive(acc_psolver)
  end function emp_active

  ! One bulk element: GetActiveElement + LocalMatrix(Element, n, nd+nb, nb,
  ! InitHandles). inith is the LOGICAL INOUT flag (1 on first call).
  subroutine emp_elem(t, inith) bind(c, name="emp_elem")
    integer(c_int), value :: t
    integer(c_int) :: inith
    character(*), parameter :: caller = 'EMPortPotential'
    integer :: n, nb, nd, m, allocstat
    logical :: inithandles
    inithandles = (inith /= 0)
    acc_pelem => GetActiveElement(t)
    n = GetElementNOFNodes(acc_pelem)
    nd = GetElementNOFDOFs(acc_pelem)
    nb = GetElementNOFBDOFs(acc_pelem)
    call emp_local_matrix(acc_pelem, n, nd + nb, nb, inithandles)
    inith = 0
  end subroutine emp_elem

  ! LocalMatrix internal of EMPortPotential: exact transcription.
  subroutine emp_local_matrix(element, n, nd, nb, inithandles)
    type(element_t), pointer, intent(in) :: element
    integer, intent(in) :: n, nd, nb
    logical, intent(inout) :: inithandles
    character(*), parameter :: caller = 'EMPortPotential'
    type(element_t), pointer :: parent
    real(c_double) :: eps0, weight
    real(c_double) :: sourceatip, epsatip, detj
    logical :: stat, found
    integer :: i, j, t, m, allocstat
    if (inithandles) then
      call ListInitElementKeyword(acc_psource_h, 'Body Force', 'Charge Density')
      call ListInitElementKeyword(acc_peps_h, 'Material', 'Relative Permittivity')
      found = .false.
      if (associated(acc_pmodel % constants)) then
        acc_peps0 = ListGetCReal(acc_pmodel % constants, 'Permittivity Of Vacuum', found)
      end if
      if (.not. found) acc_peps0 = 8.854187817e-12
      inithandles = .false.
      if (.not. allocated(acc_pbasis)) then
        m = acc_pmesh % maxelementdofs
        allocate (acc_pbasis(m), acc_pdbasis(m, 3), acc_pparentbasis(m), &
            acc_pstiff(m, m), acc_pforce(m), stat=allocstat)
        if (allocstat /= 0) call Fatal(caller, 'Local storage allocation failed')
      end if
    end if
    eps0 = acc_peps0
    acc_ip_pt = GaussPoints(element)
    call GetElementNodes(acc_nodes_pt, uelement=element)
    acc_pstiff = 0.0_c_double
    acc_pforce = 0.0_c_double
    parent => null()
    if (associated(element % boundaryinfo)) then
      parent => element % boundaryinfo % left
      if (.not. associated(parent)) then
        parent => element % boundaryinfo % right
      end if
    end if
    acc_pparent => parent
    do t = 1, acc_ip_pt % n
      stat = ElementInfo(element, acc_nodes_pt, acc_ip_pt % u(t), acc_ip_pt % v(t), &
          acc_ip_pt % w(t), detj, acc_pbasis, acc_pdbasis)
      weight = acc_ip_pt % s(t) * detj
      epsatip = ListGetElementReal(acc_peps_h, acc_pbasis, element, found, gausspoint=t)
      if (.not. found) then
        if (associated(parent)) then
          call SetParentBasis(element, n, acc_pbasis, parent, &
              parent % type % numberofnodes, acc_pparentbasis)
          epsatip = ListGetElementReal(acc_peps_h, acc_pparentbasis, parent, found)
        end if
        if (.not. found) then
          call Warn('EMPortPotential', 'Relative Permittivity not found, using 1.')
          epsatip = 1.0_c_double
        end if
      end if
      acc_pstiff(1:nd, 1:nd) = acc_pstiff(1:nd, 1:nd) + weight * &
          eps0 * epsatip * matmul(acc_pdbasis(1:nd, :), transpose(acc_pdbasis(1:nd, :)))
      sourceatip = ListGetElementReal(acc_psource_h, acc_pbasis, element, found)
      if (found) then
        acc_pforce(1:nd) = acc_pforce(1:nd) + weight * sourceatip * acc_pbasis(1:nd)
      end if
    end do
    call CondensateP(nd - nb, nb, acc_pstiff, acc_pforce)
    call DefaultUpdateEquations(acc_pstiff, acc_pforce, uelement=element)
  end subroutine emp_local_matrix

  subroutine emp_finish_bulk() bind(c, name="emp_finish_bulk")
    call DefaultFinishBulkAssembly()
  end subroutine emp_finish_bulk

  subroutine emp_finish_boundary() bind(c, name="emp_finish_boundary")
    call DefaultFinishBoundaryAssembly()
  end subroutine emp_finish_boundary

  subroutine emp_finish_assembly() bind(c, name="emp_finish_assembly")
    call DefaultFinishAssembly()
  end subroutine emp_finish_assembly

  subroutine emp_dirichlet() bind(c, name="emp_dirichlet")
    call DefaultDirichletBCs()
  end subroutine emp_dirichlet

  function emp_solve() result(n) bind(c, name="emp_solve")
    real(c_double) :: n
    n = DefaultSolve()
  end function emp_solve

  subroutine emp_finish() bind(c, name="emp_finish")
    character(*), parameter :: caller = 'EMPortPotential'
    call DefaultFinish()
    call Info(caller, '------------------------------------------------', level=7)
  end subroutine emp_finish

end module emport_acc

