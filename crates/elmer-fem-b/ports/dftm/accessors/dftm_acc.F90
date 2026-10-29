! T022-A DFT port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/DFTSolver.F90, compiled against the built
! .mod files of the target configure (T020/T021 entry pattern, ADR-057).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths — stays inside Fortran, so
! the Rust side is pure C calls. The wrapper module mirrors each
! entry's SAVE state (first-mesh scalars, table pointers, the internal
! procedures' SAVE Nodes/IP objects, and GRPULAY3's whole SAVE block)
! with identical lifetimes.
!
! All floating-point work (the three LocalMatrix integrations, the
! DOT_PRODUCT reductions, EXP/MIN/MAX/SQRT, the ABS of complex
! eigenvalue differences, MINVAL/MAXVAL/MINLOC/MAXLOC reports) and all
! formatted I/O (WRITE(Message,*)/PRINT*/INFO/WARN/ERROR/FATAL) stay in
! Fortran (same compiler), so bit-identity holds by construction
! wherever the call sequence matches. The Rust side transcribes the
! three entries' control flow (first-time blocks, assembly loops,
! relaxation dispatch, GRPULAY3 gating) statement for statement.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module dftm_acc
  use iso_c_binding
  use types
  use defutils
  use coordinatesystems
  use integration
  use matrixassembly
  use eleminfo
  use messages
  use linearalgebra
  use exchangecorrelations
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(valuelist_t), pointer, save :: acc_params => null()

  ! ChargeDensitySolver SAVE mirrors (internal LocalMatrix SAVE Nodes/IP
  ! live here with identical per-element lifetimes).
  integer(c_int), save :: cds_nofeigs = 0
  integer(c_int), save :: cds_eignum = 0
  real(c_double), pointer, save :: cds_eigw(:,:) => null()
  type(nodes_t), target, save :: cds_nodes
  type(gaussintegrationpoints_t), save :: cds_ip

  ! WaveFunctionSolver SAVE mirrors.
  integer(c_int), save :: wfs_nofigs = 0
  integer(c_int), save :: wfs_xctype = -1
  real(c_double), save :: wfs_shift = 0.0_c_double
  complex(c_double), allocatable, save :: wfs_eigstore(:)
  type(nodes_t), target, save :: wfs_nodes
  type(gaussintegrationpoints_t), save :: wfs_ip

  ! PoissonSolver SAVE mirrors, including GRPULAY3's SAVE block.
  integer(c_int), save :: ps_nofnuclei = 0
  integer(c_int), save :: ps_cdsnum = 0
  integer(c_int), save :: ps_relaxtype = 0
  integer(c_int), save :: ps_globsize = 0
  integer(c_int), save :: ps_startafter = 100000000
  integer(c_int), save :: ps_nofstored = 0
  real(c_double), save :: ps_totalnuc = 0.0_c_double
  real(c_double), save :: ps_pa = -100000.0_c_double
  real(c_double), save :: ps_pb = -100000.0_c_double
  real(c_double), save :: ps_pc = -100000.0_c_double
  real(c_double), save :: ps_nofiter = 0.0_c_double
  real(c_double), save :: ps_startrf = 1.0d8
  real(c_double), save :: ps_rf = 0.0_c_double
  logical, save :: ps_calc = .FALSE.
  logical, save :: ps_storeonly = .TRUE.
  logical, save :: ps_gfirst = .TRUE.
  logical, save :: ps_donenonrelax = .FALSE.
  real(c_double), pointer, save :: ps_nuctab(:,:) => null()
  integer(c_int), pointer, save :: ps_nucnodes(:) => null()
  integer(c_int), pointer, save :: ps_perm(:) => null()
  real(c_double), allocatable, save :: ps_old(:)
  real(c_double), allocatable, save :: ps_pulay(:,:)
  real(c_double), allocatable, save :: ps_resstor(:,:)
  real(c_double), allocatable, save :: ps_solstor(:,:)
  real(c_double), allocatable, save :: ps_rhs(:)
  real(c_double), allocatable, save :: ps_coeffs(:)
  type(nodes_t), target, save :: ps_nodes
  type(gaussintegrationpoints_t), save :: ps_ip

contains

  ! -- shared entry state -------------------------------------------
  subroutine dfm_begin(h_model, h_solver) bind(c, name="dfm_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine dfm_begin

  function dfm_max_elem_dofs() result(n) bind(c, name="dfm_max_elem_dofs")
    integer(c_int) :: n
    n = acc_solver % mesh % maxelementdofs
  end function dfm_max_elem_dofs

  function dfm_num_nodes() result(n) bind(c, name="dfm_num_nodes")
    integer(c_int) :: n
    n = acc_solver % mesh % numberofnodes
  end function dfm_num_nodes

  function dfm_nof_active() result(n) bind(c, name="dfm_nof_active")
    integer(c_int) :: n
    n = GetNOFActive()
  end function dfm_nof_active

  subroutine dfm_default_initialize() bind(c, name="dfm_default_initialize")
    call DefaultInitialize()
  end subroutine dfm_default_initialize

  subroutine dfm_default_update(ld, g, f) bind(c, name="dfm_default_update")
    integer(c_int), value :: ld
    real(c_double) :: g(ld, ld), f(ld)
    call DefaultUpdateEquations(g(1:ld, 1:ld), f(1:ld))
  end subroutine dfm_default_update

  subroutine dfm_default_update_mass(ld, g) bind(c, name="dfm_default_update_mass")
    integer(c_int), value :: ld
    real(c_double) :: g(ld, ld)
    call DefaultUpdateMass(g(1:ld, 1:ld))
  end subroutine dfm_default_update_mass

  subroutine dfm_default_finish() bind(c, name="dfm_default_finish")
    call DefaultFinishAssembly()
  end subroutine dfm_default_finish

  subroutine dfm_default_dirichlet() bind(c, name="dfm_default_dirichlet")
    call DefaultDirichletBCs()
  end subroutine dfm_default_dirichlet

  function dfm_default_solve() result(norm) bind(c, name="dfm_default_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function dfm_default_solve

  ! -- ChargeDensitySolver ------------------------------------------
  subroutine dfm_cds_info_start() bind(c, name="dfm_cds_info_start")
    call Info('ChargeDensityCalculation','...........................')
    call Info('ChargeDensityCalculation',' Charge density calculation')
    call Info('ChargeDensityCalculation','...........................')
  end subroutine dfm_cds_info_start

  subroutine dfm_cds_firsttime() bind(c, name="dfm_cds_firsttime")
    logical :: found
    integer :: k, istat
    acc_params => GetSolverParams()
    cds_nofeigs = GetInteger(acc_params, 'Number of Eigenmodes Included', found)
    if (.not. found) call Error('ChargeDensityCalculation', &
        'Number of Eigenmodes Included was not defined on the .sif file.')
    write(Message,*) 'Number of Eigenmodes included was set to ', cds_nofeigs
    call Info('ChargeDensityCalculation', Message)
    call GetConstRealArray(acc_params, cds_eigw, 'Weights of Eigen States', found)
    if (.not. found) then
      allocate(cds_eigw(cds_nofeigs, 1), stat=istat)
      if (istat /= 0) call Error('ChargeDensityCalculation', &
          'Error when allocating memory for weights of the eigen states.')
      cds_eigw = 1.0d0
      write(Message,*) 'Weights of Eigen States was not defined on', &
          'the .sif file. All weights were set to 1 .'
      call Warn('ChargeDensityCalculation', Message)
    end if
    write(Message,*) 'Weights of the Eigen States are', cds_eigw(1:cds_nofeigs, 1)
    call Info('ChargeDensityCalculation', Message)
    cds_eignum = 0
    do k = 1, acc_model % NumberOfSolvers
      acc_params => acc_model % Solvers(k) % Values
      if (GetLogical(acc_params, 'Eigen Analysis', found)) then
        if (cds_eignum == 0) then
          cds_eignum = k
        else
          write(Message,*) 'EigenSolver could not be identified.', &
              'More than one solver has Eigen Anylysis = True'
          call Error('ChargeDensityCalculation', Message)
        end if
      end if
    end do
    write(Message,*) 'Wavefunctions are variables of solver number ', cds_eignum
    call Info('ChargeDensityCalculation', Message)
  end subroutine dfm_cds_firsttime

  subroutine dfm_cds_zero_values() bind(c, name="dfm_cds_zero_values")
    acc_solver % variable % values = 0.0d0
  end subroutine dfm_cds_zero_values

  subroutine dfm_cds_begin_active(t, indexes, maxn, n, nd) &
      bind(c, name="dfm_cds_begin_active")
    integer(c_int), value :: t, maxn
    integer(c_int) :: indexes(maxn), n, nd
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes()
    nd = GetElementDOFs(indexes)
  end subroutine dfm_cds_begin_active

  subroutine dfm_cds_remap(nd, indexes, maxn) bind(c, name="dfm_cds_remap")
    integer(c_int), value :: nd, maxn
    integer(c_int) :: indexes(maxn)
    type(variable_t), pointer :: ev
    integer :: i
    ev => acc_model % Solvers(cds_eignum) % Variable
    do i = 1, nd
      indexes(i) = ev % Perm(indexes(i))
    end do
  end subroutine dfm_cds_remap

  ! The ChargeDensitySolver LocalMatrix internal procedure, transcribed
  ! statement for statement (whole-array zeroing is done by the caller,
  ! as in the original where STIFF/FORCE are zeroed on entry here).
  subroutine dfm_cds_lm(stiff, force, indexes, nmax, nd) bind(c, name="dfm_cds_lm")
    integer(c_int), value :: nmax, nd
    real(c_double) :: stiff(nmax, nmax), force(nmax)
    integer(c_int) :: indexes(nmax)
    real(c_double) :: basis(nd), detj
    real(c_double) :: localchargeatip
    logical :: stat
    integer :: t
    integer :: i, j, dim
    real(c_double) :: u, v, w
    type(variable_t), pointer :: ev
    complex(c_double), pointer :: eigv(:,:)
    dim = CoordinateSystemDimension()
    call GetElementNodes(cds_nodes)
    cds_ip = GaussPoints(acc_elem)
    do t = 1, cds_ip % n
      u = cds_ip % u(t)
      v = cds_ip % v(t)
      w = cds_ip % w(t)
      stat = ElementInfo(acc_elem, cds_nodes, u, v, w, detj, basis)
      localchargeatip = 0.0d0
      ev => acc_model % Solvers(cds_eignum) % Variable
      eigv => ev % EigenVectors
      do i = 1, cds_nofeigs
        localchargeatip = localchargeatip + cds_eigw(i, 1) * &
            sum(eigv(i, indexes(1:nd)) * basis(1:nd))**2
      end do
      do i = 1, nd
        do j = 1, nd
          stiff(i, j) = stiff(i, j) + basis(i) * basis(j) * detj * cds_ip % s(t)
        end do
        force(i) = force(i) + localchargeatip * basis(i) * detj * cds_ip % s(t)
      end do
    end do
  end subroutine dfm_cds_lm

  subroutine dfm_cds_msg_max() bind(c, name="dfm_cds_msg_max")
    integer :: n
    real(c_double) :: s
    n = acc_solver % mesh % numberofnodes
    s = maxval(acc_solver % variable % values(acc_solver % variable % perm(1:n)))
    write(Message,*) 'Greatest value of the Charge density is ', s
    call Info('ChargeDensityCalculation', Message)
  end subroutine dfm_cds_msg_max

  subroutine dfm_cds_msg_min() bind(c, name="dfm_cds_msg_min")
    integer :: n
    real(c_double) :: s
    n = acc_solver % mesh % numberofnodes
    s = minval(acc_solver % variable % values(acc_solver % variable % perm(1:n)))
    write(Message,*) 'Smallest value of the Charge density is ', s
    call Info('ChargeDensityCalculation', Message)
  end subroutine dfm_cds_msg_min

  subroutine dfm_wfs_info_start() bind(c, name="dfm_wfs_info_start")
    call Info('WavefunctionSolver','............')
    call Info('WavefunctionSolver','EigenSolver')
    call Info('WavefunctionSolver','............')
  end subroutine dfm_wfs_info_start

  subroutine dfm_wfs_firsttime() bind(c, name="dfm_wfs_firsttime")
    logical :: found
    integer :: istat
    character(len=MAX_NAME_LEN) :: str_tmp
    acc_params => GetSolverParams()
    wfs_xctype = -1
    str_tmp = GetString(acc_params, 'XC Potential type', found)
    if (trim(str_tmp) == 'perdew-zunger') then
      if (wfs_xctype /= -1) call Warn('WavefunctionSolver', &
          'Many different XC-potential types were set in .sif file.')
      if (found .and. wfs_xctype == -1) wfs_xctype = 0
    end if
    if (trim(str_tmp) == 'von barth-hedin') then
      if (wfs_xctype /= -1) call Warn('WavefunctionSolver', &
          'Many different XC-potential types were set in .sif file.')
      if (found .and. wfs_xctype == -1) wfs_xctype = 1
    end if
    if (trim(str_tmp) == 'gunnarsson-lundqvist') then
      if (wfs_xctype /= -1) call Warn('WavefunctionSolver', &
          'Many different XC-potential types were set in .sif file.')
      if (found .and. wfs_xctype == -1) wfs_xctype = 2
    end if
    if (trim(str_tmp) == 'perdew-wang') then
      if (wfs_xctype /= -1) call Warn('WavefunctionSolver', &
          'Many different XC-potential types were set in .sif file.')
      if (found .and. wfs_xctype == -1) wfs_xctype = 3
    end if
    if (trim(str_tmp) == 'none') then
      if (wfs_xctype /= -1) call Warn('WavefunctionSolver', &
          'Many different XC-potential types were set in .sif file.')
      if (found .and. wfs_xctype == -1) wfs_xctype = 4
    end if
    if (wfs_xctype == -1) call Error('WavefunctionSolver', &
        'No XC Potential type was set on the .sif file')
    if (wfs_xctype == 0) &
        call Info('WavefunctionSolver', 'Using Perdew-Zunger type XC-potential.')
    if (wfs_xctype == 1) &
        call Info('WavefunctionSolver', 'Using Von Barth-Hedin type XC-potential.')
    if (wfs_xctype == 2) &
        call Info('WavefunctionSolver', 'Using Gunnarsson-Lundqvist type XC-potential.')
    if (wfs_xctype == 3) &
        call Info('WavefunctionSolver', 'Using Perdew-Wang type XC-potential.')
    if (wfs_xctype == 4) &
        call Info('WavefunctionSolver', 'Not using any type of XC-potential.')
    wfs_shift = 0.0d0
    wfs_shift = GetConstReal(acc_params, 'Eigen System Shifting Constant', found)
    if (.not. found) then
      call Info('WavefunctionSolver', 'Eigen System Shifting constant not set. ')
    else
      print*, 'Eigen System Shifting constant was set to', wfs_shift
    end if
    wfs_nofigs = GetInteger(acc_params, 'Eigen System Values', found)
    if (.not. found) &
        call Error('WavefunctionSolver', 'Eigen System Values was not set. ')
    allocate(wfs_eigstore(wfs_nofigs), stat=istat)
    if (istat /= 0) &
        call Fatal('WavefunctionSolver', 'Memory allocation error.')
    call Info('WavefunctionSolver', 'Allocations for WavefunctionSolver are done. ')
  end subroutine dfm_wfs_firsttime

  subroutine dfm_wfs_fatal_alloc() bind(c, name="dfm_wfs_fatal_alloc")
    call Fatal('WavefunctionSolver', 'Memory allocation error.')
  end subroutine dfm_wfs_fatal_alloc

  subroutine dfm_wfs_getlocal(hart, chg, nmax) bind(c, name="dfm_wfs_getlocal")
    integer(c_int), value :: nmax
    real(c_double) :: hart(nmax), chg(nmax)
    call GetScalarLocalSolution(hart(1:nmax), 'Potential')
    call GetScalarLocalSolution(chg(1:nmax), 'Charge density')
  end subroutine dfm_wfs_getlocal

  ! The WaveFunctionSolver LocalMatrix internal procedure, transcribed
  ! statement for statement.
  subroutine dfm_wfs_lm(stiff, mass, hart, chg, nmax, nd) bind(c, name="dfm_wfs_lm")
    integer(c_int), value :: nmax, nd
    real(c_double) :: stiff(nmax, nmax), mass(nmax, nmax)
    real(c_double) :: hart(nmax), chg(nmax)
    real(c_double) :: basis(nd), dbasisdx(nd, 3), detj
    real(c_double) :: localpotentialatip, localchargeatip, localxcpotentialatip
    logical :: stat
    integer :: t
    integer :: i, j, k, dim
    real(c_double) :: u, v, w
    dim = CoordinateSystemDimension()
    call GetElementNodes(wfs_nodes)
    wfs_ip = GaussPoints(acc_elem)
    do t = 1, wfs_ip % n
      u = wfs_ip % u(t)
      v = wfs_ip % v(t)
      w = wfs_ip % w(t)
      stat = ElementInfo(acc_elem, wfs_nodes, u, v, w, detj, basis, dbasisdx)
      localpotentialatip = sum(basis(1:nd) * hart(1:nd))
      localchargeatip = sum(basis(1:nd) * chg(1:nd))
      localchargeatip = max(localchargeatip, 0.0_dp)
      if (.not. wfs_xctype == 4) then
        localxcpotentialatip = uxc(localchargeatip, 0.0d0, 1, wfs_xctype)
      else
        localxcpotentialatip = 0.0d0
      end if
      do i = 1, nd
        do j = 1, nd
          do k = 1, dim
            stiff(i, j) = stiff(i, j) &
                + 0.5 * dbasisdx(i, k) * dbasisdx(j, k) * detj * wfs_ip % s(t)
          end do
          mass(i, j) = mass(i, j) + basis(i) * basis(j) * detj * wfs_ip % s(t)
          stiff(i, j) = stiff(i, j) &
              + (localpotentialatip + localxcpotentialatip) &
              * basis(i) * basis(j) * detj * wfs_ip % s(t)
        end do
      end do
    end do
    stiff(1:nmax, 1:nmax) = 0.5d0 * (stiff(1:nmax, 1:nmax) &
        + transpose(stiff(1:nmax, 1:nmax)))
    mass(1:nmax, 1:nmax) = 0.5d0 * (mass(1:nmax, 1:nmax) &
        + transpose(mass(1:nmax, 1:nmax)))
  end subroutine dfm_wfs_lm

  subroutine dfm_wfs_maybe_shift() bind(c, name="dfm_wfs_maybe_shift")
    if (.not. abs(wfs_shift) < epsilon(wfs_shift)) &
        acc_solver % matrix % values = acc_solver % matrix % values &
        + wfs_shift * acc_solver % matrix % massvalues
  end subroutine dfm_wfs_maybe_shift

  subroutine dfm_wfs_eig_store() bind(c, name="dfm_wfs_eig_store")
    wfs_eigstore(1:wfs_nofigs) = acc_solver % variable % eigenvalues(1:wfs_nofigs)
  end subroutine dfm_wfs_eig_store

  subroutine dfm_wfs_eig_head() bind(c, name="dfm_wfs_eig_head")
    call Info('WavefunctionSolver','Absolute change of the eigenvalues:')
  end subroutine dfm_wfs_eig_head

  subroutine dfm_wfs_eig_report() bind(c, name="dfm_wfs_eig_report")
    integer :: t
    do t = 1, wfs_nofigs
      write(Message,*) 'Eigenvalue', t, 'has absolute change', &
          abs(acc_solver % variable % eigenvalues(t) - wfs_eigstore(t))
      call Info('WavefunctionSolver', Message)
    end do
  end subroutine dfm_wfs_eig_report

  ! -- PoissonSolver ------------------------------------------------
  subroutine dfm_ps_info_start() bind(c, name="dfm_ps_info_start")
    call Info('PoissonSolver','------------------')
    call Info('PoissonSolver',' Poisson equation ')
    call Info('PoissonSolver','------------------')
  end subroutine dfm_ps_info_start

  subroutine dfm_ps_firsttime() bind(c, name="dfm_ps_firsttime")
    logical :: found
    integer :: k, v, m, istat, gs
    real(c_double) :: x, y, z, tol
    character(len=MAX_NAME_LEN) :: str_tmp
    acc_params => GetSolverParams()
    ps_cdsnum = 0
    do k = 1, acc_model % NumberOfSolvers
      acc_params => acc_model % Solvers(k) % Values
      v = GetInteger(acc_params, 'Number of Eigenmodes Included', found)
      if (found) then
        if (ps_cdsnum == 0) then
          ps_cdsnum = k
        else
          write(Message,*) 'ChargeDensitySolver could not be identified.', &
              'More than one solver has Number of Eigenmodes Included', &
              'defined.'
          call Error('PoissonSolver', Message)
        end if
      end if
    end do
    if (ps_cdsnum == 0) then
      write(Message,*) 'ChargeDnsitySolver could not be identified.', &
          'No solver has Number of Eigenmodes Included', &
          'defined.'
      call Fatal('PoissonSolver', Message)
    end if
    acc_params => GetSolverParams()
    ps_perm => acc_solver % variable % perm
    ps_nofnuclei = 0
    ps_nofnuclei = GetInteger(acc_params, 'NOFnuclei', found)
    if (.not. found) call Warn('PoissonSolver', &
        'NOFnuclei was not set on the sif file.')
    write(Message,*) 'NOFnuclei was set to ', ps_nofnuclei
    call Info('PoissonSolver', Message)
    if (.not. (ps_nofnuclei == 0)) then
      allocate(ps_nuctab(ps_nofnuclei, 4), ps_nucnodes(ps_nofnuclei))
      ps_nuctab = 0.0d0
      ps_nucnodes = 0
      call GetConstRealArray(acc_params, ps_nuctab, 'NucleiTable', found)
      if (.not. found) call Warn('PoissonSolver', &
          'NucleiTable was not found from .sif file.')
      x = 0.0d0
      y = 0.0d0
      z = 0.0d0
      tol = 10 * epsilon(x)
      do k = 1, acc_solver % mesh % numberofnodes
        x = acc_solver % mesh % nodes % x(k)
        y = acc_solver % mesh % nodes % y(k)
        z = acc_solver % mesh % nodes % z(k)
        do m = 1, ps_nofnuclei
          if (abs(ps_nuctab(m, 2) - x) < tol .and. &
              abs(ps_nuctab(m, 3) - y) < tol .and. &
              abs(ps_nuctab(m, 4) - z) < tol) ps_nucnodes(m) = k
        end do
      end do
    end if
    write(Message,*) 'Nuclei are located on nodes ', ps_nucnodes
    call Info('PoissonSolver', Message)
    ps_totalnuc = 0.0d0
    do v = 1, ps_nofnuclei
      ps_totalnuc = ps_totalnuc + ps_nuctab(v, 1)
    end do
    gs = size(acc_solver % variable % values)
    ps_globsize = gs
    allocate(ps_old(gs), stat=istat)
    if (istat /= 0) call Error('PoissonSolver', &
        'Error when allocating memory for previous solution storage.')
    ps_old = 0.0d0
    ps_relaxtype = 0
    ps_pa = -100000.0d0
    ps_pb = -100000.0d0
    ps_pc = -100000.0d0
    ps_nofiter = 0
    str_tmp = GetString(acc_params, 'Relaxation Method', found)
    if (.not. found) call Warn('PoissonSolver', &
        'No Relaxation method was found the .sif file, using none.')
    ps_pa = GetConstReal(acc_params, 'Relaxation Parameter A', found)
    if (.not. found) call Info('PoissonSolver', &
        'Relaxation Parameter A was not defined' // &
        'on the .sif file.')
    ps_pb = GetConstReal(acc_params, 'Relaxation Parameter B', found)
    if (.not. found) call Info('PoissonSolver', &
        'Relaxation Parameter B was not defined' // &
        'on the .sif file.')
    ps_pc = GetConstReal(acc_params, 'Relaxation Parameter C', found)
    if (.not. found) call Info('PoissonSolver', &
        'Relaxation Parameter C was not defined' // &
        'on the .sif file.')
    if (trim(str_tmp) == 'constant mixing') then
      ps_relaxtype = 1
      write(Message,*) 'Mixing Relaxation Factor (RF) is constant', ps_pa
      call Info('PoissonSolver', Message)
      call Info('PoissonSolver', &
          'Relaxation method is r(k+1) = (1 - RF) * r(k) + RF * r(k+1)')
    end if
    if (trim(str_tmp) == 'exponential mixing') then
      ps_relaxtype = 2
      write(Message,*) 'Mixing Relaxation Factor (RF) is varying by scheme', &
          ' RF(k) = ', ps_pc, ' + 1 - ', ps_pa, '* Exp( -', ps_pb, ' * k )'
      call Info('PoissonSolver', Message)
      call Info('PoissonSolver', &
          'Relaxation method is r(k+1) = (1 - RF) * r(k) + RF * r(k+1)')
    end if
    ps_startrf = 1.0d8
    ps_startafter = 100000000
    ps_startrf = GetConstReal(acc_params, &
        'Start GRPulay if relaxation factor is more than', found)
    ps_startafter = GetInteger(acc_params, &
        'Start GRPulay after iterations', found)
    write(Message,*) 'Starting GRPulay after ', ps_startafter, 'iterations.'
    call Info('PoissonSolver', Message)
    write(Message,*) 'Starting GRPulay when RF is more than ', ps_startrf
    call Info('PoissonSolver', Message)
  end subroutine dfm_ps_firsttime

  subroutine dfm_ps_fatal_alloc() bind(c, name="dfm_ps_fatal_alloc")
    call Fatal('PoissonSolver', 'Memory allocation error.')
  end subroutine dfm_ps_fatal_alloc

  function dfm_ps_relaxtype() result(r) bind(c, name="dfm_ps_relaxtype")
    integer(c_int) :: r
    r = ps_relaxtype
  end function dfm_ps_relaxtype

  function dfm_ps_calculating() result(r) bind(c, name="dfm_ps_calculating")
    integer(c_int) :: r
    r = merge(1, 0, ps_calc)
  end function dfm_ps_calculating

  subroutine dfm_ps_save_old() bind(c, name="dfm_ps_save_old")
    if (.not. ps_relaxtype == 3) &
        ps_old(1:ps_globsize) = acc_solver % variable % values(1:ps_globsize)
  end subroutine dfm_ps_save_old

  subroutine dfm_ps_begin_active(t, n, nd) bind(c, name="dfm_ps_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n, nd
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes(acc_elem)
    nd = GetElementNOFDOFs(acc_elem)
  end subroutine dfm_ps_begin_active

  subroutine dfm_wfs_begin_active(t, n, nd) bind(c, name="dfm_wfs_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n, nd
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes(acc_elem)
    nd = GetElementNOFDOFs(acc_elem)
  end subroutine dfm_wfs_begin_active

  subroutine dfm_ps_getlocal(chg, nmax) bind(c, name="dfm_ps_getlocal")
    integer(c_int), value :: nmax
    real(c_double) :: chg(nmax)
    call GetScalarLocalSolution(chg(1:nmax), 'Charge density')
  end subroutine dfm_ps_getlocal

  ! The PoissonSolver LocalMatrix internal procedure, transcribed
  ! statement for statement. TotalCharge accumulates across elements
  ! through the caller's scalar, exactly as the original's local.
  subroutine dfm_ps_lm(stiff, force, chg, nmax, nd, total) &
      bind(c, name="dfm_ps_lm")
    integer(c_int), value :: nmax, nd
    real(c_double) :: stiff(nmax, nmax), force(nmax)
    real(c_double) :: chg(nmax), total
    real(c_double) :: basis(nd), dbasisdx(nd, 3), detj, loadatip
    real(c_double) :: localchargedensityatip
    logical :: stat
    integer :: i, j, k, t, dim
    real(c_double) :: u, v, w
    dim = CoordinateSystemDimension()
    call GetElementNodes(ps_nodes)
    ps_ip = GaussPoints(acc_elem)
    do t = 1, ps_ip % n
      u = ps_ip % u(t)
      v = ps_ip % v(t)
      w = ps_ip % w(t)
      stat = ElementInfo(acc_elem, ps_nodes, u, v, w, detj, basis, dbasisdx)
      localchargedensityatip = sum(basis(1:nd) * chg(1:nd))
      localchargedensityatip = max(localchargedensityatip, 0.0_dp)
      total = total + localchargedensityatip * detj * ps_ip % s(t)
      do i = 1, nd
        force(i) = force(i) &
            + 4 * Pi * localchargedensityatip * basis(i) * detj * ps_ip % s(t)
        do j = 1, nd
          do k = 1, dim
            stiff(i, j) = stiff(i, j) &
                + dbasisdx(i, k) * dbasisdx(j, k) * detj * ps_ip % s(t)
          end do
        end do
      end do
    end do
    stiff(1:nmax, 1:nmax) = 0.5d0 * (stiff(1:nmax, 1:nmax) &
        + transpose(stiff(1:nmax, 1:nmax)))
  end subroutine dfm_ps_lm

  subroutine dfm_ps_normalize(total) bind(c, name="dfm_ps_normalize")
    real(c_double), value :: total
    if (total > AEPS) then
      acc_solver % matrix % rhs = (ps_totalnuc / total) &
          * acc_solver % matrix % rhs
      acc_model % Solvers(ps_cdsnum) % Variable % Values = &
          (ps_totalnuc / total) &
          * acc_model % Solvers(ps_cdsnum) % Variable % Values
    end if
  end subroutine dfm_ps_normalize

  subroutine dfm_ps_add_nuclei() bind(c, name="dfm_ps_add_nuclei")
    integer :: v
    do v = 1, ps_nofnuclei
      if (ps_nucnodes(v) == 0) cycle
      if (.not. ps_perm(ps_nucnodes(v)) == 0) then
        acc_solver % matrix % rhs(ps_perm(ps_nucnodes(v))) = &
            acc_solver % matrix % rhs(ps_perm(ps_nucnodes(v))) &
            - 4 * Pi * ps_nuctab(v, 1)
      else
        acc_solver % matrix % rhs(ps_nucnodes(v)) = &
            acc_solver % matrix % rhs(ps_nucnodes(v)) &
            - 4 * Pi * ps_nuctab(v, 1)
      end if
    end do
  end subroutine dfm_ps_add_nuclei

  subroutine dfm_ps_info_norelax() bind(c, name="dfm_ps_info_norelax")
    call Info('PoissonSolver','No relaxation method applied')
  end subroutine dfm_ps_info_norelax

  subroutine dfm_ps_info_grp() bind(c, name="dfm_ps_info_grp")
    call Info('PoissonSolver','Using GRPulay.')
  end subroutine dfm_ps_info_grp

  subroutine dfm_ps_mix_const() bind(c, name="dfm_ps_mix_const")
    ps_rf = ps_pa
    acc_solver % variable % values(1:ps_globsize) = &
        (1 - ps_rf) * ps_old(1:ps_globsize) &
        + ps_rf * acc_solver % variable % values(1:ps_globsize)
    write(Message,*) 'Applied constant mixing with relaxation factor', ps_rf
    call Info('PoissonSolver', Message)
  end subroutine dfm_ps_mix_const

  subroutine dfm_ps_mix_exp() bind(c, name="dfm_ps_mix_exp")
    ps_rf = min(ps_pc + 1.0d0 - ps_pa * exp(-ps_pb * ps_nofiter), 1.0d0)
    acc_solver % variable % values(1:ps_globsize) = &
        (1 - ps_rf) * ps_old(1:ps_globsize) &
        + ps_rf * acc_solver % variable % values(1:ps_globsize)
    write(Message,*) 'Applied exponential mixing with relaxation factor', ps_rf
    call Info('PoissonSolver', Message)
  end subroutine dfm_ps_mix_exp

  ! GRPULAY3 calls plus the StoreOnly/RelaxType/NOFiter gate, exactly as
  ! the original's post-solve block.
  subroutine dfm_ps_grpulay_gate() bind(c, name="dfm_ps_grpulay_gate")
    call grpulay3_body()
    if ((ps_rf > ps_startrf) .or. (ps_nofiter == ps_startafter)) &
        ps_storeonly = .FALSE.
    if (.not. ps_storeonly) ps_relaxtype = 3
    ps_nofiter = ps_nofiter + 1
  end subroutine dfm_ps_grpulay_gate

  subroutine dfm_ps_relchange() bind(c, name="dfm_ps_relchange")
    write(Message,*) 'Relative change between relaxed solutions is', &
        sqrt(dot_product( &
        ps_old(1:ps_globsize) - acc_solver % variable % values(1:ps_globsize), &
        ps_old(1:ps_globsize) - acc_solver % variable % values(1:ps_globsize)) &
        / &
        dot_product( &
        acc_solver % variable % values(1:ps_globsize), &
        acc_solver % variable % values(1:ps_globsize)))
    call Info('PoissonSolver', Message)
    ps_old(1:ps_globsize) = acc_solver % variable % values(1:ps_globsize)
  end subroutine dfm_ps_relchange

  subroutine dfm_ps_msg_total(total) bind(c, name="dfm_ps_msg_total")
    real(c_double), value :: total
    write(Message,*) 'Total charge of electrons before normalization was', &
        total
    call Info('PoissonSolver', Message)
  end subroutine dfm_ps_msg_total

  subroutine dfm_ps_msg_nuctotal() bind(c, name="dfm_ps_msg_nuctotal")
    write(Message,*) 'Total charge of nuclei ', ps_totalnuc
    call Info('PoissonSolver', Message)
  end subroutine dfm_ps_msg_nuctotal

  subroutine dfm_ps_msg_min() bind(c, name="dfm_ps_msg_min")
    integer :: n
    real(c_double) :: s
    n = acc_solver % mesh % numberofnodes
    s = minval(acc_solver % variable % values(ps_perm(1:n)))
    write(Message,*) 'Smallest value of the potential is ', s, &
        ' at node ', minloc(acc_solver % variable % values(ps_perm(1:n)))
    call Info('PoissonSolver', Message)
  end subroutine dfm_ps_msg_min

  subroutine dfm_ps_msg_max() bind(c, name="dfm_ps_msg_max")
    integer :: n
    real(c_double) :: s
    n = acc_solver % mesh % numberofnodes
    s = maxval(acc_solver % variable % values(ps_perm(1:n)))
    write(Message,*) 'Greatest value of the potential is ', s, &
        ' at node ', maxloc(acc_solver % variable % values(ps_perm(1:n)))
    call Info('PoissonSolver', Message)
  end subroutine dfm_ps_msg_max

  ! -- GRPULAY3 (PoissonSolver CONTAINS internal) --------------------
  ! Transcribed statement for statement; serial branches only, exactly
  ! as the original (its MPI alternatives are commented out there).
  subroutine grpulay3_body()
    integer :: istat, t, p
    real(c_double) :: s
    if (ps_gfirst) then
      allocate(ps_pulay(4, 4), &
          ps_resstor(3, ps_globsize), &
          ps_coeffs(4), &
          ps_solstor(3, ps_globsize), &
          ps_rhs(4), &
          stat=istat)
      if (istat /= 0) call Error('PoissonSolver, GRPulay3', &
          'Memory allocation failed')
      ps_nofstored = 0
      ps_solstor = 0.0d0
      ps_gfirst = .FALSE.
    end if
    if (.not. ps_donenonrelax) then
      ps_solstor(2:3, 1:ps_globsize) = ps_solstor(1:2, 1:ps_globsize)
      ps_solstor(1, 1:ps_globsize) = &
          acc_solver % variable % values(1:ps_globsize)
      ps_nofstored = ps_nofstored + 1
    end if
    if (.not. ps_storeonly) then
      if (ps_nofstored < 3) call Error('PoissonSolver, GRPulay3', &
          'Less than 3 solutions stored')
      if (.not. ps_donenonrelax) then
        ps_donenonrelax = .TRUE.
        call Info('PoissonSolver, GRPulay3', &
            'Calculating the non-relaxed solution on the next step')
      else
        ps_resstor(1, 1:ps_globsize) = &
            acc_solver % variable % values(1:ps_globsize) &
            - ps_solstor(1, 1:ps_globsize)
        do t = 2, 3
          ps_resstor(t, 1:ps_globsize) = &
              ps_solstor(t - 1, 1:ps_globsize) &
              - ps_solstor(t, 1:ps_globsize)
        end do
        ps_pulay = 0.0d0
        ps_rhs = 0.0d0
        ps_rhs(4) = 1.0d0
        do t = 1, 3
          ps_pulay(t, 4) = -0.5d0
          ps_pulay(4, t) = 1.0d0
          do p = 1, 3
            ps_pulay(t, p) = dot_product( &
                ps_resstor(t, 1:ps_globsize), ps_resstor(p, 1:ps_globsize))
          end do
        end do
        call LUSolve(4, ps_pulay, ps_rhs)
        ps_coeffs(1:4) = ps_rhs(1:4)
        write(Message,*) 'Calculated relaxations coefficients are', &
            ps_coeffs(1:3)
        call Info('PoissonSolver, GRPulay3', Message)
        write(Message,*) 'Lagrange multiplier of the minimization problem', &
            ps_coeffs(4)
        call Info('PoissonSolver, GRPulay3', Message)
        write(Message,*) 'Norm of the minimized residual', &
            sqrt(dot_product( &
            acc_solver % variable % values(1:ps_globsize) &
            - ps_solstor(1, 1:ps_globsize), &
            acc_solver % variable % values(1:ps_globsize) &
            - ps_solstor(1, 1:ps_globsize)))
        call Info('PoissonSolver, GRPulay3', Message)
        acc_solver % variable % values = 0.0d0
        do t = 1, 3
          acc_solver % variable % values(1:ps_globsize) = &
              acc_solver % variable % values(1:ps_globsize) &
              + ps_coeffs(t) * ps_solstor(t, 1:ps_globsize)
        end do
        ps_donenonrelax = .FALSE.
      end if
    end if
    ps_calc = .NOT. ps_donenonrelax
  end subroutine grpulay3_body

end module dftm_acc
