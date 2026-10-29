! T022-A HelmholtzSolve port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/HelmholtzSolve.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays, C_INT
! status codes). All Fortran-object traffic — TYPE(...) handles,
! OPTIONAL/POINTER dummies, hidden lengths, assumed-shape descriptors,
! complex reductions, list-directed timing reports — stays inside
! Fortran, so the Rust side is pure C calls. The wrapper module mirrors
! HelmholtzSolver's own locals (current element, per-internal-procedure
! SAVE Nodes/IP state, entry SAVE arrays and flags) with identical
! lifetimes. All floating-point work (the Galerkin assembly, complex
! admittance/convection algebra, MATMUL/SUM reductions, the CPUTime-based
! reports) stays in Fortran (same compiler), so bit-identity holds by
! construction wherever the call sequence matches.
!
! Faithfulness notes (mirrors HelmholtzSolve.F90 statement for statement):
! - SAVE semantics replicated: the entry SAVE block (STIFF/PREC/FORCE/
!   Work/Load/SoundSpeed/Density/Damping/Impedance/ConvVelo/WallVelocity/
!   AllocationsDone) is allocated once at first-mesh sizing
!   (MaxElementDOFs) and reallocated only when `Solver % MeshChanged`,
!   exactly as the original; the three internal procedures' SAVE Nodes
!   live as wrapper-side SAVE state with identical lifetimes (consumed
!   before the next fill, as in the original).
! - Allocation failure calls the same `Fatal('HelmholtzSolve', ...)` via
!   the `STAT=istat` check, instead of aborting.
! - Pure-query touches with no observable effect are replicated as
!   touches where they could trap (`PresPerm` association,
!   `Solver % Variable % Norm` read) and documented as skipped where
!   they cannot (`GetSimulation()`, `GetEquation()`,
!   `GetElementNOFNodes()` pre-loop sizing read whose value is
!   overwritten; all are side-effect-free queries).
! - `Model`, `dt`, `TransientSimulation` are unused by the entry body
!   itself, exactly as in the original (received and ignored), except
!   for `Model % Simulation` in the `res: Frequency` list update.
! - No Rust threads: the solver entries run serially, so the module SAVE
!   state needs no guard.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module helmholtzsolve_acc
  use iso_c_binding
  use coordinatesystems
  use defutils
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(valuelist_t), pointer, save :: acc_params => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(valuelist_t), pointer, save :: acc_mat => null()
  type(valuelist_t), pointer, save :: acc_bf => null()
  type(valuelist_t), pointer, save :: acc_bc => null()
  type(valuelist_t), pointer, save :: acc_eq => null()
  integer, pointer, save :: acc_presperm(:) => null()
  ! Entry SAVE block (mirrors HelmholtzSolver's SAVE statement).
  complex(kind=dp), allocatable, save :: acc_stiff(:,:), acc_prec(:,:), acc_force(:)
  real(kind=dp), allocatable, save :: acc_work(:), acc_load(:,:), acc_snd(:)
  real(kind=dp), allocatable, save :: acc_rho(:), acc_damp(:), acc_imp(:,:)
  real(kind=dp), allocatable, save :: acc_conv(:,:)
  complex(kind=dp), allocatable, save :: acc_wall(:,:)
  ! Internal-procedure scratch (sequential lifetimes, sized 2*N like STIFF).
  real(kind=dp), allocatable, save :: acc_basis(:), acc_dbasis(:,:)
  logical, save :: acc_alloc = .false.
  ! Entry scalars (mirrors the solver body's locals).
  integer(c_int), save :: acc_niter = 0
  logical, save :: acc_bubbles = .false.
  logical, save :: acc_useconv = .false.
  logical, save :: acc_usedens = .true.
  logical, save :: acc_useprecshift = .false.
  integer(c_int), save :: acc_dim = 0
  real(kind=dp), save :: acc_omega = 0.0_dp
  complex(kind=dp), save :: acc_shift = (0.0_dp, 0.0_dp)
  logical, save :: acc_anyflow = .false.
  logical, save :: acc_anystruct = .false.
  real(kind=dp), save :: acc_norm = 0.0_dp
  real(kind=dp), save :: acc_at = 0.0_dp
  real(kind=dp), save :: acc_totat = 0.0_dp
  real(kind=dp), save :: acc_totst = 0.0_dp
  ! Flow/structure interface state (mirrors the entry locals).
  type(variable_t), pointer, save :: acc_flowsol => null()
  type(variable_t), pointer, save :: acc_dispsol => null()
  type(variable_t), pointer, save :: acc_dispsolim => null()
  real(kind=dp), pointer, save :: acc_flow(:) => null()
  real(kind=dp), pointer, save :: acc_disp(:) => null()
  real(kind=dp), pointer, save :: acc_dispim(:) => null()
  complex(kind=dp), pointer, save :: acc_dispeig(:) => null()
  integer, pointer, save :: acc_flowperm(:) => null()
  integer, pointer, save :: acc_dispperm(:) => null()
  integer(c_int), save :: acc_flowdofs = 0
  integer(c_int), save :: acc_dispdofs = 0
  integer(c_int), save :: acc_noeigen = 0
  integer(c_int), save :: acc_dismode = 0
  ! LocalMatrix context (mirrors its SAVE Nodes; IP is a fresh local
  ! there, kept module-side with identical per-element lifetime).
  type(nodes_t), target, save :: acc_nodes_lm
  type(gaussintegrationpoints_t), save :: acc_ip_lm
  ! LocalMatrixBoundary context (mirrors its SAVE Nodes likewise).
  type(nodes_t), target, save :: acc_nodes_lb
  type(gaussintegrationpoints_t), save :: acc_ip_lb
  ! LocalInterfaceMatrix context (mirrors its SAVE Nodes likewise).
  type(nodes_t), target, save :: acc_nodes_li
  type(gaussintegrationpoints_t), save :: acc_ip_li

contains

  ! -- entry state ----------------------------------------------------
  subroutine hhz_begin(h_model, h_solver) bind(c, name="hhz_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine hhz_begin

  ! Allocation block: `IF (.NOT. AllocationsDone .OR. Solver % MeshChanged)`.
  subroutine hhz_ensure_alloc() bind(c, name="hhz_ensure_alloc")
    integer :: n, istat
    if (.not. acc_alloc .or. acc_solver % meshchanged) then
      n = acc_solver % mesh % maxelementdofs
      if (acc_alloc) then
        deallocate(acc_imp, acc_work, acc_force, acc_stiff, acc_prec, &
            acc_snd, acc_rho, acc_conv, acc_damp, acc_load, acc_wall, &
            acc_basis, acc_dbasis)
      end if
      allocate(acc_imp(2,n), acc_work(n), acc_force(2*n), &
          acc_stiff(2*n,2*n), acc_prec(2*n,2*n), &
          acc_snd(n), acc_rho(n), acc_conv(3,n), acc_damp(n), &
          acc_load(2,n), acc_wall(3,n), &
          acc_basis(2*n), acc_dbasis(2*n,3), stat=istat)
      if (istat /= 0) then
        call Fatal('HelmholtzSolve', 'Memory allocation error.')
      end if
      acc_alloc = .true.
    end if
  end subroutine hhz_ensure_alloc

  ! Lines 158-296: solver params, density flag, frequency/interface
  ! setup, preconditioner shift, convection check, DefaultStart.
  ! Returns NonlinearIter (the Rust iter-loop bound).
  function hhz_pre() result(niter) bind(c, name="hhz_pre")
    integer(c_int) :: niter
    type(valuelist_t), pointer :: params
    type(variable_t), pointer :: flowsol, dispsol, dispsolim
    real(kind=dp), pointer :: flow(:), disp(:), dispim(:)
    complex(kind=dp), pointer :: dispeig(:)
    integer, pointer :: flowperm(:), dispperm(:)
    integer :: flowdofs, dispdofs, dispmode, noigen, dim, ntmp
    character(len=MAX_NAME_LEN) :: varname
    logical :: found, gotfreq, flowif, structif
    complex(kind=dp) :: shift
    real(kind=dp) :: omega
    acc_params => GetSolverParams()
    params => acc_params
    niter = GetInteger(params, 'Nonlinear System Max Iterations', found)
    if (.not. found) niter = 1
    acc_niter = niter
    acc_bubbles = GetLogical(params, 'Bubbles', found)
    acc_usedens = .true.
    call ListAddLogical(params, 'Use Density', acc_usedens)
    ntmp = GetElementNOFNodes()
    acc_dim = CoordinateSystemDimension()
    dim = acc_dim
    gotfreq = .false.
    acc_wall = 0.0_dp
    acc_presperm => acc_solver % variable % perm
    flowif = ListCheckPresentAnyBC(acc_model, 'Flow Interface')
    structif = ListCheckPresentAnyBC(acc_model, 'Structure Interface')
    acc_anyflow = flowif
    acc_anystruct = structif
    if (flowif) then
      varname = GetString(params, 'Velocity Variable Name', found)
      if (.not. found) varname = 'flow'
      flowsol => VariableGet(acc_solver % mesh % variables, varname)
      if (.not. associated(flowsol)) then
        call Fatal('HelmholtzSolver', 'No flow variable associated:'//trim(varname))
      end if
      flow => flowsol % values
      flowperm => flowsol % perm
      flowdofs = flowsol % dofs
      if (flowdofs /= 2*dim) then
        call Fatal('HelmholtzSolver', 'Harmonic flow field should have 2*dim components')
      end if
      acc_flowsol => flowsol
      acc_flow => flow
      acc_flowperm => flowperm
      acc_flowdofs = flowdofs
    end if
    if (structif) then
      varname = GetString(params, 'Displacement Variable Name', found)
      if (.not. found) varname = 'displacement'
      dispsol => VariableGet(acc_solver % mesh % variables, varname)
      if (.not. associated(dispsol)) then
        call Fatal('HelmholtzSolver', 'No displacement variable associated:'//trim(varname))
      end if
      disp => dispsol % values
      dispperm => dispsol % perm
      dispdofs = dispsol % dofs
      dispmode = -1
      noigen = GetInteger(params, 'Displacement Variable EigenMode', found)
      if (noigen > 0) then
        if (noigen > size(dispsol % eigenvalues)) then
          call Fatal('HelmholtzSolver', 'Requested eigenmode does not exist')
        end if
        dispeig => dispsol % eigenvectors(noigen,:)
        if (GetLogical(params, 'Displacement Variable Frequency', found)) then
          omega = sqrt(dispsol % eigenvalues(noigen))
          gotfreq = .true.
          acc_omega = omega
        end if
        if (dispdofs < dim) then
          call Fatal('HelmholtzSolver', 'Eigenmode displacement field should have at least 1*dim components')
        end if
        dispmode = 1
        acc_dispeig => dispeig
      else
        if (dispdofs == 2*dim) then
          dispmode = 2
          continue
        else if (dispdofs == dim) then
          dispsolim => VariableGet(acc_solver % mesh % variables, trim(varname)//' im')
          dispmode = 3
          if (associated(dispsolim)) then
            dispim => dispsolim % values
            dispmode = 4
            acc_dispsolim => dispsolim
            acc_dispim => dispim
          end if
        else
          call Fatal('HelmholtzSolver', 'Harmonic displacement field should have n*dim components')
        end if
      end if
      if (dispmode == -1) then
        call Fatal('HelmholtzSolver', 'Could not determine displacement mode!')
      end if
      acc_dispsol => dispsol
      acc_disp => disp
      acc_dispperm => dispperm
      acc_dispdofs = dispdofs
      acc_noeigen = noigen
      acc_dismode = dispmode
    end if
    if (.not. gotfreq) then
      omega = GetAngularFrequency(found=gotfreq)
      acc_omega = omega
    end if
    if (.not. gotfreq) then
      call Fatal('HelmholtzSolver', 'Could not figure out Frequency!')
    end if
    call ListAddConstReal(acc_model % simulation, 'res: Frequency', acc_omega/(2*PI))
    shift = GetCReal(params, 'Linear System Preconditioning Damp Coefficient', acc_useprecshift)
    shift = cmplx(real(shift), &
        GetCReal(params, 'Linear System Preconditioning Damp Coefficient Im', found), kind=dp)
    acc_shift = shift
    acc_useprecshift = acc_useprecshift .or. found
    if (acc_useprecshift) then
      if (ListGetString(params, 'Linear System Solver', found) /= 'iterative') then
        call Warn('HelmholtzSolver', 'Damped preconditioning makes sence for iterative methods, canceling!')
        acc_useprecshift = .false.
      else
        call Info('HelmholtzSolver', 'Generating special precondining matrix', level=12)
      end if
    end if
    acc_conv = 0.0_dp
    acc_useconv = ListCheckPrefixAnyMaterial(acc_model, 'Convection Velocity')
    acc_norm = acc_solver % variable % norm
    acc_totat = 0.0d0
    acc_totst = 0.0d0
    call DefaultStart()
  end function hhz_pre

  function hhz_need_iface() result(r) bind(c, name="hhz_need_iface")
    integer(c_int) :: r
    r = merge(1, 0, acc_anystruct .or. acc_anyflow)
  end function hhz_need_iface

  ! Iter-loop head: timing, Level=4 progress block, DefaultInitialize,
  ! StartAdvanceOutput (lines 300-318).
  subroutine hhz_iter_begin(iter) bind(c, name="hhz_iter_begin")
    integer(c_int), value :: iter
    acc_at = CPUTime()
    call Info('HelmholtzSolve', ' ', level=4)
    call Info('HelmholtzSolve', '-------------------------------------', level=4)
    write(message, *) 'Helmholtz iteration', iter
    call Info('HelmholtzSolve', message, level=4)
    write(message, *) 'Frequency (Hz): ', acc_omega/(2*PI)
    call Info('HelmholtzSolve', message, level=4)
    call Info('HelmholtzSolve', '-------------------------------------', level=4)
    call Info('HelmholtzSolve', ' ', level=4)
    call Info('HelmholtzSolve', 'Starting Assembly', level=4)
    call DefaultInitialize()
    call StartAdvanceOutput('HelmholtzSolve', 'Assembly:')
  end subroutine hhz_iter_begin

  function hhz_nactive() result(n) bind(c, name="hhz_nactive")
    integer(c_int) :: n
    n = acc_solver % numberofactiveelements
  end function hhz_nactive

  ! Bulk element body (lines 323-362): AdvanceOutput, active element,
  ! material/body-force reads, LocalMatrix, DefaultUpdateEquations.
  subroutine hhz_bulk(t) bind(c, name="hhz_bulk")
    integer(c_int), value :: t
    integer :: n, nd
    logical :: found
    call AdvanceOutput(t, acc_solver % numberofactiveelements)
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
    acc_eq => GetEquation()
    acc_mat => GetMaterial()
    acc_damp(1:n) = GetReal(acc_mat, 'Sound Damping', found)
    acc_snd(1:n) = ListGetReal(acc_mat, 'Sound Speed', &
        n, acc_elem % nodeindexes, unfoundfatal=.true.)
    acc_rho(1:n) = ListGetReal(acc_mat, 'Density', &
        n, acc_elem % nodeindexes, unfoundfatal=.true.)
    if (acc_useconv) then
      acc_conv(1,1:n) = GetReal(acc_mat, 'Convection Velocity 1', found)
      acc_conv(2,1:n) = GetReal(acc_mat, 'Convection Velocity 2', found)
      if (acc_dim == 3) then
        acc_conv(3,1:n) = GetReal(acc_mat, 'Convection Velocity 3', found)
      end if
    end if
    acc_bf => GetBodyForce()
    acc_load(1,1:n) = GetReal(acc_bf, 'Pressure Source 1', found)
    acc_load(2,1:n) = GetReal(acc_bf, 'Pressure Source 2', found)
    call hhz_local_matrix(n, nd)
    call DefaultUpdateEquations(acc_stiff, acc_force)
  end subroutine hhz_bulk

  subroutine hhz_finish_bulk() bind(c, name="hhz_finish_bulk")
    call DefaultFinishBulkAssembly()
  end subroutine hhz_finish_bulk

  function hhz_nbnd() result(n) bind(c, name="hhz_nbnd")
    integer(c_int) :: n
    n = acc_solver % mesh % numberofboundaryelements
  end function hhz_nbnd

  ! Boundary element body (lines 371-407): Neumann/Newton BCs.
  subroutine hhz_bnd(t) bind(c, name="hhz_bnd")
    integer(c_int), value :: t
    integer :: n, nd
    logical :: found, planewave
    acc_elem => GetBoundaryElement(t)
    acc_bc => GetBC()
    if (.not. associated(acc_bc)) return
    if (.not. ActiveBoundaryElement()) return
    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
    acc_load(1,1:n) = GetReal(acc_bc, 'Wave Flux 1', found)
    acc_load(2,1:n) = GetReal(acc_bc, 'Wave Flux 2', found)
    acc_imp(1,1:n) = GetReal(acc_bc, 'Wave Impedance 1', found)
    acc_imp(2,1:n) = GetReal(acc_bc, 'Wave Impedance 2', found)
    acc_rho(1:n) = GetParentMatProp('Density', acc_elem)
    planewave = GetLogical(acc_bc, 'Plane Wave BC', found)
    if (planewave) then
      acc_imp(1,1:n) = GetParentMatProp('Sound Speed', acc_elem)
      acc_imp(2,1:n) = 0.0_dp
    end if
    if (acc_useconv) then
      acc_conv(1,1:n) = GetParentMatProp('Convection Velocity 1', acc_elem, found)
      acc_conv(2,1:n) = GetParentMatProp('Convection Velocity 2', acc_elem, found)
      if (acc_dim == 3) then
        acc_conv(3,1:n) = GetParentMatProp('Convection Velocity 3', acc_elem, found)
      end if
    end if
    call hhz_local_matrix_boundary(n, nd)
    call DefaultUpdateEquations(acc_stiff, acc_force)
  end subroutine hhz_bnd

  ! Harmonic flow/structure interface element body (lines 414-469).
  subroutine hhz_iface(t) bind(c, name="hhz_iface")
    integer(c_int), value :: t
    integer :: n, j, k, l
    logical :: found, flowif, structif
    complex(kind=dp) :: cu, imunit
    acc_elem => GetBoundaryElement(t)
    if (.not. ActiveBoundaryElement()) return
    n = GetElementNOFNodes()
    if (GetElementFamily() == 1) return
    acc_bc => GetBC()
    if (.not. associated(acc_bc)) return
    flowif = ListGetLogical(acc_bc, 'Flow Interface', found)
    structif = ListGetLogical(acc_bc, 'Structure Interface', found)
    if (.not. (flowif .or. structif)) return
    imunit = cmplx(0.0d0, 1.0d0, kind=dp)
    if (flowif) then
      if (any(acc_flowperm(acc_elem % nodeindexes(1:n)) == 0)) then
        call Fatal('HelmholtzSolve', 'Flow solution is not available on boundary')
      end if
      do j = 1, n
        k = acc_flowperm(acc_elem % nodeindexes(j))
        do l = 1, acc_dim
          acc_wall(l,j) = acc_flow((k-1)*acc_flowdofs + 2*l-1) + &
              imunit * acc_flow((k-1)*acc_flowdofs + 2*l)
        end do
      end do
    else if (structif) then
      if (any(acc_dispperm(acc_elem % nodeindexes(1:n)) == 0)) then
        call Fatal('HelmholtzSolve', 'Displacement solution is not available on boundary')
      end if
      do j = 1, n
        k = acc_dispperm(acc_elem % nodeindexes(j))
        do l = 1, acc_dim
          if (acc_dismode == 1) then
            cu = acc_dispeig((k-1)*acc_dispdofs + l)
          else if (acc_dismode == 2) then
            cu = cmplx(acc_disp((k-1)*acc_dispdofs + 2*l-1), &
                acc_disp((k-1)*acc_dispdofs + 2*l), kind=dp)
          else if (acc_dismode == 3) then
            cu = cmplx(acc_disp((k-1)*acc_dispdofs + l), 0.0_dp, kind=dp)
          else
            cu = cmplx(acc_disp((k-1)*acc_dispdofs + l), &
                acc_dispim((k-1)*acc_dispdofs + l), kind=dp)
          end if
          acc_wall(l,j) = cu
        end do
      end do
      acc_wall = imunit * acc_omega * acc_wall
    end if
    call hhz_local_interface_matrix(n)
    call DefaultUpdateEquations(acc_stiff, acc_force)
  end subroutine hhz_iface

  subroutine hhz_assemble_done() bind(c, name="hhz_assemble_done")
    call DefaultFinishAssembly()
    call Info('HelmholtzSolve', 'Assembly done', level=4)
  end subroutine hhz_assemble_done

  subroutine hhz_dirichlet() bind(c, name="hhz_dirichlet")
    call DefaultDirichletBCs()
  end subroutine hhz_dirichlet

  ! Solve + Level=4 timing reports (lines 481-491). The `(s)` values are
  ! wall-clock (normalized by the proof's norm class); the `iter:` index
  ! structure is compared.
  function hhz_solve(iter) result(norm) bind(c, name="hhz_solve")
    integer(c_int), value :: iter
    real(kind=dp) :: norm, st
    acc_at = CPUTime() - acc_at
    st = CPUTime()
    norm = DefaultSolve()
    st = CPUTime() - st
    acc_totat = acc_totat + acc_at
    acc_totst = acc_totst + st
    write(message, '(a,i4,a,F8.2,F8.2)') 'iter: ', iter, ' Assembly: (s)', acc_at, acc_totat
    call Info('HelmholtzSolve', message, level=4)
    write(message, '(a,i4,a,F8.2,F8.2)') 'iter: ', iter, ' Solve:    (s)', st, acc_totst
    call Info('HelmholtzSolve', message, level=4)
    acc_norm = norm
  end function hhz_solve

  function hhz_converged() result(r) bind(c, name="hhz_converged")
    integer(c_int) :: r
    r = merge(1, 0, acc_solver % variable % nonlinconverged == 1)
  end function hhz_converged

  subroutine hhz_finish() bind(c, name="hhz_finish")
    call DefaultFinish()
  end subroutine hhz_finish

  ! -- LocalMatrix internal procedure (exact transcription) --------------
  subroutine hhz_local_matrix(n, nd)
    integer, intent(in) :: n, nd
    real(kind=dp) :: angfreq
    complex(kind=dp) :: a, b, convcoef
    real(kind=dp) :: wavenumber, d, m, l1, l2, rho, sdet
    real(kind=dp) :: diffcoef(3,3), velo(3)
    real(kind=dp) :: x, y, z, metric(3,3), sqrtmetric, symb(3,3,3), dsymb(3,3,3,3)
    integer :: i, j, k, p, q, t, dim, coordsys, nbasis
    logical :: stat
    angfreq = acc_omega
    dim = CoordinateSystemDimension()
    coordsys = CurrentCoordinateSystem()
    metric = 0.0_dp
    metric(1,1) = 1.0_dp
    metric(2,2) = 1.0_dp
    metric(3,3) = 1.0_dp
    acc_stiff = 0.0_dp
    acc_force = 0.0_dp
    if (acc_useprecshift) acc_prec = 0.0_dp
    diffcoef = 0.0_dp
    call GetElementNodes(acc_nodes_lm)
    if (acc_bubbles) then
      acc_ip_lm = GaussPoints(acc_elem, acc_elem % type % gausspoints2)
      nbasis = 2*n
    else
      nbasis = nd
      acc_ip_lm = GaussPoints(acc_elem)
    end if
    do t = 1, acc_ip_lm % n
      stat = ElementInfo(acc_elem, acc_nodes_lm, acc_ip_lm % u(t), &
          acc_ip_lm % v(t), acc_ip_lm % w(t), sdet, &
          acc_basis, acc_dbasis, bubbles=acc_bubbles)
      sdet = acc_ip_lm % s(t) * sdet
      if (coordsys /= Cartesian) then
        x = sum(acc_nodes_lm % x(1:n) * acc_basis(1:n))
        y = sum(acc_nodes_lm % y(1:n) * acc_basis(1:n))
        z = sum(acc_nodes_lm % z(1:n) * acc_basis(1:n))
        call CoordinateSystemInfo(metric, sqrtmetric, symb, dsymb, x, y, z)
        sdet = sdet * sqrtmetric
      end if
      wavenumber = angfreq / sum(acc_snd(1:n) * acc_basis(1:n))
      d = wavenumber * sum(acc_damp(1:n) * acc_basis(1:n))
      m = -wavenumber**2
      l1 = sum(acc_load(1,1:n) * acc_basis(1:n))
      l2 = sum(acc_load(2,1:n) * acc_basis(1:n))
      rho = sum(acc_rho(1:n) * acc_basis(1:n))
      do i = 1, dim
        diffcoef(i,i) = 1.0_dp
      end do
      if (acc_useconv) then
        convcoef = 2.0_dp * sqrt((-1.0_dp, 0.0_dp)) * wavenumber
        velo(1) = sum(acc_conv(1,1:n) * acc_basis(1:n))
        velo(2) = sum(acc_conv(2,1:n) * acc_basis(1:n))
        velo(3) = sum(acc_conv(3,1:n) * acc_basis(1:n))
        velo = velo / sum(acc_snd(1:n) * acc_basis(1:n))
        do i = 1, dim
          do j = 1, dim
            diffcoef(i,j) = diffcoef(i,j) - velo(i)*velo(j)
          end do
        end do
      end if
      do p = 1, nbasis
        do q = 1, nbasis
          a = cmplx(m, d, kind=dp) * acc_basis(q) * acc_basis(p)
          if (acc_useprecshift) then
            b = -m * acc_shift * acc_basis(q) * acc_basis(p)
            acc_prec(p,q) = acc_prec(p,q) + sdet * b / rho
          end if
          do i = 1, dim
            if (acc_useconv) then
              a = a + convcoef * velo(i) * acc_dbasis(q,i) * acc_basis(p)
            end if
            do j = 1, dim
              do k = 1, dim
                a = a + metric(i,j) * diffcoef(i,k) * acc_dbasis(q,k) * acc_dbasis(p,j)
              end do
            end do
          end do
          a = a / rho
          acc_stiff(p,q) = acc_stiff(p,q) + sdet * a
        end do
        b = acc_basis(p) * cmplx(l1, l2, kind=dp)
        b = b / rho
        acc_force(p) = acc_force(p) + sdet * b
      end do
    end do
    if (acc_useprecshift) then
      acc_prec = acc_stiff + acc_prec
      if (acc_bubbles) then
        call CondensateP(n, n, acc_prec)
        call DefaultUpdatePrec(acc_prec(1:n,1:n))
      else
        call DefaultUpdatePrec(acc_prec(1:nbasis,1:nbasis))
      end if
    end if
    if (acc_bubbles) then
      call CondensateP(n, n, acc_stiff, acc_force)
    end if
  end subroutine hhz_local_matrix

  ! -- LocalMatrixBoundary internal procedure (exact transcription) -----
  subroutine hhz_local_matrix_boundary(n, nd)
    integer, intent(in) :: n, nd
    complex(kind=dp) :: a, b, admittance
    real(kind=dp) :: imp1, imp2, l1, l2, sdet, rho
    real(kind=dp) :: x, y, z, normal(3), velo(3), normvelo, tangvelo(3)
    real(kind=dp) :: metric(3,3), sqrtmetric, symb(3,3,3), dsymb(3,3,3,3)
    integer :: i, p, q, t, dim, coordsys
    logical :: stat
    dim = CoordinateSystemDimension()
    coordsys = CurrentCoordinateSystem()
    acc_stiff = 0.0d0
    acc_force = 0.0d0
    call GetElementNodes(acc_nodes_lb)
    acc_ip_lb = GaussPoints(acc_elem)
    do t = 1, acc_ip_lb % n
      stat = ElementInfo(acc_elem, acc_nodes_lb, acc_ip_lb % u(t), &
          acc_ip_lb % v(t), acc_ip_lb % w(t), sdet, &
          acc_basis, acc_dbasis)
      sdet = acc_ip_lb % s(t) * sdet
      if (coordsys /= Cartesian) then
        x = sum(acc_nodes_lb % x(1:n) * acc_basis(1:n))
        y = sum(acc_nodes_lb % y(1:n) * acc_basis(1:n))
        z = sum(acc_nodes_lb % z(1:n) * acc_basis(1:n))
        call CoordinateSystemInfo(metric, sqrtmetric, symb, dsymb, x, y, z)
        sdet = sdet * sqrtmetric
      end if
      normal = Normalvector(acc_elem, acc_nodes_lb, acc_ip_lb % u(t), acc_ip_lb % v(t), .true.)
      imp1 = sum(acc_imp(1,1:n) * acc_basis(1:n))
      imp2 = sum(acc_imp(2,1:n) * acc_basis(1:n))
      if (abs(imp1) < AEPS .and. abs(imp2) < AEPS) then
        admittance = cmplx(0.0d0, 0.0d0, kind=dp)
      else
        admittance = cmplx(0.0d0, 1.0d0, kind=dp) * acc_omega / cmplx(imp1, imp2, kind=dp)
      end if
      rho = sum(acc_rho(1:n) * acc_basis(1:n))
      if (acc_useconv) then
        velo(1) = sum(acc_conv(1,1:n) * acc_basis(1:n))
        velo(2) = sum(acc_conv(2,1:n) * acc_basis(1:n))
        velo(3) = sum(acc_conv(3,1:n) * acc_basis(1:n))
        velo = velo / sum(acc_snd(1:n) * acc_basis(1:n))
        normvelo = sum(normal(1:dim) * velo(1:dim))
        tangvelo = velo - normal * normvelo
      end if
      l1 = sum(acc_load(1,1:n) * acc_basis(1:n))
      l2 = sum(acc_load(2,1:n) * acc_basis(1:n))
      do p = 1, nd
        do q = 1, nd
          a = admittance * acc_basis(q) * acc_basis(p)
          if (acc_useconv) then
            a = a - normvelo * admittance * acc_basis(q) * acc_basis(p)
            a = a + sum(tangvelo(1:dim) * acc_dbasis(q,1:dim)) * acc_basis(p)
          end if
          a = a / rho
          acc_stiff(p,q) = acc_stiff(p,q) + sdet * a
        end do
        b = acc_basis(p) * cmplx(l1, l2, kind=dp)
        b = b / rho
        acc_force(p) = acc_force(p) + sdet * b
      end do
    end do
    if (acc_useprecshift) call DefaultUpdatePrec(acc_stiff(1:nd,1:nd))
  end subroutine hhz_local_matrix_boundary

  ! -- LocalInterfaceMatrix internal procedure (exact transcription) ----
  subroutine hhz_local_interface_matrix(n)
    integer, intent(in) :: n
    complex(kind=dp) :: normvelo
    complex(kind=dp) :: imunit
    real(kind=dp) :: sdet, x, y, z, normal(3)
    real(kind=dp) :: metric(3,3), sqrtmetric, symb(3,3,3), dsymb(3,3,3,3)
    integer :: i, p, t, dim, coordsys
    logical :: stat
    dim = CoordinateSystemDimension()
    coordsys = CurrentCoordinateSystem()
    acc_stiff = 0.0d0
    acc_force = 0.0d0
    call GetElementNodes(acc_nodes_li)
    acc_ip_li = GaussPoints(acc_elem)
    do t = 1, acc_ip_li % n
      stat = ElementInfo(acc_elem, acc_nodes_li, acc_ip_li % u(t), &
          acc_ip_li % v(t), acc_ip_li % w(t), sdet, &
          acc_basis, acc_dbasis)
      sdet = acc_ip_li % s(t) * sdet
      if (coordsys /= Cartesian) then
        x = sum(acc_nodes_li % x(1:n) * acc_basis(1:n))
        y = sum(acc_nodes_li % y(1:n) * acc_basis(1:n))
        z = sum(acc_nodes_li % z(1:n) * acc_basis(1:n))
        call CoordinateSystemInfo(metric, sqrtmetric, symb, dsymb, x, y, z)
        sdet = sdet * sqrtmetric
      end if
      normal = Normalvector(acc_elem, acc_nodes_li, acc_ip_li % u(t), acc_ip_li % v(t), .true.)
      normvelo = 0.0_dp
      do i = 1, dim
        normvelo = normvelo + normal(i) * sum(acc_wall(i,1:n) * acc_basis(1:n))
      end do
      imunit = cmplx(0.0d0, 1.0d0, kind=dp)
      do p = 1, n
        acc_force(p) = acc_force(p) - sdet * acc_basis(p) * &
            imunit * acc_omega * normvelo
      end do
    end do
  end subroutine hhz_local_interface_matrix

  ! -- HelmholtzSolver_init body ----------------------------------------
  subroutine hhz_init_begin(h_solver) bind(c, name="hhz_init_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine hhz_init_begin

  subroutine hhz_init_body() bind(c, name="hhz_init_body")
    type(valuelist_t), pointer :: params
    params => GetSolverParams()
    call ListAddNewLogical(params, 'Linear System Complex', .true.)
    if (ListCheckPrefix(params, 'Linear System Preconditioning Damp Coefficient')) then
      call ListAddNewLogical(params, 'Allocate Preconditioning Matrix', .true.)
    end if
  end subroutine hhz_init_body

end module helmholtzsolve_acc
