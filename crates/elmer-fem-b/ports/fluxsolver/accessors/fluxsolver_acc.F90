! T022-A FluxSolver port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/FluxSolver.F90, compiled against the built
! .mod files of the target configure (A2 shipped mechanism, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments. All Fortran-object
! traffic — TYPE(...) handles, OPTIONAL/POINTER dummies, assumed-shape
! descriptors, hidden lengths, list-directed timing reports — stays inside
! Fortran, so the Rust side is pure C calls. The wrapper module mirrors
! FluxSolver's own locals (current element, Nodes/IP SAVE state,
! entry-level flags, FieldTable pointers, ForceVector) with identical
! lifetimes. All floating-point work (the Galerkin assembly, MATMUL/SUM
! reductions, SQRT magnitude post-processing) and the RealTime-based
! reports stay in Fortran (same compiler), so bit-identity holds by
! construction wherever the call sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module fluxsolver_acc
  use iso_c_binding
  use coordinatesystems
  use defutils
  implicit none

  type fieldtab_t
    real(kind=dp), pointer :: values(:) => null()
  end type fieldtab_t

  type(solver_t), pointer, save :: acc_solver => null()
  type(valuelist_t), pointer, save :: acc_params => null()
  type(valuelist_t), pointer, save :: acc_mat => null()
  type(element_t), pointer, save :: acc_elem => null()
  ! BulkAssembly context (mirrors its SAVE Conductivity, Nodes; IntegStuff
  ! is a fresh local there, kept module-side with identical per-element
  ! lifetime).
  real(kind=dp), pointer, save :: acc_cond(:,:,:) => null()
  type(nodes_t), target, save :: acc_nodes
  type(gaussintegrationpoints_t), save :: acc_ip
  ! Entry-level state (mirrors the solver body's locals).
  type(fieldtab_t), save :: acc_fields(8)
  real(kind=dp), allocatable, target, save :: acc_forcevec(:,:)
  real(kind=dp), pointer, save :: acc_saverhs(:) => null()
  ! LocalJumps context (mirrors its SAVE FaceNodes, P1Nodes, P2Nodes).
  type(nodes_t), save :: acc_facenodes, acc_p1nodes, acc_p2nodes
  integer(c_int), save :: acc_dim = 0
  integer(c_int), save :: acc_dofs = 0
  integer(c_int), save :: acc_firstmag = 0
  logical, save :: acc_calcflux = .false.
  logical, save :: acc_calcfluxmag = .false.
  logical, save :: acc_calcgrad = .false.
  logical, save :: acc_calcgradmag = .false.
  logical, save :: acc_usepot = .false.
  logical, save :: acc_autonorm = .false.
  logical, save :: acc_csym = .false.
  logical, save :: acc_constbulk = .false.
  logical, save :: acc_inuse = .false.
  character(len=MAX_NAME_LEN), save :: acc_varname = ''
  character(len=MAX_NAME_LEN), save :: acc_potname = ''
  character(len=MAX_NAME_LEN), save :: acc_condname = ''
  real(kind=dp), save :: acc_at0 = 0.0_dp
  real(kind=dp), save :: acc_at1 = 0.0_dp

contains

  ! -- entry state ----------------------------------------------------
  subroutine fls_begin(h_solver) bind(c, name="fls_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine fls_begin

  subroutine fls_info_start() bind(c, name="fls_info_start")
    call Info('FluxSolver', '-------------------------------------', level=4)
    call Info('FluxSolver', 'Computing the flux and/or gradient', level=4)
    call Info('FluxSolver', '-------------------------------------', level=4)
  end subroutine fls_info_start

  function fls_guard() result(code) bind(c, name="fls_guard")
    integer(c_int) :: code
    if (.not. associated(acc_solver % matrix)) then
      code = 1
      return
    end if
    if (count(acc_solver % variable % perm > 0) <= 0) then
      code = 2
      return
    end if
    code = 0
  end function fls_guard

  subroutine fls_warn_nofield() bind(c, name="fls_warn_nofield")
    call Warn('FluxSolver', 'No field computation requested, exiting...')
  end subroutine fls_warn_nofield

  subroutine fls_fatal_alloc() bind(c, name="fls_fatal_alloc")
    call Fatal('FluxSolver', 'Memory allocation error.')
  end subroutine fls_fatal_alloc

  ! -- sweep setup (entry lines up to DefaultInitialize) ----------------
  function fls_sweep_begin() result(dofs) bind(c, name="fls_sweep_begin")
    integer(c_int) :: dofs
    type(valuelist_t), pointer :: sp
    character(len=MAX_NAME_LEN) :: varname, condname, potname
    integer :: i, dim
    logical :: gotit, usepot, autonorm, constbulk
    type(variable_t), pointer :: fluxsol
    dim = CoordinateSystemDimension()
    acc_dim = dim
    sp => GetSolverParams()
    acc_params => sp
    acc_calcgradmag = GetLogical(sp, 'Calculate Grad Magnitude', gotit)
    acc_calcgrad = GetLogical(sp, 'Calculate Grad', gotit)
    acc_calcfluxmag = GetLogical(sp, 'Calculate Flux Magnitude', gotit)
    acc_calcflux = GetLogical(sp, 'Calculate Flux', gotit)
    dofs = 0
    if (acc_calcflux) dofs = dofs + dim
    if (acc_calcfluxmag) dofs = dofs + 1
    if (acc_calcgrad) dofs = dofs + dim
    if (acc_calcgradmag) dofs = dofs + 1
    acc_dofs = dofs
    if (dofs == 0) return
    autonorm = .not. ListCheckPresent(sp, 'Skip Compute Nonlinear Change')
    if (autonorm) then
      call ListAddLogical(sp, 'Skip Compute Nonlinear Change', .true.)
    end if
    acc_autonorm = autonorm
    varname = GetString(sp, 'Flux Variable', gotit)
    usepot = .false.
    if (.not. gotit) varname = GetString(sp, 'Target Variable', gotit)
    if (.not. gotit) then
      potname = GetString(sp, 'Target Expression', usepot)
      if (usepot) then
        varname = potname
      else
        call Warn('FluxSolver', '> Target Variable < not given, using Temperature')
        varname = 'Temperature'
      end if
    end if
    acc_varname = varname
    acc_potname = potname
    acc_usepot = usepot
    if (acc_calcflux .or. acc_calcfluxmag) then
      condname = ListGetString(sp, 'Flux Coefficient', gotit)
      if (.not. gotit) then
        call Warn('FluxSolver', '> Flux Coefficient < not given, using Heat Conductivity')
        condname = 'Heat Conductivity'
      end if
      acc_condname = condname
    end if
    i = 0
    if (acc_calcflux) then
      fluxsol => VariableGet(acc_solver % mesh % variables, trim(varname)//' Flux 1')
      acc_fields(1) % values => fluxsol % values
      fluxsol => VariableGet(acc_solver % mesh % variables, trim(varname)//' Flux 2')
      acc_fields(2) % values => fluxsol % values
      if (dim == 3) then
        fluxsol => VariableGet(acc_solver % mesh % variables, trim(varname)//' Flux 3')
        acc_fields(3) % values => fluxsol % values
      end if
      i = i + dim
    end if
    if (acc_calcgrad) then
      fluxsol => VariableGet(acc_solver % mesh % variables, trim(varname)//' Grad 1')
      acc_fields(i+1) % values => fluxsol % values
      fluxsol => VariableGet(acc_solver % mesh % variables, trim(varname)//' Grad 2')
      acc_fields(i+2) % values => fluxsol % values
      if (dim == 3) then
        fluxsol => VariableGet(acc_solver % mesh % variables, trim(varname)//' Grad 3')
        acc_fields(i+3) % values => fluxsol % values
      end if
      i = i + dim
    end if
    acc_firstmag = i + 1
    if (acc_calcfluxmag) then
      fluxsol => VariableGet(acc_solver % mesh % variables, trim(varname)//' Flux_mag')
      i = i + 1
      acc_fields(i) % values => fluxsol % values
    end if
    if (acc_calcgradmag) then
      fluxsol => VariableGet(acc_solver % mesh % variables, trim(varname)//' Grad_mag')
      i = i + 1
      acc_fields(i) % values => fluxsol % values
    end if
    acc_csym = CurrentCoordinateSystem() == AxisSymmetric .or. &
        CurrentCoordinateSystem() == CylindricSymmetric
    acc_at0 = RealTime()
    constbulk = GetLogical(sp, 'Constant Bulk Matrix', gotit)
    acc_constbulk = constbulk
    acc_inuse = constbulk .and. associated(acc_solver % matrix % bulkvalues)
    call DefaultInitialize(acc_solver, acc_inuse)
  end function fls_sweep_begin

  function fls_alloc_n() result(n) bind(c, name="fls_alloc_n")
    integer(c_int) :: n
    n = 2 * max(acc_solver % mesh % maxelementdofs, &
        acc_solver % mesh % maxelementnodes)
  end function fls_alloc_n

  function fls_nrhs() result(n) bind(c, name="fls_nrhs")
    integer(c_int) :: n
    n = size(acc_solver % matrix % rhs)
  end function fls_nrhs

  function fls_nof_active() result(n) bind(c, name="fls_nof_active")
    integer(c_int) :: n
    n = acc_solver % numberofactiveelements
  end function fls_nof_active

  ! -- ForceVector setup (entry ForceVector alloc + SaveRHS bind) -------
  function fls_fields_begin() result(code) bind(c, name="fls_fields_begin")
    integer(c_int) :: code
    integer :: istat
    allocate(acc_forcevec(fls_nrhs(), acc_dofs), stat=istat)
    if (istat /= 0) then
      code = 1
      return
    end if
    code = 0
    acc_forcevec = 0.0_dp
    acc_saverhs => acc_solver % matrix % rhs
  end function fls_fields_begin

  ! -- BulkAssembly element body (exact transcription) ------------------
  subroutine fls_elem(elem, nmax, ndofs, stiff, force, locpot, basis, dbasis) &
      bind(c, name="fls_elem")
    integer(c_int), value :: elem, nmax, ndofs
    real(c_double) :: stiff(nmax, nmax), force(ndofs, nmax)
    real(c_double) :: locpot(nmax), basis(nmax), dbasis(nmax, 3)
    integer :: t, i, j, k, p, q, n, nd, rank
    real(kind=dp) :: weight, coeff, detj, c(3, 3)
    real(kind=dp) :: fluxatip(3), gradatip(3)
    logical :: found
    acc_elem => GetActiveElement(elem)
    acc_mat => GetMaterial()
    call GetElementNodes(acc_nodes)
    nd = GetElementNOFDOFs()
    n = GetElementNOFNodes()
    acc_ip = GaussPoints(acc_elem)
    stiff = 0.0_c_double
    force = 0.0_c_double
    if (acc_calcflux .or. acc_calcfluxmag) then
      call GetRealArray(acc_mat, acc_cond, acc_condname, found)
      rank = 0
      if (found) rank = fls_rank(acc_cond)
      c = 0.0_dp
      fluxatip = 0.0_dp
      do i = 1, acc_dim
        c(i, i) = 1.0_dp
      end do
    end if
    if (acc_usepot) then
      locpot(1:n) = GetReal(acc_mat, acc_potname)
    else
      call GetScalarLocalSolution(locpot(1:n), acc_varname)
    end if
    do t = 1, acc_ip % n
      found = ElementInfo(acc_elem, acc_nodes, acc_ip % u(t), &
          acc_ip % v(t), acc_ip % w(t), detj, basis(1:n), dbasis(1:nd, 1:3))
      weight = acc_ip % s(t) * detj
      if (acc_csym) weight = weight * sum(basis(1:n) * acc_nodes % x(1:n))
      if (.not. acc_inuse) then
        do p = 1, nd
          do q = 1, nd
            stiff(p, q) = stiff(p, q) + weight * basis(q) * basis(p)
          end do
        end do
      end if
      gradatip(1:acc_dim) = matmul(locpot(1:nd), dbasis(1:nd, 1:acc_dim))
      if (acc_calcflux .or. acc_calcfluxmag) then
        select case (rank)
        case (0)
        case (1)
          do i = 1, acc_dim
            c(i, i) = sum(basis(1:n) * acc_cond(1, 1, 1:n))
          end do
        case (2)
          do i = 1, acc_dim
            c(i, i) = sum(basis(1:n) * acc_cond(i, 1, 1:n))
          end do
        case default
          do i = 1, acc_dim
            do j = 1, acc_dim
              c(i, j) = sum(basis(1:n) * acc_cond(i, j, 1:n))
            end do
          end do
        end select
        do i = 1, acc_dim
          fluxatip(i) = sum(c(i, 1:acc_dim) * gradatip(1:acc_dim))
        end do
      end if
      k = 0
      if (acc_calcflux) then
        do i = 1, acc_dim
          coeff = weight * fluxatip(i)
          force(i, 1:nd) = force(i, 1:nd) - coeff * basis(1:nd)
        end do
        k = k + acc_dim
      end if
      if (acc_calcgrad) then
        do i = 1, acc_dim
          coeff = weight * gradatip(i)
          force(k + i, 1:nd) = force(k + i, 1:nd) + coeff * basis(1:nd)
        end do
        k = k + acc_dim
      end if
      if (acc_calcfluxmag) then
        coeff = weight * sqrt(sum(fluxatip**2))
        k = k + 1
        force(k, 1:nd) = force(k, 1:nd) + coeff * basis(1:nd)
      end if
      if (acc_calcgradmag) then
        coeff = weight * sqrt(sum(gradatip**2))
        k = k + 1
        force(k, 1:nd) = force(k, 1:nd) + coeff * basis(1:nd)
      end if
    end do
    if (.not. acc_inuse) then
      acc_solver % matrix % rhs => acc_saverhs
      call DefaultUpdateEquations(stiff, force(1, 1:nd))
    end if
    do i = 1, acc_dofs
      acc_solver % matrix % rhs => acc_forcevec(:, i)
      call DefaultUpdateForce(force(i, 1:nd))
    end do
  end subroutine fls_elem

  ! -- face terms (AddLocalFaceTerms transcription) ---------------------
  subroutine fls_faces(nmax, ndofs, stiff, force) bind(c, name="fls_faces")
    integer(c_int), value :: nmax, ndofs
    real(c_double) :: stiff(nmax, nmax), force(ndofs, nmax)
    type(element_t), pointer :: p1, p2, face, faces(:)
    integer :: t, n, n1, n2, numberoffaces, dim
    logical :: found
    if (.not. GetLogical(acc_params, 'Discontinuous Galerkin', found)) return
    if (.not. GetLogical(acc_params, 'Average Within Materials', found)) return
    force = 0.0_c_double
    dim = CoordinateSystemDimension()
    if (dim == 2) then
      faces => acc_solver % mesh % edges
      numberoffaces = acc_solver % mesh % numberofedges
    else
      faces => acc_solver % mesh % faces
      numberoffaces = acc_solver % mesh % numberoffaces
    end if
    do t = 1, numberoffaces
      face => faces(t)
      if (.not. ActiveBoundaryElement(face)) cycle
      p1 => face % boundaryinfo % left
      p2 => face % boundaryinfo % right
      if (associated(p2) .and. associated(p1)) then
        if (.not. associated(GetMaterial(p1), GetMaterial(p2))) cycle
        n = GetElementNOFNodes(face)
        n1 = GetElementNOFNodes(p1)
        n2 = GetElementNOFNodes(p2)
        call fls_localjumps(stiff, face, n, p1, n1, p2, n2)
        call DefaultUpdateEquations(stiff, force(1, :), face)
      end if
    end do
  end subroutine fls_faces

  subroutine fls_localjumps(stiff, face, n, p1, n1, p2, n2)
    real(kind=dp) :: stiff(:, :)
    integer :: n, n1, n2
    type(element_t) :: face, p1, p2
    real(kind=dp) :: facebasis(n), p1basis(n1), p2basis(n2)
    real(kind=dp) :: jump(n1 + n2), detj, u, v, w, s
    logical :: stat
    integer :: i, j, p, q, t, nface, nparent
    type(gaussintegrationpoints_t) :: integstuff
    nface = n
    nparent = n1 + n2
    stiff = 0._dp
    call GetElementNodes(acc_facenodes, face)
    call GetElementNodes(acc_p1nodes, p1)
    call GetElementNodes(acc_p2nodes, p2)
    integstuff = GaussPoints(face)
    do t = 1, integstuff % n
      u = integstuff % u(t)
      v = integstuff % v(t)
      w = integstuff % w(t)
      s = integstuff % s(t)
      stat = ElementInfo(face, acc_facenodes, u, v, w, detj, facebasis)
      s = s * detj
      call GetParentUVW(face, n, p1, n1, u, v, w, facebasis)
      stat = ElementInfo(p1, acc_p1nodes, u, v, w, detj, p1basis)
      call GetParentUVW(face, n, p2, n2, u, v, w, facebasis)
      stat = ElementInfo(p2, acc_p2nodes, u, v, w, detj, p2basis)
      jump(1:n1) = p1basis(1:n1)
      jump(n1 + 1:n1 + n2) = -p2basis(1:n2)
      do p = 1, n1 + n2
        do q = 1, n1 + n2
          stiff(p, q) = stiff(p, q) + s * jump(q) * jump(p)
        end do
      end do
    end do
  end subroutine fls_localjumps

  ! -- assembly tail (bulk finish + assembly time report) ----------------
  subroutine fls_finish() bind(c, name="fls_finish")
    real(kind=dp) :: at1
    character(len=MAX_NAME_LEN) :: message
    if (acc_constbulk) then
      if (.not. acc_inuse) then
        call Info('FluxSolver', 'Saving the system matrix', level=6)
        call CopyBulkMatrix(acc_solver % matrix, bulkRHS=.false.)
      end if
      call DefaultFinishBulkAssembly(bulkUpdate=.false.)
    else
      call DefaultFinishBulkAssembly()
    end if
    call DefaultFinishAssembly()
    at1 = RealTime()
    acc_at1 = at1
    write (message, *) 'Assembly Time: ', at1 - acc_at0
    call Info('FluxSolver', message, level=5)
  end subroutine fls_finish

  ! -- component solve loop + magnitude post-processing + reports -------
  function fls_solve() result(totnorm) bind(c, name="fls_solve")
    real(kind=dp) :: totnorm
    integer :: i, j, k, dim
    real(kind=dp) :: unorm, val
    logical :: gotit, enforce
    type(variable_t), pointer :: fluxsol
    character(len=MAX_NAME_LEN) :: message
    real(kind=dp) :: at2
    dim = acc_dim
    enforce = GetLogical(acc_params, 'Enforce Positive Magnitude', gotit)
    totnorm = 0.0_dp
    do i = 1, acc_dofs
      acc_solver % matrix % rhs => acc_forcevec(:, i)
      if (i == acc_dofs .and. acc_autonorm) then
        call ListRemove(acc_params, 'Skip Compute Nonlinear Change')
      end if
      unorm = DefaultSolve()
      acc_fields(i) % values = acc_solver % variable % values
      if (enforce .and. i >= acc_firstmag) then
        do j = 1, size(acc_fields(i) % values)
          acc_fields(i) % values(j) = max(acc_fields(i) % values(j), 0.0_dp)
        end do
      end if
    end do
    deallocate (acc_forcevec)
    acc_solver % matrix % rhs => acc_saverhs
    totnorm = unorm
    acc_solver % variable % norm = totnorm
    fluxsol => VariableGet(acc_solver % mesh % variables, &
        trim(acc_varname)//' Flux_abs')
    if (associated(fluxsol)) then
      do j = 1, size(fluxsol % values)
        val = acc_fields(1) % values(j)**2 + acc_fields(2) % values(j)**2
        if (dim == 3) val = val + acc_fields(3) % values(j)**2
        fluxsol % values(j) = sqrt(val)
      end do
    end if
    fluxsol => VariableGet(acc_solver % mesh % variables, &
        trim(acc_varname)//' Grad_abs')
    if (associated(fluxsol)) then
      k = 0
      if (acc_calcflux) k = dim
      do j = 1, size(fluxsol % values)
        val = acc_fields(k + 1) % values(j)**2 + acc_fields(k + 2) % values(j)**2
        if (dim == 3) val = val + acc_fields(k + 3) % values(j)**2
        fluxsol % values(j) = sqrt(val)
      end do
    end if
    at2 = RealTime()
    write (message, *) 'Solution Time: ', at2 - acc_at1
    call Info('FluxSolver', message, level=5)
    write (message, *) 'Result Norm: ', totnorm
    call Info('FluxSolver', message, level=4)
    call Info('FluxSolver', 'All done', level=6)
    call Info('FluxSolver', '-------------------------------------', level=6)
  end function fls_solve

  function fls_rank(tensor) result(rank)
    real(kind=dp) :: tensor(:, :, :)
    integer :: rank
    if (size(tensor, 1) == 1) then
      rank = 1
    else if (size(tensor, 2) == 1) then
      rank = 2
    else
      rank = 3
    end if
  end function fls_rank

  ! -- FluxSolver_Init body ---------------------------------------------
  subroutine fls_init_begin(h_solver) bind(c, name="fls_init_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine fls_init_begin

  subroutine fls_init_body() bind(c, name="fls_init_body")
    type(valuelist_t), pointer :: solverpars
    integer :: dim
    character(len=MAX_NAME_LEN) :: eqname, varname, fluxname, gradname
    logical :: gotit, calcfluxmag, calcflux, calcgrad, calcgradmag
    solverpars => GetSolverParams()
    dim = CoordinateSystemDimension()
    if (dim < 2 .or. dim > 3) then
      call Fatal('FluxSolver_init', 'Flux computation makes sense only in 2D and 3D')
    end if
    calcgradmag = GetLogical(solverpars, 'Calculate Grad Magnitude', gotit)
    calcgrad = GetLogical(solverpars, 'Calculate Grad', gotit)
    calcfluxmag = GetLogical(solverpars, 'Calculate Flux Magnitude', gotit)
    calcflux = GetLogical(solverpars, 'Calculate Flux', gotit)
    if (.not. (calcflux .or. calcfluxmag .or. calcgrad .or. calcgradmag)) then
      if (.not. gotit) then
        call Warn('FluxSolver_init', 'No field computation requested, setting > Calculate Flux = True <')
      else
        call Warn('FluxSolver_init', 'No field computation requested!')
        return
      end if
      calcflux = .true.
      call ListAddLogical(solverpars, 'Calculate Flux', calcflux)
    end if
    varname = GetString(solverpars, 'Flux Variable', gotit)
    if (.not. gotit) varname = GetString(solverpars, 'Target Variable', gotit)
    if (.not. gotit) then
      if (.not. gotit) then
        varname = GetString(solverpars, 'Target Expression', gotit)
      end if
      if (.not. gotit) then
        call Warn('FluxSolver', '> Target Variable < not given, using Temperature')
        varname = 'Temperature'
      end if
    end if
    if (.not. ListCheckPresent(solverpars, 'Variable')) then
      eqname = ListGetString(solverpars, 'Equation')
      call ListAddString(solverpars, 'Variable', '-nooutput '//trim(eqname)//'_temp')
    end if
    if (calcflux) then
      fluxname = trim(varname)//' Flux'
      call Info('FluxSolver_init', 'Saving flux to: '//trim(fluxname), level=5)
      if (dim == 2) then
        call ListAddString(solverpars, &
            NextFreeKeyword('Exported Variable', solverpars), &
            trim(fluxname)//'['//trim(fluxname)//':2]')
      else if (dim == 3) then
        call ListAddString(solverpars, &
            NextFreeKeyword('Exported Variable', solverpars), &
            trim(fluxname)//'['//trim(fluxname)//':3]')
      end if
      if (GetLogical(solverpars, 'Calculate Flux Abs', gotit)) then
        fluxname = trim(varname)//' Flux_abs'
        call Info('FluxSolver_init', 'Saving flux abs to: '//fluxname)
        call ListAddString(solverpars, &
            NextFreeKeyword('Exported Variable', solverpars), trim(fluxname))
      end if
    end if
    if (calcfluxmag) then
      fluxname = trim(varname)//' Flux_mag'
      call Info('FluxSolver_init', 'Saving flux magnitude to: '//fluxname)
      call ListAddString(solverpars, &
          NextFreeKeyword('Exported Variable', solverpars), trim(fluxname))
    end if
    if (calcgrad) then
      gradname = trim(varname)//' Grad'
      call Info('FluxSolver_init', 'Saving gradient to: '//gradname)
      if (dim == 2) then
        call ListAddString(solverpars, &
            NextFreeKeyword('Exported Variable', solverpars), &
            trim(gradname)//'['//trim(gradname)//':2]')
      else if (dim == 3) then
        call ListAddString(solverpars, &
            NextFreeKeyword('Exported Variable', solverpars), &
            trim(gradname)//'['//trim(gradname)//':3]')
      end if
      if (GetLogical(solverpars, 'Calculate Grad Abs', gotit)) then
        gradname = trim(varname)//' Grad_abs'
        call Info('FluxSolver_init', 'Saving gradient abs to: '//gradname)
        call ListAddString(solverpars, &
            NextFreeKeyword('Exported Variable', solverpars), trim(gradname))
      end if
    end if
    if (calcgradmag) then
      gradname = trim(varname)//' Grad_mag'
      call Info('FluxSolver_init', 'Saving gradient magnitude to: '//gradname)
      call ListAddString(solverpars, &
          NextFreeKeyword('Exported Variable', solverpars), trim(gradname))
    end if
    call ListAddInteger(solverpars, 'Time derivative order', 0)
    call ListAddNewString(solverpars, 'Linear System Solver', 'Iterative')
    call ListAddNewString(solverpars, 'Linear System Iterative Method', 'cg')
    call ListAddNewString(solverpars, 'Linear System Preconditioning', 'ILU0')
    call ListAddNewInteger(solverpars, 'Linear System Max Iterations', 500)
    call ListAddNewInteger(solverpars, 'Linear System Residual Output', 10)
    call ListAddNewConstReal(solverpars, 'Linear System Convergence Tolerance', 1.0e-10_dp)
  end subroutine fls_init_body

end module fluxsolver_acc
