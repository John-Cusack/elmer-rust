! T022-A FreeSurfaceSolver port: statement-level BIND(C) wrappers for
! the Elmer calls in fem/src/modules/FreeSurfaceSolver.F90, compiled
! against the built .mod files of the target configure (A2 shipped
! mechanism, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, out/in-out scalars, explicit-shape
! REAL(C_DOUBLE) arrays). All Fortran-object traffic — TYPE(...) handles,
! OPTIONAL/POINTER dummies, assumed-shape descriptors, hidden lengths,
! formatted Message writes — stays inside Fortran, so the Rust side is
! pure C calls. The wrapper module mirrors the original's SAVE state
! (VariableName, PrevSize, AllocationsDone, ReAllocate, the STIFF/MASS/...
! work arrays, SubstantialSurface, NormalFlux, Bubbles, RotateFS state,
! limiter arrays, OldValues/OldRHS, FirstTime/t0 of the mesh-update
! function) with identical lifetimes, and every formatted WRITE + Info /
! Warn / Fatal call site is transcribed statement for statement so
! bit-identity holds by construction wherever the call sequence matches.
!
! The four FreeSurfaceToMeshUpdate functions share one core (fss_mu) so
! the SAVE FirstTime/t0 instance is single, exactly as in the original
! where 1/2/3 call the common internal function.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module freesurfacesolver_acc
  use iso_c_binding
  use defutils
  use elementdescription
  use differentials
  use materialmodels
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()

  ! -- FreeSurfaceSolver_RotInit state -----------------------------------
  character(len=MAX_NAME_LEN), save :: rot_sname = 'FreeSurfaceSolver_RotInit'
  character(len=MAX_NAME_LEN), save :: rot_varname = ''
  real(c_double), save :: rot_orient(3) = 0.0_c_double
  real(c_double), save :: rot_rotmat(3,3) = 0.0_c_double
  character(len=MAX_NAME_LEN), save :: rot_ovname = ''
  real(c_double), pointer, save :: rot_ovp(:) => null()
  real(c_double), pointer, save :: rot_ref(:) => null()
  integer, pointer, save :: rot_rperm(:) => null()
  real(c_double), pointer, save :: rot_fs(:) => null()
  integer, pointer, save :: rot_fperm(:) => null()

  ! -- FreeSurfaceToMeshUpdate shared SAVE -------------------------------
  logical, save :: mu_first
  real(c_double), save :: mu_t0

  ! -- FreeSurfaceSolver state (mirrors the SAVE list) -------------------
  character(len=MAX_NAME_LEN), save :: acc_sname = ''
  character(len=MAX_NAME_LEN), save :: varname = ''
  integer(c_int), save :: prevsize = 0
  type(valuelist_t), pointer, save :: sparams => null()
  type(matrix_t), pointer, save :: sysmat => null()
  real(c_double), pointer, save :: forcevec(:) => null()
  real(c_double), pointer, save :: freesurf(:) => null()
  real(c_double), pointer, save :: prefsurf(:,:) => null()
  integer, pointer, save :: fsperm(:) => null()
  type(element_t), pointer, save :: elem => null()
  type(valuelist_t), pointer, save :: equation => null()
  type(valuelist_t), pointer, save :: material => null()
  type(valuelist_t), pointer, save :: bodyforce => null()
  type(variable_t), pointer, save :: varsurfres => null()
  real(c_double), pointer, save :: surfresptr(:) => null()
  integer, pointer, save :: nidx(:) => null()

  real(c_double), allocatable, save :: stiff(:,:)
  real(c_double), allocatable, save :: mass(:,:)
  real(c_double), allocatable, save :: sourcefunc(:)
  real(c_double), allocatable, save :: force(:)
  real(c_double), allocatable, save :: timeforce(:)
  real(c_double), allocatable, save :: velo(:,:)
  real(c_double), allocatable, save :: meshvelo(:,:)
  real(c_double), allocatable, save :: flux(:,:)
  real(c_double), allocatable, save :: elemfs(:)
  real(c_double), allocatable, save :: localmaxdisp(:)
  real(c_double), allocatable, save :: lowerlim(:)
  real(c_double), allocatable, save :: upperlim(:)
  real(c_double), allocatable, save :: oldvals(:)
  real(c_double), allocatable, save :: oldrhs(:)
  real(c_double), allocatable, save :: stiffvec(:)
  real(c_double), allocatable, save :: resvec(:)
  logical, allocatable, save :: limsol(:,:)
  logical, allocatable, save :: actnode(:,:)
  type(nodes_t), save :: enodes

  integer(c_int), save :: acc_dim = 0
  integer(c_int), save :: acc_k = 0
  integer(c_int), save :: acc_l = 0
  integer(c_int), save :: acc_realloc = 0
  integer(c_int), save :: acc_nliter = 1
  integer(c_int), save :: acc_apply = 0
  integer(c_int), save :: acc_needold = 0
  integer(c_int), save :: acc_limitdisp = 0
  integer(c_int), save :: acc_complocal = 0
  integer(c_int), save :: acc_ale = 0
  integer(c_int), save :: acc_bubbles = 1
  integer(c_int), save :: acc_rotate = 0
  integer(c_int), save :: acc_reset = 0
  integer(c_int), save :: acc_allocdone = 0
  real(c_double), save :: acc_cv = 1.0_c_double
  real(c_double), save :: acc_lintol = 0.0_c_double
  real(c_double), save :: acc_nltol = 0.0_c_double
  real(c_double), save :: acc_maxdisp = 0.0_c_double
  real(c_double), save :: acc_relax = 1.0_c_double
  real(c_double), save :: acc_orient(3) = 0.0_c_double
  real(c_double), save :: acc_rotmat(3,3) = 0.0_c_double
  real(c_double), save :: acc_at = 0.0_c_double
  real(c_double), save :: acc_st = 0.0_c_double
  real(c_double), save :: acc_totat = 0.0_c_double
  real(c_double), save :: acc_totst = 0.0_c_double
  real(c_double), save :: acc_relchange = 0.0_c_double
  real(c_double), save :: acc_maxdh = 0.0_c_double
  integer(c_int), save :: acc_nsurf = 0
  ! per-element scratch (set in fss_elem_velo, consumed by fss_local_matrix)
  integer(c_int), save :: acc_nsdofs = 0
  type(variable_t), pointer, save :: acc_flowsol => null()
  integer, pointer, save :: acc_flowperm(:) => null()
  real(c_double), pointer, save :: acc_flowvals(:) => null()
  real(c_double), pointer, save :: acc_flowprev(:,:) => null()
  character(len=MAX_NAME_LEN), save :: acc_flowsolname = ''
  ! LocalMatrix scratch
  type(ElementType_t), pointer, save :: acc_savetype => null()
  integer, pointer, save :: acc_edgemap(:,:) => null()
  logical, save :: acc_normalflux = .true.
  logical, save :: acc_subst = .true.
  logical, save :: acc_bubbles_lm = .true.

contains

  ! ================= RotInit entry =================
  subroutine fsr_begin(h_model, h_solver) bind(c, name="fsr_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine fsr_begin

  subroutine fsr_params() bind(c, name="fsr_params")
    logical :: found
    integer :: i
    real(c_double), pointer :: parray(:,:)
    type(variable_t), pointer :: orientvar
    sparams => GetSolverParams()
    rot_sname = 'FreeSurfaceSolver_RotInit'
    rot_varname = ListGetString(sparams, 'Free Surface Variable Name', found)
    if (.not. found) then
      call Fatal(rot_sname, "Can't find 'Free Surface Variable Name'")
    end if
    call ListAddString(acc_model % Constants, 'Free Surface Variable Name', &
        trim(rot_varname))
    parray => ListGetConstRealArray(sparams, trim(rot_varname)//' Orientation', found)
    if (.not. found) then
      write(message, '(a,a,a)') 'Solver to Initialize Rotated Free Surface is present but &
          &can''t find ', rot_varname, ' Orientation'
      call Fatal(rot_sname, message)
    end if
    do i = 1, 3
      rot_orient(i) = parray(i,1)
    end do
    write(rot_ovname, '(a,a)') trim(rot_varname), ' Orientation'
    if (.not. associated(rot_ovp)) allocate(rot_ovp(3))
    call VariableAdd(acc_model % Mesh % Variables, acc_model % Mesh, acc_solver, &
        rot_ovname, 3, rot_ovp)
    orientvar => VariableGet(acc_model % Mesh % Variables, rot_ovname, .TRUE.)
    if (.not. associated(orientvar)) &
        call Fatal(rot_sname, 'Internal problem set/getting Orientation Variable')
    orientvar % Values = rot_orient
    rot_rotmat = ComputeRotationMatrix(rot_orient)
  end subroutine fsr_params

  subroutine fsr_bind() bind(c, name="fsr_bind")
    type(variable_t), pointer :: varfs
    rot_ref => acc_solver % Variable % Values
    if (.not. associated(rot_ref)) &
        call Fatal(rot_sname, 'Variable values not associated')
    rot_rperm => acc_solver % Variable % Perm
    varfs => VariableGet(acc_model % Mesh % Variables, trim(rot_varname))
    if (.not. associated(varfs)) then
      write(message, '(a,a)') "Can't get variable: ", rot_varname
      call Fatal(rot_sname, message)
    end if
    rot_fs => varfs % Values
    rot_fperm => varfs % Perm
  end subroutine fsr_bind

  function fsr_num_nodes() result(n) bind(c, name="fsr_num_nodes")
    integer(c_int) :: n
    n = acc_model % Mesh % NumberOfNodes
  end function fsr_num_nodes

  function fsr_ref_perm(i) result(k) bind(c, name="fsr_ref_perm")
    integer(c_int), value :: i
    integer(c_int) :: k
    k = rot_rperm(i)
  end function fsr_ref_perm

  subroutine fsr_rot_node(i) bind(c, name="fsr_rot_node")
    integer(c_int), value :: i
    real(c_double) :: holder(3)
    holder(1) = acc_model % Mesh % Nodes % x(i)
    holder(2) = acc_model % Mesh % Nodes % y(i)
    holder(3) = acc_model % Mesh % Nodes % z(i)
    holder = matmul(rot_rotmat, holder)
    rot_ref(rot_rperm(i)) = holder(3)
    rot_fs(rot_fperm(i)) = holder(3)
  end subroutine fsr_rot_node

  ! ================= MeshUpdate USF core =================
  function fss_mu(h_model, h_node, arr, axis) result(mu) bind(c, name="fss_mu")
    type(c_ptr), value :: h_model
    integer(c_int) :: h_node
    real(c_double) :: arr(2)
    integer(c_int), value :: axis
    real(c_double) :: mu
    type(model_t), pointer :: m
    type(variable_t), pointer :: orientvar, timevar
    real(c_double) :: fs, ref_fs, diff, unrot(3,3), orient(3), diffvect(3), muvect(3)
    character(len=MAX_NAME_LEN) :: usf_name, varname_l
    logical :: found
    ! NOTE: mirrors the original's SAVE (FirstTime, t0) without explicit
    ! initialization, so the zero-fill behavior of the serial build is
    ! reproduced exactly (the FirstTime branch is never taken).
    call c_f_pointer(h_model, m)
    fs = arr(1)
    ref_fs = arr(2)
    usf_name = 'FreeSurfaceToMeshUpdate'
    timevar => VariableGet(m % Mesh % Variables, 'time')
    if (mu_first) then
      if (.not. found) call Fatal(usf_name, "Can't get time variable")
      mu_t0 = timevar % Values(1)
      mu_first = .false.
    end if
    varname_l = ListGetString(m % Constants, 'Free Surface Variable Name', found)
    if (.not. found) then
      if (mu_t0 == timevar % Values(1)) then
        mu = 0
        return
      else
        call Fatal(usf_name, "Couldn't find Free Surface Variable Name")
      end if
    end if
    orientvar => VariableGet(m % Mesh % Variables, trim(varname_l)//' Orientation', found)
    if (.not. found) then
      write(message, '(a,a,a)') "Can't find ", varname_l, " Orientation"
      call Fatal(usf_name, message)
    end if
    orient = orientvar % Values
    unrot = transpose(ComputeRotationMatrix(orient))
    diff = fs - ref_fs
    diffvect = [0.0_c_double, 0.0_c_double, diff]
    muvect = matmul(unrot, diffvect)
    mu = muvect(axis)
  end function fss_mu

  ! ================= main entry =================
  subroutine fss_begin(h_model, h_solver) bind(c, name="fss_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    sparams => GetSolverParams()
    sysmat => acc_solver % Matrix
    forcevec => acc_solver % Matrix % RHS
    freesurf => acc_solver % Variable % Values
    if (.not. associated(freesurf)) &
        call Fatal(acc_sname, 'Variable values not associated')
    fsperm => acc_solver % Variable % Perm
    prefsurf => acc_solver % Variable % PrevValues
  end subroutine fss_begin

  subroutine fss_info_header() bind(c, name="fss_info_header")
    acc_sname = 'FreeSurfaceSolver ('//trim(acc_solver % Variable % Name)//')'
    call Info(acc_sname, 'Solving for free surface', level=5)
  end subroutine fss_info_header

  function fss_realloc() result(r) bind(c, name="fss_realloc")
    integer(c_int) :: r
    acc_realloc = 0
    if (varname /= trim(acc_solver % Variable % Name)) then
      varname = trim(acc_solver % Variable % Name)
      acc_realloc = 1
    end if
    if (size(acc_solver % Variable % Values) /= prevsize) then
      prevsize = size(acc_solver % Variable % Values)
      acc_realloc = 1
    end if
    r = acc_realloc
  end function fss_realloc

  function fss_free_count() result(n) bind(c, name="fss_free_count")
    integer(c_int) :: n
    n = count(fsperm /= 0)
  end function fss_free_count

  subroutine fss_warn_nofree() bind(c, name="fss_warn_nofree")
    if (parenv % PEs > 1) then
      write(message, '(A,I0,A)') 'Partition ', parenv % myPE, ' has no free surface'
      call Warn(acc_sname, message)
    else
      call Warn(acc_sname, &
          'A serial run without a free surface, but the solver switched in - weird!')
    end if
  end subroutine fss_warn_nofree

  subroutine fss_read() bind(c, name="fss_read")
    logical :: found
    character(len=MAX_NAME_LEN) :: stabflag
    acc_dim = CoordinateSystemDimension()
    acc_cv = GetCReal(sparams, 'Velocity Implicity', found)
    if (.not. found) acc_cv = 1.0_c_double
    write(message, '(a,F9.2)') 'Velocity implicity (1=fully implicit)=', acc_cv
    call Info(acc_sname, message, level=10)
    acc_lintol = GetConstReal(sparams, 'Linear System Convergence Tolerance', found)
    if (.not. found) then
      call Fatal(acc_sname, 'No > Linear System Convergence Tolerance < found')
    end if
    acc_nltol = GetConstReal(sparams, 'Nonlinear System Convergence Tolerance', found)
    acc_nliter = GetInteger(sparams, 'Nonlinear System Max Iterations', found)
    if (.not. found) acc_nliter = 1
    acc_maxdisp = GetConstReal(sparams, 'Maximum Displacement', found)
    acc_limitdisp = merge(1, 0, found)
    acc_complocal = merge(1, 0, GetLogical(sparams, &
        'Compute Local Maximum Displacement', found))
    if (found) then
      call Info(acc_sname, &
          '"Compute Local Maximum Displacement" set - ignoring global limitation, if given', &
          level=3)
      acc_limitdisp = 1
    end if
    acc_relax = GetCReal(sparams, 'Relaxation Factor', found)
    if (.not. found) acc_relax = 1.0_c_double
    if ((found .and. (acc_relax < 1.0_c_double)) .or. (acc_limitdisp /= 0)) then
      acc_needold = 1
    else
      acc_needold = 0
    end if
    if (GetLogical(sparams, 'Apply Dirichlet', found)) then
      acc_apply = 1
      call Info(acc_sname, 'Using Dirichlet method for limitation', level=6)
      if (GetLogical(sparams, 'Reset Limiter', found)) then
        acc_reset = 1
        call Info(acc_sname, 'Limiters will be reset for each nonlinear iteration', level=3)
      else
        acc_reset = 0
      end if
      if (acc_nliter < 2) then
        call Warn(acc_sname, &
            'Keyword > Apply Dirichlet < set, but > Nonlinear System Max Iterations < set to lower than 2')
      end if
    else
      acc_apply = 0
      if (.not. found) then
        call Info(acc_sname, 'No keyword > Apply Dirichlet < found. No limitation of solution', &
            level=6)
      else
        call Info(acc_sname, 'No limitation of solution', level=6)
      end if
    end if
    if (GetLogical(sparams, 'ALE Formulation', found)) then
      acc_ale = 1
      call Info(acc_sname, 'Using horizontal ALE Formulation', level=6)
    else
      acc_ale = 0
      call Info(acc_sname, 'Using horizontal Eulerian Formulation', level=6)
    end if
    stabflag = GetString(sparams, 'Stabilization Method', found)
    select case (stabflag)
    case ('stabilized')
      acc_bubbles = 0
      call Info(acc_sname, 'Using residual squared-stabilized formulation.', level=6)
    case ('bubbles')
      acc_bubbles = 1
      call Info(acc_sname, 'Using residual free bubble stabilization', level=6)
    end select
    acc_bubbles_lm = (acc_bubbles /= 0)
    acc_normalflux = .true.
    acc_subst = .true.
  end subroutine fss_read

  function fss_apply() result(r) bind(c, name="fss_apply")
    integer(c_int) :: r
    r = acc_apply
  end function fss_apply

  function fss_need_old() result(r) bind(c, name="fss_need_old")
    integer(c_int) :: r
    r = acc_needold
  end function fss_need_old

  function fss_nliter() result(r) bind(c, name="fss_nliter")
    integer(c_int) :: r
    r = acc_nliter
  end function fss_nliter

  function fss_rot_check() result(r) bind(c, name="fss_rot_check")
    integer(c_int) :: r
    type(variable_t), pointer :: rotvar
    integer :: i
    rotvar => VariableGet(acc_model % Mesh % Variables, trim(varname)//' Orientation')
    if (associated(rotvar)) then
      acc_rotate = 1
      do i = 1, 3
        acc_orient(i) = rotvar % Values(i)
      end do
      acc_rotmat = ComputeRotationMatrix(acc_orient)
      write(message, '(A,f8.2,f8.2,f8.2)') 'Rotated Free Surface defined using vector: ', acc_orient
      call Info(acc_sname, message, level=6)
    else
      acc_rotate = 0
      call Info(acc_sname, &
          'No Free Surface Orientation Vector found, assuming normal to z-axis', level=6)
    end if
    r = acc_rotate
  end function fss_rot_check

  subroutine fss_info_dims() bind(c, name="fss_info_dims")
    write(message, '(A,I0)') 'Mesh dimension: ', acc_dim
    call Info(acc_sname, message, level=8)
    write(message, '(A,I0)') 'Size of problem: ', size(acc_solver % Variable % Values)
    call Info(acc_sname, message, level=8)
  end subroutine fss_info_dims

  subroutine fss_allocate() bind(c, name="fss_allocate")
    integer :: nmax, mmax, nmatrix, istat
    if ((acc_allocdone == 0) .or. (acc_solver % MeshChanged) .or. (acc_realloc /= 0)) then
      call Info(acc_sname, 'Allocating variables for free surface', level=10)
      nmax = acc_model % MaxElementNodes
      mmax = acc_model % Mesh % NumberOfNodes
      acc_k = size(sysmat % Values)
      acc_l = size(sysmat % RHS)
      if (acc_allocdone /= 0) then
        deallocate(enodes % x, enodes % y, enodes % z, timeforce, force, stiff, &
            mass, velo, meshvelo, flux, elemfs, sourcefunc)
        if (acc_complocal /= 0) deallocate(localmaxdisp)
        if (acc_apply /= 0) deallocate(lowerlim, upperlim, limsol, actnode, &
            resvec, stiffvec, oldvals, oldrhs)
      end if
      if (acc_bubbles /= 0) then
        nmatrix = 2 * nmax
      else
        nmatrix = nmax
      end if
      allocate(enodes % x(nmax), enodes % y(nmax), enodes % z(nmax), &
          timeforce(nmatrix), force(nmatrix), stiff(nmatrix, nmatrix), &
          mass(nmatrix, nmatrix), velo(3, nmax), meshvelo(3, nmax), &
          flux(3, nmax), elemfs(nmax), sourcefunc(nmax), stat=istat)
      if (istat /= 0) call Fatal(acc_sname, 'Memory allocation error 1, Aborting.')
      elemfs = 0.0_c_double
      if (acc_complocal /= 0) then
        allocate(localmaxdisp(mmax), stat=istat)
        if (istat /= 0) call Fatal(acc_sname, 'Memory allocation error 3, Aborting.')
      end if
      if (acc_apply /= 0) then
        allocate(lowerlim(mmax), upperlim(mmax), limsol(mmax, 2), actnode(mmax, 2), &
            resvec(acc_l), stiffvec(acc_l), oldvals(acc_k), oldrhs(acc_l), stat=istat)
        if (istat /= 0) call Fatal(acc_sname, 'Memory allocation error 4, Aborting.')
        actnode = .false.
        limsol = .false.
        resvec = 0.0_c_double
      end if
      call Info(acc_sname, 'Memory allocations done')
      acc_allocdone = 1
    end if
  end subroutine fss_allocate

  subroutine fss_bind_residual() bind(c, name="fss_bind_residual")
    varsurfres => VariableGet(acc_model % Mesh % Variables, trim(varname)//' Residual')
    if (.not. associated(varsurfres)) then
      write(message, '(A)') '>'//trim(varname)//' Residual < not associated'
      call Fatal(acc_sname, message)
    end if
    surfresptr => varsurfres % Values
  end subroutine fss_bind_residual

  subroutine fss_maybe_reset() bind(c, name="fss_maybe_reset")
    if (acc_reset /= 0) actnode = .false.
  end subroutine fss_maybe_reset

  subroutine fss_iter_begin(it, nlit) bind(c, name="fss_iter_begin")
    integer(c_int), value :: it, nlit
    write(message, '(a,I0,a,I0)') 'Non-linear Iteration ', it, ' out of max. ', nlit
    call Info(acc_sname, message, level=4)
    acc_totat = 0.0_c_double
    acc_totst = 0.0_c_double
    acc_at = CPUTime()
    call Info(acc_sname, 'start assembly', level=6)
    call DefaultInitialize()
  end subroutine fss_iter_begin

  function fss_nactive() result(n) bind(c, name="fss_nactive")
    integer(c_int) :: n
    n = acc_solver % NumberOfActiveElements
  end function fss_nactive

  subroutine fss_elem_begin(t, n, skip) bind(c, name="fss_elem_begin")
    integer(c_int), value :: t
    integer(c_int) :: n, skip
    integer :: i
    real(c_double) :: holder(3)
    skip = 0
    elem => GetActiveElement(t)
    n = GetElementNOFNodes()
    if (GetElementFamily() == 1) then
      skip = 1
      return
    end if
    nidx => elem % NodeIndexes
    enodes % x(1:n) = acc_solver % Mesh % Nodes % x(nidx)
    enodes % y(1:n) = acc_solver % Mesh % Nodes % y(nidx)
    enodes % z(1:n) = acc_solver % Mesh % Nodes % z(nidx)
    if (acc_rotate /= 0) then
      do i = 1, n
        holder(1) = enodes % x(i)
        holder(2) = enodes % y(i)
        holder(3) = enodes % z(i)
        holder = matmul(acc_rotmat, holder)
        enodes % x(i) = holder(1)
        enodes % y(i) = holder(2)
        enodes % z(i) = holder(3)
      end do
    end if
    enodes % z(1:n) = 0.0_c_double
    if (acc_dim == 2) then
      enodes % y(1:n) = 0.0_c_double
    else if (acc_dim /= 3) then
      write(message, '(a,i0,a)') &
          'It is not possible to compute free-surface problems in DIM=', &
          acc_dim, ' dimensions. Aborting'
      call Fatal(acc_sname, message)
    end if
    equation => GetEquation()
    material => GetMaterial()
    bodyforce => GetBodyForce()
  end subroutine fss_elem_begin

  subroutine fss_elem_limits(n) bind(c, name="fss_elem_limits")
    integer(c_int), value :: n
    logical :: found
    lowerlim(elem % NodeIndexes(1:n)) = &
        ListGetReal(material, 'Min '//trim(varname), n, elem % NodeIndexes, found)
    limsol(elem % NodeIndexes(1:n), 1) = found
    upperlim(elem % NodeIndexes(1:n)) = &
        ListGetReal(material, 'Max '//trim(varname), n, elem % NodeIndexes, found)
    limsol(elem % NodeIndexes(1:n), 2) = found
  end subroutine fss_elem_limits

  function fss_use_local() result(r) bind(c, name="fss_use_local")
    integer(c_int) :: r
    r = acc_complocal
  end function fss_use_local

  subroutine fss_elem_maxdisp(n) bind(c, name="fss_elem_maxdisp")
    integer(c_int), value :: n
    logical :: found
    localmaxdisp(elem % NodeIndexes(1:n)) = &
        ListGetReal(bodyforce, 'Max Displacement '//trim(varname), n, &
        elem % NodeIndexes, found)
    if (.not. found) then
      write(message, *) '"Compute Local Max Displacement" set but ', &
          'Max Displacement'//trim(varname), ' not found'
      call Warn(acc_sname, message)
      acc_complocal = 0
    end if
  end subroutine fss_elem_maxdisp

  subroutine fss_elem_velo(n, transient) bind(c, name="fss_elem_velo")
    integer(c_int), value :: n, transient
    character(len=MAX_NAME_LEN) :: convflag
    logical :: found
    integer :: i, j, nsdofs
    convflag = GetString(equation, 'Convection', found)
    if (.not. found) &
        call Fatal(acc_sname, 'No string for keyword > Convection < found in Equation')
    velo = 0.0_c_double
    if (convflag == 'constant') then
      velo(1,1:n) = GetReal(material, 'Convection Velocity 1', found)
      if (.not. found) velo(1,1:n) = GetReal(equation, 'Convection Velocity 1', found)
      velo(2,1:n) = GetReal(material, 'Convection Velocity 2', found)
      if (.not. found) velo(2,1:n) = GetReal(equation, 'Convection Velocity 2', found)
      velo(3,1:n) = GetReal(material, 'Convection Velocity 3', found)
      if (.not. found) velo(3,1:n) = GetReal(equation, 'Convection Velocity 3', found)
    else if (convflag == 'computed') then
      acc_flowsolname = GetString(equation, 'Flow Solution Name', found)
      if (.not. found) then
        call Warn(acc_sname, 'Keyword > Flow Solution Name < not found in section >Equation<')
        call Warn(acc_sname, 'Taking default value > Flow Solution <')
      end if
      acc_flowsol => VariableGet(acc_solver % Mesh % Variables, acc_flowsolname)
      if (associated(acc_flowsol)) then
        acc_flowperm => acc_flowsol % Perm
        acc_nsdofs = acc_flowsol % DOFs
        acc_flowvals => acc_flowsol % Values
        acc_flowprev => acc_flowsol % PrevValues
      else
        write(message, '(A,A,A)') &
            'Convection flag set to > computed <, but no variable >', acc_flowsolname, '< found'
        call Fatal(acc_sname, message)
      end if
      if (associated(acc_flowsol)) then
        nsdofs = acc_nsdofs
        do i = 1, n
          j = nsdofs * acc_flowperm(nidx(i))
          if ((transient /= 0) .and. (abs(acc_cv - 1.0_c_double) > 0.001_c_double)) then
            if ((acc_dim == 2) .and. (nsdofs == 3)) then
              velo(1,i) = acc_cv * acc_flowvals(j-2) + (1 - acc_cv) * acc_flowprev(j-2,1)
              velo(2,i) = acc_cv * acc_flowvals(j-1) + (1 - acc_cv) * acc_flowprev(j-1,1)
              velo(3,i) = 0.0_c_double
            else if ((acc_dim == 3) .and. (nsdofs == 4)) then
              velo(1,i) = acc_cv * acc_flowvals(j-3) + (1 - acc_cv) * acc_flowprev(j-3,1)
              velo(2,i) = acc_cv * acc_flowvals(j-2) + (1 - acc_cv) * acc_flowprev(j-2,1)
              velo(3,i) = acc_cv * acc_flowvals(j-1) + (1 - acc_cv) * acc_flowprev(j-1,1)
            else if ((CurrentCoordinateSystem() == CylindricSymmetric) &
                .and. (acc_dim == 2) .and. (nsdofs == 4)) then
              velo(1,i) = acc_cv * acc_flowvals(j-3) + (1 - acc_cv) * acc_flowprev(j-3,1)
              velo(2,i) = acc_cv * acc_flowvals(j-2) + (1 - acc_cv) * acc_flowprev(j-2,1)
              velo(3,i) = acc_cv * acc_flowvals(j-1) + (1 - acc_cv) * acc_flowprev(j-1,1)
            else
              write(message, '(a,i0,a,i0,a)') 'DIM=', acc_dim, ' NSDOFs=', nsdofs, &
                  ' does not combine. Aborting'
              call Fatal(acc_sname, message)
            end if
          else
            if ((acc_dim == 2) .and. (nsdofs == 3)) then
              velo(1,i) = acc_flowvals(j-2)
              velo(2,i) = acc_flowvals(j-1)
              velo(3,i) = 0.0_c_double
            else if ((acc_dim == 3) .and. (nsdofs == 4)) then
              velo(1,i) = acc_flowvals(j-3)
              velo(2,i) = acc_flowvals(j-2)
              velo(3,i) = acc_flowvals(j-1)
            else if ((CurrentCoordinateSystem() == CylindricSymmetric) &
                .and. (acc_dim == 2) .and. (nsdofs == 4)) then
              velo(1,i) = acc_flowvals(j-3)
              velo(2,i) = acc_flowvals(j-2)
              velo(3,i) = acc_flowvals(j-1)
            else
              write(message, '(a,i0,a,i0,a)') 'DIM=', acc_dim, ' NSDOFs=', nsdofs, &
                  ' does not combine. Aborting'
              call Fatal(acc_sname, message)
            end if
          end if
        end do
      else
        velo = 0.0_c_double
      end if
    else
      velo = 0.0_c_double
    end if
    if (acc_rotate /= 0) then
      do i = 1, n
        velo(:,i) = matmul(acc_rotmat, velo(:,i))
      end do
    end if
  end subroutine fss_elem_velo

  subroutine fss_elem_meshvelo() bind(c, name="fss_elem_meshvelo")
    meshvelo = 0.0_c_double
    call GetVectorLocalSolution(meshvelo, 'Mesh Velocity', elem)
  end subroutine fss_elem_meshvelo

  subroutine fss_elem_source(n) bind(c, name="fss_elem_source")
    integer(c_int), value :: n
    logical :: found
    integer :: i
    sourcefunc = 0.0_c_double
    flux = 0.0_c_double
    acc_subst = .true.
    if (associated(bodyforce)) then
      acc_subst = .false.
      sourcefunc(1:n) = GetReal(bodyforce, trim(varname)//' Accumulation', acc_normalflux)
      if (.not. acc_normalflux) then
        flux(1,1:n) = GetReal(bodyforce, trim(varname)//' Accumulation Flux 1', found)
        if (.not. found) flux(1,1:n) = 0.0_c_double
        if (acc_dim >= 2) then
          flux(2,1:n) = GetReal(bodyforce, trim(varname)//' Accumulation Flux 2', found)
          if (.not. found) flux(2,1:n) = 0.0_c_double
        else
          flux(2,1:n) = 0.0_c_double
        end if
        if (acc_dim == 3) then
          flux(3,1:n) = GetReal(bodyforce, trim(varname)//' Accumulation Flux 3', found)
          if (.not. found) flux(3,1:n) = 0.0_c_double
        else
          flux(3,1:n) = 0.0_c_double
        end if
        sourcefunc = 0.0_c_double
        if (acc_rotate /= 0) then
          do i = 1, n
            flux(:,i) = matmul(acc_rotmat, flux(:,i))
          end do
        end if
      end if
    end if
  end subroutine fss_elem_source

  subroutine fss_elem_old(n) bind(c, name="fss_elem_old")
    integer(c_int), value :: n
    elemfs(1:n) = prefsurf(fsperm(nidx), 1)
  end subroutine fss_elem_old

  subroutine fss_local_matrix(ncoord, transient) bind(c, name="fss_local_matrix")
    integer(c_int), value :: ncoord, transient
    real(c_double) :: basis(2*ncoord), dbasis(2*ncoord,3)
    real(c_double) :: vgauss(3), vmeshgauss(3), source, gradfs(3), normg, &
        fluxgauss(3), x, y, z, u, v, w, s, sqm, su(2*ncoord), sw(2*ncoord), &
        tau, hk, unorm
    type(ElementType_t), pointer :: savetype
    integer :: lintype(2:4) = [202,303,404]
    logical :: stat, uselin
    integer :: i, j, t, p, q, n
    type(GaussIntegrationPoints_t) :: integ
    force = 0.0_c_double
    stiff = 0.0_c_double
    mass = 0.0_c_double
    if (acc_bubbles /= 0) then
      n = ncoord * 2
    else
      n = ncoord
    end if
    uselin = GetLogical(GetSolverParams(), 'Use linear elements', stat)
    if (acc_apply /= 0) then
      uselin = uselin .or. any(actnode(nidx,:))
    end if
    uselin = uselin .and. (elem % TYPE % BasisFunctionDegree == 2)
    if (uselin) then
      savetype => elem % TYPE
      elem % TYPE => GetElementType(lintype(GetElementFamily()))
    end if
    hk = ElementDiameter(elem, enodes)
    if (acc_bubbles /= 0) then
      integ = GaussPoints(elem, elem % TYPE % gausspoints2)
    else
      integ = GaussPoints(elem)
    end if
    su = 0.0_c_double
    sw = 0.0_c_double
    do t = 1, integ % n
      u = integ % u(t)
      v = integ % v(t)
      w = integ % w(t)
      s = integ % s(t)
      stat = ElementInfo(elem, enodes, u, v, w, sqm, basis, dbasis, &
          Bubbles=(acc_bubbles /= 0))
      s = s * sqm
      if (CurrentCoordinateSystem() /= Cartesian) then
        x = sum(enodes % x(1:ncoord) * basis(1:ncoord))
        y = sum(enodes % y(1:ncoord) * basis(1:ncoord))
        z = sum(enodes % z(1:ncoord) * basis(1:ncoord))
        s = s * x
      end if
      gradfs = 0.0_c_double
      vgauss = 0.0_c_double
      vmeshgauss = 0.0_c_double
      do i = 1, acc_dim - 1
        gradfs(i) = sum(dbasis(1:ncoord,i) * elemfs(1:ncoord))
      end do
      gradfs(acc_dim) = 1.0_c_double
      if (acc_ale == 0) then
        do i = 1, acc_dim
          vgauss(i) = sum(basis(1:ncoord) * (velo(i,1:ncoord)))
        end do
      else
        meshvelo(acc_dim,1:ncoord) = 0.0_c_double
        do i = 1, acc_dim
          vgauss(i) = sum(basis(1:ncoord) * (velo(i,1:ncoord) - meshvelo(i,1:ncoord)))
        end do
      end if
      if (acc_dim == 3) then
        normg = sqrt(1.0_c_double + gradfs(1)**2 + gradfs(2)**2)
      else
        normg = sqrt(1.0_c_double + gradfs(1)**2)
      end if
      unorm = sqrt(sum(vgauss(1:acc_dim-1)**2))
      if (unorm /= 0.0_c_double) then
        tau = hk / (2 * unorm)
      else
        tau = 0.0_c_double
      end if
      if (acc_bubbles == 0) then
        do p = 1, n
          su(p) = 0.0_c_double
          do i = 1, acc_dim - 1
            su(p) = su(p) + vgauss(i) * dbasis(p,i)
          end do
          sw(p) = 0.0_c_double
          do i = 1, acc_dim - 1
            sw(p) = sw(p) + vgauss(i) * dbasis(p,i)
          end do
        end do
      end if
      do p = 1, n
        do q = 1, n
          do i = 1, acc_dim - 1
            stiff(p,q) = stiff(p,q) + s * vgauss(i) * dbasis(q,i) * basis(p)
          end do
          stiff(p,q) = stiff(p,q) + s * tau * su(q) * sw(p)
        end do
      end do
      if (transient /= 0) then
        do p = 1, n
          do q = 1, n
            mass(p,q) = mass(p,q) + s * basis(q) * (basis(p) + tau * sw(p))
          end do
        end do
      end if
      if (.not. acc_subst) then
        if (acc_normalflux) then
          source = normg * sum(sourcefunc(1:ncoord) * basis(1:ncoord))
        else
          do i = 1, acc_dim
            fluxgauss(i) = sum(basis(1:ncoord) * flux(i,1:ncoord))
          end do
          source = sum(fluxgauss(1:acc_dim) * gradfs(1:acc_dim))
        end if
      else
        source = 0.0_c_double
      end if
      force(1:n) = force(1:n) + (vgauss(acc_dim) + source) * (basis(1:n) + tau * sw(1:n)) * s
    end do
    if (uselin) then
      acc_edgemap => GetEdgeMap(GetElementFamily())
      n = elem % TYPE % NumberOfNodes
      do i = n + 1, n + size(acc_edgemap, 1)
        j = acc_edgemap(i-n,1)
        q = acc_edgemap(i-n,2)
        stiff(i,:) = 0.0_c_double
        stiff(:,i) = 0.0_c_double
        mass(i,:) = 0.0_c_double
        mass(:,i) = 0.0_c_double
        stiff(i,i) = 1.0_c_double
        stiff(i,j) = -0.5_c_double
        stiff(i,q) = -0.5_c_double
        force(i) = 0.0_c_double
        elem % TYPE => savetype
      end do
    end if
  end subroutine fss_local_matrix

  subroutine fss_elem_time(transient) bind(c, name="fss_elem_time")
    integer(c_int), value :: transient
    timeforce = 0.0_c_double
    if (transient /= 0) then
      call Default1stOrderTime(mass, stiff, force)
    end if
  end subroutine fss_elem_time

  subroutine fss_elem_update(n) bind(c, name="fss_elem_update")
    integer(c_int), value :: n
    if (acc_bubbles /= 0) call Condensate(n, stiff, force, timeforce)
    call DefaultUpdateEquations(stiff, force)
  end subroutine fss_elem_update

  subroutine fss_finish() bind(c, name="fss_finish")
    call DefaultFinishBulkAssembly()
    call DefaultFinishAssembly()
    call DefaultDirichletBCs()
  end subroutine fss_finish

  subroutine fss_limits_pre() bind(c, name="fss_limits_pre")
    integer :: i, k
    oldvals = sysmat % Values
    oldrhs = forcevec
    do i = 1, acc_model % NumberOfNodes
      k = fsperm(i)
      if ((actnode(i,1) .and. actnode(i,2))) &
          call Fatal(acc_sname, 'Upper as well as lower limiter active - this is a deadlock')
      if (((actnode(i,1) .or. actnode(i,2))) .and. (k > 0)) then
        call ZeroRow(sysmat, k)
        call SetMatrixElement(sysmat, k, k, 1.0_c_double)
        if (actnode(i,1)) then
          sysmat % RHS(k) = lowerlim(i)
        else
          sysmat % RHS(k) = upperlim(i)
        end if
      end if
    end do
  end subroutine fss_limits_pre

  subroutine fss_info_assdone() bind(c, name="fss_info_assdone")
    call Info(acc_sname, 'Assembly done', level=6)
  end subroutine fss_info_assdone

  function fss_solve(relchange) result(norm) bind(c, name="fss_solve")
    real(c_double) :: relchange
    real(c_double) :: norm
    real(c_double) :: prevnorm
    acc_at = CPUTime() - acc_at
    acc_st = CPUTime()
    prevnorm = acc_solver % Variable % Norm
    norm = DefaultSolve()
    if (prevnorm + norm /= 0.0_c_double) then
      acc_relchange = 2.0_c_double * abs(prevnorm - norm) / (prevnorm + norm)
    else
      acc_relchange = 0.0_c_double
    end if
    relchange = acc_relchange
    write(message, *) 'Result Norm   : ', norm
    call Info(acc_sname, message, level=4)
    write(message, *) 'Relative Change : ', acc_relchange
    call Info(acc_sname, message, level=4)
  end function fss_solve

  subroutine fss_residual(it) bind(c, name="fss_residual")
    integer(c_int), value :: it
    integer :: i, l
    acc_nsurf = 0
    sysmat % Values = oldvals
    forcevec = oldrhs
    if (parenv % PEs > 1) then
      call ParallelInitSolve(sysmat, freesurf, forcevec, resvec)
      call ParallelMatrixVector(sysmat, freesurf, stiffvec, .TRUE.)
      resvec = stiffvec - forcevec
      call ParallelSumVector(sysmat, resvec)
    else
      call CRS_MatrixVectorMultiply(sysmat, freesurf, stiffvec)
      resvec = stiffvec - forcevec
    end if
    do i = 1, acc_model % NumberOfNodes
      l = fsperm(i)
      if (l < 1) cycle
      acc_nsurf = acc_nsurf + 1
      if ((limsol(i,1)) .and. (freesurf(l) - lowerlim(i) < 0.0_c_double)) then
        actnode(i,1) = .true.
      end if
      if ((limsol(i,2)) .and. (freesurf(l) - upperlim(i) > 0.0_c_double)) then
        actnode(i,2) = .true.
      end if
      if (limsol(i,1) .and. (resvec(l) < -acc_lintol) .and. (it > 1)) actnode(i,1) = .false.
      if (limsol(i,2) .and. (resvec(l) > acc_lintol) .and. (it > 1)) actnode(i,2) = .false.
      if (.not. actnode(i,1) .and. .not. actnode(i,2)) then
        surfresptr(varsurfres % Perm(i)) = 0.0_c_double
      else
        surfresptr(varsurfres % Perm(i)) = resvec(l)
      end if
    end do
  end subroutine fss_residual

  subroutine fss_relax() bind(c, name="fss_relax")
    integer :: i, j
    real(c_double) :: lrelax
    if (acc_limitdisp /= 0) then
      if (acc_complocal == 0) then
        acc_maxdh = -huge(acc_maxdh)
        do i = 1, acc_model % NumberOfNodes
          j = fsperm(i)
          if (j > 0) then
            acc_maxdh = max(acc_maxdh, abs(freesurf(j) - prefsurf(j,1)))
          end if
        end do
        acc_maxdh = ParallelReduction(acc_maxdh, 2)
        if (acc_maxdh > acc_maxdisp) then
          acc_relax = acc_relax * acc_maxdisp / acc_maxdh
        end if
        write(message, '(a,E9.2)') 'Maximum displacement ', acc_maxdh
        call Info(acc_sname, message, level=4)
      end if
    end if
    write(message, '(a,F9.2)') 'pp Relaxation factor', acc_relax
    call Info(acc_sname, message, level=4)
    do i = 1, acc_model % NumberOfNodes
      j = fsperm(i)
      if (j > 0) then
        lrelax = acc_relax
        if (acc_complocal /= 0) then
          acc_maxdh = abs(freesurf(j) - prefsurf(j,1))
          if (acc_maxdh > localmaxdisp(i)) lrelax = acc_relax * localmaxdisp(i) / acc_maxdh
        end if
        freesurf(j) = lrelax * freesurf(j) + (1 - lrelax) * prefsurf(j,1)
      end if
    end do
  end subroutine fss_relax

  function fss_iter_tail(it) result(done) bind(c, name="fss_iter_tail")
    integer(c_int), value :: it
    integer(c_int) :: done
    done = 0
    acc_st = CPUTIme() - acc_st
    acc_totat = acc_totat + acc_at
    acc_totst = acc_totst + acc_st
    write(message, '(a,F8.2,F8.2)') 'Assembly: (s)', acc_at, acc_totat
    call Info(acc_sname, message, level=4)
    write(message, '(a,F8.2,F8.2)') ' Solve:    (s)', acc_st, acc_totst
    call Info(acc_sname, message, level=4)
    write(message, '(a,e13.6,a,e13.6)') &
        'Max/min values surface:', maxval(freesurf(:)), '/', minval(freesurf(:))
    call Info(acc_sname, message, level=4)
    if (acc_apply /= 0) then
      write(message, '(a,i0)') 'Number of surface nodes: ', acc_nsurf
      call Info(acc_sname, message, level=4)
      write(message, '(a,i0)') 'Number of constrained points (lower limit): ', &
          count(actnode(:,1))
      call Info(acc_sname, message, level=4)
      write(message, '(a,i0)') 'Number of constrained points (upper limit): ', &
          count(actnode(:,2))
      call Info(acc_sname, message, level=4)
    end if
    if (acc_relchange < acc_nltol) then
      write(message, '(a,i0,a)') 'Converged after', it, ' iterations'
      call Info(acc_sname, message, level=4)
      done = 1
    end if
  end function fss_iter_tail

end module freesurfacesolver_acc
