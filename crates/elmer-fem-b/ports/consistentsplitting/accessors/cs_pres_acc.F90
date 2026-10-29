! ConsistentSplitting pilot: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/ConsistentSplitting/PressureSolver.F90.
! Same contract as cs_comp_acc (C-only arguments; SAVE state mirrored in
! Fortran; order-sensitive reductions stay in Fortran).
!
! The traction post-pass (boundary loop, parent-element search,
! SurfaceForceIntegration, Info/Message/ListAdd writes) runs entirely in
! one wrapper: it only touches Model/Solver state through the captured
! handles. The outflow-boundary matrix surgery runs in one wrapper per
! element. The serial tree never takes the Parallel arms; they are
! transcribed verbatim and compile-checked (parallel-tree differentials
! would drive them).
module cs_pres_acc
  use iso_c_binding
  use defutils
  use sparitersolve
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(model_t), pointer, save :: acc_model => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(variable_t), pointer, save :: acc_velovar => null()
  type(variable_t), pointer, save :: acc_divvar => null()
  type(variable_t), pointer, save :: acc_presvar => null()
  type(valuelist_t), pointer, save :: acc_bc => null()
  integer, allocatable, save :: acc_indexes(:)
  integer, allocatable, save :: acc_parentidx(:)
  ! LocalMatrix context (mirrors its SAVE Nodes and IP local).
  type(nodes_t), target, save :: acc_nodes_lm
  type(gaussintegrationpoints_t), save :: acc_ip_lm
  ! Traction context (mirrors the entry SAVE arrays used there).
  real(kind=dp), allocatable, save :: acc_cvelo(:,:), acc_ppres(:,:)
  real(kind=dp), allocatable, save :: acc_pvelo_wrk(:,:)
  real(kind=dp), allocatable, save :: acc_totalforce(:,:), acc_totalarea(:), acc_mu(:)
  type(nodes_t), target, save :: acc_nodes_tr, acc_parentnodes_tr
  character(len=512), save :: acc_message

contains

  ! -- entry context ---------------------------------------------------
  subroutine csp_begin(h_model, h_solver) bind(c, name="csp_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine csp_begin

  function csp_parallel() result(r) bind(c, name="csp_parallel")
    integer(c_int) :: r
    r = merge(1, 0, parenv % pes > 1)
  end function csp_parallel

  subroutine csp_block_precond(val, found) bind(c, name="csp_block_precond")
    integer(c_int) :: val, found
    type(valuelist_t), pointer :: sp
    logical :: v, f
    sp => GetSolverParams()
    v = GetLogical(sp, 'Block Preconditioning', f)
    if (.not. f) v = .true.
    val = merge(1, 0, v)
    found = merge(1, 0, f)
  end subroutine csp_block_precond

  function csp_dim() result(d) bind(c, name="csp_dim")
    integer(c_int) :: d
    d = CoordinateSystemDimension()
  end function csp_dim

  function csp_mesh_max_dofs() result(n) bind(c, name="csp_mesh_max_dofs")
    integer(c_int) :: n
    n = acc_solver % mesh % maxelementdofs
  end function csp_mesh_max_dofs

  function csp_matrix_rows() result(n) bind(c, name="csp_matrix_rows")
    integer(c_int) :: n
    n = acc_solver % matrix % numberofrows
  end function csp_matrix_rows

  function csp_num_bcs() result(n) bind(c, name="csp_num_bcs")
    integer(c_int) :: n
    n = acc_model % numberofbcs
  end function csp_num_bcs

  ! -- variable binding (with the source's Fatal checks) ----------------
  subroutine csp_bind_vars(block, ok_velo, ok_div, ok_pres) &
      bind(c, name="csp_bind_vars")
    integer(c_int), value :: block
    integer(c_int) :: ok_velo, ok_div, ok_pres
    type(mesh_t), pointer :: m
    m => GetMesh()
    acc_velovar => VariableGet(m % variables, "VelocityTot")
    acc_divvar => VariableGet(m % variables, "Divergence")
    ok_velo = merge(1, 0, associated(acc_velovar))
    ok_div = merge(1, 0, associated(acc_divvar))
    ok_pres = 1
    if (block /= 0) then
      acc_presvar => VariableGet(m % variables, "Flow")
      ok_pres = merge(1, 0, associated(acc_presvar))
    end if
  end subroutine csp_bind_vars

  subroutine csp_fatal_flow() bind(c, name="csp_fatal_flow")
    call Fatal('PressureSolver', 'The coupled flow variable Flow was not found')
  end subroutine csp_fatal_flow

  subroutine csp_fatal_velo() bind(c, name="csp_fatal_velo")
    call Fatal('PressureUpdate', 'Velocity Variable VelocityTot was not found')
  end subroutine csp_fatal_velo

  subroutine csp_fatal_div() bind(c, name="csp_fatal_div")
    call Fatal('PressureUpdate', 'Divergence Variable was not found')
  end subroutine csp_fatal_div

  subroutine csp_alloc(n, nbc) bind(c, name="csp_alloc")
    integer(c_int), value :: n, nbc
    if (allocated(acc_indexes)) deallocate(acc_indexes)
    if (allocated(acc_parentidx)) deallocate(acc_parentidx)
    if (allocated(acc_cvelo)) deallocate(acc_cvelo)
    if (allocated(acc_ppres)) deallocate(acc_ppres)
    if (allocated(acc_pvelo_wrk)) deallocate(acc_pvelo_wrk)
    if (allocated(acc_totalforce)) deallocate(acc_totalforce)
    if (allocated(acc_totalarea)) deallocate(acc_totalarea)
    if (allocated(acc_mu)) deallocate(acc_mu)
    allocate(acc_indexes(n), acc_parentidx(n))
    allocate(acc_cvelo(3, n), acc_ppres(n, 2), acc_pvelo_wrk(3, n))
    allocate(acc_totalforce(max(nbc, 1), 3), acc_totalarea(max(nbc, 1)), acc_mu(n))
  end subroutine csp_alloc

  subroutine csp_default_init() bind(c, name="csp_default_init")
    call DefaultInitialize()
  end subroutine csp_default_init

  ! -- active-element driver -------------------------------------------
  function csp_nof_active() result(n) bind(c, name="csp_nof_active")
    integer(c_int) :: n
    n = GetNOFActive()
  end function csp_nof_active

  subroutine csp_begin_active(t, n, nd, nb) bind(c, name="csp_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, nb
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes()
    nd = GetElementDOFs(acc_indexes)
    nb = GetElementNOFBDOFs()
  end subroutine csp_begin_active

  subroutine csp_fill_block(ntot, nd, cvelo, pvelo, ppres) bind(c, name="csp_fill_block")
    integer(c_int), value :: ntot, nd
    real(c_double) :: cvelo(*), pvelo(*), ppres(*)
    ! Block-preconditioning gather: previous velocity is divergence-free
    ! (PVelo zeroed), pressure from DivVar/Flow.
    integer :: i, d
    d = CoordinateSystemDimension()
    do i = 1, d
      call csp_set_col(cvelo, 3, i, nd, acc_velovar, 0)
      call csp_zero_col(pvelo, 3, i, nd)
    end do
    call csp_copy_div(ppres, ntot, nd, 1)
    call csp_copy_flow(ppres, ntot, nd, 2)
  end subroutine csp_fill_block

  subroutine csp_set_col(dst, stride, col, nd, var, prev)
    real(c_double) :: dst(stride, *)
    integer(c_int), value :: stride, col, nd, prev
    type(variable_t), pointer :: var
    integer :: k
    if (prev == 0) then
      do k = 1, nd
        dst(col, k) = var % values(var % dofs * (var % perm(acc_indexes(k)) - 1) + col)
      end do
    else
      do k = 1, nd
        dst(col, k) = var % prevvalues(var % dofs * (var % perm(acc_indexes(k)) - 1) + col, 1)
      end do
    end if
  end subroutine csp_set_col

  subroutine csp_zero_col(dst, stride, col, nd)
    real(c_double) :: dst(stride, *)
    integer(c_int), value :: stride, col, nd
    integer :: k
    do k = 1, nd
      dst(col, k) = 0.0d0
    end do
  end subroutine csp_zero_col

  subroutine csp_copy_div(dst, ntot, nd, col)
    real(c_double) :: dst(*)
    integer(c_int), value :: ntot, nd, col
    integer :: k
    do k = 1, nd
      dst((col - 1) * ntot + k) = acc_divvar % values(acc_divvar % perm(acc_indexes(k)))
    end do
  end subroutine csp_copy_div

  subroutine csp_copy_flow(dst, ntot, nd, col)
    real(c_double) :: dst(*)
    integer(c_int), value :: ntot, nd, col
    integer :: k, dk
    dk = CoordinateSystemDimension() + 1
    do k = 1, nd
      dst((col - 1) * ntot + k) = acc_presvar % values( &
          acc_presvar % dofs * (acc_presvar % perm(acc_indexes(k)) - 1) + dk)
    end do
  end subroutine csp_copy_flow

  subroutine csp_fill_std(ntot, nd, cvelo, pvelo, ppres) bind(c, name="csp_fill_std")
    integer(c_int), value :: ntot, nd
    real(c_double) :: cvelo(*), pvelo(*), ppres(*)
    integer :: i, k, d
    d = CoordinateSystemDimension()
    do i = 1, d
      call csp_set_col(cvelo, 3, i, nd, acc_velovar, 0)
      call csp_set_col(pvelo, 3, i, nd, acc_velovar, 1)
    end do
    call csp_copy_div(ppres, ntot, nd, 1)
    do k = 1, nd
      ppres(ntot + k) = acc_solver % variable % prevvalues( &
          acc_solver % variable % perm(acc_indexes(k)), 1)
    end do
  end subroutine csp_fill_std

  ! -- LocalMatrix context ----------------------------------------------
  subroutine csp_lm_nodes() bind(c, name="csp_lm_nodes")
    call GetElementNodes(acc_nodes_lm)
  end subroutine csp_lm_nodes

  function csp_lm_ip_begin() result(n) bind(c, name="csp_lm_ip_begin")
    integer(c_int) :: n
    acc_ip_lm = GaussPoints(acc_elem)
    n = acc_ip_lm % n
  end function csp_lm_ip_begin

  subroutine csp_lm_ip_point(t, u, v, w, s) bind(c, name="csp_lm_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip_lm % u(t)
    v = acc_ip_lm % v(t)
    w = acc_ip_lm % w(t)
    s = acc_ip_lm % s(t)
  end subroutine csp_lm_ip_point

  subroutine csp_lm_eleminfo(u, v, w, n, nd, basis, dbasis, detj, stat) &
      bind(c, name="csp_lm_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: n, nd
    real(c_double) :: basis(n), dbasis(nd, 3), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_elem, acc_nodes_lm, u, v, w, detj, &
        basis(1:n), dbasis(1:nd, 1:3))
    stat = merge(1, 0, ok)
  end subroutine csp_lm_eleminfo

  subroutine csp_lm_point(n, nd, ntot, dt, sdet, ppres, cvelo, pvelo, &
      basis, dbasis, stiff, force) bind(c, name="csp_lm_point")
    integer(c_int), value :: n, nd, ntot
    real(c_double), value :: dt, sdet
    real(c_double) :: ppres(ntot, 2), cvelo(3, ntot), pvelo(3, ntot)
    real(c_double) :: basis(ntot), dbasis(ntot, 3)
    real(c_double) :: stiff(ntot, ntot), force(ntot)
    ! Exact transcription of the per-IP body (reductions and the p/q
    ! assembly). The MATMUL pressure correction and the bubble/boundary
    ! elimination run once per element in csp_lm_correct, exactly as the
    ! source places them after the integration loop. Same compiler, same
    ! order as the original internal procedure.
    real(c_double) :: ddivu
    integer :: d, i, p, q
    d = CoordinateSystemDimension()
    ddivu = 0.0d0
    do i = 1, d
      ddivu = ddivu + (sum(cvelo(i, 1:nd) * dbasis(1:nd, i)) - &
          sum(pvelo(i, 1:nd) * dbasis(1:nd, i))) / dt
    end do
    do p = 1, n
      do q = 1, n
        stiff(p, q) = stiff(p, q) + sdet * sum(dbasis(q, 1:d) * dbasis(p, 1:d))
      end do
      force(p) = force(p) - sdet * ddivu * basis(p)
    end do
  end subroutine csp_lm_point

  subroutine csp_lm_correct(n, nd, ntot, ppres, stiff, force) &
      bind(c, name="csp_lm_correct")
    integer(c_int), value :: n, nd, ntot
    real(c_double) :: ppres(ntot, 2)
    real(c_double) :: stiff(ntot, ntot), force(ntot)
    ! Once-per-element tail: MATMUL pressure correction, then bubble
    ! elimination. Must run once (not per integration point).
    integer :: i
    force(1:n) = force(1:n) - matmul(stiff(1:n, 1:n), ppres(1:n, 1)) + &
        matmul(stiff(1:n, 1:n), ppres(1:n, 2))
    do i = n + 1, ntot
      force(i) = 0.0d0
      stiff(i, :) = 0.0d0
      stiff(:, i) = 0.0d0
      stiff(i, i) = 1.0d0
    end do
  end subroutine csp_lm_correct

  subroutine csp_lm_load_zero(nd, load) bind(c, name="csp_lm_load_zero")
    integer(c_int), value :: nd
    real(c_double) :: load(*)
    load(1:nd) = 0.0d0
  end subroutine csp_lm_load_zero

  subroutine csp_update_eq(ld, g, f) bind(c, name="csp_update_eq")
    integer(c_int), value :: ld
    real(c_double) :: g(ld, ld), f(ld)
    call DefaultUpdateEquations(g(1:ld, 1:ld), f(1:ld))
  end subroutine csp_update_eq

  subroutine csp_finish() bind(c, name="csp_finish")
    call DefaultFinishAssembly()
  end subroutine csp_finish

  function csp_solve() result(norm) bind(c, name="csp_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function csp_solve


  ! -- outflow-boundary matrix surgery ----------------------------------
  function csp_nof_boundary() result(n) bind(c, name="csp_nof_boundary")
    integer(c_int) :: n
    n = GetNOFBoundaryElements()
  end function csp_nof_boundary

  subroutine csp_begin_boundary(t, n, nd, active, family) &
      bind(c, name="csp_begin_boundary")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, active, family
    acc_elem => GetBoundaryElement(t)
    active = merge(1, 0, ActiveBoundaryElement())
    n = GetElementNOFNodes()
    nd = GetElementDOFs(acc_indexes)
    family = GetElementFamily()
  end subroutine csp_begin_boundary

  subroutine csp_bc() bind(c, name="csp_bc")
    acc_bc => GetBC()
  end subroutine csp_bc

  function csp_bc_associated() result(r) bind(c, name="csp_bc_associated")
    integer(c_int) :: r
    r = merge(1, 0, associated(acc_bc))
  end function csp_bc_associated

  subroutine csp_bc_outflow(val, found) bind(c, name="csp_bc_outflow")
    integer(c_int) :: val, found
    logical :: v, f
    v = ListGetLogical(acc_bc, 'Outflow Boundary', f)
    val = merge(1, 0, v)
    found = merge(1, 0, f)
  end subroutine csp_bc_outflow

  subroutine csp_outflow_surgery(block, nd, n) bind(c, name="csp_outflow_surgery")
    integer(c_int), value :: block, nd, n
    ! Exact transcription of the outflow-boundary body: velocity gathers
    ! plus the per-node ZeroRow/SetMatrixElement/RHS surgery.
    integer :: i, j, k, d
    real(c_double) :: ppres_j
    d = CoordinateSystemDimension()
    do i = 1, d
      call csp_set_col(acc_cvelo, 3, i, nd, acc_velovar, 0)
      call csp_set_col(acc_pvelo_wrk, 3, i, nd, acc_velovar, 1)
    end do
    do i = 1, nd
      j = acc_solver % variable % perm(acc_indexes(i))
      k = acc_divvar % perm(acc_indexes(i))
      call ZeroRow(acc_solver % matrix, j)
      call SetMatrixElement(acc_solver % matrix, j, j, 1.0d0)
      if (i > n) then
        acc_solver % matrix % rhs(j) = 0.0d0
      else
        if (block /= 0) then
          ppres_j = acc_presvar % values(acc_presvar % dofs * &
              (acc_presvar % perm(acc_indexes(i)) - 1) + d + 1) - &
              acc_divvar % values(k)
        else
          ppres_j = acc_solver % variable % prevvalues(j, 1) - acc_divvar % values(k)
        end if
        acc_solver % matrix % rhs(j) = ppres_j
      end if
    end do
  end subroutine csp_outflow_surgery


  ! -- traction post-pass (whole block, Fortran-side) -------------------
  subroutine csp_traction() bind(c, name="csp_traction")
    integer, pointer :: flowperm(:)
    integer :: t, i, j, k, n, nd, np, ndp, nlen, active, d
    logical :: found
    type(element_t), pointer :: parent, current
    type(valuelist_t), pointer :: material
    real(kind=dp) :: traction(3), area
    character(len=10) :: boundaryname
    d = CoordinateSystemDimension()
    flowperm => acc_solver % variable % perm
    acc_totalforce = 0.0d0
    acc_totalarea = 0.0d0
    active = GetNOFBoundaryElements()
    do t = 1, active
      acc_elem => GetBoundaryElement(t)
      if (.not. ActiveBoundaryElement()) cycle
      if (acc_elem % type % elementcode == 101) cycle
      if ((d > 2) .and. acc_elem % type % elementcode == 202) cycle
      n = GetElementNOFNodes()
      nd = GetElementDOFs(acc_indexes)
      do i = 1, acc_model % numberofbcs
        if (acc_elem % boundaryinfo % constraint == acc_model % bcs(i) % tag) then
          if (.not. ListGetLogical(acc_model % bcs(i) % values, &
              'Calculate Fluidic Force', found)) cycle
          parent => acc_elem % boundaryinfo % left
          found = associated(parent)
          if (found) found = all(flowperm(parent % nodeindexes(1:parent % type % numberofnodes)) > 0)
          if (.not. found) then
            parent => acc_elem % boundaryinfo % right
            found = associated(parent)
            if (found) found = all(flowperm(parent % nodeindexes(1:parent % type % numberofnodes)) > 0)
            if (.not. found) call Fatal('PressureUpdate', &
                'No parent element can be found for given boundary element')
          end if
          np = parent % type % numberofnodes
          ndp = GetElementDOFs(acc_parentidx, parent)
          k = ListGetInteger(acc_model % bodies(parent % bodyid) % values, 'Material')
          material => acc_model % materials(k) % values
          acc_mu(1:np) = ListGetReal(material, 'Viscosity', np, parent % nodeindexes(1:np))
          do j = 1, d
            acc_cvelo(j, 1:ndp) = acc_velovar % values( &
                acc_velovar % dofs * (acc_velovar % perm(acc_parentidx(1:ndp)) - 1) + j)
          end do
          acc_ppres(1:ndp, 1) = acc_solver % variable % values( &
              acc_solver % variable % perm(acc_parentidx(1:ndp)))
          traction = 0.0d0
          call csp_surface_force(acc_elem, parent, traction, area, &
              acc_cvelo, acc_ppres, acc_mu, np, ndp)
          acc_totalforce(i, 1:3) = acc_totalforce(i, 1:3) + traction(1:3)
          acc_totalarea(i) = acc_totalarea(i) + area
        end if
      end do
    end do
    do k = 1, acc_model % numberofbcs
      if (.not. ListGetLogical(acc_model % bcs(k) % values, &
          'Calculate Fluidic Force', found)) cycle
      if (acc_model % numberofbcs < 10) then
        write (boundaryname, '("bc ",I1)') k
      else
        write (boundaryname, '("bc ",I2)') k
      end if
      nlen = len_trim(boundaryname)
      call Info('ForceCompute', 'Forces on Boundary ' // boundaryname(1:nlen), level=4)
      write (acc_message, '("Fluidic Force (X,Y,Z):", 3ES17.6E2)') acc_totalforce(k, 1:3)
      call Info('ForceCompute', acc_message, level=4)
      write (acc_message, '("Contact Area:   ", ES17.6E2)') acc_totalarea(k)
      call Info('ForceCompute', acc_message, level=4)
      call ListAddConstReal(acc_model % simulation, &
          'res: contact force 1 ' // boundaryname(1:nlen), acc_totalforce(k, 1))
      call ListAddConstReal(acc_model % simulation, &
          'res: contact force 2 ' // boundaryname(1:nlen), acc_totalforce(k, 2))
      if (d > 2) call ListAddConstReal(acc_model % simulation, &
          'res: contact force 3 ' // boundaryname(1:nlen), acc_totalforce(k, 3))
      call ListAddConstReal(acc_model % simulation, &
          'res: contact force area ' // boundaryname(1:nlen), acc_totalarea(k))
    end do
  end subroutine csp_traction

  subroutine csp_surface_force(element, parent, traction, area, &
      velo, pres, viscosity, np, ndp)
    type(element_t), pointer :: element, parent
    real(kind=dp) :: traction(3), area, velo(:,:), pres(:,:), viscosity(:)
    integer :: np, ndp
    ! Verbatim transcription of SurfaceForceIntegration (dead arms kept
    ! structurally: the SELECT CASE default leaves Normal zero-filled as
    ! in the source — Normal is a fresh local, zero-initialized below).
    logical :: stat
    integer :: n_integ, t, i, d
    real(kind=dp) :: u, v, w, detj, s, basis(ndp), dbasis(ndp, 3), normal(3), &
        visc, rep, tmpmat(3, ndp), regrad(3, 3), red(3, 3)
    real(kind=dp), pointer :: u_integ(:), v_integ(:), w_integ(:), s_integ(:)
    type(gaussintegrationpoints_t), target :: integstuff
    d = CoordinateSystemDimension()
    dbasis = 0.0d0
    traction = 0.0d0
    area = 0.0d0
    normal = 0.0d0
    call GetElementNodes(acc_nodes_tr, element)
    select case (element % type % numberofnodes)
    case (2)
      normal = Normalvector(element, acc_nodes_tr, 0.0d0, 0.0d0, .true.)
    case (3)
      normal = Normalvector(element, acc_nodes_tr, 0.3d0, 0.3d0, .true.)
    case (4)
      normal = Normalvector(element, acc_nodes_tr, 0.0d0, 0.0d0, .true.)
    end select
    call GetElementNodes(acc_parentnodes_tr, parent)
    integstuff = Gausspoints(parent)
    u_integ => integstuff % u
    v_integ => integstuff % v
    w_integ => integstuff % w
    s_integ => integstuff % s
    n_integ = integstuff % n
    do t = 1, n_integ
      u = u_integ(t)
      v = v_integ(t)
      w = w_integ(t)
      stat = ElementInfo(parent, acc_parentnodes_tr, u, v, w, &
          detj, basis, dbasis)
      s = detj * s_integ(t)
      visc = sum(viscosity(1:np) * basis(1:np))
      rep = sum(pres(1:ndp, 1) * basis(1:ndp))
      tmpmat(1, 1:ndp) = velo(1, 1:ndp)
      tmpmat(2, 1:ndp) = velo(2, 1:ndp)
      if (d > 2) then
        tmpmat(3, 1:ndp) = velo(3, 1:ndp)
      else
        tmpmat(3, 1:ndp) = 0.0d0
      end if
      regrad = matmul(tmpmat, dbasis)
      red = visc * (regrad + transpose(regrad))
      traction(1:3) = traction(1:3) - rep * normal(1:3) + matmul(red, normal)
    end do
    area = ElementArea(acc_solver % mesh, element, element % type % numberofnodes)
    traction = -area / n_integ * traction
  end subroutine csp_surface_force

  subroutine csp_fatal_alloc() bind(c, name="csp_fatal_alloc")
    call Fatal('PressureUpdate', 'Memory allocation error.')
  end subroutine csp_fatal_alloc

end module cs_pres_acc
