! CutFEM advection-diffusion solver plugin: statement-level BIND(C)
! wrappers for the Elmer calls in fem/src/modules/ModelPDECutFEM.F90,
! compiled against the built .mod files of the target configure.
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module mirrors the original's locals (host-associated MatVol, current
! element handles, Nodes SAVE state per internal procedure, IP shallow
! copies) with identical lifetimes.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module modelpdecutfem_acc
  use iso_c_binding
  use types
  use defutils
  use integration
  use matrixassembly
  use eleminfo
  use messages
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(element_t), pointer, save :: acc_pelem => null()
  type(element_t), pointer, save :: acc_lmelem => null()
  type(valuelist_t), pointer, save :: acc_bf => null()
  type(valuelist_t), pointer, save :: acc_mat => null()
  type(valuelist_t), pointer, save :: acc_bc => null()
  ! Host-associated MatVol(0:10), zeroed on every solver entry.
  real(kind=dp), save :: acc_matvol(0:10)
  ! LocalMatrix context (mirrors its SAVE Nodes and IP local).
  type(nodes_t), target, save :: acc_nodes_lm
  type(gaussintegrationpoints_t), save :: acc_ip_lm
  ! LocalMatrixBC context (mirrors its SAVE Nodes and IP local).
  type(nodes_t), target, save :: acc_nodes_bc
  type(gaussintegrationpoints_t), save :: acc_ip_bc


contains
  function mcf_max_elem_dofs() result(n) bind(c, name="mcf_max_elem_dofs")
    integer(c_int) :: n
    n = acc_solver % mesh % maxelementdofs
  end function mcf_max_elem_dofs

  ! -- entry state ----------------------------------------------------
  subroutine mcf_begin(h_model, h_solver) bind(c, name="mcf_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine mcf_begin

  function mcf_cutfem_flag() result(r) bind(c, name="mcf_cutfem_flag")
    integer(c_int) :: r
    logical :: found
    r = merge(1, 0, ListGetLogical(acc_solver % values, 'CutFEM', found))
  end function mcf_cutfem_flag

  subroutine mcf_init_add_flags() bind(c, name="mcf_init_add_flags")
    call ListAddNewLogical(acc_solver % values, 'Use Global Mass Matrix', .true.)
    call ListAddNewLogical(acc_solver % values, 'CutFEM Solver', .true.)
  end subroutine mcf_init_add_flags

  function mcf_maxiter(found) result(n) bind(c, name="mcf_maxiter")
    integer(c_int) :: n, found
    logical :: f
    n = ListGetInteger(GetSolverParams(), &
        'Nonlinear System Max Iterations', f, minv=1)
    found = merge(1, 0, f)
  end function mcf_maxiter

  subroutine mcf_default_start() bind(c, name="mcf_default_start")
    call DefaultStart()
  end subroutine mcf_default_start

  subroutine mcf_default_finish() bind(c, name="mcf_default_finish")
    call DefaultFinish()
  end subroutine mcf_default_finish

  function mcf_get_nof_active() result(n) bind(c, name="mcf_get_nof_active")
    integer(c_int) :: n
    n = GetNOFActive(acc_solver)
  end function mcf_get_nof_active

  subroutine mcf_default_finish_bulk() bind(c, name="mcf_default_finish_bulk")
    call DefaultFinishBulkAssembly()
  end subroutine mcf_default_finish_bulk

  function mcf_get_nof_boundary() result(n) bind(c, name="mcf_get_nof_boundary")
    integer(c_int) :: n
    n = GetNOFBoundaryElements()
  end function mcf_get_nof_boundary

  subroutine mcf_default_finish_boundary() bind(c, name="mcf_default_finish_boundary")
    call DefaultFinishBoundaryAssembly()
  end subroutine mcf_default_finish_boundary

  subroutine mcf_default_finish_assembly() bind(c, name="mcf_default_finish_assembly")
    call DefaultFinishAssembly()
  end subroutine mcf_default_finish_assembly

  subroutine mcf_default_dirichlet() bind(c, name="mcf_default_dirichlet")
    call DefaultDirichletBCs()
  end subroutine mcf_default_dirichlet

  function mcf_default_solve() result(norm) bind(c, name="mcf_default_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function mcf_default_solve

  function mcf_default_converged() result(r) bind(c, name="mcf_default_converged")
    integer(c_int) :: r
    r = merge(1, 0, DefaultConverged())
  end function mcf_default_converged

  subroutine mcf_set_norm(norm) bind(c, name="mcf_set_norm")
    real(c_double), value :: norm
    acc_solver % variable % norm = norm
  end subroutine mcf_set_norm

  subroutine mcf_fatal_alloc() bind(c, name="mcf_fatal_alloc")
    call Fatal('AdvDiffSolver', 'Memory allocation error.')
  end subroutine mcf_fatal_alloc

  ! -- host-associated MatVol ------------------------------------------
  subroutine mcf_matvol_zero() bind(c, name="mcf_matvol_zero")
    acc_matvol = 0.0_dp
  end subroutine mcf_matvol_zero

  subroutine mcf_matvol_report() bind(c, name="mcf_matvol_report")
    integer :: n
    if (InfoActive(10)) then
      print *, 'CutFEM Interface area:', acc_matvol(0)
      do n = 1, 10
        if (acc_matvol(n) > epsilon(acc_matvol(n))) &
            print *, 'CutFEM body '//i2s(n)//' volume: ', acc_matvol(n)
      end do
    end if
  end subroutine mcf_matvol_report

  ! -- bulk active-element context --------------------------------------
  subroutine mcf_begin_active(t, n, nd, nb) bind(c, name="mcf_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, nb
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes(acc_elem)
    nd = GetElementNOFDOFs(acc_elem)
    nb = GetElementNOFBDOFs(acc_elem)
  end subroutine mcf_begin_active

  subroutine mcf_cut_check(iscut, isactive) bind(c, name="mcf_cut_check")
    integer(c_int) :: iscut, isactive
    logical :: c, a
    call CutinterfaceCheck(acc_elem, c, a)
    iscut = merge(1, 0, c)
    isactive = merge(1, 0, a)
  end subroutine mcf_cut_check

  subroutine mcf_cut_bulk(iscut, ismore) bind(c, name="mcf_cut_bulk")
    integer(c_int) :: iscut, ismore
    logical :: c, m
    acc_pelem => CutInterfaceBulk(acc_elem, c, m)
    iscut = merge(1, 0, c)
    ismore = merge(1, 0, m)
  end subroutine mcf_cut_bulk

  subroutine mcf_pelem_dims(n, nd, nb) bind(c, name="mcf_pelem_dims")
    integer(c_int) :: n, nd, nb
    n = GetElementNOFNodes(acc_pelem)
    nd = GetElementNOFDOFs(acc_pelem)
    nb = GetElementNOFBDOFs(acc_pelem)
  end subroutine mcf_pelem_dims

  function mcf_pelem_perm_ok() result(r) bind(c, name="mcf_pelem_perm_ok")
    integer(c_int) :: r
    r = merge(1, 0, &
        all(acc_solver % variable % perm(acc_pelem % nodeindexes) > 0))
  end function mcf_pelem_perm_ok

  subroutine mcf_begin_cutbc_active(t) bind(c, name="mcf_begin_cutbc_active")
    integer(c_int), value :: t
    acc_elem => GetActiveElement(t)
  end subroutine mcf_begin_cutbc_active

  subroutine mcf_cut_bc(assoc, iscut, ismore) bind(c, name="mcf_cut_bc")
    integer(c_int) :: assoc, iscut, ismore
    logical :: c, m
    acc_pelem => CutInterfaceBC(acc_elem, c, m)
    assoc = merge(1, 0, associated(acc_pelem))
    iscut = merge(1, 0, c)
    ismore = merge(1, 0, m)
  end subroutine mcf_cut_bc

  function mcf_pelem_typenodes() result(n) bind(c, name="mcf_pelem_typenodes")
    integer(c_int) :: n
    n = acc_pelem % type % numberofnodes
  end function mcf_pelem_typenodes

  subroutine mcf_set_current_to_pelem() bind(c, name="mcf_set_current_to_pelem")
    acc_model % currentelement => acc_pelem
  end subroutine mcf_set_current_to_pelem

  ! -- boundary-element context ------------------------------------------
  subroutine mcf_begin_boundary(t) bind(c, name="mcf_begin_boundary")
    integer(c_int), value :: t
    acc_elem => GetBoundaryElement(t)
  end subroutine mcf_begin_boundary

  function mcf_bnd_perm_has_zero() result(r) bind(c, name="mcf_bnd_perm_has_zero")
    integer(c_int) :: r
    r = merge(1, 0, &
        any(acc_solver % variable % perm(acc_elem % nodeindexes) == 0))
  end function mcf_bnd_perm_has_zero

  function mcf_bnd_active(n, nd) result(r) bind(c, name="mcf_bnd_active")
    integer(c_int) :: n, nd, r
    r = merge(1, 0, ActiveBoundaryElement(acc_elem))
    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
  end function mcf_bnd_active

  ! -- LocalMatrix context -----------------------------------------------
  subroutine mcf_lm_select(cut) bind(c, name="mcf_lm_select")
    integer(c_int), value :: cut
    if (cut /= 0) then
      acc_lmelem => acc_pelem
    else
      acc_lmelem => acc_elem
    end if
  end subroutine mcf_lm_select

  subroutine mcf_lm_nodes() bind(c, name="mcf_lm_nodes")
    call GetElementNodes(acc_nodes_lm, acc_lmelem)
  end subroutine mcf_lm_nodes

  function mcf_coord_dim() result(d) bind(c, name="mcf_coord_dim")
    integer(c_int) :: d
    d = CoordinateSystemDimension()
  end function mcf_coord_dim
  subroutine mcf_lm_assemble(nd, dim, ld, weight, d, c, r, rho, loadatip, &
      avec, basis, dbasis, stiff, mass, force, docut) &
      bind(c, name="mcf_lm_assemble")
    integer(c_int), value :: nd, dim, ld, docut
    real(c_double) :: stiff(ld, ld), mass(ld, ld), force(ld)
    real(c_double) :: basis(nd), dbasis(nd, 3), avec(3)
    real(c_double), value :: weight, d, c, r, rho, loadatip
    integer :: p, q
    ! Exact transcription of the update statements (same compiler, same
    ! evaluation order as the original internal procedure); the state
    ! buffers carry leading dimension ld, only the leading nd section is
    ! touched, exactly as the original's (nd,nd) automatics.
    stiff(1:nd, 1:nd) = stiff(1:nd, 1:nd) + weight * &
        d * matmul(dbasis, transpose(dbasis))
    do p = 1, nd
      do q = 1, nd
        stiff(p, q) = stiff(p, q) + weight * &
            c * sum(avec(1:dim) * dbasis(q, 1:dim)) * basis(p)
        stiff(p, q) = stiff(p, q) + weight * r * basis(q) * basis(p)
        mass(p, q) = mass(p, q) + weight * rho * basis(q) * basis(p)
      end do
    end do
    force(1:nd) = force(1:nd) + weight * loadatip * basis(1:nd)
    if (docut /= 0) &
        acc_matvol(acc_lmelem % bodyid) = acc_matvol(acc_lmelem % bodyid) + weight
  end subroutine mcf_lm_assemble
  function mcf_lm_bodyforce() result(r) bind(c, name="mcf_lm_bodyforce")
    integer(c_int) :: r
    acc_bf => GetBodyForce(acc_lmelem)
    r = merge(1, 0, associated(acc_bf))
  end function mcf_lm_bodyforce

  subroutine mcf_lm_load(n, load, found) bind(c, name="mcf_lm_load")
    integer(c_int), value :: n
    real(c_double) :: load(n)
    integer(c_int) :: found
    logical :: f
    load(1:n) = GetReal(acc_bf, 'field source', f, acc_lmelem)
    found = merge(1, 0, f)
  end subroutine mcf_lm_load

  subroutine mcf_lm_mat(n, diff, conv, react, timec) bind(c, name="mcf_lm_mat")
    integer(c_int), value :: n
    real(c_double) :: diff(n), conv(n), react(n), timec(n)
    logical :: f
    acc_mat => GetMaterial(acc_lmelem)
    diff(1:n) = GetReal(acc_mat, 'diffusion coefficient', f, acc_lmelem)
    react(1:n) = GetReal(acc_mat, 'reaction coefficient', f, acc_lmelem)
    conv(1:n) = GetReal(acc_mat, 'convection coefficient', f, acc_lmelem)
    timec(1:n) = GetReal(acc_mat, 'time derivative coefficient', f, acc_lmelem)
  end subroutine mcf_lm_mat

  subroutine mcf_lm_velo(i, n, row) bind(c, name="mcf_lm_velo")
    integer(c_int), value :: i, n
    real(c_double) :: row(n)
    logical :: f
    row(1:n) = GetReal(acc_mat, 'convection velocity '//i2s(i), f, acc_lmelem)
  end subroutine mcf_lm_velo

  function mcf_lm_ip_begin() result(n) bind(c, name="mcf_lm_ip_begin")
    integer(c_int) :: n
    acc_ip_lm = GaussPointsAdapt(acc_lmelem)
    if (acc_lmelem % elementindex == 1) then
      call Info('AdvDiffSolver', &
          'Integration points in 1st element: '//i2s(acc_ip_lm % n), level=8)
    end if
    n = acc_ip_lm % n
  end function mcf_lm_ip_begin

  subroutine mcf_lm_ip_point(t, u, v, w, s) bind(c, name="mcf_lm_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip_lm % u(t)
    v = acc_ip_lm % v(t)
    w = acc_ip_lm % w(t)
    s = acc_ip_lm % s(t)
  end subroutine mcf_lm_ip_point

  subroutine mcf_lm_eleminfo(u, v, w, n, nd, basis, dbasis, detj, stat) &
      bind(c, name="mcf_lm_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: n, nd
    real(c_double) :: basis(n), dbasis(nd, 3), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_lmelem, acc_nodes_lm, u, v, w, detj, &
        basis(1:n), dbasis(1:nd, 1:3))
    stat = merge(1, 0, ok)
  end subroutine mcf_lm_eleminfo

  function mcf_dot(n, a, b) result(s) bind(c, name="mcf_dot")
    integer(c_int), value :: n
    real(c_double) :: a(n), b(n), s
    s = sum(a(1:n) * b(1:n))
  end function mcf_dot

  subroutine mcf_lm_amatvec(n, velo, basis, a) bind(c, name="mcf_lm_amatvec")
    integer(c_int), value :: n
    real(c_double) :: velo(3, n), basis(n), a(3)
    a = matmul(velo(:, 1:n), basis(1:n))
  end subroutine mcf_lm_amatvec


  subroutine mcf_lm_first_order_time(ld, mass, stiff, force) &
      bind(c, name="mcf_lm_first_order_time")
    integer(c_int), value :: ld
    real(c_double) :: mass(ld, ld), stiff(ld, ld), force(ld)
    call Default1stOrderTime(mass, stiff, force, acc_lmelem)
  end subroutine mcf_lm_first_order_time

  subroutine mcf_lm_update_mass(n, nd, ld, mass) bind(c, name="mcf_lm_update_mass")
    integer(c_int), value :: n, nd, ld
    real(c_double) :: mass(ld, ld)
    type(matrix_t), pointer :: a
    type(variable_t), pointer :: x
    a => acc_solver % matrix
    x => acc_solver % variable
    call UpdateMassMatrix(a, mass(1:nd, 1:nd), nd, x % dofs, &
        x % perm(acc_lmelem % nodeindexes(1:n)), a % massvalues)
  end subroutine mcf_lm_update_mass

  subroutine mcf_condensate(ndmnb, nb, ld, k, f) bind(c, name="mcf_condensate")
    integer(c_int), value :: ndmnb, nb, ld
    real(c_double) :: k(ld, ld), f(ld)
    call CondensateP(ndmnb, nb, k(1:ld, 1:ld), f(1:ld))
  end subroutine mcf_condensate

  subroutine mcf_lm_update_global(nd, ld, stiff, force) &
      bind(c, name="mcf_lm_update_global")
    integer(c_int), value :: nd, ld
    real(c_double) :: stiff(ld, ld), force(ld)
    type(matrix_t), pointer :: a
    type(variable_t), pointer :: x
    a => acc_solver % matrix
    x => acc_solver % variable
    call UpdateGlobalEquations(a, stiff(1:nd, 1:nd), a % rhs, force(1:nd), &
        nd, x % dofs, x % perm(acc_lmelem % nodeindexes), uelement=acc_lmelem)
  end subroutine mcf_lm_update_global

  subroutine mcf_lm_default_update(nd, ld, stiff, force) &
      bind(c, name="mcf_lm_default_update")
    integer(c_int), value :: nd, ld
    real(c_double) :: stiff(ld, ld), force(ld)
    call DefaultUpdateEquations(stiff(1:nd, 1:nd), force(1:nd), acc_lmelem)
  end subroutine mcf_lm_default_update

  ! -- LocalMatrixBC context ----------------------------------------------
  subroutine mcf_lbc_select(cut) bind(c, name="mcf_lbc_select")
    integer(c_int), value :: cut
    if (cut /= 0) then
      acc_lmelem => acc_pelem
      acc_model % currentelement => acc_pelem
    else
      acc_lmelem => acc_elem
    end if
  end subroutine mcf_lbc_select

  function mcf_lbc_begin() result(r) bind(c, name="mcf_lbc_begin")
    integer(c_int) :: r
    acc_bc => GetBC(acc_lmelem)
    r = merge(1, 0, associated(acc_bc))
  end function mcf_lbc_begin

  subroutine mcf_lbc_nodes() bind(c, name="mcf_lbc_nodes")
    call GetElementNodes(acc_nodes_bc, acc_lmelem)
  end subroutine mcf_lbc_nodes

  subroutine mcf_lbc_vals(n, flux, coeff, ext) bind(c, name="mcf_lbc_vals")
    integer(c_int), value :: n
    real(c_double) :: flux(n), coeff(n), ext(n)
    logical :: f
    flux(1:n) = GetReal(acc_bc, 'field flux', f, acc_lmelem)
    coeff(1:n) = GetReal(acc_bc, 'robin coefficient', f, acc_lmelem)
    ext(1:n) = GetReal(acc_bc, 'external field', f, acc_lmelem)
  end subroutine mcf_lbc_vals

  function mcf_lbc_ip_begin() result(n) bind(c, name="mcf_lbc_ip_begin")
    integer(c_int) :: n
    acc_ip_bc = GaussPoints(acc_lmelem)
    n = acc_ip_bc % n
  end function mcf_lbc_ip_begin

  subroutine mcf_lbc_ip_point(t, u, v, w, s) bind(c, name="mcf_lbc_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip_bc % u(t)
    v = acc_ip_bc % v(t)
    w = acc_ip_bc % w(t)
    s = acc_ip_bc % s(t)
  end subroutine mcf_lbc_ip_point

  subroutine mcf_lbc_eleminfo(u, v, w, n, nd, basis, dbasis, detj, stat) &
      bind(c, name="mcf_lbc_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: n, nd
    real(c_double) :: basis(n), dbasis(nd, 3), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_lmelem, acc_nodes_bc, u, v, w, detj, &
        basis(1:n), dbasis(1:nd, 1:3))
    stat = merge(1, 0, ok)
  end subroutine mcf_lbc_eleminfo

  subroutine mcf_lbc_assemble(nd, ld, weight, f, c, ext, basis, stiff, force) &
      bind(c, name="mcf_lbc_assemble")
    integer(c_int), value :: nd, ld
    real(c_double), value :: weight, f, c, ext
    real(c_double) :: basis(nd), stiff(ld, ld), force(ld)
    integer :: p, q
    ! Exact transcription of the update statements; the state buffers
    ! carry leading dimension ld, only the leading nd section is touched.
    do p = 1, nd
      do q = 1, nd
        stiff(p, q) = stiff(p, q) + weight * c * basis(q) * basis(p)
      end do
    end do
    force(1:nd) = force(1:nd) + weight * (f + c * ext) * basis(1:nd)
  end subroutine mcf_lbc_assemble

  subroutine mcf_lbc_update_global(nd, ld, stiff, force, lastweight) &
      bind(c, name="mcf_lbc_update_global")
    integer(c_int), value :: nd, ld
    real(c_double) :: stiff(ld, ld), force(ld)
    real(c_double), value :: lastweight
    type(matrix_t), pointer :: a
    type(variable_t), pointer :: x
    a => acc_solver % matrix
    x => acc_solver % variable
    call UpdateGlobalEquations(a, stiff(1:nd, 1:nd), a % rhs, force(1:nd), &
        nd, x % dofs, x % perm(acc_lmelem % nodeindexes), uelement=acc_lmelem)
    acc_matvol(0) = acc_matvol(0) + lastweight
  end subroutine mcf_lbc_update_global

  subroutine mcf_lbc_default_update(nd, ld, stiff, force) &
      bind(c, name="mcf_lbc_default_update")
    integer(c_int), value :: nd, ld
    real(c_double) :: stiff(ld, ld), force(ld)
    call DefaultUpdateEquations(stiff(1:nd, 1:nd), force(1:nd), acc_lmelem)
  end subroutine mcf_lbc_default_update

end module modelpdecutfem_acc
