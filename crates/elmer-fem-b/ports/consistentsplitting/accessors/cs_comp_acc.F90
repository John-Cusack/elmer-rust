! ConsistentSplitting pilot: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/ConsistentSplitting/CompressibilitySolver.F90,
! compiled against the built .mod files of the target configure.
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths — stays inside Fortran, so the
! Rust side is pure C calls. The wrapper module mirrors the entry's own
! SAVE state (current element/solver/model, VeloVar handle, Indexes store,
! LocalMatrix Nodes/IP locals) with identical lifetimes.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module cs_comp_acc
  use iso_c_binding
  use defutils
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(model_t), pointer, save :: acc_model => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(variable_t), pointer, save :: acc_velovar => null()
  integer, allocatable, save :: acc_indexes(:)
  ! LocalMatrix context (mirrors its SAVE Nodes and IP local).
  type(nodes_t), target, save :: acc_nodes_lm
  type(gaussintegrationpoints_t), save :: acc_ip_lm

contains

  ! -- entry context ---------------------------------------------------
  subroutine csc_begin(h_model, h_solver) bind(c, name="csc_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine csc_begin

  function csc_dim() result(d) bind(c, name="csc_dim")
    integer(c_int) :: d
    d = CoordinateSystemDimension()
  end function csc_dim

  function csc_max_elem_dofs() result(n) bind(c, name="csc_max_elem_dofs")
    integer(c_int) :: n
    type(mesh_t), pointer :: m
    m => GetMesh()
    n = m % maxelementdofs
  end function csc_max_elem_dofs

  subroutine csc_const_bulk(val, found) bind(c, name="csc_const_bulk")
    integer(c_int) :: val, found
    type(valuelist_t), pointer :: sp
    logical :: v, f
    sp => GetSolverParams()
    v = GetLogical(sp, 'Constant Bulk Matrix', f)
    val = merge(1, 0, v)
    found = merge(1, 0, f)
  end subroutine csc_const_bulk

  function csc_const_bulk_in_use(cb) result(r) bind(c, name="csc_const_bulk_in_use")
    integer(c_int), value :: cb
    integer(c_int) :: r
    logical :: in_use
    in_use = (cb /= 0) .and. associated(acc_solver % matrix % bulkvalues)
    r = merge(1, 0, in_use)
  end function csc_const_bulk_in_use

  subroutine csc_velovar() bind(c, name="csc_velovar")
    type(mesh_t), pointer :: m
    m => GetMesh()
    acc_velovar => VariableGet(m % variables, "VelocityTot")
  end subroutine csc_velovar

  subroutine csc_alloc_indexes(n) bind(c, name="csc_alloc_indexes")
    integer(c_int), value :: n
    if (allocated(acc_indexes)) deallocate(acc_indexes)
    allocate(acc_indexes(n))
  end subroutine csc_alloc_indexes

  subroutine csc_default_init(cb) bind(c, name="csc_default_init")
    integer(c_int), value :: cb
    logical :: use_cb
    use_cb = (cb /= 0)
    call DefaultInitialize(acc_solver, use_cb)
  end subroutine csc_default_init

  ! -- active-element driver -------------------------------------------
  function csc_nof_active() result(n) bind(c, name="csc_nof_active")
    integer(c_int) :: n
    n = GetNOFActive()
  end function csc_nof_active

  subroutine csc_begin_active(t, n, nb, nd) bind(c, name="csc_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n, nb, nd
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes()
    nb = GetElementNOFBDOFs()
    nd = GetElementDOFs(acc_indexes)
  end subroutine csc_begin_active

  subroutine csc_fill_cvelo(dim, nd, cvelo) bind(c, name="csc_fill_cvelo")
    integer(c_int), value :: dim, nd
    real(c_double) :: cvelo(dim, *)
    integer :: i
    ! Exact transcription of the previous-iterate gather.
    do i = 1, dim
      cvelo(i, 1:nd) = acc_velovar % values( &
          acc_velovar % dofs * (acc_velovar % perm(acc_indexes(1:nd)) - 1) + i)
    end do
  end subroutine csc_fill_cvelo

  function csc_viscosity() result(v) bind(c, name="csc_viscosity")
    real(c_double) :: v
    v = ListGetConstReal(GetMaterial(), 'Viscosity')
  end function csc_viscosity

  ! -- LocalMatrix context (IP loop owned by Rust) ----------------------
  subroutine csc_lm_nodes() bind(c, name="csc_lm_nodes")
    call GetElementNodes(acc_nodes_lm)
  end subroutine csc_lm_nodes

  function csc_lm_ip_begin() result(n) bind(c, name="csc_lm_ip_begin")
    integer(c_int) :: n
    acc_ip_lm = GaussPoints(acc_elem)
    n = acc_ip_lm % n
  end function csc_lm_ip_begin

  subroutine csc_lm_ip_point(t, u, v, w, s) bind(c, name="csc_lm_ip_point")
    integer(c_int), value :: t
    real(c_double) :: u, v, w, s
    u = acc_ip_lm % u(t)
    v = acc_ip_lm % v(t)
    w = acc_ip_lm % w(t)
    s = acc_ip_lm % s(t)
  end subroutine csc_lm_ip_point

  subroutine csc_lm_eleminfo(u, v, w, ntot, basis, dbasis, detj, stat) &
      bind(c, name="csc_lm_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: ntot
    real(c_double) :: basis(ntot), dbasis(ntot, 3), detj
    integer(c_int) :: stat
    logical :: ok
    ok = ElementInfo(acc_elem, acc_nodes_lm, u, v, w, detj, &
        basis(1:ntot), dbasis(1:ntot, 1:3))
    stat = merge(1, 0, ok)
  end subroutine csc_lm_eleminfo

  subroutine csc_lm_point(dim, n, nd, ntot, cb, visc, sdet, &
      cvelo, basis, dbasis, stiff, force) bind(c, name="csc_lm_point")
    integer(c_int), value :: dim, n, nd, ntot, cb
    real(c_double), value :: visc, sdet
    real(c_double) :: cvelo(dim, *), basis(ntot), dbasis(ntot, 3)
    real(c_double) :: stiff(ntot, ntot), force(ntot)
    ! Exact transcription of the per-IP body: divergence reduction, the
    ! elemental p-loop, and the bubble/boundary-DOF elimination. Same
    ! compiler, same evaluation order as the original internal procedure.
    real(c_double) :: divu
    integer :: i, p, q
    divu = 0.0d0
    do i = 1, dim
      divu = divu + sum(cvelo(i, 1:nd) * dbasis(1:nd, i))
    end do
    do p = 1, n
      if (cb == 0) then
        do q = 1, n
          stiff(p, q) = stiff(p, q) + sdet * basis(q) * basis(p)
        end do
      end if
      force(p) = force(p) + sdet * 2.0d0 * visc * divu * basis(p)
    end do
    do i = n + 1, ntot
      force(i) = 0.0d0
      stiff(i, :) = 0.0d0
      stiff(:, i) = 0.0d0
      stiff(i, i) = 1.0d0
    end do
  end subroutine csc_lm_point

  ! -- shared assembly tail ---------------------------------------------
  subroutine csc_update_eq(ld, g, f) bind(c, name="csc_update_eq")
    integer(c_int), value :: ld
    real(c_double) :: g(ld, ld), f(ld)
    call DefaultUpdateEquations(g(1:ld, 1:ld), f(1:ld))
  end subroutine csc_update_eq

  subroutine csc_update_force(ld, f) bind(c, name="csc_update_force")
    integer(c_int), value :: ld
    real(c_double) :: f(ld)
    call DefaultUpdateForce(f(1:ld))
  end subroutine csc_update_force

  subroutine csc_finish_bulk(cb_matrix, cb_in_use) bind(c, name="csc_finish_bulk")
    integer(c_int), value :: cb_matrix, cb_in_use
    logical :: bulk_update
    bulk_update = (cb_in_use == 0)
    if (cb_matrix /= 0) then
      call DefaultFinishBulkAssembly(BulkUpdate=bulk_update, RHSUpdate=.false.)
    else
      call DefaultFinishBulkAssembly()
    end if
  end subroutine csc_finish_bulk
  subroutine csc_finish() bind(c, name="csc_finish")
    call DefaultFinishAssembly()
  end subroutine csc_finish

  subroutine csc_dirichlet() bind(c, name="csc_dirichlet")
    call DefaultDirichletBCs()
  end subroutine csc_dirichlet

  function csc_solve() result(norm) bind(c, name="csc_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function csc_solve

  subroutine csc_fatal_alloc() bind(c, name="csc_fatal_alloc")
    call Fatal('CompressibilitySolver', 'Memory allocation error.')
  end subroutine csc_fatal_alloc

end module cs_comp_acc
