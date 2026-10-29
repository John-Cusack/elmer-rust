! T022-A FreeSurfaceReduced port: statement-level BIND(C) wrappers for
! the Elmer calls in fem/src/modules/FreeSurfaceReduced.F90, compiled
! against the built .mod files of the target configure (A2 shipped
! mechanism, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, formatted Message writes —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module mirrors the original's SAVE state (FreeSurfacePoints,
! LeftNeighbours, RightNeighbours, coords/coords2 pointer aliases,
! OrigCoords, UpwindPoints, NoFreeSurfaces, DrawDirection,
! SubroutineVisited, AllocationsDone) with identical lifetimes, and every
! formatted WRITE + Info/Warn call site is transcribed statement for
! statement so bit-identity holds by construction wherever the call
! sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module freesurfacereduced_acc
  use iso_c_binding
  use types
  use defutils
  use messages
  use meshbasics, only: FindNeighbourNodes
  use coordinatesystems
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(mesh_t), pointer, save :: acc_mesh => null()
  type(variable_t), pointer, save :: acc_var => null()
  type(variable_t), pointer, save :: acc_move => null()
  type(element_t), pointer, save :: acc_elem => null()
  ! SAVE mirrors of the original's permanent storage.
  integer, pointer, save :: acc_fsp(:) => null()
  integer, pointer, save :: acc_left(:) => null()
  integer, pointer, save :: acc_right(:) => null()
  real(c_double), pointer, save :: acc_orig(:) => null()
  real(c_double), pointer, save :: acc_coords(:) => null()
  real(c_double), pointer, save :: acc_coords2(:) => null()
  integer, pointer, save :: acc_upwind(:) => null()
  integer(c_int), save :: acc_nofree = 0
  integer(c_int), save :: acc_drawdir = 0
  integer(c_int), save :: acc_visited = 0
  integer(c_int), save :: acc_alloc = 0

contains

  ! -- entry state ----------------------------------------------------
  subroutine fsa_begin(h_model, h_solver) bind(c, name="fsa_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    acc_mesh => acc_solver % mesh
    acc_move => acc_solver % variable
  end subroutine fsa_begin

  subroutine fsa_find_flow() bind(c, name="fsa_find_flow")
    acc_var => VariableGet(acc_mesh % variables, 'Flow Solution', .TRUE.)
  end subroutine fsa_find_flow

  subroutine fsa_info_header() bind(c, name="fsa_info_header")
    call Info('FreeSurfaceReduced', '-------------------------------------', level=4)
    call Info('FreeSurfaceReduced', '2D Free Surface Solver:  ', level=4)
    call Info('FreeSurfaceReduced', '-------------------------------------', level=4)
  end subroutine fsa_info_header

  ! -- model scalars --------------------------------------------------
  function fsa_max_nodes() result(n) bind(c, name="fsa_max_nodes")
    integer(c_int) :: n
    n = acc_model % maxelementnodes
  end function fsa_max_nodes

  function fsa_num_nodes() result(n) bind(c, name="fsa_num_nodes")
    integer(c_int) :: n
    n = acc_model % numberofnodes
  end function fsa_num_nodes

  function fsa_num_bcs() result(n) bind(c, name="fsa_num_bcs")
    integer(c_int) :: n
    n = acc_model % numberofbcs
  end function fsa_num_bcs

  function fsa_alloc_done() result(r) bind(c, name="fsa_alloc_done")
    integer(c_int) :: r
    r = acc_alloc
  end function fsa_alloc_done

  subroutine fsa_alloc_first(n) bind(c, name="fsa_alloc_first")
    integer(c_int), value :: n
    integer :: istat
    allocate(acc_fsp(n), acc_left(n), acc_right(n), acc_orig(n), stat=istat)
    if (istat /= 0) call Warn('FreeSurfaceReduced', 'Memory allocation error')
    acc_alloc = 1
  end subroutine fsa_alloc_first

  subroutine fsa_zero_marks() bind(c, name="fsa_zero_marks")
    acc_fsp = 0
    acc_nofree = 0
  end subroutine fsa_zero_marks

  function fsa_visited() result(n) bind(c, name="fsa_visited")
    integer(c_int) :: n
    n = acc_visited
  end function fsa_visited

  subroutine fsa_inc_visited() bind(c, name="fsa_inc_visited")
    acc_visited = acc_visited + 1
  end subroutine fsa_inc_visited

  ! -- boundary scan --------------------------------------------------
  function fsa_num_bulk() result(n) bind(c, name="fsa_num_bulk")
    integer(c_int) :: n
    n = acc_mesh % numberofbulkelements
  end function fsa_num_bulk

  function fsa_num_bnd() result(n) bind(c, name="fsa_num_bnd")
    integer(c_int) :: n
    n = acc_mesh % numberofboundaryelements
  end function fsa_num_bnd

  subroutine fsa_begin_elem(t, elemcode, n, constraint) bind(c, name="fsa_begin_elem")
    integer(c_int), value :: t
    integer(c_int) :: elemcode, n, constraint
    acc_elem => acc_mesh % elements(t)
    elemcode = acc_elem % type % elementcode
    if (elemcode /= 101) acc_model % currentelement => acc_elem
    n = acc_elem % type % numberofnodes
    constraint = acc_elem % boundaryinfo % constraint
  end subroutine fsa_begin_elem

  function fsa_bc_tag(k) result(n) bind(c, name="fsa_bc_tag")
    integer(c_int), value :: k
    integer(c_int) :: n
    n = acc_model % bcs(k) % tag
  end function fsa_bc_tag

  function fsa_bc_reduced(k, found) result(r) bind(c, name="fsa_bc_reduced")
    integer(c_int), value :: k
    integer(c_int) :: found, r
    logical :: v, f
    v = ListGetLogical(acc_model % bcs(k) % values, 'Free Surface Reduced', f)
    found = merge(1, 0, f)
    r = merge(1, 0, v)
  end function fsa_bc_reduced

  subroutine fsa_bc_surface(k, surf, found) bind(c, name="fsa_bc_surface")
    integer(c_int), value :: k
    integer(c_int) :: surf, found
    logical :: f
    surf = ListGetInteger(acc_model % bcs(k) % values, 'Free Surface Number', f)
    found = merge(1, 0, f)
  end subroutine fsa_bc_surface

  function fsa_bc_bottom(k, found) result(r) bind(c, name="fsa_bc_bottom")
    integer(c_int), value :: k
    integer(c_int) :: found, r
    logical :: v, f
    v = ListGetLogical(acc_model % bcs(k) % values, 'Free Surface Bottom', f)
    found = merge(1, 0, f)
    r = merge(1, 0, v)
  end function fsa_bc_bottom

  subroutine fsa_apply_surface(n, surf) bind(c, name="fsa_apply_surface")
    integer(c_int), value :: n, surf
    acc_fsp(acc_elem % nodeindexes(1:n)) = surf
  end subroutine fsa_apply_surface

  subroutine fsa_apply_bottom(n) bind(c, name="fsa_apply_bottom")
    integer(c_int), value :: n
    acc_fsp(acc_elem % nodeindexes(1:n)) = -1
  end subroutine fsa_apply_bottom

  subroutine fsa_note_surface(surf) bind(c, name="fsa_note_surface")
    integer(c_int), value :: surf
    acc_nofree = max(acc_nofree, surf)
  end subroutine fsa_note_surface

  function fsa_nofree() result(n) bind(c, name="fsa_nofree")
    integer(c_int) :: n
    n = acc_nofree
  end function fsa_nofree

  subroutine fsa_info_nsurf() bind(c, name="fsa_info_nsurf")
    write(message, '(A,I1)') 'Separate Free surfaces found ', acc_nofree
    call Info('FreeSurfaceReduced', message, level=5)
  end subroutine fsa_info_nsurf

  subroutine fsa_alloc_upwind() bind(c, name="fsa_alloc_upwind")
    allocate(acc_upwind(acc_nofree))
    acc_upwind = 0
  end subroutine fsa_alloc_upwind

  ! -- mesh / flow field reads ----------------------------------------
  function fsa_fsp(t) result(n) bind(c, name="fsa_fsp")
    integer(c_int), value :: t
    integer(c_int) :: n
    n = acc_fsp(t)
  end function fsa_fsp

  function fsa_node_x(t) result(x) bind(c, name="fsa_node_x")
    integer(c_int), value :: t
    real(c_double) :: x
    x = acc_mesh % nodes % x(t)
  end function fsa_node_x

  function fsa_node_y(t) result(y) bind(c, name="fsa_node_y")
    integer(c_int), value :: t
    real(c_double) :: y
    y = acc_mesh % nodes % y(t)
  end function fsa_node_y

  subroutine fsa_set_drawdir(d) bind(c, name="fsa_set_drawdir")
    integer(c_int), value :: d
    acc_drawdir = d
  end subroutine fsa_set_drawdir

  function fsa_drawdir() result(d) bind(c, name="fsa_drawdir")
    integer(c_int) :: d
    d = acc_drawdir
  end function fsa_drawdir

  function fsa_flowt(t, comp) result(v) bind(c, name="fsa_flowt")
    integer(c_int), value :: t, comp
    real(c_double) :: v
    v = acc_var % values(acc_var % dofs * (acc_var % perm(t) - 1) + comp)
  end function fsa_flowt

  function fsa_flowk(k, comp) result(v) bind(c, name="fsa_flowk")
    integer(c_int), value :: k, comp
    real(c_double) :: v
    v = acc_var % values(acc_var % dofs * (k - 1) + comp)
  end function fsa_flowk

  function fsa_var_perm(t) result(k) bind(c, name="fsa_var_perm")
    integer(c_int), value :: t
    integer(c_int) :: k
    k = acc_var % perm(t)
  end function fsa_var_perm

  function fsa_move_perm(t) result(k) bind(c, name="fsa_move_perm")
    integer(c_int), value :: t
    integer(c_int) :: k
    k = acc_move % perm(t)
  end function fsa_move_perm

  subroutine fsa_info_drawdir() bind(c, name="fsa_info_drawdir")
    write(message, '(A,I3)') 'DrawDirection was found to be', acc_drawdir
    call Info('FreeSurfaceReduced', message, level=5)
  end subroutine fsa_info_drawdir

  subroutine fsa_neighbours() bind(c, name="fsa_neighbours")
    real(c_double) :: rnorm(3)
    rnorm = 0.0_c_double
    if (abs(acc_drawdir) == 2) then
      rnorm(1) = -1.0_c_double
    else
      rnorm(2) = -1.0_c_double
    end if
    call FindNeighbourNodes(acc_mesh, rnorm, acc_left)
    call FindNeighbourNodes(acc_mesh, -1.0_c_double * rnorm, acc_right)
  end subroutine fsa_neighbours

  subroutine fsa_bind_coords() bind(c, name="fsa_bind_coords")
    if (abs(acc_drawdir) == 2) then
      acc_coords => acc_mesh % nodes % x
      acc_coords2 => acc_mesh % nodes % y
      acc_orig = acc_mesh % nodes % x
    else
      acc_coords => acc_mesh % nodes % y
      acc_coords2 => acc_mesh % nodes % x
      acc_orig = acc_mesh % nodes % y
    end if
  end subroutine fsa_bind_coords

  function fsa_get_upwind(s) result(t) bind(c, name="fsa_get_upwind")
    integer(c_int), value :: s
    integer(c_int) :: t
    t = acc_upwind(s)
  end function fsa_get_upwind

  subroutine fsa_set_upwind(s, t) bind(c, name="fsa_set_upwind")
    integer(c_int), value :: s, t
    acc_upwind(s) = t
  end subroutine fsa_set_upwind

  subroutine fsa_info_upwind(s) bind(c, name="fsa_info_upwind")
    integer(c_int), value :: s
    write(message, '(A,I6)') 'First upwinds node for surface is', acc_upwind(s)
    call Info('FreeSurfaceReduced', message, level=5)
  end subroutine fsa_info_upwind

  subroutine fsa_zero_move() bind(c, name="fsa_zero_move")
    acc_move % values = 0.0_c_double
  end subroutine fsa_zero_move

  function fsa_upwind_coord(t) result(x) bind(c, name="fsa_upwind_coord")
    integer(c_int), value :: t
    real(c_double) :: x
    if (abs(acc_drawdir) == 2) then
      x = acc_mesh % nodes % y(t)
    else
      x = acc_mesh % nodes % x(t)
    end if
  end function fsa_upwind_coord

  function fsa_coord(t) result(x) bind(c, name="fsa_coord")
    integer(c_int), value :: t
    real(c_double) :: x
    x = acc_coords(t)
  end function fsa_coord

  function fsa_coord2(t) result(x) bind(c, name="fsa_coord2")
    integer(c_int), value :: t
    real(c_double) :: x
    x = acc_coords2(t)
  end function fsa_coord2

  function fsa_left(t) result(j) bind(c, name="fsa_left")
    integer(c_int), value :: t
    integer(c_int) :: j
    j = acc_left(t)
  end function fsa_left

  function fsa_right(t) result(j) bind(c, name="fsa_right")
    integer(c_int), value :: t
    integer(c_int) :: j
    j = acc_right(t)
  end function fsa_right

  ! -- warnings (exact call sites) --------------------------------------
  subroutine fsa_warn_perm() bind(c, name="fsa_warn_perm")
    call Warn('FreeSurfaceReduced', 'Problems as Perm(i) < 0')
  end subroutine fsa_warn_perm

  subroutine fsa_warn_limits() bind(c, name="fsa_warn_limits")
    call Warn('FreeSurfaceReduced', 'Invalid integration limits')
  end subroutine fsa_warn_limits

  subroutine fsa_warn_noflow() bind(c, name="fsa_warn_noflow")
    call Warn('FreeSurfaceReduced', 'no flow solution for given point')
  end subroutine fsa_warn_noflow

  subroutine fsa_warn_dx0() bind(c, name="fsa_warn_dx0")
    call Warn('FreeSurfaceReduced', 'Integration path has dx=0')
  end subroutine fsa_warn_dx0

  ! -- solver parameters ------------------------------------------------
  function fsa_axisym() result(r) bind(c, name="fsa_axisym")
    integer(c_int) :: r
    integer :: cs
    cs = CurrentCoordinateSystem()
    if (cs == CylindricSymmetric .or. cs == AxisSymmetric) then
      r = 1
    else
      r = 0
    end if
  end function fsa_axisym

  subroutine fsa_relax(r, found) bind(c, name="fsa_relax")
    real(c_double) :: r
    integer(c_int) :: found
    logical :: f
    r = GetCReal(acc_solver % values, 'Nonlinear System Relaxation Factor', f)
    found = merge(1, 0, f)
  end subroutine fsa_relax

  subroutine fsa_mapping(m, found) bind(c, name="fsa_mapping")
    integer(c_int) :: m, found
    logical :: v, f
    v = ListGetLogical(acc_solver % values, 'Perform Mapping', f)
    found = merge(1, 0, f)
    m = merge(1, 0, v)
  end subroutine fsa_mapping

  ! -- mapping chain (exact transcription of the PerformMapping loop) ---
  subroutine fsa_map_chain(istart, iend, x0, xnew, xold, relax) &
      bind(c, name="fsa_map_chain")
    integer(c_int), value :: istart, iend
    real(c_double), value :: x0, xnew, xold, relax
    integer :: i, j
    real(c_double) :: x1
    i = istart
    do
      i = acc_right(i)
      x1 = x0 + (xnew - x0) * (acc_coords(i) - x0) / (xold - x0)
      j = acc_move % perm(i)
      acc_move % values(j) = (1.0_c_double - relax) * acc_move % values(j) + &
          relax * (x1 - acc_orig(i))
      acc_coords(i) = acc_orig(i) + acc_move % values(j)
      if (i == iend) exit
    end do
  end subroutine fsa_map_chain

  subroutine fsa_move_point(iend, xnew, relax) bind(c, name="fsa_move_point")
    integer(c_int), value :: iend
    real(c_double), value :: xnew, relax
    integer :: j
    j = acc_move % perm(iend)
    acc_move % values(j) = (1.0_c_double - relax) * acc_move % values(j) + &
        relax * (xnew - acc_orig(iend))
  end subroutine fsa_move_point

  ! -- per-surface reports (exact transcription) --------------------------
  subroutine fsa_info_surf(s, multi, fsum0, maxerror, thickmin, thickmax) &
      bind(c, name="fsa_info_surf")
    integer(c_int), value :: s, multi
    real(c_double), value :: fsum0, maxerror, thickmin, thickmax
    if (multi /= 0) then
      write(message, '(A,I1)') 'Free surface number ', s
      call Info('FreeSurfaceReduced', message, level=4)
    end if
    write(message, '(A,T30,ES17.6E2)') 'Free surface flux:', fsum0
    call Info('FreeSurfaceReduced', message, level=4)
    write(message, '(A,T30,ES17.6E2)') 'Free surface error:', maxerror
    call Info('FreeSurfaceReduced', message, level=4)
    write(message, '(A,T30,ES17.6E2)') 'Minimum thickness:', thickmin
    call Info('FreeSurfaceReduced', message, level=4)
    write(message, '(A,T30,ES17.6E2)') 'Maximum thickness:', thickmax
    call Info('FreeSurfaceReduced', message, level=4)
  end subroutine fsa_info_surf

  subroutine fsa_set_norm(n) bind(c, name="fsa_set_norm")
    real(c_double), value :: n
    acc_solver % variable % norm = n
  end subroutine fsa_set_norm

  subroutine fsa_add_results(fsum0, maxerror, thickmin, thickmax) &
      bind(c, name="fsa_add_results")
    real(c_double), value :: fsum0, maxerror, thickmin, thickmax
    call ListAddConstReal(acc_model % simulation, 'res: Free surface flux', fsum0)
    call ListAddConstReal(acc_model % simulation, 'res: Free surface max. error', maxerror)
    call ListAddConstReal(acc_model % simulation, 'res: Free thickness min', thickmin)
    call ListAddConstReal(acc_model % simulation, 'res: Free thickness max', thickmax)
  end subroutine fsa_add_results

end module freesurfacereduced_acc
