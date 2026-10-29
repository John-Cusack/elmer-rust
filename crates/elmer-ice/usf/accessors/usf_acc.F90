! ElmerIceUSF whole-library port: statement-level BIND(C) wrappers for the
! Elmer calls in the elmerice UserFunctions sources, compiled against the
! .mod files of the target ice-enabled configure.
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays, C
! strings as (ptr,len)). All Fortran-object traffic — TYPE(...) handles,
! OPTIONAL/POINTER dummies, assumed-shape descriptors, hidden lengths —
! stays inside Fortran, so the Rust side is pure C calls.
!
! The module mirrors the entry-called context in SAVE state (current
! model handle set per entry call), exactly like the cs_*_acc modules of
! the ConsistentSplitting port.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module usf_acc
  use iso_c_binding
  use defutils
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  ! SAVE'd mask perms for the GlacierMeshMetric entry (mirrors the
  ! source's SAVE'd FrontPerm/BottomPerm/SurfacePerm).
  integer, pointer, save :: acc_front_perm(:) => null()
  integer, pointer, save :: acc_bottom_perm(:) => null()
  integer, pointer, save :: acc_surface_perm(:) => null()

contains

  ! -- entry context ---------------------------------------------------
  subroutine usf_begin(h_model) bind(c, name="usf_begin")
    type(c_ptr), value :: h_model
    call c_f_pointer(h_model, acc_model)
  end subroutine usf_begin

  ! -- pure-float helpers (routed here so libm-sensitive ops evaluate in
  ! -- gfortran, bit-identical to the reference by construction) ---------
  function usf_pow_dp(base, expo) result(r) bind(c, name="usf_pow_dp")
    real(c_double), value :: base, expo
    real(c_double) :: r
    r = base ** expo
  end function usf_pow_dp

  ! -- ProjUtils projection math (ElmerIceUtils.so; identical routines) --
  subroutine usf_xy2lonlat(x, y, lon, lat) bind(c, name="usf_xy2lonlat")
    use projutils
    real(c_double), value :: x, y
    real(c_double) :: lon, lat
    call xy2LonLat(x, y, lon, lat)
  end subroutine usf_xy2lonlat

  subroutine usf_lonlat2xy(lon, lat, x, y) bind(c, name="usf_lonlat2xy")
    use projutils
    real(c_double), value :: lon, lat
    real(c_double) :: x, y
    call LonLat2xy(lon, lat, x, y)
  end subroutine usf_lonlat2xy

  ! -- model / element / material queries --------------------------------
  subroutine usf_current_element(h_elem) bind(c, name="usf_current_element")
    type(c_ptr) :: h_elem
    type(element_t), pointer :: e
    e => acc_model % currentelement
    if (associated(e)) then
      h_elem = c_loc(e)
    else
      h_elem = c_null_ptr
    end if
  end subroutine usf_current_element

  subroutine usf_get_material(h_elem, h_list) bind(c, name="usf_get_material")
    type(c_ptr), value :: h_elem
    type(c_ptr) :: h_list
    type(element_t), pointer :: e
    type(valuelist_t), pointer :: m
    call c_f_pointer(h_elem, e)
    m => GetMaterial(e)
    if (associated(m)) then
      h_list = c_loc(m)
    else
      h_list = c_null_ptr
    end if
  end subroutine usf_get_material

  function usf_get_const_real(h_list, name, namelen, found) result(v) &
      bind(c, name="usf_get_const_real")
    type(c_ptr), value :: h_list
    type(c_ptr), value :: name
    integer(c_int), value :: namelen
    integer(c_int) :: found
    real(c_double) :: v
    type(valuelist_t), pointer :: l
    character(len=:), allocatable :: key
    character(kind=c_char), pointer :: cbuf(:)
    logical :: f
    integer :: i
    call c_f_pointer(h_list, l)
    call c_f_pointer(name, cbuf, [namelen])
    allocate(character(len=namelen) :: key)
    do i = 1, namelen
      key(i:i) = cbuf(i)
    end do
    v = GetConstReal(l, key, f)
    found = merge(1, 0, f)
  end function usf_get_const_real

  subroutine usf_get_bc(h_elem, h_list) bind(c, name="usf_get_bc")
    type(c_ptr), value :: h_elem
    type(c_ptr) :: h_list
    type(element_t), pointer :: e
    type(valuelist_t), pointer :: b
    call c_f_pointer(h_elem, e)
    b => GetBC(e)
    if (associated(b)) then
      h_list = c_loc(b)
    else
      h_list = c_null_ptr
    end if
  end subroutine usf_get_bc

  subroutine usf_get_constants(h_list) bind(c, name="usf_get_constants")
    type(c_ptr) :: h_list
    type(valuelist_t), pointer :: c
    c => GetConstants()
    if (associated(c)) then
      h_list = c_loc(c)
    else
      h_list = c_null_ptr
    end if
  end subroutine usf_get_constants

  subroutine usf_get_logical(h_list, name, namelen, val, found) &
      bind(c, name="usf_get_logical")
    type(c_ptr), value :: h_list
    type(c_ptr), value :: name
    integer(c_int), value :: namelen
    integer(c_int) :: val, found
    type(valuelist_t), pointer :: l
    character(len=:), allocatable :: key
    character(kind=c_char), pointer :: cbuf(:)
    logical :: v, f
    integer :: i
    call c_f_pointer(h_list, l)
    call c_f_pointer(name, cbuf, [namelen])
    allocate(character(len=namelen) :: key)
    do i = 1, namelen
      key(i:i) = cbuf(i)
    end do
    v = GetLogical(l, key, f)
    val = merge(1, 0, v)
    found = merge(1, 0, f)
  end subroutine usf_get_logical

  subroutine usf_list_get_logical(h_list, name, namelen, val, found) &
      bind(c, name="usf_list_get_logical")
    type(c_ptr), value :: h_list
    type(c_ptr), value :: name
    integer(c_int), value :: namelen
    integer(c_int) :: val, found
    type(valuelist_t), pointer :: l
    character(len=:), allocatable :: key
    character(kind=c_char), pointer :: cbuf(:)
    logical :: v, f
    integer :: i
    call c_f_pointer(h_list, l)
    call c_f_pointer(name, cbuf, [namelen])
    allocate(character(len=namelen) :: key)
    do i = 1, namelen
      key(i:i) = cbuf(i)
    end do
    v = ListGetLogical(l, key, f)
    val = merge(1, 0, v)
    found = merge(1, 0, f)
  end subroutine usf_list_get_logical

  subroutine usf_info(caller, ncaller, msg, nmsg, level) bind(c, name="usf_info")
    type(c_ptr), value :: caller, msg
    integer(c_int), value :: ncaller, nmsg, level
    character(kind=c_char), pointer :: cbuf(:)
    character(len=:), allocatable :: fcaller, fmsg
    integer :: i
    call c_f_pointer(caller, cbuf, [ncaller])
    allocate(character(len=ncaller) :: fcaller)
    do i = 1, ncaller
      fcaller(i:i) = cbuf(i)
    end do
    call c_f_pointer(msg, cbuf, [nmsg])
    allocate(character(len=nmsg) :: fmsg)
    do i = 1, nmsg
      fmsg(i:i) = cbuf(i)
    end do
    call Info(fcaller, fmsg, level=level)
  end subroutine usf_info

  subroutine usf_fatal(caller, ncaller, msg, nmsg) bind(c, name="usf_fatal")
    type(c_ptr), value :: caller, msg
    integer(c_int), value :: ncaller, nmsg
    character(kind=c_char), pointer :: cbuf(:)
    character(len=:), allocatable :: fcaller, fmsg
    integer :: i
    call c_f_pointer(caller, cbuf, [ncaller])
    allocate(character(len=ncaller) :: fcaller)
    do i = 1, ncaller
      fcaller(i:i) = cbuf(i)
    end do
    call c_f_pointer(msg, cbuf, [nmsg])
    allocate(character(len=nmsg) :: fmsg)
    do i = 1, nmsg
      fmsg(i:i) = cbuf(i)
    end do
    call Fatal(fcaller, fmsg)
  end subroutine usf_fatal

  subroutine usf_model_variables(h_vars) bind(c, name="usf_model_variables")
    type(c_ptr) :: h_vars
    type(variable_t), pointer :: v
    v => acc_model % variables
    if (associated(v)) then
      h_vars = c_loc(v)
    else
      h_vars = c_null_ptr
    end if
  end subroutine usf_model_variables

  subroutine usf_variable_get(h_vars, name, namelen, unfoundfatal, thisonly, &
      h_var) bind(c, name="usf_variable_get")
    type(c_ptr), value :: h_vars
    type(c_ptr), value :: name
    integer(c_int), value :: namelen, unfoundfatal, thisonly
    type(c_ptr) :: h_var
    type(variable_t), pointer :: vars, v
    character(len=:), allocatable :: key
    character(kind=c_char), pointer :: cbuf(:)
    logical :: uff, to
    integer :: i
    call c_f_pointer(h_vars, vars)
    call c_f_pointer(name, cbuf, [namelen])
    allocate(character(len=namelen) :: key)
    do i = 1, namelen
      key(i:i) = cbuf(i)
    end do
    uff = (unfoundfatal /= 0)
    to = (thisonly /= 0)
    v => VariableGet(vars, key, ThisOnly=to, UnFoundFatal=uff)
    if (associated(v)) then
      h_var = c_loc(v)
    else
      h_var = c_null_ptr
    end if
  end subroutine usf_variable_get

  function usf_var_value_at(h_var, node) result(v) &
      bind(c, name="usf_var_value_at")
    type(c_ptr), value :: h_var
    integer(c_int), value :: node
    real(c_double) :: v
    type(variable_t), pointer :: var
    call c_f_pointer(h_var, var)
    v = var % values(var % perm(node))
  end function usf_var_value_at

  subroutine usf_get_string(h_list, key, nkey, buf, nbuf, outlen, found) &
      bind(c, name="usf_get_string")
    type(c_ptr), value :: h_list, key
    integer(c_int), value :: nkey, nbuf
    type(c_ptr), value :: buf
    integer(c_int) :: outlen, found
    type(valuelist_t), pointer :: l
    character(len=:), allocatable :: fkey, val
    character(kind=c_char), pointer :: cbuf(:), cout(:)
    logical :: f
    integer :: i, m
    call c_f_pointer(h_list, l)
    call c_f_pointer(key, cbuf, [nkey])
    allocate(character(len=nkey) :: fkey)
    do i = 1, nkey
      fkey(i:i) = cbuf(i)
    end do
    val = GetString(l, fkey, f)
    found = merge(1, 0, f)
    m = min(len(val), nbuf)
    if (m > 0) then
      call c_f_pointer(buf, cout, [m])
      do i = 1, m
        cout(i) = val(i:i)
      end do
    end if
    outlen = len_trim(val)
  end subroutine usf_get_string

  function usf_get_element_nof_nodes(h_elem) result(n) &
      bind(c, name="usf_get_element_nof_nodes")
    type(c_ptr), value :: h_elem
    integer(c_int) :: n
    type(element_t), pointer :: e
    call c_f_pointer(h_elem, e)
    n = GetElementNOFNodes(e)
  end function usf_get_element_nof_nodes

  subroutine usf_element_nodeindexes(h_elem, buf, n) &
      bind(c, name="usf_element_nodeindexes")
    type(c_ptr), value :: h_elem
    type(c_ptr), value :: buf
    integer(c_int), value :: n
    type(element_t), pointer :: e
    integer(c_int), pointer :: out(:)
    call c_f_pointer(h_elem, e)
    call c_f_pointer(buf, out, [n])
    out(1:n) = e % nodeindexes(1:n)
  end subroutine usf_element_nodeindexes

  function usf_max_element_nodes() result(n) &
      bind(c, name="usf_max_element_nodes")
    integer(c_int) :: n
    n = acc_model % maxelementnodes
  end function usf_max_element_nodes

  subroutine usf_get_bodyforce(h_bf) bind(c, name="usf_get_bodyforce")
    type(c_ptr) :: h_bf
    type(valuelist_t), pointer :: b
    b => GetBodyForce()
    if (associated(b)) then
      h_bf = c_loc(b)
    else
      h_bf = c_null_ptr
    end if
  end subroutine usf_get_bodyforce

  function usf_get_bodyforce_id() result(id) &
      bind(c, name="usf_get_bodyforce_id")
    integer(c_int) :: id
    id = GetBodyForceId()
  end function usf_get_bodyforce_id

  subroutine usf_bodyforce_values(bf_id, h_list) &
      bind(c, name="usf_bodyforce_values")
    integer(c_int), value :: bf_id
    type(c_ptr) :: h_list
    type(valuelist_t), pointer :: l
    l => acc_model % bodyforces(bf_id) % values
    if (associated(l)) then
      h_list = c_loc(l)
    else
      h_list = c_null_ptr
    end if
  end subroutine usf_bodyforce_values

  function usf_get_material_id(h_elem, found) result(id) &
      bind(c, name="usf_get_material_id")
    type(c_ptr), value :: h_elem
    integer(c_int) :: found
    integer(c_int) :: id
    type(element_t), pointer :: e
    logical :: f
    call c_f_pointer(h_elem, e)
    id = GetMaterialId(e, f)
    found = merge(1, 0, f)
  end function usf_get_material_id

  subroutine usf_list_get_real(h_list, key, nkey, n, idx, unfoundfatal, &
      found, out) bind(c, name="usf_list_get_real")
    type(c_ptr), value :: h_list, key, idx, out
    integer(c_int), value :: nkey, n, unfoundfatal
    integer(c_int) :: found
    type(valuelist_t), pointer :: l
    character(len=:), allocatable :: fkey
    character(kind=c_char), pointer :: cbuf(:)
    integer(c_int), pointer :: fidx(:)
    real(c_double), pointer :: fout(:)
    logical :: uff, f
    integer :: i
    call c_f_pointer(h_list, l)
    call c_f_pointer(key, cbuf, [nkey])
    allocate(character(len=nkey) :: fkey)
    do i = 1, nkey
      fkey(i:i) = cbuf(i)
    end do
    call c_f_pointer(idx, fidx, [n])
    call c_f_pointer(out, fout, [n])
    uff = (unfoundfatal /= 0)
    fout(1:n) = ListGetReal(l, fkey, n, fidx, f, UnfoundFatal=uff)
    found = merge(1, 0, f)
  end subroutine usf_list_get_real

  function usf_sum_g2_over_n(g, m, n) result(s) &
      bind(c, name="usf_sum_g2_over_n")
    type(c_ptr), value :: g
    integer(c_int), value :: m, n
    real(c_double) :: s
    real(c_double), pointer :: fg(:)
    call c_f_pointer(g, fg, [m])
    s = sum(fg(1:m) ** 2.0_dp / real(n, kind=c_double))
  end function usf_sum_g2_over_n

  subroutine usf_simulation_values(h_list) bind(c, name="usf_simulation_values")
    type(c_ptr) :: h_list
    type(valuelist_t), pointer :: s
    s => acc_model % simulation
    if (associated(s)) then
      h_list = c_loc(s)
    else
      h_list = c_null_ptr
    end if
  end subroutine usf_simulation_values

  function usf_number_of_solvers() result(n) bind(c, name="usf_number_of_solvers")
    integer(c_int) :: n
    n = acc_model % numberofsolvers
  end function usf_number_of_solvers

  subroutine usf_solver_variable_name(solver_i, buf, nbuf, outlen) &
      bind(c, name="usf_solver_variable_name")
    integer(c_int), value :: solver_i
    type(c_ptr), value :: buf
    integer(c_int), value :: nbuf
    integer(c_int) :: outlen
    character(kind=c_char), pointer :: cout(:)
    character(len=:), allocatable :: nm
    integer :: m, i
    nm = acc_model % solvers(solver_i) % variable % name
    m = min(len(nm), nbuf)
    if (m > 0) then
      call c_f_pointer(buf, cout, [m])
      do i = 1, m
        cout(i) = nm(i:i)
      end do
    end if
    outlen = len(nm)
  end subroutine usf_solver_variable_name

  subroutine usf_solver_mesh_variables(solver_i, h_vars) &
      bind(c, name="usf_solver_mesh_variables")
    integer(c_int), value :: solver_i
    type(c_ptr) :: h_vars
    type(variable_t), pointer :: v
    v => acc_model % solvers(solver_i) % mesh % variables
    if (associated(v)) then
      h_vars = c_loc(v)
    else
      h_vars = c_null_ptr
    end if
  end subroutine usf_solver_mesh_variables

  subroutine usf_calculate_nodal_weights(solver_i, varname, nvarname) &
      bind(c, name="usf_calculate_nodal_weights")
    integer(c_int), value :: solver_i, nvarname
    type(c_ptr), value :: varname
    character(kind=c_char), pointer :: cbuf(:)
    character(len=:), allocatable :: key
    integer :: i
    call c_f_pointer(varname, cbuf, [nvarname])
    allocate(character(len=nvarname) :: key)
    do i = 1, nvarname
      key(i:i) = cbuf(i)
    end do
    call CalculateNodalWeights(acc_model % solvers(solver_i), .FALSE., &
        VarName=key)
  end subroutine usf_calculate_nodal_weights

  subroutine usf_get_bodyforce_of(h_elem, h_bf) &
      bind(c, name="usf_get_bodyforce_of")
    type(c_ptr), value :: h_elem
    type(c_ptr) :: h_bf
    type(element_t), pointer :: e
    type(valuelist_t), pointer :: b
    call c_f_pointer(h_elem, e)
    b => GetBodyForce(e)
    if (associated(b)) then
      h_bf = c_loc(b)
    else
      h_bf = c_null_ptr
    end if
  end subroutine usf_get_bodyforce_of

  function usf_var_perm_at(h_var, node) result(p) &
      bind(c, name="usf_var_perm_at")
    type(c_ptr), value :: h_var
    integer(c_int), value :: node
    integer(c_int) :: p
    type(variable_t), pointer :: var
    call c_f_pointer(h_var, var)
    p = var % perm(node)
  end function usf_var_perm_at

  function usf_var_value_at_idx(h_var, idx) result(v) &
      bind(c, name="usf_var_value_at_idx")
    type(c_ptr), value :: h_var
    integer(c_int), value :: idx
    real(c_double) :: v
    type(variable_t), pointer :: var
    call c_f_pointer(h_var, var)
    v = var % values(idx)
  end function usf_var_value_at_idx

  subroutine usf_get_material_current(h_list) &
      bind(c, name="usf_get_material_current")
    type(c_ptr) :: h_list
    type(valuelist_t), pointer :: m
    m => GetMaterial()
    if (associated(m)) then
      h_list = c_loc(m)
    else
      h_list = c_null_ptr
    end if
  end subroutine usf_get_material_current

  function usf_coord_dim() result(d) bind(c, name="usf_coord_dim")
    integer(c_int) :: d
    d = CoordinateSystemDimension()
  end function usf_coord_dim

  function usf_even_random() result(r) bind(c, name="usf_even_random")
    use generalutils
    real(c_double) :: r
    r = EvenRandom()
  end function usf_even_random

  subroutine usf_write_msg_af104(buf, nbuf, prefix, nprefix, value, outlen) &
      bind(c, name="usf_write_msg_af104")
    type(c_ptr), value :: buf, prefix
    integer(c_int), value :: nbuf, nprefix
    real(c_double), value :: value
    integer(c_int) :: outlen
    character(kind=c_char), pointer :: cout(:), cpre(:)
    character(len=:), allocatable :: fpre, msg
    integer :: m, i
    call c_f_pointer(prefix, cpre, [nprefix])
    allocate(character(len=nprefix) :: fpre)
    do i = 1, nprefix
      fpre(i:i) = cpre(i)
    end do
    allocate(character(len=nbuf) :: msg)
    write(msg, '(A,F10.4)') fpre, value
    m = min(len_trim(msg), nbuf)
    if (m > 0) then
      call c_f_pointer(buf, cout, [m])
      do i = 1, m
        cout(i) = msg(i:i)
      end do
    end if
    outlen = len_trim(msg)
  end subroutine usf_write_msg_af104

  subroutine usf_write_msg_al1(buf, nbuf, prefix, nprefix, value, outlen) &
      bind(c, name="usf_write_msg_al1")
    type(c_ptr), value :: buf, prefix
    integer(c_int), value :: nbuf, nprefix
    integer(c_int), value :: value
    integer(c_int) :: outlen
    character(kind=c_char), pointer :: cout(:), cpre(:)
    character(len=:), allocatable :: fpre, msg
    logical :: v
    integer :: m, i
    call c_f_pointer(prefix, cpre, [nprefix])
    allocate(character(len=nprefix) :: fpre)
    do i = 1, nprefix
      fpre(i:i) = cpre(i)
    end do
    v = (value /= 0)
    allocate(character(len=nbuf) :: msg)
    write(msg, '(A,L1)') fpre, v
    m = min(len_trim(msg), nbuf)
    if (m > 0) then
      call c_f_pointer(buf, cout, [m])
      do i = 1, m
        cout(i) = msg(i:i)
      end do
    end if
    outlen = len_trim(msg)
  end subroutine usf_write_msg_al1

  subroutine usf_sortd(n, a, b) bind(c, name="usf_sortd")
    use generalutils
    integer(c_int), value :: n
    type(c_ptr), value :: a, b
    real(c_double), pointer :: fa(:)
    integer(c_int), pointer :: fb(:)
    call c_f_pointer(a, fa, [n])
    call c_f_pointer(b, fb, [n])
    call SortD(n, fa, fb)
  end subroutine usf_sortd

  subroutine usf_dgeev(jobvr, sig, wr, wi, vr, ldvr, work, infor) &
      bind(c, name="usf_dgeev")
    integer(c_int), value :: jobvr, ldvr
    real(c_double) :: sig(3,3), wr(3), wi(3), vr(ldvr,*), work(24)
    integer(c_int) :: infor
    character :: jv
    real(c_double) :: dumy(1)
    if (jobvr == 0) then
      jv = 'N'
    else
      jv = 'V'
    end if
    call DGEEV('N', jv, 3, sig, 3, wr, wi, dumy, 1, vr, ldvr, work, 24, &
        infor)
  end subroutine usf_dgeev

  subroutine usf_warn(caller, ncaller, msg, nmsg) bind(c, name="usf_warn")
    type(c_ptr), value :: caller, msg
    integer(c_int), value :: ncaller, nmsg
    character(kind=c_char), pointer :: cbuf(:)
    character(len=:), allocatable :: fcaller, fmsg
    integer :: i
    call c_f_pointer(caller, cbuf, [ncaller])
    allocate(character(len=ncaller) :: fcaller)
    do i = 1, ncaller
      fcaller(i:i) = cbuf(i)
    end do
    call c_f_pointer(msg, cbuf, [nmsg])
    allocate(character(len=nmsg) :: fmsg)
    do i = 1, nmsg
      fmsg(i:i) = cbuf(i)
    end do
    call Warn(fcaller, fmsg)
  end subroutine usf_warn

  subroutine usf_var_set_at(h_var, node, value) &
      bind(c, name="usf_var_set_at")
    type(c_ptr), value :: h_var
    integer(c_int), value :: node
    real(c_double), value :: value
    type(variable_t), pointer :: var
    call c_f_pointer(h_var, var)
    var % values(var % perm(node)) = value
  end subroutine usf_var_set_at

  function usf_dot(a, b, n) result(s) bind(c, name="usf_dot")
    type(c_ptr), value :: a
    type(c_ptr), value :: b
    integer(c_int), value :: n
    real(c_double) :: s
    real(c_double), pointer :: fa(:)
    real(c_double), pointer :: fb(:)
    call c_f_pointer(a, fa, [n])
    call c_f_pointer(b, fb, [n])
    s = sum(fa(1:n) * fb(1:n))
  end function usf_dot

  function usf_sum_pow2(a, n) result(s) bind(c, name="usf_sum_pow2")
    type(c_ptr), value :: a
    integer(c_int), value :: n
    real(c_double) :: s
    real(c_double), pointer :: fa(:)
    call c_f_pointer(a, fa, [n])
    s = sum(fa(1:n) ** 2.0_dp)
  end function usf_sum_pow2

  subroutine usf_get_equation_of(h_elem, h_eq, found) &
      bind(c, name="usf_get_equation_of")
    type(c_ptr), value :: h_elem
    type(c_ptr) :: h_eq
    integer(c_int) :: found
    type(element_t), pointer :: e
    type(valuelist_t), pointer :: q
    logical :: f
    call c_f_pointer(h_elem, e)
    q => GetEquation(Element=e, Found=f)
    found = merge(1, 0, f)
    if (associated(q)) then
      h_eq = c_loc(q)
    else
      h_eq = c_null_ptr
    end if
  end subroutine usf_get_equation_of

  subroutine usf_boundary_parent(h_elem, h_parent) &
      bind(c, name="usf_boundary_parent")
    type(c_ptr), value :: h_elem
    type(c_ptr) :: h_parent
    type(element_t), pointer :: be
    type(element_t), pointer :: pe
    integer :: other
    call c_f_pointer(h_elem, be)
    other = be % boundaryinfo % outbody
    if (other < 1) then
      pe => be % boundaryinfo % Right
      if (.not. associated(pe)) pe => be % boundaryinfo % Left
    else
      pe => be % boundaryinfo % Right
      if (pe % bodyid == other) pe => be % boundaryinfo % Left
    end if
    if (associated(pe)) then
      h_parent = c_loc(pe)
    else
      h_parent = c_null_ptr
    end if
  end subroutine usf_boundary_parent

  subroutine usf_write_msg_equation_default(buf, nbuf, outlen) &
      bind(c, name="usf_write_msg_equation_default")
    type(c_ptr), value :: buf
    integer(c_int), value :: nbuf
    integer(c_int) :: outlen
    character(kind=c_char), pointer :: cout(:)
    character(len=:), allocatable :: msg
    integer :: m
    integer :: i
    allocate(character(len=nbuf) :: msg)
    write(msg, '(A,I3)') 'No "Equation" found. Using default values for variables'
    m = min(len_trim(msg), nbuf)
    if (m > 0) then
      call c_f_pointer(buf, cout, [m])
      do i = 1, m
        cout(i) = msg(i:i)
      end do
    end if
    outlen = len_trim(msg)
  end subroutine usf_write_msg_equation_default

  subroutine usf_current_solver_values(h_list) &
      bind(c, name="usf_current_solver_values")
    type(c_ptr) :: h_list
    type(valuelist_t), pointer :: l
    l => acc_model % solver % values
    if (associated(l)) then
      h_list = c_loc(l)
    else
      h_list = c_null_ptr
    end if
  end subroutine usf_current_solver_values

  subroutine usf_get_real(h_list, key, nkey, n, found, out) &
      bind(c, name="usf_get_real")
    type(c_ptr), value :: h_list
    type(c_ptr), value :: key
    type(c_ptr), value :: out
    integer(c_int), value :: nkey
    integer(c_int), value :: n
    integer(c_int) :: found
    type(valuelist_t), pointer :: l
    character(len=:), allocatable :: fkey
    character(kind=c_char), pointer :: cbuf(:)
    real(c_double), pointer :: fx(:)
    real(c_double), pointer :: fout(:)
    logical :: f
    integer :: i
    call c_f_pointer(h_list, l)
    call c_f_pointer(key, cbuf, [nkey])
    allocate(character(len=nkey) :: fkey)
    do i = 1, nkey
      fkey(i:i) = cbuf(i)
    end do
    call c_f_pointer(out, fout, [n])
    fx => GetReal(l, fkey, f)
    fout(1:n) = fx(1:n)
    found = merge(1, 0, f)
  end subroutine usf_get_real

  function usf_element_body_id(h_elem) result(id) &
      bind(c, name="usf_element_body_id")
    type(c_ptr), value :: h_elem
    integer(c_int) :: id
    type(element_t), pointer :: e
    call c_f_pointer(h_elem, e)
    id = e % bodyid
  end function usf_element_body_id

  subroutine usf_model_body_values(body_id, h_list) &
      bind(c, name="usf_model_body_values")
    integer(c_int), value :: body_id
    type(c_ptr) :: h_list
    type(valuelist_t), pointer :: l
    l => acc_model % bodies(body_id) % values
    if (associated(l)) then
      h_list = c_loc(l)
    else
      h_list = c_null_ptr
    end if
  end subroutine usf_model_body_values

  subroutine usf_model_material_values(material_id, h_list) &
      bind(c, name="usf_model_material_values")
    integer(c_int), value :: material_id
    type(c_ptr) :: h_list
    type(valuelist_t), pointer :: l
    l => acc_model % materials(material_id) % values
    if (associated(l)) then
      h_list = c_loc(l)
    else
      h_list = c_null_ptr
    end if
  end subroutine usf_model_material_values

  function usf_list_get_integer(h_list, key, nkey, unfoundfatal, found, &
      use_mm, minv, maxv) result(v) bind(c, name="usf_list_get_integer")
    type(c_ptr), value :: h_list
    type(c_ptr), value :: key
    integer(c_int), value :: nkey
    integer(c_int), value :: unfoundfatal
    integer(c_int), value :: use_mm
    integer(c_int), value :: minv
    integer(c_int), value :: maxv
    integer(c_int) :: found
    integer(c_int) :: v
    type(valuelist_t), pointer :: l
    character(len=:), allocatable :: fkey
    character(kind=c_char), pointer :: cbuf(:)
    logical :: uff
    logical :: f
    integer :: i
    call c_f_pointer(h_list, l)
    call c_f_pointer(key, cbuf, [nkey])
    allocate(character(len=nkey) :: fkey)
    do i = 1, nkey
      fkey(i:i) = cbuf(i)
    end do
    uff = (unfoundfatal /= 0)
    if (use_mm /= 0) then
      v = ListGetInteger(l, fkey, f, minv=minv, maxv=maxv, UnfoundFatal=uff)
    else
      v = ListGetInteger(l, fkey, f, UnfoundFatal=uff)
    end if
    found = merge(1, 0, f)
  end function usf_list_get_integer

  function usf_get_creal(h_list, key, nkey, found) result(v) &
      bind(c, name="usf_get_creal")
    type(c_ptr), value :: h_list
    type(c_ptr), value :: key
    integer(c_int), value :: nkey
    integer(c_int) :: found
    real(c_double) :: v
    type(valuelist_t), pointer :: l
    character(len=:), allocatable :: fkey
    character(kind=c_char), pointer :: cbuf(:)
    logical :: f
    integer :: i
    call c_f_pointer(h_list, l)
    call c_f_pointer(key, cbuf, [nkey])
    allocate(character(len=nkey) :: fkey)
    do i = 1, nkey
      fkey(i:i) = cbuf(i)
    end do
    v = GetCReal(l, fkey, f)
    found = merge(1, 0, f)
  end function usf_get_creal

  function usf_nonlin_change() result(r) bind(c, name="usf_nonlin_change")
    real(c_double) :: r
    r = acc_model % solver % variable % nonlinchange
  end function usf_nonlin_change

  subroutine usf_model_mesh_variables(h_vars) &
      bind(c, name="usf_model_mesh_variables")
    type(c_ptr) :: h_vars
    type(variable_t), pointer :: v
    v => acc_model % mesh % variables
    if (associated(v)) then
      h_vars = c_loc(v)
    else
      h_vars = c_null_ptr
    end if
  end subroutine usf_model_mesh_variables

  function usf_mesh_number_of_nodes() result(n) &
      bind(c, name="usf_mesh_number_of_nodes")
    integer(c_int) :: n
    n = acc_model % mesh % numberofnodes
  end function usf_mesh_number_of_nodes

  function usf_number_of_boundary_elements() result(n) &
      bind(c, name="usf_number_of_boundary_elements")
    integer(c_int) :: n
    n = acc_model % numberofboundaryelements
  end function usf_number_of_boundary_elements

  function usf_mype() result(p) bind(c, name="usf_mype")
    use parallelutils
    integer(c_int) :: p
    p = ParEnv % MyPE
  end function usf_mype

  subroutine usf_get_boundary_element(t, h_elem) &
      bind(c, name="usf_get_boundary_element")
    integer(c_int), value :: t
    type(c_ptr) :: h_elem
    type(element_t), pointer :: e
    e => GetBoundaryElement(t)
    if (associated(e)) then
      h_elem = c_loc(e)
    else
      h_elem = c_null_ptr
    end if
  end subroutine usf_get_boundary_element

  function usf_element_part_index(h_elem) result(p) &
      bind(c, name="usf_element_part_index")
    type(c_ptr), value :: h_elem
    integer(c_int) :: p
    type(element_t), pointer :: e
    call c_f_pointer(h_elem, e)
    p = e % partindex
  end function usf_element_part_index

  function usf_list_get_string(h_list, key, nkey, buf, nbuf, outlen, found) &
      result(dummy) bind(c, name="usf_list_get_string")
    type(c_ptr), value :: h_list
    type(c_ptr), value :: key
    integer(c_int), value :: nkey
    type(c_ptr), value :: buf
    integer(c_int), value :: nbuf
    integer(c_int) :: outlen
    integer(c_int) :: found
    integer(c_int) :: dummy
    type(valuelist_t), pointer :: l
    character(len=:), allocatable :: fkey
    character(len=:), allocatable :: val
    character(kind=c_char), pointer :: cbuf(:)
    character(kind=c_char), pointer :: cout(:)
    logical :: f
    integer :: i
    integer :: m
    dummy = 0
    call c_f_pointer(h_list, l)
    call c_f_pointer(key, cbuf, [nkey])
    allocate(character(len=nkey) :: fkey)
    do i = 1, nkey
      fkey(i:i) = cbuf(i)
    end do
    val = ListGetString(l, fkey, f)
    found = merge(1, 0, f)
    m = min(len(val), nbuf)
    if (m > 0) then
      call c_f_pointer(buf, cout, [m])
      do i = 1, m
        cout(i) = val(i:i)
      end do
    end if
    outlen = len(val)
  end function usf_list_get_string

  function usf_list_get_const_real(h_list, key, nkey) result(v) &
      bind(c, name="usf_list_get_const_real")
    type(c_ptr), value :: h_list
    type(c_ptr), value :: key
    integer(c_int), value :: nkey
    real(c_double) :: v
    type(valuelist_t), pointer :: l
    character(len=:), allocatable :: fkey
    character(kind=c_char), pointer :: cbuf(:)
    logical, parameter :: uff = .TRUE.
    integer :: i
    call c_f_pointer(h_list, l)
    call c_f_pointer(key, cbuf, [nkey])
    allocate(character(len=nkey) :: fkey)
    do i = 1, nkey
      fkey(i:i) = cbuf(i)
    end do
    v = ListGetConstReal(l, fkey, UnfoundFatal=uff)
  end function usf_list_get_const_real

  subroutine usf_set_current_element(h_elem) &
      bind(c, name="usf_set_current_element")
    type(c_ptr), value :: h_elem
    type(element_t), pointer :: e
    call c_f_pointer(h_elem, e)
    acc_model % currentelement => e
  end subroutine usf_set_current_element

  subroutine usf_element_nodes_xyz(h_elem, n, x, y, z) &
      bind(c, name="usf_element_nodes_xyz")
    type(c_ptr), value :: h_elem
    integer(c_int), value :: n
    type(c_ptr), value :: x
    type(c_ptr), value :: y
    type(c_ptr), value :: z
    type(element_t), pointer :: e
    real(c_double), pointer :: fx(:)
    real(c_double), pointer :: fy(:)
    real(c_double), pointer :: fz(:)
    type(nodes_t), save, target :: acc_nodes
    call c_f_pointer(h_elem, e)
    call GetElementNodes(acc_nodes, UElement=e)
    call c_f_pointer(x, fx, [n])
    call c_f_pointer(y, fy, [n])
    call c_f_pointer(z, fz, [n])
    fx(1:n) = acc_nodes % x(1:n)
    fy(1:n) = acc_nodes % y(1:n)
    fz(1:n) = acc_nodes % z(1:n)
  end subroutine usf_element_nodes_xyz

  subroutine usf_write_msg_star_aa(pre, npre, name, nname, buf, nbuf, &
      outlen) bind(c, name="usf_write_msg_star_aa")
    type(c_ptr), value :: pre
    type(c_ptr), value :: name
    type(c_ptr), value :: buf
    integer(c_int), value :: npre
    integer(c_int), value :: nname
    integer(c_int), value :: nbuf
    integer(c_int) :: outlen
    character(kind=c_char), pointer :: cpre(:)
    character(kind=c_char), pointer :: cname(:)
    character(kind=c_char), pointer :: cout(:)
    character(len=:), allocatable :: fpre
    character(len=:), allocatable :: fname
    character(len=:), allocatable :: msg
    integer :: m
    integer :: i
    call c_f_pointer(pre, cpre, [npre])
    allocate(character(len=npre) :: fpre)
    do i = 1, npre
      fpre(i:i) = cpre(i)
    end do
    call c_f_pointer(name, cname, [nname])
    allocate(character(len=nname) :: fname)
    do i = 1, nname
      fname(i:i) = cname(i)
    end do
    allocate(character(len=nbuf) :: msg)
    write(msg, *) fpre, fname
    m = min(len_trim(msg), nbuf)
    if (m > 0) then
      call c_f_pointer(buf, cout, [m])
      do i = 1, m
        cout(i) = msg(i:i)
      end do
    end if
    outlen = len_trim(msg)
  end subroutine usf_write_msg_star_aa

  subroutine usf_write_msg_star_ar(pre, npre, value, buf, nbuf, outlen) &
      bind(c, name="usf_write_msg_star_ar")
    type(c_ptr), value :: pre
    type(c_ptr), value :: buf
    integer(c_int), value :: npre
    integer(c_int), value :: nbuf
    real(c_double), value :: value
    integer(c_int) :: outlen
    character(kind=c_char), pointer :: cpre(:)
    character(kind=c_char), pointer :: cout(:)
    character(len=:), allocatable :: fpre
    character(len=:), allocatable :: msg
    integer :: m
    integer :: i
    call c_f_pointer(pre, cpre, [npre])
    allocate(character(len=npre) :: fpre)
    do i = 1, npre
      fpre(i:i) = cpre(i)
    end do
    allocate(character(len=nbuf) :: msg)
    write(msg, *) fpre, value
    m = min(len_trim(msg), nbuf)
    if (m > 0) then
      call c_f_pointer(buf, cout, [m])
      do i = 1, m
        cout(i) = msg(i:i)
      end do
    end if
    outlen = len_trim(msg)
  end subroutine usf_write_msg_star_ar

  subroutine usf_write_msg_star_arr(pre, npre, x, y, buf, nbuf, outlen) &
      bind(c, name="usf_write_msg_star_arr")
    type(c_ptr), value :: pre
    type(c_ptr), value :: buf
    integer(c_int), value :: npre
    integer(c_int), value :: nbuf
    real(c_double), value :: x
    real(c_double), value :: y
    integer(c_int) :: outlen
    character(kind=c_char), pointer :: cpre(:)
    character(kind=c_char), pointer :: cout(:)
    character(len=:), allocatable :: fpre
    character(len=:), allocatable :: msg
    integer :: m
    integer :: i
    call c_f_pointer(pre, cpre, [npre])
    allocate(character(len=npre) :: fpre)
    do i = 1, npre
      fpre(i:i) = cpre(i)
    end do
    allocate(character(len=nbuf) :: msg)
    write(msg, *) fpre, x, y
    m = min(len_trim(msg), nbuf)
    if (m > 0) then
      call c_f_pointer(buf, cout, [m])
      do i = 1, m
        cout(i) = msg(i:i)
      end do
    end if
    outlen = len_trim(msg)
  end subroutine usf_write_msg_star_arr

  subroutine usf_write_msg_star_s(text, ntext, buf, nbuf, outlen) &
      bind(c, name="usf_write_msg_star_s")
    type(c_ptr), value :: text
    type(c_ptr), value :: buf
    integer(c_int), value :: ntext
    integer(c_int), value :: nbuf
    integer(c_int) :: outlen
    character(kind=c_char), pointer :: ctext(:)
    character(kind=c_char), pointer :: cout(:)
    character(len=:), allocatable :: ftext
    character(len=:), allocatable :: msg
    integer :: m
    integer :: i
    call c_f_pointer(text, ctext, [ntext])
    allocate(character(len=ntext) :: ftext)
    do i = 1, ntext
      ftext(i:i) = ctext(i)
    end do
    allocate(character(len=nbuf) :: msg)
    write(msg, *) ftext
    m = min(len_trim(msg), nbuf)
    if (m > 0) then
      call c_f_pointer(buf, cout, [m])
      do i = 1, m
        cout(i) = msg(i:i)
      end do
    end if
    outlen = len_trim(msg)
  end subroutine usf_write_msg_star_s

  subroutine usf_write_msg_star_asal(p1, np1, name, nname, p2, np2, lv, &
      buf, nbuf, outlen) bind(c, name="usf_write_msg_star_asal")
    type(c_ptr), value :: p1
    type(c_ptr), value :: name
    type(c_ptr), value :: p2
    type(c_ptr), value :: buf
    integer(c_int), value :: np1
    integer(c_int), value :: nname
    integer(c_int), value :: np2
    integer(c_int), value :: nbuf
    integer(c_int), value :: lv
    integer(c_int) :: outlen
    character(kind=c_char), pointer :: cp1(:)
    character(kind=c_char), pointer :: cname(:)
    character(kind=c_char), pointer :: cp2(:)
    character(kind=c_char), pointer :: cout(:)
    character(len=:), allocatable :: fp1
    character(len=:), allocatable :: fname
    character(len=:), allocatable :: fp2
    character(len=:), allocatable :: msg
    logical :: v
    integer :: m
    integer :: i
    call c_f_pointer(p1, cp1, [np1])
    allocate(character(len=np1) :: fp1)
    do i = 1, np1
      fp1(i:i) = cp1(i)
    end do
    call c_f_pointer(name, cname, [nname])
    allocate(character(len=nname) :: fname)
    do i = 1, nname
      fname(i:i) = cname(i)
    end do
    call c_f_pointer(p2, cp2, [np2])
    allocate(character(len=np2) :: fp2)
    do i = 1, np2
      fp2(i:i) = cp2(i)
    end do
    v = (lv /= 0)
    allocate(character(len=nbuf) :: msg)
    write(msg, *) fp1, fname, fp2, v
    m = min(len_trim(msg), nbuf)
    if (m > 0) then
      call c_f_pointer(buf, cout, [m])
      do i = 1, m
        cout(i) = msg(i:i)
      end do
    end if
    outlen = len_trim(msg)
  end subroutine usf_write_msg_star_asal

  subroutine usf_write_msg_star_air(pre, npre, ival, rval, buf, nbuf, &
      outlen) bind(c, name="usf_write_msg_star_air")
    type(c_ptr), value :: pre
    type(c_ptr), value :: buf
    integer(c_int), value :: npre
    integer(c_int), value :: nbuf
    integer(c_int), value :: ival
    real(c_double), value :: rval
    integer(c_int) :: outlen
    character(kind=c_char), pointer :: cpre(:)
    character(kind=c_char), pointer :: cout(:)
    character(len=:), allocatable :: fpre
    character(len=:), allocatable :: msg
    integer :: m
    integer :: i
    call c_f_pointer(pre, cpre, [npre])
    allocate(character(len=npre) :: fpre)
    do i = 1, npre
      fpre(i:i) = cpre(i)
    end do
    allocate(character(len=nbuf) :: msg)
    write(msg, *) fpre, ival, rval
    m = min(len_trim(msg), nbuf)
    if (m > 0) then
      call c_f_pointer(buf, cout, [m])
      do i = 1, m
        cout(i) = msg(i:i)
      end do
    end if
    outlen = len_trim(msg)
  end subroutine usf_write_msg_star_air



  subroutine usf_write_msg_ae148(pre, npre, value, buf, nbuf, outlen) &
      bind(c, name="usf_write_msg_ae148")
    type(c_ptr), value :: pre
    type(c_ptr), value :: buf
    integer(c_int), value :: npre
    integer(c_int), value :: nbuf
    real(c_double), value :: value
    integer(c_int) :: outlen
    character(kind=c_char), pointer :: cpre(:)
    character(kind=c_char), pointer :: cout(:)
    character(len=:), allocatable :: fpre
    character(len=:), allocatable :: msg
    integer :: m
    integer :: i
    call c_f_pointer(pre, cpre, [npre])
    allocate(character(len=npre) :: fpre)
    do i = 1, npre
      fpre(i:i) = cpre(i)
    end do
    allocate(character(len=nbuf) :: msg)
    write(msg, '(A,e14.8)') fpre, value
    m = min(len_trim(msg), nbuf)
    if (m > 0) then
      call c_f_pointer(buf, cout, [m])
      do i = 1, m
        cout(i) = msg(i:i)
      end do
    end if
    outlen = len_trim(msg)
  end subroutine usf_write_msg_ae148

  function usf_element_type_nodes(h_elem) result(n) &
      bind(c, name="usf_element_type_nodes")
    type(c_ptr), value :: h_elem
    integer(c_int) :: n
    type(element_t), pointer :: e
    call c_f_pointer(h_elem, e)
    n = e % type % numberofnodes
  end function usf_element_type_nodes

  subroutine usf_print_ice_pressure(value) &
      bind(c, name="usf_print_ice_pressure")
    real(c_double), value :: value
    print *, "Ice pressure value :", value
  end subroutine usf_print_ice_pressure

  function usf_get_element_family(h_elem) result(f) &
      bind(c, name="usf_get_element_family")
    type(c_ptr), value :: h_elem
    integer(c_int) :: f
    type(element_t), pointer :: e
    call c_f_pointer(h_elem, e)
    f = GetElementFamily(e)
  end function usf_get_element_family

  function usf_solver_dt() result(dt) bind(c, name="usf_solver_dt")
    real(c_double) :: dt
    dt = acc_model % solver % dt
  end function usf_solver_dt

  subroutine usf_current_solver(h_solver) &
      bind(c, name="usf_current_solver")
    type(c_ptr) :: h_solver
    type(solver_t), pointer :: s
    s => acc_model % solver
    if (associated(s)) then
      h_solver = c_loc(s)
    else
      h_solver = c_null_ptr
    end if
  end subroutine usf_current_solver

  function usf_mesh_tag() result(t) bind(c, name="usf_mesh_tag")
    integer(c_int) :: t
    t = acc_model % mesh % meshtag
  end function usf_mesh_tag

  function usf_mesh_changed() result(c) bind(c, name="usf_mesh_changed")
    integer(c_int) :: c
    c = merge(1, 0, acc_model % mesh % changed)
  end function usf_mesh_changed

  function usf_number_of_bodyforces() result(n) &
      bind(c, name="usf_number_of_bodyforces")
    integer(c_int) :: n
    n = acc_model % numberofbodyforces
  end function usf_number_of_bodyforces

  function usf_number_of_materials() result(n) &
      bind(c, name="usf_number_of_materials")
    integer(c_int) :: n
    n = acc_model % numberofmaterials
  end function usf_number_of_materials

  function usf_list_check_present(h_list, key, nkey) result(p) &
      bind(c, name="usf_list_check_present")
    type(c_ptr), value :: h_list
    type(c_ptr), value :: key
    integer(c_int), value :: nkey
    integer(c_int) :: p
    type(valuelist_t), pointer :: l
    character(len=:), allocatable :: fkey
    character(kind=c_char), pointer :: cbuf(:)
    integer :: i
    call c_f_pointer(h_list, l)
    call c_f_pointer(key, cbuf, [nkey])
    allocate(character(len=nkey) :: fkey)
    do i = 1, nkey
      fkey(i:i) = cbuf(i)
    end do
    p = merge(1, 0, ListCheckPresent(l, fkey))
  end function usf_list_check_present

  function usf_normal_vector(h_elem, bu, bv, out) result(dummy) &
      bind(c, name="usf_normal_vector")
    type(c_ptr), value :: h_elem
    real(c_double), value :: bu
    real(c_double), value :: bv
    real(c_double) :: out(3)
    integer(c_int) :: dummy
    type(element_t), pointer :: e
    real(c_double) :: nv(3)
    type(nodes_t), save, target :: acc_nodes_nv
    dummy = 0
    call c_f_pointer(h_elem, e)
    call GetElementNodes(acc_nodes_nv, UElement=e)
    nv = NormalVector(e, acc_nodes_nv, bu, bv, .TRUE.)
    out = nv
  end function usf_normal_vector

  function usf_consistent_normal_vector(h_solver, h_var, h_elem, node, &
      out, got) result(dummy) bind(c, name="usf_consistent_normal_vector")
    type(c_ptr), value :: h_solver
    type(c_ptr), value :: h_var
    type(c_ptr), value :: h_elem
    integer(c_int), value :: node
    real(c_double) :: out(3)
    integer(c_int) :: got
    integer(c_int) :: dummy
    type(solver_t), pointer :: s
    type(variable_t), pointer :: nv
    type(element_t), pointer :: e
    real(c_double) :: nrm(3)
    logical :: f
    dummy = 0
    call c_f_pointer(h_solver, s)
    call c_f_pointer(h_var, nv)
    call c_f_pointer(h_elem, e)
    nrm = ConsistentNormalVector(s, nv, e, f, node=node)
    out = nrm
    got = merge(1, 0, f)
  end function usf_consistent_normal_vector

  function usf_element_type_nodeu(h_elem, i) result(u) &
      bind(c, name="usf_element_type_nodeu")
    type(c_ptr), value :: h_elem
    integer(c_int), value :: i
    real(c_double) :: u
    type(element_t), pointer :: e
    call c_f_pointer(h_elem, e)
    u = e % type % nodeu(i)
  end function usf_element_type_nodeu

  function usf_element_type_nodev(h_elem, i) result(v) &
      bind(c, name="usf_element_type_nodev")
    type(c_ptr), value :: h_elem
    integer(c_int), value :: i
    real(c_double) :: v
    type(element_t), pointer :: e
    call c_f_pointer(h_elem, e)
    v = e % type % nodev(i)
  end function usf_element_type_nodev

  function usf_element_type_dim(h_elem) result(d) &
      bind(c, name="usf_element_type_dim")
    type(c_ptr), value :: h_elem
    integer(c_int) :: d
    type(element_t), pointer :: e
    call c_f_pointer(h_elem, e)
    d = e % type % dimension
  end function usf_element_type_dim

  function usf_model_node_y(node) result(y) &
      bind(c, name="usf_model_node_y")
    integer(c_int), value :: node
    real(c_double) :: y
    y = acc_model % nodes % y(node)
  end function usf_model_node_y

  function usf_model_node_z(node) result(z) &
      bind(c, name="usf_model_node_z")
    integer(c_int), value :: node
    real(c_double) :: z
    z = acc_model % nodes % z(node)
  end function usf_model_node_z

  subroutine usf_write_msg_star_al(pre, npre, v, buf, nbuf, outlen) &
      bind(c, name="usf_write_msg_star_al")
    type(c_ptr), value :: pre
    type(c_ptr), value :: buf
    integer(c_int), value :: npre
    integer(c_int), value :: nbuf
    integer(c_int), value :: v
    integer(c_int) :: outlen
    character(kind=c_char), pointer :: cpre(:)
    character(kind=c_char), pointer :: cout(:)
    character(len=:), allocatable :: fpre
    character(len=:), allocatable :: msg
    logical :: lv
    integer :: m
    integer :: i
    call c_f_pointer(pre, cpre, [npre])
    allocate(character(len=npre) :: fpre)
    do i = 1, npre
      fpre(i:i) = cpre(i)
    end do
    lv = (v /= 0)
    allocate(character(len=nbuf) :: msg)
    write(msg, *) fpre, lv
    m = min(len_trim(msg), nbuf)
    if (m > 0) then
      call c_f_pointer(buf, cout, [m])
      do i = 1, m
        cout(i) = msg(i:i)
      end do
    end if
    outlen = len_trim(msg)
  end subroutine usf_write_msg_star_al

  function usf_variable_dofs(h_var) result(d) bind(c, name="usf_variable_dofs")
    type(c_ptr), value :: h_var
    integer(c_int) :: d
    type(variable_t), pointer :: var
    call c_f_pointer(h_var, var)
    d = var % dofs
  end function usf_variable_dofs

  function usf_element_type_nodew(h_elem, i) result(w) &
      bind(c, name="usf_element_type_nodew")
    type(c_ptr), value :: h_elem
    integer(c_int), value :: i
    real(c_double) :: w
    type(element_t), pointer :: e
    call c_f_pointer(h_elem, e)
    w = e % type % nodew(i)
  end function usf_element_type_nodew

  subroutine usf_element_info(h_elem, u, v, w, detj, basis, dbasis, n, &
      stat) bind(c, name="usf_element_info")
    type(c_ptr), value :: h_elem
    real(c_double), value :: u
    real(c_double), value :: v
    real(c_double), value :: w
    real(c_double) :: detj
    integer(c_int), value :: n
    real(c_double) :: basis(n)
    real(c_double) :: dbasis(n, 3)
    integer(c_int) :: stat
    type(element_t), pointer :: e
    real(c_double) :: ddb(1, 1, 1)
    type(nodes_t), save, target :: acc_nodes_ei
    logical :: ok
    call c_f_pointer(h_elem, e)
    call GetElementNodes(acc_nodes_ei, UElement=e)
    ok = ElementInfo(e, acc_nodes_ei, u, v, w, detj, basis, dbasis, ddb, &
        .FALSE., .FALSE.)
    stat = merge(1, 0, ok)
  end subroutine usf_element_info

  subroutine usf_lgrad_dij(velo, dbasis, n, dij) &
      bind(c, name="usf_lgrad_dij")
    type(c_ptr), value :: velo
    type(c_ptr), value :: dbasis
    type(c_ptr), value :: dij
    integer(c_int), value :: n
    real(c_double), pointer :: fv(:, :)
    real(c_double), pointer :: fd(:, :)
    real(c_double), pointer :: fo(:, :)
    real(c_double) :: lg(3, 3)
    call c_f_pointer(velo, fv, [3, n])
    call c_f_pointer(dbasis, fd, [n, 3])
    call c_f_pointer(dij, fo, [3, 3])
    lg = matmul(fv(:, 1:n), fd(1:n, :))
    fo = 0.5_dp * (lg + transpose(lg))
  end subroutine usf_lgrad_dij

  subroutine usf_write_msg_aea104a(pre, npre, value, suf, nsuf, buf, &
      nbuf, outlen) bind(c, name="usf_write_msg_aea104a")
    type(c_ptr), value :: pre
    type(c_ptr), value :: suf
    type(c_ptr), value :: buf
    integer(c_int), value :: npre
    integer(c_int), value :: nsuf
    integer(c_int), value :: nbuf
    real(c_double), value :: value
    integer(c_int) :: outlen
    character(kind=c_char), pointer :: cpre(:)
    character(kind=c_char), pointer :: csuf(:)
    character(kind=c_char), pointer :: cout(:)
    character(len=:), allocatable :: fpre
    character(len=:), allocatable :: fsuf
    character(len=:), allocatable :: msg
    integer :: m
    integer :: i
    call c_f_pointer(pre, cpre, [npre])
    allocate(character(len=npre) :: fpre)
    do i = 1, npre
      fpre(i:i) = cpre(i)
    end do
    call c_f_pointer(suf, csuf, [nsuf])
    allocate(character(len=nsuf) :: fsuf)
    do i = 1, nsuf
      fsuf(i:i) = csuf(i)
    end do
    allocate(character(len=nbuf) :: msg)
    write(msg, '(a,e10.4,a)') fpre, value, fsuf
    m = min(len_trim(msg), nbuf)
    if (m > 0) then
      call c_f_pointer(buf, cout, [m])
      do i = 1, m
        cout(i) = msg(i:i)
      end do
    end if
    outlen = len_trim(msg)
  end subroutine usf_write_msg_aea104a

  subroutine usf_write_msg_ai2ai2a(p1, n1, i1, p2, n2, i2, p3, n3, buf, &
      nbuf, outlen) bind(c, name="usf_write_msg_ai2ai2a")
    type(c_ptr), value :: p1
    type(c_ptr), value :: p2
    type(c_ptr), value :: p3
    type(c_ptr), value :: buf
    integer(c_int), value :: n1
    integer(c_int), value :: n2
    integer(c_int), value :: n3
    integer(c_int), value :: nbuf
    integer(c_int), value :: i1
    integer(c_int), value :: i2
    integer(c_int) :: outlen
    character(kind=c_char), pointer :: c1(:)
    character(kind=c_char), pointer :: c2(:)
    character(kind=c_char), pointer :: c3(:)
    character(kind=c_char), pointer :: cout(:)
    character(len=:), allocatable :: f1
    character(len=:), allocatable :: f2
    character(len=:), allocatable :: f3
    character(len=:), allocatable :: msg
    integer :: m
    integer :: i
    call c_f_pointer(p1, c1, [n1])
    allocate(character(len=n1) :: f1)
    do i = 1, n1
      f1(i:i) = c1(i)
    end do
    call c_f_pointer(p2, c2, [n2])
    allocate(character(len=n2) :: f2)
    do i = 1, n2
      f2(i:i) = c2(i)
    end do
    call c_f_pointer(p3, c3, [n3])
    allocate(character(len=n3) :: f3)
    do i = 1, n3
      f3(i:i) = c3(i)
    end do
    allocate(character(len=nbuf) :: msg)
    write(msg, '(a,I2,a,I2,a)') f1, i1, f2, i2, f3
    m = min(len_trim(msg), nbuf)
    if (m > 0) then
      call c_f_pointer(buf, cout, [m])
      do i = 1, m
        cout(i) = msg(i:i)
      end do
    end if
    outlen = len_trim(msg)
  end subroutine usf_write_msg_ai2ai2a

  subroutine usf_write_msg_ai2ai2(p1, n1, i1, p2, n2, i2, buf, nbuf, &
      outlen) bind(c, name="usf_write_msg_ai2ai2")
    type(c_ptr), value :: p1
    type(c_ptr), value :: p2
    type(c_ptr), value :: buf
    integer(c_int), value :: n1
    integer(c_int), value :: n2
    integer(c_int), value :: nbuf
    integer(c_int), value :: i1
    integer(c_int), value :: i2
    integer(c_int) :: outlen
    character(kind=c_char), pointer :: c1(:)
    character(kind=c_char), pointer :: c2(:)
    character(kind=c_char), pointer :: cout(:)
    character(len=:), allocatable :: f1
    character(len=:), allocatable :: f2
    character(len=:), allocatable :: msg
    integer :: m
    integer :: i
    call c_f_pointer(p1, c1, [n1])
    allocate(character(len=n1) :: f1)
    do i = 1, n1
      f1(i:i) = c1(i)
    end do
    call c_f_pointer(p2, c2, [n2])
    allocate(character(len=n2) :: f2)
    do i = 1, n2
      f2(i:i) = c2(i)
    end do
    allocate(character(len=nbuf) :: msg)
    write(msg, '(a,I2,a,I2)') f1, i1, f2, i2
    m = min(len_trim(msg), nbuf)
    if (m > 0) then
      call c_f_pointer(buf, cout, [m])
      do i = 1, m
        cout(i) = msg(i:i)
      end do
    end if
    outlen = len_trim(msg)
  end subroutine usf_write_msg_ai2ai2

  subroutine usf_write_msg_aaai3a(p1, n1, p2, n2, p3, n3, ival, p4, n4, &
      buf, nbuf, outlen) bind(c, name="usf_write_msg_aaai3a")
    type(c_ptr), value :: p1
    type(c_ptr), value :: p2
    type(c_ptr), value :: p3
    type(c_ptr), value :: p4
    type(c_ptr), value :: buf
    integer(c_int), value :: n1
    integer(c_int), value :: n2
    integer(c_int), value :: n3
    integer(c_int), value :: n4
    integer(c_int), value :: nbuf
    integer(c_int), value :: ival
    integer(c_int) :: outlen
    character(kind=c_char), pointer :: c1(:)
    character(kind=c_char), pointer :: c2(:)
    character(kind=c_char), pointer :: c3(:)
    character(kind=c_char), pointer :: c4(:)
    character(kind=c_char), pointer :: cout(:)
    character(len=:), allocatable :: f1
    character(len=:), allocatable :: f2
    character(len=:), allocatable :: f3
    character(len=:), allocatable :: f4
    character(len=:), allocatable :: msg
    integer :: m
    integer :: i
    call c_f_pointer(p1, c1, [n1])
    allocate(character(len=n1) :: f1)
    do i = 1, n1
      f1(i:i) = c1(i)
    end do
    call c_f_pointer(p2, c2, [n2])
    allocate(character(len=n2) :: f2)
    do i = 1, n2
      f2(i:i) = c2(i)
    end do
    call c_f_pointer(p3, c3, [n3])
    allocate(character(len=n3) :: f3)
    do i = 1, n3
      f3(i:i) = c3(i)
    end do
    call c_f_pointer(p4, c4, [n4])
    allocate(character(len=n4) :: f4)
    do i = 1, n4
      f4(i:i) = c4(i)
    end do
    allocate(character(len=nbuf) :: msg)
    write(msg, '(A,A,A,i3,A)') f1, f2, f3, ival, f4
    m = min(len_trim(msg), nbuf)
    if (m > 0) then
      call c_f_pointer(buf, cout, [m])
      do i = 1, m
        cout(i) = msg(i:i)
      end do
    end if
    outlen = len_trim(msg)
  end subroutine usf_write_msg_aaai3a

  subroutine usf_list_get_real_array(h_list, key, nkey, n, idx, d1, d2, &
      got, out) bind(c, name="usf_list_get_real_array")
    type(c_ptr), value :: h_list
    type(c_ptr), value :: key
    type(c_ptr), value :: idx
    type(c_ptr), value :: out
    integer(c_int), value :: nkey
    integer(c_int), value :: n
    integer(c_int) :: d1
    integer(c_int) :: d2
    integer(c_int) :: got
    type(valuelist_t), pointer :: l
    character(len=:), allocatable :: fkey
    character(kind=c_char), pointer :: cbuf(:)
    integer(c_int), pointer :: fidx(:)
    real(c_double), pointer :: fout(:, :)
    real(c_double), pointer :: f(:, :, :) => null()
    logical :: g
    integer :: i
    integer :: j
    integer :: m
    call c_f_pointer(h_list, l)
    call c_f_pointer(key, cbuf, [nkey])
    allocate(character(len=nkey) :: fkey)
    do i = 1, nkey
      fkey(i:i) = cbuf(i)
    end do
    call c_f_pointer(idx, fidx, [n])
    call c_f_pointer(out, fout, [3, n])
    call ListGetRealArray(l, fkey, f, n, fidx, g)
    got = merge(1, 0, g)
    if (got == 1) then
      d1 = size(f, 1)
      d2 = size(f, 2)
      m = min(3, d1)
      do j = 1, m
        fout(j, 1:n) = f(j, 1, 1:n)
      end do
    else
      d1 = 0
      d2 = 0
    end if
    if (associated(f)) deallocate(f)
  end subroutine usf_list_get_real_array

  function usf_list_get_const_real_found(h_list, key, nkey, found) &
      result(v) bind(c, name="usf_list_get_const_real_found")
    type(c_ptr), value :: h_list
    type(c_ptr), value :: key
    integer(c_int), value :: nkey
    integer(c_int) :: found
    real(c_double) :: v
    type(valuelist_t), pointer :: l
    character(len=:), allocatable :: fkey
    character(kind=c_char), pointer :: cbuf(:)
    logical :: f
    integer :: i
    call c_f_pointer(h_list, l)
    call c_f_pointer(key, cbuf, [nkey])
    allocate(character(len=nkey) :: fkey)
    do i = 1, nkey
      fkey(i:i) = cbuf(i)
    end do
    v = ListGetConstReal(l, fkey, f)
    found = merge(1, 0, f)
  end function usf_list_get_const_real_found

  subroutine usf_current_solver_mesh_variables(h_vars) &
      bind(c, name="usf_current_solver_mesh_variables")
    type(c_ptr) :: h_vars
    type(variable_t), pointer :: v
    v => acc_model % solver % mesh % variables
    if (associated(v)) then
      h_vars = c_loc(v)
    else
      h_vars = c_null_ptr
    end if
  end subroutine usf_current_solver_mesh_variables

  function usf_model_node_x(node) result(x) &
      bind(c, name="usf_model_node_x")
    integer(c_int), value :: node
    real(c_double) :: x
    x = acc_model % nodes % x(node)
  end function usf_model_node_x

  subroutine usf_cubic_spline(n, x, y, r) bind(c, name="usf_cubic_spline")
    use generalutils
    integer(c_int), value :: n
    type(c_ptr), value :: x
    type(c_ptr), value :: y
    type(c_ptr), value :: r
    real(c_double), pointer :: fx(:)
    real(c_double), pointer :: fy(:)
    real(c_double), pointer :: fr(:)
    call c_f_pointer(x, fx, [n])
    call c_f_pointer(y, fy, [n])
    call c_f_pointer(r, fr, [n])
    call CubicSpline(n, fx, fy, fr)
  end subroutine usf_cubic_spline

  function usf_interpolate_curve(n, t, f, tval, c) result(v) &
      bind(c, name="usf_interpolate_curve")
    use generalutils
    integer(c_int), value :: n
    type(c_ptr), value :: t
    type(c_ptr), value :: f
    type(c_ptr), value :: c
    real(c_double), value :: tval
    real(c_double) :: v
    real(c_double), pointer :: ft(:)
    real(c_double), pointer :: ff(:)
    real(c_double), pointer :: fc(:)
    call c_f_pointer(t, ft, [n])
    call c_f_pointer(f, ff, [n])
    call c_f_pointer(c, fc, [n])
    v = InterpolateCurve(ft, ff, tval, fc)
  end function usf_interpolate_curve

  subroutine usf_print_newtime() bind(c, name="usf_print_newtime")
    print *, "USF_Zs: New time, saving previous Zs values!"
  end subroutine usf_print_newtime
  function usf_mesh_node_x(node) result(x) bind(c, name="usf_mesh_node_x")
    integer(c_int), value :: node
    real(c_double) :: x
    x = acc_model % mesh % nodes % x(node)
  end function usf_mesh_node_x

  function usf_mesh_node_y(node) result(y) bind(c, name="usf_mesh_node_y")
    integer(c_int), value :: node
    real(c_double) :: y
    y = acc_model % mesh % nodes % y(node)
  end function usf_mesh_node_y

  function usf_mesh_node_z(node) result(z) bind(c, name="usf_mesh_node_z")
    integer(c_int), value :: node
    real(c_double) :: z
    z = acc_model % mesh % nodes % z(node)
  end function usf_mesh_node_z

  subroutine usf_model_mesh(h_mesh) bind(c, name="usf_model_mesh")
    type(c_ptr) :: h_mesh
    type(mesh_t), pointer :: m
    m => acc_model % mesh
    if (associated(m)) then
      h_mesh = c_loc(m)
    else
      h_mesh = c_null_ptr
    end if
  end subroutine usf_model_mesh

  subroutine usf_mesh_element(h_mesh, i, h_elem) &
      bind(c, name="usf_mesh_element")
    type(c_ptr), value :: h_mesh
    integer(c_int), value :: i
    type(c_ptr) :: h_elem
    type(mesh_t), pointer :: m
    type(element_t), pointer :: e
    call c_f_pointer(h_mesh, m)
    e => m % elements(i)
    if (associated(e)) then
      h_elem = c_loc(e)
    else
      h_elem = c_null_ptr
    end if
  end subroutine usf_mesh_element

  function usf_mesh_bulk_elements(h_mesh) result(n) &
      bind(c, name="usf_mesh_bulk_elements")
    type(c_ptr), value :: h_mesh
    integer(c_int) :: n
    type(mesh_t), pointer :: m
    call c_f_pointer(h_mesh, m)
    n = m % numberofbulkelements
  end function usf_mesh_bulk_elements

  function usf_mesh_max_element_nodes(h_mesh) result(n) &
      bind(c, name="usf_mesh_max_element_nodes")
    type(c_ptr), value :: h_mesh
    integer(c_int) :: n
    type(mesh_t), pointer :: m
    call c_f_pointer(h_mesh, m)
    n = m % maxelementnodes
  end function usf_mesh_max_element_nodes

  function usf_element_constraint(h_elem) result(c) &
      bind(c, name="usf_element_constraint")
    type(c_ptr), value :: h_elem
    integer(c_int) :: c
    type(element_t), pointer :: e
    call c_f_pointer(h_elem, e)
    c = e % boundaryinfo % constraint
  end function usf_element_constraint

  subroutine usf_model_bc_values(bc_i, h_list) &
      bind(c, name="usf_model_bc_values")
    integer(c_int), value :: bc_i
    type(c_ptr) :: h_list
    type(valuelist_t), pointer :: l
    l => acc_model % bcs(bc_i) % values
    if (associated(l)) then
      h_list = c_loc(l)
    else
      h_list = c_null_ptr
    end if
  end subroutine usf_model_bc_values

  function usf_model_bc_tag(bc_i) result(t) bind(c, name="usf_model_bc_tag")
    integer(c_int), value :: bc_i
    integer(c_int) :: t
    t = acc_model % bcs(bc_i) % tag
  end function usf_model_bc_tag

  function usf_number_of_bcs() result(n) bind(c, name="usf_number_of_bcs")
    integer(c_int) :: n
    n = acc_model % numberofboundaries
  end function usf_number_of_bcs

  subroutine usf_make_perm(which, name, namelen, nfront) &
      bind(c, name="usf_make_perm")
    use meshbasics
    integer(c_int), value :: which
    type(c_ptr), value :: name
    integer(c_int), value :: namelen
    integer(c_int) :: nfront
    character(len=:), allocatable :: key
    character(kind=c_char), pointer :: cbuf(:)
    integer :: i
    call c_f_pointer(name, cbuf, [namelen])
    allocate(character(len=namelen) :: key)
    do i = 1, namelen
      key(i:i) = cbuf(i)
    end do
    if (which == 1) then
      if (associated(acc_front_perm)) deallocate(acc_front_perm)
      call MakePermUsingMask(acc_model, acc_model % solver, acc_model % mesh, &
          key, .FALSE., acc_front_perm, nfront, ParallelComm=.FALSE.)
    else if (which == 2) then
      if (associated(acc_bottom_perm)) deallocate(acc_bottom_perm)
      call MakePermUsingMask(acc_model, acc_model % solver, acc_model % mesh, &
          key, .FALSE., acc_bottom_perm, nfront, ParallelComm=.FALSE.)
    else
      if (associated(acc_surface_perm)) deallocate(acc_surface_perm)
      call MakePermUsingMask(acc_model, acc_model % solver, acc_model % mesh, &
          key, .FALSE., acc_surface_perm, nfront, ParallelComm=.FALSE.)
    end if
  end subroutine usf_make_perm

  function usf_perm_at(which, i) result(p) bind(c, name="usf_perm_at")
    integer(c_int), value :: which
    integer(c_int), value :: i
    integer(c_int) :: p
    if (which == 1) then
      p = acc_front_perm(i)
    else if (which == 2) then
      p = acc_bottom_perm(i)
    else
      p = acc_surface_perm(i)
    end if
  end function usf_perm_at

  subroutine usf_point_in_element(h_elem, u, v, w, lc, eg, el, en, &
      found) bind(c, name="usf_point_in_element")
    use interpolation
    type(c_ptr), value :: h_elem
    real(c_double), value :: u
    real(c_double), value :: v
    real(c_double), value :: w
    real(c_double) :: lc(3)
    real(c_double), value :: eg
    real(c_double), value :: el
    real(c_double), value :: en
    integer(c_int) :: found
    type(element_t), pointer :: e
    type(nodes_t), save, target :: acc_nodes_pi
    real(c_double) :: pnt(3)
    real(c_double) :: lcc(3)
    logical :: f
    call c_f_pointer(h_elem, e)
    call GetElementNodes(acc_nodes_pi, UElement=e)
    pnt = [u, v, w]
    f = PointInElement(e, acc_nodes_pi, pnt, lcc, eg, el, en)
    lc = lcc
    found = merge(1, 0, f)
  end subroutine usf_point_in_element

  function usf_interpolate_in_element(h_elem, vals, n, u, v, w) &
      result(r) bind(c, name="usf_interpolate_in_element")
    use elementbasis
    type(c_ptr), value :: h_elem
    type(c_ptr), value :: vals
    integer(c_int), value :: n
    real(c_double), value :: u
    real(c_double), value :: v
    real(c_double), value :: w
    real(c_double) :: r
    type(element_t), pointer :: e
    real(c_double), pointer :: fv(:)
    call c_f_pointer(h_elem, e)
    call c_f_pointer(vals, fv, [n])
    r = InterpolateInElement(e, fv, u, v, w)
  end function usf_interpolate_in_element

  function usf_mesh_boundary_elements(h_mesh) result(n) &
      bind(c, name="usf_mesh_boundary_elements")
    type(c_ptr), value :: h_mesh
    integer(c_int) :: n
    type(mesh_t), pointer :: m
    call c_f_pointer(h_mesh, m)
    n = m % numberofboundaryelements
  end function usf_mesh_boundary_elements

end module usf_acc

  ! Fortran entry shim for `glaciermeshmetricaniso_` (the single 4-arg
  ! vector entry). `TargetLength(:)` is assumed-shape: the caller
  ! (Elmer's `ElmerRealArrFn` interface) passes a descriptor, which must
  ! not cross into Rust. The shim converts it to base address + length
  ! and forwards to the Rust core; entry logic lives in `metric.rs`.
  subroutine GlacierMeshMetricAniso(Model, nodenumber, y, TargetLength)
    use iso_c_binding
    use defutils
    type(model_t), target :: Model
    integer :: nodenumber
    real(c_double) :: y
    real(c_double), target :: TargetLength(:)
    interface
      subroutine usf_metric_core(h_model, node, yy, target, ntarget) &
          bind(c, name="usf_metric_core")
        import c_ptr, c_int, c_double
        type(c_ptr), value :: h_model
        integer(c_int), value :: node
        real(c_double), value :: yy
        type(c_ptr), value :: target
        integer(c_int), value :: ntarget
      end subroutine usf_metric_core
    end interface
    type(c_ptr) :: base
    if (size(TargetLength) >= 1) then
      base = c_loc(TargetLength(1))
    else
      base = c_null_ptr
    end if
    call usf_metric_core(c_loc(Model), nodenumber, y, base, &
        size(TargetLength))
  end subroutine GlacierMeshMetricAniso

