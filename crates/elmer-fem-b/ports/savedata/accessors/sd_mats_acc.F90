! SaveData pilot: statement-level BIND(C) wrappers for the Elmer calls in
! fem/src/modules/SaveData/SaveMaterials.F90.
! Mirrors the BoundaryValues wrapper shape, plus DG-field support,
! BodyForce/Material selection, and the PrevNodes mesh-change guard.
module sd_mats_acc
  use iso_c_binding
  use defutils
  use solverbasics
  use types
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(model_t), pointer, save :: acc_model => null()
  character(len=512), save :: acc_message
  integer, save :: acc_prev_nodes = 0

contains

  subroutine sdm_begin(h_model, h_solver) bind(c, name="sdm_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine sdm_begin

  function sdm_dim() result(d) bind(c, name="sdm_dim")
    integer(c_int) :: d
    d = CoordinateSystemDimension()
  end function sdm_dim

  function sdm_local_nodes() result(n) bind(c, name="sdm_local_nodes")
    integer(c_int) :: n
    n = acc_model % numberofnodes
  end function sdm_local_nodes

  function sdm_max_elem_nodes() result(n) bind(c, name="sdm_max_elem_nodes")
    integer(c_int) :: n
    n = acc_solver % mesh % maxelementnodes
  end function sdm_max_elem_nodes

  subroutine sdm_info_creating() bind(c, name="sdm_info_creating")
    call Info('SaveMaterials', 'Creating selected material parameters as fields')
  end subroutine sdm_info_creating

  function sdm_check_mesh(local_nodes) result(changed) bind(c, name="sdm_check_mesh")
    integer(c_int), value :: local_nodes
    integer(c_int) :: changed
    changed = 0
    if (acc_prev_nodes > 0) then
      if (local_nodes /= acc_prev_nodes) changed = 1
    end if
    acc_prev_nodes = local_nodes
  end function sdm_check_mesh

  subroutine sdm_warn_remesh() bind(c, name="sdm_warn_remesh")
    call Warn('SaveMaterials', 'Solver cannot deal with changed mesh: Exiting!')
  end subroutine sdm_warn_remesh

  subroutine sdm_param(i, buf, cap, found) bind(c, name="sdm_param")
    integer(c_int), value :: i, cap
    character(kind=c_char) :: buf(*)
    integer(c_int) :: found
    type(valuelist_t), pointer :: sp
    character(len=4096) :: name, val
    logical :: got
    integer :: k, n
    sp => GetSolverParams()
    write (name, '(A,I0)') 'Parameter ', i
    val = ListGetString(sp, trim(name), got)
    found = merge(1, 0, got)
    if (got) then
      n = min(len_trim(val), cap - 1)
      do k = 1, n
        buf(k) = val(k:k)
      end do
      buf(n + 1) = c_null_char
    end if
  end subroutine sdm_param

  subroutine sdm_warn_no_params() bind(c, name="sdm_warn_no_params")
    call Warn('SaveMaterials', 'No parameters found: No fields will be created')
  end subroutine sdm_warn_no_params

  function sdm_bodyforce_params() result(n) bind(c, name="sdm_bodyforce_params")
    integer(c_int) :: n
    type(valuelist_t), pointer :: sp
    logical :: got
    sp => GetSolverParams()
    n = ListGetInteger(sp, 'Body Force Parameters', got)
  end function sdm_bodyforce_params

  subroutine sdm_create_var(cname, is_bf, local_nodes, exists, created, ndofs) &
      bind(c, name="sdm_create_var")
    character(kind=c_char) :: cname(*)
    integer(c_int), value :: is_bf, local_nodes
    integer(c_int) :: exists, created, ndofs
    character(len=4096) :: pname
    type(variable_t), pointer :: var
    type(valuelist_t), pointer :: vlist
    integer, pointer :: fperm(:)
    real(kind=dp), pointer :: field(:)
    type(element_t), pointer :: el
    integer, pointer :: nidx(:)
    integer :: i, j, n, elno, istat
    call sdm_to_fortran(cname, pname)
    var => VariableGet(acc_model % variables, trim(pname), .true.)
    exists = merge(1, 0, associated(var))
    created = 0
    ndofs = 0
    if (exists /= 0) return
    allocate(fperm(local_nodes), stat=istat)
    if (istat /= 0) call Fatal('SaveMaterials', 'Memory allocation error 2')
    fperm = 0
    do elno = 1, acc_model % numberofbulkelements + acc_model % numberofboundaryelements
      el => acc_model % elements(elno)
      n = GetElementNOFNodes(el)
      nidx => el % nodeindexes
      acc_model % currentelement => el
      if (is_bf == 0) then
        vlist => GetMaterial(el)
      else
        vlist => GetBodyForce(el)
      end if
      if (.not. associated(vlist)) cycle
      if (ListCheckPresent(vlist, trim(pname))) fperm(nidx) = 1
    end do
    j = 0
    do i = 1, local_nodes
      if (fperm(i) > 0) then
        j = j + 1
        fperm(i) = j
      end if
    end do
    ndofs = j
    if (j == 0) then
      call Warn('SaveMaterials', &
          'Parameter ' // trim(pname) // ' not present in any material')
    end if
    write (acc_message, '(A,I0,A)') 'Parameter > ' // trim(pname) // &
        ' < defined with ', j, ' dofs'
    call Info('SaveMaterials', acc_message)
    if (j > 0) then
      allocate(field(j), stat=istat)
      if (istat /= 0) call Fatal('SaveMaterials', 'Memory allocation error 3')
      field = 0.0_dp
    end if
    if (j > 0 .or. parenv % pes > 1) then
      call VariableAdd(acc_model % variables, acc_solver % mesh, acc_solver, &
          trim(pname), 1, field, fperm)
      created = 1
    end if
    if (associated(field)) nullify(field)
    nullify(fperm)
  end subroutine sdm_create_var

  subroutine sdm_update_var(cname, is_bf, localparam, maxn) bind(c, name="sdm_update_var")
    character(kind=c_char) :: cname(*)
    integer(c_int), value :: is_bf
    real(c_double) :: localparam(*)
    integer(c_int), value :: maxn
    character(len=4096) :: pname
    type(variable_t), pointer :: var
    real(kind=dp), pointer :: field(:)
    integer, pointer :: fperm(:)
    type(element_t), pointer :: el
    integer, pointer :: nidx(:)
    type(valuelist_t), pointer :: vlist
    integer :: elno, n, indexes(27)
    logical :: gotit
    call sdm_to_fortran(cname, pname)
    var => VariableGet(acc_model % variables, trim(pname), .true.)
    if (.not. associated(var)) return
    field => var % values
    if (.not. associated(field)) return
    fperm => var % perm
    if (.not. associated(fperm)) return
    if (maxval(fperm) <= 0) return
    do elno = 1, acc_model % numberofbulkelements + acc_model % numberofboundaryelements
      el => acc_model % elements(elno)
      if (el % partindex /= parenv % mype) cycle
      n = GetElementNOFNodes(el)
      nidx => el % nodeindexes
      if (var % type == variable_on_nodes) then
        indexes(1:n) = nidx(1:n)
      else if (var % type == variable_on_nodes_on_elements) then
        if (.not. associated(el % dgindexes)) cycle
        indexes(1:n) = el % dgindexes(1:n)
      else
        call Fatal('SaveMaterials', 'Only implemented for nodal and dg fields: ' // trim(pname))
      end if
      if (.not. all(fperm(indexes(1:n)) > 0)) cycle
      acc_model % currentelement => el
      if (is_bf == 0) then
        vlist => GetMaterial(el)
      else
        vlist => GetBodyForce(el)
      end if
      if (.not. associated(vlist)) cycle
      call sdm_getreal_list(vlist, pname, n, nidx, localparam, maxn, gotit)
      if (.not. gotit) cycle
      field(fperm(indexes(1:n))) = localparam(1:n)
    end do
  end subroutine sdm_update_var

  subroutine sdm_getreal_list(vlist, pname, n, nidx, vals, maxn, gotit)
    type(valuelist_t), pointer :: vlist
    character(len=*) :: pname
    integer :: n, maxn
    integer :: nidx(n)
    real(c_double) :: vals(*)
    logical :: gotit
    real(kind=dp) :: tmp(maxn)
    tmp(1:n) = ListGetReal(vlist, trim(pname), n, nidx, gotit)
    if (gotit) vals(1:n) = tmp(1:n)
  end subroutine sdm_getreal_list

  subroutine sdm_fatal_alloc1() bind(c, name="sdm_fatal_alloc1")
    call Fatal('SaveMaterials', 'Memory allocation error 1')
  end subroutine sdm_fatal_alloc1

  subroutine sdm_to_fortran(cstr, fstr)
    character(kind=c_char) :: cstr(*)
    character(len=*) :: fstr
    integer :: i
    fstr = ' '
    do i = 1, len(fstr)
      if (cstr(i) == c_null_char) exit
      fstr(i:i) = cstr(i)
    end do
  end subroutine sdm_to_fortran

end module sd_mats_acc
