! SaveData pilot: statement-level BIND(C) wrappers for the Elmer calls in
! fem/src/modules/SaveData/SaveBoundaryValues.F90.
! Same contract as the ConsistentSplitting accessors (C-only arguments;
! SAVE state mirrored in Fortran; the entry's own SAVE scalars live in
! the Rust state, element traversal stays in Fortran wrappers).
module sd_bvals_acc
  use iso_c_binding
  use defutils
  use solverbasics
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(model_t), pointer, save :: acc_model => null()
  character(len=512), save :: acc_message

contains

  subroutine sdb_begin(h_model, h_solver) bind(c, name="sdb_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine sdb_begin

  function sdb_dim() result(d) bind(c, name="sdb_dim")
    integer(c_int) :: d
    d = CoordinateSystemDimension()
  end function sdb_dim

  function sdb_local_nodes() result(n) bind(c, name="sdb_local_nodes")
    integer(c_int) :: n
    n = acc_model % numberofnodes
  end function sdb_local_nodes

  function sdb_max_elem_nodes() result(n) bind(c, name="sdb_max_elem_nodes")
    integer(c_int) :: n
    n = acc_solver % mesh % maxelementnodes
  end function sdb_max_elem_nodes

  subroutine sdb_info_creating() bind(c, name="sdb_info_creating")
    call Info('SaveBoundaryValues', 'Creating selected boundary values as fields')
  end subroutine sdb_info_creating

  subroutine sdb_param(i, buf, cap, found) bind(c, name="sdb_param")
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
  end subroutine sdb_param

  subroutine sdb_warn_no_params() bind(c, name="sdb_warn_no_params")
    call Warn('SaveBoundaryValues', 'No parameters found: No fields will be created')
  end subroutine sdb_warn_no_params

  subroutine sdb_create_var(cname, local_nodes, exists, created, valid, ndofs) &
      bind(c, name="sdb_create_var")
    character(kind=c_char) :: cname(*)
    integer(c_int), value :: local_nodes
    integer(c_int) :: exists, created, valid, ndofs
    character(len=4096) :: pname
    type(variable_t), pointer :: var
    type(valuelist_t), pointer :: vlist
    integer, pointer :: fperm(:)
    real(kind=dp), pointer :: field(:)
    type(element_t), pointer :: el
    integer, pointer :: nidx(:)
    integer :: i, j, n, elno, istat
    call sdb_to_fortran(cname, pname)
    var => VariableGet(acc_model % variables, trim(pname), .true.)
    exists = merge(1, 0, associated(var))
    created = 0
    valid = 1
    ndofs = 0
    if (exists /= 0) return
    allocate(fperm(local_nodes), stat=istat)
    if (istat /= 0) call Fatal('SaveBoundaryValues', 'Memory allocation error 2')
    fperm = 0
    do elno = acc_model % numberofbulkelements + 1, &
        acc_model % numberofbulkelements + acc_model % numberofboundaryelements
      el => acc_model % elements(elno)
      if (el % partindex /= parenv % mype) cycle
      n = GetElementNOFNodes(el)
      nidx => el % nodeindexes
      acc_model % currentelement => el
      vlist => GetBC(el)
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
      call Warn('SaveBoundaryValues', &
          'Parameter ' // trim(pname) // ' not present in any material')
      valid = 0
      deallocate(fperm)
      return
    else
      write (acc_message, '(A,I0,A)') 'Parameter > ' // trim(pname) // &
          ' < defined with ', j, ' dofs'
      call Info('SaveBoundaryValues', acc_message)
      allocate(field(j), stat=istat)
      if (istat /= 0) call Fatal('SaveBoundaryValues', 'Memory allocation error 3')
      field = 0.0_dp
      call VariableAdd(acc_model % variables, acc_solver % mesh, acc_solver, &
          trim(pname), 1, field, fperm)
      nullify(field)
      nullify(fperm)
      created = 1
    end if
  end subroutine sdb_create_var

  subroutine sdb_update_var(cname, localparam, maxn) bind(c, name="sdb_update_var")
    character(kind=c_char) :: cname(*)
    real(c_double) :: localparam(*)
    integer(c_int), value :: maxn
    character(len=4096) :: pname
    type(variable_t), pointer :: var
    real(kind=dp), pointer :: field(:)
    integer, pointer :: fperm(:)
    type(element_t), pointer :: el
    integer, pointer :: nidx(:)
    type(valuelist_t), pointer :: vlist
    integer :: elno, n
    logical :: gotit
    call sdb_to_fortran(cname, pname)
    var => VariableGet(acc_model % variables, trim(pname), .true.)
    if (.not. associated(var)) return
    field => var % values
    fperm => var % perm
    do elno = acc_model % numberofbulkelements + 1, &
        acc_model % numberofbulkelements + acc_solver % mesh % numberofboundaryelements
      el => acc_model % elements(elno)
      if (el % partindex /= parenv % mype) cycle
      n = GetElementNOFNodes(el)
      nidx => el % nodeindexes
      if (associated(fperm)) then
        if (.not. all(fperm(nidx) > 0)) cycle
      end if
      acc_model % currentelement => el
      vlist => GetBC(el)
      call sdb_getreal_list(vlist, pname, n, nidx, localparam, maxn, gotit)
      if (.not. gotit) cycle
      if (associated(fperm)) then
        field(fperm(nidx(1:n))) = localparam(1:n)
      else
        field(nidx(1:n)) = localparam(1:n)
      end if
    end do
  end subroutine sdb_update_var

  subroutine sdb_getreal_list(vlist, pname, n, nidx, vals, maxn, gotit)
    type(valuelist_t), pointer :: vlist
    character(len=*) :: pname
    integer :: n, maxn
    integer :: nidx(n)
    real(c_double) :: vals(*)
    logical :: gotit
    real(kind=dp) :: tmp(maxn)
    tmp(1:n) = ListGetReal(vlist, trim(pname), n, nidx, gotit)
    if (gotit) vals(1:n) = tmp(1:n)
  end subroutine sdb_getreal_list

  subroutine sdb_fatal_alloc1() bind(c, name="sdb_fatal_alloc1")
    call Fatal('SaveBoundaryValues', 'Memory allocation error 1')
  end subroutine sdb_fatal_alloc1

  subroutine sdb_to_fortran(cstr, fstr)
    character(kind=c_char) :: cstr(*)
    character(len=*) :: fstr
    integer :: i
    fstr = ' '
    do i = 1, len(fstr)
      if (cstr(i) == c_null_char) exit
      fstr(i:i) = cstr(i)
    end do
  end subroutine sdb_to_fortran

end module sd_bvals_acc
