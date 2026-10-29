! T021-B ScannedFieldSolver port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/ScannedFieldSolver.F90, compiled against
! the built .mod files of the target configure (ADR-057 §3, T020 Poisson
! pattern).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape integers). All Fortran-object
! traffic — TYPE(...) handles, OPTIONAL/POINTER dummies, assumed-shape
! descriptors, hidden lengths, deferred-length string results — stays inside
! Fortran, so the Rust side is pure C calls. The wrapper module mirrors the
! original's SAVE state (First flag, FieldVars pointer array, FieldNames,
! NofFields, Mesh) with identical lifetimes, and transcribes the order-
! sensitive whole-array updates, i2s-composed names, Info/Fatal text, and
! the SumFieldVar association-check shape statement for statement.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module scannedfieldsolver_acc
  use iso_c_binding
  use types
  use defutils
  implicit none

  type :: sfa_varptr_t
    type(variable_t), pointer :: var => null()
  end type sfa_varptr_t

  type(model_t), pointer, save :: sfa_model => null()
  type(solver_t), pointer, save :: sfa_solver => null()
  ! SAVE mirror of ScannedFieldSolver's internal state (First=.TRUE. plus
  ! the SAVEd FieldVars, FieldNames, NofFields, Mesh).
  logical, save :: sfa_first = .true.
  type(sfa_varptr_t), pointer, save :: sfa_fieldvars(:) => null()
  character(len=MAX_NAME_LEN), allocatable, save :: sfa_fieldnames(:)
  integer, save :: sfa_nofields = 0
  type(mesh_t), pointer, save :: sfa_mesh => null()
  ! Per-field current selections within one entry invocation.
  type(variable_t), pointer, save :: sfa_fieldvar => null()
  type(variable_t), pointer, save :: sfa_scanfieldvar => null()
  type(variable_t), pointer, save :: sfa_sumfieldvar => null()

contains

  ! -- init entry -------------------------------------------------------
  subroutine sfa_init_begin(h_model, h_solver) bind(c, name="sfa_init_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, sfa_model)
    call c_f_pointer(h_solver, sfa_solver)
  end subroutine sfa_init_begin

  ! DO WHILE loop counting pre-existing "Exported Variable " entries;
  ! returns the first free index, exactly as the original's `i`.
  function sfa_init_count_exported() result(r) bind(c, name="sfa_init_count_exported")
    integer(c_int) :: r
    type(valuelist_t), pointer :: params
    integer :: k
    params => GetSolverParams()
    k = 1
    do while (.true.)
      if (.not. ListCheckPresent(params, &
          "Exported Variable "//i2s(k)) ) exit
      k = k + 1
    end do
    r = k
  end function sfa_init_count_exported

  ! DO WHILE loop counting "Field Name " entries; returns NofFields - 1.
  function sfa_init_count_fields() result(r) bind(c, name="sfa_init_count_fields")
    integer(c_int) :: r
    type(valuelist_t), pointer :: params
    integer :: n
    params => GetSolverParams()
    n = 1
    do while (.true.)
      if (.not. ListCheckPresent(params, &
          "Field Name "//i2s(n)) ) exit
      n = n + 1
    end do
    n = n - 1
    r = n
  end function sfa_init_count_fields

  subroutine sfa_init_add_variable() bind(c, name="sfa_init_add_variable")
    type(valuelist_t), pointer :: params
    params => GetSolverParams()
    call ListAddString( params, 'Variable', '-nooutput SF_dummy' )
  end subroutine sfa_init_add_variable

  subroutine sfa_init_add_nfields(n) bind(c, name="sfa_init_add_nfields")
    integer(c_int), value :: n
    type(valuelist_t), pointer :: params
    params => GetSolverParams()
    call ListAddInteger( params, 'Number of Fields', n)
  end subroutine sfa_init_add_nfields

  function sfa_init_scan_solver() result(r) bind(c, name="sfa_init_scan_solver")
    integer(c_int) :: r
    type(valuelist_t), pointer :: params
    logical :: found
    params => GetSolverParams()
    r = GetInteger(params,'Scan Solver',Found)
    if (.not. found) call Fatal('ScannedFieldSolver_Init','Scan Solver not found.')
  end function sfa_init_scan_solver

  function sfa_init_scanmax(solverint) result(r) bind(c, name="sfa_init_scanmax")
    integer(c_int), value :: solverint
    integer(c_int) :: r
    type(valuelist_t), pointer :: params
    logical :: found
    params => GetSolverParams()
    r = GetInteger(sfa_model % solvers(solverint) % values, 'Scanning Loops', Found)
    if (.not. found) call Fatal('ScannedFieldSolver_Init','Scan Loops not found.')
    call ListAddInteger(params, 'Scan Max', r)
  end function sfa_init_scanmax

  ! One FieldInt trip of the export loop: reads 'Field Name '//i2s,
  ! appends ScanMax 'Scan i name' exports plus the '-nooutput name Dummy'
  ! export, advancing `i` exactly as the original.
  subroutine sfa_init_field_exports(fieldint, scanmax, i) bind(c, name="sfa_init_field_exports")
    integer(c_int), value :: fieldint, scanmax
    integer(c_int) :: i
    type(valuelist_t), pointer :: params
    character(len=MAX_NAME_LEN) :: fieldname
    logical :: found
    integer :: scanint
    params => GetSolverParams()
    fieldname = GetString(params,'Field Name '//i2s(fieldint),Found)
    if (.not. found) call Fatal('ScannedFieldSolver_Init','Field Name '//i2s(fieldint)//' not found.')

    do scanint = 1, scanmax
      call ListAddString( params, "Exported Variable "//i2s(i), &
        'Scan '//i2s(scanint)//" "//trim(fieldname))
      i=i+1
    end do

    call ListAddString( params, "Exported Variable "//i2s(i), &
      '-nooutput '//trim(fieldname)//' Dummy')
    i=i+1
  end subroutine sfa_init_field_exports

  ! -- main entry -------------------------------------------------------
  subroutine sfa_begin(h_model, h_solver) bind(c, name="sfa_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, sfa_model)
    call c_f_pointer(h_solver, sfa_solver)
  end subroutine sfa_begin

  ! The IF(First) block, guarded by the mirrored SAVE flag with the same
  ! lifetime as the original's internal-procedure SAVE.
  subroutine sfa_first_block() bind(c, name="sfa_first_block")
    type(valuelist_t), pointer :: sp
    integer :: fieldint, istat, n
    if (sfa_first) then
      sfa_first = .false.

      sfa_mesh => sfa_model % mesh
      n = sfa_mesh % maxelementdofs

      sfa_nofields = ListGetInteger(GetSolverParams(), 'Number of Fields')

      allocate(sfa_fieldnames(sfa_nofields), sfa_fieldvars(sfa_nofields), stat=istat)
      if ( istat /= 0 ) call Fatal('ScannedFieldSolver','Memory allocation error')

      sp => GetSolverParams()
      do fieldint = 1, sfa_nofields
        sfa_fieldnames(fieldint) = ListGetString(sp, 'Field Name '//i2s(fieldint))
        sfa_fieldvars(fieldint) % var => VariableGet( sfa_mesh % variables, sfa_fieldnames(fieldint))
        if (.not. associated(sfa_fieldvars(fieldint) % var)) then
          call Fatal('ScannedFieldSolver',trim(sfa_fieldnames(fieldint))//' not associated!')
        end if
      end do
    end if
  end subroutine sfa_first_block

  subroutine sfa_info_header() bind(c, name="sfa_info_header")
    call Info('ScannedFieldSolver','-------------------------------------------',Level=6)
    call Info('ScannedFieldSolver','Computing the scan field solutions. ',Level=5)
    call Info('ScannedFieldSolver','-------------------------------------------',Level=6)
  end subroutine sfa_info_header

  function sfa_scan_index() result(r) bind(c, name="sfa_scan_index")
    integer(c_int) :: r
    type(variable_t), pointer :: scanvar
    scanvar => VariableGet(sfa_solver % mesh % variables, 'scan')
    if (.not. associated(scanvar)) call Fatal('ScannedFieldSolver', 'Scan variable not found.')
    r = int( scanvar % values(1) )
  end function sfa_scan_index

  function sfa_scanmax() result(r) bind(c, name="sfa_scanmax")
    integer(c_int) :: r
    type(valuelist_t), pointer :: params
    logical :: found
    params => GetSolverParams()
    r = GetInteger(params, 'Scan Max', Found)
    if (.not. found) call Fatal('ScannedFieldSolver', 'The maximum number of scan loops not found.')
  end function sfa_scanmax

  function sfa_nfields() result(r) bind(c, name="sfa_nfields")
    integer(c_int) :: r
    r = sfa_nofields
  end function sfa_nfields

  subroutine sfa_field_select(fi) bind(c, name="sfa_field_select")
    integer(c_int), value :: fi
    sfa_fieldvar => sfa_fieldvars(fi) % var
    call Info('ScannedFieldSolver','Field '//trim(sfa_fieldnames(fi))//'.',Level=5)
  end subroutine sfa_field_select

  subroutine sfa_field_scan(fi, scanint) bind(c, name="sfa_field_scan")
    integer(c_int), value :: fi, scanint
    character(len=MAX_NAME_LEN) :: scanfieldname
    scanfieldname = 'scan '//i2s(scanint)//" "//trim(sfa_fieldnames(fi))
    sfa_scanfieldvar => VariableGet( sfa_mesh % variables, scanfieldname)
    if (.not. associated(sfa_scanfieldvar)) then
      call Fatal('ScannedFieldSolver', trim(scanfieldname)//' not associated!')
    end if

    call sfa_make_similar(sfa_scanfieldvar, sfa_fieldvar)
    if (size(sfa_scanfieldvar % values) /= size(sfa_fieldvar % values)) &
      call Fatal('ScannedFieldSolver','Scanned fields are of different size than &
      the defined Scan Field Variable.')
  end subroutine sfa_field_scan

  subroutine sfa_field_sum(fi) bind(c, name="sfa_field_sum")
    integer(c_int), value :: fi
    character(len=MAX_NAME_LEN) :: sumfieldname
    sumfieldname = trim(sfa_fieldnames(fi))//' Dummy'
    sfa_sumfieldvar => VariableGet( sfa_mesh % variables, sumfieldname)
    ! NOTE: the original checks ScanFieldVar here rather than SumFieldVar;
    ! replicated exactly (bug-compatible).
    if (.not. associated(sfa_scanfieldvar)) then
      call Fatal('ScannedFieldSolver', trim(sumfieldname)//' not associated')
    end if

    call sfa_make_similar(sfa_sumfieldvar, sfa_fieldvar)
    if (size(sfa_scanfieldvar % values) /= size(sfa_sumfieldvar % values)) &
      call Fatal('ScannedFieldSolver','Summed fields are of different size than &
      the defined Sum Field Variable.')
  end subroutine sfa_field_sum

  ! Whole-array updates in the original evaluation order.
  subroutine sfa_field_combine(scanint, scanmax) bind(c, name="sfa_field_combine")
    integer(c_int), value :: scanint, scanmax
    sfa_scanfieldvar % values = sfa_fieldvar % values
    if ( scanint == 1 ) then
      sfa_sumfieldvar % values = sfa_fieldvar % values
    else
      sfa_sumfieldvar % values = sfa_fieldvar % values + sfa_sumfieldvar % values
    end if

    if ( scanint == scanmax ) then
      sfa_fieldvar % values = sfa_sumfieldvar % values
    end if
  end subroutine sfa_field_combine

  ! The MakeVarSimilarToModelVar internal procedure, transcribed exactly.
  subroutine sfa_make_similar(var, modelvar)
    type(variable_t), target :: var, modelvar
    integer :: istat

    if (size(var % values) /= size(modelvar % values)) then
      deallocate(var % values)
      allocate(var % values(size(modelvar % values)), stat=istat)
      if ( istat /= 0 ) call Fatal('MakeVarSimilarToModelVar','Memory allocation error')
    end if

    if (.not. associated(var % perm,  modelvar % perm)) var % perm => modelvar % perm
    if ( var % dofs /= modelvar % dofs) var % dofs = modelvar % dofs
    if ( var % type /= modelvar % type ) var % type = modelvar % type
  end subroutine sfa_make_similar

end module scannedfieldsolver_acc
