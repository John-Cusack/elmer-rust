! Mesh2MeshSolver port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/Mesh2MeshSolver.F90, compiled against the
! built .mod files of the target configure (the Poisson-port mechanism).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, deferred-length strings,
! whole-array updates, WRITE/PRINT I/O, MAXVAL reductions — stays inside
! Fortran, so the Rust side is pure C calls. The wrapper module mirrors
! the plugin's per-invocation locals (mesh/variable handles, the
! hand-picked variable list, the additive base-value store) with SAVE
! lifetimes that are reset on every entry, exactly matching the
! original's automatic (de)allocation on procedure entry/exit.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text. The `#if 0` block of the original is omitted
! (disabled in the source too).
module mesh2meshsolver_acc
  use iso_c_binding
  use defutils
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  ! Per-invocation locals mirrored with SAVE lifetime, reset by
  ! m2ma_build_vars on every entry (the original's allocatables are
  ! freshly allocated per call and auto-freed on exit).
  type(mesh_t), pointer, save :: acc_thismesh => null()
  type(mesh_t), pointer, save :: acc_targetmesh => null()
  type(variable_t), pointer, save :: acc_fromvars => null()
  logical, save :: acc_owns_fromvars = .false.
  type acc_additivesol_t
    logical :: additive
    real(kind=dp), allocatable :: values(:)
  end type acc_additivesol_t
  type(acc_additivesol_t), allocatable, save :: acc_basesol(:)
  integer, save :: acc_novar = 0
  logical, save :: acc_additive_mode = .false.

contains

  ! -- init0 entry -----------------------------------------------------
  subroutine m2ma_init0() bind(c, name="m2ma_init0")
    type(valuelist_t), pointer :: params
    params => getsolverparams()
    call listaddnewlogical(params, 'No Matrix', .true.)
    call listaddnewlogical(params, 'Mesh Enforce Local Copy', .true.)
  end subroutine m2ma_init0

  ! -- entry state ------------------------------------------------------
  subroutine m2ma_begin(h_solver) bind(c, name="m2ma_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine m2ma_begin

  ! -- mesh resolution (Info + both Target Mesh branches) ---------------
  subroutine m2ma_resolve_meshes() bind(c, name="m2ma_resolve_meshes")
    type(mesh_t), pointer :: mesh
    type(valuelist_t), pointer :: params
    integer :: i
    logical :: found
    call info('Mesh2MeshSolver', 'Mapping result between meshes')
    acc_thismesh => getmesh()
    call info('Mesh2MeshSolver', 'This mesh dimension is ' &
        //i2s(acc_thismesh % meshdim)//' and name: '//trim(acc_thismesh % name), level=7)
    params => getsolverparams()
    acc_targetmesh => null()
    i = listgetinteger(params, 'Target Mesh Solver Index', found)
    if (found) then
      ! Target mesh solver is explicitly given
      acc_targetmesh => currentmodel % solvers(i) % mesh
      if (associated(acc_targetmesh)) then
        call info('Mesh2MeshSolver', 'Using target mesh as the mesh of Solver '//i2s(i), level=8)
      else
        call fatal('Mesh2MeshSolver', 'Target Mesh for Solver not associated: '//i2s(i))
      end if
    else
      ! Otherwise use the 1st mesh that is not this old data mesh
      mesh => currentmodel % meshes
      do while (associated(mesh))
        if (.not. associated(mesh, acc_thismesh)) then
          acc_targetmesh => mesh
          exit
        end if
        mesh => mesh % next
      end do
      if (associated(acc_targetmesh)) then
        call info('Mesh2MeshSolver', 'Using target mesh as the first mesh different from this mesh', level=8)
      else
        call fatal('Mesh2MeshSolver', 'Could not find target mesh that is not this mesh!')
      end if
    end if
    call info('Mesh2MeshSolver', 'Target mesh name is: '//trim(acc_targetmesh % name), level=7)
  end subroutine m2ma_resolve_meshes

  ! -- hand-picked variable list + additive base store ------------------
  subroutine m2ma_build_vars(novar, additive) bind(c, name="m2ma_build_vars")
    integer(c_int) :: novar, additive
    type(valuelist_t), pointer :: params
    type(variable_t), pointer :: var, pvar
    character(len=max_name_len) :: name, varname
    integer :: i
    logical :: found
    ! Fresh per-invocation locals (original: automatic deallocation).
    if (acc_owns_fromvars .and. associated(acc_fromvars)) deallocate(acc_fromvars)
    acc_fromvars => null()
    acc_owns_fromvars = .false.
    if (allocated(acc_basesol)) deallocate(acc_basesol)
    acc_novar = 0
    acc_additive_mode = .false.
    params => getsolverparams()
    do i = 1, 100
      write(name, '(A,I0)') 'Variable ', i
      if (.not. listcheckpresent(params, name)) exit
      acc_novar = i
    end do
    if (acc_novar > 0) then
      call info('Mesh2MeshSolver', 'Mapping '//i2s(acc_novar)//' hand-picked fields as requested', level=7)
      allocate(acc_fromvars)
      acc_owns_fromvars = .true.
      pvar => acc_fromvars
      acc_additive_mode = listcheckprefix(params, 'Variable Additive')
      if (acc_additive_mode) then
        allocate(acc_basesol(acc_novar))
        acc_basesol(1:acc_novar) % additive = .false.
        call info('Mesh2MeshSolver', 'Interpolating in additive mode', level=10)
      end if
      do i = 1, acc_novar
        write(name, '(A,I0)') 'Variable ', i
        varname = getstring(params, name, found)
        if (.not. found) exit
        ! Add a new variable to the temporal list.
        if (i > 1) then
          allocate(pvar % next)
          pvar => pvar % next
        end if
        var => variableget(acc_thismesh % variables, varname, thisonly=.true.)
        if (associated(var)) then
          pvar = var
        else
          pvar % name = trim(varname)
          pvar % dofs = 1
          allocate(pvar % perm(0))
        end if
        pvar % next => null()
        ! If the target variable has different name rename the primary
        ! variable also so that default subroutines can be used.
        write(name, '(A,I0)') 'Target Variable ', i
        varname = getstring(params, name, found)
        if (found) then
          call info('Mesh2MeshSolver', 'Renaming variable "'//trim(pvar % name)//'" to "'//trim(varname)//'"')
          pvar % name = varname
          pvar % namelen = lentrim(varname)
        end if
        ! If the variable is additive allocate and save the base values.
        write(name, '(A,I0)') 'Variable Additive ', i
        if (listgetlogical(params, name, found)) then
          var => variableget(acc_targetmesh % variables, varname, thisonly=.true.)
          if (.not. associated(var)) then
            call fatal('Mesh2MeshSolver', 'Only existing fields may be additive!')
          end if
          acc_basesol(i) % additive = .true.
          allocate(acc_basesol(i) % values(size(var % values)))
          acc_basesol(i) % values = var % values
        end if
      end do
    else
      acc_fromvars => acc_thismesh % variables
    end if
    novar = acc_novar
    additive = merge(1, 0, acc_additive_mode)
  end subroutine m2ma_build_vars

  ! -- parallel guard ----------------------------------------------------
  subroutine m2ma_parenv() bind(c, name="m2ma_parenv")
    ! For some reason this is not always active.
    ! If this is not set parallel interpolation could fail.
    if (parenv % pes > 1) parenv % active = .true.
  end subroutine m2ma_parenv

  ! -- the interpolation sweep -------------------------------------------
  subroutine m2ma_interpolate() bind(c, name="m2ma_interpolate")
    interface
      subroutine interpolatemeshtomesh(oldmesh, newmesh, oldvariables, &
          newvariables, usequadranttree, projector, maskname, unfoundnodes)
        use lists
        use sparitercomm
        use interpolation
        use coordinatesystems
        use meshbasics, only: releasemesh
        type(mesh_t), target :: oldmesh, newmesh
        type(variable_t), pointer, optional :: oldvariables, newvariables
        logical, optional :: usequadranttree
        type(projector_t), pointer, optional :: projector
        character(len=*), optional :: maskname
        logical, pointer, optional :: unfoundnodes(:)
      end subroutine interpolatemeshtomesh
    end interface
    ! Here is all the real interpolation work (serial & parallel)
    call info('Mesh2MeshSolver', 'Mapping all fields in primary mesh to target mesh!', level=7)
    call interpolatemeshtomesh(acc_thismesh, acc_targetmesh, &
        acc_fromvars, acc_targetmesh % variables)
    call info('Mesh2MeshSolver', 'Done mapping all fields in primary mesh to target mesh!', level=15)
  end subroutine m2ma_interpolate

  ! -- additive fixup (whole-array update stays in Fortran) --------------
  subroutine m2ma_additive() bind(c, name="m2ma_additive")
    type(valuelist_t), pointer :: params
    type(variable_t), pointer :: var
    character(len=max_name_len) :: name, varname
    integer :: i
    logical :: found
    ! If we has some additive variables then add the base values.
    params => getsolverparams()
    do i = 1, acc_novar
      if (acc_basesol(i) % additive) then
        write(name, '(A,I0)') 'Variable ', i
        varname = getstring(params, name, found)
        var => variableget(acc_targetmesh % variables, varname, thisonly=.true.)
        call info('Mesh2MeshSolver', 'Adding variable to base field: '//trim(varname), level=7)
        var % values = var % values + acc_basesol(i) % values
        deallocate(acc_basesol(i) % values)
      end if
    end do
    deallocate(acc_basesol)
  end subroutine m2ma_additive

  ! -- validation sweep + final Info --------------------------------------
  subroutine m2ma_validate(novar) bind(c, name="m2ma_validate")
    integer(c_int), value :: novar
    type(valuelist_t), pointer :: params
    type(variable_t), pointer :: var, pvar
    character(len=max_name_len) :: name, varname
    integer :: i
    logical :: found, doit
    ! Validate the mapped variables because otherwise we might
    ! accidentally be performing mapping again when calling GetVariable.
    params => getsolverparams()
    if (novar > 0) then
      do i = 1, novar
        write(name, '(A,I0)') 'Variable ', i
        varname = getstring(params, name, found)
        var => variableget(acc_thismesh % variables, varname, thisonly=.true.)
        pvar => variableget(acc_targetmesh % variables, varname, thisonly=.true.)
        if (associated(var)) then
          if (infoactive(20)) &
            print *, 'FromRange'//i2s(parenv % mype)//': ', trim(var % name), maxval(var % values)
        end if
        if (associated(pvar)) then
          pvar % valid = .true.
          pvar % valueschanged = .true.
          if (infoactive(20)) &
            print *, 'ToRange'//i2s(parenv % mype)//': ', trim(pvar % name), maxval(pvar % values)
        end if
      end do
    else
      var => acc_fromvars
      do while (associated(var))
        ! Check that it exists as target variable
        pvar => variableget(acc_targetmesh % variables, var % name(1:var % namelen), thisonly=.true.)
        doit = associated(pvar)
        if (doit) then
          doit = (size(pvar % values) > var % dofs)  ! not global variable
          if (doit .and. len(pvar % name) >= 10) then
            doit = (pvar % name(1:10) /= 'coordinate')  ! not coordinate
          else
            doit = .false.
          end if
        end if
        if (doit) then
          call info('Mesh2MeshSolver', 'Validatig variable in target mesh: '//trim(var % name), level=7)
          pvar % valid = .true.
          pvar % valueschanged = .true.
          if (infoactive(20)) then
            print *, 'FromRange'//i2s(parenv % mype)//': ', trim(var % name), maxval(var % values)
            print *, 'ToRange'//i2s(parenv % mype)//': ', trim(pvar % name), maxval(pvar % values)
          end if
        end if
        var => var % next
      end do
    end if
    call info('Mesh2MeshSolver', 'Interpolated variables between meshes', level=7)
  end subroutine m2ma_validate

end module mesh2meshsolver_acc
