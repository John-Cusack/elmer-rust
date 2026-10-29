! SaveData pilot: statement-level BIND(C) wrappers for the Elmer calls in
! fem/src/modules/SaveData/SaveProjection.F90, including its internal
! ProjectToVariable (CRS projector traversal, normalization, ToSlave /
! ToMaster paths). The internal procedure's host association is mirrored
! by module SAVE handles set per variable call.
module sd_proj_acc
  use iso_c_binding
  use defutils
  use types
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(model_t), pointer, save :: acc_model => null()
  type(variable_t), pointer, save :: acc_var => null()
  type(variable_t), pointer, save :: acc_targetvar => null()
  integer, pointer, save :: acc_activeprojs(:) => null()
  logical, save :: acc_normalize = .true.
  logical, save :: acc_toslave = .false.
  logical, save :: acc_tomaster = .false.
  character(len=512), save :: acc_message

contains

  subroutine sdp_begin(h_model, h_solver) bind(c, name="sdp_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine sdp_begin

  subroutine sdp_init() bind(c, name="sdp_init")
    call ListAddNewString(acc_solver % values, 'Variable', &
        '-nooutput -global SaveProjection_var')
  end subroutine sdp_init

  subroutine sdp_info_start() bind(c, name="sdp_info_start")
    call Info('SaveProjection', 'Creating selected projected values as fields', level=8)
  end subroutine sdp_info_start

  subroutine sdp_gen_projectors() bind(c, name="sdp_gen_projectors")
    call GenerateProjectors(acc_model, acc_solver, Nonlinear = .true.)
  end subroutine sdp_gen_projectors

  function sdp_count_vars() result(n) bind(c, name="sdp_count_vars")
    integer(c_int) :: n
    type(valuelist_t), pointer :: sp
    character(len=4096) :: vname
    type(variable_t), pointer :: var
    logical :: found
    sp => GetSolverParams()
    n = 0
    do
      n = n + 1
      vname = ListGetString(sp, 'Variable ' // trim(sdp_i2s(n)), found)
      if (.not. found) exit
      var => VariableGet(acc_model % variables, trim(vname))
      if (.not. associated(var)) then
        call Warn('SaveProjection', 'Requested variable does not exist!')
        exit
      end if
    end do
    n = n - 1
  end function sdp_count_vars

  function sdp_i2s(i) result(s)
    integer :: i
    character(len=32) :: s
    write (s, '(I0)') i
    s = adjustl(s)
  end function sdp_i2s

  function sdp_project_var(idx) result(nrm2) bind(c, name="sdp_project_var")
    integer(c_int), value :: idx
    real(c_double) :: nrm2
    type(valuelist_t), pointer :: sp
    character(len=4096) :: varname, targetname
    type(variable_t), pointer :: var, targetvar
    integer, pointer :: unitperm(:)
    integer :: i, j, dofs
    logical :: found, normalize, toslave, tomaster
    sp => GetSolverParams()
    varname = ListGetString(sp, 'Variable ' // trim(sdp_i2s(idx)), found)
    var => VariableGet(acc_model % variables, trim(varname), thisonly=.true.)
    call Info('SaveProjection', 'Doing variable: ' // trim(varname), level=8)
    targetname = ListGetString(sp, 'Target Variable ' // trim(sdp_i2s(idx)), found)
    if (.not. found) targetname = 'Projection ' // trim(varname)
    normalize = ListGetLogical(sp, 'Normalize ' // trim(sdp_i2s(idx)), found)
    if (.not. found) normalize = .true.
    toslave = ListGetLogical(sp, 'Project To Slave ' // trim(sdp_i2s(idx)), found)
    tomaster = ListGetLogical(sp, 'Project To Master ' // trim(sdp_i2s(idx)), found)
    targetvar => VariableGet(acc_model % variables, trim(targetname), thisonly=.true.)
    if (.not. associated(targetvar)) then
      if (.not. associated(var % perm)) then
        nullify(unitperm)
        allocate(unitperm(size(var % values)))
        do j = 1, size(var % values)
          unitperm(j) = j
        end do
        call VariableAddVector(acc_model % mesh % variables, acc_solver % mesh, &
            acc_solver, trim(targetname), var % dofs, perm=unitperm)
      else
        call VariableAddVector(acc_model % mesh % variables, acc_solver % mesh, &
            acc_solver, trim(targetname), var % dofs, perm=var % perm, secondary=.true.)
      end if
      targetvar => VariableGet(acc_model % variables, trim(targetname))
    end if
    acc_var => var
    acc_targetvar => targetvar
    acc_normalize = normalize
    acc_toslave = toslave
    acc_tomaster = tomaster
    if (associated(acc_activeprojs)) nullify(acc_activeprojs)
    acc_activeprojs => ListGetIntegerArray(sp, &
        'Active Projectors ' // trim(sdp_i2s(idx)), found)
    call sdp_project_to_variable()
    nrm2 = sum(targetvar % values**2)
  end function sdp_project_var

  subroutine sdp_project_to_variable()
    type(matrix_t), pointer :: a
    integer :: bc, dofs, i, j, k, pi, pj, m, n
    integer, pointer :: rows(:), cols(:)
    logical :: acti, actj
    real(kind=dp) :: r1
    real(kind=dp), pointer :: values(:)
    real(kind=dp), allocatable :: weight(:)
    logical, pointer :: isinvinvperm(:)
    type(variable_t), pointer :: var, targetvar
    var => acc_var
    targetvar => acc_targetvar
    dofs = var % dofs
    targetvar % values = 0.0_dp
    if (acc_normalize) then
      allocate(weight(size(targetvar % values)))
      weight = 0.0_dp
    end if
    do bc = 1, acc_model % numberofbcs
      if (associated(acc_activeprojs)) then
        if (.not. any(acc_activeprojs == bc)) cycle
      end if
      a => CurrentModel % BCs(bc) % pmatrix
      if (.not. associated(a)) then
        a => acc_solver % mortarbcS(bc) % projector
      end if
      if (.not. associated(a)) cycle
      n = a % numberofrows
      if (n == 0) cycle
      write (acc_message, '(A,I0,A,I0)') 'Doing projection for BC ', bc, &
          ' of size ', n
      call Info('SaveProjection', acc_message, level=20)
      rows => a % rows
      cols => a % cols
      values => a % values
      if (.not. associated(a % invperm)) then
        call Fatal('SaveProjection', 'InvPerm not associated!')
      end if
      m = maxval(cols)
      allocate(isinvinvperm(m))
      isinvinvperm = .false.
      do i = 1, n
        pi = a % invperm(i)
        if (pi > 0 .and. pi <= m) isinvinvperm(pi) = .true.
      end do
      do k = 1, dofs
        do i = 1, n
          pi = a % invperm(i)
          if (pi == 0 .or. pi > m) cycle
          acti = isinvinvperm(pi)
          if (associated(var % perm)) then
            if (pi > size(var % perm)) cycle
            pi = var % perm(pi)
            if (pi == 0) cycle
          end if
          pi = dofs * (pi - 1) + k
          r1 = 0.0_dp
          do j = rows(i), rows(i + 1) - 1
            pj = cols(j)
            if (associated(var % perm)) then
              if (pj > size(var % perm)) cycle
              pj = var % perm(pj)
              if (pj == 0) cycle
            end if
            pj = dofs * (pj - 1) + k
            actj = isinvinvperm(cols(j))
            if (acc_toslave) then
              if (.not. actj) then
                r1 = -values(j) * var % values(pj)
                targetvar % values(pi) = targetvar % values(pi) + r1
              end if
            else if (.not. acc_tomaster) then
              r1 = -values(j) * var % values(pj)
              targetvar % values(pi) = targetvar % values(pi) + r1
            end if
            if (acc_normalize .and. actj) weight(pi) = weight(pi) + values(j)
            if (acc_tomaster) then
              if (.not. actj) then
                r1 = values(j) * var % values(pi)
                targetvar % values(pj) = targetvar % values(pj) + r1
                if (acc_normalize .and. .not. actj) weight(pj) = weight(pj) + values(j)
              end if
            end if
          end do
        end do
      end do
      deallocate(isinvinvperm)
    end do
    if (acc_normalize) then
      where (abs(weight) > epsilon(r1))
        targetvar % values = targetvar % values / weight
      end where
    end if
  end subroutine sdp_project_to_variable

  subroutine sdp_set_norm(nrm) bind(c, name="sdp_set_norm")
    real(c_double), value :: nrm
    if (size(acc_solver % variable % values) == 1) then
      acc_solver % variable % values = nrm
    end if
    write (acc_message, '(A,ES12.3)') 'Combined L2 norm of all projected fields: ', nrm
    call Info('SaveProjection', acc_message)
  end subroutine sdp_set_norm

end module sd_proj_acc
