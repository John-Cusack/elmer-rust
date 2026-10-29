! T022-A FilterTimeSeries port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/FilterTimeSeries.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape C_CHAR buffers, plain
! integers/reals). All Fortran-object traffic — TYPE(...) handles,
! OPTIONAL/POINTER dummies, assumed-shape descriptors, hidden lengths,
! allocatable-string results — stays inside Fortran, so the Rust side is
! pure C calls. The wrapper module mirrors the original's SAVE state
! (acc_dt; the entry-level Model/Solver handles; the current Var/MeanVar
! pointers with identical per-iteration lifetimes).
!
! All floating-point transcendental work (EXP, SIN, COS, SQRT, `**`,
! 2*PI*freq products, INT ratio-boundary checks) and every whole-array
! update stays in Fortran wrappers (same compiler, same evaluation
! order), so bit-identity holds by construction wherever the call
! sequence matches. Plain `+ - * /` weight arithmetic lives in Rust
! (IEEE-754, identical order).
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module filtertimeseries_acc
  use iso_c_binding
  use defutils
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(variable_t), pointer, save :: acc_var => null()
  type(variable_t), pointer, save :: acc_mean => null()
  real(c_double), pointer, save :: acc_meanfield(:) => null()
  real(c_double), save :: acc_dt = 0.0_c_double

contains

  ! -- entry state ----------------------------------------------------
  subroutine fta_begin(h_model, h_solver, dt) bind(c, name="fta_begin")
    type(c_ptr), value :: h_model, h_solver
    real(c_double), value :: dt
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    acc_dt = dt
  end subroutine fta_begin

  subroutine fta_info_banner() bind(c, name="fta_info_banner")
    call Info('FilterTimeSeries', '-----------------------------------------', level=4)
    call Info('FilterTimeSeries', 'Filtering variables over time', level=4)
    call Info('FilterTimeSeries', '-----------------------------------------', level=4)
  end subroutine fta_info_banner

  ! -- small helpers --------------------------------------------------
  subroutine fta_unpack(key, keylen, kf)
    character(kind=c_char), intent(in) :: key(*)
    integer(c_int), intent(in) :: keylen
    character(len=MAX_NAME_LEN, kind=c_char), intent(out) :: kf
    integer :: i
    kf = ''
    do i = 1, min(keylen, MAX_NAME_LEN)
      kf(i:i) = key(i)
    end do
  end subroutine fta_unpack

  ! -- ValueList reads on Solver % Values ------------------------------
  subroutine fta_get_string(key, keylen, buf, buflen, outlen, found) &
      bind(c, name="fta_get_string")
    character(kind=c_char), intent(in) :: key(*)
    integer(c_int), value :: keylen, buflen
    character(kind=c_char), intent(out) :: buf(*)
    integer(c_int), intent(out) :: outlen, found
    character(len=MAX_NAME_LEN, kind=c_char) :: kf
    character(len=:), allocatable :: s
    logical :: f
    integer :: i, n
    call fta_unpack(key, keylen, kf)
    s = ListGetString(acc_solver % values, trim(kf), f)
    found = merge(1, 0, f)
    n = min(len_trim(s), buflen)
    do i = 1, n
      buf(i) = s(i:i)
    end do
    outlen = n
  end subroutine fta_get_string

  subroutine fta_get_creal(key, keylen, v, found) bind(c, name="fta_get_creal")
    character(kind=c_char), intent(in) :: key(*)
    integer(c_int), value :: keylen
    real(c_double), intent(out) :: v
    integer(c_int), intent(out) :: found
    character(len=MAX_NAME_LEN, kind=c_char) :: kf
    logical :: f
    call fta_unpack(key, keylen, kf)
    v = GetCReal(acc_solver % values, trim(kf), f)
    found = merge(1, 0, f)
  end subroutine fta_get_creal

  subroutine fta_get_int(key, keylen, v, found) bind(c, name="fta_get_int")
    character(kind=c_char), intent(in) :: key(*)
    integer(c_int), value :: keylen
    integer(c_int), intent(out) :: v
    integer(c_int), intent(out) :: found
    character(len=MAX_NAME_LEN, kind=c_char) :: kf
    logical :: f
    call fta_unpack(key, keylen, kf)
    v = GetInteger(acc_solver % values, trim(kf), f)
    found = merge(1, 0, f)
  end subroutine fta_get_int

  subroutine fta_sim_creal(key, keylen, v) bind(c, name="fta_sim_creal")
    character(kind=c_char), intent(in) :: key(*)
    integer(c_int), value :: keylen
    real(c_double), intent(out) :: v
    character(len=MAX_NAME_LEN, kind=c_char) :: kf
    call fta_unpack(key, keylen, kf)
    v = GetCReal(acc_model % simulation, trim(kf))
  end subroutine fta_sim_creal

  function fta_realtime() result(t) bind(c, name="fta_realtime")
    real(c_double) :: t
    t = RealTime()
  end function fta_realtime

  function fta_first_order() result(exists) bind(c, name="fta_first_order")
    integer(c_int) :: exists
    logical :: f, l
    l = ListGetLogical(acc_solver % values, 'First Order', f)
    if (l) then
      exists = merge(1, 0, associated(acc_var % prevvalues))
    else
      exists = 0
    end if
  end function fta_first_order

  ! -- variable context ------------------------------------------------
  subroutine fta_get_time(t, found) bind(c, name="fta_get_time")
    real(c_double), intent(out) :: t
    integer(c_int), intent(out) :: found
    type(variable_t), pointer :: tv
    tv => VariableGet(acc_model % variables, 'Time')
    found = merge(1, 0, associated(tv))
    if (associated(tv)) t = tv % values(1)
  end subroutine fta_get_time

  subroutine fta_select_var(name, namelen, assoc) bind(c, name="fta_select_var")
    character(kind=c_char), intent(in) :: name(*)
    integer(c_int), value :: namelen
    integer(c_int), intent(out) :: assoc
    character(len=MAX_NAME_LEN, kind=c_char) :: nf
    call fta_unpack(name, namelen, nf)
    acc_var => VariableGet(acc_model % variables, trim(nf))
    assoc = merge(1, 0, associated(acc_var))
  end subroutine fta_select_var

  subroutine fta_var_info(inst_dofs, nsize, has_perm) bind(c, name="fta_var_info")
    integer(c_int), intent(out) :: inst_dofs, nsize, has_perm
    inst_dofs = acc_var % dofs
    nsize = size(acc_var % values)
    has_perm = merge(1, 0, associated(acc_var % perm))
  end subroutine fta_var_info

  subroutine fta_mean_select(name, namelen, assoc, mean_dofs, nsize) &
      bind(c, name="fta_mean_select")
    character(kind=c_char), intent(in) :: name(*)
    integer(c_int), value :: namelen
    integer(c_int), intent(out) :: assoc, mean_dofs, nsize
    character(len=MAX_NAME_LEN, kind=c_char) :: nf
    call fta_unpack(name, namelen, nf)
    acc_mean => VariableGet(acc_model % variables, trim(nf))
    assoc = merge(1, 0, associated(acc_mean))
    if (associated(acc_mean)) then
      acc_meanfield => acc_mean % values
      mean_dofs = acc_mean % dofs
      nsize = size(acc_var % values) * (acc_mean % dofs / acc_var % dofs)
    else
      mean_dofs = 0
      nsize = 0
    end if
  end subroutine fta_mean_select

  subroutine fta_mean_create(name, namelen, mean_dofs, nsize, nseries) &
      bind(c, name="fta_mean_create")
    character(kind=c_char), intent(in) :: name(*)
    integer(c_int), value :: namelen, mean_dofs, nsize, nseries
    character(len=MAX_NAME_LEN, kind=c_char) :: nf
    character(len=MAX_NAME_LEN) :: tmpname
    real(c_double), pointer :: comp(:)
    integer :: j
    call fta_unpack(name, namelen, nf)
    allocate(acc_meanfield(nsize))
    acc_meanfield = 0.0_c_double
    if (associated(acc_var % perm)) then
      call VariableAdd(acc_var % PrimaryMesh % variables, acc_var % PrimaryMesh, &
          acc_var % Solver, trim(nf), mean_dofs, acc_meanfield, acc_var % perm)
      if (mean_dofs > 1) then
        do j = 1, nseries
          tmpname = ComponentName(trim(nf), j)
          comp => acc_meanfield(j:nsize-nseries+j:nseries)
          call VariableAdd(acc_var % PrimaryMesh % variables, &
              acc_var % PrimaryMesh, acc_var % Solver, tmpname, 1, comp, &
              acc_var % perm)
        end do
      end if
    else
      call VariableAdd(acc_var % PrimaryMesh % variables, acc_var % PrimaryMesh, &
          acc_var % Solver, trim(nf), mean_dofs, acc_meanfield)
      if (mean_dofs > 1) then
        do j = 1, nseries
          tmpname = ComponentName(trim(nf), j)
          comp => acc_meanfield(j:nsize-nseries+j:nseries)
          call VariableAdd(acc_var % PrimaryMesh % variables, &
              acc_var % PrimaryMesh, acc_var % Solver, tmpname, 1, comp)
        end do
      end if
    end if
    acc_mean => VariableGet(acc_model % variables, trim(nf))
  end subroutine fta_mean_create

  ! -- transcendental / order-sensitive kernels (stay in Fortran) ------
  subroutine fta_decay(mode, dt, st, relax) bind(c, name="fta_decay")
    integer(c_int), value :: mode
    real(c_double), value :: dt, st
    real(c_double), intent(out) :: relax
    if (mode == 0) then
      relax = EXP(-dt/st)
    else
      relax = EXP(-1.0_c_double/st)
    end if
  end subroutine fta_decay

  subroutine fta_reset_check(st, t0, time, dt, do_reset) &
      bind(c, name="fta_reset_check")
    real(c_double), value :: st, t0, time, dt
    integer(c_int), intent(out) :: do_reset
    if (INT((time-t0)/st) /= INT((time-t0-dt)/st)) then
      ! The original zeroes the in-progress MeanField here; on the first
      ! visit that pointer is not yet associated (reference UB). Guard the
      ! store, keep the CumWeight reset decision identical.
      if (associated(acc_meanfield)) acc_meanfield = 0.0_c_double
      do_reset = 1
    else
      do_reset = 0
    end if
  end subroutine fta_reset_check

  subroutine fta_fcoeff(freq, time, t0, fcoeff) bind(c, name="fta_fcoeff")
    real(c_double), value :: freq, time, t0
    real(c_double), intent(out) :: fcoeff
    fcoeff = 2*PI*freq*(time-t0)
  end subroutine fta_fcoeff

  subroutine fta_filter_loop(oper_code, loopsize, instdofs, ratio, cold, cnew, &
      fcoeff, nsine, ncosine, has_prev, end_ratio) bind(c, name="fta_filter_loop")
    integer(c_int), value :: oper_code, loopsize, instdofs
    integer(c_int), value :: nsine, ncosine, has_prev, end_ratio
    real(c_double), value :: ratio, cold, cnew, fcoeff
    real(c_double), pointer :: inst(:), mean(:), prev(:,:)
    real(c_double) :: val, prevval
    integer :: i, j, k
    ! Exact transcription of the per-node update loop, including the
    ! `**2` vs `**2.0_dp` exponent spellings and the (2-ratio) blending.
    inst => acc_var % values
    mean => acc_meanfield
    if (has_prev /= 0) prev => acc_var % prevvalues
    do i = 1, loopsize
      select case (oper_code)
      case (1)
        val = inst(i)**2
        if (has_prev /= 0) prevval = prev(i,1)**2
      case (2)
        val = ABS(inst(i))
      case (3)
        val = 0.0_c_double
        do k = 1, instdofs
          val = val + inst(instdofs*(i-1)+k)**2.0_dp
        end do
        val = SQRT(val)
        if (has_prev /= 0) then
          prevval = 0.0_c_double
          do k = 1, instdofs
            prevval = prevval + prev(instdofs*(i-1)+k,1)**2.0_dp
          end do
          prevval = SQRT(prevval)
        end if
      case default
        val = inst(i)
        if (has_prev /= 0) prevval = prev(i,1)
      end select
      if (has_prev /= 0) then
        if (end_ratio /= 0) then
          val = (ratio*val + (2-ratio)*prevval)/2
        else
          val = (ratio*prevval + (2-ratio)*val)/2
        end if
      end if
      if (nsine + ncosine == 0) then
        mean(i) = cold*mean(i) + cnew*val
      else if (nsine > 0) then
        do j = 1, nsine
          mean(nsine*(i-1)+j) = cold*mean(nsine*(i-1)+j) + &
              cnew*SIN(j*fcoeff)*val
        end do
      else
        do j = 1, ncosine
          mean(ncosine*(i-1)+j) = cold*mean(ncosine*(i-1)+j) + &
              cnew*COS(j*fcoeff)*val
        end do
      end if
    end do
  end subroutine fta_filter_loop

  ! -- messages --------------------------------------------------------
  subroutine fta_warn_novar(name, namelen) bind(c, name="fta_warn_novar")
    character(kind=c_char), intent(in) :: name(*)
    integer(c_int), value :: namelen
    character(len=MAX_NAME_LEN, kind=c_char) :: nf
    call fta_unpack(name, namelen, nf)
    call Warn('SaveData', 'The desired variable '//TRIM(nf)//' does not exist!')
  end subroutine fta_warn_novar

  subroutine fta_fatal_freq() bind(c, name="fta_fatal_freq")
    call Fatal('FilterTimeSeries', 'Fourier series requires frequency')
  end subroutine fta_fatal_freq

  subroutine fta_fatal_scalar() bind(c, name="fta_fatal_scalar")
    call Fatal('FilterTimeSeries', &
        'Sine and Cosine series are implemented for scalars only!')
  end subroutine fta_fatal_scalar

  subroutine fta_mean_nullify() bind(c, name="fta_mean_nullify")
    nullify(acc_meanfield)
  end subroutine fta_mean_nullify

end module filtertimeseries_acc
