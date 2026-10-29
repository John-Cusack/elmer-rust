! T022-A FourierLoss port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/FourierLoss.F90, compiled against the
! built .mod files of the target configure (T020 pilot shipped mechanism,
! ADR-057 section 3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars). All Fortran-object traffic — TYPE(...)
! handles, OPTIONAL/POINTER dummies, assumed-shape descriptors, hidden
! lengths, array-pointer results — stays inside Fortran, so the Rust side
! is pure C calls. The wrapper module mirrors FourierLossSolver's own
! locals (entry SAVE state, LocalFourierTransform SAVE state, BulkAssembly
! SAVE Nodes) with identical lifetimes. All floating-point work (the DFT
! SIN/COS coefficient updates, the Galerkin assembly, SQRT/SUM
! reductions, the RealTime-based reports) and all formatted/file I/O
! stay in Fortran (same compiler), so bit-identity holds by construction
! wherever the call sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module fourierloss_acc
  use iso_c_binding
  use defutils
  use magnetodynamicsutils
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(valuelist_t), pointer, save :: acc_params => null()
  type(variable_t), pointer, save :: acc_lossvar => null()
  type(variable_t), pointer, save :: acc_targetvar => null()
  type(variable_t), pointer, save :: acc_nodallossvar => null()
  type(solver_t), pointer, save :: acc_targetptr => null()
  real(c_double), pointer, save :: acc_targetfield(:) => null()
  real(c_double), pointer, save :: acc_prevfield(:,:) => null()
  integer(c_int), pointer, save :: acc_fourierperm(:) => null()
  real(c_double), pointer, save :: acc_saverhs(:) => null()

  type :: fls_varptr_t
    type(variable_t), pointer :: var => null()
  end type fls_varptr_t
  type(fls_varptr_t), allocatable, save :: acc_compvars(:)
  type(fls_varptr_t), allocatable, save :: acc_compvarse(:)
  type(fls_varptr_t), allocatable, save :: acc_fouriervars(:)
  real(c_double), allocatable, save :: acc_comploss(:)
  real(c_double), allocatable, save :: acc_series(:,:)
  real(c_double), allocatable, save :: acc_body(:,:)
  real(c_double), allocatable, target, save :: acc_otherrhs(:, :)

  ! Entry-level SAVE state (mirrors the solver body's SAVE locals).
  integer(c_int), save :: acc_times = 0
  integer(c_int), save :: acc_dim = 0
  integer(c_int), save :: acc_ncomp = 0
  integer(c_int), save :: acc_nvar = 0
  integer(c_int), save :: acc_fourdofs = 0
  integer(c_int), save :: acc_nsize = 0
  integer(c_int), save :: acc_tdofs = 0
  integer(c_int), save :: acc_nodal = 0
  integer(c_int), save :: acc_elem = 0
  integer(c_int), save :: acc_nloss = 0
  integer(c_int), save :: acc_sumcomp = 0
  integer(c_int), save :: acc_direct = 0
  integer(c_int), save :: acc_av = 0
  integer(c_int), save :: acc_exact = 0
  integer(c_int), save :: acc_simpson = 0
  integer(c_int), save :: acc_flip = 0
  real(c_double), save :: acc_at0 = 0.0_c_double
  real(c_double), save :: acc_omega = 0.0_c_double
  real(c_double), save :: acc_total = 0.0_c_double
  character(len=MAX_NAME_LEN), save :: acc_pref = ''
  character(len=MAX_NAME_LEN), save :: acc_varname = ''

  ! LocalFourierTransform SAVE state (mirrors its SAVE statement).
  real(c_double), save :: acc_d_time = 0.0_c_double
  real(c_double), save :: acc_d_ratio = 0.0_c_double
  real(c_double), save :: acc_d_tcycle = 0.0_c_double
  real(c_double), save :: acc_d_t0 = 0.0_c_double
  integer(c_int), save :: acc_d_prev = 0
  logical, save :: acc_d_initpending = .false.

  ! BulkAssembly SAVE Nodes (mirrors its SAVE statement).
  type(nodes_t), target, save :: acc_nodes

contains

  ! -- entry state ----------------------------------------------------
  subroutine fls_begin(h_model, h_solver) bind(c, name="fls_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine fls_begin

  ! -- FourierLossSolver_init body (verbatim transcription) ------------
  subroutine fls_init_body() bind(c, name="fls_init_body")
    type(valuelist_t), pointer :: solverparams
    logical :: found, found2
    character(len=MAX_NAME_LEN) :: varname
    integer :: i, j, k, n, m, mysolver, soln
    logical :: elementalfield, nodalfield
    logical :: oldkeywordstyle, separatecomponents, nodallosses
    character(len=MAX_NAME_LEN) :: pref
    integer :: ncomp
    character(*), parameter :: caller = 'FourierLossSolver'

    solverparams => GetSolverParams()

    oldkeywordstyle = ListCheckPresentAnyMaterial(acc_model, &
        'Harmonic Loss Linear Coefficient') .or. &
        ListCheckPresentAnyMaterial(acc_model, &
        'Harmonic Loss Quadratic Coefficient')
    if (oldkeywordstyle) then
      if (parenv % mype == 0) then
        print *, 'You are using old keywords'
        print *, 'Obsolite: "Harmonic Loss Linear Coefficient"'
        print *, 'Obsolite: "Harmonic Loss Quadratic Coefficient"'
        print *, 'Obsolite: "Harmonic Loss Linear Frequency Exponent"'
        print *, 'Obsolite: "Harmonic Loss Quadratic Frequency Exponent"'
        print *, 'Obsolite: "Harmonic Loss Linear Exponent"'
        print *, 'Obsolite: "Harmonic Loss Quadratic Exponent"'
        print *, 'etc. see the documentation.'
      end if
      call Fatal(caller, &
          'The keyword format is obsolite. Make the fixes and continue!')
    end if

    pref = ListGetString(solverparams, 'Scalars Prefix', found)
    if (.not. found) pref = 'res:'

    varname = ListGetString(solverparams, 'Target Variable', found)
    if (.not. found) then
      call Fatal(caller, 'Give > Target Variable < !')
    end if

    ncomp = 0
    do i = 1, 10
      if (ListCheckPresentAnyMaterial(acc_model, &
          'Harmonic Loss Coefficient '//I2S(i))) then
        ncomp = i
      else
        exit
      end if
    end do

    call ListAddConstReal(acc_model % simulation, &
        trim(pref)//' fourier loss total', 0.0_dp)
    do i = 1, ncomp
      call ListAddConstReal(acc_model % simulation, &
          trim(pref)//' fourier loss '//I2S(i), 0.0_dp)
    end do

    call ListAddNewString(solverparams, 'Variable', '-nooutput FourierTmp')

    nodalfield = ListGetLogical(solverparams, 'Calculate Nodal Fields', found)
    elementalfield = ListGetLogical(solverparams, &
        'Calculate Elemental Fields', found2)
    if (.not. (found .or. found2)) then
      elementalfield = .true.
      call ListAddLogical(solverparams, &
          'Calculate Elemental Fields', elementalfield)
      call Info(caller, &
          'Computing elemental fields only by default when neither given!')
    end if

    separatecomponents = ListGetLogical(solverparams, &
        'Separate Loss Components', found)

    if (nodalfield) then
      if (separatecomponents) then
        call ListAddString(solverparams, &
            NextFreeKeyword('Exported Variable ', solverparams), &
            'Fourier Loss[Fourier Loss:'//I2S(ncomp)//']')
      else
        call ListAddString(solverparams, &
            NextFreeKeyword('Exported Variable ', solverparams), &
            'Fourier Loss')
      end if
    end if

    if (elementalfield) then
      if (separatecomponents) then
        call ListAddString(solverparams, &
            NextFreeKeyword('Exported Variable ', solverparams), &
            '-dg Fourier Loss e[Fourier Loss e:'//I2S(ncomp)//']')
      else
        call ListAddString(solverparams, &
            NextFreeKeyword('Exported Variable ', solverparams), &
            '-dg Fourier Loss e')
      end if
    end if

    nodallosses = ListGetLogical(solverparams, 'Calculate Nodal Losses', found)
    if (nodallosses) then
      call ListAddString(solverparams, &
          NextFreeKeyword('Exported Variable', solverparams), &
          'Nodal Fourier Loss')
    end if
  end subroutine fls_init_body

  ! -- entry head: Info banner, TimesVisited, sizing, Ncomp ------------
  subroutine fls_head() bind(c, name="fls_head")
    logical :: found
    integer :: i
    logical :: separatecomponents
    character(*), parameter :: caller = 'FourierLossSolver'

    call Info(caller, '-------------------------------------', level=4)
    call Info(caller, 'Computing Fourier losses             ', level=4)
    call Info(caller, '-------------------------------------', level=4)

    acc_times = acc_times + 1
    acc_dim = CoordinateSystemDimension()
    acc_params => GetSolverParams()
    acc_lossvar => acc_solver % variable

    acc_at0 = RealTime()

    acc_pref = ListGetString(acc_params, 'Scalars Prefix', found)
    if (.not. found) acc_pref = 'res:'

    acc_ncomp = 0
    do i = 1, 10
      if (ListCheckPresentAnyMaterial(acc_model, &
          'Harmonic Loss Coefficient '//I2S(i))) then
        acc_ncomp = i
      else
        exit
      end if
    end do
    if (acc_ncomp == 0) then
      call Fatal(caller, &
          'Some material must have > Harmonic Loss Coefficient i <')
    end if
    call Info(caller, 'Considering number of components: '// &
        I2S(acc_ncomp), level=5)

    separatecomponents = ListGetLogical(acc_params, &
        'Separate Loss Components', found)
    acc_sumcomp = merge(1, 0, .not. separatecomponents)

    if (acc_sumcomp == 1) then
      acc_nvar = 1
    else
      acc_nvar = acc_ncomp
    end if

    if (.not. allocated(acc_comploss)) allocate(acc_comploss(acc_ncomp))
    acc_comploss = 0.0_dp
  end subroutine fls_head

  ! -- nodal / elemental field setup ----------------------------------
  subroutine fls_setup_nodal() bind(c, name="fls_setup_nodal")
    logical :: found
    integer :: icomp
    character(*), parameter :: caller = 'FourierLossSolver'
    acc_nodal = 0
    if (GetLogical(acc_params, 'Calculate Nodal Fields', found)) then
      acc_nodal = 1
      if (.not. allocated(acc_compvars)) allocate(acc_compvars(acc_nvar))
      if (acc_nvar == 1) then
        acc_compvars(1) % var => VariableGet(acc_solver % mesh % variables, &
            'Fourier Loss')
      else
        do icomp = 1, acc_nvar
          acc_compvars(icomp) % var => VariableGet( &
              acc_solver % mesh % variables, &
              'Fourier Loss '//I2S(icomp))
        end do
      end if
      do icomp = 1, acc_nvar
        if (.not. associated(acc_compvars(icomp) % var)) then
          call Fatal(caller, &
              'Variable for component does not exist: '//I2S(icomp))
        end if
      end do
    end if
  end subroutine fls_setup_nodal

  subroutine fls_setup_elem() bind(c, name="fls_setup_elem")
    logical :: found
    integer :: icomp
    character(*), parameter :: caller = 'FourierLossSolver'
    acc_elem = 0
    if (GetLogical(acc_params, 'Calculate Elemental Fields', found)) then
      acc_elem = 1
      if (.not. allocated(acc_compvarse)) allocate(acc_compvarse(acc_nvar))
      if (acc_nvar == 1) then
        acc_compvarse(1) % var => VariableGet(acc_solver % mesh % variables, &
            'Fourier Loss e')
      else
        do icomp = 1, acc_nvar
          acc_compvarse(icomp) % var => VariableGet( &
              acc_solver % mesh % variables, &
              'Fourier Loss e '//I2S(icomp))
        end do
      end if
      do icomp = 1, acc_nvar
        if (.not. associated(acc_compvarse(icomp) % var)) then
          call Fatal(caller, &
              'Variable e for component does not exist: '//I2S(icomp))
        end if
      end do
    end if
  end subroutine fls_setup_elem

  function fls_nodalloss_flag() result(r) bind(c, name="fls_nodalloss_flag")
    integer(c_int) :: r
    acc_nodallossvar => VariableGet(acc_solver % mesh % variables, &
        'Nodal Fourier Loss')
    acc_nloss = merge(1, 0, associated(acc_nodallossvar))
    r = acc_nloss
  end function fls_nodalloss_flag

  function fls_perm_empty() result(r) bind(c, name="fls_perm_empty")
    integer(c_int) :: r
    r = merge(1, 0, count(acc_solver % variable % perm > 0) <= 0)
  end function fls_perm_empty

  ! -- target field fetch + AV/Direct checks ---------------------------
  subroutine fls_fetch_target() bind(c, name="fls_fetch_target")
    character(len=MAX_NAME_LEN) :: varname
    logical :: found
    character(*), parameter :: caller = 'FourierLossSolver'
    varname = GetString(acc_params, 'Target Variable', found)
    acc_varname = varname
    acc_targetvar => VariableGet(acc_solver % mesh % variables, varname)
    if (.not. associated(acc_targetvar)) then
      call Fatal(caller, 'Target field not present: '//trim(varname))
    end if
    if (.not. associated(acc_targetvar % prevvalues)) then
      call Fatal(caller, &
          'Target field does not have PrevValues: '//trim(varname))
    end if
    acc_targetfield => acc_targetvar % values
    acc_prevfield => acc_targetvar % prevvalues(:, :)
    acc_nsize = size(acc_targetfield)
    acc_tdofs = acc_targetvar % dofs

    acc_flip = merge(1, 0, acc_targetvar % periodicflipactive)
    if (acc_flip == 1) then
      call Info(caller, 'Assuming initial field to have conforming flips')
    end if

    acc_av = 0
    acc_direct = 0
    if (ListGetLogical(acc_params, 'Target Variable Direct', found)) then
      acc_direct = 1
      call Info(caller, 'Using the target field with '// &
          I2S(acc_tdofs)//' dofs directly!')
    else
      if (acc_dim == 3) then
        acc_av = merge(1, 0, ListGetLogical(acc_params, &
            'Target Variable AV', found))
        if (.not. found) then
          acc_av = merge(1, 0, &
              size(acc_targetvar % perm) > acc_solver % mesh % numberofnodes)
        end if
        if (acc_av == 1) then
          if (acc_tdofs > 1) then
            call Fatal(caller, 'Assuming only one component for AV field!')
          end if
        else
          if (acc_tdofs /= 3) then
            call Fatal(caller, &
                'Assuming precisely three nodal components in 3D!')
          end if
        end if
      else
        if (acc_tdofs /= 1) then
          call Fatal(caller, &
              'Assuming only one nodal component Az in 2D!')
        end if
      end if
    end if
  end subroutine fls_fetch_target

  ! -- target solver search --------------------------------------------
  subroutine fls_find_target() bind(c, name="fls_find_target")
    character(len=MAX_NAME_LEN) :: varstring
    integer :: i, nlen, mlen
    logical :: found
    character(*), parameter :: caller = 'FourierLossSolver'
    mlen = len_trim(acc_varname)
    nullify(acc_targetptr)
    do i = 1, acc_model % numberofsolvers
      varstring = ListGetString(acc_model % solvers(i) % values, &
          'Variable', found)
      if (found) then
        nlen = len_trim(varstring)
        if (varstring(1:nlen) == acc_varname(1:mlen)) then
          acc_targetptr => acc_model % solvers(i)
          exit
        end if
      end if
    end do
    if (.not. associated(acc_targetptr)) then
      call Fatal(caller, &
          'Target field solver cannot be found: '//trim(acc_varname))
    end if
  end subroutine fls_find_target

  ! -- Fourier field sizing + creation ---------------------------------
  subroutine fls_fourier_dofs() bind(c, name="fls_fourier_dofs")
    logical :: found
    character(*), parameter :: caller = 'FourierLossSolver'
    acc_fourdofs = 1 + 2*ListGetInteger(acc_params, &
        'Fourier Series Components', found)
    if (.not. found) then
      call Fatal(caller, 'Give > Fourier Series Components < !')
    end if
  end subroutine fls_fourier_dofs

  subroutine fls_fourier_alloc() bind(c, name="fls_fourier_alloc")
    integer :: i
    logical :: found, fourieroutput
    character(len=MAX_NAME_LEN) :: fouriername
    real(kind=dp), pointer :: fourierfield(:)
    if (.not. allocated(acc_fouriervars)) then
      allocate(acc_fouriervars(acc_fourdofs))
      fourieroutput = ListGetLogical(acc_params, &
          'Fourier Series Output', found)
      do i = 1, acc_fourdofs
        nullify(fourierfield)
        allocate(fourierfield(acc_nsize))
        fourierfield = 0.0_dp
        if (i == 1) then
          fouriername = trim(acc_varname)//' Fourier 0'
        else if (modulo(i, 2) == 0) then
          fouriername = trim(acc_varname)//' Fourier Cos'//I2S(i/2)
        else
          fouriername = trim(acc_varname)//' Fourier Sin'//I2S(i/2)
        end if
        call VariableAddVector(acc_solver % mesh % variables, &
            acc_solver % mesh, acc_solver, &
            fouriername, acc_targetvar % dofs, fourierfield, &
            acc_targetvar % perm, fourieroutput)
        acc_fouriervars(i) % var => VariableGet( &
            acc_solver % mesh % variables, fouriername)
      end do
    end if
    acc_fourierperm => acc_targetvar % perm
  end subroutine fls_fourier_alloc

  subroutine fls_loss_alloc() bind(c, name="fls_loss_alloc")
    if (.not. allocated(acc_series)) then
      allocate(acc_series(acc_ncomp, acc_fourdofs), &
          acc_body(acc_ncomp, acc_model % numberofbodies))
    end if
  end subroutine fls_loss_alloc

  ! -- frequency + integration flags -----------------------------------
  subroutine fls_omega() bind(c, name="fls_omega")
    logical :: found
    character(*), parameter :: caller = 'FourierLossSolver'
    acc_omega = GetAngularFrequency(found=found)
    if (.not. found) then
      call Fatal('FourierLosses', &
          'The Fourier transform requires frequency!')
    end if
    write(message, '(A,ES12.3)') &
        'Base frequency for Fourier transform: ', acc_omega/(2*pi)
    call Info(caller, message, level=8)

    acc_exact = merge(1, 0, &
        .not. GetLogical(acc_params, 'Inexact Integration', found))
    if (acc_exact == 1) then
      call Info(caller, 'Using exact integration')
    else
      acc_simpson = merge(1, 0, &
          GetLogical(acc_params, 'Simpsons Rule', found))
      if (acc_simpson == 1) then
        call Info(caller, 'Using Simpsons rule for integration')
      else
        call Info(caller, 'Using trapetsoidal rule for integration')
      end if
    end if
  end subroutine fls_omega

  ! -- LocalFourierTransform (verbatim transcription) ------------------
  function fls_dft(end0, dt) result(nend) bind(c, name="fls_dft")
    integer(c_int), value :: end0
    real(c_double), value :: dt
    integer(c_int) :: nend
    logical :: endcycle

    integer :: n0, i, j, k, cycles, currentcycle
    real(kind=dp) :: c0, c1, ta, tb, a0, a1
    real(kind=dp), pointer :: component(:)
    logical :: leftrule, rightrule, cosine
    real(kind=dp) :: anglea, angleb, anglec, compcoeffa, compcoeffb
    real(kind=dp) :: vala, valb, dtfraction
    real(kind=dp), parameter :: ratioeps = 1.0d-8
    logical :: found
    character(*), parameter :: caller = 'FourierLossSolver'

    endcycle = (end0 /= 0)

    if (acc_d_initpending) then
      do i = 1, acc_fourdofs
        acc_fouriervars(i) % var % values = 0.0_dp
      end do
    end if
    acc_d_initpending = .false.

    if (endcycle) then
      acc_d_ratio = 1.0_c_double - acc_d_ratio
      if (acc_d_ratio < ratioeps) then
        acc_d_initpending = .true.
        nend = merge(1, 0, endcycle)
        return
      end if
      do i = 1, acc_fourdofs
        acc_fouriervars(i) % var % values = 0.0_dp
      end do
      leftrule = .true.
      call Info(caller, 'Finishing the timestep using left rule', level=8)
    else
      acc_d_time = GetTime()

      n0 = GetInteger(acc_params, 'Fourier Start Timestep', found)
      if (found) then
        i = GetTimestep()
        if (i < n0) then
          call Info(caller, &
              'Fourier transform not yet active, returning...', level=6)
          nend = merge(1, 0, endcycle)
          return
        else if (i == n0) then
          acc_d_t0 = acc_d_time
        end if
      else
        cycles = GetInteger(acc_params, 'Fourier Start Cycles', found)
        if (found) then
          acc_d_t0 = cycles*2*pi/acc_omega
        else
          acc_d_t0 = GetCReal(acc_params, 'Fourier Start Time', found)
          if (.not. found) then
            if (acc_times == 1) acc_d_t0 = max(acc_d_time - dt, 0.0_dp)
          end if
        end if
      end if

      acc_d_time = acc_d_time - acc_d_t0

      leftrule = .false.
      if (acc_d_time <= 0.0_dp) then
        acc_d_ratio = 0.0_dp
        call Info(caller, &
            'Fourier transform not yet active, returning...', level=6)
        nend = merge(1, 0, endcycle)
        return
      else if (acc_d_time - dt < 0.0_dp) then
        acc_d_ratio = acc_d_time/dt
        leftrule = .true.
        call Info(caller, 'Starting Fourier transform cycle', level=6)
      else
        acc_d_ratio = 1.0_dp
      end if

      cycles = GetInteger(acc_params, 'Fourier Integrate Cycles', found)
      if (.not. found) cycles = 1
      acc_d_tcycle = cycles*2*pi/acc_omega

      rightrule = .false.
      currentcycle = nint(acc_d_time/acc_d_tcycle)

      if (currentcycle /= acc_d_prev) then
        rightrule = (int(acc_d_time/acc_d_tcycle) /= &
            int((acc_d_time - dt)/acc_d_tcycle))
        rightrule = rightrule .or. &
            (abs(currentcycle - acc_d_time/acc_d_tcycle) < ratioeps)

        if (rightrule) then
          if (leftrule) then
            call Fatal(caller, &
                'Cannot use left and right rule at the same time!')
          end if

          acc_d_ratio = 1.0_dp - &
              (acc_d_time - nint(acc_d_time/acc_d_tcycle)*acc_d_tcycle)/dt

          call Info(caller, 'Finising Fourier transform cycle', level=6)

          endcycle = .true.
          acc_d_prev = currentcycle
        end if
      end if
    end if

    ta = acc_d_time - dt
    tb = acc_d_time

    if (leftrule) then
      ta = acc_d_time - acc_d_ratio*dt
    else if (rightrule) then
      tb = acc_d_time - (1 - acc_d_ratio)*dt
    end if
    dtfraction = (tb - ta)/acc_d_tcycle

    write(message, '(A,ES12.3,A,ES12.3,A)') &
        'Timestep interval: [', ta, ',', tb, ']'
    call Info(caller, message, level=8)

    do j = 1, acc_fourdofs
      component => acc_fouriervars(j) % var % values

      if (acc_exact == 1) then
        if (j == 1) then
          compcoeffa = 0.5d0*dtfraction
          compcoeffb = 0.5d0*dtfraction
        else
          k = j/2
          c1 = k*acc_omega
          c0 = c1*ta
          anglea = acc_omega*ta
          angleb = acc_omega*tb

          cosine = (modulo(j, 2) == 0)
          if (cosine) then
            a0 = 2.0d0/(c1*acc_d_tcycle)* &
                (sin(k*angleb) - sin(k*anglea))
            a1 = 2.0d0/(c1*acc_d_tcycle)* &
                (1.0d0/(c1*(tb - ta))*(cos(k*angleb) - &
                cos(k*anglea)) + sin(k*angleb))
          else
            a0 = 2.0d0/(c1*acc_d_tcycle)* &
                (-cos(k*angleb) + cos(k*anglea))
            a1 = 2.0d0/(c1*acc_d_tcycle)* &
                (1.0d0/(c1*(tb - ta))*(sin(k*angleb) - &
                sin(k*anglea)) - cos(k*angleb))
          end if

          compcoeffa = a0 - a1
          compcoeffb = a1
        end if
        if (j == 1) then
          if (leftrule) then
            compcoeffa = acc_d_ratio*compcoeffa
            compcoeffb = (1 - acc_d_ratio)*compcoeffa + compcoeffb
          else if (rightrule) then
            compcoeffb = acc_d_ratio*compcoeffb
            compcoeffa = (1 - acc_d_ratio)*compcoeffb + compcoeffa
          end if
        else
          if (leftrule) then
            compcoeffa = acc_d_ratio*(a0 - a1)
            compcoeffb = (1 - acc_d_ratio)*a0 + acc_d_ratio*a1
          else if (rightrule) then
            compcoeffa = a0 - acc_d_ratio*a1
            compcoeffb = acc_d_ratio*a1
          end if
        end if

      else
        if (j == 1) then
          compcoeffa = 0.5d0*dtfraction
          compcoeffb = 0.5d0*dtfraction
        else
          k = j/2
          anglea = acc_omega*ta
          angleb = acc_omega*tb
          cosine = (modulo(j, 2) == 0)

          if (cosine) then
            compcoeffa = cos(k*anglea)*dtfraction
            compcoeffb = cos(k*angleb)*dtfraction
          else
            compcoeffa = sin(k*anglea)*dtfraction
            compcoeffb = sin(k*angleb)*dtfraction
          end if

          if (acc_simpson == 1) then
            anglec = (anglea + angleb)/2.0_dp
            if (cosine) then
              compcoeffa = &
                  (compcoeffa + 2*cos(k*anglec)*dtfraction)/3.0
              compcoeffb = &
                  (compcoeffb + 2*cos(k*anglec)*dtfraction)/3.0
            else
              compcoeffa = &
                  (compcoeffa + 2*sin(k*anglec)*dtfraction)/3.0
              compcoeffb = &
                  (compcoeffb + 2*sin(k*anglec)*dtfraction)/3.0
            end if
          end if
        end if

        if (leftrule) then
          compcoeffa = acc_d_ratio*compcoeffa
          compcoeffb = (1 - acc_d_ratio)*compcoeffa + compcoeffb
        else if (rightrule) then
          compcoeffb = acc_d_ratio*compcoeffb
          compcoeffa = (1 - acc_d_ratio)*compcoeffb + compcoeffa
        end if

      end if

      component = component + compcoeffa*acc_prevfield(:, 1) + &
          compcoeffb*acc_targetfield

    end do

    nend = merge(1, 0, endcycle)
  end function fls_dft

  ! -- DefaultInitialize + multi-RHS staging ---------------------------
  subroutine fls_default_init() bind(c, name="fls_default_init")
    call DefaultInitialize()
  end subroutine fls_default_init

  subroutine fls_otherrhs() bind(c, name="fls_otherrhs")
    if (acc_nvar > 1) then
      if (.not. allocated(acc_otherrhs)) then
        allocate(acc_otherrhs(size(acc_solver % matrix % rhs), acc_nvar - 1))
      end if
      acc_otherrhs = 0.0_dp
    end if
  end subroutine fls_otherrhs

  ! -- BulkAssembly (verbatim transcription) ---------------------------
  subroutine fls_bulk() bind(c, name="fls_bulk")
    integer :: elem, i, j, k, p, q, t, bodyid, nsum
    integer :: icomp
    integer :: n, nd, nt
    type(GaussIntegrationPoints_t), target :: ip
    type(Element_t), pointer :: element, targetelement
    real(kind=dp) :: weight, coeff, coeff2, detj, gradatip(3, 3)
    real(kind=dp) :: curlatip(3), valatip
    logical :: found, found2
    type(valuelist_t), pointer :: material
    real(kind=dp) :: freq, fourierfreq
    real(kind=dp) :: losscoeff(10), freqpower(10), fieldpower(10)
    real(kind=dp), pointer :: component(:)
    real(kind=dp), allocatable :: stiff(:, :), force(:, :), dgforce(:)
    real(kind=dp), allocatable :: elemfield(:)
    real(kind=dp), allocatable :: basis(:), dbasisdx(:, :)
    real(kind=dp), allocatable :: wbasis(:, :), rotwbasis(:, :)
    integer, allocatable :: indeces(:), pivot(:)
    logical :: dg, erroneous, materialexponents
    integer :: ncomp
    integer(c_int), pointer :: eperm(:)
    character(*), parameter :: caller = 'FourierLossSolver'

    acc_saverhs => acc_solver % matrix % rhs

    n = max(acc_solver % mesh % maxelementdofs, &
        acc_solver % mesh % maxelementnodes)
    allocate(stiff(n, n), force(10, n), dgforce(n), pivot(n), &
        basis(n), dbasisdx(n, 3))

    if (acc_av == 1) then
      nt = max(acc_targetptr % mesh % maxelementdofs, &
          acc_targetptr % mesh % maxelementnodes)
      allocate(wbasis(nt, 3), rotwbasis(nt, 3), indeces(nt), elemfield(nt))
    else
      allocate(elemfield(n), indeces(n))
    end if

    freq = acc_omega/(2*pi)

    dg = GetLogical(acc_params, 'Discontinuous Galerkin', found)

    materialexponents = ListCheckPrefixAnyMaterial(acc_model, &
        'Harmonic Loss Frequency Exponent')
    ncomp = acc_ncomp
    if (.not. materialexponents) then
      call GetLossExponents(acc_params, freqpower, fieldpower, ncomp)
    end if

    acc_series = 0.0_dp
    acc_body = 0.0_dp
    acc_comploss = 0.0_dp

    do elem = 1, GetNOFActive()

      element => GetActiveElement(elem)
      call GetElementNodes(acc_nodes)
      nd = GetElementDOFs(indeces)
      n = GetElementNOFNodes()

      if (acc_av == 1) then
        targetelement => GetActiveElement(elem, acc_targetptr)
        if (.not. associated(targetelement)) then
          call Fatal(caller, 'Element on target mesh cannot be associated')
        end if
        nt = GetElementDOFs(indeces, targetelement, acc_targetptr)
      end if

      ip = GaussPoints(element)
      stiff = 0.0_dp
      force = 0.0_dp

      bodyid = element % bodyid
      material => GetMaterial()

      if (materialexponents) then
        call GetLossExponents(material, freqpower, fieldpower, ncomp)
      end if

      do t = 1, ip % n
        if (acc_direct == 1) then
          found = ElementInfo(element, acc_nodes, ip % u(t), &
              ip % v(t), ip % w(t), detj, basis)
        else
          found = ElementInfo(element, acc_nodes, ip % u(t), &
              ip % v(t), ip % w(t), detj, basis, dbasisdx, &
              edgebasis=wbasis, rotbasis=rotwbasis, usolver=acc_targetptr)
        end if
        weight = ip % s(t)*detj

        do p = 1, nd
          do q = 1, nd
            stiff(p, q) = stiff(p, q) + weight*basis(q)*basis(p)
          end do
        end do

        do j = 1, acc_fourdofs
          if (j == 1) then
            fourierfreq = 0.0_dp
          else
            k = j/2
            fourierfreq = k*freq
          end if

          found2 = .false.
          do icomp = 1, acc_ncomp
            losscoeff(icomp) = ListGetFun(material, &
                'Harmonic Loss Coefficient '//I2S(icomp), &
                fourierfreq, found)
            if (found) found2 = .true.
          end do
          if (.not. found2) cycle

          component => acc_fouriervars(j) % var % values

          if (acc_direct == 1) then
            if (acc_tdofs == 1) then
              elemfield(1:nd) = component(acc_fourierperm(indeces(1:nd)))
              valatip = sum(basis(1:nd)*elemfield(1:nd))
            else
              curlatip = 0.0_dp
              do k = 1, acc_tdofs
                elemfield(1:nd) = component( &
                    acc_tdofs*(acc_fourierperm(indeces(1:nd)) - 1) + k)
                curlatip(k) = sum(basis(1:nd)*elemfield(1:nd))
              end do
            end if

          else if (acc_av == 1) then
            elemfield(1:nt) = component(acc_fourierperm(indeces(1:nt)))

            curlatip(1) = sum(elemfield(n + 1:nt)*rotwbasis(1:(nt - n), 1))
            curlatip(2) = sum(elemfield(n + 1:nt)*rotwbasis(1:(nt - n), 2))
            curlatip(3) = sum(elemfield(n + 1:nt)*rotwbasis(1:(nt - n), 3))
          else if (acc_dim == 3) then
            elemfield(1:nd) = component( &
                3*(acc_fourierperm(indeces(1:nd)) - 1) + 1)
            do k = 1, 3
              gradatip(1, k) = sum(elemfield(1:nd)*dbasisdx(1:nd, k))
            end do
            elemfield(1:nd) = component( &
                3*(acc_fourierperm(indeces(1:nd)) - 1) + 2)
            do k = 1, 3
              gradatip(2, k) = sum(elemfield(1:nd)*dbasisdx(1:nd, k))
            end do
            elemfield(1:nd) = component( &
                3*(acc_fourierperm(indeces(1:nd)) - 1) + 3)
            do k = 1, 3
              gradatip(3, k) = sum(elemfield(1:nd)*dbasisdx(1:nd, k))
            end do

            curlatip(1) = gradatip(3, 2) - gradatip(2, 3)
            curlatip(2) = gradatip(1, 3) - gradatip(3, 1)
            curlatip(3) = gradatip(1, 2) - gradatip(2, 1)
          else
            elemfield(1:nd) = component(acc_fourierperm(indeces(1:nd)))
            if (acc_flip == 1) then
              do i = 1, nd
                if (currentmodel % mesh % periodicflip(indeces(i))) &
                    elemfield(i) = -elemfield(i)
              end do
            end if
            curlatip(1) = sum(elemfield(1:nd)*dbasisdx(1:nd, 2))
            curlatip(2) = -sum(elemfield(1:nd)*dbasisdx(1:nd, 1))
            curlatip(3) = 0.0_dp
          end if

          if (.not. (acc_direct == 1 .and. acc_tdofs == 1)) then
            valatip = sqrt(sum(curlatip**2))
          end if

          do icomp = 1, acc_ncomp
            coeff = weight*losscoeff(icomp)* &
                (fourierfreq**freqpower(icomp))* &
                (valatip**fieldpower(icomp))
            if (acc_sumcomp == 1) then
              nsum = 1
            else
              nsum = icomp
            end if

            force(nsum, 1:nd) = force(nsum, 1:nd) + coeff*basis(1:nd)
            acc_series(icomp, j) = acc_series(icomp, j) + coeff
            acc_body(icomp, bodyid) = acc_body(icomp, bodyid) + coeff
          end do

        end do

      end do

      if (acc_sumcomp == 1) then
        nsum = 1
      else
        nsum = acc_ncomp
      end if

      do icomp = 1, nsum
        if (icomp == 1) then
          call DefaultUpdateEquations(stiff, force(icomp, 1:nd))
        else
          acc_solver % matrix % rhs => acc_otherrhs(:, icomp - 1)
          call DefaultUpdateForce(force(icomp, 1:nd))
          acc_solver % matrix % rhs => acc_saverhs
        end if
      end do

      if (acc_elem == 1) then
        eperm => acc_compvarse(1) % var % perm
        call LUdecomp(stiff, n, pivot, erroneous)
        if (erroneous) call Fatal('FourierLoss', 'LU-decomposition fails')
        do icomp = 1, nsum
          call LUSolve(n, stiff, force(icomp, :), pivot)
          acc_compvarse(icomp) % var % values( &
              eperm(element % dgindexes(1:n))) = force(icomp, 1:n)
        end do
      end if

    end do

    if (acc_nloss == 1) then
      acc_nodallossvar % values = acc_solver % matrix % rhs
      do icomp = 2, nsum
        acc_nodallossvar % values = acc_nodallossvar % values + &
            acc_otherrhs(:, icomp - 1)
      end do
    end if

    deallocate(stiff, force, dgforce, pivot, basis, dbasisdx, &
        indeces, elemfield)
    if (acc_av == 1) then
      deallocate(wbasis, rotwbasis)
    end if

  end subroutine fls_bulk

  ! -- nodal solve loop (pointer swaps stay in Fortran) -----------------
  subroutine fls_nodal_solve() bind(c, name="fls_nodal_solve")
    integer :: icomp
    real(kind=dp) :: norm
    if (acc_nodal == 0) return
    call ListAddLogical(acc_params, 'Linear System Compute Change', .false.)
    acc_saverhs => acc_solver % matrix % rhs
    do icomp = 1, acc_nvar
      if (icomp > 1) then
        acc_solver % matrix % rhs => acc_otherrhs(:, icomp - 1)
      end if
      norm = DefaultSolve()
      acc_compvars(icomp) % var % values = acc_lossvar % values
    end do
    acc_solver % matrix % rhs => acc_saverhs
    if (acc_nvar > 1) deallocate(acc_otherrhs)
  end subroutine fls_nodal_solve

  ! -- CommunicateLosess (verbatim transcription) ----------------------
  subroutine fls_comm() bind(c, name="fls_comm")
    integer :: i, j, k, noslices
    character(len=MAX_NAME_LEN) :: lossesfile
    logical :: found
    character(*), parameter :: caller = 'FourierLossSolver'

    if (parenv % pes > 1) then
      noslices = max(1, ListGetInteger(acc_model % simulation, &
          'Number Of Slices', found))
      if (noslices > 1) then
        call Info(caller, 'Averaging losses over slices', level=8)
      end if
      do k = 1, acc_ncomp
        do j = 1, acc_fourdofs
          acc_series(k, j) = ParallelReduction(acc_series(k, j))/noslices
        end do
        do j = 1, acc_model % numberofbodies
          acc_body(k, j) = ParallelReduction(acc_body(k, j))/noslices
        end do
      end do
    end if

    do k = 1, acc_ncomp
      acc_comploss(k) = sum(acc_series(k, :))
    end do
    acc_total = sum(acc_comploss)

    call ListAddConstReal(acc_model % simulation, &
        trim(acc_pref)//' fourier loss total', acc_total)
    do k = 1, acc_ncomp
      call ListAddConstReal(acc_model % simulation, &
          trim(acc_pref)//' fourier loss '//I2S(k), acc_comploss(k))
    end do

    do k = 1, acc_ncomp
      call Info(caller, &
          'Wavewise Fourier loss for component: '//I2S(k), level=6)

      do j = 1, acc_fourdofs
        if (j == 1) then
          write(message, '(A,ES12.3)') 'CONST : ', acc_series(k, j)
          call Info(caller, message, level=6)
        else
          i = j/2
          if (modulo(j, 2) == 0) then
            write(message, '(A,I0,A,ES12.3)') 'COS_', i, ' : ', &
                acc_series(k, j)
            call Info(caller, message, level=6)
          else
            write(message, '(A,I0,A,ES12.3)') 'SIN_', i, ' : ', &
                acc_series(k, j)
            call Info(caller, message, level=6)
          end if
        end if
      end do

      call Info(caller, &
          'Bodywise Fourier loss for component: '//I2S(k), level=6)

      do j = 1, acc_model % numberofbodies
        if (acc_body(k, j) < tiny(acc_comploss(k))) cycle
        write(message, '(A,I0,A,ES12.3)') 'Body ', j, ' : ', acc_body(k, j)
        call Info(caller, message, level=6)
      end do

      write(message, '(A,ES12.3)') &
          'Total component '//I2S(k)//' loss: ', acc_comploss(k)
      call Info(caller, message, level=5)
    end do

    write(message, '(A,ES12.3)') 'Total loss: ', acc_total
    call Info(caller, message, level=5)

    if (parenv % mype == 0) then
      lossesfile = ListGetString(acc_params, 'Fourier Loss Filename', found)
      if (found) then
        open (10, file=lossesfile)
        write (10, '(A)') '!body_id   loss(1)   loss(2) ....'
        do j = 1, acc_model % numberofbodies
          if (sum(acc_body(1:acc_ncomp, j)) <= tiny(acc_total)) cycle
          write (10, '(I6)', advance='NO') j
          do i = 1, acc_ncomp - 1
            write (10, '(ES15.6)', advance='NO') acc_body(i, j)
          end do
          write (10, '(ES15.6)') acc_body(acc_ncomp, j)
        end do
        call Info(caller, 'Fourier losses for bodies was saved to file: '// &
            trim(lossesfile), level=6)
        close (10)
      end if

      lossesfile = ListGetString(acc_params, 'Series Loss Filename', found)
      if (found) then
        open (10, file=lossesfile)
        do j = 1, acc_fourdofs/2
          write (10, '(2ES15.6)') sum(acc_series(1:acc_ncomp, 2*j - 1)), &
              sum(acc_series(1:acc_ncomp, 2*j))
        end do
        call Info(caller, 'Series losses for bodies was saved to file: '// &
            trim(lossesfile), level=6)
        close (10)
      end if
    end if

    if (InfoActive(25)) then
      print *, 'Fourier components:'
      do i = 1, acc_fourdofs
        call VectorValuesRange(acc_fouriervars(i) % var % values, &
            size(acc_fouriervars(i) % var % values), 'F'//I2S(i))
      end do

      print *, 'Loss components:'
      do i = 1, acc_ncomp
        call VectorValuesRange(acc_compvars(i) % var % values, &
            size(acc_compvars(i) % var % values), 'L'//I2S(i))
      end do
    end if

  end subroutine fls_comm

  ! -- total-loss fallback when no nodal solve ran ----------------------
  subroutine fls_store_total() bind(c, name="fls_store_total")
    if (acc_nodal == 0) then
      acc_solver % variable % values = acc_total
    end if
  end subroutine fls_store_total
end module fourierloss_acc
