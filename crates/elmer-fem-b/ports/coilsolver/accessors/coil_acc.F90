! T022-A CoilSolver port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/CoilSolver.F90, compiled against the built
! .mod files of the target configure (A2 shipped mechanism, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars). All Fortran-object traffic — TYPE(...)
! handles, OPTIONAL/POINTER dummies, assumed-shape descriptors, hidden
! lengths — stays inside Fortran, so the Rust side is pure C calls. The
! wrapper module mirrors the original's host-associated state (entry
! locals shared with the CONTAINS internals, current Set/PotVar
! pointers, per-part and per-coil selections) with identical lifetimes,
! and the whole-unit internal procedures are transcribed verbatim as
! module procedures. All floating-point work (assembly, SUM/MINVAL/
! MAXVAL reductions, DSYEV, WRITE(Message) formatting) stays in Fortran
! (same compiler), so bit-identity holds by construction wherever the
! call sequence matches.
!
! Serial-proof scope: the ParEnv % PEs > 1 branches (MPI_ALLREDUCE,
! ParallelReduction, ParallelSumVectorInt) are transcribed verbatim and
! guard on the same runtime condition, so the serial path — the only one
! the proof exercises — is bit-identical by construction. The parallel
! path is preserved for completeness but not proven.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module coil_acc
  use iso_c_binding
  use defutils
  implicit none

  character(*), parameter :: Caller = 'CoilSolver'

  ! -- entry state (mirrors CoilSolver/CoilSolver_init host association) --
  type(solver_t), pointer, save :: Solver => null()
  type(model_t), pointer, save :: Model => null()
  type(mesh_t), pointer, save :: Mesh => null()
  type(valuelist_t), pointer, save :: Params => null()
  type(valuelist_t), pointer, save :: CoilList => null()
  integer, pointer, save :: TargetBodies(:) => null()
  type(matrix_t), pointer, save :: StiffMatrix => null()
  real(kind=dp), pointer, save :: ForceVector(:) => null()
  integer, pointer, save :: Perm(:) => null()
  type(variable_t), pointer, save :: SolVar => null()
  type(variable_t), pointer, save :: PotVar => null()
  type(variable_t), pointer, save :: PotVarA => null()
  type(variable_t), pointer, save :: PotVarB => null()
  type(variable_t), pointer, save :: PotSelect => null()
  integer, allocatable, target, save :: SetA(:), SetB(:)
  type(variable_t), pointer, save :: FixVar => null()
  type(variable_t), pointer, save :: LoadVar => null()
  type(variable_t), pointer, save :: CoilIndexVar => null()
  type(variable_t), pointer, save :: CoilSetVar => null()
  type(variable_t), pointer, save :: FluxVar => null()
  type(variable_t), pointer, save :: FluxVarE => null()
  type(variable_t), pointer, save :: DistVar => null()
  integer, pointer, save :: Set(:) => null()
  integer, allocatable, save :: CoilIndex(:)
  logical, allocatable, save :: CoilBodies(:)
  real(kind=dp), allocatable, save :: DesiredCoilCurrent(:)
  real(kind=dp), allocatable, save :: DesiredCurrentDensity(:)
  real(kind=dp), allocatable, save :: CoilHelicity(:)
  real(kind=dp), allocatable, save :: CoilNormals(:,:)
  logical, allocatable, save :: GotCurr(:), GotDens(:), NormalizeCoil(:)
  real(kind=dp), save :: CoilCenter(3), CoilNormal(3)
  real(kind=dp), save :: CoilTangent1(3), CoilTangent2(3)
  real(kind=dp), save :: MinCurr(3), MaxCurr(3)
  character(len=MAX_NAME_LEN), save :: CondName, EqName
  integer, save :: nsize = 0, dim = 0, CoilParts = 1
  integer, save :: MaxNonlinIter = 1, NoCoils = 0, MaxNoCoils = 1
  integer, save :: iter = 1
  logical, save :: CoilClosed = .false., CoilAnisotropic = .false.
  logical, save :: UseDistance = .false., FixConductivity = .false.
  logical, save :: UseUnityCond = .false., FitCoil = .false.
  logical, save :: SelectNodes = .false., CalcCurr = .true.
  logical, save :: OneCut = .false., TestCut = .false.
  type(element_t), pointer, save :: acc_elem => null()

  ! -- per-internal-procedure SAVE Nodes mirrors --
  type(nodes_t), save :: acc_nodes_dcc, acc_nodes_dcp
  type(nodes_t), save :: acc_nodes_pot, acc_nodes_fix, acc_nodes_flux
  type(nodes_t), save :: acc_nodes_corr, acc_nodes_corr2, acc_nodes_norm

  ! -- CoilPotential SAVE mirrors --
  logical, save :: acc_cp_visited = .false.
  type(variable_t), pointer, save :: acc_cp_pota => null()
  type(variable_t), pointer, save :: acc_cp_potb => null()
  type(variable_t), pointer, save :: acc_cp_potp => null()
  type(variable_t), pointer, save :: acc_cp_potsel => null()
  type(element_t), pointer, save :: acc_cp_prevelem => null()

  ! -- CoilPotentialNormalized SAVE mirrors --
  logical, save :: acc_cn_visited = .false.
  type(variable_t), pointer, save :: acc_cn_pota => null()
  type(variable_t), pointer, save :: acc_cn_potb => null()
  type(variable_t), pointer, save :: acc_cn_potp => null()
  type(variable_t), pointer, save :: acc_cn_potsel => null()
  type(element_t), pointer, save :: acc_cn_prevelem => null()
  type(nodes_t), save :: acc_nodes_cpotn
  real(kind=dp), allocatable, save :: acc_cn_basis(:), acc_cn_dbasis(:,:)
  real(kind=dp), allocatable, save :: acc_cn_nodalpot(:)
  real(kind=dp), save :: acc_cn_xmin = 0.0_dp
  real(kind=dp), save :: acc_cn_desdens = 1.0_dp
  real(kind=dp), save :: acc_cn_gradpot(3) = 0.0_dp
  real(kind=dp), save :: acc_cn_normcoeff = 1.0_dp

contains

  ! -- CoilSolver_init ------------------------------------------------
  subroutine coa_init_begin(h_model, h_solver) bind(c, name="coa_init_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, Model)
    call c_f_pointer(h_solver, Solver)
  end subroutine coa_init_begin

  subroutine coa_init_body() bind(c, name="coa_init_body")
    type(valuelist_t), pointer :: Pars, BC
    integer :: d, i
    logical :: Found, CalcEl, CalcNod
    integer, pointer :: ElBCs(:)

    d = CoordinateSystemDimension()
    Pars => GetSolverParams()

    CalcEl = ListGetLogical( Pars,'Calculate Elemental Fields',Found )
    CalcNod = ListGetLogical( Pars,'Calculate Nodal Fields',Found )
    if (.not. Found ) CalcNod = .true.

    call ListAddNewString( Pars,'Variable','-nooutput CoilTmp')

    call ListAddString( Pars,NextFreeKeyword('Exported Variable',Pars),&
        'CoilPot')

    if ( GetLogical( Pars,'Coil Conductivity Fix', Found ) ) then
      call ListAddString( Pars,NextFreeKeyword('Exported Variable',Pars),&
          'CoilFix')
    end if

    if ( GetLogical( Pars,'Coil Closed', Found ) .or. &
        ListGetLogicalAnyComponent( Model,'Coil Closed') ) then
      call ListAddString( Pars,NextFreeKeyword('Exported Variable',Pars),&
          'CoilPotB')
      call ListAddString( Pars,NextFreeKeyword('Exported Variable',Pars),&
          'PotSelect')

      if ( GetLogical( Pars,'Save Coil Set', Found ) ) then
        call ListAddString( Pars,NextFreeKeyword('Exported Variable',Pars),&
            'CoilSetB')
      end if
    end if

    if ( GetLogical( Pars,'Save Coil Set', Found ) ) then
      call ListAddString( Pars,NextFreeKeyword('Exported Variable',Pars),&
          'CoilSet')
    end if

    if ( GetLogical( Pars,'Save Coil Index', Found ) ) then
      call ListAddString( Pars,NextFreeKeyword('Exported Variable',Pars),&
          '-elem CoilIndex e')
    end if

    if (.not. CalcEl ) then
      if ( ListGetLogical( Pars,'Fix Input Current Density',Found ) ) then
        CalcEl = .true.
        call Info('CoilSolver_int','Computing elemental field as only that can be fixed!')
      end if
    end if

    if ( CalcEl ) then
      call ListAddString( Pars,&
          NextFreeKeyword('Exported Variable',Pars),&
          '-dg -dofs '//I2S(d)//' CoilCurrent e')
    end if
    if ( CalcNod ) then
      call ListAddString( Pars,&
          NextFreeKeyword('Exported Variable',Pars),&
          '-dofs '//I2S(d)//' CoilCurrent')
    end if

    call ListAddLogical( Pars,'Calculate Loads',.true.)

    do i=1,Model % NumberOfComponents
      Pars => Model % Components(i) % Values

      ElBCs => ListGetIntegerArray(Pars, 'Electrode Boundaries', Found)
      if (.not. Found) cycle

      if (size(ElBCs) /= 2) then
        call Warn('CoilSolver_init','For CoilSolver the size of "Electrode Boundaries must be two!')
        cycle
      end if

      BC => CurrentModel % BCs(ElBCs(1)) % Values
      if (.not. associated(BC) ) call Fatal('CoilSolver_init', 'Boundary not found!')
      if ( ListGetLogical( BC,'Coil Start', Found ) ) then
        call Fatal('CoilSolver_init',&
            '"Electrode Boundaries(1)" not consistent with given "Coil Start" in bc '//I2S(ElBCs(1)))
      end if
      call ListAddLogical( BC,'Coil End',.true.)

      BC => CurrentModel % BCs(ElBCs(2)) % Values
      if (.not. associated(BC) ) call Fatal('CoilSolver_init', 'Boundary not found!')
      if ( ListGetLogical( BC,'Coil End', Found ) ) then
        call Fatal('CoilSolver_init',&
            '"Electrode Boundaries(2)" not consistent with given "Coil End" in bc '//I2S(ElBCs(2)))
      end if
      call ListAddNewLogical( BC,'Coil Start',.true.)

      call Info('CoilSolver_init','"Coil Start" and "Coil End" set to "Electrode Boundaries"',level=5)
    end do
  end subroutine coa_init_body

  ! -- CoilSolver head: guard, header, setup, flags --------------------
  function coa_matrix_associated(h_model, h_solver) result(r) &
      bind(c, name="coa_matrix_associated")
    type(c_ptr), value :: h_model, h_solver
    integer(c_int) :: r
    call c_f_pointer(h_model, Model)
    call c_f_pointer(h_solver, Solver)
    r = merge(1, 0, associated(Solver % matrix))
  end function coa_matrix_associated

  subroutine coa_header() bind(c, name="coa_header")
    call Info(Caller,'--------------------------------------')
    call Info(Caller,'Solving current distribution in a coil')
    call Info(Caller,'--------------------------------------')
  end subroutine coa_header

  subroutine coa_setup() bind(c, name="coa_setup")
    Params => GetSolverParams()
    nsize = size( Solver % Variable % Values )
    StiffMatrix => Solver % Matrix
    ForceVector => Solver % Matrix % Rhs
    Mesh => Solver % Mesh
    SolVar => Solver % Variable
    Perm => SolVar % Perm
    dim = CoordinateSystemDimension()
  end subroutine coa_setup

  subroutine coa_read_flags() bind(c, name="coa_read_flags")
    logical :: Found
    FixConductivity = GetLogical( Params,'Coil Conductivity Fix', Found )
    CoilAnisotropic = GetLogical( Params,'Coil Anisotropic', Found )
    CondName = GetString( Params,'Electric Conductivity Name',Found )
    if ( .not. Found ) CondName = 'Electric Conductivity'
    UseUnityCond = GetLogical( Params,'Use Unity Conductivity',Found )
    UseDistance = GetLogical( Params,'Use Wall Distance',Found)
    if ( UseDistance ) then
      DistVar => VariableGet( Mesh % Variables,'Wall Distance' )
      if ( .not. associated( DistVar ) ) then
        call Fatal(Caller,'> Wall Distance < not associated!')
      end if
    end if
    CoilClosed = GetLogical( Params,'Coil Closed', Found ) .or. &
        ListGetLogicalAnyComponent( Model,'Coil Closed')
    if ( CoilClosed ) then
      call Info(Caller,'Assuming that all coils are closed!',level=5)
      CoilParts = 2
    else
      call Info(Caller,'Assuming that all coils are open!',level=5)
      if ( .not. ListGetLogicalAnyBC( Model,'Coil Start') ) then
        call Fatal(Caller,'> Coil Start < must be defined on some BC')
      end if
      if ( .not. ListGetLogicalAnyBC( Model,'Coil End') ) then
        call Fatal(Caller,'> Coil End < must be defined on some BC')
      end if
      CoilParts = 1
    end if
    CalcCurr = GetLogical( Params,'Calculate Coil Current',Found )
    if ( .not. Found ) CalcCurr = .true.
    if ( FixConductivity ) then
      if ( CoilAnisotropic ) then
        MaxNonlinIter = 3
      else
        MaxNonlinIter = 2
      end if
    else
      MaxNonlinIter = 1
    end if
  end subroutine coa_read_flags

  function coa_parts() result(n) bind(c, name="coa_parts")
    integer(c_int) :: n
    n = CoilParts
  end function coa_parts

  function coa_calc_curr() result(r) bind(c, name="coa_calc_curr")
    integer(c_int) :: r
    r = merge(1, 0, CalcCurr)
  end function coa_calc_curr

  function coa_maxnonlin() result(n) bind(c, name="coa_maxnonlin")
    integer(c_int) :: n
    n = MaxNonlinIter
  end function coa_maxnonlin

  function coa_dim() result(n) bind(c, name="coa_dim")
    integer(c_int) :: n
    n = dim
  end function coa_dim

  function coa_num_components() result(n) bind(c, name="coa_num_components")
    integer(c_int) :: n
    n = Model % NumberOfComponents
  end function coa_num_components

  function coa_test_cut() result(r) bind(c, name="coa_test_cut")
    integer(c_int) :: r
    r = merge(1, 0, TestCut)
  end function coa_test_cut

  subroutine coa_pots_setup() bind(c, name="coa_pots_setup")
    logical :: Found
    PotVarA => VariableGet( Mesh % Variables,'CoilPot' )
    allocate( SetA(nsize) )
    SetA = 0
    if ( .not. associated( PotVarA ) ) then
      call Fatal(Caller,'CoilPot not associated!')
    end if
    allocate(CoilBodies(Model % NumberOfBodies))
    CoilBodies = .false.
    TestCut = .false.
    OneCut = .false.
    if ( CoilParts == 2 ) then
      OneCut = GetLogical( Params,'Single Coil Cut',Found )
      TestCut = GetLogical( Params,'Test Coil Cut',Found )
      allocate( SetB(nsize) )
      SetB = 0
      PotVarB => VariableGet( Mesh % Variables,'CoilPotB' )
      if ( .not. associated( PotVarB ) ) then
        call Fatal(Caller,'CoilPotB not associated!')
      end if
      PotSelect => VariableGet( Mesh % Variables,'PotSelect' )
      if ( .not. associated( PotSelect ) ) then
        call Fatal(Caller,'PotSelect not associated!')
      end if
    end if
    if ( FixConductivity ) then
      FixVar => VariableGet( Mesh % Variables,'CoilFix' )
      if ( .not. associated( FixVar ) ) then
        call Fatal(Caller,'CoilFix not associated!')
      end if
      FixVar % Values = 1.0_dp
    end if
    LoadVar => VariableGet( Mesh % Variables, trim(SolVar % Name)//' Loads' )
    if ( .not. associated( LoadVar ) ) then
      call Fatal(Caller,'> '//trim(SolVar % Name)//' < Loads not associated!')
    end if
    CoilIndexVar => VariableGet( Mesh % Variables,'CoilIndex e' )
    if ( associated( CoilIndexVar ) ) then
      call Info(Caller,'Saving elemental coil index to variable >CoilIndex e<')
    end if
    MaxNoCoils = max( 1, Model % NumberOfComponents )
    allocate( DesiredCoilCurrent( MaxNoCoils ), DesiredCurrentDensity(MaxNoCoils), &
        GotCurr( MaxNoCoils ), GotDens( MaxNoCoils), NormalizeCoil(MaxNoCoils), &
        CoilHelicity(MaxNoCoils),CoilNormals(MaxNoCoils,3))
    DesiredCoilCurrent = 0.0_dp
    GotCurr = .false.
    GotDens = .false.
    NormalizeCoil = .false.
    CoilHelicity = 0.0_dp
    CoilNormals = 0.0_dp
    allocate( CoilIndex(size(LoadVar % Perm)))
    CoilIndex = 0
  end subroutine coa_pots_setup
  ! -- coil loop: head (CYCLE=0, PROCEED=1, EXIT=2) and common tail -----
  function coa_coil_try(i) result(code) bind(c, name="coa_coil_try")
    integer(c_int), value :: i
    integer(c_int) :: code
    logical :: Found, ThisIsCoil
    integer :: t, Active, j
    type(element_t), pointer :: Element
    if ( i <= Model % NumberOfComponents ) then
      CoilList => Model % Components(i) % Values
      ThisIsCoil = ListCheckPrefix( CoilList,'Coil' )
      if (.not. ThisIsCoil) then
        code = 0
        return
      end if
      TargetBodies => ListGetIntegerArray( CoilList,'Master Bodies',Found )
      if ( .not. Found ) TargetBodies => ListGetIntegerArray( CoilList,'Body',Found )
      if ( .not. Found ) call Fatal(Caller,'Coil fitting requires > Master Bodies <')
      call Info(Caller,'Treating coil in Component: '//I2S(i),level=7)
      CoilBodies(TargetBodies) = .true.
      if (i==1) then
        if ( ListCheckPresent( Params,'Coil Normal') ) &
            call Fatal(Caller,'Place "Coil Normal" also in component section')
        if ( ListCheckPresent( Params,'Coil Center') ) &
            call Fatal(Caller,'Place "Coil Center" also in component section')
      end if
      NoCoils = NoCoils + 1
      call MarkCoilNodes(TargetBodies, NoCoils, Found )
      if (.not.Found) then
        NoCoils = NoCoils - 1
        code = 0
        return
      end if
      if ( CoilClosed ) then
        SelectNodes = .true.
        call DefineCoilCenter( CoilCenter, CoilList, TargetBodies )
        call DefineCoilParameters( CoilNormal, CoilTangent1, CoilTangent2, &
            CoilList, TargetBodies )
      end if
    else
      if ( NoCoils > 0 ) then
        Active = GetNOFActive()
        do t=1,Active
          Element => GetActiveElement(t)
          j = Element % BodyId
          if (.not. CoilBodies(j)) then
            call Fatal(Caller,'CoilSolver body '//I2S(j)//' active in Equation but not in Component!')
          end if
        end do
        deallocate(CoilBodies)
        code = 2
        return
      end if
      NoCoils = 1
      CoilList => Params
      if ( CoilClosed ) then
        call DefineCoilCenter( CoilCenter, Params )
        call DefineCoilParameters( CoilNormal, CoilTangent1, CoilTangent2, Params )
      end if
      CoilNormals(NoCoils,:) = CoilNormal
    end if
    call Info(Caller,'Coil system consists of '//I2S(NoCoils)//' coils',level=7)
    code = 1
  end function coa_coil_try

  subroutine coa_coil_tail() bind(c, name="coa_coil_tail")
    logical :: Found
    real(kind=dp) :: CoilCrossSection
    if ( CoilClosed ) then
      Set => SetA
      call ChooseFixedBulkNodesNarrow(Set,1,SelectNodes)
      if ( OneCut ) call ChooseCoilCut(Set,SelectNodes)
      Set => SetB
      call ChooseFixedBulkNodesNarrow(Set,2,SelectNodes)
      if ( OneCut ) call ChooseCoilCut(Set,SelectNodes)
    else
      Set => SetA
      call ChooseFixedEndNodes(Set)
    end if
    DesiredCoilCurrent(NoCoils) = ListGetCReal( CoilList,'Desired Coil Current',Found )
    if (.not. Found ) DesiredCoilCurrent(NoCoils) = 1.0_dp
    GotCurr(NoCoils) = Found
    DesiredCurrentDensity(NoCoils) = ListGetCReal( CoilList,'Desired Current Density',Found)
    if (.not. Found ) DesiredCurrentDensity(NoCoils) = 1.0_dp
    GotDens(NoCoils) = Found
    NormalizeCoil(NoCoils) = ListGetLogical( CoilList,'Normalize Coil Current',Found )
    if (.not. Found ) then
      if ( ListGetLogical( Params,'Normalize Coil Current',Found ) ) NormalizeCoil(NoCoils) = .true.
    end if
    CoilCrossSection = ListGetCReal( CoilList,'Coil Cross Section',Found )
    if ( Found ) then
      if ( GotCurr(NoCoils) ) then
        DesiredCurrentDensity(NoCoils) = DesiredCoilCurrent(NoCoils) / CoilCrossSection
        GotDens(NoCoils) = .true.
      else if ( GotDens(NoCoils) ) then
        DesiredCoilCurrent(NoCoils) = DesiredCurrentDensity(NoCoils) * CoilCrossSection
        GotCurr(NoCoils) = .true.
      end if
    end if
    CoilHelicity(NoCoils) = ListGetCReal( CoilList,'Coil Helicity',Found )
  end subroutine coa_coil_tail

  subroutine coa_count_fixing() bind(c, name="coa_count_fixing")
    Set => SetA
    call CountFixingNodes(Set,1)
    if ( CoilClosed ) then
      call Info(Caller,'Defining 2nd set of fixing nodes',level=20)
      Set => SetB
      call CountFixingNodes(Set,2)
    end if
  end subroutine coa_count_fixing

  ! -- nonlinear iteration / potential stages --------------------------
  subroutine coa_set_iter(it) bind(c, name="coa_set_iter")
    integer(c_int), value :: it
    iter = it
  end subroutine coa_set_iter

  subroutine coa_fix_stage() bind(c, name="coa_fix_stage")
    integer :: t, Active, n, nd
    real(kind=dp) :: Norm
    call Info(Caller,'Fixing the conductivity field')
    call DefaultInitialize()
    Active = GetNOFActive()
    do t=1,Active
      acc_elem => GetActiveElement(t)
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      call LocalFixMatrix( acc_elem, n, nd )
    end do
    call ListAddLogical( Params,'Skip Compute Nonlinear Change',.true.)
    Norm = DefaultSolve()
    FixVar % Values = SolVar % Values
  end subroutine coa_fix_stage

  subroutine coa_pot_header() bind(c, name="coa_pot_header")
    call Info(Caller,'Computing the dummy potential field')
  end subroutine coa_pot_header

  subroutine coa_pot_part(part) bind(c, name="coa_pot_part")
    integer(c_int), value :: part
    integer :: t, Active, n, nd, i, sgn
    real(kind=dp) :: Norm, val
    if ( part == 1 ) then
      Set => SetA
      PotVar => PotVarA
    else
      Set => SetB
      PotVar => PotVarB
    end if
    call DefaultInitialize()
    Active = GetNOFActive()
    do t=1,Active
      acc_elem => GetActiveElement(t)
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      call LocalPotMatrix( acc_elem, n, nd )
    end do
    call DefaultFinishBulkAssembly()
    call DefaultFinishAssembly()
    do i = 1,nsize
      sgn = Set(i)
      if ( sgn == 0 ) cycle
      if ( sgn > 0 ) then
        val = 1.0_dp
      else
        val = 0.0_dp
      end if
      call UpdateDirichletDof( Solver % Matrix, i, val )
    end do
    call CutInterfaceConnections( StiffMatrix, Set )
    call DefaultDirichletBCs()
    call ListAddLogical( Params,'Skip Compute Nonlinear Change',.false.)
    Norm = DefaultSolve()
    PotVar % Values = SolVar % Values
    call ScalePotential()
  end subroutine coa_pot_part

  ! -- current stages ---------------------------------------------------
  subroutine coa_curr_begin() bind(c, name="coa_curr_begin")
    call Info(Caller,'Computing current densities',level=20)
    call ListAddLogical( Params,'Calculate Loads',.false.)
    MinCurr = 0.0_dp
    MaxCurr = 0.0_dp
  end subroutine coa_curr_begin

  subroutine coa_flux_comp(dimi) bind(c, name="coa_flux_comp")
    integer(c_int), value :: dimi
    integer :: t, Active, n, nd
    real(kind=dp) :: Norm
    call Info(Caller,'Computing current component: '//I2S(dimi),level=6)
    FluxVar => VariableGet( Mesh % Variables,'CoilCurrent '//I2S(dimi) )
    FluxVarE => VariableGet( Mesh % Variables,'CoilCurrent e '//I2S(dimi) )
    if ( .not. ( associated( FluxVarE ) .or. associated( FluxVar ) ) ) then
      call Fatal(Caller,'CoilCurrent not associated!')
    end if
    call DefaultInitialize()
    Active = GetNOFActive()
    do t=1,Active
      acc_elem => GetActiveElement(t)
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      call LocalFluxMatrix( acc_elem, n, nd, dimi )
    end do
    if (associated( FluxVar ) ) then
      call ListAddLogical( Params,'Skip Compute Nonlinear Change',.true.)
      Norm = DefaultSolve()
      MinCurr(dimi) = minval( SolVar % Values )
      MaxCurr(dimi) = maxval( SolVar % Values )
      FluxVar % Values = SolVar % Values
    else
      MinCurr(dimi) = minval( FluxVarE % Values(dimi::dim) )
      MaxCurr(dimi) = maxval( FluxVarE % Values(dimi::dim) )
    end if
  end subroutine coa_flux_comp

  subroutine coa_curr_part_tail() bind(c, name="coa_curr_part_tail")
    integer :: ierr
    real(kind=dp) :: TmpCurr(3)
    if ( ParEnv % PEs > 1 ) then
      TmpCurr = MinCurr
      call MPI_ALLREDUCE(MinCurr,TmpCurr,3,MPI_DOUBLE_PRECISION,MPI_MIN,ParEnv % ActiveComm,ierr)
      MinCurr = TmpCurr
      TmpCurr = MaxCurr
      call MPI_ALLREDUCE(MaxCurr,TmpCurr,3,MPI_DOUBLE_PRECISION,MPI_MAX,ParEnv % ActiveComm,ierr)
      MaxCurr = TmpCurr
    end if
    write( Message,'(A,3ES12.4)') 'Minimum current components: ',MinCurr
    call Info(Caller,Message,level=7)
    write( Message,'(A,3ES12.4)') 'Maximum current components: ',MaxCurr
    call Info(Caller,Message,level=7)
  end subroutine coa_curr_part_tail

  subroutine coa_curr_end() bind(c, name="coa_curr_end")
    call NormalizeCurrentDensity()
    call ListAddLogical( Params,'Calculate Loads',.true.)
  end subroutine coa_curr_end

  function coa_fixinput() result(r) bind(c, name="coa_fixinput")
    integer(c_int) :: r
    logical :: Found
    r = merge(1, 0, ListGetLogical( Params,'Fix Input Current Density',Found ))
  end function coa_fixinput

  subroutine coa_fixinput_stage() bind(c, name="coa_fixinput_stage")
    integer :: t, Active, n, nd, i
    real(kind=dp) :: Norm
    type(valuelist_t), pointer :: BC
    logical :: Found
    call Info(Caller,'Finding fixing potential for elemental current density!')
    FluxVarE => VariableGet( Mesh % Variables,'CoilCurrent E')
    if ( .not. ( associated( FluxVarE ) ) ) then
      call Fatal(Caller,'Fixing can be done only for elemental coil current!')
    end if
    call DefaultInitialize()
    Active = GetNOFActive()
    do t=1,Active
      acc_elem => GetActiveElement(t)
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      call LocalCorrMatrix( acc_elem, n, nd )
    end do
    if ( .not. CoilClosed ) then
      do i=1,Model % NumberOfBCs
        BC => Model % BCs(i) % Values
        if ( ListGetLogical( BC,'Coil Start', Found ) .or. ListGetLogical( BC,'Coil End', Found ) ) then
          call ListAddConstReal( BC,Solver % Variable % Name, 0.0_dp )
        end if
      end do
    else
      i = ParEnv % MyPe
      i = ParallelReduction(i,1)
      if (i==ParEnv % MyPe) then
        call UpdateDirichletDof(Solver % Matrix, 1, 0._dp)
      end if
    end if
    call DefaultDirichletBCs()
    Solver % Variable % Values = 0
    Norm = DefaultSolve()
    if ( .not. CoilClosed ) then
      do i=1,Model % NumberOfBCs
        BC => Model % BCs(i) % Values
        if ( ListGetLogical( BC,'Coil Start', Found ) .or. ListGetLogical( BC,'Coil End', Found ) ) then
          call ListRemove( BC,Solver % Variable % Name )
        end if
      end do
    end if
    call Info(Caller,'Fixing elemental current density to be divergence free!')
    do t=1,Active
      acc_elem => GetActiveElement(t)
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      call LocalCorrCurrent( acc_elem, n, nd )
    end do
  end subroutine coa_fixinput_stage

  ! -- label-100 tail ----------------------------------------------------
  subroutine coa_tail_100() bind(c, name="coa_tail_100")
    integer :: i, j
    CoilSetVar => VariableGet( Mesh % Variables,'CoilSet' )
    if ( associated( CoilSetVar ) ) then
      do i=1,Mesh % NumberOfNodes
        j = CoilSetVar % Perm(i)
        if ( j > 0 ) CoilSetVar % Values(j) = 1.0_dp * SetA(j)
      end do
    end if
    if ( CoilParts == 2 ) then
      CoilSetVar => VariableGet( Mesh % Variables,'CoilSetB' )
      if ( associated( CoilSetVar ) ) then
        do i=1,Mesh % NumberOfNodes
          j = CoilSetVar % Perm(i)
          if ( j > 0 ) CoilSetVar % Values(j) = 1.0_dp * SetB(j)
        end do
      end if
    end if
    Solver % Variable % Values = PotVarA % Values
    if ( InfoActive(30) ) then
      call Info(Caller,'Final current densities')
      if ( associated( FluxVar ) ) then
        call VectorValuesRange(FluxVar % Values,size(FluxVar % Values),'Current density')
      end if
      if ( associated( FluxVarE ) ) then
        call VectorValuesRange(FluxVarE % Values,size(FluxVarE % Values),'Current density e')
      end if
    end if
    call Info(Caller,'All done',level=7)
    call Info(Caller,'--------------------------------------')
  end subroutine coa_tail_100
  ! -- transcribed CONTAINS internals (verbatim bodies) ------------------
  subroutine CutInterfaceConnections( A, Set )
    type(Matrix_t) :: A
    integer, pointer :: Set(:)
    integer, pointer, contiguous :: Cols(:),Rows(:)
    real(kind=dp), pointer, contiguous :: Values(:)
    integer :: i,j,k,n,cnt,s1,s2,dia
    n = A % NumberOfRows
    Rows   => A % Rows
    Cols   => A % Cols
    if (.not. associated( A % BulkValues ) ) then
      call Fatal('CutInterfaceConnections','Do not have bulk values!')
    end if
    Values => A % BulkValues
    cnt = 0
    do i=1,n
      s1 = Set(i)
      if ( s1 == 0 ) cycle
      dia = A % Diag(i)
      do j=Rows(i),Rows(i+1)-1
        k = Cols(j)
        s2 = Set(k)
        if ( s1 * s2 < 0 ) then
          A % BulkValues(dia) = A % BulkValues(dia) + A % BulkValues(j)
          A % BulkValues(j) = 0.0_dp
          cnt = cnt + 1
        end if
      end do
    end do
    cnt = ParallelReduction(cnt)
    call Info('CutInterfaceConnections','Number of connections cut: '//I2S(cnt),level=7)
  end subroutine CutInterfaceConnections

  subroutine MarkCoilNodes(TargetBodies, Coil, CoilActive )
    integer :: TargetBodies(:)
    integer :: Coil
    logical :: CoilActive
    integer :: Active, e, n, i, j, k
    logical :: Found
    type(element_t), pointer :: Element
    CoilActive = .true.
    EqName = ListGetString( Params, 'Equation', Found)
    do i=1,size(TargetBodies)
      j = TargetBodies(i)
      k = ListGetInteger( Model % Bodies(j) % Values, 'Equation', Found )
      if (.not. ListGetLogical(Model % Equations(k) % Values,EqName,Found)) then
        CoilActive = .false.
        call Warn(Caller,'CoilSolver body '//I2S(j)//' active in Component but not in Equation!')
      end if
    end do
    if (.not. CoilActive) return
    Active = GetNOFActive()
    do e=1,Active
      Element => GetActiveElement(e)
      n  = GetElementNOFNodes()
      if ( all( TargetBodies /= Element % BodyId ) ) cycle
      CoilIndex( Element % NodeIndexes ) = NoCoils
      if ( associated( CoilIndexVar ) ) then
        n = CoilIndexVar % Perm( Element % ElementIndex )
        if ( n > 0 ) CoilIndexVar % Values(n) = 1.0_dp * NoCoils
      end if
    end do
  end subroutine MarkCoilNodes

  subroutine DefineCoilCenter(CoilCenter, Pars, TargetBodies)
    real(kind=dp) :: CoilCenter(3)
    type(valuelist_t), pointer :: Pars
    integer, optional :: TargetBodies(:)
    real(kind=dp), allocatable :: Basis(:)
    real(kind=dp) :: DetJ,r(3),s
    integer :: e,t,i,j,n,Active,ierr
    logical :: stat,Found,CoilCenterSet
    type(element_t), pointer :: Element
    type(GaussIntegrationPoints_t) :: IP
    real(kind=dp) :: Volume,Center(3),SerTmp(4),ParTmp(4)
    real(kind=dp), pointer :: HelperArray(:,:)
    real(kind=dp) :: rArray(3,1)
    n = Mesh % MaxElementNodes
    allocate( Basis(n) )
    HelperArray => ListGetConstRealArray( Pars, 'Coil Center', CoilCenterSet)
    if ( CoilCenterSet ) then
      Center(1:3) = HelperArray(1:3,1)
    else
      Center(1) = ListGetCReal( Pars,'Coil x0',CoilCenterSet)
      Center(2) = ListGetCReal( Pars,'Coil y0',Found)
      CoilCenterSet = CoilCenterSet .or. Found
      Center(3) = ListGetCReal( Pars,'Coil z0',Found)
      CoilCenterSet = CoilCenterSet .or. Found
    end if
    if ( CoilCenterSet ) then
      CoilCenter = Center
      call Info(Caller,'Coil center defined by user',level=20)
      return
    end if
    Volume = 0.0_dp
    Active = GetNOFActive()
    do e=1,Active
      Element => GetActiveElement(e)
      n  = GetElementNOFNodes()
      if ( present( TargetBodies ) ) then
        if ( all( TargetBodies /= Element % BodyId ) ) cycle
      end if
      call GetElementNodes( acc_nodes_dcc, Element )
      IP = GaussPoints(Element)
      do t=1,IP % n
        stat = ElementInfo( Element, acc_nodes_dcc, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis )
        r(1) = sum(acc_nodes_dcc % x(1:n)  *Basis(1:n))
        r(2) = sum(acc_nodes_dcc % y(1:n) *Basis(1:n))
        r(3) = sum(acc_nodes_dcc % z(1:n) *Basis(1:n))
        s = IP % s(t) * detJ
        Volume = Volume + s
        Center = Center + s * r
      end do
    end do
    if ( ParEnv % PEs > 1 ) then
      SerTmp(1:3) = Center
      SerTmp(4) = Volume
      call MPI_ALLREDUCE(SerTmp,ParTmp,4,MPI_DOUBLE_PRECISION,MPI_SUM,ParEnv % ActiveComm,ierr)
      Center = ParTmp(1:3)
      Volume = ParTmp(4)
    end if
    if ( Volume < epsilon( Volume ) ) call Fatal('DefineCoilCenter','Coil has no volume!')
    CoilCenter = Center / Volume
    write( Message,'(A,ES12.4)') 'Coil volume:',Volume
    call Info(Caller,Message,level=7)
    write( Message,'(A,3ES12.4)') 'Coil center:',CoilCenter
    call Info(Caller,Message,level=7)
    rArray(1:3,1) = CoilCenter
    call ListAddConstRealArray( Pars,'Coil Center',3,1,rArray)
  end subroutine DefineCoilCenter

  subroutine DefineCoilParameters(CoilNormal, CoilTangent1, CoilTangent2, &
      Pars, TargetBodies )
    real(kind=dp) :: CoilNormal(3), CoilTangent1(3), CoilTangent2(3)
    type(valuelist_t), pointer :: Pars
    integer, optional :: TargetBodies(:)
    real(kind=dp), pointer :: HelperArray(:,:)
    real(kind=dp), allocatable :: Basis(:)
    real(kind=dp) :: DetJ,r(3),s,CoilTangentTmp(3)
    integer :: e,t,i,j,n,Active,ierr
    logical :: stat,Found,Found2,TallCoil
    type(element_t), pointer :: Element
    type(GaussIntegrationPoints_t) :: IP
    real(kind=dp) :: Imoment(9), EigVec(3,3), EigVal(3), ParTmp(9), CP(3)
    real(kind=dp) :: EigWrk(20)
    integer :: EigInfo, Three
    real(kind=dp) :: rArray(1:3,1)
    FitCoil = GetLogical( Pars,'Fit Coil',Found )
    if (.not. Found ) FitCoil = .true.
    HelperArray => ListGetConstRealArray( Pars, 'Coil normal', Found)
    if (.not. FitCoil .or. Found) then
      Found2 = .false.
      if (Found) then
        CoilNormal(1:3) = HelperArray(1:3,1)
        HelperArray => ListGetConstRealArray( Pars, 'Coil Tangent', Found2)
        if ( Found2 ) then
          CoilTangent1(1:3) = HelperArray(1:3,1)
          CoilTangent2 = CrossProduct(CoilNormal,CoilTangent1)
          CoilTangent2 = CoilTangent2 / sqrt(sum(CoilTangent2**2))
        end if
      else
        CoilNormal = 0.0_dp
        CoilNormal(3) = 1.0_dp
      end if
      if (.not.Found2) then
        call TangentDirections(CoilNormal, CoilTangent1, CoilTangent2)
      end if
      return
    end if
    call Info('DefineCoilParameters','Fitting the coil by maximizing inertia',level=7)
    n = Mesh % MaxElementNodes
    allocate( Basis(n) )
    Imoment = 0.0_dp
    Active = GetNOFActive()
    do e=1,Active
      Element => GetActiveElement(e)
      if ( present( TargetBodies ) ) then
        if ( all( TargetBodies /= Element % BodyId ) ) cycle
      end if
      n  = GetElementNOFNodes()
      call GetElementNodes( acc_nodes_dcp, Element )
      IP = GaussPoints(Element)
      do t=1,IP % n
        stat = ElementInfo( Element, acc_nodes_dcp, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis )
        r(1) = sum(acc_nodes_dcp % x(1:n)  *Basis(1:n))
        r(2) = sum(acc_nodes_dcp % y(1:n) *Basis(1:n))
        r(3) = sum(acc_nodes_dcp % z(1:n) *Basis(1:n))
        s = IP % s(t) * detJ
        r = r - CoilCenter
        do i=1,3
          Imoment(3*(i-1)+i) = Imoment(3*(i-1)+i) + s * sum( r**2 )
          do j=1,3
            Imoment(3*(i-1)+j) = Imoment(3*(i-1)+j) - s * r(i) * r(j)
          end do
        end do
      end do
    end do
    if ( ParEnv % PEs > 1 ) then
      call MPI_ALLREDUCE(Imoment,ParTmp,9,MPI_DOUBLE_PRECISION,MPI_SUM,ParEnv % ActiveComm,ierr)
      Imoment = ParTmp
    end if
    s = sum(Imoment)
    Imoment = Imoment / s
    s = 1.0_dp
    do i=1,3
      do j=1,3
        EigVec(i,j) = Imoment(3*(i-1)+j)
      end do
      EigVec(i,i) = EigVec(i,i) - s
    end do
    EigInfo = 0
    Three = 3
    call DSYEV( 'V','U', Three, EigVec, Three, EigVal, EigWrk, size(EigWrk), EigInfo )
    if (EigInfo /= 0) then
      call Fatal( Caller, 'DSYEV cannot generate eigen basis')
    end if
    write( Message,'(A,3ES12.4)') 'Coil inertia eigenvalues:',EigVal
    call Info(Caller,Message,level=10)
    TallCoil = ListGetLogical( Pars, 'Coil Geometry Tall',Found)
    if (.not.Found) then
      TallCoil = ( abs(EigVal(2)-EigVal(3)) < 0.01 * abs(EigVal(2)-EigVal(1)) )
      if ( TallCoil ) then
      end if
    end if
    if ( TallCoil ) then
      call Info( Caller,'Assuming coil to be tall i.e. axis has minimum inertia')
      CoilNormal = EigVec(:,1)
    else
      call Info( Caller,'Assuming coil to be flat i.e. axis has maximum inertia')
      CoilNormal = EigVec(:,3)
    end if
    CP = CrossProduct( CoilCenter, CoilNormal )
    j = 1
    do i = 2, 3
      if ( abs( CP(i) ) > abs( CP(j) ) ) j = i
    end do
    if ( CP(j) < 0 ) then
      call Info(Caller,'Inverting sign of coil normal',level=10)
      CoilNormal = -CoilNormal
    end if
    call TangentDirections(CoilNormal, CoilTangent1, CoilTangent2)
    write( Message,'(A,3ES12.4)') 'Coil axis normal:',CoilNormal
    call Info(Caller,Message,level=7)
    write( Message,'(A,3ES12.4)') 'Coil tangent1:',CoilTangent1
    call Info(Caller,Message,level=10)
    write( Message,'(A,3ES12.4)') 'Coil tangent2:',CoilTangent2
    call Info(Caller,Message,level=10)
    rArray(1:3,1) = CoilNormal
    call ListAddConstRealArray( Pars,'Coil Normal',3,1,rArray)
  end subroutine DefineCoilParameters

  subroutine ChooseFixedBulkNodesNarrow( Set, SetNo, SelectNodes )
    use ParallelUtils, only: ParallelSumVectorInt
    integer :: SetNo
    integer :: Set(:)
    logical :: SelectNodes
    logical :: Mirror
    type(mesh_t), pointer :: LMesh
    real(kind=dp) :: x,y,z,x0,y0
    real(kind=dp) :: r(3),rp(3),MinCut, MaxCut, CutDist(27)
    integer :: t,i,j,k,l,n,ioffset
    logical :: Found, WrongSide
    type(element_t), pointer :: Element
    integer, pointer :: Indexes(:)
    call Info(Caller,'Choosing fixing nodes for set: '//I2S(SetNo))
    Mirror = ( SetNo == 2 )
    LMesh => Solver % Mesh
    ioffset = 10 * NoCoils
    do t=1,LMesh % NumberOfBulkElements
      Element => LMesh % Elements(t)
      Indexes => Element % NodeIndexes
      n = Element % Type % NumberOfNodes
      if ( SelectNodes ) then
        if ( any( CoilIndex(Indexes) /= NoCoils ) ) cycle
      end if
      WrongSide = .false.
      l = 0
      do k = 1, n
        i = Indexes(k)
        j = Perm(i)
        if ( j == 0 ) cycle
        r(1) = LMesh % Nodes % x(i)
        r(2) = LMesh % Nodes % y(i)
        r(3) = LMesh % Nodes % z(i)
        r = r - CoilCenter
        if ( mirror ) r = -r
        rp(1) = sum( CoilTangent1 * r )
        rp(2) = sum( CoilTangent2 * r )
        rp(3) = sum( CoilNormal * r )
        if ( SetNo == 1 ) then
          PotSelect % Values( PotSelect % Perm(i) ) = rp(1)
        end if
        if ( rp(1) < 0 ) then
          WrongSide = .true.
          exit
        end if
        l = l + 1
        CutDist(l) = rp(2)
      end do
      if (WrongSide) cycle
      if (l<2) cycle
      MaxCut = maxval( CutDist(1:l) )
      MinCut = minval( CutDist(1:l) )
      if ( MaxCut >= -epsilon( MaxCut) .and. MinCut <= epsilon( MinCut)  ) then
        l = 0
        do k=1,n
          i = Indexes(k)
          j = Perm(i)
          if (j==0) cycle
          l = l + 1
          if ( CutDist(l) > 0.0_dp ) then
            Set(j) = 2 + ioffset
          else
            Set(j) = -2 - ioffset
          end if
        end do
      end if
    end do
    i = count(Set<0)
    j = count(Set>0)
    if ( ParEnv % PEs > 1 ) then
      block
        integer, allocatable :: MinSet(:)
        allocate(MinSet(size(Set)))
        MinSet = Set
        call ParallelSumVectorInt( Solver % Matrix, Set, 2 )
        call ParallelSumVectorInt( Solver % Matrix, MinSet, 1)
        where (abs(MinSet) > Set)
          Set = MinSet
        end where
        deallocate(MinSet)
      end block
    end if
  end subroutine ChooseFixedBulkNodesNarrow

  subroutine ChooseCoilCut(Set, SelectNodes )
    use ParallelUtils, only: ParallelSumVectorInt
    integer :: Set(:)
    logical :: SelectNodes
    logical :: Ready, GotAny, Parallel
    logical :: Found
    integer :: i,j,k,l,n,t,MinIndex,MaxIndex,Loop,ParLoop,NoPieces,jmax,m,NoCand
    integer, allocatable :: MeshPiece(:)
    type(element_t), pointer :: Element
    integer, pointer :: Indexes(:)
    integer :: pIndexes(20)
    Parallel = ( ParEnv % PEs > 1 )
    m = Solver % Matrix % NumberOfRows
    allocate( MeshPiece(m) )
    MeshPiece = 0
    do t = 1, Mesh % NumberOfBulkElements
      Element => Mesh % Elements(t)
      Indexes => Element % NodeIndexes
      n = Element % TYPE % NumberOfNodes
      pIndexes(1:n) = Perm(Indexes)
      if ( any( pIndexes(1:n) == 0 ) ) cycle
      if ( SelectNodes ) then
        if ( any( CoilIndex(Indexes) /= NoCoils ) ) cycle
      end if
      do i=1,n
        if ( Parallel ) then
          if ( Set(pIndexes(i)) /= 0 ) MeshPiece(pIndexes(i)) = Mesh % Parallelinfo % GlobalDOFs(Indexes(i))
        else
          if ( Set(pIndexes(i)) /= 0 ) MeshPiece(pIndexes(i)) = Indexes(i)
        end if
      end do
    end do
    NoCand = count(MeshPiece > 0)
    call Info(Caller,&
        'Number of candidante nodes in coil '//I2S(NoCand),level=12)
    if ( Parallel ) then
      jmax = ParallelReduction(NoCand,2)
      call Info(Caller,&
          'Total number of candidante nodes in coil '//I2S(jmax),level=12)
    end if
    Ready = (NoCand == 0)
    Loop = 0
    ParLoop = 0
100  do while(.not. Ready)
      Ready = .true.
      do t = 1, Mesh % NumberOfBulkElements
        Element => Mesh % Elements(t)
        Indexes => Element % NodeIndexes
        n = Element % TYPE % NumberOfNodes
        pIndexes(1:n) = Perm(Indexes)
        GotAny = .false.
        do i=1,n
          j = Perm(Indexes(i))
          if (j==0) cycle
          k = MeshPiece(j)
          if (k==0) cycle
          if (SelectNodes) then
          end if
          if (GotAny) then
            MinIndex = min(MinIndex,k)
            MaxIndex = max(MaxIndex,k)
          else
            MinIndex = k
            MaxIndex = k
            GotAny = .true.
          end if
        end do
        if (.not. GotAny) cycle
        if (MinIndex == MaxIndex) cycle
        Ready = .false.
        do i=1,n
          j = Perm(Indexes(i))
          if (j==0) cycle
          k = MeshPiece(j)
          if (k==0) cycle
          if (SelectNodes) then
          end if
          MeshPiece(j) = MaxIndex
        end do
      end do
      Loop = Loop + 1
    end do
    call Info(Caller,'Coil cut coloring loops: '//I2S(Loop),level=12)
    if ( NoCand > 0 ) then
      MaxIndex = maxval(MeshPiece)
    else
      MaxIndex = 0
    end if
    MinIndex = MaxIndex
    if ( Parallel ) then
      i = sum(MeshPiece)
      call ParallelSumVectorInt(Solver % Matrix,MeshPiece,2)
      j = sum(MeshPiece)-i
      k = ParallelReduction(j)
      if (k > 0 ) then
        ParLoop = ParLoop + 1
        call Info(Caller,'Continuing after parallel reduction: '//I2S(k),level=6)
        Ready = .false.
        goto 100
      end if
    end if
    if ( NoCand > 0 ) then
      MaxIndex = maxval(MeshPiece)
    else
      MaxIndex = 0
    end if
    MinIndex = MaxIndex
    NoPieces = 0
    l = 0
    if ( MaxIndex > 0 ) then
      NoPieces = 1
      k = MaxIndex
200   j = 0
      do i=1,m
        if (MeshPiece(i)>j .and. MeshPiece(i)<k) then
          j=MeshPiece(i)
        end if
      end do
      if ( j>0 ) then
        NoPieces = NoPieces + 1
        k = j
        MinIndex = j
        goto 200
      end if
    end if
    NoPieces = ParallelReduction(NoPieces)
    call Info(Caller,'Number of separate cuts in mesh is '//I2S(NoPieces),level=12)
    if (NoPieces == 1 ) return
    MaxIndex = ParallelReduction(MaxIndex,2)
    if ( ListGetLogical( Solver % Values,'Select Min Coil Cut',Found ) ) then
      if ( NoCand > 0 ) then
        MinIndex = minval(MeshPiece,MeshPiece>0)
      else
        MinIndex = MaxIndex
      end if
      MinIndex = ParallelReduction(MinIndex,1)
      where ( MeshPiece /= MinIndex ) Set = 0
    else
      where ( MeshPiece /= MaxIndex) Set = 0
    end if
  end subroutine ChooseCoilCut

  subroutine ChooseFixedEndNodes( Set )
    integer :: Set(:)
    type(mesh_t), pointer :: LMesh
    integer, pointer :: Indexes(:)
    integer :: t,i,j,k,nminus,nplus,ioffset
    type(valuelist_t), pointer :: BC
    logical :: Found
    Set = 0
    ioffset = 10 * NoCoils
    LMesh => Solver % Mesh
    do t=LMesh % NumberOfBulkElements + 1, &
         LMesh % NumberOfBulkElements + LMesh % NumberOfBoundaryElements
      acc_elem => LMesh % Elements(t)
      Model % CurrentElement => acc_elem
      Indexes => acc_elem % NodeIndexes
      if ( any( Perm( Indexes ) == 0 ) ) cycle
      BC => GetBC(acc_elem)
      if ( ListGetLogical( BC,'Coil Start', Found ) ) then
        Set( Perm( Indexes ) ) = 2 + ioffset
      else if ( ListGetLogical( BC,'Coil End', Found ) ) then
        Set( Perm( Indexes ) ) = -2 - ioffset
      end if
    end do
  end subroutine ChooseFixedEndNodes

  subroutine CountFixingNodes( Set, SetNo )
    integer :: Set(:)
    integer :: SetNo
    integer :: i,j,nplus,nminus
    if ( ParEnv % PEs == 1 ) then
      nplus = count( Set > 0  )
      nminus = count( Set < 0  )
    else
      nplus = 0
      nminus = 0
      do i=1,nsize
        if ( Solver % Matrix % ParallelInfo % NeighbourList(i) % Neighbours(1) &
            /= ParEnv % MyPe ) cycle
        if ( Set(i) > 0 ) then
          nplus = nplus + 1
        else if ( Set(i) < 0 ) then
          nminus = nminus + 1
        end if
      end do
      nplus = ParallelReduction( nplus )
      nminus = ParallelReduction( nminus )
    end if
    call Info(Caller,'Set'//I2S(SetNo)//' : '&
        //I2S(nplus)//' +nodes and ' &
        //I2S(nminus)//' -nodes')
  end subroutine CountFixingNodes
  subroutine LocalPotMatrix( Element, n, nd )
    integer :: n, nd
    type(element_t), target :: Element
    real(kind=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,Weight,LoadAtIP,GradAtIp(3), &
        FixAtIp, AbsGradAtIp, CondAtIp(3), DistGradAtIp(3),AbsDistGradAtIp,AbsCondAtIp
    real(kind=dp) :: STIFF(nd,nd), FORCE(nd), NodalPot(n), &
        NodalFix(n), NodalDist(n),DotProd, ElCond(n)
    logical :: Stat,Found, GotElCond
    integer :: i,t,p,q
    type(GaussIntegrationPoints_t) :: IP
    type(valuelist_t), pointer :: Material
    call GetElementNodes( acc_nodes_pot )
    STIFF = 0._dp
    FORCE = 0._dp
    Material => GetMaterial( Element )
    if ( UseUnityCond ) then
      GotElCond = .false.
    else
      ElCond(1:n) = GetReal( Material, CondName, GotElCond )
    end if
    if ( CoilParts == 1 ) then
      NodalPot(1:n) = PotVar % Values( Perm( Element % NodeIndexes ) )
    else
      if ( minval( PotSelect % Values( PotSelect % Perm(Element % NodeIndexes)) ) > 0.0_dp ) then
        NodalPot(1:n) = PotVarB % Values( Perm( Element % NodeIndexes ) )
      else
        NodalPot(1:n) = PotVarA % Values( Perm( Element % NodeIndexes ) )
      end if
    end if
    if ( iter > 1 ) then
      NodalFix(1:n) = FixVar % Values( Perm( Element % NodeIndexes(1:n) ) )
    end if
    if ( UseDistance ) then
      NodalDist(1:n) = DistVar % Values( DistVar % Perm( Element % NodeIndexes(1:n) ) )
    end if
    IP = GaussPoints( Element )
    do t=1,IP % n
      stat = ElementInfo( Element, acc_nodes_pot, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx )
      CondAtIp = 1.0_dp
      if ( iter > 1 ) then
        FixAtIp = sum( Basis(1:n) * NodalFix(1:n) )
        if ( CoilAnisotropic ) then
          do i=1,3
            GradAtIP(i) = sum( dBasisdx(1:n,i) * NodalPot(1:n) )
          end do
          AbsGradAtIp = sqrt( sum( GradAtIp ** 2 ) )
          CondAtIp = abs( GradAtIp ) / AbsGradAtIp
        end if
      else
        FixAtIp = 1.0_dp
      end if
      if ( UseDistance ) then
        do i=1,3
          DistGradAtIP(i) = sum( dBasisdx(1:n,i) * NodalDist(1:n) )
        end do
        AbsDistGradAtIp = sqrt( sum( DistGradAtIp ** 2 ) )
        AbsCondAtIp = sqrt(sum( CondAtIP**2 ))
        DotProd = sum( CondAtIp * DistGradAtIp ) / ( AbsDistGradAtIp * sqrt(3.0) )
        if ( AbsDistGradAtIp > 1.0_dp ) then
          DistGradAtIp = DistGradAtIp / AbsDistGradAtIp
        end if
        CondAtIp = CondAtIp - DotProd * DistGradAtIp
      end if
      if ( GotElCond ) then
        CondAtIP = sum( Basis(1:n) * ElCond(1:n) ) * FixAtIp * CondAtIp
      else
        CondAtIP = FixAtIp * CondAtIp
      end if
      Weight = IP % s(t) * DetJ
      do p=1,nd
        do q=1,nd
          do i=1,dim
            STIFF(p,q) = STIFF(p,q) + Weight * &
                CondAtIp(i) * dBasisdx(p,i) * dBasisdx(q,i)
          end do
        end do
      end do
    end do
    call DefaultUpdateEquations(STIFF,FORCE)
  end subroutine LocalPotMatrix

  subroutine LocalFixMatrix( Element, n, nd )
    integer :: n, nd
    type(element_t), target :: Element
    real(kind=dp) :: diff_coeff(n), AbsGradAtIp, Weight, Dreg, AbsCondAtIp
    real(kind=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP,GradAtIp(3),DistGradAtIp(3)
    real(kind=dp) :: STIFF(nd,nd), FORCE(nd),NodalPot(nd),NodalDist(nd)
    logical :: Stat,Found
    integer :: i,t,p,q
    type(GaussIntegrationPoints_t) :: IP
    call GetElementNodes( acc_nodes_fix )
    STIFF = 0._dp
    FORCE = 0._dp
    if ( CoilParts == 1 ) then
      NodalPot(1:n) = PotVar % Values( Perm( Element % NodeIndexes ) )
    else
      if ( minval( PotSelect % Values( PotSelect % Perm(Element % NodeIndexes)) ) > 0.0_dp ) then
        NodalPot(1:n) = PotVarB % Values( Perm( Element % NodeIndexes ) )
      else
        NodalPot(1:n) = PotVarA % Values( Perm( Element % NodeIndexes ) )
      end if
    end if
    if ( UseDistance ) then
       NodalDist(1:n) = DistVar % Values( DistVar % Perm( Element % NodeIndexes ) )
    end if
    diff_coeff(1:n) = GetReal(Params,'CFix Diffusion',Found)
    IP = GaussPoints( Element )
    do t=1,IP % n
      stat = ElementInfo( Element, acc_nodes_fix, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )
      LoadAtIP = 1.0_dp
      do i=1,3
         GradAtIP(i) = sum( dBasisdx( 1:nd, i) * NodalPot(1:nd) )
      end do
      if ( UseDistance ) then
         do i=1,3
            DistGradAtIP(i) = sum( dBasisdx(1:n,i) * NodalDist(1:n) )
         end do
         GradAtIp = GradAtIp - sum( GradAtIp * DistGradAtIp ) * DistGradAtIp
      end if
      AbsGradAtIp = sqrt( sum( GradAtIp**2) )
      Dreg = sum( Basis(1:nd) * Diff_Coeff(1:nd) )
      Weight = IP % s(t) * DetJ
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + Weight * &
             Dreg * matmul( dBasisdx, transpose( dBasisdx ) )
      do p=1,nd
        do q=1,nd
          STIFF(p,q) = STIFF(p,q) + Weight * AbsGradAtIp * Basis(q) * Basis(p)
        end do
      end do
      FORCE(1:nd) = FORCE(1:nd) + Weight * LoadAtIP * Basis(1:nd)
    end do
    call DefaultUpdateEquations(STIFF,FORCE)
  end subroutine LocalFixMatrix

  subroutine LocalFluxMatrix( Element, n, nd, dimi  )
    integer :: n, nd, dimi
    type(element_t), target :: Element
    real(kind=dp) :: Basis(nd),dBasisdx(nd,3),Weight,DetJ,LoadAtIP,FixAtIp, &
         GradAtIp(3),AbsGradAtIp,CondAtIp(3),DistGradAtIp(3),AbsDistGradAtIp, &
         AbsCondAtIp,DotProd
    real(kind=dp) :: STIFF(nd,nd), FORCE(nd),NodalPot(nd),&
         NodalFix(nd),NodalDist(nd),Elcond(nd)
    logical :: Stat, GotElCond, Erroneous
    integer :: i,t,p,q
    integer :: pivot(nd)
    type(GaussIntegrationPoints_t) :: IP
    type(valuelist_t), pointer :: Material
    Material => GetMaterial( Element )
    if ( UseUnityCond ) then
      GotElCond = .false.
    else
      ElCond(1:n) = GetReal( Material, CondName, GotElCond )
    end if
    call GetElementNodes( acc_nodes_flux )
    STIFF = 0._dp
    FORCE = 0._dp
    if ( CoilParts == 1 ) then
      NodalPot(1:n) = PotVar % Values( Perm( Element % NodeIndexes ) )
    else
      if ( minval( PotSelect % Values( PotSelect % Perm(Element % NodeIndexes)) ) > 0.0_dp ) then
        NodalPot(1:n) = PotVarB % Values( Perm( Element % NodeIndexes ) )
      else
        NodalPot(1:n) = PotVarA % Values( Perm( Element % NodeIndexes ) )
      end if
    end if
    if ( FixConductivity ) then
      NodalFix(1:n) = FixVar % Values( Perm( Element % NodeIndexes ) )
    end if
    if ( UseDistance ) then
      NodalDist(1:n) = DistVar % Values( DistVar % Perm( Element % NodeIndexes ) )
    end if
    IP = GaussPoints( Element )
    do t=1,IP % n
      stat = ElementInfo( Element, acc_nodes_flux, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )
      do i=1,3
        GradAtIP(i) = sum( dBasisdx(1:n,i) * NodalPot(1:n) )
      end do
      CondAtIp = 1.0_dp
      if ( FixConductivity ) then
        FixAtIp = sum( Basis(1:n) * NodalFix(1:n) )
        if ( CoilAnisotropic ) then
          AbsGradAtIp = sqrt( sum( GradAtIp ** 2 ) )
          CondAtIp = abs( GradAtIp ) / AbsGradAtIp
        end if
      else
        FixAtIp = 1.0_dp
      end if
      if ( UseDistance ) then
        do i=1,3
          DistGradAtIP(i) = sum( dBasisdx(1:n,i) * NodalDist(1:n) )
        end do
        AbsDistGradAtIp = sqrt( sum( DistGradAtIp ** 2 ) )
        AbsCondAtIp = sqrt(sum( CondAtIP**2 ))
        DotProd = sum( CondAtIp * DistGradAtIp ) / ( AbsDistGradAtIp * sqrt(3.0) )
        if ( AbsDistGradAtIp > 1.0_dp ) then
          DistGradAtIp = DistGradAtIp / AbsDistGradAtIp
        end if
        CondAtIp = CondAtIp - DotProd * DistGradAtIp
      end if
      if ( GotElCond ) then
        CondAtIP = sum( Basis(1:n) * ElCond(1:n) ) * FixAtIp * CondAtIp
      else
        CondAtIP = FixAtIp * CondAtIp
      end if
      LoadAtIp = CondAtIp(dimi) * GradAtIp(dimi)
      Weight = IP % s(t) * DetJ
      do p=1,nd
        do q=1,nd
          STIFF(p,q) = STIFF(p,q) + Weight * Basis(q) * Basis(p)
        end do
      end do
      FORCE(1:nd) = FORCE(1:nd) + Weight * LoadAtIP * Basis(1:nd)
    end do
    call DefaultUpdateEquations(STIFF,FORCE)
    if ( associated( FluxVarE ) ) then
      call LUdecomp(STIFF,nd,pivot,Erroneous)
      if (Erroneous) call Fatal('LocalFluxMatrix', 'LU-decomposition fails')
      call LUSolve(nd,STIFF,FORCE,pivot)
      FluxVarE % Values(FluxVarE % Perm(Element % DGIndexes(1:nd))) = FORCE(1:nd)
    end if
  end subroutine LocalFluxMatrix

  subroutine LocalCorrMatrix( Element, n, nd )
    integer :: n, nd
    type(element_t), target :: Element
    real(kind=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,Weight,L(3),Lfix(3,nd)
    real(kind=dp) :: STIFF(nd,nd), FORCE(nd)
    logical :: Stat,Found
    integer :: i,j,t,p,q
    type(GaussIntegrationPoints_t) :: IP
    type(valuelist_t), pointer :: Material
    call GetElementNodes( acc_nodes_corr )
    STIFF = 0._dp
    FORCE = 0._dp
    do i=1,n
      j = FluxVarE % Perm( Element % DGIndexes(i) )
      Lfix(1:3,i) = FluxVarE % Values( 3*j-2: 3*j )
    end do
    IP = GaussPoints( Element )
    do t=1,IP % n
      stat = ElementInfo( Element, acc_nodes_corr, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx )
      Weight = IP % s(t) * detJ
      L = matmul(Lfix(1:3,1:n),Basis(1:n))
      do p=1,n
        do q=1,n
          STIFF(p,q) = STIFF(p,q) + Weight * sum(dBasisdx(q,:)*dBasisdx(p,:))
        end do
        FORCE(p) = FORCE(p) + sum(L * dBasisdx(p,:)) * Weight
      end do
    end do
    call DefaultUpdateEquations(STIFF,FORCE)
  end subroutine LocalCorrMatrix

  subroutine LocalCorrCurrent( Element, n, nd )
    integer :: n, nd
    type(element_t), target :: Element
    real(kind=dp) :: Basis(nd),dBasisdx(nd,3),NodalPot(nd),DetJ,Weight
    real(kind=dp) :: STIFF(nd,nd), FORCE(nd,3),x(nd)
    integer :: pivot(nd)
    logical :: Stat,Found,Erroneous
    integer :: i,j,k,t,p,q
    type(GaussIntegrationPoints_t) :: IP
    call GetElementNodes( acc_nodes_corr2 )
    STIFF = 0._dp
    FORCE = 0._dp
    NodalPot(1:n) = Solver % Variable % Values( Perm( Element % NodeIndexes ) )
    IP = GaussPoints( Element )
    do t=1,IP % n
      stat = ElementInfo( Element, acc_nodes_corr2, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx )
      Weight = IP % s(t) * detJ
      do p=1,n
        do q=1,n
          STIFF(p,q) = STIFF(p,q) + Weight * Basis(p) * Basis(q)
        end do
        do k=1,3
          FORCE(p,k) = FORCE(p,k) + Weight * Basis(p) * sum(dBasisdx(1:n,k)*NodalPot(1:n))
        end do
      end do
    end do
    call LUdecomp(STIFF,n,pivot,Erroneous)
    if (Erroneous) call Fatal('LocalCorrCurrent', 'LU-decomposition fails')
    do k=1,3
      x = FORCE(1:n,k)
      call LUSolve(n,STIFF,x,pivot)
      do i=1,n
        j = FluxVarE % Perm( Element % DGIndexes(i) )
        FluxVarE % Values( 3*(j-1)+k) = FluxVarE % Values( 3*(j-1)+k) - x(i)
      end do
    end do
  end subroutine LocalCorrCurrent

  subroutine ScalePotential()
    real(kind=dp) :: InitialCurrent,possum, negsum, sumerr
    integer :: i,j,k,Coil,nsize, posi,negi, sgn
    logical :: DoIt, Fail
    real(kind=dp) :: Coeff
    call Info(Caller,'Performing scaling of potential for desired current for '//I2S(NoCoils)//' coil',level=30)
    nsize = size( LoadVar % Perm )
    Fail = .false.
    do Coil = 1, NoCoils
      if ( NoCoils > 1 ) then
        call Info(Caller,'Scaling coil number: '//I2S(Coil))
      end if
      possum = 0.0_dp
      negsum = 0.0_dp
      posi = 0
      negi = 0
      do i=1,nsize
        j = LoadVar % Perm(i)
        if ( j == 0 ) cycle
        if ( ParEnv % PEs > 1 ) then
          if ( Solver % Matrix % ParallelInfo % Neighbourlist(j) % Neighbours(1) &
              /= ParEnv % MyPe ) cycle
        end if
        if ( NoCoils > 1 ) then
          if ( CoilIndex(i) /= Coil ) cycle
        end if
        sgn = Set(j)
        if ( sgn == 0 ) cycle
        if ( modulo(sgn,10) == 2  ) then
          possum = possum + LoadVar % Values(j)
          posi = posi + 1
        else if ( modulo(sgn,10)-10 == -2 ) then
          negsum = negsum + LoadVar % Values(j)
          negi = negi + 1
        end if
      end do
      if ( ParEnv % PEs > 1 ) then
        possum = ParallelReduction( possum )
        negsum = ParallelReduction( negsum )
      end if
      write (Message,'(A,I6,ES12.4)') 'Positive coil currents:',posi,possum
      call Info(Caller,Message,level=12)
      write (Message,'(A,I6,ES12.4)') 'Negative coil currents:',negi,negsum
      call Info(Caller,Message,level=12)
      DoIt = .true.
      if ( abs( possum ) < epsilon( possum ) ) then
        call Warn(Caller,'No positive current sources on coil '//I2S(Coil)//' end!')
        DoIt = .false.
      end if
      if ( abs( negsum ) < epsilon( negsum ) ) then
        call Warn(Caller,'No negative current sources on coil '//I2S(Coil)//' end!')
        DoIt = .false.
      end if
      if (.not. DoIt) then
        call Warn(Caller,'Crappy potentials in coil '//I2S(Coil))
        Fail = .true.
        cycle
      end if
      sumerr = 2.0 * abs( abs( possum ) - abs( negsum ) )  /(abs( possum ) + abs( negsum ) )
      write (Message,'(A,ES12.4)') 'Discrepancy of start and end coil currents: ',sumerr
      call Info(Caller,Message,level=7)
      InitialCurrent = ( possum - negsum ) / 2.0
      write( Message,'(A,ES12.4)') 'Initial coil current for coil '&
          //I2S(Coil)//':',InitialCurrent
      call Info(Caller,Message,level=5)
      Coeff = DesiredCoilCurrent(Coil) / InitialCurrent
      write( Message,'(A,ES12.4)') 'Coil potential multiplier:',Coeff
      call Info(Caller,Message,level=5)
      if ( NoCoils == 1 ) then
        PotVar % Values = Coeff * PotVar % Values
      else
        do i=1,nsize
          if ( CoilIndex(i) == Coil ) then
            j = PotVar % Perm(i)
            if ( j == 0 ) cycle
            PotVar % Values( j ) = Coeff * PotVar % Values( j )
          end if
        end do
      end if
    end do
    if (Fail) then
      call Fatal(Caller,'Scaling of potential failed!')
    end if
    call Info(Caller,'Scaling of potential finished',level=30)
  end subroutine ScalePotential

  subroutine NormalizeCurrentDensity()
    integer :: n, nd, dimi, Active
    type(element_t), pointer :: Element
    real(kind=dp) :: Weight,DetJ,LocalCurr(3),AbsCurr,ScaleCurr,&
        TotCurr,TotVol,TargetDensity,alpha
    real(kind=dp), allocatable :: Basis(:),NodalCurr(:,:)
    logical :: Stat
    integer :: i,j,elem,t,Coil
    type(GaussIntegrationPoints_t) :: IP
    if ( .not. any(NormalizeCoil ) ) then
      call Info(Caller,'No normalization of current requested!',level=20)
      return
    end if
    call Info(Caller,'Normalizing current density (as for stranded coil)',level=5)
    n = Mesh % MaxElementNodes
    allocate( Basis(n), NodalCurr(3,n) )
    FluxVar => VariableGet( Mesh % Variables,'CoilCurrent')
    FluxVarE => VariableGet( Mesh % Variables,'CoilCurrent E')
    if ( .not. ( associated( FluxVar ) .or. associated( FluxVarE ) ) ) then
      call Fatal(Caller,'CoilCurrent (nodal and elemental) not associated!')
    end if
    do Coil = 1, NoCoils
      if ( .not. NormalizeCoil(Coil) ) cycle
      if ( .not. GotDens(Coil) .and. GotCurr(Coil) ) then
        TotVol = 0.0_dp
        TotCurr = 0.0_dp
        NodalCurr = 0.0_dp
        LocalCurr = 0.0_dp
        Active = GetNOFActive()
          do elem=1,Active
            Element => GetActiveElement(elem)
            n  = GetElementNOFNodes()
            nd = GetElementNOFDOFs()
            if ( NoCoils > 1 ) then
              if ( any( CoilIndex( Element % NodeIndexes ) /= Coil ) ) cycle
            end if
            call GetElementNodes( acc_nodes_norm )
            if ( associated( FluxVarE ) ) then
              do dimi=1,dim
                NodalCurr(dimi,1:n) = FluxVarE % Values( dim*(FluxVarE % Perm( Element % DGIndexes )-1) + dimi )
              end do
            else
              do dimi=1,dim
                NodalCurr(dimi,1:n) = FluxVar % Values( dim*(FluxVar % Perm( Element % NodeIndexes )-1) + dimi )
              end do
            end if
            IP = GaussPoints( Element )
            do t=1,IP % n
              stat = ElementInfo( Element, acc_nodes_norm, IP % U(t), IP % V(t), &
                  IP % W(t), detJ, Basis )
              do dimi=1,dim
                LocalCurr(dimi) = sum( Basis(1:n) * NodalCurr(dimi,1:n) )
              end do
              AbsCurr = sqrt( sum( LocalCurr(1:dim)**2 ) )
              Weight = IP % s(t) * DetJ
              TotVol = TotVol + Weight
              TotCurr = TotCurr + Weight * AbsCurr
            end do
          end do
        ! end of averaging sweep
        if ( ParEnv % PEs > 1 ) then
          TotCurr = ParallelReduction( TotCurr )
          TotVol = ParallelReduction( TotVol )
        end if
        TargetDensity = TotCurr / TotVol
        write( Message,'(A,ES12.4)') 'Average current density:',TargetDensity
        call Info(Caller,Message)
      else
        TargetDensity = DesiredCurrentDensity(Coil)
        write( Message,'(A,ES12.4)') 'Desired current density:',TargetDensity
      end if
      if ( associated( FluxVar ) ) then
        do i = 1, Solver % Mesh % NumberOfNodes
          if ( NoCoils > 1 ) then
            if ( CoilIndex(i) /= Coil ) cycle
          end if
          j = FluxVar % Perm(i)
          if ( j == 0 ) cycle
          LocalCurr(1:dim) = FluxVar % Values( dim*(j-1)+1: dim*(j-1)+dim )
          alpha = CoilHelicity(Coil)
          if (abs(alpha) > tiny(alpha) ) then
            AbsCurr = sqrt( sum( LocalCurr(1:dim) ** 2 ) )
            LocalCurr = LocalCurr + alpha * AbsCurr * CoilNormals(Coil,:)
          end if
          AbsCurr = sqrt( sum( LocalCurr(1:dim) ** 2 ) )
          if ( AbsCurr > tiny( AbsCurr ) ) then
            ScaleCurr = TargetDensity / AbsCurr
            FluxVar % Values( dim*(j-1)+1: dim*(j-1)+dim ) = &
                ScaleCurr * LocalCurr(1:dim)
          end if
        end do
      end if
      if ( associated( FluxVarE ) ) then
        Active = GetNOFActive()
        do elem=1,Active
            Element => GetActiveElement(elem)
            n  = GetElementNOFNodes()
            if ( NoCoils > 1 ) then
              if ( any( CoilIndex( Element % NodeIndexes ) /= Coil ) ) cycle
            end if
            do i=1,n
              j = FluxVarE % Perm( Element % DGIndexes(i) )
              if ( j == 0 ) cycle
              LocalCurr(1:dim) = FluxVarE % Values( dim*(j-1)+1: dim*(j-1)+dim )
              alpha = CoilHelicity(Coil)
              if (abs(alpha) > tiny(alpha) ) then
                AbsCurr = sqrt( sum( LocalCurr(1:dim) ** 2 ) )
                LocalCurr = LocalCurr + alpha * AbsCurr * CoilNormals(Coil,:)
              end if
              AbsCurr = sqrt( sum( LocalCurr(1:dim) ** 2 ) )
              if ( AbsCurr > tiny( AbsCurr ) ) then
                ScaleCurr = TargetDensity / AbsCurr
                FluxVarE % Values( dim*(j-1)+1: dim*(j-1)+dim ) = &
                    ScaleCurr * LocalCurr(1:dim)
              end if
          end do
        end do
      end if
    end do
    if ( InfoActive(30) ) then
      if ( associated( FluxVar ) ) then
        call VectorValuesRange(FluxVar % Values,size(FluxVar % Values),'Current density')
      end if
      if ( associated( FluxVarE ) ) then
        call VectorValuesRange(FluxVarE % Values,size(FluxVarE % Values),'Current density e')
      end if
    end if
  end subroutine NormalizeCurrentDensity

  ! -- user-function entries --------------------------------------------
  function coa_cpot_eval(h_model, n) result(f) bind(c, name="coa_cpot_eval")
    type(c_ptr), value :: h_model
    integer(c_int), value :: n
    real(c_double) :: f
    type(element_t), pointer :: Element
    real(kind=dp) :: xmin
    call c_f_pointer(h_model, Model)
    if ( .not. acc_cp_visited ) then
      acc_cp_pota => VariableGet( Model % Mesh % Variables,'CoilPot' )
      acc_cp_potp => acc_cp_pota
      acc_cp_potb => VariableGet( Model % Mesh % Variables,'CoilPotB' )
      acc_cp_potsel => VariableGet( Model % Mesh % Variables,'PotSelect' )
      acc_cp_visited = .true.
    end if
    if ( associated( acc_cp_potsel ) ) then
      Element => Model % CurrentElement
      if ( .not. associated( Element, acc_cp_prevelem ) ) then
        xmin = minval( acc_cp_potsel % Values( acc_cp_potsel % Perm(Element % NodeIndexes ) ) )
        if ( xmin > 0.0 ) then
          acc_cp_potp => acc_cp_potb
        else
          acc_cp_potp => acc_cp_pota
        end if
      end if
    end if
    f = acc_cp_potp % Values( acc_cp_potp % Perm(n) )
  end function coa_cpot_eval

  subroutine coa_cpotn_eval(h_model, n, f) bind(c, name="coa_cpotn_eval")
    type(c_ptr), value :: h_model
    integer(c_int) :: n
    real(c_double) :: f
    logical :: Found, Stat
    integer :: i, m
    type(element_t), pointer :: Element
    type(GaussIntegrationPoints_t) :: IP
    real(kind=dp) :: u, v, w, detJ
    call c_f_pointer(h_model, Model)
    if ( .not. acc_cn_visited ) then
      acc_cn_potb => VariableGet( Model % Mesh % Variables,'CoilPotB' )
      acc_cn_pota => VariableGet( Model % Mesh % Variables,'CoilPot' )
      acc_cn_potp => acc_cn_pota
      acc_cn_potsel => VariableGet( Model % Mesh % Variables,'PotSelect' )
      do i=1,Model % NumberOfSolvers
        acc_cn_desdens = ListGetCReal( Model % Solvers(i) % Values,&
            'Desired Current Density',Found )
        if ( Found ) exit
      end do
      if (.not. Found ) acc_cn_desdens = 1.0_dp
      n = Model % Solver % Mesh % MaxElementNodes
      allocate( acc_cn_basis(n), acc_cn_dbasis(n,3), acc_cn_nodalpot(n), &
          acc_nodes_cpotn % x(n), acc_nodes_cpotn % y(n), acc_nodes_cpotn % z(n) )
      acc_cn_visited = .true.
    end if
    Element => Model % CurrentElement
    if ( .not. associated( Element, acc_cn_prevelem ) ) then
      m = GetElementNOFNodes()
      acc_nodes_cpotn % x(1:m) = Model % Mesh % Nodes % x(Element % NodeIndexes)
      acc_nodes_cpotn % y(1:m) = Model % Mesh % Nodes % y(Element % NodeIndexes)
      acc_nodes_cpotn % z(1:m) = Model % Mesh % Nodes % z(Element % NodeIndexes)
      IP = GaussPoints( Element )
      u = sum( IP % U ) / IP % n
      v = sum( IP % V ) / IP % n
      w = sum( IP % W ) / IP % n
      stat = ElementInfo( Element, acc_nodes_cpotn, u, v, w, detJ, &
          acc_cn_basis, acc_cn_dbasis )
      if ( associated( acc_cn_potsel ) ) then
        acc_cn_xmin = minval( acc_cn_potsel % Values( acc_cn_potsel % Perm(Element % NodeIndexes ) ) )
        if ( acc_cn_xmin > 0.0_dp ) then
          acc_cn_potp => acc_cn_potb
        else
          acc_cn_potp => acc_cn_pota
        end if
      end if
      acc_cn_nodalpot(1:m) = acc_cn_potp % Values( acc_cn_potp % Perm(Element % NodeIndexes) )
      do i=1,3
        acc_cn_gradpot(i) = sum( acc_cn_dbasis(1:m,i) * acc_cn_nodalpot(1:m) )
      end do
      acc_cn_normcoeff = acc_cn_desdens / sqrt( sum(acc_cn_gradpot**2) )
      acc_cn_prevelem => Element
    end if
    f = acc_cn_normcoeff * acc_cn_potp % Values( acc_cn_potp % Perm(n) )
  end subroutine coa_cpotn_eval

end module coil_acc
