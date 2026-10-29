! TopoOpt solver-entry port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/TopoOpt.F90, compiled against the built
! .mod files of the target configure (A2 shipped mechanism, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (TYPE(C_PTR) handles,
! VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All Fortran-object
! traffic — TYPE(...) handles, OPTIONAL/POINTER dummies, assumed-shape
! descriptors, hidden CHARACTER lengths, array-pointer List results — stays
! inside Fortran, so the Rust side is pure C calls. The former CONTAINS
! internals (ObjectiveGradients, UpdateDensities, ApplyFilter,
! ApplyTopologyFilter, NormalizeFilter, CreateSimpleFilter,
! CreateDistanceFilter, DefineTopologyBW, DefineTopologyZeroLevel,
! LocalMatrix, ApplyPDEFilter) are transcribed verbatim as module
! procedures sharing the module SAVE state below, which mirrors the
! original's SAVE list one for one (TimesVisited excepted: that latch
! lives in the Rust static and reaches Fortran as an explicit argument).
! Order-sensitive reductions (MATMUL, SUM), formatted I/O (PRINT/WRITE),
! whole-array updates, and the UpdateDensities bisection loop all stay in
! Fortran, so bit-identity holds by construction wherever the call
! sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module toa_acc
  use iso_c_binding
  use DefUtils
  use MeshBasics
  implicit none

  ! -- entry handles ----------------------------------------------------
  type(model_t), pointer, save :: Model => null()
  type(solver_t), pointer, save :: Solver => null()
  type(valuelist_t), pointer, save :: Params => null()

  ! -- mirror of TopoOpt's SAVE list (TimesVisited lives in Rust) -------
  type(matrix_t), pointer, save :: Fmat => null()
  real(kind=dp), pointer, save :: xTopo(:) => null()
  real(kind=dp), pointer, save :: xPhys(:) => null()
  real(kind=dp), pointer, save :: xMult(:) => null()
  real(kind=dp), pointer, save :: ce(:) => null()
  real(kind=dp), pointer, save :: dc(:) => null()
  real(kind=dp), pointer, save :: dv(:) => null()
  real(kind=dp), pointer, save :: dv0(:) => null()
  real(kind=dp), pointer, save :: bw(:) => null()
  real(kind=dp), pointer, save :: zeroset(:) => null()
  real(kind=dp), allocatable, save :: local_sol_array(:,:)
  ! Scratch `Found` flag: the original body and its internals share one
  ! host-associated `Found`; pieces that need their own declare a local
  ! shadow, exactly as the original's scoping did.
  logical, save :: Found
  real(kind=dp), allocatable, save :: local_sol(:)
  real(kind=dp), allocatable, save :: local_act(:)
  integer, pointer, save :: ElemPerm(:) => null()
  integer, allocatable, save :: SumPerm(:)
  logical, allocatable, save :: InterfaceNode(:)
  type(solver_t), pointer, save :: PhysSolver => null()
  type(mesh_t), pointer, save :: Mesh => null()
  type(variable_t), pointer, save :: uVar => null()
  type(variable_t), pointer, save :: aVar => null()
  type(nodes_t), save :: Nodes
  type(element_t), pointer, save :: Element => null()
  character(:), allocatable, save :: filterMethod, filterType
  integer, save :: Niter, dim, Nnodes, Nelems, dofs, nsize, nPer
  real(kind=dp), save :: wmin, gt, obj
  real(kind=dp), save :: Diff(3)
  logical, save :: PdeFilter, SimpleFilter, Csymmetry, SolveAdj, PhysSym
  logical, save :: SkipInterface
  ! Per-call scalars of the original body, promoted to module state so
  ! they can cross the wrapper boundaries within one entry call; they are
  ! recomputed from the same ValueList on every call, hence unobservable.
  real(kind=dp), save :: penal, volFrac, emin, efrac
  ! ApplyPDEFilter's `DoMatrix = .TRUE.` latch (implicit SAVE in the
  ! original: first call assembles, later calls only refresh the r.h.s.).
  logical, save :: acc_domatrix = .true.

  character(*), parameter :: Caller = 'TopoOpt'

contains

  ! -- entry state ------------------------------------------------------
  subroutine toa_begin(h_model, h_solver) bind(c, name="toa_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, Model)
    call c_f_pointer(h_solver, Solver)
    Params => Solver % Values
  end subroutine toa_begin

  ! -- TopoOpt_init0 body (verbatim) ------------------------------------
  subroutine toa_init0() bind(c, name="toa_init0")
    logical :: Found
    integer :: i
    Params => Solver % Values
    call ListAddNewInteger( Params,'Primary Solver Index', 1 )
    i = ListGetInteger( Params,'Primary Solver Index' )
    if (ListGetLogical( Params,'Solve Adjoint Problem', Found ) ) then
      call ListAddString( Params,NextFreeKeyword('Exported Variable ',Params),'Adjoint')
    end if
  end subroutine toa_init0

  ! -- TopoOpt_init body (verbatim) -------------------------------------
  subroutine toa_init() bind(c, name="toa_init")
    logical :: Found, HaveField
    Params => Solver % Values
    ! These automatically allocate elemental variables that are then created by library
    ! even before visiting the subroutine below.
    call ListAddString( Params,NextFreeKeyword('Exported Variable ',Params),'-elem topo rho' )
    call ListAddString( Params,NextFreeKeyword('Exported Variable ',Params),'-elem topo mult' )
    call ListAddString( Params,NextFreeKeyword('Exported Variable ',Params),'-elem topo ce' )
    call ListAddString( Params,NextFreeKeyword('Exported Variable ',Params),'-elem topo dc' )
    if ( ListGetLogical( Params,'Create BW Topology', Found ) ) then
      call ListAddString( Params,NextFreeKeyword('Exported Variable ',Params),'-elem topo bw' )
    end if
    HaveField = ( ListGetString( Params,'Filter Type', Found ) == 'pde' )
    if ( HaveField ) then
      call ListAddNewString( Params,'Variable','xNodal')
    else
      ! Add a global variable to store the norm to if no variable present.
      call ListAddNewString( Params,'Variable','-nooutput -global topoopt_nrm')
    end if
    if ( ListGetLogical( Params,'Create Zero Levelset', Found ) ) then
      if ( HaveField ) then
        call ListAddString( Params,NextFreeKeyword('Exported Variable ',Params),'topo levelset' )
      else
        call Warn('TopoOpt_init','Only PDE filter can create levelset field!')
      end if
    end if
  end subroutine toa_init

  ! -- entry open (verbatim Info pair) ----------------------------------
  subroutine toa_info_open() bind(c, name="toa_info_open")
    call Info(Caller,'-------------------------------------------')
    call Info(Caller,'Updating density for topology optimization')
  end subroutine toa_info_open

  ! -- TimesVisited == 0 block (verbatim) -------------------------------
  subroutine toa_first_touch() bind(c, name="toa_first_touch")
    integer :: i, j, n
    logical :: Found, ElemField
    real(kind=dp) :: val
    type(Variable_t), pointer :: pVar

    Mesh => Solver % Mesh
    dim = Mesh % MeshDim
    Nnodes = Mesh % NumberOfNodes
    Nelems = Mesh % NumberOfBulkElements

    CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .or. &
        CurrentCoordinateSystem() == CylindricSymmetric )

    call Info(Caller,'Number of nodes: '//I2S(Nnodes))
    call Info(Caller,'Number of bulk elements: '//I2S(Nelems))

    filterMethod = ListGetString( Params,'Filter Method',UnfoundFatal=.TRUE.)
    call Info(Caller,'Filter Method: '//trim(filterMethod))

    filterType = ListGetString( Params,'Filter Type',UnfoundFatal=.TRUE.)
    call Info(Caller,'Filter Type: '//trim(filterType))

    PdeFilter = .false.
    SimpleFilter = .false.
    select case ( filterType )
    case('simple')
      SimpleFilter = .true.
    case('pde')
      PdeFilter = .true.
    case('distance')
    case default
      call Fatal(Caller,'Unknown filter type: '//trim(filterType))
    end select

    if (PdeFilter .and. .not. associated(Solver % Matrix)) then
      call Fatal(Caller,'Pde Filter requires field variable & matrix to be associated!')
    end if

    ! This is not generic. We assume the stress solver to be the 1st solver for now.
    i = ListGetInteger( Params,'Primary Solver Index' )
    PhysSolver => Model % Solvers(i)
    if (.not. ListGetLogical(PhysSolver % Values,'Local Matrix Storage',Found ) ) then
      call Fatal(Caller,'Primary solver should have active "Local Matrix Storage"')
    end if

    uVar => PhysSolver % Variable
    dofs = uVar % dofs

    ! Get adjoint variable, or if not requested use primary variable for adjoint as well.
    SolveAdj = ListGetLogical( Params,'Solve Adjoint Problem', Found )
    if ( SolveAdj ) then
      aVar => VariableGet( Mesh % Variables,'Adjoint' )
      if (.not. associated(aVar)) call Fatal(Caller,'Could not find adjoint variable!')
    else
      aVar => uVar
    end if

    ! These fields are created also for visualization in mind!
    pVar => VariableGet( Mesh % Variables,"topo rho", UnfoundFatal = .TRUE.)
    xPhys => pVar % Values
    ElemPerm => pVar % Perm
    nsize = size(xPhys)
    call Info(Caller,'Size of elemental fields: '//I2S(nsize),level=10)

    pVar => VariableGet( Mesh % Variables,"topo mult", UnfoundFatal = .TRUE.)
    xMult => pVar % Values
    pVar => VariableGet( Mesh % Variables,"topo ce", UnfoundFatal = .TRUE.)
    ce => pVar % Values
    pVar => VariableGet( Mesh % Variables,"topo dc", UnfoundFatal = .TRUE.)
    dc => pVar % Values
    pVar => VariableGet( Mesh % Variables,"topo bw")
    if (associated(pVar)) then
      bw => pVar % Values
    else
      bw => null()
    end if

    pVar => VariableGet( Mesh % Variables,"topo levelset")
    if (associated(pVar)) then
      zeroset => pVar % Values
    else
      zeroset => null()
    end if

    ! Allocate full vectors
    allocate( xTopo(nsize) )
    xTopo = xPhys
    gt = 0.0_dp
    ce = 0.0_dp
    dc = 0.0_dp

    ! Calculate elemental volume
    allocate( dv0(nsize) )
    dv0 = 0.0_dp
    do i=1,Nelems
      j = ElemPerm(i)
      if (j==0) cycle
      Element => Mesh % Elements(i)
      call CopyElementNodesFromMesh( Nodes, Solver % Mesh, &
          Element % TYPE % NumberOfNodes, Element % NodeIndexes )
      dv0(j) = ElementSize( Element, Nodes )
      if (Csymmetry) dv0(j) = dv0(j) * sum(Nodes % x) / Element % TYPE % NumberOfNodes
    end do
    if (InfoActive(20)) then
      call VectorValuesRange(dv0,nsize,'dv0')
    end if

    nPer = ListGetInteger( Params,'Periodic PhysSolver',Found )
    if (nPer > 1) then
      allocate(SumPerm(Solver % Mesh % NumberOfBulkElements))
      PhysSym = ListGetLogical( Params,'Periodic PhysSolver Symmetric',Found )
      ElemField = .true.
      call RotationalPeriodicSumPerm(Solver, Solver % Mesh, 360.0_dp/nPer, &
          Solver % Variable % Perm, SumPerm, ElemField, PhysSym )
    end if

    ! Allocate elemental stuff
    n = Mesh % MaxElementDofs
    allocate(local_sol_array(dofs,n), local_sol(dofs*n), local_act(dofs*n))

    wmin = ListGetConstReal( Params,'Sensitivity Filter Threshold', Found )
    if (.not. Found) wmin = 1.0e-3
    SkipInterface = .false.

    if (PdeFilter ) then
      block
        real(kind=dp), pointer :: HWrk(:,:)
        integer, allocatable :: NodeCount(:)
        integer :: d1, d2, t
        Hwrk => ListGetConstRealArray( Params,'PDE Filter Diffusion Constant',UnfoundFatal = .TRUE.)

        d1 = size(Hwrk,1)
        d2 = size(Hwrk,2)

        if (d1 == 1 .and. d2 == 1 ) then
          Diff = Hwrk( 1,1)
        else if (d1 == 1 .and. d2 >= dim ) then
          Diff(1:dim) = Hwrk(1,1:dim)
        else if (d1 >= dim .and. d2 == 1 ) then
          Diff(1:dim) = Hwrk(1:dim,1)
        else
          call Fatal(Caller,'Invalid size for "PDE Filter Diffusion Constant": '//I2S(d1)//' x '//I2S(d2))
        end if

        SkipInterface = ListGetLogical( Params,'PDE Filter skip Interface',Found )
        if ( SkipInterface ) then
          allocate(NodeCount(Mesh % NumberOfNodes))
          NodeCount = 0
          do t=1,Mesh % NumberOfBulkElements
            Element => Mesh % Elements(t)
            NodeCount(Element % NodeIndexes) = NodeCount(Element % NodeIndexes) + 1
          end do
          allocate(InterfaceNode(Mesh % NumberOfNodes))
          InterfaceNode = (NodeCount < 4 )
          deallocate(NodeCount)

          t = count(InterfaceNode)
          call Info(Caller,'Number of Interface nodes: '//I2S(t),level=7)
        end if
      end block
    else
      if (ParEnv % PEs > 1 ) then
        call Fatal(Caller,'Only PDE Filter is implemented in parallel!')
      end if
      if ( SimpleFilter ) then
        FMat => CreateSimpleFilter()
        Niter = max(1,ListGetInteger( Params,'Simple Filter Iterations', Found ))
      else
        Fmat => CreateDistanceFilter()
        Niter = 1
      end if
      call NormalizeFilter(Fmat,.true.)
      val = 1.0_dp * size(Fmat % Cols) / Fmat % NumberOfRows
      write(Message,'(A,ES12.3)') 'Average number of hits in filter:',val
      call Info(Caller,Message)
    end if

    select case( FilterMethod )
    case('sensitivity')
      dv => dv0
    case('density')
      allocate(dv(nsize))
      if ( PdeFilter ) then
        call ApplyPdeFilter( dv0, dv, Diff )
      else
        call ApplyFilter( Fmat, dv0, dv, niter, Trans=.TRUE. )
      end if
      call VectorValuesRange(dv,size(dv),'dv')
    case('none')
      dv => dv0
    case default
      call Fatal(Caller,'Uknown filtering method: '//trim(FilterMethod))
    end select
  end subroutine toa_first_touch

  ! -- InfoActive(20) Kmat/uVar dump (verbatim) --------------------------
  subroutine toa_info_kmat() bind(c, name="toa_info_kmat")
    if (InfoActive(20)) then
      call VectorValuesRange(PhysSolver % Matrix % Values,size(PhysSolver % Matrix % Values),'Kmat')
      call VectorValuesRange(uVar % Values,size(uVar % Values),trim(uVar % Name))
    end if
  end subroutine toa_info_kmat

  ! -- per-call parameter reads (verbatim) ------------------------------
  subroutine toa_params() bind(c, name="toa_params")
    penal = ListGetCReal( Params,'Penalty Exponent',UnfoundFatal=.TRUE.)
    volFrac = ListGetCReal( Params,'Volume Fraction',UnfoundFatal=.TRUE.)
    emin = ListGetCReal( Params,'Minimum Relative Density',UnfoundFatal=.TRUE.)
    efrac = 1.0_dp - emin
  end subroutine toa_params

  ! -- post-params xTopo recopy (verbatim, TimesVisited == 0 only) -------
  subroutine toa_recopy_xtopo() bind(c, name="toa_recopy_xtopo")
    xTopo = xPhys
  end subroutine toa_recopy_xtopo

  ! -- GetCycleMode (verbatim logic; nT arrives from the Rust latch) -----
  function toa_cmode(nT) result(mode) bind(c, name="toa_cmode")
    integer(c_int), value :: nT
    integer(c_int) :: mode
    integer :: nCycle, nTt
    logical :: Found
    mode = 0
    nTt = nT
    nCycle = ListGetInteger( Model % Simulation,'Periodic Timesteps',Found )
    ! 1st iteration make a something to create asymmetry
    if (nCycle == 0 .or. nTt == 0 ) return
    select case( modulo(nTt,nCycle) )
    case( 0 )
      mode = 3
    case( 1 )
      mode = 1
    case default
      mode = 2
    end select
  end function toa_cmode

  ! -- obj/dc reset for cycle modes 0 and 1 (verbatim) ------------------
  subroutine toa_zero_obj_dc() bind(c, name="toa_zero_obj_dc")
    obj = 0.0_dp
    dc = 0.0_dp
  end subroutine toa_zero_obj_dc

  ! -- gradient call (verbatim) -----------------------------------------
  subroutine toa_objective_gradients() bind(c, name="toa_objective_gradients")
    call ObjectiveGradients(xPhys,ce,dc,dv0,obj)
  end subroutine toa_objective_gradients

  ! -- mid-cycle early exit note (verbatim) ------------------------------
  subroutine toa_mid_early() bind(c, name="toa_mid_early")
    call Info(Caller,'Mid of cycle, finishing early!')
  end subroutine toa_mid_early

  ! -- parallel reduction of the objective (verbatim) --------------------
  subroutine toa_parallel_obj() bind(c, name="toa_parallel_obj")
    obj = ParallelReduction( obj )
  end subroutine toa_parallel_obj

  ! -- objective report (verbatim) ---------------------------------------
  subroutine toa_report_objective() bind(c, name="toa_report_objective")
    if (InfoActive(20)) then
      call VectorValuesRange(xPhys,size(xPhys),'xPhys')
      call VectorValuesRange(ce,size(ce),'ce')
      call VectorValuesRange(dc,size(dc),'dc')
    end if
    write(Message,*) 'Objective function: ',obj
    call Info(Caller,Message)
  end subroutine toa_report_objective

  ! -- pre-filter select (verbatim) --------------------------------------
  subroutine toa_pre_filter() bind(c, name="toa_pre_filter")
    select case( FilterMethod )
    case('sensitivity')
      if ( PdeFilter ) then
        call ApplyPdeFilter( dc, dc, Diff, xTopo, wmin )
      else
        call ApplyFilter( Fmat, dc, dc, niter, xTopo, wmin )
      end if
      if (InfoActive(20)) then
        call VectorValuesRange(dc,size(dc),'dc pre')
      end if
    case('density')
      if ( PdeFilter ) then
        call ApplyPDEFilter( dc, dc, Diff )
      else
        call ApplyFilter( Fmat, dc, dc, niter, Trans=.TRUE. )
      end if
      if (InfoActive(20)) then
        call VectorValuesRange(dc,size(dc),'dc pre')
      end if
    case('none')
      call Info(Caller,'Applying no filtering')
    case default
      call Fatal(Caller,'Uknown filtering method: '//trim(FilterMethod))
    end select
  end subroutine toa_pre_filter

  ! -- density update (verbatim) -----------------------------------------
  subroutine toa_update() bind(c, name="toa_update")
    call UpdateDensities(xTopo,dc,dv,gt)
    if (InfoActive(20)) then
      call VectorValuesRange(xTopo,size(xTopo),'xTopo')
    end if
  end subroutine toa_update

  ! -- post-filter select (verbatim) -------------------------------------
  subroutine toa_post_filter() bind(c, name="toa_post_filter")
    select case( FilterMethod )
    case('density')
      if ( PdeFilter ) then
        call ApplyPDEFilter( xTopo, xPhys, Diff )
      else
        call ApplyFilter( Fmat, xTopo, xPhys, niter )
      end if
    case default
      xPhys = xTopo
    end select
  end subroutine toa_post_filter

  ! -- post physics InfoActive dumps (verbatim) ---------------------------
  subroutine toa_post_phys_info() bind(c, name="toa_post_phys_info")
    if (InfoActive(20)) then
      call VectorValuesRange(xPhys,size(xPhys),'xPhys')
      call VectorValuesRange(abs(xPhys-xTopo),size(xPhys),'dx')
    end if
  end subroutine toa_post_phys_info

  ! -- objective-as-norm assignment (verbatim) ----------------------------
  subroutine toa_assign_norm() bind(c, name="toa_assign_norm")
    if (size(Solver % Variable % Values) == 1 ) then
      Solver % Variable % Values = obj
    end if
  end subroutine toa_assign_norm

  ! -- B&W extraction (verbatim) ------------------------------------------
  subroutine toa_bw() bind(c, name="toa_bw")
    if (associated(bw)) then
      call Info(Caller,'Extracting B&W coloring',level=8)
      call DefineTopologyBW(dv0,xPhys,bw)
    end if
  end subroutine toa_bw

  ! -- zero levelset extraction (verbatim) ---------------------------------
  subroutine toa_zero_level() bind(c, name="toa_zero_level")
    if (associated(zeroset)) then
      call Info(Caller,'Extracting zero levelset function',level=8)
      call DefineTopologyZeroLevel(Solver % Variable,zeroset)
    end if
  end subroutine toa_zero_level

  ! -- label-1 tail: stiffness multiplier + periodic extension (verbatim) --
  subroutine toa_xmult() bind(c, name="toa_xmult")
    integer :: i, j
    type(Variable_t), pointer :: pVar, pVarExt
    xMult = emin + efrac * xPhys**penal
    if ( nPer > 1 ) then
      pVar => VariableGet( Mesh % Variables,"topo mult", UnfoundFatal = .TRUE.)
      pVarExt => VariableGet( Mesh % Variables,"topo mult extended", UnfoundFatal = .TRUE.)
      pVarExt % Values = 1.0_dp
      do i=1,Mesh % NumberOfBulkElements
        Element => Mesh % Elements(i)
        j = abs(SumPerm(i))
        if (j==0) cycle
        pVarExt % Values(pVarExt % Perm(i)) = pVar % Values(pVar % Perm(j))
      end do
    end if
  end subroutine toa_xmult

  !---------------------------------------------------------------------
  !> xphys: in [0,1], densities used for scaling the material properties
  !> ce: elementsize energies before scaling by density
  !> dc: gradient with respect to objective function.
  !> dv: gradient with respect to volume constraint.
  !---------------------------------------------------------------------
  subroutine ObjectiveGradients(x,ce,dc,dv,obj)
    real(kind=dp) :: x(:), ce(:), dc(:), dv(:)
    real(kind=dp) :: obj

    integer :: i,j,k,l,NoModes, NoActive, sgn, sgn1, sgn2, n, m
    logical :: UseAdjoint, Found
    type(Variable_t), pointer :: AdjSol
    real(kind=dp), allocatable :: Stiff(:,:), Force(:)
    type(Element_t), pointer :: Element
    real(kind=dp) :: spos, sneg

    n = Solver % Mesh % MaxElementNodes * dofs
    allocate(Stiff(n,n), Force(n) )

    UseAdjoint = ListGetLogical(Params,'Use Adjoint Solution',Found)
    if (UseAdjoint) then
      AdjSol => VariableGet(Solver % Mesh % Variables,trim(PhysSolver % Variable % Name)//' adjoint')
      if (.not. associated(AdjSol)) call Fatal(Caller,'Did not find Adjoint solution!')
      call Info(Caller,'Using adjoint solution: '//trim(AdjSol % Name))
      if (InfoActive(20)) then
        print *,'Adjoing interval:',minval(AdjSol % Values), maxval(AdjSol % Values), size(AdjSol % Values)
      end if
    end if

    sgn1 = LIstGetInteger(Params,'Sign A',Found )
    if (.not. Found) sgn1 = 1
    sgn2 = LIstGetInteger(Params,'Sign B',Found )
    if (.not. Found) sgn2 = 1

    NoModes = ListGetInteger(Params,'No Modes',Found )
    if (.not. Found ) then
      NoModes = PhysSolver % Variable % NumberOfConstraintModes
    end if

    if (nPer > 1) then
      NoActive = PhysSolver % NumberOfActiveElements
    else
      NoActive = Solver % NumberOfActiveElements
    end if

    !ce = 0.0_dp
    !PRINT *,'NULLIFY ce'

    spos = 0.0_dp
    sneg = 0.0_dp

    do i=1,NoActive
      if (nPer > 1 ) then
        j = PhysSolver % ActiveElements(i)
        Element => Mesh % Elements(j)
        j = SumPerm(j)
        if (j==0) cycle
        if (j<0) then
          j=-j
          sgn=-1
        else
          sgn=1
        end if
      else
        j = Solver % ActiveElements(i)
        Element => Mesh % Elements(j)
      end if
      l = ElemPerm(j)

      !IF(j>SIZE(ce)) PRINT *,'Too big j:',j,SIZE(ce)


      n = Element % TYPE % NumberOfNodes
      m = dofs*n

      ! Get the local stiffness matrix as saved by the primary solver
      call GetLocalMatrixStorage( PhysSolver, m, Stiff, Force, Found, &
          ElemInd = Element % ElementIndex )
      if (.not. Found) call Fatal(Caller,'Could not find local stiffness matrix!')

      ! Get the solution from stress solver
      if (dofs == 1) then
        call GetLocalSolution( local_sol,UElement=Element,USolver=PhysSolver)
      else
        call GetLocalSolution( local_sol_array,UElement=Element,USolver=PhysSolver)
        local_sol(1:m) = reshape( local_sol_array(1:dofs,1:n), [m] )
      end if

      local_act(1:m) = matmul( Stiff(1:m,1:m), local_sol(1:m) )

      ! Elemental energy assuming unity multiplier.
      if (UseAdjoint) then
        ! Get the solution from stress solver
        if (dofs == 1) then
          call GetLocalSolution( local_sol,UElement=Element,UVariable=AdjSol)
        else
          call GetLocalSolution( local_sol_array,UElement=Element,UVariable=AdjSol)
          local_sol(1:m) = reshape( local_sol_array(1:dofs,1:n), [m] )
        end if
        ! This sign leads to convergence of the bisection iteration.
        if (sgn>0) then
          spos = spos + sum( local_sol(1:m) * local_act(1:m) )
          ce(l) = ce(l) + sgn1 * sum( local_sol(1:m) * local_act(1:m) )
        else
          sneg = sneg + sum( local_sol(1:m) * local_act(1:m) )
          ce(l) = ce(l) + sgn2 * sum( local_sol(1:m) * local_act(1:m) )
        end if
      else if (NoModes > 0 ) then
        ce(l) = 0.0_dp
        do k=1,NoModes
          if (dofs == 1) then
            call GetLocalConsmode( local_sol,UElement=Element,USolver=PhysSolver,NoMode=k)
          else
            call GetLocalConsmode( local_sol_array,UElement=Element,USolver=PhysSolver,NoMode=k)
            local_sol(1:m) = reshape( local_sol_array(1:dofs,1:n), [m] )
          end if
          ce(l) = ce(l) + sum( local_sol(1:m) * local_act(1:m) )
        end do
      else
        ce(l) = sum( local_sol(1:m) * local_act(1:m) )
      end if
    end do

    ! Derivative of elemental energy
    dc = dc - penal*x**(penal-1) * efrac * ce

    ! Objective function
    obj = obj + sum( (emin + efrac*( x**penal ) ) * ce )

    if (InfoActive(20)) then
      print *,'Objective:',obj,spos,sneg
      print *,'ce:',sum(ce),sum(abs(ce)),minval(ce), maxval(ce)
      print *,'dc:',sum(dc),sum(abs(dc)),minval(dc), maxval(dc)
      print *,'x:',sum(x),sum(abs(x)),minval(x), maxval(x)
    end if

  end subroutine ObjectiveGradients


  !--------------------------------------------------------------------------
  !> Optimality criteria method (section 2.2 in paper) for maximum/minimum
  !> stiffness/compliance. Heuristic updating scheme for the element densities
  !> to find the Lagrangian multiplier.
  !--------------------------------------------------------------------------
  subroutine UpdateDensities(x,dc,dv,g)

    real(kind=dp) :: x(:), dc(:), dv(:)
    real(kind=dp) :: g

    real(kind=dp), allocatable :: xnew(:)
    real(kind=dp) :: Vi, l1, l2, lmid, move, err, tol, V0, damp
    integer :: k
    logical :: Visited = .FALSE.

    save move, tol

    l1 = 0.0_dp
    l2 = 1.0e9_dp
    ! maximum update of density

    move = ListGetCReal(Params,'Bisection search max change',Found )
    if (.not. Found) move = 0.2_dp

    tol = ListGetCReal(Params,'Bisection search tolerance',Found )
    if (.not. Found) tol = 1.0d-6

    damp = ListGetCReal(Params,'Bisection search damping exponent',Found )
    if (.not. Found) damp = 0.5_dp

    ! Desired total volume
    V0 = volFrac * sum(dv)
    V0 = ParallelReduction(V0)

    allocate(xnew(size(x)))
    xnew = 0.0_dp

    do k=1,1000
      lmid = 0.5_dp*(l2+l1)

      ! Note: xnew in [0,1]
      ! Suggested new density
      xnew = x*(max(1.0e-10,-dc/(dv*lmid)))**damp

      ! Regulators and limiters
      xnew = max(0.0_dp,max(x-move,min(1.0_dp,min(x+move,xnew))))

      ! Volume balance should become zero!
      Vi = sum(dv*xnew)
      Vi = ParallelReduction(Vi)

      if (Vi > V0) then
        l1 = lmid
      else
        l2 = lmid
      end if

      err = (l2-l1)/(l1+l2)
      if ( InfoActive(15)) then
        print *,'Bisection:',k,Vi,l1,l2,err
      end if

      if ( err < tol ) exit
    end do

    x = xnew
    g = Vi - V0
    call Info(Caller,'Number of bisection iterations: '//I2S(k),level=7)
    write(Message,'(A,2ES12.3)') 'Volume target and accuracy: ',V0,g
    call Info(Caller, Message, level=7)


  end subroutine UpdateDensities


  !------------------------------------------------------------------
  !> Applies a filter given by CRS matrix with rowsum scaled to unity.
  !> Optionally use a weight vector before and after scaling.
  !------------------------------------------------------------------
  subroutine ApplyFilter( Fmat, x, xf, niter, w, wmin, Trans )
    type(Matrix_t), pointer :: Fmat
    real(kind=dp), pointer :: x(:), xf(:)
    integer, optional :: niter
    real(kind=dp), pointer, optional :: w(:)
    real(kind=dp), optional :: wmin
    logical, optional :: Trans

    real(kind=dp), allocatable :: xtmp(:)
    real(kind=dp), pointer :: SValues(:)
    integer :: n, m, i, j
    logical :: DoTrans

    m = 1
    if (present(niter)) m = niter

    DoTrans = .false.
    if (present(Trans)) DoTrans = Trans

    allocate(xtmp(size(x)))

    if ( present(w)) then
      xtmp = x*w
      if (.not. present(wmin)) then
        call Fatal(Caller,'If we have weight we need "wmin" as well!')
      end if
      call Info(Caller,'Scaling filter with weight!',level=20)
    else
      xtmp = x
    end if

    do i=1,m
      if ( DoTrans ) then
        call TransposeMatrixVectorMultiply( Fmat, xtmp, xf)
      else
        call MatrixVectorMultiply( Fmat, xtmp, xf)
      end if
      if (i<m) then
        xtmp = xf
      end if
    end do

    if ( present(w)) then
      xf = xf/max(w,wmin)
    end if
    deallocate(xtmp)

  end subroutine ApplyFilter


  !------------------------------------------------------------------
  !> Applies a ave/min/max filter that is determined by the topology of
  !> the filter matrix but no values are used.
  !------------------------------------------------------------------
  subroutine ApplyTopologyFilter( Fmat, mode, x, xf, niter )
    type(Matrix_t), pointer :: Fmat
    integer :: mode
    real(kind=dp), pointer :: x(:), xf(:)
    integer, optional :: niter

    real(kind=dp), pointer :: xtmp(:)
    integer, pointer :: Cols(:), Rows(:)
    integer :: n, m, i, j

    m = 1
    if (present(niter)) m = niter

    n = Fmat % NumberOfRows
    Rows => Fmat % Rows
    Cols => Fmat % Cols

    if (m>1) then
      allocate(xtmp(n))
      xtmp = x
    else
      xtmp => x
    end if

    do j=1,m
      do i=1,n
        select case( mode )
        case( 0 )
          xf(i) = sum(xtmp(Cols(Rows(i):Rows(i+1)-1))) / (Rows(i+1)-Rows(i))
        case( 1 )
          xf(i) = minval(xtmp(Cols(Rows(i):Rows(i+1)-1)))
        case( 2 )
          xf(i) = maxval(xtmp(Cols(Rows(i):Rows(i+1)-1)))
        end select
      end do
      if (j<m) xtmp = xf
    end do

    if (m>1) deallocate(xtmp)

  end subroutine ApplyTopologyFilter



  !-------------------------------------------------------------------
  !> Normalize the entries such that the rowsum (or columnsum) is unity
  !-------------------------------------------------------------------
  subroutine NormalizeFilter(A,TransNorm)
    type(Matrix_t), pointer :: A
    logical :: TransNorm

    integer :: i,j,k,n
    real(kind=dp), allocatable :: colsum(:)
    real(kind=dp) :: rsum

    n = A % NumberOfRows

    ! Then create the standard projector normalized by rowsum
    call Info('NormalizeFilter','Normalizing filter by rowsum!')
    do i=1, A % NumberOfRows
      rsum = 0.0_dp
      do j = A % Rows(i), A % Rows(i+1)-1
        rsum = rsum + A % Values(j)
      end do
      do j = A % Rows(i), A % Rows(i+1)-1
        A % Values(j) = A % Values(j) / rsum
      end do
    end do

  end subroutine NormalizeFilter

  !----------------------------------------------------------------------------
  !> Create filter that inclues just closest neighbours associated attached by
  !> faces (3D) or edges (2D). This has rather small support and needs to be
  !> typically applied several times.
  !----------------------------------------------------------------------------
  function CreateSimpleFilter() result ( Emat )
    type(Matrix_t), pointer :: Emat

    type(Mesh_t), pointer :: Mesh
    type(Element_t), pointer :: Face, ElemA, ElemB
    integer :: i,j,k,k1,k2,NoElems,kcum(27)

    call Info(Caller,'Creating filter based on element-to-element connectivity')
    Mesh => Solver % Mesh
    if ( Mesh % MeshDim == 3 ) then
      call FindMeshFaces3D(Mesh)
      NoElems = Mesh % NumberOfFaces
    else
      call FindMeshEdges2D(Mesh)
      NoElems = Mesh % NumberOfEdges
    end if

    ! Create sparse matrix for element-to-element connectivity
    Emat => AllocateMatrix()
    Emat % FORMAT = MATRIX_LIST

    ! Add the max index first because list matrix likes this
    i = Mesh % NumberOfBulkElements
    call List_AddToMatrixElement( EMat % ListMatrix,i,i,0.0_dp )

    do i=1, NoElems
      if ( Mesh % MeshDim == 3 ) then
        Face => Mesh % Faces(i)
      else
        Face => Mesh % Edges(i)
      end if
      ElemA => Face % BoundaryInfo % Left
      ElemB => Face % BoundaryInfo % Right
      if (.not. associated(ElemA) .or. .not. associated(ElemB) ) cycle

      k1 = ElemA % ElementIndex
      k2 = ElemB % ElementIndex

      ! Each neighbour gets the same weight
      call List_AddToMatrixElement( Emat % ListMatrix,k1,k2,1.0_dp )
      call List_AddToMatrixElement( Emat % ListMatrix,k2,k1,1.0_dp )

      ! Set diagonals too. This way the filter has 0.5 weight for itself.
      call List_AddToMatrixElement( Emat % ListMatrix,k1,k1,1.0_dp )
      call List_AddToMatrixElement( Emat % ListMatrix,k2,k2,1.0_dp )
    end do

    ! Go from list matrix to more efficient CRS matrix
    call List_ToCRSMatrix(Emat)

    if (InfoActive(10)) then
      k1 = huge(k1)
      k2 = 0
      kcum = 0
      do i=1, Emat % NumberOfRows
        k = Emat % Rows(i+1) - Emat % Rows(i)
        k1 = min(k1, k)
        k2 = max(k2, k)
        kcum(k) = kcum(k) + 1
      end do
      do i=1,k2
        if (kcum(i)>0) print *,'Cumulative hits:',i,kcum(i)
      end do
    end if

    call Info(Caller,'Number of hits range for filter ['//I2S(k1)//','//I2S(k2)//']')
    call Info(Caller,'Number of rows in filter: '//trim(I2S(Emat % NumberOfRows)))
    call Info(Caller,'Number of non-zeros in filter: '//trim(I2S(size(Emat % Values))))

  end function CreateSimpleFilter


  !------------------------------------------------------------------------------------------------
  !> Create filter that includes all elements witing distance smaller than "rmin" between elements.
  !> We use the connectivity of simple filter to find the candidate elements.
  !------------------------------------------------------------------------------------------------
  function CreateDistanceFilter() result ( Rmat )
    type(Matrix_t), pointer :: Rmat

    type(Matrix_t), pointer :: Emat
    type(Mesh_t), pointer :: Mesh
    type(Element_t), pointer :: Element
    real(kind=dp), allocatable :: ElemCenters(:,:)
    real(kind=dp) :: rfilter, rfilter2, rik2
    integer :: NoElems,i,i2,j,k,k1,k2,k3,n,kmax,kmin
    integer, pointer :: NodeIndexes(:)
    integer, allocatable :: Inds(:),kcum(:)
    real(kind=dp), allocatable :: Dist(:)

    call Info(Caller,'Creating filter based on element-to-element distance')
    Mesh => Solver % Mesh

    call ResetTimer('DistanceFilter')

    dim = Mesh % MeshDim
    rfilter = ListGetCReal(Solver % Values,'Distance Filter Radius', UnfoundFatal = .TRUE.)
    rfilter2 = rfilter**2
    NoElems = Mesh % NumberOfBulkElements

    n = 1000
    allocate(Inds(n), Dist(n) )
    Inds = 0
    Dist = 0.0_dp

    Emat => CreateSimpleFilter()

    ! Compute center of elements for speedier distance computation.
    allocate(ElemCenters(dim,NoElems))
    do i=1,NoElems
      Element => Mesh % Elements(i)
      n = Element % TYPE % NumberOfNodes
      NodeIndexes => Element % NodeIndexes
      ElemCenters(1,i) = sum(Mesh % Nodes % x(NodeIndexes)) / n
      ElemCenters(2,i) = sum(Mesh % Nodes % y(NodeIndexes)) / n
      if (dim==3) ElemCenters(3,i) = sum(Mesh % Nodes % z(NodeIndexes)) / n
    end do

    ! Create sparse matrix for element-to-element connectivity
    Rmat => AllocateMatrix()
    Rmat % FORMAT = MATRIX_LIST

    ! Add the max index first because list matrix likes this
    call List_AddToMatrixElement( RMat % ListMatrix,NoElems,NoElems,0.0_dp )

    kmax = 0
    kmin = huge(kmin)

    do i=1, NoElems
      k1 = 1
      k3 = 1
      Inds(1) = i
      Dist(1) = rfilter
      do while(.true.)
        k2=k3
        do k=k1,k2
          do j=Emat % Rows(Inds(k)), Emat % Rows(Inds(k)+1)-1
            i2 = Emat % Cols(j)
            if (any(Inds(1:k3) == i2)) cycle

            ! square of distance between element centers
            rik2 = sum((ElemCenters(:,i)-ElemCenters(:,i2))**2)

            ! Check whether the centerpoints are within filter radius
            if ( rik2 < rfilter2 ) then
              k3 = k3 + 1
              if (k3 > size(Inds)) call Fatal(Caller,'Too many neighbours!')
              Inds(k3) = i2
              Dist(k3) = rfilter-sqrt(rik2)
            end if
          end do
        end do
        ! We found no new elements within radius
        if (k3 == k2) exit

        ! We have tested neighbours for 'k2' elements already
        k1 = k2+1
      end do

      ! Assemble the found ones in a new row.
      ! For List_Add the reverse ordering seems to be faster.
      call SortF(k3,Inds,Dist)

      do k1=1,k3
        do k=k1+1,k2
          if (Inds(k1) == Inds(k) ) call Fatal(Caller,'Duplicate indeces when creating distance filter!')
        end do
      end do

      do k1=k3,1,-1
        call List_AddToMatrixElement( RMat % ListMatrix,i,Inds(k1),Dist(k1))
      end do

      kmax = max(kmax,k3)
      kmin = min(kmin,k3)
    end do

    ! Go from list matrix to more efficient CRS matrix
    call List_ToCRSMatrix(Rmat)

    ! We do not need the element-to-element connectivity any more.
    call FreeMatrix(Emat)

    if (InfoActive(10)) then
      allocate(kcum(kmax))
      kcum = 0
      do i=1, Rmat % NumberOfRows
        k = Rmat % Rows(i+1) - Rmat % Rows(i)
        kcum(k) = kcum(k) + 1
      end do
      do i=1,kmax
        if (kcum(i)>0) print *,'Cumulative hits:',i,kcum(i)
      end do
    end if

    call Info(Caller,'Number of hits range for filter ['//I2S(kmin)//','//I2S(kmax)//']')
    call Info(Caller,'Number of rows in filter: '//trim(I2S(Rmat % NumberOfRows)))
    call Info(Caller,'Number of non-zeros in filter: '//trim(I2S(size(Rmat % Values))))
    call CheckTimer(Caller,Delete=.TRUE.)

  end function CreateDistanceFilter


!------------------------------------------------------------------------------
!> Given a topology xPhys create a 0/1 presentation that conserves volume.
!------------------------------------------------------------------------------
  subroutine DefineTopologyBW(dv0,xPhys,bw)
    real(kind=dp), pointer :: dv0(:), xPhys(:), bw(:)

    real(kind=dp) :: xlow, xup, xmid, h, q
    real(kind=dp), allocatable :: histv(:), cumv(:), tmp_histv(:)

    real(kind=dp) :: Vtot, Vtarget
    integer :: i,j,k,m,iter,ierr
    logical :: Hit

    m = 100
    allocate(histv(0:m),cumv(0:m))
    if (ParEnv % MyPe > 1 ) then
      allocate(tmp_histv(0:m))
    end if

    xlow = 0.0_dp
    xup = 1.0_dp

    Vtot = sum(dv0(1:nsize))
    Vtarget = volFrac * Vtot

    do iter=1,10
      h = (xup-xlow) / m
      histv = 0.0_dp
      cumv = 0.0_dp
      do i=1,nsize
        j = max(0,min(ceiling((xPhys(i)-xlow)/h),m))
        histv(j) = histv(j) + dv0(i)
      end do

      if ( ParEnv % PEs > 1 ) then
        tmp_histv(0:m) = histv(0:m)
        call MPI_ALLREDUCE( tmp_histv, histv, m+1, &
            MPI_DOUBLE_PRECISION, MPI_SUM, ELMER_COMM_WORLD, ierr )
      end if

      cumv(0) = histv(0)
      do i=1,m
        cumv(i) = cumv(i-1) + histv(i)
      end do

      Hit = .false.
      q = 1.0_dp

      do i=1,m
        if (cumv(i-1) < Vtarget .and. cumv(i) > Vtarget) then
          xlow = xlow + (i-1)*h
          xup = xlow + h

          q = (Vtarget-cumv(i-1))/(cumv(i)-cumv(i-1))

          Hit = .true.
          exit
        else if (abs(cumv(i-1)-Vtarget) < epsilon(h)) then
          xlow = xlow + (i-1)*h
          xup = xlow
          exit
        else if (abs(cumv(i)-Vtarget) < epsilon(h)) then
          xlow = xlow + i*h
          xup = xlow
          exit
        end if
      end do
      if (.not. Hit) exit
    end do

    xmid = (1-q)*xlow + q*xup

    where(xPhys > xmid )
      bw = 1.0_dp
    elsewhere
      bw = 0.0_dp
    end where

    write(Message,'(A,ES12.3)') 'Mass conserving B&W limit after '//I2S(iter)//' iters: ',xmid
    call Info(Caller,Message,level=7)

  end subroutine DefineTopologyBW


!------------------------------------------------------------------------------
!> Given a nodal topology xPhys find a zero levelset such that the volume
!> constraint is conserved as accurately as possible.
!------------------------------------------------------------------------------
  subroutine DefineTopologyZeroLevel(xPhysVar,bw)
    type(Variable_t), pointer :: xPhysVar
    real(kind=dp), pointer :: dv0(:), bw(:)

    real(kind=dp) :: xlow, xup, xmid, h, xAtIp, weight, detJ, f, q
    real(kind=dp), allocatable :: histv(:), cumv(:), tmp_histv(:), Basis(:)

    real(kind=dp) :: Vtot, Vtarget
    integer :: i,j,k,m,n,t,iter,elem,RelOrder,ierr
    logical :: Hit, Stat
    type(GaussIntegrationPoints_t) :: IP

    m = 1000
    allocate(histv(m+1),cumv(m+1))
    if ( ParEnv % MyPe > 1 ) then
      allocate(tmp_histv(m+1))
    end if

    xlow = 0.0_dp
    xup = 1.0_dp

    n = Mesh % MaxElementNodes
    allocate(Basis(n))
    Basis = 0.0_dp

    RelOrder = ListGetInteger( Solver % Values,'Levelset Integration Relative Order',Found)
    if (.not. Found) RelOrder = 1

    do iter=1,1 !0
      h = (xup-xlow) / m
      histv = 0.0_dp
      cumv = 0.0_dp

      do elem=1,Mesh % NumberOfBulkElements
        Element => Mesh % Elements(elem)
        n = Element % Type % NumberOfNodes

        IP = GaussPoints(Element, RelOrder=RelOrder)

        do t=1,IP % n
          stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), IP % W(t), detJ, Basis )
          xAtIp = sum(Basis(1:n) * xPhysVar % Values(xPhysVar % Perm(Element % NodeIndexes)))
          weight = detJ * IP % s(t)

          if (xAtIp <= xlow ) then
            histv(1) = histv(1) + weight
          else if ( xAtIp >= xup ) then
            histv(m+1) = histv(m+1) + weight
          else
            f = (xAtIp-xlow)/h
            j = ceiling(f)
            q = j-f
            histv(j+1) = histv(j+1) + (1-q) * weight
            histv(j) = histv(j) + q * weight
          end if
        end do
      end do

      if ( ParEnv % PEs > 1 ) then
        ! For parallel runs communicate the histogram here.
        tmp_histv = histv
        call MPI_ALLREDUCE( tmp_histv, histv, m+1, &
            MPI_DOUBLE_PRECISION, MPI_SUM, ELMER_COMM_WORLD, ierr )
      end if

      cumv(1) = histv(1)
      do i=2,m+1
        cumv(i) = cumv(i-1) + histv(i)
      end do

      Vtot = cumv(m+1)
      Vtarget = (1-volFrac) * Vtot

      if ( ParEnv % MyPe == 0) then
        ! We may optionally save the histogram and its cumulative sum.
        if ( ListGetLogical( Solver % Values,'Save Density Histogram',Found ) ) then
          block
            integer :: IoUnit
            open(newunit=IoUnit, FILE="dens_hist.dat")
            do i=1,m+1
              write(IoUnit,*) i-1, (i-1)*h, histv(i), cumv(i), cumv(i)/Vtot
            end do
            close(IoUnit)
          end block
        end if
      end if

      Hit = .false.
      q = 1.0_dp

      do i=1,m
        if (cumv(i) < Vtarget .and. cumv(i+1) > Vtarget) then
          xlow = xlow + (i-1)*h
          xup = xlow + h

          q = (Vtarget-cumv(i))/(cumv(i+1)-cumv(i))
          Hit = .true.
          exit
        else if (abs(cumv(i)-Vtarget) < epsilon(h)) then
          xlow = xlow + (i-1)*h
          xup = xlow
          exit
        else if (abs(cumv(i+1)-Vtarget) < epsilon(h)) then
          xlow = xlow + i*h
          xup = xlow
          exit
        end if
      end do
      if (.not. Hit) exit
    end do

    ! This is the new approximation of the mid value that gives the desider volume within (x>xmid).
    xmid = (1-q)*xlow + q*xup

    if ( ListGetLogical( Solver % Values,'Levelset Symmmetric',Found ) ) then
      ! Define levelset as simple offset from nodal density.
      bw = xPhysVar % Values-xmid
    else
      ! Map levelset between [-1,1] such that zero levelset is at desired value.
      where(xPhysVar % Values > xmid )
        bw = (xPhysVar % Values-xmid)/(1.0_dp-xmid)
      elsewhere
        bw = (xPhysVar % Values-xmid)/xmid
      end where
    end if

    write(Message,'(A,ES12.3)') 'Mass conserving zero levelset: ',xmid
    call Info(Caller,Message,level=7)

  end subroutine DefineTopologyZeroLevel



!------------------------------------------------------------------------------
!> Assembly of the matrix equation used for PDE filtering.
!> We may assembly both matrix and r.h.s., or just the r.h.s.
!------------------------------------------------------------------------------
  subroutine LocalMatrix( Element, n, DoMatrix, Diff, x )
!------------------------------------------------------------------------------
    integer :: n
    type(Element_t), target :: Element
    logical :: DoMatrix
    real(kind=dp), pointer :: x(:)
    real(kind=dp) :: Diff(3)
!------------------------------------------------------------------------------
    real(kind=dp) :: Weight, D, xi
    real(kind=dp) :: Basis(n),dBasisdx(n,3),DetJ,LoadAtIP
    real(kind=dp) :: STIFF(n,n), FORCE(n)
    logical :: Stat,Found
    integer :: i,t,p,q
    type(GaussIntegrationPoints_t) :: IP
    type(Nodes_t) :: Nodes
    save Nodes
!------------------------------------------------------------------------------

    D = 1.0_dp

    call GetElementNodes( Nodes )

    ! Separate matrix and force vector integration because we may use lower order
    ! integration scheme for the force vector.
    !----------------------------------------------------------------------------
    if (DoMatrix ) then
      STIFF = 0._dp
      IP = GaussPoints( Element )
      do t=1,IP % n
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis, dBasisdx )

        Weight = IP % s(t) * DetJ
        if (Csymmetry) weight = Weight * sum(Basis(1:n) * Nodes % x(1:n))

        do p=1,n
          do q=1,n
            STIFF(p,q) = STIFF(p,q) + Weight * Basis(p) * Basis(q)
          end do
        end do

        do p=1,n
          if (SkipInterface) then
            if (InterfaceNode(Element % NodeIndexes(p))) cycle
          end if
          do q=1,n
            STIFF(p,q) = STIFF(p,q) + Weight * &
                sum( Diff(1:dim) * dBasisdx(p,1:dim) * dBasisdx(q,1:dim) )
          end do
        end do
      end do
    end if

    FORCE = 0._dp

    if ( Csymmetry ) then
      ! We don't integrate accurately area with cylindrical coordinates.
      IP = GaussPoints( Element )
    else
      IP = GaussPoints( Element, np=1 )
    end if

    xi = x(ElemPerm(Element % ElementIndex))

    do t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis )

      Weight = IP % s(t) * DetJ
      if (Csymmetry) weight = Weight * sum(Basis(1:n) * Nodes % x(1:n))

      FORCE(1:n) = FORCE(1:n) + Weight * Basis(1:n) * xi
    end do

    if ( DoMatrix ) then
      call DefaultUpdateEquations(STIFF,FORCE)
    else
      call DefaultUpdateForce(FORCE)
    end if

!------------------------------------------------------------------------------
  end subroutine LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Solves a diffusion-reaction equation to smooth down given values "x" to "xf".
!> PDE based filtering is ideal since it can use the parallel machinery of Elmer.
!------------------------------------------------------------------------------
  subroutine ApplyPDEFilter(x, xf, Diff, w, wmin )
!------------------------------------------------------------------------------

    real(kind=dp), pointer :: x(:), xf(:)
    real(kind=dp) :: Diff(3)
    real(kind=dp), pointer, optional :: w(:)
    real(kind=dp), optional :: wmin

    real(kind=dp), pointer :: xtmp(:)
    integer :: n, t, active
    real(kind=dp) :: Norm
    type(Variable_t), pointer :: pVar

    ! Create weighted elemental field if requested.
    if ( present(w)) then
      allocate(xtmp(size(x)))
      xtmp = x*w
      if (.not. present(wmin)) then
        call Fatal(Caller,'If we have weight we need "wmin" as well!')
      end if
      call Info(Caller,'Scaling filter with weight!',level=20)
    else
      xtmp => x
    end if

    ! Assembly the matrix equation at the 1st time.
    ! Later just define the r.h.s. vector.
    if ( acc_domatrix ) then
      call DefaultInitialize()
    else
      Solver % Matrix % Rhs = 0.0_dp
    end if

    Active = GetNOFActive()
    do t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes()
      call LocalMatrix(  Element, n, acc_domatrix, Diff, xtmp )
    end do

    call DefaultDirichletBCs()

    Norm = DefaultSolve()
    pVar => Solver % Variable

    ! After solving the nodal values we need to transfer them back to elemental values.
    do t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes()
      xf(t) = sum(pVar % Values(pVar % Perm(Element % NodeIndexes)))/n
    end do

    ! If weigting was used revert back.
    if ( present(w)) then
      deallocate(xtmp)
      xf = xf/max(w,wmin)
    end if

    ! We have done the matrix. Freeze it and never touch it again.
    acc_domatrix = .false.

  end subroutine ApplyPDEFilter
#if 0
!------------------------------------------------------------------------------
!> Solves a adjoint problem of the primary problem with different rhs.
!------------------------------------------------------------------------------
  subroutine SolveAdjointProblem(x)
!------------------------------------------------------------------------------

    real(kind=dp), pointer :: x(:)

    type(Solver_t), pointer :: pSolver
    real(kind=dp), pointer :: pRhs(:)


    integer :: n, t, active
    logical :: DoMatrix = .TRUE.
    real(kind=dp) :: Norm
    real(kind=dp), pointer :: aRhs(:)


    pSolver => Model % Solver
    pRhs => PhysSolver % Matrix % rhs

    Model % Solver => PhysSolver
    PhysSolver % Variable => aVar

    if ( DoMatrix ) then
      allocate( aRhs(size(aVar % Values)) )
      PhysSolver % Matrix % Rhs => aRhs
      aRhs = 0.0_dp
      Active = GetNOFActive()
      do t=1,Mesh % NumberOfBoundaryElements
        Element => Mesh % Elements(t)
        n  = GetElementNOFNodes()
        !      call LocalMatrix(  Element, n, DoMatrix, Diff, xtmp )
      end do
      ! We have done the rhs. No need to redo.
      DoMatrix = .FALSE.
    else
      PhysSolver % Matrix % Rhs => aRhs
    end if

    ! Solve the adjoint problem with the same matrix equation is the primary matrix.
    call ListAddLogical( PhysSolver % Values,'Skip Compute Change',.TRUE.)
    Norm = DefaultSolve()
    call ListAddLogical( PhysSolver % Values,'Skip Compute Change',.FALSE.)

    ! Revert the saved pointers back
    Model % Solver => pSolver
    PhysSolver % Matrix % Rhs => pRhs
    PhysSolver % Variable => uVar

  end subroutine SolveAdjointProblem
#endif


end module toa_acc
