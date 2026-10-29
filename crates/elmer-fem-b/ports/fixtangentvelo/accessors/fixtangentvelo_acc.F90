! T022-A FixTangentVelo port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/FixTangentVelo.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments. All Fortran-object
! traffic — TYPE(...) handles, OPTIONAL/POINTER dummies, hidden lengths,
! whole-array updates, order-sensitive reductions — stays inside Fortran,
! so the Rust side is pure C calls. The wrapper module mirrors
! FixTangentVelo's own SAVE statement (Visited, Nodes, FixVar, VeloVar,
! FixVals, FixPerm, ActiveCoord, DownstreamOnly) with identical
! lifetimes; entry-scope locals (Mesh, Params, Element, BC, loop scalars)
! are re-fetched per call, which resolves to the same targets the
! original's stale pointers still address. All floating-point work (edge
! vectors, the 3x3 solve, the MAX/MIN accumulation, ComputeNorm) and the
! Info/Fatal reports stay in Fortran (same compiler), so bit-identity
! holds by construction wherever the call sequence matches. The
! ApplyTangentFix internal procedure is kept verbatim as a module
! procedure and invoked with `.FALSE.`, exactly as the original calls it.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module fixtangentvelo_acc
  use iso_c_binding
  use defutils
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  ! Per-call mesh handle (the original's `Mesh => Solver % Mesh`).
  type(mesh_t), pointer, save :: acc_mesh => null()
  ! SAVE state mirroring the original's SAVE statement.
  logical, save :: acc_visited = .false.
  type(nodes_t), target, save :: acc_nodes
  type(variable_t), pointer, save :: acc_fixvar => null()
  type(variable_t), pointer, save :: acc_velovar => null()
  real(kind=dp), pointer, save :: acc_fixvals(:) => null()
  integer, pointer, save :: acc_fixperm(:) => null()
  integer, save :: acc_activecoord = 0
  logical, save :: acc_downstreamonly = .false.

contains

  ! -- entry state ----------------------------------------------------
  subroutine ftv_begin(h_solver) bind(c, name="ftv_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
    acc_mesh => acc_solver % mesh
    call Info('FixTangentVelo','Regulating surface velocity not to go out domain')
  end subroutine ftv_begin

  function ftv_need_first() result(r) bind(c, name="ftv_need_first")
    integer(c_int) :: r
    r = merge(1, 0, .not. acc_visited)
  end function ftv_need_first

  ! -- first-visit block (original IF(.NOT. Visited) body, verbatim) ---
  subroutine ftv_first() bind(c, name="ftv_first")
    type(valuelist_t), pointer :: Params
    character(len=MAX_NAME_LEN) :: str
    type(element_t), pointer :: Element
    type(valuelist_t), pointer :: BC
    integer :: i, t, n, dim
    logical :: Found

    Params => acc_solver % Values

    dim = CoordinateSystemDimension()
    if (dim /= 3) call Fatal('FixTangentVelo','Implemented only for 3D systems!')

    acc_activecoord = ListGetInteger(Params,'Active Coordinate',Found)
    if (.not. Found) acc_activecoord = dim

    acc_downstreamonly = ListGetLogical(Params,'Fix Downstream Only',Found)

    str = ListGetString(Params,'Velocity Variable Name',Found)
    if (.not. Found) str = 'flow solution'
    acc_velovar => VariableGet(acc_mesh % Variables, str)
    if (.not. associated(acc_velovar)) then
      call Fatal('FixTangentVelo','Velocity field variable does not exist: '//trim(str))
    end if

    allocate(acc_fixperm(acc_mesh % NumberOfNodes))
    acc_fixperm = 0
    do t = acc_mesh % NumberOfBulkElements+1, &
        acc_mesh % NumberOfBulkElements + acc_mesh % NumberOfBoundaryElements
      Element => acc_mesh % Elements(t)
      BC => GetBC(Element)
      if (.not. associated(BC)) cycle
      if (ListGetLogical(BC,'Tangent Surface',Found)) then
        acc_fixperm(Element % NodeIndexes) = 1
      end if
    end do

    n = 0
    do i = 1, acc_mesh % NumberOfNodes
      if (acc_fixperm(i) > 0) then
        n = n + 1
        acc_fixperm(i) = n
      end if
    end do

    n = ParallelReduction(n)
    call Info('FixTangentVelo','Total number of nodes to fix: '//I2S(n))
    if (n == 0) call Fatal('FixTangentVelo','No nodes to fix!')

    call VariableAddVector(acc_mesh % Variables, acc_mesh, acc_solver, &
        'FixVelo', Perm = acc_fixperm)
    acc_fixvar => VariableGet(acc_mesh % Variables,'FixVelo')
    acc_fixvals => acc_fixvar % Values
    acc_visited = .true.
  end subroutine ftv_first

  ! -- per-invocation reset (original `FixVals = 0.0_dp`) ---------------
  subroutine ftv_zero() bind(c, name="ftv_zero")
    acc_fixvals = 0.0_dp
  end subroutine ftv_zero

  function ftv_num_bulk() result(n) bind(c, name="ftv_num_bulk")
    integer(c_int) :: n
    n = acc_mesh % NumberOfBulkElements
  end function ftv_num_bulk

  function ftv_num_bnd() result(n) bind(c, name="ftv_num_bnd")
    integer(c_int) :: n
    n = acc_mesh % NumberOfBoundaryElements
  end function ftv_num_bnd

  ! -- boundary-element body (original loop body, verbatim) -------------
  subroutine ftv_elem(t) bind(c, name="ftv_elem")
    integer(c_int), value :: t
    type(element_t), pointer :: Element
    type(valuelist_t), pointer :: BC
    integer :: i, j, n, i1, i2
    real(kind=dp) :: e1(3), e2(3), Amat(3,3), velem(4,4), v(3), c(3), Normal(3)
    logical :: Top, Found

    Element => acc_mesh % Elements(t)
    BC => GetBC(Element)
    if (.not. associated(BC)) return
    Top = ListGetLogical(BC,'Tangent Surface',Found)
    if (.not. Top) return

    call GetVectorLocalSolution(velem,UElement=Element,UVariable=acc_velovar)
    call GetElementNodes(acc_nodes, Element)
    Normal = NormalVector(Element, acc_nodes, Check=.true.)

    n = Element % TYPE % NumberOfNodes
    if (.not. (n == 3 .or. n == 4)) then
      call Fatal('FixTangentVelo','Only implemented for 3 and 4 nodes!')
    end if

    ! Go through each element node. So lazy coding here fetching the same stuff many times.
    do i = 1, n
      v = velem(1:3,i)
      i1 = modulo(i,n)+1
      i2 = modulo(i-2,n)+1

      ! Edges point to two other nodes
      e1(1) = acc_nodes % x(i1) - acc_nodes % x(i)
      e1(2) = acc_nodes % y(i1) - acc_nodes % y(i)
      e1(3) = acc_nodes % z(i1) - acc_nodes % z(i)

      e2(1) = acc_nodes % x(i2) - acc_nodes % x(i)
      e2(2) = acc_nodes % y(i2) - acc_nodes % y(i)
      e2(3) = acc_nodes % z(i2) - acc_nodes % z(i)

      ! The flow needs to be towards either node that we can have an issue.
      ! It is enough that either nodes gives postive dot product.
      if (acc_downstreamonly) then
        if (sum(v*e1) < 0 .and. sum(v*e2) < 0) cycle
      end if

      Amat = 0.0_dp
      Amat(:,1) = e1
      Amat(:,2) = e2
      Amat(acc_activecoord,3) = 1.0_dp

      call SolveLinSys3x3(Amat, c, v)

      j = acc_fixperm(Element % NodeIndexes(i))
      if (Normal(acc_activecoord) > 0.0_dp) then
        acc_fixvals(j) = max(acc_fixvals(j),c(3))
      else
        acc_fixvals(j) = min(acc_fixvals(j),c(3))
      end if
    end do
  end subroutine ftv_elem

  ! -- diagnostic range report (original InfoActive(20) guard) ----------
  subroutine ftv_range() bind(c, name="ftv_range")
    if (InfoActive(20)) then
      call VectorValuesRange(acc_fixvals, size(acc_fixvals),'FixVals')
    end if
  end subroutine ftv_range

  ! -- `Fix Tangent Velocity` keyword read ------------------------------
  function ftv_fix_flag() result(r) bind(c, name="ftv_fix_flag")
    integer(c_int) :: r
    type(valuelist_t), pointer :: Params
    logical :: Found
    Params => acc_solver % Values
    if (ListGetLogical(Params,'Fix Tangent Velocity',Found)) then
      r = 1
    else
      r = 0
    end if
  end function ftv_fix_flag

  ! -- ApplyTangentFix internal procedure (verbatim, module-retained) ---
  subroutine ApplyTangentFix(UndoFix)
    logical :: UndoFix
    integer :: sgn,i,jv,jf,vdofs

    call Info('FixTangentVelo','Applying the fix to current velocity field')

    sgn = -1
    if (UndoFix) sgn = 1

    vdofs = acc_velovar % dofs
    do i = 1, acc_mesh % NumberOfNodes
      jf = acc_fixvar % Perm(i)
      jv = acc_velovar % Perm(i)
      if (jf == 0 .or. jv == 0) cycle
      acc_velovar % Values(vdofs*(jv-1)+acc_activecoord) = &
          acc_velovar % Values(vdofs*(jv-1)+acc_activecoord) + &
          sgn * acc_fixvar % Values(jf)
    end do
  end subroutine ApplyTangentFix

  subroutine ftv_apply() bind(c, name="ftv_apply")
    call ApplyTangentFix(.false.)
  end subroutine ftv_apply

  ! -- norm tail (original ComputeNorm + Variable assignment) -----------
  function ftv_tail() result(nrm) bind(c, name="ftv_tail")
    real(c_double) :: nrm
    nrm = ComputeNorm(acc_solver, size(acc_fixvals), acc_fixvals)
    acc_solver % Variable % Values = nrm
    acc_solver % Variable % Norm = nrm
  end function ftv_tail

  ! -- FixTangentVelo_Init body ------------------------------------------
  subroutine ftv_init_begin(h_solver) bind(c, name="ftv_init_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine ftv_init_begin

  subroutine ftv_init_body() bind(c, name="ftv_init_body")
    call ListAddNewLogical(acc_solver % Values,'No Matrix',.true.)
    call ListAddNewString(acc_solver % Values,'Variable', &
        '-nooutput -global FixTangent_var')
  end subroutine ftv_init_body

end module fixtangentvelo_acc
