! T021-B PoissonBEM port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/PoissonBEM.F90, compiled against the built
! .mod files of the target configure (solver-entry gate shape, ADR-057
! §3; imitates the T020 Poisson pilot, whose accessor this mirrors).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls.
!
! Layout (original line numbers refer to PoissonBEM.F90):
! - `module globmat` (orig 40-45): verbatim mirror of the plugin's own
!   module. It owns the `Matrix(BoundaryNodes,BoundaryNodes)` BSS
!   descriptor, so the port's dynamic surface keeps the reference's
!   `B __globmat_MOD_matrix` symbol with identical meaning; the Rust
!   `matvec_`/`precond_` entries operate on it through the `pbb_matvec`
!   / `pbb_precond` wrappers below (same loop order, same compiler).
! - module SAVE state mirrors the host procedure's SAVE list (orig 96):
!   Load/ElementNodes/VolumeSource/VolumeForce/AllocationsDone/
!   PotentialKnown/PotentialPerm/BoundaryPerm/BoundaryNodes/Pot/Flx.
!   Per-call reassociated pointers (Potential, ForceVector) are
!   re-pointed in `pbb_begin` on every entry, exactly like the original.
! - Phase drivers (`pbb_*`) transcribe the host statement blocks
!   verbatim (same Elmer calls, same order-sensitive whole-array and
!   vector-subscript updates), so bit-identity holds by construction.
! - Internal procedures (orig 359-727: IntegrateSource,
!   IntegrateMatrix, ComputePotential, SolveFull, FullIterSolver) are
!   kept as module procedures `pbk_*` with identical bodies; host
!   association is replaced by the module SAVE state above. The only
!   deliberate deltas: `Solver`/`SolverParam` (the same object in the
!   original) is the `acc_solver` handle, and the `ALL(x == 0.0)` probe
!   reads the `acc_pot(1:acc_bn)` section the original's size-`N` dummy
!   `x` saw.
! - `USE defutils` transitively provides the `IterCall` generic
!   (via SolverBasics -> LoadMod), exactly as in the original build;
!   `pbb_*` state pointers that cross the explicit CONTIG dummies of
!   `IterCallFTNR` carry the CONTIGUOUS attribute.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module globmat
  use types
  implicit none
  real(kind=dp), allocatable :: Matrix(:,:)
end module globmat

module poissonbem_acc
  use iso_c_binding
  use defutils
  use types
  use globmat
  use huti_sfe
  use messages
  implicit none

#include "huti_fdefs.h"

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(element_t), pointer, save :: acc_elem => null()
  ! SAVE mirror of the host's ElementNodes (consumed before the next
  ! fill, exactly as in the original).
  type(nodes_t), target, save :: acc_nodes
  integer, save :: acc_bn = 0
  logical, save :: acc_alloc = .false.
  real(kind=dp), save :: acc_at = 0.0_dp, acc_st = 0.0_dp
  ! SAVE mirrors, all sized on the first call (first-mesh sizing, never
  ! resized — identical lifetimes to the original's SAVE state).
  real(kind=dp), pointer, contiguous, save :: acc_pot(:) => null()
  real(kind=dp), pointer, contiguous, save :: acc_force(:) => null()
  real(kind=dp), pointer, save :: acc_volsrc(:) => null()
  real(kind=dp), pointer, save :: acc_volforce(:) => null()
  real(kind=dp), pointer, save :: acc_load(:) => null()
  real(kind=dp), pointer, save :: acc_flx(:) => null()
  real(kind=dp), pointer, save :: acc_potp(:) => null()
  real(kind=dp), pointer, save :: acc_diag(:) => null()
  integer, pointer, save :: acc_potperm(:) => null()
  integer, pointer, save :: acc_bndperm(:) => null()
  logical, allocatable, save :: acc_potknown(:)

  ! Matvec/Precond are Rust-defined entries (same `.so`); the iterative
  ! solver below takes their addresses exactly like the original.
  external :: Matvec, Precond

contains

  ! -- entry state (orig 105-114) -------------------------------------
  subroutine pbb_begin(h_model, h_solver) bind(c, name="pbb_begin")
    type(c_ptr), value :: h_model, h_solver
    character(len=MAX_NAME_LEN) :: equation
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    if (CurrentCoordinateSystem() /= Cartesian) then
      call Fatal('PoissonBEMSolver', &
        'This solver is implemented only for cartesian coordinates!')
    end if
    acc_pot => acc_solver % variable % values
    acc_force => acc_solver % matrix % rhs
    ! Read once, like the original; the name is never used afterwards.
    equation = ListGetString(acc_model % solver % values, 'Equation')
  end subroutine pbb_begin

  ! -- first-call sizing (orig 116-170) --------------------------------
  subroutine pbb_allocate() bind(c, name="pbb_allocate")
    integer :: j, k, n, t, istat, nn
    if (acc_alloc) return
    nn = acc_solver % mesh % numberofnodes
    allocate (acc_potperm(nn), acc_bndperm(nn), stat=istat)
    if (istat /= 0) then
      call Fatal('PoissonBEMSolver', 'Memory allocation error 1.')
    end if
    acc_potperm = 0
    acc_bndperm = 0
    acc_bn = 0
    do t = 1, acc_solver % numberofactiveelements
      acc_elem => acc_solver % mesh % elements(acc_solver % activeelements(t))
      if (acc_elem % type % elementcode == 101) cycle
      if (.not. associated(acc_elem % boundaryinfo)) cycle
      do j = 1, acc_elem % type % numberofnodes
        k = acc_elem % nodeindexes(j)
        if (acc_potperm(k) == 0) then
          acc_bn = acc_bn + 1
          acc_bndperm(acc_bn) = k
          acc_potperm(k) = acc_bn
        end if
      end do
    end do
    n = acc_model % maxelementnodes
    allocate (acc_nodes % x(n), acc_nodes % y(n), acc_nodes % z(n), &
      acc_volsrc(n), acc_flx(acc_bn), acc_potp(acc_bn), &
      acc_load(acc_bn), acc_diag(acc_bn), acc_potknown(acc_bn), &
      Matrix(acc_bn, acc_bn), acc_volforce(nn), stat=istat)
    if (istat /= 0) then
      call Fatal('PoissonBEMSolver', 'Memory allocation error 2.')
    end if
    acc_alloc = .true.
  end subroutine pbb_allocate

  ! -- zero fills (orig 175-181) ---------------------------------------
  subroutine pbb_zero() bind(c, name="pbb_zero")
    acc_at = CPUTime()
    Matrix = 0.0_dp
    acc_load = 0.0_dp
    acc_diag = 0.0_dp
    acc_force = 0.0_dp
  end subroutine pbb_zero

  ! -- boundary-condition pass (orig 190-221) --------------------------
  subroutine pbb_bc_pass() bind(c, name="pbb_bc_pass")
    integer :: i, t, n
    integer, pointer :: nodeindexes(:)
    logical :: gotit
    do t = 1, acc_solver % numberofactiveelements
      acc_elem => acc_solver % mesh % elements(acc_solver % activeelements(t))
      if (.not. associated(acc_elem % boundaryinfo)) cycle
      if (acc_elem % type % elementcode == 101) cycle
      n = acc_elem % type % numberofnodes
      nodeindexes => acc_elem % nodeindexes
      do i = 1, acc_model % numberofbcs
        if (acc_elem % boundaryinfo % constraint /= acc_model % bcs(i) % tag) cycle
        acc_load(acc_potperm(nodeindexes)) = ListGetReal(acc_model % bcs(i) % values, &
          'Potential', n, nodeindexes, gotit)
        if (.not. gotit) then
          acc_load(acc_potperm(nodeindexes)) = ListGetReal(acc_model % bcs(i) % values, &
            acc_solver % variable % name(1:acc_solver % variable % namelen), &
            n, nodeindexes, gotit)
        end if
        if (.not. gotit) then
          acc_potknown(acc_potperm(nodeindexes)) = .false.
          acc_load(acc_potperm(nodeindexes)) = ListGetReal(acc_model % bcs(i) % values, &
            'Flux', n, nodeindexes, gotit)
        else
          acc_potknown(acc_potperm(nodeindexes)) = .true.
        end if
        exit
      end do
    end do
  end subroutine pbb_bc_pass

  ! -- volume-source pass (orig 228-257) -------------------------------
  subroutine pbb_volume_pass() bind(c, name="pbb_volume_pass")
    integer :: i, k, n
    integer, pointer :: nodeindexes(:)
    logical :: gotit
    acc_volforce = 0.0_dp
    do i = 1, acc_solver % numberofactiveelements
      acc_elem => acc_solver % mesh % elements(acc_solver % activeelements(i))
      if (associated(acc_elem % boundaryinfo)) cycle
      n = acc_elem % type % numberofnodes
      nodeindexes => acc_elem % nodeindexes
      k = ListGetInteger(acc_model % bodies(acc_elem % bodyid) % values, &
        'Body Force', gotit)
      if (.not. gotit) exit
      acc_volsrc(1:n) = ListGetReal(acc_model % bodyforces(k) % values, &
        'Source', n, nodeindexes, gotit)
      if (.not. gotit) exit
      acc_nodes % x(1:n) = acc_solver % mesh % nodes % x(nodeindexes)
      acc_nodes % y(1:n) = acc_solver % mesh % nodes % y(nodeindexes)
      acc_nodes % z(1:n) = acc_solver % mesh % nodes % z(nodeindexes)
      call pbk_source(acc_volsrc, acc_volforce, acc_elem, n, acc_nodes)
    end do
    do i = 1, acc_bn
      acc_force(i) = acc_force(i) + acc_volforce(acc_bndperm(i))
    end do
  end subroutine pbb_volume_pass

  ! -- matrix assembly pass (orig 262-276) -----------------------------
  subroutine pbb_assemble_pass() bind(c, name="pbb_assemble_pass")
    integer :: n, t
    integer, pointer :: nodeindexes(:)
    do t = 1, acc_solver % numberofactiveelements
      acc_elem => acc_solver % mesh % elements(acc_solver % activeelements(t))
      if (.not. associated(acc_elem % boundaryinfo)) cycle
      if (acc_elem % type % elementcode == 101) cycle
      n = acc_elem % type % numberofnodes
      nodeindexes => acc_elem % nodeindexes
      acc_nodes % x(1:n) = acc_solver % mesh % nodes % x(nodeindexes)
      acc_nodes % y(1:n) = acc_solver % mesh % nodes % y(nodeindexes)
      acc_nodes % z(1:n) = acc_solver % mesh % nodes % z(nodeindexes)
      call pbk_matrix(Matrix, acc_force, acc_load, acc_potknown, acc_elem, n, acc_nodes)
    end do
  end subroutine pbb_assemble_pass

  ! -- diagonal fix (orig 278-284) -------------------------------------
  subroutine pbb_diagonal_fix() bind(c, name="pbb_diagonal_fix")
    integer :: i
    do i = 1, acc_bn
      if (acc_potknown(i)) then
        acc_force(i) = acc_force(i) - acc_load(i)*acc_diag(i)
      else
        Matrix(i, i) = acc_diag(i)
      end if
    end do
  end subroutine pbb_diagonal_fix

  ! -- assembly timing print (orig 287-288) -----------------------------
  subroutine pbb_print_assembly() bind(c, name="pbb_print_assembly")
    acc_at = CPUTime() - acc_at
    print *, 'Assembly (s): ', acc_at
  end subroutine pbb_print_assembly

  ! -- system solve (orig 292-297) --------------------------------------
  subroutine pbb_solve() bind(c, name="pbb_solve")
    acc_st = CPUTime()
    call pbk_solvefull()
  end subroutine pbb_solve

  ! -- boundary extraction (orig 301-327) --------------------------------
  subroutine pbb_extract() bind(c, name="pbb_extract")
    type(variable_t), pointer :: var
    integer :: i, k
    do i = 1, acc_bn
      if (acc_potknown(i)) then
        acc_flx(i) = acc_pot(i)
        acc_potp(i) = acc_load(i)
      else
        acc_flx(i) = acc_load(i)
        acc_potp(i) = acc_pot(i)
      end if
    end do
    acc_pot = acc_volforce
    do i = 1, acc_bn
      k = acc_solver % variable % perm(acc_bndperm(i))
      acc_pot(k) = acc_potp(i)
    end do
    var => VariableGet(acc_solver % mesh % variables, 'Flux')
    if (associated(var)) then
      var % values = 0.0_dp
      do i = 1, acc_bn
        k = var % perm(acc_bndperm(i))
        var % values(k) = acc_flx(i)
      end do
    end if
  end subroutine pbb_extract

  ! -- volume-computing switch (orig 329) --------------------------------
  function pbb_skip_volume() result(r) bind(c, name="pbb_skip_volume")
    integer(c_int) :: r
    logical :: found
    r = merge(1, 0, GetLogical(acc_solver % values, 'Skip Volume Computing', found))
  end function pbb_skip_volume

  ! -- interior potential pass (orig 331-345) ----------------------------
  subroutine pbb_volume_compute() bind(c, name="pbb_volume_compute")
    integer :: n, t
    integer, pointer :: nodeindexes(:)
    do t = 1, acc_solver % numberofactiveelements
      acc_elem => acc_solver % mesh % elements(acc_solver % activeelements(t))
      if (.not. associated(acc_elem % boundaryinfo)) cycle
      n = acc_elem % type % numberofnodes
      nodeindexes => acc_elem % nodeindexes
      acc_nodes % x(1:n) = acc_solver % mesh % nodes % x(nodeindexes)
      acc_nodes % y(1:n) = acc_solver % mesh % nodes % y(nodeindexes)
      acc_nodes % z(1:n) = acc_solver % mesh % nodes % z(nodeindexes)
      call pbk_potential(acc_pot, acc_potp, acc_flx, acc_elem, n, acc_nodes)
    end do
  end subroutine pbb_volume_compute

  ! -- finalize: norm, invalidate, timing print (orig 350-356) ------------
  subroutine pbb_finalize() bind(c, name="pbb_finalize")
    acc_solver % variable % norm = &
      sqrt(sum(acc_pot**2)/count(acc_solver % variable % perm > 0))
    call InvalidateVariable(acc_model % meshes, acc_solver % mesh, &
      acc_solver % variable % name)
    acc_st = CPUTime() - acc_st
    print *, 'Solve (s):    ', acc_st
  end subroutine pbb_finalize

  ! -- matvec/precond bodies (called from the Rust entries) ---------------
  subroutine pbb_matvec(n, u, v) bind(c, name="pbb_matvec")
    integer(c_int), value :: n
    real(c_double) :: u(n), v(n)
    integer :: i, j
    real(kind=dp) :: s
    do i = 1, n
      s = 0.0_dp
      do j = 1, n
        s = s + Matrix(i, j)*u(j)
      end do
      v(i) = s
    end do
  end subroutine pbb_matvec

  subroutine pbb_precond(n, u, v) bind(c, name="pbb_precond")
    integer(c_int), value :: n
    real(c_double) :: u(n), v(n)
    integer :: i
    do i = 1, n
      u(i) = v(i)/Matrix(i, i)
    end do
  end subroutine pbb_precond

  ! -- internal procedures (orig 362-419) ----------------------------------
  subroutine pbk_source(Source, Force, Element, n, Nodes)
    integer :: n
    real(kind=dp) :: Source(n), Force(:)
    type(Nodes_t) :: Nodes
    type(Element_t), target :: Element
    real(kind=dp) :: Basis(n), dBasisdx(n, 3)
    logical :: Stat
    real(kind=dp) :: SqrtElementMetric, U, V, W, S, LX, LY, LZ, x, y, z, R
    integer :: i, p, t, dim
    type(GaussIntegrationPoints_t) :: IntegStuff
    dim = CoordinateSystemDimension()
    IntegStuff = GaussPoints(Element)
    do t = 1, IntegStuff % n
      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)
      S = IntegStuff % s(t)
      Stat = ElementInfo(Element, Nodes, U, V, W, SqrtElementMetric, &
        Basis, dBasisdx)
      S = S*SqrtElementMetric
      LX = sum(Nodes % x(1:n)*Basis)
      LY = sum(Nodes % y(1:n)*Basis)
      LZ = sum(Nodes % z(1:n)*Basis)
      do p = 1, acc_solver % mesh % numberofnodes
        x = LX - acc_solver % mesh % nodes % x(p)
        y = LY - acc_solver % mesh % nodes % y(p)
        z = LZ - acc_solver % mesh % nodes % z(p)
        R = sqrt(x**2 + y**2 + z**2)
        select case (dim)
        case (2)
          W = -log(R)/(2*PI)
        case (3)
          W = (1/R)/(4*PI)
        end select
        do i = 1, n
          Force(p) = Force(p) + S*Source(i)*Basis(i)*W
        end do
      end do
    end do
  end subroutine pbk_source

  ! -- internal procedures (orig 424-526) ----------------------------------
  subroutine pbk_matrix(StiffMatrix, Force, Load, PotentialKnown, Element, n, Nodes)
    real(kind=dp) :: StiffMatrix(:, :), Force(:), Load(:)
    integer :: n
    logical :: PotentialKnown(:)
    type(Nodes_t) :: Nodes
    type(Element_t), target :: Element
    real(kind=dp) :: Basis(n), dBasisdx(n, 3)
    real(kind=dp) :: R, LX, LY, LZ, x, y, z, dWdN
    logical :: Stat, CheckNormals
    real(kind=dp) :: SqrtElementMetric, U, V, W, S, gradW(3), Normal(3)
    integer :: i, p, q, t, dim
    type(GaussIntegrationPoints_t) :: IntegStuff
    dim = CoordinateSystemDimension()
    select case (Element % type % elementcode/100)
    case (2)
      IntegStuff = GaussPoints(Element, 4)
    case (3)
      IntegStuff = GaussPoints(Element, 6)
    case (4)
      IntegStuff = GaussPoints(Element, 16)
    end select
    CheckNormals = associated(Element % boundaryinfo)
    if (CheckNormals) then
      CheckNormals = associated(Element % boundaryinfo % left) .or. &
        associated(Element % boundaryinfo % right)
    end if
    do t = 1, IntegStuff % n
      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)
      S = IntegStuff % s(t)
      Stat = ElementInfo(Element, Nodes, U, V, W, SqrtElementMetric, &
        Basis, dBasisdx)
      S = S*SqrtElementMetric
      Normal = NormalVector(Element, Nodes, U, V, CheckNormals)
      LX = sum(Nodes % x(1:n)*Basis)
      LY = sum(Nodes % y(1:n)*Basis)
      LZ = sum(Nodes % z(1:n)*Basis)
      do p = 1, acc_bn
        x = LX - acc_solver % mesh % nodes % x(acc_bndperm(p))
        y = LY - acc_solver % mesh % nodes % y(acc_bndperm(p))
        z = LZ - acc_solver % mesh % nodes % z(acc_bndperm(p))
        R = sqrt(x**2 + y**2 + z**2)
        select case (dim)
        case (2)
          W = -log(R)/(2*PI)
          GradW(1) = -x/(2*PI*R**2)
          GradW(2) = -y/(2*PI*R**2)
          GradW(3) = -z/(2*PI*R**2)
        case (3)
          W = (1/R)/(4*PI)
          GradW(1) = -x/(4*PI*R**3)
          GradW(2) = -y/(4*PI*R**3)
          GradW(3) = -z/(4*PI*R**3)
        end select
        dWdN = sum(GradW*Normal)
        do i = 1, n
          q = acc_potperm(Element % nodeindexes(i))
          if (PotentialKnown(q)) then
            if (p /= q) then
              Force(p) = Force(p) - S*Load(q)*Basis(i)*dWdN
            end if
            StiffMatrix(p, q) = StiffMatrix(p, q) - S*Basis(i)*W
          else
            Force(p) = Force(p) + S*Load(q)*Basis(i)*W
            if (p /= q) then
              StiffMatrix(p, q) = StiffMatrix(p, q) + S*Basis(i)*dWdN
            end if
          end if
          if (p /= q) then
            acc_diag(p) = acc_diag(p) - S*Basis(i)*dWdN
          end if
        end do
      end do
    end do
  end subroutine pbk_matrix

  ! -- internal procedures (orig 531-619) ----------------------------------
  subroutine pbk_potential(Potential, Pot, Flx, Element, n, Nodes)
    real(kind=dp) :: Pot(:), Flx(:), Potential(:)
    integer :: n
    type(Nodes_t) :: Nodes
    type(Element_t), target :: Element
    real(kind=dp) :: Basis(n), dBasisdx(n, 3)
    real(kind=dp) :: R, LX, LY, LZ, x, y, z, dWdN, LPOT, LFLX
    logical :: Stat, CheckNormals
    real(kind=dp) :: SqrtElementMetric, U, V, W, S, gradW(3), Normal(3)
    integer :: i, j, k, q, t, dim
    type(GaussIntegrationPoints_t) :: IntegStuff
    dim = CoordinateSystemDimension()
    select case (Element % type % elementcode/100)
    case (2)
      IntegStuff = GaussPoints(Element, 4)
    case (3)
      IntegStuff = GaussPoints(Element, 6)
    case (4)
      IntegStuff = GaussPoints(Element, 16)
    end select
    CheckNormals = associated(Element % boundaryinfo)
    if (CheckNormals) then
      CheckNormals = associated(Element % boundaryinfo % left) .or. &
        associated(Element % boundaryinfo % right)
    end if
    do t = 1, IntegStuff % n
      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)
      S = IntegStuff % s(t)
      Stat = ElementInfo(Element, Nodes, U, V, W, SqrtElementMetric, &
        Basis, dBasisdx)
      S = S*SqrtElementMetric
      Normal = NormalVector(Element, Nodes, U, V, CheckNormals)
      LX = sum(Nodes % x(1:n)*Basis)
      LY = sum(Nodes % y(1:n)*Basis)
      LZ = sum(Nodes % z(1:n)*Basis)
      do i = 1, acc_solver % mesh % numberofnodes
        k = acc_solver % variable % perm(i)
        if (k <= 0) cycle
        if (acc_potperm(i) > 0) cycle
        x = LX - acc_solver % mesh % nodes % x(i)
        y = LY - acc_solver % mesh % nodes % y(i)
        z = LZ - acc_solver % mesh % nodes % z(i)
        R = sqrt(x**2 + y**2 + z**2)
        select case (dim)
        case (2)
          W = -log(R)/(2*PI)
          GradW(1) = -x/(2*PI*R**2)
          GradW(2) = -y/(2*PI*R**2)
          GradW(3) = -z/(2*PI*R**2)
        case (3)
          W = (1/R)/(4*PI)
          GradW(1) = -x/(4*PI*R**3)
          GradW(2) = -y/(4*PI*R**3)
          GradW(3) = -z/(4*PI*R**3)
        end select
        dWdN = sum(GradW*Normal)
        do j = 1, n
          q = acc_potperm(Element % nodeindexes(j))
          Potential(k) = Potential(k) - S*Basis(j)*(Pot(q)*dWdN - Flx(q)*W)
        end do
      end do
    end do
  end subroutine pbk_potential

  ! -- SolveFull (orig 624-652) ---------------------------------------------
  subroutine pbk_solvefull()
    ! Local explicit interface for the external LAPACK driver, as in the
    ! original's CONTAINS interface, but anonymous (the Elmer
    ! LinearAlgebra shape): the rank-2 module `Matrix` crosses to the
    ! rank-1 `a(n*n)` dummy by sequence association, which generic
    ! resolution rejects but an explicit specific interface accepts —
    ! the same base address the original's rank-1 dummy saw.
    interface
      subroutine SolveLapack(N, A, x)
        integer N
        double precision a(n*n), x(n)
      end subroutine SolveLapack
    end interface
    select case (ListGetString(acc_solver % values, 'Linear System Solver'))
    case ('direct')
      call SolveLapack(acc_bn, Matrix, acc_force)
      acc_pot(1:acc_bn) = acc_force(1:acc_bn)
    case ('iterative')
      call pbk_fulliter()
    case default
      call Fatal('SolveFull', 'Unknown solver type.')
    end select
  end subroutine pbk_solvefull

  ! -- FullIterSolver (orig 657-727) ------------------------------------------
  subroutine pbk_fulliter()
    integer :: n
    real(kind=dp) :: dpar(50)
    integer :: ipar(50), wsize
    real(kind=dp), allocatable :: work(:, :)
    logical :: AbortNotConverged, GotIt
    integer(kind=AddrInt) :: AddrFunc
    external :: AddrFunc
    integer(kind=addrint) :: iterProc, mvProc, pcondProc, dProc
    n = acc_bn
    ipar = 0
    dpar = 0
    dProc = 0
    HUTI_WRKDIM = HUTI_CGS_WORKSIZE
    wsize = HUTI_WRKDIM
    HUTI_NDIM = n
    allocate (work(wsize, n))
    if (all(acc_pot(1:n) == 0.0_dp)) then
      HUTI_INITIALX = HUTI_RANDOMX
    else
      HUTI_INITIALX = HUTI_USERSUPPLIEDX
    end if
    HUTI_TOLERANCE = ListGetConstReal(acc_solver % values, &
      'Linear System Convergence Tolerance')
    HUTI_MAXTOLERANCE = ListGetConstReal(acc_solver % values, &
      'Linear System Divergence Limit', GotIt)
    if (.not. GotIt) HUTI_MAXTOLERANCE = 1.0d20
    HUTI_MAXIT = ListGetInteger(acc_solver % values, &
      'Linear System Max Iterations')
    HUTI_DBUGLVL = ListGetInteger(acc_solver % values, &
      'Linear System Residual Output', GotIt)
    if (.not. GotIt) HUTI_DBUGLVL = 1
    AbortNotConverged = ListGetLogical(acc_solver % values, &
      'Linear System Abort Not Converged', GotIt)
    if (.not. GotIt) AbortNotConverged = .true.
    iterProc = AddrFunc(HUTI_D_CGS)
    mvProc = AddrFunc(Matvec)
    pcondProc = AddrFunc(Precond)
    call IterCall(iterProc, acc_pot, acc_force, ipar, dpar, work, mvProc, pcondProc, &
      dProc, dProc, dProc, dProc)
    deallocate (work)
    if (HUTI_INFO /= HUTI_CONVERGENCE) then
      if (AbortNotConverged) then
        call Fatal('IterSolve', 'Failed convergence tolerances.')
      else
        call Error('IterSolve', 'Failed convergence tolerances.')
      end if
    end if
  end subroutine pbk_fulliter

end module poissonbem_acc
