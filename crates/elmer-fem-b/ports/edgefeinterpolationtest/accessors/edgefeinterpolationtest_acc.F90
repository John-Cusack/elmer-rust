! T022-A EdgeFEInterpolationTest port: statement-level BIND(C) wrappers for
! the Elmer calls in fem/src/modules/EdgeFEInterpolationTest.F90, compiled
! against the built .mod files of the target configure (A2 shipped
! mechanism, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) /
! INTEGER(C_INT) arrays). All Fortran-object traffic — TYPE(...)
! handles, OPTIONAL/POINTER dummies, hidden lengths, assumed-shape
! descriptors, list-directed timing reports — stays inside Fortran, so the
! Rust side is pure C. The wrapper module mirrors
! EdgeFEInterpolationTest's own SAVE state (entry-level Piola/style flags,
! LocalMatrix SAVE Nodes + IP local, ComputeApproximationError SAVE Nodes
! + IP local + EK/SolNorm accumulators) with identical lifetimes.
! Order-sensitive reductions (SUM), the MATMUL-free Galerkin assembly, the
! formatted `WRITE` of the error norm, and whole-array updates stay in
! Fortran (same compiler), so bit-identity holds by construction wherever
! the call sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module edgefeinterpolationtest_acc
  use iso_c_binding
  use defutils
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(matrix_t), pointer, save :: acc_a => null()
  ! Entry-level SAVE (mirrors the solver body's locals).
  logical, save :: acc_piola = .false.
  logical, save :: acc_second = .false.
  integer(c_int), save :: acc_order = 1
  integer(c_int), save :: acc_testmode = 2
  logical, save :: acc_errest = .true.
  logical, save :: acc_curl = .false.
  integer(c_int), save :: acc_dim = 0
  ! LocalMatrix context (mirrors its SAVE Nodes and IP local).
  type(nodes_t), target, save :: acc_nodes_lm
  type(gaussintegrationpoints_t), save :: acc_ip_lm
  ! ComputeApproximationError context (mirrors its SAVE Nodes and the
  ! parent-scope EK/SolNorm accumulators).
  type(nodes_t), target, save :: acc_nodes_ce
  type(gaussintegrationpoints_t), save :: acc_ip_ce
  real(c_double), save :: acc_ek = 0.0_c_double
  real(c_double), save :: acc_solnorm = 0.0_c_double

contains

  ! -- entry state ----------------------------------------------------
  subroutine efi_begin(h_solver) bind(c, name="efi_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine efi_begin

  ! -- header: EdgeElementStyle + ValueList reads (mirrors the entry head)
  subroutine efi_header() bind(c, name="efi_header")
    type(valuelist_t), pointer :: params
    logical :: pv, sf, found, lval
    integer :: k, ival
    params => GetSolverParams()
    call EdgeElementStyle(params, pv, SecondFamily=sf, BasisDegree=k)
    acc_piola = pv
    acc_second = sf
    acc_order = k
    lval = GetLogical(params, 'Simplicial Mesh', found)
    lval = GetLogical(params, 'Error Computation', found)
    if (.not. found) lval = .true.
    acc_errest = lval
    acc_curl = GetLogical(params, 'Compute Curl Error', found)
    ival = ListGetInteger(params, 'Test Mode', found)
    if (.not. found) ival = 2
    acc_testmode = ival
    acc_dim = CoordinateSystemDimension()
  end subroutine efi_header

  function efi_error_on() result(r) bind(c, name="efi_error_on")
    integer(c_int) :: r
    r = merge(1, 0, acc_errest)
  end function efi_error_on

  subroutine efi_touch_mesh() bind(c, name="efi_touch_mesh")
    type(mesh_t), pointer :: m
    m => GetMesh()
  end subroutine efi_touch_mesh

  function efi_max_dofs() result(n) bind(c, name="efi_max_dofs")
    integer(c_int) :: n
    n = acc_solver % mesh % maxelementdofs
  end function efi_max_dofs

  subroutine efi_matrix_setup() bind(c, name="efi_matrix_setup")
    acc_solver % matrix % complex = .false.
    acc_a => GetMatrix()
  end subroutine efi_matrix_setup

  ! -- assembly driver -------------------------------------------------
  function efi_nof_active() result(n) bind(c, name="efi_nof_active")
    integer(c_int) :: n
    n = GetNOFActive()
  end function efi_nof_active

  subroutine efi_default_initialize() bind(c, name="efi_default_initialize")
    call DefaultInitialize()
  end subroutine efi_default_initialize

  subroutine efi_begin_active(t, n, nd, nb) bind(c, name="efi_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, nb
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
    nb = GetElementNOFBDOFs()
  end subroutine efi_begin_active

  ! Acoef(1:n) = 0, then Material => GetMaterial; when associated,
  ! Acoef(1:n) = GetReal(Material, 'Material Param', Found) (Found is
  ! ignored by the original too).
  subroutine efi_acoef(n, acoef) bind(c, name="efi_acoef")
    integer(c_int), value :: n
    real(c_double) :: acoef(n)
    type(valuelist_t), pointer :: mat
    logical :: found
    acoef(1:n) = 0.0_c_double
    mat => GetMaterial(acc_elem)
    if (associated(mat)) then
      acoef(1:n) = GetReal(mat, 'Material Param', found)
    end if
  end subroutine efi_acoef

  ! -- LocalMatrix internal procedure (exact transcription) --------------
  subroutine efi_local_matrix(n, nd, ld, stiff, force, acoef) &
      bind(c, name="efi_local_matrix")
    integer(c_int), value :: n, nd, ld
    real(c_double) :: stiff(ld, ld), force(ld), acoef(ld)
    real(c_double) :: matpar
    real(c_double) :: ebasis(nd, 3), curlebasis(nd, 3)
    real(c_double) :: basis(n), detj, xq, yq, zq
    logical :: stat
    integer :: t, i, j, p, q, np
    acc_ip_lm = GaussPoints(acc_elem, PReferenceElement=acc_piola, &
        EdgeBasisDegree=acc_order)
    call GetElementNodes(acc_nodes_lm)
    stiff = 0.0_c_double
    force = 0.0_c_double
    np = 0
    do t = 1, acc_ip_lm % n
      stat = ElementInfo(acc_elem, acc_nodes_lm, acc_ip_lm % u(t), &
          acc_ip_lm % v(t), acc_ip_lm % w(t), detj, basis, &
          EdgeBasis=ebasis, RotBasis=curlebasis)
      xq = sum(acc_nodes_lm % x(1:n) * basis(1:n))
      yq = sum(acc_nodes_lm % y(1:n) * basis(1:n))
      zq = sum(acc_nodes_lm % z(1:n) * basis(1:n))
      matpar = sum(acoef(1:n) * basis(1:n))
      if (np > 0) then
        do p = 1, n
          do q = 1, n
            stiff(p, q) = stiff(p, q) + basis(p) * basis(q) * detj * acc_ip_lm % s(t)
          end do
          do q = 1, nd - np
            j = np + q
            stiff(p, j) = stiff(p, j) - ebasis(q, 3) * basis(p) * detj * acc_ip_lm % s(t)
          end do
        end do
      end if
      do p = 1, nd - np
        i = np + p
        do q = 1, nd - np
          j = np + q
          stiff(i, j) = stiff(i, j) + 1.0_c_double * &
              sum(ebasis(q, 1:acc_dim) * ebasis(p, 1:acc_dim)) * detj * acc_ip_lm % s(t) + &
              matpar * sum(curlebasis(q, 1:3) * curlebasis(p, 1:3)) * detj * acc_ip_lm % s(t)
        end do
        select case (acc_testmode)
        case (1)
          force(i) = force(i) + (1.0_c_double) * ebasis(p, 1) * detj * acc_ip_lm % s(t) + &
              (1.0_c_double) * ebasis(p, 2) * detj * acc_ip_lm % s(t) + &
              (1.0_c_double) * ebasis(p, 3) * detj * acc_ip_lm % s(t)
        case (2)
          force(i) = force(i) + (1.0_c_double + zq - yq) * ebasis(p, 1) * detj * acc_ip_lm % s(t) + &
              (1.0_c_double + xq - zq) * ebasis(p, 2) * detj * acc_ip_lm % s(t) + &
              (1.0_c_double - xq + yq) * ebasis(p, 3) * detj * acc_ip_lm % s(t) + &
              matpar * (2.0_c_double) * curlebasis(p, 1) * detj * acc_ip_lm % s(t) + &
              matpar * (2.0_c_double) * curlebasis(p, 2) * detj * acc_ip_lm % s(t) + &
              matpar * (2.0_c_double) * curlebasis(p, 3) * detj * acc_ip_lm % s(t)
        case (3)
          force(i) = force(i) + (0.0_c_double) * ebasis(p, 1) * detj * acc_ip_lm % s(t) + &
              (0.0_c_double) * ebasis(p, 2) * detj * acc_ip_lm % s(t) - &
              0.5_c_double * (yq * xq**2 + xq * yq**2) * ebasis(p, 3) * detj * acc_ip_lm % s(t) + &
              matpar * (-0.5_c_double * xq**2 - yq * xq) * curlebasis(p, 1) * detj * acc_ip_lm % s(t) + &
              matpar * (0.5_c_double * yq**2 + yq * xq) * curlebasis(p, 2) * detj * acc_ip_lm % s(t) + &
              matpar * (0.0_c_double) * curlebasis(p, 3) * detj * acc_ip_lm % s(t)
        case (4)
          force(i) = force(i) + (xq * yq**2) * ebasis(p, 1) * detj * acc_ip_lm % s(t) + &
              (yq * xq**2) * ebasis(p, 2) * detj * acc_ip_lm % s(t)
        case (5)
          force(i) = force(i) + (-xq * yq**2) * ebasis(p, 1) * detj * acc_ip_lm % s(t) + &
              (yq * xq**2) * ebasis(p, 2) * detj * acc_ip_lm % s(t) + &
              matpar * (0.0_c_double) * curlebasis(p, 1) * detj * acc_ip_lm % s(t) + &
              matpar * (0.0_c_double) * curlebasis(p, 2) * detj * acc_ip_lm % s(t) + &
              matpar * (4.0_c_double * xq * yq) * curlebasis(p, 3) * detj * acc_ip_lm % s(t)
        case (6)
          force(i) = force(i) + (1.5_c_double * yq**2) * ebasis(p, 1) * detj * acc_ip_lm % s(t) + &
              (0.5_c_double * xq**2) * ebasis(p, 2) * detj * acc_ip_lm % s(t) + &
              matpar * (0.0_c_double) * curlebasis(p, 1) * detj * acc_ip_lm % s(t) + &
              matpar * (0.0_c_double) * curlebasis(p, 2) * detj * acc_ip_lm % s(t) + &
              matpar * (xq - 3.0_c_double * yq) * curlebasis(p, 3) * detj * acc_ip_lm % s(t)
        end select
      end do
    end do
  end subroutine efi_local_matrix

  ! -- shared assembly tail ---------------------------------------------
  subroutine efi_default_update(ld, g, f) bind(c, name="efi_default_update")
    integer(c_int), value :: ld
    real(c_double) :: g(ld, ld), f(ld)
    call DefaultUpdateEquations(g(1:ld, 1:ld), f(1:ld))
  end subroutine efi_default_update

  subroutine efi_default_finish_bulk() bind(c, name="efi_default_finish_bulk")
    call DefaultFinishBulkAssembly()
  end subroutine efi_default_finish_bulk

  subroutine efi_default_finish() bind(c, name="efi_default_finish")
    call DefaultFinishAssembly()
  end subroutine efi_default_finish

  subroutine efi_default_dirichlet() bind(c, name="efi_default_dirichlet")
    call DefaultDirichletBCs()
  end subroutine efi_default_dirichlet

  function efi_default_solve() result(norm) bind(c, name="efi_default_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function efi_default_solve

  subroutine efi_fatal_alloc() bind(c, name="efi_fatal_alloc")
    call Fatal('BestApproximationSolver', 'Memory allocation error.')
  end subroutine efi_fatal_alloc

  ! -- error estimation ---------------------------------------------------
  subroutine efi_error_begin() bind(c, name="efi_error_begin")
    acc_ek = 0.0_c_double
    acc_solnorm = 0.0_c_double
  end subroutine efi_error_begin

  function efi_num_active_elems() result(n) bind(c, name="efi_num_active_elems")
    integer(c_int) :: n
    n = acc_solver % numberofactiveelements
  end function efi_num_active_elems

  subroutine efi_error_elem(t, nmax, load, indices) bind(c, name="efi_error_elem")
    integer(c_int), value :: t, nmax
    real(c_double) :: load(6, nmax)
    integer(c_int) :: indices(nmax)
    integer :: n, nd
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes()
    nd = GetElementDOFs(indices)
    load(1, 1:nd) = acc_solver % variable % values( &
        acc_solver % variable % perm(indices(1:nd)))
    call efi_cerr(load, n, nd, nmax)
  end subroutine efi_error_elem

  ! ComputeApproximationError internal procedure (exact transcription).
  subroutine efi_cerr(load, n, nd, nmax)
    real(c_double) :: load(6, nmax)
    integer :: n, nd, nmax
    real(c_double) :: ebasis(nd, 3), curlebasis(nd, 3)
    real(c_double) :: basis(n), detj, xq, yq, zq
    real(c_double) :: u(3), rotu(3), sol(3), rotsol(3), e(3), rote(3)
    logical :: stat
    integer :: t, i, np
    acc_ip_ce = GaussPoints(acc_elem, PReferenceElement=acc_piola, &
        EdgeBasisDegree=acc_order)
    call GetElementNodes(acc_nodes_ce)
    np = 0
    do t = 1, acc_ip_ce % n
      stat = ElementInfo(acc_elem, acc_nodes_ce, acc_ip_ce % u(t), &
          acc_ip_ce % v(t), acc_ip_ce % w(t), detj, basis, &
          EdgeBasis=ebasis, RotBasis=curlebasis)
      xq = sum(acc_nodes_ce % x(1:n) * basis(1:n))
      yq = sum(acc_nodes_ce % y(1:n) * basis(1:n))
      zq = sum(acc_nodes_ce % z(1:n) * basis(1:n))
      u = 0.0_c_double
      do i = 1, acc_dim
        u(i) = sum(load(1, np + 1:nd) * ebasis(1:nd - np, i))
      end do
      rotu(1) = sum(load(1, np + 1:nd) * curlebasis(1:nd - np, 1))
      rotu(2) = sum(load(1, np + 1:nd) * curlebasis(1:nd - np, 2))
      rotu(3) = sum(load(1, np + 1:nd) * curlebasis(1:nd - np, 3))
      select case (acc_testmode)
      case (1)
        sol(1) = 1.0_c_double
        sol(2) = 1.0_c_double
        sol(3) = 1.0_c_double
        rotsol(1:3) = 0.0_c_double
      case (2)
        sol(1) = 1.0_c_double + zq - yq
        sol(2) = 1.0_c_double - zq + xq
        sol(3) = 1.0_c_double - xq + yq
        rotsol(1:3) = 2.0_c_double
      case (3)
        sol(1:2) = 0.0_c_double
        sol(3) = -0.5_c_double * (yq * xq**2 + xq * yq**2)
        rotsol(1) = -0.5_c_double * xq**2 - yq * xq
        rotsol(2) = 0.5_c_double * yq**2 + yq * xq
        rotsol(3) = 0.0_c_double
      case (4)
        sol(1) = xq * yq**2
        sol(2) = yq * xq**2
        sol(3) = 0.0_c_double
        rotsol = 0.0_c_double
      case (5)
        sol(1) = -xq * yq**2
        sol(2) = yq * xq**2
        sol(3) = 0.0_c_double
        rotsol(1:3) = 0.0_c_double
        rotsol(3) = 4.0_c_double * xq * yq
      case (6)
        sol(1) = 1.5_c_double * yq**2
        sol(2) = 0.5_c_double * xq**2
        sol(3) = 0.0_c_double
        rotsol(1:3) = 0.0_c_double
        rotsol(3) = xq - 3.0_c_double * yq
      end select
      e(:) = sol(:) - u(:)
      rote(:) = rotsol(:) - rotu(:)
      if (acc_curl) then
        acc_solnorm = acc_solnorm + sum(rotsol(1:3) * rotsol(1:3)) * detj * acc_ip_ce % s(t)
        acc_ek = acc_ek + sum(rote(1:3) * rote(1:3)) * detj * acc_ip_ce % s(t)
      else
        acc_solnorm = acc_solnorm + sum(sol(1:3) * sol(1:3)) * detj * acc_ip_ce % s(t)
        acc_ek = acc_ek + sum(e(1:3) * e(1:3)) * detj * acc_ip_ce % s(t)
      end if
    end do
  end subroutine efi_cerr

  subroutine efi_error_end() bind(c, name="efi_error_end")
    write (*, '(A,E16.8)') 'Error Norm = ', &
        sqrt(ParallelReduction(acc_ek)) / sqrt(ParallelReduction(acc_solnorm))
  end subroutine efi_error_end

  ! -- BestApproximationSolver_Init0 body (exact transcription) ------------
  subroutine efi_init_begin(h_solver) bind(c, name="efi_init_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine efi_init_begin

  subroutine efi_init_body() bind(c, name="efi_init_body")
    type(valuelist_t), pointer :: solverparams
    logical :: found, piolaversion, secondfamily, check
    integer :: k
    solverparams => GetSolverParams()
    if (.not. ListCheckPresent(solverparams, 'Element')) then
      call EdgeElementStyle(solverparams, piolaversion, secondfamily, &
          BasisDegree=k, Check=.true.)
      select case (k)
      case (3)
        if (secondfamily) then
          call Fatal('BestApproximationSolver', &
              'No ready support for the cubic element of the second kind')
        else
          call ListAddString(solverparams, 'Element', &
              'n:0 e:3 -tri b:6 -tetra b:6 -tri_face b:6')
        end if
      case (2)
        if (secondfamily) then
          call ListAddString(solverparams, 'Element', &
              'n:0 e:3 -tri b:3 -tri_face b:3')
        else
          call ListAddString(solverparams, 'Element', &
              'n:0 e:2 -tri b:2 -quad b:4 -brick b:6 -pyramid b:3 -prism b:2 -quad_face b:4 -tri_face b:2')
        end if
      case default
        if (secondfamily) then
          call ListAddString(solverparams, 'Element', 'n:0 e:2')
        else if (piolaversion) then
          call ListAddString(solverparams, 'Element', 'n:0 e:1 -quad b:2 -brick b:3 -quad_face b:2')
        else
          call ListAddString(solverparams, 'Element', 'n:0 e:1')
        end if
      end select
    end if
  end subroutine efi_init_body

end module edgefeinterpolationtest_acc
