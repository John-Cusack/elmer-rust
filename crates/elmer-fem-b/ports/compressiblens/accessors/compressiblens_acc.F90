! T022-A CompressibleNS port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/CompressibleNS.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 section 3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, list-directed timing
! reports — stays inside Fortran, so the Rust side is pure C calls. The
! wrapper module mirrors CompressibleNS's own locals (current element,
! LocalMatrix SAVE Nodes, entry-level flags and timings) with identical
! lifetimes. All floating-point work (the LocalMatrix Galerkin assembly,
! MATMUL/SUM reductions, the CPUTime-based reports) stays in Fortran
! (same compiler), so bit-identity holds by construction wherever the
! call sequence matches.
!
! Statement-level transcription, no upstream text.
module compressiblens_acc
  use iso_c_binding
  use coordinatesystems
  use defutils
  use solverbasics
  use elementutils
  use messages
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(valuelist_t), pointer, save :: acc_bf => null()
  type(valuelist_t), pointer, save :: acc_mat => null()
  ! LocalMatrix context (mirrors its SAVE Nodes with identical lifetime:
  ! filled by GetElementNodes, consumed before the next fill).
  type(nodes_t), target, save :: acc_nodes
  ! Entry-level state (mirrors the solver body's locals).
  integer(c_int), save :: acc_dim = 0
  integer(c_int), save :: acc_coordsys = 0
  real(c_double), save :: acc_atime = 0.0_c_double
  real(c_double), save :: acc_at0 = 0.0_c_double

contains

  ! -- entry state ----------------------------------------------------
  subroutine cns_begin(h_solver) bind(c, name="cns_begin")
    type(c_ptr), value :: h_solver
    character(len=MAX_NAME_LEN) :: itmeth
    logical :: gotit, blk
    call c_f_pointer(h_solver, acc_solver)
    acc_dim = CoordinateSystemDimension()
    acc_coordsys = CurrentCoordinateSystem()
    ! Pure keyword reads, exactly as in the original. The iteration-method
    ! choice is forced to Picard and the block preconditioning is forced
    ! off inside the assembly, as in the original.
    itmeth = ListGetString(acc_solver % values, &
        'Nonlinear Iteration Method', gotit)
    blk = ListGetLogical(acc_solver % values, 'Block Preconditioning', gotit)
    ! Both reads are pure; the forced values below match the original.
    blk = .false.
    acc_atime = CPUTime()
    acc_at0 = RealTime()
  end subroutine cns_begin

  function cns_max_dofs() result(n) bind(c, name="cns_max_dofs")
    integer(c_int) :: n
    n = acc_solver % mesh % maxelementdofs
  end function cns_max_dofs

  function cns_coord_dim() result(d) bind(c, name="cns_coord_dim")
    integer(c_int) :: d
    d = acc_dim
  end function cns_coord_dim

  function cns_done_time() result(d) bind(c, name="cns_done_time")
    integer(c_int) :: d
    d = acc_solver % donetime
  end function cns_done_time

  subroutine cns_nl_params(nliter, nltol) bind(c, name="cns_nl_params")
    integer(c_int) :: nliter
    real(c_double) :: nltol
    nliter = ListGetInteger(acc_solver % values, &
        'Nonlinear System Max Iterations', minv=0)
    nltol = ListGetConstReal(acc_solver % values, &
        'Nonlinear System Convergence Tolerance', minv=0.0_dp)
  end subroutine cns_nl_params

  subroutine cns_default_start() bind(c, name="cns_default_start")
    call DefaultStart()
  end subroutine cns_default_start

  subroutine cns_default_initialize() bind(c, name="cns_default_initialize")
    call DefaultInitialize()
  end subroutine cns_default_initialize

  function cns_nof_active() result(n) bind(c, name="cns_nof_active")
    integer(c_int) :: n
    n = acc_solver % numberofactiveelements
  end function cns_nof_active

  ! -- active-element context ------------------------------------------
  subroutine cns_begin_elem(t, n, nd, nb, idx, psz) &
      bind(c, name="cns_begin_elem")
    integer(c_int), value :: t, psz
    integer(c_int) :: n, nd, nb
    integer(c_int) :: idx(psz)
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes()
    nd = GetElementDOFs(idx)
    nb = GetElementNOFBDOFs()
  end subroutine cns_begin_elem

  function cns_bodyforce() result(r) bind(c, name="cns_bodyforce")
    integer(c_int) :: r
    acc_bf => GetBodyForce()
    r = merge(1, 0, associated(acc_bf))
  end function cns_bodyforce

  subroutine cns_load_cols(n, pld, load) bind(c, name="cns_load_cols")
    integer(c_int), value :: n, pld
    real(c_double) :: load(pld, 5)
    logical :: found
    load(1:n,1) = GetReal(acc_bf, 'Body Force 1', found)
    load(1:n,2) = GetReal(acc_bf, 'Body Force 2', found)
    load(1:n,3) = GetReal(acc_bf, 'Body Force 3', found)
  end subroutine cns_load_cols

  ! -- material parameters (mirrors the per-element material block) -----
  subroutine cns_material(n, pld, mu, cv, k, gam, lam, t0, rho0) &
      bind(c, name="cns_material")
    integer(c_int), value :: n, pld
    real(c_double) :: mu(pld)
    real(c_double) :: cv, k, gam, lam, t0, rho0
    real(c_double) :: bulk
    logical :: found
    acc_mat => GetMaterial()
    mu(1:n) = GetReal(acc_mat, 'Viscosity')
    cv = ListGetConstReal(acc_mat, 'Specific Heat')
    k = ListGetConstReal(acc_mat, 'Heat Conductivity')
    gam = ListGetConstReal(acc_mat, 'Specific Heat Ratio')
    bulk = ListGetConstReal(acc_mat, 'Bulk Viscosity', found)
    if (.not. found) bulk = 0.0_dp
    lam = bulk - 2.0_dp/3.0_dp * mu(1)
    t0 = ListGetConstReal(acc_mat, 'Equilibrium Temperature')
    rho0 = ListGetConstReal(acc_mat, 'Equilibrium Density')
  end subroutine cns_material

  ! -- previous/current iterate read (mirrors the FirstSolve branch) ---
  subroutine cns_elem_values(first, newts, iter, d1, pld, nd, velo, temp, &
      rho, idx) bind(c, name="cns_elem_values")
    integer(c_int), value :: first, newts, iter, d1, pld, nd
    real(c_double) :: velo(d1, pld), temp(pld), rho(pld)
    integer(c_int) :: idx(pld)
    integer :: i
    if (first /= 0 .and. iter == 1) then
      temp(1:nd) = 0.0_dp
      rho(1:nd) = 0.0_dp
      do i = 1, acc_dim
        velo(i,1:nd) = 0.0_dp
      end do
    else
      if (newts /= 0 .and. iter == 1) then
        do i = 1, acc_dim
          velo(i,1:nd) = acc_solver % variable % prevvalues( &
              acc_solver % variable % dofs*(acc_solver % variable % &
              perm(idx(1:nd))-1)+i,1)
        end do
        temp(1:nd) = acc_solver % variable % prevvalues( &
            acc_solver % variable % dofs*(acc_solver % variable % &
            perm(idx(1:nd))-1)+acc_dim+1,1)
        rho(1:nd) = acc_solver % variable % prevvalues( &
            acc_solver % variable % dofs*(acc_solver % variable % &
            perm(idx(1:nd))-1)+acc_dim+2,1)
      else
        do i = 1, acc_dim
          velo(i,1:nd) = acc_solver % variable % values( &
              acc_solver % variable % dofs*(acc_solver % variable % &
              perm(idx(1:nd))-1)+i)
        end do
        temp(1:nd) = acc_solver % variable % values( &
            acc_solver % variable % dofs*(acc_solver % variable % &
            perm(idx(1:nd))-1)+acc_dim+1)
        rho(1:nd) = acc_solver % variable % values( &
            acc_solver % variable % dofs*(acc_solver % variable % &
            perm(idx(1:nd))-1)+acc_dim+2)
      end if
    end if
  end subroutine cns_elem_values

  ! -- LocalMatrix internal procedure (exact transcription) --------------
  ! Dummy mapping: STIFF/Mass/FORCE/LOAD/NodalT/Nodalrho/Nodalmu/NodalVelo
  ! keep their roles; PicardIteration/Stabilization are folded to the
  ! forced values the procedure assigns on entry (.TRUE./.FALSE.).
  subroutine cns_local_matrix(nloc, pld, d1, stiff, mass, force, load, &
      tnod, rnod, mnod, vnod, n, nd, dd, gam, cv, k, la, t0, rho0) &
      bind(c, name="cns_local_matrix")
    integer(c_int), value :: nloc, pld, d1, n, nd, dd
    real(c_double), target :: stiff(nloc,nloc), mass(nloc,nloc)
    real(c_double), target :: force(nloc), load(pld,5)
    real(c_double) :: tnod(pld), rnod(pld), mnod(pld), vnod(d1,pld)
    real(c_double), value :: gam, cv, k, la, t0, rho0
    real(c_double) :: basis(nd), dbasis(nd,3), &
        detj, loadatip(dd+2), velo(dd), grad(dd,dd), ak, w3
    real(c_double), pointer :: a(:,:), f(:), m(:,:)
    logical :: stat
    integer :: t, i, j, l, p, q
    type(gaussintegrationpoints_t) :: ip
    real(c_double) :: mu, rho, s, c, c2, ch, rotterm, prevt, r, &
        gradv1(2), gradv2(2), rpos, testterm
    logical :: stab
    ! LocalMatrix forces these on entry, overriding the caller's values.
    stab = .false.
    call GetElementNodes(acc_nodes)
    stiff = 0.0_dp
    force = 0.0_dp
    mass = 0.0_dp
    ip = GaussPoints(acc_elem)
    r = (gam - 1.0_dp) * cv
    do t = 1, ip % n
      stat = ElementInfo(acc_elem, acc_nodes, ip % u(t), ip % v(t), &
          ip % w(t), detj, basis, dbasis)
      s = ip % s(t) * detj
      if (acc_coordsys == AxisSymmetric) then
        rpos = sum(basis(1:n) * acc_nodes % x(1:n))
        s = s * rpos
      end if
      velo = matmul(vnod(1:dd,1:nd), basis(1:nd))
      gradv1(1) = sum(vnod(1,1:nd) * dbasis(1:nd,1))
      gradv1(2) = sum(vnod(1,1:nd) * dbasis(1:nd,2))
      gradv2(1) = sum(vnod(2,1:nd) * dbasis(1:nd,1))
      gradv2(2) = sum(vnod(2,1:nd) * dbasis(1:nd,2))
      mu = sum(basis(1:n) * mnod(1:n))
      rho = sum(basis(1:nd) * rnod(1:nd)) + rho0
      prevt = sum(basis(1:nd) * tnod(1:nd)) + t0
      loadatip(1:dd+2) = matmul(basis(1:n), load(1:n,1:dd+2))
      do p = 1, nd
        do q = 1, nd
          i = (dd+2) * (p-1) + 1
          j = (dd+2) * (q-1) + 1
          a => stiff(i:i+dd+1,j:j+dd+1)
          m => mass(i:i+dd+1,j:j+dd+1)
          do i = 1, dd
            if ((acc_coordsys == AxisSymmetric) .and. i == 1) then
              testterm = dbasis(p,i) + basis(p)/rpos
            else
              testterm = dbasis(p,i)
            end if
            do j = 1, dd
              a(i,i) = a(i,i) + s * mu * dbasis(q,j) * dbasis(p,j)
              a(i,j) = a(i,j) + s * mu * dbasis(q,i) * dbasis(p,j)
              a(i,j) = a(i,j) + s * la * dbasis(q,j) * testterm
              if ((acc_coordsys == AxisSymmetric) .and. j == 1) then
                a(i,j) = a(i,j) + s * la * 1/rpos * basis(q) * testterm
              end if
            end do
            if ((acc_coordsys == AxisSymmetric) .and. i == 1) then
              a(i,i) = a(i,i) + s * 2.0_dp * mu * 1/rpos**2 * &
                  basis(q) * basis(p)
            end if
            m(i,i) = m(i,i) + s * rho * basis(p) * basis(q)
            if (q <= n) then
              a(i,dd+1) = a(i,dd+1) - s * r * rho * basis(q) * dbasis(p,i)
              a(i,dd+2) = a(i,dd+2) - s * r * prevt * basis(q) * dbasis(p,i)
              if ((acc_coordsys == AxisSymmetric) .and. i == 1) then
                a(i,dd+1) = a(i,dd+1) - s * r * rho * basis(q) * &
                    basis(p)/rpos
                a(i,dd+2) = a(i,dd+2) - s * r * prevt * basis(q) * &
                    basis(p)/rpos
              end if
            end if
            if (p <= n) then
              a(dd+1,i) = a(dd+1,i) + s * r * rho * dbasis(q,i) * basis(p)
              a(dd+2,i) = a(dd+2,i) + s * rho * dbasis(q,i) * basis(p)
              if ((acc_coordsys == AxisSymmetric) .and. i == 1) then
                a(dd+1,i) = a(dd+1,i) + s * r * rho * basis(q)/rpos * &
                    basis(p)
                a(dd+2,i) = a(dd+2,i) + s * rho * basis(q)/rpos * basis(p)
              end if
            end if
          end do
          m(dd+1,dd+1) = m(dd+1,dd+1) + s * rho * cv / prevt * &
              basis(p) * basis(q)
          m(dd+2,dd+2) = m(dd+2,dd+2) + s * basis(p) * basis(q)
          a(dd+1,dd+1) = a(dd+1,dd+1) + s * k / prevt * &
              sum(dbasis(p,1:dd) * dbasis(q,1:dd))
          a(dd+1,dd+1) = a(dd+1,dd+1) + s * rho * cv / prevt * &
              (velo(1) * dbasis(q,1) + velo(2) * dbasis(q,2)) * basis(p)
          a(dd+2,dd+2) = a(dd+2,dd+2) + s * &
              (velo(1) * dbasis(q,1) + velo(2) * dbasis(q,2)) * basis(p)
          if (dd == 3) then
            a(dd+1,dd+1) = a(dd+1,dd+1) + s * rho * cv / prevt * &
                velo(3) * dbasis(q,3) * basis(p)
            a(dd+2,dd+2) = a(dd+2,dd+2) + s * velo(3) * dbasis(q,3) * &
                basis(p)
          end if
          if (.false.) then
            stiff((dd+2)*(p-1)+1, (dd+2)*(q-1)+1) = &
                stiff((dd+2)*(p-1)+1, (dd+2)*(q-1)+1) + &
                rho * s * velo(2) * dbasis(q,2) * basis(p)
            stiff((dd+2)*(p-1)+1, (dd+2)*(q-1)+2) = &
                stiff((dd+2)*(p-1)+1, (dd+2)*(q-1)+2) - &
                rho * s * velo(2) * dbasis(q,1) * basis(p)
            stiff((dd+2)*(p-1)+2, (dd+2)*(q-1)+1) = &
                stiff((dd+2)*(p-1)+2, (dd+2)*(q-1)+1) - &
                rho * s * velo(1) * dbasis(q,2) * basis(p)
            stiff((dd+2)*(p-1)+2, (dd+2)*(q-1)+2) = &
                stiff((dd+2)*(p-1)+2, (dd+2)*(q-1)+2) + &
                rho * s * velo(1) * dbasis(q,1) * basis(p)
          else
            stiff((dd+2)*(p-1)+1, (dd+2)*(q-1)+1) = &
                stiff((dd+2)*(p-1)+1, (dd+2)*(q-1)+1) + &
                rho * s * velo(2) * dbasis(q,2) * basis(p)
            stiff((dd+2)*(p-1)+1, (dd+2)*(q-1)+1) = &
                stiff((dd+2)*(p-1)+1, (dd+2)*(q-1)+1) + &
                rho * s * velo(1) * dbasis(q,1) * basis(p)
            stiff((dd+2)*(p-1)+2, (dd+2)*(q-1)+2) = &
                stiff((dd+2)*(p-1)+2, (dd+2)*(q-1)+2) + &
                rho * s * velo(1) * dbasis(q,1) * basis(p)
            stiff((dd+2)*(p-1)+2, (dd+2)*(q-1)+2) = &
                stiff((dd+2)*(p-1)+2, (dd+2)*(q-1)+2) + &
                rho * s * velo(2) * dbasis(q,2) * basis(p)
            if (dd > 2) then
              stiff((dd+2)*(p-1)+1, (dd+2)*(q-1)+1) = &
                  stiff((dd+2)*(p-1)+1, (dd+2)*(q-1)+1) + &
                  rho * s * velo(3) * dbasis(q,3) * basis(p)
              stiff((dd+2)*(p-1)+2, (dd+2)*(q-1)+2) = &
                  stiff((dd+2)*(p-1)+2, (dd+2)*(q-1)+2) + &
                  rho * s * velo(2) * dbasis(q,3) * basis(p)
              stiff((dd+2)*(p-1)+3, (dd+2)*(q-1)+3) = &
                  stiff((dd+2)*(p-1)+3, (dd+2)*(q-1)+3) + &
                  rho * s * velo(1) * dbasis(q,1) * basis(p)
              stiff((dd+2)*(p-1)+3, (dd+2)*(q-1)+3) = &
                  stiff((dd+2)*(p-1)+3, (dd+2)*(q-1)+3) + &
                  rho * s * velo(2) * dbasis(q,2) * basis(p)
              stiff((dd+2)*(p-1)+3, (dd+2)*(q-1)+3) = &
                  stiff((dd+2)*(p-1)+3, (dd+2)*(q-1)+3) + &
                  rho * s * velo(3) * dbasis(q,3) * basis(p)
            end if
          end if
        end do
        i = (dd+2) * (p-1) + 1
        f => force(i:i+dd+1)
        f = f + s * loadatip * basis(p)
      end do
    end do
    do p = n+1, nd
      i = (dd+2) * p
      force(i) = 0.0_dp
      mass(:,i) = 0.0_dp
      mass(i,:) = 0.0_dp
      stiff(i,:) = 0.0_dp
      stiff(:,i) = 0.0_dp
      force(i-1) = 0.0_dp
      mass(:,i-1) = 0.0_dp
      mass(i-1,:) = 0.0_dp
      stiff(i-1,:) = 0.0_dp
      stiff(:,i-1) = 0.0_dp
      if (.not. stab) then
        stiff(i,i) = 1.0_dp
        stiff(i-1,i-1) = 1.0_dp
      end if
    end do
  end subroutine cns_local_matrix

  ! -- transient mass assembly (steady runs never call this) ------------
  subroutine cns_first_order(nloc, mass, stiff, force) &
      bind(c, name="cns_first_order")
    integer(c_int), value :: nloc
    real(c_double) :: mass(nloc,nloc), stiff(nloc,nloc), force(nloc)
    call Default1stOrderTime(mass, stiff, force)
  end subroutine cns_first_order

  ! -- shared assembly tail ---------------------------------------------
  subroutine cns_default_update(nloc, stiff, force) &
      bind(c, name="cns_default_update")
    integer(c_int), value :: nloc
    real(c_double) :: stiff(nloc,nloc), force(nloc)
    call DefaultUpdateEquations(stiff, force)
  end subroutine cns_default_update

  subroutine cns_finish_bulk() bind(c, name="cns_finish_bulk")
    call DefaultFinishBulkAssembly()
  end subroutine cns_finish_bulk

  subroutine cns_finish_assembly() bind(c, name="cns_finish_assembly")
    call DefaultFinishAssembly()
  end subroutine cns_finish_assembly

  subroutine cns_dirichlet() bind(c, name="cns_dirichlet")
    call DefaultDirichletBCs()
  end subroutine cns_dirichlet

  ! -- per-nonlinear-iteration tail (mirrors the timing/Info/Solve block)
  function cns_iter_tail() result(norm) bind(c, name="cns_iter_tail")
    real(c_double) :: norm, stime
    acc_atime = CPUTime() - acc_atime
    stime = CPUTime()
    norm = DefaultSolve()
    stime = CPUTime() - stime
    write(message, '(a,F8.2)') ' Assembly:  (s)', acc_atime
    call Info('NavierStokesSolver', message, level=4)
    write(message, '(a,F8.2)') ' Solution:  (s)', stime
    call Info('NavierStokesSolver', message, level=4)
  end function cns_iter_tail

  function cns_nonlin_converged() result(r) &
      bind(c, name="cns_nonlin_converged")
    integer(c_int) :: r
    r = merge(1, 0, acc_solver % variable % nonlinconverged > 0)
  end function cns_nonlin_converged

  subroutine cns_default_finish() bind(c, name="cns_default_finish")
    call DefaultFinish()
  end subroutine cns_default_finish

  subroutine cns_fatal_alloc() bind(c, name="cns_fatal_alloc")
    call Fatal('NavierStokesSolver', 'Memory allocation error.')
  end subroutine cns_fatal_alloc

end module compressiblens_acc
