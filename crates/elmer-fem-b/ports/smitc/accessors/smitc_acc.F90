! T021-B Smitc port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/Smitc.F90, compiled against the built .mod
! files of the target configure (T020 Poisson pattern, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module mirrors Smitc's own locals (current element, ElementNodes SAVE
! state, the internal-procedure IP state, Indexes SAVE) with identical
! lifetimes.
!
! The whole LocalMatrix internal procedure (with its PerforatedElasticity
! and CovariantInterpolation blocks inlined, including the dead
! FMIN4 IF(.FALSE.) branch omitted as never-executed) lives in the
! sma_local_matrix wrapper, so every order-sensitive reduction (SUM,
! MATMUL via StrainEnergyDensity) and the IsotropicElasticity /
! ShearCorrectionFactor / Jacobi3 / Jacobi4 calls keep the same compiler
! and evaluation order as the original. Whole-array updates stay in
! Fortran wrappers, so bit-identity holds by construction wherever the
! call sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module smitc_acc
  use iso_c_binding
  use types
  use defutils
  use solidmechanicsutils, only: strainenergydensity, &
      shearcorrectionfactor, isotropicelasticity, jacobi3, jacobi4
  use messages
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(model_t), pointer, save :: acc_model => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(valuelist_t), pointer, save :: acc_bf => null()
  type(valuelist_t), pointer, save :: acc_mat => null()
  integer, pointer, save :: acc_indexes(:) => null()
  ! LocalMatrix context (mirrors its SAVE ElementNodes and IP local).
  type(nodes_t), target, save :: acc_nodes
  type(gaussintegrationpoints_t), save :: acc_ip

contains

  ! -- init0 entry ----------------------------------------------------
  subroutine sma_init0() bind(c, name="sma_init0")
    type(valuelist_t), pointer :: params
    params => getsolverparams()
    call listaddnewstring(params, 'Variable', '-dofs 3 Deflection')
    call listaddnewinteger(params, 'Time derivative order', 2)
    call listaddnewlogical(params, 'Use Global Mass Matrix', .true.)
    call listaddnewlogical(params, 'Plate Solver', .true.)
  end subroutine sma_init0

  ! -- entry state ----------------------------------------------------
  subroutine sma_begin(h_solver, h_model) bind(c, name="sma_begin")
    type(c_ptr), value :: h_solver, h_model
    call c_f_pointer(h_solver, acc_solver)
    call c_f_pointer(h_model, acc_model)
  end subroutine sma_begin

  function sma_dofs() result(r) bind(c, name="sma_dofs")
    integer(c_int) :: r
    r = acc_solver % variable % dofs
  end function sma_dofs

  function sma_max_elem_dofs() result(n) bind(c, name="sma_max_elem_dofs")
    integer(c_int) :: n
    n = acc_solver % mesh % maxelementdofs
  end function sma_max_elem_dofs

  function sma_nof_active() result(n) bind(c, name="sma_nof_active")
    integer(c_int) :: n
    n = getnofactive()
  end function sma_nof_active

  function sma_nonlin_converged() result(r) bind(c, name="sma_nonlin_converged")
    integer(c_int) :: r
    r = merge(1, 0, acc_solver % variable % nonlinconverged == 1)
  end function sma_nonlin_converged

  subroutine sma_alloc_indexes(nmax) bind(c, name="sma_alloc_indexes")
    integer(c_int), value :: nmax
    integer :: istat
    if (.not. associated(acc_indexes)) then
      allocate(acc_indexes(nmax), stat=istat)
      if (istat /= 0) call fatal('SmitcSolver', 'Memory allocation error')
    end if
  end subroutine sma_alloc_indexes

  ! -- solver params ---------------------------------------------------
  subroutine sma_relax_alpha(alpha, found) bind(c, name="sma_relax_alpha")
    real(c_double) :: alpha
    integer(c_int) :: found
    type(valuelist_t), pointer :: params
    logical :: f
    params => getsolverparams()
    alpha = listgetcreal(params, 'Shear Relaxation Alpha', f)
    found = merge(1, 0, f)
  end subroutine sma_relax_alpha

  function sma_harmonic_mode() result(r) bind(c, name="sma_harmonic_mode")
    integer(c_int) :: r
    type(valuelist_t), pointer :: params
    logical :: f, e
    params => getsolverparams()
    e = eigenorharmonicanalysis() .or. listgetlogical(params, 'Harmonic Mode', f)
    r = merge(1, 0, e)
  end function sma_harmonic_mode

  function sma_kernel_version() result(r) bind(c, name="sma_kernel_version")
    integer(c_int) :: r
    type(valuelist_t), pointer :: params
    logical :: f
    params => getsolverparams()
    r = merge(1, 0, listgetlogical(params, 'Kernel Interpolation', f))
  end function sma_kernel_version

  function sma_use_gravity() result(r) bind(c, name="sma_use_gravity")
    integer(c_int) :: r
    type(valuelist_t), pointer :: params
    logical :: f
    params => getsolverparams()
    r = merge(1, 0, listgetlogical(params, 'Use Gravity', f))
  end function sma_use_gravity

  subroutine sma_gravity_value(grav, found) bind(c, name="sma_gravity_value")
    real(c_double) :: grav
    integer(c_int) :: found
    real(c_double), pointer :: gwork(:,:)
    logical :: f
    gwork => listgetconstrealarray(currentmodel % constants, 'Gravity', f)
    if (f) grav = abs(gwork(size(gwork,1),1))
    found = merge(1, 0, f)
  end subroutine sma_gravity_value

  subroutine sma_mat_index(idx, got) bind(c, name="sma_mat_index")
    integer(c_int) :: idx, got
    type(valuelist_t), pointer :: params
    logical :: g
    params => getsolverparams()
    idx = listgetinteger(params, 'Material Index', g, &
        minv=1, maxv=acc_model % numberofmaterials)
    got = merge(1, 0, g)
  end subroutine sma_mat_index

  subroutine sma_bf_index(idx, got) bind(c, name="sma_bf_index")
    integer(c_int) :: idx, got
    type(valuelist_t), pointer :: params
    logical :: g
    params => getsolverparams()
    idx = listgetinteger(params, 'Body Force Index', g, &
        minv=1, maxv=acc_model % numberofbodyforces)
    got = merge(1, 0, g)
  end subroutine sma_bf_index

  subroutine sma_material_select(got, mat_id) bind(c, name="sma_material_select")
    integer(c_int), value :: got, mat_id
    acc_mat => null()
    if (got /= 0) acc_mat => acc_model % materials(mat_id) % values
  end subroutine sma_material_select

  subroutine sma_bodyforce_select(got, mat_id) bind(c, name="sma_bodyforce_select")
    integer(c_int), value :: got, mat_id
    ! Faithful to the original (which reads Materials(mat_id) here too).
    acc_bf => null()
    if (got /= 0) acc_bf => acc_model % materials(mat_id) % values
  end subroutine sma_bodyforce_select

  function sma_hole_correction() result(r) bind(c, name="sma_hole_correction")
    integer(c_int) :: r
    type(valuelist_t), pointer :: params
    logical :: g
    params => getsolverparams()
    r = merge(1, 0, listgetlogical(params, 'Hole Correction', g))
  end function sma_hole_correction

  subroutine sma_max_iter(maxiter, got) bind(c, name="sma_max_iter")
    integer(c_int) :: maxiter, got
    type(valuelist_t), pointer :: params
    logical :: g
    params => getsolverparams()
    maxiter = getinteger(params, 'Nonlinear System Max Iterations', g)
    got = merge(1, 0, g)
  end subroutine sma_max_iter

  ! -- driver calls -----------------------------------------------------
  subroutine sma_info_banner() bind(c, name="sma_info_banner")
    call info('SmitcSolver', '--------------------------------------------------', level=4)
    call info('SmitcSolver', 'Solving the Reissner-Mindlin equations for plates', level=4)
    call info('SmitcSolver', '--------------------------------------------------', level=4)
  end subroutine sma_info_banner

  subroutine sma_default_start() bind(c, name="sma_default_start")
    call defaultstart()
  end subroutine sma_default_start

  subroutine sma_default_initialize() bind(c, name="sma_default_initialize")
    call defaultinitialize()
  end subroutine sma_default_initialize

  function sma_cpu() result(r) bind(c, name="sma_cpu")
    real(c_double) :: r
    r = cputime()
  end function sma_cpu

  subroutine sma_info_assembly(at) bind(c, name="sma_info_assembly")
    real(c_double), value :: at
    write(message,*) 'Assembly (s): ', at
    call info('SmitcSolver', message, level=4)
  end subroutine sma_info_assembly

  subroutine sma_info_solve(st) bind(c, name="sma_info_solve")
    real(c_double), value :: st
    write(message,*) 'Solve (s): ', st
    call info('SmitcSolver', message, level=4)
  end subroutine sma_info_solve

  ! -- active-element context --------------------------------------------
  subroutine sma_begin_active(t, n_dofs, n_nodes) bind(c, name="sma_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n_dofs, n_nodes
    acc_elem => getactiveelement(t)
    n_dofs = getelementdofs(acc_indexes)
    n_nodes = getelementnofnodes()
    call getelementnodes(acc_nodes)
  end subroutine sma_begin_active

  function sma_bodyforce() result(r) bind(c, name="sma_bodyforce")
    integer(c_int) :: r
    acc_bf => getbodyforce()
    r = merge(1, 0, associated(acc_bf))
  end function sma_bodyforce

  function sma_material() result(r) bind(c, name="sma_material")
    integer(c_int) :: r
    acc_mat => getmaterial()
    r = merge(1, 0, associated(acc_mat))
  end function sma_material

  subroutine sma_load_pressure(nn, n, load, load2) bind(c, name="sma_load_pressure")
    integer(c_int), value :: nn, n
    real(c_double) :: load(nn), load2(nn)
    logical :: gotit
    ! Exact transcription of the Pressure / Pressure B / Pressure C block.
    if (associated(acc_bf)) then
      load(1:n) = getreal(acc_bf, 'Pressure', gotit)
      if (gotit) then
        load2(1:n) = getreal(acc_bf, 'Pressure B', gotit)
        if (gotit) then
          load(1:n) = load(1:n) + load2(1:n)
          load2(1:n) = getreal(acc_bf, 'Pressure C', gotit)
          if (gotit) load(1:n) = load(1:n) + load2(1:n)
        end if
      end if
    end if
  end subroutine sma_load_pressure

  subroutine sma_material_vecs(nn, n, poisson, density, thickness, young, tension) &
      bind(c, name="sma_material_vecs")
    integer(c_int), value :: nn, n
    real(c_double) :: poisson(nn), density(nn), thickness(nn), young(nn), tension(nn)
    logical :: gotit
    poisson(1:n) = getreal(acc_mat, 'Poisson ratio')
    density(1:n) = getreal(acc_mat, 'Density')
    thickness(1:n) = getreal(acc_mat, 'Thickness')
    young(1:n) = getreal(acc_mat, 'Youngs modulus')
    tension(1:n) = getreal(acc_mat, 'Tension', gotit)
  end subroutine sma_material_vecs

  subroutine sma_dampcoef(nn, n, dampcoef) bind(c, name="sma_dampcoef")
    integer(c_int), value :: nn, n
    real(c_double) :: dampcoef(nn)
    logical :: gotit
    ! Exact transcription of the Damping BodyForce-vs-Material fallback.
    if (associated(acc_bf)) then
      dampcoef(1:n) = getreal(acc_bf, 'Damping', gotit)
    else
      gotit = .false.
    end if
    if (.not. gotit) dampcoef(1:n) = getreal(acc_mat, 'Damping', gotit)
    if (.not. gotit) dampcoef(1:n) = 0.0d0
  end subroutine sma_dampcoef

  subroutine sma_springcoef(nn, n, springcoef) bind(c, name="sma_springcoef")
    integer(c_int), value :: nn, n
    real(c_double) :: springcoef(nn)
    logical :: gotit
    ! Exact transcription of the Spring BodyForce-vs-Material fallback.
    if (associated(acc_bf)) then
      springcoef(1:n) = getreal(acc_bf, 'Spring', gotit)
    else
      gotit = .false.
    end if
    if (.not. gotit) springcoef(1:n) = getreal(acc_mat, 'Spring', gotit)
    if (.not. gotit) springcoef(1:n) = 0.0d0
  end subroutine sma_springcoef

  function sma_hole(nn, n, holesize, holefraction) result(r) &
      bind(c, name="sma_hole")
    integer(c_int), value :: nn, n
    real(c_double) :: holesize(nn), holefraction(nn)
    integer(c_int) :: r
    character(len=MAX_NAME_LEN) :: holetype
    logical :: gotholetype
    ! Exact transcription of the Hole Type / Size / Fraction block.
    holetype = getstring(acc_mat, 'Hole Type', gotholetype)
    if (gotholetype) then
      holesize(1:n) = getreal(acc_mat, 'Hole Size')
      holefraction(1:n) = getreal(acc_mat, 'Hole Fraction')
    end if
    r = merge(1, 0, gotholetype)
  end function sma_hole

  ! -- LocalMatrix (+ PerforatedElasticity + CovariantInterpolation) ------
  ! Faithful inline transcription of the SmitcSolver internal procedures.
  ! The dead FMIN4 IF(.FALSE.) branch of the original is omitted (it never
  ! executes, so omitting it is semantics-preserving).
  subroutine sma_local_matrix(ld, nn, n, dofs, stiff, damp, mass, force, &
      load, poisson, young, thickness, tension, dampcoef, springcoef, &
      holefraction, holesize, density, relax, kernel, holecorr, gotholetype, grav) &
      bind(c, name="sma_local_matrix")
    integer(c_int), value :: ld, nn, n, dofs, kernel, holecorr, gotholetype
    real(c_double) :: stiff(ld,ld), damp(ld,ld), mass(ld,ld), force(ld)
    real(c_double) :: load(nn), poisson(nn), young(nn), thickness(nn), &
        tension(nn), dampcoef(nn), springcoef(nn), holefraction(nn), holesize(nn), &
        density(nn)
    real(c_double), value :: relax, grav
    ! Locals (mirrors the internal procedures' automatics).
    integer :: i, j, p, q, t, nd
    real(c_double) :: basis(n), dbasisdx(n,3), curvature(3,3*n), &
        shearstrain(2,3*n), ematrix(3,3), gmatrix(2,2), tmatrix(2,2)
    real(c_double) :: detj, u, v, w, s, kappa, rho, h, qeff
    real(c_double) :: pressure, dampc, winkc
    real(c_double) :: euvw, puvw, guvw, tuvw, qq, aa, bb, kk, sq
    logical :: stat
    ! CovariantInterpolation locals.
    real(c_double) :: jmat(2,2), invj(2,2), shearref(2,3*n)
    real(c_double) :: tau(2), sdofs(3*n)
    real(c_double) :: a1, a2, a12, b1, b2, b12
    real(c_double) :: abf0basis(2,6)
    integer :: jj

    nd = dofs * n

    force(1:nd) = 0.0d0
    stiff(1:nd,1:nd) = 0.0d0
    damp(1:nd,1:nd) = 0.0d0
    mass(1:nd,1:nd) = 0.0d0
    curvature = 0.0d0
    shearstrain = 0.0d0
!
!   Numerical integration:
!   ----------------------
    acc_ip = gausspoints(acc_elem)

    do t = 1, acc_ip % n
      u = acc_ip % u(t)
      v = acc_ip % v(t)
      w = acc_ip % w(t)
      s = acc_ip % s(t)

!
!     Basis function values & derivatives at the integration point
!     -------------------------------------------------------------
      stat = elementinfo(acc_elem, acc_nodes, u, v, w, detj, &
          basis(1:n), dbasisdx(1:n,1:3))

      s = s * detj

      pressure = sum(load(1:n)*basis(1:n))
      h = sum(thickness(1:n)*basis(1:n))
      dampc = sum(dampcoef(1:n)*basis(1:n))
      winkc = sum(springcoef(1:n)*basis(1:n))

      if (holecorr /= 0 .and. gotholetype /= 0) then
        ! PerforatedElasticity transcription.
        euvw = sum(young(1:n)*basis(1:n))
        puvw = sum(poisson(1:n)*basis(1:n))
        tuvw = sum(thickness(1:n)*basis(1:n))
        qq = sum(holefraction(1:n)*basis(1:n))
        aa = sum(holesize(1:n)*basis(1:n))
        sq = sqrt(qq)
        bb = 2*aa/sq

        if (tuvw > bb-2*aa) then
          kk = (tuvw-0.63d0*(bb-2*aa)) * (bb-2*aa)**3 / 3.0d0
        else
          kk = ((bb-2*aa)-0.63d0*tuvw) * tuvw**3 / 3.0d0
        end if

        ematrix = 0.0d0
        ematrix(1,1) = (1.0d0-sq)/(1.0d0-puvw**2) + 0.5d0*sq*(1.0d0-sq)**2
        ematrix(2,2) = ematrix(1,1)
        ematrix(1,2) = puvw*(1.0d0-sq)/(1.0d0-puvw**2)
        ematrix(2,1) = ematrix(1,2)
        ematrix(3,3) = 0.5d0*(1.0d0-sq)/(1.0d0+puvw) + &
            1.5d0*kk*sq*(1-sq)/(bb*(1.0d0+puvw)*tuvw**3)

        guvw = euvw*ematrix(3,3)

        ematrix = ematrix*euvw*(tuvw**3)/12.0d0

        gmatrix = 0.0d0
        gmatrix(1,1) = guvw*tuvw
        gmatrix(2,2) = gmatrix(1,1)

        qeff = sum(holefraction(1:n)*basis(1:n))
        rho = (1.0d0-qeff)*sum(density(1:n)*basis(1:n))
        tmatrix = 0.0d0
        tmatrix(1,1) = sqrt(1.0d0-qeff**2)*sum(tension(1:n)*basis(1:n))*h
        tmatrix(2,2) = tmatrix(1,1)
      else
        call isotropicelasticity(ematrix, gmatrix, &
            poisson(1:nn), young(1:nn), thickness(1:nn), basis(1:n), n)
        rho = sum(density(1:n)*basis(1:n))
        tmatrix = 0.0d0
        tmatrix(1,1) = sum(tension(1:n)*basis(1:n))*h
        tmatrix(2,2) = sum(tension(1:n)*basis(1:n))*h
      end if

!
!     The degrees-of-freedom are  (u_x, u_y, u_z, r_x, r_y, r_z)
!     where u_i is the displacement [m] and r_i the rotation [1]
!     ----------------------------------------------------------

!     Bending stiffness:
!     ------------------
      do p = 1, n
        curvature(1,3*p-1) = dbasisdx(p,1)
        curvature(2,3*p  ) = dbasisdx(p,2)
        curvature(3,3*p-1) = dbasisdx(p,2)
        curvature(3,3*p  ) = dbasisdx(p,1)
      end do

      call strainenergydensity(stiff, ematrix, curvature, 3, 3*n, s)

!     In-plane stiffness:
!     -------------------


!     Shear stiffness:
!     ----------------
      ! CovariantInterpolation transcription.
      select case(n)

      ! The SMITC3 element
      ! ==================
      case(3)
        call jacobi3(jmat, invj, detj, acc_nodes % x(1:n), acc_nodes % y(1:n))
        shearref = 0.0d0
        shearstrain = 0.0d0

!       Compute the shear-dofs for edge 12:
!       ===================================
        tau(1) = 1.0d0
        tau(2) = 0.0d0
        tau = (/ 1.0d0, 0.0d0/)

        sdofs = 0.0d0
        sdofs(2) = (tau(1)*jmat(1,1)+tau(2)*jmat(2,1))/2.0d0
        sdofs(3) = (tau(1)*jmat(1,2)+tau(2)*jmat(2,2))/2.0d0
        sdofs(5) = (tau(1)*jmat(1,1)+tau(2)*jmat(2,1))/2.0d0
        sdofs(6) = (tau(1)*jmat(1,2)+tau(2)*jmat(2,2))/2.0d0

        do jj = 1, 9
          shearref(1,jj) = shearref(1,jj) + (1.0d0-v)*sdofs(jj)
          shearref(2,jj) = shearref(2,jj) + (u)*sdofs(jj)
        end do

!       Compute the shear-dofs for edge 23:
!       ===================================
        tau(1) = -1.0d0/sqrt(2.0d0)
        tau(2) =  1.0d0/sqrt(2.0d0)

        sdofs = 0.0d0
        sdofs(5) = (tau(1)*jmat(1,1)+tau(2)*jmat(2,1))/sqrt(2.0d0)
        sdofs(6) = (tau(1)*jmat(1,2)+tau(2)*jmat(2,2))/sqrt(2.0d0)
        sdofs(8) = (tau(1)*jmat(1,1)+tau(2)*jmat(2,1))/sqrt(2.0d0)
        sdofs(9) = (tau(1)*jmat(1,2)+tau(2)*jmat(2,2))/sqrt(2.0d0)

        do jj = 1, 9
          shearref(1,jj) = shearref(1,jj) + (-v)*sdofs(jj)
          shearref(2,jj) = shearref(2,jj) + (u)*sdofs(jj)
        end do

!       Compute the shear-dofs for edge 31:
!       ===================================
        tau(1) =  0.0d0
        tau(2) = -1.0d0

        sdofs = 0.0d0
        sdofs(2) = (tau(1)*jmat(1,1)+tau(2)*jmat(2,1))/2.0d0
        sdofs(3) = (tau(1)*jmat(1,2)+tau(2)*jmat(2,2))/2.0d0
        sdofs(8) = (tau(1)*jmat(1,1)+tau(2)*jmat(2,1))/2.0d0
        sdofs(9) = (tau(1)*jmat(1,2)+tau(2)*jmat(2,2))/2.0d0

        do jj = 1, 9
          shearref(1,jj) = shearref(1,jj) + (-v)*sdofs(jj)
          shearref(2,jj) = shearref(2,jj) + (-1.0d0+u)*sdofs(jj)
        end do

!       Compute the final reduced shear strain
!       ======================================
        shearstrain(1:2,1:9) = -matmul(invj,shearref(1:2,1:9))

      case(4)
        shearref = 0.0d0
        shearstrain = 0.0d0

        if (kernel /= 0) then
          ! The kernel of curl-conforming ABF_0 can be generated in terms of
          ! basis functions associated with the edges. The fourth edge (14)
          ! DOF is set such that the interpolant is in the kernel space. An
          ! equivalent constraint is to enforce the vanishing circulation
          ! around the element boundary.
          ! -----------------------------------------------------------
          abf0basis(1,1) = (1.0d0-v)/4.0d0
          abf0basis(2,1) = 0.0d0
          abf0basis(1,2) = 0.0d0
          abf0basis(2,2) = (1.0d0+u)/4.0d0
          abf0basis(1,3) = (1.0d0+v)/4.0d0
          abf0basis(2,3) = 0.0d0
          abf0basis(1,4) = 0.0d0
          abf0basis(2,4) = (1.0d0-u)/4.0d0
          ! The coefficients of the element map written as
          ! x(u,v) := a12*u*v + a1*u + a2*v + a0
          ! y(u,v) := b12*u*v + b1*u + b2*v + b0
          !
          a1 = 0.25d0*(-acc_nodes % x(1)+acc_nodes % x(2) &
              +acc_nodes % x(3)-acc_nodes % x(4))
          a2 = 0.25d0*(-acc_nodes % x(1)-acc_nodes % x(2) &
              +acc_nodes % x(3)+acc_nodes % x(4))
          a12 = 0.25d0*(acc_nodes % x(1)-acc_nodes % x(2) &
              +acc_nodes % x(3)-acc_nodes % x(4))
          b1 = 0.25d0*(-acc_nodes % y(1)+acc_nodes % y(2) &
              +acc_nodes % y(3)-acc_nodes % y(4))
          b2 = 0.25d0*(-acc_nodes % y(1)-acc_nodes % y(2) &
              +acc_nodes % y(3)+acc_nodes % y(4))
          b12 = 0.25d0*(acc_nodes % y(1)-acc_nodes % y(2) &
              +acc_nodes % y(3)-acc_nodes % y(4))
          ! The negative gradient of the deflection (no approximation
          ! introduced):
          shearref(:,1) = shearref(:,1) + abf0basis(:,1)
          shearref(:,4) = shearref(:,4) - abf0basis(:,1)
          shearref(:,4) = shearref(:,4) + abf0basis(:,2)
          shearref(:,7) = shearref(:,7) - abf0basis(:,2)
          shearref(:,7) = shearref(:,7) - abf0basis(:,3)
          shearref(:,10) = shearref(:,10) + abf0basis(:,3)
          shearref(:,10) = shearref(:,10) - abf0basis(:,4)
          shearref(:,1) = shearref(:,1) + abf0basis(:,4)
          ! Edge 12:
          ! ----------
          shearref(:,2) = shearref(:,2) + (a1-a12)*abf0basis(:,1)
          shearref(:,3) = shearref(:,3) + (b1-b12)*abf0basis(:,1)
          shearref(:,5) = shearref(:,5) + (a1-a12)*abf0basis(:,1)
          shearref(:,6) = shearref(:,6) + (b1-b12)*abf0basis(:,1)
          ! Edge 23:
          ! -----------
          shearref(:,5) = shearref(:,5) + (a2+a12)*abf0basis(:,2)
          shearref(:,6) = shearref(:,6) + (b2+b12)*abf0basis(:,2)
          shearref(:,8) = shearref(:,8) + (a2+a12)*abf0basis(:,2)
          shearref(:,9) = shearref(:,9) + (b2+b12)*abf0basis(:,2)
          ! Edge 43:
          ! -----------
          shearref(:,8) = shearref(:,8) + (a1+a12)*abf0basis(:,3)
          shearref(:,9) = shearref(:,9) + (b1+b12)*abf0basis(:,3)
          shearref(:,11) = shearref(:,11) + (a1+a12)*abf0basis(:,3)
          shearref(:,12) = shearref(:,12) + (b1+b12)*abf0basis(:,3)
          ! Edge 14:
          ! ------------
          shearref(:,2) = shearref(:,2) + (a1-a12)*abf0basis(:,4)
          shearref(:,3) = shearref(:,3) + (b1-b12)*abf0basis(:,4)
          shearref(:,5) = shearref(:,5) + (a1+a2)*abf0basis(:,4)
          shearref(:,6) = shearref(:,6) + (b1+b2)*abf0basis(:,4)
          shearref(:,8) = shearref(:,8) + (a2-a1)*abf0basis(:,4)
          shearref(:,9) = shearref(:,9) + (b2-b1)*abf0basis(:,4)
          shearref(:,11) = shearref(:,11) + (-a1-a12)*abf0basis(:,4)
          shearref(:,12) = shearref(:,12) + (-b1-b12)*abf0basis(:,4)
        else

!         The SMITC4 element
!         ==================

!         Compute the shear-dofs for edge 12:
!         ===================================
          tau(1) = 1.0d0
          tau(2) = 0.0d0

          call jacobi4(jmat, invj, detj, 0.0d0, -1.0d0, &
              acc_nodes % x(1:n), acc_nodes % y(1:n))

          sdofs = 0.0d0
          sdofs(2) = (tau(1)*jmat(1,1)+tau(2)*jmat(2,1))
          sdofs(3) = (tau(1)*jmat(1,2)+tau(2)*jmat(2,2))
          sdofs(5) = (tau(1)*jmat(1,1)+tau(2)*jmat(2,1))
          sdofs(6) = (tau(1)*jmat(1,2)+tau(2)*jmat(2,2))

          do jj = 1, 12
            shearref(1,jj) = shearref(1,jj) + (1-v)/4.0d0*sdofs(jj)
          end do

!         Compute the shear-dofs for edge 23:
!         ===================================
          tau(1) = 0.0d0
          tau(2) = 1.0d0

          call jacobi4(jmat, invj, detj, 1.0d0, 0.0d0, &
              acc_nodes % x(1:n), acc_nodes % y(1:n))

          sdofs = 0.0d0
          sdofs(5) = (tau(1)*jmat(1,1)+tau(2)*jmat(2,1))
          sdofs(6) = (tau(1)*jmat(1,2)+tau(2)*jmat(2,2))
          sdofs(8) = (tau(1)*jmat(1,1)+tau(2)*jmat(2,1))
          sdofs(9) = (tau(1)*jmat(1,2)+tau(2)*jmat(2,2))

          do jj = 1, 12
            shearref(2,jj) = shearref(2,jj) + (1+u)/4.0d0*sdofs(jj)
          end do

!         Compute the shear-dofs for edge 34:
!         ===================================
          tau(1) = -1.0d0
          tau(2) =  0.0d0

          call jacobi4(jmat, invj, detj, 0.0d0, 1.0d0, &
              acc_nodes % x(1:n), acc_nodes % y(1:n))

          sdofs = 0.0d0
          sdofs(8)  = (tau(1)*jmat(1,1)+tau(2)*jmat(2,1))
          sdofs(9)  = (tau(1)*jmat(1,2)+tau(2)*jmat(2,2))
          sdofs(11) = (tau(1)*jmat(1,1)+tau(2)*jmat(2,1))
          sdofs(12) = (tau(1)*jmat(1,2)+tau(2)*jmat(2,2))

          do jj = 1, 12
            shearref(1,jj) = shearref(1,jj) + (-1-v)/4.0d0*sdofs(jj)
          end do

!         Compute the shear-dofs for edge 41:
!         ===================================
          tau(1) =  0.0d0
          tau(2) = -1.0d0

          call jacobi4(jmat, invj, detj, -1.0d0, 0.0d0, &
              acc_nodes % x(1:n), acc_nodes % y(1:n))

          sdofs = 0.0d0
          sdofs(2)  = (tau(1)*jmat(1,1)+tau(2)*jmat(2,1))
          sdofs(3)  = (tau(1)*jmat(1,2)+tau(2)*jmat(2,2))
          sdofs(11) = (tau(1)*jmat(1,1)+tau(2)*jmat(2,1))
          sdofs(12) = (tau(1)*jmat(1,2)+tau(2)*jmat(2,2))

          do jj = 1, 12
            shearref(2,jj) = shearref(2,jj) + (-1+u)/4.0d0*sdofs(jj)
          end do

        end if

!       Compute the final reduced shear strain
!       ======================================
        call jacobi4(jmat, invj, detj, u, v, &
            acc_nodes % x(1:n), acc_nodes % y(1:n))
        shearstrain(1:2,1:12) = -matmul(invj,shearref(1:2,1:12))
      case default
        call warn('SmitcSolver', 'Illegal number of nodes for Smitc elements.')
      end select

      call shearcorrectionfactor(kappa, h, &
          acc_nodes % x(1:n), acc_nodes % y(1:n), n, relax)

      if (kernel == 0) then
        do p = 1, n
          shearstrain(1,3*p-2) = dbasisdx(p,1)
          shearstrain(2,3*p-2) = dbasisdx(p,2)
        end do
      end if

      call strainenergydensity(stiff, gmatrix, shearstrain, 2, 3*n, kappa*s)
!
!     Tensile stiffness:
!     ------------------
      shearstrain = 0.0d0
      do p = 1, n
        shearstrain(1,3*p-2) = dbasisdx(p,1)
        shearstrain(2,3*p-2) = dbasisdx(p,2)
      end do

      call strainenergydensity(stiff, tmatrix, shearstrain, 2, 3*n, s)

!     Spring Coeffficient:
!     -------------------
      do p = 1, n
        i = dofs*(p-1)+1
        do q = 1, n
          j = dofs*(q-1)+1
          stiff(i,j) = stiff(i,j) + winkc*basis(p)*basis(q)*s
        end do
      end do

!
!     Load vector:
!     ------------
      do p = 1, n
        i = dofs*(p-1)+1
        force(i) = force(i) + pressure*basis(p)*s
      end do

!     Gravity in direction down:
!     -------------------------
      do p = 1, n
        i = dofs*(p-1)+1
        force(i) = force(i) - grav*rho*h*basis(p)*s
      end do
!
!     Mass matrix:
!     ------------
      do p = 1, n
        i = dofs*(p-1)+1
        do q = 1, n
          j = dofs*(q-1)+1
          mass(i,j) = mass(i,j) + rho*h*basis(p)*basis(q)*s
        end do
      end do

!
!     Damping matrix:
!     ---------------
      do p = 1, n
        i = dofs*(p-1)+1
        do q = 1, n
          j = dofs*(q-1)+1
          damp(i,j) = damp(i,j) + dampc*basis(p)*basis(q)*s
        end do
      end do

    end do
  end subroutine sma_local_matrix

  ! -- shared assembly tail ----------------------------------------------
  subroutine sma_default_2nd_order(ld, mass, damp, stiff, force) &
      bind(c, name="sma_default_2nd_order")
    integer(c_int), value :: ld
    real(c_double) :: mass(ld,ld), damp(ld,ld), stiff(ld,ld), force(ld)
    call default2ndordertime(mass(1:ld,1:ld), damp(1:ld,1:ld), &
        stiff(1:ld,1:ld), force(1:ld))
  end subroutine sma_default_2nd_order

  subroutine sma_default_update(ld, g, f) bind(c, name="sma_default_update")
    integer(c_int), value :: ld
    real(c_double) :: g(ld, ld), f(ld)
    call defaultupdateequations(g(1:ld, 1:ld), f(1:ld))
  end subroutine sma_default_update

  subroutine sma_default_update_mass(ld, m) bind(c, name="sma_default_update_mass")
    integer(c_int), value :: ld
    real(c_double) :: m(ld, ld)
    call defaultupdatemass(m(1:ld, 1:ld))
  end subroutine sma_default_update_mass

  subroutine sma_default_update_damp(ld, b) bind(c, name="sma_default_update_damp")
    integer(c_int), value :: ld
    real(c_double) :: b(ld, ld)
    call defaultupdatedamp(b(1:ld, 1:ld))
  end subroutine sma_default_update_damp

  subroutine sma_default_finish_bulk() bind(c, name="sma_default_finish_bulk")
    call defaultfinishbulkassembly()
  end subroutine sma_default_finish_bulk

  subroutine sma_default_finish_boundary() bind(c, name="sma_default_finish_boundary")
    call defaultfinishboundaryassembly()
  end subroutine sma_default_finish_boundary

  subroutine sma_default_finish_assembly() bind(c, name="sma_default_finish_assembly")
    call defaultfinishassembly()
  end subroutine sma_default_finish_assembly

  subroutine sma_default_dirichlet() bind(c, name="sma_default_dirichlet")
    call defaultdirichletbcs()
  end subroutine sma_default_dirichlet

  function sma_default_solve() result(norm) bind(c, name="sma_default_solve")
    real(c_double) :: norm
    norm = defaultsolve()
  end function sma_default_solve

  subroutine sma_default_finish() bind(c, name="sma_default_finish")
    call defaultfinish()
  end subroutine sma_default_finish

  subroutine sma_fatal_alloc() bind(c, name="sma_fatal_alloc")
    call fatal('SmitcSolver', 'Memory allocation error')
  end subroutine sma_fatal_alloc

  subroutine sma_fatal_gravity() bind(c, name="sma_fatal_gravity")
    call fatal('SmitcSolver', 'Gravity requested but not given as constant!')
  end subroutine sma_fatal_gravity

end module smitc_acc
