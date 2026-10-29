! T022-A FDiffusion3D port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/FDiffusion3D.F90, compiled against the built .mod
! files of the target configure (A2 shipped mechanism, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (TYPE(C_PTR) handles,
! VALUE scalars). All Fortran-object traffic — TYPE(...) handles,
! OPTIONAL/POINTER dummies, hidden lengths, assumed-shape descriptor
! dummies, list-directed timing reports — stays inside Fortran, so the Rust
! side is pure C calls. The wrapper module mirrors the solver's own SAVE
! locals (entry SAVE arrays + flags, internal-procedure automatic arrays as
! per-call locals) with identical lifetimes. All floating-point work
! (complex 3-D Galerkin assembly, SUM reductions, the CPUTime/RealTime-based
! reports) stays in Fortran (same compiler), so bit-identity holds by
! construction wherever the call sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module fdiffusion3d_acc
  use iso_c_binding
  use types
  use defutils
  use lists
  use integration
  use eleminfo
  use solverbasics
  use coordinatesystems
  use messages
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(matrix_t), pointer, save :: acc_mat => null()
  type(valuelist_t), pointer, save :: acc_material => null()
  real(c_double), pointer, save :: acc_values(:) => null()
  integer(c_int), pointer, save :: acc_perm(:) => null()
  real(c_double), pointer, save :: acc_rhs(:) => null()

  ! Entry SAVE state (mirrors FourierDiffusion3DSolver's SAVE statement).
  real(c_double), pointer, save :: acc_lrf(:) => null()
  real(c_double), allocatable, save :: acc_lstiff(:,:), acc_load(:,:)
  real(c_double), allocatable, save :: acc_force(:), acc_work(:)
  real(c_double), allocatable, save :: acc_imp(:,:), acc_cond(:), acc_pm(:)
  type(nodes_t), save :: acc_nodes
  integer(c_int), save :: acc_localnodes = 0
  logical, save :: acc_allocdone = .false.
  logical, save :: acc_bubbles = .false.
  character(len=MAX_NAME_LEN), save :: acc_eqname = ''
  real(c_double), save :: acc_omega = 0.0_c_double
  real(c_double), save :: acc_at = 0.0_c_double, acc_at0 = 0.0_c_double
  real(c_double), save :: acc_totat = 0.0_c_double, acc_totst = 0.0_c_double

contains

  ! -- entry state ----------------------------------------------------
  function fd3_begin(h_model, h_solver) result(code) bind(c, name="fd3_begin")
    type(c_ptr), value :: h_model, h_solver
    integer(c_int) :: code
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    if (.not. associated(acc_solver % matrix)) then
      code = 1
      return
    end if
    acc_solver % matrix % complex = .true.
    acc_values => acc_solver % variable % values
    acc_perm => acc_solver % variable % perm
    acc_localnodes = count(acc_perm > 0)
    if (acc_localnodes <= 0) then
      code = 2
      return
    end if
    acc_mat => acc_solver % matrix
    acc_rhs => acc_mat % rhs
    code = 0
  end function fd3_begin

  subroutine fd3_alloc() bind(c, name="fd3_alloc")
    type(variable_t), pointer :: var
    type(solver_t), pointer :: sptr
    integer :: n, istat
    if (.not. acc_allocdone .or. acc_solver % meshchanged) then
      n = acc_solver % mesh % maxelementnodes
      if (acc_allocdone) then
        deallocate(acc_nodes % x, acc_nodes % y, acc_nodes % z, &
            acc_force, acc_imp, acc_work, acc_lstiff, acc_cond, acc_pm, &
            acc_load)
      end if
      allocate(acc_nodes % x(n), acc_nodes % y(n), acc_nodes % z(n), &
          acc_force(6*n), acc_imp(2,n), acc_work(n), &
          acc_lstiff(6*n,6*n), acc_cond(n), acc_pm(n), &
          acc_lrf(3*acc_model % numberofnodes), acc_load(6,n), stat=istat)
      if (istat /= 0) then
        call fatal('FourierDiffusion3DSolve', 'Memory allocation error.')
      end if
      var => variableget(acc_model % variables, 'Lorentz Force')
      if (.not. associated(var)) then
        sptr => acc_solver
        call variableaddvector(acc_solver % mesh % variables, &
            acc_solver % mesh, sptr, 'Lorentz Force', 3, acc_lrf, acc_perm)
      end if
      acc_allocdone = .true.
    end if
  end subroutine fd3_alloc

  function fd3_params() result(nl) bind(c, name="fd3_params")
    integer(c_int) :: nl
    logical :: gotit
    nl = listgetinteger(acc_solver % values, &
        'Nonlinear System Max Iterations', gotit)
    if (.not. gotit) nl = 1
    acc_eqname = listgetstring(acc_solver % values, 'Equation')
    acc_bubbles = listgetlogical(acc_solver % values, 'Bubbles', gotit)
    if (.not. gotit) acc_bubbles = .false.
    acc_omega = getangularfrequency()
  end function fd3_params

  function fd3_nbulk() result(n) bind(c, name="fd3_nbulk")
    integer(c_int) :: n
    n = acc_solver % mesh % numberofbulkelements
  end function fd3_nbulk

  ! -- nonlinear iteration head (Info banner + DefaultInitialize) -------
  subroutine fd3_iter_head(iter) bind(c, name="fd3_iter_head")
    integer(c_int), value :: iter
    acc_at = cputime()
    acc_at0 = realtime()
    call info('FourierDiffusion3DSolve', ' ', level=4)
    call info('FourierDiffusion3DSolve', ' ', level=4)
    call info('FourierDiffusion3DSolve', &
        '-------------------------------------', level=4)
    write(message, *) 'Fourier Diffusion 3D iteration', iter
    call info('FourierDiffusion3DSolve', message, level=4)
    write(message, *) 'Frequency (Hz): ', acc_omega / (2*pi)
    call info('FourierDiffusion3DSolve', message, level=4)
    call info('FourierDiffusion3DSolve', &
        '-------------------------------------', level=4)
    call info('FourierDiffusion3DSolve', ' ', level=4)
    call info('FourierDiffusion3DSolve', 'Starting Assembly...', level=4)
    call defaultinitialize()
  end subroutine fd3_iter_head

  ! -- bulk element body (exact transcription of the t-loop body) -------
  ! NOTE: the original assembles bulk elements only; there is no boundary
  ! loop (LocalMatrixBoundary below is dead code there and here).
  subroutine fd3_bulk(t) bind(c, name="fd3_bulk")
    integer(c_int), value :: t
    type(element_t), pointer :: cur
    integer, pointer :: nidx(:)
    integer :: n, k
    logical :: gotit
    cur => acc_solver % mesh % elements(t)
    if (realtime() - acc_at0 > 1.0d0) then
      write(*,'(a,i3,a)') '   Assembly: ', int(100.0d0 - 100.0d0 * &
          (acc_solver % mesh % numberofbulkelements - t) / &
          (1.0d0 * acc_solver % mesh % numberofbulkelements)), ' % done'
      acc_at0 = realtime()
    end if
    if (.not. checkelementequation(acc_model, cur, acc_eqname)) return
    acc_model % currentelement => cur
    acc_elem => cur
    n = cur % type % numberofnodes
    nidx => cur % nodeindexes
    acc_nodes % x(1:n) = acc_solver % mesh % nodes % x(nidx)
    acc_nodes % y(1:n) = acc_solver % mesh % nodes % y(nidx)
    acc_nodes % z(1:n) = acc_solver % mesh % nodes % z(nidx)
    k = listgetinteger(acc_model % bodies(cur % bodyid) % values, &
        'Equation', minv=1, maxv=acc_model % numberofequations)
    acc_work(1:n) = listgetreal(acc_model % equations(k) % values, &
        'Angular Frequency', n, nidx, gotit)
    if (gotit) then
      acc_omega = acc_work(1)
    else
      acc_work(1:n) = listgetreal(acc_model % equations(k) % values, &
          'Frequency', n, nidx, gotit)
      if (gotit) acc_omega = 2*pi*acc_work(1)
    end if
    k = listgetinteger(acc_model % bodies(cur % bodyid) % values, &
        'Material', minv=1, maxv=acc_model % numberofmaterials)
    acc_material => acc_model % materials(k) % values
    acc_cond(1:n) = listgetreal(acc_material, &
        'Electrical Conductivity', n, nidx, gotit)
    if (gotit) then
      call warn('Fdiffusion3D', &
          'Use electric conductivity instead of electrical')
    end if
    acc_cond(1:n) = listgetreal(acc_material, &
        'Electric Conductivity', n, nidx)
    acc_pm(1:n) = listgetreal(acc_material, &
        'Magnetic Permeability', n, nidx)
    k = listgetinteger(acc_model % bodies(cur % bodyid) % values, &
        'Body Force', gotit, 1, acc_model % numberofbodyforces)
    acc_load = 0.0d0
    if (k > 0) then
      acc_load(1,1:n) = listgetreal(acc_model % bodyforces(k) % values, &
          'AC Magnetic Field 1 Re', n, nidx, gotit)
      acc_load(2,1:n) = listgetreal(acc_model % bodyforces(k) % values, &
          'AC Magnetic Field 1 Im', n, nidx, gotit)
      acc_load(3,1:n) = listgetreal(acc_model % bodyforces(k) % values, &
          'AC Magnetic Field 2 Re', n, nidx, gotit)
      acc_load(4,1:n) = listgetreal(acc_model % bodyforces(k) % values, &
          'AC Magnetic Field 2 Im', n, nidx, gotit)
      acc_load(5,1:n) = listgetreal(acc_model % bodyforces(k) % values, &
          'AC Magnetic Field 3 Re', n, nidx, gotit)
      acc_load(6,1:n) = listgetreal(acc_model % bodyforces(k) % values, &
          'AC Magnetic Field 3 Im', n, nidx, gotit)
    end if
    call fd3_localmatrix(acc_lstiff, acc_force, acc_omega, &
        acc_cond(1:n)*acc_pm(1:n), acc_load, acc_bubbles, acc_elem, n, &
        acc_nodes)
    call updateglobalequations(acc_mat, acc_lstiff, acc_rhs, acc_force, &
        n, acc_solver % variable % dofs, acc_perm(nidx))
  end subroutine fd3_bulk

  ! -- nonlinear iteration tail (finish, solve, timing reports) ---------
  function fd3_iter_tail(iter) result(conv) bind(c, name="fd3_iter_tail")
    integer(c_int), value :: iter
    integer(c_int) :: conv
    real(c_double) :: st, norm, relchange
    call defaultfinishassembly()
    call defaultdirichletbcs()
    call info('FourierDiffusion3DSolve', 'Assembly done', level=4)
    acc_at = cputime() - acc_at
    st = cputime()
    norm = defaultsolve()
    relchange = acc_solver % variable % nonlinchange
    st = cputime() - st
    acc_totat = acc_totat + acc_at
    acc_totst = acc_totst + st
    write(message,'(a,i4,a,F8.2,F8.2)') 'iter: ', iter, &
        ' Assembly: (s)', acc_at, acc_totat
    call info('FourierDiffusion3DSolve', message, level=4)
    write(message,'(a,i4,a,F8.2,F8.2)') 'iter: ', iter, &
        ' Solve:    (s)', st, acc_totst
    call info('FourierDiffusion3DSolve', message, level=4)
    write(message, *) 'Result Norm    : ', norm
    call info('FourierDiffusion3DSolve', message, level=4)
    write(message, *) 'Relative Change: ', relchange
    call info('FourierDiffusion3DSolve', message, level=4)
    if (acc_solver % variable % nonlinconverged == 1) then
      conv = 1
    else
      conv = 0
    end if
  end function fd3_iter_tail

  ! -- LocalMatrix internal procedure (exact transcription) --------------
  subroutine fd3_localmatrix(stiffmatrix, force, angularfrequency, &
      nodalconductivity, load, bubbles, element, n, nodes)
    real(c_double) :: stiffmatrix(:,:), force(:), angularfrequency, &
        nodalconductivity(:), load(:,:)
    logical :: bubbles
    integer :: n
    type(nodes_t) :: nodes
    type(element_t), target :: element
    real(c_double) :: basis(2*n), dbasisdx(2*n,3)
    real(c_double) :: sqrtelem, u, v, w, s, conductivity, m, d, lre(3), lim(3)
    complex(kind=dp) :: lstiff(3*n,3*n), lforce(3*n), a(3,3)
    logical :: stat
    integer :: i, j, k, l, h, p, q, t, dim, nbasis, coordsys
    type(gaussintegrationpoints_t) :: ip
    real(c_double) :: x, y, z, metric(3,3), sqrtmetric, symb(3,3,3), &
        dsymb(3,3,3,3)
    dim = 3
    coordsys = currentcoordinatesystem()
    ! Author's evident intent (cf. CoordinateSystemInfo's Cartesian branch,
    ! which sets Symbols = 0, dSymbols = 0): the original only initializes
    ! Metric here and skips CoordinateSystemInfo on the Cartesian path, so
    ! the reference reads uninitialized Symb/dSymb there (nondeterministic
    ! across processes; see proof header). Zero them: bit-identical with
    ! the reference wherever the reference is defined (non-Cartesian calls
    ! CoordinateSystemInfo, which overwrites all four).
    metric = 0.0d0
    metric(1,1) = 1.0d0
    metric(2,2) = 1.0d0
    metric(3,3) = 1.0d0
    symb = 0.0d0
    dsymb = 0.0d0
    lstiff = 0.0d0
    lforce = 0.0d0
    a = 0.0d0
    if (bubbles) then
      print *, 'FourierDiffusion3DSolver: LocalMatrix: Cannot handle bubbles.'
      stop
      ip = gausspoints(element, element % type % gausspoints2)
      nbasis = 2*n
    else
      nbasis = n
      ip = gausspoints(element)
    end if
    do t = 1, ip % n
      u = ip % u(t)
      v = ip % v(t)
      w = ip % w(t)
      s = ip % s(t)
      stat = elementinfo(element, nodes, u, v, w, sqrtelem, &
          basis, dbasisdx, bubbles=bubbles)
      s = s * sqrtelem
      if (coordsys /= cartesian) then
        x = sum(nodes % x(1:n) * basis(1:n))
        y = sum(nodes % y(1:n) * basis(1:n))
        z = sum(nodes % z(1:n) * basis(1:n))
        call coordinatesysteminfo(metric, sqrtmetric, symb, dsymb, x, y, z)
        s = s * sqrtmetric
      end if
      conductivity = sum(nodalconductivity(1:n) * basis)
      d = -angularfrequency
      m = 0._dp
      do i = 1, dim
        lre(i) = angularfrequency * sum(load(2*i-1,1:n) * basis(1:n))
        lim(i) = angularfrequency * sum(load(2*i,1:n) * basis(1:n))
      end do
      do p = 1, nbasis
        do q = 1, nbasis
          do i = 1, dim
            a(i,i) = cmplx(m, d, kind=dp) * basis(q) * basis(p)
          end do
          do i = 1, dim
            do j = 1, dim
              do k = 1, dim
                a(i,i) = a(i,i) + metric(j,k) * dbasisdx(q,k) * &
                    dbasisdx(p,j) / conductivity
                do l = 1, dim
                  a(l,i) = a(l,i) - metric(j,k) * dbasisdx(q,k) * &
                      symb(i,j,l) * basis(p) / conductivity
                  a(i,l) = a(i,l) + metric(j,k) * symb(l,k,i) * &
                      basis(q) * dbasisdx(p,j) / conductivity
                  do h = 1, dim
                    a(h,l) = a(h,l) - metric(j,k) * basis(q) * &
                        basis(p) * symb(l,k,i) * symb(i,j,h) / &
                        conductivity
                  end do
                end do
              end do
            end do
          end do
          do i = 1, dim
            do j = 1, dim
              lstiff(dim*(p-1)+i,dim*(q-1)+j) = &
                  lstiff(dim*(p-1)+i,dim*(q-1)+j) + s*a(i,j)
            end do
          end do
        end do
        do i = 1, dim
          lforce(dim*(p-1)+i) = lforce(dim*(p-1)+i) + s * basis(p) * &
              cmplx(lre(i), lim(i), kind=dp)
        end do
      end do
    end do
    do p = 1, n
      do i = 1, dim
        force(2*dim*(p-1)+2*i-1) = real(lforce(dim*(p-1)+i))
        force(2*dim*(p-1)+2*i) = aimag(lforce(dim*(p-1)+i))
        do q = 1, n
          do j = 1, dim
            stiffmatrix(2*dim*(p-1)+2*i-1, 2*dim*(q-1)+2*j-1) = &
                real(lstiff(dim*(p-1)+i,dim*(q-1)+j))
            stiffmatrix(2*dim*(p-1)+2*i-1, 2*dim*(q-1)+2*j) = &
                -aimag(lstiff(dim*(p-1)+i,dim*(q-1)+j))
            stiffmatrix(2*dim*(p-1)+2*i, 2*dim*(q-1)+2*j-1) = &
                aimag(lstiff(dim*(p-1)+i,dim*(q-1)+j))
            stiffmatrix(2*dim*(p-1)+2*i, 2*dim*(q-1)+2*j) = &
                real(lstiff(dim*(p-1)+i,dim*(q-1)+j))
          end do
        end do
      end do
    end do
  end subroutine fd3_localmatrix

  ! -- LocalMatrixBoundary internal procedure (exact transcription) -----
  ! NOTE: dead code in the original (never called from the bulk-only
  ! assembly loop); kept here uncalled for statement fidelity.
  subroutine fd3_localmatrixboundary(stiffmatrix, force, angularfrequency, &
      impedance, load, element, n, nodes)
    real(c_double) :: stiffmatrix(:,:), force(:), impedance(:,:), load(:,:)
    real(c_double) :: angularfrequency
    integer :: n
    type(nodes_t) :: nodes
    type(element_t), target :: element
    real(c_double) :: sqrtelem, u, v, w, s, impedance1, impedance2, l1, l2
    real(c_double) :: basis(n), dbasisdx(n,3), x, y, z
    complex(kind=dp) :: lstiff(n,n), lforce(n), a
    logical :: stat
    integer :: i, j, p, q, t, dim, coordsys
    type(gaussintegrationpoints_t) :: ip
    dim = coordinatesystemdimension()
    coordsys = currentcoordinatesystem()
    lstiff = 0.0d0
    lforce = 0.0d0
    ip = gausspoints(element)
    do t = 1, ip % n
      u = ip % u(t)
      v = ip % v(t)
      w = ip % w(t)
      s = ip % s(t)
      stat = elementinfo(element, nodes, u, v, w, sqrtelem, &
          basis, dbasisdx)
      s = s * sqrtelem
      if (coordsys /= cartesian) then
        x = sum(nodes % x(1:n) * basis(1:n))
        y = sum(nodes % y(1:n) * basis(1:n))
        z = sum(nodes % z(1:n) * basis(1:n))
        s = s * coordinatesqrtmetric(x, y, z)
      end if
      impedance1 = sum(impedance(1,1:n) * basis)
      if (abs(impedance1) > aeps) &
          impedance1 = angularfrequency / impedance1
      impedance2 = sum(impedance(2,1:n) * basis)
      if (abs(impedance2) > aeps) &
          impedance2 = angularfrequency / impedance2
      l1 = sum(load(1,1:n) * basis)
      l2 = sum(load(2,1:n) * basis)
      do p = 1, n
        do q = 1, n
          a = cmplx(impedance1, impedance2, kind=dp) * basis(q) * basis(p)
          lstiff(p,q) = lstiff(p,q) + s * a
        end do
        lforce(p) = lforce(p) + s * basis(p) * cmplx(l1, l2, kind=dp)
      end do
    end do
    do i = 1, n
      force(2*(i-1)+1) = real(lforce(i))
      force(2*(i-1)+2) = aimag(lforce(i))
      do j = 1, n
        stiffmatrix(2*(i-1)+1, 2*(j-1)+1) = real(lstiff(i,j))
        stiffmatrix(2*(i-1)+1, 2*(j-1)+2) = -aimag(lstiff(i,j))
        stiffmatrix(2*(i-1)+2, 2*(j-1)+1) = aimag(lstiff(i,j))
        stiffmatrix(2*(i-1)+2, 2*(j-1)+2) = real(lstiff(i,j))
      end do
    end do
  end subroutine fd3_localmatrixboundary

  ! -- LorentzForceAve internal procedure (exact transcription) ----------
  subroutine fd3_lorentz() bind(c, name="fd3_lorentz")
    real(c_double), pointer :: lrfx(:), lrfy(:), lrfz(:)
    real(c_double), pointer :: brex(:), bimx(:), brey(:), bimy(:), &
        brez(:), bimz(:)
    type(element_t), pointer :: el
    type(nodes_t) :: nodes
    logical :: stat, gotit
    integer, pointer :: nidx(:), visited(:)
    integer :: p, q, i, t, n, k
    real(c_double) :: u, v, w, x, y, z
    real(c_double) :: b(3), dhdx(3,3), l(3)
    real(c_double) :: sqrtelem
    real(c_double) :: sqrtmetric, metric(3,3), symb(3,3,3), dsymb(3,3,3,3)
    real(c_double), allocatable :: basis(:), dbasisdx(:,:)
    real(c_double), allocatable :: permeab(:), bacre(:,:), bacim(:,:)
    lrfx => acc_lrf(1::3)
    lrfy => acc_lrf(2::3)
    lrfz => acc_lrf(3::3)
    brex => acc_values(1::6)
    bimx => acc_values(2::6)
    brey => acc_values(3::6)
    bimy => acc_values(4::6)
    brez => acc_values(5::6)
    bimz => acc_values(6::6)
    allocate(visited(currentmodel % numberofnodes))
    n = currentmodel % mesh % maxelementdofs
    allocate(nodes % x(n), nodes % y(n), nodes % z(n))
    allocate(bacre(3,n), bacim(3,n))
    allocate(permeab(n), basis(n), dbasisdx(n,3))
    visited = 0
    lrfx = 0._dp
    lrfy = 0._dp
    lrfz = 0._dp
    do t = 1, currentmodel % numberofbulkelements
      el => currentmodel % elements(t)
      n = el % type % numberofnodes
      nidx => el % nodeindexes
      nodes % x(1:n) = currentmodel % nodes % x(nidx)
      nodes % y(1:n) = currentmodel % nodes % y(nidx)
      nodes % z(1:n) = currentmodel % nodes % z(nidx)
      permeab(1:n) = listgetreal(acc_material, &
          'Magnetic Permeability', n, nidx)
      k = listgetinteger(acc_model % bodies(el % bodyid) % values, &
          'Body Force', gotit, 1, acc_model % numberofbodyforces)
      bacre = 0._dp
      bacim = 0._dp
      if (k > 0) then
        bacre(1,1:n) = listgetreal(acc_model % bodyforces(k) % values, &
            'AC Magnetic Field 1 Im', n, nidx, gotit)
        bacim(1,1:n) = -listgetreal(acc_model % bodyforces(k) % values, &
            'AC Magnetic Field 1 Re', n, nidx, gotit)
        bacre(2,1:n) = listgetreal(acc_model % bodyforces(k) % values, &
            'AC Magnetic Field 2 Im', n, nidx, gotit)
        bacim(2,1:n) = -listgetreal(acc_model % bodyforces(k) % values, &
            'AC Magnetic Field 2 Re', n, nidx, gotit)
        bacre(3,1:n) = listgetreal(acc_model % bodyforces(k) % values, &
            'AC Magnetic Field 3 Im', n, nidx, gotit)
        bacim(3,1:n) = -listgetreal(acc_model % bodyforces(k) % values, &
            'AC Magnetic Field 3 Re', n, nidx, gotit)
      end if
      if (minval(acc_perm(nidx)) > 0) then
        do p = 1, n
          b = 0._dp
          dhdx = 0._dp
          q = acc_perm(nidx(p))
          u = el % type % nodeu(p)
          v = el % type % nodev(p)
          if (el % type % dimension == 3) then
            w = el % type % nodew(p)
          else
            w = 0.0d0
          end if
          stat = elementinfo(el, nodes, u, v, w, sqrtelem, &
              basis, dbasisdx)
          x = sum(nodes % x(1:n) * basis(1:n))
          y = sum(nodes % y(1:n) * basis(1:n))
          z = sum(nodes % z(1:n) * basis(1:n))
          call coordinatesysteminfo(metric, sqrtmetric, symb, dsymb, x, y, z)
          if (currentcoordinatesystem() /= cartesian) &
              call invertmatrix(metric, 3)
          l = 0._dp
          b(1) = sum(basis(1:n)*brex(acc_perm(nidx))) + bacre(1,p)
          b(2) = sum(basis(1:n)*brey(acc_perm(nidx))) + bacre(2,p)
          b(3) = sum(basis(1:n)*brez(acc_perm(nidx))) + bacre(3,p)
          do i = 1, 3
            dhdx(1,i) = sum(dbasisdx(1:n,i)* &
                brex(acc_perm(nidx))) / permeab(p)
            dhdx(2,i) = sum(dbasisdx(1:n,i)* &
                brey(acc_perm(nidx))) / permeab(p)
            dhdx(3,i) = sum(dbasisdx(1:n,i)* &
                brez(acc_perm(nidx))) / permeab(p)
          end do
          l = fd3_computelorentz(b, dhdx, permeab(p), &
              sqrtmetric, metric, symb)
          b(1) = sum(basis(1:n)*bimx(acc_perm(nidx))) + bacim(1,p)
          b(2) = sum(basis(1:n)*bimy(acc_perm(nidx))) + bacim(2,p)
          b(3) = sum(basis(1:n)*bimz(acc_perm(nidx))) + bacim(3,p)
          do i = 1, 3
            dhdx(1,i) = sum(dbasisdx(1:n,i)* &
                bimx(acc_perm(nidx))) / permeab(p)
            dhdx(2,i) = sum(dbasisdx(1:n,i)* &
                bimy(acc_perm(nidx))) / permeab(p)
            dhdx(3,i) = sum(dbasisdx(1:n,i)* &
                bimz(acc_perm(nidx))) / permeab(p)
          end do
          l = l + fd3_computelorentz(b, dhdx, permeab(p), &
              sqrtmetric, metric, symb)
          lrfx(q) = lrfx(q) + 0.5_dp * l(1)
          lrfy(q) = lrfy(q) + 0.5_dp * l(2)
          if (currentcoordinatesystem() == cylindricsymmetric) then
            l(3) = 0._dp
          end if
          lrfz(q) = lrfz(q) + 0.5_dp * l(3)
          visited(q) = visited(q) + 1
        end do
      end if
    end do
    do i = 1, currentmodel % numberofnodes
      if (visited(i) > 1) then
        lrfx(i) = lrfx(i) / visited(i)
        lrfy(i) = lrfy(i) / visited(i)
        lrfz(i) = lrfz(i) / visited(i)
      end if
    end do
    deallocate(visited)
    deallocate(nodes % x, nodes % y, nodes % z)
    deallocate(permeab)
    deallocate(bacre, bacim, basis, dbasisdx)
  end subroutine fd3_lorentz

  ! -- ComputeLorentz internal function (exact transcription) -----------
  function fd3_computelorentz(b, dhdx, mu, sqrtmetric, metric, symb) &
      result(lf)
    real(c_double) :: b(:), dhdx(:,:), mu, lf(3), sqrtmetric, &
        metric(:,:), symb(:,:,:)
    integer :: i, j, k, l, m
    real(c_double) :: bc(3), ji(3), jc(3), s, perm(3,3,3), r
    if (currentcoordinatesystem() == cartesian) then
      ji(1) = dhdx(3,2) - dhdx(2,3)
      ji(2) = dhdx(1,3) - dhdx(3,1)
      ji(3) = dhdx(2,1) - dhdx(1,2)
      lf(1) = ji(2)*b(3) - ji(3)*b(2)
      lf(2) = ji(3)*b(1) - ji(1)*b(3)
      lf(3) = ji(1)*b(2) - ji(2)*b(1)
      return
    end if
    r = sqrtmetric
    if (currentcoordinatesystem() == cylindricsymmetric) then
      ji(1) = -dhdx(3,2)
      ji(2) = dhdx(3,1)
      if (r > 1.0d-10) then
        ji(2) = ji(2) + b(3)/(r*mu)
      else
        ji(2) = ji(2) + ji(2)
      end if
      ji(3) = dhdx(1,2) - dhdx(2,1)
      lf(1) = ji(3)*b(2) - ji(2)*b(3)
      lf(2) = ji(1)*b(3) - ji(3)*b(1)
      lf(3) = ji(2)*b(1) - ji(1)*b(2)
      return
    end if
    perm = 0
    perm(1,2,3) = -1.0d0 / sqrtmetric
    perm(1,3,2) = 1.0d0 / sqrtmetric
    perm(2,1,3) = 1.0d0 / sqrtmetric
    perm(2,3,1) = -1.0d0 / sqrtmetric
    perm(3,1,2) = -1.0d0 / sqrtmetric
    perm(3,2,1) = 1.0d0 / sqrtmetric
    bc = 0.0d0
    do i = 1, 3
      do j = 1, 3
        bc(i) = bc(i) + metric(i,j)*b(j)
      end do
    end do
    ji = 0.0d0
    do i = 1, 3
      s = 0.0d0
      do j = 1, 3
        do k = 1, 3
          if (perm(i,j,k) /= 0) then
            do l = 1, 3
              s = s + perm(i,j,k)*metric(j,l)*dhdx(l,k)
              do m = 1, 3
                s = s + perm(i,j,k)*metric(j,l)*symb(k,m,l)*b(m)/mu
              end do
            end do
          end if
        end do
      end do
      ji(i) = s
    end do
    jc = 0.0d0
    do i = 1, 3
      do j = 1, 3
        jc(i) = jc(i) + metric(i,j)*ji(j)
      end do
    end do
    lf = 0.0d0
    do i = 1, 3
      s = 0.0d0
      do j = 1, 3
        do k = 1, 3
          if (perm(i,j,k) /= 0) then
            s = s + perm(i,j,k)*jc(k)*bc(j)
          end if
        end do
      end do
      lf(i) = s
    end do
  end function fd3_computelorentz

end module fdiffusion3d_acc
