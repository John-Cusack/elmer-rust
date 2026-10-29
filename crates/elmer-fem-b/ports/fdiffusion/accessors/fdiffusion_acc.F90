! T022-A FDiffusion port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/FDiffusion.F90, compiled against the built .mod
! files of the target configure (A2 shipped mechanism, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars). All Fortran-object traffic — TYPE(...)
! handles, OPTIONAL/POINTER dummies, hidden lengths, assumed-shape
! descriptor dummies, list-directed timing reports — stays inside Fortran,
! so the Rust side is pure C calls. The wrapper module mirrors the
! solver's own SAVE locals (entry SAVE arrays + flags, internal-procedure
! automatic arrays as per-call locals) with identical lifetimes. All
! floating-point work (complex Galerkin assembly, SUM reductions, the
! CPUTime/RealTime-based reports) stays in Fortran (same compiler), so
! bit-identity holds by construction wherever the call sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module fdiffusion_acc
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

  ! Entry SAVE state (mirrors FourierDiffusionSolver's SAVE statement).
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
  function fdi_begin(h_model, h_solver) result(code) bind(c, name="fdi_begin")
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
  end function fdi_begin

  subroutine fdi_alloc() bind(c, name="fdi_alloc")
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
          acc_force(2*n), acc_imp(2,n), acc_work(n), &
          acc_lstiff(2*n,2*n), acc_cond(n), acc_pm(n), &
          acc_lrf(3*acc_model % numberofnodes), acc_load(2,n), stat=istat)
      if (istat /= 0) then
        call fatal('FourierDiffusionSolve', 'Memory allocation error.')
      end if
      var => variableget(acc_model % variables, 'Lorentz Force')
      if (.not. associated(var)) then
        sptr => acc_solver
        call variableaddvector(acc_solver % mesh % variables, &
            acc_solver % mesh, sptr, 'Lorentz Force', 3, acc_lrf, acc_perm)
      end if
      acc_allocdone = .true.
    end if
  end subroutine fdi_alloc

  function fdi_params() result(nl) bind(c, name="fdi_params")
    integer(c_int) :: nl
    logical :: gotit
    nl = listgetinteger(acc_solver % values, &
        'Nonlinear System Max Iterations', gotit)
    if (.not. gotit) nl = 1
    acc_eqname = listgetstring(acc_solver % values, 'Equation')
    acc_bubbles = listgetlogical(acc_solver % values, 'Bubbles', gotit)
    if (.not. gotit) acc_bubbles = .false.
    acc_omega = getangularfrequency()
  end function fdi_params

  function fdi_nbulk() result(n) bind(c, name="fdi_nbulk")
    integer(c_int) :: n
    n = acc_solver % mesh % numberofbulkelements
  end function fdi_nbulk

  function fdi_nbnd() result(n) bind(c, name="fdi_nbnd")
    integer(c_int) :: n
    n = acc_solver % mesh % numberofboundaryelements
  end function fdi_nbnd

  ! -- nonlinear iteration head (Info banner + DefaultInitialize) -------
  subroutine fdi_iter_head(iter) bind(c, name="fdi_iter_head")
    integer(c_int), value :: iter
    acc_at = cputime()
    acc_at0 = realtime()
    call info('FourierDiffusionSolve', ' ', level=4)
    call info('FourierDiffusionSolve', ' ', level=4)
    call info('FourierDiffusionSolve', &
        '-------------------------------------', level=4)
    write(message, *) 'Fourier Diffusion  iteration', iter
    call info('FourierDiffusionSolve', message, level=4)
    write(message, *) 'Frequency (Hz): ', acc_omega / (2*pi)
    call info('FourierDiffusionSolve', message, level=4)
    call info('FourierDiffusionSolve', &
        '-------------------------------------', level=4)
    call info('FourierDiffusionSolve', ' ', level=4)
    call info('FourierDiffusionSolve', 'Starting Assembly...', level=4)
    call defaultinitialize()
  end subroutine fdi_iter_head

  ! -- bulk element body (exact transcription of the t-loop body) -------
  subroutine fdi_bulk(t) bind(c, name="fdi_bulk")
    integer(c_int), value :: t
    type(element_t), pointer :: cur
    integer, pointer :: nidx(:)
    integer :: n, k
    logical :: gotit
    cur => acc_solver % mesh % elements(t)
    if (realtime() - acc_at0 > 1.0d0) then
      write(message,'(a,i3,a)') '   Assembly: ', int(100.0d0 - 100.0d0 * &
          (acc_solver % mesh % numberofbulkelements - t) / &
          (1.0d0 * acc_solver % mesh % numberofbulkelements)), ' % done'
      call info('FourierDiffusionSolve', message, level=5)
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
      call warn('Fdiffusion', &
          'Use electric conductivity instead of electrical')
    else
      acc_cond(1:n) = listgetreal(acc_material, &
          'Electric Conductivity', n, nidx)
    end if
    acc_pm(1:n) = listgetreal(acc_material, &
        'Magnetic Permeability', n, nidx)
    k = listgetinteger(acc_model % bodies(cur % bodyid) % values, &
        'Body Force', gotit, 1, acc_model % numberofbodyforces)
    acc_load = 0.0d0
    if (k > 0) then
      acc_load(1,1:n) = listgetreal(acc_model % bodyforces(k) % values, &
          'AC Magnetic Field Re', n, nidx, gotit)
      acc_load(2,1:n) = listgetreal(acc_model % bodyforces(k) % values, &
          'AC Magnetic Field Im', n, nidx, gotit)
    end if
    call fdi_localmatrix(n)
    call updateglobalequations(acc_mat, acc_lstiff, acc_rhs, acc_force, &
        n, acc_solver % variable % dofs, acc_perm(nidx))
  end subroutine fdi_bulk

  ! -- boundary element body (exact transcription of the t-loop body) ---
  subroutine fdi_bnd(t) bind(c, name="fdi_bnd")
    integer(c_int), value :: t
    type(element_t), pointer :: cur
    integer, pointer :: nidx(:)
    integer :: i, n
    logical :: gotit
    cur => acc_solver % mesh % elements(t)
    acc_model % currentelement => cur
    acc_elem => cur
    do i = 1, acc_model % numberofbcs
      if (cur % boundaryinfo % constraint == acc_model % bcs(i) % tag) then
        n = cur % type % numberofnodes
        nidx => cur % nodeindexes
        if (any(acc_perm(nidx) == 0)) cycle
        acc_nodes % x(1:n) = acc_solver % mesh % nodes % x(nidx)
        acc_nodes % y(1:n) = acc_solver % mesh % nodes % y(nidx)
        acc_nodes % z(1:n) = acc_solver % mesh % nodes % z(nidx)
        acc_imp(1,1:n) = listgetreal(acc_model % bcs(i) % values, &
            'Wave Impedance 1', n, nidx, gotit)
        acc_imp(2,1:n) = listgetreal(acc_model % bcs(i) % values, &
            'Wave Impedance 2', n, nidx, gotit)
        acc_load(1,1:n) = listgetreal(acc_model % bcs(i) % values, &
            'Wave Flux 1', n, nidx, gotit)
        acc_load(2,1:n) = listgetreal(acc_model % bcs(i) % values, &
            'Wave Flux 2', n, nidx, gotit)
        call fdi_localmatrixboundary(n)
        call updateglobalequations(acc_mat, acc_lstiff, acc_rhs, &
            acc_force, n, acc_solver % variable % dofs, acc_perm(nidx))
      end if
    end do
  end subroutine fdi_bnd

  ! -- nonlinear iteration tail (finish, solve, timing reports) ---------
  function fdi_iter_tail(iter) result(conv) bind(c, name="fdi_iter_tail")
    integer(c_int), value :: iter
    integer(c_int) :: conv
    real(c_double) :: st, norm, relchange
    call defaultfinishassembly()
    call defaultdirichletbcs()
    call info('FourierDiffusionSolve', 'Assembly done', level=4)
    acc_at = cputime() - acc_at
    st = cputime()
    norm = defaultsolve()
    relchange = acc_solver % variable % nonlinchange
    st = cputime() - st
    acc_totat = acc_totat + acc_at
    acc_totst = acc_totst + st
    write(message,'(a,i4,a,F8.2,F8.2)') 'iter: ', iter, &
        ' Assembly: (s)', acc_at, acc_totat
    call info('FourierDiffusionSolve', message, level=4)
    write(message,'(a,i4,a,F8.2,F8.2)') 'iter: ', iter, &
        ' Solve:    (s)', st, acc_totst
    call info('FourierDiffusionSolve', message, level=4)
    write(message, *) 'Result Norm    : ', norm
    call info('FourierDiffusionSolve', message, level=4)
    write(message, *) 'Relative Change: ', relchange
    call info('FourierDiffusionSolve', message, level=4)
    if (acc_solver % variable % nonlinconverged == 1) then
      conv = 1
    else
      conv = 0
    end if
  end function fdi_iter_tail

  ! -- LocalMatrix internal procedure (exact transcription) --------------
  subroutine fdi_localmatrix(n)
    integer, intent(in) :: n
    real(c_double) :: basis(2*n), dbasisdx(2*n,3)
    real(c_double) :: sqrtelem, u, v, w, s, cond, d, m, l1, l2
    real(c_double) :: x, y, z
    real(c_double) :: metric(3,3), sqrtmetric, symb(3,3,3), dsymb(3,3,3,3)
    real(c_double) :: nodcond(n)
    complex(kind=dp) :: lstiff(2*n,2*n), lforce(2*n), a
    logical :: stat
    integer :: i, j, p, q, tt, dim, nbasis, coordsys
    type(gaussintegrationpoints_t) :: ip
    dim = coordinatesystemdimension()
    coordsys = currentcoordinatesystem()
    metric = 0.0d0
    metric(1,1) = 1.0d0
    metric(2,2) = 1.0d0
    metric(3,3) = 1.0d0
    lstiff = 0.0d0
    lforce = 0.0d0
    if (acc_bubbles) then
      ip = gausspoints(acc_elem, acc_elem % type % gausspoints2)
      nbasis = 2*n
    else
      nbasis = n
      ip = gausspoints(acc_elem)
    end if
    nodcond(1:n) = acc_cond(1:n) * acc_pm(1:n)
    do tt = 1, ip % n
      u = ip % u(tt)
      v = ip % v(tt)
      w = ip % w(tt)
      s = ip % s(tt)
      stat = elementinfo(acc_elem, acc_nodes, u, v, w, sqrtelem, &
          basis, dbasisdx, bubbles=acc_bubbles)
      s = s * sqrtelem
      if (coordsys /= cartesian) then
        x = sum(acc_nodes % x(1:n) * basis(1:n))
        y = sum(acc_nodes % y(1:n) * basis(1:n))
        z = sum(acc_nodes % z(1:n) * basis(1:n))
        call coordinatesysteminfo(metric, sqrtmetric, symb, dsymb, x, y, z)
        s = s * sqrtmetric
      end if
      cond = sum(nodcond(1:n) * basis)
      d = -acc_omega
      m = 0._dp
      l1 = acc_omega * sum(acc_load(1,1:n) * basis(1:n))
      l2 = acc_omega * sum(acc_load(2,1:n) * basis(1:n))
      do p = 1, nbasis
        do q = 1, nbasis
          a = cmplx(m, d, kind=dp) * basis(q) * basis(p)
          do i = 1, dim
            do j = 1, dim
              a = a + metric(i,j) * dbasisdx(q,i) * dbasisdx(p,j) / cond
            end do
          end do
          lstiff(p,q) = lstiff(p,q) + s*a
        end do
        lforce(p) = lforce(p) + s * basis(p) * cmplx(l1,l2,kind=dp)
      end do
    end do
    if (acc_bubbles) then
      call condensatep(n, n, lstiff, lforce)
    end if
    do i = 1, n
      acc_force(2*(i-1)+1) = real(lforce(i))
      acc_force(2*(i-1)+2) = aimag(lforce(i))
      do j = 1, n
        acc_lstiff(2*(i-1)+1, 2*(j-1)+1) = real(lstiff(i,j))
        acc_lstiff(2*(i-1)+1, 2*(j-1)+2) = -aimag(lstiff(i,j))
        acc_lstiff(2*(i-1)+2, 2*(j-1)+1) = aimag(lstiff(i,j))
        acc_lstiff(2*(i-1)+2, 2*(j-1)+2) = real(lstiff(i,j))
      end do
    end do
  end subroutine fdi_localmatrix

  ! -- LocalMatrixBoundary internal procedure (exact transcription) -----
  subroutine fdi_localmatrixboundary(n)
    integer, intent(in) :: n
    real(c_double) :: sqrtelem, u, v, w, s, imp1, imp2, l1, l2
    real(c_double) :: basis(n), dbasisdx(n,3), x, y, z
    complex(kind=dp) :: lstiff(n,n), lforce(n), a
    logical :: stat
    integer :: i, j, p, q, tt, dim, coordsys
    type(gaussintegrationpoints_t) :: ip
    dim = coordinatesystemdimension()
    coordsys = currentcoordinatesystem()
    lstiff = 0.0d0
    lforce = 0.0d0
    ip = gausspoints(acc_elem)
    do tt = 1, ip % n
      u = ip % u(tt)
      v = ip % v(tt)
      w = ip % w(tt)
      s = ip % s(tt)
      stat = elementinfo(acc_elem, acc_nodes, u, v, w, sqrtelem, &
          basis, dbasisdx)
      s = s * sqrtelem
      if (coordsys /= cartesian) then
        x = sum(acc_nodes % x(1:n) * basis(1:n))
        y = sum(acc_nodes % y(1:n) * basis(1:n))
        z = sum(acc_nodes % z(1:n) * basis(1:n))
        s = s * coordinatesqrtmetric(x, y, z)
      end if
      imp1 = sum(acc_imp(1,1:n) * basis)
      if (abs(imp1) > aeps) imp1 = acc_omega / imp1
      imp2 = sum(acc_imp(2,1:n) * basis)
      if (abs(imp2) > aeps) imp2 = acc_omega / imp2
      l1 = sum(acc_load(1,1:n) * basis)
      l2 = sum(acc_load(2,1:n) * basis)
      do p = 1, n
        do q = 1, n
          a = cmplx(imp1, imp2, kind=dp) * basis(q) * basis(p)
          lstiff(p,q) = lstiff(p,q) + s * a
        end do
        lforce(p) = lforce(p) + s * basis(p) * cmplx(l1,l2,kind=dp)
      end do
    end do
    do i = 1, n
      acc_force(2*(i-1)+1) = real(lforce(i))
      acc_force(2*(i-1)+2) = aimag(lforce(i))
      do j = 1, n
        acc_lstiff(2*(i-1)+1, 2*(j-1)+1) = real(lstiff(i,j))
        acc_lstiff(2*(i-1)+1, 2*(j-1)+2) = -aimag(lstiff(i,j))
        acc_lstiff(2*(i-1)+2, 2*(j-1)+1) = aimag(lstiff(i,j))
        acc_lstiff(2*(i-1)+2, 2*(j-1)+2) = real(lstiff(i,j))
      end do
    end do
  end subroutine fdi_localmatrixboundary

  ! -- LorentzForceAve internal procedure (exact transcription) ---------
  subroutine fdi_lorentz() bind(c, name="fdi_lorentz")
    real(c_double), pointer :: lrfr(:), lrfz(:), lrfp(:)
    real(c_double), pointer :: bre(:), bim(:)
    type(element_t), pointer :: el
    type(nodes_t) :: nodes
    logical :: stat, gotit
    integer, pointer :: nidx(:), visited(:)
    integer :: p, q, i, tt, n, k
    real(c_double) :: u, v, w, sqrtelem
    real(c_double), allocatable :: basis(:), dbasisdx(:,:)
    real(c_double), allocatable :: permeab(:), bacre(:), bacim(:)
    lrfr => acc_lrf(1::3)
    lrfz => acc_lrf(2::3)
    lrfp => acc_lrf(3::3)
    bre => acc_values(1::2)
    bim => acc_values(2::2)
    allocate(visited(currentmodel % numberofnodes))
    n = currentmodel % mesh % maxelementdofs
    allocate(permeab(n), basis(n), dbasisdx(n,3))
    allocate(bacre(n), bacim(n))
    allocate(nodes % x(n), nodes % y(n), nodes % z(n))
    visited = 0
    lrfr = 0.0d0
    lrfz = 0.0d0
    lrfp = 0.0d0
    do tt = 1, currentmodel % numberofbulkelements
      el => currentmodel % elements(tt)
      n = el % type % numberofnodes
      nidx => el % nodeindexes
      nodes % x(1:n) = currentmodel % nodes % x(nidx)
      nodes % y(1:n) = currentmodel % nodes % y(nidx)
      nodes % z(1:n) = currentmodel % nodes % z(nidx)
      permeab(1:n) = listgetreal(acc_material, &
          'Magnetic Permeability', n, nidx)
      k = listgetinteger(acc_model % bodies(el % bodyid) % values, &
          'Body Force', gotit, 1, acc_model % numberofbodyforces)
      bacre = 0.0d0
      bacim = 0.0d0
      if (k > 0) then
        bacre(1:n) = listgetreal(acc_model % bodyforces(k) % values, &
            'AC Magnetic Field Re', n, nidx, gotit)
        bacim(1:n) = listgetreal(acc_model % bodyforces(k) % values, &
            'AC Magnetic Field Im', n, nidx, gotit)
      end if
      if (minval(acc_perm(nidx)) > 0) then
        do p = 1, n
          q = acc_perm(nidx(p))
          u = el % type % nodeu(p)
          v = el % type % nodev(p)
          if (el % type % dimension == 3) then
            w = el % type % nodew(p)
          else
            w = 0.0d0
          end if
          stat = elementinfo(el, nodes, u, v, w, sqrtelem, basis, dbasisdx)
          lrfr(q) = lrfr(q) - 0.5_dp * ( &
              sum(dbasisdx(1:n,1) * bre(acc_perm(nidx))) * &
              sum(basis(1:n) * bre(acc_perm(nidx))) + &
              sum(dbasisdx(1:n,1) * bim(acc_perm(nidx))) * &
              sum(basis(1:n) * bim(acc_perm(nidx)))) / permeab(p)
          lrfr(q) = lrfr(q) - 0.5_dp * ( &
              sum(dbasisdx(1:n,1) * bre(acc_perm(nidx))) * &
              bacim(p)) / permeab(p)
          visited(q) = visited(q) + 1
        end do
      end if
    end do
    do i = 1, currentmodel % numberofnodes
      if (visited(i) > 1) then
        lrfr(i) = lrfr(i) / visited(i)
        lrfp(i) = lrfp(i) / visited(i)
        lrfz(i) = lrfz(i) / visited(i)
      end if
    end do
    deallocate(visited)
    deallocate(nodes % x, nodes % y, nodes % z)
    deallocate(permeab, basis, dbasisdx)
    deallocate(bacre, bacim)
  end subroutine fdi_lorentz

end module fdiffusion_acc
