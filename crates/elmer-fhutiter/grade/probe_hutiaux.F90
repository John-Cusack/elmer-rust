program probe_hutiaux
  use huti_aux
  implicit none
  integer, parameter :: nd = 4
  integer :: ipar(31)
  integer :: i, j, n
  real :: sx(4,4), su(4), sv(4), sw(4)
  double precision :: dx(4,4), du(4), dv(4), dw(4)
  complex :: cx(3,3), cu(3), cv(3)
  double complex :: zx(3,3), zu(3), zv(3)
  real :: sr(4)
  double precision :: dr(4)
  complex :: cr(4)
  double complex :: zr(4)
  logical :: ok

  ipar = 0
  ipar(3) = nd

  ! ---- pcondfun: pure copy ----
  do i = 1, nd
    sw(i) = 10.0 * i + 0.5
    dw(i) = 100.0 * i + 0.25
    cr(i) = cmplx(10.0 * i + 0.5, -1.0 * i)
    zr(i) = cmplx(100.0d0 * i + 0.25d0, -2.0d0 * i)
  end do
  call huti_sdummy_pcondfun(su, sw, ipar)
  call huti_ddummy_pcondfun(du, dw, ipar)
  ipar(3) = 3
  call huti_cdummy_pcondfun(cu, cr(1:3), ipar)
  call huti_zdummy_pcondfun(zu, zr(1:3), ipar)
  ipar(3) = nd
  write(*,'(A,4(E25.17,1X))') 'S_PC:', su
  write(*,'(A,4(E25.17,1X))') 'D_PC:', du
  write(*,'(A,3(E25.17,1X,E25.17,1X))') 'C_PC:', (cu(i), i=1,3)

  write(*,'(A,3(E25.17,1X,E25.17,1X))') 'Z_PC:', (zu(i), i=1,3)

  ! ---- lusolve s/d: fixed 4x4 ----
  sx = reshape([4.0, 1.0, 2.0, 0.5, 1.0, 5.0, 3.0, 1.0, &
                2.0, 3.0, 6.0, 2.0, 0.5, 1.0, 2.0, 7.0], [4,4])
  sv = [7.0, 8.0, 9.0, 10.0]
  su = 0.0
  n = 4
  call huti_slusolve(n, sx, su, sv)
  write(*,'(A,4(E25.17,1X))') 'S_U:', su
  write(*,'(A,16(E25.17,1X))') 'S_LU:', ((sx(i,j), j=1,4), i=1,4)

  dx = reshape([4.0d0, 1.0d0, 2.0d0, 0.5d0, 1.0d0, 5.0d0, 3.0d0, 1.0d0, &
                2.0d0, 3.0d0, 6.0d0, 2.0d0, 0.5d0, 1.0d0, 2.0d0, 7.0d0], [4,4])
  dv = [7.0d0, 8.0d0, 9.0d0, 10.0d0]
  du = 0.0d0
  call huti_dlusolve(n, dx, du, dv)
  write(*,'(A,4(E25.17,1X))') 'D_U:', du
  write(*,'(A,16(E25.17,1X))') 'D_LU:', ((dx(i,j), j=1,4), i=1,4)

  ! ---- lusolve c/z: fixed 3x3 ----
  cx = reshape([(4.0,1.0), (1.0,-2.0), (2.0,0.5), &
                (1.0,2.0), (5.0,-1.0), (3.0,1.5), &
                (2.0,-0.5), (3.0,-1.5), (6.0,2.0)], [3,3])
  cv = [(7.0,-1.0), (8.0,2.0), (9.0,0.0)]
  cu = (0.0, 0.0)
  n = 3
  call huti_clusolve(n, cx, cu, cv)
  write(*,'(A,3(E25.17,1X,E25.17,1X))') 'C_U:', (cu(i), i=1,3)
  write(*,'(A,9(E25.17,1X,E25.17,1X))') 'C_LU:', ((cx(i,j), j=1,3), i=1,3)

  zx = cx
  zv = cv
  zu = (0.0d0, 0.0d0)
  call huti_zlusolve(n, zx, zu, zv)
  write(*,'(A,3(E25.17,1X,E25.17,1X))') 'Z_U:', (zu(i), i=1,3)
  write(*,'(A,9(E25.17,1X,E25.17,1X))') 'Z_LU:', ((zx(i,j), j=1,3), i=1,3)

  ! ---- pivot diagnostic: tiny positive + tiny negative pivots ----
  sx = 0.0
  sx(1,1) = 1.0e-17
  sx(2,2) = 1.0
  sx(3,3) = 1.0
  sx(4,4) = 1.0
  sv = [1.0, 2.0, 3.0, 4.0]
  su = 0.0
  n = 4
  call huti_slusolve(n, sx, su, sv)
  write(*,'(A,4(E25.17,1X))') 'S_PIVU:', su

  sx = 0.0
  sx(1,1) = -1.0e-17
  sx(2,2) = 1.0
  sx(3,3) = 1.0
  sx(4,4) = 1.0
  call huti_slusolve(n, sx, su, sv)
  write(*,'(A,4(E25.17,1X))') 'S_NPIVU:', su

  dx = 0.0d0
  dx(1,1) = 1.0d-17
  dx(2,2) = 1.0d0
  dx(3,3) = 1.0d0
  dx(4,4) = 1.0d0
  dv = [1.0d0, 2.0d0, 3.0d0, 4.0d0]
  du = 0.0d0
  call huti_dlusolve(n, dx, du, dv)
  write(*,'(A,4(E25.17,1X))') 'D_PIVU:', du

  cx = (0.0, 0.0)
  cx(1,1) = (1.0e-17, 0.0)
  cx(2,2) = (1.0, 0.0)
  cx(3,3) = (1.0, 0.0)
  cv = [(1.0,0.0), (2.0,0.0), (3.0,0.0)]
  cu = (0.0, 0.0)
  n = 3
  call huti_clusolve(n, cx, cu, cv)
  write(*,'(A,3(E25.17,1X,E25.17,1X))') 'C_PIVU:', (cu(i), i=1,3)

  zx = (0.0d0, 0.0d0)
  zx(1,1) = (1.0d-17, 0.0d0)
  zx(2,2) = (1.0d0, 0.0d0)
  zx(3,3) = (1.0d0, 0.0d0)
  zv = [(1.0d0,0.0d0), (2.0d0,0.0d0), (3.0d0,0.0d0)]
  zu = (0.0d0, 0.0d0)
  call huti_zlusolve(n, zx, zu, zv)
  write(*,'(A,3(E25.17,1X,E25.17,1X))') 'Z_PIVU:', (zu(i), i=1,3)

  ! ---- edges: n = 1 and n = 0 (s only; loop-shape identical all kinds) ----
  n = 1
  sx(1,1) = 2.0
  sv(1) = 6.0
  su(1) = 0.0
  call huti_slusolve(n, sx, su, sv)
  write(*,'(A,1(E25.17,1X))') 'S_N1U:', su(1)
  n = 0
  call huti_slusolve(n, sx, su, sv)
  write(*,'(A)') 'S_N0_OK'

  ! ---- randvec: shape/range/state checks (never the stream) ----
  ipar(3) = 4
  call huti_srandvec(sr, ipar)
  ok = .true.
  do i = 1, 4
    if (sr(i) < 0.0 .or. sr(i) >= 1.0) ok = .false.
  end do
  call huti_srandvec(sw, ipar)
  if (all(sw == sr)) ok = .false.
  if (ok) write(*,'(A)') 'S_RAND_OK'

  call huti_drandvec(dr, ipar)
  ok = .true.
  do i = 1, 4
    if (dr(i) < 0.0d0 .or. dr(i) >= 1.0d0) ok = .false.
  end do
  call huti_drandvec(dw, ipar)
  if (all(dw == dr)) ok = .false.
  if (ok) write(*,'(A)') 'D_RAND_OK'

  call huti_crandvec(cr, ipar)
  ok = .true.
  do i = 1, 4
    if (real(cr(i)) /= 0.0 .or. aimag(cr(i)) < 0.0 .or. aimag(cr(i)) >= 1.0) ok = .false.
  end do
  if (ok) write(*,'(A)') 'C_RAND_OK'

  call huti_zrandvec(zr, ipar)
  ok = .true.
  do i = 1, 4
    if (real(zr(i)) /= 0.0d0 .or. aimag(zr(i)) < 0.0d0 .or. aimag(zr(i)) >= 1.0d0) ok = .false.
  end do
  if (ok) write(*,'(A)') 'Z_RAND_OK'
end program probe_hutiaux
