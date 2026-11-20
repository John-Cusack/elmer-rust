! O7 kernel shim (report section 5.2 item 8): C entry points over the
! pristine Elmer kernels, so the oracle's Rust tests can call them.
! Oracle-only. It is compiled against one pristine arm's module files
! and linked to that arm's libelmersolver; no candidate crate links it.
!
! SPDX-License-Identifier: LGPL-2.1-or-later
MODULE o7_shim
  USE ISO_C_BINDING
  USE DefUtils
  IMPLICIT NONE
  TYPE(Element_t), POINTER, SAVE, PRIVATE :: el => NULL()
CONTAINS

  ! Element descriptions once (ELMER_HOME locates elements.def), and the
  ! model dimension ElementMetric reads through CoordinateSystemDimension.
  SUBROUTINE o7_init(cdim) BIND(C, NAME='o7_init')
    INTEGER(C_INT), VALUE :: cdim
    LOGICAL, SAVE :: done = .FALSE.
    IF (.NOT. done) THEN
      IF (.NOT. ASSOCIATED(CurrentModel)) ALLOCATE(CurrentModel)
      CALL InitializeElementDescriptions()
      el => AllocateElement()
      done = .TRUE.
    END IF
    CurrentModel % Dimension = cdim
  END SUBROUTINE o7_init

  ! One Integration.F90 rule builder. kind: 1 line, 2 tri, 3 quad, 4 tet,
  ! 5 pyr, 6 wedge, 7 brick, 8 wedge2(a,b), 9 wedge economic. Returns the
  ! point count, -1 for an unknown kind, -2 when nmax is too small.
  FUNCTION o7_rule(kind, a, b, nmax, u, v, w, s) RESULT(n) BIND(C, NAME='o7_rule')
    INTEGER(C_INT), VALUE :: kind, a, b, nmax
    REAL(C_DOUBLE) :: u(nmax), v(nmax), w(nmax), s(nmax)
    INTEGER(C_INT) :: n
    TYPE(GaussIntegrationPoints_t) :: IP
    INTEGER :: ia, ib
    ia = a
    ib = b
    SELECT CASE (kind)
    CASE (1); IP = GaussPoints1D(ia)
    CASE (2); IP = GaussPointsTriangle(ia)
    CASE (3); IP = GaussPointsQuad(ia)
    CASE (4); IP = GaussPointsTetra(ia)
    CASE (5); IP = GaussPointsPyramid(ia)
    CASE (6); IP = GaussPointsWedge(ia)
    CASE (7); IP = GaussPointsBrick(ia)
    CASE (8); IP = GaussPointsWedge2(ia, ib)
    CASE (9); IP = GaussPointsWedgeEconomic(ia)
    CASE DEFAULT
      n = -1
      RETURN
    END SELECT
    n = IP % n
    IF (n > nmax) THEN
      n = -2
      RETURN
    END IF
    u(1:n) = IP % u(1:n)
    v(1:n) = IP % v(1:n)
    w(1:n) = IP % w(1:n)
    s(1:n) = IP % s(1:n)
  END FUNCTION o7_rule

  ! ElementInfo for element type `code` with nodes (x,y,z) at (u,v,w) in a
  ! model of dimension cdim: detJ, Basis(nn), dBasisdx(nn,3) column-major.
  ! Returns 1 when ElementInfo succeeds, 0 when it fails, -1 when `code`
  ! is not a known type or nn is not its node count.
  FUNCTION o7_info(code, cdim, nn, x, y, z, u, v, w, detj, basis, dbasis) &
      RESULT(stat) BIND(C, NAME='o7_info')
    INTEGER(C_INT), VALUE :: code, cdim, nn
    REAL(C_DOUBLE) :: x(nn), y(nn), z(nn)
    REAL(C_DOUBLE), VALUE :: u, v, w
    REAL(C_DOUBLE) :: detj, basis(nn), dbasis(nn, 3)
    INTEGER(C_INT) :: stat
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp) :: uu, vv, ww, dj
    CALL o7_init(cdim)
    el % TYPE => GetElementType(INT(code), .FALSE.)
    IF (.NOT. ASSOCIATED(el % TYPE)) THEN
      stat = -1
      RETURN
    END IF
    IF (el % TYPE % NumberOfNodes /= nn) THEN
      stat = -1
      RETURN
    END IF
    ALLOCATE(Nodes % x(nn), Nodes % y(nn), Nodes % z(nn))
    Nodes % NumberOfNodes = nn
    Nodes % x = x
    Nodes % y = y
    Nodes % z = z
    uu = u
    vv = v
    ww = w
    basis = 0
    dbasis = 0
    IF (ElementInfo(el, Nodes, uu, vv, ww, dj, basis, dbasis)) THEN
      stat = 1
    ELSE
      stat = 0
    END IF
    detj = dj
    DEALLOCATE(Nodes % x, Nodes % y, Nodes % z)
  END FUNCTION o7_info

END MODULE o7_shim
