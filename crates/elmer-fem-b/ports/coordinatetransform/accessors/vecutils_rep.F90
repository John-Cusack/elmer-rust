! T022-A CoordinateTransform port: statement-for-statement transcription of
! the unit's MODULE VecUtils from fem/src/modules/CoordinateTransform.F90.
!
! Triage outcome (plain ABI): every dummy is an explicit-shape REAL(8)
! array — `a(3)`, `OldCoord(3,3)`, `NewCoord(3,3)`, `A(3,3)`, `jac(3,3)` —
! so no callee-side descriptor ABI is involved and the procedures port
! without new gate shapes. The results are array-valued
! (`norm2(3)`, `normalized(3)`, `jac(3,3)`, `transform2rank` RESULT `B(3,3)`);
! gfortran passes those as plain caller-side temporaries, exactly as the
! reference object does, because this transcription is compiled by the same
! compiler with the same flags.
!
! The evaluation order below mirrors the original term for term
! (`SUM(a**2)` broadcast, the `DOT_PRODUCT/SQRT` quotient, the `i,j,k,l`
! accumulation), so bit-identity holds by construction. The accessor
! module (`coordinatetransform_acc.F90`) USEs this module for its internal
! calls, so the exported symbols and the solver body share one source.
!
! Clean-room code: transcribed from the unit's VecUtils, no upstream text.
module vecutils
  use defutils, only: dp
  implicit none

contains

  function norm2(a)
    implicit none
    real(kind=dp) :: norm2(3)
    real(kind=dp), intent(in) :: a(3)
    norm2 = sqrt(sum(a**2))
  end function norm2

  function normalized(a)
    implicit none
    real(kind=dp) :: normalized(3)
    real(kind=dp), intent(in) :: a(3)
    normalized = a/norm2(a)
  end function normalized

  ! jacobian of the coordinate transform between the old coordinate system
  ! and the new coordinate system
  function jac(oldcoord, newcoord)
    implicit none
    real(kind=dp) :: oldcoord(3,3), newcoord(3,3)
    real(kind=dp) :: jac(3,3)
    integer :: i, j
    do i = 1, 3
      do j = 1, 3
        jac(i,j) = dot_product(newcoord(i,1:3), oldcoord(j,1:3)) &
             / sqrt(dot_product(oldcoord(j,1:3), oldcoord(j,1:3)))
      end do
    end do
  end function jac

  ! transform a 3x3 2-rank tensor with the given jacobian
  function transform2rank(a, jac) result (b)
    implicit none
    real(kind=dp) :: a(3,3), b(3,3)
    real(kind=dp) :: jac(3,3)
    integer :: i, j, k, l
    b = 0.0_dp
    do i = 1, 3
      do j = 1, 3
        do k = 1, 3
          do l = 1, 3
            b(i,j) = b(i,j) + jac(i,l)*jac(j,k) * a(k,l)
          end do
        end do
      end do
    end do
  end function transform2rank

end module vecutils
