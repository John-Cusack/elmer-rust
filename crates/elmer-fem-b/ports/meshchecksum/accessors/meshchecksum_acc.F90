! T021-B MeshChecksum port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/MeshChecksum.F90, compiled against the
! built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls.
!
! The checksum accumulation, the pseudo-norm loop, and every PRINT stay
! in Fortran as exact transcriptions, so bit-identity holds by
! construction wherever the call sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module meshchecksum_acc
  use iso_c_binding
  use types
  use defutils
  use meshbasics
  use messages
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()

contains

  ! -- entry state ----------------------------------------------------
  subroutine mca_begin(h_solver) bind(c, name="mca_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine mca_begin

  ! -- MeshChecksum_init body ------------------------------------------
  subroutine mca_init_var() bind(c, name="mca_init_var")
    call ListAddNewString(acc_solver % values, 'Variable', &
        '-nooutput -global meshcheck_var')
  end subroutine mca_init_var

  ! -- MeshChecksum body ------------------------------------------------
  subroutine mca_info() bind(c, name="mca_info")
    call Info('MeshChecksum', 'Checking for mesh consistency')
  end subroutine mca_info

  subroutine mca_checksums(cs) bind(c, name="mca_checksums")
    real(c_double) :: cs(8)
    type(element_t), pointer :: element
    type(mesh_t), pointer :: mesh
    integer :: t, ind, tag
    ! Exact transcription of the accumulation loop.
    mesh => acc_solver % mesh
    cs = 0.0_c_double
    cs(1) = mesh % numberofbulkelements + mesh % numberofbulkelements + &
        mesh % numberofnodes
    do t = 1, mesh % numberofbulkelements + mesh % numberofboundaryelements
      element => mesh % elements(t)
      if (t <= mesh % numberofbulkelements) then
        ind = 2
        tag = element % bodyid
      else
        ind = 5
        tag = 0
        if (associated(element % boundaryinfo)) then
          tag = element % boundaryinfo % constraint
        end if
      end if
      cs(ind) = cs(ind) + element % type % elementcode
      cs(ind+1) = cs(ind+1) + tag
      cs(ind+2) = cs(ind+2) + sum(element % nodeindexes)
    end do
    cs(8) = sum(mesh % nodes % x) + &
        sum(mesh % nodes % y) + sum(mesh % nodes % z)
  end subroutine mca_checksums

  subroutine mca_print_this(cs) bind(c, name="mca_print_this")
    real(c_double) :: cs(8)
    print *, 'Checksums for file output:'
    print *, 'ThisResults:', nint(cs(1:7)), cs(8)
  end subroutine mca_print_this

  function mca_get_ref(refv) result(found) bind(c, name="mca_get_ref")
    real(c_double) :: refv(8)
    integer(c_int) :: found
    real(kind=dp), pointer :: wrk(:,:)
    logical :: f
    wrk => ListGetConstRealArray(acc_solver % values, 'Reference Values', f)
    found = merge(1, 0, f)
    if (f) refv(1:8) = wrk(1:8, 1)
  end function mca_get_ref

  subroutine mca_print_ref(refv) bind(c, name="mca_print_ref")
    real(c_double) :: refv(8)
    print *, 'RefResults:', nint(refv(1:7)), refv(8)
  end subroutine mca_print_ref

  function mca_pseudonorm(cs, refv) result(p) bind(c, name="mca_pseudonorm")
    real(c_double) :: cs(8), refv(8), p
    real(kind=dp) :: c
    integer :: i, nsize
    ! Exact transcription of the pseudo-norm loop.
    nsize = 8
    p = 0.0_c_double
    do i = 1, nsize
      if (abs(refv(i)) > epsilon(c)) then
        c = cs(i) / refv(i)
        c = max(c, 1.0_dp / c)
      else
        c = 1.0_dp + cs(i)
      end if
      p = p + c / nsize
    end do
  end function mca_pseudonorm

  subroutine mca_print_pseudo(p) bind(c, name="mca_print_pseudo")
    real(c_double), value :: p
    print *, 'PseudoNorm:', p
  end subroutine mca_print_pseudo

  subroutine mca_set_pseudo(p) bind(c, name="mca_set_pseudo")
    real(c_double), value :: p
    ! Exact transcription of the result store.
    acc_solver % variable % values = p
    call ListAddNewConstReal(acc_solver % values, 'Reference Norm', 1.0_dp)
  end subroutine mca_set_pseudo

end module meshchecksum_acc
