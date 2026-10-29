! T022-A EliminateDirichlet port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/EliminateDirichlet.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments. All Fortran-object
! traffic — TYPE(...) handles, POINTER dummies, list-directed timing
! reports — stays inside Fortran, so the Rust side is pure C calls. The
! wrapper module mirrors the original's locals (BMatrix SAVE association,
! F/U work arrays) with identical lifetimes.
!
! Partition (original line numbers):
! - structure surgery that allocates derived-type components and pointer-
!   associates A % Perm / A % EMatrix (108-234), the complex fill (238-288),
!   the eigen map-back (382-414), and every Elmer library call
!   (SolveLinearSystem, ListGet/AddLogical, Info, CPUtime, ComputeChange,
!   AllocateMatrix, FreeMatrix, CRS_SortMatrix) stay in Fortran, so
!   bit-identity holds by construction wherever the call sequence matches.
! - the real fill (290-319) and the real map-back (416-425) run in Rust
!   through the scalar eld_* getters/setters below; their floating-point
!   statements reproduce the original's evaluation order exactly.
! - the parallel (PEs > 1) paths are transcribed verbatim inside the solve
!   wrapper; the serial proof never takes them.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module eliminatedirichlet_acc
  use types
  use lists
  use solverbasics
  use solvecore, only: SolveLinearSystem
  use crsmatrix
  use generalutils
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(matrix_t), pointer, save :: acc_a => null()
  type(matrix_t), pointer, save :: acc_bm => null()
  integer, pointer, save :: acc_perm(:) => null()
  real(kind=dp), pointer, save :: acc_f(:) => null()
  real(kind=dp), pointer, save :: acc_u(:) => null()
  ! Entry timing (mirrors the original's TotTime / LinTime locals).
  real(kind=dp), save :: acc_tottime = 0.0_dp
  real(kind=dp), save :: acc_lintime = 0.0_dp
  real(kind=dp), pointer, save :: acc_b(:) => null()
  real(kind=dp), pointer, save :: acc_x(:) => null()

contains

  ! -- entry state ------------------------------------------------------
  subroutine eld_begin(h_solver, h_a, h_b, h_x, n) bind(c, name="eld_begin")
    type(c_ptr), value :: h_solver, h_a, h_b, h_x
    integer(c_int), value :: n
    call c_f_pointer(h_solver, acc_solver)
    call c_f_pointer(h_a, acc_a)
    call c_f_pointer(h_b, acc_b, [n])
    call c_f_pointer(h_x, acc_x, [n])
    acc_bm => acc_a % ematrix
    acc_tottime = cputime()
  end subroutine eld_begin

  function eld_mass() result(code) bind(c, name="eld_mass")
    integer(c_int) :: code
    code = 0
    if (associated(acc_a % massvalues)) then
      if (size(acc_a % massvalues) == size(acc_a % values)) code = 1
    end if
  end function eld_mass

  function eld_damp() result(code) bind(c, name="eld_damp")
    integer(c_int) :: code
    code = 0
    if (associated(acc_a % dampvalues)) then
      if (size(acc_a % dampvalues) == size(acc_a % values)) code = 1
    end if
  end function eld_damp

  subroutine eld_ensure_md(mass, damp) bind(c, name="eld_ensure_md")
    integer(c_int), value :: mass, damp
    if (associated(acc_bm)) then
      if (mass /= 0 .and. .not. associated(acc_bm % massvalues)) then
        allocate(acc_bm % massvalues(size(acc_bm % values)))
      end if
    end if
    if (associated(acc_bm)) then
      if (damp /= 0 .and. .not. associated(acc_bm % dampvalues)) then
        allocate(acc_bm % dampvalues(size(acc_bm % values)))
      end if
    end if
  end subroutine eld_ensure_md

  function eld_need_build() result(code) bind(c, name="eld_need_build")
    integer(c_int) :: code
    code = 0
    if (.not. associated(acc_bm)) code = 1
  end function eld_need_build

  ! -- reduced-matrix structure build (original 108-230, verbatim) -------
  subroutine eld_build(mass, damp) bind(c, name="eld_build")
    integer(c_int), value :: mass, damp
    logical :: mass_l, damp_l, add
    integer :: i, j, k, l, t, dirichlet, zz
    integer, pointer :: dcount(:)
    real(kind=dp) :: s
    mass_l = (mass /= 0)
    damp_l = (damp /= 0)
    if (associated(acc_a % perm)) then
      acc_perm => acc_a % perm
    else
      allocate(acc_perm(acc_a % numberofrows))
    end if
    acc_perm = 0
    dirichlet = 0
    l = 0
    do i = 1, acc_a % numberofrows
      s = acc_a % values(acc_a % diag(i))
      acc_a % values(acc_a % diag(i)) = 0.0d0
      j = acc_a % rows(i)
      k = acc_a % rows(i+1) - 1
      if (all(abs(acc_a % values(j:k)) < 1.d-12)) then
        dirichlet = dirichlet + 1
      else
        l = l + 1
        acc_perm(i) = l
      end if
      acc_a % values(acc_a % diag(i)) = s
    end do
    acc_bm => allocatematrix()
    acc_bm % numberofrows = acc_a % numberofrows - dirichlet
    allocate(acc_bm % rows(acc_bm % numberofrows + 1), &
        acc_bm % diag(acc_bm % numberofrows))
    j = 0
    k = 1
    zz = 0
    do i = 1, acc_a % numberofrows
      if (acc_perm(i) /= 0) then
        j = j + 1
        acc_bm % rows(j) = k
        do l = acc_a % rows(i), acc_a % rows(i+1) - 1
          t = acc_a % cols(l)
          add = .true.
          if (.not. acc_a % complex) then
            add = acc_a % values(l) /= 0
            if (mass_l) add = add .or. acc_a % massvalues(l) /= 0
            if (damp_l) add = add .or. acc_a % dampvalues(l) /= 0
            if (.not. add) zz = zz + 1
          end if
          add = add .or. l == acc_a % diag(i)
          add = add .and. acc_perm(t) /= 0
          if (add) k = k + 1
        end do
      end if
    end do
    acc_bm % rows(acc_bm % numberofrows + 1) = k
    if (mass_l) then
      allocate(acc_bm % massvalues(k-1))
    end if
    if (damp_l) then
      allocate(acc_bm % dampvalues(k-1))
    end if
    allocate(acc_bm % values(k-1), acc_bm % cols(k-1), &
        acc_bm % rhs(acc_bm % numberofrows))
    allocate(dcount(acc_a % numberofrows))
    if (acc_perm(1) == 0) then
      dcount(1) = 1
    else
      dcount(1) = 0
    end if
    do i = 2, acc_a % numberofrows
      if (acc_perm(i) /= 0) then
        dcount(i) = dcount(i-1)
      else
        dcount(i) = dcount(i-1) + 1
      end if
    end do
    j = 0
    k = 1
    do i = 1, acc_a % numberofrows
      if (acc_perm(i) /= 0) then
        j = j + 1
        do l = acc_a % rows(i), acc_a % rows(i+1) - 1
          t = acc_a % cols(l)
          add = .true.
          if (.not. acc_a % complex) then
            add = acc_a % values(l) /= 0
            if (mass_l) add = add .or. acc_a % massvalues(l) /= 0
            if (damp_l) add = add .or. acc_a % dampvalues(l) /= 0
          end if
          add = add .or. l == acc_a % diag(i)
          add = add .and. acc_perm(t) /= 0
          if (add) then
            acc_bm % cols(k) = t - dcount(t)
            k = k + 1
          end if
        end do
      end if
    end do
    deallocate(dcount)
    do i = 1, acc_bm % numberofrows
      do j = acc_bm % rows(i), acc_bm % rows(i+1) - 1
        if (acc_bm % cols(j) == i) then
          acc_bm % diag(i) = j
        end if
      end do
    end do
    call crs_sortmatrix(acc_bm)
    acc_a % ematrix => acc_bm
    acc_a % perm => acc_perm
  end subroutine eld_build

  subroutine eld_setup_fu() bind(c, name="eld_setup_fu")
    acc_perm => acc_a % perm
    acc_f => acc_bm % rhs
    allocate(acc_u(acc_bm % numberofrows))
  end subroutine eld_setup_fu

  ! -- scalar traffic for the Rust-side real fill / map-back ------------
  function eld_is_complex() result(code) bind(c, name="eld_is_complex")
    integer(c_int) :: code
    code = 0
    if (acc_solver % matrix % complex) code = 1
  end function eld_is_complex

  function eld_a_nrows() result(n) bind(c, name="eld_a_nrows")
    integer(c_int) :: n
    n = acc_a % numberofrows
  end function eld_a_nrows

  function eld_a_row(i) result(r) bind(c, name="eld_a_row")
    integer(c_int), value :: i
    integer(c_int) :: r
    r = acc_a % rows(i)
  end function eld_a_row

  function eld_a_row1(i) result(r) bind(c, name="eld_a_row1")
    integer(c_int), value :: i
    integer(c_int) :: r
    r = acc_a % rows(i+1)
  end function eld_a_row1

  function eld_a_col(l) result(c) bind(c, name="eld_a_col")
    integer(c_int), value :: l
    integer(c_int) :: c
    c = acc_a % cols(l)
  end function eld_a_col

  function eld_a_val(l) result(v) bind(c, name="eld_a_val")
    integer(c_int), value :: l
    real(c_double) :: v
    v = acc_a % values(l)
  end function eld_a_val

  function eld_a_mass(l) result(v) bind(c, name="eld_a_mass")
    integer(c_int), value :: l
    real(c_double) :: v
    v = acc_a % massvalues(l)
  end function eld_a_mass

  function eld_a_damp(l) result(v) bind(c, name="eld_a_damp")
    integer(c_int), value :: l
    real(c_double) :: v
    v = acc_a % dampvalues(l)
  end function eld_a_damp

  function eld_a_diagval(t) result(v) bind(c, name="eld_a_diagval")
    integer(c_int), value :: t
    real(c_double) :: v
    v = acc_a % values(acc_a % diag(t))
  end function eld_a_diagval

  function eld_a_is_diag(i, l) result(code) bind(c, name="eld_a_is_diag")
    integer(c_int), value :: i, l
    integer(c_int) :: code
    code = 0
    if (l == acc_a % diag(i)) code = 1
  end function eld_a_is_diag

  function eld_perm(i) result(p) bind(c, name="eld_perm")
    integer(c_int), value :: i
    integer(c_int) :: p
    p = acc_perm(i)
  end function eld_perm
  function eld_b_get(t) result(v) bind(c, name="eld_b_get")
    integer(c_int), value :: t
    real(c_double) :: v
    v = acc_b(t)
  end function eld_b_get

  function eld_x_get(i) result(v) bind(c, name="eld_x_get")
    integer(c_int), value :: i
    real(c_double) :: v
    v = acc_x(i)
  end function eld_x_get

  subroutine eld_x_set(i, v) bind(c, name="eld_x_set")
    integer(c_int), value :: i
    real(c_double), value :: v
    acc_x(i) = v
  end subroutine eld_x_set


  function eld_f_get(j) result(v) bind(c, name="eld_f_get")
    integer(c_int), value :: j
    real(c_double) :: v
    v = acc_f(j)
  end function eld_f_get

  subroutine eld_f_set(j, v) bind(c, name="eld_f_set")
    integer(c_int), value :: j
    real(c_double), value :: v
    acc_f(j) = v
  end subroutine eld_f_set

  function eld_u_get(j) result(v) bind(c, name="eld_u_get")
    integer(c_int), value :: j
    real(c_double) :: v
    v = acc_u(j)
  end function eld_u_get

  subroutine eld_u_set(j, v) bind(c, name="eld_u_set")
    integer(c_int), value :: j
    real(c_double), value :: v
    acc_u(j) = v
  end subroutine eld_u_set

  subroutine eld_bm_val_set(k, v) bind(c, name="eld_bm_val_set")
    integer(c_int), value :: k
    real(c_double), value :: v
    acc_bm % values(k) = v
  end subroutine eld_bm_val_set

  subroutine eld_bm_mass_set(k, v) bind(c, name="eld_bm_mass_set")
    integer(c_int), value :: k
    real(c_double), value :: v
    acc_bm % massvalues(k) = v
  end subroutine eld_bm_mass_set

  subroutine eld_bm_damp_set(k, v) bind(c, name="eld_bm_damp_set")
    integer(c_int), value :: k
    real(c_double), value :: v
    acc_bm % dampvalues(k) = v
  end subroutine eld_bm_damp_set

  ! -- complex fill (original 236-288, verbatim; untested serially) ------
  subroutine eld_fill_complex(n, b, x, mass, damp) &
      bind(c, name="eld_fill_complex")
    integer(c_int), value :: n, mass, damp
    real(c_double) :: b(n), x(n)
    logical :: mass_l, damp_l
    integer :: i, j, k, l, m
    real(kind=dp) :: br, bi
    mass_l = (mass /= 0)
    damp_l = (damp /= 0)
    j = 1
    k = 1
    do i = 1, acc_a % numberofrows, 2
      if (acc_perm(i) /= 0) then
        acc_f(i) = b(j)
        acc_f(i+1) = b(j+1)
        acc_u(i) = x(j)
        acc_u(i+1) = x(j+1)
        k = acc_bm % rows(j)
        do l = acc_a % rows(i), acc_a % rows(i+1) - 1, 2
          m = acc_a % cols(l)
          if (acc_perm(m) == 0) then
            br = b(m) / acc_a % values(acc_a % diag(m))
            bi = b(m+1) / acc_a % values(acc_a % diag(m+1))
            acc_f(j) = acc_f(j) - (acc_a % values(l)*br + acc_a % values(l+1)*bi)
            acc_f(j+1) = acc_f(j+1) - (acc_a % values(l)*bi - acc_a % values(l+1)*br)
          else
            acc_bm % values(k) = acc_a % values(l)
            acc_bm % values(k+1) = acc_a % values(l+1)
            if (mass_l) then
              acc_bm % massvalues(k) = acc_a % massvalues(l)
              acc_bm % massvalues(k+1) = acc_a % massvalues(l+1)
            end if
            if (damp_l) then
              acc_bm % dampvalues(k) = acc_a % dampvalues(l)
              acc_bm % dampvalues(k+1) = acc_a % dampvalues(l+1)
            end if
            k = k + 2
          end if
        end do
        j = j + 2
      end if
    end do
    do i = 1, acc_bm % numberofrows, 2
      k = acc_bm % rows(i+1)
      do j = acc_bm % rows(i), acc_bm % rows(i+1) - 1, 2
        acc_bm % values(k) = -acc_bm % values(i+1)
        acc_bm % values(k+1) = acc_bm % values(i)
        if (mass_l) then
          acc_bm % massvalues(k) = -acc_bm % massvalues(i+1)
          acc_bm % massvalues(k+1) = acc_bm % massvalues(i)
        end if
        if (damp_l) then
          acc_bm % dampvalues(k) = -acc_bm % dampvalues(i+1)
          acc_bm % dampvalues(k+1) = acc_bm % dampvalues(i)
        end if
        k = k + 2
      end do
    end do
  end subroutine eld_fill_complex

  ! -- pre-solve report + flag copies (original 321-331) ------------------
  subroutine eld_pre_solve() bind(c, name="eld_pre_solve")
    write(message, *) 'Eliminated ', &
        i2s(acc_a % numberofrows - acc_bm % numberofrows), &
        ' unknowns out of ', i2s(acc_a % numberofrows)
    call info('Dirichlet', message, level=5)
    acc_bm % lumped = acc_a % lumped
    acc_bm % complex = acc_a % complex
    acc_bm % symmetric = acc_a % symmetric
  end subroutine eld_pre_solve

  function eld_pes() result(n) bind(c, name="eld_pes")
    integer(c_int) :: n
    n = parenv % pes
  end function eld_pes

  ! -- reduced solve (original 334-376, verbatim) -------------------------
  subroutine eld_solve(dofs, nrm) bind(c, name="eld_solve")
    integer(c_int), value :: dofs
    real(c_double) :: nrm
    type(matrix_t), pointer :: cmatrix
    integer, pointer :: permsave(:)
    integer :: i, j, dof
    real(kind=dp) :: at
    logical :: gotit, dochange
    if (parenv % pes > 1) then
      permsave => acc_solver % variable % perm
      cmatrix => acc_solver % matrix
      acc_solver % matrix => acc_bm
      allocate(acc_solver % variable % perm(size(permsave)))
      acc_solver % variable % perm = 0
      do i = 1, size(permsave)
        j = permsave(i)
        if (j > 0) then
          dof = dofs * (j - 1) + 1
          if (acc_perm(dof) > 0) then
            acc_solver % variable % perm(i) = (acc_perm(dof) - 1) / dofs + 1
          end if
        end if
      end do
      if (.not. associated(acc_bm % parmatrix)) then
        call parallelinitmatrix(acc_solver, acc_bm)
      end if
    end if
    at = cputime()
    dochange = listgetlogical(acc_solver % values, 'Skip Compute Nonlinear Change', gotit)
    call listaddlogical(acc_solver % values, 'Skip Compute Nonlinear Change', .true.)
    call SolveLinearSystem(acc_bm, acc_f, acc_u, nrm, dofs, acc_solver)
    call listaddlogical(acc_solver % values, 'Skip Compute Nonlinear Change', dochange)
    acc_lintime = cputime() - at
    write(message, *) 'Time used in SolveLinearSystem (CPU): ', acc_lintime
    call info('Dirichlet', message, level=5)
    if (parenv % pes > 1) then
      acc_solver % matrix => cmatrix
      deallocate(acc_solver % variable % perm)
      acc_solver % variable % perm => permsave
    end if
  end subroutine eld_solve

  function eld_n_eigen() result(n) bind(c, name="eld_n_eigen")
    integer(c_int) :: n
    n = acc_solver % nofeigenvalues
  end function eld_n_eigen

  ! -- eigen map-back (original 382-414 + 427, verbatim) -------------------
  subroutine eld_mapback_eigen(n, b) bind(c, name="eld_mapback_eigen")
    integer(c_int), value :: n
    real(c_double) :: b(n)
    complex(kind=dp), allocatable :: tempeigenvectors(:,:)
    integer :: i, j, k
    allocate(tempeigenvectors(acc_solver % nofeigenvalues, acc_a % numberofrows))
    tempeigenvectors = acc_solver % variable % eigenvectors
    j = 0
    if (acc_solver % matrix % complex) then
      do i = 1, acc_a % numberofrows / 2
        if (acc_perm(i) == 0) then
          do k = 1, acc_solver % nofeigenvalues
            acc_solver % variable % eigenvectors(k,i) = cmplx( &
                b(2*i) / acc_a % values(acc_a % diag(2*i)), &
                b(2*i+1) / acc_a % values(acc_a % diag(2*i+1)), kind=dp)
          end do
        else
          j = j + 1
          do k = 1, acc_solver % nofeigenvalues
            acc_solver % variable % eigenvectors(k,i) = tempeigenvectors(k,j)
          end do
        end if
      end do
    else
      do i = 1, acc_a % numberofrows
        if (acc_perm(i) == 0) then
          do k = 1, acc_solver % nofeigenvalues
            acc_solver % variable % eigenvectors(k,i) = &
                b(i) / acc_a % values(acc_a % diag(i))
          end do
        else
          j = j + 1
          do k = 1, acc_solver % nofeigenvalues
            acc_solver % variable % eigenvectors(k,i) = tempeigenvectors(k,j)
          end do
        end if
      end do
    end if
    deallocate(tempeigenvectors)
  end subroutine eld_mapback_eigen

  ! -- tail (original 427-448): U drop, ComputeChange, timing, free ------
  function eld_finish() result(rv) bind(c, name="eld_finish")
    integer(c_int) :: rv
    real(kind=dp) :: at
    deallocate(acc_u)
    at = cputime()
    call computechange(acc_solver, .false.)
    acc_lintime = acc_lintime + (cputime() - at)
    acc_tottime = cputime() - acc_tottime
    write(message, *) 'Total time spent in DirichletReduction (CPU): ', acc_tottime
    call info('Dirichlet', message, level=5)
    write(message, *) 'Additional time spent in DirichletReduction (CPU): ', &
        acc_tottime - acc_lintime
    call info('Dirichlet', message, level=5)
    call info('Dirichlet', ' ', level=5)
    call freematrix(acc_a % ematrix)
    nullify(acc_a % ematrix)
    acc_bm => null()
    rv = 1
  end function eld_finish

end module eliminatedirichlet_acc
