! T021-A EnergyRelease port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/EnergyRelease.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments. All Fortran-object
! traffic — TYPE(...) handles, assumed-shape descriptors, hidden
! lengths, list-directed PRINT reports — stays inside Fortran, so the
! Rust side is pure C calls. The wrapper module mirrors EnergyRelease's
! own locals (current element, SAVE Nodes state, SAVE work arrays sized
! at first-mesh sizing with identical lifetimes). All floating-point
! work (the per-integration-point tensor assembly, `MATMUL`s, `SUM`
! reductions) and every report stays in Fortran (same compiler), so
! bit-identity holds by construction wherever the call sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module energyrelease_acc
  use iso_c_binding
  use defutils
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(element_t), pointer, save :: acc_elem => null()
  ! LocalReleaseRate context (mirrors its SAVE Nodes).
  type(nodes_t), target, save :: acc_nodes
  type(gaussintegrationpoints_t), save :: acc_ip
  ! SAVE work arrays (first-mesh sizing, kept across calls).
  real(c_double), allocatable, save :: acc_displ(:, :), acc_stress(:, :), &
      acc_shape(:, :)
  logical, save :: acc_allocdone = .false.

contains

  ! -- entry state ----------------------------------------------------
  subroutine ere_begin(h_solver) bind(c, name="ere_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine ere_begin

  ! -- ReleaseRateSolver body -------------------------------------------
  subroutine ere_run() bind(c, name="ere_run")
    integer :: n, t, istat, nelem
    real(c_double) :: localgtheta, gtheta
    if (.not. acc_allocdone) then
      n = acc_solver % mesh % maxelementnodes
      allocate(acc_displ(3, n), acc_stress(6, n), &
          acc_shape(3, n), stat=istat)
      if (istat /= 0) then
        call Fatal('PoissonSolve', 'Memory allocation error.')
      end if
      acc_allocdone = .true.
    end if
    gtheta = 0.0_c_double
    nelem = acc_solver % numberofactiveelements
    do t = 1, nelem
      acc_elem => GetActiveElement(t)
      n = GetElementNOFNodes()
      call GetVectorLocalSolution(acc_displ, 'True Displacement')
      call GetVectorLocalSolution(acc_stress, 'Stress')
      call GetVectorLocalSolution(acc_shape, 'Shape Displacement')
      call ere_local_rate(localgtheta, acc_shape, acc_displ, acc_stress, n)
      gtheta = gtheta + localgtheta
    end do
    print *
    print *, '*******************************'
    print *, 'Calculated Fracture Parameters:'
    print *, '-------------------------------'
    print *, '2G = ', gtheta
    print *, '*******************************'
    print *
  end subroutine ere_run

  ! -- LocalReleaseRate internal (exact transcription) ------------------
  subroutine ere_local_rate(localgtheta, localprop, localdispl, localstress, n)
    real(c_double) :: localgtheta
    real(c_double) :: localprop(:, :)
    real(c_double) :: localdispl(:, :), localstress(:, :)
    integer :: n
    real(c_double) :: basis(n), dbasisdx(n, 3), detj, loadatip
    integer :: t, i, j
    logical :: stat
    real(c_double) :: displacement(3), graddisplacement(3, 3)
    real(c_double) :: propagationshape(3), gradpropagationshape(3, 3)
    real(c_double) :: r_tensor(3, 3), s_tensor(3, 3), stresstensor(3, 3)
    real(c_double) :: divpropagationshape
    call GetElementNodes(acc_nodes)
    localgtheta = 0.0_c_double
    acc_ip = GaussPoints(acc_elem)
    do t = 1, acc_ip % n
      stat = ElementInfo(acc_elem, acc_nodes, acc_ip % u(t), acc_ip % v(t), &
          acc_ip % w(t), detj, basis, dbasisdx)
      do i = 1, 3
        displacement(i) = sum(basis(1:n) * localdispl(i, 1:n))
        propagationshape(i) = sum(basis(1:n) * localprop(i, 1:n))
        do j = 1, 3
          graddisplacement(i, j) = sum(dbasisdx(1:n, j) * localdispl(i, 1:n))
          gradpropagationshape(i, j) = sum(dbasisdx(1:n, j) * localprop(i, 1:n))
        end do
      end do
      divpropagationshape = 0.0_c_double
      do i = 1, 3
        divpropagationshape = divpropagationshape + gradpropagationshape(i, i)
      end do
      stresstensor(1, 1) = sum(basis(1:n) * localstress(1, 1:n))
      stresstensor(2, 2) = sum(basis(1:n) * localstress(2, 1:n))
      stresstensor(3, 3) = sum(basis(1:n) * localstress(3, 1:n))
      stresstensor(1, 2) = sum(basis(1:n) * localstress(4, 1:n))
      stresstensor(2, 3) = sum(basis(1:n) * localstress(5, 1:n))
      stresstensor(1, 3) = sum(basis(1:n) * localstress(6, 1:n))
      stresstensor(2, 1) = sum(basis(1:n) * localstress(4, 1:n))
      stresstensor(3, 2) = sum(basis(1:n) * localstress(5, 1:n))
      stresstensor(3, 1) = sum(basis(1:n) * localstress(6, 1:n))
      r_tensor = -matmul(graddisplacement, gradpropagationshape)
      s_tensor = matmul(stresstensor, transpose(gradpropagationshape)) &
          - divpropagationshape * stresstensor
      localgtheta = localgtheta + (sum(s_tensor * graddisplacement) &
          - sum(r_tensor * stresstensor)) * acc_ip % s(t) * detj
    end do
  end subroutine ere_local_rate

end module energyrelease_acc
