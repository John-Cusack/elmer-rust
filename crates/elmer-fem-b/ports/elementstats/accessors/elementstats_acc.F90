! T021-A ElementStats port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/ElementStats.F90, compiled against the
! built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments. All Fortran-object
! traffic — TYPE(...) handles, POINTER components, assumed-shape
! descriptors, hidden lengths, list-directed PRINT reports — stays inside
! Fortran, so the Rust side is pure C calls. The wrapper module mirrors
! ElementStats' own locals (current element, Nodes SAVE state per
! internal study, entry-level flags) with identical lifetimes. All
! floating-point work (the three element studies, `CEILING` binning,
! `MIN`/`MAX`/`SUM` reductions, `ACOS`/`SQRT`) and every report stays in
! Fortran (same compiler), so bit-identity holds by construction wherever
! the call sequence matches. The parallel-reduction branch is transcribed
! as written; in serial builds `IsParallel` is false by the same guard as
! the original.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module elementstats_acc
  use iso_c_binding
  use defutils
  use elementdescription, only: GetEdgeMap
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(valuelist_t), pointer, save :: acc_params => null()
  type(mesh_t), pointer, save :: acc_mesh => null()
  type(element_t), pointer, save :: acc_elem => null()
  ! Internal studies' SAVE Nodes (one set per study, as in the original).
  type(nodes_t), target, save :: acc_nodes_size, acc_nodes_ratio, acc_nodes_skew
  ! Entry-level state (mirrors the solver body's locals).
  logical, save :: acc_isparallel = .false.
  logical, save :: acc_docohorts = .false.
  logical, save :: acc_applywarn = .false.
  logical, save :: acc_warnfatal = .false.
  logical, save :: acc_stopnans = .false.
  integer(c_int), save :: acc_nocohorts = 10
  integer(c_int), save :: acc_bcmax = 1
  real(c_double), save :: acc_warnsize = 0.0_c_double
  integer(c_int), save :: acc_nancount(3) = 0
  integer(c_int), allocatable, save :: acc_cohortcount(:, :)

contains

  ! -- entry state ----------------------------------------------------
  subroutine est_begin(h_model, h_solver) bind(c, name="est_begin")
    type(c_ptr), value :: h_model, h_solver
    integer :: ierr
    logical :: found
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
    call Info('ElementStats', '------------------------------', level=6)
    call Info('ElementStats', 'Studying elements of the mesh')
    acc_params => GetSolverParams()
    acc_mesh => GetMesh()
    acc_isparallel = GetLogical(acc_params, 'Reduce Parallel', found) &
        .and. (parenv % pes > 1)
    acc_docohorts = GetLogical(acc_params, 'Create Histogram', found)
    if (acc_docohorts) then
      acc_nocohorts = GetInteger(acc_params, 'Histogram Intervals', found)
      if (.not. found) acc_nocohorts = 10
      allocate(acc_cohortcount(acc_nocohorts, 3), stat=ierr)
      if (ierr /= 0) call Fatal('ElementStats', 'Memory allocation error.')
      acc_cohortcount = 0
    end if
    acc_bcmax = 1
    if (GetLogical(acc_params, 'Skip Boundaries', found)) acc_bcmax = 0
    acc_warnsize = GetConstReal(acc_params, 'Element Warn Size', acc_applywarn)
    acc_warnfatal = GetLogical(acc_params, 'Stop Below Warn Size', found)
    acc_stopnans = GetLogical(acc_params, 'Stop At Nans', found)
    acc_nancount = 0
  end subroutine est_begin

  function est_bc_max() result(n) bind(c, name="est_bc_max")
    integer(c_int) :: n
    n = acc_bcmax
  end function est_bc_max

  ! -- one is_bc pass (bulk, then boundary) ----------------------------
  subroutine est_pass(is_bc) bind(c, name="est_pass")
    integer(c_int), value :: is_bc
    integer :: t, t0, n, m, firstelem, lastelem, elmntcnt, &
        localelmntcnt, operno, i, j, ierr
    real(c_double) :: minf(3), maxf(3), avef(3), f, df, &
        localminf(3), localmaxf(3), localavef(3)
    integer :: elementcount(1000), localelementcount(1000)
    character(len=MAX_NAME_LEN) :: str, opername
    logical :: doingcohorts, havenan
    if (is_bc == 0) then
      t0 = 0
      firstelem = 1
      lastelem = acc_mesh % numberofbulkelements
      str = 'bulk'
    else
      t0 = acc_mesh % numberofbulkelements
      firstelem = 1
      lastelem = acc_mesh % numberofboundaryelements
      str = 'boundary'
    end if
    if (firstelem > lastelem) return
    call Info('ElementStats', 'Going through ' // trim(str) // ' elements')
    localminf = huge(localminf)
    localmaxf = -huge(localmaxf)
    localavef = 0.0_c_double
    doingcohorts = .false.
    localelementcount = 0
    elementcount = 0
    do ! cohort re-pass (the original's GOTO 100)
      do t = firstelem, lastelem
        acc_elem => acc_mesh % elements(t0 + t)
        acc_model % currentelement => acc_elem
        n = GetElementNOFNodes(acc_elem)
        m = GetElementCode(acc_elem)
        if (.not. doingcohorts) then
          localelementcount(m) = localelementcount(m) + 1
        end if
        do operno = 1, 3
          select case (operno)
          case (1)
            call est_size_study(acc_elem, n, f, havenan)
          case (2)
            call est_skew_study(acc_elem, n, f, havenan)
          case (3)
            call est_ratio_study(acc_elem, n, f, havenan)
          end select
          if (havenan) then
            acc_nancount(operno) = acc_nancount(operno) + 1
            cycle
          end if
          if (.not. doingcohorts) then
            localminf(operno) = min(f, localminf(operno))
            localmaxf(operno) = max(f, localmaxf(operno))
            localavef(operno) = f + localavef(operno)
          else
            i = ceiling(acc_nocohorts * (f - localminf(operno)) &
                / (localmaxf(operno) - localminf(operno)))
            i = max(1, min(i, int(acc_nocohorts)))
            acc_cohortcount(i, operno) = acc_cohortcount(i, operno) + 1
          end if
        end do
      end do
      localelmntcnt = lastelem - firstelem + 1
      if (acc_isparallel) then
        acc_nancount(1) = ParallelReduction(acc_nancount(1))
        acc_nancount(2) = ParallelReduction(acc_nancount(2))
        call MPI_ALLREDUCE(localminf, minf, 3, &
            MPI_DOUBLE_PRECISION, MPI_MIN, ELMER_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(localmaxf, maxf, 3, &
            MPI_DOUBLE_PRECISION, MPI_MAX, ELMER_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(localavef, avef, 3, &
            MPI_DOUBLE_PRECISION, MPI_SUM, ELMER_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(localelmntcnt, elmntcnt, 1, &
            MPI_INTEGER, MPI_SUM, ELMER_COMM_WORLD, ierr)
      else
        minf = localminf
        maxf = localmaxf
        avef = localavef
        elmntcnt = localelmntcnt
      end if
      if (acc_nancount(1) > 0) then
        call Info('ElementStats', &
            '***** Number of NaN sized elements: ' // i2s(acc_nancount(1)), level=3)
      end if
      if (acc_nancount(3) > 0) then
        call Info('ElementStats', &
            '***** Number of extreme element edge ratios: ' // i2s(acc_nancount(3)), &
            level=3)
      end if
      if (acc_docohorts .and. .not. doingcohorts) then
        doingcohorts = .true.
        cycle
      end if
      exit
    end do
    if (lastelem > firstelem) then
      avef = avef / elmntcnt
    end if
    do operno = 1, 3
      select case (operno)
      case (1)
        opername = 'Element size'
      case (2)
        opername = 'Element skew (degs)'
      case (3)
        opername = 'Element ratio'
      end select
      if (parenv % mype == 0) then
        print *, 'Statistics for mesh operator: ' // trim(opername) // ' on ', str
        print *, 'Minimum values ' // trim(opername) // ': ', minf(operno)
        print *, 'Maximum values ' // trim(opername) // ': ', maxf(operno)
        print *, 'Average values ' // trim(opername) // ': ', avef(operno)
        if (operno == 1) then
          print *, 'Ratio of element sizes: ', maxf(operno) / minf(operno)
        end if
        print *, ' '
        if (acc_docohorts) then
          print *, 'Histogram:'
          df = (maxf(operno) - minf(operno)) / acc_nocohorts
          do i = 1, acc_nocohorts
            if (acc_cohortcount(i, operno) > 0) then
              print *, i, ': ', acc_cohortcount(i, operno), ' (', &
                  minf(operno) + (i - 1) * df, ' to ', minf(operno) + i * df, ')'
            end if
          end do
        end if
      end if
    end do
    if (acc_applywarn .and. (minf(1) <= acc_warnsize) .and. parenv % mype == 0) then
      write(message, *) 'Minimum Element size ', minf(1), &
          ' below given threshold', acc_warnsize
      if (acc_warnfatal) then
        call Fatal('ElementStats', message)
      else
        call Warn('ElementStats', message)
      end if
    end if
  end subroutine est_pass

  subroutine est_finish() bind(c, name="est_finish")
    if (allocated(acc_cohortcount)) deallocate(acc_cohortcount)
    if (acc_stopnans .and. sum(acc_nancount) > 0) then
      call Fatal('ElementStats', 'Mesh has elements of invalid sized or proportion')
    end if
    call Info('ElementStats', '------------------------------', level=6)
  end subroutine est_finish

  ! -- internal studies (exact transcriptions) --------------------------
  subroutine est_size_study(element, n, elemsize, havenan)
    type(element_t), target :: element
    integer :: n
    real(c_double) :: elemsize
    logical :: havenan
    real(c_double) :: weight
    real(c_double) :: basis(n), detj
    logical :: stat
    integer :: t
    type(gaussintegrationpoints_t) :: ip
    call GetElementNodes(acc_nodes_size)
    ip = GaussPoints(element)
    elemsize = 0.0_c_double
    havenan = .false.
    do t = 1, ip % n
      stat = ElementInfo(element, acc_nodes_size, ip % u(t), ip % v(t), &
          ip % w(t), detj, basis)
      weight = ip % s(t) * detj
      elemsize = elemsize + weight
      if (.not. stat) exit
    end do
    if (.not. (elemsize == elemsize) .or. elemsize >= huge(elemsize) &
        .or. .not. stat) then
      havenan = .true.
      elemsize = 0.0_c_double
    end if
  end subroutine est_size_study

  subroutine est_ratio_study(element, n, q, havenan)
    type(element_t), target :: element
    integer :: n
    real(c_double) :: q
    logical :: havenan
    integer :: elemfamily
    integer, pointer :: edgemap(:, :)
    integer :: i1, i2, k
    real(c_double) :: ri(3), s2, s2min, s2max
    call GetElementNodes(acc_nodes_ratio)
    elemfamily = GetElementFamily(element)
    edgemap => GetEdgeMap(elemfamily)
    s2max = 0.0_c_double
    s2min = huge(s2min)
    havenan = .false.
    q = 0.0_c_double
    do k = 1, element % type % numberofedges
      i1 = edgemap(k, 1)
      i2 = edgemap(k, 2)
      ri(1) = acc_nodes_ratio % x(i2) - acc_nodes_ratio % x(i1)
      ri(2) = acc_nodes_ratio % y(i2) - acc_nodes_ratio % y(i1)
      ri(3) = acc_nodes_ratio % z(i2) - acc_nodes_ratio % z(i1)
      s2 = sqrt(sum(ri * ri))
      s2max = max(s2max, s2)
      s2min = min(s2min, s2)
    end do
    q = sqrt(s2max / s2min)
    if (.not. (q == q) .or. q >= huge(q)) then
      havenan = .true.
      q = 0.0_c_double
    end if
  end subroutine est_ratio_study

  subroutine est_skew_study(element, n, alpha, havenan)
    type(element_t), target :: element
    integer :: n
    real(c_double) :: alpha
    logical :: havenan
    integer :: elemfamily
    integer, pointer :: edgemap(:, :)
    integer :: noedgepairs, i1, i2, j1, j2, k
    integer :: edgepairs(24)
    real(c_double) :: beta, maxbeta, ri(3), rj(3), ri_len, rj_len
    call GetElementNodes(acc_nodes_skew)
    elemfamily = GetElementFamily(element)
    edgemap => GetEdgeMap(elemfamily)
    alpha = 0.0_c_double
    havenan = .false.
    select case (elemfamily)
    case (4)
      noedgepairs = 2
      edgepairs(1:2 * noedgepairs) = (/1, 3, 2, 4/)
    case (6)
      noedgepairs = 2
      edgepairs(1:2 * noedgepairs) = (/1, 3, 2, 4/)
    case (7)
      noedgepairs = 6
      edgepairs(1:2 * noedgepairs) = (/1, 4, 7, 8, 2, 5, 8, 9, 3, 6, 9, 7/)
    case (8)
      noedgepairs = 12
      edgepairs(1:2 * noedgepairs) = &
          (/1, 3, 2, 4, 5, 7, 6, 8, 1, 5, 9, 10, 2, 6, 10, 11, &
            3, 7, 11, 12, 4, 8, 12, 9/)
    case default
      return
    end select
    maxbeta = 0.0_c_double
    do k = 1, noedgepairs
      i1 = edgemap(edgepairs(2 * k - 1), 1)
      i2 = edgemap(edgepairs(2 * k - 1), 2)
      j1 = edgemap(edgepairs(2 * k), 1)
      j2 = edgemap(edgepairs(2 * k), 2)
      ri(1) = acc_nodes_skew % x(i2) - acc_nodes_skew % x(i1)
      ri(2) = acc_nodes_skew % y(i2) - acc_nodes_skew % y(i1)
      ri(3) = acc_nodes_skew % z(i2) - acc_nodes_skew % z(i1)
      rj(1) = acc_nodes_skew % x(j2) - acc_nodes_skew % x(j1)
      rj(2) = acc_nodes_skew % y(j2) - acc_nodes_skew % y(j1)
      rj(3) = acc_nodes_skew % z(j2) - acc_nodes_skew % z(j1)
      ri_len = sqrt(sum(ri * ri))
      rj_len = sqrt(sum(rj * rj))
      if (ri_len * rj_len > 0) then
        beta = acos(sum(ri * rj) / (ri_len * rj_len))
        beta = min(abs(beta), abs(PI - beta))
        maxbeta = max(beta, maxbeta)
      end if
    end do
    alpha = 180 * maxbeta / PI
  end subroutine est_skew_study

end module elementstats_acc
