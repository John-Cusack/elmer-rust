! Probes port: statement-level BIND(C) wrappers for the Elmer calls in
! fem/src/modules/Probes.F90, compiled against the built .mod files of
! the target configure (A2 shipped mechanism, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls.
!
! Faithfulness notes (mirrors VolumeProbe statement for statement):
! - SAVE Nodes lives here with the original's lifetime (filled per
!   active element, consumed before the next fill, as in the original).
! - The body-force probe-name string chopping + ValueHandle init block
!   is transcribed verbatim into prb_parse (same comparisons, same
!   order, same ListInitElementKeyword calls); the parse state
!   (names/biases/handles) lives in SAVE module data but is rebuilt
!   from scratch on every entry call, exactly as the original's
!   per-call locals are (previous allocation freed first — no
!   observable difference).
! - The element/gauss/probe accumulation order, the `gausspoint=ngp`
!   argument, the `(ParEnv % PEs > 1) .AND. (.NOT. SingleMesh)` guard,
!   and the `'res: volume probe ' // trim(name)` store loop all stay in
!   Fortran (prb_probe_value / prb_store), so bit-identity holds by
!   construction wherever the call sequence matches.
! - The no-names branch (keyword absent) yields zero probes here; the
!   reference build faults on SIZE() of the never-allocated handle in
!   that case (triage: SIGSEGV, exit 139), so no valid differential
!   exists for it — the proof exercises single-name vs multi-name
!   configs instead.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module probes_acc
  use iso_c_binding
  use defutils
  use messages
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(element_t), pointer, save :: acc_elem => null()
  ! SAVE Nodes mirror (the original's internal SAVE Nodes).
  type(nodes_t), target, save :: acc_nodes
  ! Integration-point context (the original's per-element IP local).
  type(gaussintegrationpoints_t), save :: acc_ip
  integer(c_int), save :: acc_ngp = 0
  ! Basis scratch (the original's per-element automatic Basis(nd)).
  real(c_double), allocatable, save :: acc_basis(:)
  ! Probe-parse state (rebuilt on every entry call, like locals).
  character(len=MAX_NAME_LEN), save :: acc_names = ''
  integer, allocatable, save :: acc_biases(:)
  type(valuehandle_t), pointer, save :: acc_h(:) => null()
  integer(c_int), save :: acc_nprobes = 0

contains

  ! -- entry state ----------------------------------------------------
  subroutine prb_begin(h_model, h_solver) bind(c, name="prb_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine prb_begin

  subroutine prb_touch_mesh() bind(c, name="prb_touch_mesh")
    type(mesh_t), pointer :: m
    m => GetMesh()
  end subroutine prb_touch_mesh

  function prb_coord_dim() result(d) bind(c, name="prb_coord_dim")
    integer(c_int) :: d
    d = CoordinateSystemDimension()
  end function prb_coord_dim

  function prb_nof_active() result(n) bind(c, name="prb_nof_active")
    integer(c_int) :: n
    n = GetNOFActive(acc_solver)
  end function prb_nof_active

  ! -- probe-name parse (verbatim transcription of the Body force -----
  ! -- probe names block; returns numsubs-1, 0 when keyword absent) ---
  function prb_parse() result(nprobes) bind(c, name="prb_parse")
    integer(c_int) :: nprobes
    character(len=MAX_NAME_LEN) :: probe_names
    logical :: found
    integer :: bias, numsubs, lenstr, counter, addbias

    if (associated(acc_h)) deallocate(acc_h)
    if (allocated(acc_biases)) deallocate(acc_biases)
    acc_names = ''
    acc_nprobes = 0

    probe_names = ListGetString(acc_solver % values, 'Body force probe names', found)

    if (found) then
      acc_names = probe_names
      lenstr = len(trim(acc_names))
      store: do while (.true.)

        bias = 1
        numsubs = 0
        counter = 0 ! to force exit at maxiter

        remove_left_ws: do while (acc_names(bias:bias) == " ")
          bias = bias + 1
        end do remove_left_ws
        if (bias >= lenstr) then ! empty string full of ws
          exit store
        end if
        numsubs = 1 ! at least 1 substring

        if (allocated(acc_biases)) then
          acc_biases(1) = bias
        end if

        scanstr: do while (.true.)
          remove_left_nws: do while (acc_names(bias:bias) /= " " .and. bias <= lenstr)
            bias = bias+1
          end do remove_left_nws

          if (bias >= lenstr) then ! no more nws
            numsubs = numsubs + 1
            if (allocated(acc_biases)) then
              acc_biases(numsubs) = lenstr+1 ! last one is length of the string
            end if
            exit scanstr
          end if

          hog_whitespace: do while (acc_names(bias:bias) == " " .and. bias <= lenstr)
            bias = bias + 1
          end do hog_whitespace

          numsubs = numsubs + 1
          if (allocated(acc_biases)) then
            acc_biases(numsubs) = bias
          end if

          ! sanity quit
          counter = counter + 1
          if (counter .ge. 10) then
            exit scanstr
          end if
        end do scanstr

        if (.not. allocated(acc_biases)) then
          allocate(acc_biases(numsubs))
          allocate(acc_h(numsubs-1))
        else
          exit store
        end if
      end do store

      do counter = 1, numsubs-1
        call ListInitElementKeyword( acc_h(counter), 'Body Force', &
            trim(acc_names(acc_biases(counter):acc_biases(counter+1)-1)) )
      end do
      acc_nprobes = numsubs - 1
    end if

    nprobes = acc_nprobes
  end function prb_parse

  ! -- active-element context ------------------------------------------
  subroutine prb_begin_active(t, n, nd, nb, ngp) bind(c, name="prb_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, nb, ngp
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
    nb = GetElementNOFBDOFs()
    acc_ip = GaussPoints(acc_elem)
    ngp = acc_ip % n
    acc_ngp = ngp
  end subroutine prb_begin_active

  subroutine prb_nodes() bind(c, name="prb_nodes")
    call GetElementNodes(acc_nodes)
  end subroutine prb_nodes

  subroutine prb_ip_point(q, u, v, w, s) bind(c, name="prb_ip_point")
    integer(c_int), value :: q
    real(c_double) :: u, v, w, s
    u = acc_ip % u(q)
    v = acc_ip % v(q)
    w = acc_ip % w(q)
    s = acc_ip % s(q)
  end subroutine prb_ip_point

  subroutine prb_eleminfo(u, v, w, nd, detj, stat) bind(c, name="prb_eleminfo")
    real(c_double), value :: u, v, w
    integer(c_int), value :: nd
    real(c_double) :: detj
    integer(c_int) :: stat
    logical :: ok
    if (.not. allocated(acc_basis) .or. size(acc_basis) /= nd) then
      if (allocated(acc_basis)) deallocate(acc_basis)
      allocate(acc_basis(nd))
    end if
    ok = ElementInfo(acc_elem, acc_nodes, u, v, w, detj, acc_basis(1:nd))
    stat = merge(1, 0, ok)
  end subroutine prb_eleminfo

  ! -- per-probe integrand (gausspoint=ngp, exactly as original) -------
  subroutine prb_probe_value(m, val, found) bind(c, name="prb_probe_value")
    integer(c_int), value :: m
    real(c_double) :: val
    integer(c_int) :: found
    logical :: f
    val = ListGetElementReal(acc_h(m), acc_basis, acc_elem, f, gausspoint=acc_ngp)
    found = merge(1, 0, f)
  end subroutine prb_probe_value

  ! -- parallel reduction + result store (verbatim tail) ---------------
  subroutine prb_store(n, vals) bind(c, name="prb_store")
    integer(c_int), value :: n
    real(c_double) :: vals(n)
    logical :: parallel
    integer :: t
    if (n <= 0) return
    parallel = (parenv % pes > 1) .and. (.not. acc_solver % mesh % singlemesh)
    if (parallel) then
      do t = 1, n
        vals(t) = parallelreduction(vals(t))
      end do
    end if
    do t = 1, n
      call ListAddConstReal( acc_model % simulation, 'res: volume probe ' // &
          trim(acc_names(acc_biases(t):acc_biases(t+1)-1)), vals(t) )
    end do
  end subroutine prb_store

  subroutine prb_fatal_alloc() bind(c, name="prb_fatal_alloc")
    call Fatal('VolumeProbe', 'Memory allocation error.')
  end subroutine prb_fatal_alloc

end module probes_acc
