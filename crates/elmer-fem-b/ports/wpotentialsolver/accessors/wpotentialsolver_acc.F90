! WPotentialSolver port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/WPotentialSolver.F90, compiled against the
! built .mod files of the target configure (solver-entry gate shape,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module mirrors the original's locals (current element/mesh/solver/model
! handles, Nodes SAVE state per internal procedure, SaveElementWSolution
! SAVE state) with identical lifetimes.
!
! Order-sensitive reductions (SUM/MATMUL), whole-array updates, string
! comparisons, and the Info call stay in Fortran, so bit-identity holds
! by construction wherever the call sequence matches. The solver-list
! surgery of Wsolve_Init0 (pointer-array copy, ListAdd* calls) lives in
! one coarse wrapper with the same statement order as the original.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module wpo_acc
  use iso_c_binding
  use types
  use defutils
  use circuitutils
  use mgdynmaterialutils
  use integration
  use eleminfo
  use messages
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(model_t), pointer, save :: acc_model => null()
  type(mesh_t), pointer, save :: acc_mesh => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(valuelist_t), pointer, save :: acc_bf => null()
  type(valuelist_t), pointer, save :: acc_bc => null()
  ! LocalMatrix context (mirrors its SAVE Nodes; IP is call-local there).
  type(nodes_t), target, save :: acc_nodes_lm
  ! AddElementWNormAndVolume context (mirrors its SAVE Nodes).
  type(nodes_t), target, save :: acc_nodes_wn
  ! SaveElementWSolution context (mirrors its SAVE Nodes plus the SAVE
  ! variable handles, First flag, and variable names).
  type(nodes_t), target, save :: acc_nodes_sv
  type(variable_t), pointer, save :: acc_wpotvar => null()
  type(variable_t), pointer, save :: acc_wvecvar => null()
  type(variable_t), pointer, save :: acc_jvecvar => null()
  logical, save :: acc_sv_first = .true.
  character(len=MAX_NAME_LEN), save :: acc_varname, acc_jvarname
  ! Wsolve_Init0 context (mirrors its SAVE visited counter).
  integer, save :: acc_visited = 0

contains

  ! -- Wsolve_Init0 ---------------------------------------------------
  subroutine wpo_init0_begin(h_model, h_solver) bind(c, name="wpo_init0_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine wpo_init0_begin

  ! Mirrors `visited = visited + 1` plus the Discontinuous Galerkin
  ! early return. Returns 1 when the caller must return immediately.
  function wpo_init0_enter() result(r) bind(c, name="wpo_init0_enter")
    integer(c_int) :: r
    logical :: found
    acc_visited = acc_visited + 1
    if (getlogical(getsolverparams(), 'Discontinuous Galerkin', found)) then
      r = 1
    else
      r = 0
    end if
  end function wpo_init0_enter

  ! The Init0 body proper: solver-list surgery creating the on-the-fly
  ! DG solver, statement for statement from the original.
  subroutine wpo_init0_extend() bind(c, name="wpo_init0_extend")
    character(len=MAX_NAME_LEN) :: varname, jvarname
    logical :: found
    integer, pointer :: active(:)
    integer :: mysolver, i, j, k, l, n, m, vdofs, soln
    type(valuelist_t), pointer :: solverparams
    type(solver_t), pointer :: solvers(:), psolver

    psolver => acc_solver
    do mysolver = 1, acc_model % numberofsolvers
      if (associated(psolver, acc_model % solvers(mysolver))) exit
    end do

    varname = getstring(getsolverparams(), 'Variable', found)
    if (.not. found) varname = 'W'

    n = acc_model % numberofsolvers
    do i = 1, acc_model % numberofequations
      active => listgetintegerarray(acc_model % equations(i) % values, &
          'Active Solvers', found)
      m = size(active)
      if (any(active == mysolver)) &
        call listaddintegerarray(acc_model % equations(i) % values, &
             'Active Solvers', m+1, [active, n+1])
    end do

    allocate(solvers(n+1))
    solvers(1:n) = acc_model % solvers
    solvers(n+1) % values => listallocate()
    solverparams => solvers(n+1) % values
    call listaddlogical(solverparams, 'Discontinuous Galerkin', .true.)
    solvers(n+1) % dg = .true.
    solvers(n+1) % procedure = c_null_funptr
    solvers(n+1) % activeelements => null()
    call listaddstring(solverparams, 'Exec Solver', 'never')
    call listaddlogical(solverparams, 'No Matrix', .true.)
    call listaddlogical(solverparams, 'Optimize Bandwidth', .false.)
    call listaddstring(solverparams, 'Equation', &
      'elementaladd'//i2s(acc_visited))
    call listaddstring(solverparams, 'Procedure', &
        'WPotentialSolver Wsolve_Dummy', .false.)
    call listaddstring(solverparams, 'Variable', '-nooutput '//trim(varname)//'_dummy')

    i = 1
    do while (.true.)
      if (listcheckpresent(solverparams, "Exported Variable "//i2s(i))) then
        i = i + 1
      else
        exit
      end if
    end do

    call listaddstring(solverparams, "Exported Variable "//i2s(i), &
          trim(varname)//" Potential")

    i = i + 1

    call listaddstring(solverparams, "Exported Variable "//i2s(i), &
            "-dg "//trim(varname)//" Vector E["//trim(varname)//" Vector E:3]")

    if (getlogical(getsolverparams(), 'Compute J Vector', found)) then
      jvarname = getstring(getsolverparams(), 'J Variable', found)
      if (.not. found) jvarname = 'J'
      i = i + 1

      call listaddstring(solverparams, "Exported Variable "//i2s(i), &
              "-dg "//trim(jvarname)//" Vector E["//trim(jvarname)//" Vector E:3]")
    end if

    deallocate(acc_model % solvers)
    acc_model % solvers => solvers
    acc_model % numberofsolvers = n + 1
  end subroutine wpo_init0_extend

  ! -- entry state ----------------------------------------------------
  subroutine wpo_begin(h_model, h_solver) bind(c, name="wpo_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine wpo_begin

  function wpo_assoc_matrix() result(r) bind(c, name="wpo_assoc_matrix")
    integer(c_int) :: r
    r = merge(1, 0, associated(acc_solver % matrix))
  end function wpo_assoc_matrix

  subroutine wpo_add_components() bind(c, name="wpo_add_components")
    call addcomponentstobodylists()
  end subroutine wpo_add_components

  function wpo_max_elem_dofs() result(n) bind(c, name="wpo_max_elem_dofs")
    integer(c_int) :: n
    n = acc_solver % mesh % maxelementdofs
  end function wpo_max_elem_dofs

  subroutine wpo_getmesh_touch() bind(c, name="wpo_getmesh_touch")
    acc_mesh => getmesh()
  end subroutine wpo_getmesh_touch

  function wpo_numberofbodies() result(n) bind(c, name="wpo_numberofbodies")
    integer(c_int) :: n
    n = acc_model % numberofbodies
  end function wpo_numberofbodies

  ! -- assembly driver -------------------------------------------------
  function wpo_nof_active() result(n) bind(c, name="wpo_nof_active")
    integer(c_int) :: n
    n = getnofactive()
  end function wpo_nof_active

  subroutine wpo_default_initialize() bind(c, name="wpo_default_initialize")
    call defaultinitialize()
  end subroutine wpo_default_initialize

  ! -- active-element context ------------------------------------------
  subroutine wpo_begin_active(t, n, nd, nb) bind(c, name="wpo_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n, nd, nb
    acc_elem => getactiveelement(t)
    n = getelementnofnodes()
    nd = getelementnofdofs()
    nb = getelementnofbdofs()
  end subroutine wpo_begin_active

  function wpo_body_id() result(b) bind(c, name="wpo_body_id")
    integer(c_int) :: b
    b = acc_elem % bodyid
  end function wpo_body_id

  function wpo_bodyforce() result(r) bind(c, name="wpo_bodyforce")
    integer(c_int) :: r
    acc_bf => getbodyforce()
    r = merge(1, 0, associated(acc_bf))
  end function wpo_bodyforce

  subroutine wpo_load_source(n, load, found) bind(c, name="wpo_load_source")
    integer(c_int), value :: n
    real(c_double) :: load(n)
    integer(c_int) :: found
    logical :: f
    load(1:n) = getreal(acc_bf, 'Source', f)
    found = merge(1, 0, f)
  end subroutine wpo_load_source

  ! Coil class of the current element: 0 = no component, 1 = stranded,
  ! 2 = massive, 3 = foil winding, -1 = unknown (the original Fatals).
  ! Mirrors the CompParams/CoilType query of the assembly loops.
  function wpo_comp_class() result(c) bind(c, name="wpo_comp_class")
    integer(c_int) :: c
    type(valuelist_t), pointer :: cp
    character(len=MAX_NAME_LEN) :: coiltype
    logical :: found
    c = 0
    cp => getcomponentparams(acc_elem)
    if (.not. associated(cp)) return
    coiltype = getstring(cp, 'Coil Type', found)
    if (.not. found) return
    select case (coiltype)
    case ('stranded')
      c = 1
    case ('massive')
      c = 2
    case ('foil winding')
      c = 3
    case default
      c = -1
    end select
  end function wpo_comp_class

  ! Mirrors `NoRotM = GetLogical(CompParams, 'Without RotM', Found)` of
  ! the stranded branch. Only called for class 1, as in the original.
  function wpo_comp_norotm() result(r) bind(c, name="wpo_comp_norotm")
    integer(c_int) :: r
    type(valuelist_t), pointer :: cp
    logical :: f, v
    cp => getcomponentparams(acc_elem)
    v = getlogical(cp, 'Without RotM', f)
    r = merge(1, 0, v)
  end function wpo_comp_norotm

  subroutine wpo_fatal_coil() bind(c, name="wpo_fatal_coil")
    call fatal('Wsolve', 'Non existent Coil Type Chosen!')
  end subroutine wpo_fatal_coil

  ! Mirrors the `RotM = 0` whole-array fill plus the SELECT CASE of the
  ! assembly loops (stranded conditional GetElementRotM, foil
  ! unconditional GetElementRotM). `norotm` is the host SAVE flag the
  ! Rust side carries, exactly like the original's SAVE local.
  subroutine wpo_rotm_fill(class, ld, rotm, n, norotm) &
      bind(c, name="wpo_rotm_fill")
    integer(c_int), value :: class, ld, n, norotm
    real(c_double) :: rotm(3, 3, ld)
    type(valuelist_t), pointer :: cp
    logical :: f, nonrotm
    rotm = 0.0_c_double
    select case (class)
    case (1)
      cp => getcomponentparams(acc_elem)
      nonrotm = getlogical(cp, 'Without RotM', f)
      if (.not. nonrotm) call getelementrotm(acc_elem, rotm, n)
    case (2)
    case (3)
      call getelementrotm(acc_elem, rotm, n)
    end select
  end subroutine wpo_rotm_fill

  ! Mirrors `Tcoef = 0` plus `Tcoef = GetElectricConductivityTensor(...)`.
  subroutine wpo_tcoef_fill(class, ld, tcoef, n) bind(c, name="wpo_tcoef_fill")
    integer(c_int), value :: class, ld, n
    real(c_double) :: tcoef(3, 3, ld)
    character(len=MAX_NAME_LEN) :: coiltype
    logical :: coilbody
    coilbody = class /= 0
    select case (class)
    case (1)
      coiltype = 'stranded'
    case (2)
      coiltype = 'massive'
    case (3)
      coiltype = 'foil winding'
    case default
      coiltype = ''
    end select
    tcoef = 0.0_c_double
    tcoef = getelectricconductivitytensor(acc_elem, n, 're', coilbody, coiltype)
  end subroutine wpo_tcoef_fill

  ! -- LocalMatrix (whole-element kernel) ------------------------------
  ! `norotm` is the host SAVE flag carried by the Rust side. `stiff` and
  ! `force` arrive pre-zeroed over the whole (ld) arrays, exactly as the
  ! original's `STIFF = 0` / `FORCE = 0` leaves them before the update
  ! statements below; only the (1:nd, 1:nd) / (1:nd) sections are
  ! touched here.
  subroutine wpo_localmatrix(ld, stiff, force, load, tcoef, rotm, n, nd, &
      class, norotm) bind(c, name="wpo_localmatrix")
    integer(c_int), value :: ld, n, nd, class, norotm
    real(c_double) :: stiff(ld, ld), force(ld), load(ld)
    real(c_double) :: tcoef(3, 3, ld), rotm(3, 3, ld)
    real(c_double) :: basis(nd), dbasisdx(nd, 3), detj, loadatip
    real(c_double) :: c(3, 3), rotmloc(3, 3)
    character(len=MAX_NAME_LEN) :: coiltype
    logical :: stat, coilbody, nonrotm
    integer :: i, j, t
    type(gaussintegrationpoints_t) :: ip

    call getelementnodes(acc_nodes_lm)

    c(1, 1:3) = 0.0_c_double
    c(2, 1:3) = 0.0_c_double
    c(3, 1:3) = 0.0_c_double

    coilbody = class /= 0
    select case (class)
    case (1)
      coiltype = 'stranded'
    case (2)
      coiltype = 'massive'
    case (3)
      coiltype = 'foil winding'
    case default
      coiltype = ''
    end select
    nonrotm = norotm /= 0

    ip = gausspoints(acc_elem)
    do t = 1, ip % n
      stat = elementinfo(acc_elem, acc_nodes_lm, ip % u(t), ip % v(t), &
        ip % w(t), detj, basis, dbasisdx)

      do i = 1, 3
        do j = 1, 3
          c(i, j) = sum(tcoef(i, j, 1:n) * basis(1:n))
          if (coilbody .and. coiltype /= 'massive' .and. .not. nonrotm) then
            rotmloc(i, j) = sum(rotm(i, j, 1:n) * basis(1:n))
          end if
        end do
      end do

      if (coilbody .and. coiltype /= 'massive' .and. .not. nonrotm) then
        c = matmul(matmul(rotmloc, c), transpose(rotmloc))
      end if

      loadatip = sum(basis(1:n) * load(1:n))

      do i = 1, nd
        do j = 1, nd
          stiff(i, j) = stiff(i, j) + &
            sum(matmul(c, dbasisdx(i, :)) * dbasisdx(j, :)) * detj * ip % s(t)
        end do
      end do

      force(1:nd) = force(1:nd) + ip % s(t) * detj * loadatip * basis(1:nd)
    end do
  end subroutine wpo_localmatrix

  ! -- shared assembly tail ---------------------------------------------
  subroutine wpo_condensate(nd, nb, ld, k, f) bind(c, name="wpo_condensate")
    integer(c_int), value :: nd, nb, ld
    real(c_double) :: k(ld, ld), f(ld)
    call condensatep(nd, nb, k(1:ld, 1:ld), f(1:ld))
  end subroutine wpo_condensate

  subroutine wpo_default_update(ld, g, f) bind(c, name="wpo_default_update")
    integer(c_int), value :: ld
    real(c_double) :: g(ld, ld), f(ld)
    call defaultupdateequations(g(1:ld, 1:ld), f(1:ld))
  end subroutine wpo_default_update

  ! -- boundary-element context ------------------------------------------
  function wpo_num_bnd_elems() result(n) bind(c, name="wpo_num_bnd_elems")
    integer(c_int) :: n
    n = acc_solver % mesh % numberofboundaryelements
  end function wpo_num_bnd_elems

  subroutine wpo_begin_boundary(t, n, active, family) &
      bind(c, name="wpo_begin_boundary")
    integer(c_int), value :: t
    integer(c_int) :: n, active, family
    acc_elem => getboundaryelement(t)
    active = merge(1, 0, activeboundaryelement())
    n = getelementnofnodes()
    family = getelementfamily()
  end subroutine wpo_begin_boundary

  subroutine wpo_bc() bind(c, name="wpo_bc")
    acc_bc => getbc()
  end subroutine wpo_bc

  subroutine wpo_bc_electrode(which, val, found) &
      bind(c, name="wpo_bc_electrode")
    integer(c_int), value :: which
    integer(c_int) :: val, found
    logical :: l, f
    if (which == 1) then
      l = getlogical(acc_bc, 'Positive Electrode', f)
    else
      l = getlogical(acc_bc, 'Negative Electrode', f)
    end if
    val = merge(1, 0, l)
    found = merge(1, 0, f)
  end subroutine wpo_bc_electrode

  subroutine wpo_default_finish() bind(c, name="wpo_default_finish")
    call defaultfinishassembly()
  end subroutine wpo_default_finish

  subroutine wpo_default_dirichlet() bind(c, name="wpo_default_dirichlet")
    call defaultdirichletbcs()
  end subroutine wpo_default_dirichlet

  function wpo_default_solve() result(norm) bind(c, name="wpo_default_solve")
    real(c_double) :: norm
    norm = defaultsolve()
  end function wpo_default_solve

  ! -- GetWNormsForBodies (whole-call kernel) ---------------------------
  ! Fills `wnorms(nofbodies)`; mirrors the function body including the
  ! ParallelReduction loop and the Info call. The AddElementWNormAndVolume
  ! internal procedure is a module subroutine below with its own SAVE
  ! Nodes mirror; `nd` is wrapper-local (host-associated in the
  ! original, but always rewritten per element on both sides).
  subroutine wpo_wnorms(nofbodies, wnorms) bind(c, name="wpo_wnorms")
    integer(c_int), value :: nofbodies
    real(c_double) :: wnorms(nofbodies)
    real(c_double) :: volumes(nofbodies)
    real(c_double) :: wnormcoeffs(nofbodies)
    integer :: active, t, n, nd, i
    logical :: coilbody, found
    character(len=MAX_NAME_LEN) :: coiltype
    type(valuelist_t), pointer :: compparams

    wnorms = 0.0_c_double
    volumes = 0.0_c_double
    wnormcoeffs = 0.0_c_double
    active = getnofactive()
    do t = 1, active
      acc_elem => getactiveelement(t)
      n = getelementnofnodes()
      nd = getelementnofdofs()

      coilbody = .false.
      compparams => getcomponentparams(acc_elem)
      coiltype = ''
      if (associated(compparams)) then
        coiltype = getstring(compparams, 'Coil Type', found)
        if (found) coilbody = .true.
      end if

      if (coilbody) then
        call wpo_addwnorm(acc_elem, n, nd, &
            wnormcoeffs(acc_elem % bodyid), &
            volumes(acc_elem % bodyid))
      end if
    end do

    do i = 1, nofbodies
      wnormcoeffs(i) = parallelreduction(wnormcoeffs(i))
      volumes(i) = parallelreduction(volumes(i))
      if (volumes(i) /= 0.0_dp) then
        wnorms(i) = wnormcoeffs(i) / volumes(i)
      else
        wnorms(i) = 0.0_c_double
      end if
    end do

    call info('GetWnormsForBodies', &
         'W norm computed in components.', level=9)
  end subroutine wpo_wnorms

  ! Mirrors AddElementWNormAndVolume (called only from above).
  subroutine wpo_addwnorm(element, n, nd, wnormcoeff, volume)
    type(element_t), target :: element
    integer :: n, nd
    real(c_double) :: wnormcoeff, volume
    integer :: active, j
    logical :: stat
    type(gaussintegrationpoints_t) :: ip
    real(c_double) :: basis(nd), dbasisdx(nd, 3), detj
    real(c_double) :: wbase(n), s, w(3)

    call getelementnodes(acc_nodes_wn)
    call getlocalsolution(wbase, 'W')

    ip = gausspoints(element)

    do j = 1, ip % n
      stat = elementinfo(element, acc_nodes_wn, ip % u(j), ip % v(j), &
          ip % w(j), detj, basis, dbasisdx)

      s = ip % s(j) * detj
      volume = volume + s

      w = -matmul(wbase(1:n), dbasisdx(1:n, :))
      wnormcoeff = wnormcoeff + sqrt(sum(w**2)) * s
    end do
  end subroutine wpo_addwnorm

  ! -- SaveElementWSolution (whole-element kernel) -----------------------
  ! Mirrors the internal procedure body; `class`/`norotm` replace the
  ! host-associated CoilBody/CoilType/NoRotM. SAVE state (variable
  ! handles, First flag, names, Nodes) lives here with the same
  ! lifetimes as the original.
  subroutine wpo_save_element(ld, tcoef, rotm, n, wnorm, class, norotm) &
      bind(c, name="wpo_save_element")
    integer(c_int), value :: ld, n, class, norotm
    real(c_double) :: tcoef(3, 3, ld), rotm(3, 3, ld)
    real(c_double), value :: wnorm
    integer :: j, k, t, i
    type(valuelist_t), pointer :: solverparams
    logical :: found, computejvec, nonrotm, coilbody
    character(len=MAX_NAME_LEN) :: coiltype
    real(c_double) :: basis(n), detj
    real(c_double) :: dbasisdx(n, 3)
    real(c_double) :: wpot(n), wvec(3), jvec(3)
    type(gaussintegrationpoints_t) :: ip
    logical :: stat
    real(c_double) :: c(3, 3), rotmloc(3, 3)
    ! Local copy of the Wnorm dummy: the original reassigns its dummy
    ! (`IF(CoilType/='stranded') Wnorm = 1`) inside the k-loops below.
    real(c_double) :: wn

    call getelementnodes(acc_nodes_sv)
    acc_varname = getstring(getsolverparams(), 'Variable', found)
    if (.not. found) acc_varname = 'W'
    computejvec = getlogical(getsolverparams(), 'Compute J Vector', found)
    if (.not. found) computejvec = .false.

    if (computejvec) then
      acc_jvarname = getstring(getsolverparams(), 'J Variable', found)
      if (.not. found) acc_jvarname = 'J'
    end if

    if (acc_sv_first) then
      acc_sv_first = .false.
      acc_wpotvar => variableget(acc_mesh % variables, &
        trim(acc_varname)//' Potential')
      if (.not. associated(acc_wpotvar)) then
        call fatal('SaveElementWSolution()', 'W Potential variable not found')
      end if
      acc_wvecvar => variableget(acc_mesh % variables, &
        trim(acc_varname)//' Vector E')
      if (.not. associated(acc_wvecvar)) then
        call fatal('SaveElementWSolution()', 'W Vector variable not found')
      end if

      if (computejvec) then
        acc_jvecvar => variableget(acc_mesh % variables, &
          trim(acc_jvarname)//' Vector E')
        if (.not. associated(acc_jvecvar)) then
          call fatal('SaveElementWSolution()', 'J Vector variable not found')
        end if
      end if
    end if

    coilbody = class /= 0
    select case (class)
    case (1)
      coiltype = 'stranded'
    case (2)
      coiltype = 'massive'
    case (3)
      coiltype = 'foil winding'
    case default
      coiltype = ''
    end select
    nonrotm = norotm /= 0

    wn = wnorm
    call getlocalsolution(wpot, acc_varname)
    ip = gausspoints(acc_elem)
    do t = 1, n
      if (associated(acc_wpotvar)) then
        do k = 1, acc_wpotvar % dofs
          if (coiltype /= 'stranded') wn = 1.0_dp
          if (wn > epsilon(wn)) then
            acc_wpotvar % values(acc_wpotvar % dofs * (acc_wpotvar % perm( &
                acc_elem % dgindexes(t)) - 1) + k) = wpot(t) / wn
          end if
        end do
      end if

      if (associated(acc_wvecvar)) then
        stat = elementinfo(acc_elem, acc_nodes_sv, ip % u(t), ip % v(t), &
                  ip % w(t), detj, basis, dbasisdx)

        wvec = matmul(wpot(1:n), dbasisdx(1:n, :))

        do k = 1, acc_wvecvar % dofs
          if (wn > epsilon(wn)) then
            if (coiltype /= 'stranded') wn = 1.0_dp
            acc_wvecvar % values(acc_wvecvar % dofs * (acc_wvecvar % perm( &
                acc_elem % dgindexes(t)) - 1) + k) = wvec(k) / wn
          end if
        end do

        if (computejvec) then
          if (associated(acc_jvecvar)) then
            do i = 1, 3
              do j = 1, 3
                c(i, j) = sum(tcoef(i, j, 1:n) * basis(1:n))
                if (coilbody .and. coiltype /= 'massive' .and. &
                    .not. nonrotm) then
                  rotmloc(i, j) = sum(rotm(i, j, 1:n) * basis(1:n))
                end if
              end do
            end do

            if (coilbody .and. coiltype /= 'massive' .and. &
                .not. nonrotm) then
              c = matmul(matmul(rotmloc, c), transpose(rotmloc))
            end if
            jvec = matmul(c, wvec)

            do k = 1, acc_wvecvar % dofs
              if (wn > epsilon(wn)) then
                acc_jvecvar % values(acc_jvecvar % dofs * &
                  (acc_jvecvar % perm( &
                  acc_elem % dgindexes(t)) - 1) + k) = jvec(k)
              end if
            end do
          end if
        end if
      end if
    end do
  end subroutine wpo_save_element

  subroutine wpo_fatal_alloc() bind(c, name="wpo_fatal_alloc")
    call fatal('Wsolve', 'Memory allocation error.')
  end subroutine wpo_fatal_alloc

end module wpo_acc
