! NodeToEdgeSolver port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/NodeToEdgeSolver.F90, compiled against the
! built .mod files of the target configure (the Poisson-port mechanism).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results,
! WRITE/PRINT I/O, MAXVAL/COUNT reductions, whole-array updates —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module mirrors the plugin's locals (current solver/model, meshes,
! variables, the LocalMatrix SAVE Nodes and IP) with SAVE lifetimes
! matching the original's automatics/SAVE, reset or refilled at the
! same points as the original.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module nodetoedgesolver_acc
  use iso_c_binding
  use defutils
  ! Keyword calls (`GaussPoints`, `ElementInfo`, `EdgeElementStyle`) need
  ! the defining modules' explicit interfaces (`USE DefUtils` alone does
  ! not re-export them for keyword use). Same procedures, same
  ! interfaces — no behavior change.
  use integration
  use eleminfo
  implicit none

  character(*), parameter :: caller_er = 'ExtrudedRestart'
  character(*), parameter :: caller_nte = 'NodeToEdgeField'

  ! -- shared entry handles -------------------------------------------
  type(solver_t), pointer, save :: acc_solver => null()
  type(model_t), pointer, save :: acc_model => null()

  ! -- ExtrudedRestart per-invocation locals (original automatics) ----
  type(mesh_t), pointer, save :: acc_thismesh => null()
  type(mesh_t), pointer, save :: acc_targetmesh => null()
  type(solver_t), pointer, save :: acc_psolver => null()
  type(variable_t), pointer, save :: acc_var => null()
  type(variable_t), pointer, save :: acc_pvar => null()
  integer, pointer, save :: acc_pperm(:) => null()
  real(kind=dp), pointer, save :: acc_pvals(:) => null()
  character(len=max_name_len), save :: acc_name
  character(len=max_name_len), save :: acc_varname
  character(len=max_name_len), save :: acc_targetname
  integer, save :: acc_solverind = 0
  integer, save :: acc_novar = 0
  integer, save :: acc_layers = 0
  integer, save :: acc_m = 0

  ! -- NodeToEdgeField per-invocation locals (original automatics) ---
  type(mesh_t), pointer, save :: acc_mesh => null()
  type(variable_t), pointer, save :: acc_nodalvar => null()
  type(variable_t), pointer, save :: acc_edgevar => null()
  type(variable_t), pointer, save :: acc_thisvar => null()
  type(element_t), pointer, save :: acc_elem => null()
  integer, pointer, save :: acc_nodeidx(:) => null()
  logical, save :: acc_constbulk = .false.
  logical, save :: acc_bulkassoc = .false.
  logical, save :: acc_ready = .false.
  logical, save :: acc_piola = .false.
  logical, save :: acc_isim = .false.
  logical, save :: acc_iscomplex = .false.
  integer, save :: acc_degree = 1
  integer, save :: acc_dofs = 0
  integer, save :: acc_dim = 0
  integer, save :: acc_activecomp = 0
  integer, save :: acc_imoffset = 0
  ! LocalMatrix context (mirrors its SAVE Nodes; IP is refilled per
  ! element exactly as the original's local is reassigned per call).
  type(nodes_t), save :: acc_nodes
  type(gaussintegrationpoints_t), save :: acc_ip

contains

  ! -- ExtrudedRestart_Init0 entry ------------------------------------
  subroutine nte_xinit0() bind(c, name="nte_xinit0")
    type(valuelist_t), pointer :: params
    params => getsolverparams()
    call listaddnewlogical(params, 'Mesh Enforce Local Copy', .true.)
    call listaddnewlogical(params, 'Mesh Output', .false.)
    call listaddnewinteger(params, 'Active Coordinate', 3)
    call listaddnewinteger(params, 'Interpolation Passive Coordinate', 3)
  end subroutine nte_xinit0

  ! -- ExtrudedRestart entry state ------------------------------------
  subroutine nte_xbegin(h_model, h_solver) bind(c, name="nte_xbegin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine nte_xbegin

  ! -- mesh resolution + variable count (original lines 94-146) -------
  subroutine nte_xstart(novar) bind(c, name="nte_xstart")
    integer(c_int) :: novar
    type(valuelist_t), pointer :: params
    integer :: i
    logical :: found
    call info(caller_er, 'Mapping result between meshes')
    acc_thismesh => getmesh()
    call info(caller_er, 'This mesh name is: '//trim(acc_thismesh % name), level=20)
    call info(caller_er, 'This mesh dimension is: '//i2s(acc_thismesh % meshdim), level=20)
    params => getsolverparams()
    acc_targetmesh => null()
    acc_solverind = listgetinteger(params, 'Target Mesh Solver Index', found)
    if (found) then
      ! Target mesh solver is explicitly given
      acc_targetmesh => currentmodel % solvers(acc_solverind) % mesh
      if (.not. associated(acc_targetmesh)) then
        call fatal(caller_er, 'No mesh associated to Solver: '//i2s(acc_solverind))
      end if
    else
      ! Otherwise use the 1st mesh that is not this old data mesh
      do i = 1, currentmodel % numberofsolvers
        acc_targetmesh => currentmodel % solvers(i) % mesh
        if (.not. associated(acc_thismesh, acc_targetmesh)) then
          acc_solverind = i
          exit
        end if
      end do
      if (acc_solverind == 0) then
        call fatal(caller_er, 'Could not find target mesh that is not this mesh!')
      end if
    end if
    call info(caller_er, 'Using target mesh from Solver index: '//i2s(acc_solverind), level=8)
    call info(caller_er, 'Target mesh name is: '//trim(acc_targetmesh % name), level=8)
    call info(caller_er, 'Target mesh dimension is: '//i2s(acc_targetmesh % meshdim), level=8)
    acc_psolver => acc_model % solvers(acc_solverind)
    acc_novar = 0
    do i = 1, 100
      write(acc_name, '(A,I0)') 'Variable ', i
      if (.not. listcheckpresent(params, acc_name)) exit
      acc_novar = i
    end do
    call info(caller_er, 'Number of variables to be mapped: '//i2s(acc_novar), level=8)
    acc_m = acc_targetmesh % numberofnodes
    if (associated(acc_psolver % variable)) then
      acc_m = max(acc_m, size(acc_psolver % variable % perm))
    end if
    acc_layers = acc_targetmesh % numberofnodes / acc_thismesh % numberofnodes
    call info(caller_er, 'Number of layers to be mapped: '//i2s(acc_layers), level=8)
    novar = acc_novar
  end subroutine nte_xstart

  ! -- one mapped variable (original lines 148-244) -------------------
  subroutine nte_xmap(i) bind(c, name="nte_xmap")
    integer(c_int), value :: i
    type(valuelist_t), pointer :: params
    integer :: j, k, l, n, dofs, maxperm
    logical :: found, createvar
    params => getsolverparams()
    write(acc_name, '(A,I0)') 'Variable ', i
    acc_varname = getstring(params, acc_name, found)
    acc_var => variableget(acc_thismesh % variables, acc_varname, thisonly=.true.)
    if (.not. associated(acc_var)) then
      if (infoactive(20)) then
        print *, 'List of variable in 2D mesh:'
        acc_var => acc_thismesh % variables
        do while (associated(acc_var))
          print *, trim(acc_var % name), size(acc_var % values), associated(acc_var % perm)
          acc_var => acc_var % next
        end do
      end if
      call fatal(caller_er, 'Could not find variable: '//trim(acc_varname))
    end if
    dofs = acc_var % dofs
    if (infoactive(20)) then
      call vectorvaluesrange(acc_var % values, size(acc_var % values), trim(acc_var % name))
    end if
    write(acc_name, '(A,I0)') 'Extruded Variable ', i
    acc_targetname = getstring(params, acc_name, found)
    if (.not. found) acc_targetname = acc_varname
    acc_pvar => variableget(acc_targetmesh % variables, acc_targetname, thisonly=.true.)
    createvar = .not. associated(acc_pvar)
    if (.not. associated(acc_var % perm)) then
      ! One intended use of this module is to extrude data from 2D electrical machine computation
      ! to 3D one. The it is often desirable also to copy the related electrical circuits that may be
      ! found in the Lagrange multiplier values not associated to any permutation. So these are
      ! copied as one-to-one from 2D mesh to 3D mesh.
      !------------------------------------------------------------------------------------------------
      n = size(acc_var % values)
      if (createvar) then
        nullify(acc_pvals)
        allocate(acc_pvals(n))
        acc_pvals = 0.0_dp
        call variableaddvector(acc_targetmesh % variables, acc_targetmesh, &
            acc_model % solvers(acc_solverind), acc_targetname, acc_var % dofs, acc_pvals)
      else
        acc_pvals => acc_pvar % values
      end if
      acc_pvals(1:n) = acc_var % values
      call info(caller_er, 'Copied variable as such from 2D mesh to 3D mesh: '//trim(acc_varname), level=8)
    else
      maxperm = maxval(acc_var % perm)
      nullify(acc_pperm, acc_pvals)
      if (createvar) then
        allocate(acc_pperm(acc_m))
        acc_pperm = 0
        allocate(acc_pvals(acc_layers*size(acc_var % values)))
      else
        acc_pperm => acc_pvar % perm
        acc_pvals => acc_pvar % values
      end if
      acc_pvals = 0.0_dp
      ! Here we assume that the fields to be mapped are nodal ones and the mesh is linear one!
      n = acc_thismesh % numberofnodes
      do j = 0, acc_layers-1
        do k = 1, acc_thismesh % numberofnodes
          if (acc_var % perm(k) == 0) cycle
          if (createvar) then
            acc_pperm(k+j*n) = acc_var % perm(k) + j * maxperm
          end if
          do l = 1, dofs
            acc_pvals(dofs*(acc_pperm(k+j*n)-1)+l) = &
                acc_var % values(dofs*(acc_var % perm(k)-1)+l)
          end do
        end do
      end do
      if (createvar) then
        call variableaddvector(acc_targetmesh % variables, acc_targetmesh, &
            acc_model % solvers(acc_solverind), &
            acc_targetname, acc_var % dofs, acc_pvals, acc_pperm, vartype=acc_var % type)
        acc_pvar => variableget(acc_targetmesh % variables, acc_varname, thisonly=.true.)
        if (associated(acc_pvar)) then
          call info(caller_er, 'Created variable: '//trim(acc_varname), level=20)
        else
          call fatal(caller_er, 'Failed to create variable: '//trim(acc_varname))
        end if
      end if
      if (.not. acc_var % periodicflipactive) then
        if (count(acc_var % perm > 0) > size(acc_var % values) / acc_var % dofs) then
          call info(caller_er, 'The variable that we read in seems to have been conforming!')
          acc_var % periodicflipactive = .true.
        end if
      end if
      ! Inherit the periodic flips, if any
      acc_pvar % periodicflipactive = acc_var % periodicflipactive
      call info(caller_er, 'Extruded variable from 2D mesh to 3D mesh: '//trim(acc_varname), level=8)
    end if
    if (infoactive(20)) then
      call vectorvaluesrange(acc_pvals, size(acc_pvals), trim(acc_targetname))
    end if
  end subroutine nte_xmap

  ! -- final tally (original line 252) --------------------------------
  subroutine nte_xfinish() bind(c, name="nte_xfinish")
    call info(caller_er, 'Transferred '//i2s(acc_novar)//' variables from 2d to 3d mesh!', level=7)
  end subroutine nte_xfinish

  ! -- NodeToEdgeField_Init0 entry (original lines 260-299) -----------
  subroutine nte_ninit0() bind(c, name="nte_ninit0")
    type(valuelist_t), pointer :: params
    logical :: found, pivolaversion, secondorder, secondkind
    params => getsolverparams()
    call listaddlogical(params, 'Linear System Refactorize', .false.)
    if (.not. listcheckpresent(params, 'Element')) then
      ! We use one place where all the edge element keywords are defined and checked.
      call edgeelementstyle(params, pivolaversion, secondkind, secondorder, check=.true.)
      if (secondorder) then
        call listaddstring(params, 'Element', &
            'n:0 e:2 -brick b:6 -pyramid b:3 -prism b:2 -quad_face b:4 -tri_face b:2')
      else if (secondkind) then
        call listaddstring(params, 'Element', &
            'n:0 e:2')
      else if (pivolaversion) then
        call listaddstring(params, 'Element', &
            'n:0 e:1 -brick b:3 -quad_face b:2')
      else
        call listaddstring(params, 'Element', 'n:0 e:1')
      end if
    end if
    call listaddnewstring(params, 'Variable', '-nooutput edgefield proj')
    call listaddinteger(params, 'Time derivative Order', 0)
    call listaddinteger(params, 'Active Coordinate', 3)
    call listaddlogical(params, 'Edge Basis', .true.)
  end subroutine nte_ninit0

  ! -- NodeToEdgeField entry state (original lines 343-345) -----------
  subroutine nte_nbegin(h_solver) bind(c, name="nte_nbegin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
    call info(caller_nte, 'Projecting 3D nodal field to Hcurl field')
    call defaultstart()
  end subroutine nte_nbegin

  ! -- resolution (original lines 347-389) ----------------------------
  subroutine nte_nsetup(dofs, dim, iscomplex) bind(c, name="nte_nsetup")
    integer(c_int) :: dofs, dim, iscomplex
    type(valuelist_t), pointer :: params
    logical :: found
    acc_dim = coordinatesystemdimension()
    params => getsolverparams()
    acc_mesh => getmesh()
    acc_psolver => acc_solver
    ! Check if accelerated assembly is desired:
    acc_constbulk = getlogical(params, 'Constant Bulk Matrix', found)
    acc_bulkassoc = associated(acc_solver % matrix % bulkvalues)
    acc_ready = acc_bulkassoc .and. acc_constbulk
    if (acc_ready) call info(caller_nte, &
        'Assuming that the system matrix is already created!', level=10)
    acc_thisvar => acc_solver % variable
    if (.not. associated(acc_thisvar)) then
      call fatal(caller_nte, 'No variable associated to solver!')
    end if
    if (acc_thisvar % dofs /= 1) call fatal(caller_nte, 'A real-valued potential expected')
    ! Find the variable which defines the nodal field to be mapped.
    ! If the nodal field is a scalar field we need to know
    ! to which component it relates to!
    !------------------------------------------------------
    acc_name = getstring(params, 'Nodal Variable', found)
    if (.not. found) acc_name = 'az'
    acc_nodalvar => variableget(acc_mesh % variables, acc_name)
    if (.not. associated(acc_nodalvar)) then
      call fatal(caller_nte, 'Nodal variable not associated: '//trim(acc_name))
    end if
    call info(caller_nte, 'Using nodal variable for projection: '//trim(acc_name), level=10)
    if (infoactive(20)) then
      call vectorvaluesrange(acc_nodalvar % values, &
          size(acc_nodalvar % values), trim(acc_nodalvar % name))
    end if
    acc_dofs = acc_nodalvar % dofs
    acc_iscomplex = (acc_dofs == 6)
    acc_isim = .false.
    acc_imoffset = 0
    if (acc_dofs < acc_dim) then
      acc_activecomp = listgetinteger(params, 'Active Coordinate')
      call info(caller_nte, 'Setting active coordinate to be: '//i2s(acc_activecomp), level=10)
    else
      acc_activecomp = 0
    end if
    dofs = acc_dofs
    dim = acc_dim
    iscomplex = merge(1, 0, acc_iscomplex)
  end subroutine nte_nsetup

  ! -- edge variable lookup (original lines 391-410) ------------------
  subroutine nte_nedge() bind(c, name="nte_nedge")
    type(valuelist_t), pointer :: params
    logical :: found
    params => getsolverparams()
    ! Find the variable which is projected:
    !---------------------------------------------------------
    acc_name = getstring(params, 'Edge Variable', found)
    if (found) then
      acc_edgevar => variableget(acc_mesh % variables, acc_name)
      if (.not. associated(acc_edgevar)) then
        call fatal(caller_nte, 'Could not find target variable for projection!')
      end if
    else
      acc_name = 'av'
      acc_edgevar => variableget(acc_mesh % variables, acc_name)
      if (.not. associated(acc_edgevar)) then
        call warn(caller_nte, 'Could not find target variable for projection!')
      end if
    end if
    if (infoactive(20)) then
      call vectorvaluesrange(acc_nodalvar % values, &
          size(acc_nodalvar % values), trim(acc_nodalvar % name))
    end if
  end subroutine nte_nedge

  ! -- edge element style (original lines 414-415) --------------------
  subroutine nte_nstyle() bind(c, name="nte_nstyle")
    type(valuelist_t), pointer :: params
    params => getsolverparams()
    ! These should be consistent with the primary solver!!
    !-------------------------------------------------------------------------------------
    call edgeelementstyle(params, acc_piola, basisdegree=acc_degree)
    if (acc_piola) call info(caller_nte, 'Using Piola-transformed finite elements', level=5)
  end subroutine nte_nstyle

  ! -- first-mesh sizing (original line 420) --------------------------
  function nte_nmaxdofs() result(n) bind(c, name="nte_nmaxdofs")
    integer(c_int) :: n
    n = acc_mesh % maxelementdofs
  end function nte_nmaxdofs

  ! -- pass state: label-1 target + imaginary redo (lines 423, 519-525)
  subroutine nte_npass(second) bind(c, name="nte_npass")
    integer(c_int), value :: second
    if (second == 0) then
      acc_isim = .false.
      acc_imoffset = 0
      acc_ready = acc_bulkassoc .and. acc_constbulk
    else
      acc_isim = .true.
      acc_imoffset = acc_dofs / 2
      acc_ready = acc_constbulk
    end if
  end subroutine nte_npass

  subroutine nte_ninit() bind(c, name="nte_ninit")
    call defaultinitialize(acc_solver, acc_ready)
  end subroutine nte_ninit

  function nte_nactive() result(n) bind(c, name="nte_nactive")
    integer(c_int) :: n
    n = getnofactive()
  end function nte_nactive

  ! -- active-element context (original lines 427-431) -----------------
  subroutine nte_nelem(t, n, nd) bind(c, name="nte_nelem")
    integer(c_int), value :: t
    integer(c_int) :: n, nd
    acc_elem => getactiveelement(t)
    n = getelementnofnodes()
    nd = getelementnofdofs()
    acc_nodeidx => acc_elem % nodeindexes
  end subroutine nte_nelem

  ! -- Anodal fill (original lines 433-444) ---------------------------
  subroutine nte_nanodal(ld, n, anodal) bind(c, name="nte_nanodal")
    integer(c_int), value :: ld, n
    real(c_double) :: anodal(3, ld)
    integer :: i
    if (acc_dofs == 1) then
      anodal = 0.0_dp
      where (acc_nodalvar % perm(acc_nodeidx(1:n)) > 0)
        anodal(acc_activecomp, 1:n) = &
            acc_nodalvar % values(acc_nodalvar % perm(acc_nodeidx(1:n)))
      end where
    else
      do i = 1, acc_dim
        where (acc_nodalvar % perm(acc_nodeidx(1:n)) > 0)
          anodal(i, 1:n) = acc_nodalvar % values( &
              acc_dofs*(acc_nodalvar % perm(acc_nodeidx(1:n))-1)+i+acc_imoffset)
        end where
      end do
    end if
  end subroutine nte_nanodal

  ! -- periodic flip (original lines 446-453) --------------------------
  subroutine nte_naflip(ld, n, anodal) bind(c, name="nte_naflip")
    integer(c_int), value :: ld, n
    real(c_double) :: anodal(3, ld)
    integer :: i, j
    if (acc_nodalvar % periodicflipactive) then
      do i = 1, acc_dofs
        if (acc_dofs == 1) j = acc_activecomp
        where (acc_mesh % periodicflip(acc_nodeidx(1:n)))
          anodal(j, 1:n) = -anodal(j, 1:n)
        end where
      end do
    end if
  end subroutine nte_naflip

  ! -- LocalMatrix head (original lines 565-574 + zeroing 567-568) ----
  function nte_nip(ld, stiff, force) result(nip) bind(c, name="nte_nip")
    integer(c_int), value :: ld
    real(c_double) :: stiff(ld, ld), force(ld)
    integer(c_int) :: nip
    call getelementnodes(acc_nodes)
    stiff = 0.0d0
    force = 0.0d0
    acc_ip = gausspoints(acc_elem, edgebasis=.true., &
        preferenceelement=acc_piola, edgebasisdegree=acc_degree)
    if (acc_dim == 2 .and. .not. acc_piola) then
      call fatal(caller_nte, '"Use Piola Transform = True" needed in 2D')
    end if
    nip = acc_ip % n
  end function nte_nip

  ! -- LocalMatrix integration point (original lines 576-597) ---------
  subroutine nte_npoint(t, ld, n, nd, anodal, stiff, force) bind(c, name="nte_npoint")
    integer(c_int), value :: t, ld, n, nd
    real(c_double) :: anodal(3, ld), stiff(ld, ld), force(ld)
    real(kind=dp) :: basis(n), dbasis(n, 3), wbasis(nd, 3), curlwbasis(nd, 3), aip(3)
    real(kind=dp) :: s, detj
    integer :: i, p, q
    logical :: stat
    stat = elementinfo(acc_elem, acc_nodes, acc_ip % u(t), acc_ip % v(t), &
        acc_ip % w(t), detj, basis, dbasis, edgebasis=wbasis, &
        rotbasis=curlwbasis, usolver=acc_psolver)
    s = detj * acc_ip % s(t)
    aip = 0.0d0
    do i = 1, acc_dim
      aip(i) = sum(anodal(i, 1:n) * basis(1:n))
    end do
    if (.not. acc_ready) then
      do p = 1, nd
        do q = 1, nd
          stiff(p, q) = stiff(p, q) + sum(wbasis(q, 1:acc_dim) * wbasis(p, 1:acc_dim)) * s
        end do
      end do
    end if
    do p = 1, nd
      force(p) = force(p) + sum(aip(1:acc_dim) * wbasis(p, 1:acc_dim)) * s
    end do
  end subroutine nte_npoint

  ! -- global update (original lines 462-466) --------------------------
  subroutine nte_nupdate(ld, stiff, force) bind(c, name="nte_nupdate")
    integer(c_int), value :: ld
    real(c_double) :: stiff(ld, ld), force(ld)
    ! Update global matrix and rhs vector from local matrix & vector:
    !---------------------------------------------------------------
    if (acc_ready) then
      call defaultupdateforce(force(1:ld))
    else
      call defaultupdateequations(stiff(1:ld, 1:ld), force(1:ld))
    end if
  end subroutine nte_nupdate

  ! -- bulk finish (original lines 469-483) ----------------------------
  subroutine nte_nbulk() bind(c, name="nte_nbulk")
    if (acc_constbulk) then
      call defaultfinishbulkassembly(bulkupdate=.not. acc_ready, rhsupdate=.false.)
    else
      call defaultfinishbulkassembly()
    end if
    call defaultfinishboundaryassembly()
    if (infoactive(20)) then
      call vectorvaluesrange(acc_solver % matrix % values, &
          size(acc_solver % matrix % values), 'A21')
    end if
    call defaultfinishassembly()
    if (infoactive(20)) then
      call vectorvaluesrange(acc_solver % matrix % values, &
          size(acc_solver % matrix % values), 'A22')
    end if
  end subroutine nte_nbulk

  ! -- dirichlet + solve (original lines 486-488) ----------------------
  function nte_nsolve() result(norm) bind(c, name="nte_nsolve")
    real(c_double) :: norm
    call defaultdirichletbcs()
    norm = defaultsolve()
  end function nte_nsolve

  ! -- pass tail (original lines 490-516) ------------------------------
  subroutine nte_ntail() bind(c, name="nte_ntail")
    integer :: i, j, k, n
    call defaultfinish()
    if (infoactive(20)) then
      call vectorvaluesrange(acc_thisvar % values, &
          size(acc_thisvar % values), trim(acc_thisvar % name))
    end if
    ! Finally, redefine the potential variable:
    !---------------------------------------------------------
    if (associated(acc_edgevar)) then
      n = size(acc_solver % variable % perm)
      if (n /= size(acc_edgevar % perm)) then
        call fatal(caller_nte, 'The variable and potential permutations differ')
      end if
      do i = 1, n
        j = acc_solver % variable % perm(i)
        if (j < 1) cycle
        k = acc_edgevar % perm(i)
        if (k < 1) call fatal(caller_nte, &
            'The variable and potential permutations are non-matching?')
        if (acc_iscomplex) then
          k = 2*k
          if (.not. acc_isim) k = k-1
        end if
        acc_edgevar % values(k) = acc_solver % variable % values(j)
      end do
    end if
  end subroutine nte_ntail

  ! -- done (original lines 527-537) -----------------------------------
  subroutine nte_ndone() bind(c, name="nte_ndone")
    if (associated(acc_edgevar)) then
      if (infoactive(20)) then
        call vectorvaluesrange(acc_edgevar % values, &
            size(acc_edgevar % values), trim(acc_edgevar % name))
      end if
    end if
    ! We should be visiting this routine only once!
    !CALL Info(Caller,'Freeing unneeded matrix structures',Level=10)
    !CALL FreeMatrix(Solver % Matrix)
    call info(caller_nte, 'Finished projection to edge basis!')
  end subroutine nte_ndone

  ! -- allocation failure (STAT=istat path) -----------------------------
  subroutine nte_nfatal_alloc() bind(c, name="nte_nfatal_alloc")
    call fatal(caller_nte, 'Memory allocation error.')
  end subroutine nte_nfatal_alloc

end module nodetoedgesolver_acc
