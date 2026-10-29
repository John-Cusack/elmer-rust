! T022-A DataToFieldSolver port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/DataToFieldSolver.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments. All Fortran-object
! traffic — TYPE(...) handles, OPTIONAL/POINTER dummies, hidden lengths,
! list-directed file I/O, whole-array updates — stays inside Fortran, so
! the Rust side is pure C calls. The wrapper module mirrors the
! original's own locals (entry-level pointers and flags, AsciiPointsToMesh
! SAVE state, BulkAssembly SAVE state, BoundaryLocalMatrix SAVE Nodes)
! with identical lifetimes. All floating-point work (basis evaluation,
! MATMUL/SUM reductions, LOG/EXP, ParallelReduction, the mass-matrix
! accumulation) and all formatted/file I/O stay in Fortran (same
! compiler), so bit-identity holds by construction wherever the call
! sequence matches.
!
! One documented repair: the original's AsciiPointsToMesh declares
! `SAVE :: AllocationsDone, ...` but never sets the flag and allocates
! unconditionally, so a second in-process call aborts (already-allocated
! variable) — undefined-behavior path, the author's evident intent being
! allocate-once. The wrapper allocates only when not yet associated (and
! grows when the mesh needs more), implementing that intent; the defined
! (first-call) path is statement-identical.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module datatoffieldsolver_acc
  use iso_c_binding
  use defutils
  use particleutils
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(valuelist_t), pointer, save :: acc_params => null()
  type(variable_t), pointer, save :: acc_var => null()
  ! Entry-level field/force/weight/mask vectors and perms (mirror the
  ! original's host-associated pointers).
  real(c_double), pointer, save :: acc_fieldvector(:) => null()
  integer, pointer, save :: acc_fieldperm(:) => null()
  real(c_double), pointer, save :: acc_forcevector(:) => null()
  integer, pointer, save :: acc_forceperm(:) => null()
  real(c_double), pointer, save :: acc_weightvector(:) => null()
  integer, pointer, save :: acc_weightperm(:) => null()
  real(c_double), pointer, save :: acc_maskvector(:) => null()
  integer, pointer, save :: acc_maskperm(:) => null()
  ! Entry-level scalars (mirror the original's locals).
  character(len=MAX_NAME_LEN), save :: acc_varname = ''
  character(len=MAX_NAME_LEN), save :: acc_filename = ''
  real(c_double), save :: acc_globalweight = 0.0_c_double
  real(c_double), save :: acc_minmaskval = 0.0_c_double
  real(c_double), save :: acc_maxmaskval = 0.0_c_double
  logical, save :: acc_found = .false.
  logical, save :: acc_found2 = .false.
  logical, save :: acc_globalnormalize = .false.
  logical, save :: acc_givennormalize = .false.
  logical, save :: acc_nodalnormalize = .false.
  logical, save :: acc_constantweightsum = .false.
  logical, save :: acc_mask = .false.
  logical, save :: acc_maskdiffusion = .false.
  logical, save :: acc_uselogscale = .false.
  logical, save :: acc_continuebc = .false.
  logical, save :: acc_revertlogscale = .false.
  ! AsciiPointsToMesh SAVE state (mirrors its SAVE clause).
  integer, save :: acc_ap_elemindex = 0
  real(c_double), pointer, save :: acc_ap_basis(:) => null()
  real(c_double), pointer, save :: acc_ap_dbasisdx(:,:) => null()
  type(nodes_t), target, save :: acc_ap_elemnodes
  ! BulkAssembly SAVE state (mirrors its SAVE clause).
  logical, save :: acc_visited = .false.
  type(nodes_t), target, save :: acc_nodes
  real(c_double), allocatable, save :: acc_stiff(:,:), acc_force(:)
  real(c_double), pointer, save :: acc_basis(:) => null()
  real(c_double), pointer, save :: acc_dbasisdx(:,:) => null()
  real(c_double), allocatable, save :: acc_nodalweight(:)
  real(c_double), allocatable, save :: acc_datadiff(:,:,:)
  integer, pointer, save :: acc_boundaryperm(:) => null()
  integer, save :: acc_boundarynodes = 0
  ! BoundaryLocalMatrix SAVE state (mirrors its SAVE Nodes).
  type(nodes_t), target, save :: acc_blm_nodes, acc_blm_parentnodes

contains

  ! -- entry state ----------------------------------------------------
  subroutine dts_begin(h_solver) bind(c, name="dts_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine dts_begin

  subroutine dts_info_start() bind(c, name="dts_info_start")
    call Info('DataToFieldSolver', '-----------------------------------------', level=4)
    call Info('DataToFieldSolver', 'Resolving field from given data', level=4)
  end subroutine dts_info_start

  ! Entry resolve block: field variable, target variable, weight and
  ! mask variables with the original's Fatal checks in the original
  ! order.
  subroutine dts_resolve() bind(c, name="dts_resolve")
    acc_var => acc_solver % variable
    acc_fieldvector => acc_var % values
    acc_fieldperm => acc_var % perm
    call Info('DataToFieldSolver', 'Fitting to variable: '//trim(acc_var % name), level=6)

    acc_params => GetSolverParams()
    acc_varname = GetString(acc_params, 'Data Field Name', acc_found)
    if (.not. acc_found) acc_varname = GetString(acc_params, 'Target Variable', acc_found)
    if (.not. acc_found) then
      call Fatal('DataToFieldSolver', '> Target Variable < must exist for the solver!')
    end if

    acc_var => VariableGet(acc_solver % mesh % variables, acc_varname)
    if (associated(acc_var)) then
      acc_forcevector => acc_var % values
      acc_forceperm => acc_var % perm
    else
      call Fatal('DataToFieldSolver', 'Variable not present:'//trim(acc_varname))
    end if

    acc_globalweight = GetCReal(acc_params, 'Weight Coefficient', acc_globalnormalize)

    acc_givennormalize = GetLogical(acc_params, 'Normalize by Given Weight', acc_found)
    if (.not. acc_givennormalize) then
      acc_givennormalize = GetLogical(acc_params, 'Normalize Data by Weight', acc_found)
    end if
    if (acc_givennormalize) then
      acc_varname = GetString(acc_params, 'Weight Field Name', acc_found)
      if (.not. acc_found) acc_varname = GetString(acc_params, 'Weight Variable', acc_found)
      if (.not. acc_found) then
        call Fatal('DataToFieldSolver', '> Weight Variable < must exist for the solver!')
      end if
      acc_var => VariableGet(acc_solver % mesh % variables, acc_varname)
      if (associated(acc_var)) then
        acc_weightvector => acc_var % values
        acc_weightperm => acc_var % perm
      else
        call Fatal('DataToFieldSolver', 'Variable not present: '//trim(acc_varname))
      end if
      call Info('DataToFieldSolver', 'Normalizing source using: '//trim(acc_varname), level=6)

      acc_constantweightsum = GetLogical(acc_params, 'Set Constant Weight Sum', acc_found)
    end if

    acc_nodalnormalize = GetLogical(acc_params, 'Normalize by Nodal Weight', acc_found)
    if (acc_nodalnormalize) then
      call Info('DataToFieldSolver', 'Normalizing source using nodal weight', level=6)
    end if

    if (acc_givennormalize .and. acc_nodalnormalize) then
      call Fatal('DataToFieldSolver', 'Normalization cannot be both given and nodal!')
    end if

    acc_maskdiffusion = .false.
    acc_varname = ListGetString(acc_params, 'Mask Field Name', acc_mask)
    if (.not. acc_mask) acc_varname = ListGetString(acc_params, 'Mask Variable', acc_mask)
    if (acc_mask) then
      acc_var => VariableGet(acc_solver % mesh % variables, acc_varname)
      if (associated(acc_var)) then
        acc_maskvector => acc_var % values
        acc_maskperm => acc_var % perm
      else
        call Fatal('DataToFieldSolver', 'Variable not present:'//trim(acc_varname))
      end if
      call Info('DataToFieldSolver', 'Masking source using: '//trim(acc_varname), level=6)

      acc_maxmaskval = ListGetCReal(acc_params, 'Max Mask Value', acc_found)
      if (.not. acc_found) acc_maxmaskval = huge(acc_maxmaskval)
      acc_minmaskval = ListGetCReal(acc_params, 'Min Mask Value', acc_found2)
      if (.not. acc_found2) then
        if (.not. acc_found) then
          acc_minmaskval = 0.0_c_double
        else
          acc_minmaskval = -huge(acc_minmaskval)
        end if
      end if
      acc_maskdiffusion = GetLogical(acc_params, 'Mask Diffusion', acc_found)
      if (acc_maskdiffusion) then
        call Info('DataToFieldSolver', 'Masking diffusion terms', level=6)
      end if
    end if
  end subroutine dts_resolve

  function dts_has_pointfile() result(r) bind(c, name="dts_has_pointfile")
    integer(c_int) :: r
    acc_filename = GetString(acc_params, 'Point Data Filename', acc_found)
    r = merge(1, 0, acc_found)
  end function dts_has_pointfile

  ! Whole AsciiPointsToMesh internal procedure (list-directed file I/O
  ! and ForceVector accumulation stay in Fortran for bit-identity).
  subroutine dts_ascii_points() bind(c, name="dts_ascii_points")
    integer :: i, j, n, dim, no, success, datacolumn, iounit
    integer, pointer :: nodeindexes(:)
    real(c_double) :: sqrtelementmetric, weight, u, v, w, localcoords(3), val, inputdata(4)
    logical :: stat
    real(c_double) :: globalcoords(3)
    type(mesh_t), pointer :: m
    type(element_t), pointer :: currentelement

    m => acc_solver % mesh

    n = m % maxelementnodes
    ! Allocate-once repair of the original's never-guarded SAVE
    ! allocation (documented in the module header): the first-call
    ! path is statement-identical.
    if (.not. associated(acc_ap_basis)) then
      allocate(acc_ap_elemnodes % x(n), acc_ap_elemnodes % y(n), &
          acc_ap_elemnodes % z(n), acc_ap_basis(n), acc_ap_dbasisdx(n, 3))
    else if (size(acc_ap_basis) < n) then
      deallocate(acc_ap_elemnodes % x, acc_ap_elemnodes % y, &
          acc_ap_elemnodes % z, acc_ap_basis, acc_ap_dbasisdx)
      allocate(acc_ap_elemnodes % x(n), acc_ap_elemnodes % y(n), &
          acc_ap_elemnodes % z(n), acc_ap_basis(n), acc_ap_dbasisdx(n, 3))
    end if

    globalcoords(3) = 0.0_c_double
    dim = CoordinateSystemDimension()

    datacolumn = ListGetInteger(acc_params, 'Point Data Column', acc_found)
    if (.not. acc_found) datacolumn = dim + 1

    open(newunit=iounit, file=acc_filename, iostat=success)
    if (success /= 0) then
      call Fatal('DataToFieldSolver', 'Could not open file for reading: '//trim(acc_filename))
    else
      call Info('DataToFieldSolver', 'Reading data points from file: '//trim(acc_filename), level=6)
    end if

    no = 0
    do while (.true.)
      read(iounit, *, iostat=success) inputdata(1:datacolumn)
      if (success /= 0) then
        call Info('DataToFieldSolver', 'End of file after '//i2s(no)//' values')
        if (no == 0) call Fatal('DataToFieldSolver', 'Could not read any points from file!')
        exit
      end if

      no = no + 1

      globalcoords(1:dim) = inputdata(1:dim)
      val = inputdata(datacolumn)

      call LocateParticleInMeshOctree(acc_ap_elemindex, globalcoords, localcoords)

      if (acc_ap_elemindex == 0) cycle

      currentelement => m % elements(acc_ap_elemindex)
      n = currentelement % type % numberofnodes
      nodeindexes => currentelement % nodeindexes
      call GetElementNodes(acc_ap_elemnodes, currentelement)

      u = localcoords(1)
      v = localcoords(2)
      w = localcoords(3)

      stat = ElementInfo(currentelement, acc_ap_elemnodes, u, v, w, sqrtelementmetric, &
          acc_ap_basis, acc_ap_dbasisdx)

      do i = 1, n
        j = acc_fieldperm(nodeindexes(i))

        if (j == 0) cycle

        ! As the weight should be proportional to the particle amount rather than
        ! element volume the weight is not multiplied with local element size!
        ! Note that the weight could be also ~1/r^2 from the nodes etc.
        !-------------------------------------------------------------------------
        weight = acc_ap_basis(i)

        acc_forcevector(j) = acc_forcevector(j) + weight * val
        if (acc_givennormalize) then
          acc_weightvector(j) = acc_weightvector(j) + weight
        end if

      end do
    end do

    close(iounit)

    call Info('DataToFieldSolver', 'Done reading data points', level=12)
  end subroutine dts_ascii_points

  function dts_use_logscale() result(r) bind(c, name="dts_use_logscale")
    integer(c_int) :: r
    acc_uselogscale = GetLogical(acc_params, 'Logarithmic Fitting', acc_found)
    if (acc_uselogscale) then
      call Info('DataToFieldSolver', 'Using logarithmic scale when fitting data!', level=6)
    end if
    r = merge(1, 0, acc_uselogscale)
  end function dts_use_logscale

  function dts_continue_bc() result(r) bind(c, name="dts_continue_bc")
    integer(c_int) :: r
    acc_continuebc = GetLogical(acc_params, 'Continue BC With Boundary Assembly', acc_found)
    if (acc_continuebc) then
      call Info('DataToFieldSolver', 'Using boundary assembly to continue the solution at BCs', level=6)
    end if
    r = merge(1, 0, acc_continuebc)
  end function dts_continue_bc

  subroutine dts_default_initialize() bind(c, name="dts_default_initialize")
    call DefaultInitialize()
  end subroutine dts_default_initialize

  ! Whole BulkAssembly internal procedure (matrix assembly, the
  ! ParallelReduction weight correction, and the diagonal/rhs update
  ! stay in Fortran for bit-identity).
  subroutine dts_bulk_assembly() bind(c, name="dts_bulk_assembly")
    integer, pointer :: indexes(:)
    integer :: i, j, p, q, j2, j3, k, t, n, istat, active, dim, maskactive
    type(element_t), pointer :: element
    type(gaussintegrationpoints_t) :: ip
    character(len=MAX_NAME_LEN) :: boundaryname, diffusivityname
    real(c_double) :: coeff, detj, weightcorr, val, wrhs, wmat, diffmatrix(3, 3)
    real(c_double), pointer :: matvalues(:), matrhs(:)
    real(c_double), pointer :: hwrk(:,:,:) => null()
    integer, pointer :: matdiag(:)
    type(matrix_t), pointer :: stiffmatrix
    logical :: stat, globaldiffuse, localdiffuse
    type(valuelist_t), pointer :: material

    ! Assembly the diffusion part used for regularization
    !----------------------------------------------------------
    coeff = GetCReal(acc_solver % values, 'Diffusion Coefficient', globaldiffuse)
    localdiffuse = .false.

    diffusivityname = GetString(acc_solver % values, 'Diffusivity Name', acc_found)
    if (.not. acc_found) diffusivityname = 'Data Diffusivity'
    localdiffuse = .false.

    active = GetNOFActive()
    stiffmatrix => acc_solver % matrix
    matvalues => stiffmatrix % values
    matrhs => stiffmatrix % rhs
    matdiag => stiffmatrix % diag
    if (.not. associated(stiffmatrix)) then
      call Fatal('DataToFieldSolver', 'StiffMatrix not associated!')
    end if

    dim = CoordinateSystemDimension()

    if (.not. acc_visited) then
      acc_visited = .true.
      n = acc_solver % mesh % maxelementnodes
      allocate(acc_basis(n), acc_dbasisdx(n, 3), acc_force(n), acc_stiff(n, n), &
          acc_datadiff(3, 3, n), stat=istat)
      if (istat /= 0) call Fatal('DataToFieldSolver', 'Allocation error 1 in BulkAssembly!')

      n = stiffmatrix % numberofrows
      allocate(acc_nodalweight(n), stat=istat)
      if (istat /= 0) call Fatal('DataToFieldSolver', 'Allocation error 2 in BulkAssembly!')

      if (.not. acc_continuebc) then
        n = acc_solver % mesh % numberofnodes
        allocate(acc_boundaryperm(n))
        acc_boundaryperm = 0
        acc_boundarynodes = 0
        boundaryname = ComponentName(acc_solver % variable)
        boundaryname = trim(boundaryname)//' continue'
        call MakePermUsingMask(CurrentModel, acc_solver, acc_solver % mesh, boundaryname, &
            .false., acc_boundaryperm, acc_boundarynodes)
        if (acc_boundarynodes > 0) then
          write(message, '(A,I0)') 'Nodes with BC > '//trim(boundaryname)//' <  true: ', acc_boundarynodes
          call Info('DataToFieldSolver', message, level=6)
        end if
      end if

    end if

    acc_nodalweight = 0.0_c_double

    do t = 1, active

      element => GetActiveElement(t)
      n = GetElementNOFNodes(element)
      indexes => element % nodeindexes

      call GetElementNodes(acc_nodes, element)
      acc_stiff = 0.0_c_double
      acc_force = 0.0_c_double

      if (.not. globaldiffuse) then
        material => GetMaterial()
        call ListGetRealArray(material, diffusivityname, hwrk, n, indexes, localdiffuse)
        if (localdiffuse) then
          acc_datadiff = 0.0_c_double
          if (size(hwrk, 1) == 1) then
            do i = 1, 3
              acc_datadiff(i, i, 1:n) = hwrk(1, 1, 1:n)
            end do
          else if (size(hwrk, 2) == 1) then
            do i = 1, min(3, size(hwrk, 1))
              acc_datadiff(i, i, 1:n) = hwrk(i, 1, 1:n)
            end do
          else
            do i = 1, min(3, size(hwrk, 1))
              do j = 1, min(3, size(hwrk, 2))
                acc_datadiff(i, j, 1:n) = hwrk(i, j, 1:n)
              end do
            end do
          end if
        end if
      end if


      ! Numerical integration:
      !----------------------
      ip = GaussPoints(element)
      do k = 1, ip % n
        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo(element, acc_nodes, ip % u(k), ip % v(k), &
            ip % w(k), detj, acc_basis, acc_dbasisdx)

        ! Finally, the elemental matrix & vector:
        !----------------------------------------
        do i = 1, n

          if (.not. acc_continuebc .and. acc_boundarynodes > 0) then
            if (acc_boundaryperm(indexes(i)) > 0) cycle
          end if

          ! Compute the rowsum that is used in the normalization
          !-----------------------------------------------------
          val = ip % s(k) * detj * acc_basis(i)
          j = indexes(i)
          acc_nodalweight(j) = acc_nodalweight(j) + val


          ! Compute the local conductivity tensor
          ! -------------------------------
          if (localdiffuse) then
            do p = 1, dim
              do q = 1, dim
                diffmatrix(p, q) = sum(acc_datadiff(p, q, 1:n) * acc_basis(1:n))
              end do
            end do
          end if

          ! This condition should remove the diffusion for proper data, and
          ! use it only for outlier data.
          !----------------------------------------------------------------
          if (acc_maskdiffusion) then
            j2 = indexes(i)
            if (associated(acc_maskperm)) j2 = acc_maskperm(j2)
            val = acc_maskvector(j2)
            if (.not. (val < acc_minmaskval .or. val > acc_maxmaskval)) cycle
          end if

          ! This condition removes the natural boundary condition that would
          ! try to fix the normal gradient of the field to zero.
          ! Does not seem to work though...
          !--------------------------------------------------------------------

          if (globaldiffuse) then
            do j = 1, n
              acc_stiff(i, j) = acc_stiff(i, j) + ip % s(k) * detj * &
                  coeff * sum(acc_dbasisdx(i, 1:dim) * acc_dbasisdx(j, 1:dim))
            end do
          else if (localdiffuse) then
            do j = 1, n
              acc_stiff(i, j) = acc_stiff(i, j) + ip % s(k) * detj * &
                  sum(matmul(diffmatrix(1:dim, 1:dim), acc_dbasisdx(j, 1:dim)) * acc_dbasisdx(i, 1:dim))
            end do
          end if

        end do
      end do

      call DefaultUpdateEquations(acc_stiff, acc_force)
    end do

    maskactive = 0
    if (acc_givennormalize .and. acc_constantweightsum) then
      !-----------------------------------------------------------------------
      ! Set the weight to the diagonal i.e. make the mass matrix contribution
      ! The data is normalized so that if it would be constant it would yield the
      ! same equation as the normal one and the weights would also be constant.
      ! This way diffusion will not depend on the amount of data, whether
      ! that is desirable, or not, I don't know.
      !-----------------------------------------------------------------------
      block
        real(c_double) :: s1, s2
        s1 = sum(acc_nodalweight)
        s2 = sum(acc_weightvector)
        s1 = ParallelReduction(s1)
        s2 = ParallelReduction(s2)
        weightcorr = s1 / s2
      end block
    else
      weightcorr = 1.0_c_double
    end if

    do i = 1, acc_solver % mesh % numberofnodes

      ! If some values are not to be trusted don't include them in the allocation process
      !----------------------------------------------------------------------------------
      if (acc_mask) then
        j = i
        if (associated(acc_maskperm)) j = acc_maskperm(i)
        if (j > 0) then
          val = acc_maskvector(j)
          if (val < acc_minmaskval .or. val > acc_maxmaskval) then
            maskactive = maskactive + 1
            cycle
          end if
        end if
      end if

      ! field to be solved for
      j = acc_fieldperm(i)
      if (j == 0) cycle

      ! force vector on the r.h.s.
      j2 = i
      if (associated(acc_forceperm)) j2 = acc_forceperm(j2)
      if (j2 > 0) then
        val = acc_forcevector(j2)
        if (acc_uselogscale) val = log(val)
      else
        val = 0.0_c_double
      end if

      if (acc_givennormalize) then
        j3 = i
        if (associated(acc_weightperm)) j3 = acc_weightperm(j3)
        if (j3 > 0) then
          wmat = weightcorr * acc_weightvector(j3)
        else
          wmat = 0.0_c_double
        end if
        wrhs = weightcorr * val

      else if (acc_nodalnormalize) then
        wmat = acc_nodalweight(i)
        wrhs = val

      else
        wmat = acc_nodalweight(i)
        wrhs = acc_nodalweight(i) * val
      end if

      if (acc_globalnormalize) then
        wmat = acc_globalweight * wmat
        wrhs = acc_globalweight * wrhs
      end if

      k = matdiag(j)
      matvalues(k) = matvalues(k) + wmat
      matrhs(j) = matrhs(j) + wrhs
    end do


    if (acc_mask) then
      write(message, '(A,I0,A,I0,A)') 'Mask is active for ', maskactive, &
          ' nodes (out of ', size(acc_maskvector), ')'
      call Info('DataToFieldSolver', message, level=4)
    end if

  end subroutine dts_bulk_assembly

  subroutine dts_default_finish_bulk() bind(c, name="dts_default_finish_bulk")
    call DefaultFinishBulkAssembly()
  end subroutine dts_default_finish_bulk

  ! Whole BoundaryAssembly internal procedure plus its BoundaryLocalMatrix
  ! (parent-element assembly stays in Fortran for bit-identity).
  subroutine dts_boundary_assembly() bind(c, name="dts_boundary_assembly")
    integer :: t, n, np
    type(element_t), pointer :: element, parentelement
    type(valuelist_t), pointer :: bc
    character(len=MAX_NAME_LEN) :: condname

    condname = ComponentName(acc_solver % variable)
    condname = trim(condname)//' continue'

    do t = 1, GetNOFBoundaryElements()

      element => GetBoundaryElement(t)
      if (.not. ActiveBoundaryElement(element)) cycle

      bc => GetBC(element)
      if (.not. GetLogical(bc, condname, acc_found)) cycle

      n = GetElementNOFNodes(element)

      parentelement => element % boundaryinfo % left
      if (.not. associated(parentelement)) then
        call Fatal('DataToFieldSolver', 'Could not find parent element!')
      else if (associated(element % boundaryinfo % right)) then
        call Fatal('DataToFieldSolver', 'This does not make sense for internal BC!')
      end if

      np = GetElementNOFNodes(parentelement)

      call BoundaryLocalMatrix(element, parentelement, n, np)
    end do
  end subroutine dts_boundary_assembly

  subroutine BoundaryLocalMatrix(element, parentelement, n, np)
    integer :: n, np
    type(element_t), target :: element, parentelement
    real(c_double), target :: stiff(np, np), force(np)
    real(c_double), pointer :: a(:,:), m(:,:)
    real(c_double) :: basis(n), dbasisdx(n, 3), s, detj, u, v, w
    real(c_double) :: parentbasis(np), parentdbasisdx(np, 3), diffmatrix(3, 3)
    real(c_double) :: nrm(3), coeff
    real(c_double) :: datadiffusivity(3, 3, np)
    real(c_double), pointer :: hwrk(:,:,:)
    logical :: stat, found
    integer :: i, j, p, q, t, dim
    type(gaussintegrationpoints_t) :: ip
    logical :: localdiffuse, globaldiffuse
    type(valuelist_t), pointer :: material

    dim = CoordinateSystemDimension()

    coeff = GetCReal(acc_solver % values, 'Diffusion Coefficient', globaldiffuse)
    localdiffuse = .false.

    call GetElementNodes(acc_blm_nodes, element)
    call GetElementNodes(acc_blm_parentnodes, parentelement)

    stiff = 0.0_c_double
    force = 0.0_c_double

    if (.not. globaldiffuse) then
      material => GetMaterial(parentelement)
      call ListGetRealArray(material, 'Data Diffusivity', hwrk, np, &
          parentelement % nodeindexes, localdiffuse)
      if (localdiffuse) then
        datadiffusivity = 0.0_c_double
        if (size(hwrk, 1) == 1) then
          do i = 1, 3
            datadiffusivity(i, i, 1:np) = hwrk(1, 1, 1:np)
          end do
        else if (size(hwrk, 2) == 1) then
          do i = 1, min(3, size(hwrk, 1))
            datadiffusivity(i, i, 1:np) = hwrk(i, 1, 1:np)
          end do
        else
          do i = 1, min(3, size(hwrk, 1))
            do j = 1, min(3, size(hwrk, 2))
              datadiffusivity(i, j, 1:np) = hwrk(i, j, 1:np)
            end do
          end do
        end if
      end if
    end if

    !Numerical integration:
    !----------------------
    ip = GaussPoints(element)
    do t = 1, ip % n
      u = ip % u(t)
      v = ip % v(t)
      w = ip % w(t)

      stat = ElementInfo(element, acc_blm_nodes, u, v, w, detj, basis, dbasisdx)

      s = ip % s(t) * detj
      nrm = NormalVector(element, acc_blm_nodes, u, v, .true.)

      call GetParentUVW(element, n, parentelement, np, u, v, w, basis)

      stat = ElementInfo(parentelement, acc_blm_parentnodes, u, v, w, detj, &
          parentbasis, parentdbasisdx)


      ! Compute the local conductivity tensor
      ! -------------------------------
      if (localdiffuse) then
        do p = 1, dim
          do q = 1, dim
            diffmatrix(p, q) = sum(datadiffusivity(p, q, 1:np) * parentbasis(1:np))
          end do
        end do
      end if


      do p = 1, np
        do q = 1, np
          if (globaldiffuse) then
            stiff(p, q) = stiff(p, q) - s * parentbasis(p)* &
                coeff * sum(parentdbasisdx(q, 1:dim) * nrm(1:dim))
          else if (localdiffuse) then
            stiff(p, q) = stiff(p, q) - s * parentbasis(p) * &
                sum(matmul(diffmatrix, dbasisdx(q, :)) * nrm(1:dim))
          end if
        end do
      end do
    end do

    call DefaultUpdateEquations(stiff, force, UElement=parentelement)
  end subroutine BoundaryLocalMatrix

  subroutine dts_default_finish() bind(c, name="dts_default_finish")
    call DefaultFinishAssembly()
  end subroutine dts_default_finish

  subroutine dts_default_dirichlet() bind(c, name="dts_default_dirichlet")
    call DefaultDirichletBCs()
  end subroutine dts_default_dirichlet

  function dts_revert_log() result(r) bind(c, name="dts_revert_log")
    integer(c_int) :: r
    acc_revertlogscale = GetLogical(acc_params, 'Revert Logarithmic Fitting', acc_found)
    r = merge(1, 0, acc_revertlogscale)
  end function dts_revert_log

  function dts_solve_revert() result(norm) bind(c, name="dts_solve_revert")
    real(c_double) :: norm
    call Info('DataToFieldSolver', 'Reverting from logarithmic scale!', level=6)
    call ListAddLogical(acc_solver % values, 'Skip Compute Nonlinear Change', .true.)
    norm = DefaultSolve()
    acc_fieldvector = exp(acc_fieldvector)
    call ComputeChange(acc_solver, .false.)
  end function dts_solve_revert

  function dts_solve() result(norm) bind(c, name="dts_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function dts_solve

  subroutine dts_info_done() bind(c, name="dts_info_done")
    call Info('DataToFieldSolver', 'All done', level=4)
    call Info('DataToFieldSolver', '-----------------------------------------', level=4)
  end subroutine dts_info_done

  subroutine dts_fatal_alloc() bind(c, name="dts_fatal_alloc")
    call Fatal('DataToFieldSolver', 'Memory allocation error.')
  end subroutine dts_fatal_alloc

  ! -- DataToFieldSolver_init body --------------------------------------
  subroutine dts_init_begin(h_solver) bind(c, name="dts_init_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine dts_init_begin

  subroutine dts_init_body() bind(c, name="dts_init_body")
    type(valuelist_t), pointer :: solverpars
    logical :: givennormalize, havefile
    character(len=MAX_NAME_LEN) :: varname

    solverpars => GetSolverParams()

    havefile = ListCheckPresent(solverpars, 'Point Data Filename')

    givennormalize = GetLogical(solverpars, 'Normalize by Given Weight', acc_found)
    if (.not. givennormalize) givennormalize = GetLogical(solverpars, 'Normalize Data by Weight', acc_found)

    ! If the file is given then the two following fields will be created internally
    ! and one can allocate them here, if not given otherwise.
    !-----------------------------------------------------------------------------
    varname = GetString(solverpars, 'Target Variable', acc_found)
    if (.not. acc_found .and. havefile) then
      call ListAddString(solverpars, &
          NextFreeKeyword('Exported Variable', solverpars), 'Point Data')
      call ListAddString(solverpars, 'Target Variable', 'Point Data')
      call Info('DataToFieldSolver_init', 'Creating > Point Data < as exported variable')
    end if

    if (givennormalize .and. havefile) then
      varname = GetString(solverpars, 'Weight Variable', acc_found)
      if (.not. acc_found) then
        call ListAddString(solverpars, &
            NextFreeKeyword('Exported Variable', solverpars), 'Point Weight')
        call ListAddString(solverpars, 'Weight Variable', 'Point Weight')
        call Info('DataToFieldSolver_init', 'Creating > Point Weight < as exported variable')
      end if
    end if

    ! If the field is not given, create it
    ! A different default name if point data set is used.
    !-----------------------------------------------------
    varname = GetString(solverpars, 'Variable', acc_found)
    if (.not. acc_found) then
      if (havefile) then
        call ListAddString(solverpars, 'Variable', 'Point Fit')
      else
        call ListAddString(solverpars, 'Variable', 'Fit')
      end if
    end if


    call ListAddInteger(solverpars, 'Time derivative order', 0)

    ! Add some cheap linear system defaults: bicgstab + none (diagonal when scaled)
    !------------------------------------------------------------------------
    call ListAddNewString(solverpars, 'Linear System Solver', 'Iterative')
    call ListAddNewString(solverpars, 'Linear System Iterative Method', 'bicgstab')
    call ListAddNewString(solverpars, 'Linear System Preconditioning', 'none')
    call ListAddNewInteger(solverpars, 'Linear System Max Iterations', 1000)
    call ListAddNewInteger(solverpars, 'Linear System Residual Output', 20)
    call ListAddNewConstReal(solverpars, 'Linear System Convergence Tolerance', 1.0e-10_dp)
  end subroutine dts_init_body

end module datatoffieldsolver_acc
