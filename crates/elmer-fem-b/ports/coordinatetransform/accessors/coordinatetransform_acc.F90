! T022-A CoordinateTransform port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/CoordinateTransform.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3), plus statement-for-statement transcriptions of the unit's
! internal procedures ComputeRotM and PolarDecomposition.
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! list-directed timing reports — stays inside Fortran, so the Rust side is
! pure C calls. The wrapper module mirrors the solver's own SAVE state
! (mesh/variable pointers, alpha_ref/beta_ref, Nodes, the polar
! decomposition settings, CoordSys_ijk) with identical lifetimes. All
! floating-point work (the MATMUL reductions, DOT_PRODUCT/SQRT, ISNAN
! checks, the polar decomposition, the PRINT/Warn reports) and the
! whole-array updates stay in Fortran (same compiler), so bit-identity
! holds by construction wherever the call sequence matches. The internal
! VecUtils calls (normalized, jac, crossproduct) resolve to the
! transcription in vecutils_rep.F90, shared with the exported symbols.
!
! Dropped as unobservable (documented, no symbol either way): the
! never-called CONTAINS internal GetElementRotM; the provably-unused
! locals of ComputeRotM (Indexes, l, m, localC, gradv3, tvec, IP).
! Where the original reads undefined automatic-array garbage (alpha/beta
! when neither 'alpha direction'/'alpha' nor 'beta direction'/'beta'
! exists), the port zero-fills instead — deterministic; the proof configs
! always provide the variables, so the defined path is unaffected.
!
! Clean-room code: statement-for-statement wrappers and transcriptions
! around Elmer library calls, no upstream text.
module coordinatetransform_acc
  use iso_c_binding
  use defutils
  use vecutils
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(mesh_t), pointer, save :: acc_mesh => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(variable_t), pointer, save :: acc_rotmvar => null()
  type(variable_t), pointer, save :: acc_alphavecvar => null()
  type(variable_t), pointer, save :: acc_betavecvar => null()
  type(variable_t), pointer, save :: acc_gammavecvar => null()
  type(variable_t), pointer, save :: acc_tmpvar => null()
  real(kind=dp), pointer, save :: acc_alpha_ref(:,:) => null()
  real(kind=dp), pointer, save :: acc_beta_ref(:,:) => null()
  real(kind=dp), save :: acc_coordsys_ijk(3,3) = 0.0_dp
  real(kind=dp), save :: acc_coordsys_ref(3,3) = 0.0_dp
  logical, save :: acc_usepdecomp = .true.
  logical, save :: acc_localsys = .false.
  ! NOTE: the original declares PDDetTol as default REAL (single
  ! precision); the GetConstReal value converts on assignment, and the
  ! tolerance comparison promotes back. Mirrored exactly.
  real(c_float), save :: acc_pddettol = 1.0e-9
  integer, save :: acc_pdmaxiter = 100
  integer, save :: acc_nd = 0
  ! ComputeRotM context (mirrors its SAVE Nodes, consumed before the next
  ! fill, as in the original) and its automatic scratch, sized per element.
  type(nodes_t), target, save :: acc_nodes
  real(kind=dp), allocatable, save :: acc_alpha(:), acc_beta(:)
  real(kind=dp), allocatable, save :: acc_gamma2(:,:)

  ! Host-scope PARAMETERs of the original (visible in ComputeRotM).
  integer, parameter :: ind1(9) = [1,1,1,2,2,2,3,3,3]
  integer, parameter :: ind2(9) = [1,2,3,1,2,3,1,2,3]

contains

  ! -- first-visit state (mirrors the IF(First) block, minus the RotM
  ! -- ALLOCATE, which the Rust side performs via try_reserve) ------------
  function cta_first_n(h_model) result(n) bind(c, name="cta_first_n")
    type(c_ptr), value :: h_model
    integer(c_int) :: n
    call c_f_pointer(h_model, acc_model)
    acc_mesh => acc_model % Mesh
    n = acc_mesh % MaxElementDOFs
  end function cta_first_n

  subroutine cta_first_vars() bind(c, name="cta_first_vars")
    type(valuelist_t), pointer :: solverpars
    logical :: found
    nullify(acc_alpha_ref)
    nullify(acc_beta_ref)
    acc_rotmvar => VariableGet(acc_mesh % Variables, 'RotM E')
    if (.not. associated(acc_rotmvar)) then
      call Fatal('RotMSolver()', 'RotM E variable not found')
    end if
    acc_alphavecvar => VariableGet(acc_mesh % Variables, 'Alpha Vector E')
    if (.not. associated(acc_alphavecvar)) then
      call Fatal('RotMSolver()', 'Alpha Vector E variable not found')
    end if
    acc_betavecvar => VariableGet(acc_mesh % Variables, 'Beta Vector E')
    if (.not. associated(acc_betavecvar)) then
      call Fatal('RotMSolver()', 'Beta Vector E variable not found')
    end if
    acc_gammavecvar => VariableGet(acc_mesh % Variables, 'Gamma Vector E')
    if (.not. associated(acc_gammavecvar)) then
      call Fatal('RotMSolver()', 'Gamma Vector E variable not found')
    end if
    solverpars => GetSolverParams()
    acc_usepdecomp = ListGetLogical(solverpars, 'Use Polar Decomposition', found)
    if (.not. found) acc_usepdecomp = .true.
    acc_pddettol = GetConstReal(solverpars, &
        'Polar Decomposition Determinant Tolerance', found)
    if (.not. found) then
      call Warn('CoordinateTransform', &
          'Polar Decomposition Determinant Tolerance not set.')
      acc_pddettol = 1e-9
    end if
    acc_pdmaxiter = GetInteger(solverpars, &
        'Polar Decomposition Max Iterations', found)
    if (.not. found) then
      call Warn('CoordinateTransform', &
          'Polar Decomposition Max Iteration not set.')
      acc_pdmaxiter = 100
    end if
    acc_coordsys_ijk(1,:) = [1, 0, 0]
    acc_coordsys_ijk(2,:) = [0, 1, 0]
    acc_coordsys_ijk(3,:) = [0, 0, 1]
  end subroutine cta_first_vars

  function cta_nactive() result(n) bind(c, name="cta_nactive")
    integer(c_int) :: n
    n = GetNOFActive()
  end function cta_nactive

  ! -- per-active-element head (main-loop body through the ComputeRotM --
  ! -- local-solution fills); result 1 requests the CYCLE skip ------------
  function cta_elem_begin(q, nn) result(code) bind(c, name="cta_elem_begin")
    integer(c_int), value :: q
    integer(c_int) :: nn
    integer(c_int) :: code
    type(valuelist_t), pointer :: bodyparams
    logical :: found
    acc_elem => GetActiveElement(q)
    nn = GetElementNOFNodes()
    acc_nd = GetElementNOFDOFs()
    bodyparams => GetBodyParams(acc_elem)
    if (.not. associated(bodyparams)) &
        call Fatal('RotMSolver', 'Body Parameters not found')
    call GetConstRealArray(bodyparams, acc_alpha_ref, 'Alpha Reference', found)
    if (.not. found) then
      code = 1
      return
    end if
    if (size(acc_alpha_ref, 1) /= 3) &
        call Fatal('RotMSolver', 'Alpha Reference should have three components!')
    acc_coordsys_ref(1,1:3) = normalized(acc_alpha_ref(1:3,1))
    call GetConstRealArray(bodyparams, acc_beta_ref, 'Beta Reference', found)
    if (.not. found) then
      code = 1
      return
    end if
    if (size(acc_beta_ref, 1) /= 3) &
        call Fatal('RotMSolver', 'Beta Reference should have three components!')
    acc_coordsys_ref(2,1:3) = normalized(acc_beta_ref(1:3,1))
    acc_localsys = ListGetLogical(bodyparams, &
        'Local Coordinate System Beta Reference and Gamma', found)
    if (.not. found) acc_localsys = .false.
    acc_coordsys_ref(3,1:3) = normalized(crossproduct( &
        acc_coordsys_ref(1,1:3), acc_coordsys_ref(2,1:3)))
    call GetElementNodes(acc_nodes)
    if (acc_localsys) then
      if (.not. allocated(acc_gamma2)) allocate(acc_gamma2(3, nn))
      if (size(acc_gamma2, 2) /= nn) then
        deallocate(acc_gamma2)
        allocate(acc_gamma2(3, nn))
      end if
      acc_gamma2 = 0.0_dp
      acc_tmpvar => VariableGet(acc_mesh % Variables, 'coilcurrent e')
      if (associated(acc_tmpvar)) then
        call GetLocalSolution(acc_gamma2, 'coilcurrent e')
      else
        call Fatal('ComputeRotM', 'Local System Beta Reference And Gamma tried but Gamma (coilcurrent e) is not found!')
      end if
    else
      if (.not. allocated(acc_alpha)) allocate(acc_alpha(nn))
      if (size(acc_alpha) /= nn) then
        deallocate(acc_alpha)
        allocate(acc_alpha(nn))
      end if
      if (.not. allocated(acc_beta)) allocate(acc_beta(nn))
      if (size(acc_beta) /= nn) then
        deallocate(acc_beta)
        allocate(acc_beta(nn))
      end if
      acc_alpha = 0.0_dp
      acc_beta = 0.0_dp
      acc_tmpvar => VariableGet(acc_mesh % Variables, 'alpha direction')
      if (associated(acc_tmpvar)) then
        call GetLocalSolution(acc_alpha, 'alpha direction')
      else
        acc_tmpvar => VariableGet(acc_mesh % Variables, 'alpha')
        if (associated(acc_tmpvar)) then
          call GetLocalSolution(acc_alpha, 'alpha')
        end if
      end if
      acc_tmpvar => VariableGet(acc_mesh % Variables, 'beta direction')
      if (associated(acc_tmpvar)) then
        call GetLocalSolution(acc_beta, 'beta direction')
      else
        acc_tmpvar => VariableGet(acc_mesh % Variables, 'beta')
        if (associated(acc_tmpvar)) then
          call GetLocalSolution(acc_beta, 'beta')
        end if
      end if
    end if
    code = 0
  end function cta_elem_begin

  ! -- ComputeRotM node body (exact transcription) ------------------------
  subroutine cta_node(j, nn) bind(c, name="cta_node")
    integer(c_int), value :: j, nn
    real(kind=dp) :: un, vn, wn, detj
    real(kind=dp) :: basis(nn), dbasis(nn, 3)
    real(kind=dp) :: coordsys(3, 3), coordsys2(3, 3), relem(3, 3)
    real(kind=dp) :: jac_ref(3, 3), jac_e(3, 3)
    integer :: k
    logical :: stat
    un = acc_elem % TYPE % NodeU(j)
    vn = acc_elem % TYPE % NodeV(j)
    wn = acc_elem % TYPE % NodeW(j)
    stat = ElementInfo(acc_elem, acc_nodes, un, vn, wn, detj, basis, dbasis)
    if (acc_localsys) then
      coordsys(3,1:3) = matmul(acc_gamma2(1:3,1:nn), basis(1:nn))
      if (any(isnan(coordsys(3,:)))) then
        print *, "Element index = ", GetElementIndex(acc_elem)
        print *, "Element aspect ratio = ", ElementAspectRatio(acc_model, acc_elem)
        call Warn('CoordinateTransform', 'Element coordinate system is NaN, this could be &
          &due to a poor mesh. Let us try to use the degenerate element normal as the local coordinate system alpha vector.')
        coordsys(3,1:3) = NormalOfDegenerateElement(acc_model, acc_elem)
        if (any(isnan(coordsys(3,:)))) &
            call Fatal('CoordinateTransform', 'Degenerate element normal did not work...')
      end if
      coordsys(2,1:3) = acc_coordsys_ref(2,1:3)
      coordsys(1,1:3) = crossproduct(coordsys(2,1:3), coordsys(3,1:3))
    else
      coordsys(1,1:3) = normalized(matmul(acc_alpha(1:nn), dbasis(1:nn,:)))
      if (any(isnan(coordsys(1,:)))) then
        print *, "Element index = ", GetElementIndex(acc_elem)
        print *, "Element aspect ratio = ", ElementAspectRatio(acc_model, acc_elem)
        call Warn('CoordinateTransform', 'Element coordinate system is NaN, this could be &
          &due to a poor mesh. Let us try to use the degenerate element normal as the local coordinate system alpha vector.')
        coordsys(1,1:3) = NormalOfDegenerateElement(acc_model, acc_elem)
        if (any(isnan(coordsys(1,:)))) &
            call Fatal('CoordinateTransform', 'Degenerate element normal did not work...')
      end if
      coordsys(2,1:3) = normalized(matmul(acc_beta(1:nn), dbasis(1:nn,:)))
      if (any(isnan(coordsys(2,:)))) then
        print *, "Element index = ", GetElementIndex(acc_elem)
        print *, "Element aspect ratio = ", ElementAspectRatio(acc_model, acc_elem)
        call Warn('CoordinateTransform', 'Element coordinate system is NaN, this could be &
          &due to a poor mesh. Let us try to use the degenerate element normal as the local coordinate system beta vector.')
        coordsys(2,1:3) = NormalOfDegenerateElement(acc_model, acc_elem)
        if (any(isnan(coordsys(2,:)))) &
            call Fatal('CoordinateTransform', 'Degenerate element normal did not work...')
      end if
      coordsys(3,1:3) = normalized(crossproduct(coordsys(1,1:3), coordsys(2,1:3)))
    end if
    coordsys2 = coordsys
    jac_ref = jac(acc_coordsys_ref, acc_coordsys_ijk)
    jac_e = jac(coordsys, acc_coordsys_ref)
    relem = matmul(jac_ref, jac_e)
    if (acc_usepdecomp) call polardecomp(relem, acc_pdmaxiter, acc_pddettol)
    coordsys2(1,1:3) = [1, 0, 0]
    coordsys2(2,1:3) = [0, 1, 0]
    coordsys2(3,1:3) = [0, 0, 1]
    coordsys2(1,1:3) = matmul(relem, coordsys2(1,1:3))
    coordsys2(2,1:3) = matmul(relem, coordsys2(2,1:3))
    coordsys2(3,1:3) = matmul(relem, coordsys2(3,1:3))
    if (associated(acc_rotmvar)) then
      do k = 1, acc_rotmvar % DOFs
        acc_rotmvar % Values(acc_rotmvar % DOFs * (acc_rotmvar % Perm( &
            acc_elem % DGIndexes(j)) - 1) + k) = relem(ind1(k), ind2(k))
      end do
    end if
    if (associated(acc_alphavecvar)) then
      do k = 1, acc_alphavecvar % DOFs
        acc_alphavecvar % Values(acc_alphavecvar % DOFs * (acc_alphavecvar % Perm( &
            acc_elem % DGIndexes(j)) - 1) + k) = coordsys2(1, k)
      end do
    end if
    if (associated(acc_betavecvar)) then
      do k = 1, acc_betavecvar % DOFs
        acc_betavecvar % Values(acc_betavecvar % DOFs * (acc_betavecvar % Perm( &
            acc_elem % DGIndexes(j)) - 1) + k) = coordsys2(2, k)
      end do
    end if
    if (associated(acc_gammavecvar)) then
      do k = 1, acc_gammavecvar % DOFs
        acc_gammavecvar % Values(acc_gammavecvar % DOFs * (acc_gammavecvar % Perm( &
            acc_elem % DGIndexes(j)) - 1) + k) = coordsys2(3, k)
      end do
    end if
  end subroutine cta_node

  ! -- PolarDecomposition internal (exact transcription; single-precision
  ! -- tolerance dummy exactly as in the original) ------------------------
  subroutine polardecomp(rotmloc, pdmaxiter, pddettol)
    real(kind=dp) :: rotmloc(3, 3), rotmlocinv(3, 3)
    real(kind=dp) :: det
    integer :: i
    real(c_float) :: pddettol
    integer :: pdmaxiter
    logical :: converged
    if (any(isnan(rotmloc))) return
    converged = .false.
    do i = 1, pdmaxiter
      det = det3x3(rotmloc)
      call invertmatrix3x3(rotmloc, rotmlocinv, det)
      rotmloc = (rotmloc + transpose(rotmlocinv)) / 2_dp
      if (abs(det - 1_dp) <= pddettol) then
        converged = .true.
        exit
      end if
    end do
    if (.not. converged) then
      print *, "i:", i, "R*R^T", matmul(rotmloc, transpose(rotmloc))
      print *, "Reference Tolerance: ", pddettol
      print *, "Tolerance: ", abs(det - 1._dp)
      call fatal('PolarDecomposition', &
          'Failed: Polar Decomposition Tolerance')
    end if
  end subroutine polardecomp

  subroutine cta_fatal_alloc() bind(c, name="cta_fatal_alloc")
    call Fatal('RotMSolver', 'Memory allocation error.')
  end subroutine cta_fatal_alloc

end module coordinatetransform_acc
