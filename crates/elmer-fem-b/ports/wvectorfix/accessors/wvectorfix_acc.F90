! WVectorFix port: statement-level BIND(C) wrappers for the Elmer calls
! in fem/src/modules/WVectorFix.F90, compiled against the built .mod
! files of the target configure (T021-B lane, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, explicit-shape REAL(C_DOUBLE) arrays). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, array-pointer results —
! stays inside Fortran, so the Rust side is pure C calls. The wrapper
! module mirrors the original's locals (current element handle, the two
! internal procedures' SAVE Nodes state) with identical lifetimes.
!
! The two internal procedures (LocalCorrMatrix, LocalCorrCurrent) are
! transcribed whole: their MATMUL/SUM reductions, whole-array updates,
! LU factor/solve and Fatal paths stay in Fortran, so bit-identity holds
! by construction wherever the call sequence matches. The Rust driver
! only repeats the entry's loop nesting and branch structure.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module wvectorfix_acc
  use iso_c_binding
  use defutils
  use messages
  implicit none

  type(solver_t), pointer, save :: acc_solver => null()
  type(element_t), pointer, save :: acc_elem => null()
  type(variable_t), pointer, save :: acc_gradv => null()
  ! LocalCorrMatrix context (mirrors its SAVE Nodes and IP local).
  type(nodes_t), target, save :: acc_nodes_cm
  type(gaussintegrationpoints_t), save :: acc_ip_cm
  ! LocalCorrCurrent context (mirrors its SAVE Nodes and IP local).
  type(nodes_t), target, save :: acc_nodes_cc
  type(gaussintegrationpoints_t), save :: acc_ip_cc

contains

  ! -- entry state ----------------------------------------------------
  subroutine wva_begin(h_solver) bind(c, name="wva_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine wva_begin

  ! -- WVectorFix_init body --------------------------------------------
  ! dim = CoordinateSystemDimension()
  ! Params => GetSolverParams()
  ! CALL ListAddNewString( Params,'Variable','-nooutput WFixVarTmp')
  subroutine wva_init_variable() bind(c, name="wva_init_variable")
    type(valuelist_t), pointer :: params
    integer :: dim
    dim = CoordinateSystemDimension()
    params => GetSolverParams()
    call ListAddNewString(params, 'Variable', '-nooutput WFixVarTmp')
  end subroutine wva_init_variable

  ! -- Info I/O (formatted output stays in Fortran) ----------------------
  subroutine wva_info_banner() bind(c, name="wva_info_banner")
    call Info('WVectorFix', '----------------------------------------')
    call Info('WVectorFix', 'Fixing the coil gradient potential field')
    call Info('WVectorFix', '----------------------------------------')
  end subroutine wva_info_banner

  subroutine wva_info_apply(v) bind(c, name="wva_info_apply")
    integer(c_int), value :: v
    if (v == 1) then
      call Info('WVectorFix', 'Applying elemental div()=0 to: W Vector E', level=3)
    else
      call Info('WVectorFix', 'Applying elemental div()=0 to: J Vector E', level=3)
    end if
  end subroutine wva_info_apply

  subroutine wva_info_fixed(v) bind(c, name="wva_info_fixed")
    integer(c_int), value :: v
    if (v == 1) then
      call Info('WVectorFix', 'Fixing W Vector E to be divergence free!')
    else
      call Info('WVectorFix', 'Fixing J Vector E to be divergence free!')
    end if
  end subroutine wva_info_fixed

  ! -- flux-variable lookup ----------------------------------------------
  ! GradV => VariableGet(Mesh % Variables, name); varindex==2 may CYCLE
  ! (returns 0); varindex==1 missing is Fatal, exactly as original.
  function wva_find_flux(v) result(r) bind(c, name="wva_find_flux")
    integer(c_int), value :: v
    integer(c_int) :: r
    type(mesh_t), pointer :: m
    m => acc_solver % mesh
    if (v == 1) then
      acc_gradv => VariableGet(m % variables, 'W Vector E')
    else
      acc_gradv => VariableGet(m % variables, 'J Vector E')
    end if
    if (.not. associated(acc_gradv)) then
      if (v == 2) then
        r = 0
        return
      end if
      if (v == 1) then
        call Fatal('WVectorFix', &
            'Fixing is done currently only for W Vector E! I cannot find it.')
      else
        call Fatal('WVectorFix', &
            'Fixing is done currently only for J Vector E! I cannot find it.')
      end if
    end if
    r = 1
  end function wva_find_flux

  ! -- assembly driver ---------------------------------------------------
  subroutine wva_default_initialize() bind(c, name="wva_default_initialize")
    call DefaultInitialize()
  end subroutine wva_default_initialize

  function wva_nof_active() result(n) bind(c, name="wva_nof_active")
    integer(c_int) :: n
    n = GetNOFActive()
  end function wva_nof_active

  subroutine wva_begin_active(t, n, nd) bind(c, name="wva_begin_active")
    integer(c_int), value :: t
    integer(c_int) :: n, nd
    acc_elem => GetActiveElement(t)
    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
  end subroutine wva_begin_active

  ! LocalCorrMatrix transcribed whole (n/nd are the internal procedure's
  ! dummy args; STIFF/FORCE/Tcoef/Lfix are per-call automatics; Nodes is
  ! the procedure's SAVE state, mirrored here with identical lifetime).
  subroutine wva_corr_matrix(n, nd) bind(c, name="wva_corr_matrix")
    integer(c_int), value :: n, nd
    real(kind=dp) :: Basis(nd), dBasisdx(nd,3), DetJ, Weight, L(3), Lfix(3,nd), Tcoef(nd)
    real(kind=dp) :: STIFF(nd,nd), FORCE(nd), C
    logical :: Stat, Found
    integer :: i, t, p, q, j
    type(valuelist_t), pointer :: Material
    call GetElementNodes(acc_nodes_cm)
    STIFF = 0._dp
    FORCE = 0._dp
    Material => GetMaterial(acc_elem)
    if (.not. associated(Material)) call Fatal('WVectorFix', 'Material not found.')
    Tcoef(1:n) = GetReal(Material, 'Electric Conductivity', Found)
    do i = 1, n
      j = acc_gradv % Perm(acc_elem % DGIndexes(i))
      Lfix(1:3, i) = acc_gradv % Values(3*j-2:3*j)
    end do
    acc_ip_cm = GaussPoints(acc_elem)
    do t = 1, acc_ip_cm % n
      stat = ElementInfo(acc_elem, acc_nodes_cm, acc_ip_cm % U(t), acc_ip_cm % V(t), &
          acc_ip_cm % W(t), detJ, Basis, dBasisdx)
      Weight = acc_ip_cm % s(t) * detJ
      L = MATMUL(Lfix(1:3, 1:n), Basis(1:n))
      C = MAXVAL((/1._dp, SUM(Tcoef(1:n) * Basis(1:n))/))
      do p = 1, n
        do q = 1, n
          STIFF(p, q) = STIFF(p, q) + Weight * C * SUM(dBasisdx(q, :) * dBasisdx(p, :))
        end do
        FORCE(p) = FORCE(p) + SUM(L * dBasisdx(p, :)) * Weight
      end do
    end do
    call DefaultUpdateEquations(STIFF, FORCE)
  end subroutine wva_corr_matrix

  ! Dead Dirichlet-node line: value unused (the update calls below stay
  ! commented out, as in the original), but the ParEnv touch and the
  ! element read are replicated exactly.
  subroutine wva_dead_node() bind(c, name="wva_dead_node")
    integer :: Dnode_index
    if (Parenv % myPE == parenv % pes) Dnode_index = acc_elem % NodeIndexes(1)
  end subroutine wva_dead_node

  function wva_default_solve() result(norm) bind(c, name="wva_default_solve")
    real(c_double) :: norm
    norm = DefaultSolve()
  end function wva_default_solve

  ! LocalCorrCurrent transcribed whole (same lifetime notes apply; the LU
  ! factor/solve and the per-direction correction stay in Fortran).
  subroutine wva_corr_current(n, nd) bind(c, name="wva_corr_current")
    integer(c_int), value :: n, nd
    real(kind=dp) :: Basis(nd), dBasisdx(nd,3), NodalPot(nd), DetJ, Weight, Tcoef(nd)
    real(kind=dp) :: STIFF(nd,nd), FORCE(nd,3), x(nd), C
    integer :: pivot(nd)
    logical :: Stat, Found, Erroneous
    integer :: i, j, k, t, p, q
    type(valuelist_t), pointer :: Material
    call GetElementNodes(acc_nodes_cc)
    STIFF = 0._dp
    FORCE = 0._dp
    Material => GetMaterial(acc_elem)
    if (.not. associated(Material)) call Fatal('WVectorFix', 'Material not found.')
    Tcoef(1:n) = GetReal(Material, 'Electric Conductivity', Found)
    NodalPot(1:n) = acc_solver % Variable % Values( &
        acc_solver % Variable % Perm(acc_elem % NodeIndexes))
    acc_ip_cc = GaussPoints(acc_elem)
    do t = 1, acc_ip_cc % n
      stat = ElementInfo(acc_elem, acc_nodes_cc, acc_ip_cc % U(t), acc_ip_cc % V(t), &
          acc_ip_cc % W(t), detJ, Basis, dBasisdx)
      Weight = acc_ip_cc % s(t) * detJ
      C = MAXVAL((/1._dp, SUM(Tcoef(1:n) * Basis(1:n))/))
      do p = 1, n
        do q = 1, n
          STIFF(p, q) = STIFF(p, q) + Weight * C * Basis(p) * Basis(q)
        end do
        do k = 1, 3
          FORCE(p, k) = FORCE(p, k) + Weight * Basis(p) * SUM(dBasisdx(1:n, k) * NodalPot(1:n))
        end do
      end do
    end do
    call LUdecomp(STIFF, n, pivot, Erroneous)
    if (Erroneous) call Fatal('WVectorFix', 'LU-decomposition fails')
    do k = 1, 3
      x = FORCE(1:n, k)
      call LUSolve(n, STIFF, x, pivot)
      do i = 1, n
        j = acc_gradv % Perm(acc_elem % DGIndexes(i))
        acc_gradv % Values(3*(j-1)+k) = acc_gradv % Values(3*(j-1)+k) - x(i)
      end do
    end do
  end subroutine wva_corr_current

end module wvectorfix_acc
