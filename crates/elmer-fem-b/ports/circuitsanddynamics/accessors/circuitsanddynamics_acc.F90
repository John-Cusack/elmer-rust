!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors:   Eelis Takala(Trafotek Oy) and Juha Ruokolainen(CSC)
! *  Emails:    eelis.takala@trafotek.fi and Juha.Ruokolainen@csc.fi
! *  Web:       http://www.trafotek.fi and http://www.csc.fi/elmer
! *  Addresses: Trafotek Oy
! *             Kaarinantie 700
! *             Turku
! *
! *             and
! *
! *             CSC - IT Center for Science Ltd.
! *             Keilaranta 14
! *             02101 Espoo, Finland 
! *
! *  Original Date: October 2015
! *
! *****************************************************************************/
 
! T022-A port of fem/src/modules/CircuitsAndDynamics.F90: statement-level
! BIND(C) wrappers plus the original CONTAINS kernels, lifted to module scope
! with their SAVE state and call order intact. The six entries live in Rust;
! every component access and every Elmer library call crosses through here,
! so the Rust side holds no layout assumption and no pinned offset.
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text beyond the lifted kernel bodies, which compile
! under the same gfortran with the same flags as the reference module.
MODULE cda_acc
  USE ISO_C_BINDING
  USE Types
  USE DefUtils
  USE CircuitUtils
  USE CircuitsMod
  USE CircMatInitMod
  USE MGDynMaterialUtils
  USE Messages
  IMPLICIT NONE
  ! Entry context: replaces the entries' host association for the lifted
  ! kernels. Set once per entry call from the Rust side.
  TYPE(Solver_t), POINTER, SAVE :: Solver => NULL()
  TYPE(Model_t), POINTER, SAVE :: Model => NULL()
  REAL(KIND=dp), SAVE :: dt = 0.0_dp
  LOGICAL, SAVE :: TransientSimulation = .FALSE.
  TYPE(CircuitModel_t), POINTER, SAVE :: Ckt => NULL()
  LOGICAL, SAVE :: EigenSystem = .FALSE.
  TYPE(Variable_t), POINTER, SAVE :: LagrangeVar => NULL()
  REAL(KIND=dp), ALLOCATABLE, SAVE :: out_crt(:)
  CHARACTER(LEN=MAX_NAME_LEN), SAVE :: out_CktPrefix = 'res:'
CONTAINS
  ! -- entry context --------------------------------------------------
  SUBROUTINE cda_begin(h_model, h_solver, dt_val, transient_val) &
      BIND(C, NAME="cda_begin")
    TYPE(C_PTR), VALUE :: h_model, h_solver
    REAL(C_DOUBLE), VALUE :: dt_val
    INTEGER(C_INT), VALUE :: transient_val
    CALL C_F_POINTER(h_model, Model)
    CALL C_F_POINTER(h_solver, Solver)
    dt = dt_val
    TransientSimulation = (transient_val /= 0)
  END SUBROUTINE cda_begin
  SUBROUTINE cda_claim_only() BIND(C, NAME="cda_claim_only")
    Ckt => GetCircuitModel(Solver)
  END SUBROUTINE cda_claim_only
  SUBROUTINE cda_claim_set() BIND(C, NAME="cda_claim_set")
    Ckt => GetCircuitModel(Solver)
    CALL SetCircuitModel(Ckt)
  END SUBROUTINE cda_claim_set
  SUBROUTINE cda_solver_bind(which) BIND(C, NAME="cda_solver_bind")
    INTEGER(C_INT), VALUE :: which
    IF (which == 0) THEN
      CALL CircuitSolverBind(Solver,'Pre Solvers','CircuitsAndDynamics_init')
    ELSE IF (which == 1) THEN
      CALL CircuitSolverBind(Solver,'Pre Solvers','CircuitsAndDynamicsHarmonic_init')
    ELSE
      CALL CircuitSolverBind(Solver,'Post Solvers','CircuitsOutput_init')
    END IF
  END SUBROUTINE cda_solver_bind
  ! -- _init value lists ----------------------------------------------
  SUBROUTINE cda_add_var(which) BIND(C, NAME="cda_add_var")
    INTEGER(C_INT), VALUE :: which
    IF (which == 0) THEN
      CALL ListAddNewString(Solver % Values,'Variable','-global ckt')
    ELSE
      CALL ListAddNewString(Solver % Values,'Variable','-global cmplxckt')
    END IF
  END SUBROUTINE cda_add_var
  SUBROUTINE cda_add_var_output() BIND(C, NAME="cda_add_var_output")
    CALL ListAddNewLogical(Solver % Values,'Variable Output',.FALSE.)
  END SUBROUTINE cda_add_var_output
  SUBROUTINE cda_add_time_order() BIND(C, NAME="cda_add_time_order")
    CALL ListAddNewInteger(Solver % Values,'Time Derivative Order',1)
  END SUBROUTINE cda_add_time_order
  SUBROUTINE cda_add_no_matrix() BIND(C, NAME="cda_add_no_matrix")
    CALL ListAddNewLogical(Solver % Values,'No Matrix',.TRUE.)
  END SUBROUTINE cda_add_no_matrix
  FUNCTION cda_sim_has_rotor_angle() RESULT(r) BIND(C, NAME="cda_sim_has_rotor_angle")
    INTEGER(C_INT) :: r
    r = MERGE(1, 0, ListCheckPresent(Model % Simulation,'Rotor Angle'))
  END FUNCTION cda_sim_has_rotor_angle
  FUNCTION cda_rot_machine() RESULT(r) BIND(C, NAME="cda_rot_machine")
    INTEGER(C_INT) :: r
    LOGICAL :: RotMachine, Found
    RotMachine = ListGetLogical(Solver % Values,'Rotating Machine',Found)
    IF (.NOT. Found) THEN
      RotMachine = ListGetLogicalAnyBC(Model,'Rotational Projector') .OR. &
          ListGetLogicalAnyBC(Model,'Anti Rotational Projector')
    END IF
    r = MERGE(1, 0, RotMachine)
  END FUNCTION cda_rot_machine
  SUBROUTINE cda_add_exported(which) BIND(C, NAME="cda_add_exported")
    INTEGER(C_INT), VALUE :: which
    IF (which == 0) THEN
      CALL ListAddString(Solver % Values, &
          NextFreeKeyword('Exported Variable ',Solver % Values), '-global Rotor Angle')
    ELSE
      CALL ListAddString(Solver % Values, &
          NextFreeKeyword('Exported Variable ',Solver % Values), '-global Rotor Velo')
    END IF
  END SUBROUTINE cda_add_exported
  FUNCTION cda_n_circuits() RESULT(n) BIND(C, NAME="cda_n_circuits")
    INTEGER(C_INT) :: n
    n = Ckt % n_Circuits
  END FUNCTION cda_n_circuits
  FUNCTION cda_hm_cm_ready() RESULT(r) BIND(C, NAME="cda_hm_cm_ready")
    INTEGER(C_INT) :: r
    TYPE(Matrix_t), POINTER :: CM
    CM => Ckt % CircuitMatrix
    r = 0
    IF (.NOT. ASSOCIATED(CM)) RETURN
    IF (.NOT. ASSOCIATED(CM % Values)) RETURN
    IF (SIZE(CM % Values) <= 0) RETURN
    r = 1
  END FUNCTION cda_hm_cm_ready
  ! -- build block ----------------------------------------------------
  SUBROUTINE cda_check_stale() BIND(C, NAME="cda_check_stale")
    CALL CircuitsCheckStale()
  END SUBROUTINE cda_check_stale
  FUNCTION cda_built_nm() RESULT(r) BIND(C, NAME="cda_built_nm")
    INTEGER(C_INT) :: r
    r = Ckt % BuiltNm
  END FUNCTION cda_built_nm
  SUBROUTINE cda_info_init(which) BIND(C, NAME="cda_info_init")
    INTEGER(C_INT), VALUE :: which
    CHARACTER(LEN=MAX_NAME_LEN) :: Caller
    IF (which == 1) THEN
      Caller = 'CircuitsAndDynamicsHarmonic'
      CALL Info(Caller,'Initializing electric circuits for harmonic simulation',Level=6)
    ELSE IF (TransientSimulation) THEN
      Caller = 'CircuitsAndDynamics'
      CALL Info(Caller,'Initializing electric circuits for transient simulation',Level=6)
    ELSE
      Caller = 'CircuitsAndDynamics'
      CALL Info(Caller,'Initializing electric circuits for steady state simulation',Level=6)
    END IF
    IF (Solver % Parallel) THEN
      CALL Info(Caller,'Assuming parallel electric circuits',Level=12)
    ELSE
      CALL Info(Caller,'Assuming serial electric circuits',Level=12)
    END IF
  END SUBROUTINE cda_info_init
  SUBROUTINE cda_set_harmonic(v) BIND(C, NAME="cda_set_harmonic")
    INTEGER(C_INT), VALUE :: v
    Ckt % Harmonic = (v /= 0)
  END SUBROUTINE cda_set_harmonic
  SUBROUTINE cda_add_comp_to_bodies() BIND(C, NAME="cda_add_comp_to_bodies")
    CALL AddComponentsToBodyLists()
  END SUBROUTINE cda_add_comp_to_bodies
  SUBROUTINE cda_circ_reset() BIND(C, NAME="cda_circ_reset")
    ! Only the total; the count is (re)read from the container after
    ! AllocateCircuitsList() populates it from the MATC definitions.
    Ckt % Circuit_tot_n = 0
  END SUBROUTINE cda_circ_reset
  SUBROUTINE cda_find_asolver() BIND(C, NAME="cda_find_asolver")
    TYPE(Solver_t), POINTER :: ASolver
    IF (Ckt % Harmonic) THEN
      ASolver => FindCircuitASolver(Solver,.TRUE.,'CircuitsAndDynamicsHarmonic')
    ELSE
      ASolver => FindCircuitASolver(Solver,.FALSE.,'CircuitsAndDynamics')
    END IF
    Ckt % ASolver => ASolver
  END SUBROUTINE cda_find_asolver
  FUNCTION cda_parallel_mismatch() RESULT(r) BIND(C, NAME="cda_parallel_mismatch")
    INTEGER(C_INT) :: r
    r = MERGE(1, 0, Solver % Parallel .NEQV. Ckt % ASolver % Parallel)
  END FUNCTION cda_parallel_mismatch
  FUNCTION cda_solver_parallel() RESULT(r) BIND(C, NAME="cda_solver_parallel")
    INTEGER(C_INT) :: r
    r = MERGE(1, 0, Solver % Parallel)
  END FUNCTION cda_solver_parallel
  SUBROUTINE cda_warn_parallel_fix() BIND(C, NAME="cda_warn_parallel_fix")
    CHARACTER(LEN=MAX_NAME_LEN) :: Caller
    IF (Ckt % Harmonic) THEN
      Caller = 'CircuitsAndDynamicsHarmonic'
    ELSE
      Caller = 'CircuitsAndDynamics'
    END IF
    CALL Warn(Caller,'Conflicting parallel status for circuit and A solver!')
    Solver % Parallel = .TRUE.
    Ckt % ASolver % Parallel = .TRUE.
  END SUBROUTINE cda_warn_parallel_fix
  SUBROUTINE cda_set_ckt_parallel() BIND(C, NAME="cda_set_ckt_parallel")
    Ckt % Parallel = Solver % Parallel
  END SUBROUTINE cda_set_ckt_parallel
  SUBROUTINE cda_check_multiplier_unique() BIND(C, NAME="cda_check_multiplier_unique")
    IF (Ckt % Harmonic) THEN
      Ckt % MultName = TRIM(LagrangeMultiplierName(Ckt % ASolver))
      CALL CheckCircuitMultiplierUnique(Ckt,'CircuitsAndDynamicsHarmonic')
    ELSE
      Ckt % MultName = TRIM(LagrangeMultiplierName(Ckt % ASolver))
      CALL CheckCircuitMultiplierUnique(Ckt,'CircuitsAndDynamics')
    END IF
  END SUBROUTINE cda_check_multiplier_unique
  SUBROUTINE cda_set_matc_prefix() BIND(C, NAME="cda_set_matc_prefix")
    IF (Ckt % Harmonic) THEN
      CALL SetCircuitMatcPrefix(Ckt,Solver,'CircuitsAndDynamicsHarmonic')
    ELSE
      CALL SetCircuitMatcPrefix(Ckt,Solver,'CircuitsAndDynamics')
    END IF
  END SUBROUTINE cda_set_matc_prefix
  SUBROUTINE cda_allocate_circuits_list() BIND(C, NAME="cda_allocate_circuits_list")
    CALL AllocateCircuitsList()
  END SUBROUTINE cda_allocate_circuits_list
  SUBROUTINE cda_set_boundary_areas() BIND(C, NAME="cda_set_boundary_areas")
    CALL SetBoundaryAreasToValueLists()
  END SUBROUTINE cda_set_boundary_areas
  SUBROUTINE cda_info_circuit(p, n) BIND(C, NAME="cda_info_circuit")
    INTEGER(C_INT), VALUE :: p, n
    CALL Info('CircuitsAndDynamics','Initializing circuit '//I2S(p)// &
        ' with '//I2S(n)//' variables!',Level=6)
  END SUBROUTINE cda_info_circuit
  FUNCTION cda_n_circ_vars(p) RESULT(n) BIND(C, NAME="cda_n_circ_vars")
    INTEGER(C_INT), VALUE :: p
    INTEGER(C_INT) :: n
    n = GetNofCircVariables(p)
  END FUNCTION cda_n_circ_vars
  SUBROUTINE cda_allocate_circuit(p) BIND(C, NAME="cda_allocate_circuit")
    INTEGER(C_INT), VALUE :: p
    CALL AllocateCircuit(p)
  END SUBROUTINE cda_allocate_circuit
  FUNCTION cda_count_comps(p, n) RESULT(nc) BIND(C, NAME="cda_count_comps")
    INTEGER(C_INT), VALUE :: p, n
    INTEGER(C_INT) :: nc
    nc = CountNofCircComponents(p, n)
  END FUNCTION cda_count_comps
  SUBROUTINE cda_alloc_comps(p, ncomp) BIND(C, NAME="cda_alloc_comps")
    INTEGER(C_INT), VALUE :: p, ncomp
    INTEGER :: astat
    ! The count is a field the later calls loop over; the original assigns
    ! it from CountNofCircComponents() directly ahead of the allocate.
    Ckt % Circuits(p) % n_comp = ncomp
    ALLOCATE(Ckt % Circuits(p) % Components(ncomp), STAT=astat)
    IF (astat /= 0) THEN
      IF (Ckt % Harmonic) THEN
        CALL Fatal('CircuitsAndDynamicsHarmonic','Memory allocation error.')
      ELSE
        CALL Fatal('CircuitsAndDynamics','Memory allocation error.')
      END IF
    END IF
  END SUBROUTINE cda_alloc_comps
  SUBROUTINE cda_setup_circuit(p, harm, par) BIND(C, NAME="cda_setup_circuit")
    INTEGER(C_INT), VALUE :: p, harm, par
    Ckt % Circuits(p) % Harmonic = (harm /= 0)
    Ckt % Circuits(p) % Parallel = (par /= 0)
    Ckt % Circuits(p) % Asolver => Ckt % ASolver
  END SUBROUTINE cda_setup_circuit
  SUBROUTINE cda_read_circuit_vars(p) BIND(C, NAME="cda_read_circuit_vars")
    INTEGER(C_INT), VALUE :: p
    CALL ReadCircuitVariables(p)
  END SUBROUTINE cda_read_circuit_vars
  SUBROUTINE cda_read_components(p) BIND(C, NAME="cda_read_components")
    INTEGER(C_INT), VALUE :: p
    CALL ReadComponents(p)
  END SUBROUTINE cda_read_components
  SUBROUTINE cda_add_comp_values(p) BIND(C, NAME="cda_add_comp_values")
    INTEGER(C_INT), VALUE :: p
    CALL AddComponentValuesToLists(p)
  END SUBROUTINE cda_add_comp_values
  SUBROUTINE cda_add_bare_vars(p) BIND(C, NAME="cda_add_bare_vars")
    INTEGER(C_INT), VALUE :: p
    CALL AddBareCircuitVariables(p)
  END SUBROUTINE cda_add_bare_vars
  SUBROUTINE cda_read_coef_matrices(p) BIND(C, NAME="cda_read_coef_matrices")
    INTEGER(C_INT), VALUE :: p
    CALL ReadCoefficientMatrices(p)
  END SUBROUTINE cda_read_coef_matrices
  SUBROUTINE cda_read_perm_vector(p) BIND(C, NAME="cda_read_perm_vector")
    INTEGER(C_INT), VALUE :: p
    CALL ReadPermutationVector(p)
  END SUBROUTINE cda_read_perm_vector
  SUBROUTINE cda_read_sources(p) BIND(C, NAME="cda_read_sources")
    INTEGER(C_INT), VALUE :: p
    CALL ReadCircuitSources(p)
  END SUBROUTINE cda_read_sources
  SUBROUTINE cda_write_coeff_vectors(p) BIND(C, NAME="cda_write_coeff_vectors")
    INTEGER(C_INT), VALUE :: p
    CALL WriteCoeffVectorsForCircVariables(p)
  END SUBROUTINE cda_write_coeff_vectors
  SUBROUTINE cda_check_component_variables() BIND(C, NAME="cda_check_component_variables")
    CALL CheckComponentVariables()
  END SUBROUTINE cda_check_component_variables
  SUBROUTINE cda_check_circuit_sources() BIND(C, NAME="cda_check_circuit_sources")
    CALL CheckCircuitSources()
  END SUBROUTINE cda_check_circuit_sources
  SUBROUTINE cda_build_elem_lists() BIND(C, NAME="cda_build_elem_lists")
    CALL BuildComponentElementLists()
  END SUBROUTINE cda_build_elem_lists
  SUBROUTINE cda_circuits_summary() BIND(C, NAME="cda_circuits_summary")
    CALL CircuitsSummary()
  END SUBROUTINE cda_circuits_summary
  SUBROUTINE cda_check_transient_components() BIND(C, NAME="cda_check_transient_components")
    CALL CheckTransientComponents()
  END SUBROUTINE cda_check_transient_components
  SUBROUTINE cda_circuits_matrix_init() BIND(C, NAME="cda_circuits_matrix_init")
    CALL Circuits_MatrixInit()
  END SUBROUTINE cda_circuits_matrix_init
  SUBROUTINE cda_realloc_crt() BIND(C, NAME="cda_realloc_crt")
    IF (ALLOCATED(Ckt % Crt)) DEALLOCATE(Ckt % Crt)
    ALLOCATE(Ckt % Crt(Ckt % Circuit_tot_n))
  END SUBROUTINE cda_realloc_crt
  FUNCTION cda_info_active(level) RESULT(r) BIND(C, NAME="cda_info_active")
    INTEGER(C_INT), VALUE :: level
    INTEGER(C_INT) :: r
    r = MERGE(1, 0, InfoActive(level))
  END FUNCTION cda_info_active
  ! -- timestep block -------------------------------------------------
  SUBROUTINE cda_set_dynamic_angle() BIND(C, NAME="cda_set_dynamic_angle")
    CALL SetDynamicAngle()
  END SUBROUTINE cda_set_dynamic_angle
  FUNCTION cda_tstep_new() RESULT(r) BIND(C, NAME="cda_tstep_new")
    INTEGER(C_INT) :: r
    IF (Ckt % Tstep /= GetTimestep()) THEN
      Ckt % Tstep = GetTimestep()
      r = 1
    ELSE
      r = 0
    END IF
  END FUNCTION cda_tstep_new
  FUNCTION cda_find_lagrange() RESULT(r) BIND(C, NAME="cda_find_lagrange")
    INTEGER(C_INT) :: r
    LagrangeVar => VariableGet(Solver % Mesh % Variables, Ckt % MultName)
    r = MERGE(1, 0, ASSOCIATED(LagrangeVar))
  END FUNCTION cda_find_lagrange
  FUNCTION cda_lag_size() RESULT(n) BIND(C, NAME="cda_lag_size")
    INTEGER(C_INT) :: n
    n = SIZE(LagrangeVar % Values)
  END FUNCTION cda_lag_size
  SUBROUTINE cda_lag_range() BIND(C, NAME="cda_lag_range")
    CALL VectorValuesRange(LagrangeVar % Values, SIZE(LagrangeVar % Values), &
        TRIM(LagrangeVar % Name))
  END SUBROUTINE cda_lag_range
  FUNCTION cda_circuit_tot_n() RESULT(n) BIND(C, NAME="cda_circuit_tot_n")
    INTEGER(C_INT) :: n
    n = Ckt % Circuit_tot_n
  END FUNCTION cda_circuit_tot_n
  SUBROUTINE cda_fatal_lagrange_small() BIND(C, NAME="cda_fatal_lagrange_small")
    CALL Fatal('CircuitsAndDynamics','Lagrange multiplier is too small for circuits!')
  END SUBROUTINE cda_fatal_lagrange_small
  FUNCTION cda_lag_is_solver_var() RESULT(r) BIND(C, NAME="cda_lag_is_solver_var")
    INTEGER(C_INT) :: r
    r = MERGE(1, 0, ASSOCIATED(LagrangeVar, Solver % Variable))
  END FUNCTION cda_lag_is_solver_var
  SUBROUTINE cda_info_assoc_lag() BIND(C, NAME="cda_info_assoc_lag")
    CALL Info('CircuitsAndDynamics','Associating circuit variable to Lagrange values!',Level=8)
  END SUBROUTINE cda_info_assoc_lag
  SUBROUTINE cda_assoc_solver_var() BIND(C, NAME="cda_assoc_solver_var")
    Solver % Variable => LagrangeVar
  END SUBROUTINE cda_assoc_solver_var
  FUNCTION cda_lag_has_prev() RESULT(r) BIND(C, NAME="cda_lag_has_prev")
    INTEGER(C_INT) :: r
    r = MERGE(1, 0, ASSOCIATED(LagrangeVar % PrevValues))
  END FUNCTION cda_lag_has_prev
  SUBROUTINE cda_info_prev() BIND(C, NAME="cda_info_prev")
    CALL Info('CircuitsAndDynamics','Add PrevValues to Lagrange multiplier!',Level=8)
  END SUBROUTINE cda_info_prev
  SUBROUTINE cda_alloc_prev() BIND(C, NAME="cda_alloc_prev")
    ALLOCATE(LagrangeVar % PrevValues(SIZE(LagrangeVar % Values),1))
  END SUBROUTINE cda_alloc_prev
  SUBROUTINE cda_rotate_lagrange() BIND(C, NAME="cda_rotate_lagrange")
    LagrangeVar % PrevValues(:,1) = LagrangeVar % Values
  END SUBROUTINE cda_rotate_lagrange
  SUBROUTINE cda_load_crt() BIND(C, NAME="cda_load_crt")
    Ckt % Crt = LagrangeVar % PrevValues(1:Ckt % Circuit_tot_n,1)
  END SUBROUTINE cda_load_crt
  SUBROUTINE cda_to_mesh_var() BIND(C, NAME="cda_to_mesh_var")
    CALL Circuits_ToMeshVariable(Solver, Ckt % Crt)
  END SUBROUTINE cda_to_mesh_var
  ! -- assembly tail --------------------------------------------------
  FUNCTION cda_max_elem_dofs() RESULT(n) BIND(C, NAME="cda_max_elem_dofs")
    INTEGER(C_INT) :: n
    n = Model % Mesh % MaxElementDOFs
  END FUNCTION cda_max_elem_dofs
  FUNCTION cda_cm_assoc() RESULT(r) BIND(C, NAME="cda_cm_assoc")
    INTEGER(C_INT) :: r
    r = MERGE(1, 0, ASSOCIATED(Ckt % CircuitMatrix))
  END FUNCTION cda_cm_assoc
  SUBROUTINE cda_td_zero_cm() BIND(C, NAME="cda_td_zero_cm")
    TYPE(Matrix_t), POINTER :: CM
    CM => Ckt % CircuitMatrix
    CM % RHS = 0._dp
    IF (ASSOCIATED(CM % Values)) CM % Values = 0._dp
  END SUBROUTINE cda_td_zero_cm
  SUBROUTINE cda_td_basic(p) BIND(C, NAME="cda_td_basic")
    INTEGER(C_INT), VALUE :: p
    CALL td_AddBasicCircuitEquations(p, Ckt % Crt, dt)
  END SUBROUTINE cda_td_basic
  SUBROUTINE cda_td_comp(p, nn) BIND(C, NAME="cda_td_comp")
    INTEGER(C_INT), VALUE :: p, nn
    CALL td_AddComponentEquationsAndCouplings(p, nn, dt, Ckt % Crt)
  END SUBROUTINE cda_td_comp
  SUBROUTINE cda_link_addmatrix() BIND(C, NAME="cda_link_addmatrix")
    TYPE(Matrix_t), POINTER :: CM
    CM => Ckt % CircuitMatrix
    Ckt % ASolver % Matrix % AddMatrix => CM
    IF (CM % FORMAT == MATRIX_LIST) CALL List_toCRSMatrix(CM)
  END SUBROUTINE cda_link_addmatrix
  SUBROUTINE cda_unlink_addmatrix() BIND(C, NAME="cda_unlink_addmatrix")
    Ckt % ASolver % Matrix % AddMatrix => NULL()
  END SUBROUTINE cda_unlink_addmatrix
  SUBROUTINE cda_cm_range() BIND(C, NAME="cda_cm_range")
    TYPE(Matrix_t), POINTER :: CM
    REAL(KIND=dp), POINTER :: px(:)
    CM => Ckt % CircuitMatrix
    px => CM % Values
    CALL VectorValuesRange(px,SIZE(px),'CircuitMatrix',.TRUE.)
    px => CM % rhs
    CALL VectorValuesRange(px,SIZE(px),'CircuitRhs',.TRUE.)
  END SUBROUTINE cda_cm_range
  SUBROUTINE cda_save_circuit_matrix() BIND(C, NAME="cda_save_circuit_matrix")
    TYPE(Matrix_t), POINTER :: CM
    LOGICAL :: Found
    CM => Ckt % CircuitMatrix
    IF (ListGetLogical(Solver % Values,'Save Circuit Matrix',Found)) THEN
      CALL SaveLinearSystem(Solver, CM,'circuit',Ckt % ASolver % Matrix % NumberOfRows)
    END IF
  END SUBROUTINE cda_save_circuit_matrix
  FUNCTION cda_free_if_empty() RESULT(r) BIND(C, NAME="cda_free_if_empty")
    INTEGER(C_INT) :: r
    TYPE(Matrix_t), POINTER :: CM
    CM => Ckt % CircuitMatrix
    IF (CM % NumberOfRows <= 0) THEN
      CALL FreeMatrix(CM)
      Ckt % ASolver % Matrix % AddMatrix => NULL()
      r = 1
    ELSE
      r = 0
    END IF
  END FUNCTION cda_free_if_empty
  SUBROUTINE cda_info_finished() BIND(C, NAME="cda_info_finished")
    IF (Ckt % Harmonic) THEN
      CALL Info('CircuitsAndDynamicsHarmonic','Finished assembly of circuit matrix',Level=12)
    ELSE
      CALL Info('CircuitsAndDynamics','Finished assembly of circuit matrix',Level=12)
    END IF
  END SUBROUTINE cda_info_finished
  ! -- harmonic extras ------------------------------------------------
  SUBROUTINE cda_hm_start() BIND(C, NAME="cda_hm_start")
    CALL DefaultStart()
  END SUBROUTINE cda_hm_start
  SUBROUTINE cda_read_eigen() BIND(C, NAME="cda_read_eigen")
    LOGICAL :: Found
    EigenSystem = GetLogical(Ckt % ASolver % Values, 'Eigen Analysis', Found)
  END SUBROUTINE cda_read_eigen
  SUBROUTINE cda_hm_zero_cm() BIND(C, NAME="cda_hm_zero_cm")
    TYPE(Matrix_t), POINTER :: CM
    CM => Ckt % CircuitMatrix
    CM % RHS = 0._dp
    CM % Values = 0._dp
    IF (ASSOCIATED(CM % MassValues)) CM % MassValues = 0._dp
  END SUBROUTINE cda_hm_zero_cm
  SUBROUTINE cda_hm_basic(p) BIND(C, NAME="cda_hm_basic")
    INTEGER(C_INT), VALUE :: p
    CALL hm_AddBasicCircuitEquations(p)
  END SUBROUTINE cda_hm_basic
  SUBROUTINE cda_hm_comp(p, nn) BIND(C, NAME="cda_hm_comp")
    INTEGER(C_INT), VALUE :: p, nn
    CALL hm_AddComponentEquationsAndCouplings(p, nn)
  END SUBROUTINE cda_hm_comp
  SUBROUTINE cda_hm_finish() BIND(C, NAME="cda_hm_finish")
    CALL DefaultFinish()
  END SUBROUTINE cda_hm_finish
  ! -- output driver pieces -------------------------------------------
  SUBROUTINE cda_out_start() BIND(C, NAME="cda_out_start")
    CALL DefaultStart()
  END SUBROUTINE cda_out_start
  SUBROUTINE cda_out_claim() BIND(C, NAME="cda_out_claim")
    Ckt => ResolveCircuitModel(Solver,'CircuitsOutput')
    CALL SetCircuitModel(Ckt)
  END SUBROUTINE cda_out_claim
  FUNCTION cda_out_asolver_assoc() RESULT(r) BIND(C, NAME="cda_out_asolver_assoc")
    INTEGER(C_INT) :: r
    r = MERGE(1, 0, ASSOCIATED(Ckt % ASolver))
  END FUNCTION cda_out_asolver_assoc
  SUBROUTINE cda_fatal_no_asolver() BIND(C, NAME="cda_fatal_no_asolver")
    CALL Fatal('CircuitsOutput','ASolver not found!')
  END SUBROUTINE cda_fatal_no_asolver
  FUNCTION cda_out_ncirc() RESULT(n) BIND(C, NAME="cda_out_ncirc")
    INTEGER(C_INT) :: n
    n = Ckt % n_Circuits
  END FUNCTION cda_out_ncirc
  SUBROUTINE cda_out_mesh_var() BIND(C, NAME="cda_out_mesh_var")
    CALL Circuits_ToMeshVariable(Solver, out_crt)
  END SUBROUTINE cda_out_mesh_var
  SUBROUTINE cda_out_finish() BIND(C, NAME="cda_out_finish")
    CALL DefaultFinish()
  END SUBROUTINE cda_out_finish
!------------------------------------------------------------------------------




    
    
!------------------------------------------------------------------------------
   SUBROUTINE td_AddBasicCircuitEquations(p,Crt,dt)
!------------------------------------------------------------------------------
     IMPLICIT NONE
    INTEGER :: p
    REAL(KIND=dp) :: Crt(:)    
    REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------     
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(CircuitVariable_t), POINTER :: Cvar
    TYPE(ValueList_t), POINTER :: BF, Params
    TYPE(Matrix_t), POINTER :: CM
    INTEGER :: i, nm, RowId, ColId, j
    REAL(KIND=dp) :: vphi
    LOGICAL :: Found
    
    Circuit => CurrentModel % CircuitModel % Circuits(p)
    nm = CurrentModel % CircuitModel % ASolver % Matrix % NumberOfRows
    BF => CurrentModel % BodyForces(1) % Values
    CM => CurrentModel % CircuitModel % CircuitMatrix
  
    Params => GetSolverParams()

    DO i=1,Circuit % n
      Cvar => Circuit % CircuitVariables(i)

      IF( Ckt % Parallel ) THEN
        IF(Cvar % Owner /= ParEnv % myPE) CYCLE
      END IF
        
      RowId = Cvar % ValueId + nm
      
      IF( LEN_TRIM( Circuit % Source(i) ) > 0 ) THEN
        vphi = GetCReal(Params, Circuit % Source(i), Found)
        IF ( .NOT. Found .AND. ASSOCIATED(BF) ) THEN
          vphi = GetCReal(BF, Circuit % Source(i), Found)
        END IF
        IF (Found) Cvar % SourceRe(i) = vphi
      ELSE
        vphi = 0.0_dp
      END IF
      
      Cvar % SourceRe(i) = vphi
      CM % RHS(RowId) = Cvar % SourceRe(i)
        
      DO j=1,Circuit % n

        ColId = Circuit % CircuitVariables(j) % ValueId + nm

        IF ( TransientSimulation ) THEN 
          ! A d/dt(x): (x could be voltage or current):
          !--------------------------------------------
          IF(Cvar % A(j) /= 0._dp) THEN
            CALL AddToMatrixElement(CM, RowId, ColId, Cvar % A(j)/dt)
            CM % RHS(RowId) = CM % RHS(RowId) + Cvar % A(j) * Crt(ColId-nm) / dt
          END IF
        END IF  
        ! B x:
        ! ------
        IF(Cvar % B(j) /= 0._dp) THEN
          CALL AddToMatrixElement(CM, RowId, ColId, Cvar % B(j))
        END IF
      END DO
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE td_AddBasicCircuitEquations
!------------------------------------------------------------------------------


   
!------------------------------------------------------------------------------
   SUBROUTINE td_AddComponentEquationsAndCouplings(p, nn, dt, crt)
!------------------------------------------------------------------------------
    USE MGDynMaterialUtils
    IMPLICIT NONE
    INTEGER :: p, CompInd, nm, nn, nd, qi
    TYPE(Solver_t), POINTER :: ASolver
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(Matrix_t), POINTER :: CM
    TYPE(Component_t), POINTER :: Comp
    TYPE(CircuitVariable_t), POINTER :: Cvar
    TYPE(Valuelist_t), POINTER :: CompParams
    TYPE(Element_t), POINTER :: Element
    INTEGER :: VvarId, IvarId, q, j, astat
    REAL(KIND=dp), ALLOCATABLE :: Tcoef(:,:,:)
    REAL(KIND=dp) :: RotM(3,3,nn)
    REAL(KIND=dp) :: val, dt, crt(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilType
    LOGICAL :: Found, IsActive

    ASolver => CurrentModel % CircuitModel % ASolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('AddComponentEquationsAndCouplings','ASolver not found!')

    Circuit => CurrentModel % CircuitModel % Circuits(p)
    
    nm = Asolver % Matrix % NumberOfRows
    CM => CurrentModel % CircuitModel % CircuitMatrix

    ALLOCATE(Tcoef(3,3,nn), STAT=astat)
    IF (astat /= 0) THEN
      CALL Fatal('AddComponentEquationsAndCouplings','Memory allocation failed!')
    END IF
    
    DO CompInd = 1, Circuit % n_comp
      Comp => Circuit % Components(CompInd)

      Comp % Resistance = 0._dp 
      Comp % Conductance = 0._dp 

      Cvar => Comp % vvar
      vvarId = Comp % vvar % ValueId + nm
      IvarId = Comp % ivar % ValueId + nm

      CompParams => CurrentModel % Components(Comp % ComponentId) % Values
      IF (.NOT. ASSOCIATED(CompParams)) CALL Fatal ('AddComponentEquationsAndCouplings', 'Component parameters not found')
      IF (Comp % CoilType == 'stranded' .OR. Comp % ComponentType == 'resistor') THEN
        Comp % Resistance = ListGetCReal(CompParams, 'Resistance', Found)
        IF (Found) THEN
          Comp % UseCoilResistance = .TRUE.
        ELSE
          Comp % UseCoilResistance = .FALSE.
        END IF
      END IF

      IsActive = .TRUE.
      IF( Circuit % Parallel ) THEN
        IsActive = ( Cvar % Owner == ParEnv % myPE )
      END IF
      
      IF ( IsActive ) THEN
        IF (Comp % ComponentType == 'resistor') THEN
            CALL Info('AddComponentEquationsAndCouplings',&
                'Writing resistor equation, component '//i2s(CompInd), Level = 7)
            CALL AddToMatrixElement(CM, VvarId, IvarId, Comp % Resistance)
            CALL AddToMatrixElement(CM, VvarId, VvarId, -1._dp)
        ELSE
          SELECT CASE (Comp % CoilType)
          CASE('stranded')
            IF (Comp % UseCoilResistance) THEN
              CALL Info('AddComponentEquationsAndCouplings',&
                  'Using coil resistance for component '//i2s(CompInd), Level = 7)
              CALL AddToMatrixElement(CM, VvarId, IvarId, Comp % Resistance)
            ELSE
              Comp % Resistance = 0._dp
            END IF
            CALL AddToMatrixElement(CM, VvarId, VvarId, -1._dp)
          CASE('massive')
            CALL AddToMatrixElement(CM, VvarId, IvarId, -1._dp)
          CASE('foil winding')
            ! Foil Winding voltage: 
            ! V + ...added next... = 0
            ! ----------------------
            CALL AddToMatrixElement(CM, VvarId, VvarId, 1._dp)
            
            DO j = 1, Cvar % pdofs 
              ! Foil Winding voltage: 
              !  ... - Nf/Lalpha * int_0^{Lalpha}(V_0+V_1*alpha+V_2*alpha**2+...) = 0
              !          => ... - Nf * (V_0*Lalpha^0 + V_1/2*Lalpha^1 + V_2/3*Lalpha^2 + ...) = 0
              ! where V_m is the mth dof of the polynomial
              ! --------------------------------------------------------------
              val = - REAL(Comp % nofturns) / REAL(j) * Comp % coilthickness**(j-1)
              CALL AddToMatrixElement(CM, VvarId, j + VvarId, val)

              ! Circuit eqns for the pdofs:
              ! - Nf/Lalpha * I * int_0^1(Vi'(alpha)) + ...added later... = 0
              ! ----------------------------------------------------------
              CALL AddToMatrixElement(CM, j + VvarId, IvarId, val)
            END DO
          END SELECT
        END IF
      END IF
      
      IF (Comp % ComponentType == 'resistor') CYCLE

      ! Walked backwards over the component's own element list, which
      ! BuildComponentElementLists() filled using ElAssocToComp(). Same elements
      ! in the same order as the old "every element, test each one" loop, so the
      ! accumulation order is unchanged.
      DO qi=SIZE(Comp % ElemIdx),1,-1
        q = Comp % ElemIdx(qi)
        Element => GetActiveElement(q)
        CompParams => GetComponentParams( Element )
        IF (.NOT. ASSOCIATED(CompParams)) &
            CALL Fatal ('AddComponentEquationsAndCouplings', 'Component parameters not found')

        CoilType = ListGetString(CompParams, 'Coil Type', UnfoundFatal=.TRUE.)
        
        nn = GetElementNOFNodes(Element)
        nd = GetElementNOFDOFs(Element,ASolver)
        
        IF (SIZE(Tcoef,3) /= nn) THEN
          DEALLOCATE(Tcoef)
          ALLOCATE(Tcoef(3,3,nn), STAT=astat)
          IF ( astat /= 0 ) THEN
            CALL Fatal('AddComponentEquationsAndCouplings', 'Memory allocation error!' )
          END IF
        END IF
        
        Tcoef = GetElectricConductivityTensor(Element, nn, 're', .TRUE., CoilType)
        SELECT CASE(CoilType)
        CASE ('stranded')
          CALL td_Add_stranded(Element,Tcoef,Comp,nn,nd,dt,CompParams)
        CASE ('massive')
          IF (.NOT. HasSupport(Element,nn)) CYCLE
          CALL td_Add_massive(Element,Tcoef,Comp,nn,nd,dt,crt)
        CASE ('foil winding')
          IF (.NOT. HasSupport(Element,nn)) CYCLE
          CALL td_Add_foil_winding(Element,Tcoef,Comp,nn,nd,dt)
        CASE DEFAULT
          CALL Fatal ('AddComponentEquationsAndCouplings', 'Non-existent Coil Type Chosen!')
        END SELECT
      END DO
    END DO

    IF( CircuitsPartitionedMesh() ) THEN
      DO CompInd = 1, Circuit % n_comp
        Comp => Circuit % Components(CompInd)
        ! Only an integrated resistance is a partial sum over the partitions. One
        ! that came from the "Resistance" keyword is already the whole value on
        ! every partition, and summing it would report it PEs times too large.
        ! UseCoilResistance is derived from a keyword, so this test is uniform
        ! over the partitions and the collective below stays collective.
        IF( .NOT. Comp % UseCoilResistance ) THEN
          Comp % Resistance = ParallelReduction(Comp % Resistance)
        END IF
        Comp % Conductance = ParallelReduction(Comp % Conductance)
        ! Note: guarded by the partitioned-mesh test at the IF above, not by
        ! Solver % Parallel - both of these are integrals over this partition's
        ! elements. See CircuitsPartitionedMesh().
      END DO
    END IF
      
    DEALLOCATE(Tcoef)
!------------------------------------------------------------------------------
   END SUBROUTINE td_AddComponentEquationsAndCouplings
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE td_Add_stranded(Element,Tcoef,Comp,nn,nd,dt,CompParams)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: Tcoef(3,3,nn),dt
    TYPE(Valuelist_t), POINTER :: CompParams
    TYPE(Component_t) :: Comp
    INTEGER :: nn, nd, nm

    INTEGER :: Indexes(nd),VvarId,IvarId
    TYPE(Solver_t), POINTER :: ASolver
    INTEGER, POINTER :: PS(:)
    TYPE(Matrix_t), POINTER :: CM
    TYPE(Nodes_t), SAVE :: Nodes
    REAL(KIND=dp) :: Basis(nd), DetJ, x,POT(nd),pPOT(nd),ppPOT(nd),tscl
    REAL(KIND=dp) :: dBasisdx(nd,3), wBase(nn), w(3)
    REAL(KIND=dp) :: localC, val, circ_eq_coeff, localR !, localL
    INTEGER :: i, j,t
    REAL(KIND=dp) :: ModelDepth
    LOGICAL :: stat

    TYPE(GaussIntegrationPoints_t) :: IP
    INTEGER :: MyGen = -1
    LOGICAL :: CSymmetry, InitHandle=.TRUE., &
               CoilUseWvec=.FALSE., CoilUseWvec0=.FALSE.,Found,Found2
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilWVecVarname, CoilType

    REAL(KIND=dp) :: WBasis(nd,3), RotWBasis(nd,3)
    INTEGER :: dim, ncdofs, q, EdgeBasisDegree
    TYPE(VariableHandle_t), SAVE :: Wvec_h
    TYPE(Variable_t), POINTER, SAVE :: Wpot
    
    LOGICAL :: PiolaVersion = .FALSE.
    
    SAVE CSymmetry, dim, MyGen, InitHandle, EdgeBasisDegree

    ASolver => CurrentModel % CircuitModel % ASolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('Add_stranded','ASolver not found!')

    IF (MyGen /= CircuitsGeneration()) THEN
      MyGen = CircuitsGeneration()
      CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )
      dim = CoordinateSystemDimension()

      ! There has been somewhat different philosophies in how to create the scalar and vector fields
      ! that span the current densities. This is an effort to enable the components to all use
      ! the current density computed by the CoilSolver without writing any additional keywords to the
      ! component sections. The idea is that the circuit then only has current densities defined
      ! by the CoilSolver. If this is not desired then also no such keywords should be used in the Solver
      ! section of this module. 
      !------------------------------------------------------------------------------------------------
      CoilUseWvec0 = GetLogical(CurrentModel % Solver % Values, 'Coil Use W Vector', Found2 ) 
      DO i=1,Model % NumberOfComponents
        CoilType = ListGetString(CurrentModel % Components(i) % Values, 'Coil Type',Found)
        IF(.NOT. Found) CYCLE
        IF(CoilType == 'stranded') THEN  ! massive, foil winding
          CoilWVecVarName = GetString(CurrentModel % Components(i) % Values,'W Vector Variable Name', Found)
          IF(Found) EXIT
        END IF
      END DO
      IF(.NOT. Found) THEN
        CoilWVecVarName = GetString(CurrentModel % Solver % Values,'W Vector Variable Name', Found)
        IF(.NOT. Found) THEN
          IF( GetLogical(CurrentModel % Solver % Values,'Use Nodal CoilCurrent',Found ) ) &
              CoilWVecVarname = 'CoilCurrent'
        END IF
         IF(.NOT. Found) THEN
          IF( GetLogical(CurrentModel % Solver % Values,'Use Elemental CoilCurrent',Found ) ) &
              CoilWVecVarname = 'CoilCurrent e'
        END IF
        IF(Found) CALL Info('Add_stranded','Setting coil current to: '//TRIM(CoilWVecVarname),Level=6)
        ! If we did not find w vector named in any component it is fair to assume that it is globally used!
        IF(.NOT. Found2) CoilUseWvec0 = Found 
      END IF
      IF(.NOT. Found) CoilWVecVarname = 'W Vector E'
      CALL ListInitElementVariable(Wvec_h, CoilWVecVarname)

      CALL EdgeElementStyle(ASolver % Values, PiolaVersion, BasisDegree = EdgeBasisDegree )
      ! Only in 3D: every use of it below is inside a dim == 3 test, and a 2D
      ! case is not expected to have a wire direction potential at all.
      Wpot => NULL()
      IF (dim == 3) CALL GetWPotentialVar(Wpot)
    END IF

    PS => Asolver % Variable % Perm

    CM => CurrentModel % CircuitModel % CircuitMatrix
    nm = CurrentModel % CircuitModel % ASolver % Matrix % NumberOfRows
    
    CALL GetElementNodes(Nodes)
    nd = GetElementDOFs(Indexes,Element,ASolver)

    CALL GetLocalSolution(pPOT,UElement=Element,USolver=ASolver,tstep=-1)

    IF(Solver % Order<2.OR.GetTimeStep()<=2) THEN 
      tscl=1.0_dp
    ELSE
      tscl=1.5_dp
      CALL GetLocalSolution(ppPOT,UElement=Element,USolver=ASolver,tstep=-2)
      pPot = 2*pPOT - 0.5_dp*ppPOT
    END IF
    
    ncdofs=nd
    IF (dim == 3) THEN
      ! We can choose the base per component.
      CoilUseWvec = GetLogical(CompParams, 'Coil Use W Vector', Found)
      IF (.NOT. Found) CoilUseWvec = CoilUseWvec0
    
      IF (.NOT. CoilUseWvec) THEN
        !CALL GetLocalSolution(Wbase, 'w')
        ! when W Potential solver is used, 'w' is not enough.
        CALL GetLocalSolution( Wbase,UElement=Element,UVariable=Wpot, Found=Found)
      END IF

      ncdofs=nd-nn
    END IF

    VvarId = Comp % vvar % ValueId + nm
    IvarId = Comp % ivar % ValueId + nm

    ! Numerical integration:
    ! ----------------------
    IF(PiolaVersion) THEN
      IP = GaussPoints(Element, PReferenceElement=PiolaVersion, EdgeBasisDegree=EdgeBasisDegree)
    ELSE
      IP = GaussPoints(Element)
    END IF
    ! Model depth is constant over the element, so read it once here. It used to
    ! be read at every integration point, and each call does two keyword lookups
    ! plus CurrentCoordinateSystem().
    ModelDepth = 1.0_dp
    IF( dim == 2 ) ModelDepth = GetCircuitModelDepth()

    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      
      circ_eq_coeff = 1._dp
      SELECT CASE(dim)
      CASE(2)
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                  IP % W(t), detJ, Basis,dBasisdx )

        w = [0._dp, 0._dp, 1._dp]
        IF( CSymmetry ) THEN
          x = SUM( Basis(1:nn) * Nodes % x(1:nn) )
          detJ = detJ * x
        END IF
        circ_eq_coeff = ModelDepth
      CASE(3)

        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), IP % W(t), &
            detJ, Basis, dBasisdx, EdgeBasis = Wbasis, RotBasis = RotWBasis, USolver = ASolver ) 

        IF (CoilUseWvec) THEN
          w = ListGetElementVectorSolution( Wvec_h, Basis, Element, dofs = dim )
        ELSE
          w = -MATMUL(WBase(1:nn), dBasisdx(1:nn,:))
        END IF
      END SELECT

      localC = SUM(Tcoef(1,1,1:nn) * Basis(1:nn))
      
      IF (.NOT. Comp % UseCoilResistance) THEN
        ! I * R, where 
        ! R = (1/sigma * js,js):
        ! ----------------------
        localR = Comp % N_j **2 * IP % s(t)*detJ*SUM(w*w)/localC*circ_eq_coeff / Comp % VoltageFactor
        Comp % Resistance = Comp % Resistance + localR
      
        CALL AddToMatrixElement(CM, VvarId, IvarId, localR)
      END IF
      
      DO j=1,ncdofs
        q=j
        IF (dim == 3) q=q+nn
        IF (Comp % N_j/=0._dp) THEN
          ! ( d/dt a,w )        

          IF ( TransientSimulation ) THEN 
            IF (dim == 2) val = Comp % N_j * IP % s(t)*detJ*Basis(j)*circ_eq_coeff/dt*w(3)
            IF (dim == 3) val = Comp % N_j * IP % s(t)*detJ*SUM(WBasis(j,:)*w)/dt
            val = val / Comp % VoltageFactor

!           localL = val
!           Comp % Inductance = Comp % Inductance + localL

            CALL AddToMatrixElement(CM, VvarId, PS(Indexes(q)), tscl * val)
            CM % RHS(vvarid) = CM % RHS(vvarid) + pPOT(q) * val
         END IF
          
          ! source: 
          ! (J, rot a'), where
          ! J = w*I, thus I*(w, rot a'):
          ! ----------------------------         
          IF (dim == 2) val = -Comp % N_j*IP % s(t)*detJ*Basis(j)*w(3)
          IF (dim == 3) val = -Comp % N_j*IP % s(t)*detJ*SUM(WBasis(j,:)*w)

          val = val * Comp % SymmetryCoeff

          CALL AddToMatrixElement(CM,PS(Indexes(q)), IvarId, val)
        END IF
      END DO
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE td_Add_stranded
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE td_Add_massive(Element,Tcoef,Comp,nn,nd,dt,crt)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t) :: Element
    REAL(KIND=dp) :: Tcoef(3,3,nn),dt, crt(:)
    TYPE(Component_t) :: Comp

    TYPE(Solver_t), POINTER :: ASolver
    INTEGER, POINTER :: PS(:)
    TYPE(Matrix_t), POINTER :: CM
    REAL(KIND=dp) :: Basis(nd), DetJ, x,POT(nd),pPOT(nd),ppPOT(nd),pVel(nd),pAcc(nd),pAcc2(nd),  &
                    tscl, ascl, bscl, alpha, beta, gamma, delta, prevV, param_a, param_b
    REAL(KIND=dp) :: dBasisdx(nd,3)
    REAL(KIND=dp) :: LondonLambda(nn)
    REAL(KIND=dp) :: LondonLambda_ip, Permittivity(nn), LocalP
    REAL(KIND=dp) :: localC, val, circ_eq_coeff, grads_coeff, localConductance !, localL
    INTEGER :: nn, nd, j, t, nm, Indexes(nd), &
               VvarId, dim
    REAL(KIND=dp) :: ModelDepth
    LOGICAL :: stat, PiolaVersion, Found
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    LOGICAL :: CSymmetry
    INTEGER :: MyGen = -1
    LOGICAL :: LondonEquations
    TYPE(ValueList_t), POINTER :: Material

    REAL(KIND=dp) :: wBase(nn), gradv(3), WBasis(nd,3), RotWBasis(nd,3)
    INTEGER :: ncdofs, q, EdgeBasisDegree
    TYPE(Variable_t), POINTER, SAVE :: Wpot
    
    SAVE CSymmetry, dim, MyGen

    IF (MyGen /= CircuitsGeneration()) THEN
      MyGen = CircuitsGeneration()
      CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )
      dim = CoordinateSystemDimension()

      ! Only in 3D: every use of it below is inside a dim == 3 test, and a 2D
      ! case is not expected to have a wire direction potential at all.
      Wpot => NULL()
      IF (dim == 3) CALL GetWPotentialVar(Wpot)
    END IF

    ASolver => CurrentModel % CircuitModel % ASolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('Add_massive','ASolver not found!')
    CALL EdgeElementStyle(ASolver % Values, PiolaVersion, BasisDegree = EdgeBasisDegree)
    
    PS => Asolver % Variable % Perm

    CM => CurrentModel % CircuitModel % CircuitMatrix
    nm = CurrentModel % CircuitModel % ASolver % Matrix % NumberOfRows

    vvarId = Comp % vvar % ValueId + nm

    Material => GetMaterial( Element )

    CALL GetElementNodes(Nodes)
    nd = GetElementDOFs(Indexes,Element,ASolver)

    ! localP below is formed at every integration point whether or not the second
    ! order branch is taken, so Permittivity has to be defined either way. It used
    ! to be read only inside that branch, leaving the other path summing
    ! uninitialized stack - harmless for the result, since localP is then unused,
    ! but it traps under -ffpe-trap and shows up in any sanitizer run.
    Permittivity = 0.0_dp

    IF (ASolver % TimeOrder==2) THEN
      CALL GetLocalSolution(pPot,UElement=Element,USolver=ASolver,tstep=-3)
      CALL GetLocalSolution(pVel,UElement=Element,USolver=ASolver,tstep=-4)
      CALL GetLocalSolution(pAcc,UElement=Element,USolver=ASolver,tstep=-5)
      CALL GetLocalSolution(pAcc2,UElement=Element,USolver=ASolver,tstep=-7)

      param_a = Asolver % alpha
      param_b = Asolver % beta
      IF(param_b>=0) THEN ! this is rho_inf for Generalized-alpha
        alpha  = (2*param_b-1)/(param_b+1)
        beta = 1/(param_b+1)**2
        gamma = (3-param_b)/(2*(param_b+1))
        delta = param_b/(param_b+1)
       ELSE !bossak
        alpha = param_a
        beta = (1-alpha)**2/4
        gamma = 0.5d0 - alpha
        delta = 0
        pAcc2(1:nd) = pAcc(1:nd)
      END IF
      prevV = Crt(vvarId-nm)
      CALL GetPermittivity(Material, Permittivity, nn)
    ELSE
      CALL GetLocalSolution(pPOT,UElement=Element,USolver=ASolver,tstep=-1)
      IF(Solver % Order<2.OR.GetTimeStep()<=2) THEN 
        tscl=1.0_dp
      ELSE
        tscl=1.5_dp
        CALL GetLocalSolution(ppPOT,UElement=Element,USolver=ASolver,tstep=-2)
        pPot = (2*pPOT - 0.5_dp*ppPOT)
      END IF
    END IF

    ncdofs=nd
    IF (dim == 3) THEN
      !CALL GetLocalSolution(Wbase, 'w')      
      CALL GetLocalSolution( Wbase,UElement=Element,UVariable=Wpot, Found=Found)
      ncdofs=nd-nn
    END IF

    LondonLambda(:) = GetReal( Material, 'London Lambda', LondonEquations, Element)

    ! Numerical integration:
    ! ----------------------
    IF (PiolaVersion) THEN
      IP = GaussPoints(Element, PReferenceElement=PiolaVersion, EdgeBasisDegree=EdgeBasisDegree)
    ELSE
      IP = GaussPoints(Element)
    END IF
    ! Model depth is constant over the element, so read it once here. It used to
    ! be read at every integration point, and each call does two keyword lookups
    ! plus CurrentCoordinateSystem().
    ModelDepth = 1.0_dp
    IF( dim == 2 ) ModelDepth = GetCircuitModelDepth()

    DO t=1,IP % n
      
      grads_coeff = -1._dp
      circ_eq_coeff = 1._dp
      SELECT CASE(dim)
      CASE(2)
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis, dBasisdx )
        
        IF( CSymmetry ) THEN
          x = SUM( Basis(1:nn) * Nodes % x(1:nn) )
          detJ = detJ * x
          grads_coeff = grads_coeff/x
        END IF
        circ_eq_coeff = ModelDepth
        grads_coeff = grads_coeff/circ_eq_coeff
      CASE(3)
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), IP % W(t), &
            detJ, Basis, dBasisdx, EdgeBasis = Wbasis, RotBasis = RotWBasis, USolver = ASolver )
        gradv = MATMUL( WBase(1:nn), dBasisdx(1:nn,:))
      END SELECT

      localC = SUM(Tcoef(1,1,1:nn) * Basis(1:nn))
      localP = SUM(Permittivity(1:nn) * Basis(1:nn))

      ! computing the source term Vi(sigma grad v0, grad si):
      ! ------------------------------------------------
      IF(dim==2) val = IP % s(t)*detJ*localC*grads_coeff**2*circ_eq_coeff
      IF(dim==3) val = IP % s(t)*detJ*localC*SUM(gradv*gradv)
      val = val * Comp % VoltageFactor

      localConductance = ABS(val)
      Comp % Conductance = Comp % Conductance + localConductance

      CALL AddToMatrixElement(CM, vvarId, vvarId, val)

      IF(Asolver % TimeOrder==2) THEN
        IF(dim==2) val = IP % s(t)*detJ*localP*grads_coeff**2*circ_eq_coeff
        IF(dim==3) val = IP % s(t)*detJ*localP*SUM(gradv*gradv)

        CALL AddToMatrixElement(CM, vvarId, vvarId, val/dt)
        CM % RHS(vvarId) = CM % RHS(vvarId) + val*prevV/dt
      END IF

      IF ( LondonEquations ) THEN
        LondonLambda_ip = SUM( Basis(1:nn) * LondonLambda(1:nn) )

        val = 0.0_dp
        IF(dim==2) val = IP % s(t)*detJ/LondonLambda_ip*grads_coeff**2*circ_eq_coeff
        val = val * Comp % VoltageFactor
        ! Phi (beta grad phi_0, grad phi')
        ! Here Phi takes the place of Vi
        ! -----------------------------------
        CALL AddToMatrixElement(CM, vvarId, vvarId, val)
      END IF

      DO j=1,ncdofs
        q=j
        IF (dim == 3) q=q+nn

        IF ( LondonEquations ) THEN
          ! Phi * ( beta * grad phi, a')
          ! where phi is the node flux scalar potential
          ! -------------------------------------------
          val = 0.0_dp
          IF(dim==2) val = IP % s(t)*detJ/LondonLambda_ip*basis(j)*grads_coeff*circ_eq_coeff
          CALL AddToMatrixElement(CM, vvarId, PS(Indexes(q)), val)

          val = 0.0_dp
          IF(dim==2) val = IP % s(t)*detJ/LondonLambda_ip*basis(j)*grads_coeff
          val = val * Comp % VoltageFactor
          CALL AddToMatrixElement(CM, PS(indexes(q)), vvarId, val)
        END IF

        ! computing the mass term (sigma d/dt a, grad si):
        ! ---------------------------------------------------------
        IF ( TransientSimulation ) THEN 
          IF(dim==2) val = IP % s(t)*detJ*basis(j)*grads_coeff*circ_eq_coeff
          IF(dim==3) val = IP % s(t)*detJ*SUM(Wbasis(j,:)*gradv)

          IF(Asolver % Timeorder==1) THEN
            CALL AddToMatrixElement(CM, vvarId, PS(Indexes(q)), tscl * val * localC/dt)
            CM % RHS(vvarid) = CM % RHS(vvarid) + pPOT(q) * val * localC/dt
          ELSE
            ! 1st time derivative
            CALL AddToMatrixElement(CM, vvarId, PS(Indexes(q)), &
                  val * localC * gamma/(beta*dt) )
            CM % RHS(vvarid) = CM % RHS(vvarid) + val * localC*( gamma/(beta*dt)*pPot(q) + &
               (gamma/beta-1)*pVel(q) - ((1-gamma)+gamma*(1-1/(2*beta)))*dt*pAcc(q) )

            ! 2nd time derivative
            CALL AddToMatrixElement(CM, vvarId, PS(Indexes(q)), &
                 val * localP * (1-alpha)/(1-delta)/(beta*dt**2) )
            CM % RHS(vvarid) = CM % RHS(vvarid) + val * localP*( &
                (1-alpha)/(1-delta)/(beta*dt**2)*pPot(q) + &
                  (1-alpha)/(1-delta)/(beta*dt)*pVel(q) + &
                    delta/(1-delta)*pAcc(q) - &
                      ((1-alpha)*(1-1/(2*beta))+alpha)/(1-delta)*pAcc2(q) )

            IF(dim==2) val = IP % s(t)*detJ*localC*basis(j)*grads_coeff
            IF(dim==3) val = IP % s(t)*detJ*localC*SUM(gradv*Wbasis(j,:))
            val = val * Comp % VoltageFactor
            CALL AddToMatrixElement(CM, PS(indexes(q)), vvarId, val /  dt)
            CM % RHS(PS(indexes(q))) = CM % RHS(PS(indexes(q))) + prevV * val / dt
          END IF
        END IF

        IF(dim==2) val = IP % s(t)*detJ*localC*basis(j)*grads_coeff
        IF(dim==3) val = IP % s(t)*detJ*localC*SUM(gradv*Wbasis(j,:))
        val = val * Comp % VoltageFactor
        CALL AddToMatrixElement(CM, PS(indexes(q)), vvarId, val)
      END DO
    END DO

!------------------------------------------------------------------------------
   END SUBROUTINE td_Add_massive
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE td_Add_foil_winding(Element,Tcoef,Comp,nn,nd,dt)
!------------------------------------------------------------------------------
    USE MGDynMaterialUtils
    IMPLICIT NONE
    INTEGER :: nn, nd
    TYPE(Element_t), TARGET :: Element
    REAL(KIND=dp) :: Tcoef(3,3,nn), C(3,3), val, dt
    TYPE(Component_t) :: Comp

    TYPE(Solver_t), POINTER :: ASolver
    INTEGER, POINTER :: PS(:)
    TYPE(Matrix_t), POINTER :: CM
    REAL(KIND=dp) :: Basis(nd), DetJ, localAlpha, localV, localVtest, &
                     x, circ_eq_coeff, grads_coeff,POT(nd),pPOT(nd),ppPOT(nd),tscl
    REAL(KIND=dp) :: dBasisdx(nd,3),alpha(nn)
    REAL(KIND=dp) :: localR !, localL
    INTEGER :: nm,p,j,t,Indexes(nd),vvarId,vpolord_tot, &
               vpolord, vpolordtest, dofId, dofIdtest, &
               dim
    REAL(KIND=dp) :: ModelDepth
    LOGICAL :: stat, PiolaVersion, Found
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    LOGICAL :: CSymmetry
    INTEGER :: MyGen = -1

    REAL(KIND=dp) :: wBase(nn), gradv(3), WBasis(nd,3), RotWBasis(nd,3), &
                     RotMLoc(3,3), RotM(3,3,nn)
    INTEGER :: i,ncdofs,q,EdgeBasisDegree
    TYPE(Variable_t), POINTER, SAVE :: Wpot
    
    SAVE CSymmetry, dim, MyGen

    IF (MyGen /= CircuitsGeneration()) THEN
      MyGen = CircuitsGeneration()
      CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )
      dim = CoordinateSystemDimension()

      ! Only in 3D: every use of it below is inside a dim == 3 test, and a 2D
      ! case is not expected to have a wire direction potential at all.
      Wpot => NULL()
      IF (dim == 3) CALL GetWPotentialVar(Wpot)
    END IF

    ASolver => CurrentModel % CircuitModel % ASolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('Add_foil_winding','ASolver not found!')
    CALL EdgeElementStyle(ASolver % Values, PiolaVersion, BasisDegree = EdgeBasisDegree)
    
    PS => Asolver % Variable % Perm

    CM => CurrentModel % CircuitModel % CircuitMatrix
    nm = CurrentModel % CircuitModel % ASolver % Matrix % NumberOfRows

    CALL GetElementNodes(Nodes)
    nd = GetElementDOFs(Indexes,Element,ASolver)
    CALL GetLocalSolution(alpha,'Alpha')

    CALL GetLocalSolution(pPOT,UElement=Element,USolver=ASolver,tstep=-1)

    IF(Solver % Order<2.OR.GetTimeStep()<=2) THEN 
      tscl=1.0_dp
    ELSE
      tscl=1.5_dp
      CALL GetLocalSolution(ppPOT,UElement=Element,USolver=ASolver,tstep=-2)
      pPot = 2*pPOT - 0.5_dp*ppPOT
    END IF

    ncdofs=nd
    IF (dim == 3) THEN
      CALL GetLocalSolution( Wbase,UElement=Element,UVariable=Wpot, Found=Found)
      !CALL GetLocalSolution(Wbase, 'w')
      CALL GetElementRotM(Element, RotM, nn)
      ncdofs=nd-nn
    END IF

    vvarId = Comp % vvar % ValueId
    vpolord_tot = Comp % vvar % pdofs - 1

    ! Numerical integration:
    ! ----------------------
    IF (PiolaVersion) THEN
      IP = GaussPoints(Element, PReferenceElement=PiolaVersion, EdgeBasisDegree=EdgeBasisDegree)
    ELSE
      IP = GaussPoints(Element)
    END IF

    ! Model depth is constant over the element, so read it once here. It used to
    ! be read at every integration point, and each call does two keyword lookups
    ! plus CurrentCoordinateSystem().
    ModelDepth = 1.0_dp
    IF( dim == 2 ) ModelDepth = GetCircuitModelDepth()

    DO t=1,IP % n
      grads_coeff = -1._dp
      circ_eq_coeff = 1._dp
      SELECT CASE(dim)
      CASE(2)
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis,dBasisdx )
        IF( CSymmetry ) THEN
          x = SUM( Basis(1:nn) * Nodes % x(1:nn) )
          detJ = detJ * x
          grads_coeff = grads_coeff/x
        END IF
        circ_eq_coeff = ModelDepth
        grads_coeff = grads_coeff/circ_eq_coeff
        C(1,1) = SUM( Tcoef(3,3,1:nn) * Basis(1:nn) )
        ! I * R, where 
        ! R = (1/sigma * js,js):
        ! ----------------------
        localR = Comp % N_j **2 * IP % s(t)*detJ/C(1,1)*circ_eq_coeff/Comp % VoltageFactor
      CASE(3)
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), IP % W(t), &
            detJ, Basis, dBasisdx, EdgeBasis = Wbasis, RotBasis = RotWBasis, USolver = ASolver )
        gradv = MATMUL( WBase(1:nn), dBasisdx(1:nn,:))
        ! Compute the conductivity tensor
        ! -------------------------------
        DO i=1,3
          DO j=1,3
            C(i,j) = SUM( Tcoef(i,j,1:nn) * Basis(1:nn) )
            RotMLoc(i,j) = SUM( RotM(i,j,1:nn) * Basis(1:nn) )
          END DO
        END DO

        ! I * R, where 
        ! R = (1/sigma * js,js):
        ! ----------------------
        localR = Comp % N_j **2 * IP % s(t)*detJ/C(3,3)/Comp % VoltageFactor
        ! Transform the conductivity tensor:
        ! ----------------------------------
        C = MATMUL(MATMUL(RotMLoc, C),TRANSPOSE(RotMLoc))
      END SELECT
      
      localAlpha = SUM(alpha(1:nn) * Basis(1:nn))
      
      ! alpha is normalized to be in [0,1] thus, 
      ! it needs to be multiplied by the thickness of the coil 
      ! to get the real alpha:
      ! ------------------------------------------------------
      localAlpha = localAlpha * Comp % coilthickness

      Comp % Resistance = Comp % Resistance + localR

      DO vpolordtest=0,vpolord_tot ! V'(alpha)
        localVtest = localAlpha**vpolordtest
        dofIdtest = vpolordtest + 1 + vvarId
        DO vpolord = 0, vpolord_tot ! V(alpha)

          localV = localAlpha**vpolord
          dofId = vpolord + 1 + vvarId
          
          ! Computing the stiff term (sigma V(alpha) grad v0, V'(alpha) grad si):
          ! ---------------------------------------------------------------------
          IF (dim == 2) val = IP % s(t)*detJ*localV*localVtest*C(1,1)*grads_coeff**2*circ_eq_coeff
          IF (dim == 3) val = IP % s(t)*detJ*localV*localVtest*SUM(MATMUL(C,gradv)*gradv)
          val = val * Comp % VoltageFactor
          CALL AddToMatrixElement(CM, dofIdtest+nm, dofId+nm, val)
        END DO


        IF ( TransientSimulation ) THEN 
          DO j=1,ncdofs
            q=j
            IF (dim == 3) q=q+nn
            ! computing the mass term (sigma * d/dt * a, V'(alpha) grad si):
            ! ---------------------------------------------------------
            IF (dim == 2) val = IP % s(t)*detJ*localVtest*C(1,1)*basis(j)*grads_coeff*circ_eq_coeff/dt
            IF (dim == 3) val = IP % s(t)*detJ*localVtest*SUM(MATMUL(C,Wbasis(j,:))*gradv)/dt
  !          localL = val
  !          Comp % Inductance = Comp % Inductance + localL
            CALL AddToMatrixElement(CM, dofIdtest+nm, PS(Indexes(q)), tscl * val)
            CM % RHS(dofIdtest+nm) = CM % RHS(dofIdtest+nm) + pPOT(q) * val
          END DO
        END IF

      END DO

      DO vpolord = 0, vpolord_tot ! V(alpha)
        localV = localAlpha**vpolord
        dofId = vpolord + 1 + vvarId

        DO j=1,ncdofs
            q=j
            IF (dim == 3) q=q+nn
            IF (dim == 2) val = IP % s(t)*detJ*localV*C(1,1)*basis(j)*grads_coeff
            IF (dim == 3) val = IP % s(t)*detJ*localV*SUM(MATMUL(C,gradv)*Wbasis(j,:))
            val = val * Comp % VoltageFactor
            CALL AddToMatrixElement(CM, PS(indexes(q)), dofId+nm, val)
        END DO
      END DO

    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE td_Add_foil_winding
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Set the rotation angle in case moment of inertia and torque are given.
!------------------------------------------------------------------------------
  SUBROUTINE SetDynamicAngle()
    TYPE(Variable_t), POINTER :: AngVar, VeloVar
    TYPE(ValueList_t), POINTER :: Simulation
    ! Note: no local "dt" here. There used to be one, which shadowed the dt
    ! argument of the host routine and was never assigned, so the rotor was
    ! integrated with an undefined timestep. dt now comes by host association.
    REAL(KIND=dp) :: ang, velo, ang0, velo0, imom, torq
    INTEGER :: tStep, tStepPrev = 0
    LOGICAL :: Found
    
    SAVE ang0, velo0, imom, tStepPrev
    
    ! Variable should already exist as it was introduced in the _init section.
    AngVar => DefaultVariableGet( 'Rotor Angle' )
    IF(.NOT. ASSOCIATED( AngVar ) ) RETURN
    
    VeloVar => DefaultVariableGet( 'Rotor Velo' )
    IF(.NOT. ASSOCIATED( VeloVar ) ) THEN
      CALL Fatal('SetDynamicAngle','Variable > Rotor Velo < does not exist!')
    END IF

    ! Start from the current state so that every exit path below reports
    ! something defined.
    ang = AngVar % Values(1)
    velo = VeloVar % Values(1)

    Simulation => GetSimulation()

    IF( ListCheckPresent( Model % Simulation,'Rotor Angle') ) THEN
      CALL Info('SetDynamicAngle','Using "Rotor Velo" from simulation section',Level=6)
    ELSE      
      CALL Info('SetDynamicAngle','Using computed torque to set rotation!',Level=6)
      tStep = GetTimestep()

      imom = GetConstReal( Simulation,'Imom',Found) ! interatial moment of the motor      
      IF(.NOT. Found ) THEN
        CALL Info('SetDynamicAngle','Moment of inertia "Imom" not giving, skipping dynamics!',Level=7)
        RETURN
      END IF
      
      IF(imom < EPSILON(imom) ) THEN
        CALL Info('SetDynamicAngle','Moment of inertia "Imom" close to zero, skipping dynamics!',Level=7)
        RETURN
      END IF
      
      ! We initiate these at the start of the timestep when they still present the previous
      ! computed values. 
      IF( tStep /= tStepPrev ) THEN
        ang0 = AngVar % Values(1)
        velo0 = VeloVar % Values(1)
        tStepPrev = tStep
      END IF

      ! Without a torque the velocity is held at the value it had at the start
      ! of the timestep and only the angle keeps integrating, as the warning
      ! says. Previously this branch left velo and ang undefined and wrote them
      ! to the variables regardless.
      velo = velo0
      torq = GetConstReal( Simulation,'res: Air Gap Torque', Found)
      IF(.NOT. Found ) THEN
        CALL Warn('SetDynamicAngle','Without torque rotor velocity stays the same!')
      ELSE
        velo = velo0 + dt * torq / imom
      END IF
      ang = ang0 + dt * velo

      VeloVar % Values(1) = velo
      AngVar % Values(1) = ang
    END IF
      
    CALL ListAddConstReal(Simulation,'res: Angle(rad)', ang)
    CALL ListAddConstReal(Simulation,'res: Speed(rpm)', velo/(2._dp*pi)*60)
      
  END SUBROUTINE SetDynamicAngle
  

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------



!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE hm_AddBasicCircuitEquations(p)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(CircuitVariable_t), POINTER :: Cvar
    TYPE(ValueList_t), POINTER :: BF, Params
    TYPE(Matrix_t), POINTER :: CM
    INTEGER :: p, i, nm, RowId, ColId, j
    REAL(KIND=dp) :: Omega, vphi
    COMPLEX(KIND=dp) :: cmplx_val
    LOGICAL :: Found
    COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
    
    Circuit => CurrentModel % CircuitModel % Circuits(p)
    nm = CurrentModel % CircuitModel % ASolver % Matrix % NumberOfRows
    BF => CurrentModel % BodyForces(1) % Values
    CM => CurrentModel % CircuitModel % CircuitMatrix

    Omega = GetAngularFrequency()
    WRITE(Message,'(A,ES12.3)') 'Angular frequency for circuit equations: ',Omega
    CALL Info('CircuitsAndDynamicsHarmonic', Message, Level=6) 
    
    Params => GetSolverParams()
    
    DO i=1,Circuit % n
      Cvar => Circuit % CircuitVariables(i)

      IF(Circuit % Parallel ) THEN
        IF(Cvar % Owner /= ParEnv % myPE) CYCLE
      END IF
        
      RowId = Cvar % ValueId + nm
      
      vphi = GetCReal(Params, TRIM(Circuit % Source(i))//" re", Found)
      IF ( .NOT.Found.AND.ASSOCIATED(BF) ) THEN
        vphi = GetCReal(BF, TRIM(Circuit % Source(i))//" re", Found)
      END IF
      IF (Found) Cvar % SourceRe(i) = vphi

      !IF(Found) PRINT *,'vphi re',Found,i,vphi,TRIM(Circuit % Source(i))
       
      vphi = GetCReal(Params, TRIM(Circuit % Source(i))//" im", Found)
      IF ( .NOT.Found.AND.ASSOCIATED(BF) ) THEN
        vphi = GetCReal(BF, TRIM(Circuit % Source(i))//" im", Found)
      END IF
      IF (Found) Cvar % SourceIm(i) = vphi
      
      !IF(Found) PRINT *,'vphi im',i,vphi,TRIM(Circuit % Source(i))

      CM % RHS(RowId) = Cvar % SourceRe(i)
      CM % RHS(RowId+1) = Cvar % SourceIm(i)
        
      DO j=1,Circuit % n

        ColId = Circuit % CircuitVariables(j) % ValueId + nm

        ! im * Omega * A x: (x could be voltage or current):
        !--------------------------------------------
        IF(Cvar % A(j) /= 0._dp) THEN
          CALL AddMatrixEntry( CM, RowId, ColId, (0._dp,0._dp), im*Omega*Cvar % A(j) )
        END IF

        ! B x:
        ! ------
        IF(Cvar % B(j) /= 0._dp) THEN
          IF (Cvar % Mre(j) /= 0._dp .OR. Cvar % Mim(j) /= 0._dp) THEN
            cmplx_val = Cvar % Mre(j) + im * Cvar % Mim(j)
            cmplx_val = cmplx_val * Cvar % B(j)
          ELSE
            cmplx_val = Cvar % B(j)
          END IF
          
          CALL AddToCmplxMatrixElement(CM, RowId, ColId, REAL(cmplx_val), AIMAG(cmplx_val))
        END IF
      END DO
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE hm_AddBasicCircuitEquations
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE hm_AddComponentEquationsAndCouplings(p, nn)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: p, nn
    INTEGER :: CompInd, nm
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(Solver_t), POINTER :: ASolver
    TYPE(Matrix_t), POINTER :: CM
    TYPE(Component_t), POINTER :: Comp
    TYPE(CircuitVariable_t), POINTER :: Cvar
    TYPE(Valuelist_t), POINTER :: CompParams
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp), ALLOCATABLE :: sigma_33(:), sigmaim_33(:)
    COMPLEX(KIND=dp), ALLOCATABLE :: Tcoef(:,:,:)
    INTEGER :: VvarId, IvarId, q, j, astat, qi
    COMPLEX(KIND=dp) :: i_multiplier, cmplx_val
    COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
    REAL(KIND=dp) :: RotM(3,3,nn)
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilType
    LOGICAL :: Found

    ASolver => Ckt % ASolver
    Circuit => CurrentModel % CircuitModel % Circuits(p)
    nm = ASolver % Matrix % NumberOfRows
    CM => CurrentModel % CircuitModel % CircuitMatrix

    DO CompInd = 1, Circuit % n_comp
      Comp => Circuit % Components(CompInd)

      Comp % Resistance = 0._dp 
      Comp % Conductance = 0._dp 

      Cvar => Comp % vvar
      vvarId = Comp % vvar % ValueId + nm
      IvarId = Comp % ivar % ValueId + nm

      CompParams => CurrentModel % Components(Comp % ComponentId) % Values
      IF (.NOT. ASSOCIATED(CompParams)) CALL Fatal ('AddComponentEquationsAndCouplings', &
               'Component parameters not found')

      IF (Comp % CoilType == 'stranded' .OR. Comp % ComponentType == 'resistor') THEN
        Comp % Resistance = ListGetCReal(CompParams, 'Resistance', Found)
        IF (Found) THEN
          Comp % UseCoilResistance = .TRUE.
        ELSE
          Comp % UseCoilResistance = .FALSE.
        END IF
      END IF

      IF ( Cvar % Owner == ParEnv % myPE .OR. .NOT. Circuit % Parallel ) THEN
        IF (Comp % ComponentType == 'resistor') THEN
          ! V = R I, mirroring the transient branch. A resistor has no coil type,
          ! so the SELECT CASE below would simply skip it and leave the row empty,
          ! which made the circuit block singular.
          CALL Info('AddComponentEquationsAndCouplings',&
              'Writing resistor equation, component '//i2s(CompInd), Level = 7)
          CALL AddToCmplxMatrixElement(CM, VvarId, IvarId, Comp % Resistance, 0._dp)
          CALL AddToCmplxMatrixElement(CM, VvarId, VvarId, -1._dp, 0._dp)
        ELSE
        SELECT CASE (Comp % CoilType)
        CASE('stranded')
          IF (Comp % UseCoilResistance) THEN
            CALL Info('AddComponentEquationsAndCouplings', &
                'Using coil resistance for component '//i2s(CompInd), Level = 7)
            CALL AddToCmplxMatrixElement(CM, VvarId, IvarId, Comp % Resistance, 0._dp)
          ELSE
            Comp % Resistance = 0._dp
          END IF
 
          CALL AddToCmplxMatrixElement(CM, VvarId, VvarId, -1._dp, 0._dp)
        CASE('massive')
          i_multiplier = Comp % i_multiplier_re + im * Comp % i_multiplier_im
          IF (i_multiplier /= 0_dp) THEN
            CALL AddToCmplxMatrixElement(CM, VvarId, IvarId, -REAL(i_multiplier), -AIMAG(i_multiplier))
          ELSE
            CALL AddToCmplxMatrixElement(CM, VvarId, IvarId, -1._dp, 0._dp)
          END IF
        CASE('foil winding')
          ! Foil Winding voltage: 
          ! V + ...added next... = 0
          ! ----------------------
          i_multiplier = Comp % i_multiplier_re + im * Comp % i_multiplier_im
          CALL AddToCmplxMatrixElement(CM, VvarId, VvarId, 1._dp, 0._dp)
          
          IF (i_multiplier == 0_dp) i_multiplier = 1.0_dp
          
          DO j = 1, Cvar % pdofs 
            ! Foil Winding voltage: 
            !  ... - Nf/Lalpha * int_0^{Lalpha}(V_0+V_1*alpha+V_2*alpha**2+...) = 0
            !          => ... - Nf * (V_0*Lalpha^0 + V_1/2*Lalpha^1 + V_2/3*Lalpha^2 + ...) = 0
            ! where V_m is the mth dof of the polynomial
            ! --------------------------------------------------------------
            cmplx_val = -i_multiplier * REAL(Comp % nofturns) / REAL(j) * Comp % coilthickness**(j-1)
            CALL AddToCmplxMatrixElement(CM, VvarId, 2*j + VvarId, &
                REAL(cmplx_val), AIMAG(cmplx_val))

            ! Circuit eqns for the pdofs:
            ! - Nf/Lalpha * I * int_0^1(Vi'(alpha)) + ...added later... = 0
            ! ----------------------------------------------------------
            CALL AddToCmplxMatrixElement(CM, 2*j + VvarId, IvarId, &
               REAL(cmplx_val), AIMAG(cmplx_val))
          END DO
        END SELECT
        END IF
      END IF

      ! A resistor has no elements of its own; as in the transient path there is
      ! nothing further to assemble for it.
      IF (Comp % ComponentType == 'resistor') CYCLE

      ! See the note in the transient driver: the component's own element lists,
      ! walked backwards to preserve the previous order.
      DO qi=SIZE(Comp % ElemIdx),1,-1
        Element => GetActiveElement(Comp % ElemIdx(qi))
        CALL AddComponentElementContributions(Element, Comp, Tcoef, &
                                              sigma_33, sigmaim_33, .False.)
      END DO

      DO qi=SIZE(Comp % BCElemIdx),1,-1
        Element => GetBoundaryElement(Comp % BCElemIdx(qi))
        CALL AddComponentElementContributions(Element, Comp, Tcoef, &
                                              sigma_33, sigmaim_33, .True.)
      END DO
    END DO

    IF( CircuitsPartitionedMesh() ) THEN
      DO CompInd = 1, Circuit % n_comp
        Comp => Circuit % Components(CompInd)
        ! Only an integrated resistance is a partial sum over the partitions. One
        ! that came from the "Resistance" keyword is already the whole value on
        ! every partition, and summing it would report it PEs times too large.
        ! UseCoilResistance is derived from a keyword, so this test is uniform
        ! over the partitions and the collective below stays collective.
        IF( .NOT. Comp % UseCoilResistance ) THEN
          Comp % Resistance = ParallelReduction(Comp % Resistance)
        END IF
        Comp % Conductance = ParallelReduction(Comp % Conductance)
        ! Note: guarded by the partitioned-mesh test at the IF above, not by
        ! Solver % Parallel - both of these are integrals over this partition's
        ! elements. See CircuitsPartitionedMesh().
      END DO
    END IF
      
    IF (ALLOCATED(Tcoef)) THEN
      DEALLOCATE(Tcoef,sigma_33,sigmaim_33)
    END IF
!------------------------------------------------------------------------------
   END SUBROUTINE hm_AddComponentEquationsAndCouplings
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE AddComponentElementContributions(Element, Comp, Tcoef, & 
                                               sigma_33, sigmaim_33, boundary)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    !------------------------------------
    TYPE(Element_t), POINTER :: Element
    TYPE(Component_t), POINTER :: Comp
    COMPLEX(KIND=dp), ALLOCATABLE :: Tcoef(:,:,:)
    REAL(KIND=dp), ALLOCATABLE :: sigma_33(:), sigmaim_33(:)
    !------------------------------------
    TYPE(Solver_t), POINTER :: ASolver
    TYPE(Valuelist_t), POINTER :: CompParams
    LOGICAL :: StrandedHomogenization
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilType
    INTEGER :: nn_elem, nd_elem
    INTEGER :: astat
    LOGICAL :: Found, FoundIm, boundary

    IF (ElAssocToComp(Element, Comp)) THEN
      ASolver => CurrentModel % CircuitModel % ASolver
      IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('AddComponentEquationsAndCouplings','ASolver not found!')

      CompParams => GetComponentParams( Element )
      IF (.NOT. ASSOCIATED(CompParams)) CALL Fatal ('AddComponentElementContributions',&
                                                    'Component parameters not found')

      StrandedHomogenization = .FALSE.
      CoilType = GetString(CompParams, 'Coil Type', Found)
      IF (.NOT. Found) CoilType = ''
      
      nn_elem = GetElementNOFNodes(Element)
      nd_elem = GetElementNOFDOFs(Element,ASolver)

      IF (.NOT. ALLOCATED(Tcoef)) THEN
        ALLOCATE(Tcoef(3,3,nn_elem), sigma_33(nn_elem), sigmaim_33(nn_elem), STAT=astat)
        IF (astat /= 0) THEN
          CALL Fatal ('AddComponentEquationsAndCouplings','Memory allocation failed')
        END IF
      ELSE IF (SIZE(Tcoef,3) /= nn_elem) THEN
        DEALLOCATE(Tcoef, sigma_33, sigmaim_33)
        ALLOCATE(Tcoef(3,3,nn_elem),sigma_33(nn_elem), sigmaim_33(nn_elem), STAT=astat)
        IF (astat /= 0) THEN
          CALL Fatal ('AddComponentEquationsAndCouplings','Memory allocation failed')
        END IF
      END IF
      
      SELECT CASE(CoilType)
      CASE ('stranded')
        StrandedHomogenization = GetLogical(CompParams, 'Homogenization Model', Found)
        IF ( StrandedHomogenization ) THEN 
          sigma_33 = GetReal(CompParams, 'sigma 33', Found)
          IF ( .NOT. Found ) sigma_33 = 0._dp
          sigmaim_33 = GetReal(CompParams, 'sigma 33 im', FoundIm)
          IF ( .NOT. FoundIm ) sigmaim_33 = 0._dp
          IF ( .NOT. Found .AND. .NOT. FoundIm ) CALL Fatal ('LocalMatrix', 'Homogenization Model sigma 33 not found!')
          IF ( .NOT. Found .AND. .NOT. FoundIm ) CALL Fatal ('AddComponentEquationsAndCouplings', &
              'Homogenization Model Sigma 33 not found!')
          Tcoef = CMPLX(0._dp, 0._dp, KIND=dp)
          Tcoef(3,3,1:nn_elem) = CMPLX(sigma_33, sigmaim_33, KIND=dp)
        ELSE
          Tcoef = GetCMPLXElectricConductivityTensor(Element, nn_elem, .TRUE., CoilType) 
        END IF
        CALL hm_Add_stranded(Element,Tcoef,Comp,nn_elem,nd_elem,CompParams)
      CASE ('massive')
        IF (HasSupport(Element,nn_elem)) THEN
          Tcoef = GetCMPLXElectricConductivityTensor(Element, nn_elem, .TRUE., CoilType) 
          CALL hm_Add_massive(Element,Tcoef,Comp,nn_elem,nd_elem)
        END IF
      CASE ('foil winding')
        IF (HasSupport(Element,nn_elem)) THEN
          Tcoef = GetCMPLXElectricConductivityTensor(Element, nn_elem, .TRUE., CoilType) 
          CALL hm_Add_foil_winding(Element,Tcoef,Comp,nn_elem,nd_elem,CompParams)
        END IF
      CASE DEFAULT
        CALL Fatal ('AddComponentEquationsAndCouplings', 'Non existent Coil Type Chosen!')
      END SELECT
    END IF
!------------------------------------------------------------------------------
   END SUBROUTINE AddComponentElementContributions
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE hm_Add_stranded(Element,Tcoef,Comp,nn,nd,CompParams)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER :: Element
    COMPLEX(KIND=dp) :: Tcoef(3,3,nn)
    TYPE(Component_t) :: Comp
    TYPE(Valuelist_t), POINTER :: CompParams
    INTEGER :: nn, nd, nm, Indexes(nd),VvarId,IvarId

    TYPE(Solver_t), POINTER :: ASolver
    INTEGER, POINTER :: PS(:)
    TYPE(Matrix_t), POINTER :: CM
    REAL(KIND=dp) :: Omega
    TYPE(Nodes_t), SAVE :: Nodes
    REAL(KIND=dp) :: Basis(nd), DetJ, x, circ_eq_coeff
    REAL(KIND=dp) :: dBasisdx(nd,3), wBase(nn), w(3)
    COMPLEX(KIND=dp) :: localC, i_multiplier, cmplx_val
    REAL(KIND=dp) :: localR !, localL
    INTEGER :: i, j,t
    REAL(KIND=dp) :: ModelDepth
    LOGICAL :: stat

    TYPE(GaussIntegrationPoints_t) :: IP
    COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
    LOGICAL :: CSymmetry
    INTEGER :: MyGen = -1

    REAL(KIND=dp) :: WBasis(nd,3), RotWBasis(nd,3)
    INTEGER :: dim, ncdofs, q, EdgeBasisDegree
    
    LOGICAL :: CoilUseWvec=.FALSE., CoilUseWvec0=.FALSE.,Found,Found2
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilWVecVarname, CoilType

    LOGICAL :: PiolaVersion = .FALSE.

    TYPE(VariableHandle_t), SAVE :: Wvec_h
    TYPE(Variable_t), POINTER, SAVE :: Wpot
    
    SAVE CSymmetry, dim, MyGen

    IF (MyGen /= CircuitsGeneration()) THEN
      MyGen = CircuitsGeneration()
      CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )
      dim = CoordinateSystemDimension()

      CoilUseWvec0 = GetLogical(CurrentModel % Solver % Values, 'Coil Use W Vector', Found2 ) 
      DO i=1,Model % NumberOfComponents
        CoilType = ListGetString(CurrentModel % Components(i) % Values, 'Coil Type',Found)
        IF(.NOT. Found) CYCLE
        IF(CoilType == 'stranded') THEN  ! massive, foil winding
          CoilWVecVarName = GetString(CurrentModel % Components(i) % Values,'W Vector Variable Name', Found)
          IF(Found) EXIT
        END IF
      END DO
      IF(.NOT. Found) THEN
        CoilWVecVarName = GetString(CurrentModel % Solver % Values,'W Vector Variable Name', Found)
        IF(.NOT. Found) THEN
          IF( GetLogical(CurrentModel % Solver % Values,'Use Nodal CoilCurrent',Found ) ) &
              CoilWVecVarname = 'CoilCurrent'
        END IF
         IF(.NOT. Found) THEN
          IF( GetLogical(CurrentModel % Solver % Values,'Use Elemental CoilCurrent',Found ) ) &
              CoilWVecVarname = 'CoilCurrent e'
        END IF
        IF(Found) CALL Info('Add_stranded','Setting coil current to: '//TRIM(CoilWVecVarname),Level=6)
        ! If we did not find w vector named in any component it is fair to assume that it is globally used!
        IF(.NOT. Found2) CoilUseWvec0 = Found 
      END IF
      IF(.NOT. Found) CoilWVecVarname = 'W Vector E'
      CALL ListInitElementVariable(Wvec_h, CoilWVecVarname)

      ! Only in 3D: every use of it below is inside a dim == 3 test, and a 2D
      ! case is not expected to have a wire direction potential at all.
      Wpot => NULL()
      IF (dim == 3) CALL GetWPotentialVar(Wpot)
    END IF

    ASolver => CurrentModel % CircuitModel % ASolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('Add_stranded','ASolver not found!')
    CALL EdgeElementStyle(ASolver % Values, PiolaVersion, BasisDegree = EdgeBasisDegree )

    PS => Asolver % Variable % Perm

    CM => CurrentModel % CircuitModel % CircuitMatrix
    nm = CurrentModel % CircuitModel % ASolver % Matrix % NumberOfRows
    Omega = GetAngularFrequency()
    
    CALL GetElementNodes(Nodes)
    nd = GetElementDOFs(Indexes,Element,ASolver)
    
    ncdofs=nd
    IF (dim == 3) THEN
      ncdofs=nd-nn

      CoilUseWvec = GetLogical(CompParams, 'Coil Use W Vector', Found)
      IF(.NOT. Found) CoilUseWvec = CoilUseWvec0 

      IF (.NOT. CoilUseWvec) THEN
        CALL GetLocalSolution( Wbase,UElement=Element,UVariable=Wpot, Found=Found)
        !CALL GetWPotential(WBase)
      END IF
    END IF

    VvarId = Comp % vvar % ValueId + nm
    IvarId = Comp % ivar % ValueId + nm

    i_multiplier = Comp % i_multiplier_re + im * Comp % i_multiplier_im


    ! Numerical integration:
    ! ----------------------
    IF(PiolaVersion) THEN
      IP = GaussPoints(Element, PReferenceElement=PiolaVersion, EdgeBasisDegree=EdgeBasisDegree)
    ELSE
      IP = GaussPoints(Element)
    END IF

    ! Model depth is constant over the element, so read it once here. It used to
    ! be read at every integration point, and each call does two keyword lookups
    ! plus CurrentCoordinateSystem().
    ModelDepth = 1.0_dp
    IF( dim == 2 ) ModelDepth = GetCircuitModelDepth()

    DO t=1,IP % n
 
      circ_eq_coeff = 1._dp
      SELECT CASE(dim)
      CASE(2)
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                  IP % W(t), detJ, Basis, dBasisdx )
        w = [0._dp, 0._dp, 1._dp]
        IF( CSymmetry ) THEN
          x = SUM( Basis(1:nn) * Nodes % x(1:nn) )
          detJ = detJ * x
        END IF
        circ_eq_coeff = ModelDepth
      CASE(3)
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), IP % W(t), &
            detJ, Basis, dBasisdx, EdgeBasis = Wbasis, RotBasis = RotWBasis, USolver = ASolver)

        IF (CoilUseWvec) THEN
          w = ListGetElementVectorSolution( Wvec_h, Basis, Element, Found = Found, dofs = dim )
          IF(.NOT. Found ) THEN
            CALL Fatal('Add_stranded','Could not find coil current density!')            
          END IF
        ELSE
          w = -MATMUL(WBase(1:nn), dBasisdx(1:nn,:))
        END IF
      END SELECT

      localC = SUM(Tcoef(3,3,1:nn) * Basis(1:nn))
      
      IF (.NOT. Comp % UseCoilResistance) THEN
        ! I * R, where 
        ! R = (1/sigma * js,js):
        ! ----------------------
        localR = REAL( Comp % N_j **2 * IP % s(t)*detJ*SUM(w*w)/localC*circ_eq_coeff &
            / Comp % VoltageFactor, KIND=dp )

        Comp % Resistance = Comp % Resistance + localR
        
        CALL AddToCmplxMatrixElement(CM, VvarId, IvarId, &
              REAL(Comp % N_j**2 * IP % s(t)*detJ*SUM(w*w)/localC*circ_eq_coeff / Comp % VoltageFactor), &
             AIMAG(Comp % N_j**2 * IP % s(t)*detJ*SUM(w*w)/localC*circ_eq_coeff / Comp % VoltageFactor))
      END IF
      
      DO j=1,ncdofs
        q=j
        IF (dim == 3) q=q+nn
        IF (Comp % N_j/=0._dp) THEN
          ! ( im * Omega a,w )
          IF (dim == 2) cmplx_val = im * Omega * Comp % N_j &
                  * IP % s(t)*detJ*Basis(j)*circ_eq_coeff*w(3)
          IF (dim == 3) cmplx_val = im * Omega * Comp % N_j &
                  * IP % s(t)*detJ*SUM(WBasis(j,:)*w)

          cmplx_val = cmplx_val / Comp % VoltageFactor
!          localL = ABS(cmplx_val)
!          Comp % Inductance = Comp % Inductance + localL

          CALL AddToCmplxMatrixElement(CM, VvarId, ReIndex(PS(Indexes(q))), &
                 REAL(cmplx_val), AIMAG(cmplx_val))
          
          IF (dim == 2) cmplx_val = -Comp % N_j*IP % s(t)*detJ*Basis(j)*w(3)
          IF (dim == 3) cmplx_val = -Comp % N_j*IP % s(t)*detJ*SUM(WBasis(j,:)*w)
          IF (i_multiplier /= 0._dp) cmplx_val = i_multiplier*cmplx_val

          cmplx_val = cmplx_val * Comp % SymmetryCoeff
          
          CALL AddToCmplxMatrixElement(CM,ReIndex(PS(Indexes(q))), IvarId, &
             REAL(cmplx_val), AIMAG(cmplx_val))

        END IF
      END DO
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE hm_Add_stranded
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE hm_Add_massive(Element,Tcoef,Comp,nn,nd)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t) :: Element
    COMPLEX(KIND=dp) :: Tcoef(3,3,nn)
    TYPE(Component_t) :: Comp

    TYPE(Solver_t), POINTER :: ASolver
    TYPE(ValueList_t), POINTER :: BC
    INTEGER, POINTER :: PS(:)
    TYPE(Matrix_t), POINTER :: CM
    REAL(KIND=dp) :: Omega, grads_coeff, circ_eq_coeff
    REAL(KIND=dp) :: Basis(nd), DetJ, x
    REAL(KIND=dp) :: dBasisdx(nd,3)
    REAL(KIND=dp) :: LondonLambda(nn)
    REAL(KIND=dp) :: LondonLambda_ip, val
    REAL(KIND=dp) :: SkinCond(nn), SkinMu(nn)
    REAL(KIND=dp) :: cond, mu, muVacuum, delta
    COMPLEX(KIND=dp) :: imu, invZs
    COMPLEX(KIND=dp) :: localC, cmplx_val=0._dp
    REAL(KIND=dp) :: localConductance !, localL
    INTEGER :: nn, nd, j, t, nm, Indexes(nd), &
               VvarId, dim
    LOGICAL :: stat, PiolaVersion, Found
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
    LOGICAL :: CSymmetry
    INTEGER :: MyGen = -1
    LOGICAL :: LondonEquations, SkinBc=.False., ElectroDynamics
    TYPE(ValueList_t), POINTER :: Material

    REAL(KIND=dp) :: wBase(nn), gradv(3), WBasis(nd,3), RotWBasis(nd,3)
    INTEGER :: ncdofs,q,EdgeBasisDegree
    REAL(KIND=dp) :: ModelDepth
    COMPLEX(KIND=dp) :: Permittivity(nn), localP
    TYPE(Variable_t), POINTER, SAVE :: Wpot

    
    SAVE CSymmetry, dim, MyGen

    IF (MyGen /= CircuitsGeneration()) THEN
      MyGen = CircuitsGeneration()
      CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )
      dim = CoordinateSystemDimension()

      ! Only in 3D: every use of it below is inside a dim == 3 test, and a 2D
      ! case is not expected to have a wire direction potential at all.
      Wpot => NULL()
      IF (dim == 3) CALL GetWPotentialVar(Wpot)
    END IF

    ASolver => CurrentModel % CircuitModel % ASolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('Add_massive','ASolver not found!')
    CALL EdgeElementStyle(ASolver % Values, PiolaVersion, BasisDegree = EdgeBasisDegree)
    
    PS => Asolver % Variable % Perm

    CM => CurrentModel % CircuitModel % CircuitMatrix
    nm = CurrentModel % CircuitModel % ASolver % Matrix % NumberOfRows
    Omega = GetAngularFrequency()

    SkinCond = 0._dp
    SkinMu = 0._dp

    CALL GetElementNodes(Nodes)
    nd = GetElementDOFs(Indexes,Element,ASolver)

    ncdofs=nd
    IF (dim == 3) THEN
      !CALL GetWPotential(WBase)     
      CALL GetLocalSolution( Wbase,UElement=Element,UVariable=Wpot, Found=Found)
      ncdofs=nd-nn
    END IF

    Material => GetMaterial( Element )

    ElectroDynamics = GetLogical( Asolver % values, 'Electrodynamics Model', Found)
    IF(ElectroDynamics) THEN
      CALL GetPermittivity(Material, Permittivity, nn)
    END IF

    LondonLambda(:) = GetReal( Material, 'London Lambda', LondonEquations, Element)

    BC => GetBC( Element )
    skinBc = .FALSE.
    IF ( ASSOCIATED(BC) ) THEN
      SkinCond = GetConstReal( BC, 'Layer Electric Conductivity', SkinBc)
      IF ( SkinBc ) THEN
        muVacuum = 4 * PI * 1d-7
        imu = CMPLX(0.0_dp, 1.0_dp, KIND=dp) 
        SkinMu = GetConstReal( BC, 'Layer Relative Permeability', Found)
      END IF
    END IF
  
    vvarId = Comp % vvar % ValueId + nm

    ! Numerical integration:
    ! ----------------------
    IF (PiolaVersion) THEN
      IP = GaussPoints(Element, PReferenceElement=PiolaVersion, EdgeBasisDegree=EdgeBasisDegree)
    ELSE
      IP = GaussPoints(Element)
    END IF

    ! Model depth is constant over the element, so read it once here. It used to
    ! be read at every integration point, and each call does two keyword lookups
    ! plus CurrentCoordinateSystem().
    ModelDepth = 1.0_dp
    IF( dim == 2 ) ModelDepth = GetCircuitModelDepth()

    DO t=1,IP % n
      grads_coeff = -1._dp
      circ_eq_coeff = 1._dp
      SELECT CASE(dim)
      CASE(2)
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis,dBasisdx )

        IF( CSymmetry ) THEN
          x = SUM( Basis(1:nn) * Nodes % x(1:nn) )
          detJ = detJ * x
          grads_coeff = grads_coeff/x
        END IF
        circ_eq_coeff = ModelDepth
        grads_coeff = grads_coeff/ModelDepth
      CASE(3)
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), IP % W(t), &
            detJ, Basis, dBasisdx, EdgeBasis = Wbasis, RotBasis = RotWBasis, USolver = ASolver )
        gradv = MATMUL( WBase(1:nn), dBasisdx(1:nn,:))
      END SELECT

      localC = SUM(Tcoef(1,1,1:nn) * Basis(1:nn))

      localP = 0._dp
      IF(ElectroDynamics) THEN
        localP = SUM( Basis(1:nn) * Permittivity(1:nn))
      END IF

      ! computing the source term Vi(sigma grad v0, grad si):
      ! ------------------------------------------------
      IF(dim==2) cmplx_val = IP % s(t)*detJ*grads_coeff**2*circ_eq_coeff * Comp % VoltageFactor

      IF(dim==3) THEN
        IF (SkinBc) THEN
        ! if SkinBC is activated:
        !   Boundary Condition: Layer Electric Conductivity
        !   Boundary Condition: Layer Relative Permeability
        !
        ! The term Vi/Z ( grad v0, grad_v0 )
        !
          cond = SUM(Basis(1:nn) * SkinCond(1:nn))
          mu  = muVacuum * SUM(Basis(1:nn) * SkinMu(1:nn))
          delta = SQRT( 2.0_dp/(cond*omega*mu))      
          invZs = (cond*delta)/(1.0_dp+imu)
          cmplx_val = IP % s(t)*detJ*invZs*SUM(gradv*gradv) * Comp % VoltageFactor
        ELSE
          cmplx_val = IP % s(t)*detJ*SUM(gradv*gradv) * Comp % VoltageFactor
        END IF
      END IF

      IF(SkinBC) THEN
        CALL AddToCmplxMatrixElement(CM, vvarId, vvarId, &
             REAL(cmplx_val), AIMAG(cmplx_val))
      ELSE
        CALL AddMatrixEntry(CM, vvarId, vvarId, &
             cmplx_val*localC, im*Omega*cmplx_val*localP )
      END IF

      localConductance = ABS(cmplx_val*localC)
      Comp % Conductance = Comp % Conductance + localConductance

      IF ( LondonEquations ) THEN
        LondonLambda_ip = SUM( Basis(1:nn) * LondonLambda(1:nn) )

        val = 0.0_dp
        IF(dim==2) val = IP % s(t)*detJ/LondonLambda_ip*grads_coeff**2*circ_eq_coeff
        val = val * Comp % VoltageFactor
        ! Phi (beta grad phi_0, grad phi')
        ! Here Phi takes the place of Vi
        ! -----------------------------------
        CALL AddToCmplxMatrixElement(CM, vvarId, vvarId, val, 0._dp)
      END IF


      DO j=1,ncdofs
        q=j
        IF (dim == 3) q=q+nn

        IF ( LondonEquations ) THEN
          ! Phi * ( beta * grad phi, a')
          ! where phi is the node flux scalar potential
          ! -------------------------------------------
          val = 0.0_dp
          IF(dim==2) val = IP % s(t)*detJ/LondonLambda_ip*basis(j)*grads_coeff*circ_eq_coeff
          CALL AddToCmplxMatrixElement(CM, vvarId, ReIndex(PS(Indexes(q))), val, 0._dp)

          val = 0.0_dp
          IF(dim==2) val = IP % s(t)*detJ/LondonLambda_ip*basis(j)*grads_coeff
          val = val * Comp % VoltageFactor
          CALL AddToCmplxMatrixElement(CM, ReIndex(PS(indexes(q))), vvarId, val, 0._dp)
        END IF

        ! computing the mass term (sigma * im * Omega * a, grad si):
        ! ---------------------------------------------------------
        IF(dim==2) cmplx_val = IP % s(t)*detJ*Basis(j)*grads_coeff*circ_eq_coeff

        IF(dim==3) THEN
          IF (SkinBc) THEN
          ! if SkinBC is activated:
          !   Boundary Condition: Layer Electric Conductivity
          !   Boundary Condition: Layer Relative Permeability
          ! Then activate
          !  (1/Z*im*Omega* a , grad v')
          !
            cmplx_val = IP % s(t)*detJ*invZs*im*Omega*SUM(Wbasis(j,:)*gradv)
          ELSE
            cmplx_val = IP % s(t)*detJ*SUM(Wbasis(j,:)*gradv)
          END IF
        END IF

        IF(SkinBc ) THEN
          CALL AddToCmplxMatrixElement(CM, vvarId, ReIndex(PS(Indexes(q))), &
                    REAL(cmplx_val), AIMAG(cmplx_val))
        ELSE
          cmplx_val = im*Omega*cmplx_val*localC - omega**2*cmplx_val*localP
          CALL AddMatrixEntry(CM, vvarId, ReIndex(PS(Indexes(q))), &
                     (0._dp, 0._dp), cmplx_val )
        END IF


!        localL = ABS(1._dp/cmplx_val)
!        Comp % Inductance = Comp % Inductance + localL
        
        IF(dim==2) cmplx_val = IP % s(t)*detJ*basis(j)*grads_coeff * Comp % VoltageFactor

        IF(dim==3) THEN
          IF (SkinBc) THEN
            ! if SkinBC is activated:
            !   Layer Electric Conductivity = Real 58e6
            !   Layer Relative Permeability = Real 1
            !
            !  + 1/Z (grad v , a') 
            !
            cmplx_val = IP % s(t)*detJ*invZs*SUM(gradv*Wbasis(j,:)) * Comp % VoltageFactor
          ELSE
            cmplx_val = IP % s(t)*detJ*SUM(gradv*Wbasis(j,:)) * Comp % VoltageFactor 
          END IF
        END IF

        IF(SkinBc) THEN
          CALL AddToCmplxMatrixElement(CM, ReIndex(PS(indexes(q))), vvarId, &
                  REAL(cmplx_val), AIMAG(cmplx_val))
        ELSE
          CALL AddMatrixEntry(CM, ReIndex(PS(indexes(q))), vvarId, &
               cmplx_val*LocalC, im*Omega*cmplx_val*LocalP )
        END IF
      END DO
    END DO

!------------------------------------------------------------------------------
   END SUBROUTINE hm_Add_massive
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE hm_Add_foil_winding(Element,Tcoef,Comp,nn,nd,CompParams)
!------------------------------------------------------------------------------
    USE MGDynMaterialUtils
    IMPLICIT NONE
    INTEGER :: nn, nd
    TYPE(Element_t), POINTER :: Element
    COMPLEX(KIND=dp) :: Tcoef(3,3,nn), C(3,3), val
    TYPE(Component_t) :: Comp
    TYPE(Valuelist_t), POINTER :: CompParams

    TYPE(Solver_t), POINTER :: ASolver
    INTEGER, POINTER :: PS(:)
    TYPE(Matrix_t), POINTER :: CM
    REAL(KIND=dp) :: Basis(nd), DetJ, Omega, localAlpha, localV, localVtest, &
                     x, circ_eq_coeff, grads_coeff
    REAL(KIND=dp) :: dBasisdx(nd,3),alpha(nn)
    INTEGER :: nm,p,j,t,Indexes(nd),vvarId,vpolord_tot, &
               vpolord, vpolordtest, dofId, dofIdtest, &
               dim
    REAL(KIND=dp) :: ModelDepth
    LOGICAL :: stat, PiolaVersion
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
    INTEGER :: MyGen = -1
    LOGICAL :: CSymmetry, InitHandle=.TRUE., &
               CoilUseWvec=.FALSE., CoilUseWvec0=.FALSE.,Found,Found2
    LOGICAL :: InitJHandle=.TRUE., FoilUseJvec=.FALSE.
    REAL(KIND=dp) :: localR
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilWVecVarname, CoilType
    CHARACTER(LEN=MAX_NAME_LEN) :: FoilJVecVarname
    TYPE(VariableHandle_t), SAVE :: Wvec_h
    TYPE(VariableHandle_t), SAVE :: Jvec_h

    REAL(KIND=dp) :: wBase(nn), gradv(3), WBasis(nd,3), RotWBasis(nd,3), &
                     RotMLoc(3,3), RotM(3,3,nn)
    ! Jvec is sigma * grad v0 and sigma may be complex ("Electric Conductivity
    ! im", homogenization models). It used to be REAL, which silently threw the
    ! imaginary part away in the V-V stiffness block and in the A-V coupling
    ! block below, while the V-A block kept the full complex C.
    COMPLEX(KIND=dp) :: Jvec(3)
    INTEGER :: i,ncdofs,q,EdgeBasisDegree
    TYPE(Variable_t), POINTER, SAVE :: Wpot

    
    SAVE CSymmetry, dim, MyGen, InitHandle, InitJHandle

    IF (MyGen /= CircuitsGeneration()) THEN
      MyGen = CircuitsGeneration()
      CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )
      dim = CoordinateSystemDimension()
      
      CoilUseWvec0 = GetLogical(CurrentModel % Solver % Values, 'Coil Use W Vector', Found2 ) 
      DO i=1,Model % NumberOfComponents
        CoilType = ListGetString(CurrentModel % Components(i) % Values, 'Coil Type',Found)
        IF(.NOT. Found) CYCLE
        IF(CoilType == 'foil winding') THEN  ! stranded, massive
          CoilWVecVarName = GetString(CurrentModel % Components(i) % Values,'W Vector Variable Name', Found)
          IF(Found) EXIT
        END IF
      END DO
      IF(.NOT. Found) THEN
        CoilWVecVarName = GetString(CurrentModel % Solver % Values,'W Vector Variable Name', Found)
        IF(.NOT. Found) THEN
          IF( GetLogical(CurrentModel % Solver % Values,'Use Nodal CoilCurrent',Found ) ) &
              CoilWVecVarname = 'CoilCurrent'
        END IF
        IF(.NOT. Found) THEN
          IF( GetLogical(CurrentModel % Solver % Values,'Use Elemental CoilCurrent',Found ) ) &
              CoilWVecVarname = 'CoilCurrent e'
        END IF
        IF(Found) CALL Info('Add_foil_winding','Setting coil current to: '//TRIM(CoilWVecVarname),Level=6)
        ! If we did not find w vector named in any component it is fair to assume that it is globally used!
        IF(.NOT. Found2) CoilUseWvec0 = Found 
      END IF
      IF(.NOT. Found) CoilWVecVarname = 'W Vector E'
      CALL ListInitElementVariable(Wvec_h, CoilWVecVarname)

      ! Only in 3D: every use of it below is inside a dim == 3 test, and a 2D
      ! case is not expected to have a wire direction potential at all.
      Wpot => NULL()
      IF (dim == 3) CALL GetWPotentialVar(Wpot)
    END IF

    ASolver => CurrentModel % CircuitModel % ASolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('Add_foil_winding','ASolver not found!')
    CALL EdgeElementStyle(ASolver % Values, PiolaVersion, BasisDegree = EdgeBasisDegree)
    
    PS => Asolver % Variable % Perm

    CM => CurrentModel % CircuitModel % CircuitMatrix
    nm = CurrentModel % CircuitModel % ASolver % Matrix % NumberOfRows
    Omega = GetAngularFrequency()

    CALL GetElementNodes(Nodes)
    nd = GetElementDOFs(Indexes,Element,ASolver)
    CALL GetLocalSolution(alpha,'Alpha')
    
    ncdofs=nd
    IF (dim == 3) THEN

      ! If we do not have a local flag then use the one from the solver section
      CoilUseWvec = GetLogical(CompParams, 'Coil Use W Vector', Found)
      IF (.NOT. Found) CoilUseWvec = CoilUseWvec0

      IF (.NOT. CoilUseWvec) THEN
        !CALL GetLocalSolution(Wbase, 'w')
        !CALL GetWPotential(WBase)
        CALL GetLocalSolution( Wbase,UElement=Element,UVariable=Wpot, Found=Found)
      END IF

      FoilUseJvec = GetLogical(CompParams, 'Foil Winding Use J Vector', Found)
      IF (.NOT. Found) FoilUseJvec = .FALSE.

      IF (FoilUseJvec) THEN
        IF( InitJHandle ) THEN
          FoilJVecVarname = GetString(CompParams, 'Foil J Vector Variable Name', Found)
          IF ( .NOT. Found) FoilJVecVarname = 'J Vector E'
          CALL ListInitElementVariable(Jvec_h, FoilJVecVarname)
          IF ( .NOT. ASSOCIATED(Jvec_h % Variable)) THEN
            CALL Fatal('Add_foil_winding','You are trying to use Foil J Vector for '//&
                'describing the component source field but I cannot find the variable')
          END IF
          InitJHandle = .FALSE.
        END IF
      END IF

      CALL GetElementRotM(Element, RotM, nn)

      ncdofs=nd-nn
    END IF

    vvarId = Comp % vvar % ValueId
    vpolord_tot = Comp % vvar % pdofs - 1

    ! Numerical integration:
    ! ----------------------
    IF (PiolaVersion) THEN
      IP = GaussPoints(Element, PReferenceElement=PiolaVersion, EdgeBasisDegree=EdgeBasisDegree)
    ELSE
      IP = GaussPoints(Element)
    END IF

    ! Model depth is constant over the element, so read it once here. It used to
    ! be read at every integration point, and each call does two keyword lookups
    ! plus CurrentCoordinateSystem().
    ModelDepth = 1.0_dp
    IF( dim == 2 ) ModelDepth = GetCircuitModelDepth()

    DO t=1,IP % n
      grads_coeff = -1._dp
      circ_eq_coeff = 1._dp
      SELECT CASE(dim)
      CASE(2)
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis,dBasisdx )
        
        IF( CSymmetry ) THEN
          x = SUM( Basis(1:nn) * Nodes % x(1:nn) )
          detJ = detJ * x
          grads_coeff = grads_coeff/x
        END IF
        circ_eq_coeff = ModelDepth
        grads_coeff = grads_coeff/circ_eq_coeff
        C(1,1) = SUM( Tcoef(3,3,1:nn) * Basis(1:nn) )
        ! I * R, where 
        ! R = (1/sigma * js,js):
        ! ----------------------
        localR = REAL( Comp % N_j **2 * IP % s(t)*detJ/C(1,1)*circ_eq_coeff &
            / Comp % VoltageFactor, KIND=dp )
      CASE(3)
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), IP % W(t), &
            detJ, Basis, dBasisdx, EdgeBasis = Wbasis, RotBasis = RotWBasis, USolver = ASolver )

        IF (CoilUseWvec) THEN
          gradv = ListGetElementVectorSolution( Wvec_h, Basis, Element, dofs = dim )
        ELSE
          gradv = MATMUL( WBase(1:nn), dBasisdx(1:nn,:))
        END IF

        ! Compute the conductivity tensor
        ! -------------------------------
        DO i=1,3
          DO j=1,3
            C(i,j) = SUM( Tcoef(i,j,1:nn) * Basis(1:nn) )
            RotMLoc(i,j) = SUM( RotM(i,j,1:nn) * Basis(1:nn) )
          END DO
        END DO

        ! I * R, where 
        ! R = (1/sigma * js,js):
        ! ----------------------
        localR = REAL( Comp % N_j **2 * IP % s(t)*detJ/C(3,3) / Comp % VoltageFactor, KIND=dp )

        C = MATMUL(MATMUL(RotMLoc, C),TRANSPOSE(RotMLoc))

        IF (FoilUseJvec) THEN
          ! A given J vector is a real field, so it has no imaginary part.
          Jvec = CMPLX( ListGetElementVectorSolution( Jvec_h, Basis, Element, dofs = dim ), &
              0.0_dp, KIND=dp )
        ELSE
          Jvec = MATMUL(C,gradv)
        END IF

        ! Transform the conductivity tensor:
        ! ----------------------------------
      END SELECT
      
      localAlpha = SUM(alpha(1:nn) * Basis(1:nn))
      
      ! alpha is normalized to be in [0,1] thus, 
      ! it needs to be multiplied by the thickness of the coil 
      ! to get the real alpha:
      ! ------------------------------------------------------
      localAlpha = localAlpha * Comp % coilthickness

      Comp % Resistance = Comp % Resistance + localR

      DO vpolordtest=0,vpolord_tot ! V'(alpha)
        localVtest = localAlpha**vpolordtest
        dofIdtest = 2*(vpolordtest + 1) + vvarId
        DO vpolord = 0, vpolord_tot ! V(alpha)

          localV = localAlpha**vpolord
          dofId = 2*(vpolord + 1) + vvarId

          ! Computing the stiff term (sigma V(alpha) grad v0, V'(alpha) grad si):
          ! ---------------------------------------------------------------------
          IF (dim == 2) val = IP % s(t)*detJ*localV*localVtest*C(1,1)*grads_coeff**2*circ_eq_coeff
          IF (dim == 3) val = IP % s(t)*detJ*localV*localVtest*SUM(Jvec*gradv)
          val = val * Comp % VoltageFactor

          CALL AddToCmplxMatrixElement(CM, dofIdtest+nm, dofId+nm, REAL(val), AIMAG(val))
        END DO

        DO j=1,ncdofs
          q=j
          IF (dim == 3) q=q+nn
          ! computing the mass term (sigma * im * Omega * a, V'(alpha) grad si):
          ! ---------------------------------------------------------
          IF (dim == 2) val = im * Omega * IP % s(t)*detJ*localVtest*C(1,1)*basis(j)*grads_coeff*circ_eq_coeff
          IF (dim == 3) val = im * Omega * IP % s(t)*detJ*localVtest*SUM(MATMUL(C,Wbasis(j,:))*gradv)
          CALL AddToCmplxMatrixElement(CM, dofIdtest+nm, ReIndex(PS(Indexes(q))), REAL(val), AIMAG(val) )
        END DO

      END DO

      DO vpolord = 0, vpolord_tot ! V(alpha)
        localV = localAlpha**vpolord
        dofId = 2*(vpolord + 1) + vvarId

        DO j=1,ncdofs
            q=j
            IF (dim == 3) q=q+nn
            IF (dim == 2) val = IP % s(t)*detJ*localV*C(1,1)*basis(j)*grads_coeff
            IF (dim == 3) val = IP % s(t)*detJ*localV*SUM(Jvec*Wbasis(j,:))
            val = val * Comp % VoltageFactor
            CALL AddToCmplxMatrixElement(CM, ReIndex(PS(indexes(q))), dofId+nm, REAL(val), AIMAG(val))
        END DO
      END DO

    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE hm_Add_foil_winding
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AddMatrixEntry( A, row, col, val, tval )
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: A
    INTEGER :: row, col
    COMPLEX(KIND=dp) :: val, tval
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp) :: cval
    REAL(KIND=dp), POINTER :: svalues(:)
!------------------------------------------------------------------------------

    IF(EigenSystem) THEN
      CALL AddToCmplxMatrixElement(A, row, col, REAL(val), AIMAG(val))

      IF(.NOT.ASSOCIATED(A % MassValues)) THEN
        ALLOCATE(A % MassValues(SIZE(A % Values)))
        A % MassValues = 0._dp
      END IF

      svalues => A % Values
      A % Values => A % MassValues

      CALL AddToCmplxMatrixElement(A, row, col, REAL(tval), AIMAG(tval))

      A % Values => svalues
    ELSE
      cval = val+tval
      CALL AddToCmplxMatrixElement(A, row, col, REAL(cval), AIMAG(cval))
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE AddMatrixEntry
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
SUBROUTINE cda_out_eec() BIND(C, NAME="cda_out_eec")
  IMPLICIT NONE
  TYPE(Solver_t), POINTER :: ASolver
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Variable_t), POINTER :: AzVar
  REAL(KIND=dp), POINTER :: Az(:)
  REAL(KIND=dp), POINTER :: Acorr(:)
  REAL(KIND=dp) :: TTime
  INTEGER :: nm
  LOGICAL :: Found
  LOGICAL, SAVE :: EEC
  INTEGER, SAVE :: MyGen = -1
  LOGICAL :: EEC_lim
  REAL(KIND=dp), SAVE :: EEC_freq, EEC_time_0
  INTEGER, SAVE :: EEC_max, EEC_cnt = 0
  REAL(KIND=dp), ALLOCATABLE, SAVE :: Az0(:)
  ASolver => Ckt % ASolver
  nm = ASolver % Matrix % NumberOfRows
  ! Re-read whenever the circuit model is rebuilt, and also whenever a different
  ! one becomes active, since generation numbers are unique over all instances.
  ! EEC's Az0 is re-sampled on such a switch; with one output solver per circuit
  ! solver that never happens in practice, but these are the last saved values
  ! in the package that are not in the container.
  IF (MyGen /= CircuitsGeneration()) THEN
    SolverParams => GetSolverParams(Solver)
    ! Reading parameter for supply frequency
    EEC_freq = GetConstReal( SolverParams, 'EEC Frequency', EEC)
    IF (EEC) THEN
      CALL Info('CircuitsOutput', "Using EEC steady state forcing.", Level=4)
      WRITE( Message,'(A,4G11.4,A)') 'EEC signal frequency: ', EEC_freq, ' Hz'
      CALL Info('CircuitsOutput', Message, Level=4)
                
      EEC_max = GetInteger( SolverParams, 'EEC Steps', EEC_lim)
      IF (.NOT. EEC_lim) EEC_max = 5 !Typically 5 correections is enough
      WRITE( Message,'(A,I5,A)') 'Applying ', EEC_max, ' halfperiod corrections'
      CALL Info('CircuitsOutput', Message, Level=4)
      
      EEC_time_0 = 0.0
      
      ! Reserve memory for storing current MVP solution
      IF(ALLOCATED(Az0)) DEALLOCATE(Az0)
      ALLOCATE(Az0(nm))
      
      ! Store MVP solution at t=0
      AzVar => Asolver % Variable 
      IF(ASSOCIATED( AzVar)) THEN
        Az => AzVar % Values
        Az0 = Az
      END IF
      
    END IF
    
    out_CktPrefix = ListGetString(SolverParams,'Scalars Prefix',Found )
    IF(.NOT. Found) out_CktPrefix = 'res:'

    MyGen = CircuitsGeneration()
  END IF


  IF (EEC .AND. (EEC_cnt < EEC_max)) THEN
    
    TTime = GetTime()
    IF(TTime >= (EEC_time_0 + 0.5/EEC_freq)) THEN
      EEC_cnt = EEC_cnt + 1
      WRITE( Message,'(A,4G11.4)') 'Performing EEC #', EEC_cnt
      CALL Info('CircuitsOutput', Message, Level=4)
      
      EEC_time_0 = EEC_time_0 + 0.5/EEC_freq

      AzVar => Asolver % Variable 
      
      IF(ASSOCIATED( AzVar)) THEN
        Az => AzVar % Values
        
        !calculate correction
        ALLOCATE(Acorr(nm))
        Acorr = -0.5*(Az0+Az)
        Az = Az+Acorr
        DEALLOCATE(Acorr)
        
        !Store corrected half-period solution
        Az0 = Az
      END IF
    END IF
  END IF
  END SUBROUTINE cda_out_eec

SUBROUTINE cda_out_gather() BIND(C, NAME="cda_out_gather")
  IMPLICIT NONE
  TYPE(Solver_t), POINTER :: ASolver
  TYPE(Variable_t), POINTER :: LagrangeVar
  TYPE(Matrix_t), POINTER :: CM
  REAL(KIND=dp), ALLOCATABLE :: crt(:), crtt(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: sname
  LOGICAL :: Parallel
  INTEGER :: i, j, nm, circuit_tot_n
  ASolver => Ckt % ASolver
  nm = ASolver % Matrix % NumberOfRows
  CM => Ckt % CircuitMatrix
  circuit_tot_n = Ckt % Circuit_tot_n
  ! Circuit variable values from previous timestep:
  ! -----------------------------------------------
  ! This is the "is the linear system parallel" question, and the library has
  ! already answered it in Solver % Parallel. It used to be recomputed here from
  ! ParEnv % PEs, SingleMesh and "Enforce Parallel" - the same rule, except that
  ! it only looked for "Enforce Parallel" in the Simulation section and missed the
  ! solver section, so the two could disagree.
  Parallel = Solver % Parallel

  ALLOCATE(crt(circuit_tot_n), crtt(circuit_tot_n))
   crt = 0._dp
   crtt = 0._dp

   sname = LagrangeMultiplierName( ASolver )
   LagrangeVar => VariableGet( ASolver % Mesh % Variables,sname)
   IF(ASSOCIATED(LagrangeVar)) THEN
     CALL Info('CircuitsOutput','Initializing Lagrange multipliers of size: '&
         //I2S(SIZE(LagrangeVar % Values)),Level=8)
     IF( Parallel ) THEN
       DO i=1,circuit_tot_n 
         IF (ASSOCIATED(Model % CircuitModel % CircuitMatrix)) THEN  
           IF( CM % RowOwner(nm+i)==Parenv%myPE) crtt(i) = LagrangeVar%Values(i)
         END IF
       END DO
       CALL MPI_ALLREDUCE(crtt,crt,circuit_tot_n, MPI_DOUBLE_PRECISION, &
                  MPI_SUM, ASolver % Matrix % Comm, j)
     ELSE
       crt(1:circuit_tot_n) = LagrangeVar % Values(1:circuit_tot_n)
    END IF
  END IF
  IF (ALLOCATED(out_crt)) DEALLOCATE(out_crt)
  ALLOCATE(out_crt(circuit_tot_n))
  out_crt = crt
  END SUBROUTINE cda_out_gather
  SUBROUTINE cda_out_time() BIND(C, NAME="cda_out_time")
  IMPLICIT NONE
  CALL ListAddConstReal(GetSimulation(),TRIM(out_CktPrefix)//' time', GetTime())
  CALL Info('CircuitsOutput','Writing Circuit Results',Level=5)
  END SUBROUTINE cda_out_time
SUBROUTINE cda_out_circuit(p) BIND(C, NAME="cda_out_circuit")
  IMPLICIT NONE
  INTEGER(C_INT), VALUE :: p
  TYPE(Circuit_t), POINTER :: Circuits(:)
  TYPE(CircuitVariable_t), POINTER :: CVar
  TYPE(Component_t), POINTER :: Comp
  TYPE(ValueList_t), POINTER :: CompParams
  COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
  COMPLEX(KIND=dp) :: Current
  REAL(KIND=dp) :: CompRealPower, p_dc_component
  LOGICAL :: Found
  INTEGER :: i, jj, j
  Circuits => Ckt % Circuits
     CALL Info('CircuitsOutput', 'Writing Circuit Variables for Circuit '//i2s(p), Level=8)
     CALL Info('CircuitsOutput', 'There are '//i2s(Circuits(p)%n)//&
       ' Circuit Variables', Level=8)
     DO i=1,Circuits(p) % n
       Cvar => Circuits(p) % CircuitVariables(i)
       
       IF (Circuits(p) % Harmonic) THEN 
         CALL SimListAddAndOutputConstReal(&
           TRIM(Circuits(p) % names(i))//' re', out_crt(Cvar % ValueId), Level=10)
         CALL SimListAddAndOutputConstReal(&
           TRIM(Circuits(p) % names(i))//' im', out_crt(Cvar % ImValueId), Level=10)

         IF (Cvar % pdofs /= 0 ) THEN
           ! The harmonic assembly places polynomial term m (m=0,...,pdofs-1) of a
           ! foil winding voltage at Cvar % ValueId + AddIndex(m+1) and its
           ! imaginary part one row further, see Add_foil_winding() and
           ! CountAndCreateFoilWinding(). Hence AddIndex/AddImIndex here.
           ! ReIndex/ImIndex used to be used instead, which are 2*jj-1 and 2*jj
           ! and so reported every dof one slot too low: "re dof 1" was the
           ! imaginary part of the voltage itself, "im dof 1" the real part of
           ! dof 1, and so on.
           DO jj = 1, Cvar % pdofs
             CALL SimListAddAndOutputConstReal(&
               TRIM(Circuits(p) % names(i))&
               //' re dof '//I2S(jj), out_crt(Cvar % ValueId + AddIndex(jj)), Level=10)
             CALL SimListAddAndOutputConstReal(&
               TRIM(Circuits(p) % names(i))&
               //' im dof '//I2S(jj), out_crt(Cvar % ValueId + AddImIndex(jj)), Level=10)
           END DO
         END IF
       ELSE
         CALL SimListAddAndOutputConstReal(&
           TRIM(Circuits(p) % names(i)), out_crt(Cvar % ValueId), Level=10)
         
         IF (Cvar % pdofs /= 0 ) THEN
           DO jj = 1, Cvar % pdofs
             CALL SimListAddAndOutputConstReal(&
               TRIM(Circuits(p) % names(i))&
               //' dof '//I2S(jj), out_crt(Cvar % ValueId + jj), Level=10)
           END DO
         END IF
       END IF

     END DO

     CALL Info('CircuitsOutput', 'Writing Component Variables for Circuit '//i2s(p), Level=8)
     DO j = 1, SIZE(Circuits(p) % Components)
         Comp => Circuits(p) % Components(j)
         IF (ABS(Comp % Resistance) < TINY(0._dp) .AND. ABS(Comp % Conductance) > TINY(0._dp)) &
             Comp % Resistance = 1._dp / Comp % Conductance

         CALL SimListAddAndOutputConstReal('r_component('//&
           i2s(Comp % ComponentId)//')', Comp % Resistance, Level=8) 

         Current = 0._dp + im * 0._dp
         Current = out_crt(Comp % ivar % ValueId) 
         IF ( Circuits(p) % Harmonic ) Current = Current + im * out_crt(Comp % ivar % ImValueId) 
              
         CompParams => CurrentModel % Components (Comp % ComponentId) % Values
         IF (.NOT. ASSOCIATED(CompParams)) CALL Fatal ('CircuitsOutput', &
           'Component parameters not found!')

         p_dc_component = ABS(Current)**2._dp * Comp % Resistance
         CALL SimListAddAndOutputConstReal('p_dc_component('//i2s(Comp % ComponentId)//')',&
           p_dc_component, Level=8) 

         CompRealPower = GetConstReal( Model % Simulation, &
             TRIM(out_CktPrefix)//' Power re in Component '//i2s(Comp % ComponentId), Found)
         IF (Found .AND. ABS(Current) > TINY(CompRealPower)) THEN
           CALL SimListAddAndOutputConstReal('p_ac_component('//&
             i2s(Comp % ComponentId)//')', CompRealPower, Level=8)
           CALL SimListAddAndOutputConstReal('r_ac_component('//&
             i2s(Comp % ComponentId)//')', CompRealPower/ABS(Current)**2._dp, Level=8)
           CALL SimListAddAndOutputConstReal('AC to DC of component '&
             //i2s(Comp % ComponentId), CompRealPower/p_dc_component, Level=8)
         END IF
          
       END DO  
  END SUBROUTINE cda_out_circuit




  
!-------------------------------------------------------------------
  SUBROUTINE SimListAddAndOutputConstReal(VariableName, VariableValue, Level)
!-------------------------------------------------------------------
  IMPLICIT NONE
  CHARACTER(LEN=MAX_NAME_LEN) :: VarVal
  CHARACTER(LEN=*) :: VariableName
  REAL(KIND=dp) :: VariableValue
  INTEGER, OPTIONAL :: Level 
  ! Not initialized in the declaration: that would give it an implicit SAVE and
  ! a call without Level would reuse whatever the previous call passed.
  INTEGER :: LevelVal

  LevelVal = 3
  IF (PRESENT(Level)) LevelVal = Level
  ! Wide enough for the longest name written here, "<name> im dof <n>".
  WRITE(Message,'(A,T34,ES15.4)') TRIM(VariableName),VariableValue
  CALL Info('CircuitsOutput',Message,Level=LevelVal)
  
  CALL ListAddConstReal(GetSimulation(),TRIM(out_CktPrefix)//' '//TRIM(VariableName), VariableValue)
!-------------------------------------------------------------------
  END SUBROUTINE SimListAddAndOutputConstReal
!-------------------------------------------------------------------

END MODULE cda_acc
