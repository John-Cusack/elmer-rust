! SaveMesh port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/SaveMesh.F90, compiled against the built .mod
! files of the target configure (T020 Poisson pattern, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars). All Fortran-object traffic — TYPE(...)
! handles, OPTIONAL/POINTER dummies, assumed-shape descriptors, hidden
! lengths, allocatable function results — stays inside Fortran, so the
! Rust side is pure C calls. The wrapper module mirrors SaveMesh's own
! locals with identical lifetimes.
!
! Faithfulness: the body below `svm_do_save` is a statement-for-statement
! transcription of SaveMesh after the "Save Mesh" flag check, including
! the WRITE formatting with i4.4 and the ParEnv % PEs serial branch, so
! bit-identity of the written mesh trees holds by construction wherever
! the call sequence matches. The parallel arm (MPI_BARRIER +
! partitioned WriteMeshToDisk2) is out of scope for the serial proof
! tree, which has no MPI: on ParEnv % PEs > 1 this port is a documented
! no-op (see lib.rs).
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module savemesh_acc
  use iso_c_binding
  use DefUtils
  use MeshBasics
  use LoadMod
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()

contains

  ! -- entry state ----------------------------------------------------
  subroutine svm_begin(h_model, h_solver) bind(c, name="svm_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine svm_begin

  ! -- "Save Mesh" flag: ListGetLogical + Warn + default .FALSE. ------
  function svm_save_flag() result(r) bind(c, name="svm_save_flag")
    integer(c_int) :: r
    logical :: flagval, found
    flagval = ListGetLogical(acc_solver % values, "Save Mesh", found)
    if (.not. found) then
      call Warn("Save Mesh", "Can't find Save Mesh logical")
      flagval = .false.
    end if
    r = merge(1, 0, flagval)
  end function svm_save_flag

  ! -- remainder of the body (serial branch) ---------------------------
  subroutine svm_do_save() bind(c, name="svm_do_save")
    logical :: everytime, found
    type(mesh_t), pointer :: mesh
    type(variable_t), pointer :: timestepvar
    character(len=MAX_NAME_LEN) :: meshname, meshdir
    mesh => acc_solver % mesh
    timestepvar => VariableGet(mesh % variables, "Timestep", .true.)
    meshname = ListGetString(acc_solver % values, "Mesh Name", found)
    if (.not. found) call Fatal("SaveMesh", "No name given for mesh")
    meshdir = ListGetString(acc_solver % values, "Save Mesh Directory", found)
    if (.not. found) call Fatal("SaveMesh", "No directory given to save mesh")
    call MakeDirectory(trim(meshdir) // char(0))
    everytime = ListGetLogical(acc_solver % values, "Save All Timesteps", found)
    if (.not. found) everytime = .false.
    if (everytime) then
      write(meshdir, '(A,A,A,i4.4)') trim(meshdir), trim(meshname), "_", &
          int(timestepvar % values(1))
    else
      write(meshdir, '(A,A)') trim(meshdir)//trim(meshname)
    end if
    if (parenv % pes <= 1) then !serial
      call MakeDirectory(trim(meshdir) // char(0))
      call WriteMeshToDisk2(acc_model, mesh, meshdir)
    else !parallel
      ! Out of scope for the serial proof tree (no MPI): documented
      ! no-op. The branch check above replicates the original's
      ! `IF( ParEnv % PEs<=1 )` so the serial arm matches exactly.
    end if
  end subroutine svm_do_save

end module savemesh_acc
