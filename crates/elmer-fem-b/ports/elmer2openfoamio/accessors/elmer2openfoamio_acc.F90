! T022-A Elmer2OpenFoamIO port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/Elmer2OpenFoamIO.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, out INTEGER(C_INT) flags). All
! Fortran-object traffic — TYPE(...) handles, OPTIONAL/POINTER dummies,
! assumed-shape descriptors, hidden lengths, fixed-length string buffers
! — stays inside Fortran, so the Rust side is pure C calls with no layout
! assumption and no pinned offset.
!
! The four CONTAINS internals of the original (OpenFOAMBlocks,
! WriteMeshCenters, CreateFOAMMesh, WriteOFField) are kept whole as
! module procedures with identical statements; only host association is
! swapped for module SAVE state (acc_params, acc_basedir, acc_nodir,
! acc_filename, acc_dirname, acc_offile, acc_mesh, acc_ofmesh,
! acc_solver, acc_varname), which has the same lifetimes as the
! original's SAVE'd host variables. Order-sensitive I/O (formatted
! reads/writes, MINVAL/MAXVAL/SUM reductions, whole-array zeroing) stays
! in Fortran, so bit-identity holds by construction wherever the call
! sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module elmer2openfoamio_acc
  use iso_c_binding
  use defutils
  use interpolation
  use meshbasics
  use elementutils
  use particleutils
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(valuelist_t), pointer, save :: acc_params => null()
  type(variable_t), pointer, save :: acc_var => null()
  type(variable_t), pointer, save :: acc_ofvar => null()
  type(mesh_t), pointer, save :: acc_mesh => null()
  type(mesh_t), pointer, save :: acc_ofmesh => null()
  character(len=MAX_NAME_LEN), save :: acc_varname = ''
  character(len=MAX_NAME_LEN), save :: acc_filename = ''
  character(len=MAX_NAME_LEN), save :: acc_dirname = ''
  character(len=MAX_NAME_LEN), save :: acc_basedir = ''
  character(len=MAX_NAME_LEN), save :: acc_offile = ''
  integer, save :: acc_nodir = 0
  integer, save :: acc_passive = 0

contains

  ! -- entry state ----------------------------------------------------
  subroutine e2o_begin(h_model, h_solver) bind(c, name="e2o_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine e2o_begin

  subroutine e2o_info_head() bind(c, name="e2o_info_head")
    CALL Info('Elmer2OpenFoamWrite','-----------------------------------------', Level=4 )
    CALL Info('Elmer2OpenFoamWrite','Projecting field to OpenFOAM cell centers',Level=4)
  end subroutine e2o_info_head

  subroutine e2o_params() bind(c, name="e2o_params")
    acc_params => GetSolverParams()
  end subroutine e2o_params

  subroutine e2o_target(found) bind(c, name="e2o_target")
    integer(c_int) :: found
    logical :: f
    acc_varname = GetString(acc_params,'Target Variable',f)
    found = merge(1, 0, f)
  end subroutine e2o_target

  subroutine e2o_fatal_no_target() bind(c, name="e2o_fatal_no_target")
    CALL Fatal('Elmer2OpenFoamWrite','> Target Variable < must exist for the solver!')
  end subroutine e2o_fatal_no_target

  subroutine e2o_proj_save(psave, pfound) bind(c, name="e2o_proj_save")
    integer(c_int) :: psave, pfound
    logical :: s, f
    s = ListGetLogical(acc_model % Simulation,'Use Mesh Projector',f)
    psave = merge(1, 0, s)
    pfound = merge(1, 0, f)
    CALL Info('Elmer2OpenFOAMWrite','Enforcing mapping without projector matrix!',Level=6)
    CALL ListAddLogical(acc_model % Simulation,'Use Mesh Projector',.FALSE.)
  end subroutine e2o_proj_save

  subroutine e2o_getmesh() bind(c, name="e2o_getmesh")
    acc_mesh => GetMesh()
  end subroutine e2o_getmesh

  function e2o_var_get() result(r) bind(c, name="e2o_var_get")
    integer(c_int) :: r
    acc_var => VariableGet(acc_mesh % Variables, acc_varname)
    r = merge(1, 0, associated(acc_var))
  end function e2o_var_get

  subroutine e2o_fatal_no_var() bind(c, name="e2o_fatal_no_var")
    CALL Fatal('Elmer2OpenFoamWrite','Variable does not exist in Elmer mesh!')
  end subroutine e2o_fatal_no_var

  subroutine e2o_release_saved() bind(c, name="e2o_release_saved")
    CALL ReleaseMesh(acc_ofmesh)
    DEALLOCATE(acc_ofmesh)
  end subroutine e2o_release_saved

  subroutine e2o_allocate() bind(c, name="e2o_allocate")
    acc_ofmesh => AllocateMesh()
  end subroutine e2o_allocate

  subroutine e2o_passive(val, found) bind(c, name="e2o_passive")
    integer(c_int) :: val, found
    logical :: f
    acc_passive = ListGetInteger(acc_solver % Values,'Passive OpenFOAM Coordinate',f)
    val = acc_passive
    found = merge(1, 0, f)
  end subroutine e2o_passive

  function e2o_meshdim() result(n) bind(c, name="e2o_meshdim")
    integer(c_int) :: n
    n = acc_mesh % MeshDim
  end function e2o_meshdim

  subroutine e2o_warn_dim() bind(c, name="e2o_warn_dim")
    CALL Warn('Elmer2OpenFOAM','Dimension of Elmer mesh is reduced, and OpenFOAM not?!')
  end subroutine e2o_warn_dim

  subroutine e2o_center_name(found) bind(c, name="e2o_center_name")
    integer(c_int) :: found
    logical :: f
    acc_filename = GetString(acc_params,'Elmer Center Filename',f)
    found = merge(1, 0, f)
  end subroutine e2o_center_name

  subroutine e2o_write_centers() bind(c, name="e2o_write_centers")
    CALL e2o_centers_sub(acc_mesh, acc_filename)
  end subroutine e2o_write_centers

  subroutine e2o_basedir(found) bind(c, name="e2o_basedir")
    integer(c_int) :: found
    logical :: f
    acc_basedir = GetString(acc_params,'OpenFOAM Directory',f)
    found = merge(1, 0, f)
  end subroutine e2o_basedir

  subroutine e2o_info_basedir() bind(c, name="e2o_info_basedir")
    CALL Info('Elmer2OpenFoamWrite','Using given > OpenFOAM Directory < : '//TRIM(acc_basedir),Level=6)
  end subroutine e2o_info_basedir

  subroutine e2o_fatal_no_basedir() bind(c, name="e2o_fatal_no_basedir")
    CALL Fatal('Elmer2OpenFoamWrite','> OpenFOAM Directory < must exist for the solver!')
  end subroutine e2o_fatal_no_basedir

  function e2o_has_mesh1() result(r) bind(c, name="e2o_has_mesh1")
    integer(c_int) :: r
    r = merge(1, 0, ListCheckPresent(acc_params,'OpenFOAM Mesh 1'))
  end function e2o_has_mesh1

  subroutine e2o_blocks() bind(c, name="e2o_blocks")
    CALL e2o_blocks_sub()
  end subroutine e2o_blocks

  subroutine e2o_reduce() bind(c, name="e2o_reduce")
    acc_nodir = ParallelReduction(acc_nodir)
  end subroutine e2o_reduce

  function e2o_nodir() result(n) bind(c, name="e2o_nodir")
    integer(c_int) :: n
    n = acc_nodir
  end function e2o_nodir

  subroutine e2o_info_nodir() bind(c, name="e2o_info_nodir")
    CALL Info('Elmer2OpenFOAMWrite','Number of active OpenFOAM blocks: '//I2S(acc_nodir),Level=5)
  end subroutine e2o_info_nodir

  function e2o_mype() result(n) bind(c, name="e2o_mype")
    integer(c_int) :: n
    n = ParEnv % MyPe
  end function e2o_mype

  subroutine e2o_info_treating(i) bind(c, name="e2o_info_treating")
    integer(c_int), value :: i
    CALL Info('Elmer2OpenFOAMWrite','Treating OpenFOAM block: '//I2S(i),Level=5)
  end subroutine e2o_info_treating

  subroutine e2o_dirname(i, found) bind(c, name="e2o_dirname")
    integer(c_int), value :: i
    integer(c_int) :: found
    logical :: f
    acc_dirname = ListGetString(acc_params,'OpenFOAM Mesh '//I2S(i),f)
    found = merge(1, 0, f)
  end subroutine e2o_dirname

  subroutine e2o_fatal_no_key() bind(c, name="e2o_fatal_no_key")
    CALL Fatal('Elmer2OpenFoamWrite','Could not find keyword: '//TRIM(acc_dirname))
  end subroutine e2o_fatal_no_key

  subroutine e2o_filename_c() bind(c, name="e2o_filename_c")
    acc_filename = TRIM(acc_dirname)//'C'
  end subroutine e2o_filename_c

  subroutine e2o_info_project() bind(c, name="e2o_info_project")
    CALL Info('Elmer2OpenFoamWrite','Projecting to OpenFOAM nodes in file: '//TRIM(acc_filename),Level=5)
  end subroutine e2o_info_project

  subroutine e2o_create() bind(c, name="e2o_create")
    CALL e2o_create_sub(acc_filename, acc_ofmesh)
  end subroutine e2o_create

  subroutine e2o_current_to_of() bind(c, name="e2o_current_to_of")
    CurrentModel % Mesh => acc_ofmesh
  end subroutine e2o_current_to_of

  subroutine e2o_zero_passive() bind(c, name="e2o_zero_passive")
    IF( acc_passive == 1 ) THEN
      acc_ofmesh % Nodes % x = 0.0_dp
    ELSE IF( acc_passive == 2 ) THEN
      acc_ofmesh % Nodes % y = 0.0_dp
    ELSE IF( acc_passive == 3 ) THEN
      acc_ofmesh % Nodes % z = 0.0_dp
    END IF
  end subroutine e2o_zero_passive

  subroutine e2o_info_mapping() bind(c, name="e2o_info_mapping")
    CALL Info('Elmer2OpenFoamWrite','Mapping data to the temporal mesh using library routines',Level=6)
  end subroutine e2o_info_mapping

  subroutine e2o_ofvar_get() bind(c, name="e2o_ofvar_get")
    acc_ofvar => VariableGet(acc_ofmesh % Variables, acc_varname)
  end subroutine e2o_ofvar_get

  subroutine e2o_solver_var() bind(c, name="e2o_solver_var")
    acc_solver % Variable => acc_ofvar
  end subroutine e2o_solver_var


  subroutine e2o_info_stats() bind(c, name="e2o_info_stats")
    REAL(KIND=dp) :: MinF, MaxF, MeanF
    MinF = MINVAL( acc_ofvar % Values )
    MaxF = MAXVAL( acc_ofvar % Values )
    MeanF = SUM( acc_ofvar % Values ) / SIZE( acc_ofvar % Values )
    WRITE( Message,'(A,ES12.5)') 'Minimum field value: ',MinF
    CALL Info('Elmer2OpenFoamWrite',Message,Level=6)
    WRITE( Message,'(A,ES12.5)') 'Maximum field value: ',MaxF
    CALL Info('Elmer2OpenFoamWrite',Message,Level=6)
    WRITE( Message,'(A,ES12.5)') 'Average field value: ',MeanF
    CALL Info('Elmer2OpenFoamWrite',Message,Level=6)
  end subroutine e2o_info_stats

  subroutine e2o_write_field() bind(c, name="e2o_write_field")
    logical :: f
    acc_offile = GetString(acc_params,'OpenFOAM file',f)
    IF(.NOT. f ) acc_offile = 'fieldSolidHS.dat'
    acc_filename = TRIM(acc_dirname)//TRIM(acc_offile)
    CALL Info('Elmer2OpenFoamWrite','Writing projected field to file: '//TRIM(acc_filename),Level=5)
    CALL e2o_write_sub(acc_filename, acc_ofmesh, acc_ofvar)
  end subroutine e2o_write_field

  subroutine e2o_release_of() bind(c, name="e2o_release_of")
    CALL ReleaseMesh(acc_ofmesh)
  end subroutine e2o_release_of

  subroutine e2o_current_to_mesh() bind(c, name="e2o_current_to_mesh")
    CurrentModel % Mesh => acc_mesh
  end subroutine e2o_current_to_mesh

  subroutine e2o_proj_restore(psave, pfound) bind(c, name="e2o_proj_restore")
    integer(c_int), value :: psave, pfound
    IF( pfound /= 0 ) THEN
      CALL ListAddLogical(acc_model % Simulation,'Use Mesh Projector',psave /= 0)
    ELSE
      CALL ListRemove(acc_model % Simulation,'Use Mesh Projector')
    END IF
  end subroutine e2o_proj_restore

  subroutine e2o_info_done() bind(c, name="e2o_info_done")
    CALL Info('Elmer2OpenFoamWrite','All done', Level=4 )
    CALL Info('Elmer2OpenFoamWrite','-----------------------------------------', Level=4 )
  end subroutine e2o_info_done

  ! -- verbatim internal procedures (host refs swapped for SAVE state)  --

  subroutine e2o_blocks_sub( )

    CHARACTER(LEN=MAX_NAME_LEN) :: DirCommand
    LOGICAL :: FileExists
    INTEGER, PARAMETER :: InFileUnit = 28
    INTEGER :: IOStatus


    acc_nodir = 0
    IF( ParEnv % MyPe /= 0 ) RETURN

#ifdef __INTEL_COMPILER
    ! Fortran standard states that inquiry for a file returns true if the queried entity is a file
    INQUIRE( Directory = TRIM(acc_basedir), Exist = FileExists )
#else
    INQUIRE( File = TRIM(acc_basedir), Exist = FileExists )
#endif
    IF(.NOT. FileExists ) THEN
      CALL Fatal('Elmer2OpenFoamWrite','OpenFOAM directory does not exist: '//TRIM(acc_basedir))
    END IF

    acc_dirname = TRIM(acc_basedir)//'/0/'
#ifdef __INTEL_COMPILER
    INQUIRE( Directory = acc_dirname, Exist = FileExists )
#else
    INQUIRE( File = acc_dirname, Exist = FileExists )
#endif
    IF(.NOT. FileExists ) THEN
      CALL Fatal('Elmer2OpenFoamWrite','OpenFOAM mesh does not exist: '//TRIM(acc_dirname))
    END IF

    acc_filename = TRIM(acc_dirname)//'C'
    CALL Info('Elmer2OpenFoamWrite','Inquire file: '//TRIM(acc_filename),Level=12)
    INQUIRE( File = acc_filename, Exist = FileExists )

    IF( FileExists ) THEN
      CALL Info('Elmer2OpenFoamWrite','Using OpenFOAM centers in: '//TRIM(acc_filename),Level=10)
      CALL ListAddString( acc_params, 'OpenFOAM Mesh 1', acc_dirname, .FALSE.)
      acc_nodir = 1
      RETURN
    END IF

    DirCommand = 'ls -d '//TRIM(acc_dirname)//'*/ > OpenFOAMBlocks.txt'
    CALL Info('Elmer2OpenFoamWrite','Performing command: '//TRIM(DirCommand),Level=12)
    CALL SystemCommand( DirCommand )

    OPEN(InFileUnit,File='OpenFOAMBlocks.txt',IOStat=IOstatus)
    IF(IOStatus /= 0 ) THEN
      CALL Fatal('Elmer2OpenFoamWrite','Could not open file: OpenFOAMBlocks.txt')
    END IF

    DO
      READ(InFileUnit,'(A)',IOStat = IOStatus) acc_dirname
      IF( IOStatus /= 0 ) EXIT
      acc_filename = TRIM(acc_dirname)//'C'
      CALL Info('Elmer2OpenFoamWrite','Inquire file: '//TRIM(acc_filename),Level=12)
      INQUIRE( File = acc_filename, Exist = FileExists )
      IF( FileExists ) THEN
        acc_nodir = acc_nodir + 1
        CALL Info('Elmer2OpenFoamWrite','Using OpenFOAM centers in: '//TRIM(acc_filename),Level=10)
        CALL ListAddString( acc_params, 'OpenFOAM Mesh '//I2S(acc_nodir), acc_dirname, .FALSE.)
      ELSE
        CALL Info('Elmer2OpenFoamWrite','No OpenFOAM center in: '//TRIM(acc_dirname),Level=12)
      END IF
    END DO
    CLOSE(InFileUnit)

    IF( acc_nodir == 0 ) THEN
      CALL Fatal('Elmer2OpenFoamWrite','No OpenFOAM mesh blocks found!')
    ELSE
      CALL Info('Elmer2OpenFoamWrite','Number of OpenFOAM blocks: '//I2S(acc_nodir),Level=10)
    END IF

  end subroutine e2o_blocks_sub


  subroutine e2o_centers_sub( Mesh, Filename )

    CHARACTER(LEN=MAX_NAME_LEN) :: FileName
    TYPE(Mesh_t), TARGET :: Mesh

    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: x,y,z
    INTEGER, PARAMETER :: OutFileUnit = 29
    INTEGER :: i,n,IOstatus


    OPEN(OutFileUnit,FILE = Filename, IOSTAT=IOStatus)
    IF( IOStatus /= 0 ) THEN
      CALL Fatal('Elmer2OpenFoamWrite','Could not open file for writing: '//TRIM(FileName))
    END IF

    CALL Info('Elmer2OpenFoamWrite','Writing Elmer element centers to file: '//TRIM(FileName),Level=6)

    DO i=1, Mesh % NumberOfBulkElements
      Element => Mesh % Elements(i)
      n = Element % TYPE % NumberOfNodes
      x = SUM( Mesh % Nodes % x( Element % NodeIndexes ) ) / n
      y = SUM( Mesh % Nodes % y( Element % NodeIndexes ) ) / n
      z = SUM( Mesh % Nodes % z( Element % NodeIndexes ) ) / n
      WRITE( OutFileUnit, * ) '(',x,y,z,')'
    END DO

    CLOSE( OutFileUnit )

  end subroutine e2o_centers_sub



  !------------------------------------------------------------------------
  !> Open file in OpenFOAM format and read the cell centers from there.
  !-------------------------------------------------------------------------
  subroutine e2o_create_sub( Filename, Mesh )

    CHARACTER(LEN=MAX_NAME_LEN) :: FileName
    TYPE(Mesh_t), TARGET :: Mesh

    INTEGER :: line,i,j,k,n
    REAL(KIND=dp) :: x,y,z
    INTEGER :: NumberOfNodes, IOStatus
    INTEGER, PARAMETER :: InFileUnit = 28
    CHARACTER(LEN=:), ALLOCATABLE :: ReadStr
    LOGICAL :: InlineCoords
    LOGICAL :: Found

    ALLOCATE( Mesh % Nodes )
    Mesh % NumberOfBulkElements = 0
    Mesh % NumberOfBoundaryElements = 0

    ! Partition zero does all the work!
    IF( ParEnv % MyPe /= 0 ) THEN
      Mesh % NumberOfNodes = 0
      GOTO 100
    END IF


    ALLOCATE(CHARACTER(MAX_STRING_LEN)::ReadStr)

    OPEN(InFileUnit,FILE = Filename, STATUS='old', IOSTAT=IOstatus)
    IF( IOStatus /= 0 ) THEN
      CALL Fatal('Elmer2OpenFoamWrite','Could not open file for reading: '//TRIM(Filename))
    END IF

    CALL Info('Elmer2OpenFoamWrite','Reading data points from file: '//TRIM(Filename),Level=7)

    j = 0
    k = 0
    DO Line = 1, 100
      READ( InFileUnit,'(A)',IOSTAT=IOStatus ) ReadStr
      IF( IOStatus /= 0 ) THEN
        CALL Warn('Elmer2OpenFoamWrite','End of file after '//I2S(Line)//' lines')
        EXIT
      END IF

      j =  INDEX( ReadStr,'internalField',.TRUE.)
      IF( j > 0 ) THEN
        ! If we have parenthesis in the same line as "internalField" then the coordinate
        ! values are in-lined.
        k = INDEX( ReadStr,'(')
        EXIT
      END IF
    END DO

    IF( j == 0 ) THEN
      CALL Fatal('Elmer2OpenFoamWrite','Could not find > internalField < in header!')
    ELSE
      CALL Info('Elmer2OpenFoamWrite','internalField found at line: '//I2S(Line),Level=10)
    END IF

    InlineCoords = ( k > 0 )
    IF( InlineCoords ) THEN
      j = INDEX( ReadStr,'<vector>')
      READ( ReadStr(j+8:k-1),*,IOSTAT=IOStatus ) NumberOfNodes
      CALL Info('Elmer2OpenFoamWrite','Reading inline coordinates',Level=10)
    ELSE
      READ(InFileUnit,*,IOSTAT=IOStatus) NumberOfNodes
    END IF

    IF( IOStatus /= 0 ) THEN
      CALL Fatal('Elmer2OpenFoamWrite','Could not read number of nodes!')
    END IF
    CALL Info('Elmer2OpenFoamWrite','Number of OpenFOAM nodes: '&
        //I2S(NumberOfNodes),Level=10)

    i = ListGetInteger(acc_params,'Number of cells',Found)
    IF( i > 0 .AND. i < NumberOfNodes ) THEN
      NumberOfNodes = i
      CALL Info('Elmer2OpenFoamWrite','Limiting number of OpenFOAM nodes: '&
          //I2S(NumberOfNodes),Level=7)
    END IF


    n = NumberOfNodes
    ALLOCATE( Mesh % Nodes % x(n), &
        Mesh % Nodes % y(n), &
        Mesh % Nodes % z(n) )

    Mesh % Nodes % x(1:n) = 0.0_dp
    Mesh % Nodes % y(1:n) = 0.0_dp
    Mesh % Nodes % z(1:n) = 0.0_dp

    Mesh % NumberOfNodes = n


    ! This is just empty left parenthesis
    IF(.NOT. InlineCoords ) THEN
      READ( InFileUnit,'(A)',IOSTAT=IOStatus ) ReadStr
    END IF

    DO i=1,n
      IF( InlineCoords ) THEN
        ReadStr = TRIM( ReadStr(k+1:) )
      ELSE
        READ( InFileUnit,'(A)',IOSTAT=IOStatus ) ReadStr
        IF( IOStatus /= 0 ) THEN
          CALL Fatal('Elmer2OpenFoamWrite','Could not read coordinate line: '//I2S(i))
        END IF
      END IF

      IF( InlineCoords ) THEN
        j =  INDEX( ReadStr,'(',.FALSE.)
        k =  INDEX( ReadStr,')',.FALSE.)
      ELSE
        j =  INDEX( ReadStr,'(',.TRUE.)
        k =  INDEX( ReadStr,')',.TRUE.)
      END IF

      IF( j == 0 ) THEN
        CALL Fatal('Elmer2OpenFoamWrite',&
            'Expecting a parenthesis at the start of OpenFOAM line: '//I2S(i))
      END IF
      IF( k == 0 ) THEN
        CALL Fatal('Elmer2OpenFoamWrite',&
            'Expecting a parenthesis at the end of OpenFOAM line: '//I2S(i))
      END IF

      READ( ReadStr(j+1:k-1),*,IOSTAT=IOStatus ) x,y,z
      IF( IOStatus /= 0 ) THEN
        CALL Fatal('Elmer2OpenFoamWrite','Could not read coordinate values: '//I2S(i))
      END IF
      Mesh % Nodes % x(i) = x
      Mesh % Nodes % y(i) = y
      Mesh % Nodes % z(i) = z
    END DO
    CLOSE( InFileUnit )


    CALL Info('Elmer2OpenFoamWrite','Creating coordinates for temporal mesh',Level=7)

100 CALL VariableAdd( Mesh % Variables, Mesh, acc_solver, &
        'Coordinate 1',1, Mesh % Nodes % x )

    CALL VariableAdd( Mesh % Variables, Mesh, acc_solver, &
        'Coordinate 2',1, Mesh % Nodes % y )

    CALL VariableAdd( Mesh % Variables, Mesh, acc_solver, &
        'Coordinate 3',1,Mesh % Nodes % z )

    CALL Info('Elmer2OpenFoamWrite','Created temporal OpenFOAM mesh just for nodes',Level=8)

  end subroutine e2o_create_sub


  !------------------------------------------------------------------------
  !> Write the data to the OpenFOAM cell centers. This assumes same order
  !> for the cells as was used in reading in the temporal mesh.
  !-------------------------------------------------------------------------
  subroutine e2o_write_sub( Filename, Mesh, Var )
    CHARACTER(LEN=MAX_NAME_LEN) :: FileName
    TYPE(Mesh_t) :: Mesh
    TYPE(Variable_t) :: Var
    INTEGER :: i,j,IOStatus

    INTEGER, PARAMETER :: OutFileUnit = 29
    CHARACTER :: NL


    NL = NEW_LINE('A')

    CALL Info('Elmer2OpenFoamWrite','Writing a file with field value for OpenFOAM: '//TRIM(FileName),Level=6)

    OPEN(OutFileUnit,FILE = Filename, IOSTAT=IOStatus)

    WRITE(OutFileUnit,'(A)') &
        '/*--------------------------------*- C++ -*----------------------------------*\ '//NL//&
        '| =========                 |                                                 | '//NL//&
        '| \\      /  F ield         | OpenFOAM: The OPEN Source CFD Toolbox           | '//NL//&
        '|  \\    /   O peration     | Version:  dev                                   | '//NL//&
        '|   \\  /    A nd           | Web:      www.OpenFOAM.org                      | '//NL//&
        '|    \\/     M anipulation  |                                                 | '//NL//&
        '|                                                                             | '//NL//&
        '|  This file written by Elmer2OpenFOAM coupling module in Elmer               | '//NL//&
        '\*---------------------------------------------------------------------------*/ '
    WRITE(OutFileUnit,'(A)') &
        'FoamFile                                '//NL//&
        '{                                       '//NL//&
        '    version     2.0;                    '//NL//&
        '    format      ascii;                  '//NL//&
        '    class       volScalarField;         '//NL//&
        '    object      fieldSolidHS;           '//NL//&
        '}'
    WRITE(OutFileUnit,'(A)') &
        'dimensions      [1 -1 -3 0 0 0 0];      '

      WRITE(OutFileUnit,'(A)') &
          'boundaryField                           '//NL//&
          '{                                       '//NL//&
          '    ".*"                                '//NL//&
          '    {                                   '//NL//&
          '        TYPE zeroGradient;              '//NL//&
          '    }                                   '//NL//&
          '}'

      WRITE(OutFileUnit,'(A)') &
          'internalField   nonuniform List<scalar> '//NL//&
          I2S(Mesh % NumberOfNodes)           //NL//&
          '('

      DO i=1,Mesh % NumberOfNodes
        IF( ASSOCIATED( Var % Perm ) ) THEN
          j = Var % Perm(i)
        ELSE
          j = i
        END IF
        IF( j < 0 .OR. j > Mesh % NumberOfNodes ) THEN
          CALL Fatal('Elmer2OpenFoamWrite','We have a troubling ENTRY: '//I2S(i)//','//I2S(j))
        END IF

        WRITE( OutFileUnit,*) Var % Values(j)
      END DO

      WRITE(OutFileUnit,'(A)') ');'
      CLOSE( OutFileUnit )
      CALL Info('Elmer2OpenFoamWrite','Created a file with OpenFOAM field values',Level=8)

  end subroutine e2o_write_sub

end module elmer2openfoamio_acc
