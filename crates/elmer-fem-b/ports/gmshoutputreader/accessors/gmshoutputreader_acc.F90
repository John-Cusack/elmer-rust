! T022-A GmshOutputReader port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/GmshOutputReader.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars, INTEGER(C_INT) codes). All Fortran-object
! traffic — TYPE(...) handles, OPTIONAL/POINTER dummies, assumed-shape
! descriptors, hidden lengths, list-directed file I/O, whole-array
! updates, MINVAL/MAXVAL reductions — stays inside Fortran, so the Rust
! side is pure C calls. The wrapper module mirrors the original's locals
! (entry-level state, the ReadSingleGmshFile SAVE GmshVer, the VarArray
! table) with identical lifetimes.
!
! The two CONTAINS internal procedures of the original (ReadSingleGmshFile,
! InterpolateFromGmshFile) become module procedures over the same
! host-associated state, which now lives as module SAVE variables; the
! Rust entry drives the outer sweep loop (the label-10 GOTO) through the
! gmr_pass / gmr_need_alloc / gmr_allocate queries, matching the
! original's call sequence exactly. All floating-point work, file reads,
! and formatted Info output stay in Fortran (same compiler), so
! bit-identity holds by construction wherever the call sequence matches.
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module gmshoutputreader_acc
  use iso_c_binding
  use defutils
  use saveutils
  implicit none

  character(*), parameter :: caller = 'GmshOutputReader'
  integer, parameter :: length = 1024

  type vararray_t
    type(variable_t), pointer :: var => null()
  end type vararray_t

  type(vararray_t) :: acc_vararray(20)

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()
  type(valuelist_t), pointer, save :: acc_params => null()
  type(mesh_t), pointer, save :: acc_tomesh => null()
  type(mesh_t), pointer, save :: acc_frommesh => null()
  type(variable_t), pointer, save :: acc_var => null()
  type(element_t), pointer, save :: acc_element => null()
  integer, pointer, save :: acc_perm(:) => null()
  integer, pointer, save :: acc_maskperm(:) => null()
  integer, pointer, save :: acc_nodeindexes(:) => null()
  integer, pointer, save :: acc_elmerindexes(:) => null()

  logical, save :: acc_found = .false.
  logical, save :: acc_usebbox = .false.
  logical, save :: acc_usequadtree = .false.
  logical, save :: acc_allocdone = .false.
  integer(c_int), save :: acc_dim = 0
  integer(c_int), save :: acc_aligncoord = 0
  integer(c_int), save :: acc_passivecoord = 0
  integer(c_int), save :: acc_inputpartitions = 1
  integer(c_int), save :: acc_dofs = 0
  integer(c_int), save :: acc_n = -1
  integer(c_int), save :: acc_nonodes = 0
  integer(c_int), save :: acc_noelems = 0
  integer(c_int), save :: acc_nobulkelems = 0
  integer(c_int), save :: acc_maxelemdim = 0
  integer(c_int), save :: acc_maxelemnodes = 0
  integer(c_int), save :: acc_meshdim = 0
  integer(c_int), save :: acc_novars = 0
  integer(c_int), save :: acc_cumnodes = 0
  integer(c_int), save :: acc_cumelems = 0
  integer(c_int), save :: acc_iostat = 0
  integer(c_int), save :: acc_fileunit = 1
  real(c_double), save :: acc_time = 0.0_c_double
  real(c_double), save :: acc_bbtol = 0.0_c_double
  real(c_double), save :: acc_dx = 0.0_c_double
  real(c_double), allocatable, save :: acc_bbox(:,:)
  real(c_double), save :: acc_mybbox(6) = 0.0_c_double
  logical, allocatable, save :: acc_activepart(:)
  character(len=length), save :: acc_txt = ''
  character(len=length), save :: acc_str = ''
  character(len=MAX_NAME_LEN), save :: acc_basefile = ''
  character(len=MAX_NAME_LEN), save :: acc_inputfile = ''
  character(len=MAX_NAME_LEN), save :: acc_varname = ''
  character(:), allocatable, save :: acc_inputdir

  interface
    subroutine interpolatemeshtomeshq( oldmesh, newmesh, oldvariables, newvariables, &
        usequadranttree, projector, maskname, foundnodes, newmaskperm, keepunfoundnodes )
      use types
      type(variable_t), pointer, optional :: oldvariables, newvariables
      type(mesh_t), target  :: oldmesh, newmesh
      logical, optional :: usequadranttree,foundnodes(:)
      character(len=*),optional :: maskname
      type(projector_t), pointer, optional :: projector
      integer, optional, pointer :: newmaskperm(:)
      logical, optional :: keepunfoundnodes
    end subroutine interpolatemeshtomeshq
  end interface

contains

  ! -- entry state ----------------------------------------------------
  subroutine gmr_begin(h_model, h_solver) bind(c, name="gmr_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine gmr_begin

  ! -- entry setup: Info start through bounding-box sizing ------------
  ! Returns 1 when the original takes the zero-masked-nodes RETURN path.
  function gmr_setup() result(code) bind(c, name="gmr_setup")
    integer(c_int) :: code
    integer :: i, nsize
    code = 0
    call info(caller,'Reading Gmsh results and interpolating to current mesh!')

    acc_params => getsolverparams()
    acc_tomesh => acc_solver % mesh

    acc_dim = coordinatesystemdimension()
    acc_aligncoord = listgetinteger(acc_params,'Align Coordinate',acc_found)

    acc_usequadtree = listgetlogical(acc_params,'Use Quadrant Tree',acc_found)

    acc_passivecoord = listgetinteger(acc_params,'Interpolation Passive Coordinate',acc_found)
    if(.not. acc_found) acc_passivecoord = abs(acc_aligncoord)

    acc_basefile = listgetstring(acc_params, 'Filename', unfoundfatal = .true.)
    if(index(acc_basefile,'.') == 0) write(acc_basefile,'(A,A)') trim(acc_basefile),".msh"
    call solveroutputdirectory(acc_solver, acc_basefile, acc_inputdir, usemeshdir = .true.)
    acc_basefile = trim(acc_inputdir)// '/' //trim(acc_basefile)

    acc_inputpartitions = listgetinteger(acc_params,'Filename Partitions', acc_found)
    if(.not. acc_found) acc_inputpartitions = 1

    acc_n = -1
    acc_maskperm => null()
    acc_str = listgetstring(acc_params,'Mask Name',acc_found)
    if(acc_found) then
      allocate(acc_maskperm(acc_tomesh % numberofnodes))
      acc_maskperm = 0
      call makepermusingmask(acc_model, acc_solver, acc_tomesh, acc_str, .false., &
          acc_maskperm, acc_n, requirelogical = .true.)
      call info(caller,'Using given mask "'//trim(acc_str)//'" for interpolation!')
      if(acc_n == 0) deallocate(acc_maskperm)
    else if(associated(acc_solver % variable)) then
      acc_maskperm => acc_solver % variable % perm
      nullify(acc_maskperm)
      if(associated(acc_maskperm)) then
        call info(caller,'Using Solver % Variable % Perm as the mask for interpolation!')
        acc_n = count(acc_maskperm(1:acc_tomesh % numberofnodes) > 0)
      end if
    end if

    if(acc_n == 0) then
      call info(caller,'Zero masked nodes, returning without interpolation!')
      code = 1
      return
    else if(acc_n > 0) then
      call info(caller,'Number of masked nodes: '//i2s(acc_n),level=7)
    end if

    acc_usebbox = .false.
    if(acc_inputpartitions > 1) then
      acc_usebbox = listgetlogical(acc_params,'Use Bounding Box',acc_found)
    end if

    ! bounding-box sizing (original lines 143-169) lives in gmr_bbox
    call gmr_bbox()
  end function gmr_setup

  subroutine gmr_bbox()
    integer :: i
    real(c_double) :: x
    if(.not. acc_usebbox) return
    if(allocated(acc_bbox)) deallocate(acc_bbox)
    allocate(acc_bbox(acc_inputpartitions,6))
    acc_bbox(:,1:3) = huge(x)
    acc_bbox(:,4:6) = -huge(x)
    if(allocated(acc_activepart)) deallocate(acc_activepart)
    allocate(acc_activepart(acc_inputpartitions))
    acc_activepart = .true.

    acc_n = acc_tomesh % numberofnodes
    if(associated(acc_maskperm)) then
      acc_mybbox(1) = minval(acc_tomesh % nodes % x(1:acc_n),acc_maskperm(1:acc_n)>0)
      acc_mybbox(2) = minval(acc_tomesh % nodes % y(1:acc_n),acc_maskperm(1:acc_n)>0)
      acc_mybbox(3) = minval(acc_tomesh % nodes % z(1:acc_n),acc_maskperm(1:acc_n)>0)
      acc_mybbox(4) = maxval(acc_tomesh % nodes % x(1:acc_n),acc_maskperm(1:acc_n)>0)
      acc_mybbox(5) = maxval(acc_tomesh % nodes % y(1:acc_n),acc_maskperm(1:acc_n)>0)
      acc_mybbox(6) = maxval(acc_tomesh % nodes % z(1:acc_n),acc_maskperm(1:acc_n)>0)
    else
      acc_mybbox(1) = minval(acc_tomesh % nodes % x(1:acc_n))
      acc_mybbox(2) = minval(acc_tomesh % nodes % y(1:acc_n))
      acc_mybbox(3) = minval(acc_tomesh % nodes % z(1:acc_n))
      acc_mybbox(4) = maxval(acc_tomesh % nodes % x(1:acc_n))
      acc_mybbox(5) = maxval(acc_tomesh % nodes % y(1:acc_n))
      acc_mybbox(6) = maxval(acc_tomesh % nodes % z(1:acc_n))
    end if

    acc_bbtol = listgetconstreal(acc_params,'Bounding box tolerance',acc_found)
    if(.not. acc_found) acc_bbtol = 1.0d-6
  end subroutine gmr_bbox

  ! -- sweep reset (label-10 prelude) ----------------------------------
  subroutine gmr_sweep_begin() bind(c, name="gmr_sweep_begin")
    acc_allocdone = .false.
    acc_maxelemdim = 0
    acc_maxelemnodes = 0
  end subroutine gmr_sweep_begin

  ! -- one full ReadPart sweep (label 10 body) --------------------------
  subroutine gmr_pass() bind(c, name="gmr_pass")
    integer :: i, readpart
    acc_cumnodes = 0
    acc_cumelems = 0

    do readpart = 1, acc_inputpartitions
      acc_nonodes = 0
      acc_noelems = 0

      if(acc_inputpartitions == 1) then
        acc_inputfile = trim(acc_basefile)
      else
        if(acc_allocdone .and. acc_usebbox) then
          if(.not. acc_activepart(readpart)) cycle
        end if
        i = modulo(readpart-1 + parenv % mype, acc_inputpartitions)
        acc_inputfile = trim(acc_basefile)//'_'//i2s(acc_inputpartitions)//'np'//i2s(i+1)
      end if

      if(acc_inputpartitions > 1 .or. .not. acc_allocdone) then
        call info(caller,'Reading mesh from file: '//trim(acc_inputfile),level=7)
        open(unit=acc_fileunit, file=acc_inputfile, status='old',iostat=acc_iostat)
        if(acc_iostat /= 0) then
          if(acc_inputpartitions > 1) then
            if(acc_usebbox) acc_activepart(readpart) = .false.
            call info(caller,'Could not open file: '//trim(acc_inputfile))
            cycle
          else
            call fatal(caller,'Could not open file: '//trim(acc_inputfile))
          end if
        end if
      end if

      call readsinglegmshfile(readpart)

      acc_cumnodes = acc_cumnodes + acc_nonodes
      acc_cumelems = acc_cumelems + acc_noelems

      if(acc_inputpartitions > 1 .or. acc_allocdone) then
        close(acc_fileunit)
      else
        rewind(acc_fileunit)
      end if
    end do
  end subroutine gmr_pass

  ! -- GOTO-10 query: 1 while the sizing pass still has to run ---------
  function gmr_need_alloc() result(code) bind(c, name="gmr_need_alloc")
    integer(c_int) :: code
    code = merge(1, 0, .not. acc_allocdone)
  end function gmr_need_alloc

  ! -- first-pass sizing tail (allocate FromMesh + Perm) ----------------
  subroutine gmr_allocate() bind(c, name="gmr_allocate")
    integer :: i, nsize
    if(acc_inputpartitions > 1) then
      call info(caller,'Number of cumulative nodes to read: '//i2s(acc_cumnodes),level=5)
      call info(caller,'Number of cumulative elements to read: '//i2s(acc_cumelems),level=5)
    end if
    if(acc_usebbox) then
      i = count(acc_activepart)
      call info(caller,'Number of active partitions '&
          //i2s(i)//' out of '//i2s(acc_inputpartitions),level=10)
    end if

    call info(caller,'Maximum element dimension: '//i2s(acc_maxelemdim),level=7)
    call info(caller,'Maximum element nodes: '//i2s(acc_maxelemnodes),level=7)
    acc_frommesh => allocatemesh(acc_cumelems,0,acc_cumnodes)

    acc_frommesh % meshdim = 3
    acc_frommesh % maxelementnodes = acc_maxelemnodes
    acc_allocdone = .true.

    call info(caller,'Creating variable structure',level=20)
    acc_frommesh % variables => null()
    allocate(acc_perm(acc_cumnodes))
    do i = 1, acc_cumnodes
      acc_perm(i) = i
    end do
  end subroutine gmr_allocate

  ! -- mesh-dimension fixup tail ----------------------------------------
  subroutine gmr_finalize() bind(c, name="gmr_finalize")
    call info(caller,'Last bulk element index: '//i2s(acc_nobulkelems),level=7)

    acc_frommesh % meshdim = max(acc_maxelemdim, acc_meshdim)
    acc_frommesh % numberofbulkelements = acc_nobulkelems
    acc_frommesh % numberofboundaryelements = acc_noelems - acc_nobulkelems

    call info(caller,'Gmsh reader complete!')
  end subroutine gmr_finalize

  subroutine gmr_interpolate() bind(c, name="gmr_interpolate")
    call interpolatefromgmshfile()
    call info(caller,'Interpolation from Gmsh format complete')
  end subroutine gmr_interpolate

  ! -- ReadSingleGmshFile (original CONTAINS procedure) ------------------
  subroutine readsinglegmshfile(readpart)
    integer, intent(in) :: readpart
    real(c_double) :: gmshver
    real(c_double) :: coord(3)
    integer, save :: gmshtoelmertype(21) = (/ 202, 303, 404, 504, 808, 706, 605, 203, 306, 409, &
        510, 827, 0, 0, 101, 408, 820, 715, 613, 0, 310 /)
    integer :: gmshindexes(27)
    integer :: i, j, k, n, elmir, gmshtype, elmdim
    save gmshver

    acc_novars = 0
    acc_meshdim = 0

    do while(.true.)
      read(acc_fileunit,'(A)',end=20,err=20) acc_str
      if(seql(acc_str, '$MeshFormat')) then
        read(acc_fileunit,'(A)',end=20,err=20) acc_str
        if(.not. acc_allocdone) then
          read(acc_str,*) gmshver
          write(message,'(A,ES12.3)') 'Gmsh file version: ',gmshver
          call info(caller,message)
        end if
      end if

      if(seql(acc_str, '$Nodes')) then
        read(acc_fileunit,'(A)',end=20,err=20) acc_str
        read(acc_str,*) acc_nonodes
        if(.not. acc_allocdone) then
          call info(caller,'Number of nodes in mesh: '//i2s(acc_nonodes),level=7)
        end if
        do i = 1, acc_nonodes
          read(acc_fileunit,'(A)',end=20,err=20) acc_str
          read(acc_str,*) j, coord
          if(i /= j) call fatal(caller,'Do node permutations!')
          if(acc_allocdone) then
            k = i + acc_cumnodes
            acc_frommesh % nodes % x(k) = coord(1)
            acc_frommesh % nodes % y(k) = coord(2)
            acc_frommesh % nodes % z(k) = coord(3)
          else if(acc_usebbox) then
            acc_bbox(readpart,1:3) = min(acc_bbox(readpart,1:3),coord)
            acc_bbox(readpart,4:6) = max(acc_bbox(readpart,4:6),coord)
          end if
        end do
        read(acc_fileunit,'(A)',end=20,err=20) acc_str

        if(.not. acc_allocdone) then
          if(acc_usebbox) then
            do i = 1, 3
              if(i == acc_passivecoord) cycle
              if(i > acc_dim) cycle
              if(acc_bbox(readpart,i) > acc_mybbox(i+3) + acc_bbtol) &
                  acc_activepart(readpart) = .false.
              if(acc_bbox(readpart,i+3) < acc_mybbox(i) - acc_bbtol) &
                  acc_activepart(readpart) = .false.
            end do
            if(.not. acc_activepart(readpart)) then
              acc_nonodes = 0
              acc_noelems = 0
              exit
            end if
          end if

          if(acc_noelems > 0 .and. acc_nonodes > 0) exit
        end if
      end if

      if(seql(acc_str, '$Elements')) then
        read(acc_fileunit,'(A)',end=20,err=20) acc_str
        read(acc_str,*) acc_noelems
        if(.not. acc_allocdone) then
          call info(caller,'Number of elements in mesh: '//i2s(acc_noelems),level=7)
        end if

        do i = 1, acc_noelems
          read(acc_fileunit,'(A)',end=20,err=20) acc_str
          read(acc_str,*) j, gmshtype
          if(i /= j) call fatal(caller,'Do element permutations!')
          elmir = gmshtoelmertype(gmshtype)

          elmdim = 0
          if(elmir > 500) then
            elmdim = 3
          else if(elmir > 300) then
            elmdim = 2
          else if(elmir > 200) then
            elmdim = 1
          end if

          if(acc_allocdone) then
            k = acc_cumelems + i
            acc_element => acc_frommesh % elements(k)
            acc_element % type => getelementtype(elmir)
            n = acc_element % type % numberofnodes
            call allocatevector(acc_element % nodeindexes, n)
            read(acc_str,*) j, gmshtype, j, j, j, acc_element % nodeindexes(1:n)
            if(acc_cumnodes > 0) then
              acc_element % nodeindexes(1:n) = acc_element % nodeindexes(1:n) + acc_cumnodes
            end if
            if(elmdim == acc_maxelemdim) acc_nobulkelems = k
          else
            acc_maxelemdim = max(acc_maxelemdim, elmdim)
            acc_maxelemnodes = max(acc_maxelemnodes, modulo(elmir, 100))
          end if
        end do
        read(acc_fileunit,'(A)',end=20,err=20) acc_str

        if(.not. acc_allocdone) then
          if(acc_noelems > 0 .and. acc_nonodes > 0) exit
        end if
      end if

      if(seql(acc_str, '$NodeData')) then
        read(acc_fileunit,'(A)',end=20,err=20) acc_str
        read(acc_fileunit,'(A)',end=20,err=20) acc_str
        n = len_trim(acc_str)

        acc_varname = acc_str(2:n-1)
        call info(caller,'Reading gmsh variable: '//trim(acc_varname))

        read(acc_fileunit,'(A)',end=20,err=20) acc_str
        read(acc_fileunit,'(A)',end=20,err=20) acc_str
        read(acc_str, *) acc_time
        read(acc_fileunit,'(A)',end=20,err=20) acc_str
        read(acc_str, *) acc_meshdim
        read(acc_fileunit,'(A)',end=20,err=20) acc_str
        read(acc_fileunit,'(A)',end=20,err=20) acc_str
        read(acc_str, *) acc_dofs
        read(acc_fileunit,'(A)',end=20,err=20) acc_str
        read(acc_str, *) n

        if(n /= acc_nonodes) then
          call info(caller,'Mismatch is vector size!')
        end if

        acc_var => null()
        if(associated(acc_frommesh % variables)) then
          call info(caller,'Reusing variable structure',level=20)
          acc_var => variableget(acc_frommesh % variables, acc_varname, thisonly = .true.)
        end if

        if(.not. associated(acc_var)) then
          call variableaddvector(acc_frommesh % variables, acc_frommesh, acc_solver, &
              acc_varname, dofs = acc_dofs, perm = acc_perm)
          acc_var => variableget(acc_frommesh % variables, acc_varname, thisonly = .true.)
        end if

        do i = 1, acc_nonodes
          read(acc_fileunit,'(A)',end=20,err=20) acc_str
          k = i + acc_cumnodes
          if(acc_dofs == 1) then
            read(acc_str, *) j, acc_var % values(k)
          else
            read(acc_str, *) j, acc_var % values(acc_dofs*(k-1)+1:acc_dofs*k)
          end if
        end do
        read(acc_fileunit,'(A)',end=20,err=20) acc_str

        acc_novars = acc_novars + 1
        acc_vararray(acc_novars) % var => acc_var

        if(infoactive(28)) then
          if(acc_novars == 1) call info(caller,'Initial field ranges:')
          call vectorvaluesrange(acc_var % values,size(acc_var % values),'From: '//trim(acc_var % name))
        end if
      end if
    end do

20  continue
  end subroutine readsinglegmshfile

  ! -- InterpolateFromGmshFile (original CONTAINS procedure) --------------
  subroutine interpolatefromgmshfile()
    real(c_double), pointer :: x1(:), x2(:)
    real(c_double) :: minx, maxx
    integer :: n1, n2, i, k

    if(acc_aligncoord /= 0) then
      k = abs(acc_aligncoord)

      if(k == 1) then
        x1 => acc_frommesh % nodes % x
        x2 => acc_tomesh % nodes % x
      else if(k == 2) then
        x1 => acc_frommesh % nodes % y
        x2 => acc_tomesh % nodes % y
      else if(k == 3) then
        x1 => acc_frommesh % nodes % z
        x2 => acc_tomesh % nodes % z
      else
        call fatal(caller,'Invalid value for "Align Coordinate": '//i2s(acc_aligncoord))
      end if

      n1 = acc_frommesh % numberofnodes
      n2 = acc_tomesh % numberofnodes

      if(acc_aligncoord > 0) then
        minx = minval(x2(1:n2))
        maxx = maxval(x1(1:n1))
        acc_dx = minx - maxx
      else
        minx = minval(x1(1:n1))
        maxx = maxval(x2(1:n2))
        acc_dx = minx - maxx
      end if

      write(message,'(A,ES12.3)') 'Aligning coordinate '//i2s(k)//' with: ',acc_dx
      call info(caller,message)

      x1 = x1 + acc_dx
    end if

    if(associated(acc_maskperm)) then
      call interpolatemeshtomeshq(acc_frommesh, acc_tomesh, acc_frommesh % variables, &
          acc_tomesh % variables, usequadranttree=acc_usequadtree,newmaskperm=acc_maskperm)
    else
      call interpolatemeshtomeshq(acc_frommesh, acc_tomesh, acc_frommesh % variables, &
          acc_tomesh % variables, usequadranttree=acc_usequadtree)
    end if

    if(infoactive(28)) then
      call info(caller,'Projected field ranges:')
      do i = 1, acc_novars
        acc_var => variableget(acc_tomesh % variables, acc_vararray(i) % var % name, &
            thisonly = .true.)
        if(associated(acc_var)) then
          call vectorvaluesrange(acc_var % values,size(acc_var % values),'To: '//trim(acc_var % name))
        end if
      end do
    end if
  end subroutine interpolatefromgmshfile

  ! -- GmshOutputReader_Init body -----------------------------------------
  subroutine gmr_init_begin(h_solver) bind(c, name="gmr_init_begin")
    type(c_ptr), value :: h_solver
    call c_f_pointer(h_solver, acc_solver)
  end subroutine gmr_init_begin

  subroutine gmr_init_body() bind(c, name="gmr_init_body")
    type(valuelist_t), pointer :: params
    params => getsolverparams()
    call listaddnewlogical(params,'No Matrix',.true.)
  end subroutine gmr_init_body

end module gmshoutputreader_acc
