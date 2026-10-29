! T021-B PartitionMesh port: statement-level BIND(C) wrappers for the Elmer
! calls in fem/src/modules/PartitionMesh.F90, compiled against the built
! .mod files of the target configure (T020 Poisson mechanism, ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars). All Fortran-object traffic — TYPE(...)
! handles, OPTIONAL/POINTER dummies, assumed-shape descriptors, hidden
! lengths, array-pointer results — stays inside Fortran, so the Rust side
! is pure C calls. The wrapper module mirrors the original's host-
! associated state (Solver/Model/Mesh/Params/ElementPart/NeighbourList
! plus the Output Directory scratch) with identical lifetimes: the
! original's CONTAINS procedures capture host locals per invocation, so
! the wrappers reset that state at each entry exactly as a fresh call
! would (deallocate-if-associated before the fresh ALLOCATE).
!
! Clean-room code: statement-for-statement wrappers around Elmer library
! calls, no upstream text.
module partitionmesh_acc
  use iso_c_binding
  use types
  use defutils
  use meshpartition
  use meshio
  use messages
  use spariterglobals, only: parenv
  implicit none

  type(solver_t), pointer, save :: pma_solver => null()
  type(model_t), pointer, save :: pma_model => null()
  type(mesh_t), pointer, save :: pma_mesh => null()
  type(valuelist_t), pointer, save :: pma_params => null()
  integer, pointer, save :: pma_elempart(:) => null()
  type(neighbourlist_t), pointer, save :: pma_neigh(:) => null()
  character(len=MAX_NAME_LEN), save :: pma_dirname = ' '

contains

  ! -- entry state ----------------------------------------------------
  subroutine pma_begin(h_solver, h_model) bind(c, name="pma_begin")
    type(c_ptr), value :: h_solver, h_model
    call c_f_pointer(h_solver, pma_solver)
    call c_f_pointer(h_model, pma_model)
    pma_params => pma_solver % values
    pma_mesh => pma_solver % mesh
  end subroutine pma_begin

  subroutine pma_info_start() bind(c, name="pma_info_start")
    call Info('PartitionMesh','Using internal mesh partitioning')
  end subroutine pma_info_start

  function pma_pes() result(n) bind(c, name="pma_pes")
    integer(c_int) :: n
    n = parenv % pes
  end function pma_pes

  subroutine pma_info_serial_only() bind(c, name="pma_info_serial_only")
    call Info('PartitionMes','Currently the mesh can only be partitioned in serial!')
  end subroutine pma_info_serial_only

  ! -- partitioning ---------------------------------------------------
  subroutine pma_partition() bind(c, name="pma_partition")
    call PartitionMeshSerial( pma_model, pma_mesh, pma_params )
    pma_elempart => pma_mesh % repartition
  end subroutine pma_partition

  ! Exact transcription of the CreateNeighbourList internal procedure
  ! (order-sensitive loops, allocation-failure Fatal calls, and the
  ! formatted I2S/Message Info tail stay in Fortran by construction).
  subroutine pma_create_neighbour_list() bind(c, name="pma_create_neighbour_list")
    integer :: i,j,k,l,n,m,partition,lsum,lmax
    integer :: TmpNeighbours(100)
    type(element_t), pointer :: element
    integer :: allocstat

    call Info('PartitionMesh','Creating neighbour list for parallel saving')

    ! Fresh-entry semantics: the original's host-local NeighbourList
    ! starts unassociated on every call; release any previous SAVE state.
    if (associated(pma_neigh)) deallocate(pma_neigh)

    n = pma_mesh % numberofnodes
    allocate( pma_neigh(n) , stat=allocstat )
    if( allocstat /= 0 ) then
      call Fatal('PartitionMesh','Allocation error for NeighbourList')
    end if

    do i=1,n
      nullify( pma_neigh(i) % neighbours )
    end do

    do i=1,pma_mesh % numberofbulkelements
      element => pma_mesh % elements(i)
      partition = pma_elempart(i)
      m = element % type % numberofnodes
      do j=1,m
        k = element % nodeindexes(j)
        if( .not. associated( pma_neigh(k) % neighbours ) ) then
          allocate( pma_neigh(k) % neighbours(1), stat = allocstat )
          if( allocstat /= 0 ) then
            call Fatal('PartitionMesh','Allocation error for Neighbours')
          end if
          pma_neigh(k) % neighbours(1) = partition
        else if( .not. any( pma_neigh(k) % neighbours == partition ) ) then
          l = size( pma_neigh(k) % neighbours )

          TmpNeighbours(1:l) = pma_neigh(k) % neighbours(1:l)
          deallocate( pma_neigh(k) % neighbours )

          allocate( pma_neigh(k) % neighbours(l+1), stat = allocstat )
          if( allocstat /= 0 ) then
            call Fatal('PartitionMesh','Allocation error for Neighbours')
          end if
          pma_neigh(k) % neighbours(1:l) = TmpNeighbours(1:l)
          pma_neigh(k) % neighbours(l+1) = partition
        end if
      end do
    end do

    lmax = 0
    lsum = 0
    do k=1,n
      l = size(pma_neigh(k) % neighbours)
      lmax = max( lmax, l )
      lsum = lsum + l
    end do

    call Info('PartitionMesh','Maximum number of partitions for a node: '//I2S(lmax))

    write(Message,'(A,F8.3)') 'Average number of partitions for a node: ',1.0_dp*lsum/n
    call Info('PartitionMesh',Message)
  end subroutine pma_create_neighbour_list

  ! -- output-directory branch (CHARACTER stays in Fortran) ------------
  function pma_output_found() result(r) bind(c, name="pma_output_found")
    integer(c_int) :: r
    logical :: found
    pma_dirname = ListGetString( pma_params,'Output Directory',found )
    r = merge(1, 0, found)
  end function pma_output_found

  subroutine pma_write_partitioned() bind(c, name="pma_write_partitioned")
    call WriteMeshToDiskPartitioned( pma_model, pma_mesh, pma_dirname, &
        pma_elempart, pma_neigh )
  end subroutine pma_write_partitioned

  ! -- partitioning-exit branch ---------------------------------------
  function pma_exit_flag() result(r) bind(c, name="pma_exit_flag")
    integer(c_int) :: r
    logical :: found
    r = merge(1, 0, ListGetLogical( pma_params,'Partitioning Exit',found ))
  end function pma_exit_flag

  subroutine pma_stop() bind(c, name="pma_stop")
    stop
  end subroutine pma_stop

  ! -- partition variable ---------------------------------------------
  subroutine pma_info_create_var() bind(c, name="pma_info_create_var")
    call Info('PartitionMesh','Create partitioning variable')
  end subroutine pma_info_create_var

  ! Exact transcription of the SetPartitionVariable internal procedure.
  subroutine pma_set_partition_variable() bind(c, name="pma_set_partition_variable")
    type(variable_t), pointer :: var
    integer :: t
    type(element_t), pointer :: element
    integer :: n

    var => pma_solver % variable
    if( .not. associated( var ) ) return

    call Info('PartitionMesh','Setting discontinuous field > Partition < ')

    var % values = 0.0_dp

    n = pma_mesh % numberofbulkelements
    do t=1,n
      element => pma_solver % mesh % elements(t)
      if( associated( element % dgindexes ) ) then
        var % values( element % dgindexes ) = 1.0_dp * pma_elempart(t)
      else
        var % values( element % nodeindexes ) = 1.0_dp * pma_elempart(t)
      end if
    end do
  end subroutine pma_set_partition_variable

  subroutine pma_info_done() bind(c, name="pma_info_done")
    call Info('PartitionMesh','All done for now')
  end subroutine pma_info_done

end module partitionmesh_acc
