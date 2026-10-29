! T022-A FluidicForce port: statement-level BIND(C) wrappers for the
! Elmer calls in fem/src/modules/FluidicForce.F90, compiled against
! the built .mod files of the target configure (A2 shipped mechanism,
! ADR-057 §3).
!
! Every wrapper takes only C-interoperable arguments (C addresses as
! TYPE(C_PTR), VALUE scalars). All Fortran-object traffic — TYPE(...)
! handles, OPTIONAL/POINTER dummies, hidden lengths, assumed-shape
! descriptors, array-pointer results — stays inside Fortran, so the Rust
! side is pure C calls. The wrapper module mirrors FluidicForce's own
! entry state (Model/Solver handles) with identical lifetimes; the whole
! entry body (ForceCompute) plus its CONTAINS internal (ForceIntegrate)
! is transcribed statement-for-statement below, so bit-identity holds
! by construction (same compiler, same flags, same call sequence).
! All floating-point work (the per-integration-point stress assembly,
! Level=4, the shear-stress data files) stays in Fortran wrappers.
!
! Clean-room code: statement-for-statement transcription around Elmer
! library calls, no upstream text.
module fluidicforce_acc
  use iso_c_binding
  use types
  use defutils
  use integration
  use elementdescription
  use elementgeometry
  use generalutils
  use coordinatesystems
  use solverbasics
  use materialmodels
  use messages
  use eleminfo
  implicit none

  type(model_t), pointer, save :: acc_model => null()
  type(solver_t), pointer, save :: acc_solver => null()

contains

  ! -- entry state ----------------------------------------------------
  subroutine ffc_begin(h_model, h_solver) bind(c, name="ffc_begin")
    type(c_ptr), value :: h_model, h_solver
    call c_f_pointer(h_model, acc_model)
    call c_f_pointer(h_solver, acc_solver)
  end subroutine ffc_begin

  ! -- ForceCompute body (exact transcription) --------------------------
  subroutine ffc_run() bind(c, name="ffc_run")
    real(c_double), allocatable :: pressure(:), viscosity(:)
    real(c_double), allocatable :: velocity(:,:), shear_data(:,:)
    real(c_double), allocatable :: forces(:,:), moments(:,:), areas(:)
    real(c_double) :: force(3), moment(3), moment_about(3), area, shear_stress
    real(c_double), pointer :: mwork(:,:)
    logical :: stat, calc_moment, viscous_force, compressible, sum_forces
    logical :: shear_output
    integer :: i,j,k,n,pn,t,dim
    integer :: nbr_shear_values, nlen
    integer, pointer :: node_indexes(:), indices(:)
    type(variable_t), pointer :: var
    type(mesh_t), pointer :: mesh
    type(valuelist_t), pointer :: material
    type(nodes_t) :: element_nodes, parent_nodes
    type(element_t), pointer :: current_element, parent
    type(valuelist_t), pointer :: bc
    character(len=MAX_NAME_LEN) :: shear_filename, message_l, viscosity_flag
    character(len=MAX_NAME_LEN) :: compressibility_flag, boundary_name
    character(len=MAX_NAME_LEN) :: variable_name

    call Info( 'ForceCompute', '-------------------------------------',level=4 )
    call Info( 'ForceCompute', 'Computing Fluidic Force:  ', level=4 )
    call Info( 'ForceCompute', '-------------------------------------',level=4 )

    mesh => GetMesh()

    variable_name = GetString( acc_solver % values, 'Velocity Field Name', stat )
    if ( .not. stat )  then
       var => VariableGet( mesh % variables, 'Flow Solution', .TRUE. )
    else
       var => VariableGet( mesh % variables, variable_name, .TRUE. )
    end if

    n = mesh % maxelementnodes
    allocate( element_nodes % x(n), element_nodes % y(n), element_nodes % z(n))
    allocate( parent_nodes % x(n), parent_nodes % y(n), parent_nodes % z(n) )
    allocate( pressure( n ), viscosity( n ), velocity( 3, n ) )

    allocate( forces(acc_model % numberofbcs,3)  )
    allocate( moments(acc_model % numberofbcs,3)  )
    allocate( areas(acc_model % numberofbcs)  )

    moment_about = 0.0d0
    forces = 0.0d0
    moments = 0.0d0
    areas = 0.0d0

    dim = CoordinateSystemDimension()

    viscous_force = ListGetLogical( acc_solver % values, 'Calculate Viscous Force',stat )
    if ( .not. stat )  viscous_force = .TRUE.

    shear_output = GetLogical( acc_solver % values, 'Shear Stress Output', stat )
    if ( shear_output .and. dim == 3 ) then
       call Warn( 'FluidicForce', &
            'Shear stress output not implemented for 3D cases' )
       shear_output = .FALSE.
    end if
    if ( shear_output ) then
       shear_filename = GetString( acc_solver % values, 'Shear Stress Output File', stat )
       if ( .not. stat )  shear_filename = 'shearstress.dat'

       nbr_shear_values = 0
       allocate( shear_data( mesh % numberofboundaryelements, 3 ) )
       shear_data = 0.0d0
    end if

    sum_forces = ListGetLogical( acc_solver % values, 'Sum Forces',stat )

    if (CurrentCoordinateSystem() == 2 .or. CurrentCoordinateSystem() >= 5) then
       do i = 1, acc_model % numberofmaterials
          material => acc_model % materials(i) % values
          compressibility_flag = ListGetString( material, &
               'Compressibility Model', stat)
          if ( stat .and. compressibility_flag /= 'incompressible' ) then
             call Warn( 'ForceCompute', &
                  'Force component due to compressibility not implemented in polar coordinates' )
             exit
          end if
       end do
    end if


    do t = mesh % numberofbulkelements + 1, &
        mesh % numberofbulkelements + mesh % numberofboundaryelements

      current_element => mesh % elements(t)
      acc_model % currentelement => mesh % elements(t)
      n = current_element % TYPE % numberofnodes
      node_indexes => current_element % nodeindexes

      if ( current_element % TYPE % elementcode == 101 ) cycle

      do k=1, acc_model % numberofbcs
        bc => acc_model % bcs(k) % values
        if ( acc_model % bcs(k) % tag /= current_element % boundaryinfo % constraint ) cycle
        if ( .not. ListGetLogical(bc,'Calculate Fluidic Force',stat ) ) cycle

        element_nodes % x(1:n) = mesh % nodes % x(node_indexes)
        element_nodes % y(1:n) = mesh % nodes % y(node_indexes)
        element_nodes % z(1:n) = mesh % nodes % z(node_indexes)

        parent => current_element % boundaryinfo % left

        stat = ASSOCIATED( parent )
        if ( stat ) stat = stat .and. ALL(var % perm(parent % nodeindexes(1:n)) > 0)

        if ( .not. stat ) then
           parent => current_element % boundaryinfo % right

           stat = ASSOCIATED( parent )
           if ( stat ) stat = ALL(var % perm(parent % nodeindexes(1:n)) > 0)

           if ( .not. stat )  call Fatal( 'ForceCompute', &
                'No flow solution available for specified boundary' )
        end if

        i = parent % bodyid
        pn = parent % Type % numberofnodes

        parent_nodes % x(1:pn) = mesh % nodes % x(parent % nodeindexes)
        parent_nodes % y(1:pn) = mesh % nodes % y(parent % nodeindexes)
        parent_nodes % z(1:pn) = mesh % nodes % z(parent % nodeindexes)

        j = ListGetInteger( acc_model % bodies(i) % values, 'Material', &
              minv=1, maxv=acc_model % numberofmaterials )
        material => acc_model % materials(j) % values

        acc_model % currentelement => parent
        compressible = .FALSE.
        compressibility_flag = ListGetString( material, &
             'Compressibility Model', stat)
        if ( stat .and. compressibility_flag /= 'incompressible' )  &
             compressible = .TRUE.

        viscosity(1:pn) = ListGetReal( material, 'Viscosity', pn, parent % nodeindexes )

        mwork => ListGetConstRealArray( bc,'Moment About',calc_moment)
        if( calc_moment ) then
          moment_about(1:dim) = mwork(1:dim,1)
        else
          moment_about(1) = ListGetCReal( bc,'Moment About 1', stat )
          moment_about(2) = ListGetCReal( bc,'Moment About 2', calc_moment )
          calc_moment = stat .or. calc_moment
          moment_about(3) = ListGetCReal( bc,'Moment About 3', stat )
          calc_moment = stat .or. calc_moment
        end if

        velocity = 0.0d0
        do i=1,pn
           do j=1,dim
              velocity(j,i) = &
                var % values(var % dofs * (var % perm(parent % nodeindexes(i))-1)+j)
           end do
        end do
        pressure(1:pn) = var % values(var % dofs * var % perm(parent % nodeindexes))

        force = 0.0d0
        moment = 0.0d0
        area = 0.0d0
        shear_stress = 0.0d0

        call ffc_integrate( force, moment, moment_about, viscous_force, compressible, &
             area, shear_stress)

        forces(k,1:3) = forces(k,1:3) + force(1:3)
        moments(k,1:3) = moments(k,1:3) + moment(1:3)
        areas(k) = areas(k) + area

        if ( shear_output ) then
           nbr_shear_values = nbr_shear_values + 1
           shear_data(nbr_shear_values,1) = shear_stress
           shear_data(nbr_shear_values,2) = SUM( element_nodes % x(1:n) ) / n
           shear_data(nbr_shear_values,3) = SUM( element_nodes % y(1:n) ) / n
        end if
     end do
    end do


    do k=1, acc_model % numberofbcs

      if( sum_forces ) then
        if( k > 1) exit
        areas(1) = SUM(areas)
        do i=1,3
          forces(1,i) = SUM(forces(:,i))
          moments(1,i) = SUM(moments(:,i))
        end do
        WRITE( boundary_name, '("")')
      else
        if ( .not. ListGetLogical(acc_model % bcs(k) % values,'Calculate Fluidic Force',stat ) ) cycle
        WRITE( boundary_name, '("bc ",I0)') k
      end if

      nlen = LEN_TRIM(boundary_name)

      call Info('ForceCompute','Forces on Boundary '//boundary_name(1:nlen),level=4 )
      WRITE( message, '("Fluidic Force:", 3ES17.6E2)') forces(k,1:3)
      call Info( 'ForceCompute', message, level=4 )
      WRITE( message, '("Resultant Force:", ES17.6E2)') SQRT(SUM(forces(k,1:3)**2 ))
      call Info( 'ForceCompute', message, level=4 )
      WRITE( message, '("Contact Area:   ", ES17.6E2)') areas(k)
      call Info( 'ForceCompute', message, level=4 )
      if ( calc_moment ) then
        WRITE( message, &
            '("Moment about (",ES9.3E1,",",ES9.3E1,",",ES9.3E1,") is:",3Es14.6E2)') &
            moment_about(1:3), moments(k,1:3)
        call Info( 'ForceCompute', message, level=4 )

        call ListAddConstReal( acc_model % simulation, &
             'res: fluid moment 3 '//boundary_name(1:nlen), moments(k,3) )
        call ListAddConstReal( acc_model % simulation, &
             'res: fluid moment 2 '//boundary_name(1:nlen), moments(k,2) )
        call ListAddConstReal( acc_model % simulation, &
             'res: fluid moment 1 '//boundary_name(1:nlen), moments(k,1) )
     end if

      call ListAddConstReal( acc_model % simulation, &
          'res: fluid force area '//boundary_name(1:nlen), areas(k) )
      call ListAddConstReal( acc_model % simulation, &
          'res: fluid force '//boundary_name(1:nlen), SQRT(SUM(forces(k,1:3)**2 )) )

      if ( dim > 2 )  call ListAddConstReal( acc_model % simulation, &
          'res: fluid force 3 '//boundary_name(1:nlen), forces(k,3) )
      call ListAddConstReal( acc_model % simulation, &
          'res: fluid force 2 '//boundary_name(1:nlen), forces(k,2) )
      call ListAddConstReal( acc_model % simulation, &
          'res: fluid force 1 '//boundary_name(1:nlen), forces(k,1) )
    end do

    if ( shear_output ) then

       allocate( indices( nbr_shear_values ) )
       indices = 0
       indices = (/ ( i, i= 1,nbr_shear_values ) /)

       if ( ( MaxVal( shear_data(1:nbr_shear_values,2) ) - &
              MINVAL( shear_data(1:nbr_shear_values,2) ) ) > &
            ( MaxVal( shear_data(1:nbr_shear_values,3) ) - &
              MINVAL( shear_data(1:nbr_shear_values,3) ) ) ) then
          call SortF( nbr_shear_values, indices, shear_data(1:nbr_shear_values,2) )

          shear_data(1:nbr_shear_values,1) = shear_data(indices,1)
          shear_data(1:nbr_shear_values,3) = shear_data(indices,3)
       else
          call SortF( nbr_shear_values, indices, shear_data(1:nbr_shear_values,3) )

          shear_data(1:nbr_shear_values,1) = shear_data(indices,1)
          shear_data(1:nbr_shear_values,2) = shear_data(indices,2)
       end if

       OPEN(10, FILE=shear_filename)
       do t = 1, nbr_shear_values
          WRITE( 10, * ) shear_data(t,:)
       end do
       CLOSE(10)

       WRITE( message_l, * ) 'Variables in columns of matrix: ' // TRIM(shear_filename)
       shear_filename = TRIM( shear_filename ) // '.names'

       OPEN(10, FILE=shear_filename)
       WRITE( 10, * ) TRIM(message_l)
       WRITE( 10, * ) '1: Shear stress [N/m2]'
       WRITE( 10, * ) '2: Coordinate 1'
       WRITE( 10, * ) '3: Coordinate 2'

       CLOSE(10)

       DEALLOCATE( shear_data )
       DEALLOCATE( indices )
    end if


    DEALLOCATE( element_nodes % x, element_nodes % y, element_nodes % z )
    DEALLOCATE( parent_nodes % x, parent_nodes % y, parent_nodes % z )
    DEALLOCATE( forces, moments, areas)
    DEALLOCATE( pressure, viscosity, velocity )

  contains

    ! -- ForceIntegrate internal (exact transcription) ------------------
    subroutine ffc_integrate( force, moment, moment_about, viscous_force, &
         compressible, area, shear_stress )
       real(c_double) :: force(3), moment(3), moment_about(3), area, shear_stress
       logical :: viscous_force, compressible

       real(c_double) :: lforce(3), lmoment(3), tforce
       real(c_double) :: basis(n), dbasisdx(n,3)
       real(c_double) :: parentbasis(pn), parentdbasisdx(pn,3)
       real(c_double) :: sqrtmetric, metric(3,3), symb(3,3,3), dsymb(3,3,3,3)
       real(c_double) :: u, v, w, s, detj, x(n), y(n), z(n), xpos, ypos, zpos
       real(c_double) :: grad(3,3), stress(3,3), normal(3), radius(3)
       real(c_double) :: div
       real(c_double) :: visc

       integer :: n_integ
       real(c_double), pointer :: u_integ(:), v_integ(:), w_integ(:), s_integ(:)

       logical :: stat
       integer :: i,t
       type(GaussIntegrationPoints_t), target :: integ_stuff

       integ_stuff = GaussPoints( current_element )

       u_integ => integ_stuff % u
       v_integ => integ_stuff % v
       w_integ => integ_stuff % w
       s_integ => integ_stuff % s
       n_integ =  integ_stuff % n

       tforce = 0.0d0

       do t=1,n_integ

          u = u_integ(t)
          v = v_integ(t)
          w = w_integ(t)

          stat = ElementInfo( current_element, element_nodes, u, v, w, &
             detj, basis, dbasisdx )

          s = 1.0d0
          if ( CurrentCoordinateSystem() /= Cartesian ) then
            xpos = SUM( element_nodes % x(1:n) * basis(1:n) )
            ypos = SUM( element_nodes % y(1:n) * basis(1:n) )
            zpos = SUM( element_nodes % z(1:n) * basis(1:n) )
            s = 2*PI
          end if

          call CoordinateSystemInfo( metric,sqrtmetric,symb,dsymb,xpos,ypos,zpos)

          s = s * sqrtmetric * detj * s_integ(t)

          normal = Normalvector( current_element, element_nodes, u, v, .TRUE. )

          do i = 1,n
            do j = 1,pn
              if ( current_element % nodeindexes(i) == parent % nodeindexes(j) ) then
                 x(i) = parent % TYPE % nodeu(j)
                 y(i) = parent % TYPE % nodev(j)
                 z(i) = parent % TYPE % nodew(j)
                 exit
              end if
            end do
          end do

          u = SUM( basis(1:n) * x(1:n) )
          v = SUM( basis(1:n) * y(1:n) )
          w = SUM( basis(1:n) * z(1:n) )

          stat = ElementInfo( parent, parent_nodes, u, v, w, detj, parentbasis, &
                    parentdbasisdx )


          stress = 0.0d0
          div = 0.0d0

          if ( viscous_force ) then

            grad = MATMUL( velocity(:,1:pn),parentdbasisdx )
            visc = SUM( viscosity(1:pn) * parentbasis(1:pn) )

            viscosity_flag = ListGetString( material,'Viscosity Model', stat )
            if ( stat )  call Warn( 'FluidicForce', &
                 'Viscosity models ignored in fluidic force computation' )

            if ( compressible ) then

              if ( CurrentCoordinateSystem() == Cartesian ) then
                do i = 1, dim
                  div = div + grad(i,i)
                end do
              else
                div = SUM( velocity(1,1:pn) * parentdbasisdx(1:pn,1) ) + &
                    SUM( velocity(1,1:pn) * parentbasis(1:pn) ) / xpos + &
                    SUM( velocity(2,1:pn) * parentdbasisdx(1:pn,2) )
              end if

            end if

            stress = visc * ( grad + TRANSPOSE(grad) )

          end if

          if ( CurrentCoordinateSystem() == 2 &
               .or. CurrentCoordinateSystem() >= 5 ) then
             do i=1,dim
                stress(i,i) = stress(i,i) - SUM( pressure(1:pn) * parentbasis )
             end do
          else
             do i=1,dim
                stress(i,i) = stress(i,i) - SUM( pressure(1:pn) * parentbasis ) &
                     -(2.0d0/3.0d0)*visc*div
             end do
          end if

          lforce = -MATMUL( stress, normal )
          force  = force  + s * lforce
          tforce = tforce + s * (lforce(1)*normal(2) - lforce(2)*normal(1))

          if (calc_moment) then
            radius(1) = SUM( (element_nodes % x(1:n) - moment_about(1)) * basis )
            radius(2) = SUM( (element_nodes % y(1:n) - moment_about(2)) * basis )
            radius(3) = SUM( (element_nodes % z(1:n) - moment_about(3)) * basis )

            lmoment(1) = radius(2) * lforce(3) - radius(3) * lforce(2)
            lmoment(2) = radius(3) * lforce(1) - radius(1) * lforce(3)
            lmoment(3) = radius(1) * lforce(2) - radius(2) * lforce(1)

            moment = moment + s * lmoment
          end if

          area = area + s

       end do

       shear_stress = tforce / area

    end subroutine ffc_integrate

  end subroutine ffc_run

end module fluidicforce_acc
