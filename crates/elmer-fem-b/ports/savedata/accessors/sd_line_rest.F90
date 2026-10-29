! SaveData pilot, SaveLine phases part 2: SavePolylines (both branches).
! Shares state with sd_line_acc via USE association.
module sd_line_rest
  use iso_c_binding
  use defutils
  use sd_line_acc
  implicit none

contains

  ! -- SavePolyLines (verbatim phase) ------------------------------------------
  subroutine sdl_save_polylines() bind(c, name="sdl_save_polylines")
    type(solver_t), pointer :: psolver
    real(kind=dp) :: linepos, tanprod(2), s, eps
    real(kind=dp) :: r0(3), r1(3), dr(3), s1(3), ds(3), localcoord(3)
    real(kind=dp) :: mincoord(3), maxcoord(3), globalcoord(3)
    real(kind=dp) :: linen(3), linet1(3), linet2(3)
    real(kind=dp) :: detj
    integer :: line, noline, nodims, nsize, i, j, t, n, m, imin, imax
    integer :: nofaces
    logical :: gotit, gotdiv
    integer, pointer :: nodeindexes(:)
    psolver => acc_solver
    eps = 1.0e-5_dp
    acc_saveaxis = .false.
    acc_saveaxis(1) = ListGetLogical(acc_params, 'Save Axis', gotit)
    if (gotit) then
      acc_saveaxis(2:3) = acc_saveaxis(1)
    else
      acc_saveaxis(1) = ListGetLogical(acc_params, 'Save Axis 1', gotit)
      acc_saveaxis(2) = ListGetLogical(acc_params, 'Save Axis 2', gotit)
      acc_saveaxis(3) = ListGetLogical(acc_params, 'Save Axis 3', gotit)
    end if
    acc_noaxis = acc_dim
    if (associated(acc_pointcoordinates)) nullify(acc_pointcoordinates)
    acc_pointcoordinates => ListGetConstRealArray(acc_params, &
        'Polyline Coordinates', gotit)
    if (gotit) then
      noline = size(acc_pointcoordinates, 1) / 2
      nodims = size(acc_pointcoordinates, 2)
      if (nodims < acc_meshdim) then
        call Warn(caller, 'Dimension of points smaller than that of mesh')
      end if
    else
      noline = 0
    end if
    acc_nolines = noline
    acc_nodims = nodims
    call Info(caller, 'Saving PolyLines into ascii table', level=8)
    gotdiv = .false.
    if (noline > 0) then
      if (associated(acc_nodivisions)) nullify(acc_nodivisions)
      acc_nodivisions => ListGetIntegerArray(acc_params, 'Polyline Divisions', gotdiv)
      if (.not. gotdiv) then
        acc_nodivisions => ListGetIntegerArray(acc_params, 'Polyline Intervals', gotdiv)
      end if
      if (gotdiv) then
        if (size(acc_nodivisions) < noline + count(acc_saveaxis)) then
          call Fatal(caller, 'Polyline divisions size too small!')
        end if
      end if
    end if
    acc_gotdiv = gotdiv
    if (noline > 0 .or. any(acc_saveaxis(1:acc_dim))) then
      acc_notests = 0
      if (.not. gotdiv) then
        if (acc_solver % timesvisited == 0) then
          call FindMeshEdges(acc_mesh, .false.)
        end if
        if (acc_dim == 2 .or. acc_intersectedge) then
          nofaces = acc_mesh % numberofedges
        else
          nofaces = acc_mesh % numberoffaces
        end if
        acc_nofaces = nofaces
      end if
      if (gotdiv) then
        t = maxval(acc_nodivisions) + 1
      else
        t = acc_mesh % numberofnodes
      end if
      if (allocated(acc_linetag)) deallocate(acc_linetag)
      allocate (acc_linetag(0:t))
      do line = 1, noline + acc_noaxis
        acc_linetag = .false.
        if (line <= noline) then
          acc_linenodes % x(1:2) = acc_pointcoordinates(2 * line - 1:2 * line, 1)
          acc_linenodes % y(1:2) = acc_pointcoordinates(2 * line - 1:2 * line, 2)
          if (acc_dim == 3) then
            acc_linenodes % z(1:2) = acc_pointcoordinates(2 * line - 1:2 * line, 3)
          else
            acc_linenodes % z(1:2) = 0.0_dp
          end if
        else
          if (.not. acc_saveaxis(line - noline)) cycle
          if (line - noline == 1) then
            acc_linenodes % x(1) = minval(acc_mesh % nodes % x)
            acc_linenodes % x(2) = maxval(acc_mesh % nodes % x)
            acc_linenodes % y(1:2) = 0.0_dp
            acc_linenodes % z(1:2) = 0.0_dp
          else if (line - noline == 2) then
            acc_linenodes % x(1:2) = 0.0_dp
            acc_linenodes % y(1) = minval(acc_mesh % nodes % y)
            acc_linenodes % y(2) = maxval(acc_mesh % nodes % y)
            acc_linenodes % z(1:2) = 0.0_dp
          else
            acc_linenodes % x(1:2) = 0.0_dp
            acc_linenodes % y(1:2) = 0.0_dp
            acc_linenodes % z(1) = minval(acc_mesh % nodes % z)
            acc_linenodes % z(2) = maxval(acc_mesh % nodes % z)
          end if
        end if
        acc_maxboundary = acc_maxboundary + 1
        if (gotdiv) then
          call sdl_poly_divided(line, noline)
        else
          call sdl_poly_faces(line, noline)
        end if
      end do
      if (acc_notests > 0) then
        call Info(caller, 'Number of candidate nodes: ' // trim(sdl_i2s(acc_notests)), level=8)
      end if
      call Info(caller, 'Number of nodes in specified lines: ' // &
          trim(sdl_i2s(acc_savenodes(2))))
      if (allocated(acc_linetag)) deallocate(acc_linetag)
    end if
  end subroutine sdl_save_polylines

  ! -- divided-polyline branch (verbatim) ---------------------------------------
  subroutine sdl_poly_divided(line, noline)
    integer :: line, noline
    real(kind=dp) :: r0(3), r1(3), dr(3), s1(3), ds(3), localcoord(3)
    real(kind=dp) :: mincoord(3), maxcoord(3), globalcoord(3)
    real(kind=dp) :: linen(3), linet1(3), linet2(3)
    real(kind=dp) :: s_len, linepos, tanprod(2), eps, detj
    logical :: stat_dummy
    integer :: nsize, t, n, m, i, j, imin, imax
    type(element_t), pointer :: el
    integer, pointer :: nodeindexes(:)
    real(kind=dp) :: dummy3(3) = 0.0_dp
    eps = 1.0e-5_dp
    r0 = [acc_linenodes % x(1), acc_linenodes % y(1), acc_linenodes % z(1)]
    r1 = [acc_linenodes % x(2), acc_linenodes % y(2), acc_linenodes % z(2)]
    dr = r1 - r0
    s_len = sqrt(sum(dr**2))
    linen = dr / s_len
    call TangentDirections(linen, linet1, linet2)
    nsize = acc_nodivisions(line)
    do t = 1, acc_mesh % numberofbulkelements + acc_mesh % numberofboundaryelements
      if (t <= acc_mesh % numberofbulkelements) then
        if (acc_intersectedge) cycle
      else
        if (.not. acc_intersectedge) cycle
      end if
      el => acc_mesh % elements(t)
      if (parenv % pes > 1) then
        if (el % partindex /= parenv % mype) cycle
      end if
      n = el % type % numberofnodes
      m = GetElementCorners(el)
      nodeindexes => el % nodeindexes
      if (acc_maskwithperm) then
        if (any(acc_solver % variable % perm(nodeindexes) == 0)) cycle
      end if
      acc_elementnodes % x(1:n) = acc_mesh % nodes % x(nodeindexes)
      acc_elementnodes % y(1:n) = acc_mesh % nodes % y(nodeindexes)
      acc_elementnodes % z(1:n) = acc_mesh % nodes % z(nodeindexes)
      do i = 1, m
        s1 = [acc_elementnodes % x(i), acc_elementnodes % y(i), &
            acc_elementnodes % z(i)]
        ds = (s1 - r0) / s_len
        localcoord(1) = sum(ds * linen)
        localcoord(2) = sum(ds * linet1)
        localcoord(3) = sum(ds * linet2)
        if (i == 1) then
          mincoord = localcoord
          maxcoord = localcoord
        else
          do j = 1, 3
            mincoord(j) = min(mincoord(j), localcoord(j))
            maxcoord(j) = max(maxcoord(j), localcoord(j))
          end do
        end if
      end do
      do i = 1, 2
        tanprod(i) = mincoord(i + 1) * maxcoord(i + 1)
      end do
      if (tanprod(1) > eps) cycle
      if (acc_dim == 3 .and. .not. acc_intersectedge) then
        if (tanprod(2) > eps) cycle
      end if
      imin = max(0, ceiling(nsize * mincoord(1)))
      imax = min(nsize, floor(nsize * maxcoord(1)))
      do i = imin, imax
        acc_notests = acc_notests + 1
        if (acc_linetag(i)) cycle
        globalcoord = r0 + i * dr / nsize
        if (PointInElement(el, acc_elementnodes, globalcoord, localcoord, &
            usolver=acc_solver, localeps=eps)) then
          if (.not. allocated(acc_basis)) then
            allocate (acc_basis(acc_mesh % maxelementnodes), &
                acc_dbasisdx(acc_mesh % maxelementnodes, 3))
          end if
          stat_dummy = ElementInfo(el, acc_elementnodes, localcoord(1), &
              localcoord(2), localcoord(3), detj, acc_basis, acc_dbasisdx)
          acc_localcoord = localcoord
          acc_linetag(i) = .true.
          acc_savenodes(2) = acc_savenodes(2) + 1
          linepos = 1.0_dp * i / nsize + 2 * (line - 1)
          n = el % type % numberofnodes
          call sdl_write_fields_impl(el % elementindex, line, i, 0, &
              .false., dummy3, .false., linepos, .true., .false.)
        end if
      end do
    end do
  end subroutine sdl_poly_divided


  ! -- faces-intersection polyline branch (verbatim) ---------------------------
  subroutine sdl_poly_faces(line, noline)
    integer :: line, noline
    real(kind=dp) :: deteps = 1.0e-10_dp
    logical :: inside
    integer :: t, n, i, k
    real(kind=dp) :: linepos
    type(element_t), pointer :: el
    integer, pointer :: nodeindexes(:)
    do t = 1, acc_nofaces
      if (acc_dim == 2 .or. acc_intersectedge) then
        el => acc_mesh % edges(t)
      else
        el => acc_mesh % faces(t)
      end if
      n = el % type % numberofnodes
      nodeindexes => el % nodeindexes
      acc_elementnodes % x(1:n) = acc_mesh % nodes % x(nodeindexes)
      acc_elementnodes % y(1:n) = acc_mesh % nodes % y(nodeindexes)
      if (acc_dim == 3) then
        acc_elementnodes % z(1:n) = acc_mesh % nodes % z(nodeindexes)
      else
        acc_elementnodes % z(1:n) = 0.0_dp
      end if
      if (.not. allocated(acc_basis)) then
        allocate (acc_basis(acc_mesh % maxelementnodes), &
            acc_dbasisdx(acc_mesh % maxelementnodes, 3))
      end if
      if (acc_intersectcoord /= 0) then
        call GlobalToLocalCoordsReduced(el, acc_elementnodes, n, acc_linenodes, &
            deteps, inside, acc_basis, i, linepos, acc_localcoord)
      else
        call GlobalToLocalCoords(el, acc_elementnodes, n, acc_linenodes, &
            deteps, inside, acc_basis, i, linepos, acc_localcoord)
      end if
      if (.not. inside) cycle
      if (1.0_dp - maxval(acc_basis(1:n)) < 1.0e-3_dp) then
        if (acc_linetag(nodeindexes(i))) cycle
        acc_linetag(nodeindexes(i)) = .true.
      end if
      acc_savenodes(2) = acc_savenodes(2) + 1
      linepos = linepos + 2 * (line - 1)
      call sdl_write_fields_impl(el % elementindex, acc_maxboundary, &
          nodeindexes(i), 0, .false., acc_dummy3, .false., &
          linepos, .true., .false.)
    end do
  end subroutine sdl_poly_faces

  ! -- SaveCircleLines (verbatim phase) -----------------------------------------
  subroutine sdl_save_circles() bind(c, name="sdl_save_circles")
    real(kind=dp) :: cylcoord(3), radius, phi, rtol, pi2
    real(kind=dp) :: r0(3), dr(3), linen(3), linet1(3), linet2(3)
    real(kind=dp) :: s1(3), localcoord(3), mincoord(3), maxcoord(3), globalcoord(3)
    real(kind=dp) :: detj, linepos
    integer :: line, noline, nodims, nsize, t, n, m, i, j, ii, imin, imax
    logical :: stat_dummy
    type(element_t), pointer :: el
    integer, pointer :: nodeindexes(:)
    type(solver_t), pointer :: psolver
    psolver => acc_solver
    pi2 = 2.0_dp * 3.14159265358979323846_dp
    if (associated(acc_pointcoordinates)) nullify(acc_pointcoordinates)
    acc_pointcoordinates => ListGetConstRealArray(acc_params, &
        'Circle Coordinates', acc_gotdiv)
    if (.not. acc_gotdiv) return
    call Info(caller, 'Saving circular lines into ascii table', level=8)
    noline = size(acc_pointcoordinates, 1)
    nodims = size(acc_pointcoordinates, 2)
    if (nodims /= 7) then
      call Fatal(caller, 'By construction the circle is defined by 7 values: ')
    end if
    if (associated(acc_nodivisions)) nullify(acc_nodivisions)
    acc_nodivisions => ListGetIntegerArray(acc_params, 'Circle Divisions', acc_gotdiv)
    if (.not. acc_gotdiv) then
      call Fatal(caller, '> Circle Coordinates < requires > Circle Divisions < !')
    end if
    if (size(acc_nodivisions) < noline) then
      call Fatal(caller, 'Polyline divisions size too small!')
    end if
    call Info(caller, 'Saving data on given circles', level=7)
    acc_notests = 0
    t = maxval(acc_nodivisions)
    if (allocated(acc_linetag)) deallocate(acc_linetag)
    allocate (acc_linetag(t))
    if (.not. allocated(acc_basis)) then
      allocate (acc_basis(acc_mesh % maxelementnodes), &
          acc_dbasisdx(acc_mesh % maxelementnodes, 3))
    end if
    do line = 1, noline
      acc_maxboundary = acc_maxboundary + 1
      acc_linetag = .false.
      nsize = acc_nodivisions(line)
      r0 = acc_pointcoordinates(line, 1:3)
      dr = acc_pointcoordinates(line, 4:6)
      radius = acc_pointcoordinates(line, 7)
      linen = dr / sqrt(sum(dr**2))
      call TangentDirections(linen, linet1, linet2)
      do t = 1, acc_mesh % numberofbulkelements
        el => acc_mesh % elements(t)
        if (parenv % pes > 1) then
          if (el % partindex /= parenv % mype) cycle
        end if
        n = el % type % numberofnodes
        m = GetElementCorners(el)
        nodeindexes => el % nodeindexes
        if (acc_maskwithperm) then
          if (any(acc_solver % variable % perm(nodeindexes) == 0)) cycle
        end if
        acc_elementnodes % x(1:n) = acc_mesh % nodes % x(nodeindexes)
        acc_elementnodes % y(1:n) = acc_mesh % nodes % y(nodeindexes)
        acc_elementnodes % z(1:n) = acc_mesh % nodes % z(nodeindexes)
        do i = 1, m
          s1(1) = acc_elementnodes % x(i)
          s1(2) = acc_elementnodes % y(i)
          s1(3) = acc_elementnodes % z(i)
          s1 = s1 - r0
          localcoord(1) = sum(s1 * linet1)
          localcoord(2) = sum(s1 * linet2)
          localcoord(3) = sum(s1 * linen)
          cylcoord(1) = sqrt(localcoord(1)**2 + localcoord(2)**2)
          cylcoord(2) = atan2(localcoord(2), localcoord(1)) / pi2
          if (cylcoord(2) < 0) cylcoord(2) = cylcoord(2) + 1.0_dp
          cylcoord(3) = localcoord(3)
          if (i == 1) then
            mincoord = cylcoord
            maxcoord = cylcoord
          else
            do j = 1, 3
              mincoord(j) = min(mincoord(j), cylcoord(j))
              maxcoord(j) = max(maxcoord(j), cylcoord(j))
            end do
          end if
        end do
        if (mincoord(3) * maxcoord(3) > 0.0_dp) cycle
        rtol = 0.5 * (maxcoord(1) - mincoord(1))
        if (mincoord(1) > radius + rtol) cycle
        if (maxcoord(1) < radius - rtol) cycle
        if (maxcoord(2) - mincoord(2) > 0.5) then
          phi = mincoord(2)
          mincoord(2) = maxcoord(2)
          maxcoord(2) = mincoord(2) + 1.0_dp
        end if
        imin = ceiling(nsize * mincoord(2) - 1.0d-20)
        imax = floor(nsize * maxcoord(2) + 1.0d-20)
        do i = imin, imax
          acc_notests = acc_notests + 1
          ii = i
          if (ii <= 0) ii = ii + nsize
          if (ii > nsize) ii = ii - nsize
          if (acc_linetag(ii)) cycle
          phi = (ii * pi2) / nsize
          globalcoord = r0 + radius * cos(phi) * linet1 + &
              radius * sin(phi) * linet2
          if (PointInElement(el, acc_elementnodes, globalcoord, localcoord)) then
            stat_dummy = ElementInfo(el, acc_elementnodes, localcoord(1), &
                localcoord(2), localcoord(3), detj, acc_basis, &
                usolver=psolver)
            acc_linetag(ii) = .true.
            acc_savenodes(3) = acc_savenodes(3) + 1
            linepos = 20 + 2 * line + 1.0_dp * ii / nsize
            acc_localcoord = localcoord
            call sdl_write_fields_impl(el % elementindex, line, ii, &
                0, .false., acc_dummy3, .false., linepos, .true., .false.)
          end if
        end do
      end do
    end do
    if (acc_notests > 0) then
      call Info(caller, 'Number of candidate nodes', level=8)
    end if
    call Info(caller, 'Number of nodes in specified circle', level=8)
    if (allocated(acc_linetag)) deallocate(acc_linetag)
  end subroutine sdl_save_circles


  ! -- SaveIsoCurves (verbatim phase, 2D only) -----------------------------------
  subroutine sdl_save_isocurves() bind(c, name="sdl_save_isocurves")
    real(kind=dp) :: f0, f1, f2, q
    integer :: line, t, n, i, i1, i2, k
    type(variable_t), pointer :: isosurfvar
    integer, pointer :: isosurfperm(:)
    real(kind=dp), pointer :: isosurf(:)
    character(len=64) :: name, varname
    type(element_t), pointer :: el
    integer, pointer :: nodeindexes(:)
    logical :: found, gotvar
    if (.not. ListGetLogical(acc_params, 'Save Isocurves', found)) return
    call Info(caller, 'Saving isocurves into ascii table', level=8)
    if (acc_dim == 3) then
      call Fatal(caller, 'Isocurves can only be saved in 2D')
    end if
    if (allocated(acc_linetag)) deallocate(acc_linetag)
    allocate (acc_linetag(0:acc_mesh % numberofnodes))
    if (acc_solver % timesvisited == 0) then
      call FindMeshEdges(acc_mesh, .false.)
    end if
    acc_nofaces = acc_mesh % numberofedges
    line = 0
    do
      line = line + 1
      acc_linetag = .false.
      write (name, '(A,I0)') 'IsoSurface Variable ', line
      varname = ListGetString(acc_params, name, gotvar)
      write (name, '(A,I0)') 'IsoSurface Value ', line
      f0 = ListGetCReal(acc_params, name, found)
      if (gotvar) then
        isosurfvar => VariableGet(acc_mesh % variables, varname)
        if (.not. associated(isosurfvar)) then
          call Warn(caller, 'Isosurface variable not given!')
          exit
        end if
        isosurfperm => isosurfvar % perm
        isosurf => isosurfvar % values
      else
        if (line == 1) then
          call Warn(caller, 'No Isosurface Variable 1 defined!')
          exit
        end if
        if (.not. found) exit
      end if
      f1 = minval(isosurf)
      f2 = maxval(isosurf)
      if (f0 <= f1 .or. f0 >= f2) then
        call Warn(caller, 'Isosurface value not within range!')
        cycle
      end if
      acc_maxboundary = acc_maxboundary + 1
      if (.not. allocated(acc_basis)) then
        allocate (acc_basis(acc_mesh % maxelementnodes), &
            acc_dbasisdx(acc_mesh % maxelementnodes, 3))
      end if
      do t = 1, acc_nofaces
        el => acc_mesh % edges(t)
        n = el % type % numberofnodes
        nodeindexes => el % nodeindexes
        acc_elementnodes % x(1:n) = acc_mesh % nodes % x(nodeindexes)
        acc_elementnodes % y(1:n) = acc_mesh % nodes % y(nodeindexes)
        acc_elementnodes % z(1:n) = 0.0_dp
        i1 = isosurfperm(nodeindexes(1))
        if (i1 == 0) cycle
        f1 = isosurf(i1) - f0
        i2 = isosurfperm(nodeindexes(2))
        if (i2 == 0) cycle
        f2 = isosurf(i2) - f0
        if (f1 * f2 >= 0.0_dp) cycle
        if (acc_maskwithperm) then
          if (any(acc_solver % variable % perm(nodeindexes(1:2)) == 0)) cycle
        end if
        q = abs(f2) / (abs(f1) + abs(f2))
        acc_basis(1:n) = 0.0_dp
        acc_basis(1) = q
        acc_basis(2) = 1 - q
        if (acc_basis(1) > acc_basis(2)) then
          k = nodeindexes(1)
        else
          k = nodeindexes(2)
          q = 1 - q
        end if
        if (q > 0.999) then
          if (acc_linetag(k)) cycle
          acc_linetag(k) = .true.
        end if
        acc_savenodes(4) = acc_savenodes(4) + 1
        call sdl_write_fields_impl(el % elementindex, acc_maxboundary, &
            k, 0, .false., acc_dummy3, .false., -1.0_dp, .true., .false.)
      end do
    end do
    if (allocated(acc_linetag)) deallocate(acc_linetag)
  end subroutine sdl_save_isocurves


  ! -- SaveVariableNames (verbatim: build + write .names file) -------------------
  subroutine sdl_save_names() bind(c, name="sdl_save_names")
    integer :: namesunit, ivar, i, j, istat, noeigenvalues
    logical :: dgvar, ipvar, elemvar, edgebasis
    type(variable_t), pointer :: var
    if (acc_solver % timesvisited == 0 .and. acc_noresults > 0 .and. &
        (.not. acc_parallel .or. parenv % mype == 0)) then
      if (allocated(acc_valuenames)) deallocate(acc_valuenames)
      allocate (acc_valuenames(acc_noresults + 5), stat=istat)
      if (istat /= 0) call Fatal(caller, 'Memory allocation error for ValueNames')
      acc_no = 0
      do ivar = -2, acc_novar
        var => sdl_getvar_n(ivar, acc_comps)
        if (associated(var % eigenvectors)) then
          noeigenvalues = size(var % eigenvalues)
          do j = 1, noeigenvalues
            do i = 1, var % dofs
              if (i == 1) then
                write (acc_valuenames(acc_no + (j - 1) * var % dofs + i), &
                    '(A,I0,A,A,I2,A,2ES20.11E3)') &
                    'Eigen ', j, ' ', trim(var % name), i, &
                    '   EigenValue = ', var % eigenvalues(j)
              else
                write (acc_valuenames(acc_no + (j - 1) * var % dofs + i), &
                    '(A,I0,A,A,I2)') &
                    'Eigen ', j, ' ', trim(var % name), i
              end if
            end do
          end do
          acc_no = acc_no + var % dofs * noeigenvalues
        else
          dgvar = (var % type == variable_on_nodes_on_elements)
          ipvar = (var % type == variable_on_gauss_points)
          elemvar = (var % type == variable_on_elements)
          edgebasis = (var % type == variable_on_edges)
          if (edgebasis) then
            acc_valuenames(acc_no + 1) = trim(var % name) // ' {e} 1'
            acc_valuenames(acc_no + 2) = trim(var % name) // ' {e} 2'
            acc_valuenames(acc_no + 3) = trim(var % name) // ' {e} 3'
            acc_no = acc_no + 3
            if (acc_avbasis) then
              acc_no = acc_no + 1
              acc_valuenames(acc_no) = trim(var % name) // ' nodal'
            end if
          else if (acc_comps > 1) then
            acc_no = acc_no + 1
            acc_valuenames(acc_no) = trim(var % name)
            if (acc_comps >= 2) then
              var => sdl_getvar_n(ivar, component=2)
              acc_no = acc_no + 1
              acc_valuenames(acc_no) = trim(var % name)
            end if
            if (acc_comps >= 3) then
              var => sdl_getvar_n(ivar, component=3)
              acc_no = acc_no + 1
              acc_valuenames(acc_no) = trim(var % name)
            end if
          else if (var % dofs == 1) then
            acc_no = acc_no + 1
            acc_valuenames(acc_no) = trim(var % name)
          else
            do i = 1, var % dofs
              acc_no = acc_no + 1
              acc_valuenames(acc_no) = ComponentName(var % name, i)
            end do
          end if
        end if
      end do
      if (acc_calculateflux) then
        acc_valuenames(acc_no + 1) = 'Flux 1'
        acc_valuenames(acc_no + 2) = 'Flux 2'
        acc_valuenames(acc_no + 3) = 'Flux 3'
        acc_no = acc_no + 3
      end if
      namesunit = 20
      open (namesunit, file=trim(acc_outputdir) // trim(acc_sidefile) // '.names', &
          status='unknown')
      do i = 1, acc_no
        write (namesunit, '(A)') trim(acc_valuenames(i))
      end do
      close (namesunit)
    end if
  end subroutine sdl_save_names

  ! -- SaveLine main driver (calls all phases in source order) -------------------
  subroutine sdl_saveline() bind(c, name="sdl_saveline")
    call sdl_save_existing()
    call sdl_save_polylines()
    call sdl_save_circles()
    call sdl_save_isocurves()
    call sdl_save_names()
  end subroutine sdl_saveline

end module sd_line_rest



