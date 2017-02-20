!code to find the separatrix surface and separators of null points
program ssfind

  use omp_lib
  use params
  use common
  use trace
  use ring

  implicit none

  character (len=8), parameter :: fmt='(I4.4)'
  character (len=5) :: fname

  integer*8 :: tstart, tstop, count_rate !to time program
  integer :: nx, ny, nz !size of grid

  real(np) :: r(3)
  real(np) :: h, h0
  real(np) :: slowdown, shift(3)
  integer :: nearflag

  integer :: iring, iline, inull, inullchk

  ! null parameters
  real(np) :: theta, phi
  real(np), allocatable, dimension(:) :: xs, ys, zs

  ! spine tracer
  real(np), dimension(3) :: rspine
  real(np), dimension(:,:), allocatable :: spine
  real(np) :: hspine
  integer :: dir

  integer :: nlines
  integer :: nrings
  integer :: nperring(0:ringsmax)

  logical :: exitcondition, out

  print*,'#######################################################################'
  print*,'#                      Separatrix Surface Finder                      #'
  print*,'#######################################################################'

  ! call omp_set_num_threads(nproc) ! have it work on 4 threads (If machine has >4 cores this should be larger, if fewer than 4 coures, this should be smaller)

  call filenames

  ! read in data
  open(unit=10, file=filein, access='stream', status='old')
    read(10) nx, ny, nz ! number of vertices
    allocate(bgrid(nx,ny,nz,3))
    allocate(x(nx), y(ny), z(nz))
    read(10) bgrid(:,:,:,1)
    read(10) bgrid(:,:,:,2)
    read(10) bgrid(:,:,:,3)
    read(10) x, y, z
  close(10)

  xmin = 1
  xmax = nx
  ymin = 1
  ymax = ny
  zmin = 1
  zmax = nz
  print*, "Data read in, grid dimensions are", nx, ny, nz

  call system_clock(tstart,count_rate) ! to time how long it takes

  ! read in null data

  open(unit=10, file=trim(fileout)//'-nullpos.dat', access='stream', status='old')
    read(10) nnulls
    allocate(rnulls(3,nnulls))
    read(10) rnulls
  close(10)
  open(unit=10, file=trim(fileout)//'-nulldata.dat', access='stream', status='old')
    read(10) nnulls
    allocate(signs(nnulls),spines(3,nnulls),fans(3,nnulls))
    read(10) signs, spines, fans
  close(10)

  rnullsalt = rnulls

#if spherical
  ! for sphericals
  do inull = 1, nnulls
    ! check whether null is at the lower phi boundary
    if (rnullsalt(3,inull) < zmin + 1) then
      rnullsalt(3,inull) = rnullsalt(3,inull) - 1
      call edgecheck(rnullsalt(:,inull))
      rnullsalt(3,inull) = rnullsalt(3,inull) + 1
    endif

    ! check whether null is at the upper phi boundary
    if (rnullsalt(3,inull) > zmax - 1) then
      rnullsalt(3,inull) = rnullsalt(3,inull) + 1
      call edgecheck(rnullsalt(:,inull))
      rnullsalt(3,inull) = rnullsalt(3,inull) - 1
    endif

    ! check whether null is at the lower theta boundary
    if (rnullsalt(2,inull) < ymin + 1) then
      rnullsalt(2,inull) = rnullsalt(2,inull) - 1
      call edgecheck(rnullsalt(:,inull))
      rnullsalt(2,inull) = rnullsalt(2,inull) - 1
    endif

    ! check whether null is at the upper theta boundary
    if (rnullsalt(2,inull) > ymax - 1) then
      rnullsalt(2,inull) = rnullsalt(2,inull) + 1
      call edgecheck(rnullsalt(:,inull))
      rnullsalt(2,inull) = rnullsalt(2,inull) + 1
    endif
  enddo
#elif cylindrical
  ! for cylindricals
  do inull = 1, nnulls
    ! check whether null is at lower phi boundary
    if (rnullsalt(2,inull) < ymin + 1) then
      rnullsalt(2,inull) = rnullsalt(2,inull) - 1
      call edgecheck(rnullsalt(:,inull))
      rnullsalt(2,inull) = rnullsalt(2,inull) + 1
    endif

    ! check whether null is at upper phi boundary
    if (rnullsalt(2,inull) > ymax - 1) then
      rnullsalt(2,inull) = rnullsalt(2,inull) + 1
      call edgecheck(rnullsalt(:,inull))
      rnullsalt(2,inull) = rnullsalt(2,inull) - 1
    endif
  enddo
#endif


  open(unit=30,file=trim(fileout)//'-spinelines.dat',status='replace',access='stream')

  !$OMP PARALLEL private(iring, iline, inull)
  do inull = 1, nnulls ! loop over all nulls
    !$OMP SINGLE
    print*, ''
    print*, 'Null number', inull, 'of', nnulls

    ! get number of start points
    nlines = nstart

    theta = acos(fans(3,inull))
    phi = atan(fans(2,inull), fans(1,inull))
    allocate(xs(nlines), ys(nlines), zs(nlines))
    call get_startpoints(theta,phi, xs,ys,zs)

    allocate(line1(3,nlines), line2(3,nlines))
    allocate(break(nlines))

    ! add start points to first ring relative to null
    do iline = 1, nlines !Go through each start point
      line1(1,iline) = rnulls(1,inull) + xs(iline)
      line1(2,iline) = rnulls(2,inull) + ys(iline)
      line1(3,iline) = rnulls(3,inull) + zs(iline)
    enddo
    line2 = line1

    write(fname,fmt) inull
    open(unit=12, file=trim(fileout)//'-separator'//trim(fname)//'.dat', access='stream', status='replace')
    open(unit=20, file=trim(fileout)//'-everything'//trim(fname)//'.dat', access='stream', status='replace')

    break = 0
    nseps = 0
    nrings = ringsmax
    nperring = 0
    nperring(0) = nlines
    slowdown = 1.0_np
    terror = 0
    out = .false.

    write(20) nrings, 0, ringsmax+1
    write(20) nperring
    write(20) [(iline,iline=1,nlines)], break, gtr(line1, x, y, z)

    exitcondition = .false.
    !$OMP END SINGLE

    do iring = 1, ringsmax ! loop over number of rings we want
#if debug
      !$OMP SINGLE
      print*, iring, nlines
      !$OMP END SINGLE
#endif

      if (signs(inull) .eq. 0) then ! skip null which is uncharacterised
        print*,'Null has zero sign'
        exit
      endif

      !$OMP SINGLE
      allocate(endpoints(nlines), association(nlines))
      !$OMP END SINGLE

      !$OMP WORKSHARE
      endpoints = 0
      !$OMP END WORKSHARE

      !$OMP SINGLE
      if (iring < 50) then
        h0 = 4d-2/slowdown
      else
        h0 = 20d-2/slowdown
      endif
      !$OMP END SINGLE

      !$OMP DO private(r, h, out, shift)
      do iline = 1, nlines ! loop over all points in ring (in parallel do)

        r(:) = line1(:,iline)
        h = h0

        call trace_line(r,signs(inull),h) !trace line by a distance of h

        shift = [0,0,0]
#if spherical
        if (abs(r(3) - line1(3,iline)) > 0.9_np*(zmax-zmin)) then
          if (r(3) - line1(3,iline) > 0) then
            shift(3) = zmax - zmin
          else
            shift(3) = -(zmax - zmin)
          endif
        elseif (abs(r(3) - line1(3,iline)) < 0.6_np*(zmax-zmin) .and. abs(r(3) - line1(3,iline)) > 0.4_np*(zmax-zmin)) then
          ! switch in theta untested
          print*, 'shifting - untested'
          if (r(3) - line1(3,iline) > 0) then
            shift(3) = (zmax - zmin)/2
          else
            shift(3) = -(zmax - zmin)/2
          endif
        endif
#elif cylindrical
        if (abs(r(2) - line1(2,iline)) > 0.9_num*(ymax-ymin)) then
          if (r(2) - line1(2,iline) > 0) then
            shift(2) = ymax - ymin
          else
            shift(2) = -(ymax - ymin)
          endif
        endif
#endif

        line2(:,iline) = line2(:,iline) + r(:) - line1(:,iline) - shift
        call edgecheck(r, out)
        ! counter to see how many points on ring have reached outer boundary
        if (out) then
          endpoints(iline) = 1
        else
          endpoints(iline) = 0
        endif
        line1(:,iline) = r(:)

        association(iline) = iline

      enddo
      !$OMP END DO

      !$OMP SINGLE
      ! print*,'Checking at null', iring, nlines, inull
      if (iring < 50) then
        maxdist = 0.1_np*h0*slowdown ! 0.0075
      elseif (iring < 100) then
        maxdist = 0.08_np*h0*slowdown ! 0.0075=0.15*0.05
      else
        maxdist = 0.6_np*h0*slowdown ! 0.15
      endif
      nulldist = 1.4_np*h0*slowdown ! 0.6
      mindist = maxdist/3
      ! print*, iring, h0, nulldist, maxdist, mindist, nlines
      !$OMP END SINGLE

      ! remove points from ring if necessary
      call remove_points(nlines, iring)

      !$OMP SINGLE
      allocate(nearnull(nlines))
      slowdown = 1.0_np
      !$OMP END SINGLE

      !$OMP WORKSHARE
      nearnull = 0
      !$OMP END WORKSHARE

      !$OMP DO private(inullchk)
      do iline = 1, nlines
        do inullchk = 1, nnulls
          if (signs(inullchk) == signs(inull)) cycle
          if (dist(rnulls(:,inullchk), line1(:, iline)) < 2.5*nulldist .or. &
            dist(rnullsalt(:,inullchk), line1(:, iline)) < 2.5*nulldist) then
            nearnull(iline) = 1
            exit
          endif
        enddo
      enddo
      !$OMP END DO

      !$OMP SINGLE
      nearflag = sum(nearnull)
      if (nearflag > 0) then
        ! print*, 'slowing down, adding points near null'
        slowdown = 2.0_np
      endif
      !$OMP END SINGLE

      ! add points to ring if necessary
      call add_points(nlines)

      !$OMP SINGLE
      deallocate(nearnull)
      !$OMP END SINGLE

      ! determine if any lines are separators
      if (nearflag > 0) call sep_detect(nlines,inull,iring)

      !$OMP SINGLE
      if (nlines > pointsmax) then
      ! exit if too many points on ring
        print*, 'Too many points on ring', nlines, iring
        exitcondition = .true.
      endif

      if (nlines == 0) then
        ! exit if all points have reached outer boundary (left box)
        print*, 'All fan points have reached the outer boundary', iring
        exitcondition = .true.
      endif

      if (terror == 1) then
        print*, 'Tracing has failed', iring
        exitcondition = .true.
      endif
      !$OMP END SINGLE

      !$OMP SINGLE
      nperring(iring) = nlines
      ! Write ring and data to file separator????.dat
      write(20) association, break, gtr(line1, x, y, z)

      deallocate(endpoints, association)
      !$OMP END SINGLE

      if (exitcondition) then
        nrings = iring
        exit
      endif

    enddo

    ! Trace spines...
    !$OMP SINGLE
    do dir = -1, 1, 2
      out = .false.
      allocate(spine(3,2))
      spine(:,1) = rnulls(:, inull)
      spine(:,2) = rnulls(:, inull) + dir*spines(:, inull)*1e-3_np ! pick a good factor
      rspine = spine(:,2)
      do while (.not. out)
        hspine = 5e-2_np
        call trace_line(rspine, -signs(inull), hspine)
        call edgecheck(rspine, out)
        call add_vector(spine, rspine) 
      enddo
      write(30) size(spine,2), gtr(spine, x, y, z)
      deallocate(spine)
    enddo

    if (nrings == ringsmax) print*, "Reached maximum number of rings"

    print*, 'number of separators=', nseps, 'number of rings', nrings

    write(12) -1
    close(12)

    write(20, pos=1) nrings, sum(nperring)
    write(20, pos=13) nperring
    close(20)

    deallocate(xs, ys, zs)
    deallocate(line1, line2, break)
    !$OMP END SINGLE

  enddo

  close(30)

  !$OMP END PARALLEL

  call system_clock(tstop, count_rate)
  print*, 'TIME = ', dble(tstop - tstart)/dble(count_rate), "(", dble(tstop - tstart)/dble(count_rate)/60, "minutes)"
  print*, 'Done!'

end program
