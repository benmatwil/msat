!code to find the separatrix surface and separators of null points
program ssfind

  use omp_lib
  use params
  use common
  use trace
  use ring

  implicit none

  integer(int64) :: tstart, tstop, count_rate !to time program
  integer(int32) :: nx, ny, nz !size of grid

  real(np) :: r(3), rmove(3)
  real(np) :: h, h0
  real(np) :: slowdown !, shift(3)
  integer(int32) :: nearflag

  integer(int32) :: iring, iline, inull, jnull

  ! null parameters
  real(np) :: theta, phi, factor
  real(np), allocatable, dimension(:) :: xs, ys, zs

  ! spine tracer
  real(np), dimension(3) :: rspine
  real(np), dimension(:,:), allocatable :: spine
  real(np) :: hspine
  integer(int32) :: dir

  integer(int32) :: nlines
  integer(int32) :: nrings
  integer(int32) :: nperring(0:ringsmax)

  logical :: exitcondition, out

  ! file writing
  integer(int64) :: ia, ib, ip
  integer(int32) :: isep, nullnum1, nullnum2, linenum, ringnum
  real(np), dimension(:,:), allocatable :: rsep, rring, linewrite
  integer(int32), dimension(:), allocatable :: brk

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
    allocate(rnulls(3,nnulls), rnullsreal(3,nnulls))
    read(10) rnulls
    read(10) rnullsreal
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

  ! stores all info regarding size of rings
  open(unit=10,file=trim(fileout)//'-ringinfo.dat',status='replace',access='stream')
  ! stores all the coordinates of rings in original coordinate system
  open(unit=20,file=trim(fileout)//'-rings.dat',status='replace',access='stream')
  ! stores all the connection info about each separator
  open(unit=40,file=trim(fileout)//'-connectivity.dat',status='replace',access='stream')
  ! stores all the coordinates of the separator lines in original coordinate system
  open(unit=50,file=trim(fileout)//'-separators.dat',status='replace',access='stream')
  ! stores all the coordinates of the spine lines in original coordinate system
  open(unit=60,file=trim(fileout)//'-spines.dat',status='replace',access='stream')

  write(10) ringsmax+1

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

    out = .false.

    ! add start points to first ring relative to null
    do iline = 1, nlines !Go through each start point
      line1(1,iline) = rnulls(1,inull) + xs(iline)
      line1(2,iline) = rnulls(2,inull) + ys(iline)
      line1(3,iline) = rnulls(3,inull) + zs(iline)
      call edgecheck(line1(:,iline))
      if (outedge(line1(:,iline))) then
        factor = 0.9_np
        do while (outedge(line1(:,iline)))
          line1(1,iline) = rnulls(1,inull) + xs(iline)*factor
          line1(2,iline) = rnulls(2,inull) + ys(iline)*factor
          line1(3,iline) = rnulls(3,inull) + zs(iline)*factor
          call edgecheck(line1(:,iline))
          factor = factor*0.9_np
        enddo
      endif
    enddo
    line2 = line1

    open(unit=90, access='stream', status='scratch')
    open(unit=95, access='stream', status='scratch')

    break = 0
    nseps = 0
    nperring = 0
    nperring(0) = nlines
    slowdown = 1.0_np
    terror = 0

    linewrite = gtr(line1, x, y, z)
    write(90) [(iline,iline=1,nlines)], break, linewrite
    write(20) break, linewrite
    deallocate(linewrite)

    exitcondition = .false.
    !$OMP END SINGLE

    do iring = 1, ringsmax ! loop over number of rings we want

      if (signs(inull) .eq. 0) then ! skip null which is uncharacterised
        print*,'Null has zero sign'
        exit
      endif
      
      !$OMP SINGLE
#if debug
      print*, iring, nlines
#endif

      allocate(endpoints(nlines), association(nlines))

      if (iring < 50) then
        h0 = 4e-2_np/slowdown
      else
        h0 = 2e-1_np/slowdown
      endif
      !$OMP END SINGLE

      !$OMP WORKSHARE
      endpoints = 0
      !$OMP END WORKSHARE

      !$OMP DO private(r, h, out)
      do iline = 1, nlines ! loop over all points in ring (in parallel do)

        r = line1(:,iline)
        h = h0

        call trace_line(r, signs(inull), h) ! trace line by a distance of h

        line2(:,iline) = line2(:,iline) + (r - line1(:,iline))
        line1(:,iline) = r

        ! counter to see how many points on ring have reached outer boundary
        if (outedge(r)) then
          endpoints(iline) = 1
        else
          endpoints(iline) = 0
        endif

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
      call remove_points(nlines)

      !$OMP SINGLE
      allocate(nearnull(nlines))
      slowdown = 1.0_np
      !$OMP END SINGLE

      !$OMP WORKSHARE
      nearnull = 0
      !$OMP END WORKSHARE

      !$OMP DO private(jnull)
      do iline = 1, nlines
        do jnull = 1, nnulls
          if (signs(jnull) == signs(inull)) cycle
          if (dist(rnulls(:,jnull), line1(:, iline)) < 2.5*nulldist .or. &
            dist(rnullsalt(:,jnull), line1(:, iline)) < 2.5*nulldist) then
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
      do iline = 1, nlines
        call edgecheck(line1(:,iline))
      enddo
      !$OMP END SINGLE

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

      if (exitcondition) deallocate(endpoints, association)
      !$OMP END SINGLE

      if (exitcondition) then
        nrings = iring
        exit
      endif

      !$OMP SINGLE
      nperring(iring) = nlines
      ! Write ring and data to file separator????.dat
      linewrite = gtr(line1, x, y, z)
      write(90) association, break, linewrite
      write(20) break, linewrite

      deallocate(endpoints, association, linewrite)
      !$OMP END SINGLE

    enddo

    !$OMP SINGLE
    write(10) nperring

    ! trace separators...
    write(95, pos=1)
    if (nseps > 0) then
      print*, "Tracing separators"
      do isep = 1, nseps
        read(95) nullnum1, nullnum2, ringnum, linenum
        allocate(rsep(3, ringnum+3))
        do iring = ringnum, 0, -1
          call file_position(iring, nperring, linenum, ia, ib, ip)
          read(90, pos=ia) linenum
          read(90, pos=ip) rsep(:,iring+2)
        enddo
        rsep(:,1) = rnullsreal(:,nullnum1)
        rsep(:,ringnum+3) = rnullsreal(:,nullnum2)
        write(50) ringnum+3
        write(50) rsep
        deallocate(rsep)
      enddo
    endif

    ! don't need scratch ring file anymore
    close(95)

    ! ! write a subset of rings to 
    ! print*, 'Writing rings'
    ! write(30) nrings/skiprings + 1
    ! do iring = 0, nrings-1, skiprings
    !   call file_position(iring, nperring, 1, ia, ib, ip)
    !   allocate(rring(3,nperring(iring)), brk(nperring(iring)))
    !   read(90, pos=ib) brk
    !   read(90, pos=ip) rring
    !   write(30) nperring(iring), rring, brk
    !   deallocate(rring, brk)
    ! enddo
    close(90)

    ! Trace spines...
    print*, "Tracing spines"
    do dir = -1, 1, 2
      out = .false.
      allocate(spine(3,1))
      spine(:,1) = rnulls(:, inull)
      rspine = rnulls(:, inull) + dir*spines(:, inull)*1e-3_np ! pick a good factor
      do while (.not. outedge(rspine))
        call add_vector(spine, rspine)
        hspine = 5e-2_np
        call trace_line(rspine, -signs(inull), hspine)
        call edgecheck(rspine)
      enddo
      write(60) size(spine,2), gtr(spine, x, y, z)
      deallocate(spine)
    enddo

    if (nrings == ringsmax) print*, "Reached maximum number of rings"

    print*, 'number of separators=', nseps, 'number of rings', nrings

    deallocate(xs, ys, zs)
    deallocate(line1, line2, break)
    !$OMP END SINGLE

  enddo

  !$OMP END PARALLEL

  close(10)
  close(20)
  write(40) -1
  close(40)
  close(50)
  close(60)

  call system_clock(tstop, count_rate)
  print*, 'TIME = ', dble(tstop - tstart)/dble(count_rate), "(", dble(tstop - tstart)/dble(count_rate)/60, "minutes)"
  print*, 'Done!'

end program

!separate outedge and edgecheck