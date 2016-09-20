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

  ! position vector of null
  double precision :: r(3)
  double precision :: h, h0
  double precision :: slowdown

  integer :: iring, iline, inull, i, inullchk

  ! null parameters
  integer :: sign
  double precision, dimension(3) :: fan, spine
  double precision :: theta, phi
  double precision, allocatable, dimension(:) :: xs, ys, zs

  integer :: nlines
  integer :: nrings, npoints
  integer :: nperring(0:ringsmax)

  logical :: exitcondition, out

  print*,'#######################################################################'
  print*,'#                      Separatrix Surface Finder                      #'
  print*,'#######################################################################'

  call omp_set_num_threads(4) ! have it work on 4 threads (If machine has >4 cores this should be larger, if fewer than 4 coures, this should be smaller)

  filename = defaultfilename
  if (command_argument_count() > 0) then
    do i = 1, command_argument_count()
      call get_command_argument(i,arg)
      if (arg(1:5) == 'data=') then
        filename = trim(arg(6:))
      endif
    enddo
  endif
  
  !read in data
  open(unit=10,file=filename,access='stream')
    read(10), nx, ny, nz !number of vertices
    allocate(bgrid(nx,ny,nz,3))
    allocate(x(nx), y(ny), z(nz))
    read(10), bgrid(:,:,:,1)
    read(10), bgrid(:,:,:,2)
    read(10), bgrid(:,:,:,3)
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

  !read in null data

  open(unit=10,file='output/nulls.dat',form='unformatted')
    read(10) nnulls
    allocate(signs(nnulls),rnulls(3,nnulls),spines(3,nnulls),fans(3,nnulls))
    read(10) signs
    read(10) rnulls, spines, fans
  close(10)
  
  rnullsalt = rnulls
  
  if (coord_type == 2) then
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
  elseif (coord_type == 3) then
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
  endif
  
  do inull = 1, nnulls ! loop over all nulls
    print*, ''
    print*, 'Null number', inull, 'of', nnulls

    r = rnulls(:, inull)
    spine = spines(:, inull)
    fan = fans(:, inull)
    sign = signs(inull)

    ! get number of start points
    nlines = nstart

    theta = acos(fan(3))
    phi = atan(fan(2), fan(1))
    allocate(xs(nlines), ys(nlines), zs(nlines))
    call get_startpoints(theta,phi, xs,ys,zs)

    allocate(line1(3,nlines), line2(3,nlines))
    allocate(break(nlines))

    ! add start points to first ring relative to null
    do iline = 1, nlines !Go through each start point
      line1(1,iline) = r(1) + xs(iline)
      line1(2,iline) = r(2) + ys(iline)
      line1(3,iline) = r(3) + zs(iline)
    enddo
    line2 = line1

    write(fname,fmt) inull
    open(unit=12,file='output/separator'//trim(fname)//'.dat',form='unformatted',access='stream',status='replace')
    open(unit=20,file='output/everything'//trim(fname)//'.dat',access='stream',status='replace')

    break = 0
    nseps = 0
    nrings = ringsmax
    nperring = 0
    nperring(0) = nlines
    slowdown = 1d0

    write(20) nrings, 0, ringsmax+1
    write(20) nperring
    write(20) [(iline,iline=1,nlines)], break, line1

    exitcondition = .false.

    do iring = 1, ringsmax ! loop over number of rings we want
      if (sign .eq. 0) then ! skip null which is uncharacterised
        print*,'Null has zero sign'
        exit
      endif

      allocate(endpoints(nlines), association(nlines))
      endpoints = 0

      if (iring < 50) then
        h0 = 5d-2/slowdown
      else
        h0 = 25d-2/slowdown
      endif

      !$OMP PARALLEL DO private(r,h)
      do iline = 1, nlines ! loop over all points in ring (in parallel do)

        r(:) = line1(:,iline)
        h = h0

        call trace_line(r,sign,h) !trace line by a distance of h

        line2(:,iline) = line2(:,iline) + r(:) - line1(:,iline)
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
      !$OMP END PARALLEL DO

      if (nlines > pointsmax) then
      ! exit if too many points on ring 
        print*, 'Too many points on ring', nlines, iring
        exitcondition = .true.
      endif

      if (sum(endpoints)/nlines == 1) then
        ! exit if all points have reached outer boundary (left box)
        print*, 'All fan points have reached the outer boundary', iring
        exitcondition = .true.
      endif
      
      if (ierror == 1) then
        print*, 'Tracing has failed', iring
        exitcondition = .true.
      endif
      
      ! print*,'Checking at null', iring, nlines, inull
      if (iring < 50) then
        maxdist = 0.1d0*h0*slowdown!0.0075d0
      elseif (iring < 100) then
        maxdist = 0.08d0*h0*slowdown!0.0075d0=0.15*0.05
      else
        maxdist = 0.6d0*h0*slowdown!0.15
      endif
      nulldist = 1.4d0*h0*slowdown !0.6
      mindist = maxdist/3
      ! print*, iring, h0, nulldist, maxdist, mindist, nlines

      ! remove points from ring if necessary
      call remove_points(nlines,iring)

      allocate(nearnull(nlines))
      nearnull = 0
      slowdown = 1d0
      do iline = 1, nlines
        do inullchk = 1, nnulls
          if (signs(inullchk) == signs(inull)) cycle
          if (dist(rnulls(:,inullchk), line1(:, iline)) < 2.5*nulldist) then
            !print*, 'close to null'
            nearnull(iline) = 1
            exit
          endif
        enddo
      enddo
      if (sum(nearnull) > 0) then
        ! print*, 'slowing down, adding points near null'
        slowdown = 2d0
      endif
      !add points to ring if necessary
      call add_points(nlines,iring)
      deallocate(nearnull)

      !determine if any lines are separators
      call sep_detect(nlines,inull,iring)

      nperring(iring) = nlines
      !Write ring and data to file separator????.dat
      write(20) association, break, line1

      deallocate(endpoints, association)

      if (exitcondition) then
        nrings = iring
        exit
      endif

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

  enddo

  call system_clock(tstop, count_rate)
  print*, 'TIME = ', dble(tstop - tstart)/dble(count_rate), "(", dble(tstop - tstart)/dble(count_rate)/60, "minutes)"
  print*, 'Done!' 

end program
