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

  double precision :: r(3)
  double precision :: h, h0
  double precision :: slowdown

  integer :: iring, iline, inull, inullchk

  ! null parameters
  integer :: sign, nearflag
  double precision, dimension(3) :: fan, spine
  double precision :: theta, phi
  double precision, allocatable, dimension(:) :: xs, ys, zs

  integer :: nlines
  integer :: nrings
  integer :: nperring(0:ringsmax)

  logical :: exitcondition, out

  print*,'#######################################################################'
  print*,'#                      Separatrix Surface Finder                      #'
  print*,'#######################################################################'

  call omp_set_num_threads(nproc) ! have it work on 4 threads (If machine has >4 cores this should be larger, if fewer than 4 coures, this should be smaller)

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

  open(unit=10, file='output/'//trim(fileout)//'-nulldata.dat', access='stream', status='old')
    read(10) nnulls
    allocate(signs(nnulls),rnulls(3,nnulls),spines(3,nnulls),fans(3,nnulls))
    read(10) rnulls
    read(10) signs, spines, fans
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
  
  !$OMP PARALLEL private(iring, iline, inull)

  do inull = 1, nnulls ! loop over all nulls
    !$OMP SINGLE
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
    open(unit=12, file='output/'//trim(fileout)//'-separator'//trim(fname)//'.dat' ,access='stream', status='replace')
    open(unit=20, file='output/'//trim(fileout)//'-everything'//trim(fname)//'.dat', access='stream', status='replace')

    break = 0
    nseps = 0
    nrings = ringsmax
    nperring = 0
    nperring(0) = nlines
    slowdown = 1d0
    terror = 0

    write(20) nrings, 0, ringsmax+1
    write(20) nperring
    write(20) [(iline,iline=1,nlines)], break, line1

    exitcondition = .false.
    !$OMP END SINGLE

    do iring = 1, ringsmax ! loop over number of rings we want
      if (sign .eq. 0) then ! skip null which is uncharacterised
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
        h0 = 5d-2/slowdown
      else
        h0 = 25d-2/slowdown
      endif
      !$OMP END SINGLE

      !$OMP DO private(r, h, out)
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
      !$OMP END DO

      !$OMP SINGLE
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
      
      if (terror == 1) then
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
      !$OMP END SINGLE
      
      ! remove points from ring if necessary
      call remove_points(nlines, iring)

      !$OMP SINGLE
      allocate(nearnull(nlines))
      slowdown = 1d0
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
        slowdown = 2d0
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
      nperring(iring) = nlines
      ! Write ring and data to file separator????.dat
      write(20) association, break, line1

      deallocate(endpoints, association)
      !$OMP END SINGLE

      if (exitcondition) then
        nrings = iring
        exit
      endif

    enddo

    !$OMP SINGLE
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

  !$OMP END PARALLEL

  call system_clock(tstop, count_rate)
  print*, 'TIME = ', dble(tstop - tstart)/dble(count_rate), "(", dble(tstop - tstart)/dble(count_rate)/60, "minutes)"
  print*, 'Done!' 

end program
