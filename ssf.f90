!code to find the separatrix surface and separators of null points
program ssfind

  use omp_lib
  use params
  use common
  use trace
  use ring

  implicit none

  character (len=8), parameter :: fmt='(I3.3)'
  character (len=5) :: fname

  integer :: nringss
  integer*8 :: tstart,tstop,count_rate !to time program
  integer :: nx, ny, nz !size of grid

  !position vector of null (and saved backup)
  double precision :: r(3)
  double precision :: a(3)
  double precision :: h

  integer :: i, j

  !null parameters
  integer :: sign
  double precision, dimension(3) :: fan, spine
  double precision :: theta, phi

  double precision, allocatable, dimension(:) :: nsepss
  double precision, allocatable, dimension(:) :: xs, ys, zs

  !number of lines
  integer :: nlines
  integer :: nnull

  integer :: nrings, npoints
  integer :: nperring(ringsmax)
  double precision :: circumference(ringsmax)

  logical :: exitcondition, out

  CALL OMP_SET_NUM_THREADS(8) !have it work on 4 threads (If machine has >4 cores this should be larger, if fewer than 4 coures, this should be smaller)
  
  !read in data
  open (unit=10,file=filename,access='stream')

    read(10), nx, ny, nz !number of vertices
    print*, nx, ny, nz
    
    allocate(bgrid(nx,ny,nz,3))
    allocate(x(nx), y(ny), z(nz))
    read(10), bgrid(:,:,:,1)
    read(10), bgrid(:,:,:,2)
    read(10), bgrid(:,:,:,3)
    read(10) x, y, z

    xmin = 1
    xmax = nx
    ymin = 1
    ymax = ny
    zmin = 1
    zmax = nz

  close(10)

  print*, 'Data read in'

  call SYSTEM_CLOCK(tstart,count_rate) !to time how long it takes

  !read in null data

  open (unit=10,file='output/nulls.dat',form='unformatted')
    read(10) nnulls
    allocate(signs(nnulls),rnulls(3,nnulls),spines(3,nnulls),fans(3,nnulls))
    read(10) signs
    read(10) rnulls, spines, fans
  close(10)
  
  rnullsalt = rnulls
  !for sphericals
  if (coord_type == 2) then 
    do i = 1, nnulls
      !check whether null is at the lower phi boundary
      if (rnullsalt(3,i) < zmin + 1) then
        rnullsalt(3,i) = rnullsalt(3,i) - 1
        call edgecheck(rnullsalt(:,i))
        rnullsalt(3,i) = rnullsalt(3,i) + 1
      endif
      
      !check whether null is at the upper phi boundary
      if (rnullsalt(3,i) > zmax - 1) then
        rnullsalt(3,i) = rnullsalt(3,i) + 1
        call edgecheck(rnullsalt(:,i))
        rnullsalt(3,i) = rnullsalt(3,i) - 1
      endif
      
      !check whether null is at the lower theta boundary
      if (rnullsalt(2,i) < ymin + 1) then
        rnullsalt(2,i) = rnullsalt(2,i) - 1
        call edgecheck(rnullsalt(:,i))
        rnullsalt(2,i) = rnullsalt(2,i) - 1
      endif
      
      !check whether null is at the upper theta boundary
      if (rnullsalt(2,i) > ymax - 1) then
        rnullsalt(2,i) = rnullsalt(2,i) + 1
        call edgecheck(rnullsalt(:,i))
        rnullsalt(2,i) = rnullsalt(2,i) + 1
      endif
    enddo
  endif
  
  allocate(nsepss(nnulls))

  !signs=-1*signs
  do nnull = 1, nnulls !loop over all nulls
    print*, ''
    print*, 'Null number', nnull, 'of', nnulls

    r = rnulls(:,nnull)
    spine = spines(:,nnull)
    fan = fans(:,nnull)
    sign = signs(nnull)

    a = r !backup null location

    !get start points
    nlines = nstart

    theta = acos(fan(3))
    phi = atan(fan(2),fan(1))

    allocate(xs(nlines),ys(nlines),zs(nlines))

    call get_startpoints(theta,phi,xs,ys,zs)

    allocate(line1(3,nlines), line2(3,nlines), add1(3,nlines), add2(3,nlines))
    allocate(association(nlines), break(nlines), remove(nlines), endpoints(nlines))

    !add start points to first ring relative to null
    do j = 1, nlines !Go through each start point
      line1(1,j) = r(1) + xs(j)
      line1(2,j) = r(2) + ys(j)
      line1(3,j) = r(3) + zs(j)
      association(j) = j
    enddo
    
    line2 = line1

    break = 0

    write(fname,fmt) nnull

    open(unit=12,file='output/separator'//trim(fname)//'.dat',form='unformatted')

    open(unit=20,file='output/everything'//trim(fname)//'.dat',access='stream',status='replace')

    nseps = 0
    nrings = 0
    nperring = 0
    nringss = 0
    circumference = 0

    write(20) nrings, npoints, ringsmax
    write(20) nperring

    exitcondition = .false.

    do i = 1, ringsmax !loop over number of rings we want
      if (sign .eq. 0) then !skip null which is uncharacterised
        print*,'Null has zero sign'
        exit
      endif

      nringss = nringss+1

      endpoints = 0
      add1 = 0
      add2 = 0
      remove = 0

      write(20) association
      nperring(i) = nlines

      ierror = 0
      !$OMP PARALLEL DO private(r,h)
        do j = 1, nlines !loop over all points in ring (in parallel do)

          r(:) = line1(:,j)

          if (i < 50) then
            h = 5d-2
          else
            h = 25d-2
          endif

          call trace_line(r,sign,h) !trace line by a distance of h

          line2(:,j) = line2(:,j) + r(:) - line1(:,j)
          call edgecheck(r, out)
          if (out) then !counter to see how many points on ring have reached outer boundary
            endpoints(j) = 1
          else
            endpoints(j) = 0
          endif
          line1(:,j) = r(:)
          
          association(j) = j

        enddo
      !$OMP END PARALLEL DO
      write(20), line1

      if (nlines > pointsmax) then
        print*, 'Too many points on ring', nlines, i
        exitcondition = .true. !exit if too many points on ring
      endif

      if (sum(endpoints)/nlines == 1) then
        print*, 'All fan points have reached the outer boundary', i
        exitcondition = .true. !exit if all points have reached outer boundary (left box)
      endif
      
      circumference(i) = circumference(i) + dist(line2(:,1),line2(:,nlines))
      do j = 2, nlines
        circumference(i) = circumference(i) + dist(line2(:,j),line2(:,j-1))
      enddo
      
      if (i > 1) then
        if (abs(circumference(i)-circumference(i-1)) .lt. 0.1*stepmin) then
          print*, 'Fan has stopped growing/shrinking', i
          exitcondition = .true. !exit if fan not changing size
        endif
      endif

      if (ierror == 1) then
        print*, 'Tracing has failed',i
        exitcondition = .true.
      endif
      
      !print*,'Checking at null', i, nlines, nnull
      call at_null(nlines,nnull,i) !determine if point is at null
      call remove_points(nlines,i) !remove points from ring if necessary
      call add_points(nlines,i) !add points to ring if necessary

      if (exitcondition) exit

    enddo
    
    print*, 'number of separators=',nseps
    nsepss(nnull) = nseps
    
    write(12) -1

    close(12)

    write(20,pos=1) nringss,sum(nperring)
    write(20,pos=13) nperring

    close(20)

    deallocate(xs,ys,zs)
    deallocate(line1,line2,remove,add1,add2,endpoints,association,break)

  enddo

  !print*,'number of separators=',int(nsepss)
  
  call SYSTEM_CLOCK(tstop,count_rate)
  print*, 'TIME =',dble(tstop-tstart)/dble(count_rate)

end program
