!code to find the separatrix surface and separators of null points

program ssfind

  use omp_lib
  use params
  use common
  use trace
  use ring

  implicit none

  !integer :: ierror

  character (len=8), parameter :: fmt='(I3.3)'
  character (len=5) :: fname

  integer :: nringss

  integer*8 :: tstart,tstop,count_rate
  !b field readin variables
  integer :: nx, ny, nz

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

  !integer :: nnulls
  !double precision, allocatable, dimension(:,:) :: rnulls

  !integer :: ringnum,index,nullnum
  integer :: nrings, npoints
  integer :: nperring(ringsmax)
  double precision :: circumference(ringsmax)

  logical :: exitcondition, out

  CALL OMP_SET_NUM_THREADS(12) !have it work on 4 threads (If machine has >4 cores this should be larger, if fewer than 4 coures, this should be smaller)
  
  !read in data
  !open (unit=10,file='output/2null-edge.dat',form='unformatted')
  open (unit=10,file=filename,access='stream')

    read(10), nx, ny, nz !number of vertices
    print*, nx, ny, nz
    
    if (1 == 0) then
      allocate(bgrid(0:nx+1,0:ny+1,0:nz+1,3))
      bgrid = 0
      allocate(x(0:nx+1), y(0:ny+1), z(0:nz+1))
      read(10), bgrid(1:nx,1:ny,1:nz,1)
      read(10), bgrid(1:nx,1:ny,1:nz,2)
      read(10), bgrid(1:nx,1:ny,1:nz,3)
      read(10) x(1:nx), y(1:ny), z(1:nz)

      !bgrid = bgrid(:,ny+1:0:-1,:,:)
      !y = y(ny+1:0:-1)
      
      !copy the grid over for periodicity of coordinate systems
      x(0) = 2*x(1) - x(2) !x,R,r don't repeat
      x(nx+1) = 2*x(nx)-x(nx-1)
      if (coord_type == 3) then !cylindrical
        !phi repeating
        bgrid(:,0,:,:) = bgrid(:,ny,:,:)
        bgrid(:,ny+1,:,:) = bgrid(:,1,:,:)
        y(0) = y(1) - (y(ny) - y(ny-1))
        y(ny+1) = y(ny) + (y(1) - y(0))
        z(0) = 2*z(1) - z(2)
        z(nz+1) = 2*z(nz)-z(nz-1)
      else if (coord_type == 2) then !spherical
        !phi repeating
        bgrid(:,:,0,:) = bgrid(:,:,nz,:)
        bgrid(:,:,nz+1,:) = bgrid(:,:,1,:)
        !have theta "periodic"
        if (2*(nz/2) == nz) then !if nz (nphi) is even
          !Top boundary
          bgrid(:,0,1:nz/2,:) = bgrid(:,2,nz/2+1:nz,:)
          bgrid(:,0,nz/2+1:nz,:) = bgrid(:,2,0:nz/2,:)
          !Bottom boundary
          bgrid(:,ny+1,1:nz/2,:) = bgrid(:,ny-1,nz/2+1:nz,:)
          bgrid(:,ny+1,nz/2+1:nz,:) = bgrid(:,ny-1,0:nz/2,:)
        else !if nz is odd
          !Top boundary
          bgrid(:,0,1:nz/2,:) = (bgrid(:,2,1+nz/2:nz-1,:) + bgrid(:,2,2+nz/2:nz,:))/2
          bgrid(:,0,nz/2+1,:) = (bgrid(:,2,nz,:) + bgrid(:,2,1,:))/2
          bgrid(:,0,nz/2+2:nz,:) = (bgrid(:,2,1:nz/2,:) + bgrid(:,2,2:nz/2+1,:))/2
          !Bottom boundary
          bgrid(:,ny+1,1:nz/2,:) = (bgrid(:,ny-1,1+nz/2:nz-1,:) + bgrid(:,ny-1,2+nz/2:nz,:))/2
          bgrid(:,ny+1,nz/2+1,:) = (bgrid(:,ny-1,nz,:) + bgrid(:,ny-1,1,:))/2
          bgrid(:,ny+1,nz/2+2:nz,:) = (bgrid(:,ny-1,1:nz/2,:) + bgrid(:,ny-1,2:nz/2+1,:))/2
        endif
        z(0) = 2*z(1) - z(2)
        z(nz+1) = 2*z(nz)-z(nz-1)
        y(0) = y(2)
        y(ny+1) = y(ny-1)
      else if (coord_type == 1) then
        z(0) = 2*z(1) - z(2)
        z(nz+1) = 2*z(nz)-z(nz-1)
      endif
    else
      allocate(bgrid(nx,ny,nz,3))
      allocate(x(nx), y(ny), z(nz))
      read(10), bgrid(:,:,:,1)
      read(10), bgrid(:,:,:,2)
      read(10), bgrid(:,:,:,3)
      read(10) x, y, z
    endif

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
      rnullsalt(3,i) = rnullsalt(3,i) - 1
      call edgecheck(rnullsalt(3,i))
      rnullsalt(3,i) = rnullsalt(3,i) + 1
      
      !check whether null is at the upper phi boundary
      rnullsalt(3,i) = rnullsalt(3,i) + 1
      call edgecheck(rnullsalt(3,i))
      rnullsalt(3,i) = rnullsalt(3,i) - 1
      
      !check whether null is at the lower theta boundary
      rnullsalt(2,i) = rnullsalt(2,i) - 1
      call edgecheck(rnullsalt(2,i))
      rnullsalt(2,i) = rnullsalt(2,i) - 1
      
      !check whether null is at the upper theta boundary
      rnullsalt(2,i) = rnullsalt(2,i) + 1
      call edgecheck(rnullsalt(2,i))
      rnullsalt(2,i) = rnullsalt(2,i) + 1
    enddo
  endif
  
  allocate(nsepss(nnulls))

  !signs=-1*signs
  do nnull = 45,45!1, nnulls !loop over all nulls
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
    allocate(association(1,nlines), break(1,nlines), remove(1,nlines), endpoints(1,nlines))
    add1 = 0
    add2 = 0

    !add start points to first ring relative to null
    do j = 1, nlines !Go through each start point
      line1(1,j) = r(1) + xs(j)
      line1(2,j) = r(2) + ys(j)
      line1(3,j) = r(3) + zs(j)
      association(1,j) = j
    enddo
    
    line2 = line1

    break = 0

    write(fname,fmt) nnull

    open(unit=12,file='output/separator'//trim(fname)//'.dat',form='unformatted')

    open(unit=20,file='output/everything'//trim(fname)//'.dat',access='stream',status='replace')

    nseps = 0

    nrings = 0

    nperring = 0

    write(20) nrings, npoints, ringsmax
    write(20) nperring
    
    circumference = 0

    exitcondition = .false.

    nringss = 0

    !$OMP PARALLEL DEFAULT(SHARED)
    do i = 1, ringsmax !loop over number of rings we want
      if (sign .eq. 0) then !skip null which is uncharacterised
        !$OMP SINGLE
          print*,'Null has zero sign'
        !$OMP END SINGLE
        exit
      endif
      !$OMP SINGLE
        nringss = nringss+1

        endpoints = 0
        add1 = 0
        add2 = 0
        remove = 0

        write(20) association
        nperring(i) = nlines

        ierror = 0
        
        if (i < 50) then
          h = 5d-2
        else
          h = 25d-2
        endif
      !$OMP END SINGLE

      !$OMP DO private(r), firstprivate(h)
        do j = 1, nlines !loop over all points in ring (in parallel do)

          r(:) = line1(:,j)

          call trace_line(r,1,sign,h) !trace line by a distance of h

          line2(:,j) = line2(:,j) + r(:) - line1(:,j)
          call edgecheck(r, out)
          if (out) then !counter to see how many points on ring have reached outer boundary
            endpoints(1,j) = 1
          else
            endpoints(1,j) = 0
          endif
          line1(:,j) = r(:)
          
          association(1,j) = j

        enddo
      !$OMP END DO

      !$OMP SINGLE
        write(20), line1

        if (nlines > pointsmax) then !exit if too many points on ring
          print*, 'Too many points on ring', nlines, i
          exitcondition = .true.
        endif

        if (sum(endpoints)/nlines > 0.99) then !exit if all points have reached outer boundary (left box)
          print*, 'All fan points have reached the outer boundary', i
          exitcondition = .true.
        endif
        
        circumference(i) = circumference(i) + dist(line2(:,1),line2(:,nlines))
        do j = 2, nlines
          circumference(i) = circumference(i) + dist(line2(:,j),line2(:,j-1))
        enddo
        
        if (i > 1)  then
          if (abs(circumference(i)-circumference(i-1)) .lt. 0.1*stepmin) then
            print*, 'Fan has stopped growing/shrinking', i
            exitcondition = .true.
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

      !$OMP END SINGLE

      if (exitcondition) exit

    enddo
    !$OMP END PARALLEL
    
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
  print*, 'TIME=',dble(tstop-tstart)/dble(count_rate)

end program
