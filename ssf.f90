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
  double precision, allocatable, dimension(:,:,:) :: bx, by, bz
  integer :: nx, ny, nz

  !position vector of null (and saved backup)
  double precision :: r(3)
  double precision :: a(3)
  double precision :: h

  integer :: i,j

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
  integer :: nrings,npoints
  integer :: nperring(ringsmax)
  double precision :: circumference(ringsmax)

  logical :: exitcondition

  CALL OMP_SET_NUM_THREADS(4) !have it work on 4 threads (If machine has >4 cores this should be larger, if fewer than 4 coures, this should be smaller)

  !read in data
  !open (unit=10,file='output/2null-edge.dat',form='unformatted')
  open (unit=10,file=filename,access='stream')

    read(10),nx,ny,nz !number of vertices

    allocate(bgrid(nx,ny,nz,3))
    !allocate(bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz))
    allocate(x(nx), y(ny), z(nz))
    read(10),bgrid(:,:,:,1)
    read(10),bgrid(:,:,:,2)
    read(10),bgrid(:,:,:,3)

    print*,nx,ny,nz

    read(10) x, y, z

    xmin=1.
    xmax=nx
    ymin=1.
    ymax=ny
    zmin=1.
    zmax=nz

    dx=(x(nx)-x(1))/nx
    dy=(y(ny)-y(1))/ny
    dz=(z(nz)-z(1))/nz

  ! deallocate(bx,by,bz)
  close(10)

  print*, 'Data read in'

  call SYSTEM_CLOCK(tstart,count_rate) !to time how long it takes

  !read in null data

  open (unit=10,file='output/nulls.dat',form='unformatted')
    read (10) nnulls
    allocate(signs(nnulls),rnulls(3,nnulls),spines(3,nnulls),fans(3,nnulls))
    read(10) signs
    read(10) rnulls,spines,fans
  close(10)

  allocate(nsepss(nnulls))

  !signs=-1*signs
  do nnull=1,nnulls !loop over all nulls
    print*,''
    print*,'Null number',nnull,'of',nnulls

    r=rnulls(:,nnull)
    spine=spines(:,nnull)
    fan=fans(:,nnull)
    sign=signs(nnull)

    a=r !backup null location

    !print*, rnulls
    !print*, r
    !print*,''

    !get start points
    nlines=nstart

    theta=acos(fan(3))
    phi=atan(fan(2),fan(1))

    !print*,'theta=',theta/dtor
    !print*,'phi=',phi/dtor

    allocate(xs(nlines),ys(nlines),zs(nlines))

    call get_startpoints(theta,phi,xs,ys,zs)

    allocate(line(3,nlines),remove(3,nlines), add(3,nlines),endpoints(3,nlines))
    allocate(association(1,nlines),break(1,nlines))
    add=0.

    !add start points to first ring
    do j=1,nlines !start point
      line(1,j)=r(1)+xs(j)
      line(2,j)=r(2)+ys(j)
      line(3,j)=r(3)+zs(j)
      association(:,j)=j
    enddo

    break=0.

    !open (unit=10,file='output/ring.dat',form='unformatted')
    !open (unit=11,file='output/association.dat',form='unformatted')

    write(fname,fmt) nnull

    open (unit=12,file='output/separator'//trim(fname)//'.dat',form='unformatted')

    open (unit=20,file='output/everything'//trim(fname)//'.dat',access='stream',status='replace')

    nseps=0

    nrings=0

    nperring=0

    write(20) nrings, npoints,ringsmax
    write(20) nperring
    
    circumference=0.

    exitcondition=.false.

    nringss=0

    !$OMP PARALLEL DEFAULT(SHARED)
    do i=1,ringsmax !loop over number of rings we want
        if (sign .eq. 0) then !skip null which is uncharacterised
        !$OMP SINGLE
        print*,'Null has zero sign'
        !$OMP END SINGLE
        exit
        endif
        !$OMP SINGLE
        nringss=nringss+1

        endpoints=0.
        add=0.
        remove=0.

        write(20) nint(association)
        nperring(i)=nlines

        ierror=0
        !call flush()
      !$OMP END SINGLE

      !$OMP DO  private(r,h)!, shared(ierror)
        do j=1,nlines !loop over all points in ring (in parallel do)

          r(:) = line(:,j)

          if (i .lt. 50) then
            h=0.05
          else
            h=0.5
          endif

          call trace_line(r,1,sign,h) !trace line by a distance of h

          if (edge(r)) then !counter to see how many points on ring have reached outer boundary
            endpoints(:,j)=1.
          else
            endpoints(:,j)=0.
          endif
          ! print*,h

          line(:,j)=r(:)

          association(:,j) = dble(j)

        enddo
      !$OMP END DO

      !stop

      !$OMP SINGLE
        write(20),line

        if (nlines .gt. pointsmax) then !exit if too many points on ring
        print*,'Too many points on ring',nlines,i
        exitcondition = .true.
        endif

        circumference(i)=circumference(i)+dist(line(:,1),line(:,nlines))
        do j=2,nlines
        circumference(i)=circumference(i)+dist(line(:,j),line(:,j-1))
        enddo

        !print*, 'ring',i,nlines!,nint(sum(endpoints(1,:))),sum(endpoints(1,:))/nlines,circumference(i)

        if (sum(endpoints(1,:))/nlines .gt. 0.99) then !exit if all points have reached outer boundary (left box)
        print*, 'All fan points have reached the outer boundary',i
        exitcondition=.true.
        endif
        if (i .gt. 1)  then
        if(abs(circumference(i)-circumference(i-1)) .lt. 0.1*stepmin) then
          print*, 'Fan has stopped growing',i
          exitcondition=.true.
        endif
        endif

        if (ierror .eq. 1) then
        print*, 'Tracing has failed',i
        exitcondition = .true.
        endif

        call at_null(nlines,nnull,i) !determine if point is at null

        call add_points(nlines,i) !add points to ring if necessary

        call remove_points(nlines,i) !remove points from ring if necessary

        !call flush()
      !$OMP END SINGLE

        if (exitcondition) exit

    enddo
    !$OMP END PARALLEL
    
    print*, 'number of separators=',nseps
    nsepss(nnull)=nseps
    
    write(12) -1

    close(12)

    write(20,pos=1) nringss,sum(nperring)
    write(20,pos=13) nperring

    close(20)

    deallocate(xs,ys,zs)
    deallocate(line,remove,add,endpoints,association,break)

  enddo

  call SYSTEM_CLOCK(tstop,count_rate)

  print*,'number of separators=',int(nsepss)

  print*, 'TIME=',dble(tstop-tstart)/dble(count_rate)

end program
