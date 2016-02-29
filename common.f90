!common functions, subroutine and variables for ssfind
module common
use params



double precision, allocatable :: bgrid(:,:,:,:)
integer, parameter :: iseed=3141592
double precision :: xmin, xmax, ymin, ymax, zmin, zmax
double precision, dimension(:), allocatable :: x, y, z
integer :: nnulls
double precision :: dx, dy, dz
double precision, allocatable, dimension(:,:) :: rnulls, spines, fans
integer , allocatable, dimension(:) :: signs
integer :: nseps

integer :: ierror



 contains




function trilinear(r,b)
    !find the value of a function, b, at (x,y,z) using the 8 vertices
    double precision :: b(:,:,:,:)
    double precision :: r(3)
    double precision :: cube(2,2,2)
    double precision :: x, y, z
    double precision :: xp, yp, zp
    double precision :: f11, f12, f21, f22
    double precision :: f1, f2
    double precision :: trilinear(3)
    integer :: nx, ny, nz
    integer :: dims
    integer :: nxpoints,nypoints,nzpoints

    nxpoints=size(bgrid,1)
    nypoints=size(bgrid,2)
    nzpoints=size(bgrid,3)

    xp=r(1)
    yp=r(2)
    zp=r(3)

    nx=floor(xp)
    ny=floor(yp)
    nz=floor(zp)



    x=xp-nx
    y=yp-ny
    z=zp-nz

    if (xp .ge. nxpoints) then
      x=nxpoints
      nx=nxpoints-1
    endif
    if (xp .le. 1.) then
      x=1.
      nx=1
    endif
    if (yp .ge. nypoints) then
      y=nypoints
      ny=nypoints-1
    endif
    if (yp .le. 1.) then
      y=1.
      ny=1
    endif
    if (zp .ge. nzpoints) then
      z=nzpoints
      nz=nzpoints-1
    endif
    if (zp .le. 1.) then
      z=1.
      nz=1
    endif


    do dims=1,3
	cube=b(nx:nx+1,ny:ny+1,nz:nz+1,dims)

	f11=(1-x)*cube(1,1,1) + x*cube(2,1,1)
	f12=(1-x)*cube(1,1,2) + x*cube(2,1,2)
	f21=(1-x)*cube(1,2,1) + x*cube(2,2,1)
	f22=(1-x)*cube(1,2,2) + x*cube(2,2,2)

	f1=(1-y)*f11 + y*f21
	f2=(1-y)*f12 + y*f22

	trilinear(dims) = (1-z)*f1 + z*f2
    enddo
end







!dot product between a and b
function dot(a,b)
  double precision, dimension(3) :: a,b
  double precision :: dot

  dot = sum(a*b)
end







!cross product between a and b
function cross(a,b)
  double precision, dimension(3) :: a,b
  double precision, dimension(3) :: cross

  cross(1) = a(2)*b(3) - a(3)*b(2)
  cross(2) = a(3)*b(1) - a(1)*b(3)
  cross(3) = a(1)*b(2) - a(2)*b(1)
end





!normalises a
function normalise(a)
  double precision, dimension(3) :: a
  double precision, dimension(3) :: normalise

  normalise = a/sqrt(dot(a,a))
end







!adds an row (val) to a nx column by ny row array at row number pos
subroutine add_element(x,val,pos)
    double precision, allocatable, dimension(:,:) :: x
    double precision :: val(:)
    double precision, allocatable, dimension(:,:) :: dummy
    integer :: pos
    integer :: nx, ny

    nx=size(x,1)
    ny=size(x,2)


    allocate(dummy(nx,ny))

    dummy=x

    deallocate(x)
    allocate(x(nx,ny+1))

    if (pos .eq. 1) then
      x(:,1)=val
      x(:,2:ny+1)=dummy
    else if (pos .eq. ny+1 .or. pos .lt. 1) then
      x(:,1:ny)=dummy
      x(:,ny+1)=val
    else
      x(:,1:pos-1) = dummy(:,1:pos-1)
      x(:,pos)=val
      x(:,pos+1:ny+1)=dummy(:,pos:ny)
    endif

    deallocate(dummy)

end

!removes row number pos from an array 
subroutine remove_element(x,pos)
    double precision, allocatable, dimension(:,:) :: x
    double precision, allocatable, dimension(:,:) :: dummy
    integer :: pos
    integer :: nx, ny

    nx=size(x,1)
    ny=size(x,2)


    allocate(dummy(nx,ny))

    dummy=x

    deallocate(x)
    allocate(x(nx,ny-1))

    if (pos .eq. 1) then
      x(:,1:ny-1)=dummy(:,2:ny)
    else if (pos .eq. ny) then
      x(:,1:ny-1)=dummy(:,1:ny-1)
    else
      x(:,1:pos-1)=dummy(:,1:pos-1)
      x(:,pos:ny-1)=dummy(:,pos+1:ny)
    endif

    deallocate(dummy)

end



!finds |a|
function modulus(a)
  double precision, dimension(3) :: a
  double precision :: modulus

  modulus=sqrt(dot(a,a))

end function




!determines if the vector r is outwith the computational box
function edge(r)
    double precision :: r(3)
    logical :: edge
    edge = .false.
    if (r(1) .gt. xmax .or. r(1) .lt. xmin) edge = .true.
    if (r(2) .gt. ymax .or. r(2) .lt. ymin) edge = .true.
    if (r(3) .gt. zmax .or. r(3) .lt. zmin) edge = .true.
end function









!from the theta,phi coordinates of a fan vector, produces ring of points in the fanplane
subroutine get_startpoints(theta,phi,xs,ys,zs)

    double precision :: theta, phi
    double precision, dimension(:) :: xs, ys, zs
    integer :: nlines
    double precision :: dtheta
    double precision :: r(3)
    double precision, parameter :: sep=0.1


    nlines=size(xs)

    dtheta=2.*pi/nlines

    !generate nlines start points in ring around equator
    do i=1,nlines
      r(1)=cos(i*dtheta)
      r(2)=sin(i*dtheta)
      r(3)=0

      r=rotate(r,theta,phi)

      xs(i)=r(1)*sep
      ys(i)=r(2)*sep
      zs(i)=r(3)*sep
    enddo



    r(1)=sin(theta)*cos(phi)
    r(2)=sin(theta)*sin(phi)
    r(3)=cos(theta)

    !print*,''
    !print*,'point number, dot product, distance from null'
    !do i=1,nlines
    !  print*,i,r(1)*xs(i)+r(2)*ys(i)+r(3)*zs(i),sqrt(xs(i)**2+ys(i)**2+zs(i)**2)
    !enddo


end subroutine



!rotates a vector (r) by a certain angle about the x-z plane (theta), then by an angle (phi) about the x-y plane
function rotate(r,theta,phi)
    double precision :: r(3)
    double precision :: theta, phi
    double precision :: rotate(3)
    double precision :: roty(3,3), rotz(3,3), rot(3,3)

    !print *, 'rotate'
    !print*,'r in', r
    !print*, theta/dtor,phi/dtor
    roty(1,1)=cos(-theta)
    roty(1,2)=0
    roty(1,3)=-sin(-theta)
    !
    roty(2,1)=0
    roty(2,2)=1
    roty(2,3)=0
    !
    roty(3,1)=sin(-theta)
    roty(3,2)=0
    roty(3,3)=cos(-theta)

    rotz(1,1)=cos(-phi)
    rotz(1,2)=sin(-phi)
    rotz(1,3)=0

    rotz(2,1)=-sin(-phi)
    rotz(2,2)=cos(-phi)
    rotz(2,3)=0

    rotz(3,1)=0
    rotz(3,2)=0
    rotz(3,3)=1

    !r=matmul(roty,r)
    !print*,r

    !r=matmul(rotz,r)
    !print*,r

    rot=matmul(rotz,roty)
    rotate=matmul(rot,r)
    r=rotate
    !print*,'r out', r
    !print*, 'theta, phi out', acos(r(3))/dtor,atan(r(2),r(1))/dtor!+180
    !rotate=r

    !print*,roty
    !print*, rotz

end



function dist(a,b)
    !calculates the distance between two points in grid units
    double precision :: dist
    double precision, dimension(3) :: a,b

    !dist=sqrt(dot(a-b,a-b))
    dist=sqrt((a(1)-b(1))*(a(1)-b(1)) + (a(2)-b(2))*(a(2)-b(2)) + (a(3)-b(3))*(a(3)-b(3)))
    end function


end module
