!common functions, subroutine and variables for ssfind
module common
use params

implicit none

real(np), allocatable :: bgrid(:,:,:,:)
real(np) :: xmin, xmax, ymin, ymax, zmin, zmax
real(np), dimension(:), allocatable :: x, y, z
integer :: nnulls
real(np) :: dx, dy, dz
real(np), allocatable, dimension(:,:) :: rnulls, rnullsalt, spines, fans
integer , allocatable, dimension(:) :: signs
integer :: nseps

contains

  function trilinear(r,b)
    ! find the value of a function, b, at (x,y,z) using the 8 vertices
    real(np), allocatable :: b(:,:,:,:)
    real(np) :: r(3)
    real(np) :: cube(2,2,2)
    real(np) :: x, y, z
    real(np) :: xp, yp, zp
    real(np) :: f11, f12, f21, f22
    real(np) :: f1, f2
    real(np) :: trilinear(3)
    integer :: nx, ny, nz
    integer :: dims

    call edgecheck(r)
    
    xp = r(1)
    yp = r(2)
    zp = r(3)

    nx = floor(xp)
    ny = floor(yp)
    nz = floor(zp)
    
    ! if point goes out of an edge which is not periodic, set point to be at boundary to trilinear
    ! it will go out and be removed soon
    if (outedge(r)) then
      ! first coordinate needs checking in all coordinate systems
      if (xp < xmin) then
        xp = xmin
        nx = nint(xp)
      elseif (xp > xmax) then
        xp = xmax
        nx = nint(xp)-1
      endif
#if cartesian
      ! for cartesians
      if (yp < ymin) then
        yp = ymin
        ny = nint(yp)
      elseif (yp > ymax) then
        yp = ymax
        ny = nint(yp)-1
      endif
      if (zp < zmin) then
        zp = zmin
        nz = nint(zp)
      elseif (zp > zmax) then
        zp = zmax
        nz = nint(zp)-1
      endif
#elif cylindrical
      ! for cylindricals
      if (zp < zmin) then
        zp = zmin
        nz = nint(zp)
      elseif (zp > zmax) then
        zp = zmax
        nz = nint(zp)-1
      endif
#endif
    endif

    x = xp-nx
    y = yp-ny
    z = zp-nz
    
    do dims = 1, 3
      cube = b(nx:nx+1,ny:ny+1,nz:nz+1,dims)

      f11 = (1-x)*cube(1,1,1) + x*cube(2,1,1)
      f12 = (1-x)*cube(1,1,2) + x*cube(2,1,2)
      f21 = (1-x)*cube(1,2,1) + x*cube(2,2,1)
      f22 = (1-x)*cube(1,2,2) + x*cube(2,2,2)

      f1 = (1-y)*f11 + y*f21
      f2 = (1-y)*f12 + y*f22
      
      trilinear(dims) = (1-z)*f1 + z*f2
    enddo

  end

  !********************************************************************************

  ! dot product between a and b
  function dot(a,b)

    real(np), dimension(3) :: a, b
    real(np) :: dot

    dot = sum(a*b)

  end

  !********************************************************************************

  ! cross product between a and b
  function cross(a,b)

    real(np), dimension(3) :: a,b
    real(np), dimension(3) :: cross

    cross(1) = a(2)*b(3) - a(3)*b(2)
    cross(2) = a(3)*b(1) - a(1)*b(3)
    cross(3) = a(1)*b(2) - a(2)*b(1)

  end

  !********************************************************************************

  ! normalises a
  function normalise(a)

    real(np), dimension(3) :: a
    real(np), dimension(3) :: normalise

    normalise = a/sqrt(dot(a,a))

  end

  !********************************************************************************

  ! adds an row (val) to a nx column by ny row array at row number pos
  subroutine add_vector(x,vec,pos)
    implicit none
    
    real(np), allocatable, dimension(:,:) :: x, dummy
    real(np), allocatable, dimension(:) :: vec1
    real(np) :: vec(:)
    integer, optional :: pos
    integer :: nx, ny, position

    vec1 = vec

    nx = size(x,1)
    ny = size(x,2)
    
    if (present(pos)) then
      position = pos
    else
      position = ny+1
    endif

    dummy = x

    deallocate(x)
    allocate(x(nx,ny+1))

    ! if (position == 1) then
    !   x(:,1) = vec1
    !   x(:,2:ny+1) = dummy
    ! elseif (position == ny+1) then
    !   x(:,1:ny) = dummy
    !   x(:,ny+1) = vec1
    ! else
      x(:,1:position-1) = dummy(:,1:position-1)
      x(:,position) = vec1
      x(:,position+1:ny+1) = dummy(:,position:ny)
    ! endif

  end

  !********************************************************************************

  subroutine add_element(x,number,pos)
    implicit none
    
    integer, allocatable, dimension(:) :: x, dummy
    integer, optional :: pos
    integer :: nx, position, number, flag
    
    flag = number

    nx = size(x,1)

    if (present(pos)) then
      position = pos
    else
      position = nx+1
    endif

    dummy = x

    deallocate(x)
    allocate(x(nx+1))

    ! if (position == 1) then
    !   x(1) = flag
    !   x(2:nx+1) = dummy
    ! elseif (position == nx+1) then
    !   x(1:nx) = dummy
    !   x(nx+1) = flag
    ! else
      x(1:position-1) = dummy(1:position-1)
      x(position) = flag
      x(position+1:nx+1) = dummy(position:nx)
    ! endif

  end

  !********************************************************************************

  ! removes row number pos from an array 
  subroutine remove_vector(x,pos)
    implicit none
    
    real(np), allocatable, dimension(:,:) :: x, dummy
    integer :: pos, nx, ny

    nx = size(x,1)
    ny = size(x,2)

    dummy = x

    deallocate(x)
    allocate(x(nx,ny-1))

    ! if (pos == 1) then
    !   x = dummy(:,2:ny)
    ! elseif (pos == ny) then
    !   x = dummy(:,1:ny-1)
    ! else
      x(:,1:pos-1) = dummy(:,1:pos-1)
      x(:,pos:ny-1) = dummy(:,pos+1:ny)
    ! endif
    
  end

  !********************************************************************************

  ! removes row number pos from an array 
  subroutine remove_element(x,pos)
    implicit none
    
    integer, allocatable, dimension(:) :: x, dummy
    integer :: pos, nx

    nx = size(x)

    dummy = x

    deallocate(x)
    allocate(x(nx-1))

    ! if (pos == 1) then
    !   x = dummy(2:nx)
    ! elseif (pos == nx) then
    !   x = dummy(1:nx-1)
    ! else
      x(1:pos-1) = dummy(1:pos-1)
      x(pos:nx-1) = dummy(pos+1:nx)
    ! endif
    
  end

  !********************************************************************************

  ! finds |a|
  function modulus(a)
    real(np), dimension(3) :: a
    real(np) :: modulus

    modulus = sqrt(dot(a,a))

  end function

  !********************************************************************************

  function outedge(r)
    real(np) :: r(3)
    logical :: outedge
    
    outedge = .false.  
    if (r(1) > xmax .or. r(1) < xmin) outedge = .true.
#if cartesian
    if (r(2) > ymax .or. r(2) < ymin) outedge = .true.
    if (r(3) > zmax .or. r(3) < zmin) outedge = .true.
#elif cylindrical
    if (r(3) > zmax .or. r(3) < zmin) outedge = .true.
#endif
    
  end function  

  !********************************************************************************

  ! determines if the vector r is outwith the computational box
  subroutine edgecheck(r, out)
    real(np) :: r(3)
    logical, optional :: out
    
    if (present(out)) out = outedge(r)
    
#if cylindrical
    if (r(2) < ymin) r(2) = r(2) + ymin - ymax
    if (r(2) > ymin) r(2) = r(2) - (ymin - ymax)
#elif spherical
    if (r(2) < ymin .or. r(2) > ymax) then
      if (r(2) < ymin) r(2) = 2*ymin - r(2)
      if (r(2) > ymax) r(2) = 2*ymax - r(2)
      if (r(3) < (zmax+zmin)/2) then
        r(3) = r(3) + (zmax-zmin)/2
      else
        r(3) = r(3) - (zmax-zmin)/2
      endif
    endif
    if (r(3) < zmin) r(3) = r(3) + zmax - zmin
    if (r(3) > zmax) r(3) = r(3) - (zmax - zmin)
#endif
        
  end subroutine

  !********************************************************************************

  ! from the theta,phi coordinates of a fan vector, produces ring of points in the fanplane
  subroutine get_startpoints(theta,phi,xs,ys,zs)

      real(np) :: theta, phi
      real(np), dimension(:) :: xs, ys, zs
      integer :: i, nlines
      real(np) :: dtheta
      real(np) :: r(3)
      real(np), parameter :: sep=0.1

      nlines = size(xs)

      dtheta = 2.0_np*pi/nlines

      ! generate nlines start points in ring around equator
      do i = 1, nlines
        r(1) = cos(i*dtheta)
        r(2) = sin(i*dtheta)
        r(3) = 0

        r = rotate(r,theta,phi)

        xs(i) = r(1)*sep
        ys(i) = r(2)*sep
        zs(i) = r(3)*sep
      enddo

      r(1) = sin(theta)*cos(phi)
      r(2) = sin(theta)*sin(phi)
      r(3) = cos(theta)

  end subroutine

  !********************************************************************************

  ! rotates a vector (r) by a certain angle about the x-z plane (theta), then by an angle (phi) about the x-y plane
  function rotate(r,theta,phi)
      real(np) :: r(3)
      real(np) :: theta, phi
      real(np) :: rotate(3)
      real(np) :: roty(3,3), rotz(3,3), rot(3,3)

      !print *, 'rotate'
      !print*,'r in', r
      !print*, theta/dtor,phi/dtor
      roty(1,1)=cos(-theta)
      roty(1,2)=0
      roty(1,3)=-sin(-theta)
      
      roty(2,1)=0
      roty(2,2)=1
      roty(2,3)=0
      
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

      rot=matmul(rotz,roty)
      rotate=matmul(rot,r)

      r=rotate

  end

  !********************************************************************************

  function dist(a,b)
    ! calculates the distance between two points in grid units
    
    real(np) :: dist
    real(np), dimension(3) :: a, b

    dist = sqrt((a(1)-b(1))**2 + (a(2)-b(2))**2 + (a(3)-b(3))**2)

  end function

end module
