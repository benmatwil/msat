!module for the spinefinder
module sf_converge_mod
  use params
  implicit none

  integer :: nx, ny, nz
  double precision, allocatable, dimension(:) :: x, y, z

  double precision, allocatable :: bgrid(:,:,:,:)

  integer, parameter :: nphi = 90
  integer, parameter :: ntheta = nphi/2

  double precision, dimension(3) :: rnull
  double precision :: thetas(ntheta), phis(nphi)

  integer, allocatable :: maxima(:,:), minima(:,:), saddle(:,:)

  contains

  !********************************************************************************

  function trilinear(r,b)
    implicit none
    !find the value of a function, b, at r=(x,y,z) using the 8 vertices
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
    
    nxpoints = size(b,1)
    nypoints = size(b,2)
    nzpoints = size(b,3)

    xp = r(1)
    yp = r(2)
    zp = r(3)

    nx = floor(xp)
    ny = floor(yp)
    nz = floor(zp)
    
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

  function dot(a,b)
    implicit none

    double precision, dimension(3) :: a, b
    double precision :: dot
    
    dot = sum(a*b)
  end function

  !********************************************************************************

  function cross(a,b)
    implicit none

    double precision, dimension(3) :: cross, a, b
    
    cross(1) = a(2)*b(3) - a(3)*b(2)
    cross(2) = a(3)*b(1) - a(1)*b(3)
    cross(3) = a(1)*b(2) - a(2)*b(1)
  end function

  !********************************************************************************

  function modulus(a)
    implicit none

    double precision :: a(3)
    double precision :: modulus

    modulus = sqrt(dot(a,a))

  end

  !********************************************************************************

  function normalise(a)
    implicit none

    double precision :: normalise(3), a(3)

    normalise = a/modulus(a)
  end

  !********************************************************************************

  !take spherical coordinates and convert to Cartesian
  function sphere2cart(r,theta,phi)
    implicit none

    double precision :: theta, phi, r
    double precision :: sphere2cart(3)

    sphere2cart(1) = r*sin(theta)*cos(phi)
    sphere2cart(2) = r*sin(theta)*sin(phi)
    sphere2cart(3) = r*cos(theta)
  end function

  !********************************************************************************

  subroutine remove_vector(x,pos)
    implicit none

    double precision, allocatable, dimension(:,:) :: x, dummy
    integer :: pos
    integer :: nx, ny

    nx = size(x,1)
    ny = size(x,2)

    dummy = x

    deallocate(x)
    allocate(x(nx,ny-1))
    
    x(:,1:pos-1) = dummy(:,1:pos-1)
    x(:,pos:ny-1) = dummy(:,pos+1:ny)
    
  end
  
  !********************************************************************************
  
  subroutine add_element(x,val)
    implicit none

    integer, allocatable, dimension(:,:) :: x, dummy
    integer :: val
    integer :: nx, ny

    nx = size(x,1)
    ny = size(x,2)

    allocate(dummy(nx,ny))

    dummy = x

    deallocate(x)
    allocate(x(nx,ny+1))

    x(:,1:ny) = dummy
    x(:,ny+1) = val

  end

  !********************************************************************************

  subroutine remove_duplicates(vecarray, accur, nclose)
    ! removes vectors with accur distance away from another vector to reduce to "unique" vectors
    implicit none
    
    double precision, allocatable :: vecarray(:,:)
    double precision :: accur
    integer, allocatable, optional :: nclose(:,:)
    integer, allocatable :: dummy(:,:)
    integer :: i, j, n, nclosei, ny
    
    n = size(vecarray,2)
    i = 1
    
    if (present(nclose)) allocate(nclose(1,1))
    
    do while (i < n)
      j = i + 1
      nclosei = 0
      do while (j < n+1)
        if (modulus(vecarray(:,i)-vecarray(:,j)) < accur) then
          call remove_vector(vecarray,j)
          n = size(vecarray,2)
          nclosei = nclosei + 1
        else
          j = j + 1
        endif
      enddo
      i = i + 1
      if (present(nclose)) call add_element(nclose,nclosei)
    enddo
    
    if (present(nclose)) then
      ny = size(nclose,2)
      dummy = nclose
      deallocate(nclose)
      nclose = dummy(:,2:ny)
      deallocate(dummy)
    endif

    end
  
  !********************************************************************************

  subroutine it_conv(rold,rnew,bnew,fact,dir)
    
    implicit none
    
    double precision, dimension(3) :: rold, rnew, bnew
    double precision :: fact
    integer :: dir
    
    rold = rnew
    rnew = rnew + dir*fact*normalise(bnew)
    rnew = rsphere*normalise(rnew)
    bnew = trilinear(rnew+rnull, bgrid)
    
  end
  
end module
