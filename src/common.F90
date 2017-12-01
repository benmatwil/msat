!common functions, subroutine and variables for ssfind
module common
  use params

  implicit none

  integer(int32) :: nx, ny, nz
  real(np), allocatable :: bgrid(:, :, :, :)
  real(np) :: xmin, xmax, ymin, ymax, zmin, zmax
  real(np), dimension(:), allocatable :: x, y, z
  integer(int32) :: nnulls
  real(np) :: dx, dy, dz
  real(np), allocatable, dimension(:,:) :: rnulls, rnullsalt, rnullsreal
  real(np), allocatable, dimension(:,:) :: spines, fans
  integer(int32), allocatable, dimension(:) :: signs
  
  real(np), parameter :: pi = acos(-1.0_np)
  real(np), parameter :: dtor = pi/180.0_np
  character(100) :: filein, fileout

  contains

    subroutine filenames

      character(100) :: arg, outname
      integer(int32) :: iarg, ic

      ic = 0
      
      if (command_argument_count() > 0) then
        do iarg = 1, command_argument_count()
          call get_command_argument(iarg,arg)
          if (trim(arg) == '-i') then
            call get_command_argument(iarg+1,arg)
            filein = trim(arg)
            ic = 1
          endif
        enddo
      endif
      if (ic == 0) stop 'No input file provided'

      outname = 'output'
      fileout = filein(1:index(filein(1:index(filein, '/', .true.)-1), '/', .true.)) &
        //trim(outname)//'/' &
        //trim(filein(index(filein, '/', .true.)+1:index(filein, '.dat', .true.)-1))
      
      ! if (oc == 0) then
      !   fileout = filein(1:index(filein(1:index(filein, '/', .true.)-1), '/', .true.)) &
      !     //trim(outname)//'/' &
      !     //trim(filein(index(filein, '/', .true.)+1:index(filein, '.dat', .true.)-1))
      ! else
      !   if (index(outname, '/', .true.) /= len(trim(outname))) outname = trim(outname)//'/'
      !   fileout = trim(outname)//trim(filein(index(filein, '/', .true.)+1:index(filein, '.dat', .true.)-1))
      ! endif

    end

    !********************************************************************************

    function trilinear(r, b)
      ! find the value of a function, b, at (x,y,z) using the 8 vertices and trilinear method

      real(np) :: trilinear(3)
      real(np), allocatable :: b(:, :, :, :)
      real(np) :: r(3)
      real(np) :: cube(2, 2, 2), square(2, 2), line(2)
      real(np) :: x, y, z, xp, yp, zp
      integer(int32) :: nx, ny, nz
      integer(int32) :: dims

#if ssf_code
#if cartesian
#else
      call edgecheck(r)
#endif
#endif

      xp = r(1)
      yp = r(2)
      zp = r(3)

      nx = floor(xp)
      ny = floor(yp)
      nz = floor(zp)

#if ssf_code
      ! if point goes out of an edge which is not periodic, set point to be at boundary to trilinear
      ! it will go out and be removed soon
      if (outedge(r)) then
        ! first coordinate needs checking in all coordinate systems
        if (xp <= xmin) then
          xp = xmin
          nx = nint(xp)
        elseif (xp >= xmax) then
          xp = xmax
          nx = nint(xp)-1
        endif
#if cartesian
        ! for cartesians
        if (yp <= ymin) then
          yp = ymin
          ny = nint(yp)
        elseif (yp >= ymax) then
          yp = ymax
          ny = nint(yp)-1
        endif
        if (zp <= zmin) then
          zp = zmin
          nz = nint(zp)
        elseif (zp >= zmax) then
          zp = zmax
          nz = nint(zp)-1
        endif
#elif cylindrical
        ! for cylindricals
        if (zp <= zmin) then
          zp = zmin
          nz = nint(zp)
        elseif (zp >= zmax) then
          zp = zmax
          nz = nint(zp)-1
        endif
#endif
      endif
#else
      if (xp >= xmax) then
        xp = xmax
        nx = nint(xp)-1
      endif
      if (yp >= ymax) then
        yp = ymax
        ny = nint(yp)-1
      endif
      if (zp >= zmax) then
        zp = zmax
        nz = nint(zp)-1
      endif
#endif

      x = xp - nx
      y = yp - ny
      z = zp - nz

      do dims = 1, 3
        cube = b(nx:nx+1, ny:ny+1, nz:nz+1, dims)
        square = (1 - z)*cube(:, :, 1) + z*cube(:, :, 2)
        line = (1 - y)*square(:, 1) + y*square(:, 2)
        trilinear(dims) = (1 - x)*line(1) + x*line(2)
      enddo

    end

    !********************************************************************************

    function dot(a, b)
    ! dot product between a and b

      real(np), dimension(3) :: a, b
      real(np) :: dot

      dot = sum(a*b)

    end

    !********************************************************************************

    function cross(a, b)
    ! cross product between a and b

      real(np), dimension(3) :: a,b
      real(np), dimension(3) :: cross

      cross(1) = a(2)*b(3) - a(3)*b(2)
      cross(2) = a(3)*b(1) - a(1)*b(3)
      cross(3) = a(1)*b(2) - a(2)*b(1)

    end

    !********************************************************************************

    function normalise(a)
    ! finds the unit vector of a

      real(np), dimension(3) :: a
      real(np), dimension(3) :: normalise

      normalise = a/sqrt(dot(a,a))

    end

    !********************************************************************************

    subroutine add_vector(x, vec, pos)
    ! adds an row (vec) to a nx column by ny row array at row number pos

      real(np), allocatable, dimension(:,:) :: x, dummy
      real(np), allocatable, dimension(:) :: vec1
      real(np) :: vec(:)
      integer(int32), optional :: pos
      integer(int32) :: nx, ny, position

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

      x(:,1:position-1) = dummy(:,1:position-1)
      x(:,position) = vec1
      x(:,position+1:ny+1) = dummy(:,position:ny)

    end

    !********************************************************************************

    subroutine add_element(x, number, pos)
      
      integer(int32), allocatable, dimension(:) :: x, dummy
      integer(int32), optional :: pos
      integer(int32) :: nx, position, number, flag
      
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

      x(1:position-1) = dummy(1:position-1)
      x(position) = flag
      x(position+1:nx+1) = dummy(position:nx)

    end

    !********************************************************************************

    subroutine remove_vector(x, pos)
    ! removes row number pos from an array

      real(np), allocatable, dimension(:,:) :: x, dummy
      integer(int32) :: pos, nx, ny

      nx = size(x,1)
      ny = size(x,2)

      dummy = x

      deallocate(x)
      allocate(x(nx,ny-1))

      x(:,1:pos-1) = dummy(:,1:pos-1)
      x(:,pos:ny-1) = dummy(:,pos+1:ny)
      
    end

    !********************************************************************************

    subroutine remove_element(x, pos)
    ! removes row number pos from an array

      integer(int32), allocatable, dimension(:) :: x, dummy
      integer(int32) :: pos, nx

      nx = size(x)

      dummy = x

      deallocate(x)
      allocate(x(nx-1))

      x(1:pos-1) = dummy(1:pos-1)
      x(pos:nx-1) = dummy(pos+1:nx)
      
    end

    !********************************************************************************

    function modulus(a)
      
      real(np), dimension(3) :: a
      real(np) :: modulus

      modulus = sqrt(dot(a,a))

    end function

    !********************************************************************************

    function outedge(r)
    ! determines if the point r is outwith the computation box across a non periodic
    ! boundary and flags true or false

      real(np) :: r(3)
      logical :: outedge
      
      outedge = .false.
      outedge = outedge .or. r(1) >= xmax .or. r(1) <= xmin
#if cartesian
      outedge = outedge .or. r(2) >= ymax .or. r(2) <= ymin .or. r(3) >= zmax .or. r(3) <= zmin
#elif cylindrical
      outedge = outedge .or. r(3) >= zmax .or. r(3) <= zmin
#endif
      
    end function  

    !********************************************************************************

    subroutine edgecheck(r)
    ! determines if the point r is outwith the computational box across a 
    ! periodic boundary and moves if necessary

      real(np) :: r(3)
      
#if cylindrical
      if (r(2) < ymin) r(2) = r(2) + (ymax - ymin)
      if (r(2) > ymax) r(2) = r(2) - (ymax - ymin)
#elif spherical
      if (r(2) < ymin .or. r(2) > ymax) then
        if (r(2) < ymin) r(2) = 2*ymin - r(2)
        if (r(2) > ymax) r(2) = 2*ymax - r(2)
        if (r(3) < (zmax + zmin)/2) then
          r(3) = r(3) + (zmax - zmin)/2
        else
          r(3) = r(3) - (zmax - zmin)/2
        endif
      endif
      if (r(3) <= zmin) r(3) = r(3) + (zmax - zmin)
      if (r(3) >= zmax) r(3) = r(3) - (zmax - zmin)
#endif
          
    end subroutine

    !********************************************************************************

    subroutine get_startpoints(theta,phi,xs,ys,zs)
    ! from the theta,phi coordinates of a fan vector, produces ring of points in the fanplane

        real(np) :: theta, phi
        real(np), dimension(:) :: xs, ys, zs
        integer(int32) :: i, nlines
        real(np) :: dtheta
        real(np) :: r(3)
        real(np), parameter :: sep = 0.05_np

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

    end subroutine

    !********************************************************************************

    function rotate(r,theta,phi)
    ! rotates a vector (r) by a certain angle about the x-z plane (theta),
    ! then by an angle (phi) about the x-y plane

      real(np), dimension(3) :: r, rotate
      real(np) :: theta, phi
      real(np) :: roty(3,3), rotz(3,3), rot(3,3)

      !print *, 'rotate'
      !print*,'r in', r
      !print*, theta/dtor,phi/dtor
      roty(1,1) = cos(-theta)
      roty(1,2) = 0
      roty(1,3) = -sin(-theta)
      
      roty(2,1) = 0
      roty(2,2) = 1
      roty(2,3) = 0
      
      roty(3,1) = sin(-theta)
      roty(3,2) = 0
      roty(3,3) = cos(-theta)

      rotz(1,1) = cos(-phi)
      rotz(1,2) = sin(-phi)
      rotz(1,3) = 0

      rotz(2,1) = -sin(-phi)
      rotz(2,2) = cos(-phi)
      rotz(2,3) = 0

      rotz(3,1) = 0
      rotz(3,2) = 0
      rotz(3,3) = 1

      rot = matmul(rotz,roty)
      rotate = matmul(rot,r)

    end

    !********************************************************************************

    function dist(a, b)
    ! calculates the distance between two points in grid units
      
      real(np) :: dist
      real(np), dimension(3) :: a, b

      dist = sqrt(sum((b - a)**2))

    end function

    !********************************************************************************

    function gtr(r, gx, gy, gz)
    ! converts a list of grid points to real coordinates

      integer(int32) :: ir
      real(np), dimension(:,:), allocatable :: r, gtr
      real(np), dimension(:), allocatable :: gx, gy, gz
      integer(int32), dimension(3) :: ig

      gtr = r

      do ir = 1, size(r, 2)
        ig = floor(r(:, ir))

        gtr(1, ir) = gx(ig(1)) + (r(1, ir) - ig(1))*(gx(ig(1)+1) - gx(ig(1)))
        gtr(2, ir) = gy(ig(2)) + (r(2, ir) - ig(2))*(gy(ig(2)+1) - gy(ig(2)))
        gtr(3, ir) = gz(ig(3)) + (r(3, ir) - ig(3))*(gz(ig(3)+1) - gz(ig(3)))

      enddo

    end function

    !********************************************************************************

    subroutine file_position(nring,nperring,index,a,b,p)
    ! calculates the position in the temporary files of the rings and specific points

      integer(int64) :: a, b, p
      integer(int64) :: uptoring
      integer(int32) :: nring, index
      integer(int32) :: nperring(0:)

      ! 1 whole ring contains 3*np vector points and np association points
      uptoring = sum(int(nperring(0:nring-1), int64))*8
      a = uptoring + int((index-1), int64)
      b = uptoring + int(nperring(nring) + (index-1), int64)
      p = uptoring + int(nperring(nring)*2 + (index-1)*6, int64)

      a = a*4 + 1
      b = b*4 + 1
      p = p*4 + 1

    end subroutine

end module
