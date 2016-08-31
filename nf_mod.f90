!module for nullfinder
module nf_mod
use params

contains

  subroutine bilin_test(cbx,cby,cbz,inull)
    !tests surface of cube using bilinear method
    implicit none

    integer :: i, inull, n
    !integer :: itest1, itest2, itest3
    !integer :: face
    !integer :: idum
    integer :: cross, sign
    integer :: edge

    double precision, dimension(2,2,2) :: cbx,cby,cbz

    double precision, dimension(2,2) :: facex, facey, facez

    integer, dimension(6,6) :: test
    integer, dimension(6) :: test2
    integer, dimension(3) :: test3
    !test is an array where the first index refers to the face number (1:6) and the second index is follows:
    !1- x-y intersection (1/0)
    !2- sign of z at this intersection (-1/1 /10 if null on face)
    !3- xz intersection (1/0)
    !4- sign of y on this intersection (-1/1 /10 if null on face)
    !5- yz intersection (1/0)
    !6- sign on x on this inersection (-1/1 /10 if null on face)

    inull = 0 !integer variable defining null. 1=null, 0=no null
    test = 0 !set test array to zero
    edge = 0 !integer variable defining null on edge. 1=null, 0=no null

    !face #1 : (bottom face)
    facex = cbx(:,:,1)
    facey = cby(:,:,1)
    facez = cbz(:,:,1)

    call face_solve(facex, facey, facez, cross, sign) !check if zeroes of x and y components of B cross on face
    test(1,1) = cross !do they cross (1/0)
    test(1,2) = sign !sign (-1/1)
    call face_solve(facey, facez, facex, cross, sign) !check if zeroes of y and z components of B cross on face
    test(1,3) = cross !do they cross (1/0)
    test(1,4) = sign !sign (-1/1)
    call face_solve(facex, facez, facey, cross, sign) !check if zeroes of x and z components of B cross on face
    test(1,5) = cross !do they cross (1/0)
    test(1,6) = sign !sign (-1/1)

    edge = edge + edge_check(facex, facey, facez) !check for null on edge (with separator lying along edge)

    !face #2 : left face
    facex = cbx(1,:,:)
    facey = cby(1,:,:)
    facez = cbz(1,:,:)

    call face_solve(facex, facey, facez, cross, sign)
    test(2,1) = cross
    test(2,2) = sign
    call face_solve(facey, facez, facex, cross, sign)
    test(2,3) = cross
    test(2,4) = sign
    call face_solve(facex, facez, facey, cross, sign)
    test(2,5) = cross
    test(2,6) = sign

    edge = edge + edge_check(facex, facey, facez)

    !face #3 : right face
    facex = cbx(2,:,:)
    facey = cby(2,:,:)
    facez = cbz(2,:,:)

    call face_solve(facex, facey, facez, cross, sign)
    test(3,1) = cross
    test(3,2) = sign
    call face_solve(facey, facez, facex, cross, sign)
    test(3,3) = cross
    test(3,4) = sign
    call face_solve(facex, facez, facey, cross, sign)
    test(3,5) = cross
    test(3,6) = sign

    edge = edge + edge_check(facex, facey, facez)

    !face 4: front face
    facex = cbx(:,1,:)
    facey = cby(:,1,:)
    facez = cbz(:,1,:)

    call face_solve(facex, facey, facez, cross, sign)
    test(4,1) = cross
    test(4,2) = sign
    call face_solve(facey, facez, facex, cross, sign)
    test(4,3) = cross
    test(4,4) = sign
    call face_solve(facex, facez, facey, cross, sign)
    test(4,5) = cross
    test(4,6) = sign

    edge = edge + edge_check(facex, facey, facez)

    !face 5: rear face
    facex = cbx(:,2,:)
    facey = cby(:,2,:)
    facez = cbz(:,2,:)

    call face_solve(facex, facey, facez, cross, sign)
    test(5,1) = cross
    test(5,2) = sign
    call face_solve(facey, facez, facex, cross, sign)
    test(5,3) = cross
    test(5,4) = sign
    call face_solve(facex, facez, facey, cross, sign)
    test(5,5) = cross
    test(5,6) = sign

    edge = edge + edge_check(facex, facey, facez)


    !face 6: top face
    facex = cbx(:,:,2)
    facey = cby(:,:,2)
    facez = cbz(:,:,2)

    call face_solve(facex, facey, facez, cross, sign)
    test(6,1) = cross
    test(6,2) = sign
    call face_solve(facey, facez, facex, cross, sign)
    test(6,3) = cross
    test(6,4) = sign
    call face_solve(facex, facez, facey, cross, sign)
    test(6,5) = cross
    test(6,6) = sign

    edge = edge + edge_check(facex, facey, facez)

    !print out test array for checking
    !do i=1,6
      !print*,test(i,:)
    !enddo

    !check for nulls within cell

    !total of each column in 'test'
    test2=sum(test,1)

    test3=0
    do i = 1, 3
      if (test2(2*i-1) >= 0) then !if there is a crossing
        if (abs(test2(2*i)) < test2(2*i-1)) test3(i) = 1
      endif
    enddo

    if (sum(test3) == 3) inull = 1 !we have a null (in the cell)

    ! check for null on surface of cell
    if (sum(test(:,2)) > 50 .or. sum(test(:,4)) > 50 .or. sum(test(:,6)) > 50) then
      !print*, 'NULL ON SURFACE'
      inull = 1
    endif


    !check for null on edge with spine/separator directed along that edge
    if (edge > 0) then
      !print*, 'NULL ON EDGE WITH SPINE ALONG CELL EDGE  :)'
      inull = 1
    endif

  end

  logical function zeropoint(a,b)
    !Checks if there is the possibility that there is a zero on the faces a and b. If not, then there cannot be a crossing of zeroes on these faces.
    implicit none
    double precision, dimension(2,2) :: a, b
    integer :: totala, totalb
    integer :: i,j

    totala = 0
    totalb = 0

    do j = 1, 2
      do i = 1, 2
        if (a(i,j) > zero) then
          totala = totala + 1
        else if (a(i,j) < -zero) then
          totala = totala - 1
        else
          totala = totala + 100
        endif

        if (b(i,j) > zero) then
          totalb = totalb + 1
        else if (b(i,j) < -zero) then
          totalb = totalb - 1
        else
          totalb = totalb + 100
        endif
      enddo
    enddo

    if (abs(totala) == 4 .and. abs(totalb) == 4) then
      zeropoint = .false.
    else
      zeropoint = .true.
    endif

  end


  function blankface(a,b,c)
    !determines if the faces a, b and c are all zero. If so, return .false.
    double precision, dimension(2,2) :: a, b, c
    double precision :: atot, btot,ctot
    logical :: blankface

    atot=sum(abs(a))
    btot=sum(abs(b))
    ctot=sum(abs(c))

    blankface = .true.

    if (atot + btot + ctot <= zero*12.d0) then
      !print*, 'BLANK FACE'
      blankface =.false.
    endif

  end


  logical function check(x,y)
    !checks of x and y are between 0 and 1
    double precision :: x, y

    check = (x >= 0.) .and. (x <= 1.) .and. (y >= 0.) .and. (y <= 1.)
  end




  subroutine face_solve(facex,facey,facez,cross,sign)
    !tests if the zeroes on facex and facey cross, and if so, what the sign of facez is at the crossing point.
    !Uses the Method described in Haynes et al. 2007

    implicit none
    double precision, dimension(2,2) :: facex,facey,facez
    integer :: cross, sign
    double precision :: a1, b1, c1, d1 !bilinear coefficients (facex)
    double precision :: a2, b2, c2, d2 !bilinear coefficients (facey)
    double precision :: a,b,c !quadratic coefficients ax^2+bx+c=0

    double precision :: x1,x2,y1,y2
    double precision :: det

    double precision :: zcomp

    ! if there could be an intersection of the zeroes on facex and face y, or if the faces are not completely zero...
    if (zeropoint(facex,facey) .or. blankface(facex,facey,facez)) then
      !get bilinear coefficients
      a1 = facex(1,1)
      b1 = facex(2,1) - facex(1,1)
      c1 = facex(1,2) - facex(1,1)
      d1 = facex(2,2) - facex(2,1) - facex(1,2) + facex(1,1)

      a2 = facey(1,1)
      b2 = facey(2,1) - facey(1,1)
      c2 = facey(1,2) - facey(1,1)
      d2 = facey(2,2) - facey(2,1) - facey(1,2) + facey(1,1)

      !get quadratic coefficients (ax^2 + bx + c = 0)
      a = b1*d2 - b2*d1
      b = (a1*d2 - a2*d1) + (b1*c2 - c1*b2)
      c = a1*c2 - a2*c1

    !determinant of quadratic
      det = b*b - 4.d0*a*c

      if (det >= 0.) then !there is a solution
        if (abs(a) < zero) then !have to solve linear
          if (b == 0.) then !no solution exists
            x1 = -1.
            y1 = -1.
          else !solution exists
            x1 = -c/b
            y1 = -1.*(a1 + b1*x1)/(c1 + d1*x1)
            x2 = -1.
            y2= - 1.
          endif
        else !have to solve quadratic
          if (det == 0.) then !one solution
            x1 = -b/(2.*a)
            y1 = -1.*(a1 + b1*x1)/(c1 + d1*x1)
            x2 = -1.
            y2 = -1.
          else !two solutions
            x1 = -b/(2*a) + sqrt(det)/(2.*a)
            y1 = -1.*(a1 + b1*x1)/(c1 + d1*x1)

            x2 = -b/(2*a) - sqrt(det)/(2.*a)
            y2 = -1.*(a1 + b1*x2)/(c1 + d1*x2)
          endif
        endif

        sign=0
        cross=0

        if (check(x1, y1)) then !if x and y are on face (between 0 and 1)
          cross = cross + 1

          zcomp = bilinear_cell(x1, y1, facez) !'z' component of field at crossing point
          if (zcomp > zero) then
            sign = sign + 1
          else if (zcomp < -zero) then
            sign = sign - 1
          else
            sign = sign + 100
          endif
        endif

        if (check(x2, y2)) then
          cross = cross + 1

          zcomp = bilinear_cell(x2, y2, facez)
          if (zcomp > zero) then
            sign = sign + 1
          else if (zcomp < -zero) then
            sign = sign - 1
          else
            sign = sign + 100
          endif
        endif
      else !b^2-4ac less than zero - no intersection
        sign = 0
        cross = 0
      endif
    else !no intersection or faces have no field on them
      cross = 0
      sign = 0
    endif

  end subroutine



  integer function edge_check(facex,facey,facez)
    !check for null along edges of a cell face
    implicit none
    double precision, dimension(2,2) :: facex,facey,facez
    double precision, dimension(2) :: linex, liney , linez
    !double precision :: x,y,z

    edge_check = 0

    linex = facex(:,1)
    liney = facey(:,1)
    linez = facez(:,1)

    edge_check = edge_check + line_check(linex, liney, linez) !check for null along this edge

    linex = facex(:,2)
    liney = facey(:,2)
    linez = facez(:,2)

    edge_check = edge_check + line_check(linex, liney, linez)

    linex = facex(1,:)
    liney = facey(1,:)
    linez = facez(1,:)

    edge_check = edge_check + line_check(linex, liney, linez)

    linex = facex(2,:)
    liney = facey(2,:)
    linez = facez(2,:)

    edge_check = edge_check + line_check(linex, liney, linez)

    if (edge_check > 0) edge_check = 1 !if nulls are found, set the counter to 1

  end


  integer function line_check(linex,liney,linez)
    !checks if null lies on an edge
    implicit none
    double precision, dimension(2) :: linex, liney, linez
    !double precision :: x, y, z
    double precision :: x1, x2, y1, y2, z1, z2

    x1 = linex(1)
    x2 = linex(2)
    y1 = liney(1)
    y2 = liney(2)
    z1 = linez(1)
    z2 = linez(2)

    if ((x1 == 0. .and. x2 == 0.) .and. (y1 == 0. .and. y2 == 0.)) then
      if (z1*z2 <= 0.) line_check = 1
    else if ((x1 == 0. .and. x2 == 0.) .and. (z1 == 0. .and. z2 == 0.)) then
      if (y1*y2 <= 0.) line_check = 1
    else if ((z1 == 0. .and. z2 == 0.) .and. (y1 == 0. .and. y2 == 0.)) then
      if (x1*x2 <= 0.) line_check = 1
    else
      line_check = 0
    endif

  end


  double precision function linear(x,line)
    !linear interpolation of two values
    double precision :: line(2), x

    linear = (1-x)*line(1) + x*line(2)

  end


  double precision function bilinear_cell(x,y,square)
    implicit none
    !interpolate the value of a function at (x,y) from the 4 corner values
    double precision :: x, y, square(2,2)
    double precision :: y1, y2

    if (.not. check(x, y)) print*, 'bilinear error!'


    if (x > 1.) x = 1.
    if (x < 0.) x = 0.
    if (y > 1.) y = 1.
    if (y < 0.) y = 0.

    y1 = (1-x)*square(1,1) + x*square(2,1)
    y2 = (1-x)*square(1,2) + x*square(2,2)

    bilinear_cell = y1*(1-y) + y2*y
  end


  double precision function trilinear_cell(x, y ,z, cube)
    implicit none
    !find the value of a function at (x,y,z) in a cell using the 8 vertices
    double precision :: cube(2,2,2)
    double precision :: x, y, z
    double precision :: f11, f12, f21, f22
    double precision :: f1, f2

    if (x > 1.) x = 1.
    if (x < 0.) x = 0.
    if (y > 1.) y = 1.
    if (y < 0.) y = 0.
    if (z > 1.) z = 1.
    if (z < 0.) z = 0.

    f11 = (1-x)*cube(1,1,1) + x*cube(2,1,1)
    f12 = (1-x)*cube(1,1,2) + x*cube(2,1,2)
    f21 = (1-x)*cube(1,2,1) + x*cube(2,2,1)
    f22 = (1-x)*cube(1,2,2) + x*cube(2,2,2)

    f1 = (1-y)*f11 + y*f21
    f2 = (1-y)*f12 + y*f22

    trilinear_cell = (1-z)*f1 + z*f2
  end

  subroutine newton_raphson(cubex,cubey,cubez,x,y,z,ierror)
    !attempts to find null within cell using Newton-Raphson
    implicit none
    double precision, dimension(2,2,2) :: cubex, cubey, cubez
    double precision :: x, y, z
    double precision, parameter :: dx=1.e-5
    double precision :: bx, by, bz
    double precision :: divbx, divby, divbz
    double precision :: xold, yold, zold
    integer :: ierror
    double precision :: b,diff
    integer :: icounter

    ! print*, 'In Newton',x,y,z
    ierror=0
    !x=0.5
    !y=0.5
    !z=0.5

    diff = 1
    b=10.

    icounter=0

    do while (b > 1.e-4 .or. diff < 1.e-4)
      icounter = icounter + 1
      xold = x
      yold = y
      zold = z

      bx = trilinear_cell(x,y,z,cubex)
      by = trilinear_cell(x,y,z,cubey)
      bz = trilinear_cell(x,y,z,cubez)

      b = sqrt(bx*bx + by*by + bz*bz)

      divbx = (trilinear_cell(x + dx/2.,y,z, cubex) - trilinear_cell(x - dx/2.,y,z, cubex))/dx
      divby = (trilinear_cell(x,y + dx/2.,z, cubey) - trilinear_cell(x,y - dx/2.,z, cubey))/dx
      divbz = (trilinear_cell(x,y,z + dx/2., cubez) - trilinear_cell(x,y,z - dx/2., cubez))/dx

      x = xold - bx/divbx
      y = yold - by/divby
      z = zold - bz/divbz

      diff = sqrt((x-xold)**2 + (y-yold)**2 + (z-zold)**2)

      if (icounter > 100000) then
        print *, 'Newton-Raphson did not converge!'
        ierror=1
        exit
      endif
    enddo
  end


  integer function switch(cube)
    !determines whether all the vertices of the cube have the same sign or not
    implicit none
    integer :: i,j,k,icount
    double precision :: cube(2,2,2)
    integer ::itest

    icount=0

    do k = 1, 2
      do j = 1, 2
        do i = 1, 2
          if (cube(i,j,k) > zero) then
            icount = icount + 1
          else if (cube(i,j,k) < -zero) then
            icount = icount - 1
          else
            icount = icount + 0
            !print *, 'ZERO'
          endif
        enddo
      enddo
    enddo

    if (abs(icount) < 8) then !don't all have same sign
      itest=1
    else !all have same sign
      itest=0
    endif

    switch = itest

  end function


  subroutine subgrid(bx,by,bz,dx,x,y,z,ierror)
    !Finds null within gridcell by subdividing it into smaller cells, and using the bilinear method on these cells. This is repeated until the required accuracy is achieved
    implicit none
    double precision, dimension(2,2,2) :: bx, by, bz
    double precision :: x, y, z
    double precision :: dx
    integer :: ierror

    integer :: i, j, k
    integer :: itestx,itesty,itestz, itest

    double precision, dimension(2,2,2) :: cbx,cby,cbz
    double precision, dimension(11,11,11) ::  cubex, cubey,cubez
    double precision :: xs,ys,zs

    !x=0.
    !y=0.
    !z=0.

    xs = x
    ys = y
    zs = z

    !print*,'dx=',dx
    !print*,'in subgrid',x,y,z

    ierror = 0
    !print*,x,y,z,'subgrid start'

    do k = 1, 11
      z = zs + (k-1)*dx
      do j = 1, 11
        y = ys + (j-1)*dx
        do i = 1, 11
          x = xs + (i-1)*dx
          cubex(i,j,k) = trilinear_cell(x, y, z, bx)
          cubey(i,j,k) = trilinear_cell(x, y, z, by)
          cubez(i,j,k) = trilinear_cell(x, y, z, bz)
          !if (k == 1 .and. j == 1) print*,x,y,z
        enddo
      enddo
    enddo

    x = xs
    y = ys
    z = zs

    outer: do k = 1, 10
      do j = 1, 10
        do i = 1, 10
          cbx = cubex(i:i+1, j:j+1, k:k+1)
          cby = cubey(i:i+1, j:j+1, k:k+1)
          cbz = cubez(i:i+1, j:j+1, k:k+1)

          itestx = switch(cbx)
          if (itestx == 0) cycle

          itesty = switch(cby)
          if (itesty == 0) cycle

          itestz = switch(cbz)
          if (itestz == 0) cycle

          call bilin_test(cbx,cby,cbz,itest) !check for null within (or on face/corner/edge of) cell

          if (itest == 1) then
            ! print*, 'Null!',xs+dx*(i-1),ys+dx*(j-1),zs+dx*(k-1)
            x = xs + dx*(i-1)
            y = ys + dx*(j-1)
            z = zs + dx*(k-1)
            exit outer
          endif
        enddo
      enddo

      if (k == 10 .and. j == 11 .and. i == 11) then
        print*,'NO NULL FOUND - ERROR :('
        ierror=1
      endif
      !print*,k,j,i

    enddo outer
  end




  subroutine brute_force(bx,by,bz,x,y,z)
  !finds a null in a cell by brite force. Namely it subdivides the cell into 10x10x10 subcells, and finds the subcell with the minimum value of |B|. This is continued until the required accuracy (sig_figs) is achieved.
    implicit none
    double precision, dimension(2,2,2) :: bx, by,bz
    double precision :: x, y, z
    integer :: i,j,k

    double precision :: dx, xs, ys,zs
    double precision :: minx, miny,minz,min
    double precision :: bbx,bby,bbz
    double precision :: b
    double precision :: range

    xs = 0.
    ys = 0.
    zs = 0.

    min = 100.

    dx = 1.
    range = 5.

    do while (range > 10.**(-sig_figs)) !num=1,sig_figs
      range = range/5.
      dx = range/10
      do k = 0, 10
      !print*,k
      z = zs + k*dx
        do j = 0, 10
          y = ys + j*dx
          do i = 0, 10
            x = xs + i*dx
            bbx = trilinear_cell(x, y, z, bx)
            bby = trilinear_cell(x, y, z, by)
            bbz = trilinear_cell(x, y, z, bz)
            b = bbx*bbx + bby*bby + bbz*bbz
            if (b < min) then
              min = b
              minx = x
              miny = y
              minz = z
            endif
          enddo
        enddo
      enddo

      !range=range/5.

      xs = minx - range/10.
      ys = miny - range/10.
      zs = minz - range/10.
    enddo
    !print*,minx,miny,minz,min
    x = minx
    y = miny
    z = minz

  end subroutine


  subroutine add_element(x,element)
    !adds an element to the end of an array
    implicit none
    double precision, allocatable :: x(:), dummy(:)
    double precision :: element
    integer :: isize

    isize = size(x)

    allocate(dummy(isize + 1))

    dummy(1:isize) = x
    dummy(isize+1) = element

    deallocate(x)
    allocate(x(isize + 1))

    x = dummy

  end


  subroutine remove_element(x,index)
    !removes element number index from an array
    double precision, allocatable :: x(:), dummy(:)
    integer :: index
    integer :: isize

    isize = size(x)

    allocate(dummy(isize-1))

    dummy(1:index-1) = x(1:index-1)
    dummy(index:isize-1) = x(index+1:isize)

    deallocate(x)
    allocate(x(isize-1))

    x = dummy

    deallocate(dummy)

  end


  subroutine normalise(a,b,c)
    double precision, dimension(2,2,2) :: a, b, c
    double precision :: scale

    scale = max(maxval(abs(a)), maxval(abs(b)), maxval(abs(c)))

    a = a/scale
    b = b/scale
    c = c/scale

  end


  double precision function div(bx,by,bz)
    !calculates div b in a cell (not used)
    double precision, dimension (2,2,2) :: bx, by,bz
    double precision, dimension(2) :: bbx, bby, bbz

    bbx(:) = 0.25*(bx(:,1,1) + bx(:,1,2) + bx(:,2,1) + bx(:,2,2))
    bby(:) = 0.25*(by(1,:,1) + by(1,:,2) + by(2,:,1) + by(2,:,2))
    bbz(:) = 0.25*(bz(1,1,:) + bz(1,2,:) + bz(2,1,:) + bx(2,2,:))

    div = bbx(2) - bbx(1) + bby(2) - bby(1) + bbz(2) - bbz(1)
  end


  subroutine remove_duplicates(x,y,z,n,xp,yp,zp)
    !removes duplicate nulls (nulls closer than 0.01 gridcells)
    implicit none
    double precision, allocatable :: x(:), y(:), z(:)
    double precision, allocatable :: xp(:), yp(:), zp(:)
    double precision, allocatable :: distances(:,:)
    integer :: n
    integer :: i,j
    integer :: counter
    double precision :: sep
    integer :: exitcondition

    if (n > 1) then
      exitcondition=0

      outer: do while (exitcondition == 0)
        !allocate(distances(n,n))
        do j = 1, n
          do i = j+1, n
            sep = (x(i) - x(j))**2 + (y(i) - y(j))**2 + (z(i) - z(j))**2
            if (sep < 1d-3) then
              print*, 'removing duplcate at index', j
              call remove_element(x, j)
              call remove_element(y, j)
              call remove_element(z, j)
              call remove_element(xp, j)
              call remove_element(yp, j)
              call remove_element(zp, j)
              n = n - 1
              cycle outer
            endif
          enddo
        enddo
        exitcondition = 1
      enddo outer

      allocate(distances(n,n))

      distances = 1000000.
      do j = 1, n
        do i = j+1, n
          !print*,nnulls,i,j
          distances(i,j) = (x(i) - x(j))**2 + (y(i) - y(j))**2 + (z(i) - z(j))**2
        enddo
      enddo

      print*,''
      print*,'Now displaying all nulls with less than 1 gridcell spacing'
      print*,'between their neighbours:'
      distances = sqrt(distances)
      counter = 0
      do j = 1, n
        do i = j+1, n
          if (distances(i,j) < 1.0) then
            print *, i, j, distances(i,j)
            counter = counter + 1
          endif
        enddo
      enddo

      print*, 'Number of close-by nulls = ', counter
      deallocate(distances)
    endif
  end subroutine

end module
