! module for nullfinder
module nf_mod
  
  use params

  implicit none

  contains

  subroutine bilin_test(cbx,cby,cbz,inull)
    ! tests surface of cube using bilinear method
    implicit none

    integer(int32) :: i, inull, n

    integer(int32) :: cross, sign
    integer(int32) :: edge

    real(np), dimension(2,2,2) :: cbx,cby,cbz

    real(np), dimension(2,2) :: facex, facey, facez

    integer(int32), dimension(6,6) :: test
    integer(int32), dimension(6) :: test2
    integer(int32), dimension(3) :: test3
    
    ! test is an array where the first index refers to the face number (1:6) and the second index is follows:
    ! 1- x-y intersection (1/0)
    ! 2- sign of z at this intersection (-1/1 /100 if null on face)
    ! 3- xz intersection (1/0)
    ! 4- sign of y on this intersection (-1/1 /100 if null on face)
    ! 5- yz intersection (1/0)
    ! 6- sign on x on this inersection (-1/1 /100 if null on face)

    inull = 0 ! variable defining null. 1=null, 0=no null
    test = 0 ! set test array to zero
    edge = 0 ! variable defining null on edge. 1=null, 0=no null

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

    ! face #2 : left face
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

    ! face #3 : right face
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

    ! face 4: front face
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

#if debug
    print*, 'Face 4: front'
    do i = 1, 2
      print*, 'x-coord', facex(1:2,i)
    enddo
    do i = 1, 2
      print*, 'z-coord', facez(1:2,i)
    enddo
#endif 

    edge = edge + edge_check(facex, facey, facez)

    ! face 5: rear face
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

#if debug
    print*, 'Face 5: rear'
    do i = 1, 2
      print*, 'x-coord', facex(1:2,i)
    enddo
    do i = 1, 2
      print*, 'z-coord', facez(1:2,i)
    enddo
#endif  

    edge = edge + edge_check(facex, facey, facez)

    ! face 6: top face
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

    ! check for nulls within cell

    ! total of each column in 'test'
    test2 = sum(test,1)

    test3 = 0
    do i = 1, 3
      if (test2(2*i-1) >= 0) then ! if there is a crossing
        if (abs(test2(2*i)) < test2(2*i-1)) test3(i) = 1
      endif
    enddo

#if debug
    ! print out test array for checking
    do i = 1, 6
      print*, test(i,:)
    enddo
    print*, 'test2'
    print*, test2
    print*, 'test3'
    print*, test3
    print*, 'edge check'
    print*, edge
    if (sum(test3) == 3) print*, 'NULL IN CELL'
    if (sum(test(:,2)) > 10 .or. sum(test(:,4)) > 10 .or. sum(test(:,6)) > 10) print*, 'NULL ON SURFACE'
    if (edge > 0) print*, 'NULL ON EDGE WITH SPINE ALONG CELL EDGE :)'
    ! stop
#endif

    ! check for null in the cell
    if (sum(test3) == 3) inull = 1

    ! check for null on surface of cell
    if (sum(test(:,2)) > 10 .or. sum(test(:,4)) > 10 .or. sum(test(:,6)) > 10) inull = 1


    ! check for null on edge with spine/separator directed along that edge
    if (edge > 0) inull = 1

  end

  !********************************************************************************

  logical function zeropoint(a,b)
    ! Checks if there is the possibility that there is a zero on the faces a and b.
    ! If not, then there cannot be a crossing of zeroes on these faces.
    implicit none
    real(np), dimension(2,2) :: a, b
    integer(int32) :: totala, totalb
    integer(int32) :: i, j

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

  !********************************************************************************

  function blankface(a,b,c)
    ! determines if the faces a, b and c are all zero. If so, return .false.
    real(np), dimension(2,2) :: a, b, c
    real(np) :: atot, btot,ctot
    logical :: blankface

    atot = sum(abs(a))
    btot = sum(abs(b))
    ctot = sum(abs(c))

    blankface = .true.

    if (atot + btot + ctot <= zero*12.0_np) then
      !print*, 'BLANK FACE'
      blankface = .false.
    endif

  end

  !********************************************************************************

  logical function check(x,y)
    !checks of x and y are between 0 and 1
    real(np) :: x, y

    check = (x >= 0.0_np) .and. (x <= 1.0_np) .and. (y >= 0.0_np) .and. (y <= 1.0_np)
  end

  !********************************************************************************

  subroutine face_solve(facex,facey,facez,cross,sign)
    !tests if the zeroes on facex and facey cross, and if so, what the sign of facez is at the crossing point.
    !Uses the Method described in Haynes et al. 2007

    implicit none
    real(np), dimension(2,2) :: facex, facey, facez
    integer(int32) :: cross, sign, i, nsol
    real(np) :: a1, b1, c1, d1 !bilinear coefficients (facex)
    real(np) :: a2, b2, c2, d2 !bilinear coefficients (facey)
    real(np) :: a, b, c !quadratic coefficients ax^2+bx+c=0

    real(np), dimension(2,2) :: x, y
    real(np) :: det

    real(np) :: zcomp

    x = -1
    y = -1

    sign = 0
    cross = 0
    nsol = 0

#if debug
    print*,'zeropoint = ', zeropoint(facex,facey)
    print*,'blankface = ', blankface(facex,facey,facez)
#endif

    ! if there could be an intersection of the zeroes on facex and face y, or if the faces are not completely zero...
    if (zeropoint(facex,facey) .or. blankface(facex,facey,facez)) then
      ! get bilinear coefficients
      a1 = facex(1,1)
      b1 = facex(2,1) - facex(1,1)
      c1 = facex(1,2) - facex(1,1)
      d1 = facex(2,2) - facex(2,1) - facex(1,2) + facex(1,1)

      a2 = facey(1,1)
      b2 = facey(2,1) - facey(1,1)
      c2 = facey(1,2) - facey(1,1)
      d2 = facey(2,2) - facey(2,1) - facey(1,2) + facey(1,1)

      ! get quadratic coefficients (ax^2 + bx + c = 0)
      a = b1*d2 - b2*d1
      b = (a1*d2 - a2*d1) + (b1*c2 - c1*b2)
      c = a1*c2 - a2*c1

      ! determinant of quadratic
      det = b**2 - 4*a*c

#if debug
    print*, 'a1, b1, c1, d1'
    print*, a1, b1, c1, d1
    print*, 'a2, b2, c2, d2'
    print*, a2, b2, c2, d2
    print*, 'ax^2 + bx + c = 0: a, b, c'
    print*, a, b, c
    print*, 'det (b**2 - 4*a*c)'
    print*, det
    print*, 'zero', zero
#endif

      if (det >= 0) then ! there is a solution
        if (abs(a) < zero) then !have to solve linear
          if (b /= 0.0_np) then
            x(1,:) = -c/b ! solution exists
            ! x(1,2) = -c/b ! solution exists
            nsol = 1
          endif
        else ! have to solve quadratic
          if (det == 0.0_np) then ! one solution
            x(1,:) = -b/(2*a)
            ! x(1,2) = -b/(2*a)
            nsol = 1
          else ! two solutions
            x(1,:) = (-b + sqrt(det))/(2*a)
            ! x(1,2) = (-b + sqrt(det))/(2*a)
            x(2,:) = (-b - sqrt(det))/(2*a)
            ! x(2,2) = (-b - sqrt(det))/(2*a)
            nsol = 2
          endif
        endif
      endif
      if (nsol > 0) then
        y(1,1) = -(a1 + b1*x(1,1))/(c1 + d1*x(1,1))
        y(1,2) = -(a2 + b2*x(1,2))/(c2 + d2*x(1,2))
        if (nsol == 2) then
          y(2,1) = -(a1 + b1*x(2,1))/(c1 + d1*x(2,1))
          y(2,2) = -(a2 + b2*x(2,2))/(c2 + d2*x(2,2))
        endif
      else
        ! if the x equation isn't solvable then try the y equation
        ! get quadratic coefficients (ay^2 + by + c = 0)
        a = c1*d2 - c2*d1
        b = (a1*d2 - a2*d1) - (b1*c2 - c1*b2)
        c = a1*b2 - a2*b1

        det = b**2 - 4*a*c

        if (det >= 0) then !there is a solution
          if (abs(a) < zero) then !have to solve linear
            if (b /= 0.0_np) then
              y(1,:) = -c/b !solution exists
              ! y(1,2) = -c/b !solution exists
              nsol = 1
            endif
          else !have to solve quadratic
            if (det == 0.0_np) then !one solution
              y(1,:) = -b/(2*a)
              ! y(1,2) = -b/(2*a)
              nsol = 1
            else !two solutions
              y(1,:) = (-b + sqrt(det))/(2*a)
              ! y(1,2) = (-b + sqrt(det))/(2*a)
              y(2,:) = (-b - sqrt(det))/(2*a)
              ! y(2,2) = (-b - sqrt(det))/(2*a)
              nsol = 2
            endif
          endif
          if (nsol > 0) then
            x(1,1) = -(a1 + c1*y(1,1))/(b1 + d1*y(1,1))
            x(1,2) = -(a2 + c2*y(1,2))/(b2 + d2*y(1,2))
            if (nsol == 2) then
              x(2,1) = -(a1 + c1*y(2,1))/(b1 + d1*y(2,1))
              x(2,2) = -(a2 + c2*y(2,2))/(b2 + d2*y(2,2))
            endif
          endif
        endif
        ! if (c1 == 0.0_np .and. c2 == 0.0_np .and. d1 == 0.0_np .and. d2 == 0.0_np .and. a /= 0.0_np .and. b /= 0.0_np .and. c /= 0.0_np) then
        ! print*, '--------------------------------------------------------------'
        !   do i = 1, 2
        !     print*, facex(:,i)
        !   enddo
        !   do i = 1, 2
        !     print*, facey(:,i)
        !   enddo

        !   print*, ''
        !   print*, det, a, b, c
        !   print*, ''
        !   print*, a1, b1, c1, d1
        !   print*, a2, b2, c2, d2
        !   print*, 'x1                        ', 'y1                          ', 'x2                         ', 'y2        '
        !   print*, x1, y1, x2, y2
        !   print*, 'c2 + d2*x1               ', 'c2 + d2*x2                 ', 'c1 + d1*x1                   ', 'c1 + d1*x2    '
        !   print*, c2 + d2*x1, c2 + d2*x2, c1 + d1*x1, c1 + d1*x2
        ! endif
      endif

      if (check(x(1,1), y(1,1))) then !if x and y are on face (between 0 and 1)
        cross = cross + 1

        zcomp = bilinear_cell(x(1,1), y(1,1), facez) !'z' component of field at crossing point
        if (zcomp > zero) then
          sign = sign + 1
        else if (zcomp < -zero) then
          sign = sign - 1
        else
          sign = sign + 100
        endif
      endif

      if (check(x(1,2), y(1,2))) then
        cross = cross + 1

        zcomp = bilinear_cell(x(1,2), y(1,2), facez)
        if (zcomp > zero) then
          sign = sign + 1
        else if (zcomp < -zero) then
          sign = sign - 1
        else
          sign = sign + 100
        endif
      endif

      if (check(x(2,1), y(2,1))) then
        cross = cross + 1

        zcomp = bilinear_cell(x(2,1), y(2,1), facez)
        if (zcomp > zero) then
          sign = sign + 1
        else if (zcomp < -zero) then
          sign = sign - 1
        else
          sign = sign + 100
        endif
      endif
    
      if (check(x(2,2), y(2,2))) then
        cross = cross + 1

        zcomp = bilinear_cell(x(2,2), y(2,2), facez)
        if (zcomp > zero) then
          sign = sign + 1
        else if (zcomp < -zero) then
          sign = sign - 1
        else
          sign = sign + 100
        endif
      endif

    endif

  end subroutine

  !********************************************************************************

  function edge_check(facex,facey,facez)
    !check for null along edges of a cell face
    implicit none
    real(np), dimension(2,2) :: facex,facey,facez
    real(np), dimension(2) :: linex, liney , linez
    integer(int32) :: edge_check

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

  !********************************************************************************

  function line_check(linex,liney,linez)
    !checks if null lies on an edge
    implicit none
    real(np), dimension(2) :: linex, liney, linez
    real(np) :: x1, x2, y1, y2, z1, z2
    integer(int32) :: line_check

    x1 = linex(1)
    x2 = linex(2)
    y1 = liney(1)
    y2 = liney(2)
    z1 = linez(1)
    z2 = linez(2)

    if ((x1 == 0.0_np .and. x2 == 0.0_np) .and. (y1 == 0.0_np .and. y2 == 0.0_np)) then
      if (z1*z2 <= 0.0_np) line_check = 1
    else if ((x1 == 0.0_np .and. x2 == 0.0_np) .and. (z1 == 0.0_np .and. z2 == 0.0_np)) then
      if (y1*y2 <= 0.0_np) line_check = 1
    else if ((z1 == 0.0_np .and. z2 == 0.0_np) .and. (y1 == 0.0_np .and. y2 == 0.0_np)) then
      if (x1*x2 <= 0.0_np) line_check = 1
    else
      line_check = 0
    endif

  end

  !********************************************************************************

  real(np) function linear(x,line)
    !linear interpolation of two values
    real(np) :: line(2), x

    linear = (1-x)*line(1) + x*line(2)

  end

  !********************************************************************************

  real(np) function bilinear_cell(x,y,square)
    implicit none
    !interpolate the value of a function at (x,y) from the 4 corner values
    real(np) :: x, y, square(2,2)
    real(np) :: y1, y2

    if (.not. check(x, y)) print*, 'bilinear error!'

    if (x > 1.0_np) x = 1
    if (x < 0.0_np) x = 0
    if (y > 1.0_np) y = 1
    if (y < 0.0_np) y = 0

    y1 = (1-x)*square(1,1) + x*square(2,1)
    y2 = (1-x)*square(1,2) + x*square(2,2)

    bilinear_cell = y1*(1-y) + y2*y
  end

  !********************************************************************************

  real(np) function trilinear_cell(x, y ,z, cube)
    implicit none
    !find the value of a function at (x,y,z) in a cell using the 8 vertices
    real(np) :: cube(2,2,2)
    real(np) :: x, y, z
    real(np) :: f11, f12, f21, f22
    real(np) :: f1, f2

    if (x > 1.0_np) x = 1
    if (x < 0.0_np) x = 0
    if (y > 1.0_np) y = 1
    if (y < 0.0_np) y = 0
    if (z > 1.0_np) z = 1
    if (z < 0.0_np) z = 0

    f11 = (1-x)*cube(1,1,1) + x*cube(2,1,1)
    f12 = (1-x)*cube(1,1,2) + x*cube(2,1,2)
    f21 = (1-x)*cube(1,2,1) + x*cube(2,2,1)
    f22 = (1-x)*cube(1,2,2) + x*cube(2,2,2)

    f1 = (1-y)*f11 + y*f21
    f2 = (1-y)*f12 + y*f22

    trilinear_cell = (1-z)*f1 + z*f2
  end

  !********************************************************************************

  function switch(cube)
    !determines whether all the vertices of the cube have the same sign or not
    implicit none
    integer(int32) :: i, j, k, icount
    real(np) :: cube(2,2,2)
    integer(int32) :: itest
    integer(int32) :: switch

    icount = 0

    do k = 1, 2
      do j = 1, 2
        do i = 1, 2
          if (cube(i,j,k) > zero) then
            icount = icount + 1
          else if (cube(i,j,k) < -zero) then
            icount = icount - 1
          else
            icount = icount + 0
          endif
        enddo
      enddo
    enddo

    if (abs(icount) < 8) then ! don't all have same sign
      itest = 1
    else ! all have same sign
      itest = 0
    endif

    switch = itest

  end function

  !********************************************************************************

  subroutine subgrid(bx,by,bz,dx,x,y,z,ierror)
    !Finds null within gridcell by subdividing it into smaller cells, and using the bilinear method on these cells. This is repeated until the required accuracy is achieved
    implicit none
    real(np), dimension(2,2,2) :: bx, by, bz
    real(np) :: x, y, z
    real(np) :: dx
    integer(int32) :: ierror

    integer(int32) :: i, j, k
    integer(int32) :: itestx,itesty,itestz, itest

    real(np), dimension(2,2,2) :: cbx,cby,cbz
    real(np), dimension(11,11,11) ::  cubex, cubey,cubez
    real(np) :: xs,ys,zs

    xs = x
    ys = y
    zs = z

    ierror = 0

    do k = 1, 11
      z = zs + (k-1)*dx
      do j = 1, 11
        y = ys + (j-1)*dx
        do i = 1, 11
          x = xs + (i-1)*dx
          cubex(i,j,k) = trilinear_cell(x, y, z, bx)
          cubey(i,j,k) = trilinear_cell(x, y, z, by)
          cubez(i,j,k) = trilinear_cell(x, y, z, bz)
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

#if debug
          print*, 'Subgrid: calling bilinear test', i, j, k, dx
#endif

          call bilin_test(cbx,cby,cbz,itest) ! check for null within (or on face/corner/edge of) cell

          if (itest == 1) then
            x = xs + dx*(i-1)
            y = ys + dx*(j-1)
            z = zs + dx*(k-1)
            exit outer
          endif
        enddo
      enddo

      if (k == 10 .and. j == 11 .and. i == 11) then
        print*,'NO NULL FOUND - ERROR :('
        ierror = 1
      endif

    enddo outer
  end

  !********************************************************************************

  subroutine add_element(x,element)
    ! adds an element to the end of an array
    implicit none
    real(np), allocatable :: x(:), dummy(:)
    real(np) :: element
    integer(int32) :: isize

    isize = size(x)

    allocate(dummy(isize + 1))

    dummy(1:isize) = x
    dummy(isize+1) = element

    deallocate(x)
    allocate(x(isize + 1))

    x = dummy

  end

  !********************************************************************************

  subroutine remove_element(x,index)
    ! removes element number index from an array
    real(np), allocatable :: x(:), dummy(:)
    integer(int32) :: index
    integer(int32) :: isize

    isize = size(x)

    allocate(dummy(isize-1))

    dummy(1:index-1) = x(1:index-1)
    dummy(index:isize-1) = x(index+1:isize)

    deallocate(x)
    allocate(x(isize-1))

    x = dummy

    deallocate(dummy)

  end

  !********************************************************************************

  subroutine normalise(a,b,c)
    
    real(np), dimension(2,2,2) :: a, b, c
    real(np) :: scale

    scale = max(maxval(abs(a)), maxval(abs(b)), maxval(abs(c)))

    a = a/scale
    b = b/scale
    c = c/scale

  end

  !********************************************************************************

  subroutine remove_duplicates(x,y,z,n,xp,yp,zp)
    ! removes duplicate nulls (nulls closer than 0.01 gridcells)
    implicit none
    real(np), allocatable :: x(:), y(:), z(:)
    real(np), allocatable :: xp(:), yp(:), zp(:)
    real(np), allocatable :: distances(:,:)
    integer(int32) :: n
    integer(int32) :: i,j
    integer(int32) :: counter
    real(np) :: sep
    integer(int32) :: exitcondition

    if (n > 1) then
      exitcondition=0

      outer: do while (exitcondition == 0)
        !allocate(distances(n,n))
        do j = 1, n
          do i = j+1, n
            sep = (x(i) - x(j))**2 + (y(i) - y(j))**2 + (z(i) - z(j))**2
            if (sep < 1e-3_np) then
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

      distances = 1e6_np
      do j = 1, n
        do i = j+1, n
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
          if (distances(i,j) < 1) then
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