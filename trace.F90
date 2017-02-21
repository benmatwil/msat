module trace
  use params
  use common

  !rkf45 parameters
  real(np), parameter :: k21 = 0.25_np
  real(np), parameter :: k31 = 3.0_np/32.0_np, k32 = 9.0_np/32.0_np
  real(np), parameter :: k41 = 1932.0_np/2197.0_np, k42 = -7200.0_np/2197.0_np, k43 = 7296.0_np/2197.0_np
  real(np), parameter :: k51 = 439.0_np/216.0_np, k52 = -8.0_np, k53 = 3680.0_np/513.0_np, k54 = -845.0_np/4104.0_np
  real(np), parameter :: k61 = -8.0_np/27.0_np, k62 = 2.0_np, k63 = -3544.0_np/2565.0_np, &
    k64 = 1859.0_np/4104.0_np, k65 = -11.0_np/40.0_np

  real(np), parameter :: y1 = 25.0_np/216.0_np, y3 = 1408.0_np/2565.0_np, y4 = 2197.0_np/4101.0_np, y5 = -1.0_np/5.0_np
  real(np), parameter :: z1 = 16.0_np/135.0_np, z3 = 6656.0_np/12825.0_np, z4 = 28561.0_np/56430.0_np, &
    z5 = -9.0_np/50.0_np, z6 = 2.0_np/55.0_np

  integer :: terror

  contains

    subroutine trace_line(r,sign,h)
    ! traces a line from 'r' for 'nsteps' integration steps in the direction along the line as specified by 'sign'. Each step is of length h
      real(np), dimension(3) :: r, r0, r1
      integer :: sign
      real(np) :: h, hdum, stepdist

      stepdist = h
      hdum = 0.0_np
      r0 = r

      do while (hdum < stepdist .and. .not. outedge(r))
        h = sign*(stepdist-hdum)
        call rk45(r, h)
        ! call edgecheck(r0)
        hdum = hdum + abs(h)
      enddo
      
      ! if (modulus(r-r0) < 0.1_np*stepdist .and. .not. outedge(r)) then
      !   print*, 'field line not tracin', modulus(r-r0), stepdist, sign
      !   print*, xmin, xmax, ymin, ymax, zmin, zmax
      !   print*, r
      !   terror = 1
      ! endif

    end

    !********************************************************************************

    subroutine rk45(r,h)
    !runge-kutta fehlberg integrator. Calculates a 4th order estimate (y) and a
    !fifth order estimate (z) and calculates the difference between them. From this it
    !determines the optimal step length with which to integrate the function by.
      real(np) :: h, hvec(3), mindist, s
      real(np), dimension(3) :: k1, k2, k3, k4, k5, k6
      real(np), dimension(3) :: r, r0, rtest, y, z, dr
      
      r0 = r
      
      !minimum physical distance corresponding to fraction of gridcell (h)
      dr = getdr(r0)
      mindist = minval(dr)*h
      
      !vector containing the grid's h for each direction
      !(so h for each direction is the same physical length, equal to mindist)
      hvec = mindist/dr

      !get rk values k1--k6
      rtest = r0
      k1 = hvec*normalise(trilinear(rtest, bgrid))
      rtest = r0 + k21*k1
      ! call edgecheck(rtest)
      k2 = hvec*normalise(trilinear(rtest, bgrid))
      rtest = r0 + k31*k1 + k32*k2
      ! call edgecheck(rtest)
      k3 = hvec*normalise(trilinear(rtest, bgrid))
      rtest = r0 + k41*k1 + k42*k2 + k43*k3
      ! call edgecheck(rtest)
      k4 = hvec*normalise(trilinear(rtest, bgrid))
      rtest = r0 + k51*k1 + k52*k2 + k53*k3 + k54*k4
      ! call edgecheck(rtest)
      k5 = hvec*normalise(trilinear(rtest, bgrid))
      rtest = r0 + k51*k1 + k52*k2 + k53*k3 + k54*k4
      ! call edgecheck(rtest)
      k6 = hvec*normalise(trilinear(rtest, bgrid))
      

      !get 4th order (y) and 5th order (z) estimates
      y = y1*k1 + y3*k3 + y4*k4 + y5*k5
      z = z1*k1 + z3*k3 + z4*k4 + z5*k5 + z6*k6

      !calculate optimum step length (s = hoptimum/h)
      s = 0.84_np*(tol/modulus(z-y))**0.25_np
      
      if (abs(s*h) < stepmin) s = stepmin/abs(h)
      if (s > 1) s = 1.0_np

      rtest = r0
      k1 = s*hvec*normalise(trilinear(rtest, bgrid))
      rtest = r0 + k21*k1
      ! call edgecheck(rtest)
      k2 = s*hvec*normalise(trilinear(rtest, bgrid))
      rtest = r0 + k31*k1 + k32*k2
      ! call edgecheck(rtest)
      k3 = s*hvec*normalise(trilinear(rtest, bgrid))
      rtest = r0 + k41*k1 + k42*k2 + k43*k3
      ! call edgecheck(rtest)
      k4 = s*hvec*normalise(trilinear(rtest, bgrid))
      rtest = r0 + k51*k1 + k52*k2 + k53*k3 + k54*k4
      ! call edgecheck(rtest)
      k5 = s*hvec*normalise(trilinear(rtest, bgrid))

      r = r0 + y1*k1 + y3*k3 + y4*k4 + y5*k5

      h = h*s

      ! if (modulus(r-r0) < 0.1*h .and. .not. outedge(r)) then
      !   print*, 'trace failure',modulus(r-r0),h
      !   stop
      ! endif

    end subroutine

    !********************************************************************************

    function getdr(r)
      !outputs length of one gridcell in 'physical' length units (essentially (dx,dy,dz))
      real(np) :: r(3), rcheck(3) !grid cell number
      real(np) :: dx, dy, dz
      real(np) :: getdr(3)
      integer :: i, j, k
      real(np) :: xc, yc, zc !x, y and z at the midpoint of the cell

      rcheck = r
      call edgecheck(rcheck)

      i = floor(rcheck(1))
      j = floor(rcheck(2))
      k = floor(rcheck(3))
      
      dx = x(i+1)-x(i)
      dy = y(j+1)-y(j)
      dz = z(k+1)-z(k)
      
      xc = x(i) + dx/2
      yc = y(j) + dy/2
      zc = z(k) + dz/2

#if cartesian
      !cartesian coordinates
      getdr = [dx, dy, dz]
#elif spherical
      !spherical coordinates
      getdr = [dx, xc*dy, xc*sin(yc)*dz]
#elif cylindrical
      !cylindrical coordinates
      getdr = [dx, xc*dy, dz]
#else
      stop 'No coordinate system selected. Please select a valid one when compiling'
#endif
      
    end function

end module
