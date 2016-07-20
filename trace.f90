module trace
use params
use common

!rkf45 parameters
double precision, parameter :: k21 = 0.25
double precision, parameter :: k31 = 3./32., k32 = 9./32.
double precision, parameter :: k41 = 1932./2197., k42 = -7200./2197., k43 = 7296./2197.
double precision, parameter :: k51 = 439./216., k52 = -8., k53 = 3680./513., k54 = -845./4104.
double precision, parameter :: k61 = -8./27., k62 = 2., k63 = -3544./2565., k64 = 1859./4104., k65 = -11./40.

double precision, parameter :: y1 = 25./216., y3 = 1408./2565., y4 = 2197./4101., y5 = -1./5.
double precision, parameter :: z1 = 16./135., z3 = 6656./12825., z4 = 28561./56430.,z5 = -9./50.,z6 = 2./55.

contains

  subroutine trace_line(r,sign,h)
  !traces a line from 'r' for 'nsteps' integration steps in the direction along the line as specified by 'sign'. Each step is of length h
    double precision :: r(3), r0(3)
    integer :: sign
    double precision :: h, hdum, stepdist

    stepdist = h
    hdum = 0
    r0 = r

    do while (hdum < stepdist .and. .not. outedge(r))
      h = sign*(stepdist-hdum)
      call rk45(r,h)
      hdum = hdum + abs(h)
    enddo

    if (modulus(r-r0) < 0.1*stepdist) then
      !print *,'field line not tracing',modulus(r-r0),stepdist
      ierror = 1
    endif

  end

  !********************************************************************************

  subroutine rk45(r,h)
  !runge-kutta fehlberg integrator. Calculates a 4th order estimate (y) and a
  !fifth order estimate (z) and calculates the difference between them. From this it
  !determines the optimal step length with which to integrate the function by.
    double precision :: h, hvec(3), mindist, s
    double precision, dimension(3) :: k1, k2, k3, k4, k5, k6
    double precision, dimension(3) :: r, r0, rtest, y, z, dl
    
    !minimum physical distance corresponding to fraction of gridcell (h)
    !dl = sign(dble([1,1,1]),getdl(r))
    !print*, 'getdl', r
    dl = getdl(r)
    mindist = minval(dl)*h
    
    !vector containing the grid's h for each direction
    !(so h for each direction is the same physical length, equal to mindist)
    !hvec = h*dl
    hvec = mindist/dl
    
    r0 = r

    !get rk values k1--k6
    rtest = r0
    !print*, 'k1'
    call edgecheck(rtest)
    k1 = hvec*normalise(trilinear(rtest, bgrid))
    rtest = r0 + k21*k1
    call edgecheck(rtest)
    !print*, 'k2'
    k2 = hvec*normalise(trilinear(rtest, bgrid))
    rtest = r0 + k31*k1 + k32*k2
    call edgecheck(rtest)
    !print*, 'k3'
    k3 = hvec*normalise(trilinear(rtest, bgrid))
    rtest = r0 + k41*k1 + k42*k2 + k43*k3
    call edgecheck(rtest)
    !print*, 'k4'
    k4 = hvec*normalise(trilinear(rtest, bgrid))
    rtest = r0 + k51*k1 + k52*k2 + k53*k3 + k54*k4
    call edgecheck(rtest)
    !print*, 'k5'
    k5 = hvec*normalise(trilinear(rtest, bgrid))
    rtest = r0 + k51*k1 + k52*k2 + k53*k3 + k54*k4
    call edgecheck(rtest)
    !print*, 'k6'
    k6 = hvec*normalise(trilinear(rtest, bgrid))
    

    !get 4th order (y) and 5th order (z) estimates
    y = y1*k1 + y3*k3 + y4*k4 + y5*k5
    z = z1*k1 + z3*k3 + z4*k4 + z5*k5 + z6*k6

    !calculate optimum step length (s = hoptimum/h)
    s = 0.84*(tol/modulus(z-y))**0.25
    
    if (abs(s*h) < stepmin) s = stepmin/abs(h)
    if (s > 1) s = 1

    rtest = r0
    call edgecheck(rtest)
    k1 = s*hvec*normalise(trilinear(rtest, bgrid))
    rtest = r0 + k21*k1
    call edgecheck(rtest)
    k2 = s*hvec*normalise(trilinear(rtest, bgrid))
    rtest = r0 + k31*k1 + k32*k2
    call edgecheck(rtest)
    k3 = s*hvec*normalise(trilinear(rtest, bgrid))
    rtest = r0 + k41*k1 + k42*k2 + k43*k3
    call edgecheck(rtest)
    k4 = s*hvec*normalise(trilinear(rtest, bgrid))
    rtest = r0 + k51*k1 + k52*k2 + k53*k3 + k54*k4
    call edgecheck(rtest)
    k5 = s*hvec*normalise(trilinear(rtest, bgrid))

    r = r0 + y1*k1 + y3*k3 + y4*k4 + y5*k5

    h = h*s

    if (modulus(r-r0) < 0.1*h) then
      print*, 'trace failure',modulus(r-r0),h
      !stop
    endif

  end subroutine

  !********************************************************************************

  function getdl(r)
    !outputs length of one gridcell in 'physical' length units (essentially (dx,dy,dz))
    double precision :: r(3), rcheck(3) !grid cell number
    double precision :: dx1, dy1, dz1
    double precision :: getdl(3)
    integer :: i, j, k
    double precision :: xh, yh, zh !x, y and z at the midpoint of the cell

    rcheck = r
    call edgecheck(rcheck)

    i = floor(rcheck(1))
    j = floor(rcheck(2))
    k = floor(rcheck(3))
    
    dx1 = x(i+1)-x(i)
    dy1 = y(j+1)-y(j)
    dz1 = z(k+1)-z(k)
    
    xh = x(i) + dx1/2
    yh = y(j) + dy1/2
    zh = z(k) + dz1/2

    if (coord_type == 1) then !Cartesian coordinates
      getdl(1) = dx1
      getdl(2) = dy1
      getdl(3) = dz1
    else if (coord_type == 2) then !spherical coordinates
      getdl(1) = dx1 ! dr
      getdl(2) = xh*dy1 ! r d(theta)
      getdl(3) = xh*sin(yh)*dz1 ! r sin(theta) d(phi)
    else if (coord_type == 3) then !cylindrical coordinates
      getdl(1) = dx1 ! dR
      getdl(2) = xh*dy1 ! R d(phi)
      getdl(3) = dz1 ! dz
    else
      !print*, 'Unknown coordinate system.'
      !print*, 'Please select a valid one'
      stop
    endif
    
  end function

end module
