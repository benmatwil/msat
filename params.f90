module params

  !basic parameters/constants
  double precision, parameter :: pi = 4.d0*atan(1.d0)
  double precision, parameter :: dtor = pi/180.d0

  !data file containing the magnetic field data
  character (len=100) :: filename = 'data/newmag.dat'

  !nullfinder parameters
  double precision, parameter :: zero = 1.d-10 !what the code treats as zero
  integer, parameter :: sig_figs = 6 !the number of significant figures of accuracy required for the null

  !spinefinder parameters
  double precision, parameter :: rsphere = 1d-4 !radius of sphere which the dot products are calculated on
  double precision, parameter :: rtraceto = 1d-1 !the radius the field line tracer traces lines to (should be larger than rsphere)


  !ssfind parameters
  double precision, parameter :: maxdist1 = 0.15 !maximum distance between points in a ring
  integer, parameter :: ringsmax = 1500 !maximum # of rings
  integer, parameter :: pointsmax = 30000 !maximum number of points in ring
  double precision, parameter :: nulldist = 0.6 !maximum distance a ring point can be for it to be treated as being 'at' a null

  double precision,parameter :: tol = 1.d-5 !tolerance of rkf45 scheme
  double precision,parameter :: stepmin = 0.0003 !minimum step length

  ! coordinate type
  ! (1) = Cartesian (x,y,z)
  ! (2) = Spherical (r,theta,phi) (theta=polar angle, phi=azimuthal angle)
  ! (3) = Cylindrical (rho,phi,z)
  integer, parameter :: coord_type = 1

end module
