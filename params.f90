module params

  ! basic parameters/constants
  double precision, parameter :: pi = 4.d0*atan(1.d0)
  double precision, parameter :: dtor = pi/180.d0

  ! data file containing the magnetic field data
  character (len=100), parameter :: defaultfilename = 'data/magfield.dat'
  character (len=100) :: filename, arg
  integer :: icommand

  ! Number of threads to be used by OpenMP
  integer, parameter :: nproc = 2

  ! nullfinder parameters
  double precision, parameter :: zero = 1.d-10 !what the code treats as zero
  integer, parameter :: sig_figs = 8 !the number of significant figures of accuracy required for the null

  ! spinefinder parameters
  double precision, parameter :: rsphere = 1d-4 !radius of sphere which the dot products are calculated on

  ! ssfind parameters
  integer, parameter :: ringsmax = 10000 !maximum # of rings
  integer, parameter :: pointsmax = 100000 !maximum number of points in ring
  double precision, parameter :: tol = 1.d-5 !tolerance of rkf45 scheme
  double precision, parameter :: stepmin = 0.0003 !minimum step length

  ! coordinate type
  ! (1) = Cartesian (x,y,z)
  ! (2) = Spherical (r,theta,phi) (theta=polar angle, phi=azimuthal angle)
  ! (3) = Cylindrical (rho,phi,z)
  integer, parameter :: coord_type = 1

end module
