module params

  implicit none

  ! basic parameters/constants
  double precision, parameter :: pi = 4.d0*atan(1.d0)
  double precision, parameter :: dtor = pi/180.d0

  ! data file containing the magnetic field data
  character(100), parameter :: defaultfilename = 'magfield.dat'
  character(100), parameter :: datadir = 'data'
  character(100), parameter :: outputdir = 'output'
  character(100) :: filein, fileout

  ! Number of threads to be used by OpenMP
  integer, parameter :: nproc = 2

  ! nullfinder parameters
  double precision, parameter :: zero = 1d-10 !what the code treats as zero
  integer, parameter :: sig_figs = 6 !the number of significant figures of accuracy required for the null

  ! spinefinder parameters
  double precision, parameter :: rsphere = 1d-4 !radius of sphere which the dot products are calculated on

  ! ssfind parameters
  integer, parameter :: nstart = 100 ! number of startpoints in ring
  integer, parameter :: ringsmax = 10000 !maximum # of rings
  integer, parameter :: pointsmax = 100000 !maximum number of points in ring
  double precision, parameter :: tol = 1d-5 !tolerance of rkf45 scheme
  double precision, parameter :: stepmin = 3d-4 !minimum step length

  ! coordinate type
  ! (1) = Cartesian (x,y,z)
  ! (2) = Spherical (r,theta,phi) (theta=polar angle, phi=azimuthal angle)
  ! (3) = Cylindrical (rho,phi,z)
  integer, parameter :: coord_type = 1

  contains
  
    subroutine filenames

      character(100) :: arg
      integer :: iarg, ichar

      fileout = defaultfilename
      if (command_argument_count() > 0) then
        do iarg = 1, command_argument_count()
          call get_command_argument(iarg,arg)
          if (arg(1:5) == 'data=') then
            fileout = trim(arg(6:))
          endif
        enddo
      endif

      filein = 'data/'//trim(fileout)
      ichar = index(fileout, '.dat')
      fileout = fileout(1:ichar-1)

    end

end module
