module params

  use iso_fortran_env, only : np => real64

  implicit none

  ! basic parameters/constants
  real(np), parameter :: pi = 4.0_np*atan(1.0_np)
  real(np), parameter :: dtor = pi/180.0_np

  ! data file containing the magnetic field data
  character(100), parameter :: defaultfilename = 'data/magfield.dat'
  character(100), parameter :: outputdir = 'output'
  character(100) :: filein, fileout

  ! Number of threads to be used by OpenMP
  integer, parameter :: nproc = 2

  ! nullfinder parameters
  real(np), parameter :: zero = 1e-10_np !what the code treats as zero
  integer, parameter :: sig_figs = 6 !the number of significant figures of accuracy required for the null

  ! spinefinder parameters
  real(np), parameter :: rsphere = 1e-4_np !radius of sphere which the dot products are calculated on

  ! ssfind parameters
  integer, parameter :: nstart = 100 ! number of startpoints in ring
  integer, parameter :: ringsmax = 10000 !maximum # of rings
  integer, parameter :: pointsmax = 100000 !maximum number of points in ring
  real(np), parameter :: tol = 1e-5_np !tolerance of rkf45 scheme
  real(np), parameter :: stepmin = 3e-4_np !minimum step length

  ! coordinate type
  ! (1) = Cartesian (x,y,z)
  ! (2) = Spherical (r,theta,phi) (theta=polar angle, phi=azimuthal angle)
  ! (3) = Cylindrical (rho,phi,z)
  integer, parameter :: coord_type = 2

  contains
  
    subroutine filenames

      character(100) :: arg, outname
      integer :: iarg, ifname, ic, oc

      ic = 0
      oc = 0
      
      outname = outputdir
      filein = trim(defaultfilename)
      
      if (command_argument_count() > 0) then
        do iarg = 1, command_argument_count()
          call get_command_argument(iarg,arg)
          if (trim(arg) == '-i') then
            call get_command_argument(iarg+1,arg)
            filein = trim(arg)
            ic = 1
          elseif (trim(arg) == '-o') then
            call get_command_argument(iarg+1,arg)
            outname = trim(arg)
          elseif (trim(arg) == '-of') then
            call get_command_argument(iarg+1,arg)
            outname = trim(arg)
            oc = 1
          endif
        enddo
      endif
      
      if (oc == 0) then
        fileout = filein(1:index(filein(1:index(filein, '/', .true.)-1), '/', .true.)) &
          //trim(outname)//'/' &
          //trim(filein(index(filein, '/', .true.)+1:index(filein, '.dat', .true.)-1))
      else
        if (index(outname, '/', .true.) /= len(trim(outname))) outname = trim(outname)//'/'
        fileout = trim(outname)//trim(filein(index(filein, '/', .true.)+1:index(filein, '.dat', .true.)-1))
      endif

    end

end module
