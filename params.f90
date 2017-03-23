module params

  use iso_fortran_env, only : np => real64, int32, int64

  implicit none

  ! data file containing the magnetic field data
  character(100), parameter :: defaultfilename = 'data/magfield.dat'
  character(100), parameter :: outputdir = 'output'

  ! Number of threads to be used by OpenMP
  integer(int32), parameter :: nproc = 4  ! set to zero to use the environment variable OMP_NUM_THREADS

  ! nullfinder parameters
  real(np), parameter :: zero = 1e-10_np ! what the code treats as zero
  integer(int32), parameter :: sig_figs = 6 ! the number of significant figures of accuracy required for the null

  ! spinefinder parameters
  real(np), parameter :: rspherefact = 50.0_np ! radius factor (in grid coordinates) of sphere which to place start points (rsphere = rspherefact*10**-sig_figs)
  integer(int32), parameter :: nphi = 90 ! number of start points in phi direction
  integer(int32), parameter :: ntheta = nphi/2 ! number of start points in theta direction

  ! ssfind parameters
  integer(int32), parameter :: nstart = 500 ! number of startpoints in ring
  integer(int32), parameter :: ringsmax = 10000 ! maximum number of rings
  integer(int32), parameter :: pointsmax = 100000 ! maximum number of points in ring
  real(np), parameter :: stepsize = 0.025_np ! 0.25_np ! step size h after 50 iterations (otherwise 5 times smaller)
  real(np), parameter :: tol = 1e-6_np ! 1e-5_np ! tolerance of rkf45 scheme
  real(np), parameter :: stepmin = 1e-5_np ! 3e-4_np ! minimum step length

  ! NO NEED TO CHANGE BEYOND HERE (basic parameters/constants)
  real(np), parameter :: pi = acos(-1.0_np)
  real(np), parameter :: dtor = pi/180.0_np
  character(100) :: filein, fileout

  contains
  
    subroutine filenames

      character(100) :: arg, outname
      integer(int32) :: iarg, ifname, ic, oc

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

      if (ic == 0) stop 'No input file provided'
      
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
