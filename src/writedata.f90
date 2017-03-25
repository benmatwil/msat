program writedata

  use params

  implicit none

  ! reading in filenames
  character(50) :: arg, outname
  integer(int32) :: iarg, oc

  ! loop indices in three directions
  integer(int32) :: ix, iy, iz
  ! three components of field
  real(np), allocatable, dimension(:,:,:) :: bx, by, bz
  !coordinates of grid
  real(np), allocatable, dimension(:) :: x, y, z

  ! grid parameters
  integer(int32), parameter :: nx = 201, ny = 201, nz = 601 ! size of grid
  real(np), parameter :: xxmin = -2.5_np, xxmax = 2.5_np  !xrange
  real(np), parameter :: yymin = -2.5_np, yymax = 2.5_np  !yrange
  real(np), parameter :: zzmin = -7.5_np, zzmax = 7.5_np  !zrange

  ! specific field parameters
  real(np), allocatable, dimension(:,:,:) :: ee
  real(np), parameter :: xc = 0.0_np, yc = 0.0_np, zc = 0.0_np, aa = 0.5_np, ll = 1.0_np
  real(np), parameter :: b0 = 1.0_np, L = 1.0_np, z0 = 5.0_np
  real(np), parameter :: b1 = 20.0_np, tt = 0.0_np
 
  ! form of current ring magnetic field (without time dependence)
  ! Bx = (b0/(L^2))*x*(z - 3*z0) 
  ! By = (b0/L^2)*y*(z + 3*z0) 
  ! Bz = (b0/L^2)*(z0^2 - z^2 + x^2 + y^2) 

  allocate(bx(nx, ny, nz), by(nx, ny, nz), bz(nx, ny, nz))
  allocate(ee(nx, ny, nz))
  allocate(x(nx), y(ny), z(nz))

  do iz = 1, nz
    z(iz) = zzmin + (zzmax - zzmin)*(iz - 1)/(nz - 1)
    do iy = 1, ny
      y(iy) = yymin + (yymax - yymin)*(iy - 1)/(ny - 1)
      do ix = 1, nx
        x(ix) = xxmin + (xxmax - xxmin)*(ix - 1)/(nx - 1)
        ee(ix, iy, iz) = exp(-((x(ix) - xc)/aa)**2-((y(iy) - yc)/aa)**2-((z(iz) - zc)/ll)**2)
        bx(ix, iy, iz) = (b0/(L**2))*x(ix)*(z(iz) - 3*z0) - &
          tt*(2*b1*(y(iy) - yc)/aa)*ee(ix, iy, iz) 
        by(ix, iy, iz) = (b0/(L**2))*y(iy)*(z(iz) + 3*z0) + &
          tt*(2*b1*(x(ix) - xc)/aa)*ee(ix, iy, iz) 
        bz(ix, iy, iz) = (b0/(L**2))*(z0**2 - z(iz)**2 + 0.5_np*x(ix)**2 + 0.5_np*y(iy)**2)
      enddo
    enddo
  enddo

  oc = 0
  if (command_argument_count() > 0) then
    do iarg = 1, command_argument_count()
      call get_command_argument(iarg,arg)
      if (trim(arg) == '-o') then
        call get_command_argument(iarg+1,arg)
        outname = trim(arg)
        oc = 1
      endif
    enddo
  endif

  if (oc == 1) then
    outname = 'data/'//trim(outname)
  else
    outname = 'data/magfield.dat'
  endif

  open(unit=10, file=trim(outname), access='stream', status='replace')
    write(10) nx, ny, nz
    write(10) bx, by, bz
    write(10) x, y, z
  close(10)

  print*, 'Done, wrote field to '//trim(outname)

end program
