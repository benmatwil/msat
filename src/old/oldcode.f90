! reading in bgrid with ghost cells
allocate(bgrid(0:nx+1,0:ny+1,0:nz+1,3))
bgrid = 0
allocate(x(0:nx+1), y(0:ny+1), z(0:nz+1))
read(10), bgrid(1:nx,1:ny,1:nz,1)
read(10), bgrid(1:nx,1:ny,1:nz,2)
read(10), bgrid(1:nx,1:ny,1:nz,3)
read(10) x(1:nx), y(1:ny), z(1:nz)

!bgrid = bgrid(:,ny+1:0:-1,:,:)
!y = y(ny+1:0:-1)

!copy the grid over for periodicity of coordinate systems
x(0) = 2*x(1) - x(2) !x,R,r don't repeat
x(nx+1) = 2*x(nx)-x(nx-1)
if (coord_type == 3) then !cylindrical
  !phi repeating
  bgrid(:,0,:,:) = bgrid(:,ny,:,:)
  bgrid(:,ny+1,:,:) = bgrid(:,1,:,:)
  y(0) = y(1) - (y(ny) - y(ny-1))
  y(ny+1) = y(ny) + (y(1) - y(0))
  z(0) = 2*z(1) - z(2)
  z(nz+1) = 2*z(nz)-z(nz-1)
else if (coord_type == 2) then !spherical
  !phi repeating
  bgrid(:,:,0,:) = bgrid(:,:,nz,:)
  bgrid(:,:,nz+1,:) = bgrid(:,:,1,:)
  !have theta "periodic"
  if (2*(nz/2) == nz) then !if nz (nphi) is even
    !Top boundary
    bgrid(:,0,1:nz/2,:) = bgrid(:,2,nz/2+1:nz,:)
    bgrid(:,0,nz/2+1:nz,:) = bgrid(:,2,0:nz/2,:)
    !Bottom boundary
    bgrid(:,ny+1,1:nz/2,:) = bgrid(:,ny-1,nz/2+1:nz,:)
    bgrid(:,ny+1,nz/2+1:nz,:) = bgrid(:,ny-1,0:nz/2,:)
  else !if nz is odd
    !Top boundary
    bgrid(:,0,1:nz/2,:) = (bgrid(:,2,1+nz/2:nz-1,:) + bgrid(:,2,2+nz/2:nz,:))/2
    bgrid(:,0,nz/2+1,:) = (bgrid(:,2,nz,:) + bgrid(:,2,1,:))/2
    bgrid(:,0,nz/2+2:nz,:) = (bgrid(:,2,1:nz/2,:) + bgrid(:,2,2:nz/2+1,:))/2
    !Bottom boundary
    bgrid(:,ny+1,1:nz/2,:) = (bgrid(:,ny-1,1+nz/2:nz-1,:) + bgrid(:,ny-1,2+nz/2:nz,:))/2
    bgrid(:,ny+1,nz/2+1,:) = (bgrid(:,ny-1,nz,:) + bgrid(:,ny-1,1,:))/2
    bgrid(:,ny+1,nz/2+2:nz,:) = (bgrid(:,ny-1,1:nz/2,:) + bgrid(:,ny-1,2:nz/2+1,:))/2
  endif
  z(0) = 2*z(1) - z(2)
  z(nz+1) = 2*z(nz)-z(nz-1)
  y(0) = y(2)
  y(ny+1) = y(ny-1)
else if (coord_type == 1) then
  z(0) = 2*z(1) - z(2)
  z(nz+1) = 2*z(nz)-z(nz-1)
endif