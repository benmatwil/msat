program writedata_ws
use params

implicit none
integer :: i,j,k
double precision :: b0, L, z0, xc, yc, zc, ll, aa, b1, tt 
double precision, allocatable, dimension(:,:,:) :: bx, by, bz, ee
double precision, allocatable, dimension(:) :: x, y, z    !coordinates of grid

!xyz grid parameters
integer, parameter :: nnx=201, nny=201, nnz=601 !size of grid
double precision, parameter :: xxmin=-2.5, xxmax=2.5  !xrange
double precision, parameter :: yymin=-2.5, yymax=2.5  !yrange
double precision, parameter :: zzmin=-7.5, zzmax=7.5  !zrange

!magnetic field parameters
!form of magnetic field
!Bx = (b0/(L^2))*x*(z - 3*z0) 
!By = (b0/L^2)*y*(z + 3*z0) 
!Bz = (b0/L^2)*(z0^2 - z^2 + x^2 + y^2) 

allocate(bx(nnx,nny,nnz),by(nnx,nny,nnz),bz(nnx,nny,nnz),ee(nnx,nny,nnz))
allocate(x(nnx),y(nny),z(nnz))

b0 = 1.0 
L = 1.0
z0 = 5.0
b1 = 20.0
xc = 0.0
yc = 0.0
zc = 0.0
ll = 1.0
aa = 0.5
tt = 0.75

do k=1,nnz
  z(k) = zzmin + (zzmax-zzmin)*(k-1)/(nnz-1)
  print*, z(k)
  do j=1,nny
    y(j) = yymin + (yymax-yymin)*(j-1)/(nny-1)
    do i=1,nnx
      x(i) = xxmin + (xxmax-xxmin)*(i-1)/(nnx-1)
      ee(i,j,k) = exp(-((x(i)-xc)/aa)**2-((y(j)-yc)/aa)**2-((z(k)-zc)/ll)**2)
      bx(i,j,k) = (b0/(L**2))*x(i)*(z(k) - 3*z0)-tt*(2*b1*(y(j)-yc)/aa)*ee(i,j,k) 
      by(i,j,k) = (b0/(L**2))*y(j)*(z(k) + 3*z0)+tt*(2*b1*(x(i)-xc)/aa)*ee(i,j,k) 
      bz(i,j,k) = (b0/(L**2))*(z0**2 - z(k)**2 + x(i)**2 + y(j)**2)
    enddo
  enddo
enddo


open (unit=10,file=filename,access='stream')
write(10) nnx, nny, nnz
write(10) bx,by,bz
write(10) x, y, z
close(10)

end
