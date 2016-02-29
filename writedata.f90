program writedata
use params

implicit none
integer :: i,j,k
double precision :: a, b, c, aa, bb, cc, dd, tt 
double precision, allocatable, dimension(:,:,:) :: bx, by, bz
double precision, allocatable, dimension(:) :: x, y, z    !coordinates of grid

!xyz grid parameters
integer, parameter :: nnx=206, nny=705, nnz=1409 !size of grid
double precision, parameter :: xxmin=1, xxmax=2.5  !xrange
double precision, parameter :: yymin=0, yymax=2*pi  !yrange
double precision, parameter :: zzmin=pi, zzmax=0  !zrange

!magnetic field parameters
!form of magnetic field
!Bx = (a_par-z)*x + b_par*y + aa_par*tt_par 
!By = c_par*x - (a_par+z)*y + bb_par*tt_par
!Bz = dd_par*tt_par*y + z*z + cc_par*tt_par
!spiral
!double precision, parameter :: a_par=1.0, b_par=1.0, c_par=-2.0 
!double precision, parameter :: aa_par=1.0, bb_par=0.0, cc_par=-2.0, dd_par=3.0
!double precision, parameter :: tt_par=0.1  
!improper
double precision, parameter :: a_par=2.0, b_par=-1.0, c_par=3.0 
double precision, parameter :: aa_par=1.0, bb_par=-1.0, cc_par=-1.0, dd_par=0.0
double precision, parameter :: tt_par=0.4

allocate(bx(nnx,nny,nnz),by(nnx,nny,nnz),bz(nnx,nny,nnz))
allocate(x(nnx),y(nny),z(nnz))

a = a_par
b = b_par
c = c_par
aa = aa_par
bb = bb_par
cc = cc_par
dd = dd_par
tt = tt_par

do k=1,nnz
  z(k) = zzmin + (zzmax-zzmin)*(k-1)/nnz 
  do j=1,nny
    y(j) = yymin + (yymax-yymin)*(j-1)/nny 
    do i=1,nnx
      x(i) = xxmin + (xxmax-xxmin)*(i-1)/nnx 
      bx(i,j,k) = (a-z(k))*x(i) + b*y(j) + aa*tt
      by(i,j,k) = c*x(i) - (a+z(k))*y(j) + bb*tt
      bz(i,j,k) = dd*tt*y(j) + z(k)*z(k) + cc*tt
    enddo
  enddo
enddo

print*,'Mag field data',a,b,c,aa,bb,cc,dd,tt

open (unit=10,file=filename,access='stream')
write(10) nnx, nny, nnz
write(10) bx,by,bz
write(10) x, y, z
close(10)

end
