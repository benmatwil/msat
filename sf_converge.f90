!spine finder with convergence method
program sf_converge
use params
use sfmod_converge

implicit none

double precision, allocatable :: rnulls(:,:), spines(:,:), fans(:,:)
integer, allocatable :: signs(:), spirals(:), warnings(:)
integer :: nnulls

integer :: pcount, ncount, ucount

integer :: sign, spiral, warning
double precision, dimension(3) :: spine,fan

integer :: i

!Read in 'null.dat'
open (unit=10,file='output/null.dat',form='unformatted')

read(10) nnulls

allocate(rnulls(3,nnulls),spines(3,nnulls),fans(3,nnulls),signs(nnulls),spirals(nnulls),warnings(nnulls))

do i = 1,nnulls
  read(10) rnulls(:,i)
enddo

close(10)

!read in bgrid
open(unit=10,file=filename,access='stream')
read(10) nx, ny, nz
allocate(bgrid(nx,ny,nz,3), x(nx), y(ny), z(nz))
read(10) bgrid
read(10) x, y, z
close(10)

dx = (x(nx)-x(1))/nx
dy = (y(ny)-y(1))/ny
dz = (z(nz)-z(1))/nz

print*, nnulls,' nulls'

!now loop over each null and characterise using get_properties
do i = 1, nnulls
  print*, 'Evaluating null', i,' of', nnulls
  rnull = rnulls(:,i)
  
  call get_properties(sign,spine,fan,spiral,warning)
  
  signs(i) = sign
  spines(:,i) = spine
  fans(:,i) = fan
  spirals(i) = spiral
  warnings(i) = warning
enddo

!now write data to nulls.dat
open(unit=10,file='output/nulls.dat',form='unformatted')
write(10) nnulls
write(10) signs
write(10) rnulls,spines,fans
close(10)

print*, ''
print*, 'Summary:'
print*, 'number, sign, spiral, warning'

pcount = 0
ucount = 0
ncount = 0
do i = 1, nnulls
print*, i, signs(i), spirals(i), warnings(i)
if (signs(i) .eq. 1) pcount = pcount+1
if (signs(i) .eq. -1) ncount = ncount+1
if (signs(i) .eq. 0) ucount = ucount+1
enddo

print*, 'Total number of nulls:', nnulls
print*, 'Positive', pcount
print*, 'Negative', ncount
print*, 'Unknown', ucount
print*, 'Warning', sum(warnings)
print*, pcount+ncount+ucount, nnulls

end program

!********************************************************************************

!characterise each null
subroutine get_properties(sign,spine,fan,spiral,warning)
use sfmod

implicit none
integer :: i, j, k, count, maxcount, n
double precision :: r(3), b(3), roldfw(3), rnewfw(3), bnewfw(3), roldbw(3), rnewbw(3), bnewbw(3), fact, acc
integer, allocatable :: fwflag(:)
double precision :: dphi, dtheta
integer :: sign, spiral, warning

double precision, dimension(3) :: spine, fan
double precision, dimension(:,:), allocatable :: rconvergefw, rconvergebw, rconverge

!set up theta and phi for sphere around null
dphi = 360.d0/dble(nphi)
dtheta = 180./dble(ntheta-1)

do i = 1,nphi
  phis(i) = (i-1)*dphi
enddo

do j = 1,ntheta
  thetas(j) = (j-1)*dtheta
enddo

phis = phis*dtor
thetas = thetas*dtor

print*, 'Null at:', rnull
print*, 'B=', trilinear(rnull, bgrid)

allocate(rconvergefw(3,nphi*ntheta), rconvergebw(3,nphi*ntheta))
allocate(fwflag(nphi*ntheta))

fact = 1d-2*rsphere
acc = 1d-10*rsphere
do j = 1, ntheta
  do i = 1, nphi
    count = 0
    r = sphere2cart(rsphere,thetas(j),phis(i))
    b = trilinear(r+rnull, bgrid)
    rnewfw = r
    rnewbw = r
    bnewfw = b
    bnewbw = b
    roldfw = [0,0,0]
    roldbw = [0,0,0]
    do while (modulus(rnewfw-roldfw) > acc .and. modulus(rnewbw-roldbw) > acc)
      
      roldfw = rnewfw
      rnewfw = rnewfw + fact*normalise(bnewfw)
      rnewfw = rsphere*normalise(rnewfw)
      bnewfw = trilinear(rnewfw+rnull, bgridfw)
      
      roldbw = rnewbw
      rnewbw = rnewbw - fact*normalise(bnewbw)
      rnewbw = rsphere*normalise(rnewbw)
      bnewbw = trilinear(rnewbw+rnull, bgridbw)
      
      count = count + 1
    enddo
    n = i+(j-1)*nphi
    rconvergefw(:,n) = normalise(rnewfw)
    rconvergebw(:,n) = normalise(rnewbw)
    if (modulus(rnewfw-roldfw) > acc) then
      fwflag(n) = 1
    else
      fwflag(n) = 0
    endif
    print*, count
  enddo
enddo

if (count(fwflag == 0) == 0) print*, "all vectors didn't converge" !spine should be contained in backward vectors
if (fwflag(n) == 0) then
  rconverge = rconvergefw
else
  rconverge = rconvergebw
endif

n = size(rconverge,2)
i = 1
do while (i < n-1)
  j = i+1
  do while (j < n)
    if (modulus(rconverge(:,i)-rconverge(:,j)) < 1d-3 .or. modulus(rconverge1(:,i)+rconverge1(:,j)) < 1d-3) remove_element(rconverge,j)
    j = j + 1
  enddo
  i = i + 1
enddo

print*, rconverge
stop
n = size(rconverge,2)
do i = 1, n
  do j = i+1, n
    print*, j, n
    print*, modulus(rconverge1(:,i)-rconverge1(:,j))
    print*, modulus(rconverge1(:,i)+rconverge1(:,j))
  enddo
enddo

print*, rconverge1

print*, ''
print*, "Determining null's properties..."

!test the estimated properties and adjust them if required. Run this twice to be sure
call test_null(spine,maxvec,minvec,sign,spiral)
call test_null(spine,maxvec,minvec,sign,spiral)
fan = normalise(cross(minvec,maxvec))

if (sign .eq. 0) then
  fan = (/0.d0,0.d0,1.d0/)
  spine = (/0.d0,0.d0,1.d0/)
endif
print*, ''
print*, 'Sign =  ', sign
print*, 'Spiral =', spiral
print*, 'Spine = ', spine
print*, 'Fan =   ', fan
print*, 'Tilt =  ', abs(90-acos(dot(fan,spine))/dtor)

warning = 0
if (abs(90-acos(dot(fan,spine))/dtor) .lt. 10.) then
  print*, 'WARNING: SPINE AND FAN STRONGLY INCLINED'
  warning = 1
endif 

end subroutine
