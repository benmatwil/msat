!module for the spinefinder
module sf_gordon_mod
use params
implicit none

!character (len=100) :: fname='output/data.dat'

double precision, parameter :: spiraltol=0.95

integer :: nx, ny, nz
double precision, allocatable, dimension(:) :: x, y, z
double precision :: dx, dy, dz

double precision, allocatable :: bgrid(:,:,:,:)

!double precision, parameter :: rsphere=0.01d0
!double precision, parameter :: rtraceto=0.1d0

integer, parameter :: nphi = 360
integer, parameter :: ntheta = nphi/2

double precision :: bmap(nphi,ntheta)
double precision :: btotal(nphi,ntheta)
double precision :: bcross(nphi,ntheta)
double precision :: modb(nphi,ntheta)
double precision :: thetas(ntheta), phis(nphi)

double precision, dimension(3) :: rnull

double precision :: thetarot, phirot

integer, allocatable :: maxima(:,:), minima(:,:), saddle(:,:)

contains

!********************************************************************************

function trilinear(r,b)
  implicit none
  !find the value of a function, b, at r=(x,y,z) using the 8 vertices
  real*8 :: b(:,:,:,:)
  real*8 :: r(3)
  real*8 :: cube(2,2,2)
  real*8 :: x, y, z
  real*8 :: xp, yp, zp
  real*8 :: f11, f12, f21, f22
  real*8 :: f1, f2
  real*8 :: trilinear(3)
  integer :: nx, ny, nz
  integer :: dims
  integer :: nxpoints,nypoints,nzpoints
  
  nxpoints=size(b,1)
  nypoints=size(b,2)
  nzpoints=size(b,3)

  xp=r(1)
  yp=r(2)
  zp=r(3)

  nx=floor(xp)
  ny=floor(yp)
  nz=floor(zp)
  
  x=xp-nx
  y=yp-ny
  z=zp-nz
  
  do dims=1,3
    cube=b(nx:nx+1,ny:ny+1,nz:nz+1,dims)

    f11=(1-x)*cube(1,1,1) + x*cube(2,1,1)
    f12=(1-x)*cube(1,1,2) + x*cube(2,1,2)
    f21=(1-x)*cube(1,2,1) + x*cube(2,2,1)
    f22=(1-x)*cube(1,2,2) + x*cube(2,2,2)

    f1=(1-y)*f11 + y*f21
    f2=(1-y)*f12 + y*f22

    trilinear(dims) = (1-z)*f1 + z*f2
  enddo
end

!********************************************************************************

function dot(a,b)
  implicit none

  real*8, dimension(3) :: a, b
  real*8 :: dot
  
  dot = sum(a*b)
end function

!********************************************************************************

function cross(a,b)
  implicit none

  real*8, dimension(3) :: cross,a,b
  
  cross(1) = a(2)*b(3) - a(3)*b(2)
  cross(2) = a(3)*b(1) - a(1)*b(3)
  cross(3) = a(1)*b(2) - a(2)*b(1)
end function

!********************************************************************************

function modulus(a)
  implicit none

  real*8 :: a(3)
  real*8 :: modulus

  modulus=sqrt(dot(a,a))

end

!********************************************************************************

function normalise(a)
  implicit none

  real*8 :: normalise(3), a(3)

  normalise=a/modulus(a)
end

!********************************************************************************

!take spherical coordinates and convert to Cartesian
function sphere2cart(r,theta,phi)
  implicit none

  double precision :: theta, phi, r
  double precision :: sphere2cart(3)

  sphere2cart(1) = r*sin(theta)*cos(phi)
  sphere2cart(2) = r*sin(theta)*sin(phi)
  sphere2cart(3) = r*cos(theta)
end function

!********************************************************************************

function rotate(r,theta,phi)
  implicit none

  real*8 :: r(3)
  real*8 :: theta, phi
  real*8 :: rotate(3)
  real*8 :: roty(3,3), rotz(3,3), rot(3,3)

  roty(1,1) = cos(-theta)
  roty(1,2) = 0
  roty(1,3) = -sin(-theta)

  roty(2,1) = 0
  roty(2,2) = 1
  roty(2,3) = 0

  roty(3,1) = sin(-theta)
  roty(3,2) = 0
  roty(3,3) = cos(-theta)

  rotz(1,1) = cos(-phi)
  rotz(1,2) = sin(-phi)
  rotz(1,3) = 0

  rotz(2,1) = -sin(-phi)
  rotz(2,2) = cos(-phi)
  rotz(2,3) = 0

  rotz(3,1) = 0
  rotz(3,2) = 0
  rotz(3,3) = 1

  rot = matmul(roty,rotz)
  rotate = matmul(rot,r)
  r = rotate

end

!********************************************************************************

!calculate angular separation between to vectors (in spherical)
function angularsep(th1,phi1,th2,phi2)
  implicit none

  double precision :: th1, th2, phi1, phi2
  double precision :: angularsep
  double precision :: r1(3), r2(3)

  !produce cartesian vectors (of length 1 for convenience)
  r1=sphere2cart(1.d0,th1,phi1)
  r2=sphere2cart(1.d0,th2,phi2)

  ! theta = acos(r1 \odot r2)
  angularsep = dot(r1,r2)
  angularsep=acos(angularsep)

end

!********************************************************************************

!find maxima/minima/saddle points in a 2d array 'map'
subroutine get_maxima(map,maxima,minima,saddle,nmax,nmin,nsaddle)
implicit none

double precision :: map(:,:)
integer, allocatable :: maxima(:,:), minima(:,:),saddle(:,:)
integer :: nmax, nmin, nsaddle
double precision :: square(3,3)
integer :: horiz(2), vert(2), diagup(2), diagdown(2), isqu(3,3)
integer :: i, j, k, ii, jj


!deallocate arrays if allocated
if (allocated(maxima)) then
  deallocate(maxima, minima, saddle)
endif

allocate(maxima(2,0),minima(2,0),saddle(2,0))

nmax = 0
nmin = 0
nsaddle = 0

do j = 1, ntheta !loop over theta
  if (j > 1 .and. j < ntheta) then !if not at poles
    do i = 1, nphi !loop over phi
      !make 3x3 'square' of points neighbouring the point (i,j)
      if (i .eq. 1) then
        square(2:3,:) = map(1:2,j-1:j+1)
        square(1,:) = map(nphi,j-1:j+1)
      else if (i .eq. nphi) then
        square(1:2,:) = map(nphi-1:nphi,j-1:j+1)
        square(3,:) = map(1,j-1:j+1)
      else
        square(:,:) = map(i-1:i+1,j-1:j+1)
      endif
      
      where (square < square(2,2))
        isqu = -1
      elsewhere (square > square(2,2))
        isqu = 1
      elsewhere (square == square(2,2))
        isqu = 0
      endwhere
      
      if (count(isqu == 1) == 0) then
        nmax = nmax + 1
        call add_element(maxima,(/i,j/),nmax)
      endif
      
      if (count(isqu == -1) == 0) then
        nmin = nmin + 1
        call add_element(minima,(/i,j/),nmin)
      endif
      ! do we want to deal with close maxima/minima
      
      !check for saddle points
      horiz(:) = (/isqu(1,2),isqu(3,2)/)
      vert(:) = (/isqu(2,1),isqu(2,3)/)
      diagup(:) = (/isqu(1,1),isqu(3,3)/)
      diagdown(:) = (/isqu(3,1),isqu(1,3)/)
      
      if ((sum(horiz)*sum(vert) == -4) .or. (sum(diagup) * sum(diagdown) == -4)) then
        nsaddle = nsaddle + 1
        call add_element(saddle,(/i,j/),nsaddle)
      endif
      
      if (i == -1) then !old code
        !subtract (i,j)th value from square to check whether surrounding values are bigger or smaller
        square = square - square(2,2)
        
        ! assign 1 to values bigger, -1 to values smaller and 0 to values equal to centre
        do jj = 1, 3
          do ii = 1, 3
            if (square(ii,jj) .gt. 0) then
              square(ii,jj) = 1
            else if (square(ii,jj) .lt. 0) then
              square(ii,jj) = -1
            else
              square(ii,jj) = 0
            endif
          enddo
        enddo
        
        square(2,2) = 0
        
        !if all surrounding points in square are greater than centre point then centre point is a minimum
        if (sum(square) .gt. 7.5) then
          nmin = nmin + 1
          call add_element(minima,(/i,j/),nmin)
        endif
        
        !centre point is a maximum
        if (sum(square) .lt. -7.5) then
          nmax = nmax + 1
          call add_element(maxima,(/i,j/),nmax)
        endif
      endif
    enddo
  else if (j == 1 .or. j == ntheta) then !if at the north pole or south pole - 'square' needs to be redefined slightly
    i = 1
    if (j == 1) k = j + 1
    if (j == ntheta) k = j - 1
    
    square(2,2) = map(i,j)
    
    square(1,1) = map(nint(0.*dble(nphi)/8.)+1,k)
    square(2,1) = map(nint(1.*dble(nphi)/8.)+1,k)
    square(3,1) = map(nint(2.*dble(nphi)/8.)+1,k)
    square(3,2) = map(nint(3.*dble(nphi)/8.)+1,k)
    square(3,3) = map(nint(4.*dble(nphi)/8.)+1,k)
    square(2,3) = map(nint(5.*dble(nphi)/8.)+1,k)
    square(1,3) = map(nint(6.*dble(nphi)/8.)+1,k)
    square(1,2) = map(nint(7.*dble(nphi)/8.)+1,k)
    
    where (square < square(2,2))
      isqu = -1
    elsewhere (square > square(2,2))
      isqu = 1
    elsewhere (square == square(2,2))
      isqu = 0
    endwhere
    
    if (count(isqu == 1) == 0) then
      nmax = nmax + 1
      call add_element(maxima,(/i,j/),nmax)
    endif
    
    if (count(isqu == -1) == 0) then
      nmin = nmin + 1
      call add_element(minima,(/i,j/),nmin)
    endif

    !check for saddle points  
    horiz(:) = (/isqu(1,2),isqu(3,2)/)
    vert(:) = (/isqu(2,1),isqu(2,3)/)
    diagup(:) = (/isqu(1,1),isqu(3,3)/)
    diagdown(:) = (/isqu(3,1),isqu(1,3)/)
    
    if ((sum(horiz)*sum(vert) == -4) .or. (sum(diagup) * sum(diagdown) == -4)) then
      nsaddle = nsaddle + 1
      call add_element(saddle,(/i,j/),nsaddle)
    endif
  endif
enddo

! remove maxima/minima/saddle points that are too close to each other (sometimes the above algorithm can give multiple nearby critical points)
call remove_duplicates(nsaddle,saddle)
call remove_duplicates(nmax,maxima)
call remove_duplicates(nmin,minima)

print*, nmax, ' Maxima'
print*, nmin, ' Minima'
print*, nsaddle, ' Saddle points'

end subroutine

!********************************************************************************

subroutine add_element(x,val,pos)
  implicit none

  integer, allocatable, dimension(:,:) :: x
  integer :: val(:)
  integer, allocatable, dimension(:,:) :: dummy
  integer :: pos
  integer :: nx, ny

  nx=size(x,1)
  ny=size(x,2)

  allocate(dummy(nx,ny))

  dummy = x

  deallocate(x)
  allocate(x(nx,ny+1))

  if (pos .eq. 1) then
    x(:,1) = val
    x(:,2:ny+1) = dummy
  else if (pos .eq. ny+1 .or. pos .lt. 1) then
    x(:,1:ny) = dummy
    x(:,ny+1) = val
  else
    x(:,1:pos-1) = dummy(:,1:pos-1)
    x(:,pos) = val
    x(:,pos+1:ny+1) = dummy(:,pos:ny)
  endif

  deallocate(dummy)

end

!********************************************************************************

subroutine remove_element(x,pos)
  implicit none

  integer, allocatable, dimension(:,:) :: x
  integer, allocatable, dimension(:,:) :: dummy
  integer :: pos
  integer :: nx, ny

  nx=size(x,1)
  ny=size(x,2)


  allocate(dummy(nx,ny))

  dummy=x

  deallocate(x)
  allocate(x(nx,ny-1))
  
  if (pos .eq. 1) then
    x(:,1:ny-1)=dummy(:,2:ny)
  else if (pos .eq. ny) then
    x(:,1:ny-1)=dummy(:,1:ny-1)
  else
    x(:,1:pos-1)=dummy(:,1:pos-1)
    x(:,pos:ny-1)=dummy(:,pos+1:ny)
  endif
  
  deallocate(dummy)

end

!********************************************************************************
!removes points if they are too close together (angular separation of <5 degrees)
subroutine remove_duplicates(n,vec) !n=number of vectors, vec = 2xn array containing vectors
implicit none

integer, allocatable, dimension(:,:) :: vec
integer :: n
double precision, allocatable, dimension(:,:) :: dist
double precision :: thi, thj, phii, phij
integer :: i, j, nstart

do while (n .gt. 1)
  allocate(dist(n,n))
  !populate distance matrix
  dist=0.
  do j=1,n
    do i=j+1,n
      thi=thetas(vec(2,i))
      thj=thetas(vec(2,j))
      phii=phis(vec(1,i))
      phij=phis(vec(1,j))
      dist(i,j) = angularsep(thi,phii,thj,phij)
    enddo
  enddo
  dist=dist/dtor !convert to degrees
  
  !locate points where distance < 5 degrees. Remove one of these points
  nstart=n
  outdo: do j=1,n
    do i=j+1,n
      if (dist(i,j) .le. 5.) then
        call remove_element(vec,j)
        n=n-1
        exit outdo
      endif
    enddo
  enddo outdo
  
  deallocate(dist)

  if (nstart .eq. n) exit !keep looping until no more points are removed

enddo 

end

!********************************************************************************

!gets the minimum and maximum value of |B| (and the vectors corresponding to their positions). This *should* give the minor and major axes of the fan.
subroutine getminmax(n,vec,minvec,maxvec)
  implicit none

  integer, dimension(:,:) :: vec
  integer :: n
  double precision, dimension(3) :: minvec,maxvec,dumvec
  double precision :: mn, mx
  integer :: i, imax, imin
  double precision :: b, angle
  double precision :: mindist

  !find minima and maxima
  modb = btotal!abs(btotal) !get rid of modb in the end
  mn = maxval(modb)
  mx = minval(modb)
  !print*, maxval(modb), minval(modb)
  do i = 1, n
    b = modb(vec(1,i),vec(2,i))
    if (b .gt. mx) then
      mx = b
      imax = i
      maxvec = sphere2cart(1.d0,thetas(vec(2,i)),phis(vec(1,i)))
    endif
    
    if (b .lt. mn) then
      mn = b
      imin = i
      minvec = sphere2cart(1.d0,thetas(vec(2,i)),phis(vec(1,i)))
    endif
  enddo

  !print*, 'max=',mx, imax,maxvec
  !print*, 'min=',mn, imin,minvec

  angle = acos(dot(minvec,maxvec))/dtor

  if (mn/mx .gt. 0.90) then !if minima and maxima are both similar
    !print*, 'Proper Null'
      
    !try to find points separated by around 90 degrees (to improve accuracy of fan vector)
    mindist = 90.d0
    do i = 1, n
      dumvec = sphere2cart(1.d0,thetas(vec(2,i)),phis(vec(1,i)))
      angle = acos(dot(dumvec,maxvec))/dtor
      if (abs(90-angle) .lt. mindist) then
        minvec = dumvec
        mindist = abs(90.-angle)
      endif
    enddo
    
  else
    !print*, 'improper null'
  endif

end

!********************************************************************************

subroutine test_null(spine,major,minor,sign,spiral)
!Test the assumptions of the null's sign and fan/spine
!first, assuming the guess is right, integrate a ring of points outwards along the expected fan plane, and inwards along it. See if the fan behaves as it should. Then also swap the major axis of the fan and the spine vectors, and try this again. See which one behaves most like a fan

implicit none

double precision, dimension(3) :: spine, major, minor, fan, dumvec
integer :: sign, spiral
integer, parameter :: nring = 1000
double precision :: dot1,dot2,sep1,sep2
double precision, dimension(3) :: v1,v2,fanold,fannew
integer :: i
double precision, dimension(3,nring) :: rings, fans1,spines1,fans2,spines2

!print*,'testing null'

spine=rotate(spine,thetarot,phirot)
major=rotate(major,thetarot,phirot)
minor=rotate(minor,thetarot,phirot)

if (sign .ne. 0) then

  !first try with the null's expected sign
  call get_ring(major,minor,nring,rings)

  ! open (unit=10,file='output/ring1.dat',form='unformatted')
  ! write(10) nring
  ! write(10) rings(1,:),rings(2,:),rings(3,:)
  ! write(10) spine(1),-1*spine(1),spine(2),-1*spine(2),spine(3),-1*spine(3)
  ! write(10) major(1),-1*major(1),major(2),-1*major(2),major(3),-1*major(3)
  ! write(10) minor(1),-1*minor(1),minor(2),-1*minor(2),minor(3),-1*minor(3)
  ! close(10)

  !integrate out along fan (hopefully)
  call integrate_rings(nring,rings,sign)
  !open(unit=90,file='output/trace1.dat',form='unformatted')
  !write(90) nring
  !write(90) rings(1,:),rings(2,:),rings(3,:)

  fans1 = rings

  !integrate back along spine (hopefully) 
  call get_ring(major,minor,nring,rings)
  call integrate_rings(nring,rings,-1*sign)
  ! write(90) nring
  ! write(90) rings(1,:),rings(2,:),rings(3,:)
  ! close(90)
  
  spines1 = rings 
  !
  !see if traced spine is where it should be
  call check_null(nring,rings,spine,dot1,sep1)

  !now try with sign and spine/major axis switched

  call get_ring(spine,minor,nring,rings)

  !open (unit=10,file='output/ring2.dat',form='unformatted')
  !write(10) nring
  !write(10) rings(1,:),rings(2,:),rings(3,:)
  !write(10) major(1),-1*major(1),major(2),-1*major(2),major(3),-1*major(3)
  !write(10) spine(1),-1*spine(1),spine(2),-1*spine(2),spine(3),-1*spine(3)
  !write(10) minor(1),-1*minor(1),minor(2),-1*minor(2),minor(3),-1*minor(3)
  !close(10)

  !integrate along fan
  call integrate_rings(nring,rings,-1*sign)
  !open(unit=90,file='output/trace2.dat',form='unformatted')
  !write(90) nring
  !write(90) rings(1,:),rings(2,:),rings(3,:)
  
  fans2 = rings
  
  !integrate along spine
  call get_ring(spine,minor,nring,rings)
  call integrate_rings(nring,rings,sign)
  !write(90) nring
  !write(90) rings(1,:),rings(2,:),rings(3,:)
  !close(90)
  
  spines2 = rings
  
  !see if traced spine is where it should be
  call check_null(nring,rings,major,dot2,sep2)
  
  !Now see if we need to switch the null's sign
  
  !print*,'Deciding null properties...'
  if ((dot1 .gt. dot2) .and. (sep1 .lt. sep2)) then
    !yay got it correct
    spine = spine

    !check to see if right axes but wrong sign
    call fansrings(fans1,spines1,sign,major,spine)

  else if ((dot2 .gt. dot1) .and. (sep2 .lt. sep1)) then
    print*, 'Swapping sign'
    !boo got it wrong, need to swap sign
    sign = sign*(-1)
    spiral = 1 !if we got it wrong the null was probably spiral
      
    dumvec = spine
    !swap spine and fan
    spine = major
    major = dumvec
      
    call fansrings(fans2,spines2,sign,major,spine)
  else
  !weird null - possibly 2d?
    print*, 'Cannot characterise null - possibly 2D null'
    sign = 0
    spiral = 0
    spine = 0
    minor = 0
    major = 0
  endif

  if (sign .ne. 0) then 
    !refine fan vector (fan is not 100% correct if the null is an inclined spiral)
    
    !trace along fan
    call get_ring(major,minor,nring,rings)
    call integrate_rings(nring,rings,sign)

    !new fan vector is cross product between traced fan vectors 
    v1 = rings(:,1)
    fanold = 0.
    do i = 2, nring
      v2 = rings(:,i)
      fannew = cross(v1,v2)
      if (modulus(fannew) .gt. modulus(fanold)) then
        fanold = fannew
      endif
    enddo

    !print*,'oldfan=',cross(major,minor)
    !print*,'newfan=',normalise(fanold)
    
    !minor=normalise(cross(major,fanold))
    !print*,'nnwfan=',normalise(cross(minor,major))
    major = normalise(v1)
    
    minor = normalise(cross(v1,fanold))
  
  endif
else
  fan = 0
  spine = 0
  major = 0
  minor = 0
endif    
        
end subroutine

!********************************************************************************

subroutine fansrings(fans,spines,sign,major,spine)
implicit none

!determines whether the traced fans/spines make a point (if so, probably the spine)
double precision, dimension (:,:) :: fans, spines
double precision, dimension(3) :: major, spine
integer :: sign
integer :: n, i, j
double precision, allocatable, dimension(:,:) :: fans2, spines2
double precision, dimension(3) :: f1, f2, s1, s2
double precision, dimension(3) :: fbar, sbar
double precision, allocatable, dimension(:,:) :: distmatrix

double precision :: spinedist, fandist

n = size(fans,2)
!print*, 'N=',n

allocate(fans2(3,n), spines2(3,n))

fans2 = fans
spines2 = spines

f1 = fans2(:,1)
s1 = spines2(:,1)

!flip spines and fans so they are all in one hemisphere
do i = 2, n
  f2 = fans2(:,i)
  
  if (dot(f1,f2) .lt. 0) then
    f2 = f2*(-1.d0)
    fans2(:,i) = f2
  endif
  
  s2 = spines2(:,i)
  if (dot(s1,s2) .lt. 0) then
    s2 = s2*(-1.d0)
    spines2(:,i) = s2
  endif
enddo

!get mean coordinate of fans and spines
do i = 1, 3
  fbar(i) = mean(fans2(i,:))
!   fsigma(i)=stddev(fans2(i,:))
  sbar(i) = mean(spines2(i,:))
!   ssigma(i)=stddev(spines2(i,:))
enddo
! 
! print*,modulus(ssigma), modulus(fsigma)

allocate(distmatrix(n,n))

!find maximum distance between spine vectors
distmatrix = 0.d0
do j = 1, n
!print*, j
  do i = j+1, n
    distmatrix(i,j) = modulus(spines2(:,i)-spines2(:,j))
  enddo
enddo
!print*, 'spinesdist=',maxval(distmatrix)
spinedist = maxval(distmatrix)

!maximum distance between fan vectors
distmatrix = 0.d0
do j = 1, n
  do i = j+1, n
    distmatrix(i,j) = modulus(fans2(:,i)-fans2(:,j))
  enddo
enddo
!print*, 'fansdist=',maxval(distmatrix)
fandist = maxval(distmatrix)

if (fandist/spinedist .lt. 2.) then
  print*, 'WARNING: FAN IS VERY POINT-LIKE'
  print*, 'NULL IS EITHER VERY IMPROPER OR 2-DIMENSIONAL'
endif

!if 'fan' seems more point-like than 'spine' then swap them
if (fandist .lt. spinedist) then
  print*,"Fan looks more like a spine. Changing null's sign"
  
  sign = sign*(-1)
  
  spine = fbar
  major = sbar
  
endif

end subroutine

!********************************************************************************

function mean(x)
implicit none

double precision :: mean
double precision :: x(:)
integer :: n

n=size(x)

mean=sum(x)/dble(n)
end function

!********************************************************************************

function stddev(x)
implicit none

double precision :: stddev
double precision :: x(:)
integer :: n

n=size(x)

stddev= sum(x*x)/dble(n) - sum(x)*sum(x)/dble(n)/dble(n)
stddev=sqrt(stddev)
end function

!********************************************************************************

subroutine check_null(nring,rings,spine,dot1,sep1)
implicit none

integer :: nring
double precision, dimension(3) :: spine
double precision, dimension(3,nring) :: rings
double precision :: r(3), dots(nring), dist(nring)
integer :: i
double precision :: dot1, sep1

spine = normalise(spine)

do i = 1, nring
  r = rings(:,i)
  r = normalise(r)
  
  dots(i) = (dot(r,spine))
  if (dots(i) .gt. 0) then
    dist(i) = modulus(r-spine)
  else
    dist(i) = modulus(r+spine)
  endif
enddo

dot1=sum(abs(dots))/nring
sep1=sum(dist)/nring

dot1=minval(abs(dots))
sep1=maxval(dist)

!print*,'dot',sum(abs(dots))/nring,minval(abs(dots))
!print*,'sep',sum(dist)/nring,maxval(dist)
!print*,'checking not complete'


end subroutine


!********************************************************************************


subroutine integrate_rings(nring,rings,sign)
  implicit none
  !integrate rings outwards (sign=sign of null) or inwards (sign= -1* sign of null)
  integer :: nring, sign
  double precision, dimension(:,:) :: rings
  double precision :: r(3)
  integer :: i
  double precision, parameter :: dr = rsphere*0.01
  integer :: icounter
  double precision :: drvec(3)
  double precision :: drreal
  !calculate dr in terms of minimum 'real' distance

  drreal = minval((/dx,dy,dz/))*dr

  ! drvec=(/(1./dx),(1./dy),(1./dz)/)
  ! drvec=normalise(drvec)
  ! drvec=drvec*dr

  drvec = drreal*(/ 1./dx , 1./dy , 1./dz/)

  do i = 1, nring
    ! print*,i,nring,sign
    r = rings(:,i) + rnull
    icounter  = 0
    do while (modulus(r-rnull) .lt. rtraceto)
      r = r + drvec*sign*normalise(trilinear(r,bgrid))
      icounter = icounter+1
      if (icounter .gt. 10000) then
        !print*,'tracer got stuck'
        !stop
        exit
      endif
    enddo
    rings(:,i) = r-rnull
  enddo

end subroutine

!********************************************************************************

!get ring from major and minor axis vectors
subroutine get_ring(major,minor,nring,rings)
  implicit none
  double precision :: major(3), minor(3)
  double precision, dimension(3) :: u, v, w
  integer :: nring
  double precision :: rings(3,nring)
  
  !integer :: nlines
  double precision :: dt, t
  !double precision :: r(3)
  double precision, parameter :: sep=rsphere
  
  integer :: i


  !nlines=size(xs)

  dt = 2.*pi/dble(nring)
  u = major
  v = minor
  ! u = major axis , v = minor axis
  ! w = unitvec((u x v) x u) - vector in plane of ring perpendicular to u
  w = cross(u,v)
  w = cross(w,u)
  w = normalise(w)
  
  !r(t) = u sin(t) + w cos(t) - parameteric description of ring
  
  !generate nlines start points in ring around equator
  do i = 1, nring
    t = (i-1)*dt
    rings(:,i) = u*cos(t)+w*sin(t)
  enddo
  
  rings = rings*sep

  !print*,''
  !print*,'point number, dot product, distance from null'
  !do i=1,nlines
  !  print*,i,r(1)*xs(i)+r(2)*ys(i)+r(3)*zs(i),sqrt(xs(i)**2+ys(i)**2+zs(i)**2)
  !enddo


end subroutine

end module
