!spine finder
program sf
use params
use sfmod

implicit none

double precision, allocatable :: rnulls(:,:), spines(:,:), fans(:,:)
integer, allocatable :: signs(:), spirals(:), warnings(:)
integer :: nnulls

integer :: pcount, ncount, ucount

integer :: sign, spiral, warning
double precision, dimension(3) :: spine,fan!,rnull

integer :: i

!Read in 'null.dat'

open (unit=10,file='output/null.dat',form='unformatted')

read(10) nnulls

allocate(rnulls(3,nnulls),spines(3,nnulls),fans(3,nnulls),signs(nnulls),spirals(nnulls),warnings(nnulls))

do i=1,nnulls
 read (10) rnulls(:,i)
enddo

close(10)

!read in bgrid

open(unit=10,file=filename,access='stream')
read(10) nx, ny, nz
allocate(bgrid(nx,ny,nz,3),x(nx), y(ny), z(nz))
read(10) bgrid
read(10) x, y, z
close(10)

dx=(x(nx)-x(1))/nx
dy=(y(ny)-y(1))/ny
dz=(z(nz)-z(1))/nz



print*,nnulls,' nulls'


!now loop over each null and characterise

do i=1,nnulls
  print*,'Evaluating null',i,' of',nnulls
  rnull=rnulls(:,i)
  
  call get_properties(sign,spine,fan,spiral,warning)
  
  signs(i)=sign
  spines(:,i)=spine
  fans(:,i)=fan
  spirals(i)=spiral
  warnings(i) = warning
enddo



!now write data to nulls.dat

open(unit=10,file='output/nulls.dat',form='unformatted')
write(10) nnulls
write(10) signs
write(10) rnulls,spines,fans
close(10)

print*,''
print*,'Summary:'
print*,'number, sign, spiral,warning'

pcount=0
ucount=0
ncount=0
do i=1,nnulls
print*, i, signs(i),spirals(i), warnings(i)
if (signs(i) .eq. 1) pcount=pcount+1
if (signs(i) .eq. -1) ncount=ncount+1
if (signs(i) .eq. 0) ucount=ucount+1
enddo

print*,'Total number of nulls:',nnulls
print*,'Positive',pcount
print*,'Negative',ncount
print*,'Unknown',ucount
print*,'Warning',sum(warnings)
print*,pcount+ncount+ucount,nnulls



end program




!characterise each null
subroutine get_properties(sign,spine,fan,spiral,warning)
use sfmod

implicit none
integer :: i, j
double precision :: r(3), b(3)
double precision :: dphi, dtheta
double precision :: mn, mx
integer :: mxloc(2), mnloc(2)
!double precision :: sep
integer :: nmax, nmin, nsaddle
double precision :: minvec(3), maxvec(3)
double precision :: flux, crossflux
integer :: sign, spiral, warning

double precision, dimension(3) :: spine, fan!, major, minor




!set up theta and phi
dphi=360.d0/dble(nphi)
dtheta=180./dble(ntheta-1)

do i=1,nphi
  phis(i) = (i-1)*dphi
enddo

do j=1,ntheta
  thetas(j) = (j-1)*dtheta
enddo

phis=phis*dtor
thetas=thetas*dtor


print*,'Null at:',rnull
print*,'B=', trilinear(rnull,bgrid)

! rotate null by arbitrary angles
! thetarot=45.*dtor
! phirot=45.*dtor
! 
! thetarot=0.
! phirot=0.



!create sphere
flux=0
crossflux=0
do j=1,ntheta
  do i=1,nphi
    r = sphere2cart(rsphere,thetas(j),phis(i))
    
   ! r=rotate(r,thetarot,phirot)
    
    b=trilinear(r+rnull,bgrid)
    
    bmap(i,j) = dot(b,r)/modulus(b)/modulus(r)
    btotal(i,j) = dot(b,r)/modulus(r)
    bcross(i,j) = modulus(cross(b,r))/modulus(b)/modulus(r)
    modb(i,j)=modulus(b)
    
    flux=flux+abs(btotal(i,j))*dphi*dtheta*sin(thetas(j))
    crossflux=crossflux + abs(bcross(i,j))*modb(i,j)*dphi*dtheta*sin(thetas(j))
    
  enddo
enddo

!print*, 'FLUX=',flux
!print*,'CROSS',crossflux
!print*,'RATIO',flux/crossflux


!open (unit=10,file='output/bmap.dat',form='unformatted')
!write(10) nphi, ntheta
!write(10) bmap, btotal,bcross,modb
!close(10)


mnloc=minloc(bmap)
mxloc=maxloc(bmap)

!print*, mxloc
!print*, mnloc

mn=minval(bmap)
mx=maxval(bmap)

!get maxima, minima and saddle points
call get_maxima(bmap, maxima,minima,saddle, nmax,nmin,nsaddle)


!open(unit=10,file='output/critical.dat',form='unformatted')
!write(10) nmax
!write(10) maxima(1,:), maxima(2,:)!

!write(10) nmin
!write(10) minima(1,:), minima(2,:)

!write(10) nsaddle
!write(10) saddle(1,:), saddle(2,:)

!close(10)



sign=0
spiral=0

!try to determine sign of nulls and vectors

if (abs(mn) .gt. spiraltol .and. abs(mx) .gt. spiraltol) then !if both signs have maxima/minima at ~1
  !probably a non-spiral null
  
  if ((nmax .eq. 2) .and. (nmin .gt. 2)) then
    !print*, 'Negative non-spiral null'
    sign=-1
    spiral=0
    
    spine=sphere2cart(1.d0,thetas(mxloc(2)),phis(mxloc(1)))
    
    call getminmax(nmin,minima,minvec,maxvec)
    
    !fan=cross(minvec,maxvec)
    
  else if ((nmin .eq. 2) .and. (nmax .gt. 2)) then
    !print*, 'Positive non-spiral null'
    sign=1
    spiral=0
    
    spine=sphere2cart(1.d0,thetas(mnloc(2)),phis(mnloc(1)))
    
    call getminmax(nmax,maxima,minvec,maxvec)

  else if ((nmin .eq. 2) .and. (nmax .eq. 2) .and. (nsaddle .eq. 2)) then
    sign=1 !pick a ramdom sign (we will test whether this sign is correct later on, and swap it if necessary)
    spiral=1
    !print*,'Spiral null - indeterminate sign'
    !try for a positive null
    
    spine=sphere2cart(1.d0,thetas(mnloc(2)),phis(mnloc(1)))
    maxvec=sphere2cart(1.d0,thetas(mxloc(2)),phis(mxloc(1)))
    minvec=sphere2cart(1.d0,thetas(saddle(2,1)),phis(saddle(1,1)))
    
  else if ( (nmax .eq. 2) .and. (nmin .eq. 1) .and. (nsaddle .eq. 1)) then
    sign=-1
    spiral=0
      spine=sphere2cart(1.d0,thetas(mxloc(2)),phis(mxloc(1)))
      maxvec=sphere2cart(1.d0,thetas(mnloc(2)),phis(mnloc(1)))
      minvec=sphere2cart(1.d0,thetas(saddle(2,1)),phis(saddle(1,1)))
  else if ( (nmax .eq. 1) .and. (nmin .eq. 2) .and. (nsaddle .eq. 1)) then
    sign=1
    spiral=0
      spine=sphere2cart(1.d0,thetas(mnloc(2)),phis(mnloc(1)))
      maxvec=sphere2cart(1.d0,thetas(mxloc(2)),phis(mxloc(1)))
      minvec=sphere2cart(1.d0,thetas(saddle(2,1)),phis(saddle(1,1)))
  else if (nmax .eq. 2 .and. nmin .eq. 2) then
    sign=1
    spiral=0
      spine=sphere2cart(1.d0,thetas(mnloc(2)),phis(mnloc(1)))
      maxvec=sphere2cart(1.d0,thetas(mxloc(2)),phis(mxloc(1)))
      minvec=normalise(cross(spine,maxvec))
  else
    sign=0
    spiral=0
    print*, "Odd null - can't characterise"
    !stop
  endif
  
else
  !print*, 'Spiral null'
  spiral=1
  
  if ((mx .gt. spiraltol) .and. (abs(mn) .lt. spiraltol)) then
   ! print*, 'Negative Spiral'
    sign=-1
    spine=sphere2cart(1.d0,thetas(mxloc(2)),phis(mxloc(1)))
    
    if ((nmax .eq. 2) .and. (nmin .eq. 2) .and. (nsaddle .eq. 2)) then
       maxvec=sphere2cart(1.d0,thetas(mnloc(2)),phis(mnloc(1)))
       minvec=sphere2cart(1.d0,thetas(saddle(2,1)),phis(saddle(1,1)))
    else if ((nmax .eq. 2) .and. (nmin .gt. 2)) then
       call getminmax(nmin,minima,minvec,maxvec)
    else if ((nmax .ge. 4) .and. (nmin .eq. 2)) then
       call getminmax(nmax, maxima,minvec,maxvec)
        if( abs(dot(spine,maxvec)) .lt. abs(dot(spine,minvec))) minvec=maxvec
        maxvec=sphere2cart(1.d0,thetas(mnloc(2)),phis(mnloc(1)))
    else if ( (nmax .eq. 2) .and. (nmin .eq. 1) .and. (nsaddle .eq. 1)) then
        maxvec=sphere2cart(1.d0,thetas(mnloc(2)),phis(mnloc(1)))
        minvec=sphere2cart(1.d0,thetas(saddle(2,1)),phis(saddle(1,1)))
    else if ( (nmax .eq. 2) .and. (nmin .eq. 2) ) then
        maxvec=sphere2cart(1.d0,thetas(mnloc(2)),phis(mnloc(1)))
        minvec=normalise(cross(spine,maxvec))
    else
      print*, 'weird null!'
      sign=0
    endif
  else if ((mx .lt. spiraltol) .and. (abs(mn) .gt. spiraltol)) then
    !print*, 'Positive spiral'
    sign=1
    spine=sphere2cart(1.d0,thetas(mnloc(2)),phis(mnloc(1)))
    
    if ((nmax .eq. 2) .and. (nmin .eq. 2) .and. (nsaddle .eq. 2)) then
       maxvec=sphere2cart(1.d0,thetas(mxloc(2)),phis(mxloc(1)))
       minvec=sphere2cart(1.d0,thetas(saddle(2,1)),phis(saddle(1,1)))
    else if ((nmax .eq. 2) .and. (nmin .gt. 2)) then
       call getminmax(nmin,minima,minvec,maxvec)
       maxvec=sphere2cart(1.d0,thetas(mxloc(2)),phis(mxloc(1)))
    else if ((nmin .ge. 4) .and. (nmax .eq. 2)) then
       call getminmax(nmax, maxima,minvec,maxvec)
    else if ( (nmax .eq. 1) .and. (nmin .eq. 2) .and. (nsaddle .eq. 1)) then
        maxvec=sphere2cart(1.d0,thetas(mxloc(2)),phis(mxloc(1)))
        minvec=sphere2cart(1.d0,thetas(saddle(2,1)),phis(saddle(1,1)))
    else if ( (nmax .eq. 2) .and. (nmin .eq. 2) ) then
        maxvec=sphere2cart(1.d0,thetas(mxloc(2)),phis(mxloc(1)))
        minvec=normalise(cross(spine,maxvec))
    else
      print*, 'weird null!'
      sign=0
    endif
    
  else
    sign=0
    spiral=1
    print *, 'strange null'
  endif
  
 ! stop
  
endif  

!if major and minor axes are almost parallel, redefine minor axis as spine x major
!print*,'CROSS',modulus(cross(maxvec,minvec))
if (modulus(cross(maxvec,minvec)) .lt. 0.1) then
  minvec=cross(spine,maxvec)
endif

print*,''
print*,"Determining null's properties..."

!test the estimated properties and adjust them if required. Run this twice to be sure
call test_null(spine,maxvec,minvec,sign,spiral)  
call test_null(spine,maxvec,minvec,sign,spiral)  
!call test_null(spine,maxvec,minvec,sign,spiral)  
fan=normalise(cross(minvec,maxvec))


if (sign .eq. 0) then
  fan = (/0.d0,0.d0,1.d0/)
  spine=(/0.d0,0.d0,1.d0/)
endif
print*, ''
print*,'Sign=',sign
print*,'Spiral=',spiral
print*, 'Spine=',spine
print*, 'Fan=  ',fan
print*, 'Tilt=',abs(90-acos(dot(fan,spine))/dtor)

warning=0

if (abs(90-acos(dot(fan,spine))/dtor) .lt. 10.) then
  print*, 'WARNING: SPINE AND FAN STRONGLY INCLINED'
  warning=1
endif 




end subroutine
