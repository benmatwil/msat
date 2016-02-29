!This code reads in the output from ssfind and writes every 'skip' rings to file to be plotted by the idl routines. It also trace the separators found, and writes them to file for readin by the IDL routine


module readin_mod

integer, parameter :: skip = 10 !number of rings to skip

integer :: nrings, npoints, nringsmax
integer, allocatable :: nperring(:)
integer, allocatable :: assoc(:)
integer :: association
integer*8 :: pos

integer :: preamble
!integer :: uptoring
!integer :: inring

double precision :: r(3)
double precision,allocatable :: ring(:,:)
double precision,allocatable, dimension(:) :: x,y,z


integer :: opt
integer :: nring, index

contains

!gets the position numbers in the file for ring number 'nring' and point number 'index' within that ring. The position numbers are:
! a - the association number (which point in ring (index-1) the point in ring (nring) came from
!p - the position of the ring vector
subroutine position(nring,index,a,p)
integer*8 :: a
integer*8 :: p
integer :: uptoring

!1 whole ring contains 3*np vector points and np association points

uptoring = preamble + (sum(nperring(1:nring-1))*7)
a=uptoring+(index-1)*1
p=uptoring + nperring(nring)*1 + (index-1)*6

a=a*4+1
p=p*4+1



end subroutine

end module



program readin
use readin_mod
use params
implicit none

integer*8 :: a, p
integer :: i, inull, nrcount

integer :: nnulls, null, nullnum_f, nullnum_t

double precision :: sep, sep1, sep2, sep3
double precision, allocatable, dimension(:,:) :: rnulls


character (len=8), parameter :: fmt='(I3.3)'
character (len=5) :: fname

!open (unit=10,file='output/nulls.dat',access='stream')
!read(10) nnulls
!print*,'num nulls',nnulls
!close(10)

open (unit=10,file='output/null.dat',form='unformatted')

read(10) nnulls
print*,'num nulls',nnulls

allocate(rnulls(3,nnulls))

do i=1,nnulls
 read (10) rnulls(:,i)
 print*,'null position',i,rnulls(:,i)
enddo

close(10)


do null=1,nnulls

write(fname,fmt) null

open(unit=12,file='output/separator'//trim(fname)//'.dat',form='unformatted')


open (unit=10,file='output/everything'//trim(fname)//'.dat',access='stream')
read(10) nrings, npoints, nringsmax
allocate(nperring(nringsmax))
read(10) nperring

open (unit=11,file='output/sep'//trim(fname)//'.dat',form='unformatted')

print*,'nrings, npoints, nringsmax'
print*,nrings, npoints, nringsmax

preamble=3+nringsmax !takes you up to end of 'header'

do
	read(12) opt
	if (opt .lt. 0) exit !if no separator in file, or reached the end of the file
	read(12) nullnum_f, nullnum_t
	print*,'Separator starts from null', nullnum_f
	print*,'Separator ends at null', nullnum_t

	read(12) nring,index !ring number and index within the ring the start point of the separator belongs to
	print*,'ring number and index within ring where separator starts'
	print*,nring,index

	allocate(x(nring+1),y(nring+1),z(nring+1))

	print*,'Tracing Separator...'

	do i=nring,1,-1 !trace back each ring to the start

	    call position(i,index,a,p) !get positions in file belonging to the current ring

	    !print*, a, p

	    read(10,pos=a) association !read association for point number index in ring 1
	    read(10,pos=p) r !read position of point index on ring i
	    if (i .eq. nring) then 
                x(nring+1) = rnulls(1,nullnum_t)
                y(nring+1) = rnulls(2,nullnum_t)
                z(nring+1) = rnulls(3,nullnum_t)
               print*,'End null, last point and start of separator'
               print*,x(nring+1),y(nring+1),z(nring+1)
               print*,r
            endif
	    !write ring position vector into x, y and z
	    x(i)=r(1)
	    y(i)=r(2)
	    z(i)=r(3)
	    index=association

	    if (i .eq. 1) print*,r

	    !print*, i,r,association

	enddo

	!print*, association
	!print*, r

	!close(10)

        nring = nring+1

	write(11) nring
	write(11) nringsmax
	write(11) x, y, z
	!close(11)


	deallocate(x,y,z)


enddo

write(11) -1
close(11)
close(12)

print*,'writing rings to file'

open (unit=11,file='output/ringidl'//trim(fname)//'.dat',form='unformatted')

print*,'nrings=',nrings
nrcount=0
do i=1,nrings-2,skip
  nrcount = nrcount+1
enddo
write(11),nrcount

do i=1,nrings-2,skip
!  print*,i
  call position(i,1,a,p)
  allocate(x(nperring(i)),y(nperring(i)),z(nperring(i)),ring(3,nperring(i)))
  read(10,pos=p) ring
  x=ring(1,:)
  y=ring(2,:)
  z=ring(3,:)

  write(11) nperring(i)
  write(11) x,y,z
  deallocate(x,y,z,ring)
enddo
write(11) -1
close(11)

close(10)

deallocate(nperring)
enddo

end
