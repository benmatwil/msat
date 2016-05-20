module ring
!this module contains parameters involved in tracing rings, and subroutines to add points or remove points from rings. Also included is the subroutine (at_null) which determines if a point on the ring has reached a null.

use params
use common
use trace

!ring drawing parameters
integer, parameter :: nstart = 100 !number of startpoints in ring

!nulldist decoupled from maxdist1 and moved into params.f90
!double precision, parameter :: nulldist=maxdist1*4. !maximum distance a ring point can be  for it to be treated as being 'at' a null

double precision, parameter :: mindist1 = maxdist1/4.d0 !minimum distance between points in a ring (defined as 1/3 of the maximum distance)

double precision, allocatable, dimension(:,:) :: line1, line2, add1, add2, remove, endpoints
double precision, allocatable, dimension(:,:) :: break, association

contains

subroutine add_points(nlines,iteration)
!adds points to rings if required
  implicit none
  
  integer :: nlines
  integer :: i, j
  double precision :: mindist, maxdist
  integer, intent(in) :: iteration
  double precision :: b1(3), b2(1)

  if (iteration < 100) then !if close-in to the starting null we want smaller max/min separations
    mindist = mindist1*0.05
    maxdist = maxdist1*0.05
  else
    mindist = mindist1
    maxdist = maxdist1
  endif

  !test for gaps. Where gaps need to be filled, put this info into 'add'
  add1 = 0
  add2 = 0

  do i = 1, nlines !loop over each point in ring
    if (break(1,i) < 0.5) then
      if (i .lt. nlines) then !if not the last point
        j = i+1
      else
        j = 1
      endif
      if (dist(line2(:,i),line2(:,j)) .gt. maxdist) then !if two adjacent points too far away
        add1(:,i) = line1(:,i) + 0.5*(line2(:,j)-line2(:,i)) !add point half way between two points
        add2(:,i) = line2(:,i) + 0.5*(line2(:,j)-line2(:,i))
      else
        add1(:,i) = 0 !don't add anything
        add2(:,i) = 0
      endif
    endif
  enddo
    
  !add new points
  !where the 'add' array has a point to be added, add this point
  i = 0
  b1 = [0d0,0d0,0d0]
  b2 = [0d0]
  do while (i .lt. nlines)
    i = i+1
    if (modulus(add1(:,i)) .gt. 0) then !if add(:,i) is not zero..
      !print*, 'point has to be added to',i
      call add_element(line1,add1(:,i),i+1) !add elements to every array
      call add_element(line2,add2(:,i),i+1)
      call add_element(add1,b1,i+1)
      call add_element(add2,b1,i+1)
      call add_element(endpoints,b2,i+1)
      call add_element(remove,b2,i+1)
      call add_element(break,b2,i+1)
      if (i .ne. nlines) then
        call add_element(association,[association(1,i)],i+1)
      else
        call add_element(association,[1.d0],i+1)
      endif
      i = i+1
      nlines = nlines+1
    endif
  enddo
end subroutine

!********************************************************************************

subroutine remove_points(nlines,iteration)
  !removes points from rings as required
  implicit none

  integer :: nlines
  integer :: i, j
  integer, intent(in) :: iteration
  double precision :: mindist, maxdist


  if (iteration < 100) then
    mindist = mindist1*0.05
    maxdist = maxdist1*0.05
  else
    mindist = mindist1
    maxdist = maxdist1
  endif

  !check for too tightly spaced points, flag points to be removed
  remove = 0

  do i = 1, nlines !loop over all points
    if (break(1,i) < 0.5) then
      if (i .lt. nlines) then !if not end point
        j = i+1
      else
        j = 1
      endif
      if (dist(line2(:,i),line2(:,j)) .lt. mindist) then !if i and i+1 are too close
        if (i /= 1) then
          if (remove(1,i-1) < 0.5) remove(:,i) = 1 !flag this point to be removed, only if adjacent hasn't been flagged
        endif
      endif
    endif
    if (endpoints(1,i) > 0.5) then !remove points that have left the simulation
      remove(:,i) = 1
      if (i .ne. 1) then
        break(:,i-1) = 1
      else
        break(:,nlines) = 1
      endif
    endif
  enddo

  !remove points
  i = 1
  if (nlines > nstart .or. sum(endpoints) > 0.5) then !if the number of points isn't too small...
    do while (i <= nlines)
      if (nlines <= nstart .and. sum(endpoints) < 0.5) exit
      if (remove(1,i) > 0.5) then !if point is flagged to be removed, then remove
        call remove_element(line1,i)
        call remove_element(line2,i)
        call remove_element(add1,i)
        call remove_element(add2,i)
        call remove_element(endpoints,i)
        call remove_element(remove,i)
        call remove_element(association,i)
        call remove_element(break,i)
        nlines = nlines-1
      else
        i = i+1
      endif
    enddo
  endif

end subroutine

!********************************************************************************

subroutine at_null(nlines,nullnum,nring)
!Determine whether a point in a ring is 'near' a null, and if so, integrate it and its neighbours forward to see if any of them diverge around the null (i.e. a separator). 

implicit none
integer :: nlines
integer :: nullnum !the null the fan is being drawn from
integer :: j, k
double precision :: sep
integer :: nring
integer :: icounter

double precision :: h

double precision :: r(3,nlines)
integer :: backindex,frontindex,index,iters
integer :: signof(nlines)
integer :: idx

double precision :: r1(3), r2(3)

double precision :: seps(nlines-1)

integer :: breakpos

nnulls = size(rnulls,2)

signof = 0

do j = 1, nlines !loop over all points in ring
  do k = 1, nnulls !loop over all nulls
    if (k .eq. nullnum) cycle !ignore the null the points belong to
    if (signs(k)*signs(nullnum) .eq. 1) cycle !ignore nulls of the same sign (but not nulls with zero/undetermined sign - just in case)
    seps = 0
    sep = dist(line1(:,j),rnulls(:,k))

    if (sep .lt. nulldist) then !point lies within nulldist
!       print*, 'AT NULL number',k,'index',j,'sep=',sep

      r(:,j) = line1(:,j) !write new dummy array (is it really needed?)

      ! check for neighbouring points that lie within 3*nulldist
      do backindex = j-1, 0, -1 !previous points
        if (dist(line1(:,backindex),rnulls(:,k)) .lt. 3*nulldist) then
          r(:,backindex) = line1(:,backindex)
!           print*,j,backindex,dist(line(:,backindex),rnulls(:,k))
        else
          exit
        endif
      enddo
      do frontindex = j+1, nlines, 1 !following points
        if (dist(line1(:,frontindex), rnulls(:,k)) .lt. 3*nulldist) then
          r(:,frontindex) = line1(:,frontindex)
!           print*,j,frontindex,dist(line(:,frontindex),rnulls(:,k))
        else
          exit
        endif
      enddo

      backindex = backindex+1
      frontindex = frontindex-1

      !integrate points forward until they have all left the null
      do index = backindex, frontindex
        sep = 0
        iters = 0
        do while (sep .lt. 3*nulldist)
          h = 0.01
          call trace_line(r(:,index),1,signs(nullnum),h)
          call edgecheck(r(:,index))
          sep = dist(r(:,index),rnulls(:,k))
          iters = iters+1
          if (iters .gt. 1000) exit
        enddo
        
        !check which side of the null the points end out on (need to know the spine vector of the null)
        if (dot(spines(:,k),r(:,index)-rnulls(:,k)) .gt. 0.) then
          signof(index) = 1
        else
          signof(index) = -1
        endif
  !     print*,j,index,signof(index),iters
      enddo


      !no knowledge of spine needed for this method of finding 'split' around null
      do index = backindex, frontindex
        r1 = normalise(r(:,backindex)-rnulls(:,k))
        r2 = normalise(r(:,index)-rnulls(:,k))

        if (dot(r1,r2) .gt. 0.) then
          signof(index) = 1
        else
          signof(index) = -1
        endif
      enddo


      !an alternative method (guaranteed to find 1 break point)
      seps = 0
      do index = backindex, frontindex-1
        r1 = normalise(r(:,index)-rnulls(:,k))
        r2 = normalise(r(:,index+1)-rnulls(:,k))
        seps(index) = modulus(r1-r2)
      enddo

      breakpos = maxloc(seps,1)

      ! signof(1:breakpos) = 1
      ! signof(breakpos+1:nlines)=-1

      !determine the index where the index+1th point goes to the other side of the null
      do index = backindex+1, frontindex
        if (signof(backindex)*signof(index) .eq. -1) then
!           print*, 'separator at index:',index-1
          break(1,index-1) = 1 !disassociate points so that new points don't get added between them as they diverge around the null
          nseps = nseps+1

          !write the point's information to the separator file
          write(12) 1
          write(12) nullnum, k
          write(12) nring, index-1

!          print*,'position of sep',line(:,index-1)
          do idx = backindex, frontindex
            line1(:,idx) = r(:,idx) !check line
          enddo
          exit
        endif
      enddo
    endif
  enddo
enddo

end subroutine

end module
