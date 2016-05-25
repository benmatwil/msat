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

double precision, allocatable, dimension(:,:) :: line1, line2, add1, add2
integer, allocatable, dimension(:) :: break, association, remove, endpoints

contains

subroutine add_points(nlines,iteration)
!adds points to rings if required
  implicit none
  
  integer :: nlines
  integer :: i, j
  double precision :: maxdist
  integer, intent(in) :: iteration
  double precision :: b1(3), b2(1)

  if (iteration < 100) then !if close-in to the starting null we want smaller max/min separations
    maxdist = maxdist1*0.05
  else
    maxdist = maxdist1
  endif

  !test for gaps. Where gaps need to be filled, put this info into 'add'
  add1 = 0
  add2 = 0

  do i = 1, nlines !loop over each point in ring
    if (break(i) == 0) then
      if (i < nlines) then !if not the last point
        j = i+1
      else
        j = 1
      endif
      if (dist(line2(:,i),line2(:,j)) > maxdist) then !if two adjacent points too far away
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
  do while (i < nlines)
    i = i+1
    if (modulus(add1(:,i)) > 0) then !if add(:,i) is not zero..
      !print*, 'point has to be added to',i
      call add_row(line1,add1(:,i),i+1) !add elements to every array
      call add_row(line2,add2(:,i),i+1)
      call add_row(add1,b1,i+1)
      call add_row(add2,b1,i+1)
      call add_element(endpoints,0,i+1)
      call add_element(remove,0,i+1)
      call add_element(break,0,i+1)
      if (i /= nlines) then
        call add_element(association,association(i),i+1)
      else
        call add_element(association,1,i+1)
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
  double precision :: mindist


  if (iteration < 100) then
    mindist = mindist1*0.05
  else
    mindist = mindist1
  endif

  !check for too tightly spaced points, flag points to be removed
  remove = 0

  do i = 1, nlines !loop over all points
    if (break(i) == 0) then
      if (i < nlines) then !if not end point
        j = i+1
      else
        j = 1
      endif
      if (dist(line2(:,i),line2(:,j)) < mindist .and. i /= 1 .and. remove(i-1) == 0) remove(i) = 1 !if i and i+1 are too close
    endif
    if (endpoints(i) == 1) then !remove points that have left the simulation
      remove(i) = 1
      if (i /= 1) then
        break(i-1) = 1
      else
        break(nlines) = 1
      endif
    endif
  enddo

  !remove points
  i = 1
  if (nlines > nstart .or. sum(endpoints) > 0) then !if the number of points isn't too small...
    do while (i <= nlines)
      if (nlines <= nstart .and. sum(endpoints) == 0) exit
      if (remove(i) == 1) then !if point is flagged to be removed, then remove
        call remove_row(line1,i)
        call remove_row(line2,i)
        call remove_row(add1,i)
        call remove_row(add2,i)
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
integer :: i, j, k
integer :: nring

double precision :: h

double precision, allocatable :: r(:,:)
integer :: close(nlines), maxcount(2)
integer :: index, count, nr
integer, allocatable :: signof(:)

close = 0

do i = 1, size(rnulls,2)
  if (i == nullnum) cycle !ignore the null the points belong to
  if (signs(i)*signs(nullnum) == 1) cycle !ignore nulls of the same sign (but not nulls with zero/undetermined sign - just in case)
  ! first find all points that lie within nulldist and note their index
  do j = 1, nlines
    if (dist(line1(:,j),rnulls(:,i)) < nulldist) close(j) = j !point lies within nulldist
  enddo

  if (maxval(close) > 0) then !if there are any points that lie within this distance
    ! find the longest consectutive number of points not near to the null (should be more not near another null)
    count = 0
    maxcount = [0,0]
    nr = nlines
    k = 0
    do while (k <= nr)
      if (close(k) == 0) count = count + 1
      if (close(k) /= 0) then
        if (count > maxcount(1)) maxcount = [count,k-1]
        count = 0
      endif
      if (k == nlines) then
        if (close(1) == 0 .and. close(nlines) == 0) then
          k = 1
          nr = 1
        else
          k = k + 1
        endif
      elseif (nr /= nlines) then
        if (close(k+1) == 0) then
          nr = nr + 1
          k = k + 1
        else
          k = k + 1
        endif
      else
        k = k + 1
      endif
    enddo
    
    !select all points not in this longest chain to be tested for change in side of the fan
    nr = size(line1,2)-maxcount(1)
    if (maxcount(2) - maxcount(1) >= 1) then !biggest gap is not contained in array
      allocate(r(3,nr))
      r(:,1:nlines-maxcount(2)) = line1(:,maxcount(2)+1:nlines)
      r(:,nlines-maxcount(2)+1:nr) = line1(:,1:maxcount(2)-maxcount(1))
    else
      r = line1(:,maxcount(2)+1:nr+maxcount(2))
    endif
    
    
    allocate(signof(nr))
    signof = 0
    do index = 1, nr
      !interpolate points along fieldlines
      count = 0
      do while (dist(r(:,index),rnulls(:,k)) < 3*nulldist .and. count < 1000)
        h = 0.01
        call trace_line(r(:,index),1,signs(nullnum),h)
        call edgecheck(r(:,index))
        count = count+1
      enddo
      
      !check which side of the null the points end out on
      if (dot(spines(:,k),r(:,index)-rnulls(:,k)) > 0) then
        signof(index) = 1
      else
        signof(index) = -1
      endif
      
      !if theres a change in sign, theres the separator
      if (index /= 1) then
        if (signof(index-1)*signof(index) == -1) then
          break(mod(index-1+maxcount(2),nlines)) = 1 !disassociate points so that new points don't get added between them as they diverge around the null
          nseps = nseps+1

          !write the point's information to the separator file
          write(12) 1
          write(12) nullnum, i
          write(12) nring, mod(index-1+maxcount(2),nlines)

          exit
        endif
      endif
    enddo
    deallocate(r, signof)
  endif
enddo

end subroutine

end module
