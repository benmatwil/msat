module ring
!this module contains parameters involved in tracing rings, and subroutines to add points or remove points from rings. Also included is the subroutine (at_null) which determines if a point on the ring has reached a null.

use params
use common
use trace

!ring drawing parameters
integer, parameter :: nstart=100 !number of startpoints in ring

!nulldist decoupled from maxdist1 and moved into params.f90
!double precision, parameter :: nulldist=maxdist1*4. !maximum distance a ring point can be  for it to be treated as being 'at' a null

double precision, parameter :: mindist1=maxdist1/4.d0 !minimum distance between points in a ring (defined as 1/3 of the maximum distance)


double precision, allocatable, dimension(:,:) :: line, add, remove, endpoints
double precision, allocatable, dimension(:,:) :: break, association

 contains


 subroutine add_points(nlines,iteration)
 !adds points to rings if required
    implicit none
    integer :: nlines
    integer :: j
    double precision :: sep
    double precision :: b(3)
    double precision :: oldj
    double precision :: mindist, maxdist
    integer, intent(in) :: iteration
   ! add=0.
    

    if (iteration .le. 100) then !if close-in to the starting null we want smaller max/min separations
      mindist=mindist1*0.05
      maxdist=maxdist1*0.05
    else
      mindist=mindist1
      maxdist=maxdist1
    endif
    
    !test for gaps. Where gaps need to be filled, put this info into 'add'
    add=0.

      do j=1,nlines !loop over each point in ring
        if (break(1,j) .lt. 0.1) then 
        if (edge(line(:,j))) cycle !if point is at the edge of the box

		if (j .lt. nlines) then !if not the last point
		  if (edge(line(:,j+1))) cycle
		  sep=dist(line(:,j),line(:,j+1)) !distance between jth and (j+1)th point
		  if (sep .gt. maxdist) then !if too far away
			add(:,j) = line(:,j) + 0.5*(line(:,j+1)-line(:,j)) !add point half way between two points
		  else
			add(:,j)=0. !don't add anything
		  endif
		else !if the last point (compare with first pint, not (j+1)th
		  if (edge(line(:,1))) cycle
		  sep=dist(line(:,j),line(:,1))
		  if (sep .gt. maxdist) then
			add(:,j) = line(:,j) + 0.5*(line(:,1)-line(:,j))
		  else
			add(:,j)=0.
		  endif
		endif
		endif
      enddo
      
      

      !add new points
      ! where the 'add' array has a point to be added, add this point
      j=0
      do while (j .lt. nlines)
		j=j+1
		b=add(:,j)
		if(modulus(b) .gt. zero) then !if add(:,j) is not zero..
		  !print*, 'point has to be added to',j
		  call add_element(line,b,j+1) !add elements to every array
		  call add_element(add,0*b,j+1)
		  call add_element(endpoints,b,j+1)
		  call add_element(remove,0.*b,j+1)
		  call add_element(break,(/0.d0/),j+1)
		  oldj=association(1,j)
		  if (j .ne. nlines) then
			call add_element(association,(/oldj/),j+1)
		  else
			call add_element(association,(/1.d0/),j+1)
		  endif
		  !print*, association(:,j),association(:,j+1)
		  j=j+1
		  nlines=nlines+1
		  !print*,'add', oldj,b
		  !stop
		endif
      enddo

 end subroutine










 subroutine remove_points(nlines,iteration)
 !removes points from rings as required
 implicit none

 integer :: nlines
 double precision :: sep, test
 integer :: j
 integer, intent(in) :: iteration
 double precision :: mindist,maxdist


 if (iteration .le. 100) then
      mindist=mindist1*0.05
      maxdist=maxdist1*0.05
    else
      mindist=mindist1
      maxdist=maxdist1
    endif

 !check for too tightly spaced points, flag points to be removed
   remove=0.
   do j=1,nlines !loop over all points
    if (break(1,j) .lt. 0.5) then
   ! if (edge(line(:,j))) cycle
    if (j .lt. nlines) then !if not end point
		sep=dist(line(:,j),line(:,j+1)) !distance between jth and j+1th point
		if (sep .lt. mindist) then !if j and j+1 are too close
		  remove(:,j)=1. !flag this point to be removed
		else
		  remove(:,j)=0.!flag point to stay
		  if (edge(line(:,j))) then !remove points that have left the simulation
			remove(:,j)=1.
			if (j .ne. 1) then
			  break(:,j-1) = 1
			else
			  break(:,nlines) = 1
			endif
		  endif
		endif
	else !if j=nlines
      sep=dist(line(:,j),line(:,1))
      if (sep .lt. mindist) then
        remove(:,j) = 1.
      else
        remove(:,j)=0.
      endif
    endif
    endif
  enddo
  
  !if 2 adjacent points are flagged to be removed, only set one of them to be removed
  do j=2,nlines
    if (remove(1,j-1) .gt. 0.5) then
      if (remove(1,j) .gt. 0.5) then
        remove(:,j) = 0.
      endif
    endif
  enddo

  !remove points

  j=1
  if (nlines .gt. nstart) then !if the number of points isn't too small...
  do while (j .lt. nlines)
    !j=j+1
    test=remove(1,j)
    if (nlines .lt. nstart) exit
    if(test .gt. zero) then !if point is flagged to be removed, then remove
      !print*,'remove', j,line(:,j)
      call remove_element(line,j)
      call remove_element(add,j)
      call remove_element(endpoints,j)
      call remove_element(remove,j)
      call remove_element(association,j)
      call remove_element(break,j)

      nlines=nlines-1
    else
      j=j+1
    endif
   enddo
   endif

end subroutine















subroutine at_null(nlines,nullnum,nring)
!Determine whether a point in a ring is 'near' a null, and if so, integrate it and its neighbours forward to see if any of them diverge around the null (i.e. a separator). 

implicit none
integer :: nlines
integer :: nullnum !the null the fan is beign drawn from)
!integer :: ringnum, index
!integer :: nnulls
integer :: j,k
!double precision :: rnull(3)
!double precision :: rnulls(:,:)
double precision :: sep
!double precision, dimension(3) :: rfront, rback, rsep
!integer :: ifront, isep, iback
integer :: nring
!integer :: nits
integer :: icounter

double precision :: h

double precision :: r(3,nlines)
integer :: backindex,frontindex,index,iters
integer :: signof(nlines)
integer :: idx

double precision :: r1(3),r2(3)

double precision :: seps(nlines-1)

integer :: breakpos

!integer ::nullinitial

nnulls=size(rnulls,2)

signof=0



do j=1,nlines !loop over all points in ring
  do k=1,nnulls !loop over all nulls
     if (k .eq. nullnum) cycle !ignore the null the points belong to
     if (signs(k)*signs(nullnum) .eq. 1) cycle !ignore nulls of the same sign (but not nulls with zero/undetermined sign - just in case)
     seps=0
     sep=dist(line(:,j),rnulls(:,k))

     if (sep .lt. nulldist) then !point lies within nulldist
!       print*, 'AT NULL number',k,'index',j,'sep=',sep

       r(:,j)=line(:,j) !write new dummy array (is it really needed?)

       ! check for neighbouring points that lie within 3*nulldist
       do backindex=j-1,0,-1 !previous points
         if (dist(line(:,backindex),rnulls(:,k)) .lt. 3*nulldist) then
           r(:,backindex)=line(:,backindex)
!           print*,j,backindex,dist(line(:,backindex),rnulls(:,k))
         else
           exit
         endif
       enddo
       do frontindex=j+1,nlines,1 !following points
         if (dist(line(:,frontindex),rnulls(:,k)) .lt. 3*nulldist) then
           r(:,frontindex)=line(:,frontindex)
!           print*,j,frontindex,dist(line(:,frontindex),rnulls(:,k))
         else
           exit
         endif
       enddo


       backindex=backindex+1
       frontindex=frontindex-1


       !integrate points forward until they have all left the null
       do index=backindex,frontindex
         sep=0.
         iters=0
         do while (sep .lt. 3*nulldist)
           h=0.01
           call trace_line(r(:,index),1,signs(nullnum),h)
           sep=dist(r(:,index),rnulls(:,k))
           iters=iters+1
           if (iters .gt. 1000) exit
         enddo
         
         !check which side of the null the points end out on (need to know the spine vector of the null)
         if (dot(spines(:,k),r(:,index)-rnulls(:,k)) .gt. 0.) then
           signof(index)=1
         else
           signof(index)=-1
         endif
    !     print*,j,index,signof(index),iters
       enddo


       !no knowledge of spine needed for this method of finding 'split' around null
       do index=backindex,frontindex
         r1=normalise(r(:,backindex)-rnulls(:,k))
         r2=normalise(r(:,index)-rnulls(:,k))

         if (dot(r1,r2) .gt. 0.) then
           signof(index)=1
         else
           signof(index)=-1
         endif
       enddo


       !an alternative method (guaranteed to find 1 break point)
       seps=0.
       do index=backindex,frontindex-1
         r1=normalise(r(:,index)-rnulls(:,k))
         r2=normalise(r(:,index+1)-rnulls(:,k))
         seps(index) = modulus(r1-r2)
       enddo

       breakpos=maxloc(seps,1)

      ! signof(1:breakpos) = 1
      ! signof(breakpos+1:nlines)=-1




       !determine the index where the index+1th point goes to the other side of the null
       do index=backindex+1,frontindex
         if (signof(backindex)*signof(index) .eq. -1) then
!           print*, 'separator at index:',index-1
           break(1,index-1) = 1. !disassociate points so that new points don't get added between them as they diverge around the null
           nseps=nseps+1

           !write the point's information to the separator file
           write(12) 1
           write(12) nullnum,k
           write(12) nring, index-1

 !          print*,'position of sep',line(:,index-1)

           do idx=backindex,frontindex
             line(:,idx)=r(:,idx)
           enddo


           exit
         endif
       enddo


     endif
   enddo
enddo




end subroutine




 end module
