module ring
!this module contains parameters involved in tracing rings, and subroutines to add points or remove points from rings. Also included is the subroutine (at_null) which determines if a point on the ring has reached a null.

use params
use common
use trace

!ring drawing parameters
integer, parameter :: nstart = 100 !number of startpoints in ring

!nulldist decoupled from maxdist1 and moved into params.f90
!double precision, parameter :: nulldist=maxdist1*4. !maximum distance a ring point can be  for it to be treated as being 'at' a null

!double precision, parameter :: mindist1 = maxdist1/4d0 !minimum distance between points in a ring (defined as 1/3 of the maximum distance)
double precision, parameter :: maxdist1, mindist, nulldist

double precision, allocatable, dimension(:,:) :: line1, line2, add1, add2
integer, allocatable, dimension(:) :: break, association, remove, endpoints

contains

  subroutine add_points(nlines,iteration)
    !adds points to rings if required
    implicit none
    
    integer :: nlines
    integer :: iline, nxtline
    double precision :: maxdist
    integer, intent(in) :: iteration
    double precision :: b(3)

    if (iteration < 100) then !if near-in to the starting null we want smaller max/min separations
      maxdist = maxdist1*0.05
    else
      maxdist = maxdist1
    endif

    !test for gaps. Where gaps need to be filled, put this info into 'add'
    add1 = 0
    add2 = 0
    do iline = 1, nlines !loop over each point in ring
      if (break(iline) == 0) then
        if (iline < nlines) then !if not the last point
          nxtline = iline+1
        else
          nxtline = 1
        endif
        if (dist(line2(:,iline),line2(:,nxtline)) > maxdist) then !if two adjacent points too far away
          add1(:,iline) = line1(:,iline) + 0.5*(line2(:,nxtline)-line2(:,iline)) !add point half way between two points
          add2(:,iline) = line2(:,iline) + 0.5*(line2(:,nxtline)-line2(:,iline))
          !if (outedge(add2(:,iline))) print*, 'adding point out'
        else
          add1(:,iline) = 0d0 !don't add anything
          add2(:,iline) = 0d0
        endif
      endif
    enddo
      
    !add new points
    !where the 'add' array has a point to be added, add this point
    iline = 0
    b = [0d0,0d0,0d0]
    do while (iline < nlines)
      iline = iline+1
      if (modulus(add1(:,iline)) > 0) then !if add(:,iline) is not zero..
        !print*, 'point has to be added to',iline
        call add_vector(line1,add1(:,iline),iline+1) !add elements to every array
        call add_vector(line2,add2(:,iline),iline+1)
        call add_vector(add1,b,iline+1)
        call add_vector(add2,b,iline+1)
        call add_element(endpoints,0,iline+1)
        call add_element(remove,0,iline+1)
        call add_element(break,0,iline+1)
        if (iline /= nlines) then
          call add_element(association,association(iline),iline+1)
        else
          call add_element(association,1,iline+1)
        endif
        iline = iline+1
        nlines = nlines+1
      endif
    enddo
  end subroutine

  !********************************************************************************

  subroutine remove_points(nlines,iteration)
    !removes points from rings as required
    implicit none

    integer :: nlines
    integer :: iline, nxtline
    integer, intent(in) :: iteration
    double precision :: mindist

    if (iteration < 100) then
      mindist = mindist1*0.05
    else
      mindist = mindist1
    endif

    !check for too tightly spaced points, flag points to be removed
    do iline = 1, nlines !loop over all points
      if (break(iline) == 0) then
        if (iline < nlines) then !if not end point
          nxtline = iline+1
        else
          nxtline = 1
        endif
        if (iline /= 1) then
          if (dist(line2(:,iline),line2(:,nxtline)) < mindist .and. remove(iline-1) == 0) remove(iline) = 1 !if iline and iline+1 are too near
        endif
      endif
      if (endpoints(iline) == 1) then !remove points that have left the simulation
        remove(iline) = 1
        if (iline /= 1) then
          break(iline-1) = 1
        else
          break(nlines) = 1
        endif
      endif
    enddo

    !do iline = 1, nlines
    !  if (iline == 1) then !if not end point
    !    nxtline = nlines
    !  else
    !    nxtline = iline-1
    !  endif
    !  if (remove(iline) == 1 .and. break(iline) == 1 .and. break(nxtline) == 0) then
    !    if (iline /= 1) then
    !      break(iline-1) = 1
    !    else
    !      break(nlines) = 1
    !    endif
    !    print*, 'going to remove bad point', iteration, iline, nxtline
    !  endif
    !enddo
    !print*, sum(remove), sum(break), iteration
    !do iline = 1, nlines
    !  if (remove(iline) == 1 .or. break(iline) == 1) print*, remove(iline), break(iline)
    !enddo

    do iline = 1, nlines-1
      if (dist(line1(:,iline),line1(:,iline+1)) > 3*maxdist1) break(iline) = 1
    enddo

    !remove points
    iline = 1
    if (nlines > nstart .or. sum(endpoints) /= 0) then !if the number of points isn't too small...
      do while (iline <= nlines)
        if (nlines <= nstart .and. sum(endpoints) == 0) exit
        if (remove(iline) >= 1) then !if point is flagged to be removed, then remove
          call remove_vector(line1,iline)
          call remove_vector(line2,iline)
          call remove_vector(add1,iline)
          call remove_vector(add2,iline)
          call remove_element(endpoints,iline)
          call remove_element(remove,iline)
          call remove_element(association,iline)
          call remove_element(break,iline)
          nlines = nlines-1
        else
          iline = iline+1
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
    integer :: i, j, k, ir
    integer :: nring

    double precision :: h

    double precision, allocatable :: r(:,:)
    integer, allocatable :: rmap(:)
    integer :: near(nlines), notnear(nlines), endgap, gapsize
    integer :: index, count, nr, nnc, nextra, n1, n2
    integer, allocatable :: signof(:)

    !$OMP PARALLEL DO default(shared) private(near, notnear, nnc, k, gapsize, endgap, nr, r, ir, rmap, signof, index, count, h)
    do i = 1, size(rnulls,2)
      if (i == nullnum) cycle !ignore the null the points belong to
      if (signs(i)*signs(nullnum) /= 1) then !ignore nulls of the same sign (but not nulls with zero/undetermined sign - just in case)
        near = 0
        notnear = 0

        ! first find all points that lie within nulldist and note their index
        do j = 1, nlines
          if (dist(line1(:,j),rnulls(:,i)) < nulldist .or. dist(line1(:,j),rnullsalt(:,i)) < nulldist) near(j) = j !point lies within nulldist
        enddo

        if (maxval(near) > 0) then !if there are any points that lie within this distance
          ! find the longest consectutive number of points not near to the null (should be more not near another null)
          nnc = 0
          do k = 1, nlines
            if (near(k) == 0) then
              nnc = nnc + 1
              notnear(k) = nnc
            else
              notnear(k) = 0
              nnc = 0
            endif
          enddo
          if (nnc /= 0) then
            k = 1
            do while (near(k) == 0)
              nnc = nnc + 1
              notnear(k) = nnc
              k = k + 1
            enddo
          endif

          endgap = maxloc(notnear,1) ! location of last non-near point
          gapsize = notnear(endgap) ! number of points in biggest gap not near null
          nextra = 20 !number of points to test either side of test points
          nr = nlines - gapsize + 2*nextra !total number to test (nlines-gapsize near null)

          if (gapsize == 0 .or. nr >= nlines) then
            n1 = 1
            n2 = nlines
            nr = nlines
          else
            !select all points not in this longest chain to be tested for change in side of the fan
            n1 = modulo(endgap + 1 - nextra - 1, nlines) + 1
            n2 = modulo(nlines + endgap - gapsize + nextra - 1, nlines) + 1
          endif
          allocate(r(3,nr), rmap(nr))

          if (n1 <= n2) then
            r = line1(:,n1:n2)
            rmap = [(i,i=n1,n2)]
          else
            r(:,1:nlines-n1+1) = line1(:,n1:nlines)
            rmap(1:nlines-n1+1) = [(i,i=n1,nlines)]
            r(:,nlines-n1+2:nr) = line1(:,1:n2)
            rmap(nlines-n1+2:nr) = [(i,i=1,n2)]
          endif

          allocate(signof(nr))
          signof = 0
          
          do index = 1, nr
            !extrapolate points along fieldlines
            count = 0
            do while (dist(r(:,index),rnulls(:,i)) < 3*nulldist .and. count < 1000)
              h = 1d-2
              call trace_line(r(:,index),signs(nullnum),h)
              call edgecheck(r(:,index))
              count = count+1
            enddo

            if (signs(i)*signs(nullnum) == -1) then              
              !check which side of the null the points end out on
              if (dot(spines(:,i),r(:,index)-rnulls(:,i)) > 0) then
                signof(index) = 1
              else
                signof(index) = -1
              endif
              
              !if theres a change in sign, theres the separator
              if (index /= 1) then
                if (signof(index-1)*signof(index) == -1) then
                  !print*, rmap(index-1), nlines, nlines-rmap(index-1)
                  print*, 'Found a separator', nring
                  break(rmap(index-1)) = 1 !disassociate points so that new points don't get added between them as they diverge around the null
                  nseps = nseps+1

                  !write the point's information to the separator file
                  !print*, 1, nullnum, i, nring, rmap(index-1), nlines
                  write(12) 1
                  write(12) nullnum, i
                  write(12) nring, rmap(index-1)
                  !exit
                endif
              endif
            else
              if (count == 1000) then !then we want to remove points as they appear to be stuck at the null
                remove(rmap(index)) = 1
                break(rmap(index-1)) = 1
                cycle
              endif
            endif
          enddo
          deallocate(r, signof, rmap)
        endif      
      endif
    enddo
    !$OMP END PARALLEL DO

  end subroutine

end module
