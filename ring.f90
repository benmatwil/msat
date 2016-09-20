module ring
!this module contains parameters involved in tracing rings, and subroutines to add points or remove points from rings. Also included is the subroutine (at_null) which determines if a point on the ring has reached a null.

use params
use common
use trace

!ring drawing parameters
integer, parameter :: nstart = 100 !number of startpoints in ring

double precision :: maxdist, mindist, nulldist

double precision, allocatable, dimension(:,:) :: line1, line2, add1, add2
integer, allocatable, dimension(:) :: break, association, remove, endpoints, nearnull

contains

  subroutine add_points(nlines,iteration)
    !adds points to rings if required
    implicit none
    
    integer :: nlines
    integer :: iline, nxtline, iadd, nadd
    integer, intent(in) :: iteration
    double precision :: b(3), maxdist0

    allocate(add1(3,nlines), add2(3,nlines))

    !test for gaps. Where gaps need to be filled, put this info into 'add'
    add1 = 0d0
    add2 = 0d0
    do iline = 1, nlines !loop over each point in ring
      maxdist0 = maxdist
      if (break(iline) == 0) then
        if (iline < nlines) then !if not the last point
          nxtline = iline+1
        else
          nxtline = 1
        endif
        if (nearnull(iline) == 1 .or. nearnull(nxtline) == 1) maxdist0 = maxdist/2d0
        if (dist(line2(:,iline),line2(:,nxtline)) > maxdist0) then !if two adjacent points too far away
          add1(:,iline) = line1(:,iline) + 0.5d0*(line2(:,nxtline)-line2(:,iline)) !add point half way between two points
          add2(:,iline) = line2(:,iline) + 0.5d0*(line2(:,nxtline)-line2(:,iline))
        else
          add1(:,iline) = 0d0 !don't add anything
          add2(:,iline) = 0d0
        endif
      endif
    enddo
      
    !add new points
    !where the 'add' array has a point to be added, add this point
    nadd = 1
!    !$OMP PARALLEL default(shared) firstprivate(iadd, nadd, iline)
    do iline = 1, nlines
      if (modulus(add1(:,iline)) > 1d-4) then !if |add(:,iline)| is not zero, all points > (1,1,1)
        iadd = iline + nadd
!        !$OMP SECTIONS
!        !$OMP SECTION
        call add_vector(line1,add1(:,iline),iadd) !add elements to every array
!        !$OMP SECTION
        call add_vector(line2,add2(:,iline),iadd)
!        !$OMP SECTION
        call add_element(break,0,iadd)
!        !$OMP SECTION
        if (iline /= nlines) then
          call add_element(association,association(iadd-1),iadd)
        else
          call add_element(association,1,iadd)
        endif
!        !$OMP END SECTIONS
        nadd = nadd + 1
        !print*, 'add', iadd, nadd, size(line1,2), nlines
      endif
    enddo
!    !$OMP END PARALLEL

    nlines = size(line1, 2)

    deallocate(add1, add2)
    
  end subroutine

  !********************************************************************************

  subroutine remove_points(nlines,iteration)
    !removes points from rings as required
    implicit none

    integer :: nlines
    integer :: iline, nxtline, iremove, nremove!, prevline
    integer, intent(in) :: iteration

    allocate(remove(nlines))
    remove = 0

    do iline = 1, nlines
      if (endpoints(iline) == 1) then !remove points that have left the simulation
        remove(iline) = 1
        if (iline /= 1) then
          break(iline-1) = 1
        else
          break(nlines) = 1
        endif
      endif
    enddo

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
    enddo

    !remove points
    nremove = 0
    if (nlines > nstart .or. sum(endpoints) /= 0) then !if the number of points isn't too small...
!      !$OMP PARALLEL default(shared) firstprivate(iremove, nremove, iline)
      do iline = 1, nlines
        if (nlines <= nstart .and. sum(endpoints) == 0) exit
        if (remove(iline) == 1) then !if point is flagged to be removed, then remove
          iremove = iline - nremove
!          !$OMP SECTIONS
!          !$OMP SECTION
          call remove_vector(line1,iremove)
!          !$OMP SECTION
          call remove_vector(line2,iremove)
!          !$OMP SECTION
          call remove_element(association,iremove)
!          !$OMP SECTION
          call remove_element(break,iremove)
!          !$OMP END SECTIONS
          nremove = nremove + 1
          !print*, 'remove', iline, iremove, nremove, size(line1, 2), nlines
        endif
      enddo
!      !$OMP END PARALLEL
    endif

    nlines = size(line1, 2)

    deallocate(remove)

  end subroutine

  !********************************************************************************

  subroutine sep_detect(nlines,nullnum,nring)
  !Determine whether a point in a ring is 'near' a null
  !If so, integrate it and its neighbours forward to see if any of them diverge around the null (i.e. a separator)

    implicit none
    integer :: nlines
    integer :: nullnum !the null the fan is being drawn from
    integer :: inull, iline, k
    integer :: nring

    double precision :: h, h0, tracedist, checkdist

    double precision, allocatable :: r(:,:)
    integer, allocatable :: rmap(:)
    integer :: near(nlines), notnear(nlines), endgap, gapsize
    integer :: index, count, nr, nnc, nextra, n1, n2
    integer, allocatable :: signof(:)

    h0 = 1d-2
    tracedist = 3*nulldist
    checkdist = tracedist + 3*h0

    !$OMP PARALLEL DO default(shared) private(near, notnear, nnc, k, gapsize, endgap, nr, r, rmap, signof, index, count, h)
    do inull = 1, size(rnulls,2)
      if (inull == nullnum) cycle !ignore the null the points belong to
      if (signs(inull)*signs(nullnum) /= 1) then !ignore nulls of the same sign (but not nulls with zero/undetermined sign - just in case)
        near = 0
        notnear = 0

        ! first find all points that lie within nulldist and note their index
        do iline = 1, nlines
          if (dist(line1(:,iline),rnulls(:,inull)) < nulldist &
            .or. dist(line1(:,iline),rnullsalt(:,inull)) < nulldist) near(iline) = iline !point lies within nulldist
        enddo

        if (maxval(near) > 0) then !if there are any points that lie within this distance
          ! find the longest consectutive number of points not near to the null (should be more not near another null)
          nnc = 0
          do iline = 1, nlines
            if (near(iline) == 0) then
              nnc = nnc + 1
              notnear(iline) = nnc
            else
              notnear(iline) = 0
              nnc = 0
            endif
          enddo
          if (nnc /= 0) then
            iline = 1
            do while (near(iline) == 0)
              nnc = nnc + 1
              notnear(iline) = nnc
              iline = iline + 1
            enddo
          endif

          endgap = maxloc(notnear,1) ! location of last non-near point
          gapsize = notnear(endgap) ! number of points in biggest gap not near null
          nextra = 10 !number of points to test either side of test points
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
            rmap = [(k,k=n1,n2)]
          else
            r(:,1:nlines-n1+1) = line1(:,n1:nlines)
            rmap(1:nlines-n1+1) = [(k,k=n1,nlines)]
            r(:,nlines-n1+2:nr) = line1(:,1:n2)
            rmap(nlines-n1+2:nr) = [(k,k=1,n2)]
          endif

          allocate(signof(nr))
          signof = 0
          
          do index = 1, nr
            !extrapolate points along fieldlines
            count = 0
            do while (dist(r(:,index),rnulls(:,inull)) < tracedist .and. count < 1000)
              h = h0
              call trace_line(r(:,index),signs(nullnum),h)
              call edgecheck(r(:,index))
              count = count+1
            enddo

            if (signs(inull)*signs(nullnum) == -1) then              
              !check which side of the null the points end out on
              if (dot(spines(:,inull),r(:,index)-rnulls(:,inull)) > 0) then
                signof(index) = 1
              else
                signof(index) = -1
              endif
              
              !if theres a change in sign, theres the separator
              if (index /= 1) then
                if (signof(index-1)*signof(index) == -1 .and. break(rmap(index-1)) /= 1 &
                  .and. dist(rnulls(:,inull), r(:,index)) < checkdist &
                  .and. dist(rnulls(:,inull), r(:,index-1)) < checkdist) then
                  print*, 'Found a separator', nring, nlines
                  break(rmap(index-1)) = 1 !disassociate points so that new points don't get added between them as they diverge around the null
                  nseps = nseps + 1
                  ! print*, nr, nulldist, nlines
                  ! print*, r(:,1)
                  ! print*, line1(:,rmap(1))
                  ! print*, r(:,index-1)
                  ! print*, line1(:,rmap(index-1))
                  ! print*, r(:,index)
                  ! print*, line1(:,rmap(index))
                  ! print*, r(:,nr)
                  ! print*, line1(:,rmap(nr))

                  !write the point's information to the separator file
                  write(12) 1
                  write(12) nullnum, inull
                  write(12) nring, rmap(index-1)
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
