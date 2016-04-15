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
  double precision, dimension(3) :: spine, fan

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
  do i = 125,125!1, nnulls
    print*, 'Evaluating null', i,' of', nnulls
    rnull = rnulls(:,i)
    
    call get_properties(sign,spine,fan,spiral,warning)
    
    signs(i) = sign
    spines(:,i) = spine
    fans(:,i) = fan
    spirals(i) = spiral
    warnings(i) = warning
    
    print*, ''
    print*, '-------------------------------------------------------------------------'
    print*, '-------------------------------------------------------------------------'
    print*, ''
  enddo
  stop
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
  use sfmod_converge

  implicit none
  integer :: i, j, count, n, imin, nfw, nbw, nspine, nfan, nfanchk, maxcount
  double precision :: r(3), b(3)
  integer :: flag, testgordon, savedata
  double precision :: dphi, dtheta
  double precision :: fact, acc, spinecheck
  double precision :: mindot, dotprod
  integer :: sign, spiral, warning
  integer, allocatable, dimension(:,:) :: densepos, denseposfw, denseposbw

  double precision, dimension(3) :: spine, fan, maxvec, minvec
  
  double precision, dimension(:,:), allocatable :: rconvergefw, rconvergebw, rconvergefw1, rconvergebw1
  double precision, dimension(:,:), allocatable :: rspine, rfan, rfanchk, rfanred, dummy, crossfan
  double precision, dimension(:), allocatable :: roldfw, rnewfw, bnewfw, roldbw, rnewbw, bnewbw

  !set up theta and phi for sphere around null
  dphi = 360.d0/dble(nphi)
  dtheta = 180.d0/dble(ntheta-1)

  do i = 1, nphi
    phis(i) = (i-1)*dphi
  enddo

  do j = 1, ntheta
    thetas(j) = (j-1)*dtheta
  enddo

  phis = phis*dtor
  thetas = thetas*dtor

  print*, 'Null at:', rnull
  print*, 'B=', trilinear(rnull, bgrid)
  print*, '--------------------------------------------------------------------------'

  allocate(rconvergefw(3,nphi*ntheta), rconvergebw(3,nphi*ntheta))
  allocate(roldfw(3), rnewfw(3), bnewfw(3), roldbw(3), rnewbw(3), bnewbw(3))

  fact = 1d-2*rsphere
  acc = 1d-10*rsphere
  ! loop over each point to find converged point
  do j = 1, ntheta
    do i = 1, nphi
      count = 0
      maxcount = 10000
      flag = 0
      r = sphere2cart(rsphere,thetas(j),phis(i))
      b = trilinear(r+rnull, bgrid)
      rnewfw = r
      rnewbw = r
      bnewfw = b
      bnewbw = b
      roldfw = [0,0,0]
      roldbw = [0,0,0]
      do while (count < maxcount) ! do enough times for spine to converge then double it
        call it_conv(roldfw,rnewfw,bnewfw,fact,1)
        call it_conv(roldbw,rnewbw,bnewbw,fact,-1)
        count = count + 1
        if (flag == 0) then
          if (modulus(rnewfw-roldfw) < acc .or. modulus(rnewbw-roldbw) < acc) then
            flag = 1
            maxcount = 2*count
          endif
        endif
      enddo
      n = i+(j-1)*nphi
      rconvergefw(:,n) = normalise(rnewfw)
      rconvergebw(:,n) = normalise(rnewbw)
    enddo
  enddo

  ! remove duplicate vectors which are a distance 1d-4 apart
  ! spine should be left with 2 vectors
  ! fan left with either 2 (minvec eigenvalue too small) or a circle of points
  rconvergefw1 = rconvergefw
  rconvergebw1 = rconvergebw
  call remove_duplicates(rconvergefw, 5d-3)
  call remove_duplicates(rconvergebw, 5d-3)

  nfw = size(rconvergefw,2)
  nbw = size(rconvergebw,2)

  if (nfw == 2 .or. nbw == 2) then ! if one of the reductions only had two vectors, can guarantee one is spine
    if (nfw < nbw) then ! if nfw is smaller then it's the spine
      ! spine going out of null
      rspine = rconvergefw
      rfan = rconvergebw
      sign = -1
    else if (nfw > nbw) then ! if nbw is smaller then it's the spine
      ! spine going into null
      rspine = rconvergebw
      rfan = rconvergefw
      sign = 1
    else if (nfw == nbw) then ! both are equal (hopefully 2=2) so just pick one and check later
      rspine = rconvergebw
      rfan = rconvergefw
      
      ! check whether this guess is correct, otherwise switch
      spinecheck = dot(trilinear(rsphere*rspine(:,1)+rnull,bgrid),rspine(:,1))
      if (abs(dot(trilinear(rsphere*rfan(:,1)+rnull,bgrid),rfan(:,1))) > abs(spinecheck)) then
        dummy = rspine
        deallocate(rspine)
        rspine = rfan
        deallocate(rfan)
        rfan = dummy
        deallocate(dummy)  
        if (spinecheck > 0) then
          sign = -1
        else
          sign = 1
        endif 
      endif
    endif
    spine = rspine(:,1)
  else ! if no reduction is two, use that the spine should be the densest accumulation of points
    call remove_duplicates(rconvergefw1, 1d-1, denseposfw)
    call remove_duplicates(rconvergebw1, 1d-1, denseposbw)
    nfw = maxval(denseposfw)
    nbw = maxval(denseposbw)
    print*, nfw, nbw
    if (nfw < nbw) then
      spine = rconvergebw1(:,maxval(maxloc(denseposbw)))
      rfan = rconvergefw
      rspine = rconvergebw
      sign = 1
    else if (nbw < nfw) then
      spine = rconvergefw1(:,maxval(maxloc(denseposfw)))
      rfan = rconvergebw
      rspine = rconvergefw
      sign = -1
    else !not sure how to to decide in this case, haven't found a case like this yet...
      print*, "Uh oh, can't decide!"
        open(unit=10, file="spinedata.dat", access="stream")
        write(10) size(rconvergefw,2), rconvergefw
        close(10)
        open(unit=10, file="fandata.dat", access="stream")
        write(10) size(rconvergebw,2), rconvergebw
        close(10)
    endif
    deallocate(rconvergefw1, rconvergebw1, denseposfw, denseposbw)
  endif

  ! We have picked rspine and rfan so can get rid of rconverges
  deallocate(rconvergebw, rconvergefw)
  
  if (size(rfan,2) /= 2) then 
    ! Check whether current fan actually will still converge to only 2 points
    rfanchk = rfan
    call remove_duplicates(rfanchk, 1d-1, densepos)
    nfanchk = size(rfanchk,2)
    if (nfanchk <= 6) then ! If so, rfan is now those two points
      print*, "Narrowed down to two fan points"
      deallocate(rfan)
      rfan = rfanchk
      deallocate(rfanchk)
      maxvec = rfan(:,maxval(maxloc(densepos)))
      minvec = normalise(cross(spine,maxvec))
    else ! have a ring/ball
      ! find the cross product of one vector with every otherwise
      ! well converged ring if all cross products are all the same point
      allocate(crossfan(3,nfanchk-1))
      do i = 2, nfanchk
        crossfan(:,i-1) = cross(rfanchk(:,1),rfanchk(:,i))
      enddo
      i = 1
      do while (i <= size(crossfan,2))
        if (modulus(crossfan(:,i)) < 2d-1) then
          call remove_element(crossfan, i)
        else
          crossfan(:,i) = normalise(crossfan(:,i))
          i = i+1
        endif
      enddo
      call remove_duplicates(crossfan,1d-1)
      print*, "Size of crossfan", size(crossfan,2)
      if (allocated(densepos)) deallocate(densepos)
      call remove_duplicates(rfan, 1d-2, densepos) ! look for maxvec in the densest area
      maxvec = rfan(:,maxval(maxloc(densepos,2)))
      if (size(crossfan,2) <= 4) then ! we have a ring, pick minvec to be vector most perpendicular
        print*, "We have a ring"
        mindot = 1
        do i = 1, size(rfan,2)
          dotprod = abs(dot(maxvec,rfan(:,i)))
          if (dotprod < mindot) then
            mindot = dotprod
            imin = i
          endif
        enddo
        minvec = rfan(:,imin)
        print*, "Perp vec is at ", imin, "out of a total of ", size(rfan)
      else ! we have a ball, find vectors approximately perpendicular and pick one in the densest area
        print*, "We have a ball"
        minvec = normalise(cross(spine,maxvec))
        i = 1
        rfanred = rfan
        do while (i <= size(rfanred,2))
          dotprod = abs(dot(maxvec,rfanred(:,i)))
          if (dotprod > 2d-1) then
            call remove_element(rfanred, i)
          else
            i = i+1
          endif
        enddo
        if (allocated(densepos)) deallocate(densepos)
        call remove_duplicates(rfanred, 5d-2, densepos)
        minvec = rfanred(:,maxval(maxloc(densepos)))
      endif
    endif
  else
    print*, "Only have two fan points"
    maxvec = rfan(:,1)
    minvec = normalise(cross(spine,maxvec))
  endif

  nspine = size(rspine,2)
  nfan = size(rfan,2)
  print*, "Number of spines:", nspine
  print*, "Number of fans:", nfan
  print*, "Number of points left in fan:", nfanchk

  ! Save data if there is a problematic null for inspection
  savedata = 1
  if (savedata == 1) then
    open(unit=10, file="spinedata.dat", access="stream")
    write(10) size(rspine,2), rspine, spine
    close(10)
    open(unit=10, file="fandata.dat", access="stream")
    write(10) size(rfan,2), rfan, maxvec, minvec
    if (allocated(crossfan)) write(10) size(crossfan,2), crossfan
    close(10)
  endif
  !stop

  print*, "Maxvec is:", maxvec
  print*, "Minvec is:", minvec

  !if (size(rfan,2) == 2) then ! print out eigenvalues
  !  print*,'-------------------------------------------------------------------------'
  !  print*,dot(trilinear(rsphere*spine+rnull,bgrid),spine)/rsphere
  !  print*,dot(trilinear(rsphere*maxvec+rnull,bgrid),maxvec)/rsphere
  !  print*,dot(trilinear(rsphere*minvec+rnull,bgrid),minvec)/rsphere
  !endif

  print*, '-------------------------------------------------------------------------'
  print*, "Final dot prod is"
  print*, "        min/max           ", "      min/spine          ", "      max/spine       "
  print*, dot(minvec,maxvec), dot(minvec,spine), dot(maxvec,spine)
  !print*, 'Pre-gordon-spine = ', spine
  !print*, 'Pre-gordon-fan =   ', cross(minvec,maxvec)
  print*, '-------------------------------------------------------------------------'

  !sign = -sign !check which sign needs to be which
  spiral = 0 ! potential so non-spiral

  testgordon = 0
  if (testgordon == 1) then
    print*, ''
    print*, "Determining null's properties..."

    !test the estimated properties and adjust them if required. Run this twice to be sure
    call test_null(spine,maxvec,minvec,sign,spiral)
    call test_null(spine,maxvec,minvec,sign,spiral)
  endif

  fan = normalise(cross(minvec,maxvec)) ! fan vector is perp to fan plane

  if (sign .eq. 0) then
    fan = (/0.d0,0.d0,1.d0/)
    spine = (/0.d0,0.d0,1.d0/)
    print*, "Warning, wasn't able to find null properties"
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
