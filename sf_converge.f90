!spine finder with convergence method
program sf_converge

  use params
  use sfmod_converge

  implicit none

  double precision, allocatable :: rnulls(:,:), spines(:,:), fans(:,:)
  integer, allocatable :: signs(:), spirals(:), warnings(:)
  
  integer :: nnulls, nullstart, nullend
  character(len=12) :: arg

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
  
  if (command_argument_count() > 0) then 
    do i = 1, command_argument_count()
      call get_command_argument(i,arg)
      if (arg(1:2) == 'n=') then
        arg = arg(3:)
        read(arg,*) nullstart
        nullend = nullstart
      endif
    enddo
  else
    nullstart = 1
    nullend = nnulls
  endif

  !now loop over each null and characterise
  do i = nullstart, nullend!1, nnulls
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
  if (nullend - nullstart /= nnulls - 1) stop
  
  !now write data to nulls.dat
  open(unit=10,file='output/nullsben.dat',form='unformatted')
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
  integer :: flag, testgordon, savedata, flagfw, flagbw
  double precision :: dphi, dtheta
  double precision :: fact, acc, spinecheck
  double precision :: mindot, dotprod
  integer :: sign, spiral, warning, signguess
  integer, allocatable, dimension(:,:) :: densepos, denseposfw, denseposbw

  double precision, dimension(3) :: spine, fan, maxvec, minvec
  
  double precision, dimension(:,:), allocatable :: rconvergefw, rconvergebw, rconvergefw1, rconvergebw1
  double precision, dimension(:,:), allocatable :: rspine, rfan, rfanchk, rfanred, dummy, crossfan, distarr
  double precision, dimension(:), allocatable :: roldfw, rnewfw, bnewfw, roldbw, rnewbw, bnewbw, distchk

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
  
  spiral = 0

  print*, 'Null at:', rnull
  print*, 'B=', trilinear(rnull, bgrid)
  print*, '--------------------------------------------------------------------------'

  allocate(rconvergefw(3,nphi*ntheta), rconvergebw(3,nphi*ntheta))
  allocate(roldfw(3), rnewfw(3), bnewfw(3), roldbw(3), rnewbw(3), bnewbw(3))

  ! Working on a sphere of size rsphere
  fact = 1d-2*rsphere
  acc = 1d-10*rsphere
  flagfw = 0
  flagbw = 0
  ! loop over each point to find converged point
  do j = 1, ntheta
    do i = 1, nphi
      count = 0
      maxcount = 5000
      flag = 0
      r = sphere2cart(rsphere,thetas(j),phis(i))
      b = trilinear(r+rnull, bgrid)
      rnewfw = r
      rnewbw = r
      bnewfw = b
      bnewbw = b
      roldfw = [0,0,0]
      roldbw = [0,0,0]
      do while (modulus(rnewfw-roldfw) > acc .and. modulus(rnewbw-roldbw) > acc .and. count < maxcount) ! do enough times for spine to converge
        call it_conv(roldfw,rnewfw,bnewfw,fact,1)
        call it_conv(roldbw,rnewbw,bnewbw,fact,-1)
        count = count + 1
        !if (flag == 0) then
        !  if (modulus(rnewfw-roldfw) < acc .or. modulus(rnewbw-roldbw) < acc) then
        !    flag = 1
        !    maxcount = 1*count ! does this factor need changing?
        !  endif
        !endif
      enddo
      if (modulus(rnewfw-roldfw) < acc) flagfw = flagfw + 1
      if (modulus(rnewbw-roldbw) < acc) flagbw = flagbw + 1
      n = i+(j-1)*nphi
      rconvergefw(:,n) = normalise(rnewfw)
      rconvergebw(:,n) = normalise(rnewbw)
    enddo
  enddo
  
  print*, 'forward stopped', flagfw, 'backward stopped', flagbw
  if (flagfw < flagbw) then
    signguess = 1
  else if (flagfw > flagbw) then
    signguess = -1
  else
    signguess = 0
  endif
  
  ! Now working on a unit sphere of points
  ! remove duplicate vectors which are a distance 1d-4 apart
  ! spine should be left with 2 vectors
  ! fan left with either 2 vectors (minvec eigenvalue too small), a circle of points or just a mess
  rconvergefw1 = rconvergefw
  rconvergebw1 = rconvergebw
  call remove_duplicates(rconvergefw, 1d-4)
  call remove_duplicates(rconvergebw, 1d-4)

  nfw = size(rconvergefw,2)
  nbw = size(rconvergebw,2)

  ! first determine the spine and which is which
  if (nfw == 2 .or. nbw == 2) then ! if one of the reductions only has two vectors, can guarantee one is spine
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
    else if (nfw == nbw) then ! both are equal (i.e. 2=2) so check which converged first and then use eigenvalues as last resort
      if (flagfw < flagbw) then ! backwards converged points stopped first
        rspine = rconvergebw
        rfan = rconvergefw
        sign = 1
      else if (flagfw > flagbw) then ! forwards converged points stopped first
        rspine = rconvergefw
        rfan = rconvergebw
        sign = -1
      else ! still don't know, check the eigenvalues as a last resort - not ideal but rare
        rspine = rconvergebw
        rfan = rconvergefw
        sign = 1
        ! check whether this guess is correct using eigenvalues, otherwise switch
        spinecheck = dot(trilinear(rsphere*rspine(:,1)+rnull,bgrid),rspine(:,1))
        print*, "Eigenvalue at spine", spinecheck
        print*, "Eigenvalue at fan", abs(dot(trilinear(rsphere*rfan(:,1)+rnull,bgrid),rfan(:,1)))
        if (abs(dot(trilinear(rsphere*rfan(:,1)+rnull,bgrid),rfan(:,1))) > abs(spinecheck)) then
          dummy = rspine
          deallocate(rspine)
          rspine = rfan
          deallocate(rfan)
          rfan = dummy
          deallocate(dummy)
          sign = -1
        endif
      endif
    endif
    spine = rspine(:,1)
  else ! if no reduction is two, use that the spine should be the densest accumulation of points
    ! also likely we have a source or a sink so div(B) is non-zero
    call remove_duplicates(rconvergefw1, 5d-2, denseposfw)
    call remove_duplicates(rconvergebw1, 5d-2, denseposbw)
    nfw = maxval(denseposfw)
    nbw = maxval(denseposbw)
    print*, "Density of densest points", nfw, nbw
    if (nfw < nbw) then
      spine = rconvergebw1(:,maxval(maxloc(denseposbw)))
      rfan = rconvergefw
      rspine = rconvergebw
      !sign = 1
    else if (nbw < nfw) then
      spine = rconvergefw1(:,maxval(maxloc(denseposfw)))
      rfan = rconvergebw
      rspine = rconvergefw
      !sign = -1
    else ! not sure how to to decide in this case, haven't found a case like this yet...
      print*, "Uh oh, can't decide!"
        open(unit=10, file="spinedata.dat", access="stream")
        write(10) size(rconvergefw,2), rconvergefw
        close(10)
        open(unit=10, file="fandata.dat", access="stream")
        write(10) size(rconvergebw,2), rconvergebw
        close(10)
    endif
    deallocate(denseposfw, denseposbw)
    sign = 0 ! check to make sure it isn't flagging any proper nulls
  endif

  ! We have picked rspine and rfan so can get rid of rconverges
  deallocate(rconvergebw, rconvergefw, rconvergefw1, rconvergebw1)
  
  print*, "Number of spines:", size(rspine,2)
  print*, "Number of fans:", size(rfan,2)
  
  ! now determine the fan
  if (size(rfan,2) /= 2) then 
    ! Check whether fan actually converges to only 2 points at lower accuracy
    rfanchk = rfan
    call remove_duplicates(rfanchk, 1d-1, densepos)
    nfanchk = size(rfanchk,2)
    print*, "Number of points left in fan:", nfanchk
    if (nfanchk <= 6) then ! If so, rfan is now those two points
      print*, "Narrowed down to two fan points"
      deallocate(rfan)
      rfan = rfanchk
      deallocate(rfanchk)
      maxvec = rfan(:,maxval(maxloc(densepos)))
      minvec = normalise(cross(spine,maxvec))
    else ! have a ring/ball (mess of points)
      ! find the cross product of one vector with every otherwise
      ! well converged ring if all cross products are all the same point
      allocate(crossfan(3,nfanchk-1))
      j = maxval(maxloc(densepos))
      do i = 1, nfanchk !perhaps change this to maxval(densepos) as chosen vector
        if (i < j) crossfan(:,i) = cross(rfanchk(:,j),rfanchk(:,i))
        if (i > j) crossfan(:,i-1) = cross(rfanchk(:,j),rfanchk(:,i))
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
      call remove_duplicates(crossfan,2d-1)
      print*, "Size of crossfan", size(crossfan,2)
      if (allocated(densepos)) deallocate(densepos)
      call remove_duplicates(rfan, 1d-2, densepos) ! look for maxvec in the densest area
      maxvec = rfan(:,maxval(maxloc(densepos,2)))
      if (size(crossfan,2) <= 8) then ! we have a ring, pick minvec to be vector most perpendicular
        print*, "We have a ring"
        ! find the most perpendicular vector to maxvec
        mindot = 1
        do i = 1, size(rfan,2)
          dotprod = abs(dot(maxvec,rfan(:,i)))
          if (dotprod < mindot) then
            mindot = dotprod
            imin = i
          endif
        enddo
        minvec = rfan(:,imin)
        print*, "Perp vec is at ", imin, "out of a total of ", size(rfan,2)
        ! Now check for a full ring for spiralling
        nfan = size(rfan,2)
        allocate(distarr(nfan,nfan))
        do i = 1, nfan
          do j = 1, nfan
            distarr(i,j) = modulus(rfan(:,i)-rfan(:,j))
          enddo
        enddo
        if (maxval(minval(distarr,2,distarr > 0)) < 1d-1) spiral = 1
        deallocate(distarr)
      else ! we have a ball, find vectors approximately perpendicular and pick one in the densest area
        print*, "We have a ball"
        ! find most set of almost perpendicular vectors and then pick densest packed one
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

  ! Save data if there is a problematic null for inspection
  savedata = 0
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

  print*, '-------------------------------------------------------------------------'
  print*, "Final dot prod is"
  print*, "        min/max           ", "      min/spine          ", "      max/spine       "
  print*, dot(minvec,maxvec), dot(minvec,spine), dot(maxvec,spine)
  print*, '-------------------------------------------------------------------------'

  !sign = -sign !check which sign needs to be which

  fan = normalise(cross(minvec,maxvec)) ! fan vector is perp to fan plane

  if (sign .eq. 0) print*, "Warning, sign = 0 and null likely a source or a sink"
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
