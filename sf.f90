!spine finder with convergence method
program sf_converge

  use params
  use sf_converge_mod

  implicit none

  double precision, allocatable :: rnulls(:,:), spines(:,:), fans(:,:)
  integer, allocatable :: signs(:), warnings(:)
  
  integer :: nnulls, nullstart, nullend, savedata = 0

  integer :: pcount, ncount, ucount

  integer :: sign, warning
  double precision, dimension(3) :: spine, fan

  integer :: inull
  
  print*,'#######################################################################'
  print*,'#                             Spinefinder                             #'
  print*,'#######################################################################'

  !Read in 'null.dat'
  open (unit=10,file='output/null.dat',form='unformatted')

  read(10) nnulls
    allocate(rnulls(3,nnulls),spines(3,nnulls),fans(3,nnulls),signs(nnulls),warnings(nnulls))

    do inull = 1,nnulls
      read(10) rnulls(:,inull)
    enddo
  close(10)

  nullstart = 1
  nullend = nnulls
  if (command_argument_count() > 0) then 
    do icommand = 1, command_argument_count()
      call get_command_argument(icommand,arg)
      if (arg(1:2) == 'n=') then
        arg = arg(3:)
        read(arg,*) nullstart
        nullend = nullstart
        savedata = 1
      endif
    enddo
  endif

  filename = defaultfilename
  if (command_argument_count() > 0) then
    do icommand = 1, command_argument_count()
      call get_command_argument(icommand,arg)
      if (arg(1:5) == 'data=') then
        filename = trim(arg(6:))
      endif
    enddo
  endif

  !read in bgrid
  open(unit=10,file=filename,access='stream')
    read(10) nx, ny, nz
    allocate(bgrid(nx,ny,nz,3), x(nx), y(ny), z(nz))
    read(10) bgrid
    read(10) x, y, z
  close(10)

  print*, 'There are ', nnulls,' nulls to analyse'
  print*, '-----------------------------------------------------------------------------'

  !now loop over each null and characterise
  do inull = nullstart, nullend!1, nnulls
    print*, 'Evaluating null', inull,' of', nnulls
    rnull = rnulls(:,inull)
    
    call get_properties(sign,spine,fan,warning,savedata)
    
    signs(inull) = sign
    spines(:,inull) = spine
    fans(:,inull) = fan
    warnings(inull) = warning
    
    print*, ''
    print*, '-----------------------------------------------------------------------------'
    print*, '-----------------------------------------------------------------------------'
    print*, ''
  enddo
  if (nullend - nullstart /= nnulls - 1) stop
  
  !now write data to nulls.dat
  open(unit=10,file='output/nulls.dat',form='unformatted')
    write(10) nnulls
    write(10) signs
    write(10) rnulls,spines,fans
  close(10)

  print*, 'Summary:'
  print*, 'number, sign, warning'

  pcount = 0
  ucount = 0
  ncount = 0
  do inull = 1, nnulls
    print*, inull, signs(inull), warnings(inull)
    if (signs(inull) .eq. 1) pcount = pcount+1
    if (signs(inull) .eq. -1) ncount = ncount+1
    if (signs(inull) .eq. 0) ucount = ucount+1
  enddo

  print*, 'Total number of nulls:', nnulls
  print*, 'Positive', pcount
  print*, 'Negative', ncount
  print*, 'Unknown', ucount
  print*, 'Warning', sum(warnings)
  print*, 'All nulls have data, should be equal', pcount+ncount+ucount, nnulls
  print*, 'Done!'

end program

!********************************************************************************

!characterise each null
subroutine get_properties(sign,spine,fan,warning,savedata)
  use sf_converge_mod

  implicit none
  integer :: itheta, iphi, itry, ifan, idense
  integer :: count, iconv, imin, nfw, nbw, maxcount
  integer :: flagfw, flagbw, savedata, angleflag

  double precision :: dphi, dtheta, angle
  double precision, dimension(3) :: roldfw, rnewfw, bnewfw, roldbw, rnewbw, bnewbw
  double precision, dimension(3) :: rold, rnew, bnew
  double precision :: fact, acc, accconv, fwvalue, bwvalue

  double precision, dimension(:,:), allocatable :: rconvergefw, rconvergebw, rmin
  double precision, dimension(:,:), allocatable :: rspine, rfan
  double precision, dimension(:,:), allocatable :: rfanperp, crossfan
  integer, dimension(:,:), allocatable :: densepos, denseposfw, denseposbw

  double precision :: mindot, dotprod
  integer :: sign, warning, signguess

  double precision, dimension(3) :: spine, fan, maxvec, minvec

  print*, 'Null position:'
  print*, rnull
  print*, 'B at null'
  print*, trilinear(rnull, bgrid)
  print*, '-----------------------------------------------------------------------------'

  allocate(rconvergefw(3,nphi*ntheta), rconvergebw(3,nphi*ntheta), rmin(3,nphi*ntheta))

  !set up theta and phi for sphere around null
  dphi = 360.d0/dble(nphi)
  dtheta = 180.d0/dble(ntheta-1)
  
  angle = 0
  main: do
    do iphi = 1, nphi
      phis(iphi) = (iphi-1)*dphi + angle
    enddo
    do itheta = 1, ntheta
      thetas(itheta) = (itheta-1)*dtheta + angle
    enddo

    phis = phis*dtor
    thetas = thetas*dtor
    
    ! Working on a sphere of size rsphere
    fact = 1d-2*rsphere
    accconv = 1d-10*rsphere
    flagfw = 0
    flagbw = 0
    maxcount = 5000
    ! loop over each point to find converged point
    do itheta = 1, ntheta
      do iphi = 1, nphi
        count = 0
        rnewfw = sphere2cart(rsphere,thetas(itheta),phis(iphi))
        rnewbw = rnewfw
        bnewfw = trilinear(rnewfw+rnull, bgrid)
        bnewbw = bnewfw
        roldfw = [0,0,0]
        roldbw = [0,0,0]
        do while (modulus(rnewfw-roldfw) > accconv .and. modulus(rnewbw-roldbw) > accconv .and. count < maxcount) ! do enough times for spine to converge
          call it_conv(roldfw,rnewfw,bnewfw,fact,1)
          call it_conv(roldbw,rnewbw,bnewbw,fact,-1)
          count = count + 1
        enddo
        if (count < 100 .and. angle < 4d0) then
          angle = angle + 1d0
          print*, 'Adjusting the initial points and starting again'
          cycle main
        endif
        if (modulus(rnewfw-roldfw) < accconv) flagfw = flagfw + 1
        if (modulus(rnewbw-roldbw) < accconv) flagbw = flagbw + 1
        iconv = iphi+(itheta-1)*nphi
        rconvergefw(:,iconv) = normalise(rnewfw)
        rconvergebw(:,iconv) = normalise(rnewbw)
      enddo
    enddo
    print*, '-----------------------------------------------------------------------------'
    exit
  enddo main

  print*, 'Number of points which stop convergence'
  print*, '    ', 'forwards: ', flagfw
  print*, '    ', 'backwards:', flagbw
  if (flagfw < flagbw) then
    signguess = 1
    print*, "Spine is backwards"
  else if (flagfw > flagbw) then
    signguess = -1
    print*, "Spine is forwards"
  else
    signguess = 0
  endif
  
  ! Now working on a unit sphere of points
  ! remove duplicate vectors which are a distance 1d-4 apart
  ! spine should be left with 2 vectors
  ! fan left with either 2 vectors (minvec eigenvalue too small), a circle of points or just a mess
  sign = 0
  do itry = 4, 2, -1
    acc = 1d1**(-itry)
    print*, 'Removing to accuracy', acc
    call remove_duplicates(rconvergefw, acc, denseposfw)
    call remove_duplicates(rconvergebw, acc, denseposbw)
    nfw = size(rconvergefw,2)
    nbw = size(rconvergebw,2)
    ! first determine the spine and which is which
    if (nfw == 2 .or. nbw == 2) then ! if one of the reductions only has two vectors, can guarantee one is spine
      print*, "Using converged to choose"
      if (nfw < nbw) then ! if nfw is smaller then it's the spine, spine going out of null
        sign = -1
      else if (nfw > nbw) then ! if nbw is smaller then it's the spine, spine going into null
        sign = 1
      else if (nfw == nbw) then ! both are equal (i.e. 2=2) so check which converged first and then use eigenvalues as last resort
        print*, "Using stopped to choose"
        if (flagfw < flagbw) then ! backwards converged points stopped first
          sign = 1
        else if (flagfw > flagbw) then ! forwards converged points stopped first
          sign = -1
        else ! still don't know, check the eigenvalues as a last resort - not ideal but rare
          print*, "Using eigenvalues to choose"
          bwvalue = dot(trilinear(rsphere*rconvergebw(:,1) + rnull, bgrid), rconvergebw(:,1))
          fwvalue = dot(trilinear(rsphere*rconvergefw(:,1) + rnull, bgrid), rconvergefw(:,1))
          print*, 'Eigenvalue on backward vector:', bwvalue
          print*, 'Eigenvalue on forward vector: ', fwvalue
          if (abs(bwvalue) > abs(fwvalue)) then
            sign = 1
          else
            sign = -1
          endif
        endif
      endif
      if (sign == 1) then
        rspine = rconvergebw
        rfan = rconvergefw
        densepos = denseposfw
      else
        rspine = rconvergefw
        rfan = rconvergebw
        densepos = denseposbw
      endif
      spine = rspine(:,1)
      exit
    else if (itry == 2 .and. sign == 0) then 
      if (nfw > 50*nbw) then
        spine = rconvergebw(:,maxval(maxloc(denseposbw)))
        rspine = rconvergebw
        rfan = rconvergefw
        densepos = denseposfw
        sign = 1
      elseif (nbw > 50*nfw) then
        spine = rconvergefw(:,maxval(maxloc(denseposfw)))
        rspine = rconvergefw
        rfan = rconvergebw
        densepos = denseposbw
        sign = -1
      else
        ! if still no reduction leaves two
        ! try using that the spine should be the densest accumulation of points
        ! also likely we have a source or a sink so div(B) is non-zero
        nfw = maxval(denseposfw)
        nbw = maxval(denseposbw)
        print*, 'Density of densest points'
        print*, '    ', 'Forwards ', maxloc(denseposfw), nfw
        print*, '    ', 'Backwards', maxloc(denseposbw), nbw
        if (2*nfw < nbw) then
          spine = rconvergebw(:,maxval(maxloc(denseposbw)))
          rfan = rconvergefw
          rspine = rconvergebw
          densepos = denseposfw
          sign = 1
          goto 100
        else if (2*nbw < nfw) then
          spine = rconvergefw(:,maxval(maxloc(denseposfw)))
          rfan = rconvergebw
          rspine = rconvergefw
          densepos = denseposbw
          sign = -1
          goto 100
        else ! not sure how to to decide in this case, haven't found a case like this yet...
          print*, "Uh oh, can't decide, error!"
          open(unit=10, file='output/spinedata.dat', access='stream')
            write(10) size(rconvergefw,2), rconvergefw
          close(10)
          open(unit=10, file='output/fandata.dat', access='stream')
            write(10) size(rconvergebw,2), rconvergebw
          close(10)
          sign = 0 ! check to make sure it isn't flagging any proper nulls
          stop
        endif
      endif
    endif
  deallocate(denseposfw, denseposbw)
  enddo
  ! We have picked rspine and rfan so can get rid of rconverges
  deallocate(rconvergebw, rconvergefw)
  print*, 'Picked a spine'

  if (size(rfan,2) == 2) then
    maxvec = rfan(:,1)
  else
    !call remove_duplicates(rfan, 1d-2, densepos)
    maxvec = rfan(:,maxval(maxloc(densepos)))
  endif
  deallocate(densepos)
  
  if (sign /= signguess) then
    sign = 0
    print*, 'The two checks do not agree, needs looking into'
    !stop
  endif
  
  print*, 'Final number of converged'
  print*, '    ','spines:', size(rspine,2)
  print*, '    ','fans:  ', size(rfan,2)
 
!   ! now determine the fan 
!   if (size(rfan,2) == 2) then
!     print*, 'Only have two fan points'
!     minvec = normalise(cross(spine,maxvec))
!   else
!     ! Check whether we have a ball or all points in a ring
!     crossfan = rfan
!     call remove_duplicates(crossfan, 1d-1, densepos)
!     idense = maxval(maxloc(densepos))
!     deallocate(densepos)
!     do ifan = 1, size(crossfan,2)
!       crossfan(:,ifan) = cross(maxvec,crossfan(:,ifan))
!     enddo
!     call remove_vector(crossfan, idense)
!     ifan = 1
!     do while (ifan <= size(crossfan,2))
!       if (modulus(crossfan(:,ifan)) < 2d-1) then
!         call remove_vector(crossfan, ifan)
!       else
!         crossfan(:,ifan) = normalise(crossfan(:,ifan))
!         ifan = ifan+1
!       endif
!     enddo
!     call remove_duplicates(crossfan,2d-1)
!     print*, 'Size of crossfan', size(crossfan,2)
!     if (size(crossfan,2) <= 4) then ! we have a ring, pick minvec to be vector most perpendicular
!       print*, 'We have a ring'
!       angleflag = 0
!       do ifan = 1, size(rfan,2)
!         if (acos(dot(rfan(:,ifan), maxvec))/dtor > 30) then
!           print*, acos(dot(rfan(:,ifan), maxvec))/dtor
!           angleflag = 1
!         endif
!       enddo
!       if (angleflag == 1) then
!         print*, 'Picking the most perpendicular point'
!         mindot = 1
!         do ifan = 1, size(rfan,2)
!           dotprod = abs(dot(maxvec,rfan(:,ifan)))
!           if (dotprod < mindot) then
!             mindot = dotprod
!             imin = ifan
!           endif
!         enddo
!         minvec = rfan(:,imin)
!         print*, 'Perp vec is at ', imin, 'out of a total of ', size(rfan,2)
!       else
!         print*, 'Taking cross product for minvec'
!         minvec = normalise(cross(spine, maxvec))
!       endif
!     else ! we have a ball
!       ifan = 1
!       rfanperp = rfan
!       do while (ifan <= size(rfanperp,2))
!         dotprod = abs(dot(maxvec,rfanperp(:,ifan)))
!         if (dotprod > 2d-1) then
!           call remove_vector(rfanperp, ifan)
!         else
!           ifan = ifan + 1
!         endif
!       enddo
!       if (allocated(densepos)) deallocate(densepos)
!       call remove_duplicates(rfanperp, 5d-2, densepos)
!       minvec = rfanperp(:,maxval(maxloc(densepos)))
!     endif
!   endif
! endif

  ! Working on a sphere of size rsphere
  accconv = accconv*2
  ! loop over each point to find converged point
  do itheta = 1, ntheta
    do iphi = 1, nphi
      count = 0
      rnew = sphere2cart(rsphere,thetas(itheta),phis(iphi))
      bnew = trilinear(rnew+rnull, bgrid)
      rold = [0,0,0]
      do while (modulus(rnew-rold) > accconv .and. count < maxcount) ! do enough times for spine to converge
        bnew = bnew - dot(bnew, normalise(maxvec))*normalise(maxvec)
        call it_conv(rold,rnew,bnew,fact,sign)
        count = count + 1
      enddo
      rmin(:,iphi+(itheta-1)*nphi) = normalise(rnew)
    enddo
  enddo
  call remove_duplicates(rmin, acc, densepos)
  minvec = rmin(:,maxval(maxloc(densepos)))

  ! Save data if there is a problematic null for inspection by map.pro
100  if (savedata == 1) then
    open(unit=10, file='output/spinedata.dat', access='stream')
    write(10) size(rspine,2), rspine, spine
    close(10)
    open(unit=10, file='output/fandata.dat', access='stream')
    write(10) size(rfan,2), rfan, maxvec, minvec
    if (allocated(crossfan)) write(10) size(crossfan,2), crossfan
    close(10)
  endif
  !stop

  print*, 'Maxvec is:'
  print*, maxvec
  print*, 'Minvec is:'
  print*, minvec

  print*, '-----------------------------------------------------------------------------'
  print*, 'Final dot products are'
  print*, '        min/max           ', '        max/spine       ', '        min/spine          '
  print*, dot(minvec,maxvec), dot(maxvec,spine), dot(minvec,spine)
  print*, '-----------------------------------------------------------------------------'

  fan = normalise(cross(minvec,maxvec)) ! fan vector is perp to fan plane

  if (sign .eq. 0) print*, 'Warning, sign = 0 and null likely a source or a sink'

  print*, ''
  print*, 'Final information for null'
  print*, '--------------------------'
  print*, 'Sign:'
  print*, sign
  print*, 'Spine:'
  print*, spine
  print*, 'Fan:'
  print*, fan
  print*, 'Tilt:'
  print*, abs(90-acos(dot(fan,spine))/dtor)

  warning = 0
  if (abs(90-acos(dot(fan,spine))/dtor) < 10.) then
    print*, 'WARNING: SPINE AND FAN STRONGLY INCLINED'
    warning = 1
  endif

end subroutine
