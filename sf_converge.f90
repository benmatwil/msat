!spine finder with convergence method
program sf_converge

  use params
  use sfmod_converge

  implicit none

  double precision, allocatable :: rnulls(:,:), spines(:,:), fans(:,:)
  integer, allocatable :: signs(:), warnings(:)
  
  integer :: nnulls, nullstart, nullend, savedata = 0
  character(len=100) :: arg

  integer :: pcount, ncount, ucount

  integer :: sign, warning
  double precision, dimension(3) :: spine, fan

  integer :: i

  nullstart = 1
  nullend = nnulls
  if (command_argument_count() > 0) then 
    do i = 1, command_argument_count()
      call get_command_argument(i,arg)
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
    do i = 1, command_argument_count()
      call get_command_argument(i,arg)
      if (arg(1:5) == 'data=')
        filename = trim(arg(6:))
      endif
    enddo
  endif

  !Read in 'null.dat'
  open (unit=10,file='output/null.dat',form='unformatted')

  read(10) nnulls
    allocate(rnulls(3,nnulls),spines(3,nnulls),fans(3,nnulls),signs(nnulls),warnings(nnulls))

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

  print*, 'There are ', nnulls,' nulls to analyse'
  print*, '-----------------------------------------------------------------------------'

  !now loop over each null and characterise
  do i = nullstart, nullend!1, nnulls
    print*, 'Evaluating null', i,' of', nnulls
    rnull = rnulls(:,i)
    
    call get_properties(sign,spine,fan,warning,savedata)
    
    signs(i) = sign
    spines(:,i) = spine
    fans(:,i) = fan
    warnings(i) = warning
    
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
  do i = 1, nnulls
    print*, i, signs(i), warnings(i)
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
subroutine get_properties(sign,spine,fan,warning,savedata)
  use sfmod_converge

  implicit none
  integer :: itheta, iphi, itry, ifan, idense
  integer :: count, iconv, imin, nfw, nbw, nspine, nfan, maxcount
  double precision :: r(3), b(3)

  double precision :: dphi, dtheta, angle
  double precision, dimension(:), allocatable :: roldfw, rnewfw, bnewfw, roldbw, rnewbw, bnewbw
  integer :: flagfw, flagbw, savedata, angleflag
  double precision :: fact, acc, fwvalue, bwvalue

  double precision :: mindot, dotprod
  integer :: sign, warning, signguess

  double precision, dimension(3) :: spine, fan, maxvec, minvec
  
  double precision, dimension(:,:), allocatable :: rconvergefw, rconvergebw
  double precision, dimension(:,:), allocatable :: rspine, rfan
  double precision, dimension(:,:), allocatable :: rfanperp, crossfan
  integer, allocatable, dimension(:,:) :: densepos, denseposfw, denseposbw

  print*, 'Null position:'
  print*, rnull
  print*, 'B at null'
  print*, trilinear(rnull, bgrid)
  print*, '-----------------------------------------------------------------------------'

  allocate(rconvergefw(3,nphi*ntheta), rconvergebw(3,nphi*ntheta))
  allocate(roldfw(3), rnewfw(3), bnewfw(3), roldbw(3), rnewbw(3), bnewbw(3))

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
    acc = 1d-10*rsphere
    flagfw = 0
    flagbw = 0
    ! loop over each point to find converged point
    do itheta = 1, ntheta
      do iphi = 1, nphi
        count = 0
        maxcount = 5000
        r = sphere2cart(rsphere,thetas(itheta),phis(iphi))
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
        enddo
        if (count < 100 .or. count == maxcount) then
          angle = angle + 1
          print*, 'Adjusting the initial points and starting again'
          cycle main
        endif
        if (modulus(rnewfw-roldfw) < acc) flagfw = flagfw + 1
        if (modulus(rnewbw-roldbw) < acc) flagbw = flagbw + 1
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
  else if (flagfw > flagbw) then
    signguess = -1
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
    call remove_duplicates(rconvergefw, acc, denseposfw)
    call remove_duplicates(rconvergebw, acc, denseposbw)
    nfw = size(rconvergefw,2)
    nbw = size(rconvergebw,2)
    ! first determine the spine and which is which
    if (nfw == 2 .or. nbw == 2) then ! if one of the reductions only has two vectors, can guarantee one is spine
      if (nfw < nbw) then ! if nfw is smaller then it's the spine, spine going out of null
        sign = -1
      else if (nfw > nbw) then ! if nbw is smaller then it's the spine, spine going into null
        sign = 1
      else if (nfw == nbw) then ! both are equal (i.e. 2=2) so check which converged first and then use eigenvalues as last resort
        if (flagfw < flagbw) then ! backwards converged points stopped first
          sign = 1
        else if (flagfw > flagbw) then ! forwards converged points stopped first
          sign = -1
        else ! still don't know, check the eigenvalues as a last resort - not ideal but rare
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
      else
        rspine = rconvergefw
        rfan = rconvergebw
      endif
      spine = rspine(:,1)
      exit
    else if (itry == 2 .and. sign == 0) then 
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
        sign = 1
      else if (2*nbw < nfw) then
        spine = rconvergefw(:,maxval(maxloc(denseposfw)))
        rfan = rconvergebw
        rspine = rconvergefw
        sign = -1
      else ! not sure how to to decide in this case, haven't found a case like this yet...
        print*, "Uh oh, can't decide, error!"
        open(unit=10, file='spinedata.dat', access='stream')
          write(10) size(rconvergefw,2), rconvergefw
        close(10)
        open(unit=10, file='fandata.dat', access='stream')
          write(10) size(rconvergebw,2), rconvergebw
        close(10)
        sign = 0 ! check to make sure it isn't flagging any proper nulls
        stop
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
    call remove_duplicates(rfan, 1d-2, densepos)
    maxvec = rfan(:,maxval(maxloc(densepos)))
    deallocate(densepos)
  endif
  
  if (sign /= signguess) then
    sign = 0
    print*, 'The two checks do not agree, needs looking into'
  endif
  
  print*, 'Final number of converged'
  print*, '    ','spines:', size(rspine,2)
  print*, '    ','fans:  ', size(rfan,2)
  
  ! now determine the fan
  nfan = size(rfan,2)
  if (nfan == 2) then
    print*, 'Only have two fan points'
    minvec = normalise(cross(spine,maxvec))
  else
    ! Check whether we have a ball or all points in a ring
    crossfan = rfan
    call remove_duplicates(crossfan, 1d-1, densepos)
    idense = maxval(maxloc(densepos))
    deallocate(densepos)
    do ifan = 1, size(crossfan,2)
      crossfan(:,ifan) = cross(maxvec,crossfan(:,ifan))
    enddo
    call remove_vector(crossfan, idense)
    ifan = 1
    do while (ifan <= size(crossfan,2))
      if (modulus(crossfan(:,ifan)) < 2d-1) then
        call remove_vector(crossfan, ifan)
      else
        crossfan(:,ifan) = normalise(crossfan(:,ifan))
        ifan = ifan+1
      endif
    enddo
    call remove_duplicates(crossfan,2d-1)
    print*, 'Size of crossfan', size(crossfan,2)
    if (size(crossfan,2) <= 4) then ! we have a ring, pick minvec to be vector most perpendicular
      print*, 'We have a ring'
      angleflag = 0
      do ifan = 1, size(rfan,2)
        if (acos(dot(rfan(:,ifan), maxvec))/dtor > 30) angleflag = 1
      enddo
      if (angleflag == 1) then
        print*, 'Picking the most perpendicular point'
        mindot = 1
        do ifan = 1, size(rfan,2)
          dotprod = abs(dot(maxvec,rfan(:,ifan)))
          if (dotprod < mindot) then
            mindot = dotprod
            imin = ifan
          endif
        enddo
        minvec = rfan(:,imin)
        print*, 'Perp vec is at ', imin, 'out of a total of ', size(rfan,2)
      else
        print*, 'Taking cross product for minvec'
        minvec = normalise(cross(spine, maxvec))
      endif
    else ! we have a ball
      ifan = 1
      rfanperp = rfan
      do while (ifan <= size(rfanperp,2))
        dotprod = abs(dot(maxvec,rfanperp(:,ifan)))
        if (dotprod > 2d-1) then
          call remove_vector(rfanperp, ifan)
        else
          ifan = ifan + 1
        endif
      enddo
      if (allocated(densepos)) deallocate(densepos)
      call remove_duplicates(rfanperp, 5d-2, densepos)
      minvec = rfanperp(:,maxval(maxloc(densepos)))
    endif
  endif

  ! Save data if there is a problematic null for inspection by map.pro
  if (savedata == 1) then
    open(unit=10, file='spinedata.dat', access='stream')
    write(10) size(rspine,2), rspine, spine
    close(10)
    open(unit=10, file='fandata.dat', access='stream')
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
