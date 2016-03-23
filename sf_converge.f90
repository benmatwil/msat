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
  do i = 1, nnulls
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
  integer :: i, j, k, count, n, imin, nfw, nbw, nspine, nfan, maxcount
  double precision :: r(3), b(3)
  double precision, dimension(:), allocatable :: roldfw, rnewfw, bnewfw, roldbw, rnewbw, bnewbw
  integer :: flag, testgordon, nfanchk
  double precision :: dphi, dtheta
  double precision :: fact, acc, spinecheck
  double precision :: mindot, dotprod
  integer :: sign, spiral, warning

  double precision, dimension(3) :: spine, fan, maxvec, minvec
  double precision, dimension(:,:), allocatable :: rconvergefw, rconvergebw, rspine, rfan, rfanchk, dummy

  !set up theta and phi for sphere around null
  dphi = 360.d0/dble(nphi)
  dtheta = 180./dble(ntheta-1)

  do i = 1,nphi
    phis(i) = (i-1)*dphi
  enddo

  do j = 1,ntheta
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
      roldbw = [0,0,0] !while limit needs changing to make accurate
      do while (count < maxcount) ! .and. modulus(rnewfw-roldfw) > acc .and. modulus(rnewbw-roldbw) > acc .and. 
        call it_conv(roldfw,rnewfw,bnewfw,fact,1)
        call it_conv(roldbw,rnewbw,bnewbw,fact,-1)
        count = count + 1
        if (modulus(rnewfw-roldfw) < acc .or. modulus(rnewbw-roldbw) < acc) flag = 1
        if (flag == 0) maxcount = 2*count
      enddo
      !print*, count
      n = i+(j-1)*nphi
      rconvergefw(:,n) = normalise(rnewfw)
      rconvergebw(:,n) = normalise(rnewbw)
    enddo
  enddo

  call remove_duplicates(rconvergefw, 1d-4) ! what do we want to do if we are still left with two vectors
  call remove_duplicates(rconvergebw, 1d-4) ! do we want to reduce rfan

  nfw = size(rconvergefw,2)
  nbw = size(rconvergebw,2)

  !with current code, this may fail to pick the correct one
  if (nfw < nbw) then !want a better way of confirming which is fan and which is spine
    ! spine going out of null
    rspine = rconvergefw
    rfan = rconvergebw
    sign = -1
  else if (nfw > nbw) then
    ! spine going into null
    rspine = rconvergebw
    rfan = rconvergefw
    sign = 1
  else if (nfw == nbw) then
    !need a way to check this
    rspine = rconvergebw
    rfan = rconvergefw
  endif

  deallocate(rconvergebw, rconvergefw)

  rfanchk = rfan
  call remove_duplicates(rfanchk,1d-1)

  nfanchk = size(rfanchk,2)
  if (nfanchk == 2) then
    deallocate(rfan)
    rfan = rfanchk
    deallocate(rfanchk)
  endif

  nspine = size(rspine,2)
  nfan = size(rfan,2)
  print*, "Number of spines:", nspine
  print*, "Number of fans:", nfan
  print*, "Number of points left in fan:", nfanchk

  if ((nfan > 2 .and. nfan < 10000) .or. nspine > 2) then
    open(unit=10, file="spinedata.dat", access="stream")
    write(10) size(rspine,2), rspine
    close(10)
    open(unit=10, file="fandata.dat", access="stream")
    write(10) size(rfan,2), rfan
    close(10)
  endif

  spine = rspine(:,1) !which should we pick if > 1 spine vector
  maxvec = rfan(:,1) !pick this more intelligently? theres a whole ring of them

  !what if nfan = 3/4?
  spinecheck = dot(trilinear(rsphere*spine+rnull,bgrid),spine)/rsphere
  if (nfan == 2 .and. abs(dot(trilinear(rsphere*maxvec+rnull,bgrid),maxvec)/rsphere) > abs(spinecheck)) then
    dummy = rspine
    deallocate(rspine)
    rspine = rfan
    deallocate(rfan)
    rfan = dummy
    
    nspine = size(rspine,2)
    nfan = size(rfan,2)
    
    spine = rspine(:,1)
    maxvec = rfan(:,1)
    
    if (spinecheck > 0) then
      sign = -1
    else
      sign = 1
    endif 
  endif

  print*, "Maxvec is:", maxvec

  if (nfan == 2) then
    minvec = cross(spine,maxvec)
  else
    mindot = 1
    do i = 2, size(rfan,2)
      dotprod = abs(dot(maxvec,rfan(:,i)))
      if (dotprod < mindot) then
        mindot = dotprod
        imin = i
      endif
    enddo
    minvec = rfan(:,imin)
    print*, "Perp vec is at ", imin, "out of a total of ", size(rfan)
  endif

  print*, "Minvec is:", minvec

  if (size(rfan,2) == 2) then
    print*,'-------------------------------------------------------------------------'
    print*,dot(trilinear(rsphere*spine+rnull,bgrid),spine)/rsphere
    print*,dot(trilinear(rsphere*maxvec+rnull,bgrid),maxvec)/rsphere
    print*,dot(trilinear(rsphere*minvec+rnull,bgrid),minvec)/rsphere
  endif

  print*, '-------------------------------------------------------------------------'
  print*, "final dot prod is"
  print*, "        min/max           ", "      min/spine          ", "      max/spine       "
  print*, dot(minvec,maxvec), dot(minvec,spine), dot(maxvec,spine)
  print*, ''
  !print*, 'Spine = ', spine
  !print*, 'Fan =   ', cross(minvec,maxvec)
  print*, '-------------------------------------------------------------------------'

  !sign = -sign !check which sign needs to be which
  spiral = 0

  testgordon = 0
  if (testgordon == 1) then
    print*, ''
    print*, "Determining null's properties..."

    !test the estimated properties and adjust them if required. Run this twice to be sure
    call test_null(spine,maxvec,minvec,sign,spiral)
    call test_null(spine,maxvec,minvec,sign,spiral)
  endif

  fan = normalise(cross(minvec,maxvec))

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
