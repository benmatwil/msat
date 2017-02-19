program make_cut

  use params
  use trace

  implicit none

  real(np) :: r0 = 1.2_np
  character(4) :: fname

  integer :: nx, ny, nz
  real(np), allocatable :: bgrid(:,:,:,:)
  real(np), dimension(:), allocatable :: x, y, z

  integer :: nnulls
  integer, allocatable :: signs(:)
  real(np), allocatable :: spines(:,:)

  integer :: iring, iline, nextline
  integer :: nrings, npoints, nringsmax, nlines
  integer(8) :: p, uptoring
  integer, dimension(:), allocatable :: nperring
  real(np), dimension(:,:), allocatable :: line1
  real(np) :: s, point(3)

  integer :: cont, inull
  real(np), dimension(:), allocatable :: x, y, z

  integer :: dir, sign
  real(np) :: rstart(:)

  call filenames

  open(unit=20, file=filein, access='stream', status='old')
    read(20) nx, ny, nz ! number of vertices
    allocate(bgrid(nx,ny,nz,3))
    allocate(x(nx), y(ny), z(nz))
    read(20) bgrid(:,:,:,1)
    read(20) bgrid(:,:,:,2)
    read(20) bgrid(:,:,:,3)
    read(20) x, y, z
  close(20)

  open(unit=10, file=trim(fileout)//'-nullpos.dat', access='stream', status='old')
    read(10) nnulls
  close(10)

  open(unit=10, file=trim(fileout)//'-nulldata.dat', access='stream', status='old')
    read(10) nnulls
    allocate(signs(nnulls), spines(3,nnulls))
    read(10) signs, spines
  close(10)


  ! Need to deal with points at the edge/breakpoints

  open(unit=12,file='output/cut_rings.dat',access='stream',status='replace')
  do inull = 1, nnulls
  print*, inull
    write(fname,'(I4.4)') inull
    open(unit=10,file=trim(fileout)//'-everything'//trim(fname)//'.dat',access='stream',status='old')
      read(10) nrings, npoints, nringsmax
      ! print*, nrings
      allocate(nperring(0:nringsmax-1))
      read(10) nperring
      ! print*, nperring(nrings-1:nrings+1)
    
    iring = 1
    do iring = 0, nrings-1
      allocate(line1(3,nperring(iring)))
      uptoring = 3 + nringsmax + sum(nperring(0:iring-1))*8
      p = uptoring + nperring(iring)*2
      p = p*4 + 1
      read(10, pos=p) line1
      do iline = 1, nperring(iring)
        if (iline == nperring(iring)) then
          nextline = 1
        else
          nextline = iline+1
        endif
        if ((line1(1,iline)-r0)*(line1(1,nextline)-r0) < 0 .and. abs(line1(1,iline) - line1(1,nextline)) < 1) then
          ! find point in between at r0 and add to line
          s = (r0 - line1(1,iline))/(line1(1,nextline) - line1(1,iline))
          point = line1(:,iline) + s*(line1(:,nextline) - line1(:,iline))
          write(12) point
          ! print*, point(1), line1(1,iline), line1(1,nextline)
        endif
      enddo
      deallocate(line1)
    enddo
    deallocate(nperring)
    close(10)
  enddo
  close(12)

  ! for each separator, check whether we cross r=r0 and add point to file
  open(unit=13, file='output/cut_seps.dat',access='stream',status='replace')
  do inull = 1, nnulls
    write(fname,'(I4.4)') inull
    open(unit=11, file=trim(fileout)//'-sep'//trim(fname)//'.dat',access='stream',status='old')
      read(11) cont
      
    do while (cont > 0)
      read(11) nringsmax
      allocate(line1(cont,3))
      read(11) line1
      do iline = 1, cont
        if (iline == cont) then
          nextline = 1
        else
          nextline = iline+1
        endif
        if ((line1(iline,1)-r0)*(line1(nextline,1)-r0) < 0) then
          ! find point in between two points (function?) and save
          s = (r0 - line1(iline,1))/(line1(nextline,1) - line1(iline,1))
          point = line1(iline,:) + s*(line1(nextline,:) - line1(iline,:))
          write(12) point
        endif
      enddo
      deallocate(line1)
      read(11) cont
    enddo

    close(11)
  enddo
  close(13)
  
  
  ! ! for each spine, check whether we cross r=r0 and add point to file
  ! do inull = 1, nnulls
  !   ! read in spine
  !   read(12) spine
  !   do dir = -1, 1, 2
  !     ! Trace out spine for each null and check for crossing the plane
  !     rstart = rnull(inull,:) + dir*spine*1e-3_np
  !     h = 0.1_np
  !     call trace_line(r, sign, h)
  !     ! if crosses, save points, interpolate
  !   enddo
  ! enddo

end