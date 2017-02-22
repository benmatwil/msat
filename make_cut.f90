program make_cut

  use params

  implicit none

  real(np) :: r0

  ! integer(int32) :: nx, ny, nz
  ! real(np), allocatable :: bgrid(:,:,:,:)
  ! real(np), dimension(:), allocatable :: x, y, z

  integer(int32) :: nnulls
  ! integer(int32), allocatable :: signs(:)
  ! real(np), allocatable :: spines(:,:)

  integer :: iarg
  character(10) :: arg

  integer(int32) :: inull, iring, iline, nextline, flag
  integer(int32) :: nrings, npoints, nlines
  ! integer(8) :: p, uptoring
  integer(int32), dimension(:), allocatable :: nperring, breaks
  real(np), dimension(:,:), allocatable :: line
  real(np) :: s, point(3)

  integer(int32) :: nspine, dir

  if (command_argument_count() > 0) then
    do iarg = 1, command_argument_count()
      call get_command_argument(iarg,arg)
      if (trim(arg) == '-r') then
        call get_command_argument(iarg+1,arg)
        read(arg,*) r0
      endif
    enddo
  else
    stop "Need to provide a radius"
  endif

  call filenames

  ! open(unit=20, file=filein, access='stream', status='old')
  !   read(20) nx, ny, nz ! number of vertices
  !   allocate(bgrid(nx,ny,nz,3))
  !   allocate(x(nx), y(ny), z(nz))
  !   read(20) bgrid(:,:,:,1)
  !   read(20) bgrid(:,:,:,2)
  !   read(20) bgrid(:,:,:,3)
  !   read(20) x, y, z
  ! close(20)

  open(unit=10, file=trim(fileout)//'-nullpos.dat', access='stream', status='old')
    read(10) nnulls
  close(10)

  ! open(unit=10, file=trim(fileout)//'-nulldata.dat', access='stream', status='old')
  !   read(10) nnulls
  !   allocate(signs(nnulls), spines(3,nnulls))
  !   read(10) signs, spines
  ! close(10)


  ! Need to deal with points at the edge/breakpoints

  open(unit=1, file=trim(fileout)//'-cut_rings.dat', access='stream', status='replace')
  open(unit=10, file=trim(fileout)//'-ringinfo.dat', access='stream', status='old')
  read(10) nrings
  open(unit=20, file=trim(fileout)//'-rings.dat', access='stream', status='old')
  
  do inull = 1, nnulls
    allocate(nperring(ringsmax+1))
    read(10) nperring
    nrings = count(nperring > 0)
    do iring = 1, nrings
      allocate(breaks(nperring(iring)), line(3,nperring(iring)))
      read(20) breaks, line
      do iline = 1, nperring(iring)
        if (iline == nperring(iring)) then
          nextline = 1
        else
          nextline = iline+1
        endif
        if ((line(1,iline)-r0)*(line(1,nextline)-r0) < 0) then
          ! find point in between at r0 and add to line
          s = (r0 - line(1,iline))/(line(1,nextline) - line(1,iline))
          point = line(:,iline) + s*(line(:,nextline) - line(:,iline))
          write(1) point
        endif
      enddo
      deallocate(line, breaks)
    enddo
    deallocate(nperring)
  enddo
  close(1)
  close(10)
  close(20)

  ! for each separator, check whether we cross r=r0 and add point to file
  open(unit=2, file=trim(fileout)//'-cut_seps.dat', access='stream', status='replace')
  open(unit=40, file=trim(fileout)//'-connectivity.dat', access='stream', status='old')
  open(unit=50, file=trim(fileout)//'-separators.dat', access='stream', status='old')
  read(40) flag
  do while (flag > 0)
    read(50) npoints
    allocate(line(3,npoints))
    read(50) line
    do iline = 1, npoints-1
      nextline = iline+1
      if ((line(1,iline)-r0)*(line(1,nextline)-r0) < 0) then
        s = (r0 - line(1,iline))/(line(1,nextline) - line(1,iline))
        point = line(:,iline) + s*(line(:,nextline) - line(:,iline))
        write(2) point
      endif
    enddo
    deallocate(line)
    read(40) flag, flag
  enddo
  close(2)
  close(50)
  
  ! for each spine, check whether we cross r=r0 and add point to file
  open(unit=3, file=trim(fileout)//'-cut_spines.dat', access='stream', status='replace')
  open(unit=30, file=trim(fileout)//'-spines.dat', access='stream', status='old')
  
  ! for each spine, check whether we cross r=r0 and add point to file
  do inull = 1, nnulls
    ! read in spine
    do dir = 1, 2
      read(30) nspine
      allocate(line(3,nspine))
      read(30) line
      do iline = 1, nspine-1
        nextline = iline+1
        if ((line(1,iline)-r0)*(line(1,nextline)-r0) < 0) then
          s = (r0 - line(1,iline))/(line(1,nextline) - line(1,iline))
          point = line(:,iline) + s*(line(:,nextline) - line(:,iline))
          write(3) point
        endif
      enddo
      deallocate(line)
    enddo
  enddo
  close(3)
  close(30)

end