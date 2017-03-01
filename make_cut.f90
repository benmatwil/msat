program make_cut

  use params
  use common

  implicit none

  real(np) :: r0
  integer(int32) :: rc = 0

  ! integer(int32) :: nnulls
  integer(int32) :: nx, ny, nz
  integer(int64) :: g
  real(np), dimension(:), allocatable :: xg, yg, zg
  real(np) :: ds

  integer(int32) :: iarg
  character(10) :: arg

  integer(int32) :: inull, iring, iline, ip, jp, nextline, flag
  integer(int32) :: imin
  integer(int32) :: nrings, npoints, nlines
  integer(int32), dimension(:), allocatable :: nperring, breaks
  real(np), dimension(:), allocatable :: dists
  real(np), dimension(:,:), allocatable :: line, points, pordered
  real(np) :: s, point(3), diff(3), mindist

  integer(int32) :: nspine, dir

  if (command_argument_count() > 0) then
    do iarg = 1, command_argument_count()
      call get_command_argument(iarg,arg)
      if (trim(arg) == '-r') then
        call get_command_argument(iarg+1,arg)
        read(arg,*) r0
        rc = 1
      endif
    enddo
  endif
  if (rc == 0) then    
    stop "Need to provide a radius"
  endif

  call filenames

  open(unit=10, file=trim(fileout)//'-nullpos.dat', access='stream', status='old')
    read(10) nnulls
  close(10)

  open(unit=10,file=filein,access='stream',status='old')
    read(10), nx, ny, nz !number of vertices
    allocate(xg(nx), yg(ny), zg(nz))
    g = 3_int64*4_int64 + 3_int64*int(nx, int64)*int(ny, int64)*int(nz, int64)*8_int64 + 1_int64
    read(10, pos=g) xg, yg, zg
  close(10)

  ds = maxval([maxval(yg(2:ny) - yg(1:ny-1)), maxval(zg(2:nz) - zg(1:nz-1))])

  ! Need to deal with points at the edge/breakpoints

  open(unit=1, file=trim(fileout)//'-cut_rings.dat', access='stream', status='replace')
  open(unit=10, file=trim(fileout)//'-ringinfo.dat', access='stream', status='old')
  read(10) nrings
  open(unit=20, file=trim(fileout)//'-rings.dat', access='stream', status='old')

  do inull = 1, nnulls
    allocate(points(3,0))
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
        if ((line(1,iline)-r0)*(line(1,nextline)-r0) < 0 .and. breaks(iline) == 0) then
          ! find point in between at r0 and add to line
          s = (r0 - line(1,iline))/(line(1,nextline) - line(1,iline))
          point = line(:,iline) + s*(line(:,nextline) - line(:,iline))
          call add_vector(points, point)
        endif
      enddo
      deallocate(line, breaks)
    enddo
    deallocate(nperring)


    do while (size(points, 2) > 0)
      allocate(pordered(3, 1))
      pordered(:, 1) = points(:, 1)
      call remove_vector(points, 1)
      ! add points in the "forward" direction
      ip = 1
      do while (size(points, 2) > 0)
        allocate(dists(size(points, 2)))
        do jp = 1, size(points, 2)
          dists(jp) = sum((points(:,jp) - pordered(:, ip))**2)
        enddo
        imin = minloc(dists, 1)
        mindist = dists(imin)
        deallocate(dists)
        if (mindist < ds/2) then
          call add_vector(pordered, points(:, imin))
          call remove_vector(points, imin)
        else
          exit
        endif
        ip = ip + 1
      enddo
      ! add points in the "backward" direction
      do while (size(points, 2) > 0)
        allocate(dists(size(points, 2)))
        do jp = 1, size(points, 2)
          dists(jp) = sum((points(:, jp) - pordered(:, 1))**2)
        enddo
        imin = minloc(dists, 1)
        mindist = dists(imin)
        deallocate(dists)
        if (mindist < ds/2) then
          call add_vector(pordered, points(:, imin), 1)
          call remove_vector(points, imin)
        else
          exit
        endif
      enddo
      write(1) size(pordered, 2), pordered
      deallocate(pordered)
    enddo
    deallocate(points)
  enddo
  write(1) -1
  close(1)
  close(10)
  close(20)

  ! for hcs
  open(unit=10, file=trim(fileout)//'-hcs.dat', access='stream', status='old')
  open(unit=5, file=trim(fileout)//'-cut_hcs.dat', access='stream', status='replace')
  read(10) npoints
  write(10, pos=(1+2*npoints*3)*4+1)
  allocate(points(3, 0))
  do ip = 1, 2*npoints
    read(10) npoints
    allocate(line(3, npoints))
    read(10) line
    do iline = 1, npoints-1
      nextline = iline+1
      if ((line(1,iline)-r0)*(line(1,nextline)-r0) < 0) then
        s = (r0 - line(1,iline))/(line(1,nextline) - line(1,iline))
        point = line(:,iline) + s*(line(:,nextline) - line(:,iline))
        call add_vector(points, point)
      endif
    enddo
    deallocate(line)
  enddo
  print*, size(points, 2)

  do while (size(points, 2) > 0)
    allocate(pordered(3, 1))
    pordered(:, 1) = points(:, 1)
    call remove_vector(points, 1)
    ! add points in the "forward" direction
    ip = 1
    do while (size(points, 2) > 0)
      allocate(dists(size(points, 2)))
      do jp = 1, size(points, 2)
        dists(jp) = sum((points(:,jp) - pordered(:, ip))**2)
      enddo
      imin = minloc(dists, 1)
      mindist = dists(imin)
      deallocate(dists)
      if (mindist < ds/2) then
        call add_vector(pordered, points(:, imin))
        call remove_vector(points, imin)
      else
        exit
      endif
      ip = ip + 1
    enddo
    ! add points in the "backward" direction
    do while (size(points, 2) > 0)
      allocate(dists(size(points, 2)))
      do jp = 1, size(points, 2)
        dists(jp) = sum((points(:, jp) - pordered(:, 1))**2)
      enddo
      imin = minloc(dists, 1)
      mindist = dists(imin)
      deallocate(dists)
      if (mindist < ds/2) then
        call add_vector(pordered, points(:, imin), 1)
        call remove_vector(points, imin)
      else
        exit
      endif
    enddo
    write(5) size(pordered, 2), pordered
    print*, size(pordered, 2)
    deallocate(pordered)
  enddo
  deallocate(points)
  write(5) -1
  close(5)
  close(10)

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
        write(2) flag, point
      endif
    enddo
    deallocate(line)
    read(40) flag, flag
  enddo
  write(2) -1
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
          write(3) inull, point
        endif
      enddo
      deallocate(line)
    enddo
  enddo
  write(3) -1
  close(3)
  close(30)

end