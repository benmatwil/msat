module make_cut_mod

  use params
  use common

  implicit none

  real(np) :: ds

  contains

  subroutine sortpoints(points, inull, unit)

    real(np), dimension(:,:), allocatable :: dists, spinedists, dots, spinedots
    real(np), dimension(:,:), allocatable :: points, pordered
    real(np) :: disttol, diff(3), mindist, mindistspine

    integer(int32), dimension(2) :: imindist, iminspinedist, imindot, iminspinedot
    integer(int32) :: ip, jp, i, npoints, unit, inull

    disttol = ds/3
    do while (size(points, 2) > 0)

      if (size(points, 2) > 1000) then
        npoints = 1000
      else
        npoints = size(points, 2)
      endif
      ! print*, size(points, 2)
      allocate(dists(npoints, npoints))
      dists = 1000 ! maybe needs improving
      do ip = 1, npoints
        do jp = ip+1, npoints
          dists(ip, jp) = dist(points(:, ip), points(:, jp))
        enddo
      enddo
      imindist = minloc(dists)
      deallocate(dists)
      allocate(pordered(3, 2))
      pordered(:, 1) = points(:, imindist(1))
      pordered(:, 2) = points(:, imindist(2))
      call remove_vector(points, imindist(1))
      if (imindist(1) < imindist(2)) then
        call remove_vector(points, imindist(2)-1)
      else
        call remove_vector(points, imindist(2))
      endif

      do i = 1, 2
        ! add points in the each direction
        do while (size(points, 2) > 0)
          allocate(dists(size(points, 2), 1))
          allocate(spinedists(size(spines, 2), 1))
          ! allocate(dots(size(points, 2), 1))
          ! allocate(spinedots(size(spines, 2), 1))
          diff = normalise(pordered(:, 1) - pordered(:, 2))
          do ip = 1, size(points, 2)
            dists(ip, 1) = dist(points(:, ip), pordered(:, 1))
            ! dots(ip, 1) = dot(diff, normalise(points(:, ip) - pordered(:, 1)))
          enddo
          do ip = 1, size(spines, 2)
            spinedists(ip, 1) = dist(spines(:, ip), pordered(:, 1))
            ! spinedots(ip, 1) = dot(diff, normalise(spines(:, ip) - pordered(:, 1)))
          enddo
          imindist = minloc(dists)
          iminspinedist = minloc(spinedists)
          mindist = dists(imindist(1), 1)
          mindistspine = spinedists(iminspinedist(1), 1)
          deallocate(dists, spinedists)
          ! deallocate(dots, spinedots)
          if (mindist < mindistspine) then
            if (dot(diff, normalise(points(:, imindist(1)) - pordered(:, 1))) < -0.98) then
              exit
            elseif (mindist < disttol) then
              call add_vector(pordered, points(:, imindist(1)), 1)
              call remove_vector(points, imindist(1))
            elseif (dot(diff, normalise(points(:, imindist(1)) - pordered(:, 1))) > 0.95 .and. &
              mindist < 3*disttol) then
              call add_vector(pordered, points(:, imindist(1)), 1)
              call remove_vector(points, imindist(1))
            else
              exit
            endif
          else
            exit
          endif
        enddo

        pordered = pordered(:, size(pordered, 2):1:-1)
      enddo

      write(unit) inull, size(pordered, 2), pordered
      deallocate(pordered)
    enddo

  end subroutine

end module

program make_cut

  use params
  use common
  use make_cut_mod

  implicit none

  real(np) :: r0
  integer(int32) :: rc = 0
  character(6) :: rstr

  ! integer(int32) :: nnulls
  integer(int32) :: nx, ny, nz
  integer(int64) :: ig
  real(np), dimension(:), allocatable :: xg, yg, zg

  integer(int32) :: iarg
  character(10) :: arg

  integer(int64) :: ia, ip, ib, uptonull
  integer(int32) :: inull, jnull, iring, iline, nextline, flag, ihcs
  integer(int32) :: nrings, npoints, nlines
  integer(int32), dimension(:), allocatable :: breaks, association
  integer(int32), dimension(0:ringsmax) :: nperring
  real(np), dimension(:,:), allocatable :: line, line2, points
  real(np) :: s, point(3)

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

  write(rstr,'(F6.4)') r0

  call filenames

  open(unit=10, file=trim(fileout)//'-nullpos.dat', access='stream', status='old')
    read(10) nnulls
  close(10)

  open(unit=10,file=filein,access='stream',status='old')
    read(10), nx, ny, nz !number of vertices
    allocate(xg(nx), yg(ny), zg(nz))
    ig = 3_int64*4_int64 + 3_int64*int(nx, int64)*int(ny, int64)*int(nz, int64)*8_int64 + 1_int64
    read(10, pos=ig) xg, yg, zg
  close(10)

  !********************************************************************************

  ! for each separator, check whether we cross r=r0 and add point to file
  open(unit=2, file=trim(fileout)//'-separators-cut_'//rstr//'.dat', access='stream', status='replace')
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
    read(40) flag, flag, flag, flag
  enddo
  write(2) -1
  close(2)
  close(50)

  !********************************************************************************

  ! for each separator, check whether we cross r=r0 and add point to file
  open(unit=2, file=trim(fileout)//'-hcs-separators-cut_'//rstr//'.dat', access='stream', status='replace')
  open(unit=40, file=trim(fileout)//'-hcs-connectivity.dat', access='stream', status='old')
  open(unit=50, file=trim(fileout)//'-hcs-separators.dat', access='stream', status='old')
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
    read(40) flag, flag, flag, flag
  enddo
  write(2) -1
  close(2)
  close(50)

  !********************************************************************************
  
  ! for each spine, check whether we cross r=r0 and add point to file
  open(unit=3, file=trim(fileout)//'-spines-cut_'//rstr//'.dat', access='stream', status='replace')
  open(unit=30, file=trim(fileout)//'-spines.dat', access='stream', status='old')

  allocate(spines(3, 0))
  
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
          call add_vector(spines, point)
        endif
      enddo
      deallocate(line)
    enddo
  enddo
  write(3) -1
  close(3)
  close(30)

  !********************************************************************************

  ds = maxval([maxval(yg(2:ny) - yg(1:ny-1)), maxval(zg(2:nz) - zg(1:nz-1))])

  open(unit=1, file=trim(fileout)//'-rings-cut_'//rstr//'.dat', access='stream', status='replace')
  open(unit=10, file=trim(fileout)//'-ringinfo.dat', access='stream', status='old')
  read(10) nrings
  open(unit=20, file=trim(fileout)//'-rings.dat', access='stream', status='old')

  uptonull = 0
  do inull = 1, nnulls
    allocate(points(3,0))
    read(10) nperring
    nrings = count(nperring > 0) - 1
    do iring = nrings, 1, -1
    
      call file_position(iring, nperring, 1, ia, ib, ip)
      allocate(association(nperring(iring)), breaks(nperring(iring)), line(3,nperring(iring)))
      read(20, pos=uptonull+ia) association, breaks, line

      call file_position(iring-1, nperring, 1, ia, ib, ip)
      allocate(line2(3,nperring(iring-1)))
      read(20, pos=uptonull+ip) line2

      do iline = 1, nperring(iring)
        if ((line(1,iline)-r0)*(line2(1,association(iline))-r0) < 0) then
          if (abs(line(3,iline) - line2(3,association(iline))) < 1) then
            ! find point in between at r0 and add to line
            s = (r0 - line(1,iline))/(line2(1,association(iline)) - line(1,iline))
            point = line(:,iline) + s*(line2(:,association(iline)) - line(:,iline))
            call add_vector(points, point)
          endif
        endif
      enddo
      deallocate(line, line2, breaks, association)
    enddo
    uptonull = uptonull + sum(int(nperring, int64))*32_int64

    ! print*, size(points, 2)
    call sortpoints(points, inull, 1)
    
    deallocate(points)
  enddo
  write(1) -1
  close(1)
  close(10)
  close(20)

  !********************************************************************************

  ! for hcs
  open(unit=10, file=trim(fileout)//'-hcs-ringinfo.dat', access='stream', status='old')
  open(unit=20, file=trim(fileout)//'-hcs-rings.dat', access='stream', status='old')
  open(unit=5, file=trim(fileout)//'-hcs-cut_'//rstr//'.dat', access='stream', status='replace')
  read(10) nrings
  
  uptonull = 0
  ihcs = 1
  do dir = 1, 2
    allocate(points(3,0))
    read(10) nperring
    nrings = count(nperring > 0) - 1
    do iring = nrings, 1, -1
    
      call file_position(iring, nperring, 1, ia, ib, ip)
      allocate(association(nperring(iring)), breaks(nperring(iring)), line(3,nperring(iring)))
      read(20, pos=uptonull+ia) association, breaks, line

      call file_position(iring-1, nperring, 1, ia, ib, ip)
      allocate(line2(3,nperring(iring-1)))
      read(20, pos=uptonull+ip) line2

      do iline = 1, nperring(iring)
        if ((line(1,iline)-r0)*(line2(1,association(iline))-r0) < 0) then
          if (abs(line(3,iline) - line2(3,association(iline))) < 1) then
            ! find point in between at r0 and add to line
            s = (r0 - line(1,iline))/(line2(1,association(iline)) - line(1,iline))
            point = line(:,iline) + s*(line2(:,association(iline)) - line(:,iline))
            call add_vector(points, point)
          endif
        endif
      enddo
      deallocate(line, line2, breaks, association)
    enddo
    uptonull = uptonull + sum(int(nperring, int64))*32_int64

    print*, 'sorting hcs points'
    print*, size(points, 2)

    call sortpoints(points, ihcs, 5)
    deallocate(points)
  enddo
  write(5) -1
  close(5)
  close(10)
  close(20)

end