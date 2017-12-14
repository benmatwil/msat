module make_cut_mod

  use params
  use common

  implicit none

  real(np) :: ds

  contains

  subroutine sortpoints(points, pordered, disttol, exit_status)

    real(np), dimension(:, :), allocatable :: dists, spinedists, dots, spinedots
    real(np), dimension(:, :), allocatable :: points, pordered
    real(np) :: disttol, diff(3), mindist, mindistspine

    integer(int32), dimension(2) :: imindist, iminspinedist, imindot, iminspinedot
    integer(int32) :: ip, jp, dir, npoints

    logical :: exit_status

    exit_status = .false.

    if (size(points, 2) > 5000) then
      npoints = 5000
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
    if (dists(imindist(1), imindist(2)) > disttol*10) then
      exit_status = .true.
    else
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

      do dir = 1, 2
        ! add points in the each direction
        do while (size(points, 2) > 0)
          allocate(dists(size(points, 2), 1))
          allocate(spinedists(size(spines, 2), 1))
          diff = normalise(pordered(:, 1) - pordered(:, 2))
          do ip = 1, size(points, 2)
            dists(ip, 1) = dist(points(:, ip), pordered(:, 1))
          enddo
          do ip = 1, size(spines, 2)
            spinedists(ip, 1) = dist(spines(:, ip), pordered(:, 1))
          enddo
          imindist = minloc(dists)
          iminspinedist = minloc(spinedists)
          mindist = dists(imindist(1), 1)
          mindistspine = spinedists(iminspinedist(1), 1)
          deallocate(dists, spinedists)
          if (mindist < mindistspine) then
            if (dot(diff, normalise(points(:, imindist(1)) - pordered(:, 1))) < -0.98) then
              exit
            elseif (mindist < disttol) then
              call add_vector(pordered, points(:, imindist(1)), 1)
              call remove_vector(points, imindist(1))
            elseif (dot(diff, normalise(points(:, imindist(1)) - pordered(:, 1))) > 0.95 .and. &
              mindist < 5*disttol) then
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
    endif

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

  integer(int64) :: ig
  real(np), dimension(:), allocatable :: xg, yg, zg

  integer(int32) :: iarg
  character(10) :: arg

  integer(int64) :: ia, ip, uptorings, uptoassoc, uptoring1, uptoring2
  integer(int32) :: inull, jnull, iring, iline, nextline, flag, ihcs, ipt, isep
  integer(int32) :: nrings, npoints, nlines, nseps, ncomp
  integer(int32), dimension(:), allocatable :: association
  integer(int32), dimension(:), allocatable :: nperring
  real(np), dimension(:, :), allocatable :: line, line2, points, pordered
  real(np) :: s, point(3)
  logical :: iexist, exit_status

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

  open(unit=10, file=filein, access='stream', status='old')
    read(10), nx, ny, nz !number of vertices
    allocate(xg(nx), yg(ny), zg(nz))
    ig = 3_int64*4_int64 + 3_int64*int(nx, int64)*int(ny, int64)*int(nz, int64)*8_int64 + 1_int64
    read(10, pos=ig) xg, yg, zg
  close(10)

  !********************************************************************************

  ! for each separator, check whether we cross r=r0 and add point to file
  print*, 'Separators'
  open(unit=80, file=trim(fileout)//'-separators-cut_'//rstr//'.dat', access='stream', status='replace')
  open(unit=20, file=trim(fileout)//'-connectivity.dat', access='stream', status='old')
  open(unit=30, file=trim(fileout)//'-separators.dat', access='stream', status='old')

  ipt = 0
  write(80) ipt
  do inull = 1, nnulls
    read(20) nseps
    do isep = 1, nseps
      read(20) flag, flag
      read(30) npoints
      allocate(line(3, npoints))
      read(30) line
      do iline = 1, npoints-1
        nextline = iline+1
        if ((line(1, iline) - r0)*(line(1, nextline) - r0) < 0) then
          s = (r0 - line(1,iline))/(line(1,nextline) - line(1,iline))
          point = line(:,iline) + s*(line(:,nextline) - line(:,iline))
          write(80) inull, flag, point
          ipt = ipt + 1
        endif
      enddo
      read(20) flag, flag
      deallocate(line)
    enddo
  enddo
  
  write(80, pos=1) ipt
  close(80)
  close(30)
  close(20)

  !********************************************************************************

  ! for each separator, check whether we cross r=r0 and add point to file
  open(unit=80, file=trim(fileout)//'-hcs-separators-cut_'//rstr//'.dat', access='stream', status='replace')
  open(unit=20, file=trim(fileout)//'-hcs-connectivity.dat', access='stream', status='old')
  open(unit=30, file=trim(fileout)//'-hcs-separators.dat', access='stream', status='old')
  
  ipt = 0
  write(80) ipt
  read(20) nseps
  do isep = 1, nseps
    read(20) inull, flag
    read(30) npoints
    allocate(line(3,npoints))
    read(30) line
    do iline = 1, npoints-1
      nextline = iline+1
      if ((line(1, iline) - r0)*(line(1, nextline) - r0) < 0) then
        s = (r0 - line(1, iline))/(line(1 ,nextline) - line(1, iline))
        point = line(:, iline) + s*(line(:, nextline) - line(:, iline))
        write(80) flag, point
        ipt = ipt + 1
      endif
    enddo
    read(20) flag, flag
    deallocate(line)
  enddo
  
  write(80, pos=1) ipt
  close(80)
  close(30)
  close(20)

  !********************************************************************************
  
  ! for each spine, check whether we cross r=r0 and add point to file
  print*, 'Spines'
  open(unit=80, file=trim(fileout)//'-spines-cut_'//rstr//'.dat', access='stream', status='replace')
  open(unit=10, file=trim(fileout)//'-spines.dat', access='stream', status='old')

  allocate(spines(3, 0))

  ipt = 0
  write(80) ipt
  
  do inull = 1, nnulls
    ! read in spine
    do dir = 1, 2
      read(10) nspine
      allocate(line(3,nspine))
      read(10) line
      do iline = 1, nspine-1
        nextline = iline + 1
        if ((line(1, iline) - r0)*(line(1, nextline) - r0) < 0) then
          s = (r0 - line(1, iline))/(line(1, nextline) - line(1, iline))
          point = line(:, iline) + s*(line(:, nextline) - line(:, iline))
          write(80) inull, point
          ipt = ipt + 1
          call add_vector(spines, point)
        endif
      enddo
      deallocate(line)
    enddo
  enddo
  write(80, pos=1) ipt
  close(80)
  close(10)

  !********************************************************************************

  ds = maxval([maxval(yg(2:ny) - yg(1:ny-1)), maxval(zg(2:nz) - zg(1:nz-1))])

  print*, 'Rings'
  open(unit=80, file=trim(fileout)//'-rings-cut_'//rstr//'.dat', access='stream', status='replace')
  open(unit=40, file=trim(fileout)//'-ringinfo.dat', access='stream', status='old')
  open(unit=50, file=trim(fileout)//'-rings.dat', access='stream', status='old')
  inquire(file=trim(fileout)//'-assocs.dat', exist=iexist)
  if (iexist) then
    open(unit=55, file=trim(fileout)//'-assocs.dat', access='stream', status='old')
  else
    stop "Associations file does not exist."
  endif

  nlines = 0

  write(80) nlines

  read(40) nrings
  allocate(nperring(0:nrings-1))
  read(40), nrings, nrings, nrings, nrings

  uptoassoc = 0
  uptorings = 0
  do inull = 1, nnulls
    if (mod(inull, 20) == 0) print*, inull
    allocate(points(3,0))
    read(40) nperring
    nrings = count(nperring > 0) - 1 ! index of nperring starts at 0
    do iring = nrings, 1, -1

      uptoring1 = sum(int(nperring(0:iring-2), int64))
      uptoring2 = uptoring1 + nperring(iring-1)
      
      allocate(association(nperring(iring)))
      ia = uptoassoc + uptoring2*4_int64 + 1
      read(55, pos=ia) association
      
      if (iring == nrings) then
        allocate(line(3, nperring(iring)))
        ip = uptorings + uptoring2*24_int64 + 1
        read(50, pos=ip) line
      endif

      allocate(line2(3, nperring(iring-1)))
      ip = uptorings + uptoring1*24_int64 + 1
      read(50, pos=ip) line2

      do iline = 1, nperring(iring)
        if ((line(1, iline) - r0)*(line2(1, association(iline)) - r0) < 0) then
          if (abs(line(3, iline) - line2(3, association(iline))) < 1) then
            ! find point in between at r0 and add to line
            s = (r0 - line(1, iline))/(line2(1, association(iline)) - line(1, iline))
            point = line(:, iline) + s*(line2(:, association(iline)) - line(:, iline))
            call add_vector(points, point)
          endif
        endif
      enddo
      call move_alloc(line2, line)
      deallocate(association)
    enddo
    deallocate(line)
    uptoring1 = sum(int(nperring, int64))
    uptoassoc = uptoassoc + uptoring1*4_int64
    uptorings = uptorings + uptoring1*24_int64

    do while (size(points, 2) > 0)
      
      call sortpoints(points, pordered, ds/3, exit_status)
      if (exit_status) then
        exit
      else
        write(80) inull, size(pordered, 2), pordered
        nlines = nlines + 1
        deallocate(pordered)
      endif

    enddo
    
    deallocate(points)
  enddo
  write(80, pos=1) nlines
  close(80)
  close(40)
  close(50)
  close(55)
  deallocate(nperring)

  print*, 'Finished null rings'

  !********************************************************************************

  ! for hcs
  open(unit=80, file=trim(fileout)//'-hcs-cut_'//rstr//'.dat', access='stream', status='replace')
  open(unit=40, file=trim(fileout)//'-hcs-ringinfo.dat', access='stream', status='old')
  open(unit=50, file=trim(fileout)//'-hcs-rings.dat', access='stream', status='old')
  inquire(file=trim(fileout)//'-hcs-assocs.dat', exist=iexist)
  if (iexist) then
    open(unit=55, file=trim(fileout)//'-hcs-assocs.dat', access='stream', status='old')
  else
    stop "HCS Associations file does not exist."
  endif

  nlines = 0
  write(80) nlines
  
  read(40) nrings
  allocate(nperring(0:nrings-1))

  read(40), nrings, nrings, nrings, nrings, ncomp
  
  uptoassoc = 0
  uptorings = 0

  ! uptonull = 0
  do ihcs = 1, ncomp
    allocate(points(3, 0))
    read(40) nperring

    nrings = count(nperring > 0) - 1
    do iring = nrings, 1, -1

      uptoring1 = sum(int(nperring(0:iring-2), int64))
      uptoring2 = uptoring1 + nperring(iring-1)

      ia = uptoassoc + uptoring2*4_int64 + 1
      ip = uptorings + uptoring2*24_int64 + 1
    
      allocate(association(nperring(iring)), line(3, nperring(iring)))
      read(55, pos=ia) association
      read(50, pos=ip) line

      ip = uptorings + uptoring1*24_int64 + 1

      allocate(line2(3,nperring(iring-1)))
      read(50, pos=ip) line2

      do iline = 1, nperring(iring)
        if ((line(1, iline) - r0)*(line2(1, association(iline)) - r0) < 0) then
          if (abs(line(3, iline) - line2(3, association(iline))) < 1) then
            ! find point in between at r0 and add to line
            s = (r0 - line(1, iline))/(line2(1, association(iline)) - line(1, iline))
            point = line(:, iline) + s*(line2(:, association(iline)) - line(:, iline))
            call add_vector(points, point)
          endif
        endif
      enddo
      deallocate(line, line2, association)
    enddo
    uptoring1 = sum(int(nperring, int64))
    uptoassoc = uptoassoc + uptoring1*4_int64
    uptorings = uptorings + uptoring1*24_int64

    print*, 'sorting hcs points'
    print*, size(points, 2)

    do while (size(points, 2) > 0)
      
      call sortpoints(points, pordered, ds/3, exit_status)
      if (exit_status) then
        exit
      else
        write(80) int(0, int32), size(pordered, 2), pordered
        nlines = nlines + 1
        deallocate(pordered)
      endif

    enddo

    deallocate(points)
  enddo
  write(80, pos=1) nlines
  close(80)
  close(40)
  close(50)
  close(55)
  deallocate(nperring)

end program