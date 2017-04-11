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

    print*, 'entering'
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
      if (dists(imindist(1), imindist(2)) > disttol*10) exit
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

      write(unit) inull, size(pordered, 2), pordered
      deallocate(pordered)
    enddo
print*, 'exit'
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
  integer(int32), dimension(:), allocatable :: association
  integer(int32), dimension(0:ringsmax) :: nperring
  real(np), dimension(:,:), allocatable :: line, line2, points
  real(np) :: s, point(3)

  integer(int32) :: nspine, dir

  integer(int32) :: ringlun, ringinfolun, ringcutlun, spinelun, seplun, conlun, sepcutlun, spinecutlun
  integer(int32) :: nulllun, fieldlun

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

  open(newunit=nulllun, file=trim(fileout)//'-nullpos.dat', access='stream', status='old')
    read(nulllun) nnulls
  close(nulllun)

  open(newunit=fieldlun,file=filein,access='stream',status='old')
    read(fieldlun), nx, ny, nz !number of vertices
    allocate(xg(nx), yg(ny), zg(nz))
    ig = 3_int64*4_int64 + 3_int64*int(nx, int64)*int(ny, int64)*int(nz, int64)*8_int64 + 1_int64
    read(fieldlun, pos=ig) xg, yg, zg
  close(fieldlun)

  !********************************************************************************

  ! for each separator, check whether we cross r=r0 and add point to file
  print*, 'Separators'
  open(newunit=sepcutlun, file=trim(fileout)//'-separators-cut_'//rstr//'.dat', access='stream', status='replace')
  open(newunit=conlun, file=trim(fileout)//'-connectivity.dat', access='stream', status='old')
  open(newunit=seplun, file=trim(fileout)//'-separators.dat', access='stream', status='old')
  
  read(conlun) flag
  do while (flag > 0)
    read(seplun) npoints
    allocate(line(3,npoints))
    read(seplun) line
    do iline = 1, npoints-1
      nextline = iline+1
      if ((line(1,iline)-r0)*(line(1,nextline)-r0) < 0) then
        s = (r0 - line(1,iline))/(line(1,nextline) - line(1,iline))
        point = line(:,iline) + s*(line(:,nextline) - line(:,iline))
        write(sepcutlun) flag, point
      endif
    enddo
    deallocate(line)
    read(conlun) flag, flag, flag, flag
  enddo
  write(sepcutlun) -1
  close(sepcutlun)
  close(seplun)

  !********************************************************************************

  ! for each separator, check whether we cross r=r0 and add point to file
  open(newunit=sepcutlun, file=trim(fileout)//'-hcs-separators-cut_'//rstr//'.dat', access='stream', status='replace')
  open(newunit=conlun, file=trim(fileout)//'-hcs-connectivity.dat', access='stream', status='old')
  open(newunit=seplun, file=trim(fileout)//'-hcs-separators.dat', access='stream', status='old')
  
  read(conlun) flag
  do while (flag > 0)
    read(seplun) npoints
    allocate(line(3,npoints))
    read(seplun) line
    do iline = 1, npoints-1
      nextline = iline+1
      if ((line(1,iline)-r0)*(line(1,nextline)-r0) < 0) then
        s = (r0 - line(1,iline))/(line(1,nextline) - line(1,iline))
        point = line(:,iline) + s*(line(:,nextline) - line(:,iline))
        write(sepcutlun) flag, point
      endif
    enddo
    deallocate(line)
    read(conlun) flag, flag, flag, flag
  enddo
  write(sepcutlun) -1
  close(sepcutlun)
  close(seplun)

  !********************************************************************************
  
  ! for each spine, check whether we cross r=r0 and add point to file
  print*, 'Spines'
  open(newunit=spinecutlun, file=trim(fileout)//'-spines-cut_'//rstr//'.dat', access='stream', status='replace')
  open(newunit=spinelun, file=trim(fileout)//'-spines.dat', access='stream', status='old')

  allocate(spines(3, 0))
  
  do inull = 1, nnulls
    ! read in spine
    do dir = 1, 2
      read(spinelun) nspine
      allocate(line(3,nspine))
      read(spinelun) line
      do iline = 1, nspine-1
        nextline = iline+1
        if ((line(1,iline)-r0)*(line(1,nextline)-r0) < 0) then
          s = (r0 - line(1,iline))/(line(1,nextline) - line(1,iline))
          point = line(:,iline) + s*(line(:,nextline) - line(:,iline))
          write(spinecutlun) inull, point
          call add_vector(spines, point)
        endif
      enddo
      deallocate(line)
    enddo
  enddo
  write(spinecutlun) -1
  close(spinecutlun)
  close(spinelun)

  !********************************************************************************

  ds = maxval([maxval(yg(2:ny) - yg(1:ny-1)), maxval(zg(2:nz) - zg(1:nz-1))])

  print*, 'Rings'
  open(newunit=ringcutlun, file=trim(fileout)//'-rings-cut_'//rstr//'.dat', access='stream', status='replace')
  open(newunit=ringinfolun, file=trim(fileout)//'-ringinfo.dat', access='stream', status='old')
  read(ringinfolun) nrings
  open(newunit=ringlun, file=trim(fileout)//'-rings.dat', access='stream', status='old')

  uptonull = 0
  do inull = 1, nnulls
    if (mod(inull, 20) == 0) print*, inull
    allocate(points(3,0))
    read(ringinfolun) nperring
    nrings = count(nperring > 0) - 1
    do iring = nrings, 1, -1
    
      call file_position(iring, nperring, 1, ia, ib, ip)
      allocate(association(nperring(iring)), line(3,nperring(iring)))
      read(ringlun, pos=uptonull+ia) association, line
      read(ringlun, pos=uptonull+ip) line

      call file_position(iring-1, nperring, 1, ia, ib, ip)
      allocate(line2(3,nperring(iring-1)))
      read(ringlun, pos=uptonull+ip) line2

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
      deallocate(line, line2, association)
    enddo
    uptonull = uptonull + sum(int(nperring, int64))*32_int64

    ! print*, size(points, 2)
    call sortpoints(points, inull, ringcutlun)
    
    deallocate(points)
    print*, 'deallocated points', inull
  enddo
  write(ringcutlun) -1
  close(ringcutlun)
  close(ringinfolun)
  close(ringlun)

  print*, 'finsihed'

  !********************************************************************************

  ! for hcs
  open(newunit=ringinfolun, file=trim(fileout)//'-hcs-ringinfo.dat', access='stream', status='old')
  open(newunit=ringlun, file=trim(fileout)//'-hcs-rings.dat', access='stream', status='old')
  open(newunit=ringcutlun, file=trim(fileout)//'-hcs-cut_'//rstr//'.dat', access='stream', status='replace')
  read(ringinfolun) nrings
  
  uptonull = 0
  ihcs = 1
  do dir = 1, 2
  print*, dir
    allocate(points(3,0))
    read(ringinfolun) nperring
    nrings = count(nperring > 0) - 1
    do iring = nrings, 1, -1
    print*, iring
    
      call file_position(iring, nperring, 1, ia, ib, ip)
      allocate(association(nperring(iring)), line(3,nperring(iring)))
      read(ringlun, pos=uptonull+ia) association, line

      call file_position(iring-1, nperring, 1, ia, ib, ip)
      allocate(line2(3,nperring(iring-1)))
      read(ringlun, pos=uptonull+ip) line2

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
      deallocate(line, line2, association)
    enddo
    uptonull = uptonull + sum(int(nperring, int64))*32_int64

    print*, 'sorting hcs points'
    print*, size(points, 2)

    call sortpoints(points, ihcs, ringcutlun)
    deallocate(points)
  enddo
  write(ringcutlun) -1
  close(ringcutlun)
  close(ringinfolun)
  close(ringlun)

end