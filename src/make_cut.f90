module make_cut_mod

  use params
  use common
  use list_mod

  implicit none

  real(np) :: ds, disttol

  real(np) :: n_pln(3), d_pln

  contains

  subroutine sortpoints(points, pordered, exit_status)

    real(np), dimension(:, :), allocatable :: dists
    real(np), dimension(:), allocatable :: ptdists, spinedists
    real(np), dimension(:, :), allocatable :: points, pordered
    real(np), dimension(3) :: diff, pt1, pt2
    real(np) :: mindist, mindistspine

    integer(int32), dimension(2) :: imindists
    integer(int32) :: imindist, iminspinedist
    integer(int32) :: ip, jp, dir, npoints

    logical :: exit_status

    type(list) :: pt_list

    exit_status = .false.

    if (size(points, 2) > 5000) then
      npoints = 5000
    else
      npoints = size(points, 2)
    endif
    allocate(dists(npoints, npoints))
    dists = 1000 ! maybe needs improving
    do ip = 1, npoints
      dists(ip, ip+1:npoints) = sqrt((points(1, ip) - points(1, ip+1:npoints))**2 &
        + (points(2, ip) - points(2, ip+1:npoints))**2 + (points(3, ip) - points(3, ip+1:npoints))**2)
    enddo
    imindists = minloc(dists)
    
    if (dists(imindists(1), imindists(2)) > disttol*10) then
      exit_status = .true.
    else
      deallocate(dists)
      call pt_list%create()
      call pt_list%append(points(:, imindists(1)))
      call pt_list%append(points(:, imindists(2)))
      call remove_vector(points, imindists(1))
      if (imindists(1) < imindists(2)) then
        call remove_vector(points, imindists(2)-1)
      else
        call remove_vector(points, imindists(2))
      endif

      do dir = 1, 2
        ! add points in the each direction
        do while (size(points, 2) > 0)
          allocate(ptdists(size(points, 2)), spinedists(size(spines, 2)))

          pt1 = pt_list%first%r
          pt2 = pt_list%first%next%r

          diff = normalise(pt1 - pt2)
          ptdists = (points(1, :) - pt1(1))**2 + (points(2, :) - pt1(2))**2 + (points(3, :) - pt1(3))**2
          spinedists = (spines(1, :) - pt1(1))**2 + (spines(2, :) - pt1(2))**2 + (spines(3, :) - pt1(3))**2

          imindist = minloc(ptdists, 1)
          iminspinedist = minloc(spinedists, 1)
          mindist = ptdists(imindist)
          mindistspine = spinedists(iminspinedist)

          deallocate(ptdists, spinedists)

          if (mindist < mindistspine) then
            if (dot(diff, normalise(points(:, imindist) - pt1)) < -0.98) then
              exit
            elseif (mindist < disttol) then
              call pt_list%prepend(points(:, imindist))
              call remove_vector(points, imindist)
            elseif (dot(diff, normalise(points(:, imindist) - pt1)) > 0.95 .and. &
              mindist < 5*disttol) then
              call pt_list%prepend(points(:, imindist))
              call remove_vector(points, imindist)
            else
              exit
            endif
          else
            exit
          endif
        enddo

        call pt_list%reverse()

      enddo

      pordered = pt_list%to_array()
      call pt_list%destroy()
    endif

  end subroutine

  function plane(r)

    real(np) :: plane, r(3)

    plane = dot(n_pln, r) - d_pln

  end function

  function check_crossing(r1, r2)

    logical :: check_crossing
    real(np), dimension(3) :: r1, r2

    check_crossing = plane(r1)*plane(r2) <= 0

  end function

  function find_crossing(r1, r2)

    real(np), dimension(3) :: find_crossing, r1, r2
    real(np) :: s

    s = (d_pln - dot(r1, n_pln))/(dot(r2 - r1, n_pln))
    find_crossing = r1 + s*(r2 - r1)

  end function

end module

program make_cut

  use params
  use common
  use make_cut_mod

  implicit none

  real(np) :: r0
  integer(int32) :: pc = 0
  character(:), allocatable :: astr, bstr, cstr, dstr, sign
  character(*), parameter :: fmt = '(A,i3.3,f5.4)', neg = '         ', pos = '        '
  character(:), allocatable :: file_pln

  integer(int64) :: ig
  real(np), dimension(:), allocatable :: xg, yg, zg

  integer(int32) :: iarg, jarg
  character(10) :: arg

  type(list) :: spine_list, ring_list
  integer(int64) :: ia, ip, uptorings, uptoassoc, uptoring1, uptoring2
  integer(int32) :: inull, jnull, iring, iline, nextline, flag, ihcs, ipt, isep
  integer(int32) :: nrings, npoints, nlines, nseps, ncomp, nskip_file
  ! integer(int32), dimension(:), allocatable :: association
  integer(int32), dimension(:), allocatable :: association, break
  integer(int32), dimension(:), allocatable :: nperring
  real(np), dimension(:, :), allocatable :: line, line2, points, pordered
  real(np) :: s, point(3), hstep
  logical :: iexist, exit_status, do_hcs = .false., do_all_ssf = .true., too_much_memory

  integer(int32) :: nspine, dir

  integer(int32) :: iprint, iorder

  integer(int64) :: num_nums, leftover_num, read_index, read_num
  integer(int64), parameter :: large_read = 500000000_int64
  real(np), dimension(:, :), allocatable :: lines
  ! integer(int32), dimension(:), allocatable :: associations
  integer(int32), dimension(:), allocatable :: associations, breaks

  call filenames

  if (command_argument_count() > 0) then
    do iarg = 1, command_argument_count()
      call get_command_argument(iarg,arg)
      if (trim(arg) == '--hcs') then
        do_hcs = .true.
      endif
      if (trim(arg) == '--only-hcs') then
        do_all_ssf = .false.
      endif
      if (trim(arg) == '-n') then
        do jarg = 1, 3
          call get_command_argument(iarg+jarg,arg)
          read(arg,*) n_pln(jarg)
        enddo
        pc = pc + 1
      endif
      if (trim(arg) == '-d') then
        call get_command_argument(iarg+1,arg)
        read(arg,*) d_pln
        pc = pc + 1
      endif
    enddo
  endif
  if (pc < 2) then
    stop "Need to provide normal to the plane -n and planar constant -d"
  endif
  
  if (n_pln(1) < 0) then
    sign = '-'
    astr = neg
  else
    sign = ''
    astr = pos
  endif
  write(astr, fmt) sign, int(abs(n_pln(1))), abs(n_pln(1))-int(abs(n_pln(1)))
  if (n_pln(2) < 0) then
    sign = '-'
    bstr = neg
  else
    sign = ''
    bstr = pos
  endif
  write(bstr, fmt) sign, int(abs(n_pln(2))), abs(n_pln(2))-int(abs(n_pln(2)))
  if (n_pln(3) < 0) then
    sign = '-'
    cstr = neg
  else
    sign = ''
    cstr = pos
  endif
  write(cstr, fmt) sign, int(abs(n_pln(3))), abs(n_pln(3))-int(abs(n_pln(3)))
  if (d_pln < 0) then
    sign = '-'
    dstr = neg
  else
    sign = ''
    dstr = pos
  endif
  write(dstr, fmt) sign, int(abs(d_pln)), abs(d_pln)-int(abs(d_pln))

  file_pln = astr // '_' // bstr // '_' // cstr // '_' // dstr

  print*, 'Making cut at in file '// filein // ' using plane:'
  print '(5X, 7A)', &
    astr, ' * x1 + ', bstr, ' * x2 + ', cstr, ' * x3 = ', dstr

  open(unit=10, file=trim(fileout)//'-nullpos.dat', access='stream', status='old')
    read(10) nnulls
  close(10)

  open(unit=10, file=filein, access='stream', status='old')
    read(10) nx, ny, nz !number of vertices
    allocate(xg(nx), yg(ny), zg(nz))
    ig = 3_int64*4_int64 + 3_int64*int(nx, int64)*int(ny, int64)*int(nz, int64)*8_int64 + 1_int64
    read(10, pos=ig) xg, yg, zg
  close(10)

  if (nnulls < 10) then
    iprint = 1
  else
    iorder = floor(log10(1.0_np*nnulls))
    iprint = (nnulls/(10**iorder) + 1)*10**(iorder-1)
  endif

  call system_clock(tstart, count_rate) ! to time how long it takes

  open(unit=40, file=trim(fileout)//'-ringinfo.dat', access='stream', status='old')
    read(40) nrings, nrings, nrings, hstep
  close(40)

  ds = maxval([maxval(xg(2:nx) - xg(1:nx-1)), &
               maxval(yg(2:ny) - yg(1:ny-1)), &
               maxval(zg(2:nz) - zg(1:nz-1))])*2*hstep
  disttol = ds/10
  print*, ds, hstep
  !********************************************************************************
    
  ! for each spine, check whether we cross r=r0 and add point to file
  print*, 'Spines'
  open(unit=80, file=trim(fileout)//'-spines-cut_'//file_pln//'.dat', access='stream', status='replace')
  open(unit=10, file=trim(fileout)//'-spines.dat', access='stream', status='old')

  call spine_list%create()

  ipt = 0
  write(80) ipt
  
  do inull = 1, nnulls
    ! read in spine
    do dir = 1, 2
      read(10) nspine
      allocate(line(3, nspine))
      read(10) line
      do iline = 1, nspine-1
        nextline = iline + 1
        if (check_crossing(line(:, iline), line(:, nextline))) then
          if (dist(line(:, iline), line(:, nextline)) < 1) then
            point = find_crossing(line(:, iline), line(:, nextline))
            write(80) inull, point
            ipt = ipt + 1
            call spine_list%append(point)
          endif
        endif
      enddo
      deallocate(line)
    enddo
  enddo
  write(80, pos=1) ipt
  close(80)
  close(10)

  spines = spine_list%to_array()

  !********************************************************************************

  if (do_all_ssf) then
    ! for each separator, check whether we cross r=r0 and add point to file
    print*, 'Separators'
    open(unit=80, file=trim(fileout)//'-separators-cut_'//file_pln//'.dat', access='stream', status='replace')
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
          if (check_crossing(line(:, iline), line(:, nextline))) then
            if (dist(line(:, iline), line(:, nextline)) < 1) then
              point = find_crossing(line(:, iline), line(:, nextline))
              write(80) inull, flag, point
              ipt = ipt + 1
            endif
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

    print*, 'Rings'
    open(unit=80, file=trim(fileout)//'-rings-cut_'//file_pln//'.dat', access='stream', status='replace')
    open(unit=40, file=trim(fileout)//'-ringinfo.dat', access='stream', status='old')
    open(unit=50, file=trim(fileout)//'-rings.dat', access='stream', status='old')
    open(unit=60, file=trim(fileout)//'-breaks.dat', access='stream', status='old')
    inquire(file=trim(fileout)//'-assocs.dat', exist=iexist)
    if (iexist) then
      open(unit=55, file=trim(fileout)//'-assocs.dat', access='stream', status='old')
    else
      stop "Associations file does not exist."
    endif

    nlines = 0

    write(80) nlines

    read(40) nrings
    read(40) nskip_file
    allocate(nperring(0:ceiling(real(nrings)/nskip_file-1)))
    read(40) nrings, nrings, nrings ! just read into unneeded variable to reach correct position in file

    uptoassoc = 0
    uptorings = 0
    do inull = 1, nnulls
      call ring_list%create()
      read(40) nperring
      nrings = count(nperring > 0) - 1 ! index of nperring starts at 0

      num_nums = sum(int(nperring, int64))
      if (num_nums > 2500000000_int64) then
        too_much_memory = .true.
      else
        too_much_memory = .false.
      endif
      ! print*, num_nums

      if (.not. too_much_memory) then
        ! 500000000 seems to be maximum number of int32 integers fortran can read at a time
        ! allocate(lines(3, num_nums), associations(num_nums))
        allocate(lines(3, num_nums), associations(num_nums), breaks(num_nums))

        ia = uptoassoc + 1
        ip = uptorings + 1
        
        if (num_nums > large_read) then
          leftover_num = num_nums
          read_index = 1
          inquire(50, pos=ip)
          inquire(55, pos=ia)
          inquire(60, pos=ia)
          do while (leftover_num > 0)
            if (leftover_num > large_read) then
              read_num = large_read
            else
              read_num = leftover_num
            endif
            read(50) lines(:, read_index:read_index+read_num-1)
            read(55) associations(read_index:read_index+read_num-1)
            read(60) breaks(read_index:read_index+read_num-1)
            read_index = read_index + read_num
            leftover_num = leftover_num - read_num
          enddo
        else
          read(50, pos=ip) lines
          read(55, pos=ia) associations
          read(60, pos=ia) breaks
        endif
      endif

      do iring = nrings, 1, -1

        uptoring1 = sum(int(nperring(0:iring-2), int64))
        uptoring2 = uptoring1 + nperring(iring-1)
        
        if (too_much_memory) then
          ! allocate(association(nperring(iring)))
          allocate(association(nperring(iring)), break(nperring(iring)))
          ia = uptoassoc + uptoring2*4_int64 + 1
          read(55, pos=ia) association
          read(60, pos=ia) break
          
          if (iring == nrings) then
            allocate(line(3, nperring(iring)))
            ip = uptorings + uptoring2*24_int64 + 1
            read(50, pos=ip) line
          endif

          allocate(line2(3, nperring(iring-1)))
          ip = uptorings + uptoring1*24_int64 + 1
          read(50, pos=ip) line2
        else
          association = associations(uptoring2+1:uptoring2+nperring(iring))
          break = breaks(uptoring2+1:uptoring2+nperring(iring))
          if (iring == nrings) line = lines(:, uptoring2+1:uptoring2+nperring(iring))
          line2 = lines(:, uptoring1+1:uptoring2)
        endif

        do iline = 1, nperring(iring)
          if (check_crossing(line(:, iline), line2(:, association(iline)))) then
            if (dist(line(:, iline), line2(:, association(iline))) < 1) then
              ! find point in between at r0 and add to line
              point = find_crossing(line(:, iline), line2(:, association(iline)))
              call ring_list%append(point)
            endif
          endif
        enddo
        do iline = 1, nperring(iring)-1
          nextline = iline + 1
          if (check_crossing(line(:, iline), line(:, nextline))) then
            if (break(iline) == 0 .and. dist(line(:, iline), line(:, nextline)) < 1) then
              ! find point in between at r0 and add to line
              point = find_crossing(line(:, iline), line(:, nextline))
              call ring_list%append(point)
            endif
          endif
        enddo
        call move_alloc(line2, line)
        ! deallocate(association)
        deallocate(association, break)
      enddo
      ! if (.not. too_much_memory) deallocate(lines, associations)
      if (.not. too_much_memory) deallocate(lines, associations, breaks)
      if (allocated(line)) deallocate(line)
      uptoring1 = sum(int(nperring, int64))
      uptoassoc = uptoassoc + uptoring1*4_int64
      uptorings = uptorings + uptoring1*24_int64

      points = ring_list%to_array()
      call ring_list%destroy()

      do while (size(points, 2) > 0)

        call sortpoints(points, pordered, exit_status)
        if (exit_status) then
          exit
        else
          write(80) inull, size(pordered, 2), pordered
          nlines = nlines + 1
          deallocate(pordered)
        endif

      enddo
      
      deallocate(points)
      if (mod(inull, iprint) == 0) print*, 'Done', inull, 'of', nnulls
    enddo
    write(80, pos=1) nlines
    close(80)
    close(40)
    close(50)
    close(55)
    deallocate(nperring)

    print*, 'Finished null rings'
  
  endif

  !********************************************************************************

  if (do_hcs) then
    ! for hcs
    open(unit=80, file=trim(fileout)//'-hcs-cut_'//file_pln//'.dat', access='stream', status='replace')
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
    read(40) nskip_file
    allocate(nperring(0:ceiling(real(nrings)/nskip_file-1)))
    read(40) nrings, nrings, nrings, ncomp
    
    uptoassoc = 0
    uptorings = 0

    ! uptonull = 0
    do ihcs = 1, ncomp
      allocate(points(3, 0))
      call ring_list%create()
      read(40) nperring

      nrings = count(nperring > 0) - 1
      num_nums = sum(int(nperring, int64))
      allocate(lines(3, num_nums), associations(num_nums))
      ia = uptoassoc + 1
      ip = uptorings + 1

      if (num_nums > large_read) then
        leftover_num = num_nums
        read_index = 1
        inquire(50, pos=ip)
        inquire(55, pos=ia)
        do while (leftover_num > 0)
          if (leftover_num > large_read) then
            read_num = large_read
          else
            read_num = leftover_num
          endif
          read(50) lines(:, read_index:read_index+read_num)
          read(55) associations(read_index:read_index+read_num)
          read_index = read_index + read_num
          leftover_num = leftover_num - read_num
        enddo
      else
        read(50, pos=ip) lines
        read(55, pos=ia) associations
      endif

      do iring = nrings, 1, -1

        uptoring1 = sum(int(nperring(0:iring-2), int64))
        uptoring2 = uptoring1 + nperring(iring-1)

        ! ia = uptoassoc + uptoring2*4_int64 + 1
        ! ip = uptorings + uptoring2*24_int64 + 1
      
        ! allocate(association(nperring(iring)), line(3, nperring(iring)))
        ! read(55, pos=ia) association
        ! read(50, pos=ip) line

        ! ip = uptorings + uptoring1*24_int64 + 1

        ! allocate(line2(3,nperring(iring-1)))
        ! read(50, pos=ip) line2

        association = associations(uptoring2+1:uptoring2+nperring(iring))
        if (iring == nrings) line = lines(:, uptoring2+1:uptoring2+nperring(iring))
        line2 = lines(:, uptoring1+1:uptoring2)

        do iline = 1, nperring(iring)
          if (check_crossing(line(:, iline), line2(:, association(iline)))) then
            if (dist(line(:, iline), line2(:, association(iline))) < 1) then
              ! find point in between at r0 and add to line
              point = find_crossing(line(:, iline), line2(:, association(iline)))
              call ring_list%append(point)
            endif
          endif
        enddo
        call move_alloc(line2, line)
        deallocate(association)
      enddo
      deallocate(lines, associations)
      if (allocated(line)) deallocate(line)
      uptoring1 = sum(int(nperring, int64))
      uptoassoc = uptoassoc + uptoring1*4_int64
      uptorings = uptorings + uptoring1*24_int64

      points = ring_list%to_array()
      call ring_list%destroy()

      print*, 'sorting hcs points'
      print*, size(points, 2)

      do while (size(points, 2) > 0)

        call sortpoints(points, pordered, exit_status)
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

  !********************************************************************************

    ! for each hcs separator, check whether we cross r=r0 and add point to file
    open(unit=80, file=trim(fileout)//'-hcs-separators-cut_'//file_pln//'.dat', access='stream', status='replace')
    open(unit=20, file=trim(fileout)//'-hcs-connectivity.dat', access='stream', status='old')
    open(unit=30, file=trim(fileout)//'-hcs-separators.dat', access='stream', status='old')
    
    ipt = 0
    write(80) ipt
    read(20) nseps
    do isep = 1, nseps
      read(20) inull, flag
      read(30) npoints
      allocate(line(3, npoints))
      read(30) line
      do iline = 1, npoints-1
        nextline = iline+1
        if (check_crossing(line(:, iline), line(:, nextline))) then
          if (dist(line(:, iline), line(:, nextline)) < 1) then
            point = find_crossing(line(:, iline), line(:, nextline))
            write(80) flag, point
            ipt = ipt + 1
          endif
        endif
      enddo
      read(20) flag, flag
      deallocate(line)
    enddo
    
    write(80, pos=1) ipt
    close(80)
    close(30)
    close(20)
  endif

  call system_clock(tstop, count_rate)
  print*, 'Time taken:', dble(tstop - tstart)/dble(count_rate), "(", dble(tstop - tstart)/dble(count_rate)/60, "minutes)"

  call spine_list%destroy()

end program