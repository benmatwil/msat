!code to find the separatrix surface and separators of null points
program ssfinder

#if _OPENMP
  use omp_lib
#endif

  use params
  use common
  use trace
  use ring

  implicit none

  integer(int64) :: tstart, tstop, count_rate !to time program

  real(np) :: r(3)
  real(np) :: h, h0
  real(np) :: slowdown, slowdist
  integer(int32) :: nearflag

  integer(int32) :: iring, iline, inull, jnull

  ! null parameters
  real(np) :: theta, phi, factor
  real(np), allocatable, dimension(:) :: xs, ys, zs

  ! spine tracer
  real(np), dimension(3) :: rspine
  real(np), dimension(:,:), allocatable :: spine
  real(np) :: hspine
  integer(int32) :: dir

  integer(int32) :: nlines
  integer(int32) :: nrings
  integer(int32) :: nperring(0:ringsmax)

  logical :: exitcondition, out

  ! file writing
  integer(int64) :: ia, ib, ip, uptonullring, uptonullconn
  integer(int32) :: isep, nullnum1, nullnum2, linenum, ringnum
  real(np), dimension(:,:), allocatable :: rsep, rring, write_ring
  integer(int32), dimension(:), allocatable :: brk
  character(:), allocatable :: tempfile

  ! for restarting
  integer(int64) :: filesize, totalpts, filepos
  integer(int32) :: istart, temp, ringsmax1

  print*,'#######################################################################'
  print*,'#                      Separatrix Surface Finder                      #'
  print*,'#######################################################################'

#if _OPENMP
  if (nproc > 0) call omp_set_num_threads(nproc)
  print*, 'Using', nproc, 'processors'
#endif

  call filenames

  ! read in data
  open(unit=10, file=filein, access='stream', status='old')
    read(10) nx, ny, nz ! number of vertices
    allocate(bgrid(nx, ny, nz, 3))
    allocate(x(nx), y(ny), z(nz))
    read(10) bgrid(:,:,:,1)
    read(10) bgrid(:,:,:,2)
    read(10) bgrid(:,:,:,3)
    read(10) x, y, z
  close(10)

  xmin = 1
  xmax = nx
  ymin = 1
  ymax = ny
  zmin = 1
  zmax = nz
  print*, "Data read in, grid dimensions are", nx, ny, nz

  call system_clock(tstart,count_rate) ! to time how long it takes

  ! read in null data

  open(unit=10, file=trim(fileout)//'-nullpos.dat', access='stream', status='old')
    read(10) nnulls
    allocate(rnulls(3,nnulls), rnullsreal(3,nnulls))
    read(10) rnulls
    read(10) rnullsreal
  close(10)
  open(unit=10, file=trim(fileout)//'-nulldata.dat', access='stream', status='old')
    read(10) nnulls
    allocate(signs(nnulls),spines(3,nnulls),fans(3,nnulls))
    read(10) signs, spines, fans
  close(10)

  rnullsalt = rnulls

#if spherical
  ! for sphericals
  do inull = 1, nnulls
    ! check whether null is at the lower phi boundary
    if (rnullsalt(3,inull) < zmin + 1) then
      rnullsalt(3,inull) = rnullsalt(3,inull) - 1
      call edgecheck(rnullsalt(:,inull))
      rnullsalt(3,inull) = rnullsalt(3,inull) + 1
    endif

    ! check whether null is at the upper phi boundary
    if (rnullsalt(3,inull) > zmax - 1) then
      rnullsalt(3,inull) = rnullsalt(3,inull) + 1
      call edgecheck(rnullsalt(:,inull))
      rnullsalt(3,inull) = rnullsalt(3,inull) - 1
    endif

    ! check whether null is at the lower theta boundary
    if (rnullsalt(2,inull) < ymin + 1) then
      rnullsalt(2,inull) = rnullsalt(2,inull) - 1
      call edgecheck(rnullsalt(:,inull))
      rnullsalt(2,inull) = rnullsalt(2,inull) - 1
    endif

    ! check whether null is at the upper theta boundary
    if (rnullsalt(2,inull) > ymax - 1) then
      rnullsalt(2,inull) = rnullsalt(2,inull) + 1
      call edgecheck(rnullsalt(:,inull))
      rnullsalt(2,inull) = rnullsalt(2,inull) + 1
    endif
  enddo
#elif cylindrical
  ! for cylindricals
  do inull = 1, nnulls
    ! check whether null is at lower phi boundary
    if (rnullsalt(2,inull) < ymin + 1) then
      rnullsalt(2,inull) = rnullsalt(2,inull) - 1
      call edgecheck(rnullsalt(:,inull))
      rnullsalt(2,inull) = rnullsalt(2,inull) + 1
    endif

    ! check whether null is at upper phi boundary
    if (rnullsalt(2,inull) > ymax - 1) then
      rnullsalt(2,inull) = rnullsalt(2,inull) + 1
      call edgecheck(rnullsalt(:,inull))
      rnullsalt(2,inull) = rnullsalt(2,inull) - 1
    endif
  enddo
#endif

  uptonullring = 0
  uptonullconn = 1

  if (restart) then
    ! stores all info regarding size of rings
    open(unit=10,file=trim(fileout)//'-ringinfo.dat',status='old',access='stream')
    ! stores all the coordinates of rings in original coordinate system
    open(unit=20,file=trim(fileout)//'-rings.dat',status='old',access='stream')
    ! stores all the connection info about each separator
    open(unit=40,file=trim(fileout)//'-connectivity.dat',status='old',access='stream')
    ! stores all the coordinates of the separator lines in original coordinate system
    open(unit=50,file=trim(fileout)//'-separators.dat',status='old',access='stream')
    ! stores all the coordinates of the spine lines in original coordinate system
    open(unit=60,file=trim(fileout)//'-spines.dat',status='old',access='stream')

    ! find number of nulls done
    read(10) ringsmax1
    if (ringsmax1 /= ringsmax) then
      print*, "Please set ringsmax to", ringsmax1, "in params.f90"
      stop
    endif
    inquire(unit=10, size=filesize)
    if (nullrestart == 0) then
      istart = (filesize-4_int64)/4_int64/ringsmax1 - 10
      if (istart < 1) istart = 1
    else
      istart = nullrestart
    endif
    print*, 'Restarting from null', istart

    ! sort out ringinfo, rings and spine files
    filepos = 1
    totalpts = 0
    do inull = 1, istart-1
      read(10) nperring
      totalpts = totalpts + sum(nperring)
      do dir = 1, 2
        ! print*, filepos
        read(60, pos=filepos) ringnum
        ! print*, inull, ringnum, filepos
        filepos = filepos + ringnum*24_int64 + 4_int64
      enddo
    enddo
    write(60, pos=filepos)
    uptonullring = totalpts*32_int64
    
    ! sort out connectivity and separator files
    filepos = 1
    inquire(unit=40, size=filesize)
    nseps = filesize/16_int64
    do isep = 1, nseps
      read(40, pos=uptonullconn) nullnum1
      ! print*, nullnum1
      if (nullnum1 > istart-1) then
        write(40, pos=uptonullconn)
        write(50, pos=filepos)
        exit
      endif
      read(50, pos=filepos) ringnum
      uptonullconn = uptonullconn + 16_int64
      filepos = filepos + 24_int64*ringnum + 4_int64
    enddo
    
  else
    ! stores all info regarding size of rings
    open(unit=10,file=trim(fileout)//'-ringinfo.dat',status='replace',access='stream')
    ! stores all the coordinates of rings in original coordinate system
    open(unit=20,file=trim(fileout)//'-rings.dat',status='replace',access='stream')
    ! stores all the connection info about each separator
    open(unit=40,file=trim(fileout)//'-connectivity.dat',status='replace',access='stream')
    ! stores all the coordinates of the separator lines in original coordinate system
    open(unit=50,file=trim(fileout)//'-separators.dat',status='replace',access='stream')
    ! stores all the coordinates of the spine lines in original coordinate system
    open(unit=60,file=trim(fileout)//'-spines.dat',status='replace',access='stream')

    write(10) ringsmax+1, nskip, bytesize, stepsize
    print*, ringsmax+1, nskip, bytesize, stepsize

    istart = 1
  endif

  tempfile = trim(fileout)//'-rings.temp'

  !$OMP PARALLEL private(iring, iline, inull)
  do inull = istart, nnulls ! loop over all nulls

    !$OMP SINGLE

    open(unit=90, file=tempfile, status='replace', access='stream')
    write(20, pos=uptonullring+1)
    write(40, pos=uptonullconn)
    print*, ''
    print*, 'Null number', inull, 'of', nnulls

    ! get number of start points
    nlines = nstart

    theta = acos(fans(3,inull))
    phi = atan(fans(2,inull), fans(1,inull))
    allocate(xs(nlines), ys(nlines), zs(nlines))
    call get_startpoints(theta,phi, xs,ys,zs)

    allocate(line1(3,nlines), line2(3,nlines))
    allocate(break(nlines))

    out = .false.

    ! add start points to first ring relative to null
    do iline = 1, nlines !Go through each start point
      line1(1,iline) = rnulls(1,inull) + xs(iline)
      line1(2,iline) = rnulls(2,inull) + ys(iline)
      line1(3,iline) = rnulls(3,inull) + zs(iline)
      call edgecheck(line1(:,iline))
      if (outedge(line1(:,iline))) then
        factor = 0.9_np
        do while (outedge(line1(:,iline)))
          line1(1,iline) = rnulls(1,inull) + xs(iline)*factor
          line1(2,iline) = rnulls(2,inull) + ys(iline)*factor
          line1(3,iline) = rnulls(3,inull) + zs(iline)*factor
          call edgecheck(line1(:,iline))
          factor = factor*0.9_np
        enddo
      endif
    enddo
    line2 = line1

    break = 0
    nseps = 0
    nperring = 0
    nperring(0) = nlines
    slowdown = 1.0_np
    terror = 0

    write_ring = gtr(line1, x, y, z)
    write(20) [(iline,iline=1,nlines)], break, real(write_ring, bytesize)
    write(90) [(iline,iline=1,nlines)], break, write_ring
    deallocate(write_ring)

    exitcondition = .false.
    !$OMP END SINGLE

    do iring = 1, ringsmax ! loop over number of rings we want

      if (signs(inull) .eq. 0) then ! skip null which is uncharacterised
        print*,'Null has zero sign'
        exit
      endif
      
      !$OMP SINGLE
      nrings = iring
#if debug
      print*, iring, nlines
#endif

      allocate(endpoints(nlines), association(nlines))

      if (iring < 50) then
        h0 = stepsize/slowdown/5
      else
        h0 = stepsize/slowdown
      endif
      !$OMP END SINGLE

      !$OMP WORKSHARE
      endpoints = 0
      !$OMP END WORKSHARE

      !$OMP DO private(r, h, out)
      do iline = 1, nlines ! loop over all points in ring (in parallel do)

        r = line1(:,iline)
        h = h0

        call trace_line(r, signs(inull), h) ! trace line by a distance of h

        line2(:,iline) = line2(:,iline) + (r - line1(:,iline))
        line1(:,iline) = r

        ! counter to see how many points on ring have reached outer boundary
        if (outedge(r)) then
          endpoints(iline) = 1
        else
          endpoints(iline) = 0
        endif

        association(iline) = iline

      enddo
      !$OMP END DO

      !$OMP SINGLE
      if (iring < 50) then
        maxdist = 0.1_np*h0*slowdown
      elseif (iring < 100) then
        maxdist = 0.08_np*h0*slowdown
      else
        maxdist = 0.5_np*h0*slowdown
      endif
      nulldist = 2.0_np*h0*slowdown
      mindist = maxdist/3
      slowdist = 2.0_np*nulldist
      !$OMP END SINGLE

      ! remove points from ring if necessary
      call remove_points(nlines)

      !$OMP SINGLE
      allocate(nearnull(nlines))
      slowdown = 1.0_np
      !$OMP END SINGLE

      !$OMP WORKSHARE
      nearnull = 0
      !$OMP END WORKSHARE

      !$OMP DO private(jnull)
      do iline = 1, nlines
        do jnull = 1, nnulls
          if (signs(jnull) == signs(inull)) cycle
          if (dist(rnulls(:,jnull), line1(:, iline)) < slowdist .or. &
            dist(rnullsalt(:,jnull), line1(:, iline)) < slowdist) then
            nearnull(iline) = 1
            exit
          endif
        enddo
      enddo
      !$OMP END DO

      !$OMP SINGLE
      nearflag = sum(nearnull)
      if (nearflag > 0) then
        slowdown = 2.0_np
      endif
      !$OMP END SINGLE

      ! add points to ring if necessary
      call add_points(nlines)

      !$OMP SINGLE
      do iline = 1, nlines
        call edgecheck(line1(:,iline))
      enddo
      !$OMP END SINGLE

      !$OMP SINGLE
      deallocate(nearnull)
      !$OMP END SINGLE

      ! determine if any lines are separators
      if (nearflag > 0) call sep_detect(nlines,inull,iring,signs(inull))

      !$OMP SINGLE

      if (nlines > pointsmax) then
      ! exit if too many points on ring
        print*, 'Reached maximum number of points'
        exitcondition = .true.
      endif

      if (nlines == 0) then
        ! exit if all points have reached outer boundary (left box)
        print*, 'All points have left the box'
        exitcondition = .true.
      endif

      if (terror == 1) then
        print*, 'Tracing has failed', iring
        exitcondition = .true.
      endif

      if (exitcondition) deallocate(endpoints, association)
      !$OMP END SINGLE

      if (exitcondition) then
        nrings = iring
        exit
      endif

      !$OMP SINGLE
      nperring(iring) = nlines
      write_ring = gtr(line1, x, y, z)
      if (modulo(iring, nskip) == 0) write(20) association, break, real(write_ring, bytesize)
      write(90) association, break, write_ring
      ! if (inull == 2) print*, association

      deallocate(endpoints, association, write_ring)
      !$OMP END SINGLE

    enddo

    !$OMP SINGLE
    write(10) nperring(::nskip)
    print*, "Tracing spines and any separators"

    ! trace separators...
    write(40, pos=uptonullconn)
    if (nseps > 0) then
      do isep = 1, nseps
        read(40) nullnum1, nullnum2, ringnum, linenum
        allocate(rsep(3, ringnum+3))
        do iring = ringnum, 0, -1
          call file_position(iring, nperring, linenum, ia, ib, ip)
          read(90, pos=ia) linenum
          read(90, pos=ip) rsep(:,iring+2)
        enddo
        rsep(:,1) = rnullsreal(:,nullnum1)
        rsep(:,ringnum+3) = rnullsreal(:,nullnum2)
        write(50) ringnum+3
        write(50) rsep
        deallocate(rsep)
      enddo
    endif

    ! Trace spines...
    do dir = -1, 1, 2
      out = .false.
      allocate(spine(3,1))
      spine(:,1) = rnulls(:, inull)
      rspine = rnulls(:, inull) + dir*spines(:, inull)*1e-3_np ! pick a good factor
      do while (.not. outedge(rspine))
        if (signs(inull) == 0) exit
        call add_vector(spine, rspine)
        hspine = 5e-2_np
        call trace_line(rspine, -signs(inull), hspine)
        call edgecheck(rspine)
      enddo
      write(60) size(spine,2), gtr(spine, x, y, z)
      deallocate(spine)
    enddo

    if (nrings == ringsmax) print*, "Reached maximum number of rings"

    print*, 'Number of separators:', nseps, 'Number of rings:', nrings

    deallocate(xs, ys, zs)
    deallocate(line1, line2, break)
    uptonullring = uptonullring + sum(int(nperring(::nskip), int64))*32_int64
    uptonullconn = uptonullconn + nseps*16_int64 ! 4*4

    close(90, status='delete')
  
    !$OMP END SINGLE

  enddo

  !$OMP END PARALLEL

  close(10)
  close(20)
  write(40) -1
  close(40)
  close(50)
  close(60)

  call system_clock(tstop, count_rate)
  print*, 'Time taken:', dble(tstop - tstart)/dble(count_rate), "(", dble(tstop - tstart)/dble(count_rate)/60, "minutes)"

end program