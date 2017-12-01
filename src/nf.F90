program nullfinder

  use params
  use nf_mod
  use common

  implicit none

  real(np), dimension(:, :, :), allocatable :: magb
  logical, dimension(:, :, :), allocatable :: candidates ! flags for candidate cells

  real(np), dimension(2, 2, 2, 3) :: cube
  real(np) :: ds
  real(np), dimension(0:10, 0:10, 0:10, 3) :: subcube
  integer(int32) :: mincube(3)
  
  real(np) :: rtri(3), rnull(3)

  integer(int32) :: ix, iy, iz, isig, ixsub, iysub, izsub

  integer(int32) :: itest

  integer(int32) :: ierror
  real(np) :: bound_dist

  integer(int32) :: i, j, counter
  real(np), allocatable :: distances(:,:)
  logical :: exitcondition

  print*, '#######################################################################'
  print*, '#                         Null Finding Code                           #'
  print*, '#######################################################################'

  call filenames

  open(unit=10, file=filein, access='stream', status='old')

    print*, ''
    print*, 'Reading in data from: '//trim(filein)

    read(10) nx, ny, nz ! number of vertices
    print*, 'Number of points in grid: ', nx, ny, nz

    allocate(bgrid(nx, ny, nz, 3), magb(nx, ny, nz))
    allocate(x(nx), y(ny), z(nz))

    read(10) bgrid(:, :, :, 1)
    read(10) bgrid(:, :, :, 2)
    read(10) bgrid(:, :, :, 3)
    read(10) x, y, z

  close(10)

  xmax = nx
  ymax = ny
  zmax = nz

  magb = sqrt(sum(bgrid(:, :, :, :)**2, 4))

  allocate(candidates(nx-1, ny-1, nz-1)) ! array where candidate cells are flagged
  candidates = .false.

  ! loop over each cell checking if bx, by and bz switch sign
  do iz = 1, nz-1
    do iy = 1, ny-1
      do ix = 1, nx-1
        
        cube = bgrid(ix:ix+1, iy:iy+1, iz:iz+1, :)
        ! check if bx changes sign
        if (allsame(cube(:, :, :, 1))) cycle
        ! check by changes sign
        if (allsame(cube(:, :, :, 2))) cycle
        ! check bz changes sign
        if (allsame(cube(:, :, :, 3))) cycle

        ! check cube isn't basically zero everywhere
        if (sum(magb(ix:ix+1, iy:iy+1, iz:iz+1)) < 8*zero) cycle

        candidates(ix, iy, iz) = .true.

      enddo
    enddo
  enddo

  nnulls = count(candidates)

  print*, ''
  print*, 'Number of candidate cells = ', nnulls

  print*, ''
  print*, '-----------------------------------------------------------------------'
  print*, ''

  allocate(rnulls(3, 0))

  do iz = 1, nz-1
    do iy = 1, ny-1
      do ix = 1, nx-1

        if (candidates(ix, iy, iz)) then
          
          rnull = [ix, iy, iz]

          ierror = 0
          sigloop: do isig = 1, sig_figs
            
            ds = 10.0_np**(-isig)
            do izsub = 0, 10
              do iysub = 0, 10
                do ixsub = 0, 10
                  ! find all field values of subgrid
                  rtri = rnull + [ixsub, iysub, izsub]*ds
                  subcube(ixsub, iysub, izsub, :) = trilinear(rtri, bgrid)
                enddo
              enddo
            enddo

            outer: do izsub = 0, 9
              do iysub = 0, 9
                do ixsub = 0, 9
                  itest = 0
                  cube = subcube(ixsub:ixsub+1, iysub:iysub+1, izsub:izsub+1, :)
                  ! check if bx changes sign
                  if (allsame(cube(:, :, :, 1))) cycle
                  ! check by changes sign
                  if (allsame(cube(:, :, :, 2))) cycle
                  ! check bz changes sign
                  if (allsame(cube(:, :, :, 3))) cycle

                  ! run bilinear test
                  call bilin_test(cube(:, :, :, 1), cube(:, :, :, 2), cube(:, :, :, 3), itest)

                  if (itest == 1) then
                    rnull = rnull + [ixsub, iysub, izsub]*ds
                    exit outer
                  endif
                enddo
              enddo
            enddo outer
            
            if (itest == 0) then
#if debug
              if (isig > 1) then
                print*, 'Candidate null in cell', ix, iy, iz
                print*, "Don't believe this null. ERROR :("
                print*, 'Point to', isig, 'sigfigs:', real(rnull, real32)
                print*, '|B| at point:', sqrt(sum(trilinear(rnull, bgrid)**2))
                print*, ''
              endif
#endif
              exit sigloop
            endif
            
          enddo sigloop

          if (itest == 1) then
            ds = 10.0_np**(-sig_figs)
            do izsub = 1, 2
              do iysub = 1, 2
                do ixsub = 1, 2
                  ! find all field values of subgrid
                  rtri = rnull + [ixsub-1, iysub-1, izsub-1]*ds
                  cube(ixsub, iysub, izsub, :) = trilinear(rtri, bgrid)
                enddo
              enddo
            enddo

            mincube = minloc(sum(cube**2, 4))

            rnull = rnull + (mincube - 1)*ds

            bound_dist = rspherefact*10**(-sig_figs)
            
            ! if (rnull(1) > 1 + bound_dist .and. rnull(1) < nx - bound_dist .and. &
            !   rnull(2) > 1 + bound_dist .and. rnull(2) < ny - bound_dist .and. &
            !   rnull(3) > 1 + bound_dist .and. rnull(3) < nz - bound_dist) then
            !   call add_vector(rnulls, rnull)
            ! endif
            if (rnull(1) > 1 + bound_dist .and. rnull(1) < nx - bound_dist) then
              call add_vector(rnulls, rnull)
#if debug
            else
              print*, 'Null on the boundary', int([ix, iy, iz], int16), 'Removing...'
#endif
            endif
          endif
        endif
      enddo
    enddo
  enddo

  nnulls = size(rnulls, 2)

  print*, 'Found', nnulls, 'nulls'

  print*, ''
  print*, '-----------------------------------------------------------------------'
  print*, ''

  do iz = 1, nz
    do iy = 1, ny
      do ix = 1, nx
        if (magb(ix, iy, iz) < zero) then ! 0 or zero?
          call add_vector(rnulls, real([ix, iy, iz], np))
        endif
      enddo
    enddo
  enddo

  nnulls = size(rnulls, 2)

  print*, 'Found', size(rnulls, 2) - nnulls, 'nulls at vertices'

  print*, ''
  print*, '-----------------------------------------------------------------------'
  print*, ''

  print*, 'Now checking for duplicate nulls:'
    
  if (nnulls > 1) then
    exitcondition = .false.

    main_while: do while (.not. exitcondition)
      !allocate(distances(n,n))
      do j = 1, nnulls
        do i = j+1, nnulls
          if (dist(rnulls(:, i), rnulls(:, j)) < 1e-4_np) then
            print*, 'removing duplcate at index', j, dist(rnulls(:, i), rnulls(:, j))
            call remove_vector(rnulls, j)
            nnulls = nnulls - 1
            cycle main_while
          endif
        enddo
      enddo
      exitcondition = .true.
    enddo main_while

    allocate(distances(nnulls, nnulls))

    distances = 1e6_np
    do j = 1, nnulls
      do i = j+1, nnulls
        distances(i, j) = dist(rnulls(:, i), rnulls(:, j))
      enddo
    enddo

    print*,''
    print*,'All null pairs with less than 1 gridcell spacing'
    counter = 0
    do j = 1, nnulls
      do i = j+1, nnulls
        if (distances(i, j) < 1) then
          print *, i, j, distances(i, j)
          counter = counter + 1
        endif
      enddo
    enddo

    deallocate(distances)
    print*, 'Number of close-by nulls =', counter

  endif

  print*, ''
  print*, '-----------------------------------------------------------------------'
  print*, ''

  print*, ''
  print*, 'Final number of nulls=', nnulls
  print*, 'Is this the same?', size(rnulls, 2)

  print*, "Writing null positions to "//trim(fileout)//"-nullpos.dat"

  open(unit=10, file=trim(fileout)//'-nullpos.dat', access='stream', status='replace')
    write(10) nnulls
    write(10) rnulls
    write(10) gtr(rnulls, x, y, z)
  close(10)

end program
