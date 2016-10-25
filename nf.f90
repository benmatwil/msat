program nullfinder
  use params 
  use nf_mod

  implicit none

  character(len=30) :: intfmt = "(a,'(',i4,',',i4,',',i4')',a)"
  character(len=40) :: dblfmt= "(a,'(',f9.5,',',f9.5,',',f9.5,')',a)"

  double precision, allocatable, dimension(:) :: xs, ys, zs !x,y and z coordinates of confirmed nulls
  double precision, allocatable, dimension(:) :: xp, yp, zp
  
  double precision, dimension(:,:,:), allocatable :: bx, by, bz

  integer, dimension(:,:,:), allocatable :: candidates !array containing coords of candidate cells (1 if candidate, 0 if not)
  double precision, dimension(2,2,2) :: cube, cbx, cby, cbz !arrays containing vertices surrounding a cell
  double precision :: x, y, z !coordinates of a null
  double precision, allocatable, dimension(:) :: xgrid, ygrid, zgrid

  integer :: nx, ny, nz !number of vertices
  integer :: nx1, ny1, nz1 !number of cells

  integer :: i, j, k
  integer :: ii, jj, kk

  integer :: nnulls, nvert
  integer :: itestx, itesty, itestz, itest

  integer :: ierror
  double precision :: dx
  double precision, allocatable, dimension(:,:) :: distances

  print*,'#######################################################################'
  print*,'#                         Null Finding Code                           #'
  print*,'#                      (Written by G P S Gibb)                        #'
  print*,'#######################################################################'

  call filenames

  open(unit=10, file=filein, access='stream', status='old')

    print*, ''
    print*, "Reading in data from: '", trim(filein), "'..."

    read(10), nx, ny, nz !number of vertices
    print*, ''
    write(*, intfmt), 'Number of points in grid: (nx,ny,nz) = ', nx, ny, nz, ''

    allocate(bx(nx,ny,nz), by(nx,ny,nz), bz(nx,ny,nz))
    allocate(xgrid(nx), ygrid(ny), zgrid(nz))

    read(10), bx
    read(10), by
    read(10), bz
    read(10), xgrid, ygrid, zgrid

  close(10)

  print*, 'Done!'
  print*, ''
  print*, '-----------------------------------------------------------------------'
  print*,  ''

  nx1 = nx-1 !number of cells (not vertices)
  ny1 = ny-1
  nz1 = nz-1

  allocate(candidates(nx1,ny1,nz1)) ! integer array where candidate cells are flagged

  nnulls = 0

  print*, 'Checking for candidate nulls...'

  !loop over each cell checking if bx, by and bz switch sign
  do k = 1, nz1
    do j = 1, ny1
      do i = 1, nx1

        candidates(i,j,k) = 0

        !check if bx changes sign
        cube(:,:,:) = bx(i:i+1, j:j+1, k:k+1)
        itestx = switch(cube)
        if (itestx == 0) cycle

        !check by changes sign
        cube(:,:,:) = by(i:i+1, j:j+1, k:k+1)
        itesty = switch(cube)
        if (itesty == 0) cycle

        !check bz changes sign
        cube(:,:,:) = bz(i:i+1, j:j+1, k:k+1)
        itestz = switch(cube)
        if (itestz == 0) cycle

        cube(:,:,:) = bx(i:i+1, j:j+1, k:k+1)**2 + by(i:i+1, j:j+1, k:k+1)**2 + bz(i:i+1, j:j+1, k:k+1)**2
        if (sum(sqrt(cube)) < 8*zero) cycle

        candidates(i,j,k) = 1
        nnulls = nnulls+1

      enddo
    enddo
  enddo

  print*, ''
  write(*,"(a,i5)"), ' Number of candidate nulls found = ', nnulls

  print*, ''
  print*, '-----------------------------------------------------------------------'
  print*, ''

  print*, 'Now looping over all candidate nulls...'
  !now loop over candidates and test for nulls using bilinear method on faces
  do k = 1, nz1
    do j = 1, ny1
      do i = 1, nx1
        if (candidates(i,j,k) == 1) then

          cbx = bx(i:i+1 ,j:j+1, k:k+1) !bx cell
          cby = by(i:i+1 ,j:j+1, k:k+1) !by cell
          cbz = bz(i:i+1 ,j:j+1, k:k+1) !bz cell
          
          call normalise(cbx, cby, cbz)

          call bilin_test(cbx, cby, cbz, itest) !check for null within (or on face/corner/edge of) cell

          if (itest == 0) then
            candidates(i,j,k) = 0 !if no null found remove this candidate
          endif

        endif
      enddo
    enddo     
  enddo

  print*, ''
  print*, '-----------------------------------------------------------------------'
  print*, ''

  allocate(xs(0), ys(0), zs(0))
  allocate(xp(0), yp(0), zp(0))

  nnulls = sum(candidates)

  !find coordinates of null
  write(*,"(a,i5)"), 'Number of confirmed nulls = ', sum(candidates)
  print*, ''
  print*, "Now determining each null's location to sub-gridcell accuracy:"
  do k = 1, nz1
    do j = 1, ny1
      do i = 1, nx1
        if (candidates(i,j,k) == 1) then

          cbx = bx(i:i+1, j:j+1, k:k+1)
          cby = by(i:i+1, j:j+1, k:k+1)
          cbz = bz(i:i+1, j:j+1, k:k+1)

          call normalise(cbx, cby, cbz)

          x = 0d0
          y = 0d0
          z = 0d0

          dx = 1.
          !Print*, 'Attempting to use trilinear method within cell...'
          do ii = 1, sig_figs
            dx = dx/10.000
            call subgrid(cbx, cby, cbz, dx, x, y, z, ierror) !chop cell into 10x10x10 subcells, and check for nulls in each subcell
            if (ierror == 1) exit
          enddo

          if (ierror == 1) then
            x = x + 5.*dx
            y = y + 5.*dx
            z = z + 5.*dx
            ! print*, 'Attemptig to use Newton Raphson...'
            ! call newton_raphson(cbx, cby, cbz, x, y, z,ierror)
            ! if (ierror .eq. 0) then
            !   Print*, 'Success!'
            !   print*, 'Null at',i+x,j+y,k+z
            ! endif

            print*, i, j, k
            print*, "Don't believe this null. Removing from list"
            nnulls = nnulls-1

            ! print*, 'Attempting to use brute force...'
            ! call brute_force(cbx,cby,cbz,x,y,z)
            ! write(*,dblfmt)'Null at ',i+x,j+y,k+z,' in gridcell coordinates'
          else
            ! write(*,dblfmt) 'Null at ', i+x, j+y, k+z, ' in gridcell coordinates'
            ! write(*,dblfmt) 'Null at ', linear(x, xgrid(i:i+1)), linear(y, ygrid(j:j+1)), linear(z, zgrid(k:k+1)), ' in real units'

            !add this null to the list of nulls (in gridcell coordinates)
            call add_element(xs, x+i)
            call add_element(ys, y+j)
            call add_element(zs, z+k)
            !add this null to the lost of nulls (in 'real' units)
            call add_element(xp, linear(x, xgrid(i:i+1)))
            call add_element(yp, linear(y, ygrid(j:j+1)))
            call add_element(zp, linear(z, zgrid(k:k+1)))

            cbx = bx(i:i+1, j:j+1, k:k+1)
            cby = by(i:i+1, j:j+1, k:k+1)
            cbz = bz(i:i+1, j:j+1, k:k+1)

            !print*, 'Sanity check:'
            ! write(*,"(a,ES10.3)"), '|B| at point is', sqrt(trilinear_cell(x, y, z, cbx)**2 + &
            ! trilinear_cell(x, y, z, cby)**2+trilinear_cell(x, y, z, cbz)**2)

          endif
        endif
      enddo
    enddo
  enddo

  print*, 'Found', nnulls, 'nulls'

  print*, ''
  print*, '-----------------------------------------------------------------------'
  print*, ''

  print*, 'Now checking for nulls of vertices:'
  nvert = 0
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        if (sqrt(bx(i,j,k)**2 + by(i,j,k)**2 + bz(i,j,k)**2) == 0) then ! 0 or zero?
          call add_element(xs, dble(i))
          call add_element(ys, dble(j))
          call add_element(zs, dble(k))
          call add_element(xp, xgrid(i))
          call add_element(yp, ygrid(j))
          call add_element(zp, zgrid(k))
          nvert = nvert + 1
        endif
      enddo
    enddo
  enddo
  print*, 'Found', nvert, 'nulls at vertices'
  nnulls = nnulls + nvert

  print*, ''
  print*, '-----------------------------------------------------------------------'
  print*, ''

  print*, 'Now checking for duplicate nulls:'
  call remove_duplicates(xs, ys, zs, nnulls, xp, yp, zp)

  print*, ''
  print*, 'Final number of nulls=', nnulls
  print*, 'Is this the same?', size(xp,1)

  print*, "Now writing null positions to "//trim(fileout)//"-nullpos.dat"

  open(unit=10, file='output/'//trim(fileout)//'-nullpos.dat', access='stream')
    write(10) nnulls

    do i = 1, nnulls
      write(10) xs(i), ys(i), zs(i)
    enddo
    do i = 1, nnulls
      write(10) xp(i), yp(i), zp(i)
    enddo
  close(10)

  print*, ''
  print*, 'Done!'

end program
