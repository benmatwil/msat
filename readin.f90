!This code reads in the output from ssfind and writes every 'skip' rings to file to be plotted by the idl routines. It also trace the separators found, and writes them to file for readin by the IDL routine
module readin_mod

integer, parameter :: skip = 30 !number of rings to skip

integer :: nrings, npoints, nringsmax
integer, allocatable :: nperring(:)
integer, allocatable :: assoc(:)
integer :: association
integer*8 :: pos

integer :: preamble

double precision :: r(3)
double precision, allocatable :: ring(:,:), brk(:)
double precision, allocatable, dimension(:) :: x, y, z

integer :: opt
integer :: nring, index

contains

  !gets the position numbers in the file for ring number 'nring' and point number 'index' within that ring. The position numbers are:
  !a - the association number (which point in ring (index-1) the point in ring (nring) came from
  !p - the position of the ring vector
  subroutine position(nring,index,a,b,p)
    integer*8 :: a, b, p
    integer :: uptoring

    !1 whole ring contains 3*np vector points and np association points
    uptoring = preamble + sum(nperring(0:nring-1))*8
    a = uptoring + (index-1)*1
    b = uptoring + nperring(nring)*1 + (index-1)*1
    p = uptoring + nperring(nring)*2 + (index-1)*6

    a = a*4 + 1
    b = b*4 + 1
    p = p*4 + 1

  end subroutine

end module

!********************************************************************************

program readin
  use readin_mod
  use params
  implicit none

  integer*8 :: a, b, p
  integer :: iring, inull, nrcount

  integer :: nnulls, nullnum_f, nullnum_t

  double precision :: sep, sep1, sep2, sep3
  double precision, allocatable, dimension(:,:) :: rnulls

  character (len=8), parameter :: fmt='(I4.4)'
  character (len=5) :: fname

  print*,'#######################################################################'
  print*,'#                      Data Read-in and Process                       #'
  print*,'#######################################################################'

  open(unit=10,file='output/null.dat',form='unformatted')
    read(10) nnulls
    print*,'num nulls',nnulls

    allocate(rnulls(3,nnulls))

    do inull=1,nnulls
      read (10) rnulls(:,inull)
      print*, 'null position', inull, rnulls(:,inull)
    enddo
  close(10)

  do inull = 1, nnulls
    write(fname,fmt) inull

    open(unit=12,file='output/separator'//trim(fname)//'.dat',access='stream')
    open(unit=10,file='output/everything'//trim(fname)//'.dat',access='stream')

    read(10) nrings, npoints, nringsmax
    allocate(nperring(0:nringsmax-1))
    read(10) nperring

    open(unit=11,file='output/sep'//trim(fname)//'.dat',access='stream')

    print*, 'nrings, npoints, nringsmax'
    print*, nrings, npoints, nringsmax-1

    preamble = 3 + nringsmax !takes you up to end of 'header'

    do
      read(12) opt
      if (opt < 0) exit !if no separator in file, or reached the end of the file
      read(12) nullnum_f, nullnum_t
      print*, 'Separator starts from null', nullnum_f
      print*, 'Separator ends at null', nullnum_t

      read(12) nring, index !ring number and index within the ring the start point of the separator belongs to
      print*, 'ring number and index within ring where separator starts'
      print*, nring, index

      allocate(x(0:nring+1), y(0:nring+1), z(0:nring+1))

      print*, 'Tracing Separator...'
      print*, 'Start ring is', nring

      do iring = nring, 0, -1 !trace back each ring to the start

        call position(iring,index,a,b,p) !get positions in file belonging to the current ring

        read(10, pos=a) association !read association for point number index in ring 1
        read(10, pos=p) r !read position of point index on ring iring
        if (iring == nring) then 
          x(nring+1) = rnulls(1,nullnum_t)
          y(nring+1) = rnulls(2,nullnum_t)
          z(nring+1) = rnulls(3,nullnum_t)
          print*, 'End null, last point and start of separator'
          print*, x(nring+1), y(nring+1), z(nring+1)
          print*, r
        endif
        !write ring position vector into x, y and z
        x(iring) = r(1)
        y(iring) = r(2)
        z(iring) = r(3)
        index = association

        if (iring == 1) print*, r

      enddo

      write(11) nring+2
      write(11) nringsmax+1
      write(11) x, y, z

      deallocate(x, y, z)
      print*, ''
    enddo

    write(11) -1
    close(11)
    close(12)

    ! Just writes a selection of rings to file i.e. every "skip"ped number of rings
    print*, 'writing rings to file'

    open(unit=11, file='output/ringidl'//trim(fname)//'.dat', form='unformatted')
    print*, 'nrings=', nrings
    nrcount = 0
    do iring = 1, nrings-2, skip
      nrcount = nrcount+1
    enddo
    write(11), nrcount

    do iring = 1, nrings-2, skip
      call position(iring,1,a,b,p)
      allocate(x(nperring(iring)), y(nperring(iring)), z(nperring(iring)), ring(3,nperring(iring)))
      allocate(brk(nperring(iring)))
      read(10, pos=b) brk
      read(10, pos=p) ring
      
      x = ring(1,:)
      y = ring(2,:)
      z = ring(3,:)

      write(11) nperring(iring)
      write(11) x, y, z
      write(11) brk
      deallocate(x, y, z, ring, brk)
    enddo
    write(11) -1
    close(11)
    close(10)

    deallocate(nperring)
  enddo

  print*, "Completed readin successfully"

end
