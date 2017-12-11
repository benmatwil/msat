! spine finder with convergence method
program signfinder

#if _OPENMP
  use omp_lib
#endif

  use params
  use sf_mod
  use common

  implicit none

  integer(int32), dimension(:), allocatable :: warnings

  integer(int32) :: nullstart, nullend

  integer(int32) :: pcount, ncount, ucount

  integer(int32) :: sign, warning
  real(np), dimension(3) :: spine, fan

  integer(int32) :: inull

  character(20) :: arg
  integer(int32) :: iarg

  print*,'#######################################################################'
  print*,'#                             Spinefinder                             #'
  print*,'#######################################################################'

#if _OPENMP
  if (nproc > 0) call omp_set_num_threads(nproc)
  print*, 'Using', nproc, 'processors'
#endif

  call filenames

  !Read in 'null.dat'
  open(unit=10, file=trim(fileout)//'-nullpos.dat', access='stream', status='old')
    read(10) nnulls
    allocate(rnulls(3,nnulls),spines(3,nnulls),fans(3,nnulls),signs(nnulls),warnings(nnulls))
    read(10) rnulls!(:,inull)
  close(10)

  nullstart = 1
  nullend = nnulls
  if (command_argument_count() > 0) then
    do iarg = 1, command_argument_count()
      call get_command_argument(iarg,arg)
      if (trim(arg) == '-n') then
        call get_command_argument(iarg+1,arg)
        read(arg,*) nullstart
        nullend = nullstart
        savedata = .true.
      endif
    enddo
  endif

  !read in bgrid
  open(unit=10, file=filein, access='stream', status='old')
    read(10) nx, ny, nz
    allocate(bgrid(nx, ny, nz, 3))
    read(10) bgrid
    allocate(x(nx), y(ny), z(nz))
    read(10) x, y, z
  close(10)

  xmax = nx
  ymax = ny
  zmax = nz

  print*, 'There are ', nnulls,' nulls to analyse'
  print*, '-----------------------------------------------------------------------------'

  !now loop over each null and characterise
  !$omp parallel do private(inull, sign, spine, fan, warning)
  do inull = nullstart, nullend!1, nnulls

    call get_properties(inull, sign, spine, fan, warning)

    signs(inull) = sign
    spines(:, inull) = spine
    fans(:, inull) = fan
    warnings(inull) = warning

    print*, ''
    print*, '-----------------------------------------------------------------------------'
    print*, ''
  enddo
  if (nullend - nullstart /= nnulls - 1) stop

  !now write data to nulls.dat
  open(unit=10, file=trim(fileout)//'-nulldata.dat', access='stream')
    write(10) nnulls
    write(10) signs, spines, fans, warnings
  close(10)

  print*, 'Summary:'
  print*, 'number, sign, warning'

  pcount = 0
  ucount = 0
  ncount = 0
  do inull = 1, nnulls
    print*, inull, signs(inull), warnings(inull)
    if (signs(inull) .eq. 1) pcount = pcount+1
    if (signs(inull) .eq. -1) ncount = ncount+1
    if (signs(inull) .eq. 0) ucount = ucount+1
  enddo

  print*, 'Total number of nulls:', nnulls
  print*, 'Positive', pcount
  print*, 'Negative', ncount
  print*, 'Unknown', ucount
  print*, 'Warning', count(warnings > 0)
  print*, 'All nulls have data, should be equal', pcount+ncount+ucount, nnulls
  print*, 'Done!'

end program
