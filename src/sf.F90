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
  print*, ''

  call filenames

  open(unit=10, file=trim(fileout)//'-nullpos.dat', access='stream', status='old')
    read(10) nnulls
    allocate(rnulls(3, nnulls))
    read(10) rnulls
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
  open(unit=20, file=filein, access='stream', status='old')
    read(20) nx, ny, nz
    allocate(bgrid(nx, ny, nz, 3))
    read(20) bgrid
    allocate(x(nx), y(ny), z(nz))
    read(20) x, y, z
  close(20)

  xmax = nx
  ymax = ny
  zmax = nz

  allocate(spines(3, nnulls), fans(3, nnulls), signs(nnulls), warnings(nnulls))
  spines = 0
  fans = 0
  signs = 0
  warnings = 0

  ! now loop over each null and characterise
  !$omp parallel do private(inull, sign, spine, fan, warning)
  do inull = nullstart, nullend!1, nnulls

    call get_properties(inull, sign, spine, fan, warning)

    signs(inull) = sign
    spines(:, inull) = spine
    fans(:, inull) = fan
    warnings(inull) = warning

  enddo
  if (nullend - nullstart /= nnulls - 1) stop

  ! now write data
  open(unit=10, file=trim(fileout)//'-nulldata.dat', access='stream')
    write(10) nnulls
    write(10) signs, spines, fans, warnings
  close(10)

  do inull = 1, nnulls
    if (warnings(inull) /= 0) print*, 'Warning on null', inull, ':', warnings(inull)
  enddo

  print*, 'Total number of nulls:', nnulls
  print*, 'Positive', count(signs == 1)
  print*, 'Negative', count(signs == -1)
  print*, 'Unknown', count(abs(signs) /= 1)
  print*, 'Warning 1:', count(warnings == 1)
  print*, 'Warning 2:', count(warnings == 2)
  print*, 'Warning 3:', count(warnings == 3)
  print*, 'Warning 4:', count(warnings == 4)

end program
