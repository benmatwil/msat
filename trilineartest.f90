program trilineartester
use sfmod

  double precision, allocatable :: field(:,:,:,:)
  integer :: nr, nt, np, check
  double precision :: r(3), tri(3), tri1(3), tri2(3), checknum(3)

  open(unit=10,file='field.dat',access='stream')
  read(10) nr, nt, np
  allocate(field(nr,nt,np,3))
  read(10) field
  close(10)
  
  r = [100d0,200d0,100d0]
  tri1 = trilinear(r,field)
  tri2 = trilinear([r(1)+1,r(2),r(3)],field)
  
  do i = 1, 3
    tri = trilinear([r(1)+0.25*i,r(2),r(3)],field)
    checknum = (tri-tri1)/(tri2-tri1)
    print*,checknum
  enddo
  
  print*,'--------------------------------------------------------------'
  
  tri1 = trilinear(r,field)
  tri2 = trilinear([r(1),r(2),r(3)+1],field)
  
  do i = 1, 3
    tri = trilinear([r(1),r(2),r(3)+0.25*i],field)
    checknum = (tri-tri1)/(tri2-tri1)
    print*,checknum
  enddo
  
  print*,trilinear([100d0,1d0,1d0],field)
  
end program