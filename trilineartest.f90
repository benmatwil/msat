program trilineartester
  use sfmod

  double precision, allocatable :: field(:,:,:,:), tripoints(:,:), tripoints1(:,:), tripoints2(:)
  integer :: nr, nt, np, check
  double precision :: r(3), tri(3), tri1(3), tri2(3), checknum(3)

  open(unit=10,file='data/testfield.dat',access='stream')
  read(10) nr, nt, np
  allocate(field(nr,nt,np,3))
  read(10) field
  close(10)
  
  print*,nr,nt,np, shape(field)
  
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
  
  print*,'---------------------------------------------------------------'
  
  allocate(tripoints(3,4))
  
  r = [100,165,236]
  do i = 0, 1
    !do j = 0, 1
      do k = 0, 1
        tripoints(:,k+1+i*2) = trilinear([r(1)+i,r(2)+j,r(3)+k],field)
        print*, i+1, k+1, tripoints(:,k+1+i*2)
      enddo
    !enddo
  enddo
  
  allocate(tripoints1(3,2))
  
  tripoints1(:,1) = tripoints(:,1)+0.2*(tripoints(:,3)-tripoints(:,1))
  tripoints1(:,2) = tripoints(:,2)+0.2*(tripoints(:,4)-tripoints(:,2))
  
  allocate(tripoints2(3))
  
  tripoints2(:) = tripoints1(:,1)+0.3*(tripoints1(:,2)-tripoints1(:,1))
  
  print*, trilinear([r(1)+0.2,r(2),r(3)+0.3],field)
  print*,tripoints2
  
end program