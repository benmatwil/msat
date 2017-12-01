!module for the spinefinder
module sf_mod
  use params
  use common

  implicit none

  real(np), parameter :: rsphere = rspherefact*(1e-1_np)**sig_figs

  contains

  !********************************************************************************

  !take spherical coordinates and convert to Cartesian
  function sphere2cart(r, theta, phi)
    implicit none

    real(np) :: theta, phi, r
    real(np) :: sphere2cart(3)

    sphere2cart(1) = r*sin(theta)*cos(phi)
    sphere2cart(2) = r*sin(theta)*sin(phi)
    sphere2cart(3) = r*cos(theta)
  end function

  !********************************************************************************

  subroutine remove_duplicates(vecarray, accur, nclose)
    ! removes vectors with accur distance away from another vector to reduce to "unique" vectors
    implicit none
    
    real(np), allocatable :: vecarray(:,:)
    real(np) :: accur
    integer(int32), allocatable, optional :: nclose(:)
    integer(int32), allocatable :: dummy(:)
    integer(int32) :: i, j, n, nclosei, ny
    
    n = size(vecarray,2)
    i = 1

    if (present(nclose)) allocate(nclose(1))
    
    do while (i < n)
      j = i + 1
      nclosei = 0
      do while (j < n+1)
        if (modulus(vecarray(:, i)-vecarray(:, j)) < accur) then
          call remove_vector(vecarray, j)
          n = size(vecarray, 2)
          nclosei = nclosei + 1
        else
          j = j + 1
        endif
      enddo
      i = i + 1
      if (present(nclose)) call add_element(nclose, nclosei)
    enddo
    
    if (present(nclose)) then
      ny = size(nclose,1)
      dummy = nclose
      deallocate(nclose)
      nclose = dummy(2:ny)
      deallocate(dummy)
    endif

    end
  
  !********************************************************************************

  subroutine it_conv(rnull, rold, rnew, bnew, fact, dir)
    
    implicit none
    
    real(np), dimension(3) :: rnull, rold, rnew, bnew
    real(np) :: fact
    integer(int32) :: dir
    
    rold = rnew
    rnew = rnew + dir*fact*normalise(bnew)
    rnew = rsphere*normalise(rnew)
    bnew = trilinear(rnew+rnull, bgrid)
    
  end
  
end module
