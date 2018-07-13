!module for the spinefinder
module sf_mod
  use params
  use common
  use trace

  implicit none

  real(np), parameter :: rsphere = rspherefact*(1e-1_np)**sig_figs
  logical :: savedata = .false.

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
    
    real(np), allocatable :: vecarray(:, :)
    real(np) :: accur
    integer(int32), allocatable, optional :: nclose(:)
    integer(int32), allocatable :: dummy(:)
    integer(int32) :: i, j, n, nclosei, ny, count_before, k
    logical, dimension(:), allocatable :: keep_point
    integer(int32), dimension(:), allocatable :: remove_pt
    
    n = size(vecarray,2)
    i = 1

    if (present(nclose)) allocate(nclose(0))
    
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
      ny = size(nclose, 1)
      dummy = nclose
      deallocate(nclose)
      nclose = dummy(2:ny)
      deallocate(dummy)
    endif

    ! n = size(vecarray,2)
    ! allocate(keep_point(n), remove_pt(n))
    ! keep_point = .true.
    ! remove_pt = 0

    ! k = 1
    ! do i = 1, n
    !   do j = i+1, n
    !     if (keep_point(j)) then
    !       if (modulus(vecarray(:, i)-vecarray(:, j)) < accur) then
    !         keep_point(j) = .false.
    !         remove_pt = k
    !       endif
    !     endif
    !   enddo
    !   if (count(keep_point(:)) < count_before) then
    !     k = k + 1
    !     count_before = count(keep_point(:))
    !   endif
    ! enddo

    ! n = count(keep_point(:))
    ! vecarray = reshape(pack(vecarray, spread(keep_point, 1, 3)), [3, n])
    ! if (present(nclose)) then
    !   allocate(nclose(n))
    !   do i = 1, n
    !     nclose(i) = count(remove_pt == i)
    !   enddo
    ! endif

    end
  
  !********************************************************************************

  subroutine it_conv(rnull, rold, rnew, bnew, it_dist, dir)
    
    implicit none
    
    real(np), dimension(3) :: rnull, rold, rnew, bnew
    real(np) :: it_dist
    integer(int32) :: dir
    
    rold = rnew
    rnew = rnew + dir*it_dist*normalise(bnew)
    rnew = rsphere*normalise(rnew)
    bnew = trilinear(rnew+rnull, bgrid)
    
  end

  !********************************************************************************

  function sign_source_sink(rnull)
    ! put a ball of points around null and trace field lines to check for source/sink

    logical :: zero_null
    integer :: sign_source_sink
    integer :: iphi, itheta, isink, isink_bw, isink_fw, itrace, dir
    real(np), dimension(3) :: rtrace, rstart, rnull

    isink = 0
    zero_null = .false.
    sign_source_sink = 0
    isink_fw = 0
    isink_bw = 0

    ! put a ball of points around null and trace field lines to check for source/sink
    do iphi = 1, 3
      do itheta = 1, 2
        rstart = sphere2cart(rsphere/10, pi/2*(itheta-0.5_np), 2*pi/3*iphi) + rnull
        do dir = -1, 1, 2
          rtrace = rstart
          do itrace = 1, 10000
            call trace_line(rtrace, dir, rsphere/1000)
            if (modulus(rtrace - rnull) > rsphere) exit
          enddo
          if (modulus(rtrace - rnull) < rsphere/10) isink = isink + 1
          if (modulus(rtrace - rnull) < rsphere/10 .and. dir == -1) isink_bw = isink_bw + 1
          if (modulus(rtrace - rnull) < rsphere/10 .and. dir == 1) isink_fw = isink_fw + 1
        enddo
      enddo
    enddo

    if (isink >= 5) zero_null = .true.
    if (isink_fw >= 5) sign_source_sink = 2
    if (isink_bw >= 5) sign_source_sink = -2
    print*, sign_source_sink, isink_fw, isink_bw
    
  end function

  !********************************************************************************

  
  subroutine get_properties(inull, sign, spine, fan, warning)
  ! characterise each null

    implicit none
    integer(int32) :: itheta, iphi, ifan, ispine
    integer(int32) :: count, iconv, imin, mincount, npts
    integer(int32) :: flagfw, flagbw, flagspine, flagfan
    integer(int32) :: nfw, nbw, nspine, nfan

    logical :: converge_minvec

    integer(int32), parameter :: nphi1 = nphi/2, ntheta1 = ntheta/2
    real(np) :: dphi, dtheta, angle, minmove
    real(np), dimension(3) :: roldfw, rnewfw, bnewfw, roldbw, rnewbw, bnewbw
    real(np), dimension(3) :: rold, rnew, bnew, rfw1, rbw1
    real(np) :: it_dist, acc, accconv
    logical :: convfw, convbw

    real(np), dimension(:,:), allocatable :: rconvergefw, rconvergebw, rmin
    real(np), dimension(:,:), allocatable :: rspine, rfan
    integer(int32), dimension(:), allocatable :: densepos, denseposfw, denseposbw, denseposspine, denseposfan
    real(np), dimension(3) :: av_spine, swap
    real(np), dimension(:), allocatable :: dotprods, dotprodsfw, dotprodsbw
    real(np) :: eigen_spine, eigen_fan
    logical :: dotprodfw_big

    integer(int32) :: inull
    real(np), dimension(3) :: rnull, spine, fan, maxvec, minvec
    integer(int32) :: sign, warning, sign_0
    character(:), allocatable :: warningmessage

#if debug
    print*, 'Evaluating null', inull,' of', nnulls
    print*, 'Null position:'
    print*, rnulls(:, inull)
    print*, 'B at null'
    print*, trilinear(rnulls(:, inull), bgrid)
    print*, '-----------------------------------------------------------------------------'
#endif

    rnull = rnulls(:, inull)

    warning = 0
    warningmessage = ''

    sign_0 = sign_source_sink(rnull)
    
    if (abs(sign_0) == 2) then
      spine = 0
      fan = 0
      sign = sign_0
    else
      allocate(rconvergefw(3, nphi*ntheta), rconvergebw(3, nphi*ntheta))

      !set up theta and phi for sphere around null
      dphi = 2*pi/nphi
      dtheta = pi/ntheta

      minmove = 2*pi*rsphere/nphi/50
      mincount = 500
      it_dist = dist_mult*rsphere
      accconv = 1e-9_np*rsphere

      angle = 0
      main: do
        
        flagfw = 0
        flagbw = 0
        ! Working on a sphere of size rsphere
        ! loop over each point to find converged point
        do itheta = 1, ntheta
          do iphi = 1, nphi
            count = 0
            rnewfw = sphere2cart(rsphere, (itheta-0.5_np)*dtheta + angle, (iphi-1)*dphi + angle)
            rnewbw = rnewfw
            rfw1 = rnewfw
            rbw1 = rnewbw
            bnewfw = trilinear(rnewfw+rnull, bgrid)
            bnewbw = bnewfw
            roldfw = [0, 0, 0]
            roldbw = [0, 0, 0]
            do count = 1, maxiter
              call it_conv(rnull, roldfw, rnewfw, bnewfw, it_dist, 1)
              call it_conv(rnull, roldbw, rnewbw, bnewbw, it_dist, -1)
              if ((modulus(rnewfw-roldfw) < accconv .or. modulus(rnewbw-roldbw) < accconv) .and. count >= mincount) exit
            enddo
            if ((modulus(rfw1-rnewfw) < minmove .or. modulus(rbw1-rnewbw) < minmove) .and. angle < dphi) then
              angle = angle + dphi/6
#if debug
              print*, 'Adjusting the initial points and starting again'
#endif
              cycle main
            endif
            convfw = modulus(rnewfw - roldfw) < accconv
            convbw = modulus(rnewbw - roldbw) < accconv
            if (convfw .and. .not. convbw) flagfw = flagfw + 1
            if (convbw .and. .not. convfw) flagbw = flagbw + 1
            iconv = iphi + (itheta-1)*nphi
            rconvergefw(:, iconv) = normalise(rnewfw)
            rconvergebw(:, iconv) = normalise(rnewbw)
          enddo
        enddo
#if debug
        print*, '-----------------------------------------------------------------------------'
#endif
        exit
      enddo main

#if debug
      print*, 'Number of points which stop convergence'
      print*, '    ', 'forwards: ', flagfw
      print*, '    ', 'backwards:', flagbw
#endif

      ! now to decide on null type and find eigenvectors
      npts = nphi*ntheta
      converge_minvec = .false.
      if (flagfw == 0 .and. flagbw == 0) then
        ! probably high current - no points converged
#if debug
        print*, 'High current, spine not converging to single point'
#endif
        allocate(dotprodsfw(npts), dotprodsbw(npts))
        do ifan = 1, npts
          dotprodsfw(ifan) = abs(dot(rconvergefw(:, 1), rconvergefw(:, ifan)))
          dotprodsbw(ifan) = abs(dot(rconvergebw(:, 1), rconvergebw(:, ifan)))
        enddo
        dotprodfw_big = minval(dotprodsfw) > minval(dotprodsbw) ! spine should be bigger, fan should have a 0

        if (minval(dotprodsfw) < minval(dotprodsbw)) then
          sign = 1
          call move_alloc(rconvergebw, rspine)
          call move_alloc(rconvergefw, rfan)
        else
          sign = -1
          call move_alloc(rconvergefw, rspine)
          call move_alloc(rconvergebw, rfan)
        endif

        deallocate(dotprodsfw, dotprodsbw)

        acc = 1e-3_np
        call remove_duplicates(rspine, acc, denseposspine)
        call remove_duplicates(rfan, acc, denseposfan)

        nspine = size(rspine, 2)
        nfan = size(rfan, 2)

        if (nspine > nfan) warning = 4

        maxvec = rfan(:, 1)

        allocate(dotprods(nfan))

        do ifan = 1, nfan
          dotprods(ifan) = abs(dot(maxvec, rfan(:, ifan)))
        enddo

        minvec = rfan(:, maxval(minloc(dotprods)))

        fan = normalise(cross(maxvec, minvec))
        
        av_spine = [0, 0, 0]
        do ispine = 1, size(rspine, 2)
          if (dot(fan, rspine(:, ispine)) > 0) av_spine = av_spine + rspine(:, ispine)
        enddo
        spine = normalise(av_spine)

      else

        if (flagfw > flagbw) then
          sign = -1
          call move_alloc(rconvergefw, rspine)
          call move_alloc(rconvergebw, rfan)
          flagspine = flagfw
          flagfan = flagbw
        elseif (flagfw < flagbw) then
          sign = 1
          call move_alloc(rconvergebw, rspine)
          call move_alloc(rconvergefw, rfan)
          flagspine = flagbw
          flagfan = flagfw
        else
          ! equal flags - attempt to use eigenvalues...
          allocate(dotprodsfw(npts), dotprodsbw(npts))
          do ifan = 1, npts
            dotprodsfw(ifan) = &
              abs(dot(rconvergefw(:, ifan), trilinear(rconvergefw(:, ifan)*rsphere + rnull, bgrid)))
            dotprodsbw(ifan) = &
              abs(dot(rconvergebw(:, ifan), trilinear(rconvergebw(:, ifan)*rsphere + rnull, bgrid)))
          enddo
          if (maxval(dotprodsfw) > maxval(dotprodsbw)) then
            sign = -1
            call move_alloc(rconvergefw, rspine)
            call move_alloc(rconvergebw, rfan)
            flagspine = flagfw
            flagfan = flagbw
          else
            sign = 1
            call move_alloc(rconvergebw, rspine)
            call move_alloc(rconvergefw, rfan)
            flagspine = flagbw
            flagfan = flagfw
          endif
        endif

        acc = 1e-6_np
        call remove_duplicates(rspine, acc, denseposspine)
        call remove_duplicates(rfan, acc, denseposfan)
        
        nspine = size(rspine, 2)
        nfan = size(rfan, 2)
        
        if (nspine == 2 .and. nfan == 2) then
          eigen_spine = dot(rspine(:, 1), trilinear(rspine(:, 1)*rsphere + rnull, bgrid))
          eigen_fan = dot(rfan(:, 1), trilinear(rfan(:, 1)*rsphere + rnull, bgrid))
          if (abs(eigen_fan) > abs(eigen_spine)) then
            rconvergefw = rspine
            call move_alloc(rfan, rspine)
            call move_alloc(rconvergefw, rfan)
            sign = -1*sign
          endif
        endif
        
#if debug
        print*, 'Number of unique points'
        print*, '    Spine:', nspine
        print*, '    Fan:  ', nfan
#endif
! print*, size(rspine)
        if (nspine <= 2) then
          spine = rspine(:, 1)
        else
          spine = rspine(:, maxval(maxloc(denseposspine)))
        endif

        if (nfan == 2) then
          maxvec = rfan(:, 1)
        else
          maxvec = rfan(:, maxval(maxloc(denseposfan)))
        endif

        if (nfan < nspine) then
          if (flagspine > 1.1*flagfan) then
            warning = 2
          else
            warning = 3
            ! perhaps need to switch them over
            swap = maxvec
            maxvec = spine
            spine = swap
            call move_alloc(rspine, rconvergefw)
            call move_alloc(rfan, rspine)
            call move_alloc(rconvergefw, rfan)
            sign = -1*sign

            nspine = size(rspine, 2)
            nfan = size(rfan, 2)
#if debug
        print*, 'Uh-oh! Swapping!'
        print*, '    Spine:', nspine
        print*, '    Fan:  ', nfan
#endif
            if (3*nspine < nfan) then
              warning = 0
            endif
          endif
        endif

        allocate(dotprods(nfan))

        do ifan = 1, nfan
          dotprods(ifan) = abs(dot(maxvec, rfan(:, ifan)))
        enddo

        if (any(dotprods < 0.1)) then
#if debug
        print*, 'Using perpendicular for minvec'
#endif

          minvec = rfan(:, maxval(minloc(dotprods)))
          
        else
#if debug
        print*, 'Converging minvec'
#endif
          converge_minvec = .true.

          allocate(rmin(3, nphi1*ntheta1))

          dphi = 2*pi/nphi1
          dtheta = pi/ntheta1

          do itheta = 1, ntheta1
            do iphi = 1, nphi1
              rnew = sphere2cart(rsphere, (itheta-0.5_np)*dtheta + angle, (iphi-1)*dphi + angle)
              bnew = trilinear(rnew+rnull, bgrid)
              rold = [0,0,0]
              do count = 1, maxiter/10
                bnew = bnew - dot(bnew, normalise(maxvec))*normalise(maxvec)
                call it_conv(rnull, rold, rnew, bnew, it_dist, sign)
                if (modulus(rnew-rold) < accconv .and. count >= mincount) exit
              enddo
              rmin(:,iphi+(itheta-1)*nphi1) = normalise(rnew)
            enddo
          enddo
          call remove_duplicates(rmin, acc, densepos)
          minvec = rmin(:, maxval(maxloc(densepos)))
          
          if (abs(dot(minvec, maxvec)) > 0.9) then
            sign = -1*sign
            rconvergefw = rfan
#if debug
            print*, 'Min-/max-vec close... switching sign'
#endif
            call move_alloc(rspine, rfan)
            call move_alloc(rconvergefw, rspine)
            
            swap = spine
            spine = maxvec
            maxvec = swap
            deallocate(rmin, densepos)
            allocate(rmin(3, nphi1*ntheta1))

            dphi = 2*pi/nphi1
            dtheta = pi/ntheta1

            do itheta = 1, ntheta1
              do iphi = 1, nphi1
                rnew = sphere2cart(rsphere, (itheta-0.5_np)*dtheta + angle, (iphi-1)*dphi + angle)
                bnew = trilinear(rnew+rnull, bgrid)
                rold = [0,0,0]
                do count = 1, maxiter/10
                  bnew = bnew - dot(bnew, normalise(maxvec))*normalise(maxvec)
                  call it_conv(rnull, rold, rnew, bnew, it_dist, sign)
                  if (modulus(rnew-rold) < accconv .and. count >= mincount) exit
                enddo
                rmin(:,iphi+(itheta-1)*nphi1) = normalise(rnew)
              enddo
            enddo
            call remove_duplicates(rmin, acc, densepos)
            minvec = rmin(:, maxval(maxloc(densepos)))
            
          endif

        endif

      endif

      fan = normalise(cross(maxvec, minvec))

    endif

    ! Save data if there is a problematic null for inspection
    if (savedata) then
      if (sign /= 0) then
        open(unit=30, file='output/spine.dat', access='stream', status='replace')
          write(30) size(rspine, 2), rspine, spine
        close(30)
        open(unit=31, file='output/maxvec.dat', access='stream', status='replace')
          write(31) size(rfan, 2), rfan, maxvec
        close(31)
        if (converge_minvec) then
          open(unit=32, file='output/minvec.dat', access='stream', status='replace')
            write(32) size(rmin, 2), rmin, minvec
          close(32)
        else
          open(unit=32, file='output/minvec.dat', access='stream', status='replace')
            write(32) 0, minvec
          close(32)
        endif
      else
        open(unit=30, file='output/spine.dat', access='stream', status='replace')
          write(30) 0, spine
        close(30)
        open(unit=31, file='output/maxvec.dat', access='stream', status='replace')
          write(31) 0, maxvec
        close(31)
        open(unit=32, file='output/minvec.dat', access='stream', status='replace')
          write(32) 0, minvec
        close(32)
      endif
    endif

    if (abs(sign) == 2) then
      warningmessage = 'NULL LIKELY A SOURCE OR SINK'
    else
      if (warning == 4) then
        warningmessage = "Check sign of high current null!"
      elseif (warning == 3) then
        warningmessage = "Convergence of fan and spine similar - switched from normal. CHECK!"
      elseif (warning == 2) then
        warningmessage = "Convergence of fan and spine similar - checks don't agree!"
      elseif (abs(90 - acos(dot(fan, spine))/dtor) < 1e1_np) then
        warningmessage = 'Spine and fan strongly inclined - high perpendicular J?'
        warning = 1
      endif
    endif

#if debug
    print*, ''
    print*, 'Final information for null'
    print*, '--------------------------'
    print*, 'Sign:   ', sign
    print*, 'Warning:', warning, trim(warningmessage)
    print*, 'Spine:  ', spine
    print*, 'Fan:    ', fan
    print*, 'Maxvec: ', maxvec
    print*, 'Minvec: ', minvec
    print*, 'Tilt:   ', abs(90 - acos(dot(fan,spine))/dtor)
    print*, 'Dot products of:'
    print*, '  min and max             ', '  max and spine           ', '  min and spine'
    print*, dot(minvec, maxvec), dot(maxvec, spine), dot(minvec, spine)
#else
    print('(A, I5, A, I5)'), ' Null', inull, ' of', nnulls
    print('(A, I5, A, I2, A, F6.2)'), ' Sign:', sign, '    Warning:', warning, '     Tilt:', abs(90 - acos(dot(fan, spine))/dtor)
    print*, '-----------------------------------------'
#endif

  end subroutine

  
end module
