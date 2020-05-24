module SpineFinder

    using ..Params
    using ..Common
    using ..Read
    using ..Trace

    function SF(bgrid::AbstractField3D)

        nnulls, rnulls, rnullsreal = Read.read_nulls(filename)

        signs = zeros(Int32, nnulls)
        spines = zeros(Common.Vector3D, nnulls)
        fans = zeros(Common.Vector3D, nnulls)
        warnings = zeros(Int32, nnulls)

        for inull in 1:nnulls
            println("Evaluating null point $inull of $nnulls")
            @debug "Null position: $(rnulls[inull])"
            @debug "B at null $(Common.trilinear(rnulls[inull], bgrid))"
            @debug "-----------------------------------------------------------------------------"
            signs[inull], spines[inull], fans[inull], warnings[inull] = get_properties(rnulls[inull], bgrid)
        end

        # now write data
        outfname = joinpath(default_output, Read.prefix(bgrid.filename) * "-nulldata.dat")
        spinefile = open(outfname, "w")
            write(spinefile, nnulls)
            write(spinefile, signs, spines, fans, warnings)
        close(spinefile)
      
        for (inull, warning) in enumerate(warnings)
            if warning != 0
                println("Warning on null point $inull with warning $warning")
            end
        end
      
        println("Total number of nulls: $nnulls")
        println("Positive: $(count(signs .== 1))")
        println("Negative: $(count(signs .== -1))")
        println("Unknown: $(count(abs.(signs) .!= 1))")
        println("Warning 1: $(count(warnings .== 1))")
        println("Warning 2: $(count(warnings .== 2))")
        println("Warning 3: $(count(warnings .== 3))")
        println("Warning 4: $(count(warnings .== 4))")

    end

    function get_properties(rnull::Vector3D, bgrid::AbstractField3D)

        rsphere = rspherefact * (1e-1)^sig_figs
        nphi1 = nphi/2
        ntheta1 = ntheta/2

        warning = 0
        warningmessage = ""

        sign_0 = sign_source_sink(rnull, bgrid)

        if abs(sign_0) == 2
            spine = zero(Vector3D)
            spine = zero(Vector3D)
            sign = sign_0
        else
            rconvergefw = Vector{Vector3D}(undef, nphi*ntheta)
            rconvergebw = Vector{Vector3D}(undef, nphi*ntheta)

            # set up theta and phi for sphere around null
            dphi = 2*pi/nphi
            dtheta = pi/ntheta
      
            minmove = 2*pi*rsphere/nphi/50
            mincount = 500
            it_dist = dist_mult*rsphere
            accconv = 1e-9*rsphere
      
            angle = 0
            while true
                global flagfw, flagbw
                flagfw = 0
                flagbw = 0
                # Working on a sphere of size rsphere
                # loop over each point to find converged point
                for itheta in 1:ntheta, iphi in 1:nphi
                    count = 0
                    rfw = sphere2cart(rsphere, (itheta-0.5)*dtheta + angle, (iphi-1)*dphi + angle)
                    rbw = copy(rfw)
                    rfw1 = copy(rfw)
                    rbw1 = copy(rbw)
                    for count in 1:maxiter
                        global oldbw, oldfw
                        oldfw, rfw = it_conv(rfw, rnull, bgrid, it_dist, 1)
                        oldbw, rbw = it_conv(rbw, rnull, bgrid, it_dist, -1)
                        # println(itheta, " ", iphi, " ", count, " ", Common.normalise(rfw))
                        # println(itheta, " ", iphi, " ", count, " ", Common.normalise(rbw))
                        if ((Common.modulus(rfw - oldfw) < accconv) | (Common.modulus(rbw - oldbw) < accconv)) & (count >= mincount)
                            break
                        end
                    end
                    if (Common.modulus(rfw1 - rfw) < minmove) | (Common.modulus(rbw1 - rbw) < minmove) & (angle < dphi)
                        angle = angle + dphi/6
                        @debug "Adjusting the initial points and starting again"
                        @goto end_main
                    end
                    convfw = Common.modulus(rfw - oldfw) < accconv
                    convbw = Common.modulus(rbw - oldbw) < accconv
                    if (convfw & !convbw)
                        flagfw += 1
                    end
                    if (convbw & !convfw)
                        flagbw += 1
                    end
                    iconv = iphi + (itheta-1)*nphi
                    rconvergefw[iconv] = Common.normalise(rfw)
                    rconvergebw[iconv] = Common.normalise(rbw)
                end
                @debug "-----------------------------------------------------------------------------"
                break
                @label end_main
            end
            println("-"^50)
      
            @debug "Number of points which stop convergence"
            @debug "    forwards:  $flagfw"
            @debug "    backwards: $flagbw"
      
            # now to decide on null type and find eigenvectors
            npts = nphi*ntheta
            converge_minvec = false
            if (flagfw == 0) & (flagbw == 0)
                # probably high current - no points converged
                @debug "High current, spine not converging to single point"
                dotprodsfw = abs.(Common.dot.(Ref(rconvergefw[1]), rconvergefw))
                dotprodsbw = abs.(Common.dot.(Ref(rconvergebw[1]), rconvergebw))
        
                if minimum(dotprodsfw) < minimum(dotprodsbw) # spine should be bigger, fan should have a 0
                    sign = 1
                    rspine, rfan = rconvergebw, rconvergefw
                else
                    sign = -1
                    rspine, rfan = rconvergefw, rconvergebw
                end
        
                acc = 1e-3
                rspine, denseposspine = remove_duplicates!(rspine, acc)
                rfan, denseposfan = remove_duplicates!(rfan, acc)
        
                nspine = size(rspine, 1)
                nfan = size(rfan, 1)
        
                if nspine > nfan
                    warning = 4
                end
        
                maxvec = rfan[1]
        
                dotprods = abs.(Common.dot.(Ref(maxvec), rfan))
        
                minvec = rfan[argmin(dotprods)]
        
                fan = normalise(cross(maxvec, minvec))

                spine = normalise(sum(rspine[Common.dot.(Ref(fan), rspine) .> 0]))
      
            else
      
                if flagfw > flagbw
                    sign = -1
                    rspine = rconvergefw
                    rfan = rconvergebw
                    flagspine = flagfw
                    flagfan = flagbw
                elseif flagfw < flagbw
                    sign = 1
                    rspine = rconvergebw
                    rfan = rconvergefw
                    flagspine = flagbw
                    flagfan = flagfw
                else
                    # equal flags - attempt to use eigenvalues...
                    dotprodsfw = abs.(Common.dot.(rconvergefw, Common.trilinear.(rconvergefw .* rsphere .+ Ref(rnull), Ref(bgrid))))
                    dotprodsbw = abs.(Common.dot.(rconvergebw, Common.trilinear.(rconvergebw .* rsphere .+ Ref(rnull), Ref(bgrid))))

                    if maximum(dotprodsfw) > maximum(dotprodsbw)
                        sign = -1
                        rspine = rconvergefw
                        rfan = rconvergebw
                        flagspine = flagfw
                        flagfan = flagbw
                    else
                        sign = 1
                        rspine = rconvergebw
                        rfan = rconvergefw
                        flagspine = flagbw
                        flagfan = flagfw
                    end
                end

                acc = 1e-6
                rspine, denseposspine = remove_duplicates!(rspine, acc)
                rfan, denseposfan = remove_duplicates!(rfan, acc)
                
                nspine = size(rspine, 1)
                nfan = size(rfan, 1)
                
                if (nspine == 2) & (nfan == 2)
                    eigen_spine = Common.dot(rspine[1], Common.trilinear(rspine[1]*rsphere + rnull, bgrid))
                    eigen_fan = Common.dot(rfan[1], Common.trilinear(rfan[1]*rsphere + rnull, bgrid))
                    if abs(eigen_fan) > abs(eigen_spine)
                        rspine, rfan = rfan, rspine
                        sign = -1*sign
                    end
                end
            
                @debug "Number of unique points"
                @debug "    Spine: $nspine"
                @debug "    Fan:   $nfan"
                if nspine <= 2
                    spine = rspine[1]
                else
                    spine = rspine[argmax(denseposspine)]
                end
        
                if nfan == 2
                    maxvec = rfan[1]
                else
                    maxvec = rfan[argmax(denseposfan)]
                end
        
                if (nfan < nspine)
                    if (flagspine > 1.1*flagfan)
                        warning = 2
                    else
                        warning = 3
                        # perhaps need to switch them over
                        spine, maxvec = maxvec, spine
                        rspine, rfan = rfan, rspine
                        sign = -1*sign
            
                        nspine = size(rspine, 2)
                        nfan = size(rfan, 2)
                        @debug "Uh-oh! Swapping!"
                        @debug "    Spine: $nspine"
                        @debug "    Fan:   $nfan"
                        if 3*nspine < nfan
                            warning = 0
                        end
                    end
                end

                dotprods = abs.(Common.dot.(Ref(maxvec), rfan))
        
                if any(dotprods .< 0.1)
                    @debug "Using perpendicular for minvec"
        
                    minvec = rfan[argmin(dotprods)]
                    
                else
                    @debug "Converging minvec"
                    converge_minvec = true
        
                    rmin = Vector{Vector3D}(undef, nphi1*ntheta1)
        
                    dphi = 2*pi/nphi1
                    dtheta = pi/ntheta1
        
                    for itheta in 1:ntheta1, iphi in 1:nphi1
                        rnew = sphere2cart(rsphere, (itheta-0.5)*dtheta + angle, (iphi-1)*dphi + angle)
                        bnew = Common.trilinear(rnew+rnull, bgrid)
                        for count in 1:maxiter/10
                            rold, rnew = it_conv(rnew, rnull, bgrid, maxvec, it_dist, sign)
                            if (Common.modulus(rnew - rold) < accconv) & (count >= mincount)
                                break
                            end
                        end
                        rmin[iphi+(itheta-1)*nphi1] = Common.normalise(rnew)
                    end
                    rmin, densepos = remove_duplicates!(rmin, acc)
                    minvec = rmin[argmax(densepos)]
                    
                    if abs(Common.dot(minvec, maxvec)) > 0.9
                        sign = -1*sign
                        @debug "Min-/max-vec close... switching sign"
                        rspine, rfan = rfan, rspine
                        spine, maxvec = maxvec, spine
                        
                        rmin = Vector{Vector3D}(undef, nphi1*ntheta1)
            
                        dphi = 2*pi/nphi1
                        dtheta = pi/ntheta1
            
                        for itheta in 1:ntheta1, iphi in 1:nphi1
                            rnew = sphere2cart(rsphere, (itheta-0.5)*dtheta + angle, (iphi-1)*dphi + angle)
                            for count = 1:maxiter/10
                                rold, rnew = it_conv(rnew, rnull, bgrid, maxvec, it_dist, sign)
                                if (modulus(rnew - rold) < accconv) & (count >= mincount)
                                    break
                                end
                            end
                            rmin[iphi+(itheta-1)*nphi1] = normalise(rnew)
                        end
                        rmin, densepos = remove_duplicates!(rmin, acc)
                        minvec = rmin[argmax(densepos)]
                    
                    end
        
                end
      
            end
      
            fan = Common.normalise(Common.cross(maxvec, minvec))      

        end

        tilt = abs(90 - rad2deg(acos(Common.dot(fan, spine))))
        println(spine)

        if abs(sign) == 2
            warningmessage = "NULL LIKELY A SOURCE OR SINK"
        else
            if warning == 4
                warningmessage = "Check sign of high current null!"
            elseif warning == 3
                warningmessage = "Convergence of fan and spine similar - switched from normal. CHECK!"
            elseif warning == 2
                warningmessage = "Convergence of fan and spine similar - checks don't agree!"
            elseif tilt < 1e1
                warningmessage = "Spine and fan strongly inclined - high perpendicular J?"
                warning = 1
            end
        end
      
        @debug """
            Final information for null
            --------------------------
            Sign:    $sign
            Warning: $warning $warningmessage
            Spine:   $spine
            Fan:     $fan
            Maxvec:  $maxvec
            Minvec:  $minvec
            Tilt:    $tilt
            Dot products of:
                min and max               max and spine           min and spine'
            dot(minvec, maxvec), dot(maxvec, spine), dot(minvec, spine)
        """

        println("Sign: $sign    Warning: $warning     Tilt: $tilt")
        println("-" ^ 25)

        return sign, spine, fan, warning

    end

    function sphere2cart(r::Float64, theta::Float64, phi::Float64)
        return r * Vector3D([sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)])
    end

    function it_conv(rold::Vector3D, rnull::Vector3D, bgrid::AbstractField3D, it_dist::Float64, dir::Integer)
        
        bvec = Common.trilinear(rold + rnull, bgrid)
        rnew = rold + dir * it_dist * Common.normalise(bvec)
        rnew = rsphere * Common.normalise(rnew)

        return rold, rnew
        
    end

    function it_conv(rold::Vector3D, rnull::Vector3D, bgrid::AbstractField3D, maxvec::Vector3D, it_dist::Float64, dir::Integer)
        
        bvec = Common.trilinear(rold + rnull, bgrid)
        bvec = bvec - Common.dot(bvec, Common.normalise(maxvec))*Common.normalise(maxvec)
        rnew = rold + dir * it_dist * Common.normalise(bvec)
        rnew = rsphere * Common.normalise(rnew)

        return rold, rnew
        
    end

    function remove_duplicates!(vecarray::Vector{Vector3D}, accur::Float64)
        
        n = size(vecarray, 1)
        i = 1

        nclose = Vector{Int32}(undef, 0)

        while i < n
            j = i + 1
            nclosei = 0
            while j < n + 1
                if (Common.modulus(vecarray[i] - vecarray[j]) < accur)
                    deleteat!(vecarray, j)
                    n = size(vecarray, 1)
                    nclosei += 1
                else
                    j += 1
                end
            end
            i += 1
            push!(nclose, nclosei)
        end

        return vecarray, nclose[2:end]

    end

    function sign_source_sink(rnull::Vector3D, bgrid)
        # put a ball of points around null and trace field lines to check for source/sink
    
        # logical :: zero_null
        # integer :: sign_source_sink
        # integer :: iphi, itheta, isink, isink_bw, isink_fw, itrace, dir
        # real(np), dimension(3) :: rtrace, rstart, rnull

        rsphere = rspherefact * (1e-1)^sig_figs
    
        isink = 0
        zero_null = false
        sign_source_sink = 0
        isink_fw = 0
        isink_bw = 0
    
        # put a ball of points around null and trace field lines to check for source/sink
        for iphi in 1:3, itheta in 1:2
            rstart = sphere2cart(rsphere/10, pi/2*(itheta-0.5), 2*pi/3*iphi) + rnull
            for dir in (-1, 1)
                rtrace = rstart
                for itrace = 1:10000
                    rtrace = Trace.trace_line(rtrace, dir, rsphere/1000, bgrid)
                    if (Common.modulus(rtrace - rnull) > rsphere)
                        break
                    end
                end
                if (Common.modulus(rtrace - rnull) < rsphere/10)
                    isink += 1
                end
                if (Common.modulus(rtrace - rnull) < rsphere/10) & (dir == -1)
                    isink_bw += 1
                end
                if (Common.modulus(rtrace - rnull) < rsphere/10) & (dir == 1)
                    isink_fw += 1
                end
            end
        end
    
        if isink >= 5
            zero_null = true
        end
        if isink_fw >= 5
            sign_source_sink = 2
        end
        if isink_bw >= 5
            sign_source_sink = -2
        end

        return sign_source_sink
        
    end

end