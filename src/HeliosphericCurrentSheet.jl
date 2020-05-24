module HeliosphericCurrentSheets

    using OffsetArrays
    using StaticArrays

    using ..Params
    import ..Common
    import ..Trace
    import ..Ring
    import ..Read

    function HCS(bgrid::Common.SphericalField3D)
        # bgrid = Read.read_field(filename, "Cartesian")
        nnulls, rnulls, rnullsreal = Read.read_nulls(bgrid.filename)
        nnulls, signs, spines, fans = Read.read_nulldata(bgrid.filename)

        rnullsalt = copy(rnulls)

        outfile_init = joinpath(default_output, Read.prefix(bgrid.filename))

        # stores all info regarding size of rings
        ringinfo_file = open(outfile_init * "-hcs-ringinfo.dat", "w+")
        # stores all the coordinates of rings in original coordinate system
        ring_file = open(outfile_init * "-hcs-rings.dat", "w+")
        # stores all the associations between rings
        if assoc_output
            assoc_file = open(outfile_init * "-hcs-assocs.dat", "w+")
        end
        # stores all the breaks on rings
        break_file = open(outfile_init * "-hcs-breaks.dat", "w+")
        # stores all the connection info about each separator
        connectivity_file = open(outfile_init * "-hcs-connectivity.dat", "w+")
        # stores all the coordinates of the separator lines in original coordinate system
        separator_file = open(outfile_init * "-hcs-separators.dat", "w+")

        hcs_components = get_hcs_lines(bgrid, 10)
        ncomps = size(hcs_components, 1)

        write(ringinfo_file, Int32(ringsmax+1), Int32(nskip), Int32(sizeof(bytesize)), stepsize, Int32(ncomps))
        write(separator_file, Int32(0))

        @info "There are $ncomps components of the hcs"

        tempfilename = outfile_init * "-hcs-rings.temp"
        uptonullconn = 0

        for ihcs in 1:ncomps
            @info "Component $ihcs of the HCS"
            for dir in (1, -1)
                if dir == 1
                    @info "Tracing forwards"
                else
                    @info "Tracing backwards"
                end
                tempfile = open(tempfilename, "w+")

                line1 = copy(hcs_components[ihcs])
                line2 = copy(line1)

                nlines = size(line1, 1)

                breaks = zeros(Int32, nlines)
                nseps = 0
                nperring = OffsetVector(zeros(Int32, ringsmax+1), 0:ringsmax)
                nperring[0] = nstart
                slowdown = 1.0
                terror = 0
                exit_condition = false

                # write initial ring to file
                write_ring = Common.gtr.(line1, Ref(bgrid))
                write(ring_file, write_ring)
                if (assoc_output)
                    write(assoc_file, Int32.(1:nlines))
                end
                write(break_file, breaks)
                write(tempfile, Int32.(1:nlines), write_ring)

                nrings = 0

                for iring in 1:ringsmax
                    if iring % 20 == 0
                        println("Ring number: $(iring) with length $(size(line1, 1))")
                    end

                    endpoints = zeros(Int32, size(line1, 1))
                    associations = zeros(Int32, size(line1, 1))

                    h0 = stepsize/slowdown

                    for (iline, r) in enumerate(line1)
                        # r = line1[iline]
                        h = h0

                        r = Trace.trace_line(r, signs[inull], h, bgrid)

                        line2[iline] = line2[iline] + (r - line1[iline])
                        line1[iline] = r

                        if Common.outedge(r, bgrid)
                            endpoints[iline] = 1
                        else
                            endpoints[iline] = 0
                        end

                        associations[iline] = iline
                    end

                    maxdist = 0.5*h0*slowdown
                    nulldist = 2.0*h0*slowdown
                    mindist = maxdist/3
                    slowdist = 2.0*nulldist

                    line1, line2, breaks, associations = Ring.remove_points(line1, line2, breaks, associations, endpoints, mindist)

                    nearnull = zeros(Int32, size(line1, 1))
                    slowdown = 1.0

                    # check for lines close to null points
                    for iline in 1:size(line1, 1)
                        for jnull in 1:nnulls
                            if signs[jnull] == dir
                                continue
                            end
                            if (Common.dist(rnulls[jnull], line1[iline]) < slowdist) | (Common.dist(rnullsalt[jnull], line1[iline]) < slowdist)
                                nearnull[iline] = 1
                                break
                            end
                        end
                    end

                    nearflag = sum(nearnull)
                    if nearflag > 0
                        slowdown = 2.0
                    end

                    Ring.add_points!(line1, line2, breaks, associations, nearnull, maxdist)

                    line1 = Common.edgecheck.(point, Ref(bgrid))

                    if size(line1, 1) > pointsmax
                        println("Reached maximum number of points")
                        exit_condition = true
                    elseif size(line1, 1) == 0
                        println("All points have left the domain")
                        exit_condition = true
                    elseif iring - samemax - 1 > 0
                        if all(nperring[iring - samemax - 1:iring-1] .== nperring[iring-1]) # maybe +1
                            println("Ring seems to have stopped growing, maybe an unfound null point: $(iring)")
                            exit_condition = true
                        end
                    end
                    
                    nrings = iring

                    if exit_condition
                        println("Exiting")
                        break
                    end

                    if nearflag > 0
                        nseps += Ring.sep_detect!(line1, rnulls, rnullsalt, signs, spines, breaks, 0, iring, dir, nulldist, bgrid, connectivity_file)
                    end

                    nperring[iring] = size(line1, 1)

                    write_ring = Common.gtr.(line1, Ref(bgrid))
                    if iring % nskip == 0
                        write(ring_file, write_ring)
                        write(break_file, breaks)
                    end
                    write(tempfile, associations, write_ring)

                end

                flush(tempfile)

                write(ringinfo_file, nperring[0:nskip:end])

                # trace separators and write to file
                println("Tracing spines and any separators and dealing with rings")
                seek(connectivity_file, uptonullconn)
                write(connectivity_file, Int32(nseps))
                if nseps > 0
                    for isep in 1:nseps
                        nullnum1, nullnum2, ringnum, linenum = read!(connectivity_file, Vector{Int32}(undef, 4)) 
                        rsep = Vector{Common.Vector3D}(undef, ringnum + 3)
                        for iring in ringnum:-1:0
                            ia, ip = Common.file_position(iring, nperring, linenum)
                            seek(tempfile, ia)
                            linenum = read(tempfile, Int32)
                            seek(tempfile, ip)
                            rsep[iring+2] = read(tempfile, Common.Vector3D)
                        end
                        rsep[1] = rnullsreal[nullnum1]
                        rsep[ringnum+3] = rnullsreal[nullnum2]
                        write(separator_file, Int32(ringnum+3), rsep)
                    end
                end
                uptonullconn = uptonullconn + nseps*16 + 4

                # write the associations to file corrected for nskip != 1
                if assoc_output
                    for iring in nskip:nskip:nskip * (nrings ÷ nskip)
                        ia, ip = Common.file_position(iring, nperring, 1)
                        # println("$ia $ip $iring $nskip $nrings $(stat(tempfile).size)")
                        seek(tempfile, ia)
                        association = read!(tempfile, Vector{Int32}(undef, nperring[iring]))
                        for iskip in 1:nskip-1
                            ia, ip = Common.file_position(iring-iskip, nperring, 1)
                            # println("In loop $ia $ip $iring $nskip $nrings $(stat(tempfile).size)")
                            seek(tempfile, ia)
                            assoc_temp = read!(tempfile, Vector{Int32}(undef, nperring[iring-iskip]))
                            association = assoc_temp[association]
                        end
                        write(assoc_file, association)
                    end
                end

                if nrings == ringsmax
                    println("Reached maximum number of rings")
                end

                println("Number of separators: $nseps\nNumber of rings: $nrings")

                # close temporary ring file and delete it
                close(tempfile)

            end

        end

        close(ringinfo_file)
        close(ring_file)
        close(assoc_file)
        close(break_file)
        close(connectivity_file)
        close(separator_file)

    end


    function get_hcs_lines(field::Common.SphericalField3D, nsplit::Integer)

        # number of theta and phi required
        nnt = nsplit * (size(field.y, 1) - 1) + 1
        nnp = nsplit * (size(field.z, 1) - 1) + 1

        tgrid = 0:1/nsplit:size(field.y, 1)
        pgrid = 0:1/nsplit:size(field.z, 1)
        tgrid[1] = tgrid[1] + 1e-6
        tgrid[end] = tgrid[nnt] - 1e-6

        points = Vector{Common.Vector3D}(undef, 0)

        rr = Float64(size(field.x, 1))
        rr1 = rr - 1e-6
        
        for theta in tgrid
            brs = [Common.trilinear(SA[rr, theta, phi], field)[1] for phi in pgrid]
            idiffs = findall(diff(brs) .< 0)
            rps = brs[idiffs]/(brs[idiffs] - brs[idiffs+1])*(pgrid[idiffs+1] - pgrid[idiffs]) + pgrid[idiffs]
            append!(points, [SA[rr1, theta, rps] for rp in rps])
        end

        for phi in pgrid
            brs = [Common.trilinear(SA[rr, theta, phi], field)[1] for theta in tgrid]
            idiffs = findall(diff(brs) .< 0)
            rts = brs[idiffs]/(brs[idiffs] - brs[idiffs+1])*(tgrid[idiffs+1] - tgrid[idiffs]) + tgrid[idiffs]
            append!(points, [SA[rr1, rt, phi] for rt in rts])
        end

        # sort points
        ds = 0.05 # need to think about this
        # allocate(pordered(3, 0))
        pordered = Vector{Vector{Common.Vector3D}}(undef, 0)
      
        while (size(points, 2) > 0)
            ind_line = [points[1]]
            deleteat!(points, 1)
            for dir in (1, 2)
                while size(points, 1) > 0
                    dists = sum.(broadcast.(^, points .- Ref(ind_line[end]), 2))
                    imin = argmin(dists)
                    if dists[imin] < ds
                        push!(ind_line, points[imin])
                        deleteat!(points, imin)
                    else
                        break
                    end
                end
                reverse!(ind_line)
            end
            push!(pordered, ind_line)
        end

        return pordered
      
    end

end