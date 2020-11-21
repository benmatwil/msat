module SeparatrixSurfaceFinder

    using OffsetArrays

    using ..Params
    import ..Common
    import ..Trace
    import ..Ring
    import ..Read

    const hspine = 0.5*stepsize
    const spinemax = 10*ringsmax

    function SSF(bgrid::Common.AbstractField3D)
        # bgrid = Read.read_field(filename, "Cartesian")
        nnulls, rnulls, rnullsreal = Read.read_nulls(bgrid.filename)
        nnulls, signs, spines, fans = Read.read_nulldata(bgrid.filename)

        rnullsalt = Common.get_rnullsalt(rnulls, bgrid)

        outfile_init = joinpath(default_output, Read.prefix(bgrid.filename))

        # stores all info regarding size of rings
        ringinfo_file = open(outfile_init * "-ringinfo.dat", "w+")
        # stores all the coordinates of rings in original coordinate system
        ring_file = open(outfile_init * "-rings.dat", "w+")
        # stores all the associations between rings
        if assoc_output
            assoc_file = open(outfile_init * "-assocs.dat", "w+")
        end
        # stores all the breaks on rings
        break_file = open(outfile_init * "-breaks.dat", "w+")
        # stores all the connection info about each separator
        connectivity_file = open(outfile_init * "-connectivity.dat", "w+")
        # stores all the coordinates of the separator lines in original coordinate system
        separator_file = open(outfile_init * "-separators.dat", "w+")
        # stores all the coordinates of the spine lines in original coordinate system
        spine_file = open(outfile_init * "-spines.dat", "w+")

        write(ringinfo_file, Int32(ringsmax+1), Int32(nskip), Int32(sizeof(bytesize)), stepsize)

        tempfilename = outfile_init * "-rings.temp"
        uptonullconn = 0

        for inull in 1:nnulls
            tempfile = open(tempfilename, "w+")

            startpoints = Common.get_startpoints(fans[inull], nstart)

            line1 = Ref(rnulls[inull]) .+ startpoints
            line2 = copy(line1)

            nlines = size(line1, 1)

            breaks = zeros(Int32, nlines)
            associations = zeros(Int32, nlines)
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
            write(connectivity_file, Int32(nseps))

            nrings = 0

            for iring in 1:ringsmax
                
                if abs(signs[inull]) != 1
                    println("Source/Sink -- skipping")
                    nrings = 1
                    break
                end

                if iring < 50
                    h0 = stepsize/slowdown/5
                else
                    h0 = stepsize/slowdown
                end

                r = Trace.trace_line.(line1, signs[inull], h0, Ref(bgrid))

                line2 .= line2 .+ (r .- line1)
                line1 .= r

                associations .= 1:length(line1)
                endpoints = Common.outedge.(line1, Ref(bgrid))

                if iring < 50
                    maxdist = 0.1*h0*slowdown
                elseif iring < 100
                    maxdist = 0.08*h0*slowdown
                else
                    maxdist = 0.5*h0*slowdown
                end
                nulldist = 2.0*h0*slowdown
                mindist = maxdist/3
                slowdist = 2.0*nulldist

                line1, line2, breaks, associations = Ring.remove_points(line1, line2, breaks, associations, endpoints, mindist)

                nearnull = zeros(Int32, size(line1, 1))
                slowdown = 1.0

                # check for lines close to null points
                for iline in 1:size(line1, 1)
                    for jnull in 1:nnulls
                        if signs[jnull] == signs[inull]
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

                line1 = Common.edgecheck.(line1, Ref(bgrid))

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
                    nseps += Ring.sep_detect!(line1, rnulls, rnullsalt, signs, spines, breaks, inull, iring, signs[inull], nulldist, bgrid, connectivity_file)
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
                        rsep[iring+2] = read(tempfile, Common.Vector3D{Float64})
                    end
                    rsep[1] = rnullsreal[nullnum1]
                    rsep[ringnum+3] = rnullsreal[nullnum2]
                    write(separator_file, Int32(ringnum+3))
                    write.(Ref(separator_file), rsep)
                end
            end
            uptonullconn = uptonullconn + nseps*16 + 4

            # trace spine lines and write file
            for dir in (-1, 1)
                out = false
                spine = [rnulls[inull]]
                rspine = rnulls[inull] + dir*spines[inull]*1e-3 # pick a good factor
                iring = 0
                while ! Common.outedge(rspine, bgrid) & (iring < spinemax)
                    iring = iring + 1
                    if abs(signs[inull]) != 1
                        break
                    end
                    push!(spine, rspine)
                    rspine = Trace.trace_line(rspine, -signs[inull], hspine, bgrid)
                    rspine = Common.edgecheck(rspine, bgrid)
                end
                write(spine_file, Int32(size(spine, 1)), Common.gtr.(spine, Ref(bgrid)))
            end

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

        close(ringinfo_file)
        close(ring_file)
        if assoc_output close(assoc_file) end
        close(break_file)
        close(connectivity_file)
        close(separator_file)
        close(spine_file)

    end

end