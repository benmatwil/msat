module SeparatrixSurfaceFinder

    using OffsetArrays

    using ..Params
    import ..Common
    import ..Trace
    import ..Ring
    import ..Read

    function SSF(filename::AbstractString)
        bgrid = Common.read_field(filename)
        nnulls, rnulls, rnullsreal = Read.read_nulls(filename)
        nnulls, signs, spines, fans = Read.read_nulldata(filename)

        rnullsalt = copy(rnulls)

        for inull in 1:nnulls
            theta = acos(fans[inull][3])
            phi = atan(fans[inull][2], fans[inull][1])
            startpoints = Common.get_startpoints(fans[inull], nstart)

            line1 = Ref(rnulls[inull]) .+ startpoints
            line2 = copy(line1)

            nlines = size(line1, 1)

            breaks = zeros(Int32, nlines)
            nseps = 0
            nperring = OffsetVector(zeros(Int32, ringsmax+1), 0:ringsmax)
            nperring[0] = nstart
            slowdown = 1.0
            terror = 0
            exit_condition = false

            for iring in 1:ringsmax
                println("Ring number: $(iring) with length $(size(line1, 1))")
                if abs(signs[inull]) != 1
                    println("Source/Sink -- skipping")
                    nrings = 1
                    break
                end

                nrings = iring
                endpoints = zeros(Int32, size(line1, 1))
                associations = zeros(Int32, size(line1, 1))

                if iring < 50
                    h0 = stepsize/slowdown/5
                else
                    h0 = stepsize/slowdown
                end

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

                for point in line1
                    Common.edgecheck(point, bgrid)
                end

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
                    nseps = Ring.sep_detect!(line1, rnulls, rnullsalt, signs, spines, breaks, inull, iring, signs[inull], nulldist, bgrid)
                end

                nperring[iring] = size(line1, 1) # maybe iring+1

                # write to file
            end

        end
    end

end