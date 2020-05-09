module Ring

    using StaticArrays

    using ..Params
    using ..Common
    using ..Trace

    function add_points!(line1::Vector{Vector3D},
                         line2::Vector{Vector3D},
                         breaks::Vector{Int32},
                         associations::Vector{Int32},
                         nearnull::Vector{Int32},
                         maxdist::Float64)

        nlines = size(line1, 1)
        add1 = zero(line1)
        add2 = zero(line2)

        for iline in 1:nlines
            maxdist0 = maxdist
            if breaks[iline] == 0
                if iline < nlines # if not the last point
                    nxtline = iline+1
                else
                    nxtline = 1
                end
                if (nearnull[iline] == 1) | (nearnull[nxtline] == 1)
                    maxdist0 = maxdist/2
                end
                if Common.dist(line2[iline], line2[nxtline]) > maxdist0 # if two adjacent points too far away
                    # add point half way between two points
                    add1[iline] = line1[iline] + 0.5 * (line2[nxtline] - line2[iline])
                    add2[iline] = line2[iline] + 0.5 * (line2[nxtline] - line2[iline])
                end
            end
        end

        # add in reverse to make indexing simple
        for iline in nlines:-1:1
            if add1[iline][1] > 1e-1
                iadd = iline + 1
                insert!(line1, iadd, add1[iline])
                insert!(line2, iadd, add2[iline])
                insert!(breaks, iadd, 0)
                if iline != nlines
                    insert!(associations, iadd, associations[iline])
                else
                    insert!(associations, iadd, 1)
                end
            end
        end

        return line1, line2, breaks, associations

    end
    
    function remove_points(line1::Vector{Vector3D},
                           line2::Vector{Vector3D},
                           breaks::Vector{Int32},
                           associations::Vector{Int32},
                           endpoints::Vector{Int32},
                           mindist::Float64)
    
        nlines = size(line1, 1)
        remove = zeros(Int32, nlines)
    
        break1 = breaks[1]
        for iline in 2:nlines
            if (endpoints[iline] == 1) | (breaks[iline] == 2) # remove points that have left the simulation or stuck
                if breaks[iline] == 2
                    endpoints[iline] = 1
                end
                remove[iline] = 1
                breaks[iline-1] = 1
            end
        end
        if (endpoints[1] == 1) | (break1 == 2) # remove points that have left the simulation or stuck
            if breaks[1] == 2
                endpoints[1] = 1
            end
            remove[1] = 1
            breaks[nlines] = 1
        end
    
        # check for too tightly spaced points, flag points to be removed
        for iline in 1:nlines # loop over all points
            if breaks[iline] == 0
                if iline != nlines # if not end point
                    nxtline = iline+1
                else
                    nxtline = 1
                end
                if (Common.dist(line2[iline], line2[nxtline]) < mindist) & (remove[iline] == 0) # if iline and iline+1 are too near
                    if breaks[nxtline] == 1
                        remove[iline] = 1
                    else
                        remove[nxtline] = 1
                    end
                end
            end
        end
    
        if (nlines > nstart) | (sum(endpoints) != 0)
            left_overs = remove .== 0
            line1 = line1[left_overs]
            line2 = line2[left_overs]
            breaks = breaks[left_overs]
            associations = associations[left_overs]
        end

        return line1, line2, breaks, associations

    end

    function sep_detect!(line1::Vector{Vector3D},
                        rnulls::Vector{Vector3D},
                        rnullsalt::Vector{Vector3D},
                        signs::Vector{Int32},
                        spines::Vector{Vector3D},
                        breaks::Vector{Int32},
                        nullnum::Integer,
                        nring::Integer,
                        sign::Integer,
                        nulldist::AbstractFloat,
                        field::Field3D)

        nlines = size(line1, 1)
        maxcount = 1000
        nseps = 0

        h0 = 5e-2 * stepsize
        tracedist = 3 * nulldist # trace points around null until tracedist away
        checkdist = tracedist + 2 * h0 # distance to check for separator after being iterated to tracedist away
    
        for inull in 1:size(rnulls, 1)
            if signs[inull] != sign # ignore nulls of the same sign (but not nulls with zero/undetermined sign - just in case)
                near = zeros(Int32, nlines)
                notnear = zeros(Int32, nlines)
        
                # first find all points that lie within nulldist and note their index
                for iline = 1:nlines
                    if (Common.dist(line1[iline], rnulls[inull]) < nulldist) | (Common.dist(line1[iline], rnullsalt[inull]) < nulldist)
                        near[iline] = iline # point lies within nulldist
                    end
                end
        
                if maximum(near) > 0 # if there are any points that lie within this distance
                    # find the longest consectutive number of points not near to the null (should be more not near another null)
                    nnc = 0
                    for iline in 1:nlines
                        if near[iline] == 0
                            nnc = nnc + 1
                            notnear[iline] = nnc
                        else
                            notnear[iline] = 0
                            nnc = 0
                        end
                    end
                    if nnc != 0
                        iline = 1
                        while near[iline] == 0
                            nnc = nnc + 1
                            notnear[iline] = nnc
                            iline = iline + 1
                        end
                    end
            
                    endgap = argmax(notnear) # location of last non-near point
                    gapsize = notnear[endgap] # number of points in biggest gap not near null
                    nextra = 10 # number of points to test either side of test points
                    nr = nlines - gapsize + 2*nextra # total number to test (nlines-gapsize near null)
            
                    if (gapsize == 0) | (nr >= nlines)
                        n1 = 1
                        n2 = nlines
                        nr = nlines
                    else
                        # select all points not in this longest chain to be tested for change in side of the fan
                        n1 = (endgap + 1 - nextra - 1) % nlines + 1
                        n2 = (nlines + endgap - gapsize + nextra - 1) % nlines + 1
                    end
            
                    if n1 <= n2
                        r = line1[n1:n2]
                        rmap = n1:n2
                    else
                        r = vcat(line1[n1:nlines], line1[1:n2])
                        rmap = vcat(n1:nlines, 1:n2)
                    end
            
                    signof = zeros(Int32, nr)
                    for index in 1:nr
                        # extrapolate points along fieldlines
                        count = 0
                        pt = r[index]
                        while (Common.dist(pt, rnulls[inull]) < tracedist) & (count < maxcount)
                            h = h0
                            pt = Trace.trace_line(pt, sign, h, field)
                            Common.edgecheck(pt, field)
                            count = count + 1
                        end
            
                        if signs[inull] * sign == -1
                            # check which side of the null the points end out on
                            if Common.dot(spines[inull], r[index] - rnulls[inull]) > 0
                                signof[index] = 1
                            else
                                signof[index] = -1
                            end
                
                            # if theres a change in sign, theres the separator
                            if index != 1
                                dist1 = Common.dist(rnulls[inull], r[index-1])
                                dist2 = Common.dist(rnulls[inull], r[index])
                                if ((signof[index-1] * signof[index] == -1) & (breaks[rmap[index-1]] != 1)
                                        & (dist1 < checkdist) & (dist2 < checkdist))
                                    println("Found a separator to null $(inull) on ring $(nring)")
                                    @debug "Found a separator on ring $nring at index $(rmap[index-1]) out of $nlines points from null number $nullnum to null number $inull"
                                    breaks[rmap[index-1]] = 1 # disassociate points so that new points don't get added between them as they diverge around the null
                                    nseps += 1
                                    # write the point's information to the separator file
                                    # if (dist1 > dist2) then
                                    #     write(40) nullnum, inull, nring, rmap(index)
                                    # else
                                    #     write(40) nullnum, inull, nring, rmap(index-1)
                                    # end
                                    if one_sep_per_ring
                                        break
                                    end
                                end
                            end
                        end
                        if count == maxcount # points appear to be stuck at the null
                            if rmap[index] != 1
                                breaks[rmap[index] - 1] = 2
                            else
                                breaks[nlines] = 2
                            end
                        end
                        # count1 = count
                    end
                end
            end
        end
        return nseps
    end

end