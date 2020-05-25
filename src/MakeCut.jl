module MakeCut

    using StaticArrays
    using OffsetArrays
    using Formatting

    using ..Params
    using ..Common
    using ..Read

    function get_ints(n::Integer)
        return Vector{Int32}(undef, n)
    end

    function MC(bgrid::AbstractField3D, plane_a::Real, plane_b::Real, plane_c::Real, plane_const::Float64; do_all_ssf=true, do_hcs=true)
        MC(bgrid, Float64[plane_a, plane_b, plane_c], plane_const; do_all_ssf=do_all_ssf, do_hcs=do_hcs)
    end

    function MC(bgrid::AbstractField3D, plane_norm::Vector{Float64}, plane_const::Float64; do_all_ssf=true, do_hcs=true)
        MC(bgrid, Vector3D(plane_norm), plane_const; do_all_ssf=do_all_ssf, do_hcs=do_hcs)
    end

    function MC(bgrid::AbstractField3D, plane_norm::Vector3D, plane_const::Float64; do_all_ssf=true, do_hcs=true)
        nnulls, _, _ = Read.read_nulls(bgrid.filename)
        
        # work out how often to print for progress update
        if nnulls < 10
            iprint = 1
        else
            iorder = floor(log10(nnulls))
            iprint = (nnulls ÷ (10^iorder) + 1) * 10^(iorder-1)
        end

        prefix = joinpath(default_output, Read.prefix(bgrid.filename))
        plane_nums = [plane_norm..., plane_const]
        file_plane = format(join([i >= 0 ? raw"{:08.4f}" : raw"{:09.4f}" for i in plane_nums], "_"), plane_nums...)

        info_file = open(prefix * "-ringinfo.dat", "r")
        skip(info_file, 12)
        hstep = read(info_file, Float64)
        close(info_file)

        ds = max(maximum(diff(bgrid.x)), maximum(diff(bgrid.y)), maximum(diff(bgrid.z))) * 2 * hstep
        disttol = ds/10

        @show ds, disttol, hstep

        large_read = 500000000

        ###############################################################################################

        # for each spine, check whether we cross desired plane and add point to file
        println("Spines")
        spines_cut_file = open(prefix * "-spines-cut_" * file_plane * ".dat", "w+")
        spines_file = open(prefix * "-spines.dat", "r")
      
        spine_list = Vector{Vector3D}(undef, 0)
      
        ipt = zero(Int32)
        write(spines_cut_file, ipt)
        
        for inull::Int32 in 1:nnulls
            # read in spines for each direction
            for dir in 1:2
                nspine = read(spines_file, Int32)
                line = read!(spines_file, Vector{Vector3D}(undef, nspine))
                for iline in 1:nspine-1
                    nextline = iline + 1
                    if check_crossing(line[iline], line[nextline], plane_norm, plane_const)
                        if Common.dist(line[iline], line[nextline]) < 1
                            point = find_crossing(line[iline], line[nextline], plane_norm, plane_const)
                            write(spines_cut_file, inull, point)
                            ipt += one(Int32)
                            push!(spine_list, point)
                        end
                    end
                end
            end
        end
        seekstart(spines_cut_file)
        write(spines_cut_file, ipt)
        close(spines_cut_file)
        close(spines_file)

        spines = spine_list
      
        # ********************************************************************************
      
        if do_all_ssf
            # for each separator, check whether we cross desired plane and add point to file
            println("Separators")
            separators_cut_file = open(prefix * "-separators-cut_" * file_plane * ".dat", "w+")
            connectivity_file = open(prefix * "-connectivity.dat", "r")
            separator_file = open(prefix * "-separators.dat", "r")
        
            ipt = zero(Int32)
            write(separators_cut_file, ipt)
            for inull in 1:nnulls
                nseps = read(connectivity_file, Int32)
                for isep in 1:nseps
                    skip(connectivity_file, 4)
                    flag = read(connectivity_file, Int32)
                    npoints = read(separator_file, Int32)
                    line = read!(separator_file, Vector{Vector3D}(undef, npoints))
                    for iline in 1:npoints-1
                        nextline = iline+1
                        if check_crossing(line[iline], line[nextline], plane_norm, plane_const)
                            if Common.dist(line[iline], line[nextline]) < 1
                                point = find_crossing(line[iline], line[nextline], plane_norm, plane_const)
                                write(separators_cut_file, Int32(inull), flag, point)
                                ipt += one(ipt)
                            end
                        end
                    end
                    skip(connectivity_file, 8)
                end
            end
            
            seekstart(separators_cut_file)
            write(separators_cut_file, ipt)
            close(separator_file)
            close(separators_cut_file)
            close(connectivity_file)
        
            # ********************************************************************************
        
            # for each adjacent ring, check whether we cross desired plane and add point to file
            println("Rings")
            ring_cut_file = open(prefix * "-rings-cut_" * file_plane * ".dat", "w+")
            ringinfo_file = open(prefix * "-ringinfo.dat", "r")
            rings_file = open(prefix * "-rings.dat", "r")
            breaks_file = open(prefix * "-breaks.dat", "r")
            try
                global assocs_file
                assocs_file = open(prefix * "-assocs.dat", "r")
            catch
                throw("Associations file does not exist.")
            end
            
            nlines = zero(Int32)
            write(ring_cut_file, nlines)
        
            nrings = read(ringinfo_file, Int32)
            nskip_file = read(ringinfo_file, Int32)
            skip(ringinfo_file, 12)

            num_rings = ceil(Int32, nrings/nskip_file)
            nperring = OffsetVector(Vector{Int32}(undef, num_rings), 0:num_rings-1)
        
            uptoassoc = 0
            uptorings = 0
            for inull in 1:nnulls
                
                read!(ringinfo_file, nperring) 
                nrings = count(nperring .> 0) - 1 # index of nperring starts at 0
        
                # two options to read in points: one fast and one slow
                # depends on if all points too big for memory/max fortran array size
                # if number of points isn't too big... uses fast method
                # and reads all points in now
                num_nums = sum(nperring)
                if num_nums > 2500000000
                    too_much_memory = true
                else
                    too_much_memory = false
                end
        
                if ! too_much_memory
                    # 500000000 seems to be maximum number of int32 integers fortran can read at a time
                    # don't have concrete evidence for this number... may need changing but seems to work
                    # allocate(lines(3, num_nums), associations(num_nums), breaks(num_nums))
                    lines = Vector{Vector3D}(undef, num_nums)
                    associations = Vector{Int32}(undef, num_nums)
                    breaks = Vector{Int32}(undef, num_nums)
            
                    ia = uptoassoc
                    ip = uptorings

                    seek(rings_file, ip)
                    seek(assocs_file, ia)
                    seek(breaks_file, ia)
                    
                    # if bigger than a certain read size then read in chunks
                    # otherwise read all of it
                    if num_nums > large_read
                        leftover_num = num_nums
                        read_index = 1
                        while leftover_num > 0
                            if leftover_num > large_read
                                read_num = large_read
                            else
                                read_num = leftover_num
                            end
                            read!(rings_file, lines[read_index:read_index+read_num-1])
                            read!(assocs_file, associations[read_index:read_index+read_num-1])
                            read!(breaks_file, breaks[read_index:read_index+read_num-1])
                            read_index += read_num
                            leftover_num -= read_num
                        end
                    else
                        read!(rings_file, lines)
                        read!(assocs_file, associations)
                        read!(breaks_file, breaks)
                    end
                end
                
                ring_list = Vector{Vector3D}(undef, 0)

                uptoring1 = sum(nperring[0:nrings-2])
                uptoring2 = uptoring1 + nperring[nrings-1]
                if too_much_memory
                    seek(rings_file, uptorings + uptoring2 * 24)
                    line = read!(rings_file, Vector{Vector3D}(undef, nperring[nrings]))
                else
                    line = lines[uptoring2+1:uptoring2+nperring[nrings]]
                end
        
                for iring in nrings:-1:1
                    # global line
        
                    uptoring1 = sum(nperring[0:iring-2])
                    uptoring2 = uptoring1 + nperring[iring-1]
                    
                    # if have to use the slow method start reading in ring by ring
                    # otherwise use the arrays read in above using fast method
                    if too_much_memory
                        ia = uptoassoc + uptoring2*4
                        seek(assocs_file, ia)
                        association = read!(assocs_file, Vector{Int32}(undef, nperring[iring]))
                        seek(breaks_file, ia)
                        break_data = read!(breaks_file, Vector{Int32}(undef, nperring[iring]))
                        seek(rings_file, uptorings + uptoring1 * 24 )
                        line2 = read!(rings_file, Vector{Vector3D}(undef, nperring[iring-1]))
                    else
                        association = associations[uptoring2+1:uptoring2+nperring[iring]]
                        break_data = breaks[uptoring2+1:uptoring2+nperring[iring]]
                        line2 = lines[uptoring1+1:uptoring2]
                    end
            
                    # find crossing of rings
                    # first by using associations from one ring to next
                    # second by going around rings in order of points
                    for iline in 1:nperring[iring]
                        if check_crossing(line[iline], line2[association[iline]], plane_norm, plane_const)
                            if Common.dist(line[iline], line2[association[iline]]) < 1
                                # find point in crossing plane and add to list of points
                                point = find_crossing(line[iline], line2[association[iline]], plane_norm, plane_const)
                                push!(ring_list, point)
                            end
                        end
                    end
                    for iline in 1:nperring[iring]-1
                        nextline = iline + 1
                        if check_crossing(line[iline], line[nextline], plane_norm, plane_const)
                            if (break_data[iline] == 0) & (Common.dist(line[iline], line[nextline]) < 1)
                                # find point in crossing plane and add to list of points
                                point = find_crossing(line[iline], line[nextline], plane_norm, plane_const)
                                push!(ring_list, point)
                            end
                        end
                    end
                    line = line2
                end

                points = ring_list
        
                # find new read file points for slow method
                uptoring1 = sum(nperring)
                uptoassoc += uptoring1 * 4
                uptorings += uptoring1 * 24
        
                # sort points until list is exhausted
                while size(points, 1) > 0
        
                    pordered, exit_status = sortpoints(points, spines, disttol)
                    if exit_status
                        break
                    else
                        write(ring_cut_file, Int32(inull), Int32(size(pordered, 1)), pordered)
                        nlines += one(Int32)
                    end
        
                end
                
                if inull % iprint == 0 println("Done $inull of $nnulls") end
            end
            seekstart(ring_cut_file)
            write(ring_cut_file, nlines)
            close(ring_cut_file)
            close(rings_file)
            close(breaks_file)
            close(assocs_file)
        
            println("Finished null rings")
        
        end
      
        # ********************************************************************************
      
        if do_hcs
            # same as the rings above but for the hcs rings - see above for comments
            hcs_ring_cut_file = open(prefix * "-hcs-cut_" * file_plane * ".dat", "w+")
            hcs_ringinfo_file = open(prefix * "-hcs-ringinfo.dat", "r")
            hcs_rings_file = open(prefix * "-hcs-rings.dat", "r")
            try
                global hcs_assocs_file
                hcs_assocs_file = open(prefix * "-hcs-assocs.dat", "r")
            catch
                throw("HCS Associations file does not exist.")
            end
        
            nlines = zero(Int32)
            write(hcs_ring_cut_file, nlines)
            
            nrings = read(hcs_ringinfo_file, Int32)
            nskip_file = read(hcs_ringinfo_file, Int32)
            num_rings = ceil(Int32, nrings/nskip_file)
            nperring = OffsetArray(Vector{Int32}(undef, num_rings), 0:num_rings-1)
            skip(hcs_ringinfo_file, 12)
            ncomp = read(hcs_ringinfo_file, Int32)
            
            uptoassoc = 0
            uptorings = 0
        
            for ihcs in 1:ncomp
                ring_list = Vector{Vector3D}(undef, 0)
                read!(hcs_ringinfo_file, nperring)
        
                nrings = count(nperring .> 0) - 1
                num_nums = sum(nperring)
                ia = uptoassoc
                ip = uptorings

                seek(hcs_rings_file, ip)
                seek(hcs_assocs_file, ia)

                lines = Vector{Vector3D}(undef, num_nums)
                associations = Vector{Int32}(undef, num_nums)
        
                if num_nums > large_read
                    leftover_num = num_nums
                    read_index = 1
                    while (leftover_num > 0)
                        if leftover_num > large_read
                            read_num = large_read
                        else
                            read_num = leftover_num
                        end
                        read!(hcs_rings_file, lines[read_index:read_index+read_num])
                        read!(hcs_assocs_file, associations[read_index:read_index+read_num])
                        read_index = read_index + read_num
                        leftover_num = leftover_num - read_num
                    end
                else
                    read!(hcs_rings_file, lines)
                    read!(hcs_assocs_file, associations)
                end

                uptoring1 = sum(nperring[0:nrings-2])
                uptoring2 = uptoring1 + nperring[nrings-1]
                line = lines[uptoring2+1:uptoring2+nperring[nrings]]
        
                for iring in nrings:-1:1
        
                    uptoring1 = sum(nperring[0:iring-2])
                    uptoring2 = uptoring1 + nperring[iring-1]
            
                    association = associations[uptoring2+1:uptoring2+nperring[iring]]
                    line2 = lines[uptoring1+1:uptoring2]
            
                    for iline in 1:nperring[iring]
                        if check_crossing(line[iline], line2[association[iline]], plane_norm, plane_const)
                            if Common.dist(line[iline], line2[association[iline]]) < 1
                                # find point in between at r0 and add to line
                                point = find_crossing(line[iline], line2[association[iline]], plane_norm, plane_const)
                                push!(ring_list, point)
                            end
                        end
                    end
                    line = line2
                end
                uptoring1 = sum(nperring)
                uptoassoc = uptoassoc + uptoring1 * 4
                uptorings = uptorings + uptoring1 * 24
        
                points = ring_list
        
                println("sorting hcs points")
                println(size(points, 1))
        
                while (size(points, 1) > 0)
        
                    pordered, exit_status = sortpoints(points, spines, disttol)
                    if exit_status
                        exit
                    else
                        write(hcs_ring_cut_file, Int32(0), Int32(size(pordered, 1)), pordered)
                        nlines += one(Int32)
                    end
        
                end
        
            end
        
            seekstart(hcs_ring_cut_file)
            write(hcs_ring_cut_file, nlines)
            close(hcs_ring_cut_file)
            close(hcs_ringinfo_file)
            close(hcs_rings_file)
            close(hcs_assocs_file)
        
            # ********************************************************************************
        
            # for each hcs separator, check whether we cross the plane and add point to file
            # need to do this after rings (can't remember why)
            hcs_separators_cut_file = open(prefix * "-hcs-separators-cut_" * file_plane * ".dat", "w+")
            hcs_connectivity_file = open(prefix * "-hcs-connectivity.dat", "r")
            hcs_separators_file = open(prefix * "-hcs-separators.dat", "r")
            
            ipt = zero(Int32)
            write(hcs_separators_cut_file, ipt)
            nseps = read(hcs_connectivity_file, Int32)
            for isep in 1:nseps
                inull, flag = read(hcs_connectivity_file, Tuple{Int32, Int32}) 
                npoints = read(hcs_separators_file, Int32) 
                line = read!(hcs_separators_file, Vector{Vector3D}(undef, npoints))
                for iline in 1:npoints-1
                    nextline = iline+1
                    if check_crossing(line[iline], line[nextline], plane_norm, plane_const)
                        if Common.dist(line[iline], line[nextline]) < 1
                            point = find_crossing(line[iline], line[nextline], plane_norm, plane_const)
                            write(hcs_separators_cut_file, flag, point)
                            ipt += one(Int32)
                        end
                    end
                end
                skip(hcs_connectivity_file, 8)
            end
            seekstart(hcs_separators_cut_file)
            write(hcs_separators_cut_file, ipt)
            close(hcs_separators_cut_file)
            close(hcs_separators_file)
            close(hcs_connectivity_file)
        end
      
        # call system_clock(tstop, count_rate)
        # print*, 'Time taken:', dble(tstop - tstart)/dble(count_rate), "(", dble(tstop - tstart)/dble(count_rate)/60, "minutes)"
      
    end

    function sortpoints(points::Vector{Vector3D}, spines::Vector{Vector3D}, disttol::Float64)

        exit_status = false

        # search for closest points in first 5000 points if large number of points
        if size(points, 1) > 5000
            npoints = 5000
        else
            npoints = size(points, 1)
        end

        # calculate point to point distances
        dists = fill(1000.0, npoints, npoints)
        for ip in 1:npoints
            dists[ip, ip+1:npoints] .= sqrt.(sum.(broadcast.(^, Ref(points[ip]) .- points[ip+1:npoints], 2)))
        end
        imindists = argmin(dists)
        
        if dists[imindists] > disttol*10
            # exit if closest points far away
            exit_status = true
        else
            # add two closest points to list to start the line
            pt_list = [points[imindists[1]], points[imindists[2]]]
            deleteat!(points, imindists[1])
            if imindists[1] < imindists[2]
                deleteat!(points, imindists[2]-1)
            else
                deleteat!(points, imindists[2])
            end

            # allocate array to store distances between current points and all spines
            nspines = size(spines, 1)
            spinedists = Vector{Float64}(undef, nspines)

            for dir in 1:2
                # add points in the each direction
                while size(points, 1) > 0
                    # set pt1/2 to be first two points in list
                    pt1 = pt_list[1]
                    pt2 = pt_list[2]

                    # vector connecting pt1/2
                    diff = Common.normalise(pt1 - pt2)

                    # required distances to first point
                    # spinedists = (spines(1, :) - pt1(1))**2 + (spines(2, :) - pt1(2))**2 + (spines(3, :) - pt1(3))**2
                    ptdists = sum.(broadcast.(^, points .- Ref(pt1), 2))
                    spinedists = sum.(broadcast.(^, spines .- Ref(pt1), 2))

                    # shortest distance from point to all other points and location in array
                    imindist = argmin(ptdists)
                    mindist = ptdists[imindist]
                    
                    # shortest distance from point to all spines
                    if nspines == 0
                        mindistspine = 1e200
                    else
                        mindistspine = minimum(spinedists)
                    end

                    # if new point is good, add to line
                    if mindist < mindistspine
                        if Common.dot(diff, Common.normalise(points[imindist] - pt1)) < -0.98 # sharp corner in line
                            break
                        elseif mindist < disttol
                            pushfirst!(pt_list, points[imindist])
                            deleteat!(points, imindist)
                        elseif (Common.dot(diff, Common.normalise(points[imindist] - pt1)) > 0.95) & (mindist < 5*disttol) # line is almost straight
                            pushfirst!(pt_list, points[imindist])
                            deleteat!(points, imindist)
                        else
                            break
                        end
                    else
                        break
                    end
                end

                # reverse line so can add to line using loop
                reverse!(pt_list)

            end

        end

        return pt_list, exit_status

    end

    function plane(r::Vector3D, plane_norm::Vector3D, plane_const::Float64)
        # function to calculate required plane

        return Common.dot(plane_norm, r) - plane_const

    end

    function check_crossing(r1::Vector3D, r2::Vector3D, plane_norm::Vector3D, plane_const::Float64)
        # checks if line between two points crosses plane

        return plane(r1, plane_norm, plane_const) * plane(r2, plane_norm, plane_const) <= 0

    end

    function find_crossing(r1::Vector3D, r2::Vector3D, plane_norm::Vector3D, plane_const::Float64)
        # finds crossing of plane using linear interpolation of two points

        s = (plane_const - Common.dot(r1, plane_norm)) / (Common.dot(r2 - r1, plane_norm))
        return r1 + s*(r2 - r1)

    end


end