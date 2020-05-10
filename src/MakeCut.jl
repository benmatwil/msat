module MakeCut

    using StaticArrays
    using OffsetArrays
    using ..Params
    using ..Common
    using ..Read

    struct FieldGrid
        x::Vector{Float64}
        y::Vector{Float64}
        z::Vector{Float64}
    end

    function MakeCut(filename::AbstractString)
        nnulls, _, _ = Read.read_nulls(filename)
        
        # read in only the grid coordinates - no need to magnetic field values here
        fieldfile = open(filename, "r")
        nx, ny, nz = read(fieldfile, Tuple{Int32, Int32, Int32})
        seek(fieldfile, 3 * 4 + 3 * nx * ny * nz * 8)
        grid = FieldGrid(read!(fieldfile, Vector{Float64}(undef, nx)),
                         read!(fieldfile, Vector{Float64}(undef, ny)),
                         read!(fieldfile, Vector{Float64}(undef, nz)))
        close(fieldfile)

        # work out how often to print for progress update
        if nnulls < 10
            iprint = 1
        else
            iorder = floor(log10(nnulls))
            iprint = (nnulls ÷ (10^iorder) + 1) * 10^(iorder-1)
        end

        prefix = joinpath("default_output", Read.prefix(filename))

        info_file = open(prefix * "-ringinfo.dat", "r")
        nrings, nrings, nrings, hstep = read(info_file, Tuple{Int32, Int32, Int32, Int32}) 
        close(info_file)

        ds = max(maximum(diff(grid.x)), maximum(diff(grid.y)), maximum(diff(grid.z))) * 2 * hstep
        disttol = ds/10

        ###############################################################################################

        # for each spine, check whether we cross desired plane and add point to file
        println("Spines")
        spines_cut_file = open(prefix * "-spines-cut_" * file_pln * ".dat", "w+")
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
                    if check_crossing(line[iline], line[nextline])
                        if Common.dist(line[iline], line[nextline]) < 1
                            point = find_crossing(line[iline], line[nextline])
                            write(spines_cut_file, inull, point)
                            ipt += one(Int32)
                            push!(spine_list, point)
                        end
                    end
                end
            end
        end
        seekstart(spines_cut_file)
        write(spines_cut_file) ipt
        close(spines_cut_file)
        close(spines_file)
      
        # ********************************************************************************
      
        if do_all_ssf
            # for each separator, check whether we cross desired plane and add point to file
            println("Separators")
            separators_cut_file = open(prefix * "-separators-cut_" * file_pln * ".dat", "w+")
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
                        if check_crossing(line[iline], line[nextline])
                            if Common.dist(line[iline], line[nextline]) < 1
                                point = find_crossing(line[iline], line[nextline])
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
            ring_cut_file = open(prefix * "-rings-cut_"//file_pln//".dat", "w+")
            ringinfo_file = open(prefix * "-ringinfo.dat", "r")
            rings_file = open(prefix * "-rings.dat", "r")
            breaks_file = open(prefix * "-breaks.dat", "r")
            try
                global assocs_file = open(prefix * "-assocs.dat", "r")
            catch
                throw("Associations file does not exist.")
            end
            
            nlines = zero(Int32)
            write(ring_cut_file, nlines)
        
            nrings = read(ringinfo_file, Int32)
            nskip_file = read(ringinfo_file, Int32)

            num_rings = ceil(nrings/nskip_file)
            nperring = OffsetVector(Vector{Int32}(undef, num_rings), 0:num_rings-1)
            skip(ringinfo_file, 12) # just read into unneeded variable to reach correct position in file
        
            uptoassoc = 0
            uptorings = 0
            for inull in 1:nnulls
                
                read!(ringinfo_file, nperring) 
                nrings = count(nperring > 0) - 1 # index of nperring starts at 0
        
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
            
                    ia = uptoassoc + 1
                    ip = uptorings + 1

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
        
                for iring in nrings:-1:1
        
                    uptoring1 = sum(nperring[0:iring-2])
                    uptoring2 = uptoring1 + nperring[iring-1]
                    
                    # if have to use the slow method start reading in ring by ring
                    # otherwise use the arrays read in above using fast method
                    if too_much_memory
                        ia = uptoassoc + uptoring2*4 + 1
                        seek(assocs_file, ia)
                        seek(breaks_file, ia)
                        association = read!(assocs_file, Vector{Int32}(undef, nperring[iring]))
                        break_data = read!(breaks_file, Vector{Int32}(undef, nperring[iring]))
                        
                        if iring == nrings
                            seek(rings_file, uptorings + uptoring2 * 24 + 1)
                            line = read!(rings_file, Vector{Vector3D}(undef, nperring[iring]))
                        end
            
                        seek(rings_file, uptorings + uptoring1 * 24 + 1)
                        line2 = read!(rings_file, Vector{Vector3D}(undef, nperring[iring-1]))
                    else
                        association = associations[uptoring2+1:uptoring2+nperring[iring]]
                        break_data = breaks[uptoring2+1:uptoring2+nperring[iring]]
                        if iring == nrings
                            line = lines[uptoring2+1:uptoring2+nperring[iring]]
                        end
                        line2 = lines[uptoring1+1:uptoring2]
                    end
            
                    # find crossing of rings
                    # first by using associations from one ring to next
                    # second by going around rings in order of points
                    for iline in 1:nperring[iring]
                        if check_crossing(line[iline], line[association[iline]])
                            if dist(line[iline], line2[association[iline]]) < 1
                                # find point in crossing plane and add to list of points
                                point = find_crossing(line[iline], line2[association[iline]])
                                push!(ring_list, point)
                            end
                        end
                    end
                    for iline in 1:nperring(iring)-1
                        nextline = iline + 1
                        if check_crossing(line[iline], line[nextline])
                            if (break_data[iline] == 0) & (dist(line[iline], line[nextline]) < 1)
                                # find point in crossing plane and add to list of points
                                point = find_crossing(line[iline], line[nextline])
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
        
                    pordered, exit_status = sortpoints(points, spines)
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
            write(ring_cut_file) nlines
            close(ring_cut_file)
            close(rings_file)
            close(breaks_file)
            close(assocs_file)
        
            println("Finished null rings")
        
        end
      
        # ********************************************************************************
      
        if do_hcs
            # same as the rings above but for the hcs rings - see above for comments
            hcs_ring_file = open(prefix * "-hcs-cut_"//file_pln//".dat", "w+")
            hcs_ringinfo_file = open(prefix * "-hcs-ringinfo.dat", "r")
            hcs_rings_file = open(prefix * "-hcs-rings.dat", "r")
            try
                hcs_assocs_file = open(prefix * "-hcs-assocs.dat", "r")
            catch
                throw("HCS Associations file does not exist.")
            end
        
            nlines = 0
            write(80) nlines
            
            read(40) nrings
            read(40) nskip_file
            allocate(nperring(0:ceiling(real(nrings)/nskip_file-1)))
            read(40) nrings, nrings, nrings, ncomp
            
            uptoassoc = 0
            uptorings = 0
        
            # uptonull = 0
            for ihcs in 1:ncomp
                allocate(points(3, 0)) # fix
                ring_list = Vector{Vector3D}(undef, 0)
                read(40) nperring
        
                nrings = count(nperring .> 0) - 1
                num_nums = sum(nperring)
                allocate(lines(3, num_nums), associations(num_nums))
                ia = uptoassoc + 1
                ip = uptorings + 1
        
                if num_nums > large_read
                    leftover_num = num_nums
                    read_index = 1
                    inquire(50, pos=ip)
                    inquire(55, pos=ia)
                    while (leftover_num > 0)
                        if (leftover_num > large_read) then
                            read_num = large_read
                        else
                            read_num = leftover_num
                        end
                        read(50) lines(:, read_index:read_index+read_num)
                        read(55) associations(read_index:read_index+read_num)
                        read_index = read_index + read_num
                        leftover_num = leftover_num - read_num
                    end
                else
                    read(50, pos=ip) lines
                    read(55, pos=ia) associations
                end
        
                for iring in nrings:-1:1
        
                    uptoring1 = sum(int(nperring(0:iring-2), int64))
                    uptoring2 = uptoring1 + nperring(iring-1)
            
                    # ia = uptoassoc + uptoring2*4_int64 + 1
                    # ip = uptorings + uptoring2*24_int64 + 1
                    
                    # allocate(association(nperring(iring)), line(3, nperring(iring)))
                    # read(55, pos=ia) association
                    # read(50, pos=ip) line
            
                    # ip = uptorings + uptoring1*24_int64 + 1
            
                    # allocate(line2(3,nperring(iring-1)))
                    # read(50, pos=ip) line2
            
                    association = associations(uptoring2+1:uptoring2+nperring(iring))
                    if (iring == nrings) line = lines(:, uptoring2+1:uptoring2+nperring(iring)) end
                    line2 = lines(:, uptoring1+1:uptoring2)
            
                    for iline in 1:nperring(iring)
                        if (check_crossing(line(:, iline), line2(:, association(iline)))) then
                            if (dist(line(:, iline), line2(:, association(iline))) < 1) then
                                # find point in between at r0 and add to line
                                point = find_crossing(line(:, iline), line2(:, association(iline)))
                                push!(ring_list, point)
                            end
                        end
                    end
                    line = line2
                    deallocate(association)
                end
                deallocate(lines, associations)
                if (allocated(line)) deallocate(line) end
                uptoring1 = sum(int(nperring, int64))
                uptoassoc = uptoassoc + uptoring1*4_int64
                uptorings = uptorings + uptoring1*24_int64
        
                points = ring_list%to_array()
                # call ring_list%destroy()
        
                println("sorting hcs points")
                println(size(points, 1))
        
                while (size(points, 2) > 0)
        
                    pordered, exit_status = sortpoints(points)
                    if (exit_status) then
                        exit
                    else
                        write(80) int(0, int32), size(pordered, 2), pordered
                        nlines = nlines + 1
                        deallocate(pordered)
                    end
        
                end
        
                deallocate(points)
            end
        
            write(80, pos=1) nlines
            close(80)
            close(40)
            close(50)
            close(55)
            deallocate(nperring)
        
            # ********************************************************************************
        
            # for each hcs separator, check whether we cross the plane and add point to file
            # need to do this after rings (can't remember why)
            open(unit=80, file=trim(fileout)//"-hcs-separators-cut_"//file_pln//".dat", "w")
            open(unit=20, file=trim(fileout)//"-hcs-connectivity.dat", "r")
            open(unit=30, file=trim(fileout)//"-hcs-separators.dat", "r")
            
            ipt = 0
            write(80) ipt
            read(20) nseps
            for isep in 1:nseps
                read(20) inull, flag
                read(30) npoints
                allocate(line(3, npoints))
                read(30) line
                for iline in 1:npoints-1
                    nextline = iline+1
                    if (check_crossing(line(:, iline), line(:, nextline))) then
                        if (dist(line(:, iline), line(:, nextline)) < 1) then
                        point = find_crossing(line(:, iline), line(:, nextline))
                        write(80) flag, point
                        ipt = ipt + 1
                        end
                    end
                end
                read(20) flag, flag # don't care about these - just reposition file marker
                deallocate(line)
            end
            
            write(80, pos=1) ipt
            close(80)
            close(30)
            close(20)
        end
      
        # call system_clock(tstop, count_rate)
        # print*, 'Time taken:', dble(tstop - tstart)/dble(count_rate), "(", dble(tstop - tstart)/dble(count_rate)/60, "minutes)"
      
        # call spine_list%destroy()

    end

    function sortpoints(points::Vector{Vector3D}, spines::Vector{Vector3D}, disttol::Float64)

        exit_status = false

        # search for closest points in first 5000 points if large number of points
        if size(points, 1) > 5000
            npoints = 5000
        else
            npoints = size(points, 2)
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
                            deleteat!(pt_list, imindist)
                        elseif (Common.dot(diff, Common.normalise(points[imindist] - pt1)) > 0.95) & (mindist < 5*disttol) # line is almost straight
                            pushfirst!(pt_list, points[imindist])
                            deleteat!(pt_list, imindist)
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

    function check_crossing(r1::Vector3D, r2::Vector3D)
        # checks if line between two points crosses plane

        return plane(r1) * plane(r2) <= 0

    end

    function find_crossing(r1::Vector3D, r2::Vector3D, plane_norm::Vector3D, plane_const::Float64)
        # finds crossing of plane using linear interpolation of two points

        s = (plane_const - Common.dot(r1, plane_norm)) / (Common.dot(r2 - r1, plane_norm))
        return r1 + s*(r2 - r1)

    end


end