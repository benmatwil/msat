module NullFinder

    using StaticArrays
    using ..Params
    using ..Common
    using ..Read

    function NF(bgrid::AbstractField3D; boundary_nulls::Bool=false, gridpoint_nulls::Bool=true)

        nx, ny, nz, ndims = size(bgrid.field)

        magb = sqrt.(sum(bgrid.field .^ 2, dims=4))

        candidates = falses((nx-1, ny-1, nz-1))
        for iz in 1:nz-1, iy in 1:ny-1, ix in 1:nx-1
            cube = bgrid.field[ix:ix+1, iy:iy+1, iz:iz+1, :]
            if allsame(cube[:, :, :, 1]) continue end
            if allsame(cube[:, :, :, 2]) continue end
            if allsame(cube[:, :, :, 3]) continue end
            if sum(magb[ix:ix+1, iy:iy+1, iz:iz+1]) < 8*zero_value continue end
            candidates[ix, iy, iz] = true
        end

        println("Number of candidate cells: $(count(candidates))")

        subcube = Array{Float64, 4}(undef, 11, 11, 11, 3)

        nullpts = Vector{Vector3D}(undef, 0)

        for iz in 1:nz-1, iy in 1:ny-1, ix in 1:nx-1
            if candidates[ix, iy, iz]
                rnull = Vector3D([ix, iy, iz])
                for isig in 1:sig_figs
                    ds = 10.0^(-isig)
                    for izsub in 1:11, iysub in 1:11, ixsub in 1:11

                        rtri = rnull + [ixsub-1, iysub-1, izsub-1]*ds
                        subcube[ixsub, iysub, izsub, :] = Common.trilinear_nf(rtri, bgrid)

                    end
                    for izsub in 1:10, iysub in 1:10, ixsub in 1:10
                        global itest

                        cube = subcube[ixsub:ixsub+1, iysub:iysub+1, izsub:izsub+1, :]
                        for i in 1:ndims
                            if allsame(cube[:, :, :, i]) continue end
                        end
                        itest = bilin_test(cube)

                        if itest == 1
                            rnull = rnull + [ixsub-1, iysub-1, izsub-1]*ds
                            @goto end_outer
                        end
                    end
                    @label end_outer

                    if itest == 0
                        if isig > 1
                            println("Candidate null in cell: [$(ix), $(iy), $(iz)]")
                            println("Don't believe this null. ERROR :(")
                            println("Point to $(isig) sigfigs:, $(Float32.(rnull))")
                            println("|B| at point: $(sqrt(sum(Common.trilinear_nf(rnull, bgrid).^2)))")
                            println("")
                        end
                        @goto end_sigloop
                    end

                end
                @label end_sigloop

                if itest == 1
                    ds = 10.0^(-sig_figs)
                    cube = Array{Float64, 4}(undef, 2, 2, 2, 3)
                    for izsub in 1:2, iysub in 1:2, ixsub in 1:2
                        # find all field values of subgrid
                        rtri = rnull + [ixsub-1, iysub-1, izsub-1]*ds
                        cube[ixsub, iysub, izsub, :] = Common.trilinear_nf(rtri, bgrid)
                    end
        
                    mincube = Tuple(argmin(dropdims(sum(cube .^ 2, dims=4), dims=4)))
        
                    rnull = rnull .+ (mincube .- 1) .* ds
        
                    bound_dist = rspherefact * 10.0^(-sig_figs)
                    
                    if boundary_nulls
                        push!(nullpts, rnull)
                    else
                        if ((rnull[1] > 1 + bound_dist) & (rnull[1] < nx - bound_dist) &
                            (rnull[2] > 1 + bound_dist) & (rnull[2] < ny - bound_dist) &
                            (rnull[3] > 1 + bound_dist) & (rnull[3] < nz - bound_dist))
                            push!(nullpts, rnull)
                        else
                            @debug "Null on the boundary: $(Int.([ix, iy, iz])), Removing..."
                        end
                    end
                end

            end
        end

        println("Found $(size(nullpts, 1)) null points")

        if gridpoint_nulls
            indices = findall(magb .< zero_value)
            append!(nullpts, Common.Vector3D.(Tuple.(indices)))
            println("Found $(size(indices, 1)) nulls at vertices")
        end

        println("Checking for duplicate nulls:")
        
        if size(nullpts, 1) > 1
            exitcondition = false

            @label main_while
            while (!exitcondition)
                for j in 1:size(nullpts, 1), i in j+1:size(nullpts, 1)
                    if Common.dist(nullpts[i], nullpts[j]) < 1e-4
                        println("removing duplcate at index $j, $(Common.dist(nullpts[i], nullpts[j]))")
                        deleteat!(nullpts, j)
                        @goto main_while
                    end
                end
                exitcondition = true
            end

            println("")
            println("All null pairs with less than 1 gridcell spacing")
            counter = 0
            for j in 1:size(nullpts, 1), i in j+1:size(nullpts, 1)
                distance = Common.dist(nullpts[i], nullpts[j])
                if distance < 1
                    println("$i, $j, $distance")
                    counter = counter + 1
                end
            end
            println("Number of close-by nulls = $(counter)")
        end

        println("Final number of nulls: $(size(nullpts, 1))")
        println("At locations:")
        for pt in nullpts
            println(pt)
        end

        for pt in Common.gtr.(nullpts, Ref(bgrid))
            println(pt)
        end

        outfname = joinpath(default_output, Read.prefix(bgrid.filename) * "-nullpos.dat")        
        println("Writing null positions to $(outfname)")
        nullfile = open(outfname, "w")
            write(nullfile, Int32(size(nullpts, 1)))
            write(nullfile, nullpts)
            write(nullfile, Common.gtr.(nullpts, Ref(bgrid)))
        close(nullfile)

        # print final run time
        # call system_clock(tstop, count_rate)
        # print*, 'Time taken:', dble(tstop - tstart)/dble(count_rate), "(", dble(tstop - tstart)/dble(count_rate)/60, "minutes)"

    end
        
    function allsame(cube::Array{Float64, 3})
        return all(cube .> zero_value) | all(cube .< -zero_value)
    end

    function bilin_test(cube)
        
        test = zeros(Int32, 6, 6)
        edge = 0
        null_test = false

        for (num, index) in enumerate(((:, :, 1),
                                       (1, :, :),
                                       (2, :, :),
                                       (:, 1, :),
                                       (:, 2, :),
                                       (:, :, 2))) # not constant type
            facex = cube[index..., 1]
            facey = cube[index..., 2]
            facez = cube[index..., 3]

            cross, sign = face_solve(facex, facey, facez)
            test[num, 1] = cross
            test[num, 2] = sign
            cross, sign = face_solve(facey, facez, facex)
            test[num, 3] = cross
            test[num, 4] = sign
            cross, sign = face_solve(facex, facez, facey)
            test[num, 5] = cross
            test[num, 6] = sign

            edge += edge_check(facex, facey, facez)
        end

        test2 = sum(test, dims=1)
        test3 = zeros(Int32, 3)
        for i in 1:3
            if test2[2*i - 1] >= 0
                if abs(test2[2*i]) < test2[2*i - 1]
                    test3[i] = 1
                end
            end
        end

        if sum(test3) == 3
            null_test = true
        end

        if sum(test[:, 2]) > 10 | sum(test[:, 4]) > 10 | sum(test[:, 6]) > 10
            null_test = true
        end

        if edge > 0
            null_test = true
        end

        return null_test

    end

    function check(x::Float64, y::Float64)
        # checks of x and y are between 0 and 1
        return (0.0 <= x <= 1.0) & (0.0 <= y <= 1.0)
    end

    function face_solve(facex, facey, facez)

        sign, cross = 0, 0
        nsol = 0

        same_sign = ((count(facex .> zero_value) == 4 ) | (count(facex .< -zero_value) == 4) |
                     (count(facey .> zero_value) == 4 ) | (count(facey .< -zero_value) == 4))

        not_zero = sum(abs.(facex) + abs.(facey) + abs.(facez)) > 12*zero_value

        if !same_sign & not_zero
            a1 = facex[1, 1]
            b1 = facex[2, 1] - facex[1, 1]
            c1 = facex[1, 2] - facex[1, 1]
            d1 = facex[2, 2] - facex[2, 1] - facex[1, 2] + facex[1, 1]

            a2 = facey[1, 1]
            b2 = facey[2, 1] - facey[1, 1]
            c2 = facey[1, 2] - facey[1, 1]
            d2 = facey[2, 2] - facey[2, 1] - facey[1, 2] + facey[1, 1]

            a = b1 * d2 - b2 * d1
            b = (a1 * d2 - a2 * d1) + (b1 * c2 - c1 * b2)
            c = a1 * c2 - a2 * c1

            det = b^2 - 4*a*c

            x, y = fill(-1.0, 2, 2), fill(-1.0, 2, 2)

            if det >= 0 # there is a solution
                if abs(a) < zero_value # have to solve linear
                    if b != 0.0
                        x[1, :] .= -c/b # solution exists
                        nsol = 1
                    end
                else # have to solve quadratic
                    if det == 0.0 # one solution
                        x[1, :] .= -b/(2*a)
                        nsol = 1
                    else # two solutions
                        x[1, :] .= (-b + sqrt(det))/(2*a)
                        x[2, :] .= (-b - sqrt(det))/(2*a)
                        nsol = 2
                    end
                end
            end

            if nsol > 0
                y[1, 1] = -(a1 + b1 * x[1, 1])/(c1 + d1 * x[1, 1])
                y[1, 2] = -(a2 + b2 * x[1, 2])/(c2 + d2 * x[1, 2])
                if nsol == 2
                    y[2, 1] = -(a1 + b1 * x[2, 1])/(c1 + d1 * x[2, 1])
                    y[2, 2] = -(a2 + b2 * x[2, 2])/(c2 + d2 * x[2, 2])
                end
            else
                # if the x equation isn't solvable then try the y equation
                # get quadratic coefficients (ay^2 + by + c = 0)
                a = c1 * d2 - c2 * d1
                b = (a1 * d2 - a2 * d1) - (b1 * c2 - c1 * b2)
                c = a1 * b2 - a2 * b1

                det = b^2 - 4*a*c

                if det >= 0 # there is a solution
                    if abs(a) < zero_value # have to solve linear
                        if b != 0.0
                            y[1, :] .= -c/b # solution exists
                            nsol = 1
                        end
                    else # have to solve quadratic
                        if det == 0.0 # one solution
                            y[1, :] .= -b/(2*a)
                            nsol = 1
                        else # two solutions
                            y[1, :] .= (-b + sqrt(det))/(2*a)
                            y[2, :] .= (-b - sqrt(det))/(2*a)
                            nsol = 2
                        end
                    end
                    if nsol > 0
                        x[1, 1] = -(a1 + c1 * y[1, 1])/(b1 + d1 * y[1, 1])
                        x[1, 2] = -(a2 + c2 * y[1, 2])/(b2 + d2 * y[1, 2])
                        if nsol == 2
                            x[2, 1] = -(a1 + c1 * y[2, 1])/(b1 + d1 * y[2, 1])
                            x[2, 2] = -(a2 + c2 * y[2, 2])/(b2 + d2 * y[2, 2])
                        end
                    end
                end
            end

            for ix in 1:2, iy in 1:2
                if check(x[ix, iy], y[ix, iy]) # if x and y are on face (between 0 and 1)
                    cross += 1
            
                    zcomp = bilinear_cell(x[ix, iy], y[ix, iy], facez) # 'z' component of field at crossing point
                    if zcomp > zero_value
                        sign += 1
                    elseif zcomp < -zero_value
                        sign += -1
                    else
                        sign += 100
                    end
                end
            end

        end

        return cross, sign

    end

    function edge_check(facex, facey, facez)
        # check for null along edges of a cell face

        edge_check = 0

        for index in ((:, 1), (:, 2), (1, :), (2, :)) # not a constant type

            linex = facex[index...]
            liney = facey[index...]
            linez = facez[index...]

            edge_check += line_check(linex, liney, linez) # check for null along this edge

        end

        if edge_check > 0
            edge_check = 1 # if nulls are found, set the counter to 1
        end

        return edge_check

    end

    function line_check(linex, liney, linez)
        # checks if null lies on an edge

        line_check = 0

        if ((linex[1] == 0.0) & (linex[2] == 0.0)) & ((liney[1] == 0.0) & (liney[2] == 0.0))
            if linez[1] * linez[2] <= 0.0
                line_check = 1
            end
        elseif ((linex[1] == 0.0) & (linex[2] == 0.0)) & ((linez[1] == 0.0) & (linez[2] == 0.0))
            if liney[1] * liney[2] <= 0.0
                line_check = 1
            end
        elseif ((linez[1] == 0.0) & (linez[2] == 0.0)) & ((liney[1] == 0.0) & (liney[2] == 0.0))
            if linex[1] * linex[2] <= 0.0
                line_check = 1
            end
        end

        return line_check

    end

    function bilinear_cell(x::Float64, y::Float64, square::Array{Float64, 2})
        # interpolate the value of a function at (x, y) from the 4 corner values
        line = (1 - x)*square[1, :] + x*square[2, :]
        return line[1]*(1 - y) + line[2]*y
    end
        
end