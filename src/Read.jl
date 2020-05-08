module Read

    using StaticArrays

    using ..Params
    import ..Common

    function prefix(filename)
        fname = splitpath(filename)[end]
        dot = findlast('.', fname)
        return fname[1:dot-1]
    end

    function read_nulls(filename::AbstractString)
        nullfile = open(joinpath(default_output, prefix(filename) * "-nullpos.dat"), "r")
		nnulls, = read!(nullfile, Array{Int32}(undef, 1))
        gridpos = convert_point_list(read!(nullfile, Array{Float64, 2}(undef, 3, nnulls)))
        pos = convert_point_list(read!(nullfile, Array{Float64, 2}(undef, 3, nnulls)))
        close(nullfile)
        return nnulls, gridpos, pos
    end

    function read_nulldata(filename::AbstractString)
        nullfile = open(joinpath(default_output, prefix(filename) * "-nulldata.dat"), "r")
		nnulls, = read!(nullfile, Array{Int32}(undef, 1))
		signs = read!(nullfile, Array{Int32, 1}(undef, nnulls))
        spines = convert_point_list(read!(nullfile, Array{Float64, 2}(undef, 3, nnulls)))
        fans = convert_point_list(read!(nullfile, Array{Float64, 2}(undef, 3, nnulls)))
		warning = read!(nullfile, Array{Int32, 1}(undef, nnulls))
        close(nullfile)
        return nnulls, signs, spines, fans, warning
    end

    function convert_point_list(list::Array{Float64, 2})
        output = Vector{Common.Vector3D}(undef, 0)
        for index in 1:size(list, 2)
            push!(output, Common.Vector3D(list[:, index]))
        end
        return output
    end


end