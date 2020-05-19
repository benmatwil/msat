module Read

    using StaticArrays

    using ..Params
    using ..Common

    function prefix(filename)
        fname = splitpath(filename)[end]
        dot = findlast('.', fname)
        return fname[1:dot-1]
    end

    function read_field(filename::String, coordinate_system::String)
        fieldfile = open(filename, "r")
		nx, ny, nz = read!(fieldfile, Array{Int32}(undef, 3))
		bgrid = read!(fieldfile, Array{Float64, 4}(undef, nx, ny, nz, 3))
		x = read!(fieldfile, Array{Float64, 1}(undef, nx))
		y = read!(fieldfile, Array{Float64, 1}(undef, ny))
		z = read!(fieldfile, Array{Float64, 1}(undef, nz))
        close(fieldfile)
        if coordinate_system == "Cartesian"
            return CartesianField3D(filename, bgrid, x, y, z)
        elseif coordinate_system == "Spherical"
            return SphericalField3D(filename, bgrid, x, y, z)
        else
            ArgumentError("Not a recognised coordinate system")
        end
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