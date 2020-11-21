module Read

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
        bgrid_read = read!(fieldfile, Array{Float64, 4}(undef, nx, ny, nz, 3))
        bgrid = [Vector3D(bgrid_read[ix, iy, iz, :]...) for ix in 1:nx, iy in 1:ny, iz in 1:nz]
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
        gridpos = read!(nullfile, Vector{Vector3D{Float64}}(undef, nnulls))
        pos = read!(nullfile, Vector{Vector3D{Float64}}(undef, nnulls))
        close(nullfile)
        return nnulls, gridpos, pos
    end

    function read_nulldata(filename::AbstractString)
        nullfile = open(joinpath(default_output, prefix(filename) * "-nulldata.dat"), "r")
		nnulls, = read!(nullfile, Array{Int32}(undef, 1))
		signs = read!(nullfile, Array{Int32, 1}(undef, nnulls))
        spines = read!(nullfile, Vector{Vector3D{Float64}}(undef, nnulls))
        fans = read!(nullfile, Vector{Vector3D{Float64}}(undef, nnulls))
		warning = read!(nullfile, Vector{Int32}(undef, nnulls))
        close(nullfile)
        return nnulls, signs, spines, fans, warning
    end

end