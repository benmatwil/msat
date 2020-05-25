module Common

	using StaticArrays
	using OffsetArrays
	using ..Params

	export Vector3D
	export AbstractField3D, CartesianField3D, SphericalField3D
	# common functions, subroutine and variables for all programs MSAT

    const Vector3D = SVector{3, Float64}

	abstract type AbstractField3D end

    struct CartesianField3D <: AbstractField3D
		filename::String
		field::Array{Float64, 4}
        x::Vector{Float64}
        y::Vector{Float64}
        z::Vector{Float64}
        xmin::Float64
        xmax::Float64
        ymin::Float64
        ymax::Float64
        zmin::Float64
        zmax::Float64
    end
	CartesianField3D(filename, field, x, y, z) = CartesianField3D(filename, field, x, y, z, 1, size(x, 1), 1, size(y, 1), 1, size(z, 1))

	struct SphericalField3D <: AbstractField3D
		filename::String
		field::Array{Float64, 4}
        x::Vector{Float64}
        y::Vector{Float64}
        z::Vector{Float64}
        xmin::Float64
        xmax::Float64
        ymin::Float64
        ymax::Float64
        zmin::Float64
        zmax::Float64
    end
	SphericalField3D(filename, field, x, y, z) = SphericalField3D(filename, field, x, y, z, 1, size(x, 1), 1, size(y, 1), 1, size(z, 1))

	# function filenames()
	# 	# sets the filenames as required for reading from and writing to

	# 	character(100) :: arg
	# 	character(:), allocatable :: outname
	# 	integer(int32) :: iarg, ic

	# 	ic = 0
	# 	outname = default_output

	# 	if (command_argument_count() > 0) then
	# 	do iarg = 1, command_argument_count()
	# 		call get_command_argument(iarg,arg)
	# 		if (trim(arg) == '-i') then
	# 		call get_command_argument(iarg+1,arg)
	# 		filein = trim(arg)
	# 		ic = 1
	# 		elseif (trim(arg) == '-o') then
	# 		call get_command_argument(iarg+1,arg)
	# 		deallocate(outname)
	# 		outname = trim(arg)
	# 		elseif (trim(arg) == '--params') then
	# 		call print_params
	# 		stop
	# 		endif
	# 	enddo
	# 	endif
	# 	if (ic == 0) stop 'No input file provided'
		
	# 	fileout = filein(1:index(filein(1:index(filein, '/', .true.)-1), '/', .true.)) &
	# 	//trim(outname)//'/' &
	# 	//trim(filein(index(filein, '/', .true.)+1:index(filein, '.dat', .true.)-1))
		
	# end

	function print_params()
		# allows the user to find out what parameters have been set in Params.jl when executable was compiled
		for name in propertynames(Params)
			println("$(name) = $(getproperty(Params, name))")
		end
	end

	# ********************************************************************************

	function trilinear(r::Vector3D, field::T) where T<:AbstractField3D
		# find the value of vector field at r using the trilinear method

		r = edgecheck(r, field)

		xp = r[1]
		yp = r[2]
		zp = r[3]

		nx = floor(Int, xp)
		ny = floor(Int, yp)
		nz = floor(Int, zp)

		# if point goes out of an edge which is not periodic, set point to be at boundary to trilinear
		# it will go out and be removed soon
		if outedge(r, field)
			if T == CartesianField3D
				if xp <= field.xmin
					xp = field.xmin
					nx = ceil(Int, xp)
				elseif xp >= field.xmax
					xp = field.xmax
					nx = ceil(Int, xp) - 1
				end
				if yp <= field.ymin
					yp = field.ymin
					ny = ceil(Int, yp)
				elseif yp >= field.ymax
					yp = field.ymax
					ny = ceil(Int, yp) - 1
				end
				if zp <= field.zmin
					zp = field.zmin
					nz = ceil(Int, zp)
				elseif zp >= field.zmax
					zp = field.zmax
					nz = ceil(Int, zp) - 1
				end
			elseif T == SphericalField3D
				if xp <= field.xmin
					xp = field.xmin
					nx = ceil(Int, xp)
				elseif xp >= field.xmax
					xp = field.xmax
					nx = ceil(Int, xp) - 1
				end
			end
		end

		x = xp - nx
		y = yp - ny
		z = zp - nz

        cube = SArray{Tuple{2, 2, 2, 3}, Float64}(field.field[nx:nx+1, ny:ny+1, nz:nz+1, :])
        square = (1 - z)*cube[:, :, 1, :] + z*cube[:, :, 2, :]
        line = (1 - y)*square[:, 1, :] + y*square[:, 2, :]
        return (1 - x)*line[1, :] + x*line[2, :]

	end

	function trilinear_nf(r::Vector3D, field::T) where T<:AbstractField3D
		# find the value of vector field at r using the trilinear method

		xp = r[1]
		yp = r[2]
		zp = r[3]

		nx = floor(Int, xp)
		ny = floor(Int, yp)
		nz = floor(Int, zp)

		if xp >= field.xmax
            xp = field.xmax
            nx = Int(xp - 1)
		end
		if yp >= field.ymax
            yp = field.ymax
            ny = Int(yp - 1)
		end
		if zp >= field.zmax
            zp = field.zmax
            nz = Int(zp - 1)
		end

		x = xp - nx
		y = yp - ny
		z = zp - nz

		# println("$nx, $ny, $nz, $xp, $yp, $zp")

        cube = SArray{Tuple{2, 2, 2, 3}, Float64}(field.field[nx:nx+1, ny:ny+1, nz:nz+1, :])
        square = (1 - z)*cube[:, :, 1, :] + z*cube[:, :, 2, :]
        line = (1 - y)*square[:, 1, :] + y*square[:, 2, :]
        return (1 - x)*line[1, :] + x*line[2, :]

	end

	# ********************************************************************************

	function dot(a::Vector3D, b::Vector3D)
		# dot product between a and b

		return sum(a .* b)

	end

	# ********************************************************************************

	function cross(a::Vector3D, b::Vector3D)
		# cross product between a and b

		return Vector3D([a[2] * b[3] - a[3] * b[2],
						 a[3] * b[1] - a[1] * b[3],
						 a[1] * b[2] - a[2] * b[1]])
        
	end

	# ********************************************************************************

	function normalise(a::Vector3D)
		# finds the unit vector of a

		return a / modulus(a)

	end

	# ********************************************************************************

	function modulus(a::Vector3D)
		# calculates the magnitude a vector a
		
		return sqrt(sum(a .^ 2))

	end

	# ********************************************************************************

	function dist(a::Vector3D, b::Vector3D)
		# calculates the distance between two points in grid units
		
		return sqrt(sum((b .- a) .^ 2))

	end

	# ********************************************************************************

	function outedge(r::Vector3D, field::T) where T<:AbstractField3D
		# determines if the point r is outwith the computation box across a non periodic
		# boundary and flags true or false

		if T == CartesianField3D
			if adjust_cartesian_periodicity
				out_cond = false
				if ! periodic_x
					out_cond = out_cond | (r[1] >= field.xmax) | (r[1] <= field.xmin)
				end
				if ! periodic_y
					out_cond = out_cond | (r[2] >= field.ymax) | (r[2] <= field.ymin)
				end
				if ! periodic_z
					out_cond = out_cond | (r[3] >= field.zmax) | (r[3] <= field.zmin)
				end
			else
				out_cond = (r[1] >= field.xmax) | (r[1] <= field.xmin) | (r[2] >= field.ymax) | (r[2] <= field.ymin) | (r[3] >= field.zmax) | (r[3] <= field.zmin)
			end
		elseif T == SphericalField3D
			out_cond = (r[1] >= field.xmax) | (r[1] <= field.xmin)
			if adjust_spherical_periodicity
				if ! periodic_theta
					out_cond = outedge | (r[2] >= field.ymax) | (r[2] <= field.ymin)
				end
				if ! periodic_phi
					out_cond = outedge | (r[3] >= field.zmax) | (r[3] <= field.zmin)
				end
			end
		end

        return out_cond
		
	end

	# ********************************************************************************

	function edgecheck(r::Vector3D, field::CartesianField3D)
		# determines if the point r is outwith the computational box across a 
		# periodic boundary and moves it if necessary

		if adjust_cartesian_periodicity
			rmut = MVector(r)
			if periodic_x
				if r[1] <= field.xmin
					rmut[1] += field.xmax - field.xmin
				end
				if r[1] >= field.xmax
					rmut[1] += field.xmin - field.xmax
				end
			end
			if periodic_y
				if r[2] <= field.ymin
					rmut[2] += field.ymax - field.ymin
				end
				if r[2] >= field.ymax
					rmut[2] += field.ymin - field.ymax
				end
			end
			if periodic_z
				if r[3] <= field.zmin
					rmut[3] += field.zmax - field.zmin
				end
				if r[3] >= field.zmax
					rmut[3] += field.zmin - field.zmax
				end
			end
			return SVector(rmut)
		end
		return r
	end

	function edgecheck(r::Vector3D, field::SphericalField3D)
		# determines if the point r is outwith the computational box across a 
		# periodic boundary and moves it if necessary

		if adjust_spherical_periodicity
			if periodic_theta
				if (r[2] < field.ymin) | (r[2] > field.ymax)
					if r[2] < field.ymin
						r[2] = 2*field.ymin - r[2]
					end
					if r[2] > field.ymax
						r[2] = 2*field.ymax - r[2]
					end
					if r[3] < (field.zmax + field.zmin)/2
						r[3] = r[3] + (field.zmax - field.zmin)/2
					else
						r[3] = r[3] - (field.zmax - field.zmin)/2
					end
				end
			end
			if periodic_phi
				if r[3] <= field.zmin
					r[3] = r[3] + (field.zmax - field.zmin)
				end
				if r[3] >= field.zmax
					r[3] = r[3] - (field.zmax - field.zmin)
				end
			end
			return r # this section needs sorting
		else
			rmut = MVector(r)
			if (r[2] < field.ymin) | (r[2] > field.ymax)
				if r[2] < field.ymin
					rmut[2] = 2*field.ymin - r[2]
				end
				if r[2] > field.ymax
					rmut[2] = 2*field.ymax - r[2]
				end
				if r[3] < (field.zmax + field.zmin)/2
					rmut[3] += (field.zmax - field.zmin)/2
				else
					rmut[3] += (field.zmin - field.zmax)/2
				end
			end
			if (r[3] <= field.zmin)
				rmut[3] += field.zmax - field.zmin
			end
			if (r[3] >= field.zmax)
				rmut[3] += (field.zmin - field.zmax)
			end
			return SVector(rmut)
		end
			
	end

	# ********************************************************************************

	function get_startpoints(fan::Vector3D, nlines::Integer)
		# from the theta, phi coordinates of a fan vector, produces ring of points in the fanplane

		theta = acos(fan[3])
        phi = atan(fan[2], fan[1])
		
		dtheta = 2*pi/nlines

		startpoints = Vector{Vector{Float64}}(undef, 0)

		# generate nlines start points in ring around equator
		for i in 1:nlines
			r = Float64[cos(i*dtheta), sin(i*dtheta), 0]
			r = rotate(r, theta, phi)
			push!(startpoints, r * start_dist)
		end

		return startpoints

	end

	# ********************************************************************************

	function rotate(r, theta::Float64, phi::Float64)
		# rotates a r by theta about the xz plane then by phi about the xy plane

		roty = Matrix{Float64}(undef, 3, 3)
		rotz = Matrix{Float64}(undef, 3, 3)

		roty[1, 1] = cos(-theta)
		roty[1, 2] = 0
		roty[1, 3] = -sin(-theta)
		
		roty[2, 1] = 0
		roty[2, 2] = 1
		roty[2, 3] = 0
		
		roty[3, 1] = sin(-theta)
		roty[3, 2] = 0
		roty[3, 3] = cos(-theta)

		rotz[1, 1] = cos(-phi)
		rotz[1, 2] = sin(-phi)
		rotz[1, 3] = 0

		rotz[2, 1] = -sin(-phi)
		rotz[2, 2] = cos(-phi)
		rotz[2, 3] = 0

		rotz[3, 1] = 0
		rotz[3, 2] = 0
		rotz[3, 3] = 1

		return (rotz * roty) * r

	end

	# ********************************************************************************

	function gtr(pt::Vector3D, field::AbstractField3D)
		# converts a list of grid points to real coordinates

		ig = floor.(Int, pt)

		return Vector3D([field.x[ig[1]] + (pt[1] - ig[1])*(field.x[ig[1] + 1] - field.x[ig[1]]),
						 field.y[ig[2]] + (pt[2] - ig[2])*(field.y[ig[2] + 1] - field.y[ig[2]]),
						 field.z[ig[3]] + (pt[3] - ig[3])*(field.z[ig[3] + 1] - field.z[ig[3]])])

	end

	# ********************************************************************************

	function file_position(nring::Integer, nperring::OffsetVector{Int32, Vector{Int32}}, index::Integer)
		# calculates the position in the temporary files of the rings and specific points

		# 1 whole ring contains 3*np vector points and np association points
		uptoring = sum(nperring[0:nring-1])*7
		a = uptoring + (index - 1)
		p = uptoring + nperring[nring] + (index - 1)*6

		a = a*4
		p = p*4

		return a, p

	end
	
end
