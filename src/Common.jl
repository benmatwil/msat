module Common

	import Base.+, Base.-, Base.*, Base./

	using OffsetArrays
	using ..Params

	export Vector3D
	export AbstractField3D, CartesianField3D, SphericalField3D, CylindricalField3D
	# common functions, subroutine and variables for all programs MSAT

	struct Vector3D{T<:Number}
		x::T
		y::T
		z::T
	end
	Vector3D(r) = Vector3D(r...)

	+(a::Vector3D, b::Vector3D) = Vector3D(a.x + b.x,
										   a.y + b.y,
										   a.z + b.z)

	-(a::Vector3D, b::Vector3D) = Vector3D(a.x - b.x,
										   a.y - b.y,
										   a.z - b.z)

	*(a::Number, b::Vector3D) = Vector3D(a * b.x,
										 a * b.y,
										 a * b.z)
	*(a::Vector3D, b::Number) = b * a
	
	/(a::Vector3D, b::Number) = Vector3D(a.x / b,
										 a.y / b,
										 a.z / b)
	/(a::Number, b::Vector3D) = Vector3D(a / b.x,
										 a / b.y,
										 a / b.z)

	*(a::Vector3D, b::Vector3D) = Vector3D(a.x * b.x,
										   a.y * b.y,
										   a.z * b.z)

	×(a::Vector3D, b::Vector3D) = Vector3D(a.y * b.z - a.z * b.y,
										   a.z * b.x - a.x * b.z,
										   a.x * b.y - a.y * b.x)

	function dot(a::Vector3D, b::Vector3D)

		return a.x * b.x + a.y * b.y + a.z * b.z

	end

	modulus(a::Vector3D) = sqrt(dot(a, a))
	normalise(a::Vector3D) = a / modulus(a)
	
	function dist(a::Vector3D, b::Vector3D)
		# calculates the distance between two points in grid units

		diff = (b - a)
		return sqrt(diff.x ^ 2 + diff.y ^ 2 + diff.z ^ 2)
	
	end
	
	Base.minimum(a::Vector3D{T}) where T = min(a.x, a.y, a.z)

	Base.zero(::Type{Vector3D{T}}) where T = Vector3D(zero(T), zero(T), zero(T))

	Base.read(io::IO, ::Type{Vector3D{T}}) where T = Vector3D{T}(read(io, T),read(io, T),read(io, T))
	Base.write(io::IO, v::Vector3D{T}) where T = write(io, v.x) + write(io, v.y) + write(io, v.z)

	abstract type AbstractField3D end
	Base.show(io::IO, a::AbstractField3D) = print(io, "$(typeof(a))($(a.filename))")

    struct CartesianField3D <: AbstractField3D
		filename::String
		field::Array{Vector3D{Float64}, 3}
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
		field::Array{Vector3D{Float64}, 3}
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

	struct CylindricalField3D <: AbstractField3D
		filename::String
		field::Array{Vector3D{Float64}, 3}
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
	CylindricalField3D(filename, field, x, y, z) = CylindricalField3D(filename, field, x, y, z, 1, size(x, 1), 1, size(y, 1), 1, size(z, 1))

	function print_params()
		# allows the user to find out what parameters have been set in Params.jl when executable was compiled
		for name in propertynames(Params)
			println("$(name) = $(getproperty(Params, name))")
		end
	end

	# ********************************************************************************

	function trilinear(r::Vector3D, field::T) where T<:AbstractField3D
		# find the value of vector field at r using the trilinear method

		xp, yp, zp = edgecheck(r.x, r.y, r.z, field)

		nx = floor(Int, xp)
		ny = floor(Int, yp)
		nz = floor(Int, zp)

		# if point goes out of an edge which is not periodic, set point to be at boundary to trilinear
		# it will go out and be removed soon
		if outedge(xp, yp, zp, field)
			if T == CartesianField3D
				if xp < field.xmin
					xp = field.xmin
					nx = Int(field.xmin)
				elseif xp >= field.xmax
					xp = field.xmax
					nx = Int(field.xmax) - 1
				end
				if yp < field.ymin
					yp = field.ymin
					ny = Int(field.ymin)
				elseif yp >= field.ymax
					yp = field.ymax
					ny = Int(field.ymax) - 1
				end
				if zp < field.zmin
					zp = field.zmin
					nz = Int(field.zmin)
				elseif zp >= field.zmax
					zp = field.zmax
					nz = Int(field.zmax) - 1
				end
			elseif T == SphericalField3D
				if xp < field.xmin
					xp = field.xmin
					nx = Int(field.xmin)
				elseif xp >= field.xmax
					xp = field.xmax
					nx = Int(field.xmax) - 1
				end
			end
		end

		x = xp - nx
		y = yp - ny
		z = zp - nz

        return tri(field, nx, ny, nz, x, y, z)

	end

	function tri(field, ix, iy, iz, x, y, z)
		
		b000 = field.field[ix, iy, iz]
		b100 = field.field[ix+1, iy, iz]
		b010 = field.field[ix, iy+1, iz]
		b110 = field.field[ix+1, iy+1, iz]
		b001 = field.field[ix, iy, iz+1]
		b101 = field.field[ix+1, iy, iz+1]
		b011 = field.field[ix, iy+1, iz+1]
		b111 = field.field[ix+1, iy+1, iz+1]

		a = b000
		b = b100 - b000
		c = b010 - b000
		d = b110 - b100 - b010 + b000
		e = b001 - b000
		f = b101 - b100 - b001 + b000
		g = b011 - b010 - b001 + b000
		h = b111 - b110 - b101 - b011 + b100 + b010 + b001 - b000        

		return a + b*x + c*y + d*x*y + e*z + f*x*z + g*y*z + h*x*y*z
	
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

	outedge(r::Vector3D, field::T) where T<:AbstractField3D = outedge(r.x, r.y, r.z, field)

	function outedge(rx, ry, rz, field::CartesianField3D)
		# determines if the point r is outwith the computation box across a non periodic
		# boundary and flags true or false
		if adjust_cartesian_periodicity
			out_cond = false
			if ! periodic_x
				out_cond = out_cond | (rx >= field.xmax) | (rx <= field.xmin)
			end
			if ! periodic_y
				out_cond = out_cond | (ry >= field.ymax) | (ry <= field.ymin)
			end
			if ! periodic_z
				out_cond = out_cond | (rz >= field.zmax) | (rz <= field.zmin)
			end
		else
			out_cond = (rx >= field.xmax) | (rx <= field.xmin) | (ry >= field.ymax) | (ry <= field.ymin) | (rz >= field.zmax) | (rz <= field.zmin)
		end
		return out_cond
	end


	function outedge(rx, ry, rz, field::SphericalField3D)
		out_cond = (rx >= field.xmax) | (rx <= field.xmin)
		if adjust_spherical_periodicity
			if ! periodic_theta
				out_cond = outedge | (ry >= field.ymax) | (ry <= field.ymin)
			end
			if ! periodic_phi
				out_cond = outedge | (rz >= field.zmax) | (rz <= field.zmin)
			end
		end
		return out_cond
	end

	# ********************************************************************************

	function edgecheck(rx, ry, rz, field::CartesianField3D)
		# determines if the point r is outwith the computational box across a
		# periodic boundary and moves it if necessary

		if adjust_cartesian_periodicity
			if periodic_x
				if rx < field.xmin
					rx += field.xmax - field.xmin
				end
				if rx >= field.xmax
					rx += field.xmin - field.xmax
				end
			end
			if periodic_y
				if ry < field.ymin
					ry += field.ymax - field.ymin
				end
				if ry >= field.ymax
					ry += field.ymin - field.ymax
				end
			end
			if periodic_z
				if rz < field.zmin
					rz += field.zmax - field.zmin
				end
				if rz >= field.zmax
					rz += field.zmin - field.zmax
				end
			end
		end
		return rx, ry, rz
	end

	edgecheck(r::Vector3D, field::AbstractField3D) = Vector3D(edgecheck(r.x, r.y, r.z, field)...)

	function edgecheck(rx, ry, rz, field::SphericalField3D)
		# determines if the point r is outwith the computational box across a
		# periodic boundary and moves it if necessary

		if adjust_spherical_periodicity
			if periodic_theta
				if (ry < field.ymin) | (ry > field.ymax)
					if ry < field.ymin
						ry = 2*field.ymin - ry
					end
					if ry > field.ymax
						ry = 2*field.ymax - ry
					end
					if rz < (field.zmax + field.zmin)/2
						rz = rz + (field.zmax - field.zmin)/2
					else
						rz = rz - (field.zmax - field.zmin)/2
					end
				end
			end
			if periodic_phi
				if rz <= field.zmin
					rz = rz + (field.zmax - field.zmin)
				end
				if rz >= field.zmax
					rz = rz - (field.zmax - field.zmin)
				end
			end
			return r # this section needs sorting
		else
			if (ry < field.ymin) | (ry > field.ymax)
				if ry < field.ymin
					ry = 2*field.ymin - ry
				end
				if ry > field.ymax
					ry = 2*field.ymax - ry
				end
				if rz < (field.zmax + field.zmin)/2
					rz += (field.zmax - field.zmin)/2
				else
					rz += (field.zmin - field.zmax)/2
				end
			end
			if rz < field.zmin
				rz += field.zmax - field.zmin
			end
			if rz >= field.zmax
				rz += field.zmin - field.zmax
			end
		end

		return rx, ry, rz
	end

	# ********************************************************************************

	function get_startpoints(fan::Vector3D, nlines::Integer)
		# from the theta, phi coordinates of a fan vector, produces ring of points in the fanplane

		theta = acos(fan.z)
        phi = atan(fan.y, fan.x)

		dtheta = 2*pi/nlines

		startpoints = Vector{Vector3D{Float64}}(undef, 0)

		# generate nlines start points in ring around equator
		for i in 1:nlines
			r = Float64[cos(i*dtheta), sin(i*dtheta), 0]
			r = Vector3D(rotate(r, theta, phi)...)
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

		igx = floor(Int, pt.x)
		igy = floor(Int, pt.y)
		igz = floor(Int, pt.z)

		return Vector3D(field.x[igx] + (pt.x - igx)*(field.x[igx + 1] - field.x[igx]),
						field.y[igy] + (pt.y - igy)*(field.y[igy + 1] - field.y[igy]),
						field.z[igz] + (pt.z - igz)*(field.z[igz + 1] - field.z[igz]))

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

	function get_rnullsalt(rnulls::Vector{Vector3D{Float64}}, bgrid::SphericalField3D)

		rnullsalt = copy(rnulls)
		
		for (inull, rnull) in enumerate(rnullsalt)

			rx = rnull.x
			ry = rnull.y
			rz = rnull.z

			# check whether null is at the lower phi boundary
			if rz < bgrid.zmin + 1
				rz = rz - 1
				rx, ry, rz = edgecheck(rx, ry, rz, bgrid)
				rz = rz + 1
			end

			# check whether null is at the upper phi boundary
			if rz > bgrid.zmax - 1
				rz = rz + 1
				rx, ry, rz = edgecheck(rx, ry, rz, bgrid)
				rz = rz - 1
			end

			# check whether null is at the lower theta boundary
			if ry < bgrid.ymin + 1
				ry = ry - 1
				rx, ry, rz = edgecheck(rx, ry, rz, bgrid)
				ry = ry - 1
			end

			# check whether null is at the upper theta boundary
			if ry > bgrid.ymax - 1
				ry = ry + 1
				rx, ry, rz = edgecheck(rx, ry, rz, bgrid)
				ry = ry + 1
			end

			rnullsalt[inull] = Vector3D(rx, ry, rz)

		end

        return rnullsalt

    end

	function get_rnullsalt(rnulls::Vector{Vector3D}, bgrid::CylindricalField3D)
		
		rnullsalt = copy(rnulls)
		
		for (inull, rnull) in enumerate(rnullalt)

			rx = rnull.x
			ry = rnull.y
			rz = rnull.z

			# check whether null is at the lower phi boundary
			if ry < bgrid.ymin + 1
				ry = ry - 1
				rx, ry, rz = edgecheck(rx, ry, rz, bgrid)
				ry = ry + 1
			end

			# check whether null is at the upper phi boundary
			if ry < bgrid.ymin + 1
				ry = ry + 1
				rx, ry, rz = edgecheck(rx, ry, rz, bgrid)
				ry = ry - 1
			end

			rnullsalt[inull] = Vector3D(rx, ry, rz)

		end
        
        return rnullsalt

    end
	
end
