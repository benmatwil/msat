module Trace
    using StaticArrays
    using ..Common
    using ..Params

	# rkf45 parameters
	const k21 = 0.25
    const k31 = 3/32
    const k32 = 9/32
    const k41 = 1932/2197
    const k42 = -7200/2197
    const k43 = 7296/2197
    const k51 = 439/216
    const k52 = -8
    const k53 = 3680/513
    const k54 = -845/4104
    const k61 = -8/27
    const k62 = 2
    const k63 = -3544/2565
    const k64 = 1859/4104
    const k65 = -11/40

    const y1 = 25/216
    const y3 = 1408/2565
    const y4 = 2197/4101
    const y5 = -1/5
    const z1 = 16/135
    const z3 = 6656/12825
    const z4 = 28561/56430
    const z5 = -9/50
    const z6 = 2/55

    function trace_line(r::Vector3D, sign::Integer, htotal::AbstractFloat, field::AbstractField3D)
    # traces a line from 'r' for 'nsteps' integration steps in the direction along the line as specified by 'sign'. Each step is of length h

        hdum = 0.0

        while hdum < htotal && !Common.outedge(r, field)
            h = sign*(htotal - hdum)
            r, h = rk45(r, h, field)
            hdum = hdum + abs(h)
        end

        return r

    end

    # ********************************************************************************

    function rk45(r::Vector3D, h::AbstractFloat, field::AbstractField3D)
    # runge-kutta fehlberg integrator. Calculates a 4th order estimate (y) and a
    # fifth order estimate (z) and calculates the difference between them. From this it
    # determines the optimal step length with which to integrate the function by.
        # real(np) :: hvec(3), mindist, s
        # real(np), dimension(3) :: k1, k2, k3, k4, k5, k6
        # real(np), dimension(3) :: r, r0, rtest, y, z, dr
        
        r0 = copy(r)
        
        # minimum physical distance corresponding to fraction of gridcell (h)
        dr = getdr(r0, field)
        mindist = minimum(dr)*h
        
        # vector containing the grid's h for each direction
        # (so h for each direction is the same physical length, equal to mindist)
        hvec = mindist ./ dr

        # get rk values k1--k6
        rtest = r0
        k1 = hvec .* Common.normalise(Common.trilinear(rtest, field))
        rtest = r0 + k21*k1
        k2 = hvec .* Common.normalise(Common.trilinear(rtest, field))
        rtest = r0 + k31*k1 + k32*k2
        k3 = hvec .* Common.normalise(Common.trilinear(rtest, field))
        rtest = r0 + k41*k1 + k42*k2 + k43*k3
        k4 = hvec .* Common.normalise(Common.trilinear(rtest, field))
        rtest = r0 + k51*k1 + k52*k2 + k53*k3 + k54*k4
        k5 = hvec .* Common.normalise(Common.trilinear(rtest, field))
        rtest = r0 + k61*k1 + k62*k2 + k63*k3 + k64*k4 + k65*k5
        k6 = hvec .* Common.normalise(Common.trilinear(rtest, field))

        # get 4th order (y) and 5th order (z) estimates
        y = y1*k1 + y3*k3 + y4*k4 + y5*k5
        z = z1*k1 + z3*k3 + z4*k4 + z5*k5 + z6*k6

        # calculate optimum step length (s = hoptimum/h)
        s = (tol ./ Common.modulus(z - y) ./ 2) .^ 0.25
        
        if (abs(s*h) < stepmin) 
            s = stepmin/abs(h)
        end
        if (s > 1)
            s = 1.0
        end

        hvec = s*hvec

        rtest = r0
        k1 = hvec .* Common.normalise(Common.trilinear(rtest, field))
        rtest = r0 + k21*k1
        k2 = hvec .* Common.normalise(Common.trilinear(rtest, field))
        rtest = r0 + k31*k1 + k32*k2
        k3 = hvec .* Common.normalise(Common.trilinear(rtest, field))
        rtest = r0 + k41*k1 + k42*k2 + k43*k3
        k4 = hvec .* Common.normalise(Common.trilinear(rtest, field))
        rtest = r0 + k51*k1 + k52*k2 + k53*k3 + k54*k4
        k5 = hvec .* Common.normalise(Common.trilinear(rtest, field))

        r = r0 + y1*k1 + y3*k3 + y4*k4 + y5*k5

        h = h*s

        return r, h

    end

    # ********************************************************************************

    function getdr(r::Vector3D, field::T) where T<:AbstractField3D
        # outputs length of one gridcell in 'physical' length units (essentially (dx,dy,dz))
        # integer(int32) :: ix, iy, iz

        rcheck = copy(r)
        Common.edgecheck(rcheck, field)

        ix = floor(Int, rcheck[1])
        iy = floor(Int, rcheck[2])
        iz = floor(Int, rcheck[3])
        
        dx = field.x[ix+1] - field.x[ix]
        dy = field.y[iy+1] - field.y[iy]
        dz = field.z[iz+1] - field.z[iz]

        if T == CartesianField3D
            return SA[dx, dy, dz]
        elseif T == SphericalField3D
            xp = field.x[ix] + (rcheck[1] - ix)*dx
            yp = field.y[iy] + (rcheck[2] - iy)*dy
            return SA[dx, xp*dy, xp*sin(yp)*dz]
        end

    end

end
