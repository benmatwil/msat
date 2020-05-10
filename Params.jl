module Params
    export MSAT_float, nproc, default_output
    export zero_value, sig_figs, boundary_nulls, gridpoint_nulls
    export rspherefact, nphi, ntheta, maxiter, dist_mult
    export nstart, start_dist, ringsmax, pointsmax, samemax, stepsize, tol, stepmin, restart, nullrestart, nskip, bytesize, assoc_output, one_sep_per_ring
    export adjust_cartesian_periodicity, adjust_cylindrical_periodicity, adjust_spherical_periodicity, periodic_x, periodic_y, periodic_z, periodic_theta, periodic_phi

    # Global Parameters
    # -----------------
    const MSAT_float = Float64
    # Number of threads to be used by OpenMP
    # set to zero to use the environment variable OMP_NUM_THREADS
    const nproc = 4
    # The default directory name to output data to
    const default_output = "output"

    # Null Finder Parameters
    # ----------------------
    # what the code treats as zero
    # for large regions of weak field, set zero to be a small (but larger) number i.e. 1e-12_np
    const zero_value = 1e-16
    # the number of significant figures of accuracy required for the null
    const sig_figs = 8
    # whether to include null points close to the boundary of vector field grids
    # if true sign finder may struggle to analyse them
    const boundary_nulls = false
    # whether to check the grid points for null points - may be useful for analytical grids
    const gridpoint_nulls = true

    # Sign Finder Parameters
    # ----------------------
    # radius factor (in grid coordinates) of sphere on which to 
    # place start points (rsphere = rspherefact*10**-sig_figs)
    const rspherefact = 100.0
    # number of start points in phi direction
    const nphi = 50
    # number of start points in theta direction
    const ntheta = nphi ÷ 2
    # maximum number of iterations for convergence method
    const maxiter = 20000
    # multiple of rsphere distance to move points at each iteration
    const dist_mult = 2e-2

    # Separatrix Surface Finder Parameters
    # (Most of these also apply to HCS)
    # ------------------------------------
    # number of startpoints in ring
    const nstart = 500
    # distance from null to start ring
    const start_dist = 0.01
    # maximum number of rings
    const ringsmax = 50000
    # maximum number of points in ring
    const pointsmax = 200000
    # number of consectutive rings which are the same size which stop ssf
    const samemax = 1000
    # step size h after 50 iterations (otherwise 5 times smaller)
    const stepsize = 0.25
    # tolerance of rkf45 scheme
    const tol = 1e-6
    # minimum step length
    const stepmin = 1e-5
    # turn on restart function of the code
    const restart = false
    # which null to restart from (make sure it's correct)
    # set to 0 to let the code decide using already written data
    const nullrestart = 0
    # number of rings to skip in output to file
    # set to 1 to write all rings to file
    const nskip = 1
    # output bitsize of floating points numbers for rings
    # default real64 (double precision)
    # also allows make_cut to read points in correctly
    const bytesize = Float64
    # whether ssf outputs associations
    # associations only currently required for make_cut
    const assoc_output = true
    # allow changes in default periodicity
    # default is no periodicity in cartesian and
    # periodic theta, phi in cylindrical/spherical
    const adjust_cartesian_periodicity = false
    const adjust_cylindrical_periodicity = false
    const adjust_spherical_periodicity = false
        # turn on and off periodicity in x, y, z, theta and phi
        # these do nothing if adjust_***_periodicities is false
        const periodic_x = false
        const periodic_y = false
        const periodic_z = false
        const periodic_theta = false
        const periodic_phi = false
    # Only allow one separator to be found per null point per ring
    # This can help to reduce the number of separators found in complex fields
    const one_sep_per_ring = false
end
