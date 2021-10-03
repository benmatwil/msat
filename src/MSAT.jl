module MSAT

    include("../Params.jl")
    using .Params
    include("Common.jl")
    using .Common
    include("Read.jl")
    include("Trace.jl")
    include("Ring.jl")
    include("NullFinder.jl")
    include("SpineFinder.jl")
    include("SeparatrixSurfaceFinder.jl")
    include("HeliosphericCurrentSheet.jl")
    include("MakeCut.jl")

    function test()
        field = Read.CartesianField3D("data/bmag0037.dat")
        middle = Common.Vector3D([10, 20, 30])
        starts = [middle + [cos(i), sin(i), 0] for i in (0:99)/100*2*pi]
        return [Trace.trace_line(start, 1, stepsize, field) for start in starts]
    end

end

if abspath(PROGRAM_FILE) == @__FILE__

    # field = MSAT.Common.read_field("data/bmag0037.dat")
    # middle = MSAT.Common.Vector3D([10, 20, 30])
    # starts = [middle + [cos(i), sin(i), 0] for i in (0:99)/100*2*pi]
    # pts = [MSAT.Trace.trace_line(start, 1, MSAT.stepsize, field) for start in starts]
    datafile = "data\\bmag0037-jl.dat"
    bgrid = MSAT.Common.CartesianField3D(datafile)
    bgrid = MSAT.Common.SphericalField3D("data\\field_20100601_0081_fft-jlv3d.dat")

    # MSAT.NullFinder.NF(datafile)
    # MSAT.SpineFinder.SF(datafile)
    MSAT.SeparatrixSurfaceFinder.SSF(bgrid)
    # MSAT.MakeCut.MC("data\\bmag0037-jl.dat", SA[0,0,1.], 2., do_hcs=false)
end