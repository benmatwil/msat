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

    function test()
        field = Common.read_field("data/bmag0037.dat")
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

    MSAT.NullFinder.NF(datafile)
    # MSAT.SeparatrixSurfaceFinder.SSF(datafile)
end