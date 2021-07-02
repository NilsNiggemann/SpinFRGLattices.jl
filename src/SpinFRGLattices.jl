module SpinFRGLattices

using Parameters,StaticArrays,Test
version() = "0.3.1"

include("GeometryEssentials.jl")

include("Polymers.jl")
include("LatticeHelpers.jl")
include("Lattices/Lattices.jl")

function LatticePlots_activate()
    @eval Main begin
        using Parameters,StaticArrays,LaTeXStrings,Plots
        include(string($(@__DIR__),"/Fourier.jl"))
        include(string($(@__DIR__),"/LatticePlot.jl"))
        # export pairsPlot,plotSystem,plotCouplings
    end
end
export LatticePlots_activate

end # module
