module SpinFRGLattices

using Parameters,StaticArrays,Plots,Test
version() = "0.2.1"

include("GeometryEssentials.jl")

include("Polymers.jl")
include("Fourier.jl")
include("LatticeHelpers.jl")
include("Lattices/Lattices.jl")

end # module
