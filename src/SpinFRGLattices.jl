module SpinFRGLattices

using Parameters,StaticArrays,Test
version() = "0.3.2"

include("GeometryEssentials.jl")

include("Polymers.jl")
include("LatticeHelpers.jl")
include("Lattices/Lattices.jl")

end # module
