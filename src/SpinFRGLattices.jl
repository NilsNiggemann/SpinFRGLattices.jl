module SpinFRGLattices

using Parameters,StaticArrays,Test
version() = "0.3.3"

include("GeometryEssentials.jl")

include("Polymers.jl")
include("LatticeHelpers.jl")
include("ReduceGeometry.jl")
include("TestFunctions.jl")
include("Lattices/Lattices.jl")

end # module
