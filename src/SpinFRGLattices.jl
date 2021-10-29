module SpinFRGLattices

using Parameters,StaticArrays,StructArrays,Test
version() = "0.3.6"

include("GeometryEssentials.jl")

include("Polymers.jl")
include("LatticeHelpers.jl")
include("ReduceGeometry.jl")
include("TestFunctions.jl")
include("Lattices/Lattices.jl")

end # module
