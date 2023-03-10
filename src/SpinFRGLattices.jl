module SpinFRGLattices

using StaticArrays,StructArrays,Test
version() = "0.3.10"

include("GeometryEssentials.jl")

include("Polymers.jl")
include("LatticeHelpers.jl")
include("ReduceGeometry.jl")
include("TestFunctions.jl")
include("Lattices/Lattices.jl")
# include("precompile.jl")
# __precompile__quiet__()
end # module
