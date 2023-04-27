module SpinFRGLattices

using StaticArrays,StructArrays,Test
version() = "0.4.0@matias"

include("GeometryEssentials.jl")
export sumElements,sitePair,Geometry,ArrayForm

include("Polymers.jl")
export getPolymer

include("LatticeHelpers.jl")
export  Basis_Struct_2D, Basis_Struct_3D,Basis_Struct, Rvec_2D, Rvec_3D, Rvec, getLatticeVec, translateToOrigin, translation, getCartesian, aboveLine, belowLine,aboveLine_strict, belowLine_strict, MirrorLine, dist, generatePairSites,generateLUnitCells,Mirror,printPairs

include("ReduceGeometry.jl")
export MapToPair, setCoupling!, setNeighborCouplings!, getSiteType, testGeometry, findSymmetryReduced, pairToInequiv_vec,generatePairToInequiv, getLatticeGeometry, getFRGComplexity

include("SpinSGeneralization.jl")
export adaptForSpinS,convertSusceptibilityToSpinS,convertSusceptibilityToSpinS!

include("TestFunctions.jl")
export testPairListSym,testPairListAdaptation


include("Lattices/Lattices.jl")

# include("precompile.jl")
# __precompile__quiet__()
end # module
