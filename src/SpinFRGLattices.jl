module SpinFRGLattices

using StaticArrays, StructArrays, Test

include("GeometryEssentials.jl")
export sumElements, sitePair, Geometry, ArrayForm

include("Polymers.jl")
export getPolymer

include("LatticeHelpers.jl")
export Basis_Struct_2D, Basis_Struct_3D, Basis_Struct, Rvec_2D, Rvec_3D, Rvec, getUnitCell, getLatticeVec, translateToOrigin, translation, getCartesian, aboveLine, belowLine, aboveLine_strict, belowLine_strict, MirrorLine, norm, dist, generatePairSites, generateLUnitCells, generateSphere, Mirror, printPairs, isInUnitCell

include("ReduceGeometry.jl")
export MapToPair, setCoupling!, setNeighborCouplings!, getSiteType, PairNumbersDict, findSymmetryReduced, pairToInequiv_vec, generatePairToInequiv, getLatticeGeometry, getFRGComplexity, generateReducedLattice, ReducedLattice

include("SpinSGeneralization.jl")
export adaptForSpinS, convertSusceptibilityToSpinS, convertSusceptibilityToSpinS!

include("TestFunctions.jl")
export testPairListSym, testPairListAdaptation, testGeometry


include("Lattices/Lattices.jl")

# include("precompile.jl")
# __precompile__quiet__()
end # module
