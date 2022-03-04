"""
Functions and data types that can be used to construct Lattices.
"""

export  Basis_Struct_2D, Basis_Struct_3D,Basis_Struct, Rvec_2D, Rvec_3D, Rvec, getLatticeVec, norm, translateToOrigin, translation, getCartesian, aboveLine, belowLine,aboveLine_strict, belowLine_strict, MirrorLine, dist, generatePairSites,generateLUnitCells


abstract type Rvec end

"""
Struct that stores the lattice indices n1, n2 as well as the basis index b = 1,2, ... 
"""
struct Rvec_2D <: Rvec
    n1::Int
    n2::Int
    b::Int
end    
"""
Struct that stores the lattice indices n1, n2 as well as the basis index b = 1,2, ... 
"""
struct Rvec_3D <: Rvec
    n1::Int
    n2::Int
    n3::Int
    b::Int
end

abstract type Basis_Struct end
"""
Basis of lattice
"""
@with_kw struct Basis_Struct_2D <:Basis_Struct
    a1::SVector{2,Float64} # lattice vector in cartesian coords
    a2::SVector{2,Float64}
    b::Vector{SVector{2, Float64}} # tuple of Basis vectors in cartesian coords
    T::SMatrix{2,2,Float64,4} =  inv([a1 a2]) # Transformation matrix as given by a1 and a2
    
    bLatt::Vector{SVector{2, Float64}} = [T * bv for bv in b] # tuple of Basis vectors in lattice coords
    NCell::Int = length(b)
    NNdist::Float64
    
    NUnique::Int = 1
    SiteType::Vector{Int} = [1 for i in 1:NCell]
    refSites::Vector{Rvec_2D} = [Rvec_2D(0,0,i) for i in 1:NUnique] # It is advisable to set the reference sites to be the first within the basis
end

@with_kw struct Basis_Struct_3D <:Basis_Struct
    a1::SVector{3,Float64} # lattice vector in cartesian coords
    a2::SVector{3,Float64}
    a3::SVector{3,Float64}
    b::Vector{SVector{3, Float64}} # tuple of Basis vectors in cartesian coords
    T::SMatrix{3,3,Float64,9} =  inv([a1 a2 a3]) # Transformation matrix as given by a1 and a2
    bLatt::Vector{SVector{3, Float64}} = [T * bv for bv in b] # tuple of Basis vectors in lattice coords

    NCell::Int = length(b)
    NNdist::Float64
    
    NUnique::Int = 1
    SiteType::Vector{Int} = [1 for i in 1:NCell]
    refSites::Vector{Rvec_3D}  =  [Rvec_3D(0,0,0,i) for i in 1:NUnique]
end

# printing form
# Base.show(io::IO, ::MIME"text/plain", R::Rvec_2D) = print(io,"[",R.n1,", ",R.n2,"]_",R.b)
# Base.show(io::IO, R::Rvec_2D) = print(io,"[",R.n1,", ",R.n2,"]_",R.b)

# Base.show(io::IO, ::MIME"text/plain", R::Rvec_3D) = print(io,"[",R.n1,", ",R.n2,", ",R.n3,"]_",R.b)
# Base.show(io::IO, R::Rvec_3D) = print(io,"[",R.n1,", ",R.n2,", ",R.n3,"]_",R.b)

"""Lazy constructor"""
function Rvec(n1,n2,b)
    return Rvec_2D(n1,n2,b)
end 

"""Lazy constructor"""
function Rvec(n1,n2,n3,b)
    return Rvec_3D(n1,n2,n3,b)
end

"""Get Type of site according to unique pairs"""
function getSiteType(R::Rvec,Basis)
    return Basis.SiteType[R.b]
end

function getCartesian(R::Rvec_2D,Basis)
    return R.n1 * Basis.a1 + R.n2 * Basis.a2 + Basis.b[R.b]
end

function getCartesian(R::Rvec_3D,Basis)
    return R.n1 * Basis.a1 + R.n2 * Basis.a2 + R.n3 * Basis.a3 + Basis.b[R.b]
end

function norm(r::AbstractArray)
    return( sqrt(sum(r .* r)))
end

function norm(R::Rvec,Basis)
    return( norm(getCartesian(R,Basis)))
end

R_isless(R1::Rvec,R2::Rvec,Basis) = norm(R1,Basis) < norm(R2,Basis) 

function dist(R1::Rvec,R2::Rvec,Basis)
    return norm(getCartesian(R1,Basis) .- getCartesian(R2,Basis))
end
    
"""Returns 2D or 3D vector in lattice coordinates"""
function getLatticeVec(r::StaticArray,Basis)::StaticArray
    # println(@allocated Basis.T*r)
    # println(@allocated SA[1. 0. 1.; 0. 1. 1.; 0. 1. 1.]*r)
    return Basis.T*r
end
"""Given a lattice vector returns a 2D or 3D vector in lattice coordinates"""
function getLatticeVec(R::Rvec,Basis)
    return getLatticeVec(getCartesian(R,Basis),Basis)
end

"""Translates Rj by -Rk """
function translateToOrigin(Rj::Rvec_2D,Rk::Rvec_2D) 
    Rvec(Rj.n1-Rk.n1 , Rj.n2-Rk.n2 , Rj.b)
end

function translateToOrigin(Rj::Rvec_3D,Rk::Rvec_3D) 
    Rvec(Rj.n1-Rk.n1 , Rj.n2-Rk.n2 , Rj.n3-Rk.n3 , Rj.b)
end

function translation(R::Rvec,T::AbstractArray,Basis)
    r = getCartesian(R,Basis) .+ T
    return getRvec(r,Basis)
end

"""
Converts cartesian coordinates to lattice struct Rvec. 
"""
function getRvec(r::StaticArray,Basis)
    r_lattice = Basis.T *r
    for (ib,b) in enumerate(Basis.bLatt)
        r_curr =  r_lattice .- b
        Intvec = round.(Int,r_curr)
        if all(abs2.(Intvec .- r_curr) .<1E-14) # r_lattice contains only Ints
            return(Rvec(Intvec...,ib))
        end
    end

    for b in Basis.bLatt
        r_curr =  r_lattice .- b
        Intvec = round.(Int,r_curr)
        println(abs2.(Intvec .- r_curr))
    end
    error("Vector not in lattice: ", r)
end

function aboveLine(f::Function,R::Rvec_2D,Basis)
    @warn "Function 'aboveLine' is deprecated, as results might be ambiguous. Please use 'aboveLine_strict' instead" maxlog=1
    r = getCartesian(R,Basis)
    y = f(r[1])
    return r[2] > y
end

function aboveLine_strict(f::Function,R::Rvec_2D,Basis)
    r = getCartesian(R,Basis)
    y = f(r[1])
    return r[2] > y && !(r[2] ≈ y)
end

function belowLine_strict(f::Function,R::Rvec_2D,Basis)
    r = getCartesian(R,Basis)
    y = f(r[1])
    return r[2] < y && !(r[2] ≈ y)
end

function MirrorLine(slope::Real,r::AbstractVector)
    l = SA[1,slope] # vector along line
    return 2* r' * l /(norm(l)^2) *l - r
end
function MirrorLine(slope::Real,R::Rvec_2D,Basis)
    mirr(R) = MirrorLine(slope,R)
    Cart(R) = getCartesian(R,Basis)
    RV(r) = getRvec(r,Basis)
    return R |> Cart |> mirr |> RV
end


"""Returns list of nearest neighbor pairs"""
function getNN(R::Rvec_2D,Basis)
    NN = []
    for n1 in -1:1,n2 in -1:1,b in 1:Basis.NCell
        Rnew = Rvec_2D(R.n1+n1,R.n2+n2,b)
        if isapprox(dist(R,Rnew,Basis),Basis.NNdist,atol = 1E-10)
            push!(NN,Rnew)
        end
    end
    return NN
end

"""Returns list of nearest neighbor pairs"""
function getNN(R::Rvec_3D,Basis)
    NN = []
    for n1 in -1:1,n2 in -1:1,n3 in -1:1, b in 1:Basis.NCell
        Rnew = Rvec_3D(R.n1+n1,R.n2+n2,R.n3+n3,b)
        if isapprox(dist(R,Rnew,Basis),Basis.NNdist,atol = 1E-10)
            push!(NN,Rnew)
        end
    end
    return NN
end

"""
Returns list of all pairs that are no more than N nearest neighbor pairs away from a specified reference site.
Symmetries may then be applied by deleting equivalent elements from PairList.
"""
function generatePairSites(N,Basis,refSite = Basis.refSites[1])
    siteList = [refSite]
    PairList = [refSite]
    CurrSiteList = []
    for _ in 1:N # procedurally generate sites going to nearest neighbors in each step
        empty!(CurrSiteList)
        for R in siteList 
            NN = getNN(R,Basis) #find all nearest neighbors
            for R_nn in NN
                if !(R_nn in PairList) #site has to be not already in the list
                    push!(CurrSiteList,R_nn)
                    push!(PairList,R_nn)
                end
            end
        end
        siteList = copy(CurrSiteList)
    end
    return PairList
end

"""
Returns list of all pairs within L unit cells.
Symmetries may then be applied by deleting equivalent elements from PairList.
"""
function generateLUnitCells(L,Basis::Basis_Struct_2D,refSite = Basis.refSites[1])
    # return
    PairList = Vector{Rvec_2D}(undef,0)
    n1ref,n2ref = refSite.n1,refSite.n2
    for n1 in -L:L, n2 in -L:L, b in 1:Basis.NCell
        push!(PairList,Rvec(n1+n1ref,n2+n2ref,b))
    end
    return PairList
end

function generateLUnitCells(L,Basis::Basis_Struct_3D,refSite = Basis.refSites[1])
    # return
    PairList = Vector{Rvec_3D}(undef,0)
    n1ref,n2ref,n3ref = refSite.n1,refSite.n2,refSite.n3
    for n1 in -L:L, n2 in -L:L, n3 in -L:L, b in 1:Basis.NCell
        push!(PairList,Rvec(n1+n1ref,n2+n2ref,n3+n3ref,b))
    end
    return PairList
end