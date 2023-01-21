"""
Includes essentials for arbitrary Geometry 
"""

export sumElements,sitePair,Geometry,getPolymer,ArrayForm
"""
Struct that contains all relevant information about the site summation.
Flow equations are of the form DV_ij = sum_k (V_ki V_kj S_k G_k) + V_ij V_ij S_i G_j 
"""


"""contains vertex indices for a single term in the site summation """
struct sumElements
    ki::Int #pair nr of left vertex
    kj::Int
    m::Int #multiplicity of this combination of vertices in site summation
    xk::Int
end

"""stores the types of site i and j in an inequivalent pair"""
struct sitePair
    xi::Int
    xj::Int
end

"""Geometry which contains all relevant information about the Lattice that is needed in the FRG
{RvecDim,SiteSumType <: StructArray{sumElements,2}}
RvecDim is type of Rvec, SiteSumType should be a struct array for performance reasons!
"""
struct Geometry{RvecDim,SiteSumType <: StructArray{sumElements,2},T}
    Name::String
    NLen::Int # System size: number of nearest neighbor pairs
    siteSum::SiteSumType #Matrix which contains all info about the site sum
    Npairs::Int
    Nsum::Vector{Int} #number of terms in the site sum for each Rij
    couplings::Vector{T} # Jij for each pair
    PairList::Vector{RvecDim} # list of inequivalent pairs
    invpairs::Vector{Int}# list of inverted pairs
    NUnique::Int # Number of unique sites in each cell
    PairTypes::Vector{sitePair} # gives pairs unique site indices according to PairList. Entries are always less or equal to NUnique
    OnsitePairs::Vector{Int} # positions of onsite pairs in PairList
end

function Geometry(;
    Name::String,
    NLen::Int = 1,
    siteSum,
    Npairs::Int = size(siteSum,2),
    Nsum::Vector{Int} = calcMaxpairs(siteSum),
    couplings::AbstractVector,
    PairList::Vector{RvecDim},
    invpairs::Vector{Int} = collect(Int,1:Npairs),
    PairTypes::Vector{sitePair} = [sitePair(1,1) for i in 1:Npairs],
    OnsitePairs::Vector{Int} = [1] ,
    NUnique::Int = length(OnsitePairs),
) where RvecDim
    return Geometry(Name,NLen,StructArray(siteSum),Npairs,Nsum,couplings,PairList,invpairs,NUnique,PairTypes,OnsitePairs)
end

"""Returns Vector containing the exact length of each row of siteSum so that site sum can be ended as soon as possible and we don't need to check if m==0 at performance critical parts."""
function calcMaxpairs(siteSum)
    Nsum,Npairs = size(siteSum)
    res = zeros(Int,Npairs)
    for j in 1:Npairs
        for k in 1:Nsum
            if siteSum[k,j].m == 0
                break
            end
            res[j] = k
        end
    end
    return res
end

Base.show(io::IO, ::MIME"text/plain", x::sumElements) = print(io,x.ki," ",x.kj," ",x.m," ",x.xk)
Base.show(io::IO, x::sumElements) = print(io,x.ki," ",x.kj," ",x.m," ",x.xk)
Base.show(io::IO, ::MIME"text/plain", x::sitePair) = print(io,x.xi,"_",x.xj)
Base.show(io::IO, x::sitePair) = print(io,x.xi,"_",x.xj)

# Base.show(io::IO, ::MIME"text/plain", siteSum::Matrix{sumElements}) = show(io,MIME"text/plain"(),ArrayForm(siteSum))

"""
For alternative printing/saving siteSum
"""
function ArrayForm(siteSum::Matrix{sumElements})
    Nsum,Npairs = size(siteSum)
    ArrayForm = Array{Int,3}(undef,4,Nsum,Npairs)
    for i in 1:Nsum, j in 1:Npairs
        x = siteSum[i,j]
        ArrayForm[:,i,j] = [x.ki,x.kj,x.m,x.xk]
    end
    ArrayForm
end

function Base.isless(x::sumElements, y::sumElements)
    if x.m == 0
        return false # move irreleveant splits to end of list when sorting
    elseif y.m == 0
        return true # move irreleveant splits to end of list when sorting
    elseif x.ki <y.ki # standard case
        return true
    elseif x.ki == y.ki && x.kj <y.kj
        return true
    else
        return false
    end
end