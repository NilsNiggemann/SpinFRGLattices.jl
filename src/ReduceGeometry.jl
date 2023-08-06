"""Gives a sorted list of pairs for all reference sites together with a list which types of sites are paired"""
function sortedPairList(N,Basis,method = generatePairSites)
    type = typeof(Basis.refSites[1])
    sortedPairs = Vector{type}(undef,0)
    sortedPairTypes = Vector{sitePair}(undef,0)
    
    refSites = Basis.refSites
    sortingfuncs = Tuple( (R1,R2) -> dist(R1,Ri,Basis) < dist(R2,Ri,Basis) for Ri in refSites)

    for (sortingfunc,refSite) in zip(sortingfuncs,refSites)
        Pairs =  method(N,Basis,refSite)
        PairTypes = [sitePair(getSiteType(refSite,Basis),getSiteType(R,Basis)) for R in Pairs]
        
        perm = sortperm(Pairs,lt= sortingfunc) 
        
        # We sort the Pairs corresponding to each ref site by norm. 
        # Be careful when switching the sorting type that siteType may be sorted accordingly using perm[Pairs]
        append!(sortedPairTypes,PairTypes[perm])
        append!(sortedPairs,Pairs[perm]) # since the array PairTypes is 
        
    end
    return sortedPairs,sortedPairTypes
end

"""Uses passed function inCorrectSector (R_ref,R) -> bool to generate a new minimal list of pairs that are not equivalent by symmetry"""
function getinequivalentPairs(PairList,PairTypes,inCorrectSector::Function,Basis)
    inds = []

    for (i,R) in enumerate(PairList)
        R_Ref = Basis.refSites[PairTypes[i].xi]
        if inCorrectSector(R_Ref,R)
            push!(inds,i)
        end
    end
    return PairList[inds],PairTypes[inds] # only select elements that are in correct subsector
end

"""
maps pair of sites to its corresponding inequivalent pair in PairList. x = 1,2 corresponds to index of reference site.
Returns 0 if no pair is found.
"""
function MapToPair(x::Integer,R::RT,PairList::AbstractVector{RT},PairTypes) where RT <:Rvec
    @inbounds for (i,b) in enumerate(PairList)
        if b === R && PairTypes[i].xi === x 
            return(i)
        end
    end
    # error("error: pair not found")
    return 0 # no pair found
    # println("cannot map: ",R," to inequivalent pair. Not found in PairList for site ",x)
end

function getMapToPairDict(PairList,PairTypes)
    Dict((PairTypes[i].xi,PairList[i]) => i for i in eachindex(PairList))
end

"""Sets coupling between ref site x and R to val. To map from x and R to the corresponding pair, we need PairList and siteTypes"""
function setCoupling!(couplings,x::Int,R::Rvec,val::Number,PairList,PairTypes)
    idx = MapToPair(x,R,PairList,PairTypes)
    if idx <=0
        println(PairList)
        println(PairTypes)
        error("problem setting coupling for (x,R,val) = ",(x,R,val))
    end
    couplings[idx] = val
end

"""Sets coupling between ref site x and R to val. Takes full Geometry as input"""

function setCoupling!(System::Geometry,x::Int,R::Rvec,val::Number)
    setCoupling!(System.couplings,x,R,val,System.PairList,System.PairTypes)
end

function getDists(SiteList,Basis)
    dists = [dist(i,j,Basis) for i in SiteList for j in SiteList] |> unique! |> sort!
    return dists
end

"""Sets couplings according to given list with the first element being the nearest neigbor coupling."""
function setNeighborCouplings!(couplings,Jparams::AbstractVector,PairList,PairTypes,Basis;tol = 1e-13)
    dists = getDists(PairList,Basis)
    filter!(>(0),dists)
    @assert length(Jparams) <= length(dists) "System is to small to set couplings"

    couplings .= 0.

    for (i,(Rj,type)) in enumerate(zip(PairList,PairTypes,couplings))
        R0 = Basis.refSites[type.xi]

        d = dist(R0,Rj,Basis)
        
        idx = findfirst(x -> abs(x - d) < tol, dists)
        if idx !== nothing && idx <= length(Jparams)
            couplings[i] = Jparams[idx]
        end
    end

    return couplings
end

setNeighborCouplings!(System,Jparams::AbstractVector,Basis) = setNeighborCouplings!(System.couplings,Jparams,System.PairList,System.PairTypes,Basis)

"""gives multiplicity of spl in list of sumelements """
function multiplicity(spl::sumElements,list)
    mult = 0
    for i in eachindex(list)
        if spl.ki == list[i].ki && spl.kj == list[i].kj && spl.xk == list[i].xk
            mult+=1
        end
    end
    return mult
end

"""gives multiplicity of spl in SORTED list of sumelements """
function multiplicity_sorted(spl::sumElements,list)
    mult = 0
    @inbounds for i in eachindex(list)
        if list[i].ki > spl.ki
            break
        end
        if spl.ki == list[i].ki && spl.kj == list[i].kj && spl.xk == list[i].xk
            mult+=1
        end
    end
    return mult
end

"""Gives number of pair in PairList corresponding to pair R1,R2.
Function pairToInequiv must map two sites to a pair containing a reference site """
function pairNumber(R1::R,R2::R,PairList::AbstractVector{R},PairTypes::AbstractVector{sitePair}, Basis::Basis_Struct,pairToInequiv::Function) where {R<:Rvec}
    R1,R2 = pairToInequiv(R1,R2)
    return MapToPair(getSiteType(R1,Basis),R2,PairList,PairTypes)
end 
function pairNumber(R1::R,R2::R,mapToPairDict::AbstractDict,nonRefSyms,refSyms,Basis::Basis_Struct) where {R<:Rvec}
    R1,R2 = pairToInequiv_vec(R1,R2,Basis,mapToPairDict,nonRefSyms,refSyms)
    x = getSiteType(R1,Basis)
    if (x,R2) in keys(mapToPairDict)
        return mapToPairDict[(x,R2)]
    end
    return 0
end

"""Given a Dict that maps a pair to the corresponding inequivalent pair index `i` returns `i` if the pair can be mapped correctly. Otherwise, return `0`
It is assumed that global translations along lattice vectors can always be performed.
"""
struct PairNumbersDict{R<:Rvec} <: AbstractDict{Tuple{R,R},Int}
    D::Dict{Tuple{R,R},Int}
end
Base.keys(P::PairNumbersDict) = keys(P.D)
Base.values(P::PairNumbersDict) = values(P.D)
Base.length(P::PairNumbersDict) = length(P.D)
Base.iterate(P::PairNumbersDict) = iterate(P.D)
Base.iterate(P::PairNumbersDict,i) = iterate(P.D,i)
Base.get(P::PairNumbersDict,i,default) = get(P.D,i,default)

function Base.getindex(P::PairNumbersDict{R},pair::Tuple{R,R}) where {R<:Rvec}
    pair = translatePairToUnitCell(pair)
    
    return get(P.D,pair,0)    
end

function pairNumbersDict(siteList::AbstractVector{R},PairList::AbstractVector{R},PairTypes,nonRefSyms,refSyms,Basis) where {R<:Rvec}
    MapToPairDict = getMapToPairDict(PairList,PairTypes)
    
    pairs = collect(Set(translatePairToUnitCell(R1,R2) for R1 in siteList for R2 in siteList))

    pairNumberVec = Vector{Pair{Tuple{R,R},Int}}(undef,length(pairs))

    Threads.@threads for i in eachindex(pairNumberVec,pairs)
        R1,R2 = pairs[i]
        pairNumberVec[i] = (R1,R2) => pairNumber(R1,R2,MapToPairDict,nonRefSyms,refSyms,Basis)
    end

    pairNumberDict = Dict(pairNumberVec)

    filter!(x-> x.second != 0,pairNumberDict) # remove all pairs that are not in PairList
    return PairNumbersDict(pairNumberDict)
end

"""This function is for backwards compatibility. It is recommended to use the function with the same name that uses the mapToPairDict instead of the pairToInequiv function"""
function pairNumbersDict(siteList,PairList,PairTypes,pairToInequiv::Function,Basis)
    pairNumber_((R1,R2)) = pairNumber(R1,R2,PairList,PairTypes,Basis,pairToInequiv)
    pairs = Set(translatePairToUnitCell(R1,R2) for R1 in siteList for R2 in siteList)
    pairNumberDict = Dict([ pair => pairNumber_(pair) for pair in pairs]) # store the pair number of all possible pairs in a dictionary
    filter!(x-> x.second != 0,pairNumberDict) # remove all pairs that are not in PairList

    return PairNumbersDict(pairNumberDict)
end


"""Taking an inequivalent pair (x,R) with x <= NUnique returns corresponding inversed pair (R,x)"""
function inversepair(R_ref::Rvec,R::Rvec,PairList,PairTypes,Basis,pairToInequiv::Function)
    return pairNumber(R,R_ref,PairList,PairTypes,Basis,pairToInequiv)
end

inversepair(R_ref::Rvec,R::Rvec,pairNumbersDict::PairNumbersDict) = pairNumbersDict[(R,R_ref)]
"""
Performs site summation for each Rj (lhs) of flow equations and returns all pairs that appear in it corresponding to pairs of vertices on rhs.
It is assumed that our reference sites are the first in the basis!
"""
function CalcSiteSum(PairList::AbstractVector{<:Rvec},PairTypes::AbstractVector{sitePair},siteList::AbstractVector{<:Rvec},pairToInequiv::Function,Basis::Basis_Struct)
    pairNumber_(R1,R2) = pairNumber(R1,R2,PairList,PairTypes,Basis,pairToInequiv)
    # pairNumberDict = Dict([ (R1,R2) => pairNumber_(R1,R2) for R1 in siteList for R2 in PairList]) # store the pair number of all possible pairs in a dictionary
    pairNumberDict = pairNumbersDict(siteList,PairList,PairTypes,pairToInequiv,Basis)
    return CalcSiteSum(PairList,PairTypes,siteList,pairNumberDict,Basis)
end

function CalcSiteSum(PairList::AbstractVector{<:Rvec},PairTypes::AbstractVector{sitePair},siteList::AbstractVector{<:Rvec},syms::Tuple,Basis::Basis_Struct)
    (nonRefSyms,refSyms) = syms
    pairNumberDict = pairNumbersDict(siteList,PairList,PairTypes,nonRefSyms,refSyms,Basis)
    return CalcSiteSum(PairList,PairTypes,siteList,pairNumberDict,Basis)
end

function CalcSiteSum(PairList::AbstractVector{<:Rv},PairTypes,siteList::AbstractVector{<:Rv},pairNumberDict::PairNumbersDict,Basis::Basis_Struct) where {Rv<:Rvec}
    Npairs = length(PairList)
    
    Nsum = length(siteList)
    
    pairs = fill(sumElements(0,0,0,0),Nsum,Npairs) # generate empty matrix

    Threads.@threads for j in eachindex(PairList) #lhs of flow eqs
        Rj = PairList[j]
        Ri = Basis.refSites[PairTypes[j].xi]  # current reference site
        for (k,Rk) in enumerate(siteList)
            
            kj_pair = pairNumberDict[(Rk,Rj)]
            kj_pair == 0 && continue # skip if pair is not in PairList

            ki_pair = pairNumberDict[(Rk,Ri)]
            ki_pair == 0 && continue # skip if pair is not in PairList

            pairs[k,j]= sumElements(ki_pair,kj_pair,1,getSiteType(Rk,Basis))
        end

    end
    return pairs
end

"""
given a list of inequivalent pair pairs, returns splits along with their multiplicity
"""
function reduceSiteSum(siteSum)
    # pairs = CalcSiteSum(L,PairList,sitetypes)

    Nsum,Npairs = size(siteSum)
    reducedSum = fill(sumElements(0,0,0,0),Nsum+1,Npairs) #add +1 to Nsum to always have a (0,0,0,0) sumElements so we can always cut off last index!
    MaxNsum = 1
    Threads.@threads for j in 1:Npairs
        jSumList = sort( @view siteSum[:,j])
        uniqueElements = unique(jSumList,dims=1)
        len = length(uniqueElements)
        if len>MaxNsum
            # determine dimension of new matrix
            MaxNsum = len
        end
        for (k,pair) in enumerate(uniqueElements)
            if pair.ki != 0 && pair.kj != 0
                m = multiplicity_sorted(pair,jSumList)
                reducedSum[k,j] = sumElements(pair.ki,pair.kj,m,pair.xk)
            end
        end
    end
    return reducedSum[1:MaxNsum-1,:] #cut off unneded part, MaxNsum-1 because every part will definitely contain an empty split 
end


"""Returns a vector of indices with symmetry inequivalent pairs by applying all symmetries in Symmetrylist.
"""
function findSymmetryReduced(PairList::AbstractVector{RV},PairTypes::AbstractVector{sitePair},refSyms::AbstractVector{<:AbstractVector}) where {RV<:Rvec}
    lt(R1,R2) = convertToTuple(R1) < convertToTuple(R2)

    uniquePairs = unique(PairList)
    NUnique = maximum(max(x.xi,x.xj) for x in PairTypes)
    UniqueSites = [ empty(PairList) for x in 1:NUnique]
    siteBuffer = empty(PairList)
    for x in 1:NUnique
        for site in uniquePairs
            push!(siteBuffer,site)
            for T in refSyms[x]
                push!(siteBuffer,T(site))
            end
            sort!(siteBuffer,lt=lt)
            push!(UniqueSites[x], siteBuffer[begin]) 
            empty!(siteBuffer)
        end
    end
    keepInds = [
        MapToPair(x,R,PairList,PairTypes)
        for x in 1:NUnique for R in UniqueSites[x]
    ]
    unique!(keepInds)
    filter!(!=(0),keepInds)
    sort!(keepInds)
    return keepInds
end



"""Given a vector of vectors return all elements in a single vector"""
function flatten_nested(v::AbstractVector{<:AbstractVector})
    return [x for vi in v for x in vi]
end

"""Fallback for only one inequiv site for backward compatibility"""
findSymmetryReduced(PairList::AbstractVector{<:Rvec_3D},Symmetrylist::AbstractVector{<:Rvec}) = findSymmetryReduced(PairList,ones(length(PairList)),[Symmetrylist,])

"""Converts a pair of sites Rk, and Rj to a pair such that corresponds to one of the reference sites. For this, translation symmetry is used as well as a list of symmetries can transforms reference sites into each other.
"""
function pairToRefSite(Rk::RV,Rj::RV,Basis::Basis_Struct,nonRefSymmetries::A,refSymmetries::AbstractVector{<:A}) where {A<:AbstractVector{T} where T,RV<:Rvec}
    refSites_b = (R.b for R in Basis.refSites)
    # check if Rk  is already in first unit cell
    Rk.b ∈ refSites_b && return translatePairToUnitCell(Rk,Rj)

    Rk´ = Rk

    for Sym in nonRefSymmetries
        Rk´ = Sym(Rk)
        Rk´.b ∈ refSites_b && return translatePairToUnitCell(Rk´,Sym(Rj))
    end
    x = getSiteType(Rk,Basis)

    for (i,SymList) in enumerate(refSymmetries) #try symmetries that do not correspond this-reference site
        i == x && continue
        for Sym in SymList
            Rk´ = Sym(Rk)
            Rk´.b ∈ refSites_b && return translatePairToUnitCell(Rk´,Sym(Rj))
        end
    end

    for s in nonRefSymmetries
        println(s(Rk))
    end

    error("Could not find reference site for pair $(Rk),$(Rj)!")
end


"""Converts a pair of sites Rk, and Rj to a symmetry inequivalent pair by first applying the symmetries in nonRefSymmetries to map Rk to a reference site and then refSymmetries to map to a symmetry reduced sector.
"""
function pairToInequiv_vec(Rk::Rvec,Rj::Rvec,Basis::Basis_Struct,MapToPairDict::AbstractDict,nonRefSymmetries,refSymmetries::AbstractVector{<:AbstractVector})
    x = getSiteType(Rk,Basis)

    Rk,Rj = pairToRefSite(Rk,Rj,Basis,nonRefSymmetries,refSymmetries)

    x = getSiteType(Rk,Basis)
    sym = refSymmetries[x]
    for s in sym
        Rjprime = s(Rj)
        (x,Rjprime) in keys(MapToPairDict) && return (Rk,Rjprime)
    end
    return Rk,Rj
end

function generatePairToInequiv(InequivPairs::AbstractVector{<:Rvec},PairTypes::AbstractVector{sitePair},Basis::Basis_Struct,nonRefSymmetries,refSymmetries)
    MapToPairDict = getMapToPairDict(InequivPairs,PairTypes)
    @inline pairToInequiv(R1::Rvec,R2::Rvec) = pairToInequiv_vec(R1,R2,Basis,MapToPairDict,nonRefSymmetries,refSymmetries)
end

"""Returns generalized Geometry struct after specification of System size, symmetry function and basis.
Function pairToInequiv must map two sites to a pair containing a reference site
inCorrectSector maps a pair to true if it is symmetry inequivalent.
"""
function getLatticeGeometry(NLen,Name,pairToInequiv::Function,inCorrectSector::Function,Basis::Basis_Struct,dtype=Float64;method= generatePairSites,kwargs...)
    sortedpairs,sortedPairTypes = sortedPairList(NLen,Basis,method) # get sorted List of pairs for all reference sites
    
    AllSites = unique(sortedpairs) # get List of all sites that are included in the system

    inequivalentPairs,PairTypes = getinequivalentPairs(sortedpairs,sortedPairTypes, inCorrectSector, Basis) # remove symmetry equivalents

    getLatticeGeometry(NLen,Name,pairToInequiv::Function,AllSites,inequivalentPairs,PairTypes,Basis,dtype;kwargs...)

end

function getLatticeGeometry(NLen,Name,pairToInequiv::Union{<:Tuple,<:Function},AllSites,inequivalentPairs,PairTypes,Basis::Basis_Struct,dtype =Float64;test = false)
    splits = CalcSiteSum(inequivalentPairs,PairTypes,AllSites,pairToInequiv,Basis)
    siteSum = reduceSiteSum(splits)
    Npairs = size(siteSum,2)
    
    invpairs = [inversepair(Basis.refSites[x.xi],R,inequivalentPairs,PairTypes,Basis,pairToInequiv) for (x,R) in zip(PairTypes,inequivalentPairs)]
    
    couplings = zeros(dtype,Npairs)

    OnsitePairs = findOnsitePairs(inequivalentPairs,PairTypes,Basis.refSites)
    
    System = Geometry(Name = Name, NLen = NLen,couplings = couplings,PairList = inequivalentPairs,invpairs = invpairs,siteSum = siteSum,PairTypes = PairTypes,NUnique = Basis.NUnique,OnsitePairs = OnsitePairs)

    if test
        testGeometry(System)
    end
    return(System)
end

function getLatticeGeometry(NLen,Name,Basis::Basis_Struct,nonRefSymmetries,refSymmetries,dtype =Float64;test = false,method = generatePairSites)
    sortedpairs,sortedPairTypes = sortedPairList(NLen,Basis,method)

    AllSites = unique(sortedpairs)
    inds = findSymmetryReduced(sortedpairs,sortedPairTypes,refSymmetries)

    inequivalentPairs = sortedpairs[inds]
    PairTypes = sortedPairTypes[inds]

    pairNumberDict = pairNumbersDict(AllSites,inequivalentPairs,PairTypes,nonRefSymmetries,refSymmetries,Basis)

    splits = CalcSiteSum(inequivalentPairs,PairTypes,AllSites,pairNumberDict,Basis)
    siteSum = reduceSiteSum(splits)
    Npairs = size(siteSum,2)
    
    invpairs = [inversepair(Basis.refSites[x.xi],R,pairNumberDict) for (x,R) in zip(PairTypes,inequivalentPairs)]
    
    couplings = zeros(dtype,Npairs)

    OnsitePairs = findOnsitePairs(inequivalentPairs,PairTypes,Basis.refSites)
    
    System = Geometry(;Name, NLen,couplings,invpairs,siteSum,PairTypes,OnsitePairs,PairList = inequivalentPairs,NUnique = Basis.NUnique)

    if test
        testGeometry(System)
    end
    return(System)
end

"""return all neccessary information about reduced lattice
"""
function generateReducedLattice(NLen,Basis::Basis_Struct,nonRefSymmetries,refSymmetries;method = generatePairSites)
    sortedpairs,sortedPairTypes = sortedPairList(NLen,Basis,method)
    AllSites = unique(sortedpairs)
    UnitCell = getUnitCell(Basis)
    inds = findSymmetryReduced(sortedpairs,sortedPairTypes,refSymmetries)
    
    PairList = sortedpairs[inds]
    PairTypes = sortedPairTypes[inds]
    pairNumberDict = pairNumbersDict(AllSites,PairList,PairTypes,nonRefSymmetries,refSymmetries,Basis)
    return (;pairNumberDict,AllSites,PairList,PairTypes,Basis)
end

function convertDictToArrays(pairNumbersDict)
    ks = collect(keys(pairNumbersDict))
    vs = collect(values(pairNumbersDict))
    R1vec = StructArray([x[1] for x in ks])
    R2vec = StructArray([x[2] for x in ks])
    return R1vec,R2vec,vs
end

function findOnsitePairs(PairList,PairTypes,refSites)
    OSP = Int[]
    for (i,(R,T)) in enumerate(zip(PairList,PairTypes))
        if R in refSites && T.xi == T.xj
            push!(OSP,i)
        end
    end
    return OSP
end

function getFRGComplexity(System::Geometry)
    s = 0
    for i in eachindex(System.Npairs)
        s += System.Npairs[i] *System.Nsum[i]
    end
    return s
end

getFRGComplexity(S1::Geometry,S2::Geometry) = getFRGComplexity(S1)/getFRGComplexity(S2)
