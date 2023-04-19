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

"""Sets couplings according to given list with the first element being the nearest neigbor coupling. Currently only works for lattices with equivalent sites and sorted PairList!"""
function setNeighborCouplings!(couplings,Jparams::AbstractVector,PairList,PairTypes,Basis)
    couplings .= 0.

    function neigborcoup(i)
        i in eachindex(Jparams) && return Jparams[i]
        return 0.
    end
    norms = similar(couplings)

    for i in eachindex(norms,PairList,PairTypes)
        xi = PairTypes[i].xi
        R_Ref = Basis.refSites[xi]
        distance(x) = dist(R_Ref,x,Basis)
        norms[i] = distance(PairList[i])
    end
    for xi in eachindex(Basis.refSites)
        # lt(x,y) = distance(x) <distance(y)
        # SortedIndices = sortperm(PairList,lt = lt)
        index = 0
        lastr = 0.
        
        pairInds = findall(x-> x.xi == xi,PairTypes)

        for i in pairInds
            r = norms[i]
            if abs(r - lastr)>1E-14
                index +=1
                lastr = r
            end
            if couplings[i] == 0. 
                couplings[i] = neigborcoup(index)
            else
                @assert couplings[i] == neigborcoup(index) "couplings are not compatible between inequivalent sites"
            end

        end
    
    end
    return couplings
end

function setNeighborCouplings!(couplings,Jparams::AbstractVector,PairList,Basis)
    couplings .= 0.
    neigborcoup = copy(couplings)

    neigborcoup[2:1+length(Jparams)] = Jparams

    R_Ref = only(Basis.refSites)
    distance(x) = dist(R_Ref,x,Basis)
    # lt(x,y) = distance(x) <distance(y)
    # SortedIndices = sortperm(PairList,lt = lt)
    norms = distance.(PairList)
    index = 1
    lastr = 0.
    for (i,r) in enumerate(norms)
        if abs(r - lastr)>1E-14
            index +=1
            lastr = r
        end
        if couplings[i] == 0. 
            couplings[i] = neigborcoup[index]
        else
            @assert couplings[i] == neigborcoup[index] "couplings are not compatible between inequivalent sites"
        end

    end
    return couplings
end

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

function generatePairNumberDict(siteList,PairList,PairTypes,nonRefSyms,refSyms,Basis)
    MapToPairDict = getMapToPairDict(PairList,PairTypes)
    pairNumberDict = Dict((R1,R2) => pairNumber(R1,R2,MapToPairDict,nonRefSyms,refSyms,Basis) for R1 in siteList for R2 in siteList)
    filter!(x-> x.second != 0,pairNumberDict) # remove all pairs that are not in PairList
    return pairNumberDict
end

"""This function is for backwards compatibility. It is recommended to use the function with the same name that uses the mapToPairDict instead of the pairToInequiv function"""
function generatePairNumberDict(siteList,PairList,PairTypes,pairToInequiv::Function,Basis)
    pairNumber_(R1,R2) = pairNumber(R1,R2,PairList,PairTypes,Basis,pairToInequiv)
    pairNumberDict = Dict([ (R1,R2) => pairNumber_(R1,R2) for R1 in siteList for R2 in siteList]) # store the pair number of all possible pairs in a dictionary
    filter!(x-> x.second != 0,pairNumberDict) # remove all pairs that are not in PairList

    return pairNumberDict
end

"""Taking an inequivalent pair (x,R) with x <= NUnique returns corresponding inversed pair (R,x)"""
function inversepair(R_ref::Rvec,R::Rvec,PairList,PairTypes,Basis,pairToInequiv::Function)
    return pairNumber(R,R_ref,PairList,PairTypes,Basis,pairToInequiv)
end

inversepair(R_ref::Rvec,R::Rvec,pairNumbersDict::AbstractDict) = pairNumbersDict[(R,R_ref)]
"""
Performs site summation for each Rj (lhs) of flow equations and returns all pairs that appear in it corresponding to pairs of vertices on rhs.
It is assumed that our reference sites are the first in the basis!
"""
function CalcSiteSum(PairList::AbstractVector{<:Rvec},PairTypes::AbstractVector{sitePair},siteList::AbstractVector{<:Rvec},pairToInequiv::Function,Basis::Basis_Struct)
    pairNumber_(R1,R2) = pairNumber(R1,R2,PairList,PairTypes,Basis,pairToInequiv)
    pairNumberDict = Dict([ (R1,R2) => pairNumber_(R1,R2) for R1 in siteList for R2 in PairList]) # store the pair number of all possible pairs in a dictionary
    filter!(x-> x.second != 0,pairNumberDict) # remove all pairs that are not in PairList
    return CalcSiteSum(PairList,PairTypes,siteList,pairNumberDict,Basis)
end

function CalcSiteSum(PairList::AbstractVector{<:Rvec},PairTypes::AbstractVector{sitePair},siteList::AbstractVector{<:Rvec},syms::Tuple,Basis::Basis_Struct)
    (nonRefSyms,refSyms) = syms
    pairNumberDict = generatePairNumberDict(siteList,PairList,PairTypes,nonRefSyms,refSyms,Basis)
    return CalcSiteSum(PairList,PairTypes,siteList,pairNumberDict,Basis)
end

function CalcSiteSum(PairList::AbstractVector{<:Rv},PairTypes,siteList::AbstractVector{<:Rv},pairNumberDict::AbstractDict,Basis::Basis_Struct) where {Rv<:Rvec}
    Npairs = length(PairList)
    
    Nsum = length(siteList)
    
    pairs = fill(sumElements(0,0,0,0),Nsum,Npairs) # generate empty matrix

    Threads.@threads for j in eachindex(PairList) #lhs of flow eqs
        Rj = PairList[j]
        Ri = Basis.refSites[PairTypes[j].xi]  # current reference site
        for (k,Rk) in enumerate(siteList)
            if (Rk,Ri) ∉ keys(pairNumberDict) || (Rk,Rj) ∉ keys(pairNumberDict)
                continue
            end

            ki_pair = pairNumberDict[(Rk,Ri)] # maps vertex V_ki to appropriate pair in PairList
            
            kj_pair = pairNumberDict[(Rk,Rj)]

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
    _tuple(R::RV) = Tuple(getproperty(R,x) for x in fieldnames(typeof(R)))
    lt(R1,R2) = _tuple(R1) < _tuple(R2)
    # orderInds = sortperm(PairList,lt = (R1,R2) -> _tuple(R1) < _tuple(R2))
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

"""Converts a pair of sites Rk, and Rj to a pair such that Rk lies in the first unit cells. For this, translation symmetry is used as well as a list of symmetries can transforms reference sites into each other.
"""
function pairToRefSite(Rk::Rvec_3D,Rj::Rvec_3D,Basis::Basis_Struct,nonRefSymmetries,refSymmetries)
    # refSites_b = getproperty.(Basis.refSites,:b)
    refSites_b = (R.b for R in Basis.refSites)
    Rk´ = Rk
    Rj´ = Rj
    x = getSiteType(Rk,Basis)

    for SymList in (nonRefSymmetries,refSymmetries...)

        SymList === refSymmetries[x] && continue

        for Sym in SymList
            Rk´.b ∈ refSites_b && return Rvec(0,0,0,Rk´.b),translateToOrigin(Rj´,Rk´)
            Rk´ = Sym(Rk)
            Rj´ = Sym(Rj)
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

    pairNumberDict = generatePairNumberDict(AllSites,inequivalentPairs,PairTypes,nonRefSymmetries,refSymmetries,Basis)

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
