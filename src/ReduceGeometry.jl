export MapToPair, setCoupling!,setNeighborCouplings!, getSiteType, testGeometry, getLatticeGeometry,getFRGComplexity
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
function pairNumber(R1,R2,PairList,PairTypes, Basis,pairToInequiv::Function)
    R1,R2 = pairToInequiv(R1,R2)
    return MapToPair(getSiteType(R1,Basis),R2,PairList,PairTypes)
end 

"""Taking an inequivalent pair (x,R) with x <= NUnique returns corresponding inversed pair (R,x)"""
function inversepair(R_ref::Rvec,R::Rvec,PairList,PairTypes,Basis,pairToInequiv::Function)
    return pairNumber(R,R_ref,PairList,PairTypes,Basis,pairToInequiv)
end

"""
Performs site summation for each Rj (lhs) of flow equations and returns all pairs that appear in it corresponding to pairs of vertices on rhs.
It is assumed that our reference sites are the first in the basis!
"""
function CalcSiteSum(PairList,siteList,PairTypes,pairToInequiv::Function,Basis)
    Npairs = length(PairList)
    
    Nsum = length(siteList)
    
    pairs = fill!(Matrix{sumElements}(undef,Nsum,Npairs),sumElements(0,0,0,0)) # generate empty matrix

    pairNumber_(R1,R2) = pairNumber(R1,R2,PairList,PairTypes,Basis,pairToInequiv)
    Threads.@threads for j in eachindex(PairList) #lhs of flow eqs
        Rj = PairList[j]
        Ri = Basis.refSites[PairTypes[j].xi]  # current reference site
        for (k,Rk) in enumerate(siteList)
            ki_pair = pairNumber_(Rk,Ri) # maps vertex V_ki to appropriate pair in PairList
            
            kj_pair = pairNumber_(Rk,Rj)
            
            if ki_pair != 0 && kj_pair != 0
                pairs[k,j]= sumElements(ki_pair,kj_pair,1,getSiteType(Rk,Basis))
            else
                pairs[k,j]= sumElements(0,0,0,0)
            end
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
    reducedSum = fill!(Matrix{sumElements}(undef,Nsum+1,Npairs),sumElements(0,0,0,0)) #add +1 to Nsum to always have a (0,0,0,0) sumElements so we can always cut off last index!
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
"""Returns generalized Geometry struct after specification of System size, symmetry function and basis.
Function pairToInequiv must map two sites to a pair containing a reference site
inCorrectSector maps a pair to true if it is symmetry inequivalent.
"""
function getLatticeGeometry(NLen,Name,pairToInequiv::Function,inCorrectSector::Function,Basis;test = false,method= generatePairSites)
    sortedpairs,sortedPairTypes = sortedPairList(NLen,Basis,method) # get sorted List of pairs for all reference sites
    
    AllSites = unique(sortedpairs) # get List of all sites that are included in the system
    Ntot = length(AllSites)

    inequivalentPairs,PairTypes = getinequivalentPairs(sortedpairs,sortedPairTypes, inCorrectSector, Basis) # remove symmetry equivalents
    println("Total Number of sites: ",Ntot, "\t Num pairs: ", length(inequivalentPairs) )

    splits = CalcSiteSum(inequivalentPairs,AllSites,PairTypes,pairToInequiv,Basis)
    siteSum = reduceSiteSum(splits)
    Npairs = size(siteSum,2)
    
    OnsitePairs = unique(i -> PairTypes[i].xi, 1:Npairs)
    
    invpairs = [inversepair(Basis.refSites[x.xi],R,inequivalentPairs,PairTypes,Basis,pairToInequiv) for (x,R) in zip(PairTypes,inequivalentPairs)]
    
    couplings = zeros(Npairs)

    System = Geometry(Name = Name, NLen = NLen,couplings = couplings,PairList = inequivalentPairs,invpairs = invpairs,siteSum = siteSum,PairTypes = PairTypes,NUnique = Basis.NUnique,OnsitePairs = OnsitePairs)

    if test
        testGeometry(System)
    end
    return(System)
end

function getFRGComplexity(System::Geometry)
    s = 0
    for i in eachindex(System.Npairs)
        s += System.Npairs[i] *System.Nsum[i]
    end
    return s
end

getFRGComplexity(S1::Geometry,S2::Geometry) = getFRGComplexity(S1)/getFRGComplexity(S2)