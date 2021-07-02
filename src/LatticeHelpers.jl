"""
Functions and data types that can be used to construct Lattices.
"""

export  Basis_Struct_2D, Basis_Struct_3D,Basis_Struct, Rvec_2D, Rvec_3D, Rvec, getLatticeVec, norm, translateToOrigin, translation, getCartesian, dist, generatePairSites,generateLUnitCells, MapToPair, setCoupling!,setNeighborCouplings!, getSiteType, testGeometry, getLatticeGeometry


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
    for n1 in -L:L, n2 in -L:L, b in 1:Basis.NCell
        push!(PairList,Rvec(n1,n2,b))
    end
    return PairList
end

function generateLUnitCells(L,Basis::Basis_Struct_3D,refSite = Basis.refSites[1])
    # return
    PairList = Vector{Rvec_3D}(undef,0)
    for n1 in -L:L, n2 in -L:L, n3 in -L:L, b in 1:Basis.NCell
        push!(PairList,Rvec(n1,n2,n3,b))
    end
    return PairList
end

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
function MapToPair(x::Integer,R::Rvec,PairList,PairTypes)
    for (i,b) in enumerate(PairList)
        if b === R && PairTypes[i].xi === x 
            return(i)
        end
    end
    # error("error: pair not found")
    return 0 # no pair found
    # println("cannot map: ",R," to inequivalent pair. Not found in PairList for site ",x)
end

"""Sets coupling between ref site x and R to val. To map from x and R to the corresponding pair, we need PairList and siteTypes"""
function setCoupling!(couplings,x,R,val::Number,PairList,PairTypes)
    idx = MapToPair(x,R,PairList,PairTypes)
    if idx <=0
        println(PairList)
        println(PairTypes)
        error("problem setting coupling for (x,R,val) = ",(x,R,val))
    end
    couplings[idx] = val
end

"""Sets couplings according to given list with the first element being the nearest neigbor coupling. Currently only works for lattices with equivalent sites and sorted PairList!"""
function setNeighborCouplings!(couplings,Jparams::AbstractVector,PairList,Basis)
    Npairs = length(PairList)
    neigborcoup = zeros(Npairs)    
    println(Npairs)
    neigborcoup[2:1+length(Jparams)] = Jparams
    
    index = 1
    lastr = 0.
    norm_B(R) = norm(R,Basis)
    norms = norm_B.(PairList)
    for (i,r) in enumerate(norms)
        if abs(r - lastr)>1E-14
            index +=1
            lastr = r
        end
        couplings[i] = neigborcoup[index]
    end
    return couplings
end

function multiplicity(spl::sumElements,list)
    N = length(list)
    mult = 0
    for i in 1:N
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
    for (j,Rj) in enumerate(PairList) #lhs of flow eqs
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
    for j in 1:Npairs
        jSumList = @view siteSum[:,j]
        uniqueElements = unique(jSumList,dims=1)
        sort!(uniqueElements)
        len = length(uniqueElements)
        if len>MaxNsum
            # determine dimension of new matrix
            MaxNsum = len
        end
        for (k,pair) in enumerate(uniqueElements)
            m = multiplicity(pair,jSumList)
            if pair.ki != 0 && pair.kj != 0
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

function testGeometry(Geometry)

    @unpack siteSum,Npairs,invpairs,OnsitePairs,NUnique,PairList,PairTypes = Geometry
    @testset "SiteSum" begin 
        for j in OnsitePairs
            for spl in siteSum[:,j]
                @test spl.ki == spl.kj
            end
        end
    end
    @testset "PairList" begin
        for (j,pair) in enumerate(PairList)
            @unpack xi = PairTypes[j]
            @test MapToPair(xi,pair,PairList,PairTypes) == j # j must correspond to pair in PairList
        end
    end

    @testset "invpairs" begin
        for i in invpairs
            @test i > 0 && i <= Npairs
        end
    end
    
    @testset "OnSitePairs" begin
        @test length(OnsitePairs) == NUnique   
        @test invpairs[OnsitePairs] == OnsitePairs # Test that OnsitePairs are their own inverse  
    end

    @testset "Allowed values" begin
        for spl in siteSum
            if spl.ki != 0 || spl.kj != 0 || spl.xk != 0 || spl.xk != 0
                @test spl.ki > 0
                @test spl.kj > 0
                @test spl.m > 0
                @test spl.xk > 0
            end
            @test spl.xk <= NUnique

            if spl.ki == 0 || spl.kj == 0 || spl.xk == 0
                @test spl.m == 0
            end
        end
    end

    symmetryWarning = "It is possible that some inequivalent pairs are missing/ double accounted for!"

    for j in 1:Npairs
        if invpairs[invpairs[j]] != j # Test that inversing twice gives the same index 
            println(symmetryWarning)
            println(invpairs[invpairs] ," != ", invpairs)
            break
        end
    end

    for j in OnsitePairs
        j1cut = siteSum[:,j] #get site sum for onsite pairs ki = kj
        inds =[s.ki for s in j1cut] #select ki elements
        deleteat!(inds, findall(x->x==0,inds)) # remove 0
        OnsiteSum = sort(invpairs[inds])  # find inverse pairs ik
        for (i,x) in enumerate(OnsiteSum)
            if x != j-1 + i # If symmetries are fully incorporated, we expect each pair to appear once
                println(symmetryWarning)
                println(OnsiteSum)
                return
            end
        end

    end
end
