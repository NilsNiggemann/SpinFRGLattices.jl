function testGeometry(Geo::Geometry)

    (;siteSum,Npairs,invpairs,couplings,OnsitePairs,NUnique,PairList,PairTypes) = Geo

    testSiteSum(siteSum,OnsitePairs)
    testAllowedValues(siteSum,NUnique)
    testOnsiteSum(siteSum,OnsitePairs)

    testPairList(PairList,PairTypes)

    @testset "invpairs" begin
        for i in invpairs
            @test i > 0 && i <= Npairs
        end
    end

    @testset "couplings" begin
        for (i,c) in zip(invpairs,couplings)
            @test isapprox(c, couplings[i],atol =eps(Float64))
        end
    end
            
    @testset "OnSitePairs" begin
        
        @test OnsitePairs == unique(i -> PairTypes[i].xi, 1:Npairs) #Convention: Onsite pairs are always the first pair for each ref site
        @test PairTypes[OnsitePairs] == [sitePair(i,i) for i in 1:NUnique]
        @test length(OnsitePairs) == NUnique   
        @test invpairs[OnsitePairs] == OnsitePairs # Test that OnsitePairs are their own inverse  
    end


    symmetryWarning = "It is possible that some inequivalent pairs are missing or doubly accounted for!"

    for j in 1:Npairs
        if invpairs[invpairs[j]] != j # Test that inversing twice gives the same index 
            @warn symmetryWarning
            println(invpairs[invpairs] ," != ", invpairs)
            break
        end
    end
    for j in OnsitePairs
        j1cut = siteSum[:,j] #get site sum for onsite pairs ki = kj
        inds = j1cut.ki #select ki elements
        deleteat!(inds, findall(x->x==0,inds)) # remove 0
        OnsiteSum = sort(invpairs[inds])  # find inverse pairs -> We expect these terms to be of the form (Vi1,Vi2,Vi3,Vi4,...) with each inequiv pair appearing once
        @test OnsiteSum == collect(OnsiteSum[1]:OnsiteSum[end])


    end

end
"""Xii = sum_k V_ki V_ki P_kk. --> In site sums corresponding to onsite pairs, the only onsite pair that can appear is (i,i). If there are more than one inequivalent sites, the other pairs can not appear in this summation!
"""
function testOnsiteSum(siteSum,OnsitePairs)
    @testset "OnSiteSum" begin 
        for j in OnsitePairs
            for ki in siteSum[:,j].ki
                if ki in OnsitePairs
                    @test ki == j
                end
            end
        end
    end
end

function testSiteSum(siteSum,OnsitePairs)

    @testset "SiteSum OnsitePairs" begin 
        for j in OnsitePairs
            for spl in siteSum[:,j]
                @test spl.ki == spl.kj
            end
        end
    end
end

function testAllowedValues(siteSum,NUnique)

    @testset "Allowed values" begin
        for spl in siteSum
            if spl.ki != 0 || spl.kj != 0 || spl.xk != 0 || spl.xk != 0
                @test spl.ki > 0
                @test spl.kj > 0
                @test spl.m > 0
                @test spl.xk > 0
            end
            @test spl.xk <= NUnique

            if 0 in (spl.ki, spl.kj, spl.xk)
                @test spl.m == 0
            end
        end
    end
    
end


function testNsum(Nsum,siteSum)
    @testset "number of sum terms" begin
        @test calcMaxpairs(siteSum) == Nsum
    end
end
testNsum(G::Geometry) = testNsum(G.Nsum,G.siteSum)

function testPairList(PairList,PairTypes)
    @testset "PairList" begin
        for (j,pair) in enumerate(PairList)
            xi = PairTypes[j].xi
            @test MapToPair(xi,pair,PairList,PairTypes) == j # j must correspond to pair in PairList
        end
    end
end

function testPairListSym(PairList::AbstractVector,Sym::Function)
    nowarnings = true
    for R in PairList
        Rn = Sym(R)
        (Rn in PairList && Rn != R) && (nowarnings = false; @warn "$R is symmetry reducible with $Rn")
    end
    return nowarnings
end

function testPairListSym(PairList::AbstractVector,Syms)
    nowarnings = true
    for Sym in Syms
        if !testPairListSym(PairList,Sym)
            nowarnings = false
        end
    end
    return nowarnings
end

function checkSymmetries(PairList::Vector{RT},Sym::Function) where RT <: Rvec
    for R in PairList
        Rn = Sym(R)
    end
    return true
end

function checkSymmetries(PairList::Vector{RT},Syms) where RT <: Rvec
    for Sym in Syms
        checkSymmetries(PairList,Sym)
    end
    return true
end


function testPairListAdaptation(G::Geometry,NCell::Int)
    K = adaptPairs(G,NCell,0.)
    @testset "test PairList Adaptation" begin
        @test K.PairList[K.OnsitePairs] == G.PairList[G.OnsitePairs]
        for i in eachindex(K.OnsitePairs,G.OnsitePairs)
            @test K.PairList[K.OnsitePairs[i]+1] == createSatellite(G.PairList[G.OnsitePairs[i]],NCell+i)
        end
    end
    @testset "test PairTypes" begin
        @test K.PairTypes[K.OnsitePairs] == G.PairTypes[G.OnsitePairs]
    end
end