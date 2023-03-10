export testPairListSym

function testGeometry(Geo::Geometry)

    (;siteSum,Npairs,invpairs,couplings,OnsitePairs,NUnique,PairList,PairTypes) = Geo
    @testset "SiteSum" begin 
        for j in OnsitePairs
            for spl in siteSum[:,j]
                @test spl.ki == spl.kj
            end
        end
    end
    @testset "PairList" begin
        for (j,pair) in enumerate(PairList)
            xi = PairTypes[j].xi
            @test MapToPair(xi,pair,PairList,PairTypes) == j # j must correspond to pair in PairList
        end
    end

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
        inds =[s.ki for s in j1cut] #select ki elements
        deleteat!(inds, findall(x->x==0,inds)) # remove 0
        OnsiteSum = sort(invpairs[inds])  # find inverse pairs ik
        for (i,x) in enumerate(OnsiteSum)
            if x != j-1 + i # If symmetries are fully incorporated, we expect each pair to appear once
                @warn symmetryWarning
                println(OnsiteSum)
                return
            end
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