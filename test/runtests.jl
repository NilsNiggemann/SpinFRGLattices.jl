using SpinFRGLattices,Test

@test norm([1,1,1]) ≈ sqrt(3)

@testset verbose = true "Testing Lattices" begin
    @testset "Testing C2SquareKagome" begin
        exemplary_PairList_i_13_23 = 
        [    
            Rvec(-1, 0, 22),
            Rvec(0, 0, 5),
            Rvec(0, 0, 15),
            Rvec(0, -1, 14),
            Rvec(0, 0, 4),
            Rvec(0, -1, 13),
            Rvec(-1, 0, 7),
            Rvec(0, 0, 13),
            Rvec(0, 0, 7),
            Rvec(-1, -1, 20),
            Rvec(-1, 0, 21),
        ]
        SK = C2SquareKagome.getC2DimerSquareKagome(6)
        @testset "Testing PairList indices 13:23" begin
            @test SK.PairList[13:23] == exemplary_PairList_i_13_23
        end

    end

    @testset "Testing LargeSquareKagome" begin
        testGeometry(LargeSquareKagome.getNoMirrDimerSquareKagome(6))
    end

    @testset "Testing Pyrochlore" begin
        testGeometry(Pyrochlore.getPyrochlore(6))
    end

    @testset "Testing SimpleCubic" begin
        testGeometry(SimpleCubic.getCubic(6))
    end

    @testset "Testing SquareKagome" begin
        exemplary_PairList_i_15_28 =
        [
            Rvec(1, 0, 6),
            Rvec(1, 1, 2),
            Rvec(1, 1, 6),
            Rvec(0, 0, 2),
            Rvec(0, 0, 6),
            Rvec(0, 0, 5),
            Rvec(0, 1, 3),
            Rvec(0, 1, 4),
            Rvec(-1, 1, 4),
            Rvec(0, 0, 4),
            Rvec(1, 0, 1),
            Rvec(0, 1, 1),
            Rvec(1, 0, 2),
            Rvec(0, 1, 2),
        ]
        SK = SquareKagome.getMirrorSquareKagome(6)
        @testset "Testing PairList indices 15:28" begin
            @test SK.PairList[15:28] == exemplary_PairList_i_15_28
        end

        @testset "Checksums" begin
            @test sum(SK.siteSum.ki) == 2052810

            @test sum(SK.siteSum.kj) == 2052810

            @test sum(SK.siteSum.xk) == 39639

            @test sum(SK.siteSum.m) == 24980

        end
        testGeometry(SK)
        testGeometry(SK)
    end

    @testset "Testing SquareKagome" begin
        testGeometry(SquareKagome.getSquareKagome(6))
    end

    @testset "Testing SquareLattice" begin
        testGeometry(SquareLattice.getSquareLattice(6))
    end

    @testset "Testing Octochlore" begin
        testGeometry(Octochlore.getOctochlore(6))
        testGeometry(Octochlore.getOctochloreGamma(6))
    end

    @testset "Testing Honeycomb" begin
        testGeometry(Honeycomb.getHoneycomb(6,[1.,2.,3.,4.]))
    end

    @testset "Testing Honeycomb" begin
        testGeometry(Honeycomb.getHoneycombGamma(6,alpha = 1., gamma = 0.3))
    end

    @testset "Testing TriangularLattice" begin
        testGeometry(TriangularLattice.getTriangularLattice(6))
    end
end

function testOnsitePairAdapt()
    @testset "adapt Onsite Pairs" begin
        @test SpinFRGLattices.adaptOnsitePairs([1]) == [1]
        @test SpinFRGLattices.adaptOnsitePairs([1,4]) == [1,5]
        @test SpinFRGLattices.adaptOnsitePairs([1,8]) == [1,9]
    end
end
testOnsitePairAdapt()

S1 = LargeSquareKagome.getDimerSquareKagome(3)
testPairListAdaptation(S1,LargeSquareKagome.Basis.NCell)