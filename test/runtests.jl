using SpinFRGLattices,Test
@test norm([1,1,1]) ≈ sqrt(3)

@testset verbose = true "Testing Lattices" begin
    @testset "Testing C2SquareKagome" begin
        testGeometry(C2SquareKagome.getC2DimerSquareKagome(6))
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
        testGeometry(SquareKagome.getMirrorSquareKagome(6))
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
    @testset "Testing TriangularLattice" begin
        testGeometry(TriangularLattice.getTriangularLattice(6))
    end
end
