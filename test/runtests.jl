using SpinFRGLattices,Test
@test norm([1,1,1]) ≈ sqrt(3)

@testset verbose = true "Testing Lattices" begin

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

    @testset "Testing Honeycomb" begin
        testGeometry(Honeycomb.getHoneycomb(6,[1.,2.,3.,4.]))
    end

    @testset "Testing TriangularLattice" begin
        testGeometry(TriangularLattice.getTriangularLattice(6))
    end
end