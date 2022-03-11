using SpinFRGLattices,Test
@test norm([1,1,1]) ≈ sqrt(3)

C2SquareKagome.getC2DimerSquareKagome(6,test = true)
println("Testing C2SquareKagome")

LargeSquareKagome.getNoMirrDimerSquareKagome(6,test = true)
println("Testing LargeSquareKagome")

Pyrochlore.getPyrochlore(6,test=true)
println("Testing Pyrochlore")

SimpleCubic.getCubic(6,test=true)
println("Testing SimpleCubic")

SquareKagome.getMirrorSquareKagome(6,test=true)
println("Testing SquareKagome")

SquareKagome.getSquareKagome(6,test=true)
println("Testing SquareKagome")

SquareLattice.getSquareLattice(6,test=true)
println("Testing SquareLattice")

Octochlore.getOctochlore(6,test=true)
println("Testing Octochlore")

Honeycomb.getHoneycomb(6,[1.,2.,3.,4.],test=true)
println("Testing Honeycomb")

Honeycomb.getHoneycombGamma(6,alpha = 1., gamma = 0.3,test=true)
println("Testing Honeycomb")

TriangularLattice.getTriangularLattice(6,test=true)
println("Testing TriangularLattice")
