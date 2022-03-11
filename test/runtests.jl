using SpinFRGLattices,Test
@test norm([1,1,1]) ≈ sqrt(3)

println("Testing C2SquareKagome")
C2SquareKagome.getC2DimerSquareKagome(6,test = true)

println("Testing LargeSquareKagome")
LargeSquareKagome.getNoMirrDimerSquareKagome(6,test = true)

println("Testing Pyrochlore")
Pyrochlore.getPyrochlore(6,test=true)

println("Testing SimpleCubic")
SimpleCubic.getCubic(6,test=true)

println("Testing SquareKagome")
SquareKagome.getMirrorSquareKagome(6,test=true)

println("Testing SquareKagome")
SquareKagome.getSquareKagome(6,test=true)

println("Testing SquareLattice")
SquareLattice.getSquareLattice(6,test=true)

println("Testing Octochlore")
Octochlore.getOctochlore(6,test=true)

println("Testing Honeycomb")
Honeycomb.getHoneycomb(6,[1.,2.,3.,4.],test=true)

println("Testing Honeycomb")
Honeycomb.getHoneycombGamma(6,alpha = 1., gamma = 0.3,test=true)

println("Testing TriangularLattice")
TriangularLattice.getTriangularLattice(6,test=true)
