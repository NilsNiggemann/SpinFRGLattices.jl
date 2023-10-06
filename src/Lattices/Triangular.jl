##
module TriangularLattice
    using StaticArrays,Test
    using ..SpinFRGLattices
    export getTriangularLattice
        
    """Constructer for a square lattice basis"""
    function TriangularBasis()
        a1 = SA[1,0] #It's best to use StaticArrays since the dimension is small and known at compile-time
        a2 = 1/2*SA[1,√3]
        
        b1 = SA[0,0]
        return Basis_Struct_2D(a1=a1,a2=a2,b=[b1], NNdist = 1.)
    end
    const Basis = TriangularBasis()
    getCartesian(R) = SpinFRGLattices.getCartesian(R,Basis) 
    dist(R1,R2) = SpinFRGLattices.dist(R1,R2,Basis) 
    getRvec(R) = SpinFRGLattices.getRvec(R,Basis)

    """performs C6 rotation on R"""
    function C6(r::AbstractVector)
        SA[1/2 -sqrt(3)/2; sqrt(3)/2 1/2] * r
    end

    """performs C6 rotation on R"""
    function C6(R::Rvec)
        return R |> getCartesian |> C6 |> getRvec
    end

    fMirr(x) = x / √3

    Mirror(x) = MirrorLine(1/sqrt(3),x,Basis)

    function inC6Sector(R::Rvec)
        (;n1,n2) = R
        n1 >= 0 && n2 >= 0
    end
    """
    Symmetry: may restrict to section with y>=0 and below mirror Line
    """
    function inCorrectSubsector(R_ref,R::Rvec)
        return(inC6Sector(R) && R.n2 <=R.n1)
    end

    inCorrectSubsector(R::Rvec) = inCorrectSubsector(Rvec(0,0,1),R)

    """gives coordinates of R' after using symmetry transformations on R
    R is arbitrary lattice vector"""
    function mapToSubsector(R::Rvec)
        i = 1
        Rold = R
        while !inC6Sector(R)
            i+=1
            i >6 && error("could not rotate to C6 sector: ", Rold)
            R = C6(R)
        end
        if !inCorrectSubsector(R)
            R = Mirror(R)
        end
        inCorrectSubsector(R) || error("Error: Vector not correctly mapped: ",Rold)
        return R
    end

    """given a pair of sites Rk and Rj from a vertex, return a corresponding inequivalent pair R_new by using lattice symmetries. This means R_new is in an equivalent position to site 0 as Rk is to Rj. (Rk,Rj) -> (0,R_new) """
    function pairToInequiv(Rk::Rvec,Rj::Rvec)
        Rj = translateToOrigin(Rj,Rk) #global translation to map Rk to first unit cell
        Rk = Rvec(0,0,1)
        return Rk, mapToSubsector(Rj) #give the vector in correct subsector
    end

    """
    returns Geometry struct with all relevant information about the Triangular lattice. We can either consider N nearest neighbor pairs or all unit cells within length L of origin.
    """
    function getTriangularLattice(NLen,J = [1.,0.5];kwargs...)
        Name = string("TriangularLattice_NLen=",NLen)
        System =  getLatticeGeometry(NLen,Name,pairToInequiv,inCorrectSubsector,Basis;kwargs...)

        (;PairList,couplings) = System
        setNeighborCouplings!(System,J,Basis)

        return(System)
    end
end