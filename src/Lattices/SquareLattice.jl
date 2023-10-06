##
module SquareLattice
    using StaticArrays,Test
    using ..SpinFRGLattices
    export getSquareLattice
    
    """Constructer for a square lattice basis"""
    function SquareBasis()
        a1 = SA[1,0] #It's best to use StaticArrays since the dimension is small and known at compile-time
        a2 = SA[0,1]
        
        b1 = SA[0,0]
        return Basis_Struct_2D(a1=a1,a2=a2,b=[b1], NNdist = 1.)
    end
    const Basis = SquareBasis()
    getCartesian(R) = SpinFRGLattices.getCartesian(R,Basis) 
    dist(R1,R2) = SpinFRGLattices.dist(R1,R2,Basis) 
    getRvec(R) = SpinFRGLattices.getRvec(R,Basis) 
    ##

    """performs C4 rotation on R"""
    function C4(R::Rvec)
        return Rvec(-R.n2,R.n1,R.b)
    end

    """mirrors R on x axis"""
    function xmirror(R::Rvec)
        return Rvec(R.n1,-R.n2,R.b)
    end

    """mirrors R on y axis"""
    function ymirror(R::Rvec)
        return Rvec(-R.n1,R.n2,R.b)
    end

    """mirrors R on y=x axis"""
    function xymirror(R::Rvec)
        return Rvec(R.n2,R.n1,R.b)
    end

    """
    Symmetry: may restrict to section with y<=x
    """
    function inCorrectSubsector(R_ref,R::Rvec)
        (;n1,n2) = R
        return(n1>= 0 && n2>=0 && n2 <=n1)
    end

    inCorrectSubsector(R::Rvec) = inCorrectSubsector(Rvec(0,0,1),R)


    """gives coordinates of R' after using symmetry transformations on R
    R is arbitrary lattice vector"""
    function mapToSubsector(R::Rvec)

        if R.n1 <0
            R= ymirror(R)
        end

        if R.n2 <0
            R= xmirror(R)
        end

        if R.n2 >R.n1
            R= xymirror(R)
        end

        if !inCorrectSubsector(R)
            error("Error: Vector not correctly mapped: ",R)
        end
        return R
    end

    """given a pair of sites Rk and Rj from a vertex, return a corresponding inequivalent pair R_new by using lattice symmetries. This means R_new is in an equivalent position to site 0 as Rk is to Rj. (Rk,Rj) -> (0,R_new) """
    function pairToInequiv(Rk::Rvec,Rj::Rvec)
        Rj = translateToOrigin(Rj,Rk) #global translation to map Rk to first unit cell
        Rk = Rvec(0,0,1)
        return Rk, mapToSubsector(Rj) #give the vector in correct subsector
    end

    """
    returns Geometry struct with all relevant information about the Square lattice. We can either consider N nearest neighbor pairs or all unit cells within length L of origin.
    """
    function getSquareLattice(NLen,J = [1.,0.5];kwargs...)
        Name = string("SquareLattice_NLen=",NLen)
        System =  getLatticeGeometry(NLen,Name,pairToInequiv,inCorrectSubsector,Basis;kwargs...)

        (;PairList,couplings) = System
        setNeighborCouplings!(System,J,Basis)

        return(System)
    end
end