##
module SimpleCubic
    using Parameters,StaticArrays,Test
    using ..SpinFRGLattices
    export getCubic
    ##
    """Constructer for a square lattice basis"""
    function CubicBasis()
        a1 = SA[1,0,0]
        a2 = SA[0,1,0]
        a3 = SA[0,0,1]
        
        b1 = SA[0,0,0]
        return Basis_Struct_3D(a1=a1,a2=a2,a3=a3,b=[b1], NNdist = 1.)
    end
    const Basis = CubicBasis()
    getCartesian(R) = SpinFRGLattices.getCartesian(R,Basis) 
    dist(R1,R2) = SpinFRGLattices.dist(R1,R2,Basis) 
    getRvec(R) = SpinFRGLattices.getRvec(R,Basis) 
    ##

    """performs C4 rotation on R around z"""
    function C4(R::Rvec)
        return Rvec(-R.n2,R.n1,R.n3,R.b)
    end


    """mirrors R on y-z plane"""
    function xmirror(R::Rvec)
        return Rvec(-R.n1,R.n2,R.n3,R.b)
    end

    """mirrors R on x-z plane"""
    function ymirror(R::Rvec)
        return Rvec(R.n1,-R.n2,R.n3,R.b)
    end

    """mirrors R on x-y plane"""
    function zmirror(R::Rvec)
        return Rvec(R.n1,R.n2,-R.n3,R.b)
    end

    """mirrors R on y=x plane"""
    function xymirror(R::Rvec)
        return Rvec(R.n2,R.n1,R.n3,R.b)
    end

    """mirrors R on z=y plane"""
    function yzmirror(R::Rvec)
        return Rvec(R.n1,R.n3,R.n2,R.b)
    end

    """mirrors R on z=x plane"""
    function xzmirror(R::Rvec)
        return Rvec(R.n3,R.n2,R.n1,R.b)
    end

    """
    Symmetry: may restrict to section with z<=y<=x
    """
    function inCorrectSubsector(R_ref,R::Rvec)
        @unpack n1,n2,n3 = R
        return(n1>= 0 && n2>=0 && n3 >=0&& n2 <=n1 && n3 <=n2 )
    end

    inCorrectSubsector(R::Rvec) = inCorrectSubsector(Rvec(0,0,0,1),R)


    """gives coordinates of R' after using symmetry transformations on R
    R is arbitrary lattice vector"""
    function mapToSubsector(R::Rvec)
        Rold = R
        if R.n1 <0
            R= xmirror(R)
        end

        if R.n2 <0
            R= ymirror(R)
        end

        if R.n3 <0
            R= zmirror(R)
        end

        if R.n3 >R.n2
            R= yzmirror(R)
        end
        
        if R.n2 >R.n1
            R= xymirror(R)
        end

        if R.n3 >R.n2
            R= yzmirror(R)
        end
        if !inCorrectSubsector(R)
            error("Error: Vector not correctly mapped: ",Rold)
        end
        return R
    end

    """given a pair of sites Rk and Rj from a vertex, return a corresponding inequivalent pair R_new by using lattice symmetries. This means R_new is in an equivalent position to site 0 as Rk is to Rj. (Rk,Rj) -> (0,R_new) """
    function pairToInequiv(Rk::Rvec,Rj::Rvec)
        Rj = translateToOrigin(Rj,Rk) #global translation to map Rk to first unit cell
        Rk = Rvec(0,0,0,1)
        return Rk, mapToSubsector(Rj) #give the vector in correct subsector
    end

    """
    returns Geometry struct with all relevant information about the Square lattice. We can either consider N nearest neighbor pairs or all unit cells within length L of origin.
    """
    function getCubic(NLen,J = [1.,0.5];kwargs...)
        Name = string("Cubic_NLen=",NLen)
        System =  getLatticeGeometry(NLen,Name,pairToInequiv,inCorrectSubsector,Basis;kwargs...)

        @unpack PairList,couplings = System
        setNeighborCouplings!(couplings,J,PairList,Basis)
        return(System)
    end
end