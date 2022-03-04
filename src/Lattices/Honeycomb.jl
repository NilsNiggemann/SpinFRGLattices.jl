module Honeycomb
    using ..SpinFRGLattices
    using StaticArrays,Parameters,StructArrays
    
    export getHoneycombLattice
    
    function HoneycombBasis()
        a1 = SA[1,0]
        a2 = SA[-1/2,√(3/4)]
        
        b0 = SA[0,0]
        b1 = SA[0,1/√3]
        return Basis_Struct_2D(a1=a1,a2=a2,b=[b0,b1],NNdist = norm(b1))
    end
    const Basis = HoneycombBasis()

    getCartesian(R) = SpinFRGLattices.getCartesian(R,Basis) 
    dist(R1,R2) = SpinFRGLattices.dist(R1,R2,Basis) 
    getRvec(R) = SpinFRGLattices.getRvec(R,Basis)

    
    yMirror(r::AbstractVector) = SA[-r[1],r[2]] 
    yMirror(R::Rvec_2D) = R |> getCartesian |> yMirror |> getRvec
    
    xMirror(r::AbstractVector) = SA[r[1],-r[2]]+Basis.b[2] 
    xMirror(R::Rvec_2D) = R |> getCartesian |> xMirror |> getRvec

    downline(x) = -1/sqrt(3)*x
    downMirror(x) = MirrorLine(-1/sqrt(3),x,Basis)

    upline(x) = 1/sqrt(3)*x
    upMirror(x) = MirrorLine(1/sqrt(3),x,Basis)

    function inRightHalf(R::Rvec)
        r = getCartesian(R)
        r[1] >= 0.
    end

    function inHexant3(R::Rvec)
        inRightHalf(R) && belowLine_strict(downline,R,Basis)
    end

    function inHexant2(R::Rvec)
       !belowLine_strict(downline,R,Basis) && belowLine_strict(upline,R,Basis)
    end

    function inHexant1(R::Rvec)
        inRightHalf(R) && !belowLine_strict(upline,R,Basis)
    end

    function inCorrectSubsector(R_ref,R::Rvec)
        inHexant1(R)
    end

    inCorrectSubsector(R::Rvec) = inHexant1(R)

    """gives coordinates of R' after using symmetry transformations on R which leave the reference site invariant.
    R is arbitrary lattice vector"""
    function mapToSubsector(R::Rvec)
        i = 1
        Rold = R
        !(inRightHalf(R)) && (R = yMirror(R))
        inHexant3(R) && (R = downMirror(R))
        inHexant2(R) && (R = upMirror(R))

        inCorrectSubsector(R) || error("Error: Vector not correctly mapped: ",Rold)
        return R
    end

    """given a pair of sites Rk and Rj from a vertex, return a corresponding inequivalent pair R_new by using lattice symmetries. This means R_new is in an equivalent position to site 0 as Rk is to Rj. (Rk,Rj) -> (0,R_new) """
    function pairToInequiv(Rk::Rvec,Rj::Rvec)
        if Rk.b != 1
            Rk = xMirror(Rk)
            Rj = xMirror(Rj)
        end
        Rj = translateToOrigin(Rj,Rk) #global translation to map Rk to first unit cell
        Rk = Rvec(0,0,1)
        return Rk, mapToSubsector(Rj) #give the vector in correct subsector
    end

    """
    returns Geometry struct with all relevant information about the Honeycomb lattice. We can either consider N nearest neighbor pairs or all unit cells within length L of origin.
    """
    function getHoneycombLattice(NLen,J = [1.,0.5];kwargs...)
        Name = string("Honeycomb_NLen=",NLen)
        System =  getLatticeGeometry(NLen,Name,pairToInequiv,inCorrectSubsector,Basis;kwargs...)

        @unpack PairList,couplings = System
        setNeighborCouplings!(couplings,J,PairList,Basis)
        return(System)
    end

end
##
