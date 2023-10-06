module Honeycomb
    using ..SpinFRGLattices
    using StaticArrays,StructArrays
    
    export getHoneycomb,getHoneycombGamma
    
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
    function getHoneycomb(NLen,J;kwargs...)
        Name = string("Honeycomb_NLen=",NLen)
        System =  getLatticeGeometry(NLen,Name,pairToInequiv,inCorrectSubsector,Basis;kwargs...)

        setNeighborCouplings!(System,J,Basis)

        return(System)
    end
    
    ## cluster Hexagon couplings
    using .Octochlore: spin, reduceCouplings,getInequivCouplings,mapCouplingsToSiteList
    import .Octochlore: getInequivCouplings
    function getInequivCouplings(alpha::Number,gamma::Number)
        a= reduceCouplings(getCouplingsToS1(alpha,gamma))
        getInequivCouplings(a,mapToSubsector)
    end

    function getHexas(refSite::Rvec,alpha = 1., gamma= 0.)
        Allsites = generateLUnitCells(2,Basis,refSite)
        center = getCartesian(refSite) + 0.5 *(Basis.a1+Basis.b[2])
        centerNorm(x) = norm(getCartesian(x) - center)
        AllNorms = sort(unique(round.(centerNorm.(Allsites),digits = 14)))
        getNeighorFromCenter(i) = filter(x->centerNorm(x) ≈ AllNorms[i],Allsites)
        firstHexa = getNeighorFromCenter(1)
        secondHexa = getNeighorFromCenter(2)
        return (StructArray(append!([spin(alpha,site) for site in firstHexa],[spin(gamma,site) for site in secondHexa])))
    end

    R(x::Integer) = Rvec(0,0,x)
    function getCouplingsToS1(alpha::T = 1., gamma::T= 0. ) where T <: Number
        L=2
        S1Terms = spin{T,Rvec_2D}[]
        for n1 in -L:L,n2 in -L:L
            hexas = getHexas(Rvec(n1,n2,1),alpha,gamma)
            for s1 in hexas
                for s2 in hexas
                    if s1.site == R(1)
                        push!(S1Terms,spin(s1.fac*s2.fac,s2.site))
                    elseif s2.site == R(1)
                        push!(S1Terms,spin(s1.fac*s2.fac,s1.site))
                    end
                end
            end
        end
        return StructArray(S1Terms)
    end

    function getHoneycombGamma(NLen;alpha=1.,gamma=0.,test = false,normalize = true)
        System = getHoneycomb(NLen,[0.,0.])
        (;PairList,couplings) = System
        IneqCouplings = getInequivCouplings(Float64(alpha),Float64(gamma))
        couplings .= 0.5 .*mapCouplingsToSiteList(IneqCouplings,PairList)
        if normalize
            couplings[System.OnsitePairs] .= 0.
            largestCoupling = maximum(abs,couplings)
            couplings ./= largestCoupling
        end
        if test 
            testGeometry(System)
        end
        return System
    end

end
##
