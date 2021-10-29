"""
Module to generate Lattice properties for the Square Kagome Lattice within FRG
"""
##
module SquareKagome
    using StaticArrays,Test,Parameters
    using ..SpinFRGLattices
    export getSquareKagome,getMirrorSquareKagome
    # "shuriken-like structure: Not used in main implementation"
    # function SquareKagomeBasis()
    #     x = SA[1,0] 
    #     y = SA[0,1]

    #     b = (√3 -1)/2
    #     c = √3 /2 * b

    #     b1 = SA[0,1/2]
    #     b2 = b1 + SA[c,b/2]
    #     b3 = b1 + SA[c,-b/2]
    #     b4 = b3 + SA[b,0]
    #     b5 = b2 + SA[b,0]
    #     b6 = SA[1/2,1.]
    #     return Basis_Struct_2D(a1=x,a2=y,b= [b1,b2,b3,b4,b5,b6],NNdist = b,NUnique = 2,SiteType= [1,2,2,2,2,1])
    # end
      
    function SquareKagomeBasis()
        x = SA[1,0] 
        y = SA[0,1]
        
        b1 = SA[0,1/2]
        b6 = SA[1/2,1.]
        
        b = sqrt(2)/4
        b2 = SA[1/4,3/4]
        b3 = SA[1/4,1/4]
        b4 = SA[3/4,1/4]
        b5 = SA[3/4,3/4]
        return Basis_Struct_2D(a1=x,a2=y,b= [b1,b2,b3,b4,b5,b6],NNdist = b,NUnique = 2,SiteType= [1,2,2,2,2,1])
    end
    const Basis = SquareKagomeBasis()
    getCartesian(R) = SpinFRGLattices.getCartesian(R,Basis) 
    dist(R1,R2) = SpinFRGLattices.dist(R1,R2,Basis) 
    getRvec(R)::Rvec_2D = SpinFRGLattices.getRvec(R,Basis) 
    getSiteType(R) = SpinFRGLattices.getSiteType(R,Basis)

    """
    C4 rotation around (1/2 ,1/2) for a cartesian vector
    """
    function C4(r::AbstractArray)
        rnew = SA[1-r[2],r[1]]
        return(rnew)
    end

    """
    Inversion at (0,1/2)
    """
    function inversion(r::AbstractArray)
        rnew = SA[-r[1],1-r[2]]
        return rnew
    end

    function inversion(R::Rvec_2D)
        return getRvec(inversion(getCartesian(R)))
    end
    
    """
    Using inversion symmetry, for the reference site b = 1 we only consider pair partners with positive values of n1.
    """
    function inSector1(R_ref::Rvec_2D,R::Rvec_2D)
        if R_ref.b == 1
            if R.n1 <0 || R.n1 == 0 && R.b == 1 && R.n2 <0 # second part: sites on n1=0 are also symmetric to each other
                return false
            end
        end
        return true # for other ref site, we are always in sector1
    end

    """
    Returns pairs of Rvec so that R is in Sector with positive x
    """
    function mapToSector1(R_ref::Rvec_2D,R::Rvec_2D)
        
        if !inSector1(R_ref,R)
            R = inversion(R)
            R_ref = inversion(R_ref)
        end

        if !inSector1(R_ref,R)
            error("could not map R to Sector 1 ",R)
        end
        return R_ref,R
    end

    function UnitCellTest(func, x = 0,y = 0)
        for b in 1:Basis.NCell
            Rstart = Rvec_2D(x,y,b)
            Rnew = getRvec(func(getCartesian(Rstart)))
            println([Rnew.n1,Rnew.n2,Rnew.b])   
            # Rnew = func(Rstart)
            # println([Rnew.n1,Rnew.n2,Rnew.b])
        end
    end

    """
    C4 rotation around (1/2 ,1/2) for Rvec
    """
    function C4(R::Rvec_2D)
        return getRvec(C4(getCartesian(R)))
    end


    function UnitTest()
        @testset "C4 cartesian check" begin 
            for x in -4:4, y in -4:4, b in 1:Basis.NCell
                Rstart = Rvec_2D(x,y,b)
                R_direct = C4(Rstart)
                Rtest = getRvec(C4(getCartesian(Rstart))) # cartesian rotation as consistency test
                @test R_direct == Rtest
                @test  (C4 ∘ C4 ∘ C4∘ C4)(Rstart) == Rstart
            end
        end
    end

    """
    Given a pair of sites Rk and Rj from a vertex, return a corresponding inequivalent pair R_new by using lattice symmetries. This means R_new is in an equivalent position to one of the reference sites x ∈ 1,2 Basis.NUnique as Rk is to Rj. (Rk,Rj) -> (Rx,R_new) 
    """
    function pairToInequiv_0(Rk::Rvec_2D,Rj::Rvec_2D)
        count = 1
        while !(Rk.b in (1,2)) # sites 1 and 2 are the inequivalent reference sites
            #gloabel C4 rotations
            Rk = C4(Rk)
            Rj = C4(Rj)
            count+=1
            if count > 4 
                error("could not map to reference site: ",(Rk,Rj))
            end
        end
        Rj = Rvec_2D(Rj.n1-Rk.n1, Rj.n2-Rk.n2, Rj.b) # global translation so that Rj is in first unit cell
        Rk = Rvec_2D(0,0, Rk.b) 
        return (Rk,Rj)
    end
    function pairToInequiv(Rk::Rvec_2D,Rj::Rvec_2D)
        Rk,Rj = pairToInequiv_0(Rk, Rj) 
        return mapToSector1(Rk, Rj) # map to correct sector according to symmetry 
    end
    """
    returns Geometry struct with all relevant information about the Pyrochlore lattice. Considers only pairs withing distane R from origin.
    """

    function getSquareKagome(NLen,J = (1.,1.,1.,0.);test = false)
        Name = string("Squagome_NLen=",NLen)
        System = getLatticeGeometry(NLen,Name,pairToInequiv,inSector1,Basis,test=test)
        @unpack couplings,PairList,PairTypes = System

        function setJ!(x,R,val)
            if val !=0.
                R_ref = Rvec(0,0,x)
                R_ref,R = mapToSector1(R_ref,R)
                return setCoupling!(couplings,x,R,val,PairList,PairTypes)
            end
        end
        J1,J2,J3,Jx = J
        setJ!(1,Rvec(0,0,2),J2)
        setJ!(1,Rvec(-1,0,4),J2)
        
        setJ!(1,Rvec(0,0,3),J3)
        setJ!(1,Rvec(-1,0,5),J3)

        setJ!(2,Rvec(0,0,3),J1)
        setJ!(2,Rvec(0,0,5),J1)

        setJ!(2,Rvec(0,0,1),J2)
        setJ!(2,Rvec(0,0,6),J3)
        setJ!(2,Rvec(0,0,4),Jx)

        return(System)
    end

    #using Mirror symmetries
    
    """Mirrors coords at x = 0 (first ref site)"""
    function Mirrorx1(R::AbstractArray)
        res = SA[-R[1],R[2]]
        return (res)
    end

    """Mirrors coords at y = 0.5 (first ref site)"""
    function Mirrory1(R::AbstractArray)
        res = SA[R[1],1-R[2]]
        return (res)
    end

    """Mirrors coords at y = -x +1 (second ref site)"""
    function Mirror2(R::AbstractArray)
        res = SA[1-R[2],1-R[1]]
        return (res)
    end 

    
    function Mirrory1(R::Rvec)
        getRvec(Mirrory1(getCartesian(R)))
    end
    function Mirrorx1(R::Rvec)
        getRvec(Mirrorx1(getCartesian(R)))
    end
    function Mirror2(R::Rvec)
        getRvec(Mirror2(getCartesian(R)))
    end
    ##
    function aboveLine(f,R)
        r = getCartesian(R)
        y = f(r[1])
        return r[2] >= y
    end

    function inMirrorSector(R_ref::Rvec_2D,R::Rvec_2D)
        reftype = Basis.SiteType[R_ref.b]
        x,y = getCartesian(R)
        if reftype ==1
            return x>=0 && aboveLine( x-> 0.5,R)
        elseif reftype == 2
            return aboveLine( x-> 1 -x,R)
        end
        error("reftype invalid: ",(R_ref,R))
    end

    """
    Returns pairs of Rvec so that R is in correct Sector according to all mirror symmetries
    """
    function mapToMirrorSector(R_ref::Rvec_2D,R::Rvec_2D)
        reftype = Basis.SiteType[R_ref.b]
        Rold = R
        x,y = getCartesian(R)
        if reftype == 1
            if x < 0.0
                R = Mirrorx1(R)
                R_ref = Mirrorx1(R_ref)
            end
            if y < 0.5
                R = Mirrory1(R)
                R_ref = Mirrory1(R_ref)
            end
        elseif reftype == 2
            if !aboveLine( x-> 1 -x,R)
                R = Mirror2(R)
                R_ref = Mirror2(R_ref)
            end
        end
        if !inMirrorSector(R_ref,R)
            error("could not map R to Sector ",R_ref,Rold)
        end
        return R_ref,R
    end

    function MirrorpairToInequiv(Rk::Rvec_2D,Rj::Rvec_2D)
        Rk,Rj = pairToInequiv_0(Rk, Rj) 
        return mapToMirrorSector(Rk, Rj) # map to correct sector according to symmetry 
    end

    function getMirrorSquareKagome(NLen,J1=1,J2_J3=1,Jx=0;test = false)
        Name = string("MSquagome_NLen=",NLen)
        System = getLatticeGeometry(NLen,Name,MirrorpairToInequiv,inMirrorSector,Basis,test=test)
        @unpack couplings,PairList,PairTypes = System

        function setJ!(x,R,val)
            if val !=0.
                R_ref = Rvec(0,0,x)
                R_ref,R = mapToMirrorSector(R_ref,R)
                return setCoupling!(couplings,x,R,val,PairList,PairTypes)
            end
        end 

        setJ!(1,Rvec(0,0,2),J2_J3)
        setJ!(2,Rvec(0,0,1),J2_J3) 
        setJ!(2,Rvec(0,0,3),J1)

        setJ!(2,Rvec(0,0,4),Jx)

        return(System)
    end
end