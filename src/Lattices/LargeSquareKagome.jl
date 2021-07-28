"""
Module to generate Lattice properties for the Square Kagome Lattice with enlarged Unit Cell within FRG
"""

module LargeSquareKagome
    using StaticArrays,Test,Parameters
    using ..SpinFRGLattices
    export getDimerSquareKagome,getNoMirrDimerSquareKagome
        
    function LargeSquareKagomeBasis()
        x = SA[1,0] 
        y = SA[0,1]
        
        # old basis: We are using these sites and shrinking the coords by a factor of two
        b1 = SA[0,1/2] / 2
        b6 = SA[1/2,1.] / 2
        
        b = sqrt(2)/4 / 2
        b2 = SA[1/4,3/4] / 2
        b3 = SA[1/4,1/4] / 2
        b4 = SA[3/4,1/4] / 2
        b5 = SA[3/4,3/4] / 2

        Cell1 = [b1,b2,b3,b4,b5,b6]
        type1 = [1,2,2,2,2,1]
        # Use originial cell to construct other sites
        Cell2 = [b + SA[1/2,0] for b in Cell1]
        type2 = [1,3,3,3,3,1]
        Cell3 = [b + SA[0,1/2] for b in Cell1]
        Cell4 = [b + SA[1/2,1/2] for b in Cell1]
        UC = vcat(Cell1,Cell2,Cell3,Cell4)
        SiteType = vcat(type1,type2,type2,type1)
        refSites = [Rvec(0,0,1),Rvec(0,0,2),Rvec(0,0,8)]
        return Basis_Struct_2D(a1=x,a2=y,b= UC,NNdist = b,NUnique = 3,SiteType= SiteType,refSites = refSites)
    end
    const Basis = LargeSquareKagomeBasis()
    getCartesian(R) = SpinFRGLattices.getCartesian(R,Basis) 
    dist(R1,R2) = SpinFRGLattices.dist(R1,R2,Basis) 
    getRvec(R)::Rvec_2D = SpinFRGLattices.getRvec(R,Basis) 
    getSiteType(R) = SpinFRGLattices.getSiteType(R,Basis)

    """
    C4 rotation around (1/4 ,1/4) for a cartesian vector
    """
    function C4(r::AbstractArray)
        rnew = SA[1/2-r[2],r[1]]
        return(rnew)
    end

    function C4(R::Rvec_2D)
        return getRvec(C4(getCartesian(R)))
    end

    """Mirrors coords at y = 1/4 (first ref site)"""
    function Mirror1(R::AbstractArray)
        res = SA[R[1],1/2-R[2]]
        return (res)
    end 

    """Mirrors coords at y = -x +1/2 (second ref site)"""
    function Mirror2(R::AbstractArray)
        res = SA[1/2-R[2],1/2-R[1]]
        return (res)
    end 

    """Mirrors coords at y = -x +1 (third ref site)"""
    function Mirror3(R::AbstractArray)
        res = SA[1-R[2],1-R[1]]
        return (res)
    end 

    function Mirror1(R::Rvec)
        getRvec(Mirror1(getCartesian(R)))
    end
    function Mirror2(R::Rvec)
        getRvec(Mirror2(getCartesian(R)))
    end
    function Mirror3(R::Rvec)
        getRvec(Mirror3(getCartesian(R)))
    end
    ##
    function aboveLine(f,R)
        r = getCartesian(R)
        y = f(r[1])
        return r[2] >= y
    end

    function inSector1(R_ref::Rvec_2D,R::Rvec_2D)
        reftype = Basis.SiteType[R_ref.b]
        mirrorLine = (x->1/4+0*x, x-> 1/2 -x, x-> 1 -x)[reftype]
        return aboveLine(mirrorLine,R)
    end

    """
    Returns pairs of Rvec so that R is in correct Sector according to all mirror symmetries
    """
    function mapToSector1(R_ref::Rvec_2D,R::Rvec_2D)
        reftype = Basis.SiteType[R_ref.b]
        mirror = (Mirror1,Mirror2,Mirror3)[reftype]
        if !inSector1(R_ref,R)
            R = mirror(R)
            R_ref = mirror(R_ref)
        end

        if !inSector1(R_ref,R)
            error("could not map R to Sector 1 ",R)
        end
        return R_ref,R
    end
    ##
    """
    Given a pair of sites Rk and Rj from a vertex, return a corresponding inequivalent pair R_new by using lattice symmetries. This means R_new is in an equivalent position to one of the reference sites x ∈ 1,2 Basis.NUnique as Rk is to Rj. (Rk,Rj) -> (Rx,R_new). Does not yet map to the correct sector!
    """
    function pairToInequiv_0(Rk::Rvec_2D,Rj::Rvec_2D)
        #Rotate until we have sites that can be translated to ref sites
        count = 1
        while !(Rk.b in (1,2,8,14,19,20)) # All these sites may be translated back to ref sites
            #gloabel C4 rotations
            Rk = C4(Rk)
            Rj = C4(Rj)
            count+=1
            if count > 4 
                error("could not map to reference site: ",(Rk,Rj))
            end
        end
        # translate so that type of Rk is in refSites
        if Rk.b == 14 #corresponds to site 3
            Rk = translation(Rk,SA[0.5,-0.5],Basis)
            Rj = translation(Rj,SA[0.5,-0.5],Basis)
        end
        if Rk.b in (19,20) #corresponds to site 1,2
            Rk = translation(Rk,SA[-0.5,-0.5],Basis)
            Rj = translation(Rj,SA[-0.5,-0.5],Basis)
        end
        ## translate to first unit cell
        Rj = translateToOrigin(Rj,Rk) # global translation so that Rj is in first unit cell
        Rk = Rvec_2D(0,0, Rk.b) 
        return Rk, Rj
    end

    function pairToInequiv(Rk::Rvec_2D,Rj::Rvec_2D)
        Rk,Rj = pairToInequiv_0(Rk, Rj)
        return mapToSector1(Rk,Rj)  # map to correct sector according to symmetry 
    end
    """
    returns Geometry struct with all relevant information about the Pyrochlore lattice. Considers only pairs withing distane R from origin.
    """

    function getDimerSquareKagome(NLen,J = (1.,1.,1.,0.),delta = 0.0;test = false)
        J1,J2,J3,Jx = J
        if J2 != J3
            error("Couplings do not fulfill assumed symmetries!")
        end
        Name = string("DimerSquagome_NLen=",NLen)
        System = getLatticeGeometry(NLen,Name,pairToInequiv,inSector1,Basis,test=test)
        @unpack couplings,PairList,PairTypes = System

        function setJ!(x,R,val)
            if val !=0.
                R_ref = Basis.refSites[x]
                R_ref,R = mapToSector1(R_ref,R)
                return setCoupling!(couplings,x,R,val,PairList,PairTypes)
            end
        end
        setJ!(1,Rvec(0,0,2),J2)
        setJ!(1,Rvec(-1,0,10),J2+delta)
        
        setJ!(1,Rvec(0,0,3),J3)
        setJ!(1,Rvec(-1, 0, 11),J3+delta)

        setJ!(2,Rvec(0,0,3),J1+delta)
        setJ!(2,Rvec(0,0,5),J1+delta)

        setJ!(2,Rvec(0,0,1),J2)
        setJ!(2,Rvec(0,0,6),J3)
        setJ!(2,Rvec(0,0,4),Jx)


        setJ!(3,Rvec(0,0,3 + 6),J1)
        setJ!(3,Rvec(0,0,5 + 6),J1)

        setJ!(3,Rvec(0,0,1 + 6),J2+delta)
        setJ!(3,Rvec(0,0,6 + 6),J3+delta)
        setJ!(3,Rvec(0,0,4 + 6),Jx)

        return(System)
    end

    ### Non-Symmetric version

    function inSector0(R_ref::Rvec_2D,R::Rvec_2D)
        return true
    end

    """
    returns Geometry struct with all relevant information about the Pyrochlore lattice. Considers only pairs withing distane R from origin.
    """
    function getNoMirrDimerSquareKagome(NLen,J = (1.,1.,1.,0.),deltaJ = (0.,0.,0.);test = false)
        J1,J2,J3,Jx = J
        deltaJ1,deltaJ2,deltaJ3 = deltaJ
        Name = string("NoMirrDimerSquagome_NLen=",NLen)
        System = getLatticeGeometry(NLen,Name,pairToInequiv_0,inSector0,Basis,test=test)
        @unpack couplings,PairList,PairTypes = System

        function setJ!(x,R,val)
            if val !=0.
                return setCoupling!(couplings,x,R,val,PairList,PairTypes)
            end
        end
        setJ!(1,Rvec(0,0,2),J2-deltaJ2)
        setJ!(1,Rvec(-1,0,10),J2+deltaJ2)
        
        setJ!(1,Rvec(0,0,3),J3-deltaJ3)
        setJ!(1,Rvec(-1, 0, 11),J3+deltaJ3)

        setJ!(2,Rvec(0,0,3),J1+deltaJ1)
        setJ!(2,Rvec(0,0,5),J1+deltaJ1)

        setJ!(2,Rvec(0,0,1),J2-deltaJ2)
        setJ!(2,Rvec(0,0,6),J3-deltaJ3)
        setJ!(2,Rvec(0,0,4),Jx)


        setJ!(3,Rvec(0,0,3 + 6),J1-deltaJ1)
        setJ!(3,Rvec(0,0,5 + 6),J1-deltaJ1)

        setJ!(3,Rvec(0,0,1 + 6),J2+deltaJ2)
        setJ!(3,Rvec(0,0,6 + 6),J3+deltaJ3)
        setJ!(3,Rvec(0,0,4 + 6),Jx)

        return(System)
    end

end

module Loop6SquareKagome
    using StaticArrays,Test,Parameters
    import ..LargeSquareKagome as LSK
    using ..SpinFRGLattices
    export getLoop6SquareKagome
    function Loop6SquareKagomeBasis()
        x = SA[1,0] 
        y = SA[0,1]
        
        # old basis: We are using these sites and shrinking the coords by a factor of two
        b1 = SA[0,1/2] / 2
        b6 = SA[1/2,1.] / 2
        
        b = sqrt(2)/4 / 2
        b2 = SA[1/4,3/4] / 2
        b3 = SA[1/4,1/4] / 2
        b4 = SA[3/4,1/4] / 2
        b5 = SA[3/4,3/4] / 2

        Cell1 = [b1,b2,b3,b4,b5,b6]
        type1 = [1,2,2,2,2,1]
        # Use originial cell to construct other sites
        Cell2 = [b + SA[1/2,0] for b in Cell1]
        type2 = [1,2,2,2,2,1]
        Cell3 = [b + SA[0,1/2] for b in Cell1]
        Cell4 = [b + SA[1/2,1/2] for b in Cell1]
        UC = vcat(Cell1,Cell2,Cell3,Cell4)
        SiteType = vcat(type1,type2,type2,type1)
        refSites = [Rvec(0,0,1),Rvec(0,0,2)]
        return Basis_Struct_2D(a1=x,a2=y,b= UC,NNdist = b,NUnique = 2,SiteType= SiteType,refSites = refSites)
    end

    const Basis = Loop6SquareKagomeBasis()
    getCartesian(R) = SpinFRGLattices.getCartesian(R,Basis) 
    dist(R1,R2) = SpinFRGLattices.dist(R1,R2,Basis) 
    getRvec(R)::Rvec_2D = SpinFRGLattices.getRvec(R,Basis) 
    getSiteType(R) = SpinFRGLattices.getSiteType(R,Basis)

    """
    C4 rotation around (1 ,1) for a cartesian vector
    """
    function C4(r::AbstractArray)
        rnew = SA[2-r[2],r[1]]
        return(rnew)
    end

    function C4(R::Rvec_2D)
        return getRvec(C4(getCartesian(R)))
    end

    function aboveLine(f,R)
        r = getCartesian(R)
        y = f(r[1])
        return r[2] >= y
    end

    function inSector(R_ref::Rvec_2D,R::Rvec_2D)
        reftype = Basis.SiteType[R_ref.b]
        x,y = getCartesian(R)
        if reftype ==1
            return aboveLine( x-> 0.25,R)
        elseif reftype == 2
            return true
        end
        error("reftype invalid: ",(R_ref,R))
    end

    function Mirror1(R::Rvec)
        getRvec(LSK.Mirror1(getCartesian(R)))
    end

    """Mirrors coords at x = 1/4 (6th basis site)"""
    function Mirror6(R::AbstractArray)
        res = SA[1/2-R[1],R[2]]
        return (res)
    end 
    function Mirror6(R::Rvec)
        getRvec(Mirror6(getCartesian(R)))
    end


    function mapToSector(R_ref::Rvec_2D,R::Rvec_2D)
        reftype = Basis.SiteType[R_ref.b]
        Rold = R
        R_refold = R_ref
        x,y = getCartesian(R)
        if reftype == 1
            if !aboveLine( x-> 0.25,R)
                R = Mirror1(R)
                R_ref = Mirror1(R_ref)
            end
        end
        if !inSector(R_ref,R)
            error("could not map R to Sector ",R_refold,Rold)
        end
        return R_ref,R
    end

    function pairToInequiv_0(Rk::Rvec_2D,Rj::Rvec_2D)
        Rk_old = Rk
        Rj_old = Rj
        #Rotate until we have sites that can be translated to ref sites
        count = 1
        while !(Rk.b in (1,2,3,4,5,7)) # All these sites may be mirrored back to ref sites
            #gloabel C4 rotations
            Rk = C4(Rk)
            Rj = C4(Rj)
            count+=1
            if count > 4 
                error("could not map to reference site: ",(Rk,Rj))
            end
        end
        Rj = translateToOrigin(Rj,Rk) # global translation so that Rj is in first unit cell
        Rk = Rvec_2D(0,0, Rk.b) 
        # mirror sites to 1,2
        if Rk.b in (4,5,7)
            Rk = Mirror6(Rk)
            Rj = Mirror6(Rj)
        end
        if Rk.b == 3 
            Rk = Mirror1(Rk)
            Rj = Mirror1(Rj)
        end
        ## translate to first unit cell
        @assert Rk.b in (1,2) "Mapping went wrong for $Rk_old , $Rj_old  mapped to $Rk , $Rj" 
        return Rk, Rj
    end

    function pairToInequiv(Rk::Rvec_2D,Rj::Rvec_2D)
        Rk,Rj = pairToInequiv_0(Rk, Rj)
        return mapToSector(Rk,Rj)  # map to correct sector according to symmetry 
    end

    function getLoop6SquareKagome(NLen,delta = 0.0,Jx = 0.;test = false)
        J = 1.
        Name = string("Loop6Squagome_NLen=",NLen)
        System = getLatticeGeometry(NLen,Name,pairToInequiv,inSector,Basis,test=test)
        @unpack couplings,PairList,PairTypes = System

        function setJ!(x,R,val)
            if val !=0.
                R_ref = Basis.refSites[x]
                R_ref,R = mapToSector(R_ref,R)
                return setCoupling!(couplings,x,R,val,PairList,PairTypes)
            end
        end
        setJ!(1,Rvec(0,0,2),J+delta)
        setJ!(1,Rvec(-1,0,10),J-delta)

        setJ!(2,Rvec(0,0,3),J-delta)
        setJ!(2,Rvec(0,0,5),J+delta)

        setJ!(2,Rvec(0,0,1),J+delta)
        setJ!(2,Rvec(0,0,6),J-delta)
        setJ!(2,Rvec(0,0,4),Jx)

        return(System)
    end

end