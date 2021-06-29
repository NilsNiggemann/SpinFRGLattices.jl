
module C2SquareKagome
    using StaticArrays,Test,Parameters
    using ..SpinFRGLattices
    export getC2DimerSquareKagome
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
        type1 = [1,2,3,2,3,4]
        # Use originial cell to construct other sites
        Cell2 = [b + SA[1/2,0] for b in Cell1]
        type2 = [1,5,6,5,6,4]
        Cell3 = [b + SA[0,1/2] for b in Cell1]
        Cell4 = [b + SA[1/2,1/2] for b in Cell1]
        UC = vcat(Cell1,Cell2,Cell3,Cell4)
        SiteType = vcat(type1,type2,type2,type1)
        refSites = [Rvec(0,0,1),Rvec(0,0,2),Rvec(0,0,3),Rvec(0,0,6),Rvec(0,0,8),Rvec(0,0,9)]
        return Basis_Struct_2D(a1=x,a2=y,b= UC,NNdist = b,NUnique = 6,SiteType= SiteType,refSites = refSites)
    end
    const Basis = LargeSquareKagomeBasis()
    getCartesian(R) = SpinFRGLattices.getCartesian(R,Basis) 
    dist(R1,R2) = SpinFRGLattices.dist(R1,R2,Basis) 
    getRvec(R)::Rvec_2D = SpinFRGLattices.getRvec(R,Basis) 
    getSiteType(R) = SpinFRGLattices.getSiteType(R,Basis)
    translation(R,T) = SpinFRGLattices.translation(R,T,Basis)

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

    function C2(R::Rvec_2D)
        return getRvec(C4(C4(getCartesian(R))))
    end
    ##

    function pairToInequiv(Rk::Rvec_2D,Rj::Rvec_2D)
        #Rotate until we have sites that can be translated to ref sites
        count = 1
        refSites = (1,2,3,6,8,9)
        translSites = (19,20,21,24,14,15)
        while !(Rk.b in refSites || Rk.b in translSites) # All these sites may be translated back to ref sites
            #gloabel C2 rotations
            Rk = C2(Rk)
            Rj = C2(Rj)
            count+=1
            if count > 2
                error("could not map to reference site: ",(Rk,Rj))
            end
        end
        # translate so that type of Rk is in refSites

        if Rk.b in (14,15) #corresponds to site 3
            Rk = translation(Rk,SA[0.5,-0.5])
            Rj = translation(Rj,SA[0.5,-0.5])
        elseif Rk.b in (19,20,21,24)  #corresponds to fourth subCell
            Rk = translation(Rk,SA[-0.5,-0.5])
            Rj = translation(Rj,SA[-0.5,-0.5])
        end
        if !(Rk.b in refSites)
            error("mapping went wrong:", Rk)
        end
        ## translate to first unit cell
        Rj = translateToOrigin(Rj,Rk) # global translation so that Rj is in first unit cell
        Rk = Rvec_2D(0,0, Rk.b) 
        return Rk, Rj
    end

    function inSector0(R_ref::Rvec_2D,R::Rvec_2D)
        if R_ref in Basis.refSites
            return true
        end
        error("called inSector0 without proper refSite: ",(R_ref,R))
    end


    """
    returns Geometry struct with all relevant information about the Pyrochlore lattice. Considers only pairs withing distane R from origin.
    """
    function getC2DimerSquareKagome(NLen,J = (1.,1.,1.,1.,0.),deltaJ = (0.,0.,0.,0.);test = false)
        J1Light,J1Dark,J2,J3,Jx = J
        deltaJ1Light,deltaJ1Dark,deltaJ2,deltaJ3 = deltaJ
        Name = string("C2DimerSquagome_NLen=",NLen)
        System = getLatticeGeometry(NLen,Name,pairToInequiv,inSector0,Basis,test=test)
        @unpack couplings,PairList,PairTypes = System

        function setJ!(x,R,val)
            #x is refSite number!
            if val !=0.
                return setCoupling!(couplings,x,R,val,PairList,PairTypes)
            end
        end

        function setJBond!(b1,T::StaticArray,val)
            #b1 is basis of refSite!
            R_Ref = Rvec(0,0,b1)
            x = getSiteType(R_Ref)
            R2 = translation(R_Ref,T)
            R_Ref,R2 = pairToInequiv(R_Ref,R2)
            if R_Ref in Basis.refSites || R2 in Basis.refSites
                setJ!(x,R2,val)
            end
        end
        
        function setJPair!(b1,T::StaticArray,val,delta)
            setJBond!(b1,T,val-delta)
            b2 = translation(Rvec(0,0,b1),SA[0.5,0.]).b # translate to second SubCell
            setJBond!(b2,T,val+delta)
        end

        UpRight = SA[1/8,1/8]
        DownRight = SA[1/8,-1/8]
        DarkGreen = SA[0,-1/4]
        LightGreen = SA[1/4,0]
        GreyUp = SA[1/4,1/4]
        GreyDown = SA[1/4,-1/4]
        
        setJPair!(1,UpRight,J2,deltaJ2)
        setJPair!(1,-UpRight,J2,-deltaJ2)

        setJPair!(1,DownRight,J3,deltaJ3)
        setJPair!(1,-DownRight,J3,-deltaJ3)

        setJPair!(2,UpRight,J3,deltaJ3)
        setJPair!(2,-UpRight,J2,deltaJ2)
        
        setJPair!(2,LightGreen,J1Light,-deltaJ1Light)
        setJPair!(2,DarkGreen,J1Dark,-deltaJ1Dark)
        
        setJPair!(2,GreyDown,Jx,0.)
        
        setJPair!(3,DownRight,J2,deltaJ2) 
        setJPair!(3,-DownRight,J3,deltaJ3)

        setJPair!(3,-DarkGreen,J1Dark,-deltaJ1Dark)
        setJPair!(3,LightGreen,J1Light,-deltaJ1Light)

        setJPair!(3,GreyUp,Jx,0.)

        setJPair!(6,UpRight,J3,-deltaJ3)
        setJPair!(6,-UpRight,J3,deltaJ3)

        setJPair!(6,DownRight,J2,deltaJ2)
        setJPair!(6,-DownRight,J2,-deltaJ2)

        return(System)
    end
end