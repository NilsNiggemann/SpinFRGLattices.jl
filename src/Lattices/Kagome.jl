module KagomeLattice
    using StaticArrays,Test
    using ..SpinFRGLattices
    """Constructer for a square lattice basis"""
    function KagomeBasis()
        a1 = SA[1,0]
        a2 = 1/2*SA[1,√3]
        
        b1 = SA[0,0]
        b2 = 1/2 *a1
        b3 = 1/2 *a2
        return Basis_Struct_2D(a1=a1,a2=a2,b=[b1,b2,b3], NNdist = 0.5)
    end
    const Basis = KagomeBasis()
    # Basis = KagomeBasis()
    getCartesian(R) = SpinFRGLattices.getCartesian(R,Basis) 
    dist(R1,R2) = SpinFRGLattices.dist(R1,R2,Basis) 
    getRvec(R) = SpinFRGLattices.getRvec(R,Basis)

    """performs C3 rotation on R"""
    function C3_Origin(r::AbstractVector)
        SA[-1/2 -sqrt(3)/2; sqrt(3)/2 -1/2] * r
    end

    """performs C3 rotation on R"""
    function C3_Origin(R::Rvec)
        return R |> getCartesian |> C3 |> getRvec
    end

    C3(r::AbstractVector) =  r|> x-> x- (Basis.b[2]+Basis.b[3])/3 |> C3_Origin |> x-> x+(Basis.b[2]+Basis.b[3])/3

    C3(R::Rvec) = R |> getCartesian |> C3 |> getRvec
    # C3(R::Rvec) =  R |> getCartesian |> x-> x-(SA[0.5,0.5])  |> C3_Origin |> getRvec
    inCorrectSubsector_0(R1::Rvec,R2::Rvec) = true
    """Pair to inequiv: Assumes no lattice symmetries"""
    function pairToInequiv_0(Rk::Rvec,Rj::Rvec)
        while Rk.b != 1
            Rk = C3(Rk)
            Rj = C3(Rj)
        end
        Rj = translateToOrigin(Rj,Rk) #global translation to map Rk to first unit cell
        Rk = Rvec(0,0,1)
        return Rk, Rj #give the vector in correct subsector
    end

    # System = getKagomeLattice(3,test = true)
    # plotSystem(System,Basis)

end

module LargeKagomeLattice

    using StaticArrays,Test
    using ..SpinFRGLattices
    export getHalfMoonDimerKagomeLattice
    """Constructer for a square lattice basis"""
    function LargeKagomeBasis()
        a1 = 4*SA[1,0] 
        a2 = 4*1/2*SA[1,√3]
        
        b1 = SA[0,0]
        b2 = 1/4 *a1
        b3 = 1/4 *a2

        b4 = 2/4 *a1
        b5 = 3/4 *a1
        b6 = 2/4 *a1 + 1/4*a2
        b7 = 2/4 *a2
        b8 = 1/4 *a1 + 2/4*a2
        b9 = 3/4 *a2
        b10 = 2/4 *a1 + 2/4*a2
        b11 = 3/4 *a1 + 2/4*a2
        b12 = 2/4 *a1 + 3/4*a2

        return Basis_Struct_2D(a1=a1,a2=a2,b=[b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12], NNdist = norm(b2),SiteType = [1,2,2,2,1,2,2,2,1,1,1,1],NUnique = 2)
    end
    const Basis = LargeKagomeBasis()

    getCartesian(R) = SpinFRGLattices.getCartesian(R,Basis) 
    dist(R1,R2) = SpinFRGLattices.dist(R1,R2,Basis) 
    getRvec(R) = SpinFRGLattices.getRvec(R,Basis)

    """performs C6 rotation on R"""
    function C6_Origin(r::AbstractVector)
        SA[1/2 -sqrt(3)/2; sqrt(3)/2 1/2] * r
    end

    """performs C6 rotation on R"""
    function C6_Origin(R::Rvec)
        return R |> getCartesian |> C6 |> getRvec
    end

    C6(r::AbstractVector) =  r|> x-> x- (Basis.a1+Basis.a2)/4 |> C6_Origin |> x-> x+(Basis.a1+Basis.a2)/4

    C6(R::Rvec) = R |> getCartesian |> C6 |> getRvec
    mapToSector(R1::Rvec,R2::Rvec) = R1,R2

    inCorrectSubsector_0(R1::Rvec,R2::Rvec) = true
    """Pair to inequiv: Assumes no lattice symmetries"""
    function pairToInequiv(Rk::Rvec,Rj::Rvec)
        while !(Rk.b in (1,2))
            Rk = C6(Rk)
            Rj = C6(Rj)
        end
        Rj = translateToOrigin(Rj,Rk) #global translation to map Rk to first unit cell
        Rk = Rvec(0,0,Rk.b)
        return Rk, Rj #give the vector in correct subsector
    end


    function getHalfMoonDimerKagomeLattice(NLen;J2, J3a, J1_2,J1_3,J1_5,J1_9,J2_3,J2_4,J2_12,kwargs...)
        Name = string("HalfMoonDimerKagomeLattice_NLen=",NLen)
        System =  getLatticeGeometry(NLen,Name,pairToInequiv,inCorrectSubsector_0,Basis;kwargs...)

        (;PairList,couplings,PairTypes) = System

        function setJ!(x,R,val)
            if val !=0.
                R_ref = Basis.refSites[x]
                R_ref,R = mapToSector(R_ref,R)
                return setCoupling!(couplings,x,R,val,PairList,PairTypes)
            end
        end

        #first and second nearest neighbors 
        setNeighborCouplings!(System,[0,J2])
        
        #third nearest neighbors  J3a
        setJ!(1,Rvec(0,0,4),J3a)
        setJ!(1,Rvec(0,0,7),J3a)
        setJ!(1,Rvec(-1,0,4),J3a)
        setJ!(1,Rvec(0,-1,7),J3a)

        setJ!(2,Rvec(0,0,5),J3a)
        setJ!(2,Rvec(-1,0,11),J3a)
        setJ!(2,Rvec(-1,0,5),J3a)
        setJ!(2,Rvec(0,-1,11),J3a)

        setJ!(1,Rvec(0,0,2),J1_2)
        setJ!(1,Rvec(0,0,3),J1_3)
        setJ!(1,Rvec(-1,0,5),J1_5)
        setJ!(1,Rvec(0,-1,9),J1_9)

        setJ!(2,Rvec(0,0,1),J1_2)

        setJ!(2,Rvec(0,0,3),J2_3)
        setJ!(2,Rvec(0,0,4),J2_4)
        setJ!(2,Rvec(0,-1,12),J2_12)

        return(System)
    end


    getHalfMoonDimerKagomeLattice(NLen,J1,J1a=J1,J2=-1,J3a = J2,delta=0.01;kwargs...) = getHalfMoonDimerKagomeLattice(NLen,
    J2 = J2,
    J3a = J3a,
    J1_2=J1+delta,
    J1_3=J1a-delta,
    J1_5=J1-delta,
    J1_9=J1-delta,
    J2_3 = J1-delta,
    J2_4 = J1-delta,
    J2_12 = J1a-delta;kwargs...)


    """Helper function for pretty plotting"""
    function bondColFunction(J::Real)
        color = :red
        label = string(round(J,digits = 2))
        if J <0
            label = "\$J_{1a}\$"
            color = :blue
        elseif J == 0.25
            label = "\$J_{3a}\$"
            color = :darkgreen
        elseif J ==0.5 
            label = "\$J_2\$"
            color = :lightgreen
        elseif J >0.5 && J<1 
            label = "\$J_1- \\delta\$"
        elseif J ==1 
            label = "\$J_1\$"
        elseif J > 1
            label = "\$J_1 + \\delta \$"
        end
        return Dict([:color => color,:label => label])
    end
end
# K = LK.getDimerKagomeLattice(4,J1_2=3,J1_3=-0.5,J1_5=1,J1_9=1,J2_3 = 1,J2_4 = 1,J2_12 = -0.5)