##
module Pyrochlore
    using Parameters,StaticArrays,Test
    using ..SpinFRGLattices
    export getPyrochlore
    
    function PyrochloreBasis()
        a1 = 1/2*SA[1,1,0]
        a2 = 1/2*SA[0,1,1]
        a3 = 1/2*SA[1,0,1]
        
        b1 = SA[0,0,0]
        b2 = 1/2*a1
        b3 = 1/2*a2
        b4 = 1/2*a3
        return Basis_Struct_3D(a1=a1,a2=a2,a3=a3,b=[b1,b2,b3,b4],NNdist = norm(b2))
    end
    const Basis = PyrochloreBasis()
    getCartesian(R) = SpinFRGLattices.getCartesian(R,Basis) 
    dist(R1,R2) = SpinFRGLattices.dist(R1,R2,Basis) 
    getRvec(R)::Rvec_3D = SpinFRGLattices.getRvec(R,Basis) 
    ##

    """performs C3 rotation on R"""
    function C3(R::Rvec_3D)
        basis = (1,3,4,2)[R.b] # rotation -> cycling permutation of sublattices

        return Rvec(R.n3,R.n1,R.n2,basis)
    end

    """performs a C2 rotation on R"""
    function C2(R::Rvec_3D)
        @unpack n1,n2,n3,b = R
        b = (2,1,4,3)[b] # C2 swaps b2 with b1 and b3 with b4
        return Rvec(-n1-n2-n3,n3,n2,b)
    end


    """gives -R in lattice coordinates"""
    function inversion(R::Rvec_3D)
        return getRvec(-getCartesian(R))
    end


    """
    inversion symmetry: may restrict to one half of space (separated by (111) plane)
    """
    function inPosHalf(R::Rvec_3D)
        @unpack n1,n2,n3 = R
        return (n1+n2+n3>=0)
    end


    """
    Returns sector according to C3 rotation symmetry around (1,1,1) axis
    """
    function whichSector(R::Rvec_3D)
        x = getLatticeVec(R,Basis)
        if x[1] >= x[2] &&  x[1] >= x[3]
            return 1
        elseif x[2] >= x[3] &&  x[2] > x[1]
            return 2
        elseif x[3] >= x[1] &&  x[3] > x[2]
            return 3
        # elseif x[3] == x[2] && x[2] == x[1]
        #     return 1
        end
        println("Vector in no sector:", R)
        error("error in sector specification")
    end

    """mirrors R at the x=y plane """
    function mirror(r::AbstractArray)
        x,y,z = r
        return SA[y,x,z]
    end

    function mirror(R::Rvec_3D)
        basis = R.b
        if basis == 3
            new_b = 4
        elseif basis == 4
            new_b = 3
        else
            new_b = basis
        end
        return Rvec(R.n1,R.n3,R.n2,new_b)
    end

    """gives subsector according to mirror symmetry"""
    function whichSubSector(R::Rvec_3D)
        x = getLatticeVec(R,Basis)
        if x[2] <= x[3]
            return 1
        else
            return 2
        end
    end

    """Checks if R is in correct Subsector as specified by symmetry. We specify R_ref to be able to use SpinFRGLattices deleteEquivalent! method"""
    function inCorrectSubsector(R_ref,R::Rvec_3D)
        return(inPosHalf(R) && whichSector(R) == 1 && whichSubSector(R) == 1)
    end

    inCorrectSubsector(R::Rvec_3D) = inCorrectSubsector(Rvec_3D(0,0,0,1),R)


    """gives coordinates of R' after using symmetry transformations on R
    R is arbitrary lattice vector"""
    function mapToSubsector(R::Rvec_3D)
        R_prime = R
        if !(inPosHalf(R))
            R_prime = inversion(R_prime)
        end

        i = 0
        while whichSector(R_prime) != 1
            R_prime = C3(R_prime)
            i+=1
            if i>2
                error("error: C3 rotation not could not map to correct sector ",R_prime)
            end
        end

        i = 0
        while whichSubSector(R_prime) != 1
            R_prime = mirror(R_prime)
            i+=1
            if i>1
                error("error: mirror symmetry not could not map to correct subsector: ",R_prime)
            end
        end
        if !inCorrectSubsector(R_prime)
            error("Error: Vector not correctly mapped: ",R_prime)
        end
        return R_prime
    end

    """given a pair of sites Rk and Rj from a vertex, return a corresponding inequivalent pair R_new by using lattice symmetries. This means R_new is in an equivalent position to site 0 as Rk is to Rj. (Rk,Rj) -> (0,R_new) """
    function pairToInequiv(Rk,Rj)
        if Rk.b == 1 #left site is already reference site 1 -> we're done
            Rj = translateToOrigin(Rj,Rk) #global translation to map Rk to first unit cell
            Rk = Rvec(0,0,0,1)
            return Rk, mapToSubsector(Rj)
        end
        # otherwise: global C3 rotation until lhs is on sublattice B
        i = 0
        while Rk.b != 2
            i+=1
            Rk = C3(Rk)
            Rj = C3(Rj)
            if i >4
                error("Problem in C3 rotation when finding pair")#
            end
        end
        Rk = C2(Rk) # global C2 maps B sublattice to A on lhs, so we now have a standard pair (0,j)
        Rj = C2(Rj)
        if Rk.b != 1
            error("Problem after C2 rotation")
        end
        Rj = translateToOrigin(Rj,Rk) #global translation to map Rk to first unit cell
        Rk = Rvec(0,0,0,1)
        return Rk, mapToSubsector(Rj) #give the vector in correct subsector
    end

    """
    returns Geometry struct with all relevant information about the Pyrochlore lattice. Considers only pairs withing distane R from origin.
    """
    function getPyrochlore(NLen,J = [1.,0.];test = false)
        Name = string("Pyrochlore_NLen=",NLen)
        System =  getLatticeGeometry(NLen,Name,pairToInequiv,inCorrectSubsector,Basis,test=test)
        @unpack PairList,couplings = System
        setNeighborCouplings!(couplings,J,PairList,Basis)

        return(System)
    end
end