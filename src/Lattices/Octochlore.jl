using SpinFRGLattices
using FRGLatticePlotting
import SpinFRGLattices as SL
using StaticArrays,Parameters
using StructArrays,Symbolics,Plots
##

function OctochloreBasis()
    a1 = SA[1,0,0]
    a2 = SA[0,1,0]
    a3 = SA[0,0,1]

    b1 = 1/2*a1
    b2 = 1/2*a2
    b3 = 1/2*a3
    
    return Basis_Struct_3D(a1=a1,a2=a2,a3=a3,b=[b1,b2,b3],NNdist = norm(b1-b2))
end
const Basis = OctochloreBasis()

import SpinFRGLattices: getCartesian,dist,norm,getRvec
getCartesian(R) = SpinFRGLattices.getCartesian(R,Basis) 
dist(R1,R2) = SpinFRGLattices.dist(R1,R2,Basis) 
norm(R) = SpinFRGLattices.norm(R,Basis) 
getRvec(R)::Rvec_3D = SpinFRGLattices.getRvec(R,Basis) 
getRefCartesian(R) = getCartesian(R) - SA[0.5,0,0]
##
R(x::Integer) = Rvec(0,0,0,x)
plotly()
##
"""mirrors x-coordinate at x = 0.5 """
xMirror(r::AbstractVector) = SA[1-r[1],r[2],r[3]]
xMirror(R::Rvec) = getRvec(xMirror(getCartesian(R)))
"""mirrors y-coordinate at y = 0 """
yMirror(r::AbstractVector) = SA[r[1],-r[2],r[3]]
yMirror(R::Rvec) = getRvec(yMirror(getCartesian(R)))
"""mirrors z-coordinate at z = 0 """
zMirror(r::AbstractVector) = SA[r[1],r[2],-r[3]]
zMirror(R::Rvec) = getRvec(zMirror(getCartesian(R)))

"""C3 rotation around (1,1,1)"""
C3(r::AbstractVector) = SA[r[3],r[1],r[2]]
C3(R::Rvec) = getRvec(C3(getCartesian(R)))
"""C4 rotation around (1,0,0)"""
C4(r::AbstractVector) = SA[r[1],-r[3],r[2]]
C4(R::Rvec) = getRvec(C4(getCartesian(R)))
##

"""
Symmetry: may restrict to section with z<=y<=x
"""
function inCorrectSubsector(R_ref,R::Rvec)
    # @unpack n1,n2,n3 = R
    x,y,z = getCartesian(R)
    return(x>= 0 && y>=0 && z >=0&& y>= z ) # todo basis!
end

inCorrectSubsector(R::Rvec) = inCorrectSubsector(Basis.refSites[1],R)

"""gives coordinates of R' after using symmetry transformations on R.
R is arbitrary lattice vector"""
function mapToSubsector(R::Rvec)
    Rold = R
    i = 0
    x,y,z = getCartesian(R)
    i = 0
    while abs(y) < abs(z)
        R = C4(R)
        x,y,z = getCartesian(R)
        i+=1
        i >=4 && error("C4 Rotation failure $Rold")
    end
    if x <0
        R= xMirror(R)
    end
    
    if y <0
        R= yMirror(R)
    end
    
    if z <0
        R= zMirror(R)
    end
    
    @assert inCorrectSubsector(R) "Error: Vector not correctly mapped: $Rold -> $R"
    return R
end

function pairToInequiv(Rk,Rj)
    i = 0
    while Rk.b != 1
        i+=1
        Rk = C3(Rk)
        Rj = C3(Rj)
        if i >=3
            error("Problem in C3 rotation when finding pair")#
        end
    end
    @assert Rk.b == 1 "Problem after C2 rotation"
    Rj = translateToOrigin(Rj,Rk) #global translation to map Rk to first unit cell
    Rk = Rvec(0,0,0,1)
    return Rk, mapToSubsector(Rj) #give the vector in correct subsector
end

## determine couplings from alpha beta gamma model
struct spin{N,R}
    fac::N
    site::R
end

##

function getOcts(refSite::Rvec,a = 1.,b=0.1,c = 0.01)
    Allsites = generateLUnitCells(3,Basis,refSite)
    # return pairsPlot(Allsites,Basis)
    center = getCartesian(refSite) -SA[0.5,0,0]
    centerNorm(x) = norm(getCartesian(x) - center)
    AllNorms = sort(unique(centerNorm.(Allsites)))
    getNeighorFromCenter(i) = filter(x->centerNorm(x) ≈ AllNorms[i],Allsites)
    firstOct = getNeighorFromCenter(1)
    secondOct = getNeighorFromCenter(2)
    thirdOct = filter(R-> sum( (R.n1-refSite.n1,R.n2-refSite.n2,R.n3-refSite.n3) .== 0) ==2 ,getNeighorFromCenter(3))
    return (StructArray(append!([spin(a,site) for site in firstOct],[spin(b,site) for site in secondOct],[spin(c,site) for site in thirdOct])))
end

function getCouplingsToS1(a = 1.,b=0.1,c = 0.01)
    L=2
    S1Terms = spin[]
    for n1 in -L:L,n2 in -L:L,n3 in -L:L
        octs = getOcts(Rvec(n1,n2,n3,1),a,b,c)
        for s1 in octs
            for s2 in octs
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

function reduceCouplings(CoupList)
    reducedList =  StructArray(spin[])
    for S in CoupList
        position = findall(x-> x==S.site,reducedList.site)
        if isempty(position)
            push!(reducedList,S)
        else
            index = only(position)
            reducedList.fac[index] += S.fac
        end
    end
    return reducedList
end
function getInequivCouplings(CoupList)
    reducedList =  StructArray(spin[])
    for S in CoupList
        S_red = mapToSubsector(S.site)
        position = findall(x-> x==S_red,reducedList.site)
        if isempty(position)
            push!(reducedList,spin(S.fac,S_red))
        else
            fac =  only(reducedList[position].fac)
            @assert S.fac - fac ≈ 0 "Coupling $(S.fac) != $fac does not match symmetry"
        end
    end
    return reducedList
end
@variables α,β,γ
AllSites = generatePairSites(4,Basis)
mapToSubsector(Rvec_3D(1, 0, 0, 3))
mapToSubsector.(AllSites)
redCoup = reduceCouplings(getCouplingsToS1(α,β,γ)) 
# redCoup = reduceCouplings(getCouplingsToS1(1.,0.5,0.1)) 
ineqCoup = getInequivCouplings(redCoup)
##
function getOctochlore(NLen,beta = 0.5,gamma = 0.1;test = false)
    Name = string("Octochlore_NLen=",NLen)
    System =  getLatticeGeometry(NLen,Name,pairToInequiv,inCorrectSubsector,Basis,test=test)
    @unpack PairList,couplings = System
    # setNeighborCouplings!(couplings,J,PairList,Basis)

    return(System)
end
System = getOctochlore(4,test = true)
# scatter!([0.5],[0.5],[0.5],color = :black,label = "center")

