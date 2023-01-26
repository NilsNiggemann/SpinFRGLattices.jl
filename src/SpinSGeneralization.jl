NumberOfSiteCopies(S::Real) = isinteger(2S) ? Int(2S) : error("S=$S ≠ n/2 is not a physical spin quantum number!")

createSatellite(R::Rvec_3D,nb::Integer) = Rvec(R.n1,R.n2,R.n3,nb)
createSatellite(R::Rvec_2D,nb::Integer) = Rvec(R.n1,R.n2,nb)

"""returns a list of indices that map pairs in PairList to pairs in PairList_new"""
function mapPairsToNewList(PairList,PairTypes,PairList_S,PairTypes_S)
    mapping = zeros(Int,length(PairList))
    for (i,(Rij,T)) in enumerate(zip(PairList,PairTypes))
        ijpr = MapToPair(T.xi,Rij,PairList_S,PairTypes_S)
        @assert ijpr != 0 "pair no found: xi = $xi , Rij = $Rij"
        mapping[i] = ijpr
    end
    return mapping
end
mapPairsToNewList(G1,G2) = mapPairsToNewList(G1.PairList,G1.PairTypes,G2.PairList,G2.PairTypes)

"""remove onsite pairs from site summation"""
function SpinSSiteSummation!(S,ijprimes,OnsitePairs,M)
    filter!(x-> !(x.ki in OnsitePairs || x.kj in OnsitePairs),S) # M Σ_(k ≠ i) V_k1i1 V_k1i1 P_kk
    S.ki .= ijprimes[S.ki]
    S.kj .= ijprimes[S.kj]
    unique!(S)
    S.m .*= M
    sort!(S)
end

"""Adapts vector containing onsite pairs by shifting the indices by one for each inequivalent satellite pair added"""
adaptOnsitePairs(OnsitePairs) = [P + (i-1) for (i,P) in enumerate(OnsitePairs)]

"""Returns vectors containing information about inequivalent pairs after inserting a new pair after each onsite pair, to keep physical properties of the system, while also including additional copies of spins on each site coupled by a couplings strenght A."""
function adaptPairs(G::Geometry,NCell::Int,A)
    (;PairList,PairTypes,OnsitePairs,invpairs,couplings) = G
    PairList_new = copy(PairList)
    PairTypes_new = copy(PairTypes)
    couplings_new = copy(couplings)
    OnsitePairs_new = adaptOnsitePairs(OnsitePairs)
    
    for (i,ij) in enumerate(OnsitePairs)
        Site = PairList[ij]
        ijpr = OnsitePairs_new[i]
        sat = createSatellite(Site,NCell+i) #for each onsite pair, we create a new satellite site
        insert!(PairList_new,ijpr+1,sat)
    end

    for (i,ij) in enumerate(OnsitePairs)
        pt = PairTypes[ij]
        ijpr = OnsitePairs_new[i]
        insert!(PairTypes_new,ijpr+1,pt)
        insert!(couplings_new,ijpr+1,A[i]) 
    end

    ijprimes = mapPairsToNewList(PairList,PairTypes,PairList_new,PairTypes_new)
    invpairs_new =  ijprimes[invpairs]
    for (i,ij) in enumerate(OnsitePairs_new)
        insert!(invpairs_new,ij+1,ij+1) # (i2,i1) = (i1,i2) from satellite permutation symmetry
    end

    return (;PairList = PairList_new,OnsitePairs = OnsitePairs_new,PairTypes = PairTypes_new,invpairs = invpairs_new,couplings = couplings_new,ijprimes)
end

adaptPairs(G::Geometry,NCell::Int,A::Real) = adaptPairs(G,NCell,fill(A,length(G.OnsitePairs)))

"""Takes a geometry and returns another geometry that is generalized for higher spin quantum number, effectively multiplying the number of sites by 2S. A is the (typically ferromagnetic) coupling by aligning the copied spins to maximize total spin. The physical Spin-S solution is recovered in the limit A → -∞. 
If there are several inequivalent sites, A may also be chosen as a vector with one component for each inequivalent site.
"""
function adaptForSpinS(G::Geometry,NCell::Integer,Spin::Real,A)
    M = NumberOfSiteCopies(Spin)
    M == 1 && return G
    if occursin("2S_",G.Name)
        @error "Invalid geometry name: '2S_' is already reserved for spin-S generalization. Please do not try to modify a geometry that was already adapted for higher spins."
    end
    (;PairList,OnsitePairs,PairTypes,invpairs,couplings,ijprimes) = adaptPairs(G,NCell,A)

    Npairs = G.Npairs + length(OnsitePairs) # one additional pair for each symmetry inequiv site
    NSumMax = size(G.siteSum,1)+2 # two terms are added
    
    zeroTerm = sumElements(0,0,0,0)

    siteSum = StructArray(fill(zeroTerm, NSumMax,Npairs))
    
    # ijprimes = mapPairsToNewList(G.PairList,G.PairTypes,PairList,PairTypes)

    for ij in eachindex(G.PairList)
        (;xi,xj) = G.PairTypes[ij]

        ijpr = ijprimes[ij]
            
        i1i1 = OnsitePairs[xi]
        i1i2 = i1i1 + 1 # PairList is arranged such that intersatellite terms are directly next to onsite bonds
        j1j1 = OnsitePairs[xj]
        j1j2 = j1j1 + 1

        i1j1 = ijpr
        j1i1 = invpairs[i1j1]

        S = G.siteSum[begin:G.Nsum[ij],ij] #allocate new vector instead of mutating geometry
        SpinSSiteSummation!(S,ijprimes,G.OnsitePairs,M)
        NumTerms = length(S)
        
        SiteSum_ij = @view siteSum[1:NumTerms,ijpr]
        S_sat_ij = @view siteSum[NumTerms+1:end,ijpr] # site sum for terms from satellite sites 

        SiteSum_ij .= S

        if ij in G.OnsitePairs
            
            # treat true onsite pair which also aquires a sum term for each of the other M-1 satellites 

            S_sat_ij[1] = sumElements(i1i1,i1i1,1,xi)
            S_sat_ij[2] = sumElements(i1i2,i1i2,M-1,xi) # (M-1) * V_i1i2 V_i1i2 P_ii 

            # treat inter-satellite pairs i1i2. These pairs are new so they are written in the next column
            @assert ijpr+1 == i1i2 "Error in spin S adaptiation: Assumption that intersatellite vertices i1i2 = $i1i2 are stored next to onsite vertices ii+1 = $(ijpr+1) is not fulfilled!" 
            
            SiteSum_ij = @view siteSum[1:NumTerms,ijpr+1]
            S_sat_ij = @view siteSum[NumTerms+1:end,ijpr+1] # site sum for terms from satellite sites 

            SiteSum_ij .= S

            S_sat_ij[1] = sumElements(i1i1,i1i2,1,xi) # V_i1i1 V_i1i2 P_ii (α = 1)
            S_sat_ij[2] = sumElements(i1i2,i1i1,1,xi) # V_i1i1 V_i1i2 P_ii (α = 2)

            M > 2 && (S_sat_ij[3] = sumElements(i1i2,i1i2,M-2,xi) ) # V_i1i1 V_i1i2 P_ii (α > 2)
        else

            S_sat_ij[1] = sumElements(i1i1,i1j1,1,xi)
            S_sat_ij[2] = sumElements(j1i1,j1j1,1,xj)

            S_sat_ij[3] = sumElements(i1i2,i1j1,M-1,xi)
            S_sat_ij[4] = sumElements(j1i1,j1j2,M-1,xj)
            
        end

    end
    for ij in axes(siteSum,2)
        S = @view siteSum[:,ij]
        sort!(S)
    end
    Name = G.Name*"_2S_$M"
    NLen = G.NLen
    NUnique = G.NUnique

    return Geometry(;Name,NLen,siteSum,Npairs,couplings,PairList,PairTypes,OnsitePairs,NUnique,invpairs)
end

"""Given a vector of susceptibilities for the spin-S generalized pairs and the indices of onsite pairs for the spin-1/2 model returns a vector of susceptibilities for the corresponding spin-S model."""
function convertSusceptibilityToSpinS!(Chinew::AbstractVector,Chi::AbstractVector,Spin::Real,OnsitePairs_S12::AbstractVector)
    M = NumberOfSiteCopies(Spin)
    OnsitePairs_S = adaptOnsitePairs(OnsitePairs_S12)
    intersitePairs = OnsitePairs_S .+1
    for i in eachindex(Chinew,Chi)
        if i in OnsitePairs_S
            Chinew[i] = M *(Chi[i] + (M-1) * Chi[i+1])
        elseif i in intersitePairs
            Chinew[i] = NaN
        else
            Chinew[i] = M^2 .*Chi[i]
        end
    end 
    deleteat!(Chinew,intersitePairs)
    return Chinew
end
convertSusceptibilityToSpinS!(Chi::AbstractVector,Spin::Real,OnsitePairs_S12::AbstractVector) = convertSusceptibilityToSpinS!(Chi,Chi,Spin,OnsitePairs_S12)
convertSusceptibilityToSpinS(Chi::AbstractVector,Spin::Real,OnsitePairs_S12::AbstractVector) = convertSusceptibilityToSpinS!(similar(Chi),Chi,Spin,OnsitePairs_S12)