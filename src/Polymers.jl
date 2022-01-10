"""Demonstration: returns polymer of Nsites equivalent sites. Useful for testing!"""
function getPolymer(Nsites;coup = [0,1,0.5])
    # Return geometry for a polymer of Nsites equivalent sites
    Npairs = Nsites÷2 +1 # includes onsite pair
    function dist(i,j,N)
        #returns distance between sites i and j. Also corresponds to pair type
        if (i < j) i,j = j,i end
        dist1(i,j) = abs(i- j)
        dist2 = dist1(i,N) + dist1(0,j) #accounts for periodic boundaries
        return(min(dist1(i,j),dist2))
    end
    pair(i,j,N) = dist(i,j,N) +1 # pair indexing starts at 1 for onsite pair

    PrepSitesum = zeros(Int,2,Nsites,Npairs) #saves summation elements pairs with  (b1,b2) for Npairs unique pairs
    for j in 1:Npairs
        for k in 1:Nsites
            ki = pair(1,k,Nsites)
            kj = pair(k,j,Nsites)
            PrepSitesum[:,k,j] = [ki,kj]
        end
    end
    #eliminate multiplicities
    siteSum = Matrix{sumElements}(undef,Nsites,Npairs)

    function numocc(obj,arr)
        m = 0
        for el in eachcol(arr)
            if (obj == el) m+=1 end
        end
        return m
    end

    for j in 1:Npairs
        for k in 1:Nsites
            VertexPair = PrepSitesum[:,k,j]
            m = numocc(VertexPair, PrepSitesum[:,:,j])
            
            xk = 1
            spl = sumElements(VertexPair...,m,xk)
            if spl in siteSum[:,j]
                spl = sumElements(0,0,0,0)
            end
            siteSum[k,j] = spl
        end
    end
    couplings = zeros(Npairs)
    if length(coup) >Npairs
        coup = coup[1:Npairs]
    end
    couplings[1:length(coup)] = coup
    Polymer = Geometry(Name = "$(Nsites)Polymer",couplings=couplings,PairList = collect(1:Npairs),siteSum=siteSum,NLen = Nsites)
    return Polymer
end