@with_kw struct LatticeInfo{BasisType,RvecType,FunctionType}
    System::Geometry
    Basis::BasisType
    NLen::Int = System.NLen
    Npairs::Int = System.Npairs
    NUnique::Int = System.NUnique
    PairList::Vector{RvecType} = System.PairList
    PairTypes::Vector{sitePair} = System.PairTypes
    SiteList::Vector{RvecType} = unique(SpinFRGLattices.sortedPairList(NLen,Basis)[1])
    UnitCell::Vector{RvecType} = [SpinFRGLattices.getRvec(b,Basis) for b in Basis.b]
    pairToInequiv::FunctionType
end

##
function FourierTransform(k::StaticArray,Chi_R, Lattice)
    @unpack PairList,SiteList,PairTypes,Basis,UnitCell,pairToInequiv = Lattice
    Chi_k = 0. +0im
    for i_site in UnitCell
        Ri = getCartesian(i_site,Basis)
        # println(Ri)
        for j_site in SiteList # site summation
            Rj = getCartesian(j_site,Basis)
            Rij = Ri - Rj
            R_Ref,ij = pairToInequiv(i_site,j_site) #Map j to correct pair so that we may use Chi_0,j'
            xi = getSiteType(R_Ref,Basis)
            pair = MapToPair(xi,ij,PairList,PairTypes)
            if pair !== 0
                Chi_k += 1/Basis.NCell * exp(1im * k' * Rij) * Chi_R[pair]
            end
            # println(j_site,Chi_R[pair])
        end
    end
    return real(Chi_k)
end

"""Returns 2D Fourier trafo in plane as specified by the "regionfunc" function. Eg for a plot in the xy plane we can use plane = (ki,kj) -> SA[ki,kj] """
function Fourier2D(Chi_R::AbstractArray,regionfunc::Function,Lattice;res=100,ext = pi)
    karray = range(-ext,stop = ext,length = res)
    Chi_k = zeros(res,res)

    Threads.@threads for i in 1:res
        ki = karray[i]
        for (j,kj) in enumerate(karray)
            Chi_k[j,i] = FourierTransform(regionfunc(kj,ki),Chi_R,Lattice)
        end
    end
    return karray,Chi_k
end
##

"""Computes χ_ij(τ=0)"""
function equalTimeChi(Chi_RNu)
    Npairs,N_nu = size(Chi_RNu)
    Chi_Tau0 = Chi_RNu[:,begin] # add static nu=0 component, appears only once in sum
    for R in eachindex(Chi_Tau0)
        sum = 0.
        for n in 2:N_nu # dynamic components
            sum += Chi_RNu[R,n]
        end
        Chi_Tau0[R] += 2*sum# Chi(nu=0)
    end
    return Chi_Tau0
end

"""Compute energy divided by temperature from spin correlations"""
function EnergyBeta(Chi_RNu, Lattice)
    @unpack PairList,SiteList,PairTypes,Basis,UnitCell,pairToInequiv,Npairs = Lattice
    J_ij = Lattice.System.couplings
    E = 0.
    Chi_Tau0 = equalTimeChi(Chi_RNu)

    for i_site in UnitCell
        # println(Ri)
        for j_site in SiteList # site summation
            R_Ref,ij = pairToInequiv(i_site,j_site) #Map j to correct pair so that we may use Chi_0,j'
            xi = getSiteType(R_Ref,Basis)
            pair = MapToPair(xi,ij,PairList,PairTypes)
            if pair !== 0
                E += 3/(2*Basis.NCell) *J_ij[pair] * Chi_Tau0[pair]
            end
            # println(j_site,Chi_R[pair])
        end
    end
    return E
end

function get_e_Chi(Chi_TRnu,Lattice)
    e_Chi = similar(Results.T)
    for (iT,T) in enumerate(Results.T)
        e_Chi[iT] = @views T*EnergyBeta(Chi_TRnu[iT,:,:],Lattice)
    end
    return e_Chi
end


##
pitick(x) = latexstring("$(round(Int,x/pi)) \\pi") 
PiMultipleTicks(ticks) = [pitick(x) for x in ticks  ]

function Chikplot(k,Chi_k;xlabel = L"k_x",ylabel= L"k_y",colorscheme = :viridis,tickfontsize = 14,labelfontsize = 20,step = 1/2, kwargs...)
    min,max = minimum(k),maximum(k)
    ticks = collect(min:max*step:max)
    ticklabels = PiMultipleTicks(ticks)
    # pl = heatmap(k,k,transpose(Chi_k),size = (640, 600),xticks=(ticks,ticklabels),yticks=(ticks,ticklabels),linewidths=0.0,xlabel=xlabel ,ylabel= ylabel,c= colorscheme,right_margin = 15 *Plots.px,aspectratio = 1,tickfontsize = tickfontsize,labelfontsize = labelfontsize;kwargs...)
    pl = heatmap(k,k,transpose(Chi_k),size = (570, 600),xticks=(ticks,ticklabels),yticks=(ticks,ticklabels),linewidths=0.0,xlabel=xlabel ,ylabel= ylabel,c= colorscheme,right_margin = 15 *Plots.px,aspectratio = 1,tickfontsize = tickfontsize,labelfontsize = labelfontsize,ylims = [min,max];kwargs...)
    return pl
end

##
function getFlow(k::StaticArray,Chi_LR,Lambdas,Lattice)
    flow = similar(Lambdas)
    FT(k,Chi) = FourierTransform(k,Chi,Lattice)
    for i in eachindex(Lambdas)
        flow[i] =  @views FT(k,Chi_LR[i,:])
    end
    return Lambdas,flow
end

function plotFlow(k::StaticArray,Chi_LR,Lambdas,Lattice,pl = plot();method = plot!,xmax=1.,kwargs...)
    flow = getFlow(k,Chi_LR,Lambdas,Lattice)
    method(pl,flow, xlims = (0.,xmax);kwargs...)
    return pl
end 

function plotMaxFlow(Chi_LR,Lambdas,Lattice,regionfunc::Function,pl = plot();  res = 30,ext = pi,xmax=1.,method = plot!,kwargs...)
    flow = similar(Lambdas)
    FT(Chi) = maximum(Fourier2D(Chi,regionfunc,Lattice,res = res,ext=ext)[2])
    for i in eachindex(Lambdas)
        flow[i] =  @views FT(Chi_LR[i,:])
    end
    method(pl,Lambdas,flow, xlims = (0.,xmax);kwargs...)
    return pl
end 
##
function pointPath(p1::StaticArray,p2::StaticArray,res)
    Path = Vector{typeof(p1)}(undef,res)
    for i in eachindex(Path)
        Path[i] = p1 + i/res*(p2 -p1)
    end
    return Path
end
"""res contains the number of points along -pi,pi"""
function fetchKPath(points,res = 100)
    Path = Vector{typeof(points[begin])}(undef,0)
    # Path = []
    PointIndices = [1]
    for i in eachindex(points[begin:end-1])
        p1 = points[i]
        p2 = points[i+1]
        append!(Path,pointPath(p1,p2,round(Int,norm(p1-p2)/2pi * res)))
        append!(PointIndices,length(Path)) # get indices corresponding to points
    end
    return PointIndices,Path
end
##

function plotKpath(Chi_R,Lattice,points,pl = plot();  PointTicks = [round.(p,digits = 2) for p in points],res = 100,xmax=1.,method = plot!,kwargs...)

    # PointTicks(ticks) = [p for p in points ]
    
    FT(k) = FourierTransform(k,Chi_R,Lattice)
    ticks,Path = fetchKPath(points,res)
    FTPath = FT.(Path)
    method(pl,1:length(Path),FTPath,xticks = (ticks,PointTicks);kwargs...)
    return pl
end 
##
# pscatter!(pArray;kwargs...) = scatter!(Tuple(p for p in pArray);kwargs...)
# pscatter(pArray;kwargs...) = scatter(Tuple(p for p in pArray);kwargs...)

pscatter!(pArray;kwargs...) = scatter!(Tuple([p[i] for p in pArray] for i in 1:length(pArray[1]));kwargs...)
pscatter(pArray;kwargs...) = scatter(Tuple([p[i] for p in pArray] for i in 1:length(pArray[1]));kwargs...)
pplot!(pArray;kwargs...) = plot!(Tuple(p for p in pArray);kwargs...)
