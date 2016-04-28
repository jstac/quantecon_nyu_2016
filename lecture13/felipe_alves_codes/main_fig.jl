using PyPlot
using QuantEcon: meshgrid



#== Value Function ==#
function fig_vf(fd::FDSpec)
    b = fd.b; a = fd.a; z = fd.z;

    an = length(a)
    bn = length(b)
    zn = length(z)
    V = vcat(fd.V...)
    V = reshape(V,bn,an,zn)
    p_args = Dict{Symbol,Float64}(:lw => 2, :alpha =>0.7)

    fig, ax = subplots()
    for ai in 1:7:an
        ax[:plot](fd.b[end-10:end], V[end-10:end,ai,1] ; p_args...)
    end
    # ax[:legend](loc= "upper right", fontsize = 11)
    ax[:set_xlabel]("Liquid asset")
    fig[:show]()

    fig, ax = subplots()
    for bi in 1:7:an
        ax[:plot](fd.a, squeeze(V[bi,:,1],1) ; p_args...)
    end
    # ax[:legend](loc= "upper right", fontsize = 11)
    ax[:set_xlabel]("Iliquid asset")
    fig[:show]()
end

"""
    Compute the figures...
"""
function fig_pol(fd::FDSpec)

    sol = fd.sol
    c = sol.c; d = sol.d; sc = sol.sc
    b = fd.b; a = fd.a; z = fd.z;

    an = length(a)
    bn = length(b)
    zn = length(z)

    netainc = fd.netainc

    p_args  = Dict{Symbol,Any}(:rstride=>2, :cstride=>2, :cmap=>ColorMap("Oranges"), :alpha=>0.7, :linewidth=>0.5)
    pp_args = Dict{Symbol,Float64}(:alpha=>0.7, :linewidth=>2)
    agrid, bgrid = meshgrid(a,b)

    #== Consumption 3D ==#
    cons_low = c[:,:,1]
    fig = figure(figsize = (8,6))
    ax = fig[:gca](projection="3d")
    ax[:plot_surface](agrid, bgrid, cons_low; p_args...)
    ax[:set_xlabel]("illiquid")
    ax[:set_ylabel]("liquid")
    ax[:legend]()
    ax[:set_title]("Consumption")
    fig[:show]()

    #== Consumption 2D ==#
    fig, ax = subplots()
    for ai in 1:7:an
        ax[:plot](fd.b, c[:,ai,1] ; pp_args...)
    end
    # ax[:legend](loc= "upper right", fontsize = 11)
    ax[:set_xlabel]("Liquid asset")
    ax[:set_title]("Consumption Low Shock ")
    fig[:show]()

    fig, ax = subplots()
    for ai in 1:7:an
        ax[:plot](fd.b, c[:,ai,2] ; pp_args...)
    end
    # ax[:legend](loc= "upper right", fontsize = 11)
    ax[:set_xlabel]("Liquid asset")
    ax[:set_title]("Consumption High Shock")
    fig[:show]()

    #== Deposit ==#
    fig = figure(figsize = (8,6))
    ax = fig[:gca](projection="3d")
    dep_high = d[:,:,end]
    ax[:plot_surface](agrid, bgrid, dep_high; p_args...)
    ax[:set_xlabel]("illiquid")
    ax[:set_ylabel]("liquid")
    ax[:legend]()
    ax[:set_title]("Deposit High")
    fig[:show]()

    #== Illiquid Savings ==#
    # adrift = zeros(d)
    # for zi in 1:zn, ai in 1:an, bi in 1:bn
    #     adrift[bi,ai,zi] = d[bi, ai, zi] + netainc[ai,zi]
    # end
    # fig = figure(figsize = (8,6))
    # ax = fig[:gca](projection="3d")
    # adrift_low = adrift[:,:,1]
    # ax[:plot_surface](agrid, bgrid, adrift_low; p_args...)
    # ax[:set_xlabel]("illiquid")
    # ax[:set_ylabel]("liquid")
    # ax[:legend]()
    # ax[:set_title]("Illiquid Saving")
    # fig[:show]()
end

function fig_dens(fd::FDSpec)

    a,b = fd.a, fd.b
    an = length(a)
    bn = length(b)
    zn = length(fd.z)
    gmat = reshape(vcat(fd.sol.g...),bn,an,zn)

    ΔTbgrid = fd.ΔTbgrid
    ΔTagrid = fd.ΔTagrid


    p_args = Dict{Symbol,Float64}(:lw => 2, :alpha =>0.7)
    #== Marginal wrt to a ==#
    adens = zeros(an,zn)
    for zi in 1:zn, ai in 1:an
        adens[ai,zi] = sum(gmat[:,ai,zi].*ΔTbgrid)
    end
    #= Check summation =#
    @printf("Check SUM   \n")
    @printf("Sum to %.2f \n",sum(adens.*repmat(fd.ΔTagrid,1,zn)))

    Ea = sum( a.* sum(adens,2).* ΔTagrid )
    fig, ax = subplots()
    for zi in 1:zn
        ax[:plot](fd.a, adens[:,zi], label="Labor shock $zi";p_args...)
    end
    ax[:legend](loc= "upper right", fontsize = 11)
    ax[:set_xlim]([0.0,2*Ea])
    ax[:set_title]("Density Illiquid")
    ax[:set_xlabel]("Illiquid asset")
    fig[:show]()

    #== Marginal wrt to b ==#
    bdens = zeros(bn,zn)
    for zi in 1:zn, bi in 1:bn
        bdens[bi,zi] = sum(squeeze(gmat[bi,:,zi],1).*ΔTagrid)
    end
    #= Check summation =#
    @printf("Check SUM   \n")
    @printf("Sum to %.2f \n",sum(bdens.*repmat(fd.ΔTbgrid,1,zn)))

    fig, ax = subplots()
    for zi in 1:zn
        ax[:plot](b, bdens[:,zi], label="Labor shock $zi";p_args...)
    end
    ax[:set_xlabel]("Liquid asset")
    ax[:set_ylim]([0,0.8*maximum(bdens[:,1])])
    ax[:set_xlim]([b[1],5.0])
    ax[:set_title]("Density Illiquid", fontsize = 20)
    ax[:legend](loc= "upper right", fontsize = 11)
    fig[:show]()

    bran_max = searchsortedfirst(b,5.0)
    # aran_max = searchsortedfirst(a,Ea)
    aran_max = searchsortedfirst(a,25.0)
    #== Joint density ==#
    arange, brange = 1:min(round(Int,aran_max),an), 1:bran_max
    # dens_low  = gmat[brange,arange,1]
    # dens_high = gmat[brange,arange,2]

    #= Construct probabilities =#
    gvec = vcat(gmat)

    ΔTab = zeros(gvec)
    bvec = zeros(gvec)
    avec = zeros(gvec)
    ijk = zero(Int)

    for zi = 1:zn, abi = 1:bn*an
        ijk += 1
        ΔTab[ijk] = ΔTagrid[afromab(abi)] * ΔTbgrid[bfromab(abi)]
        bvec[ijk] = b[bfromab(abi)]
        avec[ijk] = a[afromab(abi)]
    end

    prob = gvec .* ΔTab
    ggmat = reshape(prob,bn,an,zn)

    dens  = sum(ggmat[brange,arange,:],3)
    dens  = squeeze(dens,3)
    dens2 = ggmat[brange,arange,2]

    #== Construct 3D figure ==#
    fig = figure(figsize = (8,6))
    ax = fig[:gca](projection="3d")
    agrid, bgrid = meshgrid(a[arange],b[brange])
    ax[:plot_surface](agrid, bgrid, dens,
                    rstride=2,cstride=2, cmap=ColorMap("jet"), alpha=0.5, linewidth=0.5)
    ax[:set_xlabel]("Illiquid")
    ax[:set_ylabel]("Liquid")
    # ax[:set_xlim]([0,)
    # ax[:set_ylim]([-2,10])
    ax[:legend]()
    fig[:show]()

    fig, ax = subplots(figsize = (8,6))
    ax[:xaxis][:grid](true, zorder=0)
    ax[:yaxis][:grid](true, zorder=0)
    ax[:contourf](agrid, bgrid, dens, 5, alpha=0.6, cmap=ColorMap("jet"))
    cs1 = ax[:contour](agrid, bgrid, dens, 5, colors="black",lw=2)
    ax[:clabel](cs1, inline=1, fontsize=10)
    ax[:set_xlabel]("Illiquid")
    ax[:set_ylabel]("Liquid")
    # ax[:set_xlim]([0,)
    # ax[:set_ylim]([-2,10])
    ax[:legend]()
    fig[:show]()

# [0,0.005,0.01,0.02,0.04]

end

# function fig_pol(fd::TwoAssetsFD)
#
#     c = fd.c; d = fd.d; sc = fd.sc
#     b = fd.b; a = fd.a; z = fd.z;
#
#     an = length(a)
#     bn = length(b)
#     zn = length(z)
#
#     p_args = Dict{Symbol,Any}(:rstride=>2, :cstride=>2, :cmap=>ColorMap("Oranges"), :alpha=>0.7, :linewidth=>0.5)
#     agrid, bgrid = meshgrid(a,b)
#
#     #== Consumption ==#
#     cons_low = c[:,:,1]
#     fig = figure(figsize = (8,6))
#     ax = fig[:gca](projection="3d")
#     ax[:plot_surface](agrid, bgrid, cons_low; p_args...)
#     ax[:set_xlabel]("illiquid")
#     ax[:set_ylabel]("liquid")
#     ax[:legend]()
#     fig[:show]()
#
#     #== Deposit ==#
#     fig = figure(figsize = (8,6))
#     ax = fig[:gca](projection="3d")
#     dep_low = d[:,:,1]
#     ax[:plot_surface](agrid, bgrid, dep_low; p_args...)
#     ax[:set_xlabel]("illiquid")
#     ax[:set_ylabel]("liquid")
#     ax[:legend]()
#     ax[:set_title]("Deposit")
#     fig[:show]()
# end
