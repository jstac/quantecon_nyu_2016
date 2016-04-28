import Distances: chebyshev

##############################################################################
##
## Types and constructors
##
##############################################################################


##################################
##
## Household PROBLEM
##
##################################
"""
Consumption-saving problem with two assets and non-convex adjustment costs

## Fields

#### Parameters
- `γ::Float64`  : parameter on CRRA utility
- `ρ::Float64`  : discount factor
- `ξ::Float64`  : automatic deposit
- `σ::Float64`  : frisch elasticity of labor supply
- `ψ::Float64`  : disutility of labor

#### Prices
- `rᴬ::Float64`     : interest rate on illiquid asset
- `rᴮ::Float64`     : interest rate for savings in liquid asset
- `wedge::Float64`  : interest rate for borrowing in liquid asset
- `w::Float64`      : wage rate

#### Tax
- `τ::Float64`      : tax on income
- `T::Float64`      : Lump-sum transfer

#### Cost Function
- `χ₀::Float64` :
- `χ₁::Float64` :
"""
type TwoAssetsProb
    γ::Float64      # CRRA utility with parameter γ
    ρ::Float64      # discount rate
    ξ::Float64      # automatic deposit on illiquid asset
    σ::Float64
    ψ::Float64

    rᴬ ::Float64    # ret on illiquid asset
    rᴮ::Float64    # ret on liq asset
    wedge::Float64
    w  ::Float64    # wage rate

    τ::Float64
    T::Float64

    χ₀::Float64     # parameters on adjustment cost
    χ₁::Float64
end

function Base.show(io::IO, twoap::TwoAssetsProb)
    @printf io "    HOUSEHOLD PROBLEM\n"
    @printf io "\n"
    @printf io "    Parameters  \n"
    @printf io "   ------------ \n"
    @printf io "     γ:      %.3f\n" twoap.γ
    @printf io "     ρ:      %.3f\n" twoap.ρ
    @printf io "     ξ:      %.3f\n" twoap.ξ
    @printf io "\n"
    @printf io "    Prices  \n"
    @printf io "   -------- \n"
    @printf io "     rᴬ   :    %.3f \n" twoap.rᴬ
    @printf io "     rᴮ   :    %.3f \n" twoap.rᴮ
    @printf io "     wage :    %.3f \n" twoap.w
    @printf io "     τ    :    %.3f \n" twoap.τ
    @printf io "     T    :    %.3f \n" twoap.T
    @printf io "\n"
    @printf io "    Cost Function  \n"
    @printf io "   --------------- \n"

    m ="""
    χ₀ + χ₁ * (d/a)² * a
    """
    @printf io "     %s\n" m
end

function TwoAssetsProb2(pc::Prices; γ::Float64 = 2.0, ρ::Float64 = 0.06, ξ::Float64 = 0.1,
                            σ = 0.5, ψ = 27.0,
                            τ::Float64 = 0.0, T::Float64 = 0.,
                            χ₀::Float64 = 0.08, χ₁::Float64 = 3.0)

    # REVIEW: define some global functions that will now change on the go
    #== Substitute the functions to be used afterward ==#
    # global χ(d,a) = χ(d,a, [χ₀,χ₁])
    # global dχinv(pVb::Float64, pVa::Float64, a::Float64) = dχinv(pVb, pVa, a, χ₁)
    # global

    TwoAssetsProb(γ, ρ, ξ, σ, ψ, pc.rᴬ, pc.rᴮ, pc.wedge, pc.w, τ, T, χ₀, χ₁)
end

_unpackparams(twoap::TwoAssetsProb) =
twoap.γ, twoap.ρ, twoap.ξ, twoap.σ, twoap.ψ
_unpackprices(twoap::TwoAssetsProb) =
twoap.rᴬ, twoap.rᴮ, twoap.rᴮ+twoap.wedge, twoap.w
_unpacktax(twoap::TwoAssetsProb) =
twoap.τ, twoap.T

function utilfn(twoap::TwoAssetsProb, c::Float64, ℓutil::Float64)
    return (c-ℓutil)^(1.- twoap.γ)/(1.-twoap.γ)
end

function χ(twoap::TwoAssetsProb, d::Float64, a::Float64, abar::Float64 = 2.0)
    χ₀, χ₁ = twoap.χ₀, twoap.χ₁
    return χ₀*(abs(d)>0) + 0.5*χ₁*d^2.0/max(a,abar)
end
# χ(d::Float64,a::Float64, χval::Array{Float64,1}) = χval[1]*(abs(d)>0) + 0.5*χval[2]*d^2.0/max(a,1e-5)

function dχinv(twoap::TwoAssetsProb, ∂ratio::Float64, a::Float64, abar::Float64 = 2.0)
    χ₁ = twoap.χ₁
    return 1.0/χ₁ * ( ∂ratio-1.0 ) * max(a,abar)
end

function change_prices!(twoap::TwoAssetsProb, pc::Prices)

    twoap.rᴬ, twoap.rᴮ, twoap.w = pc.rᴬ, pc.rᴮ, pc.w
    return Void
end

function change_transfer!(twoap::TwoAssetsProb, τ::Float64, T::Float64)

    twoap.τ, twoap.T = τ, T
    return Void
end

# OLD function
# function foc_dep(∂Vb::Float64, ∂Va::Float64, a::Float64, χ₁::Float64)
#
#     d = ( ∂Va/∂Vb-1 ) * a/χ₁              # optimal d in case of adjustment
#
#     # NOTE:340 compute the change in utility and use that to decide whether to act or
#     #       not...
#
#     ∂V = ∂Va*d - ∂Vb*( d + χ(d,a)  )      # change in utility
#     d = d * (∂V > 0)
# end

##############################################################################
##
## Helping functions
##
##############################################################################

function powerspacegrid(init::Real, eend::Real, n::Int64, k::Real = 1)

    if n<=2
        println("n has to be larger than 2")
        return
    end

    x = collect(linspace(0,1,n))

    z = x.^(1.0/float(k))

    return init +(eend - init)*z
end




##################################
##
## FDSpec
##
##################################

abstract FDSol
abstract FDSpec

type SolutionExp <: FDSol

    ## Solution
    V::Vector{Vector{Float64}}      # Value function
    g::Vector{Vector{Float64}}      # Density over state space
    c::Array{Float64,3}             # Optimal consumption
    sc::Array{Float64,3}            # Savings without deposit
    d::Array{Float64,3}             # Optimal deposit flow

    ## Matrices to be filled
    A ::Vector{Base.SparseMatrix.SparseMatrixCSC{Float64, Int}}  # used on HJB to compute vⁿ⁺¹
    Au::Vector{Base.SparseMatrix.SparseMatrixCSC{Float64, Int}}  # used for KFE to compute transition dynamics

end

unpackcopy_fdsol(fdsol::SolutionExp) = deepcopy(fdsol.V), deepcopy(fdsol.g), deepcopy(fdsol.c), deepcopy(fdsol.sc), deepcopy(fdsol.d), deepcopy(fdsol.A), deepcopy(fdsol.Au)

"""
Specification for the Implicit-Explicit Finite Difference method
"""
type FDExp <: FDSpec

    ## GRID INFO
    b      ::Vector{Float64}
    Δbgrid ::Vector{Float64}
    ΔTbgrid::Vector{Float64}
    rbdrift::Vector{Float64}
    netbinc::Matrix{Float64}

    a      ::Vector{Float64}
    Δagrid ::Vector{Float64}
    ΔTagrid::Vector{Float64}
    radrift::Vector{Float64}
    netainc::Matrix{Float64}

    ΔTab::Vector{Float64}

    ## LaborSupply decisions
    ℓsupply::Vector{Float64}
    ℓutilgrid::Vector{Float64}

    ## Δ in finite difference scheme
    invΔ::Float64
    invΔᴷ::Float64

    ## Stochastic Information
    z::Vector{Float64}
    λ::Matrix{Float64}
    λdiag::Vector{Float64}
    λoff::Matrix{Float64}

    ## Solution
    sol::SolutionExp

    ## Storage
    b̃ ::Vector{Vector{Float64}}                                  # RHS of hjb
    B ::Vector{Base.SparseMatrix.SparseMatrixCSC{Float64, Int}}  # storage matrix
end

function Base.show(io::IO, fd::FDExp)
    fd.Δbgrid[2]-fd.Δbgrid[1] == 0 ? (binfo = "uniform ") : (binfo = "non-uniform ")
    fd.Δagrid[2]-fd.Δagrid[1] == 0 ? (ainfo = "uniform ") : (ainfo = "non-uniform ")
    @printf io "\n"
    @printf io "    Explicit-Implicit Finite Difference Method\n"
    @printf io "\n"
    @printf io "    Grids  \n"
    @printf io "   ------- \n"
    @printf io "     %12sb: %3d points in [% .0f, %.0f]\n" binfo length(fd.b) fd.b[1] fd.b[end]
    @printf io "     %12sa: %3d poitns in [% .0f, %.0f]\n" ainfo length(fd.a) fd.a[1] fd.a[end]
end

function FDExp(twoap::TwoAssetsProb;
z::Vector{Float64} = [.5,1.5], λ::Matrix{Float64} = [-1/3 1/3; 1/3 -1/3],
fixedΔb = true, bn::Int =100, bmin::Float64 =-1.0, bmax::Float64 = 40.0, bparam::Vector{Float64} = [0.35, 0.7],
fixedΔa = true, an::Int = 50, amin::Float64 = 1e-8, amax::Float64 = 70.0, aparam::Float64 = 0.7,
invΔ::Float64 = .025, invΔᴷ::Float64 = 0.025)
    # z = [.8,1.3];
    # λ = [-1/3 1/3;1/3 -1/3];
    #
    # fixedΔb = true;
    # bn =100;
    # bmin =-2.0;
    # bmax = 40.0;
    # bparam = [0.35;0.7];
    #
    # fixedΔa = true;
    # an = 50;
    # amin = 0.0;
    # amax = 70.0;
    # aparam = 0.7;
    #
    # invΔ = .025;
    # invΔᴷ = 0.075

    #== Household parameters and prices ==#
    γ, ρ, ξ, σ, ψ  = _unpackparams(twoap)
    rᴬ, rᴮ, rᴮ⁻, w = _unpackprices(twoap)
    τ, T           = _unpacktax(twoap)

    #== Construct the grids ==#
    if fixedΔb
        b   = powerspacegrid(bmin,bmax,bn)
    else
        ngpbPOS = Int(0.9*bn)
        ngpbNEG = bn - ngpbPOS

        bpos = powerspacegrid(1e-8, bmax, ngpbPOS, bparam[2])
        bneg = zeros(ngpbNEG)
        n = Int(ngpbNEG/2+1)
        bneg[1:n] = powerspacegrid(bmin, (bpos[1]+bmin)/2, n , bparam[1])
        b = [bneg;bpos]
        for i in n+1:ngpbNEG
            b[i] = bpos[1] - (b[ngpbNEG+2-i]-b[1])
        end
    end

    fixedΔa ?  (a = powerspacegrid(amin,amax,an)) : (a = powerspacegrid(amin,amax,an,aparam[1]))

    zn  = length(z); abn = an*bn

    #== Grid counting functions ==#
    global afromab(abi) = div(abi-1,bn)+1
    global bfromab(abi) = rem(abi-1,bn)+1
    global abfromab(ai,bi) = (ai-1)*bn + bi
    global abzfromabz(ai,bi,zi) = (zi-1)*an*bn + ( ai-1 )*bn +bi

    #== Grid info ==#
    Δbgrid = b[2:end]-b[1:end-1]
    Δagrid = a[2:end]-a[1:end-1]

    ΔTbgrid = zeros(bn); ΔTagrid = zeros(an)
    ΔTbgrid[1]       = 0.5*Δbgrid[1]
    ΔTbgrid[2:end-1] = 0.5*( Δbgrid[2:end] + Δbgrid[1:end-1] )
    ΔTbgrid[end]     = 0.5*Δbgrid[end]

    ΔTagrid[1]       = 0.5*Δagrid[1]
    ΔTagrid[2:end-1] = 0.5*( Δagrid[2:end] + Δagrid[1:end-1] )
    ΔTagrid[end]     = 0.5*Δagrid[end]

    #== Construct a vector nb*na with the terms Δb̃[bi] * Δã[ai] =#
    ΔTab = zeros(abn)
    for abi = 1:abn

        ΔTab[abi] = ΔTagrid[afromab(abi)] * ΔTbgrid[bfromab(abi)]
    end

    #== Labor grid info ==#
    ℓsupply    = zeros(zn)
    ℓsupply[:] = (1/ψ * (1-τ) * w)^σ
    ℓutilgrid   = ψ * (z .* ℓsupply.^( 1 + 1.0/σ ) ) / ( 1+1./σ )

    #== return drifts ==#
    rbdrift = (b.>= 0.0)*rᴮ + (b.<0.0)*rᴮ⁻
    rbdrift = b .* rbdrift

    # NOTE: trick on the upper part of grid
    #== impose tax on the very top part of grid ==#
    τa = 15
    τc = rᴬ * ( a[end] * 0.999 )^(1 - τa)
    radrift = rᴬ * a - τc * a .^ τa

    #== net?drifts ==#
    netbinc = zeros(bn,zn)
    netainc = zeros(an,zn)

    for zi in 1:zn
        netℓinc = (1-τ)* w * z[zi] * ℓsupply[zi] + T
        netbinc[:,zi] =  (1-ξ) * netℓinc + rbdrift
        netainc[:,zi] = radrift + ξ * netℓinc
    end

    #== Idiosyncratic shock ==#
    λdiag = diag(λ)
    λoff  = λ - diagm(λdiag)

    ### STORAGE ###
    A = SparseMatrixCSC{Float64, Int}[spdiagm(
        (ones(abn), ones(abn-1), ones(abn-1), ones( abn -bn ) , ones( abn -bn ) ),
        (0, 1, -1, bn, -bn) ) for zi =1:zn]

    # B = SparseMatrixCSC{Float64, Int}[spdiagm(
    #     (ones(abn), ones(abn-1), ones(abn-1), ones( abn -bn ) , ones( abn -bn ) ),
    #     (0, 1, -1, bn, -bn) ) for zi =1:zn]

    B  = Array(SparseMatrixCSC{Float64, Int},zn)
    Au = Array(SparseMatrixCSC{Float64, Int},zn)
    for zi=1:zn
        fill!(nonzeros(A[zi]), zero(Float64))
        B[zi]  = deepcopy(A[zi])
        Au[zi] = deepcopy(A[zi])
    end

    #== Initial Distribution ==#
    g      = zeros(abn,zn)
    bpos_ind = findfirst(b.>0)

    #REVIEW:70 change ydist assumption
    ydist = ones(zn)./zn
    for zi=1:zn
        for abi =1:abn
            ai = afromab(abi)
            bi = bfromab(abi)
            ai==1 && bi==bpos_ind && ( g[abi,zi] = ydist[zi]/(ΔTagrid[ai]*ΔTbgrid[bi]) )
        end
    end

    g  = Vector{Float64}[g[:,zi] for zi=1:zn]
    #----------------------------------------------------------------

    ### SOLUTION ###
    V = Vector{Float64}[Array(Float64, abn) for zi=1:zn]
    ij = zero(Int)
    for zi in 1:zn
        ij = Int(0)
        for ai in 1:an, bi in 1:bn
            ij += 1
            V[zi][ij] = ( (1- ξ)*w*z[zi] + rᴬ*a[ai] + rᴮ⁻*b[bi] ).^( 1- γ )/( 1- γ )/ρ
        end
    end


    # == storage arrays == #
    b̃  = Vector{Float64}[Array(Float64, abn) for zi=1:zn]         # RHS of the Bellman equation

    #== policies ==#
    c  = Array(Float64, bn, an, zn)     # consumption policy
    sc = copy(c); d  = copy(c)

    #== Create solution type ==#
    sol = SolutionExp(V, g, c, sc, d, A, Au)

FDExp(
b, Δbgrid, ΔTbgrid, rbdrift, netbinc,
a, Δagrid, ΔTagrid, radrift, netainc, ΔTab,
ℓsupply, ℓutilgrid,
invΔ, invΔᴷ,
z, λ, λdiag, λoff,
sol,
b̃, B)
end

function initialize_solexp(twoap::TwoAssetsProb, fd::FDSpec)

    #= precompute the lengths =#
    b = fd.b
    a = fd.a
    z = fd.z
    bn = length(b); an = length(a); zn = length(z); abn = bn*an

    #= extract some other from fd =#
    ΔTagrid = fd.ΔTagrid
    ΔTbgrid = fd.ΔTbgrid
    netbinc = fd.netbinc



    #== Create the A matrix ==#
    A = SparseMatrixCSC{Float64, Int}[spdiagm(
        (ones(abn), ones(abn-1), ones(abn-1), ones( abn -bn ) , ones( abn -bn ) ),
        (0, 1, -1, bn, -bn) ) for zi =1:zn]

    Au = Array(SparseMatrixCSC{Float64, Int},zn)
    for zi=1:zn
        fill!(nonzeros(A[zi]), zero(Float64))
        Au[zi] = deepcopy(A[zi])
    end

    #== Initial Distribution ==#
    g      = zeros(abn,zn)
    bpos_ind = findfirst(b.>0)

    #REVIEW:70 change ydist assumption
    ydist = ones(zn)./zn
    for zi=1:zn
        for abi =1:abn
            ai = afromab(abi)
            bi = bfromab(abi)
            ai==1 && bi==bpos_ind && ( g[abi,zi] = ydist[zi]/(ΔTagrid[ai]*ΔTbgrid[bi]) )
        end
    end

    g  = Vector{Float64}[g[:,zi] for zi=1:zn]
    #----------------------------------------------------------------

    ### SOLUTION ###
    V = Vector{Float64}[Array(Float64, abn) for zi=1:zn]
    ij = zero(Int)
    for zi in 1:zn
        ij = Int(0)
        for ai in 1:an, bi in 1:bn
            ij += 1
            V[zi][ij] = ( netbinc[bi,zi] ).^( 1- γ )/( 1- γ )/ρ
        end
    end

    #== policies ==#
    c  = Array(Float64, bn, an, zn)     # consumption policy
    sc = copy(c); d  = copy(c)

    #== Create solution type ==#
    SolutionExp(V, g, c, sc, d, A, Au)
end

"""
Update the grid for a different specification of HOUSEHOLD PROBLEM
"""
function change_grid!(fd::FDSpec, twoap::TwoAssetsProb)

    #== Parameters ==#
    _, _, ξ, σ, ψ  = _unpackparams(twoap)
    rᴬ, rᴮ, rᴮ⁻, w = _unpackprices(twoap)
    τ, T           = _unpacktax(twoap)

    b, a, z = fd.b, fd.a, fd.z
    bn = length(b)
    an = length(a)
    abn = bn * an
    zn = length(z)

    #== Labor grid info ==#
    ℓsupply    = zeros(zn)
    ℓsupply[:] = (1/ψ * (1-τ) * w)^σ
    fd.ℓsupply = ℓsupply
    fd.ℓutilgrid   = ψ * (z .* ℓsupply.^( 1 + 1.0/σ ) ) / ( 1+1./σ )

    #== Update the drifts ==#
    rbdrift    = ( b.>= 0.0)*rᴮ + (b.<0.0)*rᴮ⁻
    fd.rbdrift = b .* rbdrift

    τa = 15
    τc = rᴬ * ( a[end] * 0.999 )^(1 - τa)
    fd.radrift = rᴬ * a - τc * a .^ τa

    for zi in 1:zn
        netℓinc = (1-τ)* w *z[zi] * ℓsupply[zi] + T
        fd.netbinc[:,zi] = (1-ξ) * netℓinc + fd.rbdrift
        fd.netainc[:,zi] = fd.radrift + ξ * netℓinc
    end

end



##############################################################################
##
## Implicit scheme
##
##############################################################################

type SolutionImp <: FDSol

    ## Household Solution
    V::Vector{Float64}          # Value function
    g::Vector{Float64}          # Density over state space
    c::Array{Float64,3}         # Optimal consumption
    sc::Array{Float64,3}        # Savings without deposit
    d::Array{Float64,3}         # Optimal deposit flow

    ## matrices to be filled
    A  ::Base.SparseMatrix.SparseMatrixCSC{Float64, Int}  # intensity matrix A on vⁿ⁺¹
    Au ::Base.SparseMatrix.SparseMatrixCSC{Float64, Int}  # intensity matrix A for kfe

end
"""
Specification for the Implicit Finite Difference method
"""
type FDImp <: FDSpec

    ## GRID INFO
    b      ::Vector{Float64}
    Δbgrid ::Vector{Float64}
    ΔTbgrid::Vector{Float64}
    rbdrift::Vector{Float64}
    netbinc::Matrix{Float64}

    a      ::Vector{Float64}
    Δagrid ::Vector{Float64}
    ΔTagrid::Vector{Float64}
    radrift::Vector{Float64}
    netainc::Matrix{Float64}

    ## LaborSupply decisions
    ℓsupply  ::Vector{Float64}
    ℓutilgrid::Vector{Float64}

    ## Δ in finite difference scheme
    invΔ ::Float64
    invΔᴷ::Float64

    ## Stochastic Information
    z::Vector{Float64}
    λ::Matrix{Float64}

    ## OUTPUT
    sol::SolutionImp

    ## Storage
    b̃  ::Vector{Float64}                                  # RHS of hjb
    B  ::Base.SparseMatrix.SparseMatrixCSC{Float64, Int}  # B storage matrix
    Λ  ::Base.SparseMatrix.SparseMatrixCSC{Float64, Int}  # stochastic process terms (useful bc can fill before)
end

function Base.show(io::IO, fd::FDImp)
    fd.Δbgrid[2]-fd.Δbgrid[1] == 0 ? (binfo = "uniform ") : (binfo = "non-uniform ")
    fd.Δagrid[2]-fd.Δagrid[1] == 0 ? (ainfo = "uniform ") : (ainfo = "non-uniform ")
    @printf io "\n"
    @printf io "    Implicit Finite Difference Method\n"
    @printf io "\n"
    @printf io "    Grids  \n"
    @printf io "   ------- \n"
    @printf io "     %12sb: %3d points in [% .0f, %.0f]\n" binfo length(fd.b) fd.b[1] fd.b[end]
    @printf io "     %12sa: %3d poitns in [% .0f, %.0f]\n" ainfo length(fd.a) fd.a[1] fd.a[end]
    # @printf io "     z: %.2f\n" twoap.ξ
end

function FDImp(twoap::TwoAssetsProb;
z::Vector{Float64} = [.8,1.3], λ::Matrix{Float64} = [-1/3 1/3; 1/3 -1/3],
fixedΔb = true, bn::Int =100, bmin::Float64 =-2.0, bmax::Float64 = 40.0, bparam::Vector{Float64} = [0.35, 0.7],
fixedΔa = true, an::Int = 50, amin::Float64 = 0.0, amax::Float64 = 70.0, aparam::Float64 = 0.7,
invΔ::Float64 = .01, invΔᴷ::Float64 = .05)

    #== Parameters ==#
    γ, ρ, ξ, σ, ψ  = _unpackparams(twoap)
    rᴬ, rᴮ, rᴮ⁻, w = _unpackprices(twoap)
    τ, T           = _unpacktax(twoap)


    #== construct the grids ==#
    if fixedΔb
        b   = powerspacegrid(bmin,bmax,bn)
    else
        ngpbPOS = Int(0.8*bn)
        ngpbNEG = bn - ngpbPOS

        bpos = powerspacegrid(1e-8, bmax, ngpbPOS, bparam[2])
        bneg = zeros(ngpbNEG)
        n = Int(ngpbNEG/2+1)
        bneg[1:n] = powerspacegrid(bmin, (bpos[1]+bmin)/2, n , bparam[1])
        b = [bneg;bpos]
        for i in n+1:ngpbNEG
            b[i] = bpos[1] - (b[ngpbNEG+2-i]-b[1])
        end
    end
    fixedΔa ? (a = powerspacegrid(amin,amax,an)) : (a = powerspacegrid(amin,amax,an,aparam[1]))

    zn  = length(z)
    abn = an*bn

    #== Grid counting functions ==#
    global afromab(abi) = div(abi-1,bn)+1
    global bfromab(abi) = rem(abi-1,bn)+1
    global abfromab(ai,bi) = (ai-1)*bn + bi
    global abzfromabz(ai,bi,zi) = (zi-1)*an*bn + ( ai-1 )*bn +bi

    #== Grid info ==#
    Δbgrid = b[2:end]-b[1:end-1]
    Δagrid = a[2:end]-a[1:end-1]

    ΔTbgrid = zeros(bn); ΔTagrid = zeros(an)
    ΔTbgrid[1]       = 0.5*Δbgrid[1]
    ΔTbgrid[2:end-1] = 0.5*( Δbgrid[2:end] + Δbgrid[1:end-1] )
    ΔTbgrid[end]     = 0.5*Δbgrid[end]

    ΔTagrid[1]       = 0.5*Δagrid[1]
    ΔTagrid[2:end-1] = 0.5*( Δagrid[2:end] + Δagrid[1:end-1] )
    ΔTagrid[end]     = 0.5*Δagrid[end]

    #== Labor grid info ==#
    ℓsupply    = zeros(zn)
    ℓsupply[:] = (1/ψ * (1-τ) * w)^σ
    ℓutilgrid  = ψ * (z .* ℓsupply.^( 1 + 1.0/σ ) ) / ( 1+1./σ )

    rbdrift = (b.>= 0.0)*rᴮ + (b.<0.0)*rᴮ⁻
    rbdrift = b .* rbdrift

    #== impose tax on the very top part of grid ==#
    τa = 15
    τc = rᴬ * ( a[end] * 0.999 )^(1 - τa)
    radrift = rᴬ * a - τc * a .^ τa

    #== net_?drifts ==#
    netbinc = zeros(bn,zn)
    netainc = zeros(an,zn)

    for zi in 1:zn
        netℓinc = (1-τ)* w * z[zi] * ℓsupply[zi] + T
        netbinc[:,zi] =  (1-ξ) * netℓinc + rbdrift
        netainc[:,zi] = radrift + ξ * netℓinc
    end

    ## STORAGE

    # NOTE:250 create a sparse matrix of the RIGHT format that later will be filled with zeros
    #       here we have to be careful bc prob has 3 state variables (bi, aj, zk)
    onesλ = Array{Float64,1}[ones( (zn-i)*abn ) for i=1:zn-1]
    colλ  = Int64[(i)*abn for i=1:zn-1]

    A = spdiagm(
        (ones(abn*zn), ones(abn*zn-1), ones(abn*zn-1), ones( abn*zn -bn ) , ones( abn*zn -bn ),
        onesλ..., onesλ...),
        (0, 1, -1, bn, -bn, colλ..., -colλ...)
    )

    fill!(nonzeros(A), zero(Float64))
    Au = deepcopy(A)
    B = deepcopy(A)
    Λ = deepcopy(A)

    # == fill up matrix Λ == #
    ij = zero(Int)
    λt = λ'
    @inbounds for zi in 1:zn, abi in 1:abn
        ij += 1
        Λ[ abi:abn:(zn-1)*abn+abi , ij] = λt[:,zi]
    end

    # == storage arrays == #
    b̃  = Array(Float64, abn*zn)         # RHS of the Bellman equation
    c  = Array(Float64, bn, an, zn)     # consumption policy
    sc = copy(c); d  = copy(c)

    #== Initial value function guess ==#
    V = Array(Float64, abn*zn)
    ij = zero(Int)
    for zi in 1:zn, ai in 1:an, bi in 1:bn
        ij += 1
        V[ij] = ( (1- ξ)*w*z[zi] + rᴬ*a[ai] + rᴮ⁻*b[bi] ).^( 1- γ )/( 1- γ )/ρ
    end

    #== Initial Distribution ==#
    g = zeros(abn,zn)
    bpos_ind = findfirst(b.>0)

    #REVIEW:60 change ydist assumption
    ydist = ones(zn)./zn
    for zi=1:zn
        for abi =1:abn
            ai = afromab(abi)
            bi = bfromab(abi)
            ai==1 && bi==bpos_ind && ( g[abi,zi] = ydist[zi]/(ΔTagrid[ai]*ΔTbgrid[bi]) )
        end
    end
    g = g[:]

    sol = SolutionImp(V, g, c, sc, d, A, Au)

FDImp(
b, Δbgrid, ΔTbgrid, rbdrift, netbinc,
a, Δagrid, ΔTagrid, radrift, netainc,
ℓsupply, ℓutilgrid,
invΔ, invΔᴷ,
z, λ,
sol,
b̃, B, Λ)
end


##################################
##
## TwoAssetsFD
##
##################################
#
# type TwoAssetsFD
#     ## Parameters of the problem
#     twoap::TwoAssetsProb
#
#     ## Discretization
#     b::Vector{Float64}
#     Δb::Float64
#     rbdrift::Vector{Float64}
#     a::Vector{Float64}
#     Δa::Float64
#     radrift::Vector{Float64}
#     invΔ::Float64
#
#     ## Stochastic Information
#     λ::Matrix{Float64}
#     z::Vector{Float64}
#
#     ## Solution
#     V::Vector{Float64}          # Value function
#     gg::Vector{Float64}         # distribution
#
#     ## Storage
#     newV::Vector{Float64}                               # New Value function
#     u::Vector{Float64}                                  # utility term
#     g⁰::Vector{Float64}                                 # for KolmogorovForward
#     A::Base.SparseMatrix.SparseMatrixCSC{Float64, Int}  # rhs matrix on vⁿ⁺¹
#     B::Base.SparseMatrix.SparseMatrixCSC{Float64, Int}  # B = (indΔ + ρ)I - A
#     Λ::Base.SparseMatrix.SparseMatrixCSC{Float64, Int}  # stochastic process terms (useful bc can fill before)
# end
#
#
# function TwoAssetsFD(twoap::TwoAssetsProb,
#                      z::Vector{Float64} = [.8,1.3], λ::Matrix{Float64} = [-1/3 1/3; 1/3 -1/3],
#                      bn::Int =100, bmin::Float64 =-2.0, bmax::Float64 = 40.0,
#                      an::Int = 50, amin::Float64 = 0.0, amax::Float64 = 70.0,
#                      invΔ::Float64 = .01)
#
#     #== Parameters ==#
#     γ = twoap.γ; ρ = twoap.ρ; ξ = twoap.ξ; rᴬ = twoap.rᴬ;
#     rᴮ = twoap.rᴮ; ; w = twoap.w
#     rᴮ⁻ = twoap.rᴮ + twoap.wedge
#     b   = collect(linspace(bmin,bmax,bn))
#     a   = collect(linspace(amin,amax,an))
#     Δb = b[2]-b[1]; Δa = a[2]-a[1];
#
#     rbdrift = ( b.>= 0.0)*rᴮ + (b.<0.0)*rᴮ⁻
#     rbdrift = b .* rbdrift
#     radrift = rᴬ * a
#
#     zn  = length(z)
#     abn = an*bn
#
#     # NOTE:380 create a sparse matrix of the RIGHT format that later will be filled with zeros
#     #       here we have to be careful bc prob has 3 state variables (bi, aj, zk)
#     onesλ = Array{Float64,1}[ones( (zn-i)*abn ) for i=1:zn-1]
#     colλ  = Int64[(i)*abn for i=1:zn-1]
#
#     A = spdiagm(
#         (ones(abn*zn), ones(abn*zn-1), ones(abn*zn-1), ones( abn*zn -bn ) , ones( abn*zn -bn ),
#         onesλ..., onesλ...),
#         (0, 1, -1, bn, -bn, colλ..., -colλ...)
#     )
#
#     fill!(nonzeros(A), zero(Float64))
#     B = deepcopy(A)
#     Λ = deepcopy(A)
#
#     # == fill up matrix Λ == #
#     ij = zero(Int)
#     λt = λ'
#     @inbounds for zi in 1:zn, abi in 1:abn
#         ij += 1
#         Λ[ abi:abn:(zn-1)*abn+abi , ij] = λt[:,zi]
#     end
#
#     V = Array(Float64, abn*zn)
#     ij = zero(Int)
#     for zi in 1:zn, ai in 1:an, bi in 1:bn
#         ij += 1
#         V[ij] = ( (1- ξ)*w*z[zi] + rᴬ*a[ai] + rᴮ⁻*b[bi] ).^( 1- γ )/( 1- γ )/ρ
#     end
#
#     # == b such that Ag = b in plank
#     g⁰ = fill(zero(Float64), abn*zn)
#     i_fix = 1
#     g⁰[i_fix] = .1
#
#     # == storage arrays == #
#     gg = Array(Float64, abn*zn)
#     u  = Array(Float64, abn*zn)
#     newV = deepcopy(V)
#
# TwoAssetsFD(twoap,
#             b, Δb, rbdrift, a, Δa, radrift, invΔ,
#             λ, z,
#             V, gg,
#             newV, u, g⁰, A, B, Λ )
# end
