##################################
##
## Prices
##
##################################
"""
Type that contains information of prices

## Prices

- `rᴷ`      : rental rate of capital
- `rᴬ`      : return on the illiquid asset
- `rᴮ`      : return on the liquid asset
- `wedge`   :
- `w`       : wage rate

"""
type Prices
    rᴷ ::Float64
    rᴬ ::Float64
    rᴮ ::Float64
    wedge::Float64
    w  ::Float64
end

function update_prices!(pc::Prices, P::Vector{Float64})
    pc.rᴷ, pc.rᴬ, pc.rᴮ, pc.w = P
end

price_partial() = Prices(0.0, 0.04, 0.03, 0.09,4.0)

function Base.show(io::IO, pr::Prices)
    @printf io "\n"
    @printf io "    Prices \n"
    @printf io "      rᴷ    :    %.3f \n" pr.rᴷ
    @printf io "      rᴬ    :    %.3f \n" pr.rᴬ
    @printf io "      rᴮ    :    %.3f \n" pr.rᴮ
    @printf io "      wedge :    %.3f \n" pr.wedge
    @printf io "      wage  :    %.3f \n" pr.w
end

##################################
##
## SteadyState and Transition
##
##################################
"""
Description of the SteadyState
# Endogenous

### Prices
- `rcapital`
- `ra`
- `rb`
- `rborr`
- `wage `

### Interest rate
- `rnom`

### Aggregates to clear
- `Eb`      : aggregate liq saving from households
- `bond`    : total amount of bonds in the market
- `Ea`      : aggregate illiq saving from households
- `capital` : capital demanded by firms
- `labor`   :

### Firm Stats
- `KYratio`
- `KNratio`

- `π`
- `mc`
- `priceadjust`

- `output`
- `profit`
- `dividend`

### Fund
- `divrate`
- `investment`
- `deprec`
- `caputil`

### Government
- `G`
- `lumptransfer`
- `τ`
- `taxrev`
- `govbond`

"""

function initializeSS(KYratio = 8., rb = 0.015)

    ssvar = [:rcapital, :ra, :rb, :wedge, :wage,            # prices
             :rnom,                                         #
             :Eb, :bond, :Ea, :capital, :labor,             # aggregates
             :KYratio, :KNratio, :π, :mc, :priceadjust,     #
             :output, :profit, :dividend,                   #
             :divrate, :investment, :deprec, :caputil,      #
             :G, :lumptransfer, :τ, :taxrev, :govbond,      #
             :tfp]                                          #

    equmSS = (Symbol => Float64)[var => 0.0 for var in ssvar]

    equmSS[:tfp]      = 1.
    equmSS[:KYratio]  = KYratio
    equmSS[:τ]        = 0.25
    equmSS[:rb]       = rb    #having problem with 0.02/4

    equmSS[:wedge]    = 0.10 #0.0161134


    equmSS[:KNratio] = ( equmSS[:tfp] * equmSS[:KYratio] )^(1./(1.-α))
    equmSS[:mc] = (ϵ-1) / ϵ

    #== Prices ==#
    equmSS[:rcapital] = equmSS[:mc] * α / equmSS[:KYratio]
    equmSS[:wage]     = equmSS[:mc] * equmSS[:tfp] * (1. - α) * equmSS[:KNratio]^α

    #== Firm ==#
    equmSS[:priceadjust] = 0
    equmSS[:labor]       =  ( (1 - equmSS[:τ])*equmSS[:wage] /ψ )^σ

    equmSS[:capital] = equmSS[:KNratio] * equmSS[:labor]
    equmSS[:investment] = δbar * equmSS[:capital]

    equmSS[:profit]  = (1. - equmSS[:mc])* equmSS[:capital] / equmSS[:KYratio] - equmSS[:priceadjust]
    equmSS[:divrate] = equmSS[:profit] / equmSS[:capital]

    #== Illiquid return ==#
    equmSS[:ra] = equmSS[:rcapital] - δbar + equmSS[:divrate]

    #== Government ==#
    equmSS[:lumptransfer] = 0.20 * equmSS[:wage] * equmSS[:labor]

    ### Solution ###

    return equmSS

end

_unpackprices(equmSS::Dict{Symbol, Float64}) =
    Prices(equmSS[:rcapital], equmSS[:ra], equmSS[:rb], equmSS[:wedge] , equmSS[:wage])

##################################
## OLD AGGREGATE vars
##################################

"""
Holds the aggregate variables of the econmy

### Prices

- `prices`

### Government
- `Rev::Float64`
- `τ::Float64`
- `T::Float64`

### Firm
- `tfp::Float64`
- `KY::Float64`
- `KN::Float64`

- `output::Float64`
- `profit::Float64`
- `capital::Float64`
- `investment::Float64`

"""
type AggVar

    prices::Prices

    ##Gov
    Rev::Float64
    τ  ::Float64
    T  ::Float64
    G  ::Float64

    ## Firm
    tfp::Float64
    KY ::Float64
    KN ::Float64

    output::Float64
    profit::Float64
    capital::Float64
    investment::Float64

end

function AggVar(;tfp::Float64 = 1., KY::Float64 = 5.,
                τ::Float64 = 0.25, rᴮ::Float64 = 0.03, wedge::Float64 = 0.09)

    KN = (tfp * KY)^(1./(1.-α))

    #== Prices ==#
    rᴷ = mc * α / KY
    w  = mc * tfp * (1. - α) * KN^α
    priceadjust = 0
    ℓ = 1.0
    capital = KN * ℓ
    investment = δbar * capital
    profit  = (1. - mc)* capital / KY - priceadjust
    divr = profit / capital

    rᴬ = rᴷ - δbar + divr

    # Define this so that a fixed fraction receives more transfer than pays
    T = 0.8 * w

    output = 0.
    Rev = 0.0
    G = 0.0

    AggVar(Prices(rᴷ,rᴬ,rᴮ,wedge,w), Rev, τ, T, G, tfp, KY, KN,
    output, profit, capital, investment)
end
