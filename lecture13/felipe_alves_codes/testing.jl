#== INCLUDE files ==#
include("aggregate.jl")
include("twoassets.jl")

include("solveHJB.jl")
include("solveKFE.jl")

#== Instance of prices ==#
pr = Prices(0.0, 0.04, 0.03, 0.09, 8.0)

#== Instance of TwoAssetProblem ==#
twoap = TwoAssetsProb2(pr; γ = 2.0, ρ = 0.06, ξ = 0.10, χ₀ = 0.08, χ₁= 3.0, τ = 0.0, T = 0.0)

#== Create and FiniteDifference Structure ==#
fde = FDExp(twoap; z = [.8, 1.3], λ = [-1/3 1/3; 1/3 -1/3],
fixedΔa = true, fixedΔb = false, bn = 100, bmin = -2.0, bmax = 40.0, an =70, amax = 70.0, invΔᴷ = 0.05)

#== Solve the HJB equation ==#
solve_hjb!(twoap, fde, fde.sol)

#== Solve the KF equation ==#
solve_fp!(fde, fde.sol, maxit = 4000, tolFK = 1e-9)

include("main_fig.jl");

#== Policies ==#
# fig_pol(fde)
#== Distributions ==#
fig_dens(fde)
