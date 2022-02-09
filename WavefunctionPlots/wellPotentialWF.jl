include("../parameters.jl")

include("../underlyingFunctions.jl")

include("../PotentialEnergies/WellPotential.jl")







x = -6*A:0.01*A:15*A

potentialEnergies = finiteWellPotential.(0, 8*A, 0, 3*eV, x)

# t_11s = get_t_11.(electronEnergies)

# plot(electronEnergies/eV, t_11s)



BSWFPlot = getBoundStateWaveFunctionPlot(1e-12)

@time display(BSWFPlot)


