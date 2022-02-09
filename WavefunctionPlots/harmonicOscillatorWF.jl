include("../PotentialEnergies/HarmonicOscillator.jl")

include("../parameters.jl")

include("../underlyingFunctions.jl")



x = -20*A : 0.1*A : 20*A

potentialEnergies = getHarmonicOscillatorPotential.(0, 1.5e14, x)




HOPlot = getBoundStateWaveFunctionPlot(1e-10)

@time display(HOPlot)
savefig(HOPlot, "../../CW accompanying/Images/HOWF.png")