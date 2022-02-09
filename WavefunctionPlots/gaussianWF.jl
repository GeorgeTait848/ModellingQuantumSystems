include("../PotentialEnergies/GaussianFunction.jl")
include("../parameters.jl")
include("../underlyingFunctions.jl")

x = -10*A:0.01*A:10*A
electronEnergy = 1.75*eV
potentialEnergies = gaussianFunction.(3*eV, 0, 5*A, x)


gaussianWF = waveFunction(electronEnergy)


waveFuncPlot = getWaveFunctionPlot(electronEnergy)

@time display(waveFuncPlot)

savefig(waveFuncPlot, "../../CW accompanying/Images/gauss_wf.png")

