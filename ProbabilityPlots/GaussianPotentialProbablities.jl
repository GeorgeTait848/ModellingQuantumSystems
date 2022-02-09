include("../PotentialEnergies/GaussianFunction.jl")
include( "../underlyingFunctions.jl")
include("../parameters.jl")


electronEnergies = 1.0eV:0.01eV:5eV 

x = -16*A : 1*A : 16*A

potentialEnergies = gaussianFunction.(3.0*eV, 0, 8.0*A, x)


println(getSystemTransferMatrix(electronEnergies[1]))

@time gaussProb = getProbabilitiesPlot(electronEnergies)

display(gaussProb)

