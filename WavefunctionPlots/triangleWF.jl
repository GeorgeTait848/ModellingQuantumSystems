include("../PotentialEnergies/TriangleWave.jl")
include("../parameters.jl")
include("../underlyingFunctions.jl")


x = -5*A:0.01*A:9*A

potentialEnergies = triangleWave.(2*eV, 5*A, x, 0)



waveFuncPlot = getWaveFunctionPlot(0.75*eV)


waveFuncPlot = plot!(title = "Wave function of triangle wave potential energy")



display(waveFuncPlot)

# savefig(waveFuncPlot, "../../CW accompanying/Images/Tri_WF.png")
