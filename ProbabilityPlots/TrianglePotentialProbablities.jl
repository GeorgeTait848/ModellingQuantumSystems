include("../PotentialEnergies/TriangleWave.jl")
include( "../underlyingFunctions.jl")
include("../parameters.jl")


x = -3*A: 0.1*A:6*A


potentialEnergies = triangleWave.(2*eV, 5*A, x, 0)

electronEnergies = 0.01*eV: 0.01*eV: 4.8*eV



@time triangleProbPlot = getProbabilitiesPlot(electronEnergies)
display(triangleProbPlot)
savefig(triangleProbPlot, "../../CW accompanying/Images/Tri_Prob_120.png")
