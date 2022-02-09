include("../PotentialEnergies/WellPotential.jl")

function getCoupledWellX(startOfFirstWell, wellWidths, wellSeparations, numberOfWells, distanceFromEnd, dx)

    #distanceFromEnd was hard to name, it refers to how far from the edges of the well you want to plot. 


    #this function is to ensure youre not plotting coupled wells which are incomplete, i.e your range of x is not to small to fit all the values. 

    startOfX = startOfFirstWell - distanceFromEnd

    
    
    
    endOfX = startOfFirstWell + numberOfWells*(wellWidths + wellSeparations) + distanceFromEnd
    


    x = startOfX:dx:endOfX

    return x

    
end

x = getCoupledWellX(0.0, 5*A, 1*A, 5, 5*A, 0.01*A)



potentialEnergies = coupledWellPotentials.(0.0, 5.0*A, 0.0, 2*eV, 1*A, 5, x)




coupledWellsBSWFPlots = getBoundStateWaveFunctionPlot(1e-10)



@time display(coupledWellsBSWFPlots)
# savefig(coupledWellsBSWFPlots, "../../CW accompanying/Images/5CoupledWellsWF.png")

