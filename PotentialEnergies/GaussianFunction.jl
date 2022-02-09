include("../underlyingFunctions.jl")


function gaussianFunction(height, mean, fullWidthHalfMax, x)

    standardDeviation  = fullWidthHalfMax / (2*sqrt(2*log(2)))


    y = height*exp(-1*(x - mean)^2/(2(standardDeviation)^2))


    return y     

end




# x = -8*A:0.1*A:8*A


# y = gaussianFunction.(3*eV, 0, 5*A, x)


# gaussianPulsePlot = plot(x/A, y/eV, xlabel = "x (A)", ylabel = "U(eV)");



# display(gaussianPulsePlot)




