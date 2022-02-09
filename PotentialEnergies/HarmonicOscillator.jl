include("../parameters.jl")

include("../underlyingFunctions.jl")


function getHarmonicOscillatorPotential(offset, ang_freq, x)


    y = (1/2) * electronMass * ang_freq^2 * (x-offset)^2

    return y
    
end


x = -4*A :0.01*A : 4*A

y = getHarmonicOscillatorPotential.(0, 1.5e14, x)




plot(x/A, y/eV)

