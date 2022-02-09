include("parameters.jl")



using Plots


struct NoSolutionsError <: Exception 



end

struct EnergiesOutOfRangeError <: Exception

end

struct InvalidPotentialsError <: Exception

end


function assertValidPotentials()

    if length(x) != length(potentialEnergies)

    
        throw(InvalidPotentialsError)

    end




    
end





function getWaveNumber(E,U)

        k = sqrt(complex((E-U) * 2 * electronMass)) / hbar



    

    return k




end



function getWavenumbersArray(E, U)
    
    wavenumbersArray = getWaveNumber.(E, U)


    return wavenumbersArray
    
end






function getTransferMatrixAtBoundary(k_i, k_j, boundaryPosition)

    if k_i == 0 

        t11 = (1 - 1im*k_j*boundaryPosition)*exp(1im*k_j*boundaryPosition)

        t12 = (1 + 1im*k_j*boundaryPosition)*exp(-1im*k_j*boundaryPosition)

        t21 = 1im * k_j * exp(1im*k_j*boundaryPosition)

        t22 = -1im * k_j * exp(-1im*k_j*boundaryPosition)
        

        T = [t11 t12 ; t21 t22]


    elseif k_j == 0 

        t11 = (1/2)*exp(-1im * k_i * boundaryPosition)

        t12 = (1/(2*k_i)) * (k_i * boundaryPosition - 1im) * exp(-1im * k_i * boundaryPosition)


        t21 = (1/2)*exp(1im * k_i * boundaryPosition)

        t22 = (1/(2*k_i)) * (k_i * boundaryPosition + 1im) * exp(1im * k_i * boundaryPosition)

        

        T = [t11 t12 ; t21 t22]


    else

        t11 = (k_i + k_j) * exp(-1im * boundaryPosition * (k_i - k_j))


        t12 = (k_i - k_j) * exp(-1im * boundaryPosition * (k_i + k_j))
    
        t21 = (k_i - k_j) * exp(1im * boundaryPosition * (k_i + k_j))
    
        t22 = (k_i + k_j) * exp(1im * boundaryPosition * (k_i - k_j))

        T = [t11 t12 ; t21 t22]

        T *= (1/(2*k_i))


    end
  


    return T
    
end

function getFirstBoundaryPosition()

    return x[2]


    
end



function getSystemTransferMatrix(electronEnergy)

    dx = x[2] - x[1]

    

    wavenumbers = getWavenumbersArray(electronEnergy, potentialEnergies)

    boundaryPosition = getFirstBoundaryPosition()

    currentTransferMatrix = getTransferMatrixAtBoundary(wavenumbers[1], wavenumbers[2], boundaryPosition)
    


        for i in 2:length(wavenumbers)-1

            boundaryPosition+=dx

            boundaryMatrix = getTransferMatrixAtBoundary(wavenumbers[i], wavenumbers[i + 1], boundaryPosition)
            

            currentTransferMatrix = currentTransferMatrix * boundaryMatrix
            

        end

        
        
    return currentTransferMatrix
    
end


function getRegionCoefficients(transferMatrix, previousCoefficients)


    newCoefficients = transferMatrix * previousCoefficients


    return newCoefficients
    
end


function calculateTransmissionProbability(electronEnergy)

    assertValidPotentials()
    
    wavenumbersArray = getWavenumbersArray(electronEnergy, potentialEnergies)
    
    systemTransferMatrix = getSystemTransferMatrix(electronEnergy)
   
    t11 = systemTransferMatrix[1]

    
    finalWaveNumber = wavenumbersArray[end]

    firstWaveNumber = wavenumbersArray[1]

    

    probability = abs(1/t11)^2 * (finalWaveNumber/firstWaveNumber) 

    probability = real(probability)
    

    return probability

    
end




function calculateReflectionProbability(electronEnergy)
    

    assertValidPotentials()
    
    systemTransferMatrix = getSystemTransferMatrix(electronEnergy)

    t11 = systemTransferMatrix[1]
    t21 = systemTransferMatrix[2]
    
    probability = (abs(t21/t11))^2

    probability = real(probability)

    return probability


    
end

function getProbabilitiesPlot(electronEnergies)

    transmissionProbabilities = calculateTransmissionProbability.(electronEnergies)


    reflectionProbabilties = calculateReflectionProbability.(electronEnergies)

    plottableEnergies = electronEnergies./eV

    probabilityPlot = plot(plottableEnergies, reflectionProbabilties, label = "R(E)")


    probabilityPlot = plot!(plottableEnergies, transmissionProbabilities, label = "Tr(E)")

    probabilityPlot = plot!(xlabel = "E (eV)", ylabel = "Probability")


    return probabilityPlot




    
end




function waveFunction(electronEnergy)

    assertValidPotentials()

    boundarySeparation = x[2] - x[1]
    
    y = complex(zeros(length(x)))  
 ## Working from final region to first region (i.e right of grid to left)
    regionCounter = 0

    k_j = getWaveNumber(electronEnergy, potentialEnergies[end])

    k_i = getWaveNumber(electronEnergy, potentialEnergies[end - 1])


    firstBoundaryPosition  = getFirstBoundaryPosition()

    boundaryPosition = firstBoundaryPosition + (length(potentialEnergies) - regionCounter - 2) * boundarySeparation 

    #Since first boundary between regions 1,2 is at x[begin] + boundarySeparation: the boundary position  
    # between the nth and (n+1)th regions is firstBoundary + (n-2)boundarySeparation. where n = N - regionCounter.

    coefficients_ij = finalRegionCoefficients



    for i in 0:(length(x)-2)


        y[end-i] = coefficients_ij[1] * exp(1im*k_j*x[end - i]) + coefficients_ij[2] * exp(-1im * k_j * x[end - i])

        #for the next x, we will have changed regions, thus we must find the coefficients of the next region to the left. 


        transferMAtrix_ij = getTransferMatrixAtBoundary(k_i, k_j, boundaryPosition)

        coefficients_ij = getRegionCoefficients(transferMAtrix_ij, coefficients_ij)

        #Now reset the boundary position to be the start of this region


        boundaryPosition -= boundarySeparation

        k_j = k_i
       

        k_i = getWaveNumber(electronEnergy, potentialEnergies[end - regionCounter - 1])

        regionCounter +=1 
        

    end

    #When the for loop is termninated, it will have assigned a value to each element of y apart from y[1], but the coefficients will be 
    #correct for this region

   y[1] = coefficients_ij[1] * exp(1im*k_j*x[1]) + coefficients_ij[2] * exp(-1im * k_j * x[1])

    return y 

    
end



function getPotentialEnergyPlot()

    plottableX = (1/A)*x

    plottableU = (1/eV)*potentialEnergies

    potentialEnergyPlot = plot(plottableX, plottableU, linecolor = "black", label = "U")

    potentialEnergyPlot = plot!(xlabel = "x (Ang)", ylabel = "U (eV)")

    return potentialEnergyPlot



    
end


function getWaveFunctionPlot(electronEnergy)

    plottableU = potentialEnergies / eV

    plottableX = x / A

    y = waveFunction(electronEnergy)

    wavefuncPlot = getPotentialEnergyPlot()

    max_y = if maximum(real(y)) >= maximum(imag(y))

        maximum(real(y))

    else 

        maximum(imag(y))

    end

    max_U = maximum(plottableU)


    ratio = max_y/max_U

    y *= (1/ratio)
    

    wavefuncPlot = plot!(plottableX, real(y), label = "Re(psi)")

    wavefuncPlot = plot!(plottableX, imag(y), label = "Im(psi")


    return wavefuncPlot

    
end 




function checkIfPositive(number)

    return number >= 0
   
end


function checkForSolution(number1, number2)

    sign1 = checkIfPositive(number1)

    sign2 = checkIfPositive(number2)

    if number1 == 0 || number2 == 0

        return true

    else

        return sign1 != sign2

    end

    
end

function findFirstRootIndex(numbers)

    currentIndex = 1

    changeOfSign = false

    #check Solutions remain

    while changeOfSign == false && currentIndex < length(numbers)

        currentIndex += 1

        changeOfSign = checkForSolution(numbers[1], numbers[currentIndex])

        
        #This loop will be terminated if there is a change of sign or if 
        #the end of the array is reached 

    end

    if changeOfSign == false #i.e if the loop is terminated as reached end of array

        throw(NoSolutionsError)
    else

        return currentIndex

    end
    
end


function findAllRootIndices(numbers)

    currentNumbers = copy(numbers)

    indexInArray = 1

    indicesofSolutions= []

    currentSolutionIndex =  try
        
        findFirstRootIndex(numbers)
        
    catch

        nothing
    end

    while currentSolutionIndex !== nothing

        append!(indicesofSolutions, currentSolutionIndex)

        toBeDeletedIndices = 1:(currentSolutionIndex - indexInArray + 1)

        #To ensure when looking for the third, fourth etc solution you are removing the number 
        #of indices based on newnumbers and not numbers. 
        #i.e if 1st solution is nth element of numbers, it is the 1st element of newnumbers
        #thus the second element is the mth element of numbers and the (m + n -1)th element of new numbers.

        currentNumbers = deleteat!(currentNumbers, toBeDeletedIndices)

        indexInArray = currentSolutionIndex + 1


        currentSolutionIndex = try

            findFirstRootIndex(currentNumbers) + indexInArray - 1

            #this is to get the index relative to numbers not newnumbers, which is a shortened version. 
            
        catch

            nothing
        end 
    end

    if indicesofSolutions == []


        throw(NoSolutionsError)


    else 

        return indicesofSolutions

    end
    
end






function checkEquivalenceUnderTolerance(number1, number2, tolerance)

    
    magDiff = abs(number1 - number2)


    return magDiff < tolerance
    
end





function get_t_11(electronEnergy)


    systemTransferMatrix = getSystemTransferMatrix(electronEnergy)


    t_11 = real(systemTransferMatrix[1])

    return t_11
    
end


function refineBoundState(unrefinedSolution, de, tolerance)
    tenth_de = de/10

    solutionRange = (unrefinedSolution-de) : tenth_de : (unrefinedSolution)

    t_11_range = get_t_11.(solutionRange)

    refinedSolutionIndex = findFirstRootIndex(t_11_range)


    refined_t_11_root = t_11_range[refinedSolutionIndex]


    refinedSolution = solutionRange[refinedSolutionIndex]

    belowTolerance = checkEquivalenceUnderTolerance(refined_t_11_root, 0, tolerance)


    if belowTolerance == true 

        return refinedSolution

    else 
        return refineBoundState(refinedSolution, tenth_de, tolerance)

    end

    



    
end


function assertEnergiesInRange(electronEnergies, potentialEnergies)


    if minimum(electronEnergies) < minimum(potentialEnergies)

        minPE = minimum(potentialEnergies) / eV

        println("minimum electron energy too low. Must be greater than $minPE eV")
        throw(EnergiesOutOfRangeError)

    elseif maximum(electronEnergies) > maximum(potentialEnergies)

        maxPE = maximum(potentialEnergies) / eV

        println("maximum electron energy too high. Must be less than $maxPE eV")

        throw(EnergiesOutOfRangeError)

    end

    
end


function findBoundStateEnergies(electronEnergies, tolerance)


    assertEnergiesInRange(electronEnergies, potentialEnergies)

    t_11s = get_t_11.(electronEnergies)


    unrefinedSolutionsIndices = findAllRootIndices(t_11s)

    unrefinedSolutions = electronEnergies[unrefinedSolutionsIndices]

    de = electronEnergies[2]-electronEnergies[1]

    refiendSolutions = refineBoundState.(unrefinedSolutions, de, tolerance)


    return refiendSolutions
    
end

function getElectronEnergies(potentials)


    maxPotential = maximum(potentials)

    minPotential = minimum(potentials)

    electron_energies = (minPotential + 0.01*eV):0.01*eV:(maxPotential-0.0001*eV)
    
    return electron_energies

end



function getBoundStateWaveFunctionPlot(tolerance)
    
    electronEnergies = getElectronEnergies(potentialEnergies)

    
    
    max_U = maximum(potentialEnergies)/eV
    

    boundStateEnergies = findBoundStateEnergies(electronEnergies, tolerance)

    

    boundStateCounter = 1


    plottableX = (1/A)*x

    wavefuncPlot = getPotentialEnergyPlot()




    while boundStateCounter <= length(boundStateEnergies)

        

        currentBoundStateEnergy = boundStateEnergies[boundStateCounter]

        BSEnergy_plottable = fill(currentBoundStateEnergy /eV, length(x))

        currBoundStateEnergy_2dp = round(currentBoundStateEnergy/eV, digits = 4)

        wavefuncPlot = plot!(plottableX, BSEnergy_plottable, linecolor = "darkgray", label = "E_$boundStateCounter = $(currBoundStateEnergy_2dp) eV")

        y = waveFunction(currentBoundStateEnergy)
    
    
        max_y = maximum(real(y))
    
        ratio = max_y/max_U
    
        y *= 1/ratio

        y .+= currentBoundStateEnergy/eV
       
    
        wavefuncPlot = plot!(plottableX, real(y), label = "\U003C8_$boundStateCounter")

      
    
        
        boundStateCounter += 1 



    end


    wavefuncPlot = plot!(legend =:bottomright)
    

    return wavefuncPlot

end

function getNthBoundStateWaveFunctionPlot(tolerance, n)

    electronEnergies = getElectronEnergies(potentialEnergies)

    
    
    max_U = maximum(potentialEnergies)/eV
    

    boundStateEnergies = findBoundStateEnergies(electronEnergies, tolerance)

    

    nthBoundStateEnergy = boundStateEnergies[n]


    plottableX = (1/A)*x

    wavefuncPlot = getPotentialEnergyPlot()

    y =  waveFunction(nthBoundStateEnergy)


    max_y = maximum(real(y))
    
    ratio = max_y/max_U

    y *= 1/ratio

    y .+= nthBoundStateEnergy/eV
   
    nthBoundStateEnergy_eV = round(nthBoundStateEnergy / eV, digits=3)

    nthBoundStateEnergyPlottable = fill(nthBoundStateEnergy/eV, length(x))

    wavefuncPlot = plot!(plottableX, nthBoundStateEnergyPlottable, label = "E_$n = $nthBoundStateEnergy_eV eV ", linecolor = "darkgray")
    
    
    wavefuncPlot = plot!(plottableX, real(y), label = "\U003C8_$n")


    return wavefuncPlot

end


function plotBoundStateEnergies(tolerance)

    electronEnergies = getElectronEnergies(potentialEnergies)

    plottableX = x / A
    energiesPlot = getPotentialEnergyPlot()

    boundStateEnergies = findBoundStateEnergies(electronEnergies, tolerance)


    boundStateCounter = 1


   


    while boundStateCounter <= length(boundStateEnergies)

        currentBoundStateEnergy = round(boundStateEnergies[boundStateCounter] / eV, digits = 3)

        plottableBSEnergy = fill(currentBoundStateEnergy, length(x))

        energiesPlot = plot!(plottableX, plottableBSEnergy, label = "E_$boundStateCounter = $currentBoundStateEnergy")

        boundStateCounter += 1

        
    end


    return energiesPlot

end
