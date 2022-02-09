


include("../underlyingFunctions.jl")


using Plots
function finiteWellPotential(startOfWell, wellWidth, wellBottom, wellTop, x)


    if x < startOfWell

        y = wellTop


    elseif x < startOfWell + wellWidth


        y = wellBottom

    else

        y = wellTop

    end

    return y





    
end


function coupledWellPotentials(startOfFirstWell, wellWidths, wellBottoms, wellTops, wellSeparations, numberOfWells, x)

    endOfWell = startOfFirstWell + wellWidths + wellSeparations
    currentStart = startOfFirstWell
    

    y = nothing

    


    for _ in 1:numberOfWells

        if x < endOfWell


            y = finiteWellPotential(currentStart, wellWidths, wellBottoms, wellTops, x)

            return y

            


        else 

            currentStart = endOfWell

            endOfWell += wellWidths + wellSeparations

        

        end

    end



    if y === nothing


        y = wellTops

    end


    return y

    
end




#Plot for a single Well


# x = -4:0.01:11


# myWell = finiteWellPotential.(0, 8, 0, 5, x)

# wellPotPlot = plot(x, myWell, xlabel = "x (A)", ylabel = "U (ev)");




#plot for coupled wells 


x = -1:0.01:20


coupledWells = coupledWellPotentials.(0.0, 3.0, 0.0, 3.0, 1.0, 5.0, x)

plot(x, coupledWells)





