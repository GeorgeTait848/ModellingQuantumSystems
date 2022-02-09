
include("../underlyingFunctions.jl")

function triangleWave(peak_y, width, x, startOfTriangle)



    if x < startOfTriangle

        y = 0 

        elseif (x-startOfTriangle) < width/2.0

            y = 2*peak_y/width * (x - startOfTriangle)

        

        elseif x-startOfTriangle < width 

            gradient  = -2*peak_y/width

            y_intercept = 2*peak_y

            y = gradient * (x - startOfTriangle) + y_intercept

        

        else 


            y = 0

        
    end

   return y 
    
end





x = -3:0.01:6

y = triangleWave.(2.0, 5.0, x, 0)




trianglePlot = plot(x, y, xlabel = "x (A)", ylabel = "U (eV)")

@time display(trianglePlot)

savefig(trianglePlot, "../../CW accompanying/Images/Tri_Pot.png")

