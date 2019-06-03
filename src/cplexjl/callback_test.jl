using JuMP
using Gurobi

# We will use Gurobi, which requires that we manually set the attribute
# PreCrush to 1 if we have user cuts. We will also disable PreSolve, Cuts,
# and Heuristics so only our cut will be used
m = Model(solver=GurobiSolver(PreCrush=1, Cuts=0, Presolve=0, Heuristics=0.0))

# Define our variables to be inside a box, and integer
@variable(m, 0 <= x <= 2, Int)
@variable(m, 0 <= y <= 2, Int)

# Optimal solution is trying to go towards top-right corner (2.0, 2.0)
@objective(m, Max, x + 2y)

# We have one constraint that cuts off the top right corner
@constraint(m, y + x <= 3.5)

# Optimal solution of relaxed problem will be (1.5, 2.0)
# We can add a user cut that will cut of this fractional solution.

# We now define our callback function that takes one argument,
# the callback handle. Note that we can access m, x, and y because
# this function is defined inside the same scope
function mycutgenerator(cb)
    x_val = getvalue(x)
    y_val = getvalue(y)
    println("In callback function, x=$x_val, y=$y_val")

    # Allow for some impreciseness in the solution
    TOL = 1e-6

    # Check top right
    if y_val + x_val > 3 + TOL
        # Cut off this solution
        println("Fractional solution was in top right, cut it off")
        # Use the original variables
        return JuMP.StopTheSolver
    end
end  # End of callback function

# Tell JuMP/Gurobi to use our callback function
addheucallback(m, mycutgenerator)

# Solve the problem
solve(m)

# Print our final solution
println("Final solution: [ $(getvalue(x)), $(getvalue(y)) ]")

#################################################################

using JuMP
using Gurobi

# We will use Gurobi and disable PreSolve, Cuts, and (in-built) Heuristics so
# only our heuristic will be used
m = Model(solver=GurobiSolver(Cuts=0, Presolve=0, Heuristics=0.0))

# Define our variables to be inside a box, and integer
@variable(m, 0 <= x <= 2, Int)
@variable(m, 0 <= y <= 2, Int)

# Optimal solution is trying to go towards top-right corner (2.0, 2.0)
@objective(m, Max, x + 2y)

# We have one constraint that cuts off the top right corner
@constraint(m, y + x <= 3.5)

# Optimal solution of relaxed problem will be (1.5, 2.0)

# We now define our callback function that takes one argument,
# the callback handle. Note that we can access m, x, and y because
# this function is defined inside the same scope
function myheuristic(cb)
    x_val = getvalue(x)
    y_val = getvalue(y)
    println("In callback function, x=$x_val, y=$y_val")
    return JuMP.StopTheSolver

    # setsolutionvalue(cb, x, floor(x_val))
    # # Leave y undefined - solver should handle as it sees fit. In the case
    # # of Gurobi it will try to figure out what it should be.
    # addsolution(cb)
    #
    # # Submit a second solution
    # setsolutionvalue(cb, x, ceil(x_val))
    # addsolution(cb)
end  # End of callback function

# Tell JuMP/Gurobi to use our callback function
addheuristiccallback(m, myheuristic)

# Solve the problem
solve(m)

# Print our final solution
println("Final solution: [ $(getvalue(x)), $(getvalue(y)) ]")

#################################################################

type NodeData
    time::Float64  # in seconds since the epoch
    node::Int
    obj::Float64
    bestbound::Float64
end

# build model ``m`` up here

bbdata = NodeData[]

function infocallback(cb)
    node      = MathProgBase.cbgetexplorednodes(cb)
    obj       = MathProgBase.cbgetobj(cb)
    bestbound = MathProgBase.cbgetbestbound(cb)
    push!(bbdata, NodeData(time(),node,obj,bestbound))
    if node >= 2
        return JuMP.StopTheSolver
    end
end
addinfocallback(m, infocallback, when = :Intermediate)

solve(m)

# Save results to file for analysis later
open("bbtrack.csv","w") do fp
    println(fp, "time,node,obj,bestbound")
    for bb in bbdata
        println(fp, bb.time, ",", bb.node, ",",
                    bb.obj, ",", bb.bestbound)
    end
end
