# Simulations
using DifferentialEquations
using OrdinaryDiffEq
using ModelingToolkit 

# Figures
using Plots

# Data
using DataFrames
using CSV: CSV

# Parameters 
@parameters t β e fs fr s ft fq

# Variables
@variables S(t) A(t) I(t) R(t) Q(t)

# Differential
D = Differential(t)

# Equations 
equations =  [
    D(S) ~ -β*S*A - β*S*I - e*S,
    D(A) ~ β*S*A + β*S*I - fs*A - fr*A - s*ft*A + e*S,
    D(I) ~ fs*A - fr*I - fq*I,
    D(Q) ~ fq*I + s*ft*a - fr*Q,
    D(R) ~ fr*Q + fr*I + fr*A
]

# Problem 
@named seir = ODESystem(equations)
u0 = [S => , A => , I => , R => , Q => ]
tspan = (0., 100.)

# Simulations
for r0 in [2.5, 3.8, 4.0, 5.0, 6.0, 6.5, 8.0]
    # TODO formula for parameters from R0
    p = [ β => , e => , fs => , fr => , s => 0.9, ft => , fq => ]
    problem = ODEProblem(seir, u0, tspan, p)
    solution = solve(problem)
    plot(solution; dpi=600, frame=:box)
    savefig(joinpath("figures", "simulation_$(r0).png"))
    # TODO save data to a CSV file
end
