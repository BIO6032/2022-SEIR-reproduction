# Simulations
using DifferentialEquations
using OrdinaryDiffEq
using ModelingToolkit

# Figures
using Plots

# Data
using DataFrames
#using CSV:CSV

# Parameters 
@parameters t β e fs fr s ft fq

# Variables
@variables S(t) A(t) I(t) R(t) Q(t)

# Differential
D = Differential(t)

# Equations 
equations = [
    D(S) ~ -β * S * A - β * S * I - e * S,
    D(A) ~ β * S * A + β * S * I - fs * A - fr * A - s * ft * A + e * S,
    D(I) ~ fs * A - fr * I - fq * I,
    D(Q) ~ fq * I + s * ft * A - fr * Q,
    D(R) ~ fr * Q + fr * I + fr * A,
]

# Problem 
@named seir = ODESystem(equations)
u0 = [S => 9999, A => 1, I => 0, R => 0, Q => 0]
tspan = (0.0, 100.0)

# Simulations
for r0 in [2.5, 3.8, 4.0, 5.0, 6.0, 6.5, 8.0]
    # TODO formula for parameters from R0
    p = [
        fr => 1 / 14,
        β => r0 * fr / (S + A + I + R + Q),
        e => 4 / 10000,
        fs => 0.6,
        s => 0.9,
        ft => 2 / 7,
        fq => 0.9
    ]
    problem = ODEProblem(seir, u0, tspan, p)
    solution = solve(problem)
    plot(solution; dpi = 600, frame = :box)
    savefig(joinpath("figures", "simulation_$(r0).png"))
    # TODO save data to a CSV file
end


