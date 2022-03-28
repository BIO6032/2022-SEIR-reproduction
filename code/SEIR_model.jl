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
    D(S) ~ -β*S*A - β*S*Ifq*I - e*S,
    D(A) ~ β*S*A + β*S*I - fs*A - fr*A - s*ft*A + e*S,
    D(I) ~ fs*A - fr*I - fq*I,
    D(Q) ~ fq*I + s*ft*a - fr*Q,
    D(R) ~ fr*Q + fr*I + fr*A
]

# Problem 
@named seir = ODESystem(equations)



## Models ##



# R_0 = 8.0




# R_0 = 6.5 




# R_0 = 6.0 




# R_0 = 5.0




# R_0 = 4.0 




# R_0 = 3.8




# R_0 = 2.5 



# Figure 2A

u0 = [S => , A => , I => , R => , Q => ]
p = [ β => , e => , fs => , fr => , s => 0.9, ft => , fq => ]
tspan = (0., 100.)






# Figure 2B








# Figure 3A








# Figure 3B