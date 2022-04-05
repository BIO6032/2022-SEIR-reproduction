# Simulations
using DifferentialEquations
using OrdinaryDiffEq
using ModelingToolkit

# Figures
using Plots

# Data
using DataFrames
using CSV:CSV

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
#for r0 in [2.5, 3.8, 4.0, 5.0, 6.0, 6.5, 8.0]
    # TODO formula for parameters from R0
    r0 = 8.

    # Get Parameters values
    _fr = 1/14
    _β = r0 * _fr / (S + A + I + R + Q)
    _e = 4 / 10000
    _fs = 0.6 * _fr
    _s = 0.9
    _ft = 2/7
    _fq = 0.1
    _Ss = 0.6 #constant given in the article (#To verify)
    _Tr = 14. #given in the article (number of days of infection)


    p = [
        fr => _fr,
        β => _β,
        e => _e,
        fs => _fs,
        s => _s,
        ft => _ft,
        fq => _fq
    ]
    problem = ODEProblem(seir, u0, tspan, p)
    solution = solve(problem)

   # Visualize the solution
    plot(solution; dpi=600, frame=:box)
    
    M = solution[1:end,:] # Matrix of results 5 compartments x 31 time steps
    
    # Assign values to each compartment over time
    S = solution[1,:]
    A = solution[2,:]
    I = solution[3,:]
    Q = solution[4,:]
    R = solution[5,:]

    plot(solution.t,R; dpi = 600, frame = :box) #plot(x,y)
    #savefig(joinpath("figures", "simulation_$(r0).png"))
    # TODO save data to a CSV file
#end


# ===> NOT WORKING RN ===>
# TODO C according to time -> as a function C(t) ?

# Estimate the number of cases C(t) 
# Formula: CT(t) ≡ (Ss + sfTTR)[A(t)/Tr]
#C => [(Ss + s*ft*Tr) * (A/Tr) ] # C as vector with [] ?

C = cumsum((Ss + _s*_ft*_Tr) * (A/_Tr))

plot(solution.t,C)


# Try with different R0 values --> Figure 2A



#Figure 2B 
# Calculate Reff
#TODO: replace with parameters


Reff = r0 / (1 + Ss + s * ft * Tr)

# Find the R0 used for the numerical solution(eq 1-5)

# Parameters
R₀ = 

# Expression 10 
A = (Reff/1-Reff)*(NE₀*Tr/R₀)
A = ((NE₀*Tr)/(1-Reff))*1/(1+Ss+s*ft*Tr)

