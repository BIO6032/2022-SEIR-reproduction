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
    D(R) ~ fr * Q + fr * I + fr * A
]

# Problem 
@named seir = ODESystem(equations)
S0 = 9999.
A0 = 1.
I0 = 0.
R0 = 0.
Q0 = 0.
u0 = [S => S0, A => A0, I => I0, R => R0, Q => Q0]
tspan = (0.0, 100.0)


# Simulations
plot()
for r0 in 2:8
    # TODO formula for parameters from r0
    #r0 = 8.

    # Get Parameters values
    _fr = 1. / 14.
    #_β = r0 * _fr / _N
    _N = S0+A0+I0+R0+Q0
    _e = 4. / _N 
    _Ss = 0.6 
    _Tr = 14. 
    _fs = _Ss * _fr
    _s = 0.9
    _ft = (2. /7.)*(1. /_Tr)
    _fq = 0.9
    _β = (r0*(_fs+_fr+_s*_ft))/_N


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
  
    M = solution[1:end,:] # Matrix of results 5 compartments x 31 time steps
    
    # Assign values to each compartment over time
    S = solution[1,:]
    A = solution[2,:]
    I = solution[3,:]
    Q = solution[4,:]
    R = solution[5,:]

    
    C = cumsum((_Ss + _s*_ft*_Tr) * (A/_Tr))
    plot!(solution.t,C, lab=r0)
end
xaxis!(tspan...)

# Try with different R0 values --> Figure 2A



#Figure 2B 
# Calculate Reff
#TODO: replace with parameters


Reff = r0 / (1 + Ss + s * ft * Tr)

# Find the R0 used for the numerical solution(eq 1-5)

# Parameters
R₀ = β*N/fr

# Expression 10 
A = (Reff/1-Reff)*(NE₀*Tr/R₀)
A = ((NE₀*Tr)/(1-Reff))*1/(1+Ss+s*ft*Tr)

