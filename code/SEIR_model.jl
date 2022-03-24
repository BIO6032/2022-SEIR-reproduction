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