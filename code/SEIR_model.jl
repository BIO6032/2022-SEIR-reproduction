# Simulations
using DifferentialEquations
using OrdinaryDiffEq
#Permet de résoudre des équations différentielles 
using ModelingToolkit
#Permet de résoudre de créer des problèmes et génèrer des solutions à partir des équations imputées au modèle
# Figures
using Plots
# Data
using DataFrames
using CSV

# Parameters 
@parameters t β e fs fr s ft fq
#Les paramètres sont les valeurs qui ne changent pas au cours du temps, ce sont des taux.
#Beta = Taux de transmission 
#e = taux exogène 
#fs = taux de conversion de asymptomatique en symptomatique
#ft = taux de dépistage
#fr = taux auquel les individus deviennent recovered
#s = sensibilité du test de dépistage 
#fq = taux auquel les infectés symptomatiques et asymptomatiques sont placés en quarantaine

# Variables 
@variables S(t) A(t) I(t) R(t) Q(t)
# S(t) représente le compartiment des individus susceptibles sur le campus. 
# A(t) répresente le compartiment des individus infectieux et asymptomatiques.
# I(t) représente le compartiment des individus infectieux et symptomatiques.
# R(t) répresente le compartiment des individus qui ont récupéré de l'infection du virus.
# Q(t) représente le compartiment des individus qui ont été placées en quarantaine. 

# Differential
D = Differential(t)
# Les équations démontrent des variables qui changent selon le temps.
# Equations
equations = [
    # les susceptibles reçoivent un apport B*S*A de transmission des asymptomatiques 
    # moins B*S*I perdues vers les infectueux et e*S d'exogenou
    D(S) ~ -β * S * A - β * S * I - e * S,
    # les asymptotiques reçoivent un taux B*S*A des symptomatiques moins le taux symptomatique perdues
    # aux infectieux, un testing rate perdues aux Quarantine et un recovery rate fr perdues aux Recovery   
    D(A) ~ β * S * A + β * S * I - fs * A - fr * A - s * ft * A + e * S,
    # les infectueux reçoivent un taux symptomatiques fs des asymptomatiques et perdent un recovery rate fr et un recovery rate fq
    D(I) ~ fs * A - fr * I - fq * I,
    # les quarantined reçoivent un quarantine rate fq des infectueux et un testing rate des asymptomatiques et perdent un recovery rate
    # fq vers les quarantinedd      
    D(Q) ~ fq * I + s * ft * A - fr * Q,
    # les recovered reçoivent un recover rate des quarantined, des infectueux et des asymptotiques
    D(R) ~ fr * Q + fr * I + fr * A,
]

# Problem 
# On créé un système d'équations qui regroupe les dérivées pour chaque compartiments
@named seir = ODESystem(equations)
# On définit les valeurs initiales pour chaque variables.
S0 = 9999.0
A0 = 1.0
I0 = 0.0
R0 = 0.0
Q0 = 0.0
# On assigne à chaque compartiment/variable sa valeur initiale dans le vecteur u0.
u0 = [S => S0, A => A0, I => I0, R => R0, Q => Q0]
# On définit l'intervalle de temps (0 à 100 pas de temps)
tspan = (0.0, 100.0)

# Simulations
# On génère une fenêtre de plot ne présentant que des axes ("vide")
plot()
r0_values = collect(LinRange(0.1, 8.0, 200))
r_eff = similar(r0_values)
c_eff = similar(r0_values)
# On génère un range de valeurs pour r0 allant de 0.1 à 2 (ou autres) et prenant 5 valeurs (ou autres).
# Fonction similar permet de créer un tableau de valeurs pour un ou plusieurs éléments donnés (ici, le R0), 
# permet ainsi d'obtenir les valeurs de R eff et C eff pour les valeurs données de R0 

#La boucle qui suit a pour fonction d'effectuer les simulations du modèle définit par les équation seir à partir de différentes 
#valeurs de R_0. Cette boucle commence par définir les valeurs de R_0 qui seront utilisées pour effectuer les simulations. Ensuite, 
# les valeurs des paramètres sont définies et les paramètres sont associées à leurs différentes valeurs. 
# On définit et résout le système d'équations seir et on assigne les valeurs obtenus pour les différents compartiments à travers le temps.
# Finalement, on trace le graphique des résultats des simluations pour les différentes valeurs de R_0 (nombre de cas positifs en fonction du temps)
for (i, r0) in enumerate(r0_values)
    # La boucle utilise i iteration (partant de 1) et r0 prend successivemnet la ith valeurs de r0_values (enumerate)

    # TODO formula for parameters from r0 -> unavailable formula and unable to recover it

    # Get Parameter values - on définit les valeurs des paramètres
    _fr = 1.0 / 14.0
    _N = S0 + A0 + I0 + R0 + Q0
    _e = 4.0 / _N
    _Ss = 0.6
    _Tr = 14.0
    _fs = _Ss * _fr
    _s = 0.9
    _ft = (2.0 / 7.0)
    _fq = 0.9
    #_β = r0 * _fr / _N
    _β = (r0 * (_fs + _fr + _s * _ft)) / _N

    # On associe les paramètres aux valeurs que nous avons déterminés
    p = [
        fr => _fr,
        β => _β,
        e => _e,
        fs => _fs,
        s => _s,
        ft => _ft,
        fq => _fq,
    ]

    # Définit le problème qu'on veut résoudre en fonction du système d'équations, les valeurs initiales des paramètres, la durée des simulations
    # ainsi que les paramètres. 
    problem = ODEProblem(seir, u0, tspan, p)

    # Résout le problème définit à la ligne précédente
    solution = solve(problem)
    # Matrix of results 5 compartments * x time steps; crée une matrice qui va être utilisée pour stocker les valeurs obtenus des variables au cours du temps
    M = solution[1:end, :] 
    # Assign values to each compartment over time - stockage des valeurs des variables dans chaque compartiment à partir de la solution
    S = solution[1, :] # Dans l'objet solution prendre le 1er élément à tous les pas de temps
    A = solution[2, :]
    I = solution[3, :]
    Q = solution[4, :]
    R = solution[5, :]

    #C = cumsum((_Ss + _s*_ft*_Tr) * (A/_Tr))
    #plot!(solution.t,C, lab=r0)
    vecN = repeat([_N], outer = length(S)) # crée un vecteur N qui a la même longueur que le vecteur du compartiment S et répète la valeur de N 
    # afin de pouvoir faire le calcul qui suit

    C_alt = vecN - S # formule qui calcule le nombre de cas positifs en fonction du temps
    plot!(solution.t, C_alt, lab = "", line_z = r0, clim = (0.0, 2.0), c = :roma) # plot des différentes

    maxA = maximum(A)

    r_eff[i] = r0 / (1 + _Ss + _s * _ft * _Tr)

    c_eff[i] = (_Ss + _s * _ft * _Tr) * (maxA / _Tr)

end
xaxis!(tspan...)

plot(r_eff, c_eff)

# Figure 2B

R₀ = LinRange(0.1, 8.0, 200)
x = fill(NaN, length(R₀))
y = fill(NaN, length(R₀))

for (i, _r₀) in enumerate(R₀)
    # paramètres
    _NE₀ = 4.0
    _Ss = 0.6
    _s = 0.9
    _ft = (2.0 / 7.0)
    _Tr = 14.0

    Reff = _r₀ ./ (1 + _Ss + _s * _ft * _Tr)
    A = (Reff ./ (1.0 .- Reff)) .* (_NE₀ * _Tr ./ _r₀)

    C = (_Ss + _s * _ft * _Tr) * (A / _Tr)

    if C < 0
        break
    end

    x[i] = Reff
    y[i] = C

end

_last_valid_idx = findfirst(isnan, y)
_last_valid_idx = isnothing(_last_valid_idx) ? length(y) : _last_valid_idx
x = x[1:_last_valid_idx]
y = y[1:_last_valid_idx]

plot(x, y, legend = :bottomright, frame = :box, c = :red, lab = "Asymptotic expression")
plot!(r_eff, c_eff, c = :blue, lab = "Numerical solution")
# Try with different R0 values --> Figure 2A
# Find the R0 used for the numerical solution(eq 1-5) -> unavailable formula and unable to recover it


### Figure 3

# Importer les données
df = DataFrame(CSV.File(joinpath("data", "suffolk_county_data.csv")))

#Figure 3A pour les cas cumulatifs 
plot(
    df[:, 1],
    df[:, 5],
    frame = :box,
    lab = "cummulative cases",
    xaxis = "date",
    yaxis = "Positive cases",
    color = :black,
    linestyle = :solid,
)

#Figures 3B pour les daily cases 
plot(
    df[:, 1],
    df[:, 3],
    frame = :box,
    lab = "daily cases",
    xaxis = "date",
    yaxis = "Daily cases",
    color = :black,
    linestyle = :dash,
)

# ecrire equation #7 + p(t) + R contact tracing pour 3b
#p = P₀ * δ * (t-tp) mais on considère que p=1 ici

### Système alternatif avec D(A) = equation 7
@parameters t β e fs fr s ft fq r₀ Tr Reff NE₀ Ss
@variables S(t) A(t) I(t) R(t) Q(t)

D = Differential(t)

equations_alt = [
    D(S) ~ -β * S * A - β * S * I - e * S,
    D(A) ~ r₀ / Tr * (1 - (1 / Reff)) * A + NE₀,
    D(I) ~ fs * A - fr * I - fq * I,
    D(Q) ~ fq * I + s * ft * A - fr * Q,
    D(R) ~ fr * Q + fr * I + fr * A,
]

@named seir_alt = ODESystem(equations_alt)
S0 = 9999.0
A0 = 1.0
I0 = 0.0
R0 = 0.0
Q0 = 0.0
u0 = [S => S0, A => A0, I => I0, R => R0, Q => Q0]

## avec r₀=2.5
_fr = 1.0 / 14.0
_N = S0 + A0 + I0 + R0 + Q0
_e = 4.0 / _N
_Ss = 0.6
_Tr = 14.0
_fs = _Ss * _fr
_s = 0.9
_ft = (2.0 / 7.0)
_fq = 0.9
_r₀ = 2.5
_β = (_r₀ * (_fs + _fr + _s * _ft)) / _N
_Reff = _r₀ ./ (1 + _Ss + _s * _ft * _Tr)
_NE₀ = 4.0

p_alt = [
    fr => _fr,
    β => _β,
    e => _e,
    fs => _fs,
    Ss => _Ss,
    s => _s,
    ft => _ft,
    fq => _fq,
    r₀ => _r₀,
    NE₀ => _NE₀,
    Tr => _Tr,
    Reff => _Reff,
]

problem_alt = ODEProblem(seir_alt, u0, tspan, p_alt)

# Résout le problème définit à la ligne précédente
solution_alt = solve(problem_alt)

# Assign values to each compartment over time - stockage des valeurs des variables dans chaque compartiment à partir de la solution
S = solution_alt[1, :] # Dans l'objet solution prendre le 1er élément à tous les pas de temps
A = solution_alt[2, :]
I = solution_alt[3, :]
Q = solution_alt[4, :]
R = solution_alt[5, :]

vecN_alt = repeat([_N], outer = length(S)) 
C = vecN_alt - S # formule qui calcule le nombre de cas positifs en fonction du temps

# Positive cases fig3a (r₀=2.5 blue dash)
plot(solution_alt.t, C,      yaxis = "Positive cases (cum)",
lab = "r₀=2.5",color=:blue, linestyle=:dash) # cas cumulatifs au cours du temps


# Daily cases fig3b (r₀=2.5 black)
plot(solution_alt.t, ... , color=:black, linestyle=:solid)


# TODO: superposer les plots des fig 3