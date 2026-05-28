import DifferentialEquations as DE
import ModelingToolkit as MTK
import Plots
import ModelingToolkit: t_nounits as t, D_nounits as D, @variables, @parameters, @named, @mtkcompile, mtkcompile

# Define state variables and their intial conditions
@variables Aβ(t)=0.01 Ca(t)=0 τ(t)=0 N(t)=0 C(t)=0

# Define parameters
@parameters a1=0.25 a2=2.89
@parameters b1=11 b2=100
@parameters c1=0.25 c2=0.15 c3=0.001
@parameters d1=0.01 d2=1
@parameters e1=1.67 e2=3.83 e3=1
@parameters k1=0.35 k2=0.01 k3=0.000708 k4=1.00 k5=1.00
@parameters u1=0 u2=0 u3=0
@parameters σ=0.00027 R=1

# Define differential equations
eqs = [D(Aβ) ~ a1 + (a2 * (Ca / (Ca + σ))) - (k1 * Aβ) - (u1 * Aβ)
D(Ca) ~ b1 + (b2 * Aβ) - (k3 * Ca) - (u2 * Ca)
D(τ) ~ c1 + (c2 * Aβ) + (c3 * Ca) - (k3 * τ) - (u3 * τ)
D(N) ~ d1 + (d2 * τ) - (k4 * N)
D(C) ~ e1 + (e2 * N * R) + (e3 * τ) - (k5 * C)]

# Create ODESystem
@mtkcompile sys = MTK.ODESystem(eqs, t)

# Convert from symbolic to numeric for simulation
tspan = (0.0, 100.0)
prob = DE.ODEProblem(sys, [], tspan)

# Solve ODE
println("Solving ODE...")
sol = DE.solve(prob)
println("Done!")

# Plot Solution
p1 = Plots.plot(sol, idxs = [Aβ], title = "Aβ")
p2 = Plots.plot(sol, idxs = [Ca], title = "Ca")
p3 = Plots.plot(sol, idxs = [τ], title = "τ")
p4 = Plots.plot(sol, idxs = [N], title = "N")
p5 = Plots.plot(sol, idxs = [C], title = "C")
p = Plots.plot(p1, p2, p3, p4, p5, layout = (2, 3))
Plots.gui(p)

println("Aβ = ", sol[Aβ][end])
println("Ca = ", sol[Ca][end])
println("τ  = ", sol[τ][end])
println("N  = ", sol[N][end])
println("C  = ", sol[C][end])
println()
#println("Press Enter to Exit")
readline()
