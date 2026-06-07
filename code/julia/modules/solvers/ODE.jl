module ODE

export print_final

using DifferentialEquations
using ModelingToolkit

# --- Functions --- #

function print_final(sol)
    println("Final values at t = ", sol.t[end])
    println("Aβ = ", sol[sys.Aβ][end])
    println("Ca = ", sol[sys.Ca][end])
    println("τ  = ", sol[sys.τ][end])
    println("N  = ", sol[sys.N][end])
    println("C  = ", sol[sys.C][end])
end

end
