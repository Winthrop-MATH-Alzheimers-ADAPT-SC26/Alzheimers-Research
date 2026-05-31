# AlzheimerModel Useage

## Loading

```julia-repl
julia> include("AD_ODE.jl")
julia> using .AlzheimerModel
julia> ModelingToolkit, DifferentialEquations
```

## Basic Workflow

```julia-repl
julia> sys  = build_system()
julia> prob = make_problem(sys)
julia> sol  = solve_model(prob)

julia> plot_solution(sol,sys)
julia> print_final(sol, sys)
```

## Trying Different Timespans

```julia-repl
julia> prob2 = make_problem(sys, tspan=(0.0, 200.0))
julia> sol2  = solve_model(prob2)
```

## Tweak Parameters Without Rebuilding

```julia-repl
julia> prob3 = remake(prob, p=[sys.k1 => 0.5])
julia> sol3  = solve_model(prob3)
```
## Try Different Solvers

```julia-repl
julia> sol4 = solve_model(prob3, solver=Rodas5P())
```

If solver is specified, then underlying `solve()` will determine best solver to use.
