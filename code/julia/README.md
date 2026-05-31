# AlzheimerModel Useage

## Loading

```
julia> include("AD_ODE.jl")
```

## Basic Workflow

```
julia> sys  = AlzheimerModel.build_system()
julia> prob = AlzheimerModel.make_problem(sys)
julia> sol  = AlzheimerModel.solve_model(prob)

julia> AlzheimerModel.plot_solution(sol,sys)
julia> AlzheimerModel.print_final(sol, sys)
```

## Trying Different Timespans

```
julia> prob2 = AlzheimerModel.make_problem(sys, tspan=(0.0, 200.0))
julia> sol2  = AlzheimerModel.solve_model(prob2)
```

## Tweak Parameters Without Rebuilding

```
julia> prob3 = remake(prob, p=[sys.k1 => 0.5])
julia> sol3  = AlzheimerModel.solve_model(prob3)
```
