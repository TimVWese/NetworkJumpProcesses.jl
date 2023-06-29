# SIS

In this example, we simulate the spreading of a SIS epidemic on a contact network. Each node in the network represents a person and can be in one of two states: susceptible (`0`) or infective (`1`).
Infective nodes can infect their susceptible neighbours at a rate `β` and recover at a rate `γ`.

## Import the required packages

For this example, we will use the minimal packages `NetworkJumpProcesses`, `Graphs`, and `JumpProcesses`.
Additionally, we will use the `Plots` package to visualise the results.
The `Graphs` package provides the graphs on which the jump process will be defined.
In this case this is the `erdos_renyi` function, which generates a random graph according to the [Erdős–Rényi](https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93R%C3%A9nyi_model) model.
The `NetworkJumpProcesses` package provides the `network_jump_set` function and network jump types for defining the jump set to use in the simulation algorithms of the The `JumpProcesses` package.
Additionally, we will use the `Plots` to visualise the results.
Import the packages as as usual, or install them first using `] add NetworkJumpProcesses Graphs JumpProcesses Plots`, if not done yet.

```julia
using NetworkJumpProcesses
using Graphs
using JumpProcesses
using Plots: plot
```

## Defining the contact network

The contact network is defined using the `erdos_renyi` function. The function takes as input the number of nodes `n` and the probability `p` of an edge existing between any two nodes. In this example, we set `n = 100` and `p = 0.1`.

```julia
n = 100
g = erdos_renyi(n, 0.1)
```

## Defining the SIS model

The SIS model is defined using the `network_jump_set` function from the `JumpProcesses` package. The function takes as input the contact network `gg` and a list of vertex reactions and edge reactions.
In this example, we define one vertex reaction, the curing event, and one edge reaction, the infection event.
Each jump takes two functions, a `rate` and `effect!`.
In the case of the vertex reaction, the `rate` function takes as input the state `v`, the neighbours `nghbs`, the parameters `p`, and the time `t` and returns the rate of the reaction.
The `effect!` function takes as input the state `v`, the neighbours `nghbs`, the parameters `p`, and the time `t` and manipulates the state `v` to reflect the effect of the reaction.
The same holds for the edge, except that the functions take as input the source state `vs` and the destination state `vd`, instead of the vertex and its neighbours state.

```julia
v_IS = ConstantJumpVertex(
    (v, nghbs, p, t) -> v[1] == 1 ? p[2] : 0.0,
    (v, nghbs, p, t) -> v[1] = 0
)

e_SI = ConstantJumpEdge(
    (vs, vd, p, t) -> vs[1] == 1 ? p[1] : 0.0,
    (vs, vd, p, t) -> vd[1] = 1
)
```

Having defined the jumps, the jump set can be defined as follows:

```julia
jset = network_jump_set(gg; vertex_reactions=[v_IS], edge_reactions=[e_SI])
```

## Defining the problem

The `JumpProblem` constructor defines the jump problem.
It requires `DiscreteProblem` with initial state `u₀`, the time horizon `(t0, tf)`, and the parameters `p`, but without any other dynamics.
In this case the initial condition is such that all nodes are susceptible, except for four nodes, which are infective.
The `JumpProblem` constructor also takes as input the problem an [aggregator](https://docs.sciml.ai/JumpProcesses/stable/jump_types/#Jump-Aggregators-for-Exact-Simulation), and the jump set `jset`.
In this example, we use `Direct` aggregation and the `SSAStepper` [solver](https://docs.sciml.ai/JumpProcesses/stable/jump_solve/).

```julia
p = (0.1, 0.08)
u₀ = zeros(Int64, n)
u₀[1:4] .= 1

dprob = DiscreteProblem(u₀, (0, 40.0), p) 
jprob = JumpProblem(dprob, Direct(), jset)
sol = solve(jprob, SSAStepper())
```

## Plotting the results

We can now plot the number of infected agents over time.

```julia
using Plots

plot(sol.t, [count(sol[i].==1) for i in eachindex(sol.t)],
     xlabel="t", ylabel="I")
```
