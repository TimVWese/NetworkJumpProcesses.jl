# Reaction Diffusion

In this example, we simulate the diffusion and reaction of two chemical elements, A and B, which might react with each other and form C. C is disposed of in some sink nodes. Over the edges there happens diffusion.

## Defining the contact network

The contact network is defined using the `erdos_renyi` function from the `Graphs` package. The function takes as input the number of nodes `n` and the probability of an edge between two nodes `p`.

```julia
using Graphs

n = 12
gg = erdos_renyi(n, 0.4)
```

## Defining the vertex reactions

The vertex reactions are defined using the `ConstantJumpVertex`(@ref) constructor.
The function takes as input two anonymous functions: the first function defines the rate of the reaction, and the second function defines the update of the state.
In this example, we define three vertex reactions: A_source, B_source, and C_sink.
In the first two, respectively, A and B are produced at a constant rate, while in the third one C is disposed of.

```julia
using JumpProcesses

A_source = ConstantJumpVertex(
    (v, nghbs, p, t) -> p.rA,
    (v, nghbs, p, t) -> v[1] += 1
)

B_source = ConstantJumpVertex(
    (v, nghbs, p, t) -> p.rB,
    (v, nghbs, p, t) -> v[2] += 1
)

C_sink = ConstantJumpVertex(
    (v, nghbs, p, t) -> p.dC*v[3],
    (v, nghbs, p, t) -> v[3] -= 1
)
```

## Defining the edge reactions

The edge reactions are defined using the `ConstantJumpEdge`(@ref) function.
Once again; two functions define the rate and affect the state.
Note that the function signature sllightly differs from the vertex reactions.
It is assumed that the diffusion is the same for all components, so it can be
defined as follows.

```julia
diffusion = i -> ConstantJumpEdge(
    (vs, vd, p, t) -> max(p.D*(vs[i] - vd[i]), 0.),
    (vs, vd, p, t) -> begin
        vs[i] -= 1
        vd[i] += 1
    end
)

diffusion_A = diffusion(1)
diffusion_B = diffusion(2)
diffusion_C = diffusion(3)
```

## Defining the network jump set

Since we make use heterogenous vertex jumps, we need to define an array of
jumps for each node.
Component C can be created in any node, and it can be disposed in three of them.
A and B are also generated in 3 nodes.

```julia
A_source_reactions = [[A_source, C_creation] for i in 1:3]
B_source_reactions = [[B_source, C_creation] for i in 1:3]
C_sink_reactions = [[C_sink, C_creation] for i in 1:3]
v_reactions = vcat(A_source_reactions, B_source_reactions, C_sink_reactions)
append!(v_reactions, [[C_creation] for i in 1:(n-length(v_reactions))])
shuffle!(v_reactions)
```

The edge jumps are homogeneous for all edges, so can be defined as a single array.
Having the two sets of jumps, we can define the network jump set.

```julia
e_reactions = [diffusion_A, diffusion_B, diffusion_C]
jset = network_jump_set(gg; vertex_reactions=v_reactions,
                        edge_reactions=e_reactions, nb_states=3)
```

## Defining the problem

A discrete problem is used with initial zero state and parameter tuple `p`.
We first use the `Direct` aggregator.

```julia
p = (
    rA = 0.1,
    rB = 0.1,
    rC = 0.2,
    dC = 1.5,
    D = 0.1
)

T = 100.0
u₀ = zeros(Int64, 3*n)
dprob = DiscreteProblem(u₀, (0., T), p) 
direct_jump_prob = JumpProblem(dprob, Direct(), jset)
sol = solve(direct_jump_prob, SSAStepper())
```

## Plotting the results

The amount of each component can be computed and plotted over time.

```julia
using Plots

get_concentrations(x::Vector{Int64}) = [
    sum(x[i] for i in 1:3:length(x)),
    sum(x[i] for i in 2:3:length(x)),
    sum(x[i] for i in 3:3:length(x))
]

get_concentrations(x::ODESolution) = reduce(hcat, [get_concentrations(x[t]) for t in eachindex(x)])';

plot(sol.t, get_concentrations(sol), label=["A" "B" "C"], xlabel="t", ylabel="concentration", lw=2, legend=:topleft)
```

## Using dependency graphs

It is also possible to use more performant aggregators that require
dependency graphs, such as [RSSA]([Title](https://docs.sciml.ai/JumpProcesses/stable/jump_types/#Jump-Aggregators-Requiring-Dependency-Graphs)).
this package provides two functions to generate the dependency graphs.
Ihey require the underlying graph, the number of reactions per vertex, numer of reactions per edge and the number of states.
The first two numbers can be either a single integer or an array of integers, depending on wheter the amounts are homogeneous over the vertices and edges or not.
In this example, we use the same number of reactions for all edges, but not for the vertices.
Hence we use the following code.

```julia
nb_vertex_reacs = [length(v) for v in v_reactions]
vtj = vartojumps(gg, nb_vertex_reacs, 3, 3)
jtv = jumptovars(gg, nb_vertex_reacs, 3, 3)

rssa_jump_prob = JumpProblem(dprob, RSSA(), jset; vartojumps_map=vtj, jumptovars_map=jtv)
sol = solve(rssa_jump_prob, SSAStepper())

plot(sol.t, get_concentrations(sol), label=["A" "B" "C"], xlabel="t", ylabel="concentration", lw=2, legend=:topleft)
```
