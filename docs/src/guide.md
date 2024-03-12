# Package Guide

The idea is to define some jumps (discrete state, continuous time, stochastic variables) for the vertices and for edges of an arbitrary `Graph`, and the package functionality returns a [`JumpSet`](https://docs.sciml.ai/JumpProcesses/stable/api/#JumpProcesses.JumpSet) which can be solved with `JumpProcesses.jl` and `DifferentialEquations.jl`.

The important functions are:
1. `network_jump_set`: creates the `JumpSet` given a `Graph`, the jumps defined for vertices (`vertex_reactions`), and the jumps defined for edges, (`edge_reactions`). These last two can be each a vector of jumps, in which case all vertices or edges behave according to to the same reactions. An other option is that one or both is a vector of vectors. In this case the outer vector should have the same number of elements as vertices (`vertex_reactions`) or edges (`edge_reactions`) in the graph.
2. `dependency_graph, vartojumps, jumptovars`: these create the [dependency graphs necessary for some aggregators](https://docs.sciml.ai/JumpProcesses/stable/jump_types/#Jump-Aggregators-Requiring-Dependency-Graphs).

As in `JumpProcesses.jl`, the jumps are defined by a `rate` and `affect!` function.
The `rate` returns a single number, while `affect!` modifies the arguments directly.
In the case of a vertex the input arguments for both are
* `v`: state of the vertex itself
* `nhgbs`: state of the neighbours
* `p`: model parameter values
* `t`: current time step
In case of an edge the arguments are `(vs, vd, p, t)` where `vs` resp. `vd` are the source and destination vertex states.
