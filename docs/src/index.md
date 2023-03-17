# NetworkJumpProcesses.jl

Documentation for NetworkJumpProcesses.jl

```@docs
ConstantJumpVertex(rate, affect!)
VariableJumpVertex(rate, affect!)
ConstantJumpEdge(rate, affect!)
VariableJumpEdge(rate, affect!)
network_jump_set(graph; vertex_reactions::Vector{T}=Vector{NetworkJumpProcesses.JumpVertex}(), edge_reactions::Vector{U}=Vector{NetworkJumpProcesses.JumpEdge}()) where {T <: NetworkJumpProcesses.JumpVertex, U <: NetworkJumpProcesses.JumpEdge}
vartojumps(graph, nb_vertex_states, nb_vertex_reacs, nb_edge_reacs)
jumptovars(graph, nb_vertex_states, nb_vertex_reacs, nb_edge_reacs)
```
