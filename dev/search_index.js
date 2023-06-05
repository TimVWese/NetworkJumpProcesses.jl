var documenterSearchIndex = {"docs":
[{"location":"#NetworkJumpProcesses.jl","page":"NetworkJumpProcesses.jl","title":"NetworkJumpProcesses.jl","text":"","category":"section"},{"location":"","page":"NetworkJumpProcesses.jl","title":"NetworkJumpProcesses.jl","text":"Documentation for NetworkJumpProcesses.jl","category":"page"},{"location":"","page":"NetworkJumpProcesses.jl","title":"NetworkJumpProcesses.jl","text":"Modules = [NetworkJumpProcesses]","category":"page"},{"location":"#NetworkJumpProcesses.ConstantJumpEdge","page":"NetworkJumpProcesses.jl","title":"NetworkJumpProcesses.ConstantJumpEdge","text":"ConstantJumpEdge(rate, affect!)\n\nConstruct a constant rate jump over an edge with rate rate and affect! methods.\n\nArguments\n\nrate::Function: signature (vs, vd, p, t) -> Real\naffect!::Function: signature (vs, vd, nghbs, p, t) -> nothing\n\nSee also: Types of Jumps\n\n\n\n\n\n","category":"type"},{"location":"#NetworkJumpProcesses.ConstantJumpVertex","page":"NetworkJumpProcesses.jl","title":"NetworkJumpProcesses.ConstantJumpVertex","text":"ConstantJumpVertex(rate, affect!)\n\nConstruct a constant rate jump vertex with rate rate and affect! methods.\n\nArguments\n\nrate::Function: signature (v, nhgbs, p, t) -> Real\naffect!::Function: signature (v, nghbs, p, t)\n\nSee also: Types of Jumps\n\n\n\n\n\n","category":"type"},{"location":"#NetworkJumpProcesses.VariableJumpEdge","page":"NetworkJumpProcesses.jl","title":"NetworkJumpProcesses.VariableJumpEdge","text":"VariableJumpEdge(rate, affect!)\n\nConstruct a variable rate jump over an edge with rate rate and affect! methods.\n\nArguments\n\nrate::Function: signature (vs, vd, p, t) -> Real\naffect!::Function: signature (vs, vd, nghbs, p, t) -> nothing\n\nSee also: Types of Jumps\n\n\n\n\n\n","category":"type"},{"location":"#NetworkJumpProcesses.VariableJumpVertex","page":"NetworkJumpProcesses.jl","title":"NetworkJumpProcesses.VariableJumpVertex","text":"VariableJumpVertex(rate, affect!)\n\nConstruct a variable rate jump vertex with rate rate and affect! methods.\n\nArguments\n\nrate::Function: signature (v, nhgbs, p, t) -> Real\naffect!::Function: signature (v, nghbs, p, t) -> nothing\n\nSee also: Types of Jumps\n\n\n\n\n\n","category":"type"},{"location":"#NetworkJumpProcesses.jumptovars","page":"NetworkJumpProcesses.jl","title":"NetworkJumpProcesses.jumptovars","text":"jumptovars(graph, nb_vertex_reacs, nb_edge_reacs, nb_vertex_states=1)\n\nCreate a dependency graph that maps the vertex and edge reactions to the vertex states. This graph can be used for the RSSA and RSSACR aggregators.\n\nSee also: vartojumps, Jump Aggregators Requiring Dependency Graphs \n\n\n\n\n\n","category":"function"},{"location":"#NetworkJumpProcesses.network_jump_set-Union{Tuple{Any}, Tuple{U}, Tuple{T}} where {T<:Union{JumpVertex, Vector{var\"#s37\"} where var\"#s37\"<:JumpVertex}, U<:Union{JumpEdge, Vector{var\"#s38\"} where var\"#s38\"<:JumpEdge}}","page":"NetworkJumpProcesses.jl","title":"NetworkJumpProcesses.network_jump_set","text":"network_jump_set(\n    graph; vertex_reactions::Vector{T}=Vector{JumpVertex}(),\n    edge_reactions::Vector{U}=Vector{JumpEdge}(), nb_states=1\n) where {\n    T <: Union{JumpVertex, Vector{<:JumpVertex}},\n    U <: Union{JumpEdge, Vector{<:JumpEdge}}\n}\n\nConstruct a JumpSet from a Graph and a list of JumpVertex and JumpEdge reactions. vertex_reactions and edge_reactions can be either an vector of reactions which will all be applied to every vertex and edge respectively. The other option is a vector of vectors of reactions, where the ith vector of reactions will be applied to the ith vertex or edge. These variables may be mixed. Each vertex has nb_states variables associated.\n\nSee also: ConstantJumpVertex, ConstantJumpEdge, VariableJumpVertex, VariableJumpEdge\n\n\n\n\n\n","category":"method"},{"location":"#NetworkJumpProcesses.vartojumps-Union{Tuple{U}, Tuple{T}, Tuple{Graphs.AbstractGraph, T, U}, Tuple{Graphs.AbstractGraph, T, U, Any}} where {T<:Union{Integer, Vector{var\"#s9\"} where var\"#s9\"<:Integer}, U<:Union{Integer, Vector{var\"#s8\"} where var\"#s8\"<:Integer}}","page":"NetworkJumpProcesses.jl","title":"NetworkJumpProcesses.vartojumps","text":"vartojumps(graph::AbstractGraph, nb_vertex_reacs::T, nb_edge_reacs::U, nb_vertex_states=1) where {\n    T <: Union{Integer, Vector{<:Integer}},\n    U <: Union{Integer, Vector{<:Integer}}\n}\n\nCreate a dependency graph that maps the vertex states to the vertex and edge reactions. This graph can be used for the RSSA and RSSACR aggregators. If all vertices (edges) have the same number of reactions, then nb_vertex_reacs (nb_edge_reacs) can be an integer, equal to that number. Otherwise, nb_vertex_reacs (nb_edge_reacs) must be a vector of length nv(graph) (ne(graph)) , such that nb_vertex_reacs[v] (nb_edge_reacs[e]) is the number of reactions associated with vertex v (edge e).\n\nSee also: jumptovars, Jump Aggregators Requiring Dependency Graphs \n\n\n\n\n\n","category":"method"},{"location":"#NetworkJumpProcesses.vertex_range-Tuple{Any, Any}","page":"NetworkJumpProcesses.jl","title":"NetworkJumpProcesses.vertex_range","text":"vertex_range(n, v)\n\nReturns the range of vertices for a given vertex v, given that the number of states is n.\n\n\n\n\n\n","category":"method"},{"location":"#NetworkJumpProcesses.vertex_to_edges-Tuple{Graphs.AbstractGraph}","page":"NetworkJumpProcesses.jl","title":"NetworkJumpProcesses.vertex_to_edges","text":"vertex_to_edges(graph::AbstractGraph)\n\nCreate a dictionary that maps each vertex to the edges it is connected to.\n\n\n\n\n\n","category":"method"}]
}
