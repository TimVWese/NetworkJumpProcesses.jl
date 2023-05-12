module NetworkJumpProcesses

using Graphs
using JumpProcesses
using Logging

export JumpVertex, JumpEdge
export ConstantJumpVertex, ConstantJumpEdge
export VariableJumpVertex, VariableJumpEdge
export network_jump_set
export dependency_graph, vartojumps, jumptovars

abstract type JumpElement end
abstract type JumpVertex <: JumpElement end

"""
    ConstantJumpVertex(rate, affect!)

Construct a constant rate jump vertex with rate `rate` and `affect!` methods.

# Arguments
- `rate::Function`: signature `(v, nhgbs, p, t) -> Real`
- `affect!::Function`: signature `(v, nghbs, p, t)`

See also: [Types of Jumps](https://docs.sciml.ai/JumpProcesses/stable/jump_types/#Types-of-Jumps:-Constant-Rate,-Mass-Action,-Variable-Rate-and-Regular)
"""
Base.@kwdef struct ConstantJumpVertex{T, U} <:JumpVertex
    rate::T # signature (v, nhgbs, p, t) -> Real
    affect!::U # signature (vn, v, nghbs, p, t) -> nothing 
end

"""
    VariableJumpVertex(rate, affect!)

Construct a variable rate jump vertex with rate `rate` and `affect!` methods.

# Arguments
- `rate::Function`: signature `(v, nhgbs, p, t) -> Real`
- `affect!::Function`: signature `(v, nghbs, p, t) -> nothing`

See also: [Types of Jumps](https://docs.sciml.ai/JumpProcesses/stable/jump_types/#Types-of-Jumps:-Constant-Rate,-Mass-Action,-Variable-Rate-and-Regular)
"""
Base.@kwdef struct VariableJumpVertex{T, U} <:JumpVertex
    rate::T # signature (v, nhgbs, p, t) -> Real
    affect!::U # signature (vn, v, nghbs, p, t) -> nothing
end

abstract type JumpEdge <: JumpElement end

"""
    ConstantJumpEdge(rate, affect!)

Construct a constant rate jump over an edge with rate `rate` and `affect!` methods.

# Arguments
- `rate::Function`: signature `(vs, vd, p, t) -> Real`
- `affect!::Function`: signature `(vs, vd, nghbs, p, t) -> nothing`

See also: [Types of Jumps](https://docs.sciml.ai/JumpProcesses/stable/jump_types/#Types-of-Jumps:-Constant-Rate,-Mass-Action,-Variable-Rate-and-Regular)
"""
Base.@kwdef struct ConstantJumpEdge{T, U} <:JumpEdge
    rate::T # signature (vs, vd, p, t) -> Real
    affect!::U # signature (vsn, vdn, vs, vd, p, t) -> nothing 
end

"""
    VariableJumpEdge(rate, affect!)

Construct a variable rate jump over an edge with rate `rate` and `affect!` methods.

# Arguments
- `rate::Function`: signature `(vs, vd, p, t) -> Real`
- `affect!::Function`: signature `(vs, vd, nghbs, p, t) -> nothing`

See also: [Types of Jumps](https://docs.sciml.ai/JumpProcesses/stable/jump_types/#Types-of-Jumps:-Constant-Rate,-Mass-Action,-Variable-Rate-and-Regular)
"""
Base.@kwdef struct VariableJumpEdge{T, U} <:JumpEdge
    rate::T # signature (vs, vd, p, t) -> Real
    affect!::U # signature (vsn, vdn, vs, vd, nghbs, p, t) -> nothing
end

Base.@kwdef struct PreJumpSet
    constant::Vector{ConstantRateJump} = Vector{ConstantRateJump}()
    variable::Vector{VariableRateJump} = Vector{VariableRateJump}()
end

"""
    vertex_range(n, v)

Returns the range of vertices for a given vertex `v`, given that the number of states is `n`.
"""
function vertex_range(n, v)
    return n*(v-1)+1:n*v
end

function push_jump!(jumps::PreJumpSet, v, neighbors, vertex::ConstantJumpVertex; n=1)
    push!(jumps.constant, ConstantRateJump(
        (u, p, t) -> vertex.rate(u[vertex_range(n, v)], [view(u, vertex_range(n, nbhr)) for nbhr in neighbors], p, t),
        (integrator) -> vertex.affect!(
            view(integrator.u, vertex_range(n, v)),
            [view(integrator.u, vertex_range(n, nbhr)) for nbhr in neighbors],
            integrator.p, integrator.t
        ),
    ))
    nothing
end

function push_jump!(jumps::PreJumpSet, v, neighbors, vertex::VariableJumpVertex; n=1)
    push!(jumps.constant, VariableRateJump(
        (u, p, t) -> vertex.rate(u[vertex_range(n, v)], [view(u, vertex_range(n, nbhr)) for nbhr in neighbors], p, t),
        (integrator) -> vertex.affect!(
            view(integrator.u, vertex_range(n, v)),
            [view(integrator.u, vertex_range(n, nbhr)) for nbhr in neighbors],
            integrator.p, integrator.t
        ),
    ))
    nothing
end

function push_jump!(jumps::PreJumpSet, vs, vd, edge::ConstantJumpEdge; n=1)
    push!(jumps.constant, ConstantRateJump(
        (u, p, t) -> edge.rate(u[vertex_range(n, vs)], u[vertex_range(n, vd)], p, t),
        (integrator) -> edge.affect!(view(integrator.u, vertex_range(n, vs)), view(integrator.u, vertex_range(n, vd)), integrator.p, integrator.t),
    ))
    nothing
end

function push_jump!(jumps::PreJumpSet, vs, vd, edge::VariableJumpEdge; n=1)
    push!(jumps.variable, VariableRateJump(
        (u, p, t) -> edge.rate(u[vertex_range(n, vs)], u[vertex_range(n, vd)], p, t),
        (integrator) -> edge.affect!(view(integrator.u, vertex_range(n, vs)), view(integrator.u, vertex_range(n, vd)), integrator.p, integrator.t),
    ))
    nothing
end

"""
    network_jump_set(
        graph; vertex_reactions::Vector{T}=Vector{JumpVertex}(),
        edge_reactions::Vector{U}=Vector{JumpEdge}(), nb_states=1
    ) where {
        T <: Union{JumpVertex, Vector{<:JumpVertex}},
        U <: Union{JumpEdge, Vector{<:JumpEdge}}
    }

Construct a `JumpSet` from a `Graph` and a list of `JumpVertex` and `JumpEdge` reactions.
`vertex_reactions` and `edge_reactions` can be either an vector of reactions which
will all be applied to every vertex and edge respectively.
The other option is a vector of vectors of reactions, where the `i`th vector of reactions will be applied to the `i`th vertex or edge.
These variables may be mixed.
Each vertex has `nb_states` variables associated.

See also: [`ConstantJumpVertex`](@ref), [`ConstantJumpEdge`](@ref), [`VariableJumpVertex`](@ref), [`VariableJumpEdge`](@ref)
"""
function network_jump_set(graph;
                          vertex_reactions::Vector{T}=Vector{JumpVertex}(),
                          edge_reactions::Vector{U}=Vector{JumpEdge}(),
                          nb_states=1
                          ) where {
                            T <: Union{JumpVertex, Vector{<:JumpVertex}},
                            U <: Union{JumpEdge, Vector{<:JumpEdge}}
                          }
    jumps = PreJumpSet()
    max_nb = nv(graph)*length(vertex_reactions) + 2*ne(graph)*length(edge_reactions)
    sizehint!(jumps.constant, max_nb)
    sizehint!(jumps.variable, max_nb)

    if length(vertex_reactions) > 0
        if isa(vertex_reactions[1], Vector) && length(vertex_reactions) != nv(graph)
            throw(ArgumentError("vertex_reactions must be a vector of length nv(graph)"))
        end

        vertex_reactions_for = isa(vertex_reactions[1], Vector) ?
            v -> vertex_reactions[v] :
            v -> vertex_reactions

        for v in vertices(graph)
            nghbs = neighbors(graph, v)
            for jump_vertex in vertex_reactions_for(v)
                push_jump!(jumps, v, nghbs, jump_vertex; n=nb_states)
            end
        end
    end

    if length(edge_reactions) > 0
        if isa(edge_reactions[1], Vector) && length(edge_reactions) != ne(graph)
            throw(ArgumentError("edge_reactions must be a vector of length ne(graph)"))
        end

        edge_reactions_for = isa(edge_reactions[1], Vector) ?
            idx -> edge_reactions[idx] :
            idx -> edge_reactions

        for (idx, edge) in enumerate(edges(graph))
            for jump_edge in edge_reactions_for(idx)
                push_jump!(jumps, edge.src, edge.dst, jump_edge; n=nb_states)
                push_jump!(jumps, edge.dst, edge.src, jump_edge; n=nb_states)
            end
        end
    end

    hasConstantJumps = length(jumps.constant) > 0
    hasVariableJumps = length(jumps.variable) > 0

    if hasVariableJumps
        @warn "The use of variable rate jumps is still in an experimental phase."
    end

    if hasConstantJumps && hasVariableJumps
        return JumpSet(; constant_jumps=jumps.constant, variable_jumps=jumps.variable)
    elseif hasConstantJumps
        return JumpSet(; constant_jumps=jumps.constant)
    elseif hasVariableJumps
        return JumpSet(; variable_jumps=jumps.variable)
    else
        return JumpSet()
    end
end

function insert_vertex!(vec::AbstractArray, i)
    insert!(vec, searchsortedfirst(vec, i), i)
end

"""
    vertex_to_edges(graph::AbstractGraph)

Create a dictionary that maps each vertex to the edges it is connected to.
"""
function vertex_to_edges(graph::AbstractGraph)
    vte = Dict([v => [] for v in vertices(graph)]...)
    for (i, e) in enumerate(edges(graph))
        push!(vte[src(e)], i)
        push!(vte[dst(e)], i)
    end
    return vte
end

function dependency_map_input_preperation(graph::AbstractGraph, nb_vertex_reacs::T, nb_edge_reacs::U) where {
        T <: Union{Integer, Vector{<:Integer}},
        U <: Union{Integer, Vector{<:Integer}}
    }
    hetero_vertex = isa(nb_vertex_reacs, Vector)
    hetero_edge = isa(nb_edge_reacs, Vector)
    if hetero_vertex && length(nb_vertex_reacs) != nv(graph)
        throw(ArgumentError("nb_vertex_reacs must be a vector of length nv(graph) or an integer"))
    end
    if hetero_edge && length(nb_edge_reacs) != ne(graph)
        throw(ArgumentError("nb_edge_reacs must be a vector of length ne(graph) or an integer"))
    end

    # Number of vertex reactions
    nvr = i -> hetero_vertex ? nb_vertex_reacs[i] : nb_vertex_reacs
    # Number of edge reactions
    ner = i -> hetero_edge ? nb_edge_reacs[i] : nb_edge_reacs

    return nvr, ner
end


"""
    vartojumps(graph::AbstractGraph, nb_vertex_reacs::T, nb_edge_reacs::U, nb_vertex_states=1) where {
        T <: Union{Integer, Vector{<:Integer}},
        U <: Union{Integer, Vector{<:Integer}}
    }

Create a dependency graph that maps the vertex states to the vertex and edge reactions.
This graph can be used for the `RSSA` and `RSSACR` aggregators.
If all vertices (edges) have the same number of reactions, then `nb_vertex_reacs` (`nb_edge_reacs`)
can be an integer, equal to that number.
Otherwise, `nb_vertex_reacs` (`nb_edge_reacs`) must be a vector of length `nv(graph)` (`ne(graph)`)
, such that `nb_vertex_reacs[v]` (`nb_edge_reacs[e]`) is the number of reactions
associated with vertex `v` (edge `e`).

See also: [`jumptovars`](@ref), [`Jump Aggregators Requiring Dependency Graphs`](https://docs.sciml.ai/JumpProcesses/stable/jump_types/#Jump-Aggregators-Requiring-Dependency-Graphs) 
"""
function vartojumps(graph::AbstractGraph, nb_vertex_reacs::T, nb_edge_reacs::U, nb_vertex_states=1) where {
        T <: Union{Integer, Vector{<:Integer}},
        U <: Union{Integer, Vector{<:Integer}}
    }
    nvr, ner = dependency_map_input_preperation(graph, nb_vertex_reacs, nb_edge_reacs)
    nvs = nb_vertex_states
    vte = vertex_to_edges(graph)

    tnv = sum(v -> nvr(v), 1:nv(graph)) # total number of vertex reactions
    
    dep = Vector{Vector{Int64}}()
    for v in vertices(graph)
        nhbs = copy(neighbors(graph, v))
        insert_vertex!(nhbs, v)
        n_v_reactions = nvr(v)*length(nhbs)
        n_reactions = n_v_reactions + sum(i -> ner(i), vte[v])
        for _ in 1:nvs
            push!(dep, zeros(Int64, n_reactions))
            for i in eachindex(nhbs)
                dep[end][vertex_range(nvr(v), i)] = vertex_range(nvr(v), nhbs[i])
            end
            for i in eachindex(vte[v])
                dep[end][n_v_reactions .+ vertex_range(ner(i), i)] = tnv .+ vertex_range(ner(i), vte[v][i])
            end
        end
    end
    return dep
end

"""
    jumptovars(graph, nb_vertex_reacs, nb_edge_reacs, nb_vertex_states=1)

Create a dependency graph that maps the vertex and edge reactions to the vertex states.
This graph can be used for the `RSSA` and `RSSACR` aggregators.

See also: [`vartojumps`](@ref), [`Jump Aggregators Requiring Dependency Graphs`](https://docs.sciml.ai/JumpProcesses/stable/jump_types/#Jump-Aggregators-Requiring-Dependency-Graphs) 
"""
function jumptovars(graph, nb_vertex_reacs, nb_edge_reacs, nb_vertex_states=1)
    nvr, ner = dependency_map_input_preperation(graph, nb_vertex_reacs, nb_edge_reacs)
    nvs = nb_vertex_states
    
    dep = Vector{Vector{Int64}}()
    # Add dependences for the vertex reactions
    for v in vertices(graph)
        for _ in 1:nvr(v)
            push!(dep, vertex_range(nvs, v))
        end
    end
    for (e_idx, e) in enumerate(edges(graph)) 
        # If there is a reactiopn associated with the edge
        if ner(e_idx) > 0
            n_states = 2*nvs
            # Add the dependences for source to destination
            push!(dep, zeros(Int64, n_states))
            for (i, v) in enumerate((src(e), dst(e)))
                dep[end][vertex_range(nvs, i)] = vertex_range(nvs, v)
            end
            # Add in the other direction and for all subsequent reactions
            for _ in 1:2*ner(e_idx)-1
                push!(dep, dep[end])
            end
        end
    end

    return dep
end
end # module NetworkJumpProcesses
