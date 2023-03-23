module NetworkJumpProcesses

using Graphs
using JumpProcesses

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
- `affect!::Function`: signature `(vn, v, nghbs, p, t)`

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
- `affect!::Function`: signature `(vn, v, nghbs, p, t) -> nothing`

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
- `rate::Function`: signature `(v, nhgbs, p, t) -> Real`
- `affect!::Function`: signature `(vn, v, nghbs, p, t) -> nothing`

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

function push_jump!(jumps::PreJumpSet, v, neighbors, vertex::ConstantJumpVertex)
    push!(jumps.constant, ConstantRateJump(
        (u, p, t) -> vertex.rate(u[v], view(u, neighbors), p, t),
        (integrator) -> vertex.affect!(view(integrator.u, v), integrator.u[v], view(integrator.u, neighbors), integrator.p, integrator.t),
    ))
    nothing
end

function push_jump!(jumps::PreJumpSet, v, neighbors, vertex::VariableJumpVertex)
    push!(jumps.variable, VariableRateJump(
        (u, p, t) -> vertex.rate(u[v], view(u, neighbors), p, t),
        (integrator) -> vertex.affect!(view(integrator.u, v), integrator.u[v], view(integrator.u, neighbors), integrator.p, integrator.t),
    ))
    nothing
end

function push_jump!(jumps::PreJumpSet, vs, vd, edge::ConstantJumpEdge)
    push!(jumps.constant, ConstantRateJump(
        (u, p, t) -> edge.rate(u[vs], u[vd], p, t),
        (integrator) -> edge.affect!(view(integrator.u, vs), view(integrator.u, vd), integrator.u[vs], integrator.u[vd], integrator.p, integrator.t),
    ))
    nothing
end

function push_jump!(jumps::PreJumpSet, vs, vd, edge::VariableJumpEdge)
    push!(jumps.variable, VariableRateJump(
        (u, p, t) -> edge.rate(u[vs], u[vd], p, t),
        (integrator) -> edge.affect!(view(integrator.u, vs), view(integrator.u, vd), integrator.u[vs], integrator.u[vd], integrator.p, integrator.t),
    ))
    nothing
end

"""
    network_jump_set(graph; vertex_reactions::Vector{T}=Vector{JumpVertex}(), edge_reactions::Vector{U}=Vector{JumpEdge}()) where {T <: JumpVertex, U <: JumpEdge}

Construct a `JumpSet` from a `Graph` and a list of `JumpVertex` and `JumpEdge` reactions.

See also: [`ConstantJumpVertex`](@ref), [`ConstantJumpEdge`](@ref), [`VariableJumpVertex`](@ref), [`VariableJumpEdge`](@ref)
"""
function network_jump_set(graph;
                          vertex_reactions::Vector{T}=Vector{JumpVertex}(),
                          edge_reactions::Vector{U}=Vector{JumpEdge}()) where {T <: JumpVertex, U <: JumpEdge}

    jumps = PreJumpSet()
    max_nb = nv(graph)*length(vertex_reactions) + 2*ne(graph)*length(edge_reactions)
    sizehint!(jumps.constant, max_nb)

    if length(vertex_reactions) > 0
        for v in vertices(graph)
            nghbs = neighbors(graph, v)
            for jump_vertex in vertex_reactions
                push_jump!(jumps, v, nghbs, jump_vertex)
            end
        end
    end

    if length(edge_reactions) > 0
        for edge in edges(graph)
            for jump_edge in edge_reactions
                push_jump!(jumps, edge.src, edge.dst, jump_edge)
                push_jump!(jumps, edge.dst, edge.src, jump_edge)
            end
        end
    end

    hasConstantJumps = length(jumps.constant) > 0
    hasVariableJumps = length(jumps.variable) > 0

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

vertex_range = (n, v) -> n*(v-1)+1:n*v
insert_vertex!(vec::AbstractArray, i) = insert!(vec, searchsortedfirst(vec, i), i)

function vertex_to_edges(graph::AbstractGraph)
    vte = Dict([v => [] for v in vertices(graph)]...)
    for (i, e) in enumerate(edges(graph))
        push!(vte[src(e)], i)
        push!(vte[dst(e)], i)
    end
    return vte
end

"""
    vartojumps(graph, nb_vertex_states, nb_vertex_reacs, nb_edge_reacs)

Create a dependency graph that maps the vertex states to the vertex and edge reactions.
This graph can be used for the `RSSA` and `RSSACR` aggregators.

See also: [`jumptovars`](@ref), [`Jump Aggregators Requiring Dependency Graphs`](https://docs.sciml.ai/JumpProcesses/stable/jump_types/#Jump-Aggregators-Requiring-Dependency-Graphs) 
"""
function vartojumps(graph, nb_vertex_states, nb_vertex_reacs, nb_edge_reacs)
    nvs = number_vertex_states
    nvr = nb_vertex_reacs
    ner = nb_edge_reacs
    vte = vertex_to_edges(graph)

    tnv = nv(graph)*nvr
    
    dep = Vector{Vector{Int64}}()
    for v in vertices(graph)
        nhbs = copy(neighbors(graph, v))
        insert_vertex!(nhbs, v)
        n_v_reactions = nvr*length(nhbs)
        n_reactions = n_v_reactions + ner*(length(nhbs)-1)
        for _ in 1:nvs
            push!(dep, zeros(Int64, n_reactions))
            for i in eachindex(nhbs)
                dep[end][vertex_range(nvr, i)] = vertex_range(nvr, nhbs[i])
            end
            for i in eachindex(vte[v])
                dep[end][vertex_range(ner, n_v_reactions + i)] = vertex_range(ner, tnv + vte[v][i])
            end
        end
    end
    return dep
end

"""
    jumptovars(graph, nb_vertex_states, nb_vertex_reacs, nb_edge_reacs)

Create a dependency graph that maps the vertex and edge reactions to the vertex states.
This graph can be used for the `RSSA` and `RSSACR` aggregators.

See also: [`vartojumps`](@ref), [`Jump Aggregators Requiring Dependency Graphs`](https://docs.sciml.ai/JumpProcesses/stable/jump_types/#Jump-Aggregators-Requiring-Dependency-Graphs) 
"""
function jumptovars(graph, nb_vertex_states, nb_vertex_reacs, nb_edge_reacs)
    nvr = nb_vertex_reacs
    ner = nb_edge_reacs
    nvs = nb_vertex_states
    
    dep = Vector{Vector{Int64}}()
    if nvr > 0
        for v in vertices(graph)
            nhbs = copy(neighbors(graph, v))
            insert_vertex!(nhbs, v)
            n_states = nvs*length(nhbs)

            push!(dep, zeros(Int64, n_states))
            for i in eachindex(nhbs)
                dep[end][vertex_range(nvs, i)] = vertex_range(nvs, nhbs[i])
            end
            for i in 1:nvr-1
                push!(dep, dep[end])
            end

        end
    end
    if ner > 0
        for e in edges(graph)
            n_states = 2*nvs

            push!(dep, zeros(Int64, n_states))
            for (i, v) in enumerate((src(e), dst(e)))
                dep[end][vertex_range(nvs, i)] = vertex_range(nvs, v)
            end
            for i in 1:ner-1
                push!(dep, dep[end])
            end

        end
    end
    return dep
end
end # module NetworkJumpProcesses