using Graphs
using JumpProcesses
using Plots: plot
using Random
using Revise
using NetworkJumpProcesses

# In this expample we define source nodes for two chemical elements, A and B, which
# might react with each other and form C. C is disposed of in some sink nodes. Over
# the edges there happens diffusion.

n = 12
gg = erdos_renyi(n, 0.4)

A_source = ConstantJumpVertex(
    (v, nghbs, p, t) -> p.rA,
    (v, nghbs, p, t) -> v[1] += 1
)

B_source = ConstantJumpVertex(
    (v, nghbs, p, t) -> p.rB,
    (v, nghbs, p, t) -> v[2] += 1
)

C_creation = ConstantJumpVertex(
    (v, nghbs, p, t) -> p.rC*v[1]*v[2],
    (v, nghbs, p, t) -> begin
        v[1] -= 1
        v[2] -= 1
        v[3] += 1
    end
)

C_sink = ConstantJumpVertex(
    (v, nghbs, p, t) -> p.dC*v[3],
    (v, nghbs, p, t) -> v[3] -= 1
)

diffusion = i -> ConstantJumpEdge(
    (vs, vd, p, t) -> max(p.D*(vs[i] - vd[i]), 0.),
    (vs, vd, p, t) -> begin
        vs[i] -= 1
        vd[i] += 1
    end
)

# Define the vertex reactions
A_source_reactions = [[A_source, C_creation] for i in 1:3]
B_source_reactions = [[B_source, C_creation] for i in 1:3]
C_sink_reactions = [[C_sink, C_creation] for i in 1:3]
v_reactions = vcat(A_source_reactions, B_source_reactions, C_sink_reactions)

# Add the reactions for the nodes that are not a source or sink
append!(v_reactions, [[C_creation] for i in 1:(n-length(v_reactions))])
shuffle!(v_reactions)

# Define the edge reactions
e_reactions = [diffusion(i) for i in 1:3]
jset = network_jump_set(gg; vertex_reactions=v_reactions, edge_reactions=e_reactions, nb_states=3)

p = (
    rA = 0.1,
    rB = 0.1,
    rC = 0.2,
    dC = 1.5,
    D = 0.1
)

T = 100.0
nb_T = Int(T)
u₀ = zeros(Int64, 3*n)
dprob = DiscreteProblem(u₀, (0., T), p) 
direct_jump_prob = JumpProblem(dprob, Direct(), jset)
sol = solve(direct_jump_prob, SSAStepper())

# Let's plot the evolution of concentrations A, B and C
get_concentrations(x::Vector{Int64}) = [
    sum(x[i] for i in 1:3:length(x)),
    sum(x[i] for i in 2:3:length(x)),
    sum(x[i] for i in 3:3:length(x))
]

get_concentrations(x::ODESolution) = reduce(hcat, [get_concentrations(x[t]) for t in eachindex(x)])';

plot(sol.t, get_concentrations(sol), label=["A" "B" "C"], xlabel="t", ylabel="concentration", lw=2, legend=:topleft)

# And calculate a mean solution
function calculate_mean(prob, T, N)
    mean_SSAS = zeros(nb_T+1, 3)
    for i in 1:N
        sol = solve(prob, SSAStepper())
        for t in 0:nb_T
            mean_SSAS[t+1,:] += get_concentrations(sol(t))./N
        end 
    end
    return mean_SSAS
end

mean_direct = calculate_mean(direct_jump_prob, nb_T, 50) # This might take a while

# And plot it
plot(0:nb_T, mean_direct, label=["A" "B" "C"], xlabel="t", ylabel="concentration", lw=2, legend=:topleft, title="Direct method")

# Now do the same, but with an aggregator that makes use of dependency graphs
nb_vertex_reacs = [length(v) for v in v_reactions]
vtj = vartojumps(gg, nb_vertex_reacs, 3, 3)
jtv = jumptovars(gg, nb_vertex_reacs, 3, 3)
rssa_jump_prob = JumpProblem(dprob, RSSA(), jset; vartojumps_map=vtj, jumptovars_map=jtv)
sol = solve(rssa_jump_prob, SSAStepper())

plot(sol.t, get_concentrations(sol), label=["A" "B" "C"], xlabel="t", ylabel="concentration", lw=2, legend=:topleft)

# And calculate a mean solution
mean_RSSA = calculate_mean(rssa_jump_prob, nb_T, 250)
plot(0:nb_T, mean_RSSA, label=["A" "B" "C"], xlabel="t", ylabel="concentration", lw=2, legend=:topleft, title="RSSA")
