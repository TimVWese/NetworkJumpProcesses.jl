using NetworkJumpProcesses
using Graphs
using JumpProcesses
using Plots: plot

# 0 is susceptible state
# 1 is infected state

# Healing event (intra vertex jump)
v_IS = ConstantJumpVertex(
    (v, nghbs, p, t) -> v == 1 ? p[2] : 0.0,
    (vn, v, nghbs, p, t) -> vn .= 0
)

# Infection event (jump over an edge)
e_SI = ConstantJumpEdge(
    (vs, vd, p, t) -> vs == 1 ? p[1] : 0.0,
    (vns, vnd, vs, vd, p, t) -> vnd .= 1
)

# Define the (lattice) network and JumpSet
n = 10
gg = grid((n, n))
jset = network_jump_set(gg; vertex_reactions=[v_IS], edge_reactions=[e_SI])

# Define the problem with parameters (infection rate, helaing rate)
# and intial conditions (one infected agent)

p = (0.1, 0.08)
u₀ = zeros(Int64, n^2)
u₀[1] = 1

# Define the problem and simulate it
dprob = DiscreteProblem(u₀, (0, 40.0), p) 
jprob = JumpProblem(dprob, Direct(), jset)
sol = solve(jprob, SSAStepper())

# Plot the results
plot(sol.t, [count(sol[i].==1) for i in eachindex(sol.t)])

