@testitem "Jump set (vertices)" begin
    using Graphs, JumpProcesses

    # These dynamics should stay stable
    v_reaction1 = ConstantJumpVertex(
        (v, nghbs, p, t) -> 0.0,
        (v, nghbs, p, t) -> v[1] += 1
    )

    g = grid((2,2))
    jset = network_jump_set(g; vertex_reactions=[v_reaction1])

    @test length(jset.constant_jumps) == 4
    @test length(jset.variable_jumps) == 0
    @test isnothing(jset.massaction_jump)
    @test isnothing(jset.regular_jump)

    u0 = ones(Int64, 4)
    dprob = DiscreteProblem(u0, (0., 10.))
    jprob = JumpProblem(dprob, Direct(), jset)
    sol = solve(jprob, SSAStepper())
    @test sum(sol[:, end]) == 4

    # These dynamics should average out
    v_reaction2 = ConstantJumpVertex(
        (v, nghbs, p, t) -> 20.0,
        (v, nghbs, p, t) -> begin
            N_inv = 1/(length(nghbs) + 1)
            v[1] = N_inv*(v[1] + sum(x -> x[1], nghbs))
        end
    )

    g = complete_graph(6)
    jset = network_jump_set(g; vertex_reactions=[v_reaction2])
    @test length(jset.constant_jumps) == 6
    @test length(jset.variable_jumps) == 0
    @test isnothing(jset.massaction_jump)
    @test isnothing(jset.regular_jump)

    u0 = Float64.(0:5)
    dprob = DiscreteProblem(u0, (0., 200.))
    jprob = JumpProblem(dprob, Direct(), jset)
    sol = solve(jprob, SSAStepper())
    @test maximum(sol[:, end]) ≈ minimum(sol[:, end])
end

@testitem "Jump set (edges)" begin
    using Graphs, JumpProcesses

    # These dynamics should stay stable
    e_reaction1 = ConstantJumpEdge(
        (vs, vd, p, t) -> 0.0,
        (vs, vd, p, t) -> begin
            vs[1] += 1
            vd[1] += 1
        end
    )

    g = grid((2,2))
    jset = network_jump_set(g; edge_reactions=[e_reaction1])

    @test length(jset.constant_jumps) == 8
    @test length(jset.variable_jumps) == 0
    @test isnothing(jset.massaction_jump)
    @test isnothing(jset.regular_jump)

    u0 = ones(Int64, 4)
    dprob = DiscreteProblem(u0, (0., 10.))
    jprob = JumpProblem(dprob, Direct(), jset)
    sol = solve(jprob, SSAStepper())
    @test sum(sol[:, end]) == 4

    # These dynamics should average out
    e_reaction2 = ConstantJumpEdge(
        (vs, vd, p, t) -> 20.0,
        (vs, vd, p, t) -> begin
            avg = (vs[1] + vd[1])/2
            vs[1] = avg
            vd[1] = avg
        end
    )

    g = complete_graph(6)
    jset = network_jump_set(g; edge_reactions=[e_reaction2])
    @test length(jset.constant_jumps) == 30
    @test length(jset.variable_jumps) == 0
    @test isnothing(jset.massaction_jump)
    @test isnothing(jset.regular_jump)

    u0 = Float64.(0:5)
    dprob = DiscreteProblem(u0, (0., 200.))
    jprob = JumpProblem(dprob, Direct(), jset)
    sol = solve(jprob, SSAStepper())
    @test maximum(sol[:, end]) ≈ minimum(sol[:, end])
end

@testitem "Multiple state reactions" begin
    using Graphs, JumpProcesses

    # These dynamics should average out
    v_reaction = ConstantJumpVertex(
        (v, nghbs, p, t) -> 20.0,
        (v, nghbs, p, t) -> begin
            N_inv = 1/(length(nghbs) + 1)
            v[1] = N_inv*(v[1] + sum(x -> x[1], nghbs))
        end
    )

    e_reaction = ConstantJumpEdge(
        (vs, vd, p, t) -> 20.0,
        (vs, vd, p, t) -> begin
            avg = (vs[2] + vd[2])/2
            vs[2] = avg
            vd[2] = avg
        end
    )

    g = complete_graph(6)
    jset = network_jump_set(g; vertex_reactions=[v_reaction], edge_reactions=[e_reaction], nb_states=2)
    @test length(jset.constant_jumps) == 36
    @test length(jset.variable_jumps) == 0
    @test isnothing(jset.massaction_jump)
    @test isnothing(jset.regular_jump)

    u0 = Float64.(0:11)
    dprob = DiscreteProblem(u0, (0., 200.))
    jprob = JumpProblem(dprob, Direct(), jset)
    sol = solve(jprob, SSAStepper())
    @test maximum(sol[1:2:end, end]) ≈ minimum(sol[1:2:end, end])
    @test maximum(sol[2:2:end, end]) ≈ minimum(sol[2:2:end, end])
end

@testitem "Variable jumps" begin
    using Graphs, JumpProcesses, DifferentialEquations, Logging

    g = grid((2, 2))

    v_var = VariableJumpVertex(
        (v, nghbs, p, t) -> max(v[1], 0.0),
        (v, nghbs, p, t) -> v[1] = max(v[1] - 1.0, 0.0)
    )
    e_var = VariableJumpEdge(
        (vs, vd, p, t) -> 1.0,
        (vs, vd, p, t) -> vd[1] += 1.0
    )

    # Variable vertex jumps must end up as variable_jumps (not constant_jumps).
    # Constructing them also warns that variable rate jumps are experimental.
    jset = @test_logs (:warn,) match_mode = :any network_jump_set(g; vertex_reactions=[v_var])
    @test length(jset.constant_jumps) == 0
    @test length(jset.variable_jumps) == 4
    @test isnothing(jset.massaction_jump)
    @test isnothing(jset.regular_jump)

    # Variable edge jumps, one per direction.
    jset = with_logger(NullLogger()) do
        network_jump_set(g; edge_reactions=[e_var])
    end
    @test length(jset.constant_jumps) == 0
    @test length(jset.variable_jumps) == 8

    # Mixing constant and variable jumps populates both vectors.
    v_const = ConstantJumpVertex(
        (v, nghbs, p, t) -> 1.0,
        (v, nghbs, p, t) -> v[1] += 1.0
    )
    jset = with_logger(NullLogger()) do
        network_jump_set(g; vertex_reactions=[v_const, v_var])
    end
    @test length(jset.constant_jumps) == 4
    @test length(jset.variable_jumps) == 4

    # End-to-end: a pure-decay variable rate jump should reduce every vertex
    # towards zero without overshooting below it.
    jset = with_logger(NullLogger()) do
        network_jump_set(g; vertex_reactions=[v_var])
    end
    u0 = fill(100.0, nv(g))
    oprob = ODEProblem((du, u, p, t) -> (du .= 0.0), u0, (0.0, 10.0))
    jprob = with_logger(NullLogger()) do
        JumpProblem(oprob, Direct(), jset)
    end
    sol = solve(jprob, Tsit5())
    phys = sol[1:nv(g), end]
    @test all(phys .< 100.0)
    @test all(phys .>= 0.0)
end

@testitem "Heterogeneous and empty reactions" begin
    using Graphs, JumpProcesses

    g = grid((2, 2)) # 4 vertices, 4 edges

    v_reaction = ConstantJumpVertex(
        (v, nghbs, p, t) -> 1.0,
        (v, nghbs, p, t) -> v[1] += 1
    )
    e_reaction = ConstantJumpEdge(
        (vs, vd, p, t) -> 1.0,
        (vs, vd, p, t) -> vs[1] += 1
    )

    # Per-vertex reaction counts 1, 2, 3, 4 -> 10 constant jumps
    v_reactions = [[v_reaction for _ in 1:i] for i in 1:nv(g)]
    jset = network_jump_set(g; vertex_reactions=v_reactions)
    @test length(jset.constant_jumps) == 10
    @test length(jset.variable_jumps) == 0

    # Per-edge reaction counts 1, 2, 3, 4, doubled for both directions -> 20
    e_reactions = [[e_reaction for _ in 1:i] for i in 1:ne(g)]
    jset = network_jump_set(g; edge_reactions=e_reactions)
    @test length(jset.constant_jumps) == 20

    # A vector-of-vectors must match the number of vertices / edges
    @test_throws ArgumentError network_jump_set(g; vertex_reactions=[[v_reaction], [v_reaction]])
    @test_throws ArgumentError network_jump_set(g; edge_reactions=[[e_reaction], [e_reaction]])

    # No reactions yields an empty JumpSet
    jset = network_jump_set(g)
    @test length(jset.constant_jumps) == 0
    @test length(jset.variable_jumps) == 0
    @test isnothing(jset.massaction_jump)
    @test isnothing(jset.regular_jump)
end
