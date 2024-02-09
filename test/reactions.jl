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
