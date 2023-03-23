using Test
using NetworkJumpProcesses
using JumpProcesses
using Graphs

@testset "Jump set (vertices)" begin
    # These dynamics should stay stable
    v_reaction1 = ConstantJumpVertex(
        (v, nghbs, p, t) -> 0.0,
        (vn, v, nghbs, p, t) -> vn .+= 1
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
    @test sum(sol[end]) == 4

    # These dynamics should average out
    v_reaction2 = ConstantJumpVertex(
        (v, nghbs, p, t) -> 20.0,
        (vn, v, nghbs, p, t) -> begin
            N_inv = 1/(length(nghbs) + 1)
            vn .= N_inv*(v[1] + sum(nghbs))
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
    @test maximum(sol[end]) ≈ minimum(sol[end])
end

@testset "Jump set (edges)" begin
    # These dynamics should stay stable
    e_reaction1 = ConstantJumpEdge(
        (vs, vd, p, t) -> 0.0,
        (vns, vnd, vs, vd, p, t) -> begin
            vns[1] = vs[1] + 1
            vnd[1] = vd[1] + 1
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
    @test sum(sol[end]) == 4

    # These dynamics should average out
    e_reaction2 = ConstantJumpEdge(
        (vs, vd, p, t) -> 20.0,
        (vns, vnd, vs, vd, p, t) -> begin
            avg = (vs[1] + vd[1])/2
            vns[1] = avg
            vnd[1] = avg
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
    @test maximum(sol[end]) ≈ minimum(sol[end])
end

@testset "Var to jump map" begin
    g = grid((2, 2))
    
    # Vertex unit test
    vtj = vartojumps(g, 1, 0)
    @test vtj == [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]]
    
    # Double state unit test
    vtj = vartojumps(g, 1, 0, 2)
    @test vtj == [
        [1, 2, 3], [1, 2, 3], [1, 2, 4], [1, 2, 4],
        [1, 3, 4], [1, 3, 4], [2, 3, 4], [2, 3, 4]
    ]

    # Double vertex reaction
    vtj = vartojumps(g, 2, 0)
    @test vtj == [
        [1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 7, 8],
        [1, 2, 5, 6, 7, 8], [3, 4, 5, 6, 7, 8]
    ]

    # Edge unit test
    vtj = vartojumps(g, 0, 1)
    @test vtj == [[1, 2], [1, 3], [2, 4], [3, 4]]

    # Double edge test
    vtj = vartojumps(g, 0, 2)
    @test vtj == [
        [1, 2, 3, 4], [1, 2, 5, 6],
        [3, 4, 7, 8], [5, 6, 7, 8]
    ]

    # Double vertex and edge reaction
    vtj = vartojumps(g, 2, 1)
    @test vtj == [
        [1, 2, 3, 4, 5, 6, 9, 10],  [1, 2, 3, 4, 7, 8, 9, 11],
        [1, 2, 5, 6, 7, 8, 10, 12], [3, 4, 5, 6, 7, 8, 11, 12]
    ]
end

@testset "Jump to var map" begin
    g = grid((2, 2))

    # Vertex unit test
    jtv = jumptovars(g, 1, 0)
    @test jtv == [[1], [2], [3], [4]]

    # Double state
    jtv = jumptovars(g, 1, 0, 2)
    @test jtv == [[1, 2], [3, 4], [5, 6], [7, 8]]

    # Double vertex test
    jtv = jumptovars(g, 2, 0)
    @test jtv == [[1], [1], [2], [2], [3], [3], [4], [4]]

    # Edge unit test
    jtv = jumptovars(g, 0, 1)
    @test jtv == [
        [1, 2], [2, 1], [1, 3], [3, 1],
        [2, 4], [4, 2], [3, 4], [4, 3]
    ]

    # Double edge test
    jtv = jumptovars(g, 0, 2)
    @test jtv == [
        [1, 2], [2, 1], [1, 2], [2, 1],
        [1, 3], [3, 1], [1, 3], [3, 1],
        [2, 4], [4, 2], [2, 4], [4, 2],
        [3, 4], [4, 3], [3, 4], [4, 3]
    ]

    # Double vertex and edge reaction
    jtv = jumptovars(g, 2, 1)
    @test jtv == [
        [1], [1], [2], [2], [3], [3], [4], [4],
        [1, 2], [2, 1], [1, 3], [3, 1],
        [2, 4], [4, 2], [3, 4], [4, 3]
    ]
end