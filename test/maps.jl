@testitem "Var to jump map" begin
    using Graphs

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
    @test vtj == [
        [1, 2, 3, 4], [1, 2, 5, 6],
        [3, 4, 7, 8], [5, 6, 7, 8]
    ]

    # Double edge test
    vtj = vartojumps(g, 0, 2)
    @test vtj == [
        [1, 2, 3, 4, 5, 6, 7, 8],
        [1, 2, 3, 4, 9, 10, 11, 12],
        [5, 6, 7, 8, 13, 14, 15, 16],
        [9, 10, 11, 12, 13, 14, 15, 16]
    ]

    # Double vertex and edge reaction
    vtj = vartojumps(g, 2, 1)
    @test vtj == [
        [1, 2, 3, 4, 5, 6, 9, 10, 11, 12],
        [1, 2, 3, 4, 7, 8, 9, 10, 13, 14],
        [1, 2, 5, 6, 7, 8, 11, 12, 15, 16],
        [3, 4, 5, 6, 7, 8, 13, 14, 15, 16]
    ]

    # Multiple of everything
    vtj = vartojumps(g, 2, 3, 3)

    for i in 1:3
        @test vtj[i] == [ # Node 1
            1, 2, 3, 4, 5, 6,
            9, 10, 11, 12, 13, 14,
            15, 16, 17, 18, 19, 20
        ]
        @test vtj[i+3] == [ # Node 2
            1, 2, 3, 4, 7, 8,
            9, 10, 11, 12, 13, 14,
            21, 22, 23, 24, 25, 26
        ]
        @test vtj[i+6] == [ # Node 3
            1, 2, 5, 6, 7, 8,
            15, 16, 17, 18, 19, 20,
            27, 28, 29, 30, 31, 32
        ]
        @test vtj[i+9] == [ # Node 4
            3, 4, 5, 6, 7, 8,
            21, 22, 23, 24, 25, 26,
            27, 28, 29, 30, 31, 32
        ]
    end

    # Heterogeneous input
    @test_throws ArgumentError vartojumps(g, [1, 2], 0)
    @test_throws ArgumentError vartojumps(g, 1, [0, 1])

    vtj = vartojumps(g, [1, 2, 1, 2], 0)
    @test vtj == [
        [1, 2, 3, 4], [1, 2, 3, 5, 6],
        [1, 4, 5, 6], [2, 3, 4, 5, 6],
    ]
    vtj = vartojumps(g, [0, 2, 0, 0], 1)
    @test vtj == [
        [1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 7, 8],
        [5, 6, 9, 10], [1, 2, 7, 8, 9, 10],
    ]

    vtj = vartojumps(g, 1, [0, 1, 2, 3])
    @test vtj == [
        [1, 2, 3, 5, 6], [1, 2, 4, 7, 8, 9, 10],
        [1, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16],
        [2, 3, 4, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
    ]
    vtj = vartojumps(g, [0, 1, 2, 3], [0, 1, 1, 0], 2)
    for i in 1:2
        @test vtj[i] == [1, 2, 3, 7, 8] # Node 1
        @test vtj[i+2] == [1, 4, 5, 6, 9, 10] # Node 2
        @test vtj[i+4] == [2, 3, 4, 5, 6, 7, 8] # Node 3
        @test vtj[i+6] == [1, 2, 3, 4, 5, 6, 9, 10] # Node 4
    end
end

@testitem "Jump to var map" begin
    using Graphs

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
        [1, 2], [1, 2], [1, 3], [1, 3],
        [2, 4], [2, 4], [3, 4], [3, 4]
    ]

    # Double edge test
    jtv = jumptovars(g, 0, 2)
    @test jtv == [
        [1, 2], [1, 2], [1, 2], [1, 2],
        [1, 3], [1, 3], [1, 3], [1, 3],
        [2, 4], [2, 4], [2, 4], [2, 4],
        [3, 4], [3, 4], [3, 4], [3, 4]
    ]

    # Double vertex and edge reaction
    jtv = jumptovars(g, 2, 1)
    @test jtv == [
        [1], [1], [2], [2], [3], [3], [4], [4],
        [1, 2], [1, 2], [1, 3], [1, 3],
        [2, 4], [2, 4], [3, 4], [3, 4]
    ]

    # Multiple of everything
    jtv = jumptovars(g, 2, 3, 3)

    # Check the reactions on the vertices themselves
    # jtv[1:2] .== [1, 2, 3]
    # jtv[3:4] .== [4, 5, 6] ...
    @test length(jtv) == 32
    for i in 1:4
        @test jtv[2*(i-1)+1] == 3*(i-1) .+ [1, 2, 3]
        @test jtv[2*(i-1)+2] == 3*(i-1) .+ [1, 2, 3]
    end

    # Check the reactions on the edges
    for j in 1:6
        @test jtv[8+j] == [1, 2, 3, 4, 5, 6]
    end

    for j in 1:6
        @test jtv[14+j] == [1, 2, 3, 7, 8, 9]
    end

    for j in 1:6
        @test jtv[20+j] == [4, 5, 6, 10, 11, 12]
    end

    for j in 1:6
        @test jtv[26+j] == [7, 8, 9, 10, 11, 12]
    end
end