# NetworkJumpProcesses.jl

An interface between [Graphs.jl](https://juliagraphs.org/) and [JumpProcesses.jl](https://github.com/SciML/JumpProcesses.jl).

## Package features

- Given jumps defined on edges and vertices, the package constructs a jump process on the (undirected) graph.
- Construct a `vartojumps_map` and `jumps_tovars_map` to use as dependency graph for a the `RSSA(CR)` aggregator.

See the [Reference](@ref) for all available functions.

There is more information about installation and use of the package in the [Package Guide](@ref) and hands-on examples in the [Examples](@ref) section.

## Manual outline

```@contents
Pages = [
    "guide.md",
    "examples.md",
    "reference.md",
]
Depth = 1
```

## Package Limitations

As this is my first package and it has been developed for a specific use case, it is not as general as it could be.

- Directed graphs are not supported.
- jump-to-jump dependencies are not supported.
- `VariableJumps`` are not tested.
- Coding style and performance-wise improvements.
- ... (please let me know!)
