# NetworkJumpProcesses.jl
[![Build Status](https://github.com/TimVWese/NetworkJumpProcesses.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/TimVWese/NetworkJumpProcesses.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://timvwese.github.io/NetworkJumpProcesses.jl/dev/)

Julia package to facilitate the construction of JumpProblems on graphs.
The idea is that the relationship between this package and the [JumpProcesses.jl](https://github.com/SciML/JumpProcesses.jl) package is similar to the relationship between [NetworkDynamics.jl](https://github.com/PIK-ICoNe/NetworkDynamics.jl) and [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).

## Get started
Install with

```julia
using Pkg
Pkg.add("NetworkJumpProcesses")
```

See the [example folder](https://github.com/TimVWese/NetworkJumpProcesses.jl/tree/main/examples) and [documentation](https://timvwese.github.io/NetworkJumpProcesses.jl/dev/) to get started.

## Disclaimer
The package is in an early stage.
Reporting of encountered bugs or desired functionality is appreciated, development aid even more so.
