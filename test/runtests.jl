using Test
using ParametricLie
# Prefer `import` so Nemo's exports (dim, coefficient_ring, …) do not clash.
import AbstractAlgebra
import Nemo

@testset "ParametricLie.jl" begin
    # Layer 1 — unit tests for small components
    include("unit/types.jl")
    include("unit/bracket.jl")
    include("unit/change_of_basis.jl")
    include("unit/center.jl")

    # Layer 2 — mathematical identities (Jacobi, …)
    # include("identities/jacobi.jl")

    # Layer 3 — known examples / fixtures
    include("examples/smoke.jl")
end
