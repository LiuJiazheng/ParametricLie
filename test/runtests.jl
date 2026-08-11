using Test
using ParametricLie
# Prefer `import` so Nemo's exports (dim, coefficient_ring, …) do not clash.
import AbstractAlgebra
import Nemo

@testset "ParametricLie.jl" begin
    # --- lie (single-algebra core) -------------------------------------------
    include("unit/lie/types.jl")
    include("unit/lie/bracket.jl")
    include("unit/lie/change_of_basis.jl")
    include("unit/lie/center.jl")
    include("unit/lie/series.jl")
    include("unit/lie/killing.jl")
    include("unit/lie/levi.jl")
    include("unit/lie/ideal_decomp.jl")
    include("unit/lie/derivations.jl")
    include("unit/lie/analyze.jl")

    # --- parametric (families & stratification) ------------------------------
    include("unit/parametric/parametric.jl")
    include("unit/parametric/specialize.jl")
    include("unit/parametric/conditional.jl")
    include("unit/parametric/stratify.jl")

    # --- cohomology (exterior algebra, CE complex) ---------------------------
    include("unit/cohomology/exterior.jl")
    include("unit/cohomology/ce.jl")
    include("unit/cohomology/central_extension.jl")
end
