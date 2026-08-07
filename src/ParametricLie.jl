"""
    ParametricLie

Exact symbolic toolkit for finite-dimensional Lie algebras given by structure
constants (parametric families, stratification, cohomology — see
`docs/POSITIONING.md`). Complements LieGroups.jl; does not replace it.

v0.1: coefficient field via AbstractAlgebra/Nemo, structure, Jacobi, invariants, `analyze`.
"""
module ParametricLie

using LinearAlgebra
using SparseArrays
using AbstractAlgebra
using Nemo

# Include order = dependency order. Future: radical.jl, quotient.jl, levi.jl.
include("types.jl")
include("bracket.jl")
include("change_of_basis.jl")
include("subspace.jl")
include("center.jl")
include("series.jl")
include("killing.jl")
include("derivations.jl")
include("analyze.jl")

export LieAlgebra, LieAlgebraElem, dim, coefficient_ring, base_ring, structure_constants
export coords, basis_elem
export lie_bracket, lie_bracket!, bracket, ad, check_jacobi, check_antisymmetry, JacobiCertificate
export center, derived_series, lower_central_series, derived_algebra
export is_solvable, is_nilpotent, derivations, change_of_basis
export killing_form, killing_rank, killing_radical, killing_orthogonal
export cartan_orthogonal, is_cartan_solvable, is_semisimple
export LieSubspace, basis_matrix, basis_elems, full_space, zero_space
export commutator_span, LieSeries, SeriesKind, terms, layers
export analyze

end # module ParametricLie
