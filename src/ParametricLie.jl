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

# Core building blocks for v0.1 — implement these next.
include("types.jl")
include("bracket.jl")
include("invariants.jl")
include("analyze.jl")

export LieAlgebra, LieAlgebraElem, dim, coefficient_ring, base_ring, structure_constants
export coords, basis_elem
export lie_bracket, lie_bracket!, bracket, ad, check_jacobi, check_antisymmetry, JacobiCertificate
export center, derived_series, lower_central_series
export is_solvable, is_nilpotent, killing_form, derivations, change_of_basis
export analyze

end # module ParametricLie
