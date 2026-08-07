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

# Include order = dependency order.
include("types.jl")
include("bracket.jl")
include("change_of_basis.jl")
include("subspace.jl")
include("center.jl")
include("series.jl")
include("killing.jl")
include("radical.jl")
include("quotient.jl")
include("levi.jl")
include("ideal_decomp.jl")
include("derivations.jl")
include("analyze.jl")

export LieAlgebra, LieAlgebraElem, dim, coefficient_ring, base_ring, structure_constants
export coords, basis_elem
export lie_bracket, lie_bracket!, bracket, ad, check_jacobi, check_antisymmetry, JacobiCertificate
export center, derived_series, lower_central_series, derived_algebra
export is_solvable, is_nilpotent, derivations, Derivations, basis_matrices
export is_derivation, apply_derivation, change_of_basis
export killing_form, killing_rank, killing_radical, killing_orthogonal
export cartan_orthogonal, is_cartan_solvable, is_semisimple
export LieSubspace, basis_matrix, basis_elems, full_space, zero_space
export commutator_span, complement, is_subalgebra, is_ideal
export LieSeries, SeriesKind, terms, layers
export radical, radical_derived_series
export QuotientAlgebra, quotient_algebra
export LeviDecomposition, levi_decomposition, levi_subalgebra
export adjoint_commutant, ideal_decomposition, is_simple
export analyze, LieAlgebraReport, jacobi, levi_kind, simple_factor_dims

end # module ParametricLie
