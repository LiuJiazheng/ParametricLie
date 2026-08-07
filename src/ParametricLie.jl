"""
    ParametricLie

Exact symbolic toolkit for finite-dimensional Lie algebras given by structure
constants (parametric families, stratification, cohomology — see
`docs/POSITIONING.md`). Complements LieGroups.jl; does not replace it.

Layout:

- `src/lie/` — single-algebra core (structure, invariants, `analyze`)
- `src/parametric/` — parameterized families (`specialize` / generic analyze)
"""
module ParametricLie

using LinearAlgebra
using SparseArrays
using AbstractAlgebra
using Nemo

include("lie/types.jl")
include("lie/bracket.jl")
include("lie/change_of_basis.jl")
include("lie/subspace.jl")
include("lie/center.jl")
include("lie/series.jl")
include("lie/killing.jl")
include("lie/radical.jl")
include("lie/quotient.jl")
include("lie/levi.jl")
include("lie/ideal_decomp.jl")
include("lie/derivations.jl")
include("lie/analyze.jl")
include("parametric/specialize.jl")

export LieAlgebra, LieAlgebraElem, dim, coefficient_ring, base_ring, structure_constants
export parameters, domain_denominators
export coords, basis_elem
export lie_bracket, lie_bracket!, bracket, ad, check_jacobi, check_antisymmetry
export JacobiCertificate, JacobiResidual
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
export specialize, generic_algebra, analyze_generic

end # module ParametricLie
