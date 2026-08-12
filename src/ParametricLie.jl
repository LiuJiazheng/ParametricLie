"""
    ParametricLie

Exact symbolic toolkit for finite-dimensional Lie algebras given by structure
constants: invariants, parameterized families / stratification, Chevalley–Eilenberg
cohomology, central extensions, and formal deformations. Complements LieGroups.jl;
does not replace group numerics.

Layout:

- `src/lie/` — single-algebra core (structure, invariants, `analyze`, isomorphism)
- `src/parametric/` — families, conditional LA, stratification, jump causes
- `src/cohomology/` — exterior algebra, modules, CE complex, central extensions
- `src/deformations/` — NR bracket, truncated Maurer–Cartan, gauge equivalence

Docs: `README.md`, `docs/README.md`.
"""
module ParametricLie

using LinearAlgebra
using SparseArrays
using AbstractAlgebra
using Nemo

# Include order = dependency order.
# --- lie (single algebra) ----------------------------------------------------
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
include("lie/isomorphism.jl")

# --- parametric (families & stratification) ----------------------------------
include("parametric/specialize.jl")
include("parametric/assumptions.jl")
include("parametric/conditional_linalg.jl")
include("parametric/cond_tree.jl")
include("parametric/stratify.jl")
include("parametric/strata_compare.jl")

# --- cohomology (exterior algebra, modules, CE complex, central extensions) --
include("cohomology/exterior.jl")
include("cohomology/module.jl")
include("cohomology/ce.jl")
include("cohomology/central_extension.jl")

# --- deformations (NR, MC, gauge equivalence) --------------------------------
include("deformations/nr_bracket.jl")
include("deformations/maurer_cartan.jl")
include("deformations/gauge.jl")

# jump explanation depends on stratify + deformations + isomorphism
include("parametric/jump_explain.jl")

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
export IsoCertificate, isomorphism, isomorphic, is_structure_isomorphism
export specialize, generic_algebra, analyze_generic
export PivotStatus, PIVOT_ZERO, PIVOT_NONZERO, PIVOT_UNKNOWN
export AssumptionSet, status, normalize_poly, assume_zero, assume_nonzero
export empty_assumptions, assumptions_from_domain
export PivotCertificate, BranchBudget, CondLeaf, CondResult, CondTree
export conditional_rank, conditional_nullspace, conditional_rref, conditional_solve
export conditional_center, conditional_killing, conditional_invariants
export cond_tree, refine, refine_leaves, leaves, algebra
export center_dim, center_basis, killing_rank_of, killing_radical_dim
export derived_dim_of, derived_profile_of, radical_dim_of, der_dim_of
export is_solvable_of, is_nilpotent_of
export analyze_conditional, default_conditional_suite, is_complete
export unresolved_invariants, invariant_signature
export Stratum, JumpEntry, JumpReport, Stratification
export stratify, jump_table
export StratumComparison, compare, validate_stratum, validate_stratification
export JumpCause, explain_jump, explain_jumps!, explain_jumps
export sample_stratum_point, wall_cocycle
export ExteriorAlgebra, ExteriorElem, exterior_algebra, ambient_dim
export multi_indices, coord_index, exterior_generator, exterior_monomial
export wedge, homogeneous_part, support_degrees, is_homogeneous, exterior_degree
export interior_product, form_eval
export LieModule, trivial_module, adjoint_module, action_matrices, act
export CEComplex, ce_complex, coefficient_module, cochain_dim
export ce_differential, cocycles, coboundaries, cohomology, CohomologyGroup
export central_extension, is_trivial_cocycle
export adjoint_bracket_cochain, nr_circle, nr_bracket
export FormalDeformation, MCCertificate, formal_deformation, extend!
export deformation_term, mc_rhs, obstruction_cochain, mc_certificate
export is_integrable, stalled_at, max_order, filled_order, seed, is_rigid
export equivalent, equivalent_with_gauge
export gauge_normal_form, gauge_normal_form!, gauge_normal_form_with_gauge

end # module ParametricLie
