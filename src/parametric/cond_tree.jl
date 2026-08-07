# v0.2.3 — Conditional computation tree + incremental invariant refine.
#
# Semantic model: one shared tree; each invariant pass reads leaf Σ, attaches
# results, and subdivides only that leaf when an UNKNOWN pivot blocks progress.

"""
    CondTree

Assumption-aware computation state for a parametric Lie algebra (generic Frac
view). Leaves carry accumulated `Σ`, pivot trails, and invariant bags.
"""
struct CondTree{P,C<:FieldElem}
    algebra::LieAlgebra{C}
    result::CondResult{P}
    budget::BranchBudget
end

leaves(T::CondTree) = leaves(T.result)
Base.length(T::CondTree) = length(T.result)
algebra(T::CondTree) = T.algebra

"""
    cond_tree(L; max_branches, max_depth, max_expression_size, sigma) -> CondTree

Root tree: single leaf with domain-denominator nonvanishings (or given `sigma`).
"""
function cond_tree(
    L::LieAlgebra;
    sigma = nothing,
    max_branches::Int = 64,
    max_depth::Int = 16,
    max_expression_size::Int = 10_000,
)
    Lg = generic_algebra(L)
    F = coefficient_ring(Lg)
    Σ = if sigma !== nothing
        sigma
    else
        assumptions_from_domain(F, collect(domain_denominators(Lg)))
    end
    P = typeof_param(Σ)
    leaf = CondLeaf(Σ)
    budget = BranchBudget(; max_branches, max_depth, max_expression_size)
    return CondTree{P,elem_type(F)}(Lg, CondResult{P}([leaf], false), budget)
end

"""
    refine_leaves(res, f!) -> CondResult

Map over complete leaves; `f!(invs, leaf)` mutates a copy of `leaf.invariants`.
"""
function refine_leaves(res::CondResult{P}, f!) where {P}
    out = CondLeaf{P}[]
    for leaf in res.leaves
        if !leaf.complete
            push!(out, leaf)
            continue
        end
        invs = copy(leaf.invariants)
        f!(invs, leaf)
        push!(out, CondLeaf{P}(
            leaf.sigma, leaf.rank, leaf.nullspace, leaf.rref, leaf.solution,
            leaf.trail, leaf.complete, leaf.stop_reason, invs,
        ))
    end
    return CondResult{P}(out, res.incomplete)
end

# --- merge LA children into a parent leaf -------------------------------------

function _merge_la_child(parent::CondLeaf{P}, child::CondLeaf{P}, attach!) where {P}
    invs = copy(parent.invariants)
    if child.complete
        attach!(invs, child)
    end
    return CondLeaf{P}(
        child.sigma,
        child.rank,
        child.nullspace,
        child.rref,
        child.solution,
        child.trail,  # already includes parent trail when passed into LA
        child.complete,
        child.stop_reason,
        invs,
    )
end

function _expand_leaf_with_la(parent::CondLeaf{P}, la::CondResult{P}, attach!) where {P}
    return CondLeaf{P}[_merge_la_child(parent, ch, attach!) for ch in la.leaves]
end

# --- matrix builders ----------------------------------------------------------

function _bracket_generator_matrix(L::LieAlgebra{C}) where {C<:FieldElem}
    n = dim(L)
    F = coefficient_ring(L)
    cols = Vector{C}[]
    for i in 1:n, j in (i + 1):n
        v = lie_bracket(L, basis_elem(L, i), basis_elem(L, j)).coords
        all(iszero, v) || push!(cols, v)
    end
    if isempty(cols)
        return zero_matrix(F, n, 0)
    end
    m = length(cols)
    M = zero_matrix(F, n, m)
    for j in 1:m, i in 1:n
        M[i, j] = cols[j][i]
    end
    return M
end

function _commutator_generator_matrix(
    L::LieAlgebra{C}, BA::MatElem, BB::MatElem,
) where {C<:FieldElem}
    n = dim(L)
    F = coefficient_ring(L)
    dA, dB = ncols(BA), ncols(BB)
    cols = Vector{C}[]
    for j in 1:dA, k in 1:dB
        aj = LieAlgebraElem{C}(L, [BA[i, j] for i in 1:n])
        bk = LieAlgebraElem{C}(L, [BB[i, k] for i in 1:n])
        v = lie_bracket(L, aj, bk).coords
        all(iszero, v) || push!(cols, v)
    end
    if isempty(cols)
        return zero_matrix(F, n, 0)
    end
    m = length(cols)
    M = zero_matrix(F, n, m)
    for j in 1:m, i in 1:n
        M[i, j] = cols[j][i]
    end
    return M
end

"""
Column-span basis via conditional rref on generator rows (same layout as `_column_span`).
"""
function _conditional_column_span_la(
    L::LieAlgebra{C},
    gens_as_columns::MatElem{C},
    parent::CondLeaf{P},
    budget::BranchBudget,
) where {C,P}
    n = dim(L)
    F = coefficient_ring(L)
    m = ncols(gens_as_columns)
    if m == 0
        N = zero_matrix(F, n, 0)
        ch = CondLeaf(parent.sigma; rank = 0, nullspace = N, trail = copy(parent.trail),
            invariants = copy(parent.invariants))
        return CondResult{P}([ch], false)
    end
    # rows = generators (m × n), as in `_column_span`
    A = zero_matrix(F, m, n)
    for i in 1:m, j in 1:n
        A[i, j] = gens_as_columns[j, i]
    end
    return conditional_rref(
        A;
        sigma = parent.sigma,
        budget,
        trail = copy(parent.trail),
        invariants = copy(parent.invariants),
    )
end

function _basis_from_rref(L::LieAlgebra{C}, leaf::CondLeaf) where {C}
    n = dim(L)
    F = coefficient_ring(L)
    R = leaf.rref
    r = leaf.rank === nothing ? 0 : leaf.rank
    B = zero_matrix(F, n, r)
    for j in 1:r, i in 1:n
        B[i, j] = R[j, i]
    end
    return B
end

# --- per-invariant leaf expanders ---------------------------------------------

function _refine_center_leaf(L::LieAlgebra{C}, leaf::CondLeaf{P}, budget::BranchBudget) where {C,P}
    n = dim(L)
    F = coefficient_ring(L)
    if n == 0
        N = zero_matrix(F, 0, 0)
        invs = copy(leaf.invariants)
        invs[:center_dim] = 0
        invs[:center_basis] = N
        return [CondLeaf(leaf.sigma; trail = copy(leaf.trail), invariants = invs)]
    end
    if all(iszero, structure_constants(L))
        N = identity_matrix(F, n)
        invs = copy(leaf.invariants)
        invs[:center_dim] = n
        invs[:center_basis] = N
        return [CondLeaf(leaf.sigma; trail = copy(leaf.trail), invariants = invs)]
    end
    E = _center_equation_matrix(L)
    la = conditional_nullspace(
        E; sigma = leaf.sigma, budget, trail = copy(leaf.trail),
        invariants = copy(leaf.invariants),
    )
    return _expand_leaf_with_la(leaf, la, (invs, ch) -> begin
        N = ch.nullspace
        invs[:center_dim] = N === nothing ? nothing : ncols(N)
        invs[:center_basis] = N
    end)
end

function _refine_killing_rank_leaf(L::LieAlgebra{C}, leaf::CondLeaf{P}, budget::BranchBudget) where {C,P}
    n = dim(L)
    if n == 0
        invs = copy(leaf.invariants)
        invs[:killing_rank] = 0
        return [CondLeaf(leaf.sigma; trail = copy(leaf.trail), invariants = invs)]
    end
    K = killing_form(L)
    la = conditional_rank(
        K; sigma = leaf.sigma, budget, trail = copy(leaf.trail),
        invariants = copy(leaf.invariants),
    )
    return _expand_leaf_with_la(leaf, la, (invs, ch) -> begin
        invs[:killing_rank] = ch.rank
    end)
end

function _refine_killing_radical_leaf(L::LieAlgebra{C}, leaf::CondLeaf{P}, budget::BranchBudget) where {C,P}
    n = dim(L)
    F = coefficient_ring(L)
    if n == 0
        N = zero_matrix(F, 0, 0)
        invs = copy(leaf.invariants)
        invs[:killing_radical_dim] = 0
        invs[:killing_radical_basis] = N
        invs[:killing_rank] = 0
        return [CondLeaf(leaf.sigma; trail = copy(leaf.trail), invariants = invs)]
    end
    K = killing_form(L)
    la = conditional_nullspace(
        K; sigma = leaf.sigma, budget, trail = copy(leaf.trail),
        invariants = copy(leaf.invariants),
    )
    return _expand_leaf_with_la(leaf, la, (invs, ch) -> begin
        N = ch.nullspace
        invs[:killing_radical_basis] = N
        invs[:killing_radical_dim] = N === nothing ? nothing : ncols(N)
        invs[:killing_rank] = ch.rank
    end)
end

function _refine_derived_dim_leaf(L::LieAlgebra{C}, leaf::CondLeaf{P}, budget::BranchBudget) where {C,P}
    G = _bracket_generator_matrix(L)
    if ncols(G) == 0
        F = coefficient_ring(L)
        N = zero_matrix(F, dim(L), 0)
        invs = copy(leaf.invariants)
        invs[:derived_dim] = 0
        invs[:derived_basis] = N
        return [CondLeaf(leaf.sigma; trail = copy(leaf.trail), invariants = invs)]
    end
    la = _conditional_column_span_la(L, G, leaf, budget)
    out = CondLeaf{P}[]
    for ch in la.leaves
        invs = copy(leaf.invariants)
        if ch.complete
            B = _basis_from_rref(L, ch)
            invs[:derived_dim] = ch.rank
            invs[:derived_basis] = B
        end
        push!(out, CondLeaf{P}(
            ch.sigma, ch.rank, ch.nullspace, ch.rref, ch.solution,
            ch.trail, ch.complete, ch.stop_reason, invs,
        ))
    end
    return out
end

"""
Iteratively refine derived (or lower-central) dimension profile on a leaf,
allowing splits at each step.
"""
function _refine_series_profile_leaf(
    L::LieAlgebra{C},
    leaf::CondLeaf{P},
    budget::BranchBudget;
    kind::Symbol,   # :derived or :lower_central
    profile_key::Symbol,
    solvable_key::Union{Symbol,Nothing} = nothing,
) where {C,P}
    n = dim(L)
    F = coefficient_ring(L)
    # Worklist of (leaf, basis_of_current_term, profile_so_far)
    B0 = identity_matrix(F, n)
    work = Tuple{CondLeaf{P},MatElem,Vector{Int}}[(leaf, B0, Int[n])]
    done = CondLeaf{P}[]

    while !isempty(work)
        cur, Bold, profile = pop!(work)
        d_old = ncols(Bold)
        if d_old == 0
            invs = copy(cur.invariants)
            invs[profile_key] = profile
            if solvable_key !== nothing
                invs[solvable_key] = true
            end
            push!(done, CondLeaf(cur.sigma; trail = copy(cur.trail), invariants = invs,
                complete = cur.complete, stop_reason = cur.stop_reason))
            continue
        end

        G = if kind === :derived
            _commutator_generator_matrix(L, Bold, Bold)
        else
            I = identity_matrix(F, n)
            _commutator_generator_matrix(L, I, Bold)
        end

        if ncols(G) == 0
            invs = copy(cur.invariants)
            push!(profile, 0)
            invs[profile_key] = profile
            if solvable_key !== nothing
                invs[solvable_key] = true
            end
            if kind === :derived
                # dim[L,L] is the second profile entry when present; do not
                # overwrite with the terminal zero of the series.
                if length(profile) >= 2
                    invs[:derived_dim] = get(invs, :derived_dim, profile[2])
                else
                    invs[:derived_dim] = 0
                    invs[:derived_basis] = zero_matrix(F, n, 0)
                end
            end
            push!(done, CondLeaf(cur.sigma; trail = copy(cur.trail), invariants = invs))
            continue
        end

        la = _conditional_column_span_la(L, G, cur, budget)
        for ch in la.leaves
            if !ch.complete
                invs = copy(cur.invariants)
                invs[profile_key] = profile
                push!(done, CondLeaf{P}(
                    ch.sigma, ch.rank, ch.nullspace, ch.rref, ch.solution,
                    ch.trail, false, ch.stop_reason, invs,
                ))
                continue
            end
            Bnew = _basis_from_rref(L, ch)
            d_new = ch.rank === nothing ? 0 : ch.rank
            prof2 = copy(profile)
            push!(prof2, d_new)
            if d_new == 0 || d_new == d_old
                invs = copy(cur.invariants)
                invs[profile_key] = prof2
                if solvable_key !== nothing
                    invs[solvable_key] = (d_new == 0)
                end
                if kind === :derived && length(prof2) >= 2
                    invs[:derived_dim] = get(invs, :derived_dim, prof2[2])
                    if length(prof2) == 2
                        invs[:derived_dim] = d_new
                        invs[:derived_basis] = Bnew
                    end
                end
                push!(done, CondLeaf{P}(
                    ch.sigma, ch.rank, ch.nullspace, ch.rref, ch.solution,
                    ch.trail, true, nothing, invs,
                ))
            else
                invs = copy(cur.invariants)
                if kind === :derived && length(prof2) == 2
                    invs[:derived_dim] = d_new
                    invs[:derived_basis] = Bnew
                end
                push!(work, (CondLeaf{P}(
                    ch.sigma, ch.rank, ch.nullspace, ch.rref, ch.solution,
                    ch.trail, true, nothing, invs,
                ), Bnew, prof2))
            end
        end
    end
    return done
end

function _refine_derived_profile_leaf(L, leaf, budget)
    return _refine_series_profile_leaf(
        L, leaf, budget;
        kind = :derived,
        profile_key = :derived_profile,
        solvable_key = :is_solvable,
    )
end

function _refine_lower_central_profile_leaf(L, leaf, budget)
    return _refine_series_profile_leaf(
        L, leaf, budget;
        kind = :lower_central,
        profile_key = :lower_central_profile,
        solvable_key = :is_nilpotent,
    )
end

function _refine_solvability_leaf(L, leaf, budget)
    # Prefer reusing profile if present
    if haskey(leaf.invariants, :derived_profile)
        invs = copy(leaf.invariants)
        prof = invs[:derived_profile]
        invs[:is_solvable] = (!isempty(prof) && last(prof) == 0)
        return [CondLeaf(leaf.sigma; trail = copy(leaf.trail), invariants = invs)]
    end
    return _refine_derived_profile_leaf(L, leaf, budget)
end

function _refine_nilpotency_leaf(L, leaf, budget)
    if haskey(leaf.invariants, :lower_central_profile)
        invs = copy(leaf.invariants)
        prof = invs[:lower_central_profile]
        invs[:is_nilpotent] = (!isempty(prof) && last(prof) == 0)
        return [CondLeaf(leaf.sigma; trail = copy(leaf.trail), invariants = invs)]
    end
    return _refine_lower_central_profile_leaf(L, leaf, budget)
end

function _refine_radical_leaf(L::LieAlgebra{C}, leaf::CondLeaf{P}, budget::BranchBudget) where {C,P}
    # Need derived basis under current Σ
    leaves1 = if haskey(leaf.invariants, :derived_basis)
        CondLeaf{P}[leaf]
    else
        _refine_derived_dim_leaf(L, leaf, budget)
    end
    out = CondLeaf{P}[]
    n = dim(L)
    F = coefficient_ring(L)
    for lf in leaves1
        if !lf.complete
            push!(out, lf)
            continue
        end
        D = lf.invariants[:derived_basis]
        if D === nothing
            push!(out, lf)
            continue
        end
        if ncols(D) == 0
            invs = copy(lf.invariants)
            N = n == 0 ? zero_matrix(F, 0, 0) : identity_matrix(F, n)
            invs[:radical_dim] = n
            invs[:radical_basis] = N
            push!(out, CondLeaf(lf.sigma; trail = copy(lf.trail), invariants = invs))
            continue
        end
        K = killing_form(L)
        E = transpose(D) * K
        la = conditional_nullspace(
            E; sigma = lf.sigma, budget, trail = copy(lf.trail),
            invariants = copy(lf.invariants),
        )
        append!(out, _expand_leaf_with_la(lf, la, (invs, ch) -> begin
            N = ch.nullspace
            invs[:radical_basis] = N
            invs[:radical_dim] = N === nothing ? nothing : ncols(N)
        end))
    end
    return out
end

function _refine_derivations_leaf(L::LieAlgebra{C}, leaf::CondLeaf{P}, budget::BranchBudget) where {C,P}
    n = dim(L)
    F = coefficient_ring(L)
    if n == 0
        invs = copy(leaf.invariants)
        invs[:der_dim] = 0
        invs[:der_basis] = MatElem[]
        return [CondLeaf(leaf.sigma; trail = copy(leaf.trail), invariants = invs)]
    end
    E = _derivation_equation_matrix(L)
    la = conditional_nullspace(
        E; sigma = leaf.sigma, budget, trail = copy(leaf.trail),
        invariants = copy(leaf.invariants),
    )
    return _expand_leaf_with_la(leaf, la, (invs, ch) -> begin
        N = ch.nullspace
        if N === nothing
            invs[:der_dim] = nothing
            invs[:der_basis] = nothing
            return
        end
        d = ncols(N)
        mats = MatElem[]
        for j in 1:d
            M = zero_matrix(F, n, n)
            for q in 1:n, p in 1:n
                M[p, q] = N[_der_unknown_index(n, p, q), j]
            end
            push!(mats, M)
        end
        invs[:der_dim] = d
        invs[:der_basis] = mats
    end)
end

const _REFINE_DISPATCH = Dict{Symbol,Function}(
    :center => _refine_center_leaf,
    :killing_rank => _refine_killing_rank_leaf,
    :killing_radical => _refine_killing_radical_leaf,
    :derived_dim => _refine_derived_dim_leaf,
    :derived_profile => _refine_derived_profile_leaf,
    :lower_central_profile => _refine_lower_central_profile_leaf,
    :solvability => _refine_solvability_leaf,
    :nilpotency => _refine_nilpotency_leaf,
    :is_solvable => _refine_solvability_leaf,
    :is_nilpotent => _refine_nilpotency_leaf,
    :radical => _refine_radical_leaf,
    :derivations => _refine_derivations_leaf,
    :der => _refine_derivations_leaf,
)

"""
    refine(T::CondTree, inv::Symbol; kwargs...) -> CondTree
    refine(T::CondTree, invs::Symbol...) -> CondTree

Incrementally refine every leaf by the named invariant(s). Each pass may split
leaves further; subsequent invariants inherit the refined `Σ`.

Supported symbols (must-have / nice-to-have for v0.2.3):

| Symbol | Attaches |
|--------|----------|
| `:center` | `center_dim`, `center_basis` |
| `:killing_rank` | `killing_rank` |
| `:killing_radical` | `killing_rank`, `killing_radical_dim`, `killing_radical_basis` |
| `:derived_dim` | `derived_dim`, `derived_basis` |
| `:derived_profile` | `derived_profile`, `is_solvable`, also derived_dim/basis |
| `:solvability` / `:is_solvable` | `is_solvable` (via derived profile if needed) |
| `:nilpotency` / `:is_nilpotent` | `is_nilpotent` (via lower-central profile) |
| `:radical` | `radical_dim`, `radical_basis` |
| `:derivations` / `:der` | `der_dim`, `der_basis` |
"""
function refine(T::CondTree, invs::Symbol...)
    for inv in invs
        T = _refine_one(T, inv)
    end
    return T
end

function refine(T::CondTree, invs::AbstractVector{Symbol})
    return refine(T, invs...)
end

"""
    conditional_invariants(L, invs...; kwargs...) -> CondTree

`cond_tree(L)` followed by sequential [`refine`](@ref) passes.
"""
function conditional_invariants(L::LieAlgebra, invs::Symbol...; kwargs...)
    T = cond_tree(L; kwargs...)
    return refine(T, invs...)
end

# --- convenience one-shot clients ---------------------------------------------

"""
    conditional_center(L; ...) -> CondTree

Equivalent to `refine(cond_tree(L), :center)`.
"""
conditional_center(L::LieAlgebra; kwargs...) = refine(cond_tree(L; kwargs...), :center)

"""
    conditional_killing(L; ...) -> CondTree

Refine Killing radical (includes Killing rank).
"""
conditional_killing(L::LieAlgebra; kwargs...) =
    refine(cond_tree(L; kwargs...), :killing_radical)

# --- leaf accessors -----------------------------------------------------------

center_dim(leaf::CondLeaf) = get(leaf.invariants, :center_dim, nothing)
center_basis(leaf::CondLeaf) = get(leaf.invariants, :center_basis, nothing)
killing_rank_of(leaf::CondLeaf) = get(leaf.invariants, :killing_rank, nothing)
killing_radical_dim(leaf::CondLeaf) = get(leaf.invariants, :killing_radical_dim, nothing)
derived_dim_of(leaf::CondLeaf) = get(leaf.invariants, :derived_dim, nothing)
derived_profile_of(leaf::CondLeaf) = get(leaf.invariants, :derived_profile, nothing)
radical_dim_of(leaf::CondLeaf) = get(leaf.invariants, :radical_dim, nothing)
der_dim_of(leaf::CondLeaf) = get(leaf.invariants, :der_dim, nothing)
is_solvable_of(leaf::CondLeaf) = get(leaf.invariants, :is_solvable, nothing)
is_nilpotent_of(leaf::CondLeaf) = get(leaf.invariants, :is_nilpotent, nothing)

# --- v0.2.4: suite completeness & analyze_conditional -------------------------

"""
Keys that must be present for an invariant symbol to count as certified on a leaf.
"""
const _CERTIFIED_KEYS = Dict{Symbol,Vector{Symbol}}(
    :center => [:center_dim],
    :killing_rank => [:killing_rank],
    :killing_radical => [:killing_radical_dim],
    :derived_dim => [:derived_dim],
    :derived_profile => [:derived_profile],
    :lower_central_profile => [:lower_central_profile],
    :solvability => [:is_solvable],
    :is_solvable => [:is_solvable],
    :nilpotency => [:is_nilpotent],
    :is_nilpotent => [:is_nilpotent],
    :radical => [:radical_dim],
    :derivations => [:der_dim],
    :der => [:der_dim],
)

"""
    default_conditional_suite() -> Vector{Symbol}

Standard v0.2.4 matrix-based invariant suite.
"""
default_conditional_suite() = Symbol[
    :center, :killing_radical, :derived_profile, :nilpotency,
]

function _leaf_has_invariant(leaf::CondLeaf, inv::Symbol)
    keys = get(_CERTIFIED_KEYS, inv, Symbol[inv])
    return all(k -> haskey(leaf.invariants, k) && leaf.invariants[k] !== nothing, keys)
end

"""
    unresolved_invariants(leaf, suite) -> Vector{Symbol}

Members of `suite` not yet certified on `leaf`.
"""
function unresolved_invariants(leaf::CondLeaf, suite)
    return Symbol[inv for inv in suite if !_leaf_has_invariant(leaf, inv)]
end

"""
    is_complete(leaf; invariants=…) -> Bool

Whether `leaf` is complete for the suite (certified values, no budget stop).
"""
function is_complete(leaf::CondLeaf; invariants = default_conditional_suite())
    leaf.complete || return false
    return isempty(unresolved_invariants(leaf, invariants))
end

"""
    is_complete(T::CondTree; invariants=…) -> Bool

Every leaf is complete for the suite and the tree is not marked incomplete.
"""
function is_complete(T::CondTree; invariants = default_conditional_suite())
    T.result.incomplete && return false
    return all(leaf -> is_complete(leaf; invariants), leaves(T))
end

function _refine_one(T::CondTree{P,C}, inv::Symbol) where {P,C}
    haskey(_REFINE_DISPATCH, inv) || throw(ArgumentError(
        "unknown refine invariant :$inv; known = $(sort(collect(keys(_REFINE_DISPATCH))))",
    ))
    expander = _REFINE_DISPATCH[inv]
    L = T.algebra
    budget = T.budget
    new_leaves = CondLeaf{P}[]
    for leaf in T.result.leaves
        if !leaf.complete
            push!(new_leaves, leaf)
            continue
        end
        # Inheritance: skip recomputation when already certified on this leaf
        if _leaf_has_invariant(leaf, inv)
            push!(new_leaves, leaf)
            continue
        end
        append!(new_leaves, expander(L, leaf, budget))
    end
    incomplete = any(!l.complete for l in new_leaves)
    return CondTree{P,C}(L, CondResult{P}(new_leaves, incomplete), budget)
end

"""
    analyze_conditional(L; invariants, kwargs...) -> CondTree

Run [`cond_tree`](@ref) then sequential [`refine`](@ref) for the requested suite.
Partial suites are allowed (e.g. only `:center` and `:derived_dim`).
"""
function analyze_conditional(
    L::LieAlgebra;
    invariants = default_conditional_suite(),
    kwargs...,
)
    suite = invariants isa Symbol ? Symbol[invariants] : collect(Symbol, invariants)
    T = cond_tree(L; kwargs...)
    return refine(T, suite...)
end

"""
    invariant_signature(leaf) -> Dict{Symbol,Any}

Basis-invariant summary bag used by stratification / comparison.
"""
function invariant_signature(leaf::CondLeaf)
    d = Dict{Symbol,Any}()
    for (k, acc) in (
        (:center_dim, center_dim),
        (:derived_dim, derived_dim_of),
        (:killing_rank, killing_rank_of),
        (:killing_radical_dim, killing_radical_dim),
        (:radical_dim, radical_dim_of),
        (:der_dim, der_dim_of),
        (:is_solvable, is_solvable_of),
        (:is_nilpotent, is_nilpotent_of),
        (:derived_profile, derived_profile_of),
    )
        v = acc(leaf)
        v === nothing || (d[k] = v)
    end
    return d
end

function Base.show(io::IO, ::MIME"text/plain", T::CondTree)
    println(io, "CondTree(dim=$(dim(T.algebra)), leaves=$(length(T))",
        T.result.incomplete ? ", incomplete)" : ")")
    for (i, leaf) in enumerate(leaves(T))
        sig = invariant_signature(leaf)
        parts = String[]
        for k in (:center_dim, :derived_dim, :killing_rank, :is_solvable, :is_nilpotent)
            haskey(sig, k) && push!(parts, "$k=$(sig[k])")
        end
        st = leaf.complete ? "" : " [incomplete]"
        println(io, "  [$i] Σ=$(leaf.sigma)$st")
        isempty(parts) || println(io, "      ", join(parts, ", "))
    end
end

