# v0.2.3 — Assumption-aware rank / nullspace / rref / solve (lazy, generic-first).

"""
    BranchBudget

Limits for conditional linear algebra. Exceeding a limit yields an incomplete leaf
(no guessing, no silent truncation).
"""
mutable struct BranchBudget
    max_branches::Int
    max_depth::Int
    max_expression_size::Int
    branches_used::Int
end

function BranchBudget(; max_branches::Int = 64, max_depth::Int = 16, max_expression_size::Int = 10_000)
    return BranchBudget(max_branches, max_depth, max_expression_size, 0)
end

"""
    PivotCertificate{P}

Record of an explicit split on a parameter-dependent pivot.
"""
struct PivotCertificate{P}
    poly::P
    row::Int
    col::Int
    branch::Symbol  # :zero or :nonzero
end

function Base.show(io::IO, c::PivotCertificate)
    rel = c.branch === :zero ? "=" : "≠"
    print(io, "PivotCertificate($(c.poly) $rel 0 at ($(c.row),$(c.col)))")
end

"""
    CondLeaf

A completed (or incomplete) leaf of a conditional computation.
`invariants` holds Lie-algebra results attached by [`refine`](@ref).
LA payloads (`rank`, `nullspace`, `rref`, `solution`) refer to the last matrix
operation that produced this leaf (when applicable).
"""
struct CondLeaf{P}
    sigma::AssumptionSet{P}
    rank::Union{Int,Nothing}
    nullspace::Any
    rref::Any
    solution::Any   # nothing | (:inconsistent,) | (:unique, x) | (:affine, x0, N)
    trail::Vector{PivotCertificate{P}}
    complete::Bool
    stop_reason::Union{Nothing,String}
    invariants::Dict{Symbol,Any}
end

function CondLeaf(
    sigma::AssumptionSet{P};
    rank = nothing,
    nullspace = nothing,
    rref = nothing,
    solution = nothing,
    trail = PivotCertificate{P}[],
    complete::Bool = true,
    stop_reason = nothing,
    invariants::Dict{Symbol,Any} = Dict{Symbol,Any}(),
) where {P}
    return CondLeaf{P}(
        sigma, rank, nullspace, rref, solution, trail, complete, stop_reason, invariants,
    )
end

function Base.show(io::IO, leaf::CondLeaf)
    st = leaf.complete ? "complete" : "incomplete($(leaf.stop_reason))"
    r = leaf.rank === nothing ? "?" : string(leaf.rank)
    keys = isempty(leaf.invariants) ? "" : ", inv=$(join(string.(sort(collect(Base.keys(leaf.invariants)))), ","))"
    print(io, "CondLeaf($st, rank=$r$keys, Σ=$(leaf.sigma))")
end

"""
    CondResult

Flat list of leaves from a conditional matrix operation or a refined tree.
"""
struct CondResult{P}
    leaves::Vector{CondLeaf{P}}
    incomplete::Bool
end

function Base.show(io::IO, ::MIME"text/plain", res::CondResult)
    n = length(res.leaves)
    println(io, "CondResult with $n leaf$(n == 1 ? "" : "s")",
        res.incomplete ? " (incomplete)" : "")
    for (i, leaf) in enumerate(res.leaves)
        println(io, "  [$i] ", leaf)
    end
end

Base.length(res::CondResult) = length(res.leaves)
Base.iterate(res::CondResult, args...) = iterate(res.leaves, args...)
leaves(res::CondResult) = res.leaves

# --- materialize & GE helpers -------------------------------------------------

function _materialize(A::MatElem{C}, Σ::AssumptionSet) where {C}
    B = deepcopy(A)
    m, n = nrows(B), ncols(B)
    z = zero(AbstractAlgebra.base_ring(B))
    for i in 1:m, j in 1:n
        if status(Σ, B[i, j]) === PIVOT_ZERO
            B[i, j] = z
        end
    end
    return B
end

function _expr_too_big(x, limit::Int)
    pay = _pivot_payload(x)
    pay[1] !== :poly && return false
    p = pay[2]
    try
        return length(p) > limit
    catch
        return false
    end
end

function _pick_candidate(cands)
    best = nothing
    best_score = Inf
    for c in cands
        s = pivot_complexity(c.entry)
        if s < best_score
            best_score = s
            best = c
        end
    end
    return best
end

function _eliminate_column!(A::MatElem{C}, prow::Int, col::Int) where {C}
    F = AbstractAlgebra.base_ring(A)
    m = nrows(A)
    n = ncols(A)
    piv = A[prow, col]
    iszero(piv) && return
    invp = inv(piv)
    for j in 1:n
        A[prow, j] = A[prow, j] * invp
    end
    for i in 1:m
        i == prow && continue
        factor = A[i, col]
        iszero(factor) && continue
        for j in 1:n
            A[i, j] = A[i, j] - factor * A[prow, j]
        end
        A[i, col] = zero(F)
    end
    A[prow, col] = one(F)
    return A
end

struct _Cand
    row::Int
    col::Int
    entry
end

# --- finish helpers -----------------------------------------------------------

function _extract_solve(Aaug::MatElem{C}, Σ::AssumptionSet, nvars::Int) where {C}
    B = _materialize(Aaug, Σ)
    F = AbstractAlgebra.base_ring(B)
    m = nrows(B)
    # Consistency: row with all-zero coeffs but nonzero rhs
    for i in 1:m
        allz = true
        for j in 1:nvars
            if status(Σ, B[i, j]) !== PIVOT_ZERO && !iszero(B[i, j])
                allz = false
                break
            end
        end
        if allz && status(Σ, B[i, nvars + 1]) !== PIVOT_ZERO && !iszero(B[i, nvars + 1])
            return nothing, (:inconsistent,)
        end
    end
    # Use AA nullspace / solve on materialized system A x = b
    A = zero_matrix(F, m, nvars)
    b = zero_matrix(F, m, 1)
    for i in 1:m
        for j in 1:nvars
            A[i, j] = B[i, j]
        end
        b[i, 1] = B[i, nvars + 1]
    end
    ok, sol = AbstractAlgebra.Solve.can_solve_with_solution(A, b; side = :right)
    if !ok
        return Int(rank(A)), (:inconsistent,)
    end
    _nu, N = nullspace(A)
    r = Int(rank(A))
    if ncols(N) == 0
        return r, (:unique, sol)
    end
    return r, (:affine, sol, N)
end

function _finish_leaf!(
    leaves_out,
    A::MatElem{C},
    Σ::AssumptionSet{P},
    trail,
    mode::Symbol,
    nvars::Int,
    invariants::Dict{Symbol,Any},
) where {C,P}
    if mode === :rank
        B = _materialize(A, Σ)
        r = Int(rank(B))
        push!(leaves_out, CondLeaf(Σ; rank = r, trail = copy(trail), invariants = copy(invariants)))
    elseif mode === :nullspace
        B = _materialize(A, Σ)
        _nu, N = nullspace(B)
        r = Int(rank(B))
        push!(leaves_out, CondLeaf(Σ; rank = r, nullspace = N, trail = copy(trail),
            invariants = copy(invariants)))
    elseif mode === :rref
        B = _materialize(A, Σ)
        r, R = rref(B)
        push!(leaves_out, CondLeaf(Σ; rank = Int(r), rref = R, trail = copy(trail),
            invariants = copy(invariants)))
    elseif mode === :solve
        r, sol = _extract_solve(A, Σ, nvars)
        push!(leaves_out, CondLeaf(Σ; rank = r, solution = sol, trail = copy(trail),
            invariants = copy(invariants)))
    else
        error("unknown conditional LA mode $mode")
    end
    return nothing
end

# --- core recursion -----------------------------------------------------------

function _conditional_la(
    A::MatElem{C},
    Σ::AssumptionSet{P};
    mode::Symbol = :nullspace,
    nvars::Int = ncols(A),
    budget::BranchBudget = BranchBudget(),
    trail::Vector{PivotCertificate{P}} = PivotCertificate{P}[],
    invariants::Dict{Symbol,Any} = Dict{Symbol,Any}(),
    depth::Int = 0,
) where {C,P}
    leaves_out = CondLeaf{P}[]
    _cond_la_rec!(leaves_out, deepcopy(A), Σ, budget, trail, mode, nvars, invariants, depth)
    incomplete = any(!l.complete for l in leaves_out)
    return CondResult{P}(leaves_out, incomplete)
end

function _cond_la_rec!(
    leaves_out,
    A::MatElem{C},
    Σ::AssumptionSet{P},
    budget::BranchBudget,
    trail::Vector{PivotCertificate{P}},
    mode::Symbol,
    nvars::Int,
    invariants::Dict{Symbol,Any},
    depth::Int,
) where {C,P}
    if depth > budget.max_depth
        push!(leaves_out, CondLeaf(Σ; trail = copy(trail), complete = false,
            stop_reason = "max_depth", invariants = copy(invariants)))
        return
    end

    m, n = nrows(A), ncols(A)
    # For solve, only pivot in the coefficient part (not the rhs column)
    pivot_nmax = mode === :solve ? nvars : n
    row_has_pivot = falses(m)
    pivot_cols = Int[]

    while true
        nz_cands = _Cand[]
        un_cands = _Cand[]
        oversized = false

        for col in 1:pivot_nmax
            col in pivot_cols && continue
            for row in 1:m
                row_has_pivot[row] && continue
                x = A[row, col]
                if _expr_too_big(x, budget.max_expression_size)
                    oversized = true
                    break
                end
                st = status(Σ, x)
                st === PIVOT_ZERO && continue
                c = _Cand(row, col, x)
                if st === PIVOT_NONZERO
                    push!(nz_cands, c)
                else
                    push!(un_cands, c)
                end
            end
            oversized && break
        end

        if oversized
            push!(leaves_out, CondLeaf(Σ; trail = copy(trail), complete = false,
                stop_reason = "max_expression_size", rank = length(pivot_cols),
                invariants = copy(invariants)))
            return
        end

        if !isempty(nz_cands)
            c = _pick_candidate(nz_cands)
            _eliminate_column!(A, c.row, c.col)
            row_has_pivot[c.row] = true
            push!(pivot_cols, c.col)
            continue
        end

        if !isempty(un_cands)
            if budget.branches_used >= budget.max_branches
                push!(leaves_out, CondLeaf(Σ; trail = copy(trail), complete = false,
                    stop_reason = "max_branches", rank = length(pivot_cols),
                    invariants = copy(invariants)))
                return
            end
            c = _pick_candidate(un_cands)
            pay = _pivot_payload(c.entry)
            pay[1] === :poly || error("internal: UNKNOWN pivot without poly")
            p = pay[2]

            Σ_nz = assume_nonzero(Σ, p)
            if Σ_nz !== nothing
                budget.branches_used += 1
                trail_nz = copy(trail)
                push!(trail_nz, PivotCertificate{P}(p, c.row, c.col, :nonzero))
                _cond_la_rec!(
                    leaves_out, deepcopy(A), Σ_nz, budget, trail_nz,
                    mode, nvars, invariants, depth + 1,
                )
            end

            Σ_z = assume_zero(Σ, p)
            if Σ_z !== nothing
                budget.branches_used += 1
                trail_z = copy(trail)
                push!(trail_z, PivotCertificate{P}(p, c.row, c.col, :zero))
                _cond_la_rec!(
                    leaves_out, deepcopy(A), Σ_z, budget, trail_z,
                    mode, nvars, invariants, depth + 1,
                )
            end
            return
        end

        _finish_leaf!(leaves_out, A, Σ, trail, mode, nvars, invariants)
        return
    end
end

function _budget_kwargs(; max_branches = 64, max_depth = 16, max_expression_size = 10_000)
    return BranchBudget(; max_branches, max_depth, max_expression_size)
end

function _default_sigma(M::MatElem, sigma)
    F = AbstractAlgebra.base_ring(M)
    return sigma === nothing ? empty_assumptions(F) : sigma
end

"""
    conditional_rank(M; sigma, ...) -> CondResult
"""
function conditional_rank(
    M::MatElem;
    sigma::Union{AssumptionSet,Nothing} = nothing,
    max_branches::Int = 64,
    max_depth::Int = 16,
    max_expression_size::Int = 10_000,
    budget::Union{BranchBudget,Nothing} = nothing,
    trail = nothing,
    invariants::Dict{Symbol,Any} = Dict{Symbol,Any}(),
)
    Σ = _default_sigma(M, sigma)
    b = budget === nothing ? _budget_kwargs(; max_branches, max_depth, max_expression_size) : budget
    P = typeof_param(Σ)
    tr = trail === nothing ? PivotCertificate{P}[] : trail
    return _conditional_la(M, Σ; mode = :rank, budget = b, trail = tr, invariants)
end

"""
    conditional_nullspace(M; sigma, ...) -> CondResult
"""
function conditional_nullspace(
    M::MatElem;
    sigma::Union{AssumptionSet,Nothing} = nothing,
    max_branches::Int = 64,
    max_depth::Int = 16,
    max_expression_size::Int = 10_000,
    budget::Union{BranchBudget,Nothing} = nothing,
    trail = nothing,
    invariants::Dict{Symbol,Any} = Dict{Symbol,Any}(),
)
    Σ = _default_sigma(M, sigma)
    b = budget === nothing ? _budget_kwargs(; max_branches, max_depth, max_expression_size) : budget
    P = typeof_param(Σ)
    tr = trail === nothing ? PivotCertificate{P}[] : trail
    return _conditional_la(M, Σ; mode = :nullspace, budget = b, trail = tr, invariants)
end

"""
    conditional_rref(M; sigma, ...) -> CondResult

Each complete leaf stores the reduced row-echelon form in `leaf.rref`.
"""
function conditional_rref(
    M::MatElem;
    sigma::Union{AssumptionSet,Nothing} = nothing,
    max_branches::Int = 64,
    max_depth::Int = 16,
    max_expression_size::Int = 10_000,
    budget::Union{BranchBudget,Nothing} = nothing,
    trail = nothing,
    invariants::Dict{Symbol,Any} = Dict{Symbol,Any}(),
)
    Σ = _default_sigma(M, sigma)
    b = budget === nothing ? _budget_kwargs(; max_branches, max_depth, max_expression_size) : budget
    P = typeof_param(Σ)
    tr = trail === nothing ? PivotCertificate{P}[] : trail
    return _conditional_la(M, Σ; mode = :rref, budget = b, trail = tr, invariants)
end

"""
    conditional_solve(A, b; sigma, ...) -> CondResult

Solve `A x = b` under assumptions. Each leaf's `solution` is one of:

- `(:inconsistent,)`
- `(:unique, x)` with `x` an `n×1` matrix
- `(:affine, x0, N)` particular solution plus nullspace basis columns
"""
function conditional_solve(
    A::MatElem{C},
    b;
    sigma::Union{AssumptionSet,Nothing} = nothing,
    max_branches::Int = 64,
    max_depth::Int = 16,
    max_expression_size::Int = 10_000,
    budget::Union{BranchBudget,Nothing} = nothing,
    trail = nothing,
    invariants::Dict{Symbol,Any} = Dict{Symbol,Any}(),
) where {C}
    F = AbstractAlgebra.base_ring(A)
    m, n = nrows(A), ncols(A)
    bv = _as_rhs_matrix(F, m, b)
    Aug = zero_matrix(F, m, n + 1)
    for i in 1:m, j in 1:n
        Aug[i, j] = A[i, j]
    end
    for i in 1:m
        Aug[i, n + 1] = bv[i, 1]
    end
    Σ = _default_sigma(A, sigma)
    bud = budget === nothing ? _budget_kwargs(; max_branches, max_depth, max_expression_size) : budget
    P = typeof_param(Σ)
    tr = trail === nothing ? PivotCertificate{P}[] : trail
    return _conditional_la(Aug, Σ; mode = :solve, nvars = n, budget = bud, trail = tr, invariants)
end

function _as_rhs_matrix(F, m::Int, b::MatElem)
    nrows(b) == m || throw(ArgumentError("rhs rows $(nrows(b)) ≠ $m"))
    ncols(b) == 1 || throw(ArgumentError("rhs must be a single column"))
    return b
end

function _as_rhs_matrix(F, m::Int, b::AbstractVector)
    length(b) == m || throw(ArgumentError("rhs length $(length(b)) ≠ $m"))
    M = zero_matrix(F, m, 1)
    for i in 1:m
        M[i, 1] = F(b[i])
    end
    return M
end
