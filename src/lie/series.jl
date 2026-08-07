# Derived series, lower central series, solvability / nilpotency.

"""
Series kind for [`LieSeries`](@ref).
"""
@enum SeriesKind begin
    DerivedSeriesKind      # D₁ = g, Dₖ₊₁ = [Dₖ, Dₖ]
    LowerCentralSeriesKind # C₁ = g, Cₖ₊₁ = [g, Cₖ]
end

"""
    LieSeries{C}

Lazy derived or lower-central series of a Lie algebra.

- Construct with [`derived_series`](@ref) / [`lower_central_series`](@ref)
  (no heavy work yet).
- `terms(S, k)` / `S[1:k]` compute (and cache) terms on demand.
- [`layers`](@ref)`(S)`, [`is_solvable`](@ref)`(S)`, and [`is_nilpotent`](@ref)`(S)`
  force a full evaluation and then cache:

  1. the sequence until `0` or a stable nonzero (“irreducible”) term,
  2. the layer count ([`layers`](@ref): index of first zero, or `Inf`),
  3. whether the series reaches zero (solvability / nilpotency for the matching kind).

Indexing is **1-based** with term `1` equal to the whole algebra `g`.
"""
mutable struct LieSeries{C<:FieldElem}
    parent::LieAlgebra{C}
    kind::SeriesKind
    cache::Vector{LieSubspace{C}}
    reached_zero::Bool
    stabilized::Bool                 # stopped at nonzero dim (never → 0)
    layers::Union{Nothing, Int, typeof(Inf)}  # nothing until fully evaluated
end

Base.parent(S::LieSeries) = S.parent

function Base.show(io::IO, S::LieSeries)
    k = length(S.cache)
    ly = S.layers === nothing ? "lazy" : string(S.layers)
    print(io, "LieSeries($(S.kind), cached=$(k), layers=$(ly), reached_zero=$(S.reached_zero))")
end

function _empty_series(L::LieAlgebra{C}, kind::SeriesKind) where {C<:FieldElem}
    return LieSeries{C}(L, kind, LieSubspace{C}[], false, false, nothing)
end

function _finalize_layers!(S::LieSeries)
    if S.reached_zero
        S.layers = length(S.cache)   # 1-based index of first zero term
    elseif S.stabilized
        S.layers = Inf
    end
    return S
end

function _ensure_term!(S::LieSeries{C}, k::Int) where {C<:FieldElem}
    k >= 1 || throw(ArgumentError("series index must be ≥ 1, got $k"))
    L = parent(S)

    # Seed term 1 = g
    if isempty(S.cache)
        push!(S.cache, full_space(L))
        if dim(L) == 0
            S.reached_zero = true
            _finalize_layers!(S)
        end
    end

    while length(S.cache) < k && !S.reached_zero && !S.stabilized
        prev = S.cache[end]
        if dim(prev) == 0
            S.reached_zero = true
            _finalize_layers!(S)
            break
        end
        next = if S.kind === DerivedSeriesKind
            commutator_span(L, prev, prev)
        else
            commutator_span(L, full_space(L), prev)
        end
        push!(S.cache, next)
        if dim(next) == 0
            S.reached_zero = true
            _finalize_layers!(S)
        elseif dim(next) == dim(prev)
            # Nested nonincreasing dims: no further drop ⇒ never reaches 0
            S.stabilized = true
            _finalize_layers!(S)
        end
    end
    return S
end

"""
    _evaluate!(S::LieSeries) -> LieSeries

Force full evaluation until the series hits zero or stabilizes nonzero.
Caches the terminal sequence, [`layers`](@ref), and `reached_zero`.
"""
function _evaluate!(S::LieSeries)
    S.layers !== nothing && return S
    # At most dim strict dimension drops before zero or a plateau
    _ensure_term!(S, dim(parent(S)) + 2)
    if S.layers === nothing
        # Safety: cannot drop further without having marked a terminal state
        S.stabilized = !S.reached_zero
        _finalize_layers!(S)
    end
    return S
end

"""
    terms(S::LieSeries, k::Int) -> Vector{LieSubspace}

Compute and return the first `k` terms of the series (1-based, term 1 = `g`).
If the series reaches zero (or stabilizes nonzero) earlier, the returned vector
may be shorter than `k`.
"""
function terms(S::LieSeries, k::Int)
    k >= 0 || throw(ArgumentError("k must be nonnegative"))
    k == 0 && return empty(S.cache)
    _ensure_term!(S, k)
    n = min(k, length(S.cache))
    return S.cache[1:n]
end

"""
    terms(S::LieSeries) -> Vector{LieSubspace}

Fully evaluate and return the cached sequence through `0` or the first
stable nonzero term.
"""
function terms(S::LieSeries)
    _evaluate!(S)
    return S.cache
end

function Base.getindex(S::LieSeries, k::Int)
    ts = terms(S, k)
    length(ts) >= k || throw(BoundsError(S, k))
    return ts[k]
end

function Base.getindex(S::LieSeries, r::UnitRange{Int})
    r.start >= 1 || throw(ArgumentError("series range must start at ≥ 1"))
    ts = terms(S, r.stop)
    length(ts) >= r.stop || throw(BoundsError(S, r))
    return ts[r]
end

"""
    layers(S::LieSeries) -> Int or Inf

Number of terms until the series reaches zero: the **1-based index of the first
zero term**. Forces full evaluation and caches the result.

If the series never reaches zero (stabilizes at a nonzero subspace), returns
`Inf`.

Classical derived length / nilpotency class (with `D⁰ = g`) equal `layers(S) - 1`
when finite.
"""
function layers(S::LieSeries)
    _evaluate!(S)
    return S.layers
end

"""
    derived_series(L::LieAlgebra) -> LieSeries

Lazy derived series `D₁ = g`, `Dₖ₊₁ = [Dₖ, Dₖ]`.
"""
derived_series(L::LieAlgebra{C}) where {C<:FieldElem} =
    _empty_series(L, DerivedSeriesKind)

"""
    lower_central_series(L::LieAlgebra) -> LieSeries

Lazy lower central series `C₁ = g`, `Cₖ₊₁ = [g, Cₖ]`.
"""
lower_central_series(L::LieAlgebra{C}) where {C<:FieldElem} =
    _empty_series(L, LowerCentralSeriesKind)

"""
    derived_algebra(L::LieAlgebra) -> LieSubspace

`[g,g]`, i.e. the second term of the derived series.
"""
function derived_algebra(L::LieAlgebra)
    dim(L) == 0 && return zero_space(L)
    return derived_series(L)[2]
end

"""
    is_solvable(S::LieSeries) -> Bool

For a derived series: fully evaluate and return whether some term is zero.
Caches sequence, [`layers`](@ref), and the result on `S`.
"""
function is_solvable(S::LieSeries)
    S.kind === DerivedSeriesKind ||
        throw(ArgumentError("is_solvable requires a derived series; got $(S.kind)"))
    _evaluate!(S)
    return S.reached_zero
end

"""
    is_solvable(L::LieAlgebra) -> Bool

Whether the derived series of `L` reaches zero.
"""
is_solvable(L::LieAlgebra) = is_solvable(derived_series(L))

"""
    is_nilpotent(S::LieSeries) -> Bool

For a lower-central series: fully evaluate and return whether some term is zero.
Caches sequence, [`layers`](@ref), and the result on `S`.
"""
function is_nilpotent(S::LieSeries)
    S.kind === LowerCentralSeriesKind ||
        throw(ArgumentError("is_nilpotent requires a lower-central series; got $(S.kind)"))
    _evaluate!(S)
    return S.reached_zero
end

"""
    is_nilpotent(L::LieAlgebra) -> Bool

Whether the lower central series of `L` reaches zero.
"""
is_nilpotent(L::LieAlgebra) = is_nilpotent(lower_central_series(L))
