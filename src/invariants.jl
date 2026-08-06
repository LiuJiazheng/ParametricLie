# Structural invariants (v0.1)
# center, series, change_of_basis; killing/derivations still placeholders.

function killing_form end
function derivations end

# --- subspaces / center -----------------------------------------------------

"""
    LieSubspace{C}

A linear subspace of a [`LieAlgebra`](@ref), given by a basis matrix whose
**columns** are coordinates in the ambient `F^n`.
"""
struct LieSubspace{C<:FieldElem}
    parent::LieAlgebra{C}
    basis_matrix::MatElem  # n × d over coefficient_ring(parent)
end

Base.parent(S::LieSubspace) = S.parent
dim(S::LieSubspace) = ncols(S.basis_matrix)
basis_matrix(S::LieSubspace) = S.basis_matrix

function basis_elems(S::LieSubspace{C}) where {C<:FieldElem}
    L = parent(S)
    n = dim(L)
    d = dim(S)
    B = S.basis_matrix
    return [LieAlgebraElem{C}(L, [B[i, j] for i in 1:n]) for j in 1:d]
end

function Base.show(io::IO, S::LieSubspace)
    print(io, "LieSubspace(dim=$(dim(S)) ⊆ LieAlgebra(dim=$(dim(parent(S)))))")
end

"""
    full_space(L::LieAlgebra) -> LieSubspace

The whole ambient space `F^n` as a subspace (standard basis).
"""
function full_space(L::LieAlgebra{C}) where {C<:FieldElem}
    n = dim(L)
    F = coefficient_ring(L)
    B = n == 0 ? zero_matrix(F, 0, 0) : identity_matrix(F, n)
    return LieSubspace{C}(L, B)
end

"""
    zero_space(L::LieAlgebra) -> LieSubspace

The zero subspace.
"""
function zero_space(L::LieAlgebra{C}) where {C<:FieldElem}
    n = dim(L)
    F = coefficient_ring(L)
    return LieSubspace{C}(L, zero_matrix(F, n, 0))
end

"""
    _column_span(L, vectors) -> LieSubspace

Exact column-space basis of the given coordinate vectors (via `rref` on rows).
"""
function _column_span(L::LieAlgebra{C}, vectors::Vector{Vector{C}}) where {C<:FieldElem}
    n = dim(L)
    F = coefficient_ring(L)
    # drop exact zeros
    gens = filter(v -> !all(iszero, v), vectors)
    if isempty(gens)
        return zero_space(L)
    end
    m = length(gens)
    A = zero_matrix(F, m, n)  # rows = generators
    for i in 1:m, j in 1:n
        A[i, j] = gens[i][j]
    end
    r, R = rref(A)
    B = zero_matrix(F, n, r)
    for j in 1:r, i in 1:n
        B[i, j] = R[j, i]
    end
    return LieSubspace{C}(L, B)
end

"""
    commutator_span(L, A::LieSubspace, B::LieSubspace) -> LieSubspace

Subspace bracket
`[A,B] = span{ [a_i, b_j] }` over bases of `A` and `B` (exact column reduction).
"""
function commutator_span(L::LieAlgebra{C}, A::LieSubspace{C}, B::LieSubspace{C}) where {C<:FieldElem}
    parent(A) === L && parent(B) === L ||
        throw(ArgumentError("subspaces must belong to L"))
    vecs = Vector{C}[]
    for a in basis_elems(A), b in basis_elems(B)
        push!(vecs, lie_bracket(L, a, b).coords)
    end
    return _column_span(L, vecs)
end

commutator_span(A::LieSubspace{C}, B::LieSubspace{C}) where {C<:FieldElem} =
    commutator_span(parent(A), A, B)

"""
    _center_equation_matrix(L) -> MatElem

Build the `n² × n` matrix `E` such that `z = ∑ a_j e_j` is central iff `E * a = 0`,
where row `(i-1)*n + k` encodes `[e_i, z]_k = ∑_j a[i,j,k] a_j`.
"""
function _center_equation_matrix(L::LieAlgebra{C}) where {C<:FieldElem}
    n = dim(L)
    F = coefficient_ring(L)
    a = structure_constants(L)
    E = zero_matrix(F, n * n, n)
    for i in 1:n, k in 1:n
        row = (i - 1) * n + k
        for j in 1:n
            E[row, j] = a[i, j, k]
        end
    end
    return E
end

"""
    center(L::LieAlgebra) -> LieSubspace

Center `Z(L) = { z | [x,z] = 0 ∀ x }`, computed as the exact nullspace of the
structure-constant equation matrix over `coefficient_ring(L)`.

Returns a [`LieSubspace`](@ref) whose columns form a basis of `Z` (certificate).
For an abelian algebra the equation matrix is zero and `dim(center(L)) == dim(L)`.
"""
function center(L::LieAlgebra{C}) where {C<:FieldElem}
    n = dim(L)
    F = coefficient_ring(L)

    # Empty algebra
    if n == 0
        return LieSubspace{C}(L, zero_matrix(F, 0, 0))
    end

    # Abelian / all brackets zero: ker = whole F^n (avoid relying only on nullspace edge cases)
    a = structure_constants(L)
    if all(iszero, a)
        return LieSubspace{C}(L, identity_matrix(F, n))
    end

    E = _center_equation_matrix(L)
    # AbstractAlgebra/Nemo: nullspace(E) -> (nullity, N) with E*N == 0 exactly
    _nullity, N = nullspace(E)
    return LieSubspace{C}(L, N)
end

# --- derived / lower central series (lazy) ----------------------------------

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

# --- change of basis --------------------------------------------------------
#
# Convention: columns of P are the new basis vectors in old coordinates,
#   e'_j = ∑_i e_i P[i,j].
# Then old coords ↔ new coords by  x_old = P * x_new,  x_new = P^{-1} * x_old.
# Structure constants a[i,j,k] = a^k_{ij} transform as the mixed tensor:
#   a'[p,q,r] = ∑_{i,j,k} (P^{-1})[r,k] * a[i,j,k] * P[i,p] * P[j,q].

function _matrix_n(F::Field, P, n::Int)
    if P isa AbstractAlgebra.MatElem
        size(P) == (n, n) || throw(ArgumentError("P must be $n×$n, got $(size(P))"))
        M = zero_matrix(F, n, n)
        for i in 1:n, j in 1:n
            M[i, j] = F(P[i, j])
        end
        return M
    elseif P isa AbstractMatrix
        size(P) == (n, n) || throw(ArgumentError("P must be $n×$n, got $(size(P))"))
        M = zero_matrix(F, n, n)
        for i in 1:n, j in 1:n
            M[i, j] = F(P[i, j])
        end
        return M
    else
        throw(ArgumentError("P must be an n×n matrix over the coefficient field"))
    end
end

"""
    change_of_basis(L::LieAlgebra, P) -> LieAlgebra

Return a Lie algebra isomorphic to `L` whose structure constants are expressed
in the new basis defined by `P`.

**Convention.** Columns of `P` are new basis vectors in old coordinates:

    e'_j = ∑_i e_i P[i,j]

so `x_old = P * x_new`. With `a[i,j,k] = a^k_{ij}`,

    a'[p,q,r] = ∑_{i,j,k} (P^{-1})[r,k] * a[i,j,k] * P[i,p] * P[j,q].

`P` may be an AbstractAlgebra matrix or a Julia `AbstractMatrix` of coefficients.
"""
function change_of_basis(L::LieAlgebra{C}, P) where {C<:FieldElem}
    n = dim(L)
    F = coefficient_ring(L)
    M = _matrix_n(F, P, n)
    is_unit(det(M)) || throw(ArgumentError("change-of-basis matrix P must be invertible"))
    Minv = inv(M)

    a = structure_constants(L)
    a_new = fill(zero(F), n, n, n)
    for p in 1:n, q in 1:n, r in 1:n
        s = zero(F)
        for i in 1:n, j in 1:n, k in 1:n
            aik = a[i, j, k]
            iszero(aik) && continue
            s += Minv[r, k] * aik * M[i, p] * M[j, q]
        end
        a_new[p, q, r] = s
    end
    return LieAlgebra{C}(F, a_new)
end

"""
    change_of_basis(x::LieAlgebraElem, P) -> LieAlgebraElem

Express `x` in the new basis: builds `L' = change_of_basis(parent(x), P)` and
returns the element with coordinates `P^{-1} * coords(x)`.
"""
function change_of_basis(x::LieAlgebraElem{C}, P) where {C<:FieldElem}
    L = parent(x)
    n = dim(L)
    F = coefficient_ring(L)
    M = _matrix_n(F, P, n)
    is_unit(det(M)) || throw(ArgumentError("change-of-basis matrix P must be invertible"))
    Minv = inv(M)
    Lnew = change_of_basis(L, M)

    # x_new = P^{-1} x_old
    xcol = matrix(F, n, 1, x.coords)
    ycol = Minv * xcol
    coords_new = [ycol[i, 1] for i in 1:n]
    return LieAlgebraElem{C}(Lnew, coords_new)
end
