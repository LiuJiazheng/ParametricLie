# Ideal decomposition of semisimple Lie algebras via the adjoint commutant.
#
# Pipeline (characteristic 0, exact over F):
#   Aᵢ = ad(eᵢ)
#   → M stacks (Aᵢᵀ ⊗ I − I ⊗ Aᵢ)
#   → C = { T | T Aᵢ = Aᵢ T ∀i }  ≅  unvec(ker M)
#   → nontrivial idempotent P ∈ C  (via minpoly factorization)
#   → s = im(P) ⊕ ker(P)  (both ideals)
#   → recurse until blocks are simple (dim C = 1)

"""
    adjoint_commutant(L::LieAlgebra) -> Vector{MatElem}

Basis of the commutant of the adjoint representation:

    C = { T ∈ Mₙ(F) | T ad(eᵢ) = ad(eᵢ) T  for all basis vectors eᵢ }.

Computed as the nullspace of the stacked Kronecker system

    M = [ A₁ᵀ⊗I − I⊗A₁ ; … ; Aₙᵀ⊗I − I⊗Aₙ ],

with column-major `vec(T)`. For semisimple `L` over a char‑0 field, `dim(C)`
equals the number of simple ideals whenever each factor has endomorphism ring `F`.
"""
function adjoint_commutant(L::LieAlgebra{C}) where {C<:FieldElem}
    n = dim(L)
    F = coefficient_ring(L)
    n == 0 && return MatElem[]

    I_n = identity_matrix(F, n)
    # Stack n blocks of size n²×n² → n³×n²
    M = zero_matrix(F, n * n * n, n * n)
    for i in 1:n
        Ai = ad(L, basis_elem(L, i))
        Mi = kronecker_product(transpose(Ai), I_n) - kronecker_product(I_n, Ai)
        r0 = (i - 1) * n * n
        for r in 1:(n*n), c in 1:(n*n)
            M[r0+r, c] = Mi[r, c]
        end
    end

    _nu, N = nullspace(M)
    d = ncols(N)
    Ts = MatElem[zero_matrix(F, n, n) for _ in 1:d]
    for j in 1:d
        Ts[j] = _unvec_matrix(F, n, view_col(N, j))
    end
    return Ts
end

"""
Column `j` of matrix `N` as a dense vector of length `nrows(N)`.
"""
function view_col(N::MatElem, j::Int)
    F = AbstractAlgebra.base_ring(N)
    return [N[i, j] for i in 1:nrows(N)]
end

"""
Column-major unvec: index `(p,q)` ↔ `(q-1)*n + p`.
"""
function _unvec_matrix(F::Field, n::Int, v::AbstractVector)
    T = zero_matrix(F, n, n)
    length(v) == n * n || throw(ArgumentError("expected length n²=$(n*n), got $(length(v))"))
    for q in 1:n, p in 1:n
        T[p, q] = F(v[(q-1)*n+p])
    end
    return T
end

function _is_scalar_matrix(T::MatElem)
    n = nrows(T)
    n == 0 && return true
    F = AbstractAlgebra.base_ring(T)
    λ = T[1, 1]
    for i in 1:n, j in 1:n
        expected = i == j ? λ : zero(F)
        T[i, j] == expected || return false
    end
    return true
end

function _is_zero_matrix(T::MatElem)
    for i in 1:nrows(T), j in 1:ncols(T)
        iszero(T[i, j]) || return false
    end
    return true
end

function _is_identity_matrix(T::MatElem)
    n = nrows(T)
    F = AbstractAlgebra.base_ring(T)
    for i in 1:n, j in 1:n
        expected = i == j ? one(F) : zero(F)
        T[i, j] == expected || return false
    end
    return true
end

"""
Evaluate a univariate polynomial at a square matrix (Horner), returning an
`n×n` matrix over the same coefficient field.
"""
function _poly_at_matrix(p, T::MatElem)
    F = AbstractAlgebra.base_ring(T)
    n = nrows(T)
    S = zero_matrix(F, n, n)
    # Horner: ((a_d x + a_{d-1}) x + … )
    for k in degree(p):-1:0
        S = S * T
        ak = coeff(p, k)
        for i in 1:n
            S[i, i] += ak
        end
    end
    return S
end

"""
Integer coefficient tuples for generic elements of an `d`-dimensional space.
"""
function _generic_coeff_tuples(d::Int)
    tuples = Vector{Vector{Int}}()
    # Basis directions
    for i in 1:d
        v = zeros(Int, d)
        v[i] = 1
        push!(tuples, v)
    end
    d == 1 && return tuples
    # All-ones and arithmetic progression
    push!(tuples, ones(Int, d))
    push!(tuples, collect(1:d))
    # A few mixed small combinations
    if d >= 2
        v = zeros(Int, d)
        v[1] = 1
        v[2] = 1
        push!(tuples, v)
        v2 = zeros(Int, d)
        v2[1] = 1
        v2[2] = -1
        push!(tuples, v2)
    end
    if d >= 3
        push!(tuples, [1, 2, 3, zeros(Int, d - 3)...])
    end
    return tuples
end

"""
    _splitting_idempotent(Ts) -> Union{MatElem, Nothing}

Find a nontrivial idempotent `P` in `span(Ts)` with `P² = P`, `P ≠ 0`, `P ≠ I`,
by trying generic linear combinations and splitting the minimal polynomial into
coprime factors (CRT / spectral projection).

Returns `nothing` if no split is found (block treated as simple over `F`).
"""
function _splitting_idempotent(Ts::Vector)
    isempty(Ts) && return nothing
    T0 = Ts[1]
    F = AbstractAlgebra.base_ring(T0)
    n = nrows(T0)
    d = length(Ts)
    n <= 1 && return nothing

    R, x = polynomial_ring(F, :x; cached = false)

    for coeffs in _generic_coeff_tuples(d)
        T = zero_matrix(F, n, n)
        for j in 1:d
            cj = F(coeffs[j])
            iszero(cj) && continue
            T += cj * Ts[j]
        end
        _is_scalar_matrix(T) && continue

        μ = minpoly(R, T)
        degree(μ) <= 1 && continue

        fac = factor(μ)
        factors = [p for (p, _e) in fac]
        length(factors) < 2 && continue

        # Split first irreducible vs the rest
        f = factors[1]
        g = divexact(μ, f)   # μ = f * g, gcd(f,g)=1 when μ square-free
        # Semisimple commutant ⇒ minpoly of generic element is square-free
        δ, a, b = gcdx(f, g)
        isone(δ) || continue

        # P = (b g)(T)  — idempotent, im(P) = ker(f(T))
        P = _poly_at_matrix(b * g, T)
        if !_is_zero_matrix(P) && !_is_identity_matrix(P)
            # Numerical check P² = P
            if P * P == P
                return P
            end
        end
    end
    return nothing
end

"""
Induce the Lie algebra structure of an ideal/subalgebra `S ⊆ L` on coordinates
of a basis of `S`. Returns `(L_S, embed)` where `embed` is `n×dim(S)` and
ambient coords = `embed * coords_S`.
"""
function _induced_lie_algebra(
    L::LieAlgebra{C}, S::LieSubspace{C}
) where {C<:FieldElem}
    parent(S) === L || throw(ArgumentError("subspace must belong to L"))
    F = coefficient_ring(L)
    n = dim(L)
    d = dim(S)
    B = basis_matrix(S)
    f = fill(zero(F), d, d, d)
    if d == 0
        return LieAlgebra(F, f), B
    end
    lifts = [C[B[i, j] for i in 1:n] for j in 1:d]
    for a in 1:d, b in 1:d
        br = lie_bracket(L, lifts[a], lifts[b])
        rhs = matrix(F, n, 1, br)
        ok, sol = AbstractAlgebra.Solve.can_solve_with_solution(B, rhs; side = :right)
        ok || throw(ArgumentError("bracket left the subspace; not a subalgebra"))
        for c in 1:d
            f[a, b, c] = sol[c, 1]
        end
    end
    return LieAlgebra(F, f), B
end

"""
Map a subspace of the induced algebra (coords in `F^d`) back to ambient `L`
via the embedding matrix `embed` (`n×d`).
"""
function _embed_subspace(
    L::LieAlgebra{C}, embed::MatElem, S_ind::LieSubspace{C}
) where {C<:FieldElem}
    Bind = basis_matrix(S_ind)
    if ncols(Bind) == 0
        return zero_space(L)
    end
    B = embed * Bind
    return LieSubspace{C}(L, B)
end

function _image_as_subspace(Lind::LieAlgebra{C}, P::MatElem) where {C<:FieldElem}
    F = coefficient_ring(Lind)
    n = dim(Lind)
    cols = Vector{C}[]
    for j in 1:n
        v = C[P[i, j] for i in 1:n]
        all(iszero, v) || push!(cols, v)
    end
    return _column_span(Lind, cols)
end

function _kernel_as_subspace(Lind::LieAlgebra{C}, P::MatElem) where {C<:FieldElem}
    _nu, N = nullspace(P)
    return LieSubspace{C}(Lind, N)
end

function _ideal_decomp_block(
    L::LieAlgebra{C}, S::LieSubspace{C}
) where {C<:FieldElem}
    d = dim(S)
    d == 0 && return LieSubspace{C}[]
    d == 1 && return LieSubspace{C}[S]  # 1-dim semisimple ⇒ simple (over char ≠ 2,3 quirks aside)

    Lind, embed = _induced_lie_algebra(L, S)
    Ts = adjoint_commutant(Lind)
    if length(Ts) <= 1
        return LieSubspace{C}[S]
    end

    P = _splitting_idempotent(Ts)
    if P === nothing
        return LieSubspace{C}[S]
    end

    Im = _image_as_subspace(Lind, P)
    Ker = _kernel_as_subspace(Lind, P)
    # Degenerate split — bail out as simple
    if dim(Im) == 0 || dim(Ker) == 0 || dim(Im) + dim(Ker) != d
        return LieSubspace{C}[S]
    end

    Im_amb = _embed_subspace(L, embed, Im)
    Ker_amb = _embed_subspace(L, embed, Ker)
    return vcat(_ideal_decomp_block(L, Im_amb), _ideal_decomp_block(L, Ker_amb))
end

"""
    ideal_decomposition(L::LieAlgebra) -> Vector{LieSubspace}

Decompose a **semisimple** Lie algebra into simple ideals (characteristic 0).

Uses the adjoint-commutant method: solve `T ad(x) = ad(x) T`, find a nontrivial
idempotent in the commutant by minimal-polynomial splitting, and recurse on
`im(P) ⊕ ker(P)`.

Each returned [`LieSubspace`](@ref) is an ideal of `L`; they form a direct sum
equal to `L`. Requires [`is_semisimple`](@ref)`(L)`.
"""
function ideal_decomposition(L::LieAlgebra{C}) where {C<:FieldElem}
    is_semisimple(L) ||
        throw(ArgumentError("ideal_decomposition requires a semisimple Lie algebra"))
    dim(L) == 0 && return LieSubspace{C}[]
    return _ideal_decomp_block(L, full_space(L))
end

"""
    ideal_decomposition(lev::LeviDecomposition) -> Vector{LieSubspace}

Decompose the semisimple Levi factor / quotient algebra of a Levi decomposition.
"""
function ideal_decomposition(lev::LeviDecomposition{C}) where {C<:FieldElem}
    return ideal_decomposition(lev.quotient)
end

"""
    is_simple(L::LieAlgebra) -> Bool

Whether `L` is simple: semisimple with adjoint commutant of dimension 1
(only scalar endomorphisms).
"""
function is_simple(L::LieAlgebra{C}) where {C<:FieldElem}
    is_semisimple(L) || return false
    dim(L) == 0 && return false
    return length(adjoint_commutant(L)) == 1
end
