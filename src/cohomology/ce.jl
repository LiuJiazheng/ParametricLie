# Chevalley–Eilenberg cochain complex: module-induced d, lazy Z / B / H.
#
# C^k(g, M) ≅ Hom(Λ^k g, M), dim = m * binom(n,k).
# Coordinates: multi-indices I in lex order; for each I, m coordinates of
# ω(e_{i1},…,e_{ik}). When m = 1 this is the degree-k slice of Λ^•.
#
# Lazy prefix (cf. LieSeries): requesting degree k fills and caches d, Z, B, H
# for 0,…,k (or 0,…,n if k > n).

"""
    CohomologyGroup{C}

`H^k(g, M)` as a certificate: columns of [`basis_matrix`](@ref) are cocycle
representatives in `C^k` coordinates (a complement of `B^k` inside `Z^k`).
"""
struct CohomologyGroup{C<:FieldElem}
    parent::Any  # CEComplex{C}; filled after CEComplex is defined
    degree::Int
    basis_matrix::MatElem  # dim(C^k) × dim(H^k)
end

dim(H::CohomologyGroup) = ncols(H.basis_matrix)
basis_matrix(H::CohomologyGroup) = H.basis_matrix
Base.parent(H::CohomologyGroup) = H.parent

function Base.show(io::IO, H::CohomologyGroup)
    print(io, "H^$(H.degree)(dim=$(dim(H)))")
end

"""
    CEComplex{C}

Lazy Chevalley–Eilenberg complex of `L` with coefficients in `M`.

Construct with [`ce_complex`](@ref) (no linear algebra yet). [`cocycles`](@ref),
[`coboundaries`](@ref), [`cohomology`](@ref), and [`ce_differential`](@ref) at
degree `k` compute (and cache) all of `d, Z, B, H` for degrees `0,…,k`.
"""
mutable struct CEComplex{C<:FieldElem}
    algebra::LieAlgebra{C}
    module_::LieModule{C}
    filled::Int  # highest degree with d/Z/B/H cached; -1 = nothing
    d_cache::Vector{Union{Nothing, MatElem}}
    Z_cache::Vector{Union{Nothing, MatElem}}
    B_cache::Vector{Union{Nothing, MatElem}}
    H_cache::Vector{Union{Nothing, CohomologyGroup{C}}}
end

Base.parent(C::CEComplex) = C.algebra
coefficient_module(C::CEComplex) = C.module_
coefficient_ring(C::CEComplex) = coefficient_ring(C.algebra)

function Base.show(io::IO, C::CEComplex)
    n = dim(C.algebra)
    print(io, "CEComplex(g dim=$(n), M dim=$(dim(C.module_)), filled=$(C.filled))")
end

"""
    ce_complex(L, M = trivial_module(L)) -> CEComplex

Lazy CE complex. No matrices are built until a degree is requested.
"""
function ce_complex(L::LieAlgebra{C}, M::LieModule{C}) where {C<:FieldElem}
    parent(M) === L || throw(ArgumentError("module must be a module over L"))
    n = dim(L)
    slots = n + 1
    return CEComplex{C}(
        L,
        M,
        -1,
        Union{Nothing, MatElem}[nothing for _ in 1:slots],
        Union{Nothing, MatElem}[nothing for _ in 1:slots],
        Union{Nothing, MatElem}[nothing for _ in 1:slots],
        Union{Nothing, CohomologyGroup{C}}[nothing for _ in 1:slots],
    )
end

ce_complex(L::LieAlgebra{C}) where {C<:FieldElem} = ce_complex(L, trivial_module(L))

"""
    cochain_dim(L, M, k) -> Int
    cochain_dim(C, k) -> Int

Dimension of `C^k(g, M)`. Zero for `k < 0` or `k > dim(g)`.
"""
function cochain_dim(L::LieAlgebra, M::LieModule, k::Integer)
    kk = Int(k)
    (kk < 0 || kk > dim(L)) && return 0
    return dim(M) * binomial(dim(L), kk)
end

cochain_dim(C::CEComplex, k::Integer) = cochain_dim(C.algebra, C.module_, k)

# --- CE differential on coordinates -----------------------------------------

function _cochain_block(ω::Vector{C}, m::Int, rank::Int) where {C}
    off = (rank - 1) * m
    return view(ω, off + 1:off + m)
end

# ω(e_t, e_{K1}, …) with K increasing; 0 if t ∈ K.
function _alt_front(ω::Vector{C}, n::Int, k::Int, m::Int, t::Int, K::Vector{Int}) where {C}
    z = C[zero(ω[1]) for _ in 1:m]
    for κ in K
        t == κ && return z
    end
    args = Vector{Int}(undef, k)
    args[1] = t
    args[2:k] = K
    perm = sortperm(args)
    J = args[perm]
    s = _perm_sign(perm)
    blk = _cochain_block(ω, m, _lex_rank(n, J))
    if s == 1
        return C[blk[i] for i in 1:m]
    else
        return C[-blk[i] for i in 1:m]
    end
end

"""
Apply `d^k: C^k → C^{k+1}` to a coordinate vector of length `m * binom(n,k)`.
"""
function _apply_ce_d(
    L::LieAlgebra{C}, M::LieModule{C}, k::Int, ω::Vector{C}
) where {C<:FieldElem}
    n = dim(L)
    m = dim(M)
    F = coefficient_ring(L)
    a = structure_constants(L)
    nk1 = m * binomial(n, k + 1)
    out = fill(zero(F), nk1)
    k < 0 && return out
    isempty(ω) && return out
    Is = multi_indices(n, k + 1)
    for (rI, I) in enumerate(Is)
        val = fill(zero(F), m)
        # ∑_p (-1)^{p-1} e_{I[p]} · ω(I without p)
        for p in 1:k + 1
            J = deleteat!(copy(I), p)
            s = isodd(p) ? one(F) : -one(F)
            ωJ = C[_cochain_block(ω, m, _lex_rank(n, J))[t] for t in 1:m]
            av = act(M, I[p], ωJ)
            for t in 1:m
                val[t] += s * av[t]
            end
        end
        # ∑_{p<q} (-1)^{p+q} ω([e_{Ip}, e_{Iq}], e_{I without p,q})
        if k >= 1
            for p in 1:k + 1, q in (p + 1):k + 1
                s2 = iseven(p + q) ? one(F) : -one(F)
                K = deleteat!(deleteat!(copy(I), q), p)
                ip, iq = I[p], I[q]
                for t in 1:n
                    c = a[ip, iq, t]
                    iszero(c) && continue
                    ωalt = _alt_front(ω, n, k, m, t, K)
                    for u in 1:m
                        val[u] += s2 * c * ωalt[u]
                    end
                end
            end
        end
        off = (rI - 1) * m
        for u in 1:m
            out[off + u] = val[u]
        end
    end
    return out
end

function _ce_differential_matrix(L::LieAlgebra{C}, M::LieModule{C}, k::Int) where {C<:FieldElem}
    F = coefficient_ring(L)
    ncols_d = cochain_dim(L, M, k)
    nrows_d = cochain_dim(L, M, k + 1)
    D = zero_matrix(F, nrows_d, ncols_d)
    (ncols_d == 0 || nrows_d == 0) && return D
    ω = fill(zero(F), ncols_d)
    for j in 1:ncols_d
        ω[j] = one(F)
        dω = _apply_ce_d(L, M, k, ω)
        ω[j] = zero(F)
        for i in 1:nrows_d
            D[i, j] = dω[i]
        end
    end
    return D
end

function _image_basis(A::MatElem)
    F = AbstractAlgebra.base_ring(A)
    n = nrows(A)
    c = ncols(A)
    (n == 0 || c == 0) && return zero_matrix(F, n, 0)
    At = transpose(A)
    r, R = rref(At)
    B = zero_matrix(F, n, r)
    for j in 1:r, i in 1:n
        B[i, j] = R[j, i]
    end
    return B
end

function _kernel_basis(A::MatElem)
    F = AbstractAlgebra.base_ring(A)
    c = ncols(A)
    c == 0 && return zero_matrix(F, 0, 0)
    if nrows(A) == 0
        return identity_matrix(F, c)
    end
    _nullity, N = nullspace(A)
    return N
end

"""
Complement of the column space of `B` inside the column space of `Z`
(assumes `im B ⊆ im Z`). Columns of the result are a basis of `Z/B`.
"""
function _cochain_quotient_basis(Z::MatElem, B::MatElem)
    F = AbstractAlgebra.base_ring(Z)
    n = nrows(Z)
    nrows(B) == n || throw(ArgumentError("Z and B must have the same number of rows"))
    z = ncols(Z)
    z == 0 && return zero_matrix(F, n, 0)
    kept = Vector{elem_type(F)}[]
    for j in 1:z
        v = [Z[i, j] for i in 1:n]
        span = if isempty(kept)
            B
        elseif ncols(B) == 0
            _matrix_from_coord_cols(F, n, kept)
        else
            hcat(B, _matrix_from_coord_cols(F, n, kept))
        end
        if ncols(span) == 0
            all(iszero, v) || push!(kept, v)
        elseif !_in_column_span(span, v)
            push!(kept, v)
        end
    end
    return _matrix_from_coord_cols(F, n, kept)
end

function _zero_cohomology(C::CEComplex{Celt}, k::Int) where {Celt<:FieldElem}
    F = coefficient_ring(C)
    nk = cochain_dim(C, k)
    return CohomologyGroup{Celt}(C, k, zero_matrix(F, nk, 0))
end

function _fill_degree!(C::CEComplex{Celt}, k::Int) where {Celt<:FieldElem}
    L = C.algebra
    M = C.module_
    F = coefficient_ring(C)
    idx = k + 1
    D = _ce_differential_matrix(L, M, k)
    C.d_cache[idx] = D
    Z = _kernel_basis(D)
    C.Z_cache[idx] = Z
    B = if k == 0
        zero_matrix(F, cochain_dim(C, 0), 0)
    else
        _image_basis(C.d_cache[k])  # d^{k-1}
    end
    C.B_cache[idx] = B
    Hmat = _cochain_quotient_basis(Z, B)
    C.H_cache[idx] = CohomologyGroup{Celt}(C, k, Hmat)
    return C
end

"""
    _ensure_degree!(C, k)

Fill caches for degrees `0,…,min(k,n)` (prefix evaluation, like `LieSeries`).
"""
function _ensure_degree!(C::CEComplex, k::Int)
    k < 0 && throw(ArgumentError("cohomology degree must be ≥ 0, got $k"))
    n = dim(C.algebra)
    target = min(k, n)
    while C.filled < target
        _fill_degree!(C, C.filled + 1)
        C.filled += 1
    end
    return C
end

# --- public accessors -------------------------------------------------------

"""
    ce_differential(C, k) -> MatElem

Matrix of `d^k: C^k → C^{k+1}` (fills the prefix through `k`).
"""
function ce_differential(C::CEComplex, k::Integer)
    kk = Int(k)
    kk < 0 && throw(ArgumentError("degree must be ≥ 0, got $kk"))
    n = dim(C.algebra)
    if kk > n
        F = coefficient_ring(C)
        return zero_matrix(F, 0, 0)
    end
    _ensure_degree!(C, kk)
    return C.d_cache[kk + 1]
end

"""
    cocycles(C, k) -> MatElem

Basis matrix of `Z^k = ker d^k` (columns). Fills degrees `0,…,k`.
"""
function cocycles(C::CEComplex, k::Integer)
    kk = Int(k)
    kk < 0 && throw(ArgumentError("degree must be ≥ 0, got $kk"))
    n = dim(C.algebra)
    if kk > n
        F = coefficient_ring(C)
        return zero_matrix(F, 0, 0)
    end
    _ensure_degree!(C, kk)
    return C.Z_cache[kk + 1]
end

"""
    coboundaries(C, k) -> MatElem

Basis matrix of `B^k = im d^{k-1}` (columns). `B^0 = 0`. Fills degrees `0,…,k`.
"""
function coboundaries(C::CEComplex, k::Integer)
    kk = Int(k)
    kk < 0 && throw(ArgumentError("degree must be ≥ 0, got $kk"))
    n = dim(C.algebra)
    if kk > n
        F = coefficient_ring(C)
        return zero_matrix(F, 0, 0)
    end
    _ensure_degree!(C, kk)
    return C.B_cache[kk + 1]
end

"""
    cohomology(C::CEComplex, k) -> CohomologyGroup
    cohomology(L, k)
    cohomology(L, M, k)

`H^k = Z^k / B^k` with cocycle representatives. Requesting `k` computes and
caches `Z^j`, `B^j`, `H^j` for all `j ≤ k` (up to `dim(L)`).
"""
function cohomology(C::CEComplex{Celt}, k::Integer) where {Celt<:FieldElem}
    kk = Int(k)
    kk < 0 && throw(ArgumentError("degree must be ≥ 0, got $kk"))
    n = dim(C.algebra)
    if kk > n
        _ensure_degree!(C, n)
        return _zero_cohomology(C, kk)
    end
    _ensure_degree!(C, kk)
    return C.H_cache[kk + 1]
end

function cohomology(L::LieAlgebra{C}, M::LieModule{C}, k::Integer) where {C<:FieldElem}
    return cohomology(ce_complex(L, M), k)
end

function cohomology(L::LieAlgebra{C}, k::Integer) where {C<:FieldElem}
    return cohomology(L, trivial_module(L), k)
end
