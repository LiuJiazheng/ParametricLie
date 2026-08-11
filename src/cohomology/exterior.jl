# Exterior algebra Λ^• V of an n-dimensional F-vector space.
#
# Coordinate convention (graded lex, matching future cochain layout):
#   degree 0, then 1, …, n; within degree k, strictly increasing multi-indices
#   I = (i1 < … < ik) in lexicographic order.
# Index of I is 1-based:  offset(k) + lex_rank(I).

"""
    ExteriorAlgebra{C}

Graded-commutative algebra `Λ^• V` over a field `F`, with `dim V = n` and
dual generators `e^1, …, e^n`. Vector-space dimension is `2^n`.
"""
struct ExteriorAlgebra{C<:FieldElem}
    F::Field
    n::Int
    function ExteriorAlgebra{C}(F::Field, n::Int) where {C<:FieldElem}
        elem_type(F) === C || throw(ArgumentError("elem_type(F)=$(elem_type(F)) ≠ C=$C"))
        n >= 0 || throw(ArgumentError("ambient dimension must be nonnegative, got $n"))
        return new{C}(F, n)
    end
end

"""
    exterior_algebra(F::Field, n::Int) -> ExteriorAlgebra
    exterior_algebra(L::LieAlgebra) -> ExteriorAlgebra

Exterior algebra of an `n`-dimensional space over `F`, or of the underlying
space of a Lie algebra (field coefficients only).
"""
function exterior_algebra(F::Field, n::Int)
    return ExteriorAlgebra{elem_type(F)}(F, n)
end

function exterior_algebra(L::LieAlgebra{C}) where {C<:FieldElem}
    return exterior_algebra(coefficient_ring(L), dim(L))
end

coefficient_ring(Λ::ExteriorAlgebra) = Λ.F
base_ring(Λ::ExteriorAlgebra) = Λ.F
ambient_dim(Λ::ExteriorAlgebra) = Λ.n
dim(Λ::ExteriorAlgebra) = 2^Λ.n
dim(Λ::ExteriorAlgebra, k::Integer) = _binom(Λ.n, Int(k))

function _binom(n::Int, k::Int)
    (k < 0 || k > n) && return 0
    return binomial(n, k)
end

function _degree_offset(n::Int, k::Int)
    k <= 0 && return 0
    s = 0
    for j in 0:k-1
        s += _binom(n, j)
    end
    return s
end

"""
    multi_indices(n, k) -> Vector{Vector{Int}}

Strictly increasing `k`-tuples in `{1,…,n}`, lexicographic order.
"""
function multi_indices(n::Int, k::Int)
    k < 0 && throw(ArgumentError("k must be nonnegative, got $k"))
    k > n && return Vector{Int}[]
    k == 0 && return [Int[]]
    out = Vector{Vector{Int}}()
    c = collect(1:k)
    while true
        push!(out, copy(c))
        i = k
        while i >= 1 && c[i] == n - k + i
            i -= 1
        end
        i < 1 && break
        c[i] += 1
        for j in (i + 1):k
            c[j] = c[j - 1] + 1
        end
    end
    return out
end

multi_indices(Λ::ExteriorAlgebra, k::Integer) = multi_indices(Λ.n, Int(k))

# 1-based lex rank of increasing I among k-subsets of {1,…,n}.
function _lex_rank(n::Int, I::Vector{Int})
    k = length(I)
    k == 0 && return 1
    r = 1
    prev = 0
    for t in 1:k
        for x in (prev + 1):(I[t] - 1)
            r += _binom(n - x, k - t)
        end
        prev = I[t]
    end
    return r
end

function _unrank_lex(n::Int, k::Int, r::Int)
    k == 0 && return Int[]
    I = Vector{Int}(undef, k)
    x = 1
    for t in 1:k
        while true
            c = _binom(n - x, k - t)
            if r > c
                r -= c
                x += 1
            else
                I[t] = x
                x += 1
                break
            end
        end
    end
    return I
end

"""
    coord_index(n, I) -> Int

1-based index of multi-index `I` in the graded-lex basis of `Λ^•` (`dim = 2^n`).
"""
function coord_index(n::Int, I::Vector{Int})
    k = length(I)
    return _degree_offset(n, k) + _lex_rank(n, I)
end

coord_index(Λ::ExteriorAlgebra, I::Vector{Int}) = coord_index(Λ.n, I)

function _index_to_multi(n::Int, idx::Int)
    rem = idx
    k = 0
    while k <= n
        b = _binom(n, k)
        rem <= b && break
        rem -= b
        k += 1
    end
    k > n && throw(ArgumentError("index $idx out of range for n=$n"))
    return _unrank_lex(n, k, rem)
end

# --- elements ----------------------------------------------------------------

"""
    ExteriorElem{C}

Element of [`ExteriorAlgebra`](@ref), stored as a dense coordinate vector of
length `2^n` in graded-lex order.
"""
struct ExteriorElem{C<:FieldElem}
    parent::ExteriorAlgebra{C}
    coords::Vector{C}
    function ExteriorElem{C}(Λ::ExteriorAlgebra{C}, coords::Vector{C}) where {C<:FieldElem}
        length(coords) == dim(Λ) ||
            throw(ArgumentError("expected $(dim(Λ)) coordinates, got $(length(coords))"))
        return new{C}(Λ, coords)
    end
end

function ExteriorElem(Λ::ExteriorAlgebra{C}, coords::AbstractVector) where {C<:FieldElem}
    F = coefficient_ring(Λ)
    N = dim(Λ)
    length(coords) == N ||
        throw(ArgumentError("expected $N coordinates, got $(length(coords))"))
    v = Vector{C}(undef, N)
    for i in 1:N
        v[i] = F(coords[i])
    end
    return ExteriorElem{C}(Λ, v)
end

Base.parent(α::ExteriorElem) = α.parent
coords(α::ExteriorElem) = α.coords
coefficient_ring(α::ExteriorElem) = coefficient_ring(parent(α))
ambient_dim(α::ExteriorElem) = ambient_dim(parent(α))
dim(α::ExteriorElem) = dim(parent(α))

function Base.zero(Λ::ExteriorAlgebra{C}) where {C<:FieldElem}
    F = coefficient_ring(Λ)
    return ExteriorElem{C}(Λ, fill(zero(F), dim(Λ)))
end

Base.zero(α::ExteriorElem) = zero(parent(α))

function Base.one(Λ::ExteriorAlgebra{C}) where {C<:FieldElem}
    F = coefficient_ring(Λ)
    v = fill(zero(F), dim(Λ))
    v[1] = one(F)  # empty multi-index, degree 0
    return ExteriorElem{C}(Λ, v)
end

Base.one(α::ExteriorElem) = one(parent(α))

function Base.:(==)(α::ExteriorElem{C}, β::ExteriorElem{C}) where {C<:FieldElem}
    return α.parent === β.parent && α.coords == β.coords
end

function Base.iszero(α::ExteriorElem)
    return all(iszero, α.coords)
end

function Base.isone(α::ExteriorElem)
    iszero(α) && return false
    F = coefficient_ring(α)
    isone(α.coords[1]) || return false
    return all(i -> iszero(α.coords[i]), 2:length(α.coords))
end

function Base.show(io::IO, Λ::ExteriorAlgebra)
    print(io, "ExteriorAlgebra(dim V=$(ambient_dim(Λ)), dim Λ=$(dim(Λ)) over $(coefficient_ring(Λ)))")
end

function Base.show(io::IO, α::ExteriorElem)
    print(io, "ExteriorElem(n=$(ambient_dim(α)), coords=$(α.coords))")
end

function _check_same_exterior(α::ExteriorElem{C}, β::ExteriorElem{C}) where {C<:FieldElem}
    α.parent === β.parent || throw(ArgumentError("exterior elements have different parents"))
    return α.parent
end

# --- generators / monomials --------------------------------------------------

"""
    exterior_generator(Λ, i) -> ExteriorElem

Dual generator `e^i` (degree 1).
"""
function exterior_generator(Λ::ExteriorAlgebra{C}, i::Int) where {C<:FieldElem}
    n = ambient_dim(Λ)
    (1 <= i <= n) || throw(ArgumentError("generator index $i out of range 1:$n"))
    F = coefficient_ring(Λ)
    v = fill(zero(F), dim(Λ))
    v[coord_index(n, Int[i])] = one(F)
    return ExteriorElem{C}(Λ, v)
end

"""
    exterior_monomial(Λ, I) -> ExteriorElem

Basis monomial `e^{i1} ∧ … ∧ e^{ik}` for increasing (or unsorted) index list `I`.
"""
function exterior_monomial(Λ::ExteriorAlgebra{C}, I) where {C<:FieldElem}
    n = ambient_dim(Λ)
    F = coefficient_ring(Λ)
    J = collect(Int, I)
    if length(unique(J)) != length(J)
        return zero(Λ)
    end
    for i in J
        (1 <= i <= n) || throw(ArgumentError("index $i out of range 1:$n"))
    end
    perm = sortperm(J)
    Js = J[perm]
    s = _perm_sign(perm)
    v = fill(zero(F), dim(Λ))
    v[coord_index(n, Js)] = F(s)
    return ExteriorElem{C}(Λ, v)
end

function _perm_sign(perm::Vector{Int})
    # sign of permutation sending 1:k to perm (as a listing of positions)
    k = length(perm)
    s = 1
    seen = fill(false, k)
    for i in 1:k
        seen[i] && continue
        j = i
        len = 0
        while !seen[j]
            seen[j] = true
            j = perm[j]
            len += 1
        end
        iseven(len) && (s = -s)  # cycle length ℓ contributes (-1)^{ℓ-1}
    end
    return s
end

# --- arithmetic --------------------------------------------------------------

function Base.:+(α::ExteriorElem{C}, β::ExteriorElem{C}) where {C<:FieldElem}
    Λ = _check_same_exterior(α, β)
    return ExteriorElem{C}(Λ, α.coords .+ β.coords)
end

function Base.:-(α::ExteriorElem{C}, β::ExteriorElem{C}) where {C<:FieldElem}
    Λ = _check_same_exterior(α, β)
    return ExteriorElem{C}(Λ, α.coords .- β.coords)
end

function Base.:-(α::ExteriorElem{C}) where {C<:FieldElem}
    return ExteriorElem{C}(parent(α), map(-, α.coords))
end

function _scale_exterior(s, α::ExteriorElem{C}) where {C<:FieldElem}
    F = coefficient_ring(α)
    return ExteriorElem{C}(parent(α), [F(s) * c for c in α.coords])
end

Base.:*(s::Number, α::ExteriorElem) = _scale_exterior(s, α)
Base.:*(α::ExteriorElem, s::Number) = _scale_exterior(s, α)
Base.:*(s::C, α::ExteriorElem{C}) where {C<:FieldElem} = _scale_exterior(s, α)
Base.:*(α::ExteriorElem{C}, s::C) where {C<:FieldElem} = _scale_exterior(s, α)

# --- degree / homogeneous projection -----------------------------------------

"""
    support_degrees(α) -> Vector{Int}

Degrees with a nonzero component (sorted). The zero element has empty support.
"""
function support_degrees(α::ExteriorElem)
    n = ambient_dim(α)
    degs = Int[]
    idx = 1
    for k in 0:n
        nk = _binom(n, k)
        any(!iszero, view(α.coords, idx:idx+nk-1)) && push!(degs, k)
        idx += nk
    end
    return degs
end

function is_homogeneous(α::ExteriorElem)
    d = support_degrees(α)
    return length(d) <= 1
end

"""
    exterior_degree(α::ExteriorElem) -> Int

Unique degree of a homogeneous element. The zero element has degree 0.
Throws if `α` is inhomogeneous.
"""
function exterior_degree(α::ExteriorElem)
    d = support_degrees(α)
    isempty(d) && return 0
    length(d) == 1 || throw(ArgumentError("element is not homogeneous (degrees $d)"))
    return d[1]
end

"""
    homogeneous_part(α, k) -> ExteriorElem

Degree-`k` component of `α`.
"""
function homogeneous_part(α::ExteriorElem{C}, k::Integer) where {C<:FieldElem}
    Λ = parent(α)
    n = ambient_dim(Λ)
    kk = Int(k)
    (0 <= kk <= n) || return zero(Λ)
    F = coefficient_ring(Λ)
    v = fill(zero(F), dim(Λ))
    off = _degree_offset(n, kk)
    nk = _binom(n, kk)
    copyto!(v, off + 1, α.coords, off + 1, nk)
    return ExteriorElem{C}(Λ, v)
end

# --- wedge -------------------------------------------------------------------

function _wedge_sign_and_union(I::Vector{Int}, J::Vector{Int})
    for i in I, j in J
        i == j && return (0, Int[])
    end
    invs = 0
    for i in I, j in J
        i > j && (invs += 1)
    end
    U = sort!(vcat(I, J))
    return (iseven(invs) ? 1 : -1, U)
end

"""
    wedge(α, β) -> ExteriorElem

Associative graded-commutative product `α ∧ β`. Also available as `α * β`.
"""
function wedge(α::ExteriorElem{C}, β::ExteriorElem{C}) where {C<:FieldElem}
    Λ = _check_same_exterior(α, β)
    n = ambient_dim(Λ)
    F = coefficient_ring(Λ)
    out = fill(zero(F), dim(Λ))
    for ia in 1:length(α.coords)
        ca = α.coords[ia]
        iszero(ca) && continue
        I = _index_to_multi(n, ia)
        for ib in 1:length(β.coords)
            cb = β.coords[ib]
            iszero(cb) && continue
            s, U = _wedge_sign_and_union(I, _index_to_multi(n, ib))
            s == 0 && continue
            out[coord_index(n, U)] += F(s) * ca * cb
        end
    end
    return ExteriorElem{C}(Λ, out)
end

Base.:*(α::ExteriorElem{C}, β::ExteriorElem{C}) where {C<:FieldElem} = wedge(α, β)

# --- contraction / evaluation ------------------------------------------------

"""
    interior_product(v, α) -> ExteriorElem

Interior product `ι_v α`. `v` is a coordinate vector in `F^n` (length `n`).
Lowers degree by 1.
"""
function interior_product(v::AbstractVector, α::ExteriorElem{C}) where {C<:FieldElem}
    Λ = parent(α)
    n = ambient_dim(Λ)
    length(v) == n || throw(ArgumentError("expected length $n, got $(length(v))"))
    F = coefficient_ring(Λ)
    out = fill(zero(F), dim(Λ))
    for ia in 1:length(α.coords)
        ca = α.coords[ia]
        iszero(ca) && continue
        I = _index_to_multi(n, ia)
        k = length(I)
        k == 0 && continue
        for p in 1:k
            coeff = F(v[I[p]])
            iszero(coeff) && continue
            J = deleteat!(copy(I), p)
            s = isodd(p) ? 1 : -1  # (-1)^{p-1}
            out[coord_index(n, J)] += F(s) * coeff * ca
        end
    end
    return ExteriorElem{C}(Λ, out)
end

"""
    form_eval(α, vs) -> F

Evaluate a homogeneous `k`-form on `k` vectors (`vs` a vector of length-`n`
coordinate vectors). The zero form evaluates to 0.
"""
function form_eval(α::ExteriorElem{C}, vs::Vector) where {C<:FieldElem}
    Λ = parent(α)
    n = ambient_dim(Λ)
    F = coefficient_ring(Λ)
    k = length(vs)
    for v in vs
        length(v) == n || throw(ArgumentError("each vector must have length $n"))
    end
    if iszero(α)
        return zero(F)
    end
    is_homogeneous(α) || throw(ArgumentError("form_eval requires a homogeneous element"))
    exterior_degree(α) == k ||
        throw(ArgumentError("degree $(exterior_degree(α)) form evaluated on $k vectors"))
    s = zero(F)
    for I in multi_indices(n, k)
        c = α.coords[coord_index(n, I)]
        iszero(c) && continue
        M = AbstractAlgebra.zero_matrix(F, k, k)
        for a in 1:k, b in 1:k
            M[a, b] = F(vs[b][I[a]])
        end
        s += c * AbstractAlgebra.det(M)
    end
    return s
end
