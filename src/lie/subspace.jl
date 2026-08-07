# Linear subspaces of a Lie algebra (certificates via basis matrices).

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

# --- membership / complement / ideals ---------------------------------------

"""
    _in_column_span(B, v) -> Bool

Whether coordinate vector `v` lies in the column space of `B` (exact).
"""
function _in_column_span(B::MatElem, v::AbstractVector)
    n = nrows(B)
    length(v) == n || return false
    F = AbstractAlgebra.base_ring(B)
    if ncols(B) == 0
        return all(iszero, v)
    end
    vv = matrix(F, n, 1, [F(v[i]) for i in 1:n])
    ok, _ = AbstractAlgebra.Solve.can_solve_with_solution(B, vv; side = :right)
    return ok
end

function Base.in(x::LieAlgebraElem{C}, S::LieSubspace{C}) where {C<:FieldElem}
    parent(x) === parent(S) || return false
    return _in_column_span(S.basis_matrix, x.coords)
end

function Base.in(v::AbstractVector, S::LieSubspace{C}) where {C<:FieldElem}
    return _in_column_span(S.basis_matrix, v)
end

"""
    complement(L::LieAlgebra, S::LieSubspace) -> LieSubspace

A vector-space complement of `S` in `L`: the (standard) orthogonal of `S`,
i.e. `nullspace(Dᵀ)` where `D = basis_matrix(S)`. Not in general a subalgebra.
"""
function complement(L::LieAlgebra{C}, S::LieSubspace{C}) where {C<:FieldElem}
    parent(S) === L || throw(ArgumentError("subspace must belong to L"))
    n = dim(L)
    n == 0 && return zero_space(L)
    D = basis_matrix(S)
    ncols(D) == 0 && return full_space(L)
    ncols(D) == n && return zero_space(L)
    _nullity, N = nullspace(transpose(D))
    return LieSubspace{C}(L, N)
end

complement(S::LieSubspace) = complement(parent(S), S)

"""
    is_subalgebra(L::LieAlgebra, S::LieSubspace) -> Bool

Whether `[S,S] ⊆ S`.
"""
function is_subalgebra(L::LieAlgebra{C}, S::LieSubspace{C}) where {C<:FieldElem}
    parent(S) === L || throw(ArgumentError("subspace must belong to L"))
    for a in basis_elems(S), b in basis_elems(S)
        lie_bracket(L, a, b) ∈ S || return false
    end
    return true
end

is_subalgebra(S::LieSubspace) = is_subalgebra(parent(S), S)

"""
    is_ideal(L::LieAlgebra, S::LieSubspace) -> Bool

Whether `[L,S] ⊆ S`.
"""
function is_ideal(L::LieAlgebra{C}, S::LieSubspace{C}) where {C<:FieldElem}
    parent(S) === L || throw(ArgumentError("subspace must belong to L"))
    n = dim(L)
    for i in 1:n, s in basis_elems(S)
        lie_bracket(L, basis_elem(L, i), s) ∈ S || return false
    end
    return true
end

is_ideal(S::LieSubspace) = is_ideal(parent(S), S)

"""
    _matrix_from_coord_cols(F, n, cols) -> MatElem

Build an `n × length(cols)` matrix with the given coordinate columns.
"""
function _matrix_from_coord_cols(F::Field, n::Int, cols::Vector{Vector{C}}) where {C<:FieldElem}
    m = length(cols)
    B = zero_matrix(F, n, m)
    for j in 1:m, i in 1:n
        B[i, j] = cols[j][i]
    end
    return B
end

"""
    _complement_in(L, V, W) -> MatElem

Basis matrix (`n × (dim V − dim W)`) of a complement of `W` inside `V`
(assumes `W ⊆ V`).
"""
function _complement_in(
    L::LieAlgebra{C}, V::LieSubspace{C}, W::LieSubspace{C}
) where {C<:FieldElem}
    parent(V) === L && parent(W) === L ||
        throw(ArgumentError("subspaces must belong to L"))
    F = coefficient_ring(L)
    n = dim(L)
    BW = basis_matrix(W)
    kept = Vector{C}[]
    for v in basis_elems(V)
        span_cols = if isempty(kept)
            BW
        else
            ncols(BW) == 0 ? _matrix_from_coord_cols(F, n, kept) :
                hcat(BW, _matrix_from_coord_cols(F, n, kept))
        end
        if ncols(span_cols) == 0
            all(iszero, v.coords) || push!(kept, copy(v.coords))
        elseif !_in_column_span(span_cols, v.coords)
            push!(kept, copy(v.coords))
        end
    end
    return _matrix_from_coord_cols(F, n, kept)
end

"""
    _u_coords(v, BU, BW) -> Vector

Coordinates of `v` along columns of `BU` in the direct sum decomposition
`span(BU) ⊕ span(BW)` (requires `v ∈ span(BU) + span(BW)`).
"""
function _u_coords(v::AbstractVector{C}, BU::MatElem, BW::MatElem) where {C<:FieldElem}
    F = AbstractAlgebra.base_ring(BU)
    n = nrows(BU)
    d = ncols(BU)
    d == 0 && return C[]
    vv = matrix(F, n, 1, [F(v[i]) for i in 1:n])
    P = ncols(BW) == 0 ? BU : hcat(BU, BW)
    ok, sol = AbstractAlgebra.Solve.can_solve_with_solution(P, vv; side = :right)
    ok || throw(ArgumentError("vector is not in the span of [BU|BW]"))
    return C[sol[t, 1] for t in 1:d]
end
