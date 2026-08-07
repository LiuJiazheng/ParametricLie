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
