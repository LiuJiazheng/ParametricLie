# Killing form and Cartan criteria (char 0).
#
# K_ij = tr(ad_ei ∘ ad_ej).  For coords-as-columns, B(u,v) = uᵀ K v.
# Cartan solvability (char 0): B(L,[L,L]) = 0 iff τ = L, where
#   τ = { x | B(x,y) = 0 ∀ y ∈ [L,L] } = nullspace(Dᵀ K),
#   D = basis_matrix(derived_algebra(L)).
# Cartan semisimplicity (char 0): K nondegenerate iff rad(B) = 0.

"""
    killing_form(L::LieAlgebra) -> MatElem

Killing form matrix `K` over `coefficient_ring(L)`, with

    K[i,j] = tr(ad(eᵢ) ∘ ad(eⱼ)).

Symmetric by construction (`K[j,i]` filled from `K[i,j]`).
"""
function killing_form(L::LieAlgebra{C}) where {C<:FieldElem}
    n = dim(L)
    F = coefficient_ring(L)
    K = zero_matrix(F, n, n)
    n == 0 && return K
    ads = [ad(L, basis_elem(L, i)) for i in 1:n]
    for i in 1:n
        for j in i:n
            tij = tr(ads[i] * ads[j])
            K[i, j] = tij
            if i != j
                K[j, i] = tij
            end
        end
    end
    return K
end

"""
    killing_rank(L::LieAlgebra) -> Int

Rank of the Killing form matrix.
"""
killing_rank(L::LieAlgebra) = Int(rank(killing_form(L)))

"""
    killing_radical(L::LieAlgebra) -> LieSubspace

Killing radical `rad(B) = { x | B(x,y) = 0 ∀ y ∈ L } = ker(K)`.

Cartan semisimplicity (char 0): `L` is semisimple iff this is `{0}`.
"""
function killing_radical(L::LieAlgebra{C}) where {C<:FieldElem}
    n = dim(L)
    F = coefficient_ring(L)
    n == 0 && return zero_space(L)
    K = killing_form(L)
    if rank(K) == 0
        return full_space(L)
    end
    _nullity, N = nullspace(K)
    return LieSubspace{C}(L, N)
end

"""
    killing_orthogonal(L::LieAlgebra, S::LieSubspace) -> LieSubspace

Orthogonal of `S` w.r.t. the Killing form:

    S^⊥ = { x ∈ L | B(x,y) = 0 ∀ y ∈ S } = nullspace(Dᵀ K),

where `D = basis_matrix(S)` (columns = basis of `S`) and `K = killing_form(L)`.
"""
function killing_orthogonal(L::LieAlgebra{C}, S::LieSubspace{C}) where {C<:FieldElem}
    parent(S) === L || throw(ArgumentError("subspace must belong to the same Lie algebra"))
    n = dim(L)
    n == 0 && return zero_space(L)
    D = basis_matrix(S)
    if ncols(D) == 0
        return full_space(L)
    end
    K = killing_form(L)
    E = transpose(D) * K   # d × n; E x = 0 ⇔ Dᵀ K x = 0
    _nullity, N = nullspace(E)
    return LieSubspace{C}(L, N)
end

"""
    cartan_orthogonal(L::LieAlgebra) -> LieSubspace

Cartan space

    τ = { x ∈ L | B(x,y) = 0 ∀ y ∈ [L,L] } = killing_orthogonal(L, [L,L]).

Computed as `nullspace(Dᵀ K)` with `D = basis_matrix(derived_algebra(L))`.
"""
cartan_orthogonal(L::LieAlgebra) = killing_orthogonal(L, derived_algebra(L))

"""
    is_cartan_solvable(L::LieAlgebra) -> Bool

Cartan's solvability criterion (characteristic 0): `B(L, [L,L]) = 0`,
equivalently `dim(cartan_orthogonal(L)) == dim(L)`.

Agrees with [`is_solvable`](@ref) (derived series) over char‑0 fields such as `QQ`.
"""
is_cartan_solvable(L::LieAlgebra) = dim(cartan_orthogonal(L)) == dim(L)

"""
    is_semisimple(L::LieAlgebra) -> Bool

Cartan's semisimplicity criterion (characteristic 0): Killing form nondegenerate,
i.e. `killing_rank(L) == dim(L)` (equivalently `dim(killing_radical(L)) == 0`).
"""
is_semisimple(L::LieAlgebra) = killing_rank(L) == dim(L)
