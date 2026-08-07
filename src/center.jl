# Center Z(L) = { z | [x,z] = 0 ∀ x }.

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
