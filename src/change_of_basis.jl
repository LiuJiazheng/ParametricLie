# Change of basis for structure constants and elements.
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
