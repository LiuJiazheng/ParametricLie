# Structural invariants (v0.1)
# center + change_of_basis implemented; other invariants still placeholders.

function derived_series end
function lower_central_series end
function is_solvable end
function is_nilpotent end
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
