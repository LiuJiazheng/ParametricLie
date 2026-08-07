# Derivations Der(L).
#
# A linear map D: L → L is a derivation iff it satisfies the Leibniz rule
#   D([x,y]) = [D(x), y] + [x, D(y)].
# In the standard basis, D is an n×n matrix (columns = images of basis vectors).
# The identity on all basis pairs becomes a homogeneous linear system in the n²
# entries of D; a basis of Der(L) is the nullspace of that system.

"""
    Derivations{C}

Certificate for the derivation algebra `Der(L)`: a basis of `n×n` matrices over
`coefficient_ring(L)` in the current coordinates of `L`.

The span of [`basis_matrices`](@ref) is exactly `Der(L)`; these matrices generate
`Der` as a vector space (and as a Lie algebra under the commutator of endomorphisms).
"""
struct Derivations{C<:FieldElem}
    parent::LieAlgebra{C}
    matrices::Vector{MatElem}  # each n×n
end

Base.parent(D::Derivations) = D.parent
dim(D::Derivations) = length(D.matrices)
basis_matrices(D::Derivations) = D.matrices

function Base.show(io::IO, D::Derivations)
    print(io, "Derivations(dim=$(dim(D)) of LieAlgebra(dim=$(dim(parent(D)))))")
end

"""
Column-major index of matrix entry `D[p,q]` among `n²` unknowns (1-based).
"""
@inline _der_unknown_index(n::Int, p::Int, q::Int) = (q - 1) * n + p

"""
    _derivation_equation_matrix(L) -> MatElem

Build the `(n³)×(n²)` matrix `E` such that `vec(D)` (column-major) is a
derivation iff `E * vec(D) = 0`.

For each triple `(i,j,m)`, the Leibniz identity on `[eᵢ,eⱼ]` at coordinate `m`
reads

    ∑ₖ a[i,j,k] D[m,k] − ∑ᵣ a[r,j,m] D[r,i] − ∑ₛ a[i,s,m] D[s,j] = 0.
"""
function _derivation_equation_matrix(L::LieAlgebra{C}) where {C<:FieldElem}
    n = dim(L)
    F = coefficient_ring(L)
    a = structure_constants(L)
    E = zero_matrix(F, n * n * n, n * n)
    n == 0 && return E

    for i in 1:n, j in 1:n, m in 1:n
        row = ((i - 1) * n + (j - 1)) * n + m

        # ∑ₖ a[i,j,k] D[m,k]
        for k in 1:n
            aik = a[i, j, k]
            iszero(aik) && continue
            E[row, _der_unknown_index(n, m, k)] += aik
        end

        # − ∑ᵣ a[r,j,m] D[r,i]
        for r in 1:n
            arjm = a[r, j, m]
            iszero(arjm) && continue
            E[row, _der_unknown_index(n, r, i)] -= arjm
        end

        # − ∑ₛ a[i,s,m] D[s,j]
        for s in 1:n
            aism = a[i, s, m]
            iszero(aism) && continue
            E[row, _der_unknown_index(n, s, j)] -= aism
        end
    end
    return E
end

"""
Reshape a length-`n²` coordinate vector (column-major) into an `n×n` matrix.
"""
function _unvec_derivation(F::Field, n::Int, v::AbstractVector)
    D = zero_matrix(F, n, n)
    for q in 1:n, p in 1:n
        D[p, q] = F(v[_der_unknown_index(n, p, q)])
    end
    return D
end

"""
    is_derivation(L::LieAlgebra, D::MatElem) -> Bool

Whether the `n×n` matrix `D` satisfies the Leibniz identity on all basis pairs.
"""
function is_derivation(L::LieAlgebra{C}, D::MatElem) where {C<:FieldElem}
    n = dim(L)
    size(D) == (n, n) || return false
    F = coefficient_ring(L)
    for i in 1:n, j in 1:n
        ei = basis_elem(L, i)
        ej = basis_elem(L, j)
        br = lie_bracket(L, ei, ej)
        lhs = D * matrix(F, n, 1, br.coords)
        Dei = [D[p, i] for p in 1:n]
        Dej = [D[p, j] for p in 1:n]
        rhs1 = lie_bracket(L, Dei, ej.coords)
        rhs2 = lie_bracket(L, ei.coords, Dej)
        for m in 1:n
            lhs[m, 1] == rhs1[m] + rhs2[m] || return false
        end
    end
    return true
end

"""
    apply_derivation(L, D, x) -> LieAlgebraElem

Apply derivation matrix `D` to an element: coordinates `D * coords(x)`.
"""
function apply_derivation(
    L::LieAlgebra{C}, D::MatElem, x::LieAlgebraElem{C}
) where {C<:FieldElem}
    x.parent === L || throw(ArgumentError("element must belong to L"))
    size(D) == (dim(L), dim(L)) || throw(ArgumentError("D must be n×n"))
    F = coefficient_ring(L)
    ycol = D * matrix(F, dim(L), 1, x.coords)
    return LieAlgebraElem{C}(L, [ycol[i, 1] for i in 1:dim(L)])
end

function apply_derivation(L::LieAlgebra{C}, D::MatElem, x::AbstractVector) where {C<:FieldElem}
    return apply_derivation(L, D, LieAlgebraElem(L, x))
end

"""
    _matrix_in_span(D, Ms) -> Bool

Whether matrix `D` lies in the F-span of the matrices `Ms` (exact).
"""
function _matrix_in_span(D::MatElem, Ms::Vector)
    isempty(Ms) && return _is_zero_matrix_der(D)
    F = AbstractAlgebra.base_ring(D)
    n = nrows(D)
    d = length(Ms)
    # Stack vec(M_j) as columns; solve A ξ = vec(D)
    A = zero_matrix(F, n * n, d)
    for j in 1:d, q in 1:n, p in 1:n
        A[_der_unknown_index(n, p, q), j] = Ms[j][p, q]
    end
    rhs = zero_matrix(F, n * n, 1)
    for q in 1:n, p in 1:n
        rhs[_der_unknown_index(n, p, q), 1] = D[p, q]
    end
    ok, _ = AbstractAlgebra.Solve.can_solve_with_solution(A, rhs; side = :right)
    return ok
end

function _is_zero_matrix_der(T::MatElem)
    for i in 1:nrows(T), j in 1:ncols(T)
        iszero(T[i, j]) || return false
    end
    return true
end

"""
    derivations(L::LieAlgebra) -> Derivations

Compute a basis of `Der(L)` by solving the Leibniz linear system on matrix
entries of `D` in the current basis of `L`.

Returns a [`Derivations`](@ref) certificate: `dim(derivations(L))` matrices that
span `Der(L)`. Special cases:
- abelian `L`: `Der(L) ≅ 𝔤𝔩ₙ`, dimension `n²`;
- semisimple (char 0): `Der(L) = Inn(L) ≅ L / Z(L)`.
"""
function derivations(L::LieAlgebra{C}) where {C<:FieldElem}
    n = dim(L)
    F = coefficient_ring(L)

    if n == 0
        return Derivations{C}(L, MatElem[])
    end

    # Abelian: every endomorphism is a derivation
    a = structure_constants(L)
    if all(iszero, a)
        mats = MatElem[]
        for q in 1:n, p in 1:n
            Dmat = zero_matrix(F, n, n)
            Dmat[p, q] = one(F)
            push!(mats, Dmat)
        end
        return Derivations{C}(L, mats)
    end

    E = _derivation_equation_matrix(L)
    _nullity, N = nullspace(E)
    d = ncols(N)
    mats = MatElem[
        _unvec_derivation(F, n, [N[r, j] for r in 1:(n * n)]) for j in 1:d
    ]
    return Derivations{C}(L, mats)
end
