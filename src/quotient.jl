# Semisimple quotient L / rad(L) via a vector-space complement.

"""
    QuotientAlgebra{C}

Certificate for the quotient Lie algebra `L / rad(L)`.

- `ambient` — original Lie algebra `L`
- `radical` — solvable radical as a subspace of `L`
- `complement` — chosen lifts of a basis of `L/rad(L)` (vector-space complement;
  **not** necessarily a Levi subalgebra)
- `algebra` — abstract Lie algebra on `F^{n-k}` with the induced structure constants
"""
struct QuotientAlgebra{C<:FieldElem}
    ambient::LieAlgebra{C}
    radical::LieSubspace{C}
    complement::LieSubspace{C}
    algebra::LieAlgebra{C}
end

Base.parent(Q::QuotientAlgebra) = Q.ambient
dim(Q::QuotientAlgebra) = dim(Q.algebra)
coefficient_ring(Q::QuotientAlgebra) = coefficient_ring(Q.algebra)

function Base.show(io::IO, Q::QuotientAlgebra)
    print(
        io,
        "QuotientAlgebra(dim=$(dim(Q.algebra)) = ambient=$(dim(Q.ambient)) / radical=$(dim(Q.radical)))",
    )
end

"""
    _structure_constants_mod(
        L, lifts::Vector{Vector{C}}, R::LieSubspace
    ) -> Array{C,3}

Structure constants `f[a,b,c]` of the span of `lifts` modulo `R`: write

    [x_a, x_b] = ∑_c f[a,b,c] x_c + (element of R)

in the splitting `span(lifts) ⊕ R`.
"""
function _structure_constants_mod(
    L::LieAlgebra{C}, lifts::Vector{Vector{C}}, R::LieSubspace{C}
) where {C<:FieldElem}
    F = coefficient_ring(L)
    n = dim(L)
    m = length(lifts)
    f = fill(zero(F), m, m, m)
    m == 0 && return f

    BX = _matrix_from_coord_cols(F, n, lifts)
    BR = basis_matrix(R)
    P = ncols(BR) == 0 ? BX : hcat(BX, BR)
    size(P) == (n, n) ||
        throw(ArgumentError("lifts ⊕ radical must be a direct-sum basis of L"))
    is_unit(det(P)) ||
        throw(ArgumentError("lifts ⊕ radical must form a basis of L"))

    for a in 1:m, b in 1:m
        br = lie_bracket(L, lifts[a], lifts[b])
        rhs = matrix(F, n, 1, br)
        ok, sol = AbstractAlgebra.Solve.can_solve_with_solution(P, rhs; side = :right)
        ok || throw(ErrorException("bracket not in lifts ⊕ radical"))
        for c in 1:m
            f[a, b, c] = sol[c, 1]
        end
    end
    return f
end

"""
    quotient_algebra(L::LieAlgebra; radical = radical(L), complement = nothing)
        -> QuotientAlgebra

Build the quotient Lie algebra `L / rad(L)`.

1. Take a basis of the radical (`dim = k`).
2. Choose a complementary basis `{x₁,…,x_{n-k}}` (default: [`complement`](@ref)).
3. Induce structure constants of the quotient from brackets modulo the radical.

The returned `complement` is only a **vector-space** section of `L → L/rad(L)`.
It need not be a subalgebra; use [`levi_decomposition`](@ref) to correct it to a
Levi factor.
"""
function quotient_algebra(
    L::LieAlgebra{C};
    radical::Union{Nothing,LieSubspace{C}} = nothing,
    complement::Union{Nothing,LieSubspace{C}} = nothing,
) where {C<:FieldElem}
    # Default via cartan_orthogonal to avoid kwarg name shadowing `radical(L)`.
    R = radical === nothing ? cartan_orthogonal(L) : radical
    parent(R) === L || throw(ArgumentError("radical must belong to L"))
    X = complement === nothing ? ParametricLie.complement(L, R) : complement
    parent(X) === L || throw(ArgumentError("complement must belong to L"))
    dim(X) + dim(R) == dim(L) ||
        throw(ArgumentError("complement ⊕ radical must equal L as vector spaces"))

    lifts = [copy(e.coords) for e in basis_elems(X)]
    f = _structure_constants_mod(L, lifts, R)
    Qalg = LieAlgebra(coefficient_ring(L), f)
    return QuotientAlgebra{C}(L, R, X, Qalg)
end
