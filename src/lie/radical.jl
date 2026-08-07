# Solvable radical via Cartan's criterion (characteristic 0).

"""
    radical(L::LieAlgebra) -> LieSubspace

Solvable radical of `L` (largest solvable ideal).

In characteristic 0 this equals the Killing orthogonal of the derived algebra:

    rad(L) = [L,L]^⊥ = cartan_orthogonal(L).

Returns a [`LieSubspace`](@ref) certificate (basis matrix).
"""
radical(L::LieAlgebra) = cartan_orthogonal(L)

"""
    radical_derived_series(L::LieAlgebra, R::LieSubspace = radical(L))
        -> Vector{LieSubspace}

Derived series of the radical as subspaces of `L`:

    R = R⁽⁰⁾ ⊇ R⁽¹⁾ = [R,R] ⊇ R⁽²⁾ ⊇ ⋯ ⊇ 0.

Used by the Levi lifting algorithm (each successive quotient is abelian).
"""
function radical_derived_series(
    L::LieAlgebra{C}, R::LieSubspace{C} = radical(L)
) where {C<:FieldElem}
    parent(R) === L || throw(ArgumentError("radical subspace must belong to L"))
    terms = LieSubspace{C}[R]
    # Solvable ⇒ at most dim(R) strict drops before hitting 0
    for _ in 1:(dim(R) + 1)
        prev = terms[end]
        dim(prev) == 0 && break
        nxt = commutator_span(L, prev, prev)
        push!(terms, nxt)
        dim(nxt) == dim(prev) && break  # safety (should not occur for rad)
    end
    return terms
end
