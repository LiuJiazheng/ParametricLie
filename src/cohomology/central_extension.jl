# Central extensions from trivial-module 2-cocycles.
#
# A class [ω] ∈ H²(g, F) classifies (equivalence classes of) central extensions
#   0 → F → ĝ → g → 0,
# with bracket [x,y]_ĝ = [x,y]_g + ω(x,y)·z  and  [ĝ, z] = 0.
# This file is the H² application layer (not a separate “extensions” module).

"""
    central_extension(L::LieAlgebra, ω) -> LieAlgebra

Central extension `ĝ = L ⊕ F·z` by a 2-cochain `ω ∈ C²(L, F)`.

`ω` is a coordinate vector of length `binom(n,2)` in the graded-lex layout of
[`cochain_dim`](@ref) (same as columns of `basis_matrix(cohomology(L, 2))`).

Requires `dω = 0` (otherwise the new bracket fails Jacobi). The zero cocycle
gives the split extension `L ⊕ F` (abelian direct summand). A nontrivial class
gives a nonsplit extension (e.g. Heisenberg from abelian `ℚ²`).
"""
function central_extension(L::LieAlgebra{C}, ω::AbstractVector) where {C<:FieldElem}
    n = dim(L)
    F = coefficient_ring(L)
    nk = binomial(n, 2)
    length(ω) == nk ||
        throw(ArgumentError("2-cochain must have length binom($n,2)=$nk, got $(length(ω))"))
    ωv = C[F(ω[i]) for i in 1:nk]

    # Jacobi of the extension ⇔ d²ω = 0 over the trivial module
    dω = _apply_ce_d(L, trivial_module(L), 2, ωv)
    all(iszero, dω) ||
        throw(ArgumentError("ω is not a 2-cocycle (dω ≠ 0); central extension would fail Jacobi"))

    a = structure_constants(L)
    c = fill(zero(F), n + 1, n + 1, n + 1)
    # copy original brackets into the first n basis vectors
    for i in 1:n, j in 1:n, k in 1:n
        c[i, j, k] = a[i, j, k]
    end
    # central brackets: [e_i, e_j] += ω(e_i,e_j) e_{n+1}
    Is = multi_indices(n, 2)
    for (r, I) in enumerate(Is)
        w = ωv[r]
        iszero(w) && continue
        i, j = I[1], I[2]
        c[i, j, n + 1] += w
        c[j, i, n + 1] -= w
    end
    # e_{n+1} is central: all brackets involving it stay zero
    return LieAlgebra(F, c)
end

"""
    is_trivial_cocycle(L, ω) -> Bool

Whether the 2-cochain `ω` is a coboundary (`[ω] = 0` in `H²(L, F)`),
i.e. the corresponding central extension is split.
"""
function is_trivial_cocycle(L::LieAlgebra{C}, ω::AbstractVector) where {C<:FieldElem}
    n = dim(L)
    F = coefficient_ring(L)
    nk = binomial(n, 2)
    length(ω) == nk ||
        throw(ArgumentError("2-cochain must have length binom($n,2)=$nk, got $(length(ω))"))
    ωv = C[F(ω[i]) for i in 1:nk]
    B = coboundaries(ce_complex(L, trivial_module(L)), 2)
    return _in_column_span(B, ωv)
end
