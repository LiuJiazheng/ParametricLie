# Levi decomposition via successive linear corrections along the radical's
# derived series (abelian layers). Ideal decomposition of the semisimple factor
# is reserved — see `ideal_decomposition`.

"""
    LeviDecomposition{C}

Levi decomposition certificate for a Lie algebra `L` (characteristic 0):

    L = levi ⋉ radical   (vector-space direct sum; `levi` a semisimple subalgebra).

Fields
- `ambient` — original `L`
- `radical` — solvable radical
- `levi` — Levi subalgebra as a subspace of `L`
- `quotient` — abstract algebra `L/rad(L)` (≅ the Levi factor)
- `lifts` — corrected basis `{s_a}` of the Levi factor with
  `[s_a,s_b] = ∑_c f_{ab}^c s_c` for the quotient structure constants `f`
"""
struct LeviDecomposition{C<:FieldElem}
    ambient::LieAlgebra{C}
    radical::LieSubspace{C}
    levi::LieSubspace{C}
    quotient::LieAlgebra{C}
    lifts::Vector{LieAlgebraElem{C}}
end

Base.parent(lev::LeviDecomposition) = lev.ambient
dim(lev::LeviDecomposition) = dim(lev.levi)

function Base.show(io::IO, lev::LeviDecomposition)
    print(
        io,
        "LeviDecomposition(levi=$(dim(lev.levi)), radical=$(dim(lev.radical)), ambient=$(dim(lev.ambient)))",
    )
end

"""
    _levi_correct_layer!(L, s_coords, f, V, W)

One abelian-layer correction: find `u_a ∈ V` (unique mod `W`) such that with
`s_a ← s_a + u_a`,

    [s_a, s_b] − ∑_c f[a,b,c] s_c  ∈  W.

Working modulo `W = [V,V]` (when `V` is a derived term of the radical) kills the
quadratic `[u_a,u_b]` term, so the system is linear in the coefficients of `u_a`.
"""
function _levi_correct_layer!(
    L::LieAlgebra{C},
    s_coords::Vector{Vector{C}},
    f::Array{C,3},
    V::LieSubspace{C},
    W::LieSubspace{C},
) where {C<:FieldElem}
    F = coefficient_ring(L)
    n = dim(L)
    m = length(s_coords)
    dV = dim(V)
    dW = dim(W)
    d = dV - dW
    (d == 0 || m == 0) && return s_coords

    BU = _complement_in(L, V, W)
    ncols(BU) == d ||
        throw(ErrorException("failed to split radical layer V/W"))
    BW = basis_matrix(W)

    # Precompute ad(s_a) * BU  (n × d)
    AdBU = Vector{typeof(BU)}(undef, m)
    for a in 1:m
        AdBU[a] = ad(L, s_coords[a]) * BU
    end

    npairs = m * (m - 1) ÷ 2
    nrows = npairs * d
    ncols_sys = m * d
    A = zero_matrix(F, nrows, ncols_sys)
    rhs = zero_matrix(F, nrows, 1)

    row0 = 0
    for a in 1:m-1, b in (a+1):m
        # error e_ab = [s_a,s_b] − ∑ f s_c  (∈ V)
        e = lie_bracket(L, s_coords[a], s_coords[b])
        for c in 1:m
            fc = f[a, b, c]
            iszero(fc) && continue
            @inbounds for i in 1:n
                e[i] -= fc * s_coords[c][i]
            end
        end

        αe = _u_coords(e, BU, BW)
        for t in 1:d
            rhs[row0+t, 1] = -αe[t]
        end

        # residual += Ad_sa BU ξ_b − Ad_sb BU ξ_a − BU (∑ f_ab^c ξ_c)
        for j in 1:d
            # ξ_a : contribution − (Ad_sb BU)[:,j]
            col_a = (a - 1) * d + j
            va = Vector{C}(undef, n)
            @inbounds for i in 1:n
                va[i] = -AdBU[b][i, j]
            end
            αa = _u_coords(va, BU, BW)
            for t in 1:d
                A[row0+t, col_a] += αa[t]
            end

            # ξ_b : contribution + (Ad_sa BU)[:,j]
            col_b = (b - 1) * d + j
            vb = Vector{C}(undef, n)
            @inbounds for i in 1:n
                vb[i] = AdBU[a][i, j]
            end
            αb = _u_coords(vb, BU, BW)
            for t in 1:d
                A[row0+t, col_b] += αb[t]
            end
        end

        for c in 1:m
            fc = f[a, b, c]
            iszero(fc) && continue
            for j in 1:d
                col_c = (c - 1) * d + j
                vc = Vector{C}(undef, n)
                @inbounds for i in 1:n
                    vc[i] = -fc * BU[i, j]
                end
                αc = _u_coords(vc, BU, BW)
                for t in 1:d
                    A[row0+t, col_c] += αc[t]
                end
            end
        end

        row0 += d
    end

    ok, sol = AbstractAlgebra.Solve.can_solve_with_solution(A, rhs; side = :right)
    ok || throw(ErrorException(
        "Levi layer correction has no solution (unexpected in characteristic 0)",
    ))

    for a in 1:m, j in 1:d
        ξ = sol[(a-1)*d+j, 1]
        iszero(ξ) && continue
        @inbounds for i in 1:n
            s_coords[a][i] += ξ * BU[i, j]
        end
    end
    return s_coords
end

"""
    levi_decomposition(L::LieAlgebra; radical = radical(L), complement = nothing)
        -> LeviDecomposition

Compute a Levi subalgebra of `L` (characteristic 0).

Algorithm
1. Form a vector-space complement `{x_a}` to `rad(L)` and the quotient structure
   constants `f_{ab}^c` ([`quotient_algebra`](@ref)).
2. Walk the derived series of the radical
   `R ⊇ R⁽¹⁾ ⊇ ⋯ ⊇ 0`. At each layer the quotient is abelian, so the correction
   equation for `u_a ∈ R⁽ᵏ⁾` is **linear**.
3. After the last layer the corrected lifts `{s_a}` satisfy
   `[s_a,s_b] = ∑ f_{ab}^c s_c` and span a Levi subalgebra.

The Levi factor is unique up to conjugation by `exp(rad)`; different complements
may yield different (isomorphic) Levi subalgebras.
"""
function levi_decomposition(
    L::LieAlgebra{C};
    radical::Union{Nothing,LieSubspace{C}} = nothing,
    complement::Union{Nothing,LieSubspace{C}} = nothing,
) where {C<:FieldElem}
    Q = quotient_algebra(L; radical = radical, complement = complement)
    R = Q.radical
    X = Q.complement
    f = structure_constants(Q.algebra)

    s_coords = [copy(e.coords) for e in basis_elems(X)]
    series = radical_derived_series(L, R)
    for k in 1:(length(series)-1)
        V = series[k]
        W = series[k+1]
        dim(V) == dim(W) && continue
        _levi_correct_layer!(L, s_coords, f, V, W)
    end

    F = coefficient_ring(L)
    n = dim(L)
    m = length(s_coords)
    B = _matrix_from_coord_cols(F, n, s_coords)
    levi_space = LieSubspace{C}(L, B)
    lifts = [LieAlgebraElem{C}(L, s_coords[a]) for a in 1:m]

    return LeviDecomposition{C}(L, R, levi_space, Q.algebra, lifts)
end

"""
    levi_subalgebra(L::LieAlgebra; kwargs...) -> LieSubspace

Levi subalgebra of `L` as a subspace (see [`levi_decomposition`](@ref)).
"""
function levi_subalgebra(L::LieAlgebra{C}; kwargs...) where {C<:FieldElem}
    return levi_decomposition(L; kwargs...).levi
end
