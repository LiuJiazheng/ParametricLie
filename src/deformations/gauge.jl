# Gauge equivalence of truncated formal deformations (partial gauge theory).
#
# Ordinary comparison: two FormalDeformation objects are equivalent through
# order N iff there exist α_k ∈ C¹ (k=1..N) such that, acting on D₁'s terms,
#   φ′_k = φ_k + d α_k + ∑_{i=1}^{k-1} [φ_i, α_{k-i}]_NR.
# Order 1 reduces to φ′₁ − φ₁ ∈ B² = im d¹ (almost free with CE cache).
#
# Canonical representatives: at each order, kill the B²-component so that
# φ_k lies in a fixed complement W of B² in C² (greedy on the standard basis,
# same style as CE quotients). Not a classification pipeline — just a
# deterministic section of each gauge orbit.

"""
    equivalent(D1, D2; order=nothing) -> Bool

Whether two truncated formal deformations of the **same** Lie algebra are
gauge-equivalent through `order`.

- Default `order` is `min(filled_order(D1), filled_order(D2))`.
- Both sides are `extend!`ed to `order` when needed; if either stalls before
  that order, returns `false`.
- Uses cached `φ_k`; order 1 is a single membership test in `B²`.
"""
function equivalent(
    D1::FormalDeformation{C}, D2::FormalDeformation{C}; order::Union{Nothing, Integer} = nothing
) where {C<:FieldElem}
    D1.algebra === D2.algebra ||
        throw(ArgumentError("equivalent requires deformations of the same LieAlgebra instance"))
    N = order === nothing ? min(filled_order(D1), filled_order(D2)) : Int(order)
    N >= 1 || throw(ArgumentError("order must be ≥ 1, got $N"))

    extend!(D1, N)
    extend!(D2, N)
    (is_integrable(D1) && filled_order(D1) >= N) || return false
    (is_integrable(D2) && filled_order(D2) >= N) || return false

    L = D1.algebra
    Comp = D1.complex
    n1 = cochain_dim(Comp, 1)
    n2 = cochain_dim(Comp, 2)
    F = coefficient_ring(L)
    Dd = ce_differential(Comp, 1)   # d¹: C¹ → C²

    alphas = Vector{C}[]
    for k in 1:N
        φ = deformation_term(D1, k)
        φ′ = deformation_term(D2, k)
        rhs = C[φ′[t] - φ[t] for t in 1:n2]
        for i in 1:(k - 1)
            br = nr_bracket(L, deformation_term(D1, i), 2, alphas[k - i], 1)
            for t in 1:n2
                rhs[t] -= br[t]
            end
        end
        ok, αk = _solve_d1(Dd, rhs, n1, n2, F)
        ok || return false
        push!(alphas, αk)
    end
    return true
end

equivalent(c1::MCCertificate, c2::MCCertificate; kwargs...) =
    equivalent(c1.deformation, c2.deformation; kwargs...)

function _solve_d1(Dd::MatElem, rhs::Vector{C}, n1::Int, n2::Int, F) where {C<:FieldElem}
    length(rhs) == n2 || throw(ArgumentError("RHS length $(length(rhs)) ≠ dim C²=$n2"))
    if n2 == 0
        return true, fill(zero(F), n1)
    end
    if n1 == 0
        return all(iszero, rhs), fill(zero(F), 0)
    end
    if nrows(Dd) == 0
        return all(iszero, rhs), fill(zero(F), n1)
    end
    b = matrix(F, n2, 1, rhs)
    ok, sol = AbstractAlgebra.Solve.can_solve_with_solution(Dd, b; side = :right)
    ok || return false, fill(zero(F), n1)
    return true, C[sol[i, 1] for i in 1:n1]
end

"""
    equivalent_with_gauge(D1, D2; order=nothing) -> NamedTuple

Like [`equivalent`](@ref), but on success also returns the gauge cochains
`alphas::Vector` with `alphas[k] ∈ C¹`. On failure `ok=false` and `alphas=[]`.
"""
function equivalent_with_gauge(
    D1::FormalDeformation{C}, D2::FormalDeformation{C}; order::Union{Nothing, Integer} = nothing
) where {C<:FieldElem}
    D1.algebra === D2.algebra ||
        throw(ArgumentError("equivalent requires deformations of the same LieAlgebra instance"))
    N = order === nothing ? min(filled_order(D1), filled_order(D2)) : Int(order)
    N >= 1 || throw(ArgumentError("order must be ≥ 1, got $N"))

    extend!(D1, N)
    extend!(D2, N)
    if !(is_integrable(D1) && filled_order(D1) >= N && is_integrable(D2) && filled_order(D2) >= N)
        return (ok = false, order = N, alphas = Vector{C}[])
    end

    L = D1.algebra
    Comp = D1.complex
    n1 = cochain_dim(Comp, 1)
    n2 = cochain_dim(Comp, 2)
    F = coefficient_ring(L)
    Dd = ce_differential(Comp, 1)

    alphas = Vector{C}[]
    for k in 1:N
        φ = deformation_term(D1, k)
        φ′ = deformation_term(D2, k)
        rhs = C[φ′[t] - φ[t] for t in 1:n2]
        for i in 1:(k - 1)
            br = nr_bracket(L, deformation_term(D1, i), 2, alphas[k - i], 1)
            for t in 1:n2
                rhs[t] -= br[t]
            end
        end
        ok, αk = _solve_d1(Dd, rhs, n1, n2, F)
        if !ok
            return (ok = false, order = N, alphas = Vector{C}[])
        end
        push!(alphas, αk)
    end
    return (ok = true, order = N, alphas = alphas)
end

# --- gauge normal form ------------------------------------------------------

"""
Fixed splitting `C² = B² ⊕ W` with `B = coboundaries(·,2)` and `W` the greedy
complement of `B` inside the standard basis of `C²` (same style as CE `H^k`).
"""
function _c2_splitting(Comp::CEComplex{C}) where {C<:FieldElem}
    F = coefficient_ring(Comp)
    n2 = cochain_dim(Comp, 2)
    B = coboundaries(Comp, 2)
    W = _cochain_quotient_basis(identity_matrix(F, n2), B)
    return B, W
end

"""
Unique decomposition `u = b + w` with `b ∈ im B`, `w ∈ im W` for the splitting
`C² = B² ⊕ W`.
"""
function _split_along_b2(u::Vector{C}, B::MatElem, W::MatElem) where {C<:FieldElem}
    F = AbstractAlgebra.base_ring(B)
    n2 = length(u)
    nb = ncols(B)
    nw = ncols(W)
    if n2 == 0
        return C[], C[]
    end
    if nb == 0
        return fill(zero(F), n2), C[u[t] for t in 1:n2]
    end
    if nw == 0
        return C[u[t] for t in 1:n2], fill(zero(F), n2)
    end
    M = hcat(B, W)
    rhs = matrix(F, n2, 1, u)
    ok, sol = AbstractAlgebra.Solve.can_solve_with_solution(M, rhs; side = :right)
    ok || throw(ArgumentError("C² ≠ B² ⊕ W splitting failed (inconsistent complement)"))
    b = fill(zero(F), n2)
    w = fill(zero(F), n2)
    for j in 1:nb
        cj = sol[j, 1]
        iszero(cj) && continue
        for i in 1:n2
            b[i] += B[i, j] * cj
        end
    end
    for j in 1:nw
        cj = sol[nb + j, 1]
        iszero(cj) && continue
        for i in 1:n2
            w[i] += W[i, j] * cj
        end
    end
    return C[b[t] for t in 1:n2], C[w[t] for t in 1:n2]
end

"""
Compute gauge-canonical terms `φ̃_k ∈ W` and the realizing `α_k`, using the
original (pre-gauge) `φ_k` in the NR corrections.
"""
function _gauge_canonical_terms(D::FormalDeformation{C}, N::Int) where {C<:FieldElem}
    L = D.algebra
    Comp = D.complex
    n1 = cochain_dim(Comp, 1)
    n2 = cochain_dim(Comp, 2)
    F = coefficient_ring(L)
    Dd = ce_differential(Comp, 1)
    Bmat, Wmat = _c2_splitting(Comp)

    phi_orig = Vector{C}[deformation_term(D, k) for k in 1:N]
    alphas = Vector{C}[]
    phi_can = Dict{Int, Vector{C}}()
    for k in 1:N
        u = C[phi_orig[k][t] for t in 1:n2]
        for i in 1:(k - 1)
            br = nr_bracket(L, phi_orig[i], 2, alphas[k - i], 1)
            for t in 1:n2
                u[t] += br[t]
            end
        end
        b, w = _split_along_b2(u, Bmat, Wmat)
        # d α_k = −b  ⇒  u + dα_k = w ∈ W
        negb = C[-b[t] for t in 1:n2]
        ok, αk = _solve_d1(Dd, negb, n1, n2, F)
        ok || throw(ErrorException("gauge normal form: d¹ α_$k = −b failed (b ∉ B²)"))
        push!(alphas, αk)
        phi_can[k] = w
    end
    return phi_can, alphas
end

function _formal_deformation_from_phi(
    L::LieAlgebra{C}, Comp::CEComplex{C}, phi::Dict{Int, Vector{C}}, N::Int
) where {C<:FieldElem}
    φ1 = phi[1]
    rhs = Dict{Int, Vector{C}}()
    for k in 2:N
        rhs[k] = _mc_rhs(L, phi, k)
    end
    return FormalDeformation{C}(L, Comp, φ1, N, N, phi, rhs, nothing, nothing)
end

"""
    gauge_normal_form(D; order=nothing) -> FormalDeformation

Gauge-canonical representative of `D` through `order` (default:
`filled_order(D)`): each `φ_k` lies in a fixed complement `W` of `B²` in `C²`.

Returns a **new** [`FormalDeformation`](@ref) with caches filled from the
normalized terms (same adjoint CE complex). Requires integrability through
`order`. The realizing gauge is available via [`gauge_normal_form_with_gauge`](@ref).
"""
function gauge_normal_form(
    D::FormalDeformation{C}; order::Union{Nothing, Integer} = nothing
) where {C<:FieldElem}
    N = order === nothing ? filled_order(D) : Int(order)
    N >= 1 || throw(ArgumentError("order must be ≥ 1, got $N"))
    extend!(D, N)
    (is_integrable(D) && filled_order(D) >= N) ||
        throw(ArgumentError("cannot normalize: deformation not integrable through order $N"))
    phi_can, _ = _gauge_canonical_terms(D, N)
    return _formal_deformation_from_phi(D.algebra, D.complex, phi_can, N)
end

"""
    gauge_normal_form_with_gauge(D; order=nothing) -> NamedTuple

Like [`gauge_normal_form`](@ref), also returning `alphas` with `D`'s terms
mapped to the canonical ones by the usual truncated gauge formula.
"""
function gauge_normal_form_with_gauge(
    D::FormalDeformation{C}; order::Union{Nothing, Integer} = nothing
) where {C<:FieldElem}
    N = order === nothing ? filled_order(D) : Int(order)
    N >= 1 || throw(ArgumentError("order must be ≥ 1, got $N"))
    extend!(D, N)
    (is_integrable(D) && filled_order(D) >= N) ||
        throw(ArgumentError("cannot normalize: deformation not integrable through order $N"))
    phi_can, alphas = _gauge_canonical_terms(D, N)
    Dn = _formal_deformation_from_phi(D.algebra, D.complex, phi_can, N)
    return (deformation = Dn, order = N, alphas = alphas)
end

"""
    gauge_normal_form!(D; order=nothing) -> FormalDeformation

In-place version of [`gauge_normal_form`](@ref): replaces cached `φ_k` / `ψ_k`
by the gauge-canonical ones through `order`.
"""
function gauge_normal_form!(
    D::FormalDeformation{C}; order::Union{Nothing, Integer} = nothing
) where {C<:FieldElem}
    N = order === nothing ? filled_order(D) : Int(order)
    N >= 1 || throw(ArgumentError("order must be ≥ 1, got $N"))
    extend!(D, N)
    (is_integrable(D) && filled_order(D) >= N) ||
        throw(ArgumentError("cannot normalize: deformation not integrable through order $N"))
    phi_can, _ = _gauge_canonical_terms(D, N)
    empty!(D.phi)
    empty!(D.rhs)
    for (k, v) in phi_can
        D.phi[k] = v
    end
    D.φ1 = phi_can[1]
    for k in 2:N
        D.rhs[k] = _mc_rhs(D.algebra, D.phi, k)
    end
    D.order = max(D.order, N)
    D.filled = N
    D.stalled_at = nothing
    D.obstruction = nothing
    return D
end

gauge_normal_form(c::MCCertificate; kwargs...) = gauge_normal_form(c.deformation; kwargs...)
gauge_normal_form!(c::MCCertificate; kwargs...) = gauge_normal_form!(c.deformation; kwargs...)
