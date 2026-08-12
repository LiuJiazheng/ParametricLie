# Nijenhuis–Richardson bracket on adjoint cochains C^•(g, g).
#
# Same coordinates as CE: for ω ∈ C^k, length n * binom(n,k); lex multi-index
# I, then n coordinates of ω(e_I) ∈ g ≅ F^n.
#
# For φ ∈ C^p, ψ ∈ C^q (p,q ≥ 1):
#   (φ ∘ ψ)(x₁,…,x_{p+q-1})
#     = ∑_{σ ∈ Sh(q,p-1)} sgn(σ)  φ( ψ(x_{σ(1)},…,x_{σ(q)}), x_{σ(q+1)},… )
#   [φ,ψ]_NR = φ∘ψ − (−1)^{(p−1)(q−1)} ψ∘φ  ∈ C^{p+q−1}.

"""
    adjoint_bracket_cochain(L::LieAlgebra) -> Vector

The Lie bracket of `L` as an element `μ ∈ C²(L, L)` in adjoint cochain
coordinates: `μ(e_i, e_j) = [e_i, e_j]`.
"""
function adjoint_bracket_cochain(L::LieAlgebra{C}) where {C<:FieldElem}
    n = dim(L)
    F = coefficient_ring(L)
    a = structure_constants(L)
    μ = fill(zero(F), n * binomial(n, 2))
    for (r, I) in enumerate(multi_indices(n, 2))
        i, j = I[1], I[2]
        off = (r - 1) * n
        for k in 1:n
            μ[off + k] = a[i, j, k]
        end
    end
    return μ
end

"""
Evaluate `ω ∈ C^k(g,g)` on an ordered list of `k` basis indices (not necessarily
increasing). Returns the coordinate vector in `F^n`.
"""
function _cochain_eval_basis(ω::Vector{C}, n::Int, k::Int, idxs::Vector{Int}) where {C}
    length(ω) >= 1 || throw(ArgumentError("empty cochain"))
    F = parent(ω[1])
    z = fill(zero(F), n)
    length(idxs) == k || throw(ArgumentError("expected $k indices, got $(length(idxs))"))
    k == 0 && return C[ω[t] for t in 1:n]
    length(unique(idxs)) != k && return z
    for i in idxs
        (1 <= i <= n) || throw(ArgumentError("basis index $i out of range 1:$n"))
    end
    perm = sortperm(idxs)
    J = idxs[perm]
    s = _perm_sign(perm)
    blk = _cochain_block(ω, n, _lex_rank(n, J))
    return s == 1 ? C[blk[t] for t in 1:n] : C[-blk[t] for t in 1:n]
end

"""
Evaluate `φ ∈ C^p` on `(v, e_{j1}, …, e_{j_{p-1}})` with `v ∈ F^n`.
"""
function _cochain_eval_vec_first(
    φ::Vector{C}, n::Int, p::Int, v::Vector{C}, rest::Vector{Int}
) where {C}
    F = parent(first(φ))
    out = fill(zero(F), n)
    length(rest) == p - 1 ||
        throw(ArgumentError("expected $(p - 1) trailing indices, got $(length(rest))"))
    p == 0 && throw(ArgumentError("C⁰ has no vec-first evaluation in NR circle"))
    if p == 1
        # φ ∈ Hom(g,g): φ(v) = ∑_t v_t φ(e_t)
        for t in 1:n
            iszero(v[t]) && continue
            φt = _cochain_eval_basis(φ, n, 1, Int[t])
            for u in 1:n
                out[u] += v[t] * φt[u]
            end
        end
        return out
    end
    for t in 1:n
        iszero(v[t]) && continue
        args = Vector{Int}(undef, p)
        args[1] = t
        args[2:p] = rest
        φt = _cochain_eval_basis(φ, n, p, args)
        for u in 1:n
            out[u] += v[t] * φt[u]
        end
    end
    return out
end

"""
Shuffles `Sh(q, r)`: permutations `σ` of `1:(q+r)` with
`σ(1)<⋯<σ(q)` and `σ(q+1)<⋯<σ(q+r)`. Yields `(sign, σ)` where `σ` is the
array `σ[i] = σ(i)`.
"""
function _nr_shuffles(q::Int, r::Int)
    m = q + r
    out = Tuple{Int, Vector{Int}}[]
    if m == 0
        push!(out, (1, Int[]))
        return out
    end
    if q == 0
        σ = collect(1:m)
        push!(out, (1, σ))
        return out
    end
    if r == 0
        σ = collect(1:m)
        push!(out, (1, σ))
        return out
    end
    # Choose q positions among 1:m for the first block (ψ arguments)
    # Generate combinations iteratively
    c = collect(1:q)
    while true
        rest = Int[]
        j = 1
        for x in 1:m
            if j <= q && c[j] == x
                j += 1
            else
                push!(rest, x)
            end
        end
        σ = vcat(c, rest)
        # sign of σ as a permutation (σ[i] = image of i)
        push!(out, (_perm_sign(σ), σ))
        # next combination
        i = q
        while i >= 1 && c[i] == m - q + i
            i -= 1
        end
        i < 1 && break
        c[i] += 1
        for j in (i + 1):q
            c[j] = c[j - 1] + 1
        end
    end
    return out
end

function _check_adjoint_cochain(L::LieAlgebra{C}, ω::AbstractVector, k::Int) where {C<:FieldElem}
    n = dim(L)
    nk = n * binomial(n, k)
    length(ω) == nk ||
        throw(ArgumentError("C^$k(g,g) cochain must have length $nk, got $(length(ω))"))
    F = coefficient_ring(L)
    return C[F(ω[i]) for i in 1:nk]
end

"""
    nr_circle(L, φ, p, ψ, q) -> Vector

Nijenhuis–Richardson circle product `φ ∘ ψ ∈ C^{p+q-1}(g,g)`.

`φ` is a length-`n binom(n,p)` adjoint cochain, `ψ` length-`n binom(n,q)`.
Requires `p ≥ 1`, `q ≥ 1`.
"""
function nr_circle(
    L::LieAlgebra{C}, φ::AbstractVector, p::Integer, ψ::AbstractVector, q::Integer
) where {C<:FieldElem}
    pp, qq = Int(p), Int(q)
    pp >= 1 || throw(ArgumentError("nr_circle requires p ≥ 1, got $pp"))
    qq >= 1 || throw(ArgumentError("nr_circle requires q ≥ 1, got $qq"))
    n = dim(L)
    F = coefficient_ring(L)
    φv = _check_adjoint_cochain(L, φ, pp)
    ψv = _check_adjoint_cochain(L, ψ, qq)
    deg = pp + qq - 1
    out = fill(zero(F), n * binomial(n, deg))
    deg > n && return out

    shuf = _nr_shuffles(qq, pp - 1)
    for (rJ, J) in enumerate(multi_indices(n, deg))
        # args x_i = e_{J[i]}
        acc = fill(zero(F), n)
        for (sgn, σ) in shuf
            ψ_idxs = [J[σ[t]] for t in 1:qq]
            rest = [J[σ[qq + t]] for t in 1:(pp - 1)]
            v = _cochain_eval_basis(ψv, n, qq, ψ_idxs)
            all(iszero, v) && continue
            w = _cochain_eval_vec_first(φv, n, pp, v, rest)
            sF = sgn == 1 ? one(F) : -one(F)
            for u in 1:n
                acc[u] += sF * w[u]
            end
        end
        off = (rJ - 1) * n
        for u in 1:n
            out[off + u] = acc[u]
        end
    end
    return out
end

"""
    nr_bracket(L, φ, p, ψ, q) -> Vector

Nijenhuis–Richardson bracket

    [φ, ψ]_NR = φ∘ψ − (−1)^{(p−1)(q−1)} ψ∘φ  ∈ C^{p+q−1}(g,g).

Same coordinate convention as [`nr_circle`](@ref) / CE adjoint cochains.
"""
function nr_bracket(
    L::LieAlgebra{C}, φ::AbstractVector, p::Integer, ψ::AbstractVector, q::Integer
) where {C<:FieldElem}
    pp, qq = Int(p), Int(q)
    circ = nr_circle(L, φ, pp, ψ, qq)
    # graded sign for desuspended degrees (p-1), (q-1)
    if isodd((pp - 1) * (qq - 1))
        return circ .+ nr_circle(L, ψ, qq, φ, pp)
    else
        return circ .- nr_circle(L, ψ, qq, φ, pp)
    end
end

nr_circle(L::LieAlgebra, φ::AbstractVector, ψ::AbstractVector; p::Integer, q::Integer) =
    nr_circle(L, φ, p, ψ, q)

nr_bracket(L::LieAlgebra, φ::AbstractVector, ψ::AbstractVector; p::Integer, q::Integer) =
    nr_bracket(L, φ, p, ψ, q)
