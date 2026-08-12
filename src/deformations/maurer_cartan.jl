# Truncated Maurer–Cartan integration for adjoint deformations.
#
# Convention (matched to our CE / NR signs): for φ ∈ C²,
#   dφ = −[μ, φ]_NR,
# and [μ+φ, μ+φ] = 0 expands to
#   dφ = ½[φ, φ]_NR.
# With φ = ∑_{k≥1} ε^k φ_k this is, order by order,
#   d φ_k = ½ ∑_{i+j=k} [φ_i, φ_j]_NR   (k ≥ 1),
# with the k=1 case reducing to d φ_1 = 0.
#
# Computation is truncated at `order`; results are cached so φ_k (C²) and the
# order-k RHS / obstruction (C³) can be read without recomputing.

"""
    FormalDeformation{C}

Truncated formal deformation of a Lie bracket along a seed cocycle `φ₁ ∈ Z²(g,g)`.

Caches successful terms `φ_k ∈ C²` and the Maurer–Cartan right-hand sides
`ψ_k ∈ C³`. Extending the truncation reuses the cache. The
[`MCCertificate`](@ref) is a thin view of the same caches.
"""
mutable struct FormalDeformation{C<:FieldElem}
    algebra::LieAlgebra{C}
    complex::CEComplex{C}                 # adjoint CE complex (owns d / Z / B / H caches)
    φ1::Vector{C}
    order::Int                            # requested truncation
    filled::Int                           # largest k with φ_k successfully obtained
    phi::Dict{Int, Vector{C}}             # k => φ_k ∈ C²
    rhs::Dict{Int, Vector{C}}             # k => ψ_k ∈ C³  (k ≥ 2)
    stalled_at::Union{Nothing, Int}
    obstruction::Union{Nothing, Vector{C}}  # ψ_k at failure (∈ Z³ \\ B³)
end

Base.parent(D::FormalDeformation) = D.algebra
coefficient_ring(D::FormalDeformation) = coefficient_ring(D.algebra)
ce_complex(D::FormalDeformation) = D.complex
seed(D::FormalDeformation) = D.φ1
max_order(D::FormalDeformation) = D.order
filled_order(D::FormalDeformation) = D.filled
is_integrable(D::FormalDeformation) = D.stalled_at === nothing && D.filled >= D.order
stalled_at(D::FormalDeformation) = D.stalled_at

function Base.show(io::IO, D::FormalDeformation)
    st = D.stalled_at === nothing ? "ok" : "stalled@$(D.stalled_at)"
    print(io, "FormalDeformation(order=$(D.order), filled=$(D.filled), $st)")
end

"""
    MCCertificate

Snapshot view of a [`FormalDeformation`](@ref): integrability up to the
requested order, with accessors that read the deformation caches.
"""
struct MCCertificate{C<:FieldElem}
    deformation::FormalDeformation{C}
end

deformation(c::MCCertificate) = c.deformation
is_integrable(c::MCCertificate) = is_integrable(c.deformation)
stalled_at(c::MCCertificate) = stalled_at(c.deformation)
max_order(c::MCCertificate) = max_order(c.deformation)
filled_order(c::MCCertificate) = filled_order(c.deformation)

function Base.show(io::IO, c::MCCertificate)
    D = c.deformation
    if is_integrable(D)
        print(io, "MCCertificate(integrable to order $(D.order))")
    else
        print(io, "MCCertificate(obstruction at order $(D.stalled_at))")
    end
end

# --- helpers ----------------------------------------------------------------

function _require_char_not_two(F::Field)
    c = characteristic(F)
    c == 2 && throw(ArgumentError("Maurer–Cartan integration requires char ≠ 2"))
    return inv(F(2))
end

function _as_c2_cocycle(L::LieAlgebra{C}, φ1::AbstractVector) where {C<:FieldElem}
    φ = _check_adjoint_cochain(L, φ1, 2)
    Comp = ce_complex(L, adjoint_module(L))
    D2 = ce_differential(Comp, 2)
    if ncols(D2) > 0 && nrows(D2) > 0
        v = matrix(coefficient_ring(L), length(φ), 1, φ)
        dφ = D2 * v
        all(iszero, dφ) ||
            throw(ArgumentError("seed φ₁ must be a 2-cocycle (dφ₁ = 0)"))
    elseif nrows(D2) > 0 && ncols(D2) == 0
        # empty domain
    end
    # When C³=0, d² is 0×(dim C²) and every φ₁ is a cocycle.
    return φ, Comp
end

"""
RHS `ψ_k = ½ ∑_{i=1}^{k-1} [φ_i, φ_{k-i}]_NR ∈ C³`.
"""
function _mc_rhs(L::LieAlgebra{C}, phi::Dict{Int, Vector{C}}, k::Int) where {C<:FieldElem}
    F = coefficient_ring(L)
    half = _require_char_not_two(F)
    n = dim(L)
    acc = fill(zero(F), n * binomial(n, 3))
    for i in 1:(k - 1)
        br = nr_bracket(L, phi[i], 2, phi[k - i], 2)
        for t in eachindex(acc)
            acc[t] += br[t]
        end
    end
    return C[half * acc[t] for t in eachindex(acc)]
end

function _solve_d2(Comp::CEComplex{C}, ψ::Vector{C}) where {C<:FieldElem}
    F = coefficient_ring(Comp)
    n2 = cochain_dim(Comp, 2)
    n3 = cochain_dim(Comp, 3)
    length(ψ) == n3 || throw(ArgumentError("RHS length $(length(ψ)) ≠ dim C³=$n3"))
    if n3 == 0
        return true, fill(zero(F), n2)
    end
    D2 = ce_differential(Comp, 2)
    if n2 == 0
        return all(iszero, ψ), fill(zero(F), 0)
    end
    rhs = matrix(F, n3, 1, ψ)
    ok, sol = AbstractAlgebra.Solve.can_solve_with_solution(D2, rhs; side = :right)
    ok || return false, fill(zero(F), n2)
    return true, C[sol[i, 1] for i in 1:n2]
end

# --- ensure / extend --------------------------------------------------------

"""
    _ensure_order!(D::FormalDeformation, N)

Fill caches through order `N` (or until an obstruction). Idempotent: already
computed orders are not recomputed.
"""
function _ensure_order!(D::FormalDeformation{C}, N::Int) where {C<:FieldElem}
    N >= 1 || throw(ArgumentError("order must be ≥ 1, got $N"))
    D.order = max(D.order, N)
    D.stalled_at !== nothing && return D
    L = D.algebra
    while D.filled < N
        k = D.filled + 1
        if k == 1
            # seed already validated
            D.filled = 1
            continue
        end
        ψ = _mc_rhs(L, D.phi, k)
        D.rhs[k] = ψ
        ok, φk = _solve_d2(D.complex, ψ)
        if !ok
            D.stalled_at = k
            D.obstruction = ψ
            return D
        end
        D.phi[k] = φk
        D.filled = k
    end
    return D
end

"""
    extend!(D::FormalDeformation, N) -> FormalDeformation

Raise the truncation to order `N`, reusing cached terms. No-op if `N ≤ order`
and already filled (or previously stalled).
"""
function extend!(D::FormalDeformation, N::Integer)
    _ensure_order!(D, Int(N))
    return D
end

# --- constructors -----------------------------------------------------------

"""
    formal_deformation(L, φ1; order=2) -> FormalDeformation

Start a truncated Maurer–Cartan integration along the adjoint 2-cocycle `φ1`,
computing through `order` (default **2**, the minimal nontrivial slice).

Requires characteristic ≠ 2. Results are cached on the returned object.
"""
function formal_deformation(
    L::LieAlgebra{C}, φ1::AbstractVector; order::Integer = 2
) where {C<:FieldElem}
    N = Int(order)
    N >= 1 || throw(ArgumentError("order must be ≥ 1, got $N"))
    _require_char_not_two(coefficient_ring(L))
    φ, Comp = _as_c2_cocycle(L, φ1)
    D = FormalDeformation{C}(
        L,
        Comp,
        φ,
        N,
        0,
        Dict{Int, Vector{C}}(1 => φ),
        Dict{Int, Vector{C}}(),
        nothing,
        nothing,
    )
    D.filled = 1
    _ensure_order!(D, N)
    return D
end

"""
    formal_deformation(L; order=2) -> Vector{FormalDeformation}

One deformation problem per column of `H²(L, L)` (adjoint).
"""
function formal_deformation(L::LieAlgebra{C}; order::Integer = 2) where {C<:FieldElem}
    H2 = cohomology(L, adjoint_module(L), 2)
    B = basis_matrix(H2)
    n2 = size(B, 1)
    out = FormalDeformation{C}[]
    for j in 1:ncols(B)
        φ1 = [B[i, j] for i in 1:n2]
        push!(out, formal_deformation(L, φ1; order = order))
    end
    return out
end

# --- cache accessors --------------------------------------------------------

"""
    deformation_term(D, k) -> Vector

Cached `φ_k ∈ C²`. Extends the truncation if needed and still integrable.
Throws if the integration stalled before `k`.
"""
function deformation_term(D::FormalDeformation{C}, k::Integer) where {C<:FieldElem}
    kk = Int(k)
    kk >= 1 || throw(ArgumentError("k must be ≥ 1"))
    if D.stalled_at !== nothing && kk >= D.stalled_at
        throw(ArgumentError("deformation stalled at order $(D.stalled_at); no φ_$kk"))
    end
    if kk > D.filled
        _ensure_order!(D, kk)
        D.stalled_at !== nothing && kk >= D.stalled_at &&
            throw(ArgumentError("deformation stalled at order $(D.stalled_at); no φ_$kk"))
    end
    return D.phi[kk]
end

deformation_term(c::MCCertificate, k::Integer) = deformation_term(c.deformation, k)

"""
    mc_rhs(D, k) -> Vector

Cached Maurer–Cartan right-hand side `ψ_k ∈ C³` for `k ≥ 2`
(`ψ_k = ½ ∑_{i+j=k}[φ_i,φ_j]`). Available for successful steps and for the
failing order (same vector as the obstruction cochain).
"""
function mc_rhs(D::FormalDeformation{C}, k::Integer) where {C<:FieldElem}
    kk = Int(k)
    kk >= 2 || throw(ArgumentError("mc_rhs requires k ≥ 2, got $kk"))
    if haskey(D.rhs, kk)
        return D.rhs[kk]
    end
    if D.stalled_at !== nothing && kk > D.stalled_at
        throw(ArgumentError("no RHS beyond stalled order $(D.stalled_at)"))
    end
    # Need φ_1..φ_{k-1}
    for j in 1:(kk - 1)
        deformation_term(D, j)
    end
    if !haskey(D.rhs, kk)
        # compute without solving (e.g. user asks RHS at filled+1 after stall path)
        ψ = _mc_rhs(D.algebra, D.phi, kk)
        D.rhs[kk] = ψ
    end
    return D.rhs[kk]
end

mc_rhs(c::MCCertificate, k::Integer) = mc_rhs(c.deformation, k)

"""
    obstruction_cochain(D) -> Union{Nothing, Vector}

The `C³` cocycle `ψ_k` at the first failing order, or `nothing` if integrable
through the requested order.
"""
obstruction_cochain(D::FormalDeformation) = D.obstruction
obstruction_cochain(c::MCCertificate) = obstruction_cochain(c.deformation)

"""
    mc_certificate(D) -> MCCertificate

Certificate view of the (cached) truncated integration.
"""
mc_certificate(D::FormalDeformation) = MCCertificate(D)

"""
    is_rigid(L) -> Bool

Infinitesimal rigidity: `dim H²(L, L) == 0` (adjoint).
"""
is_rigid(L::LieAlgebra) = dim(cohomology(L, adjoint_module(L), 2)) == 0
