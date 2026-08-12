# Jump explanation: attach H² / wall-cocycle / MC / isomorphism causes to jumps.
#
# Stratification already certifies *what* jumped. This module optionally explains
# *why* in the local deformation sense, and whether fibers are isomorphic.

"""
    JumpCause

Deformation / isomorphism annotation of a confirmed [`JumpReport`](@ref).

`verdict` is one of:

| Symbol | Meaning |
|--------|---------|
| `:isomorphic_degeneration` | fibers isomorphic (`P` found) — coordinate / iso-type degeneration |
| `:gauge_like` | wall class ∈ `B²` |
| `:integrable_deformation` | nontrivial wall class + MC integrable through `order` |
| `:obstructed` | nontrivial class but MC stalls |
| `:infinitesimal_moduli` | nontrivial `H²` class, MC not conclusive |
| `:invariant_only` | signature jump only (no wall / deformation data) |
| `:unknown` | incomplete |
"""
struct JumpCause{C<:FieldElem}
    verdict::Symbol
    point_generic
    point_special
    fiber_generic::LieAlgebra{C}
    fiber_special::LieAlgebra{C}
    h2_dim_generic::Int
    h2_dim_special::Int
    rigid_generic::Bool
    rigid_special::Bool
    wall_cocycle::Union{Nothing, Vector{C}}
    wall_class::Symbol
    mc::Union{Nothing, FormalDeformation{C}}
    iso::Union{Nothing, IsoCertificate{C}}
    notes::String
end

function Base.show(io::IO, ::MIME"text/plain", c::JumpCause)
    println(io, "JumpCause(verdict=$(c.verdict), wall_class=$(c.wall_class))")
    println(io, "  H² dims: generic=$(c.h2_dim_generic), special=$(c.h2_dim_special)")
    println(io, "  rigid:   generic=$(c.rigid_generic), special=$(c.rigid_special)")
    if c.iso !== nothing
        println(io, "  iso:     ", c.iso)
    end
    isempty(c.notes) || println(io, "  notes:   ", c.notes)
end

function Base.show(io::IO, c::JumpCause)
    print(io, "JumpCause($(c.verdict))")
end

# --- stratum sampling -------------------------------------------------------

"""
    sample_stratum_point(L, Σ; trials=40) -> Dict

Heuristic integer point on the base field satisfying assumption set `Σ`.
Throws if none found in `trials` attempts.
"""
function sample_stratum_point(L::LieAlgebra, Σ::AssumptionSet; trials::Int = 40)
    gens = collect(parameters(L))
    isempty(gens) && throw(ArgumentError("family has no free parameters"))
    k = _coefficient_base_field(coefficient_ring(generic_algebra(L)))
    m = length(gens)

    # Seed: set lone-variable zeros to 0; set differences x-y to equal values.
    function seed_from_zeros()
        vals = Dict{Any,Any}(g => k(1) for g in gens)
        for f in Σ.zeros
            # single generator
            for g in gens
                if f == g || f == -g
                    vals[g] = k(0)
                end
            end
            # x - y
            if AbstractAlgebra.total_degree(f) == 1
                for g in gens
                    d = AbstractAlgebra.derivative(f, g)
                    AbstractAlgebra.total_degree(d) == 0 || continue
                    c = AbstractAlgebra.constant_coefficient(d)
                    if c == 1 || c == -1
                        for h in gens
                            h === g && continue
                            dh = AbstractAlgebra.derivative(f, h)
                            AbstractAlgebra.total_degree(dh) == 0 || continue
                            ch = AbstractAlgebra.constant_coefficient(dh)
                            if c + ch == 0
                                vals[h] = vals[g]
                            end
                        end
                    end
                end
            end
        end
        return vals
    end

    function ok_point(vals)
        try
            assign = Dict(gens[i] => k(vals[gens[i]]) for i in 1:m)
            _assignment_from_point(L, Σ, assign)
            return assign
        catch
            return nothing
        end
    end

    base = seed_from_zeros()
    got = ok_point(base)
    got !== nothing && return got

    # Brute small integer grid
    candidates = (-2, -1, 0, 1, 2, 3)
    count = 0
    function rec(i, cur)
        count >= trials && return nothing
        if i > m
            count += 1
            return ok_point(cur)
        end
        g = gens[i]
        for t in candidates
            cur[g] = k(t)
            got2 = rec(i + 1, cur)
            got2 !== nothing && return got2
        end
        return nothing
    end
    got = rec(1, Dict{Any,Any}(g => k(0) for g in gens))
    got === nothing && throw(ArgumentError("could not sample a point on Σ = $Σ"))
    return got
end

# --- wall cocycle -----------------------------------------------------------

"""
Directional derivative of a polynomial / fraction at `vals` along `dir`.
`vals` / `dir` are vectors aligned with `gens` (the family's parameter gens).
"""
function _directional_deriv_at(f, gens, vals, dir, k)
    # Lift Frac(k[a]) elements to numerator/denominator in the poly ring of `gens`.
    if f isa AbstractAlgebra.FracElem || (parent(f) isa AbstractAlgebra.FracField)
        num = AbstractAlgebra.numerator(f)
        den = AbstractAlgebra.denominator(f)
        # quotient rule: (n' d - n d') / d^2
        np = _directional_deriv_at(num, gens, vals, dir, k)
        dp = _directional_deriv_at(den, gens, vals, dir, k)
        n0 = k(AbstractAlgebra.evaluate(num, vals))
        d0 = k(AbstractAlgebra.evaluate(den, vals))
        iszero(d0) && throw(ArgumentError("denominator vanishes at sample point"))
        return (np * d0 - n0 * dp) * inv(d0)^2
    end
    # Constant field element (no parameters)
    try
        pr = parent(f)
        if pr isa Field && !(pr isa AbstractAlgebra.FracField)
            return zero(k)
        end
        AbstractAlgebra.ngens(pr)
    catch
        return zero(k)
    end
    s = zero(k)
    for (i, g) in enumerate(gens)
        iszero(dir[i]) && continue
        df = AbstractAlgebra.derivative(f, g)
        s += k(dir[i]) * k(AbstractAlgebra.evaluate(df, vals))
    end
    return s
end

"""
    wall_cocycle(family, point, direction) -> Vector

C² cochain on the specialized fiber: directional derivative of the structure
constants of `family` at `point` along parameter `direction` (Dict / Vector in
`parameters(family)` order), packed as an adjoint cochain.
"""
function wall_cocycle(family::LieAlgebra, point, direction)
    gens = collect(parameters(family))
    k = _coefficient_base_field(coefficient_ring(generic_algebra(family)))
    vals = _assignment_values(gens, point, k)
    dir = _assignment_values(gens, direction, k)
    n = dim(family)
    a = structure_constants(family)
    # Build C² length n*binom(n,2)
    φ = fill(zero(k), n * binomial(n, 2))
    for (r, I) in enumerate(multi_indices(n, 2))
        i, j = I[1], I[2]
        off = (r - 1) * n
        for t in 1:n
            φ[off + t] = _directional_deriv_at(a[i, j, t], gens, vals, dir, k)
        end
    end
    return φ
end

"""
Wall normal in parameter space from new equalities on the target stratum.
Evaluated at `point_special`.
"""
function _wall_direction(family::LieAlgebra, Σ_ref::AssumptionSet, Σ_tgt::AssumptionSet, point_special)
    gens = collect(parameters(family))
    k = _coefficient_base_field(coefficient_ring(generic_algebra(family)))
    vals = _assignment_values(gens, point_special, k)
    dir = fill(zero(k), length(gens))
    for f in Σ_tgt.zeros
        status(Σ_ref, f) === PIVOT_ZERO && continue
        for (i, g) in enumerate(gens)
            df = AbstractAlgebra.derivative(f, g)
            dir[i] += k(AbstractAlgebra.evaluate(df, vals))
        end
    end
    if all(iszero, dir)
        dir[1] = one(k)
    end
    return Dict(gens[i] => dir[i] for i in eachindex(gens))
end

function _classify_c2_class(L::LieAlgebra{C}, φ::Vector{C}) where {C<:FieldElem}
    Comp = ce_complex(L, adjoint_module(L))
    n2 = cochain_dim(Comp, 2)
    length(φ) == n2 || throw(ArgumentError("wall cocycle length $(length(φ)) ≠ dim C²=$n2"))
    D2 = ce_differential(Comp, 2)
    if nrows(D2) > 0 && ncols(D2) > 0
        v = matrix(coefficient_ring(L), n2, 1, φ)
        dφ = D2 * v
        all(iszero, dφ) || return :not_cocycle
    end
    B = coboundaries(Comp, 2)
    _in_column_span(B, φ) && return :coboundary
    return :nontrivial
end

function _project_to_cocycle_rep(L::LieAlgebra{C}, φ::Vector{C}) where {C<:FieldElem}
    # If already a cocycle, return φ (caller may still want H² section).
    # If nontrivial, leave as-is for MC seed (must be cocycle).
    return φ
end

function _verdict_from(; iso, wall_class, mc, order)
    if iso !== nothing && iso.isomorphic
        return :isomorphic_degeneration
    end
    if wall_class === :coboundary
        return :gauge_like
    end
    if wall_class === :nontrivial
        if mc !== nothing && is_integrable(mc) && filled_order(mc) >= order
            return :integrable_deformation
        elseif mc !== nothing && stalled_at(mc) !== nothing
            return :obstructed
        else
            return :infinitesimal_moduli
        end
    end
    if iso !== nothing && !iso.isomorphic && iso.reason === :invariants
        return :invariant_only
    end
    return :unknown
end

"""
    explain_jump(family, jump; point_generic=nothing, point_special=nothing,
                 order=2, direction=nothing) -> JumpCause

Build a [`JumpCause`](@ref) for one [`JumpReport`](@ref): specialize fibers,
compute adjoint `H²`, wall cocycle class, optional MC along the wall class, and
[`isomorphism`](@ref) between the two fibers.
"""
function explain_jump(
    family::LieAlgebra,
    jump::JumpReport;
    point_generic = nothing,
    point_special = nothing,
    order::Integer = 2,
    direction = nothing,
)
    N = Int(order)
    pg = point_generic === nothing ? sample_stratum_point(family, jump.reference.sigma) : point_generic
    ps = point_special === nothing ? sample_stratum_point(family, jump.target.sigma) : point_special
    Lg = specialize(family, pg)
    Ls = specialize(family, ps)
    C = elem_type(coefficient_ring(Lg))

    H2g = cohomology(Lg, adjoint_module(Lg), 2)
    H2s = cohomology(Ls, adjoint_module(Ls), 2)
    rg = is_rigid(Lg)
    rs = is_rigid(Ls)

    dir = direction === nothing ? _wall_direction(family, jump.reference.sigma, jump.target.sigma, ps) : direction
    φ = wall_cocycle(family, ps, dir)
    # φ lives over base field of specialization — same as coefficient_ring(Ls)
    φC = C[φ[t] for t in eachindex(φ)]
    wclass = _classify_c2_class(Ls, φC)

    mc = nothing
    notes = ""
    if wclass === :nontrivial
        try
            mc = formal_deformation(Ls, φC; order = N)
            notes = is_integrable(mc) ? "MC integrable to order $N on special fiber" :
                "MC stalled at $(stalled_at(mc))"
        catch e
            notes = "MC failed: $e"
            wclass = :not_cocycle  # seed rejected
            mc = nothing
        end
    elseif wclass === :coboundary
        notes = "wall cocycle is a coboundary on the special fiber"
    elseif wclass === :not_cocycle
        notes = "wall directional derivative is not a 2-cocycle on the special fiber"
    end

    iso = isomorphism(Lg, Ls)
    verdict = _verdict_from(iso = iso, wall_class = wclass, mc = mc, order = N)
    if verdict === :invariant_only && wclass === :nontrivial && mc !== nothing && is_integrable(mc)
        verdict = :integrable_deformation
    end

    return JumpCause{C}(
        verdict, pg, ps, Lg, Ls,
        dim(H2g), dim(H2s), rg, rs,
        φC, wclass, mc, iso, notes,
    )
end

"""
    explain_jumps!(S::Stratification; order=2, points=nothing) -> Stratification

Annotate every confirmed jump in `S` with a [`JumpCause`](@ref).
`points` may be a `Dict` mapping stratum index or `AssumptionSet` → assignment.
"""
function explain_jumps!(
    S::Stratification;
    order::Integer = 2,
    points = nothing,
)
    for J in jump_table(S)
        pg = nothing
        ps = nothing
        if points !== nothing
            for (key, pt) in pairs(points)
                if key isa Integer
                    st = S.strata[Int(key)]
                    if st === J.reference || st.sigma == J.reference.sigma
                        pg = pt
                    end
                    if st === J.target || st.sigma == J.target.sigma
                        ps = pt
                    end
                elseif key isa AssumptionSet
                    if key == J.reference.sigma
                        pg = pt
                    end
                    if key == J.target.sigma
                        ps = pt
                    end
                end
            end
        end
        J.cause = explain_jump(S.family, J; point_generic = pg, point_special = ps, order = order)
    end
    return S
end

function explain_jumps(S::Stratification; kwargs...)
    # shallow: mutate jump reports in place (they are the same objects)
    return explain_jumps!(S; kwargs...)
end
