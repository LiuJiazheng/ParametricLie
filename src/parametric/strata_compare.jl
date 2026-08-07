# v0.2.6 — Stratum comparison and specialization validation.

"""
    StratumComparison

Result of [`compare`](@ref) between two strata.
`same_signature` does **not** imply isomorphism.
"""
struct StratumComparison
    left::Stratum
    right::Stratum
    changed::Vector{JumpEntry}
    unchanged::Vector{Symbol}
    same_signature::Bool
end

function Base.show(io::IO, ::MIME"text/plain", C::StratumComparison)
    println(io, "StratumComparison(same_signature=$(C.same_signature))")
    println(io, "  left:  ", C.left.sigma)
    println(io, "  right: ", C.right.sigma)
    if isempty(C.changed)
        println(io, "  changed: (none)")
    else
        println(io, "  changed:")
        for c in C.changed
            println(io, "    ", c)
        end
    end
    isempty(C.unchanged) || println(io, "  unchanged: ", join(string.(C.unchanged), ", "))
end

"""
    compare(S1::Stratum, S2::Stratum) -> StratumComparison

Compare certified invariant signatures. Never claims isomorphism.
"""
function compare(S1::Stratum, S2::Stratum)
    J = _jump_report(S1, S2)
    return StratumComparison(S1, S2, J.changed, J.unchanged, isempty(J.changed))
end

function compare(S::Stratification, i::Int, j::Int)
    return compare(S.strata[i], S.strata[j])
end

# --- specialization validation -----------------------------------------------

"""
Extract a concrete assignment Dict for parameters that are forced to zero in Σ,
and check that suggested values satisfy all nonvanishings.
"""
function _assignment_from_point(L::LieAlgebra, Σ::AssumptionSet, point)
    gens = collect(parameters(L))
    isempty(gens) && throw(ArgumentError("family has no free parameters"))
    k = _coefficient_base_field(coefficient_ring(generic_algebra(L)))
    vals = _assignment_values(gens, point, k)
    # Check zeros: each equality poly evaluates to 0
    for f in Σ.zeros
        ev = AbstractAlgebra.evaluate(f, vals)
        iszero(ev) || throw(ArgumentError(
            "point does not satisfy equality $f = 0 (got $ev)",
        ))
    end
    for g in Σ.nonzeros
        ev = AbstractAlgebra.evaluate(g, vals)
        iszero(ev) && throw(ArgumentError(
            "point violates nonvanishing $g ≠ 0",
        ))
    end
    return Dict(gens[i] => vals[i] for i in eachindex(gens))
end

"""
    validate_stratum(L, stratum, point) -> NamedTuple

Specialize `L` at `point` (must satisfy `stratum.sigma`) and compare v0.1
[`analyze`](@ref) dims/flags to the stratum signature.

Returns `(ok, fiber_report, mismatches)`.
"""
function validate_stratum(L::LieAlgebra, stratum::Stratum, point)
    assign = _assignment_from_point(L, stratum.sigma, point)
    Lf = specialize(L, assign)
    r = analyze(Lf)
    mismatches = Pair{Symbol,Any}[]
    sig = stratum.signature

    function check(key::Symbol, got)
        haskey(sig, key) || return
        exp = sig[key]
        exp == got || push!(mismatches, key => (; expected = exp, got = got))
    end

    check(:center_dim, dim(center(r)))
    check(:derived_dim, dim(derived_algebra(Lf)))
    check(:killing_rank, killing_rank(r))
    check(:is_solvable, is_solvable(r))
    check(:is_nilpotent, is_nilpotent(r))
    if haskey(sig, :killing_radical_dim)
        check(:killing_radical_dim, dim(killing_radical(r)))
    end
    if haskey(sig, :radical_dim)
        check(:radical_dim, dim(radical(r)))
    end

    return (ok = isempty(mismatches), fiber_report = r, mismatches = mismatches)
end

"""
    validate_stratification(L, S::Stratification, points) -> Vector

`points` is a vector of `(stratum_index, assignment)` pairs (1-based index into
`S.strata`). Returns one validation result per pair.
"""
function validate_stratification(L::LieAlgebra, S::Stratification, points)
    return [validate_stratum(L, S.strata[i], pt) for (i, pt) in points]
end
