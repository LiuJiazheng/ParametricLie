# v0.2.5 — Parameter stratification and jump analysis.
#
# Interprets a completed conditional analysis tree as strata + jump reports.
# Does not introduce a new linear-algebra engine.

"""
    Stratum

One certified region of parameter space: assumptions, invariant signature,
optional witnesses, and the source leaf certificate trail.
"""
struct Stratum{P}
    sigma::AssumptionSet{P}
    signature::Dict{Symbol,Any}
    witnesses::Dict{Symbol,Any}
    trail::Vector{PivotCertificate{P}}
    complete::Bool
end

function Base.show(io::IO, S::Stratum)
    parts = String[]
    for k in (:center_dim, :derived_dim, :killing_rank, :is_solvable, :is_nilpotent)
        haskey(S.signature, k) && push!(parts, "$k=$(S.signature[k])")
    end
    print(io, "Stratum($(S.sigma); ", join(parts, ", "), ")")
end

"""
    JumpEntry

One invariant change relative to a reference (usually generic) stratum.
"""
struct JumpEntry
    invariant::Symbol
    generic_value
    special_value
end

function Base.show(io::IO, j::JumpEntry)
    print(io, "$(j.invariant): $(j.generic_value) → $(j.special_value)")
end

"""
    JumpReport

Confirmed jumps from a reference stratum to an exceptional stratum.
Unchanged invariants are listed separately.
"""
struct JumpReport
    reference::Stratum
    target::Stratum
    changed::Vector{JumpEntry}
    unchanged::Vector{Symbol}
end

function Base.show(io::IO, ::MIME"text/plain", J::JumpReport)
    println(io, "Jump  ", J.reference.sigma, "  →  ", J.target.sigma)
    if isempty(J.changed)
        println(io, "  (no confirmed invariant change)")
    else
        for c in J.changed
            println(io, "  ", c)
        end
    end
end

"""
    Stratification

User-facing parameter stratification produced by [`stratify`](@ref).
Certified but not necessarily geometrically minimal/canonical.
"""
struct Stratification{P}
    family::LieAlgebra
    tree::CondTree
    generic::Union{Stratum{P},Nothing}
    strata::Vector{Stratum{P}}
    jumps::Vector{JumpReport}
    invariants::Vector{Symbol}
end

function _signature_line(sig::Dict{Symbol,Any})
    parts = String[]
    for k in (:center_dim, :derived_dim, :killing_rank, :radical_dim,
              :is_solvable, :is_nilpotent)
        haskey(sig, k) && push!(parts, "$k=$(sig[k])")
    end
    return join(parts, ", ")
end

function _jump_for_target(S::Stratification, st::Stratum)
    for j in S.jumps
        (j.target === st || j.target.sigma == st.sigma) && return j
    end
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", S::Stratification)
    n = length(S.strata)
    n_jump = length(jump_table(S))
    println(io, "Stratification: $n certified region$(n == 1 ? "" : "s")",
        " (suite=$(S.invariants))")
    if S.generic === nothing
        println(io, "  (no generic region identified)")
    end
    for (i, st) in enumerate(S.strata)
        is_gen = S.generic !== nothing &&
            (st === S.generic || st.sigma == S.generic.sigma)
        if is_gen
            println(io, "  [$i] GENERIC  Σ=$(st.sigma)")
            line = _signature_line(st.signature)
            isempty(line) || println(io, "         ", line)
            continue
        end
        J = _jump_for_target(S, st)
        if J !== nothing && !isempty(J.changed)
            println(io, "  [$i] JUMP from generic  Σ=$(st.sigma)")
            for c in J.changed
                println(io, "         ", c)
            end
        else
            println(io, "  [$i] exceptional (no confirmed jump)  Σ=$(st.sigma)")
            line = _signature_line(st.signature)
            isempty(line) || println(io, "         ", line)
        end
    end
    n_jump > 0 && println(io, "  ($n_jump confirmed jump$(n_jump == 1 ? "" : "s") vs generic)")
end

function _is_confirmed_jump(S::Stratification, st::Stratum)
    any(j -> j.target === st && !isempty(j.changed), S.jumps) && return true
    # identity compare by sigma string if === fails after copy
    for j in S.jumps
        !isempty(j.changed) || continue
        j.target.sigma == st.sigma && return true
    end
    return false
end

function _stratum_from_leaf(leaf::CondLeaf{P}) where {P}
    wit = Dict{Symbol,Any}()
    for k in (:center_basis, :derived_basis, :killing_radical_basis, :radical_basis, :der_basis)
        haskey(leaf.invariants, k) && (wit[k] = leaf.invariants[k])
    end
    return Stratum{P}(
        leaf.sigma,
        invariant_signature(leaf),
        wit,
        leaf.trail,
        leaf.complete,
    )
end

"""
Pick a generic stratum: maximize number of nonvanishings, then minimize equalities.
"""
function _pick_generic(strata::Vector{Stratum{P}}) where {P}
    isempty(strata) && return nothing
    best = strata[1]
    best_score = (-length(best.sigma.nonzeros), length(best.sigma.zeros))
    for st in strata[2:end]
        score = (-length(st.sigma.nonzeros), length(st.sigma.zeros))
        if score < best_score
            best = st
            best_score = score
        end
    end
    return best
end

const _COMPARE_KEYS = (
    :center_dim, :derived_dim, :killing_rank, :killing_radical_dim,
    :radical_dim, :der_dim, :is_solvable, :is_nilpotent,
)

function _jump_report(ref::Stratum, tgt::Stratum)
    changed = JumpEntry[]
    unchanged = Symbol[]
    present = Set{Symbol}()
    for k in keys(ref.signature)
        push!(present, k)
    end
    for k in keys(tgt.signature)
        push!(present, k)
    end
    for k in _COMPARE_KEYS
        k in present || continue
        vg = get(ref.signature, k, nothing)
        vs = get(tgt.signature, k, nothing)
        (vg === nothing && vs === nothing) && continue
        if vg == vs
            push!(unchanged, k)
        else
            push!(changed, JumpEntry(k, vg, vs))
        end
    end
    return JumpReport(ref, tgt, changed, unchanged)
end

"""
    stratify(L; invariants, kwargs...) -> Stratification
    stratify(T::CondTree; invariants) -> Stratification

Convert a (possibly freshly computed) conditional analysis tree into strata and
confirmed jump reports relative to the identified generic stratum.
"""
function stratify(
    L::LieAlgebra;
    invariants = default_conditional_suite(),
    kwargs...,
)
    suite = invariants isa Symbol ? Symbol[invariants] : collect(Symbol, invariants)
    T = analyze_conditional(L; invariants = suite, kwargs...)
    return stratify(T; invariants = suite)
end

function stratify(T::CondTree{P,C}; invariants = default_conditional_suite()) where {P,C}
    suite = invariants isa Symbol ? Symbol[invariants] : collect(Symbol, invariants)
    strata = Stratum{P}[_stratum_from_leaf(leaf) for leaf in leaves(T) if leaf.complete]
    gen = _pick_generic(strata)
    jumps = JumpReport[]
    if gen !== nothing
        for st in strata
            st.sigma == gen.sigma && continue
            push!(jumps, _jump_report(gen, st))
        end
    end
    return Stratification{P}(T.algebra, T, gen, strata, jumps, suite)
end

"""
    jump_table(S::Stratification) -> Vector{JumpReport}

Confirmed jumps only (nonempty `changed`). Algorithmic exceptional conditions
with identical signatures are omitted.
"""
function jump_table(S::Stratification)
    return JumpReport[j for j in S.jumps if !isempty(j.changed)]
end
