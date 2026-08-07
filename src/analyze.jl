# High-level analysis entry point (v0.1).
#
# `analyze(L)` builds a [`LieAlgebraReport`](@ref) that caches structural
# certificates (series, Killing, Levi, Der, …). `show` prints a compact
# summary of basis-invariant facts; detail queries reuse the cached objects.

"""
    LieAlgebraReport{C}

Structured analysis of a [`LieAlgebra`](@ref). Constructed by [`analyze`](@ref).

**Summary (`show`)** — dimensions, series profiles, solvability / semisimplicity,
Levi type (`simple` / `composite` / empty), simple-factor dimensions, `dim Der`.

**Detail queries** (reuse cached certificates; bases/matrices are basis-dependent):

| Query | Returns |
|-------|---------|
| [`jacobi`](@ref) / `report.jacobi` | [`JacobiCertificate`](@ref) |
| [`center`](@ref) | [`LieSubspace`](@ref) |
| [`derived_series`](@ref) / [`lower_central_series`](@ref) | [`LieSeries`](@ref) |
| [`killing_form`](@ref) / [`killing_radical`](@ref) | matrix / subspace |
| [`radical`](@ref) | solvable radical |
| [`levi_decomposition`](@ref) | [`LeviDecomposition`](@ref) |
| [`ideal_decomposition`](@ref) | simple ideals of the Levi quotient |
| [`derivations`](@ref) | [`Derivations`](@ref) |
"""
struct LieAlgebraReport{C<:FieldElem}
    algebra::LieAlgebra{C}
    jacobi::JacobiCertificate
    center::LieSubspace{C}
    derived::LieSeries{C}
    lower_central::LieSeries{C}
    solvable::Bool
    nilpotent::Bool
    killing::MatElem
    killing_rank::Int
    killing_radical::LieSubspace{C}
    cartan_solvable::Bool
    semisimple::Bool
    radical::LieSubspace{C}
    levi::LeviDecomposition{C}
    simple_factor_dims::Vector{Int}
    simple_ideals::Vector{LieSubspace{C}}  # ideals of `levi.quotient`
    derivations::Derivations{C}
end

Base.parent(r::LieAlgebraReport) = r.algebra
dim(r::LieAlgebraReport) = dim(r.algebra)
coefficient_ring(r::LieAlgebraReport) = coefficient_ring(r.algebra)

# --- detail accessors (cached) ----------------------------------------------

jacobi(r::LieAlgebraReport) = r.jacobi
center(r::LieAlgebraReport) = r.center
derived_series(r::LieAlgebraReport) = r.derived
lower_central_series(r::LieAlgebraReport) = r.lower_central
is_solvable(r::LieAlgebraReport) = r.solvable
is_nilpotent(r::LieAlgebraReport) = r.nilpotent
killing_form(r::LieAlgebraReport) = r.killing
killing_rank(r::LieAlgebraReport) = r.killing_rank
killing_radical(r::LieAlgebraReport) = r.killing_radical
is_cartan_solvable(r::LieAlgebraReport) = r.cartan_solvable
is_semisimple(r::LieAlgebraReport) = r.semisimple
radical(r::LieAlgebraReport) = r.radical
levi_decomposition(r::LieAlgebraReport) = r.levi
levi_subalgebra(r::LieAlgebraReport) = r.levi.levi
ideal_decomposition(r::LieAlgebraReport) = r.simple_ideals
derivations(r::LieAlgebraReport) = r.derivations
simple_factor_dims(r::LieAlgebraReport) = r.simple_factor_dims

"""
    levi_kind(r::LieAlgebraReport) -> Symbol

`:empty` if the Levi factor is zero, `:simple` if one simple ideal,
`:composite` if two or more.
"""
function levi_kind(r::LieAlgebraReport)
    ds = r.simple_factor_dims
    isempty(ds) && return :empty
    length(ds) == 1 && return :simple
    return :composite
end

function is_simple(r::LieAlgebraReport)
    return r.semisimple && levi_kind(r) === :simple && dim(r.radical) == 0
end

# --- summary display --------------------------------------------------------

function _series_profile(S::LieSeries)
    return join(string.(dim.(terms(S))), " → ")
end

function _field_label(F)
    # Prefer short names for common rings
    s = string(F)
    occursin("Rational", s) && return "QQ"
    return s
end

function _format_factor_dims(ds::Vector{Int})
    isempty(ds) && return "[]"
    return "[" * join(string.(ds), ", ") * "]"
end

function Base.show(io::IO, ::MIME"text/plain", r::LieAlgebraReport)
    L = r.algebra
    F = coefficient_ring(L)
    println(io, "LieAlgebra dim=$(dim(L)) over $(_field_label(F))")
    println(io)
    jac = r.jacobi.ok ? "OK" : "FAIL"
    println(io, "Jacobi:           ", jac)
    println(io)
    println(io, "center:           dim $(dim(r.center))")
    println(io, "derived:          ", _series_profile(r.derived))
    println(io, "lower_central:    ", _series_profile(r.lower_central))
    println(io, "solvable:         ", r.solvable)
    println(io, "nilpotent:        ", r.nilpotent)
    println(io)
    println(io, "Killing rank:     ", r.killing_rank)
    println(io, "semisimple:       ", r.semisimple)
    println(io)
    println(io, "radical:          dim $(dim(r.radical))")
    kind = levi_kind(r)
    if kind === :empty
        println(io, "Levi:             dim 0")
    else
        println(io, "Levi:             dim $(dim(r.levi)), ", string(kind))
    end
    println(io, "simple factors:   dims ", _format_factor_dims(r.simple_factor_dims))
    println(io)
    print(io, "Der:              dim $(dim(r.derivations))")
end

function Base.show(io::IO, r::LieAlgebraReport)
    print(
        io,
        "LieAlgebraReport(dim=$(dim(r)), Levi=$(levi_kind(r)), Der=$(dim(r.derivations)))",
    )
end

# --- analyze ----------------------------------------------------------------

"""
    analyze(L::LieAlgebra) -> LieAlgebraReport

Compute classical structural invariants of `L` and cache the main certificates
(Levi decomposition, derivations, series, Killing form, …).

The printed summary emphasizes basis-invariant facts (dimensions, series
profiles, whether the Levi factor is simple or composite and the dimensions of
its simple summands). Call the detail accessors on the returned report — or the
standalone APIs — for bases and matrices.
"""
function analyze(L::LieAlgebra{C}) where {C<:FieldElem}
    jac = check_jacobi(L)
    Z = center(L)

    Dser = derived_series(L)
    Cser = lower_central_series(L)
    solv = is_solvable(Dser)
    nilp = is_nilpotent(Cser)

    K = killing_form(L)
    krank = Int(rank(K))
    Krad = killing_radical(L)
    cartan_solv = is_cartan_solvable(L)
    ss = is_semisimple(L)

    lev = levi_decomposition(L)
    R = lev.radical

    # Simple ideals of the semisimple quotient (empty if Levi = 0)
    ideals = if dim(lev.quotient) == 0
        LieSubspace{C}[]
    else
        ideal_decomposition(lev.quotient)
    end
    factor_dims = sort(Int[dim(I) for I in ideals]; rev = true)

    Der = derivations(L)

    return LieAlgebraReport{C}(
        L,
        jac,
        Z,
        Dser,
        Cser,
        solv,
        nilp,
        K,
        krank,
        Krad,
        cartan_solv,
        ss,
        R,
        lev,
        factor_dims,
        ideals,
        Der,
    )
end
