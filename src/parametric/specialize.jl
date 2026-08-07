# v0.2.2 — Specialization and generic (fraction-field) views of a family.
#
# Generic analysis runs ordinary field algorithms over Frac(k[a]).
# Results hold on a nonempty Zariski-open set, but this milestone does **not**
# discover or certify the exceptional locus (see v0.2.3–0.2.4).

"""
    _coefficient_base_field(R::Ring) -> Field

Base field `k` of a polynomial ring `k[a,…]` or its fraction field.
"""
function _coefficient_base_field(R::Ring)
    if R isa AbstractAlgebra.FracField
        return _coefficient_base_field(AbstractAlgebra.base_ring(R))
    end
    try
        B = AbstractAlgebra.base_ring(R)
        B isa Field && return B
    catch
    end
    R isa Field && return R
    throw(ArgumentError("cannot determine base field of coefficient ring $R"))
end

"""
    _is_polynomial_coeff_ring(R) -> Bool

Whether `R` is a polynomial ring (not a field / fraction field).
"""
function _is_polynomial_coeff_ring(R::Ring)
    R isa AbstractAlgebra.FracField && return false
    R isa Field && return false
    try
        AbstractAlgebra.ngens(R)
        return true
    catch
        return false
    end
end

"""
    generic_algebra(L::LieAlgebra) -> LieAlgebra

View of `L` suitable for **generic** field algorithms:

- already over a field (incl. `Frac(R)`) → return `L`
- over a polynomial ring `R = k[a,…]` → coerce structure constants into `Frac(R)`
"""
function generic_algebra(L::LieAlgebra{C}) where {C<:RingElem}
    R = coefficient_ring(L)
    R isa Field && return L
    _is_polynomial_coeff_ring(R) ||
        throw(ArgumentError("generic_algebra: unsupported coefficient ring $R"))

    K = AbstractAlgebra.fraction_field(R)
    n = dim(L)
    CK = elem_type(K)
    a = structure_constants(L)
    cK = Array{CK,3}(undef, n, n, n)
    for i in 1:n, j in 1:n, k in 1:n
        cK[i, j, k] = K(a[i, j, k])
    end
    return LieAlgebra(K, cK)
end

"""
Parse specialization assignments into a value vector aligned with `gens`.
"""
function _assignment_values(gens::Vector, assignments, k::Field)
    m = length(gens)
    m == 0 && return elem_type(k)[]

    if assignments isa AbstractVector
        length(assignments) == m ||
            throw(ArgumentError("expected $m specialization values (gens order), got $(length(assignments))"))
        return elem_type(k)[k(assignments[i]) for i in 1:m]
    end

    # Dict / NamedTuple — match by indeterminate, Symbol, or name string
    vals = Vector{elem_type(k)}(undef, m)
    for i in 1:m
        g = gens[i]
        vals[i] = k(_lookup_assignment(assignments, g))
    end
    return vals
end

function _lookup_assignment(assignments, g)
    if assignments isa NamedTuple
        sg = Symbol(string(g))
        haskey(assignments, sg) && return assignments[sg]
        throw(ArgumentError("missing specialization for parameter $g"))
    end

    # AbstractDict
    if haskey(assignments, g)
        return assignments[g]
    end
    sg = Symbol(string(g))
    if haskey(assignments, sg)
        return assignments[sg]
    end
    strg = string(g)
    for (key, val) in pairs(assignments)
        string(key) == strg && return val
    end
    throw(ArgumentError("missing specialization for parameter $g"))
end

"""
Specialize one coefficient at `vals` (images of the parameter gens).
Throws if a denominator vanishes (pole).
"""
function _specialize_coeff(x, vals, k::Field)
    P = parent(x)
    if P isa AbstractAlgebra.FracField
        num = AbstractAlgebra.evaluate(numerator(x), vals)
        den = AbstractAlgebra.evaluate(denominator(x), vals)
        iszero(den) && throw(ArgumentError(
            "specialization undefined: denominator vanishes (pole in structure constants)",
        ))
        # num, den live in k
        return k(num) * inv(k(den))
    end
    # polynomial (or already in k)
    if parent(x) === k || x isa typeof(one(k))
        return k(x)
    end
    return k(AbstractAlgebra.evaluate(x, vals))
end

"""
    specialize(L::LieAlgebra, assignments) -> LieAlgebra

Exact specialization of a parametric family to a concrete fiber over the base
field `k`.

`assignments` may be:

- a `Dict` / iterable of pairs: parameter → value (`a => 0`, `:a => 0`, …)
- a `NamedTuple`: `(a = 0,)`
- a `Vector` of values in [`parameters`](@ref)`(L)` order

Fails if any structure-constant denominator vanishes at the point.

Does **not** search for exceptional parameters — the caller chooses the fiber.
"""
function specialize(L::LieAlgebra{C}, assignments) where {C<:RingElem}
    gens = collect(parameters(L))
    isempty(gens) && throw(ArgumentError(
        "specialize: algebra has no free parameters (already over a parameter-free ring)",
    ))

    R = coefficient_ring(L)
    k = _coefficient_base_field(R)
    vals = _assignment_values(gens, assignments, k)

    n = dim(L)
    a = structure_constants(L)
    Ck = elem_type(k)
    c = Array{Ck,3}(undef, n, n, n)
    for i in 1:n, j in 1:n, t in 1:n
        c[i, j, t] = _specialize_coeff(a[i, j, t], vals, k)
    end
    return LieAlgebra(k, c)
end

"""
    analyze_generic(L::LieAlgebra) -> LieAlgebraReport

Run [`analyze`](@ref) on the **generic** view of `L` (fraction field of the
parameter ring when needed).

The report's invariants are valid on a nonempty Zariski-open subset of parameter
space (where denominators and all pivots used by the v0.1 algorithms are
nonzero). This does **not** certify the exceptional locus — see v0.2.3+.
"""
function analyze_generic(L::LieAlgebra{C}) where {C<:RingElem}
    return analyze(generic_algebra(L))
end
