# v0.2.3 — Assumption sets Σ and conservative three-valued pivot status.
#
# Core rule: prove before split. status returns ZERO / NONZERO only when proven.

"""
    PivotStatus

Three-valued judgment of a pivot polynomial under an assumption set.
"""
@enum PivotStatus::Int8 begin
    PIVOT_ZERO = 0
    PIVOT_NONZERO = 1
    PIVOT_UNKNOWN = 2
end

const _FACTOR_DEG_BUDGET = 8
const _FACTOR_TERM_BUDGET = 24

"""
    AssumptionSet{P}

Conjunction of polynomial equalities `f = 0` and nonvanishings `g ≠ 0`.
Polynomials are stored in a normalized form (primitive, monic, square-free when factoring succeeds).
"""
struct AssumptionSet{P}
    zeros::Vector{P}
    nonzeros::Vector{P}
end

AssumptionSet{P}() where {P} = AssumptionSet{P}(P[], P[])

function Base.show(io::IO, Σ::AssumptionSet)
    z = isempty(Σ.zeros) ? "∅" : join(string.(Σ.zeros), ", ") * " = 0"
    nz = isempty(Σ.nonzeros) ? "∅" : join(string.(Σ.nonzeros), ", ") * " ≠ 0"
    print(io, "AssumptionSet($z; $nz)")
end

Base.isempty(Σ::AssumptionSet) = isempty(Σ.zeros) && isempty(Σ.nonzeros)

function Base.:(==)(A::AssumptionSet{P}, B::AssumptionSet{P}) where {P}
    return A.zeros == B.zeros && A.nonzeros == B.nonzeros
end


# --- polynomial helpers -------------------------------------------------------

"""
    _poly_ring_of(x) -> Ring

Parent polynomial ring of a polynomial or of the numerator of a fraction.
"""
function _poly_ring_of(x)
    P = parent(x)
    if P isa AbstractAlgebra.FracField
        return AbstractAlgebra.base_ring(P)
    end
    return P
end

"""
    _numerator_poly(x)

Numerator polynomial of a field/ring element for vanishing tests.
Constants in a base field map to the constant polynomial in a dummy sense:
return `nothing` with a separate constant path via `_pivot_payload`.
"""
function _numerator_poly(x)
    P = parent(x)
    if P isa AbstractAlgebra.FracField
        return numerator(x)
    end
    return x
end

function _is_field_zero(x)
    return iszero(x)
end

"""
    _leading_coeff_inv(p)

Inverse of the leading coefficient in the coefficient field, or `nothing`.
"""
function _leading_coeff_inv(p)
    iszero(p) && return nothing
    lc = leading_coefficient(p)
    iszero(lc) && return nothing
    try
        return inv(lc)
    catch
        return nothing
    end
end

"""
    normalize_poly(p) -> p′

Primitive, monic when possible, and square-free when factorization is within budget.
Returns `0` unchanged.
"""
function normalize_poly(p)
    iszero(p) && return p
    R = parent(p)
    # content → primitive
    try
        c = content(p)
        if !iszero(c) && !isone(c)
            p = divexact(p, c)
        end
    catch
    end
    # monic
    invlc = _leading_coeff_inv(p)
    if invlc !== nothing && !isone(invlc)
        p = invlc * p
    end
    # square-free part via factor (bounded)
    if _within_factor_budget(p)
        try
            p = _square_free_part(p)
            invlc = _leading_coeff_inv(p)
            if invlc !== nothing && !isone(invlc)
                p = invlc * p
            end
        catch
        end
    end
    return p
end

function _within_factor_budget(p)
    try
        AbstractAlgebra.total_degree(p) > _FACTOR_DEG_BUDGET && return false
    catch
    end
    try
        length(p) > _FACTOR_TERM_BUDGET && return false
    catch
    end
    return true
end

function _square_free_part(p)
    iszero(p) && return p
    fac = factor(p)
    r = one(parent(p))
    for (f, _e) in fac
        r *= f
    end
    return r
end

"""
    _factor_atoms(p) -> Vector

Irreducible (or atomic) factors of `p` for assumption reasoning.
Falls back to `[normalize_poly(p)]` when factoring is skipped or fails.
"""
function _factor_atoms(p)
    iszero(p) && return typeof(p)[]
    p = normalize_poly(p)
    iszero(p) && return typeof(p)[]
    if _is_nonzero_constant_poly(p)
        return typeof(p)[]
    end
    if !_within_factor_budget(p)
        return typeof(p)[p]
    end
    try
        fac = factor(p)
        atoms = typeof(p)[]
        for (f, _e) in fac
            fn = normalize_poly(f)
            (_is_nonzero_constant_poly(fn) || iszero(fn)) && continue
            push!(atoms, fn)
        end
        return isempty(atoms) ? typeof(p)[p] : atoms
    catch
        return typeof(p)[p]
    end
end

function _is_nonzero_constant_poly(p)
    iszero(p) && return false
    try
        return is_constant(p) && !iszero(p)
    catch
        # degree-based fallback
        try
            return iszero(AbstractAlgebra.total_degree(p)) && !iszero(p)
        catch
            return false
        end
    end
end

function _poly_in(list, p)
    for q in list
        q == p && return true
    end
    return false
end

function _divides_poly(a, b)
    iszero(a) && return false
    try
        return divides(b, a)[1]  # a divides b ?
    catch
        return false
    end
end

# AbstractAlgebra.divides(a, b) means "does b divide a"? Check docs...
# In AA: divides(a, b) -> (true, q) if a = b*q, i.e. b divides a.
function _strict_divides(divisor, dividend)
    iszero(divisor) && return false
    iszero(dividend) && return true
    try
        flag, _ = divides(dividend, divisor)
        return flag
    catch
        return false
    end
end

"""
    _pivot_payload(x) -> (:zero,) | (:constant_nonzero,) | (:poly, p)

Classify a matrix entry for pivot logic. Vanishing depends on the numerator only.
"""
function _pivot_payload(x)
    _is_field_zero(x) && return (:zero,)
    num = _numerator_poly(x)
    if _is_nonzero_constant_poly(num) || (try
        parent(num) isa Field
    catch
        false
    end)
        # constant in poly ring, or element of a parameter-free field
        if parent(num) isa Field
            return iszero(num) ? (:zero,) : (:constant_nonzero,)
        end
        return _is_nonzero_constant_poly(num) ? (:constant_nonzero,) : (:poly, normalize_poly(num))
    end
    p = normalize_poly(num)
    iszero(p) && return (:zero,)
    _is_nonzero_constant_poly(p) && return (:constant_nonzero,)
    return (:poly, p)
end

# --- status -------------------------------------------------------------------

"""
    status(Σ, p) -> PivotStatus

Conservative judgment whether `p` vanishes under `Σ`.
Accepts a polynomial or a matrix entry (uses numerator).
"""
function status(Σ::AssumptionSet, x)
    pay = _pivot_payload(x)
    pay[1] === :zero && return PIVOT_ZERO
    pay[1] === :constant_nonzero && return PIVOT_NONZERO
    return _status_poly(Σ, pay[2])
end

function _status_poly(Σ::AssumptionSet, p)
    iszero(p) && return PIVOT_ZERO
    p = normalize_poly(p)
    iszero(p) && return PIVOT_ZERO
    _is_nonzero_constant_poly(p) && return PIVOT_NONZERO

    # Proven zero: matches an equality, or an equality divides p
    for f in Σ.zeros
        f == p && return PIVOT_ZERO
        _strict_divides(f, p) && return PIVOT_ZERO
    end

    atoms = _factor_atoms(p)
    isempty(atoms) && return PIVOT_NONZERO  # unit / constant

    # Product is ZERO if any atom is proven ZERO
    any_unknown = false
    all_nonzero = true
    for f in atoms
        st = _status_atom(Σ, f)
        if st === PIVOT_ZERO
            return PIVOT_ZERO
        elseif st === PIVOT_UNKNOWN
            any_unknown = true
            all_nonzero = false
        end
    end
    all_nonzero && return PIVOT_NONZERO
    return PIVOT_UNKNOWN
end

function _status_atom(Σ::AssumptionSet, f)
    f = normalize_poly(f)
    iszero(f) && return PIVOT_ZERO
    _is_nonzero_constant_poly(f) && return PIVOT_NONZERO
    for z in Σ.zeros
        z == f && return PIVOT_ZERO
        _strict_divides(z, f) && return PIVOT_ZERO
    end
    for g in Σ.nonzeros
        g == f && return PIVOT_NONZERO
        # g ≠ 0 and g divides f does NOT prove f ≠ 0
        # but f ≠ 0 if f divides g? If f|g and g≠0 then not necessarily f≠0 (g=f*h).
        # If f == g → nonzero. If g divides f and f/g is unit → same.
        if _strict_divides(f, g) && _strict_divides(g, f)
            return PIVOT_NONZERO
        end
    end
    return PIVOT_UNKNOWN
end

# --- assume / consistency -----------------------------------------------------

"""
    assume_zero(Σ, p) -> Union{AssumptionSet, Nothing}

Adjoin `p = 0` with factor-aware reduction. Returns `nothing` if inconsistent.
"""
function assume_zero(Σ::AssumptionSet{P}, x) where {P}
    pay = _pivot_payload(x)
    pay[1] === :constant_nonzero && return nothing
    pay[1] === :zero && return Σ
    return _assume_zero_poly(Σ, pay[2])
end

function _assume_zero_poly(Σ::AssumptionSet{P}, p) where {P}
    p = normalize_poly(p)
    iszero(p) && return Σ
    _is_nonzero_constant_poly(p) && return nothing

    atoms = _factor_atoms(p)
    remaining = P[]
    for f in atoms
        st = _status_atom(Σ, f)
        st === PIVOT_NONZERO && continue
        st === PIVOT_ZERO && return Σ  # already implied
        push!(remaining, f)
    end
    isempty(remaining) && return nothing  # all factors nonzero ⇒ p ≠ 0

    # Single conjunction: product of remaining unknown factors vanishes
    q = remaining[1]
    for i in 2:length(remaining)
        q *= remaining[i]
    end
    q = normalize_poly(q)

    # Conflict with nonvanishings: some g ≠ 0 is forced to 0
    for g in Σ.nonzeros
        if g == q || _strict_divides(q, g)
            return nothing
        end
        # also each atom
        for f in remaining
            if g == f || _strict_divides(f, g)
                return nothing
            end
        end
    end

    new_zeros = copy(Σ.zeros)
    _poly_in(new_zeros, q) || push!(new_zeros, q)
    return AssumptionSet{P}(new_zeros, copy(Σ.nonzeros))
end

"""
    assume_nonzero(Σ, p) -> Union{AssumptionSet, Nothing}

Adjoin `p ≠ 0` (all atomic factors). Returns `nothing` if inconsistent.
"""
function assume_nonzero(Σ::AssumptionSet{P}, x) where {P}
    pay = _pivot_payload(x)
    pay[1] === :zero && return nothing
    pay[1] === :constant_nonzero && return Σ
    return _assume_nonzero_poly(Σ, pay[2])
end

function _assume_nonzero_poly(Σ::AssumptionSet{P}, p) where {P}
    p = normalize_poly(p)
    iszero(p) && return nothing
    _is_nonzero_constant_poly(p) && return Σ

    atoms = _factor_atoms(p)
    isempty(atoms) && return Σ

    new_nz = copy(Σ.nonzeros)
    for f in atoms
        st = _status_atom(Σ, f)
        st === PIVOT_ZERO && return nothing
        st === PIVOT_NONZERO && continue
        _poly_in(new_nz, f) || push!(new_nz, f)
    end
    # Conflict: existing equality forced nonzero
    for z in Σ.zeros
        for f in atoms
            if z == f || _strict_divides(f, z) || _strict_divides(z, f)
                # z=0 and f≠0 with f|z or z|f or equal → conflict when equal or f divides z
                if z == f || _strict_divides(f, z)
                    return nothing
                end
            end
        end
    end
    return AssumptionSet{P}(copy(Σ.zeros), new_nz)
end

"""
    empty_assumptions(R) -> AssumptionSet

Empty Σ with polynomial element type inferred from ring `R` (poly or Frac).
"""
function empty_assumptions(R::Ring)
    if R isa AbstractAlgebra.FracField
        P = elem_type(AbstractAlgebra.base_ring(R))
        return AssumptionSet{P}()
    end
    if R isa Field
        # parameter-free field: use the field element type (constants only)
        return AssumptionSet{elem_type(R)}()
    end
    return AssumptionSet{elem_type(R)}()
end

"""
    assumptions_from_domain(dens) -> AssumptionSet

Initial Σ with domain denominators as nonvanishing conditions.
"""
function assumptions_from_domain(R::Ring, dens::Vector)
    Σ = empty_assumptions(R)
    P = typeof_param(Σ)
    for d in dens
        dP = _coerce_poly(P, d, R)
        Σ2 = assume_nonzero(Σ, dP)
        Σ2 === nothing && throw(ArgumentError("domain denominators inconsistent"))
        Σ = Σ2
    end
    return Σ
end

typeof_param(::AssumptionSet{P}) where {P} = P

function _coerce_poly(::Type{P}, d, R) where {P}
    d isa P && return normalize_poly(d)
    if R isa AbstractAlgebra.FracField
        BR = AbstractAlgebra.base_ring(R)
        return normalize_poly(BR(d))
    end
    # `d` should already live in the parameter polynomial ring
    return normalize_poly(d)
end

"""
    pivot_complexity(p) -> Float64

Symbolic complexity score for pivot selection (lower is better).
"""
function pivot_complexity(x)
    pay = _pivot_payload(x)
    pay[1] === :zero && return Inf
    pay[1] === :constant_nonzero && return 0.0
    p = pay[2]
    deg = try
        float(AbstractAlgebra.total_degree(p))
    catch
        0.0
    end
    nterms = try
        float(length(p))
    catch
        1.0
    end
    nfac = try
        float(length(_factor_atoms(p)))
    catch
        1.0
    end
    return deg + 0.25 * nterms + 0.5 * nfac
end
