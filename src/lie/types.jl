# Scalar domains: reuse AbstractAlgebra's parent/element model.
#
# Supported coefficient domains (via AA/Nemo parents):
#   Fields:  Nemo.QQ, Nemo.GF(p), Frac(R), number fields, …
#   Rings:   multivariate (and uni) polynomial rings k[a₁,…,aₘ]  (v0.2.1)
#
# Parameters in a polynomial ring are treated as **algebraically independent**.
# Equality is polynomial (resp. rational) identity — not solving zero sets.
#
# Storage vs input:
#   - Internally we keep a dense n×n×n tensor (simple algorithms).
#   - Dense ctor:  LieAlgebra(R, c::Array{C,3})
#   - DOF ctor:    LieAlgebra(R, n, brackets)
#   - Abelian:     LieAlgebra(R, n)

"""
    LieAlgebra{C<:RingElem}

Finite-dimensional Lie algebra over a coefficient **ring** `R` (often a field),
presented on the coordinate space `R^n` with the standard basis `e₁,…,eₙ`.

The algebraic structure is the structure-constant tensor
`struct_consts[i,j,k] = aᵢⱼₖ` with

    [eᵢ, eⱼ] = ∑ₖ aᵢⱼₖ eₖ.

Optional metadata (v0.2.1):

- `parameters` — free indeterminates when `R = k[a₁,…]` or `Frac(k[a₁,…])`
- `domain_denominators` — polynomials that must be nonzero for a rational-function
  family to be defined (denominators appearing in structure constants)

Field-only algorithms (center nullspace, Killing, Der, Levi, …) require
`C <: FieldElem` and simply do not apply to plain polynomial rings.
"""
struct LieAlgebra{C<:RingElem}
    F::Ring
    struct_consts::Array{C,3}
    parameters::Vector
    domain_denominators::Vector
    function LieAlgebra{C}(
        F::Ring,
        struct_consts::Array{C,3};
        parameters::Vector = Any[],
        domain_denominators::Vector = Any[],
    ) where {C<:RingElem}
        elem_type(F) === C || throw(ArgumentError("elem_type(F)=$(elem_type(F)) ≠ C=$C"))
        n1, n2, n3 = size(struct_consts)
        (n1 == n2 == n3) ||
            throw(ArgumentError("structure constants must be n×n×n, got $(size(struct_consts))"))
        n1 >= 0 || throw(ArgumentError("dimension must be nonnegative"))
        new{C}(F, struct_consts, parameters, domain_denominators)
    end
end

# --- parametric metadata helpers --------------------------------------------

"""
    _parameter_gens(R::Ring) -> Vector

Indeterminates of a polynomial ring, or of `base_ring` of a fraction field;
empty for ordinary fields such as `QQ`.
"""
function _parameter_gens(R::Ring)
    if R isa AbstractAlgebra.FracField
        return _parameter_gens(AbstractAlgebra.base_ring(R))
    end
    try
        return collect(AbstractAlgebra.gens(R))
    catch
        return Any[]
    end
end

"""
    _collect_denominators(R, struct_consts) -> Vector

Unique nontrivial denominators of rational structure constants (for domain metadata).
Empty when coefficients are not in a fraction field.
"""
function _collect_denominators(R::Ring, struct_consts::Array{C,3}) where {C}
    R isa AbstractAlgebra.FracField || return Any[]
    dens = Any[]
    seen = Set{Any}()
    for x in struct_consts
        d = denominator(x)
        isone(d) && continue
        s = string(d)
        if !(s in seen)
            push!(seen, s)
            push!(dens, d)
        end
    end
    return dens
end

function _finish_lie_algebra(F::Ring, c::Array{C,3}) where {C}
    params = _parameter_gens(F)
    dens = _collect_denominators(F, c)
    return LieAlgebra{C}(F, c; parameters = params, domain_denominators = dens)
end

# --- dense path -------------------------------------------------------------

"""
    LieAlgebra(R::Ring, struct_consts::Array{C,3}) where {C<:RingElem}

**Dense** constructor: full `n×n×n` tensor (includes antisymmetric duplicates).
`R` may be a field or a polynomial / fraction ring.
"""
function LieAlgebra(R::Ring, struct_consts::Array{C,3}) where {C<:RingElem}
    return _finish_lie_algebra(R, struct_consts)
end

"""
    LieAlgebra(R::Ring, n::Int)

Abelian Lie algebra of dimension `n` over `R` (zero structure constants).
"""
function LieAlgebra(R::Ring, n::Int)
    n >= 0 || throw(ArgumentError("dimension must be nonnegative"))
    C = elem_type(R)
    return _finish_lie_algebra(R, fill(zero(R), n, n, n))
end

# --- DOF / sparse-pair path -------------------------------------------------

"""
    LieAlgebra(R::Ring, n::Int, brackets)

**Degree-of-freedom** constructor: specify only independent brackets.

`brackets` is an iterator of `pair => coords` (e.g. `Dict` or vector of `Pair`s):

```julia
# Parametric Heisenberg: [e₁,e₂] = a·e₃ over QQ[a]
R, (a,) = polynomial_ring(Nemo.QQ, [:a])
LieAlgebra(R, 3, Dict((1, 2) => [0, 0, a]))
```

- Keys `(i,j)` are 1-based basis indices; if `i > j`, the entry is treated as `[eᵢ,eⱼ]`
  and stored accordingly (you may also pass only `i < j`).
- `coords` is a length-`n` vector: coefficients of `[eᵢ,eⱼ]` in the standard basis.
- Diagonal brackets default to zero; opposite pairs are filled by antisymmetry.
- Skipped pairs default to zero (sparse-friendly).
"""
function LieAlgebra(R::Ring, n::Int, brackets)
    n >= 0 || throw(ArgumentError("dimension must be nonnegative"))
    c = _struct_consts_from_brackets(R, n, brackets)
    return _finish_lie_algebra(R, c)
end

function _struct_consts_from_brackets(R::Ring, n::Int, brackets)
    C = elem_type(R)
    c = fill(zero(R), n, n, n)
    seen = Dict{Tuple{Int,Int},Vector{C}}()

    for br in brackets
        br isa Pair || throw(ArgumentError("each bracket must be (i,j) => coords, got $br"))
        ij, coords = br
        i, j = _pair_indices(ij, n)
        v = _coords_in_R(R, n, coords)

        if i == j
            iszero(v) || throw(ArgumentError("[e_i,e_i] must be zero; got nonzero for i=$i"))
            continue
        end

        key = i < j ? (i, j) : (j, i)
        v_low = i < j ? v : .-v

        if haskey(seen, key)
            seen[key] == v_low || throw(ArgumentError("contradictory brackets for pair $key"))
            continue
        end
        seen[key] = v_low

        a, b = key
        for k in 1:n
            c[a, b, k] = v_low[k]
            c[b, a, k] = -v_low[k]
        end
    end
    return c
end

function _pair_indices(ij, n::Int)
    if ij isa Tuple && length(ij) == 2
        i, j = Int(ij[1]), Int(ij[2])
    else
        throw(ArgumentError("bracket key must be (i,j), got $ij"))
    end
    (1 <= i <= n && 1 <= j <= n) || throw(ArgumentError("indices ($i,$j) out of range 1:$n"))
    return i, j
end

function _coords_in_R(R::Ring, n::Int, coords)
    length(coords) == n || throw(ArgumentError("coords length must be n=$n, got $(length(coords))"))
    C = elem_type(R)
    v = Vector{C}(undef, n)
    for k in 1:n
        v[k] = R(coords[k])
    end
    return v
end

# --- accessors --------------------------------------------------------------

dim(L::LieAlgebra) = size(L.struct_consts, 1)

"""
    coefficient_ring(L::LieAlgebra)

Coefficient ring parent of `L` (e.g. `Nemo.QQ`, `QQ[a,b]`, `Frac(QQ[a])`, …).
"""
coefficient_ring(L::LieAlgebra) = L.F

base_ring(L::LieAlgebra) = coefficient_ring(L)  # AA-style alias

"""
    parameters(L::LieAlgebra)

Free parameter indeterminates attached to `L` (empty for concrete fields).
"""
parameters(L::LieAlgebra) = L.parameters

"""
    domain_denominators(L::LieAlgebra)

Polynomials that must be nonzero for a rational-function family to be defined.
Empty for polynomial rings and ordinary fields.
"""
domain_denominators(L::LieAlgebra) = L.domain_denominators

"""
    structure_constants(L::LieAlgebra)

The `n×n×n` tensor `aᵢⱼₖ` for the current (standard) basis of `R^n`.
"""
structure_constants(L::LieAlgebra) = L.struct_consts

# --- elements ---------------------------------------------------------------

"""
    LieAlgebraElem{C}

Element of a [`LieAlgebra`](@ref), stored as coordinates in `R^n` with a parent.
"""
struct LieAlgebraElem{C<:RingElem}
    parent::LieAlgebra{C}
    coords::Vector{C}
    function LieAlgebraElem{C}(L::LieAlgebra{C}, coords::Vector{C}) where {C<:RingElem}
        length(coords) == dim(L) ||
            throw(ArgumentError("expected $(dim(L)) coordinates, got $(length(coords))"))
        new{C}(L, coords)
    end
end

function LieAlgebraElem(L::LieAlgebra{C}, coords::AbstractVector) where {C<:RingElem}
    R = coefficient_ring(L)
    v = Vector{C}(undef, dim(L))
    length(coords) == dim(L) ||
        throw(ArgumentError("expected $(dim(L)) coordinates, got $(length(coords))"))
    for i in eachindex(v)
        v[i] = R(coords[i])
    end
    return LieAlgebraElem{C}(L, v)
end

Base.parent(x::LieAlgebraElem) = x.parent
coords(x::LieAlgebraElem) = x.coords
dim(x::LieAlgebraElem) = dim(x.parent)
coefficient_ring(x::LieAlgebraElem) = coefficient_ring(x.parent)

function Base.zero(L::LieAlgebra{C}) where {C<:RingElem}
    return LieAlgebraElem{C}(L, fill(zero(coefficient_ring(L)), dim(L)))
end

Base.zero(x::LieAlgebraElem) = zero(parent(x))

"""
    basis_elem(L::LieAlgebra, i::Int)

The `i`-th standard basis vector as a [`LieAlgebraElem`](@ref).
"""
function basis_elem(L::LieAlgebra{C}, i::Int) where {C<:RingElem}
    n = dim(L)
    (1 <= i <= n) || throw(ArgumentError("basis index $i out of range 1:$n"))
    R = coefficient_ring(L)
    v = fill(zero(R), n)
    v[i] = one(R)
    return LieAlgebraElem{C}(L, v)
end

function Base.:(==)(x::LieAlgebraElem{C}, y::LieAlgebraElem{C}) where {C<:RingElem}
    return x.parent === y.parent && x.coords == y.coords
end

function Base.show(io::IO, x::LieAlgebraElem)
    print(io, "LieAlgebraElem(dim=$(dim(x)), coords=$(x.coords))")
end
