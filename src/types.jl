# Scalar fields: reuse AbstractAlgebra's parent/element model.
# Coefficient field is named `F` (not `R`) to avoid reading it as ℝ.
#
# Supported coefficient domains (via AA/Nemo parents):
#   Nemo.QQ, Nemo.GF(p), fraction fields, number fields,
#   Nemo.RealField / Nemo.ComplexField (arb), …
#
# Out of v0.1 scope as coefficient *fields*:
#   quaternions / other division rings → AA NCRing; revisit later.
#
# Storage vs input:
#   - Internally we keep a dense n×n×n tensor (simple algorithms).
#   - Dense ctor:  LieAlgebra(F, c::Array{C,3})
#   - DOF ctor:    LieAlgebra(F, n, brackets)  # only i < j (or sparse pairs)
#   - Abelian:     LieAlgebra(F, n)

"""
    LieAlgebra{C<:FieldElem}

Finite-dimensional Lie algebra over a coefficient field `F`, presented on the
coordinate space `F^n` with the standard basis `e₁,…,eₙ`.

The algebraic structure is the structure-constant tensor
`struct_consts[i,j,k] = aᵢⱼₖ` with

    [eᵢ, eⱼ] = ∑ₖ aᵢⱼₖ eₖ.

(The basis itself is not stored: it is the standard basis of `F^n`.)
"""
struct LieAlgebra{C<:FieldElem}
    F::Field
    struct_consts::Array{C,3}
    function LieAlgebra{C}(F::Field, struct_consts::Array{C,3}) where {C<:FieldElem}
        elem_type(F) === C || throw(ArgumentError("elem_type(F)=$(elem_type(F)) ≠ C=$C"))
        n1, n2, n3 = size(struct_consts)
        (n1 == n2 == n3) || throw(ArgumentError("structure constants must be n×n×n, got $(size(struct_consts))"))
        n1 >= 0 || throw(ArgumentError("dimension must be nonnegative"))
        new{C}(F, struct_consts)
    end
end

# --- dense path -------------------------------------------------------------

"""
    LieAlgebra(F::Field, struct_consts::Array{C,3}) where {C<:FieldElem}

**Dense** constructor: full `n×n×n` tensor (includes antisymmetric duplicates).
"""
function LieAlgebra(F::Field, struct_consts::Array{C,3}) where {C<:FieldElem}
    return LieAlgebra{C}(F, struct_consts)
end

"""
    LieAlgebra(F::Field, n::Int)

Abelian Lie algebra of dimension `n` over `F` (zero structure constants).
"""
function LieAlgebra(F::Field, n::Int)
    n >= 0 || throw(ArgumentError("dimension must be nonnegative"))
    C = elem_type(F)
    return LieAlgebra{C}(F, fill(zero(F), n, n, n))
end

# --- DOF / sparse-pair path -------------------------------------------------

"""
    LieAlgebra(F::Field, n::Int, brackets)

**Degree-of-freedom** constructor: specify only independent brackets.

`brackets` is an iterator of `pair => coords` (e.g. `Dict` or vector of `Pair`s):

```julia
# Heisenberg: [e₁,e₂] = e₃  (only i < j needed)
LieAlgebra(Nemo.QQ, 3, Dict(
    (1, 2) => [0, 0, 1],
))
```

- Keys `(i,j)` are 1-based basis indices; if `i > j`, the entry is treated as `[eᵢ,eⱼ]`
  and stored accordingly (you may also pass only `i < j`).
- `coords` is a length-`n` vector: coefficients of `[eᵢ,eⱼ]` in the standard basis.
- Diagonal brackets default to zero; opposite pairs are filled by antisymmetry.
- Skipped pairs default to zero (sparse-friendly).

Internally this still builds a dense tensor for v0.1 algorithms.
"""
function LieAlgebra(F::Field, n::Int, brackets)
    n >= 0 || throw(ArgumentError("dimension must be nonnegative"))
    C = elem_type(F)
    c = _struct_consts_from_brackets(F, n, brackets)
    return LieAlgebra{C}(F, c)
end

function _struct_consts_from_brackets(F::Field, n::Int, brackets)
    C = elem_type(F)
    c = fill(zero(F), n, n, n)
    # Track which unordered pairs were set, to catch contradictory duplicates.
    seen = Dict{Tuple{Int,Int},Vector{C}}()

    for br in brackets
        br isa Pair || throw(ArgumentError("each bracket must be (i,j) => coords, got $br"))
        ij, coords = br
        i, j = _pair_indices(ij, n)
        v = _coords_in_F(F, n, coords)

        if i == j
            iszero(v) || throw(ArgumentError("[e_i,e_i] must be zero; got nonzero for i=$i"))
            continue
        end

        key = i < j ? (i, j) : (j, i)
        v_low = i < j ? v : .-v   # store canonical [e_min, e_max]

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

function _coords_in_F(F::Field, n::Int, coords)
    length(coords) == n || throw(ArgumentError("coords length must be n=$n, got $(length(coords))"))
    C = elem_type(F)
    v = Vector{C}(undef, n)
    for k in 1:n
        v[k] = F(coords[k])
    end
    return v
end

# --- accessors --------------------------------------------------------------

dim(L::LieAlgebra) = size(L.struct_consts, 1)

"""
    coefficient_ring(L::LieAlgebra)

Coefficient field parent `F` of `L` (e.g. `Nemo.QQ`, a `FracField`, …).
"""
coefficient_ring(L::LieAlgebra) = L.F

base_ring(L::LieAlgebra) = coefficient_ring(L)  # AA-style alias

"""
    structure_constants(L::LieAlgebra)

The `n×n×n` tensor `aᵢⱼₖ` for the current (standard) basis of `F^n`.
"""
structure_constants(L::LieAlgebra) = L.struct_consts

# --- elements ---------------------------------------------------------------

"""
    LieAlgebraElem{C}

Element of a [`LieAlgebra`](@ref), stored as coordinates in `F^n` with a parent.
"""
struct LieAlgebraElem{C<:FieldElem}
    parent::LieAlgebra{C}
    coords::Vector{C}
    function LieAlgebraElem{C}(L::LieAlgebra{C}, coords::Vector{C}) where {C<:FieldElem}
        length(coords) == dim(L) ||
            throw(ArgumentError("expected $(dim(L)) coordinates, got $(length(coords))"))
        new{C}(L, coords)
    end
end

function LieAlgebraElem(L::LieAlgebra{C}, coords::AbstractVector) where {C<:FieldElem}
    F = coefficient_ring(L)
    v = Vector{C}(undef, dim(L))
    length(coords) == dim(L) ||
        throw(ArgumentError("expected $(dim(L)) coordinates, got $(length(coords))"))
    for i in eachindex(v)
        v[i] = F(coords[i])
    end
    return LieAlgebraElem{C}(L, v)
end

Base.parent(x::LieAlgebraElem) = x.parent
coords(x::LieAlgebraElem) = x.coords
dim(x::LieAlgebraElem) = dim(x.parent)
coefficient_ring(x::LieAlgebraElem) = coefficient_ring(x.parent)

function Base.zero(L::LieAlgebra{C}) where {C<:FieldElem}
    return LieAlgebraElem{C}(L, fill(zero(coefficient_ring(L)), dim(L)))
end

Base.zero(x::LieAlgebraElem) = zero(parent(x))

"""
    basis_elem(L::LieAlgebra, i::Int)

The `i`-th standard basis vector as a [`LieAlgebraElem`](@ref).
"""
function basis_elem(L::LieAlgebra{C}, i::Int) where {C<:FieldElem}
    n = dim(L)
    (1 <= i <= n) || throw(ArgumentError("basis index $i out of range 1:$n"))
    v = fill(zero(coefficient_ring(L)), n)
    v[i] = one(coefficient_ring(L))
    return LieAlgebraElem{C}(L, v)
end

function Base.:(==)(x::LieAlgebraElem{C}, y::LieAlgebraElem{C}) where {C<:FieldElem}
    return x.parent === y.parent && x.coords == y.coords
end

function Base.show(io::IO, x::LieAlgebraElem)
    print(io, "LieAlgebraElem(dim=$(dim(x)), coords=$(x.coords))")
end
