# Bracket, adjoint, and Jacobi.
# Naming: lie_bracket / lie_bracket! aligned with LieGroups.jl.
# Math: [x,y]_k = ∑_{i,j} x_i y_j a_{ijk};  ad_x has (ad_x)[k,j] = ∑_i x_i a_{ijk}.
#
# v0.2.1: Jacobi over free-parameter rings is exact identity testing —
# identity certificate if all residuals vanish; otherwise explicit nonzero
# residuals (no root-finding / stratification).

"""
    JacobiResidual

Explicit nonzero Jacobi residual for one ordered basis triple `(i,j,k)`.

- `residual` — coordinate vector of the Jacobi sum after clearing denominators
  (polynomial identity form)
- `domain` — denominators that must be nonzero for the (rational) expression to
  be defined; empty over polynomial rings / fields without denominators
"""
struct JacobiResidual{C}
    triple::Tuple{Int,Int,Int}
    residual::Vector{C}
    domain::Vector
end

function Base.show(io::IO, r::JacobiResidual)
    print(io, "JacobiResidual(triple=$(r.triple), residual=$(r.residual)")
    isempty(r.domain) || print(io, ", domain≠0: $(r.domain)")
    print(io, ")")
end

"""
    JacobiCertificate

Result of [`check_jacobi`](@ref).

- `ok` — antisymmetry and Jacobi hold identically over the coefficient ring
- `residuals` — explicit nonzero Jacobi residuals (empty when `ok`);
  does **not** solve their zero loci
"""
struct JacobiCertificate
    ok::Bool
    antisym_failing::Vector{Tuple{Int,Int,Int}}
    jacobi_failing::Vector{Tuple{Int,Int,Int}}
    residuals::Vector{JacobiResidual}
end

Base.Bool(c::JacobiCertificate) = c.ok

# Back-compat alias used in early tests / docs
function Base.getproperty(c::JacobiCertificate, s::Symbol)
    s === :failing && return getfield(c, :jacobi_failing)
    return getfield(c, s)
end

function Base.show(io::IO, c::JacobiCertificate)
    if c.ok
        print(io, "JacobiCertificate(identity)")
    else
        print(
            io,
            "JacobiCertificate(nonidentity; antisym=$(length(c.antisym_failing)), ",
            "jacobi=$(length(c.jacobi_failing)), residuals=$(length(c.residuals)))",
        )
    end
end

# --- residual normalization (poly identity / Frac numerators) ---------------

"""
Clear a single coefficient to polynomial residual form.

Returns `(poly_or_elem, domain_denom_or_nothing)`.
Over `Frac(R)`, uses `numerator` / `denominator` and records `denom ≠ 0`.
Over a polynomial ring or ordinary field, returns the element itself.
"""
function _identity_form(x)
    R = parent(x)
    if R isa AbstractAlgebra.FracField
        num = numerator(x)
        den = denominator(x)
        return num, (isone(den) ? nothing : den)
    end
    return x, nothing
end

function _normalize_residual_coords(coords::Vector{C}) where {C}
    n = length(coords)
    # Determine output element type from first cleared numerator
    out = Vector(undef, n)
    domain = Any[]
    seen = Set{String}()
    for i in 1:n
        poly, den = _identity_form(coords[i])
        out[i] = poly
        if den !== nothing
            s = string(den)
            if !(s in seen)
                push!(seen, s)
                push!(domain, den)
            end
        end
    end
    C1 = typeof(out[1])
    residual = C1[out[i] for i in 1:n]
    return residual, domain
end

# --- coordinate helpers -----------------------------------------------------

function _coords_vector(L::LieAlgebra{C}, x::AbstractVector) where {C<:RingElem}
    n = dim(L)
    length(x) == n || throw(ArgumentError("expected length $n, got $(length(x))"))
    R = coefficient_ring(L)
    v = Vector{C}(undef, n)
    for i in 1:n
        v[i] = R(x[i])
    end
    return v
end

_coords_vector(::LieAlgebra{C}, x::Vector{C}) where {C<:RingElem} = x

function _check_same_parent(x::LieAlgebraElem, y::LieAlgebraElem)
    x.parent === y.parent || throw(ArgumentError("elements must belong to the same LieAlgebra"))
    return x.parent
end

# --- lie_bracket ------------------------------------------------------------

"""
    lie_bracket(L::LieAlgebra, x, y)

Lie bracket `[x, y]` via structure constants (triple loop):

    [x,y]_k = ∑_{i,j} x_i y_j a_{ijk}.

Accepts `LieAlgebraElem` or coordinate vectors. Elem inputs return an elem.
Works over fields and free-parameter polynomial / fraction rings.
"""
function lie_bracket(L::LieAlgebra{C}, x::LieAlgebraElem{C}, y::LieAlgebraElem{C}) where {C<:RingElem}
    x.parent === L || throw(ArgumentError("x does not belong to L"))
    y.parent === L || throw(ArgumentError("y does not belong to L"))
    return LieAlgebraElem{C}(L, lie_bracket(L, x.coords, y.coords))
end

function lie_bracket(x::LieAlgebraElem{C}, y::LieAlgebraElem{C}) where {C<:RingElem}
    L = _check_same_parent(x, y)
    return lie_bracket(L, x, y)
end

function lie_bracket(L::LieAlgebra{C}, x::AbstractVector, y::AbstractVector) where {C<:RingElem}
    xv = _coords_vector(L, x)
    yv = _coords_vector(L, y)
    n = dim(L)
    z = fill(zero(coefficient_ring(L)), n)
    lie_bracket!(z, L, xv, yv)
    return z
end

"""
    lie_bracket!(z, L, x, y)

In-place bracket into preallocated coordinate vector `z`.
"""
function lie_bracket!(
    z::AbstractVector{C}, L::LieAlgebra{C}, x::AbstractVector{C}, y::AbstractVector{C}
) where {C<:RingElem}
    n = dim(L)
    length(z) == n && length(x) == n && length(y) == n ||
        throw(ArgumentError("dimension mismatch"))
    a = structure_constants(L)
    R = coefficient_ring(L)
    for k in 1:n
        z[k] = zero(R)
    end
    for i in 1:n
        iszero(x[i]) && continue
        for j in 1:n
            iszero(y[j]) && continue
            xiyj = x[i] * y[j]
            for k in 1:n
                aik = a[i, j, k]
                iszero(aik) && continue
                z[k] += xiyj * aik
            end
        end
    end
    return z
end

# --- adjoint ----------------------------------------------------------------

"""
    ad(L::LieAlgebra, x) -> MatElem

Adjoint matrix of `x`: AbstractAlgebra matrix `A_x` with
`(A_x * y)_k = [x,y]_k`, i.e. `A_x[k,j] = sum_i x_i * a[i,j,k]`.

Computed on demand from structure constants (not cached on `L`).
"""
function ad(L::LieAlgebra{C}, x::AbstractVector) where {C<:RingElem}
    xv = _coords_vector(L, x)
    n = dim(L)
    R = coefficient_ring(L)
    a = structure_constants(L)
    A = zero_matrix(R, n, n)
    for j in 1:n, k in 1:n
        s = zero(R)
        for i in 1:n
            iszero(xv[i]) && continue
            s += xv[i] * a[i, j, k]
        end
        A[k, j] = s
    end
    return A
end

ad(L::LieAlgebra{C}, x::LieAlgebraElem{C}) where {C<:RingElem} = ad(L, x.coords)

function ad(x::LieAlgebraElem{C}) where {C<:RingElem}
    return ad(x.parent, x)
end

# --- structure axioms -------------------------------------------------------

"""
    check_antisymmetry(L::LieAlgebra)

Check `a[i,j,k] + a[j,i,k] == 0` for all indices (includes `[e_i,e_i]=0`).

Returns failing triples `(i,j,k)`. Dense constructors do not enforce this;
DOF constructors usually do, but we still verify the stored tensor.
"""
function check_antisymmetry(L::LieAlgebra{C}) where {C<:RingElem}
    n = dim(L)
    a = structure_constants(L)
    failing = Tuple{Int,Int,Int}[]
    for i in 1:n, j in 1:n, k in 1:n
        if !iszero(a[i, j, k] + a[j, i, k])
            push!(failing, (i, j, k))
        end
    end
    return failing
end

"""
    check_jacobi(L::LieAlgebra) -> JacobiCertificate

Validate Lie bracket axioms by **exact identity** over the coefficient ring:

1. **Antisymmetry:** `a[i,j,k] + a[j,i,k] == 0`
2. **Jacobi** on all basis triples.

Over a free-parameter polynomial ring `k[a₁,…]`, parameters are algebraically
independent: `f = 0` means `f` is the zero polynomial. Over `Frac(R)`,
`p/q = 0` means `p = 0` (with domain condition `q ≠ 0` recorded on residuals).

On success returns an **identity** certificate (`ok = true`, empty residuals).
On failure returns **explicit nonzero residuals** — it does **not** solve their
zero sets or stratify parameter space (v0.2.3+).
"""
function check_jacobi(L::LieAlgebra{C}) where {C<:RingElem}
    antisym_failing = check_antisymmetry(L)

    n = dim(L)
    jacobi_failing = Tuple{Int,Int,Int}[]
    residuals = JacobiResidual[]

    for i in 1:n, j in 1:n, k in 1:n
        ei = basis_elem(L, i)
        ej = basis_elem(L, j)
        ek = basis_elem(L, k)
        s = lie_bracket(L, ei, lie_bracket(L, ej, ek)).coords .+
            lie_bracket(L, ej, lie_bracket(L, ek, ei)).coords .+
            lie_bracket(L, ek, lie_bracket(L, ei, ej)).coords
        if !all(iszero, s)
            push!(jacobi_failing, (i, j, k))
            res_coords, domain = _normalize_residual_coords(s)
            # Skip if clearing made everything zero (should not happen if iszero worked)
            if !all(iszero, res_coords)
                push!(residuals, JacobiResidual( (i, j, k), res_coords, domain ))
            end
        end
    end

    ok = isempty(antisym_failing) && isempty(jacobi_failing)
    return JacobiCertificate(ok, antisym_failing, jacobi_failing, residuals)
end
