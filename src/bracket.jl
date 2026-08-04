# Bracket, adjoint, and Jacobi (v0.1)
# Naming: lie_bracket / lie_bracket! aligned with LieGroups.jl.
# Math: [x,y]_k = ∑_{i,j} x_i y_j a_{ijk};  ad_x has (ad_x)[k,j] = ∑_i x_i a_{ijk}.

"""
    JacobiCertificate

Result of [`check_jacobi`](@ref): antisymmetry and Jacobi witnesses.
"""
struct JacobiCertificate
    ok::Bool
    antisym_failing::Vector{Tuple{Int,Int,Int}}  # (i,j,k) with a[i,j,k] + a[j,i,k] ≠ 0
    jacobi_failing::Vector{Tuple{Int,Int,Int}}   # basis triples with nonzero Jacobi sum
end

Base.Bool(c::JacobiCertificate) = c.ok

# Back-compat alias used in early tests / docs
function Base.getproperty(c::JacobiCertificate, s::Symbol)
    s === :failing && return getfield(c, :jacobi_failing)
    return getfield(c, s)
end

# --- coordinate helpers -----------------------------------------------------

function _coords_vector(L::LieAlgebra{C}, x::AbstractVector) where {C<:FieldElem}
    n = dim(L)
    length(x) == n || throw(ArgumentError("expected length $n, got $(length(x))"))
    F = coefficient_ring(L)
    v = Vector{C}(undef, n)
    for i in 1:n
        v[i] = F(x[i])
    end
    return v
end

_coords_vector(::LieAlgebra{C}, x::Vector{C}) where {C<:FieldElem} = x

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
"""
function lie_bracket(L::LieAlgebra{C}, x::LieAlgebraElem{C}, y::LieAlgebraElem{C}) where {C<:FieldElem}
    x.parent === L || throw(ArgumentError("x does not belong to L"))
    y.parent === L || throw(ArgumentError("y does not belong to L"))
    return LieAlgebraElem{C}(L, lie_bracket(L, x.coords, y.coords))
end

function lie_bracket(x::LieAlgebraElem{C}, y::LieAlgebraElem{C}) where {C<:FieldElem}
    L = _check_same_parent(x, y)
    return lie_bracket(L, x, y)
end

function lie_bracket(L::LieAlgebra{C}, x::AbstractVector, y::AbstractVector) where {C<:FieldElem}
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
function lie_bracket!(z::AbstractVector{C}, L::LieAlgebra{C}, x::AbstractVector{C}, y::AbstractVector{C}) where {C<:FieldElem}
    n = dim(L)
    length(z) == n && length(x) == n && length(y) == n ||
        throw(ArgumentError("dimension mismatch"))
    a = structure_constants(L)
    F = coefficient_ring(L)
    for k in 1:n
        z[k] = zero(F)
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
function ad(L::LieAlgebra{C}, x::AbstractVector) where {C<:FieldElem}
    xv = _coords_vector(L, x)
    n = dim(L)
    F = coefficient_ring(L)
    a = structure_constants(L)
    A = zero_matrix(F, n, n)
    for j in 1:n, k in 1:n
        s = zero(F)
        for i in 1:n
            iszero(xv[i]) && continue
            s += xv[i] * a[i, j, k]
        end
        A[k, j] = s
    end
    return A
end

ad(L::LieAlgebra{C}, x::LieAlgebraElem{C}) where {C<:FieldElem} = ad(L, x.coords)

function ad(x::LieAlgebraElem{C}) where {C<:FieldElem}
    return ad(x.parent, x)
end

# --- structure axioms -------------------------------------------------------

"""
    check_antisymmetry(L::LieAlgebra)

Check `a[i,j,k] + a[j,i,k] == 0` for all indices (includes `[e_i,e_i]=0`).

Returns failing triples `(i,j,k)`. Dense constructors do not enforce this;
DOF constructors usually do, but we still verify the stored tensor.
"""
function check_antisymmetry(L::LieAlgebra{C}) where {C<:FieldElem}
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
    check_jacobi(L::LieAlgebra)

Validate Lie bracket axioms on the stored structure constants:

1. **Antisymmetry** (and skew-diagonals): `a[i,j,k] + a[j,i,k] == 0`
   i.e. `[e_i,e_j] = -[e_j,e_i]`, including `[e_i,e_i] = 0`.
2. **Jacobi** on all basis triples:
   `[e_i,[e_j,e_k]] + [e_j,[e_k,e_i]] + [e_k,[e_i,e_j]] = 0`.

Returns a [`JacobiCertificate`](@ref). Important for **dense** input, which does
not auto-enforce antisymmetry; also catches corrupted / inconsistent tensors.
"""
function check_jacobi(L::LieAlgebra{C}) where {C<:FieldElem}
    antisym_failing = check_antisymmetry(L)

    n = dim(L)
    jacobi_failing = Tuple{Int,Int,Int}[]
    # Jacobi is only meaningful if we treat the tensor as a bracket; still report
    # even when antisymmetry fails (helps debugging dense tables).
    for i in 1:n, j in 1:n, k in 1:n
        ei = basis_elem(L, i)
        ej = basis_elem(L, j)
        ek = basis_elem(L, k)
        s = lie_bracket(L, ei, lie_bracket(L, ej, ek)).coords .+
            lie_bracket(L, ej, lie_bracket(L, ek, ei)).coords .+
            lie_bracket(L, ek, lie_bracket(L, ei, ej)).coords
        if !all(iszero, s)
            push!(jacobi_failing, (i, j, k))
        end
    end

    ok = isempty(antisym_failing) && isempty(jacobi_failing)
    return JacobiCertificate(ok, antisym_failing, jacobi_failing)
end
