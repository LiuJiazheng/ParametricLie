# Isomorphism helpers: invariant obstruction + search for P ∈ GL(n).
#
# Not a classification engine. `change_of_basis(L, P)` applies a known P;
# here we try to find P or prove none exists via basis-invariant screens /
# explicit matrix search (small n).

"""
    IsoCertificate

Result of [`isomorphism`](@ref) / [`isomorphic`](@ref).

- `isomorphic == true` with `matrix = P` means `change_of_basis(L, P)` matches `L′`
  as structure-constant tensors.
- `isomorphic == false` with `reason == :invariants` is a hard non-isomorphism
  certificate (signature mismatch).
- `reason == :not_found` means no obstruction and no matrix found (incomplete search).
"""
struct IsoCertificate{C<:FieldElem}
    isomorphic::Bool
    matrix::Union{Nothing, MatElem}
    reason::Symbol
    mismatches::Vector{Pair{Symbol,Any}}
end

function Base.show(io::IO, ::MIME"text/plain", c::IsoCertificate)
    if c.isomorphic
        println(io, "IsoCertificate(isomorphic=true, reason=$(c.reason))")
        c.matrix === nothing || println(io, "  P = ", c.matrix)
    else
        println(io, "IsoCertificate(isomorphic=false, reason=$(c.reason))")
        for (k, v) in c.mismatches
            println(io, "  $k: $v")
        end
    end
end

function Base.show(io::IO, c::IsoCertificate)
    print(io, "IsoCertificate(isomorphic=$(c.isomorphic), reason=$(c.reason))")
end

# --- invariant screen -------------------------------------------------------

function _iso_signature(L::LieAlgebra)
    r = analyze(L)
    return Dict{Symbol,Any}(
        :dim => dim(L),
        :center_dim => dim(r.center),
        :derived_dim => dim(derived_algebra(L)),
        :killing_rank => r.killing_rank,
        :radical_dim => dim(r.radical),
        :der_dim => dim(r.derivations),
        :is_solvable => r.solvable,
        :is_nilpotent => r.nilpotent,
        :semisimple => r.semisimple,
        :simple_factor_dims => copy(r.simple_factor_dims),
        :levi_kind => levi_kind(r),
    )
end

function _iso_invariant_mismatches(L::LieAlgebra, L′::LieAlgebra)
    s = _iso_signature(L)
    t = _iso_signature(L′)
    mism = Pair{Symbol,Any}[]
    for k in keys(s)
        haskey(t, k) || continue
        s[k] == t[k] || push!(mism, k => (s[k], t[k]))
    end
    return mism
end

function _structure_constants_equal(a::Array, b::Array)
    size(a) == size(b) || return false
    n = size(a, 1)
    for i in 1:n, j in 1:n, k in 1:n
        a[i, j, k] == b[i, j, k] || return false
    end
    return true
end

"""
    is_structure_isomorphism(L, L′, P) -> Bool

Whether invertible `P` realizes `L ≅ L′` via [`change_of_basis`](@ref).
"""
function is_structure_isomorphism(L::LieAlgebra{C}, L′::LieAlgebra{C}, P) where {C<:FieldElem}
    dim(L) == dim(L′) || return false
    try
        M = _matrix_n(coefficient_ring(L), P, dim(L))
        is_unit(det(M)) || return false
        Lp = change_of_basis(L, M)
        return _structure_constants_equal(structure_constants(Lp), structure_constants(L′))
    catch
        return false
    end
end

# --- polynomial search for P ------------------------------------------------

function _iso_build_equations(L::LieAlgebra{C}, L′::LieAlgebra{C}) where {C<:FieldElem}
    n = dim(L)
    F = coefficient_ring(L)
    R, gensp = polynomial_ring(F, ["p$(i)_$(j)" for i in 1:n for j in 1:n])
    P = zero_matrix(R, n, n)
    for i in 1:n, j in 1:n
        P[i, j] = gensp[(i - 1) * n + j]
    end
    μ = structure_constants(L)
    ν = structure_constants(L′)
    # change_of_basis(L, P) = L′  iff for all p,q,s:
    #   ∑_r P[s,r] ν[p,q,r] = ∑_{i,j} μ[i,j,s] P[i,p] P[j,q]
    eqs = elem_type(R)[]
    for p in 1:n, q in 1:n, s in 1:n
        lhs = zero(R)
        for r in 1:n
            c = ν[p, q, r]
            iszero(c) || (lhs += P[s, r] * R(c))
        end
        rhs = zero(R)
        for i in 1:n, j in 1:n
            c = μ[i, j, s]
            iszero(c) || (rhs += R(c) * P[i, p] * P[j, q])
        end
        e = lhs - rhs
        iszero(e) || push!(eqs, e)
    end
    return R, gensp, eqs
end

function _poly_constant_in_field(e, F)
    iszero(e) && return true, zero(F)
    AbstractAlgebra.total_degree(e) == 0 || return false, zero(F)
    return true, F(AbstractAlgebra.constant_coefficient(e))
end

"""
If `e` is linear in exactly one variable `x ∈ vars_left` with coefficients free of
`vars_left`, return `(x, a, b)` for `a*x + b = 0` with `a,b ∈ F`, `a ≠ 0`.
"""
function _linear_univariate_form(e, vars_left, F)
    R = parent(e)
    iszero(e) && return nothing
    for x in vars_left
        AbstractAlgebra.degree(e, x) == 1 || continue
        a_poly = AbstractAlgebra.derivative(e, x)
        # require deg_x(a) = 0 already from degree(e,x)==1 for monic-in-x sense;
        # also a,b free of all remaining vars
        AbstractAlgebra.degree(a_poly, x) == 0 || continue
        b_poly = AbstractAlgebra.evaluate(e, [x], [zero(R)])
        ok = true
        for y in vars_left
            if AbstractAlgebra.degree(a_poly, y) != 0 || AbstractAlgebra.degree(b_poly, y) != 0
                ok = false
                break
            end
        end
        ok || continue
        _, ac = _poly_constant_in_field(a_poly, F)
        _, bc = _poly_constant_in_field(b_poly, F)
        iszero(ac) && continue
        return (x, ac, bc)
    end
    return nothing
end

function _assignment_to_matrix(gensp, assign, n, F)
    P = zero_matrix(F, n, n)
    for i in 1:n, j in 1:n
        v = gensp[(i - 1) * n + j]
        haskey(assign, v) || return nothing
        P[i, j] = assign[v]
    end
    return P
end

const _ISO_TRIAL_VALUES = (-2, -1, 0, 1, 2, 3)

function _iso_search(L::LieAlgebra{C}, L′::LieAlgebra{C}; max_branch::Int = 4000) where {C<:FieldElem}
    n = dim(L)
    F = coefficient_ring(L)
    R, gensp, eqs0 = _iso_build_equations(L, L′)
    counter = Ref(0)

    function dfs(eqs, assign)
        counter[] += 1
        counter[] > max_branch && return nothing
        neq = elem_type(R)[]
        for e in eqs
            iszero(e) && continue
            okc, cval = _poly_constant_in_field(e, F)
            if okc
                iszero(cval) || return nothing
                continue
            end
            push!(neq, e)
        end
        left = [g for g in gensp if !haskey(assign, g)]
        if isempty(left)
            isempty(neq) || return nothing
            M = _assignment_to_matrix(gensp, assign, n, F)
            M === nothing && return nothing
            is_unit(det(M)) || return nothing
            is_structure_isomorphism(L, L′, M) || return nothing
            return M
        end
        # forced linear
        for e in neq
            lin = _linear_univariate_form(e, left, F)
            lin === nothing && continue
            x, a, b = lin
            val = -b * inv(a)
            assign2 = copy(assign)
            assign2[x] = val
            eqs2 = [AbstractAlgebra.evaluate(f, [x], [R(val)]) for f in neq]
            return dfs(eqs2, assign2)
        end
        # branch
        x = left[1]
        for t in _ISO_TRIAL_VALUES
            assign2 = copy(assign)
            assign2[x] = F(t)
            eqs2 = [AbstractAlgebra.evaluate(f, [x], [R(F(t))]) for f in neq]
            got = dfs(eqs2, assign2)
            got !== nothing && return got
        end
        return nothing
    end

    return dfs(eqs0, Dict{Any,C}())
end

"""
    isomorphism(L, L′) -> IsoCertificate

Try to decide whether `L ≅ L′` over the common coefficient field.

1. Basis-invariant screen (`analyze` dims/flags) — hard negative certificate.
2. Identical structure constants → identity matrix.
3. Small-n polynomial search for `P ∈ GL(n)`.

Incomplete search yields `reason = :not_found` (not a proof of non-isomorphism).
"""
function isomorphism(L::LieAlgebra{C}, L′::LieAlgebra{C}) where {C<:FieldElem}
    if dim(L) != dim(L′)
        return IsoCertificate{C}(false, nothing, :invariants, [:dim => (dim(L), dim(L′))])
    end
    coefficient_ring(L) == coefficient_ring(L′) ||
        throw(ArgumentError("isomorphism requires the same coefficient field"))

    mism = _iso_invariant_mismatches(L, L′)
    isempty(mism) || return IsoCertificate{C}(false, nothing, :invariants, mism)

    F = coefficient_ring(L)
    n = dim(L)
    if _structure_constants_equal(structure_constants(L), structure_constants(L′))
        return IsoCertificate{C}(true, identity_matrix(F, n), :identity, Pair{Symbol,Any}[])
    end

    n > 6 && return IsoCertificate{C}(false, nothing, :not_found, Pair{Symbol,Any}[])

    P = _iso_search(L, L′)
    if P !== nothing
        return IsoCertificate{C}(true, P, :found, Pair{Symbol,Any}[])
    end
    return IsoCertificate{C}(false, nothing, :not_found, Pair{Symbol,Any}[])
end

"""
    isomorphic(L, L′) -> Bool

`true` only with a positive [`IsoCertificate`](@ref). Invariant mismatches give
`false`. Incomplete search also gives `false` (see `isomorphism` for the reason).
"""
isomorphic(L::LieAlgebra{C}, L′::LieAlgebra{C}) where {C<:FieldElem} =
    isomorphism(L, L′).isomorphic
