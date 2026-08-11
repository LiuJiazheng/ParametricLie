# Coefficient modules for Chevalley–Eilenberg cohomology.
# A module M is given by matrices ρ(e_i) ∈ End(M), i = 1,…,n.

"""
    LieModule{C}

Finite-dimensional left module over a [`LieAlgebra`](@ref) `L`, presented on
`F^m` with action matrices `ρ(e_i)` (each `m×m`, columns = images of the
standard basis of `M`).
"""
struct LieModule{C<:FieldElem}
    parent::LieAlgebra{C}
    m::Int
    actions::Vector{MatElem}  # length n; each m×m
    function LieModule{C}(L::LieAlgebra{C}, m::Int, actions::Vector) where {C<:FieldElem}
        n = dim(L)
        m >= 0 || throw(ArgumentError("module dimension must be nonnegative, got $m"))
        length(actions) == n ||
            throw(ArgumentError("expected $n action matrices, got $(length(actions))"))
        F = coefficient_ring(L)
        acts = MatElem[]
        for (i, A) in enumerate(actions)
            size(A) == (m, m) ||
                throw(ArgumentError("action[$i] must be $m×$m, got $(size(A))"))
            M = zero_matrix(F, m, m)
            for r in 1:m, c in 1:m
                M[r, c] = F(A[r, c])
            end
            push!(acts, M)
        end
        return new{C}(L, m, acts)
    end
end

function LieModule(L::LieAlgebra{C}, m::Int, actions::Vector) where {C<:FieldElem}
    return LieModule{C}(L, m, actions)
end

Base.parent(M::LieModule) = M.parent
dim(M::LieModule) = M.m
action_matrices(M::LieModule) = M.actions
coefficient_ring(M::LieModule) = coefficient_ring(parent(M))

function Base.show(io::IO, M::LieModule)
    print(io, "LieModule(dim=$(dim(M)) over LieAlgebra(dim=$(dim(parent(M)))))")
end

"""
    trivial_module(L::LieAlgebra) -> LieModule

One-dimensional trivial module `F` (all actions zero).
"""
function trivial_module(L::LieAlgebra{C}) where {C<:FieldElem}
    n = dim(L)
    F = coefficient_ring(L)
    z = zero_matrix(F, 1, 1)
    return LieModule{C}(L, 1, [deepcopy(z) for _ in 1:n])
end

"""
    adjoint_module(L::LieAlgebra) -> LieModule

The adjoint module: `M = L` with `x · y = [x,y]`.
"""
function adjoint_module(L::LieAlgebra{C}) where {C<:FieldElem}
    n = dim(L)
    F = coefficient_ring(L)
    if n == 0
        return LieModule{C}(L, 0, MatElem[])
    end
    acts = MatElem[ad(L, basis_elem(L, i)) for i in 1:n]
    return LieModule{C}(L, n, acts)
end

"""
    act(M::LieModule, i::Int, v) -> Vector

Action of the basis vector `e_i` on a coordinate vector in `M`.
"""
function act(M::LieModule{C}, i::Int, v::AbstractVector) where {C<:FieldElem}
    m = dim(M)
    length(v) == m || throw(ArgumentError("expected length $m, got $(length(v))"))
    n = dim(parent(M))
    (1 <= i <= n) || throw(ArgumentError("basis index $i out of range 1:$n"))
    F = coefficient_ring(M)
    A = M.actions[i]
    out = Vector{C}(undef, m)
    for r in 1:m
        s = zero(F)
        for c in 1:m
            s += A[r, c] * F(v[c])
        end
        out[r] = s
    end
    return out
end
