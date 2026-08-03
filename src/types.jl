# Scalar fields: reuse AbstractAlgebra's parent/element model.
# Do not invent a parallel Field trait — see AbstractAlgebra.Field / FieldElem.
#
# Supported coefficient domains (via AA/Nemo parents):
#   Nemo.QQ, Nemo.GF(p), fraction fields, number fields,
#   Nemo.RealField / Nemo.ComplexField (arb), …
# Julia Base: Rational / AbstractFloat via FieldElement (ad-hoc; prefer AA parents).
#
# Out of v0.1 scope as coefficient *fields*:
#   quaternions / other division rings → AA NCRing (noncommutative); revisit later.

"""
    LieAlgebra{C<:FieldElem}

Finite-dimensional Lie algebra over a coefficient field whose elements have type `C`.

The field itself is the AbstractAlgebra/Nemo **parent** stored in `coefficient_ring`.
Structure constants live in that field (dense table for now; sparse later if needed).
"""
struct LieAlgebra{C<:FieldElem}
    R::Field
    dim::Int
    # struct_consts::Array{C,3}  # [i,j,k] = a_{ijk} with [e_i,e_j] = ∑_k a_{ijk} e_k
    function LieAlgebra{C}(R::Field, n::Int) where {C<:FieldElem}
        elem_type(R) === C || throw(ArgumentError("elem_type(R)=$(elem_type(R)) ≠ C=$C"))
        n >= 0 || throw(ArgumentError("dimension must be nonnegative"))
        new{C}(R, n)
    end
end

"""
    LieAlgebra(R::Field, n::Int)

Lie algebra of dimension `n` over coefficient field `R` (AbstractAlgebra/Nemo parent).
"""
LieAlgebra(R::Field, n::Int) = LieAlgebra{elem_type(R)}(R, n)

dim(L::LieAlgebra) = L.dim

"""
    coefficient_ring(L::LieAlgebra)

Coefficient field parent of `L` (e.g. `Nemo.QQ`, a `FracField`, …).
"""
coefficient_ring(L::LieAlgebra) = L.R

base_ring(L::LieAlgebra) = coefficient_ring(L)  # AA-style alias
