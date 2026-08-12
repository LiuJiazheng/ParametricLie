# example/stratify_ab_family.jl
#
# User-facing parametric workflow:
#   family → stratify (classification + jumps)
#         → pick a stratum → specialize → lie/analyze (Levi, …)
#
# CondTree is the internal IR behind stratify; prefer S = stratify(L).
#
# Run from the package root:
#   julia --project=. example/stratify_ab_family.jl

using ParametricLieAlgebras
import AbstractAlgebra
import Nemo

const QQ = Nemo.QQ

println("="^60)
println("ParametricLieAlgebras — stratify → specialize → analyze")
println("="^60)

R, (a, b) = AbstractAlgebra.polynomial_ring(QQ, [:a, :b])
L = LieAlgebra(R, 3, Dict(
    (1, 2) => [R(0), a, R(0)],
    (1, 3) => [R(0), R(0), a - b],
))

println()
println("Family over QQ[a,b]:")
println("  [e1,e2] = a·e2")
println("  [e1,e3] = (a−b)·e3")
println("  [e2,e3] = 0")
println("Jacobi identity: ", check_jacobi(L).ok)

# --- 1. Classification + jumps (user view) -----------------------------------

suite = [:center, :derived_dim, :solvability, :nilpotency]
println()
println("-"^60)
println("1. stratify  — certified parameter regions + jumps vs generic")
println("-"^60)
S = stratify(L; invariants = suite)
show(stdout, MIME("text/plain"), S)
println()
println("IR behind this (not the user table):")
show(stdout, MIME("text/plain"), S.tree)
println()

# --- 2. Pick a stratum and specialize ----------------------------------------

# Region a−b = 0, a ≠ 0  (first confirmed jump in this family)
idx = findfirst(s -> status(s.sigma, a - b) === PIVOT_ZERO &&
                     status(s.sigma, a) === PIVOT_NONZERO, S.strata)
idx === nothing && error("expected a−b=0, a≠0 stratum")
st = S.strata[idx]

point = Dict(a => QQ(1), b => QQ(1))  # satisfies Σ
println("-"^60)
println("2. specialize onto a chosen stratum")
println("-"^60)
println("  stratum Σ = ", st.sigma)
println("  point     = (a,b) = (1,1)")
println("  certified signature: ", join(
    ["$k=$(st.signature[k])" for k in (:center_dim, :derived_dim, :is_solvable, :is_nilpotent)
     if haskey(st.signature, k)], ", "))

Lf = specialize(L, point)
println()
println("  specialize(L, point) → LieAlgebra over ", coefficient_ring(Lf))

# --- 3. Full fiber analysis (Levi, …) on the specialized algebra --------------

println()
println("-"^60)
println("3. analyze(fiber)  — full lie toolkit on the concrete algebra")
println("-"^60)
r = analyze(Lf)
show(stdout, MIME("text/plain"), r)
println()

lev = levi_decomposition(Lf)
println("Levi decomposition on this fiber:")
println("  radical dim = ", dim(lev.radical))
println("  Levi dim    = ", dim(lev.levi))
println("  (this family is solvable on every stratum → Levi is trivial)")

# --- 4. Sanity: fiber matches stratum certificate ----------------------------

println()
println("-"^60)
println("4. validate_stratum  — fiber analyze vs stratum signature")
println("-"^60)
v = validate_stratum(L, st, point)
println("  ok=$(v.ok)  mismatches=$(v.mismatches)")

# Spot-check other regions quickly
for (label, pt) in (
    ("generic (1,2)", Dict(a => QQ(1), b => QQ(2))),
    ("a=0,b=1", Dict(a => QQ(0), b => QQ(1))),
    ("abelian (0,0)", Dict(a => QQ(0), b => QQ(0))),
)
    i = findfirst(S.strata) do s
        try
            return validate_stratum(L, s, pt).ok
        catch
            return false
        end
    end
    i === nothing && continue
    vv = validate_stratum(L, S.strata[i], pt)
    println("  $label  ok=$(vv.ok)")
end

println()
println("="^60)
println("Done.  User path: stratify → specialize → analyze / Levi / …")
println("       CondTree remains available as S.tree when debugging.")
println("  docs: docs/parametric.md")
println("="^60)
