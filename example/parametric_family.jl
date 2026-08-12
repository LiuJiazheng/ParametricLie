# example/parametric_family.jl
#
# v0.2.1 + v0.2.2 use case: free-parameter families.
#
#   v0.2.1  — Does the parametric bracket identically satisfy Jacobi?
#   v0.2.2  — What does the generic fiber look like, and what about a
#             concrete fiber the user specializes to?
#
# Run from the package root:
#   julia --project=. example/parametric_family.jl

using ParametricLieAlgebras
import AbstractAlgebra
import Nemo

const QQ = Nemo.QQ

println("="^60)
println("ParametricLieAlgebras example — v0.2.1 / v0.2.2 parameter families")
println("="^60)

R, (a,) = AbstractAlgebra.polynomial_ring(QQ, [:a])
println()
println("Parameter ring: QQ[a]  (a free / algebraically independent)")

# =============================================================================
# v0.2.1 — Exact Jacobi identity over the polynomial ring
# =============================================================================

println()
println("="^60)
println("v0.2.1 — Jacobi as polynomial identity (no root-finding)")
println("="^60)

# Identity case: parametric Heisenberg  [e1,e2] = a·e3
H_a = LieAlgebra(R, 3, Dict((1, 2) => [R(0), R(0), a]))
println()
println("Family H_a:  [e1, e2] = a·e3   over QQ[a]")
println("parameters(H_a) = ", parameters(H_a))
jac_ok = check_jacobi(H_a)
println("check_jacobi → ", jac_ok)
println("  ok = ", jac_ok.ok, "  (identity certificate: holds for all a)")

# Non-identity case: explicit nonzero residual
L_bad = LieAlgebra(R, 3, Dict(
    (1, 2) => [a, R(0), R(0)],
    (1, 3) => [R(0), a, R(0)],
    (2, 3) => [R(0), R(0), R(1)],
))
jac_bad = check_jacobi(L_bad)
println()
println("Broken family (for contrast): brackets that fail Jacobi identically")
println("check_jacobi → ", jac_bad)
if !isempty(jac_bad.residuals)
    r0 = first(jac_bad.residuals)
    println("  first residual triple ", r0.triple)
    println("  residual coords      ", r0.residual)
    println("  (returned as a polynomial — we do NOT solve its zeros)")
end

# =============================================================================
# v0.2.2 — Generic analysis vs user-chosen specialization
# =============================================================================

println()
println("="^60)
println("v0.2.2 — generic fiber vs specialize(a => …)")
println("="^60)

# Core demo family: [e1, e2] = a·e2
L_a = LieAlgebra(R, 2, Dict((1, 2) => [R(0), a]))
println()
println("Family L_a:  [e1, e2] = a·e2   over QQ[a]")
println()
println("Tree:")
println("  family L_a")
println("   ├── analyze_generic  →  invariants over QQ(a)")
println("   └── specialize(a=>α) →  concrete fiber over QQ")
println()
println("v0.2.2 does NOT auto-discover that a=0 is exceptional.")

# --- generic ---
println()
println("-"^60)
println("analyze_generic(L_a)  — view over QQ(a)")
println("-"^60)
println()
r_gen = analyze_generic(L_a)
show(stdout, MIME("text/plain"), r_gen)
println()
println()
println("Generic summary (basis-invariant):")
println("  center dim     = ", dim(center(r_gen)))
println("  derived dim    = ", dim(derived_algebra(generic_algebra(L_a))))
println("  solvable       = ", is_solvable(r_gen))
println("  nilpotent      = ", is_nilpotent(r_gen))
println("  (same as non-abelian 2D algebra when a ≠ 0, e.g. a=1)")

# --- fiber a = 0 (user chose this point) ---
println()
println("-"^60)
println("specialize(L_a, a => 0)  — fiber chosen by the user")
println("-"^60)
L0 = specialize(L_a, Dict(a => QQ(0)))
println()
println("L0 structure constants zero? ", all(iszero, structure_constants(L0)))
r0 = analyze(L0)
println()
show(stdout, MIME("text/plain"), r0)
println()
println()
println("Fiber a=0 summary:")
println("  center dim     = ", dim(center(r0)))
println("  derived dim    = ", dim(derived_algebra(L0)))
println("  solvable       = ", is_solvable(r0))
println("  nilpotent      = ", is_nilpotent(r0))

# --- fiber a = 1 ---
println()
println("-"^60)
println("specialize(L_a, a => 1)  — another concrete fiber")
println("-"^60)
L1 = specialize(L_a, Dict(:a => 1))
r1 = analyze(L1)
println()
println("Fiber a=1: center=$(dim(center(r1))), nilpotent=$(is_nilpotent(r1))")
println("  (matches generic qualitative invariants)")

# --- contrast table ---
println()
println("-"^60)
println("What we demonstrated")
println("-"^60)
println()
println("  view              center   derived   nilpotent")
println("  ----------------  ------   -------   ---------")
println("  generic QQ(a)        ", dim(center(r_gen)), "        ",
        dim(derived_algebra(generic_algebra(L_a))), "       ", is_nilpotent(r_gen))
println("  specialize a=0       ", dim(center(r0)), "        ",
        dim(derived_algebra(L0)), "       ", is_nilpotent(r0))
println("  specialize a=1       ", dim(center(r1)), "        ",
        dim(derived_algebra(L1)), "       ", is_nilpotent(r1))
println()
println("generic report ≠ fiber a=0  →  invariants can jump")
println("finding that jump automatically  →  v0.2.3+ (not this milestone)")

println()
println("="^60)
println("Done.")
println("  docs: docs/parametric.md")
println("="^60)
