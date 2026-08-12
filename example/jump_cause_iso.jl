# example/jump_cause_iso.jl
#
# End-to-end narrative:
#   stratify  →  Jump (e.g. center 0→1)
#            →  explain_jumps!  (H² / wall / MC cause)
#            →  isomorphism on fiber representatives
#
# Run from the package root:
#   julia --project=. example/jump_cause_iso.jl

using ParametricLie
import AbstractAlgebra
import Nemo

const QQ = Nemo.QQ

println("="^64)
println("ParametricLie — Jump → cause (H²/MC) → isomorphism")
println("="^64)

R, (a, b) = AbstractAlgebra.polynomial_ring(QQ, [:a, :b])
L = LieAlgebra(R, 3, Dict(
    (1, 2) => [R(0), a, R(0)],
    (1, 3) => [R(0), R(0), a - b],
))

println()
println("Family over QQ[a,b]:")
println("  [e1,e2] = a·e2")
println("  [e1,e3] = (a−b)·e3")

# --- 1. Stratify ------------------------------------------------------------

S = stratify(L; invariants = [:center, :derived_dim, :solvability, :nilpotency])
println()
println("-"^64)
println("1. stratify — certified regions + jumps")
println("-"^64)
show(stdout, MIME("text/plain"), S)
println()

# Focus on the wall a−b = 0, a ≠ 0  (center 0→1)
jumps = jump_table(S)
jidx = findfirst(j -> status(j.target.sigma, a - b) === PIVOT_ZERO &&
                      status(j.target.sigma, a) === PIVOT_NONZERO, jumps)
jidx === nothing && error("expected jump across a−b = 0")
J = jumps[jidx]

println("Selected jump:")
show(stdout, MIME("text/plain"), J)
println()

# --- 2. Explain with H² / wall cocycle / MC ---------------------------------

pg = Dict(a => QQ(1), b => QQ(2))   # generic fiber
ps = Dict(a => QQ(1), b => QQ(1))   # special fiber on the wall
ig = findfirst(s -> s.sigma == S.generic.sigma, S.strata)
is = findfirst(s -> s.sigma == J.target.sigma, S.strata)

println("-"^64)
println("2. explain_jumps! — deformation cause of the jump")
println("-"^64)
println("  generic point  (a,b) = (1,2)")
println("  special point  (a,b) = (1,1)")

explain_jumps!(S; order = 2, points = Dict(ig => pg, is => ps))
cause = J.cause
show(stdout, MIME("text/plain"), cause)
println()

# --- 3. Isomorphism of the two fibers ---------------------------------------

println("-"^64)
println("3. isomorphism — are the fibers the same GL(n) orbit?")
println("-"^64)
Lg = cause.fiber_generic
Ls = cause.fiber_special
iso = cause.iso
show(stdout, MIME("text/plain"), iso)
println()

if iso.isomorphic
    println("Found P ⇒ isomorphic degeneration / coordinate artefact")
    println("  (aligns with a gauge_like wall class when the jump is only apparent)")
else
    println("No P (here: invariant obstruction) ⇒ genuine change of isomorphism class")
    println("  Jump center 0→1 already separates the fibers; wall class is")
    println("  nontrivial + MC-integrable ⇒ local moduli explain the transition.")
end

println()
println("Narrative:")
println("  Jump: center 0→1 (stratify)")
println("  → cause: $(cause.verdict) (wall=$(cause.wall_class), MC order 2)")
println("  → isomorphism: isomorphic=$(iso.isomorphic), reason=$(iso.reason)")

println()
println("="^64)
println("Done.  Path: stratify → explain_jumps! → isomorphism")
println("  docs: docs/README.md, docs/parametric.md")
println("="^64)
