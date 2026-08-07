# example/stratify_ab_family.jl
#
# End-to-end v0.2.1 → v0.2.6 story for the 3D family
#   [e1,e2] = a e2,  [e1,e3] = (a-b) e3  over QQ[a,b].
#
# Run from the package root:
#   julia --project=. example/stratify_ab_family.jl

using ParametricLie
import AbstractAlgebra
import Nemo

const QQ = Nemo.QQ

println("="^60)
println("ParametricLie — 3D family stratification (a, a−b)")
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
println()
println("v0.2.1 Jacobi identity: ", check_jacobi(L).ok)

println()
println("-"^60)
println("v0.2.2 analyze_generic")
println("-"^60)
rg = analyze_generic(L)
println("  center_dim  = ", dim(center(rg)))
println("  derived_dim = ", dim(derived_algebra(generic_algebra(L))))
println("  solvable    = ", is_solvable(rg))
println("  nilpotent   = ", is_nilpotent(rg))

suite = [:center, :derived_dim, :solvability, :nilpotency]
println()
println("-"^60)
println("v0.2.4 analyze_conditional  suite=$suite")
println("-"^60)
T = analyze_conditional(L; invariants = suite)
show(stdout, MIME("text/plain"), T)
println()
println("complete? ", is_complete(T; invariants = suite))

println()
println("-"^60)
println("v0.2.5 stratify + jump_table")
println("-"^60)
S = stratify(T; invariants = suite)
show(stdout, MIME("text/plain"), S)
println()
println("Confirmed jumps:")
for j in jump_table(S)
    show(stdout, MIME("text/plain"), j)
    println()
end

println()
println("-"^60)
println("v0.2.6 validate fibers vs v0.1 analyze")
println("-"^60)
checks = [
    ("generic (1,2)", S.generic, Dict(a => QQ(1), b => QQ(2))),
]
for st in S.strata
    cd = get(st.signature, :center_dim, nothing)
    if cd == 3
        push!(checks, ("abelian (0,0)", st, Dict(a => QQ(0), b => QQ(0))))
    elseif cd == 1 && status(st.sigma, a) === PIVOT_ZERO
        push!(checks, ("a=0,b=1", st, Dict(a => QQ(0), b => QQ(1))))
    elseif cd == 1 && status(st.sigma, a - b) === PIVOT_ZERO
        push!(checks, ("a=b=1", st, Dict(a => QQ(1), b => QQ(1))))
    end
end

for (label, st, pt) in checks
    st === nothing && continue
    v = validate_stratum(L, st, pt)
    println("  $label  ok=$(v.ok)  mismatches=$(v.mismatches)")
end

println()
println("="^60)
println("Done.  docs: docs/parametric.md")
println("="^60)
