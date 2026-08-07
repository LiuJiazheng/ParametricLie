# example/heisenberg.jl
#
# First ParametricLie use case: the 3-dimensional Heisenberg algebra over ℚ.
#
#   [e₁, e₂] = e₃,   all other brackets zero.
#
# Run from the package root:
#   julia --project=. example/heisenberg.jl

using ParametricLie
import Nemo

const F = Nemo.QQ

println("="^60)
println("ParametricLie example — Heisenberg algebra")
println("="^60)

# --- construct --------------------------------------------------------------

# Basis e1, e2, e3 with sole nonzero relation [e1, e2] = e3.
H = LieAlgebra(F, 3, Dict((1, 2) => [0, 0, 1]))

println()
println("Lie algebra H: dim=$(dim(H)) over QQ")
println("Structure: [e1, e2] = e3  (Heisenberg)")
println()
println("Jacobi check: ", check_jacobi(H).ok ? "OK" : "FAIL")

# --- summary report ---------------------------------------------------------

println()
println("-"^60)
println("analyze(H)  — structural summary")
println("-"^60)
println()

r = analyze(H)
show(stdout, MIME("text/plain"), r)
println()

# --- detail: series ---------------------------------------------------------

println()
println("-"^60)
println("Detail: derived & lower-central series (cached on report)")
println("-"^60)

D = derived_series(r)
C = lower_central_series(r)

println()
println("derived_series dims:       ", join(dim.(terms(D)), " → "))
println("lower_central_series dims: ", join(dim.(terms(C)), " → "))
println("layers(derived) = ", layers(D), "   (index of first zero term)")
println("layers(lower_central) = ", layers(C))
println()
println("solvable:  ", is_solvable(r))
println("nilpotent: ", is_nilpotent(r))

println()
println("derived algebra [H,H] basis (coords):")
for v in basis_elems(derived_algebra(H))
    println("  ", v.coords)
end

# --- detail: center & radical -----------------------------------------------

println()
println("-"^60)
println("Detail: center & solvable radical")
println("-"^60)

Z = center(r)
R = radical(r)
println()
println("center dim=$(dim(Z)); basis coords:")
for z in basis_elems(Z)
    println("  ", z.coords)
end
println()
println("radical dim=$(dim(R))  (here radical = H, since H is solvable)")
println("is_ideal(H, radical)? ", is_ideal(H, R))

# --- detail: Levi -----------------------------------------------------------

println()
println("-"^60)
println("Detail: Levi decomposition (cached)")
println("-"^60)

lev = levi_decomposition(r)
println()
println(lev)
println("levi_kind(r)          = ", levi_kind(r))
println("simple_factor_dims(r) = ", simple_factor_dims(r))
println("dim(Levi subalgebra)  = ", dim(lev.levi))
println("dim(radical)          = ", dim(lev.radical))
println("dim(quotient H/rad)   = ", dim(lev.quotient))
println()
println("Interpretation: H is nilpotent ⇒ Levi factor is {0};")
println("the whole algebra is the radical. No simple factors.")

# --- detail: Killing & Der --------------------------------------------------

println()
println("-"^60)
println("Detail: Killing form & Der(H)")
println("-"^60)

K = killing_form(r)
println()
println("Killing matrix:")
println(K)
println("rank = ", killing_rank(r))

Der = derivations(r)
println()
println(Der)
println("dim Der(H) = ", dim(Der), "  (classical: 6 for Heisenberg₃)")
println()
println("First derivation matrix (one basis element of Der):")
println(first(basis_matrices(Der)))

# --- optional: inspect a bracket --------------------------------------------

println()
println("-"^60)
println("Sanity: bracket from structure constants")
println("-"^60)
e1, e2, e3 = basis_elem(H, 1), basis_elem(H, 2), basis_elem(H, 3)
println()
println("[e1, e2] = ", lie_bracket(H, e1, e2).coords, "  (expect e3)")
println("[e1, e3] = ", lie_bracket(H, e1, e3).coords, "  (expect 0)")
println("[e2, e3] = ", lie_bracket(H, e2, e3).coords, "  (expect 0)")

println()
println("="^60)
println("Done. For bases of Levi/Der on other algebras, call the same")
println("detail APIs on analyze(L) — matrices depend on the chosen basis.")
println("="^60)
