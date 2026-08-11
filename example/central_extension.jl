# example/central_extension.jl
#
# Nontrivial H² application: the Heisenberg algebra as a central extension of
# the abelian plane ℚ² by a generator of H²(ℚ², ℚ) ≅ ℚ.
#
#   0 → ℚ → Heisenberg₃ → ℚ² → 0
#
# Run from the package root:
#   julia --project=. example/central_extension.jl

using ParametricLie
import AbstractAlgebra
import Nemo

const F = Nemo.QQ

println("="^60)
println("ParametricLie example — central extension from H²")
println("="^60)

# --- base algebra: abelian ℚ² ----------------------------------------------

A = LieAlgebra(F, 2)
println()
println("Base algebra A = ℚ² (abelian), dim=$(dim(A))")
println("Jacobi: ", check_jacobi(A).ok ? "OK" : "FAIL")

# --- classify central extensions: H²(A, F) ---------------------------------

H2 = cohomology(A, 2)
println()
println("-"^60)
println("H²(A, F) classifies central extensions of A by F")
println("-"^60)
println()
println("dim H²(A, F) = ", dim(H2), "   (expect 1 = dim Λ²(ℚ²)*)")
println("cocycle representative ω (coords on [e1∧e2]):")
ω = [basis_matrix(H2)[i, 1] for i in 1:binomial(2, 2)]
println("  ω = ", ω)
println("trivial class? ", is_trivial_cocycle(A, ω), "   (expect false)")

# --- build the extension ----------------------------------------------------

Ê = central_extension(A, ω)
println()
println("-"^60)
println("Central extension Ê = A ⊕ F·z  with  [e1,e2] = ω·z")
println("-"^60)
println()
println("dim(Ê) = ", dim(Ê))
println("Jacobi: ", check_jacobi(Ê).ok ? "OK" : "FAIL")
println()
println("[e1, e2] = ", lie_bracket(Ê, basis_elem(Ê, 1), basis_elem(Ê, 2)).coords)
println("[e1, e3] = ", lie_bracket(Ê, basis_elem(Ê, 1), basis_elem(Ê, 3)).coords)
println("[e2, e3] = ", lie_bracket(Ê, basis_elem(Ê, 2), basis_elem(Ê, 3)).coords)

r = analyze(Ê)
println()
show(stdout, MIME("text/plain"), r)
println()

println()
println("Interpretation: Ê is the 3D Heisenberg algebra (up to scaling z).")
println("Nonsplit because [ω] ≠ 0: derived algebra = center ≠ 0.")

# --- contrast: split extension ---------------------------------------------

println()
println("-"^60)
println("Split extension from the zero cocycle")
println("-"^60)
E0 = central_extension(A, [F(0)])
println()
println("dim(E0) = ", dim(E0), "  derived dim = ", dim(derived_algebra(E0)),
        "  center dim = ", dim(center(E0)))
println("⇒ abelian A ⊕ F (direct product).")

# --- sl₂ has no nontrivial central extensions ------------------------------

println()
println("-"^60)
println("sl₂: H² = 0 ⇒ only the split extension exists")
println("-"^60)
sl2 = LieAlgebra(F, 3, Dict(
    (1, 2) => [0, 2, 0],
    (1, 3) => [0, 0, -2],
    (2, 3) => [1, 0, 0],
))
println()
println("dim H²(sl₂, F) = ", dim(cohomology(sl2, 2)))
println("="^60)
println("Done.")
println("="^60)
