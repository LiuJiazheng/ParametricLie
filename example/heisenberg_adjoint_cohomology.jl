# example/heisenberg_adjoint_cohomology.jl
#
# Nontrivial coefficient module: Chevalley–Eilenberg cohomology of the
# 3D Heisenberg algebra with values in the adjoint module.
#
#   H = ⟨e1, e2, e3 | [e1,e2]=e3⟩
#   M = H  (adjoint),   x · y = [x,y]
#
# Classical dictionary:
#   H⁰(H,H) ≅ Z(H)          center
#   H¹(H,H) ≅ Out(H)        outer derivations Der/Inn
#   H²(H,H)                 infinitesimal deformations
#   H³(H,H)                 obstruction space
#
# Run from the package root:
#   julia --project=. example/heisenberg_adjoint_cohomology.jl

using ParametricLie
import Nemo

const F = Nemo.QQ

println("="^60)
println("ParametricLie example — Heisenberg adjoint cohomology")
println("="^60)

# --- algebra & module -------------------------------------------------------

H = LieAlgebra(F, 3, Dict((1, 2) => [0, 0, 1]))
M = adjoint_module(H)

println()
println("H = Heisenberg₃ over QQ,  [e1,e2]=e3")
println("Jacobi: ", check_jacobi(H).ok ? "OK" : "FAIL")
println()
println("Module M = adjoint_module(H),  dim(M) = ", dim(M))
println("dim C^k = dim(M) * binom(3,k):")
for k in 0:3
    println("  C^$k  dim = ", cochain_dim(H, M, k))
end

# --- lazy complex: ask for H³, get H⁰…H³ cached -----------------------------

C = ce_complex(H, M)
println()
println("-"^60)
println("ce_complex(H, M)  — lazy; filled=$(C.filled) before any query")
println("-"^60)

H3 = cohomology(C, 3)   # fills degrees 0,1,2,3 and caches
println()
println("After cohomology(C, 3): filled = ", C.filled, "  (expect 3)")

# --- Betti table ------------------------------------------------------------

println()
println("-"^60)
println("Betti numbers  H^k(H, H)")
println("-"^60)
println()
println(lpad("k", 4), lpad("dim C", 10), lpad("dim Z", 10),
        lpad("dim B", 10), lpad("dim H", 10))
for k in 0:3
    Zk = cocycles(C, k)
    Bk = coboundaries(C, k)
    Hk = cohomology(C, k)
    println(lpad(string(k), 4),
            lpad(string(cochain_dim(C, k)), 10),
            lpad(string(size(Zk, 2)), 10),
            lpad(string(size(Bk, 2)), 10),
            lpad(string(dim(Hk)), 10))
end
println()
println("Expected Betti:  1, 4, 5, 2")

# --- interpret H⁰ and H¹ ----------------------------------------------------

println()
println("-"^60)
println("Interpretation")
println("-"^60)

H0 = cohomology(C, 0)
H1 = cohomology(C, 1)
H2 = cohomology(C, 2)

println()
println("H⁰ ≅ center(H):")
println("  dim H⁰ = ", dim(H0), "   dim center = ", dim(center(H)))

println()
println("H¹ ≅ Out = Der / Inn:")
println("  dim Der  = ", dim(derivations(H)))
println("  dim Inn  = dim(H) − dim(center) = ", dim(H) - dim(center(H)))
println("  dim Out  = ", dim(derivations(H)) - (dim(H) - dim(center(H))))
println("  dim H¹   = ", dim(H1), "   (match)")

println()
println("H² = infinitesimal deformations of the bracket (dim ", dim(H2), ")")
println("H³ = obstruction space (dim ", dim(H3), ")")

# --- cocycle representatives ------------------------------------------------

println()
println("-"^60)
println("Cocycle representatives  (columns of basis_matrix)")
println("-"^60)
println()
println("Coordinates for C^k: lex multi-indices I, then m=3 coords of ω(e_I).")
println("For k=1: blocks (ω(e1), ω(e2), ω(e3)), each in F³  → length 9.")
println("For k=2: blocks on (e1∧e2, e1∧e3, e2∧e3)           → length 9.")
println("For k=3: single block ω(e1,e2,e3)                   → length 3.")

for k in 1:3
    Hk = cohomology(C, k)
    println()
    println("H^$k  (", dim(Hk), " classes):")
    println(basis_matrix(Hk))
end

# --- contrast with trivial coefficients ------------------------------------

println()
println("-"^60)
println("Contrast: trivial module H^k(H, F)")
println("-"^60)
println()
Ct = ce_complex(H)   # trivial
print("Betti trivial: ")
println(join([string(dim(cohomology(Ct, k))) for k in 0:3], ", "))
print("Betti adjoint: ")
println(join([string(dim(cohomology(C, k))) for k in 0:3], ", "))
println()
println("Trivial H² classifies central extensions by F;")
println("adjoint H² classifies deformations of the Lie bracket.")
println()
println("See also: example/central_extension.jl  (trivial H² → Heisenberg)")
println("="^60)
println("Done.")
println("="^60)
