# Deformations (`src/deformations/`)

Formal deformations of a Lie bracket in the adjoint Chevalley–Eilenberg complex:
Nijenhuis–Richardson bracket, truncated Maurer–Cartan integration, and gauge
equivalence. Builds on [cohomology.md](cohomology.md) (same \(C^\bullet(\mathfrak{g},\mathfrak{g})\)
coordinates).

```julia
using ParametricLieAlgebras
import Nemo

L  = LieAlgebra(Nemo.QQ, 3, Dict((1, 2) => [0, 0, 1]))
H2 = cohomology(L, adjoint_module(L), 2)
φ1 = [basis_matrix(H2)[i, 1] for i in 1:size(basis_matrix(H2), 1)]
D  = formal_deformation(L, φ1; order = 2)
```

---

## Nijenhuis–Richardson bracket

On adjoint cochains \(C^p(\mathfrak{g},\mathfrak{g})\):

```julia
μ  = adjoint_bracket_cochain(L)     # bracket as μ ∈ C²
φ  = basis_matrix(cohomology(L, adjoint_module(L), 2))[:, 1]
br = nr_bracket(L, φ, 2, φ, 2)      # [φ,φ]_NR ∈ C³
```

| API | Role |
|-----|------|
| `adjoint_bracket_cochain(L)` | \(\mu\in C^2\) from structure constants |
| `nr_circle(L, φ, p, ψ, q)` | \(\varphi\circ\psi\in C^{p+q-1}\) |
| `nr_bracket(L, φ, p, ψ, q)` | \([\varphi,\psi]_{\mathrm{NR}}\) |

With our CE sign convention,

\[
d\varphi = (-1)^{p+1}\,[\mu,\varphi]_{\mathrm{NR}}.
\]

---

## Truncated Maurer–Cartan

Integrate \(\phi=\sum_{k\ge 1}\varepsilon^k\phi_k\) order by order up to a truncation
`order` (default **2**). Convention (char \(\neq 2\)):

\[
d\phi = \tfrac12[\phi,\phi]_{\mathrm{NR}}.
\]

```julia
D = formal_deformation(L, φ1; order = 2)
extend!(D, 4)
deformation_term(D, 2)     # φ₂ ∈ C²
mc_rhs(D, 2)               # ψ₂ ∈ C³
mc_certificate(D)
```

| API | Role |
|-----|------|
| `formal_deformation(L, φ1; order=N)` | MC along a 2-cocycle |
| `formal_deformation(L; order=N)` | one problem per \(H^2(L,L)\) basis column |
| `extend!(D, N)` | raise truncation; reuse cache |
| `deformation_term(D, k)` / `mc_rhs(D, k)` | cached terms |
| `obstruction_cochain(D)` | failing \(\psi_k\), or `nothing` |
| `is_rigid(L)` | \(\dim H^2(L,L)=0\) |

Always terminates: at most `order` linear solves.

---

## Gauge equivalence and normal form

Same-order comparison of two [`FormalDeformation`](@ref) objects:

| API | Role |
|-----|------|
| `equivalent(D1, D2; order=N)` | gauge-equivalent through order \(N\)? |
| `equivalent_with_gauge(...)` | same, plus witness \(\alpha_k\in C^1\) |
| `gauge_normal_form(D; order=N)` | section: each \(\phi_k\) in a fixed complement of \(B^2\) |
| `gauge_normal_form!` / `_with_gauge` | in-place / with realizing \(\alpha_k\) |

Order 1 reduces to \(\phi'_1-\phi_1\in B^2\). Higher orders solve recursively with
NR corrections. Normal form uses a fixed splitting \(C^2=B^2\oplus W\) (greedy on
the standard basis)—a section of each gauge orbit, not a classification engine.

---

## Source layout

```text
src/deformations/
  nr_bracket.jl
  maurer_cartan.jl
  gauge.jl
```

See also: [cohomology.md](cohomology.md), [parametric.md](parametric.md)
(jump causes use wall cocycles + MC + fiber isomorphism).
