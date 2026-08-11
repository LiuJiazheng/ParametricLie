# Cohomology (`src/cohomology/`)

Compute Chevalley–Eilenberg cohomology of a finite-dimensional Lie algebra over
a field (primarily `Nemo.QQ`): pick a coefficient module, build the complex, ask
for \(H^k\) on demand. Certificates are bases of cocycle representatives — not
dimensions alone.

```julia
using ParametricLie
import Nemo

L = LieAlgebra(Nemo.QQ, 3, Dict((1, 2) => [0, 0, 1]))  # Heisenberg
H2 = cohomology(L, adjoint_module(L), 2)
dim(H2)                 # 5
basis_matrix(H2)        # columns = cocycle reps in C²
```

---

## What you can do

| Goal | Entry point |
|------|-------------|
| Betti numbers / cocycle bases | `cohomology(L, k)` or `cohomology(L, M, k)` |
| Reuse one complex, many degrees | `C = ce_complex(L, M)` then `cohomology(C, k)` |
| Inspect \(Z^k\), \(B^k\), \(d^k\) | `cocycles`, `coboundaries`, `ce_differential` |
| Trivial coeffs \(H^\bullet(\mathfrak{g}, F)\) | default / `trivial_module(L)` |
| Adjoint coeffs \(H^\bullet(\mathfrak{g}, \mathfrak{g})\) | `adjoint_module(L)` |
| Custom module | `LieModule(L, m, action_matrices)` |
| Central extension from a 2-cocycle | `central_extension(L, ω)` |
| Graded forms / wedge product | `exterior_algebra`, `wedge`, `form_eval` |

**Lazy evaluation.** `ce_complex` builds nothing. Requesting degree \(k\) fills and
**caches** \(d^j, Z^j, B^j, H^j\) for all \(j \le k\) (same idea as `LieSeries`).
Asking for \(H^3\) therefore also computes \(H^0,H^1,H^2\).

**Finite.** \(\dim C^k = m\binom{n}{k}\), and \(C^k = 0\) for \(k > n\). Exact
linear algebra over the coefficient field; always terminates.

---

## Quick start

### 1. Trivial coefficients (scalar cohomology)

```julia
L  = LieAlgebra(Nemo.QQ, 3, Dict((1, 2) => [0, 0, 1]))
C  = ce_complex(L)                    # = trivial_module
for k in 0:3
    Hk = cohomology(C, k)
    println("H^$k dim = ", dim(Hk))
end
# Heisenberg: 1, 2, 2, 1
```

Interpretation highlights:

- \(H^0(\mathfrak{g}, F) \cong F\) always (constants), for connected picture
- \(H^1(\mathfrak{g}, F) \cong (\mathfrak{g}/[\mathfrak{g},\mathfrak{g}])^*\)
- \(H^2(\mathfrak{g}, F)\) classifies **central extensions** by \(F\)

### 2. Adjoint coefficients (outer derivations / deformations)

```julia
M  = adjoint_module(L)
C  = ce_complex(L, M)
H0 = cohomology(C, 0)   # ≅ center
H1 = cohomology(C, 1)   # ≅ Out = Der / Inn
H2 = cohomology(C, 2)   # infinitesimal deformations
```

Runnable walkthrough for Heisenberg:
[`example/heisenberg_adjoint_cohomology.jl`](../example/heisenberg_adjoint_cohomology.jl).

### 3. Central extension from \(H^2(-, F)\)

```julia
A  = LieAlgebra(Nemo.QQ, 2)           # abelian plane
ω  = basis_matrix(cohomology(A, 2))[:, 1]
Ê  = central_extension(A, ω)          # ≅ Heisenberg
is_trivial_cocycle(A, ω)              # false → nonsplit
```

See [`example/central_extension.jl`](../example/central_extension.jl).

---

## API reference

### Coefficient modules

| API | Role |
|-----|------|
| `LieModule` | Left \(L\)-module on \(F^m\) via matrices \(\rho(e_i)\) |
| `trivial_module(L)` | \(m=1\), all actions zero |
| `adjoint_module(L)` | \(M = L\), \(x\cdot y = [x,y]\) |
| `LieModule(L, m, actions)` | Custom; `actions` length \(n\), each \(m\times m\) |
| `dim(M)`, `action_matrices(M)`, `act(M, i, v)` | Inspect / apply \(\rho(e_i)\) |

### CE complex and \(H^k\)

| API | Role |
|-----|------|
| `ce_complex(L)` / `ce_complex(L, M)` | Lazy complex (no matrices yet) |
| `cohomology(C, k)` | \(H^k = Z^k/B^k\) as `CohomologyGroup` |
| `cohomology(L, k)` | Same, trivial module |
| `cohomology(L, M, k)` | Same, given module |
| `dim(Hk)`, `basis_matrix(Hk)` | Dimension and cocycle representatives |
| `Hk.degree` | The index \(k\) |
| `cocycles(C, k)` | Basis matrix of \(Z^k = \ker d^k\) |
| `coboundaries(C, k)` | Basis matrix of \(B^k = \mathrm{im}\, d^{k-1}\) |
| `ce_differential(C, k)` | Matrix of \(d^k : C^k \to C^{k+1}\) |
| `cochain_dim(C, k)` | \(\dim C^k = m\binom{n}{k}\) |
| `coefficient_module(C)` | The module used by `C` |

### Central extensions (trivial-module \(H^2\) application)

| API | Role |
|-----|------|
| `central_extension(L, ω)` | \(\hat L = L \oplus F\cdot z\) with \([x,y] += \omega(x,y)\, z\) |
| `is_trivial_cocycle(L, ω)` | \([\omega] = 0\) in \(H^2(L,F)\) (extension splits) |

`ω` is a length-\(\binom{n}{2}\) coordinate vector (same layout as a column of
`basis_matrix(cohomology(L, 2))`). Requires \(d\omega = 0\).

### Exterior algebra (forms / graded product)

Used as the graded language behind cochains when \(m=1\); also useful on its own.

| API | Role |
|-----|------|
| `exterior_algebra(F, n)` / `exterior_algebra(L)` | \(\Lambda^\bullet V\) |
| `exterior_generator(Λ, i)`, `exterior_monomial(Λ, I)` | Basis elements |
| `wedge` / `*` | Graded-commutative product |
| `homogeneous_part`, `exterior_degree`, `support_degrees` | Degree helpers |
| `interior_product(v, α)`, `form_eval(α, vs)` | Contraction / evaluation |
| `ambient_dim(Λ)`, `dim(Λ)`, `dim(Λ, k)` | \(n\), \(2^n\), \(\binom{n}{k}\) |
| `multi_indices(n, k)`, `coord_index(n, I)` | Lex indexing |

---

## Coordinates (when reading matrices)

Cochain space \(C^k(\mathfrak{g}, M) \cong \mathrm{Hom}(\Lambda^k\mathfrak{g}, M)\).

1. List increasing multi-indices \(I = (i_1 < \cdots < i_k)\) in **lex** order
   (`multi_indices(n, k)`).
2. For each \(I\), store \(m = \dim M\) coordinates of \(\omega(e_{i_1},\ldots,e_{i_k})\).

So \(\dim C^k = m\binom{n}{k}\). When \(m=1\), this is exactly the degree-\(k\)
slice of `ExteriorAlgebra`. Columns of `basis_matrix(Hk)` live in these
coordinates.

Differential (module-induced):

```
(dω)(x₀,…,xₖ) = ∑ᵢ (-1)ⁱ xᵢ · ω(…ˆxᵢ…)
              + ∑_{i<j} (-1)^{i+j} ω([xᵢ,xⱼ], …ˆxᵢ…ˆxⱼ…)
```

---

## Classical dictionary

| Object | Cohomology |
|--------|------------|
| Center | \(H^0(\mathfrak{g}, \mathfrak{g})\) |
| Outer derivations | \(H^1(\mathfrak{g}, \mathfrak{g})\) |
| Central extensions by \(F\) | \(H^2(\mathfrak{g}, F)\) |
| Infinitesimal deformations | \(H^2(\mathfrak{g}, \mathfrak{g})\) |
| Obstructions (next order) | \(H^3(\mathfrak{g}, \mathfrak{g})\) |

Known smoke values used in tests:

| Algebra | Module | Betti \((H^0,\ldots)\) |
|---------|--------|-------------------------|
| Heisenberg₃ | trivial | \(1,2,2,1\) |
| Heisenberg₃ | adjoint | \(1,4,5,2\) |
| \(\mathfrak{sl}_2\) | trivial | \(1,0,0,1\) |
| \(\mathfrak{sl}_2\) | adjoint | \(0,0,0,\ldots\) (rigid, centerless) |

---

## Source layout

```
src/cohomology/
  exterior.jl            # Λ^• V
  module.jl              # LieModule, trivial / adjoint
  ce.jl                  # d, lazy Z / B / H
  central_extension.jl   # H²(—, F) → central extension
```

Field coefficients only (`FieldElem`) for the CE layer.
