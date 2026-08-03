# ParametricLie.jl

Exact symbolic toolkit for **finite-dimensional Lie algebras** defined by structure constants — parameterized families, structural stratification, cohomology, extensions, and deformations.

> Complements [LieGroups.jl](https://github.com/JuliaManifolds/LieGroups.jl) (smooth geometry / numerics on concrete groups).  
> **Not** a competitor for `exp`/`log`/`hat`/`vee` on SO(3)/SE(3).

Full positioning and locked decisions: [docs/POSITIONING.md](docs/POSITIONING.md)  
Roadmap & testing philosophy: [ParametricLie_Project_Overview.md](ParametricLie_Project_Overview.md)

**Brand:** Exact · Parametric · Cohomological · Structure-aware

## Status

**v0.1 (in progress)** — **symbolic core we own**: structure-constant `LieAlgebra` over an AA/Nemo field (tests on `Nemo.QQ`), `lie_bracket`, Jacobi certificates, classical invariants, `analyze`. Parametric / cohomology / LieGroups adapters come after.

### Scalar fields

We do **not** define a custom Field trait. Coefficient domains are [AbstractAlgebra.jl](https://nemocas.github.io/AbstractAlgebra.jl/stable/field_interface/) parents (`Field`) with elements `FieldElem`:

| Domain | Parent (examples) |
|--------|-------------------|
| ℚ | `Nemo.QQ` |
| 𝔽ₚ | `Nemo.GF(p)` |
| ℝ / ℂ (arb) | `Nemo.RealField()`, `Nemo.ComplexField()` |
| rational functions | `fraction_field(polynomial_ring(...))` |

### Ecosystem stance

| Layer | Responsibility |
|-------|----------------|
| **LieGroups.jl / Manifolds.jl** | Groups, exp/log, tangent geometry, numerical ops |
| **This package** | Structure constants, exact algebra, parameters, cohomology, deformations |
| **AbstractAlgebra / Nemo (/ Oscar)** | Fields, polynomials, exact linear algebra |
| **Adapters (later, soft deps)** | Bridge concrete `LieAlgebra(G)` ↔ our view; specialize symbols → numerics |

Shared concept names follow LieGroups where possible (`lie_bracket`). Structure-constant brackets are still implemented here; concrete-group kernels are not.

## Requirements

- Julia ≥ 1.10 ([juliaup](https://github.com/JuliaLang/juliaup) recommended)
- Dependencies: `AbstractAlgebra`, `Nemo` (plus stdlibs `LinearAlgebra`, `SparseArrays`)

## Setup

```bash
cd /path/to/LieGroupJulia
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Develop

```bash
julia --project=.
```

```julia
using ParametricLie
import Nemo   # use `import` to avoid clashing exports (dim, …)

L = LieAlgebra(Nemo.QQ, 3)
```

## Test

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

| Path | Layer |
|------|--------|
| `test/unit/` | small components |
| `test/identities/` | Jacobi and other identities |
| `test/examples/` | known algebras + fixtures |
| `test/property/` | random basis changes, etc. |
| `test/differential/` | compare with Sage/GAP/Oscar |
| `fixtures/` | expected invariants / structure data |

## Layout

```
src/
  ParametricLie.jl
  types.jl
  bracket.jl       # lie_bracket (+ bracket alias)
  invariants.jl
  analyze.jl
docs/POSITIONING.md
test/
fixtures/
```
