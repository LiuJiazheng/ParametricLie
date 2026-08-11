# ParametricLie.jl

Exact symbolic toolkit for **finite-dimensional Lie algebras** defined by structure constants — parameterized families, structural stratification, cohomology, extensions, and deformations.

> Complements [LieGroups.jl](https://github.com/JuliaManifolds/LieGroups.jl) (smooth geometry / numerics on concrete groups).  
> **Not** a competitor for `exp`/`log`/`hat`/`vee` on SO(3)/SE(3).

| Doc | Contents |
|-----|----------|
| [docs/POSITIONING.md](docs/POSITIONING.md) | Ecosystem stance and brand |
| [docs/lie.md](docs/lie.md) | Single-algebra core — API & algorithms (`src/lie/`) |
| [docs/parametric.md](docs/parametric.md) | Families, conditional LA, stratification (`src/parametric/`) |
| [docs/cohomology.md](docs/cohomology.md) | CE cohomology — what you can do, API, examples (`src/cohomology/`) |
| [ParametricLie_Project_Overview.md](ParametricLie_Project_Overview.md) | Roadmap |

**Brand:** Exact · Parametric · Cohomological · Structure-aware

## Status

- **`lie`** — structure-constant `LieAlgebra` over AA/Nemo fields, Jacobi certificates, classical invariants, `analyze`.
- **`parametric`** — polynomial/rational families, assumption-aware linear algebra, conditional analysis trees, `stratify` + fiber validation.
- **`cohomology`** — exterior algebra \(\Lambda^\bullet V\), coefficient modules, lazy CE \(H^\bullet\).

### Scalar fields

Coefficient domains are [AbstractAlgebra.jl](https://nemocas.github.io/AbstractAlgebra.jl/stable/field_interface/) parents (`Field` / polynomial / fraction rings):

| Domain | Parent (examples) |
|--------|-------------------|
| ℚ | `Nemo.QQ` |
| 𝔽ₚ | `Nemo.GF(p)` |
| rational functions | `fraction_field(polynomial_ring(...))` |

### Ecosystem stance

| Layer | Responsibility |
|-------|----------------|
| **LieGroups.jl / Manifolds.jl** | Groups, exp/log, tangent geometry, numerical ops |
| **This package** | Structure constants, exact algebra, parameters, stratification |
| **AbstractAlgebra / Nemo** | Fields, polynomials, exact linear algebra |

## Requirements

- Julia ≥ 1.10
- Dependencies: `AbstractAlgebra`, `Nemo` (plus stdlibs `LinearAlgebra`, `SparseArrays`)

## Setup

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Quick start

```julia
using ParametricLie
import AbstractAlgebra, Nemo

F = Nemo.QQ
H = LieAlgebra(F, 3, Dict((1, 2) => [0, 0, 1]))  # Heisenberg
r = analyze(H)

R, (a, b) = AbstractAlgebra.polynomial_ring(F, [:a, :b])
L = LieAlgebra(R, 3, Dict(
    (1, 2) => [R(0), a, R(0)],
    (1, 3) => [R(0), R(0), a - b],
))
S = stratify(L; invariants = [:center, :derived_dim, :solvability, :nilpotency])
```

## Examples

```bash
julia --project=. example/heisenberg.jl
julia --project=. example/parametric_family.jl
julia --project=. example/stratify_ab_family.jl
```

## Test

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Unit tests live under `test/unit/lie/` and `test/unit/parametric/`.

## Layout

```
src/
  ParametricLie.jl
  lie/                 # single-algebra core
  parametric/          # families, CondTree, stratify
docs/
  lie.md
  parametric.md
  POSITIONING.md
test/unit/{lie,parametric}/
example/
```
