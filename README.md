# ParametricLie.jl

Exact symbolic toolkit for **finite-dimensional Lie algebras** given by structure
constants: classical invariants, parameterized families with certified strata,
Chevalley–Eilenberg cohomology, central extensions, and formal deformations.

---

## Motivation

Many Lie-algebra questions are algebraic and exact:

- What is the structure of this algebra (center, radical, Levi, …)?
- In a family `L(a)`, where do invariants jump?
- Is a jump a genuine change of isomorphism class, or only a change of coordinates?
- Which H² classes extend, and which formal deformations are obstructed?

ParametricLie answers these with **certificates** (bases, cocycles, pivot trails,
Maurer–Cartan caches)—not dimensions alone—over AbstractAlgebra/Nemo fields
(ℚ, 𝔽ₚ, polynomial and rational-function rings).

---

## Features

| Area | Capability |
|------|------------|
| **Single algebra** | Exact Jacobi check; compute center, derived/lower-central series, Killing form, radical, Levi decomposition, and derivations, with bases and matrices as certificates |
| **Parameter families** | Structure constants over polynomial / rational rings; branch safely when pivots may vanish; partition parameter space into certified regions where classical invariants are constant |
| **Phase / jump analysis** | Report where invariants change across regions; optionally explain a jump via the wall class in H², truncated Maurer–Cartan integrability, and whether two fibers are isomorphic |
| **Cohomology** | Chevalley–Eilenberg cohomology (lazy evaluation); trivial and adjoint coefficients; construct central extensions from 2-cocycles |
| **Deformations** | Infinitesimal and truncated formal deformations of the bracket (Nijenhuis–Richardson / Maurer–Cartan); compare deformations up to gauge; optional gauge normal form |
| **Isomorphism helpers** | Obstruct non-isomorphism via invariants; search for a change-of-basis matrix in small dimension (not a full classification) |

**Out of scope:** general CAS / own Gröbner engine, numerical Lie-group ops, root
systems, complete isomorphism classification of all Lie algebras.

---

## LLM disclaimer

Parts of this repository were drafted with assistance from large language models.
Every design choice and implementation was reviewed by a human maintainer.
Correctness is checked against a test suite whose cases are taken from classical
Lie-algebra textbooks and standard examples (Heisenberg, sl₂, known Betti numbers,
Jacobi identities, etc.).

---

## Installation

```julia
# from the repo
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Requirements: Julia ≥ 1.10; `AbstractAlgebra`, `Nemo`.

Coefficient rings are AA/Nemo parents, e.g. `Nemo.QQ`, `Nemo.GF(p)`, or
`fraction_field(polynomial_ring(...))`.

---

## Quick start

```julia
using ParametricLie
import AbstractAlgebra, Nemo

F = Nemo.QQ
H = LieAlgebra(F, 3, Dict((1, 2) => [0, 0, 1]))  # Heisenberg
r = analyze(H)

# Parameter family → stratification
R, (a, b) = AbstractAlgebra.polynomial_ring(F, [:a, :b])
L = LieAlgebra(R, 3, Dict(
    (1, 2) => [R(0), a, R(0)],
    (1, 3) => [R(0), R(0), a - b],
))
S = stratify(L; invariants = [:center, :derived_dim, :solvability, :nilpotency])

# Adjoint H² and a truncated deformation
H2 = cohomology(H, adjoint_module(H), 2)
D  = formal_deformation(H, basis_matrix(H2)[:, 1]; order = 2)
```

End-to-end narrative (jump → deformation cause → fiber isomorphism):

```bash
julia --project=. example/jump_cause_iso.jl
```

---

## Documentation

| Document | Contents |
|----------|----------|
| [docs/README.md](docs/README.md) | Documentation map |
| [docs/lie.md](docs/lie.md) | Single-algebra core (`src/lie/`) |
| [docs/parametric.md](docs/parametric.md) | Families, stratification, jump causes |
| [docs/cohomology.md](docs/cohomology.md) | CE cohomology & central extensions |
| [docs/deformations.md](docs/deformations.md) | NR bracket, Maurer–Cartan, gauge |
| [example/README.md](example/README.md) | Runnable scripts |

---

## Package layout

```text
src/
  ParametricLie.jl
  lie/              # structure constants, invariants, analyze, isomorphism
  parametric/       # specialize, CondTree, stratify, jump_explain
  cohomology/       # exterior algebra, modules, CE complex, central extensions
  deformations/     # NR bracket, Maurer–Cartan, gauge
docs/               # module references (start at docs/README.md)
example/            # demos
test/unit/          # mirrored by layer: lie/, parametric/, cohomology/, deformations/
```

---

## Examples

```bash
julia --project=. example/heisenberg.jl
julia --project=. example/parametric_family.jl
julia --project=. example/stratify_ab_family.jl
julia --project=. example/central_extension.jl
julia --project=. example/heisenberg_adjoint_cohomology.jl
julia --project=. example/jump_cause_iso.jl
```

See [example/README.md](example/README.md) for what each script demonstrates.

---

## Testing

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Continuous integration runs on pushes and pull requests
([`.github/workflows/ci.yml`](.github/workflows/ci.yml)).

---

## Citation

If you use ParametricLie in academic work, please cite it. A machine-readable
record is in [`CITATION.cff`](CITATION.cff). APA-style:

> Liu, J. (2026). *ParametricLie.jl* (Version 0.1.0) [Computer software].
> https://github.com/LiuJiazheng/ParametricLie

---

## Future work

| Direction | Intent |
|-----------|--------|
| Gauge theory | Residual H¹ / Aut in normal forms; gauge beyond finite truncation; NF-based ID of MC solutions |
| Parametric cohomology | Conditional CE: certify jumps of dim Hᵏ(L(a), M) under assumption sets Σ |
| MC → families | Realize truncated deformations as polynomial bracket families, then `stratify` |
| Isomorphism of fibers | Stronger search for `P` in GL(`n`) (adapted flags); not a full classification |

**Non-goals:** own Gröbner engine; numerical `exp`/`log` on SO(n)/SE(n); root-system
representation theory; exhaustive iso classification.

---

## Related packages

| Package | Role vs ParametricLie |
|---------|------------------------|
| [AbstractAlgebra.jl](https://nemocas.github.io/AbstractAlgebra.jl/) / [Nemo.jl](https://nemocas.github.io/Nemo.jl/) | Exact fields, polynomials, linear algebra (our scalar backend) |
| [LieGroups.jl](https://github.com/JuliaManifolds/LieGroups.jl) | Concrete groups and numerics; different starting point (optional later bridge) |
| Oscar / GAP / Sage Lie algebras | Broader CAS; we focus on parametric strata + CE witnesses + deformations |

---

## License

[MIT](LICENSE) © 2026 Jiazheng Liu
