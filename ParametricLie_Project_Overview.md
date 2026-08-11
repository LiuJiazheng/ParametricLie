# ParametricLie.jl --- Project Vision & Design

> **Mission:** Build a Julia-first symbolic computation library for
> finite-dimensional Lie algebras, emphasizing exact algebra,
> parameterized families, cohomology, and mathematically verifiable
> computation.

> **Differentiation (locked):** complements LieGroups.jl — we start from
> structure constants and symbolic structure space; they start from concrete
> smooth groups and numerics. See [docs/POSITIONING.md](docs/POSITIONING.md).

------------------------------------------------------------------------

# 1. Vision

ParametricLie.jl is **not** intended to be another computer algebra
system.

Instead, it aims to become a **research-oriented Lie algebra toolkit**
capable of answering questions such as:

-   What is the structure of this Lie algebra?
-   What changes when parameters vary?
-   Where do symmetries appear?
-   Which deformations are genuine?
-   Which extensions are nontrivial?
-   Can these results be certified mathematically?

The long-term goal is to bridge

> **Abstract Lie theory → Symbolic algorithms → Research exploration.**

------------------------------------------------------------------------

# 2. Why Julia?

Julia provides an ideal balance between mathematical expressiveness and
performance.

Rather than reimplementing everything, the project will build upon the
existing ecosystem.

## Core dependencies

-   LinearAlgebra
-   SparseArrays
-   AbstractAlgebra.jl
-   Nemo.jl
-   Symbolics.jl
-   Test

Future integrations:

-   DifferentiationInterface.jl
-   ChainRulesCore.jl
-   ForwardDiff.jl
-   SciML ecosystem

Early versions focus on **exact symbolic computation**, not numerical
approximation.

------------------------------------------------------------------------

# 3. Scope

## Included

-   Finite-dimensional Lie algebras
-   Exact rational arithmetic
-   Polynomial/rational-function parameters
-   Structure constants
-   Structural invariants
-   Cohomology
-   Central extensions
-   Infinitesimal deformations
-   Parameter stratification

## Explicitly out of scope (initially)

-   General CAS
-   Gröbner basis implementation
-   Numerical Lie-group optimization
-   Root systems & representation theory
-   Infinite-dimensional Lie algebras
-   Complete isomorphism classification

------------------------------------------------------------------------

# 4. Roadmap

## v0.1 --- Symbolic Core (owned by us)

**Goal (one sentence):**

> Own a reliable **structure-constant Lie algebra** over an AbstractAlgebra/Nemo
> field (primary tests on \(\mathbb{Q}\)): verify Jacobi, compute classical
> invariants, and return **certificates** — not dimensions alone.

This is the layer we must implement ourselves. LieGroups.jl cannot supply
arbitrary \(c_{ij}^k\) over \(\mathbb{Q}\) / \(\mathbb{F}_p\); adapters come later.

### In scope

-   `LieAlgebra{C}` + `coefficient_ring` (AA/Nemo `Field`, not a custom trait)
-   Structure constants (dense ok; sparse optional)
-   `lie_bracket` / `lie_bracket!` (structure-constant; name aligned with LieGroups)
-   `check_jacobi` → witness / residual certificate
-   Change of basis (and transformed structure constants)
-   Exact linear algebra over the coefficient field (via AA/Nemo matrices)
-   Adjoint representation (needed by Killing / many invariants)
-   Center (with basis)
-   Derived series / lower central series
-   Solvability / nilpotency
-   Killing form (matrix + rank / radical hooks)
-   Derivations (basis of \(\mathrm{Der}\))
-   `analyze(L)` → structured report of the above
-   Fixture suite: abelian, Heisenberg, \(\mathfrak{sl}_2\), … with expected.toml
-   Tests: unit + identities + examples + basis-change invariance

### Explicitly out of v0.1

-   Parameter families / stratification → **v0.2**
-   CE cohomology / central extensions → **v0.3**
-   Deformations / \(H^2(\mathfrak{g},\mathfrak{g})\) → **v0.4**
-   LieGroups.jl PackageExtension / numerical adapters → **post–v0.1 soft bridge**
-   Group exp/log/hat/vee, SO(3)/SE(3) numerics → **never our job**
-   Root systems, full isomorphism classification, Gröbner engine of our own

### Deliverable

> Exact, tested analysis of concrete finite-dimensional Lie algebras over a
> field (esp. \(\mathbb{Q}\)), with reproducible certificates and a clear API
> we own. Parametric / cohomological differentiation starts *after* this core.

**Single-algebra reference (`src/lie/`):** [docs/lie.md](docs/lie.md)

**Parametric reference (`src/parametric/`):** [docs/parametric.md](docs/parametric.md)

------------------------------------------------------------------------

## v0.2 --- Parameterized Lie Algebras

**Goal (one sentence):**

> Compute classical invariants of a **parameter family** \(L_{\mathbf{a}}\)
> generically, and **certify exceptional loci** where those invariants jump
> (parameter stratification).

The hard core is not “polynomial \(c_{ij}^k\)” alone, but **conditional exact
linear algebra**: Gaussian elimination must branch when a pivot such as
\(a-b\) may vanish.

**Reference:** [docs/parametric.md](docs/parametric.md)

Delivered as a continuous chain:

| Layer | Product |
|-------|---------|
| Conditional LA | pivot certificates under assumption set \(\Sigma\) |
| Analysis tree | shared `CondTree` + incremental `refine` / `analyze_conditional` |
| Stratification | `stratify` → strata + confirmed jump reports |
| Validation | `compare`, specialize ↔ `analyze` fiber checks |

------------------------------------------------------------------------

## v0.3 --- Cohomology

Engine: full exterior algebra, general coefficient modules, induced CE \(d\),
and \(H^k = Z^k/B^k\) with cocycle representatives. Central extensions are a
trivial-module \(H^2\) **test**, not a library module.

In progress:

-   [x] Exterior algebra \(\Lambda^\bullet V\) (`src/cohomology/exterior.jl`)
-   [x] Coefficient modules (trivial + adjoint)
-   [x] Induced \(d\) and lazy cached \(Z/B/H\) (`ce_complex`)
-   [x] Central extensions from trivial-module \(H^2\) (application / test)

Goal: turn cohomology into a practical computational tool. See
[`docs/cohomology.md`](docs/cohomology.md).

------------------------------------------------------------------------

## v0.4 --- Deformation Theory

Implement

-   Infinitesimal deformations
-   NR bracket
-   Obstruction
-   Cohomology over parameter families
-   Isomorphism helpers

------------------------------------------------------------------------

# 5. Core Modules

Current source layout:

-   `src/lie/` — single-algebra structure & invariants (`docs/lie.md`)
-   `src/parametric/` — families, conditional LA, stratification (`docs/parametric.md`)
-   `src/cohomology/` — exterior algebra, modules, CE complex (`docs/cohomology.md`)

------------------------------------------------------------------------

# 6. Distinctive Features

## Exact computation

Everything should be exact whenever possible.

## Parameter awareness

Instead of computing one Lie algebra,

compute **families**.

Example output:

    Generic region:

    center dimension = 0

    a = 0:

    center dimension = 1
    new deformation directions appear

## Explicit certificates

Algorithms should return

-   bases
-   cocycles
-   derivations
-   witnesses

instead of only dimensions.

## Research workflow

The library should help discover mathematics rather than merely execute
formulas.

------------------------------------------------------------------------

# 7. Why Cohomology?

Cohomology unifies several seemingly unrelated computations.

  Object                       Interpretation
  ---------------------------- ----------------
  Center                       H\^0
  Outer derivations            H\^1
  Central extensions           H\^2(K)
  Infinitesimal deformations   H\^2(g)
  Obstructions                 H\^3

Without cohomology:

> "What is this Lie algebra?"

With cohomology:

> "How can this Lie algebra change?"

------------------------------------------------------------------------

# 8. Testing Philosophy

Every mathematical feature must satisfy four requirements.

1.  Hand-computable examples
2.  Mathematical identities
3.  Independent verification
4.  Basis-change invariance

------------------------------------------------------------------------

## Five layers of testing

### Unit tests

Small components.

### Mathematical identities

Jacobi

d²=0

Representation identities

Euler-Poincaré

etc.

### Known examples

Abelian

Heisenberg

sl₂

Bianchi

...

### Differential testing

Compare with

-   SageMath
-   GAP
-   Oscar

### Property-based tests

Random basis changes

Random direct sums

Random extensions

Random parameter specializations

------------------------------------------------------------------------

# 9. Testing Dataset

Classical examples (abelian, Heisenberg, \(\mathfrak{sl}_2\), parametric
\(a\) / \(a-b\) families, …) are exercised directly in
`test/unit/lie/` and `test/unit/parametric/`. Runnable demos live under
`example/`.

------------------------------------------------------------------------

# 10. Future Lean Integration

Long-term idea:

Julia performs computation.

Lean verifies certificates.

Examples:

-   Jacobi
-   Center basis
-   Rank certificates
-   Cohomology representatives
-   Isomorphism matrices

This provides

> fast computation + trusted verification.

------------------------------------------------------------------------

# 11. North-Star Example

``` julia
L = LieAlgebra(...)

check_jacobi(L)

analyze(L)

stratify(L)

H2 = cohomology(L, adjoint_module(L), 2)
```

Possible output:

    Generic region:

    center dimension = 0

    H² = 0

    infinitesimally rigid

    --------------------------------

    a = 0

    center dimension = 1

    H² = 2

    two deformation directions

    higher symmetry detected

------------------------------------------------------------------------

# 12. Guiding Principles

-   Learn by implementing.
-   Every theorem should become an algorithm.
-   Every algorithm should have mathematical tests.
-   Prefer exact algebra over floating point.
-   Reuse Julia's ecosystem instead of reinventing it.
-   Return certificates, not only numbers.
-   Parameterized symbolic computation is the project's defining
    feature.
