# `lie` — Single-algebra core

Reference for the **concrete** Lie-algebra toolkit in `src/lie/`: structure
constants over an AbstractAlgebra/Nemo field, classical invariants, and
[`analyze`](@ref).

Parametric families, conditional branching, and stratification live in
[`parametric`](parametric.md).

---

## Role

Given a finite-dimensional Lie algebra by structure constants over a field
\(F\) (tests use \(F = \mathbb{Q}\)), verify the Lie axioms and compute classical
invariants with **exact certificates** (bases / matrices, not only dimensions).

```text
F + structure constants
        │
        ▼
   LieAlgebra{C}          (types.jl)
        │
        ├─ lie_bracket / ad / check_jacobi
        ├─ change_of_basis
        ├─ LieSubspace / commutator_span
        ├─ center, series, Killing, radical
        ├─ quotient, Levi, ideal decomposition
        ├─ derivations
        ├─ analyze(L) → LieAlgebraReport
        └─ isomorphism / IsoCertificate
```

---

## Module map (`src/lie/`)

| File | Responsibility |
|------|----------------|
| `types.jl` | `LieAlgebra`, `LieAlgebraElem`, `coefficient_ring`, parameters metadata |
| `bracket.jl` | `lie_bracket`, `ad`, `check_jacobi` / `JacobiCertificate` |
| `change_of_basis.jl` | change of basis for algebras and elements |
| `subspace.jl` | `LieSubspace`, span, complement, ideal/subalgebra tests |
| `center.jl` | center as nullspace of the structure-constant equation matrix |
| `series.jl` | derived / lower-central series, solvability, nilpotency |
| `killing.jl` | Killing form, radical, Cartan criteria |
| `radical.jl` | solvable radical (= Cartan orthogonal in char 0) |
| `quotient.jl` | quotient by an ideal |
| `levi.jl` | Levi decomposition (lifting along radical derived series) |
| `ideal_decomp.jl` | adjoint-commutant ideal decomposition / simplicity |
| `derivations.jl` | `Der(L)` via Leibniz linear system |
| `analyze.jl` | `LieAlgebraReport` summary + cached detail accessors |
| `isomorphism.jl` | `isomorphism` / `IsoCertificate` (helpers, not full classification) |

Exact linear algebra uses AbstractAlgebra/Nemo matrices (never `Float64`).

---

## Construction

```julia
using ParametricLieAlgebras
import Nemo

F = Nemo.QQ

# Abelian
L = LieAlgebra(F, 3)

# Degree-of-freedom brackets (recommended)
H = LieAlgebra(F, 3, Dict((1, 2) => [0, 0, 1]))   # [e1,e2]=e3

# Dense n×n×n tensor
# L = LieAlgebra(F, c::Array{C,3})
```

Constructors are distinguished by argument types. Internally both paths store a
dense structure-constant tensor.

---

## Core operations

| API | Meaning |
|-----|---------|
| `lie_bracket(L, x, y)` | structure-constant bracket |
| `ad(L, x)` | adjoint matrix |
| `check_jacobi(L)` | identity certificate; residuals on failure |
| `change_of_basis(L, P)` | isomorphic copy in a new basis |
| `center(L)` | `LieSubspace` certificate |
| `derived_series` / `lower_central_series` | lazy `LieSeries` |
| `derived_algebra(L)` | \([L,L]\) |
| `is_solvable` / `is_nilpotent` | series reaches zero |
| `killing_form` / `killing_rank` / `killing_radical` | Killing data |
| `radical(L)` | solvable radical (char 0) |
| `quotient_algebra(L, I)` | \(L/I\) |
| `levi_decomposition(L)` | radical + Levi factor |
| `ideal_decomposition` / `is_simple` | simple ideals of the Levi quotient |
| `derivations(L)` | basis of derivation matrices |

---

## Algorithms (sketch)

**Center.** Build the \(n^2\times n\) matrix whose nullspace is \(\{z : [x,z]=0\ \forall x\}\);
call field `nullspace`.

**Series.** Iterate commutator spans (`commutator_span`) until dimension drops
to zero or stabilizes; cache terms and layer count.

**Killing.** \(K_{ij}=\mathrm{tr}(\mathrm{ad}_{e_i}\mathrm{ad}_{e_j})\). Radical /
Cartan orthogonal are further nullspaces involving \(K\) and \([L,L]\).

**Derivations.** Encode Leibniz on basis pairs as a homogeneous system in the
\(n^2\) entries of \(D\); nullspace → `Derivations`.

**Levi.** In char 0, \(\mathrm{rad}(L)=[L,L]^\perp\); lift a Levi factor of the
semisimple quotient along the radical derived series by solving linear
correction equations layer by layer.

**Ideal decomposition.** Adjoint commutant via Kronecker nullspace; split by
idempotents from the minimal polynomial; recurse.

---

## `analyze(L)`

Primary UX for a concrete algebra. Summary prints **basis-invariant** facts
(dims, series profiles, flags, Levi kind, Der dim). Bases and matrices stay in
cached certificates:

```julia
r = analyze(L)
center(r)
derived_series(r)
killing_form(r)
levi_decomposition(r)
derivations(r)
```

Example summary shape:

```text
LieAlgebra dim=3 over QQ
Jacobi:           OK
center:           dim 1
derived:          3 → 1 → 0
…
Der:              dim 6
```

---

## Design rules

- Certificates over bare numbers (return subspaces / matrices when meaningful).
- Summary displays stay basis-invariant.
- No LieGroups dependency; no parameter stratification (see [`parametric`](parametric.md)).
- Coefficient domain is an AA/Nemo field (`FieldElem`); parametric rings are
  accepted on `LieAlgebra` for Jacobi identity checks, but field algorithms
  require a field (use `generic_algebra` / `specialize` from the parametric layer).

See also: [parametric.md](parametric.md), [docs/README.md](README.md).
