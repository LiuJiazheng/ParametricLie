# Positioning & Differentiation

> **One line:** LieGroups.jl studies geometry and numerics on *known* Lie groups;
> we study the *symbolic structure* of finite-dimensional Lie algebras given by
> structure constants — especially parameterized families, cohomology, and deformations.

Brand keywords: **Exact · Parametric · Cohomological · Structure-aware**

---

## Elevator pitch

**EN:** ParametricLie is an exact symbolic toolkit for finite-dimensional Lie
algebras defined by structure constants, with emphasis on parameterized
families, structural stratification, cohomology, extensions, and deformations.
It complements LieGroups.jl (smooth Lie-group geometry and numerical group ops),
rather than replacing it.

**中文:** 由结构常数驱动的有限维李代数精确符号计算库，重点支持参数族、结构分层、
上同调、扩张与形变。它是 LieGroups.jl 的代数与符号补充，不是替代品。

---

## Not competitors

### vs [LieGroups.jl](https://github.com/JuliaManifolds/LieGroups.jl) (JuliaManifolds)

| | LieGroups.jl | This project |
|--|--------------|--------------|
| Starting point | Concrete group \(G\) | Structure constants \(c_{ij}^k\) over a field/ring |
| Lie algebra means | \(T_e G\) (tangent space) | Algebra \((V,[\cdot,\cdot])\) |
| Core ops | compose, inv, exp, log, hat/vee, group actions | Jacobi, center, series, Killing, Der, \(H^\bullet\), stratify |
| Scalars | Float / manifold coords | Exact: \(\mathbb{Q}\), \(\mathbb{F}_p\), \(\mathbb{Q}(a,b)\), … |
| Typical question | How do I move on \(G\)? | Where does structure jump in the family \(g_a\)? |

**Overlap is thin:** bracket, bases, coordinates, a few classical examples.
We do **not** reimplement `exp`/`log`/`hat`/`vee`/`compose` for SO(3)/SE(3).

### vs [LieTypes.jl](https://github.com/JoanAguilar/LieTypes.jl)

Lightweight robot/pose types (`SO2`/`SE2`/`SO3`/`SE3`, dual quaternions) with
multiply / inv / exp / log. Almost no overlap with parametric abstract Lie
algebras or cohomology.

### Real overlap risk: CAS Lie algebras

Closer peers: **Oscar.jl**, GAP Lie packages, SageMath Lie algebras.
Our durable differentiator must stay:

> **parameterized structure analysis + cohomology witnesses + deformations +
> automatic stratum detection**

—not merely “another Julia Lie algebra library.”

---

## Meaning of “Lie algebra”

| Ecosystem | Meaning |
|-----------|---------|
| LieGroups.jl | Object attached to a group; matrices / tangent reps / `ArrayPartition` |
| This project | Base + coefficient field + structure constants; e.g. `(3//2)*e₁ + (a²-1)*e₃` |

Same name, different center of gravity.

---

## Engineering value (honest)

Not “make your robot 30% faster.” Upstream **math infrastructure**:

```
discover / verify structure  →  generate executable model  →  numerical stack
```

Concrete uses:

1. **Model verifier** — Jacobi, center, series, hidden symmetry for custom/extended state algebras  
2. **Symbolic front-end** — validated ad / bracket / Jacobians / BCH stubs → LieGroups / SciML / codegen  
3. **Parameter degeneration** — exact exceptional strata (`a=0`, `a=b`, …) before numerics blow up  
4. **Isomorphism / convention bridge** — invariants + change-of-basis between papers/codebases  
5. **Central extensions** — \(H^2\) candidates for bias/drift/gauge states (algebraically complete)  
6. **Deformation filter** — coboundary vs nontrivial \(H^2(\mathfrak{g},\mathfrak{g})\)  
7. **Research infrastructure** — certificates, fixtures, future Lean; not one-off paper scripts  

**False claims to avoid:** replace LieGroups.jl; general ML speedups; AD as the product; cohomology solves robotics by itself.

---

## Locked design decisions

### D1 — Complement, do not compete on numerics

Out of scope (use LieGroups.jl / Manifolds.jl):

- group multiply / inverse  
- exp / log on concrete groups  
- manifold metrics, group actions  
- SO(n) / SE(n) / SU(n) numerical kernels  

In scope: structure-constant algebras over AA/Nemo fields; invariants; parameters; cohomology; deformations.

### D2 — API names align with LieGroups where meanings match

Prefer LieGroups naming for shared concepts:

| Prefer | Alias (optional) |
|--------|------------------|
| `lie_bracket(L, x, y)` | `bracket` |
| `lie_bracket!(…)` when mutating | — |

We still **implement** structure-constant `lie_bracket` ourselves (LieGroups
only brackets `LieAlgebra(G)` for concrete \(G\)). Alignment is **interface**,
not “delete our bracket.”

### D3 — Own the symbolic core; bridge later

**v0.1 owns** structure-constant `lie_bracket`, Jacobi, invariants, certificates.
Naming may follow LieGroups; implementation does not.

**After v0.1**, optional soft bridge:

```
LieGroups.jl  ←→  (adapter / PackageExtension)  ←→  ParametricLie
                         ↓
              AbstractAlgebra / Nemo / (Oscar)
```

- `structure_constant_view(LieAlgebra(G); basis=…)` → our `LieAlgebra`  
- specialize exact tensors → Float64 / arb for SciML / LieGroups  
- best-effort `realize_as_lie_group` only when recognition is justified  

No hard dependency on LieGroups.jl in v0.1.

### D4 — Coefficient fields via AbstractAlgebra (not a custom Field trait)

See README / `src/types.jl`. Primary tests on `Nemo.QQ`; parametric families via
fraction fields of polynomial rings.

### D5 — Package name (deferred)

Working name: `ParametricLie`. Stronger candidate later:
`ParametricLieAlgebras.jl` (makes “algebras, not groups” explicit).
Rename after v0.1 prototype stabilizes.

### D6 — Certificates over bare numbers

Algorithms return bases, cocycles, matrices, strata — not only dimensions —
so results can feed adapters, codegen, and (eventually) Lean.

---

## North-star workflow

```julia
L = LieAlgebra(F, struct_consts)   # F = QQ or QQ(a); n×n×n tensor

check_jacobi(L)
analyze(L)
stratify(L)                        # v0.2+

H2 = cohomology(L, adjoint_module(L), 2)   # v0.3+
```

Optional later:

```julia
using LieGroups
g = LieAlgebra(SpecialOrthogonalGroup(3))
L = structure_constant_view(g)     # adapter
analyze(L)

# or: specialize symbolic tensors → numerical kernels for ODE / optimization
```

---

## Decision checklist (when adding features)

1. Is this SO(3)/SE(3)/exp/log/hat/vee numerics? → **LieGroups**, not us.  
2. Does a shared concept already have a LieGroups name? → **reuse that name**.  
3. Can Oscar/GAP/Sage already do it for a *single* algebra over ℚ? → OK to overlap
   only if we add **parameters / strata / witnesses / deformations**.  
4. Does the feature need Manifolds at compile time? → **PackageExtension / adapter**,
   not a hard dep in v0.1.  
5. Exact first; floats only at the adapter boundary.
