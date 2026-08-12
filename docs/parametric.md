# `parametric` — Families, conditional analysis, stratification

Reference for the **parametric layer** in `src/parametric/`: coefficient rings
with free parameters, assumption-aware linear algebra, shared conditional
analysis trees, stratification, and fiber validation against the
[`lie`](lie.md) core.

---

## Role in the stack

```text
structure constants in k[a] or Frac(k[a])
        │
        ▼
  Jacobi identity (exact)
        │
        ▼
  stratify → regions + jumps vs generic     ← user view
        │         ▲
        │         │  CondTree / refine      ← internal IR
        │         │  conditional LA under Σ
        ▼
  specialize a chosen stratum → fiber over k
        │
        ▼
  analyze / Levi / Der / …                  ← src/lie on the fiber
        │
        ▼
  validate_stratum (optional sanity check)
```

**Core rule:** never divide by a parameter-dependent expression silently.
Every such division is either proven nonzero on the current branch or produces
an explicit \(p=0\) / \(p\neq 0\) split with a pivot certificate.

---

## Module map (`src/parametric/`)

| File | Responsibility |
|------|----------------|
| `specialize.jl` | `generic_algebra`, `specialize`, `analyze_generic` |
| `assumptions.jl` | \(\Sigma\), `status`, normalize, `assume_zero` / `assume_nonzero` |
| `conditional_linalg.jl` | conditional `rank` / `nullspace` / `rref` / `solve` |
| `cond_tree.jl` | `CondTree`, `refine`, `analyze_conditional`, invariant clients |
| `stratify.jl` | `stratify`, `Stratum`, `JumpReport` |
| `strata_compare.jl` | `compare`, `validate_stratum` |
| `jump_explain.jl` | `explain_jumps!`, `JumpCause`, `wall_cocycle` |

---

## Coefficient domains and Jacobi

Parameters are algebraically independent indeterminates. Over \(R=k[a_i]\) or
\(K=\mathrm{Frac}(R)\), `check_jacobi` tests the **polynomial / rational
identity** — it does not solve zero loci. Non-identity residuals are returned
explicitly (`JacobiResidual`).

`LieAlgebra` may carry `parameters` and `domain_denominators` (nonvanishing
conditions for rational structure constants).

---

## Specialization and generic analysis

| API | Meaning |
|-----|---------|
| `generic_algebra(L)` | view over a field (`Frac(R)` when needed) |
| `analyze_generic(L)` | `analyze` on that view (Zariski-open validity) |
| `specialize(L, assignments)` | exact fiber over the base field \(k\) |

Generic analysis **does not** certify the exceptional locus; that requires the
conditional tree.

---

## Assumption set \(\Sigma\)

\(\Sigma\) is a conjunction of equalities \(f_i=0\) and nonvanishings \(g_j\neq 0\).

```text
status(Σ, p) ∈ { ZERO, NONZERO, UNKNOWN }
```

Conservative three-valued logic: return ZERO/NONZERO only when **proven**.
Normalization (content, monic, square-free when cheap) and factor-aware
reasoning reduce duplicate or vacuous branches. Inconsistent branches are
pruned when emptiness is proven.

---

## Conditional linear algebra

Primitives under \(\Sigma\):

- `conditional_rank`
- `conditional_nullspace`
- `conditional_rref`
- `conditional_solve`

Behaviour:

1. Prove before split (`status`).
2. On `UNKNOWN`, split **lazily** (generic-first: \(p\neq 0\) then \(p=0\)).
3. Prefer symbolically simple pivots.
4. Respect `BranchBudget` (`max_branches`, `max_depth`, `max_expression_size`);
   on exhaustion return **incomplete** leaves — never silent truncation.

Each split records a `PivotCertificate` (polynomial, row/col, branch kind).

Vanishing of a rational pivot depends on its **numerator**; denominators must
already be nonzero on the branch.

---

## User path vs internal IR

**User path** (what `show` emphasizes):

```julia
S = stratify(L; invariants = [:center, :derived_dim, :solvability, :nilpotency])
# pick a Stratum → specialize → full lie toolkit on the fiber
Lf = specialize(L, Dict(a => QQ(1), b => QQ(1)))
analyze(Lf)                 # center, derived, solvability, …
levi_decomposition(Lf)      # and any other src/lie API
validate_stratum(L, st, point)
```

`show(S)` prints the certified parameter regions: one **GENERIC** signature,
then each **JUMP** as a diff against generic (the two ends of the change).

**Internal IR:** `CondTree` / `analyze_conditional` / `refine` is the shared
assumption tree behind `stratify`. Prefer not to present it as a second
classification table; inspect `S.tree` when debugging incomplete branches or
certificates.

```julia
T = cond_tree(L)
T = refine(T, :center)
T = refine(T, :derived_dim)
# or
T = analyze_conditional(L; invariants = …)
S = stratify(T)
```

Each leaf stores \(\Sigma\), pivot trail, and an invariant bag. Refinement:

1. Read current \(\Sigma\).
2. If the invariant is already certified → **inherit**, do not recompute.
3. Else run the matrix client; on UNKNOWN pivots, subdivide **only that leaf**.
4. Attach certificates to children.

A leaf is **complete for a suite** when every requested invariant is certified
and the leaf is not budget-incomplete (`is_complete`).

### Supported refine symbols

| Symbol | Attaches (among others) |
|--------|-------------------------|
| `:center` | `center_dim`, `center_basis` |
| `:killing_rank` / `:killing_radical` | Killing rank / radical basis |
| `:derived_dim` / `:derived_profile` | derived dim / series profile |
| `:solvability` / `:nilpotency` | boolean flags (via series profiles) |
| `:radical` | solvable radical under \(\Sigma\) |
| `:derivations` | `Der` dimension / matrices |

Levi / ideal decomposition are **fiber** APIs (`src/lie`): run them after
`specialize`, not as parametric `refine` clients (yet).

---

## Stratification and jumps

```julia
S = stratify(L; invariants = …)   # or stratify(T; …)
```

Converts complete leaves into `Stratum` objects (assumptions, signature,
witnesses, trail). Identifies a **generic** stratum (prefer many nonvanishings,
few equalities).

Keep three concepts distinct:

| Concept | Meaning |
|---------|---------|
| Generic validity condition | \(p\neq 0\) obligations used on the generic path |
| Algorithmic exceptional condition | locus where a generic certificate *may* fail |
| Confirmed invariant jump | after analyzing \(p=0\), some invariant actually differs |

`jump_table(S)` lists only confirmed jumps (nonempty signature changes); the
default `show(S)` already inlines those diffs.

The stratification is **certified but not necessarily geometrically minimal**
(no primary decomposition / canonical constructible simplification required).

---

## Comparison and validation

```julia
compare(S1, S2)           # same_signature ≠ isomorphic
validate_stratum(L, stratum, point)
```

Fiber check: if `point` satisfies \(\Sigma\), then
`analyze(specialize(L, point))` must agree with the stratum signature
(center dim, derived dim, Killing rank, solvability, nilpotency, …).

This closes the loop from parametric certificates back to the trusted
[`lie`](lie.md) core.

---

## End-to-end example

Family over \(\mathbb{Q}[a,b]\):

\[
[e_1,e_2]=a e_2,\qquad [e_1,e_3]=(a-b)e_3,\qquad [e_2,e_3]=0.
\]

`stratify` yields four certified regions; jumps vs generic:

| \(\Sigma\) | vs generic |
|------------|------------|
| \(a\neq 0,\ a-b\neq 0\) | GENERIC (`center=0`, `derived=2`) |
| \(a\neq 0,\ a-b=0\) | `center 0→1`, `derived 2→1` |
| \(a=0,\ a-b\neq 0\) | same dim jump |
| \(a=a-b=0\) | dims + becomes nilpotent |

Then pick e.g. \((a,b)=(1,1)\), `specialize`, and run `analyze` /
`levi_decomposition` on the fiber (this family is solvable everywhere → Levi
trivial). Fiber signatures match the strata under `validate_stratum`.

### Jump cause + isomorphism

After `stratify`, annotate confirmed jumps with deformation / iso certificates:

```julia
S = stratify(L; invariants = [:center, :derived_dim])
explain_jumps!(S; order = 2, points = Dict(1 => pg, 2 => ps))
J = jump_table(S)[1]
J.cause.verdict          # e.g. :integrable_deformation
J.cause.wall_class       # :nontrivial / :coboundary / …
J.cause.iso              # IsoCertificate between fiber representatives
```

| API | Role |
|-----|------|
| `explain_jump` / `explain_jumps!` | wall cocycle + `H²` + MC + fiber `isomorphism` |
| `wall_cocycle(family, point, dir)` | `∂_t μ` as adjoint `C²` cochain |
| `isomorphism(L, L′)` | find `P` or invariant obstruction (`IsoCertificate`) |

Runnable narrative: `example/jump_cause_iso.jl`
(Jump `center 0→1` → nontrivial integrable wall class → no fiber isomorphism).

Runnable script: `example/stratify_ab_family.jl`.

---

## Design rules

- Branch explicitly; certificates over samples.
- One shared tree; inherit parent results when still valid.
- Specialize to `lie`/`analyze` for validation.
- Budgets fail loudly (`incomplete`), never invent answers.
- Prefer one-parameter / hypersurface demos; multivariate geometry is opt-in
  complexity.

See also: [lie.md](lie.md), [deformations.md](deformations.md), [docs/README.md](README.md).
