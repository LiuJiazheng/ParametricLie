# Examples

Runnable demos (not part of the test suite). From the package root:

```bash
julia --project=. example/<script>.jl
```

| Script | Layer | What it shows |
|--------|-------|----------------|
| [`heisenberg.jl`](heisenberg.jl) | `lie` | Heisenberg over ℚ — `analyze` |
| [`parametric_family.jl`](parametric_family.jl) | `parametric` | Jacobi over `QQ[a]`; generic vs specialize |
| [`stratify_ab_family.jl`](stratify_ab_family.jl) | `parametric` | `stratify` → specialize → `analyze` / Levi |
| [`central_extension.jl`](central_extension.jl) | `cohomology` | Heisenberg as central extension via trivial \(H^2\) |
| [`heisenberg_adjoint_cohomology.jl`](heisenberg_adjoint_cohomology.jl) | `cohomology` | Adjoint Betti \(1,4,5,2\); Out / deformations |
| [`jump_cause_iso.jl`](jump_cause_iso.jl) | `parametric` + deformations | Jump → `explain_jumps!` → fiber `isomorphism` |

Module docs: [docs/README.md](../docs/README.md).
