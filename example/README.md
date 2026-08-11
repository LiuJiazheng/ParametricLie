# Examples

Runnable use cases (not part of the test suite).

```bash
julia --project=. example/heisenberg.jl
julia --project=. example/parametric_family.jl
julia --project=. example/stratify_ab_family.jl
julia --project=. example/central_extension.jl
julia --project=. example/heisenberg_adjoint_cohomology.jl
```

| Script | Layer | What it shows |
|--------|-------|----------------|
| [`heisenberg.jl`](heisenberg.jl) | `lie` | 3D Heisenberg over ℚ — `analyze` + detail queries |
| [`parametric_family.jl`](parametric_family.jl) | `parametric` | Jacobi over `QQ[a]`; `analyze_generic` vs `specialize` |
| [`stratify_ab_family.jl`](stratify_ab_family.jl) | `parametric` | `stratify` → specialize a stratum → `analyze` / Levi; CondTree as IR |
| [`central_extension.jl`](central_extension.jl) | `cohomology` | Heisenberg as nonsplit central extension of abelian \(\mathbb{Q}^2\) via trivial \(H^2\) |
| [`heisenberg_adjoint_cohomology.jl`](heisenberg_adjoint_cohomology.jl) | `cohomology` | Heisenberg \(H^\bullet(H,H)\): Betti \(1,4,5,2\), Out / deformations |
