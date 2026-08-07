# Examples

Runnable use cases (not part of the test suite).

```bash
julia --project=. example/heisenberg.jl
julia --project=. example/parametric_family.jl
julia --project=. example/stratify_ab_family.jl
```

| Script | Layer | What it shows |
|--------|-------|----------------|
| [`heisenberg.jl`](heisenberg.jl) | `lie` | 3D Heisenberg over ℚ — `analyze` + detail queries |
| [`parametric_family.jl`](parametric_family.jl) | `parametric` | Jacobi over `QQ[a]`; `analyze_generic` vs `specialize` |
| [`stratify_ab_family.jl`](stratify_ab_family.jl) | `parametric` | Full conditional tree → strata → jump table → fiber validation |
