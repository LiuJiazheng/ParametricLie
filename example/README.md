# Examples

Runnable use cases for ParametricLie.jl (not part of the test suite).

```bash
julia --project=. example/heisenberg.jl
julia --project=. example/parametric_family.jl
```

| Script | What it shows |
|--------|----------------|
| [`heisenberg.jl`](heisenberg.jl) | v0.1: 3D Heisenberg over `ℚ` — `analyze` + series / Levi / Der details |
| [`parametric_family.jl`](parametric_family.jl) | v0.2.1 + v0.2.2: Jacobi identity over `QQ[a]`, then `analyze_generic` vs `specialize` on `[e1,e2]=a e2` |
