# Test fixtures

Reusable Lie algebra examples with expected invariants.

Each subdirectory should eventually contain:

| File | Purpose |
|------|---------|
| `structure.json` or `.jl` | Structure constants / constructor |
| `expected.toml` | Invariants: center dim, Killing rank, series, etc. |

## v0.1 fixtures

- `abelian` — trivial bracket
- `heisenberg` — 3D Heisenberg
- `sl2` — \(\mathfrak{sl}_2\)
- `affine` — affine algebras (e.g. \(\mathfrak{aff}(1)\))
- `bianchi` — Bianchi classification samples
- `nilpotent` — additional nilpotent examples

## Later

- `known_cohomology` — Betti numbers / cocycle witnesses (v0.3+)
