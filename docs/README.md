# Documentation

Module references for ParametricLieAlgebras. Start from the root [README](../README.md)
for motivation, install, and quick start.

| Guide | Source tree | Topics |
|-------|-------------|--------|
| [lie.md](lie.md) | `src/lie/` | Structure constants, invariants, `analyze`, isomorphism |
| [parametric.md](parametric.md) | `src/parametric/` | Families, conditional LA, `stratify`, jump causes |
| [cohomology.md](cohomology.md) | `src/cohomology/` | Exterior algebra, CE complex, central extensions |
| [deformations.md](deformations.md) | `src/deformations/` | NR bracket, Maurer–Cartan, gauge equivalence |

## Suggested reading order

1. **Concrete algebra** — [lie.md](lie.md) (`analyze`, certificates)
2. **Families & strata** — [parametric.md](parametric.md) (`stratify`, `validate_stratum`)
3. **Cohomology** — [cohomology.md](cohomology.md) (Betti numbers, central extensions)
4. **Deformations** — [deformations.md](deformations.md) (MC, gauge, rigidity)
5. **Glue** — jump causes in [parametric.md](parametric.md#jump-cause--isomorphism) and
   `example/jump_cause_iso.jl`

## Examples

Runnable scripts: [example/README.md](../example/README.md).
