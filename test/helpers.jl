"""
Shared test helpers for ParametricLie.

Use for exact comparisons, random basis changes, and fixture loading.
"""

const FIXTURES_DIR = joinpath(dirname(@__DIR__), "fixtures")

"""
    fixture_path(name)

Absolute path to `fixtures/<name>/`.
"""
fixture_path(name::AbstractString) = joinpath(FIXTURES_DIR, name)

"""
    approx_zero(x; atol)

Placeholder for exact/near-exact zero checks once fields are wired up.
"""
approx_zero(x; atol=0) = iszero(x)
