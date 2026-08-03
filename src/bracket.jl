# Bracket and Jacobi (v0.1)
# Name `lie_bracket` matches LieGroups.jl; we still implement the
# structure-constant case ourselves (see docs/POSITIONING.md D2).

"""
    lie_bracket(L, x, y)

Lie bracket `[x, y]` for algebra `L` (coordinate vectors or algebra elements).

Aligned with LieGroups.jl naming. Structure-constant implementation lives here;
concrete-group brackets stay in LieGroups.jl.
"""
function lie_bracket end

"""
    bracket(L, x, y)

Alias for [`lie_bracket`](@ref). Prefer `lie_bracket` in new code.
"""
bracket(L, x, y) = lie_bracket(L, x, y)

"""
    check_jacobi(L; atol=nothing)

Verify the Jacobi identity. Returns a certificate / report (not only a Bool).
"""
function check_jacobi end
