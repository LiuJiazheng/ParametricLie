@testset "jump explanation" begin
    F = Nemo.QQ
    R, (a, b) = AbstractAlgebra.polynomial_ring(F, [:a, :b])
    L = LieAlgebra(R, 3, Dict(
        (1, 2) => [R(0), a, R(0)],
        (1, 3) => [R(0), R(0), a - b],
    ))
    S = stratify(L; invariants = [:center, :derived_dim, :solvability, :nilpotency])
    @test coefficient_ring(S.family) isa AbstractAlgebra.MPolyRing ||
          coefficient_ring(S.family) == R

    Jidx = findfirst(j -> status(j.target.sigma, a - b) === PIVOT_ZERO &&
                          status(j.target.sigma, a) === PIVOT_NONZERO,
                     jump_table(S))
    @test Jidx !== nothing
    J = jump_table(S)[Jidx]
    @test any(e -> e.invariant === :center_dim && e.generic_value == 0 && e.special_value == 1,
              J.changed)

    pg = Dict(a => F(1), b => F(2))
    ps = Dict(a => F(1), b => F(1))
    cause = explain_jump(L, J; point_generic = pg, point_special = ps, order = 2)
    @test cause.wall_class === :nontrivial
    @test cause.verdict === :integrable_deformation
    @test cause.iso !== nothing
    @test !cause.iso.isomorphic
    @test cause.iso.reason === :invariants
    @test !all(iszero, cause.wall_cocycle)

    explain_jumps!(S; order = 2, points = Dict(
        findfirst(s -> s.sigma == S.generic.sigma, S.strata) => pg,
        findfirst(s -> s.sigma == J.target.sigma, S.strata) => ps,
    ))
    @test J.cause !== nothing
    @test J.cause.verdict === :integrable_deformation
end
