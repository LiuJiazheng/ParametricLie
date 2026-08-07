@testset "v0.2.2 specialize + analyze_generic" begin
    R, (a,) = AbstractAlgebra.polynomial_ring(Nemo.QQ, [:a])
    F = Nemo.QQ

    # Core family: [e1, e2] = a e2
    L_R = LieAlgebra(R, 2, Dict((1, 2) => [R(0), a]))
    @test check_jacobi(L_R).ok

    # --- generic over QQ(a) ---
    L_K = generic_algebra(L_R)
    @test coefficient_ring(L_K) isa AbstractAlgebra.FracField
    @test check_jacobi(L_K).ok

    r_gen = analyze_generic(L_R)
    @test r_gen isa LieAlgebraReport
    @test dim(center(r_gen)) == 0
    @test dim(derived_algebra(generic_algebra(L_R))) == 1
    @test is_solvable(r_gen)
    @test !is_nilpotent(r_gen)
    # Same as analyzing the Frac algebra directly
    r_K = analyze(L_K)
    @test dim(center(r_K)) == dim(center(r_gen))
    @test is_solvable(r_K) == is_solvable(r_gen)
    @test is_nilpotent(r_K) == is_nilpotent(r_gen)

    # --- fiber a = 0: abelian ---
    L0 = specialize(L_R, Dict(a => F(0)))
    @test coefficient_ring(L0) === F
    @test check_jacobi(L0).ok
    @test all(iszero, structure_constants(L0))

    r0 = analyze(L0)
    @test dim(center(r0)) == 2
    @test dim(derived_algebra(L0)) == 0
    @test is_solvable(r0) && is_nilpotent(r0)

    # The point of v0.2.2: generic ≠ special fiber (user chose a=0)
    @test dim(center(r_gen)) != dim(center(r0))
    @test is_nilpotent(r_gen) != is_nilpotent(r0)

    # Fiber a = 1: non-abelian 2D (same qualitative invariants as generic)
    L1 = specialize(L_R, Dict(:a => 1))
    @test check_jacobi(L1).ok
    r1 = analyze(L1)
    @test dim(center(r1)) == 0
    @test dim(derived_algebra(L1)) == 1
    @test is_solvable(r1) && !is_nilpotent(r1)

    # Vector assignment in parameter order
    L0v = specialize(L_R, [F(0)])
    @test structure_constants(L0v) == structure_constants(L0)

    # Specialize from Frac family
    L0f = specialize(L_K, Dict(a => 0))
    @test dim(center(analyze(L0f))) == 2

    # Pole: denominator vanishes
    K = AbstractAlgebra.fraction_field(R)
    Lpole = LieAlgebra(K, 2, Dict((1, 2) => [K(0), a // (a - 1)]))
    @test_throws ArgumentError specialize(Lpole, Dict(a => 1))
    # Defined fiber still works
    L2 = specialize(Lpole, Dict(a => 2))
    @test structure_constants(L2)[1, 2, 2] == F(2) * inv(F(1))  # 2/(2-1) = 2
end
