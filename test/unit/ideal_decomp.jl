@testset "ideal decomposition / commutant" begin
    F = Nemo.QQ

    sl2 = LieAlgebra(F, 3, Dict(
        (1, 2) => [0, 2, 0],
        (1, 3) => [0, 0, -2],
        (2, 3) => [1, 0, 0],
    ))
    H = LieAlgebra(F, 3, Dict((1, 2) => [0, 0, 1]))

    # --- commutant of simple sl₂: scalars only ---
    Csl = adjoint_commutant(sl2)
    @test length(Csl) == 1
    T = Csl[1]
    @test T[1, 1] == T[2, 2] == T[3, 3]
    @test iszero(T[1, 2]) && iszero(T[1, 3]) && iszero(T[2, 3])
    @test !iszero(T[1, 1])
    @test is_simple(sl2)
    @test is_semisimple(sl2)

    ideals_sl = ideal_decomposition(sl2)
    @test length(ideals_sl) == 1
    @test dim(ideals_sl[1]) == 3
    @test is_ideal(sl2, ideals_sl[1])

    # --- non-semisimple rejected ---
    @test_throws ArgumentError ideal_decomposition(H)
    @test !is_simple(H)

    # --- sl₂ ⊕ sl₂ ---
    sl2sl2 = LieAlgebra(F, 6, Dict(
        (1, 2) => [0, 2, 0, 0, 0, 0],
        (1, 3) => [0, 0, -2, 0, 0, 0],
        (2, 3) => [1, 0, 0, 0, 0, 0],
        (4, 5) => [0, 0, 0, 0, 2, 0],
        (4, 6) => [0, 0, 0, 0, 0, -2],
        (5, 6) => [0, 0, 0, 1, 0, 0],
    ))
    @test check_jacobi(sl2sl2).ok
    @test is_semisimple(sl2sl2)
    @test !is_simple(sl2sl2)

    C2 = adjoint_commutant(sl2sl2)
    @test length(C2) == 2

    ideals2 = ideal_decomposition(sl2sl2)
    @test length(ideals2) == 2
    @test sort(dim.(ideals2)) == [3, 3]
    for I in ideals2
        @test is_ideal(sl2sl2, I)
        @test is_subalgebra(sl2sl2, I)
    end
    # Direct sum: basis matrices together invertible
    P = hcat(basis_matrix(ideals2[1]), basis_matrix(ideals2[2]))
    @test AbstractAlgebra.is_unit(AbstractAlgebra.det(P))
    # Cross brackets vanish (direct sum of ideals)
    for a in basis_elems(ideals2[1]), b in basis_elems(ideals2[2])
        @test all(iszero, lie_bracket(sl2sl2, a, b).coords)
    end
    # Each factor is simple
    for I in ideals2
        Lind, _ = ParametricLie._induced_lie_algebra(sl2sl2, I)
        @test is_simple(Lind)
        @test is_semisimple(Lind)
    end

    # --- change-of-basis invariance ---
    Pmix = AbstractAlgebra.matrix(F, 6, 6, [
        1, 0, 0, 1, 0, 0,
        0, 1, 0, 0, 1, 0,
        0, 0, 1, 0, 0, 1,
        0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 1,
    ])
    Lmix = change_of_basis(sl2sl2, Pmix)
    @test is_semisimple(Lmix)
    ideals_m = ideal_decomposition(Lmix)
    @test length(ideals_m) == 2
    @test sort(dim.(ideals_m)) == [3, 3]
    for I in ideals_m
        @test is_ideal(Lmix, I)
    end
    @test AbstractAlgebra.is_unit(AbstractAlgebra.det(hcat(
        basis_matrix(ideals_m[1]), basis_matrix(ideals_m[2])
    )))

    # --- three copies ---
    sl2_3 = LieAlgebra(F, 9, Dict(
        (1, 2) => [0, 2, 0, 0, 0, 0, 0, 0, 0],
        (1, 3) => [0, 0, -2, 0, 0, 0, 0, 0, 0],
        (2, 3) => [1, 0, 0, 0, 0, 0, 0, 0, 0],
        (4, 5) => [0, 0, 0, 0, 2, 0, 0, 0, 0],
        (4, 6) => [0, 0, 0, 0, 0, -2, 0, 0, 0],
        (5, 6) => [0, 0, 0, 1, 0, 0, 0, 0, 0],
        (7, 8) => [0, 0, 0, 0, 0, 0, 0, 2, 0],
        (7, 9) => [0, 0, 0, 0, 0, 0, 0, 0, -2],
        (8, 9) => [0, 0, 0, 0, 0, 0, 1, 0, 0],
    ))
    @test check_jacobi(sl2_3).ok
    @test is_semisimple(sl2_3)
    @test length(adjoint_commutant(sl2_3)) == 3
    ideals3 = ideal_decomposition(sl2_3)
    @test length(ideals3) == 3
    @test sort(dim.(ideals3)) == [3, 3, 3]
    for I in ideals3
        @test is_ideal(sl2_3, I)
    end

    # --- via Levi quotient: sl₂ ⊕ Heisenberg ---
    sl2H = LieAlgebra(F, 6, Dict(
        (1, 2) => [0, 2, 0, 0, 0, 0],
        (1, 3) => [0, 0, -2, 0, 0, 0],
        (2, 3) => [1, 0, 0, 0, 0, 0],
        (4, 5) => [0, 0, 0, 0, 0, 1],
    ))
    lev = levi_decomposition(sl2H)
    @test is_semisimple(lev.quotient)
    ideals_q = ideal_decomposition(lev)
    @test length(ideals_q) == 1
    @test dim(ideals_q[1]) == 3
    @test is_simple(lev.quotient)

    # --- dim 0 ---
    A0 = LieAlgebra(F, 0)
    @test is_semisimple(A0)
    @test isempty(ideal_decomposition(A0))
    @test !is_simple(A0)
end
