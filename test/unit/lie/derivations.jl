@testset "derivations" begin
    F = Nemo.QQ

    # --- abelian: Der ≅ gl(n), dim = n² ---
    for n in (0, 1, 2, 3)
        A = LieAlgebra(F, n)
        DA = derivations(A)
        @test DA isa Derivations
        @test dim(DA) == n * n
        @test length(basis_matrices(DA)) == n * n
        for D in basis_matrices(DA)
            @test is_derivation(A, D)
            @test size(D) == (n, n)
        end
    end

    # Standard matrix units span Der for abelian n=2
    A2 = LieAlgebra(F, 2)
    mats2 = basis_matrices(derivations(A2))
    E11 = AbstractAlgebra.matrix(F, 2, 2, [1, 0, 0, 0])
    E12 = AbstractAlgebra.matrix(F, 2, 2, [0, 1, 0, 0])
    E21 = AbstractAlgebra.matrix(F, 2, 2, [0, 0, 1, 0])
    E22 = AbstractAlgebra.matrix(F, 2, 2, [0, 0, 0, 1])
    for E in (E11, E12, E21, E22)
        @test ParametricLieAlgebras._matrix_in_span(E, mats2)
        @test is_derivation(A2, E)
    end

    # --- Heisenberg: dim Der = 6 ---
    H = LieAlgebra(F, 3, Dict((1, 2) => [0, 0, 1]))
    @test check_jacobi(H).ok
    DH = derivations(H)
    @test dim(DH) == 6
    for D in basis_matrices(DH)
        @test is_derivation(H, D)
    end
    # Inner derivations ad(e_i) ⊆ Der
    for i in 1:3
        Ai = ad(H, basis_elem(H, i))
        @test is_derivation(H, Ai)
        @test ParametricLieAlgebras._matrix_in_span(Ai, basis_matrices(DH))
    end
    # apply_derivation sanity: D(e1) coords = column 1 of D
    D0 = first(basis_matrices(DH))
    y = apply_derivation(H, D0, basis_elem(H, 1))
    @test y.coords == [D0[1, 1], D0[2, 1], D0[3, 1]]

    # --- sl₂: Der = Inn ≅ sl₂, dim 3 ---
    sl2 = LieAlgebra(F, 3, Dict(
        (1, 2) => [0, 2, 0],
        (1, 3) => [0, 0, -2],
        (2, 3) => [1, 0, 0],
    ))
    @test check_jacobi(sl2).ok
    Dsl = derivations(sl2)
    @test dim(Dsl) == 3
    for D in basis_matrices(Dsl)
        @test is_derivation(sl2, D)
    end
    # Every ad(e_i) lies in Der, and Inn spans Der
    ads = [ad(sl2, basis_elem(sl2, i)) for i in 1:3]
    for Ai in ads
        @test ParametricLieAlgebras._matrix_in_span(Ai, basis_matrices(Dsl))
    end
    for D in basis_matrices(Dsl)
        @test ParametricLieAlgebras._matrix_in_span(D, ads)
    end

    # Non-derivation witness: identity is not a derivation of sl₂
    I3 = AbstractAlgebra.identity_matrix(F, 3)
    @test !is_derivation(sl2, I3)

    # --- change-of-basis: dimension invariant ---
    P = AbstractAlgebra.matrix(F, 3, 3, [1, 1, 0, 0, 1, 1, 0, 0, 1])
    for L in (H, sl2)
        Lp = change_of_basis(L, P)
        @test dim(derivations(Lp)) == dim(derivations(L))
    end

    # --- sl₂ ⊕ ℚ: Der larger than Inn ---
    # e1..e3 = sl₂, e4 central abelian factor
    sl2a = LieAlgebra(F, 4, Dict(
        (1, 2) => [0, 2, 0, 0],
        (1, 3) => [0, 0, -2, 0],
        (2, 3) => [1, 0, 0, 0],
    ))
    @test check_jacobi(sl2a).ok
    Da = derivations(sl2a)
    # Inn has dim 3 (center trivial on sl2 factor; ad(e4)=0); Der includes
    # scaling of the abelian summand and maps into/out of it → dim ≥ 4
    @test dim(Da) >= 4
    for D in basis_matrices(Da)
        @test is_derivation(sl2a, D)
    end
end
