@testset "killing / Cartan" begin
    F = Nemo.QQ

    # --- abelian: K = 0, τ = L, solvable, not semisimple (n>0) ---
    A = LieAlgebra(F, 3)
    KA = killing_form(A)
    @test size(KA) == (3, 3)
    @test all(iszero, KA)
    @test killing_rank(A) == 0
    @test dim(killing_radical(A)) == 3
    @test dim(cartan_orthogonal(A)) == 3
    @test is_cartan_solvable(A)
    @test is_solvable(A)
    @test !is_semisimple(A)

    A0 = LieAlgebra(F, 0)
    @test size(killing_form(A0)) == (0, 0)
    @test killing_rank(A0) == 0
    @test is_cartan_solvable(A0)
    @test is_semisimple(A0)  # vacuous

    # --- Heisenberg: K = 0, [L,L] = Z = span{e3}, τ = L ---
    H = LieAlgebra(F, 3, Dict((1, 2) => [0, 0, 1]))
    @test check_jacobi(H).ok
    @test killing_rank(H) == 0
    @test all(iszero, killing_form(H))
    DH = derived_algebra(H)
    @test dim(DH) == 1
    τH = cartan_orthogonal(H)
    @test dim(τH) == 3
    @test is_cartan_solvable(H)
    @test is_solvable(H) == is_cartan_solvable(H)
    @test !is_semisimple(H)

    # Explicit algorithm: nullspace(Dᵀ K) with K = 0 → full space
    KH = killing_form(H)
    DH_mat = basis_matrix(DH)
    EH = transpose(DH_mat) * KH
    _nH, NH = AbstractAlgebra.nullspace(EH)
    @test size(NH) == (3, 3)

    # --- sl₂: nondegenerate Killing, not solvable, semisimple ---
    # [H,X]=2X, [H,Y]=-2Y, [X,Y]=H  with e1=H, e2=X, e3=Y
    sl2 = LieAlgebra(F, 3, Dict(
        (1, 2) => [0, 2, 0],
        (1, 3) => [0, 0, -2],
        (2, 3) => [1, 0, 0],
    ))
    @test check_jacobi(sl2).ok
    Ksl = killing_form(sl2)
    # Classical matrix: diag(8,0,0) with off-diagonal K[2,3]=K[3,2]=4
    @test Ksl[1, 1] == F(8)
    @test Ksl[2, 3] == F(4) && Ksl[3, 2] == F(4)
    @test iszero(Ksl[1, 2]) && iszero(Ksl[1, 3]) && iszero(Ksl[2, 2]) && iszero(Ksl[3, 3])
    @test killing_rank(sl2) == 3
    @test dim(killing_radical(sl2)) == 0
    @test is_semisimple(sl2)
    @test !is_cartan_solvable(sl2)
    @test is_solvable(sl2) == is_cartan_solvable(sl2)
    # [sl₂,sl₂] = sl₂ ⇒ τ = rad(B) = 0
    @test dim(derived_algebra(sl2)) == 3
    @test dim(cartan_orthogonal(sl2)) == 0

    # Symmetry of K
    for L in (A, H, sl2)
        K = killing_form(L)
        n = dim(L)
        for i in 1:n, j in 1:n
            @test K[i, j] == K[j, i]
        end
    end

    # Center ⊆ Killing radical
    for L in (A, H, sl2)
        Z = center(L)
        R = killing_radical(L)
        # every central basis vector is in ker(K)
        K = killing_form(L)
        for z in basis_elems(Z)
            v = AbstractAlgebra.matrix(F, dim(L), 1, z.coords)
            @test all(iszero, K * v)
        end
        @test dim(Z) <= dim(R)
    end

    # Change-of-basis: rank / Cartan / semisimplicity invariant
    P = AbstractAlgebra.matrix(F, 3, 3, [1, 1, 0, 0, 1, 1, 0, 0, 1])
    for L in (H, sl2)
        Lp = change_of_basis(L, P)
        @test killing_rank(Lp) == killing_rank(L)
        @test is_cartan_solvable(Lp) == is_cartan_solvable(L)
        @test is_semisimple(Lp) == is_semisimple(L)
        @test dim(cartan_orthogonal(Lp)) == dim(cartan_orthogonal(L))
    end
end
