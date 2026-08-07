@testset "center" begin
    F = Nemo.QQ

    # --- Abelian / commutative: Z = whole algebra ---
    for n in (0, 1, 2, 5)
        A = LieAlgebra(F, n)
        Z = center(A)
        @test dim(Z) == n
        @test dim(Z) == dim(A)
        @test size(basis_matrix(Z)) == (n, n)
        # brackets with every basis vector vanish
        for z in basis_elems(Z), i in 1:n
            @test iszero(lie_bracket(A, basis_elem(A, i), z).coords)
        end
    end

    # Explicit abelian via dense zero tensor (same corner: E = 0)
    A3 = LieAlgebra(F, fill(F(0), 3, 3, 3))
    @test dim(center(A3)) == 3
    @test check_jacobi(A3).ok

    # --- Heisenberg: Z = span{e3} ---
    H = LieAlgebra(F, 3, Dict((1, 2) => [0, 0, 1]))
    ZH = center(H)
    @test dim(ZH) == 1
    @test check_jacobi(H).ok
    z = only(basis_elems(ZH))
    # must be multiple of e3
    @test iszero(z.coords[1]) && iszero(z.coords[2]) && !iszero(z.coords[3])
    for i in 1:3
        @test iszero(lie_bracket(H, basis_elem(H, i), z).coords)
    end

    # --- sl₂: trivial center ---
    # [H,X]=2X, [H,Y]=-2Y, [X,Y]=H  with e1=H, e2=X, e3=Y
    sl2 = LieAlgebra(F, 3, Dict(
        (1, 2) => [0, 2, 0],
        (1, 3) => [0, 0, -2],
        (2, 3) => [1, 0, 0],
    ))
    @test check_jacobi(sl2).ok
    @test dim(center(sl2)) == 0
    @test size(basis_matrix(center(sl2))) == (3, 0)

    # --- change-of-basis invariance of dim Z ---
    P = AbstractAlgebra.matrix(F, 3, 3, [1, 1, 0, 0, 1, 1, 0, 0, 1])
    Hp = change_of_basis(H, P)
    @test dim(center(Hp)) == dim(center(H)) == 1
end
