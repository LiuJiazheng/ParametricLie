@testset "isomorphism helpers" begin
    F = Nemo.QQ
    H = LieAlgebra(F, 3, Dict((1, 2) => [0, 0, 1]))
    sl2 = LieAlgebra(F, 3, Dict(
        (1, 2) => [0, 2, 0],
        (1, 3) => [0, 0, -2],
        (2, 3) => [1, 0, 0],
    ))

    @testset "invariant obstruction" begin
        c = isomorphism(H, sl2)
        @test !c.isomorphic
        @test c.reason === :invariants
        @test !isempty(c.mismatches)
        @test !isomorphic(H, sl2)
    end

    @testset "identity and change_of_basis search" begin
        c0 = isomorphism(H, H)
        @test c0.isomorphic
        @test c0.reason === :identity

        L = LieAlgebra(F, 3, Dict((1, 2) => [0, 1, 0], (1, 3) => [0, 0, 1]))
        P = AbstractAlgebra.matrix(F, 3, 3, [2, 1, 0, 0, 1, 0, 0, 0, 3])
        Lp = change_of_basis(L, P)
        @test is_structure_isomorphism(L, Lp, P)
        c = isomorphism(L, Lp)
        @test c.isomorphic
        @test c.matrix !== nothing
        @test is_structure_isomorphism(L, Lp, c.matrix)
    end
end
