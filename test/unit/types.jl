@testset "types / coefficient field" begin
    # Primary v0.1 target: exact rationals via Nemo (abelian = zero tensor)
    L = LieAlgebra(Nemo.QQ, 3)
    @test dim(L) == 3
    @test coefficient_ring(L) === Nemo.QQ
    @test base_ring(L) === Nemo.QQ
    @test L isa LieAlgebra{Nemo.QQFieldElem}
    @test size(structure_constants(L)) == (3, 3, 3)
    @test all(iszero, structure_constants(L))

    # Explicit structure constants (Heisenberg-style [e1,e2]=e3) — dense
    C = Nemo.QQFieldElem
    c = fill(Nemo.QQ(0), 3, 3, 3)
    c[1, 2, 3] = Nemo.QQ(1)
    c[2, 1, 3] = Nemo.QQ(-1)
    H = LieAlgebra(Nemo.QQ, c)
    @test dim(H) == 3
    @test structure_constants(H)[1, 2, 3] == Nemo.QQ(1)

    # Same Heisenberg via DOF / sparse-pair constructor
    H2 = LieAlgebra(Nemo.QQ, 3, Dict((1, 2) => [0, 0, 1]))
    @test structure_constants(H2) == structure_constants(H)

    # i > j is accepted and antisymmetrized
    H3 = LieAlgebra(Nemo.QQ, 3, [(2, 1) => [0, 0, -1]])
    @test structure_constants(H3) == structure_constants(H)

    # Finite field — same API, different parent
    F5 = Nemo.GF(5)
    L5 = LieAlgebra(F5, 2)
    @test dim(L5) == 2
    @test coefficient_ring(L5) === F5
    @test L5 isa LieAlgebra{Nemo.elem_type(F5)}

    # Parameter-ready: frac field of a polynomial ring (v0.2 substrate)
    P, _ = AbstractAlgebra.polynomial_ring(Nemo.QQ, :t)
    Kt = AbstractAlgebra.fraction_field(P)
    Lt = LieAlgebra(Kt, 1)
    @test coefficient_ring(Lt) === Kt
end
