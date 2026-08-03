@testset "types / coefficient field" begin
    # Primary v0.1 target: exact rationals via Nemo
    L = LieAlgebra(Nemo.QQ, 3)
    @test dim(L) == 3
    @test coefficient_ring(L) === Nemo.QQ
    @test base_ring(L) === Nemo.QQ
    @test L isa LieAlgebra{Nemo.QQFieldElem}

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
