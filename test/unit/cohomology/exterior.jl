@testset "exterior algebra" begin
    QQ = Nemo.QQ

    @testset "dims / graded-lex indexing" begin
        Λ = exterior_algebra(QQ, 3)
        @test Λ isa ExteriorAlgebra{Nemo.QQFieldElem}
        @test ambient_dim(Λ) == 3
        @test dim(Λ) == 8
        @test dim(Λ, 0) == 1
        @test dim(Λ, 1) == 3
        @test dim(Λ, 2) == 3
        @test dim(Λ, 3) == 1
        @test dim(Λ, 4) == 0
        @test sum(dim(Λ, k) for k in 0:3) == dim(Λ)

        @test multi_indices(3, 0) == [Int[]]
        @test multi_indices(3, 2) == [[1, 2], [1, 3], [2, 3]]
        @test coord_index(3, Int[]) == 1
        @test coord_index(3, Int[1]) == 2
        @test coord_index(3, Int[2]) == 3
        @test coord_index(3, Int[3]) == 4
        @test coord_index(3, Int[1, 2]) == 5
        @test coord_index(3, Int[1, 3]) == 6
        @test coord_index(3, Int[2, 3]) == 7
        @test coord_index(3, Int[1, 2, 3]) == 8

        # lex rank / unrank round-trip via monomials
        for n in 0:5, k in 0:n
            Is = multi_indices(n, k)
            @test length(Is) == binomial(n, k)
            Λn = exterior_algebra(QQ, n)
            for (r, I) in enumerate(Is)
                @test coord_index(n, I) == ParametricLieAlgebras._degree_offset(n, k) + r
                α = exterior_monomial(Λn, I)
                @test coords(α)[coord_index(n, I)] == QQ(1)
                @test count(!iszero, coords(α)) == 1
            end
        end
    end

    @testset "unit / generators / from LieAlgebra" begin
        Λ = exterior_algebra(QQ, 2)
        @test isone(one(Λ))
        @test iszero(zero(Λ))
        @test one(Λ) + zero(Λ) == one(Λ)
        e1 = exterior_generator(Λ, 1)
        e2 = exterior_generator(Λ, 2)
        @test exterior_degree(e1) == 1
        @test is_homogeneous(e1)
        @test support_degrees(e1 + one(Λ)) == [0, 1]
        @test homogeneous_part(e1 + one(Λ), 0) == one(Λ)
        @test homogeneous_part(e1 + one(Λ), 1) == e1

        L = LieAlgebra(QQ, 4)
        ΛL = exterior_algebra(L)
        @test ambient_dim(ΛL) == dim(L)
        @test dim(ΛL) == 16
    end

    @testset "wedge: signs, associativity, graded commutativity" begin
        Λ = exterior_algebra(QQ, 4)
        e = [exterior_generator(Λ, i) for i in 1:4]

        # antisymmetry of generators
        @test wedge(e[1], e[2]) == -wedge(e[2], e[1])
        @test iszero(wedge(e[1], e[1]))
        @test exterior_monomial(Λ, [2, 1]) == -exterior_monomial(Λ, [1, 2])
        @test iszero(exterior_monomial(Λ, [1, 1]))

        # odd square vanishes; even need not
        α = e[1] + e[2]
        @test iszero(wedge(α, α))
        β = wedge(e[1], e[2])  # even
        @test wedge(β, β) == zero(Λ)  # 4-form e1∧e2∧e1∧e2 = 0 by repeat

        # associativity on generators
        lhs = wedge(wedge(e[1], e[2]), e[3])
        rhs = wedge(e[1], wedge(e[2], e[3]))
        @test lhs == rhs
        @test lhs == exterior_monomial(Λ, [1, 2, 3])

        # general associativity sample
        x = QQ(2) * e[1] + e[3]
        y = e[2] - QQ(3) * wedge(e[1], e[4])
        z = one(Λ) + e[4]
        @test wedge(wedge(x, y), z) == wedge(x, wedge(y, z))
        @test x * y == wedge(x, y)

        # graded commutativity for homogeneous elements
        a = wedge(e[1], e[2])  # deg 2
        b = e[3]               # deg 1
        @test wedge(a, b) == (-1)^(exterior_degree(a) * exterior_degree(b)) * wedge(b, a)

        # scalar linearity
        @test wedge(QQ(2) * e[1], e[2]) == QQ(2) * wedge(e[1], e[2])
        @test 3 * e[1] == QQ(3) * e[1]
    end

    @testset "interior product and form_eval" begin
        Λ = exterior_algebra(QQ, 3)
        e1, e2, e3 = exterior_generator(Λ, 1), exterior_generator(Λ, 2), exterior_generator(Λ, 3)
        ω = wedge(e1, e2)

        @test interior_product([1, 0, 0], ω) == e2
        @test interior_product([0, 1, 0], ω) == -e1
        @test iszero(interior_product([0, 0, 1], ω))
        @test interior_product([1, 0, 0], e1) == one(Λ)
        @test iszero(interior_product([1, 0, 0], one(Λ)))

        F = QQ
        v1 = [F(1), F(0), F(0)]
        v2 = [F(0), F(1), F(0)]
        @test form_eval(ω, [v1, v2]) == F(1)
        @test form_eval(ω, [v2, v1]) == F(-1)
        @test form_eval(ω, [v1, v1]) == F(0)
        vol = wedge(wedge(e1, e2), e3)
        @test form_eval(vol, [v1, v2, [F(0), F(0), F(1)]]) == F(1)
    end

    @testset "form_eval degree 0" begin
        Λ = exterior_algebra(QQ, 2)
        @test form_eval(one(Λ), Vector[]) == QQ(1)
        @test form_eval(zero(Λ), Vector[]) == QQ(0)
    end

    @testset "finite field" begin
        F5 = Nemo.GF(5)
        Λ = exterior_algebra(F5, 2)
        e1, e2 = exterior_generator(Λ, 1), exterior_generator(Λ, 2)
        @test iszero(wedge(e1, e1))
        @test wedge(e1, e2) == -wedge(e2, e1)
        @test ambient_dim(Λ) == 2
        @test dim(Λ) == 4
    end
end
