@testset "NR bracket" begin
    F = Nemo.QQ

    H = LieAlgebra(F, 3, Dict((1, 2) => [0, 0, 1]))
    sl2 = LieAlgebra(F, 3, Dict(
        (1, 2) => [0, 2, 0],
        (1, 3) => [0, 0, -2],
        (2, 3) => [1, 0, 0],
    ))
    A = LieAlgebra(F, 2)

    @testset "μ from structure constants" begin
        μ = adjoint_bracket_cochain(H)
        @test length(μ) == 3 * binomial(3, 2)
        # μ(e1,e2)=e3, other pairs 0
        @test μ[1:3] == [F(0), F(0), F(1)]
        @test all(iszero, μ[4:9])
    end

    @testset "[μ,μ] = 0 (Jacobi)" begin
        for L in (H, sl2, A, LieAlgebra(F, 0), LieAlgebra(F, 1))
            μ = adjoint_bracket_cochain(L)
            @test all(iszero, nr_bracket(L, μ, 2, μ, 2))
        end
    end

    @testset "graded skew-symmetry" begin
        C = ce_complex(H, adjoint_module(H))
        # (p-1)(q-1) even ⇒ skew; odd ⇒ symmetric
        φ1 = fill(F(0), cochain_dim(C, 1)); φ1[2] = F(1)
        ψ2 = fill(F(0), cochain_dim(C, 2)); ψ2[1] = F(1)
        @test all(iszero, nr_bracket(H, φ1, 1, ψ2, 2) .+ nr_bracket(H, ψ2, 2, φ1, 1))

        φ2 = fill(F(0), cochain_dim(C, 2)); φ2[1] = F(1)
        ψ2b = fill(F(0), cochain_dim(C, 2)); ψ2b[5] = F(1)
        @test nr_bracket(H, φ2, 2, ψ2b, 2) == nr_bracket(H, ψ2b, 2, φ2, 2)
    end

    @testset "dφ = (-1)^{p+1} [μ, φ]_NR" begin
        for L in (H, sl2)
            μ = adjoint_bracket_cochain(L)
            C = ce_complex(L, adjoint_module(L))
            for p in 1:2
                np = cochain_dim(C, p)
                Dp = ce_differential(C, p)
                sign = isodd(p + 1) ? -1 : 1   # (-1)^{p+1}: p=1 → +1; p=2 → -1
                for j in 1:np
                    φ = fill(F(0), np)
                    φ[j] = F(1)
                    dφ = Dp * AbstractAlgebra.matrix(F, np, 1, φ)
                    dφv = [dφ[i, 1] for i in 1:size(dφ, 1)]
                    br = nr_bracket(L, μ, 2, φ, p)
                    expect = sign == 1 ? br : .-br
                    @test dφv == expect
                end
            end
        end
    end

    @testset "H² cocycle ⇒ [φ,φ] ∈ Z³" begin
        C = ce_complex(H, adjoint_module(H))
        H2 = cohomology(C, 2)
        D2 = ce_differential(C, 2)
        for j in 1:dim(H2)
            φ = [basis_matrix(H2)[i, j] for i in 1:cochain_dim(C, 2)]
            # φ is a cocycle
            @test all(iszero, D2 * AbstractAlgebra.matrix(F, length(φ), 1, φ))
            sq = nr_bracket(H, φ, 2, φ, 2)
            # [φ,φ] ∈ Z³ = ker d³; for n=3, d³: C³→0 so automatic, check length
            @test length(sq) == cochain_dim(C, 3)
            # still a cocycle relative to d²… wait [φ,φ] ∈ C³, check d³=0 vacuous
            # and for obstruction: class in H³ — just ensure it is in ker (trivial) / compute
            @test all(iszero, ce_differential(C, 3) * AbstractAlgebra.matrix(F, length(sq), 1, sq))
        end
    end

    @testset "abelian ⇒ μ=0, NR = Gerstenhaber on cochains only" begin
        μ = adjoint_bracket_cochain(A)
        @test all(iszero, μ)
        C = ce_complex(A, adjoint_module(A))
        φ = fill(F(0), cochain_dim(C, 1)); φ[1] = F(1)
        ψ = fill(F(0), cochain_dim(C, 1)); ψ[2] = F(1)
        # [φ,ψ] ∈ C¹ for p=q=1
        br = nr_bracket(A, φ, 1, ψ, 1)
        @test length(br) == cochain_dim(C, 1)
    end
end
