@testset "Maurer–Cartan truncation" begin
    F = Nemo.QQ

    H = LieAlgebra(F, 3, Dict((1, 2) => [0, 0, 1]))
    sl2 = LieAlgebra(F, 3, Dict(
        (1, 2) => [0, 2, 0],
        (1, 3) => [0, 0, -2],
        (2, 3) => [1, 0, 0],
    ))

    @testset "rigid ⇒ only zero seed / empty H² list" begin
        @test is_rigid(sl2)
        @test isempty(formal_deformation(sl2; order = 2))
        φ0 = fill(F(0), 3 * binomial(3, 2))
        D0 = formal_deformation(sl2, φ0; order = 3)
        @test is_integrable(D0)
        @test filled_order(D0) == 3
        @test all(iszero, deformation_term(D0, 1))
        @test all(iszero, deformation_term(D0, 2))
        @test all(iszero, deformation_term(D0, 3))
    end

    @testset "seed must be a cocycle" begin
        bad = fill(F(0), 9)
        # find a non-cocycle if any
        C = ce_complex(H, adjoint_module(H))
        Z2 = cocycles(C, 2)
        found = false
        for j in 1:9
            v = fill(F(0), 9); v[j] = F(1)
            if !ParametricLieAlgebras._in_column_span(Z2, v)
                @test_throws ArgumentError formal_deformation(H, v; order = 2)
                found = true
                break
            end
        end
        # Heisenberg: dim Z² = 8, dim C² = 9 ⇒ exists
        @test found
    end

    @testset "order-2 slice + cache" begin
        H2 = cohomology(H, adjoint_module(H), 2)
        @test dim(H2) == 5
        φ1 = [basis_matrix(H2)[i, 1] for i in 1:9]
        D = formal_deformation(H, φ1; order = 2)
        @test D isa FormalDeformation
        @test max_order(D) == 2
        @test filled_order(D) >= 1
        @test deformation_term(D, 1) == φ1
        # RHS / φ₂ cached after order-2 attempt
        ψ2 = mc_rhs(D, 2)
        @test length(ψ2) == 3 * binomial(3, 3)
        cert = mc_certificate(D)
        @test cert isa MCCertificate
        @test is_integrable(cert) == is_integrable(D)
        @test deformation_term(cert, 1) === deformation_term(D, 1)
        if is_integrable(D)
            @test length(deformation_term(D, 2)) == 9
            @test stalled_at(D) === nothing
            @test obstruction_cochain(D) === nothing
            # d φ₂ = ψ₂
            D2 = ce_differential(ce_complex(D), 2)
            φ2 = deformation_term(D, 2)
            lhs = D2 * AbstractAlgebra.matrix(F, 9, 1, φ2)
            @test [lhs[i, 1] for i in 1:3] == ψ2
        else
            @test stalled_at(D) == 2
            @test obstruction_cochain(D) == ψ2
        end
    end

    @testset "all H² directions to order 2" begin
        dirs = formal_deformation(H; order = 2)
        @test length(dirs) == 5
        for D in dirs
            @test filled_order(D) >= 1
            @test haskey(D.phi, 1)
            @test haskey(D.rhs, 2) || stalled_at(D) == 2 || is_integrable(D)
        end
    end

    @testset "extend! reuses cache" begin
        H2 = cohomology(H, adjoint_module(H), 2)
        φ1 = [basis_matrix(H2)[i, 1] for i in 1:9]
        D = formal_deformation(H, φ1; order = 2)
        φ1_id = deformation_term(D, 1)
        extend!(D, 3)
        @test max_order(D) == 3
        @test deformation_term(D, 1) === φ1_id  # same cached object
        if is_integrable(D)
            @test filled_order(D) == 3
            @test haskey(D.phi, 3)
            @test haskey(D.rhs, 3)
        else
            @test stalled_at(D) !== nothing
            @test stalled_at(D) <= 3
            @test obstruction_cochain(D) !== nothing
        end
    end

    @testset "zero deformation integrates to any order" begin
        φ0 = fill(F(0), 9)
        D = formal_deformation(H, φ0; order = 4)
        @test is_integrable(D)
        for k in 1:4
            @test all(iszero, deformation_term(D, k))
        end
        for k in 2:4
            @test all(iszero, mc_rhs(D, k))
        end
    end
end
