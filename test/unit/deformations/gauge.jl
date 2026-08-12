@testset "gauge equivalence" begin
    F = Nemo.QQ
    H = LieAlgebra(F, 3, Dict((1, 2) => [0, 0, 1]))
    C = ce_complex(H, adjoint_module(H))
    _, W2 = ParametricLieAlgebras._c2_splitting(C)

    @testset "order-1: same class via dα" begin
        H2 = cohomology(H, adjoint_module(H), 2)
        φ1 = [basis_matrix(H2)[i, 1] for i in 1:9]
        α = fill(F(0), cochain_dim(C, 1))
        α[1] = F(1)
        dα = ce_differential(C, 1) * AbstractAlgebra.matrix(F, length(α), 1, α)
        dαv = [dα[i, 1] for i in 1:9]
        φ1p = φ1 .+ dαv
        D = formal_deformation(H, φ1; order = 1)
        Dp = formal_deformation(H, φ1p; order = 1)
        @test equivalent(D, Dp; order = 1)
        @test equivalent(Dp, D; order = 1)
        w = equivalent_with_gauge(D, Dp; order = 1)
        @test w.ok
        @test length(w.alphas) == 1
    end

    @testset "order-1: distinct H² directions inequivalent" begin
        dirs = formal_deformation(H; order = 1)
        @test length(dirs) >= 2
        @test equivalent(dirs[1], dirs[1]; order = 1)
        @test !equivalent(dirs[1], dirs[2]; order = 1)
    end

    @testset "higher order: self-equivalent; extend uses cache" begin
        H2 = cohomology(H, adjoint_module(H), 2)
        φ1 = [basis_matrix(H2)[i, 1] for i in 1:9]
        D = formal_deformation(H, φ1; order = 2)
        if is_integrable(D)
            @test equivalent(D, D; order = 2)
            extend!(D, 3)
            if is_integrable(D)
                @test equivalent(D, D; order = 3)
            end
            α = fill(F(0), 9); α[2] = F(1)
            dα = ce_differential(C, 1) * AbstractAlgebra.matrix(F, 9, 1, α)
            φ1p = φ1 .+ [dα[i, 1] for i in 1:9]
            Dp = formal_deformation(H, φ1p; order = 1)
            @test equivalent(D, Dp; order = 1)
        end
    end

    @testset "gauge_normal_form: orbit section" begin
        H2 = cohomology(H, adjoint_module(H), 2)
        φ1 = [basis_matrix(H2)[i, 1] for i in 1:9]
        α = fill(F(0), cochain_dim(C, 1))
        α[1] = F(1)
        dα = ce_differential(C, 1) * AbstractAlgebra.matrix(F, length(α), 1, α)
        φ1p = φ1 .+ [dα[i, 1] for i in 1:9]

        D = formal_deformation(H, φ1; order = 1)
        Dp = formal_deformation(H, φ1p; order = 1)
        Dn = gauge_normal_form(D; order = 1)
        Dpn = gauge_normal_form(Dp; order = 1)
        @test deformation_term(Dn, 1) == deformation_term(Dpn, 1)
        @test ParametricLieAlgebras._in_column_span(W2, deformation_term(Dn, 1))
        @test equivalent(D, Dn; order = 1)
        @test equivalent(Dp, Dn; order = 1)

        # idempotent
        Dn2 = gauge_normal_form(Dn; order = 1)
        @test deformation_term(Dn2, 1) == deformation_term(Dn, 1)

        # distinct H² directions stay distinct after NF
        dirs = formal_deformation(H; order = 1)
        n1 = gauge_normal_form(dirs[1]; order = 1)
        n2 = gauge_normal_form(dirs[2]; order = 1)
        @test deformation_term(n1, 1) != deformation_term(n2, 1)
        @test !equivalent(n1, n2; order = 1)

        # with_gauge witness + in-place
        g = gauge_normal_form_with_gauge(Dp; order = 1)
        @test deformation_term(g.deformation, 1) == deformation_term(Dn, 1)
        @test length(g.alphas) == 1
        Dc = formal_deformation(H, φ1p; order = 1)
        gauge_normal_form!(Dc; order = 1)
        @test deformation_term(Dc, 1) == deformation_term(Dn, 1)
    end

    @testset "gauge_normal_form higher order" begin
        H2 = cohomology(H, adjoint_module(H), 2)
        φ1 = [basis_matrix(H2)[i, 1] for i in 1:9]
        D = formal_deformation(H, φ1; order = 2)
        if is_integrable(D)
            Dn = gauge_normal_form(D; order = 2)
            @test equivalent(D, Dn; order = 2)
            @test ParametricLieAlgebras._in_column_span(W2, deformation_term(Dn, 1))
            @test ParametricLieAlgebras._in_column_span(W2, deformation_term(Dn, 2))
            # re-normalize is fixed point
            Dn2 = gauge_normal_form(Dn; order = 2)
            @test deformation_term(Dn2, 1) == deformation_term(Dn, 1)
            @test deformation_term(Dn2, 2) == deformation_term(Dn, 2)
        end
    end
end
