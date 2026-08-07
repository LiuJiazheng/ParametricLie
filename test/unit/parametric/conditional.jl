@testset "v0.2.3 conditional linear algebra" begin
    QQ = Nemo.QQ

    # -------------------------------------------------------------------------
    # Layer 1 — pure matrix, pivot a
    # -------------------------------------------------------------------------
    @testset "matrix pivot a" begin
        R, (a,) = AbstractAlgebra.polynomial_ring(QQ, [:a])
        K = AbstractAlgebra.fraction_field(R)
        M = AbstractAlgebra.matrix(K, [K(a) K(0); K(0) K(1)])

        res = conditional_rank(M)
        @test !res.incomplete
        @test length(res) == 2

        ranks = sort([leaf.rank for leaf in leaves(res)])
        @test ranks == [1, 2]

        branches = Set(c.branch for leaf in leaves(res) for c in leaf.trail)
        @test :zero in branches && :nonzero in branches
        polys = Set(string(c.poly) for leaf in leaves(res) for c in leaf.trail)
        @test "a" in polys

        for leaf in leaves(res)
            @test leaf.complete
            if any(c.branch === :nonzero for c in leaf.trail)
                @test leaf.rank == 2
                @test status(leaf.sigma, a) === PIVOT_NONZERO
            else
                @test leaf.rank == 1
                @test status(leaf.sigma, a) === PIVOT_ZERO
            end
        end

        # rref / solve smoke on the same matrix family
        rr = conditional_rref(M)
        @test length(rr) == 2
        @test all(leaf -> leaf.rref !== nothing || !leaf.complete, leaves(rr))

        A = AbstractAlgebra.matrix(K, [K(a) K(0); K(0) K(1)])
        b = AbstractAlgebra.matrix(K, 2, 1, [K(1), K(0)])
        sol = conditional_solve(A, b)
        @test !sol.incomplete
        # a≠0 → unique; a=0 → inconsistent (rhs first row becomes 1 with zero row)
        kinds = Set(s.solution[1] for s in leaves(sol) if s.complete && s.solution !== nothing)
        @test :unique in kinds
        @test :inconsistent in kinds
    end

    # -------------------------------------------------------------------------
    # Layer 2 — pure matrix, pivot a − b
    # -------------------------------------------------------------------------
    @testset "matrix pivot a-b" begin
        R, (a, b) = AbstractAlgebra.polynomial_ring(QQ, [:a, :b])
        K = AbstractAlgebra.fraction_field(R)
        M = AbstractAlgebra.matrix(K, [K(a - b) K(0); K(0) K(0)])

        res = conditional_rank(M)
        @test !res.incomplete
        @test length(res) == 2

        trail_polys = [normalize_poly(c.poly) for leaf in leaves(res) for c in leaf.trail]
        @test any(p -> p == normalize_poly(a - b), trail_polys)

        ranks = sort([leaf.rank for leaf in leaves(res)])
        @test ranks == [0, 1]

        for leaf in leaves(res)
            @test leaf.complete
            has_eq = any(c.branch === :zero for c in leaf.trail)
            has_nz = any(c.branch === :nonzero for c in leaf.trail)
            @test xor(has_eq, has_nz)
            if has_nz
                @test leaf.rank == 1
                @test status(leaf.sigma, a - b) === PIVOT_NONZERO
            else
                @test leaf.rank == 0
                @test status(leaf.sigma, a - b) === PIVOT_ZERO
            end
        end

        Σ = assume_nonzero(empty_assumptions(K), a)
        @test Σ !== nothing
        @test status(Σ, a * (a - b)) === PIVOT_UNKNOWN
        Σ0 = assume_zero(Σ, a * (a - b))
        @test Σ0 !== nothing
        @test status(Σ0, a - b) === PIVOT_ZERO
        @test status(Σ0, a) === PIVOT_NONZERO
        @test all(z -> z != normalize_poly(a), Σ0.zeros)
    end

    # -------------------------------------------------------------------------
    # Layer 3 — 2D Lie family + incremental refine
    # -------------------------------------------------------------------------
    @testset "2D family refine center/Killing/derived" begin
        R, (a,) = AbstractAlgebra.polynomial_ring(QQ, [:a])
        L = LieAlgebra(R, 2, Dict((1, 2) => [R(0), a]))
        @test check_jacobi(L).ok

        T0 = cond_tree(L)
        @test length(T0) == 1

        T = refine(T0, :center, :killing_radical, :derived_profile)
        @test !T.result.incomplete
        @test length(T) == 2

        dims = sort([center_dim(leaf) for leaf in leaves(T)])
        @test dims == [0, 2]

        for leaf in leaves(T)
            @test leaf.complete
            @test center_basis(leaf) !== nothing
            @test size(center_basis(leaf), 2) == center_dim(leaf)
            @test killing_rank_of(leaf) !== nothing
            @test derived_dim_of(leaf) !== nothing
            @test derived_profile_of(leaf) !== nothing
            @test is_solvable_of(leaf) === true
            if status(leaf.sigma, a) === PIVOT_NONZERO
                @test center_dim(leaf) == 0
                @test derived_dim_of(leaf) == 1
                @test killing_rank_of(leaf) == 1   # K = diag(a², 0)
            elseif status(leaf.sigma, a) === PIVOT_ZERO
                @test center_dim(leaf) == 2
                @test derived_dim_of(leaf) == 0
                @test killing_rank_of(leaf) == 0
            end
        end

        # nilpotency as separate pass inherits Σ from center
        Tn = refine(conditional_center(L), :nilpotency)
        for leaf in leaves(Tn)
            if status(leaf.sigma, a) === PIVOT_NONZERO
                @test is_nilpotent_of(leaf) === false
            elseif status(leaf.sigma, a) === PIVOT_ZERO
                @test is_nilpotent_of(leaf) === true
            end
        end

        # Der nice-to-have
        Td = refine(cond_tree(L), :derivations)
        @test !Td.result.incomplete
        @test all(leaf -> der_dim_of(leaf) !== nothing, leaves(Td))
    end

    # -------------------------------------------------------------------------
    # Layer 4 — 3D family with a and a−b
    # -------------------------------------------------------------------------
    @testset "3D family center (a and a-b)" begin
        R, (a, b) = AbstractAlgebra.polynomial_ring(QQ, [:a, :b])
        L = LieAlgebra(R, 3, Dict(
            (1, 2) => [R(0), a, R(0)],
            (1, 3) => [R(0), R(0), a - b],
        ))
        @test check_jacobi(L).ok

        T = conditional_invariants(L, :center, :killing_rank, :derived_dim, :radical)
        @test !T.result.incomplete
        @test length(T) >= 2

        dim_set = Set(center_dim(leaf) for leaf in leaves(T))
        @test 0 in dim_set
        @test 1 in dim_set
        @test 3 in dim_set

        found_generic = false
        found_ab = false
        for leaf in leaves(T)
            sa = status(leaf.sigma, a)
            sab = status(leaf.sigma, a - b)
            if sa === PIVOT_NONZERO && sab === PIVOT_NONZERO
                @test center_dim(leaf) == 0
                found_generic = true
            end
            if sa === PIVOT_ZERO && sab === PIVOT_ZERO
                @test center_dim(leaf) == 3
                found_ab = true
            end
            B = center_basis(leaf)
            @test B !== nothing
            @test size(B, 1) == 3
            @test size(B, 2) == center_dim(leaf)
            @test killing_rank_of(leaf) !== nothing
            @test derived_dim_of(leaf) !== nothing
            @test radical_dim_of(leaf) !== nothing
        end
        @test found_generic
        @test found_ab
    end
end
