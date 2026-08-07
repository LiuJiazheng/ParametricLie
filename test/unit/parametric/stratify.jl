@testset "v0.2.4–0.2.6 conditional analysis / stratify / validate" begin
    QQ = Nemo.QQ
    R, (a, b) = AbstractAlgebra.polynomial_ring(QQ, [:a, :b])
    # Running example: [e1,e2]=a e2, [e1,e3]=(a-b) e3
    L = LieAlgebra(R, 3, Dict(
        (1, 2) => [R(0), a, R(0)],
        (1, 3) => [R(0), R(0), a - b],
    ))
    @test check_jacobi(L).ok

    suite = [:center, :derived_dim, :solvability, :nilpotency]

    @testset "v0.2.4 analyze_conditional" begin
        T = analyze_conditional(L; invariants = suite)
        @test !T.result.incomplete
        @test is_complete(T; invariants = suite)
        @test length(T) >= 3

        # Inheritance: refining center again does not grow the tree
        n0 = length(T)
        T2 = refine(T, :center)
        @test length(T2) == n0
        @test all(leaf -> isempty(unresolved_invariants(leaf, [:center])), leaves(T2))

        dims = Set(center_dim(leaf) for leaf in leaves(T))
        @test 0 in dims
        @test 1 in dims
        @test 3 in dims

        # Generic-ish leaf: a≠0 and a-b≠0 → center 0, derived 2
        found = false
        for leaf in leaves(T)
            if status(leaf.sigma, a) === PIVOT_NONZERO &&
               status(leaf.sigma, a - b) === PIVOT_NONZERO
                @test center_dim(leaf) == 0
                @test derived_dim_of(leaf) == 2
                @test is_solvable_of(leaf) === true
                @test is_nilpotent_of(leaf) === false
                found = true
            end
        end
        @test found
    end

    @testset "v0.2.5 stratify + jumps" begin
        S = stratify(L; invariants = suite)
        @test S.generic !== nothing
        @test length(S.strata) >= 3
        @test S.generic.signature[:center_dim] == 0
        @test S.generic.signature[:derived_dim] == 2

        jt = jump_table(S)
        @test !isempty(jt)
        # At least one jump in center_dim
        @test any(j -> any(e -> e.invariant === :center_dim, j.changed), jt)

        # Abelian corner should appear as a confirmed nilpotency jump
        found_ab = false
        for j in jt
            if get(j.target.signature, :center_dim, nothing) == 3
                @test get(j.target.signature, :derived_dim, nothing) == 0
                @test get(j.target.signature, :is_nilpotent, nothing) === true
                found_ab = true
            end
        end
        @test found_ab
    end

    @testset "v0.2.6 compare + validate fibers" begin
        S = stratify(L; invariants = suite)
        # Find indices for generic / a=b=0 / a=0 b≠0 style leaves
        igen = findfirst(st -> st.sigma == S.generic.sigma, S.strata)
        @test igen !== nothing

        i_ab = findfirst(st -> get(st.signature, :center_dim, -1) == 3, S.strata)
        @test i_ab !== nothing
        C = compare(S.strata[igen], S.strata[i_ab])
        @test C.same_signature === false
        @test any(e -> e.invariant === :center_dim && e.special_value == 3, C.changed)
        @test any(e -> e.invariant === :is_nilpotent, C.changed)

        # Fiber validations
        v_gen = validate_stratum(L, S.generic, Dict(a => QQ(1), b => QQ(2)))
        @test v_gen.ok

        # exceptional a=0, b=1
        i_a0 = findfirst(S.strata) do st
            get(st.signature, :center_dim, -1) == 1 &&
                status(st.sigma, a) === PIVOT_ZERO &&
                status(st.sigma, a - b) === PIVOT_NONZERO
        end
        if i_a0 !== nothing
            v0 = validate_stratum(L, S.strata[i_a0], Dict(a => QQ(0), b => QQ(1)))
            @test v0.ok
        end

        # a=b≠0
        i_abeq = findfirst(S.strata) do st
            get(st.signature, :center_dim, -1) == 1 &&
                status(st.sigma, a) === PIVOT_NONZERO &&
                status(st.sigma, a - b) === PIVOT_ZERO
        end
        if i_abeq !== nothing
            v1 = validate_stratum(L, S.strata[i_abeq], Dict(a => QQ(1), b => QQ(1)))
            @test v1.ok
        end

        v_ab = validate_stratum(L, S.strata[i_ab], Dict(a => QQ(0), b => QQ(0)))
        @test v_ab.ok
        @test dim(center(v_ab.fiber_report)) == 3
        @test is_nilpotent(v_ab.fiber_report)
    end
end
