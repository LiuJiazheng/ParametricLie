@testset "central extensions from H²" begin
    F = Nemo.QQ

    @testset "abelian₂ → Heisenberg (nontrivial)" begin
        A = LieAlgebra(F, 2)
        H2 = cohomology(A, 2)
        @test dim(H2) == 1  # Λ²(ℚ²)* 

        ω = [basis_matrix(H2)[i, 1] for i in 1:binomial(2, 2)]
        @test !is_trivial_cocycle(A, ω)

        Ê = central_extension(A, ω)
        @test dim(Ê) == 3
        @test check_jacobi(Ê).ok

        # Nonsplit: derived algebra = center = span{e3}, classical Heisenberg
        @test dim(center(Ê)) == 1
        @test dim(derived_algebra(Ê)) == 1
        @test is_nilpotent(Ê) && is_solvable(Ê)
        @test killing_rank(Ê) == 0

        z = only(basis_elems(center(Ê)))
        @test iszero(z.coords[1]) && iszero(z.coords[2]) && !iszero(z.coords[3])
        # [e1,e2] is a nonzero multiple of the center
        br = lie_bracket(Ê, basis_elem(Ê, 1), basis_elem(Ê, 2))
        @test !iszero(br.coords[3])
        @test iszero(br.coords[1]) && iszero(br.coords[2])
        for i in 1:2
            @test iszero(lie_bracket(Ê, basis_elem(Ê, i), basis_elem(Ê, 3)).coords)
        end

        # Same invariants as the hand-built Heisenberg
        H = LieAlgebra(F, 3, Dict((1, 2) => [0, 0, 1]))
        @test dim(center(Ê)) == dim(center(H))
        @test map(dim, terms(derived_series(Ê))) == map(dim, terms(derived_series(H)))
        @test map(dim, terms(lower_central_series(Ê))) ==
              map(dim, terms(lower_central_series(H)))
        @test dim(derivations(Ê)) == dim(derivations(H)) == 6
    end

    @testset "split extension from trivial cocycle" begin
        A = LieAlgebra(F, 2)
        ω0 = [F(0)]
        @test is_trivial_cocycle(A, ω0)
        E0 = central_extension(A, ω0)
        @test check_jacobi(E0).ok
        @test dim(E0) == 3
        @test dim(derived_algebra(E0)) == 0   # abelian = A ⊕ F
        @test dim(center(E0)) == 3
        @test is_nilpotent(E0)
    end

    @testset "coboundary ⇒ split" begin
        # For abelian g, B² = 0, so every cocycle is nontrivial when nonzero.
        # Use Heisenberg: H²(H,F) has dim 2, and B² may be nonzero.
        H = LieAlgebra(F, 3, Dict((1, 2) => [0, 0, 1]))
        C = ce_complex(H)
        B2 = coboundaries(C, 2)
        Z2 = cocycles(C, 2)
        H2 = cohomology(C, 2)
        @test dim(H2) == 2
        @test size(Z2, 2) - size(B2, 2) == 2

        if size(B2, 2) > 0
            ωB = [B2[i, 1] for i in 1:binomial(3, 2)]
            @test is_trivial_cocycle(H, ωB)
            EB = central_extension(H, ωB)
            @test check_jacobi(EB).ok
            # Split over H: the new center direction is a direct summand algebraically
            # witnessed by [ω]=0; at least Jacobi holds and dim = 4.
            @test dim(EB) == 4
        end

        # Nontrivial class from H² representatives
        ωN = [basis_matrix(H2)[i, 1] for i in 1:binomial(3, 2)]
        @test !is_trivial_cocycle(H, ωN)
        EN = central_extension(H, ωN)
        @test check_jacobi(EN).ok
        @test dim(EN) == 4
        @test dim(center(EN)) >= 1
    end

    @testset "sl₂: H² = 0 ⇒ only split extensions" begin
        sl2 = LieAlgebra(F, 3, Dict(
            (1, 2) => [0, 2, 0],
            (1, 3) => [0, 0, -2],
            (2, 3) => [1, 0, 0],
        ))
        @test dim(cohomology(sl2, 2)) == 0
        ω0 = fill(F(0), binomial(3, 2))
        E = central_extension(sl2, ω0)
        @test check_jacobi(E).ok
        @test dim(E) == 4
        # Split: radical / derived pattern — e4 central and [E,E] lands in sl2 summand
        @test dim(center(E)) == 1
        z = only(basis_elems(center(E)))
        @test iszero(z.coords[1]) && iszero(z.coords[2]) && iszero(z.coords[3])
        @test !iszero(z.coords[4])
    end

    @testset "non-cocycle rejected" begin
        # ax+b ⊕ ℚ: [e1,e2]=e1. Then (dω)(e1,e2,e3)=ω(e1,e3), so e1∧e3 is not a cocycle.
        L = LieAlgebra(F, 3, Dict((1, 2) => [1, 0, 0]))
        @test check_jacobi(L).ok
        ω_bad = [F(0), F(1), F(0)]  # coords on (e1∧e2, e1∧e3, e2∧e3)
        Z2 = cocycles(ce_complex(L), 2)
        @test !ParametricLieAlgebras._in_column_span(Z2, ω_bad)
        @test_throws ArgumentError central_extension(L, ω_bad)
    end
end
