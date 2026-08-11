@testset "CE cohomology" begin
    F = Nemo.QQ

    sl2 = LieAlgebra(F, 3, Dict(
        (1, 2) => [0, 2, 0],
        (1, 3) => [0, 0, -2],
        (2, 3) => [1, 0, 0],
    ))
    H3 = LieAlgebra(F, 3, Dict((1, 2) => [0, 0, 1]))

    @testset "modules" begin
        A = LieAlgebra(F, 2)
        T = trivial_module(A)
        @test T isa LieModule
        @test dim(T) == 1
        @test parent(T) === A
        @test all(iszero, action_matrices(T)[1])
        @test act(T, 1, [F(7)]) == [F(0)]

        Ad = adjoint_module(sl2)
        @test dim(Ad) == 3
        # e1 · e2 = [H,X] = 2X
        v = act(Ad, 1, [F(0), F(1), F(0)])
        @test v == [F(0), F(2), F(0)]
    end

    @testset "lazy prefix cache" begin
        C = ce_complex(H3)
        @test C isa CEComplex
        @test C.filled == -1
        H2 = cohomology(C, 2)
        @test H2 isa CohomologyGroup
        @test H2.degree == 2
        @test C.filled == 2
        @test C.Z_cache[1] !== nothing  # Z^0
        @test C.Z_cache[2] !== nothing  # Z^1
        @test C.Z_cache[3] !== nothing  # Z^2
        @test C.B_cache[3] !== nothing
        @test C.H_cache[3] !== nothing
        @test C.Z_cache[4] === nothing  # Z^3 not yet
        cohomology(C, 3)
        @test C.filled == 3
        @test C.Z_cache[4] !== nothing
    end

    @testset "abelian trivial: d = 0, H^k = Λ^k" begin
        A = LieAlgebra(F, 3)
        C = ce_complex(A)
        for k in 0:3
            @test cochain_dim(C, k) == binomial(3, k)
            @test all(iszero, ce_differential(C, k))
            @test size(cocycles(C, k), 2) == binomial(3, k)
            @test size(coboundaries(C, k), 2) == 0
            @test dim(cohomology(C, k)) == binomial(3, k)
        end
        @test dim(cohomology(C, 4)) == 0
        A0 = LieAlgebra(F, 0)
        @test dim(cohomology(A0, 0)) == 1
        @test dim(cohomology(A0, 1)) == 0
    end

    @testset "Heisenberg trivial Betti 1,2,2,1" begin
        C = ce_complex(H3)
        @test dim(cohomology(C, 0)) == 1
        @test dim(cohomology(C, 1)) == 2  # (g/[g,g])*
        @test dim(cohomology(C, 2)) == 2
        @test dim(cohomology(C, 3)) == 1
        @test dim(derived_algebra(H3)) == 1
        @test dim(cohomology(H3, 1)) == dim(H3) - dim(derived_algebra(H3))
    end

    @testset "sl2 trivial Betti 1,0,0,1" begin
        C = ce_complex(sl2, trivial_module(sl2))
        @test dim(cohomology(C, 0)) == 1
        @test dim(cohomology(C, 1)) == 0
        @test dim(cohomology(C, 2)) == 0
        @test dim(cohomology(C, 3)) == 1
    end

    @testset "adjoint: H^0 = Z, H^1 = Out" begin
        @test dim(cohomology(H3, adjoint_module(H3), 0)) == dim(center(H3)) == 1
        @test dim(cohomology(sl2, adjoint_module(sl2), 0)) == 0
        @test dim(cohomology(sl2, adjoint_module(sl2), 1)) == 0
        @test dim(cohomology(sl2, adjoint_module(sl2), 2)) == 0
        # Heisenberg: Der 6, Inn = 3 − 1 = 2 ⇒ Out = 4
        @test dim(derivations(H3)) == 6
        @test dim(cohomology(H3, adjoint_module(H3), 1)) == 4
    end

    @testset "d² = 0 and dim H = dim Z − dim B" begin
        for L in (LieAlgebra(F, 2), H3, sl2)
            for M in (trivial_module(L), adjoint_module(L))
                C = ce_complex(L, M)
                n = dim(L)
                cohomology(C, n)  # fill all
                for k in 0:n-1
                    D0 = ce_differential(C, k)
                    D1 = ce_differential(C, k + 1)
                    if size(D0, 2) > 0 && size(D1, 1) > 0
                        @test all(iszero, D1 * D0)
                    end
                end
                chiC = 0
                chiH = 0
                for k in 0:n
                    Zk = cocycles(C, k)
                    Bk = coboundaries(C, k)
                    Hk = cohomology(C, k)
                    @test dim(Hk) == size(Zk, 2) - size(Bk, 2)
                    @test size(basis_matrix(Hk), 1) == cochain_dim(C, k)
                    # B ⊆ Z: d^k ∘ B^k = 0
                    if size(Bk, 2) > 0 && size(ce_differential(C, k), 1) > 0
                        @test all(iszero, ce_differential(C, k) * Bk)
                    end
                    chiC += (isodd(k) ? -1 : 1) * cochain_dim(C, k)
                    chiH += (isodd(k) ? -1 : 1) * dim(Hk)
                end
                @test chiC == chiH
            end
        end
    end

    @testset "change-of-basis Betti invariant" begin
        P = AbstractAlgebra.matrix(F, 3, 3, [1, 1, 0, 0, 1, 1, 0, 0, 1])
        for L in (H3, sl2)
            Lp = change_of_basis(L, P)
            for k in 0:3
                @test dim(cohomology(Lp, k)) == dim(cohomology(L, k))
                @test dim(cohomology(Lp, adjoint_module(Lp), k)) ==
                      dim(cohomology(L, adjoint_module(L), k))
            end
        end
    end
end
