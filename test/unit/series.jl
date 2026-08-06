@testset "series / commutator_span" begin
    F = Nemo.QQ

    # --- Cache order: is_solvable → layers → terms (Heisenberg) ---
    H = LieAlgebra(F, 3, Dict((1, 2) => [0, 0, 1]))
    @test check_jacobi(H).ok

    D = derived_series(H)
    @test isempty(D.cache) && D.layers === nothing   # still lazy
    @test is_solvable(D)                             # first call triggers eval
    @test !isempty(D.cache) && D.layers == 3 && D.reached_zero
    cached = D.cache
    @test layers(D) == 3                             # reads cache
    @test D.cache === cached                         # no recompute / replace
    @test map(dim, terms(D)) == [3, 1, 0]            # full cached sequence
    @test terms(D) === D.cache

    C = lower_central_series(H)
    @test isempty(C.cache) && C.layers === nothing
    @test is_nilpotent(C)                            # first call triggers eval
    @test C.layers == 3 && C.reached_zero
    cached_c = C.cache
    @test layers(C) == 3
    @test C.cache === cached_c
    @test map(dim, terms(C)) == [3, 1, 0]
    @test terms(C) === C.cache

    # --- Same order on sl₂: !is_solvable first, then Inf layers ---
    sl2 = LieAlgebra(F, 3, Dict(
        (1, 2) => [0, 2, 0],
        (1, 3) => [0, 0, -2],
        (2, 3) => [1, 0, 0],
    ))
    @test check_jacobi(sl2).ok
    Dsl = derived_series(sl2)
    @test isempty(Dsl.cache) && Dsl.layers === nothing
    @test !is_solvable(Dsl)                          # first call triggers eval
    @test Dsl.layers == Inf && Dsl.stabilized && !Dsl.reached_zero
    cached_sl = Dsl.cache
    @test layers(Dsl) == Inf
    @test Dsl.cache === cached_sl
    @test dim(terms(Dsl)[end]) == 3
    @test terms(Dsl) === Dsl.cache

    # --- Abelian: both series hit 0 at step 2; solvable + nilpotent ---
    A = LieAlgebra(F, 3)
    DA = derived_series(A)
    @test DA isa LieSeries
    @test isempty(DA.cache)          # lazy: nothing computed yet
    @test DA.layers === nothing
    dterms = terms(DA, 5)
    @test length(dterms) == 2       # g, then 0 (shorter than 5)
    @test dim(dterms[1]) == 3
    @test dim(dterms[2]) == 0
    @test DA.layers == 2             # finalized while filling terms
    @test layers(DA) == 2
    @test is_solvable(DA)
    @test is_solvable(A)
    Cser = lower_central_series(A)
    @test layers(Cser) == 2
    @test is_nilpotent(Cser)
    @test is_nilpotent(A)
    @test dim(derived_algebra(A)) == 0

    cterms = terms(Cser)            # full cached sequence
    @test length(cterms) == 2
    @test dim(cterms[1]) == 3 && dim(cterms[2]) == 0

    # --- Heisenberg: derived algebra certificate ---
    Hg = derived_algebra(H)
    @test dim(Hg) == 1
    z = only(basis_elems(Hg))
    @test iszero(z.coords[1]) && iszero(z.coords[2]) && !iszero(z.coords[3])
    @test is_solvable(H)
    @test is_nilpotent(H)

    # commutator_span reuse: [g, span e3] = 0
    Z = center(H)
    @test dim(commutator_span(H, full_space(H), Z)) == 0

    # --- sl₂ algebra-level flags ---
    @test dim(derived_algebra(sl2)) == 3
    @test !is_solvable(sl2)
    @test !is_nilpotent(sl2)
    Csl = lower_central_series(sl2)
    @test layers(Csl) == Inf
    @test !is_nilpotent(Csl)

    # --- change-of-basis: dimension sequences invariant ---
    P = AbstractAlgebra.matrix(F, 3, 3, [1, 1, 0, 0, 1, 1, 0, 0, 1])
    Hp = change_of_basis(H, P)
    @test map(dim, terms(derived_series(Hp))) == [3, 1, 0]
    @test map(dim, terms(lower_central_series(Hp))) == [3, 1, 0]
    @test layers(derived_series(Hp)) == 3
    @test layers(lower_central_series(Hp)) == 3
end
