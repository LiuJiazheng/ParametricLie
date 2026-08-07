@testset "analyze report" begin
    F = Nemo.QQ

    # --- Heisenberg: solvable, Levi empty, Der dim 6 ---
    H = LieAlgebra(F, 3, Dict((1, 2) => [0, 0, 1]))
    rH = analyze(H)
    @test rH isa LieAlgebraReport
    @test dim(rH) == 3
    @test jacobi(rH).ok
    @test dim(center(rH)) == 1
    @test is_solvable(rH) && is_nilpotent(rH)
    @test killing_rank(rH) == 0
    @test !is_semisimple(rH)
    @test dim(radical(rH)) == 3
    @test levi_kind(rH) === :empty
    @test simple_factor_dims(rH) == Int[]
    @test dim(derivations(rH)) == 6
    # Cached objects are the real certificates
    @test levi_decomposition(rH) isa LeviDecomposition
    @test dim(levi_decomposition(rH).levi) == 0
    @test derivations(rH) === rH.derivations
    @test killing_form(rH) === rH.killing
    @test map(dim, terms(derived_series(rH))) == [3, 1, 0]
    summary_H = sprint(show, MIME("text/plain"), rH)
    @test occursin("Levi:             dim 0", summary_H)
    @test occursin("Der:              dim 6", summary_H)
    @test occursin("simple factors:   dims []", summary_H)

    # --- sl₂: simple semisimple ---
    sl2 = LieAlgebra(F, 3, Dict(
        (1, 2) => [0, 2, 0],
        (1, 3) => [0, 0, -2],
        (2, 3) => [1, 0, 0],
    ))
    rS = analyze(sl2)
    @test jacobi(rS).ok
    @test is_semisimple(rS) && is_simple(rS)
    @test !is_solvable(rS)
    @test dim(radical(rS)) == 0
    @test levi_kind(rS) === :simple
    @test simple_factor_dims(rS) == [3]
    @test dim(derivations(rS)) == 3
    @test dim(levi_decomposition(rS).levi) == 3
    @test length(ideal_decomposition(rS)) == 1
    summary_S = sprint(show, MIME("text/plain"), rS)
    @test occursin("Levi:             dim 3, simple", summary_S)
    @test occursin("dims [3]", summary_S)

    # --- sl₂ ⊕ sl₂: composite ---
    sl2sl2 = LieAlgebra(F, 6, Dict(
        (1, 2) => [0, 2, 0, 0, 0, 0],
        (1, 3) => [0, 0, -2, 0, 0, 0],
        (2, 3) => [1, 0, 0, 0, 0, 0],
        (4, 5) => [0, 0, 0, 0, 2, 0],
        (4, 6) => [0, 0, 0, 0, 0, -2],
        (5, 6) => [0, 0, 0, 1, 0, 0],
    ))
    r2 = analyze(sl2sl2)
    @test is_semisimple(r2) && !is_simple(r2)
    @test levi_kind(r2) === :composite
    @test simple_factor_dims(r2) == [3, 3]
    @test length(ideal_decomposition(r2)) == 2
    summary_2 = sprint(show, MIME("text/plain"), r2)
    @test occursin("composite", summary_2)
    @test occursin("dims [3, 3]", summary_2)

    # --- sl₂ ⊕ Heisenberg: nontrivial radical + simple Levi ---
    sl2H = LieAlgebra(F, 6, Dict(
        (1, 2) => [0, 2, 0, 0, 0, 0],
        (1, 3) => [0, 0, -2, 0, 0, 0],
        (2, 3) => [1, 0, 0, 0, 0, 0],
        (4, 5) => [0, 0, 0, 0, 0, 1],
    ))
    rLH = analyze(sl2H)
    @test !is_semisimple(rLH)
    @test dim(radical(rLH)) == 3
    @test levi_kind(rLH) === :simple
    @test simple_factor_dims(rLH) == [3]
    @test dim(levi_decomposition(rLH).levi) == 3
    @test derivations(rLH) isa Derivations

    # --- abelian ---
    A = LieAlgebra(F, 2)
    rA = analyze(A)
    @test is_solvable(rA) && is_nilpotent(rA)
    @test levi_kind(rA) === :empty
    @test dim(derivations(rA)) == 4
end
