@testset "v0.2.1 parametric coefficient domains" begin
    R, (a, b) = AbstractAlgebra.polynomial_ring(Nemo.QQ, [:a, :b])

    # --- identity over QQ[a,b]: parametric Heisenberg ---
    H = LieAlgebra(R, 3, Dict((1, 2) => [R(0), R(0), a]))
    @test parameters(H) == collect(AbstractAlgebra.gens(R))
    @test isempty(domain_denominators(H))
    jacH = check_jacobi(H)
    @test jacH.ok
    @test isempty(jacH.residuals)
    @test sprint(show, jacH) == "JacobiCertificate(identity)"

    e1, e2 = basis_elem(H, 1), basis_elem(H, 2)
    br = lie_bracket(H, e1, e2)
    @test iszero(br.coords[1]) && iszero(br.coords[2]) && br.coords[3] == a

    # --- nonidentity: explicit polynomial residual (no root-finding) ---
    Lbad = LieAlgebra(R, 3, Dict(
        (1, 2) => [a, R(0), R(0)],
        (1, 3) => [R(0), b, R(0)],
        (2, 3) => [R(0), R(0), R(1)],
    ))
    jac_bad = check_jacobi(Lbad)
    @test !jac_bad.ok
    @test !isempty(jac_bad.residuals)
    @test occursin("nonidentity", sprint(show, jac_bad))
    r0 = first(jac_bad.residuals)
    @test r0.triple == (1, 2, 3) || r0.triple isa Tuple{Int,Int,Int}
    @test any(!iszero, r0.residual)
    @test all(t -> parent(t) === R, r0.residual)
    @test isempty(r0.domain)
    # Residual is returned as a polynomial expression — not solved
    @test any(t -> !iszero(t), r0.residual)

    # --- Frac(R): numerator residual + domain metadata ---
    K = AbstractAlgebra.fraction_field(R)
    Lfrac_ok = LieAlgebra(K, 3, Dict((1, 2) => [K(0), K(0), K(a)]))
    @test check_jacobi(Lfrac_ok).ok
    @test length(parameters(Lfrac_ok)) == 2
    @test isempty(domain_denominators(Lfrac_ok))

    Lfrac_bad = LieAlgebra(K, 3, Dict(
        (1, 2) => [(a - b^2) // (a + 1), K(0), K(0)],
        (1, 3) => [K(0), K(b), K(0)],
        (2, 3) => [K(0), K(0), K(1)],
    ))
    @test !isempty(domain_denominators(Lfrac_bad))
    @test any(d -> iszero(d - (a + 1)), domain_denominators(Lfrac_bad))

    jac_f = check_jacobi(Lfrac_bad)
    @test !jac_f.ok
    @test !isempty(jac_f.residuals)
    rf = first(jac_f.residuals)
    @test all(t -> parent(t) === R, rf.residual)
    @test any(!iszero, rf.residual)
    # Domain condition from clearing / structure constants
    @test !isempty(rf.domain) || !isempty(domain_denominators(Lfrac_bad))

    # --- concrete QQ regression ---
    Hqq = LieAlgebra(Nemo.QQ, 3, Dict((1, 2) => [0, 0, 1]))
    @test check_jacobi(Hqq).ok
    @test isempty(parameters(Hqq))
end
