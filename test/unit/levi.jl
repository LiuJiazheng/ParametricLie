@testset "radical / quotient / levi" begin
    F = Nemo.QQ

    # Shared examples
    sl2 = LieAlgebra(F, 3, Dict(
        (1, 2) => [0, 2, 0],
        (1, 3) => [0, 0, -2],
        (2, 3) => [1, 0, 0],
    ))
    H = LieAlgebra(F, 3, Dict((1, 2) => [0, 0, 1]))
    A = LieAlgebra(F, 2)

    # --- radical = cartan_orthogonal; ideal ---
    @test dim(radical(sl2)) == dim(cartan_orthogonal(sl2))
    @test dim(radical(sl2)) == 0
    @test is_ideal(sl2, radical(sl2))

    @test dim(radical(H)) == 3
    @test is_ideal(H, radical(H))
    @test is_solvable(H)

    @test dim(radical(A)) == 2

    # --- complement ⊕ S = L ---
    Rsl = radical(sl2)
    Xsl = complement(sl2, Rsl)
    @test dim(Xsl) == 3
    @test is_subalgebra(sl2, Xsl)  # whole algebra

    RH = radical(H)
    XH = complement(H, RH)
    @test dim(XH) == 0

    # --- quotient: sl₂ ---
    Qsl = quotient_algebra(sl2)
    @test dim(Qsl) == 3
    @test dim(Qsl.radical) == 0
    @test check_jacobi(Qsl.algebra).ok
    @test is_semisimple(Qsl.algebra)
    # Structure matches sl₂ (same basis)
    fQ = structure_constants(Qsl.algebra)
    fL = structure_constants(sl2)
    @test fQ == fL

    # --- quotient: solvable → 0 ---
    QH = quotient_algebra(H)
    @test dim(QH) == 0
    QA = quotient_algebra(A)
    @test dim(QA) == 0

    # --- direct sum sl₂ ⊕ ℚ (abelian radical) ---
    # e1,e2,e3 = sl₂; e4 central
    sl2_plus = LieAlgebra(F, 4, Dict(
        (1, 2) => [0, 2, 0, 0],
        (1, 3) => [0, 0, -2, 0],
        (2, 3) => [1, 0, 0, 0],
    ))
    @test check_jacobi(sl2_plus).ok
    R4 = radical(sl2_plus)
    @test dim(R4) == 1
    @test is_ideal(sl2_plus, R4)
    # radical = span{e4}
    z = only(basis_elems(R4))
    @test iszero(z.coords[1]) && iszero(z.coords[2]) && iszero(z.coords[3])
    @test !iszero(z.coords[4])

    Q4 = quotient_algebra(sl2_plus)
    @test dim(Q4) == 3
    @test check_jacobi(Q4.algebra).ok
    @test is_semisimple(Q4.algebra)
    @test killing_rank(Q4.algebra) == 3

    lev4 = levi_decomposition(sl2_plus)
    @test dim(lev4.radical) == 1
    @test dim(lev4.levi) == 3
    @test is_subalgebra(sl2_plus, lev4.levi)
    @test is_semisimple(lev4.quotient)
    # Levi ⊕ radical = ambient as vector spaces
    @test dim(lev4.levi) + dim(lev4.radical) == dim(sl2_plus)
    # Bracket relations among lifts match quotient structure constants
    f4 = structure_constants(lev4.quotient)
    m = length(lev4.lifts)
    for a in 1:m, b in 1:m
        br = lie_bracket(sl2_plus, lev4.lifts[a], lev4.lifts[b])
        expected = zero(sl2_plus)
        for c in 1:m
            fc = f4[a, b, c]
            iszero(fc) && continue
            for i in 1:4
                expected.coords[i] += fc * lev4.lifts[c].coords[i]
            end
        end
        @test br == expected
    end

    # --- mixed basis: complement not yet a subalgebra ---
    # Columns of P = new basis in old coords:
    # e1'=e1+e4, e2'=e2+e4, e3'=e3+e4, e4'=e4
    Pmix = AbstractAlgebra.matrix(F, 4, 4, [
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        1, 1, 1, 1,
    ])
    Lmix = change_of_basis(sl2_plus, Pmix)
    @test check_jacobi(Lmix).ok
    @test dim(radical(Lmix)) == 1
    @test !is_semisimple(Lmix)

    # Raw complement need not be a subalgebra
    Xmix = complement(Lmix, radical(Lmix))
    # After Levi correction it must be
    levm = levi_decomposition(Lmix)
    @test dim(levm.levi) == 3
    @test dim(levm.radical) == 1
    @test is_subalgebra(Lmix, levm.levi)
    @test is_semisimple(levm.quotient)
    @test check_jacobi(levm.quotient).ok
    # Direct sum of subspaces
    Blev = basis_matrix(levm.levi)
    Brad = basis_matrix(levm.radical)
    Pfull = hcat(Blev, Brad)
    @test AbstractAlgebra.is_unit(AbstractAlgebra.det(Pfull))

    # Verify structure-constant relations for corrected lifts
    fm = structure_constants(levm.quotient)
    mm = length(levm.lifts)
    for a in 1:mm, b in 1:mm
        br = lie_bracket(Lmix, levm.lifts[a], levm.lifts[b])
        expected_coords = fill(zero(F), 4)
        for c in 1:mm
            fc = fm[a, b, c]
            iszero(fc) && continue
            for i in 1:4
                expected_coords[i] += fc * levm.lifts[c].coords[i]
            end
        end
        @test br.coords == expected_coords
    end

    # --- non-abelian radical: sl₂ ⊕ Heisenberg ---
    # dim 6: e1..e3 = sl₂, e4,e5,e6 = Heisenberg ([e4,e5]=e6)
    sl2H = LieAlgebra(F, 6, Dict(
        (1, 2) => [0, 2, 0, 0, 0, 0],
        (1, 3) => [0, 0, -2, 0, 0, 0],
        (2, 3) => [1, 0, 0, 0, 0, 0],
        (4, 5) => [0, 0, 0, 0, 0, 1],
    ))
    @test check_jacobi(sl2H).ok
    R6 = radical(sl2H)
    @test dim(R6) == 3
    @test is_ideal(sl2H, R6)
    # Derived series of radical: 3 → 1 → 0 (Heisenberg)
    rser = radical_derived_series(sl2H, R6)
    @test map(dim, rser) == [3, 1, 0]

    Q6 = quotient_algebra(sl2H)
    @test dim(Q6) == 3
    @test is_semisimple(Q6.algebra)

    lev6 = levi_decomposition(sl2H)
    @test dim(lev6.levi) == 3
    @test dim(lev6.radical) == 3
    @test is_subalgebra(sl2H, lev6.levi)
    @test is_semisimple(lev6.quotient)

    # Mix radical into complement, require multi-layer lift
    P6 = AbstractAlgebra.matrix(F, 6, 6, [
        1, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0,
        1, 0, 0, 1, 0, 0,
        0, 1, 0, 0, 1, 0,
        0, 0, 1, 0, 0, 1,
    ])
    L6m = change_of_basis(sl2H, P6)
    lev6m = levi_decomposition(L6m)
    @test dim(lev6m.levi) == 3
    @test is_subalgebra(L6m, lev6m.levi)
    @test AbstractAlgebra.is_unit(AbstractAlgebra.det(hcat(basis_matrix(lev6m.levi), basis_matrix(lev6m.radical))))
    f6 = structure_constants(lev6m.quotient)
    for a in 1:3, b in 1:3
        br = lie_bracket(L6m, lev6m.lifts[a], lev6m.lifts[b])
        expected_coords = fill(zero(F), 6)
        for c in 1:3
            fc = f6[a, b, c]
            iszero(fc) && continue
            for i in 1:6
                expected_coords[i] += fc * lev6m.lifts[c].coords[i]
            end
        end
        @test br.coords == expected_coords
    end

    # --- pure semisimple / pure solvable Levi edge cases ---
    lev_sl = levi_decomposition(sl2)
    @test dim(lev_sl.levi) == 3
    @test dim(lev_sl.radical) == 0
    @test is_subalgebra(sl2, lev_sl.levi)

    lev_H = levi_decomposition(H)
    @test dim(lev_H.levi) == 0
    @test dim(lev_H.radical) == 3

    @test dim(levi_subalgebra(sl2_plus)) == 3

    # --- ideal_decomposition: see unit/ideal_decomp.jl ---
    @test_throws ArgumentError ideal_decomposition(H)  # not semisimple
end
