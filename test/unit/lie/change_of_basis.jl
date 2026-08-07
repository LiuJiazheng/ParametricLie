@testset "change_of_basis" begin
    F = Nemo.QQ
    H = LieAlgebra(F, 3, Dict((1, 2) => [0, 0, 1]))

    # Identity: unchanged structure constants
    Id = AbstractAlgebra.identity_matrix(F, 3)
    H2 = change_of_basis(H, Id)
    @test structure_constants(H2) == structure_constants(H)
    @test check_jacobi(H2).ok

    # Generic GL(3,Q) change
    P = AbstractAlgebra.matrix(F, 3, 3, [1, 1, 0,
                                         0, 1, 1,
                                         0, 0, 1])
    Hp = change_of_basis(H, P)
    @test check_jacobi(Hp).ok
    @test isempty(check_antisymmetry(Hp))

    # Bracket compatibility: P [x',y'] = [P x', P y'] in old algebra
    xp = LieAlgebraElem(Hp, [1, 0, 0])
    yp = LieAlgebraElem(Hp, [0, 1, 0])
    zp = lie_bracket(Hp, xp, yp)

    # old coords = P * new coords
    function apply_P(P, v)
        n = length(v)
        col = AbstractAlgebra.matrix(F, n, 1, v)
        w = P * col
        return [w[i, 1] for i in 1:n]
    end
    x_old = apply_P(P, xp.coords)
    y_old = apply_P(P, yp.coords)
    z_old = lie_bracket(H, x_old, y_old)
    @test z_old == apply_P(P, zp.coords)

    # Element transform
    e1 = basis_elem(H, 1)
    e1p = change_of_basis(e1, P)
    @test apply_P(P, e1p.coords) == e1.coords
end
