@testset "lie_bracket / ad / Jacobi" begin
    F = Nemo.QQ

    # Abelian: all brackets zero, Jacobi OK
    A = LieAlgebra(F, 3)
    x = LieAlgebraElem(A, [1, 2, 3])
    y = LieAlgebraElem(A, [4, 5, 6])
    @test iszero(lie_bracket(A, x, y).coords)
    @test check_jacobi(A).ok

    # Heisenberg: [e1,e2]=e3
    H = LieAlgebra(F, 3, Dict((1, 2) => [0, 0, 1]))
    e1, e2, e3 = basis_elem(H, 1), basis_elem(H, 2), basis_elem(H, 3)

    z = lie_bracket(H, e1, e2)
    @test z == e3
    @test lie_bracket(H, e2, e1) == LieAlgebraElem(H, [0, 0, -1])
    @test iszero(lie_bracket(H, e1, e1).coords)
    @test iszero(lie_bracket(H, e1, e3).coords)

    # vector API
    @test lie_bracket(H, [1, 0, 0], [0, 1, 0]) == [F(0), F(0), F(1)]

    # antisymmetry on general vectors
    u = LieAlgebraElem(H, [1, 2, 0])
    v = LieAlgebraElem(H, [3, 4, 0])
    @test lie_bracket(H, u, v).coords == map(-, lie_bracket(H, v, u).coords)

    # ad(x) * y == [x,y]
    Ax = ad(H, e1)
    @test Ax isa AbstractAlgebra.MatElem
    ycol = AbstractAlgebra.matrix(F, 3, 1, e2.coords)
    Ay = Ax * ycol
    @test [Ay[k, 1] for k in 1:3] == lie_bracket(H, e1, e2).coords

    jac = check_jacobi(H)
    @test jac.ok
    @test isempty(jac.antisym_failing)
    @test isempty(jac.jacobi_failing)

    # Structure that violates Jacobi (but is antisymmetric via DOF ctor)
    B = LieAlgebra(F, 3, Dict(
        (1, 2) => [1, 0, 0],
        (1, 3) => [0, 1, 0],
    ))
    jacB = check_jacobi(B)
    @test !jacB.ok
    @test isempty(jacB.antisym_failing)
    @test !isempty(jacB.jacobi_failing)

    # Dense tensor that breaks antisymmetry
    c = fill(F(0), 2, 2, 2)
    c[1, 2, 1] = F(1)
    c[2, 1, 1] = F(0)   # should be -1
    Bad = LieAlgebra(F, c)
    jacBad = check_jacobi(Bad)
    @test !jacBad.ok
    @test !isempty(jacBad.antisym_failing)
    @test (1, 2, 1) in jacBad.antisym_failing || (2, 1, 1) in jacBad.antisym_failing
end
