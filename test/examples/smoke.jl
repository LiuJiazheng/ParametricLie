@testset "examples smoke" begin
    # Fixture loader will live here once algebras are constructible.
    # For now, ensure the package loads and fixture paths exist.
    root = dirname(dirname(@__DIR__))
    fixtures = joinpath(root, "fixtures")
    @test isdir(fixtures)
    for name in ("abelian", "heisenberg", "sl2", "affine", "bianchi", "nilpotent")
        @test isdir(joinpath(fixtures, name))
    end
end
