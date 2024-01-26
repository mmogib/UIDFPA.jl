using UIDFPA
using Test

@testset "UIDFPA.jl" begin
    # Write your tests here.

    @test 3 == firstfun(1, 2)
    @test -3 == firstfun(-5, 2)
end
