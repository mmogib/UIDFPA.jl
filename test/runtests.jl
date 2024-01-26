using UIDFPA
using Test

@testset "UIDFPA.jl" begin
    # Write your tests here.

    @test 3 == firstfun(1, 2)
    @test 1 == firstfun(1, 2)
end
