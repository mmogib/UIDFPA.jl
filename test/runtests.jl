using UIDFPA
using Test

@testset "UIDFPA.jl" begin
    # Write your tests here.
    n = 1_000
    F(x) = exp.(x) .- 1
    A = ones(1, n)
    b = [n]
    lb = -ones(n)
    ub = n * ones(n)
    Omega = Polyhedral(A, b, lb, ub)
    problem = NECProblem(F, Omega)
    params = UIDFPAParams(0.25, 0.5)
    x0 = ones(n)

    options = UIDFPAOptions(params, :default, :default, 1e-6, 2000, true, false)
    solution = uidfpa(problem, x0; options=options)
    @test solution.solution.Fnorm <= 1e-6
    @test maximum(abs.(solution.solution.x)) <= 1e-6
end
