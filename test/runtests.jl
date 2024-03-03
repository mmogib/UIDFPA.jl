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
    x0 = zeros(n)
    solution = uidfpa(problem, x0, params, ϵ=1e-6, mxitrs=2000)
    @test solution.solution.Fnorm <= 1e-6
    @test solution.solution.x ≈ zeros(n)
end
