
function successiveProject(P, w, z, ϵ, projections)
    for i in 1:300
        try
            zbar = P(z - w, z)
            if isnothing(zbar)
                return nothing, projections + i
            end
            FvalAprox = dot(z - w, zbar - z)
            if (FvalAprox >= -ϵ)
                return z, projections + i
            end
            pcg = -FvalAprox
            qcg = norm(zbar - z)^2
            β = min(1, pcg / qcg)
            z = (1 - β) * z + β * zbar
        catch
            return nothing, projections + i
        end
    end
    return z, 300
end
# β, ρ, μ, γ,
function uidfpa(Problem::NECProblem, x0::Vector{Float64}, params::UIDFPAParams; LS::Function=compute_alpha_k, ϵ::Float64=1e-6, mxitrs::Int64=2000)
    # ζk in our paper is μk in Goncalvves (theta in MATLAB code /ϵCG) = 0.25

    θ = params.θ
    ρ = params.ρ
    σ = params.σ
    ζ = params.ζ
    ζk = ζ()
    F = Problem.F
    P = Problem.Ω.project
    Pcheck = Problem.Ω.check
    T(x, i=0) = (F(x), i + 1)
    compute_alpha = LS(F, ρ, σ, mxitrs)
    x1 = w0 = w1 = copy(x0)
    Fvalue, Fevals = T(x0)
    Fw0 = Fw1 = Fvalue
    normF = norm(Fvalue)
    k = 1
    LSitrs = PjItrs = 0
    dk = -Fvalue
    while true
        if normF <= ϵ
            curret_solution = SolutionInfo(x1, Fvalue, normF, k, PjItrs, LSitrs, Fevals)
            success_solution = SuccessSolution(curret_solution, "Solution Found")
            return success_solution
        end

        if PjItrs >= 10000
            curret_solution = SolutionInfo(x1, Fvalue, normF, k, PjItrs, LSitrs, Fevals)
            failure_solution = FailureMessage(curret_solution, "Maximum number of projections allowed reached")
            return failure_solution
        end


        if k >= mxitrs
            curret_solution = SolutionInfo(x1, Fvalue, normF, k, PjItrs, LSitrs, Fevals)
            failure_solution = FailureMessage(curret_solution, "Maximum number of iterations allowed reached")
            return failure_solution
        end

        αk, LSitrs = compute_alpha(dk, w1; params=Dict(:k => k))
        k += 1
        Fevals += LSitrs
        if isnothing(αk)
            curret_solution = SolutionInfo(x1, Fvalue, normF, k, PjItrs, LSitrs, Fevals)
            failure_solution = FailureMessage(curret_solution, "Algorithm Failed during line search.")
            return failure_solution
        end

        z1 = w1 + αk * dk
        Fzvalue, Fevals = T(z1, Fevals)
        normFz = norm(Fzvalue)
        z0 = copy(w1)
        if Pcheck(z1)
            if (normFz <= ϵ)
                curret_solution = SolutionInfo(z1, Fzvalue, normFz, k, PjItrs, LSitrs, Fevals)
                success_solution = SuccessSolution(curret_solution, "Solution Found")
                return success_solution
            else
                z0 = copy(z1)
            end
        end

        λk = dot(Fzvalue, -αk * dk) / normFz^2
        ϵk = max((ζk * λk * normFz)^2, 1e-2)
        x0, (x1, PjItrs) = x1, successiveProject(P, w1 - λk * Fzvalue, z0, ϵk, PjItrs)
        if isnothing(x1)
            curret_solution = SolutionInfo(z1, Fzvalue, normFz, k, PjItrs, LSitrs, Fevals)
            failure_solution = FailureMessage(curret_solution, "Erro in linrprog in the second test")
            return failure_solution
        end

        θk = norm(x1 - x0) ≈ 0 ? θ : min(1 / (k^2 * norm(x1 - x0)), θ)
        w0, w1 = w1, x1 + θk * (x1 - x0)
        Fw0 = copy(Fw1)
        Fw1, Fevals = T(w1, Fevals)

        dk = compute_search_direction(Fw0, Fw1, dk, w0, w1)
        Fvalue = Fw1
        ζk = ζ(ζk)
    end
end

function compute_alpha_k(F::Function, β::Float64, σ::Float64, maxitrs::Int)
    ϵ = 1e-6
    η = 0.001
    ξ = 0.6
    return (d::Vector{Float64}, u::Vector{Float64}; params::Dict{Symbol,<:Number}) -> begin
        max_iters = 10_000
        for i in 0:max_iters
            tk = β^i
            arg = F(u + tk .* d)
            lhs = -dot(arg, d)
            χ = norm(arg)
            Pηξχ = min(ξ, max(χ, η))
            rhs = σ * tk * Pηξχ * norm(d)^2
            # rhs = σ * tk * norm(d)^2
            if lhs >= rhs || tk <= ϵ
                return tk, i
            end
        end
        nothing, max_iters
    end

end


function compute_search_direction(Fu0, Fu1, d0, u0, u1, r=0.1, ψ=0.2, αmin=1.0e-10, αmax=Inf64)
    y0 = Fu1 - Fu0
    s0 = u1 - u0 + r * y0
    # v1 = max(ψ * norm(d0) * norm(y0), dot(d0, y0), norm(Fu0)^2)
    v1 = max(ψ * norm(d0) * norm(y0), norm(Fu0)^2)
    β1 = dot(Fu1, y0) / v1
    α12 = dot(Fu1, d0) / v1
    α11 = min(αmax, max(αmin, dot(s0, y0) / dot(y0, y0)))

    return -α11 * Fu1 + β1 * d0 - α12 * y0

end