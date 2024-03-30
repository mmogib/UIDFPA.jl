
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
function successiveProject2(P, w, z, ϵ, projections)
    return P(w), 1
end
# β, ρ, μ, γ,
function uidfpa(Problem::NECProblem, x0::Vector{Float64}; options::UIDFPAOptions)
    # ζk in our paper is μk in Goncalvves (theta in MATLAB code /ϵCG) = 0.25
    @debug "Starting UIDFPA..."
    params = options.params
    ϵ = options.ϵ
    mxitrs = options.mxitrs
    LineSearch = options.LineSearch == :default ? compute_alpha_k : options.LineSearch
    SearchDirection = options.SearchDirection == :default ? compute_search_direction : options.SearchDirection
    withInertia = options.Inertia
    θ = params.θ
    ρ = params.ρ
    σ = params.σ
    ζ = params.ζ
    ζk = ζ()
    F = Problem.F
    P = options.Approximate ? Problem.Ω.approximate_project : Problem.Ω.exact_project
    Pcheck = Problem.Ω.check
    T(x, i=0) = (F(x), i + 1)
    compute_alpha = LineSearch(F, ρ, σ, mxitrs)
    compute_direction = SearchDirection()
    x0 = Pcheck(x0) ? x0 : P(x0)
    x1 = w0 = w1 = copy(x0)
    Fvalue, Fevals = T(x0)
    Fw0 = Fw1 = Fvalue
    normF = norm(Fvalue)
    k = 0
    LSitrs = PjItrs = 0
    dk = -Fvalue
    while true
        @debug "Function Norm" k normF norm_dk = norm(dk) norm_xk = norm(x1) norm_wk = norm(w1)
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

        αk, LSitrs = withInertia ? compute_alpha(dk, w1; params=Dict(:k => k)) : compute_alpha(dk, x1; params=Dict(:k => k))
        k += 1
        Fevals += LSitrs
        if isnothing(αk)
            curret_solution = SolutionInfo(x1, Fvalue, normF, k, PjItrs, LSitrs, Fevals)
            failure_solution = FailureMessage(curret_solution, "Algorithm Failed during line search.")
            return failure_solution
        end

        z1 = withInertia ? w1 + αk * dk : x1 + αk * dk
        Fzvalue, Fevals = T(z1, Fevals)
        normF = norm(Fzvalue)
        z0 = withInertia ? copy(w1) : copy(x1)
        if Pcheck(z1)
            if (normF <= ϵ)
                curret_solution = SolutionInfo(z1, Fzvalue, normF, k, PjItrs, LSitrs, Fevals)
                success_solution = SuccessSolution(curret_solution, "Solution Found")
                return success_solution
            else
                z0 = copy(z1)
            end
        end

        λk = dot(Fzvalue, -αk * dk) / normF^2
        ϵk = max((ζk * λk * normF)^2, 1e-2)
        wkfz = withInertia ? w1 - λk * Fzvalue : x1 - λk * Fzvalue
        x0, (x1, PjItrs) = x1, successiveProject2(P, wkfz, z0, ϵk, PjItrs)
        if isnothing(x1)
            curret_solution = SolutionInfo(z1, Fzvalue, normF, k, PjItrs, LSitrs, Fevals)
            failure_solution = FailureMessage(curret_solution, "Erro in linrprog in the second test")
            return failure_solution
        end


        Fw0 = copy(Fw1)
        if withInertia
            θk = norm(x1 - x0) ≈ 0 ? θ : min(1 / (k^2 * norm(x1 - x0)), θ)
            w0, w1 = w1, x1 + θk * (x1 - x0)
            Fw1, Fevals = T(w1, Fevals)
            dk = compute_direction(Fw0, Fw1, dk, w0, w1, k)
        else
            Fw1, Fevals = T(x1, Fevals)
            k == 1 && println("ji")
            dk = compute_direction(Fw0, Fw1, dk, x0, x1, k)
        end
        Fvalue = Fw1
        ζk = ζ(ζk)
        @debug "Same k" norm_zk = norm(z1)
    end
end
