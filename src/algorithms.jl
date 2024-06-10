
function successiveProject(P, w, z, ϵ, projections, B=I, is_approximate=true)
    if !is_approximate
        return P(w), projections + 1
    end
    for i in 1:300
        try
            wp = B * z - w
            zbar = P(wp)
            if isnothing(zbar)
                return nothing, projections + i
            end
            FvalAprox = dot(wp, zbar - z)
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
    return P(w), projections + 1
end
# β, ρ, μ, γ,
function uidfpa(Problem::NECProblem, x0::Vector{Float64}; options::UIDFPAOptions)
    @debug "Starting UIDFPA..."
    # ζk in our paper is μk in Goncalvves (theta in MATLAB code /ϵCG) = 0.25

    # Unpacking options and parameters
    params = options.params
    ϵ = options.ϵ
    mxitrs = options.mxitrs
    timelimit = options.timelimit
    gather_norms = options.GatherNorms
    withInertia = options.Inertia
    θ = params.θ
    ρ = params.ρ
    σ = params.σ
    ζ = params.ζ
    ζk = ζ()

    # Determine the line search and search direction methods
    LineSearch = options.LineSearch == :default ? compute_alpha_k : options.LineSearch
    SearchDirection = options.SearchDirection == :default ? compute_search_direction : options.SearchDirection

    # Function and projection operators
    F = Problem.F
    is_approximate = options.Approximate
    P = is_approximate ? Problem.Ω.approximate_project : Problem.Ω.exact_project
    Pcheck = Problem.Ω.check
    T(x, i=0) = (F(x), i + 1)

    # Initial computations
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
    lsparams = LSParams(1, Inf64, Inf64, Inf64)

    # Setting timers
    time_start = Dates.now()
    FNORMS = if gather_norms
        -1.0 * ones(mxitrs)
    else
        nothing
    end

    # Main loop
    while true
        @debug "Iteration: $k, Function Norm: $normF"
        current_time = Dates.now()
        if round(current_time - time_start, Second) > Second(timelimit)
            curret_solution = SolutionInfo(x1, Fvalue, normF, k, PjItrs, LSitrs, Fevals, FNORMS)
            failure_solution = FailureMessage(curret_solution, "Maximum time limit: $timelimit seconds reached")
            return failure_solution
        end
        if normF <= ϵ
            curret_solution = SolutionInfo(x1, Fvalue, normF, k, PjItrs, LSitrs, Fevals, FNORMS)
            success_solution = SuccessSolution(curret_solution, "Solution Found")
            return success_solution
        end

        if PjItrs >= mxitrs
            curret_solution = SolutionInfo(x1, Fvalue, normF, k, PjItrs, LSitrs, Fevals, FNORMS)
            failure_solution = FailureMessage(curret_solution, "Maximum number of projections allowed reached")
            return failure_solution
        end


        if k >= mxitrs
            curret_solution = SolutionInfo(x1, Fvalue, normF, k, PjItrs, LSitrs, Fevals, FNORMS)
            failure_solution = FailureMessage(curret_solution, "Maximum number of iterations allowed reached")
            return failure_solution
        end

        k += 1
        lsout = withInertia ? compute_alpha(dk, w1; params=lsparams) : compute_alpha(dk, x1; params=lsparams)
        αk, LSitrs, Δ, αstab = lsout.α, lsout.iters, lsout.Δ, lsout.αstab
        @debug "Lineserach out" αk
        Fevals += LSitrs

        if isnothing(αk)
            curret_solution = SolutionInfo(x1, Fvalue, normF, k, PjItrs, LSitrs, Fevals, FNORMS)
            failure_solution = FailureMessage(curret_solution, "Algorithm Failed during line search.")
            return failure_solution
        end

        z1 = withInertia ? w1 + αk * dk : x1 + αk * dk
        Fzvalue, Fevals = T(z1, Fevals)
        normF = norm(Fzvalue)
        if gather_norms
            FNORMS[k] = normF
        end
        # z0 = withInertia ? copy(w1) : copy(x1)

        if Pcheck(z1) && (normF <= ϵ)
            curret_solution = SolutionInfo(z1, Fzvalue, normF, k, PjItrs, LSitrs, Fevals, FNORMS)
            success_solution = SuccessSolution(curret_solution, "Solution Found")
            return success_solution
        end

        z0 = copy(z1)
        λk = dot(Fzvalue, -αk * dk) / normF^2
        ϵk = max((ζk * λk * normF)^2, 1e-2)
        wkfz = withInertia ? w1 - λk * Fzvalue : x1 - λk * Fzvalue
        x0, (x1, PjItrs) = x1, successiveProject(P, wkfz, z0, ϵk, PjItrs, I, is_approximate)

        if isnothing(x1)
            curret_solution = SolutionInfo(z1, Fzvalue, normF, k, PjItrs, LSitrs, Fevals, FNORMS)
            failure_solution = FailureMessage(curret_solution, "Erro in Projection")
            return failure_solution
        end

        Fw0 = copy(Fw1)

        if withInertia
            θk = norm(x1 - x0) ≈ 0 ? θ : min(1 / (k^2 * norm(x1 - x0)), θ)
            w0, w1 = w1, x1 + θk * (x1 - x0)
            lsparams = LSParams(k, Δ, norm(w1 - w0), αstab)
            Fw1, Fevals = T(w1, Fevals)
            dk = compute_direction(Fw0, Fw1, dk, w0, w1, k)
            @debug "Search direction with inertia: θk = $θk"
        else
            Fw1, Fevals = T(x1, Fevals)
            lsparams = LSParams(k, Δ, norm(x1 - x0), αstab)
            dk = compute_direction(Fw0, Fw1, dk, x0, x1, k)
        end

        Fvalue = Fw1
        ζk = ζ(ζk)
    end
end
