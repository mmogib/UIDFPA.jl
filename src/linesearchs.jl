struct LSParams
    k::Int
    Δ::Float64
    sk::Float64
    αstab::Float64
end
struct LSOutput
    α::Union{Nothing,Float64}
    iters::Int
    Δ::Float64
    αstab::Float64
end
function LSI(F::Function, β::Float64, σ::Float64, maxitrs::Int)
    ϵ = 1e-6
    return (dk::Vector{Float64}, xk::Vector{Float64}; params::LSParams) -> begin
        max_iters = 10_000
        for i in 0:max_iters
            αk = β^i
            arg = F(xk + αk * dk)
            lhs = -dot(arg, dk)
            rhs = σ * αk * norm(dk)^(2)
            if lhs >= rhs || αk <= ϵ
                return LSOutput(αk, i, Inf64, Inf64)
            end
        end
        LSOutput(nothing, max_iters, Inf64, Inf64)
    end
end
function LSII(F::Function, β::Float64, σ::Float64, maxitrs::Int)
    ϵ = 1e-6
    return (dk::Vector{Float64}, xk::Vector{Float64}; params::LSParams) -> begin
        max_iters = 10_000
        for i in 0:max_iters
            αk = β^i
            arg = F(xk + αk * dk)
            normFk = norm(arg)
            γk = normFk
            lhs = -dot(arg, dk)
            rhs = σ * αk * γk * norm(dk)^(2)
            if lhs >= rhs || αk <= ϵ
                return LSOutput(αk, i, Inf64, Inf64)
            end
        end
        LSOutput(nothing, max_iters, Inf64, Inf64)
    end
end

function LSIII(F::Function, β::Float64, σ::Float64, maxitrs::Int)
    ϵ = 1e-6
    return (dk::Vector{Float64}, xk::Vector{Float64}; params::LSParams) -> begin
        max_iters = 10_000
        for i in 0:max_iters
            αk = β^i
            arg = F(xk + αk * dk)
            normFk = norm(arg)
            γk = normFk / (1 + normFk)
            lhs = -dot(arg, dk)
            rhs = σ * αk * γk * norm(dk)^(2)
            if lhs >= rhs || αk <= ϵ
                return LSOutput(αk, i, Inf64, Inf64)
            end
        end
        LSOutput(nothing, max_iters, Inf64, Inf64)
    end
end


function LSIV(F::Function, β::Float64, σ::Float64, maxitrs::Int)
    b = zeros(1:(maxitrs+1))
    foreach(1:(maxitrs+1)) do i
        if 1 == i
            b[i] = 1
        elseif i == 2
            b[i] = 0.5
        else
            b[i] = 0.5 * sum(b[i-2:i-1])
        end
    end
    λ = b[2:end]
    ϵ = 1e-6
    return (dk::Vector{Float64}, xk::Vector{Float64}; params::LSParams) -> begin
        k = params.k
        max_iters = 10_000
        for i in 0:max_iters
            αk = β^i
            arg = F(xk + αk * dk)
            normFk = norm(arg)
            γk = λ[k] + (1 - λ[k]) * normFk
            lhs = -dot(arg, dk)
            rhs = σ * αk * γk * norm(dk)^(2)
            if lhs >= rhs || αk <= ϵ
                return LSOutput(αk, i, Inf64, Inf64)
            end
        end
        LSOutput(nothing, max_iters, Inf64, Inf64)
    end
end

function LSV(F::Function, β::Float64, σ::Float64, maxitrs::Int)
    ϵ = 1e-6
    return (dk::Vector{Float64}, xk::Vector{Float64}; params::LSParams) -> begin
        max_iters = 10_000
        for i in 0:max_iters
            αk = β^i
            arg = F(xk + αk * dk)
            normFk = norm(arg)
            γk = min(1, normFk)
            lhs = -dot(arg, dk)
            rhs = σ * αk * γk * norm(dk)^(2)
            if lhs >= rhs || αk <= ϵ
                return LSOutput(αk, i, Inf64, Inf64)
            end
        end
        LSOutput(nothing, max_iters, Inf64, Inf64)
    end
end

function compute_alpha_k(F::Function, β::Float64, σ::Float64, maxitrs::Int)
    ϵ = 1e-6
    η = 0.001
    ξ = 0.6
    return (d::Vector{Float64}, u::Vector{Float64}; params::LSParams) -> begin
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
                return LSOutput(tk, i, params.Δ, params.αstab)
            end
        end
        LSOutput(nothing, max_iters, params.Δ, params.αstab)
    end

end

function LSVI(F::Function, β::Float64, σ::Float64, maxitrs::Int)
    compute_alpha_k(F, β, σ, maxitrs)
end

function LSVII(F::Function, β::Float64, σ::Float64, maxitrs::Int)
    sd_fun = compute_alpha_k(F, β, σ, maxitrs)
    return (d::Vector{Float64}, u::Vector{Float64}; params::LSParams) -> begin
        k = params.k
        out1 = sd_fun(d, u; params)
        α, iters, Δ, αstab = out1.α, out1.iters, params.Δ, params.αstab
        if k <= 3
            Δ = min(Δ, params.sk)
            αstab = Δ / norm(d)
            return LSOutput(α, iters, Δ, αstab)
        end

        return LSOutput(0.25 * min(α, αstab), iters, Δ, αstab)

    end

end
"""
Line search in 
Yu, Z., Lin, J., Sun, J., Xiao, Y., Liu, L., & Li, Z. (2009). 
Spectral gradient projection method for monotone nonlinear equations with convex constraints. Applied Numerical Mathematics, 59(10), 2416–2423. https://doi.org/10.1016/j.apnum.2009.04.004
"""
function LSVIII(evalFun::Function, β::Float64, σ::Float64, maxitrs::Int)

    return (d::Vector{Float64}, u::Vector{Float64}; params::LSParams) -> begin
        It = 0
        α = 1.0
        normd = norm(d)
        xtrial = u + α * d
        Ftrial = evalFun(xtrial)
        Ftd = dot(Ftrial, d)

        while true
            if -Ftd >= σ * α * (normd^2)
                break
            else
                α /= 2.0
            end

            It += 1

            if It >= maxitrs
                return LSOutput(nothing, maxitrs, params.Δ, params.αstab)
            end

            xtrial = u .+ α .* d
            Ftrial = evalFun(xtrial)
            Ftd = dot(Ftrial, d)
        end


        return LSOutput(α, It, params.Δ, params.αstab)
    end
end

export LSI, LSII, LSIII, LSIV, LSV, LSVI, compute_alpha_k, LSVII, LSVIII