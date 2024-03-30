function LSI(F::Function, β::Float64, σ::Float64, maxitrs::Int)
    ϵ = 1e-6
    return (dk::Vector{Float64}, xk::Vector{Float64}; params::Dict{Symbol,<:Number}) -> begin
        max_iters = 10_000
        for i in 0:max_iters
            αk = β^i
            arg = F(xk + αk * dk)
            lhs = -dot(arg, dk)
            rhs = σ * αk * norm(dk)^(2)
            if lhs >= rhs || αk <= ϵ
                return αk, i
            end
        end
        nothing, max_iters
    end
end
function LSII(F::Function, β::Float64, σ::Float64, maxitrs::Int)
    ϵ = 1e-6
    return (dk::Vector{Float64}, xk::Vector{Float64}; params::Dict{Symbol,<:Number}) -> begin
        max_iters = 10_000
        for i in 0:max_iters
            αk = β^i
            arg = F(xk + αk * dk)
            normFk = norm(arg)
            γk = normFk
            lhs = -dot(arg, dk)
            rhs = σ * αk * γk * norm(dk)^(2)
            if lhs >= rhs || αk <= ϵ
                return αk, i
            end
        end
        nothing, max_iters
    end
end

function LSIII(F::Function, β::Float64, σ::Float64, maxitrs::Int)
    ϵ = 1e-6
    return (dk::Vector{Float64}, xk::Vector{Float64}; params::Dict{Symbol,<:Number}) -> begin
        max_iters = 10_000
        for i in 0:max_iters
            αk = β^i
            arg = F(xk + αk * dk)
            normFk = norm(arg)
            γk = normFk / (1 + normFk)
            lhs = -dot(arg, dk)
            rhs = σ * αk * γk * norm(dk)^(2)
            if lhs >= rhs || αk <= ϵ
                return αk, i
            end
        end
        nothing, max_iters
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
    return (dk::Vector{Float64}, xk::Vector{Float64}; params::Dict{Symbol,<:Number}) -> begin
        k = get(params, :k, 1)
        max_iters = 10_000
        for i in 0:max_iters
            αk = β^i
            arg = F(xk + αk * dk)
            normFk = norm(arg)
            γk = λ[k] + (1 - λ[k]) * normFk
            lhs = -dot(arg, dk)
            rhs = σ * αk * γk * norm(dk)^(2)
            if lhs >= rhs || αk <= ϵ
                return αk, i
            end
        end
        nothing, max_iters
    end
end

function LSV(F::Function, β::Float64, σ::Float64, maxitrs::Int)
    ϵ = 1e-6
    return (dk::Vector{Float64}, xk::Vector{Float64}; params::Dict{Symbol,<:Number}) -> begin
        max_iters = 10_000
        for i in 0:max_iters
            αk = β^i
            arg = F(xk + αk * dk)
            normFk = norm(arg)
            γk = min(1, normFk)
            lhs = -dot(arg, dk)
            rhs = σ * αk * γk * norm(dk)^(2)
            if lhs >= rhs || αk <= ϵ
                return αk, i
            end
        end
        nothing, max_iters
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

function LSVI(F::Function, β::Float64, σ::Float64, maxitrs::Int)
    compute_alpha_k(F, β, σ, maxitrs)
end



export LSI, LSII, LSIII, LSIV, LSV, LSVI, compute_alpha_k