struct Polyhedral
    A::Matrix{Float64}
    b::Vector{Float64}
    lb::Vector{Float64}
    ub::Vector{Float64}
    approximate_project::Function
    exact_project::Function
    check::Function
    function Polyhedral(A::Matrix{Float64}, b::Vector{Float64}, lb::Vector{Float64}, ub::Vector{Float64})
        n = size(A, 2)
        solver = optimizer_with_attributes(
            Ipopt.Optimizer,
            MOI.Silent() => true,
            "sb" => "yes",
            "max_iter" => 10_000,
        )
        function exact_project(w, z=zeros(n))
            model = Model(solver)
            set_silent(model)
            @variable(model, u[1:n])
            @constraint(model, A * u .<= b)
            @constraint(model, lb .<= u .<= ub)
            @objective(model, Min, 0.5 * dot(u, u) - dot(w, u))
            try

                optimize!(model)
                if termination_status(model) == MOI.LOCALLY_SOLVED
                    return value.(u)
                else
                    return nothing
                end
            catch err
                # println(err)
                return nothing
            end
        end
        function approx_project(w, z=zeros(n))
            model = Model(HiGHS.Optimizer)
            set_silent(model)
            @variable(model, u[1:n])
            @constraint(model, A * u .<= b)
            @constraint(model, lb .<= u .<= ub)
            @objective(model, Min, dot(w, u - z))
            try

                optimize!(model)
                if termination_status(model) == MOI.OPTIMAL
                    return value.(u)
                else
                    return nothing
                end
            catch err
                # println(err)
                return nothing
            end
        end
        function check(x, ϵ=1e-8)
            all(lb - x .<= ϵ) && all(x - ub .<= ϵ) && all(A * x - b .<= ϵ)
        end
        new(A, b, lb, ub, approx_project, exact_project, check)
    end
end
Polyhedral(A::Matrix{<:Real}, b::Vector{<:Real}, lb::Vector{<:Real}, ub::Vector{<:Real}) = Polyhedral(Float64.(A), Float64.(b), Float64.(lb), Float64.(ub))
Polyhedral(A::Matrix{<:Real}, b::Vector{<:Real}) = begin
    _, n = size(A)
    lb = Vector{Float64}(-Inf64, n)
    ub = Vector{Float64}(Inf64, n)

    Polyhedral(Float64.(A), Float64.(b), lb, ub)
end

Polyhedral(A::Matrix{<:Real}, b::Vector{<:Real}, lb::Vector{<:Real}) = begin
    _, n = size(A)
    ub = Vector{Float64}(Inf64, n)
    Polyhedral(Float64.(A), Float64.(b), Float64.(lb), ub)
end


struct SolutionInfo
    x::Union{Nothing,Vector{Float64}}
    Fvalue::Union{Nothing,Vector{Float64}}
    Fnorm::Union{Nothing,Float64}
    Iterations::Union{Nothing,Int64}
    ProjectionIterations::Union{Nothing,Int64}
    LineseatchIterations::Union{Nothing,Int64}
    FunctionEvaluations::Union{Nothing,Int64}
end
abstract type ResultMessage end
struct FailureMessage <: ResultMessage
    solution::SolutionInfo
    message::String
end

struct SuccessSolution <: ResultMessage
    solution::SolutionInfo
    message::String
end

Base.show(io::IO, ::MIME"text/plain", s::T where {T<:ResultMessage}) = begin
    msg = "message: $(s.message), iterations: $(s.solution.Iterations), projections: $(s.solution.ProjectionIterations), linesearch_iterations: $(s.solution.LineseatchIterations), evaluations: $(s.solution.FunctionEvaluations)"
    if isa(s, FailureMessage)
        print(io, "Failure: \n $msg")
    else
        print(io, "Success: \n $msg")
    end
end

# Nonlinear Equations with Constraints
struct NECProblem
    F::Function
    Ω::Polyhedral
end

struct UIDFPAParams
    ρ::Float64
    θ::Float64
    σ::Float64
    ζ::Function
end
UIDFPAParams(ρ::Float64, θ::Float64) = UIDFPAParams(ρ, θ, 0.01, (x = 0.5) -> x)
UIDFPAParams(ρ::Float64, θ::Float64, σ::Float64) = UIDFPAParams(ρ, θ, σ, (x = 0.5) -> x)

struct UIDFPAOptions
    params::UIDFPAParams
    LineSearch::Union{Symbol,Function}
    SearchDirection::Union{Symbol,Function}
    ϵ::Float64
    mxitrs::Int
    Inertia::Bool
    Approximate::Bool
end

UIDFPAOptions(p::UIDFPAParams) = UIDFPAOptions(p, :default, :default, 1e-6, 2000, false, false)
UIDFPAOptions(p::UIDFPAParams, ls::Function) = UIDFPAOptions(p, ls, :default, 1e-6, 2000, false, false)
UIDFPAOptions(p::UIDFPAParams, ls::Union{Symbol,Function}, sd::Function) = begin
    if isa(ls, Symbol)
        UIDFPAOptions(p, :default, sd, 1e-6, 2000, false, false)
    else
        UIDFPAOptions(p, ls, sd, 1e-6, 2000, false, false)
    end

end
UIDFPAOptions(p::UIDFPAParams, inertia::Bool) = UIDFPAOptions(p, :default, :default, 1e-6, 2000, inertia, false)
UIDFPAOptions(p::UIDFPAParams, ϵ::Float64) = UIDFPAOptions(p, :default, :default, ϵ, 2000, inertia, false)
UIDFPAOptions(p::UIDFPAParams, mxitrs::Int) = UIDFPAOptions(p, :default, :default, ϵ, mxitrs, inertia, false)
