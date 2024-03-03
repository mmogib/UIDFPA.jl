struct Polyhedral
    A::Matrix{Float64}
    b::Vector{Float64}
    lb::Vector{Float64}
    ub::Vector{Float64}
    project::Function
    check::Function
    function Polyhedral(A::Matrix{Float64}, b::Vector{Float64}, lb::Vector{Float64}, ub::Vector{Float64})
        n = size(A, 2)
        function project(w, z=zeros(n))
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
        new(A, b, lb, ub, project, check)
    end
end

Polyhedral(A::Matrix{Float64}, b::Vector{Float64}) = begin
    _, n = size(A)
    lb = Vector{Float64}(-Inf64, n)
    ub = Vector{Float64}(Inf64, n)

    Polyhedral(A, b, lb, ub)
end
Polyhedral(A::Matrix{Float64}, b::Vector{Float64}, lb::Vector{Float64}) = begin
    _, n = size(A)
    ub = Vector{Float64}(Inf64, n)
    Polyhedral(A, b, lb, ub)
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
