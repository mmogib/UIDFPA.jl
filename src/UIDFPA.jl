module UIDFPA
using LinearAlgebra
using JuMP, HiGHS
include("types.jl")
include("algorithms.jl")

# Write your package code here.

export uidfpa,
    Polyhedral,
    SolutionInfo,
    ResultMessage,
    SuccessSolution,
    FailureMessage,
    NECProblem,
    UIDFPAParams,
    compute_alpha_k

end
