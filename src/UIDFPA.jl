module UIDFPA
using LinearAlgebra, Dates
using JuMP, HiGHS, Ipopt
# import Clarabel
include("types.jl")
include("searchdirections.jl")
include("linesearchs.jl")
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
    compute_alpha_k,
    UIDFPAOptions
end
