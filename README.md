# UIDFPA

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mmogib.github.io/UIDFPA.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mmogib.github.io/UIDFPA.jl/dev/)
[![Build Status](https://github.com/mmogib/UIDFPA.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/mmogib/UIDFPA.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Build Status](https://app.travis-ci.com/mmogib/UIDFPA.jl.svg?branch=master)](https://app.travis-ci.com/mmogib/UIDFPA.jl)
[![Coverage](https://codecov.io/gh/mmogib/UIDFPA.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mmogib/UIDFPA.jl)
[![Coverage](https://coveralls.io/repos/github/mmogib/UIDFPA.jl/badge.svg?branch=master)](https://coveralls.io/github/mmogib/UIDFPA.jl?branch=master)


## Install
```julia
add https://github.com/mmogib/UIDFPA.jl
```

## Usage
Let $F:\mathbb{R}^n\to \mathbb{R}^n$. ``UIDFPA`` algorithm solves the problem
```math
\text{Find } x^*\in\Omega\subseteq \mathbb{R}^n \quad \text{such}\quad F(x^*)=0 
```
where 
```math
\Omega =\{x\in\mathbb{R}^n \;|\; Ax\leq b, \quad lb\leq x \leq ub\}.
```
Here is an example 
```math
F(x)=e^x-1 (\text{componentwise}), \quad \sum_{i}^n x_i \leq n, \quad  x_i\in [-1,n],\;\text{for }i\in \{1,\cdots,n\}
```
So $A=[1\;\; 1\;\; \cdots\;\; 1]$, $lb=[-1\;\; -1\;\; \cdots\;\; -1]$ and $ub=[n\;\; n\;\; \cdots\;\; n]$. To solve this example for $n=1000$.
```julia
using UIDFPA
n = 1_000;
F(x) = exp.(x) .- 1;
A = ones(1,n);
b =[1000.0];
lb = -ones(n);
ub = n*ones(n);
Omega = Polyhedral(A,b,lb,ub)
problem = NECProblem(F,Omega)
params = UIDFPAParams(0.25,0.5)
x0 = zeros(n)
solution = uidfpa(problem, x0, params, ϵ=1e-6, mxitrs=2000)
```
Here 
```julia
# The algorithm parameters
struct UIDFPAParams
    ρ::Float64 
    θ::Float64
    σ::Float64 # default 0.1
    ζ::Function #  g(initial) default 0.5
end
```
The solution is either `SuccessSolution` or `FailureMessage`.
1. `SuccessSolution`: has two fields: `solution` and `message`
   - `solution`: has the following fields
     - `x`: solution vector, i.e. $F(x)=0$.
     - Fvalue: $F(x)$.
     - Fnorm: $\|F(x)\|$
     - Iterations: number of iterations
     - ProjectionIterations: number of projections
     - LineseatchIterations: number of linesearch iterations
     - FunctionEvaluations: number of function evaluations.
   - `message`: a message indicating at what stage the solution was found.
2. `FailureMessage`: has alsw the same fields: `solution` and `message`. However the `solution` is for the last iteration before the algorithm stops. The `message` indicates the reason for failure.