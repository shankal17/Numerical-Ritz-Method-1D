module RitzMethod

import Optim

using ArgCheck
using ForwardDiff: derivative
using FastGaussQuadrature: gausslegendre
using LinearAlgebra: dot
using ADTypes: AutoForwardDiff

export Interval, minimizefunctional, RitzResult, summarize

"""
    Interval(low, high)

Closed interval `[low, high]` used to define a domain.

# Fields
- `low::Float64`: Lower bound
- `high::Float64`: Upper bound

# Throws
- `ArgumentError` (via `@argcheck`) if `low >= high`.
"""
struct Interval
    low::Float64
    high::Float64
    function Interval(low, high)
        @argcheck low < high
        new(low, high)
    end

end

"""
    minimizefunctional(F, domain, phi_s; phi_0, num_quad_nodes, initial_coeffs, method, options)

Approximate the minimizer of the functional

```math
J[y] = \\int_{a}^{b} F(x, y(x), y'(x))\\, dx
```

over `domain` `= [a, b]` using the Ritz method, where the trial solution is

```math
y(x) = \\phi_0(x) + \\sum_{i=1}^{n} c_i \\phi_i(x)
```

The integral is evaluated via Gauss-Legendre quadrature, and the coefficients
`c_i` are found by minimizing the resulting objective with BFGS.

# Arguments
- `F::Function`: the integrand, called as `F(x, y, y')`.
- `domain::Interval`: the interval `[a, b]` of integration.
- `phi_s::Vector`: basis functions `\\phi_i` (callable) used to build the trial solution.

# Keyword Arguments
- `phi_0::Function = x -> 0.0`: fixed term `\\phi_0(x)` added to the trial solution, typically used to satisfy boundary conditions.
- `num_quad_nodes::Int = 15`: number of Gauss-Legendre quadrature nodes.
- `initial_coeffs::Vector{Float64} = zeros(length(phi_s))`: initial guess for the coefficients `c_i`.
- `method::Optim.AbstractOptimizer = BFGS()`: optimization algorithm used to minimize the objective. See the [Optim.jl documentation](https://julianlsolvers.github.io/Optim.jl/stable/) for available algorithms.
- `options::Optim.Options = Optim.Options()`: convergence and iteration settings for `method`.

# Returns
A [`RitzResult`](@ref) containing:
- `coeffs`: the optimized coefficients `c_i`.
- `minimized_value`: the minimized value of the functional `J[y]`.
- `domain`: the interval `[a, b]` the functional was minimized over.
- `converged`: whether the optimizer reported convergence.
- `y`: the trial solution as a callable function, `y(x) = phi_0(x) + sum(c_i * phi_i(x))`.

For example, `result.coeffs` gives the coefficient vector and `result.y(0.5)` evaluates
the approximate solution at `x = 0.5`.

# Throws
- `ArgumentError` (via `@argcheck`) if `length(phi_s) != length(initial_coeffs)` or `num_quad_nodes <= 0`.
"""
function minimizefunctional(
    F::Function,
    domain::Interval,
    phi_s::Vector;
    phi_0::Function = x-> 0.0,
    num_quad_nodes::Int = 15,
    initial_coeffs::Vector{Float64} = zeros(length(phi_s)),
    method::Optim.AbstractOptimizer = Optim.BFGS(),
    options::Optim.Options = Optim.Options(),
)::RitzResult
    @argcheck length(phi_s) == length(initial_coeffs)
    @argcheck num_quad_nodes > 0

    nodes, weights = gausslegendre(num_quad_nodes)

    # Change of interval
    coeff = (domain.high - domain.low) / 2 # Also coefficient that is applied during quadrature
    bias = (domain.low + domain.high) / 2
    new_nodes = map(x_i -> coeff * x_i + bias, nodes)

    # Precompute basis functions and their derivatives at nodes
    phi_0_vals = map(x -> phi_0(x), new_nodes)
    phi_i_vals = [phi(x) for x in new_nodes, phi in phi_s]

    phi_0_prime_vals = map(x -> derivative(phi_0, x), new_nodes)
    phi_i_prime_vals = [derivative(phi, x) for x in new_nodes, phi in phi_s]

    function objective(basis_coeffs)
        y_vals = phi_0_vals .+ phi_i_vals * basis_coeffs
        y_prime_vals = phi_0_prime_vals .+ phi_i_prime_vals * basis_coeffs
        F_vals = F.(new_nodes, y_vals, y_prime_vals)

        return coeff * dot(weights, F_vals)
    end

    # Optimize
    res = Optim.optimize(objective, initial_coeffs, method, options; autodiff = AutoForwardDiff())

    return postprocessoptimization(res, domain, phi_s, phi_0)
end

"""
    RitzResult

Result of a Ritz method optimization, as returned by [`minimizefunctional`](@ref).

# Fields
- `coeffs::Vector{Float64}`: the optimized coefficients `c_i`.
- `minimized_value::Float64`: the minimized value of the functional `J[y]`.
- `domain::Interval`: the interval `[a, b]` over which the functional was minimized.
- `converged::Bool`: whether the underlying optimizer reported convergence.
- `y::Function`: the trial solution `y(x) = phi_0(x) + sum(c_i * phi_i(x))`, callable as `y(x)`.
"""
struct RitzResult
    coeffs::Vector{Float64}
    minimized_value::Float64
    domain::Interval
    converged::Bool
    y::Function
end

"""
    summarize(result::RitzResult)::String

Build a human-readable summary of a [`RitzResult`](@ref).

# Arguments
- `result::RitzResult`: the result to summarize.

# Returns
A multi-line `String` describing the domain, optimized coefficients,
minimized functional value, and convergence status.
"""
function summarize(result::RitzResult)::String
    return """
    RitzResult:
      Domain:          [$(result.domain.low), $(result.domain.high)]
      Coefficients:    $(result.coeffs)
      Minimized value: $(result.minimized_value)
      Converged:       $(result.converged)"""
end

Base.show(io::IO, result::RitzResult) = print(io, summarize(result))

"""
    postprocessoptimization(res, domain, phi_s, phi_0)

Build a [`RitzResult`](@ref) from the raw output of `Optim.optimize`.

# Arguments
- `res::Optim.OptimizationResults`: the result returned by `Optim.optimize`.
- `domain::Interval`: the interval `[a, b]` the functional was minimized over.
- `phi_s::Vector`: basis functions `\\phi_i` used to build the trial solution.
- `phi_0::Function`: fixed term `\\phi_0(x)` added to the trial solution.

# Returns
A [`RitzResult`](@ref) containing the optimized coefficients, minimized value,
domain, convergence flag, and the resulting trial solution `y(x) = phi_0(x) + sum(c_i * phi_i(x))`.
"""
function postprocessoptimization(
    res::Optim.OptimizationResults,
    domain::Interval,
    phi_s::Vector,
    phi_0::Function,
)::RitzResult
    coeffs = copy(Optim.minimizer(res))
    minimized_value = Optim.minimum(res)
    conv = Optim.converged(res)

    y(x) = phi_0(x) + sum(coeffs[i] * phi_s[i](x) for i in eachindex(phi_s); init = 0.0)

    return RitzResult(coeffs, minimized_value, domain, conv, y)
end

end # module ritz_method
