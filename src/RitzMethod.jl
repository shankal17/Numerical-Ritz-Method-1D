module RitzMethod

using ArgCheck
using ForwardDiff: derivative
using FastGaussQuadrature: gausslegendre
using LinearAlgebra: dot
using ADTypes: AutoForwardDiff
using Optim: optimize, BFGS, AbstractOptimizer, Options

export Interval, minimizefunctional

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
An `Optim.OptimizationResults` object; the minimizing coefficients are available via `Optim.minimizer(result)`.

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
    method::AbstractOptimizer = BFGS(),
    options::Options = Options(),
)
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
    res = optimize(objective, initial_coeffs, method, options; autodiff = AutoForwardDiff())

    return res
end


end # module ritz_method
