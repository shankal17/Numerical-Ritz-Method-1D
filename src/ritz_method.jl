module ritz_method

using ArgCheck
using ForwardDiff: derivative
using FastGaussQuadrature: gausslegendre
using LinearAlgebra: dot
using ADTypes: AutoForwardDiff
using Optim: optimize, BFGS

export Interval, ritz

# TODO: Add docstrings
struct Interval
    low::Float64
    high::Float64
    function Interval(low, high)
        @argcheck low < high
        new(low, high)
    end

end

# TODO: Add docstrings
function ritz(
    F::Function,
    domain::Interval,
    phi_s::Vector{<:Function};
    phi_0::Function = x-> 0.0,
    num_quad_nodes::Int = 15,
    initial_coeffs::Vector{Float64} = zeros(length(phi_s))
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
    res = optimize(objective, initial_coeffs, BFGS(); autodiff = AutoForwardDiff())

    return res
end


end # module ritz_method
