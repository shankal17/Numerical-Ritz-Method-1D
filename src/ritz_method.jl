module ritz_method

using ArgCheck
using ForwardDiff: derivative
using FastGaussQuadrature: gausslegendre
import LinearAlgebra

export integrate_function, Interval

struct Interval
    low::Float64
    high::Float64
    function Interval(low, high)
        @argcheck low < high
        new(low, high)
    end

end

function integrate_function(f::Function, domain::Interval, num_nodes::Int)
    @argcheck num_nodes > 0
    nodes, weights = gausslegendre(num_nodes)

    # Change of interval
    coeff = (domain.high - domain.low) / 2
    bias = (domain.low + domain.high) / 2
    new_nodes = map(x_i -> coeff * x_i + bias, nodes)

    # Quadrature
    res = coeff * LinearAlgebra.dot(weights, f.(new_nodes))

    return res
end

# TODO: I'm pretty sure that we don't need to evaluate phi's so much, cause we're optimizing over coeffs, and the values of basis
# functions don't change per opimization pass
function eval_y_bar(x::Float64, coeffs::Vector{Float64}, phi_0::Function, phi_s::Vector{Function})::Float64
    val = phi_0(x)
    for i in eachindex(c)
        val += coeffs[i] * phi_s[i](x)
    end

    return val
end

function eval_y_bar_prime(x::Float64, coeffs::Vector{Float64}, phi_0::Function, phi_s::Vector{Function})::Float64
    return ForwardDiff.derivative(x_val -> eval_y_bar(x_val, coeffs, phi_0, phi_s), x)
end

function create_objective(
    integrand::Function,
    domain::Interval,
    phi_s::Vector{Function},
    phi_0::Function = x-> 0.0,
    num_quad_nodes::Int = 15,
    initial_coeffs::Vector{Float64} = zeros(length(phi_s))
)

end

function ritz(
    integrand::Function,
    domain::Interval,
    phi_s::Vector{Function},
    phi_0::Function = x-> 0.0,
    num_quad_nodes::Int = 15,
    initial_coeffs::Vector{Float64} = zeros(length(phi_s))
)

end

end # module ritz_method
