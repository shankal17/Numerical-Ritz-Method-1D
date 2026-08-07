using Revise
using RitzMethod
using Optim: minimizer, BFGS, Options

sample_function(x, y, y_prime) = y_prime^2 - y^2 - 2*x*y

my_domain = Interval(0.0, 1.0)
basis = [x->x*(1-x), x->x^2*(1-x), x->x^3*(1-x)]

# Customize the optimizer algorithm and convergence settings
res_custom = minimizefunctional(
    sample_function, my_domain, basis;
    num_quad_nodes=5,
    method=BFGS(),
    options=Options(iterations=50, g_abstol=1e-20, x_reltol=1e-8),
)
minimizer(res_custom)
