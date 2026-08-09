using Revise
using RitzMethod
import Optim

sample_function(x, y, y_prime) = y_prime^2 + y^2 + 2*x*y

my_domain = Interval(0.0, 2.0)

basis = [x->x*(2-x), x->x^2*(2-x), x->x^3*(2-x)]

# Customize the optimizer algorithm and convergence settings
res = minimizefunctional(
    sample_function, my_domain, basis;
    num_quad_nodes=5,
    method=Optim.BFGS(),
    options=Optim.Options(iterations=50, g_abstol=1e-20, x_reltol=1e-8),
)
Optim.minimizer(res)
