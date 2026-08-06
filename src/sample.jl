using Revise
using ritz_method
using Optim: minimizer
sample_function(x, y, y_prime) = y_prime^2 - y^2 - 2*x*y

my_domain = Interval(0.0, 1.0)

res = ritz(sample_function, my_domain, [x->x*(1-x), x->x^2*(1-x)]; num_quad_nodes=5)
minimizer(res)
