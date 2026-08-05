using Revise
using ritz_method

sample_function(x) = x^3

my_domain = Interval(-3.0, 2.0)

res = integrate_function(sample_function, my_domain, 3)

println(res)
