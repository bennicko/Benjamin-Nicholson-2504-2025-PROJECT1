using Pkg
Pkg.activate(".")

include("poly_factorization_project.jl")

x = x_poly(PolynomialDense)

p1 = 2x^2
p2 = 3x^2
p3 = 4x^2
p4 = 5x^2
p5 = 10^32*x^2
@show p1+p2+p3+p4+p5
