using Pkg
Pkg.activate(".")

include("poly_factorization_project.jl")

x = x_poly(PolynomialDense)
p1 = 100x^2
p2 = 50x^2
p3 = 10^10*x^2
p4 = 10^25*x^2
p5 = 10^32*x^2;


@show p1+p2
@show p2+p3
@show p3+p4
@show p4+p5
@show p5+p1