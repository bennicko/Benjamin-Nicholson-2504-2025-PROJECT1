using Pkg
Pkg.activate(".")

include("poly_factorization_project.jl")


#q2 p1

x = x_poly(PolynomialDense)
p1 = 10x^2
p2 = 20x^2
p3 = 30x^2
p4 = 40x^2
p5 = (10^64)*x^2

@show p1+p2+p3+p4+p5

