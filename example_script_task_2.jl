using Pkg
Pkg.activate(".")

include("poly_factorization_project.jl")


x1 = x_poly(PolynomialDense{BigInt,Int})
x2 = x_poly(PolynomialDense{Int,Int})

p1intint = (10*10^18)*x2^2
p2intint = (10*10^17)*x2^2

p1bigintint = big(10)*10^18*x1^2
p2bigintint = big(10)*10^17*x1^2

@show p2intint+p1intint
@show p1bigintint+p2bigintint
