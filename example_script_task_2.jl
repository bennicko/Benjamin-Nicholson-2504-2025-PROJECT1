using Pkg
Pkg.activate(".")

include("poly_factorization_project.jl")


a = Term{BigInt,Int}(BigInt(1),1)
b = Term{BigInt,Int}(BigInt(2),2)
p = PolynomialDense{BigInt,Int}([a])



x = x_poly(PolynomialDense)
p1 = 100x^2
p2 = 50x^2
p3 = 10^10*x^2
p4 = 10^25*x^2
p5 = 10^32*x^2


# @show p1+p2
# @show p2+p3
# @show p3+p4
# @show p4+p5
# @show p5+p1;


p6 = Term{BigInt,Int}(BigInt(100),2)
p7 = Term{BigInt,Int}(BigInt(50),2)
p8 = Term{BigInt,Int}(BigInt(10^10),2)
p9 = Term{BigInt,Int}(BigInt(10^25),2)
p10 = Term{BigInt,Int}(BigInt(10^32),2)

# p1 = Term{Int,Int}(Int(100),2)
# p2 = Term{Int,Int}(Int(50),2)
# p3 = Term{Int,Int}(Int(10^10),2)
# p4 = Term{Int,Int}(Int(10^25),2)
# p5 = Term{Int,Int}(Int(10^32),2)


@show PolynomialDense{BigInt, Int}([p6])+PolynomialDense{BigInt, Int}([p7])
@show p6+p7
# @show p7+p8
# @show p8+p9
@show PolynomialDense{BigInt, Int}([p9])+PolynomialDense{BigInt, Int}([p10])
@show p9+p10
# @show p10+p6;
