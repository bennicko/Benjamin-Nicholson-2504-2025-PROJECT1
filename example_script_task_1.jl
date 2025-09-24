using Pkg
Pkg.activate(".")

include("poly_factorization_project.jl")

x = x_poly(PolynomialDense{Int,Int})
f = 2x^3 + 4x^2 - 3x
g = 2x^4 - 4x^2 - 3x + 3
h = 3x^2 + 13x^4 - 8x

@show f+g 
@show f*h
@show g*h

@show derivative(f*g)
@show derivative(f)*g + f*derivative(g)


mod5q, mod5r = div_rem_mod_p(f*h,h,5)
mod17q, mod17r = div_rem_mod_p(f*h,h,17)
mod101q, mod101r = div_rem_mod_p(f*h,h,101)
fmod5, fmod17, fmod101 = mod(f,5), mod(f,17), mod(f,101)


println("f * h ÷ h mod 5 = ", mod5q)
println("fmod5 = ", fmod5)
@show fmod5 == mod5q

println("f * h ÷ h mod 17 = ", mod17q)
println("fmod17 = ", fmod17)
println("fmod101 = ", fmod101)


println("f * h ÷ h mod 101 = ", mod101q)
println("fmod101 = ", fmod101)




println("q1 p5")
p = [5, 11, 13]
for i in p
    println("for p = ", i, ", gcd_mod_p(f*h, g*h, ",i,") = ", gcd_mod_p(f*h, g*h, i))
end