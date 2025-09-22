using Pkg;
# To be able to run this, have the worksheet inside the project repository, it should not be part of your assignment submission. 
Pkg.activate(".");
# Pkg.instantiate(); #This should be uncommented when the project is first run
include("poly_factorization_project.jl");
     
struct RationalFunction
    numerator::Polynomial
    denominator::Polynomial
end

x = x_poly(PolynomialDense)
r1 = RationalFunction(5x^2-3x+4, 6x^4-2x+5)
r2 = RationalFunction(-3x+4, 2x^2-2x+5)



function show(io::IO, r::RationalFunction)
    println(io, r.numerator)
    num_chars = max(length(string(r.numerator)),length(string(r.denominator)))
    println(io, "-"^num_chars)
    println(io,r.denominator)
end

*(rf1::RationalFunction, rf2::RationalFunction) = 
        RationalFunction(rf1.numerator * rf2.numerator, rf1.denominator * rf2.denominator)

function derivative(r::RationalFunction)
    n = r.numerator
    d = r.denominator

    RationalFunction(derivative(n)*d - n*derivative(d), d^2)
end

function +(rf1::RationalFunction, rf2::RationalFunction)
    a, b = rf1.numerator, rf1.denominator
    c, d = rf2.numerator, rf2.denominator
    common = b*d
    return RationalFunction(a*d + c*b, common)
end

clean(n::Int) = abs(2*(n÷2))
[(n, clean(n)) for n=-5:5] |> println


clean(t::Term) = Term(clean(t.coeff),t.degree)


cleaned = clean(Term(1,3))
     
function clean(p::PolynomialDense)
    p_out = PolynomialDense()
    terms = deepcopy(p).terms
    for t in terms
        clean_t = clean(t)
        !iszero(clean_t) && push!(p_out,clean(t))
    end
    return p_out
end

     

using Random; Random.seed!(0)
p = rand(PolynomialDense) + 1x^50


clean(p)
     
clean(r::RationalFunction) = RationalFunction(clean(r.numerator), clean(r.denominator))

struct PointFloat_attempt1
    x::Float64
    y::Float64
end

struct PointInt_attempt1
    x::Int
    y::Int
end

# @show norm2.(PointFloats)
# @show norm2.(IntFloats)

abstract type Point2 end

struct PointFloat_attempt2 <: Point2 #error because we're changing the type structure
    x::Float64
    y::Float64
end

struct PointInt_attempt2  <: Point2
    x::Int
    y::Int

end



PointFloats = [PointFloat_attempt1(3rand(),3rand()) for _ in 1:10]
IntFloats = [PointInt_attempt1(i,j) for i in 1:3 for j in 1:3]

FloatPoints2 = [PointFloat_attempt2(3rand(), 3rand()) for _ in 1:10] #10 random floating points in [0,3]×[0,3]
IntPoints2 = [PointInt_attempt2(i,j) for i in 1:3 for j in 1:3] #9 integer points at x=(1,2,3), y=(1,2,3)

function norm(point::PointFloat_attempt1)
    return sqrt(point.x^2 + point.y^2)
end

function norm(point::PointInt_attempt1)
    return sqrt(point.x^2 + point.y^2)
end

function norm2(point::T) where {T <: Union{PointFloat_attempt1, PointInt_attempt1}}
    return sqrt(point.x^2 + point.y^2)
end


function norm(point::T) where {T <: Point2}
    println("Finding norm of an integer point")
    return sqrt(point.x^2 + point.y^2)
end
norm(IntPoints2[1])

using Plots
function plotpoints(points::Vector{T}; plt=plot()) where {T <: Point2}
    plot!(plt, xlabel="x", ylabel = "y", legend=:none)
    scatter!(
        [p.x for p in points],
        [p.y for p in points],
        color = (T==PointInt_attempt2) ? :lightblue : :orange
    )
end
plt = plotpoints(FloatPoints2)
plt = plotpoints(IntPoints2; plt=plt)



struct Point_attempt3{T <: Real} <: Point2 #With this signature we're creating a new concrete type Point_attempt3{T} for every subtype of Real. All the new concrete types will be subtypes of Point2
    x::T
    y::T
end

FloatPoints3 = [Point_attempt3{Float64}(3rand(), 3rand()) for _ in 1:10] #10 random floating points in [0,3]×[0,3]
IntPoints3 = [Point_attempt3{Int}(i,j) for i in 1:3 for j in 1:3] #9 integer points at x=(1,2,3), y=(1,2,3)

function plotpoints(points::Point_attempt3{Int}; plt=plot())
    plot!(plt, xlabel="x", ylabel="y", legend=:none)
    scatter!(
        [p.x for p in points], #unpack all x values from points
        [p.y for p in points], #unpack all y values from points
        color = :lightblue 
    )
end

plt = plotpoints(FloatPoints3)
plotpoints(IntPoints3; plt=plt) 