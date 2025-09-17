using Pkg


Pkg.activate(".")

include("poly_factorization_project.jl")

function pretty_print(poly_str::String)
    # Step 1: Parse the terms
    terms = split(poly_str, r"\s*\+\s*")  # split on '+', allowing spaces
    parsed_terms = []

    for term in terms
        term = replace(term, " " => "")
        # Extract coefficient and power
        m = match(r"([+-]?\d+)⋅x\^(\d+)", term)
        if m !== nothing
            coeff = parse(Int, m.captures[1])
            power = parse(Int, m.captures[2])
            push!(parsed_terms, (power, coeff))
        else
            # Handle constants like "5"
            if occursin("x", term)
                error("Unrecognized term format: $term")
            else
                coeff = parse(Int, term)
                push!(parsed_terms, (0, coeff))
            end
        end
    end

    # Step 2: Sort by descending degree
    sorted_terms = sort(parsed_terms, by = x -> -x[1])

    # Step 3: Build pretty string
    result = ""
    for (i, (power, coeff)) in enumerate(sorted_terms)
        if coeff == 0
            continue
        end

        # Determine sign
        sign_str = ""
        if coeff < 0
            sign_str = i == 1 ? "-" : " - "
        elseif coeff > 0 && i > 1
            sign_str = " + "
        end

        abs_coeff = abs(coeff)

        # Format term
        term_str = ""
        if power == 0
            term_str = string(abs_coeff)
        elseif power == 1
            term_str = abs_coeff == 1 ? "x" : "$(abs_coeff)x"
        else
            term_str = abs_coeff == 1 ? "x^$power" : "$(abs_coeff)x^$power"
        end

        result *= sign_str * term_str
    end

    return(result)  
end




#q1 p1
x = x_poly(PolynomialDense)
f = x^4-3x^2
g = x^2 + 5x^3
h = -3x^3 + 2x + (-5)

println("q1 p2")
println("f + g: ", pretty_print(string(f+g)))
println("f * h: ", pretty_print(string(f*h)))
println("g * h: ", pretty_print(string(g*h)))


println("q1 p3")
println("derivative f * g = ", pretty_print(string(derivative(f * g))))
println("product rule; f'g + g'f = ", pretty_print(string(derivative(f) * g + f * derivative(g))))

println("q1 p4")

mod5q, mod5r = div_rem_mod_p(f*h,h,5)
mod17q, mod17r = div_rem_mod_p(f*h,h,17)
mod101q, mod101r = div_rem_mod_p(f*h,h,101)
fmod5, fmod17, fmod101 = mod(f,5), mod(f,17), mod(f,101)


println("f * h ÷ h mod 5 = ", pretty_print(string(mod5q)))
println("fmod3 = ", pretty_print(string(fmod5)))
@show fmod5 == mod5q

println("f * h ÷ h mod 17 = ", pretty_print(string(mod17q)))
println("fmod17 = ", pretty_print(string(fmod17)))


println("f * h ÷ h mod 101 = ", pretty_print(string(mod101q)))
println("fmod101 = ", pretty_print(string(fmod101)))




println("q1 p5")
p = [5, 11, 13]
for i in p
    println("for p = ", i, ", gcd_mod_p(f*h, g*h, ",i,") = ", pretty_print(string(gcd_mod_p(f*h, g*h, i))))
end