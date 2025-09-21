#############################################################################
#############################################################################
#
# This file defines the `PolynomialDense` type with several operations 
#                                                                               
#############################################################################
#############################################################################

#########################################
# PolynomialDense type and construction #
#########################################

"""
A Polynomial type - designed to be for polynomials with integer coefficients.

This type utilises a dense representation for a polynomial. This means that
zero terms with degree less than the degree of the polynomial are stored.

E.g, for x^3 + 2x we store in memory:
    [Term(0, 0), Term(2, 1), Term(0, 2), Term(1, 3)]
"""
struct PolynomialDense{C, D} <: Polynomial{C, D}
    # FIXME in future make storing correct degree for zero terms an invariant - Mittun
    #A zero packed vector of terms
    #Terms are assumed to be in order with first term having degree 0, second degree 1, and so fourth
    #until the degree of the polynomial. The leading term (i.e. last) is assumed to be non-zero except 
    #for the zero polynomial where the vector is of length 1.
    terms::Vector{Term{C, D}}
    
    #Inner constructor of 0 polynomial
    PolynomialDense() = new{C,D}([zero(Term{Int, Int})])

    #Inner constructor of polynomial based on arbitrary list of terms
    function PolynomialDense(vt::Vector{Term{Int, Int}})

        #Filter the vector so that there is not more than a single zero term
        vt = filter((t)->!iszero(t), vt)
        if isempty(vt)
            vt = [zero(Term{Int, Int})]
        end

        max_degree = maximum((t)->t.degree, vt)
        terms = [zero(Term{Int, Int}) for i in 0:max_degree] #First set all terms with zeros

        #now update based on the input terms
        for t in vt
            terms[t.degree + 1] = t #+1 accounts for 1-indexing
        end
        return new{C,D}(terms)
    end
end


###########
# Display #
###########

# We don't need to override the `show` function from this section.
# The implementation from `Polynomial` is sufficient for our needs.`

##############################################
# Iteration over the terms of the polynomial #
##############################################

"""
Allows to do iteration over the non-zero terms of the polynomial. This implements the iteration interface.
"""
function iterate(p::PolynomialDense{C,D}, state=1) where {C,D}
    return iterate(p.terms, state)
end


##############################
# Queries about a polynomial #
##############################

"""
The number of terms of the polynomial.
"""
length(p::PolynomialDense{C,D}) where {C,D} = length(p.terms)


"""
The leading term of the polynomial.
"""
function leading(p::PolynomialDense{C,D}) where {C,D}
    isempty(p.terms) ? zero(Term{Int, Int}) : last(p.terms)
end


"""
The term of smallest degree in this polynomial.
"""
function last(p::PolynomialDense{C,D}) where {C,D} 
    iszero(p) && return leading(p)  # zero Term
    p.terms[findfirst(t -> !iszero(t), p.terms)]
end



################################
# Pushing and popping of terms #
################################

"""
Push a new leading term into the polynomial (note - a constant can be pushed onto the zero polynomial).
"""
function push!(p::PolynomialDense{C,D}, t::Term{C,D}) where {C,D}
    if t.degree < degree(p) || (t.degree == degree(p) && !iszero(p))
        error("Cannot push a term $(t) that is not a new leading term (the polynomial had degree $(degree(p)))")
    elseif iszero(p) && iszero(t.degree)  # New constant polynomial
        p.terms[1] = t
    else
        append!(p.terms, [zero(Term{C,D}) for _ in 1:(t.degree - degree(p) - 1)])
        push!(p.terms, t)
    end
    return p        
end


"""
Pop the leading term out of the polynomial. When polynomial is 0, keep popping out 0.
"""
function pop!(p::PolynomialDense{C,D}) where {C,D}
    popped_term = pop!(p.terms)  # last element popped is leading term

    while !isempty(p.terms) && iszero(last(p.terms))
        pop!(p.terms)
    end

    if isempty(p.terms)
        push!(p.terms, zero(Term{C,D}))
    end

    return popped_term
end


#################################################################
# Transformation of the polynomial to create another polynomial #
#################################################################

# We again can rely on the Julia type system - the implementations for `Polynomial` will work for `PolynomialDense`

#################################
# Queries about two polynomials #
#################################

"""
Check if two polynomials are the same.

Note - even though this is done for `Polynomial`, we can override it for `PolynomialDense`
to leverage Julia's speed with vectors.
"""
==(p1::PolynomialDense{C,D}, p2::PolynomialDense{C,D}) where {C,D} = p1.terms == p2.terms




##################################################################
# Operations with two objects where at least one is a polynomial #
##################################################################

# Again, technically the abstract implementations for `Polynomial` will work correctly for `PolynomialDense`.
# However, this doesn't mean they are particularly efficient - if you wish you can override any particular
# operation and re-implement it for `PolynomialDense`.`