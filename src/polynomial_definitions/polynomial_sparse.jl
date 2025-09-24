#############################################################################
#############################################################################
#
# This file defines the `PolynomialSparse` type with several operations 
#                                                                               
#############################################################################
#############################################################################

#########################################
# PolynomialSparse type and construction #
#########################################

"""
A Polynomial type - designed to be for polynomials with integer coefficients.

This type utilises a dense representation for a polynomial. This means that
zero terms with degree less than the degree of the polynomial are stored.

E.g, for x^3 + 2x we store in memory:
    [Term(0, 0), Term(2, 1), Term(0, 2), Term(1, 3)]
"""
# struct PolynomialSparse{C,D} <: Polynomial{C,D}
#     terms::Vector{Term{C,D}}

#     # zero polynomial
#     PolynomialSparse{C, D}() where {C, D} = new{C, D}([zero(Term{C, D})])

#     # construct from a vector of terms
#     function PolynomialSparse{C,D}(vt::Vector{Term{C,D}}) where {C,D}
#         # Remove zero terms
#         vt = filter(t -> !iszero(t), vt)
#         isempty(vt) && (vt = [zero(Term{C,D})])

#         max_degree = maximum(t -> t.degree, vt)
#         terms = [zero(Term{C,D}) for i in 0:max_degree]

#         for t in vt
#             terms[t.degree + 1] = t
#         end

#         return new{C,D}(terms)
#     end
# end

struct PolynomialSparse{C,D} <: Polynomial{C,D}
    terms::Heap{Term{C,D}}

    # zero polynomial
    PolynomialSparse{C, D}() where {C, D} = new{C, D}(Heap([zero(Term{C, D})]))

    # construct from a heap of terms
    function PolynomialSparse{C,D}(vt::Vector{Term{C,D}}) where {C,D}
        # Remove zero terms
        vt = filter(t -> !iszero(t), vt)
        isempty(vt) && (vt = [zero(Term{C,D})])

        max_degree = maximum(t -> t.degree, vt)
        terms = [zero(Term{C,D}) for i in 0:max_degree]

        for t in vt
            terms[t.degree + 1] = t
        end

        return new{C,D}(Heap(vt))
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
iterate(p::PolynomialSparse, state=1) = iterate(p.terms, state)

##############################
# Queries about a polynomial #
##############################

"""
The number of terms of the polynomial.
"""
length(p::PolynomialSparse) = length(p.terms) 

"""
The leading term of the polynomial.
"""
# function leading(p::PolynomialSparse)::Term
#     isempty(p.terms) ? zero(Term{Int, Int}) : last(p.terms) 
# end

function leading(p::PolynomialSparse{C,D})::Term{C,D} where {C,D}
    isempty(p.terms) ? zero(Term{C, D}) : last(p.terms) 
end

"""
The term of smallest degree in this polynomial.
"""
# function last(p::PolynomialSparse) 
#     iszero(p) && return leading(p) # zero Term
#     p.terms[findfirst(t -> !iszero(t), p.terms)]
# end
function last(p::PolynomialSparse{C,D}) where {C,D}
    iszero(p) && return leading(p) # zero Term
    p.terms[findfirst(t -> !iszero(t), p.terms)]
end


################################
# Pushing and popping of terms #
################################

"""
Push a new leading term into the polynomial (note - a constant can be pushed onto the zero polynomial).
"""
# function push!(p::PolynomialSparse, t::Term)
#     if t.degree < degree(p) || (t.degree == degree(p) && !iszero(p))
#         error("Cannot push a term $(t) that is not a new leading term (the polynomial had degree $(degree(p)))")
#     elseif iszero(p) && iszero(t.degree) # New constant polynomial
#          p.terms[1] = t
#     else
#         append!(p.terms, zeros(Term{Int, Int}, t.degree - degree(p)-1))
#         push!(p.terms, t)
#     end
#     return p        
# end

function push!(p::PolynomialSparse{C,D}, t::Term{C,D}) where {C,D}
    if t.degree < degree(p) || (t.degree == degree(p) && !iszero(p))
        error("Cannot push a term $(t) that is not a new leading term (the polynomial had degree $(degree(p)))")
    elseif iszero(p) && iszero(t.degree) # New constant polynomial
         p.terms[1] = t
    else
        append!(p.terms, zeros(Term{C, D}, t.degree - degree(p)-1))
        push!(p.terms, t)
    end
    return p        
end


"""
Pop the leading term out of the polynomial. When polynomial is 0, keep popping out 0.
"""
# function pop!(p::PolynomialSparse)::Term 
#     popped_term = pop!(p.terms) #last element popped is leading coefficient

#     while !isempty(p.terms) && iszero(last(p.terms))
#         pop!(p.terms)
#     end

#     if isempty(p.terms)
#         push!(p.terms, zero(Term{Int, Int}))
#     end

#     return popped_term
# end

function pop!(p::PolynomialSparse{C,D})::Term{C,D} where {C,D}
    popped_term = pop!(p.terms) #last element popped is leading coefficient

    while !isempty(p.terms) && iszero(last(p.terms))
        pop!(p.terms)
    end

    if isempty(p.terms)
        push!(p.terms, zero(Term{C, D}))
    end

    return popped_term
end

#################################################################
# Transformation of the polynomial to create another polynomial #
#################################################################

# We again can rely on the Julia type system - the implementations for `Polynomial` will work for `PolynomialSparse`

#################################
# Queries about two polynomials #
#################################

"""
Check if two polynomials are the same.

Note - even though this is done for `Polynomial`, we can override it for `PolynomialSparse`
to leverage Julia's speed with vectors.
"""
# ==(p1::PolynomialSparse, p2::PolynomialSparse)::Bool = p1.terms == p2.terms

function ==(p1::PolynomialSparse{C,D}, p2::PolynomialSparse{C,D}) where {C,D}
    return p1.terms == p2.terms
end



##################################################################
# Operations with two objects where at least one is a polynomial #
##################################################################

# Again, technically the abstract implementations for `Polynomial` will work correctly for `PolynomialSparse`.
# However, this doesn't mean they are particularly efficient - if you wish you can override any particular
# operation and re-implement it for `PolynomialSparse`.`


"""
reworking show because it was bricking my code
"""

function show(io::IO, p::PolynomialSparse{C,D}) where {C,D}
    if iszero(p)
        print(io,"0")
    else
        n = length(p)
        for (i,t) in enumerate(p)
            if !iszero(t)
                print(io, t, i != n ? " + " : "")
            end
        end
    end
end