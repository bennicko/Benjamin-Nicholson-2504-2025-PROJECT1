#############################################################################
#############################################################################
#
# This file implements polynomial addition for dense polynomials.
#                                                                               
#############################################################################
#############################################################################

"""
Add a polynomial and a term.
"""
function +(p::PolynomialDense{C,D}, t::Term{C,D}) where {C,D}
    p = deepcopy(p)
    if t.degree > degree(p)
        push!(p, t)
    else
        # +1 accounts for 1-based indexing
        if !iszero(p.terms[t.degree + 1])
            p.terms[t.degree + 1] = Term(p.terms[t.degree + 1].coeff + t.coeff, t.degree)
        else
            p.terms[t.degree + 1] = t
        end
    end
    trim!(p)
    return p
end


# We won't re-implement any of these functions for dense polynomials, the abstract versions will 
# produce the correct result.

# function +(p1::PolynomialDense, p2::PolynomialDense)::PolynomialDense
# +(p::PolynomialDense, n::Int) = p + Term(n,0)
# +(n::Int, p::PolynomialDense) = p + Term(n,0)