#############################################################################
#############################################################################
#
# This file implements polynomial addition for abstract polynomials.
#                                                                               
#############################################################################
#############################################################################

"""
Add a polynomial and a term.

This must be overridden by concrete subtypes.
"""
function +(p::Polynomial{C,D}, t::Term{C,D}) where {C,D}
    not_implemented_error(p, "Polynomial + Term")
end

+(t::Term{C,D}, p::Polynomial{C,D}) where {C,D} = p + t


"""
Add two polynomials of the same concrete subtype.

Note: This operation may be slow for some concrete subtypes. You may wish to override this to factor
in the details of your polynomial representation when implementing your concrete subtype.
"""
function +(p1::Polynomial{C,D}, p2::Polynomial{C,D}) where {C,D}
    p = deepcopy(p1)
    for t in p2.terms
        p += t
    end
    return p
end


"""
Add a polynomial and an integer.
"""
+(p::Polynomial{C,D}, n::C) where {C,D} = p + Term(n, zero(D))
+(n::C, p::Polynomial{C,D}) where {C,D} = p + Term(n, zero(D))
