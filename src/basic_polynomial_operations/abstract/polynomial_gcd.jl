#############################################################################
#############################################################################
#
# This file implements polynomial GCD for abstract polynomials.
#                                                                               
#############################################################################
#############################################################################

"""
The extended euclid algorithm for polynomials (of the same concrete subtype) modulo prime.

Note, when working with polynomials over Z[x] we cannot compute the extended euclidean algorithm (EEA) directly.
This is because we do not have division in Z[x]. Hence, we must drop to the field Zp[x] and compute the EEA
with respect to some prime. 
"""
function extended_euclid_alg_mod_p(a::P, b::P, prime::Integer) where {P <: Polynomial}
    old_r, r = mod(a, prime), mod(b, prime)
    old_s, s = one(P), zero(P)
    old_t, t = zero(P), one(P)

    while !iszero(mod(r,prime))
        q = div_mod_p(old_r, r, prime)
        old_r, r = r, mod(old_r - q*r, prime)
        old_s, s = s, mod(old_s - q*s, prime)
        old_t, t = t, mod(old_t - q*t, prime)
    end
    g, s, t = old_r, old_s, old_t

    @assert mod(s*a + t*b - g, prime) == 0
    return g, s, t  
end

function extended_euclid_alg_mod_p(a::P, b::P, prime::Integer) where {C,D,P <: Polynomial{C,D}}
    primeC = C(prime) #avoid typing issues
    old_r, r = mod(a, primeC), mod(b, primeC)
    old_s, s = one(P), zero(P)
    old_t, t = zero(P), one(P)

    while !iszero(mod(r,primeC))
        q = div_mod_p(old_r, r, primeC)
        old_r, r = r, mod(old_r - q*r, primeC)
        old_s, s = s, mod(old_s - q*s, primeC)
        old_t, t = t, mod(old_t - q*t, primeC)
    end
    g, s, t = old_r, old_s, old_t

    @assert mod(s*a + t*b - g, primeC) == 0
    return g, s, t  
end

"""
The GCD of two polynomials (of the same concrete subtype) modulo prime.
"""
gcd_mod_p(a::P, b::P, prime::Integer) where {P <: Polynomial} = extended_euclid_alg_mod_p(a, b, prime) |> first

gcd_mod_p(a::P, b::P, prime::Integer) where {C,D,P <: Polynomial{C,D}} = extended_euclid_alg_mod_p(a, b, C(prime)) |> first
