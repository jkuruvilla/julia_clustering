import Dierckx.Spline1D, Dierckx.evaluate
import Cubature.hquadrature


π = 3.141592653589793

function Hubble(a::Float64, H0::Float64, OmegaM::Float64, OmegaL::Float64)
    H = H0 * sqrt(OmegaM/a^3 + OmegaL)
    return H
end

function growth_factor(a::Float64, OmegaM::Float64, OmegaL::Float64)
    integrand = 1/(OmegaM/a + OmegaL*a^2)^(3/2)
    return integrand
end

Hubblezero = 67.74
MatterDensity = 0.3089
LambdaDensity = 0.6911

scalefactor50 = 1/51
scalefactor0 = 1.0

f = x -> growth_factor(x, MatterDensity, LambdaDensity)
growthfactor50, e50 = hquadrature(f, 0.0, scalefactor50)

growthfactor50 *= (5/2 * MatterDensity * Hubble(scalefactor50, Hubblezero, MatterDensity, LambdaDensity) / Hubblezero)
println(growthfactor50)

f = x -> growth_factor(x, MatterDensity, LambdaDensity)
growthfactor0, e0 = hquadrature(f, 0.0, 1)

growthfactor0 *= (5/2 * MatterDensity * Hubble(scalefactor0, Hubblezero, MatterDensity, LambdaDensity) / Hubblezero)
println(growthfactor0)

println((growthfactor0/growthfactor50))
println((growthfactor0/growthfactor50)^2)
