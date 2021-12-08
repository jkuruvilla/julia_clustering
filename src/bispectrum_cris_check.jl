import Dierckx.Spline1D, Dierckx.evaluate
import Cubature.hquadrature
import DelimitedFiles.readdlm, DelimitedFiles.writedlm
import LegendrePolynomials.Pl

π = 3.141592653589793

function b_multipole(mu::Float64, ell::Real, bk_int)
    integr = bk_int(mu) * Pl(mu, ell)
    return integr
end

function bispec(k1::Float64, k2::Float64, pk_int)
    mu = collect(-1.0:0.05:1.0)
    k1x, k1y, k1z = 0,0,k1
    len_mu = length(mu)
    B_averaged = zeros((len_mu))

    for i in 1:len_mu
        k2x, k2y, k2z = 0, k2*sqrt(1-mu[i]^2), k2*mu[i]
        k3x, k3y, k3z = 0, -k2*sqrt(1-mu[i]^2), -(k2*mu[i] + k1)
        k3 = sqrt((k3x^2 + k3y^2 + k3z^2))
        println(k3)
        mu12 = (k1x*k2x+k1y*k2y+k1z*k2z) / (k1*k2)
        mu23 = (k2x*k3x+k2y*k3y+k2z*k3z) / (k2*k3)
        mu31 = (k3x*k1x+k3y*k1y+k3z*k1z) / (k3*k1)
        println(mu[i],mu12)

        B_first = B_kernel(k1, k2, mu12, pk_int)
        B_second = B_kernel(k2, k3, mu23, pk_int)
        B_third = B_kernel(k3, k1, mu31, pk_int)

        B_averaged[i] = B_first + B_second + B_third
    end

    return mu, B_averaged
end

function B_kernel(k1::Float64, k2::Float64, mu::Float64, pk_int)
	return (10/7 + (k1/k2 + k2/k1)*mu + 4/7*mu^2) * pk_int(k1)*pk_int(k2)
end
aa = readdlm("../input/flagship_linear_cb_hr_matterpower_z0p0.dat", skipstart=1)
kkk = aa[:,1]
pkkk = aa[:,2]
#plot(wavenumber, pk)

len_pk = length(pkkk)
damped_pk = zeros((len_pk))

for i in 1:len_pk
    damped_pk[i] = pkkk[i] * exp(-(kkk[i]/10)^2)
end

pk_inter = Spline1D(kkk, damped_pk, k=3)

muu, bii = bispec(1.0, 1.0, pk_inter)
println(bii)

open("../results/B_full_k1_1_k2_1.txt", "w") do io
    writedlm(io, [muu,bii])
end
