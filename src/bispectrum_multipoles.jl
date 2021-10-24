import Dierckx.Spline1D, Dierckx.evaluate
import Cubature.hquadrature
import DelimitedFiles.readdlm, DelimitedFiles.writedlm
import LegendrePolynomials.Pl

π = 3.141592653589793

function b_multipole(mu::Float64, ell::Real, bk_int)
    integr = bk_int(mu) * Pl(mu, ell)
    return integr
end

aa = readdlm("../../Anna_Checks/results/bispectrum_2d.txt")
mu = aa[:,1]
bk = aa[:,2]
counts = aa[:,3]

aa = readdlm("../input/flagship_linear_cb_hr_matterpower_z0p0.dat", skipstart=1)
k = aa[:,1]
pkkk = aa[:,2]
#plot(wavenumber, pk)

len_mu = 100#length(mu)
len_k = 1000 #261 #400


llim = mu[1]
ulim = mu[99]

B_m = zeros((len_k))
_ell = 0
for i in 1:len_k
    first_in = (len_mu * (i-1)) + 1
    last_in = len_mu * i - 1
    bk_inter = Spline1D(mu[first_in:last_in], bk[first_in:last_in], k=3)
    println(i)
    f = x -> b_multipole(x, _ell, bk_inter)
    I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-8, maxevals=10^7)#, abstol=1e-8)
    B_m[i] = I * (2*_ell +1)/2
    #err_corr[i] = e / (2*_ell +1)/2
end

open("../results/B_multipole_0_k2_0p001.txt", "w") do io
    writedlm(io, [B_m])
end

B_m = zeros((len_k))
_ell = 1
for i in 1:len_k
    first_in = len_mu * (i-1) + 1
    last_in = len_mu * i - 1
    bk_inter = Spline1D(mu[first_in:last_in], bk[first_in:last_in], k=3)
    println(i)
    f = x -> b_multipole(x, _ell, bk_inter)
    I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-8, maxevals=10^7)#, abstol=1e-8)
    B_m[i] = I * (2*_ell +1)/2
    #err_corr[i] = e / (2*_ell +1)/2
end

open("../results/B_multipole_1_k2_0p001.txt", "w") do io
    writedlm(io, [B_m])
end

B_m = zeros((len_k))
_ell = 2
for i in 1:len_k
    first_in = len_mu * (i-1) + 1
    last_in = len_mu * i - 1
    bk_inter = Spline1D(mu[first_in:last_in], bk[first_in:last_in], k=3)
    println(i)
    f = x -> b_multipole(x, _ell, bk_inter)
    I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-8, maxevals=10^7)#, abstol=1e-8)
    B_m[i] = I * (2*_ell +1)/2
    #err_corr[i] = e / (2*_ell +1)/2
end
open("../results/B_multipole_2_k2_0p001.txt", "w") do io
    writedlm(io, [B_m])
end


B_m = zeros((len_k))
_ell = 3
for i in 1:len_k
    first_in = len_mu * (i-1) + 1
    last_in = len_mu * i - 1
    bk_inter = Spline1D(mu[first_in:last_in], bk[first_in:last_in], k=3)
    println(i)
    f = x -> b_multipole(x, _ell, bk_inter)
    I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-8, maxevals=10^7)#, abstol=1e-8)
    B_m[i] = I * (2*_ell +1)/2
    #err_corr[i] = e / (2*_ell +1)/2
end
open("../results/B_multipole_3_k2_0p001.txt", "w") do io
    writedlm(io, [B_m])
end

B_m = zeros((len_k))
_ell = 4
for i in 1:len_k
    first_in = len_mu * (i-1) + 1
    last_in = len_mu * i - 1
    bk_inter = Spline1D(mu[first_in:last_in], bk[first_in:last_in], k=3)
    println(i)
    f = x -> b_multipole(x, _ell, bk_inter)
    I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-8, maxevals=10^7)#, abstol=1e-8)
    B_m[i] = I * (2*_ell +1)/2
    #err_corr[i] = e / (2*_ell +1)/2
end
open("../results/B_multipole_4_k2_0p001.txt", "w") do io
    writedlm(io, [B_m])
end

B_m = zeros((len_k))
_ell = 5
for i in 1:len_k
    first_in = len_mu * (i-1) + 1
    last_in = len_mu * i - 1
    bk_inter = Spline1D(mu[first_in:last_in], bk[first_in:last_in], k=3)
    println(i)
    f = x -> b_multipole(x, _ell, bk_inter)
    I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-8, maxevals=10^7)#, abstol=1e-8)
    B_m[i] = I * (2*_ell +1)/2
    #err_corr[i] = e / (2*_ell +1)/2
end
open("../results/B_multipole_5_k2_0p001.txt", "w") do io
    writedlm(io, [B_m])
end
