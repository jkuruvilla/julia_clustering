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

        mu12 = (k1x*k2x+k1y*k2y+k1z*k2z) / (k1*k2)
        mu23 = (k2x*k3x+k2y*k3y+k2z*k3z) / (k2*k3)
        mu31 = (k3x*k1x+k3y*k1y+k3z*k1z) / (k3*k1)

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

#
# aa = readdlm("../../Anna_Checks/results/bispectrum_2d.txt")
# mu = aa[:,1]
# bk = aa[:,2]
# counts = aa[:,3]

aa = readdlm("../input/flagship_linear_cb_hr_matterpower_z0p0.dat", skipstart=1)
kkk = aa[:,1]
pkkk = aa[:,2]
#plot(wavenumber, pk)

len_pk = length(pkkk)
#damped_pk = zeros((len_pk))

#for i in 1:len_pk
#    damped_pk[i] = pkkk[i] * exp(-(kkk[i]/10)^2)
#end

#pk_inter = Spline1D(kkk, damped_pk, k=3)
pk_inter = Spline1D(kkk, pkkk, k=3)

#
# len_mu = 100#length(mu)
# len_k = 1000 #261 #400


# llim = mu[1]
# ulim = mu[99]

k1_ = 0.0001
#k2_ = collect(0.0001:0.5:100)
cc = readdlm("../results/multi-test-all_Cris.dat")
k2_ = cc[:,1]
b1_cris = cc[:,2]
len_k = length(k2_)

#B_m = zeros((len_k),6) # Array{Float64}(undef, 143, 6)
B_m = zeros((len_k)) # Array{Float64}(undef, 143, 6)

B_full = zeros((len_k))

for i in 1:len_k
    println(i)
    mu_, bk_ = bispec(k1_::Float64, k2_[i]::Float64, pk_inter)
    #B_full[i] = bk_

    bk_inter = Spline1D(mu_, bk_, k=3)
    len_mu_ = length(mu_)
    #for j in 1:6
    _ell = 3 #j - 1
    f = x -> b_multipole(x, _ell, bk_inter)
    llim = mu_[1]
    ulim = mu_[len_mu_]
    #println(llim,ulim)
    I,e = hquadrature(f, llim, ulim, reltol=1e-10, maxevals=10^7)#, abstol=1e-8)

    #B_m[i,j] = I * (2*_ell +1)/2
    B_m[i] = I * (2*_ell +1)/2
    #end
end

# open("../results/B_full_k1_0p0001.txt", "w") do io
#     writedlm(io, [B_full])
# end

open("../results/B_multipole_3_k1_0p0001.txt", "w") do io
    writedlm(io, [k2_, B_m])
end
#
# open("../results/B_multipole_0_k2_0p001.txt", "w") do io
#     writedlm(io, [B_m])
# end
#
# B_m = zeros((len_k))
# _ell = 1
# for i in 1:len_k
#     first_in = len_mu * (i-1) + 1
#     last_in = len_mu * i - 1
#     bk_inter = Spline1D(mu[first_in:last_in], bk[first_in:last_in], k=3)
#     println(i)
#     f = x -> b_multipole(x, _ell, bk_inter)
#     I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-8, maxevals=10^7)#, abstol=1e-8)
#     B_m[i] = I * (2*_ell +1)/2
#     #err_corr[i] = e / (2*_ell +1)/2
# end
#
# open("../results/B_multipole_1_k2_0p001.txt", "w") do io
#     writedlm(io, [B_m])
# end
#
# B_m = zeros((len_k))
# _ell = 2
# for i in 1:len_k
#     first_in = len_mu * (i-1) + 1
#     last_in = len_mu * i - 1
#     bk_inter = Spline1D(mu[first_in:last_in], bk[first_in:last_in], k=3)
#     println(i)
#     f = x -> b_multipole(x, _ell, bk_inter)
#     I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-8, maxevals=10^7)#, abstol=1e-8)
#     B_m[i] = I * (2*_ell +1)/2
#     #err_corr[i] = e / (2*_ell +1)/2
# end
# open("../results/B_multipole_2_k2_0p001.txt", "w") do io
#     writedlm(io, [B_m])
# end
#
#
# B_m = zeros((len_k))
# _ell = 3
# for i in 1:len_k
#     first_in = len_mu * (i-1) + 1
#     last_in = len_mu * i - 1
#     bk_inter = Spline1D(mu[first_in:last_in], bk[first_in:last_in], k=3)
#     println(i)
#     f = x -> b_multipole(x, _ell, bk_inter)
#     I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-8, maxevals=10^7)#, abstol=1e-8)
#     B_m[i] = I * (2*_ell +1)/2
#     #err_corr[i] = e / (2*_ell +1)/2
# end
# open("../results/B_multipole_3_k2_0p001.txt", "w") do io
#     writedlm(io, [B_m])
# end
#
# B_m = zeros((len_k))
# _ell = 4
# for i in 1:len_k
#     first_in = len_mu * (i-1) + 1
#     last_in = len_mu * i - 1
#     bk_inter = Spline1D(mu[first_in:last_in], bk[first_in:last_in], k=3)
#     println(i)
#     f = x -> b_multipole(x, _ell, bk_inter)
#     I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-8, maxevals=10^7)#, abstol=1e-8)
#     B_m[i] = I * (2*_ell +1)/2
#     #err_corr[i] = e / (2*_ell +1)/2
# end
# open("../results/B_multipole_4_k2_0p001.txt", "w") do io
#     writedlm(io, [B_m])
# end
#
# B_m = zeros((len_k))
# _ell = 5
# for i in 1:len_k
#     first_in = len_mu * (i-1) + 1
#     last_in = len_mu * i - 1
#     bk_inter = Spline1D(mu[first_in:last_in], bk[first_in:last_in], k=3)
#     println(i)
#     f = x -> b_multipole(x, _ell, bk_inter)
#     I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-8, maxevals=10^7)#, abstol=1e-8)
#     B_m[i] = I * (2*_ell +1)/2
#     #err_corr[i] = e / (2*_ell +1)/2
# end
# open("../results/B_multipole_5_k2_0p001.txt", "w") do io
#     writedlm(io, [B_m])
# end
