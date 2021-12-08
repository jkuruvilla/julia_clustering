import Dierckx.Spline1D, Dierckx.evaluate, Dierckx.derivative
import Cubature.hquadrature
import DelimitedFiles.readdlm, DelimitedFiles.writedlm
import LegendrePolynomials.Pl

π = 3.141592653589793

function correlation_integrand(kk::Float64, rr::Float64, pk_int)
    pk_ev = evaluate(pk_int, kk)
    _va = kk * sin(kk*rr)/(rr)
    integr = pk_ev * _va
    return integr
end


function xi_oneplus(kk::Float64, rr::Float64, pk_int)
    arg_ = kk*rr
    pk_ev = evaluate(pk_int, kk)
    _va = kk^3 * ( (sin(arg_)/arg_^2) - (cos(arg_)/arg_) )
    integr = pk_ev * _va
    return integr
end

function xi_oneminus(kk::Float64, rr::Float64, pk_int)
    arg_ = kk*rr
    pk_ev = evaluate(pk_int, kk)
    _va = kk * ( (sin(arg_)/arg_^2) - (cos(arg_)/arg_) )
    integr = pk_ev * _va
    return integr
end

function xi_two(kk::Float64, rr::Float64, pk_int)
    arg_ = kk*rr
    pk_ev = evaluate(pk_int, kk)
    _va = kk^2 * ( ((3/arg_^2 - 1) *(sin(arg_)/arg_)) - (3*cos(arg_)/arg_^2) )
    integr = pk_ev * _va
    return integr
end

function Slepian_precyclic(r1::Float64, r2::Float64, _xi, _xim, _xip, _xitwo)
    pc0 = 34/21 * _xi(r1) * _xi(r2)
    pc1 = -1 * ((_xim(r1)*_xip(r2)) + (_xim(r2)*_xip(r1)))
    pc2 = 8/21 * (_xitwo(r1) * _xitwo(r2))
    return pc0, pc1, pc2
end


aa = readdlm("../input/flagship_linear_cb_hr_matterpower_z0p0.dat", skipstart=1)
wavenumber = aa[:,1]
pk = aa[:,2]

#plot(wavenumber, pk)

len_pk = length(pk)
damped_pk = zeros((len_pk))
for i in 1:len_pk
    damped_pk[i] = pk[i] * exp(-(wavenumber[i]/10)^2)
end


pk_inter = Spline1D(wavenumber, damped_pk, k=3)

rrr = collect(0.01:0.1:160)
# rrr = readdlm("rbins_anna.txt")
# rrr = readdlm("rbins_fftlog_anna.txt")
len_r00 = length(rrr) # 300
println(len_r00)

llim = wavenumber[1]
ulim = wavenumber[len_pk] #100

corr = zeros((len_r00))
corr_oneplus = zeros((len_r00))
corr_oneminus = zeros((len_r00))
corr_two = zeros((len_r00))

err_corr = zeros((len_r00))
err_corr_oneplus = zeros((len_r00))
err_corr_oneminus = zeros((len_r00))
err_corr_two = zeros((len_r00))

_rrr = zeros((len_r00))

rr = 0.0

for i in 1:len_r00
    f = x -> correlation_integrand(x, rrr[i], pk_inter)
    I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-8, maxevals=10^7)#, abstol=1e-8)
    corr[i] = I / (2*pi^2)
    err_corr[i] = e / (2*pi^2)

    f = x -> xi_oneplus(x, rrr[i], pk_inter)
    I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-8, maxevals=10^7)#, abstol=1e-8)
    corr_oneplus[i] = I / (2*pi^2) #* (4*pi)
    err_corr_oneplus[i] = e / (2*pi^2) #* (4*pi)

    f = x -> xi_oneminus(x, rrr[i], pk_inter)
    I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-8, maxevals=10^7)#, abstol=1e-8)
    corr_oneminus[i] = I / (2*pi^2) #* (4*pi)
    err_corr_oneminus[i] = e / (2*pi^2) #* (4*pi)

    f = x -> xi_two(x, rrr[i], pk_inter)
    I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-8, maxevals=10^7)#, abstol=1e-8)
    corr_two[i] = I / (2*pi^2) #* (4*pi)
    err_corr_two[i] = e / (2*pi^2) #* (4*pi)

    _rrr[i] = rrr[i]
    #println(i)
end

xi_inter = Spline1D(_rrr, corr, k=3)
xi_onep_inter = Spline1D(_rrr, corr_oneplus, k=3)
xi_onem_inter = Spline1D(_rrr, corr_oneminus, k=3)
xi_two_inter = Spline1D(_rrr, corr_two, k=3)

mu = collect(-1:0.1:1)
lenmu = length(mu)
llim_mu = mu[1]
ulim_mu = mu[lenmu]

tri_rr = readdlm("../input/triangles_binSize_5.dat", skipstart=1)
_r12 = tri_rr[:,1]
_r23 = tri_rr[:,2]
_r31 = tri_rr[:,3]

lenr12_ = length(_r12)
_numell = 20
len_multipoles = (lenr12_ * _numell) + _numell
println(len_multipoles)

slepian_multi = zeros((len_multipoles))
err_slepian_multi = zeros((len_multipoles))

for i in 1:lenr12_
    #println(i)
    pc0, pc1, pc2 = Slepian_precyclic(_r12[i], _r23[i], xi_inter, xi_onem_inter, xi_onep_inter, xi_two_inter)
    pc_sum12 = zeros((lenmu))
    pc_sum23 = zeros((lenmu))
    pc_sum31 = zeros((lenmu))
    for j in 1:lenmu
        rrr31 = sqrt(_r12[i]^2 + _r23[i]^2 - (2*mu[j]*_r12[i]*_r23[i]))
        mu1 = mu[j]
        mu2 = (_r23[i]^2+rrr31^2-_r12[i]^2)/(2*_r23[i]*rrr31)
        mu3 = (rrr31^2+_r12[i]^2-_r23[i]^2)/(2*rrr31*_r12[i])
        pc_sum12[j] = pc0 *  Pl(mu1, 0) + pc1 *  Pl(mu1, 1) + pc2 *  Pl(mu1, 2)
        pc_sum23[j] = pc0 *  Pl(mu2, 0) + pc1 *  Pl(mu2, 1) + pc2 *  Pl(mu2, 2)
        pc_sum31[j] = pc0 *  Pl(mu3, 0) + pc1 *  Pl(mu3, 1) + pc2 *  Pl(mu3, 2)
    end
    pc_sum12_inter = Spline1D(mu, pc_sum12, k=3)
    pc_sum23_inter = Spline1D(mu, pc_sum23, k=3)
    pc_sum31_inter = Spline1D(mu, pc_sum31, k=3)

    for k in 1:_numell
        _ell = k-1
        f = x -> Slepian_multipole_ell(x, _r12[i], _r23[i], _ell, pc_sum12_inter, pc_sum23_inter, pc_sum31_inter)
        I,e = hquadrature(f, llim_mu::Real, ulim_mu::Real, reltol=1e-10, maxevals=10^7)#, abstol=1e-8)
        _idx = (i-1)*_numell + k
        slepian_multi[_idx] = I * (2*_ell+1)/2 #* (4*pi)
        err_slepian_multi[_idx] = e * (2*_ell+1)/2 #* (4*pi)
    end
end

open("../results/Euclid/Slepian_julia_flagship_damped_higher_qmax10_ell20.txt", "w") do io
    writedlm(io, [slepian_multi])
end

slepian_corr = zeros((lenr12_))
for i in 1:lenr12_
    mu11 = (_r12[i]^2+_r23[i]^2-_r31[i]^2)/(2*_r12[i]*_r23[i])
    for j in 1:_numell
        slepian_corr[i]  += (slepian_multi[((i-1)*_numell)+j] * Pl(mu11, j-1))
    end
end


open("../results/Euclid/Slepian_corr_julia_flagship_damped_higher_qmax10.txt", "w") do io
    writedlm(io, [_r12, _r23, _r31, slepian_corr])
end
