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
    _va = k^3 * ( (sin(arg_)/arg_^2) - (cos(arg_)/arg_) )
    integr = pk_ev * _va
    return integr
end

function xi_oneminus(kk::Float64, rr::Float64, pk_int)
    arg_ = kk*rr
    pk_ev = evaluate(pk_int, kk)
    _va = k * ( (sin(arg_)/arg_^2) - (cos(arg_)/arg_) )
    integr = pk_ev * _va
    return integr
end

function xi_two(kk::Float64, rr::Float64, pk_int)
    arg_ = kk*rr
    pk_ev = evaluate(pk_int, kk)
    _va = k^2 * ( ((3/arg_^2 - 1) *(sin(arg_)/arg_)) - (3*cos(arg_)/arg_^2) )
    integr = pk_ev * _va
    return integr
end

function Slepian_precyclic(r1::Float64, r2::Float64, _xi, _xim, _xip, _xitwo)
    pc0 = 34/21 * _xi(r1) * _xi(r2)
    pc1 = -1 * ((_xim(r1)*_xip(r2)) + (_xim(r2)*_xip(r1)))
    pc2 = 8/21 * (_xitwo(r1) * _xitwo(r2))
    _sum = pc0 + pc1 + pc2
    return _sum
end

function Slepian_multipole_ell(_mu::Float64, _ell::Real, _prec_sum)
    integr = _prec_sum(_mu) * Pl(_mu, _ell)
    return integr
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

rrr = collect(0.2:0.05:160)
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
    I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-10, maxevals=10^7)#, abstol=1e-8)
    corr[i] = I / (2*pi^2)
    err_corr[i] = e / (2*pi^2)

    f = x -> xi_oneplus(x, rrr[i], pk_inter)
    I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-10, maxevals=10^7)#, abstol=1e-8)
    corr_oneplus[i] = I / (2*pi^2) #* (4*pi)
    err_corr_oneplus[i] = e / (2*pi^2) #* (4*pi)

    f = x -> xi_oneminus(x, rrr[i], pk_inter)
    I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-10, maxevals=10^7)#, abstol=1e-8)
    corr_oneminus[i] = I / (2*pi^2) #* (4*pi)
    err_corr_oneminus[i] = e / (2*pi^2) #* (4*pi)

    f = x -> xi_two(x, rrr[i], pk_inter)
    I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-10, maxevals=10^7)#, abstol=1e-8)
    corr_two[i] = I / (2*pi^2) #* (4*pi)
    err_corr_two[i] = e / (2*pi^2) #* (4*pi)

    _rrr[i] = rrr[i]
    println(i)
end

xi_inter = Spline1D(_rrr, corr, k=3)
xi_onep_inter = Spline1D(_rrr, corr_oneplus, k=3)
xi_onem_inter = Spline1D(_rrr, corr_oneminus, k=3)
xi_two_inter = Spline1D(_rrr, corr_two, k=3)

mu = collect(-1.0:0.05:1)

open("../results/Euclid/linear_correlation_julia_flagship_damped_bg.txt", "w") do io
    writedlm(io, [_rrr, corr])
end

open("../results/Euclid/xi_prime_julia_flagship_damped_bg.txt", "w") do io
    writedlm(io, [_rrr, xi_prime])
end


open("../results/Euclid/phi_julia_flagship_damped_higher.txt", "w") do io
    writedlm(io, [_rrr, phi])
end

open("../results/Euclid/phi_prime_julia_flagship_damped_higher.txt", "w") do io
    writedlm(io, [_rrr, phi_prime, phi_der_inter])
end


tri_rr = readdlm("../input/triangles_binSize_5.dat", skipstart=1)
_r12 = tri_rr[:,1]
_r23 = tri_rr[:,2]
_r31 = tri_rr[:,3]

len_r12 = length(_r12)
bg_ = zeros((len_r12))
for i in 1:len_r12
    _precyclic  = BarrigaGatzanaga_connected(_r12[i], _r23[i], (_r12[i]^2+_r23[i]^2-_r31[i]^2)/(2*_r12[i]*_r23[i]), xi_inter, phi_inter)
    _cyclic_one = BarrigaGatzanaga_connected(_r23[i], _r31[i], (_r23[i]^2+_r31[i]^2-_r12[i]^2)/(2*_r23[i]*_r31[i]), xi_inter, phi_inter)
    _cyclic_two = BarrigaGatzanaga_connected(_r31[i], _r12[i], (_r31[i]^2+_r12[i]^2-_r23[i]^2)/(2*_r31[i]*_r12[i]), xi_inter, phi_inter)

    bg_[i] = _precyclic + _cyclic_one + _cyclic_two
end

open("../results/Euclid/bg_julia_flagship_damped_higher.txt", "w") do io
    writedlm(io, [_r12, _r23, _r31, bg_])
end
