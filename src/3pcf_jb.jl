import Dierckx.Spline1D, Dierckx.evaluate
import Cubature.hquadrature
import DelimitedFiles.readdlm, DelimitedFiles.writedlm

π = 3.141592653589793

function correlation_integrand(kk::Float64, rr::Float64, pk_int)
    pk_ev = evaluate(pk_int, kk)
    _va = kk * sin(kk*rr)/(rr)
    integr = pk_ev * _va
    return integr
end

function eta_ell_integrand(kk::Float64, rr::Float64, ell::Real, pk_int)
    arg_ = kk*rr
    pk_ev = evaluate(pk_int, kk)
    k_ell = (1/(kk^ell)) * kk^2
    _va =  (arg_*cos(arg_) - sin(arg_)) / (kk*rr^3)
    integr = pk_ev * k_ell * _va
    return integr
end

function epsilon_ell_integrand(kk::Float64, rr::Float64, ell::Real, pk_int)
    arg_ = kk*rr
    pk_ev = evaluate(pk_int, kk)
    k_ell = (1/(kk^ell)) * kk^2
    _va =  (3*(sin(arg_) - arg_*cos(arg_)) - arg_^2*sin(arg_)) / (kk*rr^5)
    integr = pk_ev * k_ell * _va
    return integr
end

function jingborner_connected(r1::Float64, r2::Float64, coschi::Float64, xi, eta_zero_interpolated, eta_two_interpolated, epsilon_two_interpolated)
    first_term  = 10.0/7.0 * xi(r1) * xi(r2)
    second_term = - (eta_two_interpolated(r1) * eta_zero_interpolated(r2) * coschi *r1 * r2)
    third_term  = - (eta_zero_interpolated(r1) * eta_two_interpolated(r2) * coschi *r1 * r2)
    fourth_term = 4.0/7.0 * epsilon_two_interpolated(r1) * epsilon_two_interpolated(r2) * coschi^2 * r1^2 * r2^2
    fifth_term  = 4.0/7.0 * epsilon_two_interpolated(r1) * eta_two_interpolated(r2) * r1^2
    sixth_term  = 4.0/7.0 * eta_two_interpolated(r1) * epsilon_two_interpolated(r2) * r2^2
    seven_term  = 4.0/7.0 * 3.0 * eta_two_interpolated(r1) * eta_two_interpolated(r2)
    _sum = first_term + second_term + third_term + fourth_term + fifth_term + sixth_term + seven_term
    return _sum
end

#aa = readdlm("Linear_power_spectrum_sim.txt", skipstart=3)
aa = readdlm("../input/flagship_linear_cb_hr_matterpower_z0p0.dat", skipstart=1)
wavenumber = aa[:,1]
pk = aa[:,2]

#plot(wavenumber, pk)

len_pk = length(pk)
damped_pk = zeros((len_pk))
for i in 1:len_pk
    damped_pk[i] = pk[i] * exp(-(wavenumber[i]/10)^2)
end

open("../results/Euclid/damped_pk_j.txt", "w") do io
    writedlm(io, [wavenumber, damped_pk])
end


pk_inter = Spline1D(wavenumber, damped_pk, k=3)

#rrr = collect(0.1:0.1:160)
rrr = collect(0.5:1.0:160)
# rrr = readdlm("rbins_anna.txt")
# rrr = readdlm("rbins_fftlog_anna.txt")
len_r00 = length(rrr) # 300
println(len_r00)

llim = wavenumber[1]
ulim = wavenumber[len_pk] #100

corr = zeros((len_r00))
eta_zero = zeros((len_r00))
eta_two = zeros((len_r00))
# epsilon_zero = zeros((len_r00))
epsilon_two = zeros((len_r00))

err_corr = zeros((len_r00))
err_eta_zero = zeros((len_r00))
err_eta_two = zeros((len_r00))
# err_epsilon_zero = zeros((len_r00))
err_epsilon_two = zeros((len_r00))

_rrr = zeros((len_r00))

rr = 0.0

for i in 1:len_r00
    f = x -> correlation_integrand(x, rrr[i], pk_inter)
    I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-10, maxevals=10^7)#, abstol=1e-8)
    corr[i] = I / (2*pi^2)
    err_corr[i] = e / (2*pi^2)

    f = x -> eta_ell_integrand(x, rrr[i], 0, pk_inter)
    I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-9, maxevals=10^7)#, abstol=1e-8)
    eta_zero[i] = I / (2*pi^2)
    err_eta_zero[i] = e / (2*pi^2)
    println(i)

    f = x -> eta_ell_integrand(x, rrr[i], 2, pk_inter)
    I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-9, maxevals=10^7)#, abstol=1e-8)
    eta_two[i] = I / (2*pi^2)
    err_eta_two[i] = e / (2*pi^2)

    # f = x -> epsilon_ell_integrand(x, rrr[i], 0, pk_inter)
    # I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-8, maxevals=10^7)#, abstol=1e-8)
    # epsilon_zero[i] = I / (2*pi^2)
    # err_epsilon_zero[i] = e / (2*pi^2)

    f = x -> epsilon_ell_integrand(x, rrr[i], 2, pk_inter)
    I,e = hquadrature(f, llim::Real, ulim::Real,reltol=1e-9, maxevals=10^7)#, abstol=1e-8)
    epsilon_two[i] = I / (2*pi^2)
    err_epsilon_two[i] = e / (2*pi^2)

    _rrr[i] = rrr[i]

end

open("../results/Euclid/linear_correlation_julia_flagship_damped_higher_qmax10.txt", "w") do io
    writedlm(io, [_rrr, corr])
end

open("../results/Euclid/eta_zero_julia_flagship_damped_higher_qmax10.txt", "w") do io
    writedlm(io, [_rrr, eta_zero])
end
open("../results/Euclid/eta_two_julia_flagship_damped_higher_qmax10.txt", "w") do io
    writedlm(io, [_rrr, eta_two])
end

# open("epsilon_zero_julia.txt", "w") do io
#     writedlm(io, [rrr, epsilon_zero])
# end
open("../results/Euclid/epsilon_two_julia_flagship_damped_higher_qmax10.txt", "w") do io
    writedlm(io, [_rrr, epsilon_two])
end

xi_inter = Spline1D(_rrr, corr, k=3)
eta_zero_inter = Spline1D(_rrr, eta_zero, k=3)
eta_two_inter = Spline1D(_rrr, eta_two, k=3)
epsilon_two_inter = Spline1D(_rrr, epsilon_two, k=3)

tri_rr = readdlm("../input/triangles_binSize_5.dat", skipstart=1)
_r12 = tri_rr[:,1]
_r23 = tri_rr[:,2]
_r31 = tri_rr[:,3]

len_r12 = length(_r12)
jb_ = zeros((len_r12))
for i in 1:len_r12
    _precyclic  = jingborner_connected(_r12[i], _r23[i], (_r12[i]^2+_r23[i]^2-_r31[i]^2)/(2*_r12[i]*_r23[i]), xi_inter, eta_zero_inter, eta_two_inter, epsilon_two_inter)
    _cyclic_one = jingborner_connected(_r23[i], _r31[i], (_r23[i]^2+_r31[i]^2-_r12[i]^2)/(2*_r23[i]*_r31[i]), xi_inter, eta_zero_inter, eta_two_inter, epsilon_two_inter)
    _cyclic_two = jingborner_connected(_r31[i], _r12[i], (_r31[i]^2+_r12[i]^2-_r23[i]^2)/(2*_r31[i]*_r12[i]), xi_inter, eta_zero_inter, eta_two_inter, epsilon_two_inter)

    jb_[i] = _precyclic + _cyclic_one + _cyclic_two
end

open("../results/Euclid/JingBorner_Joseph_binsize5.txt", "w") do io
    writedlm(io, [_r12, _r23, _r31, jb_])
end
