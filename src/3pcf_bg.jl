import Dierckx.Spline1D, Dierckx.evaluate, Dierckx.derivative
import Cubature.hquadrature
import DelimitedFiles.readdlm, DelimitedFiles.writedlm

π = 3.141592653589793

function correlation_integrand(kk::Float64, rr::Float64, pk_int)
    pk_ev = evaluate(pk_int, kk)
    _va = kk * sin(kk*rr)/(rr)
    integr = pk_ev * _va
    return integr
end

function phi_integrand(kk::Float64, rr::Float64, pk_int)
    arg_ = kk*rr
    pk_ev = evaluate(pk_int, kk)
    _va =  ( sin(arg_) ) / arg_
    integr = pk_ev * _va
    return integr
end

function phi_prime_integrand(kk::Float64, rr::Float64, pk_int)
    arg_ = kk*rr
    pk_ev = evaluate(pk_int, kk)
    _va =  ( ( arg_*cos(arg_) )-sin(arg_) ) / arg_^2
    integr = pk_ev * _va * kk
    return integr
end

function BarrigaGatzanaga_connected(r1::Float64, r2::Float64, coschi::Float64, _xi, _phi)
    first_term  = 10.0/7.0 * _xi(r1) * _xi(r2)
    second_term = - (coschi * (derivative(_xi,r1) * derivative(_phi,r2) + derivative(_xi,r2) * derivative(_phi,r1)))
    third_term = - ((4.0/7.0) * 3.0 * derivative(_phi,r1) * derivative(_phi,r2) / (r1*r2))
    fourth_term  = - ((4.0/7.0) * _xi(r1) * derivative(_phi,r2) / r2)
    fifth_term  = - ((4.0/7.0) * _xi(r2) * derivative(_phi,r1) / r1)
    sixth_term  = 4.0/7.0 * coschi^2 * (_xi(r1) + ((3*derivative(_phi,r1))/r1)) * (_xi(r2) + ((3*derivative(_phi,r2))/r2))
    _sum = first_term + second_term + third_term + fourth_term + fifth_term + sixth_term
    return _sum
end

function BarrigaGatzanaga_bu_connected(r1::Float64, r2::Float64, coschi::Float64, _xi, _phi_prime)
    first_term  = 10.0/7.0 * _xi(r1) * _xi(r2)
    second_term = - (coschi * (derivative(_xi,r1) * _phi_prime(r2) + derivative(_xi,r2) * _phi_prime(r1)))
    third_term = - ((4.0/7.0) * 3.0 * _phi_prime(r1) * _phi_prime(r2) / (r1*r2))
    fourth_term  = - ((4.0/7.0) * _xi(r1) * _phi_prime(r2) / r2)
    fifth_term  = - ((4.0/7.0) * _xi(r2) * _phi_prime(r1) / r1)
    sixth_term  = 4.0/7.0 * coschi^2 * (_xi(r1) + ((3*_phi_prime(r1))/r1)) * (_xi(r2) + ((3*_phi_prime(r2))/r2))
    _sum = first_term + second_term + third_term + fourth_term + fifth_term + sixth_term
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


pk_inter = Spline1D(wavenumber, damped_pk, k=3)

rrr = collect(0.5:0.1:160)
# rrr = collect(0.2:0.1:160)
# rrr = readdlm("rbins_anna.txt")
# rrr = readdlm("rbins_fftlog_anna.txt")
len_r00 = length(rrr) # 300
println(len_r00)

llim = wavenumber[1]
ulim = wavenumber[len_pk] #100

corr = zeros((len_r00))
phi = zeros((len_r00))
# phi_prime = zeros((len_r00))

err_corr = zeros((len_r00))
err_phi = zeros((len_r00))
# err_phi_prime = zeros((len_r00))

_rrr = zeros((len_r00))

rr = 0.0

for i in 1:len_r00
    f = x -> correlation_integrand(x, rrr[i], pk_inter)
    I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-10, maxevals=10^7)#, abstol=1e-8)
    corr[i] = I / (2*pi^2)
    err_corr[i] = e / (2*pi^2)

    f = x -> phi_integrand(x, rrr[i], pk_inter)
    I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-10, maxevals=10^7)#, abstol=1e-8)
    phi[i] = I / (2*pi^2) #* (4*pi)
    err_phi[i] = e / (2*pi^2) #* (4*pi)

    # f = x -> phi_prime_integrand(x, rrr[i], pk_inter)
    # I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-11, maxevals=10^7)#, abstol=1e-8)
    # phi_prime[i] = I / (2*pi^2) #* (4*pi)
    # err_phi_prime[i] = e / (2*pi^2) #* (4*pi)

    _rrr[i] = rrr[i]
    println(i)
end

xi_inter = Spline1D(_rrr, corr, k=3)
xi_prime = derivative(xi_inter, _rrr)
phi_inter = Spline1D(_rrr, phi, k=3)
phi_der_inter = derivative(phi_inter, _rrr)
#phi_prime_inter = Spline1D(_rrr, phi_prime, k=3)



open("../results/Euclid/linear_correlation_julia_flagship_damped_bg_qmax10.txt", "w") do io
    writedlm(io, [_rrr, corr])
end

open("../results/Euclid/xi_prime_julia_flagship_damped_bg_qmax10.txt", "w") do io
    writedlm(io, [_rrr, xi_prime])
end


open("../results/Euclid/phi_julia_flagship_damped_higher_qmax10.txt", "w") do io
    writedlm(io, [_rrr, phi])
end

# open("../results/Euclid/phi_prime_julia_flagship_damped_higher.txt", "w") do io
#     writedlm(io, [_rrr, phi_prime, phi_der_inter])
# end

open("../results/Euclid/phi_prime_julia_flagship_damped_higher_qmax10.txt", "w") do io
    writedlm(io, [_rrr, phi_der_inter])
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

open("../results/Euclid/bg_julia_flagship_damped_higher_qmax10.txt", "w") do io
    writedlm(io, [_r12, _r23, _r31, bg_])
end
