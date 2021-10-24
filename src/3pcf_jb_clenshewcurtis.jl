import Dierckx.Spline1D, Dierckx.evaluate
import Cubature.pquadrature
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

aa = readdlm("Linear_power_spectrum_sim.txt", skipstart=3)
wavenumber = aa[:,1]
pk = aa[:,2]

#plot(wavenumber, pk)

len_pk = length(pk)

pk_inter = Spline1D(wavenumber, pk, k=3)

rrr = collect(0.1:0.5:200.1)
len_r00 = length(rrr)
println(len_r00)

llim = wavenumber[1]
ulim = wavenumber[len_pk]

corr = zeros((len_r00))
eta_zero = zeros((len_r00))
eta_two = zeros((len_r00))
epsilon_zero = zeros((len_r00))
epsilon_two = zeros((len_r00))

err_corr = zeros((len_r00))
err_eta_zero = zeros((len_r00))
err_eta_two = zeros((len_r00))
err_epsilon_zero = zeros((len_r00))
err_epsilon_two = zeros((len_r00))

rr = 0.0

for i in 1:len_r00
    println(i)
    # f = x -> correlation_integrand(x, rrr[i], pk_inter)
    # I,e = pquadrature(f, llim::Real, ulim::Real, maxevals=10^6)#, abstol=1e-8)
    # corr[i] = I / (2*pi^2)
    # err_corr[i] = e / (2*pi^2)

    f = x -> eta_ell_integrand(x, rrr[i], 0, pk_inter)
    I,e = pquadrature(f, llim::Real, ulim::Real, reltol=1e-8, maxevals=10^8)#, abstol=1e-8)
    eta_zero[i] = I / (2*pi^2)
    err_eta_zero[i] = e / (2*pi^2)

    # f = x -> eta_ell_integrand(x, rrr[i], 2, pk_inter)
    # I,e = pquadrature(f, llim::Real, ulim::Real, maxevals=10^6)#, abstol=1e-8)
    # eta_two[i] = I / (2*pi^2)
    # err_eta_two[i] = e / (2*pi^2)

    # f = x -> epsilon_ell_integrand(x, rrr[i], 0, pk_inter)
    # I,e = pquadrature(f, llim::Real, ulim::Real, maxevals=10^6)#, abstol=1e-8)
    # epsilon_zero[i] = I / (2*pi^2)
    # err_epsilon_zero[i] = e / (2*pi^2)

    # f = x -> epsilon_ell_integrand(x, rrr[i], 2, pk_inter)
    # I,e = pquadrature(f, llim::Real, ulim::Real, maxevals=10^6)#, abstol=1e-8)
    # epsilon_two[i] = I / (2*pi^2)
    # err_epsilon_two[i] = e / (2*pi^2)

end

# open("linear_correlation_julia_clenshaw.txt", "w") do io
#     writedlm(io, [rrr, corr])
# end

open("eta_zero_julia_clenshaw_reltolem8.txt", "w") do io
    writedlm(io, [rrr, eta_zero])
end
# open("eta_two_julia_clenshaw.txt", "w") do io
#     writedlm(io, [rrr, eta_two])
# end

# open("epsilon_zero_julia_clenshaw.txt", "w") do io
#     writedlm(io, [rrr, epsilon_zero])
# end
# open("epsilon_two_julia_clenshaw.txt", "w") do io
#     writedlm(io, [rrr, epsilon_two])
# end
