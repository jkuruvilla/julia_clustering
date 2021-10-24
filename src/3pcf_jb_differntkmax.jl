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

aa = readdlm("Linear_power_spectrum_sim.txt", skipstart=3)
wavenumber = aa[:,1]
pk = aa[:,2]

#plot(wavenumber, pk)

len_pk = length(pk)

pk_inter = Spline1D(wavenumber, pk, k=3)

rrr = collect(0.1:0.5:200.1)
# rrr = readdlm("rbins_anna.txt")
# rrr = readdlm("rbins_fftlog_anna.txt")
len_r00 = length(rrr)
println(len_r00)

llim = wavenumber[1]

corr = zeros((len_r00))
eta_zero = zeros((len_r00))
eta_two = zeros((len_r00))
epsilon_two = zeros((len_r00))

err_corr = zeros((len_r00))
err_eta_zero = zeros((len_r00))
err_eta_two = zeros((len_r00))
err_epsilon_two = zeros((len_r00))

rr = 0.0

ulim_list = [30, 50, 100, 150, 200, 300]

for j in 1:6
    ulim = ulim_list[j]
    for i in 1:len_r00
        println(i)
        f = x -> correlation_integrand(x, rrr[i], pk_inter)
        I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-8, maxevals=10^7)#, abstol=1e-8)
        corr[i] = I / (2*pi^2)
        err_corr[i] = e / (2*pi^2)

        f = x -> eta_ell_integrand(x, rrr[i], 0, pk_inter)
        I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-8, maxevals=10^7)#, abstol=1e-8)
        eta_zero[i] = I / (2*pi^2)
        err_eta_zero[i] = e / (2*pi^2)

        f = x -> eta_ell_integrand(x, rrr[i], 2, pk_inter)
        I,e = hquadrature(f, llim::Real, ulim::Real, reltol=1e-8, maxevals=10^7)#, abstol=1e-8)
        eta_two[i] = I / (2*pi^2)
        err_eta_two[i] = e / (2*pi^2)

        f = x -> epsilon_ell_integrand(x, rrr[i], 2, pk_inter)
        I,e = hquadrature(f, llim::Real, ulim::Real,reltol=1e-8, maxevals=10^7)#, abstol=1e-8)
        epsilon_two[i] = I / (2*pi^2)
        err_epsilon_two[i] = e / (2*pi^2)

    end

    open("jbelements_r_xi_eta0_eta2_epsilon2_kmax$(ulim_list[j]).txt", "w") do io
        writedlm(io, [rrr, corr, eta_zero, eta_two, epsilon_two])
    end
end
