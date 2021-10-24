import Dierckx.Spline1D, Dierckx.evaluate
import Cubature.hquadrature, Cubature.hcubature
import Plots.plot, Plots.savefig, Plots.plot!
import DelimitedFiles.readdlm
import Printf.@printf

using CPUTime

π = 3.141592653589793

function correlation_integrand(kk::Float64, rr::Float64, pk_int)
    pk_ev = evaluate(pk_int, kk)
    _va = kk * sin(kk*rr)/(rr)
    integr = pk_ev * _va
    return integr
end

function Cone(rr::Float64, kk::Float64, pk_int)
    arg1 = kk*rr
    ev = evaluate(pk_int, arg1)
    rr2 = rr^2
    factor = 12/rr2 + 10 + 100*rr2 - 42*rr2^2 + 3/rr^3 * (rr2 - 1)^3 * (7*rr2 + 2) * log(abs((1+rr)/(1-rr)))
    integr = ev * factor
    return integr
end

function Ctwo(var, kk::Float64, pk_int)
    arg1 = kk*var[1]
    #print(arg1)
    ev1 = evaluate(pk_int, arg1)
    #println((1+var[1]^2 - 2*var[1]*var[2]))
    #print("\r",(1+var[1]^2 - 2*var[1]*var[2])^0.5)
    arg2 = kk * sqrt(1+var[1]^2 - 2*var[1]*var[2])
    ev2 = evaluate(pk_int, arg2)
    factor = (3*var[1] + 7*var[2] - 10*var[1]*var[2]^2)^2 / (1 + var[1]^2 - 2*var[1]*var[2])^2
    #println(factor)
    integr = ev1 * ev2 * factor
    return integr
end


aa = readdlm("Linear_power_spectrum_Planck15_LRT.txt", skipstart=3)
wavenumber = aa[:,1]
pk = aa[:,2]
len_pk = length(pk)
pk_inter = Spline1D(wavenumber, pk, k=3)

knl = hquadrature(x -> evaluate(pk_inter, x), wavenumber[1], wavenumber[len_pk])
#plot(wavenumber, pk)

new_k = wavenumber[220:15:4000]
len_nk = length(new_k)
print(len_nk)


pk_lrt = zeros((len_nk))
err_pk_lrt = zeros((len_nk))

CPUtic()
for i in 1:len_nk
    u_limit = 10 #50/(1+exp(-1 * (new_k[i])))
    pkk = evaluate(pk_inter, new_k[i])
    fone = y -> Cone(y, new_k[i], pk_inter)
    I_one,e_one = hquadrature(fone, 0, u_limit)
    c1 = 1/252 * new_k[i]^3/(4*pi^2) * pkk * I_one
    ftwo = x -> Ctwo(x, new_k[i], pk_inter)
    I,e = hcubature(ftwo, Vector([0,-1]), Vector([u_limit, 1]))
    #println(I)
    c2 = I * new_k[i]^3/(4*pi^2) * 1/98
    damp = exp(-new_k[i]^2/(6*pi^2)*knl[1])
    pk_lrt[i] = (pkk + c1 + c2) * damp
    err_pk_lrt[i] = e
    println(pk_lrt[i])
end
CPUtoc()


#plot(rrr, corr_r2, label="ξ", xlabel="r", ylabel="r^2 ξ")
#savefig("sim_corr.pdf")
