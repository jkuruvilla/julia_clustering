import Printf.@printf
import Dierckx.Spline1D, Dierckx.evaluate
import Cubature.hquadrature
import Plots.plot, Plots.savefig, Plots.plot!

using HDF5
using CPUTime

π = 3.141592653589793

function correlation_integrand(kk::Float64)
    pk_ev = evaluate(pk_int, kk)
    _va = kk * sin(kk*rr)/(rr)
    integr = pk_ev * _va
    return integr
end

fid = h5open("/mnt/Edu/Quijote_Pk_env/Particles/Fiducial/0_real.h5", "r")
pk = fid["lambdaRef/AutoSpectra/Monopole/AA"][:]
close(fid)

header = h5readattr("/mnt/Edu/Quijote_Pk_env/Particles/Fiducial/0_real.h5", "Header")
wavenumber = header["ks"]

plot(wavenumber, pk)

len_pk = length(pk)

pk_int = Spline1D(wavenumber, pk, k=3)

rrr = collect(0.1:0.5:140.1)
len_r00 = length(rrr)

corr = zeros((len_r00))

rr = 0.0
CPUtic()
for i in 1:len_r00
    #@printf("%i\n", i)
    #global rr = rrr[i]
    global rr = rrr[i]
    I,e = hquadrature(correlation_integrand, wavenumber[1]::Real, wavenumber[len_pk]::Real)
    corr[i] = I / (2*pi^2)
end
CPUtoc()

corr_r2 = zeros((len_r00))
for i in 1:len_r00
   corr_r2[i] = rrr[i]^2 * corr[i]
end

plot(rrr, corr_r2, label="ξ", xlabel="r", ylabel="r^2 ξ")
savefig("sim_corr.pdf")
