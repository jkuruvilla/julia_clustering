import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

plt.rcParams['xtick.labelsize']=16
plt.rcParams['ytick.labelsize']=16

_k1,bk = np.loadtxt("../Anna_Checks/input/Bispectrum_kmin5e-05_kmax100.0_N400_Nmu5000_mu-1_k20.07197604397698391.txt", skiprows=2, unpack=True)

_k1 = _k1[80:341]

b0 = np.loadtxt("B_multipole_0.txt")
b1 = np.loadtxt("B_multipole_1.txt")
b2 = np.loadtxt("B_multipole_2.txt")
b3 = np.loadtxt("B_multipole_3.txt")
b4 = np.loadtxt("B_multipole_4.txt")
b5 = np.loadtxt("B_multipole_5.txt")

plt.plot(_k1, np.abs(b0), "C0-", label="$B_0$")
plt.plot(_k1, np.abs(b1), "C1-", label="$B_1$")
plt.plot(_k1, np.abs(b2), "C2-", label="$B_2$")
plt.plot(_k1, np.abs(b3), "C3-", label="$B_3$")
plt.plot(_k1, np.abs(b4), "C4-", label="$B_4$")
plt.plot(_k1, np.abs(b5), "C5-", label="$B_5$")
plt.legend(fontsize=14, loc=1)
plt.yscale('log')
plt.xscale('log')
plt.title("$k_2 = 0.07197$", fontsize=16)
plt.ylabel("$\\left|B_{\\ell}(k_1, k_2)\\right|$", fontsize=16)
plt.xlabel("$k_1\ (h\,\\mathrm{Mpc}^{-1})$", fontsize=16)
plt.tick_params(top=True, right=True, length=6)
plt.savefig("../Anna_Checks/plots/bispectrum_multipoles_check_2.png", format='png', dpi=300, bbox_inches='tight')
plt.show()
