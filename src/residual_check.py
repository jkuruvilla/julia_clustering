__license__   = "GNU GPLv3 <https://www.gnu.org/licenses/gpl.txt>"
__copyright__ = "2021, Joseph Kuruvilla"
__author__    = "Joseph Kuruvilla <joseph.kuruvilla@universite-paris-saclay.fr>"
__version__   = "1.0"

'''
Program to plot the residual between JB and BG results

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

# ------------------
# Importing Modules
# ------------------

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = r"\usepackage{bm} \usepackage{amsmath}"
plt.rcParams['xtick.labelsize']=16
plt.rcParams['ytick.labelsize']=16

# --------------------
#   Program  start
# --------------------

if __name__ == "__main__":

    r1, r2, r3, bg = np.loadtxt("../results/Euclid/bg_julia_flagship_damped_higher_qmax50.txt")
    r1, r2, r3, bg10 = np.loadtxt("../results/Euclid/bg_julia_flagship_damped_higher_qmax10.txt")

    r1, r2, r3, jb = np.loadtxt("../results/Euclid/jb_julia_flagship_damped_higher_qmax50.txt")
    r1, r2, r3, jb10 = np.loadtxt("../results/Euclid/jb_julia_flagship_damped_higher_qmax10.txt")

    plt.plot((bg10-jb10)/jb10)
    #plt.legend(fontsize=14)
    plt.title("$q_\\mathrm{max} = 10$", fontsize=16, y=1.02)
    plt.ylabel("$\\frac{\\zeta_{\\mathrm{BG}}-\\zeta_{\\mathrm{JB}}}{\\zeta_{\\mathrm{JB}}}$", fontsize=16)
    plt.xlabel("$\\mathrm{Triangle\,\,Configurations}$", fontsize=16)
    plt.tick_params(top=True, right=True, length=6)
    plt.savefig("../plots/qmax10_jb_bg.png", format='png', dpi=300, bbox_inches='tight')
    plt.show()

    plt.plot((bg-jb)/jb)
    #plt.legend(fontsize=14)
    plt.title("$q_\\mathrm{max} = 50$", fontsize=16, y=1.02)
    plt.ylabel("$\\frac{\\zeta_{\\mathrm{BG}}-\\zeta_{\\mathrm{JB}}}{\\zeta_{\\mathrm{JB}}}$", fontsize=16)
    plt.xlabel("$\\mathrm{Triangle\,\,Configurations}$", fontsize=16)
    plt.tick_params(top=True, right=True, length=6)
    plt.savefig("../plots/qmax50_jb_bg.png", format='png', dpi=300, bbox_inches='tight')
    plt.show()
