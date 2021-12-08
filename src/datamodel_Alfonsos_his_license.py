#!/opt/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table

class ThreePCF_DataModel:

    def __init__(self):
        self.author_name = None
        self.author_surname = None
        self.author_abbr = None
        self.method_3PCF = None
        self.LMAX = None

        self.r12 = None
        self.r13 = None
        self.zeta = None
        self.variable_list = ["author_name", "author_surname", "author_abbr", "method_3PCF", "r12", "r13", "zeta"]

    def check_filled (self):

        error = False

        for vl in self.variable_list:
            if getattr(self, vl) is None:
                error = True
                print("%s is not set"%vl)

        if error:
            raise Exception("Variables not set")

    def write (self, table, output_dir, output_file, overwrite=True):

        table.header.set('CHALL', self.author_name+" "+self.author_surname)
        table.header.set('ABBR', self.author_abbr)
        table.header.set('3PCFMETH', self.method_3PCF)
        table.header.set('LMAX', self.LMAX)
        table.writeto(output_dir+output_file, overwrite=overwrite)

        return table

class ThreePCF_Triangles(ThreePCF_DataModel):

    def __init__(self):
        super().__init__()
        self.r23 = None
        self.variable_list += ["r23", ]
        print(self.variable_list)

    def write (self, output_dir, output_file, overwrite=True):

        self.check_filled()

        cr12 = fits.Column(name='r12', array=np.array(self.r12), format='D')
        cr13 = fits.Column(name='r13', array=np.array(self.r13), format='D')
        cr23 = fits.Column(name='r23', array=np.array(self.r23), format='D')
        czeta = fits.Column(name='zeta', array=np.array(self.zeta), format='D')

        table = fits.BinTableHDU.from_columns([cr12, cr13, cr23, czeta])
        table = super().write(table, output_dir, output_file, overwrite)

class ThreePCF_Multipoles(ThreePCF_DataModel):

    def __init__(self):
        super().__init__()
        self.ell = None
        self.variable_list += ["ell" ,]

    def write (self, output_dir, output_file, overwrite=True):

        self.check_filled()

        cr12 = fits.Column(name='r12', array=np.array(self.r12), format='D')
        cr13 = fits.Column(name='r13', array=np.array(self.r13), format='D')
        cell = fits.Column(name='ell', array=np.array(self.ell), format='K')
        czeta = fits.Column(name='zeta', array=np.array(self.zeta), format='D')

        table = fits.BinTableHDU.from_columns([cr12, cr13, cell, czeta])
        table = super().write(table, output_dir, output_file, overwrite)

def main ():

    zeta_j = np.loadtxt("../results/Euclid/BarrigaGaztanaga_Joseph_binsize5_multipoles_new.txt")
    r1, r2, zeta_a = np.loadtxt("../input/Barriga_Gaztanaga_l0_Anna_binSize5.txt", skiprows=1, unpack=True)

    _ellmax = 10
    new_r12 = np.zeros(len(zeta_j))
    new_r13 = np.zeros(len(zeta_j))
    new_ell = np.zeros(len(zeta_j))

    for i in range(len(r1)):
        for j in range(_ellmax):
            _idx = i+j
            new_r12[_idx] = r1[i]
            new_r13[_idx] = r2[i]
            new_ell[_idx] = j

    threepcf = ThreePCF_Multipoles()
    threepcf.author_name = "Joseph"
    threepcf.author_surname = "Kuruvilla"
    threepcf.author_abbr = "JosephK"
    threepcf.method_3PCF = "BarrigaGaztanaga"
    threepcf.r12 = new_r12
    threepcf.r13 = new_r13
    threepcf.ell = new_ell
    threepcf.zeta = zeta_j

    threepcf.write("./", "BarrigaGaztanaga_multipoles_lmax10.fits")

if __name__ == "__main__":
    main()
