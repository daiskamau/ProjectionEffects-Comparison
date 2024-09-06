#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

file_path =  r"/global/u1/k/kamau/SE-CLMM-LSSTDESC/project-1/Data/skysim-members/"

pmem_cut = 0.9
z1, g_r, grlow, grhigh = np.loadtxt(file_path + 'g_r.dat',skiprows=1, unpack=True)
g_r_vs_redshift = interp1d(z1, g_r)
sigma_g_r_vs_redshift = interp1d(z1, (grhigh-grlow)/2.)

z3, r_i, rilow, rihigh = np.loadtxt(file_path + 'r_i.dat',skiprows=1, unpack=True)
r_i_vs_redshift = interp1d(z3, r_i)
sigma_r_i_vs_redshift = interp1d(z3, (rihigh-rilow)/2.)


z5, i_z, izlow, izhigh = np.loadtxt(file_path + 'i_z.dat',skiprows=1, unpack=True)
i_z_vs_redshift = interp1d(z5, i_z)
sigma_i_z_vs_redshift = interp1d(z5, (izhigh-izlow)/2.)


z7, z_y, zylow, zyhigh = np.loadtxt(file_path + 'z_y.dat',skiprows=1, unpack=True)
z_y_vs_redshift = interp1d(z7, z_y)
sigma_z_y_vs_redshift = interp1d(z7, (zyhigh-zylow)/2.)

redshift, z = np.loadtxt(file_path + 'zband_max.dat',skiprows=1, unpack=True)
z_vs_redshift = interp1d(redshift, z)



if __name__ == "__main__":
    for ipanel in range(3):
        if ipanel == 0:
            ylabel = 'g-r'
        if ipanel == 1:
            ylabel = 'r-i'
        if ipanel == 2:
            ylabel = 'i-z'
        z, med, low, high = np.loadtxt(file_path + '%s_p.dat'%(ylabel), skiprows=1, unpack=True)
        line = plt.plot(z, med, label=ylabel)
        plt.fill_between(z, low, high, interpolate=True, facecolor=line[0].get_c(), alpha=0.1)
        plt.legend()
        plt.xlabel('z')
        plt.ylabel('color')
        plt.title(r'$\rm p_{mem} > %g$'%(pmem_cut))
    plt.savefig('member_color_vs_z_p%g.png'%(pmem_cut))

    plt.show()
