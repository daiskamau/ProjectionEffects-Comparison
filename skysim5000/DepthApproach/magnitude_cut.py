#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py

def mag_i_lim_Rykoff14(z): # Rykoff14, Eq9, mstar+1.75
    if z <= 0.5:
        mag_star = 22.44 + 3.36*np.log(z) + 0.273 * np.log(z)**2 - 0.0618 * np.log(z)**3 - 0.0227 * np.log(z)**4
    else:
        mag_star = 22.94 + 3.08*np.log(z) - 11.22 * np.log(z)**2 -  27.11 * np.log(z)**3 -  18.02 * np.log(z)**4
    return mag_star + 1.75


# z, i = np.loadtxt(file_path + 'skysim-members/iband_max.dat',skiprows=1, unpack=True)
# i_vs_redshift = interp1d(z, i)


if __name__ == "__main__":
    z_list = np.linspace(0.1,1)
    plt.plot(z_list, mag_i_lim_Rykoff14(z_list))

    plt.xlabel('z')
    plt.ylabel('mag i cut')
    plt.savefig('/magnitude_cut.png')
    plt.show()
