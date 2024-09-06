#!/usr/bin/env python
import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
import time 
# import fitsio
from scipy import spatial
from multiprocessing import Pool, cpu_count
import argparse

from astropy.table import Table
import astropy.coordinates
from astropy.cosmology import FlatLambdaCDM
astropy_cosmo = FlatLambdaCDM(H0=70, Om0=0.286,Ob0 = 0.047)
h = 0.7

from magnitude_cut import mag_i_lim_Rykoff14
from member_color_interp import *
from Extraction import galaxies

filepath = r"/bsuscratch/shuleicao/Cardinalv3/"

file_path = r"/bsuhome/gladyskamau/BSU-Research/Cardinal/Data/"

# file_path1 = r"/bsuhome/gladyskamau/BSU-Research/Cardinal/Data/DepthData/"

parser = argparse.ArgumentParser()
parser.add_argument("output_file_path", type=str, help="Path to the output file directory")
parser.add_argument("depth", type=int, help="Depth value")

args = parser.parse_args()

output_file_path = args.output_file_path

depth = args.depth
chisq_cut = 9 #1e10 # TODO: no color cut for now  # 9

## Read in All the  halos
sky_02_1_all = Table(np.load( file_path + "halos_cardinal_w_lensing_1687.npy"))

### Check whether cylinder richness has already been computed to avoid repetition
# lambda_60 = Table.read(file_path + 'RichnessDepth-Data/RichnessDepth_1_035_05.dat', format='ascii')
# pixelss = list(np.unique(sky_02_1_all['pixel_id'][~np.isin(sky_02_1_all['pixel_id'],np.unique(lambda_60['pixelid']))]))
# sky_02_1_all = sky_02_1_all[np.isin(sky_02_1_all['pixel_id'],pixelss)]

# sky_02_1_all = sky_02_1_all[np.isin(sky_02_1_all['pixel_id'],np.unique(sky_02_1_all['pixel_id'])[:5])]

## Read all the galaxies (gold)
### getting halos
fname = filepath + 'Cardinal-3_v2.0_Y6a_gold.h5'
f = h5py.File(fname,'r')
catalog = f['catalog']
gold = catalog['gold']

mag_g = gold['mag_g'][:]
mag_r = gold['mag_r'][:]
mag_i = gold['mag_i'][:]
mag_z = gold['mag_z'][:]

mag_g_mask = mag_g < 24
mag_r_mask = mag_r < 24
mag_i_mask = mag_i < 24
mag_z_mask = mag_z < 24

# Combine the masks using logical AND to select values that satisfy all conditions
combined_mask =  mag_g_mask & mag_r_mask &mag_i_mask &mag_z_mask 

# Use the combined mask to filter the arrays
filtered_mag_g = mag_g[combined_mask]
filtered_mag_r = mag_r[combined_mask]
filtered_mag_i = mag_i[combined_mask]
filtered_mag_z = mag_z[combined_mask]


#### getting redshift of the galaxies
fname_z = filepath + 'Cardinal-3_v2.0_Y6a_bpz.h5'
f_z = h5py.File(fname_z,'r')
catalog_z = f_z['catalog']
bpz = catalog_z['bpz']
z_gal = bpz['redshift_cos'][:]
chi_gal = astropy_cosmo.comoving_distance(z_gal).value*h

ra_g = gold['ra'][:][combined_mask]
ra_filtered = ra_g - 360 * (ra_g>180)
dec_filtered = gold['dec'][:][combined_mask]
z_filtered = z_gal[combined_mask]
chi_filtered = chi_gal[combined_mask]

############################################################################################
class AlternativeRichness(object):
    def __init__(self, zmin, zmax, pixel_id):#, sample, method, radius, depth, chisq_cut=chisq_cut):
        self.zmin = zmin
        self.zmax = zmax
        self.pixel_id = pixel_id
        self.zmid = 0.5 * (self.zmin + self.zmax)
        self.save_name = 'z_%g_%g_%i'%(self.zmin, self.zmax,self.pixel_id)


        self.richness_file = output_file_path+'ngal_%s.dat'%(self.save_name)

        self.halos = galaxies(pixel_id,sky_02_1_all)
        
    def get_galaxies(self):
        ra_all = ra_filtered  
        dec_all = dec_filtered 
        z_all = z_filtered 
        chi_all = chi_filtered 
        mag_g_all = filtered_mag_g 
        mag_r_all = filtered_mag_r 
        mag_i_all = filtered_mag_i 
        mag_z_all = filtered_mag_z 
        print('len_ra_gall_all',len(ra_all))

      #### Step 2: cut redshift and magnitude
        sel = (z_all >= self.zmin-0.1)&(z_all <= self.zmax+0.1)  & (mag_i_all > 10) & (mag_i_all < i_vs_redshift(self.zmid)) 

        # print('done selection')
        self.chi_gal = chi_all[sel]
        self.ra_gal = ra_all[sel] # deg
        self.dec_gal = dec_all[sel] # deg

        #### Step 3: cut color chisq based on the color templates 
        g_r_mean = g_r_vs_redshift(self.zmid)
        g_r_std = sigma_g_r_vs_redshift(self.zmid)

        r_i_mean = r_i_vs_redshift(self.zmid)
        r_i_std = sigma_r_i_vs_redshift(self.zmid)

        i_z_mean = i_z_vs_redshift(self.zmid)
        i_z_std = sigma_i_z_vs_redshift(self.zmid)

        g_r = mag_g_all[sel] - mag_r_all[sel]
        r_i = mag_r_all[sel] - mag_i_all[sel]
        i_z = mag_i_all[sel] - mag_z_all[sel]
        # print(len(g_r))

        chisq = (g_r - g_r_mean)**2 / g_r_std**2
        chisq += (r_i - r_i_mean)**2 / r_i_std**2
        chisq += (i_z - i_z_mean)**2 / i_z_std**2

        sel2 = (chisq < chisq_cut)

        self.chi_gal = self.chi_gal[sel2]
        self.ra_gal = self.ra_gal[sel2] * np.pi / 180.
        self.dec_gal = self.dec_gal[sel2] * np.pi / 180.

        self.gal_taken = np.zeros(len(self.chi_gal)) # for percolation


    def get_halos(self):
        hid, redshift_halo, mass, ra_halo, dec_halo, pixelid, chi_halo, DA_halo = self.halos

         ## Step 2: filter redshift 
        sel = (redshift_halo > self.zmin)&(redshift_halo < self.zmax) &(mass > 0)
        hid = np.array(hid[sel])
        ra_halo = np.array(ra_halo[sel])
        dec_halo = np.array(dec_halo[sel])
        redshift_halo = np.array(redshift_halo[sel])  
        mass = np.array(mass[sel])
        pixelid = np.array(pixelid[sel])
        chi_halo = np.array(chi_halo[sel])
        DA_halo = np.array(DA_halo[sel])
        
        ## Step 3: sort by mass
        sort = np.argsort(-mass)
        self.hid = hid[sort]
        self.ra_halo = ra_halo[sort] * np.pi / 180.
        self.dec_halo = dec_halo[sort] * np.pi / 180.
        self.redshift_halo = redshift_halo[sort]
        self.mass = mass[sort]
        self.pixelid = pixelid[sort]
        self.chi_halo = chi_halo[sort]
        self.DA_halo = DA_halo[sort]
#         print(self.DA_halo,self.redshift_halo)

    def build_trees(self):
        x_gal = (0.5*np.pi - self.dec_gal) # like theta
        y_gal = np.cos(self.dec_gal) * self.ra_gal # like phi sin(theta)
        x_halo = (0.5*np.pi - self.dec_halo)
        y_halo = np.cos(self.dec_halo) * self.ra_halo

        gal_position = np.dstack([x_gal, y_gal])[0]
        gal_tree = spatial.cKDTree(gal_position)

        halo_position = np.dstack([x_halo, y_halo])[0]
        halo_tree = spatial.cKDTree(halo_position)

        DA_min = min(self.DA_halo)
        rmax_tree = 4 / DA_min #
        self.indexes_tree = halo_tree.query_ball_tree(gal_tree, r=rmax_tree)
        self.gal_taken = np.zeros(len(self.chi_gal)) # for percolation


    def get_richness(self, i_halo):
        #### step 0: get the tree
        # print(i_halo)
        gal_ind = self.indexes_tree[i_halo]
        #### step 1: cut the LOS ####
        d_los = self.chi_gal[gal_ind] - self.chi_halo[i_halo]
        sel_z = (np.abs(d_los) < depth)
        sel_z = sel_z & (self.gal_taken[gal_ind] < 1e-4)

        #### steo 2: calculate angular separation & projected radius
        d_ra = self.ra_gal[gal_ind][sel_z] - self.ra_halo[i_halo]
        d_dec = self.dec_gal[gal_ind][sel_z] - self.dec_halo[i_halo]
        ang_sep = d_ra**2 * np.cos(self.dec_halo[i_halo])**2 + d_dec**2
        ang_sep = np.sqrt(ang_sep)
#         DA = self.DA_halo # Mpc/h
        r = self.DA_halo[i_halo] * ang_sep # physical Mpc/h

        #### step 3: iteratively calculating r_lambda
        rlam_ini = 2 #1
        rlam = rlam_ini
        for iteration in range(100):
            sel = (r < rlam)
            ngal = len(r[sel])
            rlam_old = rlam
            rlam = (ngal/100.)**0.2 # r is already in physical Mpc/h
            if abs(rlam - rlam_old) < 1e-5 or rlam < 1e-6:
                break

        #### Step 4: do percolation ####
        if rlam > 0:
            sel_mem = (r < rlam)
            self.gal_taken[np.array(gal_ind)[sel_z][sel_mem]] = 1
        return rlam, ngal


    def measure_richness(self, pixel_id):
        nh = len(self.ra_halo)
        #### save the richness file:  ####
        outfile1 = open(self.richness_file, 'w')
        outfile1.write('#haloid, mass, ra, dec, redshift, rlam, lam, pixelid \n')
        for ih in range(nh):
            rlam, lam = self.get_richness(ih)
            # print(lam>0,rlam, lam)
            if (lam > 0) and (self.pixelid[ih] == pixel_id): #and (self.pixelid == pixel_id)
                outfile1.write('%12i %15e %12g %12g %12g %12g %12g %12i \n'%(self.hid[ih], self.mass[ih], self.ra_halo[ih], 
                                                                        self.dec_halo[ih], self.redshift_halo[ih], rlam,
                                                                             lam, self.pixelid[ih]))
        outfile1.close()

        print('DONE!!')


def run_parallel(pixels):
    for i in range(len(pixels)):
        hpix = pixels[i]
        zmin = 0.35 #0.2 #0.3
        zmax = 0.5 #0.32
        car = AlternativeRichness(zmin=zmin, zmax=zmax, pixel_id = hpix)
        car.get_galaxies()
        car.get_halos()
        car.build_trees()
        car.measure_richness(hpix)
    
     
procs = 3    #cpu_count()
df_chunks = np.array_split(np.unique(sky_02_1_all['pixel_id']) ,procs)
if __name__ == '__main__':
    start_time = time.time() 
    with Pool(processes=procs) as p: 
        p.map(run_parallel, df_chunks)
        end_time = time.time() 
        elapsed_time = (end_time - start_time)/60
        print(f"Elapsed time: {elapsed_time:.2f} minutes")

