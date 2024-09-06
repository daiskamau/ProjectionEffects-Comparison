#!/usr/bin/env python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time 
# import fitsio
# import h5py
from scipy import spatial
from multiprocessing import Pool

# import GCRCatalogs as gcr

from astropy.table import Table
import astropy.coordinates
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=71, Om0=0.265)

# skysim = gcr.load_catalog('skysim5000_v1.1.1')
# skycosmology = skysim.cosmology # Cosmological parameters
# cosmo = FlatLambdaCDM(H0=skycosmology.H0.value, Om0= skycosmology.Om0)
# cosmo = FlatLambdaCDM(H0=100, Om0=0.3) # TODO for Gladys, please change it it back to skyssim when you run on nersc

from magnitude_cut import mag_i_lim_Rykoff14
from member_color_interp import *
from Extraction import galaxies
file_path = r"/global/u1/k/kamau/SE-CLMM-LSSTDESC/project-1/"


output_file_path = file_path + 'Notebooks/skysim5000/Depth/calc-richness/Data/one/' # TODO: Gladys, please set it to your output directory
input_file_path = file_path + './'

depth = 
chisq_cut = 20 #1e10 # TODO: no color cut for now

## Read in All the gcr halos
sky_02_1_all = Table.read(file_path + 'Data/halos/skysim_0.2-1-2152757_comv&DA.dat', format='ascii')

# ### Check whether cylinder richness has already been computed to avoid repetition
# lambda_1 = Table.read(file_path + 'Data/halos/RichnessDepth_1.dat', format='ascii')
# pixelss = list(np.unique(sky_02_1_all['pixel_id'][~np.isin(sky_02_1_all['pixel_id'],np.unique(lambda_1['pixelid']))]))
# sky_02_1_all = sky_02_1_all[np.isin(sky_02_1_all['pixel_id'],pixelss)]

############################################################################################
class AlternativeRichness(object):
    def __init__(self, zmin, zmax, pixel_id):#, sample, method, radius, depth, chisq_cut=chisq_cut):
        self.zmin = zmin
        self.zmax = zmax
        self.pixel_id = pixel_id
        self.zmid = 0.5 * (self.zmin + self.zmax)
        self.save_name = 'z_%g_%g_%i'%(self.zmin, self.zmax,self.pixel_id)

        # self.mag_i_cut = mag_i_lim_Rykoff14(zmin)

        self.richness_file = output_file_path+'ngal_%s.dat'%(self.save_name)

        self.galaxies, self.halos = galaxies(pixel_id,sky_02_1_all)
        
    def get_galaxies(self):

        ra_all,dec_all,z_all,mag_g_all,mag_r_all,mag_i_all,mag_z_all,mag_y_all,chi_all = self.galaxies
  

        # print('ngal = ', len(z_all))
        # print(f'galaxies zmin = {min(z_all)}, zmax = {max(z_all)}')
        

      #### Step 2: cut redshift and magnitude
        ## TODO: magnitude cut not working 
        sel = (z_all >= self.zmin-0.1)&(z_all <= self.zmax+0.1)  & (mag_i_all > 10) & (mag_i_all < i_vs_redshift(self.zmid)) # & mag_i_all > 10)&(mag_i_all < self.mag_i_cut) & (i_vs_redshift(self.zmid))

        # print('done selection')
        self.chi_gal = chi_all[sel]
        self.ra_gal = ra_all[sel] # deg
        self.dec_gal = dec_all[sel] # deg

        #### Step 3: cut color chisq based on the color templates 
        ## TODO for Heidi & Gladys: there might be a bug in color templates
        g_r_mean = g_r_vs_redshift(self.zmid)
        g_r_std = sigma_g_r_vs_redshift(self.zmid)

        r_i_mean = r_i_vs_redshift(self.zmid)
        r_i_std = sigma_r_i_vs_redshift(self.zmid)

        i_z_mean = i_z_vs_redshift(self.zmid)
        i_z_std = sigma_i_z_vs_redshift(self.zmid)

        z_y_mean = z_y_vs_redshift(self.zmid)
        z_y_std = sigma_z_y_vs_redshift(self.zmid)

        g_r = mag_g_all[sel] - mag_r_all[sel]
        r_i = mag_r_all[sel] - mag_i_all[sel]
        i_z = mag_i_all[sel] - mag_z_all[sel]
        z_y = mag_z_all[sel] - mag_y_all[sel]

        chisq = (g_r - g_r_mean)**2 / g_r_std**2
        chisq += (r_i - r_i_mean)**2 / r_i_std**2
        chisq += (i_z - i_z_mean)**2 / i_z_std**2
        chisq += (z_y - z_y_mean)**2 / z_y_std**2

        sel2 = (chisq < chisq_cut)

        self.chi_gal = self.chi_gal[sel2]
        self.ra_gal = self.ra_gal[sel2] * np.pi / 180.
        self.dec_gal = self.dec_gal[sel2] * np.pi / 180.

        self.gal_taken = np.zeros(len(self.chi_gal)) # for percolation

        # print('done selecting galaxies, ngals = ', len(self.ra_gal))

    def get_halos(self):
        
        hid,redshift_halo,mass,ra_halo,dec_halo, pixelid = self.halos #galaxies(pixel_id,sky_02_1_all)[1]
        
        # print('ngal = ', len(hid))

        ## Step 2: filter redshift 
        sel = (redshift_halo > self.zmin)&(redshift_halo < self.zmax)&(mass > 0)
        hid = np.array(hid[sel])
        ra_halo = np.array(ra_halo[sel])
        dec_halo = np.array(dec_halo[sel])
        redshift_halo = np.array(redshift_halo[sel])  # astropy needs arrays
        mass = np.array(mass[sel])
        pixelid = np.array(pixelid[sel])

        ## Step 3: sort by mass
        sort = np.argsort(-mass)
        self.hid = hid[sort]
        self.ra_halo = ra_halo[sort] * np.pi / 180.
        self.dec_halo = dec_halo[sort] * np.pi / 180.
        self.redshift_halo = redshift_halo[sort]
        self.mass = mass[sort]
        self.pixelid = pixelid[sort]

        # Step 4: calculate LOS distance
        self.chi_halo = cosmo.comoving_distance(self.redshift_halo).value

        # print('after cutting, halo number = ', len(redshift_halo))
        # print(f'halo min mass = {min(mass):e}, max mass = {max(mass):e}')
        # print(f'halo zmin = {min(redshift_halo)}, zmax = {max(redshift_halo)}')


    def build_trees(self):
        x_gal = (0.5*np.pi - self.dec_gal) # like theta
        y_gal = np.cos(self.dec_gal) * self.ra_gal # like phi sin(theta)

        x_halo = (0.5*np.pi - self.dec_halo)
        y_halo = np.cos(self.dec_halo) * self.ra_halo

        gal_position = np.dstack([x_gal, y_gal])[0]
        gal_tree = spatial.cKDTree(gal_position)

        halo_position = np.dstack([x_halo, y_halo])[0]
        halo_tree = spatial.cKDTree(halo_position)

        DA_min = min(self.chi_halo)/(1+self.zmin)
        rmax_tree = 4 / DA_min #
        self.indexes_tree = halo_tree.query_ball_tree(gal_tree, r=rmax_tree)
        self.gal_taken = np.zeros(len(self.chi_gal)) # for percolation


    def get_richness(self, i_halo):
        #### step 0: get the tree
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
        DA = self.chi_halo[i_halo] / (1+self.redshift_halo[i_halo]) # Mpc/h
        r = DA * ang_sep # physical Mpc/h

        #### step 3: iteratively calculating r_lambda
        rlam_ini = 1
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
        #print(ngal, len(self.gal_taken[self.gal_taken>1e-4]))
        return rlam, ngal


    def measure_richness(self, pixel_id):
        nh = len(self.ra_halo)
        # print('nh =', nh)
        #### save the richness file:  ####
        outfile1 = open(self.richness_file, 'w')
        outfile1.write('#hid, mass, ra, dec, redshift, rlam, lam, pixelid \n')
        ## TODO: Gladys, please only save the halos in the central pixels (6080) because we need percolation
        for ih in range(nh):
            rlam, lam = self.get_richness(ih)
            if (lam > 0) and (self.pixelid[ih] == pixel_id): #and (self.pixelid == pixel_id)
                outfile1.write('%12i %15e %12g %12g %12g %12g %12g %12i \n'%(self.hid[ih], self.mass[ih], self.ra_halo[ih], 
                                                                        self.dec_halo[ih], self.redshift_halo[ih], rlam,
                                                                             lam, self.pixelid[ih]))
        outfile1.close()

        print('DONE!!')

## TODO: Heidi & Gladys: what to do for halos very close to survey edge??

def run_parallel(pixels):
    for i in range(len(pixels)):
        hpix = pixels[i]
        zmin = 0.2 #0.3
        zmax = 0.35 #0.32
        car = AlternativeRichness(zmin=zmin, zmax=zmax, pixel_id = hpix)
        car.get_galaxies()
        car.get_halos()
        car.build_trees()
        car.measure_richness(hpix)
    
     
procs = 32
df_chunks = np.array_split(np.unique(sky_02_1_all['pixel_id']) ,procs)
if __name__ == '__main__':
    start_time = time.time() 
    with Pool(processes=procs) as p:
        p.map(run_parallel, df_chunks)
        end_time = time.time() 
        elapsed_time = (end_time - start_time)/60
        print(f"Elapsed time: {elapsed_time:.2f} minutes")