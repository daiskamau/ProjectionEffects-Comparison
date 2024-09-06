import numpy as np
import pandas as pd
import healpy as hp

import GCRCatalogs as gcr
from astropy.table import Table
from astropy.coordinates import SkyCoord
from multiprocessing import Pool, cpu_count


c = 3e5
G = 4.302e-9
nside = 32
h=0.71
skysim = gcr.load_catalog('skysim5000_v1.1.1')
skycosmology = skysim.cosmology # Cosmological parameters
# cosmoskysim = FlatLambdaCDM(H0=skycosmology.H0.value, Om0= skycosmology.Om0)


## Get the all galaxies in the central and neighbouring pixels
def galaxies(hpix,halos):
    hpix_list = skysim.available_healpix_pixels
    hpix_neighbors = hp.pixelfunc.get_all_neighbours(nside, hpix)
    hpix_toread = list(set(hpix_list) & set(hpix_neighbors))
    hpix_toread.append(hpix)
    # print(hpix_toread)

    ## Getting the galaxies
    ra_all = []
    dec_all = []
    z_all = []
    mag_g_all = []
    mag_r_all = []
    mag_i_all = []
    mag_z_all = []
    mag_y_all = []

    ## getting galaxies
    for hpix_i in hpix_toread:
        galaxy_data = skysim.get_quantities(['redshift','ra', 'dec','mag_true_g','mag_true_r',
                                             'mag_true_i','mag_true_z','mag_true_y'], 
                                            native_filters=f'healpix_pixel == {hpix_i}')

        ra_all.extend(galaxy_data['ra'])
        dec_all.extend(galaxy_data['dec'])
        z_all.extend(galaxy_data['redshift'])
        mag_g_all.extend(galaxy_data['mag_true_g'])
        mag_r_all.extend(galaxy_data['mag_true_r'])
        mag_i_all.extend(galaxy_data['mag_true_i'])
        mag_z_all.extend(galaxy_data['mag_true_z'])
        mag_y_all.extend(galaxy_data['mag_true_y'])

    chi_all = skycosmology.comoving_distance(z_all).value * h #Mpc/h
    galaxies_all = np.array([ra_all,dec_all,z_all,mag_g_all,mag_r_all,mag_i_all,mag_z_all,mag_y_all,chi_all])



    ### Halos within the central pixel_id
    pixel_halos = halos.to_pandas()[pd.Series(halos['pixel_id']).isin(hpix_toread)]

    ## Get the halo properties
    hid = pixel_halos['halo_id']
    redshift_halo = pixel_halos['redshift'] ## redshift (with LOS velocity) instead of redshift_true to be consistent with lensing
    mass = pixel_halos['baseDC2/sod_halo_mass'] # M200c  # TODO: some central halos have mass = -1 (ask on SLACK)
    ra_halo = pixel_halos['ra']
    dec_halo = pixel_halos['dec']
    pixelid = pixel_halos['pixel_id']
    chi_halo = pixel_halos['comoving_dis']
    DA_halo = pixel_halos['Angular_dis']

    pixel_halos = np.array([hid,redshift_halo,mass,ra_halo,dec_halo,pixelid,chi_halo,DA_halo])

    return galaxies_all, pixel_halos
