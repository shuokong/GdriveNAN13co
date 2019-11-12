# from https://pyregion.readthedocs.io/en/latest/getting_started.html

import pyregion
from astropy.io import fits
import numpy as np

region_name = "ds9.reg"
r = pyregion.open(region_name)

hdu1 = fits.open('nofreq_mom1_SE_mask_13co_pix_2_Tmb.fits')[0]
mymask = r.get_mask(hdu=hdu1)
hdu1.data[~mymask] = 0
hdu1.data[mymask] = 1
fits.writeto('SE_filaments.fits', hdu1.data, hdu1.header, clobber=True)

