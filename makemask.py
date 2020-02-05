# from https://pyregion.readthedocs.io/en/latest/getting_started.html

import pyregion
from astropy.io import fits
import numpy as np

region_name = "F6polygon.reg"
r = pyregion.open(region_name)

temphdu = fits.open('peak_SE_mask_13co_pix_2_Tmb.fits')[0]
temphdu.data = temphdu.data[0,:,:]
temphdu.header['NAXIS'] = 2
del temphdu.header['CRPIX3']
del temphdu.header['CDELT3']
del temphdu.header['CRVAL3']
del temphdu.header['CTYPE3']
#del temphdu.header['NAXIS3']
mymask = r.get_mask(hdu=temphdu)
temphdu.data[mymask] = 1
temphdu.data[~mymask] = 0
fits.writeto('F6polygon.fits', temphdu.data, temphdu.header, clobber=True)

