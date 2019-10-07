import numpy as np
import pyfits

hdulist = pyfits.open('pmo13_cube_13CO.fits')
hdulist[0].data=hdulist[0].data[0,:,:,:]/0.48
hdulist[0].header['NAXIS'] = 3
hdulist[0].header['BMAJ'] = 0.015222
hdulist[0].header['BMIN'] = 0.015222
hdulist[0].header['BUNIT'] = 'K'
hdulist[0].header['CTYPE3'] = 'VELO-LSR'
del hdulist[0].header['DATAMIN'] 
del hdulist[0].header['DATAMAX'] 
#del hdulist[0].header['NAXIS4'] 
del hdulist[0].header['CRPIX4'] 
del hdulist[0].header['CDELT4'] 
del hdulist[0].header['CRVAL4'] 
del hdulist[0].header['CTYPE4'] 
del hdulist[0].header['CROTA4'] 
del hdulist[0].header['DATE']
del hdulist[0].header['ORIGIN']
del hdulist[0].header[''] 
hdulist.writeto('cube_13CO.fits',output_verify='exception',clobber=True,checksum=False)


