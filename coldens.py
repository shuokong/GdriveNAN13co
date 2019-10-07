import sys
from astropy.io import fits
import numpy as np

rregion = 'SE'
rregion = 'NW'

hdu1 = fits.open('mom0_'+rregion+'_mask_13co_pix_2_Tmb.fits')[0]
hdu2 = fits.open('peak_'+rregion+'_mask_13co_pix_2_Tmb.fits')[0]
hdu3 = fits.open(rregion+'_tex12.fits')[0]
hdu4 = fits.open('peak_'+rregion+'_mask_13co_pix_2_Tmb.fits')[0]
hdu5 = fits.open(rregion+'_mask_13co_pix_2_Tmb.fits')[0]

mom013data = hdu1.data
peak13data = hdu2.data
tex12data = hdu3.data
cube13 = hdu5.data[0,:,:,:]
print cube13.shape
n1,n2,n3 = cube13.shape
dv = 0.16 # km/s
cuberms = 0.7 # K
mom0sens = cuberms * dv * (n1)**0.5 # integrate over total channels
sensbool = np.copy(mom013data)
sensbool[mom013data < 3.*mom0sens] = np.nan
sensbool[mom013data >= 3.*mom0sens] = 1
mask = 1
#sys.exit()

tau13peak = -np.log(1. - (peak13data / 5.3) / (1. / (np.exp(5.3 / tex12data) - 1) - 0.16))
hdu4.data = -np.log(1. - (peak13data / 5.3) / (1. / (np.exp(5.3 / tex12data) - 1) - 0.16))
if mask == 1:
    hdu4.data = hdu4.data * sensbool
hdu4.writeto(rregion+'tau13peak.fits', output_verify='exception', overwrite=True, checksum=False)

## coldens13_tauthin.fits
coldens13_thin = 3.0e14 * mom013data / (1. - np.exp(-5.3 / tex12data))
hdu4.data = 3.0e14 * mom013data / (1. - np.exp(-5.3 / tex12data))
if mask == 1:
    hdu4.data = hdu4.data * sensbool
hdu4.writeto(rregion+'coldens13_thin.fits', output_verify='exception', overwrite=True, checksum=False)

## coldens13_taupeak.fits
hdu4.data = tau13peak / (1-np.exp(-tau13peak)) * coldens13_thin
if mask == 1:
    hdu4.data = hdu4.data * sensbool
hdu4.writeto(rregion+'coldens13_taupeak.fits', output_verify='exception', overwrite=True, checksum=False)

## coldens13_tauinte.fits
tex12cube = np.repeat(tex12data,n1,axis=0)
print tex12cube.shape
tau13 = -np.log(1. - (cube13 / 5.3) / (1. / (np.exp(5.3 / tex12cube) - 1) - 0.16))
print tau13.shape
tau13_inte = np.nansum(tau13,axis=0,keepdims=True) * dv
print tau13.shape
one_minus_tau13_inte = np.nansum(1.-np.exp(-tau13),axis=0,keepdims=True)
print one_minus_tau13_inte.shape
#hdu4.data = tau13_inte / one_minus_tau13_inte * coldens13_thin
hdu4.data = 2.42e14 * (tex12data + 0.88) / (1. - np.exp(-5.3 / tex12data)) * tau13_inte
if mask == 1:
    hdu4.data = hdu4.data * sensbool
hdu4.writeto(rregion+'coldens13_tauinte.fits', output_verify='exception', overwrite=True, checksum=False)

