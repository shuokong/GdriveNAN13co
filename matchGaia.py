import numpy as np
import os

################## query based on ICRS coordinates 2 arcsec search radius

os.system('cp rebull_raw_ICRS.reg rebull_ICRS_Gaia.reg') # template region file to show Gaia results

rebull_raw_ICRS_Gaia = np.loadtxt('rebull_raw_ICRS_Gaia.txt',dtype='string')
rows, cols = rebull_raw_ICRS_Gaia.shape
distdic = {}
for rr in range(rows):
    kk = str(rebull_raw_ICRS_Gaia[rr,0])
    dd = int(float(rebull_raw_ICRS_Gaia[rr,5]))
    distdic[kk] = str(dd)

rebull_ICRS_Gaia = open('rebull_ICRS_Gaia.reg', 'r')
newregfile = open('ds9.reg','w') 
reglines = rebull_ICRS_Gaia.readlines()
reghead = reglines[:3]
regdata = reglines[3:]
for ll in reghead:
    newregfile.write(ll)
for lll in regdata:
    if lll[6:28] in distdic.keys():
        lll_1 = lll.split('{')[0]
        newregfile.write(lll_1+'{'+distdic[lll[6:28]]+'}\n')

rebull_ICRS_Gaia.close() 
newregfile.close()
os.system('cp ds9.reg rebull_ICRS_Gaia.reg')

################## query based on 2MASS name

os.system('cp rebull_raw_2MASS.reg rebull_2MASS_Gaia.reg') # template region file to show Gaia results

rebull_2MASS_Gaia = np.loadtxt('rebull_2MASS_Gaia.txt',dtype='string')
rows, cols = rebull_2MASS_Gaia.shape
distdic = {}
for rr in range(rows):
    kk = str(rebull_2MASS_Gaia[rr,1]) # Jxxxxxxxx+xxxxxxx
    dd = int(float(rebull_2MASS_Gaia[rr,6])) # rest distance in pc
    distdic[kk] = str(dd)

rebull_2MASS_Gaia = open('rebull_2MASS_Gaia.reg', 'r')
newregfile = open('ds9.reg','w') 
reglines = rebull_2MASS_Gaia.readlines()
reghead = reglines[:3]
regdata = reglines[3:]
for ll in reghead:
    newregfile.write(ll)
for lll in regdata:
    if lll[-19:-2] in distdic.keys(): # Jxxxxxxxx+xxxxxxx
        lll_1 = lll.split('{')[0]
        newregfile.write(lll_1+'{'+distdic[lll[-19:-2]]+'}\n')

rebull_2MASS_Gaia.close() 
newregfile.close()
os.system('cp ds9.reg rebull_2MASS_Gaia.reg')
