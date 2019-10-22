import numpy as np
import os

os.system('cp rebull_raw_ICRS.reg rebull_ICRS_Gaia.reg')

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
