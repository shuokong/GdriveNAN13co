import numpy as np
import os

fk5gaiacoord,gaiadist = np.loadtxt('rebull_FK5_Gaia_results.txt',delimiter=';',dtype='string',unpack=True,usecols=(0,5))
gaiadict = dict(zip(fk5gaiacoord,gaiadist))
cluster,coord = np.loadtxt('rebull_raw_cluster.txt',delimiter=';',dtype='string',unpack=True,usecols=(1,2))

clusterdict = {'Pelican':[],'PelHat':[],'Gulf':[],'NGC6997':[]}
for nn,cc in enumerate(cluster):
    if cc.strip() == '': continue 
    if coord[nn] in fk5gaiacoord:
        dist = float(gaiadict[coord[nn]]) 
        clusterdict[cc.strip()].append(dist) 
for kk in clusterdict.keys():
     print kk,clusterdict[kk]

