import aplpy
import astropy.wcs as wcs
import numpy as np
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib import rc
rc('text', usetex=True)
font = {'weight' : 'normal','size':20,'family':'sans-serif','sans-serif':['Helvetica']}
rc('font', **font)

cellsize = 2. # arcsec
distance = 700. # parsec
cellpc = cellsize * distance / 206265. # cell size in pc
sehdu1 = fits.open('mom0_SE_mask_13co_pix_2_Tmb.fits')[0]
nwhdu1 = fits.open('mom0_NW_mask_13co_pix_2_Tmb.fits')[0]

makeYSOmap = 0
SEmom0 = 1
NWmom0 = 1

if makeYSOmap:

    wse = wcs.WCS(sehdu1.header)
    maskse = np.full_like(sehdu1.data, np.nan, dtype=np.double)
    
    wnw = wcs.WCS(nwhdu1.header)
    masknw = np.full_like(nwhdu1.data, np.nan, dtype=np.double)
    
    rebull = np.loadtxt('rebull_plusO.txt')
    nyso, ncoords = rebull.shape
    worldcoord = np.pad(rebull, ((0, 0), (0, 1)), 'constant', constant_values=((0, 0), (0, 0)))
    
    pixcoordse = wse.all_world2pix(worldcoord,0) # [[x,y],[x,y]]
    pixcoordnw = wnw.all_world2pix(worldcoord,0)
    
    for ind, vv in np.ndenumerate(masknw):
        if np.isnan(nwhdu1.data[ind]): continue
        zz, yy, xx = ind
        dists = np.nansum((pixcoordnw[:,:2] - [xx,yy])**2, axis=1) 
        dist5 = sorted(dists)[4]**0.5 # in cells 
        masknw[ind] = 5./(np.pi*(dist5*cellpc)**2) # stars per pc2 
    
    nw_mask_hdu = fits.PrimaryHDU(masknw, nwhdu1.header)
    nw_mask_hdu.writeto('NW_SigmaYSO.fits',output_verify='exception',clobber=True,checksum=False)
    
    for ind, vv in np.ndenumerate(maskse):
        if np.isnan(sehdu1.data[ind]): continue
        zz, yy, xx = ind
        dists = np.nansum((pixcoordse[:,:2] - [xx,yy])**2, axis=1) 
        dist5 = sorted(dists)[4]**0.5 # in cells 
        maskse[ind] = 5./(np.pi*(dist5*cellpc)**2) # stars per pc2 
    
    se_mask_hdu = fits.PrimaryHDU(maskse, sehdu1.header)
    se_mask_hdu.writeto('SE_SigmaYSO.fits',output_verify='exception',clobber=True,checksum=False)
    sys.exit()


if SEmom0 == 1:
    mask_sehdu1 = fits.open('SE_SigmaYSO.fits')[0]
    xcenter = 314.0150412
    ycenter = 43.70085501
    wid = 1.1022223
    hei = 0.7527778
    xpanels = 1
    ypanels = 1
    fig=plt.figure(figsize=(3*xpanels*1.1*(wid/(wid+hei))*10.,3*ypanels/1.1*(hei/(wid+hei))*10.))
    ff = aplpy.FITSFigure(sehdu1, figure=fig)
    ff.recenter(xcenter,ycenter,width=wid,height=hei) 
    ff.set_theme('publication')
    ff.set_system_latex(True)
    maxcolor = 80.
    mincolor = 5.
    ff.show_colorscale(cmap='afmhot', vmin=mincolor, vmax=maxcolor, stretch='sqrt')
    ff.show_contour(data=mask_sehdu1, levels=np.arange(50.,1200.,50.), colors='white', linewidths=0.5) 
    #ff.show_regions('olay.reg')
    #ff.show_regions('olay1.reg')
    ff.add_colorbar() 
    ff.colorbar.set_pad(0.5)
    ff.colorbar.set_axis_label_text('K km s$^{-1}$')
    ff.add_scalebar(0.1,corner='top right',pad=1) # degree for 1pc at 550 pc
    ff.scalebar.set_label('1 pc') 
    beamx = xcenter + wid/2.*1.2
    beamy = ycenter - hei/2.*9./10.
    bmaj = sehdu1.header['BMAJ']
    bmin = sehdu1.header['BMIN']
    beamangle = sehdu1.header['BPA'] 
    ff.show_ellipses(beamx,beamy,bmaj,bmin,angle=beamangle-90,facecolor='black',edgecolor='black')
    textx = xcenter - wid/2.*9./10.
    texty = ycenter - hei/2.*9./10.
    ff.add_label(textx,texty,'0th-moment $^{13}$CO(1-0)',weight='bold')
    #ff.tick_labels.set_xformat('dd')
    #ff.tick_labels.set_yformat('dd')
    pdfname = 'SE_yso_13co.pdf'
    os.system('rm '+pdfname)
    plt.savefig(pdfname,bbox_inches='tight')
    os.system('open '+pdfname)
    os.system('cp '+pdfname+os.path.expandvars(' /Users/shuokong/GoogleDrive/imagesCARMANan/'))

if NWmom0 == 1:
    mask_nwhdu1 = fits.open('NW_SigmaYSO.fits')[0]
    xcenter = 312.5676799
    ycenter = 44.5018964
    wid = 0.3561111
    hei = 0.5444445
    xpanels = 1
    ypanels = 1
    fig=plt.figure(figsize=(3*xpanels*1.1*(wid/(wid+hei))*10.,3*ypanels/1.1*(hei/(wid+hei))*10.))
    ff = aplpy.FITSFigure(nwhdu1, figure=fig)
    ff.recenter(xcenter,ycenter,width=wid,height=hei) 
    ff.set_theme('publication')
    #ff.set_system_latex(True)
    maxcolor = 80.
    mincolor = 5.
    ff.show_colorscale(cmap='afmhot', vmin=mincolor, vmax=maxcolor, stretch='sqrt')
    ff.show_contour(data=mask_nwhdu1, levels=np.arange(50.,1200.,50.), colors='white', linewidths=0.5) 
    #ff.show_regions('olay.reg')
    #ff.show_regions('olay1.reg')
    ff.add_colorbar() 
    ff.colorbar.set_pad(0.5)
    ff.colorbar.set_axis_label_text('K km s$^{-1}$')
    ff.add_scalebar(0.05,corner='bottom right',pad=1) # degree for 0.5 pc at 550 pc
    ff.scalebar.set_label('0.5 pc') 
    beamx = xcenter + wid/2.*1.3
    beamy = ycenter - hei/2.*9./10.
    bmaj = nwhdu1.header['BMAJ']
    bmin = nwhdu1.header['BMIN']
    beamangle = nwhdu1.header['BPA'] 
    ff.show_ellipses(beamx,beamy,bmaj,bmin,angle=beamangle-90,facecolor='black',edgecolor='black')
    textx = xcenter + wid/2.*4./5.
    texty = ycenter + hei/2.*9./10.
    ff.add_label(textx,texty,'0th-moment $^{13}$CO(1-0)',weight='bold')
    #ff.tick_labels.set_xformat('dd')
    #ff.tick_labels.set_yformat('dd')
    pdfname = 'NW_yso_13co.pdf'
    os.system('rm '+pdfname)
    plt.savefig(pdfname,bbox_inches='tight')
    os.system('open '+pdfname)
    os.system('cp '+pdfname+os.path.expandvars(' /Users/shuokong/GoogleDrive/imagesCARMANan/'))


