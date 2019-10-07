import aplpy
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib import rc
import sys
rc('text', usetex=True)
font = {'weight' : 'normal','size':20,'family':'sans-serif','sans-serif':['Helvetica']}
rc('font', **font)

## SE
hdu1 = fits.open('mom2_SE_mask_13co_pix_2_Tmb.fits')[0]
xcenter = 314.0150412
ycenter = 43.70085501
wid = 1.1022223
hei = 0.7527778
xpanels = 1
ypanels = 1
fig=plt.figure(figsize=(3*xpanels*1.1*(wid/(wid+hei))*10.,3*ypanels/1.1*(hei/(wid+hei))*10.))
ff = aplpy.FITSFigure(hdu1, figure=fig)
ff.recenter(xcenter,ycenter,width=wid,height=hei) 
ff.set_theme('publication')
#ff.set_system_latex(True)
maxcolor = np.nanmax(hdu1.data)
ff.show_colorscale(cmap='jet', vmin=0, vmax=10, stretch='linear')
#ff.show_regions('olay2.reg')
#ff.show_contour(mask_hdu, levels=1, colors='yellow', linewidths=0.1)
ff.add_colorbar() 
ff.colorbar.set_pad(0.5)
ff.colorbar.set_axis_label_text('km s$^{-1}$')
ff.add_scalebar(0.1,corner='top right',pad=1) # degree for 1pc at 550 pc
ff.scalebar.set_label('1 pc') 
beamx = xcenter + wid/2.*1.3
beamy = ycenter - hei/2.*9./10.
bmaj = hdu1.header['BMAJ']
bmin = hdu1.header['BMIN']
beamangle = hdu1.header['BPA'] 
ff.show_ellipses(beamx,beamy,bmaj,bmin,angle=beamangle-90,facecolor='black',edgecolor='black')
textx = xcenter - wid/2.*9./10.
texty = ycenter - hei/2.*9./10.
ff.add_label(textx,texty,'1st-moment $^{13}$CO(1-0)',weight='bold')
pdfname = 'mom2_SE_13co_pix_2_Tmb.pdf'
os.system('rm '+pdfname)
plt.savefig(pdfname,bbox_inches='tight')
os.system('open '+pdfname)
os.system('cp '+pdfname+os.path.expandvars(' /Users/shuokong/GoogleDrive/imagesCARMANan/'))

## NW
hdu1 = fits.open('mom2_NW_mask_13co_pix_2_Tmb.fits')[0]
xcenter = 312.5676799
ycenter = 44.5018964
wid = 0.3561111
hei = 0.5444445
xpanels = 1
ypanels = 1
fig=plt.figure(figsize=(3*xpanels*1.1*(wid/(wid+hei))*10.,3*ypanels/1.1*(hei/(wid+hei))*10.))
ff = aplpy.FITSFigure(hdu1, figure=fig)
ff.recenter(xcenter,ycenter,width=wid,height=hei) 
ff.set_theme('publication')
#ff.set_system_latex(True)
maxcolor = np.nanmax(hdu1.data)
ff.show_colorscale(cmap='jet', vmin=0, vmax=10, stretch='linear')
#ff.show_regions('olay2.reg')
#ff.show_contour(mask_hdu, levels=1, colors='yellow', linewidths=0.1)
ff.add_colorbar() 
ff.colorbar.set_pad(0.5)
ff.colorbar.set_axis_label_text('km s$^{-1}$')
ff.add_scalebar(0.05,corner='bottom right',pad=1) # degree for 0.5 pc at 550 pc
ff.scalebar.set_label('0.5 pc') 
beamx = xcenter + wid/2.*1.2
beamy = ycenter - hei/2.*9./10.
bmaj = hdu1.header['BMAJ']
bmin = hdu1.header['BMIN']
beamangle = hdu1.header['BPA'] 
ff.show_ellipses(beamx,beamy,bmaj,bmin,angle=beamangle-90,facecolor='black',edgecolor='black')
textx = xcenter + wid/2.*4./5.
texty = ycenter + hei/2.*9./10.
ff.add_label(textx,texty,'1st-moment $^{13}$CO(1-0)',weight='bold')
pdfname = 'mom2_NW_13co_pix_2_Tmb.pdf'
os.system('rm '+pdfname)
plt.savefig(pdfname,bbox_inches='tight')
os.system('open '+pdfname)
os.system('cp '+pdfname+os.path.expandvars(' /Users/shuokong/GoogleDrive/imagesCARMANan/'))
