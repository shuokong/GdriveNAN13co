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
SEmom0 = 0
NWmom0 = 0
SEplot = 1
NWplot = 1

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

if NWplot == 1:
    mask_nwhdu1 = fits.open('NW_SigmaYSO.fits')[0]
    ########################
    plot_xlist = [nwhdu1.data.flatten()]
    plot_ylist = [mask_nwhdu1.data.flatten()]
    
    plot_xlabel_list = [r'$\Sigma_{\rm total}~\rm (M_\odot~pc^{-2}$']
    plot_ylabel_list = [r'$\Sigma_{\rm YSO}~\rm (stars~pc^{-2}$']
    
    xpanels = 1
    ypanels = 1
    xpanelwidth = 5
    ypanelwidth = 5
    
    datafiles = {}
    pdfname = 'NW_Sigma_yso_13co.pdf'
    for i in range(0,xpanels):
        for j in range(0,ypanels):
            panel = i+j*xpanels+1
            print 'panel',panel 
            velocity = plot_xlist[panel-1]
            rawintens = plot_ylist[panel-1]
            datafiles['panel'+str(panel)] = {'title':'','lines':{'1':{'x':velocity,'y':rawintens,'legends':'data','drawsty':'default',}},'ylim':[1.e1,1.e3],'xlim':[1.,1.e2],'xscale':'log','yscale':'log','xlabel':plot_xlabel_list[panel-1],'ylabel':plot_ylabel_list[panel-1],'text':''}
    fig, axes = plt.subplots(ncols=xpanels, nrows=ypanels, figsize=(5*xpanels,5*ypanels))
    plt.subplots_adjust(wspace=0.1,hspace=0.1)
    for i in range(0,xpanels):
        for j in range(0,ypanels):
            panelnum = i+j*xpanels+1
            #ax = axes[j,i]
            ax = axes
            ax.set_xscale(datafiles['panel'+str(panelnum)]['xscale']) 
            ax.set_yscale(datafiles['panel'+str(panelnum)]['yscale']) 
            for datafilenum in range(len(datafiles['panel'+str(panelnum)]['lines'].keys())): 
                x = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['x']
                y = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['y']
                legend = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['legends']
                drawsty = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['drawsty']
                ax.scatter(x,y,label=legend)
            #ax.legend(frameon=False,prop={'size':14},labelspacing=0.1) 
            #if j == 0:
            #    ax.set_title(datafiles['panel'+str(panelnum)]['title'])
            #ax.text(0.1, 0.9,datafiles['panel'+str(panelnum)]['title']+' '+datafiles['panel'+str(panelnum)]['text'],horizontalalignment='left',verticalalignment='center',transform = ax.transAxes,fontsize=12)
            #ax.text(0.05, 0.95,'('+lletter[panelnum-1]+')',horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
            xlabel = datafiles['panel'+str(panelnum)]['xlabel']
            ylabel = datafiles['panel'+str(panelnum)]['ylabel']
            #ax.set_xticks(np.arange(datafiles['panel'+str(panelnum)]['xlim'][0],datafiles['panel'+str(panelnum)]['xlim'][1],0.1),minor=True)
            #for vl in vertlinex:
            #    ax.vlines(vl,ydown,yup,linestyles='dashed')
            #ax.hlines(2,datafiles['panel'+str(panelnum)]['xlim'][0],datafiles['panel'+str(panelnum)]['xlim'][1],linestyle='dashed')
            if j != ypanels-1:
                ax.set_yticks(ax.get_yticks()[1:])
                ax.set_xticklabels(ax.get_xlabel(),visible=False)
            else: 
                ax.set_xlabel(xlabel)
            if i != 0:
                ax.set_yticklabels(ax.get_ylabel(),visible=False) 
                ax.set_xticks(ax.get_xticks()[1:]) 
            else: 
                ax.set_ylabel(ylabel)
            if datafiles['panel'+str(panelnum)]['ylim']:
                ydown = datafiles['panel'+str(panelnum)]['ylim'][0]
                yup   = datafiles['panel'+str(panelnum)]['ylim'][1]
                ax.set_ylim(ydown,yup)
            #else:
            #    ydown,yup = ax.get_ylim()
            if datafiles['panel'+str(panelnum)]['xlim']:
                ax.set_xlim(datafiles['panel'+str(panelnum)]['xlim'][0],datafiles['panel'+str(panelnum)]['xlim'][1])
                #ax.hlines(0,datafiles['panel'+str(panelnum)]['xlim'][0],datafiles['panel'+str(panelnum)]['xlim'][1],linestyle='dotted')
            
    os.system('rm '+pdfname)
    plt.savefig(pdfname,bbox_inches='tight')
    plt.close(fig)
    os.system('open '+pdfname)
    os.system('cp '+pdfname+os.path.expandvars(' /Users/shuokong/GoogleDrive/imagesCARMANan/'))



