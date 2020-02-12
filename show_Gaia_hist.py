import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rc
rc('text', usetex=True)
font = {'weight' : 'normal','size':20,'family':'sans-serif','sans-serif':['Helvetica']}
rc('font', **font)
lletter = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']

rebull_ICRS_Gaia = open('rebull_ICRS_Gaia.reg', 'r')
reglines = rebull_ICRS_Gaia.readlines()
reghead = reglines[:3]
regdata = reglines[3:]
allysos = []
for lll in regdata:
    lll_1 = lll.split('{')[1]
    lll_2 = lll_1.split('}')[0]
    allysos.append(int(lll_2))
rebull_ICRS_Gaia.close() 

Pelican = [2434.82178815343, 905.059003320025, 721.636767331034, 749.226180632678, 859.194536987336, 799.190905091087, 503.385643498088, 754.023834174072, 723.031980179476, 1112.82916214578, 881.734886337228, 776.795358396845, 2004.35980118232, 844.395490081145, 691.86095236827, 1906.8920934061, 3786.3162416268, 1621.20279665251, 912.029427018509, 764.500865172174, 1276.44189948567, 760.859339418918, 4815.92561991346, 1976.84505147026, 867.985558166753, 1857.48644346684, 4294.62884562402, 683.454660225107, 2072.06301598032, 723.589234164306, 1407.52398440353, 767.749242809037, 884.941510280575, 625.71852008682, 799.274156847239, 947.906635604635, 1836.35881258639, 767.859569170012, 777.875014853869, 838.7365552612, 811.271278140984, 3316.71997563507, 571.499135106475, 761.789278572869, 852.464752160486, 746.56779297015, 868.356079643741, 975.507445277943, 803.085869514186, 250.575126481165, 872.705823556434, 2886.54497003378, 1729.23356024887, 816.946539703473, 922.625248503084, 915.241252838128, 1744.67715661054, 1202.42423131564, 825.961796720554, 842.362114614054, 2702.05942495823, 752.653489941106, 662.708929034271, 861.393179237551, 1689.5969389148, 905.719974246922, 771.121758631907, 763.295172589652, 3563.44414247642, 952.227550060169, 1166.21750505641, 2873.45026752559, 778.198043080506, 769.115383401083, 667.586655493889, 697.198182147022, 783.292381232644, 2609.32973883258, 734.453543456201, 1894.23793871673, 2820.27045864306, 783.153545148421, 1116.9048077341, 779.822688055989, 891.965991458161, 4235.85479184503, 1381.95024716679, 792.707595665285, 751.644769756222, 725.283061764409, 898.398469553683, 888.211145556609, 1614.8505711268, 852.749153286709, 1971.73194478772, 698.004738962798, 1914.3872783759, 746.235265614942, 1940.07054777389, 731.702595937933, 805.315168356481, 725.781888457524, 816.060742167335, 910.139591795079, 468.481219866894, 408.905972253502, 863.140191847476, 854.500492363945, 851.685110466478, 731.466299108073, 651.354419894784, 803.436474601046, 1413.26370047616, 755.015499445039, 1518.58221982549, 387.745226340284, 768.029159424976, 790.859947126273, 791.763841124086, 696.569654699851, 750.168968429683, 522.486475485046, 969.746826475215, 379.471880648502, 1011.54258529349, 883.874804473825, 496.918542803599, 348.068414058764, 820.158228581386, 705.374807929477, 1006.32419140855, 813.160739405869, 725.533319840291, 1856.23948618869, 773.803892256781, 729.938541863254, 1878.69197479626, 1529.33821022423, 800.332653081625, 819.28584329989, 3801.92635113712, 866.518362616397, 1491.80346727631, 813.407947346681, 973.49472654731, 1768.42607842712, 799.971428749532, 1978.61611965194, 823.660960882639, 813.504648748231, 2316.65195121524, 506.973840186841, 901.252906294015, 1527.20566197042, 1370.39651181901, 852.56395685007, 519.745495926029, 734.299485818984, 552.51441984638, 4015.64353312615, 502.327363198938, 834.077244783787, 904.117685510665, 860.396721264234, 889.43392037579, 479.301138340759, 798.953714765902, 3399.63324446406, 448.441109886709, 651.256153827257, 727.256963451992, 961.34755989208, 740.026482583995, 879.906501680387]
NGC6997 = [860.140030977463, 892.10949622255, 902.514505400935, 847.439810238413, 881.271633610038, 842.57143489398]
PelHat = [939.227210255998, 965.221091334315, 844.948483492299, 1525.86512482428, 767.777274789607, 1790.71330762667, 888.933273391326, 2046.89948760523, 562.646293570125]
Gulf = [1103.99415377816, 1863.35141353043, 1111.28517301647, 1611.91434721785, 3565.55158356839, 2457.71518647829, 750.160710341511, 495.573978317199, 1423.03126492669, 661.633415573617, 397.55708709332, 2755.64599505009, 2618.06423243466, 681.150687550381, 844.587985759287, 383.443796359857, 4272.83322578446, 874.668047171683, 739.297626106859, 735.493266075061, 1403.08712213707, 696.852438516762, 1461.8024259132, 503.260569475898, 420.95242199612, 653.875641700789, 1234.41772734148, 737.695310637663, 907.298133807293, 659.314530495693, 765.205417913019, 673.771117269154, 1445.38590759965, 1659.01834881985, 754.449442919191, 741.186882205338, 912.289823263289, 734.938692934406, 1000.70468256195, 759.291454723858, 1020.17613692611, 1683.48266654352, 697.94215429709, 567.702706482062, 1534.29426126875, 779.438243738859]

clusterlist = [allysos,Gulf,Pelican,PelHat]
titlelist = ['All','Gulf','Pelican','PelHat']

xpanels = 4
ypanels = 1
tkin_mean = np.nanmean(allysos)
tkin_median = np.nanmedian(allysos)
tkin_min = np.nanmin(allysos)
tkin_max = np.nanmax(allysos)
print 'total num',len(allysos),'Tkin, mean',tkin_mean,'median',tkin_median,'min',tkin_min,'max',tkin_max
tlow,thigh = (100,1500)
datafiles = {}
for i in range(0,xpanels):
    for j in range(0,ypanels):
        panel = i+j*xpanels+1
        print 'panel',panel 
        hist, bin_edges = np.histogram(clusterlist[panel-1],bins='auto',range=(tlow,thigh))
        bincenter = (bin_edges[:-1] + bin_edges[1:]) / 2.
        datafiles['panel'+str(panel)] = {'title':titlelist[panel-1],'lines':{
                                                         '1':{'x':bincenter,'y':hist,'peaksnr':[],'legends':'data','linestyles':'k-','drawsty':'steps-mid'},
                                                         #'2':{'x':bincenter,'y':gaus(bincenter,*popt),'peaksnr':[],'legends':'fit '+r'$\sigma$='+'{0:.3f}'.format(popt[2])+r'$\pm$'+'{0:.3f}'.format(perr[2]),'linestyles':'b-','drawsty':'default'},
                                                            },
                                          'xlim':[tlow,thigh],'ylim':[0,np.nanmax(hist)*1.1],'xscale':'linear','yscale':'linear',
                                        'xlabel':r'$\rm Gaia~Distance~\rm (pc)$','ylabel':r'$\rm number$','text':'','vertlines':[],
                                        }

xpanelwidth = 7
ypanelwidth = 6
fig=plt.figure(figsize=(xpanelwidth*xpanels,ypanelwidth*ypanels))
plt.subplots_adjust(wspace=0.1,hspace=0.1)
pdfname='Rebull_Gaia_hist.pdf'
for i in range(0,xpanels):
    for j in range(0,ypanels):
        panelnum = i+j*xpanels+1
        ax = fig.add_subplot(ypanels,xpanels,panelnum)
        if 'panel'+str(panelnum) not in datafiles.keys(): continue
        ax.set_xscale(datafiles['panel'+str(panelnum)]['xscale']) 
        ax.set_yscale(datafiles['panel'+str(panelnum)]['yscale']) 
        if datafiles['panel'+str(panelnum)]['ylim'] != []:
            ydown = datafiles['panel'+str(panelnum)]['ylim'][0]
            yup   = datafiles['panel'+str(panelnum)]['ylim'][1]
            ax.set_ylim(ydown,yup)
        if datafiles['panel'+str(panelnum)]['xlim'] != []:
            xmin,xmax = datafiles['panel'+str(panelnum)]['xlim']
            ax.set_xlim(xmin,xmax)
            #ax.hlines(0,xmin,xmax,linestyle='dotted')
        for datafilenum in range(len(datafiles['panel'+str(panelnum)]['lines'].keys())): 
            x = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['x']
            y = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['y']
            #xx = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['velocity']
            #print 'weird cores',corenames[(xx<-2)|(xx>2)]
            #ax.hist(x,bins='auto',range=(xmin,xmax))
            linestyle = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['linestyles']
            legend = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['legends']
            drawsty = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['drawsty']
            ax.plot(x,y,linestyle,label=legend,drawstyle=drawsty)
            #ax.text(peakvelocity+0.8, yup*0.9, '%.1f' % peakvelocity + ',' + '%.1f' % peaksnr,horizontalalignment='left',verticalalignment='center',fontsize=12)
        #ax.legend(frameon=False,prop={'size':14},labelspacing=0.2) 
        if j == 0:
            ax.set_title(r'')
        ax.text(0.05, 0.9,datafiles['panel'+str(panelnum)]['title'],horizontalalignment='left',verticalalignment='center',transform = ax.transAxes)
        #ax.text(0.1, 0.9,datafiles['panel'+str(panelnum)]['title']+' '+datafiles['panel'+str(panelnum)]['text'],horizontalalignment='left',verticalalignment='center',transform = ax.transAxes,fontsize=12)
        ax.text(0.9, 0.9,'('+lletter[panelnum-1]+')',horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
        xlabel = datafiles['panel'+str(panelnum)]['xlabel']
        ylabel = datafiles['panel'+str(panelnum)]['ylabel']
        vertlinex = datafiles['panel'+str(panelnum)]['vertlines']
        for vl in vertlinex:
            ax.vlines(vl,ydown,yup,linestyles='dashed',colors='k')
        if j != ypanels-1:
            ax.set_yticks(ax.get_yticks()[1:])
            ax.set_xticklabels(ax.get_xlabel(),visible=False)
        else: 
            ax.set_xlabel(xlabel)
        if i != 0:
            #ax.set_yticklabels(ax.get_ylabel(),visible=False) 
            ax.set_xticks(ax.get_xticks()[1:]) 
        else: 
            ax.set_ylabel(ylabel)
        #minor_locator = AutoMinorLocator(5)
        #ax.xaxis.set_minor_locator(minor_locator)
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.tick_params(axis='both',direction='in',length=5,which='major',top=True,right=True)
        ax.tick_params(axis='both',direction='in',length=3,which='minor',top=True,right=True)
        
os.system('rm '+pdfname)
plt.savefig(pdfname,bbox_inches='tight')
plt.close(fig)
os.system('open '+pdfname)
os.system('cp '+pdfname+os.path.expandvars(' ~/GoogleDrive/imagesCARMANan/'))
