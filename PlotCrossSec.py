
import matplotlib.pyplot as plt
import numpy as np
import matplotlib

def PlotCrossSec(cgrid,ind):

    matplotlib.rcParams['image.cmap'] = 'jet'

    fig,ax = plt.subplots(3,1,figsize=(6,9), constrained_layout=True,
        sharex=True, sharey=True)
    x = cgrid['alongx'][ind]
    x0=x
    x = x0[:-1]+np.diff(x)/2.
    x = np.concatenate(([x[0]-3.],x,[x[-1]+3.]))
    denc = np.arange(21, 28, 0.25)
    dencl = np.arange(21, 28, 1.)
    slim = [29, 33.]
    tlim = [8, 14.]
    olim = [0, 60.]

    pc=ax[0].pcolormesh(x,cgrid['depths'],cgrid['sal'][:,ind],rasterized=True,vmin=slim[0],vmax=slim[1])
    ax[0].contour(x0,cgrid['depths'],cgrid['pden'][:,ind]-1000.,denc,colors='0.5',linewidths=0.2)
    cs=ax[0].contour(x0,cgrid['depths'],cgrid['pden'][:,ind]-1000.,dencl,colors='0')
    for nn, xx in enumerate(cgrid['alongx']):
        print( cgrid['id'][nn].data, xx)
        ax[0].text(xx, 280, cgrid['id'][nn].data, ha='center')
    ax[0].clabel(cs,fmt='%1.0f', fontsize=9)
    ax[0].set_ylim(280,0)
    ax[0].set_title('S [psu]')
    ax[0].set_ylabel('DEPTH [m]')
    fig.colorbar(pc,ax=ax[0],shrink=0.6,extend='both')

    pc=ax[1].pcolormesh(x,cgrid['depths'],cgrid['temp'][:,ind],rasterized=True,vmin=tlim[0],vmax=tlim[1])
    ax[1].contour(x0,cgrid['depths'],cgrid['pden'][:,ind]-1000.,denc,colors='0.5',linewidths=0.2)
    ax[1].contour(x0,cgrid['depths'],cgrid['pden'][:,ind]-1000.,dencl,colors='0')
    ax[1].set_ylim(280,0)
    ax[1].set_title('$T\ [^oC]$')
    fig.colorbar(pc,ax=ax[1],shrink=0.6,extend='both')

    pc=ax[2].pcolormesh(x,cgrid['depths'],np.ma.masked_invalid(cgrid['O2'][:,ind]),rasterized=True,vmin=olim[0],vmax=olim[1])
    ax[2].contour(x0,cgrid['depths'],cgrid['pden'][:,ind]-1000., denc, colors='0.1',linewidths=0.2)
    ax[2].contour(x0,cgrid['depths'],cgrid['pden'][:,ind]-1000., dencl, colors='0')
    ax[2].set_ylim(280,0)
    ax[2].set_title('$O2 $')
    ax[2].set_xlabel('X [km] from S4')
    ax[2].set_xlim([-3,35])
    fig.colorbar(pc,ax=ax[2],shrink=0.6,extend='both')

    return fig,ax
#fig.savefig('2016CrossSec.pdf')
