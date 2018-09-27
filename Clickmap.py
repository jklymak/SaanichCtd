import numpy as np
import matplotlib.pyplot as plt
import matfile
import seawater.eos80 as sw
import seawater.constants as swcons

# plot topo
cc=matfile.loadmatbunch('201710/CtdGrid.mat')
cgrid = cc['cgrid']

dat=matfile.loadmatbunch('../../topo/SouthVanIsle.mat')
topo = dat['VanIsleTopo']
print(topo.keys())

fig,ax = plt.subplots()
ax.plot(cgrid.lon,cgrid.lat,'rd')
ax.get_xaxis().get_major_formatter().set_useOffset(False)
#ax.set_aspect(1./np.cos(48.7*np.pi/180.))
ax.contour(topo['Lon'],topo['Lat'],topo['z'],[-1000,0],colors='k')
pc=ax.pcolormesh(topo['Lon'],topo['Lat'],topo['z'],vmin=-1000,vmax=0,cmap=plt.get_cmap('ocean'),rasterized=True)
ax.set_xlim([-123.59,-123.1])
ax.set_ylim([48.51,48.86])
plt.colorbar(pc,shrink=0.7,extend='both')

ax.set_xlabel('Longitude $[^o E]$')
ax.set_ylabel('Latitude $[^o N]$')
#fig.savefig('2016CruiseTrack.pdf')

def on_click(event):
    # get the x and y coords, flip y from top to bottom
    x, y = event.x, event.y
    if event.button == 1:
        if event.inaxes is not None:
            print('[%f, %f]' % (event.xdata, event.ydata))
plt.connect('button_press_event', on_click)

plt.show()
