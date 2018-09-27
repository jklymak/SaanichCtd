import numpy as np

def dist(x1,y1, x2,y2, x3,y3): # x3,y3 is the point
    px = x2-x1
    py = y2-y1

    dsq = px*px + py*py

    u =  ((x3 - x1) * px + (y3 - y1) * py) / dsq

    if u > 1:
        u = 1
    elif u < 0:
        u = 0

    x = x1 + u * px
    y = y1 + u * py

    dx = x - x3
    dy = y - y3

    dist = np.sqrt(dx*dx + dy*dy)

    return dist

def getInletX(lon,lat):
    dat=[[-123.534430, 48.523672],
    [-123.543075, 48.536432],
    [-123.544309, 48.551016],
    [-123.509733, 48.560130],
    [-123.497384, 48.569245],
    [-123.496149, 48.585651],
    [-123.499854, 48.639427],
    [-123.498619, 48.679531],
    [-123.496149, 48.695937],
    [-123.456633, 48.714167],
    [-123.420822, 48.726927],
    [-123.393654, 48.743333],
    [-123.370192, 48.757005],
    [-123.344259, 48.767031],
    [-123.326971, 48.765208],
    [-123.286220, 48.738776],
    [-123.259052, 48.718724],
    [-123.245469, 48.712344],
    [-123.261522, 48.693203],
    [-123.223241, 48.584740]]

    lon = np.atleast_1d(lon)
    lat = np.atleast_1d(lat)
    lon0=-123.5;
    lat0=48.65;

    mpernm = 1852;


    xx = np.zeros(len(dat))+1j*np.zeros(len(dat))
    latt = 0.*np.real(xx)
    lonn = 0.*np.real(xx)

    for n,d in enumerate(dat):
        xx[n]=(d[0]-lon0)*np.cos(lat0*np.pi/180)*mpernm*60.+1j*(d[1]-lat0)*mpernm*60.
        latt[n]=d[1]
        lonn[n]=d[0]
    alongx = np.real(xx)*0.
    alongx[1:]=np.real(np.cumsum(np.abs(np.diff(xx))))
    alongxnew = np.linspace(0.,60.e3,10000)
    latnew = np.interp(alongxnew,alongx,latt)
    lonnew = np.interp(alongxnew,alongx,lonn)

    xr = np.interp(alongxnew,alongx,np.real(xx))
    xi = np.interp(alongxnew,alongx,np.imag(xx))

    xnew = xr+1j*xi

    x0 = 15457.545754575458/1e3  # this is S4...
    # OK, have new high def version.  Now get lat and lon as x
    x=(lon-lon0)*np.cos(lat0*np.pi/180)*mpernm*60.+1j*(lat-lat0)*mpernm*60.
    xalong = 0*np.real(x)
    xacross = 0*np.real(x)
    for ind,xx in enumerate(x):
        closest = np.argmin(np.abs(xx - xnew))
        xalong[ind]  = alongxnew[closest]/1.e3 - x0
        xacross[ind] = np.abs(xx - xnew[closest])/1.e3

    return xalong, xacross

def getInletXNew(linelons, linelats, lon, lat, lon0, lat0, anchorind=0):

    '''
    get potison along an inlet line   Anchor Ind should be S4
    '''

    xx = (linelons - lon0)*np.cos(lat0 * np.pi / 180.) * 60. * 1.85  # km
    yy = (linelats - lat0)*60.*1.85     # km

    distline = np.cumsum(np.sqrt(np.diff(xx)**2 + np.diff(yy)**2))
    distline = np.append([0.], distline)

    x = (np.atleast_1d(lon) - lon0)*np.cos(lat0 * np.pi / 180.) * 60. * 1.85  # km
    y = (np.atleast_1d(lat) - lat0)*60.*1.85     # km

    ind = np.zeros(len(x))
    for j in range(len(x)):
        thedist = np.Inf
        for i in range(len(xx)-1):
            dd = dist(xx[i], yy[i], xx[i+1], yy[i+1], x[j], y[j])
            if dd < thedist:
                thedist = dd
                ind[j] = i
    alongx = x * 0.
    for i in range(len(x)):
        # get the distance along the line....
        indd = int(ind[i])
        x0 = xx[indd]
        x1 = xx[indd+1]
        y0 = yy[indd]
        y1 = yy[indd+1]

        xp = x[i] - x0
        yp = y[i] - y0

        dot = xp * (x1 - x0) + yp * (y1 - y0)
        dot = dot / np.sqrt((x1 - x0)**2 + (y1 - y0)**2)
        print(dot)
        alongx[i] = dot + distline[indd]

    alongx = distline[anchorind] - alongx

    return alongx
