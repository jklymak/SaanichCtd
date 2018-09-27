import logging
logging.basicConfig(level=logging.INFO)
import CtdHexTrans as ctdtrans
from datetime import datetime, timedelta
import scipy.io as sio
import xarray as xr
import numpy as np
import seawater.eos80 as sw
import seawater
import glob as glob
import os.path
import subprocess
from jmkdata import bindata1d
from timeit import default_timer as timer
import getInletX
import sys

todo = sys.argv[1]

coffset = 1.7

lat0 = 48.7
lon0 = -123.5
points = np.array([[-123.241147, 48.720417],
                    [-123.338702, 48.766345],
                    [-123.418352, 48.722311],
                    [-123.469599, 48.705739],
                    [-123.5, 48.700530],
                    [-123.5            , 48.6418],
                    [-123.5, 48.5]])

innames = glob.glob(todo+'/*.hex')
cgrid = dict()

for name in innames:
    basename = os.path.splitext(os.path.basename(name))[0]
    dirname = os.path.dirname(name)
    logging.info(dirname)
    ncname = dirname+'/'+basename+'.nc'
    if not os.path.exists(ncname) or (os.path.getctime(name) >
            os.path.getctime(ncname)):
        logging.info(name)
        ctd = ctdtrans.CtdHex2mat(name)

        # parse some header stuff
        for l in ctd['headinfo'].splitlines():
            try:

                if (l.find('** Lon') == 0):
                    lon, minute = (l.split(':')[1]).split(' ')
                    ctd['lon'] = -float(lon) - float(minute)/60.
                if (l.find('** Lat') == 0):
                    lat, minute = (l.split(':')[1]).split(' ')
                    ctd['lat'] = float(lat) + float(minute)/60.
                if (l.find('** Sta') == 0):
                    ctd['id'] = l.split(':')[1]
                if (l.find('* cast') == 0 ):
                    timest = l[11:28]
                    ctd['time'] = datetime.strptime(timest, '%d %b %Y %H:%M')
                    ctd['matlabtime'] = ctdtrans.datetime2matlab(ctd['time'])
            except:
                logging.info(basename)
                pass

        #
        # get salt and pden...
        ctd['c0'] = ctd['c']
        N = len(ctd['c'])
        ctd['coffset'] = coffset
        ctd['c'] = np.interp(np.arange(N)-coffset, np.arange(N), ctd['c0'])
        ctd['sal'] = sw.salt(ctd['c'] * 10. / seawater.constants.c3515,
                ctd['t'], ctd['p'])
        ctd['pden']=sw.dens(ctd['sal'], ctd['t'], ctd['p'])

        ctd = ctdtrans.getO2(ctd,'v1');
        ctd = ctdtrans.getFlu(ctd,'v2');
        ctd = ctdtrans.getPar(ctd,'v3');

        alongx, acrossx = getInletX.getInletX(ctd['lon'], ctd['lat'])

        ctd['alongx'] = alongx
        ctd['acrossx'] = acrossx

        ds = xr.Dataset(
            {
            'cond': (['scan'], ctd['c'], {'units':'S/m', 'offset': '{} scans'.format(coffset)}),
            'cond0': (['scan'], ctd['c0'], {'units':'S/m'}),
            'temp': (['scan'],ctd['t'], {'units':'deg C'}),
            'pres': (['scan'],ctd['p'], {'units':'dbar'}),
            'O2':  (['scan'],ctd['O2'], {'units':'mmol/kg'}),
            'O2sat':  (['scan'],ctd['O2sat'], {'units':'percent saturation'}),
            'Par':  (['scan'],ctd['Par'], {'units':'E/m^2'}),
            'Flu':  (['scan'],ctd['Flu'], {'units':'Flu'}),
            'sal':  (['scan'],ctd['sal'], {'units':'psu'}),
            'pden':  (['scan'],ctd['pden'], {'units':'kg/m^2'}),
            'time': ctd['time'],
            },
            coords={'scan': (['scan'], np.arange(len(ctd['c'])))},
            attrs={
            'alongx':ctd['alongx'],
            'acrossx':ctd['acrossx'],
            'lat':ctd['lat'],
            'lon':ctd['lon'],
            'id':ctd['id'],
            'header':ctd['headinfo']}
        )
        logging.debug('ds:')
        logging.debug(ds)

        ds.to_netcdf(ncname)
        logging.info('saved nc')

        ctdout = ctd.copy()
        ctdout.pop('time')
        ctdout.pop('matlabtime')
        ctdout['den'] = ctdout['pden']
        ctdout.pop('pden')
        ctdout['time'] = ctd['matlabtime']

        sio.savemat(dirname+'/'+basename+'New.mat', ctdout, format='5')
        subprocess.call(['octave', 'saveasstruct.m', dirname+'/'+basename+'New.mat', 'ctd'])

        logging.info(ctd['lon'])
        logging.info(ctd['lat'])
        logging.info(ctd['time'])
        logging.info(ctd['matlabtime'])

        logging.debug(cgrid.keys())

# make cgrid

innames = glob.glob(todo+'/20*.nc')
zbins = np.arange(325)
cgrid = dict()
cgrid['depths'] = zbins[1:] - 0.5
cgrid['name'] = []
for name in innames:
    print(f'Doing {name}')
    ds = xr.open_dataset(name)
    logging.debug(ds)
    logging.info('Starting %s', name)
    # find downcast...
    ind = np.where(ds.pres > 20)[0][0]
    start = ind
    stop = ind
    p = ds.pres.values
    while (p[start - 10] < p[start]) and (start > 10):
        start = start - 10

    while ((p[stop + 10] > p[stop])
            and (stop + 10 < len(p))):
        stop = stop + 10
    # OK stop will have to be deepest
    stop = np.argmax(p)
    inds = range(start, stop)
    for tobin in ds.keys():
        logging.debug('Starting %s', tobin)
        start = timer()
        if ds[tobin].dims == ('scan',):
            dat, vvv, vv = bindata1d(zbins,
                    ds.pres.values[inds], ds[tobin].values[inds])
        else:
            dat = ds[tobin]
        end = timer()
        logging.debug('Done %s %1.2f', tobin, end-start)
        start = timer()
        if tobin in cgrid:
            logging.debug('append variable %s', tobin)
            if len(np.shape(dat)):
                ax = 1
                cgrid[tobin] = np.append(cgrid[tobin],
                        dat[:, np.newaxis], axis=ax)
            else:
                ax = 0
                dat = np.array([np.array(dat)])
                cgrid[tobin] = np.append(cgrid[tobin],
                                    dat, axis=ax)

        else:
            if len(np.shape(dat)):
                cgrid[tobin] = dat[:, np.newaxis]
            else:
                cgrid[tobin] = np.array([np.array(dat)])
        stop = timer()
        logging.debug('Concat time %1.2f', end-start)
    for towrite in ds.attrs.keys():
        logging.debug('Starting %s', towrite)
        if towrite in cgrid:
            cgrid[towrite] = np.append(cgrid[towrite], ds.attrs[towrite])
        else:
            cgrid[towrite] = ds.attrs[towrite]
    cgrid['name'] += [os.path.splitext(os.path.basename(name))[0]]
    logging.info('Done %s', name)

# sort by alongx:

print(cgrid)

ind = np.argsort(cgrid['alongx'])
for key in cgrid.keys():
    print(key)
    if isinstance(cgrid[key], np.ndarray):
        if len(cgrid[key].shape) == 2:
            cgrid[key] = cgrid[key][:, ind]
        elif cgrid[key].shape[0] == len(ind):
            cgrid[key] = cgrid[key][ind]
    else:
        if len(ind) > 1:
            cgrid[key] = np.asarray([cgrid[key][i] for i in ind])
        else:
            cgrid[key] = np.asarray([cgrid[key]])


logging.debug(cgrid['alongx'])

# save to netcdf...
ctd = cgrid.copy()
ds = xr.Dataset(
    {
    'cond': (['depths', 'time'], ctd['cond'], {'units':'S/m', 'offset': '{} scans'.format(coffset)}),
    'cond0': (['depths', 'time'], ctd['cond0'], {'units':'S/m'}),
    'temp': (['depths', 'time'],ctd['temp'], {'units':'deg C'}),
    'pres': (['depths', 'time'],ctd['pres'], {'units':'dbar'}),
    'O2':  (['depths', 'time'],ctd['O2'], {'units':'mmol/kg'}),
    'O2sat':  (['depths', 'time'],ctd['O2sat'], {'units':'percent saturation'}),
    'Par':  (['depths', 'time'],ctd['Par'], {'units':'E/m^2'}),
    'Flu':  (['depths', 'time'],ctd['Flu'], {'units':'Flu'}),
    'sal':  (['depths', 'time'],ctd['sal'], {'units':'psu'}),
    'pden':  (['depths', 'time'],ctd['pden'], {'units':'kg/m^2'}),
    'lat': (['time'], ctd['lat'], {'units':'deg N'}),
    'lon': (['time'], ctd['lon'], {'units':'deg W'}),
    'alongx': (['time'], ctd['alongx'], {'units':'dist from S4 [km]'}),
    'acrossx': (['time'], ctd['acrossx'], {'units':'dist from S4 [km]'}),
    'id': (['time'], ctd['id'])
    },
    coords={'depths': (['depths'], ctd['depths']),
            'time': (['time'], ctd['time'])},
    attrs={}
)
ds.to_netcdf(dirname+'/CtdGrid.nc', 'w')
logging.debug(ds)
logging.info('Saved %s', dirname+'/CtdGrid.nc')

# save cgrid to matfile
cgridout = cgrid.copy()
cgridout['time'] = np.array([ctdtrans.datetime2matlab(
                    datetime.utcfromtimestamp(t.tolist()/1.e9))
        for t in cgrid['time']])
cgridout['den'] = cgridout['pden'].data
cgridout['t'] = cgridout['temp'].data
cgridout['c'] = cgridout['cond'].data
cgridout['sal'] = cgridout['sal'].data
cgridout['c0'] = cgridout['cond0'].data
cgridout['flu'] = cgridout['Flu'].data
cgridout['p'] = cgridout['pres'].data
cgridout['O2'] = cgridout['O2'].data
cgridout['O2sat'] = cgridout['O2sat'].data
cgridout['par'] = cgridout['Par'].data
for topop in ('pden', 'temp', 'cond', 'cond0', 'Flu'):
    cgridout.pop(topop)
sio.savemat(dirname+'/CtdGridNew.mat', cgridout, format='5')
logging.info('Saved %s', dirname+'/CtdGridNew.mat')
ans = subprocess.call(['octave', 'saveasstruct.m', dirname+'/CtdGridNew.mat', 'cgrid'])
print(ans)
logging.info('Saved %s structured matlab', dirname+'/CtdGrid.mat')
