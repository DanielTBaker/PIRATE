

import numpy as np
from astropy.io import fits
import astropy.units as u


def axis_array(header,axID):
    try:
        nAx = header['NAXIS'+str(axID)]
    except:
        nAx = header['NAXIS']
    try:
        dAx = header['CDELT'+str(axID)]
    except:
        dAx = header['CDELT']            
    try:
        refVal = header['CRVAL'+str(axID)]
    except:
        refVal = header['CRVAL']            
    try:
        refID = int(header['CRPIX'+str(axID)])-1
    except:
        refID = int(header['CRPIX'])-1            
    ax = np.linspace(0,nAx-1,nAx)*dAx
    ax += refVal - ax[refID]
    return(ax)

def dspec_from_fits(fname):
    hdul = fits.open(fname)
    dspec = hdul['PRIMARY'].data.T
    try:
        time = axis_array(hdul['PRIMARY'].header,2)*u.s
        freq = axis_array(hdul['PRIMARY'].header,1)*u.MHz
    except Exception as e:
        print(e)
        print(fname)
        print(hdul['PRIMARY'].header)
        return None
    hdul.close()
    return(dspec,time,freq)
   
    