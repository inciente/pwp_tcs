'''
Module meant to streamline loading of the Roemmich and Gilson 
ARGO climatology. Picks up pieces of code from multiple places.
'''
import numpy as np; import xarray as xr; 
import gsw, sys; import pandas as pd; 
from dateutil.relativedelta import relativedelta
from datetime import datetime, timedelta

sys.path.append( '/home/noelgb/repositories/EastPac/' )
import intercomparison

#sys.path.append( modview_path )
#import loader, timetools, mapper
#sys.path.append( '/home/noel/Documents/PISTON/pyscripts')
#import heat_funcs

argo_RG_path = '/home/noelgb/data/ARGO/RG_ArgoClim_Merged.json'

monthly_paths = []; # paths to monthly updates
tmean = 'ARGO_TEMPERATURE_MEAN';
tanom = 'ARGO_TEMPERATURE_ANOMALY';
smean = 'ARGO_SALINITY_MEAN';
sanom = 'ARGO_SALINITY_ANOMALY'; # variable names in files

def to_real_time( xr_obj ):
    # Change time from months since 2004 Jan 15
    basetime = datetime( 2004, 1, 15 ); 
    new_time = [ basetime + relativedelta( months =+ \
            np.floor( nmonths) ) for nmonths in xr_obj['TIME'].values]
    xr_obj['TIME'] = xr.DataArray( data = new_time, dims=('TIME') )
    #xr_obj = xr_obj.swap_dims( {'TIME':'date'} )
    return xr_obj

def make_single( xr_obj, annual = False ):
    if annual:
        tan = 'ARGO_TEMPERATURE_ANNUAL_ANOMALY';
        san = 'ARGO_SALINITY_ANNUAL_ANOMALY';
    else:
        tan = tanom; san = sanom;
    gnames = ['TEMP','SALT'];
    cnames = [[tmean, tan], [smean, san]]
    # Take anomalies and means and add them together
    for jj in [0,1]:
        set_as = gnames[jj]
        get_from = cnames[jj]
        # Check what variables are present in xr_obj to add
        if get_from[0] in xr_obj.keys():
            xr_obj[set_as] = xr_obj[get_from[0]] + xr_obj[get_from[1]]
            xr_obj = xr_obj.drop_vars( get_from )
    return xr_obj 


def make_single2( ds ):
    ds[tmean] = ds[tmean].isel( TIME = 0 )
    ds[smean] = ds[smean].isel( TIME = 0 )
    ds['TEMP'] = ds[tmean] + ds['ARGO_TEMPERATURE_ANOMALY']
    ds['SALT'] = ds[smean] + ds['ARGO_SALINITY_ANOMALY']
    ds.drop_vars( ['ARGO_TEMPERATURE_ANOMALY','ARGO_SALINITY_ANOMALY'] )
    return ds 

def load_all( ):
    # Load entire RG climatology including monthly updates
    ds = intercomparison.json_to_xr( argo_RG_path , decode_times = False )
    ds = make_single2( ds );
    ds = to_real_time( ds )
    return ds

def area_slicer( xr_obj , xslice, yslice ):
    return xr_obj.sel( {'LONGITUDE':xslice, 'LATITUDE':yslice } )
