import xarray as xr; import numpy as np; import pandas as pd; 
import gsw

from argopy import DataFetcher
from argopy import ArgoIndex

# Pre-load to speed up future calls
#ArgoIndex.load()

'''
Tools to fetch ARGO profiles and test for thermal inversions within them. 
'''

class argo_search:

    def __init__( self, region, plims, dates ):
        # use argopy to find argo profiles within specified parameters
        # Input: 
        # ----- region is a list [xmin, xmax, ymin, ymax]
        # ----- plims is a list [pmin, pmax]
        # ----- dates is a list of two dates as strings 

        self.region = region;
        self.plims = plims;
        self.dates = dates; 

        # Use parameters to find individual argo profiles
        self.ds, self.directory = search_region( self.region + self.plims , dates )
        self.ds = self.ds.to_xarray()

    def single_profile( self, N ):
        # Search self.ds and find data for profile in row N of self.directory
        # get profile metadata
        cyc = self.directory['cyc'][N]; wmo = self.directory['wmo'][N]
        
        # search self.ds for matching entries
        nancheck = self.ds['PLATFORM_NUMBER'] == wmo;
        nancheck = nancheck * ( self.ds['CYCLE_NUMBER'] == cyc )

        # isolate results
        prof = self.ds.isel( N_POINTS = nancheck ).swap_dims( {'N_POINTS':'PRES'} )
        # add CT and RHO to profile
        prof = prepare_profile( prof ); 
        return prof


# To be used with argopy.DataFetcher.region()
def search_region( region, period ):
    # region is a list of 4 coordinates [xmin, xmax, ymin, ymax, pmin, pmax ]
    # period is a list of two strings w/ format ['2011-01-01','2011-04-01']
    
    # instantiate finder
    fetched = DataFetcher( mode = 'research' , parallel = True , 
                           progress = True ).region( region + period  )
    # make dataframe of individual profiles found
    info = fetched.to_index(); # use wmo and cyc columns to explore dataset
    return fetched, info


# and the profiles that are stored within the output

def prepare_profile( xr_prof ):
    # Take in a profile,  pass through GSW and return gridded CT, dens
    
    # data into xarray with easier format
    #xr_prof = profile.to_xarray(); # easier format
    #xr_prof = xr_prof.swap_dims( { 'N_POINTS':'PRES' } )
        
    # save conservative temperature and sigma0
    xr_prof['CT'] = gsw.CT_from_t( xr_prof['PSAL'], xr_prof['TEMP'], 
                                   xr_prof['PRES'] )
    xr_prof['RHO'] = gsw.density.rho( xr_prof['PSAL'], xr_prof['CT'], 
                                   xr_prof['PRES'] )
    # Many variables have values repeated for all pressures
    # clean those up
    xr_prof = clean_argo_prof( xr_prof );    

    return xr_prof 

def clean_argo_prof( xr_prof ):
    # Many variables have values repeated for all pressures.
    # Instead of keeping repeated values along PRES, create
    # store somewhere else
    vars2clean = ['LATITUDE','LONGITUDE','TIME','CYCLE_NUMBER',
        'DIRECTION','PLATFORM_NUMBER','PRES_ERROR', 'PSAL_ERROR',
        'TEMP_ERROR']
    for var in vars2clean: 
        pass
    return xr_prof

















