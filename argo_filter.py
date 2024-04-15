import xarray as xr; import numpy as np; import pandas as pd; 
import gsw

from argopy import DataFetcher
from argopy import ArgoIndex

# Pre-load to speed up future calls
ArgoIndex.load()

'''
Tools to fetch ARGO profiles and test for thermal inversions within them. 
'''

# Define sample regions
WNP = [ 100, 160, 10, 35, 0, 250 ]

# To be used with argopy.DataFetcher.region()
def search_region( region, period ):
    # region is a list of 4 coordinates [xmin, xmax, ymin, ymax, pmin, pmax ]
    # period is a list of two strings w/ format ['2011-01-01','2011-04-01']
    
    # instantiate finder
    fetched = DataFetcher( mode = 'research' , parallel = True , progress = True ).region( region + period  )
    # make dataframe of individual profiles found
    info = fetched.to_index(); # use wmo and cyc columns to explore dataset
    return fetched, info


# and the profiles that are stored within the output

def prepare_profile( profile ):
    # Take in a profile,  pass through GSW and return gridded CT, dens
    
    # data into xarray with easier format
    profile_dat = profile.to_xarray(); # easier format
    profile_dat = profile_dat.swap_dims( { 'N_POINTS':'PRES' } )
        
    # save conservative temperature and sigma0
    xr_prof['CT'] = gsw.CT_from_t( xr_prof['PSAL'], xr_prof['TEMP'], xr_prof['PRES'] )
    xr_prof['RHO'] = gsw.density.rho( xr_prof['PSAL'], xr_prof['CT'], xr_prof['PRES'] )
    
    return xr_prof 


