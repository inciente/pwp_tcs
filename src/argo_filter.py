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
 
    def distribute_task( self, func ):
        # Applies func to all profiles in the search class
        N = len( self.directory ); 
        print('Will apply function to ' + str(N) + ' profiles.')
        
        # Now access individual profiles and store output in list
        results = []; # make sure that func preserves some metadata
        failed = []; # store rows where execution fails 

        for jj in range( N ):
            # Extract data
            profile = self.single_profile( jj )
            
            # Argopy loads some broken data, so func() will not always run successfully
            try: 
                # Run func and save output only if useful
                res = func( profile );
                if res is not None:
                    # useful output. simplify drops redundant info
                    results.append( simplify_profile( res ) ) 
            except:
                # Store position of failure
                failed.append( jj )

        print( 'Execution failed on ' + str( len( failed ) ) + ' profiles.' )
        return results, failed 


def simplify_profile( profile ):
    # Reduce the data stored in xr.Datasets of argo profiles that come out of argopy

    # List of variables that will be switched to nprof 
    vars2clean = ['LATITUDE','LONGITUDE','TIME','CYCLE_NUMBER',
        'DIRECTION','PLATFORM_NUMBER','PRES_ERROR', 'PSAL_ERROR',
        'TEMP_ERROR']

    for var in vars2clean:
        # Get variable as single-element 1-D array
        dat = np.array( [ profile[var][0].values ] ).flatten()
        # Save it to replace extended version
        profile[ var ] = xr.DataArray( data = dat, dims = ('nprof') )
    
    return profile.drop_vars( 'N_POINTS' )

 
# To be used with argopy.DataFetcher.region()
def search_region( region, period, opt_dict = None ):
    # region is a list of 4 coordinates [xmin, xmax, ymin, ymax, pmin, pmax ]
    # period is a list of two strings w/ format ['2011-01-01','2011-04-01']
    
    # instantiate finder
    if opt_dict is None:
        # This is where options for fetcher are set, enables BGC usage
        opt_dict = { 'mode' : 'research', 'parallel' : True, 'progress' : True }

    fetched = DataFetcher( opt_dict ).region( region + period  )
    # make dataframe of individual profiles found
    info = fetched.to_index(); # use wmo and cyc columns to explore dataset
    return fetched, info
















