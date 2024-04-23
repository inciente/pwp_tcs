import numpy as np; import xarray as xr; import pandas as pd; 
import sys, gsw

''' 
Functions to analyze characteristics of inversions in argo profiles.
Storing them within a class is useful because otherwise we'd be calculating dz and heaviside again and again. 
'''

class inversion:

    def __init__( self, profile, world = None ):
        # Store profile and save basic quantities needed to describe it.
        self.profile = profile;
        self.sst = profile['temp'].isel( z = 0 ); 
        # Boolean mask to mark the inversion itself
        self.where = heaviside( profile['temp'] - self.sst , delta = 0.05 )
        # dz, preferrably is constant
        self.dz = np.gradient( profile['z'].values )

    def thickness( self ):
        return ( self.dz * self.where ).sum( 'z' )

    def heat_content( self ):
        # Total heat above SST available within inversion
        rel_temp = ( self.profile['temp'] - self.sst ).where( self.where )
        # multiply relative temp by rho, cp, and integrate vertically
        heat_inv = 1024 * 4000 * ( self.dz * rel_temp ).sum( 'z' )
        return heat_inv

    def max_temp( self ):
        # Maximum temperature (relative to sst) within inversion
        rel_temp = ( self.profile['temp'] - self.sst ).where( self.where )
        return rel_temp.max( 'z' )

    def nsquared( self ): 
        # Return mean value of N2 inside the inversion
        N2 = 9.81 / 1024 * self.profile['dens'].differentiate( 'z' )
        # Save mean value within inversion
        N2 = N2.where( self.where ).mean( 'z' )
        return N2

def heaviside( xr_obj, delta = 0 ):
    # 0 when dat < 0, 1 when dat >= 0;
    nu_obj = np.zeros( xr_obj.shape );
    nu_obj[ xr_obj >= delta ] = 1; 
     
    # reformat as xr with same coords
    nu_obj = xr.DataArray( data = nu_obj , coords = xr_obj.coords )
    return nu_obj 



