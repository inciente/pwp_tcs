import xarray as xr; import pandas as pd; 
from datetime import datetime, timedelta; 
import numpy as np; import seawater, gsw; 
from dataclasses import dataclass


''' 
This module is largely refactored from Earlew's pwp_python_00 package. Here, PWP functions have been broken into multiple subfunctions, and relevant classes were to improve organization and modularity with other applications, such as coupling with atmospheric models.

@earlew describes the algorithm of his pwp.run() as 

        1) Set model parameters (see set_params function in PWP_helper.py). 
        2) Read in forcing and initial profile data.
        3) Prepare forcing and profile data for model run (see prep_data in PWP_helper.py).
            3.1) Interpolate forcing data to prescribed time increments.
            3.2) Interpolate profile data to prescribed depth increments.
            3.3) Initialize model output variables.
        4) Iterate the PWP model specified time interval:
            4.1) apply heat and salt fluxes
            4.2) rotate, adjust to wind, rotate
            4.3) apply bulk Richardson number mixing
            4.4) apply gradient Richardson number mixing
            4.5) apply drag associated with internal wave dissipation
            4.5) apply diapycnal diffusion       
        5) Save results to output file

'''

# So many steps use the series of profiles t, s, d, u, v (d is density). This can be simplified to an xr.Dataset() that includes all. WIll learn more about implementation to figure out requirements.

class World: 
    ''' 
    This class stores PWP simulation parameters and some major methods involved in timestepping.
    Methods described here become available to any other functions once parameters are set. 
    All parameters are set to default values unless user decides otherwise. 
    This class exists to ensure consistency between steps and experiments.
    '''

    # parameters necessary for grids and vectors
    lat : float; 
    dz : float = 0.5; # in meters
    dt : float = 300; # in seconds
    zmax : float = 150; # in meters
    # parameters involved in mixing
    mld_thres : float = 1e-4; # density threshold for MLD
    Ri_g : float = 0.25; # threshold gradient richardson number
    Ri_b : float = 0.65; # threshold bulk richardson number
    rkz : float = 1e-6; # backgroudn vertical diffusivity
    # thermodynamics and optics
    beta_1 : float = 0.6; # longwave extinction coefficient (meters)
    beta_2 : float = 20; # shortwave extinction coefficient
    cpw : float = 4183; # heat capacity of seawater
    drag_coef : float = 0.1; # internal wave drag, to be multiplied by f
    # set flags to decouple aspects of air_sea interaction
    wind_ON : bool = True; # momentum fluxes
    heat_ON : bool = True; # heat fluxes
    drag_ON : bool = True; # rayleigh friction on flow
    emp_ON : bool = True; # evap minus precip freshwater fluxes


    def f( self ):
        return 4 * np.pi * np.sin( self.lat / 180 * np.pi ) / 24 / 3600

    def make_absorption( self, profile ):
        rs1 = 0.6; # fraction of heat contained in shortwave
        ztop = profile['z'] - self.dz / 2; # depth of top of cells
        zbot = profile['z'] + self.dz / 2; # depth of bottom of cells
        # treat shortwave and longwave separately
        absorb_sw = rs1 * ( np.exp( - ztop / self.beta_1 ) - np.exp( - zbot / self.beta_1 ) ) 
        absorb_lw = ( 1 - rs1 ) * ( np.exp( -ztop / self.beta_2 ) - np.exp( zbot / self.beta_2 ) );
        # save total absorbption
        tot_absorb = xr.DataArray( data = abosrb_sw + absorb_lw , 
                            coords = { 'z' : profile['z'] } )
        return tot_absorb


    def prepare_profile( self , profile ):
        # Profile is xr.dataset with temp, sal, dens, u, v
        zvals = np.arange( self.dz / 2 , self.zmax + 0.0001 , self.dz ); # depth of cell centers
        profile = profile.interp( z = zvals );
        # add absorption profile
        profile['absorb'] = self.make_absorption( profile )
        return profile

    def flag_forcing( self, forcing ):
        ''' 
        Apply forcing flags. This might fit better in a separate class that handles forcing exclusive and 
        that must be instantiated in relation to an instance of World. TBD as project evolves.
        '''
        if ~ wind_ON :
            forcing['taux'] = 0; 
            forcing['tauy'] = 0; 
        if ~ heat_ON :
            forcing['q_in'] = 0; 
            forcing['q_out'] = 0;
        if ~ emp_ON : 
            forcing['emp'] = 0; 
        return forcing

    def rotate( self, profile ):
        # Apply rotation that results from Coriolis over dt / 2
        f = self.f(); 
        profile['u'] = profile['u'] + f * profile['v'] * self.dt / 2
        profile['v'] = profile['v'] - f * profile['u'] * self.dt / 2 
        return profile 

    def find_MLD( self, profile ):
        # find ml depth of profile and its index in z
        rho_diff = profile['dens'] - profile['dens'].isel( z = 0 ); 
        # find first index for which diff is above threshold
        mld_idx = np.nonflatzero( rho_diff > self.mld_thres )[0]
        assert mld_idx.size != 0, 'Error: ML depth is undefined'
        # get numerical value of MLD
        mld = profile['z'].isel( z = mld_idz )
        return mld, mld_idx

    def wind_on_ML( self, profile , forcing ):
        # apply momentum forcing to mixed layer
        mld, mld_idx = self.find_MLD( profile );
        mass = mld * profile['dens'][0]; 
        dU = forcing['taux'] / mass * self.dt 
        dV = forcing['tauy'] / mass * self.dt
        # ------ now add velocities to old profile
        profile['u'][ : mld_idx ] = profile['u'][ : mld_idx ] + dU
        profile['v'][ : mld_idx ] = profile['v'][ : mld_idx ] + dV;        
        return profile

    def update_surface( self, profile , forcing ):
        # Apply heat and salinity fluxes to uppermost level of a profile
        dT = ( forcing['q_in'] * profile['absorb'].isel( z = 0 )  + forcing['q_out'] ) \
                 * self.dt / ( self.dz * profile['dens'].isel( z = 0 ) * self.cpw )
        dS = forcing['emp'] * profile['sal'].isel( z = 0 ) * self.dt / self.dz
        # set new values
        profile['temp'][0] += dT
        profile['sal'][0] += dS
        return profile

    def subsurface_sw( self, profile, forcing ):
        # apply downwelling sw radiation to subsurface layers
        dT = ( forcing['q_in'] * profile['absorb'][1:] ) * self.dt 
        dT = dT / ( self.dz * self.cpw * profile['dens'][1:] )
        profile['temp'][1:] = profile['temp'][1:] + dT
        return profile 


def pwp_step( world, profile, forcing ):
    # Apply PWP algorithm 
    # Heat and salinity fluxes at the surface and then below
    profile = world.update_surface( profile, forcing );
    profile = world.subsurface_sw( profile, forcing );          
    # --- might want to add a point checking for freezing here
    profile['dens'] = sw.dens0( profile['sal'], profile['temp'] ); # update density
    # --- relieve static instability
    # apply momentum flux and coriolis rotation
    profile = world.rotate( profile );
    profile = world.wind_on_ML( profile, forcing ); 
    profile = world.rayleigh_friction( profile );
    profile = world.rotate( profile ); # coriolis for dt / 2
    # time to apply mixing parameterizations
    profile = world.bulk_mix( profile );


def stir( profile, rc, r, j):
    # This mixes cells around j (depth-wise) to ensure that the 
    # Ri after (r_new) is r_new > rc (stands for critical Ri)
    pass 

def interpolator( xr_obj, x_locs ):
    # interpolate xr_obj onto locations specified by x_locs (dict)
    xr_obj = xr_obj.interp( x_locs )
    return xr_obj

def update_uv( ocn , met, params ):
    # This is the function that takes care of changing ocean u,v
    # Will return a new ocn object
    pass 

def update_TS( ocn , met, params ):
    # Function to implement a single timestep on ocn T,S
    # Will return a new ocn object
    pass 

def SW_radiation( ocn, met, params ):
    # Given shortwave flux, determine exponential profile of absorption. Might not need to know ocn state. 
    pass 

def bulk_Ri( ocn , params ):
    # Compute bulk Richardson Number
    pass

def get_MLD( profile , params ):
    # Return ML depth given profile and parameters.
    pass 

# And many other functions that the existing functions in Earlew's code will be broken down into. 
