import xarray as xr; import pandas as pd; 
from datetime import datetime, timedelta; 
import numpy as np; import seawater, gsw; 

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
