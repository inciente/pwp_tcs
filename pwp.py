import xarray as xr; import pandas as pd; 
from datetime import datetime, timedelta; 
import numpy as np; import seawater, gsw; 

''' 
This module is largely refactored from Earlew's pwp_python_00 package. Here, PWP functions have been broken into multiple subfunctions, and relevant classes were to improve organization and modularity with other applications, such as coupling with atmospheric models.
'''

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

# And many other functions that the existing functions in Earlew's code will be broken down into. 
