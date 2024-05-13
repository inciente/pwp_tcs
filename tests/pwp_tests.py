import unittest

import xarray as xr; import pandas as pd; 
import numpy as np; 
import os

'''
Module to test features of pwp code.
'''

# Load src components of pwp package
os.path.abspath( os.path.join( os.path.dirname( __file__ ) , 
                 '..', 'src' ) )
import pwp


# Load an ARGO profile to use while running tests
profile_src = 'data/float-6902855-cyc-87-Oct-2019.nc'
test_profile = xr.open_dataset( profile_src )

# Create instance of world to call methods
test_world = pwp.World( lat = 20 )








