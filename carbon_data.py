import xarray as xr; import pandas as pd; 
from datetime import datetime, timedelta
import numpy as np; import sys; 

'''
Tools to analyze variations in ocean carbon in response to tropical cyclones. 

'''

    
def numberfy( column ):
    # Turn column entries into floats. Use with dfs that are imported as str
    nucol = [ float( jj ) for jj in column ]
    return nucol

def decimal_to_date( decimal ):
    # Turn decimal dates into datetime
    # e.g. 1980.003 is January 2 of 1980
    dates = []
    for jj in range( len( decimal ) ):
        yrdate = decimal.loc[jj]
        date = datetime( int( np.floor( yrdate ) ), 1, 1, 0 ) \
            + timedelta( days = round( 365 * np.mod( yrdate , 1 ) ) ) 
        dates.append( date )
    return dates
        
class carbon_station:
    
    '''
    Interface to load and process observations of oceanic carbon from 
    BERM, BATS, and HAWI timeseries (Scripps CO2 program, R. Keeling)
    '''
    # Names of data columns in csv files. 
    # A dict format can be used to facilitate renaming.
    varnames = [' Decimal',' Depth', '    Temp', ' d13C-DIC',
                '  DIC', '   Salinity', '   ALK']            
    
    def __init__( self , fname ):
        # Find data and load with a reasonable format
        self.path = fname; 
        self.table = self.prepare_table()
        
    def prepare_table( self ):
        table = pd.read_csv( self.path , skiprows = 34 ) 
        table = table.iloc[1:].reset_index()
        # ensure that numerical variables are floats
        for var in self.varnames:
            table[var] = numberfy( table[var] )
        # Turn dates into datetime
        table[' Sample'] = decimal_to_date( table[' Decimal'] )
        table['Month'] = [ time.month for time in table[' Sample']]
        return table 
    
    def var_as_xr( self, var ):
        # Take var column and turn into xarray
        data = self.table[ var ];
        xr_obj = xr.DataArray( data = data , coords = \
                              {'time':self.table[' Sample']})
        return xr_obj 
    
    def seasonal_cycle( self, var ):
        # Take var column, turn into xarray
        xr_obj = self.var_as_xr( var )
        seas = xr_obj.groupby( 'time.month' ).mean()
        return seas
    
    def table_to_tz( self ):
        # Return an xr.dataset with time, depth dimensions
        zvec = np.sort( np.array( self.table[' Depth'].unique() ) )
        tvec = np.sort( np.array( self.table[' Sample'].unique() ) )
        # Iterate through rows and place data in correct spots
        nufields = ['Temp','DIC','d13C-DIC','Salinity','ALK']
        fieldsas = ['   Temp','    DIC', ' d13C-DIC' , '   Salinity', '    ALK']
        # Probably better to ask chat gpt help with this one. 
        # Might not be necessary. 
    
        
    

