import xarray as xr; import numpy as np; import pandas as pd; 

'''
Module containing supporting functions to estimate metrics of ocean-atmosphere
interactions under tropical cyclones.
'''

def get_layer_thickness( pressure , z0 = 0 ):
    # Take in xr_obj with pressure coordinates.
    # Assuming that data goes up to z0, and that pressure levels are given at the 
    # center of each cell, return the thickness of each cell
    
    dz = [None] * len( pressure ); # thickness of each cell
    bottoms = [None] * len( pressure ); # bottom of each cell
    tops = [None] * ( len( pressure ) ); 
    tops[0] = z0; 

    # Now proceed iteratively
    for jj in range( len( pressure ) ):
        if jj >= 1:
            tops[jj] = tops[jj-1] + dz[jj-1]; 
        dz[jj] = ( pressure[jj] - tops[jj] ) * 2 ; # twice distance from center to top
        bottoms[jj] = tops[jj] + dz[jj] 

    # Store everything as xr dataset
    zgrid = xr.Dataset(  )
    zgrid['dz'] = xr.DataArray( data = np.array( dz ) , dims = pressure.dims , coords = pressure.coords ); 
    zgrid['tops'] = xr.DataArray( data = np.array( tops ), dims = pressure.dims, coords = pressure.coords ); 
    zgrid['bottoms'] = xr.DataArray( data = np.array( bottoms ) , dims = pressure.dims , coords = pressure.coords )
    
    return zgrid

def make_profile_conscious( profile, vert_coord = 'PRESSURE' ):
    # Teach a vertical profile about its own vertical grid
    zgrid = get_layer_thickness( profile[ vert_coord ] )

    profile['dz'] = zgrid['dz']; profile['tops'] = zgrid['tops']
    profile['bottoms'] = zgrid['bottoms']
    return profile


def mix_profile( profile, zmix , vert_coord = 'PRESSURE' ):
    # Take in a temperature with ['TEMP'], ['SALT'], ['dz'], ['bottoms'], and ['tops'] 
    # and mix its properties down to depth zmix. 
    # Return transformed profile, and assume no mixing below zmix.

    # Find all cells fully contained between the surface and zmix
    fully_contained = profile['bottoms'] < zmix
    final_cell = ( profile['tops'] < zmix ) * ( profile['bottoms'] > zmix ); 

    # Weighted average temp and salinity in layers fully above zmix
    def weighted_mean( xr_obj, var ):
        full_part = xr_obj[ var ].isel( { vert_coord : fully_contained } ).weighted( xr_obj['dz'][fully_contained] \
                                            ).mean( vert_coord )
        full_part = full_part * ( profile['dz'][fully_contained].sum() / zmix )
        return full_part 
    # Now we need a factor for the layer that includes zmix
    final_factor = ( zmix - profile['tops'].isel( { vert_coord : final_cell } ) ) / zmix
    
    Tmean = weighted_mean( profile, 'TEMP' ) + final_factor * profile['TEMP'].isel( { vert_coord : final_cell } ).values
    Smean = weighted_mean( profile, 'SALT' ) + final_factor * profile['SALT'].isel( { vert_coord : final_cell } ).values

    # Now update TEMP and SALT in all cells (fully and partially contained)
    nu_profile = xr.Dataset()
    nu_temp = profile['TEMP'].values.copy(); nu_salt = profile['SALT'].values.copy();
    
    nu_temp[ fully_contained.values ] = Tmean ; 
    nu_salt[ fully_contained.values ] = Smean ;

    nu_temp[ final_cell.values ] = profile['TEMP'][final_cell].values*(1-final_factor) + final_factor * Tmean;
    nu_salt[ final_cell.values ] = profile['SALT'][final_cell].values*(1-final_factor) + final_factor * Smean;

    nu_profile['TEMP'] = xr.DataArray( data = nu_temp, dims = profile.dims , coords = profile.coords ); 
    nu_profile['SALT'] = xr.DataArray( data = nu_salt, dims = profile.dims, coords = profile.coords )
    
    return nu_profile


def PE_change( prof_before, prof_after, vert_coord = 'PRESSURE' ):
    # Take in two profiles, compute their PEs, and return difference (integer)
    for prof in [prof_before, prof_after]:
        prof['DENS'] = gsw.rho( prof['SALT'], prof['TEMP'], prof[ vert_coord ] )
        prof = prof.sel( { vert_coord : slice(0, 300) } )

    def get_PE( profile ):
        pe_dens = ( profile['DENS'] * 9.81 * ( profile[ vert_coord ] + 300 ) ).integrate( vert_coord )
        return pe_dens
                                                                            
    diff = get_PE( prof_after ) - get_PE( prof_before )
    return diff.values

  
