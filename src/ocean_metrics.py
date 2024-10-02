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

def get_zmix( profile, dT, vert_coord = 'PRESSURE' ):
    # Gets conscious vertical profile and returns depth above which mean T is SST - dT
    # Temperature averaged between surface and bottom of cells
    noise = ( np.random.rand( len( profile['TEMP'] ) ) - 0.5 ) * 1e-3
    noise = np.expand_dims( noise, 1 )
    mean_at_bottom = ( ( profile['TEMP'] + noise ) \
                   * profile['dz'] ).cumsum( vert_coord ) \
                       / profile['dz'].cumsum( vert_coord )
    #mean_at_bottom = mean_at_bottom.where( profile['TEMP'].isel( { vert_coord : 0 } ) > 26 )
    target_temp = profile['TEMP'].isel( { vert_coord : 0 } ) - dT;  # sst - dT
    diff = ( profile['TEMP'] - target_temp ); # zmix is wherever this is 0 
    
    # Find temperature at bottom of cells surrounding zmix
    above_temp = mean_at_bottom.where( diff > 0 ).min( \
                             dim = vert_coord )
    below_temp = mean_at_bottom.where( diff < 0 ).max( \
                             dim = vert_coord )
    
    # Now find the depth of those cells 
    above_p = mean_at_bottom[ vert_coord ].where( \
           mean_at_bottom == above_temp ).mean( dim = vert_coord )
    below_p = mean_at_bottom[ vert_coord ].where( \
           mean_at_bottom == below_temp ).mean( dim = vert_coord )
    # Get dTdz to approximate zmix by taylor expansions
    dTdz = ( above_temp - below_temp ) / ( below_p - above_p ); # this should be positive
    zmix = above_p + ( above_temp - target_temp ) / dTdz # first order taylor
    
    return zmix



def CI_index( profile , dT = 2, vert_coord = 'PRESSURE' ):
    # Compute the cooling inhibition index given SST cooling of magnitude dT induced by mixing
    profile = make_profile_conscious( profile ).persist()
    # ----------- Find zmix
    zmix = get_zmix( profile , dT = dT, vert_coord = vert_coord );

    # Generate mixed profile
    after = mix_profile( profile, zmix, vert_coord = vert_coord )

    # PE difference
    CI = PE_change( profile, after , vert_coord = vert_coord )
    return CI, zmix
     
