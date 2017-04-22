'''
Module containing different interpolation functions in 3D to be easily imported into other modules
'''
from __future__ import division
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import scipy.interpolate as si

def cubic_3D_splines_fn(U, V, X, Y, TIME):
    '''

    '''
    def interpolator(y, derivs):
        '''

        '''


def regular_grid_interpolator_fn(U, V, X, Y, TIME):
    '''
    Function to input data into regular_grid_interpolator without inputting it directly into the function
    '''
    def regular_grid_interpolator_scalar_t(coords, current_time, bound_interpolation = False):
        '''
        Function for interpolating velocity
        ~~~~~~~~~
        Inputs:
        coords = coordinates to find the interpolated velocities at 2*ny*nx*4
        current_time = time to find the interpolated velocity at (scalar)
        U, V = original velocity data that is used for the interpolation (len(TIME)*len(Y)*len(X) shaped arrays)
        X,Y,TIME = original 1D arrays of coordinates that velocity data is on
        bound_interpolation = if set to True the interpolator will raise error when used outside of data range
            if set to False the interpolator will extrapolate outside of the given data
        ~~~~~~~~~~
        Outputs:
        Interpolated velocity at coordinates and current_time inputted
        '''
        #t_repeat = np.repeat(current_time, 4, axis=-1) #grid repeated to 1*ny*nx*4 shape for broadcasting
        t_repeat = np.zeros((1,np.shape(coords)[1],np.shape(coords)[2]))+current_time
        #print np.shape(t_repeat)
        mesh_coords = np.concatenate((t_repeat,coords)) # 3*ny*nx shape now (T,X,Y) indexing  (but X = Y here - so doesn't matter)
        mesh_coords1 = mesh_coords[::-1] # Flips to (Y,X,T) indexing
        mesh_coords2 = np.roll(mesh_coords1, 1, axis=0) #Rolls to (T,Y,X) indexing
        #Need to meshgrid coordinates into a N*3 form to input into regulargridinterpolator
        mc3 = np.rollaxis(mesh_coords2,0,3)  #ny*nx*3
        # print np.shape(coords)
        # print np.shape(mc3)
        shape = np.int(np.prod(np.shape(coords))/2)
        # print shape
        mc4 = np.reshape(mc3, (shape,3)) #(ny*nx)*3
        # New attempt at mesh coords


        #print np.shape(time)
        Uint = RegularGridInterpolator((TIME,Y,X),U,bounds_error=bound_interpolation,fill_value=None)
        Vint = RegularGridInterpolator((TIME,Y,X),V,bounds_error=bound_interpolation,fill_value=None)  #fill_value = None extrapolates
        U1 = Uint(mc4).reshape(np.shape(coords)[1],np.shape(coords)[2])  #insert N*3 return N shaped arrays
        V1 = Vint(mc4).reshape(np.shape(coords)[1],np.shape(coords)[2])
        return U1, V1

    def regular_grid_interpolator_array_t(coords, current_time, bound_interpolation = False):
        '''
        Function for interpolating velocity
        ~~~~~~~~~
        Inputs:
        coords = coordinates to find the interpolated velocities at 2*ny*nx*4
        current_time = time to find the interpolated velocity at (array)
        U, V = original velocity data that is used for the interpolation (len(TIME)*len(Y)*len(X) shaped arrays)
        X,Y,TIME = original 1D arrays of coordinates that velocity data is on
        bound_interpolation = if set to True the interpolator will raise error when used outside of data range
            if set to False the interpolator will extrapolate outside of the given data
        ~~~~~~~~~~
        Outputs:
        Interpolated velocity at coordinates and current_time inputted
        '''
        #print np.shape(coords)
        #t_repeat = np.repeat(current_time, 4, axis=-1) #grid repeated to 1*ny*nx*4 shape for broadcasting
        t_repeat = np.zeros((1,np.shape(coords)[1],np.shape(coords)[2],4))+current_time
        #print t_repeat
        #print np.shape(current_time), np.shape(t_repeat)
        mesh_coords = np.concatenate((t_repeat,coords)) # 3*ny*nx*n4 shape now (T,X,Y) indexing  (but X = Y here - so doesn't matter)
        #join = np.concatenate((current_time[np.newaxis],coords))
        mesh_coords1 = mesh_coords[::-1] # Flips to (Y,X,T) indexing
        mesh_coords2 = np.roll(mesh_coords1, 1, axis=0) #Rolls to (T,Y,X) indexing (moves the final axes to position before 1)
        #Need to meshgrid coordinates into a N*3 form to input into regulargridinterpolator
        mc3 = np.swapaxes(mesh_coords2 , 0, -1)  #4*ny*nx*3
        shape = np.prod(np.shape(coords))/2
        mc4 = np.reshape(mc3, (shape,3)) #(4*ny*nx)*3
        # New attempt at mesh coords
        #print np.shape(coords), np.shape(current_time)

        #print np.shape(time)
        Uint = RegularGridInterpolator((TIME,Y,X),U,bounds_error=bound_interpolation,fill_value=None)
        Vint = RegularGridInterpolator((TIME,Y,X),V,bounds_error=bound_interpolation,fill_value=None)  #fill_value = None extrapolates
        U1 = Uint(mc4).reshape(4,np.shape(coords)[1],np.shape(coords)[2])  #insert N*3 return N shaped arrays
        V1 = Vint(mc4).reshape(4,np.shape(coords)[1],np.shape(coords)[2])

        return np.rollaxis(U1,0,3),np.rollaxis(V1,0,3)



    return regular_grid_interpolator_scalar_t, regular_grid_interpolator_array_t
