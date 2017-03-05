'''
Module for calculating the FTLE field of the simulated hydrodynamc and dynamo flow data
'''
from __future__ import division
import numpy as np
import general_functions as fn
from scipy.io import netcdf
from scipy.interpolate import RegularGridInterpolator
import plotting_code as plot

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
        t_repeat = np.zeros((1,np.shape(coords)[1],np.shape(coords)[2],4))+current_time
        #print np.shape(t_repeat)
        mesh_coords = np.concatenate((t_repeat,coords)) # 3*ny*nx*n4 shape now (T,X,Y) indexing  (but X = Y here - so doesn't matter)
        mesh_coords1 = mesh_coords[::-1] # Flips to (Y,X,T) indexing
        mesh_coords2 = np.roll(mesh_coords1, 1, axis=0) #Rolls to (T,Y,X) indexing
        #Need to meshgrid coordinates into a N*3 form to input into regulargridinterpolator
        mc3 = np.swapaxes(mesh_coords2 , 0, -1)  #4*ny*nx*3
        shape = np.prod(np.shape(coords))/2
        mc4 = np.reshape(mc3, (shape,3)) #(4*ny*nx)*3
        # New attempt at mesh coords


        #print np.shape(time)
        Uint = RegularGridInterpolator((TIME,Y,X),U,bounds_error=bound_interpolation,fill_value=None)
        Vint = RegularGridInterpolator((TIME,Y,X),V,bounds_error=bound_interpolation,fill_value=None)  #fill_value = None extrapolates
        U1 = Uint(mc4).reshape(4,np.shape(coords)[1],np.shape(coords)[2])  #insert N*3 return N shaped arrays
        V1 = Vint(mc4).reshape(4,np.shape(coords)[1],np.shape(coords)[2])
        #print U1
        return np.rollaxis(U1,0,3),np.rollaxis(V1,0,3)

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



hyd = netcdf.netcdf_file('data/dynamo10x10Re470.nc', 'r', mmap=False)
print hyd.variables.keys()
U_hyd = hyd.variables['U'][:]  #T x Y x X sized array (assume its Y x X order)
V_hyd = hyd.variables['V'][:]   # 76*512*512 here
TIME_hyd = hyd.variables['TIME'][:]
hyd.close()

# dyn = netcdf.netcdf_file('dynamo10x10Re470.nc', 'r', mmap=False)
# print dyn.variables.keys()
# U_dyn = fh.variables['U'][:]  #T x Y x X sized array (assume its Y x X order)
# V_dyn = fh.variables['V'][:]  # 60*512*512
# TIME_dyn = fh.variables['TIME'][:]
# dyn.close()
# print U_hyd
# print TIME_hyd[-1], TIME_hyd[0]
#~~~~~~~~~~~~~~ INITIALISE PARAMETERS ~~~~~~~~~~~~~~~~~~~~~
nx = 100
ny = 100
t_0 = TIME_hyd[0]                  # Initial time
aux_grid_spacing = 10**-3
int_time  = 40 #hydro data goes ~211-264
dt_min = np.sign(int_time)*0.01
dt_max = np.sign(int_time)*5
adap_error_tol = 4

# Compute nx*ny grid of coordinates
X = np.linspace(0.,10.,512)
Y = np.linspace(0.,10.,512)
xx = np.linspace(X[0]+1/(2.*nx),X[-1]-1/(2.*nx), nx)
yy = np.linspace(Y[0]+1/(2.*ny),Y[-1]-1/(2.*ny), ny)
coord_grid = np.array(np.meshgrid(xx,yy,indexing='xy'))
# Compute auxiliary grid for differentiation
aux_grid = fn.generate_auxiliary_grid(coord_grid, aux_grid_spacing)

# Perform RKF45 scheme on aux_grid
# final_positions = fn.rkf45_loop_fixed_step(
#     derivs=regular_grid_interpolator_fn(U_hyd, V_hyd, X, Y, TIME_hyd)[1], aux_grid=aux_grid,
#     adaptive_error_tol=adap_error_tol, t_0=t_0,
#     int_time=int_time, dt_min=dt_min, dt_max = dt_max)
final_positions = fn.rk4_loop(
    derivs=regular_grid_interpolator_fn(U_hyd, V_hyd, X, Y, TIME_hyd)[1], aux_grid=aux_grid,
    t_0=t_0,
    int_time=int_time, dt = 1)
jac = fn.jacobian_matrix_aux(final_positions,aux_grid_spacing=aux_grid_spacing)
cgst = fn.cauchy_green_tensor(jac)
ev = np.linalg.eigvalsh(cgst)
ev_max = np.amax(ev,-1)
ftle = np.log(ev_max)/(2.*np.abs(int_time))

#
# Plotting code for plot of eigenvalues/FTLE field
plot.FTLE_plot(ftle, X[0], X[-1], Y[0], Y[-1], int_time, t_0, adap_error_tol)#, colour_range=(0,0.8))