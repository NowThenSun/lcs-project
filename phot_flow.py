'''
Module to compute and plot FTLE field of photospheric flow data
'''
from __future__ import division
import numpy as np
import general_functions as fn
from scipy.io import netcdf
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
#import matplotlib.animation as anim
import matplotlib as mpl



def regular_grid_interpolator_fn(U, V, X, Y, TIME):
    '''
    Function to input data into regular_grid_interpolator without inputting it directly into the function
    '''
    def regular_grid_interpolator(coords, current_time, bound_interpolation = False):
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

    return regular_grid_interpolator


# Code to access netcdf data file
fh = netcdf.netcdf_file('data/phot_flow_welsch.nc', 'r', mmap=False)

print fh.variables.keys()

X = fh.variables['X'][:]
Y = fh.variables['Y'][:]
U = fh.variables['U'][:]  #T x Y x X sized array (assume its Y x X order)
V = fh.variables['V'][:]
TIME = fh.variables['TIME'][:]

fh.close()

#~~~~~~~~~~~~~~ INITIALISE PARAMETERS ~~~~~~~~~~~~~~~~~~~~~
nx = 400
ny = 400
t_0 = TIME[0]                  # Initial time
aux_grid_spacing = 1
int_time  = 14400 # in seconds (21600s = 6 hrs)
dt_min = np.sign(int_time)*200
dt_max = np.sign(int_time)*2000
adap_error_tol = 20


# Compute nx*ny grid of coordinates
xx = np.linspace(X[0]+1/(2.*nx),X[-1]-1/(2.*nx), nx)
yy = np.linspace(Y[0]+1/(2.*ny),Y[-1]-1/(2.*ny), ny)
coord_grid = np.array(np.meshgrid(xx,yy,indexing='xy'))
# Compute auxiliary grid for differentiation
aux_grid = fn.generate_auxiliary_grid(coord_grid, aux_grid_spacing)

# Perform RKF45 scheme on aux_grid
final_positions = fn.rkf45_loop_fixed_step(
    derivs=regular_grid_interpolator_fn(U, V, X, Y, TIME), aux_grid=aux_grid,
    adaptive_error_tol=adap_error_tol, t_0=t_0,
    int_time=int_time, dt_min=dt_min, dt_max = dt_max)

jac = fn.jacobian_matrix_aux(final_positions,aux_grid_spacing=aux_grid_spacing)
cgst = fn.cauchy_green_tensor(jac)
ev = np.linalg.eigvalsh(cgst)
ev_max = np.amax(ev,-1)
ftle = np.log(ev_max)/(2.*np.abs(int_time))

#
# Plotting code for plot of eigenvalues/FTLE field
fig = plt.figure()
ax = plt.axes()
#mpl.cm.register_cmap(name='mycubehelix',data=mpl._cm.cubehelix(g,s,r,sat)) ##way to add a colourmap directly (from segmentdata) into mpl.cm or something
#cmap = mpl.colors.LinearSegmentedColormap(name='abab', segmentdata =mpl._cm.cubehelix(g,s,r,sat))  ##way to create colormap from colour dictionary
def cubehelix_cmap(g=1.0, s=0.5, r = -1.5, sat = 1.0):
    '''
    ~~~~~~~~~~
    Inputs:
    g : gamma value (can increase intensity of high valued colors for g>1, increase intensity for low values for g<1)
    s : starting color
    r : number of rotations through B -> G -> R
    sat : saturation value (0 for grayscale)
    ~~~~~~~~~~~
    Outputs:
    cubehelix colourmap
    reverse cubehelix colourmap
    '''
    cdict = mpl._cm.cubehelix(g,s,r,sat)
    def rev_fn(f):
        def reverse_f(x):
            return f(1-x)
        return reverse_f
    b_r = rev_fn(cdict['blue'])
    g_r = rev_fn(cdict['green'])
    r_r = rev_fn(cdict['red'])
    cdict_r = {u'blue' : b_r, u'green' : g_r, u'red' : r_r}
    cmap = mpl.colors.LinearSegmentedColormap(name='ch', segmentdata=cdict)
    cmap_r = mpl.colors.LinearSegmentedColormap(name='ch_r', segmentdata=cdict_r)
    return cmap, cmap_r


im = ax.imshow(ftle, interpolation='none', origin='lower', extent=(X[0],X[-1],Y[0],Y[-1]),
    cmap=cubehelix_cmap(g=1.0,s=-1.2,r=-0.85,sat=1.0)[1]) #,aspect='auto' vmin=-0.0001,vmax=0.0001,

ax.text(0.8,1.02,'T = %.1f' %int_time, transform=ax.transAxes)
ax.text(-0.1,1.02,'t_0 = %.1f' %t_0, transform=ax.transAxes)
#ax.text(0.3,1.02,'average dt = %.2e' %np.average(dt), transform=ax.transAxes)
ax.text(0.6,-0.17,'error tol in dt= %r' %adap_error_tol, transform=ax.transAxes)
cbar_ax = fig.add_axes([0.855, 0.15, 0.025, 0.75])
#cbar_ax.set_title('title',fontsize=11,y=1.02,x=1.005)
#ax1.text(0.8,0.9,r'$t$ = %d $\mu$s' %t[T],fontsize=13,transform=ax1.transAxes, color='Azure')
ax.set_xlabel('x')
ax.set_ylabel('y')
cb = fig.colorbar(im, cax=cbar_ax)
plt.show()
