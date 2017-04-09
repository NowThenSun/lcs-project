'''
Module to compute and plot FTLE field of photospheric flow data
'''
from __future__ import division
import numpy as np
import general_functions as fn
from scipy.io import netcdf
from scipy.interpolate import RegularGridInterpolator
import plotting_code as plot
import ODE_solvers as ODE
import time
import interpolating_functions as interp


# Code to access netcdf data file
fh = netcdf.netcdf_file('data/phot_flow_welsch.nc', 'r', mmap=False)

print fh.variables.keys()

X = fh.variables['X'][:]
Y = fh.variables['Y'][:]
U = fh.variables['U'][:]  #T x Y x X sized array (assume its Y x X order)
V = fh.variables['V'][:]
TIME = fh.variables['TIME'][:]

fh.close()
def plot_phot(int_time, t_0 = TIME[0]):
    #~~~~~~~~~~~~~~ INITIALISE PARAMETERS ~~~~~~~~~~~~~~~~~~~~~
    nx = 150
    ny = 150
    # t_0 = TIME[0]                  # Initial time
    # int_time  = 400 # in seconds (21600s = 6 hrs)
    dt_min = np.sign(int_time)*10
    dt_max = np.sign(int_time)*250
    etol = 0.0001

    grid_lim_step = 2
    X_min,X_max, Y_min, Y_max = (X[grid_lim_step],X[-grid_lim_step-1],Y[grid_lim_step],Y[-grid_lim_step-1])  # Limit initial grid size
    aux_grid_spacing = ((X_max-X_min)/nx)*0.08
    # Compute nx*ny grid of coordinates
    xx = np.linspace(X_min,X_max, nx)
    yy = np.linspace(Y_min,Y_max, ny)
    coord_grid = np.array(np.meshgrid(xx,yy,indexing='xy'))
    # Compute auxiliary grid for differentiation
    aux_grid = fn.generate_auxiliary_grid(coord_grid, aux_grid_spacing)
    # Perform RKF45 scheme on aux_grid
    t_beforeloop = time.time()
    final_positions = ODE.dp45_loop(
        derivs=interp.regular_grid_interpolator_fn(U, V, X, Y, TIME)[1], aux_grid=aux_grid,
        t_0=t_0,
        int_time=int_time, dt_min=dt_min, dt_max = dt_max,maxiters = 1000, atol=etol, rtol =etol)#00001)
    # final_positions = fn.rk4_loop(
    #     derivs=regular_grid_interpolator_fn(U, V, X, Y, TIME)[1], aux_grid=aux_grid,
    #     t_0=t_0,
    #     int_time=int_time, dt = 250,return_data=False)
    t_afterloop = time.time()
    print "Time taken to integrate ODE:", t_afterloop - t_beforeloop

    jac = fn.jacobian_matrix_aux(final_positions,aux_grid_spacing=aux_grid_spacing)
    cgst = fn.cauchy_green_tensor(jac)
    ev = np.linalg.eigvalsh(cgst)
    ev_max = np.amax(ev,-1)
    ftle = np.log(ev_max)/(2.*np.abs(int_time))

    #
    # Plotting code for plot of eigenvalues/FTLE field
    #labtop path
    # Google drive path
    plot_name = "C:/Users/Harry/Google Drive/Project IV Lagrangian Coherent Structures/plots/nx%r_T%r_etol%r_t0rel%r.png" %(nx,int_time,etol,t_0-TIME[0])
    plot.FTLE_plot(ftle, X_min, X_max, Y_min, Y_max, int_time, t_0-TIME[0], adap_error_tol=0, colour_range=(-0.0001,0.0001), save_name = plot_name,g=1,s=0.8,r=1.2,sat=1)


int_times_array = [14400]#np.arange(3600,7200,600)
t_0_array = [0,5000,10000,15000]
t_0_array += TIME[0]
for k in xrange(len(t_0_array)):
    print "Starting FTLE calculation for time", t_0_array[k]
    plot_phot(int_time = int_times_array[0], t_0=t_0_array[k])
