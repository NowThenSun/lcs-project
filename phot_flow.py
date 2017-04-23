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



def plot_phot(int_time, t_0 = TIME[0],g=1,s=0.7,r=1.2,sat=1):
    #~~~~~~~~~~~~~~ INITIALISE PARAMETERS ~~~~~~~~~~~~~~~~~~~~~
    nx = 100
    ny = 100
    # t_0 = TIME[0]                  # Initial time
    # int_time  = 400 # in seconds (21600s = 6 hrs)
    dt_min = np.sign(int_time)*10
    dt_max = np.sign(int_time)*200
    etol = 0.01
    # ~~~~~~~~~~~~~~  SETUP UNITS ~~~~~~~~~~~~~~
    lenunits = 'km' #units of distance for axis labels
    INIT_TIME = 1165933440.0  # Epoch time (rel to Jan 1 1970) in seconds of Dec 12 2006 (initial time of magnetograms)
    # print time.strftime('%d-%b-%Y %H:%M GMT', time.gmtime(INIT_TIME))
    REF_TIME = INIT_TIME - TIME[0]
    # ~~~~~~~~~~~~~~  RESTRICT GRID DOMAIN ~~~~~~~~~~~~~~
    grid_lim_step = 4
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
    #     derivs=interp.regular_grid_interpolator_fn(U, V, X, Y, TIME)[1], aux_grid=aux_grid,
    #     t_0=t_0,
    #     int_time=int_time, dt = 250,return_data=False)
    t_afterloop = time.time()
    print "Time taken to integrate ODE:", t_afterloop - t_beforeloop

    jac = fn.jacobian_matrix_aux(final_positions,aux_grid_spacing=aux_grid_spacing)
    cgst = fn.cauchy_green_tensor(jac)
    ev = np.linalg.eigvalsh(cgst)
    ev_max = np.amax(ev,-1)
    ftle = np.log(ev_max)/(2.*np.abs(int_time))

    time_label = 'T = %+.1f h' %(int_time/3600.)
    # Plotting code for plot of eigenvalues/FTLE field
    #labtop path
    # Google drive path
    plot_name = "C:/Users/Harry/Google Drive/Project IV Lagrangian Coherent Structures/plots/phot/nx%r_t0rel%.1f_T%r_etol%r_dp45.png" %(nx,t_0-TIME[0],int_time,etol)
    plot.FTLE_plot(ftle, X_min, X_max, Y_min, Y_max, int_time=0, t_0=0, adap_error_tol=0,
        colour_range=(-0.0001,0.0001),colour_rescale=0, save_name = plot_name,g=g,s=s,r=r,sat=sat, lenunits=lenunits, label1 = time.strftime('%d-%b-%Y %H:%M UT', time.gmtime((t_0 + REF_TIME))),label2 = time_label)


int_times_array = [-1800]#np.arange(3600,7200,600)
t_0_array = [TIME[-1]]
# t_0_array += TIME[0]
for k in xrange(len(t_0_array)):
    print "Relative starting time for FTLE calculation", t_0_array[k] - TIME[0]
    plot_phot(int_time = int_times_array[0], t_0=t_0_array[k],g=1,s=0.9,r=1.2,sat=1)

# ,g=1,s=0.7,r=1.2,sat=1  #green backwards
# ,g=1,s=-0.9,r=0.9,sat=1 # blue purple forwards
