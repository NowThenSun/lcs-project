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
import matplotlib.pyplot as plt
from plotting_code import cubehelix_cmap

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
    nx = 400
    ny = nx
    # t_0 = TIME[0]                  # Initial time
    # int_time  = 400 # in seconds (21600s = 6 hrs)
    dt_min = np.sign(int_time)*1
    dt_max = np.sign(int_time)*200
    etol = 0.1
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
    plot_name = "C:/Users/Harry/Google Drive/Project IV Lagrangian Coherent Structures/plots/phot/nx%r_t0rel%.1f_T%r_etol%r_dp45.pdf" %(nx,t_0-TIME[0],int_time,etol)
    plot.FTLE_plot(ftle, X_min, X_max, Y_min, Y_max, int_time=0, t_0=0, adap_error_tol=0,
        colour_range=(-0.0001,0.0001),colour_rescale=0, save_name = plot_name,g=g,s=s,r=r,sat=sat, lenunits=lenunits, label1 = time.strftime('%d-%b-%Y %H:%M UT', time.gmtime((t_0 + REF_TIME))),label2 = time_label)

def ftle_phot(nx,ny,int_time,t_0):
    #~~~~~~~~~~~~~~ INITIALISE PARAMETERS ~~~~~~~~~~~~~~~~~~~~~
    dt_min = np.sign(int_time)*1
    dt_max = np.sign(int_time)*200
    etol = 0.1

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

    return ftle, X_min, X_max, Y_min, Y_max, time_label

def subplot1x2_phot(ftle, X_min,X_max,Y_min,Y_max, main_label=0, subplot_labels=0, lenunits=False, colour_range=0,
    g=1,s=0.5,r=0.5,sat=1):
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey='row')#
    #
    cmap = cubehelix_cmap(g=g,s=s,r=r,sat=sat)[1]
    ax1.set(adjustable='box-forced')
    ax2.set(adjustable='box-forced')
    ax1.set_ylim([Y_min,Y_max])
    ax1.set_xlim([X_min,X_max])
    ax2.set_ylim([Y_min,Y_max])
    ax2.set_xlim([X_min,X_max])

    vmin = colour_range[0]
    vmax = colour_range[1]
    im = ax1.imshow(ftle[0], cmap = cmap, origin = 'lower', extent=(X_min,X_max,Y_min,Y_max), vmin=vmin,vmax=vmax)
    ax2.imshow(ftle[1], cmap = cmap, origin = 'lower', extent=(X_min,X_max,Y_min,Y_max), vmin=vmin,vmax=vmax)
    #
    ax1.axes.autoscale(True)
    ax2.axes.autoscale(True)

    if lenunits == False:
        fig.text(0.47, 0.05, 'x (%s)' %lenunits, ha='center', va='center', fontsize=10)  #A different way to add axes labels onto plots
        fig.text(0.05, 0.5, 'y (%s)' %lenunits, ha='center', va='center', rotation='vertical',fontsize=10)
    else:
        fig.text(0.47, 0.23, 'x (%s)' %lenunits, ha='center', va='center', fontsize=10)  #A different way to add axes labels onto plots
        fig.text(0.04, 0.5, 'y (%s)' %lenunits, ha='center', va='center', rotation='vertical',fontsize=10)

    if not main_label == 0:
        fig.text(0.47, 0.77, main_label, ha='center', va='center', fontsize=10)  #A different way to add axes labels onto plots

    if not subplot_labels == 0:
        ax1.set_title(subplot_labels[0],fontsize=10)
        ax2.set_title(subplot_labels[1],fontsize=10)

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.25, 0.03, 0.5])
    cbar_ax.set_title(r'FTLE (s$^{-1}$)',fontsize=10,y=1.02,x=0.6)
    cbar = fig.colorbar(im,cbar_ax,format='%.1e')

    plt.savefig("C:/Users/Harry/Google Drive/Project IV Lagrangian Coherent Structures/plots/phot/phot_fwd_ftle_subplot_comparison.pdf",bbox_inches='tight')
    # plt.savefig('dg_ftle_plot_T_var_A0-1_eps0-2_t00_vmax0-5.pdf', bbox_inches='tight')
    plt.show()

# ~~~~~~~~~~~~~~  SETUP UNITS ~~~~~~~~~~~~~~
lenunits = 'km' #units of distance for axis labels
INIT_TIME = 1165933440.0  # Epoch time (rel to Jan 1 1970) in seconds of Dec 12 2006 (initial time of magnetograms)
# print time.strftime('%d-%b-%Y %H:%M GMT', time.gmtime(INIT_TIME))
REF_TIME = INIT_TIME - TIME[0]

t0 = TIME[0]+6*3600
int_time=[4*3600,6*3600]
nx_res = 400
ftle_1 = ftle_phot(nx=nx_res,ny=nx_res,int_time=int_time[0], t_0 = t0)
ftle_2 = ftle_phot(nx=nx_res,ny=nx_res,int_time=int_time[1], t_0 = t0)
# print a[1:-1]
subplot1x2_phot((ftle_1[0],ftle_2[0]),*ftle_1[1:-1], main_label=time.strftime('%d-%b-%Y %H:%M UT', time.gmtime(t0+REF_TIME)),
    subplot_labels=(ftle_1[-1],ftle_2[-1]),lenunits=lenunits, colour_range=(-0.0001,0.0001),g=1,s=-0.9,r=0.9,sat=1)
#
# int_times_hrs = [-12]
# int_times_array = np.array(int_times_hrs)*3600
# t_0_array = [TIME[0]+12*3600]
# # t_0_array += TIME[0]
# for k in xrange(len(t_0_array)):
#     print "Relative starting time for FTLE calculation", t_0_array[k] - TIME[0]
#     for j in xrange(len(int_times_array)):
#         plot_phot(int_time = int_times_array[j], t_0=t_0_array[k],g=1,s=0.7,r=1.2,sat=1)

# ,g=1,s=0.7,r=1.2,sat=1  #green backwards
# ,g=1,s=-0.9,r=0.9,sat=1 # blue purple forwards
