'''
Module for calculating the FTLE field of the simulated hydrodynamc and dynamo flow data
'''
from __future__ import division
import numpy as np
import general_functions as fn
import ODE_solvers as ODE

from scipy.io import netcdf
import plotting_code as plot
import time
import interpolating_functions as interp

hyd = netcdf.netcdf_file('data/hydro10x10Re470.nc', 'r', mmap=False)
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
nx = 300
ny = 300
t_0 = TIME_hyd[10]                  # Initial time
int_time  = 5#hydro data goes ~211-264
dt_min = np.sign(int_time)*0.005
dt_max = np.sign(int_time)*0.01
adaptive_error_tol = 10**-3

# Compute nx*ny grid of coordinates
X = np.linspace(0.,10.,512)
Y = np.linspace(0.,10.,512)
X_min,X_max, Y_min, Y_max = (2.,8.,2.,8.)  # Limit initial grid size
aux_grid_spacing = ((X_max-X_min)/nx)*0.08


xx = np.linspace(X_min, X_max, nx)
yy = np.linspace(Y_min, Y_max, ny)

coord_grid = np.array(np.meshgrid(xx,yy,indexing='xy'))
# Compute auxiliary grid for differentiation
aux_grid = fn.generate_auxiliary_grid(coord_grid, aux_grid_spacing)

t_beforeloop = time.time()
# Perform RKF45 scheme on aux_grid
# final_positions = fn.rkf45_loop_fixed_step(
#     derivs=regular_grid_interpolator_fn(U_hyd, V_hyd, X, Y, TIME_hyd)[1], aux_grid=aux_grid,
#     adaptive_error_tol=adap_error_tol, t_0=t_0,
#     int_time=int_time, dt_min=dt_min, dt_max = dt_max)
# final_positions = fn.rk4_loop(
#     derivs=regular_grid_interpolator_fn(U_hyd, V_hyd, X, Y, TIME_hyd)[1], aux_grid=aux_grid,
#     t_0=t_0,
#     int_time=int_time, dt = 0.2,return_data=False)
# NEW RKF45 SCHEME
final_positions = ODE.dp45_loop(derivs=interp.regular_grid_interpolator_fn(U_hyd, V_hyd, X, Y, TIME_hyd)[1], aux_grid=aux_grid,
    t_0=t_0, int_time=int_time, dt_min=dt_min,dt_max= dt_max,
    maxiters = 1000, atol = adaptive_error_tol, rtol = adaptive_error_tol)

t_afterloop = time.time()
print "Time taken to integrate ODE:", t_afterloop - t_beforeloop

#print final_positions
jac = fn.jacobian_matrix_aux(final_positions,aux_grid_spacing=aux_grid_spacing)
cgst = fn.cauchy_green_tensor(jac)
ev = np.linalg.eigvalsh(cgst)
ev_max = np.amax(ev,-1)
#print ev_max, np.average(ev_max)
ftle = np.log(ev_max)/(2.*np.abs(int_time))


#
# Plotting code for plot of eigenvalues/FTLE field
plot.FTLE_plot(ftle, X_min, X_max, Y_min, Y_max, int_time, t_0, adap_error_tol=adaptive_error_tol)#, colour_range=(0,0.1))
