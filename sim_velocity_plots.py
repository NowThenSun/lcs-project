'''
Simulated convection velocity plots
'''
import numpy as np
from plotting_code import velocity_2x2_subplot, velocity_singleplot
from scipy.io import netcdf
import interpolating_functions as interp
import time

hyd = netcdf.netcdf_file('data/hydro10x10Re470.nc', 'r', mmap=False)
print hyd.variables.keys()
U = hyd.variables['U'][:]  #T x Y x X sized array (assume its Y x X order)
V = hyd.variables['V'][:]   # 76*512*512 here
TIME = hyd.variables['TIME'][:]
hyd.close()

# dyn = netcdf.netcdf_file('dynamo10x10Re470.nc', 'r', mmap=False)
# print dyn.variables.keys()
# U_dyn = fh.variables['U'][:]  #T x Y x X sized array (assume its Y x X order)
# V_dyn = fh.variables['V'][:]  # 60*512*512
# TIME_dyn = fh.variables['TIME'][:]
# dyn.close()
nx = 50
ny = 50
nx_dense = 200
ny_dense = 200

t_diff = [0,10,20,30]
t_base = TIME[0]
t = t_diff + t_base
X = np.linspace(0.,10.,512)
Y = np.linspace(0.,10.,512)
X_min,X_max, Y_min, Y_max = (2.,8.,2.,8.)  # Limit initial grid size
aux_grid_spacing = ((X_max-X_min)/nx)*0.08

xx = np.linspace(X_min, X_max, nx)
yy = np.linspace(Y_min, Y_max, ny)
xx_de = np.linspace(X_min,X_max, nx_dense)
yy_de = np.linspace(Y_min,Y_max, ny_dense)
coord_grid = np.array(np.meshgrid(xx,yy,indexing='xy'))
coord_grid_de = np.array(np.meshgrid(xx_de,yy_de,indexing='xy'))

time_single = TIME[10]
vel_single = interp.regular_grid_interpolator_fn(U, V, X, Y, TIME)[0](coord_grid, time_single)
vel_high_res_single =interp.regular_grid_interpolator_fn(U, V, X, Y, TIME)[0](coord_grid_de, time_single)
velocity_singleplot(coords=coord_grid, velocity=vel_single,
    velocity_high_res=vel_high_res_single, xlower=X_min, xupper=X_max, ylower=Y_min, yupper=Y_max, save_name="sim_velocity_plot_TIME[10].pdf", lenunits = False)

# subplot_labels = []
# vel =[]
# vel_high_res = []
#
# for k in xrange(4):
#     # calculate velocities
#     vel.append(interp.regular_grid_interpolator_fn(U, V, X, Y, TIME)[0](coord_grid, t[k]))
#     vel_high_res.append(interp.regular_grid_interpolator_fn(U, V, X, Y, TIME)[0](coord_grid_de,t[k]))
#     subplot_labels.append('t = %.1d' %t[k])

# velocity_2x2_subplot(coords=coord_grid,labels=subplot_labels, velocity=vel,
#     velocity_high_res=vel_high_res, xlower=X_min, xupper=X_max, ylower=Y_min, yupper=Y_max, save_name=False, lenunits = False)
