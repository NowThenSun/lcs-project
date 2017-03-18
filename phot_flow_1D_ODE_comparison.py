'''
Module for comparing the 1D plots of different ODE methods for photospheric flow data
'''
from __future__ import division
import numpy as np
import general_functions as fn
import matplotlib.pyplot as plt
import ODE_solvers as ODE
import interpolating_functions as interp
from scipy.io import netcdf

# Code to access netcdf data file
fh = netcdf.netcdf_file('data/phot_flow_welsch.nc', 'r', mmap=False)

print fh.variables.keys()

X = fh.variables['X'][:]
Y = fh.variables['Y'][:]
U = fh.variables['U'][:]  #T x Y x X sized array (assume its Y x X order)
V = fh.variables['V'][:]
TIME = fh.variables['TIME'][:]

fh.close()

#~~~~~ 1D PLOTTING PARAMETERS ~~~~~
y = 8000                                # fixed y-value the FTLE is evaluated at
x0 = 4000                                  # initial x value
x1 = 6000                                # final x value
res = np.array([50,300,500]) # Resolutions tested at (number of grid points between x0 and x1)               0.2/100 ~ 1000x500 res
grid_spacing = (x1-x0)/res
# double gyre parameters - Note: dg.gyre_global_params's global parameters have dg. in front of them so primarily exist within the dg module (although can still be called upon here with e.g. dg.amplitude_g)
dg_params = dg.gyre_global_params(amplitude=0.1, epsilon=0.1, omega=2*np.pi/10.)
# other parameters
aux_grid_spacing = grid_spacing*0.08
rkf45_error_tol = 1*10**-3
dp45_error_tol = rkf45_error_tol
t_0 = 2.
int_time = 10.
dt_min = np.sign(int_time)*10
dt_max = np.sign(int_time)*250
dt_fixed = 250      # timestep used for RK4 integration


# Bools of which methods are compared
RK4 = True
RKF45 = True
DP45 = True
no_methods = [RK4,RKF45,DP45].count(True)
#list to put in final data of calculated FTLE fields

FTLE_list_RK4 = []
FTLE_list_RKF45 = []
FTLE_list_DP45 = []


for j in xrange(len(res)):
    print "Grid spacing =", grid_spacing[j]
    # Generate initial grid and aux. grid
    initial_grid = np.array(np.meshgrid(np.linspace(x0,x1,res[j]),y,indexing='xy'))  # 2*1*res[j] grid
    a_grid = fn.generate_auxiliary_grid(initial_grid, aux_grid_spacing=aux_grid_spacing)
    print np.shape(a_grid)

    if RKF45 == True:
        # Use rkf45 integration scheme
        final_pos = fn.rkf45_loop(derivs=dg.analytic_velocity, aux_grid=a_grid, adaptive_error_tol=adaptive_error_tol, t_0=t_0, int_time=int_time, dt_min=dt_min, dt_max = dt_max)
        jac = fn.jacobian_matrix_aux(final_pos,aux_grid_spacing=aux_grid_spacing)
        cgst = fn.cauchy_green_tensor(jac)
        ev = np.linalg.eigvalsh(cgst)
        ev_max = np.amax(ev,-1)
        ftle = np.log(ev_max)/(2.*np.abs(int_time))
        # Append FTLE data to FTLE_list
        FTLE_list_RKF45.append(ftle)

    if RK4 == True:
        # Now use rk4 scheme for
        final_pos = fn.rk4_loop(derivs=interp.regular_grid_interpolator_fn(U, V, X, Y, TIME)[1], aux_grid=a_grid, t_0=t_0, int_time=int_time, dt=dt_fixed)
        jac = fn.jacobian_matrix_aux(final_pos,aux_grid_spacing=aux_grid_spacing)
        cgst = fn.cauchy_green_tensor(jac)
        ev = np.linalg.eigvalsh(cgst)
        ev_max = np.amax(ev,-1)
        ftle = np.log(ev_max)/(2.*np.abs(int_time))
        # Append FTLE data to FTLE_list
        FTLE_list_RK4.append(ftle)

    if DP45 == True:
        #DP45 scheme
        final_pos = fn.dp45_loop(
            derivs=interp.regular_grid_interpolator_fn(U, V, X, Y, TIME)[1], aux_grid=aux_grid_spacing,
            t_0=t_0, int_time=int_time,
            dt_min=dt_min, dt_max = dt_max,maxiters = 1000, atol=0.0001, rtol =0.0001
            )
        cgst = fn.cauchy_green_tensor(jac)
        ev = np.linalg.eigvalsh(cgst)
        ev_max = np.amax(ev,-1)
        ftle = np.log(ev_max)/(2.*np.abs(int_time))
        FTLE_list_DP45.append(ftle)

#Plotting code
print len(FTLE_list)

fig = plt.figure()
ax = plt.axes()
ax.set_xlabel('x')
ax.set_ylabel('FTLE')
legends = []
for j in xrange(len(res)):
    plt.plot(np.linspace(x0,x1,res[j]),FTLE_list[2*j][0], '--')
    legends.append('RKF45 %ix%i' %(actual_res[j],actual_res[j]/2))
    plt.plot(np.linspace(x0,x1,res[j]),FTLE_list[2*j+1][0], '-')
    legends.append('RK4 %ix%i dt = %.2f' %(actual_res[j],actual_res[j]/2, dt_fixed))
plt.legend(legends)
#plt.savefig('testfig23.pdf',transparent=True)
plt.show()
