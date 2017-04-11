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
y = 5000                                # fixed y-value the FTLE is evaluated at
x0 = X[0]                                  # initial x value
x1 = 5500                                # final x value
res = np.array([200,500]) # Resolutions tested at (number of grid points between x0 and x1)               0.2/100 ~ 1000x500 res
grid_spacing = (x1-x0)/res
#True resolution calculation
# Main parameter (need to make sure this matches with grid_lim_step used in data calculation)
grid_lim_step = 2
X_min,X_max, Y_min, Y_max = (X[grid_lim_step],X[-grid_lim_step-1],Y[grid_lim_step],Y[-grid_lim_step-1])  # Limit initial grid size
actual_res = (X_max-X_min)*res/(x1-x0)

# other parameters
aux_grid_spacing = grid_spacing*0.08
rkf45_error_tol = 1*10**2
dp45_error_tol = rkf45_error_tol
t_0 = TIME[-1]
int_time = -14400
dt_min = np.sign(int_time)*10
dt_max = np.sign(int_time)*300
dt_fixed = np.sign(int_time)*300      # timestep used for RK4 integration


# Bools of which methods are compared
Truez=False
RK4 = True
RKF45 = True
DP45 = True
no_methods = [RK4,RKF45,DP45].count(True)
#list to put in final data of calculated FTLE fields

FTLE_list_RK4 = []
FTLE_list_RKF45 = []
FTLE_list_DP45 = []
fig = plt.figure()
ax = plt.axes()
ax.set_xlabel('x')
ax.set_ylabel('FTLE')
legends = []
line_alpha = 0.66
dp_line_alpha = 1.

for j in xrange(len(res)):
    print "Grid spacing =", grid_spacing[j]
    # Generate initial grid and aux. grid
    initial_grid = np.array(np.meshgrid(np.linspace(x0,x1,res[j]),y,indexing='xy'))  # 2*1*res[j] grid
    print np.shape(initial_grid)
    a_grid = fn.generate_auxiliary_grid(initial_grid, aux_grid_spacing=aux_grid_spacing[j])
    print np.shape(a_grid)

    if RK4 == True:
        # Now use rk4 scheme for
        final_pos = fn.rk4_loop(derivs=interp.regular_grid_interpolator_fn(U, V, X, Y, TIME)[1], aux_grid=a_grid, t_0=t_0, int_time=int_time, dt=dt_fixed)

        jac = fn.jacobian_matrix_aux(final_pos,aux_grid_spacing=aux_grid_spacing[j])
        cgst = fn.cauchy_green_tensor(jac)
        ev = np.linalg.eigvalsh(cgst)
        ev_max = np.amax(ev,-1)
        ftle = np.log(ev_max)/(2.*np.abs(int_time))
        # Append FTLE data to FTLE_list
        FTLE_list_RK4.append(ftle)

        plt.plot(np.linspace(x0,x1,res[j]),FTLE_list_RK4[j][0],'-', alpha=line_alpha)
        legends.append('RK4 nx%i dt = %.2f' %(actual_res[j], dt_fixed))

    if RKF45 == True:
        # Use rkf45 integration scheme
        final_positions = ODE.rkf45_loop(
        derivs=interp.regular_grid_interpolator_fn(U, V, X, Y, TIME)[1], aux_grid=a_grid,
        t_0=t_0,
        int_time=int_time, dt_min=dt_min, dt_max = dt_max,maxiters = 1000, atol=rkf45_error_tol, rtol =rkf45_error_tol)

        jac = fn.jacobian_matrix_aux(final_positions,aux_grid_spacing=aux_grid_spacing[j])
        cgst = fn.cauchy_green_tensor(jac)
        ev = np.linalg.eigvalsh(cgst)
        ev_max = np.amax(ev,-1)
        ftle = np.log(ev_max)/(2.*np.abs(int_time))
        # Append FTLE data to FTLE_list
        FTLE_list_RKF45.append(ftle)
        plt.plot(np.linspace(x0,x1,res[j]),FTLE_list_RKF45[j][0],'--', alpha=line_alpha)
        legends.append('RKF45 nx%i etol = %r' %(actual_res[j], rkf45_error_tol))

    if DP45 == True:
        #DP45 scheme
        final_positions = ODE.dp45_loop(
            derivs=interp.regular_grid_interpolator_fn(U, V, X, Y, TIME)[1], aux_grid=a_grid,
            t_0=t_0,
            int_time=int_time, dt_min=dt_min, dt_max = dt_max,maxiters = 1000, atol=dp45_error_tol, rtol =dp45_error_tol)

        jac = fn.jacobian_matrix_aux(final_positions,aux_grid_spacing=aux_grid_spacing[j])
        cgst = fn.cauchy_green_tensor(jac)
        ev = np.linalg.eigvalsh(cgst)
        ev_max = np.amax(ev,-1)
        ftle = np.log(ev_max)/(2.*np.abs(int_time))
        FTLE_list_DP45.append(ftle)
        plt.plot(np.linspace(x0,x1,res[j]),FTLE_list_DP45[j][0],':', alpha=dp_line_alpha, marker='x')
        legends.append('DP45 nx%i etol = %r' %(actual_res[j], dp45_error_tol))



#
plt.legend(legends)
ax.text(0.8,1.02,'y = %r' %y, transform=ax.transAxes)
ax.text(0.2,1.02,'t_0 (relative) = %r' %(t_0-TIME[0]), transform=ax.transAxes)

#plt.savefig('testfig23.pdf',transparent=True)
plt.show()
