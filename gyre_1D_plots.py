'''
Module that creates 1D plots over ridges in the double gyre field to compare errors in peaks at different resolutions, methods and different time steps
compare to RKF45 that alternates all of the timesteps
Can also compare dgyre analytic to different interpolation functions
'''
from __future__ import division
import numpy as np
import general_functions as fn
import double_gyre_functions as dg
import matplotlib.pyplot as plt



#~~~~~ 1D PLOTTING PARAMETERS ~~~~~
y = 0.3                                  # fixed y-value the FTLE is evaluated at
x0 = 1.0                                # first x value
x1 = 1.2                                # final x value
res = np.array([50,300,1000]) # Resolutions tested at (number of grid spacing between x0 and x1)               0.2/100 ~ 1000x500 res
actual_res = res*2./(x1-x0)     # True resolution on a full [0,2]x[0,1] double gyre domain
# double gyre parameters - Note: dg.gyre_global_params's global parameters have dg. in front of them so primarily exist within the dg module (although can still be called upon here with e.g. dg.amplitude_g)
dg_params = dg.gyre_global_params(amplitude=0.1, epsilon=0.1, omega=2*np.pi/10.)
# other parameters
aux_grid_spacing = 1*10**-5
adaptive_error_tol = 1*10**-2
t_0 = 2.
int_time = 10.
dt_min = np.sign(int_time)*1*10**-3
dt_max = 1.
dt_fixed = 0.5      # timestep used for RK4 integration

#list to put in final data of calculated FTLE fields
FTLE_list = []

for j in xrange(len(res)):
    print "RESOLUTION =", actual_res[j]
    # Generate initial grid and aux. grid
    initial_grid = np.array(np.meshgrid(np.linspace(x0,x1,res[j]),y,indexing='xy'))  # 2*1*res[j] grid
    a_grid = fn.generate_auxiliary_grid(initial_grid, aux_grid_spacing=aux_grid_spacing)
    print np.shape(a_grid)
    # Use rkf45 integration scheme
    final_pos = fn.rkf45_loop(derivs=dg.analytic_velocity, aux_grid=a_grid, adaptive_error_tol=adaptive_error_tol, t_0=t_0, int_time=int_time, dt_min=dt_min, dt_max = dt_max)
    jac = fn.jacobian_matrix_aux(final_pos,aux_grid_spacing=aux_grid_spacing)
    cgst = fn.cauchy_green_tensor(jac)
    ev = np.linalg.eigvalsh(cgst)
    ev_max = np.amax(ev,-1)
    ftle = np.log(ev_max)/(2.*np.abs(int_time))
    # Append FTLE data to FTLE_list
    FTLE_list.append(ftle)

    # Now use rk4 scheme for
    final_pos = fn.rk4_loop(derivs=dg.analytic_velocity, aux_grid=a_grid, t_0=t_0, int_time=int_time, dt=dt_fixed)
    jac = fn.jacobian_matrix_aux(final_pos,aux_grid_spacing=aux_grid_spacing)
    cgst = fn.cauchy_green_tensor(jac)
    ev = np.linalg.eigvalsh(cgst)
    ev_max = np.amax(ev,-1)
    ftle = np.log(ev_max)/(2.*np.abs(int_time))
    # Append FTLE data to FTLE_list
    FTLE_list.append(ftle)

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
plt.show()
