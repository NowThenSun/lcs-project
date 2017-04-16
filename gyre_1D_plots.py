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
import ODE_solvers as ODE
import matplotlib


#~~~~~ 1D PLOTTING PARAMETERS ~~~~~
y = 0.3                                  # fixed y-value the FTLE is evaluated at
x0 = 1.0                                # first x value
x1 = 1.2                                # final x value
res = np.array([500,5000]) # Resolutions tested at (number of grid spacing between x0 and x1)               0.2/100 ~ 1000x500 res
actual_res = res*2./(x1-x0)     # True resolution on a full [0,2]x[0,1] double gyre domain
grid_spacing = (x1-x0)/res

# double gyre parameters - Note: dg.gyre_global_params's global parameters have dg. in front of them so primarily exist within the dg module (although can still be called upon here with e.g. dg.amplitude_g)
dg_params = dg.gyre_global_params(amplitude=0.1, epsilon=0.1, omega=2*np.pi/10.)
# other parameters
aux_grid_spacing = 0.08*grid_spacing
rkf45_error_tol = 1*10**-4
dp45_error_tol = rkf45_error_tol
t_0 = 2.
int_time = 10.
dt_min = np.sign(int_time)*1*10**-3
dt_max = 1.
# dt_fixed = 0.5      # timestep used for RK4 integration
dt_fixed_array = [0.5]        # array of timesteps used for RK4 integration

# Bools of which methods are compared
Truez=False
RK4 = True
RKF45 = True
DP45 = True
no_methods = [RK4,RKF45,DP45].count(True)
#list to put in final data of calculated FTLE fields

fig = plt.figure()
ax = plt.axes()

ax.set_xlabel('x')
ax.set_ylabel('FTLE')
legends = []
line_alpha = 0.66
dp_line_alpha = 0.7

for j in xrange(len(res)):
    print "Grid spacing =", grid_spacing[j]
    # Generate initial grid and aux. grid
    initial_grid = np.array(np.meshgrid(np.linspace(x0,x1,res[j]),y,indexing='xy'))  # 2*1*res[j] grid
    print np.shape(initial_grid)
    a_grid = fn.generate_auxiliary_grid(initial_grid, aux_grid_spacing=aux_grid_spacing[j])
    print np.shape(a_grid)

    if RK4 == True:
        for k in xrange(len(dt_fixed_array)):
            # Now use rk4 scheme for
            final_pos = fn.rk4_loop(derivs=dg.analytic_velocity, aux_grid=a_grid, t_0=t_0, int_time=int_time, dt=dt_fixed_array[k])

            jac = fn.jacobian_matrix_aux(final_pos,aux_grid_spacing=aux_grid_spacing[j])
            cgst = fn.cauchy_green_tensor(jac)
            ev = np.linalg.eigvalsh(cgst)
            ev_max = np.amax(ev,-1)
            ftle = np.log(ev_max)/(2.*np.abs(int_time))
            # Append FTLE data to FTLE_list
            # FTLE_list_RK4.append(ftle)

            plt.plot(np.linspace(x0,x1,res[j]),ftle[0],'-', alpha=line_alpha)
            legends.append(r'RK4: $\mathregular{h_\mathrm{grid}}$=%.1e  , $\mathregular{\Delta t}$ = %.1f' %(grid_spacing[j], dt_fixed_array[k]))

    if RKF45 == True:
        # Use rkf45 integration scheme
        final_positions = ODE.rkf45_loop(
        derivs=dg.analytic_velocity, aux_grid=a_grid,
        t_0=t_0,
        int_time=int_time, dt_min=dt_min, dt_max = dt_max,maxiters = 1000, atol=rkf45_error_tol, rtol =rkf45_error_tol)

        jac = fn.jacobian_matrix_aux(final_positions,aux_grid_spacing=aux_grid_spacing[j])
        cgst = fn.cauchy_green_tensor(jac)
        ev = np.linalg.eigvalsh(cgst)
        ev_max = np.amax(ev,-1)
        ftle = np.log(ev_max)/(2.*np.abs(int_time))
        # Append FTLE data to FTLE_list
        # FTLE_list_RKF45.append(ftle)
        plt.plot(np.linspace(x0,x1,res[j]),ftle[0],'--', alpha=line_alpha)
        legends.append(r'RKF45: $\mathregular{h_{grid}}$=%.1e  , $\mathregular{\tau}$ =%.0e' %(grid_spacing[j], rkf45_error_tol))

    if DP45 == True:
        #DP45 scheme
        final_positions = ODE.dp45_loop(
            derivs=dg.analytic_velocity, aux_grid=a_grid,
            t_0=t_0,
            int_time=int_time, dt_min=dt_min, dt_max = dt_max,maxiters = 1000, atol=dp45_error_tol, rtol =dp45_error_tol)

        jac = fn.jacobian_matrix_aux(final_positions,aux_grid_spacing=aux_grid_spacing[j])
        cgst = fn.cauchy_green_tensor(jac)
        ev = np.linalg.eigvalsh(cgst)
        ev_max = np.amax(ev,-1)
        ftle = np.log(ev_max)/(2.*np.abs(int_time))
        # FTLE_list_DP45.append(ftle)
        plt.plot(np.linspace(x0,x1,res[j]),ftle[0],':', alpha=dp_line_alpha, marker='x')
        legends.append(r'DOPRI54: $\mathregular{h_{grid}}$=%.1e  , $\mathregular{\tau}$ =%.0e' %(grid_spacing[j], dp45_error_tol))



#

plt.legend(legends)
ax.text(0.8,1.02,'y = %r' %y, transform=ax.transAxes)
ax.text(0.2,1.02,'$\mathregular{t_{0}}$= %r' %(t_0), transform=ax.transAxes)
ax.text(0.4,1.02,'Resolution = %r' %actual_res, transform=ax.transAxes) # remove later

#plt.savefig('testfig23.pdf',transparent=True)
plt.show()
