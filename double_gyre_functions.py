'''
Module that includes double gyre specific functions for FTLE computation.
'''
from __future__ import division
import numpy as np
import general_functions as fn
import ODE_solvers as ODE

def main(amplitude, epsilon, omega, nx, ny, aux_grid_spacing, t_0, int_time, adaptive_error_tol, dt_min, dt_max, method ='rkf45'):
	'''
	Main function that puts together other functions to compute the FTLE field of the double gyre
	~~~~~~~~~~
	Inputs:
	==
	double-gyre parameters:
	amplitude = amplitude of oscillations
	epsilon = amplitude of time dependent perturbations
	omega = angular frequency of double gyre system
	==
	general parameters:
	nx = # points of original grid in x-direction
	ny = # points of original grid in y-direction
	aux_grid_spacing = spacing of the auxiliary grid from the original grid used for finite differencing
	t_0 = initial time to start integration at
	int_time = time integrated over
	adaptive_error_tol = error tolerance used in the adaptive timestep integration scheme
	method = method used for numerical integration -- currently only supports rkf45 fixed
	~~~~~~~~~~
	Outputs:
	ftle = the finite-time Lyapunov exponent for the double gyre in a ny*nx array
	ev = eigenvalues of the Cauchy-Green strain tensor
	'''
	# Globalise double gyre parameters (for convenience of using analytic_velocity in rkf45 loop)
	dg_params = gyre_global_params(amplitude=amplitude, epsilon=epsilon, omega=omega)
	#~~gen_params = fn.general_params(nx=nx, ny=ny, aux_grid_spacing=aux_grid_spacing, adaptive_error_tol=adaptive_error_tol, t_0=t_0, int_time=int_time, dt_i=dt_i)

	a_grid = fn.generate_auxiliary_grid(generate_grid(nx,ny), aux_grid_spacing=aux_grid_spacing)
	#~~print np.shape(a_grid)

	if method == 'rkf45_fixed':
		# Perform rkf45 loop for fixed step on all grid points
		final_pos = fn.rkf45_loop_fixed_step(derivs=analytic_velocity, aux_grid=a_grid,
			adaptive_error_tol=adaptive_error_tol, t_0=t_0, int_time=int_time, dt_min=dt_min,dt_max= dt_max)
		print "RKF45 fixed time step method used"

	if method == 'rk4':
		# Perform rk4 loop procedure

		print "RK4 method used"

	if method == 'rkf45':
		print "OLD VERSION OF RKF45 USED WHERE NO RESTEPPING OCCURS"
		final_pos = fn.rkf45_loop(derivs=analytic_velocity, aux_grid=a_grid,
			adaptive_error_tol=adaptive_error_tol, t_0=t_0, int_time=int_time, dt_min=dt_min,dt_max= dt_max)
		print "OLD VERSION OF RKF45 USED WHERE NO RESTEPPING OCCURS"

		print "RKF45 full adaptive grid method used"

	if method == 'rkf45error':
		final_pos = ODE.rkf45_loop(derivs=analytic_velocity, aux_grid=a_grid,
			t_0=t_0, int_time=int_time, dt_min=dt_min,dt_max= dt_max,
			maxiters = 1000, atol = adaptive_error_tol, rtol = adaptive_error_tol)

	if method == 'dp45':
		final_pos = ODE.dp45_loop(derivs=analytic_velocity, aux_grid=a_grid,
			t_0=t_0, int_time=int_time, dt_min=dt_min,dt_max= dt_max,
			maxiters = 1000, atol = adaptive_error_tol, rtol = adaptive_error_tol)
	# Calculate jacobian matrix
	jac = fn.jacobian_matrix_aux(final_pos,aux_grid_spacing=aux_grid_spacing)
	cgst = fn.cauchy_green_tensor(jac)

	ev = np.linalg.eigvalsh(cgst)
	ev_max = np.amax(ev,-1)
	ftle = np.log(ev_max)/(2.*np.abs(int_time))

	return ftle

def gyre_global_params(amplitude, epsilon, omega):
	'''
	Function that returns and globalises the parameters used in the analytic double gyre example
	'''
	global amplitude_g
	global epsilon_g
	global omega_g
	amplitude_g, epsilon_g, omega_g = amplitude, epsilon, omega
	return amplitude, epsilon, omega

def generate_grid(nx, ny):
	'''
	Generates grid on the domain (x,y)=[0,2]x[0,1]
	~~~~~~~~~~~~
	Inputs:
	nx = # grid points in x
	ny = # grid points in y
	~~~~~~~~~~~~
	Outputs:
	coords = 2*ny*nx array of coordinates with [0] component returning the x-component and [1] returning the y-component
	'''

	xlower = 0.
	xupper = 2.
	ylower = 0.
	yupper = 1.

	#Generate equally spaced points just from the edge of the grid
	X = np.linspace(xlower+1/(2.*nx), xupper-1/(2.*nx), nx)
	Y = np.linspace(ylower + 1/(2.*ny), yupper-1/(2.*ny), ny)

	coords = np.meshgrid(X, Y, indexing='xy') 			# 2*ny*nx array with [0] component giving x values

	return np.array(coords)#, X_, Y_

def analytic_velocity(coordinates, time_array):
	'''
	Function that calculates the analytic velocity for the double gyre for an array of coordinates and array of times
	~~~~~~~~~~
	Inputs:
	coordinates = 2*ny*nx(*x) sized array of input coordinates (where x can be any integer)
	time_array = ny*nx(*x) sized array of times for each point in the grid
	amplitude, epsilon, omega = the double gyre parameters
	~~~~~~~~~~
	Outputs:
	u, v = analytic velocity fields of the double gyre (same shape as the input coordinates)
	'''

	a = epsilon_g*np.sin(omega_g*time_array)

	b = 1. - 2.*epsilon_g*np.sin(omega_g*time_array)

	f = a*coordinates[0]**2 + b*coordinates[0]
	# df = df/dx
	df = 2.*a*coordinates[0] + b

	u = -np.pi*amplitude_g*np.sin(np.pi*f)*np.cos(np.pi*coordinates[1])

	v = np.pi*amplitude_g*np.cos(np.pi*f)*np.sin(np.pi*coordinates[1])*df

	return u, v  #Note this returns (u,v) as a tuple


def analytic_velocity_noglobal(epsilon_loc, amplitude_loc, omega_loc, coordinates, time_array):
	'''
	Function that calculates the analytic velocity for the double gyre for an array of coordinates and array of times, without requiring global parameters
	~~~~~~~~~~
	Inputs:
	coordinates = 2*ny*nx(*x) sized array of input coordinates (where x can be any integer)
	time_array = ny*nx(*x) sized array of times for each point in the grid
	amplitude, epsilon, omega = the double gyre parameters
	~~~~~~~~~~
	Outputs:
	u, v = analytic velocity fields of the double gyre (same shape as the input coordinates)
	'''

	a = epsilon_loc*np.sin(omega_loc*time_array)

	b = 1. - 2.*epsilon_loc*np.sin(omega_loc*time_array)

	f = a*coordinates[0]**2 + b*coordinates[0]
	# df = df/dx
	df = 2.*a*coordinates[0] + b

	u = -np.pi*amplitude_loc*np.sin(np.pi*f)*np.cos(np.pi*coordinates[1])
	v = np.pi*amplitude_loc*np.cos(np.pi*f)*np.sin(np.pi*coordinates[1])*df
	return u, v  #Note this returns (u,v) as a tuple

#
# nx = 500
# ftle = main(amplitude=0.3, epsilon=0.2, omega=2*np.pi/10.,
# 	nx=nx, ny=nx/2, aux_grid_spacing=0.08*2/nx,
# 	t_0=2., int_time=100, adaptive_error_tol=1.*10**-3,
# 	dt_min=0.001, dt_max=0.2, method = 'dp45')
# 
#
# import matplotlib.pyplot as plt
# import matplotlib as mpl
# # Plotting code for plot of eigenvalues/FTLE field
# fig = plt.figure()
# ax = plt.axes()
# #mpl.cm.register_cmap(name='mycubehelix',data=mpl._cm.cubehelix(g,s,r,sat)) ##way to add a colourmap directly (from segmentdata) into mpl.cm or something
# #cmap = mpl.colors.LinearSegmentedColormap(name='abab', segmentdata =mpl._cm.cubehelix(g,s,r,sat))  ##way to create colormap from colour dictionary
# def cubehelix_cmap(g=1.0, s=0.5, r = -1.5, sat = 1.0):
#     '''
#     ~~~~~~~~~~
#     Inputs:
#     g : gamma value (can increase intensity of high valued colors for g>1, increase intensity for low values for g<1)
#     s : starting color
#     r : number of rotations through B -> G -> R
#     sat : saturation value (0 for grayscale)
#     ~~~~~~~~~~~
#     Outputs:
#     cubehelix colourmap
#     reverse cubehelix colourma	p
#     '''
#     cdict = mpl._cm.cubehelix(g,s,r,sat)
#     def rev_fn(f):
#         def reverse_f(x):
#             return f(1-x)
#         return reverse_f
#     b_r = rev_fn(cdict['blue'])
#     g_r = rev_fn(cdict['green'])
#     r_r = rev_fn(cdict['red'])
#     cdict_r = {u'blue' : b_r, u'green' : g_r, u'red' : r_r}
#     cmap = mpl.colors.LinearSegmentedColormap(name='ch', segmentdata=cdict)
#     cmap_r = mpl.colors.LinearSegmentedColormap(name='ch_r', segmentdata=cdict_r)
#     return cmap, cmap_r
#
#
# im = ax.imshow(ftle, interpolation='none', origin='lower', extent=(0,2,0,1),
#     cmap=cubehelix_cmap(g=2.0,s=-1.2,r=-0.85,sat=1.0)[1],vmax=0.5,vmin=0) #,aspect='auto' vmin=-0.0001,vmax=0.0001,
# #
# # ax.text(0.8,1.02,'T = %.1f' %int_time, transform=ax.transAxes)
# # ax.text(-0.1,1.02,'t_0 = %.1f' %t_0, transform=ax.transAxes)
# # #ax.text(0.3,1.02,'average dt = %.2e' %np.average(dt), transform=ax.transAxes)
# # ax.text(0.6,-0.17,'error tol in dt= %r' %adap_error_tol, transform=ax.transAxes)
# cbar_ax = fig.add_axes([0.900, 0.23, 0.025, 0.54])
# pos1 = ax.get_position() # get the original position
# pos2 = [pos1.x0 - 0.04, pos1.y0 ,  pos1.width , pos1.height ]
# ax.set_position(pos2) # set a new position
# #print pos1
# #cbar_ax.set_title('title',fontsize=11,y=1.02,x=1.005)
# #ax1.text(0.8,0.9,r'$t$ = %d $\mu$s' %t[T],fontsize=13,transform=ax1.transAxes, color='Azure')
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# cb = fig.colorbar(im, cax=cbar_ax ) #ticks=np.arange(8.)/10.
# # plt.savefig('DG_FTLE_-1A_-15e_15T_5t0_1000x500_v4.pdf',transparent=True,bbox_inches='tight')
# plt.show()
