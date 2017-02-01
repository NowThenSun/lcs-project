'''
Module that includes double gyre specific functions for FTLE computation.
'''
from __future__ import division
import numpy as np
import general_functions as fn

def main(amplitude, epsilon, omega, nx, ny, aux_grid_spacing, t_0, int_time, adaptive_error_tol, dt_i):
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
	# Perform rkf45 loop
	final_pos = fn.rkf45_loop_fixed_step(derivs=analytic_velocity, aux_grid=a_grid, adaptive_error_tol=adaptive_error_tol, t_0=t_0, int_time=int_time, dt_i=dt_i)
	
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
	


#main(amplitude=0.1, epsilon=0.1, omega=2*np.pi/10., nx=200, ny=100, aux_grid_spacing=1.*10**-5, t_0=0., int_time=10., adaptive_error_tol=1.*10**-4, dt_i=0.0001)

