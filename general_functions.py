'''
Module containing general functions for FTLE computation. Created for use in Project IV in Lagrangian Coherent Structures.
'''
from __future__ import division
import numpy as np

def general_params(nx, ny, aux_grid_spacing, adaptive_error_tol, t_0, int_time, dt_i):
	'''
	Function for generating general parameters for FTLE computation
	'''
	return nx, ny, aux_grid_spacing, t_0, int_time, adaptive_error_tol, dt_i

def generate_auxiliary_grid(original_grid, aux_grid_spacing):
	'''
	Function that generates auxiliary grid based on an original coordinate with 4 extra points around each original grid point.
	~~~~~~~~~
	Inputs:
	original_grid = 2*ny*nx array of original coordinate points
	aux_grid_spacing = the distance to space the auxiliary grid points away from the original points
	~~~~~~~~~
	Outputs:
	aux_grid = 2*ny*nx*4 array of auxiliary grid. [...,0] component referring to negative shift in the x-direction, [...,1] component +ve shift in the x-direction. [...,2] for -ve shift in y-direction, [...,3] for +ve shift in the y-direction
	'''
	#Make sure aux_grid_spacing is a +ve values
	h = np.abs(aux_grid_spacing)

	#Initialise aux_grid same size as original grid but with extra axis of size 4 at the end
	aux_grid = np.zeros(np.append(np.shape(original_grid),4))

	aux_grid[:] = original_grid[...,np.newaxis]

	aux_grid[0,...,0] -= aux_grid_spacing
	aux_grid[0,...,1] += aux_grid_spacing
	aux_grid[1,...,2] -= aux_grid_spacing
	aux_grid[1,...,3] += aux_grid_spacing

	return aux_grid

def rk4(y, time, dt, derivs):
	'''
	Function that moves value of y forward by a single step of size dt by the 4th order Runge-Kutta method.
	y is the whole coord grid in this case
	'''
	k0 = dt*np.array(derivs(y, time))
	k1 = dt*np.array(derivs(y + 0.5*k0, time + 0.5*dt))
	k2 = dt*np.array(derivs(y + 0.5*k1, time + 0.5*dt))
	k3 = dt*np.array(derivs(y + k2, time + dt))
	y_next = y + (k0+2.*k1+2.*k2+k3)/6.
	return y_next

def rk4_loop(derivs, aux_grid, t_0, int_time, dt,return_data=False):
	'''
	Function that performs the rk4 loop over an auxiliary grid.
	~~~~~~~~~
	Inputs:
	derivs = function that returns the derivatives of the positions (velocity)
	aux_grid = 2*ny*nx(*4) initial grid of coordinatesd
	t_0 = initial time before timestepping begins
	int_time = time integrated over
	dt = fixed timestep
	~~~~~~~~~
	Outputs:
	positions = final array of positions
	'''
	# Initialise positions array to aux. grid
	positions = aux_grid
	# Initialise time
	time = t_0

	if return_data==True:
		np.append(np.shape(positions),np.int(np.abs(int_time/dt)))
		full_positions = np.zeros((np.append(np.shape(positions),np.int(np.abs(int_time/dt)))))
		for k in xrange(np.int(np.abs(int_time/dt))):
			if np.abs(time - t_0<np.abs(int_time)):   # want to generalise so negative dt works too - hence the use of mod(k*dt) < mod(intT)
				positions = rk4(positions, time, dt, derivs)
				time += dt
				full_positions[...,k] = positions
				print k, "current time =", np.average(time), "~~~~~~~~dt =", np.average(dt)

		return positions, full_positions

	for k in xrange(np.int(np.abs(int_time/dt))):
		if np.abs(time - t_0<np.abs(int_time)):   # want to generalise so negative dt works too - hence the use of mod(k*dt) < mod(intT)
			positions = rk4(positions, time, dt, derivs)
			time += dt
			print k, "current time =", np.average(time), "~~~~~~~~dt =", np.average(dt)


	return positions

# def rk4_loop_alt(derivs, aux_grid, t_0, int_time, dt, dt_i):
	# '''
	# Quick alternative RK4 test code (don't use for -ve integration times)
	# '''
	# positions = aux_grid
	# for k in xrange(np.int(int_time/dt)+1):
		# positions = rk4(positions, t_0 + k*dt, dt, derivs)

	# return positions

def rkf45(y, time, dt, derivs, adaptive_error_tol):
	'''
	Function that moves value of y forward by a single step of size dt by the 4/5th order Runge-Kutta-Fehlberg method that uses adaptive step sizing.
	y is the whole coord grid in this case
	etol : the error tolerance in dt
	'''
	# Initialise Runge-Kutta-Fehlberg coefficients
	a2 = [1/4.]
	a3 = [3/32. , 9/32. ]
	a4 = [1932/2197. , -7200/2197. , 7296/2197. ]
	a5 = [439/216. , -8. , 3680/513. , -845/4104. ]
	a6 = [-8/27. , 2. , -3544/2565. , 1859/4104. , -11/40. ]

	k1 = dt*np.array(derivs(y, time))
	#print np.shape(k1), np.shape(dt)
	k2 = dt*np.array(derivs(y + a2*k1, time + np.sum(a2)*dt))
	k3 = dt*np.array(derivs(y + a3[0]*k1 + a3[1]*k2, time + np.sum(a3)*dt))
	k4 = dt*np.array(derivs(y + a4[0]*k1 + a4[1]*k2 + a4[2]*k3 ,time + np.sum(a4)*dt))
	k5 = dt*np.array(derivs(y + a5[0]*k1 + a5[1]*k2 + a5[2]*k3 + a5[3]*k4, time + np.sum(a5)*dt))
	k6 = dt*np.array(derivs(y + a6[0]*k1 + a6[1]*k2 + a6[2]*k3 + a6[3]*k4 + a6[4]*k5, time + np.sum(a6)*dt))

	#4th order method
	y_next = y + 25/216.*k1 + 1408/2565.*k3 + 2197/4101.*k4 - 1/5.*k5
	#print "Y_NEXT IS:", y_next-y
	#5th order method
	z_next = y  + 16/135.*k1 + 6656/12825.*k3 + 28561/56430.*k4 - 9/50.*k5 + 2/55.*k6
	#print "Z_NEXT IS:", z_next
	#scaling factor s (for dt)
	#s = np.ones(np.shape(dt))
	#If y_next = z_next (e.g. for dt = 0) replace s with s = 0
	# condition = (y_next[0]==z_next[0])*(y_next[1]==z_next[1]) # Condition where both components of y_next and z_next are 0
	# s = np.where(condition, 0., s)
	# new_s = (adaptive_error_tol/(2.*np.sqrt((z_next[0]-y_next[0])**2+(z_next[1]-y_next[1])**2)))**0.25
	# s = np.where(s==1, new_s, s)
	#print y_next - z_next
	#Alternative strategy for s
	s = np.zeros_like(dt)

	s =(adaptive_error_tol*dt/(2.*
	np.sqrt((z_next[0]-y_next[0])**2+(z_next[1]-y_next[1])**2))
	)**0.25
	s = np.where(dt==0, 0., s)

	#print np.abs(z_next-y_next)
	#print "s values......", s

	return y_next, s

def rkf45_loop_fixed_step(derivs, aux_grid, adaptive_error_tol, t_0, int_time, dt_min, dt_max):
	'''
	Function that performs the looping/timestepping part of the RKF45 scheme for an auxiliary grid with a fixed timestep for ALL points in the grid
	~~~~~~~~~
	Inputs:
	derivs = function that returns the derivatives of the positions (velocity)
	aux_grid = 2*ny*nx(*4) initial grid of coordinates
	adaptive_error_tol = error tolerance in rkf45 method
	t_0 = initial time before timestepping begins
	int_time = time integrated over
	dt_min = impose minimum timestep allowed for adaptive method to choose
	dt_max = impose maximum timestep allowed for adaptive method to choose
	~~~~~~~~~
	Outputs:
	positions = final array of positions
	'''
	# Initialise scalar of dt and times
	time = t_0	# Doesn't have to be an array for time steps that don't vary spatially
	print "timeshape", np.shape(time)
	# set initial dt -- use a scalar dt for fixed step rkf45
	dt = dt_min

	# Set initial position to aux. grid
	positions = np.zeros_like(aux_grid)
	positions[:] = aux_grid[:]

	for k in xrange(np.int(np.abs(int_time/dt_min))):
		# Condition that if k*dt = int_time the integration is complete
		if np.any(np.isclose(np.abs(int_time),np.abs(time-t_0))):
			break

		# Condition that if integration time left is less than a full step of dt a new dt is defined so the integration ends exactly after stepping through int_time
		elif np.any(np.abs(int_time)-np.abs(time-t_0)<np.abs(dt)):
			dt_final = np.sign(dt)*(np.abs(int_time)-np.abs(time-t_0))
			#print "final dt step:", dt_final
			step = rkf45(positions, time, dt_final, derivs, adaptive_error_tol = adaptive_error_tol)
			#dt = step[1]
			positions = step[0]
			time += dt_final
			#print "final time (not to be evaluated at):", np.average(time)
			break

		# If neither of above conditions met code will do the normal time-stepping method below
		else:
			step = rkf45(positions, time, dt, derivs, adaptive_error_tol = adaptive_error_tol)
			#print np.shape(step[0]),np.shape(step[1])
			# Update time to be used in the next step (the dt used in rkf45 has been used already for time)
			positions = step[0] #positions = rk4(positions,time, dt,gyre_analytic)
			time += dt

			# Update dt
			dt *= np.min(step[1])
			# Set dt = dt_min if dt too small
			if np.abs(dt) < np.abs(dt_min): # assuming scalar dt here
				dt = dt_min

			# Set = dt_max if dt too large
			if np.abs(dt) > np.abs(dt_max):
				dt = dt_max


			print k, "average time =", np.average(time), "~~~~~~~~dt =", np.average(dt)


	return positions



def rkf45_loop(derivs, aux_grid, adaptive_error_tol, t_0, int_time, dt_min, dt_max):
	'''
	Function that performs the looping/timestepping part of the RKF45 scheme for an auxiliary grid with timesteps varying for different points in the grid.
	~~~~~~~~~
	Inputs:
	derivs = function that returns the derivatives of the positions (velocity)
	aux_grid = 2*ny*nx(*4) initial grid of coordinates
	adaptive_error_tol = error tolerance in rkf45 method
	t_0 = initial time before timestepping begins
	int_time = time integrated over
	dt_min = impose minimum timestep allowed for adaptive method to choose
	dt_max = impose maximum timestep allowed for adaptive method to choose
	~~~~~~~~~
	Outputs:
	positions = final array of positions
	'''
	print "RKF45 LOOP INITIATED"
	# Initialise scalar of dt and times
	time = np.zeros((np.shape(aux_grid)[1],np.shape(aux_grid)[2],np.shape(aux_grid)[3]))+t_0	# ny*nx*4 array of times
	print "timeshape", np.shape(time)
	# set initial dt -- use an array of dt's the same size as the grid
	dt = np.zeros_like(time) + dt_min

	# Set initial position to aux. grid
	positions = np.zeros_like(aux_grid)
	positions[:] = aux_grid[:]

	for k in xrange(np.int(np.abs(int_time/dt_min))):
		# Condition that if for all grid points k*dt = int_time the integration is complete
		if np.all(np.isclose(np.abs(int_time),np.abs(time-t_0))):
			break

		# Perform rkf45 stepping until all grid points reach the total integration time
		# Use various conditions to keep dt bounded and not overstep past the final integration time
		else:
			#Perform step through dt first
			step = rkf45(positions, time, dt, derivs, adaptive_error_tol)
			# Update time used for the NEXT step
			time += dt
			positions = step[0] # Update position

			# Work out different conditions to impose onto dt
			time_left = (int_time)-(time-t_0)
			# Modify dt for each point in the grid by the recommended scalar s
			dt *= step[1]
			# Set dt = dt_min if dt too small
			dt = np.where(np.abs(dt)<np.abs(dt_min), dt_min, dt)
			# Set = dt_max if dt too large
			dt = np.where(np.abs(dt)>np.abs(dt_max), dt_max, dt)
			# Finally set dt values which are less than a full step away to the size left
			dt = np.where(np.abs(dt)>np.abs(time_left), time_left, dt)
			# Note order that these conditions are taken in is important
			# If time_left = 0, dt = 0 which keeps the time fixed at the final time once it reaches it

			print np.count_nonzero(time-t_0!=int_time)
			print k, "average time =", np.average(time), "~~~~~~~~dt =", np.average(dt)
			#break
			#print positions-aux_grid
	return positions




def jacobian_matrix_aux(aux_grid_f, aux_grid_spacing):
	'''
	Function that calculates the jacobian matrix for the flow of an auxiliary grid of points using a 2nd order central finite differences method.
	~~~~~~~~~
	Inputs:
	aux_grid_f = 2*ny*nx*4 array of final auxiliary grid points
	aux_grid_spacing = spacing of the initial auxiliary grid from original grid
	(normally would use initial grid to calculate Jacobian - but this is not necessary for auxiliary grids)
	~~~~~~~~~
	Outputs:
	jacobian = ny*nx*2*2 array representing the 2x2 Jacobian matrix of the flow at each point of the original grid
	'''
	# Initialise Jacobian matrix shape ny*nx*2*2
	jacobian = np.zeros(np.shape(aux_grid_f)[1:])
	jacobian = np.reshape(jacobian, np.append(np.shape(aux_grid_f)[1:3],(2,2)))

	# Finite differences to calculate the values of the matrix
	jacobian[...,0,0] = (aux_grid_f[0,...,1] - aux_grid_f[0,...,0])	# Variation in x wrt x
	jacobian[...,0,1] = (aux_grid_f[0,...,3] - aux_grid_f[0,...,2])	# Variation in x wrt y
	jacobian[...,1,0] = (aux_grid_f[1,...,1] - aux_grid_f[1,...,0])	# Variation in y wrt x
	jacobian[...,1,1] = (aux_grid_f[1,...,3] - aux_grid_f[1,...,2])	# Variation in y wrt y
	jacobian /= 2.*aux_grid_spacing		# Central differences factor same here since aux. grid spacing is the same in x and y

	return jacobian

def cauchy_green_tensor(jac):
	'''
	Function that takes in the Jacobian flow matrix and outputs the Cauchy-Green strain tensor (using CG = J^T * J).
	~~~~~~~~~~
	Inputs:
	jacobian =  (ny*nx*)2*2 Jacobian matrix
	~~~~~~~~~~
	Outputs:
	cgt = (ny*nx*)2*2 Cauchy-Green strain tensor
	'''
	# Initialise cgt array
	cgt = np.zeros_like(jac)
	# Input values using symmetry of Cauchy-Green tensor
	cgt[...,0,0] = jac[...,0,0]**2 + jac[...,1,0]**2
	cgt[...,1,0] = jac[...,0,0]*jac[...,0,1] + jac[...,1,0]*jac[...,1,1]
	cgt[...,0,	1] = jac[...,0,0]*jac[...,0,1] + jac[...,1,0]*jac[...,1,1]
	cgt[...,1,1] = jac[...,0,1]**2 + jac[...,1,1]**2

	return cgt
