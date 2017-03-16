'''
Module containing methods for integrating ordinary differential equations
'''
from __future__ import division
import numpy as np


	# h_try - trial step
	# Take a step
	# Evaluate accuracy
	# Step rejected - try again with reduced h given by step controller

	# Step success
	# Reuse last derivative for next step (First Same As Last)

def dp45(y, time, dt, derivs, atol, rtol):
	'''
	Subroutine that calculates one step of Dormand-Prince 4-5 method.
	'''
	# Initialise Runge-Kutta-Fehlberg coefficients
	a2 = [1/5]
	a3 = [3/40. , 9/40. ]
	a4 = [44/45. , -56/15. , 32/9. ]
	a5 = [19372/6561. , -25360/2187. , 64448/6561. , -212/729. ]
	a6 = [9017/3168. , -355/33. , 46732/5247. , 49/176. , -5103/18656. ]
	a7 = [ 35/384. , 0. , 500/1113. , 125/192. , -2187/6784. , 11/84.]
	#4th order method coefficients
	b1, b2, b3, b4, b5, b6, b7 = (5179/57600. , 0. , 7571/16695. , 393/640. , -92097/339200. , 187/2100. , 1/40.)
	#5th order methd coefficients
	B1, B2, B3, B4, B5, B6, B7 = (35/384. , 0. , 500/1113. , 125/192. , -2187/6784. , 11/84. , 0.)

	k1 = dt*np.array(derivs(y, time))
	k2 = dt*np.array(derivs(y + a2*k1, time + np.sum(a2)*dt))
	k3 = dt*np.array(derivs(y + a3[0]*k1 + a3[1]*k2, time + np.sum(a3)*dt))
	k4 = dt*np.array(derivs(y + a4[0]*k1 + a4[1]*k2 + a4[2]*k3 ,time + np.sum(a4)*dt))
	k5 = dt*np.array(derivs(y + a5[0]*k1 + a5[1]*k2 + a5[2]*k3 + a5[3]*k4, time + np.sum(a5)*dt))
	k6 = dt*np.array(derivs(y + a6[0]*k1 + a6[1]*k2 + a6[2]*k3 + a6[3]*k4 + a6[4]*k5, time + np.sum(a6)*dt))
	k7 = dt*np.array(derivs(y + a7[0]*k1 + a7[1]*k2 + a7[2]*k3 + a7[3]*k4 + a7[4]*k5 + a7[5]*k6, time + np.sum(a7)*dt)) # First same as last (FSAL) step in DP
	# So k7 is the 5th order estimate for the next step in y
	# 4th order estimate for y_next
	Y_next = y + b1*k1 + b2*k2 + b3*k3 + b4*k4 + b5*k5 + b6*k6 + b7*k7
	y_next = y + B1*k1 + B2*k2 + B3*k3 + B4*k4 + B5*k5 + B6*k6 + B7*k7
	
	delta = y_next - Y_next
	scale = atol + np.maximum(np.abs(y),np.abs(y_next))*rtol
	err = np.sqrt(0.5*((delta/scale)[0]**2+(delta/scale)[1]**2))
	print "shape of error:", np.shape(err)
	#s =(atol*dt/(2.*delta))**0.25
	safety_factor = 0.97
	s = np.where(dt==0, 0., safety_factor*(1./err)**0.2)

	return y_next, s, err


def dp45_loop(derivs, aux_grid, t_0, int_time, dt_min, dt_max, maxiters, atol, rtol):
	'''
	Function that performs the looping/timestepping part of the RKF45 scheme for an auxiliary grid with timesteps varying for different points in the grid.
	~~~~~~~~~
	Parameters:				Inputs:
	derivs 					function that returns the derivatives of the positions (velocity)
	aux_grid  				2*ny*nx(*4) initial grid of coordinates (y(t0) = y0 = aux_grid)
	t_0 					initial time
	int_time 				time integrated over
	dt_min 					impose minimum timestep allowed for adaptive method to choose
	dt_max 					impose maximum timestep allowed for adaptive method to choose
	maxiters				Maximum number of iterations
	atol 					absolute error tolerance
	rtol 					relative error tolerance (scaled by position values)
	~~~~~~~~~
	Outputs:
	pos 						final array of positions
	flag						error flag
								0   no errors
								1   hmin exceeded (not implemented at the moment)
								2   maximum iterations exceeded
	'''
	print "DP45 LOOP INITIATED"
	# Initialise scalar of dt and times
	time = np.zeros((np.shape(aux_grid)[1],np.shape(aux_grid)[2],np.shape(aux_grid)[3]))+t_0	# ny*nx*4 array of times
	print "timeshape", np.shape(time)
	# set initial dt -- use an array of dt's the same size as the grid
	dt = np.zeros_like(time) + dt_min

	# Set initial position to aux. grid
	pos = np.zeros_like(aux_grid)
	pos[:] = aux_grid[:]


	for k in xrange(maxiters):
		# Condition that if for all grid points k*dt = int_time the integration is complete
		if np.all(np.isclose(np.abs(int_time),np.abs(time-t_0))):
			flag = 0
			break

		# Perform rkf45 stepping until all grid points reach the total integration time
		# Use various conditions to keep dt bounded and not overstep past the final integration time
		else:
			#Perform step through dt first
			pos_new, s, err = dp45(pos, time, dt, derivs, atol, rtol)
			# Total Error tolerance
			# If delta < tol step is successful so update times and positions
			time = np.where(err<1,time+dt,time)
			pos = np.where(err<1,pos_new,pos)
			# Keep s restricted between 0.1-4 of the previous value
			s[s<0.1] = 0.1
			s[s>4] = 4

			# Work out different conditions to impose onto dt
			time_left = (int_time)-(time-t_0)
			# Modify dt for each point in the grid by the recommended scalar s
			dt *= s
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
	return pos





def rkf45(y, time, dt, derivs, atol, rtol):
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
	#4th order method coefficients
	b1, b2, b3, b4, b5, b6, b7 = (25/216, 0, 1408/2565, 2197/4104, -1/5, 0, 0)
	#5th order methd coefficients
	B1, B2, B3, B4, B5, B6, B7 = (16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55, 0)


	k1 = dt*np.array(derivs(y, time))
	#print np.shape(k1), np.shape(dt)
	k2 = dt*np.array(derivs(y + a2*k1, time + np.sum(a2)*dt))
	k3 = dt*np.array(derivs(y + a3[0]*k1 + a3[1]*k2, time + np.sum(a3)*dt))
	k4 = dt*np.array(derivs(y + a4[0]*k1 + a4[1]*k2 + a4[2]*k3 ,time + np.sum(a4)*dt))
	k5 = dt*np.array(derivs(y + a5[0]*k1 + a5[1]*k2 + a5[2]*k3 + a5[3]*k4, time + np.sum(a5)*dt))
	k6 = dt*np.array(derivs(y + a6[0]*k1 + a6[1]*k2 + a6[2]*k3 + a6[3]*k4 + a6[4]*k5, time + np.sum(a6)*dt))
	k7 = 0
	#4th order method
	y_next = y + b1*k1 + b2*k2 + b3*k3 + b4*k4 + b5*k5 + b6*k6 + b7*k7
	#5th order method
	Y_next = y + B1*k1 + B2*k2 + B3*k3 + B4*k4 + B5*k5 + B6*k6 + B7*k7
	# difference between different orderered solutions - taking the Euclidean norm
	delta = (np.sqrt( ((b1-B1)*k1 + (b2-B2)*k2 + (b3-B3)*k3 + (b4-B4)*k4 + (b5-B5)*k5 + (b6-B6)*k6 + (b7-B7)*k7)[0]**2
	+((b1-B1)*k1 + (b2-B2)*k2 + (b3-B3)*k3 + (b4-B4)*k4 + (b5-B5)*k5 + (b6-B6)*k6 + (b7-B7)*k7)[1]**2 ))
	# full error tolerance
	scale = atol + np.sqrt(np.maximum(y,y_next)[0]**2+np.maximum(y,y_next)[1]**2)*rtol

	# Calculate scaling factor s for the stepsize
	s = np.zeros_like(dt)
	#s =(atol*dt/(2.*delta))**0.25
	s = np.where(dt==0, 0., 0.84*(scale/delta)**0.2)


	return y_next, s, delta, scale




def rkf45_loop(derivs, aux_grid, t_0, int_time, dt_min, dt_max, maxiters, atol, rtol):
	'''
	Function that performs the looping/timestepping part of the RKF45 scheme for an auxiliary grid with timesteps varying for different points in the grid.
	~~~~~~~~~
	Parameters:				Inputs:
	derivs 					function that returns the derivatives of the positions (velocity)
	aux_grid  				2*ny*nx(*4) initial grid of coordinates (y(t0) = y0 = aux_grid)
	t_0 					initial time
	int_time 				time integrated over
	dt_min 					impose minimum timestep allowed for adaptive method to choose
	dt_max 					impose maximum timestep allowed for adaptive method to choose
	maxiters				Maximum number of iterations
	atol 					absolute error tolerance
	rtol 					relative error tolerance (scaled by position values)
	~~~~~~~~~
	Outputs:
	pos 						final array of positions
	flag						error flag
								0   no errors
								1   hmin exceeded (not implemented at the moment)
								2   maximum iterations exceeded
	'''
	print "RKF45 LOOP INITIATED"
	# Initialise scalar of dt and times
	time = np.zeros((np.shape(aux_grid)[1],np.shape(aux_grid)[2],np.shape(aux_grid)[3]))+t_0	# ny*nx*4 array of times
	print "timeshape", np.shape(time)
	# set initial dt -- use an array of dt's the same size as the grid
	dt = np.zeros_like(time) + dt_min

	# Set initial position to aux. grid
	pos = np.zeros_like(aux_grid)
	pos[:] = aux_grid[:]


	for k in xrange(maxiters):
		# Condition that if for all grid points k*dt = int_time the integration is complete
		if np.all(np.isclose(np.abs(int_time),np.abs(time-t_0))):
			flag = 0
			break

		# Perform rkf45 stepping until all grid points reach the total integration time
		# Use various conditions to keep dt bounded and not overstep past the final integration time
		else:
			#Perform step through dt first
			pos_new, s, delta,scale = rkf45(pos, time, dt, derivs, atol, rtol)
			# Total Error tolerance
			#scale = atol + np.sqrt(np.maximum(pos,pos_new)[0]**2+np.maximum(pos,pos_new)[1]**2)*rtol
			# If delta < tol step is successful so update times and positions
			time = np.where(delta<scale,time+dt,time)
			pos = np.where(delta<scale,pos_new,pos)
			# time = np.where(err<1,time+dt,time)
			# pos = np.where(err<1,pos_new,pos)
			# Keep s restricted between 0.1-4 of the previous value
			s[s<0.1] = 0.1
			s[s>4] = 4

			# Work out different conditions to impose onto dt
			time_left = (int_time)-(time-t_0)
			# Modify dt for each point in the grid by the recommended scalar s
			dt *= s
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
	return pos
































#
#
# def rkf45_loop(derivs, aux_grid, t0, int_time, dt_min, dt_max, adaptive_error_tol, maxiters):
# 	'''
# 	Function that performs the looping/timestepping part of the RKF45 scheme for an auxiliary grid with timesteps varying for different points in the grid.
# 	~~~~~~~~~
# 	Parameters:				Inputs:
# 	derivs 					function that returns the derivatives of the positions (velocity)
# 	aux_grid  				2*ny*nx(*4) initial grid of coordinates (y(t0) = y0 = aux_grid)
# 	t0 						initial time
# 	int_time 				time integrated over
# 	dt_min 					impose minimum timestep allowed for adaptive method to choose
# 	dt_max 					impose maximum timestep allowed for adaptive method to choose
# 	adaptive_error_tol 		error tolerance in rkf45 method
# 	maxiters				Maximum number of iterations
# 	~~~~~~~~~
# 	Outputs:
# 	positions 				 final array of positions
# 	flag						error flag
# 								0   no errors
# 								1   hmin exceeded
# 								2   maximum iterations exceeded
# 	'''
# 	print "RKF45 LOOP INITIATED"
# 	# Initialise scalar of dt and times
# 	time = np.zeros((np.shape(aux_grid)[1],np.shape(aux_grid)[2],np.shape(aux_grid)[3]))+t_0	# ny*nx*4 array of times
# 	#print "timeshape", np.shape(time)
# 	# set initial dt -- use an array of dt's the same size as the grid
# 	dt = np.zeros_like(time) + dt_min
#
# 	# Set initial position to aux. grid
# 	positions = np.zeros_like(aux_grid)
# 	positions[:] = aux_grid[:]
# 	for i in xrange(maxiters):
# 		rkf45()
# 		#Check error size --- base this on an absolute tolerance and relative tolerance
# 		#If error <= 1 : Step succeeded
# 			#Compute dt_next
#
# 			#Keep scaling factor within some bounds of previous scaling factor s
#
#
#
#
# 		if maxiters - i == 0:
# 			#raise flag of maximum iterations reached
# 			flag = 2
# 	return positions
