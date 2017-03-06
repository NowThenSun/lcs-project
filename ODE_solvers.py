'''
Module containing methods for integrating ordinary differential equations
'''
from __future__ import division
import numpy as np

def dp45():
	'''
	Subroutine that calculates one step of Dormand-Prince 4-5 method.
	'''



# h_try - trial step
# Take a step
# Evaluate accuracy
# Step rejected - try again with reduced h given by step controller

# Step success
# Reuse last derivative for next step (First Same As Last)

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

	#4th order method
	y_next = y + b1*k1 + b2*k2 + b3*k3 + b4*k4 +b5*k5 +b6*k6 + b7*k7
	# difference between different orderered solutions
	delta =np.abs ((b1-B1)*k1 + (b2-B2)*k2 + (b3-B3)*k3 + (b4-B4)*k4 + (b5-B5)*k5 + (b6-B6)*k6 + (b7-B7)*k7)

	return y_next, delta
