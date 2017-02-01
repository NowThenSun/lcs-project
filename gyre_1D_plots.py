'''
Module that creates 1D plots over ridges in the double gyre field to compare errors in peaks at different resolutions, methods and different time steps
compare to RKF45 that alternates all of the timesteps
Can also compare dgyre analytic to different interpolation functions
'''
from __future__ import division
import numpy as np
import general_functions as fn
import double_gyre_functions as dg

def generate_1D_grid():
	pass

#~~~~~ 1D PLOTTING PARAMETERS ~~~~~
y = 0.3  # fixed y-value the FTLE is evaluated at
x0 = 1.0
x1 = 1.2
res = np.array([50,300,1000]) # Resolutions tested at (number of grid spacing)               0.2/100 ~ 1000x500 res
actual_res = res*2./(x1-x0)

#list to put in final data of calculated FTLE fields
FTLE_list = [] 

for j in xrange(len(res)):
	dg.main(amplitude=0.1, epsilon=0.1, omega=2*np.pi/10., nx=200, ny=100, aux_grid_spacing=1.*10**-5, t_0=0., int_time=10., adaptive_error_tol=1.*10**-4, dt_i=0.0001)