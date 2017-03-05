'''
Module for simple analytic examples of dynamical systems
'''
from __future__ import division
import numpy as np
import general_functions as fn
import matplotlib.pyplot as plt

def ex1(coordinates,time_array):
    '''
    Function that calculates the analytic velocity for simple time-independent analytic examples for an array of coordinates
    ~~~~~~~~~~
    Inputs:
    coordinates = 2*ny*nx(*x) sized array of input coordinates (where x can be any integer)
    ~~~~~~~~~~
    Outputs:
    u, v = analytic velocity fields (same shape as the input coordinates)
    '''
    u = coordinates[0]
    v = -coordinates[1]-coordinates[1]**3
    return u, v  #Note this returns (u,v) as a tuple

nx = 7
ny = 7
xlower = -1.
xupper = 1.
ylower = -1
yupper = 1

X = np.linspace(xlower+1/(2.*nx), xupper-1/(2.*nx), nx)
Y = np.linspace(ylower + 1/(2.*ny), yupper-1/(2.*ny), ny)
coords = np.meshgrid(X, Y, indexing='xy')
U,V = ex1(coords,"doesnt matter")
# Want to seed initial points
#Define a dense fluid square blob
density_blob=200
corner_pos_blob = (-.05,0.9)
centre_pos_blob = (corner_pos_blob[0]+diam_blob/2,corner_pos_blob[1]+diam_blob/2)
diam_blob = 0.1
xline = np.linspace(corner_pos_blob[0], corner_pos_blob[0]+diam_blob, density_blob)
yline = np.linspace(corner_pos_blob[1], corner_pos_blob[1]+diam_blob, density_blob)
square_blob = np.meshgrid(xline, yline, indexing='xy')

#Advect with rk4
advected_pos = fn.rk4_loop(derivs=ex1, aux_grid=square_blob, t_0=0, int_time=3, dt=0.1)
print np.shape(advected_pos)
plt.figure()
plt.streamplot(X,Y,U,V)
plt.quiver(X,Y,U,V)
plt.plot(square_blob[0],square_blob[1],marker='o',markeredgewidth=0, color='b', alpha=1)
plt.plot(advected_pos[0],advected_pos[1],marker='o',markeredgewidth=0, color='r', alpha=1)
plt.show()
