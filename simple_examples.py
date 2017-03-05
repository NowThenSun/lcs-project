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

def ex2(coordinates,time_array):
    '''
    Function that calculates the analytic velocity for simple time-independent analytic examples for an array of coordinates
    ~~~~~~~~~~
    Inputs:
    coordinates = 2*ny*nx(*x) sized array of input coordinates (where x can be any integer)
    ~~~~~~~~~~
    Outputs:
    u, v = analytic velocity fields (same shape as the input coordinates)
    '''
    u = 1 + (np.tanh(coordinates[0]))**2
    v = -coordinates[1]
    return u, v  #Note this returns (u,v) as a tuple


nx = 7
ny = 7
xlower = -2.
xupper = 2.
ylower = -2
yupper = 2
velocity_example = ex2

X = np.linspace(xlower+1/(2.*nx), xupper-1/(2.*nx), nx)
Y = np.linspace(ylower + 1/(2.*ny), yupper-1/(2.*ny), ny)
coords = np.meshgrid(X, Y, indexing='xy')
U,V = velocity_example(coords,"doesnt matter")
# Want to seed initial points
#Define a dense fluid square blob
density_blob=100

diam_blob = 0.5
#corner_pos_blob = (-diam_blob/2,0.85)#
corner_pos_blob=  (-2, -diam_blob/2)
centre_pos_blob = (corner_pos_blob[0]+diam_blob/2,corner_pos_blob[1]+diam_blob/2)

xline = np.linspace(corner_pos_blob[0], corner_pos_blob[0]+diam_blob, density_blob)
yline = np.linspace(corner_pos_blob[1], corner_pos_blob[1]+diam_blob, density_blob)
square_blob = np.meshgrid(xline, yline, indexing='xy')
circle_bool = np.zeros_like(square_blob[0])
circle_bool[((square_blob[0]-centre_pos_blob[0])**2+(square_blob[1]-centre_pos_blob[1])**2)<(diam_blob/2)**2] = 1
circle_blob = np.zeros((2,np.int(np.shape(np.nonzero(circle_bool))[1])))
circle_blob[0] = np.array(square_blob[0])[circle_bool==1.]
circle_blob[1] = np.array(square_blob[1])[circle_bool==1.]
#Advect with rk4
final_advected_pos, full_positions = fn.rk4_loop(derivs=velocity_example, aux_grid=circle_blob, t_0=0, int_time=2.4, dt=0.06,return_data=True)
print np.shape(final_advected_pos)
fig= plt.figure()
ax = plt.axes()

#plt.quiver(X,Y,U,V)
# plt.plot(square_blob[0],square_blob[1],marker='o',markeredgewidth=0, color='b', alpha=1)
print np.shape(full_positions)
# plt.plot(full_positions[0,...,10],full_positions[1,...,10],marker='o',markeredgewidth=0, color='b', alpha=1)
# plt.plot(full_positions[0,...,20].flatten(),full_positions[1,...,20].flatten(),marker='o',markeredgewidth=0, color='b', alpha=1)
# ax.fill(final_advected_pos[0].flatten(),final_advected_pos[1].flatten(),'r')
blob_color = (0.1,0.4,0.7,0.5)
# Set time indices for intermediate blobs
blob_timestep1 = 5
blob_timestep2 = 14
blob_timestep3 = 25

ax.scatter(circle_blob[0],circle_blob[1], c=blob_color,
    alpha=0.3, edgecolors='none')
ax.scatter(full_positions[0,...,blob_timestep1],full_positions[1,...,blob_timestep1], c=blob_color,
    alpha=0.3, edgecolors='none')
ax.scatter(full_positions[0,...,blob_timestep2],full_positions[1,...,blob_timestep2], c=blob_color,
    alpha=0.3, edgecolors='none')
ax.scatter(full_positions[0,...,blob_timestep3],full_positions[1,...,blob_timestep3], c=blob_color,
    alpha=0.3, edgecolors='none')
ax.scatter(full_positions[0,...,-1],full_positions[1,...,-1], c=blob_color,
    alpha=0.3, edgecolors='none')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.streamplot(X,Y,U,V, density=1, color='cornflowerblue')#, color='0.8', density=0.5)
#plt.plot(final_advected_pos[0],final_advected_pos[1],marker='o',markeredgewidth=0, color='r', alpha=1)
plt.show()