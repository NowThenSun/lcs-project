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

def ex3(coordinates,time_array):
    '''
    Function that calculates the analytic velocity for simple time-independent analytic examples for an array of coordinates
    ~~~~~~~~~~
    Inputs:
    coordinates = 2*ny*nx(*x) sized array of input coordinates (where x can be any integer)
    ~~~~~~~~~~
    Outputs:
    u, v = analytic velocity fields (same shape as the input coordinates)
    '''
    u = coordinates[0] - coordinates[0]**3
    v = -coordinates[1]
    return u, v  #Note this returns (u,v) as a tuple

def ex4(coordinates,time_array):
    u = coordinates[0] + 2*coordinates[0]**3
    v = -coordinates[1]
    return u, v  #Note this returns (u,v) as a tuple


def plot_fluid_parcel(nx = 7, ny = 7, xlower = -2., xupper = 2.,ylower = -2, yupper = 2,
    velocity_example = ex2,
    density_blob=100,diam_blob = 1.5 , corner_pos_blob=  (-3, -1.5/2),
    blob_timestep1 = -1, blob_timestep2 = -1, blob_timestep3 = 20,
    t_0=0, int_time=2.4, dt=0.06
    ):
    '''
    Plotting code to plot streamlines and a fluid parcel at several different points.
    The fluid parcel is advected using Runge-Kutta 4th order method
    '''
    X = np.linspace(xlower, xupper, nx)
    Y = np.linspace(ylower, yupper, ny)
    coords = np.meshgrid(X, Y, indexing='xy')
    U,V = velocity_example(coords,"doesnt matter")
    # Want to seed initial points
    #Define a dense fluid square blob



    #corner_pos_blob = (-diam_blob/2,0.85)#

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
    final_advected_pos, full_positions = fn.rk4_loop(derivs=velocity_example, aux_grid=circle_blob, t_0=t_0, int_time=int_time, dt=dt,return_data=True)
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


# plot_fluid_parcel(diam_blob=1.2, corner_pos_blob=(-2,-0.6))  # Default parameters for example2
#
# plot_fluid_parcel(nx = 20, ny = 20, xlower = -1.3, xupper = 1.3,ylower = -1.3, yupper = 1.3,
#     velocity_example = ex1,
#     density_blob=200,diam_blob = 0.2 , corner_pos_blob=  (-0.1, 0.85),
#     blob_timestep1 = 5, blob_timestep2 = 10, blob_timestep3 = -1,
#     t_0=0, int_time=2.4, dt=0.06
#     )  # parameters for example1

# plot_fluid_parcel(nx = 20, ny = 20, xlower = -2.5, xupper = 2.5,ylower = -2.5, yupper = 2.5,
#     velocity_example = ex3,
#     density_blob=200,diam_blob = 0.2 , corner_pos_blob=  (-0.5, 0.85),
#     blob_timestep1 = 5, blob_timestep2 = 15, blob_timestep3 = -1,
#     t_0=0, int_time=2.4, dt=0.06
#     )  # parameters for example1


def plot_streamlines(nx = 7, ny = 7, xlower = -2., xupper = 2.,ylower = -2, yupper = 2,
    velocity_example = ex2,
    ):
    '''
    Plotting code to plot streamlines and dynamical systems like plots.

    '''
    X = np.linspace(xlower, xupper, nx)
    Y = np.linspace(ylower, yupper, ny)
    coords = np.meshgrid(X, Y, indexing='xy')
    U,V = velocity_example(coords,"doesnt matter (in autonomous case)")
    # Want to seed initial points
    #Define a dense fluid square blob


    fig= plt.figure()
    ax = plt.axes()
    #start_points =  [[0,0],[1,1],[1,-1],[1,0],[-1,0],[-1,1],[-1,-1],[0,1],[0,-1]]
    ax.streamplot(X,Y,U,V, color='cornflowerblue', density=0.5)


    mlw = 1.5
    hl = 0.06
    hw = 0.1
    oh = 0.3


    plt.axhline(y=0,xmin=0,xmax=1, lw= mlw, c='Black' )
    plt.axvline(x=0, lw = mlw, c= 'Black')
    # plt.axvline(x=1, lw = mlw, c= 'Black')
    # plt.axvline(x=-1,lw = mlw, c= 'Black')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.xlim(xlower,xupper)
    plt.ylim(ylower,yupper)

    plt.xticks([0.])
    plt.yticks([0.])

    plt.arrow(0.0,0,0.5,0.0, lw = mlw, head_width = hw, head_length=hl, overhang=oh, ec='none', fill=True, fc='Black', ls='--', zorder=3)
    plt.arrow(0.0,0,-0.5,0.0, lw = mlw, head_width = hw, head_length=hl, overhang=oh, ec='none', fill=True, fc='Black', ls='--', zorder=3)
    # plt.arrow(-2,0,0.5,0.0, lw = mlw, head_width = hw, head_length=hl, overhang=oh, ec='none', fill=True, fc='Black', ls='--', zorder=3)
    # plt.arrow(2,0,-0.5,0.0, lw = mlw, head_width = hw, head_length=hl, overhang=oh, ec='none', fill=True, fc='Black', ls='--', zorder=3)

    hw = 0.06
    hl = 0.1
    # plt.arrow(1.,1.,0,-0.5, lw = mlw, head_width = hw, head_length=hl, overhang=oh, ec='none', fill=True, fc='Black', ls='--', zorder=3)
    plt.arrow(0,1.,0,-0.5, lw = mlw, head_width = hw, head_length=hl, overhang=oh, ec='none', fill=True, fc='Black', ls='--', zorder=3)
    # plt.arrow(-1.,1.,0,-0.5, lw = mlw, head_width = hw, head_length=hl, overhang=oh, ec='none', fill=True, fc='Black', ls='--', zorder=3)

    # plt.arrow(1 , -1 , 0, 0.5, lw = mlw, head_width = hw, head_length=hl, overhang=oh, ec='none', fill=True, fc='Black', ls='--', zorder=3)
    plt.arrow(0 , -1 , 0, 0.5, lw = mlw, head_width = hw, head_length=hl, overhang=oh, ec='none', fill=True, fc='Black', ls='--', zorder=3)
    # plt.arrow(-1 , -1 , 0, 0.5, lw = mlw, head_width = hw, head_length=hl, overhang=oh, ec='none', fill=True, fc='Black', ls='--', zorder=3)
    fprad = 0.03
    fp1 = plt.Circle((0, 0), fprad, color='Black', zorder=3)
    fp2 = plt.Circle((-1, 0), fprad, color='Black', zorder=3)
    fp3 = plt.Circle((1, 0), fprad, color='Black',zorder=3)
    ax.add_artist(fp1)
    # ax.add_artist(fp2)
    # ax.add_artist(fp3)


    plt.savefig('streamplot_Saddle.pdf')#,transparent=True)

    plt.show()

plot_streamlines(nx= 20, ny = 20,xlower=-1.2,xupper=1.2,ylower=-1.2,yupper=1.2, velocity_example=ex4)
