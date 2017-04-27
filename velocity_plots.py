'''
Code for plotting the double gyre velocities in a subplot arrangement
'''
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from plotting_code import velocity_2x2_subplot
from scipy.io import netcdf
import interpolating_functions as interp
import time
# Code to access netcdf data file
fh = netcdf.netcdf_file('data/phot_flow_welsch.nc', 'r', mmap=False)

print fh.variables.keys()

X = fh.variables['X'][:]
Y = fh.variables['Y'][:]
U = fh.variables['U'][:]  #T x Y x X sized array (assume its Y x X order)
V = fh.variables['V'][:]
TIME = fh.variables['TIME'][:]


fh.close()

lenunits = 'km' #units of distance for axis labels
INIT_TIME = 1165933440.0  # Epoch time (rel to Jan 1 1970) in seconds of Dec 12 2006 (initial time of magnetograms)
print time.strftime('%d-%b-%Y %H:%M GMT', time.gmtime(INIT_TIME))
REF_TIME = INIT_TIME - TIME[0]
nx = 40
ny = 40
nx_dense = 260
ny_dense = 260

t_hrs = [5.5,6,6.5,7]
t_hrs = np.array(t_hrs)
t = t_hrs*60*60 + TIME[0]

grid_lim_step = 0
X_min,X_max, Y_min, Y_max = (X[grid_lim_step],X[-grid_lim_step-1],Y[grid_lim_step],Y[-grid_lim_step-1])  # Limit initial grid size
# Compute nx*ny grid of coordinates
xx = np.linspace(X_min,X_max, nx)
yy = np.linspace(Y_min,Y_max, ny)
coord_grid = np.array(np.meshgrid(xx,yy,indexing='xy'))
xx_de = np.linspace(X_min,X_max, nx_dense)
yy_de = np.linspace(Y_min,Y_max, ny_dense)
coord_grid_de = np.array(np.meshgrid(xx_de,yy_de,indexing='xy'))

subplot_labels = []
vel =[]
vel_high_res = []

for k in xrange(4):
    # calculate velocities
    vel.append(interp.regular_grid_interpolator_fn(U, V, X, Y, TIME)[0](coord_grid, t[k]))
    vel_high_res.append(interp.regular_grid_interpolator_fn(U, V, X, Y, TIME)[0](coord_grid_de,t[k]))
    subplot_labels.append(time.strftime('%d-%b-%Y %H:%M UT', time.gmtime((t[k] + REF_TIME))))

# print np.shape(vel)
velocity_2x2_subplot(coords=coord_grid,labels=subplot_labels, velocity=vel,
    velocity_high_res=vel_high_res, xlower=X_min, xupper=X_max, ylower=Y_min, yupper=Y_max, save_name="phot_flow_velocity_subplots_chelixclr.pdf", lenunits = lenunits,g=0.7,s=1.2,r=2,sat =1)






# import double_gyre_functions as dg
# # Double gyre velocity parameters
# eps = [0.2,0.2,0.2,0.2]
# amp = [0.1,0.1,0.1,0.1]
# period = np.array([10,10,10,10])
# om = 2.*np.pi/period
# nx = 20
# ny = 10
# t_0 = [2,4,6,8]
# dg_vel = []
# dg_vel_high_res = []
# coords = dg.generate_grid(nx,ny)
# for k in xrange(4):
#     #Calculate double gyre velocities
#     dg_vel.append( dg.analytic_velocity_noglobal(eps[k],amp[k],om[k],
#         coords, time_array = (t_0[k] + np.zeros_like(dg.generate_grid(nx,ny))[0] ))
#         )
#     dg_vel_high_res.append( dg.analytic_velocity_noglobal(eps[k],amp[k],om[k],
#         dg.generate_grid(200,100), time_array = (t_0[k] + np.zeros_like(dg.generate_grid(200,100))[0] ))
#         )
#
# fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')#
#
# cmap = plt.cm.plasma
#
# vel_sq = np.sqrt(np.array(dg_vel_high_res)[:,0]**2+np.array(dg_vel_high_res)[:,1]**2)
# vel_sq_rel = vel_sq/np.max(vel_sq)
# # ax5 = fig.add_axes([0.85, 0.15, 0.03, 0.7])
#
# # ax5.cla()
#
#
# fig.text(0.47, 0.04, 'x', ha='center', va='center', fontsize=10)  #A different way to add axes labels onto plots
# fig.text(0.06, 0.5, 'y', ha='center', va='center', rotation='vertical',fontsize=10)
#
#
# # print np.shape(dg_vel)
# # print np.shape(coords)
#
# scl = 6
# ax1.quiver(coords[0],coords[1], np.array(dg_vel[0][0]), np.array(dg_vel[0][1]), scale = scl, scale_units = 'inches')
# ax2.quiver(coords[0],coords[1], np.array(dg_vel[1][0]), np.array(dg_vel[1][1]), scale = scl, scale_units = 'inches')
# ax3.quiver(coords[0],coords[1], np.array(dg_vel[2][0]), np.array(dg_vel[2][1]), scale = scl, scale_units = 'inches')
# ax4.quiver(coords[0],coords[1], np.array(dg_vel[3][0]), np.array(dg_vel[3][1]), scale = scl, scale_units = 'inches')
#
#
# ax1.set_ylim([0,1])
# ax1.set_xlim([0,2])
# # ax3.set_xticklabels([0.0,0.5,1.0,1.5,2.0])
# ax4.set_ylim([0,1])
# ax4.set_xlim([0,2])
#
#
#
# alpha_im = 0.5
#
# im = ax1.imshow(vel_sq_rel[0], cmap = cmap, origin = 'lower', extent=(0,2,0,1), vmin=0,vmax=1,aspect='auto', alpha=alpha_im)
# ax2.imshow(vel_sq_rel[1], cmap = cmap, origin = 'lower', extent=(0,2,0,1), vmin=0,vmax=1,aspect='auto', alpha=alpha_im)
# ax3.imshow(vel_sq_rel[2], cmap = cmap, origin = 'lower', extent=(0,2,0,1), vmin=0,vmax=1,aspect='auto', alpha=alpha_im)
# ax4.imshow(vel_sq_rel[3], cmap = cmap, origin = 'lower', extent=(0,2,0,1), vmin=0,vmax=1,aspect='auto', alpha=alpha_im)
#
# ax1.axes.autoscale(False)
# ax2.axes.autoscale(False)
# ax3.axes.autoscale(False)
# ax4.axes.autoscale(False)
# ax1.set_title('t = %d ' %t_0[0],fontsize=10)
# ax2.set_title('t = %d ' %t_0[1],fontsize=10)
# ax3.set_title('t = %d ' %t_0[2],fontsize=10)
# ax4.set_title('t = %d ' %t_0[3],fontsize=10)
# # ax1.set_title('$\mathregular{A}$ = %.1f ' %amp[0],fontsize=10)
# # ax2.set_title('$\mathregular{A}$ = %.1f ' %amp[1],fontsize=10)
# # ax3.set_title('$\mathregular{A}$ = %.1f ' %amp[2],fontsize=10)
# # ax4.set_title('$\mathregular{A}$ = %.1f ' %amp[3],fontsize=10)
# # ax1.set_title('$\mathregular{\epsilon}$ = %.1f ' %eps[0],fontsize=10)
# # ax2.set_title('$\mathregular{\epsilon}$ = %.1f ' %eps[1],fontsize=10)
# # ax3.set_title('$\mathregular{\epsilon}$ = %.1f ' %eps[2],fontsize=10)
# # ax4.set_title('$\mathregular{\epsilon}$ = %.1f ' %eps[3],fontsize=10)
#
# fig.subplots_adjust(right=0.8)
# cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
# cbar_ax.set_title('$|\mathbf{{v}}_{\mathregular{{rel}}}|$',fontsize=10,y=1.02,x=0.9)
# cbar = fig.colorbar(im,cbar_ax)
#
# cbar.solids.set_rasterized(True)
# # # cbar.solids.set_edgecolor("face")
# # cbar.set_alpha(1)
# # cbar.draw_all()
# # cbar.set_alpha(0.7)
# # cbar.draw_all()
#
# #
# # ax1.text(-10,9.8,'t = %d $\mu$s' %(t[T]),color='DimGray', fontsize=12)
# # plt.savefig('dg_velocity_plot_t_var_A0-3_eps0-2.pdf', bbox_inches='tight', dpi=1000)
# plt.show()
