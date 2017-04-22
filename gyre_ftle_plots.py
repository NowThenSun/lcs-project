from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import double_gyre_functions as dg
from plotting_code import cubehelix_cmap

# Double gyre velocity parameters
eps = [0.2,0.2,0.2,0.2]
amp = [0.1,0.1,0.1,0.1]
period = np.array([10,10,10,10])
om = 2.*np.pi/period
<<<<<<< HEAD
nx = 100
=======
nx = 50
>>>>>>> 5640c3f5761a508dce09d72c14619e823b120695
ny = nx/2
aux_grid_spacing = 0.08*2/nx
t_0 = [0,0,0,0]
int_time = [5,10,15,20]
etol = 10**-3
dt_min = 10**-3
dt_max = 0.2

ftle = []

for k in xrange(4):
    #Calculate double gyre velocities

    ftle.append(dg.main(amplitude=amp[k], epsilon=eps[k], omega=om[k],
    	nx=nx, ny=ny, aux_grid_spacing=aux_grid_spacing,
    	t_0=t_0[k], int_time=int_time[k], adaptive_error_tol=etol,
    	dt_min=np.sign(int_time[k])*dt_min, dt_max=np.sign(int_time[k])*dt_max, method = 'dp45'))



#
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')#
#

cmap = cubehelix_cmap(g=2.0,s=-1.2,r=-0.85,sat=1.0)[1]

fig.text(0.47, 0.04, 'x', ha='center', va='center', fontsize=10)  #A different way to add axes labels onto plots
fig.text(0.06, 0.5, 'y', ha='center', va='center', rotation='vertical',fontsize=10)

ax1.set_ylim([0,1])
ax1.set_xlim([0,2])
# ax3.set_xticklabels([0.0,0.5,1.0,1.5,2.0])
ax4.set_ylim([0,1])
ax4.set_xlim([0,2])
ax3.set_ylim([0,1])
ax3.set_xlim([0,2])
ax2.set_ylim([0,1])
ax2.set_xlim([0,2])

vmin = 0
vmax = 0.5
im = ax1.imshow(ftle[0], cmap = cmap, origin = 'lower', extent=(0,2,0,1), vmin=vmin,vmax=vmax,aspect='auto')
ax2.imshow(ftle[1], cmap = cmap, origin = 'lower', extent=(0,2,0,1), vmin=vmin,vmax=vmax,aspect='auto')
ax3.imshow(ftle[2], cmap = cmap, origin = 'lower', extent=(0,2,0,1), vmin=vmin,vmax=vmax,aspect='auto')
ax4.imshow(ftle[3], cmap = cmap, origin = 'lower', extent=(0,2,0,1), vmin=vmin,vmax=vmax,aspect='auto')

ax1.axes.autoscale(False)
ax2.axes.autoscale(False)
ax3.axes.autoscale(False)
ax4.axes.autoscale(False)
# ax1.set_title('$\mathregular{t_0}$ = %d ' %t_0[0],fontsize=10)
# ax2.set_title('$\mathregular{t_0}$ = %d ' %t_0[1],fontsize=10)
# ax3.set_title('$\mathregular{t_0}$ = %d ' %t_0[2],fontsize=10)
# ax4.set_title('$\mathregular{t_0}$ = %d ' %t_0[3],fontsize=10)
ax1.set_title('$\mathregular{T}$ = %d ' %int_time[0],fontsize=10)
ax2.set_title('$\mathregular{T}$ = %d ' %int_time[1],fontsize=10)
ax3.set_title('$\mathregular{T}$ = %d ' %int_time[2],fontsize=10)
ax4.set_title('$\mathregular{T}$ = %d ' %int_time[3],fontsize=10)
# ax1.set_title('$\mathregular{\epsilon}$ = %.1f ' %eps[0],fontsize=10)
# ax2.set_title('$\mathregular{\epsilon}$ = %.1f ' %eps[1],fontsize=10)
# ax3.set_title('$\mathregular{\epsilon}$ = %.1f ' %eps[2],fontsize=10)
# ax4.set_title('$\mathregular{\epsilon}$ = %.1f ' %eps[3],fontsize=10)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
# if np.sign(int_time[0]) == +1:
#     cbar_ax.set_title('f-FTLE',fontsize=10,y=1.02,x=0.6)
# else:
#     cbar_ax.set_title('b-FTLE',fontsize=10,y=1.02,x=0.6)
cbar_ax.set_title('FTLE',fontsize=10,y=1.02,x=0.6)

cbar = fig.colorbar(im,cbar_ax)

# cbar.solids.set_rasterized(True)
# # cbar.solids.set_edgecolor("face")
# cbar.set_alpha(1)
# cbar.draw_all()
# cbar.set_alpha(0.7)
# cbar.draw_all()

#
# ax1.text(-10,9.8,'t = %d $\mu$s' %(t[T]),color='DimGray', fontsize=12)
plt.savefig('dg_ftle_plot_T_var_A0-1_eps0-2_t00_vmax0-5.pdf', bbox_inches='tight')
plt.show()
