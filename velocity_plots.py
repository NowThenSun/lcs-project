'''
Code for plotting the double gyre velocities in a subplot arrangement
'''
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import double_gyre_functions as dg

# Double gyre velocity parameters
eps = [0.,0.1,0.2,0.5]
amp = [0.1,0.1,0.1,0.1]
period = [10,10,10,10]
om = 2.*np.pi/period
coord_grid = 


for k in xrange(4):
    #Calculate double gyre velocities
    dg.analytic_velocity_noglobal(eps[k],amp[k],om[k], coord_grid, time_grid)
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
fig.text(0.5, 0.04, 'x (cm)', ha='center', va='center', fontsize=13)  #A different way to add axes labels onto plots
fig.text(0.06, 0.5, 'y (cm)', ha='center', va='center', rotation='vertical',fontsize=13)

cmap = plt.cm.coolwarm

#ax2.set_xlim([-12,12])
#ax3.set_xlim([-12,12])
#ax4.set_xlim([-12,12])

ax1.quiver(dg_vel1)


ax1.set_ylim([-12,12])
ax1.set_xlim([-12,12])

ax2.imshow(JzRel[1], cmap = cmap, origin = 'lower', extent=(-12,12,-12,12), vmin=-1,vmax=1,aspect='auto')
ax3.imshow(JzRel[2], cmap = cmap, origin = 'lower', extent=(-12,12,-12,12), vmin=-1,vmax=1,aspect='auto')
ax4.imshow(JzRel[3], cmap = cmap, origin = 'lower', extent=(-12,12,-12,12), vmin=-1,vmax=1,aspect='auto')
ax1.axes.autoscale(False)
ax2.axes.autoscale(False)
ax3.axes.autoscale(False)
ax4.axes.autoscale(False)
ax1.set_title(r'$z$ = %d cm' %z[Z[0]],fontsize=12)
ax2.set_title(r'$z$ = %d cm' %z[Z[1]],fontsize=12)
ax3.set_title(r'$z$ = %d cm' %z[Z[2]],fontsize=12)
ax4.set_title(r'$z$ = %d cm' %z[Z[3]],fontsize=12)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
cbar_ax.set_title(r'$J_{z}$ rel.',fontsize=12,y=1.02,x=1.005)
fig.colorbar(im, cax=cbar_ax)

ax1.text(-10,9.8,r't = %d $\mu$s' %(t[T]),color='DimGray', fontsize=12)

plt.show()
