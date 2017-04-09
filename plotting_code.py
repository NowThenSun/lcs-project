'''
Module that contains various code for plotting
'''
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


#mpl.cm.register_cmap(name='mycubehelix',data=mpl._cm.cubehelix(g,s,r,sat)) ##way to add a colourmap directly (from segmentdata) into mpl.cm or something
#cmap = mpl.colors.LinearSegmentedColormap(name='abab', segmentdata =mpl._cm.cubehelix(g,s,r,sat))  ##way to create colormap from colour dictionary
def cubehelix_cmap(g=1.0, s=0.5, r = -1.5, sat = 1.0):
    '''
    ~~~~~~~~~~
    Inputs:
    g : gamma value (can increase intensity of high valued colors for g>1, increase intensity for low values for g<1)
    s : starting color
    r : number of rotations through B -> G -> R
    sat : saturation value (0 for grayscale)
    ~~~~~~~~~~~
    Outputs:
    cubehelix colourmap
    reverse cubehelix colourmap
    '''
    cdict = mpl._cm.cubehelix(g,s,r,sat)
    def rev_fn(f):
        def reverse_f(x):
            return f(1-x)
        return reverse_f
    b_r = rev_fn(cdict['blue'])
    g_r = rev_fn(cdict['green'])
    r_r = rev_fn(cdict['red'])
    cdict_r = {u'blue' : b_r, u'green' : g_r, u'red' : r_r}
    cmap = mpl.colors.LinearSegmentedColormap(name='ch', segmentdata=cdict)
    cmap_r = mpl.colors.LinearSegmentedColormap(name='ch_r', segmentdata=cdict_r)
    return cmap, cmap_r


def FTLE_plot(ftle, xlower, xupper, ylower, yupper, int_time, t_0, adap_error_tol, colour_range =(0,0), save_name=False,g=1.0,s=-1.2,r=-0.85,sat=1.0):
    '''
    Function that plots a colourmap of the FTLE field
    '''
    fig = plt.figure()
    ax = plt.axes()
    if colour_range == (0,0):
        # Automatic colour bar range
        im = ax.imshow(ftle, interpolation='none', origin='lower', extent=(xlower,xupper, ylower, yupper),
            cmap=cubehelix_cmap(g,s,r,sat)[1])
    else:
        im = ax.imshow(ftle, interpolation='none', origin='lower', extent=(xlower,xupper, ylower, yupper),
            cmap=cubehelix_cmap(g,s,r,sat)[1],
            vmin=colour_range[0],vmax=colour_range[1]) #,aspect='auto'

    ax.text(0.8,1.02,'T = %.1f' %int_time, transform=ax.transAxes)
    ax.text(-0.1,1.02,'t_0 = %.1f' %t_0, transform=ax.transAxes)
    #ax.text(0.3,1.02,'average dt = %.2e' %np.average(dt), transform=ax.transAxes)

    if not adap_error_tol == 0:
        ax.text(0.3,1.02,'error tol in dt= %r' %adap_error_tol, transform=ax.transAxes)
    cbar_ax = fig.add_axes([0.855, 0.15, 0.025, 0.75])
    #cbar_ax.set_title('title',fontsize=11,y=1.02,x=1.005)
    #ax1.text(0.8,0.9,r'$t$ = %d $\mu$s' %t[T],fontsize=13,transform=ax1.transAxes, color='Azure')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    cb = fig.colorbar(im, cax=cbar_ax)
    if save_name == False:
        plt.show()
    else:
        plt.savefig(save_name, bbox_inches='tight')
