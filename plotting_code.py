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


def FTLE_plot(ftle, xlower, xupper, ylower, yupper, int_time=0, t_0=0, adap_error_tol=0, colour_range =(0,0), colour_rescale = 0,
    save_name=False,g=1.0,s=-1.2,r=-0.85,sat=1.0, lenunits = False, label1=0, label2=0):
    '''
    Function that plots a colourmap of the FTLE field
    '''
    fig = plt.figure()
    ax = plt.axes()


    fig.subplots_adjust(right=0.97)
    cbar_ax = fig.add_axes([0.9, 0.12, 0.03, 0.7])

    if colour_range == (0,0):
        # Automatic colour bar range
        im = ax.imshow(ftle, interpolation='none', origin='lower', extent=(xlower,xupper, ylower, yupper),
            cmap=cubehelix_cmap(g,s,r,sat)[1])
        cb = fig.colorbar(im, cax=cbar_ax,format='%.1e')

    else:
        if not colour_rescale == 0:
            im = ax.imshow(ftle/colour_rescale, interpolation='none', origin='lower', extent=(xlower,xupper, ylower, yupper),
                cmap=cubehelix_cmap(g,s,r,sat)[1],
                vmin=colour_range[0]/colour_rescale,vmax=colour_range[1]/colour_rescale) #,aspect='auto'
            cb = fig.colorbar(im, cax=cbar_ax,format='%.1e')
        else:
            im = ax.imshow(ftle, interpolation='none', origin='lower', extent=(xlower,xupper, ylower, yupper),
                cmap=cubehelix_cmap(g,s,r,sat)[1],
                vmin=colour_range[0],vmax=colour_range[1]) #,aspect='auto'
            cb = fig.colorbar(im, cax=cbar_ax,format='%.1e')

    if not t_0 == 0:
        ax.text(0.1,1.02,'t_0 = %.1f' %t_0, transform=ax.transAxes)

    if not int_time == 0:
        ax.text(0.7,1.02,'T = %.1f' %int_time, transform=ax.transAxes)
    #ax.text(0.3,1.02,'average dt = %.2e' %np.average(dt), transform=ax.transAxes)

    if not label1 == 0:
        ax.text(0.1,1.02,label1, transform=ax.transAxes)

    if not label2 == 0:
        ax.text(0.7,1.02,label2, transform=ax.transAxes)

    if not adap_error_tol == 0:
        ax.text(0.3,1.02,'error tol in dt= %r' %adap_error_tol, transform=ax.transAxes)


    # cbar_ax = fig.add_axes([0.855, 0.15, 0.025, 0.75])
    #ax1.text(0.8,0.9,r'$t$ = %d $\mu$s' %t[T],fontsize=13,transform=ax1.transAxes, color='Azure')
    if lenunits == False:
        ax.set_xlabel('x')
        ax.set_ylabel('y')
    else:
        ax.set_xlabel('x (%s)' %lenunits)
        ax.set_ylabel('y (%s)' %lenunits)

    if save_name == False:
        plt.show()
    else:
        plt.savefig(save_name, bbox_inches='tight')

def velocity_2x2_subplot(coords,labels, velocity, velocity_high_res, xlower, xupper, ylower, yupper, save_name=False, lenunits = False):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')#

    cmap = plt.cm.plasma
    print np.shape(velocity_high_res)
    vel_sq = np.sqrt(np.array(velocity_high_res)[:,0]**2+np.array(velocity_high_res)[:,1]**2)
    vel_sq_rel = vel_sq/np.max(vel_sq)
    # ax5 = fig.add_axes([0.85, 0.15, 0.03, 0.7])

    # ax5.cla()

    if lenunits ==False:
        fig.text(0.47, 0.04, 'x', ha='center', va='center', fontsize=10)  #A different way to add axes labels onto plots
        fig.text(0.06, 0.5, 'y', ha='center', va='center', rotation='vertical',fontsize=10)
    else:
        fig.text(0.47, 0.04, 'x (%s)' %lenunits, ha='center', va='center', fontsize=10)  #A different way to add axes labels onto plots
        fig.text(0.06, 0.5, 'y (%s)' %lenunits, ha='center', va='center', rotation='vertical',fontsize=10)

    # print np.shape(velocity)
    # print np.shape(coords)

    scl = 6
    ax1.quiver(coords[0],coords[1], np.array(velocity[0][0]), np.array(velocity[0][1]), scale = scl, scale_units = 'inches')
    ax2.quiver(coords[0],coords[1], np.array(velocity[1][0]), np.array(velocity[1][1]), scale = scl, scale_units = 'inches')
    ax3.quiver(coords[0],coords[1], np.array(velocity[2][0]), np.array(velocity[2][1]), scale = scl, scale_units = 'inches')
    ax4.quiver(coords[0],coords[1], np.array(velocity[3][0]), np.array(velocity[3][1]), scale = scl, scale_units = 'inches')


    ax1.set_ylim([ylower,yupper])
    ax1.set_xlim([xlower,xupper])
    # ax3.set_xticklabels([0.0,0.5,1.0,1.5,2.0])
    ax4.set_ylim([ylower,yupper])
    ax4.set_xlim([xlower,xupper])



    alpha_im = 0.5
    vel_max = np.max(vel_sq)
    vel_min = np.min(vel_sq)
    im = ax1.imshow(vel_sq[0], cmap = cmap, origin = 'lower', extent=(xlower,xupper,ylower,yupper), vmin=vel_min,vmax=vel_max,aspect='auto', alpha=alpha_im)
    ax2.imshow(vel_sq[1], cmap = cmap, origin = 'lower', extent=(xlower,xupper,ylower,yupper), vmin=vel_min,vmax=vel_max,aspect='auto', alpha=alpha_im)
    ax3.imshow(vel_sq[2], cmap = cmap, origin = 'lower', extent=(xlower,xupper,ylower,yupper), vmin=vel_min,vmax=vel_max,aspect='auto', alpha=alpha_im)
    ax4.imshow(vel_sq[3], cmap = cmap, origin = 'lower', extent=(xlower,xupper,ylower,yupper), vmin=vel_min,vmax=vel_max,aspect='auto', alpha=alpha_im)
    print vel_min,vel_max
    ax1.axes.autoscale(False)
    ax2.axes.autoscale(False)
    ax3.axes.autoscale(False)
    ax4.axes.autoscale(False)
    ax1.set_title(labels[0], fontsize=10)
    ax2.set_title(labels[1], fontsize=10)
    ax3.set_title(labels[2], fontsize=10)
    ax4.set_title(labels[3],fontsize=10)
    # ax1.set_title('$\mathregular{A}$ = %.1f ' %amp[0],fontsize=10)
    # ax2.set_title('$\mathregular{A}$ = %.1f ' %amp[1],fontsize=10)
    # ax3.set_title('$\mathregular{A}$ = %.1f ' %amp[2],fontsize=10)
    # ax4.set_title('$\mathregular{A}$ = %.1f ' %amp[3],fontsize=10)
    # ax1.set_title('$\mathregular{\epsilon}$ = %.1f ' %eps[0],fontsize=10)
    # ax2.set_title('$\mathregular{\epsilon}$ = %.1f ' %eps[1],fontsize=10)
    # ax3.set_title('$\mathregular{\epsilon}$ = %.1f ' %eps[2],fontsize=10)
    # ax4.set_title('$\mathregular{\epsilon}$ = %.1f ' %eps[3],fontsize=10)

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
    if lenunits == False:
        cbar_ax.set_title('$|\mathbf{v}|$',fontsize=10,y=1.02,x=0.5)
    else:
        cbar_ax.set_title('$|\mathbf{v}|$ (%s/s)' %lenunits,fontsize=9,y=1.02,x=0.9)

    cbar = fig.colorbar(im,cbar_ax)

    cbar.solids.set_rasterized(True)
    # # cbar.solids.set_edgecolor("face")
    # cbar.set_alpha(1)
    # cbar.draw_all()
    # cbar.set_alpha(0.7)
    # cbar.draw_all()

    #
    # ax1.text(-10,9.8,'t = %d $\mu$s' %(t[T]),color='DimGray', fontsize=12)
    # plt.savefig('velocityocity_plot_t_var_A0-3_eps0-2.pdf', bbox_inches='tight', dpi=1000)
    if save_name == False:
        plt.show()
    else:
        plt.savefig(save_name, bbox_inches='tight', dpi=1000)


def velocity_singleplot(coords, velocity, velocity_high_res, xlower, xupper, ylower, yupper, save_name=False, lenunits = False):
    fig, ax1 = plt.subplots(1, 1, sharex='col', sharey='row')#

    cmap = plt.cm.inferno
    vel_sq = np.sqrt(np.array(velocity_high_res)[0]**2+np.array(velocity_high_res)[1]**2)
    vel_sq_rel = vel_sq/np.max(vel_sq)
    print np.shape(velocity_high_res)
    print np.shape(velocity)
    if lenunits ==False:
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
    else:
        fig.text(0.47, 0.04, 'x (%s)' %lenunits, ha='center', va='center', fontsize=10)  #A different way to add axes labels onto plots
        fig.text(0.06, 0.5, 'y (%s)' %lenunits, ha='center', va='center', rotation='vertical',fontsize=10)

    scl = 6
    ax1.quiver(coords[0],coords[1], np.array(velocity[0]), np.array(velocity[1]), scale = scl, scale_units = 'inches')

    ax1.set_ylim([ylower,yupper])
    ax1.set_xlim([xlower,xupper])

    alpha_im = 0.5
    vel_max = np.max(vel_sq)
    vel_min = np.min(vel_sq)
    im = ax1.imshow(vel_sq, cmap = cmap, origin = 'lower', extent=(xlower,xupper,ylower,yupper), vmin=vel_min,vmax=vel_max, alpha=alpha_im) #aspect='auto',
    # ax1.axes.autoscale(False)

    fig.subplots_adjust(right=0.97)
    cbar_ax = fig.add_axes([0.9, 0.12, 0.03, 0.7])
    if lenunits == False:
        cbar_ax.set_title('$|\mathbf{v}|$',fontsize=10,y=1.02,x=0.5)
    else:
        cbar_ax.set_title('$|\mathbf{v}|$ (%s/s)' %lenunits,fontsize=9,y=1.02,x=0.9)

    cbar = fig.colorbar(im,cbar_ax)
    cbar.solids.set_rasterized(True)

    if save_name == False:
        plt.show()
    else:
        plt.savefig(save_name, bbox_inches='tight', dpi=1000)
        plt.show()
