#! /usr/bin/env python
from astropy.io import fits
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator


def plot_hist(ax, mags, bins, color, alpha=0.5):
    '''Plot histograms of SDSS magnitudes for a given WISP field 
       using predetermined magnitude bins (so all magnitudes
       and errors use the same bin)
    '''
    hist,bins = np.histogram(mags, bins=bins)
    bincenters = 0.5 * (bins[1:] + bins[:-1])
    width = 1.01 * (bins[1] - bins[0])
    ax.bar(bincenters, hist, align='center', width=width, linewidth=1,
        fill=False, edgecolor=color, alpha=alpha)
    return hist


def bin_errors(ax, mag, emag, magbins, color, alpha=0.5):
    '''Bin all errors for a given filter in predetermined 
       magnitude bins (so all filters use the same bins)
    '''
    binwidth = magbins[1] - magbins[0]
    # consider only sources with reasonable errors
    mag = mag[emag < 5]
    emag = emag[emag < 5]
    # get the median error and standard deviation 
    # of the errors in each bin
    bin_median = np.zeros(magbins.shape, dtype=float)
    std = np.zeros(magbins.shape, dtype=float)
    for i,b in enumerate(magbins):
        # get all errors in this magnitude bin
        errors = emag[(mag >= (b-binwidth/2.)) & (mag < (b+binwidth/2.))]
        if errors.shape[0] != 0:
            bin_median[i] = np.median(errors)
            std[i] = np.std(errors)
    # plot median error for each bin   
    ax.scatter(magbins[bin_median != 0], bin_median[bin_median != 0], 
        marker='o', s=10, edgecolor='none', color=color, alpha=0.7)
    # fill in the errors within 1 stddev of the median
    ax.fill_between(magbins[bin_median != 0], 
        bin_median[bin_median != 0]-std[bin_median != 0], 
        bin_median[bin_median != 0]+std[bin_median != 0],
        edgecolor='none', alpha=alpha, color=color)
    

def plot_limit(wispfield):
    '''Plot the distribution of all SDSS magnitudes for the given 
       wispfield. Plot also the distribution of sources with S/N >= 10.

       Output is sdss_hist.pdf and is useful for checking the 
       SDSS catalog used for calibrating the Palomar photmetry
    '''
    # set up plot
    fig = plt.figure(figsize=(8,6))
    gs = gridspec.GridSpec(2,2)
    gs.update(left=0.1, right=0.95, top=0.9, bottom=0.1)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[1,0])
    ax4 = fig.add_subplot(gs[1,1])
    
    # use same magnitude bins for all binning
    bins = np.arange(12, 28.5, 0.2)
    ebins = np.arange(0, 2, 0.05)

    # read in SDSS catalog
    cat = fits.open(os.path.join(wispfield, 'result.fits'))
    data = cat[1].data
    cat.close()

    # plot the histogram of all magnitudes in purple
    purple = [0.5, 0., 1., 1.]
    ghist1 = plot_hist(ax1, data['g'], bins, purple, alpha=0.2)
    ihist1 = plot_hist(ax2, data['i'], bins, purple, alpha=0.2)
    # plot the histogram of all magnitudes for which S/N >= 10. in blue
    ghist2 = plot_hist(ax1, data['g'][1.0875/(data['Err_g']) >= 10.],
                       bins, 'b', alpha=0.6)
    ihist2 = plot_hist(ax2, data['i'][1.0875/(data['Err_i']) >= 10.],
                       bins, 'b', alpha=0.6)

    # plot the median error in each bin and +/-1 stddev for each bin
    bin_errors(ax3, data['g'], data['Err_g'], bins, purple, alpha=0.2)
    bin_errors(ax4, data['i'], data['Err_i'], bins, purple, alpha=0.2)
    # and for sources with S/N >= 10
    bin_errors(ax3, data['g'][1.0875/(data['Err_g']) >= 10.],
               data['Err_g'][1.0875/(data['Err_g']) >= 10.], bins, 
               'b', alpha=0.6)
    bin_errors(ax4, data['i'][1.0875/(data['Err_i']) >= 10.],
               data['Err_i'][1.0875/(data['Err_i']) >= 10.], bins, 
               'b', alpha=0.6)

    # plot SDSS's magnitude limit in g
    ax1.plot([22.2,22.2], [0,500], 'k--', linewidth=2)
    ax2.plot([21.3,21.3], [0,500], 'k--', linewidth=2)
    ax3.plot([22.2,22.2], [-0.2,3], 'k--', linewidth=2)
    ax4.plot([21.3,21.3], [-0.2,3], 'k--', linewidth=2)

    ax3.plot([13,14.4], [0.42,0.42], 'k', linewidth=1.5)
    ax3.text(14.9, 0.41, r'Mag limit, 22.2', fontsize=10) 
    ax4.plot([13,14.4], [0.42,0.42], 'k', linewidth=1.5)
    ax4.text(14.9, 0.41, r'Mag limit, 21.3', fontsize=10)
    ax3.text(13, 0.36, r'All Magnitudes', color=purple, fontsize=10)
    ax3.text(13, 0.32, r'Mags (S/N > 10)', color='b', fontsize=10)

    ax3.set_xlabel(r'SDSS $g$')
    ax4.set_xlabel(r'SDSS $i$')

    for ax in [ax1,ax2,ax3,ax4]:
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.set_xlim(12,30)
    ax1.set_ylim(0,1.1*np.max([ghist1, ghist2]))
    ax2.set_ylim(0,1.1*np.max([ihist1, ihist2]))
    ax3.set_ylim(-0.01,0.5)
    ax4.set_ylim(-0.01,0.5)
    ax1.set_title(r'SDSS $g$')
    ax2.set_title(r'SDSS $i$')

    fig.suptitle(wispfield, fontsize=20)
    fig.savefig(os.path.join(wispfield,'sdss_hist.pdf'))


def main():
    pass


if __name__ == '__main__':
    main()
