#! /usr/bin/env python
#######################################################################
## Micaela Bagley
## November 2014
##
## usage: calibrate_sdss.py [-h] wispfield
##
## Calibrate Palomar photometry with the SDSS catalog
##
## positional arguments:
##   wispfield   WISP field for which to calibrate photometry. Must 
##               match name of directory containing all relevant files
##
## optional arguments:
##   -h, --help  show this help message and exit
#######################################################################
import argparse
import os
from glob import glob
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
from match_cats import match_cats
import sdss_histograms
from run_SE import run_SE
from maglim import find_maglim


def read_cat(catfile):
    '''Read in fits tables'''
    f = fits.open(catfile)
    cat = f[1].data
    f.close()
    return cat


def calc_zp(palomar, sdss):
    '''Calculate the zero point shift for a Palomar catalog
       between SDSS magnitudes and Palomar instrumental magnitudes.
    '''
    # calculate difference between Palomar and SDSS photometry
    diff = (sdss - palomar)
    zp = np.median(diff)
    print '\nZeropoint:   %f' % zp
    print 'Stddev:        %f' % np.std(diff)
    return zp,diff


def fit_line(xx, yy):
    '''Fit a best-fit line to the data.'''
    # fit a line
    alpha = np.polyfit(xx, yy, 1)
    xvals = np.arange(np.min(xx)-2, np.max(xx)+2, 0.1)
    line = np.zeros((2,xvals.shape[0]), dtype=float)
    line[0,:] = xvals
    line[1,:] = alpha[0]*line[0,:] + alpha[1]
    return alpha,line


def y_distance(xx, yy, line):
    '''Calculate the y-distance of each point off of a line.
       Cut sources that are too far away
    '''
    dist = np.zeros(xx.shape, dtype=float)
    for i,v in enumerate(xx):
        dist[i] = np.abs(line[1,:][np.argmin(np.abs(line[0,:] - v))] - yy[i])
    return dist


def colorterms(diff, color, ax, cutoff):
    '''Calculate color terms by:
       1) Fitting a line to the arrays of zero point shifts (diff*)
       2) Removing outliers that are farther from the line than cutoff
       3) Re-fit a line to the remaining points

       Plot the zero point shifts as a function of color. 
       Indicate points identified as outliers in red.
       Plot both lines for reference. Better line is in black.

       Return the parameters for the better lines for both filters
       and array indices for the good (non-outlier) sources
    '''
    # fit lines 
    a,line = fit_line(color, diff)
    # get distance of each point off of the line
    dist = y_distance(color, diff, line)
    # fit second line with outliers removed
    good = np.where(dist < cutoff) 
    bad = np.where(dist > cutoff)
    a2,line2 = fit_line(color[good], diff[good])

    # plot all points 
    ax.scatter(color, diff, marker='o', c='b', s=20, edgecolor='none')
    # plot the points that are outliers
    ax.scatter(color[bad], diff[bad], marker='o', c='r', s=25,edgecolor='none')
    # plot lines 
    ax.plot(line[0,:], line[1,:], 'r', linewidth=1.5)
    # lines fit to points excluding outliers
    ax.plot(line2[0,:], line2[1,:], 'k', linewidth=1.5)

    return a2,good


def make_reg(RA, Dec, reg, radius, color, width=1):
    '''Add circular regions in fk5 coords to a ds9 region file'''
    for r,d in zip(RA, Dec):
        reg.write('circle(%.6f,%.6f,%f") # color=%s width=%i\n' % \
            (r,d,radius,color,width))


def calibrate(Palcat, filt, threshold, wispfield, cutoff=0.2):
    '''Calibrate Palomar photometry by:
        1) Match catalog to SDSS
        2) Calculate 0th order zeropoint (no color term)
        3) Consider only sources within 1sigma of median zero point
        4) Calculate color terms from array of zero points. Remove outliers
        5) Calibrate with color term:
                m_cal = m_Pal + alpha[0]*instr_color + alpha[1]
                alpha[0] = slope of line (color term)
                alpha[1] = offset (zero point)
        6) Plot colorterms and residuals 
    '''

    # region file of Palomar objects
    reg = open(os.path.join(wispfield,'Palomar-SDSS.reg'), 'w')
    reg.write('global color=green dashlist=8 3 width=1 font="helvetica 10 '
            'normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 '
            'move=1 delete=1 include=1 source=1 \n')
    reg.write('fk5 \n')

    # read in Palomar catalogs
    pal = read_cat(Palcat)
    palRA = pal['X_WORLD']
    palDec = pal['Y_WORLD']
    palMag = pal['MAG_AUTO']
    epalMag = pal['MAGERR_AUTO']
    # add all sources to region file
    make_reg(palRA, palDec, reg, 2, 'blue')

    # read in SDSS 
    sdss = read_cat(os.path.join(wispfield,'result.fits'))
    # take only the SDSS sources with a S/N >= 10 in both bands
    wSN = np.where((1.0875/(sdss['Err_g']) >= 10.) & 
                   (1.0875/(sdss['Err_i']) >= 10.) & (sdss['g'] >= 14))
    sdssRA = sdss['ra'][wSN]
    sdssDec = sdss['dec'][wSN]
    sdss_g = sdss['g'][wSN]
    sdss_i = sdss['i'][wSN]
    # make plot of SDSS photometry and show the S/N >= 10 sources
    # output is called sdss_hist.pdf
    sdss_histograms.plot_limit(wispfield)
    # add all SDSS sources to region file
    make_reg(sdssRA, sdssDec, reg, 1, 'red')

    # match Palomar to SDSS 
    idx,separc = match_cats(palRA, palDec, sdssRA, sdssDec)
    match = (separc.value*3600. <= threshold)  
    print '\n%i Palomar objs matched to SDSS objs\n'%idx[match].shape[0]
    # add successfully matched sources to region file with width=4
    make_reg(sdssRA[idx[match]], sdssDec[idx[match]], reg, 1, 'red', width=4)
    
    # resize photometry arrays
    # we only need sources that are matched to SDSS
    palMag = palMag[match]
    epalMag = epalMag[match]
    sdss_g = sdss_g[idx[match]]
    sdss_i = sdss_i[idx[match]]

    ###################
    ''' CALIBRATION '''
    ###################
    # set up plots
    fig = plt.figure(figsize=(4,5.5))
    gs = gridspec.GridSpec(2,1)
    gs.update(left=0.2, right=0.9, top=0.9, bottom=0.1, hspace=0.7) 
    ax1 = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[1,0])
    fig.suptitle(wispfield, fontsize=20)
    plt.figtext(0.2, 0.45, 'Calibrated Palomar Photometry')

    # calculate zero point shifts
    print '%s band:' % filt
    if filt == 'g':
        zp,diff = calc_zp(palMag, sdss_g)
    if filt == 'i':
        zp,diff = calc_zp(palMag, sdss_i)
    palMag = palMag + zp
    # plot median zero point shift for filter
    ax1.plot([-1,4], [zp,zp], 'k:', linewidth=1.5)

    # calculate percentiles
    # consider only sources within 1sigma of median zero point
    lowsig,upsig = np.percentile(diff,15.9),np.percentile(diff,84.1)
    dist = upsig - lowsig
    w = np.where(np.abs(diff-zp) <= dist) 
    # resize arrays to drop large outliers
    palMag = palMag[w]
    sdss_g = sdss_g[w]
    sdss_i = sdss_i[w]
    diff = diff[w]

    # plot zero points as a function of SDSS colors
    # to check for any obvious color terms we are missing
    sdss_color = sdss_g - sdss_i
    alpha,good = colorterms(diff, sdss_color, ax1, cutoff)

    # resize arrays once more to remove outliers based on their
    # y-distances from the best fit line
    palMag = palMag[good]
    sdss_g = sdss_g[good]
    sdss_i = sdss_i[good]
    sdss_color = sdss_color[good]

    # calibrate 
    cal_inst = palMag #- (zp - alpha[1])

    
    if filt == 'g':
        # plot calibrated Palomar photometry against SDSS
        ax2.scatter(sdss_g, cal_inst, marker='o', edgecolor='none')
        # fit and plot lines of best fit
        alpha2,line2 = fit_line(sdss_g, cal_inst)
    if filt == 'i':
        ax2.scatter(sdss_i, cal_inst, marker='o', edgecolor='none')
        alpha2,line2 = fit_line(sdss_i, cal_inst)

    # plot the best fit line
    ax2.plot(line2[0,:], line2[1,:], 'k')
    # plot the one-to-one line
    ax2.plot([-5,30], [-5,30], 'k:')

    # plot axes labels and limits 
    ax1.set_xlabel(r'$(g-i)_{SDSS}$', fontsize=15)
    ax1.set_ylabel(r'$%s_{SDSS}-%s_{Pal}$'%(filt,filt), fontsize=15)
    ax1.set_xticks([-1,0,1,2,3,4])
    ax1.set_xlim(np.min(sdss_color)-0.5, np.max(sdss_color)+0.5)
#    ax1.set_ylim(np.min(diff[good])-0.5, np.max(diff[good])+0.5)

    ax2.set_xlabel('$%s_{SDSS}$'%filt, fontsize=15)
    ax2.set_ylabel('$%s_{Pal,cal}$'%filt, fontsize=15)
    ax2.set_xlim(np.min(sdss_g)-1,np.max(sdss_g)+1)
    ax2.set_ylim(np.min(sdss_g)-1,np.max(sdss_g)+1)

    for ax in [ax1,ax2]:
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    fig.savefig(os.path.join(wispfield, 'sdss_calibration.pdf'))

    # calculate limiting magnitude
    segmap = os.path.join(wispfield, '%s_%s_calib_seg.fits'%(wispfield,filt))
    maglim,sig = find_maglim(os.path.join(wispfield,'%s_%s.fits'%
                             (wispfield,filt)), segmap, zp, wispfield)

    # print out info for calibration
    print '\nCalibration information for %s'%wispfield
    print '\n%i sources used in fit'%sdss_g.shape[0]
    if filt == 'g':
        print 'stddev(SDSSg - Palg) = %f' % np.std(sdss_g - cal_inst)
    if filt == 'i':
        print 'stddev(SDSSi - Pali) = %f' % np.std(sdss_i - cal_inst)
    print '\nCalibration: '
    print '   m_cal = m_Pal + alpha[0]*(g-1)_Pal + alpha[1],  where:'
    print '      alpha[0] = 0.0, alpha[1] = %f' % zp
    print '\nLimiting Magnitudes: '
    print '   %s: %f'%(filt,maglim)
    # print to file
    with open(os.path.join(wispfield,'sdss_calibration.dat'), 'w') as catfile:
        catfile.write('# Calibration information for %s\n'%wispfield)
        catfile.write('# \n')
        catfile.write('# %i sources used in fit\n'%sdss_g.shape[0])
        if filt == 'g':
            catfile.write('# stddev(SDSSg-Palg) = %f\n' % \
                np.std(sdss_g - cal_inst))
        if filt == 'i':
            catfile.write('# stddev(SDSSi-Pali) = %f\n' % \
                np.std(sdss_i - cal_inst))
        catfile.write('# \n')
        catfile.write('# Calibration: \n')
        catfile.write('#    m_cal = m_Pal + alpha[0]*(g-i)_Pal + alpha[1]\n')
        catfile.write('# \n')
        catfile.write('# Filter   alpha[0]   alpha[1]    maglim    sigma \n')
        catfile.write('  %s  %f  %f    %f    %f\n' % (filt,0.0,zp,maglim,sig))
   

def main():
    parser = argparse.ArgumentParser(description=
        'Calibrate Palomar photometry with the SDSS catalog')
    parser.add_argument('wispfield', type=str, nargs=1,
        help='WISP field for which to calibrate photometry. ' +\
             'Must match name of directory containing all relevant files')
    args = parser.parse_args()
    wispfield = args.wispfield[0].strip('/')

    # get list of images for this field
    images = [x for x in [os.path.join(wispfield,'%s_g.fits'%wispfield), \
               os.path.join(wispfield,'%s_i.fits'%wispfield)] \
               if x in glob(os.path.join(wispfield,'%s_*.fits'%wispfield))]
    images.sort()

    # if both SDSS filters are present, quit and use calibrate_sdss.py
    if (os.path.join(wispfield,'%s_g.fits'%wispfield) in images) & \
       (os.path.join(wispfield,'%s_i.fits'%wispfield) in images):
        print 'Both SDSS filters are present. Use calibrate_sdss.py instead'
        exit()

    print '\nCalibrating:   %s\n' % images

    filt = fits.getheader(images[0])['FILTER'].rstrip("'")

    # Run SE on Palomar image
    run_SE([images[0]], 'Calibration', mode='single')
    
    Palcat = glob(os.path.join(wispfield, '%s_*_calib_cat.fits' % wispfield))

    # calibrate photometry
    threshold = 1  # arcsec for matching
    calibrate(Palcat[0], filt, threshold, wispfield)



if __name__ == '__main__':
    main()


