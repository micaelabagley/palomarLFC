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
import pyfits
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
    f = pyfits.open(catfile)
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


def colorterms(diff_g, diff_i, color, axg, axi, cutoff):
    '''Calculate color terms by:
       1) Fitting a line to the arrays of zero point shifts (diff_*)
       2) Removing outliers that are farther from the line than cutoff
       3) Re-fit a line to the remaining points

       Plot the zero point shifts as a function of color. 
       Indicate points identified as outliers in red.
       Plot both lines for reference. Better line is in black.

       Return the parameters for the better lines for both filters
       and array indices for the good (non-outlier) sources
    '''
    # fit lines 
    a_g,line_g = fit_line(color, diff_g)
    a_i,line_i = fit_line(color, diff_i)
    # get distance of each point off of the line
    dist_g = y_distance(color, diff_g, line_g)
    dist_i = y_distance(color, diff_i, line_i)
    # fit second line with outliers removed
    good = np.where((dist_g < cutoff) & (dist_i < cutoff))
    bad = np.where((dist_g > cutoff) | (dist_i > cutoff))
    a2_g,line2_g = fit_line(color[good], diff_g[good])
    a2_i,line2_i = fit_line(color[good], diff_i[good])

    # plot all points 
    axg.scatter(color, diff_g, marker='o', c='b', s=20, edgecolor='none')
    axi.scatter(color, diff_i, marker='o', c='b', s=20, edgecolor='none')
    # plot the points that are outliers
    axg.scatter(color[bad], diff_g[bad], marker='o', c='r', 
        s=25, edgecolor='none')
    axi.scatter(color[bad], diff_i[bad], marker='o', c='r', 
        s=25, edgecolor='none')
    # plot lines 
    axg.plot(line_g[0,:], line_g[1,:], 'r', linewidth=1.5)
    axi.plot(line_i[0,:], line_i[1,:], 'r', linewidth=1.5)
    # lines fit to points excluding outliers
    axg.plot(line2_g[0,:], line2_g[1,:], 'k', linewidth=1.5)
    axi.plot(line2_i[0,:], line2_i[1,:], 'k', linewidth=1.5)

    return a2_g,a2_i,good


def make_reg(RA, Dec, reg, radius, color, width=1):
    '''Add circular regions in fk5 coords to a ds9 region file'''
    for r,d in zip(RA, Dec):
        reg.write('circle(%.6f,%.6f,%f") # color=%s width=%i\n' % \
            (r,d,radius,color,width))


def calibrate(Palcats, threshold, wispfield, cutoff=0.2):
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
    if len(Palcats) == 2:
        # 2 filters used
        # g in dual image mode
        pal1 = read_cat(Palcats[0])
        pal1RA = pal1.field('X_WORLD')
        pal1Dec = pal1.field('Y_WORLD')
        AUTO_g = pal1.field('MAG_AUTO')
        eAUTO_g = pal1.field('MAGERR_AUTO')
        # i in single image mode
        pal2 = read_cat(Palcats[1])
        pal2RA = pal2.field('X_WORLD')
        pal2Dec = pal2.field('Y_WORLD')
        AUTO_i = pal2.field('MAG_AUTO')
        eAUTO_i = pal2.field('MAGERR_AUTO')
    else:
        print 'Not sure what to do with %i Palomar SE catalogs'%len(Palcats)
        exit()
    # add all sources to region file
    make_reg(pal1RA, pal1Dec, reg, 2, 'blue')

    # read in SDSS 
    sdss = read_cat(os.path.join(wispfield,'result.fits'))
    # take only the SDSS sources with a S/N >= 10 in both bands
    wSN = np.where((1.0875/(sdss['Err_g']) >= 10.) & 
                   (1.0875/(sdss['Err_i']) >= 10.))
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
    idx,separc = match_cats(pal1RA, pal1Dec, sdssRA, sdssDec)
    match = (separc.value*3600. <= threshold)  
    print '\n%i Palomar objs matched to SDSS objs\n'%idx[match].shape[0]
    # add successfully matched sources to region file with width=4
    make_reg(sdssRA[idx[match]], sdssDec[idx[match]], reg, 1, 'red', width=4)
    reg.close()

    # resize photometry arrays
    # we only need sources that are matched to SDSS
    AUTO_g = AUTO_g[match]
    eAUTO_g = eAUTO_g[match]
    AUTO_i = AUTO_i[match]
    eAUTO_i = eAUTO_i[match]
    sdss_g = sdss_g[idx[match]]
    sdss_i = sdss_i[idx[match]]

    ###################
    ''' CALIBRATION '''
    ###################
    # set up plots
    fig = plt.figure(figsize=(8,11))
    gs1 = gridspec.GridSpec(2,2)
    gs1.update(left=0.15, right=0.9, top=0.93, bottom=0.66, 
        hspace=0.45, wspace=0.6)
    ax1 = plt.subplot(gs1[0,0])
    ax2 = plt.subplot(gs1[0,1])
    ax3 = plt.subplot(gs1[1,0])
    ax4 = plt.subplot(gs1[1,1])
    plt.figtext(0.34, 0.59, 'Calibrated Palomar Photometry', fontsize=15)
    gs2 = gridspec.GridSpec(1,2)
    gs2.update(left=0.12, right=0.95, top=0.57, bottom=0.44, 
        hspace=0.4, wspace=0.5)
    ax5 = plt.subplot(gs2[0,0])
    ax6 = plt.subplot(gs2[0,1])
    gs3 = gridspec.GridSpec(1,1)
    gs3.update(left=0.37, right=0.7, top=0.38, bottom=0.25)
    ax7 = plt.subplot(gs3[0])
    gs4 = gridspec.GridSpec(1,2)
    gs4.update(left=0.12, right=0.95, top=0.19, bottom=0.06,
        hspace=0.4, wspace=0.5)
    ax8 = plt.subplot(gs4[0,0])
    ax9 = plt.subplot(gs4[0,1])

    # calculate zero point shifts
    print 'g band:'
    zp_g,diff_g = calc_zp(AUTO_g, sdss_g)
    print 'i band:'
    zp_i,diff_i = calc_zp(AUTO_i, sdss_i)
    # plot median zero point shift for filter
    for ax in [ax1,ax3]:
        ax.plot([-1,4], [zp_g,zp_g], 'k:', linewidth=1.5)
    for ax in [ax2,ax4]:
        ax.plot([-1,4], [zp_i,zp_i], 'k:', linewidth=1.5)

    # calculate percentiles
    # consider only sources within 1sigma of median zero point
    g_lowsig,g_upsig = np.percentile(diff_g, 15.9),np.percentile(diff_g, 84.1)
    i_lowsig,i_upsig = np.percentile(diff_i, 15.9),np.percentile(diff_i, 84.1)
    dist_g = g_upsig - g_lowsig
    dist_i = i_upsig - i_lowsig
    w = np.where( (np.abs(diff_g-zp_g) <= dist_g) & 
                  (np.abs(diff_i-zp_i) <= dist_i) &
                  (AUTO_i != 99.0))
    # resize arrays to drop large outliers
    AUTO_g = AUTO_g[w]
    sdss_g = sdss_g[w]
    AUTO_i = AUTO_i[w]
    sdss_i = sdss_i[w]
    diff_g = diff_g[w]
    diff_i = diff_i[w]

    # get color terms
    # first get color terms from instrumental Palomar colors
    instr_color = AUTO_g - AUTO_i
    
    alpha_g,alpha_i,good = colorterms(diff_g, diff_i, instr_color, 
                                      ax1, ax2, cutoff)
    # now get color terms from SDSS colors (for checking)
    sdss_color = sdss_g - sdss_i
    sdss_alpha_g,sdss_alpha_i,sdss_good = colorterms(diff_g,diff_i, sdss_color,
                                                     ax3, ax4, cutoff)

    # set xaxes so they are the same for all plots as a fn of color
    for ax in [ax1,ax2,ax3,ax4,ax7,ax8,ax9]:
        ax.set_xticks([-1,0,1,2,3,4])
        ax.set_xlim(np.min([np.min(sdss_color)-0.5,np.min(instr_color)-0.1]),
                    np.max([np.max(sdss_color)+0.5,np.max(instr_color)+0.1]))
        if ax == ax7:
            ax.set_yticks([-1,0,1,2,3,4])
            ax.set_ylim(
                np.min([np.min(sdss_color)-0.5,np.min(instr_color)-0.1]),
                np.max([np.max(sdss_color)+0.5,np.max(instr_color)+0.1]))

    # resize arrays once more to remove outliers based on their
    # y-distances from the best fit line
    AUTO_g = AUTO_g[good]
    AUTO_i = AUTO_i[good]
    sdss_g = sdss_g[good]
    sdss_i = sdss_i[good]
    instr_color = AUTO_g - AUTO_i
    sdss_color = sdss_g - sdss_i

    # plot calibrated photometry
    fig.suptitle(wispfield, fontsize=20)

    # calibrate with color terms
    g_cal_inst = AUTO_g + alpha_g[0]*instr_color + alpha_g[1]
    i_cal_inst = AUTO_i + alpha_i[0]*instr_color + alpha_i[1]
    
    # color-code plots by SDSS colors
    w1 = np.where(sdss_g-sdss_i <= 0.5)
    w2 = np.where((sdss_g-sdss_i > 0.5) & (sdss_g-sdss_i <=1.5))
    w3 = np.where((sdss_g-sdss_i > 1.5) & (sdss_g-sdss_i <=2.5))
    w4 = np.where((sdss_g-sdss_i > 2.5) & (sdss_g-sdss_i <=3.5))
    w5 = np.where(sdss_g-sdss_i > 3.5)

    for wcolor,color in zip([w1,w2,w3,w4,w5],['r','#ff7d40','g','b','k']):
        ax5.scatter(sdss_g[wcolor], g_cal_inst[wcolor], marker='o', 
                    c=color, edgecolor='none')
        ax6.scatter(sdss_i[wcolor], i_cal_inst[wcolor], marker='o', 
                    c=color, edgecolor='none')
        ax7.scatter(sdss_color[wcolor], g_cal_inst[wcolor]-i_cal_inst[wcolor],
                    marker='o', c=color, edgecolor='none')
        ax8.scatter(sdss_color[wcolor], sdss_g[wcolor]-g_cal_inst[wcolor], 
                    marker='o', c=color, edgecolor='none')
        ax9.scatter(sdss_color[wcolor], sdss_i[wcolor]-i_cal_inst[wcolor], 
                    marker='o', c=color, edgecolor='none')

    # fit and plot lines of best fit
    alpha5,line5 = fit_line(sdss_g, g_cal_inst)
    alpha6,line6 = fit_line(sdss_i, i_cal_inst)
    alpha7,line7 = fit_line(sdss_color, g_cal_inst-i_cal_inst)
    alpha8,line8 = fit_line(sdss_color, sdss_g-g_cal_inst)
    alpha9,line9 = fit_line(sdss_color, sdss_i-i_cal_inst)
    for ax,line in zip([ax5,ax6,ax7,ax8,ax9],[line5,line6,line7,line8,line9]):
        ax.plot(line[0,:], line[1,:], 'k')

    # plot the one-to-one line
    for ax in [ax5,ax6,ax7]:
        ax.plot([-5,30], [-5,30], 'k:')
    # plot a line at y=0
    for ax in [ax8,ax9]:
        ax.plot([-1,5], [0,0], 'k:')

    # plot axes labels and limits 
    for ax in [ax1,ax2]:
        ax.set_xlabel(r'$(g-i)_{Pal}$', fontsize=15)
    for ax in [ax3,ax4,ax7,ax8,ax9]:
        ax.set_xlabel(r'$(g-i)_{SDSS}$', fontsize=15)
    for ax in [ax1,ax3]:
        ax.set_ylabel(r'$g_{SDSS} - g_{Pal}$', fontsize=15)
    for ax in [ax2,ax4]:
        ax.set_ylabel(r'$i_{SDSS} - i_{Pal}$', fontsize=15)
    ax5.set_xlabel('$g_{SDSS}$', fontsize=15)
    ax5.set_ylabel('$g_{Pal,cal}$', fontsize=15)
    ax6.set_xlabel('$i_{SDSS}$', fontsize=15)
    ax6.set_ylabel('$i_{Pal,cal}$', fontsize=15)
    ax7.set_ylabel(r'$(g-i)_{Pal,cal}$', fontsize=15)
    ax8.set_ylabel(r'$g_{SDSS} - g_{Pal, cal}$', fontsize=15)
    ax9.set_ylabel(r'$i_{SDSS} - i_{Pal,cal}$', fontsize=15)

    ax5.set_xlim(np.min(sdss_g)-1,np.max(sdss_g)+1)
    ax5.set_ylim(np.min(sdss_g)-1,np.max(sdss_g)+1)
    ax6.set_xlim(np.min(sdss_i)-1,np.max(sdss_i)+1)
    ax6.set_ylim(np.min(sdss_i)-1,np.max(sdss_i)+1)
    ax8.set_ylim(-(np.abs(np.max(sdss_g-g_cal_inst)+0.05)), 
                 np.abs(np.max(sdss_g-g_cal_inst)+0.05))
    ax9.set_ylim(-(np.abs(np.max(sdss_i-i_cal_inst)+0.05)), 
                 np.abs(np.max(sdss_i-i_cal_inst)+0.05))

    for ax in [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9]:
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    fig.savefig(os.path.join(wispfield, 'sdss_calibration.pdf'))
    
    # calculate limiting magnitude for each image
    segmap = os.path.join(wispfield,'%s_i_calib_seg.fits'%wispfield)
    maglim_g,sig_g = find_maglim(os.path.join(wispfield,'%s_g.fits'%wispfield), 
                                 segmap, alpha_g[1], wispfield)
    maglim_i,sig_i = find_maglim(os.path.join(wispfield,'%s_i.fits'%wispfield),
                                 segmap, alpha_i[1], wispfield)

    # print out info for calibration
    print '\nCalibration information for %s'%wispfield
    print '\n%i sources used in fit'%sdss_g.shape[0]
    print 'stddev(SDSSg - Palg) = %f' % np.std(sdss_g - g_cal_inst)
    print 'stddev(SDSSi - Pali) = %f' % np.std(sdss_i - i_cal_inst)
    print '\nResiduals (slopes of linear fits post-calibration): '
    print '   SDSSg - Palg residuals: %f'%alpha8[0]
    print '   SDSSi - Pali residuals: %f'%alpha9[0]
    print '\nCalibration: '
    print '   m_cal = m_Pal + alpha[0]*(g-i)_Pal + alpha[1],  where:'
    print '   g: alpha_g[0] = %f, alpha_g[1] = %f'%(alpha_g[0],alpha_g[1])
    print '   i: alpha_i[0] = %f, alpha_i[1] = %f'%(alpha_i[0],alpha_i[1])
    print '\nLimiting Magnitudes: '
    print '   g: %f'%maglim_g
    print '   i: %f'%maglim_i
    # print to file
    with open(os.path.join(wispfield,'sdss_calibration.dat'), 'w') as catfile:
        catfile.write('# Calibration information for %s\n'%wispfield)
        catfile.write('# \n')
        catfile.write('# %i sources used in fit\n'%sdss_g.shape[0])
        catfile.write('# stddev(SDSSg-Palg) = %f\n'%np.std(sdss_g - g_cal_inst))
        catfile.write('# stddev(SDSSi-Pali) = %f\n'%np.std(sdss_i - i_cal_inst))
        catfile.write('# \n')
        catfile.write('# Residuals (slopes of linear fits post-calibration):\n')
        catfile.write('#   SDSSg - Palg: %f\n'%alpha8[0])
        catfile.write('#   SDSSi - Pali: %f\n'%alpha9[0])
        catfile.write('# \n')
        catfile.write('# Calibration: \n')
        catfile.write('#    m_cal = m_Pal + alpha[0]*(g-i)_Pal + alpha[1]\n')
        catfile.write('# \n')
        catfile.write('# Filter   alpha[0]   alpha[1]    maglim    sigma \n')
        catfile.write('  g   %f   %f   %f    %f\n' % 
                        (alpha_g[0],alpha_g[1],maglim_g,sig_g))
        catfile.write('  i   %f   %f   %f    %f' % 
                        (alpha_i[0],alpha_i[1],maglim_i,sig_i))
    

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

    print 'Calibrating:'
    print images
    print ' '

    # Run SE on all Palomar images in dual image mode 
    run_SE(images, 'Calibration', mode='dual')    

    Palcats = glob(os.path.join(wispfield, '%s_*_calib_cat.fits' % wispfield))
    Palcats.sort()

    # calibrate photometry
    threshold = 1 # arcsec for matching
    calibrate(Palcats, threshold, wispfield)



if __name__ == '__main__':
    main()


