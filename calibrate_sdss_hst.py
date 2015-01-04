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
import numpy as np
from glob import glob
import pyfits
import subprocess,os,shutil
from astropy.coordinates import SkyCoord,match_coordinates_sky
import astropy.units as u
from astropy.table import Table
from astropy.io import ascii
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
from utils.match_cats import match_cats
from scipy.optimize import curve_fit
import sdss_histograms
from align_images import run_wregister


def read_cat(catfile):
    '''Read in fits tables'''
    f = pyfits.open(catfile)
    cat = f[1].data
    f.close()
    return cat


def filter_combo(Palim, UVISims):
    '''Return the filter combination that should be used for 
       determination of the color term.
    '''
    # get Palomar SDSS filters
    hdr = pyfits.getheader(Palim)
    Pal_filt = hdr['FILTER'].rstrip("'")
    
    # get UVIS filters
    UVIS_filt = []
    for im in UVISims:
        hdr = pyfits.getheader(im)
        UVIS_filt.append(hdr['FILTER'])

    # are both UVIS filters available?
    if len(UVIS_filt) > 1:
        print '\nNOTE: Both UVIS filters -- %s -- are available for '+\
            'this field\n' % [x for x in UVIS_filt]
   
    # possible filter combinations for photometric calibration
    if Pal_filt == 'g':
        use_filt = [x for x in UVIS_filt if x == 'F600LP' or x == 'F814W']
        use_im = [x for x in UVISims if \
                  os.path.basename(x).split('_')[0] == 'F600LP'
    if Pal_filt == 'i':
        use_filt = [x for x in UVIS_filt if x == 'F475X' or x == 'F606W']
        use_im = [x for x in UVISims if \
                  os.path.basename(x).split('_')[0] == 'F475X'
    
    use_filt.append(Pal_filt)
    return use_filt, use_im
 
   
def estimate_exptime(base):
    '''Estimate the effective exposure time of a combined image based on
       how many images were combined and roughly accounting for 
       rejected pixels.
    '''
    # add up exptimes from all images used to create combined image
    data = np.genfromtxt(base+'.lst', dtype=[('l','S50')])
    ims = data['l']
    exptime = 0.
    for i,v in enumerate(ims):
        expt = pyfits.getheader(ims[i])['EXPTIME']
        exptime += expt
    # remove 1 or 2 images' worth of time to account for rejected pixels
    no_ims = len(ims)
    avg_time = exptime / no_ims
    if no_ims > 5:
        exptime = exptime - 2.*avg_time
    if no_ims <= 5:
        exptime = exptime - 1.*avg_time
    return exptime


def make_UVIS_rms(UVISim):
    '''Make UVIS weight map and rms map for use with SExtractor'''
    UVISbase = os.path.splitext(UVISim)[0]
    whtim = UVISbase + '_wht.fits'
    rmsim = UVISbase + '_rms.fits'

    HDUList = pyfits.open(UVISim)
    prihdr =  HDUList[0].header
    sci = HDUList['SCI'].data
    scihdr = HDUList['SCI'].header
    wht = HDUList['WHT'].data
    whthdr = HDUList['WHT'].header
    HDUList.close()

    # wht image
    # create new HDUList
    newpri = pyfits.PrimaryHDU(header=prihdr)
    whtdata = pyfits.PrimaryHDU(header=whthdr, data=wht)
    whtlist = pyfits.HDUList(hdus=newpri)
    # append SCI data and header
    whtlist.append(whtdata)
    whtlist.update_extend()
    whtlist.writeto(whtim, clobber=True)

    # rms image
    rms = np.where(wht != 0, 1/np.sqrt(wht), np.nan)
    rmsdata = pyfits.PrimaryHDU(header=whthdr, data=rms)
    rmslist = pyfits.HDUList(hdus=newpri)
    rmslist.append(rmsdata)
    rmslist.update_extend()
    rmslist.writeto(rmsim, clobber=True)
  

def get_zp(UVISim):
    '''Get zero point for WFC3 UVIS filter'''
    # get DATE-OBS from header
    obsdate = pyfits.getheader(UVISim)['DATE-OBS']
    date = datetime.strptime(obsdate, '%Y-%m-%d')

    zp = {}
    if date.date() >= HSTdate.date():
        # new zero points
        z['F475X'] = 26.1579
        z['F600LP'] = 25.8746
        z['F606W'] = 26.0691
        z['F814W'] = 25.0985

    if date.date() < HSTdate.date():
        # old zero points
        z['F475X'] = 26.15
        z['F600LP'] = 25.85
        z['F606W'] = 26.08
        z['F814W'] = 25.09

    return z


def run_SE(UVISim, Palim, wispfield):
    '''Run SExtractor on Palomar and UVIS images. '''
    # UVIS filter
    filt = pyfits.getheader(UVISim)['FILTER']
    # pixscales
    UVISps = pyfits.getheader(UVISim)['D001SCAL'] 
    PALps = pyfits.getheader(Palim)['SECPIX1']
    # file names
    UVISbase = os.path.splitext(UVISim)[0]
    Palbase = os.path.splitext(Palim)[0]
    UVISrms = UVISbase + '_rms.fits'
    Palcat = Palbase + '_calib_cat.fits'
    UVIScat = UVISbase + '_calib_cat.fits'
    Palseg = Palbase + '_calib_seg.fits'
    UVISseg = UVISbase + '_calib_seg.fits'

    # exposure times
    exptime_UVIS = pyfits.getheader(UVISim)['EXPTIME']
    # rough estimate for the exptime of the Palomar combined images
    exptime_Pal = estimate_exptime(Palbase)

    # make UVIS weight map
    make_UVIS_rms(UVISim)
    
    '''MAG ZP'''
    UVISzp = get_zp(UVISim)
    filtzp = UVISzp[filt]
    Palzp = 0.0

    # detect thresh 2.2
    cmd_hst = ('sex %s -c config.sex -CATALOG_NAME %s '% (UVISim, UVIScat) +\
               '-THRESH_TYPE RELATIVE -DETECT_MINAREA 9 -DETECT_THRESH ' +\
               '1.0 -ANALYSIS_THRESH 1.0 -WEIGHT_TYPE MAP_RMS ' +\
               'WEIGHT_IMAGE %s ' % UVISrms +\
               '-PHOT_APERTURES 30 -BACK_SIZE 512 -GAIN %f ' % exptime_UVIS +\
               '-BACK_FILTERSIZE 6,6 -BACKPHOTO_THICK 37.5 ' +\
               '-FILTER_NAME gauss_uvis.conv ' +\
               '-CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME %s '% UVISseg +\
               '-PIXEL_SCALE %f -MAG_ZEROPOINT %f' % (UVISps, filtzp))
    subprocess.call(cmd_i, shell=True)

    cmd_pal = ('sex %s -c config.sex -CATALOG_NAME %s ' % (Palim, Palcat) +\
               '-THRESH_TYPE RELATIVE -DETECT_MINAREA 5 -DETECT_THRESH ' +\
               '1.0 -ANALYSIS_THRESH 1.0 -WEIGHT_TYPE NONE ' +\
               '-PHOT_APERTURES 30 -BACK_SIZE 64 -GAIN %f ' % exptime_Pal +\
               '-CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME %s '%Palseg +\
               '-PIXEL_SCALE %f -MAG_ZEROPOINT %f' % (pixscale, Palzp))
    subprocess.call(cmd, shell=True)


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
        4) Match catalog to WISP UVIS catalog
        5) Calculate color terms from array of zero points. Remove outliers
        6) Calibrate with color term:
                m_cal = m_Pal + alpha[0]*instr_color + alpha[1]
                alpha[0] = slope of line (color term)
                alpha[1] = offset (zero point)
        7) Plot colorterms and residuals 
    '''

    # region file of Palomar objects
    reg = open('Palomar-SDSS-UVIS.reg', 'w')
    reg.write('global color=green dashlist=8 3 width=1 font="helvetica 10 '
            'normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 '
            'move=1 delete=1 include=1 source=1 \n')
    reg.write('fk5 \n')

    # read in Palomar catalog
    if len(Palcats) == 1:
        pal = read_cat(Palcats[0])
        palRA = pal.field('X_WORLD')
        palDec = pal.field('Y_WORLD')
        palMag = pal.field('MAG_AUTO')
        palEMag = pal.field('MAGERR_AUTO')
    else:
        print 'Not sure what to do with %i Palomar SE catalogs'%len(Palcats)
        exit()
    # add all sources to region file
    make_reg(palRA, palDec, reg, 2, 'blue')

    # read in UVIS
    uvis = read_cat(UVIScats[0])
    uvisRA = uvis['X_WORLD']
    uvisDec = uvis['Y_WORLD']
    uvisMag = uvis['MAG_AUTO']
    uvisEMag = uvis['MAGERR_AUTO']
    # add all sources to region file
    make_reg(uvisRA, uvisDec, reg, 1.5, 'green')

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
    idx,separc = match_cats(palRA, palDec, sdssRA, sdssDec)
    match = (separc.value*3600. <= threshold)  
    print '\n%i Palomar objs matched to SDSS objs\n'%idx[match].shape[0]
    # add successfully matched sources to region file with width=4
    make_reg(sdssRA[idx[match]], sdssDec[idx[match]], reg, 1, 'red', width=4)

    # resize photometry arrays
    # we only need sources that are matched to SDSS
    palRA = palRA[match]
    palDec = palDec[match]
    palMag = palMag[match]
    palEMag = palEMag[match]
    sdss_g = sdss_g[idx[match]]
    sdss_i = sdss_i[idx[match]]

    # match Palomar to WISP UVIS
    idx,separc = match_cats(palRA, palDec, uvisRA, uvisDec)
    match = (separc.value*3600. <= threshold)
    print '\n%i Palomar objs matched to UVIS objs\n'%idx[match].shape[0]
    # add successfully matched sources to region file with width=4
    make_reg(uvisRA[idx[match]], uvisDec[idx[match]],reg,1.5,'green',width=4)

    # resize again
    palMag = palMag[match]
    palEMag = palEMag[match]
    sdss_g = sdss_g[match]
    sdss_i = sdss_i[match]
    uvisMag = uvisMag[idx[match]]
    uvisEMag = uvisEMag[idx[match]]

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
#    AUTO_g = AUTO_g + zp_g
#    AUTO_i = AUTO_i + zp_i
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
                  (np.abs(diff_i-zp_i) <= dist_i))
    # resize arrays to drop large outliers
    AUTO_i = AUTO_i[w]
    sdss_i = sdss_i[w]
    AUTO_g = AUTO_g[w]
    sdss_g = sdss_g[w]
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
    ax5.plot(line5[0,:], line5[1,:], 'k')
    ax6.plot(line6[0,:], line6[1,:], 'k')
    ax7.plot(line7[0,:], line7[1,:], 'k')
    ax8.plot(line8[0,:], line8[1,:], 'k')
    ax9.plot(line9[0,:], line9[1,:], 'k')

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
        catfile.write('# Filter  alpha[0]  alpha[1] \n')
        catfile.write('  g   %f   %f\n'%(alpha_g[0],alpha_g[1]))
        catfile.write('  i   %f   %f'%(alpha_i[0],alpha_i[1]))
    


def main():
    parser = argparse.ArgumentParser(description=
        'Calibrate Palomar photometry with the SDSS catalog')
    parser.add_argument('wispfield', type=str, nargs=1,
        help='WISP field for which to calibrate photometry. ' +\
             'Must match name of directory containing all relevant files')
    args = parser.parse_args()
    wispfield = args.wispfield[0].strip('/')

    # get list of images for this field
    Palim = [x for x in [os.path.join(wispfield,'%s_g.fits'%wispfield), \
             os.path.join(wispfield,'%s_i.fits'%wispfield)] \
             if x in glob(os.path.join(wispfield,'%s_*.fits'%wispfield))]
    Palim.sort()

    print '\nCalibrating:   %s\n' % Palim

    # if both SDSS filters are present, quit and use calibrate_sdss.py
    if (os.path.join(wispfield,'%s_g.fits'%wispfield) in Palim) & \
       (os.path.join(wispfield,'%s_i.fits'%wispfield) in Palim):
        print 'Both SDSS filters are present. Use calibrate_sdss.py instead'
        exit()

    # copy the corresponding UVIS image to wispfield directory
    parnum = wispfield.split('WISP')[1]
    WISPdir = '/data/highzgal/PUBLICACCESS/WISPS/data/Par'+parnum+ \
              '/DATA/UVIS/'
    cpfiles = [os.path.join(WISPdir,'F600LP_drz.fits'),
               os.path.join(WISPdir,'F475X_drz.fits')]
    UVISims = []
    for cpf in cpfiles:
        if os.path.exists(cpf):
            shutil.copy(cpf, wispfield)
            UVISims.append(os.path.join(wispfield, os.path.basename(cpf)))

    # UVIS and Palomar images are too different in size. 
    # Alignment won't work.

    # which filters are available for determination of the color term?
    color_filters,use_im = filter_combo(Palim, UVISims)

    # Run SE on Palomar image in dual image mode with the UVIS image
    run_SE(use_im[0], Palim[0], wispfield)    


    ## START HERE
    # calibrate photometry
    threshold = 0.5  # arcsec for matching
    calibrate(Palomar_catalogs, threshold, wispfield)



if __name__ == '__main__':
    main()


