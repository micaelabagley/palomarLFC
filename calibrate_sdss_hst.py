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
from datetime import datetime

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
    hdr = pyfits.getheader(Palim[0])
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
                  os.path.basename(x).split('_')[0] == 'F600LP']
    if Pal_filt == 'i':
        use_filt = [x for x in UVIS_filt if x == 'F475X' or x == 'F606W']
        use_im = [x for x in UVISims if \
                  os.path.basename(x).split('_')[0] == 'F475X']
    
    use_filt.append(Pal_filt)
    use_filt.sort()
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
    # obs date for new photometry
    photdate = '2012-03-06'
    HSTdate = datetime.strptime(photdate, '%Y-%m-%d')

    # get DATE-OBS from header
    obsdate = pyfits.getheader(UVISim)['DATE-OBS']
    date = datetime.strptime(obsdate, '%Y-%m-%d')

    zp = {}
    if date.date() >= HSTdate.date():
        # new zero points
        zp['F475X'] = 26.1579
        zp['F600LP'] = 25.8746
        zp['F606W'] = 26.0691
        zp['F814W'] = 25.0985

    if date.date() < HSTdate.date():
        # old zero points
        zp['F475X'] = 26.15
        zp['F600LP'] = 25.85
        zp['F606W'] = 26.08
        zp['F814W'] = 25.09

    return zp


def run_SE(UVISim, Palim, wispfield):
    '''Run SExtractor on Palomar and UVIS images. '''
    # UVIS filter
    filt = pyfits.getheader(UVISim)['FILTER']
    # pixscales
    UVISps = pyfits.getheader(UVISim)['D001SCAL'] 
    Palps = pyfits.getheader(Palim)['SECPIX1']
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
               '-WEIGHT_IMAGE %s ' % UVISrms +\
               '-PHOT_APERTURES 30 -BACK_SIZE 512 -GAIN %f ' % exptime_UVIS +\
               '-BACK_FILTERSIZE 6,6 -BACKPHOTO_THICK 37.5 ' +\
               '-FILTER_NAME gauss_uvis.conv ' +\
               '-CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME %s '% UVISseg +\
               '-PIXEL_SCALE %f -MAG_ZEROPOINT %f' % (UVISps, filtzp))
    subprocess.call(cmd_hst, shell=True)

    cmd_pal = ('sex %s -c config.sex -CATALOG_NAME %s ' % (Palim, Palcat) +\
               '-THRESH_TYPE RELATIVE -DETECT_MINAREA 5 -DETECT_THRESH ' +\
               '1.0 -ANALYSIS_THRESH 1.0 -WEIGHT_TYPE NONE ' +\
               '-PHOT_APERTURES 30 -BACK_SIZE 64 -GAIN %f ' % exptime_Pal +\
               '-CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME %s '%Palseg +\
               '-PIXEL_SCALE %f -MAG_ZEROPOINT %f' % (Palps, Palzp))
    subprocess.call(cmd_pal, shell=True)


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
       1) Fitting a line to the array of zero point shifts (diff)
       2) Removing outliers that are farther from the line than cutoff
       3) Re-fit a line to the remaining points

       Plot the zero point shifts as a function of color. 
       Indicate points identified as outliers in red.
       Plot both lines for reference. Better line is in black.

       Return the parameters for the better line
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


def calibrate(Palcat, UVIScat, filters, threshold, wispfield, cutoff=0.4):
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
    reg = open(os.path.join(wispfield,'Palomar-SDSS-UVIS.reg'), 'w')
    reg.write('global color=green dashlist=8 3 width=1 font="helvetica 10 '
            'normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 '
            'move=1 delete=1 include=1 source=1 \n')
    reg.write('fk5 \n')

    # read in Palomar catalog
    pal = read_cat(Palcat)
    palRA = pal.field('X_WORLD')
    palDec = pal.field('Y_WORLD')
    palMag = pal.field('MAG_AUTO')
    palEMag = pal.field('MAGERR_AUTO')
    # add all sources to region file
    make_reg(palRA, palDec, reg, 2, 'blue')

    # read in UVIS
    uvis = read_cat(UVIScat)
    uvisRA = uvis['X_WORLD']
    uvisDec = uvis['Y_WORLD']
    uvisMag = uvis['MAG_AUTO']
    uvisEMag = uvis['MAGERR_AUTO']
    # add all sources to region file
    make_reg(uvisRA, uvisDec, reg, 1.5, 'green')

    # read in SDSS 
    sdss = read_cat(os.path.join(wispfield,'result.fits'))
    # take only the SDSS sources with a S/N >= 10 in both bands
#    wSN = np.where((1.0875/(sdss['Err_g']) >= 10.) & 
#                   (1.0875/(sdss['Err_i']) >= 10.))
    sdssRA = sdss['ra']#[wSN]
    sdssDec = sdss['dec']#[wSN]
    # filters is a sorted list: 
    # 1st element is UVIS filter, 2nd is Palomar filter
    if filters[1] == 'g':
        sdssMag = sdss['g']#[wSN]
    if filters[1] == 'i':
        sdssMag = sdss['i']#[wSN]
    
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
    sdssMag = sdssMag[idx[match]]

    # match Palomar to WISP UVIS
    idx,separc = match_cats(palRA, palDec, uvisRA, uvisDec)
    match = (separc.value*3600. <= threshold)
    print '\n%i Palomar objs matched to UVIS objs\n'%idx[match].shape[0]
    # add successfully matched sources to region file with width=4
    make_reg(uvisRA[idx[match]], uvisDec[idx[match]],reg,1.5,'green',width=4)

    # resize again
    palMag = palMag[match]
    palEMag = palEMag[match]
    sdssMag = sdssMag[match]
    uvisMag = uvisMag[idx[match]]
    uvisEMag = uvisEMag[idx[match]]

    ###################
    ''' CALIBRATION '''
    ###################
    # set up plots
    fig = plt.figure(figsize=(8,11))
    gs1 = gridspec.GridSpec(1,2)
    gs1.update(left=0.12, right=0.94, top=0.9, bottom=0.75,
        hspace=0.4, wspace=0.5)
    ax1 = plt.subplot(gs1[0,0])
    ax2 = plt.subplot(gs1[0,1])
    plt.figtext(0.34, 0.65, 'Calibrated Palomar Photometry', fontsize=15)
    gs2 = gridspec.GridSpec(1,2)
    gs2.update(left=0.12, right=0.94, top=0.6, bottom=0.45,
        hspace=0.4, wspace=0.5)
    ax3 = plt.subplot(gs2[0,0])
    ax4 = plt.subplot(gs2[0,1])

    gs3 = gridspec.GridSpec(1,1)
    gs3.update(left=0.37, right=0.7, top=0.38, bottom=0.25)
    ax5 = plt.subplot(gs3[0])


    # calculate zero point shift
    print '%s band:'%filters[1]
    zp,diff = calc_zp(palMag, sdssMag)
    palMag = palMag + zp
    # plot median zero point shift for filter
    for ax in [ax1,ax2]:
        ax.plot([-1,4], [zp,zp], 'k:', linewidth=1.5)

    # calculate percentiles
    # consider only sources within 1sigma of median zero point
    lowsig,upsig = np.percentile(diff,5.9),np.percentile(diff,84.1)
    dist = upsig - lowsig
    w = np.where((np.abs(diff-zp) <= dist) & (sdssMag-uvisMag > -2))
    # resize arrays to drop large outliers
    palMag = palMag[w]
    sdssMag = sdssMag[w]
    uvisMag = uvisMag[w]
    diff = diff[w]

    # get color terms
    if filters[1] == 'g':
        instr_color = palMag - uvisMag
        sdss_color = sdssMag - uvisMag
        instr_color_pstr = r'$%s_{Pal}-%s$'%(filters[1], filters[0])
        instr_color_str = '%s_Pal-%s'%(filters[1], filters[0])
        sdss_color_str = r'$%s_{SDSS}-%s$'%(filters[1], filters[0])
    if filters[1] == 'i':
        instr_color = uvisMag - palMag
        sdss_color = uvisMag - sdssMag
        instr_color_pstr = r'$%s-%s_{Pal}$'%(filters[0], filters[1])
        instr_color_str = '%s-%s_Pal'%(filters[0], filters[1])
        sdss_color_str = r'$%s-%s_{SDSS}$'%(filters[0], filters[1])
    
    alpha,good = colorterms(diff, instr_color, ax1, cutoff)
    # now get color terms from SDSS colors (for checking)
    sdss_alpha,sdss_good = colorterms(diff, sdss_color, ax2, cutoff)

    # set xaxes so they are the same for all plots as a fn of color
    for ax in [ax1,ax2,ax4,ax5]:
        ax.set_xticks([-1,0,1,2,3,4])
        ax.set_xlim(np.min([np.min(sdss_color)-0.5,np.min(instr_color)-0.1]),
                    np.max([np.max(sdss_color)+0.5,np.max(instr_color)+0.1]))
        if ax == ax4:
            ax.set_yticks([-1,0,1,2,3,4])
            ax.set_ylim(
                np.min([np.min(sdss_color)-0.5,np.min(instr_color)-0.1]),
                np.max([np.max(sdss_color)+0.5,np.max(instr_color)+0.1]))

    # resize arrays once more to remove outliers based on their
    # y-distances from the best fit line
    palMag = palMag[good]
    sdssMag = sdssMag[good]
    instr_color = instr_color[good]
    sdss_color = sdss_color[good]

    # plot calibrated photometry
    fig.suptitle(wispfield, fontsize=20)

    # calibrate with color terms
    #cal_inst = palMag + alpha[0]*instr_color #+ alpha[1]
    cal_inst = palMag + alpha[0]*instr_color - (zp - alpha[1])
    if filters[1] == 'g':
        cal_inst_color = cal_inst - uvisMag[good]
        
    if filters[1] == 'i':
        cal_inst_color = uvisMag[good] - cal_inst
    
    # color-code plots by SDSS colors
    w1 = np.where(sdss_color <= 0.5)
    w2 = np.where((sdss_color > 0.5) & (sdss_color <=1.5))
    w3 = np.where((sdss_color > 1.5) & (sdss_color <=2.5))
    w4 = np.where((sdss_color > 2.5) & (sdss_color <=3.5))
    w5 = np.where(sdss_color > 3.5)

    for wcolor,color in zip([w1,w2,w3,w4,w5],['r','#ff7d40','g','b','k']):
        ax3.scatter(sdssMag[wcolor], cal_inst[wcolor], marker='o', 
                    c=color, edgecolor='none')
        ax4.scatter(sdss_color[wcolor], cal_inst_color[wcolor],
                    marker='o', c=color, edgecolor='none')
        ax5.scatter(sdss_color[wcolor], sdssMag[wcolor]-cal_inst[wcolor], 
                    marker='o', c=color, edgecolor='none')

    # fit and plot lines of best fit
    alpha3,line3 = fit_line(sdssMag, cal_inst)
    alpha4,line4 = fit_line(sdss_color, cal_inst_color)
    alpha5,line5 = fit_line(sdss_color, sdssMag-cal_inst)
    ax3.plot(line3[0,:], line3[1,:], 'k')
    ax4.plot(line4[0,:], line4[1,:], 'k')
    ax5.plot(line5[0,:], line5[1,:], 'k')

    # plot the one-to-one line
    for ax in [ax3,ax4]:
        ax.plot([-5,30], [-5,30], 'k:')
    # plot a line at y=0
    ax5.plot([-1,5], [0,0], 'k:')

    # plot axes labels and limits 
    ax1.set_xlabel(instr_color_pstr, fontsize=15)
    for ax in [ax2, ax4, ax5]:
        ax.set_xlabel(sdss_color_str, fontsize=15)
    ax3.set_xlabel(r'$%s_{SDSS}$'%filters[1], fontsize=15)
 
    for ax in [ax1, ax2]:
        ax.set_ylabel(r'$%s_{SDSS}-%s_{Pal}$'%(filters[1],filters[1]), 
            fontsize=15)
    ax3.set_ylabel(r'$%s_{Pal,cal}$'%filters[1], fontsize=15)
    ax4.set_ylabel(r'(%s)$_{cal}$'%instr_color_pstr, fontsize=15)
    ax5.set_ylabel(r'$%s_{SDSS}-%s_{Pal,cal}$'%(filters[1],filters[1]),
        fontsize=15)

    ax3.set_xlim(np.min(sdssMag)-1,np.max(sdssMag)+1)
    ax3.set_ylim(np.min(sdssMag)-1,np.max(sdssMag)+1)

    for ax in [ax1,ax2,ax3,ax4,ax5]:
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    fig.savefig(os.path.join(wispfield, 'sdss_calibration.pdf')) 
   
    # print out info for calibration
    print '\nCalibration information for %s'%wispfield
    print '\n%i sources used in fit'%sdssMag.shape[0]
    print 'stddev(SDSS%s - Pal%s) = %f' % \
        (filters[1], filters[1], np.std(sdssMag - cal_inst))
    print '\nResiduals (slope of linear fit post-calibration): '
    print '   SDSS%s - Pal%s residuals: %f'%(filters[1],filters[1],alpha5[0])
    print '\nCalibration: '
    print '   m_cal = m_Pal + alpha[0]*%s + alpha[1],  where:'%instr_color_str
    print '   g: alpha[0] = %f, alpha[1] = %f'%(alpha[0],alpha[1])
    # print to file
    with open(os.path.join(wispfield,'sdss_calibration.dat'), 'w') as catfile:
        catfile.write('# Calibration information for %s\n'%wispfield)
        catfile.write('# \n')
        catfile.write('# %i sources used in fit\n'%sdssMag.shape[0])
        catfile.write('# stddev(SDSS%s-Pal%s) = %f\n' % \
            (filters[1], filters[1], np.std(sdssMag - cal_inst)))
        catfile.write('# \n')
        catfile.write('# Residuals (slope of linear fit post-calibration):\n')
        catfile.write('#   SDSS%s - Pal%s: %f\n' % 
            (filters[1],filters[1],alpha5[0]))
        catfile.write('# \n')
        catfile.write('# Calibration: \n')
        catfile.write('#    m_cal = m_Pal + alpha[0]*%s + alpha[1]\n' % 
            instr_color_str)
        catfile.write('# Filter  alpha[0]  alpha[1] \n')
        catfile.write('  g   %f   %f\n'%(alpha[0],alpha[1]))


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

    # which filters are available for determination of the color term?
    use_filts,use_im = filter_combo(Palim, UVISims)

    # Run SE on Palomar image in dual image mode with the UVIS image
#    run_SE(use_im[0], Palim[0], wispfield)    

    ## START HERE
    # find palcat, find UVIScat names
    Palcats = glob(os.path.join(wispfield, '%s_*_calib_cat.fits'%wispfield))
    if len(Palcats) != 1:
        print 'Not sure what to do with %i Palomar SE catalogs'%len(Palcats)
        exit()
    UVIScats = glob(os.path.join(wispfield, '*_drz_calib_cat.fits'))
    if len(UVIScats) != 1:
        print 'Not sure what to do with %i UVIS SE catalogs'%len(UVIScats)
        exit()

    # use_filts is already sorted, so first is UVIS filter, second is Palomar
    # calibrate photometry
    threshold = 0.75  # arcsec for matching
    calibrate(Palcats[0], UVIScats[0], use_filts, threshold, wispfield)



if __name__ == '__main__':
    main()


