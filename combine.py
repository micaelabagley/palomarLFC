#! /usr/bin/env python
#######################################################################
## Micaela Bagley
## June 2014
##
## usage: combine.py [-h]
##                   [--rejection {none,minmax,ccdclip,crreject,sigclip,
##                                 avsigclip,pclip}]
##                   [--combine {average,median,sum}]
##                   wispfield
## 
## Stack images from each filter using IRAF's imcombine task
## 
## positional arguments:
##   wispfield             WISP field for which to stack images. Must 
##                         match name of directory containing all 
##                         relevant FITS files
##
## optional arguments:
##   -h, --help            show this help message and exit
##   --rejection {none,minmax,ccdclip,crreject,sigclip,avsigclip,pclip}
##                         Type of rejection algorithm to use. 
##                         Default is crreject
##   --combine {average,median,sum}
##                         Type of combining operation to perform. 
##                         Default is median.
#######################################################################
import argparse
import numpy as np
import pyfits
import subprocess
import os
from glob import glob
from astropy.table import Table
from scipy.optimize import curve_fit
from scipy.stats import norm
from pyraf import iraf
from iraf import imcombine,imalign

def get_filter(image):
    hdr = pyfits.getheader(image)
    return hdr['FILTER'].rstrip("'")


def write_file(value_list, filename):
    '''Write files to be used with IMCOMBINE and IMALIGN'''
    f = open(filename, 'w')
    for x in value_list:
        f.write('%s\n' % x) 
    f.close()


def run_SE(image, wispfield):
    '''Run SExtractor on the images'''
    # names of output files
    base = os.path.splitext(image)[0]
    cat = base + '_cat.fits'
    seg = base + '_seg.fits'

    # pixscale
    pixscale = pyfits.getheader(image)['SECPIX1']

    # run SE
    cmd = ('sex ' +image+ ' -c config.sex -CATALOG_NAME ' +cat+\
           ' -THRESH_TYPE RELATIVE -DETECT_MINAREA 5 -DETECT_THRESH ' +\
           '2.2 -ANALYSIS_THRESH 2.2 -WEIGHT_TYPE NONE ' +\
           '-BACK_SIZE 64 -PIXEL_SCALE ' +str(pixscale) +\
           ' -CHECKIMAGE_TYPE SEGMENTATION ' +\
           '-CHECKIMAGE_NAME ' +seg)
    subprocess.call(cmd, shell=True)


def get_rdnoise(wispfield):
    '''Roughly calculate the read noise of the images for use
    with IMCOMBINE's sigma clipping algorithm'''
    nbins = 20
    biaslist = glob(os.path.join(wispfield,'Zero*.fits'))
    rd = np.zeros((len(biaslist)), dtype=float)
    for i,bias in enumerate(biaslist):
        im = pyfits.getdata(bias)
        im = im.flatten()
        # inefficient way to sigma clip the bias to remove crazy outliers
        im = im[np.where((im < np.median(im)+3*np.std(im)) & 
            (im > np.median(im)-3*np.std(im)))]
        im = im[np.where((im < np.median(im)+3*np.std(im)) & 
            (im > np.median(im)-3*np.std(im)))]
        im = im[np.where((im < np.median(im)+3*np.std(im)) & 
            (im > np.median(im)-3*np.std(im)))]

        # make a histogram of the data
        hist, bins = np.histogram(im, bins=nbins)

        width = 0.8 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2

        # best fit of data
        gaussfunc = lambda x,a,b,c,d: a * np.exp(-(x-b)**2/(2*c**2)) + d
        p0 = [np.max(hist), bins[len(bins)/2], 1., 0]
        popt,pcov = curve_fit(gaussfunc, center, hist, p0=p0)
        sigma = popt[2]
        fwhm = sigma * 2.355
        gain = 1.1
        rd[i] = fwhm * gain

        """
        plt.bar(center, hist, width=width)
        plt.plot(center, gaussfunc(center,*popt), 'r')
        plt.xlabel('flux')
        plt.ylabel('number of pixels')

        plt.show()
        plt.close()
        """
    rdnoise = np.median(rd)
    return rdnoise


def run_imcombine(wispfield, rejection, combine):
    '''Set up parameters for IRAF's task IMCOMBINE and run 
       via PYRAF
    '''
    # determine which filters were used this night
    all_images = glob(os.path.join(wispfield, '*_solved.fits'))
    filters = []
    for im in all_images:
        filters.append(get_filter(im))
    unique_filts = np.unique(filters)

    # airmass coefficients
    # http://www.ifa.hawaii.edu/~rgal/science/dposs/ccdproc/ccdproc_node3.html
    gcoeff = -0.152
    rcoeff = -0.094
    icoeff = -0.07

    for f in unique_filts:
        # airmass coefficient for filter
        if f == 'g':
            coeff = gcoeff
        elif f == 'r':
            coeff = rcoeff
        elif f == 'i':
            coeff = icoeff
    
        # get a list of all images with this filter
        filtlist = [x for x in all_images if get_filter(x) == f]
        
        # determine which images to throw out due to bad seeing
        # or other problems
        print 'Check all %s band images for bad seeing or other issues.' % f
        print '\n%s images:' % f
        for im in filtlist:
            print '   %s' % im
        print (
            "\nUse IRAF's IMEXAM to check the FWHM of "
            "a few point sources across each image.\n"
            "If the FWHM's are drastically different "
            "for one of the images, remove it from \n"
            "the directory. Make sure to check the "
            "observing logs as well and remove any \n"
            "images that have problems indicated there.\n ")
        raw_input('Press ENTER to continue...')

        new_list = [x for x in glob(os.path.join(wispfield, '*_solved.fits')) \
                    if x in filtlist]
        nim = len(new_list)

        # file names
        image_list = os.path.join(wispfield, '%s_%s.lst'%(wispfield,f))
        sky_list = os.path.join(wispfield, '%s_%s.sky'%(wispfield,f))
        scale_list = os.path.join(wispfield, '%s_%s.scale'%(wispfield,f))
        output = os.path.join(wispfield, '%s_%s.fits'%(wispfield,f))
        sigmas = os.path.join(wispfield, '%s_%s_sig.fits'%(wispfield,f))
        logfile = os.path.join(wispfield, '%s_imcombine.log'%wispfield)
        ''' 
        # These output masks are produced as IRAF's pixel lists
        bpmask = os.path.join(wispfield, '%s_%s_bp.fits'%(wispfield,f))
        rejmask = os.path.join(wispfield, '%s_%s_rej.fits'%(wispfield,f))
        #'''

        # get the gain and read noise to use for all images
        gain = pyfits.getheader(new_list[0])['GAIN']
        rdnoise = get_rdnoise(wispfield)

        # get median sky values and airmasses 
        median_array = np.zeros((nim), dtype=float)
        airmass = np.zeros((nim), dtype=float)
        rd = np.zeros((nim), dtype=float)
        for i,image in enumerate(new_list):
            # get median 'sky' value
            median_array[i] = np.median(pyfits.getdata(image))
            # get airmass
            airmass[i] = pyfits.getheader(image)['AIRMASS']

        # sort the list by sky value
        #   image with lowest sky value should be listed first
        #   all zero 'corrections' are done wrt 1st image
        #   (makes more sense to subtract off a sky value)
        sort = np.argsort(median_array)
        # make a list of the sorted images as an input to IMCOMBINE
        write_file(np.array(new_list)[sort], image_list)
        # make a list of sorted sky values
        write_file(median_array[sort], sky_list)
        # make a list of sorted airmass corrections
        scale_factor = 10**(0.4 * coeff * airmass[sort])        
        write_file(scale_factor, scale_list)

        # parameters for IMCOMBINE
        params = {'input':'@%s'%image_list, 'output':'@%s'%output,
                  'sigmas':'@%s'%sigmas, 'combine':combine,
                  'reject':rejection, 'scale':'@%s'%scale_list, 
                  'zero':'median', 'mclip':'yes', 'rdnoise':rdnoise,
                  'gain':gain, 'logfile':logfile, 'offsets':'wcs'}
        # choice of how many pixels to keep and/or reject depends on
        # choice of rejection algorithm
        # for minmax rejection
        if rejection == 'minmax':
            # how many images are there? How many high/low pixels to reject?
            print ' '
            print 'There are %i images in the %s band' % (nim, f)
            nrej = raw_input('How many high/low pixels should be rejected? ')
            params['nlow'] = int(nrej)
            params['nhigh'] = int(nrej)
            params['nkeep'] = nim-2*int(nrej)            
        elif rejection != 'none':
            # (the clipping algorithms)
            # how many images are there? Minimum # of pixels to keep
            # (given as positive: min number to keep
            #  given as negative: max number to reject)
            print ' '
            print 'There are %i images in the %s band' % (nim, f)
            nkeep = raw_input('What is the minimum number of pixels to keep? ')
            params['nkeep'] = int(nkeep)

        # run imcombine
        imcombine(**params)        


def run_imalign(images):
    '''Set up parameters for IRAF's task IMALIGN and run via PYRAF'''
    # run SE on images to get x coordinates
    for image in images:
        run_SE(image)

    # read in Palomar catalogs
    if len(Palcats) == 2:
        # 2 filters used
        # g in dual image mode
        pal1 = read_cat(Palcats[0])
        pal1RA = pal1.field('X_WORLD')
        pal1Dec = pal1.field('Y_WORLD')
        pal1x = pal1.field('X_IMGE')
        pal1y = pal1.field('Y_IMAGE')
        # i in single image mode
        pal2 = read_cat(Palcats[1])
        pal2RA = pal2.field('X_WORLD')
        pal2Dec = pal2.field('Y_WORLD')
        pal2x = pal2.field('X_IMAGE')
        pal2y = pal2.field('Y_IMAGE')
    else:
        print 'Not sure what to do with %i Palomar SE catalogs'%len(Palcats)
        exit()

    # match catalogs
    idx,separc = match_cats(pal1RA, pal1Dec, pal2RA, pal2Dec)
    match = (separc.value*3600. <= 0.5)
    print '\n%i i band objs matched to g band objs\n'%idx[match].shape[0]

    # calculate shifts
    pal1x = pal1x[match]
    pal1y = pal1x[match]
    pal2x = pal2x[idx[match]]
    pal2y = pal2y[idx[match]]
    shiftx = pal1x - pal2x  
    shifty = pal1y - pal2y
    # write to file
    shift = ['%.2f %.2f' % (shiftx, shifty)]
    write_file(shift, 'imalign.shifts')
    
    # write x coords to file for reference image
    coords = ['%.3f %.3f'%(xx,yy) for xx,yy in zip(pal1x, pal1y)]
    write_file(coords, 'imalign.coords')

    params = {'input':'%s,%s'%(images[0],images[1]), 'reference':images[0],
              'coords':'imalign.coords', 'shifts':'imalign.shifts',
              'output':'%s,%s'%(images[0],images[1]), 'negative':'no',
              'shiftimages':'yes', 'interp_type':'spline3',
              'boundary_type':'constant', 'constant':0, 'trimimages':'yes'}
    imalign(**params)


def main():
    parser = argparse.ArgumentParser(description=
        "Stack images from each filter using IRAF's imcombine task")
    parser.add_argument('wispfield', type=str, nargs=1, 
        help='WISP field for which to stack images. Must match ' +\
             'name of directory containing all relevant FITS files')
    parser.add_argument('--rejection', type=str, default='crreject',
        choices=['none', 'minmax', 'ccdclip', 'crreject', 'sigclip', 
                 'avsigclip', 'pclip'],
        help='Type of rejection algorithm to use. Default is crreject')
    parser.add_argument('--combine', type=str, default='median',
        choices=['average','median','sum'],
        help='Type of combining operation to perform. Default is median.')
    args = parser.parse_args()
    # if user added a trailing slash to wispfield, 
    # remove for creating file names
    wispfield = args.wispfield[0].strip('/')
    rejection = args.rejection
    combine = args.combine

    # get all scripts and files set up for running imcombine in IRAF
    run_imcombine(wispfield, rejection, combine)

    # get list of combined images for this field
    images = [x for x in [os.path.join(wispfield,'%s_g.fits'%wispfield), \
               os.path.join(wispfield,'%s_i.fits'%wispfield)] \
               if x in glob(os.path.join(wispfield,'%s_*.fits'%wispfield))]
    images.sort()
    if len(images) > 1:
        # align images with IMALIGN
        run_imalign(images)

if __name__ == '__main__':
    main()

