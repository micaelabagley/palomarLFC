#! /usr/bin/env python
import os
import subprocess
from astropy.io import fits
import ConfigParser

import palomarLFC

def estimate_effective_gain(image):
    '''Estimate the effective gain of a combined image
       based on how many images were combined
    
       All images are normalized by their exptimes before combining
       so they are all in units of counts/s. The SE parameter GAIN
       should be set to gain * total_exptime (Roy Gal; SE for Dummies)
    '''
    wispfield = os.path.dirname(image)
    Palomar_gain = fits.getheader(image)['GAIN']
    ims = fits.getheader(image)['IMCMB*']
    exptime = 0.
    Nim = len(ims)
    for i in range(Nim):
        if os.path.splitext(ims[i])[1] != '.fits':
            ims[i] = ims[i] + '.fits'
        exptime += fits.getheader(os.path.join(wispfield,ims[i]))['exptime']
    gain = Palomar_gain * Nim
    return gain


def single_SE(images, outstr, params):
    '''Run SE in single image mode and output a segmentation map'''
    for image in images:
        base = os.path.splitext(image)[0]
        cat = '%s_%s_cat.fits' % (base, outstr)
        seg = '%s_%s_seg.fits' % (base, outstr)
        params['-catalog_name'] = cat
        params['-checkimage_type'] = 'segmentation'
        params['-checkimage_name'] = seg
        # find effective gain of the combined image
        effgain = estimate_effective_gain(image)
        params['-gain'] = '%.1f'%effgain
        # set up SE parameters
        args = ['sex', image]
        for key,value in params.iteritems():
            args.append(key)
            args.append(value)
        # run SE
        subprocess.check_call(args)


def dual_SE(image_g, image_i, outstr, params):
    '''Run SE in dual image mode using the i band for detection'''
    base_g = os.path.splitext(image_g)[0]
    cat_g = '%s_%s_cat.fits' % (base_g, outstr)
    params['-catalog_name'] = cat_g
    params['-checkimage_type'] = 'none'
    # find effective gains of the combined images
    effgain_g = estimate_effective_gain(image_g)
    params['-gain'] = '%.1f'%effgain_g
    # set up SE parameters
    args = ['sex', image_i, image_g]
    for key,value in params.iteritems():
        args.append(key)
        args.append(value)
    # run SE
    subprocess.check_call(args)


def run_SE(images, section, mode='single', updates={}):
    '''Run SExtractor on Palomar images using the parameters
       in SE_parameters.cfg 

       If section = 'Calibration', the thresholds are set low
       to detect as much light as possible from the sources.
       If section = 'Catalog', the thresholds are set at 2.2 sigma

       Additional parameters that will be set:
               pixel_scale - from header (SECPIX1)
                      gain - from an estimate of the effective exposure time
                             of the combined image
              catalog_name - [base]_calib_cat.fits if section='Calibration'
                             [base]_final_cat.fits if section='Catalog'
           checkimage_type - segmentation if single image mode
                             none if dual image mode
        
       SE may be run in either single or dual-image mode
    '''
    # directory containing SE files
    dirSE = os.path.join(palomarLFC.__path__[0], 'SE_parameters')

    Config = ConfigParser.ConfigParser()
    Config.read(os.path.join(dirSE,'SE_parameters.cfg'))
    # the sections are 
    #   Calibration - for calibrating the images
    #   Catalog - for creating the final catalog               
    options = Config.options(section)
    params = {}
    for option in options:
        params[option] = Config.get(section, option)

    # add absolute path to necessary SE files
    params['-c'] = os.path.join(dirSE, Config.get(section, '-c'))
    params['-parameters_name'] = os.path.join(dirSE, 'default.param')
    params['-starnnw_name'] = os.path.join(dirSE, 'default.nnw')
    params['-filter_name'] = os.path.join(dirSE, 'default.conv')
    
    # add some parameters
    # pixscale
    pixscale = fits.getheader(images[0])['SECPIX1']
    params['-pixel_scale'] = '%f'%pixscale

    # override any parameters or add new ones?
    params.update(updates)

    if section == 'Calibration':
        outstr = 'calib'
    if section == 'Catalog':
        outstr = 'final'
    if section == 'Astrometry':
        outstr = 'astr'

    if mode == 'single':
        single_SE(images, outstr, params)

    if mode == 'dual':
        images.sort()
        # g band is 1st image, i band is 2nd
        single_SE([images[1]], outstr, params)
        dual_SE(images[0], images[1], outstr, params)

