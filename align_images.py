#! /usr/bin/env python
#######################################################################
## Micaela Bagley
## November 2014
##
## usage: align_images.py [-h] wispfield
##
## Align combined images using IRAF's wregister task
##
## positional arguments:
##   wispfield   WISP field for which to align images. Must match name of
##               directory containing all relevant FITS files
##
## optional arguments:
##   -h, --help  show this help message and exit
#######################################################################
import argparse
import subprocess
import os
from glob import glob
import numpy as np
from astropy.io import fits
from pyraf.iraf import wregister


def run_wregister(image, reference):
    '''Align image to reference image using the WCS in 
       the header and IRAF's WREGISTER task
    '''
    # rename original image so it is not overwritten
    old = '%s_old.fits' % os.path.splitext(image)[0]
    os.rename(image, old)
    
    params = {'input':'%s[0]'%old, 
              'reference':'%s[0]'%reference,
              'output':image, 
              'xmin':'INDEF', 'xmax':'INDEF', 'ymin':'INDEF', 'ymax':'INDEF', 
              'wcs':'world', 'transpose':'no', 'fitgeom':'general',
              'function':'polynomial', 'calctype':'double', 
              'geometry':'geometric', 'interpolant':'spline3',
              'boundary':'constant', 'constant':0.0, 'fluxconserve':'yes',
              'wcsinherit':'yes'}

    wregister(**params)


def main():
    parser = argparse.ArgumentParser(description=
        "Align combined images using IRAF's wregister task")
    parser.add_argument('wispfield', type=str, nargs=1,
        help='WISP field for which to align images. Must match ' +\
             'name of directory containing all relevant FITS files')
    args = parser.parse_args()
    # if user added a trailing slash to wispfield, 
    # remove for creating file names
    wispfield = args.wispfield[0].strip('/')

    # get list of combined images for this field
    images = [x for x in [os.path.join(wispfield,'%s_g.fits'%wispfield), \
               os.path.join(wispfield,'%s_i.fits'%wispfield)] \
               if x in glob(os.path.join(wispfield,'%s_*.fits'%wispfield))]
    images.sort()
    if len(images) == 1:
        # if there is only 1 combined image, no need to align
        print '\nOnly 1 combined image for %s (%s). Exiting.' % \
              (wispfield, images[0])
        exit()
    image = os.path.join(wispfield,'%s_g.fits'%wispfield)
    reference = os.path.join(wispfield,'%s_i.fits'%wispfield)

    # run wregister in IRAF
    run_wregister(image, reference)


if __name__ == '__main__':
    main()
