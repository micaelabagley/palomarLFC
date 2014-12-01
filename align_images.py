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
import pyfits
from pyraf.iraf import wregister


def run_wregister(wispfield):
    '''Align combined i image to combined g image using the WCS in 
       the header and IRAF's WREGISTER task
    '''
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

    # rename old i band image so it is not overwritten
    os.rename(os.path.join(wispfield,'%s_i.fits'%wispfield), 
              os.path.join(wispfield,'%s_i_old.fits'%wispfield))
    images.append(os.path.join(wispfield,'%s_i_old.fits'%wispfield))
    
    params = {'input':'%s[0]'%images[2],
              'reference':'%s[0]'%images[0],
              'output':images[1], 
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

    # get all scripts and files set up for running wregister in IRAF
    run_wregister(wispfield)


if __name__ == '__main__':
    main()
