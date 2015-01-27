#! /usr/bin/env python
#######################################################################
## Micaela Bagley
## June 2014
## 
## usage: astrometry.py [-h] [--useSE] wispfield
##
## Run all images in directory wispfield through astrometry.net
##
## positional arguments:
##   wispfield   WISP field to be fit with a WCS, must match name of directory
##               containing all relevant FITS files
##
## optional arguments:
##   -h, --help  show this help message and exit
##   --useSE    Use SExtractor instead of image2xy.py
#######################################################################
import argparse
import numpy as np
from glob import glob
import pyfits
import subprocess
import os
import shutil

def get_images(wispfield):
    '''Get list of images that have yet to be solved.'''
    # check that astrometry directory exists
    if not os.path.isdir(os.path.join(wispfield, 'astrometric_solns')):
        # create dir for astrometry.net outputs if it does not exist
        os.mkdir(os.path.join(wispfield, 'astrometric_solns'))
        # then all files in directory need to be solved
        images = [x for x in glob(os.path.join(wispfield, '*.fits')) if not \
                  os.path.basename(x).startswith(('result', 'Bias', 'Zero'))]

    else:
        # any files that have already been solved
        complete = [os.path.splitext(os.path.basename(x))[0]+'_solved' for x \
            in glob(os.path.join(wispfield, 'astrometric_solns', '*.solved'))]
        # get list of images that have yet to be solved
        images = [x for x in glob(os.path.join(wispfield, '*.fits')) if not \
                  os.path.basename(x).startswith(('result','Bias','Zero')) and \
                  os.path.splitext(os.path.basename(x))[0] not in complete]
    images.sort()
    return images


def sex_to_deg(RA, Dec):
    '''Convert a sexigessimal RA,Dec into decimal RA,Dec'''
    RA = map(float, RA.split(':'))
    Dec = map(float, Dec.split(':'))
    decimalRA = 15.*(RA[0] + (RA[1] + RA[2]/60.)/60.)
    decimalDec = np.sign(Dec[0])*(np.abs(Dec[0]) + (Dec[1] + Dec[2]/60.)/60.)
    return decimalRA, decimalDec


def get_RADec(image):
    '''Get RA,Dec from image header'''
    hdr = pyfits.getheader(image)
    RA = hdr['RA']
    Dec = hdr['DEC']
    decimalRA,decimalDec = sex_to_deg(RA, Dec)
    return decimalRA, decimalDec


def solve_image(image, RA, Dec, useSE=False):
    '''Run astrometry.net on an image using the RA,Dec from the
       header as an initial guess.

       --no-fits2fits skips running the FITS files through a
       "sanitizer" which tries to clean up non-standards-compliant
       images. Our FITS files are already compliant.

       A low and high pixel scale is given (0.1 - 0.45"/pix). 
       This should cover both binned and unbinned images.

       Set --useSE to use sextractor to detect sources rather than
       astrometry.net's bundled "images2xy" program. 
       Note: if sextractor is used, the xyls output files may be
       incorrectly all (0,0). 

       Solved files are called [base]_solved.fits. Several other
       output files are created. See below for filenames and descriptions
    '''
    new = ''.join([os.path.splitext(image)[0], '_solved.fits'])
    if useSE:
        cmd = 'solve-field --no-plots --no-fits2fits --use-sextractor ' + \
              '--scale-units app -L 0.1 -H 0.45 --crpix-center --radius 1 ' + \
              '--ra %.4f --dec %.4f --new-fits %s %s' % (RA, Dec, new, image)
    else:
        cmd = 'solve-field --no-plots --no-fits2fits --scale-units app ' + \
              '-L 0.1 -H 0.45 --crpix-center --radius 1 ' + \
              '--ra %.4f --dec %.4f --new-fits %s %s' % (RA, Dec, new, image)
    print cmd
    subprocess.call(cmd, shell=True)
    return



def main():
    parser = argparse.ArgumentParser(description=
        'Run all images in directory wispfield through astrometry.net')
    parser.add_argument('wispfield', type=str, nargs=1,
        help='WISP field to be fit with a WCS, must match name of ' +\
             'directory containing all relevant FITS files')
    parser.add_argument('--useSE', action="store_true",
        help='Use SExtractor instead of image2xy.py')
    args = parser.parse_args()
    wispfield = args.wispfield[0].strip('/')

    # get list of images that need to be solved
    images = get_images(wispfield)
    for image in images:
        RA,Dec = get_RADec(image)
        solve_image(image, RA, Dec, useSE=args.useSE)
        if os.path.isfile(os.path.splitext(image)[0]+'.solved'):
            # move all but solved image to astrometric_solns directory
            soln_files = [x for x in glob(os.path.join(wispfield, '*')) if \
                          os.path.splitext(x)[0] == os.path.splitext(image)[0]]
            soln_files.append(os.path.splitext(image)[0] + '-indx.xyls')
            for f in soln_files:
                shutil.move(f, os.path.join(wispfield, 'astrometric_solns'))
    
    print '\nCheck your /tmp directory'
    print 'Delete any tmp.sanitized.xxxx files in your /tmp directory.\n'


if __name__ == '__main__':
    main()


'''
   Astrometry.net Default Outputs:
    <base>.wcs       - a FITS WCS header for the solution
    <base>.new       - a new FITS file containing the WCS header
    <base>-indx.xyls - a FITS BINTABLE with the pixel locations of
                       stars from the indx
    <base>.rdls      - a FITS BINTABLE with the Ra,Dec of sources we
                       extracted from the image
    <base>.axy       - a FITS BINTABLE of the sources we extracted, plus
                       headers that describe the job (how the image is
                       going to be solved)
    <base>.solved    - exists and contains (binary) 1 if the field solved
    <base>.match     - a FITS BINTABLE describing the quad match that 
                       solved the image
'''
