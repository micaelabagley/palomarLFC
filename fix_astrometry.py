#! /usr/bin/env python
#######################################################################
## Micaela Bagley
## September 2014
##
## usage: fix_astrometry.py [-h] wispfield
## 
## Fix the astrometry on Palomar images using IRAF's ccmap task
##
## positional arguments:
##   wispfield   WISP field for which to construct a catalog. Must 
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
from astropy.table import Table
from match_cats import match_cats
from run_SE import run_SE
from pyraf import iraf


def read_cat(catfile):
    '''Read in fits tables'''
    f = fits.open(catfile)
    cat = f[1].data
    f.close()
    return cat


def run_ccmap(Palcat, ims, threshold, wispfield):
    '''Match Palomar and SDSS catalogs for fixing astrometry
       with IRAF task CCMAP
    '''
    # read in SE catalog for Palomar field
    pal = read_cat(Palcat)
    palRA = pal.field('X_WORLD')
    palDec = pal.field('Y_WORLD')
    palx = pal.field('X_IMAGE')
    paly = pal.field('Y_IMAGE')

    # read in SDSS catalog 
    sdss = read_cat(os.path.join(wispfield,'result.fits'))
    sdssRA = sdss.field('ra')
    sdssDec = sdss.field('dec')

    # match Palomar to SDSS 
    idx,separc = match_cats(palRA, palDec, sdssRA, sdssDec)
    match = (separc.value*3600. <= threshold)  
    print '\n%i Palomar objs matched to SDSS objs\n' % idx[match].shape[0]
    
    # create coordinate file for ccmap
    t = Table([palx[match], paly[match], 
              sdssRA[idx[match]], sdssDec[idx[match]]],
              names=('xcolumn', 'ycolumn', 'lngcolumn', 'latcolumn'))

    # print matched catalog to file
    t.write('SDSS.match', format='ascii.no_header')

    # set up parameters for ccmap
    params = {'input':'SDSS.match', 'database':'SDSS.db', 'images':ims, \
        'xcol':1, 'ycol':2, 'lngcol':3, 'latcol':4, 'lngunit':'deg', \
        'latunit':'deg', 'insyste':'j2000', 'refpoin':'coords', \
        'lngref':'CRVAL1',
        'latref':'CRVAL2', 'refsyst':'j2000', 'lngrefu':'CUNIT1', \
        'latrefu':'CUNIT2', 'project':'tan', 'fitgeom':'general', \
        'xxorder':5, 'xyorder':5, 'yxorder':5, 'yyorder':5, 'update':'yes'}
        #'xref':'CRPIX1', 'yref':'CRPIX2', 'lngref':'CRVAL1', \
    # run ccmap
    from iraf import images
    images.imcoords.ccmap(**params)


def main():
    parser = argparse.ArgumentParser(description=
        "Fix the astrometry on Palomar images using IRAF's ccmap task")
    parser.add_argument('wispfield', type=str, nargs=1,
        help='WISP field for which to fix the astrometry. ' +\
             'Must match name of directory containing all relevant files')
    args = parser.parse_args()
    wispfield = args.wispfield[0].strip('/')

    # get list of images for this field
    images = [x for x in glob(os.path.join(wispfield, '*_solved.fits'))]
    images.sort()

    print '\nFixing astrometry on images: %s\n' % images
    # Run SE on all Palomar images
    for image in images:
        run_SE([image], 'Astrometry')
    
    catalogs = glob(os.path.join(wispfield, '*astr_cat.fits'))
    catalogs.sort()
   
    # set up and run ccmap to update astrometry in all headers
    threshold = 3.5  # arcsec for matching radius
    for i,v in enumerate(catalogs):
        print images[i], catalogs[i]
        run_ccmap(catalogs[i], images[i], threshold, wispfield)



if __name__ == '__main__':
    main()


