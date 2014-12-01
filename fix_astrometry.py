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
import numpy as np
from glob import glob
import pyfits
import subprocess,os,shutil
from astropy.coordinates import SkyCoord,match_coordinates_sky
import astropy.units as u
from astropy.table import Table
from astropy.io import ascii
from utils.match_cats import match_cats
from pyraf import iraf


def run_SE(image, wispfield):
    '''Run SExtractor on the images'''
    # names of output files
    base = os.path.splitext(image)[0]
    cat = base + '_astr_cat.fits'
    seg = base + '_astr_seg.fits'

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


def read_cat(catfile):
    '''Read in fits tables'''
    f = pyfits.open(catfile)
    cat = f[1].data
    f.close()
    return cat


def make_reg(RA, Dec, reg, radius, color, width=1):
    '''Add circlular regions in fk5 coords to a ds9 region file'''
    for r,d in zip(RA, Dec):
        reg.write('circle(%.6f,%.6f,%f") # color=%s width=%i\n' % \
                 (r,d,radius,color,width))


def run_ccmap(Palcat, ims, threshold, wispfield):
    '''Match Palomar and SDSS catalogs for fixing astrometry
       with IRAF task CCMAP
    '''
    # region file of Palomar, SDSS, and matched objects
    reg = open('Palomar-SDSS.reg', 'w')
    reg.write('global color=green dashlist=8 3 width=1 font="helvetica 10 '
              'normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 '
              'move=1 delete=1 include=1 source=1 \n')
    reg.write('fk5 \n')

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

    # add all Palomar sources to region file
    make_reg(palRA, palDec, reg, 2, 'blue')
    # add all SDSS sources to region file
    make_reg(sdssRA, sdssDec, reg, 1, 'red')

    # match Palomar to SDSS 
    idx,separc = match_cats(palRA, palDec, sdssRA, sdssDec)
    match = (separc.value*3600. <= threshold)  
    print '\n%i Palomar objs matched to SDSS objs\n' % idx[match].shape[0]
    
    # add successfully matched sources to region file with width=4
    make_reg(sdssRA[idx[match]], sdssDec[idx[match]], reg, 1, 'red', width=4)
    reg.close()

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
        'xref':'CRPIX1', 'yref':'CRPIX2', 'lngref':'CRVAL1', \
        'latref':'CRVAL2', 'refsyst':'j2000', 'lngrefu':'CUNIT1', \
        'latrefu':'CUNIT2', 'project':'tan', 'fitgeom':'general', \
        'xxorder':5, 'xyorder':5, 'yxorder':5, 'yyorder':5, 'update':'yes'}
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
        run_SE(image, wispfield)
    
    catalogs = glob(os.path.join(wispfield, '*astr_cat.fits'))
    catalogs.sort()
   
    # set up and run ccmap to update astrometry in all headers
    threshold = 3.5  # arcsec for matching radius
    for i,v in enumerate(catalogs):
        print images[i], catalogs[i]
        run_ccmap(catalogs[i], images[i], threshold, wispfield)



if __name__ == '__main__':
    main()


