#! /usr/bin/env python
import argparse
import pyfits
import os
from glob import glob
import numpy as np
import time


def dark_combine(darks, binning):
    '''Median combine the darks. Requires >3 darks.

       Keep track of image binning. Darks are either binned 2x2 
       or unbinned (1x1). 
       
       Output Master dark frame will be Dark_[binning].fits
    '''
    # are there enough darks for median combining?
    ndarks = len(darks)
    if ndarks < 3:
        print 'There are not enough darks to median combine ' +\
            '(fewer than 3 darks) '
        exit()

    # filename of output Master Dark frame   
    output = os.path.join(os.path.dirname(darks[0]), 'Dark_%s.fits' % binning)


    # get basic header information
    dhd = pyfits.getheader(darks[0])
    nx = dhd['NAXIS1']
    ny = dhd['NAXIS2']
    median_array = np.zeros((ndarks,ny,nx), dtype=float)

    for i,im in enumerate(darks):
        dd,dhd = pyfits.getdata(im, header=True)
        exptime = dhd['EXPTIME']
        median_array[i,:,:] = dd / exptime

    Dark = np.median(median_array, axis=0)

    # add keywords to header
    dhd['OBJECT'] = 'Master Dark'
    dhd['NCOMBINE'] = (ndarks, 'Number of darks combined to form master dark')
    # write out master dark
    pyfits.writeto(output, Dark, header=dhd, clobber=True)

    return output


def dark_subtract(imlist, MasterDark, SaveSteps=False):
    '''Dark subtract the images in imlist using the previously
       combined Master Dark frame. Scale the Master Dark up to 
       the exptime of each image before subtraction.

       Set SaveSteps = True to write each dark-subtracted image 
       to a new file named [input_image].ds.fits. This option is 
       good for checking the reduction at each step. 
       Default is to overwrite input image files.
    '''
    # read in Master Dark frame
    Dark = pyfits.getdata(MasterDark)

    # get the date and time
    now = time.strftime('%c')
    for image in imlist:
        im,hdr = pyfits.getdata(image, header=True)
        # write a keyword to header 
        darkstr = 'Dark subtracted: %s' % now
        hdr['DARKSUB'] = darkstr

        # output file name
        if SaveSteps:
            output = ''.join([os.path.splitext(image)[0], '.ds.fits'])
        else:            
            output = image

        # scale Dark up to exptime of image
        exptime = hdr['EXPTIME']
        ScaledDark = Dark * exptime
        new = im - ScaledDark
        pyfits.writeto(output, new, header=hdr, clobber=True)
    
    return



def main():
    pass


if __name__ == '__main__':
    main()

