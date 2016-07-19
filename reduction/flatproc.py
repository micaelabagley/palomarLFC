#! /usr/bin/env python
import argparse
from astropy.io import fits
import os
from glob import glob
import numpy as np
import time


def flat_combine(flats, filt, binning):
    '''Median combine the flats for filter 'filt'. Requires >3 flats.

       Keep track of image binning. Flats are either binned 2x2 
       or unbinned (1x1). 
       
       Output Master flat frame will be Flat[filt]_[binning].fits
    '''
    # are there enough flats for median combining?
    nflats = len(flats)
    print flats
    print nflats
    if nflats < 3:
        print 'There are not enough flats to median combine ' +\
            '(fewer than 3 flats) '
        exit()

    # filename of output Master flat frame   
    output = os.path.join(os.path.dirname(flats[0]), 
                          'Flat%s_%s.fits' % (filt,binning))

    # get basic header information
    fhd = fits.getheader(flats[0])
    nx = fhd['NAXIS1']
    ny = fhd['NAXIS2']
    median_array = np.zeros((nflats,ny,nx), dtype=float)

    for i,im in enumerate(flats):
        ff,fhd = fits.getdata(im, header=True)
        # check that median of image is != 0
        med = np.median(ff)
        if med == 0:
            print 'Median value for ' + im + ' is 0!'
            exit()
        # normalize by median value
        median_array[i,:,:] = ff / med

    Flat = np.median(median_array, axis=0)

    # add keywords to header
    fhd['OBJECT'] = 'Master %s Band Flat'%filt
    fhd['NCOMBINE'] = (nflats, 'Number of flats combined to form master flat')
    # write out master flat
    fits.writeto(output, Flat, header=fhd, clobber=True)

    return output


def flat_field(imlist, MasterFlat, SaveSteps=False):
    '''Flat field the images in imlist using the previously
       combined Master Flat frames.

       Set SaveSteps = True to write each flat-fielded image 
       to a new file named [input_image].bs.ff.fits. This option is 
       good for checking the reduction at each step. 
       Default is to overwrite input image files.
    '''
    # read in Master Flat frame
    Flat = fits.getdata(MasterFlat)

    # check the flat for zeros, replace with median of flat
    Flat[Flat == 0.] = np.median(Flat)
    Flat = Flat / np.median(Flat)

    # get the date and time
    now = time.strftime('%c')
    for image in imlist:
        im,hdr = fits.getdata(image, header=True)
        # write a keyword to header 
        flatstr = 'Flat fielded: %s' % now
        hdr['FLATFLD'] = flatstr

        # output file name
        if SaveSteps:
            output = ''.join([os.path.splitext(image)[0], '.ff.fits'])
        else:            
            output = image

        new = im / Flat
        fits.writeto(output, new, header=hdr, clobber=True)
    
    return



def main():
    pass


if __name__ == '__main__':
    main()

