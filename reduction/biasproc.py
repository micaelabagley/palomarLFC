#! /usr/bin/env python
import argparse
from astropy.io import fits
import os
from glob import glob
import numpy as np
import time


def bias_combine(biases, binning):
    '''Median combine the biases. Requires >3 biases.

       Keep track of image binning. Biases are either binned 2x2 
       or unbinned (1x1). 
       
       Output Master bias frame will be Bias_[binning].fits
    '''
    # are there enough biases for median combining?
    nbias = len(biases)
    if nbias < 3:
        print 'There are not enough biases to median combine ' +\
            '(fewer than 3 biases) '
        exit()

    # filename of output Master Bias frame   
    output = os.path.join(os.path.dirname(biases[0]), 'Bias_%s.fits' % binning)

    # get basic header information
    bhd = fits.getheader(biases[0])
    nx = bhd['NAXIS1']
    ny = bhd['NAXIS2']
    median_array = np.zeros((nbias,ny,nx), dtype=float)

    for i,im in enumerate(biases):
        bb,bhd = fits.getdata(im, header=True)
        median_array[i,:,:] = bb

    Bias = np.median(median_array, axis=0)

    # add keywords to header
    bhd['OBJECT'] = 'Master Bias'
    bhd['NCOMBINE'] = (nbias, 'Number of biases combined to form master bias')
    # write out master bias
    fits.writeto(output, Bias, header=bhd, clobber=True)
    
    return output


def bias_subtract(imlist, MasterBias, SaveSteps=False):
    '''Bias subtract the images in imlist using the previously
       combined Master Bias frame.

       Set SaveSteps = True to write each bias-subtracted image 
       to a new file named [input_image].bs.fits. This option is 
       good for checking the reduction at each step. 
       Default is to overwrite input image files.
    '''
    # read in Master Bias frame
    Bias = fits.getdata(MasterBias)

    # get the date and time
    now = time.strftime('%c')
    for image in imlist:
        im,hdr = fits.getdata(image, header=True)
        # write a keyword to header 
        biasstr = 'Bias subtracted: %s' % now
        hdr['BIASSUB'] = biasstr

        # output file name
        if SaveSteps:
            output = ''.join([os.path.splitext(image)[0], '.bs.fits'])
        else:            
            output = image

        new = im - Bias
        fits.writeto(output, new, header=hdr, clobber=True)
    
    return



def main():
    pass


if __name__ == '__main__':
    main()

