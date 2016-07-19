#! /usr/bin/env python
import argparse
from astropy.table import Table
from astropy.io import fits
from glob import glob
import os

def extract_header(imlist):
    '''Extract image type, target, exposure time, filter, airmass, and 
       the binning of the array'''
 
    outfile = os.path.join(os.path.dirname(imlist[0]), 'list.txt')
    imlist.sort()

    # set up table
    names = ('image', 'imtype', 'obj', 'exptime', 'filt', 'airmass', 'binning')
    datatype = ('S40', 'S20', 'S20', 'float', 'S5', 'float', 'S10')
    t = Table(data=None, names=names, dtype=datatype)

    for i,image in enumerate(imlist):
        hdr = fits.getheader(image)
        imtype = hdr['IMAGETYP']
        obj = hdr['OBJECT']
        exptime = hdr['EXPTIME']
        filt = hdr['FILTER']
        try:
            airmass = hdr['AIRMASS']
        except:
            airmass = 1.0
        nx1 = hdr['NAXIS1']
        nx2 = hdr['NAXIS2']
        bin1 = hdr['CCDBIN1']
        bin2 = hdr['CCDBIN2']
        if (bin1,bin2) == (2,2):
            binning = 'binned'
        elif (bin1,bin2) == (1,1):
            binning = 'unbinned'
        else:
            print image + ' does not match known size for either binned ' +\
                'or unbinned data '
            exit()
        im = os.path.basename(image)
        t.add_row([im, imtype, obj, exptime, filt, airmass, binning])
    t.write(outfile, format='ascii.tab') #, delimiter='\t')
    return t

    
def sort_images(imlist):
    '''Sort images into biases, flats, science images'''
    


def sort_fields(imlist):
    '''Sort the reduced data by 
    '''
'''
# how to find field numbers form header info where obj name may be
'wisps33', 'wisp33', 'par33', etc.
regular expressions:
import re
check = re.search('\d+', objname)
# \ either escapes special characters or signals a special sequence
# d special sequence: matches decimal digit
# +  causes the resulting RE to match 1 or more repetitions 
parnum = check.group(0)

'''
def main():
    parser = argparse.ArgumentParser(description='Sort images by image type.')
    parser.add_argument('--dir', type=str, default='.', 
        help='Directory in which      images are located')
    args = parser.parse_args()
    print args.dir

    imlist = glob(os.path.join(args.dir, '*.fits'))
    # extract info from headers
    t = extract_header(imlist)


if __name__ == '__main__':
    main()




