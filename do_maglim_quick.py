#! /usr/bin/env python
import argparse
import numpy as np


def run_maglim(wispfield):
    # get list of images for this field
    images = [x for x in [os.path.join(wispfield,'%s_g.fits'%wispfield), \
               os.path.join(wispfield,'%s_i.fits'%wispfield)] \
               if x in glob(os.path.join(wispfield,'%s_*.fits'%wispfield))]
    images.sort()

    # read in calibration info for this field
    #   if maglim has already been calculated, there will be 4 columns
    test = np.genfromtxt(os.path.join(wispfield, 'sdss_calibration.dst'))
    if test[0].shape[0] == 3:
        # open file to get lines of comments
        f = open(os.path.join(wispfield, 'sdss_calibration.dat'), 'r')
        lines = f.readlines()
        f.close()
        # open file for writing
        f = open(os.path.join(wispfield, 'sdss_calibration.dat'), 'w')
        # write out all but the last 3 lines
        for line in lines[:-3]:
            f.write(line)
        f.write('# \n')
        f.write('# Filter   alpha[0]   alpha[1]    maglim \n')
        # get filters, zero points and color terms
        cal = np.genfromtxt(os.path.join(wispfield, 'sdss_calibration.dat'),
            dtype=[('filts','S10'), ('cterm',float), ('zp',float)])
        filts = cal['filts']
        cterm = cal['cterm']
        zp = cal['zp']
        # add in magnitude limit
        for image in images:
            # image filter
            filt = fits.getheader(image)['FILTER'].rstrip("'")
            # get 1 sigma limiting magnitude
            if len(images) == 1:
                segmap = os.path.splitext(image)[0] + '_calib_seg.fits'
            elif len(image) == 2:
                segmap = os.path.join(wispfield,'%s_i_calib_seg.fits'%wispfield)
            maglim = find_maglim(image, segmap, zp[filts==filt], wispfield)
            f.write('  %s   %f   %f   %f\n'%(filt, cterm[filts==filt], 
                                             zp[filts == filt], maglim))
        f.close()

    elif test[0].shape[0] == 4:
        print 'Maglim already done. Skipping'    
        exit()
    else:
        print '%i columns in sdss_calibration.dat?'%test[0].shape[0]
        

def main():
    parser = argparse.ArgumentParser(description=
        'Get the maglimits for a field and add them to sdss_calibration.dat')
    parser.add_argument('wispfield', type=str, nargs=1,
        help='WISP field for which to construct a catalog. ' +\
             'Must match name of directory containing all relevant files')
    args = parser.parse_args()
    wispfield = args.wispfield[0].strip('/')

    run_maglim(wispfield)

if __name__ == '__main__':
    main()
