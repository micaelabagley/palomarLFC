#! /usr/bin/env python
import pyfits
import pywcs
import numpy as np

def trim_image(image, cra, cdec, xsize, ysize, outname):
    '''Trim images and maintain WCS.

       Provide the center of the field and desired size
    '''
    im,hdr = pyfits.getdata(image, header=True)
    hdr_wcs = pywcs.WCS(hdr)

    # convert size to pixels
    xs = xsize / hdr['SECPIX1']
    ys = ysize / hdr['SECPIX2']

    # update header to have correct WCS
    # from Adam Ginsburg's cutout.py
    hdr['CRPIX1'] = xs / 2.
    hdr['CRPIX2'] = ys / 2.
    hdr['CRVAL1'] = cra
    hdr['CRVAL2'] = cdec
    hdr['NAXIS1'] = xs
    hdr['NAXIS2'] = ys

    # find center of WISP field in Palomar pixel coordinates
    pix = hdr_wcs.wcs_sky2pix(np.array([cra,cdec]),1)
    x1 = pix[0] - xs/2.
    x2 = pix[0] + xs/2.
    y1 = pix[1] - ys/2.
    y2 = pix[1] + ys/2.

    new = im[y1:y2,x1:x2]

    pyfits.writeto(outname, new, header=hdr, clobber=True)

