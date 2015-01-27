#! /usr/bin/env python
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
import astropy.units as u



def match_cats(RA, Dec, refRA, refDec):
    '''Match a catalog of RA's and Dec's to a reference catalog
    (refRA and refDec).
    
    Return the indices of the reference catalog that match each
    source in the input catalog, and the on-sky
    separation between each source's closest match in arcsec'''
    # create SkyCoord objects to use with the matching
    SC_cat = SkyCoord(ra=RA, dec=Dec, frame='icrs', unit=(u.deg,u.deg))
    SC_refcat = SkyCoord(ra=refRA, dec=refDec, frame='icrs',unit=(u.deg,u.deg))

    #    idx - indices of matched sources in reference cat
    #  sep2d - on-sky angular separation between closest match
    # dist3d - 3D distance between closest matches
    idx,sep2d,dist3d = match_coordinates_sky(SC_cat, SC_refcat)
    # convert separation to arcsecs
    separc = sep2d * u.arcsec
    return (idx,separc)


def main():
    match_cats(RA, Dec, refRA, refDec)


if __name__ == '__main__':
    main()
