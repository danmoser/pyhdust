#!/usr/bin/env python
# -*- coding:utf-8 -*-

""" Program to normalize a FITS (*.fit*) or a 2-columns text file (other 
extensions) spectrum based on non-parametric fitting.
"""

from glob import glob
import numpy as np
import os
import sys
import pyfits as pf
import pyhdust.spectools as spt
from pyhdust.tabulate import tabulate as tab
from argparse import ArgumentParser

__version__ = "0.93"
__author__ = "Daniel Moser"
__email__ = "dmfaes@gmail.com"


class MyParser(ArgumentParser): 
    def error(self, message):
        sys.stderr.write('# ERROR! %s\n' % message)
        self.print_help()
        sys.exit(2)

parser = MyParser(description=__doc__)
parser.add_argument('--version', action='version', 
    version='%(prog)s {0}'.format(__version__))
parser.add_argument("INPUT", help=("String to retrive the spectrum(a) "
    "filename(s) (wildcards accepted)"))
# parser.add_argument("-l", "--line", action="store", dest="line", 
#     help=("Wavelength of the reference line [default: %(default)s]"), 
#     type=float, default=6562.79)

args = parser.parse_args()


def norm_text(fname):
    f0 = np.loadtxt(fname)
    ny = spt.normalize_spec(f0[:, 0], f0[:, 1])
    outname = '{0}_norm{1}'.format(*os.path.splitext(fname))
    # f0 = np.savetxt(outname, np.column_stack((f0[:, 0], ny)))
    f1 = open(outname, 'w')
    f1.writelines( tab(np.column_stack((f0[:, 0], ny)), tablefmt='plain') )
    print('# {0} saved!'.format(outname))
    f1.close()
    return


def norm_spec(fname):
    imfits = pf.open(fname)
    hdr = imfits[0].header
    flux = imfits[0].data
    wl = np.arange(len(flux)) * hdr['CDELT1'] + hdr['CRVAL1']
    nflux = spt.normalize_spec(wl, flux.astype(float))
    outname = '{0}_norm{1}'.format(*os.path.splitext(fname))
    hdu = pf.PrimaryHDU(nflux)
    hdulist = pf.HDUList([hdu])
    hdulist[0].header = hdr
    hdu.writeto(outname, clobber=True)
    print('# {0} saved!'.format(outname))
    imfits.close()
    return

if __name__ == '__main__':

    linput = args.INPUT

    for specfile in glob(linput):
        if specfile.find('_norm.') > 0:
            print('# Warning! Remove "_norm" from {0} if you want to normalize'
                ' it!'.format(specfile))
            continue
        if specfile.find('.fit') > 0:
            norm_spec(specfile)
        else:
            norm_text(specfile)

    if len(glob(linput)) == 0:
        print('# ERROR! {0} files were not found!'.format(linput))
