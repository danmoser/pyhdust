#!/usr/bin/env python
# -*- coding:utf-8 -*-

""" Program to extract 2D images from a (3D) FITS cube
"""

import sys
from argparse import ArgumentParser
import os
import pyfits as pf

__version__ = "0.91"
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
parser.add_argument("INPUT", help=("Filename of the 3D FITS cube"), type=str)

args = parser.parse_args()


if __name__ == '__main__':

    fname = args.INPUT
    fits = pf.open(fname)
    nimgs = len(fits[0].data)
    for i in range(nimgs):
        hdu = pf.PrimaryHDU(fits[0].data[0])
        hdulist = pf.HDUList([hdu])
        hdulist[0].header = fits[0].header
        path, oname = os.path.split(fname)
        pref, ext = os.path.splitext(oname)
        hdu.writeto(os.path.join(path, pref+"_{0:04d}".format(i)+ext), 
            clobber=True)
    print('# Saved {0} 2D images from {1}'.format(nimgs, oname))
