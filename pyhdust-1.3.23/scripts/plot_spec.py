#!/usr/bin/env python
# -*- coding:utf-8 -*-

""" Program to export as figure line profiles from standard FITS spectra.
"""

import sys
from glob import glob
import pyhdust.spectools as spt
import pyhdust.phc as phc
from argparse import ArgumentParser
import os
import matplotlib.pyplot as plt

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
parser.add_argument("INPUT", help=("String to retrive the spectrum(a) "
    "filename(s) (wildcards accepted wrapped with \"\")"), type=str)
parser.add_argument("-l", "--line", action="store", dest="line", 
    help=("Wavelength of the reference line [default: %(default)s]"), 
    type=float, default=6562.79)
parser.add_argument("-w", "--hwidth", action="store", dest="hwidth", 
    help=("Half-width of the line profile (in km/s) [default: %(default)s]"), 
    type=float, default=1000.)

# group2 = parser.add_argument_group('other arguments (overwrite power)')
# group2.add_argument("-r", "--remove", action="store_true", dest="rm", 
#     help=("If this flag is enabled, it restores the spectrum(a) to the "
#         "original value (i.e., VELSHIFT = 0)"), default=False)
# group2.add_argument("-f", "--force-dv", action="store", dest="vs", 
#     help=("Force VELSHIFT to the non-zero VS value (in km/s) "
#         "[default: None]"), type=float, default=0.)

args = parser.parse_args()


if __name__ == '__main__':

    lbc = args.line
    linput = args.INPUT
    hwidth = args.hwidth

    plev = 0.
    lfiles = glob(linput)
    lfiles.sort()
    if len(lfiles):
        fig, ax = plt.subplots(figsize=(8, 6*(1+len(lfiles)/25.)))

    for fitsfile in lfiles:
        sdata = spt.loadfits(fitsfile)
        vl, fx = spt.lineProf(sdata[0], sdata[1], lbc=lbc, hwidth=hwidth)
        ax.plot(vl, fx+plev)
        ax.text(0.95*hwidth, 1.+plev+0.05/2, os.path.splitext(fitsfile)[0], 
            horizontalalignment='right',  # transform=ax.transAxes,
        verticalalignment='center', fontsize=8)
        if plev == 0:
            ymin = min(fx)
        plev += 0.04

    if len(glob(linput)) == 0:
        print('# ERROR! {0} files were not found!'.format(linput))
    else:
        # ax.legend(loc='best', fancybox=True, framealpha=0.5, fontsize=8, 
        #         labelspacing=0.05)
        ax.set_ylim([ymin, plev+max(fx)])
        ax.set_xlim([-hwidth, hwidth])
        phc.savefig(fig, figname='plot_spec_'+phc.dtflag())
