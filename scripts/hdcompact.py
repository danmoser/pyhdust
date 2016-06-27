#!/usr/bin/env python
# -*- coding:utf-8 -*-

""" Program to recursively compact *.tau* and *.map* files (from current 
folder).

If the compacted files have more than 40 KB (i.e., the size of compacted 
header), then the original file is automatically excluded. 
"""
import bz2
# import gzip
import os
import sys
from multiprocessing import Pool
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
parser.add_argument("-n", action="store", dest="n", 
    help=("Number of cores to be used [default: %(default)s]"), 
    type=float, default=2)

group2 = parser.add_argument_group('other arguments')
group2.add_argument("-c", "--check", action="store_true", dest="chk", 
    help=("If this flag is enabled, it only removed the *.map* file if the "
        "*.map.bz2 file exists"), default=False)

args = parser.parse_args()


def compact(fp):
    fpnew = fp+'.bz2'
    # fpnew = fp+'.gz'
    print('# Creating '+fpnew)
    # with gzip.open(fpnew, 'wb') as fout:
    with bz2.BZ2File(fpnew, 'wb', compresslevel=1) as fout:
        fout.write(open(fp).read())
    if os.path.exists(fpnew):
        if os.path.getsize(fpnew) > 40000:
            os.remove(fp)
    return '# Removing '+fp

if __name__ == '__main__':
    extensions = ['.tau', '.map']

    path = "."
    if len(sys.argv) == 2:
        path = sys.argv[1]

    clist = []
    for r, d, fnames in os.walk(path):
        for f in fnames:
            fp = os.path.join(r, f)
            ext = os.path.splitext(fp)[1]
            if any(ext.startswith(e) for e in extensions):
                if not args.chk:
                    clist.append(fp)
                else:
                    fpnew = fp+'.bz2'
                    if os.path.exists(fpnew):
                        if os.path.getsize(fpnew) > 40000:
                            os.remove(fp)
                            print('# Removing '+fp)

    if not args.chk:
        pool = Pool(processes=args.n)
        for result in pool.imap_unordered(compact, clist):
            print(result)
