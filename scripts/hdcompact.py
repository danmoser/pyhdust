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

extensions = ['.tau', '.map']

path = "."
if len(sys.argv) == 2:
    path = sys.argv[1]

for r, d, fnames in os.walk(path):
    for f in fnames:
        fp = os.path.join(r, f)
        ext = os.path.splitext(fp)[1]
        if any(ext.startswith(e) for e in extensions):
            fpnew = fp+'.bz2'
            # fpnew = fp+'.gz'
            print('# Creating '+fpnew)
            # with gzip.open(fpnew, 'wb') as fout:
            with bz2.BZ2File(fpnew, 'wb', compresslevel=1) as fout:
                fout.write(open(fp).read())
            if os.path.exists(fpnew):
                if os.path.getsize(fpnew) > 40000:
                    print('# Removing '+fp)
                    os.remove(fp)
