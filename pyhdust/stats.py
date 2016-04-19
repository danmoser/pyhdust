#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
TO compara SUMMARY:

for i in range(8):
    a = _np.random.randn(10**i)+2
    print(_np.average(a), _np.std(a), summary(a))

"""

import numpy as _np


def mad(data, axis=None):
    """ Return 1.48xMAD (median absolute deviation) """
    return 1.48*_np.median(_np.abs(data - _np.median(data, axis)), axis)


def summary(x, verbose=False):
    """ Returns the summary of the variable: "median", "minus sigma" and 
    "plus sigma" ROBUST values (i.e., median and [15.9, 84.1] percentiles). 
    """
    data = _np.hstack((_np.median(x), _np.percentile(x, (15.87, 84.13))))
    if verbose:
        print('# median and [15.9, 84.1] percentiles: ')
        print('# {0} {1} {2}'.format(*data))
    return data
