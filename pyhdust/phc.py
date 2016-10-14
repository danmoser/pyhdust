# -*- coding:utf-8 -*-

"""PyHdust *phc* module: physical constants and general use functions

Includes functions for:
- Data manipulation (average, bin, fitting)
- Manipulation of 3D coordinates 
- Convolution functions
- Manipulation of strings and lists
- Manipulation of angular coordinates and dates 
- Files or directories manipulation
- Plotting 
- Physical functions
- `List of constants`_
- Python-version issues

.. _`List of constants`: phc_list.html


:license: GNU GPL v3.0 https://github.com/danmoser/pyhdust/blob/master/LICENSE
"""
from __future__ import print_function
import sys as _sys
import os as _os
import re as _re
import numpy as _np
import datetime as _dt
from dateutil.relativedelta import relativedelta as _dtdelta
import gzip as _gzip
from bz2 import BZ2File as _BZ2File
from glob import glob as _glob
from itertools import product as _product
from collections import Iterable as _It
import pyhdust.jdcal as _jdcal
from pyhdust.tabulate import tabulate as _tab
from six import string_types as _strtypes
import warnings as _warn
import struct as _struct

try:
    import matplotlib.pyplot as _plt
    from scipy import optimize as _optimize
except ImportError:
    ('matplotlib and/or scipy module not installed!!!')

__author__ = "Daniel Moser"
__email__ = "dmfaes@gmail.com"


# Data manipulation (includes fitting)
def find_nearest_pt(x0, y0, x, y, z, case=1):
    dx = _np.max(x) - _np.min(x)
    dy = _np.max(y) - _np.min(y)
    if case == 1:
        idx = _np.where( (x < x0) & (y < y0))
    elif case == 2:
        idx = _np.where( (x > x0) & (y < y0))
    elif case == 3:
        idx = _np.where( (x < x0) & (y > y0))
    elif case == 4:
        idx = _np.where( (x > x0) & (y > y0))
    else:
        idx = _np.argsort(x)
    z = z[idx]
    if len(z) > 0:
        x = x[idx]
        y = y[idx]
        d = _np.sqrt( ((x0-x)/dx)**2 + ((y0-y)/dy)**2 )
        idx = _np.where(d == _np.min(d))
        return _np.average(x[idx]), _np.average(y[idx]), _np.average(z[idx]), \
            _np.average(d[idx])
    else:
        return _np.NaN, _np.NaN, _np.NaN, _np.NaN


def baricent_calc(x0, y0, x, y, z, fullrange=False):
    x1, y1, z1, d1 = find_nearest_pt(x0, y0, x, y, z, case=1)
    x2, y2, z2, d2 = find_nearest_pt(x0, y0, x, y, z, case=2)
    x3, y3, z3, d3 = find_nearest_pt(x0, y0, x, y, z, case=3)
    x4, y4, z4, d4 = find_nearest_pt(x0, y0, x, y, z, case=4)
    if _np.sum(_np.isnan([x1, x2, x3, x4])) == 0:
        # DO THE BARICENTER CALCULATIONS HERE
        idx = _np.argsort([d1, d2, d3, d4])
        xs = _np.array([x1, x2, x3, x4])[idx]
        ys = _np.array([y1, y2, y3, y4])[idx]
        zs = _np.array([z1, z2, z3, z4])[idx]
        inside = False
        seq = _np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
        for i in seq:
            # http://mathworld.wolfram.com/TriangleInterior.html    
            # http://stackoverflow.com/questions/13300904/
            # determine-whether-point-lies-inside-triangle
            if not inside:
                X, Y = (xs[i], ys[i])
                alpha = ((Y[1] - Y[2])*(x0 - X[2]) + (X[2] - X[1])*(y0 - 
                    Y[2])) / ((Y[1] - Y[2])*(X[0] - X[2]) + (X[2] - 
                    X[1])*(Y[0] - Y[2]))
                beta = ((Y[2] - Y[0])*(x0 - X[2]) + (X[0] - X[2])*(y0 - 
                    Y[2])) / ((Y[1] - Y[2])*(X[0] - X[2]) + (X[2] - 
                    X[1])*(Y[0] - Y[2]))
                gamma = 1.0 - alpha - beta
                if alpha > 0 and beta > 0 and gamma > 0:
                    Z = zs[i]
                    inside = True
        # http://www.mathopenref.com/coordtrianglearea.html
        a = _np.zeros(3.)
        a[2] = _np.abs(0.5 * (X[0]*(Y[1]-y0) + X[1]*(y0-Y[0]) + 
            x0*(Y[0]-Y[1])))
        a[1] = _np.abs(0.5 * (X[0]*(Y[2]-y0) + X[2]*(y0-Y[0]) + 
            x0*(Y[0]-Y[2])))
        a[0] = _np.abs(0.5 * (X[1]*(Y[2]-y0) + X[2]*(y0-Y[1]) + 
            x0*(Y[1]-Y[2])))
        return _np.sum(Z*a)/_np.sum(a)
    elif _np.sum(_np.isnan([x1, x2, x3, x4])) == 1 and fullrange:
        idx = _np.argsort([d1, d2, d3, d4])
        zs = _np.array([z1, z2, z3, z4])[idx][:3]
        ds = _np.array([d1, d2, d3, d4])[idx][:3]
        return _np.sum(zs/ds) / _np.sum(1/ds)
    else:
        return _np.NaN


def baricent_map(x, y, z, res=100, fullrange=False):
    X = _np.linspace(_np.min(x), _np.max(x), res)
    Y = _np.linspace(_np.min(y), _np.max(y), res)
    X, Y = _np.meshgrid(X, Y[::-1])
    img = _np.zeros((res, res))
    for i, j in _product(range(res), range(res)):
        img[i, j] = baricent_calc(X[i, j], Y[i, j], x, y, z, 
            fullrange=fullrange)
    return img


def flatten(list_, outlist=None):
    """ Return a flatten list. 

    If ``outlist`` is provided, the output is extended to it.
    """
    if outlist is None:
        outlist = []
    for item in list_:
        if isinstance(item, _It) and not isinstance(item, _strtypes):
            flatten(item, outlist)
        else:
            outlist.append(item)
    return outlist


def wg_avg_and_std(values, sigma):
    """ Return the weighted average and standard deviation.

    This IS NOT the problem of Wikipedia > Weighted_arithmetic_mean >
        Weighted_sample_variance.

    average = _np.average(values, weights=weights)
    #Fast and numerically precise:
    variance = _np.average((values-average)**2, weights=weights)
    return (average, _np.sqrt(variance))

    INPUT: values, sigma -- arrays with the same shape.

    OUTPUT: average, avg_sig (float, float)
    """
    avg = _np.average(values, weights=1 / sigma)
    return (avg, _np.sqrt(_np.sum(sigma**2)) / len(values) )


def bindata(x, y, nbins=20, yerr=None, xlim=None, perc=0):
    """
    Return the weighted binned data.

    if ``perc > 0``, then it returns the percentile value of the interval.

    if yerr is None:
        yerr = _np.ones(shape=_np.shape(x))

    INPUT: x, y, err - arrays with the same shape (they don't need to be
    sorted); nbins=int, xlim=[xmin, xmax]

    OUTPUT: xvals, yvals, new_yerr (arrays)
    """
    x = _np.array(x)
    y = _np.array(y)
    if yerr is None:
        yerr = _np.ones(len(x))
    else:
        yerr = _np.array(yerr)
    nans, tmp = nan_helper(y)
    x = x[~nans]
    y = y[~nans]
    yerr = yerr[~nans]
    if xlim is None:
        xmax = _np.max(x)
        xmin = _np.min(x)
    else:
        xmin, xmax = xlim
    shift = (xmax - xmin) / (nbins - 1)
    tmpx = _np.arange(nbins) * shift + xmin
    tmpy = _np.zeros(nbins)
    tmpyerr = _np.zeros(nbins)
    idx = _np.zeros(nbins, dtype=bool)
    for i in range(nbins):
        selx = _np.where( abs(x - tmpx[i]) <= shift / 2. )
        if len(selx[0]) >= 1:
            if perc <= 0:
                tmpx[i] = _np.average(x[selx])
                tmpy[i], tmpyerr[i] = wg_avg_and_std(y[selx], yerr[selx])
            else:
                tmpy[i] = _np.percentile(y[selx], perc)
                pid = find_nearest(y[selx], tmpy[i], idx=True)
                tmpx[i] = x[selx][pid]
                _, tmpyerr[i] = wg_avg_and_std(y[selx], yerr[selx])
            idx[i] = True
    if _np.sum(yerr) / len(x) == 1:
        return tmpx[idx], tmpy[idx]
    else:
        return tmpx[idx], tmpy[idx], tmpyerr[idx]


def chi2calc(mod, obs, sig_obs=None, npar=1):
    """ Calculate the chi2 """
    nans, tmp = nan_helper(obs)
    obs = _np.array(obs)[~nans]
    mod = _np.array(mod)[~nans]
    if sig_obs is None:
        sig_obs = _np.ones(len(obs))
    else:
        sig_obs = _np.array(sig_obs)[~nans]
    return _np.sum( (mod - obs)**2 / sig_obs**2 ) / (len(sig_obs) - npar - 1)


def splitequal(n, N):
    """ Split `N` in approx. `n` igual parts 

    `N` must be integer.

    Suggestion: *phc.splitequal(N/8., N)* split N into sequences of 
    approx. 8 itens.
    """
    n = int(round(n))
    idx = []
    for i in range(n):
        idx.append([i * N / n, (i + 1) * N / n])
    if n == 0:
        idx = [[0, N]]
    return idx


def interLinND(X, X0, X1, Fx, disablelog=False):
    """
    N-dimensional linear interpolation in LOG space!!

    Pay attention: Fx must always be > 0. If not, put disablelog=True.

    Other important thing: Fx must be regularly spaced (i.e., [Fx0y0, Fx0y1,
    Fx1y0, Fx1y1] if X0=[x0,y0] and X1=[x1,y1]). 

    If it is not the case, see *interBar?D* function.

    | INPUT:
    | X = position in with the interpolation is desired;
    | X0 = minimal values of the interval;
    | X1 = maximum values of the inveral
    | Fx = function values along the interval, ORDERED BY DIMENSTION.
    | Example: Fx = [F00, F01, F10, F11]

    OUTPUT: interpolated value (float)"""
    X = _np.array(X)
    X0 = _np.array(X0)
    X1 = _np.array(X1)
    Xd = (X - X0) / (X1 - X0)
    DX = _np.array([ [(1 - x), x] for x in Xd ])
    #
    i = 0
    F = 0
    for prod in _product(*DX):
        if disablelog:
            F += Fx[i] * _np.product(prod)
        else:
            F += _np.log(Fx[i]) * _np.product(prod)
        i += 1
    #
    if not disablelog:
        return _np.exp(F)
    else:
        return F


def optim(p0, x, y, yerr, func, errfunc=None):
    """ Do scipy.optimize.leastsq minimization. 
    Default error function is the chi2 function.

    Requirements:
        - func(p, x) is a previously user defined function.
        - len(x)=len(y)=len(yerr)

    Output:
        best parameters (p), chi2_red value
    """
    if errfunc is None:
        def errfunc(p, x, y, yerr, func):
            """ error function """
            # return _np.sum( ((y-func(p,x))/yerr)**2 )
            return ((y - func(p, x)) / yerr)**2
    bestp, tmp = _optimize.leastsq(errfunc, p0, args=(x, y, yerr, func))
    c2red = chi2calc(func(bestp, x), y, yerr, npar=len(p0))
    return bestp, c2red


def optim2(p0, x, y, yerr, func):
    """ Do scipy.optimize.curve_fit minimization. 
    It returns errors to the parameters fitting!!!

    Requirements:
        - func(p, x) is a previously user defined function.
        - len(x)=len(y)=len(yerr)

    Output:
        best params (p), params errors (perr), chi2_red value
     """
    bestp, cov = _optimize.curve_fit(func, x, y, p0=p0, sigma=yerr)
    c2red = chi2calc(func(x, bestp), y, yerr, npar=len(p0))
    return bestp, _np.sqrt(_np.diag(cov)), c2red


def bin_ndarray(ndarray, new_shape, operation='avg'):
    """
    Bins an ndarray in all axes based on the target shape, by summing or
        averaging.
    Number of output dimensions must match number of input dimensions.

    Example
    -------
    >>> m = _np.arange(0,100,1).reshape((10,10))
    >>> n = bin_ndarray(m, new_shape=(5,5), operation='sum')
    >>> print(n)
    [[ 22  30  38  46  54]
     [102 110 118 126 134]
     [182 190 198 206 214]
     [262 270 278 286 294]
     [342 350 358 366 374]]
    """
    if not operation.lower() in ['sum', 'mean', 'average', 'avg']:
        raise ValueError("Operation {} not supported.".format(operation))
    if ndarray.ndim != len(new_shape):
        raise ValueError("Shape mismatch: {} -> {}".format(ndarray.shape,
                                                           new_shape))
    compression_pairs = [(d, c//d) for d, c in zip(new_shape,
                                                   ndarray.shape)]
    flattened = [l for p in compression_pairs for l in p]
    # print flattened
    ndarray = ndarray.reshape(flattened)
    for i in range(len(new_shape)):
        if operation.lower() == "sum":
            ndarray = ndarray.sum(-1*(i+1))
        elif operation.lower() in ["mean", "average", "avg"]:
            ndarray = ndarray.mean(-1*(i+1))
    return ndarray


# Convolution functions
def normgauss(sig, x=None, xc=0.):
    """Normalized Gaussian function.

    INPUT: sigma value (float), x (array), xc (float)

    OUTPUT: yvals (array)
    """
    if x is None:
        x = _np.linspace(-5 * sig, 5 * sig, 101)
    return 1. / (2*_np.pi*sig**2)**.5 * _np.exp(-(x - xc)**2. / (2 * sig**2.))


def normbox(hwidth, x=None, xc=0.):
    """Normalized Box function.

    INPUT: half-width (number), x (array), xc (float)

    OUTPUT: yvals (array)
    """
    if x is None:
        x = _np.linspace(-hwidth * 2 + xc, hwidth * 2 + xc, 101)
    y = _np.zeros(len(x))
    idx = _np.where(_np.abs(x - xc) <= hwidth)
    y[idx] = _np.zeros(len(y[idx])) + 2. / hwidth
    return y


def convnorm(x, arr, pattern):
    """Do the convolution of arr with pattern.
    Vector x is required for normalization. Its length must be odd!

    INPUT: units space (x, array), original array (arr), pattern (array)

    OUTPUT: array
    """
    if (len(x) != len(pattern)):
        _warn.warn('Wrong format of x and/or pattern arrays!')
        return None
    if (len(x) % 2 == 0):
        _warn.warn('Even length of arrays. Interpolating n-1 '
            'dimension!')
        idx = _np.argsort(x)
        x0 = x[idx]
        pattern0 = pattern[idx]
        x = _np.linspace(x0[0], x0[-1], len(x0)-1)
        pattern = _np.interp(x, x0, pattern0)
    dx = (x[-1] - x[0]) / (len(x) - 1.)
    cut = len(x) / 2
    return _np.convolve(pattern, arr)[cut:-cut] * dx


# 3D Coordinates manipulation (rotation)
def cart2sph(x, y, z=None):
    """ Cartesian to spherical coordinates.

    INPUT: arrays of same length

    OUTPUT: arrays """
    if z is None:
        z = _np.zeros(len(x))
    hxy = _np.hypot(x, y)
    r = _np.hypot(hxy, z)
    el = _np.arctan2(z, hxy)
    az = _np.arctan2(y, x)
    return r, az, el


def sph2cart(r, az, el=None):
    if el is None:
        el = _np.zeros(len(r))
    rcos_theta = r * _np.cos(el)
    x = rcos_theta * _np.cos(az)
    y = rcos_theta * _np.sin(az)
    z = r * _np.sin(el)
    return x, y, z


def cart_rot(x, y, z, ang_xy=0., ang_yz=0., ang_zx=0.):
    """ Apply rotation in Cartesian coordinates.

    INPUT: 3 arrays of same length, 3 angles (float, in radians).

    OUTPUT: arrays """
    rotmtx = _np.array([ [_np.cos(ang_zx) * _np.cos(ang_xy), -_np.cos(ang_yz) *
        _np.sin(ang_xy) + _np.sin(ang_yz) * _np.sin(ang_zx) * _np.cos(ang_xy),
        _np.sin(ang_yz) * _np.sin(ang_xy) + _np.cos(ang_yz) * _np.sin(ang_zx) *
        _np.cos(ang_xy)],
        [_np.cos(ang_zx) * _np.sin(ang_xy), _np.cos(ang_yz) * _np.cos(ang_xy) +
        _np.sin(ang_yz) * _np.sin(ang_zx) * _np.sin(ang_xy), -_np.sin(ang_yz) *
        _np.cos(ang_xy) + _np.cos(ang_yz) + _np.sin(ang_zx) * _np.sin(ang_xy)],
        [-_np.sin(ang_zx), _np.sin(ang_yz) * _np.cos(ang_zx), _np.cos(ang_yz) *
        _np.cos(ang_zx)] ])
    vec = _np.array([x, y, z])
    return _np.dot(rotmtx, vec)


# Lists and strings manipulation
def readpck(n, tp, ixdr, f):
    """ Read XDR 

    - n: length
    - tp: type ('i', 'l', 'f', 'd')
    - ixdr: counter
    - f: file-object

    :returns: ixdr (counter), _np.array
    """    
    sz = dict(zip(['i', 'l', 'f', 'd'], [4, 4, 4, 8]))
    s = sz[tp]
    upck = '>{0}{1}'.format(n, tp)
    return ixdr+n*s, _np.array(_struct.unpack(upck, f[ixdr:ixdr+n*s]))


def reshapeltx(ltxtb, ncols=2, latexfmt=True):
    """ Reshape a latex table. 

    :param ltxtb: latex table (only contents)
    :type ltxtb: string, or list of strings
    :rtype: string, or _np.array (str)
    :returns: latex formatted table, or matrix
    """
    t = ltxtb
    if not isinstance(t, _strtypes):
        t = ''.join(flatten(t))
    t = t.replace(' ', '').replace(r'\\', '&').replace('\n', '')
    t = t.split('&')
    t = _np.array(t)
    t = _np.delete(t, _np.where(t == ''))
    nadd = ncols - len(t) % ncols
    if nadd > 0:
        t = _np.append(t, _np.tile('', nadd))
    if latexfmt:
        return _tab(t.reshape((-1, ncols)), tablefmt='latex')
    else:
        return t.reshape((-1, ncols))


# def find_after(id, ref):
#     """ Returns the line where `id` if found after `ref`.
#     """
#     return 


def log_norm(x, vmax=1., vmin=0.01, autoscale=True, clip=True):
    """ Renormalize ``x`` (value or vector) making a correspondance of [0-1] 
    to [vmin, vmax] in log scale.

    If ``autoscale`` and x is a vector, ``vmax`` and ``vmin`` are automatically 
    set (``x`` must have at least a positive value).

    ``clip`` force the range to be between 0 and 1.
    """
    x = _np.array(x)
    if autoscale:
        (vmin, vmax) = (_np.min(x[_np.where(x > 0)]), _np.max(x))
    if vmin <= 0:
        _warn.warn('vmin <= 0 at phc.log_norm!')
        return x
    if vmax <= vmin:
        _warn.warn('vmax <= vmin at phc.log_norm!')
        return x
    if clip:
        x = x.clip(vmin, vmax)
    return _np.log10(x/vmin)/_np.log10(vmax/vmin)


def renormvals(xlist, xlims, ylims):
    """ Renormalize ``xlist`` according to ``_lims`` values.

    len(xlims)=len(ylims)=2
    """
    a = _np.diff(ylims)/_np.diff(xlims)
    b = ylims[0] - a*xlims[0]
    return a*_np.array(xlist)+b


def keys_values(keys, text, delimiter='_'):
    """ Return the values in a string coded as *KeyValueDelimiter*. The keys 
    do not need to be in order in the text. 

    Important! The text can not start with a keyword (if it is the case, add 
    a delimiter first) and it must have a delimiter after the last key. 

    TODO: Relax the important message considering starting with the first key 
    and a sign as last delimiter.

    Example: 

    .. code:: python

        keys = ['M', 'ob', 'H', 'Z', 'b']
        text = 'Be_M04.20_ob1.30_H0.30_Z0.014_bE_Ell.txt'

        print( phc.keys_values(keys, text) )
    """
    d = delimiter
    vals = []
    for k in keys:
        afterk = text[text.find(d+k)+len(d+k):]
        vals.append( afterk[:afterk.find(d)] )
    if len(vals) == 1:
        vals = vals[0]
    return vals


def fltTxtOccur(s, lines, n=1, seq=1, after=True, asstr=False):
    """ Return the seq-th float of the line after the n-th
    occurrence of `s` in the array `lines`.

    INPUT: s=string, lines=array of strings, n/seq=int (starting at 1)

    OUTPUT: float"""
    if isinstance(lines, _strtypes):
        lines = [lines]
    fltregex = r'[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?'
    if after:
        occur = [x[x.find(s) + len(s):] for x in lines if x.find(s) > -1]
    else:
        occur = [x for x in lines if x.find(s) > -1]
    out = _np.NaN
    if len(occur) >= n:
        occur = occur[n - 1]
        out = _re.findall(fltregex, occur)[seq - 1]
    if not asstr:
        out = float(out)
    return out


def strrep(seq, n, newseq):
    """ Insert `newseq` at position `n` of the string `seq`.

    seq[n] = seq[:n]+newseq+seq[n+1:]

    Note: the string at the position `n` is replaced!!!

    OUTPUT: string
    """
    return seq[:n] + newseq + seq[n + 1:]


def find_nearest(array, value, bigger=None, idx=False):
    """ Find nearest VALUE in the array and return it. 

    INPUT: array, value

    OUTPUT: closest value (array dtype)
    """
    if bigger is None:
        array = _np.array(array)
        i = (_np.abs(array - value)).argmin()
        found = array[i]
    elif bigger:
        found = _np.min([x for x in array if x > value])
        i = _np.where(array == found)
    elif not bigger:
        found = _np.max([x for x in array if x < value])
        i = _np.where(array == found)
    else:
        _warn.warn('# ERROR at bigger!!')
    # return
    if not idx:
        return found
    else:
        return i


def find_nearND(matrix, array, idx=False, bigger=None, outlen=1):
    """ Find nearest array values in a MATRIX.

    INPUT: array, value

    OUTPUT: closest value (array dtype)
    """
    if len(array) == 1:
        array = array[0]
    if _np.shape(matrix)[1] != len(array):
        raise ValueError('shape(matrix)[1] must be == len(array)') 
    dists = []
    for i in range(len(array)):
        nfact = _np.max(matrix[:, i]) - _np.min(matrix[:, i])
        if nfact <= 0:
            continue
        dists.append(_np.abs(matrix[:, i] - array[i])/nfact)

    dists = _np.sum(dists, axis=0)
    idsort = _np.argsort(dists)
    if bigger:
        valid = []
        for nid in idsort:
            chk = [True for i in range(len(array)) if (
                matrix[nid, i] > array[i])]
            if len(chk) == len(array):
                valid.extend([nid])
            if len(valid) == outlen:
                idsort = valid
                break
        if valid == []:
            _warn.warn('# Invalid values for bigger == {0}'.format(bigger))
    elif bigger is False:
        valid = []
        for nid in idsort:
            chk = [True for i in range(len(array)) if (
                matrix[nid, i] < array[i])]
            if len(chk) == len(array):
                valid.extend([nid])
            if len(valid) == outlen:
                idsort = valid
                break
        if valid == []:
            _warn.warn('# Invalid values for bigger == {0}'.format(bigger))
    # return
    if 'valid' in locals():
        if len(valid) == 0:
            print(idsort)
            raise ValueError('Invalid values for bigger == {0}'.format(bigger))
    if not idx:
        if outlen <= 1:
            return matrix[idsort[0]]
        else:
            return matrix[idsort[:outlen]]
    else:
        if outlen <= 1:
            return idsort[0]
        else:
            return idsort[:outlen]


def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x = nan_helper(y)
        >>> # x is a lambda "sequence" function to interpolation purposes
        >>> y[nans]= _np.interp(x(nans), x(~nans), y[~nans])
    """
    return _np.isnan(y), lambda z: z.nonzero()[0]


# Angular coordinates and dates manipulation
def dtflag(ms=False):
    """ Return a "datetime" flag, i.e., a string the the current date and time
    formated as yyyymmdd-hhMM."""
    now = _dt.datetime.now()
    if not ms:
        return '{0}{1:02d}{2:02d}-{3:02d}{4:02d}{5:02d}'.format(now.year, 
            now.month, now.day, now.hour, now.minute, now.second)
    else:
        return '{0}{1:02d}{2:02d}-{3:02d}{4:02d}{5:02d}{6:02.0f}'.format(
            now.year, now.month, now.day, now.hour, now.minute, now.second, 
            now.microsecond/1e4)


def longdate2MJD(ldate):
    """ FROM YYYY-MM-HHThh:mm:ss.sss to MJD (float). """
    ldate, hms = ldate.split('T')
    ldate = _np.array(ldate.split('-'), dtype=int)
    mjd = _jdcal.gcal2jd(*ldate)[1]
    return mjd + hms2fracday(hms)


def hms2fracday(hms):
    """ Enter hour:min:sec (string) and return fraction of a day (float) """
    hms = _np.array(hms.split(':'), dtype='float')
    return (hms[0] + hms[1]/60. + hms[2]/3600) / 24.


def fracday2hms(frac):
    """Enter fraction of a day (e.g., MJD) and return integers of hour, min,
    sec.

    INPUT: float

    OUTPUT: hour, min, sec (int, int, int)
    """
    hh = frac * 24
    if int(hh) > 0:
        mm = hh % int(hh)
        hh = int(hh)
    else:
        mm = hh
        hh = 0
    mm = mm * 60
    if int(mm) > 0:
        ss = mm % int(mm)
        mm = int(mm)
    else:
        ss = mm
        mm = 0
    ss = int(round(ss * 60))
    return hh, mm, ss


def ra2degf(rastr):
    """ RA to degrees (decimal). Input is string. """
    rastr = rastr.replace('::', ':')
    rastr = rastr.replace(',', '.')
    vals = _np.array(rastr.split(':')).astype(float)
    return (vals[0] + vals[1] / 60. + vals[2] / 3600.) * 360. / 24


def dec2degf(decstr, delimiter=":"):
    """ Sexagesimal to decimal. Input is string. """
    vals = _np.array(decstr.split(delimiter)).astype(float)
    if vals[0] < 0:
        vals[1:] *= -1
    return vals[0] + vals[1] / 60. + vals[2] / 3600.


def gentkdates(mjd0, mjd1, fact, step, dtstart=None):
    """ Generates round dates between ``mjd0`` and ``mjd1`` in a given step.
    Valid ``steps`` are:

        'd/D/dd/DD' for days;
        'w/W/ww/WW' for weeks;
        'm/M/mm/MM' for months;
        'y/Y/yy/YY/yyyy/YYYY' for years.

    dtstart (optional) is expected to be in _datetime._datetime.date() format
    [i.e., datetime.date(yyyy, m, d)].

    ``fact`` must be an integer.

    INPUT: float, float, float, int, step (see above), dtstart (see above)

    OUTPUT: list of _datetime.date
    """
    # check sanity of dtstart
    if dtstart is None:
        dtstart = _dt.datetime(*_jdcal.jd2gcal(_jdcal.MJD_0, mjd0)[:3]).date()
        mjdst = _jdcal.gcal2jd(dtstart.year, dtstart.month, dtstart.day)[1]
    else:
        mjdst = _jdcal.gcal2jd(dtstart.year, dtstart.month, dtstart.day)[1]
        if mjdst < mjd0 - 1 or mjdst > mjd1:
            _warn.warn('Invalid "dtstart". Using mjd0.')
            dtstart = _dt.datetime(
                *_jdcal.jd2gcal(_jdcal.MJD_0, mjd0)[:3]).date()
    # define step 'position' and vector:
    # if dtstart.day > 28:
    #     dtstart = dtstart.replace(day=28)
    basedata = dtstart
    dates = []
    mjd = mjdst
    while mjd < mjd1 + 1:
        dates += [basedata]
        if step.upper() in ['Y', 'YY', 'YYYY']:
            basedata += _dtdelta(years=fact)
        elif step.upper() in ['M', 'MM']:
            basedata += _dtdelta(months=fact)
        elif step.upper() in ['W', 'WW']:
            basedata += _dtdelta(weeks=fact)
        elif step.upper() in ['D', 'DD']:
            basedata += _dtdelta(days=fact)
        else:
            _warn.warn('# ERROR! Invalid step')
            raise SystemExit(1)
        mjd = _jdcal.gcal2jd(basedata.year, basedata.month, basedata.day)[1]
    return dates


def deg2rad(deg=1.):
    r""" :math:`1^\circ=\frac{\pi}{180}` rad."""
    return deg*_np.pi/180.


def rad2deg(rad=1.):
    r""" :math:`1^\circ=\frac{\pi}{180}` rad."""
    return rad*180./_np.pi


def arcmin2rad(arm=1.):
    r""" :math:`1'=\frac{\pi}{10800}` rad."""
    return arm*_np.pi/10800.


def rad2arcmin(rad=1.):
    r""" :math:`1'=\frac{\pi}{10800}` rad."""
    return rad*10800./_np.pi


def arcsec2rad(ars=1.):
    r""" :math:`1"=\frac{\pi}{648000}` rad."""
    return ars*_np.pi/648000.


def rad2arcsec(rad=1.):
    r""" :math:`1"=\frac{\pi}{648000}` rad. """
    return rad*648000./_np.pi


def mas2rad(mas=1.):
    r""" 1 mas :math:`=\frac{\pi}{648000000}` rad."""
    return mas*_np.pi/648000000.


def rad2mas(rad=1.):
    r""" 1 mas :math:`=\frac{\pi}{648000000}` rad."""
    return rad*648000000./_np.pi


def deg2mas(deg=1.):
    r""" :math:`1^\circ = 3600000` mas."""
    return deg*3600000.


def mas2deg(mas=1.):
    r""" :math:`1^\circ = 3600000` mas."""
    return mas/3600000.


# Files manipulation
def fileread(fname, bin=False):
    rmode = 'r'
    if bin:
        rmode = 'rb'
    ext = _os.path.splitext(fname)[1]
    if ext == '.gz':
        return _gzip.open(fname, rmode).read()
    elif ext == '.bz2':
        return _BZ2File(fname, rmode).read()
    else:
        return open(fname, rmode).read()


def readfixwd(fname, wlims, chklims=True):
    """ ``chklims`` : read line from beginning until the end
    """
    lines = fileread(fname).split('\n')
    lsz = [len(l) for l in lines]
    if wlims[-1] is None or wlims[-1] < 0:
        chklims = True
    wlims = [int(i-1) for i in wlims if (i is not None and i >= 0)]
    if wlims[0] < 0:
        wlims[0] = 0
    if chklims and wlims[0] != 0:
        wlims = [0] + wlims
    if chklims and wlims[-1] < max(lsz)-1:
        wlims += [max(lsz)]
    out = _np.chararray( (len(lines), len(wlims)-1), 
        itemsize=_np.max([ _np.max(_np.diff(wlims)), wlims[0] ]) )
    for j in range(len(lines)):
        out[j] = [ lines[j][wlims[i-1]:wlims[i]].strip().replace(',', '^') 
            for i in range(1, len(wlims))]
    return out


def readtextable(fname, commentchars='#', splitchars=r'&|\pm', 
    delchars=r'\$'+'\r '):
    lines = fileread(fname).split('\n')
    out = []
    for il in lines:
        if len(il) == 0:
            continue
        if il[0] in commentchars:
            continue
        for c in delchars:
            il = il.replace(c, '')
        out.append( _re.split('|'.join(splitchars.split('|')), il) )
    return out


def sortfile(file, quiet=False):
    """ Sort the file. """
    f0 = open(file, 'r')
    lines = f0.readlines()
    f0.close()
    lines.sort()
    f0 = open(file, 'w')
    f0.writelines(lines)
    f0.close()
    if not quiet:
        print('# File {0} sorted!'.format(file))
    return


def outfld(fold='hdt'):
    """
    Check and create (if necessary) an (sub)folder - generally used for output.

    INPUT: *fold=string

    OUTPUT: *system [folder creation]
    """
    _warn.warn('# Deprecated! Use `os` tools directly')
    if not _os.path.exists(fold):
        _os.system('mkdir {0}'.format(fold))
    return


def trimpathname(file):
    """Trim the full path string to return path and filename.

    INPUT: full file path

    OUTPUT: folder path, filename (strings)"""
    return [ file[:file.rfind('/') + 1], file[file.rfind('/') + 1:] ]


def rmext(name):
    """Remove the extension of a filename.
    Criteria: last `.` sets the extension.

    INPUT: filename (string)

    OUTPUT: filename without extension (string) """
    i = name.rfind('.')
    if i == -1:
        return name
    return name[:name.rfind('.')]


def readrange(file, i0, ie):
    """ Read a specific range of lines of a file with minimal memory use.

    Note that i == n-1 for the n-th line.

    INPUT: string, int, int

    OUTPUT: list of strings """
    _warn.warn('Update this with linecache')
    lines = []
    fp = open(file)
    for i, line in enumerate(fp):
        if i >= i0:
            lines += [line]
        if i >= ie:
            break
    fp.close()
    return lines


def splitfilelines(n, file):
    """ Break the *file* into *n* files. It also erases the expression "qsub ".

    OUTPUT: `file_##.txt` """
    f0 = open(file)
    lines = f0.readlines()
    f0.close()
    lines.sort()
    lines = [line.replace('qsub ', '') for line in lines]
    outname = trimpathname(file)[1].replace('.sh', '')
    N = len(lines)
    for i in range(n):
        f0 = open('{0}_{1:02d}.txt'.format(outname, i), 'w')
        f0.writelines(lines[i * N / n:(i + 1) * N / n])
        f0.close()
    print('# {0} files created!'.format(n))
    return


def recsearch(root='./', fstr=[''], allfstr=True, fullpath=False):
    """
    Do a recursive search in `root` looking for files containg `fstr`.

    `fstr` must be a list!

    If `allfstr` is True, *all* itens in fstr must be in the path+file names 
    structure. If False, *any* of the itens must be.

    INPUT: root (string), fstr (list of strings)

    OUTPUT: list of strings
    """
    outflist = []
    for i in range(len(fstr)):
        fstr[i] = fstr[i].replace('*', '')
    for root, subfolders, files in _os.walk(root):
        for f in files:
            fcomp = f
            if fullpath:
                fcomp = _os.path.join(root, f)
            if allfstr:
                if all(x in fcomp for x in fstr):
                    outflist.append(root+'/'+f)
            else:
                if any(x in fcomp for x in fstr):
                    outflist.append(root+'/'+f)
    return outflist


def renlist(root, newr):
    """ The routine changes each A_STR_B to A_NEW_B inside the running folder.
    """
    files = _glob('*{0}*'.format(root))
    files.sort()
    for i in range(len(files)):
        _os.system(
            'mv "' + files[i] + '" "' + files[i].replace(root, newr) + '"' )
        print("# " + files[i] + " renamed to: " +
              files[i].replace(root, newr) )
    return


def repl_fline_val(f, iline, oldval, newval):
    """ Replace `oldval` by `newval` in `f[iline]`

    return `f` replaced.
    """
    f[iline] = f[iline].replace(str(oldval), str(newval))
    return f


# Python-version stuff
def user_input(arg):
    if _sys.version_info[0] > 2:
        return input(arg)
    else: 
        return raw_input(arg)


# Plot-related
def civil_ticks(ax, civcfg=[1, 'm'], civdt=None, tklab=True, label="%y %b %d"):
    """ Add the civil ticks in the axis.

    :param civcfg: forces a given timestep between the ticks [`n`, interval].
    Interval can be day (`d`), month (`m`) or year (`y`). 
    :param civdt: sets the initial tick date. Format: [Y, M, D]
    :param tklab: if False, the civil date labels are not written.
    :param label: define the label for the ``date.strftime`` to be displayed in 
    the tick labels.

    :EXAMPLE:

        fig, ax = plt.subplots()
        ax.plot(x, y)
        ax = phc.civil_ticks(ax)
    """
    if civdt is not None:
        civdt = _dt.datetime(civdt[0], civdt[1], civdt[2]).date()
    mjd0, mjd1 = ax.get_xlim()
    dtticks = gentkdates(mjd0, mjd1, civcfg[0], civcfg[1], dtstart=civdt)
    mjdticks = [_jdcal.gcal2jd(date.year, date.month, date.day)[1] for 
        date in dtticks]
    ax2 = ax.twiny()
    ax2.set_xlim( ax.get_xlim() )
    ax2.set_xticks(mjdticks)
    if tklab:
        ax2.set_xticklabels([date.strftime(label) for date in dtticks])
    else:
        ax2.set_xticklabels([])
    return ax


def savefig(fig, figname=None, fmt=['png'], keeppt=False, dpi=80, transp=True):
    """ Standard way of saving a figure in PyHdust. """
    if figname is None or figname == "":
        figname = dtflag(ms=True)
    elif not keeppt:
        figname = figname.replace('.', 'p')
    if _os.path.basename(figname) == figname:
        figname = _os.getcwd() + '/' + figname
    for f in fmt:
        fig.savefig(figname+'.{0}'.format(f), transparent=transp, 
            bbox_inches='tight', dpi=dpi)
        print('# Saved {1}.{0}'.format(f, figname))
    _plt.close(fig)
    return


def normGScale(val, min=None, max=None, log=False):
    """ Return the normalized values of a given array to 0 and 255 (gray 
    scale).

    If `log` then the normalization is done in this scale.

    If `min` and `max` are not set, it is assumed that the values are from the 
    list.

    .. code::

        >>> phc.normGScale(np.linspace(0,10,10), 0, 10)
        array([  0,  28,  57,  85, 113, 142, 170, 198, 227, 255])
        >>> phc.normGScale(np.linspace(0,10,10))
        array([  0,  28,  57,  85, 113, 142, 170, 198, 227, 255])
        >>> phc.normGScale(np.linspace(0,10,10), 0, 10, log=True)
        array([  0,   0,  80, 128, 161, 187, 208, 226, 241, 255])
        >>> phc.normGScale(np.logspace(0,1,10), 0, 10, log=True)
        array([  0,  28,  57,  85, 113, 142, 170, 198, 227, 255])
        >>> phc.normGScale(np.logspace(0,1,10), 0, 10)
        array([ 26,  33,  43,  55,  71,  92, 118, 153, 197, 255])
        >>> phc.normGScale(np.logspace(0,1,10))
        array([  0,   8,  19,  33,  51,  73, 103, 142, 191, 255])
    """
    if len(val) == 1 and (min is None or max is None):
        _warn.warn('Wrong normGScale call!!')
        return 127
    #
    val = _np.array(val).astype(float)
    if min is None:
        min = _np.min(val)
    if max is None:
        max = _np.max(val)
    # 
    if not log:
        val = (val - min) / (max - min) * 255
    else:
        if min <= 0:
            min = _np.min(val[_np.where(val > 0)])
            val[_np.where(val <= 0)] = min
        val = (_np.log(_np.array(val).astype(float)) - _np.log(min)) / \
            (_np.log(max) - _np.log(min)) * 255
    return _np.round(val).astype(int)


def gradColor(val, cmapn='jet', min=None, max=None, log=False):
    """ Return the corresponding value(s) color of a given colormap.

    Good options, specially for lines, are 'jet', 'gnuplot', 'brg', 
    'cool' and 'gist_heat' (attention! Here max is white!). 

    .. code-block:: python

        cor = phc.gradColor(arange(10), cmapn='gist_heat')
        for i in range(0,10):
            cor = phc.gradColor(arange(10), cmapn='gist_heat')[i]
            print cor
            plt.plot(arange(5)+i, color=cor, label='GB='+('{0:4.2f},'*3).
                format(*cor)[:-1])

        plt.legend(fontsize=8)

    .. image:: _static/phc_gradColor.png
        :width: 512px
        :align: center
        :alt: phc.gradColor example
    """
    val = normGScale(val, min=min, max=max, log=log)
    cmap = _plt.get_cmap(cmapn)
    return cmap(val)


def cycles(i=0, ctype='cor'):
    """ Cycle between values of the phc.colors, phc.line_styles and 
    phc.filled_markers lists. 

    INPUT: valid ctypes are: ['cor','ls','mk']

    OUTPUT: the corresponding value of the list. """
    if ctype not in ['cor', 'ls', 'mk']:
        _warn.warn('Invalid ctype calling phc.cycles!')
        return
    elif ctype == 'cor':
        return colors[_np.mod(i, len(colors))]
    elif ctype == 'ls':
        return line_styles[_np.mod(i, len(line_styles))]
    elif ctype == 'mk':
        return filled_markers[_np.mod(i, len(filled_markers))]


def dashes(i=0):
    """ Dashes scheme for plot 
    """
    i = int(_np.mod(i, 7))
    if i == 0:
        return []
    if i < 4:
        return (32/2**i, 4, 4, 4)
    if i < 7:
        return (32/2**(i-2), 4, 2, 4)


colors = ['Black', 'Blue', 'Green', 'red', 'orange', 'brown', 'purple', 'gray',
    'dodgerblue', 'lightgreen', 'tomato', 'yellow', 'peru', 'MediumVioletRed',
    'LightSteelBlue', 'cyan', 'darkred', 'olive']

filled_markers = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 
    'd']

line_styles = ['-', '--', '-.', ':']


# Physical functions
def BBlbd(T, lbd=None):
    """ Black body radiation as function of lambda. CGS units (erg s-1 sr−1 cm−3).

    INPUT: lambda vector in cm. If None, lbd = _np.arange(1000, 10000, 100) *
    1e-8 #Angs -> cm"""
    if lbd is None:
        lbd = _np.arange(1000, 10000., 100) * 1e-8  # Angs -> cm
    ft = h.cgs * c.cgs / (lbd * kB.cgs * T)
    return 2 * h.cgs * c.cgs**2 / (lbd**5 * (_np.exp(ft) - 1))


def fBBcor(T):
    """ Correction as appendix B Vieira+2015. Stellar atmospheric models
    systematically have a LOWER flux than a BB of a given Teff temperature in 
    IR (at least for Early-type stars).

    fBBcor(T)*BBlbd(Teff) == Kurucz(Teff) 

    INPUT: Teff (Kurucz). log g = 4.0

    OUTPUT: fBBcor(T) < 1.0 """
    return 1.015 - 0.301 * (T / 1e4) + 0.064 * (T / 1e4)**2.


def gbf(T, lbd):
    r""" Gaunt factors from Vieira+2015. 

    INPUT: T (K) and lbd (:math:`\mu`m, array)

    log(T /K) G0 G1 G2 B0 B1 B2
    """
    vals = _np.array([
        3.70, 0.0952, 0.0215, 0.0145, 2.2125, -1.5290, 0.0563,
        3.82, 0.1001, 0.0421, 0.0130, 1.6304, -1.3884, 0.0413,
        3.94, 0.1097, 0.0639, 0.0111, 1.1316, -1.2866, 0.0305,
        4.06, 0.1250, 0.0858, 0.0090, 0.6927, -1.2128, 0.0226,
        4.18, 0.1470, 0.1071, 0.0068, 0.2964, -1.1585, 0.0169,
        4.30, 0.1761, 0.1269, 0.0046, -0.0690, -1.1185, 0.0126,
    ]).reshape((6, -1))
    if T < 5000 or T > 22500:
        _warn.warn('# ERROR! Invalid temperature for Gaunt factors '
            'calculation!')
        return _np.zeros(len(lbd)), _np.zeros(len(lbd))
    elif T >= 5000 and T < 10**vals[0, 0]:
        _warn.warn('Extrapolated Gaunt factors!!')
        g0, g1, g2, b0, b1, b2 = vals[0, 1:]
    elif T <= 22500 and T > 10**vals[-1, 0]:
        _warn.warn('Extrapolated Gaunt factors!!')
        g0, g1, g2, b0, b1, b2 = vals[-1, 1:]
    else:
        i = _np.where(vals[:, 0] == find_nearest(
            vals[:, 0], _np.log10(T), bigger=False))[0]
        # print i, vals[i,0]
        g0 = interLinND(
            [_np.log10(T)], [vals[i, 0]], [vals[i + 1, 0]], vals[i:i + 2, 1])
        g1 = interLinND(
            [_np.log10(T)], [vals[i, 0]], [vals[i + 1, 0]], vals[i:i + 2, 2])
        g2 = interLinND(
            [_np.log10(T)], [vals[i, 0]], [vals[i + 1, 0]], vals[i:i + 2, 3])
        b0 = interLinND(
            [_np.log10(T)], [vals[i, 0]], [vals[i + 1, 0]], vals[i:i + 2, 4], 
            disablelog=True)
        b1 = interLinND(
            [_np.log10(T)], [vals[i, 0]], [vals[i + 1, 0]], vals[i:i + 2, 5], 
            disablelog=True) 
        b2 = interLinND(
            [_np.log10(T)], [vals[i, 0]], [vals[i + 1, 0]], vals[i:i + 2, 6]) 
    return _np.exp(g0 + g1 * _np.log(lbd) + g2 * _np.log(lbd)**2), _np.exp(b0 + 
        b1 * _np.log(lbd) + b2 * _np.log(lbd)**2)


def lawkep(M=None, m=None, P=None, a=None):
    """ Kepler law calc. Kepp `None` on what you what to calc.

    Units are in *Solar System* one, what is, masses in Msun and `P` in years.

    `a` is the distance between the two bodies, measured in AU. 

    One can also use the Centre-of-Mass equation a1*M1 = a2*M2 to relate the 
    two masses. 
    """
    Gc = G.cgs
    pi2_4 = 4*_np.pi**2
    if M is None:
        m *= Msun.cgs
        P *= yr.cgs
        a *= au.cgs
        return (pi2_4*a**3/P**2/Gc-m)/Msun.cgs
    elif m is None:
        M *= Msun.cgs
        P *= yr.cgs
        a *= au.cgs
        return (pi2_4*a**3/P**2/Gc-M)/Msun.cgs
    elif P is None:
        M *= Msun.cgs
        m *= Msun.cgs
        a *= au.cgs
        return _np.sqrt( pi2_4*a**3/(M+m)/Gc )/yr.cgs
    elif a is None:
        M *= Msun.cgs
        m *= Msun.cgs
        P *= yr.cgs
        return ( P**2*Gc*(M+m)/pi2_4 )**(1./3)/au.cgs
    else:
        print('# Wrong call of phc.lawkep! Keep `None` to calc that qtt.')
        return None


class Constant(object):

    """ Class for a physical/astronomical constant
    """

    def __init__(self, cgs, SI, unitscgs='', info='No available description'):
        self.cgs = cgs
        self.SI = SI
        self.unitscgs = unitscgs
        self.info = info

    def __repr__(self):
        return str('{0:.7e} in {1} (cgs)'.format(self.cgs, self.unitscgs))


# From CODATA/NIST in May/2016 (related to Mohr+2015)
#  http://arxiv.org/abs/1507.07956
nA = Constant(6.022140857e23, 6.022140857e23, '', "Avogadro's constant")
amu = Constant(1.66053904e-24, 1.66053904e-27, 'g', 'Unified atomic mass unit')
me = Constant(9.10938356e-28, 9.10938356e-31, 'g', 'Mass of electron')
mp = Constant(1.672621898e-24, 1.672621898e-27, 'g', 'Mass of proton')
mn = Constant(1.00137841898*mp.cgs, 1.00137841898*mp.SI, 'g', 
    'Mass of neutron')
h = Constant(6.62607004e-27, 6.62607004e-34, 'erg s-1', 'Planck constant')
eV = Constant(1.6021766208e-12, 1.6021766208e-19, 'erg', 'Electron volt')
kB = Constant(1.38064852e-16, 1.38064852e-23, 'erg K-1', 'Boltzmann constant')
alpha = Constant(7.2973525664e-3, 7.2973525664e-3, '', 
    'Fine structure constant')
Rinf = Constant(109737.31568, 10973731.568, 'cm-1', 'Rydberg constant')
sigT = Constant(6.6524587158e-25, 6.6524587158e-29, 'cm2', 
    'Thomson cross section')
# From Luzum et al., 2011
au = Constant(1.49597870700e13, 1.49597870700e11, 'cm', 'Astronomical unit')
# From Prsa & Harmanec, 2012
G = Constant(6.67384e-8, 6.67384e-11, 'cm3 g-1 s-2', 
    'Gravitational constant')
Mea = Constant(398600.4418e15/G.cgs, 398600.4418e9/G.SI, 'g', 'Earth mass')
Rea = Constant(6371.e5, 6371.e3, 'cm', 'Earth radius')
Mju = Constant(126686535e15/G.cgs, 126686535e9/G.SI, 'g', 'Jupiter mass')
Rju = Constant(71492.e5, 71492.e3, 'cm', 'Jupiter radius')
c = Constant(2.99792458e10, 299792458., 'cm s-1', 'speed of light in vacuum')
sigma = Constant(5.670400e-5, 5.670400e-8, 'erg cm-2 K-4 s-1', 
    'Stefan-Boltzmann constant')
Msun = Constant(1.988547e33, 1.988547e30, 'g', 'Solar mass')
Rsun = Constant(6.95508e10, 695508e3, 'cm', 'Solar radius')
Lsun = Constant(3.846e33, 3.846e26, 'erg s-1', 'Solar luminosity')
Tsun = Constant(5779.57, 5779.57, 'K', 'Solar Temperature')
# Derived quantities
# In astronomy, the Julian year is defined as 365.25 days of exactly 86400 SI 
#  seconds each, totalling exactly 31557600 SI (IAU style manual, 1989, 
#  Wilkins)
# For the Gregorian calendar the average length of the calendar year (the mean 
#  year) across the complete leap cycle of 400 years is 365.2425 days!
yr = Constant(60*60*24*365.25, 60*60*24*365.25, 'sec', 'year')
ly = Constant(yr.cgs*c.cgs, yr.SI*c.SI, 'cm', 'Light year')
pc = Constant(au.cgs*60*60*180/_np.pi, au.SI*60*60*180/_np.pi, 'cm', 'Parsec')
e = Constant(10*c.cgs*1.6021766208e-19, 1.6021766208e-19, 'esu', 
    'Elementary charge')
a = Constant(4*sigma.cgs/c.cgs, 4*sigma.SI/c.SI, 'erg cm-3 K-4', 
    'Radiation density constant')
hbar = Constant(h.cgs/2/_np.pi, h.SI/2/_np.pi, 'erg s', 
    'Planck constant/(2*pi)')
# From Wikipedia
ep0 = Constant(1., 8.854187187e-12, '', 'Permittivity of Free Space')
mH = Constant(1.00794*amu.cgs, 1.00794*amu.SI, 'g', 'Mass of hydrogen')


# MAIN ###
if __name__ == "__main__":
    pass
