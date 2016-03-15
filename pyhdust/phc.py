# -*- coding:utf-8 -*-

"""
PyHdust *phc* module: physical constants and general use functions

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

.. _`List of constants`: phc_list.html


:license: GNU GPL v3.0 https://github.com/danmoser/pyhdust/blob/master/LICENSE
"""
import os as _os
import re as _re
import numpy as _np
import datetime as _dt
from glob import glob as _glob
from itertools import product as _product
import pyhdust.jdcal as _jdcal

try:
    import matplotlib.pyplot as _plt
    from scipy import optimize as _optimize
except:
    print('# Warning! matplotlib and/or scipy module not installed!!!')

__author__ = "Daniel Moser"
__email__ = "dmfaes@gmail.com"


# Data manipulation (includes fitting)
def wg_avg_and_std(values, sigma):
    """
    Return the weighted average and standard deviation.

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


def bindata(x, y, nbins=20, yerr=None, xrange=None):
    """
    Return the weighted binned data.

    if yerr is None:
        yerr = _np.ones(shape=_np.shape(x))

    INPUT: x, y, err - arrays with the same shape (they don't need to be
    sorted); nbins=int, xrange=[xmin, xmax]

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
    if xrange is None:
        xmax = _np.max(x)
        xmin = _np.min(x)
    else:
        xmin, xmax = xrange
    shift = (xmax - xmin) / (nbins - 1)
    tmpx = _np.arange(nbins) * shift + xmin
    tmpy = _np.zeros(nbins)
    tmpyerr = _np.zeros(nbins)
    idx = _np.zeros(nbins, dtype=bool)
    for i in range(nbins):
        selx = _np.where( abs(x - tmpx[i]) <= shift / 2. )
        if len(selx[0]) >= 1:
            tmpx[i] = _np.average(x[selx])
            tmpy[i], tmpyerr[i] = wg_avg_and_std(y[selx], yerr[selx])
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


def interBarND(point, grid, Zgrid):
    """ TODO """
    d = _np.sqrt( (point[0] - grid[:, 0])**2 + (point[1] - grid[:, 1])**2)
    print(d)
    return


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
            # return np.sum( ((y-func(p,x))/yerr)**2 )
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
    >>> m = np.arange(0,100,1).reshape((10,10))
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
    print flattened
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
    return 1. / (2 * _np.pi)**.5 / sig * _np.exp(-(x - xc)**2. / (2 * sig**2.))


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
    if (len(x) % 2 == 0) or (len(x) != len(pattern)):
        print('# Wrong format of x and/or pattern arrays!')
        return None
    dx = (x[-1] - x[0]) / (len(x) - 1.)
    cut = len(x) / 2
    return _np.convolve(pattern, arr)[cut:-cut] * dx


# 3D Coordinates manipulation (rotation)
def cart2sph(x, y, z):
    """ Cartesian to spherical coordinates.

    INPUT: arrays of same length

    OUTPUT: arrays """
    r = _np.sqrt(x**2 + y**2 + z**2)
    idx = _np.where(r == 0)
    r[idx] = 1e-9
    th = _np.arccos(z / r)
    phi = _np.arctan2(y, x)
    # ind = (y<0) & (x<0)
    # phi[ind] = phi+_np.pi
    return r, th, phi


def sph2cart(r, th, phi):
    """  Spherical to Cartesian coordinates.

    INPUT: arrays of same length

    OUTPUT: arrays """
    x = r * _np.sin(th) * _np.cos(phi)
    y = r * _np.sin(th) * _np.sin(phi)
    z = r * _np.cos(th)
    return x, y, z


def cart_rot(x, y, z, ang_xy, ang_yz, ang_zx):
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


def rotate_coords(x, y, theta, ox, oy):
    """Rotate arrays of coordinates x and y by theta radians about the
    point (ox, oy).

    This routine was inspired on a http://codereview.stackexchange.com post.
    """
    s, c = _np.sin(theta), _np.cos(theta)
    x, y = _np.asarray(x) - ox, _np.asarray(y) - oy
    return x * c - y * s + ox, x * s + y * c + oy


def rotate_image(src, theta, ox=None, oy=None, fill=0):
    """Rotate the image src by theta radians about (ox, oy).
    Pixels in the result that don't correspond to pixels in src are
    replaced by the value fill.

    This routine was inspired on a http://codereview.stackexchange.com post.
    """
    # Images have origin at the top left, so negate the angle.
    theta = -theta

    # Dimensions of source image. Note that scipy.misc.imread loads
    # images in row-major order, so src.shape gives (height, width).
    sh, sw = src.shape
    if ox is None:
        ox = sw/2
    if oy is None:
        oy = sh/2

    # Rotated positions of the corners of the source image.
    cx, cy = rotate_coords([0, sw, sw, 0], [0, 0, sh, sh], theta, ox, oy)

    # Determine dimensions of destination image.
    dw, dh = (int(_np.ceil(c.max() - c.min())) for c in (cx, cy))

    # Coordinates of pixels in destination image.
    dx, dy = _np.meshgrid(_np.arange(dw), _np.arange(dh))

    # Corresponding coordinates in source image. Since we are
    # transforming dest-to-src here, the rotation is negated.
    sx, sy = rotate_coords(dx + cx.min(), dy + cy.min(), -theta, ox, oy)

    # Select nearest neighbour.
    sx, sy = sx.round().astype(int), sy.round().astype(int)

    # Mask for valid coordinates.
    mask = (0 <= sx) & (sx < sw) & (0 <= sy) & (sy < sh)

    # Create destination image.
    dest = _np.empty(shape=(dh, dw), dtype=src.dtype)

    # Copy valid coordinates from source image.
    dest[dy[mask], dx[mask]] = src[sy[mask], sx[mask]]

    # Fill invalid coordinates.
    dest[dy[~mask], dx[~mask]] = fill

    return dest


# Lists and strings manipulation
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
        # 
    """
    d = delimiter
    vals = []
    for k in keys:
        afterk = text[text.find(d+k)+len(d+k):]
        vals.append( afterk[:afterk.find(d)] )
    return vals


def fltTxtOccur(s, lines, n=1, seq=1, after=False, asstr=False):
    """ Return the seq-th float of the line after the n-th
    occurrence of `s` in the array `lines`.

    INPUT: s=string, lines=array of strings, n/seq=int (starting at 1)

    OUTPUT: float"""
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
        idx = (_np.abs(array - value)).argmin()
        found = array[idx]
    elif bigger:
        found = _np.min([x for x in array if x > value])
        idx = _np.where(array == found)
    elif not bigger:
        found = _np.max([x for x in array if x < value])
        idx = _np.where(array == found)
    else:
        print('# ERROR at bigger!!')
    # return
    if not idx:
        return found
    else:
        return idx[0]


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
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """
    return _np.isnan(y), lambda z: z.nonzero()[0]


# Angular coordinates and dates manipulation
def dtflag():
    """ Return a "datetime" flag, i.e., a string the the current date and time
    formated as yyyymmdd-hhMM."""
    now = _dt.datetime.now()
    return '{0}{1:02d}{2:02d}-{3:02d}{4:02d}{5:02d}'.format(now.year, 
        now.month, now.day, now.hour, now.minute, now.second)


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
    """ Generates round dates between > mjd0 and < mjd1 in a given step.
    Valid steps are:

        'd/D/dd/DD' for days;
        'm/M/mm/MM' for months;
        'y/Y/yy/YY/yyyy/YYYY' for years.

    dtstart (optional) is expected to be in _datetime._datetime.date() format
    [i.e., _datetime.date(yyyy, m, d)].

    fact must be an integer.

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
            print('# Warning! Invalid "dtstart". Using mjd0.')
            dtstart = _dt.datetime(
                *_jdcal.jd2gcal(_jdcal.MJD_0, mjd0)[:3]).date()
    # define step 'position' and vector:
    if dtstart.day > 28:
        dtstart = dtstart.replace(day=28)
    basedata = [dtstart.year, dtstart.month, dtstart.day]
    dates = []
    mjd = mjdst
    if step.upper() in ['Y', 'YY', 'YYYY']:
        i = 0
        while mjd < mjd1 + 1:
            dates += [_dt.datetime(*basedata).date()]
            basedata[i] += fact
            mjd = _jdcal.gcal2jd(*basedata)[1]
    elif step.upper() in ['M', 'MM']:
        i = 1
        while mjd < mjd1 + 1:
            dates += [_dt.datetime(*basedata).date()]
            basedata[i] += fact
            while basedata[i] > 12:
                basedata[0] += 1
                basedata[1] -= 12
            mjd = _jdcal.gcal2jd(*basedata)[1]
    elif step.upper() in ['D', 'DD']:
        i = 2
        # daysvec = _np.arange(1, 29, fact)
        # if basedata[i] not in daysvec:
        #     j = 0
        #     while daysvec[j + 1] < basedata[i]:
        #         j += 1
        #     daysvec += basedata[i] - daysvec[j]
        #     idx = _np.where(daysvec < 29)
        #     daysvec = daysvec[idx]
        # else:
        #     j = _np.where(daysvec == basedata[i])[0]
        j = 0
        daysvec = _np.arange(1, 29)[j::fact]
        while basedata[i] not in daysvec:
            j += 1
            daysvec = _np.arange(1, 29)[j::fact]
        j = 0
        while mjd < mjd1 + 1:
            dates += [_dt.datetime(*basedata).date()]
            j += 1
            if j == len(daysvec):
                j = 0
                basedata[1] += 1
                if basedata[1] == 13:
                    basedata[1] = 1
                    basedata[0] += 1
            basedata[i] = daysvec[j]
            mjd = _jdcal.gcal2jd(*basedata)[1]
    else:
        print('# ERROR! Invalid step')
        raise SystemExit(1)
    return dates


def deg2rad(deg=1.):
    """ :math:`1^\circ=\frac{pi}{180}` """
    return deg*_np.pi/180.


def rad2deg(rad=1.):
    """ :math:`1^\circ=\frac{pi}{180}` """
    return rad*180./_np.pi


def arcmin2rad(arm=1.):
    """ :math:`1'=\frac{pi}{10800}` """
    return arm*_np.pi/10800.


def rad2arcmin(rad=1.):
    """ :math:`1^\circ=\frac{pi}{10800}` """
    return rad*10800./_np.pi


def arcsec2rad(ars=1.):
    """ :math:`1'=\frac{pi}{648000}` """
    return ars*_np.pi/648000.


def rad2arcsec(rad=1.):
    """ :math:`1^\circ=\frac{pi}{648000}` """
    return rad*648000./_np.pi


def mas2rad(mas=1.):
    """ :math:`1^\circ=\frac{pi}{648000000}` """
    return mas*_np.pi/648000000.


def rad2mas(rad=1.):
    """ :math:`1^\circ=\frac{pi}{648000000}` """
    return rad*648000000./_np.pi


def deg2mas(deg=1.):
    """ :math:`1^\circ = 3600000` """
    return deg*3600000.


def mas2deg(mas=1.):
    """ :math:`1^\circ = 3600000` """
    return mas/3600000.


# Files manipulation
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


def recsearch(root='./', fstr=[''], allfstr=True):
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
            if allfstr:
                if all(x in f for x in fstr):
                    outflist.append(root+'/'+f)
            else:
                if any(x in f for x in fstr):
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


# Plot-related
def civil_ticks(ax, civcfg=[1, 'm'], civdt=None, tklab=True):
    """ Add the civil ticks in the axis.  """
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
        ax2.set_xticklabels([date.strftime("%d %b %y") for date in dtticks])
    else:
        ax2.set_xticklabels([])
    return ax


def savefig(fig, figname=None, fmt=['png'], keeppt=False):
    """ Standard way of saving a figure in PyHdust. """
    if figname is None or figname == "":
        figname = dtflag()
    elif not keeppt:
        figname = figname.replace('.', 'p')
    if _os.path.basename(figname) == figname:
        figname = _os.getcwd() + '/' + figname
    for f in fmt:
        print('# Saved {1}.{0}'.format(f, figname))
        fig.savefig(figname+'.{0}'.format(f), transparent=True)
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
        print('# Warning! Wrong normGScale call!!')
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
        print('# Warning! Invalid ctype calling phc.cycles!')
        return
    elif ctype == 'cor':
        return colors[_np.mod(i, len(colors))]
    elif ctype == 'ls':
        return line_styles[_np.mod(i, len(line_styles))]
    elif ctype == 'mk':
        return filled_markers[_np.mod(i, len(filled_markers))]


colors = ['Black', 'Blue', 'Green', 'red', 'orange', 'brown', 'purple', 'gray',
    'dodgerblue', 'lightgreen', 'tomato', 'yellow', 'peru', 'MediumVioletRed',
    'LightSteelBlue', 'cyan', 'darkred', 'olive']

filled_markers = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 
    'd']

line_styles = ['-', '--', ':', '-.']


# Physical functions
def BBlbd(T, lbd=None):
    """ Black body radiation as function of lambda. CGS units (erg s-1 sr−1 cm−3).

    INPUT: lambda vector in cm. If None, lbd = _np.arange(1000, 10000, 100) *
    1e-8 #Angs -> cm"""
    if lbd is None:
        lbd = _np.arange(1000, 10000, 100) * 1e-8  # Angs -> cm
    ft = h.cgs * c.cgs / (lbd * kB.cgs * T)
    return 2 * h.cgs * c.cgs**2 / (lbd**5 * (_np.exp(ft) - 1))


def lawkep(M=None, m=None, P=None, a=None):
    """ Kepler law calc. But `None` on what you what to calc.

    Units are in *Solar System* one, what is, masses in Msun and `P` in years.

    `a` is the distance between the two bodies, measured in AU. 

    One can also use the Centre-of-Mass equation a1*M1 = a2*M2 to relate the 
    two masses. 
    """
    if M is None:
        m *= Msun.cgs
        P *= yr.cgs
        a *= au.cgs
        return 4*_np.pi**2*a**3/P**2/G.cgs-m
    elif m is None:
        M *= Msun.cgs
        P *= yr.cgs
        a *= au.cgs
        return 4*_np.pi**2*a**3/P**2/G.cgs-M
    elif P is None:
        M *= Msun.cgs
        m *= Msun.cgs
        a *= au.cgs
        return _np.sqrt( 4*_np.pi**2*a**3/(M+m)/G.cgs )
    elif a is None:
        M *= Msun.cgs
        m *= Msun.cgs
        P *= yr.cgs
        return ( P**2*G.cgs*(M+m)/4/_np.pi**2 )**(1./3)
    else:
        print('# Wrong call of phc.lawkep! Put `None` calc that qtt.')
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


# Constants
c = Constant(2.99792458e10, 299792458., 'cm s-1', 'speed of light in vacuum')
h = Constant(6.6260755e-27, 6.62606957e-34, 'erg s-1', 'Planck constant')
hbar = Constant(
    1.05457266e-27, 1.05457172534e-34, 'erg s', 'Planck constant/(2*pi)')
G = Constant(
    6.674253e-8, 6.674253e-11, 'cm3 g-1 s-2', 'Gravitational constant')
e = Constant(4.8032068e-10, 1.60217657e-19, 'esu', 'Elementary charge')
ep0 = Constant(1., 8.8542e-12, '', 'Permittivity of Free Space')
me = Constant(9.1093897e-28, 9.10938291e-31, 'g', 'Mass of electron')
mp = Constant(1.6726231e-24, 1.67262178e-27, 'g', 'Mass of proton')
mn = Constant(1.6749286e-24, 1.674927351e-27, 'g', 'Mass of neutron')
mH = Constant(1.6733e-24, 1.6737236e-27, 'g', 'Mass of hydrogen')
amu = Constant(1.6605402e-24, 1.66053892e-27, 'g', 'Atomic mass unit')
nA = Constant(6.0221367e23, 6.0221413e23, '', "Avagadro's number")
kB = Constant(1.3806504e-16, 1.3806488e-23, 'erg K-1', 'Boltzmann constant')
eV = Constant(1.6021764871e-12, 1.602176565e-19, 'erg', 'Electron volt')
a = Constant(
    7.5646e-15, 7.5646e-16, 'erg cm-3 K-4', 'Radiation density constant')
sigma = Constant(
    5.67051e-5, 5.670373e-8, 'erg cm-2 K-4 s-1', 'Stefan-Boltzmann constant')
alpha = Constant(7.29735308e-3, 7.2973525698e-3, '', 'Fine structure constant')
Rinf = Constant(109737.316, 10973731.6, 'cm', 'Rydberg constant')
sigT = Constant(
    6.65245854533e-25, 6.65245854533e-29, 'cm2', 'Thomson cross section')

au = Constant(1.49597870691e13, 1.49597870691e11, 'cm', 'Astronomical unit')
pc = Constant(3.08567758e18, 3.08567758e16, 'cm', 'Parsec')
ly = Constant(9.4605284e17, 9.4605284e15, 'cm', 'Light year')
Mea = Constant(5.9722e27, 5.9722e24, 'g', 'Earth mass')
Rea = Constant(6371.e5, 6371.e3, 'cm', 'Earth radius')
Msun = Constant(1.9891e33, 1.9891e30, 'g', 'Solar mass')
Rsun = Constant(6.961e10, 696100e3, 'cm', 'Solar radius')
Lsun = Constant(3.846e33, 3.846e26, 'erg s-1', 'Solar luminosity')
Tsun = Constant(5778., 5778., 'K', 'Solar Temperature')

yr = Constant(3.15569e7, 3.15569e7, 'sec', 'year')

bestars = [
    # The numbers below are based on Harmanec 1988
    # SpType     Tpole    Mass    Rp      Lum     Rp2
    ['B0',   29854, 14.57, 05.80, 27290, 6.19],
    ['B0.5', 28510, 13.19, 05.46, 19953, 5.80],
    ['B1',   26182, 11.03, 04.91, 11588, 5.24],
    ['B2',   23121, 08.62, 04.28,  5297, 4.55],
    ['B3',   19055, 06.07, 03.56,  1690, 3.78],
    ['B4',   17179, 05.12, 03.26,   946, 3.48],
    ['B5',   15488, 04.36, 03.01,   530, 3.21],
    ['B6',   14093, 03.80, 02.81,   316, 2.99],
    ['B7',   12942, 03.38, 02.65,   200, 2.82],
    ['B8',   11561, 02.91, 02.44,   109, 2.61],
    ['B9',   10351, 02.52, 02.25,   591, 2.39],
    ['B9.5',  9886, 02.38, 02.17,    46, 2.32]]


# MAIN ###
if __name__ == "__main__":
    pass
