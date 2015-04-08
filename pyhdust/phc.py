#-*- coding:utf-8 -*-

""" 
PyHdust *phc* module: physical constants and general use functions

`List of constants`_

.. _`List of constants`: phc_list.html

:license: GNU GPL v3.0 (https://github.com/danmoser/pyhdust/blob/master/LICENSE)
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
except:
    print('# Warning! matplotlib module not installed!!!')

__author__ = "Daniel Moser"
__email__ = "dmfaes@gmail.com"


def dtflag():
    """ Return a "datetime" flag, i.e., a string the the current date and time
    formated as yyyy-mm-dd_hh-MM."""
    now = _dt.datetime.now()
    return '{0}-{1:02d}-{2:02d}_{3:02d}-{4:02d}'.format(now.year, now.month,
    now.day, now.hour, now.minute)
    

class Constant(object):
    """ Class for a physical/astronomical constant
    """
    def __init__(self, cgs, SI, unitscgs='', info='No available description'):
        self.cgs = cgs
        self.SI  = SI
        self.unitscgs = unitscgs
        self.info = info

    def __repr__(self):
        return str('{:.7e} in {} (cgs)'.format(self.cgs, self.unitscgs))


def fltTxtOccur(s, lines, n=1, seq=1, after=False):
    """ Return the seq-th float of the line after the n-th
    occurrence of `s` in the array `lines`. 

    INPUT: s=string, lines=array of strings, n/seq=int (starting at 1)

    OUTPUT: float"""
    fltregex = r'[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?'
    if after:
        occur = [x[x.find(s)+len(s):] for x in lines if x.find(s) > -1]
    else:
        occur = [x for x in lines if x.find(s) > -1]
    out = _np.NaN
    if len(occur) >= n:
        occur = occur[n-1]
        out = _re.findall(fltregex, occur)[seq-1]
    return float(out)

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
    if _os.path.exists(fold) == False:
        _os.system('mkdir {}'.format(fold))
    return

def strrep(seq, n, newseq):
    """ Insert `newseq` at position `n` of the string `seq`.

    seq[n] = seq[:n]+newseq+seq[n+1:]

    Note: the string at the position `n` is replaced!!!

    OUTPUT: string
    """
    return seq[:n]+newseq+seq[n+1:]

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
    avg = _np.average(values, weights=1/sigma)
    return (avg, _np.sqrt(_np.sum(sigma**2))/len(values) )

def find_nearest(array,value, bigger=None):
    """ Find nearest VALUE in the array and return it.

    INPUT: array, value

    OUTPUT: closest value (array dtype)
    """
    if bigger == None:
        array = _np.array(array)
        idx = (_np.abs(array-value)).argmin()
        found = array[idx]
    elif bigger == True:
        found = min(x for x in array if x > value)        
    elif bigger == False:
        found = max(x for x in array if x < value)
    else:
        print('# ERROR at bigger!!')
    return found

def bindata(x, y, yerr, nbins, xrange=None):
    """
    Return the weighted binned data.

    if yerr == None:
        yerr = _np.ones(shape=_np.shape(x))

    INPUT: x, y, err - arrays with the same shape (they don't need to be sorted);
    nbins=int, xrange=[xmin, xmax]

    OUTPUT: xvals, yvals, new_yerr (arrays)
    """
    if xrange == None:
        max = _np.max(x)
        min = _np.min(x)
    else:
        min,max = xrange
    shift = (max-min) / (nbins-1)
    tmpx = _np.arange(nbins)*shift+min
    tmpy = _np.zeros(nbins)
    tmpyerr = _np.zeros(nbins)
    idx = _np.zeros(nbins, dtype=bool)
    for i in range(nbins):
        selx = _np.where( abs(x-tmpx[i])<=shift/2. )
        if len(selx[0]) >= 1:
            tmpy[i],tmpyerr[i] = wg_avg_and_std(y[selx], yerr[selx])
            idx[i] = True
    return (tmpx[idx], tmpy[idx], tmpyerr[idx])

def trimpathname(file):
    """Trim the full path string to return path and filename.

    INPUT: full file path

    OUTPUT: folder path, filename (strings)"""
    return [ file[:file.rfind('/')+1], file[file.rfind('/')+1:] ]

def rmext(name):
    """Remove the extension of a filename.
    Criteria: last `.` sets the extension.

    INPUT: filename (string)

    OUTPUT: filename without extension (string) """
    i = name.rfind('.')
    if i == -1:
        return name
    return name[:name.rfind('.')]

def normgauss(sig, x=None, xc=0.):
    """Normalized Gaussian function.

    INPUT: sigma value (float), x (array), xc (float)

    OUTPUT: yvals (array)
    """
    if x == None:
        x = _np.linspace(-5*sig,5*sig,101)
    return 1./(2*_np.pi)**.5/sig*_np.exp(-(x-xc)**2./(2*sig**2.))

def normbox(hwidth, x=None, xc=0.):
    """Normalized Box function.

    INPUT: half-width (number), x (array), xc (float)

    OUTPUT: yvals (array)
    """
    if x == None:
        x = _np.linspace(-hwidth*2+xc, hwidth*2+xc, 101)
    y = _np.zeros(len(x))
    idx = _np.where(_np.abs(x-xc) <= hwidth)
    y[idx] = _np.zeros(len(y[idx]))+2./hwidth 
    return y

def convnorm(x, arr, pattern):
    """Do the convolution of arr with pattern.
    Vector x is required for normalization. Its length must be odd!

    INPUT: units space (x, array), original array (arr), pattern (array)

    OUTPUT: array
    """
    if (len(x)%2 == 0) or (len(x) != len(pattern)):
        print('# Wrong format of x and/or pattern arrays!')
        return None
    dx = (x[-1]-x[0])/(len(x)-1.)
    cut = len(x)/2
    return _np.convolve(pattern, arr)[cut:-cut]*dx

def BBlbd(T,lbd=None):
    """ Black body radiation as function of lambda. CGS units (erg s-1 sr−1 cm−3).

    INPUT: lambda vector in cm. If None, lbd = 1000-10000 Angs."""
    if lbd==None:
        lbd = _np.arange(1000, 10000, 100)*1e-8 #Angs -> cm
    ft = h.cgs*c.cgs/(lbd*kB.cgs*T)
    return 2*h.cgs*c.cgs**2/(lbd**5*(_np.exp(ft)-1))

def recsearch(path, star, fstr):
    """
    Do a recursive search in PATH, looking inside `*star*` folders for 
    matching `fstr` (fullpath/night/star/files structure).

    Alternatively, one can use _os.walk(path).

    INPUT: path (string), folder (star, string), files (fstr, string)

    OUTPUT: list of strings
    """
    outfilelist = []
    nights = [o for o in _os.listdir(path) if _os.path.isdir('{}/{}'.format(path,\
    o))]
    for night in nights:
        targets = [o for o in _os.listdir('%s/%s' % (path,night)) if
            _os.path.isdir('%s/%s/%s' % (path,night,o))]
        for target in targets:
            if target.find(star)>-1:
                files = _glob('%s/%s/%s/*%s' % (path,night,target,fstr))
                for file in files:
                    outfilelist += [file]
    return outfilelist

def fracday2hms(frac):
    """Enter fraction of a day (e.g., MJD) and return integers of hour, min,
    sec.

    INPUT: float

    OUTPUT: hour, min, sec (int, int, int)
    """
    hh = frac*24
    if int(hh) > 0:
        mm = hh%int(hh)
        hh = int(hh)
    else:
        mm = hh
        hh = 0
    mm = mm*60
    if int(mm) > 0:
        ss = mm%int(mm)
        mm = int(mm)
    else:
        ss = mm
        mm = 0
    ss = int(round(ss*60))
    return hh,mm,ss

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
    #check sanity of dtstart
    if dtstart is None:
        dtstart = _dt.datetime(*_jdcal.jd2gcal(_jdcal.MJD_0,mjd0)[:3]).date()
        mjdst = _jdcal.gcal2jd(dtstart.year,dtstart.month,dtstart.day)[1]
    else:
        mjdst = _jdcal.gcal2jd(dtstart.year,dtstart.month,dtstart.day)[1]
        if mjdst < mjd0-1 or mjdst > mjd1:
            print('# Warning! Invalid "dtstart". Using mjd0.')
            dtstart = _dt.datetime(*_jdcal.jd2gcal(_jdcal.MJD_0,mjd0)[:3]).date()
    #define step 'position' and vector:
    basedata = [dtstart.year, dtstart.month, dtstart.day]
    dates =  []
    mjd = mjdst
    if step.upper() in ['Y','YY','YYYY']:
        i = 0
        while mjd < mjd1+1:
            dates +=  [_dt.datetime(*basedata).date()]
            basedata[i] += fact
            mjd = _jdcal.gcal2jd(*basedata)[1]
    elif step.upper() in ['M','MM']:
        i = 1
        while mjd < mjd1+1:
            dates += [_dt.datetime(*basedata).date()]
            basedata[i] += fact
            while basedata[i] > 12:
                basedata[0] += 1
                basedata[1] -= 12
            mjd = _jdcal.gcal2jd(*basedata)[1]
    elif step.upper() in ['D','DD']:
        i = 2
        daysvec = _np.arange(1,29,fact)
        if basedata[i] not in daysvec:
            j = 0
            while daysvec[j+1] < basedata[i]:
                j += 1
            daysvec += basedata[i]-daysvec[j]
            idx = _np.where(daysvec < 29)
            daysvec = daysvec[idx]
        else: 
            j = _np.where(daysvec == basedata[i])[0]
        while mjd < mjd1+1:
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

def cart2sph(x,y,z):
    """ Cartesian to spherical coordinates.

    INPUT: arrays of same length

    OUTPUT: arrays """
    r = _np.sqrt(x**2+y**2+z**2)
    idx = _np.where(r == 0)
    r[idx] = 1e-9
    th = _np.arccos(z/r)
    phi = _np.arctan2(y,x)
    #ind = (y<0) & (x<0)
    #phi[ind] = phi+_np.pi
    return r,th,phi
    
def sph2cart(r,th,phi):
    """  Spherical to Cartesian coordinates.

    INPUT: arrays of same length

    OUTPUT: arrays """
    x = r*_np.sin(th)*_np.cos(phi)
    y = r*_np.sin(th)*_np.sin(phi)
    z = r*_np.cos(th)
    return x,y,z

def cart_rot(x,y,z,ang_xy,ang_yz,ang_zx):
    """ Apply rotation in Cartesian coordinates.

    INPUT: 3 arrays of same length, 3 angles (float, in radians).

    OUTPUT: arrays """
    rotmtx = _np.array([ 
    [_np.cos(ang_zx)*_np.cos(ang_xy), -_np.cos(ang_yz)*_np.sin(ang_xy)+_np.sin(ang_yz)*_np.sin(ang_zx)*_np.cos(ang_xy),    _np.sin(ang_yz)*_np.sin(ang_xy)+_np.cos(ang_yz)*_np.sin(ang_zx)*_np.cos(ang_xy) ],
    [_np.cos(ang_zx)*_np.sin(ang_xy),    _np.cos(ang_yz)*_np.cos(ang_xy)+_np.sin(ang_yz)*_np.sin(ang_zx)*_np.sin(ang_xy), -_np.sin(ang_yz)*_np.cos(ang_xy)+_np.cos(ang_yz)+_np.sin(ang_zx)*_np.sin(ang_xy) ],
    [-_np.sin(ang_zx), _np.sin(ang_yz)*_np.cos(ang_zx), _np.cos(ang_yz)*_np.cos(ang_zx) ]
    ])
    vec = _np.array([x,y,z])
    return _np.dot(rotmtx,vec)

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

def interLinND(X, X0, X1, Fx, disablelog=False):
    """
    | N-dimensional linear interpolation in LOG space!!
    | Pay attention: Fx must always be > 0. If not, put disablelog=True.

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
    Xd = (X-X0)/(X1-X0)
    DX = _np.array([ [(1-x),x] for x in Xd ])
    #
    i = 0
    F = 0
    for prod in _product(*DX):
        if disablelog:
            F+= Fx[i]*_np.product(prod)
        else:
            F+= _np.log(Fx[i])*_np.product(prod)
        i+= 1
    #
    if not disablelog:
        return _np.exp(F)
    else:
        return F

ra2deg = _np.double(180./_np.pi)

ls = ['-','--',':','-.']

#Constants
c = Constant(2.99792458e10, 299792458., 'cm s-1', 'speed of light in vacuum') #:
h = Constant(6.6260755e-27, 6.62606957e-34, 'erg s-1', 'Planck constant') #:
hbar = Constant(1.05457266e-27, 1.05457172534e-34, 'erg s', 'Planck constant/(2*pi)') #:
G = Constant(6.674253e-8, 6.674253e-11, 'cm3 g-1 s-2', 'Gravitational constant') 
e = Constant(4.8032068e-10, 1.60217657e-19,	'esu', 'Elementary charge')
ep0 = Constant(1., 8.8542e-12, '', 'Permittivity of Free Space')
me = Constant(9.1093897e-28, 9.10938291e-31, 'g', 'Mass of electron')
mp = Constant(1.6726231e-24, 1.67262178e-27, 'g', 'Mass of proton')
mn = Constant(1.6749286e-24, 1.674927351e-27, 'g', 'Mass of neutron')
mH = Constant(1.6733e-24, 1.6737236e-27, 'g', 'Mass of hydrogen')	
amu = Constant(1.6605402e-24, 1.66053892e-27, 'g', 'Atomic mass unit')
nA = Constant(6.0221367e23, 6.0221413e23, '', "Avagadro's number")
kB = Constant(1.3806504e-16, 1.3806488e-23, 'erg K-1', 'Boltzmann constant')
eV = Constant(1.6021764871e-12, 1.602176565e-19, 'erg', 'Electron volt')
a = Constant(7.5646e-15, 7.5646e-16, 'erg cm-3 K-4', 'Radiation density constant') 	
sigma = Constant(5.67051e-5, 5.670373e-8, 'erg cm-2 K-4 s-1', 'Stefan-Boltzmann constant')
alpha = Constant(7.29735308e-3, 7.2973525698e-3, '', 'Fine structure constant')
Rinf = Constant(109737.316, 10973731.6, 'cm', 'Rydberg constant') 	
sigT = Constant(6.65245854533e-25, 6.65245854533e-29, 'cm2', 'Thomson cross section')

au = Constant(1.49597870691e13, 149597871, 'cm', 'Astronomical unit')
pc = Constant(3.08567758e18, 3.08567758e16, 'cm', 'Parsec')
ly = Constant(9.4605284e17, 9.4605284e15, 'cm', 'Light year')		
Msun = Constant(1.9891e33, 1.9891e30, 'g', 'Solar mass')	
Rsun = Constant(6.961e10, 696100, 'cm', 'Solar radius')
Lsun = Constant(3.846e33,	3.846e26, 'erg s-1', 'Solar luminosity')
Tsun = Constant(5778., 5778., 'K', 'Solar Temperature')

yr =  Constant(3.15569e7, 3.15569e7, 'sec', 'year')

colors = ['Black','Blue','Green','red','orange','brown','purple','gray',
    'dodgerblue','lightgreen','tomato','yellow','peru','MediumVioletRed',
    'LightSteelBlue','cyan','darkred','olive']

### MAIN ###
if __name__ == "__main__":
    pass
