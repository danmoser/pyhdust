##!/usr/bin/env python
#-*- coding:utf-8 -*-
#Modified by D. Moser in 2014-10-17

"""
Physical constants

includes *basicfuncs*
"""

import numpy as np
import pyhdust.jdcal as jdcal
import datetime as dt
import os, re
from glob import glob

class Constant(object):
    """ Class for a physical/astronomical constant
    """
    def __init__(self, cgs, SI, cgsunits='', info='No available description'):
        self.cgs = cgs
        self.SI  = SI
        self.cgsunits = cgsunits
        self.info = info

    def __repr__(self):
        return str('{:.7e} in {} (cgs)'.format(self.cgs, self.cgsunits))


def fltTxtOccur(s,block,n=1):
    """ Return the first float of the line in 'block' with the n-th
    occurrence of 's' """
    regex=r'[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?'
    occur = [x for x in block if x.find(s) > -1]
    out = np.NaN
    if len(occur) >= n:
        occur = occur[n-1]
        out = re.findall(regex, occur)[0]
    return float(out)

def outfld(fold='hdt'):
    """
    Check and create (if necessary) an (sub)folder - generally used for output.
    """
    if os.path.exists(fold) == False:
        os.system('mkdir {}'.format(fold))

def strrep(seq, n, newseq):
    """ seq[n] = seq[:n]+newseq+seq[n+1:]

    Note: the string at the position `n` is replaced!!!

    @param PARAM: DESCRIPTION
    @return RETURN: DESCRIPTION
    """
    return seq[:n]+newseq+seq[n+1:]

def wg_avg_and_std(values, sigma):
    """
    Return the weighted average and standard deviation.

    values, sigma -- Numpy ndarrays with the same shape.

    This IS NOT the problem of Wikipedia > Weighted_arithmetic_mean >
        Weighted_sample_variance.

    average = np.average(values, weights=weights)
    #Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, np.sqrt(variance))
    """
    avg = np.average(values, weights=1/sigma)
    return (avg, np.sqrt(np.sum(sigma**2))/len(values) )

def find_nearest(array,value):
    """ Find nearest VALUE in the array and return it
    """
    array = np.array(array)
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def bindata(x, y, yerr, nbins, xrange=None):
    """
    Return the weighted binned data.

    if yerr == None:
        yerr = np.ones(shape=np.shape(x))

    x, y, err - Numpy arrays with the same shape.
    They don't need to be sorted.
    nbins - integer
    """
    if xrange == None:
        max = np.max(x)
        min = np.min(x)
    else:
        min,max = xrange
    shift = (max-min) / (nbins-1)
    tmpx = np.arange(nbins)*shift+min
    tmpy = np.zeros(nbins)
    tmpyerr = np.zeros(nbins)
    idx = np.zeros(nbins, dtype=bool)
    for i in range(nbins):
        selx = np.where( abs(x-tmpx[i])<=shift/2. )
        if len(selx[0]) >= 1:
            tmpy[i],tmpyerr[i] = wg_avg_and_std(y[selx], yerr[selx])
            idx[i] = True
    return (tmpx[idx], tmpy[idx], tmpyerr[idx])

def trimpathname(file):
    """Trim the full path string to return path and filename."""
    return [ file[:file.rfind('/')+1], file[file.rfind('/')+1:] ]

def rmext(name):
    """Remove the extension of a filename
    Criteria: last `.` sets the extension"""
    i = name.rfind('.')
    if i == -1:
        return name
    return name[:name.rfind('.')]

def normgauss(sig, x=None, xc=0.):
    """Normalized Gaussian function
    """
    if x == None:
        x = np.linspace(-5*sig,5*sig,101)
    return 1./(2*np.pi)**.5/sig*np.exp(-(x-xc)**2./(2*sig**2.))

def normbox(width, x=None, xc=0.):
    """Normalized Box function
    """
    if x == None:
        x = np.linspace(-width, width, 101)
        xc = 0.
    y = np.zeros(len(x))
    idx = np.where(np.abs(x-xc) <= width/2.)
    y[idx] = np.zeros(len(y[idx]))+1./width 
    return y

def convnorm(x, arr, pattern):
    """Do the convolution of arr with pattern.
    Vector x is required for normalization. Its length must be odd!
    """
    if (len(x)%2 == 0) or (len(x) != len(pattern)):
        print('# Wrong format of x and/or pattern arrays!')
        return None
    dx = (x[-1]-x[0])/(len(x)-1.)
    cut = len(x)/2
    return np.convolve(pattern, arr)[cut:-cut]*dx

def recsearch(path, star, fstr):
    """
    Do a recursive search in PATH, looking inside `*star*` folders for 
    matching `*files` (fullpath/night/star/files structure).

    Alternatively, one can use os.walk(path).
    """
    outfilelist = []
    nights = [o for o in os.listdir(path) if os.path.isdir('{}/{}'.format(path,\
    o))]
    for night in nights:
        targets = [o for o in os.listdir('%s/%s' % (path,night)) if
            os.path.isdir('%s/%s/%s' % (path,night,o))]
        for target in targets:
            if target.find(star)>-1:
                files = glob('%s/%s/%s/*%s' % (path,night,target,fstr))
                for file in files:
                    outfilelist += [file]
    return outfilelist

def fracday2hms(frac):
    """Enter fraction of a day (e.g., MJD) and return integers of hour, min,
    sec.
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

    dtstart (optional) is expected to be in datetime.datetime.date() format
    [i.e., datetime.date(yyyy, m, d)]

    fact must be an integer
    """
    #check sanity of dtstart
    if dtstart is None:
        dtstart = dt.datetime(*jdcal.jd2gcal(jdcal.MJD_0,mjd0)[:3]).date()
        mjdst = jdcal.gcal2jd(dtstart.year,dtstart.month,dtstart.day)[1]
    else:
        mjdst = jdcal.gcal2jd(dtstart.year,dtstart.month,dtstart.day)[1]
        if mjdst < mjd0-1 or mjdst > mjd1:
            print('# Warning! Invalid "dtstart". Using mjd0.')
            dtstart = dt.datetime(*jdcal.jd2gcal(jdcal.MJD_0,mjd0)[:3]).date()
    #define step 'position' and vector:
    basedata = [dtstart.year, dtstart.month, dtstart.day]
    dates =  []
    mjd = mjdst
    if step.upper() in ['Y','YY','YYYY']:
        i = 0
        while mjd < mjd1+1:
            dates +=  [dt.datetime(*basedata).date()]
            basedata[i] += fact
            mjd = jdcal.gcal2jd(*basedata)[1]
    elif step.upper() in ['M','MM']:
        i = 1
        while mjd < mjd1+1:
            dates += [dt.datetime(*basedata).date()]
            basedata[i] += fact
            while basedata[i] > 12:
                basedata[0] += 1
                basedata[1] -= 12
            mjd = jdcal.gcal2jd(*basedata)[1]
    elif step.upper() in ['D','DD']:
        i = 2
        daysvec = np.arange(1,29,fact)
        if basedata[i] not in daysvec:
            j = 0
            while daysvec[j+1] < basedata[i]:
                j += 1
            daysvec += basedata[i]-daysvec[j]
            idx = np.where(daysvec < 29)
            daysvec = daysvec[idx]
        else: 
            j = np.where(daysvec == basedata[i])[0]
        while mjd < mjd1+1:
            dates += [dt.datetime(*basedata).date()]
            j += 1
            if j == len(daysvec):
                j = 0
                basedata[1] += 1
                if basedata[1] == 13:
                    basedata[1] = 1
                    basedata[0] += 1
            basedata[i] = daysvec[j]
            mjd = jdcal.gcal2jd(*basedata)[1]
    else:
        print('# ERROR! Invalid step')
        raise SystemExit(1)
    return dates

def cart2sph(x,y,z):
    """ ### GEOMETRY ### """
    r = np.sqrt(x**2+y**2+z**2)
    th = np.arccos(z/r)
    phi = np.arctan2(y,x)
    #ind = (y<0) & (x<0)
    #phi[ind] = phi+np.pi
    return r,th,phi
    
def sph2cart(r,th,phi):
    """ ### GEOMETRY ### """
    x = r*np.sin(th)*np.cos(phi)
    y = r*np.sin(th)*np.sin(phi)
    z = r*np.cos(th)
    return x,y,z

def cart_rot(x,y,z,ang_xy,ang_yz,ang_zx):
    """ ### GEOMETRY ### """
    rotmtx = np.array([ 
    [np.cos(ang_zx)*np.cos(ang_xy), -np.cos(ang_yz)*np.sin(ang_xy)+np.sin(ang_yz)*np.sin(ang_zx)*np.cos(ang_xy),    np.sin(ang_yz)*np.sin(ang_xy)+np.cos(ang_yz)*np.sin(ang_zx)*np.cos(ang_xy) ],
    [np.cos(ang_zx)*np.sin(ang_xy),    np.cos(ang_yz)*np.cos(ang_xy)+np.sin(ang_yz)*np.sin(ang_zx)*np.sin(ang_xy), -np.sin(ang_yz)*np.cos(ang_xy)+np.cos(ang_yz)+np.sin(ang_zx)*np.sin(ang_xy) ],
    [-np.sin(ang_zx), np.sin(ang_yz)*np.cos(ang_zx), np.cos(ang_yz)*np.cos(ang_zx) ]
    ])
    vec = np.array([x,y,z])
    return np.dot(rotmtx,vec)

def readrange(file, i0, ie):
    """ Read a specific range of lines with minimal memory use.

    Note that i == n-1 for the nth line."""
    lines = []
    fp = open(file)
    for i, line in enumerate(fp):
        if i >= i0:
            lines += [line]
        if i >= ie:
            break
    fp.close()
    return lines

#Physical constants
c = Constant(2.99792458e10, 299792458., 'cm s-1', 'speed of light in vacuum')
h = Constant(6.6260755e-27, 6.62606957e-34, 'erg s-1', 'Planck constant')
hbar = Constant(1.05457266e-27, 1.05457172534e-34, 'erg s', 'Planck constant/(2*pi)')
G = Constant(6.67259e-8, 6.67384e-11, 'cm3 g-1 s-2', 'Gravitational constant') 
e = Constant(4.8032068e-10, 1.60217657e-19,	'esu', 'Elementary charge')
ep0 = Constant(1., 8.8542e-12, '', 'Permittivity of Free Space')
me = Constant(9.1093897e-28, 9.10938291e-31, 'g', 'Mass of electron')
mp = Constant(1.6726231e-24, 1.67262178e-27, 'g', 'Mass of proton')
mn = Constant(1.6749286e-24, 1.674927351e-27, 'g', 'Mass of neutron')
mH = Constant(1.6733e-24, 1.6737236e-27, 'g', 'Mass of hydrogen')	
amu = Constant(1.6605402e-24, 1.66053892e-27, 'g', 'Atomic mass unit')
nA = Constant(6.0221367e23, 6.0221413e23, '', "Avagadro's number")
kB = Constant(1.380658e-16, 1.3806488e-23, 'erg K-1', 'Boltzmann constant')
eV = Constant(1.6021772e-12, 1.602176565e-19, 'erg', 'Electron volt')
a = Constant(7.5646e-15, 7.5646e-16, 'erg cm-3 K-4', 'Radiation density constant') 	
sigma = Constant(5.67051e-5, 5.670373e-8, 'erg cm-2 K-4 s-1', 'Stefan-Boltzmann constant')
alpha = Constant(7.29735308e-3, 7.2973525698e-3, '', 'Fine structure constant')
Rinf = Constant(109737.316, 10973731.6, 'cm', 'Rydberg constant') 	
sigT = Constant(6.65245854533e-25, 6.65245854533e-29, 'cm2', 'Thomson cross section')

au = Constant(1.496e13, 149597871, 'cm', 'Astronomical unit')
pc = Constant(3.086e18, 3.08567758e16, 'cm', 'Parsec')
ly = Constant(9.463e17, 9.4605284e15, 'cm', 'Light year')		
Msun = Constant(1.99e33, 1.9891e30, 'g', 'Solar mass')	
Rsun = Constant(6.96e10, 695800, 'cm', 'Solar radius')
Lsun = Constant(3.9e33,	3.846e26, 'erg s-1', 'Solar luminosity')
Tsun = Constant(5.780e3, 5778., 'K', 'Solar Temperature')

yr =  Constant(3.15569e7, 3.15569e7, 'sec', 'year')

ra2deg = np.double(180./np.pi)

colors = ['Black','Blue','Green','red','orange','brown','purple','gray',
    'dodgerblue','lightgreen','tomato','yellow','peru','MediumVioletRed',
    'LightSteelBlue','cyan','darkred','olive']

ls = ['-','--',':','-.']

