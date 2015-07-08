# -*- coding:utf-8 -*-

"""
PyHdust main module: Hdust tools.

:license: GNU GPL v3.0 (https://github.com/danmoser/pyhdust/blob/master/LICENSE)
"""
import os as _os
import time as _time
import datetime as _datetime
import math as _math
import numpy as _np
import re as _re
import struct as _struct
from glob import glob as _glob
import pyhdust.phc as _phc
import pyhdust.jdcal as _jdcal

try:
    import matplotlib.pyplot as _plt
    from scipy import interpolate as _interpolate
except:
    print('# Warning! matplotlib and/or scipy module not installed!!!')

__version__ = 0.961
__release__ = "Beta"
__author__ = "Daniel Moser"
__email__ = "dmfaes@gmail.com"


def readscr(file):
    '''
    Read source generated with `ref_estrela.txt`.

    OUTPUT: M, Req and TP (2*solar units and K).
    '''
    f0 = open(file)
    lines = f0.readlines()
    f0.close()

    n = int(_phc.fltTxtOccur('STAR =',lines,n=1))
    M = _phc.fltTxtOccur('M =',lines,n=n)
    Rp = _phc.fltTxtOccur('R_pole =',lines,n=n)
    if n == 2:
        ob = _phc.fltTxtOccur('R_eq/R_pole =',lines,n=1)
        Tp = _phc.fltTxtOccur('Teff_pole =',lines,n=1)
    else:
        W = _phc.fltTxtOccur('W =',lines,n=1)
        bet = _phc.fltTxtOccur('Beta_GD =',lines,n=1)
        L = _phc.fltTxtOccur('L =',lines,n=n)
        wfrac = _np.sqrt(27./8*(1+0.5*W**2)**-3*W**2)
        ob, Tp, A = rotStar(Tp=L, M=M, rp=Rp, beta=bet, wfrac=wfrac, quiet=True,
        LnotTp=True)
    #print M,Rp*ob,Tp
    return M,Rp*ob,Tp


def n0toSigma0(n0, M, Req, f, Tp, mu):
    """ VDD Steady-State conversion between `n0` to `Sigma0`.

    INPUT: n0 (float), M, Req (mass and equatorial radius, Solar units),
    fraction and polar temperature (0-1, Kelvin) and mu molecular weight
    (0.5-1.0)

    OUTPUT: float (g cm-2) """
    rho0 = n0 * mu * _phc.mH.cgs
    a = (_phc.kB.cgs * f * Tp / mu / _phc.mH.cgs) ** .5
    sig0 = (2 * _np.pi) ** .5 * a / (_phc.G.cgs * M * _phc.Msun.cgs / (Req * _phc.Rsun.cgs)) ** .5 * \
           Req * _phc.Rsun.cgs * rho0
    return sig0


def n0toMdot(n0, M, Req, f, Tp, mu, alpha, R0):
    """ VDD Steady-State conversion between `n0` to `Mdot`.

    INPUT: n0 (float), M, Req (mass and equatorial radius, Solar units),
    fraction and polar temperature (0-1, Kelvin), mu molecular weight
    (0.5-1.0), alpha (viscous parameter), R0 ("truncation" radius, Solar unit).

    OUTPUT: float (Msun yr-1)"""
    Req = Req * _phc.Rsun.cgs
    R0 = R0 * _phc.Rsun.cgs
    M = M * _phc.Msun.cgs
    rho0 = n0 * mu * _phc.mH.cgs
    a = (_phc.kB.cgs * f * Tp / mu / _phc.mH.cgs) ** .5
    Mdot = 3 * _np.pi * (2 * _np.pi) ** .5 * alpha * a ** 3. / (_phc.G.cgs * M / Req) * rho0 * Req ** 2. * \
           ((R0 / Req) ** .5 - 1) ** -1.
    return Mdot / _phc.Msun.cgs * _phc.yr.cgs


def plotMJDdates(spec=None, pol=None, interf=None, limits=None):
    """
    Plot dates from spec (Class), pol (routines) and interf (ESO query)

    This need to be polished !!!!
    """
    fig, ax = _plt.subplots()
    spec = 'data_aeri_splot.txt'
    if spec is not None:
        spJD = _np.loadtxt(spec)
        spJD = spJD[:, 0]
        y = [0. for JD in spJD]
        ax.plot(spJD, y, marker='d', color='lightgray', ls='')
        # yerr = [ [1. for JD in spJD], [1. for JD in spJD] ]
        #ax.errorbar(spJD, y, yerr, marker='o', color='blue', ls='')

    pol = 'pol_aeri.log'
    if pol is not None:
        polJD = _np.loadtxt(pol, dtype=str)
        polJD = polJD[:, 9]
        polJD = _np.array(polJD, dtype=float) - 2400000.5
        y = [-.5 for JD in polJD]
        ax.plot(polJD, y, marker='o', color='gray', ls='')
        # yerr = [ [.5 for JD in polJD], [1.5 for JD in polJD] ]
        #ax.errorbar(polJD, y, yerr, marker='x', color='green', ls='')

    interf = 'interf_aeri.txt'
    if interf is not None:
        intJD = _np.loadtxt(interf, dtype=str, delimiter=',', skiprows=1)
        intJD = _np.array(intJD[:, -2], dtype=float)
        y = [.5 for JD in intJD]
        ax.plot(intJD, y, marker='s', color='darkgrey', ls='')
        # yerr = [ [1.5 for JD in intJD], [.5 for JD in intJD] ]
        #ax.errorbar(intJD, y, yerr, marker='s', color='red', ls='')

    limits = (56100., 56750.)
    if limits is None:
        mjd0, mjd1 = ax.get_xlim()
    else:
        mjd0, mjd1 = limits
        ax.set_xlim(limits)
    ticks = _phc.gentkdates(mjd0, mjd1, 3, 'm', dtstart=dt.datetime(2012, 7, 1). \
                            date())
    mjdticks = [_jdcal.gcal2jd(date.year, date.month, date.day)[1] for date in \
                ticks]
    ax2 = ax.twiny()
    ax2.set_xlim(limits)
    ax2.set_xticks(mjdticks)
    ax2.set_xlabel('Civil date')
    ax2.set_xticklabels([date.strftime("%d %b %y") for date in ticks])
    _plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45)
    ax.set_yticklabels([])
    ax.set_xlabel('MJD')
    return


def hdtpath():
    """
    Return the path os the module.

    >>> hdt.hdtpath()
    /home/user/Scripts/pyhdust/
    """
    fulldir = __file__[:__file__.rfind('/') + 1]
    return fulldir[:fulldir[:-1].rfind('/') + 1]


def doFilterConv(x0, y0, filter, pol=False):
    """
    Return the convolved filter total flux for a given flux profile y0,
    at wavelengths x0 (um).

    INPUT: x0 lambda array (um), y0 flux array, filter (string)

    OUTPUT: summed flux (y0 units)
    """
    fdat = _np.loadtxt('{0}/filters/{1}.dat'.format(hdtpath(), filter.lower()), \
                       skiprows=1)
    fdat[:, 0] /= 10000.  # from Angs to microns
    idx = _np.where((x0 >= fdat[0, 0]) & (x0 <= fdat[-1, 0]))
    x0 = x0[idx]
    y0 = y0[idx]
    # interpfunc = interpolate.interp1d(fdat[:,0], fdat[:,1], kind='linear')
    interpfunc = _interpolate.InterpolatedUnivariateSpline(fdat[:, 0], fdat[:, 1])

    if not pol:
        y = interpfunc(x0) * y0
        return _np.trapz(y, x0)
    else:
        return phc.wg_avg_and_std(y0, 1/interpfunc(x0))[0]


def doPlotFilter(pref, obs, filter, fsed2data, pol=False):
    """
    pref = output prefix; obs = integer; filter = single string
    """
    x0 = fsed2data[obs, :, 2]
    if pol:
        y0 = fsed2data[obs, :, 7]
        savename = 'pol_{0}_{1}_{2}'.format(pref, obs, filter)
    else:
        y0 = fsed2data[obs, :, 3] / x0
        savename = '{0}_{1}_{2}'.format(pref, obs, filter)

    fdat = _np.loadtxt('{0}filters/{1}.dat'.format(hdtpath(), filter.lower()), \
                       skiprows=1)
    fdat[:, 0] /= 10000.
    # interpfunc = interpolate.interp1d(fdat[:,0], fdat[:,1], kind='linear')
    interpfunc = _interpolate.InterpolatedUnivariateSpline(fdat[:, 0], fdat[:, 1])

    idx = _np.where((x0 >= fdat[0, 0]) & (x0 <= fdat[-1, 0]))
    x0 = x0[idx]
    y0 = y0[idx]
    y = interpfunc(x0) * y0  #/_np.sum( interpfunc(x0) )

    fig, ax = _plt.subplots()
    ax.plot(x0, y0, label='SED')
    #ax.plot(fdat[:,0], fdat[:,1], label='Filter')
    ax.plot(x0, interpfunc(x0) * y0, label='Convolved')
    ax.plot(fdat[:, 0], fdat[:, 1] * _np.max(y0), label='Filter (eff.)')
    ax.set_title(savename)
    ax.legend()
    fig.savefig(savename + '.png', transparent=True)
    _plt.close()
    print('# Saved {0}.png !'.format(savename))
    return


def sed2info(file):
    """
    Read info from SED2 file.

    INPUT: file (path string)

    OUTPUT: nlbd, nobs, Rstar, Rwind (as floats)
    """
    f0 = open(file, 'r')
    fcont = f0.readlines()
    f0.close()
    info = ''
    i = 0
    while info == '':
        if fcont[i][0] != '%':
            info = _np.array(fcont[i].split(), dtype=float)
        i += 1
    return info


def readfullsed2(file):
    """
    Read data from FULLSED2 file.

    INPUT: file (path string)

    OUTPUT: array[nobs,nlbd,-1]
        number of columns from SED2file replaces "-1"
    """
    nlbd, nobs, Rstar, Rwind = sed2info(file)
    sed2data = _np.loadtxt(file, skiprows=5)
    sed2data = sed2data.reshape((nobs, nlbd, -1))
    return sed2data


def readsed2(file):
    """
    Read data from SED2 file.

    INPUT: file (path string)

    OUTPUT: array[nobs*nlbd,-1]
        number of columns from SED2file replaces "-1".
        NOTE: this format is different of `readfullsed2`.
    """
    nlbd, nobs, Rstar, Rwind = sed2info(file)
    sed2data = _np.loadtxt(file, skiprows=1)
    return sed2data


def chkObsLog(path=None, nights=None, badweath=None):
    """ Check if there is data for all nights with observations.

    If not, check if the night is in the list of night lost due to bad weather.

    If no data and no bad weather info is registered, prints an error.

    If the night is included as bad weather and is not in night list, prints a
    warning.
    """
    if path == None:
        path = _os.getcwd()
    if nights == None:
        nights = '{0}refs/noites.txt'.format(hdtpath())
    lnights = _np.loadtxt(nights, dtype=str)
    if badweath == None:
        badweath = '{0}refs/maltempo.txt'.format(hdtpath())
    lbadweath = _np.loadtxt(badweath, dtype=str)
    for night in lnights:
        if night in lbadweath:
            pass
        elif not _os.path.exists(night):
            print('# ERROR! {0} has no data and was not lost for bad weather!'. \
                  format(night))
    flds = [fld for fld in _os.listdir('{0}'.format(path)) if \
            _os.path.isdir(_os.path.join('{0}'.format(path), fld))]
    for fld in flds:
        if fld not in lnights:
            print('# Warning! Night {0} is not recorded as OPD night!'. \
                  format(fld))
            print('# Update the file {0}'.format(nights))
    for night in lbadweath:
        if night not in lnights:
            print('# Warning! Bad weather {0} is not recorded as OPD night!'. \
                  format(night))
            print('# Probably it is a spec night.')
    return


def readtemp(tfile, quiet=False):
    """ Read *.temp file

    - ncr = número de células da simulação na coordenada radial
    - ncmu = número de células da simulação na coordenada latitudinal
    - ncphi = número de células da simulação na coordenada azimutal
    - nLTE = number of atomic LTE levels
    - nNLTE = number of atomic NLTE levels
    - beta = flare disk parameter
    - Rstar = raio da estrela
    - Ra = raio máximo da simulação, onde acabam todas as poeiras
    - pcrc = coordenadas das células em raio
    - pcmuc = coordenadas das células em mu
    - pcphic = coordenadas das células em phi
    - pcr = distância entre as células em raio
    - pcmu = distância entre as células em mu
    - pcphi = distância entre as células em phi

    OUTPUT = ncr,ncmu,ncphi,nLTE,nNLTE,Rstar,Ra,beta,data,pcr,pcmu,pcphi
    """
    f = open(tfile).read()
    ixdr=0
    ncr, ncmu, ncphi, nLTE, nNLTE = _struct.unpack('>5l', f[ixdr:ixdr+4*5])
    ixdr+=4*5
    Rstar, Ra, beta = _struct.unpack('>3f', f[ixdr:ixdr+4*3])
    ixdr+=4*3
    #~ 
    rlen = (nLTE+6)*ncr*ncmu*ncphi
    data = _struct.unpack('>{0}f'.format(rlen), f[ixdr:ixdr+4*rlen])
    ixdr+=4*rlen
    data = _np.reshape(data, (nLTE+6,ncr,ncmu,ncphi), order='F')
    #~
    #this will check if the XDR is finished.
    if ixdr == len(f):
        if not quiet:
            print('# XDR {0} completely read!'.format(tfile))
    else:
        print('# Warning: XDR {0} not completely read!'.format(tfile))
        print('# length difference is {0}'.format( (len(f)-ixdr)/4 ) )
    #~
    pcrc = data[0,:,0,0]
    pcr = _np.zeros(ncr+1)
    pcr[0] = Rstar
    for icr in range(1,ncr+1):
        pcr[icr] = pcr[icr-1] + 2*(pcrc[icr-1]-pcr[icr-1])
    pcmu = _np.zeros((ncmu+1,ncr))
    #~ 
    for icr in range(0, ncr):
        pcmuc = data[1,icr,:,0]
        pcmu[0,icr] = -1.
        for icmu in range(1,ncmu+1):
            pcmu[icmu,icr] = pcmu[icmu-1,icr] + 2.*(pcmuc[icmu-1]-pcmu[icmu-1,icr])
    #~
    pcphic = data[2,0,0,:]
    pcphi = _np.zeros(ncphi+1)
    pcphi[0] = 0.
    for icphi in range(1,ncphi+1):
        pcphi[icphi] = pcphi[icphi-1] + 2*(pcphic[icphi-1]-pcphi[icphi-1])
    #~ 
    return ncr,ncmu,ncphi,nLTE,nNLTE,Rstar,Ra,beta,data,pcr,pcmu,pcphi


def readdust(tfile):
    """ TBD!!

    - ntip, = número de tipos de poeira determinados para a simulação (composição)
    - na, = número do tipo da poeira (tamanhos)
    - NdustShells, = número de camadas de poeiras da simulação
    - Rdust, = raio onde está(ão) a(s) camada(s)
    - Tdestruction, = temperatura na qual os grãos são evaporados
    - Tdust, = temperatura da poeira numa da posição da grade da simulação (r,phi,mu)
    - lacentro, = controla os tipos e tamanhos das poeiras

    """
    return


def plotdust(tfile):
    """ TBD!!

    For more info, see `readdust` help. 
    
    """
    return


def gentemplist(tfile, tfrange=None, avg=True):
    """
    `tfile` = file name or file prefix to be plotted. If `tfrange` (e.g.,
    tfrange=[20,24] is present, it you automatically plot the interval.

    `avg` = include tfile_avg.temp file (if exists).

    OUTPUT = file list, label list
    """
    if tfrange is None:
        tfrange = [int(tfile.replace('.temp','')[-2:])]
    ltfile = []
    ltlabel = []
    for tn in range(int(tfrange[0]), int(tfrange[-1])+1):
        if tfile.find('.temp') > 0:
            ltfile+= [tfile.replace(tfile[-7:-5], '{0:02d}'.format(tn))]
        else:
            ltfile+= [tfile+'{0:02d}.temp'.format(tn)]
        ltlabel += [str(tn)]
    tfile = ltfile[-1].replace('.temp', '_avg.temp')
    if _os.path.exists(tfile) and avg:
        ltfile+= [tfile]
        ltlabel += ['Avg.']
    return ltfile, ltlabel


def plottemp(tfiles, philist=[0], interpol=False, xax=0, fmts=['png'],
    outpref=None, tlabels=None, title=None):
    """
    .. code::

        >>> import pyhdust as hdt
        >>> tfiles, tlabels = hdt.gentemplist('bestar2.02/mod01/mod01b33.temp', tfrange=[30,33])
        >>> hdt.plottemp(tfiles, tlabels=tlabels)

    .. image:: _static/hdt_plottemp.png
        :width: 512px
        :align: center
        :alt: hdt.plottemp example

    `tfile` = filenames list.

    `xax` = 0: log10(r/R-1), 1: r/R; 2: 1-R/r

    `outpref` = prefix of the output images

    `fmts` = format os the output images.

    If interpol==True, what will be plotted is the population along rays of
    a given latitude. The latitudes are defined in array muplot below.

    If interpola==False, what will be plotted is the population for a given mu
    index as a function of radius, starting with index ncmu/2(midplane) + plus

    OUTPUT = ...
    """
    if isinstance(tfiles, basestring):
        tfiles = [tfiles]
    if title is None:
        title = tfiles[-1]
    if interpol:
        print('# Interpol option currently is not available.')
        raise SystemExit(0)
    if xax not in [0,1,2]:
        print('# Invalid `xax` option. Try again.')
        raise SystemExit(0)
    #~
    fig, axs = _plt.subplots(1, 1)#, figsize=(21./3,29.7/3), sharex=True)
    lev = 0
    np = 0
    rplus = 0
    for i in range(len(tfiles)):
        rtfile = tfiles[i]
        ncr,ncmu,ncphi,nLTE,nNLTE,Rstar,Ra,beta,data,pcr,pcmu,pcphi = \
        readtemp(rtfile, quiet=True)
        for phiidx in range(0,len(philist)):
            icphi = philist[phiidx]
            x = data[0,:,0,icphi]
            if (xax == 0): x = _np.log10(x/Rstar-1.); xtitle = r'$\log_{10}(r/R_*-1)$'
            elif (xax == 1): x = x/Rstar; xtitle = r'$r/R_*$'
            elif (xax == 2): x = 1.-Rstar/x; xtitle = r'$1-R_*/r$'
            y = data[3+lev,:,ncmu/2+np+rplus,icphi]
            y = y/1000.
            fmt = 'o:'
            if rtfile.find('avg') > 0:
                fmt = 'o-'
            if tlabels is not None:
                axs.plot(x, y, fmt, label=tlabels[i])
            else:
                axs.plot(x, y, fmt)
    #~
    axs.legend()
    axs.set_title(title)
    axs.set_xlabel(xtitle)
    axs.set_ylabel(r'Temperature (10$^3$ K)')
    if outpref == None:
        outpref = _phc.dtflag()
    for fmt in fmts:
        _plt.savefig('{0}.{1}'.format(outpref,fmt), transparent=True)
    _plt.close()
    return


def mergesed2(models, Vrots, path=None):
    """
    Merge all mod#/*.sed2 files into the fullsed file.
    
    It will check if all sed2info are the same (i.e., nlbd, nobs, Rstar, Rwind).
    If not, it ask if you want to continue (and receive and error).
    
    The presence of the SED file is not required anymore.
    The structure is set by the first broadband sed2 found.
    
    NO AVERAGE is coded (yet).

    INPUT: models lists (*.txt or *.inp), Vrots (array).

    OUTPUT: *files written (status printed).
    """
    sufbands = ['SED', 'UV', 'IR', 'NIR', 'BALMER', 'PASCHEN', 'CM', 'MM',  # 0-7
                'J', 'H', 'K', 'L', 'M', 'N', 'Q1', 'Q2']  # 8-14
    # wavelength in microns
    suflines = {'H12': .372300, 'H11': .373543, 'H10': .375122, 'H9': .377170,
    'H8': .379899, 'H7': .383649, 'H6': .389017, 'H5': .397120, 'Hd': .410289,
    'Hg': .434169, 'Hb': .486271, 'Ha': .656461, 'Br13':1.61137, 'Br12':1.6416,
    'Brg':2.166}

    for model in models:
        model = model.replace('.inp', '.txt')
        modfld, modelname = _phc.trimpathname(model)
        path = _phc.trimpathname(modfld[:-1])[0]
        if not _os.path.exists('{0}fullsed'.format(path)):
            _os.system('mkdir {0}fullsed'.format(path))
        #
        #~ modelname = modelname.replace('.txt','.inp')
        sed2data = _np.empty(0)
        sfound = []
        #Process broad-bands
        for suf in sufbands:
            file = modfld + '{0}_{1}.sed2'.format(suf, modelname.replace(".txt", ""))
            if _os.path.exists(file):
                sfound += [suf]
                newdata = readsed2(file)
                if len(sed2data) == 0:
                    sed2data = newdata.copy()
                    nlbd, nobs, Rstar, Rwind = sed2info(file)
                else:
                    #Check if the SED2 file has the info as the first file
                    if _np.product((nobs, Rstar, Rwind) == sed2info(file)[1:]) == 0:
                        key = ''
                        while key.upper() != Y:
                            print('# WARNING: {0} has different HDUST output!!!'. \
                                  format(modelname))
                            key = raw_input('Do you want do proceed? (y/other): ')
                    nlbd += sed2info(file)[0]
                    sed2data = _np.vstack((sed2data, newdata))
        #Process lines
        Vrot = Vrots[models.index(model)]
        for suf in suflines:
            file = modfld + '{0}_{1}_SEI.sed2'.format(suf, modelname.replace(".txt", ""))
            if _os.path.exists(file):
                sfound += [suf]
                newdata = readsed2(file)
                #print("# TRIMMING INPUT SPECTRUM {}".format(modelname))
                #print("# TO ACCOUNT FOR THE NON-ZERO ROTATION VELOCITY.")
                deltalbd = Vrot / _phc.c.cgs / 1e-5 * suflines[suf]
                mini = newdata[0, 2] + deltalbd
                maxi = newdata[-1, 2] - deltalbd
                #print deltalbd, Vrot
                idx = _np.where((newdata[:, 2] >= mini) & (newdata[:, 2] <= maxi))
                #print len(newdata), len(newdata[idx])
                ncut = len(newdata) - len(newdata[idx])
                newdata = newdata[idx]
                if len(sed2data) == 0:
                    sed2data = newdata.copy()
                    nlbd, nobs, Rstar, Rwind = sed2info(file)
                else:
                    #Check if the SED2 file has the info as the first file
                    if _np.product((nobs, Rstar, Rwind) == sed2info(file)[1:]) == 0:
                        key = ''
                        while key.upper() != Y:
                            print('# WARNING: {0} has different HDUST input!!!'. \
                                  format(modelname))
                            key = raw_input('Do you want do proceed? (y/other): ')
                    idx = _np.where((sed2data[:, 2] < mini) | (sed2data[:, 2] > maxi))
                    ncut += len(sed2data) - len(sed2data[idx])
                    sed2data = sed2data[idx]
                    nlbd += sed2info(file)[0] - ncut / sed2info(file)[1]
                    sed2data = _np.vstack((sed2data, newdata))

        if len(sfound) > 0:
            print('# PROCESSED: {0} with {1}'.format(model, sfound))
            fullsed2 = _np.zeros((len(sed2data), 16))
            fullsed2[:, 0:3 + 1] = sed2data[:, 0:3 + 1]
            fullsed2[:, 4] = sed2data[:, 11]
            fullsed2[:, 5] = sed2data[:, 19]
            fullsed2[:, 6] = sed2data[:, 27]
            fullsed2[:, 7:8 + 1] = sed2data[:, 4:5 + 1]

            #a = _np.arange(9).reshape(3,3)
            #_np.core.records.fromarrays(a.transpose(), names='a,b,c', formats='f4,f4,f4')
            #fullsed2 = fullsed2[fullsed2[:,2].argsort()]
            #fullsed2 = fullsed2[fullsed2[:,0].argsort()]
            if _np.max(fullsed2[_np.isfinite(fullsed2)]) < 100000:
                fmt = '%13.6f'
            elif _np.max(fullsed2[_np.isfinite(fullsed2)]) < 1000000:
                fmt = '%13.5f'
            else:
                print('# ERROR at max values of fullsed2!!!!!!!!')
                raise SystemExit(0)
            fullsed2 = _np.core.records.fromarrays(fullsed2.transpose(), names= \
                'MU,PHI,LAMBDA,FLUX,SCT FLUX,EMIT FLUX,TRANS FLUX,Q,U,Sig FLUX,\
                Sig SCT FLUX,Sig EMIT FLUX,Sig TRANS FLU,Sig Q,Sig U',
                                                   formats='f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8')
            idx = _np.argsort(fullsed2, order=('MU', 'LAMBDA'))
            fullsed2 = fullsed2[idx]

            hd = '%CONTAINS: {0}\n'.format(' + '.join(sfound))
            hd += '%CREATED: {0}\n'.format(_time.asctime(_time.localtime(_time.time())))
            hd += '%{0:>7s}{1:>8s}{2:>13s}{3:>13s}'.format('nlbd', 'nobs', 'Rstar', 'Rwind') + '\n'
            hd += '{0:8d}{1:8d}{2:13.4f}{3:13.2f}\n'.format(int(nlbd), int(nobs), Rstar, Rwind)
            hd += '%{0:>12s}'.format('MU') + (15 * '{:>13s}').format('PHI', 'LAMBDA', 'FLUX', \
                                                                      'SCT FLUX', 'EMIT FLUX', 'TRANS FLUX', 'Q', 'U',
                                                                      'Sig FLUX', 'Sig FLUX', \
                                                                      'SigSCTFLX', 'SigEMITFLX', 'SigTRANSFLX',
                                                                      'Sig Q', 'Sig U')
            _np.savetxt(path + 'fullsed/fullsed_' + modelname.replace('.txt', '.sed2'), \
                        fullsed2, header=hd, comments="", fmt=fmt, delimiter='')
        else:
            print('# WARNING: No SED2 found for {0}'.format(model))
    return


def calcTeff(lum, size, M=None):
    """
    Calculate Teff for the non-rotating case.

    INPUT: Lum, Radius (or Mass in Solar units). `size` variable is assumed to
    be  the stellar radius (i.e., M==None). If M is given, size is assumed to
    be log(g) (cgs units).

    OUTPUT: float (Kelvin)
    """
    # ~ M, Rp, Lum = fundline
    if M == None:
        Rp = size * _phc.Rsun.cgs
    else:
        Rp = (M * _phc.Msun.cgs * _phc.G.cgs / 10 ** size) ** .5
    L = lum * _phc.Lsun.cgs
    #Lum = 4*_np.pi*Rp**2*_phc.sigma*Teff**4
    Teff = (L / (4 * _np.pi * Rp ** 2 * _phc.sigma.cgs)) ** .25
    return Teff


def calclogg(M, R):
    """
    Calculate log(g) for the non-rotating case. 

    INPUT: Radius and Mass in Solar units.

    OUTPUT: log(g) in cgs units.
    """
    R = R * _phc.Rsun.cgs
    logg = _np.log10(_phc.G.cgs * M * _phc.Msun.cgs / R ** 2)
    return logg


def genlog(path=None, extrainfo=None):
    """Gen. log of the calculated models of the project.

    ppath = Project's path. If it is not given, it assumes the local pwd.

    | extrainfo = {\
    | 'mod01':'i=60',\
    | 'mod02':'i=60+70',\
    | 'mod03':'i=60+70',\
    | 'mod04':'i=60+70/-source NO ROT',\
    | 'mod05':'i=60+70/?',\
    | 'mod06':'i=60+70/?'}

    INPUT: *path (string), *extrainfo (dictionary with modn number as index)

    OUTPUT: file written.
    """
    if path == None:
        path = _os.getcwd()
    modfld = _glob('{0}/mod*'.format(path))
    # while len(modfld) == 0:
    #    proj = raw_input('Type the project name: ')
    #    modfld = _glob('{0}/mod*'.format(proj))
    modfld.sort()

    #MODN, steps, sed2, maps, extrainfo
    tab = _np.zeros((5000, 5), dtype='|S127')
    i = 0

    for modn in modfld:
        modnn = _phc.trimpathname(modn)[1]
        modglob = _glob('{0}/*'.format(modn))
        mods = [x for x in modglob if (x.find('.txt') > -1 and x.find('{0}_'.\
        format(modnn)) >-1 )]
        print('# Catalogue of {0}'.format(modn))
        for mod in mods:
            suf = mod[mod.rfind('/') + 1:-4]

            #~ step1 = _glob('{0}/{1}??.temp'.format(modn, suf))
            sufglob = [x for x in modglob if (x.find(suf) > -1)]
            step1 = [x for x in sufglob if (x.find('.temp') > -1 and
            x.find('/{0}_'.format(modnn)) >-1 and x.find('_avg') == -1)]
            if len(step1) == 0:
                step1 = ['0']
            step1.sort()

            #~ sed2 = _glob('{0}/*{1}_*.sed2'.format(modn, suf))
            #~ sed2 += _glob('{0}/*{1}.sed2'.format(modn, suf))
            sed2 = [x for x in sufglob if x.find('.sed2') > -1]
            sed2.sort()
            s2out = ''
            for sedi in sed2:
                s2out += sedi[sedi.rfind('/') + 1:sedi.find('_')] + '+'

            #~ maps = _glob('{0}/*{1}_*.maps'.format(modn, suf))
            #~ maps += _glob('{0}/*{1}.maps'.format(modn, suf))
            maps = [x for x in sufglob if x.find('.map') > -1]
            maps.sort()
            mout = ''
            for mapi in maps:
                mout += mapi[mapi.rfind('/') + 1:mapi.find('_')] + '+'

            if extrainfo != None:
                if modnn in extrainfo:
                    extra = extrainfo[modnn]
            else:
                extra = ''

            tab[i] = (suf, step1[-1][-7:-5], s2out[:-1], mout[:-1], extra)
            i += 1

        if len(mods) == 0:
            print('# NO model found in {0}'.format(modn))

    #tab = tab[:i,:]
    #tab = tab[tab[:,0].argsort()]
    #_np.savetxt('log.csv', tab, fmt='%s', delimiter=',')
    _np.savetxt('{0}/log.csv'.format(path), tab[:i, :], fmt='%s', delimiter=',')
    #tab = _np.loadtxt('log.csv', dtype=str)
    #tab = tab[tab[:,0].argsort()]
    #_np.savetxt('log.csv', tab, fmt='%s', delimiter=',')
    print('# Generated {0}/log.csv !'.format(path))
    return


def printN0(modn):
    """
    Print the n0 of the model inside modn folder. It does a grep of `n_0` of
    all *.txt files of the folder.

    INPUT: string

    OUTPUT: figures written. 
    """
    import matplotlib.pyplot as plt

    # generate and print n_0
    _os.system('grep n_0 mod{0}/*.txt > n0_mod{0}.txt'.format(modn))
    #n0file = _np.genfromtxt('n0_mod01.txt', delimiter = (6,3,3,4,4,2,3,3,5,\
    #5,5,3,4), usecols=(4,10,12), dtype=None) 
    #cols = (0#, 1type, 2typeVal, 3#, 4sig, 5#, 6h, 7#, 8Rd, 9#, 10M, 11#,
    # 12ob, 13#)
    n0file = _np.genfromtxt('n0_mod{0}.txt'.format(modn), delimiter= \
        (22, 4, 18, 5, 3, 4, 35, 8), usecols=(1, 3, 5, 7), dtype=None)
    #cols = (sig0, M, ob, n0)

    sig0s = []
    for item in n0file[:, 0]:
        if item not in sig0s:
            sig0s += [item]

    Ms = []
    for item in n0file[:, 1]:
        if item not in Ms:
            Ms += [item]

    obs = [1.1, 1.2, 1.3, 1.4, 1.45]

    _plt.clf()
    for item in obs:
        idx = _np.where(n0file[:, 2] == item)
        _plt.title('Sig0 range: [{0}, {1}] (g cm-2)'.format(sig0s[0], sig0s[-1]))
        _plt.plot(n0file[:, 0][idx], n0file[:, 3][idx], 'o')
        _plt.yscale('log')
        _plt.xscale('log')
        _plt.ylabel('n0 (# cm-3)')
        _plt.xlabel('sig0 (g cm-2)')
    _plt.savefig('n0_vs_sig0.png')
    _plt.savefig('n0_vs_sig0.pdf')

    Btp = ['B0.5', 'B1', 'B1.5', 'B2', 'B2.5', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9']
    _plt.clf()
    for item in obs:
        idx = _np.where(n0file[:, 2] == item)
        _plt.title('Sig0 range: [{0}, {1}] (g cm-2)'.format(sig0s[0], sig0s[-1]))
        _plt.plot(n0file[:, 1][idx], n0file[:, 3][idx], 'o')
        _plt.plot([3.4, 14.6], [1e13, 1e14], 'r-')
        _plt.plot([3.4, 14.6], [1e12, 1e12], 'r-')
        _plt.xlim([16, 2])
        _plt.yscale('log')
        _plt.ylabel('n0 (# cm-3)')
        _plt.xlabel('mass (Solar units)')
    for i in range(len(Btp)):
        _plt.text(Ms[::-1][i], 4e14, Btp[i])
    _plt.savefig('n0_vs_M.png')
    _plt.savefig('n0_vs_M.pdf')

    _plt.clf()
    for item in obs:
        idx = _np.where(n0file[:, 2] == item)
        _plt.title('Sig0 range: [{0}, {1}] (g cm-2)'.format(sig0s[0], sig0s[-1]))
        _plt.plot(n0file[:, 2][idx], n0file[:, 3][idx], 'o')
        _plt.yscale('log')
        _plt.ylabel('n0 (# cm-3)')
        _plt.xlabel('radii ratio (Req/Rp)')
        _plt.xlim([1.0, 1.5])
    _plt.savefig('n0_vs_ob.png')
    _plt.savefig('n0_vs_ob.pdf')
    return


def diskcalcs(M, R, Tpole, T, alpha, R0, mu, rho0, Rd):
    """ Do the equivalence of disk density for different quantities.

    Note that they all depend of specific stellar quantities!!!

    INPUT: 
    M = 10.3065*Msun,
    R = 7*Rsun,
    Tpole = 26025.,
    T = 0.72*Tpole,
    alpha = 1.,
    R0 = 1e14*R,
    mu = 0.5,
    rho0 = 5e12 #in particles per cubic centimeter,
    Rd = 18.6*R.

    OUTPUT: printed status
    """

    def rho2sigp(R, rho0, a, M):
        sig = (2 * _np.pi) ** .5 * a / (_phc.G.cgs * M / R) ** .5 * R * rho0
        return sig

    def rho2Mdot(R, alpha, a, M, rho0, R0):
        Mdot = 3 * _np.pi * (2 * _np.pi) ** .5 * alpha * a ** 3. / (_phc.G.cgs * M / R) * rho0 * R ** 2. * ((
                                                                                                                R0 / R) ** .5 - 1) ** -1.
        return Mdot

    def Mdot2sig(R, Mdot, alpha, a, M, R0):
        sig = Mdot * (_phc.G.cgs * M / R) ** .5 / (3. * _np.pi * alpha * a ** 2 * R) * ((R0 / R) ** .5 - 1)
        return sig

    a = (_phc.kB.cgs * T / mu / _phc.mH.cgs) ** .5
    rho0 = rho0 * mu * _phc.mH.cgs
    sigp = rho2sigp(R, rho0, a, M)
    rho0p = sigp / (2 * _np.pi) ** .5 / a * (_phc.G.cgs * M / R) ** .5 / R
    Mdot = rho2Mdot(R, alpha, a, M, rho0, R0)
    sig = Mdot2sig(R, Mdot, alpha, a, M, R0)
    sigl = Mdot * (_phc.G.cgs * M) ** .5 / 3 / _np.pi
    Mdisk0 = (2 * _np.pi) ** 1.5 * rho0 * R ** 2. * (Rd - R) * a / (_phc.G.cgs * M / R) ** .5
    Mdisk = 2 * _np.pi * Mdot * (
                                    _phc.G.cgs * M / R) ** .5 * R ** .5 / 3 / _np.pi / alpha / a ** 2. * R0 ** .5 * _np.log(
        Rd / R)
    if Mdisk == 0:
        Mdisk = 4 * _np.pi * Mdot * (_phc.G.cgs * M / R) ** .5 * R ** .5 / 3 / _np.pi / alpha / a ** 2. * (
            Rd ** .5 - R ** .5)
    MdiskG = 2 * Mdot * (_phc.G.cgs * M / R) ** .5 * R ** .5 / 3 / alpha / a ** 2 * (
        R0 ** .5 * _np.log(Rd / R) + 2 * R ** .5 - 2 * Rd ** .5)

    print('R0/R  = {0:.1f}'.format(R0 / R))
    print('Valid sigma (1)?: {0}'.format(round(sigp / sig) == 1))
    print('Valid sigma (2)?: {0}'.format(round(rho0 / rho0p) == 1))
    print('rho0  = {0:.2e} g/cm3'.format(rho0))
    print('sigma0= {0:.2e} g/cm2'.format(sig / alpha / a ** 2))
    print('Mdot  = {0:.2e} Msun/yr'.format(Mdot / _phc.Msun.cgs * _phc.yr.cgs))
    print('Mdisk0= {0:.2e} Msun [#from rho0]'.format(Mdisk0 / _phc.Msun.cgs))
    print('Mdisk = {0:.2e} Msun [#approx.]'.format(Mdisk / _phc.Msun.cgs))
    print('MdiskG= {0:.2e} Msun'.format(MdiskG / _phc.Msun.cgs))
    print('PS: Mdisk for both isothermal Sigma(r) and H(r)')
    return


def rotStar(Tp=20000., M=10.3065, rp=5.38462, star='B', beta=0.25, wfrac=0.8,
            th_res=5001, quiet=False, LnotTp=False):
    """ Return the photospheric parameters of a rotating star.

    `LnotTp`: the value of "Tp" is the Luminosity (in solar units).

    Calculation of Von Zeipel's Beta parameter as function of W: see math...

    INPUT: th_res (theta resolution, integer)...

    OUTPUT: printed status + (ob, Tp values, Area[cm2])
    """
    Rsun = _phc.Rsun.cgs
    Msun = _phc.Msun.cgs
    Lsun = _phc.Lsun.cgs
    G = _phc.G.cgs
    AU = _phc.au.cgs
    pc = _phc.pc.cgs
    sigma = _phc.sigma.cgs
    M = M * Msun
    rp = rp * Rsun
    if wfrac == 0.:
        wfrac = 1e-9
    if LnotTp:
        Tp = (Tp*Lsun/4./_np.pi/rp**2/sigma)**.25

    ### DEFS ###
    def rt(th, wfrac):
        if th == 0:
            r = 1.
        else:
            r = (-3. * _np.cos((_np.arccos(wfrac * _np.sin(th)) + 4 * _np.pi) / 3)) / (wfrac * _np.sin(th))
        return r

    def area(wfrac):
        ths = _np.linspace(_np.pi / 2, 0, th_res)
        a = 0.
        for i in range(len(ths)):
            a = a + 2 * _np.pi * rt(ths[i], wfrac) ** 2 * _np.sin(ths[i])
        return 2 * a * ths[-2]

    def g(wfrac, M, rp, th):
        wcrit = _np.sqrt(8 * G * M / (27 * rp ** 3))
        g = (wcrit * wfrac) ** 2 * rp * rt(th, wfrac) * _np.sin(th) ** 2 - G * M / (rp * rt(th, wfrac)) ** 2
        return g

    def lum(wfrac, Tp, rp, M, C, beta):
        ths = _np.linspace(_np.pi / 2, 0, th_res)
        l = 0.
        for i in range(len(ths)):
            l = l + rt(ths[i], wfrac) ** 2 * _np.sin(ths[i]) * (abs(g(wfrac, M, rp, ths[i]))) ** (4 * beta)
        return 2 * 2 * _np.pi * ths[-2] * sigma * rp ** 2 * C ** (4 * beta) * l

    def lumf(wfrac, Tp, rp, M, beta):
        ths = _np.linspace(_np.pi / 2, 0, th_res)
        l = 0.
        for i in range(len(ths)):
            l = l + rt(ths[i], wfrac) ** 2 * _np.sin(ths[i]) * abs(g(wfrac, M, rp, ths[i])) ** (4 * beta)
        return l * ths[-2] * rp ** 2

    Bstars = _np.array(_phc.bestars, dtype=str)
    if star in Bstars:
        i = _np.where(Bstars[:, 0] == star)
        i = i[0][0]
        print Bstars[i][0]
        Tp = float(Bstars[i][1])
        M = float(Bstars[i][2]) * Msun
        rp = float(Bstars[i][3]) * Rsun
        # comentar linha abaixo se 1a. rodada:
        #Tp = 27438.63 #K

    wcrit = _np.sqrt(8 * G * M / (27 * rp ** 3))
    C = Tp ** (1. / beta) / abs(G * M / rp ** 2)

    vrot = wcrit * wfrac * rp * rt(_np.pi / 2, wfrac)
    lum0 = 4 * _np.pi * rp ** 2 * sigma * Tp ** 4 / Lsun

    ###a = rp**2*Tp**4*abs(g(wfrac,M,rp,0.))**(4*beta)
    ###print('Teff_pol* = %.2f' % ( (a/b)**beta ) )
    b = lumf(wfrac, Tp, rp, M, beta)
    c = lumf(0.0001, Tp, rp, M, beta)
    Cw = (c / b) ** (1. / (4. * beta)) * C
    ob = rt(_np.pi / 2, wfrac)#/(rp / Rsun)

    ### OUTPUT ###
    if not quiet:
        print('# Parameters:')
        print('wfrac     = %.4f' % (wfrac))
        print('W         = %.4f' % (_np.sqrt(2*(ob-1))))
        print('Star Mass = %.2f Msun' % (M / Msun))
        print('Rpole     = %.2f Rsun' % (rp / Rsun))
        print('Req       = %.2f Rpole' % (rt(_np.pi / 2, wfrac)))
        print('Teff_pol  = %.1f' % (Tp))
    
        print('Star Area = %.2f' % (area(wfrac)))
        print('Star Lum. = %.1f' % (lum(wfrac, Tp, rp, C, M, beta) / Lsun))
        print('Star Lum.*= %.1f' % (lum0))
    
        print('vrot(km/s)= %.1f' % (vrot / 1e5))
        print('vorb(km/s)= %.1f' % (_np.sqrt(G * M / rp / rt(_np.pi / 2, wfrac)) / 1e5) )
        print('vcrt(km/s)= %.1f' % (wcrit * rp * rt(_np.pi / 2, 1.) / 1e5))
    
        print('log(g)pole= %.2f' % (_np.log10(abs(g(wfrac, M, rp, 0.))) ))
        print('log(g)eq  = %.2f' % (_np.log10(abs(g(wfrac, M, rp, _np.pi / 2))) ))
        print('Teff_eq   = %.1f' % ( (C * abs(g(wfrac, M, rp, _np.pi / 2))) ** beta) )
        print('Teff_eq*  = %.1f' % ( (Cw * abs(g(wfrac, M, rp, _np.pi / 2))) ** beta) )
    
        print('Teff_pol* = %.2f' % ( (Cw * abs(g(wfrac, M, rp, 0.))) ** beta) )
        print(
            'T_pol/eq* = %.4f' % (
            (Cw * abs(g(wfrac, M, rp, 0.))) ** beta / (Cw * abs(g(wfrac, M, rp, _np.pi / 2))) ** beta) )
    
        print('# \"*\" == case where L is constant!')
    return ob, (Cw * abs(g(wfrac, M, rp, 0.))) ** beta, area(wfrac)*(rp**2)


def obsCalc():
    """ Obs. Calculation

    INPUT:

    OUTPUT:
    """

    def base60_to_decimal(xyz, delimiter=None):
        """Decimal value from numbers in sexagesimal system.

        The input value can be either a floating point number or a string
        such as "hh mm ss.ss" or "dd mm ss.ss". Delimiters other than " "
        can be specified using the keyword ``delimiter``.
        """
        divisors = [1, 60.0, 3600.0]
        xyzlist = str(xyz).split(delimiter)
        sign = -1 if xyzlist[0].find("-") != -1 else 1
        xyzlist = [abs(float(x)) for x in xyzlist]
        decimal_value = 0

        for i, j in zip(xyzlist, divisors):  # if xyzlist has <3 values then
            # divisors gets clipped.
            decimal_value += i / j

        decimal_value = -decimal_value if sign == -1 else decimal_value
        return decimal_value

    def julian_date(year, month, day, hour, minute, second):
        """Given year, month, day, hour, minute and second return JD.

        ``year``, ``month``, ``day``, ``hour`` and ``minute`` are integers,
        truncates fractional part; ``second`` is a floating point number.
        For BC year: use -(year-1). Example: 1 BC = 0, 1000 BC = -999.
        """
        MJD0 = 2400000.5  # 1858 November 17, 00:00:00 hours

        year, month, day, hour, minute = \
            int(year), int(month), int(day), int(hour), int(minute)

        if month <= 2:
            month += 12
            year -= 1

        modf = _math.modf
        # Julian calendar on or before 1582 October 4 and Gregorian calendar
        # afterwards.
        if ((10000L * year + 100L * month + day) <= 15821004L):
            b = -2 + int(modf((year + 4716) / 4)[1]) - 1179
        else:
            b = int(modf(year / 400)[1]) - int(modf(year / 100)[1]) + \
                int(modf(year / 4)[1])

        mjdmidnight = 365L * year - 679004L + b + int(30.6001 * (month + 1)) + day

        fracofday = base60_to_decimal( \
            " ".join([str(hour), str(minute), str(second)])) / 24.0

        return MJD0 + mjdmidnight + fracofday

    def obs_ver(hmin, hnas, hmax, hpoe):
        debug = False
        if hnas - hmin < -12:
            cor = 24
        elif hnas - hmin > 12:
            cor = -24
        else:
            cor = 0
        if hmin > hnas + cor and hpoe > hmin:
            hnas = hmin
            if debug:
                print 'ok0'
        elif hmin > hnas + cor and hmax > hpoe and (hpoe - hmin) < -12:
            hnas = hmin
            if debug:
                print 'ok0a'
        elif hmin > hnas + cor and hmax < hpoe and hnas < hmax:
            hnas = hmin
            if debug:
                print 'ok1'
        elif hmin < hnas + cor and hmax < hpoe and hnas < hmax:
            hpoe = hmax
            if debug:
                print 'ok2'
        elif hmin < hnas + cor and hmax < hpoe and (hnas - hmax) > 12:
            hpoe = hmax
            if debug:
                print 'ok2a'
        elif hmin < hnas + cor and hmax > hpoe and (hpoe - hmin) < -12:
            if debug:
                print 'ok3'
        elif hmin < hnas + cor and hmax > hpoe and hpoe > hmin:
            if debug:
                print 'ok4'
        else:
            if debug:
                print(hmin, hnas, hmax, hpoe)
            hnas = hpoe = _np.NaN
        return (hnas, hpoe)

    # equinocio set 2011 (djsol)
    djsol = 2455827.87778

    #Input: dia gregoriano (dg), para dia juliano (dj) as 0h local
    #+3 do UT #OPD@LNA
    ut = -3.

    print("# Programa de planejamento de alvos BEACON\n#")
    print("# Horas em UT=%d ! (horario de 'inverno' de Brasilia)" % (ut))
    print("# Digite a Data de Observacao, ou 'ENTER' para hoje...\n#")
    dg = []
    dg = dg + [str(raw_input("Digite o Ano (xxxx): "))]
    try:
        j = float(dg[-1])
        dg = dg + [str(raw_input("Digite o Mes (1-12): "))]
        dg = dg + [str(raw_input("Digite o Dia (1-31): "))]
    except:
        now = _datetime.datetime.now()
        dg = []
        dg = dg + [now.year]
        dg = dg + [now.month]
        dg = dg + [now.day]

    dg = _np.array(dg, dtype=float)
    dj = julian_date(dg[0], dg[1], dg[2] + 1, 0 + ut, 0, 0.)
    #dj = julian_date(now.year, now.month, now.day, now.hour, now.minute, now.second)

    #dias/ano passados do solsticio de verao (ndays)
    ndays = dj - djsol
    while ndays > 365:
        ndays = ndays - 365
    while ndays < 0:
        ndays = ndays + 365

    #Hora sideral 'as 0h local (hz) em horas
    hz = ndays * (0.065711111)

    #Reducao do tempo de observacao (rt) em minutos
    if ndays <= 91.25:
        rt = 75 + 15 / 91.25 * ndays
    elif ndays <= 273.75:
        rt = 90 - (ndays - 91.25) * 180 / 365
    else:
        rt = (ndays - 273.75) * 75 / 91.25
    #rt de min. para horas
    rt = rt / 60

    #hora maximo (hmax) e minimo (hmin) de observacao
    hmin = hz - 6 + rt / 2
    if hmin < 0:
        hmin = hmin + 24
    hmax = hz + 6 - rt / 2
    if hmax > 24:
        hmax = hmax - 24

    ##observacao minima e maxima em dj (dmin, dmax)
    ##rt de horas para segs.
    ##rt = rt*3600.
    ##dmin = julian_date(dg[0],dg[1],dg[2]+1,3-6,0,rt/2)
    ##dmax = julian_date(dg[0],dg[1],dg[2]+1,3+6,0,-rt/2) #seg. so' >0!!!

    #carrega lista de alvos
    alvos = _np.loadtxt('{0}refs/obs_alvos.txt'.format(hdtpath()), dtype=str, \
                        delimiter='\t')
    #carrega tempo das declinacoes
    obsdec = _np.loadtxt('{0}refs/obs_dec.txt'.format(hdtpath()), delimiter='\t')
    #carrega efemerides
    if _os.path.exists('{0}refs/obs_ef.txt'):
        ef_alvos = _np.loadtxt('{0}refs/obs_ef.txt'.format(hdtpath()), \
                               delimiter='\t', dtype=str)
        ef_alvos = ef_alvos.T

    print("\n# Info. da noite: %2d %2d %4d" % ( dg[2], dg[1], dg[0]) )
    print(  "Tempo sideral no anoitecer : %2d %2d" % ( int(hmin), int((hmin - int(hmin)) * 60) )  )
    print(  "Tempo sideral 'a meia noite: %2d %2d" % ( int(hz), int((hz - int(hz)) * 60) )  )
    print(  "Tempo sideral no amanhecer : %2d %2d" % ( int(hmax), int((hmax - int(hmax)) * 60) )  )
    print(  "Dif. de observ. ao inverno : -%d %2d" % ( int(rt), int((rt - int(rt)) * 60) )  )
    print("\nALVO\tINICIO\tFIM\tMERID.\tF_INI\tF_FIM")

    outfile = "ALVO\tINICIO\tFIM\tMERID.\tF_INI\tF_FIM\n"
    for i in range(len(alvos)):
        #calcula a ascencao reta (ra) e declinacao (dec)
        ra = float(alvos[i][2][:2]) + float(alvos[i][2][3:5]) / 60
        dec = float(alvos[i][3][:3]) + float(alvos[i][3][0] + alvos[i][3][4:6]) / 60
        #calcula tempo observavel na declinacao (to)
        j = odmin = odmax = 0
        while odmax == 0:
            decmin = obsdec[j][0]
            odmin = obsdec[j][1]
            if obsdec[j + 1][0] < dec:
                decmax = obsdec[j + 1][0]
                odmax = obsdec[j + 1][1]
            j = j + 1

        to = (abs(dec - decmax) * odmin + abs(dec - decmin) * odmax) / (abs(dec - decmin) + abs(dec - decmax))

        #calcula para cada alvo o tempo de observacao
        #hora de inicio da observacao ('nascer')
        hnas = ra - to / 2
        if hnas < 0:
            hnas = hnas + 24
        #hora de termino da observacao ('poente')
        hpoe = ra + to / 2
        if hpoe > 24:
            hpoe = hpoe - 24
        #observavel ao anoitecer?
        (hnas, hpoe) = obs_ver(hmin, hnas, hmax, hpoe)
        #estara' observavel por mais de um 3/4 de hora? (dt)
        dt = hpoe - hnas
        if dt < 0:
            dt = dt + 24
        if dt > 12:
            dt = dt - 24
        if dt < 3 / 4.:
            hnas = _np.NaN
            hpoe = _np.NaN

        #procura posicao nas efemerides (pef)
        if _os.path.exists('database/obs_ef.txt'):
            pef = [j for j, x in enumerate(ef_alvos[1]) if x.find(alvos[i][1]) > -1]
        else:
            pef = []
        if len(pef) == 1 and hnas > 0:
            pef = pef[0]
            per = float(ef_alvos[2][pef])
            T0 = float(ef_alvos[3][pef])
            #calcula fases
            #dj dia juliano 'a meia-noite local
            #dj = julian_date(dg[0],dg[1],dg[2]+1,0+ut,0,0.)
            #observacao ANTES da 0h UT
            cor = j = 0
            if (hnas - hz) > 12:
                cor = -24
            if hnas - hz + ut + cor < 0:  ##or (hz-hnas < 0 and hz-nas<-12):
                j = 1
                hour = hnas - hz + ut + cor
                minute = (hour - int(hour)) * 60
                hour = int(hour)
                djnas = julian_date(dg[0], dg[1], dg[2], hour, minute, 0.)
            else:
                j = 2
                hour = hnas - hz + ut + cor
                minute = (hour - int(hour)) * 60
                hour = int(hour)
                djnas = julian_date(dg[0], dg[1], dg[2] + 1, hour, minute, 0.)
            ##print(j,hour,djnas)
            j = 0
            if hpoe - hz + ut + cor < 0:
                j = 1
                hour = hpoe - hz + ut + cor
                minute = (hour - int(hour)) * 60
                hour = int(hour)
                djpoe = julian_date(dg[0], dg[1], dg[2], hour, minute, 0.)
            else:
                j = 2
                hour = hpoe - hz + ut + cor
                minute = (hour - int(hour)) * 60
                hour = int(hour)
                djpoe = julian_date(dg[0], dg[1], dg[2] + 1, hour, minute, 0.)
            ##print(j,hour,djpoe)
            #Warnings
            ##if T0 == 0.:
            ##  print("WARNING: No 'T0' is available for this star!!!")
            ##if per == 1.:
            ##  print("WARNING: No Period is available for this star!!!")
            #print("# A fase atual e' "+str( (j-T0)%Per/Per )+"\n" )
            fase_inc = ((djnas - T0) % per) / per
            fase_fim = ((djpoe - T0) % per) / per
        else:
            fase_inc = _np.NaN
            fase_fim = _np.NaN

        #tempos locais
        hlmin = hnas - hz
        if hlmin < 0:
            hlmin = hlmin + 24
        hlmax = hpoe - hz
        if hlmax < 0:
            hlmax = hlmax + 24
        hmer = ra - hz
        if hmer < 0:
            hmer = hmer + 24

        if hlmin > 0:
            print("%s\t%2d %2d\t%2d %2d\t%2d %2d\t%.2f\t%.2f" % (
                alvos[i][0][:7], int(hlmin), int((hlmin - int(hlmin)) * 60), int(hlmax), int((hlmax - int(hlmax)) * 60),
                int(hmer), int((hmer - int(hmer)) * 60), fase_inc, fase_fim) )
            outfile = outfile + ("%s\t%2d %2d\t%2d %2d\t%2d %2d\t%.2f\t%.2f" % (
                alvos[i][0][:7], int(hlmin), int((hlmin - int(hlmin)) * 60), int(hlmax), int((hlmax - int(hlmax)) * 60),
                int(hmer), int((hmer - int(hmer)) * 60), fase_inc, fase_fim) ) + '\n'
        else:
            print("%s\t-- --\t-- --\t-- --\t-.--\t-.--" % (alvos[i][0][:7]) )
            outfile = outfile + ("%s\t-- --\t-- --\t-- --\t-.--\t-.--" % (alvos[i][0][:7]) ) + '\n'

    wfile = raw_input("\nDeseja salvar a lista?(Sim/outro): ")
    if wfile in ['s', 'sim', 'Sim', 'S', 'y', 'yes', 'Yes', 'Y']:
        f0 = open('obs_%2d_%2d_%4d.txt' % (dg[2], dg[1], dg[0]), 'w')
        f0.writelines(outfile)
        f0.close()
    return


def readSingleBe(sBfile):
    """ Read the singleBe output

    OUTPUT = rgrid, lsig_r, nsnaps, simdays, alpha
    """

    def readSBeBlock(lines):
        """ """
        tau = _np.array(lines[0]).astype(float)              # tauintval in rad
        tausec = _np.array(lines[1]).astype(float)           # tausec in rad
        sinject = _np.array(lines[2]).astype(float)          # `sinject` ?
        alpha_r =  _np.array(lines[3].split()).astype(float) # alpha(r)
        s1_r =  _np.array(lines[4].split()).astype(float)    # s1(r) = sig/sig0 ?
        sig_r =  _np.array(lines[5].split()).astype(float)   # sigma(r)
        maxr =  _np.array(lines[6]).astype(float)            # maxr = maximmum
                                                            #non-zero cell
        vr_cs  =  _np.array(lines[7].split()).astype(float)  # vel_rad/cs ?
        decrr =  _np.array(lines[8].split()).astype(float)   # Decretion rate
                                                            #(units?)
        return

    hs = 15                     # header size
    f0 = open(sBfile).read().split('\n')
    f0 = f0[:-1]
    nsnaps = (len(f0)-hs+1)/9   # number of snapshots
    line0 = f0[0].split()
    
    alpha=float(line0[0])       # constant alpha parameter
    teff=float(line0[1])        # stellar effective temperature in K
    tdisk=float(line0[3])       # disk temperature in K
    cs=float(line0[4])          # "sound speed" in cm/s
    mstar=float(line0[5])       # mass of the star, in solar masses
    req=float(line0[6])         # equatorial radius, in solar units
    omega0=float(line0[8])      # disk angular velocity at equator in rad/s?
    rho0=float(line0[13])       # g cm-3
    sigma0=float(line0[14])     # g cm-2
    rin=float(line0[15])        # internal radius of the disk (in req?)
    rout=float(line0[16])       # external radius of the disk (in req?)
    rinject=float(line0[17])    # radius of injection in the disk (in req?)
    n=int(line0[18])            # number of radial cells
    kinject=int(line0[19])      # cell of mass injection
    dt=float(line0[24])         # ?
    tauintval=float(line0[26])  # time steps (in rad)
    
    #BLOCKS
    ltau = _np.array(f0[hs+0::9]).astype(float)     # tauintval in rad
    ltausec = _np.array(f0[hs+1::9]).astype(float)  # tausec in seconds
    lsinject = _np.array(f0[hs+2::9]).astype(float) # `sinject` ?
    # alpha(r)
    lalpha_r =  _np.array([l.split() for l in f0[hs+3::9]]).astype(float) 
    # s1(r) = sig/sig0 ?
    ls1_r =  _np.array([l.split() for l in f0[hs+4::9]]).astype(float)
    # sigma(r)
    lsig_r =  _np.array([l.split() for l in f0[hs+5::9]]).astype(float)
    # maxr = maximmum non-zero cell
    lmaxr =  _np.array(f0[hs+6::9]).astype(float)            
    # vel_rad/cs ?  ##VARIABLE SIZE = not read
    #~ lvr_cs  =  _np.array([l.split() for l in f0[hs+7::9]]).astype(float)  
    # Decretion rate (units?)  ##VARIABLE SIZE = not read
    #~ ldecrr =  _np.array([l.split() for l in f0[hs+8::9]]).astype(float)   
    
    tauintvaldays = tauintval/omega0/24/3600       # time steps in days
    rgrid = _np.array(f0[4].split()).astype(float)  # radial grid values
    simdays = ltausec[-1]/24/3600                  # simulation total time in days
    
    return rgrid, lsig_r, nsnaps, simdays, alpha


### MAIN ###
if __name__ == "__main__":
    pass
