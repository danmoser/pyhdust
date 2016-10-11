# -*- coding:utf-8 -*-

"""PyHdust main module: Hdust tools.

This module contains:

- PyHdust package routines
- Hdust I/O functions
- Hdust useful routines and plots
- Be stars quantities conversions
- Astronomic useful and plotting functions
- Astronomic filters tools

:co-author: Despina Panoglou
:license: GNU GPL v3.0 https://github.com/danmoserp/yhdust/blob/master/LICENSE
"""
from __future__ import print_function
import os as _os
import time as _time
import datetime as _dt
import math as _math
import numpy as _np
import struct as _struct
from glob import glob as _glob
from itertools import product as _itprod
import pyhdust.phc as _phc
import pyhdust.jdcal as _jdcal
from pyhdust.tabulate import tabulate as _tab
from six import string_types as _strtypes
import warnings as _warn

try:
    import matplotlib.pyplot as _plt
    import matplotlib.patches as _mpatches
    from scipy import interpolate as _interpolate
    import pyfits as _pf
except ImportError:
    _warn.warn('# matplotlib, pyfits, six and/or scipy module not installed!!')

__version__ = '1.0.8'
__release__ = "Stable"
__author__ = "Daniel Moser"
__email__ = "dmfaes@gmail.com"


# Package tools
def hdtpath():
    """ 
    :rtype: str
    :returns: The module path.

    :Example:

    >>> hdt.hdtpath()
    /home/user/Scripts/pyhdust/
    """
    fulldir = _os.path.split(__file__)[0] + _os.path.sep
    # fulldir = _os.path.split(fulldir)[0] + _os.path.sep
    end = 'pyhdust' + _os.path.sep
    if not fulldir.endswith(end):
        fulldir += end
    return fulldir


# Hdust I/O
def sed2info(sfile):
    """ Read info from a HDUST SED2 file.

    :type  sfile: str
    :param sfile: HDUST SED2 file path

    :rtype: tuple of floats
    :returns: ``nlbd, nobs, Rstar, Rwind``, the HDUST parameters
    """
    f0 = open(sfile, 'r')
    fcont = f0.readlines()
    f0.close()
    info = ''
    i = 0
    while info is '':
        if fcont[i][0] != '%':
            info = _np.array(fcont[i].split(), dtype=float)
        i += 1
    return info


def readsed2(sfile):
    """ Read data from HDUST SED2 file.

    Note: this format is different of **readfullsed2**.

    :type sfile: str
    :param sfile: HDUST SED2 file path

    :rtype: ``np.array(nobs*nlbd, -1)``, where number of columns from SED2 file 
    replaces "-1"
    :returns: SED2 header info
    """
    sed2data = _np.loadtxt(sfile, skiprows=1)
    return sed2data


def readfullsed2(sfile):
    """ Read data from HDUST fullSED2 file.

    :type sfile: str
    :param sfile: HDUST fullSED2 file path

    :rtype: ``np.array(nobs, nlbd, -1)``, where the number of columns from 
        fullSED2 file replaces "-1"
    :returns: HDUST fullSED2 file content
    """
    sed2data = _np.loadtxt(sfile, skiprows=5)
    nlbd, nobs, Rstar, Rwind = sed2info(sfile)
    # print(_np.shape(sed2data), int(nobs), int(nlbd))
    sed2data = sed2data.reshape((int(nobs), int(nlbd), -1))
    for i, j in _itprod([3, 7], range(int(nobs))):
        isnan = _np.isnan(sed2data[j, :, i])
        if any(isnan):
            sed2data[j, :, i][isnan] = _np.interp(sed2data[0, :, 2][isnan], 
                sed2data[0, :, 2][~isnan], sed2data[j, :, i][~isnan])
    return sed2data


def readtemp(tfile, quiet=False):
    """ Read HDUST temp file

    - ncr = número de células da simulação na coordenada radial
    - ncmu = número de células da simulação na coordenada latitudinal
    - ncphi = número de células da simulação na coordenada azimutal
    - nLTE = number of atomic LTE levels
    - nNLTE = number of atomic NLTE levels
    - beta = flare disk parameter
    - Rstar = raio da estrela (Rsun)
    - Ra = raio máximo da simulação, onde acabam todas as poeiras (Rsun)
    - pcrc = coordenadas das células em raio
    - pcmuc = coordenadas das células em mu
    - pcphic = coordenadas das células em phi
    - pcr = distância entre as células em raio
    - pcmu = distância entre as células em mu
    - pcphi = distância entre as células em phi

    `data` format is: `data[nLTE+6, ncr, ncmu, ncphi]`

    .. seealso::

        Temperature are in `data[3, ...]`. More info. see :py:func:`plottemp` 
        function.

    :type sfile: str
    :param sfile: HDUST fullSED2 file path.

    :rtype: ``np.array(nobs, nlbd, -1)``, number of columns from fullSED2 file 
        replaces "-1".)
    :returns: HDUST fullSED2 file content.

    OUTPUT = ncr,ncmu,ncphi,nLTE,nNLTE,Rstar,Ra,beta,data,pcr,pcmu,pcphi
    """
    f = open(tfile, 'rb').read()
    ixdr = 0
    ncr, ncmu, ncphi, nLTE, nNLTE = _struct.unpack('>5l', f[ixdr:ixdr + 4 * 5])
    ixdr += 4 * 5
    Rstar, Ra, beta = _struct.unpack('>3f', f[ixdr:ixdr + 4 * 3])
    ixdr += 4 * 3
    # 
    rlen = (nLTE + 6) * ncr * ncmu * ncphi
    data = _struct.unpack('>{0}f'.format(rlen), f[ixdr:ixdr + 4 * rlen])
    ixdr += 4 * rlen
    data = _np.reshape(data, (nLTE + 6, ncr, ncmu, ncphi), order='F')
    #
    # this will check if the XDR is finished.
    if ixdr == len(f):
        if not quiet:
            print('# XDR {0} completely read!'.format(tfile))
    else:
        _warn.warn( '# XDR {0} not completely read!\n# length difference is '
            '{1}'.format(tfile, (len(f)-ixdr)/4 ) )
    #
    pcrc = data[0, :, 0, 0]
    pcr = _np.zeros(ncr + 1)
    pcr[0] = Rstar
    for icr in range(1, ncr + 1):
        pcr[icr] = pcr[icr - 1] + 2 * (pcrc[icr - 1] - pcr[icr - 1])
    pcmu = _np.zeros((ncmu + 1, ncr))
    # 
    for icr in range(0, ncr):
        pcmuc = data[1, icr, :, 0]
        pcmu[0, icr] = -1.
        for icmu in range(1, ncmu + 1):
            pcmu[icmu, icr] = pcmu[icmu - 1, icr] + 2. * \
                (pcmuc[icmu - 1] - pcmu[icmu - 1, icr])
    #
    pcphic = data[2, 0, 0, :]
    pcphi = _np.zeros(ncphi + 1)
    pcphi[0] = 0.
    for icphi in range(1, ncphi + 1):
        pcphi[icphi] = pcphi[icphi - 1] + 2 * \
            (pcphic[icphi - 1] - pcphi[icphi - 1])
    # 
    return ncr, ncmu, ncphi, nLTE, nNLTE, Rstar, Ra, beta, data, pcr, pcmu,\
        pcphi


def readdust(dfile):
    """ Read *.dust files

    - ncr = número de células da simulação na coordenada radial
    - ncmu = número de células da simulação na coordenada latitudinal
    - ncphi = número de células da simulação na coordenada azimutal
    - Rstar = raio da estrela (Rsun)
    - Ra = raio máximo da simulação, onde acabam todas as poeiras (Rsun)
    - pcrc = coordenadas das células em raio
    - pcmuc = coordenadas das células em mu
    - pcphic = coordenadas das células em phi
    - pcr = distância entre as células em raio
    - pcmu = distância entre as células em mu
    - pcphi = distância entre as células em phi

    - ntip = número de tipos de poeira determinados para a simulação 
    (composição)
    - na = número de tamanho por tipo da poeira 
    - NdustShells = número de camadas de poeiras da simulação
    - Rdust = raio onde está(ão) a(s) camada(s)
    - Tdestruction = temperatura na qual os grãos são evaporados
    - Tdust = temperatura da poeira numa da posição da grade da simulação 
    (r,phi,mu)
    - lacentro = controla os tipos e tamanhos das poeiras
    """
    f0 = open(dfile).read()
    f0 = f0.split('\n')
    # Header
    tmp = f0[0].split()
    ncr, ncmu, ncphi, NdustShells = _np.array(tmp[:4], dtype='int')
    Rstar, Ra = _np.array(tmp[4:], dtype='float')
    tmp = f0[1].split()
    ncrdust = int(tmp[0])
    Rdust, Tdestruction = _np.array(tmp[1:], dtype='float')
    ntip = int(f0[2])
    # 
    pcrc = _np.zeros(ncr)
    pcmuc = _np.zeros((ncmu, ncr))
    pcphic = _np.zeros(ncphi)
    i = 3
    for stip in range(ntip):
        na = int(f0[i])
        # 
        if stip is 0:
            Tdust = _np.zeros((ntip, na, ncr, ncmu, ncphi))
            lacentro = _np.zeros((ntip, na))
        i += 1
        lacentro[stip] = _np.array(f0[i])
        for icphi, icmu, icr in _itprod(range(ncphi), range(ncmu), range(ncr)):
            i += 1 
            tmp = _np.array(f0[i])
            pcrc[icr] = tmp[0]
            pcmuc[icmu, icr] = tmp[1]
            pcphic[icphi] = tmp[2]
            Tdust[stip, :, icr, icmu, icphi] = tmp[3:3 + na]
    # 
    pcr = _np.zeros(ncr + 1)
    pcr[0] = Rdust[0]
    for icr in range(1, ncr + 1):
        pcr[icr] = pcr[icr - 1] + 2 * (pcrc[icr - 1] - pcr[icr - 1])
    # 
    pcmu = _np.zeros((ncmu + 1, ncr))
    for icr in range(ncr):
        pcmu[0, icr] = -1.
        for icmu in range(1, ncmu + 1):
            pcmu[icmu, icr] = pcmu[icmu - 1, icr] + 2. * (pcmuc[icmu - 1] -
                pcmu[icmu - 1, icr])
    # 
    pcphi = _np.zeros(ncphi + 1)
    pcphi[0] = 0.
    for icphi in range(1, ncphi + 1):
        pcphi[icphi] = pcphi[icphi - 1] + 2 * \
            (pcphic[icphi - 1] - pcphi[icphi - 1])
    return ncr, ncmu, ncphi, ntip, na, Rstar, Ra, NdustShells, Rdust,\
        Tdestruction, Tdust, pcrc, pcmuc, pcphic, pcr, pcmu, pcphi, lacentro


def mergesed2(models, Vrots, path=None, checklineval=False, onlyfilters=None):
    """
    Merge all mod#/*.sed2 files into the fullsed file.

    It will check if all sed2info are the same (i.e., nobs, Rstar, Rwind). 
    That's the reason why SED was the first band in previous versions.
    If not, it ask if you want to continue (and receive and error).

    The presence of the SED file is not required anymore.
    The structure is set by the first broadband sed2 found.

    NO AVERAGE is coded (yet).

    `onlyfitlters` is a iterable of desired filters (``None`` for all).

    IMPORTANT: the line rest wavelength is assumed to be the center of the SEI 
    band!

    INPUT: models lists (*.txt or *.inp), Vrots (array).

    OUTPUT: *files written (status printed).
    """
    if isinstance(models, _strtypes):
        models = [models]
    if isinstance(Vrots, (int, long, float)):
        Vrots = [Vrots]
        if len(Vrots) < len(models):
            Vrots = list(Vrots)+(len(models)-len(Vrots))*Vrots[-1:]

    for i in range(len(models)):
        models[i] = models[i].replace('.inp', '.txt')

    for model in models:
        modfld, modelname = _phc.trimpathname(model)
        path = _phc.trimpathname(modfld[:-1])[0]
        if not _os.path.exists('{0}fullsed'.format(path)):
            _os.system('mkdir {0}fullsed'.format(path))
        sed2data = _np.empty(0)
        sfound = []
        # Get all *.sed2 and choose if it is a broad-band or a line
        if onlyfilters is None:
            lsed2 = _glob('{0}*{1}.sed2'.format(modfld, modelname[:-4]))
            lsed2.extend(_glob('{0}*{1}_SEI.sed2'.format(modfld, 
                modelname[:-4])))
        else:
            lsed2 = []
            for f in onlyfilters:
                pattern = _os.path.join(modfld, '{0}*{1}.sed2'.format(f, 
                    modelname[:-4]))
                lsed2.extend(_glob(pattern))
                pattern = _os.path.join(modfld, '{0}*{1}_SEI.sed2'.format(f, 
                    modelname[:-4]))
                lsed2.extend(_glob(pattern))
        # print lsed2, '{0}*{1}*.sed2'.format(modfld, modelname[:-4])
        for file in lsed2:
            suf = _phc.trimpathname(file)[-1].split('_')[0]
            sfound += [suf]
            newdata = readsed2(file)
            if file.find('_SEI.') == -1 or file.find('SED_') > -1:
                # Process broad-band
                if len(sed2data) == 0:
                    sed2data = newdata.copy()
                    nlbd, nobs, Rstar, Rwind = sed2info(file)
                else:
                    # Check if the SED2 file has the info as the first file
                    if _np.product((nobs, Rstar, Rwind) == sed2info(file)[1:])\
                        == 0:
                        key = ''
                        while key.upper() != 'Y':
                            _warn.warn('# {0} has different HDUST '
                                'output!'.format(modelname))
                            key = _phc.user_input(
                                'Do you want do proceed? (y/other): ')
                    nlbd += sed2info(file)[0]
                    sed2data = _np.vstack((sed2data, newdata))
            else:
                # Process lines
                Vrot = Vrots[models.index(model)]
                nlbdSEI, tmp, tmp, tmp = sed2info(file)
                # Rest wavelength
                # print nlbd, nlbdSEI, newdata[:nlbdSEI, 2]
                lbrest = (newdata[int(nlbdSEI)-1, 2] - newdata[0, 2]) / 2. +\
                    newdata[0, 2]
                if checklineval:
                    print('# I found {0} as wavelength for {1}'.format(lbrest,
                        suf))
                    print('# Type a number to change it (no sci. notation/' +
                        'or empty to continue): ')
                    loop = True
                    while loop:
                        userinp = _phc.user_input('# lbd: ')
                        if len(userinp) == 0:
                            loop = False
                        elif userinp.replace('.', '', 1).isdigit():
                            lbrest = float(userinp)
                            loop = False
                deltalbd = Vrot / _phc.c.cgs / 1e-5 * lbrest
                mini = newdata[0, 2] + deltalbd
                maxi = newdata[-1, 2] - deltalbd
                # print deltalbd, Vrot
                idx = _np.where(
                    (newdata[:, 2] >= mini) & (newdata[:, 2] <= maxi))
                # print len(newdata), len(newdata[idx])
                ncut = len(newdata) - len(newdata[idx])
                newdata = newdata[idx]
                if len(sed2data) == 0:
                    sed2data = newdata.copy()
                    nlbd, nobs, Rstar, Rwind = sed2info(file)
                else:
                    # Check if the SED2 file has the info as the first file
                    if _np.product((nobs, Rstar, Rwind) == sed2info(file)[1:])\
                        == 0:
                        key = ''
                        while key.upper() != 'Y':
                            _warn.warn('# {0} has different HDUST '
                                'input!!'.format(modelname))
                            key = _phc.user_input(
                                'Do you want do proceed? (y/other): ')
                    idx = _np.where(
                        (sed2data[:, 2] < mini) | (sed2data[:, 2] > maxi))
                    ncut += len(sed2data) - len(sed2data[idx])
                    sed2data = sed2data[idx]
                    nlbd += sed2info(file)[0] - ncut / sed2info(file)[1]
                    sed2data = _np.vstack((sed2data, newdata))
        # 
        if len(sfound) > 0:
            print('# PROCESSED: {0} with {1}'.format(model, sfound))
            fullsed2 = _np.zeros((len(sed2data), 16))
            fullsed2[:, 0:3 + 1] = sed2data[:, 0:3 + 1]
            fullsed2[:, 4] = sed2data[:, 11]
            fullsed2[:, 5] = sed2data[:, 19]
            fullsed2[:, 6] = sed2data[:, 27]
            fullsed2[:, 7:8 + 1] = sed2data[:, 4:5 + 1]
            # Old 
            # if _np.max(fullsed2[_np.isfinite(fullsed2)]) < 100000:
                # fmt = '%13.6f'
            # elif _np.max(fullsed2[_np.isfinite(fullsed2)]) < 1000000:
                # fmt = '%13.5f'
            # else:
                # print('# ERROR at max values of fullsed2!!!!!!!!')
                # print(_np.max(fullsed2[_np.isfinite(fullsed2)]))
                # raise SystemExit(0)
            fullsed2 = _np.core.records.fromarrays(fullsed2.transpose(), 
                names='MU,PHI,LAMBDA,FLUX,SCT FLUX,EMIT FLUX,TRANS FLUX,Q,U,' +
                'Sig FLUX,Sig SCT FLUX,Sig EMIT FLUX,Sig TRANS FLU,Sig Q,' +
                'Sig U', formats='f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,' +
                'f8,f8')
            idx = _np.argsort(fullsed2, order=('MU', 'PHI', 'LAMBDA'))
            fullsed2 = fullsed2[idx]
            hd = '%CONTAINS: {0}\n'.format(' + '.join(sfound))
            hd += '%CREATED: {0}\n'.format(
                _time.asctime(_time.localtime(_time.time())))
            hd += '%{0:>7s}{1:>8s}{2:>13s}{3:>13s}\n'.format(
                'nlbd', 'nobs', 'Rstar', 'Rwind')
            hd += '{0:8d}{1:8d}{2:13.4f}{3:13.2f}\n'.format(
                int(nlbd), int(nobs), Rstar, Rwind)
            # hd += '%{0:>12s}'.format('MU') + (15 * '{:>13s}').format('PHI', 
                # 'LAMBDA', 'FLUX', \
                # 'SCT FLUX', 'EMIT FLUX', 'TRANS FLUX', 'Q', 'U',
                # 'Sig FLUX', 'Sig FLUX', \
                # 'SigSCTFLX', 'SigEMITFLX', 'SigTRANSFLX',
                # 'Sig Q', 'Sig U')
            # _np.savetxt(path + 'fullsed/fullsed_' + modelname.replace('.txt', 
                # '.sed2'), fullsed2, header=hd, comments="", fmt=fmt, 
                # delimiter='')
            # New
            header = ['% PHI', 'LAMBDA', 'FLUX', 'SCT FLUX', 'EMIT FLUX', 
                'TRANS FLUX', 'Q', 'U', 'Sig FLUX', 'Sig FLUX', 'SigSCTFLX', 
                'SigEMITFLX', 'SigTRANSFLX', 'Sig Q', 'Sig U']
            outfile = hd + _tab(fullsed2, headers=header, tablefmt="plain")
            f0 = open(
                path + 'fullsed/fullsed_' + modelname.replace('.txt', '.sed2'),
                'w')
            f0.writelines(outfile)
            f0.close()
        else:
            _warn.warn('# No SED2 found for {0}'.format(model))
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
        alpha_r = _np.array(lines[3].split()).astype(float)  # alpha(r)
        s1_r = _np.array(lines[4].split()).astype(
            float)    # s1(r) = sig/sig0 ?
        sig_r = _np.array(lines[5].split()).astype(float)   # sigma(r)
        maxr = _np.array(lines[6]).astype(float)            # maxr = maximmum
                                                            # non-zero cell
        vr_cs = _np.array(lines[7].split()).astype(float)  # vel_rad/cs ?
        decrr = _np.array(lines[8].split()).astype(float)   # Decretion rate
                                                            #(units?)
        return

    hs = 15                     # header size
    f0 = open(sBfile).read().split('\n')
    f0 = f0[:-1]
    nsnaps = (len(f0) - hs + 1) / 9   # number of snapshots
    line0 = f0[0].split()

    alpha = float(line0[0])       # constant alpha parameter
    teff = float(line0[1])        # stellar effective temperature in K
    tdisk = float(line0[3])       # disk temperature in K
    cs = float(line0[4])          # "sound speed" in cm/s
    mstar = float(line0[5])       # mass of the star, in solar masses
    req = float(line0[6])         # equatorial radius, in solar units
    omega0 = float(line0[8])      # disk angular velocity at equator in rad/s?
    rho0 = float(line0[13])       # g cm-3
    sigma0 = float(line0[14])     # g cm-2
    rin = float(line0[15])        # internal radius of the disk (in req?)
    rout = float(line0[16])       # external radius of the disk (in req?)
    rinject = float(line0[17])    # radius of injection in the disk (in req?)
    n = int(line0[18])            # number of radial cells
    kinject = int(line0[19])      # cell of mass injection
    dt = float(line0[24])         # ?
    tauintval = float(line0[26])  # time steps (in rad)

    # BLOCKS
    ltau = _np.array(f0[hs + 0::9]).astype(float)     # tauintval in rad
    ltausec = _np.array(f0[hs + 1::9]).astype(float)  # tausec in seconds
    lsinject = _np.array(f0[hs + 2::9]).astype(float)  # `sinject` ?
    # alpha(r)
    lalpha_r = _np.array([l.split() for l in f0[hs + 3::9]]).astype(float) 
    # s1(r) = sig/sig0 ?
    ls1_r = _np.array([l.split() for l in f0[hs + 4::9]]).astype(float)
    # sigma(r)
    lsig_r = _np.array([l.split() for l in f0[hs + 5::9]]).astype(float)
    # maxr = maximmum non-zero cell
    lmaxr = _np.array(f0[hs + 6::9]).astype(float)            
    # vel_rad/cs ?  ##VARIABLE SIZE = not read
    # lvr_cs  =  _np.array([l.split() for l in f0[hs+7::9]]).astype(float)  
    # Decretion rate (units?)  ##VARIABLE SIZE = not read
    # ldecrr =  _np.array([l.split() for l in f0[hs+8::9]]).astype(float)   

    tauintvaldays = tauintval / omega0 / 24 / 3600  # time steps in days
    rgrid = _np.array(f0[4].split()).astype(float)  # radial grid values
    # simulation total time in days
    simdays = ltausec[-1] / 24 / 3600

    return rgrid, lsig_r, nsnaps, simdays, alpha


# Hdust useful
def plotdust(tfile):
    """ TBD!!

    .. seealso:: py:func:`readdust`.

    """
    return


def gentemplist(tfile, tfrange=None, avg=True):
    """ Generate a list of *.temp files.

    `tfile` = file name or file prefix to be plotted. If `tfrange` (e.g.,
    tfrange=[20,24] is present, it you automatically plot the interval.

    `avg` = include tfile_avg.temp file (if exists).

    OUTPUT = file list, label list
    """
    if tfrange is None:
        tfrange = [int(tfile.replace('.temp', '')[-2:])]
    ltfile = []
    ltlabel = []
    for tn in range(int(tfrange[0]), int(tfrange[-1]) + 1):
        if tfile.find('.temp') > 0:
            ltfile += [tfile.replace(tfile[-7:-5], '{0:02d}'.format(tn))]
        else:
            ltfile += [tfile + '{0:02d}.temp'.format(tn)]
        ltlabel += [str(tn)]
    tfile = ltfile[-1].replace('.temp', '_avg.temp')
    if _os.path.exists(tfile) and avg:
        ltfile += [tfile]
        ltlabel += ['Avg.']
    return ltfile, ltlabel


def genlog(proj=None, mods=None):
    """ Generate a log of the runned simulations. 
    """
    if proj is None:
        proj = _os.getcwd()
    modfld = mods
    if mods is None:
        modfld = [f for f in _os.listdir(proj) if (f.startswith('mod') and 
            _os.path.isdir(f))]
    modfld.sort()

    basedir = _os.getcwd()
    # name, temp, sed2, map
    allmods = [[], [], [], []]
    for modn in modfld:
        _os.chdir(_os.path.join(proj, modn))
        modn_list = _glob(modn+'*.txt')
        allmods[0].extend(modn_list)
        for imodn in modn_list:
            suf = _os.path.splitext(imodn)[0]
            #
            tmp = sorted(_glob(suf+'[0-9][0-9].temp'))
            if len(tmp) > 0:
                allmods[1].append( _os.path.splitext(tmp[-1])[0][-2:] )
            else:
                allmods[1].append( '' )
            #
            tmp = [s2[:s2.find(suf)] for s2 in 
                sorted(_glob('*'+suf+'*.sed2'))]
            if len(tmp) > 0:
                allmods[2].append( ' '.join(tmp) )
            else: 
                allmods[2].append( '' )
            #
            tmp = [mp[:mp.find(suf)] for mp in 
                sorted(_glob('*'+suf+'*.map*'))]
            if len(tmp) > 0:
                allmods[3].append( ' '.join(tmp) )
            else:
                allmods[3].append( '' )
        _os.chdir(basedir)

    glout = 'genlog_{0}.txt'.format(_time.strftime('%y%m%d'))
    f0 = open(_os.path.join(proj, glout), 'w')
    f0.writelines(_tab(_np.array(allmods).T, tablefmt="tsv"))
    print('# Generated {0} !'.format(_os.path.join(proj, glout)))
    f0.close()
    return _np.array(allmods)


def printN0(modn):
    """ Print the n0 of the model inside modn folder. It does a grep of `n_0` of
    all *.txt files of the folder.

    INPUT: string

    OUTPUT: figures written. 
    """
    # generate and print n_0
    _os.system('grep n_0 mod{0}/*.txt > n0_mod{0}.txt'.format(modn))
    # n0file = _np.genfromtxt('n0_mod01.txt', delimiter = (6,3,3,4,4,2,3,3,5,\
    # 5,5,3,4), usecols=(4,10,12), dtype=None) 
    # cols = (0#, 1type, 2typeVal, 3#, 4sig, 5#, 6h, 7#, 8Rd, 9#, 10M, 11#,
    # 12ob, 13#)
    n0file = _np.genfromtxt('n0_mod{0}.txt'.format(modn), delimiter=(
        22, 4, 18, 5, 3, 4, 35, 8), usecols=(1, 3, 5, 7), dtype=None)
    # cols = (sig0, M, ob, n0)

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
        _plt.title(
            'Sig0 range: [{0}, {1}] (g cm-2)'.format(sig0s[0], sig0s[-1]))
        _plt.plot(n0file[:, 0][idx], n0file[:, 3][idx], 'o')
        _plt.yscale('log')
        _plt.xscale('log')
        _plt.ylabel('n0 (# cm-3)')
        _plt.xlabel('sig0 (g cm-2)')
    _plt.savefig('n0_vs_sig0.png')
    _plt.savefig('n0_vs_sig0.pdf')

    Btp = ['B0.5', 'B1', 'B1.5', 'B2', 'B2.5',
        'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9']
    _plt.clf()
    for item in obs:
        idx = _np.where(n0file[:, 2] == item)
        _plt.title(
            'Sig0 range: [{0}, {1}] (g cm-2)'.format(sig0s[0], sig0s[-1]))
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
        _plt.title(
            'Sig0 range: [{0}, {1}] (g cm-2)'.format(sig0s[0], sig0s[-1]))
        _plt.plot(n0file[:, 2][idx], n0file[:, 3][idx], 'o')
        _plt.yscale('log')
        _plt.ylabel('n0 (# cm-3)')
        _plt.xlabel('radii ratio (Req/Rp)')
        _plt.xlim([1.0, 1.5])
    _plt.savefig('n0_vs_ob.png')
    _plt.savefig('n0_vs_ob.pdf')
    return


def fs2rm_nan(fsed2, cols=[3], refcol=None):
    """ Remove ``nan`` values present in columns of a matrix file.
    In a *fullsed2* file, ``cols=[3]`` and ``refcol=[2]``.

    Input=filename, cols

    Output=file overwritten.

    TDB: refcol apparently is not working...
    """
    f2mtx = _np.loadtxt(fsed2, skiprows=5)
    for col in cols:
        nans, x = _phc.nan_helper(f2mtx[:, col])
        if refcol is None:
            f2mtx[:, col][nans] = _np.interp(
                x(nans), x(~nans), f2mtx[:, col][~nans])
        else:
            _warn.warn('# The output must be checked!')
            f2mtx[:, col][nans] = _np.interp(f2mtx[:, refcol][nans], 
                f2mtx[:, refcol][~nans], f2mtx[:, col][~nans])
    #
    oldf = open(fsed2).read().split('\n')
    if _np.max(f2mtx[_np.isfinite(f2mtx)]) < 100000:
        fmt = '%13.6f'
    elif _np.max(f2mtx[_np.isfinite(f2mtx)]) < 1000000:
        fmt = '%13.5f'
    else:
        raise ValueError('# ERROR at max values of fullsed2 {0}!!!!!!'.format(
            fsed2))
    _np.savetxt(fsed2, f2mtx, header='\n'.join(
        oldf[:5]), comments="", fmt=fmt, delimiter='')
    print('# {0} file updated!'.format(fsed2))
    return


def plottemp(tfiles, philist=[0], interpol=False, xax=0, fmt=['png'],
    figname=None, tlabels=None, title=None):
    """
    :Example:

    >>> import pyhdust as hdt
    >>> tfiles, tlabels = hdt.gentemplist('bestar2.02/mod01/mod01b33.temp', 
        tfrange=[30,33])
    >>> hdt.plottemp(tfiles, tlabels=tlabels)

    .. image:: _static/hdt_plottemp.png
        :width: 512px
        :align: center
        :alt: hdt.plottemp example

    `tfile` = filenames list.

    `xax` = 0: log10(r/R-1), 1: r/R; 2: 1-R/r

    `figname` = prefix of the output images

    `fmts` = format os the output images.

    If interpol==True, what will be plotted is the population along rays of
    a given latitude. The latitudes are defined in array muplot below.

    If interpola==False, what will be plotted is the population for a given mu
    index as a function of radius, starting with index ncmu/2(midplane) + plus

    OUTPUT = ...
    """
    if isinstance(tfiles, _strtypes):
        tfiles = [tfiles]
    if title is None:
        title = tfiles[-1]
    if interpol:
        raise NotImplementedError('# Interpol option not yet available.')
    if xax not in [0, 1, 2]:
        raise ValueError('# Invalid `xax` option. Try again.')
    #
    fig, axs = _plt.subplots(1, 1)  # , figsize=(21./3,29.7/3), sharex=True)
    lev = 0
    np = 0
    rplus = 0
    for i in range(len(tfiles)):
        rtfile = tfiles[i]
        ncr, ncmu, ncphi, nLTE, nNLTE, Rstar, Ra, beta, data, pcr, pcmu, pcphi\
            = readtemp(rtfile, quiet=True)
        for phiidx in range(0, len(philist)):
            icphi = philist[phiidx]
            x = data[0, :, 0, icphi]
            if (xax == 0):
                x = _np.log10(x / Rstar - 1.)
                xtitle = r'$\log_{10}(r/R_*-1)$'
            elif (xax == 1):
                x = x / Rstar
                xtitle = r'$r/R_*$'
            elif (xax == 2):
                x = 1. - Rstar / x
                xtitle = r'$1-R_*/r$'
            y = data[3 + lev, :, ncmu / 2 + np + rplus, icphi]
            y = y / 1000.
            mk = 'o:'
            if rtfile.find('avg') > 0:
                mk = 'o-'
            if tlabels is not None:
                axs.plot(x, y, mk, label=tlabels[i])
            else:
                axs.plot(x, y, mk)
    #
    axs.legend(loc='best', fancybox=True, framealpha=0.5)
    axs.set_title(title)
    axs.set_xlabel(xtitle)
    axs.set_ylabel(r'Temperature (10$^3$ K)')
    _phc.savefig(fig, figname=figname, fmt=fmt)
    return


def plotdens(tfiles, philist=[0], interpol=False, xax=0, fmt=['png'],
    outpref=None, tlabels=None, title=None):
    """
    VARIABLES:
        interpol: 
            - True: What will be plotted is the population along rays of a
                given latitude. The latitudes are defined in array muplot 
                below.
            - False: What will be plotted is the population for a given mu
                index as a function of radius, starting with index
                ncmu/2(midplane) + plus
        tfile:   file names' list
        xax:     0: log10(r/R-1), 1: r/R; 2: 1-R/r
        figname: prefix of the output images
        fmts:    format of the output images.

    :Example:
        >>> import pyhdust as hdt
        >>> tfiles, tlabels = hdt.gentemplist('bestar2.02/mod01/mod01b33.temp', 
                tfrange=[30,33])
        >>> hdt.plottemp(tfiles, tlabels=tlabels)
            .. image:: _static/hdt_plottemp.png
                :width: 512px
                :align: center
                :alt: hdt.plottemp example

    OUTPUT = ...
    """
    # despo SP 160908: copied from pyhdust's plottemp
    if isinstance(tfiles, _strtypes):
        tfiles = [tfiles]
    if title is None:
        title = tfiles[-1]
    if interpol:
        raise NotImplementedError('# Interpol option not available.')
    if xax not in [0, 1, 2, 3]:
        raise ValueError('# Invalid `xax` option. Try again.')
    #
    fig, axs = _plt.subplots(1, 1)  # , figsize=(21./3,29.7/3), sharex=True)
    np = 0
    rplus = 0
    for i in range(len(tfiles)):
        rtfile = tfiles[i]
        ncr, ncmu, ncphi, nLTE, nNLTE, Rstar, Ra, beta, data, pcr, pcmu, \
            pcphi = readtemp(rtfile, quiet=True)
        # print(rtfile, nLTE, nNLTE)
        for phiidx in range(0, len(philist)):
            icphi = philist[phiidx]
            x = data[0, :, 0, icphi]
            if (xax == 0):
                x = _np.log10(x / Rstar - 1.)
                xtitle = r'$\log_{10}(r/R_*-1)$'
            elif (xax == 1):
                x = x / Rstar
                xtitle = r'$r/R_*$'
            elif (xax == 2):
                x = 1. - Rstar / x
                xtitle = r'$1-R_*/r$'
            elif (xax == 3):
                x = _np.log10(x / Rstar)
                xtitle = r'$\log_{10}(r/R_*)$'
            # despo: number density numbers appear at column 5
            #        temperature is at lev+3
            y = data[5+nLTE, :, ncmu / 2 + np + rplus, icphi]
            y = _np.log10(y)  # /.6/mp)
            mk = 'o:'
            if rtfile.find('avg') > 0:
                mk = 'o-'
            if tlabels is not None:
                axs.plot(x, y, mk, label=tlabels[i])
            else:
                axs.plot(x, y, mk)
    #
    axs.legend(loc='best', fancybox=True, framealpha=0.5)
    axs.set_title(title)
    axs.set_xlabel(xtitle)
    axs.set_ylabel(r'$\log_{10}(n)$')
    _phc.savefig(fig, figname=outpref, fmt=fmt)
    return


def plottemp2d(tfile, figname=None, fmt=['png'], icphi=0, itype='linear',
    nimg=128, Rmax=None, trange=None, hres=False):
    """ itype = 'nearest', linear' or 'cubic'

    ``xax`` = 0: log10(r/R-1), 1: r/R; 2: 1-R/r
    """
    xax = 1
    if xax not in [0, 1, 2]:
        raise ValueError('# Invalid `xax` option. Try again.')
    lev = 0
    ncr, ncmu, ncphi, nLTE, nNLTE, Rstar, Ra, beta, data, pcr, pcmu, pcphi\
        = readtemp(tfile)
    #
    fig, ax = _plt.subplots()
    # pcmu[0, :] = 2*pcmu[1, :] - pcmu[2, :]
    # pcmu[-1, :] = 2*pcmu[-2, :] - pcmu[-3, :]
    ccmu = _np.arccos(pcmu[:-1] + _np.diff(pcmu, axis=0)/2.).flatten()
    ccr = _np.tile(pcr[:-1] + _np.diff(pcr)/2., ncmu)
    tmp = _phc.sph2cart(ccr, ccmu)
    x = tmp[1]
    y = tmp[0]
    #
    if (xax == 0):
        idx = _np.where(y > Rstar)
        x = _np.log10(x[idx]/Rstar - 1.)
        y = _np.log10(y[idx]/Rstar - 1.)
        idata = data[3+lev, :, :, icphi].T.flatten()[idx]/1e3
        xtitle = r'$\log_{10}(r/R_*-1)$'
    elif (xax == 1):
        idata = data[3+lev, :, :, icphi].T.flatten()/1e3
        idx = _np.where(idata < data[3+lev, 0, 0, icphi]/1e3)
        idata = idata[idx]
        x = x[idx]/Rstar
        y = y[idx]/Rstar
        xtitle = r'$r/R_*$'
    elif (xax == 2):
        idx = _np.where(y > Rstar)
        x = 1.-Rstar/x[idx]
        y = 1.-Rstar/y[idx]
        idata = data[3+lev, :, :, icphi].T.flatten()[idx]/1e3
        xtitle = r'$1-R_*/r$'
    if Rmax is None:
        xmax = _np.max(x)
    else:
        if xax == 0: 
            xmax = _np.log10(Rmax - 1.)
        elif xax == 1:
            xmax = Rmax
        elif xax == 2:
            xmax = 1-1./Rmax
    #
    ymax = (xmax-_np.min(x))/2.
    #
    coords = _np.column_stack((x, y))
    ax.set_xlim([_np.min(x), xmax])
    ax.set_ylim([-ymax, ymax])
    # xo, yo = _np.meshgrid( _np.linspace(_np.min(x), _np.max(x), nimg), 
    #     _np.linspace(_np.min(y), _np.max(y), nimg) )
    xo, yo = _np.meshgrid( _np.linspace(_np.min(x), xmax, nimg), 
        _np.linspace(-ymax, ymax, nimg) )
    msgri = _np.column_stack((xo.flatten(), yo.flatten()))
    # xo = _np.linspace(_np.min(ccr), _np.max(ccr), 21)
    # yo = _np.linspace(_np.min(ccmu), _np.max(ccmu), 21) 
    vmin = _np.min(idata)
    vmax = _np.max(idata)
    if trange is not None:
        vmin = trange[0]
        vmax = trange[-1]        
    if hres:
        img = _interpolate.griddata(coords, idata, msgri, 
            method=itype).reshape((nimg, nimg))
    else:
        img = _phc.baricent_map(x, y, idata)
    cax = ax.imshow(img, origin='lower', extent=[_np.min(x), xmax, -ymax, 
        ymax], vmin=vmin, vmax=vmax, cmap='gist_heat', 
        interpolation='bilinear')
    ax.set_title(_os.path.basename(tfile))
    ax.set_xlabel(xtitle)
    ax.set_ylabel(xtitle)
    cbar = fig.colorbar(cax, label=r'Temp. (10$^3$ K)')  
    # , orientation='horizontal')
    # cbar.ax.set_yticklabels(['< -1', '0', '> 1'])  # vertically colorbar
    # if Rmax is None:
    #     Rmax = Ra
    # else:
    #     Rmax *= Rstar
    # ymax = abs(Rmax/Rstar-1)/2
    # ax.set_xlim([1, Rmax/Rstar])
    # ax.set_ylim([-ymax, ymax])
    ax.set_aspect(abs(xmax-_np.min(x))/abs(ymax-(-ymax)))
    _phc.savefig(fig, figname=figname, fmt=fmt)
    return


def plotdens2d(tfile, figname=None, fmt=['png'], icphi=0, itype='linear',
    nimg=128, Rmax=None, trange=None, hres=False, zlim=None):
    """ itype = 'nearest', linear' or 'cubic'

    :param hres: high-resolution mode
    :param zlim: limits the z-axis (color) scale. Ex.: `[6, 13]` sets the log
    scale to :math:`10^6-10^{13}` .
    """
    xax = 1
    if xax not in [0, 1, 2]:
        raise ValueError('# Invalid `xax` option. Try again.')
    ncr, ncmu, ncphi, nLTE, nNLTE, Rstar, Ra, beta, data, pcr, pcmu, pcphi\
        = readtemp(tfile)
    lev = nLTE+2
    #
    fig, ax = _plt.subplots()
    # pcmu[0, :] = 2*pcmu[1, :] - pcmu[2, :]
    # pcmu[-1, :] = 2*pcmu[-2, :] - pcmu[-3, :]
    ccmu = _np.arccos(pcmu[:-1] + _np.diff(pcmu, axis=0)/2.).flatten()
    ccr = _np.tile(pcr[:-1] + _np.diff(pcr)/2., ncmu)
    tmp = _phc.sph2cart(ccr, ccmu)
    x = tmp[1]
    y = tmp[0]
    #
    if (xax == 0):
        idx = _np.where(y > Rstar)
        x = _np.log10(x[idx]/Rstar - 1.)
        y = _np.log10(y[idx]/Rstar - 1.)
        idata = data[3+lev, :, :, icphi].T.flatten()[idx]
        xtitle = r'$\log_{10}(r/R_*-1)$'
    elif (xax == 1):
        idata = ( data[3+lev, :, :, icphi].T.flatten() )
        idx = _np.where(idata > data[3+lev, 0, 0, icphi])
        idata = _np.log10(idata[idx])
        x = x[idx]/Rstar
        y = y[idx]/Rstar
        xtitle = r'$r/R_*$'
    elif (xax == 2):
        idx = _np.where(y > Rstar)
        x = 1.-Rstar/x[idx]
        y = 1.-Rstar/y[idx]
        idata = data[3+lev, :, :, icphi].T.flatten()[idx]
        xtitle = r'$1-R_*/r$'
    if Rmax is None:
        xmax = _np.max(x)
    else:
        if xax == 0: 
            xmax = _np.log10(Rmax - 1.)
        elif xax == 1:
            xmax = Rmax
        elif xax == 2:
            xmax = 1-1./Rmax
    #
    ymax = (xmax-_np.min(x))/2.
    #
    coords = _np.column_stack((x, y))
    ax.set_xlim([_np.min(x), xmax])
    ax.set_ylim([-ymax, ymax])
    # xo, yo = _np.meshgrid( _np.linspace(_np.min(x), _np.max(x), nimg), 
    #     _np.linspace(_np.min(y), _np.max(y), nimg) )
    xo, yo = _np.meshgrid( _np.linspace(_np.min(x), xmax, nimg), 
        _np.linspace(-ymax, ymax, nimg) )
    msgri = _np.column_stack((xo.flatten(), yo.flatten()))
    # xo = _np.linspace(_np.min(ccr), _np.max(ccr), 21)
    # yo = _np.linspace(_np.min(ccmu), _np.max(ccmu), 21) 
    vmin = _np.min(idata)
    vmax = _np.max(idata)
    if zlim is not None:
        vmin, vmax = zlim
    if trange is not None:
        vmin = trange[0]
        vmax = trange[-1]
    if hres:
        img = _interpolate.griddata(coords, idata, msgri, 
            method=itype).reshape((nimg, nimg))
    else:
        img = _phc.baricent_map(x, y, idata, fullrange=False)
    cax = ax.imshow(img, origin='lower', 
        extent=[_np.min(x), xmax, -ymax, ymax], vmin=vmin, 
        vmax=vmax, cmap='gist_heat', interpolation='bilinear')
    ax.set_title(_os.path.basename(tfile))
    ax.set_xlabel(xtitle)
    ax.set_ylabel(xtitle)
    cbar = fig.colorbar(cax, label=r'$\log_{10}(d)$ (cm$^{-3}$)')  
    # , orientation='horizontal')
    # cbar.ax.set_yticklabels(['< -1', '0', '> 1'])  # vertically colorbar
    # if Rmax is None:
    #     Rmax = Ra
    # else:
    #     Rmax *= Rstar
    # ymax = abs(Rmax/Rstar-1)/2
    # ax.set_xlim([1, Rmax/Rstar])
    # ax.set_ylim([-ymax, ymax])
    ax.set_aspect(abs(xmax-_np.min(x))/abs(ymax-(-ymax)))
    _phc.savefig(fig, figname=figname, fmt=fmt)
    return


# Be quantities convertions
def rho2sigp(R, rho0, a, M):
    sig = (2 * _np.pi) ** .5 * a / (_phc.G.cgs * M / R) ** .5 * R * rho0
    return sig


def rho2Mdot(R, alpha, a, M, rho0, R0):
    Mdot = 3 * _np.pi * (2 * _np.pi) ** .5 * alpha * a ** 3. / \
        (_phc.G.cgs * M / R) * rho0 * R ** 2. * ((R0 / R) ** .5 - 1) ** -1.
    return Mdot


def Mdot2sig(R, Mdot, alpha, a, M, R0):
    sig = Mdot * (_phc.G.cgs * M / R) ** .5 / \
        (3. * _np.pi * alpha * a ** 2 * R) * ((R0 / R) ** .5 - 1)
    return sig


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
    a = (_phc.kB.cgs * T / mu / _phc.mH.cgs) ** .5
    rho0 = rho0 * mu * _phc.mH.cgs
    sigp = rho2sigp(R, rho0, a, M)
    rho0p = sigp / (2 * _np.pi) ** .5 / a * (_phc.G.cgs * M / R) ** .5 / R
    Mdot = rho2Mdot(R, alpha, a, M, rho0, R0)
    sig = Mdot2sig(R, Mdot, alpha, a, M, R0)
    # sigl = Mdot * (_phc.G.cgs * M) ** .5 / 3 / _np.pi
    Mdisk0 = (2 * _np.pi) ** 1.5 * rho0 * R ** 2. * \
        (Rd - R) * a / (_phc.G.cgs * M / R) ** .5
    Mdisk = 2 * _np.pi * Mdot * \
        (_phc.G.cgs * M / R) ** .5 * R ** .5 / 3 / _np.pi / \
        alpha / a ** 2. * R0 ** .5 * _np.log( Rd / R)
    if Mdisk == 0:
        Mdisk = 4 * _np.pi * Mdot * \
            (_phc.G.cgs * M / R) ** .5 * R ** .5 / 3 / \
            _np.pi / alpha / a ** 2. * (Rd ** .5 - R ** .5)
    MdiskG = 2 * Mdot * (_phc.G.cgs * M / R) ** .5 * R ** .5 / 3 / alpha /\
        a ** 2 * (R0 ** .5 * _np.log(Rd / R) + 2 * R ** .5 - 2 * Rd ** .5)

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


def n0toSigma0(n0, M, Req, f, Tp, mu):
    """ VDD Steady-State conversion between `n0` (particles[ionized]/volume) 
    to `Sigma0`.

    INPUT: n0 (float), M, Req (mass and equatorial radius, Solar units),
    fraction and polar temperature ([0-1], Kelvin) and mu molecular weight
    [0.5-2.0].

    OUTPUT: float (g cm-2) """
    rho0 = n0 * mu * _phc.mH.cgs
    a = (_phc.kB.cgs * f * Tp / mu / _phc.mH.cgs) ** .5
    sig0 = (2 * _np.pi) ** .5 * a / (_phc.G.cgs * M * _phc.Msun.cgs / (Req *
        _phc.Rsun.cgs)) ** .5 * Req * _phc.Rsun.cgs * rho0
    return sig0


def n0toMdot(n0, M, Req, f, Tp, mu, alpha, R0):
    """ VDD Steady-State conversion between `n0` to `Mdot`.

    INPUT: n0 (float), M, Req (mass and equatorial radius, Solar units),
    fraction and polar temperature ([0-1], Kelvin), mu molecular weight
    [0.5-2.0], alpha (viscous parameter), R0 ("truncation" radius, Solar unit).

    OUTPUT: float (Msun yr-1)"""
    Req = Req * _phc.Rsun.cgs
    R0 = R0 * _phc.Rsun.cgs
    M = M * _phc.Msun.cgs
    rho0 = n0 * mu * _phc.mH.cgs
    a = (_phc.kB.cgs * f * Tp / mu / _phc.mH.cgs) ** .5
    Mdot = 3 * _np.pi * (2 * _np.pi) ** .5 * alpha * a ** 3. / (_phc.G.cgs *
        M / Req) * rho0 * Req ** 2. * ((R0 / Req) ** .5 - 1) ** -1.
    return Mdot / _phc.Msun.cgs * _phc.yr.cgs


def calcTeff(lum, size, M=None):
    """
    Calculate Teff for the non-rotating case.

    INPUT: Lum, Radius (or Mass in Solar units). `size` variable is assumed to
    be  the stellar radius (i.e., M==None). If M is given, size is assumed to
    be log(g) (cgs units).

    OUTPUT: float (Kelvin)
    """
    # M, Rp, Lum = fundline
    if M is None:
        Rp = size * _phc.Rsun.cgs
    else:
        Rp = (M * _phc.Msun.cgs * _phc.G.cgs / 10 ** size) ** .5
    L = lum * _phc.Lsun.cgs
    # Lum = 4*_np.pi*Rp**2*_phc.sigma*Teff**4
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


# Astro Useful and Plots
def loadfits(filename, hdu=0):
    """ Load a generic FITS. This function is different of the function as 
    pyhdust.spectools.

    if `hdu` is not a integer, it returns a list with all headers and data.

    OUTPUT: header, data 
    """
    fits = _pf.open(filename)
    if isinstance(hdu, int):
        header = fits[hdu].header
        data = fits[hdu].data
    else:
        header, data = ([], [])
        for fitshdu in fits:
            header.append(fitshdu.header)
            data.append(fitshdu.data)
    return header, data


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
        if ((long(10000) * year + long(100) * month + day) <= long(15821004)):
            b = -2 + int(modf((year + 4716) / 4)[1]) - 1179
        else:
            b = int(modf(year / 400)[1]) - int(modf(year / 100)[1]) + \
                int(modf(year / 4)[1])

        mjdmidnight = long(365) * year - long(679004) + \
            b + int(30.6001 * (month + 1)) + day

        fracofday = base60_to_decimal(
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
                print('ok0')
        elif hmin > hnas + cor and hmax > hpoe and (hpoe - hmin) < -12:
            hnas = hmin
            if debug:
                print('ok0')
        elif hmin > hnas + cor and hmax < hpoe and hnas < hmax:
            hnas = hmin
            if debug:
                print('ok1')
        elif hmin < hnas + cor and hmax < hpoe and hnas < hmax:
            hpoe = hmax
            if debug:
                print('ok2')
        elif hmin < hnas + cor and hmax < hpoe and (hnas - hmax) > 12:
            hpoe = hmax
            if debug:
                print('ok2a')
        elif hmin < hnas + cor and hmax > hpoe and (hpoe - hmin) < -12:
            if debug:
                print('ok3')
        elif hmin < hnas + cor and hmax > hpoe and hpoe > hmin:
            if debug:
                print('ok4')
        else:
            if debug:
                print(hmin, hnas, hmax, hpoe)
            hnas = hpoe = _np.NaN
        return (hnas, hpoe)

    # equinocio set 2011 (djsol)
    djsol = 2455827.87778

    # Input: dia gregoriano (dg), para dia juliano (dj) as 0h local
    # +3 do UT #OPD@LNA
    ut = -3.

    print("# Programa de planejamento de alvos BEACON\n#")
    print("# Horas em UT=%d ! (horario de 'inverno' de Brasilia)" % (ut))
    print("# Digite a Data de Observacao, ou 'ENTER' para hoje...\n#")
    dg = []
    dg = dg + [str(_phc.user_input("Digite o Ano (xxxx): "))]
    try:
        j = float(dg[-1])
        dg = dg + [str(_phc.user_input("Digite o Mes (1-12): "))]
        dg = dg + [str(_phc.user_input("Digite o Dia (1-31): "))]
    except:
        now = _dt.datetime.now()
        dg = []
        dg = dg + [now.year]
        dg = dg + [now.month]
        dg = dg + [now.day]

    dg = _np.array(dg, dtype=float)
    dj = julian_date(dg[0], dg[1], dg[2] + 1, 0 + ut, 0, 0.)
    # dj = julian_date(now.year, now.month, now.day, now.hour, now.minute, 
        # now.second)

    # dias/ano passados do solsticio de verao (ndays)
    ndays = dj - djsol
    while ndays > 365:
        ndays = ndays - 365
    while ndays < 0:
        ndays = ndays + 365

    # Hora sideral 'as 0h local (hz) em horas
    hz = ndays * (0.065711111)

    # Reducao do tempo de observacao (rt) em minutos
    if ndays <= 91.25:
        rt = 75 + 15 / 91.25 * ndays
    elif ndays <= 273.75:
        rt = 90 - (ndays - 91.25) * 180 / 365
    else:
        rt = (ndays - 273.75) * 75 / 91.25
    # rt de min. para horas
    rt = rt / 60

    # hora maximo (hmax) e minimo (hmin) de observacao
    hmin = hz - 6 + rt / 2
    if hmin < 0:
        hmin = hmin + 24
    hmax = hz + 6 - rt / 2
    if hmax > 24:
        hmax = hmax - 24

    # observacao minima e maxima em dj (dmin, dmax)
    # rt de horas para segs.
    # rt = rt*3600.
    # dmin = julian_date(dg[0],dg[1],dg[2]+1,3-6,0,rt/2)
    # dmax = julian_date(dg[0],dg[1],dg[2]+1,3+6,0,-rt/2) #seg. so' >0!!!

    # carrega lista de alvos
    alvos = _np.loadtxt('{0}refs/obs_alvos.txt'.format(hdtpath()), 
        dtype=str, delimiter='\t')
    # carrega tempo das declinacoes
    obsdec = _np.loadtxt(
        '{0}refs/obs_dec.txt'.format(hdtpath()), delimiter='\t')
    # carrega efemerides
    if _os.path.exists('{0}refs/obs_ef.txt'.format(hdtpath())):
        ef_alvos = _np.loadtxt('{0}refs/obs_ef.txt'.format(hdtpath()),
                               delimiter='\t', dtype=str)
        ef_alvos = ef_alvos.T

    print("\n# Info. da noite: %2d %2d %4d" % ( dg[2], dg[1], dg[0]) )
    print(  "Tempo sideral no anoitecer : %2d %2d" %
          ( int(hmin), int((hmin - int(hmin)) * 60) )  )
    print(  "Tempo sideral 'a meia noite: %2d %2d" %
          ( int(hz), int((hz - int(hz)) * 60) )  )
    print(  "Tempo sideral no amanhecer : %2d %2d" %
          ( int(hmax), int((hmax - int(hmax)) * 60) )  )
    print(  "Dif. de observ. ao inverno : -%d %2d" %
          ( int(rt), int((rt - int(rt)) * 60) )  )
    print("\nALVO\tINICIO\tFIM\tMERID.\tF_INI\tF_FIM")

    outfile = "ALVO\tINICIO\tFIM\tMERID.\tF_INI\tF_FIM\n"
    for i in range(len(alvos)):
        # calcula a ascencao reta (ra) e declinacao (dec)
        ra = float(alvos[i][2][:2]) + float(alvos[i][2][3:5]) / 60
        dec = float(alvos[i][3][:3]) + \
            float(alvos[i][3][0] + alvos[i][3][4:6]) / 60
        # calcula tempo observavel na declinacao (to)
        j = odmin = odmax = 0
        while odmax == 0:
            decmin = obsdec[j][0]
            odmin = obsdec[j][1]
            if obsdec[j + 1][0] < dec:
                decmax = obsdec[j + 1][0]
                odmax = obsdec[j + 1][1]
            j = j + 1

        to = (abs(dec - decmax) * odmin + abs(dec - decmin) * odmax) / \
            (abs(dec - decmin) + abs(dec - decmax))

        # calcula para cada alvo o tempo de observacao
        # hora de inicio da observacao ('nascer')
        hnas = ra - to / 2
        if hnas < 0:
            hnas = hnas + 24
        # hora de termino da observacao ('poente')
        hpoe = ra + to / 2
        if hpoe > 24:
            hpoe = hpoe - 24
        # observavel ao anoitecer?
        (hnas, hpoe) = obs_ver(hmin, hnas, hmax, hpoe)
        # estara' observavel por mais de um 3/4 de hora? (dt)
        dt = hpoe - hnas
        if dt < 0:
            dt = dt + 24
        if dt > 12:
            dt = dt - 24
        if dt < 3 / 4.:
            hnas = _np.NaN
            hpoe = _np.NaN

        # procura posicao nas efemerides (pef)
        if _os.path.exists('{0}refs/obs_ef.txt'.format(hdtpath())):
            pef = [k for k, x in enumerate(
                ef_alvos[1]) if x.find(alvos[i][1]) > -1]
        else:
            pef = []
        if len(pef) == 1 and hnas > 0:
            pef = pef[0]
            per = float(ef_alvos[2][pef])
            T0 = float(ef_alvos[3][pef])
            # calcula fases
            # dj dia juliano 'a meia-noite local
            # dj = julian_date(dg[0],dg[1],dg[2]+1,0+ut,0,0.)
            # observacao ANTES da 0h UT
            cor = j = 0
            if (hnas - hz) > 12:
                cor = -24
            if hnas - hz + ut + cor < 0:  # or (hz-hnas < 0 and hz-nas<-12):
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
            # print(j,hour,djnas)
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
            # print(j,hour,djpoe)
            # Warnings
            # if T0 == 0.:
                # print("WARNING: No 'T0' is available for this star!!!")
            # if per == 1.:
                # print("WARNING: No Period is available for this star!!!")
            # print("# A fase atual e' "+str( (j-T0)%Per/Per )+"\n" )
            fase_inc = ((djnas - T0) % per) / per
            fase_fim = ((djpoe - T0) % per) / per
        else:
            fase_inc = _np.NaN
            fase_fim = _np.NaN

        # tempos locais
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
                alvos[i][0][:7], int(hlmin), int((hlmin - int(hlmin)) * 60), 
                int(hlmax), int((hlmax - int(hlmax)) * 60), int(hmer), 
                int((hmer - int(hmer)) * 60), fase_inc, fase_fim) )
            outfile = outfile + ("%s\t%2d %2d\t%2d %2d\t%2d %2d\t%.2f\t%.2f" %
                (alvos[i][0][:7], int(hlmin), int((hlmin - int(hlmin)) * 60), 
                int(hlmax), int((hlmax - int(hlmax)) * 60), int(hmer), 
                int((hmer - int(hmer)) * 60), fase_inc, fase_fim) ) + '\n'
        else:
            print("%s\t-- --\t-- --\t-- --\t-.--\t-.--" % (alvos[i][0][:7]) )
            outfile = outfile + ("%s\t-- --\t-- --\t-- --\t-.--\t-.--" %
                 (alvos[i][0][:7]) ) + '\n'

    wfile = _phc.user_input("\nDeseja salvar a lista?(Sim/outro): ")
    if wfile in ['s', 'sim', 'Sim', 'S', 'y', 'yes', 'Yes', 'Y']:
        f0 = open('obs_%2d_%2d_%4d.txt' % (dg[2], dg[1], dg[0]), 'w')
        f0.writelines(outfile)
        f0.close()
    return


def plot_obs(observ_dates=[], legend=[], civcfg=[1, 'm'], civdt=None,  
    fmt=['png'], addsuf=None, colors=None, ls=None, markers=None, 
    mjdlims=None, alpha=1.0, graymjds=[], dolines=False, mincivcfg=[10, 'd']):
    """ Plot observations as done in Faes+2016. 

    INPUT: `observ_dates` in a list of arrays. Each array contains the dates of 
    a given observation. Dates, in principle (TODO), must be in MJD.

    `legend` is a list of names, one for each array.

    `civcfg` = [step, 'd'|'m'|'y'] 

    `civdt` = starting date at the Civil Date axis (MJD?).

    if `addsuf` is None, a date-time string will be added to the filename.

    `colors` = colors for each observ_dates array. If `None`, then 
    phc.cycles(ctype='cor').

    `ls` same as above, for linestyle.

    `markers` same as above, for markers.

    `graymjds`, gray areas; example `graymjds =[[56200, 56260], [57000, 57100]]

    OUTPUT: File saved.
    """
    if civdt is not None:
        civdt = _dt.datetime(civdt[0], civdt[1], civdt[2]).date()
    k = len(observ_dates)
    if colors is None:
        colors = [_phc.cycles(c, 'cor') for c in range(k)]
    if ls is None:
        ls = [_phc.cycles(c, 'ls') for c in range(k)]
    if markers is None:
        mk = [_phc.cycles(c, 'mk') for c in range(k)]
    #
    fig, ax = _plt.subplots(figsize=(8, 3))
    ys = _np.linspace(0, 1, len(observ_dates) + 2)
    for i in range(k):
        for j in range(len(observ_dates[i])):
            dt = observ_dates[i][j]
            if dolines:
                ax.plot([dt, dt], [0, 1], color=colors[i], ls=ls[i], alpha=0.1)
            if j == 0:
                ax.scatter([dt], [ys[k-i]], color=colors[i], 
                    marker=mk[i], s=12, label=legend[i], alpha=alpha)
            else:
                ax.scatter([dt], [ys[k-i]], color=colors[i], 
                    marker=mk[i], s=12, alpha=alpha)
    for g in graymjds:
        rect = _mpatches.Rectangle([g[0], 0], g[1]-g[0], 1., ec="gray", 
            fc='gray', alpha=0.5, zorder=1)
        ax.add_patch(rect)
    # 
    flatdataes = [item for sublist in observ_dates for item in sublist]
    mjd0, mjd1 = (_np.min(flatdataes), _np.max(flatdataes))
    civvals = {'Y': 365, 'M': 30, 'D': 1}
    extratick = civcfg[0]*civvals[civcfg[1][0].upper()]
    ax.set_ylim([0, 1])
    # ax.set_title('Title')
    ax.legend(fontsize=10, loc='lower left', fancybox=True, scatterpoints=1)
    ax.set_xlabel('Julian date - 2400000.5')
    dtticks = _phc.gentkdates(mjd0, mjd1+extratick, civcfg[0], civcfg[1], 
        dtstart=civdt)
    mjdticks = [_jdcal.gcal2jd(date.year, date.month, date.day)[1] for date in 
        dtticks]
    for pair in zip(dtticks, mjdticks):
        print(pair)
    # ax.plot([mjdticks[-1], mjdticks[-1]], [0, 1], alpha=0)
    # xlim = [mjdticks[0], ax.get_xlim()[-1]]
    if mjdlims is None:
        xlim = [mjd0, mjd1+extratick]
    else:
        xlim = mjdlims
    # print xlim
    ax.minorticks_on()
    # ax.set_xlim(limits)
    ax2 = ax.twiny()
    ax2.spines['top'].set_position(('axes', 1.1))
    ax2.minorticks_on()
    ax2.set_xlabel('Civil date')
    dtminticks = _phc.gentkdates(xlim[0], xlim[1], mincivcfg[0], mincivcfg[1])
    i = 1
    while dtticks[0] not in dtminticks:
        dtminticks = _phc.gentkdates(xlim[0]+i, xlim[1], mincivcfg[0], 
            mincivcfg[1])
        i += 1
    minjdticks = [_jdcal.gcal2jd(date.year, date.month, date.day)[1] for date 
        in dtminticks]
    ax2.set_xticks(mjdticks)
    xlabs = [date.strftime('%Y-%m-%d') for date in dtticks]
    xlabs[1::2] = ['']*len(xlabs[1::2])
    ax2.set_xticklabels(xlabs)  
        # , fontsize=10)
    ax2.set_xticks(minjdticks, minor=True)
    # ax2.set_xticklabels([date.strftime("%Y/%M/%d") for date in dtticks])
    ax.set_yticks([])
    ax.set_xlim(xlim)
    ax2.set_xlim(xlim)
    # ax.xaxis.set_minor_locator(_ML(50))
    ax.xaxis.set_ticks_position('both')
    ax.xaxis.set_tick_params(length=8, width=1.5)
    ax.xaxis.set_tick_params(length=6, which='minor')
    ax2.xaxis.set_tick_params(length=4, which='minor')
    ax2.xaxis.set_tick_params(length=8, width=1.5)
    _plt.setp( ax2.xaxis.get_majorticklabels(), rotation=0 )
    # ax2.fmt_xdata = mdates.DateFormatter('%Y-%m-%d')
    # fig.autofmt_xdate()
    _plt.subplots_adjust(left=0.02, right=0.98, top=0.67, bottom=0.16, 
        hspace=.15)
    if addsuf is None:
        addsuf = _phc.dtflag()
    # for f in fmt:
    #     print('# Saved plot_obs{0}.{1}'.format(addsuf, f))
    #     _plt.savefig('plot_obs{0}.{1}'.format(addsuf, f), transparent=True)
    _phc.savefig(fig, fmt=fmt)
    _plt.close()
    return


def plotMJDdates(spec=None, pol=None, interf=None, limits=None):
    """
    Plot dates from spec (Class), pol (routines) and interf (ESO query)

    TODO: This need to be polished !!!!

    spec = 'data_aeri_splot.txt'
    pol = 'pol_aeri.log'
    interf = 'interf_aeri.txt'
    """
    fig, ax = _plt.subplots()

    if spec is not None:
        spJD = _np.loadtxt(spec)
        spJD = spJD[:, 0]
        y = [0. for JD in spJD]
        ax.plot(spJD, y, marker='d', color='lightgray', ls='')
        # yerr = [ [1. for JD in spJD], [1. for JD in spJD] ]
        # ax.errorbar(spJD, y, yerr, marker='o', color='blue', ls='')

    if pol is not None:
        polJD = _np.loadtxt(pol, dtype=str)
        polJD = polJD[:, 9]
        polJD = _np.array(polJD, dtype=float) - 2400000.5
        y = [-.5 for JD in polJD]
        ax.plot(polJD, y, marker='o', color='gray', ls='')
        # yerr = [ [.5 for JD in polJD], [1.5 for JD in polJD] ]
        # ax.errorbar(polJD, y, yerr, marker='x', color='green', ls='')

    if interf is not None:
        intJD = _np.loadtxt(interf, dtype=str, delimiter=',', skiprows=1)
        intJD = _np.array(intJD[:, -2], dtype=float)
        y = [.5 for JD in intJD]
        ax.plot(intJD, y, marker='s', color='darkgrey', ls='')
        # yerr = [ [1.5 for JD in intJD], [.5 for JD in intJD] ]
        # ax.errorbar(intJD, y, yerr, marker='s', color='red', ls='')

    limits = (56100., 56750.)
    if limits is None:
        mjd0, mjd1 = ax.get_xlim()
    else:
        mjd0, mjd1 = limits
        ax.set_xlim(limits)
    ticks = _phc.gentkdates(mjd0, mjd1, 3, 'm', dtstart=_dt.datetime(2012, 7,
        1).date())
    mjdticks = [_jdcal.gcal2jd(date.year, date.month, date.day)[1] for date in
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


def plot_gal(ra, dec, addsuf=None, fmt=['png']):
    """ Plot in "Galatic Coordinates" (i.e., Mollweide projection).

    Input: RA and DEC arrays in RADIANS (fraction)

    Output: Saved images """

    fig = _plt.figure()  # figsize=(8,6))
    ax = fig.add_subplot(111, projection="mollweide")
    ax.set_xticklabels(['14h', '16h', '18h', '20h', '22h', '0h', '2h', '4h', 
        '6h', '8h', '10h'])
    ax.scatter(ra, dec)

    if addsuf is None:
        addsuf = _phc.dtflag()
    for f in fmt:
        print('# Saved plot_gal{0}.{1}'.format(addsuf, f))
        _plt.savefig('plot_gal{0}.{1}'.format(addsuf, f), transparent=True)
    _plt.close()
    return


# Filters tools
def doFilterConv(x0, y0, filt, zeropt=False):
    r""" Return the convolved filter total flux for a given flux profile y0,
    at wavelengths x0.

    .. math::

        \text{zero}_{\rm pt} = \frac{\int F(\lambda) Sp(\lambda)\,d\lambda}
        {\int F(\lambda)\,d\lambda} 

    INPUT: x0 lambda array (**Angstroms**), y0 flux array, filter (string, 
    uppercase for standard filters)

    OUTPUT: summed flux (y0 units; default) or *zero point level* 
    (**polarimetry**)
    """
    fpath = _os.path.join(hdtpath(), 'refs', 'filters', filt+'.dat')
    fdat = _np.loadtxt(fpath, skiprows=1)
    #
    # idx = _np.where((x0 >= fdat[0, 0]) & (x0 <= fdat[-1, 0]))
    # x0 = x0[idx]
    # y0 = y0[idx]
    # interpfunc = interpolate.interp1d(fdat[:,0], fdat[:,1], 
    #     kind='linear')
    # interpfunc = _interpolate.InterpolatedUnivariateSpline(
    #     fdat[:, 0], fdat[:, 1])
    # dfintp = interpfunc(x0)
    fdintp = _np.interp(x0, fdat[:, 0], fdat[:, 1], left=0, right=0)
    if not zeropt:
        return _np.trapz(fdintp * y0, x0)
    else:
        return _np.trapz(fdintp * y0, x0)/_np.trapz(fdintp, x0)


def doPlotFilter(obs, filter, fsed2data, pol=False, addsuf=None, fmt=['png']):
    """
    obs = integer; filter = single string

    TODO = lambda standard in .../filters/*.dat
    """
    x0 = fsed2data[obs, :, 2]
    if pol:
        y0 = fsed2data[obs, :, 7]
    else:
        y0 = fsed2data[obs, :, 3] / x0

    fdat = _np.loadtxt('{0}filters/{1}.dat'.format(hdtpath(), filter.lower()),
                       skiprows=1)
    fdat[:, 0] /= 10000.
    # interpfunc = interpolate.interp1d(fdat[:,0], fdat[:,1], kind='linear')
    interpfunc = _interpolate.InterpolatedUnivariateSpline(
        fdat[:, 0], fdat[:, 1])

    idx = _np.where((x0 >= fdat[0, 0]) & (x0 <= fdat[-1, 0]))
    x0 = x0[idx]
    y0 = y0[idx]
    # y = interpfunc(x0) * y0  # /_np.sum( interpfunc(x0) )

    fig, ax = _plt.subplots()
    ax.plot(x0, y0, label='SED')
    # ax.plot(fdat[:,0], fdat[:,1], label='Filter')
    ax.plot(x0, interpfunc(x0) * y0, label='Convolved')
    ax.plot(fdat[:, 0], fdat[:, 1] * _np.max(y0), label='Filter (eff.)')
    if addsuf is not None:
        ax.set_title(addsuf)
    ax.legend(loc='best', fancybox=True, framealpha=0.5)

    if addsuf is None:
        addsuf = _phc.dtflag()
        if pol:
            addsuf += '_pol'
    for f in fmt:
        print('# Saved doPlotFilter{0}.{1}'.format(addsuf, f))
        _plt.savefig('doPlotFilter{0}.{1}'.format(addsuf, f), transparent=True)
    _plt.close()

    return


def plot_hdt_filters(outname=None):
    "Plot all filters available in PyHdust"
    filters = _glob( _os.path.join(hdtpath(), 'refs', 'filters', 
        '*.dat') )
    fig, axs = _plt.subplots(3, 1)
    for f in filters:
        data = _np.loadtxt(f)
        if data[0, 0] < 9000:
            i = 0
        elif data[0, 0] < 50000:
            i = 1
        else:
            i = 2
        axs[i].plot(data[:, 0], data[:, 1], label=_os.path.split(f)[1].
            replace('.dat', ''))
    for i in range(len(axs)):
        axs[i].legend(fontsize=6)
    axs[i].set_xlabel(r'Wavelength ($\AA$)')
    _phc.savefig(fig, figname=outname)
    # print('# Filtres ploted in file {0}'.format(outname))
    return 


def readphotxdr(xdrfile='kur_ap00k0.xdr', quiet=False):
    """ doc """
    f = open(xdrfile, 'rb').read()
    ixdr = 0
    ixdr, arr = _phc.readpck(4, 'l', ixdr, f)
    _, nlin, nlbd, nmod = arr

    ixdr, _ = _phc.readpck(1, 'l', ixdr, f)
    ixdr, lbdcentr = _phc.readpck(nlin, 'f', ixdr, f)

    ixdr, _ = _phc.readpck(1, 'l', ixdr, f)
    ixdr, ltemp = _phc.readpck(nmod, 'f', ixdr, f)

    ixdr, _ = _phc.readpck(1, 'l', ixdr, f)
    ixdr, llogg = _phc.readpck(nmod, 'f', ixdr, f)

    lambdas = _np.zeros((nlin, nlbd))
    for i in range(nlin):
        _ = _struct.unpack('>l', f[ixdr:ixdr+ 4])
        ixdr += 4
        istr = '>{0:d}f'.format(nlbd)
        lambdas[i] = _struct.unpack(istr, f[ixdr:ixdr+ 4*nlbd])
        ixdr += 4*nlbd

    profiles = _np.zeros((nmod, nlin, nlbd))
    for i, j in _itprod(range(nmod), range(nlin)):
        _ = _struct.unpack('>l', f[ixdr:ixdr+ 4])
        ixdr += 4
        istr = '>{0:d}f'.format(nlbd)
        profiles[i, j] = _struct.unpack(istr, f[ixdr:ixdr+ 4*nlbd])
        ixdr += 4*nlbd

    if ixdr == len(f):
        if not quiet:
            print('# XDR {0} completely read!'.format(xdrfile))
    else:
        _warn.warn('# XDR {0} not completely read!\n# length difference is '
            '{1}'.format(xdrfile, (len(f) - ixdr) / 4 ) )

    return lbdcentr, ltemp, llogg, lambdas, profiles

# def chkObsLog(path=None, nights=None, badweath=None):
    # """ Check if there is data for all nights with observations.
# 
    # If not, check if the night is in the list of night lost due to bad 
        # weather.
# 
    # If no data and no bad weather info is registered, prints an error.
# 
    # If the night is included as bad weather and is not in night list, prints
    # a warning.
    # """
    # if path == None:
        # path = _os.getcwd()
    # if nights == None:
        # nights = '{0}pyhdust/refs/noites.txt'.format(hdtpath())
    # lnights = _np.loadtxt(nights, dtype=str)
    # if badweath == None:
        # badweath = '{0}pyhdust/refs/maltempo.txt'.format(hdtpath())
    # lbadweath = _np.loadtxt(badweath, dtype=str)
    # for night in lnights:
        # if night in lbadweath:
            # pass
        # elif not _os.path.exists(night):
            # print('# ERROR! {0} has no data and was not lost for bad' +\
                # 'weather!'.format(night))
    # flds = [fld for fld in _os.listdir('{0}'.format(path)) if \
            # _os.path.isdir(_os.path.join('{0}'.format(path), fld))]
    # for fld in flds:
        # if fld not in lnights:
            # print('# Warning! Night {0} is not recorded as OPD night!'. \
                # format(fld))
            # print('# Update the file {0}'.format(nights))
    # for night in lbadweath:
        # if night not in lnights:
            # print('# Warning! Bad weather {0} is not recorded as OPD ' +\
                # 'night!'.format(night))
            # print('# Probably it is a spec night.')
    # return


def tefflum_dJN(s, b):
    """ Calculate the Teff and Lum from *s* and *b* variables from de Jager & 
    Niewuwenhuijzen (1987).

    return log10_T, log10_L (in Solar and Kelvin units)

    ====== ======== ==============
    Spec   s-value  s-step/0.1-Sp
    ====== ======== ==============
    01-09  0.1-0.9  0.1
    09-B2  0.9-1.8  0.3
    B2-A0  1.8-3.0  0.15
    A0-F0  3.0-4.0  0.1
    F0-G0  4.0-5.0  0.1
    G0-K0  5.0-5.5  0.05
    K0-M0  5.5-6.5  0.1
    M0-M10 6.5-8.5  0.2
    ====== ======== ==============

    ========== ========
    L class    b-value
    ========== ========
    V          5.0
    IV         4.0
    III        3.0
    II         2.0
    Ib         1.4
    Iab (or I) 1.0
    Ia         0.6
    Ia+        0.0
    ========== ========

    The calcs are made based on Chebyshev polynomials. 
    """
    cij = [
        [+3.82573, -2.13868, -0.46357, +0.02076, -0.11937],
        [-1.55607, -1.89216, -0.96916, -0.08869, -0.20423],
        [+1.05165, +0.42330, -0.94379, -0.07438],
        [-0.01663, -0.20024, -0.18552],
        [-0.07576, -0.10934],
        [+0.11008]]

    dij = [
        [+3.96105, +0.03165, -0.02963, +0.01307, -0.01172],
        [-0.62945, +0.02596, -0.06009, +0.01881, -0.01121],
        [+0.14370, -0.00977, -0.03265, +0.01649],
        [+0.00791, +0.00076, -0.03006],
        [+0.00723, -0.02621],
        [+0.02755]]

    def Tk(k, x):
        if not isinstance( k, ( int, long ) ):
            _warn.warn('# Wrong Tk call! Invalid k')
            return
        if k == 0:
            return 1
        if k == 1:
            return x
        return 2*x*Tk(k-1, x) - Tk(k-2, x)

    logL = 0
    for n in range(5):
        for i in range(n+1):
            j = n-i
            logL += cij[i][j]*Tk(i, (s-4.25)/4.25)*Tk(j, (b-2.5)/2.5)
    logT = 0
    for n in range(5):
        for i in range(n+1):
            j = n-i
            logT += dij[i][j]*Tk(i, (s-4.25)/4.25)*Tk(j, (b-2.5)/2.5)
    return logT, logL


def masslum_OB(logL, classL='V'):
    """ Mass determination based on luminosity and class (Claret 2004). 

    Warning! Valid only for late O to late B range (ie., Be stars).

    return log10_M (solar units)
    """
    cs = [0.1997, 0.0844, 0.0312]
    if classL.upper() == 'IV':
        cs = [0.2055, 0.0746, 0.0322]
    elif classL.upper() == 'III':
        cs = [0.2084, 0.0643, 0.0325]
    elif classL.upper() != 'V':
        _warn.warn('# Invalid "classL". Assuming "V" for masslum_OB...')
    return cs[0] + cs[1]*logL + cs[2]*logL**2

# MAIN ###
if __name__ == "__main__":
    pass
