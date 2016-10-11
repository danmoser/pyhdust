# -*- coding:utf-8 -*-

"""PyHdust *beatlas* module: BeAtlas specific variables and functions.

Module contains:
- BAstar class
- BAmod class

:license: GNU GPL v3.0 https://github.com/danmoser/pyhdust/blob/master/LICENSE
"""
from __future__ import print_function
import re as _re
import os as _os
import numpy as _np
import struct as _struct
from glob import glob as _glob
from itertools import product as _product
import pyhdust.phc as _phc
import pyhdust as _hdt
import warnings as _warn

try: 
    from scipy.interpolate import griddata as _griddata
except ImportError:
    _warn.warn('# scipy module not installed!')

__author__ = "Daniel Moser"
__email__ = "dmfaes@gmail.com"


class BAstar(object):

    """ BeAtlas source star filename structure.

    The filename must follow this structure (the sequence in filename is not 
    importante):
    - keys: ['M', 'ob', 'H', 'Z', 'b']
    - last-key (no identifier): 'Ell'

    See BAmod, phc.keys_values.
    """

    def __init__(self, f0):
        vals = _phc.keys_values(['M', 'ob', 'H', 'Z', 'b'], f0)
        self.M, self.ob, self.H, self.Z, self.beta = vals
        self.shape = f0.split('_')[-1].replace('.txt', '')
        self._f0 = f0

    def __repr__(self):
        return self._f0


class BAmod(BAstar):

    """ BeAtlas disk model filename structure.

    It could be f0.split('_'), but the f0.find('_X') way was chosen.

    See that the parameters sequence is not important for this reading (this
    may not be the case of other routines). And, by definition, the source star
    has a specific name added at the end of disk model name, starting with
    'Be_'. """

    def __init__(self, f0):
        """ Class initialiser """
        BAstar.__init__(self, f0[f0.find('Be_'):])
        self.param = False
        if f0.find('_PL') > -1:
            self.param = True
            self.n = _phc.keys_values(['PLn'], f0)
        self.sig, self.h, self.Rd = _phc.keys_values(['sig', 'h', 'Rd'], f0)

    def build(self, ctrlarr, listpars):
        """ Set full list of parameters. """
        for i in range(len(ctrlarr)):
            if i == 0:
                self.M = _phc.find_nearest(listpars[i], ctrlarr[i])
            if i == 1:
                self.ob = _phc.find_nearest(listpars[i], ctrlarr[i])
            if i == 2:
                self.Z = _phc.find_nearest(listpars[i], ctrlarr[i])
            if i == 3:
                self.H = _phc.find_nearest(listpars[i], ctrlarr[i])
            if i == 4:
                self.sig = _phc.find_nearest(listpars[i], ctrlarr[i])
            if i == 5:
                self.Rd = _phc.find_nearest(listpars[i], ctrlarr[i])
            if i == 6:
                self.h = _phc.find_nearest(listpars[i], ctrlarr[i])
            if len(listpars) == 9:
                if i == 7:
                    # print(_phc.find_nearest(listpars[i], ctrlarr[i]))
                    self.n = _phc.find_nearest(listpars[i], ctrlarr[i])
                    self.param = True
                if i == 8:
                    self.cosi = _phc.find_nearest(listpars[i], ctrlarr[i])
            else:
                if i == 7:
                    self.cosi = _phc.find_nearest(listpars[i], ctrlarr[i])

    def getidx(self, minfo):
        """ Find index of current model in minfo array. """
        if len(minfo[0]) == 9:
            self.idx = (minfo[:, 0] == self.M) & (minfo[:, 1] == self.ob) &\
                (minfo[:, 2] == self.Z) & (minfo[:, 3] == self.H) &\
                (minfo[:, 4] == self.sig) & (minfo[:, 5] == self.Rd) &\
                (minfo[:, 6] == self.h) & (minfo[:, 7] == self.n) &\
                (minfo[:, -1] == self.cosi)
        else:
            self.idx = (minfo[:, 0] == self.M) & (minfo[:, 1] == self.ob) &\
                (minfo[:, 2] == self.Z) & (minfo[:, 3] == self.H) &\
                (minfo[:, 4] == self.sig) & (minfo[:, 5] == self.Rd) &\
                (minfo[:, 6] == self.h) & (minfo[:, -1] == self.cosi)
        return self.idx

# Only for H=0.30
vrots = [
    [259.759, 354.834, 417.792, 464.549, 483.847],
    [252.050, 346.163, 406.388, 449.818, 468.126],
    [245.127, 336.834, 399.983, 448.076, 467.806],
    [239.522, 329.496, 388.734, 432.532, 450.806],
    [234.301, 321.139, 379.297, 423.241, 441.122],
    [228.538, 313.797, 370.343, 412.488, 429.914],
    [219.126, 299.656, 354.547, 395.821, 413.008],
    [211.544, 288.840, 341.081, 380.426, 396.978],
    [203.438, 279.328, 328.666, 365.697, 380.660],
    [197.823, 268.964, 316.901, 353.568, 368.506],
    [192.620, 262.688, 308.208, 341.963, 356.410],
    [187.003, 255.125, 299.737, 332.511, 346.043]]

obs = [1.1, 1.2, 1.3, 1.4, 1.45]
ms = [14.6, 12.5, 10.8, 9.6, 8.6, 7.7, 6.4, 5.5, 4.8, 4.2, 3.8, 3.4]
xc = [0.08, 0.30, 0.42, 0.54, 0.64, 0.77]
ns = [3.0, 3.5, 4.0, 4.5]
nia = [0.5]
z = [0.014]
h = [72]
Rd = [50.0]
sig0 = _np.logspace(_np.log10(0.02), _np.log10(4.0), 7)

# Only for H=0.30
Tp11 = _np.array([28905.8, 26945.8, 25085.2, 23629.3, 22296.1, 20919.7,
18739.3, 17063.8, 15587.7, 14300.3, 13329.9, 12307.1])
Ms = _np.array([14.6, 12.5, 10.8, 9.6, 8.6, 7.7, 6.4, 5.5, 4.8, 4.2, 3.8, 3.4],
dtype=str)
Sig0 = ['{0:.2f}'.format(x) for x in sig0]


def rmMods(modn, Ms, clusters=['job']):
    """
    Remove the *.inp models of models `modn` according to the list structure
    below.

    | Masses list ans sig0 POSITION do be excluded
    | Ms = [
    | ['14.6', [0]],
    | ['12.5', [0,-1]],
    | ['10.8', [0,-1]],
    | ['09.6',  [0,-2,-1]],
    | ['08.6',  [0,-2,-1]],
    | ['07.7',  [0,-2,-1]],
    | ['06.4',  [0,-3,-2,-1]],
    | ['05.5',  [0,-3,-2,-1]],
    | ['04.8',  [-4,-3,-2,-1]],
    | ['04.2',  [-4,-3,-2,-1]],
    | ['03.8',  [-4,-3,-2,-1]],
    | ['03.4',  [-4,-3,-2,-1]],]

    INPUT: string, structured list

    OUTPUT: *files removed
    """
    # Create sig0 list
    sig0s = Sig0
    project = _phc.trimpathname(_os.getcwd())[1]
    for cl in clusters:
        file = open('{0}s/{0}s_{1}_mod{2}.sh'.format(cl, project, modn))
        lines = file.readlines()
        file.close()
        for item in Ms:
            M = item[0]
            exsig = item[1]
            for rm in exsig:
                _os.system('rm mod{0}/mod{0}*_sig{1}*_M{2}*.inp'.format(modn,
                sig0s[rm], M))
                print('# Deleted mod{0}/mod{0}*_sig{1}*_M{2}*.inp'.format(modn,
                sig0s[rm], M))
                _os.system('rm {3}s/mod{0}*_sig{1}*_M{2}*.{3}'.format(modn,
                sig0s[rm], M, cl))
                lines = [line for line in lines if (line.find('_sig{0}'.format(
                    sig0s[rm])) == -1 or line.find('_M{0}'.format(M)) == -1)]
        file = open('{0}s/{0}s_{1}_mod{2}.sh'.format(cl, project, modn), 'w')
        file.writelines(lines)
        file.close()
    # End prog
    return


def fsedList(fsedlist, param=True):
    """ Return the total of models and the parameters values in the fullsed list.

    The len of fsedlist is 9 (param=True) for the parametric case and 8
    to the VDD-ST one.

    The sequence is: M, ob(W), Z, H, sig, Rd, h, *n*, cos(i).

    It is assumed that all models have the same `observers` configuration."""
    nq = 9
    if not param:
        nq = 8
    listpar = [[] for i in range(nq)]
    nm = 0
    for sed in fsedlist:
        mod = BAmod(sed)
        if mod.param == param:
            nm += 1
            if mod.M not in listpar[0]:
                listpar[0].append(mod.M)
            if mod.ob not in listpar[1]:
                listpar[1].append(mod.ob)
            if mod.Z not in listpar[2]:
                listpar[2].append(mod.Z)
            if mod.H not in listpar[3]:
                listpar[3].append(mod.H)
            if mod.sig not in listpar[4]:
                listpar[4].append(mod.sig)
            if mod.Rd not in listpar[5]:
                listpar[5].append(mod.Rd)
            if mod.h not in listpar[6]:
                listpar[6].append(mod.h)
            if param:
                if mod.n not in listpar[7]:
                    listpar[7].append(mod.n)
            if listpar[-1] == []:
                sed2data = _hdt.readfullsed2(sed)
                listpar[-1] = list(sed2data[:, 0, 0])
    #
    for vals in listpar:
        vals.sort()
    return nm * len(listpar[-1]), listpar


def createBAsed(fsedlist, xdrpath, lbdarr, param=True, savetxt=False,
    ignorelum=False, pol=False, saveextra=None):
    """ Create the BeAtlas SED XDR release.

    WARNING: The file names must be in this format: 
    `mod01_PLn3.5_sig0.00_h072_Rd000.0_Be_M14.60_ob1.45_H0.77_Z0.014_bE_Ell`

    | The file structure:
    | -n_quantities, n_lbd, n_models,
    | -n_qt_vals1, n_qt_vals2, .. n_qt_valsn
    | -quantities values =  M, ob(W), Z, H, sig, Rd, h, *n*, cos(i).
    | -(Unique) lbd array
    | -Loop:
    |   *model values
    |   *model SED

    | Definitions:
    | -photospheric models: sig0 = 0.00
    | -Parametric disk model default (`param` == True)
    | -VDD-ST models: n excluded (alpha and R0 fixed. Confirm?)
    | -The flux will be given in ergs/s/um2/um. If ignorelum==True, the usual
    |   F_lbda/F_bol unit will be given.

    Since the grid is not symmetric, there is no index to jump directly to the
    desired model. So the suggestion is to use the index matrix, or read the
    file line by line until find the model (if exists).

    :Example: 

    def genxdr(xdrname='PL.xdr', param=True, pol=False):
        fs2l = glob('fullsed/*.sed2')
        print('# Using {0} as reference!'.format(fs2l[0]))
        lbdarr = hdt.readfullsed2(fs2l[0])
        lbdarr = lbdarr[0, :, 2]
        nm, listpar = bat.fsedList(fs2l)
        bat.createBAsed(fs2l, xdrname, lbdarr, param=param, savetxt=False, 
            pol=pol, saveextra=xdrname.replace('xdr', 'txt'))
        return

    genxdr(xdrname='Yudin_PL.xdr')
    """
    fsedlist.sort()
    nq = 9
    if not param:
        nq = 8
    nm, listpar = fsedList(fsedlist, param=param)
    header2 = []
    for vals in listpar:
        header2 += [len(vals)]
    nlb = len(lbdarr)
    header1 = [nq, nlb, nm]
    models = _np.zeros((nm, nlb))
    minfo = _np.zeros((nm, nq))
    k = 0
    iflx = 3
    if pol:
        iflx = 7
    for i in range(len(fsedlist)):
        mod = BAmod(fsedlist[i])
        # Select only `param` matching cases:
        if mod.param == param:
            sed2data = _hdt.readfullsed2(fsedlist[i])
            iL = 1.
            dist = 1/_np.sqrt(4 * _np.pi)
            if not ignorelum and not pol:
                j = fsedlist[i].find('fullsed_mod')
                # modn = fsedlist[i][j + 11:j + 13]
                modn = _re.match(r'.*mod(\d+)_', fsedlist[i]).group(1)
                log = fsedlist[i].replace('fullsed_mod', '../mod{0}/mod'.
                    format(modn)).replace('.sed2', '.log')
                if not _os.path.exists(log):
                    log = _glob(log.replace('../mod{0}/mod'.format(modn),
                    '../mod{0}/*mod'.format(modn)))
                    if len(log) >= 1:
                        log = log[0]
                    else:
                        raise LookupError('# No log file found for {0}'.
                            format(fsedlist[i]))
                f0 = open(log)
                lines = f0.readlines()
                f0.close()
                iL = _phc.fltTxtOccur('L =', lines, seq=2) * _phc.Lsun.cgs
                if saveextra is not None:
                    R_pole = _phc.fltTxtOccur('R_pole =', lines, seq=2)
                    Vrot = _phc.fltTxtOccur('Vrot', lines)
                    f0 = open(saveextra, 'a')
                    f0.writelines('{0}\t{1}\t{2}\n'.format(R_pole, Vrot, iL))
                    f0.close()
                dist = 10. * _phc.pc.cgs
            for j in range(header2[-1]):
                #  M, ob(W), Z, H, sig, Rd, h, *n*, cos(i).
                if param:
                    minfo[k * header2[-1] + j] = _np.array([ mod.M, mod.ob, 
                    mod.Z, mod.H, mod.sig, mod.Rd, mod.h, mod.n, 
                    listpar[-1][j] ]).astype(float)
                else:
                    minfo[k * header2[-1] + j] = _np.array([ mod.M, mod.ob, 
                    mod.Z, mod.H, mod.sig, mod.Rd, mod.h,
                    listpar[-1][j] ]).astype(float)
                if len(sed2data[j, :, 2]) != nlb:
                    models[k * header2[-1] + j] = _np.interp(lbdarr, 
                        sed2data[j, :, 2], sed2data[j, :, iflx]) * iL / 4 / \
                        _np.pi / dist**2
                else:
                    models[k * header2[-1] + j] = sed2data[j, :, iflx] * \
                        iL / 4 / _np.pi / dist**2
                if _np.sum(_np.isnan(models[k * header2[-1] + j])) > 0:
                    nans, x = _phc.nan_helper(models[k * header2[-1] + j]) 
                    models[k * header2[-1] + j][nans] = _np.interp(x(nans), 
                        x(~nans), models[k * header2[-1] + j][~nans])
            k += 1
    #
    f0 = open(xdrpath, 'wb')
    stfmt = '>{0}l'.format(3)
    f0.write(_struct.pack(stfmt, *header1))
    stfmt = '>{0}l'.format(nq)
    f0.write(_struct.pack(stfmt, *header2))
    for vals in listpar:
        stfmt = '>{0}f'.format(len(vals))
        f0.write(_struct.pack(stfmt, *_np.array(vals).astype(float)))
    stfmt = '>{0}f'.format(nlb)
    f0.write(_struct.pack(stfmt, *_np.array(lbdarr).astype(float)))
    for i in range(nm):
        stfmt = '>{0}f'.format(nq)
        f0.write(_struct.pack(stfmt, *minfo[i]))
        stfmt = '>{0}f'.format(nlb)
        f0.write(_struct.pack(stfmt, *_np.array(models[i]).astype(float)))
    f0.close()
    print('# XDR file {0} saved!'.format(xdrpath))

    if savetxt:
        f0 = open(xdrpath + '.txt', 'w')
        f0.writelines('{0} \n'.format(header1))
        f0.writelines('{0} \n'.format(header2))
        for vals in listpar:
            f0.writelines('{0} \n'.format(vals))
        f0.writelines('{0} \n'.format(lbdarr))
        for i in range(nm):
            f0.writelines('{0} \n'.format(minfo[i]))
            f0.writelines('{0} \n'.format(models[i]))
        f0.close()
        print('# TXT file {0} saved!'.format(xdrpath + '.txt'))
    return


def createXDRmap(maplist, xdrpath, refclass, lbdlim=None, npix=128):
    """ Is this possible? 

    ``lbdlim``: discard images that are out of this limits (``None`` keeps all;
    units of microns).

    Criteria:
    - all the images are converted the same size in pixels (``npix``).
    - only images with zoom = `renv` are used (checked by `*.log` file)
    - the pixel values are converted to flux units (erg/s/Ang)
    - the pixel scale is converted to length unit (at d = 10 pc)

    Example:
    maspp = maspp * 10/5  # transform the pixel scale for d = 5 pc
    """
    # maspp = milli arcsec per pixel
    raise NotImplementedError("To be done!")
    dval, renv, maspp, cubeimages = range(4)
    listpar = [dval, renv, maspp]
    return listpar, cubeimages


def createXDRsed(fsedlist, xdrpath, refclass, lbdarr, ignorelum=False, 
    pol=False):
    """ Create the generic SED XDR release.

    nob = (individual) number of observers
    listpar = parameters of each model
    nq = number of parameters
    nmod = number of models

    output units: 10**12 erg/s/Ang/cm2
    """
    if pol:
        ignorelum = True
    fsedlist.sort()
    ifact = 1.
    nq = len(refclass.vdict.keys())+1
    nlbd = len(lbdarr)
    listpar = _np.zeros( nq )
    models = _np.zeros( nlbd )
    for fs in fsedlist:
        print('# Processing {0}'.format(fs))
        m = refclass(fs)
        nob = m.get_nob()
        m.readfs2ob(pol=pol)
        if not ignorelum:
            ifact = m.get_lum() * _phc.Lsun.cgs / (4 * _np.pi * 
                (10 * _phc.pc.cgs)**2) * 1e8
        for j in range(nob):
            listpar = _np.vstack(( listpar, [ getattr(m, it) for it in 
                refclass.vdict.keys()]+[ m.obs[j] ] ))
            models = _np.vstack(( models, 
                _np.interp(lbdarr, m.arr_lbd, m.arr_flx[j])*ifact ))

    models = models[1:]
    listpar = listpar[1:]
    nmod = len(models)

    chk = []
    for i in range( nq ):
        if len(_np.unique(listpar[:, i])) == 1:
            chk.append(i)
    if len(chk) > 0:
        listpar = _np.delete(listpar, chk, axis=1)
    nq = len(listpar[0])
    #
    intervals = _np.zeros(( nq, 2 ))
    intervals[:, 0] = _np.min(listpar, axis=0)
    intervals[:, 1] = _np.max(listpar, axis=0)
    #
    f0 = open(xdrpath, 'wb')
    stfmt = '>{0}l'.format(3)
    f0.write( _struct.pack(stfmt, nq, nlbd, nmod) )
    stfmt = '>{0}f'.format(nq*2)
    f0.write( _struct.pack(stfmt, *intervals.flatten()) )
    stfmt = '>{0}f'.format(nlbd)
    f0.write( _struct.pack(stfmt, *lbdarr))
    stfmt = '>{0}f'.format(nq*nmod)
    f0.write( _struct.pack(stfmt, *listpar.flatten()) )
    stfmt = '>{0}f'.format(nlbd*nmod)
    f0.write( _struct.pack(stfmt, *models.flatten()) )
    f0.close()
    print('# XDR file {0} saved!'.format(xdrpath))
    return


def readXDRsed(xdrpath, quiet=False):
    """ Doc
    """
    ixdr = 0
    f = open(xdrpath, 'rb').read()
    ixdr, ninfo = _phc.readpck(3, 'l', ixdr, f)
    nq, nlbd, nm = ninfo
    ixdr, intervals = _phc.readpck(nq*2, 'f', ixdr, f)
    ixdr, lbdarr = _phc.readpck(nlbd, 'f', ixdr, f)
    ixdr, listpar = _phc.readpck(nq*nm, 'f', ixdr, f)
    ixdr, models = _phc.readpck(nlbd*nm, 'f', ixdr, f)
    #
    if ixdr == len(f):
        if not quiet:
            print('# XDR {0} completely read!'.format(xdrpath))
    else:
        _warn.warn('# XDR {0} not completely read!\n# length '
            'difference is {1} /4'.format(xdrpath), (len(f)-ixdr) )
    # 
    return ( ninfo, intervals.reshape((nq, 2)), lbdarr, 
        listpar.reshape((nm, nq)), models.reshape((nm, nlbd)) )


def readBAsed(xdrpath, quiet=False):
    """ Read the BeAtlas SED release.

    | Definitions:
    | -photospheric models: sig0 (and other quantities) == 0.00
    | -Parametric disk model default (`param` == True)
    | -VDD-ST models: n excluded (alpha and R0 fixed. Confirm?)
    | -The models flux are given in ergs/s/cm2/um. If ignorelum==True in the
    |   XDR creation, F_lbda/F_bol unit will be given.

    INPUT: xdrpath

    | OUTPUT: listpar, lbdarr, minfo, models 
    | (list of mods parameters, lambda array (um), mods index, mods flux)
    """
    f = open(xdrpath, 'rb').read()
    ixdr = 0
    # 
    npxs = 3
    upck = '>{0}l'.format(npxs)
    header = _np.array(_struct.unpack(upck, f[ixdr:ixdr + npxs * 4]) )
    ixdr += npxs * 4
    nq, nlb, nm = header
    # 
    npxs = nq
    upck = '>{0}l'.format(npxs)
    header = _np.array(_struct.unpack(upck, f[ixdr:ixdr + npxs * 4]) )
    ixdr += npxs * 4
    # 
    listpar = [[] for i in range(nq)]
    for i in range(nq):
        npxs = header[i]
        upck = '>{0}f'.format(npxs)
        listpar[i] = _np.array(_struct.unpack(upck, f[ixdr:ixdr + npxs * 4]) )
        ixdr += npxs * 4
    # 
    npxs = nlb
    upck = '>{0}f'.format(npxs)
    lbdarr = _np.array(_struct.unpack(upck, f[ixdr:ixdr + npxs * 4]) )
    ixdr += npxs * 4
    # 
    npxs = nm * (nq + nlb)
    upck = '>{0}f'.format(npxs)
    models = _np.array(_struct.unpack(upck, f[ixdr:ixdr + npxs * 4]) )
    ixdr += npxs * 4
    models = models.reshape((nm, -1))
    # this will check if the XDR is finished.
    if ixdr == len(f):
        if not quiet:
            print('# XDR {0} completely read!'.format(xdrpath))
    else:
        _warn.warn('# XDR {0} not completely read!\n# length '
            'difference is {0}'.format(xdrpath, (len(f)-ixdr)/4) )
    # 
    return listpar, lbdarr, models[:, 0:nq], models[:, nq:]


def parnorm(dvals, vmax, vmin_non0, issig0=True, s_non0=0):
    r""" Converts density in normalized range [0-1], and vice-versa.

    If ``issig0``, treats ``r01`` as :math:`\Sigma_0`; otherwise, use it
    as [0-1] value.
    """
    dvals = _np.array(dvals)
    if vmin_non0 <= 0 or _np.min(dvals) < 0:
        raise ValueError('`vmin_non0` > 0 and `dvals` >=0 must be True')
    if s_non0 >= 1:
        raise ValueError('`s_non0` <1 must be True')
    if s_non0 > 0:
        vmin_non0 /= (vmax/vmin_non0)**((1/s_non0-1)**-1.)
    if issig0:
        dvals[_np.where(dvals < vmin_non0)] = vmin_non0
        return _np.log(dvals/vmin_non0)/_np.log(vmax/vmin_non0)
    else:
        return _np.exp(dvals*_np.log(vmax/vmin_non0))*vmin_non0


def densBAnorm(r01, M, issig0=True):
    r""" Converts density in normalized range [0-1], and vice-versa for the 
    BeAtlas.

    If ``issig0``, treats ``r01`` as :math:`\Sigma_0`; otherwise, use it
    as [0-1] value.

    """
    if M < 3.8 or M > 14.6:
        raise ValueError('# Wrong M at bat.normdens() !')
    vmin = 0.02/2.5
    r01 = _np.array(r01)
    if issig0:
        r01[_np.where(r01 < vmin)] = vmin

    # Completo, convergido etapa1, com DEPENDENCIA do vizinho inferior
    x = [4.2, 4.8, 5.5, 6.4, 7.7, 8.6, 9.6, 10.8]
    y = [0.05, 0.12, 0.28, 0.28, 0.68, 1.65, 1.65, 4]

    # Completo, convergido etapa2, com DEPENDENCIA do vizinho inferior
    # x = [4.2, 4.8, 5.5, 6.4, 7.7, 8.6, 9.6, 10.8]
    # y = [0.05, 0.12, 0.28, 0.40, 0.68, 1.65, 2.46, 4]

    d = 7
    x = _np.array(x) + 10**-d
    vmax = _np.interp(M, x, y)
    if issig0:
        return _np.round( _np.log(r01/vmin)/_np.log(vmax/vmin), d)
    else:
        return _np.round( _np.exp(r01*_np.log(vmax/vmin))*vmin, d)


def interpolBA2(params, ctrlarr, minfo, models):
    ctrlarr[_np.isnan(ctrlarr)] = params
    return _griddata(minfo, models, ctrlarr)


def interpolBA(params, ctrlarr, lparams, minfo, models, param=True):
    """ Interpola os `modelos` para os parametros `params` 

    | -params = from emcee minimization
    | -ctrlarr = the fixed value of M, ob(W), Z, H, sig, Rd, h, *n*, cos(i).
    |            If it is not fixed, use np.NaN.
    | -Parametric disk model default (`param` == True).

    This function always returns a valid result (i.e., extrapolations from the
    nearest values are always on).

    If it is a 'Non-squared grid' (asymmetric), it will return a zero array if
    a given model is not found.
    """
    nq = 9
    if not param:
        nq = 8
    if len(ctrlarr) != nq:
        raise ValueError('# Wrong ctrlarr format!!')
    params = params[:_np.sum(_np.isnan(ctrlarr))]
    nlb = len(models[0])
    outmodels = _np.empty((2**len(params), nlb))
    mod = BAmod('')
    parlims = _np.zeros((len(params), 2))
    j = 0
    for i in range(nq):
        if ctrlarr[i] is _np.NaN:
            parlims[j] = [_phc.find_nearest(lparams[i], params[j], 
            bigger=False), _phc.find_nearest(lparams[i], params[j], 
            bigger=True)]
            j += 1
    j = 0
    for prod in _product(*parlims):
        allpars = _np.array(ctrlarr)
        idx = _np.isnan(allpars)
        allpars[idx] = prod
        mod.build(allpars, lparams)
        idx = mod.getidx(minfo)
        if _np.sum(idx) == 0:
            return _np.zeros(nlb)
        outmodels[j] = models[idx]
        j += 1
    X0 = parlims[:, 0]
    X1 = parlims[:, 1]
    return _phc.interLinND(params, X0, X1, outmodels)


# MAIN ###
if __name__ == "__main__":
    pass
