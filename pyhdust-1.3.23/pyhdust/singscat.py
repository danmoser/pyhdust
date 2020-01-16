# -*- coding:utf-8 -*-

"""PyHdust *singscat* module: Single Scatering Dumbbell+disk model. 

This release contains a simplified version of the one applied to Sigma Ori E in 
Carciofi+2013.

Conditions:
- Calculations for a given observer position: observer at z = +infinity
- Single scattering by electrons (Thomson)
- Star is a point-source (*or at least, a spherical emitting source)
- :math:`\theta` is latitude. 0 = North pole; 180 = South pole.
- :math:`\phi` is longitude
- Spherical coordinates seguence (:math:`r`, :math:`\phi`, :math:`\theta`)

Remeber:
- phi=0 is x direction (perpendicular to the observer)
- (x, y) are not the same for (N, -E) astronomical images

:license: GNU GPL v3.0  https://github.com/danmoser/pyhdust/blob/master/LICENSE 
"""
from __future__ import print_function
# import os as _os
import numpy as _np
# import time as _time
# import pyhdust as _hdt
# import pyhdust.jdcal as _jdcal
# import pyhdust.poltools as _polt
import pyhdust.phc as _phc
# import pyhdust.triangle as _triangle
import warnings as _warn

try:
    import matplotlib.pyplot as _plt
    import matplotlib.gridspec as _gridspec
    from mpl_toolkits.mplot3d import Axes3D as _Axes3D
    # import emcee as _emcee
    # from scipy.stats import percentileofscore as _perct
    # from scipy import optimize as _optimize
    from scipy import interpolate as _interpolate

except ImportError:
    _warn.warn('matplotlib, scipy and/or emcee module not installed!!!')

__author__ = "Daniel Moser"
__email__ = "dmfaes@gmail.com"


# sing scat mod base
def sph_hyper_rej(res):
    """ Hypercube rejection method
    """
    x0 = _np.linspace(-1, 1, res)
    y0 = _np.linspace(-1, 1, res)
    z0 = _np.linspace(-1, 1, res)
    xx0, yy0, zz0 = _np.meshgrid(x0, y0, z0)
    rr = _np.sqrt(xx0**2+yy0**2+zz0**2)
    idx = _np.where(rr <= 1) 
    return xx0[idx], yy0[idx], zz0[idx]


def blobs_coords(db=2, rb=1./3, res=3, Req=1., phi0=0, beta=0.):
    """ 
    :param db: distance from stellar center to blob center (Req)
    :param rb: radius of the blob (sphere; Req) 
    :param res: **linear** resolution used in the *Hypercube rejection* method.
    :param Req: stellar equatorial radius (in cm)
    :param phi0: central position of the blob (in rad) 
    """
    x, y, z = sph_hyper_rej(res)
    x, y, z = _phc.cart_rot(x, y, z, ang_xy=0., ang_yz=beta, ang_zx=0.)
    r, phi, th = _phc.cart2sph((x*rb+db), y*rb, z*rb)
    phi += phi0
    return r*Req, phi, th


def cil_slice_rej(res, dd, rd):
    """ NO *z* dimension
    """
    x0 = _np.linspace(-1, 1, res)*(rd+dd)
    y0 = _np.linspace(-1, 1, res)*(rd+dd)
    xx0, yy0 = _np.meshgrid(x0, y0)
    rr = _np.sqrt(xx0**2+yy0**2)
    idx = _np.where((rr <= dd+rd) & (rr >= dd-rd))
    return xx0[idx], yy0[idx], _np.zeros(_np.shape(xx0[idx]))


def disk_coords(dd=2, rd=1./3, res=11., Req=1., beta=0., phi0=0.):
    """  Cilinder slice rejection method

    :param beta: obliquity, in rad
    :param phi0: **only** effective when beta != 0.
    """
    x, y, z = cil_slice_rej(res, dd, rd)
    # The disk is invariant with :math:`\delta\phi`. So, instead of appling
    #  phi0=_np.pi/2 and ang_yz, one could directly use ang_zx.
    # x, y, z = _phc.cart_rot(x, y, z, ang_xy=0., ang_yz=0., ang_zx=beta)
    x, y, z = _phc.cart_rot(x, y, z, ang_xy=0., ang_yz=beta, ang_zx=0.)
    r, phi, th = _phc.cart2sph(x, y, z)
    phi += phi0
    return r*Req, phi, th


def idx_phi_blobdisk_coords(phib, rd, phid, thd):
    # idx = _np.where( (_np.cos(phid) < _np.min(_np.cos(phib))) | 
    #     (_np.cos(phid) > _np.max(_np.cos(phib))) )
    # return rd[idx], phid[idx], thd[idx]
    # Function standard: always between 0 to 2pi
    # if len(radcut)!= 2:
    #     raise ValueError("len(`radcut`) != 2: {}".format(radcut))
    # elif _np.min(radcut) == _np.max(radcut):
    #     raise ValueError("delta-radcut > 0: {}".format(radcut))
    tpi = 2*_np.pi
    phib = _np.array([_np.min(phib), _np.max(phib)])
    phib = _np.where(phib < 0, phib+tpi, phib)
    phib = _np.where(phib > tpi, phib-tpi, phib)
    if _np.diff(phib) > _np.pi:
        phib = phib[::-1]
    phid = _np.where(phid < 0, phid+tpi, phid)
    phid = _np.where(phid > tpi, phid-tpi, phid)
    if phib[0] < phib[1]:
        # 350, 10
        idx = _np.where( (phid < phib[0]) | (phid > phib[1]) )
    else:
        idx = _np.where( (phid < phib[0]) & (phid > phib[1]) )
    return rd[idx], phid[idx], thd[idx]


def idx_phi_disklimit_coords(radcut, rd, phid, thd):
    """ Enter the coordinates of a 2*pi disk and cuts it with the limits of 
    radcut

    :param radcut: vector as [-np.pi/6, np.pi/6]
    """
    # Function standard: always between 0 to 2pi
    if len(radcut)!= 2:
        raise ValueError("len(`radcut`) != 2: {}".format(radcut))
    elif _np.min(radcut) == _np.max(radcut):
        raise ValueError("delta-radcut > 0: {}".format(radcut))
    tpi = 2*_np.pi
    radcut = _np.array(radcut)
    radcut = _np.where(radcut < 0, radcut+tpi, radcut)
    radcut = _np.where(radcut > tpi, radcut-tpi, radcut)
    phid = _np.where(phid < 0, phid+tpi, phid)
    phid = _np.where(phid > tpi, phid-tpi, phid)
    if radcut[0] < radcut[1]:
        # 350, 10
        idx = _np.where( (phid < radcut[0]) | (phid > radcut[1]) )
    else:
        idx = _np.where( (phid < radcut[0]) & (phid > radcut[1]) )
    return rd[idx], phid[idx], thd[idx]


def stokes(r, phi, th, irad, dV, ne, Req, occult=True):
    """ 
    :param irad: in rad
    :param dV: in cm3

    TODO: intensity calc.

    output in fraction (multiple by 100 to have percentage
    """
    x, y, z = _phc.sph2cart(r, phi, th)
    x, y, z = _phc.cart_rot(x, y, z, 0., irad, 0.)
    norm = 3*_phc.sigT.cgs/8/_np.pi*(1./r)**2*_np.cos(_np.arctan(Req/r))*dV*ne
    cos_chi = _np.cos(th)*_np.cos(irad) + _np.cos(phi-0)*_np.sin(th)*_np.sin(
        irad)
    sI = norm/2. * ( 1 + cos_chi**2 )
    sQ = norm/2. * ( (x / r)**2 - (y / r)**2 )
    sU = norm/2. * ( 2. * x * y / r**2 )
    if occult is True: 
        ind = _np.where((x**2 + y**2 < Req**2) & (z < 0))
        sI[ind] = 0.
        sQ[ind] = 0.
        sU[ind] = 0.
        # ind = (x**2 + y**2 < Req**2) & (z > 0)
        # sI[ind] = sI[ind]
    return sI, sQ, sU, _np.zeros(len(r))


# single scat higher level
def blobsdiskmodel_geo(phi0, Req=_phc.Rsun.cgs, rb=1/3., db=2., bres=3, 
        nb=8e11, rd=1/3., dd=2., dres=11, H=.1, beta=0., nd=8e11, 
        overlap=False):
    """ Doc

    :param rb: in Req units
    :param overlap: disk and blob overlap each other? Default is `False`.
    Warning: only :math:`\phi` direction is checked (ie., coplanar overlap).
    :param beta: obliquity, in rad
    """
    r1, phi1, th1 = blobs_coords(phi0=phi0, rb=rb, Req=Req, db=db, res=bres, 
        beta=beta)
    xb1, yb1, zb1 = _phc.sph2cart(r1, phi1, th1)
    r2, phi2, th2 = blobs_coords(phi0=phi0+_np.pi, rb=rb, Req=Req, db=db, 
        res=bres, beta=beta)
    xb2, yb2, zb2 = _phc.sph2cart(r2, phi2, th2)
    Vb = 4./3*_np.pi*(rb*Req)**3 
    dVb = _np.zeros(len(r1))+Vb/len(r1)
    nb = _np.zeros(len(r1))+nb

    rD, phid, thd = disk_coords(rd=rd, Req=Req, dd=dd, res=dres, beta=beta,
        phi0=phi0)
    area_d = _np.pi*( (dd+rd)**2 - (dd-rd)**2 )
    Vd = area_d*H*Req**3
    if not overlap:
        # Remove the overlap of the disk over the blobs
        rD, phid, thd = idx_phi_blobdisk_coords(phi1, rD, phid, thd)
        rD, phid, thd = idx_phi_blobdisk_coords(phi2, rD, phid, thd)
        area_bd = _np.pi*rd**2
        Vd = (area_d-2*area_bd)*H*Req**3
    dVd = _np.zeros(len(rD))+Vd/len(rD)
    nd = _np.zeros(len(rD))+nd
    xd, yd, zd = _phc.sph2cart(rD, phid, thd)

    return (_np.concatenate((xb1, xb2, xd)), _np.concatenate((yb1, yb2, yd)), 
        _np.concatenate((zb1, zb2, zd)), _np.concatenate((r1, r2, rD)),
        _np.concatenate((phi1, phi2, phid)), _np.concatenate((th1, th2, thd)),
        _np.concatenate((dVb, dVb, dVd)), _np.concatenate((nb, nb, nd)))


def blobsmodel_geo(phi0, Req=_phc.Rsun.cgs, rb=1/3., db=2., bres=3, 
        nb=8e11, beta=0.):
    """ Doc

    :param overlap: disk and blob overlap each other? Default is `False`.
    Warning: only :math:`\phi` direction is checked (ie., coplanar overlap).
    :param beta: obliquity, in rad
    """
    r1, phi1, th1 = blobs_coords(phi0=phi0, rb=rb, Req=Req, db=db, res=bres, 
        beta=beta)
    xb1, yb1, zb1 = _phc.sph2cart(r1, phi1, th1)
    r2, phi2, th2 = blobs_coords(phi0=phi0+_np.pi, rb=rb, Req=Req, db=db, 
        res=bres, beta=beta)
    xb2, yb2, zb2 = _phc.sph2cart(r2, phi2, th2)
    Vb = 4./3*_np.pi*(rb*Req)**3 
    dVb = _np.zeros(len(r1))+Vb/len(r1)
    nb = _np.zeros(len(r1))+nb

    return (_np.concatenate((xb1, xb2)), _np.concatenate((yb1, yb2)), 
        _np.concatenate((zb1, zb2)), _np.concatenate((r1, r2)),
        _np.concatenate((phi1, phi2)), _np.concatenate((th1, th2)),
        _np.concatenate((dVb, dVb)), _np.concatenate((nb, nb)))


def diskmodel_geo(phi0, Req=_phc.Rsun.cgs, rd=1/3., dd=2., dres=11, H=.1, 
    beta=0., nd=8e11, radcut=[]):
    """ Doc

    :param rd: in Req units
    :param beta: obliquity, in rad
    :param radcut: vector as [[1*np.pi/4, 3*np.pi/4], [5*np.pi/4, 7*np.pi/4]]. 
        It must be around pi/2 and 3pi/2 for phase consistency. ``[]`` returns 
        a full (2*pi) disk (Radian units).
    """
    rD, phid, thd = disk_coords(rd=rd, Req=Req, dd=dd, res=dres, beta=beta,
        phi0=phi0)
    total_phi_disk = 2*_np.pi
    for philim in radcut:
        # Remove an interval from the disk model
        philim = _np.array(philim)+phi0
        rD, phid, thd = idx_phi_disklimit_coords(philim, rD, phid, thd)
        total_phi_disk -= _np.max(philim)-_np.min(philim)
    area_d = (total_phi_disk/2.)*( (dd+rd)**2 - (dd-rd)**2 )    
    Vd = area_d*H*Req**3
    dVd = _np.zeros(len(rD))+Vd/len(rD)
    nd = _np.zeros(len(rD))+nd
    xd, yd, zd = _phc.sph2cart(rD, phid, thd)

    return xd, yd, zd, rD, phid, thd, dVd, nd


def blobs_cicle(phi_ar, irad, Req=_phc.Rsun.cgs, rb=1/3., db=2., bres=3, 
        nb=8e11):
    """ Calculates the stokes parameters of blobsmodel over a full cycle

    :param phi_ar: phi array (in rad; [0, 2pi))
    """
    stk = _np.zeros((len(phi_ar), 4))
    for i in range(len(phi_ar)):
        phi0 = phi_ar[i]
        x, y, z, r, phi, th, dV, ne = blobsmodel_geo(phi0, Req=Req, rb=rb, 
            db=db, bres=bres)
        istk = stokes(r, phi, th, dV=dV, ne=ne, irad=irad, Req=Req)
        stk[i] = _np.array(istk).sum(axis=1)
    return stk.T


def blobsdisk_cicle(phi_ar, irad, Req=_phc.Rsun.cgs, rb=1/3., db=2., bres=3, 
        nb=8e11, rd=1/3., dd=2., dres=11, H=.1, beta=0., nd=8e11, 
        overlap=False):
    """ Calculates the stokes parameters of blobs+disk model over a full cycle

    :param phi_ar: phi array (in rad; [0, 2pi))
    """
    stk = _np.zeros((len(phi_ar), 4))
    for i in range(len(phi_ar)):
        phi0 = phi_ar[i]
        x, y, z, r, phi, th, dV, ne = blobsdiskmodel_geo(phi0, Req=Req, rb=rb, 
            db=db, bres=bres, nb=nb, rd=rd, dd=dd, dres=dres, H=H, beta=beta, 
            nd=nd, overlap=overlap)
        istk = stokes(r, phi, th, dV=dV, ne=ne, irad=irad, Req=Req)
        stk[i] = _np.array(istk).sum(axis=1)
    return stk.T


def disk_cicle(phi_ar, irad, Req=_phc.Rsun.cgs, rd=1/3., dd=2., dres=11, 
        H=.1, beta=0., nd=8e11, radcut=[]):
    """ Calculates the stokes parameters of disk model over a full cycle

    :param phi_ar: phi array (in rad; [0, 2pi))
    :param radcut: vector as [[-np.pi/4, np.pi/4], [-np.pi/4+np.pi, 
        np.pi/4+np.pi]]. ``[]`` returns a full (2*pi) disk (Radian units).
    """
    stk = _np.zeros((len(phi_ar), 4))
    for i in range(len(phi_ar)):
        phi0 = phi_ar[i]
        x, y, z, r, phi, th, dV, ne = diskmodel_geo(phi0, Req=Req, rd=rd, 
            dd=dd, dres=dres, H=H, beta=beta, nd=nd, radcut=radcut)
        istk = stokes(r, phi, th, dV=dV, ne=ne, irad=irad, Req=Req)
        stk[i] = _np.array(istk).sum(axis=1)
    return stk.T


# general 
def angQU(Q, U, filter=True):
    """ Calculate [Q, U] angles with filters* (avg-180 < degs > avg+180) """
    Q = _np.array(Q)
    U = _np.array(U)
    ind = _np.where(Q == 0)
    Q[ind] = 1e-34
    ang = _np.arctan(U / Q)
    #
    ind = _np.where(Q <= 0.)
    ang[ind] = ang[ind] + _np.pi
    ind = _np.where((Q > 0) & (U < 0))
    ang[ind] = ang[ind] + 2 * _np.pi
    ang = ang / 2.
    # ind = _np.where(ang >= _np.pi)
    # ang[ind] = ang[ind] - _np.pi
    if filter:
        avg = _np.median(ang)
        avg = _phc.find_nearest(
            [0, _np.pi / 4, _np.pi / 2, _np.pi * 3. / 4], avg)
        ind = _np.where((ang - avg) > 1. / 2 * _np.pi)
        ang[ind] = ang[ind] - _np.pi
        ind = _np.where((ang - avg) < -1. / 2 * _np.pi)
        ang[ind] = ang[ind] + _np.pi
    return ang


def mod2obs(Qmod, Umod, Qis, Uis, ths):
    """ MODEL TO OBSERV

    :param ths: theta sky in rad
    :param Qis: IS Q
    :param Uis: IS U
    """
    ang = angQU(Qmod, Umod)
    Qobc = _np.sqrt(Qmod**2 + Umod**2) * _np.cos(2 * (ang + ths))
    Uobc = _np.sqrt(Qmod**2 + Umod**2) * _np.sin(2 * (ang + ths))         
    Qobc = Qobc + Qis
    Uobc = Uobc + Uis
    return Qobc, Uobc


def obs2mod(Qobs, Uobs, Qis, Uis, ths):
    """ ### MODEL TO OBSERV ### 

    :param ths: theta sky"""
    Qcob = Qobs - Qis
    Ucob = Uobs - Uis
    ang = angQU(Qcob, Ucob)
    Qcob = _np.sqrt(Qcob**2 + Ucob**2) * _np.cos(2 * (ang - ths))
    Ucob = _np.sqrt(Qcob**2 + Ucob**2) * _np.sin(2 * (ang - ths))         
    return Qcob, Ucob


def calc_phase(MJDar, period, MJD0=0.):
    """ Calc phase [0,1) 
    """
    phase = _np.modf( (_np.array(MJDar)-MJD0)/period )[0]
    return _np.where(phase>0, phase, phase+1)


# beacon
def loadpol(txt, old=False, p_sigP=0):
    """ Load polarization txt file. 

    :param p_sigP: if > 0, filter values P/sigP
    """
    dtb = _np.genfromtxt(txt, dtype=str)
    if not old:
        dtb = _np.core.records.fromarrays(dtb.transpose(), 
            names='MJD,night,ccd,filt,calc,stdstars,dth,sigdth,P,Q,U,th,sigP,'
                'sigth,outfile,star,flag,tags', 
            formats='f8,{0},{0},{0},f8,{0},f8,f8,f8,f8,f8,f8,f8,f8,{0},{0},{0}'
                ',{0}'.format(dtb.dtype))
    else:
        dtb = _np.core.records.fromarrays(dtb.transpose(), names='MJD,night,'
            'filt,calc,stdstars,dth,devdth,P,Q,U,th,sigP,sigth', 
            formats='f8,{0},{0},f8,{0},f8,f8,f8,f8,f8,f8,f8,f8'.format(
                dtb.dtype))
    if p_sigP > 0:
        idx = _np.where(dtb['P']/dtb['sigP'] >= p_sigP)
        return dtb[idx]
    return dtb


# plotting
def plot_coords(xb1, yb1, zb1, Req, slim=3):
    """ :param slim: Scale limit [-slim, slim] 
    """
    fig = _plt.figure()
    lins, cols = (2, 2)
    gs = _gridspec.GridSpec(lins, cols)

    ax = _plt.subplot(gs[0, 0], projection='3d')  # first line, first col
    ax10 = _plt.subplot(gs[1, 0])  # second line, first col
    ax01 = _plt.subplot(gs[0, 1])  
    ax11 = _plt.subplot(gs[1, 1])

    ax.scatter(xb1/Req, yb1/Req, zb1/Req)
    ax.set_xlim3d(-slim, slim)
    ax.set_ylim3d(-slim, slim)
    ax.set_zlim3d(-slim, slim)

    ax10.scatter(yb1/Req, zb1/Req)
    ax10.set_xlabel('y')
    ax10.set_ylabel('z')
    ax10.set_xlim(-slim, slim)
    ax10.set_ylim(-slim, slim)

    ax01.scatter(xb1/Req, yb1/Req)
    ax01.set_xlabel('x')
    ax01.set_ylabel('y')
    ax01.set_xlim(-slim, slim)
    ax01.set_ylim(-slim, slim)

    ax11.scatter(xb1/Req, zb1/Req)
    ax11.set_xlabel('x')
    ax11.set_ylabel('z')
    ax11.set_xlim(-slim, slim)
    ax11.set_ylim(-slim, slim)
    return fig, [ax, ax01, ax10, ax11]


def plot_QU_gd(x, y, Q, U, irad, Req):
    """ using griddata 
    """
    fig = _plt.figure()
    lins, cols = (1, 2)
    gs = _gridspec.GridSpec(lins, cols)

    axq = _plt.subplot(gs[0, 0])  
    axu = _plt.subplot(gs[0, 1])  

    xmin = _np.min(x)/Req
    xmax = _np.max(x)/Req
    ymin = _np.min(y)/Req
    ymax = _np.max(y)/Req
    xx, yy = _np.meshgrid(_np.linspace(xmin, xmax, 32), 
        _np.linspace(ymin, ymax, 32)[::-1])
    yo = y*_np.cos(irad)
    q = _interpolate.griddata( _np.array([x, yo]).T/Req, Q, 
        _np.array([xx.flatten(), yy.flatten()]).T )
    u = _interpolate.griddata( _np.array([x, yo]).T/Req, U, 
        _np.array([xx.flatten(), yy.flatten()]).T )

    axq.imshow(q.reshape(32, 32), origin='lower', extent=[xmin, xmax, 
        ymin, ymax])
    axu.imshow(u.reshape(32, 32), origin='lower', extent=[xmin, xmax, 
        ymin, ymax])
    return fig, [axq, axu]


def plot_QU_bc(x, y, Q, U, irad, Req):
    """ using baricenter 

    TODO: discover why it is not working
    """ 
    fig2 = _plt.figure()
    lins, cols = (1, 2)
    gs = _gridspec.GridSpec(lins, cols)

    axq2 = _plt.subplot(gs[0, 0])  
    axu2 = _plt.subplot(gs[0, 1])  

    xmin = _np.min(x)/Req
    xmax = _np.max(x)/Req
    ymin = _np.min(y)/Req
    ymax = _np.max(y)/Req
    axq2.imshow(_phc.baricent_map(x, y, Q, res=2*32, yfact=_np.cos(irad)), 
        origin='lower', extent=[xmin, xmax, ymin, ymax])
    axu2.imshow(_phc.baricent_map(x, y, U, res=2*32, yfact=_np.cos(irad)), 
        origin='lower', extent=[xmin, xmax, ymin, ymax])
    return fig2, [axq2, axu2]


# MAIN ###
if __name__ == "__main__":
    pass
