# -*- coding:utf-8 -*-

"""PyHdust *singscat* module: Single Scatering Dumbbell+disk model. 

This release contains a simplified version of the one applied to Sigma Ori E in 
Carciofi+2013.

| VARIABLES:
|    BLOB
|    - i
|    - Rstar
|    - Blob diam.
|    - Blob dist. (center)
|    - n_e
|
|    DISK
|    - n_d
|    - Disk length (assumed centered at blob dist.)
|    - Disk height
|    - Angle between stellar rotation and magnetic field (disk inclination)
|
|    OBSERVATIONAL
|    - Pis = Qis + Uis
|    - Theta_sky
|    - Delta phase from photometry
|
|    CALCULATIONS
|    - n0 = blob resolution (n0**3)
|    - phases
|    - dh = disk height step
|    - ddr = disk radius steps
|    - dphi = disk angle step

:license: GNU GPL v3.0  https://github.com/danmoser/pyhdust/blob/master/LICENSE 
"""
from __future__ import print_function
import os as _os
import numpy as _np
import time as _time
import pyhdust as _hdt
import pyhdust.jdcal as _jdcal
import pyhdust.poltools as _polt
import pyhdust.phc as _phc
import pyhdust.triangle as _triangle
import warnings as _warn

try:
    import matplotlib.pyplot as _plt
    import emcee as _emcee
    from scipy.stats import percentileofscore as _perct
except ImportError:
    _warn.warn('matplotlib, scipy and/or emcee module not installed!!!')

__author__ = "Daniel Moser"
__email__ = "dmfaes@gmail.com"

# PHYSICAL CONSTANTS VARS ###
_sigT = _phc.sigT.cgs  # cm^2 = Thomson cross section


def diskcoords(pdk):
    """ ### DISK geom. ### """
    # Qis,Uis,ne,phi0,ths,iang,fact = pob
    rdi, rdf, H, alpha, ned, dh, ddr, dphi = pdk

    dvarr = (rdf - rdi) / ddr
    varrD = _np.zeros((ddr, dh, dphi))
    thD = _np.zeros((ddr, dh, dphi))
    phiD = _np.zeros((ddr, dh, dphi))
    for i in range(ddr):
        for j in range(dh):
            for k in range(dphi):
                varrD[i][j][k] = rdi + dvarr * (_np.arange(ddr)[i] + .5)
                thD[i][j][k] = 90. * _np.pi / 180  # TODO
                phiD[i][j][k] = _np.linspace(
                    0, (2 - 1. / dphi * 2) * _np.pi, dphi)[k]

    varrD = varrD[~_np.isnan(varrD)]
    phiD = phiD[~_np.isnan(phiD)]
    thD = thD[~_np.isnan(thD)]

    (xD, yD, zD) = _phc.sph2cart(varrD, thD, phiD)
    (xD, yD, zD) = _phc.cart_rot(xD, yD, zD, 0., 0., alpha)
    (varrD, thD, phiD) = _phc.cart2sph(xD, yD, zD)

    return varrD, thD, phiD


def stokesD(varrD, thD, phiD, pst, pob, pdk):
    """ ### STOKES CALC given DISK ### """
    rs, diamb, distb, n0, occult = pst
    Qis, Uis, ne, phi0, ths, iang, fact = pob
    rdi, rdf, H, alpha, ned, dh, ddr, dphi = pdk

    (x, y, z) = _phc.sph2cart(varrD, thD, phiD)
    (x, y, z) = _phc.cart_rot(x, y, z, 0., iang, 0.)

    dV = H * varrD * 2 * _np.pi / dphi * H / dh  # TODO
    rl = (rdf - rdi) / ddr
    P0 = 3. / 8 * (ned * _sigT * rl) * 1. / 2 * (rs / varrD)**2 * _np.sqrt(1 -
        (rs / varrD)**2) * (2 * _np.pi / dphi * varrD * H / dh) / (varrD**2)
    cosX = _np.sin(iang) * _np.sin(thD) * _np.cos(phiD) + \
        _np.cos(iang) * _np.cos(thD)
    # Ii = _np.zeros(len(varr))+1./(len(varr)*2)
    Ii = P0 * (1 + cosX**2)
    Qi = P0 * ( (x / varrD)**2 - (y / varrD)**2 )
    Ui = P0 * ( 2. * x * y / varrD**2 )
    Pi = _np.sqrt(Qi**2 + Ui**2)
    # Occultation ##
    if occult is True: 
        ind = (x**2 + y**2 < rs**2) & (z < 0)
        Qi[ind] = 0
        Ui[ind] = 0
        Pi[ind] = 0
        Ii[ind] = 0.
        ind = (x**2 + y**2 < rs**2) & (z > 0)
        Ii[ind] = Ii[ind] - (rs / varrD[ind])**2 * (_np.exp(-ned * _sigT * 
            rl)) * (2 * _np.pi / dphi * varrD[ind] * H / dh) / \
            (varrD[ind]**2) * 1 / 4. / _np.pi
    return Ii.sum(), fact * Pi.sum(), fact * Qi.sum(), fact * Ui.sum()    


def modcycleF(varr, th, phi_0, n3, varrD, thD, phiD_0, phiT, pst, pob, pdk):
    """ modcycleF """
    # pst = [rs,diamb,distb,n0,occult]
    rs, diamb, distb, n0, occult = pst
    # pob = [Qis,Uis,ne,phi0,ths,iang,fact]
    Qis, Uis, ne, phi0, ths, iang, fact = pob
    rdi, rdf, H, alpha, ned, dh, ddr, dphi = pdk
    #
    nT = len(phiT)
    I = _np.zeros(nT)
    P = _np.zeros(nT)
    Q = _np.zeros(nT)
    U = _np.zeros(nT)
    # (varr,th,phi_0,n3) = geogen(pst)
    # (varrD,thD,phiD_0) = diskcoords(pdk)
    xD, yD, zD = _phc.sph2cart(varrD, thD, phiD_0)
    #
    # phib blobs in position phibi ##
    ind = (_np.sqrt((xD - 0)**2 + (yD - distb)**2 + (zD - 0)**2) < diamb / 2.)\
        | (_np.sqrt((xD - 0)**2 + (yD + distb)**2 + (zD - 0)**2) < diamb / 2.)
    varrD[ind] = _np.nan
    thD[ind] = _np.nan
    phiD_0[ind] = _np.nan
    varrD = varrD[~_np.isnan(varrD)]
    phiD_0 = phiD_0[~_np.isnan(phiD_0)]
    thD = thD[~_np.isnan(thD)]
    #
    phib = [_np.pi / 2, -_np.pi / 2]
    nT = len(phiT)
    for i in range(nT):
        for phibi in phib:
            phi = phi_0 + phibi + phiT[i] - phi0
            (Ii, Pi, Qi, Ui) = stokes(varr, th, phi, n3, pst, pob)
            I[i] = I[i] + Ii
            P[i] = P[i] + Pi
            Q[i] = Q[i] + Qi
            U[i] = U[i] + Ui
            phiD = phiD_0 + phiT[i] - phi0
            (Ii, Pi, Qi, Ui) = stokesD(varrD, thD, phiD, pst, pob, pdk)
            I[i] = I[i] + Ii
            P[i] = P[i] + Pi
            Q[i] = Q[i] + Qi
            U[i] = U[i] + Ui
    Q = Q0check(Q)
    return I, P, Q, U, n3


def modcycleD(phiT, pst, pob, pdk):
    """ modcycleD """
    # pst = [rs,diamb,distb,n0,occult]
    rs, diamb, distb, n0, occult = pst
    # pob = [Qis,Uis,ne,phi0,ths,iang,fact]
    Qis, Uis, ne, phi0, ths, iang, fact = pob
    rdi, rdf, H, alpha, ned, dh, ddr, dphi = pdk
    #
    nT = len(phiT)
    I = _np.zeros(nT)
    P = _np.zeros(nT)
    Q = _np.zeros(nT)
    U = _np.zeros(nT)
    (varr, th, phi_0, n3) = geogen(pst)
    (varrD, thD, phiD_0) = diskcoords(pdk)
    xD, yD, zD = _phc.sph2cart(varrD, thD, phiD_0)
    #
    # phib blobs in position phibi ##
    ind = (_np.sqrt((xD - 0)**2 + (yD - distb)**2 + (zD - 0)**2) < diamb / 2.)\
        | (_np.sqrt((xD - 0)**2 + (yD + distb)**2 + (zD - 0)**2) < diamb / 2.)
    varrD[ind] = _np.nan
    thD[ind] = _np.nan
    phiD_0[ind] = _np.nan
    varrD = varrD[~_np.isnan(varrD)]
    phiD_0 = phiD_0[~_np.isnan(phiD_0)]
    thD = thD[~_np.isnan(thD)]
    #
    phib = [0.]
    nT = len(phiT)
    for i in range(nT):
        for phibi in phib:
            phiD = phiD_0 + phiT[i] - phi0
            (Ii, Pi, Qi, Ui) = stokesD(varrD, thD, phiD, pst, pob, pdk)
            I[i] = I[i] + Ii
            P[i] = P[i] + Pi
            Q[i] = Q[i] + Qi
            U[i] = U[i] + Ui
    Q = Q0check(Q)
    return I, P, Q, U, n3


def geogen(pst):
    """ ### BLOB dVs ### """
    # pst = [rs,diamb,distb,n0,occult]
    rs, diamb, distb, n0, occult = pst
    #
    dr = 2. * diamb / n0  # 2xblobdiameter
    dx = distb + dr * (_np.arange(-n0 / 2., n0 / 2) + .5)
    dy = dr * (_np.arange(-n0 / 2., n0 / 2) + .5)
    dz = dr * (_np.arange(-n0 / 2., n0 / 2) + .5)
    #
    varr = _np.zeros((n0, n0, n0))
    th = _np.zeros((n0, n0, n0))
    phi = _np.zeros((n0, n0, n0))
    for i in range(n0):
        for j in range(n0):
            for k in range(n0):
                if _np.sqrt((dx[i] - distb)**2 + dy[j]**2 + dz[k]**2) <= \
                    diamb / 2.:
                    varr[i][j][k] = _np.sqrt(dx[i]**2 + dy[j]**2 + dz[k]**2)
                    th[i][j][k] = _np.arccos(dz[k] / varr[i][j][k])                 
                    phi[i][j][k] = _np.arctan(dy[j] / dx[i])
                else:
                    varr[i][j][k] = _np.nan
                    th[i][j][k] = _np.nan
                    phi[i][j][k] = _np.nan
    #
    varr = varr[~_np.isnan(varr)]
    th = th[~_np.isnan(th)]
    phi = phi[~_np.isnan(phi)]
    n3 = len(varr)  # n**3 is final blob V division
    return varr, th, phi, n3


def stokes(varr, th, phi, n3, pst, pob):
    """ ### STOKES CALC given BLOB dVs ### """
    # pst = [rs,diamb,distb,n0,occult]
    rs, diamb, distb, n0, occult = pst
    # pob = [Qis,Uis,ne,phi0,ths,iang,fact]
    Qis, Uis, ne, phi0, ths, iang, fact = pob
    #
    # x = varr*_np.sin(th)*_np.cos(phi)
    # y = varr*(_np.cos(iang)*_np.sin(th)*_np.sin(phi)+_np.sin(iang)*_
        # np.cos(th))
    # z = varr*(-_np.sin(iang)*_np.sin(th)*_np.sin(phi)+_np.cos(iang)*_
        # np.cos(th))
    (x, y, z) = _phc.sph2cart(varr, th, phi)
    (x, y, z) = _phc.cart_rot(x, y, z, 0., iang, 0.)
    #
    rl = (4 * _np.pi / 3 / n3)**(1. / 3) * diamb / 2.
    # P0 = (rs/varr)**2*sigT*ne*rl*3./16./_np.pi*(rl**2./4/_np.pi/rs**2.)
    P0 = 3. / 8 * (ne * _sigT * rl) * 1. / 2 * (rs / varr)**2 * \
        _np.sqrt(1 - (rs / varr)**2) * (rl / varr)**2
    cosX = _np.sin(iang) * _np.sin(th) * _np.cos(phi) + \
        _np.cos(iang) * _np.cos(th)
    # Ii = _np.zeros(len(varr))+1./(len(varr)*2)
    Ii = P0 * (1 + cosX**2)
    Qi = P0 * ( (x / varr)**2 - (y / varr)**2 )
    Ui = P0 * ( 2. * x * y / varr**2 )
    Pi = _np.sqrt(Qi**2 + Ui**2)
    # Occultation ##
    if occult is True: 
        ind = (x**2 + y**2 < rs**2) & (z < 0)
        Qi[ind] = 0
        Ui[ind] = 0
        Pi[ind] = 0
        Ii[ind] = 0.
        ind = (x**2 + y**2 < rs**2) & (z > 0)
        Ii[ind] = Ii[ind] - (rs / varr[ind])**2 * (_np.exp(-
            ne * _sigT * rl)) * (rl**2. / 4 / _np.pi / varr[ind]**2.)
        # ind = _np.where(z<0) and _np.where(x**2+y**2 < rs**2)
        # Qi[ind] = 0
        # Ui[ind] = 0
        # Pi[ind] = 0
        # Ii[ind] = 0
        # ind = _np.where(z>0) and _np.where(x**2+y**2 < rs**2)
        # Ii[ind] = Ii[ind] - (rs/varr[ind])**2*(_np.exp(-ne*sigT*rl))*
            # (rl**2./4/_np.pi/varr[ind]**2.)
    return .5 + Ii.sum(), fact * Pi.sum(), fact * Qi.sum(), fact * Ui.sum()


def modcycle(phiT, pst, pob):
    """ ### ONE MODEL CALC ### """
    # pst = [rs,diamb,distb,n0,occult]
    # rs,diamb,distb,n0,occult = pst
    # pob = [Qis,Uis,ne,phi0,ths,iang,fact]
    Qis, Uis, ne, phi0, ths, iang, fact = pob
    #
    nT = len(phiT)
    I = _np.zeros(nT)
    P = _np.zeros(nT)
    Q = _np.zeros(nT)
    U = _np.zeros(nT)
    # phib blobs in position phibi ##
    phib = [_np.pi / 2, -_np.pi / 2]
    (varr, th, phi_0, n3) = geogen(pst)
    nT = len(phiT)
    for i in range(nT):
        for phibi in phib:
            phi = phi_0 + phibi + phiT[i] - phi0
            (Ii, Pi, Qi, Ui) = stokes(varr, th, phi, n3, pst, pob)
            I[i] = I[i] + Ii
            P[i] = P[i] + Pi
            Q[i] = Q[i] + Qi
            U[i] = U[i] + Ui
    Q = Q0check(Q)
    return I, P, Q, U, n3


def errorcalc(pst, pob, dobs, pdk, phiobs):
    """ error calc """
    rdi, rdf, H, alpha, ned, dh, ddr, dphi = pdk
    Qis, Uis, ne, phi0, ths, iang, fact = pob
    Qobs, Uobs, sigobs = dobs
    (varr, th, phi_0, n3) = geogen(pst)
    (varrD, thD, phiD_0) = diskcoords(pdk)
    (I, P, Q, U, n3) = modcycleF(
        varr, th, phi_0, n3, varrD, thD, phiD_0, phiobs, pst, pob, pdk)
    #
    (Qobc, Uobc) = mod2obs(Q, U, pob)
    (Qchi2, Uchi2) = chi2calc(Qobc, Uobc, dobs)
    chi2min = Qchi2 + Uchi2
    #
    steps = [-1, 1]
    for step in steps:
        Qis0 = Qis
        Uis0 = Uis
        ne0 = ne
        ths0 = ths
        phi00 = phi0
        ned0 = ned
        alpha0 = alpha
        #
        (varrD, thD, phiD_0) = diskcoords(pdk)
        (I, P, Q, U, n3) = modcycleF(
            varr, th, phi_0, n3, varrD, thD, phiD_0, phiobs, pst, pob, pdk)
        chi2 = chi2min / 2
        i = 1 * step
        while chi2 < 2 * chi2min:
            Qis = Qis * (1 + i * .005)
            pob = [Qis, Uis, ne, phi0, ths, iang, fact]
            (Qobc, Uobc) = mod2obs(Q, U, pob)
            (Qchi2, Uchi2) = chi2calc(Qobc, Uobc, dobs)
            chi2 = Qchi2 + Uchi2
            i = i + 1 * step
        print('Finished Qis! %s, %d' % (Qis, i))
        Qis = Qis0
        pob = [Qis, Uis, ne, phi0, ths, iang, fact]
        #
        chi2 = chi2min / 2
        i = 1 * step
        while chi2 < 2 * chi2min:
            Uis = Uis * (1 + i * .005)
            pob = [Qis, Uis, ne, phi0, ths, iang, fact]
            (Qobc, Uobc) = mod2obs(Q, U, pob)
            (Qchi2, Uchi2) = chi2calc(Qobc, Uobc, dobs)
            chi2 = Qchi2 + Uchi2
            i = i + 1 * step
        print('Finished Uis! %s, %d' % (Uis, i))
        Uis = Uis0
        pob = [Qis, Uis, ne, phi0, ths, iang, fact]
        #
        chi2 = chi2min / 2
        i = 1 * step
        while chi2 < 2 * chi2min / 1.17741:
            ne = ne * (1 + i * .05)
            pob = [Qis, Uis, ne, phi0, ths, iang, fact]
            (I, P, Q, U, n3) = modcycleF(
                varr, th, phi_0, n3, varrD, thD, phiD_0, phiobs, pst, pob, pdk)
            (Qobc, Uobc) = mod2obs(Q, U, pob)
            (Qchi2, Uchi2) = chi2calc(Qobc, Uobc, dobs)
            chi2 = Qchi2 + Uchi2
            i = i + 1 * step
        print('Finished ne! %s, %d' % (ne, i))
        ne = ne0
        pob = [Qis, Uis, ne, phi0, ths, iang, fact]
        #
        chi2 = chi2min / 2
        i = 1 * step
        while chi2 < 2 * chi2min / 1.17741:
            ths = ths * (1 + i * .005)
            pob = [Qis, Uis, ne, phi0, ths, iang, fact]
            (I, P, Q, U, n3) = modcycleF(
                varr, th, phi_0, n3, varrD, thD, phiD_0, phiobs, pst, pob, pdk)
            (Qobc, Uobc) = mod2obs(Q, U, pob)
            (Qchi2, Uchi2) = chi2calc(Qobc, Uobc, dobs)
            chi2 = Qchi2 + Uchi2
            i = i + 1 * step
        print('Finished ths! %.1f, %d' % (ths * 180 / _np.pi, i))
        ths = ths0
        pob = [Qis, Uis, ne, phi0, ths, iang, fact]
        #
        chi2 = chi2min / 2
        i = 1 * step
        while chi2 < 2 * chi2min / 1.17741:
            phi0 = phi0 * (1 + i * .005)
            pob = Qis, Uis, ne, phi0, ths, iang, fact
            (I, P, Q, U, n3) = modcycleF(
                varr, th, phi_0, n3, varrD, thD, phiD_0, phiobs, pst, pob, pdk)
            (Qobc, Uobc) = mod2obs(Q, U, pob)
            (Qchi2, Uchi2) = chi2calc(Qobc, Uobc, dobs)
            chi2 = Qchi2 + Uchi2
            i = i + 1 * step
        print('Finished phi0! %.3f, %d' % (phi0 / 2. / _np.pi, i))
        phi0 = phi00
        pob = [Qis, Uis, ne, phi0, ths, iang, fact]
        #
        chi2 = chi2min / 2
        i = 1 * step
        while chi2 < 2 * chi2min / 1.17741:
            ned = ned * (1 + i * .005)
            pdk = [rdi, rdf, H, alpha, ned, dh, ddr, dphi]
            (I, P, Q, U, n3) = modcycleF(
                varr, th, phi_0, n3, varrD, thD, phiD_0, phiobs, pst, pob, pdk)
            (Qobc, Uobc) = mod2obs(Q, U, pob)
            (Qchi2, Uchi2) = chi2calc(Qobc, Uobc, dobs)
            chi2 = Qchi2 + Uchi2
            i = i + 1 * step
        print('Finished ned! %.2e, %d' % (ned, i))
        ned = ned0
        pdk = [rdi, rdf, H, alpha, ned, dh, ddr, dphi]
        #
        chi2 = chi2min / 2
        i = 1 * step
        while chi2 < 2 * chi2min / 1.17741:
            alpha = alpha * (1 + i * .05)
            pdk = [rdi, rdf, H, alpha, ned, dh, ddr, dphi]
            (varrD, thD, phiD_0) = diskcoords(pdk)
            (I, P, Q, U, n3) = modcycleF(
                varr, th, phi_0, n3, varrD, thD, phiD_0, phiobs, pst, pob, pdk)
            (Qobc, Uobc) = mod2obs(Q, U, pob)
            (Qchi2, Uchi2) = chi2calc(Qobc, Uobc, dobs)
            chi2 = Qchi2 + Uchi2
            i = i + 1 * step
        print('Finished alpha! %.3f, %d' % (alpha * 180 / _np.pi, i))
        alpha = alpha0
        pdk = [rdi, rdf, H, alpha, ned, dh, ddr, dphi]
    return 


def QUang(Q, U, filter=True):
    """ ### Q,U angles ### """
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
        ind = _np.where((ang - avg) > 2. / 4 * _np.pi)
        ang[ind] = ang[ind] - _np.pi
        ind = _np.where((ang - avg) < -2. / 4 * _np.pi)
        ang[ind] = ang[ind] + _np.pi
    return ang


def mod2obs(Q, U, pob):
    """ ### MODEL TO OBSERV ### """
    # pst = [rs,diamb,distb,n0,occult]
    # rs,diamb,distb,n0,occult = pst
    # pob = [Qis,Uis,ne,phi0,ths,iang,fact]
    Qis, Uis, ne, phi0, ths, iang, fact = pob
    # dobs = [phiobs,Qobs,Uobs,sigobs]
    # phiobs,Qobs,Uobs,sigobs = dobs

    ang = QUang(Q, U)
    Qobc = _np.sqrt(Q**2 + U**2) * _np.cos(2 * (ang + ths))
    Uobc = _np.sqrt(Q**2 + U**2) * _np.sin(2 * (ang + ths))         
    Qobc = Qobc + Qis
    Uobc = Uobc + Uis
    return Qobc, Uobc


def obs2mod(pob, dobs):
    """ ### MODEL TO OBSERV ### """
    # pst = [rs,diamb,distb,n0,occult]
    # rs,diamb,distb,n0,occult = pst
    # pob = [Qis,Uis,ne,phi0,ths,iang]
    Qis, Uis, ne, phi0, ths, iang, fact = pob
    # dobs = [phiobs,Qobs,Uobs,sigobs]
    Qobs, Uobs, sigobs = dobs

    Qcob = Qobs - Qis
    Ucob = Uobs - Uis
    ang = QUang(Qcob, Ucob)
    Qcob = _np.sqrt(Qcob**2 + Ucob**2) * _np.cos(2 * (ang - ths))
    Ucob = _np.sqrt(Qcob**2 + Ucob**2) * _np.sin(2 * (ang - ths))         
    return Qcob, Ucob


def chi2calc_old(Qobc, Uobc, dobs):
    """ ### CHI2 CALC given Pobs and Pobc ### """
    Qobs, Uobs, sigobs = dobs
    Qchi2 = (Qobc - Qobs)**2 / sigobs**2
    Uchi2 = (Uobc - Uobs)**2 / sigobs**2
    return (Qchi2.sum() + Uchi2.sum()) / (2 * len(sigobs) - 7. - 1.) / 2,\
        (Qchi2.sum() + Uchi2.sum()) / (2 * len(sigobs) - 7. - 1.) / 2


def chi2calc(Qobc, Uobc, dobs):
    """ ### CHI2 CALC given Pobs and Pobc ### """
    Qobs, Uobs, sigobs = dobs
    Qchi2 = (Qobc - Qobs)**2 / sigobs**2 + (Uobc - Uobs)**2 / sigobs**2 
    return Qchi2.sum() / (2 * len(sigobs) - 7. - 1.) / 2,\
        Qchi2.sum() / (2 * len(sigobs) - 7. - 1.) / 2


def Q0check(Q):
    """ ### Q0 check ### """
    ind = _np.where(Q == 0)
    if len(ind) > 0:
        Q[ind] = _np.zeros(len(ind[0])) + 1e-34
    return Q


def extend(phiT, Q, U):
    """ extend """
    ind = (phiT > 0.8 * 2 * _np.pi)
    phiT = _np.hstack((phiT[ind] - 1. * 2 * _np.pi, phiT))
    Q = _np.hstack((Q[ind], Q))
    U = _np.hstack((U[ind], U))
    ind2 = (phiT < 0.2 * 2 * _np.pi) & (phiT > 0.)
    phiT = _np.hstack((phiT, phiT[ind2] + 1. * 2 * _np.pi))
    Q = _np.hstack((Q, Q[ind2]))
    U = _np.hstack((U, U[ind2]))
    #
    ang = QUang(Q, U) * 180. / _np.pi - 90.
    #
    P = _np.sqrt(Q**2 + U**2)
    #
    extmtx = _np.vstack(( (phiT) / 2. / _np.pi, -Q, U, P, ang )).T
    _np.savetxt('mod.txt', extmtx)
    return


def poly_curve_output(best_params, errors):
    '''
    Terminal fit output
    '''
    print("=" * 60)
    print('''        MCMC fitting:       ''')
    print('*' * 15 + ' Fitted parameters ' + '*' * 15)
    for value, sig in zip(best_params, errors):
            # print("{:e} +- {:e}" .format(value, sig))
        print(value, sig)
    return    


def phsort(ph, P):
    """ sort by phases """
    idx = _np.argsort(ph)
    return ph[idx], P[idx]


def extdata(ph, P):
    """ extend data for phases """
    ph, P = phsort(ph, P)
    idx = _np.where(ph > .8)
    phadd0 = ph[idx] - 1.
    Padd0 = P[idx]
    idx = _np.where(ph < .2)
    phadd = ph[idx] + 1.
    Padd = P[idx]
    return _np.hstack((phadd0, ph, phadd)), _np.hstack((Padd0, P, Padd))


class BlobDiskMod(object):

    """ Class BlobDiskMod doc

    Definicoes:

        -fases sempre entre 0-1
        -phiin = phiobs-phi0 = phimod
        -phiobs = phimod+phi0
        -Q1 = Qobs-Qis
        -Qobs = Qmod+Qis
        -Q2 = Pobs**2*_np.cos(2*(ang-ths))
        -Q2obs = Pmod**2*_np.cos(2*(angmod+ths))

        *Pobs = Raw Data
        *P1 = Raw Data - P_IS
        *P2 = P1 * rot(ths)
        *Pmod = intrinsic mod (phi = linspace)
        *Pmodobs = observed mod (phi = linspace)
        *Psetobs = observed mod (phi = phiobs)
    """
    #

    def __init__(self, tgt='sori', rs=4.28, diamb=2 / 3., distb=2.4, n0=5,
        ne=1e12, iang=75., dlt0=-0.17, ths=150.3, Qis=-0.350, Uis=0.025,
        rdi=0., rdf=0., Hd=0.01, alpha=28, dh=1, ddr=3, dphi=180, ned=2.7e12, 
        nbins=0, outname=None, hasdata=True, path=None):
        """ Class initialiser """
        self.tgt = tgt
        if outname is None:
            self.outname = tgt
        else:
            self.outname = outname
        self.path = path
        # star-fixed parameters
        self.rs = rs * _phc.Rsun.cgs
        self.diamb = diamb * self.rs
        self.distb = distb * self.rs
        self.n0 = n0
        # star-variable parameters
        self.ne = ne
        self.iang = iang
        self.irad = self.iang * _np.pi / 180
        self.dlt0 = dlt0
        self.phi0 = self.dlt0 / 2 / _np.pi
        self.ths = ths
        self.trad = self.ths * _np.pi / 180
        self.Qis = Qis
        self.Uis = Uis
        # disk
        if rdi == 0:
            self.rdi = self.distb - self.diamb / 2.
        else:
            self.rdi = rdi
        if rdf == 0:
            self.rdf = self.distb + self.diamb / 2.
        else:
            self.rdf = rdf
        self.Hd = Hd * self.rs
        self.alpha = alpha
        self.alphad = -self.alpha * _np.pi / 180.
        self.dh = dh
        self.ddr = ddr
        self.dphi = dphi
        self.ned = ned
        self.nbins = nbins
        if hasdata is True:
            self.setobs(path=path)
            self.setRmPis()
            self.setRot()
        else:
            self.phiobs = _np.linspace(0, 2 * _np.pi, 80)[:-1]
        self.setmod()
        if hasdata is True:
            self.setbin(nbins)

    def setobs(self, path=None, r=0):
        """ It reads the `tgt.log` file inside path and sets the observation
        variables ph0 and Period (from light-curve database) and phiobs, Pobs,
        Qobs, Uobs, sigP, sigth and angobs from polarimetry. 

        `r` applies the filtraobs with this ratio. """

        if path is None:
            path = _os.getcwd()
        lmags = _np.loadtxt(
            '{0}/refs/pol_mags.txt'.format(_hdt.hdtpath()), dtype=str)
        if not _os.path.exists('{0}/{1}.log'.format(path, self.tgt)):
            self.phiobs = _np.empty(0)
            self.Qobs = _np.empty(0)
            self.Uobs = _np.empty(0)
            self.sigP = _np.empty(0)
            self.sigth = _np.empty(0)
            _warn.warn('Invalid {0} log file. Nothing done.'.format(self.tgt))
        #
        data = _np.loadtxt('{0}/{1}.log'.format(path, self.tgt), dtype=str)
        data = _np.core.records.fromarrays(data.transpose(), names='MJD,night,'
            'filt,calc,ang.ref,dth,P,Q,U,th,sigP,sigQU,sigth', formats='f8,a7,'
            'a1,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8')
        idx = _np.where(data['filt'] == 'v')
        data = data[idx]
        if r > 0:
            data = _polt.filtraobs(data, r=r)
        #
        idx = _np.where(lmags[:, 0] == self.tgt)
        Period, ph0 = lmags[idx][0][1:]
        self.Period = Period
        ph0 = float(ph0) - _jdcal.MJD_0
        self.ph0 = ph0
        #
        phase = data['MJD'] - ph0
        phase /= float(Period)
        phase = _np.modf(phase)[0]
        idx = _np.where(phase < 0)
        if len(idx[0]) > 0:
            raise ValueError('BlobDiskMod calc phases error')
        #
        self.phiobs = phase
        self.Pobs = data['P']
        self.Qobs = data['Q']   
        self.Uobs = data['U']
        self.angobs = QUang(self.Qobs, self.Uobs)
        self.sigP = data['sigP']
        self.sigth = data['sigth']
        return

    def getobs(self):
        """ print obs. info """
        return self.phiobs, self.phiin, self.Pobs, self.Qobs, self.Uobs, \
            self.angobs

    def setRmPis(self):
        """ Remove P_IS """
        self.Q1 = self.Qobs - self.Qis
        self.U1 = self.Uobs - self.Uis
        self.P1 = _np.sqrt(self.Q1**2 + self.U1**2)
        self.ang1 = QUang(self.Q1, self.U1)
        return

    def getRmPis(self):
        """ print rmP_IS """
        return self.phiobs, self.P1, self.Q1, self.U1, self.ang1

    def setRot(self):
        """ setRot """
        self.Q2 = self.P1 * _np.cos(2 * (self.ang1 - self.trad))
        self.U2 = self.P1 * _np.sin(2 * (self.ang1 - self.trad))
        self.P2 = _np.sqrt(self.Q2**2 + self.U2**2)
        self.ang2 = QUang(self.Q2, self.U2)
        # print(self.tgt, _np.average(self.U2), self.ths)
        return

    def getRot(self):
        """ setRot """
        return self.phiobs, self.P2, self.Q2, self.U2, self.ang2

    def setmod(self):
        """ set mod """
        self.phimod = _np.linspace(0, 1., 80)[:-1]
        self.phimodobs = self.phimod + self.phi0
        self.pmodrad = self.phimod * _np.pi * 2
        pst = [self.rs, self.diamb, self.distb, self.n0, True]
        pdk = [self.rdi, self.rdf, self.Hd, self.alphad,
            self.ned, self.dh, self.ddr, self.dphi]
        (varr, th, phi_0, n3) = geogen(pst)
        (varrD, thD, phiD_0) = diskcoords(pdk)
        pob = [0., 0., self.ne, 0., 0., self.irad, 1.]
        (I, Pmod, Qmod, Umod, n3) = modcycleF(varr, th, phi_0, n3, varrD, thD, 
            phiD_0, self.pmodrad, pst, pob, pdk)
        self.Qmod = Qmod * 100
        self.Umod = Umod * 100
        self.Pmod = _np.sqrt(self.Qmod**2 + self.Umod**2)
        self.angmod = QUang(self.Qmod, self.Umod)
        pob = [self.Qis, self.Uis, self.ne,
            self.dlt0, self.trad, self.irad, 1.]
        # (I,Pmod,Qmod,Umod,n3) = modcycleF(varr,th,phi_0,n3,varrD,thD,phiD_0,
            # self.pmodrad,pst,pob,pdk)
        # (Qobc,Uobc) = mod2obs(Qmod,Umod,pob)
        (Qobc, Uobc) = mod2obs(self.Qmod, self.Umod, pob)
        self.Qmodobs = Qobc
        self.Umodobs = Uobc
        self.Pmodobs = _np.sqrt(self.Qmodobs**2 + self.Umodobs**2)
        self.angmodobs = QUang(self.Qmodobs, self.Umodobs)
        #
        self.phisetin = self.phiobs - self.phi0
        self.psetrad = self.phisetin * _np.pi * 2
        pst = [self.rs, self.diamb, self.distb, self.n0, True]
        pdk = [self.rdi, self.rdf, self.Hd, self.alphad,
            self.ned, self.dh, self.ddr, self.dphi]
        (varr, th, phi_0, n3) = geogen(pst)
        (varrD, thD, phiD_0) = diskcoords(pdk)
        pob = [0., 0., self.ne, 0., 0., self.irad, 1.]
        (I, Pmod, Qtmp, Utmp, n3) = modcycleF(varr, th, phi_0, n3, varrD, thD, 
            phiD_0, self.psetrad, pst, pob, pdk)
        self.Qsetin = Qtmp * 100
        self.Usetin = Utmp * 100
        self.Psetin = _np.sqrt(self.Qsetin**2 + self.Usetin**2)
        pob = [self.Qis, self.Uis, self.ne,
            self.dlt0, self.trad, self.irad, 1.]
        # (I,Pmod,Qmod,Umod,n3) = modcycleF(varr,th,phi_0,n3,varrD,thD,phiD_0,
            # self.pmodrad,pst,pob,pdk)
        # (Qobc,Uobc) = mod2obs(Qmod,Umod,pob)
        (Qobc, Uobc) = mod2obs(self.Qsetin, self.Usetin, pob)
        self.Qsetobs = Qobc
        self.Usetobs = Uobc
        self.Psetobs = _np.sqrt(self.Qsetobs**2 + self.Usetobs**2)
        self.angsetobs = QUang(self.Qsetobs, self.Usetobs)        
        return

    def getmod(self):
        """ print mod """
        return self.phimodobs, self.phimod, self.Pmod, self.Qmod, self.Umod, \
            self.angmod

    def setbin(self, nbins):
        """ set bin """
        self.nbins = nbins
        if self.nbins > 0:
            self.outname += '_binned'
            phiobsB, self.Pobs, sigPB = _phc.bindata(
                self.phiobs, self.Pobs, self.nbins, yerr=self.sigP)
            self.Qobs = _phc.bindata(
                self.phiobs, self.Qobs, self.nbins, yerr=self.sigP)[1]
            self.Uobs = _phc.bindata(
                self.phiobs, self.Uobs, self.nbins, yerr=self.sigP)[1]
            self.sigP = sigPB
            self.sigth = _phc.bindata(
                self.phiobs, self.sigth, self.nbins, yerr=self.sigth)[2]
            self.phiobs = phiobsB
            self.setRmPis()
            self.setRot()
            self.setmod()
        # else:
            # print('# Warning! Invalid `nbins` value. Nothing done.')
        return

    def calcavgU0(self):
        """ calcavgU0 """
        lths = _np.linspace(0, 180, 1801) * _np.pi / 180.
        # lths = _np.linspace(0,90,9001)*_np.pi/180.
        avgU = _np.inf
        thmin = 0.
        for ths in lths:
            Qrot = self.P1 * _np.cos(2 * (self.ang1 - ths))
            Urot = self.P1 * _np.sin(2 * (self.ang1 - ths))
            avg = _np.abs(_phc.wg_avg_and_std(Urot, self.sigP)[0])
            if avg < avgU:
                avgU = avg
                thmin = ths
        # Qrot = _np.sqrt(Qin**2+Uin**2)*_np.cos(2*(ang-thmin))
        # Urot = _np.sqrt(Qin**2+Uin**2)*_np.sin(2*(ang-thmin))
        # Prot = _np.sqrt(Qrot**2+Urot**2)
        # angrot = QUang(Qrot,Urot)
        self.ths = thmin * 180 / _np.pi
        self.trad = thmin
        return


def plotPolISrot(a):
    """ plot pol IS Rot """
    fig, (ax0, ax1, ax2) = _plt.subplots(
        3, 1, figsize=(5, 9), sharex=True)  # , sharey=True)
    ax0.plot([-0.2, 1.2], [0, 0], color='k', ls=':')

    # a.phiobs += -a.phi0
    # a.phimodobs += -a.phi0
    df = 0.

    ax0.errorbar(*extdata(a.phiobs + df, a.Pobs), yerr=extdata(a.phiobs +
                 df, a.sigP)[1], label='P', color='r', marker='o', ls='')
    ax0.errorbar(*extdata(a.phiobs + df, a.Qobs), yerr=extdata(a.phiobs +
                 df, a.sigP)[1], label='Q', color='k', marker='o', ls='')
    ax0.errorbar(*extdata(a.phiobs + df, a.Uobs), yerr=extdata(a.phiobs +
                 df, a.sigP)[1], label='U', color='b', marker='o', ls='')
    ax0b = ax0.twinx()
    ax0b.errorbar(*extdata(a.phiobs + df, a.angobs * 180 / _np.pi), 
        yerr=extdata(a.phiobs + df, a.sigth)[1], color='g', marker='^', ls='')

    ax1.plot([-0.2, 1.2], [0, 0], color='k', ls=':')
    ax1.errorbar(*extdata(a.phiobs + df, a.P1), yerr=extdata(a.phiobs +
                 df, a.sigP)[1], label='P', color='r', marker='o', ls='')
    ax1.errorbar(*extdata(a.phiobs + df, a.Q1), yerr=extdata(a.phiobs +
                 df, a.sigP)[1], label='Q', color='k', marker='o', ls='')
    ax1.errorbar(*extdata(a.phiobs + df, a.U1), yerr=extdata(a.phiobs +
                 df, a.sigP)[1], label='U', color='b', marker='o', ls='')
    ax1b = ax1.twinx()
    ax1b.errorbar(*extdata(a.phiobs + df, a.ang1 * 180 / _np.pi),
        yerr=extdata(a.phiobs + df, a.sigth)[1], color='g', marker='^', ls='')

    # lths = _np.linspace(0,180,1801)*_np.pi/180.
    # lths = _np.linspace(0,90,901)*_np.pi/180.
    # avgU = _np.inf
    # thmin = 0.
    # for ths in lths:
        # Qrot = _np.sqrt(Qin**2+Uin**2)*_np.cos(2*(ang-ths))
        # Urot = _np.sqrt(Qin**2+Uin**2)*_np.sin(2*(ang-ths))
        # avg = _np.abs(phc.wg_avg_and_std(Urot, sigP)[0])
        # if avg < avgU:
            # avgU = avg
            # thmin = ths
    # Qrot = _np.sqrt(Qin**2+Uin**2)*_np.cos(2*(ang-thmin))
    # Urot = _np.sqrt(Qin**2+Uin**2)*_np.sin(2*(ang-thmin))
    # Prot = _np.sqrt(Qrot**2+Urot**2)
    # angrot = QUang(Qrot,Urot)
    # idx = _np.where(angrot > 3./4*_np.pi)
    # angrot[idx] = angrot[idx]-_np.pi
    # idx = _np.where(angrot < -3./4*_np.pi)
    # angrot[idx] = angrot[idx]+_np.pi

    ax2.plot([-0.2, 1.2], [0, 0], color='k', ls=':')
    ax2.errorbar(*extdata(a.phiobs + df, a.P2), yerr=extdata(a.phiobs +
                 df, a.sigP)[1], label='P', color='r', marker='o', ls='')
    ax2.errorbar(*extdata(a.phiobs + df, a.Q2), yerr=extdata(a.phiobs +
                 df, a.sigP)[1], label='Q', color='k', marker='o', ls='')
    ax2.errorbar(*extdata(a.phiobs + df, a.U2), yerr=extdata(a.phiobs +
                 df, a.sigP)[1], label='U', color='b', marker='o', ls='')
    ax2b = ax2.twinx()
    ax2b.errorbar(*extdata(a.phiobs + df, a.ang2 * 180 / _np.pi),
        yerr=extdata(a.phiobs + df, a.sigth)[1], color='g', marker='^', ls='')
    # ax2b.plot([-0.2,1.2], [0,0], color='k', ls=':')

    ax0.plot(*extdata(a.phimodobs + df, a.Pmodobs),
             label='P', color='r', marker='', ls='-')
    ax0.plot(*extdata(a.phimodobs + df, a.Qmodobs),
             label='Q', color='k', marker='', ls='-')
    ax0.plot(*extdata(a.phimodobs + df, a.Umodobs),
             label='U', color='b', marker='', ls='-')
    ax0.set_xlim([-.2, 1.2])
    ax0b.plot(*extdata(a.phimodobs + df, a.angmodobs * 180 / _np.pi),
              color='g', marker='', ls='-')

    # ax1.plot(*extdata(a.phimodobs+df, a.Pmodobs), label='P', color='r', 
        # marker='', ls='-')
    # ax1.plot(*extdata(a.phimodobs+df, a.Qmodobs), label='Q', color='k', 
        # marker='', ls='-')
    # ax1.plot(*extdata(a.phimodobs+df, a.Umodobs), label='U', color='b', 
        # marker='', ls='-')
    # ax1.set_xlim([-.2,1.2])
    # ax1b.plot(*extdata(a.phimodobs+df, a.angmodobs*180/_np.pi), color='g', 
        # marker='', ls='-')

    ax2.plot(*extdata(a.phimodobs + df, a.Pmod),
             label='P', color='r', marker='', ls='-')
    ax2.plot(*extdata(a.phimodobs + df, a.Qmod),
             label='Q', color='k', marker='', ls='-')
    ax2.plot(*extdata(a.phimodobs + df, a.Umod),
             label='U', color='b', marker='', ls='-')
    ax2.set_xlim([-.2, 1.2])
    ax2b.plot(*extdata(a.phimodobs + df, a.angmod * 180 / _np.pi),
              color='g', marker='', ls='-')

    # ax2.legend()
    _plt.subplots_adjust(top=.96, bottom=.06, hspace=.008)
    if a.nbins < 1:
        ax0.set_title(a.tgt)
        figname = '{0}_PisRot.png'.format(a.outname)
    else:
        ax0.set_title(a.tgt + ' (binned)')
        figname = '{0}_binned_PisRot.png'.format(a.outname)
    #
    fig.savefig(figname, transparent=True)
    # plt.close()
    # print a.tgt, thmin*180/_np.pi, 180-thmin*180/_np.pi, 
        # _np.average(ang)*180/_np.pi, _np.average(angin)*180/_np.pi, 
        # _np.average(angrot)*180/_np.pi
    return


def QUplots(a):
    """ QU plots """
    fig, (ax0, ax1) = _plt.subplots(2, 1, figsize=(5, 7))  # , sharey=True)

    ax0.errorbar(
        [a.Qis], [a.Uis], yerr=[.033], xerr=[.033], marker='D', color='gray')
    ax0.errorbar(
        a.Qobs, a.Uobs, yerr=a.sigP / 2, xerr=a.sigP / 2, marker='o', ls='')
    ax0.errorbar(a.Qmodobs, a.Umodobs, marker='x')
    ax0.set_xlabel('$Q$ (%)')
    ax0.set_ylabel('$U$ (%)')
    ax0.axis('equal')

    ax1.errorbar(
        a.Q2, a.U2, yerr=a.sigP / 2, xerr=a.sigP / 2, marker='o', ls='')
    ax1.errorbar(a.Qmod, a.Umod, marker='x')
    ax1.set_xlabel('$Q$ (%)')
    ax1.set_ylabel('$U$ (%)')
    ax1.axis('equal')

    if a.nbins < 1:
        ax0.set_title(a.tgt)
        figname = '{0}_QU.png'.format(a.outname)
    else:
        ax0.set_title(a.tgt + ' (binned)')
        figname = '{0}_binned_QU.png'.format(a.outname)
    fig.savefig(figname, transparent=True)
    return


def chi2f(params, p_info, tgt):
    """ chi2f """
    ind = _np.where(p_info[:, 3] == 0)
    p_info[:, 1][ind] = params
    #
    for i in range( len(p_info) ):
        if p_info[i][1] < p_info[i][0] or p_info[i][1] > p_info[i][2]:
            c2nfact = 1
            c2Pisfact = 1
            chi2 = _np.inf
            return -0.5 * (chi2) - c2nfact - c2Pisfact
    #
    iang, ne, phi0, ths, Qis, Uis, alpha, ned = p_info[:, 1]
    thmin = ths
    lths = _np.linspace(0, 180, 1801) * _np.pi / 180.
    avgU = _np.inf
    P1 = _np.sqrt((tgt.Qobs - Qis)**2 + (tgt.Uobs - Uis)**2)
    ang1 = QUang((tgt.Qobs - Qis), (tgt.Uobs - Uis))
    for ths in lths:
        Urot = P1 * _np.sin(2 * (ang1 - ths))
        avg = _np.abs(_phc.wg_avg_and_std(Urot, tgt.sigP)[0])
        if avg < avgU:
            avgU = avg
            thmin = ths
    # if thmin < p_info[3][0] or thmin > p_info[3][2]:
        # c2nfact = 1
        # c2Pisfact = 1
        # chi2 = _np.inf
        # print thmin
        # return -0.5*(chi2)-c2nfact-c2Pisfact
    #
    pst = [tgt.rs, tgt.diamb, tgt.distb, tgt.n0, True]
    pdk = [tgt.rdi, tgt.rdf, tgt.Hd, tgt.alphad,
        tgt.ned, tgt.dh, tgt.ddr, tgt.dphi]
    (varr, th, phi_0, n3) = geogen(pst)
    (varrD, thD, phiD_0) = diskcoords(pdk)
    # pob = [0.,0.,ne,0.,0.,iang,1.]
    # (I,P,Qsetin,Usetin,n3) = modcycleF(varr,th,phi_0,n3,varrD,thD,phiD_0,
        # tgt.phiobs*2*_np.pi-phi0,pst,pob,pdk)
    pob = [Qis, Uis, ne, phi0, thmin, iang, 1.]
    (I, P, Qsetobs, Usetobs, n3) = modcycleF(varr, th, phi_0, n3, varrD, thD, 
        phiD_0, tgt.phiobs * 2 * _np.pi, pst, pob, pdk)
    Qsetobs *= 100.
    Usetobs *= 100.
    (Qobc, Uobc) = mod2obs(Qsetobs, Usetobs, pob)
    dobs = [tgt.Qobs, tgt.Uobs, tgt.sigP]
    (Qchi2, Uchi2) = chi2calc(Qobc, Uobc, dobs)
    chi2 = Qchi2 + Uchi2
    #
    # dobs = [tgt.Qsetin, tgt.Uobs, tgt.sigP]
    # (Qchi2,Uchi2) = chi2calc(Qobc,Uobc,dobs)
    # chi2 = Qchi2+Uchi2i2
    #
    # 3/8.*(0.5*(_np.log10(p_info[7,2]/ned))+(_np.log10(p_info[1,2]/ne)))
    c2nfact = 0.
    # .5*(_np.abs(Qis-Qis0)/_np.abs(Qis0)+_np.abs(Uis-Uis0)/_np.abs(Uis0))
    c2Pisfact = 0.
    # print -0.5*(chi2), c2nfact, c2Pisfact
    return -0.5 * (chi2) - c2nfact - c2Pisfact


def chi2fin(params, p_info, tgt):
    """ chi2f """
    ind = _np.where(p_info[:, 3] == 0)
    p_info[:, 1][ind] = params
    #
    for i in range( len(p_info) ):
        if p_info[i][1] < p_info[i][0] or p_info[i][1] > p_info[i][2]:
            c2nfact = 1
            c2Pisfact = 1
            chi2 = _np.inf
            return -0.5 * (chi2) - c2nfact - c2Pisfact
    #
    iang, ne, phi0, ths, Qis, Uis, alpha, ned = p_info[:, 1]
    thmin = ths
    lths = _np.linspace(0, 180, 1801) * _np.pi / 180.
    avgU = _np.inf
    P1 = _np.sqrt((tgt.Qobs - Qis)**2 + (tgt.Uobs - Uis)**2)
    ang1 = QUang((tgt.Qobs - Qis), (tgt.Uobs - Uis))
    for ths in lths:
        Urot = P1 * _np.sin(2 * (ang1 - ths))
        avg = _np.abs(_phc.wg_avg_and_std(Urot, tgt.sigP)[0])
        if avg < avgU:
            avgU = avg
            thmin = ths
    # if thmin < p_info[3][0] or thmin > p_info[3][2]:
        # c2nfact = 1
        # c2Pisfact = 1
        # chi2 = _np.inf
        # print thmin
        # return -0.5*(chi2)-c2nfact-c2Pisfact
    #
    pst = [tgt.rs, tgt.diamb, tgt.distb, tgt.n0, True]
    pdk = [tgt.rdi, tgt.rdf, tgt.Hd, tgt.alphad,
        tgt.ned, tgt.dh, tgt.ddr, tgt.dphi]
    (varr, th, phi_0, n3) = geogen(pst)
    (varrD, thD, phiD_0) = diskcoords(pdk)
    pob = [0., 0., ne, 0., 0., iang, 1.]
    (I, P, Qsetin, Usetin, n3) = modcycleF(varr, th, phi_0, n3, varrD,
     thD, phiD_0, tgt.phiobs * 2 * _np.pi - phi0, pst, pob, pdk)
    Qsetin *= 100
    Usetin *= 100
    U2 = P1 * _np.sin(2 * (ang1 - thmin))
    Q2 = P1 * _np.cos(2 * (ang1 - thmin))
    dobs = [Q2, U2, tgt.sigP]
    (Qchi2, Uchi2) = chi2calc(Qsetin, Usetin, dobs)
    chi2 = Qchi2 + Uchi2
    #
    # pob = [Qis,Uis,ne,phi0,thmin,iang,1.]
    # (I,P,Qsetobs,Usetobs,n3) = modcycleF(varr,th,phi_0,n3,varrD,thD,phiD_0,
        # tgt.phiobs*2*_np.pi,pst,pob,pdk)
    # Qsetobs *= 100.
    # Usetobs *= 100.
    # (Qobc,Uobc) = mod2obs(Qsetobs,Usetobs,pob)
    # dobs = [tgt.Qobs, tgt.Uobs, tgt.sigP]
    # (Qchi2,Uchi2) = chi2calc(Qobc,Uobc,dobs)
    # chi2 = Qchi2+Uchi2
    #
    # 3/8.*(0.5*(_np.log10(p_info[7,2]/ned))+(_np.log10(p_info[1,2]/ne)))
    c2nfact = 0.
    # .5*(_np.abs(Qis-Qis0)/_np.abs(Qis0)+_np.abs(Uis-Uis0)/_np.abs(Uis0))
    c2Pisfact = 0.
    # print -0.5*(chi2), c2nfact, c2Pisfact
    return -0.5 * (chi2) - c2nfact - c2Pisfact


def plotResiduals(a):
    """ Residuals """
    ms = 5
    fig = _plt.figure()
    ax0 = _plt.subplot2grid((3, 2), (0, 0), rowspan=2)
    ax2 = _plt.subplot2grid((3, 2), (0, 1), rowspan=2)
    ax1 = _plt.subplot2grid((3, 2), (2, 0))
    ax3 = _plt.subplot2grid((3, 2), (2, 1))
#
    ax0.plot([-.2, 1.2], [0, 0], ls=':', color='k')
    ax0.errorbar(*extdata(a.phiobs, a.Pobs), yerr=extdata(a.phiobs, a.sigP)
                 [1], label='P', color='r', marker='o', ls='', markersize=ms)
    ax0.errorbar(*extdata(a.phiobs, a.Qobs), yerr=extdata(a.phiobs, a.sigP)
                 [1], label='Q', color='k', marker='o', ls='', markersize=ms)
    ax0.errorbar(*extdata(a.phiobs, a.Uobs), yerr=extdata(a.phiobs, a.sigP)
                 [1], label='U', color='b', marker='o', ls='', markersize=ms)
    ax0.set_ylabel('Pol. (%)')
    # ax0.set_xlabel('$\phi$')
    ax0.set_xticklabels([])
#
    ax0.plot(*extdata(a.phiobs, a.Psetobs), label='P',
             color='r', marker='x', ls='', markersize=ms)
    ax0.plot(*extdata(a.phiobs, a.Qsetobs), label='Q',
             color='k', marker='x', ls='', markersize=ms)
    ax0.plot(*extdata(a.phiobs, a.Usetobs), label='U',
             color='b', marker='x', ls='', markersize=ms)
#
    ax1.plot([-.2, 1.2], [0, 0], ls=':', color='k')
    ax1.plot(*extdata(a.phiobs, a.Pobs - a.Psetobs), label='P',
             color='r', marker='^', ls='', markersize=ms)
    ax1.plot(*extdata(a.phiobs, a.Qobs - a.Qsetobs), label='Q',
             color='k', marker='^', ls='', markersize=ms)
    ax1.plot(*extdata(a.phiobs, a.Uobs - a.Usetobs), label='U',
             color='b', marker='^', ls='', markersize=ms)
    ax1.set_ylabel('Residual (Obs.-Model)')
    ax1.set_xlabel('$\phi$')
#
    ax2.plot([-.2, 1.2], [0, 0], ls=':', color='k')
    ax2.errorbar(*extdata(a.phiobs, a.P2), yerr=extdata(a.phiobs, a.sigP)
                 [1], label='P', color='r', marker='o', ls='', markersize=ms)
    ax2.errorbar(*extdata(a.phiobs, a.Q2), yerr=extdata(a.phiobs, a.sigP)
                 [1], label='Q', color='k', marker='o', ls='', markersize=ms)
    ax2.errorbar(*extdata(a.phiobs, a.U2), yerr=extdata(a.phiobs, a.sigP)
                 [1], label='U', color='b', marker='o', ls='', markersize=ms)
    ax2.set_ylabel('Pol. (%)')
    # ax0.set_xlabel('$\phi$')
    ax2.set_xticklabels([])
#
    ax2.plot(*extdata(a.phiobs, a.Psetin), label='P',
             color='r', marker='x', ls='', markersize=ms)
    ax2.plot(*extdata(a.phiobs, a.Qsetin), label='Q',
             color='k', marker='x', ls='', markersize=ms)
    ax2.plot(*extdata(a.phiobs, a.Usetin), label='U',
             color='b', marker='x', ls='', markersize=ms)
#
    ax3.plot([-.2, 1.2], [0, 0], ls=':', color='k')
    ax3.plot(*extdata(a.phiobs, a.P2 - a.Psetin), label='P',
             color='r', marker='^', ls='', markersize=ms)
    ax3.plot(*extdata(a.phiobs, a.Q2 - a.Qsetin), label='Q',
             color='k', marker='^', ls='', markersize=ms)
    ax3.plot(*extdata(a.phiobs, a.U2 - a.Usetin), label='U',
             color='b', marker='^', ls='', markersize=ms)
    ax3.set_ylabel('Residual (Obs.-Model)')
    ax3.set_xlabel('$\phi$')
#
    # ax2.set_xlim([-.2,1.2])
    ax0.set_title('Raw data')
    ax2.set_title('Intrinsic polarization')
    _plt.subplots_adjust(
        top=.95, bottom=.096, hspace=.01, left=.09, right=.95, wspace=.3)
    fig.savefig('{0}_press2.png'.format(a.outname), transparent=True)
    # _plt.close()
    return
    return


def plotCaribe(a):
    """ Plot Caribe """
    ms = 2
    # fig, ((ax0, ax2), (ax1, ax3))  = _plt.subplots(2,2, figsize=(5,7), 
        # sharex=True)
    fig = _plt.figure()
    ax0 = _plt.subplot2grid((3, 2), (0, 0), rowspan=2)
    ax2 = _plt.subplot2grid((3, 2), (0, 1), rowspan=2)
    ax1 = _plt.subplot2grid((3, 2), (2, 0))
    ax3 = _plt.subplot2grid((3, 2), (2, 1))
#
    ax0.plot([-.2, 1.2], [0, 0], ls=':', color='k')
    ax0.errorbar(*extdata(a.phiobs, a.Pobs), yerr=extdata(a.phiobs, a.sigP)
                 [1], label='P', color='r', marker='o', ls='', markersize=ms)
    ax0.errorbar(*extdata(a.phiobs, a.Qobs), yerr=extdata(a.phiobs, a.sigP)
                 [1], label='Q', color='k', marker='o', ls='', markersize=ms)
    ax0.errorbar(*extdata(a.phiobs, a.Uobs), yerr=extdata(a.phiobs, a.sigP)
                 [1], label='U', color='b', marker='o', ls='', markersize=ms)
    ax0.set_ylabel('Pol. (%)')
    # ax0.set_xlabel('$\phi$')
    ax0.set_xticklabels([])
#
    ax1.errorbar(*extdata(a.phiobs, a.angobs * 180 / _np.pi), 
        yerr=extdata(a.phiobs, a.sigth)[1], label='PA', color='g', marker='^', 
        ls='', markersize=ms)
    ax1.set_ylabel('PA (deg.)')
    ax1.set_xlabel('$\phi$')
#
    ax2.plot([-.2, 1.2], [0, 0], ls=':', color='k')
    ax2.errorbar(*extdata(a.phiobs, a.P2), yerr=extdata(a.phiobs, a.sigP / 2)
                 [1], label='P', color='r', marker='o', ls='', markersize=ms)
    ax2.errorbar(*extdata(a.phiobs, a.Q2), yerr=extdata(a.phiobs, a.sigP / 2)
                 [1], label='Q', color='k', marker='o', ls='', markersize=ms)
    ax2.errorbar(*extdata(a.phiobs, a.U2), yerr=extdata(a.phiobs, a.sigP / 2)
                 [1], label='U', color='b', marker='o', ls='', markersize=ms)
    ax2.set_ylim([_np.min(a.U2), _np.max(a.P2)])
#
    ax2.plot(*extdata(a.phimodobs, a.Pmod), label='P',
             color='r', marker='', ls='-', markersize=ms)
    ax2.plot(*extdata(a.phimodobs, a.Qmod), label='Q',
             color='k', marker='', ls='-', markersize=ms)
    ax2.plot(*extdata(a.phimodobs, a.Umod), label='U',
             color='b', marker='', ls='-', markersize=ms)
#
    ax2.set_xticklabels([])
    ax2.set_ylabel('Pol. (%)')
    # ax2.set_xlabel('$\phi$')
#
    ax3.plot([-.2, 1.2], [0, 0], ls=':', color='k')
    ax3.plot(*extdata(a.phimodobs, a.angmod * 180 / _np.pi),
             color='g', marker='', ls='-', markersize=ms)
#
    ax3.errorbar(*extdata(a.phiobs, a.ang2 * 180 / _np.pi), 
        yerr=extdata(a.phiobs, a.sigth * 45 / _np.pi)[1], label='PA', 
        color='g', marker='^', ls='', markersize=ms)
    ax3.set_ylabel('PA (deg.)')
    ax3.set_xlabel('$\phi$')
    ax3.set_ylim(
        [_np.min(a.ang2 * 180 / _np.pi), _np.max(a.ang2 * 180 / _np.pi)])
#
    # ax2.set_xlim([-.2,1.2])
    ax0.set_title('Raw data')
    ax2.set_title('Intrinsic polarization')
    _plt.subplots_adjust(
        top=.95, bottom=.096, hspace=.01, left=.09, right=.95, wspace=.3)
    fig.savefig('{0}_press1.png'.format(a.outname), transparent=True)
    # _plt.close()
    return


def plotClean(a, angopt1=False):
    """ Plot Clean """
    ms = 2
    # fig, ((ax0, ax2), (ax1, ax3))  = _plt.subplots(2,2, figsize=(5,7), 
        # sharex=True)
    fig = _plt.figure()
    ax0 = _plt.subplot2grid((3, 2), (0, 0), rowspan=2)
    ax2 = _plt.subplot2grid((3, 2), (0, 1), rowspan=2)
    ax1 = _plt.subplot2grid((3, 2), (2, 0))
    ax3 = _plt.subplot2grid((3, 2), (2, 1))
#
    ax0.plot([-.2, 1.2], [0, 0], ls=':', color='k')
    ax0.errorbar(*extdata(a.phiobs, a.Pobs), yerr=extdata(a.phiobs, a.sigP)
                 [1], label='P', color='r', marker='o', ls='', markersize=ms)
    ax0.errorbar(*extdata(a.phiobs, a.Qobs), yerr=extdata(a.phiobs, a.sigP)
                 [1], label='Q', color='k', marker='o', ls='', markersize=ms)
    ax0.errorbar(*extdata(a.phiobs, a.Uobs), yerr=extdata(a.phiobs, a.sigP)
                 [1], label='U', color='b', marker='o', ls='', markersize=ms)
    ax0.set_ylabel('Pol. (%)')
    # ax0.set_xlabel('$\phi$')
    ax0.set_xticklabels([])
#
    ax1.errorbar(*extdata(a.phiobs, a.angobs * 180 / _np.pi), 
        yerr=extdata(a.phiobs, a.sigth)[1], label='PA', color='g', marker='^', 
        ls='', markersize=ms)
    ax1.set_ylabel('PA (deg.)')
    ax1.set_xlabel('$\phi$')
#
    ax2.plot([-.2, 1.2], [0, 0], ls=':', color='k')
    ax2.errorbar(*extdata(a.phiobs, a.P2), yerr=extdata(a.phiobs, a.sigP / 2)
                 [1], label='P', color='r', marker='o', ls='', markersize=ms)
    ax2.errorbar(*extdata(a.phiobs, a.Q2), yerr=extdata(a.phiobs, a.sigP / 2)
                 [1], label='Q', color='k', marker='o', ls='', markersize=ms)
    ax2.errorbar(*extdata(a.phiobs, a.U2), yerr=extdata(a.phiobs, a.sigP / 2)
                 [1], label='U', color='b', marker='o', ls='', markersize=ms)
    # ax2.set_ylim([_np.min(a.U2),_np.max(a.P2)])
#
    ax2.set_xticklabels([])
    ax2.set_ylabel('Pol. (%)')
    # ax2.set_xlabel('$\phi$')
#
    ang = a.ang2
    if angopt1:
        idx = _np.where(ang > _np.pi / 2)
        ang[idx] -= _np.pi
    ax3.plot([-.2, 1.2], [0, 0], ls=':', color='k')
    ax3.errorbar(*extdata(a.phiobs, ang * 180 / _np.pi), yerr=extdata(a.phiobs,
        a.sigth * 45 / _np.pi)[1], label='PA', color='g', marker='^', ls='', 
        markersize=ms)
    ax3.set_ylabel('PA (deg.)')
    ax3.set_xlabel('$\phi$')
    # ax3.set_ylim([_np.min(a.ang2*180/_np.pi),_np.max(a.ang2*180/_np.pi)])
#
    # ax2.set_xlim([-.2,1.2])
    ax0.set_title('Raw data')
    ax2.set_title('Intrinsic polarization')
    _plt.subplots_adjust(
        top=.95, bottom=.096, hspace=.01, left=.09, right=.95, wspace=.3)
    fig.savefig('{0}_press3.png'.format(a.outname), transparent=True)
    # _plt.close()
    return


def mcmc(tgt, p_info, birun=50, emrun=500, nd_walk=10, a=2, threads=4,
        savesteps=False):
    """ 
    | tgt = target (sst.BlobDiskMod class)
    | p_info = parameters info matrix 
    |
    | birun = burn-in steps
    | emrun = emcee run steps 
    | nd_walk = nwalkers = nd_walk*ndim

    Returns a sst.BlodDiskMod as best-fit and the emcee smampler:
        bestfit, sampler 
    """
    # Build the parameters to be minimized
    ind = _np.where(p_info[:, 3] == 0)
    p0val = p_info[:, 1][ind]
    # Choose an initial set of positions for the walkers.
    ndim = len(p0val)
    # We'll sample with X (nwalkers; 250) walkers.
    nwalkers = ndim * nd_walk
    # Choose an initial set of positions for the walkers.
    if False:
        p0sig = (p_info[:, 2][ind] - p_info[:, 0][ind]) / 6.
        p0 = _np.array([_np.random.normal(p0val, p0sig, ndim)
                       for i in range(nwalkers)])
    else:
        p0 = _np.array([_np.random.uniform(p_info[:, 0][ind], 
            p_info[:, 2][ind], ndim) for i in range(nwalkers)])
    # Initialize the sampler with the chosen specs.
    sampler = _emcee.EnsembleSampler(
        nwalkers, ndim, chi2f, args=[p_info, tgt], a=a, threads=threads)
    # Run a burn-in.
    print("Burning-in ...")
    pos, prob, state = sampler.run_mcmc(p0, birun)
    # Reset the chain to remove the burn-in samples.
    sampler.reset()
    # Starting from the final position in the burn-in chain, sample for 1500
    # steps. (rstate0 is the state of the internal random number generator)
    print("Running MCMC ...")
    localtime = _time.localtime(_time.time())
    if savesteps:
        filename2 = 'r_mcmc%d%02d%02d_%02d%02d.txt' % (localtime[0:5])
        f0 = open(filename2, 'w')
        f0.writelines('# R_MCMC walkers=%d iterations=%d parameters=\n# ' %
            (nwalkers, emrun))
        f0.writelines( '%d\t' * len(ind[0]) % tuple(ind[0]) ) 
        f0.writelines( 'lnprob\n' ) 
        f0.close()
    # 
    for result in sampler.sample(pos, iterations=emrun, rstate0=state):
        pos = result[0]
        prob = result[1]
        state = result[2]
        if savesteps:
            f_handle = file(filename2, 'a')
            _np.savetxt( f_handle, _np.hstack((pos, _np.array([prob]).T)) )
            f_handle.close()        
    # 
    # Print out the mean acceptance fraction. 
    af = sampler.acceptance_fraction
    print("Mean acceptance fraction:", _np.mean(af))
    # Get the index with the highest probability
    maxprob_index = _np.argmax(prob)
    # Get the best parameters and their respective errors
    params_fit = pos[maxprob_index]
    errors_fit = [sampler.flatchain[:, i].std() for i in range(ndim) ]
    # ~ print("# Best fit ...")
    # ~ print('# 1 == Fixed; ',p_info[:,-1])
    # ~ print('# chi2_red = ',_np.max(prob)*-2)
    p_info[:, 1][ind] = params_fit
    err = _np.zeros(len(p_info[:, 1]))
    err[ind] = errors_fit
    fact = _np.array(
        [180. / _np.pi, 1, 2 * _np.pi, 180. / _np.pi, 1, 1, -180. / _np.pi, 1])
    outmtx = _np.hstack(( (_np.array(p_info[:, 1]) * fact).reshape(-1, 1),
    (_np.array(err) * fact).reshape(-1, 1) ))
    for i in range(len(outmtx)):
        print(outmtx[i])

    bestfit = BlobDiskMod(tgt=tgt.tgt, rs=tgt.rs / _phc.Rsun.cgs, 
        diamb=tgt.diamb / tgt.rs, distb=tgt.distb / tgt.rs, n0=tgt.n0, 
        ne=p_info[:, 1][1], iang=p_info[:, 1][0] * 180 / _np.pi, 
        dlt0=p_info[:, 1][2], ths=p_info[:, 1][3] * 180 / _np.pi,
        Qis=p_info[:, 1][4], Uis=p_info[:, 1][5], rdi=0., rdf=0., Hd=0.01, 
        alpha=p_info[:, 1][6] * -1 * 180 / _np.pi, dh=1, ddr=3, dphi=180, 
        ned=p_info[:, 1][7], path=tgt.path)

    bestfit.outname = tgt.outname

    distrib = sampler.chain[:, :, :].reshape((-1, ndim))
    labels = list(_np.array([r'$i$ (rad)', r'$n_e$ (cm$^{-3}$)', 
        r'$\delta$ (rad)', r'$\theta$ (rad)', r'$Q_{IS}$ (%)', r'$U_{IS}$ (%)', 
        r'$\alpha$ (rad)', r'$n_e^d$ (cm$^{-3}$)'])[ind])

    figure = _triangle.corner(distrib, labels=labels, 
        truths=_np.median(distrib, axis=0), quantiles=[0.16, 0.5, 0.84],
        show_titles=True, title_kwargs={"fontsize": 16}, 
        label_kwargs={"fontsize": 16}, title_fmt='.3f')
    print("{0}, median($\chi^2_{{red}}$)= {1:.2f}".format(tgt.tgt, 
        _np.median(prob) * -2))
    # figure.gca().annotate("{0}, median($\chi^2_{{red}}$)= {1:.2f}".
        # format(tgt.tgt, _np.median(prob)*-2), xy=(0.5, 1.0), 
        # xycoords="figure fraction",xytext=(0, -5), 
        # textcoords="offset points",ha="center", va="top", fontsize=20)
    figure.savefig(
        '{0}_triangle{1}.png'.format(bestfit.outname, 0), transparent=False)
    figure.savefig(
        '{0}_triangle{1}.eps'.format(bestfit.outname, 0), transparent=False)
    _plt.close()

    dif = _np.zeros(ndim)
    for i in range(ndim):
        dif[i] = _perct(distrib[:, i], params_fit[i]) / 100. - .5
    figure = _triangle.corner(distrib, labels=labels, truths=params_fit, 
        quantiles=[0.16, 0.5, 0.84], qtdiff=dif, show_titles=True, 
        title_kwargs={"fontsize": 16}, label_kwargs={"fontsize": 16}, 
        title_fmt='.3f')
    print("{0}, $\chi^2_{{red}}$= {1:.2f}".format(tgt.tgt, _np.max(prob) * -2))
    # figure.gca().annotate("{0}, $\chi^2_{{red}}$= {1:.2f}".format(tgt.tgt, 
        # _np.max(prob)*-2), xy=(0.5, 1.0), xycoords="figure fraction",
        # xytext=(0, -5), textcoords="offset points",ha="center", va="top", 
        # fontsize=20)
    figure.savefig(
        '{0}_triangle{1}.png'.format(bestfit.outname, 1), transparent=False)
    figure.savefig(
        '{0}_triangle{1}.eps'.format(bestfit.outname, 1), transparent=False)
    _plt.close()

    f0 = open('{0}_triangle.txt'.format(bestfit.outname), 'w')
    for i in range(len(outmtx)):
        f0.writelines(str(outmtx[i]) + '\n')
    f0.close()

    return bestfit, sampler, params_fit, prob, pos


# MAIN ###
if __name__ == "__main__":
    pass
