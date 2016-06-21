#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
Program to normalize spectra from OPERA Echelle format based on derivative 
variations.
"""
import os
from glob import glob
import numpy as np
import pyhdust.spectools as spt
from scipy import interpolate
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

__version__ = "0.9"
__author__ = "Daniel Moser"
__email__ = "dmfaes@gmail.com"


def plot_orders(ax, fwl, fflx, mode='o'):
    for i in range(len(fwl)):
        ax.plot(fwl[i], fflx[i], mode, ms=1)
    return ax


def do_norm(wl, flx):
    xout, yout = spt.max_func_pts(wl, flx, ws=0.01, avgbox=5)
    tck = interpolate.splrep(xout, yout)
    ynew = interpolate.splev(xout, tck, der=1)
    dmed = np.median( np.abs(ynew) )
    l = len(wl)
    # ConditionA: small derivative 
    conda = ( np.abs(ynew) < dmed )  
    # ConditionB: in the limits
    condb = ( wl[int(l*0.15)] > xout ) | ( xout > wl[int(l*(1-.15))] )
    # xmax = phc.find_nearest(yout, np.max(yout), idx=True)
    # xmax = xout[xmax]
    # xmax = wl[l/2]
    # ConditionC: xout < max(flx) AND ynew > 0
    # condc = ( xout <= xmax ) & ( ynew >= 0 )
    # ConditionD: xout > max(flx) AND ynew < 0
    # condd = ( xout >= xmax ) & ( ynew <= 0 )
    conde = ( ( xout > 658.45 ) | ( xout < 654.1 ) )
    # idx = np.where( (( conda & (condd | condc) ) | condb) & (conde) )  
    idx = np.where( (conda | condb) & (conde) )  
    # & (condd | condc))
    xout = xout[idx]
    yout = yout[idx]
    tck = interpolate.splrep(xout, yout, k=3)
    wl = wl[2:-2]
    flx = flx[2:-2]
    ynew = interpolate.splev(wl, tck, der=0)
    return wl, flx, ynew


path = '.'
opsps = glob('*.spc.gz')
print opsps

avgspecs = []
for sp in opsps:
    os.system('gunzip -c {0} > tmp.txt'.format(sp))
    spec = np.loadtxt('tmp.txt', skiprows=11)

    fig, ax = plt.subplots()
    ax.plot(spec[:, 5], spec[:, 10], label='Leg.')
    # ax.set_title('Title')
    # ax.legend()
    plt.savefig(sp.replace('.spc.gz', '_ori'))
    plt.close(fig)

    # i = 0
    onorm = []
    fwl = []
    oflx = []
    ocont = []
    orders = np.unique(spec[:, 0])
    for o in orders:
        idx = np.where(spec[:, 0] == o)
        wl = spec[:, 5][idx]
        flx = spec[:, 8][idx]
        try:
            wl, flx, ynew = do_norm(wl, flx)
            fwl.append(wl)
            onorm.append(flx/ynew)
            ocont.append(ynew)
            oflx.append(flx)        
        except:
            print('# Error in order {0}'.format(o))

    fig, ax = plt.subplots()
    ax = plot_orders(ax, fwl, oflx, 'o')
    ax = plot_orders(ax, fwl, ocont, '-')
    plt.savefig(sp.replace('.spc.gz', '.orders'))
    plt.close(fig)

    swl, sflx = spt.sum_ec(fwl, oflx)
    swl, scont = spt.sum_ec(fwl, ocont)
    swl, snorm = spt.sum_ec(fwl, onorm)

    avgspecs.append(sp.replace('.spc.gz', '.rv.fits'))
    nans, x = (np.isnan(scont), lambda z: z.nonzero()[0])
    sflx[nans] = 1.
    scont[nans] = 1.
    idx = np.where(scont <= 1e-3)
    scont[idx] = 1.
    sflx[idx] = 1.
    final = sflx/scont
    idx = np.where((final > 0.986) & (final < 1.014))
    pflx = final[idx]
    pwl = swl[idx]
    pflx = np.interp(swl, pwl, savgol_filter(pflx, 3, 1, mode='nearest'))
    final = final/pflx
    spt.writeFits(final, swl*10, 
        savename=sp.replace('.spc.gz', '.rv.fits'), quiet=False, 
        externhd='../../{0}'.format(sp.replace('.spc.gz', '.fits')) )

    fig, ax = plt.subplots()
    ax.plot(swl, final)
    ax.set_ylim([0, 3])
    plt.savefig(figname=sp.replace('.spc.gz', ''))
    plt.close(fig)

if os.path.exists('tmp.txt'):
    os.system('rm tmp.txt')

spt.averagespecs(avgspecs)
