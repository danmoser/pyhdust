##!/usr/bin/env python
#-*- coding:utf-8 -*-

"""
Tools for field stars

:author: D. Bednarski
:license: GNU GPL v3.0 (https://github.com/danmoser/pyhdust/blob/master/LICENSE)
"""

import os
import re
import csv
import numpy as np
import matplotlib.pyplot as plt
from itertools import product
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection


#from glob import glob
#import pyhdust.phc as phc
#import pyhdust.jdcal as jdcal
#import datetime as dt
#from pyhdust import hdtpath
#from itertools import product

filters = ['u','b','v','r','i']
fonts = [20, 17, 17, 14]  # Font sizes for titles, axes labels, axes values, key label of graphs

# Dictionary for filter colors for the graphs
colrs = {   'u' : 'Orchid',
            'b' : 'MediumBlue',
            'v' : 'Green',
            'r' : 'Red',
            'i' : 'GoldenRod',
        }

# Dictionary for column numbers of csv table file
idx =  {'be'    : 0,       # Be star name
        'tgt'   : 1,       # target name
        'dRA'   : 2,       # delta RA (in degree)
        'dDEC'  : 3,       # delta DEC (in degree)
        'diang' : 4,       # angular distance (in degree)
        'r'     : 5,       # radial distance (per cent of Be distance)
        'sr'    : 6,       # radial distance error (idem)
        'type'  : 7,       # target star type
        'date'  : 13,      # observation date
        'flt'   : 14,      # filter
        'lbd'   : 15,      # lambda_0 of filter
        'q'     : 16,      # Stokes Q (%)
        'u'     : 17,      # Stokes U (%)
        'p'     : 18,      # pol (%)
        's'     : 19,      # pol error  (%)
        'thet'  : 20,      # pol angle
        'sthet' : 21,      # pol angle error
        'p/s'   : 22,      # ratio p/s above
        'status': 23,      # validation status (ok/ws/wc/on/od/du)
        'notes' : 24,      # notes
        'out'   : 40,      # target out file
        'std1'  : 50,      # standard name #1
        'tgtHD' : 63,      # HD/BD/CD target number
        'plx'   : 64,      # target parallax (pc)
        'splx'  : 65,      # target parallax error (pc)
        'plxbe' : 66,      # be parallax (pc)
        'splxbe': 67,      # be parallax error (pc)
        'RA'    : 71,      # RA (hours, decimal)
        'DEC'   : 72,      # DEC (degree, decimal)
        'std2'  : 77,      # standard name #2
        }
       
# Dictionary for Be names
bes =  {'aeri'  : 'alpha Eri',      '28tau'   : '28 Tau',         'leri'   : 'lambda Eri',
        'oori'  : 'omega Ori',      'acol'    : 'alpha Col',      'kcma'   : 'kappa CMa',
        '27cma' : '27 CMa',         '28cma'   : '28 CMa',         'bcmi'   : 'beta CMi',
        'opup'  : 'omicron Pup',    'ocar'    : 'omega Car',      'mcen'   : 'mu Cen',
        'ecen'  : 'eta Cen',        '48lib'   : '48 Lib',         'dsco'   : 'delta Sco',
        'coph'  : 'chi Oph',        'hd148937': 'HD 148937',      'tsco'   : 'tau Sco',
        'iara'  : 'iota Ara',       '51oph'   : '51 Oph',         'aara'   : 'alpha Ara',
        'lpav'  : 'lambda Pav',     'hr7355'  : 'HR 7355',        'oaqr'   : 'omicron Aqr',
        'paqr'  : 'pi Aqr',         'hr5907'  : 'HR 5907',        '228eri' : '228 Eri',
        'dori'  : 'delta Ori C',    'ztau'    : 'zeta Tau',       'v725tau': 'V725 Tau',
        'pcar'  : 'p Car',          'dcen'    : 'delta Cen',      '39cru'  : '39 Cru',
        'slup'  : 'sigma Lupi',     'klup'    : 'kappa Lupi',     'gara'   : 'gamma Ara',
        'epscap': 'epsilon Cap',    'ecap'    : 'epsilon Cap',    '31peg'  : '31 Peg',
        'epspsa': 'epsilon PsA',    'epsa'    : 'epsilon PsA',    'bpsc'   : 'beta Psc',
        'epstuc': 'epsilon Tuc',    'etuc'    : 'epsilon Tuc',
        }




def readcsv(csvfile, be):
    """
    Read lines of Be star 'be' from 'cvsfile' and return a list with all components
    """
    data = []
    tgt_curr = 'VOID'
    f0 = open(csvfile, 'ro')
    for line in csv.reader(f0, delimiter=';'):
        if line == '':
            continue
        elif line[idx['be']] == be:
            tgt = line[idx['tgt']]
            if tgt != tgt_curr:
                tgt_curr = tgt
                data += [[]]
            # Copy only if is correct data
            if line[idx['status']] not in ('wc', 'ws', 'du'):
                data[-1] += [line]
            
    return data



def getTable(data, x, y, z=None, sx=None, sy=None, sz=None, \
                                vfilter=['du','wc','ws'], nobin=False):
    """
      Receive the list 'data' with many collumns and return the collumns
    concerning to the x, y, z quantities (and their errors sx, sy, sz).

      x, y, z are the label (key) of 'idx' dictionary.
      Returns also 'objarr', list with object names, filters and validation status.

      vfilter is a list whose elements are the status label to be filtered of
    output. The list of status is:

        ws	Without standard
        wc	Without catalog
        on	Standard from other night
        od	Other delta theta
        ok	Ok
        du  Don't use the observation
    """
    objarr  = [[],[],[]]
    xarr = [[],[]]
    yarr = [[],[]]
    if z != None:
        zarr = [[],[]]

    # Verify keys
    for param in (x, y, z, sx, sy, sz):
        if param != None and param not in idx:
            print('Error: key {0} not found'.format(param))
            return 1

    for block in data:
        for line in block:
            # Test the filter
            validate = True
            for vfilt in vfilter:
                if line[idx['status']] == vfilt:
                    validate = False
                    break
            if validate:
                objarr[0] +=  [line[idx['tgt']]]
                objarr[1] += [line[idx['flt']]]
                objarr[2] += [line[idx['status']]]
                xarr[0] += [line[idx[x]]]
                yarr[0] += [line[idx[y]]]
                if sx != None:
                    xarr[1] += [line[idx[sx]]]
                else:
                    xarr[1] += ['']
                if sy != None:
                    yarr[1] += [line[idx[sy]]]
                else:
                    yarr[1] += ['']
                if z != None:
                    zarr[0] += [line[idx[z]]]
                if sz != None:
                    zarr[1] += [line[idx[sz]]]

    try: xarr = [[float(xelem) for xelem in xarr[0]],[float(xelem) for xelem in xarr[1]]]
    except: pass
    try: yarr = [[float(yelem) for yelem in yarr[0]],[float(yelem) for yelem in yarr[1]]]
    except: pass
    try: zarr = [[float(zelem) for zelem in zarr[0]],[float(zelem) for zelem in zarr[1]]]
    except: pass

    if z == None:
        if not nobin:
            bin_data(objarr, xarr, yarr)
        return objarr, xarr, yarr
    else:
        if not nobin:
            bin_data(objarr, xarr, yarr, zarr=zarr)
        return objarr, xarr, yarr, zarr
    



def graf_field(csvfile, be, vfilter=['du','wc','ws'], save=False, nobin=False,\
                                                    extens='pdf'):
    """
    Plot a field graph with polarization directions for Be star 'be'
    through 'csvfile' data table.
    
    """
    plt.close('all')

    # Subtask to calculate the polarization vector coordinates
    def gen_polar(x, y, l, thet):
        thet_rad = np.deg2rad(thet)
        xmin=x-l*np.cos(thet_rad)
        xmax=x+l*np.cos(thet_rad)
        ymin=y-l*np.sin(thet_rad)
        ymax=y+l*np.sin(thet_rad)
        return [xmin, xmax], [ymin, ymax]

    # Subtask to calculate the coordinates of 'error shade'
    def gen_spolar(x, y, l, thet, sthet):
        ang1 = np.deg2rad(thet+sthet)
        ang2 = np.deg2rad(thet-sthet)
        x1=x-l*np.cos(ang1)
        y1=y-l*np.sin(ang1)
        x2=x+l*np.cos(ang1)
        y2=y+l*np.sin(ang1)
        x3=x+l*np.cos(ang2)
        y3=y+l*np.sin(ang2)
        x4=x-l*np.cos(ang2)
        y4=y-l*np.sin(ang2)
        return [[x1,y1], [x2,y2], [x3,y3], [x4,y4]]

    # Read file and generate table
    data = readcsv(csvfile, be)
    if data == []:
        print('No {0} data!'.format(be))
        return
    objarr, raarr, decarr, tharr = getTable(data, 'dRA', 'dDEC', z='thet', sz='sthet', \
                                                        vfilter=vfilter, nobin=nobin)
    if objarr == []:
        print('No {0} valid data!'.format(be))
        return

    fig = plt.figure(1)
    ax = plt.subplot(1, 1, 1)
    ax.set_title('Field Stars - {0}'.format(bes[be]), fontsize=fonts[0], verticalalignment='bottom')
#    ax.set_xlabel(r'RA $-$ RA${}_{Be}$ (degree)', size=fonts[1])
#    ax.set_ylabel(r'DEC $-$ DEC${}_{Be}$ (degree)', size=fonts[1])
    ax.set_xlabel(r'$\Delta$ RA  (degree)', size=fonts[1])
    ax.set_ylabel(r'$\Delta$ DEC  (degree)', size=fonts[1])

    leg = []
    delt = 0.1    # Variable for label names displacement
    lsize = 0.16  # Variable for vectors size
    ymin = 100.  # 100. to pass in first iteration
    # Do the subplots
    for i in range(len(objarr[0])):

        obj = objarr[0][i]
        filt = objarr[1][i]
        x = float(raarr[0][i])
        y = float(decarr[0][i])
        thet = float(tharr[0][i])
        sthet = float(tharr[1][i])
        color = colrs[filt]
        if y < ymin:
            ymin = y
            
        # Plot vectors
        xvert, yvert = gen_polar(x, y, lsize, thet)
        plt.plot(xvert, yvert, color=color)

        # Plot errors
        coords = gen_spolar(x,y,lsize*0.8,thet,sthet)
        polygon = Polygon(coords, True, color=color, alpha=0.3)
        ax.add_patch(polygon)

        # Print object names
        if i==0 or obj not in [ileg[0] for ileg in leg]:
            n=0
            scale = abs(min([float(idecarr) for idecarr in decarr[0]]) - \
                        max([float(idecarr) for idecarr in decarr[0]]))
            for ileg in leg:
                if abs(y - ileg[2]) < delt and abs(x - ileg[1]) < delt:
                    if n==0:
                        x1 = ileg[1]
                        y1 = ileg[2]
                    n += 1
            if n == 0:
                x1=x
                yleg = y-lsize*1.5
                objleg = obj
            else:
                yleg = y1-lsize*1.5-n*(0.04*scale+0.02)
                objleg = '+ '+obj
            ax.text(x1, yleg, '{0}'.format(objleg), horizontalalignment='center', verticalalignment='baseline', fontsize=fonts[3], color='black')
            leg += [[obj,x,y]]
            if yleg < ymin:
                ymin = yleg

    # Print point for Be star
    ax.plot([0,0], [0,0], 'o', color='black', markersize=7)
    # Set limits to show the Be star
    '''
    if ax.get_xlim()[1] < 0:
        ax.set_xlim([ax.get_xlim()[0],0])
    elif ax.get_xlim()[0] > 0:
        ax.set_xlim([0,ax.get_xlim()[1]])
    if ax.get_ylim()[1] < 0:
        ax.set_ylim([ax.get_ylim()[0],0])
    elif ax.get_ylim()[0] > 0:
        ax.set_ylim([0,ax.get_ylim()[1]])
    '''
    ax.set_ylim([ymin-0.1, ax.get_ylim()[1]])
  
    # Fix limits
    ax.autoscale(False)
    ax.plot(ax.get_xlim(), [0,0], 'k--')
    ax.plot([0,0], ax.get_ylim(), 'k--')

    # Setting sizes
    ax.xaxis.label.set_fontsize(fonts[1])
    ax.yaxis.label.set_fontsize(fonts[1])
    for item in (ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fonts[2])

    if save:
        plt.savefig('{0}_field.{1}'.format(be,extens), bbox_inches='tight')
    else:
        plt.show()



# bin same object+filter data
def bin_data(objarr, xarr, yarr, zarr=None):
    """
    Bin data

    objarr == list  [[objects], [filters]]
      xarr == list [[x values], [sx values]]  (idem for yarr and zarr)

      If zarr is void, bin the y values and calculate error by propagation only
    (no stddev). Otherwise, do the same over z values.
      CAUTION: only bins where the contents of objarr and xarr (and yarr,
    case zarr != None) are the same!
    """

    # Check sizes of lists
    for elem in product([len(lista) for lista in (objarr[0], objarr[1], \
                        xarr[0], xarr[1], yarr[0], yarr[1])],repeat=2):
        if elem[0] != elem[1]:
            print('# ERROR: Data binning not processed because the lists have distinct sizes')
            return

    i=0
    # Loop on lines
    while i < len(objarr[0]):
        obj=objarr[0][i]
        fil=objarr[1][i]
        x=xarr[0][i]
        sx=xarr[1][i]
        if zarr == None:
            yarr[1][i] = yarr[1][i]**2
        else:
            y=yarr[0][i]
            sy=yarr[1][i]
            zarr[1][i] = zarr[1][i]**2
        n=1
        j=i+1
        # looking for the same object/filter from i-esim line until the end
        while True:
            # break if another object or end of list
            if j >= len(objarr[0]) or objarr[0][j] != obj:
                break
            # If all j-esim components are equal, except the z values (or y values,
            #   case zarr is void), calculate the mean
            elif objarr[0][j] == obj and objarr[1][j] == fil and \
                             xarr[0][j] == x and xarr[1][j] == sx and \
                       (zarr == None or (yarr[0][j] == y and yarr[1][j] == sy)):
                n+=1
                if zarr == None:
                    yarr[0][i] += yarr[0][j]
                    yarr[1][i] += yarr[1][j]**2
                else:
                    zarr[0][i] += zarr[0][j]
                    zarr[1][i] += zarr[1][j]**2
                    del(zarr[0][j])
                    del(zarr[1][j])
                del(objarr[0][j])
                del(objarr[1][j])
                del(xarr[0][j])
                del(xarr[1][j])
                del(yarr[0][j])
                del(yarr[1][j])
            # skip to next
            else:
                j += 1
        # Conclude the computation of average and propagated error
        if zarr == None:
            yarr[0][i] = yarr[0][i]/n
            yarr[1][i] = np.sqrt(yarr[1][i])/n
        else:
            zarr[0][i] = zarr[0][i]/n
            zarr[1][i] = np.sqrt(zarr[1][i])/n
        i+=1




def graf_theta(csvfile, be, vfilter=['du','wc','ws'], save=False, nobin=False, \
                                                            extens='pdf'):

    plt.close('all')
    data = readcsv(csvfile, be)
    if data == []:
        print('No {0} data!'.format(be))
        return

    objarr, plxarr, tharr = getTable(data, 'plx', 'thet', sx='splx', sy='sthet', \
                                                     vfilter=vfilter, nobin=nobin)
    if objarr == []:
        print('No {0} valid data!'.format(be))
        return

    thmean = meanAngle([[float(elem) for elem in tharr[0]], \
                        [float(elem) for elem in tharr[1]]])
    fig = plt.figure(1)
    ax = plt.subplot(1, 1, 1)

    ax.set_title('Field Stars - {0}'.format(bes[be]), fontsize=fonts[0], verticalalignment='bottom')
    ax.set_xlabel('star/filter', size=fonts[1])
    ax.set_ylabel(r'$\theta$ (degree)', size=fonts[1])

    j = 0
    pts=[[],[],[]]
    # Do the subplots
    while True:
        pts[0] = [i+1 for i in range(len(objarr[0])) if objarr[0][j] == objarr[0][i]]
        pts[1] = [float(tharr[0][i]) for i in range(len(objarr[0])) if objarr[0][j] == objarr[0][i]]
        pts[2] = [float(tharr[1][i]) for i in range(len(objarr[0])) if objarr[0][j] == objarr[0][i]]

#        pts[0].sort(key=lambda x: x[0])
        ax.errorbar(pts[0], pts[1], yerr=pts[2], label = objarr[0][j], linestyle='', marker='o')
        j += len(pts[0])
        if j >= len(objarr[0]):
            break

    ax.set_xlim([0,j+1])
    
    # Setting legend
    ax.legend(loc='upper left', borderaxespad=0., numpoints=1, prop={'size':fonts[3]})
    ax.plot(ax.get_xlim(), [thmean[0], thmean[0]], 'k--')

    # Setting sizes
    ax.xaxis.label.set_fontsize(fonts[1])
    ax.yaxis.label.set_fontsize(fonts[1])
    for item in (ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fonts[2])

    if save:
        plt.savefig('{0}_tht.{1}'.format(be,extens), bbox_inches='tight')
    else:
        plt.show()



def meanAngle(tht):
    """
    Return the mean angle, error propagated and stddev from
        list tht=[[theta], [stheta]]
    """
    mean_err = np.sqrt(sum([elem**2 for elem in tht[1]]))/len(tht[1])

    return np.mean(tht[0]), mean_err, np.std(tht[0])/np.sqrt(len(tht[0]))
  
    
    

def graf_p(csvfile, be, vfilter=['du','wc','ws'], save=False, nobin=False, extens='pdf'):
    """
    Plot  P x wavelength  for star 'be' from 'cvsfile'
    """

    plt.close('all')
    data = readcsv(csvfile, be)
    if data == []:
        print('No {0} data!'.format(be))
        return

    objarr, lbdarr, parr = getTable(data, 'lbd', 'p', sy='s', \
                                    vfilter=vfilter, nobin=nobin)
    if objarr == []:
        print('No {0} valid data!'.format(be))
        return

    fig = plt.figure(1)
    ax = plt.subplot(1, 1, 1)

    ax.set_title('Field Stars - {0}'.format(bes[be]), fontsize=fonts[0], verticalalignment='bottom')
    ax.set_xlabel(r'$\lambda\ (\AA)$', size=fonts[1])
    ax.set_ylabel('P (%)', size=fonts[1])
    ax.set_xlim([2500, 8500])

    j = 0
    pts=[[],[],[]]
    # Do the subplots
    while True:
        pts[0] = [float(lbdarr[0][i]) for i in range(len(objarr[0])) if objarr[0][j] == objarr[0][i]]
        pts[1] = [float(parr[0][i]) for i in range(len(objarr[0])) if objarr[0][j] == objarr[0][i]]
        pts[2] = [float(parr[1][i]) for i in range(len(objarr[0])) if objarr[0][j] == objarr[0][i]]

        # sort lists pts[0],pts[1],pts[2] based on values of pts[0]
        for i in range(len(pts[0])):
            idx = pts[0].index(min(pts[0][i:]))
            if idx != i:
                tmp0 = pts[0][i]
                tmp1 = pts[1][i]
                tmp2 = pts[2][i]
                pts[0][i] = pts[0][idx]
                pts[1][i] = pts[1][idx]
                pts[2][i] = pts[2][idx]
                pts[0][idx] = tmp0
                pts[1][idx] = tmp1
                pts[2][idx] = tmp2

        ax.errorbar(pts[0], pts[1], yerr=pts[2], label='{0:.9s}'.format(objarr[0][j]), \
                                                linestyle='--', marker='o')
        j += len(pts[0])
        if j >= len(objarr[0]):
            break

    # Setting legend
    ax.legend(loc='upper left', borderaxespad=0., numpoints=1, prop={'size':fonts[3]})
#        bbox_to_anchor=(1., .5))

    # Setting sizes
    ax.xaxis.label.set_fontsize(fonts[1])
    ax.yaxis.label.set_fontsize(fonts[1])
    for item in (ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fonts[2])

    if save:
        plt.savefig('{0}_p.{1}'.format(be,extens), bbox_inches='tight')
    else:
        plt.show()




def graf_pradial(csvfile, be, filt, vfilter=['du','wc','ws'], save=False,
                                                    nobin=False, extens='pdf'):
    """
    Plot a field graph with polarization versus distance.
    
    """
    plt.close('all')

    # Read file and generate table
    data = readcsv(csvfile, be)
    if data == []:
        print('No {0} data!'.format(be))
        return
    objarr, plxarr, bearr, parr = getTable(data, 'plx', 'plxbe', sx='splx',
                    sy='splxbe', z='p', sz='s', vfilter=vfilter, nobin=nobin)
    if objarr == []:
        print('No {0} valid data!'.format(be))
        return

    fig = plt.figure(1)
    ax = plt.subplot(1, 1, 1)
    ax.set_title('Field Stars - {0} - filter {1}'.format(bes[be], filt), \
                                        fontsize=fonts[0], verticalalignment='bottom')
    ax.set_xlabel(r'r (pc)', size=fonts[1])
    ax.set_ylabel(r'P (%)', size=fonts[1])

    # Extract data only in filter filt
    obj, x, y = [], [[],[]], [[],[]]
    for i in range(len(objarr[0])):
        if objarr[1][i] == filt:
            try:
                # Skip if plx is less than 0
                if float(plxarr[0][i]) > 0:
                    x[0] += [float(plxarr[0][i])]
                    x[1] += [float(plxarr[1][i])]
                    y[0] += [float(parr[0][i])]
                    y[1] += [float(parr[1][i])]
                    obj += [objarr[0][i]]
            except:
                continue

    '''
    # sort lists based on values of x[0]
    for i in range(len(obj)):
        idx = x[0].index(min(x[0][i:]))
        if idx != i:
            tmp0 = obj[i]
            tmp1 = x[0][i]
            tmp2 = x[1][i]
            tmp3 = y[0][i]
            tmp4 = y[1][i]
            obj[i] = obj[idx]
            x[0][i] = x[0][idx]
            x[1][i] = x[1][idx]
            y[0][i] = y[0][idx]
            y[1][i] = y[1][idx]
            obj[idx] = tmp0
            x[0][idx] = tmp1
            x[1][idx] = tmp2
            y[0][idx] = tmp3
            y[1][idx] = tmp4
    '''
            
    for i in range(len(obj)):
        ax.errorbar(x[0][i], y[0][i], xerr=x[1][i], yerr=y[1][i], label='{0:.9s}'.\
                            format(obj[i]), linestyle='--', marker='o')#, color=colrs[filt])

    '''
    leg = []
    delt = 0.1    # Variable for label names displacement
    lsize = 0.16  # Variable for vectors size
    ymin = 100.  # 100. to pass in first iteration

    for i in range(len(objarr[0])):

        obj = objarr[0][i]
        objfilt = objarr[1][i]
        y = float(parr[0][i])
        sy = float(parr[1][i])
        color = colrs[objfilt]
        try:
            x = float(plxarr[0][i])
            sx = float(plxarr[1][i])
        except:
            continue

        if objfilt != filt:
            continue

        if y < ymin:
            ymin = y

        ax.errorbar(x, y, xerr=sx, yerr=sy, label='{0:.9s}'.format(obj), \
                                                linestyle='--', marker='o', color=color)

        # Print object names
        if i==0 or obj not in [ileg[0] for ileg in leg]:
            n=0
            scale = abs(min([float(iparr) for iparr in parr[0]]) - \
                        max([float(iparr) for iparr in parr[0]]))
            for ileg in leg:
                if abs(y - ileg[2]) < delt and abs(x - ileg[1]) < delt:
                    if n==0:
                        x1 = ileg[1]
                        y1 = ileg[2]
                    n += 1
            if n == 0:
                x1=x
                yleg = y-lsize*1.5
                objleg = obj
            else:
                yleg = y1-lsize*1.5-n*(0.04*scale+0.02)
                objleg = '+ '+obj
            ax.text(x1, yleg, '{0}'.format(objleg), horizontalalignment='center', verticalalignment='baseline', fontsize=fonts[3], color='black')
            leg += [[obj,x,y]]
            if yleg < ymin:
                ymin = yleg

    ax.set_ylim([-0.1*ax.get_ylim()[1], ax.get_ylim()[1]])
    '''

    # Fix limits
    if ax.get_xlim()[1] < bearr[0][0]:
        ax.set_xlim([ax.get_xlim()[0], bearr[0][0]+bearr[1][0]])
    elif ax.get_xlim()[0] > bearr[0][0]:
        ax.set_xlim([bearr[0][0]-bearr[1][0],ax.get_xlim()[1]])

    ax.autoscale(False)
    ax.plot([bearr[0][0],bearr[0][0]], ax.get_ylim(), linestyle ='--', color='gray')
    ax.plot([bearr[0][0]-bearr[1][0],bearr[0][0]-bearr[1][0]], ax.get_ylim(), \
                                    linestyle =':', color='gray')
    ax.plot([bearr[0][0]+bearr[1][0],bearr[0][0]+bearr[1][0]], ax.get_ylim(), \
                                    linestyle =':', color='gray')
    ax.legend(loc='best', borderaxespad=0., numpoints=1, prop={'size':fonts[3]})

    # Setting sizes
    ax.xaxis.label.set_fontsize(fonts[1])
    ax.yaxis.label.set_fontsize(fonts[1])
    for item in (ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fonts[2])

    if save:
        plt.savefig('{0}_field.{1}'.format(be,extens), bbox_inches='tight')
    else:
        plt.show()

        


def genAll(csvfile):

    nobin=False
    save=True
    extens='png'
    if extens == 'eps':
        print('Warning: field graphs will lose the transparency at .eps format')
    for star in bes.keys():
        graf_p(csvfile, star, nobin=nobin, save=save, extens=extens)
        graf_field(csvfile, star, nobin=nobin, save=save, extens=extens)
        graf_theta(csvfile, star, nobin=nobin, save=save, extens=extens)
