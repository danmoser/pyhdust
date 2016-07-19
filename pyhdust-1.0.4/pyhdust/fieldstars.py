##!/usr/bin/env python
#-*- coding:utf-8 -*-

"""
Tools for field stars

:author: D. Bednarski
:license: GNU GPL v3.0 (https://github.com/danmoser/pyhdust/blob/master/LICENSE)
"""
from __future__ import print_function
import os
import re
import csv
import copy
import numpy as np
#import matplotlib.cm as mplcm
from itertools import product
from glob import glob
from pyhdust import hdtpath
import pyhdust.poltools as polt
import pyhdust.phc as phc
import sys as _sys


def eprint(*args, **kwargs):
    print(*args, file=_sys.stderr, **kwargs)
    return

try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize, hsv_to_rgb
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection
    mpl.rcParams['pdf.fonttype']=42     # Fix the % Symbol in pdf images
except ImportError:
    eprint('# Warning! matplotlib module not installed!!!')

#from glob import glob
#import pyhdust.phc as phc
#import pyhdust.jdcal as jdcal
#import datetime as dt
#from pyhdust import hdtpath
#from itertools import product

filters = ['u','b','v','r','i'] 
fonts = [20, 17, 17, 14, 13]          # Font sizes for titles, axes labels, axes values, key label of graphs, subplot labels


# Dictionary for filter colors for the graphs
colrs = { 'u' : 'Orchid',
          'b' : 'MediumBlue',
          'v' : 'Green',
          'r' : 'Red',
          'i' : 'GoldenRod',
        }

# Dictionary for column numbers of csv table infile
idx1 = {'be?'   : 0,       # line concerning to the Be star?
        'tgt'   : 2,       # target name
        'tgtHD' : 3,       # HD/BD/CD target number
        'be'    : 4,       # Be star name
        'sel?'  : 5,       # Standard selected? ('Y'/'')
        'diang' : 7,       # angular distance (in degree)
        'r'     : 8,       # radial distance (per cent of Be distance)
        'sr'    : 9,       # radial distance error (idem)
        'type'  : 10,      # target star type
        'stype' : 11,      # target star spectral type
        'mplx'  : 13,      # target parallax (")
        'msplx' : 14,      # target parallax error (")
        'plx'   : 15,      # target parallax (pc)
        'splx'  : 16,      # target parallax error (pc)
        'magu'  : 17,      # magnitude filter u
        'magb'  : 18,      # magnitude filter b
        'magv'  : 19,      # magnitude filter v
        'magr'  : 20,      # magnitude filter r
        'magi'  : 21,      # magnitude filter i
        'coor'  : 23,      # coordinates
        'RA'    : 37,      # RA (hours, decimal)
        'DEC'   : 38,      # DEC (degree, decimal)
        }

# Dictionary for column numbers of csv table outfile
idx2 = {'MJD'   : 0,       # MJD
        'date'  : 1,       # date
        'ccd'   : 2,       # CCD name
        'filt'  : 3,       # filter
        'calc'  : 4,       # calcite
        'std'   : 5,       # standard names
        'dth'   : 6,       # delta theta
        'sdth'  : 7,       # delta theta error
        'p'     : 8,       # pol (%)
        'q'     : 9,       # Stokes Q (%)
        'u'     : 10,      # Stokes U (%)
        'thet'  : 11,      # pol angle
        's'     : 12,      # pol error  (%)
        'sthet' : 13,      # pol angle error
        'out'   : 14,      # target out file
        '#star' : 15,      # star number inside the outfile
        'flag'  : 16,      # star flag ('OK'/'W'/'E')
        'tags'  : 17,      # star tags
        'tgt'   : 18,      # target name
        'tgtHD' : 19,      # HD/BD/CD target number
        'be'    : 20,      # Be star name
        'sel?'  : 21,      # Field star selected? ('Y'/'')
        'magu'  : 22,      # magnitude filter u
        'magb'  : 23,      # magnitude filter b
        'magv'  : 24,      # magnitude filter v
        'magr'  : 25,      # magnitude filter r
        'magi'  : 26,      # magnitude filter i
        'type'  : 27,      # target star type
        'stype' : 28,      # target star spectral type
        'diang' : 29,      # angular distance (in degree)
        'plx'   : 30,      # target parallax (pc)
        'splx'  : 31,      # target parallax error (pc)
        'plxbe' : 32,      # Be parallax (pc)
        'splxbe': 33,      # Be parallax error (pc)
        'coor'  : 34,      # target coordinates
        'dRA'   : 35,      # delta RA (target-Be) (degree, decimal)
        'dDEC'  : 36,      # delta DEC (target-Be) (degree, decimal)
        'leff'  : 37,      # lambda_eff (Angstrom)
        'sleff' : 38,      # lambda_eff error (Angstrom)
        }



def readcsv(csvfile, be):
    """
    Read lines of Be star *be* from *cvsfile* and return a
    list with all components:

    [
     [[star 1, filter 1],
      [star 1, filter 2]],
                          [[star 2, filter 1],
                           [star 2, filter 2]],
                                               ...
                                                  ]

    Ignores the observations with 'E' flag.
    
    """
    data = []
    tgt_curr = 'VOID'
    f0 = open(csvfile, 'ro')
    reader = csv.reader(f0, delimiter=';')
    sortedcsv = sorted(reader, key=lambda x: [x[idx2['be']],x[idx2['tgt']],x[idx2['filt']]])

    
    for line in sortedcsv:
        if line == '':
            continue
        elif line[idx2['be']] == be:
            tgt = line[idx2['tgt']]
            if tgt != tgt_curr:
                tgt_curr = tgt
                data += [[]]
            # Copy only if is correct data
            if line[idx2['flag']] != 'E':
                data[-1] += [line]
            
    return data



def gencsv(csvin, path=None, skipdth=False, delta=3.5, epssig=2.0):
    """
    Generate a csvfile with every observations for the standard stars
    listed in pyhdust/refs/pol_hip.txt.

    Compute lambda_eff using color index and an mean airmass X=1.35.
    The error is the (lambda_eff(X=1.7)-lambda_eff(X=1.))/2

    INPUT
        csvin        The csv file with the informations
                     about the standard stars. The columns are
                     indentified through idx1[] dictionary.
        path         The path to the ``red`` directory, inside
                     which are the nights with the observations
        skipdth      print all target values anyway, even when
                     there are no standard star?
        epssig       sigP/P max for unpolarized target (sigP/P
                     up to epssig doesn't need standard star
                     case skipdth==False)
        delta        tolerance for the angle between the two beams
                     of calcite.

    """

    if path == None:
        path = os.getcwd()

    try:
        objs = np.loadtxt('{0}/refs/pol_hip.txt'.format(hdtpath()), dtype=str)
    except:
        print('# ERROR: Can\'t read files pyhdust/refs/pol_hip.txt.')
        raise SystemExit(1)

    fout = open('{0}/dados.csv'.format(path), 'w')
    csvout = csv.writer(fout, delimiter=';')#, quoting=csv.QUOTE_NONE, quotechar='')
    
    # Main loop
    for obj in objs:

        # Generate the table file for each field star
        polt.genTarget(obj, path=path, skipdth=skipdth, delta=delta, epssig=epssig)

        if os.path.exists('{0}/{1}.log'.format(path,obj)):
            # Read informations about the field star and its Be
            with open(csvin, 'ro') as fin:
                beline, tgtline, linetmp = [],[],[]
                for line in csv.reader(fin, delimiter=';'):
                    if line == []:
                        continue
                    # Case the line is for some Be
                    elif line[idx1['be?']] != '':
                        linetmp = line[:]
                    # Case it was matched the field star line
                    elif line[idx1['tgt']] == obj:
                        tgtline += [line[:]]
                        beline += [linetmp[:]]
                # At this point we have two lists: tgtline and beline, which hold the
                # informations about the field star 'obj' and the respective Be star.
                # These variables are lists because more the same field star can be associated
                # with more than a single Be.

            semilines = []
            for i in range(len(tgtline)):
                semilines += [[]]
                # preparar uma semi-linha para ser copiada para todas as linhas com as informações referentes à estrela de campo e suas relações com a Be de referência. Compilar na linha abaixo apenas as informações que são importantes.
                for tag in ('tgt', 'tgtHD', 'be', 'sel?', 'magu', 'magb', 'magv', 'magr', 'magi', 'type', 'stype', 'diang', 'plx', 'splx'):
                    semilines[i] += [tgtline[i][idx1[tag]]]
                semilines[i] += [beline[i][idx1['plx']], beline[i][idx1['splx']]]
                semilines[i] += [tgtline[i][idx1['coor']]]
                semilines[i] += ['{0:.7f}'.format((float(tgtline[i][idx1['RA']])-float(beline[i][idx1['RA']]))*360./24.)]
                semilines[i] += ['{0:.7f}'.format(float(tgtline[i][idx1['DEC']])-float(beline[i][idx1['DEC']]))]
#                print semilines
#                semilines += [tgtline[i]+beline[i]]

#            print '############### AQUI!!!!!!!!!!!!! : ' + str(len(tgtline)) + ' linhas\n'
#            if len(tgtline) > 1:
#                print '##########################################################'
#                print '################### MAIS QUE 1 LINHA #####################'
#                for i in range(len(tgtline)):
#                    print tgtline[i][idx1['tgt']], tgtline[i][idx1['be']]
#                    print beline[i][idx1['be']]
#                print semilines
#                lixo = raw_input('Clique qualquer coisa para continuar...')

            fobj = np.loadtxt('{0}/{1}.log'.format(path,obj), dtype=str, comments='#')
            # Test if is needed to reshape
            if len(fobj) > 0 and type(fobj[0]) != np.ndarray:
                fobj = fobj.reshape(-1,18)
#            print fobj, semilines
            for semiline in semilines:
                for line in fobj.tolist():
                    # Compute the lambda_effective and errors
                    color = ''
                    if line[3] == 'u' and semiline[4] not in ('~','') and semiline[5] not in ('~',''):
                        color = float(semiline[4])-float(semiline[5])
                    elif line[3] in filters[1:] and semiline[5] not in ('~','') and semiline[6] not in ('~',''):
                        color = float(semiline[5])-float(semiline[6])

                    if color != '':
                        leff = polt.lbds(color, line[3], line[2], airmass=1.35)
                        leff1 = polt.lbds(color, line[3], line[2], airmass=1.)
                        leff2 = polt.lbds(color, line[3], line[2], airmass=1.7)
                        sleff = abs(leff2-leff1)/2
#                        print leff, abs(leff1-leff), abs(leff2-leff)
                    else:
                        leff = phc.lbds[line[3]]
                        sleff = 0.
                    csvout.writerow(line+semiline+[leff,sleff])
               
    # Process the found star with the substring 'field'
    polt.genTarget('field', path=path, skipdth=skipdth, delta=delta, epssig=epssig)
    if os.path.exists('{0}/field.log'.format(path)):
        fobj = np.loadtxt('{0}/field.log'.format(path), dtype=str, comments='#')
        # Test if is needed to reshape
        if len(fobj) > 0 and type(fobj[0]) != np.ndarray:
            fobj = fobj.reshape(-1,18)
        for line in fobj.tolist():
            leff = phc.lbds[line[3]]
            csvout.writerow(line[:-1]+[line[-1],'',line[-1].split('_')[0],'Y']+['']*15+[leff])

    fout.close()

    
    return




def getTable(data, x, y, z=None, sx=None, sy=None, sz=None, \
                        vfilter=['no-std'], bin_data=True, onlyY=False, unbias='wk'):
    """
      Receive the list *data* with many collumns and return
      the columns concerning to the *x*, *y*, *z* quantities
      (and their errors *sx*, *sy*, *sz*). Returns also 'objarr',
      a list with object names, filters and validation status.

      IMPORTANT: 1) Propagates error of delta theta factor for
                    'q', 'u', 'p' and 'thet' labels
                 2) Apply unbias correction for P and theta values
                    when specified through 'unbias' variable.
                 3) Case some label is 'q', 'u', 'p' or 'thet',
                    'no-std' is not in vfilter list and
                    bin_data==True, don't bin the data of those
                    which have 'no-std' flag to prevent to compute
                    values WITH NO MEANING. Return various lines
                    for every unbinnable data and a line for the
                    others else binned data.

    INPUT
        x,y,z,sx,sy,sz    Label (key) of *idx2* dictionary.
        bin_data          Bin the data in the observations for
                          the same object+filter+x variable? Note:
                          Case y and/or z are 'p' or 'thet',
                          compute the bins *over the Stokes
                          parameters* and then return the values
                          to 'p' or 'thet'!
        onlyY             True = only selects the list with
                          a "Y" marked
        unbias            Estimator to unbias the data when 'y' or 'z'
                          is equals to 'p' or 'thet'. Three options:
                             a) 'ml': Maximum Likelihood (K=1.41)
                             b) 'wk': Wardle & Kronberg (K=1.0)
                             c) ''  : None (K=0.)

                          where K is the factor to unbias (read the
                          description of routines 'unbiasData' and
                          'meanAngle'.
        vfilter           List whose elements are the labels (tags)
                          to be filtered from the output. It can be
                          'full', 'prob' or 'comp' for these
                          pre-defined lists in polt.vfil dictionary.

        To mount your own list vfilter, select the current tags:

      # General tags for target/standard observation
      - 'bad-mod'         bad modulation
      - 'very-bad-mod'    very bad modulation
      - 'incomp-mods'     some modulations are incompatible
      - 'obs-prob'        some observational problem/error
      - 'iagpol-prob'     polarimeter problem suspected
      - 'other-prob'      another relevant problem
      - 'obs!=pub'        std incompatible with the published

      # Tags for standard status
      - 'no-std'          no standard in the night
      - 'oth-day-std'     standard from another day
      - 'oth-dth'         delta theta estimated from another filter
    """

    def pro_arr (xx,yy,zz=None,sxx=None,syy=None,szz=None):

        objarr  = [[],[],[]]
        xarr = [[],[]]
        yarr = [[],[]]
        zarr = [[],[]]

        for block in data:
            for line in block:

                # Skip the line: 1) case the observation has some tag to be filtered;
                #                2) case the xx, yy or zz has not value in line (i.e.,
                #                   if it is equal to '' or '~');
                #                3) case the observation has not the mark 'Y' that indicates
                #                   it was selected for THAT Be star, since onlyY==True.
                if (onlyY==True and line[idx2['sel?']]!='Y') or \
                                     line[idx2[xx]] in ('','~') or \
                                        line[idx2[yy]] in ('','~') or \
                                        (zz != None and line[idx2[zz]] in ('','~')):
                    validate = False
                else:
                    validate = True
                    for vfilt in vfilter:
                        if vfilt in line[idx2['tags']]:
                            validate = False
#                           print line[idx2['tgt']] + 'FALSE: sel? = ' + line[idx2['sel?']]
                            break

                # IF VALIDATE, proceed
                if validate:
                    objarr[0] +=  [line[idx2['tgt']]]
                    objarr[1] += [line[idx2['filt']]]
                    objarr[2] += [line[idx2['tags']]]
                    xarr[0] += [line[idx2[xx]]]
                    yarr[0] += [line[idx2[yy]]]
                    if sxx != None:
                        xarr[1] += [line[idx2[sxx]]]
                    else:
                        xarr[1] += ['']
                    if syy != None:
                        yarr[1] += [line[idx2[syy]]]
                    else:
                        yarr[1] += ['']
                    if zz != None:
                        zarr[0] += [line[idx2[zz]]]
                    if szz != None:
                        zarr[1] += [line[idx2[szz]]]

        try: xarr[0] = [float(xelem) for xelem in xarr[0]]
        except: pass
        try: xarr[1] = [float(xelem) for xelem in xarr[1]]
        except: pass
        try: yarr[0] = [float(yelem) for yelem in yarr[0]]
        except: pass
        try: yarr[1] = [float(yelem) for yelem in yarr[1]]
        except: pass
        try: zarr[0] = [float(zelem) for zelem in zarr[0]]
        except: pass
        try: zarr[1] = [float(zelem) for zelem in zarr[1]]
        except: pass

        return objarr, xarr, yarr, zarr



    # Verify keys
    for param in (x, y, z, sx, sy, sz):
        if param != None and param not in idx2:
            print('Error: key {0} not found'.format(param))
            return 1

    # Verify if vfilter is a special filter
    if vfilter in polt.vfil.keys():
        vfilter = polt.vfil[vfilter]

    if bin_data:
        if x in ('p', 'thet'):
            print('Error: key \'p\' or \'thet\' can\'t be used as x values when bin_data==True. Try to pass the data that must be binned through the y values or set bin_data=False.')
            return 1
        elif z != None and y in ('p', 'thet'):
            print('Error: keys \'p\' and/or \'thet\' can\'t be used as x and/or y values when bin_data==True and z!=None. Try to pass the data that must be binned through the z values or set bin_data=False.')
            return 1

    if x == 'leff': ignx = True
    else: ignx = False
    if y == 'leff': igny = True
    else: igny = False

    objarr, xarr, yarr, zarr = pro_arr(x, y, z, sx, sy, sz)

    if bin_data and y not in ('p','thet','q','u') and z not in ('p','thet','q','u'):
        if z == None:
            binData(objarr, xarr, yarr, ignx=ignx)
        else:
            binData(objarr, xarr, yarr, zarr=zarr, ignx=ignx, igny=igny)

    # Case y or z is the polarization P, theta, Q or U, bin the data on the *Stokes parameters QU* and
    # compute the new P (or theta, Q, U) values.
    if (y in ('p','thet','q','u') or z in ('p','thet','q','u')):

        objarr, parr, thetarr, dtharr = pro_arr(xx='p', yy='thet', sxx='s', zz='dth', szz='sdth')
        objarr_aux, qarr, uarr, lixo = pro_arr(xx='q', yy='u', sxx='s', syy='s')

        thetarr[1], qarr[1], uarr[1] = polt.propQU(parr[0], thetarr[0], parr[1], dtharr[1])

        if bin_data:
            if z == None:
                xarr_aux = copy.deepcopy(xarr)
                binData(objarr, xarr, qarr, prevent=True, ignx=ignx)
                binData(objarr_aux, xarr_aux, uarr, prevent=True, ignx=ignx, igny=igny)
            else:
                xarr_aux = copy.deepcopy(xarr)
                yarr_aux = copy.deepcopy(yarr)
                binData(objarr, xarr, yarr, zarr=qarr, prevent=True, ignx=ignx, igny=igny)
                binData(objarr_aux, xarr_aux, yarr_aux, zarr=uarr, prevent=True, ignx=ignx, igny=igny)


        if 'thet' in (y, z):

            # Call meanAngle over the qarr and uarr lists, one by one element, just to compute theta correctly
            thetarr[0] = [meanAngle([qi],[ui],[sqi],[sui], estim=unbias)[0] for qi,ui,sqi,sui in zip(qarr[0],uarr[0],qarr[1],uarr[1])]
            thetarr[1] = [meanAngle([qi],[ui],[sqi],[sui], estim=unbias)[1] for qi,ui,sqi,sui in zip(qarr[0],uarr[0],qarr[1],uarr[1])]

            if y=='thet': yarr = thetarr
            else: zarr = thetarr

        if 'p' in (y, z):

            parr[0] = [np.sqrt(qi**2+ui**2) for qi,ui in zip(qarr[0],uarr[0])]
            parr[1] = [np.sqrt((qi*sqi)**2 + (ui*sui)**2)/pi if pi!=0 else 0. for qi,ui,sqi,sui,pi in zip(qarr[0],uarr[0],qarr[1],uarr[1],parr[0])]
            if unbias != '':
                unbiasData(parr[0], parr[1], estim=unbias)

            if y=='p': yarr = parr
            else: zarr = parr

        if y=='q': yarr = qarr
        elif z=='q': zarr = qarr
        if y=='u': yarr = uarr
        elif z=='u': zarr = uarr


    # Return the values
    if z == None:
        return objarr, xarr, yarr
    else:
        return objarr, xarr, yarr, zarr

   


def graf_field(csvfile, be, vfilter=['no-std'], squared=True, save=False, bin_data=True, \
                                                    onlyY=False, extens='pdf', unbias='wk'):
    """
    Plot a field graph with polarization directions for Be star *be*
    through *csvfile* data table.

    squared: use the same scale in both x and y axes?
    
    
    """

    # Subtask to calculate the polarization vector coordinates
    def gen_polar(x, y, l, thet):

        thet_rad = -np.deg2rad(thet)+np.pi/2

        xmin=x-np.cos(thet_rad)*l
        xmax=x+np.cos(thet_rad)*l
        ymin=y-np.sin(thet_rad)*l
        ymax=y+np.sin(thet_rad)*l

        return [xmin, xmax], [ymin, ymax]


    # Subtask to calculate the coordinates of 'error shade'
    def gen_spolar(x, y, l, thet, sthet):

        thet_rad = -np.deg2rad(thet)+np.pi/2
        s_rad = np.deg2rad(sthet)
        l = 0.8*l

        x1=x-l*np.cos(thet_rad-s_rad)/np.cos(s_rad)
        y1=y-l*np.sin(thet_rad-s_rad)/np.cos(s_rad)
        x2=x+l*np.cos(thet_rad-s_rad)/np.cos(s_rad)
        y2=y+l*np.sin(thet_rad-s_rad)/np.cos(s_rad)
        x3=x+l*np.cos(thet_rad+s_rad)/np.cos(s_rad)
        y3=y+l*np.sin(thet_rad+s_rad)/np.cos(s_rad)
        x4=x-l*np.cos(thet_rad+s_rad)/np.cos(s_rad)
        y4=y-l*np.sin(thet_rad+s_rad)/np.cos(s_rad)

#        return [[x1,y1], [x2,y2], [x3,y3], [x4,y4]]
        return [[x1,y1], [x2,y2], [x3,y3], [x4,y4]]


    plt.close('all')
    # Verify if vfilter is a special filter
    if vfilter in polt.vfil.keys():
        vfilter = polt.vfil[vfilter]

    # Read file and generate table
    data = readcsv(csvfile, be)
    if data == []:
        print('No {0} data!'.format(be))
        return
    objarr, raarr, decarr, tharr = getTable(data, 'dRA', 'dDEC', z='thet', sz='sthet', \
                                                vfilter=vfilter, bin_data=bin_data, onlyY=onlyY, unbias=unbias)
#    print tharr
    if objarr == [] or objarr == [[],[],[]]:
        print('No {0} valid data!'.format(be))
        return

    fig = plt.figure(1)
    ax = plt.subplot(1, 1, 1)
    ax.set_title('Field Stars - {0}'.format(phc.bes[be]), fontsize=fonts[0], verticalalignment='bottom')
    ax.set_xlabel(r'$\Delta$ RA  (degree)', size=fonts[1])
    ax.set_ylabel(r'$\Delta$ DEC  (degree)', size=fonts[1])

    leg = []
    vec = [float(veci) for veci in raarr[0]+decarr[0]]
    lsize = 0.075*(max(vec)-min(vec))  # Variable for vectors size
    delt = 1.5*lsize    # Variable for label names displacement

    if squared:
        ax.set_xlim([min(vec)-2*lsize, max(vec)+2*lsize])
        ax.set_ylim([min(vec)-2*lsize, max(vec)+2*lsize])
        ax.autoscale(False)


#    lsize = 0.06*(max(ra)-min(ra))  # Variable for vectors size
#    lysize = 0.1*(max(dec)-min(dec))/(max(dec)+min(dec))  # Variable for vectors size
#    lsize = np.sqrt(lxsize**2 + lysize**2)
#    lsize = 0.1
#    ratio = (max(dec)-min(dec)+2*lsize)/(max(ra)-min(ra)+2*lsize)  # Variable for vectors size
#    ratio=0.65
    ymin = 100.  # 100. to pass in first iteration
    # Do the subplots
#    vec = np.arange(0,180.,180./len(objarr[0]))
#    tharr[0] = vec.tolist()
#    raarr[0] = [0 for i in range(len(raarr[0]))]
#    decarr[0] = [0 for i in range(len(raarr[0]))]
    
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
        coords = gen_spolar(x, y, lsize, thet, sthet)
        polygon = Polygon(coords, True, color=color, alpha=0.18)
        ax.add_patch(polygon)

        # Print object names
        if i==0 or obj not in [ileg[0] for ileg in leg]:
            n=0

            # Compute scale to fix the label positions
            if squared:
                scale = 0.035*abs(ax.get_ylim()[1] - ax.get_ylim()[0])
            else:
                scale = 0.06*abs(min([float(idecarr) for idecarr in decarr[0]]) - \
                            max([float(idecarr) for idecarr in decarr[0]]))

            # Verify if another label shall be closer to the current label
            for ileg in leg:
                if abs(y - ileg[2]) < delt and abs(x - ileg[1]) < delt:
                    if n==0:
                        x1 = ileg[1]
                        y1 = ileg[2]
                    n += 1
            if n == 0:
                x1=x
                yleg = y-lsize*1.5
                objleg = fixName(obj)
            else:
                yleg = y1-lsize*1.5-n*scale
                objleg = '+ ' + fixName(obj)
            ax.text(x1, yleg, '{0}'.format(objleg), horizontalalignment='center', verticalalignment='baseline', fontsize=fonts[3], color='black')
            leg += [[obj,x,y]]
            if yleg < ymin:
                ymin = yleg

    # Print point for Be star
    ax.plot([0,0], [0,0], 'o', color='black', markersize=7)

    # Fix y limits if squared==False
    if not squared:
        ax.set_ylim([ymin-2*lsize, ax.get_ylim()[1]])
  
    # Fix limits
    ax.autoscale(False)
    ax.plot(ax.get_xlim(), [0,0], 'k--')
    ax.plot([0,0], ax.get_ylim(), 'k--')
    plt.gca().invert_xaxis()

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
def binData(objarr, xarr, yarr, zarr=None, ignx=False, igny=False, prevent=False):
    """
    Bin data

    objarr is a list [[objects], [filters], [flags]]
      xarr is a list [[x values], [sx values]]  (idem for yarr and zarr)

    If zarr is void, bin just the y values and calculate error by
    propagation only (no stddev). Otherwise, bin only the z values.

    - prevent:   prevent to bin when flags contain 'no-std'? If True,
                 case there are more than one point able to be binned,
                 bin only the data which don't have 'no-std' flag. The
                 output lists will have a copy of every line from those
                 in input lists with 'no-std' flag, as like as the
                 binned lines among those that just don't have 'no-std'
                 flag.

    CAUTION: only bins when the contents of objarr AND xarr are the
    same (AND yarr also, case zarr != None)! Only if ignx or/and igny
    parameter have value 'True' the routine ignores when xarr or/and
    yarr have not the same values. In this case, the x or/and y value
    to be returned are taken from the first occurrence.

    CONSIDERATIONS:

    - The error of binned data is just the propagated error,
      sqrt(sum(sigma_i^2))/n, and doesn't consider the stddev
    
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
        tags=objarr[2][i]
        x=xarr[0][i]
        sx=xarr[1][i]
        # If it is to prevent, preserve the line without to bin
        # when 'no-std' is in tags
        if prevent and 'no-std' in tags:
            i += 1
            continue
        if zarr == None:
            try:
                yarr[1][i] = yarr[1][i]**2
            except:
                yarr[1][i] = ''
        else:
            y=yarr[0][i]
            sy=yarr[1][i]
            try:
                zarr[1][i] = zarr[1][i]**2
            except:
                zarr[1][i] = ''
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
                       (ignx or (xarr[0][j] == x and xarr[1][j] == sx)) and \
                       (zarr == None or igny or (yarr[0][j] == y and yarr[1][j] == sy)):
                n+=1
                if zarr == None:
                    yarr[0][i] += yarr[0][j]
                    try:
                        yarr[1][i] += yarr[1][j]**2
                    except:
                        pass
                else:
                    zarr[0][i] += zarr[0][j]
                    try:
                        zarr[1][i] += zarr[1][j]**2
                    except:
                        pass
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
            if yarr[1][i] != '':
                yarr[1][i] = np.sqrt(yarr[1][i])/n
        else:
            zarr[0][i] = zarr[0][i]/n
            if zarr[1][i] != '':
                zarr[1][i] = np.sqrt(zarr[1][i])/n
        i+=1



def graf_theta(csvfile, be, vfilter=['no-std'], save=False, bin_data=True, onlyY=False, \
                                                                       extens='pdf', unbias='wk'):

    # Verify if vfilter is a special filter
    if vfilter in polt.vfil.keys():
        vfilter = polt.vfil[vfilter]

    plt.close('all')
    data = readcsv(csvfile, be)
    if data == []:
        print('No {0} data!'.format(be))
        return

    objarr, filtarr, tharr = getTable(data, 'filt', 'thet', sy='sthet', \
                                            vfilter=vfilter, bin_data=bin_data, onlyY=onlyY, unbias=unbias)

    if objarr == [] or objarr == [[],[],[]]:
        print('No {0} valid data!'.format(be))
        return

    # Compute the mean angle
    # A new objarr is needed because the bin_data==False can do the objarr below have a larger length
    objarr_qu, qarr, uarr = getTable(data, 'q', 'u', sx='s', sy='s', \
                                    vfilter=vfilter, bin_data=False, onlyY=onlyY, unbias=unbias)
    thmean = meanAngle(qarr[0], uarr[0], qarr[1], uarr[1], estim=unbias)

        
    fig = plt.figure(1)
    ax = plt.subplot(1, 1, 1)

    ax.set_title('Field Stars - {0}'.format(phc.bes[be]), fontsize=fonts[0], verticalalignment='bottom')
    ax.set_xlabel('star/filter', size=fonts[1])
    ax.set_ylabel(r'$\theta$ (degree)', size=fonts[1])

    j = 0
    pts=[[],[],[]]
    longname = False
    
    # Do the subplots
    while True:
        pts[0] = [i+1 for i in range(len(objarr[0])) if objarr[0][j] == objarr[0][i]]
        pts[1] = [float(tharr[0][i]) for i in range(len(objarr[0])) if objarr[0][j] == objarr[0][i]]
        pts[2] = [float(tharr[1][i]) for i in range(len(objarr[0])) if objarr[0][j] == objarr[0][i]]

#        pts[0].sort(key=lambda x: x[0])
        color = gen_color(csvfile, be, objarr[0][j], onlyY=onlyY)
        nome = fixName(objarr[0][j])
        if len(nome) > 13:
            longname = True
        
        ax.errorbar(pts[0], pts[1], yerr=pts[2], label=nome, linestyle='', marker='o', color=color)
        j += len(pts[0])
        if j >= len(objarr[0]):
            break

    ax.set_xlim([0,j+j/3])
    
    # Setting legend
    leg = ax.legend(loc='best', borderaxespad=0., numpoints=1, prop={'size':fonts[3]})
    leg.get_frame().set_alpha(0.5)

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




def meanAngle(q,u,sq,su,estim='wk'):
    """
    Return the mean angle and the error propagated

    The computation is over QU parameters.


    Unbias theta error using 'estim' estimator:

        if p/sp <= K,   s_theta = psi
        otherwise,      s_theta = propagated error

      where K is given by the estimator related to the
      'estim' variable:

         a) 'ml' : Maximum Likelihood (K=1.41, psi=51.96)
         b) 'wk' : Wardle & Kronberg (K=1.0, psi=51.96)
         c) ''   : None (K=0, psi=51.96)
         d) 'mts': Maier, Tenzer & Santangelo (estimates
                   from Bayesian analysis, psi=61.14)
             
    """

    if estim=='wk':
        k=1.
    elif estim=='ml':
        k=1.41
    elif estim=='':
        k=0.
    elif estim!='mts':
        print('# ERROR: estimation type `{0}` not valid!.'.format(estim))
        raise SystemExit(1)

    qq, uu, sqq, suu = q[:], u[:], sq[:], su[:]
    # Use binData to compute the mean Q and U values
    # espn variables are artifices to bin all data from qq and uu lists using binData
    esp1 = [['' for i in range(len(q))], ['' for i in range(len(q))], ['' for i in range(len(q))]]
    esp2 = copy.deepcopy(esp1)
    esp3 = copy.deepcopy(esp1)
    esp4 = copy.deepcopy(esp1)
    binData(esp1, esp2, [qq,sqq])
    binData(esp3, esp4, [uu,suu])

    # Prevent division by 0
    if qq[0] == 0:
        if uu[0] > 0:
            tht = 45.0
        elif uu[0] < 0:
            tht = 135.0
    else:
#        print qq, uu
#        print sqq, suu
        tht = np.arctan(uu[0]/qq[0])*90/np.pi

    p = np.sqrt(qq[0]**2 + uu[0]**2)
    if p!=0:
        sp = np.sqrt((qq[0]*sqq[0])**2 + (uu[0]*suu[0])**2)/p
        saux = np.sqrt((qq[0]*suu[0])**2 + (uu[0]*sqq[0])**2)/p
        if estim!='mts':
            if sp!=0 and p/saux > k:
                stht = 28.65*saux/p
            else:
                stht = 51.96
        else:
            # We need to use 'saux' below, not 'sp'! It is because if U=0, the error that will
            # influence the theta uncertainty is the sigma_U, because sigma_Q will be along the
            # P direction. Comparing the formulas for sp and saux, the one that do it correctly
            # is saux.
            if sp!=0 and p/saux > 6:
                stht = 28.65*saux/p
            elif saux!=0:
                a=32.50
                b=1.350
                c=0.739
                d=0.801
                e=1.154
                stht = (a*(b+np.tanh(c*(d-p/saux))) - e*p/saux)
            else:
                stht = 61.14
    else:
        tht = 0.
        if estim!='mts':
            stht = 51.96
        else:
            stht = 61.14

    # Fix the angle to the correct in QU diagram
    if qq[0] < 0:
        tht += 90
    if tht < 0:
        tht += 180

    return [tht, stht]



def meanAngle_star(csvfile, be, obj, filts='ubvri', vfilter=['no-std'], onlyY=False, estim='wk'):
    """
        Return the meanAngle for star 'obj' in field of Be
        'be' computed in all filters specified in string
        format in the variable 'filts'.

        For example, if filts='ubv', this routine will
        compute the mean angle among UBV filters.


        Unbias theta error using 'estim' estimator:

            if p/sp <= K,   s_theta = psi
            otherwise,      s_theta = propagated error

          where K is given by the estimator related to the
          'estim' variable:

             a) 'ml' : Maximum Likelihood (K=1.41, psi=51.96)
             b) 'wk' : Wardle & Kronberg (K=1.0, psi=51.96)
             c) ''   : None (K=0, psi=51.96)
             d) 'mts': Maier, Tenzer & Santangelo (estimates
                       from Bayesian analysis, psi=61.14)


        CONSIDERATIONS:

        - Operations over QU parameters

        - Error from standard star is propagated

        - If estim != '' and the final polarization p is
          smaller than its uncertainty multiplied by the
          estimator factor K, the error of theta has the
          value 51.96 or 61.14, depending of the model

    """


    # Verify if vfilter is a special filter
    if vfilter in polt.vfil.keys():
        vfilter = polt.vfil[vfilter]

    # Read tables and unbias data
    data = readcsv(csvfile, be)
    if data == []:
        print('No {0} data!'.format(be))
        return

    # Get tables
    objarr, qarr, uarr = getTable(data, 'q', 'u', sx='s', sy='s', \
                                    vfilter=vfilter, bin_data=False, onlyY=onlyY)
    lixo, parr, arr_aux, = getTable(data, 'p', 'thet', sx='s', sy='sdth', \
                                    vfilter=vfilter, bin_data=False, onlyY=onlyY)

    # Propagate error from standard star
    lixo, qarr[1], uarr[1] = polt.propQU(parr[0], arr_aux[0], parr[1], arr_aux[1])

    qarr_filt, uarr_filt = [[],[]], [[],[]]
    # Create new qarr and uarr with the data of star 'obj' in filters 'filts'
    for i, (obji, flti) in enumerate(zip(objarr[0], objarr[1])):
        if obji == obj and flti in filts:
            qarr_filt[0] += [qarr[0][i]]
            qarr_filt[1] += [qarr[1][i]]
            uarr_filt[0] += [uarr[0][i]]
            uarr_filt[1] += [uarr[1][i]]

    # Compute the thmean
    if qarr_filt != [[],[]]:
        thmean = meanAngle(qarr_filt[0], uarr_filt[0], qarr_filt[1], uarr_filt[1], estim=estim)
    else:
        thmean = [0,0]

    return thmean


    


def gen_color(csvfile, be, obj, onlyY=False):
    """
    Receives a string 'csvfile' pointing to the csvfile,
    a Be star 'be' and the field star named 'obj'. 
    Returns a smart computed np array for the color
    of such field star.

    Note that if you run this task with distinct values
    for 'onlyY', the color returned may be distinct also.

    'be' parameter is necessary because some field stars
    are used for more than one Be star.
    """


    data = readcsv(csvfile, be)
#    print data

    if data == []:
        print('No {0} data!'.format(be))
        return

    stars = []
    for block in data:
#        print block
        line = block[0]
        target = line[idx2['tgt']]
        if target not in stars and ((onlyY and line[idx2['sel?']] == 'Y') or not onlyY):
            stars += [target]

    if len(stars) != 0:
        if len(stars) != 1:
            norm = np.arange(0.04,1.,4.775/(5*(len(stars)-1)))
        else:
            norm = [0.5]
        cm = plt.cm.spectral(norm)
#        norm = Normalize(vmin=0., vmax=len(stars))
#        cm = plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.get_cmap(cmap)).to_rgba(norm)

        if obj in stars:
            idx = stars.index(obj)
            fcm = cm[idx]
        else:
            fcm = np.array([])

    else:
        fcm = np.array([])

    return fcm




def fixName(star):
    """
    Fix the name of star 'star', returning the
    printable name.

    'hd'      ->  'HD '
    '2MASS-'  ->  ''
    'hip'     ->  'HIP '
    'TYC-'    ->  'TYC '

    Return identical 'star' name if the name was not
    identified
    
    """

    def fixx(starr):
        if bool(re.match('^2MASS-J[0-9-+]*$',starr)):
            nome = starr[6:]
        elif bool(re.match('^(hip|HIP)[0-9]*$',starr)):
            nome = 'HIP '+ starr[3:]
        elif bool(re.match('^(tyc|TYC)-[0-9-]*$',starr)):
            nome = 'TYC ' + starr[4:]
        elif bool(re.match('^(hd|HD)[0-9]*$',starr)):
            nome = 'HD ' + starr[2:]
        elif bool(re.match('^h[0-9]*$',starr)):
            nome = 'HD ' + starr[1:]
        elif bool(re.match('^hr[0-9]*$',starr)):
            nome = 'HR ' + starr[2:]
        elif bool(re.match('^bd-[0-9-+]*$',starr)):
            nome = 'BD-' + starr[3:5] + ' ' + starr[6:]
        elif starr in phc.bes:
            nome = phc.bes[starr]
        else:
            nome = starr

        return nome


    if bool(re.match('.*_field',star)):
        if star.split('_')[0] in phc.bes:
            nome = 'Star #' + star[star.index('_field')+6:]
        else:
            nome = fixx(star.split('_')[0]) + ' (Star #' + star[star.index('_field')+6:] + ')'
    elif bool(re.match('.*_.$',star)):
            nome = fixx(star.split('_')[0]) + ' ' + star[star.index('_')+1:].upper()
            
    else:
        nome = fixx(star)


    return nome



def unbiasData(p, s, estim='wk'):
    """
    Unbias P data, appling the operation over input
    lists/arrays:

    p^2 -> p^2 - K^2 s^2

    INPUT

          p: % polarization
          s: pol error
      estim: estimative model:
             a) 'ml': Maximum Likelihood (K=1.41)
             b) 'wk': Wardle & Kronberg (K=1.0)
    
    """

    if estim=='wk':
        k=1.
    elif estim=='ml':
        k=1.41
    elif estim=='mts':
        print('# Warning: changing estimation type from `mts` to `wk` to the % of pol...')
        k=1.
    else:
        print('# ERROR: estimation type `{0}` not valid!.'.format(estim))
        raise SystemExit(1)


    for i in range(len(p)):
        if p[i]<k*s[i]:
            p[i] = 0.
        else:
            p[i] = np.sqrt(p[i]**2-(k*s[i])**2)

    return



def graf_p(csvfile, be, thetfile=None, path=None, vfilter=[], save=False, \
           bin_data=True, onlyY=False, every=False, rotate=False, fit=True, unbias='wk', law='w82', extens='pdf'):
    """
    Plot  P x wavelength for star 'be' and operate over a
    /'be'_is.csv file. The field stars are read from 'csvfile'.

    Fix the polarization bias using Wardle & Kronberg formula.


    'csvfile' : location of dados.csv (field stars data).

    'path' : the path where is located the log file for star 'be'
             (out from polt.genTarget). If None, it is supposed
             inside the current directory (.).

    'bin_data': bin data in graphs?
    
    'onlyY'   : use only the field stars with 'Y' marks (field
                stars originally selected for Be star 'be')?

    'thetfile' : file with the intrinsic angles (oufile from
                 fs.genInt). If 'thetfile' != None, plot the
                 ISP component parallel to the disk of the Be
                 star. It is done by rotating the data in QU
                 diagram at the angle corresponding to the disk
                 inclination.

    'unbias'   :   Estimator to unbias the data when 'y' or 'z'
                   is equals to 'p' or 'thet'. Three options:
                       a) 'ml': Maximum Likelihood (K=1.41)
                       b) 'wk': Wardle & Kronberg (K=1.0)
                       c) ''  : None (K=0.)

                   where K is the factor to unbias (read the
                   description of routines 'unbiasData' and
                   'meanAngle'.

    'every'    : use one intrinsic angle for each one filter to
                 obtain the // component? If every=False makes
                 all data to use a mean value at the -4:-2 collums
                 (22th to 24th) from 'thetfile'

    'fit': fit the ISP using MCMC? Case True, generate the
           graphs and a file ./'be'_is.csv with the best values.
           Write the mean polarization angles inside this file
           also. Case fit==True and there exists this file, this
           routine will ask to the user if he wants to read the
           previous fitted parameters and plot the Serkowski curves
           using them, or if he wants to run MCMC again,
           overwriting this file. If there exists a ./'be'_is.csv
           file and some standard star was not fitted yet, this
           routine will do that and append a line to the csv file.

    'rotate' : DON'T WORKING PROPERLY. Rotate the polarization of
               field stars in QU diagram by the mean angle? A
               option to be explored to replace the unbias procedure.


    CONSIDERATIONS:

    - Error for parallel component is the combination between
    the propagated error and stddev.
    - Error of standard stars is propagated to the field stars
    - The polarization angles (PA) inside 'be'_is.csv are binned
      over Stokes parameters ALWAYS. The errors are computed by the
      simple propagation.
    - The mean PA for each star is computed using QU parameters
      in all filters, even when p/sigma_p<1 (when the individual
      PA error is equals to 51.96), because we are actually
      operating over QU space!
    - p/sigma_p inside 'be'_is.csv is computed by the ratio of
      28.65 and the PA error.

    """


    def copyFit(star, csvwriter, secondcol=None):
        """
        Copy the line of star 'star' inside 'be'_is.csv file
        to the writer csvwriter
        If 'secondcol' != None, do the copy only when the 2nd
        column inside 'be'_is.csv has the value 'secondcol'.

        Return pmax, lmax arrays with the values and
        errors (+ and -).
        """

        pmax, lmax = [], []
        fr = open('{0}_is.csv'.format(be), 'r')
        csvread = csv.reader(fr, delimiter=';')#, quoting=csv.QUOTE_NONE, quotechar='')
        for i, line in enumerate(csvread):
            if line[0] == star and (secondcol==None or line[1]==secondcol):
                csvwriter.writerow(line)
                pmax, lmax = map(lambda v: float(v), line[2:5]), map(lambda v: float(v), line[5:8])
#                break

        fr.close()

        return pmax, lmax


        

    if path == None:
        path = os.getcwd()


    # Verify if vfilter is a special filter
    if vfilter in polt.vfil.keys():
        vfilter = polt.vfil[vfilter]

    plt.close('all')


    # Read tables and unbias data
    data = readcsv(csvfile, be)
    if data == []:
        print('No {0} data!'.format(be))
        return

    objarr, larr, parr = getTable(data, 'leff', 'p', sy='s', \
                                    vfilter=vfilter, bin_data=bin_data, onlyY=onlyY, unbias=unbias)
    if objarr == [] or objarr == [[],[],[]]:
        print('No {0} valid data!'.format(be))
        return

    # Done inside getTable routine now
#    if unbias in ('ml','wk'):
#        unbiasData(parr[0], parr[1], unbias)


    # Try to find a current 'be'_is.csv file, request to the user and initialize
    # the file for the output of the fitted parameters
    if fit:
        if os.path.exists('{0}_is.csv'.format(be)):
            opt = ''
            while opt not in ('1','2'):
                print(('# File with ISP fits {0}_is.csv already exists. Type the option:\n' +\
                       '  1) Use this saved values to plot the fitted Serkowski curve;\n' +\
                       '  2) Run the MCMC to fit all again.').format(be))
                opt = raw_input('Type the option: ')
        
            if opt == '1':
                usePrevious=True
            else:
                usePrevious=False
        else:
            usePrevious=False

        fout = open('{0}_is_tmp.csv'.format(be), 'w')
        csvout = csv.writer(fout, delimiter=';')#, quoting=csv.QUOTE_NONE, quotechar='')
        csvout.writerow(['#obj','#name']+['Pmax','sPmax_+','sPmax_-', 'lmax','slmax_+','slmax_-']+['chi','n']+\
#        csvout.writerow(['#obj','#name']+\
                        ['<th>', 's<th>','<p/sp>']+\
                        ['plot point?', 'use in fit?', 'comments','']+\
                        ['th_u', 'sth_u','th_b', 'sth_b','th_v', 'sth_v','th_r', 'sth_r','th_i', 'sth_i']+\
                        ['p/sp_u', 'p/sp_b', 'p/sp_v', 'p/sp_r', 'p/sp_i'])


    # Get the table for theta values, propagating errors from standard and computing the mean angle by filter.
    # IMPORTANT: Allways bin data below to compute the mean angle in all cases! Allways include 'no-std' to
    # filter a incorrect theta value without standard calibration, even if vfilter doesn't contain this flag.
    if rotate or fit:
        if 'no-std' not in vfilter:
            print('Warning: no tag `no-std` in vfilter variable. Adding this tag in computation of theta values...')
            vfilter_thet = vfilter + ['no-std']
        else:
            vfilter_thet = vfilter
        objarr_tht, lixo, thtarr = getTable(data, 'filt', 'thet', sy='sthet',  \
                                         vfilter=vfilter_thet, bin_data=True, onlyY=onlyY, unbias=unbias)
                            

    # Initialize the graph
    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)

    ax.set_title('Field Stars - {0}'.format(phc.bes[be]), fontsize=fonts[0], verticalalignment='bottom')
    ax.set_xlabel(r'$\lambda\ (\AA)$', size=fonts[1])
    if rotate:
        ax.set_ylabel('Q` (%)', size=fonts[1])
    else:
        ax.set_ylabel('P (%)', size=fonts[1])
    ax.set_xlim([2500, 8500])


    # Initialize Variables
    j = 0
    pts=[[],[],[]]
    tht = [[],[]]
    longname = False


    ##############################
    ####       MAIN LOOP      ####
    # Do the subplots and fits
    while True:
        star = objarr[0][j]
        pts[0] = [float(larr[0][i]) for i in range(len(objarr[0])) if star == objarr[0][i]]
        pts[1] = [float(parr[0][i]) for i in range(len(objarr[0])) if star == objarr[0][i]]
        pts[2] = [float(parr[1][i]) for i in range(len(objarr[0])) if star == objarr[0][i]]

#        print objarr_qu[0][j:]
#        print objarr_qu[1][j:]

        # Get the mean angles
        if rotate or fit:
            thmean, psp = [], []
            for k, ffil in enumerate(filters):
                ver = False
                ni = 0
                # thmean is a plain list with the angles and errors.
                # By the bin performed some lines above, the first matched element is
                # the only one, and concerning to the mean angle computed.
                for i in range(len(objarr_tht[0])):
                    if objarr_tht[0][i] == star and objarr_tht[1][i] == ffil:
                        thmean += ['{0:.2f}'.format(float(thtarr[0][i]))]
                        thmean += ['{0:.2f}'.format(float(thtarr[1][i]))]
                        if float(thtarr[1][i]) not in (0, 51.96, 61.14):
                            psp += [28.65/float(thtarr[1][i])]
                        else:
                            psp += [0.]
                        ver = True
                        break
                if not ver:
                    thmean += [0.,0.]
                    psp += [0.]
                """
                # It is needed to compute the mean p/sp because bin_data can be False
                for i in range(len(objarr[0])):
                    print objarr[0][i], star
                    print objarr[1][i], ffil
                    print parr[0][i], parr[1][i]
                    print float(parr[0][i])/float(parr[1][i])
                    if objarr[0][i] == star and objarr[1][i] == ffil and float(parr[1][i]) != 0.:
                        print '*********ENTROU******************'
                        psp[k] += float(parr[0][i])/float(parr[1][i])
                        ni += 1
                    elif objarr[0][i] == star and objarr[1][i] == ffil:
                        print '*********NÃOOOOOOO ENTROU******************'
                    print ''
                        
                if ni > 1:
                    psp[k] = psp[k]/ni
                """
            means = meanAngle_star(csvfile, be, star, filts='ubvri', vfilter=vfilter_thet, onlyY=onlyY, estim=unbias)
            means += [np.mean(filter(lambda v: v!=0, psp))]

        # Store the new rotated values
        # It is needed TO FIX above because thmean[0] is not the mean theta, but a list with the mean theta and its error for filter u.
        if rotate:
            q[0], u[0], q[1], u[1] = rotQU(q[0], u[0], q[1], u[1], thmean[0], thmean[1])

        # Sort lists pts[0],pts[1],pts[2] based on values of pts[0]
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
                if rotate:
                    tmp3 = q[0][i]
                    tmp4 = u[0][i]
                    tmp5 = q[1][i]
                    tmp6 = u[1][i]
                    q[0][i] = q[0][idx]
                    u[0][i] = u[0][idx]
                    q[1][i] = q[1][idx]
                    u[1][i] = u[1][idx]
                    q[0][idx] = tmp3
                    u[0][idx] = tmp4
                    q[1][idx] = tmp5
                    u[1][idx] = tmp6

        # Get the color and the the common name to plot
        color = gen_color(csvfile, be, star, onlyY=onlyY)
        nome = fixName(star)
        if len(nome) > 13:
            longname = True

        if not rotate:
            ax.errorbar(pts[0], pts[1], yerr=pts[2], label=nome, linestyle='', marker='o', color=color)

            # Fit data HERE or copy from the previous csv file
            if fit:
                if usePrevious:
#                    means = map(lambda v: '{0:.2f}'.format(v), means[:2]) + ['{:.1f}'.format(means[2])]
#                    psp = map(lambda v: '{0:.1f}'.format(v), psp)
                    pmax_fit, lmax_fit = copyFit(star, csvout)
                else:
                    pmax_fit, lmax_fit = [], []
#                print usePrevious
#                print pmax_fit, lmax_fit

                # Run MCMC again if not usePrevious or if there is no 'star' inside the
                # current csv file.
                if (pmax_fit, lmax_fit) == ([],[]):
                    # Convert to microns and fit
                    lb = [lbi/10000 for lbi in pts[0]]
                    pmax_fit, lmax_fit, chi2 = fitSerk(lb, pts[1], pts[2], star=star, extens=extens)

                    # Fix the format to the lists
                    pmax_fit_str = map(lambda v: '{0:.5f}'.format(v), list(pmax_fit))
                    lmax_fit_str = map(lambda v: '{0:.5f}'.format(v), list(lmax_fit))
                    means = map(lambda v: '{0:.2f}'.format(v), means[:2]) + ['{:.1f}'.format(means[2])]
                    psp = map(lambda v: '{0:.1f}'.format(v), psp)
                    chi2 = '{0:.3f}'.format(chi2)

                    csvout.writerow([star, nome] + pmax_fit_str + lmax_fit_str + [chi2, len(pts[1])] +\
                                     means + ['1','1','---',''] + thmean + psp)
#                    csvout.writerow([star, nome] + means + ['1','1','---',''] + thmean + psp)

                ll = np.linspace(3000.,8300.,100)
                pp = np.array([])
                for lli in ll:
                    pp = np.append(pp, polt.serkowski(pmax_fit[0], lmax_fit[0]*10000, lli, law=law, mode=2))

                # Only plot the graph if there are more than one data, because with an only point
                # the curve is not defined! But the emcee was runned to generate the covariance map
                if len(pts[0]) > 1:
                     ax.plot(ll, pp, color=color)

        else:
            ax.errorbar(pts[0], q[0], yerr=q[1], label=nome+'_q', linestyle='-', marker='', color=color)
#            ax.errorbar(pts[0], u[0], yerr=u[1], label=nome+'_u', linestyle='-.', marker='', color=color)
#            ax.errorbar(pts[0], pts[1], yerr=pts[2], label=nome+'_p', linestyle='--', marker='', color=color)
            ax.plot(pts[0], pts[1], label=nome+'_p', linestyle='--', marker='', color=color)


        # Refresh the lines to the next iteration
        j += len(pts[0])
        if j >= len(objarr[0]):
            break

    # Setting legend
#    fig.subplots_adjust(right=0.8)
#    ax.legend(loc='center left', borderaxespad=0., numpoints=1, prop={'size':fonts[4]}, bbox_to_anchor=(1.02,.5))

    # Plot ISP parallel to the disk direction of Be star
    if thetfile != None:

        lbd, pBe, sBe, devBe = rotQUBe(be, thetfile, path=path, every=every, vfilter=vfilter)
        # Using combined error for points to // pol to Be star
        s = [np.sqrt(sBe[i]**2+devBe[i]**2) for i in range(len(sBe))]

        if lbd != []:
            ax.errorbar(lbd, pBe, yerr=s, label=r'$//$ comp.', linestyle='', \
                                                            marker='s', color='black')

            if fit:

                if usePrevious:
                    pmax_fit, lmax_fit = copyFit(be, csvout, secondcol=fixName(be)+' (//)')
                else:
                    pmax_fit, lmax_fit = [], []

                if (pmax_fit,lmax_fit) == ([],[]):
                    # Convert to microns
                    lb = [lbi/10000 for lbi in lbd]
                    print('')
                    pmax_fit, lmax_fit, chi2 = fitSerk(lb, pBe, s, star=be, extens=extens)

                    csvout.writerow([be, fixName(be)+' (//)'] + list(pmax_fit) + list(lmax_fit) + [chi2, len(pBe)] + ['0']*3 + ['0','0','---',''] + ['0']*15)

                ll = np.linspace(3000.,8300.,100)
                pp = np.array([])
                for lli in ll:
                    pp = np.append(pp, polt.serkowski(pmax_fit[0], lmax_fit[0]*10000, lli, law=law, mode=2))
                # Only plot the graph if there are more than one data, because with an only point
                # the curve is not defined! But the emcee was runned to generate the covariance map
                if len(pBe) > 1:
                    ax.plot(ll, pp, color='black')



    # Copy the dummy lines beginning with the 'be' and fixed Be name
    if fit and usePrevious:
        copyFit(be, csvout, secondcol=fixName(be))

    if longname:
        ax.set_xlim([1500, 8500])

    leg = ax.legend(loc='best', borderaxespad=0., numpoints=1, prop={'size':fonts[4]})
    leg.get_frame().set_alpha(0.5)
    
    # Setting sizes
    ax.xaxis.label.set_fontsize(fonts[1])
    ax.yaxis.label.set_fontsize(fonts[1])
    for item in (ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fonts[2])

    if save:
        fig.savefig('{0}_p.{1}'.format(be,extens), bbox_inches='tight')
    else:
        fig.show()

    print('\nDone!\n')
    if fit:
        if os.path.exists('{0}_is.csv'.format(be)):
            os.remove('{0}_is.csv'.format(be))
        os.rename('{0}_is_tmp.csv'.format(be),'{0}_is.csv'.format(be))
        fout.close()
        
    return




def graf_pradial(csvfile, be, filt='pmax', vfilter=[], isfile=None, fit=False, \
                            bin_data=True, onlyY=True, save=False, extens='pdf', unbias='wk'):
    """
    Plot a field graph with polarization versus distance.

    filt : is the filter to plot in y axes - 'pmax','u','b',
           'v','r','i'. If 'pmax' use the Pmax values
           from ./'be'_is.csv file to plot.

    isfile : 'be'_is.csv file location (out from fs.graf_p).
              If None, it is suposed ./'be'_is.csv

    fit   : (only for filt=='pmax') fit a line in graph? This
            routine will not use in fitting the points whose
            values in 22th column inside isfile have a '0'
            character (this column defines the points to be
            used in the adjust). The points not used will be
            marked with a 'x' in the graph.

    csvfile : location of dados.csv.



    Considerations:

    - If filt=='pmax', don't plot the data whose values in 21th
    column inside isfile have a '0' character.

    - If filt=='pmax', don't use in fitting the points whose
    values in 22th column inside isfile have a '0' character
    when fit=True.  The points not used will be marked with a
    'x' in the graph.

    - Skip the lines inside isfile commented with a '#', or with
    a void first element
    
    """

    # Verify if vfilter is a special filter
    if vfilter in polt.vfil.keys():
        vfilter = polt.vfil[vfilter]

    plt.close('all')

    # Read csvfile and get the table
    data = readcsv(csvfile, be)
    if data == []:
        print('No {0} data!'.format(be))
        return
    objarr, plxarr, bearr, parr = getTable(data, 'plx', 'plxbe', sx='splx',
                    sy='splxbe', z='p', sz='s', vfilter=vfilter, bin_data=bin_data, onlyY=onlyY, unbias=unbias)
    if objarr == [] or objarr == [[],[],[]]:
        print('No {0} valid data!'.format(be))
        return

    # Verify isfile
    if filt == 'pmax':
        if isfile == None:
            isfile = '{0}_is.csv'.format(be)
        if not os.path.exists(isfile):
            print('No {0} file found!'.format(isfile))
            return

    # Initialize graphs
    fig = plt.figure(1)
    ax = plt.subplot(1, 1, 1)
    ax.set_title('Field Stars - {0}'.format(phc.bes[be]), \
                                        fontsize=fonts[0], verticalalignment='bottom')
    ax.set_xlabel(r'r (pc)', size=fonts[1])
    if filt=='pmax':
        ax.set_ylabel(r'$P_{max}$ (%)', size=fonts[1])
    else:
        ax.set_ylabel(r'$P_{0}$ (%)'.format(filt.upper()), size=fonts[1])

    # Extract data only in filter filt
    obj, x, y, y0, colors = [], [[],[]], [[],[],[]], [[],[],[]], []
    obj_filt, x_filt, y_filt, colors_filt = [], [[],[]], [[],[],[]], []
    x0 = [float(bearr[0][0]), float(bearr[1][0])]
    longname, useinfit = False, False

    
    # CASE 1)  filt=='pmax'
    if filt=='pmax':

        with open(isfile, 'r') as fr:
            csvread = csv.reader(fr, delimiter=';')#, quoting=csv.QUOTE_NONE, quotechar='')
            for i, line in enumerate(csvread):

                plotpointi = line[14]
                useinfiti = line[15]
                if line[0] == '' or line[0][0] == '#':
                    continue

                if line[0] != be and plotpointi=='1':
                    # 1) Filtered data (only works when fit==True)
                    if fit and useinfiti == '0':
                        # read paralaxes from csvfile
                        for j, obs in enumerate(objarr[0]):
                            if line[0]==obs:
                                x_filt[0] += [float(plxarr[0][j])]
                                x_filt[1] += [float(plxarr[1][j])]
                                break
                        obj_filt += [line[1]]
                        y_filt[0] += [float(line[2])]
                        y_filt[1] += [float(line[3])]
                        y_filt[2] += [float(line[4])]
                        colors_filt += [gen_color(csvfile, be, line[0], onlyY=onlyY)]
                    # b) Main data
                    else:
                        # read paralaxes from csvfile
                        for j, obs in enumerate(objarr[0]):
                            if line[0]==obs:
                                x[0] += [float(plxarr[0][j])]
                                x[1] += [float(plxarr[1][j])]
                                break
                        obj += [line[1]]
                        y[0] += [float(line[2])]
                        y[1] += [float(line[3])]
                        y[2] += [float(line[4])]
                        colors += [gen_color(csvfile, be, line[0], onlyY=onlyY)]

                    if len(line[1]) > 13:
                        longname = True

                # Don't plot Be's // component, but plot Be entire polarization
                # if it was writen in a line and the line[2] element is just
                # the Be fixed name + there are a '1' character in 'plot point'
                # column
                elif line[0] == be and plotpointi=='1':
                    if useinfiti == '1':
                        useinfit = True
                    obj0 = line[1]
                    x0 = [float(bearr[0][0]), float(bearr[1][0])]
                    y0 = [float(line[2]), float(line[3]), float(line[4])]
                    if len(obj0) > 13:
                        longname = True


    else:
        for i in range(len(objarr[0])):
            if objarr[1][i] == filt:
                try:
                    # Skip if plx is less than 0
                    if float(plxarr[0][i]) > 0:
                        x[0] += [float(plxarr[0][i])]
                        x[1] += [float(plxarr[1][i])]
                        y[0] += [float(parr[0][i])]
                        y[1] += [float(parr[1][i])]
                        y[2] += [float(parr[1][i])]
                        obj += [fixName(objarr[0][i])]
                        colors += [gen_color(csvfile, be, objarr[0][i], onlyY=onlyY)]
                        if len(obj[-1]) > 13:
                            longname = True
                except:
                    continue


    # Plot data
    for i in range(len(obj)):
        ax.errorbar(x[0][i], y[0][i], xerr=x[1][i], yerr=[[y[2][i]],[y[1][i]]], label=obj[i], \
                                            linestyle='--', marker='o', color=colors[i])
    for i in range(len(obj_filt)):
        ax.errorbar(x_filt[0][i], y_filt[0][i], xerr=x_filt[1][i], yerr=[[y_filt[2][i]],[y_filt[1][i]]], \
                                     label=obj_filt[i], linestyle='--', marker='x', color=colors_filt[i],\
                                     markersize=10, markeredgewidth=1.5)
    if filt=='pmax' and y0 != [[],[],[]]:
        ax.errorbar(x0[0], y0[0], xerr=x0[1], yerr=[[y0[2]],[y0[1]]], label=obj0, \
                                            linestyle='--', marker='s', color='black')

    # Fit the line
    if fit and filt=='pmax':
        # If it is to use the point of Be data in the ajust:
        if useinfit:
            x[0] += [x0[0]]
            x[1] += [x0[1]]
            y[0] += [y0[0]]
            y[1] += [y0[1]]
            y[2] += [y0[2]]

        # Fit by the total least squares method (orthogonal distance regression) without clipping
        param, sparam, cov, chi2, niter,bolfilt = phc.fit_linear(x[0], y[0], x[1], [(y[1][i]+y[2][i])/2 for i in range(len(y[0]))], clip=False)

        if len(y[0]) > 2:
            rchi2 = chi2[0]/(len(y[0])-2)
        else:
            rchi2 = 0

        pmaxfit = param[0]*x0[0]+param[1]
        spmaxfit = np.sqrt((param[0]*x0[1])**2 + (sparam[0]*x0[0])**2 + sparam[1]**2)
        
        print(55*'-')
        print('  Total least squares fit  (y = a*x+b):')
        print(55*'-')
        print('             a = {0:.3f} +- {1:.3f}'.format(param[0], sparam[0]))
        print('             b = {0:.3f} +- {1:.3f}'.format(param[1], sparam[1]))
        print('')
        print('                 N = {0:d}'.format(len(y[0])))
        print('         red chi^2 = {0:2f}'.format(rchi2))
        print('')
        print('')
        print('  Extrapolated Pmax: ')
        print('')
        print('          Pmax = {0:.4f} +- {1:.4f}'.format(pmaxfit, spmaxfit))
        print(55*'-')
        print('')

        xadj = np.linspace(ax.get_xlim()[0],ax.get_xlim()[1],3)
        yadj = param[0]*xadj+param[1]
        ax.plot(xadj, yadj, '-', color='dimgray', linewidth=1.5)



    # Fix limits
    if ax.get_xlim()[1] < x0[0]:
        ax.set_xlim([ax.get_xlim()[0], x0[0]+x0[1]])
    elif ax.get_xlim()[0] > x0[0]:
        ax.set_xlim([x0[0]-x0[1],ax.get_xlim()[1]])
    ax.autoscale(False)

    # Plot marks for Be star
    ax.plot([x0[0], x0[0]], ax.get_ylim(), linestyle ='--', color='gray')
    ax.plot([x0[0]-x0[1], x0[0]-x0[1]], ax.get_ylim(), \
                                    linestyle =':', color='gray')
    ax.plot([x0[0]+x0[1], x0[0]+x0[1]], ax.get_ylim(), \
                                    linestyle =':', color='gray')

    leg = ax.legend(loc='best', borderaxespad=0., numpoints=1, prop={'size':fonts[3]})
    if leg != None:
        leg.get_frame().set_alpha(0.5)

    # Setting sizes
    ax.xaxis.label.set_fontsize(fonts[1])
    ax.yaxis.label.set_fontsize(fonts[1])
    for item in (ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fonts[2])

    if save:
        plt.savefig('{0}_radial_{1}.{2}'.format(be,filt,extens), bbox_inches='tight')
    else:
        plt.show()





def graf_inst(logfile, mode=1, vfilter=['no-std'], save=False, extens='pdf'):
    """
    Plot a QU diagram for unpolarized standard stars in the
    logfile .log file (the outfile from polt.genTarget, appended by
    the name of the object at each line). Propagates error of
    standard star.

    mode=1 plot BVRI graphs in the same figure;
    mode=2 plot UBVRI graphs in separeted figures.

    Return 3 lists:

    1) qarr = [[mean Q], [propagated Q error], [Q stddev]]
    2) uarr = [[mean U], [propagated U error], [U stddev]]
    3) narr = [n] (the number of data used to compute the averages)

    Where [mean Q], [propagated Q error], ..., [n] are lists with
    the values for each filter.
    
    """


    # IMPLEMENTAR PESOS PARA CADA ALVO CONFORME ALGUM CRITÉRIO
    def meanQU(lines, filt):
        """
        Return the mean QU lists and the number of lines (n):
        [mean Q, propagated Q error, Q stddev],
        [mean U, propagated U error, U stddev],
        n

        lines: the lines readed from .logfile to be computed
        
        """

        p, q, u, s, thet, sdth = [],[],[],[],[],[]
        for line in lines:
            if line[3] == filt:
                p += [float(line[8])]
                q += [float(line[9])]
                u += [float(line[10])]
                thet += [float(line[11])]
                s += [float(line[12])]
                sdth += [float(line[7])]

        if q == []:
            print('No {0} data!'.format(filt))
            return [0,0,0],[0,0,0],0

        # Propagate error of delta_theta to Stokes QU
        lixo, sq, su = polt.propQU(p, thet, s, sdth)

        # List with mean QU, the propagated error and the stddev
        meanq = [np.mean(q), np.sqrt(np.dot(sq,sq))/len(sq), np.std(q)/np.sqrt(len(q))]
        meanu = [np.mean(u), np.sqrt(np.dot(su,su))/len(su), np.std(u)/np.sqrt(len(u))]

        return meanq, meanu, len(q)



    def plotQU(filt, fig, ax):
        """
        Receive figure and axes objects and do the plot for filter
        filt WITHOUT show or save the image.
        Return the same than meanQU().
        """

        # Factor to fix the font sizes
        if mode==1:
            factor=0.7
        else:
            factor=1.
            
        try:
            lines = np.loadtxt(logfile, dtype=str)
        except:
            print('# ERROR: Can\'t read file {0}.'.format(logfile))
            raise SystemExit(1)
        
#        ax.set_title('{0} filter'.format(filt.upper()), fontsize=fonts[0]*factor, verticalalignment='bottom')
        ax.text(0.98, 0.9, '{0} filter'.format(filt.upper()), horizontalalignment='right', \
                 verticalalignment='bottom', transform=ax.transAxes, fontsize=fonts[1]*factor)
    #    ax.set_xlabel(r'RA $-$ RA${}_{Be}$ (degree)', size=fonts[1])
    #    ax.set_ylabel(r'DEC $-$ DEC${}_{Be}$ (degree)', size=fonts[1])
        ax.set_xlabel(r'Q (%)', size=fonts[1]*factor)
        ax.set_ylabel(r'U (%)', size=fonts[1]*factor)


        # Do the subplots when flag is not 'E', tag is not in vfilter and filter is equal to filt
        j = 0
        while True:
            pts, spts = [[],[]], [[],[]]
            pts[0] = [float(line[9]) for line in lines if lines[j][-1] == line[-1] and \
                    line[3] == filt and line[16] != 'E' and not any(sub in line[17] for sub in vfilter)]
            pts[1] = [float(line[10]) for line in lines if lines[j][-1] == line[-1] and \
                       line[3] == filt and line[16] != 'E' and not any(sub in line[17] for sub in vfilter)]
            ptmp = [float(line[8]) for line in lines if lines[j][-1] == line[-1] and \
                    line[3] == filt and line[16] != 'E' and not any(sub in line[17] for sub in vfilter)]
            thettmp = [float(line[11]) for line in lines if lines[j][-1] == line[-1] and \
                    line[3] == filt and line[16] != 'E' and not any(sub in line[17] for sub in vfilter)]
            s = [float(line[12]) for line in lines if lines[j][-1] == line[-1] and \
                      line[3] == filt and line[16] != 'E' and not any(sub in line[17] for sub in vfilter)]
            sdth = [float(line[7]) for line in lines if lines[j][-1] == line[-1] and \
                    line[3] == filt and line[16] != 'E' and not any(sub in line[17] for sub in vfilter)]

            lixo, spts[0], spts[1] = polt.propQU(ptmp, thettmp, s, sdth)

            #print pts
            #print j, len(lines)
            
            if pts != [[],[]]:
                ax.errorbar(pts[0], pts[1], xerr=spts[0], yerr=spts[1], label=lines[j][-1], \
                        elinewidth=0.5, markersize=4, linestyle='', color='black', marker='o', alpha=0.7)

            nn = 0
            for line in lines[j:]:
                if line[-1] == lines[j][-1]:
                    nn += 1
            j += nn
            if j >= len(lines):
                break

        meanq, meanu, n = meanQU(lines, filt)
        if n == 0:
            return meanq, meanu, n

        print('FILTER {0}  ->  N = {1}'.format(filt.upper(), n))
        print('FILTER {0}  ->  Q (%): mean, error, stddev  =  {1:.7f}, {2:.7f}, {3:.7f}'\
                                        .format(filt.upper(), meanq[0], meanq[1], meanq[2]))
        print('FILTER {0}  ->  U (%): mean, error, stddev  =  {1:.7f}, {2:.7f}, {3:.7f}'\
                                        .format(filt.upper(), meanu[0], meanu[1], meanu[2]))

        coords = [[meanq[0]-meanq[2], meanu[0]-meanu[2]]]
        coords += [[meanq[0]+meanq[2], meanu[0]-meanu[2]]]
        coords += [[meanq[0]+meanq[2], meanu[0]+meanu[2]]]
        coords += [[meanq[0]-meanq[2], meanu[0]+meanu[2]]]
        polygon = Polygon(coords, True, color='blue', alpha=0.6, visible=True, fill='wheat')
        ax.add_patch(polygon)
    #    ax.errorbar(meanq[0], meanu[0], xerr=meanq[2], yerr=meanu[2], linestyle='', marker='o', \
    #                  elinewidth=2, fillstyle='full', markersize=8, color='black')
        ax.errorbar(meanq[0], meanu[0], xerr=meanq[2], yerr=meanu[2], linestyle='', marker='o', \
                      elinewidth=0.6, fillstyle='full', markersize=3, color='black')

        # Fix limits
        ax.autoscale(False)
        ax.set_xlim([-0.05, 0.05])
        ax.set_ylim([-0.05, 0.05])
        ax.plot(ax.get_xlim(), [0,0], 'k--')
        ax.plot([0,0], ax.get_ylim(), 'k--')
        ax.plot(ax.get_xlim(), [meanu[0],meanu[0]], ':', color='grey')
        ax.plot([meanq[0],meanq[0]], ax.get_ylim(), ':', color='grey')

        # Setting sizes
        ax.xaxis.label.set_fontsize(fonts[1]*factor)
        ax.yaxis.label.set_fontsize(fonts[1]*factor)
        for item in (ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(fonts[2]*factor)

        return meanq, meanu, n



    # Verify if vfilter is a special filter
    if vfilter in polt.vfil.keys():
        vfilter = polt.vfil[vfilter]

    plt.close('all')
    qarr, uarr, narr = [[],[],[]], [[],[],[]], []

    if mode==1:
        fig = plt.figure(1)
        axs = [plt.subplot(2, 2, 1)]
        axs += [plt.subplot(2, 2, 2, sharey=axs[0])]
        axs += [plt.subplot(2, 2, 3, sharex=axs[0])]
        axs += [plt.subplot(2, 2, 4, sharex=axs[1], sharey=axs[2])]
        plt.subplots_adjust(hspace=0.05, wspace=0.05)

        nax = 0
        for filt in ('b','v','r','i'):
            meanq, meanu, n = plotQU(filt, fig, axs[nax])
            for i in range(3):
                qarr[i] += [meanq[i]]
                uarr[i] += [meanu[i]]
            narr += [n]
            nax += 1

            print('')

#            axs[0].set_yticks(axs[0].get_yticks()[1:])
#            axs[3].set_xticks(axs[3].get_xticks()[1:])
#            print axs[3].get_xticks()
        plt.setp(axs[0].get_xticklabels(), visible=False)
        plt.setp(axs[1].get_xticklabels(), visible=False)
        plt.setp(axs[1].get_yticklabels(), visible=False)
        plt.setp(axs[3].get_yticklabels(), visible=False)
        axs[0].set_xlabel('')
        axs[1].set_xlabel('')
        axs[1].set_ylabel('')
        axs[3].set_ylabel('')

        if save:
            plt.savefig('qu_inst.{0}'.format(extens), bbox_inches='tight')
        else:
            plt.show(block=False)
            
    elif mode==2:
        nfig = 1
        for filt in ('u','b','v','r','i'):
            fig = plt.figure(nfig)
            ax = plt.subplot(1, 1, 1)
            meanq, meanu, n = plotQU(filt, fig, ax)
            for i in range(3):
                qarr[i] += [meanq[i]]
                uarr[i] += [meanu[i]]
            narr += [n]
            nfig += 1

            print('')
            if n == 0:
                continue

            if save:
                plt.savefig('qu_inst_{0}.{1}'.format(filt,extens), bbox_inches='tight')
            else:
                plt.show(block=False)


    return qarr, uarr, narr




def genAll(csvfile, path=None, genlogs=True, genint=True, vfilter=['no-std'], vfilter_graf_p=[], extens='pdf'):

    bin_data=True
    onlyY=True
    rotate=False  # Here
    every=False
    mcmc=True    # Here
    odr=True
    save=True

    if path == None or path == '.':
        path = os.getcwd()
        
    if extens == 'eps':
        print('Warning: field graphs will lose the transparency at .eps format')

    try:
        objs = np.loadtxt('{0}/refs/pol_alvos.txt'.format(hdtpath()), dtype=str)
    except:
        print('# ERROR: Can\'t read files pyhdust/refs/pol_alvos.txt.')
        raise SystemExit(1)

    # Generating logfiles for all Be stars
    if genlogs:
        for star in objs:
            print('Generating logfile for star {0}...'.format(star))
            polt.genTarget(star, path=path, ispol=None, skipdth=False, delta=3.5)

    # Generating thet_int.csv file and QU graphs
    if genint:
        for star in objs:
            print('Generating QU graphs for star {0}...'.format(be))
            genInt(star, path, vfilter=vfilter, extens=extens)
        
    for star in objs:
        print('='*50)
        print('Generating graphs for star {0}...'.format(star))
        graf_p(csvfile, star, rotate=rotate, path=path, bin_data=bin_data, onlyY=onlyY,
               save=save, fit=mcmc, extens=extens, vfilter=vfilter_graf_p)
        graf_pradial(csvfile, star, 'v', bin_data=bin_data, onlyY=onlyY, save=save, extens=extens, vfilter=vfilter)
        graf_field(csvfile, star, bin_data=bin_data, onlyY=onlyY, save=save, extens=extens)
        graf_theta(csvfile, star, bin_data=bin_data, onlyY=onlyY, save=save, extens=extens)
        if os.path.exists(path+'/'+star+'.log'):
            if not genint:
                arr_u, arr_b, arr_v, arr_r, arr_i = polt.graf_qu('{0}/{1}.log'.format(path,star),
                                                    mcmc=mcmc, odr=odr, save=True, extens=extens)
            polt.graf_t(path+'/'+star+'.log', save=save, extens=extens, vfilter=vfilter)
        print('\n\n')




def genInt(be, path=None, vfilter=['no-std'], extens='pdf'):
    """
    Call polt.graf_qu() for Be star 'be' and save the intrinsic
    angles inside thet_int.csv.

    path: path inside which the logfile 'be'.log (out from
          polt.genTarget) is located. The code will try to open
          thet_int.csv inside the current directory (.). If it
          was not found, it will be created. Otherwise, the
          routine will append a new line inside it, asking about
          to overwrite an eventual existing line for star 'be'.


    CONSIDERATIONS:

    - Propagates errors from standard star.
    
    """

    if path == None or path == '.':
        path = os.getcwd()

    if not os.path.exists('{0}/{1}.log'.format(path,be)):
        print('# No {0}/{1}.log file found.'.format(path,be))
        return 

    fw = open('thet_int_tmp.csv', 'w')
    csvwrite = csv.writer(fw, delimiter=';')#, quoting=csv.QUOTE_NONE, quotechar='')

    if os.path.exists('thet_int.csv'):
        fr = open('thet_int.csv', 'r')
        csvread = csv.reader(fr, delimiter=';')#, quoting=csv.QUOTE_NONE, quotechar='')
        opt = ''
        for i, line in enumerate(csvread):
            if line[0] == be:
                while opt not in ('y','Y','n','N'):
                    print('# Star {0} is already present inside file thet_int.csv. Run again and overwrite the current values?'.format(be))
                    opt = raw_input('(y/n): ')
                if opt in ('n','N'):
                    fr.close()
                    os.remove('thet_int_tmp.csv')
                    return
                else:
                    print('\n')
                    continue
            else:
                csvwrite.writerow(line)
        fr.close()
    else:
        csvwrite.writerow(['#obj','name']+['th_int', 'sth_int_+', 'sth_int_-', 'n']*5 +\
                          ['th_int', 'sth_int_+',' sth_int_-', 'comment'] +\
                          ['b*costh', 'sb*costh_+', 'sb*costh_-']*5 +\
                          ['Second peaks: ', 'repeated 7by7', 'cells following','these labels in', 'subheader below'])
        csvwrite.writerow(['#','']+['u']*4+['b']*4+['v']*4+['r']*4+['i']*4+[' ']*4 +\
                          ['u']*3+['b']*3+['v']*3+['r']*3+['i']*3 +\
                          ['filter', 'th_int', 'sth_int_+', 'sth_int_-','b*costh', 'sb*costh_+', 'sb*costh_-'])#+['to use']*4)
                
    arr_u, arr_b, arr_v, arr_r, arr_i = polt.graf_qu('{0}/{1}.log'.format(path,be),
                                            mcmc=True, odr=False, save=True, extens=extens, Vb_ran=[0., 1.])

#        arr_u[0][0] = np.arctan(arr_u[0][0])*90/np.pi
#        arr_u[1][0] = (90*arr_u[1][0])/(np.pi*(arr_u[1][0]**2+1))
#        arr_u[2][0] = arr_u[1][0]

    # Reshape lists. Copy the second peak informations to 'addpeak' plain list
    # and keep only the first peek informations inside arr_u, arr_b, etc.
    addpeak=[]
    arrs = (arr_u,arr_b,arr_v,arr_r,arr_i)
    for i,arr in enumerate(arrs):
        # Case one of filters didn't have been fitted
        if arr[0] == []:
            for j in (0,1,2):
                arrs[i][j] = [0,0,0,0,0]
        if len(arr[0])==2 and type(arr[0][0])==list:
            addpeak += [filters[i]]
            for j in (0,1,2):
                addpeak += [arr[j][1][0]/2]
            for j in (0,1,2):
                addpeak += [arr[j][1][1]]
                arrs[i][j] = arrs[i][j][0]

    # The operation '/2' below is because the intrinsic angle is the inclination angle of the
    # line in in QU diagram diveded by 2! (because PA = 1/2*arctan(U/Q))
    csvwrite.writerow([be, fixName(be)]+[arr_u[0][0]/2]+[arr_u[1][0]/2]+[arr_u[2][0]/2]+[arr_u[3]]+\
                                        [arr_b[0][0]/2]+[arr_b[1][0]/2]+[arr_b[2][0]/2]+[arr_b[3]/2]+\
                                        [arr_v[0][0]/2]+[arr_v[1][0]/2]+[arr_v[2][0]/2]+[arr_v[3]/2]+\
                                        [arr_r[0][0]/2]+[arr_r[1][0]/2]+[arr_r[2][0]/2]+[arr_r[3]/2]+\
                                        [arr_i[0][0]/2]+[arr_i[1][0]/2]+[arr_i[2][0]/2]+[arr_i[3]/2]+\
                                        ['','','','---']+\
                                        [arr_u[0][1]]+[arr_u[1][1]]+[arr_u[2][1]]+\
                                        [arr_b[0][1]]+[arr_b[1][1]]+[arr_b[2][1]]+\
                                        [arr_v[0][1]]+[arr_v[1][1]]+[arr_v[2][1]]+\
                                        [arr_r[0][1]]+[arr_r[1][1]]+[arr_r[2][1]]+\
                                        [arr_i[0][1]]+[arr_i[1][1]]+[arr_i[2][1]]+\
                                        addpeak)

    # Refresh the csv file
    print('\nDone!\n')
    fw.close()
    try:
        os.remove('thet_int.csv')
    except:
        pass
    os.rename('thet_int_tmp.csv','thet_int.csv')

    return



def rotQU(q, u, sq, su, ang, sang):
    """
    Rotates lists/arrays q and u in QU diagram at an angle
    2*(ang +- sang) (clockwise).

    Look that if 'ang' is the mean polarization angle,
    this rotation will transpose all polarization to
    Q parameter. U parameter will have residual variations
    with respect to the 0.

    Return q_rot, u_rot, sq_rot, su_rot

    Todo: sometimes it's better don't use any error value
    for ang?

    """

#    fig = plt.figure(1)
#    ax = plt.subplot(1, 1, 1)

    # Rotates each QU value
    qRot, sqRot, uRot, suRot = [],[],[],[]
    for i in range(len(q)):
        rad = 2*ang*np.pi/180.
        srad = 2*sang*np.pi/180.
        qRot += [q[i]*np.cos(rad)+u[i]*np.sin(rad)]
        uRot += [-q[i]*np.sin(rad)+u[i]*np.cos(rad)]
        sqRot += [np.sqrt((uRot[i]*srad)**2+(sq[i]*np.cos(rad))**2+(su[i]*np.sin(rad))**2)]
        suRot += [np.sqrt((qRot[i]*srad)**2+(sq[i]*np.sin(rad))**2+(su[i]*np.cos(rad))**2)]

#    ax.errorbar(q,u,xerr=sq,yerr=su)
#    ax.errorbar(qRot,uRot,xerr=sqRot,yerr=suRot)
#    plt.show(fig)

    return qRot, uRot, sqRot, suRot




def rotQUBe(be, thetfile, path=None, every=False, vfilter=['no-std']):
    """
    Rotates the QU values for Be star 'be' in the intrinsic
    angle specified inside thetfile and computes the <Q'>
    and <U'> of parameters rotated, returning four lists:
    lbd (the wavelength values), <U'>, sigma U' and stddev U',
    with the values in each one of the filters UBVRI.


    'thetfile' : the location (path+filename) of thet_int.csv
                 file with the intrinsic angles (out from fs.gen).
    
    'path'     : the path where is located the log file for star
                 'be'. If None, is supposed inside '.'.

    'every'    : use one intrinsic angle for each one filter to
                 rotate them? If every=False makes all data to use
                 a mean value at the -4:-2 collums (22th to 24th)
                 from 'thetfile'.


    CONSIDERATIONS:

    - Propagates errors from standard star.
    

    The rotated parameters are just the polarization
    components perpendicular and parallel to the
    orientation of the disk.

    If every==True, uses the intrinsic angle of each filter to
    rotate its QU values (if some angle==0, skip this filter!);
    otherwise, use the same value specified in last 4 columns
    for every one of the 5 filters.

    """

    if path == None:
        path = os.getcwd()

    try:
        f0 = open(thetfile, 'ro')
        reader = csv.reader(f0, delimiter=';')
    except:
        print('# ERROR: Can\'t read file {0}.'.format(thetfile))
        raise SystemExit(1)

    try:
        lines = np.loadtxt('{0}/{1}.log'.format(path, be), dtype=str)
    except:
        print('# ERROR: {0}.log file not found inside {1}.'.format(be, path))
        return [],[],[],[]

#    fig = plt.figure(1)
#    ax = plt.subplot(1, 1, 1)

    # Read the intrinsic angles
    thet, sthet, comment = [],[],[]
    for line in reader:
        if line[0] == be:
            if not every:
                try:
                    thet = [float(line[22])]*5
                    sthet = [(float(line[23])+float(line[24]))/2]*5
                    comment = [line[25]]*5
                except:
                    print('# ERROR: Components [-4],[-3],[-2] in lines of {0} file seem don\'t exist or not be float type for Be star {1}.'.format(thetfile, be))
                    return [],[],[],[]
            else:
                try:
                    for i in range(2,19,4):
                        thet += [float(line[i])]
                        sthet += [(float(line[i+1])+float(line[i+2]))/2]
                        comment += ['']
                        if thet[-1] == 0:
                            print('# WARNING: No intrinsic angle defined to filter {0}.'.format(filters[i/4]))
                except:
                    print('# ERROR: Components in lines of {0} file seem not be float type.'.format(thetfile))
                    return [],[],[],[]

    if thet == []:
        print('# ERROR: Star {0} not found in {1} file.'.format(be, thetfile))
        return [],[],[],[]

    j = 0
    lbd, qqRot, uuRot, sqqRot, suuRot = [],[],[], [[],[]], [[],[]]

    # Getting the values for each filter
    for filt in filters:

        if thet[j] == 0:
            j += 1
            continue

        JD, p, q, u, s, thet, sdth = [],[],[],[],[],[],[]
        for line in lines:
            if line[3] == filt and line[16] != 'E' and not any(sub in line[17] for sub in vfilter):
                JD += [float(line[0])]
                p += [float(line[9])]
                q += [float(line[9])]
                u += [float(line[10])]
                thet += [float(line[11])]
                s += [float(line[12])]
                sdth += [float(line[7])]

        lixo, sq, su = polt.propQU(p, thet, s, sdth)

        if len(q) != 0:
            # Rotates each QU value
            qRot, uRot, sqRot, suRot = rotQU(q,u,sq,su,thet[j],sthet[j])

            # Computes the mean of rotated QU parameters, propagates the errors and compute stddev
            if len(q) != 0:
                lbd += [phc.lbds[filt]]
                qqRot += [np.mean(qRot)]
                uuRot += [np.mean(uRot)]
                sqqRot[0] += [np.sqrt(sum([el**2 for el in sqRot]))/len(sqRot)]
                sqqRot[1] += [np.std(qRot)/np.sqrt(len(qRot))]
                suuRot[0] += [np.sqrt(sum([el**2 for el in suRot]))/len(suRot)]
                suuRot[1] += [np.std(uRot)/np.sqrt(len(uRot))]

        j += 1

#    plt.show()

#    fig = plt.figure(1)
#    ax = plt.subplot(1, 1, 1)
#    ax.errorbar(lbd, uuRot,yerr=suuRot[1])
#    plt.show()


    return lbd, uuRot, suuRot[0], suuRot[1]




def fitSerk(larr, parr, sarr, star='', law='w82', n_burnin=300, n_mcmc=800, \
                                                n_walkers=120, extens='pdf'):
    """
        Fit Serkowski law using Markov Chain Monte Carlo
        from emcee code.
        The likelihood function (L) supposes Gaussian
        errors around the Pmax values:

        log(L) = -0.5*chi2 -0.5*sum(ln(2*pi*sy^2))


      INPUT:
        
          larr: array/list with lambda values
          parr: array/list with P values
          sarr: array/list with the sigma_P values
          star: star name to be printed in the graph and
                its filename. If it's a void str '', this
                routine give a random number to prevent
                overwriting data.
           law: what K value in Serkowski's law use?
                (see polt.serkowski).
      n_burnin: number of iterations for burning-in
        n_mcmc: number of iterations to run emcee
     n_walkers: number of walkers to map the posterior
                probabilities.
        extens: extension for the graph file


      OUTPUT: sorted like "pmax_fit, lmax_fit, chi2"

        pmax_fit: [Pmax, sPmax_+, sPmax_-], the Pmax value and
                  its errors (at right and left from it). Pmax
                  is the median of distribution probability and
                  sPmax_+, sPmax_- are the range within which
                  there are 68.3% of the points in such
                  distribution.
        lmax_fit: Idem, for lmax.
            chi2: Reduced chi2.

    """

    import emcee
    import triangle.nov
    from matplotlib.ticker import MaxNLocator


    def lnprob(params, x, y, sy):
        """
        Return the log of posterior probability (p_pos) in
        bayesian statistics for the parameters 'params'
        ([Pmax,lmax]) and the data poits x, y and sy
        (y error values).

        p_pos = L*p_prior (unless by a normalization constant),
        where L is the likelihood function and p_prior is the
        prior probability function.

        In our case, for gaussian and independent uncertainies,
        the log of likelihood:

        log(L) = -0.5*chi2 -0.5*sum(ln(2*pi*sy^2))

        Now, p_prior = constant for 'params' values inside the
        range defined by 'intervalos'; otherwise, it is 0.
        That is the only determination that we can do.

        So, p_pos = -0.5*chi2 -0.5*sum(ln(2*pi*sy^2)) or
        -inf case 'params' are out from the allowed range.
        """

        Pmax, lmax = params

        # Set prior ln prob
        if Pmax < intervalos[0][0] or Pmax > intervalos[0][1] or \
                     lmax < intervalos[1][0] or lmax > intervalos[1][1] :
            lnprior = -np.inf
        else:
            lnprior = 0.

        # Return posterior prob
        if not np.isfinite(lnprior):
            return -np.inf
        else:
            return lnprior -0.5*np.sum(((polt.serkowski(Pmax, lmax*10000, x*10000, law=law, mode=2) - y)/sy)**2 + np.log(2*np.pi*(sy**2)))



    def run_emcee(sampler, p0):
        """
        Run emcee.

        p0 is the initial positions for the walkers
        """

        print("Burning-in ...")
        pos, prob, state = sampler.run_mcmc(p0, n_burnin)
        sampler.reset()

        print("Running MCMC ...")
        pos, prob, state = sampler.run_mcmc(pos, n_mcmc, rstate0=state)

        #~ Print out the mean acceptance fraction. 
        af = sampler.acceptance_fraction
        print("Mean acceptance fraction:", np.mean(af))


        # The lines below were to compute the best fit parameters using the maximum value
        """
        #~ Get the index with the highest probability
        maxprob_idx = np.argmax(prob)
        minprob_idx = np.argmin(prob)

        #~ Get the best parameters and their respective stddev + chi2
        params_fit = pos[maxprob_idx]
        stddev_fit = [sampler.flatchain[:,i].std() for i in xrange(ndim)]
        samples_aux = map(lambda v: [v[0]-params_fit[0], v[1]-params_fit[1]], samples)
        if len(larr) == 2:
            chi1 = 1.
        else:
            chi1 = np.sum(((polt.serkowski(params_fit[0], params_fit[1]*10000, larr*10000, law=law, mode=2) - parr)/sarr)**2)/(len(larr)-2)

        # Computing errors accordind to the distance to the maximum value inside which
        # there are 68.3% of the data
        p_fit1 = np.array([], dtype=float)
        p_fit2 = np.array([], dtype=float)
        l_fit1 = np.array([], dtype=float)
        l_fit2 = np.array([], dtype=float)

        print('Please wait, computing errors...')
        for elem in samples_aux:
            if elem[0] >= 0:
                p_fit1 = np.append(p_fit1, [elem[0]])
            else:
                p_fit2 = np.append(p_fit2, [elem[0]])
            if elem[1] >= 0:
                l_fit1 = np.append(l_fit1, [elem[1]])
            else:
                l_fit2 = np.append(l_fit2, [elem[1]])

        # Caution, 65.85 = 100-34.15 !!!
        p_error = [np.percentile(p_fit1, 68.3), -np.percentile(p_fit2, 31.7)]
        l_error = [np.percentile(l_fit1, 68.3), -np.percentile(l_fit2, 31.7)]

        print(74*'-')
        print('Output')
        print(74*'-')
        print '  P_max = {0:.4f}  +{1:.4f}  -{2:.4f}'.format(params_fit[0],p_error[0],p_error[1],stddev_fit[0])
        print 'lbd_max = {0:.4f}  +{1:.4f}  -{2:.4f}'.format(params_fit[1],l_error[0],l_error[1],stddev_fit[1])
        print 'reduced chi2 = {0:.4f}'.format(chi1)
        print(74*'-')

        """

        ### 1) Compute the results using all interval
        print('Please wait, computing errors...')
        samples = sampler.chain[:, :, :].reshape((-1, ndim))
        p_mcmc, l_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(samples, [16.075, 50, 83.925], axis=0)))

        if len(larr) == 2:
            chi = 0.
        else:
            chi = np.sum(((polt.serkowski(p_mcmc[0], l_mcmc[0]*10000, larr*10000, law=law, mode=2) - parr)/sarr)**2)/(len(larr)-2)

        #~ Plot the graphs -- histogram, corner and convergence map
#        fighists = plot_samples_hist(sampler)
        plot_conv(sampler, [p_mcmc[0], l_mcmc[0]])
        fig = triangle.nov.corner(samples, title=fixName(star), \
#                            truths=[p_mcmc[0], l_mcmc[0]], \
#                            extents=[(p_range[0],l_range[0]),(p_range[1],l_range[1])], \
                             quantiles=[0.16075, 0.50, 0.83925], \
                             labels=['$P_{max}\,($%$)$', '$\lambda_{max}\,(\mu m)$'], \
                             verbose=False)
        fig.savefig('{0}_correl.{1}'.format(star,extens))
        fig1 = triangle.nov.corner(samples, title=fixName(star), \
#                            truths=[p_mcmc[0], l_mcmc[0]], \
#                            extents=[(p_range[0],l_range[0]),(p_range[1],l_range[1])], \
                             quantiles=[0.16075, 0.50, 0.83925], \
                             labels=['$P_{max}\,($%$)$', '$\lambda_{max}\,(\mu m)$'], \
                             verbose=False)
        fig1.show()


        #~ Print the output using all interval
        """ TBD """
        print(74*'-')
        print('Output')
        print(74*'-')
        print('  P_max = {0:.4f}  +{1:.4f}  -{2:.4f}'.format(p_mcmc[0],p_mcmc[1],p_mcmc[2]))
        print('lbd_max = {0:.4f}  +{1:.4f}  -{2:.4f}'.format(l_mcmc[0],l_mcmc[1],l_mcmc[2]))
        print('reduced chi2 = {0:.4f}'.format(chi))
        print(74*'-')



        ### 2) NEW: Requests if the user want to use some specific interval
        opt = ''
        while opt not in ('y','Y','n','N'):
            print('These values were calculated using all Pmax and lbdmax data.\nDo you want' +\
                  ' to select specific ranges to use to compute the uncertainties?')
            opt = raw_input('(y/n): ')
        if opt in ('y','Y'):
            while True:
                print('')
                while True:
                    try:
                        petr = raw_input('Pmax: specify Pmax in format `Pmax_min,Pmax_max`: ')
#                        p_int = [float(ei)-params_fit[0] for ei in petr.split(',')]
                        p_range = [float(ei) for ei in petr.split(',')]
                        if len(p_range) == 2:
                            if p_range[1] > p_range[0]:
                                break
                            else:
                                print('Error: Pmax_max must be greather than Pmax_min!')
                        else:
                            print('Invalid input!')
                    except:
                        print('Invalid input!')

                while True:
                    try:
                        letr = raw_input('lbdmax: specify lbdmax in format `lbdmax_min,lbdmax_max`: ')
#                        l_int = [float(ei)-params_fit[1] for ei in letr.split(',')]
                        l_range = [float(ei) for ei in letr.split(',')]
                        if len(l_range) == 2:
                            if l_range[1] > l_range[0]:
                                break
                            else:
                                print('Error: lbdmax_max must be grather than lbdmax_min!')
                        else:
                            print('Invalid input!')
                    except:
                        print('Invalid input!')

                opt = ''
                while opt not in ('y','Y','n','N'):
                    print('\nIs it correct?     Pmax_min,Pmax_max = ' + petr + '\n' +\
                            '               lbdmax_min,lbdmax_max = ' + letr)
                    opt = raw_input('(y/n): ' )
                if opt in ('y','Y'):
                    break

            # The lines below were to compute the best fit parameters using the maximum value
            """
            p_fit1 = np.array([], dtype=float)
            p_fit2 = np.array([], dtype=float)
            l_fit1 = np.array([], dtype=float)
            l_fit2 = np.array([], dtype=float)
            print('Please wait, computing errors...')

            # fit1: arrays for the right of the best values
            # fit2: arrays for the left of the best values
            for elem in samples_aux:
                if elem[0] >= 0 and elem[0] < p_int[1]:
                    p_fit1 = np.append(p_fit1, [elem[0]])
                elif elem[0] < 0 and elem[0] > p_int[0]:
                    p_fit2 = np.append(p_fit2, [elem[0]])
                if elem[1] >= 0 and elem[0] < l_int[1]:
                    l_fit1 = np.append(l_fit1, [elem[1]])
                elif elem[1] < 0 and elem[1] > l_int[0]:
                    l_fit2 = np.append(l_fit2, [elem[1]])

            # Caution, 100-68.3 = 31.7!
            p_error = [np.percentile(p_fit1, 68.3), -np.percentile(p_fit2, 31.7)]
            l_error = [np.percentile(l_fit1, 68.3), -np.percentile(l_fit2, 31.7)]

            print(74*'-')
            print('Output')
            print(74*'-')
            print '  P_max = {0:.4f}  +{1:.4f}  -{2:.4f}'.format(params_fit[0],p_error[0],p_error[1],stddev_fit[0])
            print 'lbd_max = {0:.4f}  +{1:.4f}  -{2:.4f}'.format(params_fit[1],l_error[0],l_error[1],stddev_fit[1])
            print 'reduced chi2 = {0:.4f}'.format(chi1)
            print(74*'-')
            """

            # Filtering 'samples' array
            print('Please wait, computing errors...')
            samples_new = np.empty(shape=[0, 2])
            for elem in samples:
                if elem[0] >= p_range[0] and elem[0] <= p_range[1] and \
                   elem[1] >= l_range[0] and elem[1] <= l_range[1]:
                    samples_new = np.vstack([samples_new, elem])

            # Computing NEW medians and errors
            p_mcmc, l_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                                 zip(*np.percentile(samples_new, [16.075, 50, 83.925], axis=0)))

            if len(larr) == 2:
                chi = 0.
            else:
                chi = np.sum(((polt.serkowski(p_mcmc[0], l_mcmc[0]*10000, larr*10000, law=law, mode=2) - parr)/sarr)**2)/(len(larr)-2)

            
            #~ Print the output using the specific range
            """ TBD """
            print(74*'-')
            print('Output')
            print(74*'-')
            print('  P_max = {0:.4f}  +{1:.4f}  -{2:.4f}'.format(p_mcmc[0],p_mcmc[1],p_mcmc[2]))
            print('lbd_max = {0:.4f}  +{1:.4f}  -{2:.4f}'.format(l_mcmc[0],l_mcmc[1],l_mcmc[2]))
            print('reduced chi2 = {0:.4f}'.format(chi))
            print(74*'-')

            # Save the new triangle graph
            fig = triangle.nov.corner(samples_new, title=fixName(star), \
#                                 truths=[p_mcmc[0], l_mcmc[0]], \
#                                 extents=[(p_range[0],l_range[0]),(p_range[1],l_range[1])], \
                                  quantiles=[0.16075, 0.50, 0.83925], \
                                  labels=['$P_{max}\,($%$)$', '$\lambda_{max}\,(\mu m)$'], \
                                  verbose=False)
            fig.savefig('{0}_correl_cut.{1}'.format(star,extens))

        else:
            try:
                os.remove('{0}_correl_cut.{1}'.format(star,extens))
            except:
                pass

#        plt.close(fighists[0])
#        plt.close(fighists[1])
        plt.close(fig1)

        return sampler, p_mcmc, l_mcmc, chi



    def plot_samples_hist(sampler):
        """
        Plot two figures with the histograms
        """
        samples = [sampler.flatchain[:,i] for i in (0,1)]
        par = ['$P_{max}$', '$\lambda_{max}$']

        fig = []
        for i, sample in enumerate(samples):
            fig += [plt.figure()]
            plt.hist(sample, 100)
            plt.title('Sample of parameter {0}'.format(par[i]))
            fig[-1].show()
            
        return fig



    def plot_conv(sampler, param):
        """
        Plot convergence map. 'param' are the values to be highlighted
        """

        fig, axes = plt.subplots(2, 1, sharex=True, figsize=(8, 6))
        axes[0].plot(sampler.chain[:, :, 0].T, color="k", alpha=0.4)
        axes[0].yaxis.set_major_locator(MaxNLocator(5))
        axes[0].axhline(param[0], color="#888888", lw=2)
        axes[0].set_ylabel("$P_{max}$")

        axes[1].plot(sampler.chain[:, :, 1].T, color="k", alpha=0.4)
        axes[1].yaxis.set_major_locator(MaxNLocator(5))
        axes[1].axhline(param[1], color="#888888", lw=2)
        axes[1].set_ylabel("$\lambda_{max}$")
        axes[1].set_xlabel("Step number")

        fig.tight_layout(h_pad=0.0)
        fig.savefig('{0}_conv.{1}'.format(star,extens))

        return



    # Setting parameters and limits
    Pmax_max = 9.
    Pmax_min = 0.
    lmax_max = 1.
    lmax_min = 0.
    intervalos = np.array([[Pmax_min, Pmax_max], [lmax_min, lmax_max]])
    ndim = 2

    # Converting lists to np.array
    if type(parr) == list:
        parr = np.array(parr)
    if type(larr) == list:
        larr = np.array(larr)
    if type(sarr) == list:
        sarr = np.array(sarr)

    # If 'star' was not specified, generate a random number to append to the graph name to be saved
    if star == '':
        star = 'rand' + str(int(np.random.rand(1)[0]*10000))

#    larr = np.array([0.3650, 0.4450, 0.5510, 0.6580, 0.8060])
#    parr = np.array([0.26685, 0.34856, 0.36904, 0.33956, 0.29710])
#    sigma = np.array([0.09753, 0.00881, 0.00749, 0.01132, 0.01120])

#    larr = np.array([0.4450, 0.5510, 0.6580, 0.8060])
#    parr = np.array([0.34856, 0.36904, 0.33956, 0.29710])
#    sarr = np.array([0.00881, 0.00749, 0.01132, 0.01120])

    # Define random values to be used as priori numbers within the interval
    p0 = np.array( [np.random.rand(ndim) for n in xrange(n_walkers)] )
    for k in range(ndim):
        p0[:,k] = intervalos[k][0]+p0[:,k]*(intervalos[k][1]-intervalos[k][0])

    # Initialize the sampler and run mcmc
    sampler = emcee.EnsembleSampler(n_walkers, ndim, lnprob, args=[larr, parr, sarr], a=3)#, threads=2)
    sampler, pmax_fit, lmax_fit, chi = run_emcee(sampler, p0)


    return pmax_fit, lmax_fit, chi


