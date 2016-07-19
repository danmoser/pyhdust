# -*- coding:utf-8 -*-

"""PyHdust *poltools* module: polarimetry tools

History:
-grafpol working for *_WP1110....log files!
-grafpol working for log/out files with more than a single star

:co-author: Daniel Bednarski
:license: GNU GPL v3.0 https://github.com/danmoser/pyhdust/blob/master/LICENSE 
"""
from __future__ import print_function
import os as _os
import re as _re
import pwd as _pwd
import time as _time
from glob import glob as _glob
import numpy as _np
import datetime as _dt
import shutil as _shutil
# import csv as _csv
# from itertools import product as _product
# from glob import glob as _glob
from inspect import getouterframes as _getouterframes
from inspect import currentframe as _currentframe
import pyhdust.phc as _phc
import pyhdust.jdcal as _jdcal
from pyhdust import hdtpath as _hdtpath
# from sys import _argv
# from matplotlib import rc as _rc
import sys as _sys


def eprint(*args, **kwargs):
    print(*args, file=_sys.stderr, **kwargs)
    return

try:
    import matplotlib.pyplot as _plt
    from matplotlib.transforms import offset_copy as _offset_copy
    import pyfits as _pyfits
    from scipy.optimize import curve_fit as _curve_fit
    from scipy.integrate import simps as _simps
    from scipy.interpolate import interp1d as _interp1d
    import matplotlib as _mpl
    _mpl.rcParams['pdf.fonttype']=42
except ImportError:
    eprint('# Warning! matplotlib and/or pyfits module not installed!!!')

__author__ = "Daniel Moser"
__email__ = "dmfaes@gmail.com"


filters = ['u','b','v','r','i']
fonts = [20, 17, 17, 14, 13]  # Font sizes for titles, axes labels, axes values, key label of graphs, subplot labels

# Setting an "initial value" for ccd
ccd = '---'

# Dictionary for the tags entered by the user
dictags = {0: ['bad modulation','bad-mod'],
           1: ['very bad modulation','very-bad-mod'],
           2: ['the pol values have values incompatible from each other','incomp-mods'],
           3: ['some observational problem/error','obs-prob'],
           4: ['polarimeter problem suspected','iagpol-prob'],
           5: ['another relevant problem','other-prob'],
          }

# Dictionary for the tags assigned automatically
# If you want to add another value, add inside verout routin also.
dictests = {0: ['std incompatible with the published','obs!=pub', 'W'],
            1: ['sig >> theorical_sig','s>>theor_s', 'OK'],
            2: ['no standard in the night','no-std', 'W'],
            3: ['standard from another day','oth-day-std', 'W'],
            4: ['delta theta estimated from another filter','oth-dth', 'W'],
           }

# Dictionary for pre-defined vfilters
vfil = { 'comp' : ['no-std','iagpol-prob','oth-day-std'],
         'prob' : ['no-std','iagpol-prob','incomp-mods','obs-prob','other-prob','oth-day-std'],
         'full' : ['no-std','iagpol-prob','incomp-mods','obs-prob','other-prob','oth-day-std','bad-mod','very-bad-mod'],
       }


#################################################
#################################################
#################################################
def stdchk(stdname):
    """
    Check if the standard star name contains a known name, and return
    its position in `padroes.txt`.
    """
    lstds = list(_np.loadtxt('{0}/refs/pol_padroes.txt'.format(_hdtpath()), dtype=str,\
    usecols=[0]))
    chk = False
    i = -1
    for std in lstds:
        if stdname.find(std) > -1:
            chk = True
            i = lstds.index(std)
    return chk, i



#################################################
#################################################
#################################################
def countStars(objdir, f):
    """
    Count how many stars there are inside outfiles in 'objdir'
    and filter f. Return 0 if there are no outfiles
    """
    louts = _glob('{0}/*_{1}_*.out'.format(objdir,f))
    if len(louts) == 0:
        counts = 0
    else:
        file0 = _np.loadtxt(louts[0], dtype=str, delimiter='\n', comments=None)
        counts = len(file0)-1    # -1 because the header line

    return counts



#################################################
#################################################
#################################################
def thtFactor(MJD):
    """
    Return the factor for polarization angle 'theta'. This factor
    indicates when theta must be taken as 180-theta (factor==-1)
    or +theta (factor==+1).

    It's based on engine of IAGPOL polarimeter.
    If MJD < 57082.5, return -1 (including MJD=-1, for no assigned
    value); otherwise, return +1.

    Theta out from polrap is correct for a WP rotating in
    'counter-clockwise' direction, then:
    
    factor = -1 when WP rotating in clockwise
    factor = +1 when WP rotating in counter-clockwise

    """

    if float(MJD) < 54101.5:  # Eu desconfio que antes de 2007. Confirmar a data exata
        factor=1.
    elif float(MJD) < 57082.5: # before 2015, March 1st
        factor=-1.
    else:
        factor=1.

    return factor



#################################################
#################################################
#################################################
def readout(out, nstar=1):
    """
    Read the *.out file from IRAF reduction and return a float array

    Q   U        SIGMA   P       THETA  SIGMAtheor.  APERTURE  STAR

    'nstar' == star number inside 'out' file (usefull when there are
               more than a single star inside .out)
    """
    f0 = open(out)
    data = f0.readlines()
    f0.close()
    data = data[nstar].split()
    return [float(x) for x in data]



#################################################
#################################################
#################################################
def readoutMJD(out, nstar=1):
    """
    Read the 'out' file from IRAF reduction in a float array (fout),
    appending the MJD date and the angle of the beams from
    calcite.

    'nstar' == star number inside 'out' file. PS: calcice angle
            is allways evaluated using the first star coordinates.
    """

    path = _phc.trimpathname(out)[0]
    outn = _phc.trimpathname(out)[1]
    try:
        data = readout(out, nstar=nstar)
    except:
        print('# ERROR: Can\'t open/read file {0}. Verify and run again.\n'.format(out))
        raise SystemExit(1)
        
    WP = False
    if '_WP' in outn:
        outn = outn[:outn.rfind('_')]
        WP = True

    i = outn.rfind('_')+1
    seq = int(outn[i:i+2])
    npos = int(outn[i+2:i+5])
    f = outn[outn.find('_')+1:outn.rfind('_')]
    JD = _glob('{0}/JD_*_{1}'.format(path,f))
    try:
        f0 = open(JD[0])
        date = f0.readlines()
        f0.close()
        datei = float(date[npos-1].split()[-1])-2400000.5
        datef = float(date[npos-1+seq-1].split()[-1])-2400000.5
    except:
        print(('# ERROR: Found *_{0}_*.out files, but none JD file found as {1}/JD_*_{2}. '+\
                            'Verify and run again.\n').format(f,path,f))
        raise SystemExit(1)

    if WP:
        ver = outn[-1]
    else:
        i = outn[:-4].rfind('.')+1
        ver = outn[i:i+1]

    coords = _glob('{0}/coord_*_{1}.{2}.ord'.format(path,f,ver))
    if len(coords) == 0 and ccd not in ('301','654'):
        coords = _glob('{0}/coord_*_{1}.ord'.format(path,f))
        if len(coords) == 0:
            coords = _glob('{0}/coord_*_{1}_[0-9]*.ord'.format(path,f))
            if len(coords) == 0:
                print(('# ERROR: Found *_{0}_*.out files, but none COORD file found as '+\
                        '{1}/coord_*_{2}_*.ord. Verify and run again.\n').format(f,path,f))
                raise SystemExit(1)

    try:
        if ccd not in ('301','654'):
            coords = _np.loadtxt(coords[0])
            ang = _np.arctan( (coords[1,1]-coords[0,1])/(coords[1,0]-coords[0,0]) )*180/_np.pi
        else:
            coords = _np.array([[0.,0.],[0.,0.]])
            ang = 0.
        while ang < 0:
            ang += 180.
    except:
        print('# ERROR: Can\'t open coords file {0}/coord_*_{1}_*.ord. Verify and run again.\n'.format(path,f))
        raise SystemExit(1)

    if date != -1:
        if datei == datef:
            print('# Strange JD file for '+out)
        date = (datef+datei)/2


    return [float(x) for x in data]+[date]+[ang]



#################################################
#################################################
#################################################
def chooseout(objdir, obj, f, nstar=1, sigtol=lambda sig: 1.4*sig):
    """
    Olha na noite, qual(is) *.OUT(s) de um filtro que tem o menor erro.

    Retorna um mais valores para toda a sequencia (i.e., pasta).
    Criterios definidos no anexo de polarimetria.

    minerror == True: recebe o out de menor erro em qualquer caso.

    sigtol: condicao de tolerancia para pegar o agrupamento de menor erro.
    Se o sigma do agrupamento com todas N posicoes for menor que a funcao
    sigtol sobre o sigma do agrupamento de menor erro, entao usa o agrupamento
    com as N posicoes; do contrario usa o de menor erro.
    Exemplo com 16 posicoes: erro do grupo de 16 posicoes == 0.000230;
        menor erro == 0.000200. Se sigtol(sig) = 1.4*sig, como
        sigtol(0.000200) == 0.000240 > 0.000230, usa o agrupamento de 16 posicoes.

    O numero de posicoes eh baseado no numero de arquivos *.fits daquele
    filtro.

    'nstar' == star number inside 'out' file (usefull when there are
               more than a single star inside .out)
    """


    def minErrBlk16(serie='16001'):
        """
        Calculate the out with best error out of type *_f_*serie.?.out
        for star number 'nstar' (where 'f' is the filter, 'serie' is
        the five-numbers concerning to the WP positions (like 16001)
        and ? is some char.

        Return err, out. If not found, return 1000.,''.
        """

        err = 1000.
        out = ''
        ls = [objdir+'/'+fl for fl in _os.listdir('{0}'.format(objdir)) if _re.search(r'_{0}'. \
                    format(f) + r'_.*_?{0}\..\.out'.format(serie), fl)]

        if len(ls) > 0:
            err = float(readout(ls[0],nstar=nstar)[2])
            out = ls[0]
            for outi in ls:
                if float(readout(outi,nstar=nstar)[2]) < err:
    #                print float(readout(outi,nstar=nstar)[2])
                    err = float(readout(outi,nstar=nstar)[2])
                    out = outi

        return err, out


    npos = len(_glob('{0}/*_{1}_*.fits'.format(objdir,f)))
    if npos == 0:
        npos = len(_glob('{0}/{1}/p??0'.format(objdir,f)))

    louts = _glob('{0}/*_{1}_*.out'.format(objdir,f))

    # Check reduction
    if len(louts) == 0 and npos != 0:
        print(('# ERROR: There are observations not reduced for {0}/{1}_{2}_*.fits. ' +\
                            'Verify and run again.\n').format(objdir,obj,f))
        raise SystemExit(1)

    # Calculate the number of outfiles to be returned.
    n=npos/16   # operacao em formato int!
    rest = npos%16
    if n!=0:
        if rest == 0:
            nlast = 16
        elif rest >= 8:
            n += 1
            nlast = rest
        elif rest < 8:
            nlast = 16+rest
    elif rest > 0:
        n = 1
        nlast = rest

#    print n, rest, nlast
    err = [1000.]*n
    outs = ['']*n
    # n contem o numero de outs que serao obtidos
    # nlast contem o numero de posicoes para o ultimo out

    # Loop to get the n outfiles
    for i in range(n):

        # Get the best outfile with all WP positions.
        if i+1 < n or (i+1==n and nlast >= 16):
            serie='{0:02d}{1:03d}'.format(16,i*16+1)
        else:
            serie='{0:02d}{1:03d}'.format(nlast,i*16+1)

        err[i], outs[i] = minErrBlk16(serie)
        errtmp = err[i]   # errtmp==1000 if no outfiles were found by minErrBlk16

        # Tests if there is some better group, with smaller error
        for outi in louts:
            if outi.find('_WP') == -1:
                indx = outi.rfind('_')+1
            else:
                indx = outi[:outi.find('_WP')].rfind('_')+1
            combi = outi[indx:indx+5]
            n1= int(combi[:2]) # First part of '16001', etc 
            n2= int(combi[2:]) # Last part of '16001', etc 
#            if i+1==n:
#                print n1, n2
            # Default case
            if i+1 != n or (i+1 == n and nlast == 16):
                # Get only the groups with independent data
                if n2  >= 16*i+1 and n2 <= 16*i+1 + (16-n1):
                    if float(readout(outi,nstar=nstar)[2]) < errtmp:
                        errtmp = float(readout(outi,nstar=nstar)[2])
                        outtmp = outi
            # Case i==n (and nlast!=16)
            else:
#                print 'entrou1'
#                print n1,n2,16*i+1
                if n2  >= 16*i+1:
                    if float(readout(outi,nstar=nstar)[2]) < errtmp:
                        errtmp = float(readout(outi,nstar=nstar)[2])
                        outtmp = outi

        if errtmp != err[i] and err[i] > sigtol(errtmp):
            outs[i] = outtmp

    # if some element of outs is '', chooseout has failed to find the best out in such 16-position serie.
    # But don't panic. It can happen due some espurious .fits file
#    print [out for out in outs if out != '']
    return [out for out in outs if out != '']



#################################################
#################################################
#################################################
def verout(out, obj, f, nstar=1, verbose=True, delta=3.5):
    """
    Function to do tests on outfile 'out' concerning to
    star number 'nstar', object name 'obj' in filter 'f'.
    Tests: (1) test if P_pub for standards is compatible with
               P_obs value within 10 sigmas (10 because
               std pol can changes with the time).
           (2) test if sig < 3*sig_theorical.
           (3) test if there are some standard star (only if
                'obj' is a target star)

    - If verbose==True, show warnings in screen
    - In objdir==None, outfile is supposed in current dir

    Return a boolean list with three components concerning
    to the tests (1)-(3) above + log string. If some test has failed,
    the concerning value is assigned as \"True\"; otherwise,
    \"False\".
    """

    tests = [False]*len(dictests)
    loglines = ''
    # The complex lines below is to extract 'path' from 'out' (considering) fixing consecutives '//'
    if out[0] == '/':
        path = '/'+'/'.join(s for s in [s for s in out.split('/') if s][:-2])
    else:
        path = '/'.join(s for s in [s for s in out.split('/') if s][:-2])

    [Q,U,sig,P,th,sigT,ap,star,MJD,calc] = readoutMJD(out, nstar=nstar)
    sig_ratio = float(sig)/float(sigT)
    ztest = verStdPol(obj, f, float(P)*100, float(sig*100))
    
    # Some tests.
    if ztest > 10.0:     # Case the object is not a standard, ztest==-1 and tests[0] remains False.
        tests[0] = True
    if sig_ratio > 6.:
        tests[1] = True
    if not stdchk(obj)[0]:  # Only if object is not a standard star, tests if there exists some standard star for it
        tests[2] = not chkStdLog(f, calc, path=path, delta=delta, verbose=False)
    
    # Print tests
    if tests[0]:
        loglines += ('# WARNING: {0}_{1}: The standard has polarization only compatible '+\
                                  'within {2:.1f} sigma with the published value.\n').format(obj, f, ztest)
    if tests[1]:
        loglines += ('# WARNING: {0}_{1}: Polarization has sig >> theorical_sig ' +\
                        '({2:.4f} >> {3:.4f}).\n').format(obj, f, sig*100, sigT*100)
    if tests[2]:
        loglines += ('# WARNING: {0}_{1}: Standard star not found yet '+\
                                             '(calc. {2:.1f})\n').format(obj, f, calc)

    if verbose and loglines != '':
        print('\n'+loglines)

    
    return tests, loglines



#################################################
#################################################
#################################################
def queryout(objdir, obj, f, nstar=1, sigtol=lambda sig: 1.4*sig):
    """
    Call chooseout for 'obj' at filter 'f' (in 'objdir'),
    print the graphs and query to the user if he wants
    to use the selected outs. If not, he musts answer
    what to use.

    'nstar' == star number inside 'out' file (usefull when there are
               more than a single star inside .out)

    """

    _plt.close('all')
    outs = chooseout(objdir, obj, f, nstar=nstar, sigtol=sigtol)
    if outs == ['']:
        return outs, None, None

    # Initialize the components for each outfile
    tags = [[]]*len(outs)
    flag = ['']*len(outs)

    for i in range(len(outs)):

        sortout = []
        while True:

            _plt.close('all')
            # Only in the first iteration stack the values in sortout
            if sortout == []:
                sortout = grafall(objdir, f, n=i+1 ,nstar=nstar, bestouts=[outs[i]], shortmode=True)
            else: 
                lixo = grafall(objdir, f, n=i+1 ,nstar=nstar, bestouts=[outs[i]], shortmode=True)

            print('\n'+'_'*80)
            print(' {0:<10s} {1:<5s} {2:<6s}  {3:<7s} {4:<10s} {5:<7s} {6:<s}'.format('Obj', 'Filt', \
                    'Pol (%)', '', 'sig/ThSig', 'ztest', 'out/num'))
            try:
                [Q,U,sig,P,th,sigT,ap,star,MJD,calc] = readoutMJD(outs[i], nstar=nstar)
                sig_ratio = float(sig)/float(sigT)
                z = verStdPol(obj, f, float(P)*100, float(sig*100))
                numout = '(#{0})'.format(sortout.index(outs[i]))
            except:
                print('# ERROR: It shouldn\'t enter here in queryout!\n')
                raise SystemExit(1)

            # Reassigns ztest value for printing
            if z == -1:
                zstr = '-----'
            else:
                zstr = '{0:.1f}'.format(z)

            # Prints the values
            print(' {0:<10s} {1:<5s} {2:<6.4f}+-{3:<7.4f} {4:<10.1f} {5:<7s} {6:<s} {7:<s}'.\
                     format(obj.upper(), f.upper(), float(P)*100, float(sig)*100, \
                            sig_ratio, zstr, _phc.trimpathname(outs[i])[1], numout))
            print('_'*80+'\n')

            # Test the out file to print tests
            testout, logout = verout(outs[i], obj, f, nstar=nstar, verbose=True)

            while True:
                verif = raw_input('Use this out? (y/n): ')
                if verif not in ('y','Y','n','N'):
                    print('Invalid choise!')
                else:
                    break

            if verif in ('y', 'Y'):
                break
            elif verif in ('n', 'N'):
                opt=''  # for the first iteration
                while True:
                    opt = raw_input('Type the out number: ')
                    # If opt is a valid value, assign the input number with the concerning out file
                    if opt in [str(j) for j in range(1,len(sortout)) if sortout[j] != '']:
                        outs[i] = sortout[int(opt)]
                        break
                    else:
                        print('Wrong value!')
                        opt=''

        # Request what tags to assign (flexible through the dictags global dictionary)
        print('\n# TAGS LIST\n  0: none')
        for j in dictags.keys():
            print('  {0}: {1}'.format(j+1, dictags[j][0]))
        print('')

        while True:
            verif = True
            tags[i] = [False for j in dictags.keys()]
            strin = raw_input('Select all tags that apply separated by commas (\'0\' for none): ')
            if strin == '0':
                flag[i]='OK'
                break
            opts = strin.split(',')
            for opt in opts:
                if opt in [str(j+1) for j in dictags.keys()]:
                    opt = int(opt)-1
                    tags[i][opt] = True
                else:
                    print('Invalid choise!')
                    verif = False
                    break

            # If some tag was selected, request a flag below
            if verif:
                verif2=''
                while verif2 not in ('y','Y','n','N'):
                    verif2 = raw_input('For you, this data should appear as usable? (y/n): ')
                if verif2 in ('y','Y'):
                    flag[i] = 'W'
                else:
                    flag[i] = 'E'
                break

    _plt.close('all')

    return outs, tags, flag



#################################################
#################################################
#################################################
# Falta alterar para novos índices das colunas dos arquivos std.dat e obj.dat
# das mudanças que fiz. Bednarski.
def plotfrompollog(path, star, filters=None, colors=None):
    """ Plot default including civil dates
    """
    tab = _np.genfromtxt('{0}/{1}.log'.format(path,star), dtype=str, autostrip=True)

    MJD = tab[:,0].astype(float)
    nights = tab[:,1]
    filt = tab[:,2]
    calc = tab[:,3]
    ang_ref = tab[:,4].astype(float)
    dth = tab[:,5].astype(float)
    P = tab[:,7].astype(float)
    Q = tab[:,8].astype(float)
    U = tab[:,9].astype(float)
    th = tab[:,10].astype(float)
    sigP = tab[:,11].astype(float)
    sigth = tab[:,12].astype(float)

    if colors == None:
        colors = _phc.colors
    if filters == None:
        filters = ['b','v','r','i']
        colors = ['b','y','r','brown']
    leg = ()
    fig, ax = _plt.subplots()
    for f in filters:
        i = [i for i,x in enumerate(filters) if x == f][0]
        leg += (f.upper()+' band',)
        ind = _np.where(filt == f)
        x = MJD[ind]
        y = P[ind]
        yerr = sigP[ind]
        ax.errorbar(x, y, yerr, marker='o', color=colors[i], fmt='--')
    ax.legend(leg,'upper left')#, fontsize='small')    
    #ax.legend(leg,'lower left', fontsize='small')  
    ax.set_ylabel('Polarization (%)')
    ax.plot(ax.get_xlim(),[0,0],'k--')
    ylim = ax.get_ylim()
    ax.set_xlabel('MJD')
    xlim = ax.get_xlim()
    ticks = _phc.gentkdates(xlim[0], xlim[1], 3, 'm',\
    dtstart=_dt.datetime(2012,7,1).date())
    mjdticks = [_jdcal.gcal2jd(date.year,date.month,date.day)[1] for date in \
    ticks]
    ax2 = ax.twiny()
    ax2.set_xlabel('Civil date')
    ax2.set_xlim(xlim)
    ax2.set_xticks(mjdticks)
    ax2.set_xticklabels([date.strftime("%d %b %y") for date in ticks])
    _plt.setp( ax2.xaxis.get_majorticklabels(), rotation=45 )
    _plt.subplots_adjust(left=0.1, right=0.95, top=0.84, bottom=0.1)
    _plt.savefig('{0}/{1}.png'.format(path,star))
    _plt.savefig('{0}/{1}.eps'.format(path,star))
    _plt.close()
    
    bres = 20
    _plt.clf()
    leg = ()
    fig, ax = _plt.subplots()
    for f in filters:
        ind = _np.where(filt == f)
        x, y, yerr = _phc.bindata(MJD[ind], P[ind], sigP[ind], bres)
        leg += (f.upper()+' band',)
        ax.errorbar(x, y, yerr, marker='o', color=colors[filters.index(f)], fmt='-')
    ax.legend(leg,'upper left', fontsize='small')            
    ax.set_ylabel('Polarization (%) (binned)')
    ax.plot(ax.get_xlim(),[0,0],'k--')
    ax.set_ylim(ylim)
    ax.set_xlabel('MJD')
    xlim = ax.get_xlim()
    ticks = _phc.gentkdates(xlim[0], xlim[1], 3, 'm',\
    dtstart=_dt.datetime(2012,7,1).date())
    mjdticks = [_jdcal.gcal2jd(date.year,date.month,date.day)[1] for date in \
    ticks]
    ax2 = ax.twiny()
    ax2.set_xlabel('Civil date')
    ax2.set_xlim(xlim)
    ax2.set_xticks(mjdticks)
    ax2.set_xticklabels([date.strftime("%d %b %y") for date in ticks])
    _plt.setp( ax2.xaxis.get_majorticklabels(), rotation=45 )
    _plt.subplots_adjust(left=0.1, right=0.95, top=0.84, bottom=0.1)
    _plt.savefig('{0}/{1}_binned.png'.format(path,star))
    _plt.savefig('{0}/{1}_binned.eps'.format(path,star))
    _plt.close()       
    
    for f in filters:
        ind = _np.where(filt == f)
        avg,sigm = _phc.wg_avg_and_std(P[ind], sigP[ind])
        print('# Averaged {} band is {:.3f} +/- {:.3f} %'.format(f.upper(),avg,\
        sigm))
    return



#################################################
#################################################
#################################################
# Bednarski: This function generates the graphs for all filter "filt" logfiles found inside "objdir"
def grafall(objdir, filt, nstar=1, n=1, bestouts=[], shortmode=False):
    """
    Multiple plot modulations for the object inside 'objdir'
    in filter 'filt'.
    Return a list 'sortout' with the log files sorted by plot
    number. Allways sortout[0]=='' and sortout[i] is the plot #i.

    Optional:
        bestouts - list of outs which have smaller error to
                   highlight in figures
        shortmode - case True, plot only the groups 16001, 08001,
                    08009 (ver .1 and .2) and an eventual 7th group.
        n - is concerning to the n-esim 16-sequence to use.
            Exemple, if there are 40 WP positions, n=1 is to plot
            the graphs for positions 1-16; n=2 is for 17-32;
            n=3 is for 33-40.
    """

    # Receive a list of logfiles and return lists for each group of WP/version:
    # sublogs == list of sorted logfiles.
    # groups == list with informations about lamina groups.
    # ver == list with the reduction versions of files in sublogs.
    #
    # Ex, sublogs == [[*16001.1.log], [*16001.2.log], [*0800[1-9].1.log], [*0800[1-9].2.log]]
    #      groups == [   [16, 1, .1],    [16, 1, .2],         [8, 9, .1],         [8, 9, .2]]
    #         ver == [          [.1],         [best],           [.1 ...],          [.2, ...]]
    def combineout(logs, n):

        logsplit, groups, sublogs, ver = [], [], [[]], [[]]
        
        # Split the name files at last "_" character in each .log name. Working for _WP1111110110 ... files
        if len(logs) > 0:
            direc = _phc.trimpathname(logs[0])[0]
        else:
            print("ERROR: no log files found to plot.")
            return
        for log in [_phc.trimpathname(li)[1] for li in logs]:
            if log.find('_WP') == -1:
                indx = log.rfind('_')+1
            else:
                indx = log[:log.find('_WP')].rfind('_')+1
            combi = log[indx:indx+5]
            n1= int(combi[:2])
            n2= int(combi[2:])

            if n2  >= 16*(n-1)+1 and n2 <= 16*(n-1)+1 + (16-n1):
                logsplit += [[log[:indx], log[indx:]]]

        # Sort by lamina (high to low) -> version (including versions with _WP111100 ...)
        logsplit.sort(key=lambda x: [x[1][6:8],x[1][:]])
        logsplit.sort(key=lambda x: x[1][:1], reverse=True)

        j=0
        # Separate the lamina groups
        for i in range(len(logsplit)):
            sublogs[j] += [direc+logsplit[i][0]+logsplit[i][1]]
            if sublogs[j][-1][:-4]+'.out' in bestouts:
                ver[j] += ['best']
            else:
                ver[j] += [logsplit[i][1][5:-4]]
                      
            if (i != len(logsplit)-1 and (logsplit[i][1][0:2] != logsplit[i+1][1][0:2] or \
                        logsplit[i][1][6:8] != logsplit[i+1][1][6:8])) or i == len(logsplit)-1:
                groups += [[int(logsplit[i][1][0:2]), len(sublogs[j]), \
                                                        logsplit[i][1][5:-4]]]
                if i != len(logsplit)-1:
                    j+=1
                    sublogs[:] += [[]]
                    ver[:] += [[]]
            
        return groups, sublogs, ver
        

    # Generate background color for the graphs, depending the reduction version
    def gencolour(ver):

        if ver == 'best': bkg = '#d3dcf9'
        elif ver == '.1' or ver[:4] == '.1_WP': bkg = '#f0ffe0'
        elif ver == '.2' or ver[:4] == '.2_WP': bkg = '#f0fff0'
        else: bkg = '#f5f5f5'

        return bkg
        

    # Plot graphs for 'shortmode' and return 'sortout' list. Shortmode consists in the
    # processing of groups nn001.1, nn001.2 (nn is the number of WP), 08001.1, 08001.2,
    # 08009.1, 08009.2 and an eventual 7th element in 'bestouts' variable.
    # The variable returned is a list with the outfiles displayed, sorted in same order
    # that the showed.
    # The input variables are exactly the outputs "sublogs" and "groups" from combineout subroutine.
    def gengraphshort (sublogs, groups):

        # Variables below are to mark the positions in lists
        pos8ver1, pos8ver2, posnver1, posnver2 = -1,-1,-1,-1
        maxver1, maxver2 = 0,0

        # Find index for 16 and 08 groups positions
        for i in range(len(groups)):
            # Only shows the groups with 8 positions if there are 08001 to 08009 files
            if groups[i][0] == 8 and groups[i][1] == 9 and groups[i][2] == '.1':
                pos8ver1 = i
                if maxver1 < 8:
                    maxver1 = 8
            elif groups[i][0] > maxver1 and groups[i][2] == '.1':
                maxver1 = groups[i][0]
                posnver1 = i
            if groups[i][0] == 8 and groups[i][1] == 9 and groups[i][2] == '.2':
                pos8ver2 = i
                if maxver2 < 8:
                    maxver2 = 8
            elif groups[i][0] > maxver2 and groups[i][2] == '.2':
                maxver2 = groups[i][0]
                posnver2 = i

        # Set the logfiles in a first time
        tlogs = ['']*6
        tver = ['.1','.2']*3
        if posnver1 != -1:
            tlogs[0] = sublogs[posnver1][0]
        if posnver2 != -1:
            tlogs[1] = sublogs[posnver2][0]
        if pos8ver1 not in (-1,posnver1):
            tlogs[2], tlogs[4] = sublogs[pos8ver1][0], sublogs[pos8ver1][8]
        if pos8ver2 not in (-1,posnver2):
            tlogs[3], tlogs[5] = sublogs[pos8ver2][0], sublogs[pos8ver2][8]

        # Set the bestout
        if bestouts != []:
            if len(bestouts) > 1:
                print('# WARNING: grafall: more than one value of best .out passed as' + \
                ' parameter in short mode. Only using the first one...\n')
            if bestouts[0][:-4]+'.log' not in tlogs:
                tlogs += [bestouts[0][:-4]+'.log']
                tver += ['best']
            else:
                tver[tlogs.index(bestouts[0][:-4]+'.log')] = 'best'

        # Set the logfiles once, erasing the void components
        if (posnver1 == -1 and pos8ver1 == -1) or (posnver2 == -1 and pos8ver2 == -1):
            logs = [tlogs[i] for i in range(len(tlogs)) if tlogs[i] != '']
            ver = [tver[i] for i in range(len(tlogs)) if tlogs[i] != '']
            mode='lin'
        else:
            mode='col'
            logs = []
            ver = []
            if posnver1 != -1 or posnver2 != -1:
                logs += tlogs[:2]
                ver += tver[:2]
            if pos8ver1 != -1 or pos8ver2 != -1:
                logs += tlogs[2:6]
                ver += tver[2:6]
            if len(tlogs) == 7:
                logs += tlogs[6:]
                ver += tver[6:]

        # Run gengraphl4!
        gengraphl4(logs,ver,1,align=mode)

        # Return the 'sortlog' file (sorted outfiles with extension '.log')
        return [''] + [log[:-4]+'.out' for log in logs if log != '']
        
    
    # Generate graphs for the cases nlog <= 8
    # align: align by lines or columns? 1234//5678 ('lin') or 1357//2468 ('col')?
    def gengraphl4(logs, ver, count, align='lin'):

        if align=='lin':
            if len(logs) < 4:
                nlin, ncol = 1, len(logs)
            else:
                nlin, ncol = 2, 4
        elif align=='col':
            if len(logs) == 1:
                nlin, ncol = 1, 1
            else:
                nlin, ncol = 2, len(logs)/2
                if len(logs)%2 != 0:
                    ncol += 1
        else:
            print('# ERROR: align mode {0} is not valid in grafall! Graphs not displayed!'.format(align))
            return            
        if ncol > 4:
            print('# ERROR: {0} figure(s) was(were) not displayed by grafall'.format(len(logs)-8))
            return

        if   ncol == 1: linit=0.15
        elif ncol == 2: linit=0.08
        else:           linit=0.05

        if   nlin == 1: binit=0.12; tinit=0.88
        elif nlin == 2: binit=0.19; tinit=0.94

        fig = _plt.figure(figsize=(4*ncol,3.4*nlin))

        # Estabelece os grids e cria todos os eixos
        grids = [ _plt.GridSpec(2*nlin, ncol, hspace=0, wspace=0.35, \
                        top=tinit, bottom=binit, left=linit, right=0.95) ]
        if nlin == 2:
            grids += [ _plt.GridSpec(2*nlin, ncol, hspace=0, wspace=0.35, \
                        top=0.81, bottom=0.06, left=linit, right=0.95) ]
        ax = []
        if align == 'lin':
            for j in range(ncol):
                if logs[j] != '':
                    ax += [ fig.add_subplot(grids[0][0,j]),\
                            fig.add_subplot(grids[0][1,j]) ]
                else:
                    ax += [None,None]
            for j in range(len(logs)-ncol):
                if logs[j+ncol] != '':
                    ax += [ fig.add_subplot(grids[1][2,j]),\
                            fig.add_subplot(grids[1][3,j]) ]
                else:
                    ax += [None,None]
        elif align == 'col':
            for j in range(ncol):
                if logs[2*j] != '':
                    ax += [ fig.add_subplot(grids[0][0,j]),\
                            fig.add_subplot(grids[0][1,j]) ]
                else:
                    ax += [None,None]
                if 2*j+1 < len(logs) and logs[2*j+1] != '':
                    ax += [ fig.add_subplot(grids[1][2,j]),\
                            fig.add_subplot(grids[1][3,j]) ]
                else:
                    ax += [None,None]

        k=0
        for j in range(len(logs)):
            # Case of even logs, breaks after the last one
#            print logs[j]
            if len(logs) <= j:
                break
            elif logs[j] != '':
                grafpol(logs[j], nstar, fig, ax[2*j], ax[2*j+1])
                ax[2*j].text(0.85, 0.85, '#{0:<2d}'.format(count+k), \
                    horizontalalignment='left', verticalalignment='center', style='italic', \
                    transform=ax[2*j].transAxes, fontsize=20, color='red')
                ax[2*j].set_axis_bgcolor(gencolour(ver[j]))
                ax[2*j+1].set_axis_bgcolor(gencolour(ver[j]))
                k += 1
            
        _plt.show(block=False)


    # Generate graphs for the cases nlog > 8
    def gengraphm4(logs, ver, count):

        nwin = len(logs)/12 + 1

        for i in range(nwin):

            # set nlin/ncol values
            nlog = len(logs)-i*12
            ncol = 4
            if   nlog <= 4:  nlin=1; ncol=nlog
            elif nlog <= 8:  nlin=2
            elif nlog <= 12: nlin=3
            else:            nlin=3; nlog=12

            # set left/right parametes
            if   nlog == 1: linit=0.15
            elif nlog == 2: linit=0.08
            else:           linit=0.05

            # set top/bottom parameters
            if   nlog <= 4: delt=0.10; tinit=0.86; binit=0.00
            elif nlog <= 8: delt=0.13; tinit=0.90; binit=0.00
            else:           delt=0.10; tinit=0.95; binit=0.05
            
            fig = _plt.figure(figsize=(4*ncol,3*nlin))

            # Creates the axes of first row
            grids = [ _plt.GridSpec(2*nlin, ncol, hspace=0, wspace=0.35, \
                            top=tinit, bottom=binit+2*delt, left=linit, right=0.95) ]
            ax = []
            for j in range(0,ncol):
                ax += [ fig.add_subplot(grids[0][0,j]),\
                        fig.add_subplot(grids[0][1,j]) ]

            # Creates the axes of second row
            if nlin > 1:
                grids += [ _plt.GridSpec(2*nlin, ncol, hspace=0, wspace=0.35, \
                            top=tinit-delt, bottom=binit+delt, left=linit, right=0.95) ]
                for j in range(0,ncol):
                    if j+ncol >= nlog:
                        break
                    ax += [ fig.add_subplot(grids[1][2,j]),\
                            fig.add_subplot(grids[1][3,j]) ]

            # Creates the axes of third row
            if nlin > 2:
                grids += [ _plt.GridSpec(2*nlin, ncol, hspace=0, wspace=0.35, \
                            top=tinit-2*delt, bottom=binit, left=linit, right=0.95) ]
                for j in range(0,ncol):
                    if j+2*ncol >= nlog:
                        break
                    ax += [ fig.add_subplot(grids[2][4,j]),\
                            fig.add_subplot(grids[2][5,j]) ]

            # Generates the plots on the axes
            for j in range(0,nlog):
                grafpol(logs[j+12*i], nstar, fig, ax[2*j], ax[2*j+1])
                ax[2*j].set_axis_bgcolor(gencolour(ver[j]))
                ax[2*j+1].set_axis_bgcolor(gencolour(ver[j]))
                ax[2*j].text(0.85, 0.85, '#{0:<2d}'.format(count+j+i*12), \
                        horizontalalignment='left', verticalalignment='center', \
                        style='italic', transform=ax[2*j].transAxes, fontsize=20, color='red')

            _plt.show(block=False)


    logs = _glob('{0}/*_{1}_*.log'.format(objdir, filt))
    if logs == []:
        print('# ERROR: log files not found to plot. May the file names \
              {0}/*_{1}_*.log are wrong!'.format(objdir, filt))
        return 1
    gps, sublogs, ver = combineout(logs, n)

    # 1) Case short mode
    if shortmode:
        sortout = gengraphshort(sublogs, gps)

    # 2) Case long mode
    else:
        nlog = sum([len(subb) for subb in sublogs])
        # If a few logfiles, tries to use only one window
        if nlog <= 8:
            test=True
            # Test if all groups have two reduction versions
            if len(gps)%2 == 0:
                for i in range(0,len(gps),2):
                    # [:2] and not [:1] because gps is a list of lists
                    if gps[i][:2] != gps[i+1][:2]:
                        test=False
                        break
            else:
                test=False
            if test:
                tver = []
                tlogs = []
                for i in range(len(sublogs)):
                    tlogs += sublogs[i]
                    tver += ver[i]
                gengraphl4(tlogs, tver, 1)
        # Otherwise, loop on lamina groups
        else:
            i = 0
            count = 1       # variable to print the graph number
            while i < len(gps):
                nout = 0
                iver = []
                ilogs = []
    #            print i
                for j in range(i,len(gps)):
                    if nout <= 8 and gps[j][:2] == gps[i][:2]:
                        nout += gps[j][1]
                        iver += ver[j]
                        ilogs += sublogs[j]
                    else:
                        break
                # if isn't last gps element, nout<=8 and number of versions is 2, 4, 6, etc
                if j != len(gps)-1 and nout <= 8 and (j-i)%2 == 0:
                    gengraphl4(ilogs,iver,count)
                    count += nout
                    i = j
    #                print('entrou 1')
                else:
                    gengraphm4(sublogs[i],ver[i],count)
                    count += len(sublogs[i])
    #                print('entrou 2')
                    i += 1

        sortout = ['']
        for ilogs in sublogs:
            sortout += [log[:-4]+'.out' for log in ilogs]

    # returns the sorted logs/outs, with '.log' changed to '.out'
    return sortout



#################################################
#################################################
#################################################
def grafpol(filename, nstar=1, fig=None, ax1=None, ax2=None, save=False, extens='png'):
    """
    Program to plot the best adjust and its residuals of IRAF reduction.
    'filename' is the path to the .log output file from reduction.
    nstar is the star number inside out/log files to be plotted.

    NEW: Working for *_WP1110....log files!
    NEW (2): Working for logfiles with more than a single star!

    Two working modes:

    1) If only filename is given: displays the graph or, case save=True,
       save it into the logfile directory.
       'extens' parameter changes the extension.
    2) If filename, one figure and two axes are given: changes these axes,
       adding plots in them. Doesn't display, either save at the end,
       only adds a plot to the axes. Usefull for subplots in same figure.
  
    Authors: Moser and Bednarski
    Current version: May, 2015
    """
    
    def readlog(filename):

        try:
            # CAUTION! BLANK LINES WILL BE SKIPPED!
            file0 = _np.loadtxt(filename, dtype=str, delimiter='\n', comments=None)
        except:
            print('# ERROR: File {0} not found!\n'.format(filename))
            raise SystemExit(1)

        [lixo,lixo,lixo,lixo,lixo,lixo,lixo,lixo,MJD,lixo] = \
                readoutMJD(filename.replace('.log','.out'), nstar=nstar)
        MJD = float(MJD)


        # npts: number of WP valid to be plotted.
        # totpts: used for the total number of WP (for the *.2_WP11110...1.out)
        nstars = int(file0[6].split()[-1])
        npts = int(file0[9].split()[-1])    # Blank lines were SKIPPED!
        tnpts = int(file0[8].split()[-1])
        delta = float(file0[14].split()[-1])
        sigma = 1.
        isinstar=False

        if nstars < nstar:
            print('# ERROR: File {0} has {1} stars (you have selected star #{2}).'.\
                                                        format(filename, nstars, nstar))
            raise SystemExit(1)
            
        # Bednarski: corrected (25 -> 19, because the blank lines had been ignorated by
        #            np.loadtext function)
        for i in range(19, len(file0)):
            if 'STAR # {0:d}'.format(nstar) in file0[i]:
                isinstar=True
            elif 'STAR # {0:d}'.format(nstar+1) in file0[i]:
                break
            if isinstar and 'APERTURE' in file0[i]:
                sig = float(file0[i+2].split()[2])
                if sig < sigma:
                    sigma = sig

                    fator = thtFactor(MJD)
                    thet = fator * float(file0[i+2].split()[4])
                    while thet >= 180:
                        thet -= 180
                    while thet < 0:
                        thet += 180
                        
                    '''
                    # Bed: Os Q e U sao os abaixo, conforme copiei da rotina graf.cl
                    if float(file0[i+2].split()[4]) < 0:
                        thet = - float(file0[i+2].split()[4])
                    else:
                        thet = 180. - float(file0[i+2].split()[4])
                    '''

                    # Recalculating the new QU parameters
                    Q = float(file0[i+2].split()[3])*_np.cos(2.*thet*_np.pi/180.)
                    U = float(file0[i+2].split()[3])*_np.sin(2.*thet*_np.pi/180.)
#                    print Q, U, thet, float(file0[i+2].split()[3])
                    n = npts/4 
                    if npts%4 != 0:
                        n = n+1
                    P_pts = []
                    for j in range(n):
                        P_pts += file0[i+4+j].split()
                    # I think the P values are in reverse order inside Pereyra's .log files:
                    # Uncomment the two lines below if you want to show the ascending x-axes
                    # (and not descending x-axes) and comment the next two lines.
#                    P_pts = _np.array(P_pts, dtype=float)[::-1]
#                    th_pts = fator*(22.5*_np.arange(1,tnpts+1)+delta/2.)
                    P_pts = _np.array(P_pts, dtype=float)
                    th_pts = -fator*(22.5*_np.arange(tnpts)-delta/2.)
                    j = filename.find('.')
                    delta2 = int(filename[-2+j:j])-1
                    # Bed: Funcionando para nlam >= 10  para impressão correta
#                    str_pts = map(str, _np.arange(1,tnpts+1)+delta2)[::-1]
                    str_pts = map(str, _np.arange(1,tnpts+1)+delta2)

                    # Case _WP11110...1.log file
                    if npts != tnpts:
                        refs = file0[9].split()[3:-2]
                        rm = [j for j in range(tnpts) if refs[j] == '0']
                        th_pts = _np.delete(th_pts, rm)
                        str_pts = _np.delete(str_pts, rm)

        if sigma == 1.:
            print('# ERROR reading the file %s !' % filename)
            Q = U = 0
            P_pts = th_pts = _np.arange(1)
            str_pts = ['0','0']

        return(Q, U, sigma, P_pts, th_pts, str_pts, nstars, fator)


    def plotlog(ax1, ax2, Q,U,sigma,P_pts,th_pts,str_pts,filename,fator):

        # extract the group number
        WPpos = filename.find('_WP')
        if WPpos == -1:
            suff = filename[filename.rfind('_')+1:]
        else:
            suff = filename[filename.rfind('_', 0, WPpos)+1:]

        ax1.set_title(r'Q={0:.3f}, U={1:.3f}, $\sigma$={2:.3f}'.format(Q*100,U*100,sigma*100),
                       fontsize=14, verticalalignment='bottom')
        ax1.text(0.98, 0.01, '{0}'.format(suff), horizontalalignment='right', \
                 verticalalignment='bottom', transform=ax1.transAxes, fontsize=9)
        ax1.set_ylabel('p (%)', size=9)
        
        ysigma = _np.zeros(len(th_pts))+sigma
        ax1.errorbar(th_pts,P_pts*100,yerr=ysigma*100)
    
        th_det = _np.linspace(th_pts[0]*.98,th_pts[-1]*1.02,100)
        P_det = Q*_np.cos(4*th_det*_np.pi/180)+U*_np.sin(4*th_det*_np.pi/180)

        ax1.plot(th_det, P_det*100)
        ax1.plot([th_det[0],th_det[-1]], [0,0], 'k--')    
        ax1.set_xlim([th_pts[0]+fator*4,th_pts[-1]*1.02-fator*1.5])
#        ax1.set_xlim([th_pts[0]-4,th_pts[-1]*1.02+3][::-1])
#        ax1.set_ylim([min(P_pts*100)*1.1, max(P_pts*100)*1.1])

        ax2.set_xlabel('WP position', size=9)
        ax2.set_ylabel('Residuals', size=9)
        ax2.set_xlim([th_pts[0]+fator*4,th_pts[-1]*1.02-fator*1.5])
#        ax2.set_xlim([th_pts[0]-4,th_pts[-1]*1.02+3][::-1])
        ax2.plot([th_det[0],th_det[-1]], [0,0], 'k--')

        _plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.set_yticks(ax1.get_yticks()[1:])
        transOffset = _offset_copy(ax2.transData, fig=fig, x=0.00, y=0.10, units='inches')
        
        P_fit = Q*_np.cos(4*th_pts*_np.pi/180)+U*_np.sin(4*th_pts*_np.pi/180)
        
        # Bed: Agora plota os residuos relativos (residuos divididos por sigma)
        ax2.errorbar(th_pts, (P_pts-P_fit)/sigma, yerr=1.)
    
        for i in range(len(th_pts)):
            ax2.text(th_pts[i], (P_pts-P_fit)[i]/sigma, str_pts[i], transform=transOffset)  

        if int(ax2.get_yticks()[0])-int(ax2.get_yticks()[-1]) > 5:
            passo = 1
        else:
            passo = 2

        ax2.set_yticks(range(int(ax2.get_yticks()[0]),int(ax2.get_yticks()[-1]+1), passo))
        ax1.set_xticklabels(ax1.get_xticks(), size=7)
        ax1.set_yticklabels(ax1.get_yticks(), size=7)
        ax2.set_xticklabels([int(ax2.get_xticks()[i]) for i in range(len(ax2.get_xticks()))], size=7)
        ax2.set_yticklabels(ax2.get_yticks(), size=7)

        return



    if fig == None or ax1 == None or ax2 == None:
        _plt.close('all')
        fig = _plt.figure(1)
        ax1 = _plt.subplot(2, 1, 1)
        ax2 = _plt.subplot(2, 1, 2, sharex=ax1)
        _plt.subplots_adjust(hspace = 0)
        Q, U, sigma, P_pts, th_pts, str_pts, nstars, fator = readlog(filename)
        plotlog(ax1,ax2, Q,U,sigma,P_pts,th_pts,str_pts,filename,fator)
        if save:
            if nstars == 1:
                _plt.savefig(filename.replace('.log','.'+extens))
            else:
                _plt.savefig(filename.replace('.log','_star{0}.{1}'.format(nstar,extens)))
        else:
            _plt.show()
    else:
        Q, U, sigma, P_pts, th_pts, str_pts, nstars, fator = readlog(filename)
        plotlog(ax1,ax2, Q,U,sigma,P_pts,th_pts,str_pts,filename,fator)
            
    return



#################################################
#################################################
#################################################
def verStdPol(std, filt, p, sig):
    """
    Calculate z test for standard 'std', filter 'filt', comparing
    the observed polarization 'p' +- 'sig' and the published value.

    Return z = abs(ppub-p)/sqrt(sigpub^2+sig^2) or -1 if there is
    no such object or filter.
    """
    lstds = _np.loadtxt('{0}/refs/pol_padroes.txt'.format(_hdtpath()), dtype=str, usecols=range(0,22))

    # Get P_pub value
    i = stdchk(std)[1]
    if i == -1:
        return -1
    else:
        j = filters.index(filt)+11 # +11 devido as colunas a serem puladas pra
                                   # chegar a coluna das polarizacoes
        ppub = float(lstds[i,j])
        sppub = float(lstds[i,j+5])

    if ppub==0. or (sppub==0. and sig==0.):
        return -1
#    if ztest > 2.5
#       loglines += ('# CAUTION! Standard {0}, {1}, has polarization only compatible '+\
#                                  'within {2:.1f} sigma.\n').format(std, f, ztest)

    return abs(ppub-p)/_np.sqrt(sppub**2+sig**2)



#################################################
#################################################
#################################################
def readTests(tests, tags=None, flag=None):
    """
    Read boolean list 'tests' concerning to dictests
    dictionary and return a string list with tags and
    the flag value ('OK','E','W')

    'tags' and 'flag' are optional and are concerning to
    another tags already assigned and the current flag. The flag
    returned is the worst flag found between them
    (e.g., if input 'flag' is 'W' and tests results on flag
    'OK', return flag 'W'; if tests results in flag 'E',
    return 'E'). Also, if they are given as input, return
    'tags'+tags concerning to 'tests' list.
    """

    tagstr = ''

    # Generate a string for such tests tags
    for i in dictests.keys():
        if tests[i]:
            if tagstr != '':
                tagstr += ','
            tagstr += dictests[i][1]

    # Generate a string for such tags
    if tags != None:
        for i in dictags.keys():
            if tags[i]:
                if tagstr != '':
                    tagstr += ','
                tagstr += dictags[i][1]

    if flag == None:
        flag = 'OK'

    # Get the worst case for the flag
    if 'E' in [dictests[j][2] for j in range(len(tests)) if tests[j]]+[flag]:
        flag2 = 'E'
    elif 'W' in [dictests[j][2] for j in range(len(tests)) if tests[j]]+[flag]:
        flag2 = 'W'
    else:
        flag2 = 'OK'

    if tagstr == '':
        tagstr = '---'
    
    return tagstr, flag2

    

    
#################################################
#################################################
#################################################
# Bednarski: I added delta variable to (calc-calcst) tolerance
def chkStdLog(f, calc, path=None, delta=3.5, verbose=True):
    """
    Verify if there are standards for filter `f` and
    calcite `calc` inside path/std.dat. Return True if
    successful, unless the standard has been reduced
    and marked with `E` flag.

    delta is the allowed variation for the angles between the two
    beams for one same calcite.
    """

    loglines = ''
    if path == None or path == '.':
        path = _os.getcwd()

    # Read `obj.dat` and `std.dat`. If there are errors, assigns [''] to get inside
    # ifs below and print error messages
    try:
        std = _np.loadtxt('{0}/std.dat'.format(path), dtype=str)
    except:
        std = _np.array([], dtype=str)

    # Verify if std.dat has more than one line. Caso no, do reshape (transform
    # list type [] in [[]] for further compatibility)
    if _np.size(std) != 0:
        # If std is of type [] -- one line with 9 elements:
        if type(std[0]) != _np.ndarray and _np.size(std) == 9:
            std = std.reshape(-1,9)
        elif (type(std[0]) == _np.ndarray and _np.size(std[0]) != 9) \
                                or (type(std[0]) != _np.ndarray and _np.size(std) != 8):
            # Save loglines only if this function was called by genLog
            if _getouterframes(_currentframe(), 2)[1][3] == 'genLog':
                writeLog(path, '# ERROR: polt.chkStdLog() not runned! Incompatible number '+ \
                                                            'of columns in `std.dat`.\n')
            else:
                print('# ERROR: Incompatible number of columns in `std.dat`.\n')

    foundstd = False
    for stdi in std:
    # Skip if stdi is not to use ('Error' flag)
        if stdi[7] == 'E':
            continue
        fst = stdi[3]
        calcst = float(stdi[4])
        if f == fst and abs(calc-calcst) < delta:
            foundstd = True
            break
    if not foundstd and verbose:
        print(('# WARNING: Standard star not found for filt. {0} and '+\
            'calc. {1:.1f}\n').format(f, calc))

    return foundstd



#################################################
#################################################
#################################################
def writeLog(path, strin):
    """
    Append 'strin' string into 'path'/polt.log file

    """
    
    f0 = open('{0}/polt.log'.format(path), 'a')
    f0.writelines(strin)
    f0.close()

    return



#################################################
#################################################
#################################################
def genLog(path, subdirs, tgts, fileout, sigtol=lambda sigm: 1.4*sigm, \
                    autochoose=False, delta=3.5):
    """
    Generate the .dat file with data of objects 'tgts[:]' inside
    'path'/'subdirs[:]' directories

    Save the results in 'path'/'fileout'
    Usable to generate target and standard lists (obj.dat and std.dat)

    delta: tolerance for the angle between the two beams of calcite.
           If abs(angle1 - angle2) < delta, both observations 1 and 2
           are assigned as the same calcite.
    sigtol: tolerance to use the outfiles with all WP instead the
            out with best error. Must be a 'function' that receives
            a pol sigma value (in decimal system and NOT in per cent,
            i.e., value from 0. to 1., where 1 is a 100% polarized
            source) and return the maximum sigma for which to ignore
            the best out. Its format must be a 'python's lambda function'!
            The default values is sigtol=lambda sigm: 1.4*sigm, while
            the old value was sigtol=lambda sigm: 1.1*sigm + 0.00005. If
            you want to take just the groups with all WP, and none other,
            you can specify sigtol=lambda sigm: 1000.*sigm, for example.
    autochoose: choose best outfiles automatically, without
                interaction?
    """

    if fileout.split('.')[0] == 'std':
        typ = 'standards'
    elif fileout.split('.')[0] == 'obj':
        typ = 'targets'
    else:
        typ = fileout

#    if mode not in ('std','obj'):
#        print('\n# ERROR: mode \'{0}\' not valid (it\'s only valid \'std\' and \'obj\')'.format(mode))
#        return 1

    if len(tgts) != len(subdirs):
        print('\n# ERROR: polt.genLog() NOT RUNNED for {0} (len(tgts) != len(subdirs))'.format(typ))
        writeLog(path, '# ERROR: polt.genLog() NOT RUNNED for {0}! (len(tgts) != len(subdirs))\n'.format(typ))
        return 1

    continuerun = False
    # Checking if there exists a previous run and if it has generated unless one line.
    if _os.path.exists('{0}/{1}.tmp'.format(path,fileout)) and len(_np.loadtxt('{0}/{1}.tmp'.format(path,fileout), dtype=str)) != 0:
        opt = ''
        while opt not in ('y','Y','n','N'):
            opt = raw_input(('There exists one file concerning to a uncompleted previous run for {0}. ' +\
                            'Do you want to continue where it was stopped? (y/n): ').format(typ))
            if opt in ('y','Y'):
                continuerun = True

    if not continuerun:
        f0 = open('{0}/{1}.tmp'.format(path,fileout), 'w')
        f0.writelines('{:12s} {:>7s} {:>10s} {:4s} {:>5s} {:<s} {:>4s} {:>5s}  {:<s}\n'.format('#MJD','ccd',\
                    'target','filt','calc',':::outfile:::','star','flag','tags'))
        f0.close()
    # Case continuing a previous run, identify the stars already runned
    else:
        ftemp = _np.loadtxt('{0}/{1}.tmp'.format(path,fileout), dtype=str)
        odone=[]  # odone and fdone is lists that contains subdirectories, star number and the
        fdone=[]                                    # filters already done by the previous run
        # If there is just one line, transform np array type [] for [[]]
        if len(ftemp) > 0 and len(ftemp[-1]) != 9:
            ftemp = ftemp.reshape(-1,9)
        for line in ftemp:
            # [3:] because the firsts characters are ':::'
            objct = [line[5].split('/')[0][3:], line[6]]
            if objct in odone:
                indx = odone.index(objct)
                fdone[indx] += [line[3]]
            else:
                odone += [objct]
                fdone += [[line[3]]]
#        print odone
#        print fdone
        
    # Loop on list of objects
    for i in range(len(tgts)):

        obj = tgts[i]
        objdir = subdirs[i]
        if obj == '':
            continue

        # Loop on filters
        for f in filters:

            nstars = countStars('{0}/{1}'.format(path,objdir), f)

            # Check if there exist fits files for object/filter, but not .out files (target not reduced)
            if nstars == 0 and (len(_glob('{0}/{1}/*_{2}_*.fits'.format(path,objdir,f))) > 0 \
                            or len(_glob('{0}/{1}/{2}/p??0'.format(path,objdir,f))) > 0):
                print(('\n# ERROR: {0}_{1}: Fits files found, but the object was not reduced! ' +\
                        'Reduce and run again...\n\n - HINT: if these fits files compose some ' +\
                        'non-valid serie but need be kept in, move them for a subdir {2}/tmp, ' +\
                        'and hence, the path will not be sweept by routine.\n').format(objdir,f,objdir))
                raise SystemExit(1)
            # Check if there exist some .out file for such object/filter, but not the fits files
            elif nstars != 0 and (len(_glob('{0}/{1}/*_{2}_*.fits'.format(path,objdir,f))) == 0 \
                                    and len(_glob('{0}/{1}/{2}/p??0'.format(path,objdir,f))) == 0):
                print(('\n# ERROR: {0}_{1}: Fits files not found, but were found *_{2}_* files. ' +\
                        'It can be by three reasons:\n'+\
                        '  1) Fits files missing (in this case, search by them and add in such directory);\n' +\
                        '  2) The found *_{3}_* files can be \'spurious files\' (in this case, delete them);\n' +\
                        '  3) The preffix of *_{4}_*.fits files can be different of the rest of *_{5}_* ' +\
                        'files (in this case, rename them).\n').format(objdir,f,f,f,f,f))
                raise SystemExit(1)

            # Working for more than a single star inside .out files
            for nstar in range(1,nstars+1):

                loglines = ''
                # Skip if the object/filter was done already in a previous run
                if continuerun and [objdir,str(nstar)] in odone and \
                                                f in fdone[odone.index([objdir,str(nstar)])]:
                    continue
                elif nstars == 1:
                    obj = tgts[i]   # It is needed this line here again
                else:
                    while True:
                        obj = raw_input(('Type a name for star #{0:d} (of {1:d}) ' +\
                               'inside {2} dir, filter {3}: '). format(nstar, nstars, objdir, f))
                        if obj != '' and ' ' not in obj and '#' not in obj and ':' not in obj:
                            loglines += '# WARNING: {0}_{1}: There are more than one star inside .out files.'\
                                     .format(objdir, f)+' Star #{0:d} was named as {1}.\n'.format(nstar, obj)
                            break

                if autochoose:
                    outs = chooseout('{0}/{1}'.format(path,objdir), obj, f, nstar=nstar, sigtol=sigtol)
                    tags=None
                    flag=None
                else:
                    outs, tags, flag = queryout('{0}/{1}'.format(path,objdir), obj, f, nstar=nstar, sigtol=sigtol)
                    print('')
                        
                # Loop on outfiles
                lines = ''
                for j in range(len(outs)):
                    if outs[j] != '':
                        [Q,U,sig,P,th,sigT,ap,star,MJD,calc] = readoutMJD(outs[j], nstar=nstar)
                        tests, logs = verout(outs[j], obj, f, nstar=nstar, verbose=False, delta=delta)
                        loglines += logs
                        if tags!=None and flag!=None:
                            tagstr, flagout = readTests(tests, tags=tags[j], flag=flag[j])
                        else:
                            tagstr, flagout = readTests(tests, tags=None, flag=None)
                        
                        # It is needed to open and close in each object
                        lines += ('{0:12.6f} {1:>7s} {2:>10s} {3:>4s} {4:>5.1f} {5:<s} {6:>4d} ' \
                                    '{7:>5s}  {8:<s}\n').format(MJD, ccd, obj, f, float(calc), \
                                    ':::'+_os.path.relpath(outs[j], path)+':::', nstar, flagout, tagstr)

                # Write lines after process all outfiles for object in one filter
                if lines != '':
                    writeLog(path, loglines)
                    f0 = open('{0}/{1}.tmp'.format(path,fileout), 'a')
                    f0.writelines(lines)
                    f0.close()


    # Read fileout+'.tmp', realign columns to fileout and delete fileout+'.tmp'
    fin = open('{0}/{1}.tmp'.format(path,fileout), 'r')
    lines = [li.split(':::') for li in fin.readlines()]
    fin.close()
    if len(lines) == 1:
        writeLog(path, '# WARNING: No valid {0} were found by polt.genLog().\n'.format(typ))
    else:
        maxsize = max([len(li[1]) for li in lines])
        linesout = []
        for i in range(len(lines)):
            if lines[i] in ('','\n'):
                continue        
            linesout += [lines[i][0]+lines[i][1].rjust(maxsize+2)+lines[i][2]]
        fout = open('{0}/{1}'.format(path,fileout), 'w')
        fout.writelines(linesout)
        fout.close()
    try:
        _os.unlink('{0}/{1}.tmp'.format(path,fileout))
    except:
        pass


    return 0




#################################################
#################################################
#################################################
def genAllLog(path=None, sigtol=lambda sigm: 1.4*sigm, autochoose=False, delta=3.5):
    """
    Generate the std.dat/obj.dat for one reduced night
    
    path: path of night
    delta: tolerance for the angle between the two beams of calcite.
           If abs(angle1 - angle2) < delta, both observations 1 and 2
           are assigned as the same calcite.
    sigtol: tolerance to use the outfiles with all WP instead the
            out with best error. Must be a 'function' that receives
            a pol sigma value (in decimal system and NOT in per cent,
            i.e., value from 0. to 1., where 1 is a 100% polarized
            source) and return the maximum sigma for which to ignore
            the best out. Its format must be a 'python's lambda function'!
            The default values is sigtol=lambda sigm: 1.4*sigm, while
            the old value was sigtol=lambda sigm: 1.1*sigm + 0.00005. If
            you want to take just the groups with all WP, and none other,
            you can specify sigtol=lambda sigm: 1000.*sigm, for example.
    autochoose: choose best outfiles automatically, without
                interaction?
    """

    if path == None or path == '.':
        path = _os.getcwd()

    # Verifies if std.dat and obj.dat files exist. Case True, queries to delete
    # Doesn't delete polt.log because the aiming is keep every information about the previous run
    if _os.path.exists('{0}/std.dat'.format(path)):
        if _os.path.exists('{0}/obj.dat'.format(path)):
            while True:
                verif = raw_input('Caution: obj.dat and std.dat already exists. Are you sure to overwrite it/them (y/n): ')
                if verif in ('n','N'):
                    print('Aborted!')
                    return
                elif verif in ('y','Y'):
                    break
            for arq in (path+'/obj.dat', path+'/std.dat'): #, path+'/polt.log'):
                try:
                    _os.unlink(arq)
                except:
                    pass
        elif not _os.path.exists('{0}/obj.dat'.format(path)) and _os.path.exists('{0}/obj.dat.tmp'.format(path)):
            print('# WARNING: keeping file std.dat and processing only obj.dat...\n')

    # Generates lists
    try:
        ltgts = _np.loadtxt('{0}/refs/pol_alvos.txt'.format(_hdtpath()), dtype=str)
        if _os.path.exists('{0}/refs/pol_hip.txt'.format(_hdtpath())):
            try:
                ltgts = _np.concatenate((ltgts,_np.loadtxt('{0}/refs/pol_hip.txt'.\
                                                    format(_hdtpath()), dtype=str)))
            except:
                pass
        if _os.path.exists('{0}/refs/pol_unpol.txt'.format(_hdtpath())):
            try:
                ltgts = _np.concatenate((ltgts,_np.loadtxt('{0}/refs/pol_unpol.txt'.\
                                                    format(_hdtpath()), dtype=str)))
            except:
                pass
        lstds = _np.loadtxt('{0}/refs/pol_padroes.txt'.format(_hdtpath()), \
                                                           dtype=str, usecols=[0])
    except:
        print('# ERROR: Can\'t read files pyhdust/refs/pol_alvos.txt and/or pyhdust/refs/pol_padroes.txt.\n')
        raise SystemExit(1)

    subdirs = [fld for fld in _os.listdir('{0}'.format(path)) if \
            _os.path.isdir(_os.path.join('{0}'.format(path), fld))]
    tgts=[]    # variable for real target names, not the subdirectory names
    stds=[]    # variable for real standard names, not the subdirectory names
    lines = ''

    # set ccd name from first .fits file
    try:
        setCCD(_glob('{0}/{1}/*.fits'.format(path, next(k for j,k in enumerate(subdirs)\
                                                                            if k!='')))[0])
    except:
        setCCD('')     # set manually inside the function

    # Verifies if object is a standard star or target
    # (Works on directories with suffix also (like 'dsco_a0'))
    for obj in [elem.split('_')[0] for elem in subdirs]:
        obj_curr = obj
        while obj not in _np.hstack((ltgts,lstds,_np.array(['calib','flat','dark','bias']))):
            if obj_curr == obj:
                print('\nObject {0} is not a known target or standard!!'.format(obj_curr))
            else:
                print('\nObject {0} (and {1}) is not a known target or standard!!'.format(obj_curr, obj))
            obj = raw_input('Type the common name (you can add a new target inside pyhdust/ref/pol_*.txt files, but be careful!): ')
        if obj in lstds:
            tgts.append('')
            # Only assigns standard's name if there is no std.dat file. Otherwise,
            # it's because the actual run will process only the targets, not standards.
            if not _os.path.exists('{0}/std.dat'.format(path)):
                stds.append(obj)
        elif obj in ltgts:
            tgts.append(obj)
            stds.append('')
        elif obj in ('calib','flat','dark','bias'):
            tgts.append('')
            stds.append('')

    print('')
    writeLog(path, '#### BEGIN\n')
    if not _os.path.exists('{0}/std.dat'.format(path)):
        genLog(path, subdirs, stds, fileout='std.dat', delta=delta, sigtol=sigtol, autochoose=autochoose)
    genLog(path, subdirs, tgts, fileout='obj.dat', delta=delta, sigtol=sigtol, autochoose=autochoose)

    # Write user name and date+time
    username = _pwd.getpwuid(_os.getuid())[4]
    if username.find != -1:
        username = username[:username.find(',')]
    loglines = _time.strftime("\nGenerated at: %Y/%m/%d - %I:%M %p\n")
    loglines += '          by: ' + _pwd.getpwuid(_os.getuid())[0] + ' (' + username + ')\n\n'
    writeLog(path, loglines)
    
    with open('{0}/polt.log'.format(path), 'r') as fl:
        print('\n{0}\nPOLTOOLS LOG (content of {1}/polt.log)\n\n{2}'.format('='*40, path, fl.read()))

    return




#################################################
#################################################
#################################################
def corObjStd(night, f, calc, path=None, delta=3.5, verbose=True):
    """
    Find the correction factor delta theta for filter 'f'
    inside night 'night', for calcite 'calc' (the last variable
    must be the angle between the ord. and extraord. beams).
    Returns the values for matching standard stars inside
    'night'/std.dat file, except when marked with an 'E' flag.

    NEW: In case of missing data in filter 'f', this routine
         tries to compute delta theta by using the data from
         another filter, making use of the relations among the
         values of delta theta for each filter). But only the
         estimate presenting the smaller error will be taken.

    verbose: auxiliary variable used as False inside
             polt.genTarget routine. Keep it as True
             otherwise.

    Formulas:
        * dth = fact*th^std_measured - th^std_published
        * fact = -1  for observations before 2015 March 1st
                 +1  otherwise

    delta: tolerance, in degree, of angle of two beams for
           one same calcite (default: +/-3.5 degree)

    Output, in order: stdnames, mdth, smdth, flag, tags
      - stdnames: string with standard star names separated
                  by commas (,)
      - mdth: mean delta theta (corrected by the output factor
              +-1 of routine polt.thtFactor())
      - smdth: the error of the mean delta theta
      - flag: flag concerning to the standards.
      - tags: string with the tags concerning to the standards,
              separated by commas (,), or '---' if none.
    """

    ########################
    def computeDth():
        """
        Main routine
        """
        
        try:
            dthref = _np.loadtxt('{0}/refs/dths.txt'.format(_hdtpath()), dtype=str)
        except:
            print('# ERROR: Can\'t read files pyhdust/refs/dths.txt')
            raise SystemExit(1)

        # Try to use the standard star observation at filter f
        stdnames, dth, sdth, flag, tags = readFilter(f)
#        print('STD REPORT: filter {0} runned...\t stdnames={1}\t dth={2}\t sdth={3}\t flag={4}\t tags={5}'.format(f, stdnames, dth, sdth, flag, tags))

        # Case successful on find
        if stdnames != '---':
            mdth = sum(dth)/len(dth) # Mean dth
            devth = _np.std(dth)/_np.sqrt(len(dth))  # Std dev of the mean
            smdth = _np.sqrt(devth**2 + _np.dot(sdth,sdth)/(len(sdth)**2))  # Combined final error
#            print('STD REPORT: filter {0} has...\t dth={1}\t mdth=mean(dth)={2}\t smdth={3}'.format(f, dth, mdth, smdth))
        # Otherwise
        else:
 #           print('STD REPORT: filter {0} has not standard... trying compute dth from another filter...'.format(f))
#            print('{0:<10s} Trying compute delta_theta from another filter...'.format(f))
            mdth = 0.
            smdth = 10000.
 
            # To identify the calcite I propose the if below, because is just like the beams
            # seem to behave within the CCD.
#            if calc == 0:
#                calcite = ''
#                print('{0:<12s} WARNING! Calcite name not identified for an angle of 0 degrees.'.format(night+', '+f+':'))
#                while calcite not in ('a0','a2'):
#                    calcite = raw_input('      Type the calcite name (a0/a2): ')
            if (calc < 12 and calc >= 0) or (calc > 78 and calc < 102) or (calc > 168 and calc < 180):
                calcite = 'a2'
            elif calc >= 0 and calc < 180:
                calcite = 'a0'
            # Useful because sometimes I set the calc value manually inside std.dat/obj.dat for special
            # cases, putting a negative value.
            else:
                calcite = ''
                print('{0:<12s} WARNING! Calcite name not identified for an angle of {1} degrees.'.format(night+', '+f+':', calc))
                while calcite not in ('a0','a2'):
                    calcite = raw_input('      Type the calcite name (a0/a2): ')

            for filt in filters:

                if filt == f[0]:
                    continue

                ddth = 0.
                sddth = 0.
                stdnamesaux, dth, sdth, flagaux, tagsaux = readFilter(filt)
#                print('STD REPORT: filter {0} runned to correction of filter {1}...\t stdnames={2}\t dth={3}\t sdth={4}\t flag={5}\t tags={6}'.format(filt, f, stdnamesaux, dth, sdth, flagaux, tagsaux))
                if stdnamesaux == '---':
#                    print('STD REPORT: filter {0} (runned to correction of filter {1}) has not standard... trying compute dth from another filter...'.format(filt, f))
                    continue

                # Try to find the 'delta delta theta'
                for dthi in dthref:

                    if dthi[0] == '{0}-{1}'.format(f[0],filt) and dthi[1] == calcite:

                        ddth = float(dthi[2])
                        sddth = float(dthi[3])

                        dth = [di+ddth for di in dth]               # Refresh the new dth list
                        mdthaux  = sum(dth)/len(dth)                # Mean dth
                        devthaux = _np.std(dth)/_np.sqrt(len(dth))  # Std dev of the mean
                        smdthaux = _np.sqrt(devthaux**2 + _np.dot(sdth,sdth)/(len(sdth)**2) + sddth**2)  # Combined final error

#                        print('STD REPORT: filter {0} (runned to correction of filter {1}) has the identified {2}-{3} value for calcite {4} ({5}): ddth={6}\t sddth={7}'.format(filt, f, f, filt, calc, calcite, ddth, sddth))

                        # Case there exists some observation in filter filt, as well as the study
                        # of delta that inside sths.txt file, compare if it have a best error
                        if smdthaux < smdth:
                            stdnames = stdnamesaux
                            flag = flagaux
                            tags = tagsaux
                            mdth = mdthaux
                            devth = devthaux
 #                           print('STD REPORT: filter {0} (runned to correction of filter {1}) has a best error in mdth...\t dth={2}\t ddth={3}\t mdth=mean(dth)+ddth={4}\t smdth={5}<{6}. Refreshing...'.format(filt, f, [di-ddth for di in dth], ddth, mdth, smdthaux, smdth))
                            smdth = smdthaux
 #                       else:
#                            print('STD REPORT: filter {0} (runned to correction of filter {1}) has NOT a best error in mdth: smdth={2}>={3}. Skipping this filter...'.format(filt, f, smdthaux, smdth))
                        break

#                if ddth == 0 and sddth == 0:
#                    print('STD REPORT: filter {0} (runned to correction of filter {1}) doesn\'t have {2}-{3} ddth value found for calcite {4} ({5}). Skipping this filter...'.format(filt, f, f, filt, calc, calcite))

            if stdnames == '---':
                if verbose:
                    print(('{0:<12s} WARNING! None standard found in `std.dat` for filter {1}.').format(night+', '+f+':', f))
                smdth = 0.
            else:
                if flag == 'OK':
                    flag = 'W'
                if tags == '---':
                    tags = 'oth-dth'
                else:
                    tags += ',oth-dth'
                if verbose:
                    print(('{0:<12s} WARNING! Delta theta for filter {1} was computed using an delta theta from another filter.').format(night+', '+f+':', f))
                else:
                    print(('{0:<12s} WARNING! Delta theta for filter {1} was computed using an delta theta from another filter.').format('', f))

 #       print('STD REPORT: filter {0} - final values...\t stdnames={1}\t mdth={2}\t smdth={3}\t flag={4}\t tags={5}\n'.format(f, stdnames, mdth, smdth, flag, tags))

        return stdnames, mdth, smdth, flag, tags
        ########################


    ########################
    def readFilter(filt):
        """
        Receive a filter 'filt' and return two lists with the delta theta values
        and the errors. Return also a string with the names of standards and
        the flag and tags string.
        """

        try:
            stdref = _np.loadtxt('{0}/refs/pol_padroes.txt'.format(_hdtpath()), dtype=str, usecols=range(0,22))
        except:
            print('# ERROR: Can\'t read files pyhdust/refs/pol_padroes.txt')
            raise SystemExit(1)

        dth = []
        sdth = []
        stdnames = '---'
        tags = '---'
        flag = 'OK'

        if _os.path.exists('{0}/{1}/std.dat'.format(path,night)):
            stds = _np.loadtxt('{0}/{1}/std.dat'.format(path,night), dtype=str)

            if len(stds) > 0 and len(stds[-1]) != 9:
                stds = stds.reshape(-1,9)

            for stdinf in stds:
                if stdinf[7] == 'E' and stdchk(stdinf[2])[0] and stdinf[3] == filt and \
                                                        abs(float(stdinf[4])-calc) <= delta:
                    if f == filt and verbose:
                        print(('{0:<12s} WARNING! Standard `{1}` ({2}) wasn\'t used because it ' +\
                                    'had `E` flag. Skipping this standard data...').format(night+', '+f+':',stdinf[2],f))
#                        else:
#                            print(('{0:<12s} WARNING! Standard `{1}` ({2}) wasn\'t used because it ' +\
#                                    'had `E` flag. Skipping this obs serie...').format('',stdinf[2],f))
                    continue
                elif stdinf[7] != 'E' and stdchk(stdinf[2])[0] and stdinf[3] == filt and \
                                                        abs(float(stdinf[4])-calc) <= delta:
                    # Bednarski: Show error message now
                    try:
                        nstar = int(stdinf[6])
                        Q, U, sig, P, th, sigT, tmp, tmp2 = readout('{0}/{1}'.\
                                                format(path+'/'+night,stdinf[5]), nstar=nstar)
                    except:
                        if f == filt and verbose:
                            print(('{0:<12s} WARNING! Standard `{1}` ({2}) wasn\'t used because' +\
                                                ' can\'t open/read {3}. Skipping this standard data...').\
                                                format(night+', '+f+':', stdinf[2], filt, stdinf[5]))
 #                           else:
 #                               print(('{0:<12s} WARNING! Standard `{1}` ({2}) wasn\'t used because' +\
 #                                               ' can\'t open/read {3}. Skipping this obs serie...').\
  #                                              format('', stdinf[2], filt, stdinf[5]))
                        continue
                    if stdref[stdchk(stdinf[2])[1],filters.index(filt[0])+1] == '0':
                        if f == filt and verbose:
                            print(('{0:<12s} WARNING! Standard `{1}` ({2}) wasn\'t used because' +\
                                ' there is no published value in such filter. Skipping this standard data...').\
                                                format(night+', '+f+':', stdinf[2], filt))
                        continue
                        
                    # Refresh string for std names
                    if stdnames == '---':
                        stdnames = stdinf[2]
                    elif stdinf[2] not in stdnames:
                        stdnames += ','+stdinf[2]

                    # Receive the published theta and its error
                    # (Let filt[0] because sometimes filter can be 'v2' for example)
                    i = stdchk(stdinf[2])[1]
                    angref = float(stdref[i,filters.index(filt[0])+1])
                    sangref = float(stdref[i,filters.index(filt[0])+6])

                    # Calculate dth and its error
                    dth += [thtFactor(float(stdinf[0]))*float(th)-angref]
                    while dth[-1] >= 180:
                        dth[-1] -= 180
                    while dth[-1] < 0:
                        dth[-1] += 180
                    sdth += [_np.sqrt((28.65*float(sig)/float(P))**2 + sangref**2)]

#                    print('STD:    dth({0}) = factor*th_obs - th_pub = {1} * {2:.2f} - {3:.2f} = {4:.3f}'.format(len(dth),thtFactor(float(stdinf[0])),float(th), angref, dth[-1]))
 #                   print('STD:    s_th_pub={0:.3f}'.format(sangref))
                    
                    # Receive the flag
                    if flag == 'OK' and stdinf[7] == 'W':
                        flag = 'W'

                    # Refresh the tag list
                    for tagi in stdinf[8].split(','):
                        if tagi not in tags+',---':
                            if tags == '---':
                                tags = ''
                            else:
                                tags += ','
                            tags += tagi

            # Fixes the cases where dth is close to 0 (for instance, if dth[0]=0.3,
            # dth[1] must be -1.0 and not 179.0
            for i in range(len(dth)):
                if min(dth) < 10 and dth[i] > 170:
                    dth[i] -= 180

            if stdnames == '---':
                flag = 'W'
                tags = 'no-std'
                   
        return stdnames, dth, sdth, flag, tags
        ########################


    if path == None or path == '.':
        path = _os.getcwd()

    calc = float(calc)

    if not _os.path.exists('{0}/{1}/std.dat'.format(path,night)):
#        print('{0:<12s} WARNING! `std.dat` file not found.'.format(night+', '+f+':'))
        return '---', 0., 0., 'W', 'no-std'
    else:
        return computeDth()




#################################################
#################################################
#################################################
def genTarget(target, path=None, path2=None, ispol=None, skipdth=False, delta=3.5, epssig=2.0):
    """ Gen. target

    Generate a table with all observations found for 'target',
    unless those which have `E` flag.

    The error in delta theta (factor from observation of standard
    star) IS NOT propagated to the theta of Be star.

    skipdth: skip the observations without estimates for delta
             theta (no standard star in none filter and no
             standard star in previous and next nights)? Case
             True, the exception is only when the data is
             'unpolarized' (defined by 'epssig' variable).
    epssig: sigP/P max for unpolarized target (sigP/P up to
            epssig doesn't need standard star when skipdth=True)
    ispol: the Serkowski parameters [P_max, lambda_max, theta_IS]
            to correct IS polarization (P_max ein % and lambda_max in
            Angstrom). If ispol==None, don't make the correction
            of ISP.
    
    Syntax of out tags:   tags1:tags2:tags3, where tags1 is concerning
                          to the method to calculate delta_theta
                          correction factor, tags2 is the tag list of
                          the observation of object and tags3 is the
                          tags of the standard star that has been used.


    If no standard star was found in some night for one filter, this
    routine tries to use the standard from another filter to compute
    the delta theta in missing filter. Case there are no standard star
    in the other filters also, the routine tries to read a file named
    as std.link, whose lines content must be [calc   night]. It points to
    a standard from another night `night` in calcite `calc` (an example
    of sintax: [calc   night] = [136.1   15out22]). Caso this file was
    not found, the routine queries directly of what night to use the
    standard stars.

    Formulas:
        * th_eq = fact*th_measured - fact*th^std_measured + th^std_published
        * th_eq = fact*th_measured - dth
        * dth = fact*th^std_measured - th^std_published
        * fact = -1 for observations before 2015 March 1st
                 +1 otherwise
        * Q_eq = P*cos(2*th_eq*pi/180)
        * U_eq = P*sin(2*th_eq*pi/180)
        * sigth = 28.65*sig/P

    """

    verbose = ''
    nlines = 0
    
    if path == None or path == '.':
        path = _os.getcwd()
    if path2 == None or path2 == '.':
        path2 = _os.getcwd()

    print(target, path)

    # Read lists and verify if target is a valid target
    try:
        obj = _np.loadtxt('{0}/refs/pol_alvos.txt'.format(_hdtpath()), dtype=str)
        if _os.path.exists('{0}/refs/pol_hip.txt'.format(_hdtpath)):
            obj = _np.concatenate((obj,_np.loadtxt('{0}/refs/pol_hip.txt'.format(_hdtpath()), dtype=str)))
        if _os.path.exists('{0}/refs/pol_unpol.txt'.format(_hdtpath)):
            obj = _np.concatenate((obj,_np.loadtxt('{0}/refs/pol_unpol.txt'.format(_hdtpath()), dtype=str)))
        std = _np.loadtxt('{0}/refs/pol_padroes.txt'.format(_hdtpath()), dtype=str, usecols=range(0,22))
    except:
        print('# ERROR: Can\'t read files pyhdust/refs/pol_alvos.txt and/or pyhdust/refs/pol_padroes.txt.')
        raise SystemExit(1)

    if target in std[:,0]:
        ftype='std'
    else:
        ftype='obj'

    if target not in _np.hstack((std[:,0],obj)) and 'field' not in target:
        print('\nWARNING: Target {0} is not a default target or standard!'.\
        format(target))

    print('\n'+'='*30+'\n')

    nights = [fld for fld in _os.listdir(path) if _os.path.isdir(_os.path.join(path, fld))]
    lines = ''

    for night in nights:

        # Check obj.dat/std.dat for the night
        if _os.path.exists('{0}/{1}/{2}.dat'.format(path,night,ftype)):
            try:
                objs = _np.loadtxt('{0}/{1}/{2}.dat'.format(path,night,ftype), dtype=str)
            except:
                print('{0:<12s} ERROR! Can\'t read {1}.dat file. Ignoring this night...\n'.format(night+':',ftype))
                continue

            # Verify if std has more than one line. Case not, do the reshape
            if _np.size(objs) == 9:
                objs = objs.reshape(-1,9)
            elif _np.size(objs) % 9 != 0:
                print('{0:<12s} ERROR! Wrong column type in {1}.dat file. Ignoring this night...\n'.format(night+':',ftype))
                continue

            valc = True
            
            # Loop on found nights
            for objinf in objs:

                if objinf[2] == target or ('field' in objinf[2] and target == 'field'):
                    tags = ['---','---','---']
                    MJD, ccd, obj, f, calc, out, nstar, flag, tags[1] = objinf
                    if flag == 'E':
                        print(('{0:<12s} WARNING! Star found ({1}), but with `E` flag ' +\
                                        'and tags `{2}`. Ignoring this data...').format(night+', '+f+':',f,tags[1]))
                        continue
                    try:
                        # Fator is a var to indicate when polarization angle must be taken as 180-theta or +theta
                        fator = thtFactor(float(MJD))
                        Q, U, sig, P, th, sigT, tmp, tmp2 = readout('{0}/{1}'.\
                                                    format(path+'/'+night,out), nstar=int(nstar))
                    except:
                        print('{0:<12s} ERROR! Can\'t open/read out file {1}. Ignoring this data...\n'.format(night+', '+f+':',out))
                        continue

                    P = float(P)*100
                    th = float(th)
                    sig = float(sig)*100
                    sigth = 28.65*sig/P
#                    print objinf, objinf[2], target, tags[1]

                    # Try to get the night's standard
                    if ftype == 'obj':
                        # Print below the warning message and only one time by night
                        if valc and not _os.path.exists('{0}/{1}/std.dat'.format(path,night)):
                            print('{0:<12s} WARNING! `std.dat` file not found.'.format(night+':'))
                            valc = False
                        stdnames, mdth, smdth, flagstd, tags[2] = corObjStd(night, f, calc, path=path, delta=delta)
                    else:
                        stdnames = '---'
                        mdth, smdth = 0, 0
                        flagstd, tags[2] = 'OK', '---'

#                    if flagstd == 'E':
#                        mdth = 0.
#                        smdth = 0.
#                        stdnames = '---'
#                        flagstd = 'W'
#                        tags[2] = 'no-std'

                    vald=True
                    # APPLY ALTERNATIVE METHOD TO COMPUTE DTHETA IN CASES WHERE THERE IS NO NIGHT'S STANDARD
                    while stdnames == '---' and ftype == 'obj':

                        night_alt=''
#                        print '{0}/{1}/std.link'.format(path,night)
#                        print _os.path.exists('{0}/{1}/std.link'.format(path,night))
                        if vald and _os.path.exists('{0}/{1}/std.link'.format(path,night)):
#                            print 'entrou'
                            file0 = _np.loadtxt('{0}/{1}/std.link'.format(path,night), dtype=str)
                            if type(file0[0]) != _np.ndarray and _np.size(file0) == 2:
                                file0 = file0.reshape(-1,2)
                            for line0 in file0:
                                if abs(float(line0[0])-float(calc)) <= delta:
                                    night_alt = line0[1]
                                    break
                            # if temporary for me.
                        if night_alt == 's' or _os.path.exists('{0}/{1}/skipstd'.format(path,night)):
                            print(('{0:<12s} WARNING! No standard correction as specified inside std.link.\n').format(night+', '+f+':', night_alt))
                            break
                        if night_alt=='':
                            night_alt = raw_input('\n{0:<12s} Do you want to select some standard from another day?\n{0:<12s} #Type the date or `s` to skip: '.format('','#'))
                            print('')
                            if night_alt in ('s','S'):
                                break

                        if night_alt!='' and _os.path.exists('{0}/{1}'.format(path,night_alt)):
                            stdnames, mdth, smdth, flagstd, tags[2] = corObjStd(night_alt, f, calc, path=path, delta=delta, verbose=False)
                            valc = False
                            if stdnames != '---':
                                if flagstd == 'OK':
                                    flagstd = 'W'
                                if tags[2] == '---':
                                    tags[2] = 'oth-day-std'
                                else:
                                    tags[2] += ',oth-day-std'
                                if vald:
                                    print(('{0:<12s} WARNING! Using standard from another night ({1}) as specified inside std.link.\n').format(night+', '+f+':', night_alt))
                            elif vald: 
                                print(('\n{0:<12s} ERROR! Standard not found inside the alternative night {1}!').format(night+', '+f+':', night_alt))
                                vald = False
                            else:
                                print(('\n{0:<12s} ERROR! Standard not found inside the alternative night {1}!').format('', night_alt))
                        elif vald: 
                            print(('\n{0:<12s} ERROR! Standard not found inside the alternative night {1}!').format(night+', '+f+':', night_alt))
                            vald = False
                        else:
                            print(('\n{0:<12s} ERROR! Standard not found inside the alternative night {1}!').format('', night_alt))
#                        print stdname, thstd, angref, flagstd, tags[2]

                    # Set the tags concerning to the standard
                    if stdnames == '---' and ftype == 'obj':
                        tags[0] = 'no-std'
                    elif ftype == 'obj':
                        if 'oth-day-std' in tags[2]:
                            tags[0] = 'oth-day-std'
                        if 'oth-dth' in tags[2] and tags[0] == '---':
                            tags[0] = 'oth-dth'
                        elif 'oth-dth' in tags[2]:
                            tags[0] += ',oth-dth'

                    # Refresh tags and flags
                    for specialtag in ('no-std','oth-day-std','oth-dth'):
                        for i in (1,2):
                            if tags[i] == specialtag:
                                tags[i] = '---'
                                if i == 1:
                                    flag = 'OK'
                            elif tags[i][0:7] == specialtag+',':
                                tags[i] = tags[i].replace(specialtag+',','')
                            else:
                                tags[i] = tags[i].replace(','+specialtag,'')

                    # Set the "global" flag (for object+standard)
                    if flag == 'E' or flagstd == 'E':
                        flag = 'E'
                    elif flag == 'W' or flagstd == 'W':
                        flag = 'W'
                    else:
                        flag = 'OK'

                    # Applying the correction of standard star
                    th = fator*th-mdth

                    # Fixing the angle value and computing QU parameters
                    while th >= 180:
                        th-= 180
                    while th < 0:
                        th+= 180
                    Q = P*_np.cos(2*th*_np.pi/180)
                    U = P*_np.sin(2*th*_np.pi/180)

                    # Correction of IS polarization
                    if ispol != None:
                        QIS, UIS = serkowski(ispol[0], ispol[1], str(f), mode=1, pa=ispol[2])
                        Q = Q - QIS
                        U = U - UIS
                        P = _np.sqrt(Q**2 + U**2)
                        th = _np.arctan(Q/U)*90/_np.pi
#                        sigth = 28.65*sig/P

                        # Fix the angle to the correct in QU diagram
                        if Q < 0:
                            th += 90
                        elif Q >= 0 and U < 0:
                            th += 180
                

                    # Write the line
                    if stdnames != '---' or (not skipdth) or P/sig <= epssig or ftype == 'std':
                        if out.find('_WP') == -1:
                            outn = out[-11:]
                        else:
                            outn = out[out.find('_WP')-7:]
                        lines += ('{:12s} {:>7s} {:>7s} {:>4s} {:>5s} {:>12s} {:>6.1f} {:>6.1f}'+
                                ' {:>8.4f} {:>8.4f} {:>8.4f} {:>7.2f} {:>7.4f} '+
                                '{:>6.2f} {:>13s} {:>4s} {:>5s} {:>s}').format(MJD, night, ccd, f, \
                                calc, stdnames, mdth, smdth, P, Q, U, th, sig, sigth, outn, nstar, \
                                flag, ';'.join(tags))
                        if target == 'field':
                            lines += '   {0}\n'.format(obj)
                        else:
                            lines += '\n'
                            
                        nlines += 1
                    else:
                        print(('{0:<12s} ERROR! No valid delta_theta value estimated in filter {1}.' +\
                                        ' Ignoring this data...\n').format(night+', '+f+':', f))
        else:
            if verbose != '':
                verbose += ', '
            verbose += night

    # Print "no obj/std.dat found" message
    if verbose != '':
        print('{0:<12s} WARNING! No `{1}.dat` file found for the following nights: {2}. Ignoring these nights...'.format('-------',ftype,verbose))

    print('\n'+'='*30+'\n')

    # Write the output
    if lines != '':
        if ispol==None:
            f0 = open('{0}/{1}.log'.format(path2,target),'w')
        else:
            f0 = open('{0}/{1}_iscor.log'.format(path2,target),'w')

        if target == 'field':
            lines = ('#{:>11s} {:>7s} {:>7s} {:>4s} {:>5s} {:>12s} {:>6s} {:>6s}' +\
                        ' {:>8s} {:>8s} {:>8s} {:>7s} {:>7s} {:>6s} {:>13s}' +\
                        ' {:>4s} {:>5s} {:>s}   {:s}\n').format('MJD', 'night',\
                        'ccd', 'filt', 'calc', 'stdstars', 'dth', 'sigdth', 'P', 'Q', 'U',\
                        'th', 'sigP', 'sigth', 'outfile', 'star', 'flag', 'tags', 'obj_name')+lines
        else:
            lines = ('#{:>11s} {:>7s} {:>7s} {:>4s} {:>5s} {:>12s} {:>6s} {:>6s}' +\
                        ' {:>8s} {:>8s} {:>8s} {:>7s} {:>7s} {:>6s} {:>13s}' +\
                        ' {:>4s} {:>5s} {:>s}\n').format('MJD', 'night',\
                        'ccd', 'filt', 'calc', 'stdstars', 'dth', 'sigdth', 'P', 'Q', 'U',\
                        'th', 'sigP', 'sigth', 'outfile', 'star', 'flag', 'tags')+lines

        if ispol==None:
            ispol = [0,0,0]
        lines = ('# ISP parameters used:\n#\n# Pmax (%)   lmax (A)     PA\n# {0:>8.4f} {1:>9.2f} {2:>7.2f}\n#\n')\
                                                .format(ispol[0],ispol[1],ispol[2]) + lines
        f0.writelines(lines)
        f0.close()
        
        if ispol==[0,0,0]:
            print('DONE! {0} lines written in {1}/{2}.log.'.format(nlines,path2,target))
        else:
            print('DONE! {0} lines written in {1}/{2}_iscor.log.'.format(nlines,path2,target))
    else:
        print('NOT DONE! No valid observation was found for target `{0}`.'.format(target))
      
    return




#################################################
#################################################
#################################################
def fixISP(logfile, ispol, path2=None):
    """
    Correct the interstellar polarization in file logfile

    logfile: outfile from genTarget.
    ispol: the Serkowski parameters [P_max, lambda_max, theta_IS]
            to correct IS polarization (P_max ein % and lambda_max in
            Angstrom). If ispol==None, don't make the correction
            of ISP.
    path2:   path where to save the data file. If None, it is
             supposed the same of logfile.


    """

    star = _phc.trimpathname(logfile)[1].split('.')[0].split('_')[0]
    if star in _phc.bes:
        be = _phc.bes[star]
    else:
        be = star

    if path2 == None or path2 == '.':
        path2 = _os.getcwd()

#    if path2 == None or path2 == '.':
#        path2 = _phc.trimpathname(logfile)[0]
#        if path2=='':
#            path2 = _os.getcwd()


    try:
        lines = _np.loadtxt(logfile, dtype=str)
    except:
        print('# ERROR: Can\'t read file {0}.'.format(logfile))
        raise SystemExit(1)

    if type(lines[0]) != _np.ndarray and _np.size(lines) == 18:
        lines = lines.reshape(-1,18)

    linesout=''

    for i,line in enumerate(lines):

        # read filter
        f = line[3]

        # Correction of IS polarization
        QIS, UIS = serkowski(ispol[0], ispol[1], str(f), mode=1, pa=ispol[2])
        Q = float(line[9]) - QIS
        U = float(line[10]) - UIS
        P = _np.sqrt(Q**2 + U**2)
        th = _np.arctan(Q/U)*90/_np.pi
#                        sigth = 28.65*sig/P

        # Fix the angle to the correct in QU diagram
        if Q < 0:
            th += 90
        elif Q >= 0 and U < 0:
            th += 180

        linesout += ('{:12s} {:>7s} {:>7s} {:>4s} {:>5s} {:>12s} {:>6.1f} {:>6.1f}'+
                                ' {:>8.4f} {:>8.4f} {:>8.4f} {:>7.2f} {:>7.4f} '+
                                '{:>6.2f} {:>13s} {:>4s} {:>5s} {:>s}\n').format(line[0], line[1], \
                                line[2], line[3], line[4], line[5], float(line[6]), float(line[7]), \
                                P, Q, U, th, float(line[12]), float(line[13]), line[14], line[15], \
                                line[16], line[17])

#        lines[i][8] = P
#        lines[i][9] = Q
#        lines[i][10] = U
#        lines[i][11] = th

    if linesout != '':
        f0 = open('{0}/{1}_iscor.log'.format(path2,star),'w')
        linesout = ('# ISP parameters used:\n#\n# Pmax (%)   lmax (A)     PA\n# {0:>8.4f}'+\
                    ' {1:>9.2f} {2:>7.2f}\n#\n').format(ispol[0],ispol[1],ispol[2]) + \
                   ('#{:>11s} {:>7s} {:>7s} {:>4s} {:>5s} {:>12s} {:>6s} {:>6s}' +\
                        ' {:>8s} {:>8s} {:>8s} {:>7s} {:>7s} {:>6s} {:>13s}' +\
                        ' {:>4s} {:>5s} {:>s}\n').format('MJD', 'night',\
                        'ccd', 'filt', 'calc', 'stdstars', 'dth', 'sigdth', 'P', 'Q', 'U',\
                        'th', 'sigP', 'sigth', 'outfile', 'star', 'flag', 'tags')+linesout

        f0.writelines(linesout)
        f0.close()

        print('DONE! File written in {0}/{1}_iscor.log.'.format(path2,star))
    else:
        print('NOT DONE! No observation for target `{0}`.'.format(star))


    return




#################################################
#################################################
#################################################
def serkowski(pmax, lmax, wlen, mode, pa=None, law='w82'):
    """

    Mode==1
    Receive ISP parameters 'pmax', 'lmax' e 'pa' and return
    the Stokes QU parameters concerning to the value
    of Serkowski's law at wavelenght 'wlen'.

    Mode==2
    Receive ISP parameters 'pmax', 'lmax' and return
    the P value computed by the Serkowski's law at
    wavelenght 'wlen'.
    

    'wlen' and 'lmax' must be in Angstrom.
    'wlen' can be 'u'/'b'/'v'/'r'/'i' also.

    'law' defines what value use to K parameter:

    w82  -  Wilking (1982)
        K = 1.86*lmax - 0.1
    w80  -  Wilking (1980)
        K = 1.68*lmax - 0.002
    serk  -  Serkowski
        K = 1.15

    Serkowski's Law:
        P = pmax*np.exp(-K*np.log(lmax/wlen)**2)
    
    """

    if law=='w82':
        K = 1.86*lmax/10000 - 0.1     # Wilking (1982)
    elif law=='w80':
        K = 1.68*lmax/10000 - 0.002   # Wilking (1980)
    elif law=='serk':
        K = 1.15                # Serkowski

    if pmax==0 and lmax==0:
        P = 0.
    elif mode==1 and wlen in _phc.lbds:
        P = pmax*_np.exp(-K*_np.log(lmax/_phc.lbds[wlen])**2)
    else:
        P = pmax*_np.exp(-K*_np.log(lmax/wlen)**2)

    if mode == 1:
#    print wlen, P, pa
        Q = P*_np.cos(2*pa*_np.pi/180)
        U = P*_np.sin(2*pa*_np.pi/180)
        return Q, U
    elif mode == 2:
        return P
    else:
        return




#################################################
#################################################
#################################################
### IMPLEMENTAR CORLOWPOL QUE SERÁ A CORREÇÃO
### DAS INCERTEZAS PARA POLARIZAÇÕES BAIXAS
def propQU(p, th, sp, sdth, estim='wk'):
    """
    Propagate the delta theta error over the polarization
    angle, computing the new errors for theta, Q and U.

    Input:

    - p, th: lists holding the P and theta values.
    - sp, sdth: lists holding the P and delta theta errors.

    Return lists containing the new errors for theta,
    Q and U.: sth, sq, su.

    Formulas:
    
    sth = sqrt( sth0^2 + sdth^2 )
    sq = sqrt( (cos(2*th)*sp)^2 + (p*sin(2*th)*sth)**2 )
    su = sqrt( (sin(2*th)*sp)^2 + (p*cos(2*th)*sth)**2 )


    Unbias theta error using 'estim' estimator:

        if p/sp <= K,   sth0 = psi
        otherwise,      sth0 = propagated error

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

    sth,sq,su = [],[],[]
    sth0=0
    
    for i in range(len(p)):

        if p[i] != 0:
            if estim!='mts':
                if sp[i]!=0 and p[i]/sp[i] > k:
                    sth0 = 28.65*sp[i]/p[i]
                else:
                    sth0 = 51.96
            else:
                if sp[i]!=0 and p[i]/sp[i] > 6:
                    sth0 = 28.65*sp[i]/p[i]
                elif sp[i]!=0:
                    a=32.50
                    b=1.350
                    c=0.739
                    d=0.801
                    e=1.154
                    sth0 = a*(b+_np.tanh(c*(d-p[i]/sp[i]))) - e*p[i]/sp[i]
                else:
                    sth0 = 61.14

            sth += [_np.sqrt(sth0**2 + sdth[i]**2)]
        else:
            if estim!='mts':
                sth += [51.96]
            else:
                sth += [61.14]

        sq += [_np.sqrt( (_np.cos(th[i]*_np.pi/90)*sp[i])**2 + (p[i]*_np.sin(th[i]*_np.pi/90)*sth[i]*_np.pi/90)**2 )]
        su += [_np.sqrt( (_np.sin(th[i]*_np.pi/90)*sp[i])**2 + (p[i]*_np.cos(th[i]*_np.pi/90)*sth[i]*_np.pi/90)**2 )]
       
 
    return sth, sq, su

    


#################################################
#################################################
#################################################
def genJD(path=None):
    """Generate de JD file for the fits inside the folder
    """
    if path == None or path == '.':
        path = _os.getcwd()
    for f in filters:
        lfits = _glob('*_{0}_*.fits'.format(f))
        if len(lfits) > 0:
            lfits.sort()
            i = lfits[0].find('_{0}_'.format(f))
            pref = lfits[0][:i+2]
            lfits = _glob('{0}_*.fits'.format(pref))
            lfits.sort()
            if len(lfits)%8 != 0:
                print('# Warning! Strange number of fits files!')
                print(lfits)
            JDout = ''
            i = 0
            for fits in lfits:
                i += 1
                imfits = _pyfits.open(fits)
                dtobs = imfits[0].header['DATE']
                if 'T' in dtobs:
                    dtobs, tobs = dtobs.split('T')
                    dtobs = dtobs.split('-')
                    tobs = tobs.split(':')
                    tobs = float(tobs[0])*3600+float(tobs[1])*60+float(tobs[2])
                    tobs /= (24*3600)
                else:
                    print('# ERROR! Wrong DATE-OBS in header! {0}'.format(fits))
                    raise systemExit(1)
                JD = _np.sum(_jdcal.gcal2jd(*dtobs))+tobs
                JDout += 'WP {0}  {1:.7f}\n'.format(i,JD)
            f0 = open('JD_{0}'.format(pref),'w')
            f0.writelines(JDout)
            f0.close()
    return




#################################################
#################################################
#################################################
def listNights(path, tgt):
    """
    List Nights
    """
    ltgts = _np.loadtxt('{0}/refs/pol_alvos.txt'.format(_hdtpath()), dtype=str)
    lstds = _np.loadtxt('{0}/refs/pol_padroes.txt'.format(_hdtpath()), dtype=str,\
    usecols=[0])
    if tgt not in _np.hstack((ltgts,lstds)):
        print('# Warning! Target {0} is not a default target or standard!!'.\
        format(tgt))

    lnights = []
    nights = [fld for fld in _os.listdir(path) if \
    _os.path.isdir(_os.path.join(path, fld))]

    for night in nights:
        tgts = [fld for fld in _os.listdir('{0}/{1}'.format(path,night)) if \
    _os.path.isdir(_os.path.join('{0}/{1}'.format(path,night), fld))]
        if tgt in tgts:
            lnights += [night]
        else:
            out = [obj for obj in tgts if obj.find(tgt) > -1]
            if len(out) > 0:
                lnights += [night]

    return lnights




#################################################
#################################################
#################################################
# Falta alterar para novos índices das colunas dos arquivos std.dat e obj.dat
# das mudanças que fiz. Bednarski.
def plotMagStar(tgt, path=None):
    """ Function doc

    @param PARAM: DESCRIPTION
    @return RETURN: DESCRIPTION
    """
    if path == None or path == '.':
        path = _os.getcwd()
    lmags = _np.loadtxt('{0}/refs/pol_mags.txt'.format(_hdtpath()), dtype=str)

    if tgt not in lmags[:,0]:
        print('# ERROR! {0} is not a valid mag. star!'.format(tgt))
        return

    data = _np.loadtxt('{0}/{1}.log'.format(path,tgt), dtype=str)
    data = _np.core.records.fromarrays(data.transpose(), names='MJD,night,filt,\
    calc,ang.ref,dth,P,Q,U,th,sigP,sigth', formats='f8,a7,a1,f8,f8,f8,f8,\
    f8,f8,f8,f8,f8')
    
    if False:
        fig = _plt.figure()#figsize=(5.6,8))
        ax0 = fig.add_subplot(311)
        ax0.errorbar(data['MJD'], data['P'], data['sigP'], color='black')
        ax1 = fig.add_subplot(312)
        ax1.errorbar(data['MJD'], data['Q'], data['sigP'], color='blue')
        ax2 = fig.add_subplot(313)
        ax2.errorbar(data['MJD'], data['U'], data['sigP'], color='red')

    idx = _np.where(lmags[:,0] == tgt)
    P, ph0 = lmags[idx][0][1:]
    ph0 = float(ph0) - _jdcal.MJD_0
    
    phase = data['MJD']-ph0
    phase /= float(P)
    phase = _np.modf(phase)[0]
    idx = _np.where(phase < 0)
    phase[idx] = phase[idx]+1

    fig2, (ax0, ax1, ax2, ax3) = _plt.subplots(4,1, sharex=True)
    ax0.errorbar(phase, data['P'], yerr=data['sigP'], fmt='ok')
    ax1.errorbar(phase, data['Q'], yerr=data['sigP'], fmt='or')
    ax2.errorbar(phase, data['U'], yerr=data['sigP'], fmt='ob')
    ax3.errorbar(phase, data['th'], yerr=data['sigth'], fmt='ok')
    ax3.set_xlabel('Phase')
    ax0.set_ylabel(u'P (%)')
    ax1.set_ylabel(u'Q (%)')
    ax2.set_ylabel(u'U (%)')
    ax3.set_ylabel(r'$\theta$ (deg.)')
    ax0.set_title('{0} ; P={1} days, ph0={2:.3f}'.format(tgt,P,ph0+_jdcal.MJD_0))
    _plt.savefig('{0}/{1}.png'.format(path,tgt))

    bphase, bP, bsigP = _phc.bindata(phase, data['P'], data['sigP'], 30)
    bphase, bQ, bsigP = _phc.bindata(phase, data['Q'], data['sigP'], 30)
    bphase, bU, bsigP = _phc.bindata(phase, data['U'], data['sigP'], 30)
    bphase, bth, bsigth = _phc.bindata(phase, data['th'], data['sigth'], 30)
    fig3, (ax0, ax1, ax2, ax3) = _plt.subplots(4,1, sharex=True)
    ax0.errorbar(bphase, bP, yerr=bsigP, fmt='ok')
    ax1.errorbar(bphase, bQ, yerr=bsigP, fmt='or')
    ax2.errorbar(bphase, bU, yerr=bsigP, fmt='ob')
    ax3.errorbar(bphase, bth, yerr=bsigth, fmt='ok') 
    ax3.set_xlabel('Phase')
    ax0.set_ylabel(u'P (%)')
    ax1.set_ylabel(u'Q (%)')
    ax2.set_ylabel(u'U (%)')
    ax3.set_ylabel(r'$\theta$ (deg.)')
    ax0.set_title('{0} (binned); P={1} days, ph0={2:.3f}'.format(tgt,P,ph0+_jdcal.MJD_0))
    _plt.savefig('{0}/{1}_bin.png'.format(path,tgt))
    #_plt.show()
    return




#################################################
#################################################
#################################################
def sortLog(filename):
    """ Sort the *.out file """
    f0 = open(filename)
    lines = f0.readlines()
    f0.close()
    log = _np.loadtxt(filename, dtype=str)
    log = log[log[:,0].argsort()]
    fmt = '%12s %7s %1s %5s %5s %6s %5s %6s %6s %6s %5s %5s'
    _np.savetxt(filename.replace('.log','.txt'), log, fmt=fmt, header=lines[0])
    return




#################################################
#################################################
#################################################
def filtra_obs(n,obs):
    """ ### FILTER OBSERV. ### """
    nobs = [ ]
    for i in range(len(obs)):
        if obs[i][5]/obs[i][3] > n:
            nobs = nobs+[obs[i]]
    return _np.array(nobs)




#################################################
#################################################
#################################################
def filtraobs(data, r=20):
    """ filtro! """
    idx = data['P']/data['sigP'] > r
    return data[idx]




#################################################
#################################################
#################################################
def setCCD(fitsfile):
    """
    Set CCD name in global variable 'ccd'.

    The CCD name can be: 'ikon', 'ixon', '301' or
    'ikon-14912' (the last one is the Ikon CCD with
    Deep Depletion).
    """
    global ccd
    ccd = ''

    if fitsfile != '':
        try:
            fits = _pyfits.open(fitsfile)
            instrume = '{0}'.format(fits[0].header['SERNO'])
            if instrume.find('4335') != -1:
                ccd = 'ixon'
            elif instrume.find('4269') != -1:
                ccd = 'ixon'
            elif instrume.lower().find('10127') != -1:
                ccd = 'ikon'
            elif instrume.lower().find('9867') != -1:
                ccd = 'ikon'
            elif instrume.lower().find('14912') != -1:
                ccd = 'ikon-14912'
            else:
                ccd = ''
        except:
            pass

    while ccd not in ('ikon', 'ixon', '301', '654', 'ikon-14912'):
        ccd = raw_input('Type the CCD name (301/654/ikon/ikon-14912/ixon): ')




#################################################
#################################################
#################################################
def splitData(night, path_raw='raw', path_red='red'):
    """
    Split the raw files and reduced files for a night.

    Parameters:
        night: path to the night (this directory will be fully preserved)
        path_raw: directory with the raw data of the nights
        path_red: directory with the output files of reduction
    """

    print('')
    
    if not _os.path.exists(path_raw) or not _os.path.exists(path_red):
        print('Error: Directory \'{0}\' and/or \'{1}\' doesn\'t exist!\n'.format(path_raw, path_red))
        return 1
    elif not _os.path.exists(night):
        print('Error: Directory \'{0}\' doesn\'t exist!\n'.format(night))
        return 1

    if night[len(night)-1] == '/':
        night = night[:-1]

    if night.find('/') == -1:
        path = '.'
    else:
        path = _phc.trimpathname(night)[0]
        night = _phc.trimpathname(night)[1]
    
    #print path, night

    # Verify if the splitted directories exist
    if _os.path.exists('{0}/{1}'.format(path_raw, night)) or  \
                        _os.path.exists('{0}/{1}'.format(path_red, night)):
        while True:
            verif = raw_input('CAUTION: Directory \'{0}/{1}\' and/or \'{2}/{1}\' already exists! '.\
                           format(path_raw, night, path_red) + 'Are you sure to continue, ' + \
                           'overwriting all data inside these directories? (y/n): ')
            print('')
            if verif in ['y', 'Y']:
                
                if _os.path.exists('{0}/{1}'.format(path_raw, night)):
                    try:
                        _shutil.rmtree('{0}/{1}'.format(path_raw, night))
                    except:
                        print('ERROR when deleting directory \'{0}/{1}\'. Check the permissions.\n'.format(path_raw, night))
                        return 2

                if _os.path.exists('{0}/{1}'.format(path_red, night)):
                    try:
                        _shutil.rmtree('{0}/{1}'.format(path_red, night))
                    except:
                        print('ERROR when deleting directory \'{0}/{1}\'. Check the permissions.\n'.format(path_red, night))
                        return 2
                break
            elif verif in ['n', 'N']:
                print('Aborted!\n')
                return 1
            else:
               print('Value not valid!\n')

    # Loop for splitting
    for (direc, _, files) in _os.walk('{0}/{1}'.format(path, night)):
        for f in files:
            file_old = _os.path.join(direc, f)
           
            # Case a raw file
            if _re.search(r'[0-9][0-9].fits$', f) or _re.search(r'.doc$', f):
                file_new = path_raw + '/' + _os.path.relpath(file_old, path)
            # Case a file from reduction
            else:
                file_new = path_red + '/' + _os.path.relpath(file_old, path)

            # Create all subdirectories
            if not _os.path.exists(_phc.trimpathname(file_new)[0]):
                try:
                    _os.makedirs(_phc.trimpathname(file_new)[0])
                except:
                    print('ERROR when copying files. Check the permissions to directories' + \
                                ' \'{0}\' and \'{1}\'\n'.format(path_raw, path_red))
                    return 2

            _shutil.copy2(file_old, file_new)

    print('Done!\n')
    return 0




#################################################
#################################################
#################################################
def splitData301(night, path_raw='raw', path_red='red'):
    """
    Split the raw files and reduced files for a night, renaming according
       CCD iXon nomenclature.

    Parameters:
        night: path to the night (this directory will be fully preserved)
        path_raw: directory with the raw data of the nights
        path_red: directory with the output files of reduction
    """

    print('')

    # Verify if raw and red directories exist
    if not _os.path.exists(path_raw) or not _os.path.exists(path_red):
        print('Error: Directory \'{0}\' and/or \'{1}\' doesn\'t exist!\n'.format(path_raw, path_red))
        return 1
    elif not _os.path.exists(night):
        print('Error: Directory \'{0}\' doesn\'t exist!\n'.format(night))
        return 1

    if night[len(night)-1] == '/':
        night = night[:-1]
    if night.find('/') == -1:
        path = '.'
    else:
        path = _phc.trimpathname(night)[0]
        night = _phc.trimpathname(night)[1]

    # Verify if the splitted directories exist
    if _os.path.exists('{0}/{1}'.format(path_raw, night)) or  \
                        _os.path.exists('{0}/{1}'.format(path_red, night)):

        while True:
            verif = raw_input('CAUTION: Directory \'{0}/{1}\' and/or \'{2}/{1}\' already exists!\n'.\
                           format(path_raw, night, path_red) + 'Are you sure to continue, ' + \
                           'overwriting all data inside these directories? (y/n): ')
            print('')
            if verif in ['y', 'Y']:
                if _os.path.exists('{0}/{1}'.format(path_raw, night)):
                    try:
                        _shutil.rmtree('{0}/{1}'.format(path_raw, night))
                    except:
                        print('ERROR when deleting directory \'{0}/{1}\'. Check the permissions.\n'.format(path_raw, night))
                        return 2

                if _os.path.exists('{0}/{1}'.format(path_red, night)):
                    try:
                        _shutil.rmtree('{0}/{1}'.format(path_red, night))
                    except:
                        print('ERROR when deleting directory \'{0}/{1}\'. Check the permissions.\n'.format(path_red, night))
                        return 2
                break
            elif verif in ['n', 'N']:
                print('Aborted!\n')
                return 1
            else:
               print('Value not valid!\n')


    # Loop for splitting
    for (direc, _, files) in _os.walk('{0}/{1}'.format(path, night)):

        jds = [f for f in files if _re.search(r'^jd[0-9]*', f)]

        for f in files:
            file_old = _os.path.join(direc, f)

            # If filename is 'logfile'
            if f == 'logfile':
                file_new = path_red + '/' + _os.path.relpath(file_old, path)
            # If file is inside dark/flat/bias path
            elif file_old.find('/dark/') != -1 or file_old.find('/flat/') != -1 or \
                                            file_old.find('/bias/') != -1:
                if _re.search(r'.fits$', f):
                    file_new = path_red + '/' + _os.path.relpath(file_old, path)
                else:
                    file_new = path_raw + '/' + _os.path.relpath(file_old, path)
            # If filename is type p010001, fb01001, s001001, d001001, b01001
            elif _re.search(r'p[0-9]*$', f):
                file_new = path_raw + '/' + _os.path.relpath(file_old, path)
            # The other cases
            else:

                # Tries to get object and filter names (it can be a file of none object)
                elem = filter(None, file_old.split('/{0}/'.format(night))[1].split('/'))
                if len(elem) > 0:
                    obj = elem[0]

                if len(elem) > 1 and elem[1] in ('u','b','v','r','i'):
                    filt = elem[1]
                else:
                    filt = ''

                # Case sum*.dat files
                if _re.search(r'^sum', f):
#                    print f
                    if filt == '':
                        print(('\nERROR: No filter was identified for file {0}.\nPut this files ' +\
                        ' inside a subdir named with the filter name and run again!').format(file_old))
                        return 3
                    file_new = path_red + '/' + night + '/' + obj + '/sum_' + obj + '_' +\
                                                                      filt + '_0' + f[-8:]
                # Case w*.dat, w*.log, w*.out files
                elif _re.search(r'\.dat|\.log|\.out$', f):
                    if filt == '':
                        print(('# ERROR: No filter was identified for file {0}.\nPut this files ' +\
                        ' inside a subdir named with the filter letter and run again!').format(file_old))
                        return 3
                    file_new = path_red + '/' + night + '/' + obj + '/w' + obj + '_' +\
                                                                        filt + '_' + f[1:]
                # Case coord.*.ord files
                elif _re.search(r'^coord', f):
                    if filt == '':
                        print(('# ERROR: No filter was identified for file {0}.\nPut this files ' +\
                        ' inside a subdir named with the filter letter and run again!').format(file_old))
                        return 3
                    file_new = path_red + '/' + night + '/' + obj + '/coord_' + obj + '_' +\
                                                                        filt + f[-6:]
                # Case JD* files, pass and concatenate only at the end of this loop
                elif _re.search(r'^jd[0-9]*', f):
                    continue
                else:
                    file_new = path_red + '/' + _os.path.relpath(file_old, path)

            # Create all subdirectories
            if not _os.path.exists(_phc.trimpathname(file_new)[0]):
                try:
                    _os.makedirs(_phc.trimpathname(file_new)[0])
                except:
                    print('ERROR when copying files. Check the permissions to directories' + \
                                ' \'{0}\' and \'{1}\'\n'.format(path_raw, path_red))
                    return 2

            print('OLD:' + file_old)
            print('NEW:' + file_new + '\n')
            _shutil.copy2(file_old, file_new)

        # Concatenate all JD files now
        if len(jds) > 0:
            jds.sort()
            with open('{0}/{1}/{2}/JD_{3}_{4}'.format(path_red,night,obj,obj,filt), 'a') as f_out:
                for f in jds:
                    with open('{0}/{1}'.format(direc,f), 'r') as f_in:
                        f_out.write(f_in.readline())

    print('Done!\n')
    return 0




#################################################
#################################################
#################################################
def graf_t(logfile, path2=None, vfilter=['no-std'], save=False, extens='pdf', filt='v'):
    """
    Plot a P_V x t, theta_V x t and P_B/P_I x t graphs for the
    Be star in the logfile .log file (the outfile from
    polt.genTarget). Propagates error of standard star.

    If 'no-std' is in vfilter, no data with 'no-std' tag will be
    displayed, but the others filtered data will be showed
    with a 'x' symbol.
    If 'no-std' is not in vfilter, the data with 'no-std' will
    be displayed normally and the other filtered will be
    showed with a 'x' symbol.

    """

    ###########
    ## FUNC
    def polColor(dayp, jdp, p, s, dayp_filt, jdp_filt, p_filt, s_filt):
        """
        Return pb/pi from input lists
        """
       
        jd, pbpi, spbpi = [],[],[]
        for i in range(len(dayp[1])):
            for j in range(len(dayp[4])):
                # p[4][j] != 0 is to prevent division by 0.
                if dayp[1][i] == dayp[4][j] and p[4][j] != 0:
 #                   print dayp[1][i], dayp[4][j]
                    jd += [(jdp[1][i] + jdp[4][j])/2]
                    pbpi += [p[1][i]/p[4][j]]
                    spbpi += [pbpi[-1]*_np.sqrt((s[1][i]/p[1][i])**2 + (s[4][j]/p[4][j])**2)]
#                    print pbpi[-1], p[1][i], p[4][j]
#                    print ''

                    # This line below prevent to take the same point more once
                    p[4][j] = 0.
                    break

        return jd, pbpi, spbpi
        ###########


    ###########
    ## FUNC
    def plot(fig, axs):
        """
        Receive figure and an axes list (with three axes objects)and do the plots
        filt WITHOUT show or save the image.
        """

        cm = _plt.cm.gist_rainbow     # Setting the color map
        factor=0.7                   # Factor to fix the font sizes
            
        try:
            lines = _np.loadtxt(logfile, dtype=str)
        except:
            print('# ERROR: Can\'t read file {0}.'.format(logfile))
            raise SystemExit(1)

        if type(lines[0]) != _np.ndarray and _np.size(lines) == 18:
            lines = lines.reshape(-1,18)

#        ax.set_title('{0} filter'.format(filt.upper()), fontsize=fonts[0]*factor, verticalalignment='bottom')
#        ax.text(0.98, 0.9, '{0} filter'.format(filt.upper()), horizontalalignment='right', \
#                 verticalalignment='bottom', transform=ax.transAxes, fontsize=fonts[1]*factor)
        axs[2].set_xlabel(r'MJD', size=fonts[1]*factor)
        axs[0].set_ylabel(r'$P_{0}$ (%)'.format(filt.upper()), size=fonts[1]*factor)
        axs[1].set_ylabel(r'$P_B / P_I$', size=fonts[1]*factor)
        axs[2].set_ylabel(r'$\theta_{0}$ (degree)'.format(filt.upper()), size=fonts[1]*factor)


        dayp, jdp, p, s = [[],[],[],[],[]], [[],[],[],[],[]], [[],[],[],[],[]], [[],[],[],[],[]]
        daythet, jdthet, thet, sthet = [[],[],[],[],[]], [[],[],[],[],[]], [[],[],[],[],[]], [[],[],[],[],[]]
        dayp_filt, jdp_filt, p_filt, s_filt = [[],[],[],[],[]], [[],[],[],[],[]], [[],[],[],[],[]], [[],[],[],[],[]]
        daythet_filt, jdthet_filt, thet_filt, sthet_filt = [[],[],[],[],[]], [[],[],[],[],[]], [[],[],[],[],[]], [[],[],[],[],[]]
        image = []
        limJD=[10000000000000.,0]

        # Getting the values to fit and plot
        for line in lines:

            if line[16] != 'E':
                idx = filters.index(line[3])

                # Refresh the limits
                if float(line[0]) < limJD[0]:
                    limJD[0] = float(line[0])
                if float(line[0]) > limJD[1]:
                    limJD[1] = float(line[0])

                # Not filtered data
                if not any(sub in line[17] for sub in vfilter): # if sub != 'no-std'):
                    jdp[idx] += [float(line[0])]
                    dayp[idx] += [line[1]]
                    p[idx] += [float(line[8])]
                    s[idx] += [float(line[12])]
                    if 'no-std' not in line[17]:
                        jdthet[idx] += [float(line[0])]
                        daythet[idx] += [line[1]]
                        thet[idx] += [float(line[11])]
                        sthet[idx] += [_np.sqrt(float(line[7])**2 + float(line[13])**2)]

                # Filtered data
                elif 'no-std' not in line[17]:
                    jdp_filt[idx] += [float(line[0])]
                    dayp_filt[idx] += [line[1]]
                    p_filt[idx] += [float(line[8])]
                    s_filt[idx] += [float(line[12])]
                    jdthet_filt[idx] += [float(line[0])]
                    daythet_filt[idx] += [line[1]]
                    thet_filt[idx] += [float(line[11])]
                    sthet_filt[idx] += [_np.sqrt(float(line[7])**2 + float(line[13])**2)]

        jdpbpi, pbpi, spbpi = polColor(dayp, jdp, p, s, dayp_filt, jdp_filt, p_filt, s_filt)

#        print jdp, p, s
#        print jdp_filt, p_filt, s_filt


        ##############
        ### 1ST GRAPH
        idx = filters.index(filt)

        # Plot data
        if jdp[idx] != []:
            if limJD[0] != limJD[1]:
                col = _plt.cm.gist_rainbow([(jdi-limJD[0])/(limJD[1]-limJD[0]) for jdi in jdp[idx]])
                image += [axs[0].scatter(jdp[idx], p[idx], marker='o', c=jdp[idx], vmin=limJD[0], \
                                    vmax=limJD[1], edgecolors='white', s=60, cmap=cm)]
            else:
                col = _plt.cm.gist_rainbow([0.5])
                lixo = [axs[0].scatter(jdp[idx], p[idx], marker='o', c=col, \
                                    edgecolors='white', s=60, cmap=cm)]

            for i in range(len(jdp[idx])):
                axs[0].errorbar(jdp[idx][i], p[idx][i], yerr=s[idx][i], linestyle='', \
                                       elinewidth=0.6, marker='', c=col[i], alpha=0.7)

        # Plot ignored data
        if jdp_filt[idx] != []:
            if limJD[0] != limJD[1]:
                col = _plt.cm.gist_rainbow([(jdi-limJD[0])/(limJD[1]-limJD[0]) for jdi in jdp_filt[idx]])
                image += [axs[0].scatter(jdp_filt[idx], p_filt[idx], marker='x', c=jdp_filt[idx], vmin=limJD[0], \
                                    vmax=limJD[1], s=50, linewidths=1.8, cmap=cm)]
            else:
                col = _plt.cm.gist_rainbow([0.5])
                lixo = [axs[0].scatter(jdp_filt[idx], p_filt[idx], marker='x', c=col, \
                                                    s=50, linewidths=1.8, cmap=cm)]

            for i in range(len(jdp_filt[idx])):
                axs[0].errorbar(jdp_filt[idx][i], p_filt[idx][i], yerr=s_filt[idx][i], linestyle='', \
                                    elinewidth=0.6, marker='', color=col[i], alpha=0.7)

        ##############
        ### 2ND GRAPH

        # Plot data 
        if jdpbpi != []:
            if limJD[0] != limJD[1]:
                col = _plt.cm.gist_rainbow([(jdi-limJD[0])/(limJD[1]-limJD[0]) for jdi in jdpbpi])
                image += [axs[1].scatter(jdpbpi, pbpi, marker='o', c=jdpbpi, vmin=limJD[0], \
                                    vmax=limJD[1], edgecolors='white', s=60, cmap=cm)]
            else:
                col = _plt.cm.gist_rainbow([0.5])
                lixo = [axs[1].scatter(jdpbpi, pbpi, marker='o', c=col, \
                                     edgecolors='white', s=60, cmap=cm)]

            for i in range(len(jdpbpi)):
                axs[1].errorbar(jdpbpi[i], pbpi[i], yerr=spbpi[i], linestyle='', \
                                    elinewidth=0.6, marker='', color=col[i], alpha=0.7)

        # Plot ignored data
#        image += [axs[1].scatter(jdp_filt[idx], p_filt[idx], marker='x', c=jdp_filt[idx], vmin=limJD[0], \
 #                                   vmax=limJD[1], s=50, linewidths=1.8, cmap=cm)]
 #       axs[1].errorbar(jdp_filt[idx], p_filt[idx], yerr=s_filt[idx], linestyle='', \
 #                                   elinewidth=0.6, marker='', color='black', alpha=0.7)

        ##############
        ### 3RD GRAPH
        idx = filters.index(filt)

        # Plot data 
        if jdthet[idx] != []:
            if limJD[0] != limJD[1]:
                col = _plt.cm.gist_rainbow([(jdi-limJD[0])/(limJD[1]-limJD[0]) for jdi in jdthet[idx]])
                image += [axs[2].scatter(jdthet[idx], thet[idx], marker='o', c=jdthet[idx], vmin=limJD[0], \
                                    vmax=limJD[1], edgecolors='white', s=60, cmap=cm)]
            else:
                col = _plt.cm.gist_rainbow([0.5])
                lixo = [axs[2].scatter(jdthet[idx], thet[idx], marker='o', c=col, \
                                       edgecolors='white', s=60, cmap=cm)]


            for i in range(len(jdthet[idx])):
                axs[2].errorbar(jdthet[idx][i], thet[idx][i], yerr=sthet[idx][i], linestyle='', \
                                    elinewidth=0.6, marker='', color=col[i], alpha=0.7)

        # Plot ignored data
        if jdthet_filt[idx] != []:
            if limJD[0] != limJD[1]:
                col = _plt.cm.gist_rainbow([(jdi-limJD[0])/(limJD[1]-limJD[0]) for jdi in jdthet_filt[idx]])
                image += [axs[2].scatter(jdthet_filt[idx], thet_filt[idx], marker='x', c=jdthet_filt[idx], vmin=limJD[0], \
                                    vmax=limJD[1], s=50, linewidths=1.8, cmap=cm)]
            else:
                col = _plt.cm.gist_rainbow([0.5])
                lixo = [axs[2].scatter(jdthet_filt[idx], thet_filt[idx], marker='x', c=col, \
                                       s=50, linewidths=1.8, cmap=cm)]
                            
            for i in range(len(jdthet_filt[idx])):
                axs[2].errorbar(jdthet_filt[idx][i], thet_filt[idx][i], yerr=sthet_filt[idx][i], linestyle='', \
                                    elinewidth=0.6, marker='', color=col[i], alpha=0.7)

        # Fix limits
#        ax.autoscale(False)
#        ax.plot(ax.get_xlim(), [0,0], 'k--')
#        ax.plot([0,0], ax.get_ylim(), 'k--')

        # Setting sizes
        axs[2].xaxis.label.set_fontsize(fonts[1]*factor)
        axs[0].yaxis.label.set_fontsize(fonts[1]*factor)
        axs[1].yaxis.label.set_fontsize(fonts[1]*factor)
        axs[2].yaxis.label.set_fontsize(fonts[1]*factor)

        for item in axs[2].get_xticklabels():
            item.set_fontsize(fonts[2]*factor)
        for item in axs[0].get_yticklabels():
            item.set_fontsize(fonts[2]*factor)
        for item in axs[1].get_yticklabels():
            item.set_fontsize(fonts[2]*factor)
        for item in axs[2].get_yticklabels():
            item.set_fontsize(fonts[2]*factor)

        return image
        ###########

   

    _plt.close('all')
    star = _phc.trimpathname(logfile)[1].split('.')[0].split('_')[0]
    if star in _phc.bes:
        be = _phc.bes[star]
    else:
        be = star
    images = []

    if path2 == None or path2 == '.':
        path2 = _os.getcwd()

    # Verify if vfilter is a special filter
    if vfilter in vfil.keys():
        vfilter = vfil[vfilter]
    elif type(vfilter) != list:
        vfilter = []

    # Generate the four axes (sorted as BVRI)
    fig = _plt.figure(1)
    axs = [_plt.subplot(3, 1, 1)]
    axs += [_plt.subplot(3, 1, 2, sharex=axs[0])]
    axs += [_plt.subplot(3, 1, 3, sharex=axs[0])]
    
    # Fix the spacing among the subgraphs and set the QU limits
    _plt.subplots_adjust(hspace=0.06, wspace=0.06)

    # Do the graphs
    images += [plot(fig, axs)]

    # Unset the ticklabels
    _plt.setp(axs[0].get_xticklabels(), visible=False)
    _plt.setp(axs[1].get_xticklabels(), visible=False)
    axs[0].set_xlabel('')
    axs[1].set_xlabel('')
    fig.subplots_adjust(right=0.8)

    # Plot colormap
    if images != [[]]:
        cax = fig.add_axes([0.85, 0.3, 0.02, 0.5])
        cb = _plt.colorbar(images[0][0], cax=cax, orientation='vertical')
        cb.set_label('MJD')
#        cb.ColorbarBase(cax, orientation='vertical', cmap=_plt.cm.gist_rainbow)
#        cb.set_ticklabels(range(int(limJD[0]),int(limJD[1]),50))

    if save:
        _plt.savefig('{0}/{1}_t.{2}'.format(path2,star,extens), bbox_inches='tight')
    else:
        _plt.show(block=False)

    return




#################################################
#################################################
#################################################
def graf_qu(logfile, path2=None, mode=1, thetfile=None, isp=[], odr=True, mcmc=False, \
             nn=[120, 200, 600], thet_ran=[0., 180.], b_ran=[-1., 1.], Pb_ran=[0., 1.], \
             Yb_ran=[-1., 1.], Vb_ran=[0., 1.], clip=True, sclip=4.5, nmax=5, \
             vfilter=['no-std'], save=False, extens='pdf', limQ=None, limU=None, limJD=None):
    """
    Plot a QU diagram for the Be star in the logfile .log
    file (the outfile from polt.genTarget) and fit a line
    if specified. Propagates error of standard star.

         ODR: dot-line
        MCMC: continuous line
        
    mode=1 plot BVRI graphs in the same figure + U in another;
    mode=2 plot UBVRI graphs in separated figures.

    INPUT

        logfile: Logfile with QU data
          path2: Path to save the output graphs
           mode: 1) Plot a figure to filter U and another for
                 BVRI filters; 2) Plot a figure for filter
        thetfile: thet_int.csv file (out from fs.genInt) to
                 to plot the lines using the values inside.
                 In this case, mcmc variable don`t matters in
                 the running.
            isp: interstellar polarization to plot direction
                 in QU diagram
            odr: Run phc.fit_linear to fit a line?
           mcmc: Run fitMCMCline to fit a line?
             nn: fitMCMCline: [n_walkers, n_burnin, n_mcmc]
       thet_ran: fitMCMCline: [thet_min, thet_max]
          b_ran: fitMCMCline: [b_min, b_max]
         Pb_ran: fitMCMCline: [Pb_min, Pb_max], with Pb_min >= 0
         Yb_ran: fitMCMCline: [Yb_min, Yb_max]
         Vb_ran: fitMCMCline: [Vb_min, Vb_max], with Vb_min >= 0
           clip: phc.fit_linear: apply sigma clipping?
          sclip: phc.fit_linear: sigma value to clip
           nmax: phc.fit_linear: sigma clipping max number of
                 iterations
        vfilter: list of flags to filter (will be marked with
                 a 'x' symbol and won't considered in odr
                 fitting). The observations with flag 'no-std'
                 never are shown. mcmc fitting uses the filtered
                 observations, except those with 'no-std' flag.
           save: Save the graphs? If False, just shows
         extens: Extension for the graphs

    OUTPUT

        1) None if odr and mcmc are False
        2) When mcmc==True:
           [[[param_u], [sparam_+_u], [sparam_-_u], n_u, obj_u],
             ...
            [[param_i], [sparam_+_i], [sparam_-_i], n_i, obj_i]],

           where param, sparam_+ and sparam_- are arrays with
           the five parameters for MCMC fitting (thet, b, Pb,
           Yb and Vb) and its errors (at right and left), n is
           the number of points and obj is one graphical dummy object.
        3) When odr==True AND mcmc==False:
           Idem, but param, sparam_+ and sparam_- are arrays with the
           two parameters for ODR fitting (a, b) and its errors.


    CAUTION:

      - The angle returned IS NOT the PA angle obtained from the slop,
        but the own inclination angle of the fitted line! PA is this
        angle divided by 2.
      - PA angle is indefined by a factor of +-90 degrees
      - This routine NEVER shows the data with 'no-std' tag,
        independently of vfilter variable!

    """

    def plotQU(filt, fig, ax, vfilter, odr_fit, mcmc_fit, limq=None, limu=None, limjd=None):
        """
        Receive figure and axes objects and do the plot for filter
        filt WITHOUT show or save the image.

        limq=[qmin, qmax] and limu=[umin, umax] are lists for the min
        and max limits for the axes. Case no specified, they are
        automatically chosen.
        
        Return [param, sparam_+, sparam_-, n, image]

        param, sparam_+ and sparam_- are lists when two peaks are selected
        from chi2 distribution.

        """

        # Setting the color map
        cm = _plt.cm.gist_rainbow
        
        # Factor to fix the font sizes
        if mode==1:
            factor=0.7
        else:
            factor=1.
            
        try:
            lines = _np.loadtxt(logfile, dtype=str)
        except:
            print('# ERROR: Can\'t read file {0}.'.format(logfile))
            raise SystemExit(1)

#        ax.set_title('{0} filter'.format(filt.upper()), fontsize=fonts[0]*factor, verticalalignment='bottom')
        ax.text(0.98, 0.9, '{0} filter'.format(filt.upper()), horizontalalignment='right', \
                 verticalalignment='bottom', transform=ax.transAxes, fontsize=fonts[1]*factor)
        ax.set_xlabel(r'Q (%)', size=fonts[1]*factor)
        ax.set_ylabel(r'U (%)', size=fonts[1]*factor)

        if type(lines[0]) != _np.ndarray and _np.size(lines) == 18:
            lines = lines.reshape(-1,18)

        JD, p, q, u, s, thet, sdth = [],[],[],[],[],[],[]
        JD_filt, p_filt, q_filt, u_filt, s_filt, thet_filt, sdth_filt = [],[],[],[],[],[],[]
        sq, su, sq_filt, su_filt = [],[],[],[]
        image = []

        # Getting the values of the points and filtered points
        for line in lines:
            if line[3] == filt and line[16] != 'E' and not any(sub in line[17] for sub in vfilter):
                JD += [float(line[0])]
                p += [float(line[8])]
                q += [float(line[9])]
                u += [float(line[10])]
                thet += [float(line[11])]
                s += [float(line[12])]
                sdth += [float(line[7])]
            elif line[3] == filt and line[16] != 'E' and 'no-std' not in line[17]:
                JD_filt += [float(line[0])]
                p_filt += [float(line[8])]
                q_filt += [float(line[9])]
                u_filt += [float(line[10])]
                thet_filt += [float(line[11])]
                s_filt += [float(line[12])]
                sdth_filt += [float(line[7])]

        # Propagate errors
        lixo, sq, su = propQU(p, thet, s, sdth)
        lixo, sq_filt, su_filt = propQU(p_filt, thet_filt, s_filt, sdth_filt)
#        print 'original:'
#        print q
#        print q_filt

        # Case some valid data was found
        if [q, u] != [[],[]]:

            # Fitting and plotting by least squares (must be found at least 2 points)
            if odr_fit and len(q) > 1:
                print('='*60)
                print('='*6 + '  FILTER ' + filt.upper())
                print('='*60)
                print('')
                tht, stht, param, sparam1, num = fitodr(q,u,sq,su,JD,q_filt,u_filt,sq_filt,su_filt,JD_filt,filt)
                sparam2 = sparam1
                delt = (max(q+q_filt)-min(q+q_filt))/8
                xadj = _np.linspace(min(q+q_filt)-delt,max(q+q_filt)+delt,3)
                yadj = param[0]*xadj+param[1]
                ax.plot(xadj, yadj, ':', color='dimgray', linewidth=1.7*factor, label='odr')

            # Fitting and plotting by MCMC
            if mcmc_fit or thetfile != None:
                if thetfile != None:
                    param, sparam1, sparam2, n = [],[],[], 0
                    if _os.path.exists(thetfile):
                        fr = open(thetfile, 'r')
                        csvread = csv.reader(fr, delimiter=';')#, quoting=csv.QUOTE_NONE, quotechar='')
                        for i, line in enumerate(csvread):
                            if line[0] == star:
                                # These variables below contain informations about the column numbers where the
                                # filter 'filt' begins inside thetfile file
                                posfilt = 2 + filters.index(filt)*4
                                posfilt2 = 26 + filters.index(filt)*3

                                param = [[float(line[posfilt])*2,float(line[posfilt2]),0,0,0]]
                                sparam1 = [[float(line[posfilt+1])*2,float(line[posfilt2+1]),0,0,0]]
                                sparam2 = [[float(line[posfilt+2])*2,float(line[posfilt2+2]),0,0,0]]

                                # Case there exists a second solution
#                                print line[41:]
                                if filt in line[41:]:
                                    posfilt = 41 + line[41:].index(filt)
                                    param = [param[0]] + [[float(line[posfilt+1])*2,float(line[posfilt+4]),0,0,0]]
                                    sparam1 = [sparam1[0]] + [[float(line[posfilt+2])*2,float(line[posfilt+5]),0,0,0]]
                                    sparam2 = [sparam2[0]] + [[float(line[posfilt+3])*2,float(line[posfilt+6]),0,0,0]]
                        if param == []:
                            print('# WARNING: No star named {0} found inside thetfile file.  The lines won`t be plotted.'.format(be))
                    else:
                        print('# WARNING: No thetfile file found. The lines won`t be plotted.')
                    
                else:
                    if not odr_fit or len(q) <= 1:
                        print('='*60)
                        print('='*6 + '  FILTER ' + filt.upper())
                        print('='*60)
                        print('')
                    param, sparam1, sparam2 = fitmcmc(q+q_filt,u+u_filt,sq+sq_filt,su+su_filt,filt)

                num=len(q+q_filt)
#                print q, q_filt
#                print q+q_filt
#                print param
                # Only plot the curve when the number of points is > 1
                if len(q+q_filt) > 1:
#                    print 'new', param
                    delt = (max(q+q_filt)-min(q+q_filt))/8
                    xadj = _np.linspace(min(q+q_filt)-delt,max(q+q_filt)+delt,3)
                    for i,parami in enumerate(param):
#                        print parami
                        b0 = parami[1]/_np.cos(parami[0]*_np.pi/180)
                        a0 = _np.tan(parami[0]*_np.pi/180)
                        yadj = a0*xadj+b0
#                        print xadj, yadj
                        if i==0:
                            ax.plot(xadj, yadj, '-', color='dimgray', linewidth=1.7*factor, label='mcmc')
                        else:
                            ax.plot(xadj, yadj, '-.', color='dimgray', linewidth=1.7*factor, label='mcmc')
                # Reshape the lists when there is just one peak selected from mcmc.
                if len(param) == 1:
                    param = param[0]
                    sparam1 = sparam1[0]
                    sparam2 = sparam2[0]
            elif not odr or len(q) == 1:
                param, sparam1, sparam2 = [], [], []
                num = 0


            # Specify the colors according to the JD if limjd is None
            if limjd == None or limjd == []:
                limjd = [min(JD+JD_filt), max(JD+JD_filt)]

            # Plot data 
            if len(q) > 0:
#                image += [ax.scatter(pts[0], pts[1], marker='o', edgecolors=colors, facecolors=colors, s=50)]
                if limjd[0] != limjd[1]:
                    col = _plt.cm.gist_rainbow([(jdi-limjd[0])/(limjd[1]-limjd[0]) for jdi in JD])
                    image += [ax.scatter(q, u, marker='o', c=JD, vmin=limjd[0], \
                                            vmax=limjd[1], edgecolors='white', s=72*factor, cmap=cm)]
                else:
                    col = _plt.cm.gist_rainbow([0.5])
                    lixo = ax.scatter(q, u, marker='o', c=col, \
                                      edgecolors='white', s=72*factor, cmap=cm)

                for i in range(len(q)):
                    ax.errorbar(q[i], u[i], xerr=sq[i], yerr=su[i], linestyle='', \
                                                elinewidth=1.05*factor, capsize=4*factor, marker='', c=col[i], alpha=0.7)

            # Plot ignored data
            if len(q_filt) > 0:
#                image += [ax.scatter(pts_filt[0], pts_filt[1], marker='x', edgecolors=colors_filt, facecolors=colors_filt, s=60, linewidths=2)]
#                print "Entrou IF"
                if limjd[0] != limjd[1]:
                    col = _plt.cm.gist_rainbow([(jdi-limjd[0])/(limjd[1]-limjd[0]) for jdi in JD_filt])
                    image += [ax.scatter(q_filt, u_filt, marker='x', c=JD_filt, \
                                vmin=limjd[0], vmax=limjd[1], s=72*factor, linewidths=1.4, cmap=cm)]
                else:
                    col = _plt.cm.gist_rainbow([0.5])
                    lixo = ax.scatter(q_filt, u_filt, marker='x', c=col, \
                                      s=72*factor, linewidths=1.4, cmap=cm)

                for i in range(len(q_filt)):
                    ax.errorbar(q_filt[i], u_filt[i], xerr=sq_filt[i], yerr=su_filt[i], \
                                        linestyle='', elinewidth=1.05*factor, capsize=4*factor, marker='', c=col[i], alpha=0.7)

            # If mode==2, print the colorbar here
            if mode == 2:
                if image != []:
                    fig.subplots_adjust(right=0.8)
                    cax = fig.add_axes([0.85, 0.3, 0.02, 0.5])
                    cb = _plt.colorbar(image[0], cax=cax, orientation='vertical')
                    cb.set_label('MJD')

        else:
            param, sparam1, sparam2 = [], [], []
            num = 0
            image = []


        # Fix limits
        ax.autoscale(False)
        if limq != None and limq != []:
            ax.set_xlim(limq)
        if limu != None and limu != []:
            ax.set_ylim(limu)
        ax.plot(ax.get_xlim(), [0,0], 'k--')
        ax.plot([0,0], ax.get_ylim(), 'k--')

        # plot the diretions of ISP
        for ispi in isp:
#            delt = (max(q+q_filt)-min(q+q_filt))/8
            if ((2*ispi > 270 and 2*ispi <= 360) or (2*ispi< 90 and 2*ispi>= 0)):
                if ax.get_xlim()[1]>0:
                    xisp = _np.array([0,ax.get_xlim()[1]])
                else:
                    print('# WARNING: Angle {0} degree don`t displayed inside the graph area.'.format(ispi))
                    continue
            elif (2*ispi > 90 and 2*ispi < 270):
                if ax.get_xlim()[0]<0:
                    xisp = _np.array([ax.get_xlim()[0],0])
                else:
                    print('# WARNING: Angle {0} degree don`t displayed inside the graph area.'.format(ispi))
                    continue
            else:
                print('# WARNING: This routine can`t plot the ISP on an angle of {0} degree.'.format(ispi))
                continue
            yisp = _np.tan(2*ispi*_np.pi/180)*xisp
            ax.plot(xisp, yisp, ':', color='red', linewidth=1.7*factor)


        # Setting sizes
        ax.xaxis.label.set_fontsize(fonts[1]*factor)
        ax.yaxis.label.set_fontsize(fonts[1]*factor)
        for item in (ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(fonts[2]*factor)


        return [param, sparam1, sparam2, num, image]



    def fitodr(q,u,sq,su,JD,q_filt,u_filt,sq_filt,su_filt,JD_filt,filt):
        """
        Fit ODR
        """

        jd, jd_filt, pts, spts, pts_filt, spts_filt = [],[],[[],[]],[[],[]],[[],[]],[[],[]]
        jd_filt = JD_filt[:]
        pts_filt[0] = q_filt[:]
        pts_filt[1] = u_filt[:]
        spts_filt[0] = sq_filt[:]
        spts_filt[1] = su_filt[:]

        # Fit the simple least squares (only considering y errors) to find
        # initial parameters to the next adjust
        param0, cov = _curve_fit(lambda x,a,b: a*x + b, q, u, sigma=su)
        sparam0 = [_np.sqrt(cov[0][0]), _np.sqrt(cov[1][1])]
        tht0 = _np.arctan(param0[0])*180/_np.pi
        stht0 = (180*sparam0[0])/(_np.pi*(param0[0]**2+1))

        # Fit by the total least squares method (orthogonal distance regression) with clipping
        param, sparam, cov, chi2, niter,bolfilt = _phc.fit_linear(q, u, sq, su, param0=param0,
                                                     clip=clip, sclip=sclip, nmax=nmax)
        tht = _np.arctan(param[0])*180/_np.pi
        stht = (180*sparam[0])/(_np.pi*(param[0]**2+1))

        # Splitting the data into the selected data and the filtered by the clipping
        for i in range(len(q)):
            if bolfilt[i] == 1:
                jd += [JD[i]]
                pts[0] += [q[i]]
                pts[1] += [u[i]]
                spts[0] += [sq[i]]
                spts[1] += [su[i]]
            else:
                jd_filt += [JD[i]]
                pts_filt[0] += [q[i]]
                pts_filt[1] += [u[i]]
                spts_filt[0] += [sq[i]]
                spts_filt[1] += [su[i]]

        # Calculate the reduced chi-squared
        if len(pts[0]) > 2:
            rchi2 = chi2[0]/(len(pts[0])-2)
        else:
            rchi2 = 0

        # Print informations only if mcmc==False (to prevent to print twice)
        if not mcmc:
            print(55*'-')
            print('  Total least squares fit  (y = a*x+b):')
            print(55*'-')
            print('             a = {0:.3f} +- {1:.3f}'.format(param[0], sparam[0]))
            print('             b = {0:.3f} +- {1:.3f}'.format(param[1], sparam[1]))
            print('         theta = {0:.2f} +- {1:.2f} (+- n*90)'.format(tht, stht))
            print('             N = {0:d}'.format(len(pts[0])))
            print('')
            print('             red chi^2 = {0:2f}'.format(rchi2))
            print(55*'-')
            print('')
            
        return tht, stht, param, sparam, len(pts[0])


    def fitmcmc(q, u, sq, su, filt):
        """
        Fit MCMC
        The returned variable are "lists"
        """

        dictvar = ['theta','b', 'P_b', 'Y_b', 'V_b']
        ranges = [thet_ran, b_ran, Pb_ran, Yb_ran, Vb_ran]

        opt2 = ''
        while True:
            thet_mcmc, b_mcmc, Pb_mcmc, Yb_mcmc, Vb_mcmc, fig1, fig2 = fitMCMCline(q, u, sq, su, \
                                star=star+'_'+filt, plot_adj=True, margin=True, n_burnin=nn[1], \
                                n_mcmc=nn[2], n_walkers=nn[0], thet_ran=ranges[0],   \
                                b_ran=ranges[1], Pb_ran=ranges[2], Yb_ran=ranges[3], \
                                Vb_ran=ranges[4], extens=extens)
                        
            if opt2 in ('y','Y'):
                param = [param[0]] + [[thet_mcmc[0], b_mcmc[0], Pb_mcmc[0], Yb_mcmc[0], Vb_mcmc[0]]]
                sparam1 = [sparam1[0]] + [[thet_mcmc[1], b_mcmc[1], Pb_mcmc[1], Yb_mcmc[1], Vb_mcmc[1]]]
                sparam2 = [sparam2[0]] + [[thet_mcmc[2], b_mcmc[2], Pb_mcmc[2], Yb_mcmc[2], Vb_mcmc[2]]]
            else:
                param = [[thet_mcmc[0], b_mcmc[0], Pb_mcmc[0], Yb_mcmc[0], Vb_mcmc[0]]]
                sparam1 = [[thet_mcmc[1], b_mcmc[1], Pb_mcmc[1], Yb_mcmc[1], Vb_mcmc[1]]]
                sparam2 = [[thet_mcmc[2], b_mcmc[2], Pb_mcmc[2], Yb_mcmc[2], Vb_mcmc[2]]]

            # NEW: Requests if the user want to use another interval
            opt = ''
            while opt not in ('y','Y','n','N'):
                print('Do you want to prune the limits and run again MCMC?')
                opt = raw_input('(y/n): ')

            if opt in ('y','Y'):
                while True:
                    ranges = [[],[],[],[],[]]
                    print('')
                    for i, var in enumerate(dictvar):
                        while True:
                            try:
                                etr = raw_input('{0}: specify in format `{0}_min,{0}_max`: '.format(var))
            #                        p_int = [float(ei)-params_fit[0] for ei in petr.split(',')]
                                ranges[i] = [float(ei) for ei in etr.split(',')]
                                if len(ranges[i]) == 2:
                                    if ranges[i][1] > ranges[i][0]:
                                        break
                                    else:
                                        print('Error: {0}_max must be greather than {0}_min!'.format(var))
                                else:
                                    print('Invalid input!')
                            except:
                                print('Invalid input!')

                    opt = ''
                    while opt not in ('y','Y','n','N'):
                        print('\nIs it correct?')
                        for i, var in enumerate(dictvar):
                            print('              {0}_min,{0}_max: {1},{2}'.format(var, ranges[i][0], ranges[i][1]))

                        opt = raw_input('(y/n): ' )
                    if opt in ('y','Y'):
                        _plt.close(fig1)
                        _plt.close(fig2)
                        break
            else:
                # To precent a 'third' peak
                if opt2 in ('y','Y'):
                    opt2 = 'N'
                while opt2 not in ('y','Y','n','N'):
                    print('Do you want select a second peak?')
                    opt2 = raw_input('(y/n): ')

                _plt.close(fig1)
                _plt.close(fig2)
                if opt2 in ('n','N'):
                    break


        return param, sparam1, sparam2



    def fixLimits():
        """
        Found the min and max Q, U and JD values to set the limits of plot
        (excluding the observations for U filter).
        Return 3 lists: [qmin, qmax], [umin, umax], [JDmin, JDmax]
        Return [],[],[] if there is none observations.
        """

        try:
            lines = _np.loadtxt(logfile, dtype=str)
        except:
            print('# ERROR: Can\'t read file {0}.'.format(logfile))
            raise SystemExit(1)

        if type(lines[0]) != _np.ndarray and _np.size(lines) == 18:
            lines = lines.reshape(-1,18)

        q = [float(line[9]) for line in lines if line[3] != 'u' and line[16] != 'E' and 'no-std' not in line[17]]
#                                                 not any(sub in line[17] for sub in vfilter)]
        u = [float(line[10]) for line in lines if line[3] != 'u' and line[16] != 'E' and 'no-std' not in line[17]]
#                                                 not any(sub in line[17] for sub in vfilter)]
        JD = [float(line[0]) for line in lines if line[3] != 'u' and line[16] != 'E' and 'no-std' not in line[17]]
#                                                 not any(sub in line[17] for sub in vfilter)]

        if q==[]:
            return [],[],[]
        
        # A scaled value to shift
        deltq = (max(q)-min(q))/8
        deltu = (max(u)-min(u))/8

        return [min(q)-deltq, max(q)+deltq], [min(u)-deltu, max(u)+deltu], [min(JD), max(JD)]



    _plt.close('all')
    star = _phc.trimpathname(logfile)[1].split('.')[0].split('_')[0]
    if star in _phc.bes:
        be = _phc.bes[star]
    else:
        be = star
    arr, images = [],[]

    if path2 == None or path2 == '.':
        path2 = _os.getcwd()


    # Verify if vfilter is a special filter
    if vfilter in vfil.keys():
        vfilter = vfil[vfilter]
    elif type(vfilter) != list:
        vfilter = []

    ######
    ## 1) Mode 1 plots QU diagram for BVRI filters in the same image and for U in another
    if mode==1:

        ### 1.1 Do the graph for U filter
        fig = _plt.figure()
        ax = _plt.subplot(1, 1, 1)
        arr += [plotQU('u', fig, ax, vfilter, odr, mcmc)]
#        _plt.close(fig_aux)
        
        if save:
            fig.savefig('{0}/{1}_qu_u.{2}'.format(path2,star,extens), bbox_inches='tight')
#            _plt.close(fig)
        else:
            fig.show()
        if odr:
            print('\n')
        
        # Generate the four axes (sorted as BVRI)
        fig = _plt.figure()
        axs = [_plt.subplot(2, 2, 1)]
        axs += [_plt.subplot(2, 2, 2, sharey=axs[0])]
        axs += [_plt.subplot(2, 2, 3, sharex=axs[0])]
        axs += [_plt.subplot(2, 2, 4, sharex=axs[1], sharey=axs[2])]

        for ax in axs:
#            ax.locator_params(axis='x', nbins=6)
            xloc = _plt.MaxNLocator(6)
            ax.xaxis.set_major_locator(xloc)
        
        # Fix the spacing among the subgraphs and set the QU limits
        _plt.subplots_adjust(hspace=0.05, wspace=0.05)
        limq, limu, limjd = fixLimits()
        if limQ != None: limq=limQ
        if limU != None: limu=limU
        if limJD != None: limjd=limJD


        ### 1.2 Do the graphs fo BVRI
        nax = 0
        for filt in ('b','v','r','i'):
            arr += [plotQU(filt, fig, axs[nax], vfilter, odr, mcmc, limq=limq, limu=limu, limjd=limjd)]
            nax += 1

            if arr[-1][-1] != []:
                images += [arr[-1][-1]]
            if odr:
                print('\n')

        # Unset the ticklabels
        _plt.setp(axs[0].get_xticklabels(), visible=False)
        _plt.setp(axs[1].get_xticklabels(), visible=False)
        _plt.setp(axs[1].get_yticklabels(), visible=False)
        _plt.setp(axs[3].get_yticklabels(), visible=False)
        axs[0].set_xlabel('')
        axs[1].set_xlabel('')
        axs[1].set_ylabel('')
        axs[3].set_ylabel('')
        fig.subplots_adjust(right=0.8)

        # Plot colormap
        if images != []:
            cax = fig.add_axes([0.85, 0.3, 0.02, 0.5])
            cb = _plt.colorbar(images[0][0], cax=cax, orientation='vertical')
            cb.set_label('MJD')
#        cb.ColorbarBase(cax, orientation='vertical', cmap=_plt.cm.gist_rainbow)
#        cb.set_ticklabels(range(int(limjd[0]),int(limjd[1]),50))

        if save:
            fig.savefig('{0}/{1}_qu.{2}'.format(path2,star,extens), bbox_inches='tight')
#            _plt.close(fig)
        else:
            fig.show()


    ######
    ## 2) Mode 2 plots QU diagram for UBVRI filters in different images
    elif mode==2:
        for filt in ('u','b','v','r','i'):
            fig = _plt.figure()
            ax = _plt.subplot(1, 1, 1)
            arr += [plotQU(filt, fig, ax, vfilter, odr, mcmc, limq=limQ, limu=limU, limjd=limJD)]
#            _plt.close(fig_aux)

            if save:
                fig.savefig('{0}/{1}_qu_{2}.{3}'.format(path2,star,filt,extens), bbox_inches='tight')
#                _plt.close(fig)
            else:
                _plt.show()

    if odr or mcmc or thetfile != None:
        return arr
    else:
        return




def sintLeff(ccdn='ixon', step=5., save=True, extens='pdf'):
    """
    Sintetizes the response curve, considering the
    CCD Quantum Efficience (QE) and filter transmitance
    from OPD and the stellar models from Pickles (1998).
    Interpolations are made by using cubic splines.

    This code DON'T use none curve for sky transmitance!

    Creates two data files:

    leff_stars_[ccdn].dat : table with l_eff calculated
                            for each star type
    leff_[ccdn].dat : table with the parameters for the
                      adjusted cubic function for l_eff(b-v)
                      (or l_eff(u-b) to the filter U). The
                      M-type stars were excluded from the
                      adjusts, because the molecular lines
                      were affecting these curves.

    'ccdn':  CCD to use in QE curve
                 ixon: CON, Frame Transfer
                 ikon: High Sensitive, Frame Transfer
           ikon-14912: High Sensitive, Frame Transfer
                  301: not avaible

        New eventual CCD files must sample the QE, at least,
        inside the range [2800,11000] Angstrom.

    'step': step, in angstrom, used for the integration
            (Simpson's method) to calculate lambda_eff.
            Allowed values are 5, 10, 15, ..., 500.


    FORMULAS:

    The lbd_eff is computed as **the flux-weighted mean
    wavelenght in terms of photons**:

    lbd_eff = \int(lbd*phi(lbd) * d_lbd) / \int(phi(lbd) * d_lbd)


    where  phi(lbd) = lbd * F_lbd * QE(lbd) * T(lbd)
           F_lbd: stellar flux in erg cm-2 s-1 A-1
           QE(lbd): curve for quantum efficience
           T(lbd):  curve for filter transmitance
           d_lbd: infinitesimal element of wavelength
    
    """

    # List the files for standard stars models
    stars = _glob('{0}/stars/uk*.dat'.format(_hdtpath()))

    # Open file with informations about the standard stars models
    dstars = _np.loadtxt('{0}/stars/synphot.dat'.format(_hdtpath()),usecols=[4,6,7], dtype=str)
    lbds = _np.arange(2800.,11000.001,step)

    # Open file with CCD Quantum Efficience (QE)
    try:
        fqe = _np.loadtxt('{0}/refs/QE_{1}.dat'.format(_hdtpath(),ccdn),skiprows=1, dtype=float, unpack=True) # unpack is to get the transposed array
    except:
        print('ERROR: CCD name \'{0}\' not identified!'.format(ccdn))
        return

    if step not in range(5, 505, 5):
        print('ERROR: step value not valid! Put some value among 5, 10, 15, ..., 500')
        return        

    # Interpolate QE
    qe = _interp1d(fqe[0], fqe[1], kind='cubic')

    # Delete the old .dat files
    if _os.path.exists('leff_stars_{0}.dat'.format(ccdn)):
        _os.unlink('leff_stars_{0}.dat'.format(ccdn))
    if _os.path.exists('leff_{0}.dat'.format(ccdn)):
        _os.unlink('leff_{0}.dat'.format(ccdn))

    with open('leff_stars_{0}.dat'.format(ccdn), 'w') as f0:
        f0.write('{0:5s} {1:>6s} {2:>7s} {3:>7s} {4:>10s}\n'.format('#filt','stype','u-b','b-v','leff'))
    with open('leff_{0}.dat'.format(ccdn), 'w') as f0:
        f0.write('# For U filter:      leff = l0 + k1*(u-b) + k2*(u-b)^2 + k3*(u-b)^3\n')
        f0.write('# For BVRI filters:  leff = l0 + k1*(b-v) + k2*(b-v)^2 + k3*(b-v)^3\n#\n')
        f0.write(('{0:5s} {1:>8s} {2:>10s} {3:>8s} {4:>8s} {5:>7s} {6:>7s} {7:>7s}' +\
                  '{8:>7s} {9:>7s}\n').format('#filt','adj_col','l0','k1','k2','k3',\
                  'sl0','sk1','sk2','sk3'))

    for filt in filters:

        # Open file with Filter Transmitance
        ftr = _np.loadtxt('{0}/filters/T{1}_POL.dat'.format(_hdtpath(),filt.upper()),skiprows=1)

        # Interpolate Filter Transmitance
        if filt=='u':
            tr = _interp1d(_np.arange(2800., 11000.001, 50.), ftr, kind='cubic')
        else:
            tr = _interp1d(_np.arange(2800., 11000.001, 100.), ftr, kind='cubic')

        for star in stars:

            # Open file with flux for star model (genfromtxt allows skip lines in both header and footer)
            fspec = _np.genfromtxt(star, usecols=[0,1], unpack=True, skip_header=330, skip_footer=2800)

            # Mount the array for star spectrum
            spec = _np.array([], dtype=float)
            for i in range(len(fspec[0])):
                if fspec[0][i] in lbds:
                    spec = _np.append(spec, [fspec[1][i]])

            # Interpolate the star spectrum
#            spec = _interp1d(fspec[0], fspec[1], kind='cubic')
            stype = star.split('/')[-1][2:-4].upper()

            # read the color index
            for di in dstars:
                if di[0] == stype:
                    bv = float(di[2])
                    ub = float(di[1])-float(di[2])
                    break

            # Convolute all the curves in the response function
            resp = _np.array([], dtype=float)
            for i in range(len(lbds)):
                resp = _np.append(resp, [lbds[i]*spec[i]*tr(lbds[i])*qe(lbds[i])])

            # Integrate reponse*lambda and response to compute lambda_eff
            l_on = _simps(lbds*lbds*resp, lbds)
            l_under = _simps(lbds*resp, lbds)
            leff = l_on/l_under

            line = '{0:5s} {1:>6s} {2:>7.3f} {3:>7.3f} {4:>10.3f}'.format(filt,stype,ub,bv,leff)
            print(line)
            with open('leff_stars_{0}.dat'.format(ccdn), 'a') as f0:
                f0.write(line+'\n')

            if False:
                _plt.figure()
                ax = _plt.axes()
                _plt.xlabel(r'$\lambda\ (\AA)$', size=fonts[1])
                _plt.ylabel(r'Curves', size=fonts[1])

                _plt.plot(lbds,resp/max(resp), '-', c='black', label='Combined')
                _plt.plot(lbds,spec/max(spec), 'r-.', label='{0} spec'.format(stype))
                tr2 = [tr(lbd) for lbd in lbds]
                _plt.plot(lbds,tr2/max(tr2), 'g--', label='Transm')
                qe2 = [qe(lbd) for lbd in lbds]
                _plt.plot(lbds,qe2/max(qe2), 'b--', label='QE')

                _plt.autoscale(False)
                _plt.ylim([-0.1,1.1])
                _plt.plot([leff,leff], [-0.1,1.1], 'k--', label=r'$\lambda_{eff}$')
                _plt.legend(loc='best', prop={'size':fonts[3]})

                # Setting sizes
                ax.xaxis.label.set_fontsize(fonts[1])
                ax.yaxis.label.set_fontsize(fonts[1])
                for item in (ax.get_xticklabels() + ax.get_yticklabels()):
                    item.set_fontsize(fonts[2])

                _plt.savefig('{0}_{1}_{2}.{3}'.format(stype, ccdn, filt, extens), bbox_inches='tight')



    # Once concluded, we need compute lambda_eff as function of u-b and b-v
    print('\n\n\n# GENERAL ADJUST\n')
    leffs = _np.loadtxt('leff_stars_{0}.dat'.format(ccdn), dtype=str, unpack=True)
    for filt in filters:
        
        ub, bv, leff = _np.array([], dtype=float), _np.array([], dtype=float), _np.array([], dtype=float)

        for i in range(len(leffs[0])):
            # Excluding 'M' type because it is problematic, due to the molecular lines
            if leffs[0][i] == filt and leffs[1][i][0] != 'M':
                ub = _np.append(ub, [float(leffs[2][i])])
                bv = _np.append(bv, [float(leffs[3][i])])
                leff = _np.append(leff, [float(leffs[4][i])])

        _plt.figure()
        ax = _plt.axes()
        _plt.xlabel(r'Color', size=fonts[1])
        _plt.ylabel(r'$\lambda_{eff}\ (\AA)$', size=fonts[1])
        _plt.plot(bv, leff, 'o', c='grey', label='B-V')
        _plt.plot(ub, leff, 's', c='blue', label='U-B')

        if filt == 'u':
            color = ub
            colorstr = 'U-B'
        else:
            color = bv
            colorstr = 'B-V'

        # Fit the lambdas/colors
        param, cov = _curve_fit(lambda x,l0,k1,k2,k3: l0 + k1*x + k2*(x**2) + k3*(x**3), color, leff)
        sparam = _np.array([_np.sqrt(cov[0][0]), _np.sqrt(cov[1][1]),_np.sqrt(cov[2][2]), _np.sqrt(cov[3][3])])
        x = _np.linspace(-1,2,100) 
        y = param[0] + param[1]*x + param[2]*(x**2) + param[3]*(x**3)
        _plt.plot(x, y, '--', c='black', label='{0} fit'.format(colorstr))

        _plt.legend(loc='best', prop={'size':fonts[3]})

        # Setting sizes
        ax.xaxis.label.set_fontsize(fonts[1])
        ax.yaxis.label.set_fontsize(fonts[1])
        for item in (ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(fonts[2])

        line = ('{0:5s} {1:>8s} {2:>10.2f} {3:>8.2f} {4:>8.2f} {5:>7.2f} {6:>7.2f} {7:>7.2f}' +\
                '{8:>7.2f} {9:>7.2f}').format(filt,colorstr,param[0],param[1],param[2],param[3],\
                         sparam[0],sparam[1],sparam[2],sparam[3])
        print(line)
        with open('leff_{0}.dat'.format(ccdn), 'a') as f0:
            f0.write(line+'\n')

        if save:
            _plt.savefig('leff_{0}_{1}.{2}'.format(ccdn, filt, extens), bbox_inches='tight')
        else:
            _plt.show()
    return




def lbds(color, filt, ccdn, airmass=1.3, skiperror=False):
    """
    Return the lambda_eff in angstrom for star with
    color index 'color', in filter 'filt', CCD 'ccdn'
    and airmass 'airmass'. The default airmass is 1.3
    for an average value.

    If skiperror==True, tries to use lambda_eff as the value
    from phc.lbds[] in case of missing information for
    'filt' or 'ccd'.

    CAUTION:
    If filt=='u', the color must be U-B
       Otherwise, the color must be B-V

    FORMULAS:
      - Atm. reddening:   redn = 2.5*_np.log10(_np.e)*
                                    *airmass*(taub-tauv)
      - Redd. color:  color_av = color + redn
      - lambda_eff:      l_eff = l0 + k1*color_av +
                                    + k2*(color_av**2) +
                                    + k3*(color_av**3)

       where tauu, taub, tau are the optical deepth of
       atmosphere (values used from Kepler de Oliveira
       et al, Astronomia e Astrofisica).
    
    """
    
    data = _np.loadtxt('{0}/filters/leff.dat'.format(_hdtpath()), dtype=str)

    # Optical deepth according to Kepler de Oliveira et al (Astronomia
    # e Astrofisica), for altitude above 2000m
    tauu=1.36
    taub=0.52
    tauv=0.37

    if filt=='u':
        redn = 2.5*_np.log10(_np.e)*airmass*(tauu-taub)
    else:
        redn = 2.5*_np.log10(_np.e)*airmass*(taub-tauv)

    l0=0
    for line in data:
        if line[0] == filt and line[2] == ccdn:
            l0 = float(line[3])
            k1 = float(line[4])
            k2 = float(line[5])
            k3 = float(line[6])
            break

    if l0==0:
        if skiperror:
            leff = _phc.lbds[filt]
        else:
            print('# ERROR: parameters to calculate lambda_eff in filter {0} and CCD {1} not found.'.format(filt,ccdn))
            raise SystemExit(1)
    else:
        color = color + redn
        leff = l0 + k1*color + k2*(color**2) + k3*(color**3)

    return leff



#################################################
#################################################
#################################################
### MAIN ###
if __name__ == "__main__":
    pass





def fitMCMCline(x, y, sx, sy, star='', margin=False, plot_adj=True, fig=None, ax=None, \
                                            n_burnin=350, n_mcmc=600, \
                                    n_walkers=120, thet_ran=[0., 180.], \
                                    b_ran=[-1., 1.], Pb_ran=[0., 1.], \
                                   Yb_ran=[-1., 1.], Vb_ran=[0., 1.], extens='pdf'):
    """
        Fit a line using Markov Chain Monte Carlo for data
        with both x and y errors and with bad points
        from emcee code.
        
        The model is tha sum of two models:
        a) a line model for the good points, y = ax + b, with
        data displaced ortogonally by a gaussian diplacentment
        (covariance is suposed null!);
        b) a gaussian distribution orthogonal to the line for
        the bad points with amplitude, mean and variance equal
        to Pb, Yb and Vb (see Hoog, Bovy and Lang,
        ``Data analysis recipes: Fitting a model to data'')

        The MCMC does a ``mixture'' among both model for each
        point. So, it is not needed know about the bad and
        good points necessairly! The five parameters to be
        found are theta=arctan(a), b, Pb, Yb and Vb.


      INPUT:
        
      x/y/sx/sy: array/list with the data
           star: star name to be printed in the graph and
                 its filename. If it's a void str '', this
                 routine give a random number to prevent
                 overwriting data.
         margin: marginalize the corner graphs over Pb,
                 Yb, Vb (generating the graph only for
                 theta and b)?
       plot_adj: show a plot of data+fit?
            fig: Figure to append the plots for data+fit.
                 If None, a new figure is generated.
             ax: Axes, as like fig above.
       n_burnin: number of iterations for burning-in
         n_mcmc: number of iterations to run emcee
      n_walkers: number of walkers to map the posterior
                 probabilities.
       thet_ran: [thet_min, thet_max]
          b_ran: [b_min, b_max]
         Pb_ran: [Pb_min, Pb_max], with Pb_min >= 0
         Yb_ran: [Yb_min, Yb_max]
         Vb_ran: [Vb_min, Vb_max], with Vb_min >= 0
         extens: extension for the graph file


      OUTPUT:

         theta_fit: [theta, theta_+, theta_-], the theta value
                    and its errors (at right and left from it).
                    theta is the median of distribution
                    probability and theta_+, theta_- are the
                    range within which there are 68.3% of the
                    points in such distribution.
      b*cos(theta): Idem, for b*cos(theta).
                Pb: Idem, for Pb.
                Yb: Idem, for Yb.
                Vb: Idem, for Vb.
              fig1: Figure pointer to the corner graph ([]
                    if show==False).
              fig2: Figure pointer to the data+fit graph ([]
                    if plot_adj==False).
        

      FORMULAS:

        Supposing i as each data point, the log of Likelihood
        function (L) takes form:

          log(L) = sum(log(p_good_i+p_bad_i))

        with

          p_good_i = (1-Pb)/sqrt(2*pi*var_i)*
                                exp(-0.5*(disp_i**2/var_i)),

           p_bad_i = Pb/sqrt(2*pi*(Vb+var)*
                            exp(-0.5*(disp_i-Yb)**2)/(Vb+var_i))

        where disp_i is the total projection of the (x_i,y_i)
        values over the line and var_i, the projected variance:

          disp_i = v*Z_i -b cos(thet)
           var_i = v*S_i*v

        with
        
            v = (-sin(theta), cos(theta)) (versor orthogonal
                                                 to the line)
          Z_i = (x_i, y_i)    (data point)
          S_i = | sx_i^2  sxy_i^2| = |sx_i^2     0 | = (covariance 
                |syx_i^2   sy_i^2|   | 0     sy_i^2|       matrix)


    """

    import emcee
    import triangle.nov
    from matplotlib.ticker import MaxNLocator


    def lnprob(params, xx, yy, sxx, syy):
        """
        Return the log of posterior probability (p_pos) in
        bayesian statistics for the parameters 'params' and the
        data poits xx, yy, sxx and syy.

        p_pos = L*p_prior (unless by a normalization constant),
        where L is the likelihood function and p_prior is the
        prior probability function.


        a) Likelihood
        
        In our case, for gaussian and independent uncertaities,
        in both x and y axes and with bad points:

        log(L) = sum(log(p_good_i+p_bad_i))

        with

        p_good_i = (1-Pb)/sqrt(2*pi*var_i)*exp(-0.5*(disp_i**2/var_i)),
        p_bad_i = Pb/sqrt(2*pi*(Vb+var)*exp(-0.5*(disp_i-Yb)**2)/(Vb+var_i))

        where disp_i is the total projection of the (x_i,y_i) values
        over the line and var_i, the projected variance; Pb, Yb, Vb are
        the gaussian model for the bad points - the amplitude, mean and
        variance (see Hoog, Bovy and Lang, ``Data analysis recipes:
        Fitting a model to data'')

        Taking the model for the line y = ax + b, where a = tan(thet),
        let

        v = (-sin(thet), cos(thet)) (versor orthogonal to the line)
        Z_i = (x_i, y_i)
        S_i = | sx_i^2   sxy_i^2|  (covariance matrix)
              |syx_i^2    sy_i^2|

        The formulas for disp_i and var_i are:
        
        disp_i = v*Z_i -b cos(thet)
        var_i = v*S_i*v


        b) p_prior
        
        Now, p_prior = constant for 'params' values inside the
        range defined by 'intervalos'; otherwise, it is 0.
        That is the only determination that we can do.

        So, p_pos = log(L) or -inf case 'params' are out from
        the allowed range.
        """

        thet, b, Pb, Yb, Vb = params
        b0 = b/_np.cos(thet*_np.pi/180)

        # Set prior ln prob
        lnprior = 0
        for i, interv in enumerate(intervalos):
            if params[i] < interv[0] or params[i] > interv[1]:
                lnprior = -_np.inf

        # Return posterior prob
        if not _np.isfinite(lnprior):
            return -_np.inf
        else:
            sin = _np.sin(thet*_np.pi/180)
            cos = _np.cos(thet*_np.pi/180)
#            cov = sin*cos*(b0**2)
            disp = -x*sin + y*cos - b
#            var = (1-cov)*(sin*sxx)**2 + (1-1/cov)*(cos*syy)**2
            # Projected variance WITHOUT covariance terms
            var = (sin*sxx)**2 + (cos*syy)**2

            prob_good = (1-Pb)/(_np.sqrt(2*_np.pi*var))*_np.exp(-0.5*(disp**2/var))
            prob_bad = Pb/_np.sqrt(2*_np.pi*(Vb+var))*_np.exp(-0.5*((disp-Yb)**2)/(Vb+var))
#            print '-'*40
#            print prob_good, prob_bad, disp, var#, cov

            for i in range(len(prob_good)):
                if prob_good[i]+prob_bad[i] <= 0:
                    return -_np.inf
            
            return lnprior + sum(_np.log(prob_good + prob_bad))



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
        print("Mean acceptance fraction:", _np.mean(af))

        # Compute the results using all interval
        fig1, thet_mcmc, b_mcmc, Pb_mcmc, Yb_mcmc, Vb_mcmc = gen_results(sampler, intervalos, save=True, show=True)

        # Plot convergence map
        plot_conv(sampler, [thet_mcmc[0], b_mcmc[0], Pb_mcmc[0], Yb_mcmc[0], Vb_mcmc[0]])

        if plot_adj:
            ax2.errorbar(x,y,xerr=sx,yerr=sy, linestyle='', marker='o')

            if len(x) > 1:
                xadj = _np.linspace(min(x),max(x),2)
                yadj = podr[0]*xadj+podr[1]
                ax2.plot(xadj, yadj, '-.', color='dimgray', label='odr')

            xadj = _np.linspace(min(x),max(x),2)
            b0 = b_mcmc[0]/_np.cos(thet_mcmc[0]*_np.pi/180)
            a0 = _np.tan(thet_mcmc[0]*_np.pi/180)
            yadj = a0*xadj+b0
            ax2.plot(xadj, yadj, '-',color='dimgray', label='mcmc')
            ax2.legend(loc='best')
            fig2.show()

        # NEW: Requests if the user want to use some specific interval
        opt = ''
        while opt not in ('y','Y','n','N'):
            print('Do you want to select specific ranges to compute the values?')
            opt = raw_input('(y/n): ')
        if opt in ('y','Y'):

            while True:
                ranges = [[],[],[],[],[]]
                print('')
                for i, var in enumerate(dictvar[0]):
                    while True:
                        try:
                            etr = raw_input('{0}: specify in format `{0}_min,{0}_max`: '.format(var))
    #                        p_int = [float(ei)-params_fit[0] for ei in petr.split(',')]
                            ranges[i] = [float(ei) for ei in etr.split(',')]
                            if len(ranges[i]) == 2:
                                if ranges[i][1] > ranges[i][0]:
                                    break
                                else:
                                    print('Error: {0}_max must be greather than {0}_min!'.format(var))
                            else:
                                print('Invalid input!')
                        except:
                            print('Invalid input!')

                opt = ''
                while opt not in ('y','Y','n','N'):
                    print('\nIs it correct?')
                    for i, var in enumerate(dictvar[0]):
                        print('              {0}_min,{0}_max: {1},{2}'.format(var, ranges[i][0], ranges[i][1]))

                    opt = raw_input('(y/n): ' )
                if opt in ('y','Y'):
                    print('')
                    break

            _plt.close(fig1)
            fig1, thet_mcmc, b_mcmc, Pb_mcmc, Yb_mcmc, Vb_mcmc = gen_results(sampler, ranges, save=True, show=True)
#        else:
#            _plt.close(fig2)


        return fig1, thet_mcmc, b_mcmc, Pb_mcmc, Yb_mcmc, Vb_mcmc, opt



    def gen_results(sampler, ranges, show=True, save=True):

        # Read the sampler
        samples = sampler.chain[:, :, :].reshape((-1, ndim))
        samples_new = _np.empty(shape=[0, ndim])

        # Filtering 'samples' array
        print('Please wait, computing values from samples...')
        new=False
        for i, rang in enumerate(ranges):
            if intervalos[i][0] != rang[0] or intervalos[i][1] != rang[1]:
                new = True
                break
        if new:
            for elem in samples:
                ver = True
                for i in range(ndim):
                    if elem[i] > ranges[i][1] or elem[i] < ranges[i][0]:
                        ver = False
                        break
                if ver:
                    samples_new = _np.vstack([samples_new, elem])
        else:
            samples_new = samples

        # Ploting corner graph
        fig1 = triangle.nov.corner(samples_new, title=star, \
#                            truths=[p_mcmc[0], l_mcmc[0]], \
#                            extents=[(p_range[0],l_range[0]),(p_range[1],l_range[1])], \
                             quantiles=[0.16075, 0.50, 0.83925], \
                             labels=dictvar[1], \
                             verbose=False)
        if save:
            fig1.savefig('{0}_correl.{1}'.format(star,extens))
        if show:
            fig1.show()
        else:
            fig1 = []
        if margin:
            # Ploting corner graph
            fig3 = triangle.nov.corner(samples_new[:,0:2], title=star, \
#                            truths=[p_mcmc[0], l_mcmc[0]], \
#                            extents=[(p_range[0],l_range[0]),(p_range[1],l_range[1])], \
                             quantiles=[0.16075, 0.50, 0.83925], \
                             labels=dictvar[1], \
                             verbose=False)
            fig3.savefig('{0}_correl_m.{1}'.format(star,extens))
        
        # Computing the medians and errors according to the range from median
        # inside which there are 68.3% of the data
        thet_mcmc, b_mcmc, Pb_mcmc, Yb_mcmc, Vb_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                                                zip(*_np.percentile(samples_new, [16.075, 50, 83.925], axis=0)))

        #~ Print the output
        """ TBD """
        print('')
        print(55*'-')
        print('  2) MCMC best values  (y = tan(theta)*x + b*cos(theta)):')
        print(55*'-')
        print('          theta = {0:9.4f}  +{1:.4f}  -{2:.4f}'.format(thet_mcmc[0],thet_mcmc[1],thet_mcmc[2]))
        print('   b*cos(theta) = {0:9.4f}  +{1:.4f}  -{2:.4f}'.format(b_mcmc[0],b_mcmc[1],b_mcmc[2]))
        print('             Pb = {0:9.4f}  +{1:.4f}  -{2:.4f}'.format(Pb_mcmc[0],Pb_mcmc[1],Pb_mcmc[2]))
        print('             Yb = {0:9.4f}  +{1:.4f}  -{2:.4f}'.format(Yb_mcmc[0],Yb_mcmc[1],Yb_mcmc[2]))
        print('             Vb = {0:9.4f}  +{1:.4f}  -{2:.4f}'.format(Vb_mcmc[0],Vb_mcmc[1],Vb_mcmc[2]))
#        print 'reduced chi2 = {0:.4f}'.format(chi)
        print(55*'-')
        print('')


        return fig1, thet_mcmc, b_mcmc, Pb_mcmc, Yb_mcmc, Vb_mcmc


    def plot_conv(sampler, param):
        """
        Plot convergence map. 'param' are the values to be highlighted
        """

        fig4, axes = _plt.subplots(ndim, 1, sharex=True, figsize=(8, 15))

        for i in range(ndim):
            axes[i].plot(sampler.chain[:, :, i].T, color="k", alpha=0.4)
            axes[i].yaxis.set_major_locator(MaxNLocator(5))
            axes[i].axhline(param[i], color="#888888", lw=2)
            axes[i].set_ylabel(dictvar[1][i])

        axes[4].set_xlabel("Step number")

        fig4.tight_layout(h_pad=0.0)
        fig4.savefig('{0}_conv.{1}'.format(star,extens))

        return



    def fitodr():
        """
        Fit by Least Squares
        """

        # Fit the simple least squares (only considering y errors) to find
        # initial parameters to the next adjust
        param0, cov = _curve_fit(lambda t,a,b: a*t + b, x, y, sigma=sy)
#        sparam0 = [_np.sqrt(cov[0][0]), _np.sqrt(cov[1][1])]
#        tht0 = _np.arctan(param0[0])*180/_np.pi
#        stht0 = (180*sparam0[0])/(_np.pi*(param0[0]**2+1))

        # Fit by the total least squares method (orthogonal distance regression) with clipping
        param, sparam, cov, chi2, niter,bolfilt = _phc.fit_linear(x, y, sx, sy, param0=param0,
                                                                        clip=False)
        tht = _np.arctan(param[0])*180/_np.pi
        if tht < 0:
            tht += 180.
        elif tht >= 180:
            tht -= 180
        stht = (180*sparam[0])/(_np.pi*(param[0]**2+1))
        nn = sum(bolfilt)

        # Calculate the reduced chi-squared
        if nn > 2:
            rchi2 = chi2[0]/(nn-2)
        else:
            rchi2 = 0

        # Print informations
        print(55*'-')
        print('  1) Total least squares fit  (y = a*x+b):')
        print(55*'-')
        print('             a = {0:.3f} +- {1:.3f}'.format(param[0], sparam[0]))
        print('             b = {0:.3f} +- {1:.3f}'.format(param[1], sparam[1]))
        print('         theta = {0:.2f} +- {1:.2f}'.format(tht, stht))
        print('             N = {0:d}'.format(nn))
        print('')
        print('             red chi^2 = {0:2f}'.format(rchi2))
        print(55*'-')
        print('')

        return param, sparam




    # Dictionary for the graph labels
    dictvar =      [['theta',
              'b*cos(theta)',
                       'P_b',
                       'Y_b',
                       'V_b'],
                  [r'$\theta$',
          r'$b\,\cos \theta $',
                      r'$P_b$',
                      r'$Y_b$',
                      r'$V_b$']]


#    try:
#        _plt.close(fig2)
#    except:
#        pass
#    x = _np.array([201., 244., 47., 287., 203., 58., 210., 202., 198., 158., 165., 201., 157., 131., 166., 160., 186., 125., 218., 146.])
#    y = _np.array([592., 401., 583., 402., 495., 173., 479., 504., 510., 416., 393., 442., 317., 311., 400., 337., 423., 334., 533., 344.])
#    sx = _np.array([9.,4.,11.,7.,5.,9.,4.,4.,11.,7.,5.,5.,5.,6.,6.,5.,9.,8.,6.,5.])
#    sy = _np.array([61.,25.,38.,15.,21.,15.,27.,14.,30.,16.,14.,25.,52.,16.,34.,31.,42.,26.,16.,22.])

    # thet, b, Pb, Yb, Vb

    # Setting parameters and limits
    intervalos = _np.array([thet_ran, b_ran, Pb_ran, Yb_ran, Vb_ran])
    ndim = 5

    # Converting lists to np.array
    if type(x) == list:
        x = _np.array(x)
    if type(y) == list:
        y = _np.array(y)
    if type(sx) == list:
        sx = _np.array(sx)
    if type(sy) == list:
        sy = _np.array(sy)

    # Fit by Least Squares
    if len(x) > 1:
        podr, spodr = fitodr()
    if plot_adj:
        if ax==None or fig==None:
            fig2 = _plt.figure()
            ax2 = _plt.subplot(1, 1, 1)
        else:
            fig2 = fig
            ax2 = ax

    # If 'star' was not specified, generate a random number to append to the graph name to be saved
    if star == '':
        star = 'rand' + str(int(_np.random.rand(1)[0]*10000))

    # Define random values to be used as priori numbers within the interval
    p0 = _np.array( [_np.random.rand(ndim) for n in xrange(n_walkers)] )
    for k in range(ndim):
        p0[:,k] = intervalos[k][0]+p0[:,k]*(intervalos[k][1]-intervalos[k][0])

    # Initialize the sampler and run mcmc
    sampler = emcee.EnsembleSampler(n_walkers, ndim, lnprob, args=[x, y, sx, sy], a=3)#, threads=2)
    fig1, thet_mcmc, b_mcmc, Pb_mcmc, Yb_mcmc, Vb_mcmc, opt = run_emcee(sampler, p0)


    # Plot only if plot_adj==True or a new computation was done
    if plot_adj and opt in ('y', 'Y'):

        xadj = _np.linspace(min(x), max(x), 2)
        b0 = b_mcmc[0]/_np.cos(thet_mcmc[0]*_np.pi/180)
        a0 = _np.tan(thet_mcmc[0]*_np.pi/180)
        yadj = a0*xadj+b0
        ax2.plot(xadj, yadj, '--', color='dimgray', label='mcmc_new')
        ax2.legend(loc='best')
        fig2.show()
    elif not plot_adj:
        fig2 = []

    return thet_mcmc, b_mcmc, Pb_mcmc, Yb_mcmc, Vb_mcmc, fig1, fig2


def loadpol(txt):
    """ Load polarization txt file. """
    dtb = _np.loadtxt(txt, dtype=str)
    dtb = _np.core.records.fromarrays(dtb.transpose(), names='MJD,night,filt,\
    calc,stdstars,dth,devdth,P,Q,U,th,sigP,sigth', formats='f8,{0},{0},f8,{0},\
    f8,f8,f8,f8,f8,f8,f8,f8'.format(dtb.dtype))
    return dtb


if __name__ == '__main__':
    pass
