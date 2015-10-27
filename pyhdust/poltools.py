#-*- coding:utf-8 -*-

"""
PyHdust *poltools* module: polarimetry tools

History:
-grafpol working for *_WP1110....log files!
-grafpol working for log/out files with more than a single star

:co-author: Daniel Bednarski
:license: GNU GPL v3.0 (https://github.com/danmoser/pyhdust/blob/master/LICENSE)
"""

import os as _os
import re as _re
import pwd as _pwd
import time as _time
import glob as _glob
import numpy as _np
import datetime as _dt
import shutil as _shutil
#from itertools import product as _product
from glob import glob as _glob
from inspect import getouterframes as _getouterframes
from inspect import currentframe as _currentframe
import pyhdust.phc as _phc 
import pyhdust.jdcal as _jdcal
from pyhdust import hdtpath as _hdtpath
#from sys import _argv
#from matplotlib import rc as _rc

try:
    import matplotlib.pyplot as _plt
    from matplotlib.transforms import offset_copy as _offset_copy
    import pyfits as _pyfits
except:
    print('# Warning! matplotlib and/or pyfits module not installed!!!')

__author__ = "Daniel Moser"
__email__ = "dmfaes@gmail.com"


filters = ['u','b','v','r','i']
# Setting an "initial value" for ccd
ccd = '---'

# Dictionary for the tags entered by the user
dictags = {0: ['bad modulation','bad-mod'],
           1: ['very bad modulation','very-bad-mod'],
           2: ['some modulations are incompatible','incomp-mods'],
           3: ['some observational problem/error','obs-prob'],
           4: ['polarimeter problem suspected','iagpol-prob'],
           5: ['another relevant problem','other-prob'],
          }

# Dictionary for the tags assigned automatically
# If you want to add another value, add inside verout routin also.
dictests = {0: ['std incompatible with the published','obs!=pub', 'W'],
            1: ['sig >> theorical_sig','s>>theor_s', 'W'],
            2: ['no standard in the night','no-std', 'W'],
           }



#################################################
#################################################
#################################################
def stdchk(stdname):
    """
    Check if the standard star name contains an known name, and return
    its position in `padroes.txt`.
    """
    lstds = list(_np.loadtxt('{0}/pyhdust/refs/pol_padroes.txt'.format(_hdtpath()), dtype=str,\
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

    """

    if float(MJD) < 57082.5: # before 2015, March 1st
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
    if len(coords) == 0:
        coords = _glob('{0}/coord_*_{1}.ord'.format(path,f))
        if len(coords) == 0:
            coords = _glob('{0}/coord_*_{1}_[0-9]*.ord'.format(path,f))
            if len(coords) == 0:
                print(('# ERROR: Found *_{0}_*.out files, but none COORD file found as '+\
                        '{1}/coord_*_{2}_*.ord. Verify and run again.\n').format(f,path,f))
                raise SystemExit(1)

    if len(coords) != 0:
        try:
            coords = _np.loadtxt(coords[0])
            ang = _np.arctan( (coords[1,1]-coords[0,1])/(coords[1,0]-coords[0,0]) )*180/_np.pi
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
def chooseout(objdir, obj, f, nstar=1, sigtol=lambda sig: 1.2*sig):
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
        menor erro == 0.000200. Se sigtol(sig) = 1.2*sig, como
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

        Return err, out. If not found, return 100.,''.
        """

        err = 100.
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
    err = [100.]*n
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
        errtmp = err[i]   # errtmp==100 if no outfiles were found by minErrBlk16

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
def verout(out, obj, f, nstar=1, verbose=True, delta=2.5):
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
    if sig_ratio > 3.:
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
def queryout(objdir, obj, f, nstar=1, sigtol=lambda sig: 1.2*sig):
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
            print(' {0:<10s} {1:<5s} {2:<7s} {3:<8s} {4:<10s} {5:<7s} {6:<s}'.format('Obj', 'Filt', \
                    'P(%)', 'sig(%)', 'sig/ThSig', 'ztest', 'out/num'))
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
            print(' {0:<10s} {1:<5s} {2:<7.4f} {3:<8.4f} {4:<10.3f} {5:<7s} {6:<s} {7:<s}'.\
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
        print('# Averaged {0} band is {1:.3f} +/- {2:.3f} %'.format(f.upper(),avg,\
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
    # The input variables are exactly the outuput sublogs and groups from combineout subroutine.
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

        # Set the logfiles once, erasing some void components
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
            print('# ERROR: {0} figure was not displayed by grafall'.format(len(logs)-8))
            return

        if   ncol == 1: linit=0.15
        elif ncol == 2: linit=0.08
        else:           linit=0.05

        fig = _plt.figure(figsize=(4*ncol,3.4*nlin))

        # Estabelece os grids e cria todos os eixos
        grids = [ _plt.GridSpec(2*nlin, ncol, hspace=0, wspace=0.35, \
                        top=0.94, bottom=0.19, left=linit, right=0.95) ]
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

        return(Q, U, sigma, P_pts, th_pts, str_pts, nstars)


    def plotlog(ax1, ax2, Q,U,sigma,P_pts,th_pts,str_pts,filename):

        # extract the group number
        WPpos = filename.find('_WP')
        if WPpos == -1:
            suff = filename[filename.rfind('_')+1:]
        else:
            suff = filename[filename.rfind('_', 0, WPpos)+1:]

        ax1.set_title('P(%) = {0:.4f}+-{1:.4f}'.format(_np.sqrt(Q**2+U**2)*100,
                              sigma*100), fontsize=14, verticalalignment='bottom')
        ax1.text(0.98, 0.01, '{0}'.format(suff), horizontalalignment='right', \
                 verticalalignment='bottom', transform=ax1.transAxes, fontsize=9)
        ax1.set_ylabel('P (%)', size=9)
        
        ysigma = _np.zeros(len(th_pts))+sigma
        ax1.errorbar(th_pts,P_pts*100,yerr=ysigma*100)
    
        th_det = _np.linspace(th_pts[0]*.98,th_pts[-1]*1.02,100)
        P_det = Q*_np.cos(4*th_det*_np.pi/180)+U*_np.sin(4*th_det*_np.pi/180)

        ax1.plot(th_det, P_det*100)
        ax1.plot([th_det[0],th_det[-1]], [0,0], 'k--')    
        ax1.set_xlim([th_pts[0]-4,th_pts[-1]*1.02+3])
#        ax1.set_xlim([th_pts[0]-4,th_pts[-1]*1.02+3][::-1])
#        ax1.set_ylim([min(P_pts*100)*1.1, max(P_pts*100)*1.1])

        ax2.set_xlabel('WP position', size=9)
        ax2.set_ylabel('Residuals', size=9)
        ax2.set_xlim([th_pts[0]-4,th_pts[-1]*1.02+3])
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
        Q, U, sigma, P_pts, th_pts, str_pts, nstars = readlog(filename)
        plotlog(ax1,ax2, Q,U,sigma,P_pts,th_pts,str_pts,filename)
        if save:
            if nstars == 1:
                _plt.savefig(filename.replace('.log','.'+extens))
            else:
                _plt.savefig(filename.replace('.log','_star{0}.{1}'.format(nstar,extens)))
        else:
            _plt.show()
    else:
        Q, U, sigma, P_pts, th_pts, str_pts, nstars = readlog(filename)
        plotlog(ax1,ax2, Q,U,sigma,P_pts,th_pts,str_pts,filename)
            
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
    lstds = _np.loadtxt('{0}/pyhdust/refs/pol_padroes.txt'.format(_hdtpath()), dtype=str)

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
    another tags already assigned and current flag. The flag
    returned is the worst flag between that and those found
    (e.g., if input 'flag' is 'W' and tests results on flag
    'OK', return flag 'W'; if tests results on flag 'E',
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
def chkStdLog(f, calc, path=None, delta=2.5, verbose=True):
    """
    Verify if there are standards for filter `f` and
    calcite `calc` inside path/std.dat. Return True if
    successful, unless the standard has been reduced
    and marked with `E` flag.

    delta is the allowed variation for the angles between the two
    beams for one same calcite.
    """

    loglines = ''
    if path == None:
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
def genLog(path, subdirs, tgts, fileout, sigtol=lambda sigm: 1.2*sigm, \
                    autochoose=False, delta=2.5):
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
            The default values is sigtol=lambda sigm: 1.2*sigm, while
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
            if nstars == 0 and len(_glob('{0}/{1}/*_{2}_*.fits'.format(path,objdir,f))) > 0:
                print(('\n# ERROR: {0}_{1}: Fits files found, but the object was not reduced! ' +\
                        'Reduce and run again...\n\n - HINT: if these fits files compose some ' +\
                        'non-valid serie but need be kept in, move them for a subdir {2}/tmp, ' +\
                        'and hence, the path will not be sweept by routine.\n').format(objdir,f,objdir))
                raise SystemExit(1)
            # Check if there exist some .out file for such object/filter, but not the fits files
            elif nstars != 0 and len(_glob('{0}/{1}/*_{2}_*.fits'.format(path,objdir,f))) == 0:
                print(('\n# ERROR: {0}_{1}: Fits files not found, but were found another *_{2}_* files. ' +\
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
def genAllLog(path=None, sigtol=lambda sigm: 1.2*sigm, autochoose=False, delta=2.5):
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
            The default values is sigtol=lambda sigm: 1.2*sigm, while
            the old value was sigtol=lambda sigm: 1.1*sigm + 0.00005. If
            you want to take just the groups with all WP, and none other,
            you can specify sigtol=lambda sigm: 1000.*sigm, for example.
    autochoose: choose best outfiles automatically, without
                interaction?
    """

    if path == None:
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
        ltgts = _np.loadtxt('{0}/pyhdust/refs/pol_alvos.txt'.format(_hdtpath()), dtype=str)
        if _os.path.exists('{0}/pyhdust/refs/pol_hip.txt'.format(_hdtpath())):
            try:
                ltgts = _np.concatenate((ltgts,_np.loadtxt('{0}/pyhdust/refs/pol_hip.txt'.\
                                                    format(_hdtpath()), dtype=str)))
            except:
                pass
        lstds = _np.loadtxt('{0}/pyhdust/refs/pol_padroes.txt'.format(_hdtpath()), \
                                                           dtype=str, usecols=[0])
    except:
        print('# ERROR: Can\'t read files pyhdust/pyhdust/refs/pol_alvos.txt and/or pyhdust/pyhdust/refs/pol_padroes.txt.\n')
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
        while obj not in _np.hstack((ltgts,lstds,_np.array(['calib']))):
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
        elif obj == 'calib':
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
def corObjStd(night, f, calc, path=None, delta=2.5):
    """
    Correlate a target observed at filter 'f' and calcite
    'calc' (expected the angle between ord. and extraord. beams)
    and return the values for matching standard stars inside
    'night'/std.dat file, except by those marked with 'E' flag.

    delta: tolerance, in degree, of angle of two beams for
           one same calcite (default: +/-2.5 degree)

    Output, in order: stdname, thstd, angref, flagstd
      - stdname: list with standard names
      - thstd: list with theta measured values (the angle
               returned is the same from .out, NOT 180-theta!)
      - angref: list with theta published values
      - flastd: list with flags concerning to the target data.
    """
    if path == None:
        path = _os.getcwd()
    try:
        stdref = _np.loadtxt('{0}/pyhdust/refs/pol_padroes.txt'.format(_hdtpath()), dtype=str)
#        f0 = open('{0}/pyhdust/refs/pol_padroes.txt'.format(_hdtpath()))
#        data = f0.readlines()
#        f0.close()
#        stdref = []
#        for datai in data:
#            stdref += [datai.split()]
    except:
        print('# ERROR: Can\'t read files pyhdust/pyhdust/refs/pol_padroes.txt')
        raise SystemExit(1)

    calc = float(calc)
    angref = []
    thstd = []
    stdname = []
    flagstd = []
    tagstd = []

    if _os.path.exists('{0}/{1}/std.dat'.format(path,night)):
        stds = _np.loadtxt('{0}/{1}/std.dat'.format(path,night), dtype=str)
#        print stds
        if len(stds) > 0 and len(stds[-1]) != 9:
            stds = stds.reshape(-1,9)
#        print stds
        k = 0
        sigs = []
        for stdinf in stds:
            if stdinf[7] == 'E' and stdchk(stdinf[2])[0] and stdinf[3] == f and \
                                                    abs(float(stdinf[4])-calc) <= delta:
                print(('{0:<10s} WARNING! Standard `{1}` ({2}) wasn\'t used because it ' +\
                                'had `E` flag. Skipping this standard...').format(night+':',stdinf[2],f))
                continue
            elif stdinf[7] != 'E' and stdchk(stdinf[2])[0] and stdinf[3] == f and \
                                                    abs(float(stdinf[4])-calc) <= delta:
                # Bednarski: Show error message now
                try:
                    nstar = int(stdinf[6])
                    Q, U, sig, P, th, sigT, tmp, tmp2 = readout('{0}/{1}'.\
                                            format(path+'/'+night,stdinf[5]), nstar=nstar)
                except:
                    print(('{0:<10s} WARNING! Standard `{1}` ({2}) wasn\'t used because' +\
                                            ' can\'t open/read {3}. Skipping this standard...').\
                                            format(night+':', stdinf(2), f, stdinf[5]))
                    continue
                stdname += [ stdinf[2] ]
                thstd += [ float(th) ]
                flagstd += [stdinf[7]]
                tagstd += [stdinf[8]]
                if thstd == 0.:
                    thstd = 0.01
                sigs += [ float(sig)*100 ]
                sigth = 28.65*sig/P 
                i = stdchk(stdinf[2])[1]
                j = filters.index(f)+1 #+1 devido aos nomes na 1a coluna
                angref += [ float(stdref[i,j]) ]

#                print stdname, thstd, angref, flagstd, tagstd
                
        if stdname == []:
            print(('{0:<10s} WARNING! None standard found in `std.dat` for filter {1}.').format(night+':', f))

    else:
        print('{0:<10s} ERROR! `std.dat` file not found (filter {1}).'.format(night+':', f))
#    print stdname, thstd, angref

#    print stdname, thstd, angref, flagstd, tagstd
    return stdname, thstd, angref, flagstd, tagstd



#################################################
#################################################
#################################################
# Bednarski: I CHANGED MANY THINGS TO BE COMPATIBLE WITH THE REMAINING CHANGES.
#            I TESTED ONLY A LITTLE BIT. COMPLETE TESTS MUST BE DONE.
def genTarget(target, path=None, PAref=None, skipdth=False, delta=2.5, epssig=2.0):
    """ Gen. target

    Generate a table with all observations found for 'target', unless
    those which have `E` flag or haven't standard star (except when
    some data is unpolarized or if skipdth=True):

    epssig: sigP/P max for unpolarized target (sigP/P up to epssig
            doesn't need standard star)
    skipdth: print all target values anyway, even when there are no
             standard star?

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
    
    if path == None or path == '.':
        path = _os.getcwd()

    # Carregar angulos de referencia para todos os filtros de um conjunto de padroes
    if PAref is not None:
        for line in PAref:
            if len(line) < 21:
                print('# ERROR: Wrong PAref matrix format.')
                raise SystemExit(1)
    else:
        PAref = _np.loadtxt('{0}/pyhdust/refs/pol_padroes.txt'.format(_hdtpath()), dtype=str)

    # Read lists and verify if target is a valid target

    try:
        obj = _np.loadtxt('{0}/pyhdust/refs/pol_alvos.txt'.format(_hdtpath()), dtype=str)
        if _os.path.exists('{0}/pyhdust/refs/pol_hip.txt'.format(_hdtpath)):
            obj = _np.concatenate((obj,v_np.loadtxt('{0}/pyhdust/refs/pol_hip.txt'.format(_hdtpath()), dtype=str)))
        std = _np.loadtxt('{0}/pyhdust/refs/pol_padroes.txt'.format(_hdtpath()), dtype=str)
    except:
        print('# ERROR: Can\'t read files pyhdust/pyhdust/refs/pol_alvos.txt and/or pyhdust/pyhdust/refs/pol_padroes.txt.')
        raise SystemExit(1)
        
    if target not in _np.hstack((std[:,0],obj)):
        print('# WARNING: Target {0} is not a default target or standard!'.\
        format(target))
        tmp = raw_input('Type something to continue...')

    print('\n'+'='*30+'\n')

    nights = [fld for fld in _os.listdir(path) if _os.path.isdir(_os.path.join(path, fld))]
    lines = ''

    for night in nights:

        # Check obj.dat for the night
        if _os.path.exists('{0}/{1}/obj.dat'.format(path,night)):
            try:
                objs = _np.loadtxt('{0}/{1}/obj.dat'.format(path,night), dtype=str)
            except:
                print('{0:<10s} ERROR! Can\'t read obj.dat file. Ignoring this night...'.format(night+':'))
                continue

            # Verify if std has more than one line. Case not, do the reshape
            if _np.size(objs) == 9:
                objs = objs.reshape(-1,9)
            elif _np.size(objs) % 9 != 0:
                print('{0:<10s} ERROR! Wrong column type in obj.dat file. Ignoring this night...'.format(night+':'))
                continue
#                raise SystemExit(1)

            # Loop on found nights
            for objinf in objs:
                dth = []
                stdnames = ''
                if objinf[2] == target:
                    MJD, ccd, obj, f, calc, out, nstar, flag, tags = objinf
                    if flag == 'E':
                        print(('{0:<10s} WARNING! Star found ({1}), but with `E` flag ' +\
                                        'and tags `{2}`. Ignoring this data...').format(night+':',f,tags))
                        continue
                    try:
                        # Fator is a var to indicate when polarization angle must be taken as 180-theta or +theta
                        fator = thtFactor(float(MJD))
                        Q, U, sig, P, th, sigT, tmp, tmp2 = readout('{0}/{1}'.\
                                                    format(path+'/'+night,out), nstar=int(nstar))
                    except:
                        print('{0:<10s} ERROR! Can\'t open/read out file {1}. Ignoring this data...'.format(night+':',out))
                        continue

                    P = float(P)*100
                    th = float(th)
                    sig = float(sig)*100
                    sigth = 28.65*sig/P
#                    print objinf, objinf[2], target, tags
                    if 'no-std' not in tags:
#                        print('entrou1')
                        stdname, thstd, angref, flagstd, tagstd = corObjStd(night, f, calc, path=path, delta=delta)
#                        print stdname, thstd, angref, flagstd, tagstd
                    else:
#                       FURTHER: APPLY ALTERNATIVE METHOD TO COMPUTE DTHETA HERE AND DELETE LIKE BELLOW.
#                                CAUTION because the line "if stdname == [] and 'no-std' not in tags:" must be changed also
#                        print('entrou2')
                        stdname, thstd, angref, flagstd, tagstd = corObjStd(night, f, calc, path=path, delta=delta)
#                        print stdname, thstd, angref, flagstd, tagstd

                    # Refresh tags and flags
                    if stdname == [] and 'no-std' not in tags:
                        if tags == '---':
                            tags = 'no-std'
                        else:
                            tags += ',no-std'
                        if flag == 'OK':
                            flag = 'W'
                    elif stdname != [] and 'no-std' in tags:
                        if tags == 'no-std':
                            tags = '---'
                        elif tags[0:6] == 'no-std,':
                            tags = tags.replace('no-std,','')
                        else:
                            tags = tags.replace(',no-std','')

                    # Bednarski: working for more than one standard star
                    # 1) Case there is no standard star, set the variables
                    if 'no-std' in tags: # or stdname==[]:
                        mdth = 0.
                        devth = 0.
                        stdnames = '---'
                        thstd=[]
                        angref=[]
                        flagstd=[]
                        tagstd=[]
                    # 2) Case there is some standard star, correct the angle and set variables
                    else:
                        for i in range(len(thstd)):
                            if flag == 'OK' and flagstd[i] == 'W':
                                flag = 'W'
                            # Refresh tag list
                            for tagi in tagstd[i].split(','):
                                if tagi not in tags+',---':
                                    if tags == '---':
                                        tags = ''
                                    else:
                                        tags += ','
                                    tags += tagi
                            dth += [ fator*thstd[i]-angref[i] ]
                            while dth[i] >= 180:
                                dth[i]-= 180
                            while dth[i] < 0:
                                dth[i]+= 180
                            if i != 0:
                                stdnames += ','
                            stdnames += stdname[i]
                        mdth=sum(dth)/len(dth)  # evalute the mean dth
                        devth=_np.std(dth)

                    th = fator*th-mdth
                    # Bednarski: I changed 'if' for 'while' because it can be necessary more than once
                    while th >= 180:
                        th-= 180
                    while th < 0:
                        th+= 180
                    Q = P*_np.cos(2*th*_np.pi/180)
                    U = P*_np.sin(2*th*_np.pi/180)
                    if thstd != [] or skipdth or P/sig <= epssig:
                        if out.find('_WP') == -1:
                            outn = out[-11:]
                        else:
                            outn = out[out.find('_WP')-7:]
                        lines += ('{:12s} {:>7s} {:>7s} {:>4s} {:>5s} {:>12s} {:>6.1f} {:>6.1f}'+
                        ' {:>8.4f} {:>8.4f} {:>8.4f} {:>7.2f} {:>7.4f} '+
                        '{:>6.2f} {:>13s} {:>4s} {:>5s} {:>s}\n').format(MJD, night, ccd, f, calc, stdnames, \
                                    mdth, devth, P, Q, U, th, sig, sigth, outn, nstar, flag, tags)
                        print('{0:<10s} One line added to {1}.log...'.format(night+':', obj))
                    else:
                        print(('{0:<10s} WARNING! No valid delta_theta value estimated in filter {1}.' +\
                                        ' Ignoring this data...').format(night+':', f))

#            print('')
        else:
            if verbose != '':
                verbose += ', '
            verbose += night

    # Print "no obj.dat found" messages
    if verbose != '':
        print('{0:<10s} WARNING! No `obj.dat` file found for the following nights: {1}. Ignoring these nights...'.format('-------',verbose))

    print('\n'+'='*30+'\n')
    #arquivo de saida
    if lines != '':
        print('DONE! Output written in {0}/{1}.log.'.format(path,target))
        f0 = open('{0}/{1}.log'.format(path,target),'w')
        lines = ('{:12s} {:>7s} {:>7s} {:>4s} {:>5s} {:>12s} {:>6s} {:>6s}' +\
                        ' {:>8s} {:>8s} {:>8s} {:>7s} {:>7s} {:>6s} {:>13s}' +\
                        ' {:>4s} {:>5s} {:>s}\n').format('#MJD', 'night',\
                        'ccd', 'filt', 'calc', 'stdstars', 'dth', 'devdth', 'P', 'Q', 'U',\
                        'th', 'sigP', 'sigth', 'outfile', 'star', 'flag', 'tags')+lines
        f0.writelines(lines)
        f0.close()
    else:
        print('NOT DONE! No valid observation was found for target `{0}`.'.format(target))
      
    return



#################################################
#################################################
#################################################
def pur2red():
    """
    Pure to Reduced
    """
    return



#################################################
#################################################
#################################################
def red2pur():
    """
    Reduced to Pure
    """
    return



#################################################
#################################################
#################################################
def genJD(path=None):
    """Generate de JD file for the fits inside the folder
    """
    if path == None:
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
    ltgts = _np.loadtxt('{0}/pyhdust/refs/pol_alvos.txt'.format(_hdtpath()), dtype=str)
    lstds = _np.loadtxt('{0}/pyhdust/refs/pol_padroes.txt'.format(_hdtpath()), dtype=str,\
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
    if path == None:
        path = _os.getcwd()
    lmags = _np.loadtxt('{0}/pyhdust/refs/pol_mags.txt'.format(_hdtpath()), dtype=str)

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
    Set CCD name in global variable 'ccd'
    """
    global ccd
    ccd = ''

    if fitsfile != '':
        try:
            fits = _pyfits.open(fitsfile)
            instrume = '{0}'.format(fits[0].header['SERNO'])
            if instrume.find('4335') != -1:
                ccd = 'ixon'
            elif instrume.lower().find('10127') != -1:
                ccd = 'ikon'
            elif instrume.lower().find('9867') != -1:
                ccd = 'ikon'
            else:
                ccd = ''
        except:
            pass

    while ccd not in ('ikon', 'ixon', '301'):
        ccd = raw_input('Type the CCD name (301/ikon/ixon): ')



#################################################
#################################################
#################################################
def graf_dtheta(fname, save=False):
    """
    Plot graph for delta theta study.

    fname  --  is the data file. Its lines must be of type
           "date standard_name lambda dtheta", separeted by a
            simple space.
    Blank lines will be skiped by the algorithm.
    """

    _plt.close('all')
    #_rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

    # Read file
    try:
        file0 = _np.loadtxt(fname, dtype=str, delimiter=' ')
    except:
        print('Error when reading the file or file doesn\'t exist!')
        return
    flist = file0.tolist()
    flist.sort(key=lambda x: [x[1],x[0],x[2]])

    # Delete void lines
    i=0
    while i < len(flist):
        if flist[i][0] == '':
            del(flist[i])
        else:
            i += 1
    flist += [['','',0,0]]    # To work for last line

    # Count number of plots
    num_plots = 0
    for i in range(len(flist)-1):
        if flist[i][0]!=flist[i-1][0] or flist[i][0]!=flist[i-1][0]:
            num_plots += 1

    # Generate figure and axes
    fig = _plt.figure(1)
    ax = _plt.subplot(1, 1, 1)
    #ax.set_title('Estudo delta theta', fontsize=16, verticalalignment='bottom')
    ax.set_xlabel(r'$\lambda \ (\AA)$', size=15)
    ax.set_ylabel(r'$\Delta\theta\ -\ <\Delta\theta>$ (degree)', size=15)
    ax.set_xlim([3000, 9000])
    colormap = _plt.cm.spectral
    _plt.gca().set_color_cycle([colormap(i) for i in _np.linspace(0, 0.9, num_plots)])

    # Plot the graphs
    for i in range(len(flist)):
        if i==0 or flist[i][0]!=flist[i-1][0] or flist[i][0]!=flist[i-1][0]:
            if i!=0:
                ax.plot(vals[0], vals[1], label = '{0} ({1})'.format(flist[i-1][1], \
                                    flist[i-1][0]), linestyle='-', marker='o')
            vals = [[],[]]
        vals[0]+=[float(flist[i][2])]
        vals[1]+=[float(flist[i][3])]

    # Final changes
    ax.autoscale(False)
    ax.plot(ax.get_xlim(), [0,0], 'k--')
    ax.legend(loc='lower left', borderaxespad=0., numpoints=1, prop={'size':11})

    #ax.xaxis.label.set_fontsize(14)
    #ax.yaxis.label.set_fontsize(14)
    for item in (ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(15)

    if save:
        _plt.savefig('dtheta.pdf', bbox_inches='tight')
    else:
        _plt.show()



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

def loadpol(txt):
    """ Load polarization txt file. """
    dtb = _np.loadtxt(txt, dtype=str)
    dtb = _np.core.records.fromarrays(dtb.transpose(), names='MJD,night,filt,\
    calc,stdstars,dth,devdth,P,Q,U,th,sigP,sigth', formats='f8,{0},{0},f8,{0},\
    f8,f8,f8,f8,f8,f8,f8,f8'.format(dtb.dtype))
    return dtb

#################################################
#################################################
#################################################
## AINDA FALTA TESTAR!!!
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

            # If filename is 'logfile'
            if f == 'logfile':
                file_new = path_red + '/' + _os.path.relpath(file_old, path)
            # If file is inside dark/flat/bias path
            elif file_old.find('/dark/') != -1 or file_old.find('/flat/') != -1 or \
                                            file_old.find('/bias/') != -1:
                if _re.search(r'.fits$', f):
                    file_new = path_red + '/' + night + '/' + calib
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
                elif len(elem) > 1:
                    filt = elem[1]
                else:
                    filt = 'wf'

                # Case sum*.dat files
                if _re.search(r'^sum', f):
                    file_new = path_red + '/' + obj + '/sum_' + obj + '_' + filt + '_' + f[-6:]
                # Case w*.dat, w*.log, w*.out files
                elif _re.search(r'.[dat|log|out]$', f):
                    file_new = path_red + '/' + obj + '/w' + obj + '_' + filt + '_' + f[1:]
                # Case coord.*.ord files
                elif _re.search(r'^coord', f):
                    file_new = path_red + '/' + obj + '/coord_' + obj + '_' + filt + f[-6:]
                # Case JD* files
                elif _re.search(r'^jd[0-9]*', f):
                    # PAREI AQUI. CONCATENAR TODOS AQUIVOS JD EM UM SÓ
                    pass
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
### MAIN ###
if __name__ == "__main__":
    pass

