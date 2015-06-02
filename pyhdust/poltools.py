#-*- coding:utf-8 -*-

"""
PyHdust *poltools* module: polarimetry tools

History:
-Added options `force` and `chknames` to genStdLog and genObjLog

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
           2: ['wrong object suspected','wrong-obj?'],
           3: ['polarimeter problem suspected','iagpol-prob?'],
           4: ['another relevant problem','other-prob'],
          }

# Dictionary for the tags assigned automatically
# If you want to add another value, add inside verout routin also.
dictests = {0: ['no reduced','no-red', 'E'],
            1: ['std incompatible with the published','wrong-std?', 'W'],
            2: ['sig < theorical_sig','s<theor_s', 'W'],
            3: ['no standard in the night','no-std', 'W'],
           }



#################################################
#################################################
#################################################
def stdchk(stdname):
    """ Check it the standard star name contains an known name, and return
    its position in `padroes.txt` """
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
def readout(out):
    """
    Read the *.out file from IRAF reduction and return a float array

    Q   U        SIGMA   P       THETA  SIGMAtheor.  APERTURE  STAR
    """
    f0 = open(out)
    data = f0.readlines()
    f0.close()
    data = data[1].split()
    return [float(x) for x in data]



#################################################
#################################################
#################################################
def readoutMJD(out, obj=None, verbose=True):
    """
    Read the 'out' file from IRAF reduction in a float array (fout),
    appending the MJD date and the angle of the beams from
    calcite. Stack error/warning messages in a string (log).

    Return: fout, log    
    - If unable to read 'out' file, assign fout as []
    - If unable to read JD file, assign MJD value as -1
    - If unable to read coords file, assign beams angle as -1

    If verbose==True, print log also.
    'obj' (object name) is optional and used only in error messages
    """

    loglines = ''
    path = _phc.trimpathname(out)[0]
    outn = _phc.trimpathname(out)[1]
    try:
        data = readout(out)
    except:
        if obj != None:
            loglines += '# ERROR: {0}: Can\'t open/read file {1}. Verify and write the values manually.\n'.format(obj,out)
        else:
            loglines += 'ERROR: Can\'t open/read file {0}. Verify and write the values manually.\n'.format(out)
        if verbose:
            print(loglines)
        return [], loglines
        
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
#        print obj, f, _glob('{0}/JD_*_{1}'.format(path,f))
        if obj != None:
            loglines += '# CAUTION: {0}_{1}: No JD file found. Edit the value manually.\n'.format(obj,f)
        else:
            loglines += '# CAUTION: No JD file found as {0}/JD_*_{1}. Edit the value manually.\n'.format(path,f)           
        date = -1

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
                if obj != None:
                    loglines += '# CAUTION: {0}_{1}: No COORDS file found. Edit the value manually.\n'.format(obj,f)
                else:
                    loglines += '# CAUTION: No COORDS file found as {0}/coord_*_{1}_*.ord. Edit the value manually.\n'.format(path,f)
                ang = -1

    if len(coords) != 0:
        try:
            coords = _np.loadtxt(coords[0])
            ang = _np.arctan( (coords[1,1]-coords[0,1])/(coords[1,0]-coords[0,0]) )*180/_np.pi
            while ang < 0:
                ang += 180.
        except:
            if obj != None:
                loglines += '# CAUTION: {0}_{1}: Can\'t open coords file. Edit the value manually.\n'.format(obj,f)
            else:
                loglines += '# CAUTION: Can\'t open coords file at {0}/coord_*_{1}_*.ord. Edit the value manually.\n'.format(path,f)
            ang = -1

    if date != -1:
        if datei == datef:
            print('# Strange JD file for '+out)
        date = (datef+datei)/2

    if verbose and loglines != '':
        print(loglines)

    return [float(x) for x in data]+[date]+[ang], loglines



#################################################
#################################################
#################################################
def minErrBlk16(night,f,i):
    """Retorna o menor erro num bloco de 16 posicoes polarimetricas
    (i.e., arquivo de menor erro entre [08(i),08(i+9)] e [16(i)]).
    """

    err = _np.ones(3)
    out = _np.zeros(3, dtype='|S256')

    # First: Get 1st 8WP group
    # Bednarski: substituí pra usar regexp e funcionar com o sufixo depois do filtro nos nomes.
    ls = [night+'/'+fl for fl in _os.listdir('{0}'.format(night)) if _re.search(r'_{0}'. \
                format(f,night) + r'_.*_?08' + r'{0:03d}\..\.out'.format(i), fl)]
    if len(ls) > 0:
        err[0] = float(readout(ls[0])[2])
        out[0] = ls[0]
        for outi in ls:
            if float(readout(ls[0])[2]) < err[0]:
                err[0] = float(readout(outi)[2])
                out[0] = outi
    
    # Second: Get 2nd 8WP group
    ls = [night+'/'+fl for fl in _os.listdir('{0}'.format(night)) if _re.search(r'_{0}'. \
                format(f) + r'_.*_?08' + r'{0:03d}\..\.out'.format(i+8), fl)]
    if len(ls) > 0:
        err[1] = float(readout(ls[0])[2])
        out[1] = ls[0]
        for outi in ls:
            if float(readout(ls[0])[2]) < err[0]:
                err[1] = float(readout(outi)[2])
                out[1] = outi

    # Third: Get 16WP group
    ls = [night+'/'+fl for fl in _os.listdir('{0}'.format(night)) if _re.search(r'_{0}'. \
                format(f) + r'_.*_?16' + r'{0:03d}\..\.out'.format(i), fl)]
    if len(ls) > 0:
        err[2] = float(readout(ls[0])[2])
        out[2] = ls[0]
        for outi in ls:
            if float(readout(ls[0])[2]) < err[0]:
                err[2] = float(readout(outi)[2])
                out[2] = outi

    j = _np.where(err == _np.min(err))[0]

    return err[j], out[j]



#################################################
#################################################
#################################################
def chooseout(objdir, obj, f, sigtol=lambda sig: 1.1*sig +0.00005):
    """
    Olha na noite, qual(is) *.OUT(s) de um filtro que tem o menor erro.

    Retorna um mais valores para toda a sequencia (i.e., pasta).
    Criterios definidos no anexo de polarimetria.

    minerror == True: recebe o out de menor erro em qualquer caso.
                False: so ocorre se numero de posicoes <= 16.

    sigtol: tolerancia para o sigma (válido apenas se houver < 16 posicoes
    de lamina). Se ha N<=16 pos de lâmina e o erro do agrupamento de N
    posicoes for menor que a funcao sigtol sobre o erro do agrupamento de menor erro,
    entao usa o agrupamento com as N posicoes; do contrario usa o de menor erro.
    Exemplo com 16 posicoes: erro do grupo de 16 posicoes == 0.000140;
        menor erro == 0.000100. Se sigtol(sig) = 1.1*sig+0.00005, como
        sigtol(0.000100) == 0.000160 > 0.000140, usa o agrupamento de 16 posicoes.

    O numero de posicoes eh baseado no numero de arquivos *.fits daquele
    filtro.
    """

    # Verify if queryout was the function that has called this chooseout
    if _getouterframes(_currentframe(), 2)[1][3] == 'queryout':
        queryverif = True
    else:
        queryverif = False

    loglines = ''
    npos = len(_glob('{0}/*_{1}_*.fits'.format(objdir,f)))
    louts = _glob('{0}/*_{1}_*.out'.format(objdir,f))

    # Verify if there are more than one star inside outfiles (working only to the star #1)
    for fout in louts:
        tmpfile = open(fout)
        count = sum(1 for line in tmpfile if line.rstrip('\n'))
        tmpfile.close()
        if count > 2:
            print('\n# CAUTION: {0}_{1}: There are more than one star inside .out files.'\
               .format(obj, f)+' Only star #1 was used! Check and add the others manually.')
            loglines += '# CAUTION: {0}_{1}: There are more than one star inside .out files.'\
               .format(obj, f)+' Only star #1 was used! Check and add the others manually.\n'
            break

    # Check reduction
    if len(louts) == 0:
#        print obj, objdir, f
        outs = []
        if npos != 0:
            tests, loglines = verout('', obj, f, verbose=True)
            return [''], tests, loglines
    #Se ateh 16 npos, pega o de menor erro
    elif npos <= 16:
        err0 = 10.
        err1 = 10.
        for outi in louts:
            # Calculate min error for .out file with all WP positions
            if outi.find('_{0:2d}'.format(npos)) != -1:
                if float(readout(outi)[2]) < err0:
                    err0 = float(readout(outi)[2])
                    out0 = outi
            # Calculate min error for the others .out files
            else:
                if float(readout(outi)[2]) < err1:
                    err1 = float(readout(outi)[2])
                    out1 = outi
        if err0 == 10. and err1 == 10.:
            print npos,louts
            print('# ERROR: Something went wrong with *.out sigma values!')
            raise SystemExit(1)
        elif err1 == 10. or err0 <= sigtol(err1): outs = [out0]
        elif err0 == 10. or err0 > sigtol(err1): outs = [out1]
        else:
            print('# ERROR: It shouldn\'t enter here in chooseout!')
            raise SystemExit(1)
        
    # Se a partir de 16, faz o seguinte criterio: dentro de npos%8, ve qual
    # posicao inicial do block de 16 (ou 8+8) tem o menor erro, e joga no valor
    # de i.
    else:
        i = 1
        err1 = minErrBlk16(objdir,f,i)[0]
        #outs = list(minErrBlk16(objdir,f,i)[1])
        #print 'b',outs
        # Case lamina position is among 17,18,...,23 , 25,26,...,31 , 33,34,...,39 , ...
        for j in range(1,npos%8+1):
            if minErrBlk16(objdir,f,j+1)[1] < err1:
                err1 = minErrBlk16(objdir,f,j+1)[0]
                i = j+1
                #outs = list(minErrBlk16(objdir,f,j+1)[1])
                #print 'c',outs
        outs = []
        while i+16-1 <= npos:
            outi = minErrBlk16(objdir,f,i)[1][0]
            if outi.find('_{0}_16'.format(f)) > -1:
                outs += [outi]
            else:
                for j in [i,i+8]:
#                    ls = _glob('{0}/*_{1}_08{2:03d}*.out'.format(objdir,f,j))
                    ls = [objdir+'/'+fl for fl in _os.listdir('{0}'.format(objdir)) \
                            if _re.search(r'_{0}'.format(f) + r'_.*_?08' + \
                                        r'{0:03d}\..\.out'.format(j), fl)]
                    if len(ls) != 2:
                        print(('# Warning! Check the *.out 2 '+\
                        'versions for filter {0} of {1}').format(f,obj))
                        if len(ls) == 1:
                            out = ls[0]
                        #else:
                        #    raise SystemExit(1)
                    else:
                        if float(readout(ls[0])[2]) < float(readout(ls[0])[2]):
                            out = ls[0]
                        else:
                            out = ls[1]
                        outs += [out]
            i += 16
        #Se apos o bloco de 16 ha 8 pontos independentes, ve dentre eles qual
        #tem o menor erro.
        if i <= npos-8+1:
            outi = [objdir+'/'+fl for fl in _os.listdir('{0}'.format(objdir)) \
                            if _re.search(r'_{0}'.format(f) + r'_.*_?08' +  \
                                r'{0:03d}\..\.out'.format(i), fl)]
            if len(outi) == 0:
                print('# ERROR! Strange *.out selection! Case 1')
                print('# Probably 08pos files missing.')
                print('{0}/*_{1}_[suf_]08{2:03d}*.out'.format(objdir,f,i))
                print i,npos,npos-8+1,objdir,f
                raise SystemExit(1)
            else:
                err1 = float(readout(outi[0])[2])
                if len(outi) == 1:
                    outi = outi[0]
                elif len(outi) == 2:
                    if float(readout(outi[1])[2]) < err1:
                        err1 = float(readout(outi[1])[2])
                        outi = outi[1]
                    else:
                        err1 = float(readout(outi[0])[2])
                        outi = outi[0]                
                else:
                    print('# ERROR! Strange *.out selection! Case 2')
                    print('# Probably there is position-excluded POLRAP *.out...')
                    print('{0}/*_{1}_[suf_]08{2:03d}*.out'.format(objdir,f,i))
                    print i,npos,npos-8+1,outi,objdir,f
                    raise SystemExit(1)
            for j in range(i+1,npos+1-8+1):
                tmp = [objdir+'/'+fl for fl in _os.listdir('{0}'.format(objdir)) \
                                if _re.search(r'_{0}'.format(f) + r'_.*_?08' + \
                                r'{0:03d}\..\.out'.format(j), fl)]
                #print j, npos+1-8+1, len(tmp), tmp
                if len(tmp) == 1:
                    tmp = tmp[0]
                elif len(tmp) == 2:
                    if float(readout(tmp[1])[2]) < float(readout(tmp[0])[2]):
                        tmp = tmp[1]
                    else:
                        tmp = tmp[0]
                else:
                    print('# ERROR! Strange *.out selection! Case 3')
                    print('# Probably there is position-excluded POLRAP *.out...')
                    print('{0}/*_{1}_[suf_]08{2:03d}*.out'.format(objdir,f,j))
                    print i,npos,npos-8+1,tmp,outi
                    raise SystemExit(1)                
                if float(readout(tmp)[2]) < err1:
                    err1 = float(readout(tmp)[2])
                    outi = tmp
            outs += [outi]

    # If there is no fits files for such object/filter
    if outs == []:
#        print('********* HEY: Here didn\'t reduced? ')
        return [''], [[]], loglines
    else:
        tests = []
        for out in outs:
            testout, logout = verout(out, obj, f, verbose=False)
            tests += [testout]
            if not queryverif:
                loglines += logout

    return outs, tests, loglines



#################################################
#################################################
#################################################
def verout(out, obj, f, verbose=True):
    """
    Function to do tests on outfile 'out' concerning to
    object 'obj' in filter 'f'.
    Tests: (1) test if file 'out' exists.
           (2) test if P_pub for standards is compatible with
               P_obs value.
           (3) test if sig > sig_theorical.

    - If verbose==True, show warnings in screen
    - In objdir==None, outfile is supposed in current dir

    Return a boolean list with three components concerning
    to the tests (1)-(3) above + log string. If some test has failed,
    the concerning value is assigned as \"True\"; otherwise,
    \"False\".
    """
    tests = [False]*len(dictests)
    loglines = ''
#    if objdir != None:
#        out = objdir+'/'+out

    try:
        [Q,U,sig,P,th,sigT,ap,star,MJD,calc],loglines = readoutMJD(out, obj=obj, verbose=False)
        sig_ratio = float(sig)/float(sigT)
        ztest = verStdPol(obj, f, float(P)*100, float(sig*100))
    except:
        tests[0] = True
        loglines += '# ERROR: {0}_{1}: Reduction files were not found.\n'.format(obj,f)
        if verbose:
            print(loglines)
        return tests, loglines
    
    # Some tests.
    if ztest > 2.5:     # Case the object is not a standard, ztest==-1 and tests[0] remains False.
        tests[1] = True
    if sig_ratio <= 1.:
        tests[2] = True
    
    # Print tests
    if tests[1]:
        loglines += ('# WARNING: {0}_{1}: The standard has polarization only compatible '+\
                                  'within {2:.1f} sigma with the published value.\n').format(obj, f, ztest)
    if tests[2]:
        loglines += ('# WARNING: {0}_{1}: Polarization has sig < theorical_sig ' +\
                        '({2:.4f} < {3:.4f}).\n').format(obj, f, sig*100, sigT*100)
    if verbose and loglines != '':
        print('\n'+loglines)
    
    return tests, loglines



#################################################
#################################################
#################################################
def queryout(obj, objdir, f, sigtol=lambda sig: 1.1*sig +0.005):
    """
    Call chooseout for 'obj' at filter 'f' (in 'objdir'),
    print the graphs and query to the user if he wants
    to use the selected outs. If not, he musts answer
    what to use.

    """

    loglines = ''
    _plt.close('all')
    outs, tests, loglines = chooseout(objdir, obj, f, sigtol=sigtol)

    if outs == ['']:
        return outs, tests, None, None, loglines

    sortout = grafall(objdir, f, bestouts=outs)

    while True:
        print('\n'+'_'*80)
        print(' {0:<10s} {1:<5s} {2:<7s} {3:<8s} {4:<10s} {5:<7s} {6:<s}'.format('Obj', 'Filt', \
                'P(%)', 'sig(%)', 'sig/ThSig', 'ztest', 'out/num'))

        for out in outs:
            if out != '':
                try:
                    [Q,U,sig,P,th,sigT,ap,star,MJD,calc],loglixo = readoutMJD(out, obj=obj, verbose=False)
                    sig_ratio = float(sig)/float(sigT)
                    z = verStdPol(obj, f, float(P)*100, float(sig*100))
                    numout = '(#{0})'.format(sortout.index(out))
                except:
                    print('ERROR! It shouldn\'t enter here in queryout!')
                    raise SystemExit(1)

                # Reassigns ztest value for printing
                if z == -1:
                    zstr = '-----'
                else:
                    zstr = '{0:.1f}'.format(z)
                print(' {0:<10s} {1:<5s} {2:<7.4f} {3:<8.4f} {4:<10.3f} {5:<7s} {6:<s} {7:<s}'.\
                         format(obj.upper(), f.upper(), float(P)*100, float(sig)*100, \
                                sig_ratio, zstr, _phc.trimpathname(out)[1], numout))
        print('_'*80+'\n')

        # Get tests for outfiles
        tests = []
        logs = ''
        # All out files shouldn't be '', unless the object was not reduced
        # and the task should has enter in a previous 'if'
        for out in outs:
            if out == '':  # Case is not to use such out
                tests += [[True] + [False]*(len(dictests)-1)]
            else:          # Otherwise, test the out file
                testout, logout = verout(out, obj, f, verbose=True)
                tests += [testout]
                logs += logout

        while True:
            verif = raw_input('Use this(ese) out(s)? (y/n): ')
            if verif not in ('y','Y','n','N'):
                print('Invalid choise!')
            else:
                break

        if verif in ('y', 'Y'):
            break
        elif verif in ('n', 'N'):
            opt=''  # for the first iteration
            for i in range(len(outs)):
                while True:
                    if len(outs) == 1:
                        # At least one out file shall be attributed
                        opt = raw_input('Type the out number: '.format(i, len(outs)))
                    else:
                        if opt != '0':
                            # At least one out file shall be attributed
                            if i == 0: opt = raw_input('Type the 1st out number: ')
                            elif i == 1: opt = raw_input('Type the 2nd out number or \'0\' to finish: ')
                            elif i == 2: opt = raw_input('Type the 3rd out number or \'0\' to finish: ')
                            else: opt = raw_input('Type the {0}th out number or \'0\' to finish: '.format(i+1))
                    # if opt == '0', don't use the outfiles for current and next i values of loop
                    if opt == '0' and i != 0:
                        outs[i] = ''
                        break
                    # If opt is a valid value, assign the input number with the concerning out file
                    elif opt in [str(j) for j in range(1,len(sortout))]:
                        outs[i] = sortout[int(opt)]
                        break
                    else:
                        print('Wrong value!')
                        opt=''

    # Request what tags assign (flexible by means of dictags global dictionary)
    print('\n# TAGS LIST\n  0: none')
    for i in dictags.keys():
        print('  {0}: {1}'.format(i+1, dictags[i][0]))
    print('')

    while True:
        verif = True
        tags = [False for i in dictags.keys()]
        strin = raw_input('Select all tags that apply separated by commas (\'0\' for none): ')
        if strin == '0':
            flag='OK'
            break
        opts = strin.split(',')
        for opt in opts:
            if opt in [str(j+1) for j in dictags.keys()]:
                opt = int(opt)-1
                tags[opt] = True
            else:
                print('Invalid choise!')
                verif = False
                break

        # If some tag was selected, request a flag below
        if verif:
            flag=''
            while flag not in ('E','e','W','w','OK','ok','Ok'):
                flag = raw_input('Select one flag for usage: [W]arning, [E]rror or [OK]: ')
            flag = flag.upper()
            break

    loglines += logs
    _plt.close('all')

    return outs, tests, tags, flag, loglines



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
def grafall(objdir, filt, bestouts=[], onlyhigh=False):
    """
    Multiple plot modulations for the object inside \'objdir\'
    in filter \'filt\'.
    Return a list 'sortout' with the log files sorted by plot
    number. Allways sortout[0]=='' and sortout[i] is the plot #i.

    Optional:
        bestouts - list of outs which have least error to
                   highlight in figures
        onlyhigh - plot only the group with the highest number
                   of WP? (True/False)
    """

    # Receive a list of logfiles and return lists for each group of WP/version:
    # sublogs == list of sorted logfiles.
    # groups == list with informations about lamina groups.
    # ver == list with the reduction versions of files in sublogs.
    #
    # Ex, sublogs == [[*16001.1.log], [*16001.2.log], [*0800[1-9].1.log], [*0800[1-9].2.log]]
    #      groups == [   [16, 1, .1],    [16, 1, .2],         [8, 9, .1],         [8, 9, .2]]
    #         ver == [          [.1],         [best],           [.1 ...],          [.2, ...]]
    def combineout(logs):

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
                    if onlyhigh and logsplit[i][1][0:2] != logsplit[i+1][1][0:2]:
                        break
                    j+=1
                    sublogs[:] += [[]]
                    ver[:] += [[]]
                    
        return groups, sublogs, ver
        

    # Generate background color for the graphs, depending the reduction version
    def gencolour(ver):

        if ver == 'best': bkg = '#ffbbbb'
        elif ver == '.1': bkg = '#f3ffd8'
        elif ver == '.2': bkg = '#e3ecf9'
        elif ver == '.3': bkg = '#ffeee5'
        elif ver[:4] == '.1_WP': bkg = '#f3ffd3'
        elif ver[:4] == '.2_WP': bkg = '#e3ecf4'
        elif ver[:4] == '.3_WP': bkg = '#faeee7'
        else: bkg = '#f5f5f5'

        return bkg

    
    # Generate graphs for the cases nlog <= 8
    def gengraphl4(logs, ver, count):

        nlin = 2
        ncol = len(logs)/2
        if ncol > 4:
            print("Error: one figure was not displayed")
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
        for j in range(ncol):
            ax += [ fig.add_subplot(grids[0][0,j]),\
                    fig.add_subplot(grids[0][1,j]) ]
        for j in range(ncol):
            ax += [ fig.add_subplot(grids[1][2,j]),\
                    fig.add_subplot(grids[1][3,j]) ]
        for j in range(len(logs)):
            grafpol(logs[j], fig, ax[2*j], ax[2*j+1])
            ax[2*j].set_axis_bgcolor(gencolour(ver[j]))
            ax[2*j+1].set_axis_bgcolor(gencolour(ver[j]))
            ax[2*j].text(0.85, 0.85, '#{0:<2d}'.format(count+j), \
                horizontalalignment='left', verticalalignment='center', style='italic', \
                transform=ax[2*j].transAxes, fontsize=20, color='red')
            
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
                grafpol(logs[j+12*i], fig, ax[2*j], ax[2*j+1])
                ax[2*j].set_axis_bgcolor(gencolour(ver[j]))
                ax[2*j+1].set_axis_bgcolor(gencolour(ver[j]))
                ax[2*j].text(0.85, 0.85, '#{0:<2d}'.format(count+j+i*12), \
                        horizontalalignment='left', verticalalignment='center', \
                        style='italic', transform=ax[2*j].transAxes, fontsize=20, color='red')

            _plt.show(block=False)


    logs = _glob('{0}/*_{1}_*.log'.format(objdir, filt))
    if logs == []:
        print('ERROR: log files not found to plot. May the file names \
              {0}/*_{1}_*.log are wrong!'.format(objdir, filt))
        return 1
    gps, sublogs, ver = combineout(logs)
    nlog = sum([len(subb) for subb in sublogs])
#    print(gps, sublogs, ver)

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
def grafpol(filename, fig=None, ax1=None, ax2=None, save=False, extens='png'):
    """
    Program to plot the best adjust and its residuals of IRAF reduction.
    'filename' is the path to the .log output file from reduction.

    NEW: Working for *_WP1110....log files!

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

        file = _np.loadtxt(filename, dtype=str, delimiter='\n')

        # npts: number of WP valid to be plotted.
        # totpts: used for the total number of WP (for the *.2_WP11110...1.out)
        npts = int(file[9].split()[-1])
        tnpts = int(file[8].split()[-1])
        delta = float(file[14].split()[-1])
        sigma = 1.

        # Bednarski: corrected (25 -> 20, because the blank lines had been ignorated by
        #            np.loadtext function)
        for i in range(20, len(file)):
            if 'APERTURE' in file[i]:
                sig = float(file[i+2].split()[2])
                if sig < sigma:
                    sigma = sig

                    # Bed: Os Q e U sao os abaixo, conforme copiei da rotina graf.cl
                    if float(file[i+2].split()[4]) < 0:
                        thet = - float(file[i+2].split()[4])
                    else:
                        thet = 180. - float(file[i+2].split()[4])

                    # Recalculating the new QU parameters
                    Q = float(file[i+2].split()[3])*_np.cos(2.*thet*_np.pi/180.)
                    U = float(file[i+2].split()[3])*_np.sin(2.*thet*_np.pi/180.)
                    n = npts/4 
                    if npts%4 != 0:
                        n = n+1
                    P_pts = []
                    for j in range(n):
                        P_pts += file[i+4+j].split()
#                    print P_pts
                    P_pts = _np.array(P_pts, dtype=float)
                    th_pts = 22.5*_np.arange(tnpts)-delta/2.
                    j = filename.find('.')
                    delta2 = int(filename[-2+j:j])-1
                    # Bed: Funcionando para nlam >= 10  para impressão correta
                    str_pts = map(str, _np.arange(1,tnpts+1)+delta2)

                    # Case _WP11110...1.log file
                    if npts != tnpts:
                        refs = file[9].split()[3:-2]
                        rm = [j for j in range(tnpts) if refs[j] == '0']
                        th_pts = _np.delete(th_pts, rm)
                        str_pts = _np.delete(str_pts, rm)

        if sigma == 1.:
            print('# ERROR reading the file %s !' % filename)
            Q = U = 0
            P_pts = th_pts = _np.arange(1)
            str_pts = ['0','0']

        return(Q, U, sigma, P_pts, th_pts,str_pts)


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
#        ax1.set_ylim([min(P_pts*100)*1.1, max(P_pts*100)*1.1])

        ax2.set_xlabel('WP position', size=9)
        ax2.set_ylabel('Residuals', size=9)
        ax2.set_xlim([th_pts[0]-4,th_pts[-1]*1.02+3])
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
        Q, U, sigma, P_pts, th_pts, str_pts = readlog(filename)
        plotlog(ax1,ax2, Q,U,sigma,P_pts,th_pts,str_pts,filename)
        if save:
            _plt.savefig(filename.replace('.log','.'+extens))
        else:
            _plt.show()
    else:
        Q, U, sigma, P_pts, th_pts, str_pts = readlog(filename)
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
    lstds = _np.loadtxt('{0}/refs/pol_padroes.txt'.format(_hdtpath()), dtype=str)

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
    another tags already assigned and current flag.
    If they are given as input, return 'tags'+tags
    concerning to 'tests' list and the resulting flag.
    """

    tagstr = ''

    # Generate a string for such tags
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
def genStdLog(path, subdirs, stds, sigtol=lambda sigm: 1.1*sigm + 0.00005, \
                    autochoose=False):
    """
    Generate Standards Log
    """

    if len(stds) != len(subdirs):
        print('\nERROR: polt.genStdLog() NOT RUNNED! (len(stds) != len(subdirs))')
        writeLog(path, '# ERROR: polt.genStdLog() NOT RUNNED! (len(stds) != len(subdirs))\n')
        return 1

    lstds = _np.loadtxt('{0}/refs/pol_padroes.txt'.format(_hdtpath()), dtype=str)

    i=0
    lines=''
    loglines = ''
    maxsize = 7  # var to align the columns (initial value concerning to 'outfile' header)

    # Loop on list of standards
    for i in range(len(stds)):
        obj = stds[i]
        objdir = subdirs[i]
        if obj == '':
            continue
#        print('\n\n**{0}**: outfile(s) with least error '.format(obj.upper()))
        for f in filters:
            if autochoose:
                outs, tests, logs = chooseout('{0}/{1}'.format(path,objdir), obj, f, sigtol=sigtol)
                tags=None
                flag=None
            else:
                outs, tests, tags, flag, logs = queryout(obj, '{0}/{1}'.format(path,objdir), f, sigtol=sigtol)
            loglines += logs

            logs=''
            # Loop on outfiles
            for j in range(len(outs)):
                if outs[j] != '':
                    [Q,U,sig,P,th,sigT,ap,star,MJD,calc],loglixo = readoutMJD(outs[j],obj=obj,verbose=False)
                    tagstr, flagout = readTests(tests[j], tags=tags, flag=flag)
                    lines += '{0:12.6f} {1:>7s} {2:>10s} {3:>4s} {4:>5.1f} {5:<s} {6:>5s}  {7:<s}\n'.\
                                format(MJD, ccd, obj, f, float(calc), ':::'+_os.path.relpath(outs[j], path)+':::', flagout, tagstr)
                    if len(_os.path.relpath(outs[j], path)) > maxsize:
                        maxsize = len(_os.path.relpath(outs[j], path))

            # Only the last 'logs' (concerning to the outs[-1]), not each one
#            if logs != '':
#                loglines += logs

            # Case all components of outs are the void string '' (object not reduced)
            # tests != [[]] is to filter objects+filter not observed
            if outs.count('') == len(outs) and tests != [[]]:
                tagstr, flagout = readTests([True] + [False]*(len(dictests)-1), tags=tags, flag=flag)
                lines += '{0:12.6f} {1:>7s} {2:>10s} {3:>4s} {4:>5.1f} {5:<s} {6:>5s}  {7:<s}\n'.\
                        format(-1, ccd, obj, f, -1, ':::'+objdir+'/-----:::', flagout, tagstr)
                if len(objdir)+6 > maxsize:
                    maxsize = len(objdir)+6

    if lines == '':
        loglines += '# WARNING! No valid standards were found by polt.genStdLog().\n'
    else:
        lines = '{:12s} {:>7s} {:>10s} {:4s} {:>5s} {:<s} {:>5s}  {:<s}\n'.format('#MJD','ccd',\
        'target','filt','calc',':::outfile:::','flag','tags')+lines

    # Realign columns
    lines_tmp = lines.split('\n')
    lines = ''
    for line in lines_tmp:
        if line == '':
            continue
#        print line
        line = line.split(':::')
#        print line
        lines += line[0][:]
        lines += line[1][:].rjust(maxsize+2)
        lines += line[2][:]
        lines += '\n'
    
    f0 = open('{0}/std.dat'.format(path), 'w')
    f0.writelines(lines)
    f0.close()
    writeLog(path, loglines)

    return 0


    
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

    # If calcite is not known, standard star can't be applied
    if calc == -1:
        return False

    loglines = ''
    if path == None:
        path = _os.getcwd()

    # Read `obj.dat` and `std.dat`. If there are errors, assigns [''] to get inside
    # ifs below and print error messages
    try:
        std = _np.loadtxt('{0}/std.dat'.format(path), dtype=str)
    except:
        std = _np.array([''], dtype=str)

    # Verify if std.dat has more than one line. Caso no, do reshape (transform
    # list type [] in [[]] for further compatibility)
    if _np.size(std) != 0:
        # If std is of type [] -- one line with 6 elements (5 collums + 6th element for data type):
        if type(std[0]) != _np.ndarray and _np.size(std) == 8:
            std = std.reshape(-1,8)
        elif (type(std[0]) == _np.ndarray and _np.size(std[0]) != 8) \
                                or (type(std[0]) != _np.ndarray and _np.size(std) != 8):
            verif = False
            loglines += '# ERROR: polt.chkStdLog() not runned! Incompatible number '+ \
                                                            'of collumns in `std.dat`.\n'

    foundstd = False
    for stdi in std:
    # If stdi is reduced and is not to use ('Error' flag)
        if stdi[6] == 'E' and 'no-red' not in stdi[7]:
            continue
        fst = stdi[3]
        calcst = float(stdi[4])
        if f == fst and abs(calc-calcst) < delta:
            foundstd = True
            break
    if not foundstd and verbose:
        print(('# WARNING: Standard star not found for filt. {0} and '+\
            'calc. {1:.1f}\n').format(f, calc))

    if loglines != '':
        writeLog(path, loglines)

    return foundstd



#################################################
#################################################
#################################################
def writeLog(path, strin):
    
    f0 = open('{0}/polt.log'.format(path), 'a')
    f0.writelines(strin)
    f0.close()

    return



#################################################
#################################################
#################################################
def genObjLog(path, subdirs, tgts, sigtol=lambda sigm: 1.1*sigm + 0.00005, \
                    autochoose=False, delta=2.5):
    """
    Generate Objects Log.
    Must be runned after genStdLog, because uses its results.

    delta is the allowed variation for the angles between
    the two beams for one same calcite.
    """

    if len(tgts) != len(subdirs):
        print('\n# ERROR: polt.genObjLog() NOT RUNNED! (len(tgts) != len(subdirs))')
        writeLog(path, '# ERROR: polt.genObjLog() NOT RUNNED! (len(tgts) != len(subdirs))\n')
        return 1

    loglines = ''
    lines=''
    i=0
    maxsize = 7 

    # Loop on list of standards
    for i in range(len(tgts)):
        obj = tgts[i]
        objdir = subdirs[i]
        if obj == '':
            continue
#        print('\n\n**{0}**: outfile(s) with least error '.format(obj.upper()))
        for f in filters:
            if autochoose:
                outs, tests, logs = chooseout('{0}/{1}'.format(path,objdir), obj, f, sigtol=sigtol)
                tags=None
                flag=None
            else:
                outs, tests, tags, flag, logs = queryout(obj, '{0}/{1}'.format(path,objdir), f, sigtol=sigtol)
            loglines += logs

            logs=''
            # Loop on outfiles
            for j in range(len(outs)):
                if outs[j] != '':
                    [Q,U,sig,P,th,sigT,ap,star,MJD,calc],loglixo = readoutMJD(outs[j],obj=obj,verbose=False)
                    # Verify if there are standard star and assigns the 3rd element of tests list
                    tests[j][3] = not chkStdLog(f, calc, path=path, delta=delta, verbose=False)
                    if tests[j][3] == True:
                        loglines += ('# WARNING: {0}_{1}: Standard star not found '+\
                                        '(calc. {2:.1f})\n').format(obj, f, calc)
                    tagstr, flagout = readTests(tests[j], tags=tags, flag=flag)
                    lines += '{0:12.6f} {1:>7s} {2:>10s} {3:>4s} {4:>5.1f} {5:<s} {6:>5s}  {7:<s}\n'.\
                                format(MJD, ccd, obj, f, float(calc), ':::'+_os.path.relpath(outs[j], path)+':::', flagout, tagstr)
                    if len(_os.path.relpath(outs[j], path)) > maxsize:
                        maxsize = len(_os.path.relpath(outs[j], path))

            # Only the last 'logs' (concerning to the outs[-1]), not each one
#            if logs != '':
#                loglines += logs

            # Case all components of outs are the void string '' (object not reduced)
            # tests != [[]] is to filter objects+filter not observed
            if outs.count('') == len(outs) and tests != [[]]:
                tagstr, flagout = readTests([True] + [False]*(len(dictests)-1), tags=tags, flag=flag)
                lines += '{0:12.6f} {1:>7s} {2:>10s} {3:>4s} {4:>5.1f} {5:<s} {6:>5s}  {7:<s}\n'.\
                        format(-1, ccd, obj, f, -1, ':::'+objdir+'/-----:::', flagout, tagstr)
                if len(objdir)+6 > maxsize:
                    maxsize = len(objdir)+6

    if lines == '':
        loglines += '# WARNING: No valid targets were found by polt.genObjLog().\n'
    else:
        lines = '{:12s} {:>7s} {:>10s} {:4s} {:>5s} {:<s} {:>5s}  {:<s}\n'.format('#MJD','ccd',\
        'target','filt','calc',':::outfile:::','flag','tags')+lines

    # Realign columns
    lines_tmp = lines.split('\n')
    lines = ''
    for line in lines_tmp:
        if line == '':
            continue
#        print line
        line = line.split(':::')
#        print line
        lines += line[0][:]
        lines += line[1][:].rjust(maxsize+2)
        lines += line[2][:]
        lines += '\n'

    f0 = open('{0}/obj.dat'.format(path), 'w')
    f0.writelines(lines)
    f0.close()

    writeLog(path, loglines)

    return 0



#################################################
#################################################
#################################################
def genAllLog(path=None, delta=2.5, sigtol=lambda sigm: 1.1*sigm + 0.00005, autochoose=False):
    """
    Generate the std.dat/obj.dat for one reduced night
    
    path: path of night
    delta: allowed variation for the angles between the two
           beams for one same calcite.
    sigtol: tolerance to use the outfiles with all WP instead
            the out with best error. Must be a 'function' that
            receives a sigma value and return the maximum
            sigma for which ignore the best out. Python's
            lambda function can be usefull (e.g., default is
            sigtol=lambda sig: 1.1*sig + 0.00005).
    autochoose: choose best outfiles automatically, without
                interaction?

    """

    if path == None:
        path = _os.getcwd()
    try:
        _os.remove(path+'/polt.log')
    except:
        pass
    
    # Verifies if any log exists
    if _os.path.exists('{0}/std.dat'.format(path)) or _os.path.exists('{0}/obj.dat'.format(path)):
        
        print('\nCAUTION: obj.dat and/or std.dat file already exists for {0} night. '.\
        format(path) + 'Are you sure to continue, overwriting all manual' + \
        'changes on these files? (y/n):')
        while True:
            verif = raw_input('  ')
            if verif in ['y', 'Y']:
                for arq in (path+'/obj.dat', path+'/std.dat', path+'/polt.log'):
                    try:
                        _os.rename(arq, arq+'.back')
                    except:
                        try:
                            _os.unlink(arq+'.back')
                        except:
                            pass
                break
            elif verif in ['n', 'N']:
                print('\nAborted!')
                return
            else:
               print('Value not valid. Please, type \'y\' ou \'n\':')

    # Generates lists
    ltgts = _np.loadtxt('{0}/refs/pol_alvos.txt'.format(_hdtpath()), dtype=str)
    lstds = _np.loadtxt('{0}/refs/pol_padroes.txt'.format(_hdtpath()), dtype=str,\
            usecols=[0])
    subdirs = [fld for fld in _os.listdir('{0}'.format(path)) if \
            _os.path.isdir(_os.path.join('{0}'.format(path), fld))]
    tgts=[]    # variable for real target names, not the subdirectory names
    stds=[]    # variable for real standard names, not the subdirectory names
    lines = ''
    loglines = ''

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
        while obj not in _np.hstack((ltgts,lstds,_np.array(['calib','']))):
            print('\nObject {0} is not a known target or standard!!'.format(obj))
            obj = raw_input('Type the common name or press -Enter- to ignore it: ')
        if obj in lstds:
            tgts.append('')
            stds.append(obj)
        elif obj in ltgts:
            tgts.append(obj)
            stds.append('')
        else:
            if obj != 'calib':
                loglines += ('# CAUTION: Object inside {0} directory '.format(obj_curr) + \
                                'was left out of obj.dat/std.dat files.\n')
            tgts.append('')
            stds.append('')

#    print subdirs
#    print tgts
#    print stds
    print('')
    writeLog(path, loglines)
    genStdLog(path, subdirs, stds, sigtol=sigtol, autochoose=autochoose)
    genObjLog(path, subdirs, tgts, delta=delta, sigtol=sigtol, autochoose=autochoose)

    # Write user name and date+time
    username = _pwd.getpwuid(_os.getuid())[4]
    if username.find != -1:
        username = username[:username.find(',')]
    loglines = _time.strftime("\nGenerated at: %Y/%m/%d - %I:%M %p\n")
    loglines += '          by: ' + _pwd.getpwuid(_os.getuid())[0] + ' (' + username + ')'
    writeLog(path, loglines)
    
    with open('{0}/polt.log'.format(path), 'r') as fl:
        print('\n{0}\nPOLTOOLS LOG (content of {1}/polt.log)\n\n{2}'.format('='*40, path, fl.read()))

    return

    

#################################################
#################################################
#################################################
def corObjStd(night, f, calc, path=None, delta=2.5):
    """
    Correlate a target observer at one filter and calcite
    (`f` and `calc`) and return the values for matching
    standard stars inside `night`'s std.dat of those which
    are not marked with `E` flag.

    delta: tolerance, in degree, of angle of two beams for
           one same calcite (default: +/-2.5 degree)

    Output, in order: stdname, thstd, angref, flagstd
      - stdname: list with standard names
      - thstd: list with theta measured values (the angle
               returned is the same from .out, NOT 180-theta!)
      - angref: list with theta published values
      - flastd: list with flags concerning to the data.
    """
    if path == None:
        path = _os.getcwd()
    try:
        std = _np.loadtxt('{0}/refs/pol_padroes.txt'.format(_hdtpath()), dtype=str)
    except:
        print('ERROR: Can\'t read files pyhdust/refs/pol_padroes.txt')
        raise SystemExit(1)

    calc = float(calc)
    angref = []
    thstd = []
    stdname = []
    flagstd = []
    tagstd = []

    if _os.path.exists('{0}/{1}/std.dat'.format(path,night)):
        stds = _np.loadtxt('{0}/{1}/std.dat'.format(path,night), dtype=str)
        if len(stds) > 0 and len(stds[-1]) != 8:
            stds = stds.reshape(-1,8)
        k = 0
        sigs = []
        for stdinf in stds:
            if stdinf[6] == 'E' and 'no-red' in stdinf[7] and stdinf[3] == f:
                print('WARNING: {0} night: standard {1}_{2} not reduced and unknown calcite.'.format(night,stdinf[2],f))
                continue
            elif stdinf[6] == 'E' and stdchk(stdinf[2])[0] and stdinf[3] == f and abs(float(stdinf[4])-calc) <= delta:
                print('WARNING: {0} night: standard {1}_{2} wasn\'t used because has `E` flag.'.format(night,stdinf[2],f))
                continue
            elif stdinf[6] != 'E' and stdchk(stdinf[2])[0] and stdinf[3] == f and abs(float(stdinf[4])-calc) <= delta:
                # Bednarski: Show error message now
                try:
                    Q, U, sig, P, th, sigT, tmp, tmp2 = readout('{0}/{1}'.\
                                            format(path+'/'+night,stdinf[5]))
                except:
                    print('WARNING: {0} night: standard {1}_{2} wasn\'t used because can\'t open/read {3}.'.format(night, stdinf(2), f, stdinf[5]))
                    continue
                stdname += [ stdinf[2] ]
                thstd += [ float(th) ]
                flagstd += [stdinf[6]]
                tagstd += [stdinf[7]]
                if thstd == 0.:
                    thstd = 0.01
                sigs += [ float(sig)*100 ]
                sigth = 28.65*sig/P 
                i = stdchk(stdinf[2])[1]
                j = filters.index(f)+1 #+1 devido aos nomes na 1a coluna
                angref += [ float(std[i,j]) ]
        if stdname == []:
            print('# WARNING: {0} night: no standard found for filt. {1} and calc. {2:.1f}.'.format(night, f, calc))
    else:
        print('# ERROR: {0} night: `std.dat` not found. Ignoring its standards...'.format(night))
#    print stdname, thstd, angref

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
        * th_eq = -th_measured + th^std_measured + th^std_published
        * th_eq = 180 - th_measured - dth
        * dth = 180 - th^std_measured - th^std_published
        * Q_eq = P*cos(2*th_eq*pi/180)
        * U_eq = P*sin(2*th_eq*pi/180)
        * sigth = 28.65*sig/P
    """

    if path == None or path == '.':
        path = _os.getcwd()

    # Carregar angulos de referencia para todos os filtros de um conjunto de padroes
    if PAref is not None:
        for line in PAref:
            if len(line) < 21:
                print('# ERROR! Wrong PAref matrix format')
                raise SystemExit(1)
    else:
        PAref = _np.loadtxt('{0}/refs/pol_padroes.txt'.format(_hdtpath()), dtype=str)

    # Read lists and verify if target is a valid target

    try:
        obj = _np.loadtxt('{0}/refs/pol_alvos.txt'.format(_hdtpath()), dtype=str)
        std = _np.loadtxt('{0}/refs/pol_padroes.txt'.format(_hdtpath()), dtype=str)
    except:
        print('ERROR: Can\'t read files pyhdust/refs/pol_alvos.txt and/or pyhdust/refs/pol_padroes.txt.')
        raise SystemExit(1)
        
    if target not in _np.hstack((std[:,0],obj)):
        print('# Warning: Target {0} is not a default target or standard!!'.\
        format(target))
        tmp = raw_input('# Type something to continue...')

    nights = [fld for fld in _os.listdir(path) if _os.path.isdir(_os.path.join(path, fld))]
    lines = ''

    for night in nights:

        # Check obj.dat for the night
        if _os.path.exists('{0}/{1}/obj.dat'.format(path,night)):
            try:
                objs = _np.loadtxt('{0}/{1}/obj.dat'.format(path,night), dtype=str)
            except:
                print('# ERROR! {0} night: Can\'t read obj.dat file. Ignoring it...'.format(night))
                continue
#                raise SystemExit(1)

            # Verify if std has more than one line. Case not, do the reshape
            if _np.size(objs) == 8:
                objs = objs.reshape(-1,8)
            elif _np.size(objs) % 8 != 0:
                print('# ERROR! {0} night: Wrong column type in obj.dat file. Ignoring it...'.format(night))
                continue
#                raise SystemExit(1)

            # Loop on found nights
            for objinf in objs:
                dth = []
                stdnames = ''

                if objinf[2] == target:
                    MJD, ccd, obj, f, calc, out, flag, tags = objinf
                    if flag == 'E':
                        print('Warning: {0} night: {1}_{2} with `E` flag and tags `{3}`. Ignoring it...'.format(night,obj,f,tags))
                        continue
                    try:
                        Q, U, sig, P, th, sigT, tmp, tmp2 = readout('{0}/{1}'.\
                                                    format(path+'/'+night,out))
                    except:
                        print('ERROR: {0} night: Can\'t open/read out file {1}. Ignoring it...'.format(night,out))
                        continue

                    P = float(P)*100
                    th = float(th)
                    sig = float(sig)*100
                    sigth = 28.65*sig/P
                    if 'no-std' not in tags:
                        stdname, thstd, angref, flagstd, tagstd = corObjStd(night, f, calc, path=path, delta=delta)
                    else:
#                        FURTHER: APPLY ALTERNATIVE METHOD TO COMPUTE DTHETA HERE
                        pass

                    # Bednarski: working for more than one standard star
                    if 'no-std' in tags or stdname==[]:
                        mdth = 0.
                        devth = 0.
                        stdnames = '---'
                        thstd=[]
                        angref=[]
                        flagstd=[]
                        tagstd=[]
                        # Case, for some reason, target observation was not marked with 'no-std' tag and flag=='OK'
                        if 'no-std' not in tags:
                            tags += ',no-std'
                        if flag == 'OK':
                            flag = 'W'
                        ##print('# ERROR with thstd ou angref!')
                        ##print(night, f, calc)
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
                            dth += [ -thstd[i]-angref[i] ]
                            while dth[i] >= 180:
                                dth[i]-= 180
                            while dth[i] < 0:
                                dth[i]+= 180
                            if i != 0:
                                stdnames += ','
                            stdnames += stdname[i]
                        mdth=sum(dth)/len(dth)  # evalute the mean dth
                        devth=_np.std(dth)
                    th = -th-mdth
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
                        '{:>6.2f} {:>13s} {:>5s} {:>s}\n').format(MJD, night, ccd, f, calc, stdnames, \
                                    mdth, devth, P, Q, U, th, sig, sigth, outn, flag, tags)
        else:
            print('# Warning: {0} night: `obj.dat` not found. Ignoring it...'.format(night))
    #arquivo de saida
    if lines != '':
        print('\n# Output written: {0}/{1}.log'.format(path,target))
        f0 = open('{0}/{1}.log'.format(path,target),'w')
        lines = ('#MJD            night     ccd filt  calc     stdstars    dth devdth' +\
                '        P        Q        U      th    sigP  sigth       outfile flag   tags\n')+lines
        f0.writelines(lines)
        f0.close()
    else:
        print('\n# ERROR! No valid observation was found for target {0}'.format(target))
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
    if path == None:
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
def sortLog(file):
    """ Sort the *.out file """
    f0 = open(file)
    lines = f0.readlines()
    f0.close()
    log = _np.loadtxt(file, dtype=str)
    log = log[log[:,0].argsort()]
    fmt = '%12s %7s %1s %5s %5s %6s %5s %6s %6s %6s %5s %5s'
    _np.savetxt(file.replace('.log','.txt'), log, fmt=fmt, header=lines[0])
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
        file = _np.loadtxt(fname, dtype=str, delimiter=' ')
    except:
        print('Error when reading the file or file doesn\'t exist!')
        return
    flist = file.tolist()
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

