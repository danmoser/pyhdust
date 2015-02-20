##!/usr/bin/env python
#-*- coding:utf-8 -*-
#Modified by D. Moser in 2014-10-17

"""
HDUST tools

includes *stars, *filters
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from glob import glob
import os, time
import pyhdust.phc as phc
import pyhdust.jdcal as jdcal

def n0toSigma0(n0,M,Req,f,Tp,mu):
    """ From n0 to sigma0 """
    rho0 = n0*mu*phc.mH.cgs
    a = (phc.kB.cgs*f*Tp/mu/phc.mH.cgs)**.5
    sig0 = (2*np.pi)**.5*a/(phc.G.cgs*M/Req)**.5*Req*rho0
    return sig0

def n0toMdot(n0,M,Req,f,Tp,mu,alpha,R0):
    """ From n0 to Mdot """
    rho0 = n0*mu*phc.mH.cgs
    a = (phc.kB.cgs*f*Tp/mu/phc.mH.cgs)**.5
    Mdot = 3*np.pi*(2*np.pi)**.5*alpha*a**3./(phc.G.cgs*M/Req)*rho0*Req**2.*((R0/Req)**.5-1)**-1.
    return Mdot/phc.Msun.cgs*phc.yr.cgs

def plotMJDdates(spec=None, pol=None, interf=None, limits=None):
    """
    Plot dates from spec (Class), pol (routines) and interf (ESO query)

    This need to be polished !!!!
    """
    fig, ax = plt.subplots()
    spec = 'data_aeri_splot.txt'
    if spec is not None:
        spJD = np.loadtxt(spec)
        spJD = spJD[:,0]
        y = [ 0. for JD in spJD ]
        ax.plot(spJD, y, marker='d', color='lightgray', ls='')
        #yerr = [ [1. for JD in spJD], [1. for JD in spJD] ]
        #ax.errorbar(spJD, y, yerr, marker='o', color='blue', ls='')

    pol = 'pol_aeri.log'
    if pol is not None:
        polJD = np.loadtxt(pol, dtype=str)
        polJD = polJD[:,9]
        polJD = np.array(polJD, dtype=float)-2400000.5
        y = [ -.5 for JD in polJD ]
        ax.plot(polJD, y, marker='o', color='gray', ls='')
        #yerr = [ [.5 for JD in polJD], [1.5 for JD in polJD] ]
        #ax.errorbar(polJD, y, yerr, marker='x', color='green', ls='')

    interf = 'interf_aeri.txt'
    if interf is not None:
        intJD = np.loadtxt(interf, dtype=str, delimiter=',', skiprows=1)
        intJD = np.array(intJD[:,-2], dtype=float)
        y = [ .5 for JD in intJD ]
        ax.plot(intJD, y, marker='s', color='darkgrey', ls='')
        #yerr = [ [1.5 for JD in intJD], [.5 for JD in intJD] ]
        #ax.errorbar(intJD, y, yerr, marker='s', color='red', ls='')

    limits = (56100., 56750.)
    if limits is None:
        mjd0, mjd1 = ax.get_xlim()
    else:
        mjd0, mjd1 = limits
        ax.set_xlim(limits)
    ticks = phc.gentkdates(mjd0, mjd1, 3, 'm', dtstart=dt.datetime(2012,7,1).\
    date())
    mjdticks = [jdcal.gcal2jd(date.year,date.month,date.day)[1] for date in \
    ticks]
    ax2 = ax.twiny()
    ax2.set_xlim(limits)
    ax2.set_xticks(mjdticks)
    ax2.set_xlabel('Civil date')
    ax2.set_xticklabels([date.strftime("%d %b %y") for date in ticks])
    plt.setp( ax2.xaxis.get_majorticklabels(), rotation=45 )
    ax.set_yticklabels([])
    ax.set_xlabel('MJD')
    return
    

def hdtpath():
    """
    Return the path os the module.

    >>> hdt.hdtpath()
    /home/user/Scripts/pyhdust/
    """
    return __file__[:__file__.rfind('/')+1]

def doFilterConv(x0, y0, filter):
    """
    Return the convolved filter total flux for a given flux y0,
    at wavelengths x0.
    """
    fdat = np.loadtxt('{}/filters/{}.dat'.format(hdtpath(),filter.lower()),\
    skiprows=1)
    fdat[:,0] /= 10000 #from Angs to microns
    #interpfunc = interpolate.interp1d(fdat[:,0], fdat[:,1], kind='linear')
    interpfunc = interpolate.InterpolatedUnivariateSpline(fdat[:,0], fdat[:,1])

    idx = np.where( (x0 >= fdat[0,0]) & (x0 <= fdat[-1,0]) )
    x0 = x0[idx]
    y0 = y0[idx]
    y = interpfunc(x0)*y0/np.sum( interpfunc(x0) )

    return np.sum(y)

def doPlotFilter(file, obs, filter, sed2data, pol=False):
    """
    doPlotFilter
    """
    x0 = sed2data[j,:,2]
    if pol:
        y0 = sed2data[j,:,7]
        savename = 'pol_{}_{}_{}'.format(file[file.rfind('/')+1:-5],obs,\
        filter)
    else:
        y0 = sed2data[j,:,3]
        savename = '{}_{}_{}'.format(file[file.rfind('/')+1:-5],obs,filter)

    fdat = np.loadtxt('filters/{}.dat'.format(filter),\
    skiprows=1)
    fdat[:,0] /= 10000
    #interpfunc = interpolate.interp1d(fdat[:,0], fdat[:,1], kind='linear')
    interpfunc = interpolate.InterpolatedUnivariateSpline(fdat[:,0], fdat[:,1])

    idx = np.where( (x0 >= fdat[0,0]) & (x0 <= fdat[-1,0]) )
    x0 = x0[idx]
    y0 = y0[idx]
    y = interpfunc(x0)*y0/np.sum( interpfunc(x0) )

    fig, ax = plt.subplots()
    ax.plot(x0, y0, label='SED')
    #ax.plot(fdat[:,0], fdat[:,1], label='Filter')
    ax.plot(x0, interpfunc(x0)*y0, label='Convolved')
    ax.set_title(savename)
    ax.legend()
    fig.savefig(savename+'.png', transparent=True)
    fig.clf()
    return    

def sed2info(file):
    """
    Read info from SED2 file.
    input = file
    output = nlbd, nobs, Rstar, Rwind (as floats)
    """
    f0 = open(file, 'r')
    fcont = f0.readlines()
    f0.close()
    info = ''
    i=0
    while info == '':
        if fcont[i][0] != '%':
            info = np.array(fcont[i].split(), dtype=float)
        i+=1
    return info

def readfullsed2(file):
    """
    Read data from FULLSED2 file.
    input = file
    output = ndarray((nobs,nlbd,-1))
        number of columns from SED2file replaces "-1"
    """
    nlbd, nobs, Rstar, Rwind = sed2info(file)
    sed2data = np.loadtxt(file, skiprows=5)
    sed2data = sed2data.reshape((nobs,nlbd,-1))
    return sed2data

def readsed2(file):
    """
    Read data from SED2 file.
    input = file
    output = ndarray((nobs,nlbd,-1))
        number of columns from SED2file replaces "-1"
    """
    nlbd, nobs, Rstar, Rwind = sed2info(file)
    sed2data = np.loadtxt(file, skiprows=1)
    #sed2data = sed2data.reshape((nobs,nlbd,-1))
    return sed2data

def chkObsLog(path=None, nights=None, badweath=None):
    """ Check if there is data for all nights with observations.

    If not, check if the night is in the list of night lost due to bad weather.

    If no data and no bad weather info is registered, prints an error.

    If the night is included as bad weather and is not in night list, prints a
    warning.
    """
    if path == None:
        path = os.getcwd()
    if nights == None:
        nights = '{0}/refs/noites.txt'.format(hdtpath())
    lnights = np.loadtxt(nights, dtype=str)
    if badweath == None:
        badweath = '{0}/refs/maltempo.txt'.format(hdtpath())
    lbadweath = np.loadtxt(badweath, dtype=str)
    for night in lnights:
        if night in lbadweath:
            pass
        elif os.path.exists(night) == False:
            print('# ERROR! {0} has no data and was not lost for bad weather!'.\
            format(night))
    flds = [fld for fld in os.listdir('{0}'.format(path)) if \
    os.path.isdir(os.path.join('{0}'.format(path), fld))]
    for fld in flds:
        if fld not in lnights:
            print('# Warning! Night {0} is not recorded as OPD night!'.\
            format(fld))
            print('# Update the file {0}'.format(nights))
    for night in lbadweath:
        if night not in lnights:
            print('# Warning! Bad weather {0} is not recorded as OPD night!'.\
            format(night))
            print('# Probably it is a spec night.')            
    return

def mergesed2(models, Vrots, path=None):
    """
    Merge all mod#/*.sed2 files into the fullsed file.
    
    It will check if all sed2info are the same (i.e., nlbd, nobs, Rstar, Rwind).
    If not, it ask if you want to continue (and receive and error).
    
    The presence of the SED file is not required anymore.
    The structure is set by the first broadband sed2 found.
    
    NO AVERAGE is coded.

    Input: *.txt ou *.inp lists.
    Output: files written.
    Print status.
    """

    sufbands = ['SED', 'UV', 'IR', 'NIR', 'BALMER', 'PASCHEN', 'CM', 'MM', #0-7
    'J', 'H', 'K', 'L', 'M', 'N', 'Q1', 'Q2'] #8-14
    #wavelength in microns
    suflines = {'Hb':.486268, 'Ha':.656461, 'Br15':1.641, 'Brg':2.166}
    
    for model in models:
        modfld, modelname = phc.trimpathname(model)
        path = phc.trimpathname(modfld[:-1])[0]
        if os.path.exists('{0}/fullsed'.format(path)) == False:
            os.system('mkdir {0}/fullsed'.format(path))
        #
        modelname = modelname.replace('.txt','.inp')
        sed2data = np.empty(0)
        sfound = []
        #Process broad-bands
        for suf in sufbands:
            file = modfld+'{}_{}.sed2'.format(suf,modelname.replace(".inp",""))
            if os.path.exists(file):
                sfound += [suf]
                newdata = readsed2(file)
                if len(sed2data) == 0:
                    sed2data = newdata.copy()
                    nlbd, nobs, Rstar, Rwind = sed2info(file)
                else:
                    #Check if the SED2 file has the info as the first file
                    if np.product( (nobs, Rstar, Rwind) == sed2info(file)[1:] )==0:
                        key = ''
                        while key.upper() != Y:
                            print('# WARNING: {} has different HDUST output!!!'.\
                            format(modelname))
                            key = raw_input('Do you want do proceed? (y/other): ')
                    nlbd += sed2info(file)[0]
                    sed2data = np.vstack((sed2data,newdata))
        #Process lines
        Vrot = Vrots[models.index(model)]
        for suf in suflines:
            file = modfld+'{}_{}_SEI.sed2'.format(suf,modelname.replace(".inp",""))
            if os.path.exists(file):
                sfound += [suf]
                newdata = readsed2(file)
                #print("# TRIMMING INPUT SPECTRUM {}".format(modelname))
                #print("# TO ACCOUNT FOR THE NON-ZERO ROTATION VELOCITY.")
                deltalbd = Vrot/phc.c.cgs/1e-5*suflines[suf]
                mini = newdata[0,2]+deltalbd
                maxi = newdata[-1,2]-deltalbd
                #print deltalbd, Vrot
                idx = np.where((newdata[:,2] >= mini) & (newdata[:,2] <= maxi))
                #print len(newdata), len(newdata[idx])
                ncut = len(newdata)-len(newdata[idx])
                newdata = newdata[idx]
                if len(sed2data) == 0:
                    sed2data = newdata.copy()
                    nlbd, nobs, Rstar, Rwind = sed2info(file)
                else:
                    #Check if the SED2 file has the info as the first file
                    if np.product((nobs, Rstar, Rwind) == sed2info(file)[1:]) == 0:
                        key = ''
                        while key.upper() != Y:
                            print('# WARNING: {} has different HDUST input!!!'.\
                            format(modelname))
                            key = raw_input('Do you want do proceed? (y/other): ')
                    idx = np.where((sed2data[:,2] < mini) | (sed2data[:,2] > maxi))
                    ncut += len(sed2data)-len(sed2data[idx])
                    sed2data = sed2data[idx]
                    nlbd += sed2info(file)[0]-ncut/sed2info(file)[1]
                    sed2data = np.vstack((sed2data,newdata))

        if len(sfound) > 0:
            print('# PROCESSED: {} with {}'.format(model, sfound))
            fullsed2 = np.zeros((len(sed2data),16))
            fullsed2[:,0:3+1] = sed2data[:,0:3+1]
            fullsed2[:,4] = sed2data[:,11]
            fullsed2[:,5] = sed2data[:,19]
            fullsed2[:,6] = sed2data[:,27]
            fullsed2[:,7:8+1] = sed2data[:,4:5+1]

            #a = np.arange(9).reshape(3,3)
            #np.core.records.fromarrays(a.transpose(), names='a,b,c', formats='f4,f4,f4')
            #fullsed2 = fullsed2[fullsed2[:,2].argsort()]
            #fullsed2 = fullsed2[fullsed2[:,0].argsort()]
            fullsed2 = np.core.records.fromarrays(fullsed2.transpose(), names=\
            'MU,PHI,LAMBDA,FLUX,SCT FLUX,EMIT FLUX,TRANS FLUX,Q,U,Sig FLUX,\
            Sig SCT FLUX,Sig EMIT FLUX,Sig TRANS FLU,Sig Q,Sig U',
            formats='f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8')
            idx = np.argsort(fullsed2, order=('MU','LAMBDA'))
            fullsed2 = fullsed2[idx]
        
            hd = '%CONTAINS: {}\n'.format(' + '.join(sfound))
            hd+= '%CREATED: {}\n'.format(time.asctime( time.localtime(time.time()) ))
            hd+= '%{:>7s}{:>8s}{:>13s}{:>13s}'.format('nlbd','nobs','Rstar','Rwind')+'\n'
            hd+= '{:8d}{:8d}{:13.4f}{:13.2f}\n'.format(int(nlbd),int(nobs),Rstar,Rwind)
            hd+= '%{:>12s}'.format('MU')+(15*'{:>13s}').format('PHI','LAMBDA','FLUX',\
            'SCT FLUX','EMIT FLUX','TRANS FLUX','Q','U','Sig FLUX','Sig FLUX',\
            'Sig SCT FLUX','Sig EMIT FLUX','Sig TRANS FLUX','Sig Q','Sig U')
        
            np.savetxt(path+'/fullsed/fullsed_'+modelname.replace('.inp','.sed2'),\
            fullsed2, header=hd, comments="", fmt='%13.6f', delimiter='')
        else:
            print('# WARNING: No SED2 found for {}'.format(model))
    return


### STARS ###
def calcTeff(Lum, size, M=None):
    """
    Calculate Teff for the non-rotating case. `size` variable is assumed to be
    the stellar radius (M==None). If M is given, size is assumed to be log(g).

    log(g) in cgs units.

    Lum, Radius and Mass in Solar units.
    """
    #~ M, Rp, Lum = fundline
    if M == None:
        Rp = size*phc.Rsun.cgs
    else:
        Rp = (M*phc.Msun.cgs*phc.G.cgs/10**size)**.5
    L = Lum*phc.Lsun.cgs
    #Lum = 4*np.pi*Rp**2*phc.sigma*Teff**4
    Teff = (L/(4*np.pi*Rp**2*phc.sigma.cgs))**.25
    return Teff


def calclogg(fundline):
    """
    STARS:
    Calculate Teff
    """
    M, Rp, Lum = fundline
    Rp = Rp*phc.Rsun
    logg =np.log10(phc.G*M*phc.Msun/Rp**2)
    return logg

def genlog(path=None, extrainfo=None):
    """Gen. log of the calculated models of the project.

    ppath = Project's path. If it is not given, it assumes the local pwd.

    extrainfo = {\
    'mod01':'i=60',\
    'mod02':'i=60+70',\
    'mod03':'i=60+70',\
    'mod04':'i=60+70/-source NO ROT',\
    'mod05':'i=60+70/?',\
    'mod06':'i=60+70/?'}
    """
    if path == None:
        path = os.getcwd()
    modfld = glob('{0}/mod*'.format(path))
    #while len(modfld) == 0:
    #    proj = raw_input('Type the project name: ')
    #    modfld = glob('{0}/mod*'.format(proj))
    modfld.sort()

    #MODN, steps, sed2, maps, extrainfo
    tab = np.zeros((5000,5), dtype='|S127')
    i = 0
    
    for modn in modfld:
        modnn = phc.trimpathname(modn)[1]
        mods = glob('{0}/{1}_*.txt'.format(modn,modnn))
        print('# Catalogue of {0}'.format(modn))
        for mod in mods:
            #if mod.find('aeri') > 0:
                suf = mod[mod.rfind('/')+1:-4]
        
                step1 = glob('{0}/{1}??.temp'.format(modn,suf))
                if len(step1) == 0:
                    step1=['0']
                step1.sort()
        
                sed2 = glob('{0}/*{1}_*.sed2'.format(modn,suf))
                sed2+= glob('{0}/*{1}.sed2'.format(modn,suf))
                sed2.sort()
                s2out = ''
                for sedi in sed2:
                    s2out+=sedi[sedi.rfind('/')+1:sedi.find('_')]+'+'
                
                maps = glob('{0}/*{1}_*.maps'.format(modn,suf))
                maps+= glob('{0}/*{1}.maps'.format(modn,suf))
                maps.sort()
                mout = ''
                for mapi in maps:
                    mout+=mapi[mapi.rfind('/')+1:mapi.find('_')]+'+'

                if extrainfo != None:
                    if modnn in extrainfo:
                        extra = extrainfo[modnn]
                    else:
                        extra = ''
    
                tab[i] = (suf, step1[-1][-7:-5], s2out[:-1], mout[:-1], extra)
                i+=1
    
        if len(mods) == 0:
            print('# NO model found in {0}'.format(modn))
    
    #tab = tab[:i,:]
    #tab = tab[tab[:,0].argsort()]
    #np.savetxt('log.csv', tab, fmt='%s', delimiter=',')
    np.savetxt('{0}/log.csv'.format(path), tab[:i,:], fmt='%s', delimiter=',')
    #tab = np.loadtxt('log.csv', dtype=str)
    #tab = tab[tab[:,0].argsort()]
    #np.savetxt('log.csv', tab, fmt='%s', delimiter=',')
    return

def printN0(modn):
    """
    modn = '02'
    
    """
    import matplotlib.pyplot as plt

    #generate and print n_0
    os.system('grep n_0 mod{0}/*PLn3.5*.txt > n0_mod{0}.txt'.format(modn))
    #n0file = np.genfromtxt('n0_mod01.txt', delimiter = (6,3,3,4,4,2,3,3,5,\
    #5,5,3,4), usecols=(4,10,12), dtype=None) 
    #cols = (0#, 1type, 2typeVal, 3#, 4sig, 5#, 6h, 7#, 8Rd, 9#, 10M, 11#,
    # 12ob, 13#)
    n0file = np.genfromtxt('n0_mod{}.txt'.format(modn), delimiter=\
    (22,4,18,5,3,4,35,8), usecols=(1,3,5,7), dtype=None) 
    #cols = (sig0, M, ob, n0)
    
    sig0s = []
    for item in n0file[:,0]:
        if item not in sig0s:
            sig0s += [item]
    
    Ms = []
    for item in n0file[:,1]:
        if item not in Ms:
            Ms += [item]
    
    obs = [1.1,1.2,1.3,1.4,1.45]
    
    plt.clf()
    for item in obs:
        idx = np.where(n0file[:,2] == item)
        plt.title('Sig0 range: [{}, {}] (g cm-2)'.format(sig0s[0], sig0s[-1]))    
        plt.plot(n0file[:,0][idx], n0file[:,3][idx], 'o')
        plt.yscale('log')
        plt.xscale('log')
        plt.ylabel('n0 (# cm-3)')
        plt.xlabel('sig0 (g cm-2)')
    plt.savefig('n0_vs_sig0.png')
    plt.savefig('n0_vs_sig0.pdf')
    
    Btp = ['B0.5','B1','B1.5','B2','B2.5','B3','B4','B5','B6','B7','B8','B9']
    plt.clf()
    for item in obs:
        idx = np.where(n0file[:,2] == item)
        plt.title('Sig0 range: [{}, {}] (g cm-2)'.format(sig0s[0], sig0s[-1]))
        plt.plot(n0file[:,1][idx], n0file[:,3][idx], 'o')
        plt.plot([3.4,14.6], [1e13,1e14], 'r-')
        plt.plot([3.4,14.6], [1e12,1e12], 'r-')
        plt.xlim([16,2])
        plt.yscale('log')
        plt.ylabel('n0 (# cm-3)')
        plt.xlabel('mass (Solar units)')
    for i in range(len(Btp)):
        plt.text(Ms[::-1][i], 4e14, Btp[i])
    plt.savefig('n0_vs_M.png')
    plt.savefig('n0_vs_M.pdf')  
        
    
    plt.clf()
    for item in obs:
        idx = np.where(n0file[:,2] == item)
        plt.title('Sig0 range: [{}, {}] (g cm-2)'.format(sig0s[0], sig0s[-1]))
        plt.plot(n0file[:,2][idx], n0file[:,3][idx], 'o')
        plt.yscale('log')
        plt.ylabel('n0 (# cm-3)')
        plt.xlabel('radii ratio (Req/Rp)')
        plt.xlim([1.0,1.5])
    plt.savefig('n0_vs_ob.png')
    plt.savefig('n0_vs_ob.pdf')    
    return 

def rmMods(modn, Ms):
    """
    #MODEL NUMBER
    modn        = '02'
    #Masses list ans sig0 POSITION do be excluded
    Ms = [\
    ['14.6', [0]],\
    ['12.5', [0,-1]],\
    ['10.8', [0,-1]],\
    ['09.6',  [0,-2,-1]],\
    ['08.6',  [0,-2,-1]],\
    ['07.7',  [0,-2,-1]],\
    ['06.4',  [0,-3,-2,-1]],\
    ['05.5',  [0,-3,-2,-1]],\
    ['04.8',  [-4,-3,-2,-1]],\
    ['04.2',  [-4,-3,-2,-1]],\
    ['03.8',  [-4,-3,-2,-1]],\
    ['03.4',  [-4,-3,-2,-1]],\
    ]
    """
    #Create sig0 list
    mods = glob('mod{}/*NI*.txt'.format(modn))
    n0file = np.genfromtxt(mods, delimiter = (22,4), usecols=(1), dtype=None)
    sig0s = []
    for item in n0file:
        if item not in sig0s:
            sig0s += [item]
    sig0s.sort()
    
    for item in Ms:
        M = item[0]
        exsig = item[1]
        for rm in exsig:
            os.system('rm mod{0}/mod{0}*_sig{1}*_M{2}*.txt'.format(modn,sig0s[rm],\
            M))
            print('# Deleted mod{0}/mod{0}*_sig{1}*_M{2}*.txt'.format(modn,\
            sig0s[rm],M))
    #End prog
    return 

def splitKurucz(path=None):
    """ Function doc

    @param PARAM: DESCRIPTION
    @return RETURN: DESCRIPTION
    """
    if path == None:
        path = os.getcwd()
    allk = np.loadtxt('ap00k0.dat', dtype=str, delimiter='\n')
    
    for i in range(0,len(allk)-1):
        if 'EFF' in allk[i]:
            iref = i
            teff = int(allk[i].split()[1][:-1])
            logg = float(allk[i].split()[3][:-3])
        elif 'DECK6 72' in allk[i]:
            allk[i] = allk[i].replace('DECK6 72','DECK6 71')
        elif 'EFF' in allk[i+1]:
            np.savetxt('ap00k0tef%05dg%.1f.dat' % (teff,logg), allk[iref:i+1], fmt='%s')
    
    np.savetxt('ap00k0tef%05dg%.1f.dat' % (teff,logg), allk[iref:], fmt='%s')
    return

def diskcalcs(M, R, Tpole, T, alpha, R0, mu, rho0, Rd):
    """ Do the equivalence of disk density for different quantities.

    Note that they all depend of specific stellar quantities!!!

    @param PARAM: DESCRIPTION
    M = 10.3065*Msun    
    R = 7*Rsun
    #R = 5.38462*Rsun
    #R = 2*5.38462*Rsun
    Tpole = 26025.
    #Tpole = 23168.
    #Tpole = 16382.
    T = 0.72*Tpole
    alpha = 1.
    R0 = 1e14*R
    R0 = 0.*R
    mu = 0.5
    rho0 = 5e12 #in particles per cubic centimeter
    rho0 = 2.35e13
    Rd = 18.6*R

    @return RETURN: DESCRIPTION
    """
    def rho2sigp(R,rho0,a,M):
        sig = (2*np.pi)**.5*a/(phc.G.cgs*M/R)**.5*R*rho0
        return sig
        
    def rho2Mdot(R,alpha,a,M,rho0,R0):
        Mdot = 3*np.pi*(2*np.pi)**.5*alpha*a**3./(phc.G.cgs*M/R)*rho0*R**2.*((R0/R)**.5-1)**-1.
        return Mdot
        
    def Mdot2sig(R,Mdot,alpha,a,M,R0):
        sig = Mdot*(phc.G.cgs*M/R)**.5/(3.*np.pi*alpha*a**2*R)*((R0/R)**.5-1)
        return sig

    a = (phc.kB.cgs*T/mu/phc.mH.cgs)**.5
    rho0 = rho0*mu*phc.mH.cgs
    sigp = rho2sigp(R,rho0,a,M)
    rho0p = sigp/(2*np.pi)**.5/a*(phc.G.cgs*M/R)**.5/R
    Mdot = rho2Mdot(R,alpha,a,M,rho0,R0)
    sig = Mdot2sig(R,Mdot,alpha,a,M,R0)
    sigl= Mdot*(phc.G.cgs*M)**.5/3/np.pi
    Mdisk0 = (2*np.pi)**1.5*rho0*R**2.*(Rd-R)*a/(phc.G.cgs*M/R)**.5
    Mdisk = 2*np.pi*Mdot*(phc.G.cgs*M/R)**.5*R**.5/3/np.pi/alpha/a**2.*R0**.5*np.log(Rd/R) 
    if Mdisk == 0:
        Mdisk = 4*np.pi*Mdot*(phc.G.cgs*M/R)**.5*R**.5/3/np.pi/alpha/a**2.*(Rd**.5-R**.5)
    MdiskG= 2*Mdot*(phc.G.cgs*M/R)**.5*R**.5/3/alpha/a**2*(R0**.5*np.log(Rd/R)+2*R**.5-2*Rd**.5)
    
    print('R0/R  = {:.1f}'.format(R0/R))
    print('Valid sigma (1)?: {}'.format(round(sigp/sig)==1))
    print('Valid sigma (2)?: {}'.format(round(rho0/rho0p)==1))
    print('rho0  = {:.2e} g/cm3'.format(rho0))
    print('sigma0= {:.2e} g/cm2'.format(sig/alpha/a**2))
    print('Mdot  = {:.2e} Msun/yr'.format(Mdot/phc.Msun.cgs*phc.yr.cgs))
    print('Mdisk0= {:.2e} Msun [#from rho0]'.format(Mdisk0/phc.Msun.cgs))
    print('Mdisk = {:.2e} Msun [#approx.]'.format(Mdisk/phc.Msun.cgs))
    print('MdiskG= {:.2e} Msun'.format(MdiskG/phc.Msun.cgs))
    print('PS: Mdisk for both isothermal Sigma(r) and H(r)')
    return

### MAIN ###
if __name__ == "__main__":
    pass
