#-*- coding:utf-8 -*-

import sys
import numpy as np
import pyhdust.lrr as lrr
import pyhdust.phc as phc

try:
    import matplotlib.pyplot as plt
except:
    print('# Warning! matplotlib module not installed!!!')

__author__ = "Leandro Rimulo"
__email__ = "lrrimulo@gmail.com"



def read_singlebe_outputs(lista_fullpath,path_to_singlebe_outputs):
    """
    Reading SINGLEBE files and creating 'todo_list' and 'outputs_list', 
    to be used by the other routines that work on the SINGLEBE output.
    
    INPUT: 
    . lista_fullpath: Name of the file containing the names of the 
    SINGLEBE output files + "instructions".
    . path_to_singlebe_outputs: 
    OUTPUT: 
    . todo_list: each element of the list contains the name of the 
    SINGLEBE file further "instructions" (which should be interpreted 
    by another program made by the user). 
    . outputs_list: each element of the list contains the outputs of 
    the specific SINGLEBE file from the list 'todo_list'.
    """
    
    func_name = sys._getframe().f_code.co_name
    todo_list=[]; outputs_list=[]
    
    ### Composing the todo_list based on the file lista_fullpath.
    f = open(lista_fullpath, "r"); nr_of_lines = sum(1 for line in f); f.close() 
    f = open(lista_fullpath, "r")
    for line in xrange(0,nr_of_lines):
        linha=f.readline(); linha=linha.split()
        ### The first word must be "FILE" and there must be at least one 
        ### word after it (supposed to be the name of a SINGLEBE output 
        ### file).
        if len(linha) > 1 and linha[0] == "FILE":
            todo_list.append([linha[i] for i in xrange(1,len(linha))])
    f.close()
    
    todelete=[] ### to be used in removing problematic files from todo_list.
    for i in xrange(0,len(todo_list)):
        try:
            f = open(path_to_singlebe_outputs+todo_list[i][0], "r")
            chave_elem=1
        except:
            print "<<",func_name,">>"
            print "I couldn't find "+path_to_singlebe_outputs+todo_list[i][0]
            todelete.append(i)
            chave_elem=0
            
        if chave_elem == 1: ### If the SINGLEBE output file was found, 
                            ### then continue composing the outputs_list.
            nr_of_lines = sum(1 for line in f); f.close(); ite=(nr_of_lines-6)/9;
            f = open(path_to_singlebe_outputs+todo_list[i][0], "r")
            ### Creating one list for the present SINGLEBE output file
            outputs_list.append([todo_list[i][0],\
                                    np.nan,\
                                    np.nan,np.nan,np.nan,np.nan,np.nan,\
                                    np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan])
            
            ### Reading the many constants from the beginning of the 
            ### SINGLEBE output file.
            lixo=f.readline(); lixo=lixo.split()
            lixo[0]=float(lixo[0])      ### [1][0] constant alpha parameter
            lixo[1]=float(lixo[1])      ### [1][1] stellar effective temperature [K]
            lixo[2]=float(lixo[2])      ### [1][2] tfac = Tdisk/Teff
            lixo[3]=float(lixo[3])      ### [1][3] disk temperature [K]
            lixo[4]=float(lixo[4])      ### [1][4] isothermal sound speed [cm/s]
            lixo[5]=float(lixo[5])      ### [1][5] stellar mass [Msun]
            lixo[6]=float(lixo[6])      ### [1][6] equatorial radius [Rsun]
            lixo[7]=float(lixo[7])      ### [1][7] vorb [cm/s]
            lixo[8]=float(lixo[8])      ### [1][8] orbital angular velocity at 
                                        ###     stellar equator [rad/s]
            lixo[9]=float(lixo[9])      ### [1][9] -dot_Minj0 [Msun/yr]
            lixo[10]=float(lixo[10])    ### [1][10] -dot_Minj0 [g/s]
            lixo[11]=float(lixo[11])    ### [1][11] cs/vorb
            lixo[12]=float(lixo[12])    ### [1][12] scale height at Rinj [cm]
            lixo[13]=float(lixo[13])    ### [1][13] rho_inj_0 [g/cm^3]
            lixo[14]=float(lixo[14])    ### [1][14] Sigma_inj_0 [g/cm^2]
            lixo[15]=float(lixo[15])    ### [1][15] inner boundary [stellar radii]
            lixo[16]=float(lixo[16])    ### [1][16] outer boundary [stellar radii]
            lixo[17]=float(lixo[17])    ### [1][17] injection radius [stellar radii]
            lixo[18]=int(lixo[18])      ### [1][18] number of grid points N
            lixo[19]=int(lixo[19])      ### [1][19] index of injection radius
            lixo[20]=float(lixo[20])    ### [1][20] dzeta (for the meshgrid)
            lixo[21]=float(lixo[21])    ### [1][21] dtau0
            lixo[22]=float(lixo[22])    ### [1][22] C factor for CFL condition
            lixo[23]=float(lixo[23])    ### [1][23] dtau
            lixo[24]=float(lixo[24])    ### [1][24] dt [s]
            lixo[25]=float(lixo[25])    ### [1][25] taumax
            lixo[26]=float(lixo[26])    ### [1][26] tauintval
            lixo[27]=int(lixo[27])      ### [1][27] intval=INT(tauintval/dtau)
            lixo[28]=int(lixo[28])      ### [1][28] ionoff (number of )
            lixo[29]=float(lixo[29])    ### [1][29] (injection of S at Rinj 
                                        ###     during the interval dt)_0
            outputs_list[-1][1]=[x for x in lixo]
            ### [2] tonoff's
            lixo=f.readline(); lixo=lixo.split()
            outputs_list[-1][2]=np.array([float(x) for x in lixo])
            ### [3] emdot's [g/s]
            lixo=f.readline(); lixo=lixo.split()
            outputs_list[-1][3]=np.array([float(x) for x in lixo])
            ### [4] zeta's (usually 400)
            lixo=f.readline(); lixo=lixo.split()
            outputs_list[-1][4]=np.array([float(x) for x in lixo])
            ### [5] r's (usually 400)
            lixo=f.readline(); lixo=lixo.split()
            outputs_list[-1][5]=np.array([float(x) for x in lixo])
            ### [6] d's (usually 400)
            lixo=f.readline(); lixo=lixo.split()
            outputs_list[-1][6]=np.array([float(x) for x in lixo])
            
            N=outputs_list[-1][1][18] ### (see above)
            ### [7] tau's
            array07=np.zeros(ite); array07[:]=np.nan
            ### [8] time's [s]
            array08=np.zeros(ite); array08[:]=np.nan
            ### [9] injections of S at Rinj during the interval dt
            array09=np.zeros(ite); array09[:]=np.nan
            ### [10] alpha's (usually 400)
            array10=np.zeros((N,ite)); array10[:,:]=np.nan
            ### [11] S (usually 400)
            array11=np.zeros((N,ite)); array11[:,:]=np.nan
            ### [12] Sigma's [g/cm^2] (usually 400)
            array12=np.zeros((N,ite)); array12[:,:]=np.nan
            ### [13] number of grid points for the vectors in [14] and [15]
            array13=np.zeros(ite); array13[:]=np.nan
            ### [14] amach's (usually 399 or less)
            array14=np.zeros((N-1,ite)); array14[:,:]=np.nan
            ### [15] decrate's [Msun/yr] (usually 399 or less)
            array15=np.zeros((N-1,ite)); array15[:,:]=np.nan
            
            for j in xrange(0,ite):
                lixo=f.readline(); array07[j]=float(lixo)
                lixo=f.readline(); array08[j]=float(lixo)
                lixo=f.readline(); array09[j]=float(lixo)
                lixo=f.readline(); lixo=lixo.split()
                for ii in xrange(0,len(lixo)):
                    array10[ii,j]=float(lixo[ii])
                lixo=f.readline(); lixo=lixo.split()
                for ii in xrange(0,len(lixo)):
                    array11[ii,j]=float(lixo[ii])
                lixo=f.readline(); lixo=lixo.split()
                for ii in xrange(0,len(lixo)):
                    array12[ii,j]=float(lixo[ii])
                lixo=f.readline(); array13[j]=float(lixo)
                lixo=f.readline(); lixo=lixo.split()
                for ii in xrange(0,len(lixo)):
                    array14[ii,j]=float(lixo[ii])
                lixo=f.readline(); lixo=lixo.split()
                for ii in xrange(0,len(lixo)):
                    array15[ii,j]=float(lixo[ii])
            array13.astype(int)
            outputs_list[-1][7]=array07
            outputs_list[-1][8]=array08
            outputs_list[-1][9]=array09
            outputs_list[-1][10]=array10
            outputs_list[-1][11]=array11
            outputs_list[-1][12]=array12
            outputs_list[-1][13]=array13
            outputs_list[-1][14]=array14
            outputs_list[-1][15]=array15
            
            
            f.close()
            print("The file "+outputs_list[-1][0]+" was read.")
    
    ### deleting problematic files from todo_list.
    todo_list=[todo_list[j] for j in range(0,len(todo_list)) \
                    if j not in todelete]
            
    return todo_list,outputs_list




### Now, some functions that return important constants associated with
### the SINGLEBE output file

def f_vorb(outputs_list_element):
    """
    Returns vorb in cm/s.
    """
    massa=outputs_list_element[1][5]
    relReq=outputs_list_element[1][6]
    return (phc.G.cgs*massa*phc.Msun.cgs/relReq/phc.Rsun.cgs)**0.5

def f_alphatau(outputs_list_element,units="cgs"):
    """
    Returns alphatau in s or in days.
    """
    vorb=f_vorb(outputs_list_element)
    massa=outputs_list_element[1][5]
    relReq=outputs_list_element[1][6]
    cs=outputs_list_element[1][11]*vorb
    if units=="cgs":
        return vorb**2./cs**2.*((relReq*phc.Rsun.cgs)**3./\
                (phc.G.cgs*massa*phc.Msun.cgs))**0.5
    if units=="days":
        return vorb**2./cs**2.*((relReq*phc.Rsun.cgs)**3./\
                (phc.G.cgs*massa*phc.Msun.cgs))**0.5/3600./24.

def f_asymp0(outputs_list_element):
    relRout=outputs_list_element[1][16]
    relRin=outputs_list_element[1][15]
    relRinj=outputs_list_element[1][17]
    Xi=(relRout**0.5-relRinj**0.5)/(relRout**0.5-relRin**0.5)
    return relRinj*relRinj*outputs_list_element[1][14]/Xi
        
def f_typ_decrate0(outputs_list_element,units="cgs"):
    Req=outputs_list_element[1][6]*phc.Rsun.cgs
    asymp0=f_asymp0(outputs_list_element)
    alphatau=f_alphatau(outputs_list_element)
    alpha0=outputs_list_element[1][0]
    if units=="cgs":
        return 2.*np.pi*Req*asymp0*alpha0*Req/alphatau
    if units=="astro":
        return 2.*np.pi*Req*asymp0*alpha0*Req/alphatau/phc.Msun.cgs*phc.yr.cgs

def f_std_decrate0(outputs_list_element,units="cgs"):
    relRout=outputs_list_element[1][16]
    relRin=outputs_list_element[1][15]
    typ_decrate0=f_typ_decrate0(outputs_list_element)
    if units=="cgs":
        return typ_decrate0/(relRout**0.5-relRin**0.5)
    if units=="astro":
        return typ_decrate0/(relRout**0.5-relRin**0.5)/phc.Msun.cgs*phc.yr.cgs
    
def f_std_mdJdt0(outputs_list_element):
    massa=outputs_list_element[1][5]
    relRout=outputs_list_element[1][16]
    lout=(phc.G.cgs*massa*phc.Msun.cgs*relRout*phc.Rsun.cgs)**0.5
    std_decrate0=f_std_decrate0(outputs_list_element)
    return std_decrate0*lout
    
def f_cs(outputs_list_element):
    return outputs_list_element[1][11]*f_vorb(outputs_list_element)

def f_typ_flow_speed(outputs_list_element):
    Req=outputs_list_element[1][6]*phc.Rsun.cgs
    alphatau=f_alphatau(outputs_list_element)
    alpha0=outputs_list_element[1][0]
    return alpha0*Req/alphatau

# 

def f_M_disk(outputs_list_element):
    M_disk=[]
    Req2=(outputs_list_element[1][6]*phc.Rsun.cgs)**2.
    ln_relR=outputs_list_element[4]
    d_ln_relR=np.array([outputs_list_element[4][j+1]-outputs_list_element[4][j] \
                    for j in xrange(0,len(outputs_list_element[4])-1)])
    relR=outputs_list_element[5]
    # over iterations...
    for i in xrange(0,len(outputs_list_element[8])):
        Sigmas=np.array([elem for elem in outputs_list_element[12][:,i]])
        M_disk.append(2.*np.pi*Req2*lrr.integrate_trapezia(relR**2.*Sigmas,d_ln_relR))
    return np.array([elem for elem in M_disk])

def f_J_disk(outputs_list_element):
    J_disk=[]
    Req=outputs_list_element[1][6]*phc.Rsun.cgs
    Req2=(outputs_list_element[1][6]*phc.Rsun.cgs)**2.
    massa=outputs_list_element[1][5]
    leq=(phc.G.cgs*massa*phc.Msun.cgs*Req)**0.5
    ln_relR=outputs_list_element[4]
    d_ln_relR=np.array([outputs_list_element[4][j+1]-outputs_list_element[4][j] \
                    for j in xrange(0,len(outputs_list_element[4])-1)])
    relR=outputs_list_element[5]
    # over iterations...
    for i in xrange(0,len(outputs_list_element[8])):
        Sigmas=np.array([elem for elem in outputs_list_element[12][:,i]])
        J_disk.append(2.*np.pi*Req2*leq*lrr.integrate_trapezia(relR**2.5*Sigmas,d_ln_relR))
    return np.array([elem for elem in J_disk])

def f_dMdt_disk(outputs_list_element):
    times=outputs_list_element[8]
    M_disk=f_M_disk(outputs_list_element)
    dMdt_disk=np.zeros(len(times))
    for i in xrange(0,len(times)):
        if i == 0:
            dMdt_disk[i]=(M_disk[i+1]-M_disk[i])/(times[i+1]-times[i])
        elif i == len(times)-1:
            dMdt_disk[i]=(M_disk[i]-M_disk[i-1])/(times[i]-times[i-1])
        else:
            dMdt_disk[i]=(M_disk[i+1]-M_disk[i-1])/(times[i+1]-times[i-1])
    return dMdt_disk
    
def f_dJdt_disk(outputs_list_element):
    times=outputs_list_element[8]
    J_disk=f_J_disk(outputs_list_element)
    dJdt_disk=np.zeros(len(times))
    for i in xrange(0,len(times)):
        if i == 0:
            dJdt_disk[i]=(J_disk[i+1]-J_disk[i])/(times[i+1]-times[i])
        elif i == len(times)-1:
            dJdt_disk[i]=(J_disk[i]-J_disk[i-1])/(times[i]-times[i-1])
        else:
            dJdt_disk[i]=(J_disk[i+1]-J_disk[i-1])/(times[i+1]-times[i-1])
    return dJdt_disk






def f_m_exponent(outputs_list_element):
    
    ln_relR=outputs_list_element[4]
    Sigma=outputs_list_element[12]
    N=len(Sigma[:,0]); ite=len(Sigma[0,:])
    m_exponent=np.zeros((N,ite)); m_exponent[:,:]=np.nan
    for iite in xrange(0,ite):
        for iN in xrange(0,N):
            if iN == 0:
                if Sigma[iN+1,iite] > 0. and Sigma[iN,iite] > 0.:
                    m_exponent[iN,iite]=-(np.log(Sigma[iN+1,iite])-\
                        np.log(Sigma[iN,iite]))/(ln_relR[iN+1]-ln_relR[iN])
            elif iN == N-1:
                if Sigma[iN,iite] > 0. and Sigma[iN-1,iite] > 0.:
                    m_exponent[iN,iite]=-(np.log(Sigma[iN,iite])-\
                        np.log(Sigma[iN-1,iite]))/(ln_relR[iN]-ln_relR[iN-1])
            else:
                if Sigma[iN+1,iite] > 0. and Sigma[iN-1,iite] > 0.:
                    m_exponent[iN,iite]=-(np.log(Sigma[iN+1,iite])-\
                        np.log(Sigma[iN-1,iite]))/(ln_relR[iN+1]-ln_relR[iN-1])                              
    return m_exponent
    

def f_massflux(outputs_list_element):
    ln_relR=outputs_list_element[4]
    relR=outputs_list_element[5]
    vorb=f_vorb(outputs_list_element)
    vorb2=vorb*vorb
    cs2=(outputs_list_element[1][11]*vorb)**2.
    massa=outputs_list_element[1][5]
    relReq=outputs_list_element[1][6]
    omegaorb=(phc.G.cgs*massa*phc.Msun.cgs/(relReq*phc.Rsun.cgs)**3.)**0.5
    C=-4.*np.pi*(relReq*phc.Rsun.cgs)**2.*omegaorb*cs2/vorb2
    Sigma=outputs_list_element[12]
    alpha=outputs_list_element[10]
    N=len(Sigma[:,0]); ite=len(Sigma[0,:])
    massflux=np.zeros((N,ite)); massflux[:,:]=np.nan
    for iite in xrange(0,ite):
        for iN in xrange(0,N):
            if iN == 0:
                massflux[iN,iite]=C*relR[iN]**-0.5*(alpha[iN+1,iite]*relR[iN+1]**2.*\
                    Sigma[iN+1,iite]-alpha[iN,iite]*relR[iN]**2.*Sigma[iN,iite])/\
                    (ln_relR[iN+1]-ln_relR[iN])
            elif iN == N-1:
                massflux[iN,iite]=C*relR[iN]**-0.5*(alpha[iN,iite]*relR[iN]**2.*\
                    Sigma[iN,iite]-alpha[iN-1,iite]*relR[iN-1]**2.*Sigma[iN-1,iite])/\
                    (ln_relR[iN]-ln_relR[iN-1])
            else:
                massflux[iN,iite]=C*relR[iN]**-0.5*(alpha[iN+1,iite]*relR[iN+1]**2.*\
                    Sigma[iN+1,iite]-alpha[iN-1,iite]*relR[iN-1]**2.*Sigma[iN-1,iite])/\
                    (ln_relR[iN+1]-ln_relR[iN-1])
    return massflux

def f_AMflux(outputs_list_element):
    ln_relR=outputs_list_element[4]
    relR=outputs_list_element[5]
    vorb=f_vorb(outputs_list_element)
    cs2=(outputs_list_element[1][11]*vorb)**2.
    relReq=outputs_list_element[1][6]
    C=-4.*np.pi*(relReq*phc.Rsun.cgs)**2.*cs2
    Sigma=outputs_list_element[12]
    alpha=outputs_list_element[10]
    N=len(Sigma[:,0]); ite=len(Sigma[0,:])
    AMflux=np.zeros((N,ite)); AMflux[:,:]=np.nan
    for iite in xrange(0,ite):
        for iN in xrange(0,N):
            if iN == 0:
                AMflux[iN,iite]=C*relR[iN]**0.5*(alpha[iN+1,iite]*relR[iN+1]**1.5*\
                    Sigma[iN+1,iite]-alpha[iN,iite]*relR[iN]**1.5*Sigma[iN,iite])/\
                    (ln_relR[iN+1]-ln_relR[iN])
            elif iN == N-1:
                AMflux[iN,iite]=C*relR[iN]**0.5*(alpha[iN,iite]*relR[iN]**1.5*\
                Sigma[iN,iite]-alpha[iN-1,iite]*relR[iN-1]**1.5*Sigma[iN-1,iite])/\
                (ln_relR[iN]-ln_relR[iN-1])
            else:
                AMflux[iN,iite]=C*relR[iN]**0.5*(alpha[iN+1,iite]*relR[iN+1]**1.5*\
                Sigma[iN+1,iite]-alpha[iN-1,iite]*relR[iN-1]**1.5*Sigma[iN-1,iite])/\
                (ln_relR[iN+1]-ln_relR[iN-1])
    return AMflux

def f_stelmassvarrate(outputs_list_element,ii="standard"):
    dMdt_disk=f_dMdt_disk(outputs_list_element)
    massflux=f_massflux(outputs_list_element)
    N=len(massflux[:,0]); ite=len(massflux[0,:])
    if ii == "standard":
        ii=N-2
    if ii == "last":
        ii=N-1
    stelmassvarrate=-(massflux[ii,:]+dMdt_disk[:])
    return stelmassvarrate

def f_stelAMvarrate(outputs_list_element,ii="standard"):
    dJdt_disk=f_dJdt_disk(outputs_list_element)
    AMflux=f_AMflux(outputs_list_element)
    N=len(AMflux[:,0]); ite=len(AMflux[0,:])
    if ii == "standard":
        ii=N-2
    if ii == "last":
        ii=N-1
    stelAMvarrate=-(AMflux[ii,:]+dJdt_disk[:])
    return stelAMvarrate

def f_stelmassvar(outputs_list_element,ii="standard"):
    jj=ii
    times=outputs_list_element[8]
    stelmassvarrate=f_stelmassvarrate(outputs_list_element,ii=jj)
    stelmassvar=np.zeros(len(times))
    for i in xrange(0,len(times)):
        if i == 0:
            stelmassvar[0]=0.
        else:
            stelmassvar[i]=stelmassvar[i-1]+0.5*(times[i]-times[i-1])*\
                (stelmassvarrate[i-1]+stelmassvarrate[i])
    return stelmassvar

def f_stelAMvar(outputs_list_element,ii="standard"):
    jj=ii
    times=outputs_list_element[8]
    stelAMvarrate=f_stelAMvarrate(outputs_list_element,ii=jj)
    stelAMvar=np.zeros(len(times))
    for i in xrange(0,len(times)):
        if i == 0:
            stelAMvar[0]=0.
        else:
            stelAMvar[i]=stelAMvar[i-1]+0.5*(times[i]-times[i-1])*\
                (stelAMvarrate[i-1]+stelAMvarrate[i])
    return stelAMvar

#



def f_dotMinjs(outputs_list_element):
    m_emdot0=-outputs_list_element[1][10]
    times=outputs_list_element[8]
    Req=outputs_list_element[1][6]*phc.Rsun.cgs
    massa=outputs_list_element[1][5]
    omegaorb=(phc.G.cgs*massa*phc.Msun.cgs/(Req)**3.)**0.5
    timesonoff=outputs_list_element[2]/omegaorb
    ratios=outputs_list_element[3]/outputs_list_element[1][9]
    dotMinjs=[]
    for itime in xrange(0,len(times)):
        for ionoffs in xrange(0,len(timesonoff)):
            if ionoffs == 0 and times[itime] >= 0. and \
                    times[itime] < timesonoff[ionoffs]:
                dotMinjs.append(m_emdot0*ratios[ionoffs])
            if ionoffs > 0 and ionoffs < len(timesonoff)-1 and \
                    times[itime] >= timesonoff[ionoffs-1] and \
                    times[itime] < timesonoff[ionoffs]:
                dotMinjs.append(m_emdot0*ratios[ionoffs])
            if ionoffs == len(timesonoff)-1 and \
                    times[itime] >= timesonoff[ionoffs-1]:
                dotMinjs.append(m_emdot0*ratios[ionoffs])
    return np.array(dotMinjs)   



def f_typdecrates(outputs_list_element):
    if 1==2:
        massa=outputs_list_element[1][5]
        relReq=outputs_list_element[1][6]
        omegaorb=(phc.G.cgs*massa*phc.Msun.cgs/(relReq*phc.Rsun.cgs)**3.)**0.5
        alphatau=f_alphatau(outputs_list_element)
        asymp0=f_asymp0(outputs_list_element)
        typdecrate_auxi=2.*np.pi*(relReq*phc.Rsun.cgs)**2.*asymp0/alphatau
        timesonoff=outputs_list_element[2]/omegaorb
        ratios=outputs_list_element[3]/outputs_list_element[1][9]
        #
        times=outputs_list_element[8]
        kinj=outputs_list_element[1][19]
        alphas_inj=outputs_list_element[10][kinj-1,:]
        #
        typdecrates=[]
        for itime in xrange(0,len(times)):
            for ionoffs in xrange(0,len(timesonoff)):
                if ionoffs == 0 and times[itime] >= 0. and \
                        times[itime] < timesonoff[ionoffs]:
                    typdecrates.append(typdecrate_auxi*ratios[ionoffs]*alphas_inj[itime])
                if ionoffs > 0 and ionoffs < len(timesonoff)-1 and \
                        times[itime] >= timesonoff[ionoffs-1] and \
                        times[itime] < timesonoff[ionoffs]:
                    typdecrates.append(typdecrate_auxi*ratios[ionoffs]*alphas_inj[itime])
                if ionoffs == len(timesonoff)-1 and \
                        times[itime] >= timesonoff[ionoffs-1]:
                    typdecrates.append(typdecrate_auxi*ratios[ionoffs]*alphas_inj[itime])
        return np.array(typdecrates)
    if 1==1:
        m_emdot0=-outputs_list_element[1][10]
        Req=outputs_list_element[1][6]*phc.Rsun.cgs
        massa=outputs_list_element[1][5]
        omegaorb=(phc.G.cgs*massa*phc.Msun.cgs/(Req)**3.)**0.5
        relRinj=outputs_list_element[1][17]
        relRin=outputs_list_element[1][15]
        #
        times=outputs_list_element[8]
        timesonoff=outputs_list_element[2]/omegaorb
        ratios=outputs_list_element[3]/outputs_list_element[1][9]
        C=m_emdot0*(relRinj**0.5-relRin**0.5)
        #
        typdecrates=[]
        for itime in xrange(0,len(times)):
            for ionoffs in xrange(0,len(timesonoff)):
                if ionoffs == 0 and times[itime] >= 0. and \
                        times[itime] < timesonoff[ionoffs]:
                    typdecrates.append(C*ratios[ionoffs])
                if ionoffs > 0 and ionoffs < len(timesonoff)-1 and \
                        times[itime] >= timesonoff[ionoffs-1] and \
                        times[itime] < timesonoff[ionoffs]:
                    typdecrates.append(C*ratios[ionoffs])
                if ionoffs == len(timesonoff)-1 and \
                        times[itime] >= timesonoff[ionoffs-1]:
                    typdecrates.append(C*ratios[ionoffs])
        return np.array(typdecrates)                
        
        
        

def f_asymps(outputs_list_element):
    typ_decrate0=f_typ_decrate0(outputs_list_element)
    massa=outputs_list_element[1][5]
    relReq=outputs_list_element[1][6]
    omegaorb=(phc.G.cgs*massa*phc.Msun.cgs/(relReq*phc.Rsun.cgs)**3.)**0.5
    alphatau=f_alphatau(outputs_list_element)
    ratios=outputs_list_element[3]/outputs_list_element[1][9]
    times=outputs_list_element[8]
    timesonoff=outputs_list_element[2]/omegaorb
    kinj=outputs_list_element[1][19]
    alphas_inj=outputs_list_element[10][kinj-1,:]
    C=2.*np.pi*(relReq*phc.Rsun.cgs)**2./alphatau
    asymps=[]
    for itime in xrange(0,len(times)):
        for ionoffs in xrange(0,len(timesonoff)):
            if ionoffs == 0 and times[itime] >= 0. and \
                    times[itime] < timesonoff[ionoffs]:
                asymps.append(typ_decrate0*ratios[ionoffs]/(C*alphas_inj[itime]))
            if ionoffs > 0 and ionoffs < len(timesonoff)-1 and \
                    times[itime] >= timesonoff[ionoffs-1] and times[itime] < timesonoff[ionoffs]:
                asymps.append(typ_decrate0*ratios[ionoffs]/(C*alphas_inj[itime]))
            if ionoffs == len(timesonoff)-1 and \
                    times[itime] >= timesonoff[ionoffs-1]:
                asymps.append(typ_decrate0*ratios[ionoffs]/(C*alphas_inj[itime]))
    return asymps

def f_std_decrates(outputs_list_element):
    relRin=outputs_list_element[1][15]
    relRout=outputs_list_element[1][16]
    typdecrates=f_typdecrates(outputs_list_element)
    C=(relRout**0.5-relRin**0.5)
    return typdecrates/C

def f_std_mdJdt(outputs_list_element):
    massa=outputs_list_element[1][5]
    Req=outputs_list_element[1][6]*phc.Rsun.cgs
    relRout=outputs_list_element[1][16]
    std_decrates=f_std_decrates(outputs_list_element)
    C=(phc.G.cgs*massa*phc.Msun.cgs*Req*relRout)**0.5
    return C*std_decrates

def f_int_std_decrates(outputs_list_element):
    times=outputs_list_element[8]
    std_decrates=f_std_decrates(outputs_list_element)
    int_std_decrates=np.zeros(len(times))
    for i in xrange(0,len(times)):
        if i == 0:
            int_std_decrates[0]=0.
        else:
            int_std_decrates[i]=int_std_decrates[i-1]+\
                0.5*(times[i]-times[i-1])*(std_decrates[i-1]+std_decrates[i]) 
    return int_std_decrates

def f_int_std_mdJdt(outputs_list_element):
    times=outputs_list_element[8]
    std_mdJdt=f_std_mdJdt(outputs_list_element)
    int_std_mdJdt=np.zeros(len(times))
    for i in xrange(0,len(times)):
        if i == 0:
            int_std_mdJdt[0]=0.
        else:
            int_std_mdJdt[i]=int_std_mdJdt[i-1]+\
                0.5*(times[i]-times[i-1])*(std_mdJdt[i-1]+std_mdJdt[i]) 
    return int_std_mdJdt

def f_timepars(outputs_list_element):
    """
    According to relation dtimepar=dt/tau(t), 
    obtains the vector 'timepars', from 0 to maximum time parameter 
    (from reading the vector 'times', from 0 to maximum calculated time [s]).
    """
    
    times=outputs_list_element[8] ### (N)-shaped array with times in s
    kinj=outputs_list_element[1][19]    ### index of the injection radius 
                                        ### as defined in SINGLEBE
    alphas_inj=outputs_list_element[10][kinj-1,:]   ### (ite)-shaped array
                                                    ### with the values of 
                                                    ### alpha_inj
    alphatau=f_alphatau(outputs_list_element)
    
    timepars=np.zeros(len(times))
    for i in range(0,len(times)):
        if i == 0:
            timepars[0]=0.
        else:
            timepars[i]=timepars[i-1]+0.5*(times[i]-times[i-1])*\
                (alphas_inj[i-1]+alphas_inj[i])/alphatau
    return timepars






def f_convertion_to_mohammad_paper(outputs_list_element):
    relRinj=outputs_list_element[1][17]
    relRout=outputs_list_element[1][16]
    Lambda=1./(1.-relRout**-0.5)
    Req=outputs_list_element[1][6]*phc.Rsun.cgs
    massa=outputs_list_element[1][5]
    C=Lambda*(phc.G.cgs*massa*phc.Msun.cgs*Req)**0.5*(relRinj**0.5-1.)*\
        phc.Msun.cgs/phc.yr.cgs
    return C













def filesinfo(outputs_list,filesinfo_output):
    ### needs to be improved!
    """
    
    """


    extfile=open(filesinfo_output,'w')


    for ifile in range(0,len(outputs_list)):
    
        extfile.write("FILE #: "+str(ifile)+"\n")
    
        extfile.write("FILE NAME: "+outputs_list[ifile][0]+"\n")
        extfile.write("   STELLAR PARAMETERS"+"\n")
        extfile.write("mass [Msun]                       = "+\
            str(outputs_list[ifile][1][5])+"\n")
        extfile.write("equatorial radius [Rsun]          = "+\
            str(outputs_list[ifile][1][6])+"\n")
        extfile.write("stellar effective temperature [K] = "+\
            str(outputs_list[ifile][1][1])+"\n")

        extfile.write("   DISK PARAMETERS"+"\n")
        extfile.write("disk temperature [K]                 = "+\
            str(outputs_list[ifile][1][3])+"\n")
        extfile.write("mean molecular weight                = ?"+"\n")
        extfile.write("asymptotic surface density [g/cm^2]  = "+\
            str(f_asymp0(outputs_list[ifile]))+"\n")
        extfile.write("alpha0                               = "+\
            str(outputs_list[ifile][1][0])+"\n")
        extfile.write("alphatau [days]                      = "+\
            str(f_alphatau(outputs_list[ifile],units="days"))+"\n")
        extfile.write("isothermal sound speed [km/s]        = "+\
            str(f_cs(outputs_list[ifile])/100000.)+"\n")
        extfile.write("typical flow speed [km/s]            = "+\
            str(f_typ_flow_speed(outputs_list[ifile])/100000.)+"\n")
        extfile.write("typical decretion rate [Msun/yr]     = "+\
            str(f_typ_decrate0(outputs_list[ifile],units="astro"))+"\n")
        extfile.write("steady-state decrate [Msun/yr]       = "+\
            str(f_std_decrate0(outputs_list[ifile],units="astro"))+"\n")
        extfile.write("steady-state AMloss rate [cgs]       = "+\
            str(f_std_mdJdt0(outputs_list[ifile]))+"\n")

        extfile.write("   SOURCE OF MASS AND RADIAL GRID"+"\n")
        extfile.write("rin     = "+str(outputs_list[ifile][1][15])+\
            " "+str(1)+"\n")
        extfile.write("rinject = "+str(outputs_list[ifile][1][17])+\
            " "+str(outputs_list[ifile][1][19])+"\n")
        extfile.write("rout    = "+str(outputs_list[ifile][1][16])+\
            " "+str(outputs_list[ifile][1][18])+"\n")
#
#       extfile.write("   REFERENCE VALUES"+"\n")
#       extfile.write("rho00 (g/cm^3)       = "+str(rho00[ifile])+"\n")
#       extfile.write("ndens0 (1/cm^3)      = "+str(ndens0[ifile])+"\n")
#       extfile.write("sigma0 (g/cm^2)      = "+str(sigma0[ifile])+"\n")
#       extfile.write("emdot0 (g/s)         = "+str(emdot0[ifile]/year*solmas)+"\n")
#       extfile.write("emdot0 (solmas/year) = "+str(emdot0[ifile])+"\n")
#
#       extfile.write("   TIMESTEP"+"\n")
#       extfile.write("timestep (days) =                   "+str(dt[ifile]/day)+"\n")
#       deltatimeprt=(timesecond-timefirst)/day
#       extfile.write("output rough time interval (days) = "+str(deltatimeprt)+"\n")
#
## if len(emdot[:,ifile]) < tresmil2:
        extfile.write("   MASS INJECTION HISTORY"+"\n")
        extfile.write("mass injection rate [units of reference value]  ;  until time [yrs]"+\
            "\n")
        for iiii in range(0,len(outputs_list[ifile][2])):
#      extfile.write(str(emdot[iiii,ifile]*year/solmas)+' '+str(tonoff[iiii,ifile]/(year*omega0[ifile]))+'\n')
            extfile.write(str(outputs_list[ifile][3][iiii]/\
                outputs_list[ifile][1][9])+" "+str(outputs_list[ifile][2][iiii]/\
                (phc.yr.cgs*outputs_list[ifile][1][8]))+"\n")

        extfile.write("   "+"\n")


    extfile.close()

    return







