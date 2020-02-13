#-*- coding:utf-8 -*-

import sys
import numpy as np
import pyhdust.lrr as lrr
import pyhdust.phc as phc

__author__ = "Leandro Rimulo"
__email__ = "lrrimulo@gmail.com"

def rocheparams(param,tipo):
    """
    Reference: 2013A&ARv..21...69R
    
    Input 'param' as the real number corresponding to 'tipo' (where 'tipo' 
    can be 'W', 'omega', 'oblateness' or 'gamma'). 
    
    Output: ob,omega,gamma,W
    """
    
    func_name = sys._getframe().f_code.co_name
    
    if tipo=='W':
        if param < 0.0 or param > 1.0: 
            print('<<',func_name,'>>')
            print('Please, give me <<W>> between 0 and 1!')
            print('')
            ob=np.nan; omega=np.nan; gamma=np.nan; W=np.nan
        else:
            W=param
            ob=1.+0.5*W*W
            if W != 0.:
                omega=np.sqrt(27./8.*ob**-3.*W**2.)
            else:
                omega=0.
            gamma=2./3.*ob*omega
    elif tipo=='omega':
        if param < 0.0 or param > 1.0: 
            print('<<',func_name,'>>')
            print('Please, give me <<omega>> between 0 and 1!')
            print('')
            ob=np.nan; omega=np.nan; gamma=np.nan; W=np.nan
        else:
            omega=param
            if omega != 0.:
                ob=3./omega*np.cos(np.pi/3.+1./3.*np.arccos(omega))
            else:
                ob=1.0
            gamma=2./3.*ob*omega
            if ob != 1.0:
                W=np.sqrt(2.*(ob-1.0))
            else:
                W=0.
    elif tipo=='oblateness':
        if param < 1.0 or param > 1.5: 
            print('<<',func_name,'>>')
            print('Please, give me <<oblateness>> between 1 and 1.5!')
            print('')
            ob=np.nan; omega=np.nan; gamma=np.nan; W=np.nan
        else:
            ob=param
            if ob != 1.0:
                omega=np.sqrt(27./4.*ob**-3.0*(ob-1.0))
            else:
                omega=0.
            if omega != 0.:
                gamma=2.0*np.cos(np.pi/3.+1./3.*np.arccos(omega))
            else:
                gamma=0.
            if ob != 1.0:
                W=np.sqrt(2.*(ob-1.0))
            else:
                W=0.
    elif tipo=='gamma':
        if param < 0.0 or param > 1.0: 
            print('<<',func_name,'>>')
            print('Please, give me <<gamma>> between 0 and 1!')
            print('')
            ob=np.nan; omega=np.nan; gamma=np.nan; W=np.nan
        else:
            gamma=param
            if gamma != 0.0:
                omega=np.cos(3.*np.arccos(gamma/2.)-np.pi)
            else:
                omega=0.
            if omega != 0.:
                ob=1.5*gamma/omega
            else:
                ob=1.0
            if ob != 1.0:
                W=np.sqrt(2.*(ob-1.0))
            else:
                W=0.
    else:
        print('<<',func_name,'>>')
        print(tipo+' ???')
        print('I do not know the meaning of this parameter you gave me!')
        print('')
        ob=np.nan; omega=np.nan; gamma=np.nan; W=np.nan
    
    return ob,omega,gamma,W


###########################################

def f_ctes1(r0,M,par,tp):
    """
    Calculates the values of the constants associated with the Roche 
    model for a single star rotating as a solid body, given r0, M and 
    one of the three (F0, T0 or L0).

    | INPUT:
    |. r0 [Rsun]
    |. M [Msun]
    |. par (should receive any of three parameters: F0 [erg/s cm2], 
    | T0 [K] and L0 [Lsun]) 
    |. tp (should receive any of the corresponding lables: 
    | 'flux', 'temp' and 'lum')
    NOTE: If not all of the input parameters are available, the NaN 
    should be provided instead. This will result in some NaNs being 
    returned, because some of the parameters hadn't the necessary input 
    parameters to be calculated.
    
    | OUTPUT:
    |. A0 [cm2]
    |. V0 [cm3]
    |. Psi0 [erg/g]
    |. Omega0 [1/s]
    |. g0 [cm/s2]
    |. F0 [erg/s cm2]
    |. T0 [K]
    |. L0 [Lsun]

    """

    func_name = sys._getframe().f_code.co_name
    
    if np.isnan(r0) == False:
        A0=4.*np.pi*(r0*phc.Rsun.cgs)**2.
        V0=4.*np.pi*(r0*phc.Rsun.cgs)**3./3.
        if np.isnan(M) == False:
            Psi0=-phc.G.cgs*M*phc.Msun.cgs/(r0*phc.Rsun.cgs)
            Omega0=np.sqrt((8./27.)*phc.G.cgs*M*phc.Msun.cgs/(r0*phc.Rsun.cgs)**3.)
            g0=phc.G.cgs*M*phc.Msun.cgs/(r0*phc.Rsun.cgs)**2.
        else:
            Psi0=np.nan; Omega0=np.nan; g0=np.nan 

    else:
        V0=np.nan; A0=np.nan; Psi0=np.nan; Omega0=np.nan; g0=np.nan 	

    if np.isnan(par) == False:
        if tp == "flux":
            F0=par
            T0=(F0/phc.sigma.cgs)**0.25
            if np.isnan(A0) == False:
                L0=F0*A0/phc.Lsun.cgs
            else: L0=np.nan
        elif tp == "temp":
            T0=par
            F0=phc.sigma.cgs*T0**4.
            if np.isnan(A0) == False:
                L0=F0*A0/phc.Lsun.cgs
            else: L0=np.nan
        elif tp == "lum":
            L0=par
            if np.isnan(A0) == False:
                F0=L0*phc.Lsun.cgs/A0
                T0=(F0/phc.sigma.cgs)**0.25
            else: 
                F0=np.nan; T0=np.nan
        else:
            print('<<',func_name,'>>')
            print('WARNING: Type of parameter not recognized: '+tp)
            F0=np.nan; T0=np.nan; L0=np.nan
    else:
        F0=np.nan; T0=np.nan; L0=np.nan
        
    return V0,A0,Psi0,Omega0,g0,F0,T0,L0




def cte_veq(r0,M,omega,psi):
    """
    Returns the velocity in the equator in km/s.
    """
    V0,A0,Psi0,Omega0,g0,F0,T0,L0=f_ctes1(r0,M,np.nan,"lum")
    ob,omega,gamma,W=rocheparams(omega,"omega")

    return Omega0*(r0*phc.Rsun.cgs)*psi**0.5*omega*ob*1e-5





###########################################

def psi_Fremat2005(W=np.nan,mass=np.nan):
    """Reference: 2005A&A...440..305F"""
    func_name = sys._getframe().f_code.co_name
    psi=1.0
    
    if not np.isnan(mass) and not np.isnan(W):
        tau=(0.0072+0.008*W)*W
        pm=5.66+9.43/mass**2.0
        psi=1.0/(1.0-pm*tau)
    else:
        print('<<',func_name,'>>')
        print('I did not receive meaningful values of the mass of the star or W.')
        print('Therefore, I considered the polar radius unvariable with W.')
        print('')
        psi=1.0
            
    return psi

def lum_fac_Fremat2005(W=np.nan,mass=np.nan):
    """Reference: 2005A&A...440..305F"""
    func_name = sys._getframe().f_code.co_name
    lum_fac=1.0
    
    if not np.isnan(mass) and not np.isnan(W):
        a=0.675+0.046*mass**0.5
        b=52.71+20.63*mass**0.5
        tau=(0.0072+0.008*W)*W
        lum_fac=a+(1.0-a)*np.exp(-b*tau)
    else:
        print('<<',func_name,'>>')
        print('I did not receive meaningful values of the mass of the star or W.')
        print('Therefore, I considered the luminosity to be unvariable with W.')
        print('')
        lum_fac=1.0
        
    return lum_fac

###########################################

def go_first_quadrant(theta):
    theta=abs(theta%(2.*np.pi))
    if theta > 0.5*np.pi and theta < 1.0*np.pi: theta=np.pi-theta
    if theta >= 1.0*np.pi and theta < 1.5*np.pi: theta=-np.pi+theta
    if theta >= 1.5*np.pi and theta < 2.0*np.pi: theta=2.*np.pi-theta
    return theta

def f_rocheradius(theta,omega,psi):
    """
    Returns the function s_psi(theta) for the Roche surface of a 
    single star rotating with parameters omega and psi.
    
    Input: 
    . 'theta' is the angle theta in radians. 
    . 'omega'
    . 'psi'
    """
    if omega != 0.0:
        theta=go_first_quadrant(theta)
        if theta!=0.0*np.pi:
            rocheradius=3./(psi*omega*np.sin(theta))*\
                np.cos(np.pi/3.0+np.arccos(omega*np.sin(theta))/3.0)
        else:
            rocheradius=1.0/psi
    else:
        rocheradius=1.0/psi
    return rocheradius


def f_tilde_r_roche(omega,psi,r):


    ob,omega,gamma,W=rocheparams(omega,"omega")
    return r*psi/ob


    
def f_hdustquadratic(theta,omega,psi):
    """
    Returns the function s_psi(theta) for the surface of the 
    ellipsoid used by HDUST.

    Input: 
    . 'theta' is the angle theta in radians. 
    . 'omega'
    . 'psi'    
    """
    rochepole=f_rocheradius(0.0,omega,psi)
    rocheeq=f_rocheradius(0.5*np.pi,omega,psi)
    r_hdustquad=rocheeq*rochepole/np.sqrt(rocheeq**2.*np.cos(theta)**2.+rochepole**2.*np.sin(theta)**2.)
    return r_hdustquad















def f_eff_gravity_vector(theta,omega,psi,r):
    """
    

    Input: 
    . 'theta' is the angle theta in radians. 
    . 'omega'
    . 'psi'  
    . 'r' is the scaled position.
    """
    eff_gravity=np.sqrt((r**-2.0-(8./27.)*omega**2.0*psi**3.0*\
                            r*np.sin(theta)**2.0)**2.0\
                            +((8./27.)*omega**2.0*psi**3.0*\
                            r*np.sin(theta)*np.cos(theta))**2.0)
                            
    r_component=-r**-2.0+(8./27.)*omega**2.0*psi**3.0\
        *r*np.sin(theta)**2.0
    theta_component=(8./27.)*omega**2.0*psi**3.0*\
        r*np.sin(theta)*np.cos(theta)
    return np.array([r_component,theta_component])

def f_eff_gravity(theta,omega,psi,r):
    """
    Returns the absolute value of the scaled effective gravity for 
    every point (r,theta), for a Roche single star with parameters 
    omega and psi.

    Input: 
    . 'theta' is the angle theta in radians. 
    . 'omega'
    . 'psi'  
    . 'r' is the scaled position.
    """
    eff_gravity_vector=f_eff_gravity_vector(theta,omega,psi,r)
    r_component=eff_gravity_vector[0]
    theta_component=eff_gravity_vector[1]
    eff_gravity=np.sqrt(r_component**2.0+theta_component**2.0)
    return eff_gravity

def f_cosepsilon(theta,omega,psi,r):
    """
	
    """
    eff_gravity_vector=f_eff_gravity_vector(theta,omega,psi,r)
    r_component=eff_gravity_vector[0]
    theta_component=eff_gravity_vector[1]

    return 1./np.sqrt(1.+(theta_component/r_component)**2.)



def f_normal_roche(theta,omega,psi,r,unity=True):

    r_component=1.
    eff_gravity_vector=f_eff_gravity_vector(theta,omega,psi,r)
    theta_component=eff_gravity_vector[1]/eff_gravity_vector[0]
    if unity == True:
        cosepsilon=f_cosepsilon(theta,omega,psi,r)
        r_component=r_component*cosepsilon
        theta_component=theta_component*cosepsilon

    return np.array([r_component,theta_component])














def f_drochedtheta(theta,omega,psi):
    """
    Returns dr/dtheta(theta), for a point in the single star Roche 
    surface.

    Input: 
    . 'theta' is the angle theta in radians. 
    . 'omega'
    . 'psi'  
    """
    theta=go_first_quadrant(theta)
    if omega == 1.0 and theta == 0.5*np.pi:
        drochedtheta=0.0 # Removable discontinuity corrected.
    else:
        rocheradius=f_rocheradius(theta,omega,psi)
        drochedtheta=(8./27.)*omega**2.0*psi**3.0*\
                rocheradius**2.*np.sin(theta)*np.cos(theta)*\
                (rocheradius**(-2.0)-(8./27.)*omega**2.0*psi**3.0*\
                rocheradius*np.sin(theta)**2.0)**(-1.)
    return drochedtheta

def f_deffgravdtheta_roche(theta,omega,psi):
    """
    Returns d(effgrav)/dtheta(theta), for a point in the single star 
    Roche surface.

    Input: 
    . 'theta' is the angle theta in radians. 
    . 'omega'
    . 'psi'  
    """
    rocheradius=f_rocheradius(theta,omega,psi)
    drochedtheta=f_drochedtheta(theta,omega,psi)
    g1=-rocheradius**(-2.0)+(8./27.)*omega**2.0*psi**3.0*\
        rocheradius*np.sin(theta)**2.0
    g2=(8./27.)*omega**2.0*psi**3.0*rocheradius*np.sin(theta)*\
        np.cos(theta)
    dg1=2.0*rocheradius**-3.*drochedtheta+8./27.*psi**3.*\
        omega**2.*(np.sin(theta)**2.*drochedtheta+\
        2.*rocheradius*np.sin(theta)*np.cos(theta))
    dg2=8./27.*psi**3.*omega**2.*(np.sin(theta)*np.cos(theta)*\
        drochedtheta+rocheradius*np.cos(2.*theta))
    deffgravdtheta_roche=(g1*dg1+g2*dg2)/np.sqrt(g1*g1+g2*g2)
    return deffgravdtheta_roche
    




    

def f_roche_dafac(rocheradius,drochedtheta):
    """Returns the factor multiplied to the element of area of the 
    equipotential surface"""
    roche_dafac=np.sqrt(rocheradius**4.+rocheradius**2.*\
        drochedtheta**2.)
    return roche_dafac
    
def f_ds(theta,dtheta,dphi,dafac):
    ds=dafac*np.sin(theta)*dtheta*dphi/4./np.pi
    return ds
    
def f_integral_surfaceroche(theta,omega,psi,x):
    """
    Returns the integral of x over the equipotential surface.
    $I = \int_S x(\theta) \mathrm{d}S$
    
    
    """
    if len(theta) < 2 or len(x) < 2 or len(theta) != len(x):
        func_name = sys._getframe().f_code.co_name
        print('<<',func_name,'>>')
        print('There is something wrong with the integration!')
        print('')
        return np.nan
    for i in xrange(0,len(theta)):
        if not (theta[i] >= 0.0*np.pi and theta[i] <= 0.5*np.pi):
            func_name = sys._getframe().f_code.co_name
            print('<<',func_name,'>>')
            print('The thetas you gave me are out of the desired range!')
            print('')
            return np.nan
        if i > 0 and theta[i] <= theta[i-1]:
            func_name = sys._getframe().f_code.co_name
            print('<<',func_name,'>>')
            print('The thetas you gave me are not in ascending order!')
            print('')
            return np.nan
    npts=len(theta)
    dtheta=np.array([theta[i+1]-theta[i] \
        for i in xrange(0,len(theta)-1)])
    rocheradius=np.array([f_rocheradius(theta[i],omega,psi) \
        for i in xrange(0,len(dtheta))])
    drochedtheta=np.array([f_drochedtheta(theta[i],omega,psi) \
        for i in xrange(0,len(dtheta))])
    roche_dafac=np.array([f_roche_dafac(rocheradius[i],\
        drochedtheta[i]) for i in xrange(0,len(dtheta))])
    ds=np.array([f_ds(theta[i],dtheta[i],2.*np.pi,roche_dafac[i]) \
        for i in xrange(0,len(dtheta))])

    integral_roche=2.*lrr.integrate_trapezia(x,ds)
    return integral_roche








def f_dalinha(theta,omega,psi,dtheta):
    """
    
    """
    rocheradius=f_rocheradius(theta,omega,psi)
    cosepsilon=f_cosepsilon(theta,omega,psi,rocheradius)
    return 0.5/cosepsilon*rocheradius**2.*np.sin(theta)*dtheta


def f_integral_surfaceroche_v2(theta,omega,psi,x):
    """
    Returns the integral of x over the equipotential surface.
    $I = \int_S x(\theta) \mathrm{d}S$
    
    
    """
    if len(theta) < 2 or len(x) < 2 or len(theta) != len(x):
        func_name = sys._getframe().f_code.co_name
        print('<<',func_name,'>>')
        print('There is something wrong with the integration!')
        print('')
        return np.nan
    for i in xrange(0,len(theta)):
        if not (theta[i] >= 0.0*np.pi and theta[i] <= 0.5*np.pi):
            func_name = sys._getframe().f_code.co_name
            print('<<',func_name,'>>')
            print('The thetas you gave me are out of the desired range!')
            print('')
            return np.nan
        if i > 0 and theta[i] <= theta[i-1]:
            func_name = sys._getframe().f_code.co_name
            print('<<',func_name,'>>')
            print('The thetas you gave me are not in ascending order!')
            print('')
            return np.nan
    npts=len(theta)
    dtheta=np.array([theta[i+1]-theta[i] \
        for i in xrange(0,len(theta)-1)])
    rocheradius=np.array([f_rocheradius(theta[i],omega,psi) \
        for i in xrange(0,len(dtheta))])
    ds=np.array([f_dalinha(theta[i],omega,psi,dtheta[i]) \
        for i in xrange(0,len(dtheta))])

    integral_roche=2.*lrr.integrate_trapezia(x,ds)
    return integral_roche







####
def f_flux_roche_betalaw(theta,omega,psi,lum_fac,r,betacte,npts=18001,\
        conservation=True):
    """
    
    """
    effgrav=f_eff_gravity(theta,omega,psi,r)
    #
    if conservation == True:
        thetamin=0.0; thetamax=0.5*np.pi
        thetas=np.array([thetamin+(thetamax-thetamin)*\
            float(i)/float(npts-1) for i in xrange(0,npts)])
        rocheradiuss=np.array([f_rocheradius(thetas[i],omega,psi) \
            for i in xrange(0,npts)])
        effgravs=np.array([f_eff_gravity(thetas[i],omega,psi,\
            rocheradiuss[i]) for i in xrange(0,npts)])
        cte=lum_fac/f_integral_surfaceroche_v2(thetas,omega,psi,\
            effgravs**(4.*betacte))
        flux=cte*effgrav**(4.*betacte)
    else: flux=effgrav**(4.*betacte)
    return flux

    
def f_vw(omega,psi,npts):
    """Returns the volume inside the equipotential surface."""
    thetamin=0.0; thetamax=0.5*np.pi
    theta=np.array([thetamin+(thetamax-thetamin)*\
        float(i)/float(npts-1) for i in xrange(0,npts)])
    rocheradius=np.array([f_rocheradius(theta[i],omega,psi) \
        for i in xrange(0,len(theta))])
    integ=np.array([0.5*np.sin(theta[i])*rocheradius[i]**3. \
        for i in xrange(0,len(theta))])
    dtheta=np.array([theta[i+1]-theta[i] \
        for i in xrange(0,len(theta)-1)])
    vw=2.*lrr.integrate_trapezia(integ,dtheta)
    return vw
    
def f_meanflux(omega,psi,lum_fac,npts=18001):
    thetamin=0.0; thetamax=0.5*np.pi
    theta=np.array([thetamin+(thetamax-thetamin)*\
        float(i)/float(npts-1) for i in xrange(0,npts)])
    x=np.array([1. for i in xrange(0,npts)])
    meanflux=lum_fac/f_integral_surfaceroche(theta,omega,psi,x)
    return meanflux

def Tmean4_qu(omega,psi,lum_fac,npts=18001):
    Tmean4_qu=f_meanflux(omega,psi,lum_fac,npts)**0.25
    return Tmean4_qu


def f_meanflux_v2(omega,psi,lum_fac,npts=18001):
    thetamin=0.0; thetamax=0.5*np.pi
    theta=np.array([thetamin+(thetamax-thetamin)*\
        float(i)/float(npts-1) for i in xrange(0,npts)])
    x=np.array([1. for i in xrange(0,npts)])
    meanflux=lum_fac/f_integral_surfaceroche_v2(theta,omega,psi,x)
    return meanflux

def Tmean4_qu_v2(omega,psi,lum_fac,npts=18001):
    Tmean4_qu=f_meanflux_v2(omega,psi,lum_fac,npts)**0.25
    return Tmean4_qu



    

###########################################
# Symmetric quartics

def f_symmquartic(r,theta,A,B,C,D,E,F):
	x=r*np.sin(theta)
	y=r*np.cos(theta)
	return A+B*x**2.+C*y**2.+D*x**2.*y**2.+E*x**4.+F*y**4.
	
def f_symmquartic_x(r,theta,A,B,C,D,E,F):
	x=r*np.sin(theta)
	y=r*np.cos(theta)
	return 2.*B*x+2.*D*x*y**2.+4.*E*x**3.	

def f_symmquartic_y(r,theta,A,B,C,D,E,F):
	x=r*np.sin(theta)
	y=r*np.cos(theta)
	return 2.*C*y+2.*D*x**2.*y+4.*F*y**3.


def f_symmquartic_prod(r,theta,A1,B1,C1,D1,E1,F1,A2,B2,C2,D2,E2,F2):
	Q=f_symmquartic(r,theta,A1,B1,C1,D1,E1,F1)
	G=f_symmquartic(r,theta,A2,B2,C2,D2,E2,F2)
	return Q*G
	
def f_symmquartic_prod_x(r,theta,A1,B1,C1,D1,E1,F1,A2,B2,C2,D2,E2,F2):
	Q=f_symmquartic(r,theta,A1,B1,C1,D1,E1,F1)
	Q_x=f_symmquartic_x(r,theta,A1,B1,C1,D1,E1,F1)
	G=f_symmquartic(r,theta,A2,B2,C2,D2,E2,F2)
	G_x=f_symmquartic_x(r,theta,A2,B2,C2,D2,E2,F2)
	return Q*G_x+Q_x*G

def f_symmquartic_prod_y(r,theta,A1,B1,C1,D1,E1,F1,A2,B2,C2,D2,E2,F2):
	Q=f_symmquartic(r,theta,A1,B1,C1,D1,E1,F1)
	Q_y=f_symmquartic_y(r,theta,A1,B1,C1,D1,E1,F1)
	G=f_symmquartic(r,theta,A2,B2,C2,D2,E2,F2)
	G_y=f_symmquartic_y(r,theta,A2,B2,C2,D2,E2,F2)
	return Q*G_y+Q_y*G

def f_symmquartic_prod_drdtheta(r,theta,A1,B1,C1,D1,E1,F1,\
        A2,B2,C2,D2,E2,F2):
	dHdx=f_symmquartic_prod_x(r,theta,A1,B1,C1,D1,E1,F1,\
        A2,B2,C2,D2,E2,F2)
	dHdy=f_symmquartic_prod_y(r,theta,A1,B1,C1,D1,E1,F1,\
        A2,B2,C2,D2,E2,F2)
	return -r*(dHdx*np.cos(theta)-dHdy*np.sin(theta))/\
        (dHdx*np.sin(theta)+dHdy*np.cos(theta))




###########################################


def f_v_ELR(theta,omega,psi,r):
    """Reference: 2011A&A...533A..43E"""
    toler=0.0000001 ### This precision usually requires 
                    ### 4 or 5 iterations.
    theta=go_first_quadrant(theta)
    if theta == 0.0: v=theta; return v
    elif theta == 0.5*np.pi: v=theta; return v
    else:
        auxi=8./81.*psi**3.*omega**2.*r**3.*np.cos(theta)**3.+\
            np.cos(theta)+np.log(abs(np.tan(theta/2.)))
        xprev=theta
        yprev=-np.cos(xprev)-np.log(abs(np.tan(xprev/2.)))+auxi
        yyprev=np.sin(xprev)-1./np.sin(xprev)
        xnew=xprev-yprev/yyprev
        N=1
        while abs(xnew-xprev) > toler:
            yprev=-np.cos(xnew)-np.log(abs(np.tan(xnew/2.)))+auxi
            yyprev=yyprev=np.sin(xnew)-1./np.sin(xnew)
            xprev=xnew
            xnew=xprev-yprev/yyprev
            N=N+1
        v=xnew; return v
    
def f_F_ELR(theta,v,omega,psi,r):
    """Reference: 2011A&A...533A..43E"""
    theta=go_first_quadrant(theta)
    if theta == 0.0: F=np.exp(16./81.*psi**3.*omega**2.*r**3.)
    elif theta == 0.5*np.pi:
        if omega != 1.0: F=(1.-8./27.*psi**3.*omega**2.*r**3.)**(-2./3.)
        else: F=0.0
    else: F=(np.tan(v)/np.tan(theta))**2.
    return F
    

def f_flux_vector_ELR(theta,omega,psi,lum_fac,r):
    """Reference: 2011A&A...533A..43E"""
    v=f_v_ELR(theta,omega,psi,r)
    F=f_F_ELR(theta,v,omega,psi,r)
    effgrav_vec=f_eff_gravity_vector(theta,omega,psi,r)
    flux_r=-lum_fac*F*effgrav_vec[0]
    flux_theta=-lum_fac*F*effgrav_vec[1]
    return np.array([flux_r,flux_theta])

def f_flux_ELR(theta,omega,psi,lum_fac,r):
    """Reference: 2011A&A...533A..43E"""
    #v=f_v_ELR(theta,omega,psi,r)
    #F=f_F_ELR(theta,v,omega,psi,r)
    #effgrav=f_eff_gravity(theta,omega,psi,r)
    #flux=lum_fac*F*effgrav
    flux_vec=f_flux_vector_ELR(theta,omega,psi,lum_fac,r)
    
    return np.sqrt(flux_vec[0]**2.+flux_vec[1]**2.)

def f_Teff_ELR(theta,omega,psi,lum_fac,r):
    """Reference: 2011A&A...533A..43E"""
    flux=f_flux_ELR(theta,omega,psi,lum_fac,r)
    
    return flux**0.25






####
def Tmean_ELR(omega,psi,lum_fac,npts=18001):

    thetamin=0.0*np.pi; thetamax=0.5*np.pi
    theta=np.array([thetamin+(thetamax-thetamin)*\
        float(i)/float(npts-1) for i in xrange(0,npts)])
    rocheradius=np.array([f_rocheradius(theta[i],omega,psi) \
        for i in xrange(0,len(theta))])
    oness=np.array([1. for i in xrange(0,len(theta))])
    flux_nc=np.array([f_flux_ELR(theta[i],omega,psi,lum_fac,\
        rocheradius[i],conservation=False) \
        for i in xrange(0,len(theta))])
    sw=f_integral_surfaceroche(theta,omega,psi,oness)
    Tmean=lum_fac**0.25/sw*f_integral_surfaceroche(theta,omega,psi,\
        flux_nc**0.25)/f_integral_surfaceroche(theta,omega,psi,\
        flux_nc)**0.25

    return Tmean




def f_dvdtheta_roche_ELR(theta,v,omega,psi):
    """Reference: 2011A&A...533A..43E"""
    theta=go_first_quadrant(theta)
    if theta == 0.0: dvdtheta=1.; return dvdtheta
    elif theta == 0.5*np.pi: dvdtheta=np.nan; return dvdtheta
    rocheradius=f_rocheradius(theta,omega,psi)
    drochedtheta=f_drochedtheta(theta,omega,psi)
    dvdtheta_roche=(8./27.*psi**3.*omega**2.*rocheradius**2.*\
        np.cos(theta)**2.*\
        (np.cos(theta)*drochedtheta-\
        rocheradius*np.sin(theta))-np.sin(theta)+\
        1./np.sin(theta))/\
        (-np.sin(v)+1./np.sin(v))
    return dvdtheta_roche

def f_dFdtheta_roche_ELR(theta,v,omega,psi):
    """Reference: 2011A&A...533A..43E"""
    theta=go_first_quadrant(theta)
    if theta == 0.0: 
        dFdtheta_roche=np.nan; return dFdtheta_roche
    elif theta == 0.5*np.pi: 
        dFdtheta_roche=np.nan; return dFdtheta_roche
    else:
        rocheradius=f_rocheradius(theta,omega,psi)
        F=f_F_ELR(theta,v,omega,psi,rocheradius)
        dvdtheta_roche=f_dvdtheta_roche_ELR(theta,v,omega,psi)
        dFdtheta_roche=2.*(dvdtheta_roche/np.sin(v)/np.cos(v)-\
            1./np.sin(theta)/np.cos(theta))*F
        return dFdtheta_roche

def f_beta_ELR(theta,v,omega,psi,estimate=True):
    """Reference: 2011A&A...533A..43E"""
    theta=go_first_quadrant(theta)
    rocheradius=f_rocheradius(theta,omega,psi)
    eff_gravity=f_eff_gravity(theta,omega,psi,rocheradius)
    F=f_F_ELR(theta,v,omega,psi,rocheradius)
    if (theta != 0.0 and theta != 0.5*np.pi):
        dFdtheta=f_dFdtheta_roche_ELR(theta,v,omega,psi)
        deffgravdtheta=f_deffgravdtheta_roche(theta,omega,psi)
        beta=0.25*(1.+(dFdtheta/F)/(deffgravdtheta/eff_gravity))
        return beta
    else:
        if estimate == True:
            deltinha=0.001
            theta2=theta+deltinha; theta3=theta2+deltinha
            theta2=go_first_quadrant(theta2)
            theta3=go_first_quadrant(theta3)
            rocheradius2=f_rocheradius(theta2,omega,psi)
            eff_gravity2=f_eff_gravity(theta2,omega,psi,rocheradius2)
            v2=f_v_ELR(theta2,omega,psi,rocheradius2)
            F2=f_F_ELR(theta2,v2,omega,psi,rocheradius2)
            dFdtheta2=f_dFdtheta_roche_ELR(theta2,v2,omega,psi)
            deffgravdtheta2=f_deffgravdtheta_roche(theta2,omega,psi)
            beta2=0.25*(1.+(dFdtheta2/F2)/\
                (deffgravdtheta2/eff_gravity2))        
            rocheradius3=f_rocheradius(theta3,omega,psi)
            eff_gravity3=f_eff_gravity(theta3,omega,psi,rocheradius3)
            v3=f_v_ELR(theta3,omega,psi,rocheradius3)
            F3=f_F_ELR(theta3,v3,omega,psi,rocheradius3)
            dFdtheta3=f_dFdtheta_roche_ELR(theta3,v3,omega,psi)
            deffgravdtheta3=f_deffgravdtheta_roche(theta3,omega,psi)
            beta3=0.25*(1.+(dFdtheta3/F3)/\
                (deffgravdtheta3/eff_gravity3))
            beta=beta2+(beta3-beta2)/(theta3-theta2)*(theta-theta2)        
            return beta
        else:
            beta=np.nan
            return beta









###########################################
# Rodrigo Vieira's functions

# CONVERT OBLATENESS TO W
def oblat2w(oblat):
    '''
    Converts oblateness into wc=Omega/Omega_crit
    Ekstrom et al. 2008, Eq. 9

    Usage:
    w = oblat2w(oblat)
    '''
    if (np.min(oblat) < 1.) or (np.max(oblat) > 1.5):
        print('Warning: values out of allowed range')

    oblat = np.array([oblat]).reshape((-1))
    nw = len(oblat)
    w = np.zeros(nw)

    for iw in xrange(nw):
        if oblat[iw] <= 1.:
            w[iw] = 0.
        elif oblat[iw] >= 1.5:
            w[iw] = 1.
        else:
            w[iw] = (1.5**1.5) * np.sqrt(2.*(oblat[iw] - 1.) / \
                oblat[iw]**3.)

    if nw == 1:
        w = w[0]

    return w

# CONVERT W TO WK
def w2wk(w):
    '''
    Expressions for beta are in terms of w_k, so
    w_c has to be converted

    In can be done according the relations at
    Espinosa-Lara lecture, slide 18

    w=Omega/Omega_c
    w=w_k=Omega/Omega_c
    w^2=(3/(2+w^2))^3*w^2

    Solving for w with Newton's method
    f=w^2-w_c^2

    Usage:
    wk = w2wk(w)
    '''
    # limit cases
    if w == 0:
        wk = 0.
    elif w == 1.:
        wk = 1.
    elif (w < 0.) or (w > 1.):
        print('w out of allowed range')
        wk = w
    else:
        # initial conditions
        delt = 1.
        wk1 = 0.
        wk = w
        # iteration loop
        while delt >= 1e-5:
            f = (3./(2.+wk**2.))**3. * wk**2. - w**2.
            df = -108. * wk * (wk**2.-1.) / (wk**2.+2.)**4.
            wk1 = wk - f/df
            delt = np.abs(wk1 - wk)/wk
            wk = wk1

    return wk

# BETA ESPINOSA-LARA
def beta_espinosa(w, oblat=False):
    '''
    Beta Espinosa-Lara (VLTI School lecture)

    par = w = Omega/Omega_crit

    OR

    par = oblateness = Re/Rp
    (with /oblat option on)

    Usage:
    bet = beta_espinosa(w, oblat=False)

    w can be either a scalar or an array
    '''
    if oblat:
        w = oblat2w(w)

    # make sure it is a numpy array
    w = np.array([w]).reshape((-1))

    # for w<0.022, poly_fit at the end writes
    # '% POLY_FIT: Warning: Invert detected a small pivot element.'
    w[np.where(w <= 0.025)] = 0.

    # convert wc to wk, used in the beta expressions
    # w=Omega/Omega_crit
    # wk=Omega/Omega_kep
    # (see Espinosa-Lara lecture, slide 18)
    nw = len(w)
    omega = np.array([w2wk(w[iw]) for iw in xrange(nw)])

    # initializing variables
    pi = np.pi
    nthe = 100
    theta = np.linspace(0., pi/2., nthe+1)
    theta = 0.5*(theta[0:-1] + theta[1:]) # to avoid singularities at 0 and 90
    grav = np.zeros(nthe * nw).reshape(nthe, nw)
    teff = np.zeros(nthe * nw).reshape(nthe, nw)
    corr = np.zeros(nthe * nw).reshape(nthe, nw)
    beta = np.zeros(nw)

    # main loop
    for iw in xrange(nw):
        if omega[iw] == 0.:
            beta[iw] = 0.25
        else:
            for ithe in xrange(nthe):

                # R_roche(theta)
                delt = 1.
                r1 = 0.
                r = 1.
                while delt >= 1e-5:
                    f = omega[iw]**2. * r**3. * np.sin(theta[ithe])**2.\
                        -(2.+omega[iw]**2.) * r + 2.
                    df = 3. * omega[iw]**2. * r**2. * \
                        np.sin(theta[ithe])**2. \
                        -(2. + omega[iw]**2.)
                    r1 = r - f/df
                    delt = np.abs(r1 - r) / r
                    r = r1
       
                # eta (Espinosa Lara lecture, slide 30)
                delt = 1.
                n1 = 0.
                ftheta = 1./3. * omega[iw]**2.*r**3. * \
                        np.cos(theta[ithe])**3. \
                        + np.cos(theta[ithe]) + \
                        np.log(np.tan(theta[ithe]/2.))
                n = theta[ithe]
                while delt >= 1e-5:
                    f = np.cos(n) + np.log(np.tan(n/2.)) - ftheta
                    df = -np.sin(n) + 1./np.sin(n)
                    n1 = n - f/df
                    delt = np.abs(n1 - n) / n
                    n = n1
       
                # parameters
                grav[ithe, iw] = np.sqrt(1./r**4. + omega[iw]**4. * \
                                r**2. \
                                * np.sin(theta[ithe])**2. - \
                                2. * omega[iw]**2.
                                * np.sin(theta[ithe])**2. / r)
       
                corr[ithe, iw] = np.sqrt(np.tan(n) / \
                    np.tan(theta[ithe]))
       
                teff[ithe, iw] = corr[ithe, iw] * grav[ithe,iw]**0.25

            # fitting effective beta to the obtained result
            coef = np.polyfit(np.log(grav[:, iw]), \
                np.log(teff[:, iw]), 1)
            beta[iw]=coef[0]

    if nw == 1:
        beta = beta.reshape((-1))[0]

    return beta


# 

def f_beta_ELR_fit(omega,psi,lum_fac):


    # for w<0.022, poly_fit at the end writes
    # '% POLY_FIT: Warning: Invert detected a small pivot element.'
    if omega <= 0.025:
        omega=0.


    if omega == 0.:
        beta=0.25
    else:
        # initializing variables
        nthe = 100
        theta = np.linspace(0., np.pi/2., nthe+1)
        theta = 0.5*(theta[0:-1] + theta[1:]) # to avoid singularities at 0 and 90
    
        rocheradius=np.array([f_rocheradius(theta[i],omega,psi) \
            for i in xrange(0,nthe)])
        v=np.array([f_v_ELR(theta[i],omega,psi,rocheradius[i]) \
            for i in xrange(0,nthe)])

        grav=np.array([f_eff_gravity(theta[i],omega,psi,\
            rocheradius[i]) for i in xrange(0,nthe)])
        teff=np.array([f_Teff_ELR(theta[i],omega,psi,lum_fac,\
            rocheradius[i]) for i in xrange(0,nthe)])

        # fitting effective beta to the obtained result
        coef = np.polyfit(np.log(grav), np.log(teff), 1)
        beta=coef[0]

    return beta





###########################################



def ctes1(mass,r0,lum):
    
    g0=np.nan; angularvel0=np.nan; potential0=np.nan; surface0=np.nan; T0=np.nan
    
    if not np.isnan(mass) and not np.isnan(r0):
        g0=10.0**(np.log10(phc.G.cgs)+np.log10(mass*phc.Msun.cgs)-2.0*np.log10(r0*phc.Rsun.cgs))
        angularvel0=np.sqrt(8./27.*10.0**(np.log10(phc.G.cgs)+np.log10(mass*phc.Msun.cgs)-3.0*np.log10(r0*phc.Rsun.cgs)))
        potential0=-10.0**(np.log10(phc.G.cgs)+np.log10(mass*phc.Msun.cgs)-1.0*np.log10(r0*phc.Rsun.cgs))
        surface0=4.*np.pi*(r0*phc.Rsun.cgs)**2.
    
    if not np.isnan(lum) and not np.isnan(r0):
        T0=(10.**(np.log10(lum/lum_fac*phc.Lsun.cgs)-np.log10(phc.sigma.cgs)-np.log10(surface0)))**0.25
        
    return g0,angularvel0,potential0,surface0,T0
    
    
###########################################






#READ model parameters
#ob,omega,gamma,W=rocheparams(0.7,'W')
#psi=1.0
#lum_fac=1.0
#
#if 1==2:    # 
#            #
#    thetamin=0.0*np.pi; thetamax=0.5*np.pi; npts=20
#    theta=np.array([thetamin+(thetamax-thetamin)*float(i)/float(npts-1) for i in xrange(0,npts)])
#    rocheradius=np.array([f_rocheradius(theta[i],omega,psi) for i in xrange(0,len(theta))])
#    eff_gravity=np.array([f_eff_gravity(theta[i],omega,psi,rocheradius[i]) for i in xrange(0,len(theta))])
#    for i in xrange(0,len(theta)):
#        print theta[i]/np.pi*180.,rocheradius[i],eff_gravity[i]
#    print ''
#
#if 1==2:
#    v=np.array([f_v_ELR(theta[i],omega,psi,rocheradius[i]) for i in xrange(0,len(theta))])
#    F=np.array([f_F_ELR(theta[i],v[i],omega,psi,rocheradius[i]) for i in xrange(0,len(theta))])
#    flux_ELR=np.array([f_flux_ELR(theta[i],omega,psi,lum_fac,rocheradius[i]) for i in xrange(0,len(theta))])
#    efftemp_ELR=np.array([flux_ELR[i]**0.25 for i in xrange(0,len(theta))])
#    #efftemp_betacte=np.array([f_efftemp_betacte(theta[i],omega,psi,0.25) for i in xrange(0,len(theta))])
#    beta_ELR=np.array([f_beta_ELR(theta[i],v[i],omega,psi) for i in xrange(0,len(theta))])
#    for i in xrange(0,len(theta)):
#        #print theta[i]/np.pi*180.,F[i],eff_gravity[i],flux_ELR[i],efftemp_ELR[i],efftemp_ELR[i]*((1.+0.5*W*W)/psi)**0.5
#        print theta[i]/np.pi*180.,F[i],eff_gravity[i],flux_ELR[i],efftemp_ELR[i],beta_ELR[i]
#    print ''


#import sys; sys.exit()
#
#import pyhdust.phc as phc
#
#
#
#
#f0 = open('./Be_M15.00_ob1.33_H0.30_Z0.002_bE_Ell.txt','r')
#lines = f0.readlines()
#f0.close()
#
#psi=1.0
#lum_fac=1.0
#mass=float(lines[3].split()[2])
#r0=float(lines[4].split()[2])
#W=float(lines[5].split()[2])
#lum=float(lines[6].split()[2])
#betaa=float(lines[7].split()[2])
#
#ob,omega,gamma,W=rocheparams(W,'W')
#
#    
#def ctes(mass,r0,lum):
#    
#    g0=np.nan; angularvel0=np.nan; potential0=np.nan; surface0=np.nan; T0=np.nan
#    
#    if not np.isnan(mass) and not np.isnan(r0):
#        g0=10.0**(np.log10(phc.G.cgs)+np.log10(mass*phc.Msun.cgs)-2.0*np.log10(r0*phc.Rsun.cgs))
#        angularvel0=np.sqrt(8./27.*10.0**(np.log10(phc.G.cgs)+np.log10(mass*phc.Msun.cgs)-3.0*np.log10(r0*phc.Rsun.cgs)))
#        potential0=-10.0**(np.log10(phc.G.cgs)+np.log10(mass*phc.Msun.cgs)-1.0*np.log10(r0*phc.Rsun.cgs))
#        surface0=4.*np.pi*(r0*phc.Rsun.cgs)**2.
#    
#    if not np.isnan(lum) and not np.isnan(r0):
#        T0=(10.**(np.log10(lum/lum_fac*phc.Lsun.cgs)-np.log10(phc.sigma.cgs)-np.log10(surface0)))**0.25
#        
#    return g0,angularvel0,potential0,surface0,T0
#
#
#
#
#g0,angularvel0,potential0,surface0,T0=ctes(mass,r0,lum)
#vpar=angularvel0*r0*(27./8.)**0.5
#vorb=vpar*f_rocheradius(0.5*np.pi,omega,psi)**-0.5
#vcrit=vorb*(2./3.*ob)**0.5
#veq=W*vorb
#Prot=2.*np.pi/omega*psi**-1.5*angularvel0**-0.5
#print vpar,vorb,vcrit,veq,Prot
#
#
#
#
#print psi
#print lum_fac
#print mass
#print r0
#print W
#print lum
#print betaa
#b=Tmean_ELR(omega,psi,lum_fac)
#cte=ctes(mass,r0,lum)
#c=Tmean4_qu(omega,psi,lum_fac)
#print 0.6*c*cte[4]
#print 0.6*b*cte[4]
#print ''
#
#
#
#
#import sys; sys.exit()
#
#
#npts=18001; thetamin=0.0; thetamax=1.0*np.pi
#theta[i]=np.array([thetamin+(thetamax-thetamin)*float(i)/float(npts-1) for i in xrange(0,npts)])
#rocheradius[i]=np.array([f_rocheradius(theta[i],omega,psi) for i in xrange(0,npts)])
#gravityroche[i]=np.array([f_eff_gravity(theta[i],omega,psi,rocheradius[i]) for i in xrange(0,npts)])
#
#
#
#
#
#vcrit=10.0**(0.5*(np.log10(gconst)+np.log10(mass*solmas)-np.log10(1.5*rocheradius[0])-np.log10(r0*solrad)))
#vorb =10.0**(0.5*(np.log10(gconst)+np.log10(mass*solmas)-np.log10(rocheradius[int(float(npts-1)/2.0)])-np.log10(r0*solrad)))
#prot =10.0**(np.log10(2.0*np.pi)-np.log10(omega)-0.5*np.log10(8./27.*gconst)-0.5*np.log10(mass*solmas)+1.5*np.log10(XXXX*rocheradius[int(float(npts-1)/2.0)]*solrad))
#
#da=zeros(npts-1)
#for i in range(0,npts-1):
# if theta[i]!=0.5*pi:
#  da[i]=2.0*pi*sqrt((rocheradius[i]**2.0*sin(theta[i]))**2.0+(rocheradius[i]*sin(theta[i])*((8./27.)*s_omega**2.0*s_psi**3.0*rocheradius[i]**2.0*sin(theta[i])*cos(theta[i]))/(rocheradius[i]**(-2.0)-(8./27.)*s_omega**2.0*s_psi**3.0*rocheradius[i]*sin(theta[i])**2.0))**2.0)*(theta[i+1]-theta[i])
# if theta[i]==0.5*pi:
#  da[i]=2.0*pi*sqrt((rocheradius[i]**2.0*sin(theta[i]))**2.0)*(theta[i+1]-theta[i])
#da[:]=10.0**(log10(da[:])+2.0*log10(src[iz,imass,itau,iob,1]*solrad))
#
#b_theta=zeros(npts)
#if b_law==0:
# for i in range(0,npts):
#  b_theta[i]=src[iz,imass,itau,iob,4]
#
#areag=0.0; area=0.0
#for i in range(0,npts-1):
# areag=areag+10.0**((4.0*b_theta[i])*(log10(gp0)+log10(gravity[i]))+log10(da[i]))
# area=area+10.0**(log10(da[i]))
#teff[iz,imass,itau,iob]=10.**(0.25*log10(src[iz,imass,itau,iob,3]*sollum)-0.25*log10(sigmac)-0.25*log10(area))
#print teff[iz,imass,itau,iob]











