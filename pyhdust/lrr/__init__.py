# -*- coding:utf-8 -*-

import sys
import os
import numpy as np
import itertools as it
import pyhdust.phc as phc

__author__ = "Leandro Rimulo"
__email__ = "lrrimulo@gmail.com"

cwd = os.path.dirname(os.path.realpath(__file__))
cwd+='/'


def makeitso(x=None):
    """
    A good function to run at the beginning of your program! :)
    """
    
    N = 5 
    if x == None:
        x = np.random.random_integers(1,N)
    
    if x == 1:
        print("\n MAKE IT SO! \n")
        print("                                _____")
        print("                       __...---'-----`---...__")
        print("                  _===============================")
        print(" ______________,/'      `---..._______...---'")
        print("(____________LL). .    ,--'")
        print(" /    /.---'       `. /")
        print("'--------_  - - - - _/")
        print("          `~~~~~~~~'")
        print("")
    if x == 2:
        print("\n ENGAGE! \n")
        print("               .")
        print("              .:.")
        print("             .:::.")
        print("            .:::::.")
        print("        ***.:::::::.***")
        print("   *******.:::::::::.*******")       
        print(" ********.:::::::::::.********")     
        print("********.:::::::::::::.********")    
        print("*******.::::::'***`::::.*******")    
        print("******.::::'*********`::.******")    
        print(" ****.:::'*************`:.****")
        print("   *.::'*****************`.*")
        print("   .:'  ***************    .")
        print("  .")
        print("")
    if x == 3:
        print("\n I'M GIVIN' HER ALL SHE'S GOT, CAPTAIN! \n")
        print("__________________           __")
        print("\_________________|)____.---'--`---.____")
        print("              ||    \----.________.----/")
        print("              ||     / /    `--'")
        print("            __||____/ /_")
        print("           |___         \ ")
        print("               `--------'")
        print("")
    if x == 4:    
        print("\n WARP FACTOR 9, MISTER SULU! \n")
        print("       .----------.___")
        print("     / |        |||\  ~~~--_               ____")
        print("  __|-------------| |~~~~~~~~----______--~~    ~~-_")
        print(" | _|-_           |_|                  =============")
        print(" |_|   ~~---------'__\_____------~~~~~~--.______.-~")
        print("  |__|~~ ~-__/   ~-_~-_")
        print("     ~--__   ~-__/  ~-_~-_")
        print("          ~--__  ~-__/ ~-_~-_")
        print("               ~--__ ~-___/  ~-_")
        print("                     ~--__       ~-_")
        print("                        [__________]=====[[")
        print("")
    if x == 5:    
        print("\n THE HUMAN ADVENTURE IS JUST BEGINNING! \n")        
        print("           ______")
        print("        _-' .   .`-_")
        print("    |/ /  .. . '   .\ \|")
        print("   |/ /            ..\ \|")
        print(" \|/ |: .   ._|_ .. . | \|/")
        print("  \/ |   _|_ .| . .:  | \/")
        print(" \ / |.   |  .  .    .| \ /")
        print("  \||| .  . .  _|_   .|||/")
        print(" \__| \  . :.  .|.  ./ |__/")
        print("   __| \_  .    .. _/ |__")
        print("    __|  `-______-'  |__")
        print("       -,____  ____,-")
        print("         ---'  `---")
        print("")


    return



def integrate_trapezia(f,dg):
    """
    Integration by trapezia, with differentials of possibly 
    different size:
             ____
             \
    Returns: /___ f dg
    
    INPUT: lists/arrays of values of f and dg (Obs: len(f)=len(dg)+1)
    """
    if len(f) < 2 or len(dg) != len(f)-1:
        func_name = sys._getframe().f_code.co_name
        print '<<',func_name,'>>'
        print 'There is something wrong with the integration!'
        print ''
        return np.nan
    else: 
        soma = 0.0
        for i in xrange(0,len(f) - 1):
            termo = float((f[i] + f[i + 1]) * dg[i])/2.
            if not np.isnan(termo):
                soma += termo
        return soma
        
def dec_2_binary(number,dim=None):
    """Returns the number in its binary form.
    
    Ex. 
    >>> lrr.dec_2_binary(5)
    array([1, 0, 1])
    >>> lrr.dec_2_binary(5,5)
    array([0, 0, 1, 0, 1])
    >>> lrr.dec_2_binary(5,2)
    << dec_2_binary >>
    I need more dimensions to write this binary number!
    
    nan
    """

    if number != int(number) or number < 0:
        func_name = sys._getframe().f_code.co_name
        print '<<',func_name,'>>'
        print 'The number must be a natural number.'
        print ''
        return np.nan

    addreverse=[]; res=number%2; quo=number/2
    addreverse.append(res)
    while quo != 0: 
        res=quo%2; quo=quo/2; addreverse.append(res)
    
    if dim != None and dim < len(addreverse):
        func_name = sys._getframe().f_code.co_name
        print '<<',func_name,'>>'
        print 'I need more dimensions to write this binary number!'
        print ''
        return np.nan
    elif dim != None and dim >= len(addreverse): 
        addn=np.array([0 for k in xrange(0,dim)])
    else: 
        addn=np.array([0 for k in xrange(0,len(addreverse))])
    
    for j in xrange(0,len(addreverse)): 
        addn[len(addn)-1-j]=addreverse[j]

    return addn
    
def dice(minx=None,maxx=None,p='no'):
    """Returns a random integer number between 
    two boundaries (inclusive)."""
    
    between=[np.nan,np.nan]


    if minx==None or maxx==None:
        while True:
            try:
                between[0] = int(raw_input("I will give you an integer number between: "))
                between[1] = int(raw_input("and: "))
                teste = np.sqrt(float(between[1])-float(between[0]))
                break
            except ValueError:
                print "Oops!  That was no valid dominium.  Try again..."
        x=np.random.randint(between[0],between[1]+1)
        print " ==> Here is your number: "+str(x)
        return x

    else:
        while True:
            try:
                between[0] = int(minx)
                between[1] = int(maxx)
                teste = np.sqrt(float(between[1])-float(between[0]))
                break
            except ValueError:
                print "Oops!  That was no valid dominium.  Try again..."
        x=np.random.randint(between[0],between[1]+1)
        if p == 'yes': print " ==> Here is your number: "+str(x)
        return x      
        
def find_interval(X,a,interval,ind):
    """ Finds the interval, among the elements of a, 
    in which X is located.
    X is a number. a is a list of numbers.
    >>> Please, insert a in ascending order! <<<
    Returns: "X belongs to [ a[ind],a[ind+1] )" & ind
    """

    # Empty and unitary list
    if len(a)==0:
        interval=[-np.inf,np.inf]
        ind=np.nan
        return interval,ind
    elif len(a)==1: 
        if X < a[0]:
            interval=[-np.inf,a[0]]
        else:
            interval=[a[0],np.inf]
        ind=np.nan
        return interval,ind
    # Now, test if X is outside the boundaries of the list.
    elif X < a[0]: # 
        interval=[-np.inf,a[0]]
        ind=-np.inf
        return interval,ind
    elif X >= a[len(a)-1]:
        interval=[a[len(a)-1],np.inf]
        ind=len(a)-1
        return interval,ind
    # Now that I know that X is between the boundaries of the list: 
    else:
        i=0
        while not (X >= a[i] and X < a[i+1]):
            i=i+1
        interval=[a[i],a[i+1]]
        ind=i
        return interval,ind        




def logsumexp_trick(expon):
    """
    Returns the exponent of the single exponential that corresponds 
    to the sum of exponentials:
    logexpsum = log( sum( e^expon[i] ) ) 
    
    ...using the log-sum-exp trick, to avoid rounding problems.
    """
    
    M=len(expon)
    if M > 0:
        expmax=np.nanmax(expon)
        soma=0.
        for elem in expon:
            if not np.isnan(elem):
                soma+=np.exp(elem-expmax)
        logsumexp=expmax+np.log(soma)
        return logsumexp
    else: np.nan



def round_sig(x, sig=3):
    """
    Returns the value of 'x' with 'sig' significant numbers.
    """
    return round(x, sig-int(np.floor(np.log10(abs(x))))-1)





##########################################

def scale_log10(x,k1,k2,m="normal"):

    if m != "inverse":
        if (not np.isnan(x)) and x>0.: return k1*np.log10(k2*x)
        else: return np.nan
    else: return 1./k2*10.**(x/k1)

def scale_log(x,k1,k2,m="normal"):

    if m != "inverse":
        if (not np.isnan(x)) and x>0.: return k1*np.log(k2*x)
        else: return np.nan
    else: return 1./k2*np.exp(x/k1)

def scale_two_propto(x,up1,down1,m="normal"):

    if m != "inverse":
        if x >= 0: return up1*x
        if x < 0: return down1*x
    else:
        if x >= 0: return x/up1
        if x < 0: return x/down1

def scale_two_arcsinh(x,up1,up2,down1,down2,m="normal"):
    """
    Transforms a real number x into the real number y given by
    y = A1*arcsinh(A2*x) = A1*ln[A2*x + ((A2*x)^2 + 1)^0.5]
    
    Good guidelines for choosing the quantity A2 are the following properties:
    . If |A2*x| >> 1: y ~ sign(A2*x)*A1*ln[2*|A2*x|]
    . If |A2*x| << 1: y ~ A1*A2*x

    The inverse of the transformation is given by 
    x = 1/A2 * sinh(y/A1) = 1/A2 * 0.5*[e^(y/A1)-e^(-y/A1)]
    """
    
    if m != "inverse":
        if x >= 0: return up1*np.arcsinh(x*up2)
        if x < 0: return down1*np.arcsinh(x*down2)
    else:
        if x >= 0: return 1./up2*np.sinh(x/up1)
        if x < 0: return 1./down2*np.sinh(x/down1)

def scale_arctan(x,k1,k2,m="normal"):
    
    if m != "inverse": return 2./np.pi*k1*np.arctan(x*k2)
    else: return 1./k2*np.tan(np.pi/2.*x/k1)

def scale_powerlaw(x,up1,up2,down1,down2,m="normal"):
    
    if m != "inverse":
        if (not np.isnan(x)) and not (x == 0. and up2<0.):
            if x >= 0: return up1*x**up2
            if x < 0: return down1*x**down2
    else:
        if (not np.isnan(x)) and not (x == 0. and up2<0.):
            if x >= 0: return (x/up1)**(1./up2)
            if x < 0: return (x/down1)**(1./down2)





    


##########################################

def interLinND(X, X0, X1, Fx, tp="linear"):
    """
    N-dimensional linear interpolation. (Based on Moser's function
    with the same name.)

    | INPUT:
    | X = list containing the position in with the interpolation 
    |     is desired;
    | X0 = list containing minimal values of the interval;
    | X1 = list containing maximum values of the inveral
    | Fx = list containing function values along the interval, 
    |     ORDERED BY DIMENSION.
    |   Example: Fx = [F00, F01, F10, F11]
    |   Example: Fx = [F000, F001, F010, F011, F100, F101, F110, F111]

    OUTPUT: interpolated value (float)"""
    X = np.array(X)     ### an array of N positions
    X0 = np.array(X0)   ### an array of N backward positions
    X1 = np.array(X1)   ### an array of N forward positions
    Xd = (X-X0)/(X1-X0) ### an array of N normalized positions
    DX = np.array([ [(1-x),x] for x in Xd ]) ### an array containing N 
                                             ### "2-arrays" 
                                             ### of backward and forward 
                                             ### normalized intervals
    ### OBS: Fx is an array of 2**N elements correctly arranged
    i = 0
    F = 0
    for prod in it.product(*DX):    ### itertools.product(*DX) creates 
                                    ### 2**N tuples of 
                                    ### N elements => the set of tuples 
                                    ### encompasses all the 
                                    ### possibilities of combinations 
                                    ### of normalized intervals.
        ### np.product returns the product 
        ### of the N elements of each
        ### of the 2**N tuples
        if tp == "linear":
            F+= Fx[i]*np.product(prod) 
        elif tp == "ln":
            F+= np.log(Fx[i])*np.product(prod) 
        elif tp == "arcsinh":
            F+= np.arcsinh(Fx[i])*np.product(prod) 

        i+= 1
    #
    if tp == "linear":
        return F
    elif tp == "ln":
        return np.exp(F)
    elif tp == "arcsinh":
        return np.sinh(F)

    
def low_high(X,axis):

    low=[]; high=[]; ind=[]
    extrapolated="no"
    for j in range(0,len(axis)):
        chave=0
        ### If it is to extrapolate to a value before a minimum
        if X[j]<axis[j][0]:
            ind.append(0)
            low.append(axis[j][0])
            high.append(axis[j][1])
            chave+=1
            extrapolated="yes"
        ### If it is to extrapolate to a value after a maximum
        if X[j]>axis[j][len(axis[j])-1]:
            ind.append(len(axis[j])-2)
            low.append(axis[j][len(axis[j])-2])
            high.append(axis[j][len(axis[j])-1])
            chave+=1
            extrapolated="yes"
        ### If it is to interpolate
        i=0
        while chave==0:
            if X[j] >= axis[j][i] \
                    and X[j] <= axis[j][i+1] \
                    and i < len(axis[j][:])-1:
                ind.append(i)
                low.append(axis[j][i])
                high.append(axis[j][i+1])
                chave+=1
            i+=1

    return np.array(low),np.array(high),np.array(ind),extrapolated

def find_index(ind,axis):

    leng=[len(x) for x in axis]
    indd=0
    for i in xrange(0,len(axis)-1):
        prod=1
        for j in xrange(i+1,len(axis)):
            prod=prod*leng[j]
        indd+=ind[i]*prod
    indd+=ind[len(axis)-1]
    
    return indd
    

def build_Fx(axis,values,ind):
    """
    This function builds the Fx array, to be used in the interpolation 
    function 'interLinND'. 
    
    It uses the grid's domain defined by 'axis' 
    and the 'values' contained in each element of the grid. With the 
    index 'ind', it finds the 2**len(ind) elements of Fx in the correct 
    order. 
    """

    Fx=[]
    N_Fx=2**len(ind)
    for i in xrange(0,N_Fx):
        addn=dec_2_binary(i,len(ind))
        #print("addn = ",addn)
        newind=np.array([ind[j]+addn[j] for j in xrange(len(ind))])
        #print("newind = ",newind)
        indd=find_index(newind,axis)
        #print("indd = ",indd)
        Fx.append(values[indd])
    Fx=np.array(Fx)
    return Fx
    
def interpLinND(X,axis,values,tp="linear",allow_extrapolation="yes"):
    """N-dimensional linear interpolation.
    
    | INPUT:
    | 'X' = the position in which the interpolation is desired
    | 'axis' = list of vectors that create domain of the grid
    | 'values' = vector with all the values associated with each 
    |           element of the grid domain.
    |
    | OUTPUT:
    | 'interp' = interpolated value
    """
    
    if len(X) != len(axis):
        func_name = sys._getframe().f_code.co_name
        print '<<',func_name,'>>'
        print 'X and axis must have the same dimension!'
        print ''
        return np.nan
    low,high,ind,extrapolated=low_high(X,axis)
    if allow_extrapolation != "yes" and extrapolated == "yes":
        return np.nan
    else:
        Fx=build_Fx(axis,values,ind)
        interp=interLinND(X, low, high, Fx, tp)
        return interp








######################################################################



def photosystem_k(filtername):
    if filtername == 'STMAG': return 0.00
    elif filtername == 'white': return 0.00
    ### Reference: 1998A&A...333..231B
    elif filtername == 'bess-u': return 0.03 
    ### Reference: 1998A&A...333..231B
    elif filtername == 'bess-b': return 0.03 
    ### Reference: 1998A&A...333..231B
    elif filtername == 'bess-v': return 0.03 
    ### Reference: 1998A&A...333..231B
    elif filtername == 'bess-r': return 0.03 
    ### Reference: 1998A&A...333..231B
    elif filtername == 'bess-i': return 0.03 
    ### Reference: 1998A&A...333..231B
    elif filtername == 'bess-j': return 0.03 
    ### Reference: 1998A&A...333..231B
    elif filtername == 'bess-h': return 0.03 
    ### Reference: 1998A&A...333..231B
    elif filtername == 'bess-k': return 0.03 
    ### Reference: 1998A&A...333..231B
    elif filtername == 'bess-l': return 0.03 
    ### Reference: 1998A&A...333..231B
    elif filtername == 'bess-ll': return 0.03 
    ### Reference: 1998A&A...333..231B
    elif filtername == 'bess-m': return 0.03 
    ### Reference: 2003AJ....126.1090C
    elif filtername == 'cohen-j': return 0.00 
    ### Reference: 2003AJ....126.1090C
    elif filtername == 'cohen-h': return 0.00 
    ### Reference: 2003AJ....126.1090C
    elif filtername == 'cohen-k': return 0.00 
    ### Reference: 1998AJ....116..482G
    elif filtername == 'crawford-stromgren-u': return 1.445 
    ### Reference: 1998AJ....116..482G
    elif filtername == 'crawford-stromgren-v': return 0.195 
    ### Reference: 1998AJ....116..482G
    elif filtername == 'crawford-stromgren-b': return 0.034 
    ### Reference: 1998AJ....116..482G
    elif filtername == 'crawford-stromgren-y': return 0.030 
    elif filtername == 'int_wfc-hbetan': return 9999.
    elif filtername == 'int_wfc-hbetaw': return 9999.
    elif filtername == 'int_wfc-halpha': return 9999.
    ### Reference: http://www.ctio.noao.edu/soar/content/filters-available-soar
    elif filtername == 'filt_Ha': return 9999.
    else: 
        func_name = sys._getframe().f_code.co_name
        print '<<',func_name,'>>'
        print filtername,' ?? I do not know which filter is that!'
        print ''
        return np.nan


def filter_passband(filtername):
    """
    Returns, for the filter 'filtername'
    * array of lambdas [Angstroms]
    * array of the function R_X
    """
    
    if filtername == 'STMAG':
        func_name = sys._getframe().f_code.co_name
        #print '<<',func_name,'>>'
        #print 'Warning! A passband for STMAG was requested! I am returning NaN to the calculation.'
        #print ''
        return np.nan,np.nan
    elif filtername == 'white':
        return np.array([0.,10.**99.]),np.array([1.,1.])
    elif filtername == 'bess-u': ### Reference: 1990PASP..102.1181B (Table 2)
        f0 = open(cwd+'../refs/filters/Bessell_U.dat','r')
        lines = f0.readlines()
        f0.close()
        ### Eliminating text and separating columns
        lines=[lines[i].split() for i in xrange(1,len(lines))]
        ### lambda (Angstroms)
        lamb=np.array([float(lines[i][0]) for i in xrange(0,len(lines))])
        ### Reponse function
        R=np.array([float(lines[i][1]) for i in xrange(0,len(lines))])
        return lamb,R
    elif filtername == 'bess-b': ### Reference: 1990PASP..102.1181B (Table 2)
        f0 = open(cwd+'../refs/filters/Bessell_B.dat','r')
        lines = f0.readlines()
        f0.close()
        ### Eliminating text and separating columns
        lines=[lines[i].split() for i in xrange(1,len(lines))]
        ### lambda (Angstroms)
        lamb=np.array([float(lines[i][0]) for i in xrange(0,len(lines))])
        ### Reponse function
        R=np.array([float(lines[i][1]) for i in xrange(0,len(lines))])
        return lamb,R
    elif filtername == 'bess-v': ### Reference: 1990PASP..102.1181B (Table 2)
        f0 = open(cwd+'../refs/filters/Bessell_V.dat','r')
        lines = f0.readlines()
        f0.close()
        ### Eliminating text and separating columns
        lines=[lines[i].split() for i in xrange(1,len(lines))]
        ### lambda (Angstroms)
        lamb=np.array([float(lines[i][0]) for i in xrange(0,len(lines))])
        ### Reponse function
        R=np.array([float(lines[i][1]) for i in xrange(0,len(lines))])
        return lamb,R
    elif filtername == 'bess-r': ### Reference: 1990PASP..102.1181B (Table 2)
        f0 = open(cwd+'../refs/filters/Bessell_R.dat','r')
        lines = f0.readlines()
        f0.close()
        ### Eliminating text and separating columns
        lines=[lines[i].split() for i in xrange(1,len(lines))]
        ### lambda (Angstroms)
        lamb=np.array([float(lines[i][0]) for i in xrange(0,len(lines))])
        ### Reponse function
        R=np.array([float(lines[i][1]) for i in xrange(0,len(lines))])
        return lamb,R
    elif filtername == 'bess-i': ### Reference: 1990PASP..102.1181B (Table 2)
        f0 = open(cwd+'../refs/filters/Bessell_I.dat','r')
        lines = f0.readlines()
        f0.close()
        ### Eliminating text and separating columns
        lines=[lines[i].split() for i in xrange(1,len(lines))]
        ### lambda (Angstroms)
        lamb=np.array([float(lines[i][0]) for i in xrange(0,len(lines))])
        ### Reponse function
        R=np.array([float(lines[i][1]) for i in xrange(0,len(lines))])
        return lamb,R
    elif filtername == 'bess-j': ### Reference: 1988PASP..100.1134B (Table IV)
        f0 = open(cwd+'../refs/filters/Bessell_J.dat','r')
        lines = f0.readlines()
        f0.close()
        ### Eliminating text and separating columns
        lines=[lines[i].split() for i in xrange(1,len(lines))]
        ### lambda (Angstroms)
        lamb=np.array([float(lines[i][0]) for i in xrange(0,len(lines))])
        ### Reponse function
        R=np.array([float(lines[i][1]) for i in xrange(0,len(lines))])
        return lamb,R
    elif filtername == 'bess-h': ### Reference: 1988PASP..100.1134B (Table IV)
        f0 = open(cwd+'../refs/filters/Bessell_H.dat','r')
        lines = f0.readlines()
        f0.close()
        ### Eliminating text and separating columns
        lines=[lines[i].split() for i in xrange(1,len(lines))]
        ### lambda (Angstroms)
        lamb=np.array([float(lines[i][0]) for i in xrange(0,len(lines))])
        ### Reponse function
        R=np.array([float(lines[i][1]) for i in xrange(0,len(lines))])
        return lamb,R
    elif filtername == 'bess-k': ### Reference: 1988PASP..100.1134B (Table IV)
        f0 = open(cwd+'../refs/filters/Bessell_K.dat','r')
        lines = f0.readlines()
        f0.close()
        ### Eliminating text and separating columns
        lines=[lines[i].split() for i in xrange(1,len(lines))]
        ### lambda (Angstroms)
        lamb=np.array([float(lines[i][0]) for i in xrange(0,len(lines))])
        ### Reponse function
        R=np.array([float(lines[i][1]) for i in xrange(0,len(lines))])
        return lamb,R
    elif filtername == 'bess-l': ### Reference: 1988PASP..100.1134B (Table IV)
        f0 = open(cwd+'../refs/filters/Bessell_L.dat','r')
        lines = f0.readlines()
        f0.close()
        ### Eliminating text and separating columns
        lines=[lines[i].split() for i in xrange(1,len(lines))]
        ### lambda (Angstroms)
        lamb=np.array([float(lines[i][0]) for i in xrange(0,len(lines))])
        ### Reponse function
        R=np.array([float(lines[i][1]) for i in xrange(0,len(lines))])
        return lamb,R
    elif filtername == 'bess-ll': ### Reference: 1988PASP..100.1134B (Table IV)
        f0 = open(cwd+'../refs/filters/Bessell_LL.dat','r')
        lines = f0.readlines()
        f0.close()
        ### Eliminating text and separating columns
        lines=[lines[i].split() for i in xrange(1,len(lines))]
        ### lambda (Angstroms)
        lamb=np.array([float(lines[i][0]) for i in xrange(0,len(lines))])
        ### Reponse function
        R=np.array([float(lines[i][1]) for i in xrange(0,len(lines))])
        return lamb,R
    elif filtername == 'bess-m': ### Reference: 1988PASP..100.1134B (Table IV)
        f0 = open(cwd+'../refs/filters/Bessell_M.dat','r')
        lines = f0.readlines()
        f0.close()
        ### Eliminating text and separating columns
        lines=[lines[i].split() for i in xrange(1,len(lines))]
        ### lambda (Angstroms)
        lamb=np.array([float(lines[i][0]) for i in xrange(0,len(lines))])
        ### Reponse function
        R=np.array([float(lines[i][1]) for i in xrange(0,len(lines))])
        return lamb,R
    elif filtername == 'cohen-j': ### Reference: 2003AJ....126.1090C
        f0 = open(cwd+'../refs/filters/Cohen-J.pass','r')
        lines = f0.readlines()
        f0.close()
        ### Eliminating text and separating columns
        lines=[lines[i].split() for i in xrange(1,len(lines))]
        ### lambda (Angstroms)
        lamb=np.array([float(lines[i][0]) for i in xrange(0,len(lines))])
        ### Reponse function
        R=np.array([float(lines[i][1]) for i in xrange(0,len(lines))])
        return lamb,R
    elif filtername == 'cohen-h': ### Reference: 2003AJ....126.1090C
        f0 = open(cwd+'../refs/filters/Cohen-H.pass','r')
        lines = f0.readlines()
        f0.close()
        ### Eliminating text and separating columns
        lines=[lines[i].split() for i in xrange(1,len(lines))]
        ### lambda (Angstroms)
        lamb=np.array([float(lines[i][0]) for i in xrange(0,len(lines))])
        ### Reponse function
        R=np.array([float(lines[i][1]) for i in xrange(0,len(lines))])
        return lamb,R
    elif filtername == 'cohen-k': ### Reference: 2003AJ....126.1090C
        f0 = open(cwd+'../refs/filters/Cohen-K.pass','r')
        lines = f0.readlines()
        f0.close()
        ### Eliminating text and separating columns
        lines=[lines[i].split() for i in xrange(1,len(lines))]
        ### lambda (Angstroms)
        lamb=np.array([float(lines[i][0]) for i in xrange(0,len(lines))])
        ### Reponse function
        R=np.array([float(lines[i][1]) for i in xrange(0,len(lines))])
        return lamb,R
    elif filtername == 'crawford-stromgren-u': ### Reference: 1970AJ.....75..978C
        f0 = open(cwd+'../refs/filters/Crawford-Stromgren-U.dat','r')
        lines = f0.readlines()
        f0.close()
        ### Eliminating text and separating columns
        lines=[lines[i].split() for i in xrange(2,len(lines))]
        ### lambda (Angstroms)
        lamb=np.array([float(lines[i][0]) for i in xrange(0,len(lines))])
        ### Reponse function
        R=np.array([float(lines[i][1]) for i in xrange(0,len(lines))])
        return lamb,R
    elif filtername == 'crawford-stromgren-v': ### Reference: 1970AJ.....75..978C
        f0 = open(cwd+'../refs/filters/Crawford-Stromgren-V.dat','r')
        lines = f0.readlines()
        f0.close()
        ### Eliminating text and separating columns
        lines=[lines[i].split() for i in xrange(2,len(lines))]
        ### lambda (Angstroms)
        lamb=np.array([float(lines[i][0]) for i in xrange(0,len(lines))])
        ### Reponse function
        R=np.array([float(lines[i][1]) for i in xrange(0,len(lines))])
        return lamb,R
    elif filtername == 'crawford-stromgren-b': ### Reference: 1970AJ.....75..978C
        f0 = open(cwd+'../refs/filters/Crawford-Stromgren-B.dat','r')
        lines = f0.readlines()
        f0.close()
        ### Eliminating text and separating columns
        lines=[lines[i].split() for i in xrange(2,len(lines))]
        ### lambda (Angstroms)
        lamb=np.array([float(lines[i][0]) for i in xrange(0,len(lines))])
        ### Reponse function
        R=np.array([float(lines[i][1]) for i in xrange(0,len(lines))])
        return lamb,R
    elif filtername == 'crawford-stromgren-y': ### Reference: 1970AJ.....75..978C
        f0 = open(cwd+'../refs/filters/Crawford-Stromgren-Y.dat','r')
        lines = f0.readlines()
        f0.close()
        ### Eliminating text and separating columns
        lines=[lines[i].split() for i in xrange(2,len(lines))]
        ### lambda (Angstroms)
        lamb=np.array([float(lines[i][0]) for i in xrange(0,len(lines))])
        ### Reponse function
        R=np.array([float(lines[i][1]) for i in xrange(0,len(lines))])
        return lamb,R
    elif filtername == 'int_wfc-hbetan':
        f0 = open(cwd+'../refs/filters/INT_WFC-Hbetan.dat','r')
        lines = f0.readlines()
        f0.close()
        ### Eliminating text and separating columns
        lines=[lines[i].split() for i in xrange(0,len(lines))]
        ### lambda (Angstroms)
        lamb=np.array([float(lines[i][0]) for i in xrange(0,len(lines))])
        ### Reponse function
        R=np.array([float(lines[i][1]) for i in xrange(0,len(lines))])
        return lamb,R
    elif filtername == 'int_wfc-hbetaw':
        f0 = open(cwd+'../refs/filters/INT_WFC-Hbetaw.dat','r')
        lines = f0.readlines()
        f0.close()
        ### Eliminating text and separating columns
        lines=[lines[i].split() for i in xrange(0,len(lines))]
        ### lambda (Angstroms)
        lamb=np.array([float(lines[i][0]) for i in xrange(0,len(lines))])
        ### Reponse function
        R=np.array([float(lines[i][1]) for i in xrange(0,len(lines))])
        return lamb,R
    elif filtername == 'int_wfc-halpha':
        f0 = open(cwd+'../refs/filters/INT_WFC-Halpha.dat','r')
        lines = f0.readlines()
        f0.close()
        ### Eliminating text and separating columns
        lines=[lines[i].split() for i in xrange(0,len(lines))]
        ### lambda (Angstroms)
        lamb=np.array([float(lines[i][0]) for i in xrange(0,len(lines))])
        ### Reponse function
        R=np.array([float(lines[i][1]) for i in xrange(0,len(lines))])
        return lamb,R
    elif filtername == 'filt_Ha': ### Reference: http://www.ctio.noao.edu/soar/content/filters-available-soar
        f0 = open(cwd+'../refs/filters/filt_Ha.dat','r')
        lines = f0.readlines()
        f0.close()
        ### Eliminating text and separating columns
        lines=[lines[i].split() for i in xrange(5,len(lines))]
        ### lambda (Angstroms)
        lamb=np.array([float(lines[i][0])*10. for i in xrange(0,len(lines))])
        ### Reponse function
        R=np.array([float(lines[i][1])/100. for i in xrange(0,len(lines))])
        return lamb,R
    else: 
        func_name = sys._getframe().f_code.co_name
        print '<<',func_name,'>>'
        print filtername,' ?? I do not know which filter is that!'
        print ''
        return np.nan,np.nan













def photonflux(lamb,flambda,filtername,npts=50):
    
    ### Obtaining the passband associated with the filter 'filtername'
    lamb_r,R=filter_passband(filtername)
    
    if filtername == 'white':
        
        dens=np.array(lamb)*np.array(flambda)/\
            (phc.h.cgs*phc.c.cgs)   ### lamb*flamb/hc [cgs units]
        
        ### domain of lambdas [Angstroms] for interpolation
        l=np.array([np.nanmin(lamb)+(np.nanmax(lamb)-np.nanmin(lamb))/\
            float(npts-1)*float(i) for i in xrange(0,npts)])
        
        ### variation of lambda [Angstroms] for integration
        dl=np.array([l[i+1]-l[i] for i in xrange(0,npts-1)])
        
        ### interpolated lamb*flamb/hc [cgs units]
        d=np.array([interpLinND([l[i]],[lamb],dens) \
            for i in xrange(0,npts)])
        ### photon flux [s^-1 cm^-2]
        photflux=integrate_trapezia(d,dl*10.**-8.)
        return photflux
    
    elif np.nanmin(lamb_r) >= np.nanmin(lamb) and \
            np.nanmax(lamb_r) <= np.nanmax(lamb):
        
        dens=np.array(lamb)*np.array(flambda)/\
            (phc.h.cgs*phc.c.cgs)   ### lamb*flamb/hc [cgs units]
        
        ### domain of lambdas [Angstroms] for interpolation
        l=np.array([np.nanmin(lamb_r)+(np.nanmax(lamb_r)-\
            np.nanmin(lamb_r))/float(npts-1)*float(i) \
            for i in xrange(0,npts)])
        ### variation of lambda [Angstroms] for integration
        dl=np.array([l[i+1]-l[i] for i in xrange(0,npts-1)])
        
        ### interpolated lamb*flamb/hc [cgs units]
        d=np.array([interpLinND([l[i]],[lamb],dens) \
            for i in xrange(0,npts)])
        ### Interpolated passband
        r=np.array([interpLinND([l[i]],[lamb_r],R) \
            for i in xrange(0,npts)])
        ### photon flux [s^-1 cm^-2]
        photflux=integrate_trapezia(d*r,dl*10.**-8.)
        return photflux
        
    else: return np.nan



def pogson(X,zp):
    return -2.5*np.log10(X)+zp



def VEGA_spct(spct_name):
    """
    Returns the SED of VEGA:
    * array of lambdas [Angstroms]
    * array of flambda [erg cm^-2 s^-1 A^-1]
    """

    if spct_name == 'spct1': ### Reference: 1994A&A...281..817C
        ### Reading Vega's spectrum
        spct1=cwd+'../refs/stars/fm05t9550g395k2odfnew.dat'
        f0 = open(spct1,'r')
        lines = f0.readlines()
        f0.close()

        ### Eliminating text and separating columns
        lines=[lines[i].split() for i in xrange(12,len(lines))]

        dist_fac=1./(1.62*10.**16.) ### Reference: 1994A&A...281..817C
        ### lambda [angstroms]
        lamb=np.array([float(lines[i][2])*10. for i in xrange(0,len(lines))])
        ### flambda [erg cm^-2 s^-1 A^-1]
        flambda=np.array([10.**(np.log10(4.*np.pi)+np.log10(phc.c.cgs)+\
                            np.log10(float(lines[i][4]))\
                            -2.*np.log10(float(lines[i][2])*10.**-7.)) \
                            for i in xrange(0,len(lines))])/10.**8.*\
                            dist_fac
        return lamb,flambda



def obtain_pogson_zp(spct_name,filtername,npts=50):
    """
    Returns the zp for a pogson magnitude system, whose filter is
    'filtername'.
    The 'spct_name' contains the name of the spectrum (usually of Vega)
    to be used as calibrator.
    """


    #lamb,flambda=VEGA_spct(spct_name)
    if filtername == 'STMAG':
        #flb=interpLinND([5480.],[lamb],flambda)
        #k=photosystem_k(filtername)
        #zp=2.5*np.log10(flb)+k
        zp = -21.1
        return zp
    else:
        lamb,flambda=VEGA_spct(spct_name)
        photflux=photonflux(lamb,flambda,filtername,npts)
        k=photosystem_k(filtername)
        zp=2.5*np.log10(photflux)+k
        return zp





#####

def fullsed2photonflux(fullsed,source,filtername,npts=50,dist=10.):
    """
    Reads a HDUST's fullsed and source files + 1 'filtername'
    and returns lists of 
    * mu (cosi), 
    * arrays of lambda [in Angs](for each cosi), 
    * arrays of flambda [erg/cm^2 s A](for each cosi), 
    * photon fluxes [1/cm^2 s](for each cosi) in the filter 'filtername'
    """
    ### only the total flux (dont care about polarization)
    ### I am not taking phi into account. 
    ### Currently, this function is only usable for axisymmetric 
    ### simulations.

    ### Reading the fullsed file
    f0 = open(fullsed,'r')
    lines = f0.readlines()
    f0.close()
    nlbd=int(lines[3].split()[0]); nobs=int(lines[3].split()[1])
    lines=[lines[i].split() for i in xrange(5,len(lines))]
    
    ### Reading the source file
    f0 = open(source,'r')
    src = f0.readlines()
    f0.close()
    lum=float(src[6].split()[2])
    normflx=(lum*phc.Lsun.cgs)/(4.*np.pi*(dist*phc.pc.cgs)**2.)

    ### Calculating the photon fluxes for each inclination
    mu=[];lamb=[];flambda=[];photflux=[]
    for iob in xrange(0,nobs):
        mu_aux=float(lines[0+nlbd*iob][0]) ### cosi
        lamb_aux=np.array([float(lines[ilamb+nlbd*iob][2])*10000. \
            for ilamb in xrange(0,nlbd)]) ### lambda [Angstroms]
        flambda_aux=np.array([float(lines[ilamb+nlbd*iob][3])*\
            normflx/10000. \
            for ilamb in xrange(0,nlbd)]) ### flambda [erg/cm^2 s A]
        photflux_aux=photonflux(lamb_aux,flambda_aux,filtername,npts)
        ### Appending results for this specific inclination
        mu.append(mu_aux)
        lamb.append(lamb_aux)
        flambda.append(flambda_aux)
        photflux.append(photflux_aux)
    return mu,lamb,flambda,photflux



    








