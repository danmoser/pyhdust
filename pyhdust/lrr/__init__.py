# -*- coding:utf-8 -*-

import sys
import os
import numpy as np
import itertools as it
import pyhdust.phc as phc
from scipy.optimize import curve_fit

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
        print('<<',func_name,'>>')
        print('There is something wrong with the integration!')
        print('')
        return np.nan
    else: 
        soma = 0.0
        for i in range(0,len(f) - 1):
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
        print('<<',func_name,'>>')
        print('The number must be a natural number.')
        print('')
        return np.nan

    addreverse=[]; res=number%2; quo=int(number/2)
    addreverse.append(res)
    while quo != 0: 
        res=quo%2; quo=int(quo/2); addreverse.append(res)
    
    if dim != None and dim < len(addreverse):
        func_name = sys._getframe().f_code.co_name
        print('<<',func_name,'>>')
        print('I need more dimensions to write this binary number!')
        print('')
        return np.nan
    elif dim != None and dim >= len(addreverse): 
        addn=np.array([0 for k in range(0,dim)])
    else: 
        addn=np.array([0 for k in range(0,len(addreverse))])
    
    for j in range(0,len(addreverse)): 
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
                print("Oops!  That was no valid dominium.  Try again...")
        x=np.random.randint(between[0],between[1]+1)
        print(" ==> Here is your number: "+str(x))
        return x

    else:
        while True:
            try:
                between[0] = int(minx)
                between[1] = int(maxx)
                teste = np.sqrt(float(between[1])-float(between[0]))
                break
            except ValueError:
                print("Oops!  That was no valid dominium.  Try again...")
        x=np.random.randint(between[0],between[1]+1)
        if p == 'yes': print(" ==> Here is your number: "+str(x))
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

def scale_two_propto(x,A1,A2,m="normal"):

    if m != "inverse":
        if x >= 0: return A1*x
        if x < 0: return A2*x
    else:
        if x >= 0: return x/A1
        if x < 0: return x/A2

def scale_two_arcsinh(x,A1,B1,A2,B2,m="normal"):
    """
    Transforms a real number x into the real number y given by
    y = A*arcsinh(B*x) = A*ln[B*x + ((B*x)^2 + 1)^0.5]
    Also, the derivative of y is 
    y' = A*B*((B*x)^2 + 1)^(-0.5)
    
    Good guidelines for choosing the quantity B are the following properties:
    . If |B*x| >> 1: y ~ sign(B*x)*A*ln[2*|B*x|] and y' ~ A*sign(B)*|x|^(-1)
    . If |B*x| << 1: y ~ A*B*x and y' ~ A*B

    The inverse of the transformation is given by 
    x = 1/B * sinh(y/A) = 1/B * 0.5*[e^(y/A)-e^(-y/A)]
    """
    if np.isnan(x):
        return np.nan
    
    if m == "normal":
        if x >= 0: return A1*np.arcsinh(x*B1)
        if x < 0: return A2*np.arcsinh(x*B2)
    elif m == "inverse":
        if x >= 0: return 1./B1*np.sinh(x/A1)
        if x < 0: return 1./B2*np.sinh(x/A2)
    elif m == "deriv":
        if x >= 0: return A1*B1*((B1*x)*(B1*x)+1.)**(-0.5)
        if x < 0: return A2*B2*((B2*x)*(B2*x)+1.)**(-0.5)
    else:
        return np.nan

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
    for i in range(0,len(axis)-1):
        prod=1
        for j in range(i+1,len(axis)):
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
    for i in range(0,N_Fx):
        addn=dec_2_binary(i,len(ind))
        #print("addn = ",addn)
        newind=np.array([ind[j]+addn[j] for j in range(len(ind))])
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
        print('<<',func_name,'>>')
        print('X and axis must have the same dimension!')
        print('')
        return np.nan
    low,high,ind,extrapolated=low_high(X,axis)
    if allow_extrapolation != "yes" and extrapolated == "yes":
        return np.nan
    else:
        Fx=build_Fx(axis,values,ind)
        interp=interLinND(X, low, high, Fx, tp)
        return interp






def indexes_of_arrows(axis,index_edge00):
    """
    This is an auxiliary function of the function 'interpLinNDpowerful'.
    
    With the N-dimensional 'axis', containing the domain of the 
    function to be interpolated, 
    and with the indexes of an N-dimensional point in this domain, 
    it creates ALL possible "forward arrows" for that point.
    
    The point and the arrows are like this:
                      |
      o ------->      |
      |               |
      |               |
      v               |
                      | * limits of hyperrectangular domain
                      |
    __________________| 

    and the idea is that, with N normal arrows, we define an N-dimensional
    hyperrectangle.
    
    In the above drawing, the specific 'indx_arrows' will be given by (7,2).
    """
    ### 
    auxi_indx_arrows = []
    ### 
    for idim in range(0,len(axis)):
        auxi = []
        j0 = int(np.nanmax([0,index_edge00[idim]+1]))
        index = j0
        for j in range(j0,len(axis[idim])):
            auxi.append(index-index_edge00[idim])
            index+=1
        auxi_indx_arrows.append(auxi)
    list_indx_arrows = list(it.product(*auxi_indx_arrows))
    return list_indx_arrows
    

def axis_for_hyperrectangle(index_edge00,list_indx_arrows,axis):
    """
    This is an auxiliary function of the function 'interpLinNDpowerful'.
    
    This function generates a list for every hyperrectangle defined by 
    'list_indx_arrows' (generated by function 'indexes_of_arrows'). 
    The elements of this list are lists N 2-arrays 
    (where N is the dimension of the problem) which are 
    the axis that generate each hyperrectangle.
    """
    ### 
    list_new_axis = []
    list_new_indexes = []
    ### 
    for i_ia in range(0,len(list_indx_arrows)):
        auxi_new_axis = []
        auxi_new_indexes = []
        ### Loop over the dimensions of the problem.
        for idim in range(0,len(axis)):
            auxi_new_indexes.append((
                index_edge00[idim],\
                index_edge00[idim]+list_indx_arrows[i_ia][idim]\
                ))
                
            auxi_new_axis.append((
                axis[idim][index_edge00[idim]],\
                axis[idim][index_edge00[idim]+list_indx_arrows[i_ia][idim]]\
                ))
        list_new_indexes.append(auxi_new_indexes)
        list_new_axis.append(auxi_new_axis)
        
    return list_new_axis, list_new_indexes

def find_extrap(point,list_new_axis):
    foundextrap = []
    for i_ia in range(0,len(list_new_axis)):
        foundextrap.append(0)
        for idim in range(0,len(list_new_axis[i_ia])):
            if not (list_new_axis[i_ia][idim][0] <= point[idim] <= \
                    list_new_axis[i_ia][idim][1]):
                foundextrap[-1] = 1
    return foundextrap    

    
    
def get_hypervolumes(list_new_axis):
    
    
    hypervolumes = []
    for i_ia in range(0,len(list_new_axis)):
        vol = 1.
        for idim in range(0,len(list_new_axis[i_ia])):
            vol = vol*(list_new_axis[i_ia][idim][1]-list_new_axis[i_ia][idim][0])
        hypervolumes.append(vol)
    
    return hypervolumes
       






def get_list_values(list_new_axis,list_new_indexes,axis,values):
    """
    This is an auxiliary function of the function 'interpLinNDpowerful'.
    
    This function generates the list of 2**N values associated with the 
    elements of 'list_new_axis'. 
    """
    
    ### 'list_vertexes_new_axis' will contain a list composed of the 2**N 
    ### vertexes for each N-dimensional element of 'list_new_axis'.
    list_vertexes_new_axis = []
    for i_ia in range(0,len(list_new_axis)):
        list_vertexes_new_axis.append(\
                list(it.product(*list_new_axis[i_ia]))\
                )
    
    ### 'list_indvertexes_new_axis' will contain a list composed of the 2**N 
    ### indexes of vertexes for each N-dimensional element of 'list_new_axis'.
    list_indvertexes_new_axis = []
    for i_ia in range(0,len(list_new_axis)):
        list_indvertexes_new_axis.append(\
                list(it.product(*list_new_indexes[i_ia]))\
                )
    
    list_values = []
    ### 'listao_axis' will contain every point of the domain defined by
    ### 'axis'.
    listao_axis = list(it.product(*axis))
    ### Loop over every collection of vertexes associated with an element of 
    ### 'list_new_axis'. 
    leng = len(axis)
    leng_vec = [len(x) for x in axis]
    for i_ia in range(0,len(list_indvertexes_new_axis)):
        auxi_vals = []
        ### Loop over every vertex. 
        ### The index of the correspondent element of 'listao_axis'
        ### will be collected. It will be used to select the elements of 
        ### 'values' and compose the list 'list_values'
        for ivert in range(0,len(list_indvertexes_new_axis[i_ia])):
            indx_vals = 0.
            ### Loop over the dimensions of the problem:
            ### (This is the heaviest part of this function.)
            for idim in range(0,leng):
                prod = 1.
                for iax in range(1,idim+1):
                    prod = prod*leng_vec[leng-iax]
                indx_vals += list_indvertexes_new_axis[i_ia][ivert][leng-1-idim]*prod
            auxi_vals.append(values[int(indx_vals)])
        list_values.append(np.array(auxi_vals))
    
    
    return list_values



def findNaN(list_values):
    foundNaN = []
    for i_ia in range(0,len(list_values)):
        foundNaN.append(0)
        #print(list_values[i_ia])
        for elem in list_values[i_ia]:
            if np.isnan(elem):
                foundNaN[-1] = 1
    return foundNaN


def eliminate_foundNaNs(lista,foundNaN):
        
    new_lista = []
    if len(lista) == len(foundNaN):
        for i in range(0,len(lista)):
            if foundNaN[i] == 0:
                new_lista.append(lista[i])
    else:
        new_lista = [x for x in lista]
        
    return new_lista
    

def hyperrectangles_for_interp(axis, values, prints="no"):
    """
    
    """
    
    ### Listing the indexes of all the "00-vertexes".
    ### From each of these 00-edges, the "forward hyperrectangles"
    ### will be built.
    list_index_edge00 = []
    xs = []
    for iaxis in range(0,len(axis)):
        xs.append([0+ix for ix in range(0,len(axis[iaxis]))])
    list_index_edge00 = list(it.product(*xs))
    
    ### 
    llist_new_axis = []
    llist_new_indexes = []
    llist_values = []
    if prints != "no":
        import sys
    for ilie in range(0,len(list_index_edge00)):
        ### One specific edge's index
        index_edge00 = list_index_edge00[ilie]
        ### for this specific 'index_edge00',
        ### obtaining the indexes of the arrows that build the 
        ### "forward hyperrectangles".
        list_indx_arrows = indexes_of_arrows(axis,index_edge00)
        ### For each "forward hyperrectangle" (of dimension N), 
        ### find the domain of N 2-arrays that define the hyperrectangle.
        auxi_list_new_axis, auxi_list_new_indexes = \
            axis_for_hyperrectangle(index_edge00,list_indx_arrows,axis)
        llist_new_axis.append(auxi_list_new_axis)
        llist_new_indexes.append(auxi_list_new_indexes)
        ### For each "forward hyperrectangle" (of dimension N), 
        ### find the 2**N array of values in each of the vertexes.
        ### (This is the heaviest part of this function!!)
        auxi_list_values = get_list_values(auxi_list_new_axis,\
            auxi_list_new_indexes,axis,values)
        llist_values.append(auxi_list_values)
        if prints != "no":
            sys.stdout.write("\rFinding all hyperrectangles: {:2.0%}".\
                    format(float(ilie)/\
                    (len(list_index_edge00)-1))+"     ")
            sys.stdout.flush()
    if prints != "no":
        sys.stdout.write("\rFinding all hyperrectangles: DONE \n")

    ### Joining the lists within each of the two lists from the 
    ### above loop in one big list (for each of the two lists).
    list_new_axis = []
    list_values = []
    for ilie in range(0,len(list_index_edge00)):
        for ielem in range(0,len(llist_new_axis[ilie])):
            list_new_axis.append(llist_new_axis[ilie][ielem])
        for ielem in range(0,len(llist_values[ilie])):
            list_values.append(llist_values[ilie][ielem])

    ### Now, comes the selection of the best hyperrectangle to perform
    ### the interpolation. It will be the one that satisfies a certain
    ### criterion of distance of the hyperrectangle from the point, 
    ### volume of the hyperrectangle, and if extrapolation is allowed.

    ### For every hyperrectangle, check if there is at least a NaN
    ### in the values at the vertexes
    if prints != "no":
        sys.stdout.write("Finding NaNs        ")
        sys.stdout.flush()
    foundNaN = findNaN(list_values)
    if prints != "no":
        sys.stdout.write("\rFinding NaNs: DONE \n")
    
    return list_new_axis, list_values, foundNaN
    


    
    


def get_quaddistances(X,list_new_axis):
        
    ### Obtaining the coordinates of the centers of the hyperrectangles
    centers = []
    for ilist in range(0,len(list_new_axis)):
        auxi_centers = []
        for idim in range(0,len(list_new_axis[ilist])):
            mean = 0.5*(list_new_axis[ilist][idim][0]+\
                        list_new_axis[ilist][idim][1])
            auxi_centers.append(mean)
        centers.append(auxi_centers)
        
    ### Obtaining quadratic distances of every center of hyperrectangle
    ### and the point 'X'
    quaddistances = []
    for ilist in range(0,len(centers)):
        sum_qd = 0.
        for idim in range(0,len(centers[ilist])):
            sum_qd += (centers[ilist][idim]-X[idim])*\
                            (centers[ilist][idim]-X[idim])
        quaddistances.append(sum_qd)
        
    return quaddistances    




def selecting_by_merit(X,list_new_axis,tominimize,\
            list_values,allow_extrapolation):   
    """
    
    """




    ### The criterion is: the best hyperrectangle is the one that 
    ### minimizes its volume multiplied by the sum of quadratic distances 
    ### of its vertexes to the point.
    lower_hyp_d2 = np.inf
    indx_min = np.nan
    if allow_extrapolation != "yes":
        for i in range(0,len(list_new_axis)):
            if tominimize[i] < lower_hyp_d2:
                ### For every hyperrectangle, check if the point 'X' is inside or 
                ### outside the hyperrectangle. (In the letter case, extrapolation is 
                ### required.)
                foundextrap = find_extrap(X,[list_new_axis[i]])
                #print(tominimize[i],foundextrap)
                if not (1 in foundextrap):
                    lower_hyp_d2 = tominimize[i]
                    indx_min = i
    else:
        for i in range(0,len(list_new_axis)):
            if tominimize[i] < lower_hyp_d2:
                lower_hyp_d2 = tominimize[i]
                indx_min = i                    
    
    ### Finding the list of N 2-arrays for the vertexes of the best 
    ### hyperrectangle, and
    ### the 2**N array of values in each of the vertexes.
    if ~np.isnan(indx_min):
        the_new_axis = [x for x in list_new_axis[indx_min]]
        the_new_values = np.array([x for x in list_values[indx_min]])
    else:
        the_new_axis = []
        the_new_values = []
        
    return the_new_axis, the_new_values

    




    

def interpLinNDpowerful(X,axis,values,tp="linear",allow_extrapolation="yes"):
    """
    N-dimensional linear interpolation (in the same way as 'interpLinND').
    
    This function, however, is to be used when there are "holes" in the 
    'values' (containing NaNs). 
    
    This function searches for the "closest" hyperrectangle in the domain, 
    whose vertexes all contain values different of NaN. Then it performs 
    the interpolation using this hyperrectangle.
    """
    
    ### 
    list_new_axis, list_values, foundNaN = \
        hyperrectangles_for_interp(axis, values, prints="no")
    ### 
    hypervolumes = get_hypervolumes(list_new_axis)
    
    
    e_list_new_axis = eliminate_foundNaNs(list_new_axis,foundNaN)
    e_list_values = eliminate_foundNaNs(list_values,foundNaN)
    e_hypervolumes = eliminate_foundNaNs(hypervolumes,foundNaN)
    
    e_quaddistances = get_quaddistances(X,e_list_new_axis)
    e_hyp_d2 = [e_quaddistances[i]*e_hypervolumes[i] \
        for i in range(0,len(e_hypervolumes))]    
    
    
    
    the_new_axis, the_new_values = \
                    selecting_by_merit(X,\
                    e_list_new_axis,e_hyp_d2,e_list_values,\
                    allow_extrapolation)
    
    
    ### Finally, performing the interpolation:
    if len(the_new_axis) > 0 and len(the_new_values) > 0:
        interp = interpLinND(X,the_new_axis,the_new_values,\
                    tp,allow_extrapolation)
    else:
        interp = np.nan
    
    return interp





def fill_NaNs_interp(axis,values,tp="linear",allow_extrapolation="yes",prints="no"):
    """
    This function will fill the NaNs in 'values' with the powerful
    interpolation. 
    
    
    """

    if prints != "no":
        import sys
    
    
    ### 
    list_new_axis, list_values, foundNaN = \
            hyperrectangles_for_interp(axis, values, prints)
    ### 
    hypervolumes = get_hypervolumes(list_new_axis)
    ### 
    allpoints = list(it.product(*axis))
    ### 
    not_foundNaN = 0 in foundNaN

    ### TODO: make new 'list_new_axis', 'list_values', 'hypervolumes', 
    ### which are clean of all occurencies of 'foundNaN[i]'=1.
    ### They will enter into 'selecting_by_merit_v2'.
    ### DONE
    e_list_new_axis = eliminate_foundNaNs(list_new_axis,foundNaN)
    e_list_values = eliminate_foundNaNs(list_values,foundNaN)
    e_hypervolumes = eliminate_foundNaNs(hypervolumes,foundNaN)
    
    
    new_values = []
    for ivalue in range(0,len(allpoints)):
        if np.isnan(values[ivalue]) and not_foundNaN:
            ### TODO: here, call function to calculate quadratic distances
            ### to 'allpoints[ivalue]'.
            ### DONE
            e_quaddistances = get_quaddistances(allpoints[ivalue],\
                    e_list_new_axis)
            e_hyp_d2 = [e_quaddistances[i]*e_hypervolumes[i] \
                for i in range(0,len(e_hypervolumes))]
            ### 
            the_new_axis, the_new_values = \
                    selecting_by_merit(allpoints[ivalue],\
                    e_list_new_axis,e_hyp_d2,e_list_values,\
                    allow_extrapolation)
            if len(the_new_axis) > 0 and len(the_new_values) > 0:
                newval = interpLinND(allpoints[ivalue],the_new_axis,\
                        the_new_values,tp,allow_extrapolation)
            else:
                newval = np.nan
    
        else:
            newval = values[ivalue]
        new_values.append(newval)
        if prints != "no":
            sys.stdout.write("\rTrying to fill NaNs: {:2.3%}".\
                    format(float(ivalue)/(len(allpoints)-1))+"     ")
            sys.stdout.flush()
    if prints != "no":
        sys.stdout.write("\rTrying to fill NaNs: DONE      \n")
    
    
    return new_values







######################################################################
### Photometric routines


def photosystem_k(filtername):
    """
    Returns the magnitude of Vega in the photometric system of the 
    filter 'filtername'.
    
    * Non-standard filters receive 9999 as the magnitude of Vega.
    """
    if filtername == 'STMAG': return 0.00
    elif filtername == 'white': return 0.00
    ### Reference: 1998A&A...333..231B
    elif filtername == 'bess-u': return 0.03 
    elif filtername == 'bess-b': return 0.03 
    elif filtername == 'bess-v': return 0.03 
    elif filtername == 'bess-r': return 0.03 
    elif filtername == 'bess-i': return 0.03 
    elif filtername == 'bess-j': return 0.03 
    elif filtername == 'bess-h': return 0.03 
    elif filtername == 'bess-k': return 0.03 
    elif filtername == 'bess-l': return 0.03 
    elif filtername == 'bess-ll': return 0.03 
    elif filtername == 'bess-m': return 0.03 
    ### Reference: 2003AJ....126.1090C
    elif filtername == 'cohen-j': return 0.00 
    elif filtername == 'cohen-h': return 0.00 
    elif filtername == 'cohen-k': return 0.00 
    ### Reference: 1998AJ....116..482G
    elif filtername == 'crawford-stromgren-u': return 1.445 
    elif filtername == 'crawford-stromgren-v': return 0.195 
    elif filtername == 'crawford-stromgren-b': return 0.034 
    elif filtername == 'crawford-stromgren-y': return 0.030
    ### 
    elif filtername == 'int_wfc-hbetan': return 9999.
    elif filtername == 'int_wfc-hbetaw': return 9999.
    elif filtername == 'int_wfc-halpha': return 9999.
    ### Reference: http://www.ctio.noao.edu/soar/content/filters-available-soar
    elif filtername == 'filt_Ha': return 9999.
    ### Reference: 2008PASP..120.1233H
    elif filtername == 'irac_ch1': return 0.00
    elif filtername == 'irac_ch2': return 0.00
    elif filtername == 'irac_ch3': return 0.00
    elif filtername == 'irac_ch4': return 0.00
    ### Reference: 2011ApJ...735..112J
    elif filtername == 'wise1': return 0.00
    elif filtername == 'wise2': return 0.00
    elif filtername == 'wise3': return 0.00
    elif filtername == 'wise4': return 0.00
    else: 
        func_name = sys._getframe().f_code.co_name
        print('<<',func_name,'>>')
        print(filtername,' ?? I do not know which filter is that!')
        print('')
        return np.nan


def filter_auxi(filename,lambfac,Rfac,lininit): 

    f0 = open(filename,'r')
    lines = f0.readlines()
    f0.close()
    ### Eliminating text and separating columns
    lines=[lines[i].split() for i in range(lininit,len(lines))]
    ### lambda (Angstroms)
    lamb=np.array([float(lines[i][0])*lambfac for i in range(0,len(lines))])
    ### Reponse function
    R=np.array([float(lines[i][1])*Rfac for i in range(0,len(lines))])
    return lamb,R


def filter_passband(filtername):
    """
    Returns, for the filter 'filtername'
    * array of lambdas [Angstroms]
    * array of the function R_X (I don't know if it is normalized!)
    """
    
    if filtername == 'STMAG':
        return np.nan,np.nan
    elif filtername == 'white':
        return np.array([0.,1e300]),np.array([1.,1.])
    elif filtername == 'bess-u': ### Reference: 1990PASP..102.1181B (Table 2)
        lamb, R = filter_auxi(cwd+'../refs/filters/Bessell_U.dat',1.,1.,1)
        return lamb,R
    elif filtername == 'bess-b': ### Reference: 1990PASP..102.1181B (Table 2)
        lamb, R = filter_auxi(cwd+'../refs/filters/Bessell_B.dat',1.,1.,1)
        return lamb,R
    elif filtername == 'bess-v': ### Reference: 1990PASP..102.1181B (Table 2)
        lamb, R = filter_auxi(cwd+'../refs/filters/Bessell_V.dat',1.,1.,1)
        return lamb,R
    elif filtername == 'bess-r': ### Reference: 1990PASP..102.1181B (Table 2)
        lamb, R = filter_auxi(cwd+'../refs/filters/Bessell_R.dat',1.,1.,1)
        return lamb,R
    elif filtername == 'bess-i': ### Reference: 1990PASP..102.1181B (Table 2)
        lamb, R = filter_auxi(cwd+'../refs/filters/Bessell_I.dat',1.,1.,1)
        return lamb,R
    elif filtername == 'bess-j': ### Reference: 1988PASP..100.1134B (Table IV)
        lamb, R = filter_auxi(cwd+'../refs/filters/Bessell_J.dat',1.,1.,1)
        return lamb,R
    elif filtername == 'bess-h': ### Reference: 1988PASP..100.1134B (Table IV)
        lamb, R = filter_auxi(cwd+'../refs/filters/Bessell_H.dat',1.,1.,1)
        return lamb,R
    elif filtername == 'bess-k': ### Reference: 1988PASP..100.1134B (Table IV)
        lamb, R = filter_auxi(cwd+'../refs/filters/Bessell_K.dat',1.,1.,1)
        return lamb,R
    elif filtername == 'bess-l': ### Reference: 1988PASP..100.1134B (Table IV)
        lamb, R = filter_auxi(cwd+'../refs/filters/Bessell_L.dat',1.,1.,1)
        return lamb,R
    elif filtername == 'bess-ll': ### Reference: 1988PASP..100.1134B (Table IV)
        lamb, R = filter_auxi(cwd+'../refs/filters/Bessell_LL.dat',1.,1.,1)
        return lamb,R
    elif filtername == 'bess-m': ### Reference: 1988PASP..100.1134B (Table IV)
        lamb, R = filter_auxi(cwd+'../refs/filters/Bessell_M.dat',1.,1.,1)
        return lamb,R
    elif filtername == 'cohen-j': ### Reference: 2003AJ....126.1090C
        lamb, R = filter_auxi(cwd+'../refs/filters/Cohen-J.pass',1.,1.,1)
        return lamb,R
    elif filtername == 'cohen-h': ### Reference: 2003AJ....126.1090C
        lamb, R = filter_auxi(cwd+'../refs/filters/Cohen-H.pass',1.,1.,1)
        return lamb,R
    elif filtername == 'cohen-k': ### Reference: 2003AJ....126.1090C
        lamb, R = filter_auxi(cwd+'../refs/filters/Cohen-K.pass',1.,1.,1)
        return lamb,R
    elif filtername == 'crawford-stromgren-u': ### Reference: 1970AJ.....75..978C
        lamb, R = filter_auxi(cwd+'../refs/filters/Crawford-Stromgren-U.dat',1.,1.,2)
        return lamb,R
    elif filtername == 'crawford-stromgren-v': ### Reference: 1970AJ.....75..978C
        lamb, R = filter_auxi(cwd+'../refs/filters/Crawford-Stromgren-V.dat',1.,1.,2)
        return lamb,R
    elif filtername == 'crawford-stromgren-b': ### Reference: 1970AJ.....75..978C
        lamb, R = filter_auxi(cwd+'../refs/filters/Crawford-Stromgren-B.dat',1.,1.,2)
        return lamb,R
    elif filtername == 'crawford-stromgren-y': ### Reference: 1970AJ.....75..978C
        lamb, R = filter_auxi(cwd+'../refs/filters/Crawford-Stromgren-Y.dat',1.,1.,2)
        return lamb,R
    elif filtername == 'int_wfc-hbetan':
        lamb, R = filter_auxi(cwd+'../refs/filters/INT_WFC-Hbetan.dat',1.,1.,0)
        return lamb,R
    elif filtername == 'int_wfc-hbetaw':
        lamb, R = filter_auxi(cwd+'../refs/filters/INT_WFC-Hbetaw.dat',1.,1.,0)
        return lamb,R
    elif filtername == 'int_wfc-halpha':
        lamb, R = filter_auxi(cwd+'../refs/filters/INT_WFC-Halpha.dat',1.,1.,0)
        return lamb,R
    elif filtername == 'filt_Ha': ### Reference: http://www.ctio.noao.edu/soar/content/filters-available-soar
        lamb, R = filter_auxi(cwd+'../refs/filters/filt_Ha.dat',10.,0.01,5)
        return lamb,R
    elif filtername == 'irac_ch1': ### Reference: 2008PASP..120.1233H
        lamb, R = filter_auxi(cwd+'../refs/filters/irac_ch1.txt',1e4,1.,3)
        return lamb,R
    elif filtername == 'irac_ch2': ### Reference: 2008PASP..120.1233H
        lamb, R = filter_auxi(cwd+'../refs/filters/irac_ch2.txt',1e4,1.,3)
        return lamb,R
    elif filtername == 'irac_ch3': ### Reference: 2008PASP..120.1233H
        lamb, R = filter_auxi(cwd+'../refs/filters/irac_ch3.txt',1e4,1.,3)
        return lamb,R
    elif filtername == 'irac_ch4': ### Reference: 2008PASP..120.1233H
        lamb, R = filter_auxi(cwd+'../refs/filters/irac_ch4.txt',1e4,1.,3)
        return lamb,R
    elif filtername == 'wise1': ### Reference: 2011ApJ...735..112J
        lamb, R = filter_auxi(cwd+'../refs/filters/RSR-W1.EE.txt',1e4,1.,2)
        return lamb,R
    elif filtername == 'wise2': ### Reference: 2011ApJ...735..112J
        lamb, R = filter_auxi(cwd+'../refs/filters/RSR-W2.EE.txt',1e4,1.,2)
        return lamb,R
    elif filtername == 'wise3': ### Reference: 2011ApJ...735..112J
        lamb, R = filter_auxi(cwd+'../refs/filters/RSR-W3.EE.txt',1e4,1.,2)
        return lamb,R
    elif filtername == 'wise4': ### Reference: 2011ApJ...735..112J
        lamb, R = filter_auxi(cwd+'../refs/filters/RSR-W4.EE.txt',1e4,1.,2)
        return lamb,R

    else: 
        func_name = sys._getframe().f_code.co_name
        print('<<',func_name,'>>')
        print(filtername,' ?? I do not know which filter is that!')
        print('')
        return np.nan,np.nan







def mean_RX(filtername,m_wavelength="no"):
    """
    Returns [Angstroms]:
    * bandwidth of filter 'filtername', if 'bandwidth' == "yes"
    * mean wavelength and bandwidth of filter 'filtername', if 'bandwidth' != "yes"
    """

    ### Obtaining the passband associated with the filter 'filtername'
    lamb_r,R=filter_passband(filtername)
        
    ### variation of lambda [Angstroms] for integration
    dl=np.array([lamb_r[i+1]-lamb_r[i] for i in range(0,len(lamb_r)-1)])
    ### Obtaining bandwidth [Angstroms]
    bdwidth = integrate_trapezia(R,dl)

    if m_wavelength == "no":
        return bdwidth
    else:
        mean = integrate_trapezia(R*lamb_r,dl)/bdwidth  ### mean wavelength [Angs]
        return mean, bdwidth






def color_from_alpha(alphavec,filter1,filter2,zp1,zp2):
    """
    |OUTPUT: 
    | * a vector of colors associated with a vector of alphas (spectral indexes).
    | * a vector of multiplicative factors, to multiply on the error of color and 
    | get the error on alpha.
    """
    
    ### Obtaining the passbands associated with the filter 'filter1' and 'filter2'
    lamb_r1,R1=filter_passband(filter1)
    lamb_r2,R2=filter_passband(filter2)

    ### variation of lambda [Angstroms] for integration
    dl1=np.array([lamb_r1[i+1]-lamb_r1[i] for i in range(0,len(lamb_r1)-1)])
    dl2=np.array([lamb_r2[i+1]-lamb_r2[i] for i in range(0,len(lamb_r2)-1)])
    
    ### 
    color = []
    err_factor = []
    for alpha in alphavec:
        integral1 = integrate_trapezia(R1*lamb_r1**alpha,dl1)
        integral_ln1 = integrate_trapezia(R1*lamb_r1**alpha*np.log(lamb_r1),dl1)
        integral2 = integrate_trapezia(R2*lamb_r2**alpha,dl2)
        integral_ln2 = integrate_trapezia(R2*lamb_r2**alpha*np.log(lamb_r2),dl2)
        color.append(-2.5*np.log10(integral1/integral2)+zp1-zp2)
        err_factor.append(\
            1./abs(-2.5/np.log(10.)*(integral_ln1/integral1-integral_ln2/integral2))\
            )
    color = np.array([elem for elem in color])
    err_factor = np.array([elem for elem in err_factor])
    
    return color, err_factor
    
    



def photonflux(lamb,flambda,filtername,npts=50,forcedred = 0,pf = "yes"):
    """
    | INPUT: 
    | * 'lamb': array of lambda [Angstroms]
    | * 'flambda': flux density [erg cm^-2 s^-1 A^-1]
    | * 'filtername'
    | * 'forcedred': number of last points used in the extrapolation of 
    | 'flambda', if the filter extends to redder wavelengths.
    | * 'pf': * if "yes", then calculates the photon flux
    |         * if not "yes", then calculates the mean flux density

    | OUTPUT: 
    | * If 'pf' == "yes": photon flux [photons/cm2 s]
    | * If 'pf' != "yes": mean flux density*bandwidth [erg cm^-2 s^-1]
    """
    
    
    ### Obtaining the passband associated with the filter 'filtername'
    lamb_r,R=filter_passband(filtername)
    
    if filtername == 'white':
        
        if pf == "yes":
            dens=np.array(lamb)*np.array(flambda)/\
                (phc.h.cgs*phc.c.cgs)   ### lamb*flamb/hc [cgs units]
        else:
            dens=np.array(flambda)*1e8      ### flamb [erg cm^-2 s^-1 cm^-1]
        
        ### domain of lambdas [Angstroms] for interpolation
        l=np.array([np.nanmin(lamb)+(np.nanmax(lamb)-np.nanmin(lamb))/\
            float(npts-1)*float(i) for i in range(0,npts)])
        
        ### variation of lambda [Angstroms] for integration
        dl=np.array([l[i+1]-l[i] for i in range(0,npts-1)])
        
        ### interpolated lamb*flamb/hc [cgs units]
        d=np.array([interpLinND([l[i]],[lamb],dens) \
            for i in range(0,npts)])
        ### photon flux [s^-1 cm^-2] or 
        ### mean flux density*bandwidth [erg cm^-2 s^-1]
        photflux=integrate_trapezia(d,dl*1e-8)
        return photflux
    
    elif np.nanmin(lamb_r) >= np.nanmin(lamb) and \
            (np.nanmax(lamb_r) <= np.nanmax(lamb) or forcedred >= 2):
        
        ### Extending 'lamb' and 'flambda', if forcedred >= 2:
        if np.nanmax(lamb_r) > np.nanmax(lamb):
            def powerlaw(x,C,alpha):
                """
                Assuming a power-law for the fitting of a portion of the SED.
                """
                return C*x**alpha
            
            ### Collecting the last elements of the SED
            auxi_lamb = []
            auxi_flambda = []
            for ii in range(0,forcedred):
                auxi_lamb.append(lamb[len(lamb)-1-(forcedred-1-ii)])
                auxi_flambda.append(flambda[len(flambda)-1-(forcedred-1-ii)])
            auxi_lamb = np.array([elem for elem in auxi_lamb])
            auxi_flambda = np.array([elem for elem in auxi_flambda])
            ### Fitting power-law to last portion of the SED
            auxi_alpha0 = []
            for ii in range(0,forcedred-1):
                auxi_alpha0.append(np.log(auxi_flambda[-1-ii]/auxi_flambda[-2-ii])/\
                            np.log(auxi_lamb[-1-ii]/auxi_lamb[-2-ii]))
            alpha0 = np.mean(np.array(auxi_alpha0))
            C0 = auxi_flambda[-1]/auxi_lamb[-1]**alpha0

            ### Adding extra wavelengths to 'lamb' and extra flambdas to 'flambda'
            lamb = [elem for elem in lamb]
            flambda = [elem for elem in flambda];
            factor = np.nanmin([1.,(np.nanmax(lamb_r)-np.nanmin(lamb_r))/\
                        (np.nanmax(lamb_r)-np.nanmax(lamb))])
            Nlamb = int(np.trunc(factor*npts)+1.)
            auxi_lamb = np.linspace(np.nanmax([np.nanmin(lamb_r),np.nanmax(lamb)]),\
                            np.nanmax(lamb_r),Nlamb+1)
            auxi_lamb = np.array([auxi_lamb[ii] for ii in range(1,len(auxi_lamb))])
            for ii in range(0,Nlamb):
                lamb.append(auxi_lamb[ii])
                flambda.append(powerlaw(auxi_lamb[ii],C0,alpha0))
            lamb = np.array([elem for elem in lamb])
            flambda = np.array([elem for elem in flambda])


        if pf == "yes":
            dens=np.array(lamb)*np.array(flambda)/\
                (phc.h.cgs*phc.c.cgs)   ### lamb*flamb/hc [cgs units]
        else:
            dens=np.array(flambda)*1e8      ### flamb [erg cm^-2 s^-1 cm^-1]        
        
        ### domain of lambdas [Angstroms] for interpolation
        l=np.array([np.nanmin(lamb_r)+(np.nanmax(lamb_r)-\
            np.nanmin(lamb_r))/float(npts-1)*float(i) \
            for i in range(0,npts)])
        ### variation of lambda [Angstroms] for integration
        dl=np.array([l[i+1]-l[i] for i in range(0,npts-1)])
        
        ### interpolated lamb*flamb/hc [cgs units]
        d=np.array([interpLinND([l[i]],[lamb],dens) \
            for i in range(0,npts)])
        ### Interpolated passband
        r=np.array([interpLinND([l[i]],[lamb_r],R) \
            for i in range(0,npts)])
        ### photon flux [s^-1 cm^-2] or 
        ### mean flux density*bandwidth [erg cm^-2 s^-1]
        photflux=integrate_trapezia(d*r,dl*1e-8)
        return photflux
        
    else: 
        print("WARNING: Spectrum out of bounds of filter "+filtername+"!")
        print("Writing NaN to photon flux.")
        return np.nan




def iso_wavelength(lamb,flambda,U_meanfluxdensity,filtername):
    """
    Returns the isophotal wavelength [same unit as 'lamb']
    
    
    """
    
    U_mfd = U_meanfluxdensity   ### [Angstrom*'flambda']
    U = mean_RX(filtername,m_wavelength="no")
    mfd = U_mfd/U   ### mean flux density [same unit as 'flambda']
    
    y = flambda - mfd   ### y ['flambda']
    
    
    ### y = y0 +R*(lamb-lamb0)
    ### Hence, y = 0 ---> lamb = lamb0 - y0/R
    R = np.array([(y[i+1]-y[i])/(lamb[i+1]-lamb[i]) for i in range(0,len(y)-1)])
    roots = np.array([lamb[i]-y[i]/R[i] for i in range(0,len(y)-1)])
    
    rootsforiso = []
    for i in range(0,len(roots)):
        if lamb[i] <= roots[i] <= lamb[i+1]:
            rootsforiso.append(roots[i])
    
    if len(rootsforiso) > 0:
        lamb_iso = 0.
        for elem in rootsforiso:
            lamb_iso += elem
        lamb_iso = lamb_iso/float(len(rootsforiso))
        
        return lamb_iso
    else:
        return np.nan






def pogson(X,zp):
    return -2.5*np.log10(X)+zp



def VEGA_spct(spct_name):
    """
    Returns the SED of VEGA (9 nm - 160 microns):
    * array of lambdas [Angstroms]
    * array of Flambda [erg cm^-2 s^-1 A^-1]
    """

    if spct_name == 'spct1': ### Reference: 1994A&A...281..817C
        ### Reading Vega's spectrum
        spct1=cwd+'../refs/stars/fm05t9550g395k2odfnew.dat'
        f0 = open(spct1,'r')
        lines = f0.readlines()
        f0.close()

        ### Eliminating text and separating columns
        lines=[lines[i].split() for i in range(12,len(lines))]

        dist_fac=1./(1.62*10.**16.) ### Reference: 1994A&A...281..817C
        ### lambda [angstroms]
        lamb=np.array([float(lines[i][2])*10. for i in range(0,len(lines))])
        ### flambda [erg cm^-2 s^-1 A^-1]
        flambda=np.array([10.**(np.log10(4.*np.pi)+np.log10(phc.c.cgs)+\
                            np.log10(float(lines[i][4]))\
                            -2.*np.log10(float(lines[i][2])*1e-7)) \
                            for i in range(0,len(lines))])/1e8*\
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

def fullsed2photonflux(fullsed,source,filtername,npts=50,dist=10.,\
                        forcedred = 0, doit="yes"):
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
    ### Currently, this function is only usable for axisymmetric simulations.

    ### Reading the fullsed file
    f0 = open(fullsed,'r')
    lines = f0.readlines()
    f0.close()
    nlbd=int(lines[3].split()[0]); nobs=int(lines[3].split()[1])
    lines=[lines[i].split() for i in range(5,len(lines))]
    
    ### Reading the source file
    f0 = open(source,'r')
    src = f0.readlines()
    f0.close()
    lum=float(src[6].split()[2])
    normflx=(lum*phc.Lsun.cgs)/(4.*np.pi*(dist*phc.pc.cgs)**2.)

    ### Calculating the photon fluxes for each inclination
    mu=[];lamb=[];flambda=[];photflux=[]
    for iob in range(0,nobs):
        mu_aux=float(lines[0+nlbd*iob][0]) ### cosi
        lamb_aux=np.array([float(lines[ilamb+nlbd*iob][2])*1e4 \
            for ilamb in range(0,nlbd)]) ### lambda [Angstroms]
        flambda_aux=np.array([float(lines[ilamb+nlbd*iob][3])*\
            normflx/1e4 \
            for ilamb in range(0,nlbd)]) ### flambda [erg/cm^2 s A]
        if doit == "yes":
            photflux_aux=photonflux(lamb_aux,flambda_aux,filtername,npts,forcedred)
        else: 
            photflux_aux=np.nan
        ### Appending results for this specific inclination
        mu.append(mu_aux)
        lamb.append(lamb_aux)
        flambda.append(flambda_aux)
        photflux.append(photflux_aux)
    return mu,lamb,flambda,photflux



    








