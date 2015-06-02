#-*- coding:utf-8 -*-

"""
PyHdust *poltools* module: polarimetry tools

Em reuniões em nov/2014 com Moser, Bednarski e Alex, ficou decidido:

- Cada estrela (alvo ou padrão) terá um PREFIXO único de observação. Essa lista é gerada da seguinte forma: copia-se a coluna no arquivo *planejamento_LNA.xls* num arquivo de texto, que será lido pelas rotinas (arquivo *pyhdust/pol/alvos.txt*). Destaca-se que a identificação do alvo é feita desta forma, não importando o prefixo dos arquivos *fits*.

- A estrutura de diretórios "local" é *noite/prefixoalvo_n*. Assim, todo sub-diretório da noite é considerado um alvo, exceto *calib*. Caso uma estrela tenha sido observada em mais de uma sequencia, adicionar *_n*, onde pode ser um número (e.g., 2) ou um descritivo (e.g., *a2*). Tudo o que está em subdiretórios dos alvos é ignorado.

- A estrutura de diretórios será *dadospuros/noite/alvo* e *reduzidos/noite/alvo*. O script *pur2red* limpa o resultado da redução dos dados puros. O script *red2pur* copia os resultados da redução em reduzidos para a pasta dos dados puros.

-  O critério para a escolha do(s) arquivos *out(s)* será o de minimizar os erros em blocos de 16 ou 8+8 posições, mantendo o máximo de pontos independentes possível. Os detalhes estarão no manual das rotinas.

-  As rotinas python propostas são as seguintes: (i) *genStdLog*, (ii) *genObjLog*, (iii) *genStdDat*, (iv) *genTarget*, (v) *pur2red*, (vi) *red2pur* e (vii) *listNights*. Consultar o manual das rotinas para a lista completa e detalhes.

-  É possivel demonstrar que :math:`P=\sqrt{Q^2+P^2}` e :math:`\sigma_P=\sigma_{Q,U}`.

Passos da redução: 

-  Obtem-se a lista de noites do objeto (ou via página Beacon, ou via rotina *listNights* caso tenha acesso a pasta *puros*).

-  Gera-se os arquivos de calibração-padrão (rotina *iraf* *calib*).

-  Reduz-se as padrões da noite no *iraf*. Roda-se  *genStdLog*. A rotina imprime alertas de dados de padrões não reduzidos. Colunas da saída de *genStdLog*:

MJD|target|filt|calc|out|flags| (um *out* por filtro/padrão). 

-  Reduz-se os alvos da noite (*iraf*). Roda-se *genObjLog* e imprime alertas de dados não reduzidos/estrelas fora da nomenclatura. Colunas da saída de *genObjLog*: MJD|target|filt|calc|out\_n|flags| (pode haver mais de um *out* por filtro/alvo; registro em nova linha).

-  Roda-se *genStdDat*. A rotina faz: (i) olha na saída de *genObjLog* os filtros/calcitas dos alvos utilizados, (ii) e espera encontrar as padrões correspondentes; (iii) caso não encontre uma padrão, faz interpolação (de filtros). Se faltou para toda a calcita, pede para o usuário fazer a referência a saída de outro *genStdLog*.

-  Reduz-se todas as noites (de interesse). Roda-se *genTarget* para um, múltiplos ou todos os alvos. Procura noite à noite o alvo na saída de *genObjLog*, calibrado com *genStdDat* da noite correspondente (detalhe: nesta configuração, os arquivos *outs* serão lindos neste momento. Ou seja, o usuário depois dos passos 1, 2 e 3, ainda precisa ter acesso a estes arquivos). Neste momento, o usuário confirma qual ângulo de referência da padrão será utilizado (padrão é o "publicado"). Caso exista mais de uma padrão para o filtro/calcita, faz uma média do :math:`\Delta\theta`. Colunas da saída de *genTarget* (saída individual para cada estrela):

MJD|noite|filt|calc|ang.ref|del.th|P|Q|U|th|sigP|sigth|flags|.

Outras definições sugeridas:

-  No caso de um alvo observado com duas calcitas (como padrões), não colocar no mesmo diretório. Isso complica muito o desflexibiliza muito as rotinas que gerenciam os artivos *\*.out*. A regra é criar uma nova pasta com o sufixo "_a2" contendo os arquivos fits. A raiz destes arquivos é arbitrária.

-  Além disso, outra situação é necessária para se abrir outra pasta: ao se iniciar uma nova sequência temporal.

-  Nomenclatura de observação *target_a0_g1_f.fits* onde *a0* e *g1* são argumentos opcionais (calcita e ganho/config\*.\* do CCD) e *f* é o filtro. Na ausência destes parâmetros, assume-se valor padrão (ou lê-se no *header* dos fits).

-  Os avisos em tela *Warning* avisa o usuário que uma informação foi registrada, mas ela não é a otimizada. Por exemplo, na ausência de uma das versões de redução (.1 ou .2) ou do agrupamento de posições da lâmina (08 ou 16).

-  Os avisos em tela *Error* avisa o usuário que uma informação foi perdida. Dependendo o erro, o script será interrompido (não gerando nenhum output parcial).

Definições das flags.

- [0] As flags são cumulativas. Se estiver 0, significa que tudo ocorreu bem (nenhuma decisão ou aproximação no processo de redução foi feita). Se ocorreram, terá uma ou mais das flags abaixo.

- [pA] Estrela sem padrão para calibração de *dth* (output é salvo mesmo assim).

- [pB] A correção *dth* foi feita usando uma estrela de outra noite.

- [pC] A correção *dth* foi feita usando duas ou mais estrelas.

:co-author: Daniel Bednarski
:license: GNU GPL v3.0 (https://github.com/danmoser/pyhdust/blob/master/LICENSE)
"""
import os as _os
import glob as _glob
import numpy as _np
import datetime as _dt
from itertools import product as _product
from glob import glob as _glob
import pyhdust.phc as _phc 
import pyhdust.jdcal as _jdcal
from pyhdust import hdtpath as _hdtpath

try:
    import matplotlib.pyplot as _plt
    from matplotlib.transforms import offset_copy as _offset_copy
    import pyfits as _pyfits
except:
    print('# Warning! matplotlib and/or pyfits module not installed!!!')

__author__ = "Daniel Moser"
__email__ = "dmfaes@gmail.com"


filters = ['u','b','v','r','i']

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


def readoutMJD(out):
    """
    Read the *.out file from IRAF reduction and return a float array + str(MJD
    date) + float(angle) of the beams from calcite on the CDD. 

    readoutMJD
    """
    data = readout(out)
    WP = False
    if '_WP' in out:
        out = out[:out.rfind('_')]
        WP = True
    i = out.rfind('_')+1
    seq = int(out[i:i+2])
    npos = int(out[i+2:i+5])
    f = out[out.rfind('_')-1:out.rfind('_')]
    path = _phc.trimpathname(out)[0]
    JD = _glob('{0}/JD_*_{1}'.format(path,f))
    try:
        f0 = open(JD[0])
    except:
        print("# ERROR! No JD file found... {0}/JD_*_{1}".format(path,f))
        print('Script stops!!')
        raise SystemExit(1)
    date = f0.readlines()
    f0.close()
    datei = float(date[npos-1].split()[-1])-2400000.5
    datef = float(date[npos-1+seq-1].split()[-1])-2400000.5
    if WP:
        ver = out[-1]
    else:
        i = out[:-4].rfind('.')+1
        ver = out[i:i+1]
    coords = _glob('{0}/coord_*_{1}.{2}.ord'.format(path,f,ver))
    if len(coords) == 0:
        coords = _glob('{0}/coord_*_{1}.ord'.format(path,f))
        if len(coords) == 0:    
            print("# ERROR! No COORDS file found... {0}".format(out))
            print('Script stops!! {0}/coord_*_{1}.{2}.ord'.format(path,f,ver))
            raise SystemExit(1)
    coords = _np.loadtxt(coords[0])
    ang = _np.arctan( (coords[1,1]-coords[0,1])/(coords[1,0]-coords[0,0]) )*\
    180/_np.pi
    while ang < 0:
        ang += 180.
    if datei == datef:
        print('# Strange JD file for '+out)
    date = '%.7f' % ((datef+datei)/2)
    return [float(x) for x in data]+[date]+[ang]


def minErrBlk16(night,f,i):
    """Retorna o menor erro num bloco de 16 posicoes polarimetricas
    (i.e., arquivos [08(i),08(i+9)] e [16(i)]).
    """
    err = _np.ones(3)
    out = _np.zeros(3, dtype='|S256')
    ls = _glob('{0}/*_{1}_08{2:03d}*.out'.format(night,f,i))
    if len(ls) > 0:
        err[0] = float(readout(ls[0])[2])
        out[0] = ls[0]
        for outi in ls:
            if float(readout(ls[0])[2]) < err[0]:
                err[0] = float(readout(outi)[2])
                out[0] = outi
    else:
        print('# Warning! No 08pos *.out found at {0} for {1} band!'.format(\
        _phc.trimpathname(night)[1],f))
    ls = _glob('{0}/*_{1}_08{2:03d}*.out'.format(night,f,i+8))
    if len(ls) > 0:
        err[1] = float(readout(ls[0])[2])
        out[1] = ls[0]
        for outi in ls:
            if float(readout(ls[0])[2]) < err[0]:
                err[1] = float(readout(outi)[2])
                out[1] = outi
    ls = _glob('{0}/*_{1}_16{2:03d}*.out'.format(night,f,i))
    if len(ls) > 0:
        err[2] = float(readout(ls[0])[2])
        out[2] = ls[0]
        for outi in ls:
            if float(readout(ls[0])[2]) < err[0]:
                err[2] = float(readout(outi)[2])
                out[2] = outi
    #Soma na igualdade, pois pode procurar em i != 1 para npos=16.
    else:
        if len(_glob('{0}/*_{1}_*.fits'.format(night,f))) >= 16-1+i:
            print('# Warning! Some *_16*.out found at {0} for {1} band!'.format(\
            _phc.trimpathname(night)[1],f))
            #print i, ls, len(_glob('{0}/*_{1}_*.fits'.format(night,f)))
    j = _np.where(err == _np.min(err))[0]
    return err[j], out[j]


def chooseout(night,f):
    """
    Olha na noite, qual(is) *.OUT(s) de um filtro que tem o menor erro.

    Retorna um mais valores para toda a sequencia (i.e., pasta).
    Criterios definidos no anexo de polarimetria.

    O numero de posicoes eh baseado no numero de arquivos *.fits daquele
    filtro.
    """
    npos = len(_glob('{0}/*_{1}_*.fits'.format(night,f)))
    louts = _glob('{0}/*_{1}_*.out'.format(night,f))
    #check reducao
    if len(louts) == 0:
        outs = []
        if npos != 0:
            print('# Reducao nao feita em {0} para filtro {1}!'.format(\
            _phc.trimpathname(night)[1],f))
    #Se ateh 16 npos, pega o de menor erro (16 nao incluso)
    elif npos < 16:
        err1 = 1.
        for outi in louts:
            if float(readout(outi)[2]) < err1:
                err1 = float(readout(outi)[2])
                outs = [outi]
                #print 'a',outs
        if err1 == 1:
            print npos,louts
            print('# ERROR! Something went wrong with *.out Sigma values!!')
    #Se a partir de 16, faz o seguinte criterio: dentro de npos%8, ve qual
    #posicao inicial do block de 16 (ou 8+8) tem o menor erro, e joga no valor
    #de i.
    else:
        i = 1
        err1 = minErrBlk16(night,f,i)[0]
        #outs = list(minErrBlk16(night,f,i)[1])
        #print 'b',outs
        for j in range(1,npos%8+1):
            if minErrBlk16(night,f,j+1)[1] < err1:
                err1 = minErrBlk16(night,f,j+1)[0]
                i = j+1
                #outs = list(minErrBlk16(night,f,j+1)[1])
                #print 'c',outs
        outs = []
        while i+16-1 <= npos:
            outi = minErrBlk16(night,f,i)[1][0]
            if outi.find('_{0}_16'.format(f)) > -1:
                outs += [outi]
            else:
                for j in [i,i+8]:
                    ls = _glob('{0}/*_{1}_08{2:03d}*.out'.format(night,f,j))
                    if len(ls) != 2:
                        print(('# Warning! Check the *.out 2 '+\
                        'versions for filter {0} of {1}').\
                        format(f,_phc.trimpathname(night)[1]))
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
            outi = _glob('{0}/*_{1}_08{2:03d}*.out'.format(night,f,i))
            if len(outi) == 0:
                print('# ERROR! Strange *.out selection! Case 1')
                print('# Probably 08pos files missing.')
                print('{0}/*_{1}_08{2:03d}*.out'.format(night,f,i))
                print i,npos,npos-8+1,night,f
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
                    print('{0}/*_{1}_08{2:03d}*.out'.format(night,f,i))
                    print i,npos,npos-8+1,outi,night,f
                    raise SystemExit(1)
            for j in range(i+1,npos+1-8+1):
                tmp = _glob('{0}/*_{1}_08{2:03d}*.out'.format(night,f,j))
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
                    print('{0}/*_{1}_08{2:03d}*.out'.format(night,f,j))
                    print i,npos,npos-8+1,tmp,outi
                    raise SystemExit(1)                
                if float(readout(tmp)[2]) < err1:
                    err1 = float(readout(tmp)[2])
                    outi = tmp
            outs += [outi]
    return outs
    

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
    dtstart=dt.datetime(2012,7,1).date())
    mjdticks = [jdcal.gcal2jd(date.year,date.month,date.day)[1] for date in \
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
    dtstart=dt.datetime(2012,7,1).date())
    mjdticks = [jdcal.gcal2jd(date.year,date.month,date.day)[1] for date in \
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


def grafpol(argv):
    """
    Program to plot the best adjust and its residuals of IRAF reduction
    """    
    tela = False
    todos = False
    for i in range(len(argv)):
        if argv[i] == '--screen' or argv[i] == '-s':
            tela = True
        if argv[i] == '--all' or argv[i] == '-a':
            todos = True
    
    def arraycheck(refs,values):
        if len(refs) != len(values):
            print('# PROBLEM HERE!')
            print refs, values
            raise SystemExit(0)
        outvec = []
        for i in range(len(refs)):
            if refs[i] != '0':
                outvec += [i]
        return values[outvec]
    
    #ler *.log
    def readlog(filename,sigma):
        file = _np.loadtxt(filename, dtype=str, delimiter='\n')
        npts = int(file[8].split()[-1])
        delta = float(file[14].split()[-1])
        sigma = 1.
        for i in range(25, len(file)):
            if 'APERTURE' in file[i]:
                sig = float(file[i+2].split()[2])
                if sig < sigma:
                    sigma = sig
                    # Bed: Os Q e U sao os abaixo, conforme copiei da rotina graf.cl
                    if float(file[i+2].split()[4]) < 0:
                        thet = - float(file[i+2].split()[4])
                    else:
                        thet = 180. - float(file[i+2].split()[4])
                    Q = float(file[i+2].split()[3])*_np.cos(2.*thet*_np.pi/180.)
                    U = float(file[i+2].split()[3])*_np.sin(2.*thet*_np.pi/180.)
                    n = npts/4 
                    if npts%4 != 0:
                        n = n+1
                    P_pts = []
                    for j in range(n):
                        P_pts += file[i+4+j].split()
                    P_pts = _np.array(P_pts, dtype=float)
                    th_pts = 22.5*_np.arange(npts)-delta/2.
                    j = filename.find('.')
                    delta2 = int(filename[-2+j:j])-1
                    # Bed: Modifiquei abaixo para as pos lam >= 10 terem o num impresso corretamente no grafico
                    str_pts = map(str, _np.arange(1,npts+1)+delta2)
                    if int(file[9].split()[-1]) != npts:
                        refs = file[9].split()[3:-2]
                        #P_pts = arraycheck(refs,P_pts)
                        th_pts = arraycheck(refs,th_pts)
                        str_pts = arraycheck(refs,str_pts)
        if sigma == 1.:
            print('# ERROR reading the file %s !' % filename)
            Q = U = 0
            P_pts = th_pts = _np.arange(1)
            str_pts = ['0','0']
        return(Q, U, sigma, P_pts, th_pts,str_pts)
    
    def plotlog(Q,U,sigma,P_pts,th_pts,str_pts,filename):    
        
        fig = _plt.figure(1)
        fig.clf()
        ax1 = _plt.subplot(2, 1, 1)
        _plt.title('Ajuste do arquivo '+filename)
        _plt.ylabel('Polarizacao')
        
        ysigma = _np.zeros(len(th_pts))+sigma
        _plt.errorbar(th_pts,P_pts,yerr=ysigma, linewidth=0.7)
    
        th_det = _np.linspace(th_pts[0]*.98,th_pts[-1]*1.02,100)
        P_det = Q*_np.cos(4*th_det*_np.pi/180)+U*_np.sin(4*th_det*_np.pi/180)
    
        _plt.plot(th_det, P_det)
        _plt.plot([th_det[0],th_det[-1]], [0,0], 'k--')    
        ax1.set_xlim([th_pts[0]*.98,th_pts[-1]*1.02])
        _plt.setp( ax1.get_xticklabels(), visible=False)
        
        # Bed: Retirei o compartilhamento do eixo y com ax1, pois as escalas
        #devem ser independentes
        ax2 = _plt.subplot(2, 1, 2, sharex=ax1)
        _plt.xlabel('Posicao lamina (graus)')
        _plt.ylabel('Residuo (sigma)')
        
        P_fit = Q*_np.cos(4*th_pts*_np.pi/180)+U*_np.sin(4*th_pts*_np.pi/180)
        
        transOffset = offset_copy(ax2.transData, fig=fig, x = 0.00, y=0.10, units='inches')
        # Bed: Agora plota os residuos relativos (residuos divididos por sigma)
        _plt.errorbar(th_pts, (P_pts-P_fit)/sigma, yerr=1)
    
        for i in range(len(th_pts)):
            _plt.text(th_pts[i], (P_pts-P_fit)[i]/sigma, str_pts[i], transform=transOffset)  
        _plt.plot([th_det[0],th_det[-1]], [0,0], 'k--')  
        
        if tela:
            _plt.show()
        else:
            _plt.savefig(filename.replace('.log','.png'))
        return
    
    fileroot = raw_input('Digite o nome do arquivo (ou parte dele): ')
    print('# Lembre-se de que se houver mais de um arquivo com os caracteres acima')
    print("  sera' o grafico do que tiver menor incerteza.")
    print("# Para gerar todos os graficos, use a flag '-a'; para exibir na tela '-s'")
    
    fileroot = fileroot.replace('.log','')
    
    files = _glob('*'+fileroot+'*.log')
    
    if len(files) == 0:
        print('# Nenhum arquivo encontrado com *'+fileroot+'*.log !!')
        raise SystemExit(0) 
    
    sigmaf = 1.
    sigma = 1.
    for filename in files:
        Q, U, sigma, P_pts, th_pts, str_pts = readlog(filename,sigma)
        if todos:
            plotlog(Q,U,sigma,P_pts,th_pts,str_pts,filename)
            
        if sigma < sigmaf:
            Qf=Q
            Uf=U
            sigmaf=sigma
            P_ptsf=P_pts
            th_ptsf=th_pts
            str_ptsf=str_pts
            filenamef=filename
    
    if not todos:
        plotlog(Qf,Uf,sigmaf,P_ptsf,th_ptsf,str_ptsf,filenamef)
    
    print('\n# Fim da rotina de plot! Ver arquivo '+filenamef.replace('.log','.png'))
    return

   
def genStdLog(path=None, force=False, chknames=True):
    """
    Generate Standards Log

    If `force`==True, it overrides security condition.

    If `chknames`==True, unknown targets (i.e., folder name) will be manually
    checked. 
    """
    if path == None:
        path = _os.getcwd()

    # bednarski: security condition
    if not force:
        if not securCondit('{0}/std.log'.format(path)):
            return

    ltgts = _np.loadtxt('{0}/refs/pol_alvos.txt'.format(_hdtpath()), dtype=str)
    lstds = _np.loadtxt('{0}/refs/pol_padroes.txt'.format(_hdtpath()), dtype=str,\
    usecols=[0])
    tgts = [fld for fld in _os.listdir('{0}'.format(path)) if \
    _os.path.isdir(_os.path.join('{0}'.format(path), fld))]
    lines = ''
    rtgts=[]
    
    # Bednarski: I added below to verify and query correct standard names.
    # I changed to work on directories with suffix also (like 'eaql_a0').
    for tgt in tgts:
        if tgt.split('_')[0] not in _np.hstack((ltgts,lstds,_np.array(['calib'])))\
        and chknames:
            print('Target {0} is not a known standard or target!!\n'.\
            format(tgt.split('_')[0]) + '  If it\'s a standard, type the ' + \
            'common name. Otherwise, press -Enter- to skip:')
            while True:
                tgt_real = raw_input('  ')
                if tgt_real == '':
                    rtgts.append(tgt.split('_')[0])
                    break
                elif tgt_real in _np.hstack((lstds,_np.array(['calib']))):
                    rtgts.append(tgt_real)
                    break
                else:
                    print('{0} is not a known standard. '.format(tgt_real) + \
                    '(Wouldn\'t it be a target?)\n  Write again or just ' + \
                    'press -Enter- to skip:')
        else:
            rtgts.append(tgt.split('_')[0])

    i=0
    for prod in _product(rtgts,lstds,filters):
        tgt,std,f = prod
#        print tgt, i, i/(len(lstds)*len(filters))
        # num var is necessary to correct work when a target has more
        # than one directory in the same night
        tgtindex=i/(len(lstds)*len(filters))
        i+=1
#        if tgt.find(std) != 0:
#            print('teste {0}'.format(std))
        if tgt.find(std) == 0:
            outs = chooseout('{0}/{1}'.format(path,tgts[tgtindex]),f)
            for out in outs:
                Q,U,sig,P,th,sigT,ap,star,MJD,calc = readoutMJD(out)
                # Bednarski: I replaced the column concerning to the .out
                #    by the .out PATH.
                lines += '{0:12.6f} {1:>10s} {2:1s} {3:>5.1f} {4:>40s}\n'.\
                format(float(MJD), tgt, f, float(calc), \
                _os.path.relpath(out, path))
                
    if lines == '':
        print('# ERROR! No valid standards were found!')
    else:
        lines = '{:12s} {:>8s} {:1s} {:5s} {:>35s}\n'.format('#MJD','target',\
        'filt','calc','out')+lines
    f0 = open('{0}/std.log'.format(path), 'w')
    f0.writelines(lines)
    f0.close()
    return


# Bednarski: I added delta variable to (calc-calcst) tolerance
def chkStdLog(path=None, delta=2.5):
    """
    Generate Standards Data

    delta is the allowed variation for the angles between the two
    beams for one same calcite.
    """
    allinfo = True
    if path == None:
        path = _os.getcwd()
    #le `obj.log` e `std.log`
    std = _np.loadtxt('{0}/std.log'.format(path), dtype=str)
    obj = _np.loadtxt('{0}/obj.log'.format(path), dtype=str)
    #verifica se std tem mais de uma linha. Caso nao, aplica reshape
    if _np.size(std) == 5:
        std = std.reshape(-1,5)
    elif _np.size(std) < 5:
        print('# ERROR! Check `std.log` file!!')
    if _np.size(obj) == 5:
        obj = obj.reshape(-1,5)
    elif _np.size(obj) < 5:
        print('# ERROR! Check `std.log` file!!')
        #raise SystemExit(1)
    # verifica se todas as observacoes contem padroes,
    for obji in obj:
        foundstd = False
        f = obji[2]
        calc = float(obji[3])
        for stdi in std:
            fst = stdi[2]
            calcst = float(stdi[3])
            if f == fst and abs(calc-calcst) < delta:
                foundstd = True
        if not foundstd:
            #implementar ajuste de filtro e apontamento para outro std.log
            print(('# ERROR! Standard star info not found for filt. {0} and '+\
            'calc. {1} in {2}!').format(f, calc, path))
            allinfo = False
    return allinfo


def securCondit(chkfile):
    """ Ask the user if he wants to continue since the file already exists...

    OUTPUT: False (do not) and True (continue)
    """
    flag = True
    if _os.path.exists(chkfile):
        print('\n# CAUTION: the file {0} already exists. '.format(chkfile) +
        'Are you sure to continue, overwriting all manual ' + \
        'changes on this file? (y/n):')
        verif = ''
        while verif not in ['y','n']:
            verif = raw_input('  ')
            verif = verif.lower()
            if verif == 'n':
                print('# Aborted!')
                flag = False
            elif verif not in ['y','n']:
               print('Value not valid. Please, type \'y\' ou \'n\':')
    return flag
    

# Bednarski: I added delta variable to (calc-calcst) tolerance
def genObjLog(path=None, delta=2.5, force=False, chknames=True):
    """
    Generate Objects Log.
    Must be runned after genStdLog.

    delta is the allowed variation for the angles between the two
    beams for one same calcite.

    If `force`==True, it overrides security condition.

    If `chknames`==True, unknown targets (i.e., folder name) will be manually
    checked. 
    """
    if path == None:
        path = _os.getcwd()

    # bednarski: security condition
    if not force:
        if not securCondit('{0}/obj.log'.format(path)):
            return

    ltgts = _np.loadtxt('{0}/refs/pol_alvos.txt'.format(_hdtpath()), dtype=str)
    tgts = [fld for fld in _os.listdir('{0}'.format(path)) if \
    _os.path.isdir(_os.path.join('{0}'.format(path), fld))]
    lstds = _np.loadtxt('{0}/refs/pol_padroes.txt'.format(_hdtpath()), dtype=str,\
    usecols=[0])
    lines = ''
    rtgts=[]
    
    # Bednarski: I added below to verify and query correct target names.
    # I changed to work on directories with suffix also (like 'dsco_a0').
    for tgt in tgts:
        if tgt.split('_')[0] not in _np.hstack((ltgts,lstds,_np.array(['calib'])))\
        and chknames:
            print('Target {0} is not a known target either standard!!\n'.\
            format(tgt.split('_')[0]) + '  Type the common name ' + \
            'case it be a target or press -Enter- to skip:')
            while True:
                tgt_real = raw_input('  ')
                if tgt_real == '':
                    rtgts.append(tgt.split('_')[0])
                    break
                elif tgt_real in _np.hstack((ltgts,_np.array(['calib']))):
                    rtgts.append(tgt_real)
                    break
                else:
                    print('{0} is not a known target. '.format(tgt_real) + \
                    '(Wouldn\'t it be a standard?)\n  Write again or just ' + \
                    'press -Enter- to skip:')
        else:
            rtgts.append(tgt.split('_')[0])

    i=0
    for prod in _product(rtgts,ltgts,filters):
        tgt,star,f = prod
#        print tgt, i, i/(len(lstds)*len(filters))
        # num var is necessary to correct work when a target has more
        # than one directory in the same night
        tgtindex=i/(len(ltgts)*len(filters))
        i+=1
        if tgt.find(star) == 0:
            outs = chooseout('{0}/{1}'.format(path,tgts[tgtindex]),f)
            for out in outs:
                if out != '': #Erro bizarro... Nao deveria acontecer
                    Q,U,sig,P,th,sigT,ap,nstar,MJD,calc = readoutMJD(out)
                    # Bednarski: I replaced the column concerning to the .out
                    #    by the .out PATH.
                    lines += '{0:12.6f} {1:>10s} {2:1s} {3:>5.1f} {4:>40s}\n'.\
                    format(float(MJD), tgt, f, float(calc), \
                    _os.path.relpath(out, path))
    if lines == '':
        print('# ERROR! No valid targets were found!')
    else:
        lines = '{:12s} {:>8s} {:1s} {:5s} {:>35s}\n'.format('#MJD','target',\
        'filt','calc','out')+lines
    #check std info. Only if it is okay the `obj.log` will be written
    stdchk = False

    f0 = open('{0}/obj.log'.format(path), 'w')
    f0.writelines(lines)
    f0.close()
    if not chkStdLog(path=path, delta=delta):
        print('# ERROR! Missing information in the `std.log` file!!')

    # Bednarski: changed here  tgt  ->  tgt.split('_')[0]
    for tgt in rtgts: #implementar aqui
        if tgt.split('_')[0] not in _np.hstack((ltgts,lstds,_np.array(['calib']))):
            print('# Warning! Target {0} is not a known target!! '.\
            format(tgt) + 'Change manually if it isn\'t a standard.')
    return


# Bednarski: I added delta variable to (calc-calcst) tolerance
# It's really necessary the parameter target here?
def corObjStd(target, night, f, calc, path=None, delta=2.5):
    """
    Correlaciona um alvo (`target`) num dado filtro e calcita (`f` e `calc`) e
    retorna o valor da(s) padrao(oes) correspondente(s) em `std.log`.

    O angulo retornado eh o mesmo theta do .out e nao 180-theta!
    
    Tolerancia da calcita default: +/-2.5 graus

    Input: target, f, calc, stds
    Output: angref, thstd
    """
    if path == None:
        path = _os.getcwd()
    std = _np.loadtxt('{0}/refs/pol_padroes.txt'.format(_hdtpath()), dtype=str)
    calc = float(calc)
    angref = 0.
    thstd = 0.
    stdname = ''
    if _os.path.exists('{0}/{1}/std.log'.format(path,night)):
        stds = _np.loadtxt('{0}/{1}/std.log'.format(path,night), dtype=str)
        if len(stds) > 0 and len(stds[-1]) != 5:
            stds = stds.reshape(-1,5)
        k = 0
        angref = []
        thstd = []
        stdname = []
        sigs = []
        for stdinf in stds:
            if stdchk(stdinf[1])[0] and stdinf[2] == f and \
            abs(float(stdinf[3])-calc) <= delta:
                # Bednarski: I changed below to work correctly
                Q, U, sig, P, th, sigT, tmp, tmp2 = readout('{0}/{1}'.\
                    format(path+night,stdinf[4]))
                stdname += [ stdinf[1] ]
                thstd += [ float(th) ]
                if thstd == 0.:
                    thstd = 0.01
                sigs += [ float(sig)*100 ]
                sigth = 28.65*sig/P   # Nao faltava dividir por P?
                i = stdchk(stdinf[1])[1]
                j = filters.index(f)+1 #+1 devido aos nomes na 1a coluna
                angref += [ float(std[i,j]) ]
# Bednarski: Comentei abaixo pra retornar todas padrões
#        if angref != []:
#            if len(angref) == 1 and len(thstd) == 1:
#                thstd = thstd[0]
#                angref = angref[0]
#            else:
#            if len(angref) == len(thstd):
#                idx = sigs.index(np.min(sigs))
#                thstd = thstd[idx]
#                angref = angref[idx]
                #print('# ERROR! More than one std for {0}!'.format(night))
                #raise SystemExit(1)
                # AQUI fazer algoritmo pra quando houverem mais que uma padrão
    else:
        print('# ERROR! `std.log` not found for {0}'.format(night))
#    print stdname, thstd, angref
    return stdname, thstd, angref


# Bednarski: epssig define dentro de quantos sigmas uma polarizacao medida pode ser
# compativel com zero (e nao depender de dados de padrao polarizada)
def genTarget(target, path=None, PAref=None, skipdth=False, delta=2.5, epssig=2.0):
    """ Gen. target

    PAref deve ser do formato `pyhdust/refs/pol_padroes.txt`.

    Calculos/definicoes:

        * th_eq = -th_medido + th^pad_medido + th^pad_publicado
        * th_eq = 180 - th_medido - dth
        * dth = 180 - th^pad_medido - th^pad_publicado
        * Q = P*cos(2*th_eq*pi/180)
        * U = P*sin(2*th_eq*pi/180)
        * sigth = 28.6*sig

    `skipdth` ignores the delta-theta correction!
    """
    if path == None:
        path = _os.getcwd()
    if path[-1] != '/':
        path += '/'
    #Carregar angulos de referencia para todos os filtros de um conjunto de
    # padroes
    if PAref is not None:
        for line in PAref:
            if len(line) != 6 and len(line) != 10:
                print('# ERROR! Wrong PAref matrix format')
                raise SystemExit(1)
    else:
        PAref = _np.loadtxt('{0}/refs/pol_padroes.txt'.format(_hdtpath()), dtype=str)
    #le listas de padroes e alvos
    std = _np.loadtxt('{0}/refs/pol_padroes.txt'.format(_hdtpath()), dtype=str)
    obj = _np.loadtxt('{0}/refs/pol_alvos.txt'.format(_hdtpath()), dtype=str)
    #verifica se teu objeto eh um objeto valido
    if target not in _np.hstack((std[:,0],obj)):
        print('# Warning! Target {0} is not a default target or standard!!'.\
        format(target))
        tmp = raw_input('# Type something to continue...')
    #define noites
    nights = [fld for fld in _os.listdir(path) if \
    _os.path.isdir(_os.path.join(path, fld))]
    lines = ''
    for night in nights:
        #check obj.log for the night
        if _os.path.exists('{0}/{1}/obj.log'.format(path,night)):
            objs = _np.loadtxt('{0}/{1}/obj.log'.format(path,night), dtype=str)
            #verifica se std tem mais de uma linha. Caso nao, aplica reshape
            if _np.size(objs) == 5:
                objs = objs.reshape(-1,5)
            elif _np.size(objs) % 5 != 0:
                print('# ERROR! Check {0}/{1}/obj.log'.format(path,night))
                raise SystemExit(1)
            for objinf in objs:
                #~ thstd = 0.
                dth = []
                stdnames = ''
                #~ Existe algum motivo abaixo para encontrar se objinf[1] CONTEM target e nao
                #~ se eh de fato IGUAL a target? Comentei e mudei.
                #~ if objinf[1].find(target) > -1:
                if objinf[1] == target:
                    MJD, tmp, f, calc, out = objinf
                    # Bednarski: I changed below to work correctly
                    Q, U, sig, P, th, sigT, tmp, tmp2 = readout('{0}/{1}'.\
                    format(path+night,out))
                    P = float(P)*100
                    th = float(th)
                    sig = float(sig)*100
                    sigth = 0.
                    if P > 0:
                        sigth = 28.65*sig/P  # Nao precisava dividir por P ainda?
                    stdname, thstd, angref = corObjStd(target, night, f, calc, path=path, delta=delta)
                    #print f,  thstd, angref
                    # Bednarski: reorganizei esses ifs abaixo pra funcionar agora com
                    #     mais de uma padrao
                    if thstd==[] or angref==[]:
                        thstd = []
                        angref = []
                        mdth = 0.
                        devth = 0.
                        stdnames = 'no-stdstar'
                        ##print('# ERROR with thstd ou angref!')
                        ##print(night, f, calc)
                    else:
                        for i in range(0,len(thstd)):
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
                        lines += ('{:12s} {:>7s} {:1s} {:>5s} {:>12s} {:>6.1f} {:>4.1f}'+
                        '{:>7.3f} {:>6.3f} {:>6.3f} {:>6.1f} {:>5.3f} '+
                        '{:>5.1f}\n').format(MJD,\
                        night, f, calc, stdnames, mdth, devth, P, Q, U, th, sig, sigth)
                    else:
                        print('# ERROR! No std found for non-null polarization of ' + \
                        '{0}, filter {1}, in {2}'.format(target,f,night))
        else:
            print('# ERROR! `obj.log` not found for {0}'.format(night))
    #arquivo de saida
    if lines != '':
        print('# Output written: {0}/{1}.log'.format(path,target))
        f0 = open('{0}/{1}.log'.format(path,target),'w')
        lines = ('#MJD       night filt calc  stdstars  dth devdth   P      Q      U'+\
        '    th    sigP  sigth\n')+lines
        f0.writelines(lines)
        f0.close()
    else:
        print('# ERROR! No observation was found for target {0}'.format(target))
    return


def pur2red():
    """
    Pure to Reduced
    """
    return


def red2pur():
    """
    Reduced to Pure
    """
    return


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
                JD = _np.sum(jdcal.gcal2jd(*dtobs))+tobs
                JDout += 'WP {0}  {1:.7f}\n'.format(i,JD)
            f0 = open('JD_{0}'.format(pref),'w')
            f0.writelines(JDout)
            f0.close()
    return


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
    ph0 = float(ph0) - jdcal.MJD_0
    
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
    ax0.set_title('{0} ; P={1} days, ph0={2:.3f}'.format(tgt,P,ph0+jdcal.MJD_0))
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
    ax0.set_title('{0} (binned); P={1} days, ph0={2:.3f}'.format(tgt,P,ph0+jdcal.MJD_0))
    _plt.savefig('{0}/{1}_bin.png'.format(path,tgt))
    #_plt.show()
    return

def sortLog(file):
    """ Sort the *.out file """
    f0 = open(file)
    lines = f0.readlines()
    f0.close()
    log = _np.loadtxt(file, dtype=str)
    log = log[log[:,0].argsort()]
    fmt = '%12s %7s %1s %5s %5s %6s %5s %6s %6s %6s %5s %5s'
    #~ Format error
    #~ _np.savetxt(file.replace('.log','.txt'), log, fmt=fmt, header=lines[0])
    _np.savetxt(file.replace('.log','.txt'), log, fmt='%s', header=lines[0])
    return

def filtra_obs(n,obs):
    """ ### FILTER OBSERV. ### """
    nobs = [ ]
    for i in range(len(obs)):
        if obs[i][5]/obs[i][3] > n:
            nobs = nobs+[obs[i]]
    return _np.array(nobs)

def filtraobs(data, r=20):
    """ filtro! """
    idx = data['P']/data['sigP'] > r
    return data[idx]

def loadpol(txt):
    """ Load polarization txt file. """
    dtb = _np.loadtxt(txt, dtype=str)
    dtb = _np.core.records.fromarrays(dtb.transpose(), names='MJD,night,filt,\
    calc,stdstars,dth,devdth,P,Q,U,th,sigP,sigth', formats='f8,{0},{0},f8,{0},\
    f8,f8,f8,f8,f8,f8,f8,f8'.format(dtb.dtype))
    return dtb
    

### MAIN ###
if __name__ == "__main__":
    pass

