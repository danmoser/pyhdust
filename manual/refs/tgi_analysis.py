#!/usr/bin/env python
# -*- coding: utf-8 -*-

from subprocess import call
import matplotlib.pyplot as mpl
import numpy as np
import os
import re

deltaVelocity = 1000		# [km/s]
lightSpeed = 299792.458		# [km/s]
numberOfPoints = 120

# Linhas espectrais que iremos analisar. Use a segunda linha para testes.
#specLines = {'ha': 6562.8, 'he': 6678.2, 'si2': 6347.1, 'si3': 4553}
specLines = {'si3': 4553}	# para testes

# Regex's usadas para identificar parâmetros, como temperatura e gravidade su-
# perficial, nos nomes dos arquivos de modelos atmosféricos e separar colunas
# ao ler os arquivos de dados..
atmosphere = re.compile(r'BG([\d]{5})g([\d]{3})v2')
twoColumns = re.compile(r'[\s]*([\d,.]*)[\s]*([\d,.]*)[\s]*')

###############################################################################
# calc_limits() - calcula os extremos do intervalo válidos para a normalização,
# ou seja, os primeiros pontos de cada lado que não pertencem a nenhuma linha
# espectral. - OBSOLETA!
###############################################################################
def calc_limits(x, y, delta_y):
    size = len(y)
    min, max = (0, size - 1)
    mean = np.median(y)
    lsuccess = False
    rsuccess = False

    for i in range(size/4):
        ldelta = y[i]/mean
        rdelta = y[size-i-1]/mean
        if ldelta > 1 - delta_y/2 and ldelta < 1 + delta_y/2 and not lsuccess:
            min = i
            lsuccess = True
        if rdelta > 1 - delta_y/2 and rdelta < 1 + delta_y/2 and not rsuccess:
            max = size - i - 1
            rsuccess = True

    return min, max

###############################################################################
# check_limits() - checa se os extremos do intervalo são válidos para a norma-
# lização obtida por linfit(). Caso não sejam, calcula os primeiros pontos fo-
# ra da linha espectral. - OBSOLETA!
###############################################################################
def check_limits(x, y, iMin, iMax, delta_y):
    jMin, jMax = (iMin, iMax)
    size = iMax + 1
    for i in range(size/4):
        if y[i] > 1 - delta_y and y[i] < 1 + delta_y:
            jMin = i
        if y[size-i-1] > 1 - delta_y and y[size-i-1] < 1 + delta_y:
            jMax = size - i - 1
    if jMin != iMin and jMax != iMax:
        normalizedFlux = linfit(x, y, jMin, jMax)
        print('New limits!')
    else:
        print('Limits are good!')

###############################################################################
# idl() - executa um comando no idl.
###############################################################################
def idl(command):
    NULL = open('/dev/null', 'w')
    call(['idl', '-quiet', '-e', command], stdout=NULL, stderr=NULL)
    NULL.close()

###############################################################################
# genspec() - gera um espectro para o modelo e linha espectral passados como
# argumentos. Os dados são salvos na pasta 'spectra' no diretório atual.
###############################################################################
def genspec(model, specLine, deltaWl, resolution):
    spec = "synplot49,0,0,0,atmos='BGmodels_v2/%s',wstart=%d,wend=%d,wdist=%d,x,y" % (model, specLines[specLine] - deltaWl, specLines[specLine] + deltaWl, resolution)
    save = "openw,lun,'spectra/%s_%s',/get_lun & printf,lun,[x,y] & free_lun,lun" % (model, specLine)
    idl(spec + ' & ' + save)

###############################################################################
# linfit() - retorna um array (y) normalizado. - OBSOLETA!
###############################################################################
def linfit(x, y, iMin, iMax):
    new_y = y[iMin] + (y[iMax] - y[iMin]) * (x - x[iMin]) / (x[iMax] - x[iMin])
    avgValue = (new_y[iMin] + new_y[iMax]) / 2
    return new_y / avgValue

###############################################################################
# normalize() - normaliza um array de acordo com a mediana de seus valores. Não
# sei justificar o uso da mediana, mas acho que funciona ou é um bom ponto de
# partida.
###############################################################################
def normalize(array):
    new_array = array / np.median(array)
    return new_array

###############################################################################
# load_data() - carrega os espectros gerados anteriormente (da pasta 'spectra')
# em arrays do numpy. Retorna os dois arrays (x, y).
###############################################################################
def load_data(model, specLine):
    data = open('spectra/' + model + '_' + specLine, 'rt')

    wl_ = []			# vetores temporários, já que não é possível
    flux_ = []			# inicializar um array (numpy) vazio

    for line in data:
        x, y = twoColumns.match(line).groups()
        wl_.append(float(x))
        flux_.append(float(y))
    data.close()

    wl = np.array(wl_)
    flux = np.array(flux_)

    return wl, flux

###############################################################################
# analyze() - para determinado modelo de atmosfera e linha, devolve a largura
# equivalente e a profundidade da mesma.
###############################################################################
def analyze(model, specLine, deltaWl, resolution):
    if not os.path.isfile('spectra/' + model + '_' + specLine):
        genspec(model, specLine, deltaWl, resolution)

    wl, flux = load_data(model, specLine)

# As linhas comentadas a seguir eram usadas com as funções agora obsoletas.
#    size = len(wl)
#    iMin, iMax = calc_limits(wl, flux, 0.05)
#    normalizedFlux = linfit(wl, flux, iMin, iMax)
#    check_limits(wl, normalizedFlux, 0, size - 1, 0.05)

    size = len(wl)
    normalizedFlux = normalize(flux)

    temp, grav = atmosphere.match(model).groups()
    temp = int(temp)
    grav = int(grav)/100.0

    depth = 1 - min(normalizedFlux[int(size*0.25):int(size*0.75)])
#    width = eqwidth()		# TODO

    return temp, grav, depth


# Criação de uma lista com todos os modelos de atmosferas do Tlusty. A segunda linha apenas
# restringe a lista completa para fins de testes.
models_ = [model.strip('.5') for model in os.listdir('BGmodels_v2') if model.endswith('.5')]
models = [model for model in models_ if model.endswith('g250v2')]	# para testes
models.sort()

depth = {'ha': [], 'he': [], 'si2': [], 'si3': []}

#
# Comente as linhas abaixo se quiser trabalhar com o interpretador Python.
#

for specLine in specLines:
    deltaWl = deltaVelocity * specLines[specLine] / lightSpeed
    resolution = 2 * deltaWl / numberOfPoints
    temp = []
    for model in models:
        temp_, grav_, depth_ = analyze(model, specLine, deltaWl, resolution)
        temp.append(temp_)
        depth[specLine].append(depth_)
        print('%d, %.2f, %f') % (temp_, grav_, depth_)

#
# TODO!
# As linhas abaixo plotam os dados gerados usando a biblioteca matplotlib. Ainda faltam ajustes.
#

#mpl.plot(temp, depth['ha'], '', temp, depth['he'], '', temp, depth['si2'], '', temp, depth['si3'], '')
#mpl.plot(temp, depth['si3'], '')
#mpl.show()
