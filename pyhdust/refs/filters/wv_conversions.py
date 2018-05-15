#!/usr/bin/env python
# -*- coding:utf-8 -*-

import numpy as np
from pyhdust.tabulate import tabulate as tab

filt = 'Q'
fact = 1e4

data = np.loadtxt(filt+'.dat')
data[:, 0] *= fact

outf = tab(data, tablefmt="plain")

f0 = open(filt+'.dat', 'w')
f0.writelines(outf)
f0.close()
